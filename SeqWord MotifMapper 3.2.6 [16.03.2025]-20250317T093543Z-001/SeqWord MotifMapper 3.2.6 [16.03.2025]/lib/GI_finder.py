import sys, os, math, time, re, subprocess, tools, copy, operator
import blast, seq_io, reports, progressbar
from svg import SVG_linear as svg

########################################################################
class Main:
    def __init__(self, task_list, options, modified_sites=[], immutable_tasks=[], flg_exclude_strand="off", max_sites_for_verification=3000, 
        motif_setting="", completed_tasks=None, graph_title=""):
        # ATTRIBUTES
        self.IO = seq_io.IO()
        self.graph_title = graph_title
        self.task_list = self.tasks = task_list
        self.oCompletedTasks = completed_tasks
        self.sample_size = len(modified_sites)
        self.motif_setting = motif_setting
        self.frame = options["-l"]
        self.bigstep = options["-b"]
        self.mediumstep = options["-m"]
        self.smallstep = options["-s"]
        self.reference_step = options['-r']
        self.input_path = options["-i"]
        self.output_path = options["-o"]
        self.promoter_length = int(options["-p"])
        self.flg_saveFastaFile = options["-f"]
        self.flg_saveSVGimage = options["-v"]
        self.flg_doBLAST = options["-u"]
        self.flg_search_tRNA = options["-n"]
        self.scenario = options["-c"]
        self.flg_coding_sequences = False
        self.flg_distribution_stat = "Z-score" # LD | Z-score
        self.flg_exclude_strand = flg_exclude_strand
        self.max_sites_for_verification = int(max_sites_for_verification)
        
        self.modified_sites = modified_sites
        self.immutable_tasks = immutable_tasks
        self.silant_tasks = ["contigs",""]
        
        # BLAST bin and tmp folders
        self.BLAST_bin = os.path.join("lib","bin")
        self.BLAST_tmp = os.path.join(self.BLAST_bin,"tmp")
        
        self.echo = None
        if options["-e"] == "Contrasting/Iteration":
            self.flg_contrasting = True
            self.flg_iteration = True
        elif options["-e"] == "No":
            self.flg_contrasting = False
            self.flg_iteration = False
        elif options["-e"] == "Contrasting":
            self.flg_contrasting = True
            self.flg_iteration = False
        elif options["-e"] == "Iteration":
            self.flg_contrasting = False
            self.flg_iteration = True
            
        if self.flg_doBLAST.upper() == "YES" and self.scenario.find("MGE") > -1:
            self.flg_doBLAST = True
        else:
            self.flg_doBLAST = None
        
        self.StartTime = None
        self.DataSet = None
        self.StandardPatterns = {}
        self.oSVG = None
        self.oStatReport = None
        self.svg_reports = []
        self.stat_reports = []
        
        #### TMP
        self.calculate_GI = True

        # Validator
        self.oValidator = Validator()
        for task in self.tasks:
            if self.tasks[task]['scenario']:
                self.tasks[task]['scenario'] = self.oValidator.validateTasks(task_list=self.tasks[task]['scenario'], frame=self.frame, keep_initial_settings=True)

        # CONSTANTS
        self.symbols = {'GC-content':' %;',
                        'G/C-skew':';',
                        'A/T-skew':';',
                        'Pattern deviation':' %;',
                        'Absolute deviation':' %;',
                        'Pattern skew':' %;',
                        'Variance':';'}
        # IMPLEMENTATION
        self.initiateDataSet()

    # METHODS

    def initiateDataSet(self):
        self.DataSet = {'Sequence name':'',
                        'Sequence link':'',
                        'Sequence description':'',
                        'Total sequence length':0,
                        'Locus length':0,
                        'Left border':0,
                        'Frame':0,
                        'Step':0,
                        'Time':'',
                        'Tasks':{}
                        }
        
    # process source data files
    def process(self, gbk_files=[], immutable_tasks=[], modified_sites=[], top_indend=50, flg_exclude_strand=None, verified=False):
        # Set flag to exclude strans
        if flg_exclude_strand != None:
            self.flg_exclude_strand = flg_exclude_strand
        # Set immutable tasks and modified sites, if passed
        if immutable_tasks:
            self.immutable_tasks = immutable_tasks
        
        if modified_sites:
            self.modified_sites = modified_sites
        
        # Check if the list of input sequence files corresponds to the list of modified files
        if len(gbk_files) != len(self.modified_sites):
            raise ValueError(f"Length of input sequence files {len(GBK_files)} must equal the list of modified sites {len(self.modified_sites)}!")

        if not self.tasks:
            raise ValueError("No tasks were specified for this run!")
            
        # get file list
        if not gbk_files:
            filelist = self.get_fileList(self.input_path)
        else:
            filelist = gbk_files
        if len(filelist) == 0:
            raise ValueError("No input files were provided!")
            
        # loop files in the list
        outdata = {}
        for i in range(len(filelist)):
            modified_sites = self.modified_sites[i]
            self.sample_size = len(modified_sites)
            StartTime = time.perf_counter()  # Use perf_counter for high-resolution timing
            fname = filelist[i]
            # Input file is read and converted to seqlist dictionary:
            # {oGBK.description:{"sequence":str(oGBK.seq),"dataset":{'Gene map':dict(zip(cds_titles,cds)),"Accession":oGBK.accession}}}
            if fname[fname.rfind('.'):].upper() in (".GBK",".GB"):
                seqlist = self.getSequenceFromGBK(fname)
            elif fname[fname.rfind('.'):].upper() in (".GBF",".GBFF"):
                seqlist = self.getSeqFromGBFF(fname)
            elif fname[fname.rfind('.'):].upper() in (".FA",".FASTA",".FST",".FSA",".FAS",".FNA"):
                seqlist = self.getSeqFromFASTA(fname)
            else:
                raise ValueError(f"Input file extension {fname[fname.rfind('.')+1:]} is not recognised.\n" +
                "Use files with the following extensions:\n\tGBK, GB - for single record genbank files;" +
                "\n\tGBF, GBFF - for multiple record genbank files;" +
                "FA, FASTA, FST, FSA, FAS, FNA - for fasta files.")
            if not seqlist:
                raise ValueError(f"File {fname} has wrong format or corrupted!")
            
            # Process list of sequences
            seqlist_output = self.process_seqlist(seqlist=seqlist, step=self.bigstep, top_indend=top_indend, modified_sites=modified_sites, verified=verified)
            outdata[fname] = seqlist_output
        return outdata            
    
    # Process sequences listed in dictionary seqlist
    def process_seqlist(self, seqlist, step=None, frame=None, tasks=None, reference_step=None,
        top_indend=150, modified_sites=[], verified=True, flg_chatty=True):
        if tasks == None:
            tasks = self.tasks
            
        if step == None:
            step = self.bigstep
            
        if frame == None:
            frame = self.frame

        if reference_step == None:
            reference_step = self.reference_step
            
        # modified_sites are in the format [[start,end,strand,{data}],...]
        sites = [ls[:3] if not ls[3] else [int(ls[3]['Site']), int(ls[3]['Site']), ls[2]] for ls in modified_sites]
        
        if self.max_sites_for_verification and len(sites) > self.max_sites_for_verification:
            verified = False

        output = {}

        for seqname in seqlist:
            # Start time of sequence processing
            StartTime = time.perf_counter()
            # Gene coordinates can be used to calculate reference pattern based on only coding sequences
            genes = None
            if seqlist[seqname]['dataset']['Gene map']:
                genes = seqlist[seqname]['dataset']['Gene map']
            contigs = seqlist[seqname]['dataset']['contigs']
            seq = seqlist[seqname]["sequence"]
            if len(seq) < 5*frame:
                tools.alert(f"Sequence {seqname} of {len(seq)} bp is to short for slidong window of {frame} bp!")
                continue
                
            if seqlist[seqname]["dataset"] and seqlist[seqname]["dataset"]["Accession"]:
                acc = seqlist[seqname]["dataset"]["Accession"]
            else:
                acc = f"GI{tools.random_id(4)}"
                
            if flg_chatty:
                tools.msg(f"Processing sequence {seqname}")
            
            bar = None
            if flg_chatty:
                bar = progressbar.indicator(len(seq),"Process...")

            if len(seq) < 2 * reference_step:
                self.calculate(subseq=seq, full_length=len(seq), name=seqname, tasks=tasks, step=step, frame=frame, flg_update_references=True, 
                    flg_renew_dataset=True, genes=genes, progress_bar=bar)
            else:
                start = 0
                stop = reference_step
                while stop < len(seq):
                    locus = seq[start:stop]
                    self.calculate(subseq=locus, full_length=len(seq), name=seqname, tasks=tasks, step=step, frame=frame, shift=start, flg_update_references=True, 
                        flg_renew_dataset = True if start == 0 else False, genes=genes, progress_bar=bar)
                    start = stop
                    stop += reference_step

                # Remaining part of the sequence
                remaining_seqlength = len(seq) - stop + reference_step
                if remaining_seqlength > 5 * frame: # if remaining part of the sequence is smaller that 5 sliding windows, ignore this part
                    locus = seq[len(seq) - reference_step:]
                    start = len(seq) - remaining_seqlength
                    self.calculate(subseq=locus, full_length=len(seq), name=seqname, tasks=tasks, step=step, frame=frame, shift=start, 
                        flg_update_references=True, flg_renew_dataset = False, genes=genes, progress_bar=bar)
            if bar:
                bar.stop()
                        
            # Link calculated values with the respective tasks
            self.charge_tasks(modified_sites=sites)
            
            # Prform calculations for all tasks except 'MOD' and 'contigs' if available 
            if self.oCompletedTasks:
                for task in tasks:
                    if task in ('MOD', 'contigs'):
                        continue
                    tasks[task] = {}
                    tasks[task].update(self.oCompletedTasks[task])
                
            # Calculate Kendal correlation with immutable tasks
            self.calculate_correlation(tasks=tasks)
            
            # Calculate distribution across contigs
            if contigs:
                self.calcualte_distribution_across_contigs(contigs=contigs, modified_sites=sites, seqlength=len(seq), tasks=tasks)

            if "MGE" in tasks and not self.oCompletedTasks:
                # Identified genomic islands and overlapped genes are recorded to tasks[task]['frames'] as [['start-end':statistics],...] 
                # and to tasks[task]['genes'] as [{'start':start,'end':end,'tag':tag,'name':name,'product':product},...], respectively
                self.identify_GenomicIslands(seq=seq, task="MGE", tasks=tasks, genes=genes)
                # Calculate statistics of distribution of modified nucleotides acccross genomic islands
                self.calculate_distribution_statistics_for_GIs(seqlength=len(seq), modified_sites=sites, tasks=tasks)
                
            # Statistics of distribution over coding, promoter and non-coding regions
            # Performed if the annotation data is available and locations of modified sites are verified
            if genes and "MOD" in tasks and verified:
                results = self.calculate_CDS_nonCDS_distribution(seqlength=len(seq), modified_sites=modified_sites, genes=genes)
                if results:
                    tasks["MOD"][self.flg_distribution_stat] = {"coding":results['coding'],
                        "noncoding":results['noncoding'],"promoter":results['promoter']}
           
            # Prepare report
            self.report_output(seq=seq, seqname=seqname, acc=acc, tasks=tasks, genes=genes, top_indend=top_indend, verified=verified)

            elapsed_time = time.perf_counter() - StartTime  # Measure elapsed time
            print(f"has been done in {elapsed_time:.2f} sec.")
            print()
            
            # Mark tasks as completed
            for task in self.tasks:
                self.tasks[task]['completed'] = True
            output[seqname] = tasks

            return output
    
    def calculate_CDS_nonCDS_distribution(self, seqlength, modified_sites, genes):
        tools.msg("\nCalculate distribution of modified sites across coding, promoter, and non-coding regions...")
        if not len(modified_sites):
            raise ValueError("List of modified sites is empty!")
        
        sites = [ls[:3] if (len(ls) < 4 or not ls[3]) else 
            # either [start, stop, strand] or [location, location, strand]
            [int(ls[3]['Site']), int(ls[3]['Site']), ls[2]] 
            for ls in modified_sites]
        # Annotation data is either {} or {'Annotation': 'non-coding', 'Site location relative to TSC': 'n/a'}
        annotation_data = modified_sites[0][3]
        
        # Convert genes to a list of gene objects
        gene_list = list(genes.values())
        # This parameter doubles the genome length for both strands
        glength_correction = 2
        if self.flg_exclude_strand.upper() in ("LEADING","+"):
            gene_list = [g for g in gene_list if g.strand == 1]
            # Only one strand is used
            glength_correction = 1
        if self.flg_exclude_strand.upper() in ("LAGGING","-"):
            gene_list = [g for g in gene_list if g.strand == -1]
            # Only one strand is used
            glength_correction = 1

        # Calculate sequence lengths
        coding_seqlength = sum(abs(gene.location.end - gene.location.start) for gene in gene_list)
        promoter_locations, promoter_overlap = tools.select_promoter_regions(genes=gene_list, 
            promoter_length=self.promoter_length, seqlength=seqlength)
        promoter_seqlength = sum([ls[1] - ls[0] for ls in promoter_locations]) - promoter_overlap
        noncoding_seqlength = glength_correction * seqlength - coding_seqlength - promoter_seqlength

        # Count modified loci
        # In coding regions
        gene_locations = [[gene.location.start, gene.location.end, gene.strand] for gene in gene_list]
        if annotation_data:
            coding_counts = len([ls for ls in modified_sites if ls[3] and ls[3]['Annotation'] != 'non-coding' and 
                int(ls[3]['Site location relative to TSC']) > 0])
        else:
            coding_counts = tools.count_loci_in_genes(genes=gene_locations, modified_sites=sites)

        # In promoters\
        if annotation_data:
            promoter_counts = len([ls for ls in modified_sites if ls[3] and ls[3]['Annotation'] != 'non-coding' \
                and -self.promoter_length <= int(ls[3]['Site location relative to TSC']) <= 0])
        else:
            promoter_counts = tools.count_loci_in_genes(genes=promoter_locations, modified_sites=sites, 
                loci_to_filter=gene_locations, echo=False)
            
        noncoding_counts = len(modified_sites) - coding_counts - promoter_counts

        # Calculate total modified loci
        total_counts = len(modified_sites)
        
        # Calculate expectations
        coding_expectation = total_counts * coding_seqlength/glength_correction/seqlength
        promoter_expectation = total_counts * promoter_seqlength/glength_correction/seqlength
        noncoding_expectation = total_counts * noncoding_seqlength/glength_correction/seqlength

        if self.flg_distribution_stat == "Z-score":
            # Z-values and p-values
            coding_v, coding_error, coding_p, coding_expected, coding_observed = tools.calculate_zsore(target_seqlength=coding_seqlength, total_seqlength=(glength_correction * seqlength), 
                target_counts=coding_counts, total_counts=total_counts)
    
            promoter_v, promoter_error, promoter_p, promoter_expected, promoter_observed = tools.calculate_zsore(target_seqlength=promoter_seqlength, total_seqlength=(glength_correction * seqlength), 
                target_counts=promoter_counts, total_counts=total_counts)
    
            noncoding_v, noncoding_error, noncoding_p, noncoding_expected, noncoding_observed = tools.calculate_zsore(target_seqlength=noncoding_seqlength, total_seqlength=(glength_correction * seqlength), 
                target_counts=noncoding_counts, total_counts=total_counts)
    
        elif self.flg_distribution_stat == "LD":
            # LD-values and p-values
            coding_v, coding_error, coding_p = tools.contingency_table_statistics(seq_length = glength_correction * seqlength, region_length=coding_seqlength, 
                modified_in_regions = coding_counts, modified_outside_regions = len(modified_sites) - coding_counts, 
                unmodified_in_regions = coding_seqlength - coding_counts, unmodified_outside_regions = glength_correction * seqlength - len(modified_sites) - coding_seqlength + coding_counts)
    
            promoter_v, promoter_error, promoter_p = tools.contingency_table_statistics(seq_length = glength_correction * seqlength, region_length=promoter_seqlength, 
                modified_in_regions = promoter_counts, modified_outside_regions = len(modified_sites) - promoter_counts, 
                unmodified_in_regions = promoter_seqlength - promoter_counts, unmodified_outside_regions = glength_correction * seqlength - len(modified_sites) - promoter_seqlength + promoter_counts)
    
            noncoding_v, noncoding_error, noncoding_p = tools.contingency_table_statistics(seq_length = glength_correction * seqlength, region_length=noncoding_seqlength, 
                modified_in_regions = noncoding_counts, modified_outside_regions = len(modified_sites) - noncoding_counts, 
                unmodified_in_regions = noncoding_seqlength - noncoding_counts, unmodified_outside_regions = glength_correction * seqlength - len(modified_sites) - noncoding_seqlength + noncoding_counts)
                
        else:
            raise ValueError(f"Unknown statistics setting {self.flg_distribution_stat}!")
    
        if any([v == None for v in [coding_v, promoter_v, noncoding_v]]):
            return None
        # [value (Z-score or LD), std. error, p-value, expected counts, observed counts]
        return {
            'coding': [coding_v, coding_error, coding_p, coding_expectation, coding_counts],
            'promoter': [promoter_v, promoter_error, promoter_p, promoter_expectation, promoter_counts],
            'noncoding': [noncoding_v, noncoding_error, noncoding_p, noncoding_expectation, noncoding_counts]
        }    
        
    def report_output(self, seq, seqname, acc, top_indend=150, tasks=None, genes=None, verified=False):
        if tasks == None:
            tasks = self.tasks
            
        # Create report object
        self.oStatReport = reports.StatPlotReport(sample_size=self.sample_size, title=seqname, motif_setting=self.motif_setting, verified=verified)
        
        # Create SVG
        self.oSVG = svg(seqname=seqname, seqlength=len(seq), task_list=self.task_list, top_indend=top_indend, title=acc, graph_title=self.graph_title)
        
        # Adding tasks to SVG and StatReport
        for task in tasks:
            # Record CDS/nonCDS distribution
            if task == "MOD":
                if self.flg_distribution_stat not in tasks[task]:
                    raise ValueError(f"The requested {self.flg_distribution_stat} statistics has not been calculated")
                    
                self.oStatReport.set_CDS_nonCDS_stat(algorithm=self.flg_distribution_stat, 
                    statistics=tasks[task][self.flg_distribution_stat])   # self.flg_distribution_stat = Z-score | LD

            # Record contig distribution
            if task == "contigs":
                self.oStatReport.set_contig_statistics(contigs=tasks[task]['contig borders'],
                    anova_p_value=tasks[task]['p-value'])
                
            if task in self.immutable_tasks or task in self.silant_tasks:
                continue

            if task == "MGE" and tasks["MGE"]["frames"]:
                # Add MGE statistics to report
                self.oStatReport.set_MGE_stat(**dict(zip(['Z_score', 'std_error', 'p_value', 'expected_number', 'observed_number'], tasks["MGE"]['Z-score'])))
                self.oStatReport.set_MGE_loci(locations = tasks["MGE"]["frames"], genes = tasks["MGE"]['genes'], seq = seq)
                # add_gi(self,lb,rb,color="",title="")
                for i in range(len(tasks["MGE"]["frames"])):
                    gi = tasks["MGE"]["frames"][i]
                    lb, rb = self.oValidator.parse_location(gi[0])
                    title = f"GI:{i+1} [{lb}..{rb}]"
                    self.oSVG.add_gi(lb=lb, rb=rb, title=title)
            elif tasks[task]['frames']:
                statistics={}
                if 'correlation' in tasks[task] and tasks[task]['correlation'] and self.immutable_tasks:
                    statistics = {"correlation":tasks[task]['correlation'][self.immutable_tasks[0]]}
                self.oSVG.add_task(task=task, task_description=tasks[task]['description'], windows=tasks[task]['frames'],
                    statistics=statistics)
                self.oStatReport.set_task(title=task, statistics=statistics)
            else:
                tools.alert(f"\t\tTask {task} has not been added!")
        self.oSVG.set_svg(verified=verified) 
        return True
        
    def get_svg(self):
        return self.oSVG

    def get_svg_code(self):
        return self.oSVG.get()
        
    def get_stat_report(self):
        return self.oStatReport.get()

    def get_GISeq(self,lb,rb,seq,flg_circular=True):
        before = after = ""
        start = int(lb)
        end = int(rb)
        if lb < 1 and not flg_circular:
            lb = start = 1
        elif lb < 1 and flg_circular:
            start = 1
            before = seq[len(seq)+lb:]
        elif lb==0:
            lb = start = 1
            
        if rb > len(seq) and not flg_circular:
            rb = end = len(seq)
        elif rb > len(seq) and flg_circular:
            after = seq[:rb-len(seq)]
            end = len(seq)
        GISeq = "".join([before,seq[start-1:end],after])
        return lb,rb,GISeq
        
    # Calculate kendal correlation between immutable tasks and other tasks
    def calculate_correlation(self, tasks=None, method=None):
        '''
        Calculate associations between 2 arrays of values: values_1 
        Methods:
            RAS - rank association analysis using Kendall tau correlation and Mann-Whitney U Test test through defined ranks (rank_number = 4 by default)
            U - Mann-Whitney U Test (default);
            KS - Kolmogorov-Smirnov Test;
            TAU - Kendall tau correlation;
            SP - Spearman correlation;
            MI - Mutual Information (k-nearest neighbors)
        '''
        # Set or check correlation method
        methods = {
            "RAS":tools.calculate_rank_association,
            "U":tools.u_correlation,
            "KS":tools.ks_correlation,
            "TAU":tools.kendal_correlation,
            "SP":tools.spearman_correlation,
            "MI":tools.calculate_mutual_info
        }
        if not method:
            method = "SP"
        if method not in list(methods.keys()):
            raise ValueError(f"Unknown correlation method abbreviation {method}!")
            
        # Check if any immutable tesks, like "MOD", were set
        if not self.immutable_tasks:
            return False

        if tasks == None:
            tasks = self.tasks
            
        for qtask in tasks:
            if qtask in self.immutable_tasks or qtask in ("MGE", ):
                continue
            for stask in self.immutable_tasks:
                values_1 = [float(ls[1]) for ls in sorted(tasks[qtask]['frames'], 
                    key=lambda ls: self.oValidator.parse_location(ls[0]))]  # Sliding window properties like GC-content or GC-skew
                values_2 = [float(ls[1]) for ls in sorted(tasks[stask]['frames'], 
                    key=lambda ls: self.oValidator.parse_location(ls[0]))]  # Counts of modified nucleotides across sliding windows
                results = methods[method](values_1, values_2, parameter=qtask)
                if results:
                    correlation = results['correlation']
                    error = results['correlation error']
                    p_value = results['p-value']
                    if 'correlation' not in tasks[qtask]:
                        tasks[qtask]['correlation'] = {}
                    tasks[qtask]['correlation'][stask] = [correlation, error, p_value, method]
                
    # Calculate statistics of distribution of modified sites genomic islands
    def calculate_distribution_statistics_for_GIs(self, seqlength, modified_sites, tasks=None):
        if not self.immutable_tasks or "MOD" not in self.immutable_tasks:
            return False
        if tasks == None:
            tasks = self.tasks
            
        sites = [ls[:3] for ls in modified_sites]
        
        regions_of_interest = [self.oValidator.parse_location(location) for location in [ls[0] for ls in tasks["MGE"]["frames"]]]
        results = tools.calculate_distribution_statistics(seq_length=seqlength, 
            modified_nucleotides=sites, 
            regions_of_interest=regions_of_interest,
            mode=self.flg_distribution_stat)

        if results and self.flg_distribution_stat == "LD":
            LD = results['LD']
            error = results['LD error']
            p_value = results['p-value']
            tasks["MGE"]["LD"] = [LD, error, p_value, 0, 0]
        elif results and self.flg_distribution_stat == "Z-score":
            Z_score = results['Z-score']
            error = results['Z error']
            p_value = results['p-value']
            expected_count = results['expected']
            observed_count = results['observed']
            tasks["MGE"]["Z-score"] = [Z_score, error, p_value, expected_count, observed_count]
    
    # Calculate statistics of distribution of modified sites across contigs
    def calcualte_distribution_across_contigs(self, contigs, modified_sites, seqlength, tasks=None):
        if tasks == None:
            tasks = self.tasks
        if not contigs:
            return False
        contig_borders = [[contig.location.start,contig.location.end] for contig in contigs]
        
        sites = [ls[:2] for ls in modified_sites]
        
        nonrandom_site_distribution_across_contigs = tools.compare_sites_across_contigs(contigs=contig_borders, modified_sites=sites)
        if "contigs" not in tasks:
            tasks['contigs'] = {'contig borders':[], 'p-value':nonrandom_site_distribution_across_contigs, 'whole genome length':seqlength}
        tasks['contigs']['contig borders'] = [{"title":contig.qualifiers['name'][0] if 'name' in contig.qualifiers else "",
            "start":contig.location.start,
            "end":contig.location.end,
            "Z-score":tools.calculate_zsore(target_seqlength=abs(contig.location.end-contig.location.start),    # [z_score, std_error, p_value, expectation, observation]
                total_seqlength=seqlength, 
                target_counts = len([site for site in modified_sites if site[0] > contig.location.start and site[1] <= contig.location.end]), 
                total_counts=len(modified_sites))} 
            for contig in contigs]

    # Calculate values for requested tasks
    def charge_tasks(self, modified_sites=[], tasks=None):
        if tasks == None:
            tasks = self.tasks
            
        def charge_task_with_data(D, task, subtask=""):
            # Mapping operand symbols to operator functions
            operations = {
                "/": operator.truediv,
                "*": operator.mul,
                "+": operator.add,
                "-": operator.sub
            }            
            if subtask:
                D = D[subtask]
                task = subtask
            D['frames'] = {}
            # Parse task
            subject, symbol, operand = self.oValidator.parse_task(task)

            if operand == None:
                # Not a complex task, simply link calculated values in self.DataSet with the respective task
                if subject not in self.DataSet['Tasks']:
                    raise ValueError(f"There are no values for the task {subject}!")
                D['frames'] = [[ls[0],float(ls[1]['value']) if isinstance(ls[1],dict) else float(ls[1])] for ls in self.DataSet['Tasks'][subject].items()]
            elif isinstance(operand,float):
                # Complex task like val / 100
                if subject not in self.DataSet['Tasks']:
                    raise ValueError(f"There are no values for the task {subject}!")
                D['frames'] = [[ls[0],operations[symbol](float(ls[1]['value']) if isinstance(ls[1],dict) else float(ls[1]), operand)] for ls in self.DataSet['Tasks'][subject].items()]
            elif isinstance(subject,float):
                # Complex task like 100 * val
                if operand not in self.DataSet['Tasks']:
                    raise ValueError(f"There are no values for the task {operand}!")
                D['frames'] = [[ls[0],operations[symbol](subject, float(ls[1]['value']) if isinstance(ls[1],dict) else float(ls[1]))] for ls in self.DataSet['Tasks'][operand].items()]
            elif isinstance(operand,str):
                if subject not in self.DataSet['Tasks'] or operand not in self.DataSet['Tasks']:
                    raise ValueError(f"There are no values either for the task {subject} or for the task {operand}!")
                D['frames'] = [[ls1[0],operations[symbol](float(ls1[1]['value']) if isinstance(ls1[1],dict) else float(ls1[1]), float(ls2[1]['value']) if isinstance(ls2[1],dict) else float(ls2[1]))] for ls1, ls2 in zip(self.DataSet['Tasks'][subject].items(), self.DataSet['Tasks'][operand].items())]
        
        for task in tasks:
            if task == "MOD":
                print("\nCalculate density of modified nucleotides")
                sliding_windows = sorted(list(self.DataSet['Tasks'][task].values()), key=lambda ls: int(ls[0]))
                counts, fdr_cutoff_p = tools.density_distribution(frames=sliding_windows, 
                    modified_sites=modified_sites)
                # Convert modified_site_distribution to format [['start-end':value],...]
                tasks[task]['frames'] = [[f"{sliding_windows[i][0]}-{sliding_windows[i][1]}", counts[i]] for i in range(len(counts))]
                tasks[task]['fdr_cutoff_p'] = fdr_cutoff_p
            
            elif tasks[task]['scenario']:
                tasks[task]['frames'] = {}
                for subtask in tasks[task]['scenario']:
                    charge_task_with_data(tasks[task]['scenario'], task, subtask)
                    # Link two possible locations of data access
                    tasks[task]['frames'][subtask] = tasks[task]['scenario'][subtask]['frames']
                    
            else:
                charge_task_with_data(tasks[task], task)
    
    def getAnnotation(self,data,lb,rb):
        def overlap(gene,lb,rb):
            coords = gene.split(" | ")[-1][1:-1]
            start,end = [int(v) for v in (coords.split("..") if coords.find("..") > -1 else coords.split("-"))]
            if (start >= int(lb) and start < int(rb)) or (end >= int(lb) and end < int(rb)):
                return True
            return False
            
        genes = [gene for gene in data if overlap(gene, lb, rb)]
        return "\t" + ";\n\t".join(genes), lb, rb

    def getSeqFromFASTA(self,fname,names=[]):
        oFASTA = self.IO.parse(fname,"fasta")
        return {oFASTA.description:{"sequence":oFASTA.seq,"dataset":None}}

    def getSequenceFromGBK(self,fname):
        def get_gene_description(g):
            tags = ['locus_tag','gene','product']
            values = [g.qualifiers[tag][0].strip().strip("'") if tag in g.qualifiers else "." for tag in tags]
            values = [s for s in values if s] + [f"[{g.location.start}..{g.location.end}]"]
            return " | ".join(values)
             
        oGBK = self.IO.parse_seq(fname,"genbank",concatenate=True)
        if not oGBK:
            raise ValueError(f"File {fname} is empty or corrupted!")
        contigs = [feature for feature in oGBK.features if feature.type=="contig"]
        cds = [feature for feature in oGBK.features if feature.type=="CDS"]
        if not cds:
            cds = [feature for feature in oGBK.features if feature.type=="gene"]
        cds_titles = [get_gene_description(g) for g in cds]
        return {oGBK.description:{"sequence":str(oGBK.seq),"dataset":{'contigs':contigs,'Gene map':dict(zip(cds_titles,cds)),"Accession":oGBK.accession}}}

    def getSeqFromGBFF(self,fname):
        return getSequenceFromGBK(fname)

    def parse_seqname(self,line,seq_names):
        if line > 3 and line[:3]=="gi|":
            seqname_elements = str.split(line,"|")
            if len(seqname_elements)==5:
                seqname = seqname_elements[4]
                pos = str.find(seqname,", complete sequence")
                if pos > -1:
                    seqname = seqname[:pos]
                if str.find(seqname_elements[3],"NC_")==0:
                    accession = seqname_elements[3]
                    pos = str.find(accession,".")
                    if pos > -1:
                        accession = accession[:pos]
                    seqname += " [" + accession + "]"
            else:
                seqname = seqname_elements[0]
        else:
            seqname = line
        seqname = self.check_seqname(seqname,seq_names)
        return seqname
    
    def check_seqname(self,seqname,names=[]):
        mark = str.rfind(seqname,", complete")
        if mark > 0:
            seqname = seqname[:mark]
        for symb in ("\\","/",":","*","?","\"","<",">","|"):
            seqname = str.replace(seqname,symb," ")
        while seqname and seqname[0] == " ":
            seqname = seqname[1:]
        pos = str.find(seqname,", complete")
        if pos != -1:
            seqname = seqname[:pos]
        while seqname and (seqname[-1] == "." or seqname[-1]) == " ":
            seqname = seqname[:-1]
        if not seqname:
            seqname = "#1"
        if not names:
            return seqname
        while str.upper(seqname) in names:
            sep = str.rfind(seqname,"#")+1
            if sep == 0:
                seqname += " #1"
            else:
                try:
                    n = int(seqname[sep:])
                except:
                    return seqname + " #1"
                seqname = seqname[:sep] + str(n+1)
        return seqname
       
    def identify_GenomicIslands(self, seq, window=None, bigstep=None, smallstep=None, ref_pattern=None, tasks=None,
        genes=None, task="MGE", mode="relaxed", main_task='n0_4mer:D'):    # Mode cane be relaxed | average | stringent
        if bigstep == None:
            bigstep = self.bigstep
        if smallstep == None:
            smallstep = self.smallstep
        if tasks == None:
            tasks = self.tasks
        if window == None:
            window = self.frame - bigstep
            
        def get_precise_borders(loci, seq, window, bigstep, smallstep, main_task, mode="relaxed", ref_pattern=None):
            # Parse k-mer pattern description to paraemeters
            task_title, taskID, ptype, norm, wlength = self.oValidator.parse_task_title(main_task)
            subseq_number = 0
            for i in range(len(loci)-1,-1,-1):
                coords, stat = loci[i]
                start, end = new_coords = self.oValidator.parse_location(coords)
                locus_seq = seq[start:end]
                lb, rb = extend_borders(start, end, len(seq))
                extended_locus_seq = seq[lb:rb]
                # Set or reset a reference pattern
                subseq = ""
                if self.reference_step and len(seq) > 3*self.reference_step:
                    new_subseq_number = lb//self.reference_step
                    if new_subseq_number != subseq_number:
                        subseq = (seq[new_subseq_number*self.reference_step:(new_subseq_number+1)*self.reference_step] if (new_subseq_number+1)*self.reference_step <=len(seq) else 
                            seq[new_subseq_number*self.reference_step - ((new_subseq_number+1)*self.reference_step -len(seq)):len(seq)])
                    elif ref_pattern == None:
                        subseq[:self.reference_step]
                elif ref_pattern == None:
                    subseq = seq
                if subseq:
                    ref_pattern = self.getPattern(seq=subseq, wlength=wlength, norm=norm, ptype=ptype)

                gMin, gMax, gAvr = stat[main_task]

                # Check left border
                success = False
                values = []
                begin = 0
                finish = begin + window
                while finish + 1 < len(extended_locus_seq):
                    value = ref_pattern - self.getPattern(seq=extended_locus_seq[begin:finish+1], 
                        wlength=wlength, norm=norm, ptype=ptype)
                    values.append([lb+begin, lb+finish, int(value)])
                    if ((mode == "relaxed" and value >= gMin) or (mode == "average" and value >= gAvr) or 
                        (mode == "stringent" and value >= gMax)):
                            new_coords[0] = (begin + finish + 2*lb)//2
                            success = True
                            break
                    begin += smallstep
                    finish = begin + window

                # Check right border
                if success:
                    begin = len(extended_locus_seq) - window
                    finish = len(extended_locus_seq) - 1
                    while begin > 0:
                        value = ref_pattern - self.getPattern(seq=extended_locus_seq[begin:finish+1], 
                            wlength=wlength, norm=norm, ptype=ptype)
                        if ((mode == "relaxed" and value >= gMin) or (mode == "average" and value >= gAvr) or 
                            (mode == "stringent" and value >= gMax)):
                                new_coords[1] = (begin + finish + 2*lb)//2
                                break
                        begin -= smallstep
                        finish = begin + window
                    loci[i][0] = f"{new_coords[0]}-{new_coords[1]}"
                else:
                    loci.pop(i)
            return loci
        
        def extend_borders(start, end, seqlength):
            lb, rb = sorted([int(start), int(end)])
            length = rb - lb
            if not length:
                return lb, rb
            window = length // 2
            lb -= window
            if lb < 0:
                lb = 0
            rb += window
            if rb >= seqlength:
                rb = seqlength-1
            return lb, rb
        
        # Identify GIs
        print("\nIdentification of genomic islands...")
        loci = self.identify_GI_Borders(task=task)
        print(f"\nIn total, {len(loci)} candidate genomic islands were identified.")
        
        # Precise identification of borders
        print("Checking genomic island borders...")
        loci = get_precise_borders(loci=loci, seq=seq, window=window, bigstep=bigstep, smallstep=smallstep, 
            main_task=main_task, mode=mode)
        print(f"\nIn total, {len(loci)} genomic islands were identified.")
        
        if self.flg_doBLAST:
            loci = self.filter_rrn(loci=loci, seq=seq)
            print(f"\nAfter rrn filtering, {len(loci)} genomic islands have remained.")
        
        gi_genes = None
        if genes:
            print("\nSelect GI's genes and align GI coordinates with overlapped genes.")
            gene_list = sorted(self.oValidator.gene_strings_to_dicts(genes), key=lambda gene: gene['start'])
            loci, gi_genes = self.select_gi_genes(loci, gene_list)
        
        # GI loci are added as [['start-end',statistics],...]
        tasks[task]['frames'] = [[f"{ls[0]}-{ls[1]}",ls[2]] for ls in loci]
        # Genes for each GI are added as [{'start':start,'end':end,'tag':tag,'name':name,'product':product},...]
        tasks[task]['genes'] = gi_genes
        return True
        
    def filter_rrn(self, loci, seq, frame=None, p_cutoff=0.0001):
        if not frame:
            frame = self.frame
        # Create BLAST object
        oBlast = blast.BLAST("dna", binpath=self.BLAST_bin, source_path=self.BLAST_tmp)
        if not os.path.exists(self.BLAST_tmp):
            os.mkdir(self.BLAST_tmp)
        for i in range(len(loci)-1,-1,-1):
            coords, stat = loci[i]
            start, end = self.oValidator.parse_location(coords)
            # Extract GI sequene
            query_seq = seq[start:end]
            # Create temporary fasta file for the sequence
            while True:
                # Select temporary file name
                query_fname = f"{tools.random_id(6)}.tmp"
                if os.path.exists(os.path.join(self.BLAST_tmp,query_fname)):
                    continue
                break
            self.IO.save(f">tmp\n{query_seq}", os.path.join(self.BLAST_tmp,query_fname))
            # Perform balst search through 16S database
            oBlast.execute(query_fname,"16S")
            blast_result = oBlast.get_top_alignment()
            # Parse blast result
            # If rrn is found at the end of GI, this part is removed
            # If rrn takes the central part of GI, GI is removed
            if blast_result:
                p, score, sbjct_title, hit_list = blast_result
                # Select the largest hit coordiantes
                lb = hit_list[0].query_start
                rb = hit_list[0].query_end
                if float(p) < p_cutoff:
                    if lb < frame and end - start - rb < frame:
                        del loci[i]
                        continue
                    if lb < frame:
                        loci[i] = [start - rb + 1, end, stat]
                    elif end - start - rb < frame:
                        loci[i] = [start, end - rb + lb, stat]
            # Remove temporary file
            if os.path.exists(os.path.join(self.BLAST_tmp,query_fname)):
                os.remove(os.path.join(self.BLAST_tmp,query_fname))
        return loci
        
    def select_gi_genes(self, loci, gene_list):
        gi_genes = []
        for i in range(len(loci)):
            coords, stat = loci[i]
            start, end = self.oValidator.parse_location(coords)
            loci[i] = [start, end, stat]
            selected_genes = []
            for gene in gene_list:
                # Get gene start and end
                if (gene['start'] >= start and gene['start'] <= end) or (gene['end'] >= start and gene['end'] <= end):
                    selected_genes.append(gene)
            if selected_genes:
                # Set GI start and end if previous borders were within overlapped genes
                loci[i] = [selected_genes[0]['start'] if selected_genes[0]['start'] < start else start, 
                    selected_genes[-1]['end'] if selected_genes[-1]['end'] > end else end, stat]
            gi_genes.append(selected_genes)
        return loci, gi_genes
    
    def join_regions(self,first,second):
        box = [first[0],second[1],{}]
        for key in first[2]:
            a = float(first[1]-first[0])
            b = float(second[1]-second[0])
            c = float(first[1]-first[0]+second[1]-second[0])
            box[2][key] = (first[2][key]*a + second[2][key]*b)/c
        return box

    def getExactValue(self,left,seq,mode):
        if len(seq) < self.frame:
            return left
        self.initiateDataSet()
        self.calculate(seq,"",self.smallstep,None,False)
        key = list(self.DataSet["Tasks"].keys())[0]
        coordinates = list(self.DataSet["Tasks"][key].keys())
        coordinates.sort(key=lambda s: int(s.split("..")[0]) if s.find("..") > -1 else int(s.split("-")[0]))
        for i in range(1,len(coordinates)):
            isSignal = 1
            for task in self.tasks:
                subtr = self.tasks[task]['subtr']
                if subtr:
                    subtrahend = self.DataSet["Tasks"][subtr][coordinates[i]]['value']
                else:
                    subtrahend = 0
                divisorId = self.tasks[task]['divisor']
                if divisorId:
                    divisor = self.DataSet["Tasks"][divisorId][coordinates[i]]['value']
                else:
                    divisor = 1.0
                taskId = self.tasks[task]['ID']
                value = self.DataSet["Tasks"][taskId][coordinates[i]]['value']
                if divisor==0 or not self.check_condition(self.task_list[task],(value-subtrahend)/divisor):
                    isSignal = 0
                    break
            if (mode and not isSignal) or (not mode and isSignal):
                return left + i*self.smallstep - self.frame/2
            elif mode:
                return left + len(coordinates)*self.smallstep - self.frame/2
            else:
                return left + self.frame/2
                
    def check_scenario_availability(self, task, tasks=None):
        if tasks == None:
            tasks = self.tasks
        # Check if all scenario parameters were calculated
        if task not in tasks or "scenario" not in tasks[task] or "frames" not in tasks[task] or not tasks[task]['frames']:
            print(task, tasks[task])
            raise ValueError(f"Scenario for {task} identification has not been set!")
        if any([subtask not in tasks[task]['scenario'] for subtask in tasks[task]['frames']]):
            raise ValueError(f"Not all {task} scenario parameters were calculated!")
        # Check if all scenario parameters are asigned with lists of values of the same length
        lengths = [len(tasks[task]['frames'][subtask]) for subtask in tasks[task]['frames']]
        if lengths[0] != sum(lengths)//len(lengths):
            raise ValueError(f"{task} scenario parameter datasets have different lengths!")
        return True
                     
    def meet_scenario(self, task, index, tasks=None, scenario_tasks=[], flg_perform_scenario_checking=False):
        if tasks == None:
            tasks = self.tasks
            
        def meet_conditions(i, task, subtask, tasks):
            settings = tasks[task]['scenario'][subtask]
            # Every locus contain data-point value in the format ['start-end',value]
            value = tasks[task]['frames'][subtask][i][1]
            # Sigma setting means that the value must be converted to (value - mean)/stdev
            if settings['mode'] == 'sigma':
                if 'stat' not in settings:
                    # Calculate mean and standard deviations for all values of the given parameter accross the whole genome
                    settings['stat'] = tools.getMeanAndStDev([float(ls[1]) for ls in tasks[task]['frames'][subtask]])
                mean,stdev = settings['stat']
                if not stdev:
                    return False
                value = (value - mean)/stdev
            if settings["condition"] == "bigger than" and value > settings["val1"]:
                return True
            elif settings["condition"] == "smaller than" and value < settings["val1"]:
                return True
            elif settings["condition"] == "between" and value >= settings["val1"] and value <= settings["val2"]:
                return True
            else:
                return False
            
        if flg_perform_scenario_checking:
            # Check scenario
            self.check_scenario_availability(task=task, tasks=tasks)
        if not scenario_tasks:
            scenario_tasks = list(tasks[task]['frames'].keys())
            
        # Check scenario conditions
        for subtask in scenario_tasks:
            if not meet_conditions(index, task, subtask, tasks):
                return False
        return True
            
    def identify_GI_Borders(self, task="MGE", tasks=None, flg_report=False):
        if tasks == None:
            tasks = self.tasks
        # Sort list by coordinates
        def parse_coordinates(coord):
            if "-" in coord:
                return list(map(int, coord.split("-")))
            elif ".." in coord:
                return list(map(int, coord.split("..")))
            else:
                return [int(coord)]
    
        def combine(ls, distance=0):
            # Convert format [['start-end', index], ...] to [['start-end', [index]], ...]
            ls = [[l[0], [l[1]]] for l in ls]
            if len(ls) < 2:
                return ls
        
            ls.sort(key=lambda x: parse_coordinates(x[0]))
        
            # Combine overlapping regions
            for i in range(len(ls) - 1, 0, -1):
                start1, end1 = parse_coordinates(ls[i][0])
                start2, end2 = parse_coordinates(ls[i - 1][0])
        
                if start1 - end2 <= distance:
                    # Merge regions
                    ls[i - 1] = [f"{start2}-{end1}", ls[i - 1][1] + ls[i][1]]
                    ls.pop(i)  # Remove the combined element        
            return ls
        
        # Format list [['start-end',[indices]],...) to [['start-end',{'task':[min,max,avr],...}],...)
        def replace_indices_to_values(gi, scenario_tasks, task="MGE"):  # gi = ['start-end',[indices]]; scenario_tasks - list of scenario tasks (parameters)
            values = dict(zip(scenario_tasks,[[] for i in range(len(scenario_tasks))]))
            for subtask in scenario_tasks:
                for i in gi[1]:
                    # Data-point values in the format ['start-end',value]
                    values[subtask].append(self.tasks[task]['frames'][subtask][i][1])
            for subtask in scenario_tasks:
                values[subtask] = tools.get_MinMaxAvr(values[subtask])
            return [gi[0],values]

        # Check scenario
        self.check_scenario_availability(task=task, tasks=tasks)
        
        # List of scenario tasks
        scenario_tasks = list(self.tasks[task]['frames'].keys())
        
        # Sort all data points in scenario value lists by coordinates
        for subtask in scenario_tasks:
            self.tasks[task]['frames'][subtask].sort(key=lambda ls: parse_coordinates(ls[0]))
            
        # Identify genomic loci which meet scenario criteria
        # Every locus contain data-point value in the format ['start-end',value]
        loci = self.tasks[task]['frames'][scenario_tasks[0]]
        GI_loci = combine([[loci[i][0],i] for i in range(len(loci)) if self.meet_scenario(task=task, index=i, scenario_tasks=scenario_tasks)], distance=self.bigstep)
        
        # Replace sliding window indices with max-min-avr values of task parameters
        GI_loci = [replace_indices_to_values(gi, scenario_tasks, task) for gi in GI_loci]
        return GI_loci
    
    # return list of ".fst" files in the home directory
    def get_fileList(self,home_folder):
        extensions = ('.FNA','.FAS','.FST','.FASTA',".GBK",".GBFF")
        l = os.listdir(home_folder)
        filelist = []
        for fname in l:
            if os.path.isfile(os.path.join(home_folder,fname)):
                if fname[fname.rfind("."):].upper() in extensions:
                    filelist.append(os.path.join(home_folder,fname))
                else:
                    continue
            elif os.path.isdir(os.path.join(home_folder,fname)):
                local_files = self.get_fileList(os.path.join(home_folder,fname))
                for item in local_files:
                    filelist.append(item)
            else:
                pass
        return filelist
    
    # Set new dataset for given sequence
    def set_dataset(self, seqname, seqLen, step):
        self.DataSet['Sequence name'] = seqname
        self.DataSet['Total sequence length'] = seqLen
        self.DataSet['Locus length'] = seqLen
        self.DataSet['Left border'] = 0
        self.DataSet['Frame'] = self.frame
        self.DataSet['Step'] = step
        for currTask in self.tasks:
            if currTask in self.DataSet["Tasks"]:
                continue
            # Prepare dataset for the tasks with scanario, like MGE
            if self.tasks[currTask]['scenario']:
                scenario_tasks = self.tasks[currTask]['scenario']
                for task in scenario_tasks:
                    for subtask in [s.strip() for s in (task.split("-") if task.find("-") > -1 else task.split("/"))]:
                        if subtask in self.DataSet["Tasks"]:
                            continue
                        self.DataSet["Tasks"][subtask] = {}
            else:
                self.DataSet["Tasks"][currTask] = {}
                for subtask in (self.tasks[currTask]["modifiers"]["subtrahend"],self.tasks[currTask]["modifiers"]['divisor']):
                    if not subtask or subtask in self.DataSet["Tasks"]:
                        continue
                    self.DataSet["Tasks"][subtask] = {}
    
    # Calculate basic values like GC-content, which do not require k-mer patterns
    def calculate_basic_parameters(self, currTask, key, locus):
        task = currTask
        if task not in self.DataSet["Tasks"] or not isinstance(self.DataSet["Tasks"][task],dict):
            raise ValueError(f"Programrun error!\nTask {task} was not added to DataSet.")
        if task == "GC":
            self.DataSet["Tasks"][task][key] = tools.calculate_GC(locus)
        elif task == "GCS":
            self.DataSet["Tasks"][task][key] = tools.calculate_GC_skew(locus)
        elif task == "AT":
            self.DataSet["Tasks"][task][key] = tools.calculate_AT(locus)
        elif task == "ATS":
            self.DataSet["Tasks"][task][key] = tools.calculate_AT_skew(locus)
    
    # Calculate scenario subtask paraemeters
    def calculate_scenario_subtask_parameters(self, currTask, key, locus, seq, tasks=None, genes=None):
        if tasks == None:
            tasks = self.tasks
        for subtask in tasks[currTask]['scenario']:
            if self.is_complex_task(subtask):
                # Process 'subtr' and 'divisor'
                for task_title in [s.strip() for s in (subtask.split("-") if subtask.find("-") > -1 else subtask.split("/"))]:
                    pt_name, taskID, pt_type, norm, wlength = self.oValidator.parse_task_title(task_title)
                    # Check if this value was already calculated
                    if task_title in self.DataSet["Tasks"] and key in self.DataSet["Tasks"][task_title]:
                        continue
                    self.evaluate_pattern(subtask=task_title, seq=seq, locus=locus, taskID=taskID, key=key,
                        pt_name=pt_name, wlength=wlength, norm=norm, pt_type=pt_type, genes=genes)
            else:
                task = tasks[currTask]['scenario'][subtask]
                self.evaluate_pattern(subtask=subtask, seq=seq, locus=locus, taskID=task['task'], key=key,
                    pt_name=task['k-pattern'], wlength=task['wlength'], norm=task['norm'], pt_type=task['type'], genes=genes)
            
    # Calculate k-mer pattern parameters
    def calculate_pattern_parameters(self, currTask, key, locus):        
        # Process 'subtr' and 'divisor'
        if self.is_complex_task(subtask):
            for task_title in [s.strip for s in (subtask.split("-") if subtask.find("-") > -1 else subtask.split("/"))]:
                pt_name, taskID, pt_type, norm, wlength = self.oValidator.parse_task_title(task_title)
                # Check if this value was already calculated
                if task_title in self.DataSet["Tasks"] and key in self.DataSet["Tasks"][task_title]:
                    continue
                self.evaluate_pattern(subtask=task_title, seq=seq, locus=locus, taskID=taskID, key=key,
                    pt_name=pt_name, wlength=wlength, norm=norm, pt_type=pt_type, genes=genes)
        else:
            pt_name, taskID, pt_type, norm, wlength = self.oValidator.parse_task_title(task_title)
            self.evaluate_pattern(subtask=currTask, seq=seq, locus=locus, taskID=task['task'], key=key,
                pt_name=task['k-pattern'], wlength=task['wlength'], norm=task['norm'], pt_type=task['type'], genes=genes)
        
    # Create k-mer patterns and record respective parameters to sliding windows
    def evaluate_pattern(self, subtask, seq, locus, taskID, key, pt_name, wlength, norm, pt_type, genes=None):
        #task_title = f"MGE:{subtask}"
        task_title = subtask
        # create reference pattern for a new task
        if pt_name not in self.StandardPatterns:
            self.StandardPatterns[pt_name] = self.getPattern(seq, wlength, norm, pt_type, genes)
        # create local pattern for a frame
        # Check if the pattern was alredy calculated for another task, otherwise create a new local pattern
        oPattern = None
        if (task_title in self.DataSet and key in self.DataSet["Tasks"][task_title] and self.DataSet["Tasks"][task_title][key] and 
                pt_name in self.DataSet["Tasks"][task_title][key]['oup']):
            oPattern = self.DataSet["Tasks"][task_title][key]['oup'][pt_name]
        value,oPattern = self.get_value(locus=locus,
                               mode=self.oValidator.tasks[taskID],
                               normalization=norm,
                               wlength=wlength,
                               pattern_type=pt_type,
                               oStdPattern=self.StandardPatterns[pt_name],
                               oCurrPattern=oPattern)
        # Set dataset
        start,stop = [int(v) for v in key.split("-")]
        DataSet = {'value':value,
                   'start':start,
                   'stop':stop,
                   'oup':{pt_name:oPattern}}
        if task_title not in self.DataSet["Tasks"]:
            self.DataSet["Tasks"][task_title] = {}
        self.DataSet["Tasks"][task_title][key] = {}
        self.DataSet["Tasks"][task_title][key].update(DataSet)

    # return requested data table for the sequence in the given file
    def calculate(self, subseq, full_length, name, step=None, frame=None, shift=0, tasks= None, genes=None, progress_bar=None, flg_update_references=True, flg_renew_dataset=False):
        if tasks == None:
            tasks = self.tasks
            
        if step == None:
            step = self.bigstep
            
        if frame == None:
            frame = self.frame

        # Setting dataset
        if flg_renew_dataset:
            self.set_dataset(name, len(subseq), self.bigstep)
        
        if flg_update_references:
            self.StandardPatterns = {}
            
        # maximal number of digits in coordinates
        maxnumlen = len(str(full_length))
        start = 0
        stop = frame-1
        
        while stop < len(subseq) and stop + shift < full_length:
            if self.echo:
                print([start,stop])
            if stop > len(subseq):
                locus = subseq[start:]
                stop = len(subseq)
            else:
                locus = subseq[start:stop]
                          
            key = (maxnumlen - len(str(start+shift)))*" " + str(start+shift) + '-' + str(stop+shift)
            
            # Calculate values of tasks for the selected locus
            for currTask in tasks:
                if currTask == "MOD":
                    self.DataSet["Tasks"][currTask][key] = [start+shift,stop+shift]
                # if Tasks were already calculated in another process, they can be taken from the object self.oCompletedTasks
                if not self.oCompletedTasks:
                    if currTask != "MGE" and not tasks[currTask]['k-pattern']:
                        # Calculate basic values like GC-content, which do not require k-mer patterns
                        self.calculate_basic_parameters(currTask=currTask, key=key, locus=locus)
                    else:
                        # Calculate GIs and parameters, which require k-mer patterns
                        if currTask == "MGE":
                            self.calculate_scenario_subtask_parameters(tasks=tasks, currTask=currTask, key=key, locus=locus, seq=subseq, genes=genes)
                        else:
                            self.calculate_pattern_parameters(currTask=currTask, key=key, locus=locus, seq=subseq, genes=genes)
            
            # Show progress bar
            if progress_bar:
                try:
                    progress_bar(shift + start + 1)
                except:
                    pass

            start += step
            stop += step
            if self.echo:
                print("has been done in " + str(time.perf_counter() - TaskStartTime) + " sec.")
                print()
    
    # Return requested values for provided reference and local k-mer patterns
    def get_value(self, locus, mode, normalization, wlength, pattern_type, oStdPattern, oCurrPattern=None):
        # create local pattern
        if not oCurrPattern:
            oCurrPattern = self.getPattern(locus,wlength,normalization,pattern_type)
        oSubtrPattern = None
        if mode == 'GC-content':
            return '%s' % oCurrPattern.getPercentage("GC")
        elif mode == 'G/C-skew':
            return '%s' % oCurrPattern.getGCskew()
        elif mode == 'A/T-skew':
            return '%s' % oCurrPattern.getATskew()
        elif mode in ('Generalized distance',
                    'Generalized pattern skew',
                    'Generalized relative pattern skew',
                    'Generalized variance',
                    'Generalized relative variance') and normalization:
            oNormalizationTable = oStdPattern.getNormalizationTable()
            oCurrPattern.setNormalizationTable(oNormalizationTable)
        else:
            pass
            
        if mode in ('Distance', 'Generalized distance'):
            return oCurrPattern-oStdPattern,oCurrPattern
        
        elif mode in ('Pattern skew',
                      'Relative pattern skew',
                      'Generalized pattern skew',
                      'Generalized relative pattern skew'):
            return '%s' % oCurrPattern.getPS(),oCurrPattern
        
        elif mode in ('Variance',
                      'Generalized variance',
                      'Relative variance',
                      'Generalized relative variance'):
            return '%s' % oCurrPattern.getOUV(),oCurrPattern
        
        else:
            print('Error mode ' + mode)
            return None
    
    # Create requested k-mer pattern
    def getPattern(self,seq,wlength,norm,ptype,genes=None):
        if self.flg_coding_sequences:
            seq = self.getCodingSequence(seq,genes)
        oPattern = Pattern(wlength)
        oPattern.setPattern(seq, norm, ptype)
        return oPattern
        
    # Return the completed tasks
    def getCompletedTasks(self):
        return self.tasks
    
    # Return coding sequences of a replicon according to gene locations provided
    def getCodingSequence(self,seq,genes):
        if not self.flg_contrasting or not genes:
            return seq
        cseq = ""
        for gene in genes:
            description = ""
            for key in ['name','description','remark']:
                description += genes[gene].qualifiers['name'][0] if key in genes[gene].qualifiers else ""
            if description.find('hypothetical') > -1:
                continue
            coords = gene.split(" | ")[-1][1:-1]
            lb,rb = [int(v) for v in (coords.split("..") if coords.find("..") > -1 else coords.split("-"))]
            cseq += seq[lb:rb]
        if len(cseq) < 10000:
            return seq
        return cseq
    
    def updateStdPatterns(self,seq):
        for key in self.StandardPatterns:
            ptype = key.split(":")[0]
            norm = int(ptype[1])
            wlength = int(ptype[3])
            ptype = ptype[0]
            self.StandardPatterns[key] = self.getPattern(seq,wlength,norm,ptype)
    
    def _generic_name(self,seqname,acc):
        cgr = seqname.rfind(", complete")
        if cgr > -1:
            seqname = seqname[:cgr]
        seqname = seqname.replace(" ","_")
        return "%s_[%s]" % (seqname,acc)
        
    def getCurrTime(self):
        CurrTime = time.localtime()
        CurrTime = (str(CurrTime[3]) + ":" + str(CurrTime[4]) + ":" + str(CurrTime[5]) + " " +
                str(CurrTime[2]) + "." + str(CurrTime[1]) + "." + str(CurrTime[0]))
        return CurrTime
        
    def is_complex_task(self, task):
        task, symbol, operand = self.oValidator.parse_task(task)
        if operand == None:
            return False
        return True

    def save_svg(self,seqname):
        oIO = seq_io.IO()
        oIO.save(self.oSVG.get_svg(),os.path.join(self.output_path,seqname+".svg"))

#############################################################################################################
class Validator:
    def __init__(self):
        pass
        # ATTRIBUTES
        self.MinimalLength = {
            7:295000,
            6:74000,
            5:18500,
            4:4600,
            3:1200,
            2:300
            }
        self.tasks = {'GC':'GC-content',
                      'GCS':'G/C-skew',
                      'ATS':'A/T-skew',
                      'D':'Distance',
                      'GD':'Generalized distance',
                      'PS':'Pattern skew',
                      'RPS':'Relative pattern skew',
                      'GPS':'Generalized pattern skew',
                      'GRPS':'Generalized relative pattern skew',
                      'V':'Variance',
                      'GV':'Generalized variance',
                      'RV':'Relative variance',
                      'GRV':'Generalized relative variance'}

    # METHODS
    def select_optimal_k(self, seqlength):
        length_borders = sorted(list(self.MinimalLength.items()))
        for k, length in length_borders:
            if seqlength < length:
                return k - 1
        return length_borders[-1][0]
        
    # Parse pattern names like n0_4mer:D
    def parse_task_title(self,task_title):
        taskID = ""
        if task_title.find(":"):
            task_title,taskID = task_title.split(":")
        else:
            return [task_title]
        pt_type, norm, wlength = self.parse_pattern(task_title)
        return task_title, taskID, pt_type, norm, wlength
        
    def parse_pattern(self,pattern):
        pt_type = pattern[0]
        norm = int(pattern[1:pattern.find("_")])
        wlength = int(pattern[pattern.find("_")+1:-3])
        return pt_type, norm, wlength
        
    def parse_gene_string(self, gene, flg_location=False):
        try:
            elements = [s.strip() for s in gene.split("|")]
            if len(elements) == 4:
                locus_tag, gene_name, product, location = elements
            elif len(elements) == 3:
                gene_name = "."
                locus_tag, product, location = elements
            start, end = self.parse_location(location)
            if flg_location:
                return start, end
            return locus_tag, gene_name, product, start, end
        except:
            raise ValueError(f"Wrong gene string format {gene}.\nMust be: tag|name|product|[start..end].\n")
            
    def gene_strings_to_dicts(self, genes):
        #gene_list = list(map(lambda gene: self.parse_gene_string(gene), genes))
        gene_list = [self.parse_gene_string(gene) for gene in genes]
        return [dict(zip(['locus_tag','name','product','start','end'],item)) for item in gene_list]
        
    # Parse complex tasks with operands
    def parse_task(self, task):
        subject = task
        operand = operator = None
        for symbol in ("*", "/", "+", "-"):
            p = task.find(symbol)
            if p > -1:
                subject = task[:p].strip()
                operand = task[p+1:].strip()
                # Check if operand is a number
                try:
                    operand = float(operand)
                except:
                    pass
                # Check if subject is a number
                try:
                    subject = float(subject)
                except:
                    pass
                operator = task[p]
                break
        return subject, operator, operand
        
    # Parse coords in formats 'start-end', or 'start..end', or 'sart:end' to [int(start),int(end)]
    def parse_location(self, coords):
        coords = coords.strip().replace("[","").replace("]","")
        separators = ("-", "..", ":")
        for symbol in separators:
            if coords.find(symbol) > -1:
                coords = coords.split(symbol)
                if len(coords) != 2:
                    raise ValueError(f"Separator symbold {symbol} was found more then 1 time!")
                try:
                    coords = [int(s) for s in coords]
                except:
                    raise ValueError("Coordinates must be 'start' and 'end' integers!")
                return coords
        raise ValueError(f"Separators {'; '.join(separators)} were not found in the string {coords}!")
        
    def validateTasks(self,task_list,frame,keep_initial_settings=False):
        tasks = {}
        # check tasks
        LengthLimits = []
        for currTask in task_list:
            SuperTask = str.split(currTask,'-')
            subtrahend = None
            divisor = None
            if len(SuperTask) > 2:
                showError("Wrong task: " + currTask)
                return None
            elif len(SuperTask) == 2:
                subtrahend = SuperTask[1]
                divisor = None
            else:
                SuperTask = str.split(SuperTask[0],'/')
                if len(SuperTask) > 2:
                    showError("Wrong task: " + currTask)
                    return None
                elif len(SuperTask) == 2:
                    divisor = SuperTask[1]
                    subtrahend = None
                else:
                    pass
            TaskId = SuperTask[0]
            
            Values = str.split(TaskId,":")
            if len(Values) == 1:
                Values = [None,Values[0]]
            elif len(Values) == 2:
                pass
            else:
                showError("Wrong task: '" + str(TaskId) + "'.\nMust be like 'n0_4mer:D'")
                return None
            
            Task = Values[1]
            if Values[0]:
                Type = str.split(str(Values[0]),"_")
                if len(Type) != 2:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                
                Norm = str(Type[0])
                if len(Norm) != 2 or Norm[0] != "n":
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                try:
                    pattern_type = Norm[0]
                    Norm = int(Norm[1])
                except:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                
                Wlength = str(Type[1])
                if len(Wlength) != 4 or Wlength[1:] != "mer":
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                try:
                    Wlength = int(Wlength[0])
                except:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
            else:
                Task = Values[1]
                Norm = 0
                Wlength = 4
                pattern_type = "n"
            
            # add task to the task list
            tasks[currTask] = {'task':Task,
                            'norm':Norm,
                            'type':pattern_type,
                            'wlength':Wlength,
                            'subtr':subtrahend,
                            'divisor':divisor,
                            'ID':TaskId,
                            'k-pattern':f'n{Norm}_{Wlength}mer'}
            subtrahend = None
            divisor = None
            LengthLimits.append(self.get_MinimalLength(Wlength))

            # check frame sizes
            if frame < min(LengthLimits):
                result = showError("The sliding window size must be at least " + str(min(LengthLimits)) + " bp.")
                return None
            else:
                pass
            if keep_initial_settings:
                tasks[currTask].update(task_list[currTask])
        return tasks

    # ACCESSIONS
    def get_MinimalLength(self, n):
        return self.MinimalLength[n]

    def get_TaskCategory(self, task):
        return self.tasks[task]

###################################################################################################################
class Pattern:
    def __init__(self,wlength):
        # ATTRIBUTES
            # Number of words in collection
        self.TotalWordNumber = 0
            # Word length
        self.wlength = wlength
            # Collection of words: word name like {'AAAA' : [real number, expected number, word spreading]}
            # diviation is measured in values of expected standard diviation calculated by regression
            # equation for given word length and sequence length
        self.words = {}
            # AGCT-box
        self.boxAGCT = {'A':0,'G':0,'C':0,'T':0}
            # Normalization type
        self.normalization = None
            # object NormalizationTable
        self.oNormalizationTable = None
            # reverse complementation
        self.RevComplTable = {'A':'T','T':'A','G':'C','C':'G'}
            # pattern type
        self.PatternType = None
            # word statistics
        self.flg_addWordStatistics = 0
        
        self.initiate()

    # METHODS
    # initiate words
    def initiate(self):
        # Fill in the generic pattern with words
        pmlib = PatternMethod()
        WordList = pmlib.getWordList(self.wlength)
        for i in range(len(WordList)):
            # word data [real,expected,regularity]
            self.words[WordList[i]] = [0,0,None]
            

    # Methods
    def setPattern(self, strSeq, normalization, pattern_type):
        if not len(strSeq):
            print(861)
        self.PatternType = pattern_type
        self.normalization = normalization
        if self.PatternType in ("d","s"):
            flg_addWordStatistics = 1
        else:
            flg_addWordStatistics = 0
        if normalization == 0:
            normalization = 1
        
        # set mononucleotite matrix
        strSeq = str.upper(strSeq)
        self.setMatrix(strSeq)
        # Calculate expected and real numbers of words: [real number, expected number]
        self.setWordStatistics(strSeq, flg_addWordStatistics)
        if self.PatternType in ("d","s"):
            self.processStatData()
        # set normalization table
        if self.oNormalizationTable == None:
            WordFrequencies = {}
            if normalization == 1:
                WordFrequencies.update(self.boxAGCT)
            else:
                oNormPattern = Pattern(normalization)
                oNormPattern.setPattern(strSeq)
                for word in list(oNormPattern.words.keys()):
                    WordFrequencies[word] = oNormPattern.getWordFrequency(word)
                del oNormPattern
            oNormalizationTable = NormalizationTable(WordFrequencies,normalization)
            self.setNormalizationTable(oNormalizationTable)
            self.setExpectation()
        else:
            self.oNormalizationTable = None

    def setExpectation(self):
        wordlist = list(self.words.keys())
        pVal = float(len(wordlist)/self.TotalWordNumber)
        for word in wordlist:
            if self.oNormalizationTable:
                pVal = self.oNormalizationTable.getWordLikelihood(word)
            self.words[word][1] = pVal*self.TotalWordNumber
    
    def setMatrix(self, strSeq):
        # Calculate numbers of A, G, C and T in the sequence and assign values to the 4-cell table boxAGCT
        letters = list(self.boxAGCT.keys())
        for l in letters:
            self.boxAGCT[l] = float(str.count(strSeq,l))/len(strSeq)
        return 1

    def setWordStatistics(self, strSeq, flg_addWordStatistics, option=None):
        keys = list(self.words.keys())
        # set real numbers of words and total number of words depending on allowing of word overlapping
        badWords = ''
        # total number of words
        self.TotalWordNumber = len(strSeq)-self.wlength+1
        # real number of words
        WordPosition = {}
        # number of lots for intrinsic PS
        lNum = self.getLotLength(self.TotalWordNumber)
        lot = self.TotalWordNumber/lNum
        k = 1
        for i in range(self.TotalWordNumber):
            if k and i >= k*lot:
                k = k + 1
                if k > lNum:
                    k = None
            word = strSeq[i:i+self.wlength]
            try:
                self.words[word][0] += 1
                if flg_addWordStatistics:
                    # set statistics for intrinsic pattern skew
                    if self.words[word][2]:
                        distance = i - WordPosition[word]
                        WordPosition[word] = i
                        self.words[word][2]['Sum 1'] = self.words[word][2]['Sum 1'] + distance
                        self.words[word][2]['Sum 2'] = self.words[word][2]['Sum 2'] + distance*distance
                        if k:
                            self.words[word][2]['PartialFrequ'][k-1] = self.words[word][2]['PartialFrequ'][k-1] + 1
                    else:
                        WordPosition[word] = i
                        self.words[word][2] = {'Sum 1':0,
                                               'Sum 2':0,
                                               'Start':i,
                                               'PartialFrequ':[]}
                        for j in range(lNum):
                            self.words[word][2]['PartialFrequ'].append(0)
            except:
                if str.find(word,'*') == -1:
                    badWords = badWords + word + '; '
        return badWords

    def processStatData(self):
        tmp_words = []
        for word in self.words:
            try:
                RevComplWord = tools.getReverseComplement(word)
            except:
                tmp_words.append([word,None])
            try:
                lNum = len(self.words[word][2]['PartialFrequ'])
            except:
                tmp_words.append([word,None])
            if word == RevComplWord:
                tmp_words.append([word,0])

            words = list(self.words.keys())
            sum = 0
            
            # distance in ranks
            for i in range(lNum):
                data = []
                for w in words:
                    if self.words[w][2]:
                        data.append([self.words[w][2]['PartialFrequ'][i],w])
                    else:
                        data.append([0,w])
                data.sort(key=lambda ls: [ls[0],ls[1]], reverse=True)
                r = 1
                for j in range(len(data)):
                    if data[j][1] == word or data[j][1] == RevComplWord:
                        if r:
                            rank = j
                            r = r - 1
                        else:
                            sum = sum + abs(rank - j)
                            break
            tmp_words.append([word,float(sum)/float(lNum)])

        for item in tmp_words:
            if self.words[item[0]][2]:
                self.words[item[0]][2]['PartialFrequ'] = item[1]
            else:
                self.words[item[0]][2] = {'Sum 1':0,
                                               'Sum 2':0,
                                               'Start':0,
                                               'PartialFrequ':0}

    def getLotLength(self, slength):
        MinimalSequencesLength = {
            7:295000,
            6:74000,
            5:18500,
            4:4600,
            3:1200,
            2:300
            }
        minLength = MinimalSequencesLength[self.wlength]
        lLen = int(2.0*math.sqrt(float(slength+self.wlength)/float(minLength)))
        if lLen:
            return lLen
        else:
            return 1

    def getWordLength(self):
        return self.wlength
    
    def getPatternType(self):
        return self.PatternType
    
    def getPatternName(self):
        return self.getPatternType() + str(self.getNormalization()) + "_" + str(self.getWordLength()) + "mer"
    
    def getNormalization(self):
        return self.normalization

    def getSeqLength(self):
        return self.TotalWordNumber + self.wlength
    
    def getWordRanks(self):
        oWordList = self.getWordList()
        return oWordList.getWordList()

    # take pattern, standard list of words and mode, return ['word','pattern type'],
    # if take None, return current WordList
    def getWordList(self):
        oWordList = WordList(self.getWordLength())
        #normstdiv = 0.14+(4**self.getWordLength())/self.getSeqLength()
        normal_expect = self.getMeanWordNumber()
        for item in oWordList.wordList:
            value = self.getWordAbundance(item[0])
            # set word frequencies
            item[3] = self.getWordFrequency(item[0])
            if self.PatternType == 'n' and self.normalization == 0:
                item[2] = self.getDeviation(value,normal_expect,normal_expect)
            elif self.PatternType == 's':
                # item[2] = self.getWordSkewness(item[0])/stdivIntrinsicPS
                item[2] = self.getWordSkewness(item[0])/10.0
                if item[2] == None:
                    return None
            elif self.PatternType == 'd':
                IPS = self.getWordSkewness(item[0])/10.0
                if IPS == None:
                    return None
                if self.normalization:
                    # normalized patterns
                    # expected number of word
                    expect = self.getWordExpectation(item[0])
                    # correction for the letter inequivalence
                    if expect == 0:
                        expect = 1
                    if normstdiv == 0:
                        normstdiv = 1
                    item[2] = [self.getDeviation(value,expect,normal_expect),IPS]
                else:
                    item[2] = [self.getDeviation(value,normal_expect,normal_expect),IPS]
            else:
                # normalized patterns
                # expected number of word
                expect = self.getWordExpectation(item[0])
                # correction for the letter inequivalence
                if expect == 0:
                    expect = 1
                item[2] = self.getDeviation(value,expect,normal_expect)
        return oWordList
    
    def getDeviation(self,value,expect,stdev,flg_lognormal=0):
        if flg_lognormal:
            if value == 0:
                value = 1
            numerator = math.log((value**2)*math.sqrt(stdev**2 + expect**2)/(expect**2)/math.sqrt(stdev**2 + value**2))
            denumerator = math.sqrt(math.log((stdev/expect)**2 + 1))
            return 6.0*numerator/denumerator
        else:
            normstdiv = 0.14+(4**self.getWordLength())/self.getSeqLength()
            return float((value - expect)/stdev)/normstdiv

    def getNormalizationTable(self):
        return self.oNormalizationTable
    
    def getNormalizationTableParametr(self):
        if self.oNormalizationTable:
            return self.oNormalizationTable.getWordLength()
        else:
            return None

    def setNormalizationTable(self, oNormalizationTable):
        del self.oNormalizationTable
        self.oNormalizationTable = oNormalizationTable
        self.setExpectation()

    def getWordAbundance(self, word):
        return self.words[word][0]

    def getWordFrequency(self, word):
        return float(self.words[word][0])/float(self.TotalWordNumber)

    def getWordExpectation(self, word):
        if self.normalization == 0:
            return self.getMeanWordNumber()
        else:
            return self.words[word][1]

    def getWordSkewness(self, word):
        try:
            return float(self.words[word][2]['PartialFrequ'])
        except:
            return None

    def getMeanWordNumber(self):
        return self.TotalWordNumber/len(list(self.words.keys()))
    
    def getOUV(self):
        oWordList1 = self.getWordList()
        wordList = oWordList1.getWordList(self.PatternType)[0]
        sum = 0
        sum_sq = 0
        N = len(wordList)
        for i in range(N):
            sum = sum + wordList[i][2]
            sum_sq = sum_sq + wordList[i][2] * wordList[i][2]
        val = (sum_sq - sum * sum/N)/(N - 1)
        #return math.sqrt(val)
        return val
    
    def getPS(self):
        # set variables
        oWordList1 = self.getWordList()
        wl1 = oWordList1.getWordList(self.PatternType)
        wl2 = self.getReverseComplement(wl1)
        if self.PatternType == "d":
            wl1 = [wl1[0],]
            wl2 = [wl2[0],]
        D = self.getDistance(wl1,wl2,self.getMinimalDistance(),0,0)[0]
        return D
    
    def getPercentage(self,symbol):
        if len(symbol)==1:
            val = self.getNucleotideFrequency(symbol)
        elif symbol == "GC":
            val = self.getNucleotideFrequency("G") + self.getNucleotideFrequency("C")
        elif symbol == "AT":
            val = self.getNucleotideFrequency("A") + self.getNucleotideFrequency("T")
        else:
            val = 0
        return val*100
    
    def getGCskew(self):
        G = self.getNucleotideFrequency("G")
        C = self.getNucleotideFrequency("C")
        return (float(G-C)/float(G+C))
        
    def getATskew(self):
        A = self.getNucleotideFrequency("A")
        T = self.getNucleotideFrequency("T")
        return (float(A-T)/float(A+T))

    def getMinimalDistance(self):
        if self.wlength%2 == 0:
            return (4**self.wlength - 2**self.wlength)/2.0
        else:
            return (4**self.wlength)/2.0

    def getNucleotideFrequency(self,letter):
        return self.boxAGCT[letter]
    
    def getReverseComplement(self, wl):
        NewWordList = []
        for d in range(len(wl)):
            WordList = copy.deepcopy(wl[d])
            WordDic = dict(zip([item[0] for item in WordList],[item[1:] for item in WordList]))
            RevCompWordDic = dict(zip(list(WordDic.keys()),[[v for v in WordDic[tools.reverse_complement(key)]] for key in list(WordDic.keys())]))
            # Restor order numbers
            for key in list(RevCompWordDic.keys()):
                RevCompWordDic[key][0] = WordDic[key][0]
            RevCompWordList = [[key] + RevCompWordDic[key] for key in list(RevCompWordDic.keys())]
            NewWordList.append(RevCompWordList)
        return NewWordList

    # return a copy of the current pattern object
    def copy(self):
        oPattern = Pattern(self.wlength)
        oPattern.TotalWordNumber = self.TotalWordNumber
        oPattern.PatternType = self.PatternType
        oPattern.words = self.words
        oPattern.boxAGCT = self.boxAGCT
        oPattern.normalization = self.normalization
        oPattern.flg_addWordStatistics = self.flg_addWordStatistics
        try:
            oPattern.setNormalizationTable(self.getNormalizationTable())
        except:
            oPattern.setNormalizationTable(None)
        return oPattern
    
    def convert(self, pattern_type=None, normalization=None):
        oPattern = self.copy()
        if pattern_type == None:
            pattern_type = self.getPatternType()
        if normalization == None:
            normalization = self.getNormalization()
        if pattern_type == self.getPatternType():
            if normalization == self.getNormalization():
                return oPattern
            if normalization == 0 or normalization == self.getNormalizationTableParametr():
                oPattern.normalization = normalization
            else:
                return None
            oPattern.setExpectation()
            return oPattern
        elif pattern_type == "n" or oPattern.flg_addWordStatistics:
            oPattern.PatternType = pattern_type
        else:
            return None
        
    def compare(self,other):
        if type(other) != type(Pattern(2)):
            return None
        if (self.normalization != other.normalization) or (self.PatternType != other.PatternType):
            return None
        else:
            return 1
    
    # calculate distance between two pattern lists 
    def getDistance(self,first_list,second_list,minDist=0,flg_BestHit=0,flg_Normolize=0):
        wl1 = first_list
        dm = min((len(first_list),len(second_list)))
        total_sum = 0
        wlength = len(first_list[0][0][0])
        if dm > 1:
            WordNumber = int(4**wlength//2)
        else:
            WordNumber = int(4**wlength)
        maxDist = float(WordNumber*(WordNumber+1))
        StdDev = 0
        dist = maxDist
        StdDev = 0
        complement = 0
        if flg_BestHit:
            count = 3
        else:
            count = 0
        # calculate distance
        while count >= 0:
            s = []
            values = []
            for i in range(dm):
                if count==1:
                    wl1 = first_list
                    wl1[i].sort(key=lambda ls: ls[1])
                    wl2 = self.getReverseComplement(second_list)
                    wl2[i].sort(key=lambda ls: ls[1])
                elif count==2:
                    wl1 = self.getReverseComplement(first_list)
                    wl1[i].sort(key=lambda ls: ls[1])
                    wl2 = second_list
                    wl2[i].sort(key=lambda ls: ls[1])
                elif count==3:
                    wl1 = self.getReverseComplement(first_list)
                    wl1[i].sort(key=lambda ls: ls[1])
                    wl2 = self.getReverseComplement(second_list)
                    wl2[i].sort(key=lambda ls: ls[1])
                else:
                    wl1 = first_list
                    wl1[i].sort(key=lambda ls: ls[1])
                    wl2 = second_list
                    wl2[i].sort(key=lambda ls: ls[1])
                s.append(0)
                values.append([])
                for j in range(WordNumber):
                    val = abs(wl1[i][j][3]-wl2[i][j][3])/2
                    s[i] += val
                    values[i].append(val)
            
            subtotal = math.sqrt(sum([v*v for v in s])) if len(s) else 0
            if  subtotal < dist:
                dist = subtotal
                StdDev = self.getMeanAndStDev(values)[1]
            if count == 0:
                break
            count -= 1
        # normolize distance
        if flg_Normolize:
            OUV1 = getOUV(wl1[0])
            OUV2 = getOUV(wl2[0])
            m = -.25*math.sqrt(OUV1*OUV1+OUV2*OUV2)+4
        else:
            m = 1.0
        if dist < minDist:
            dist = minDist
        dist = 400.0*m*(dist-minDist)/(math.sqrt(dm)*maxDist-minDist)/dm
        return [dist,StdDev,complement]

    # Get list of values and return [mean value,standard deviation]
    def getMeanAndStDev(self,DataList,ind=2,binome=None):
        s = 0.0
        s_sq = 0.0
        s_cube = 0.0
        # if there is a number of lists, recalculate stdDev for all lists and return the average and the skew for the last list
        if type(DataList[0]) == type([1,2]):
            dm = len(DataList)
            N = len(DataList[0])
        else:
            dm = 1
            N = len(DataList)
            DataList = [DataList,]
        StdDev = 0
        for k in range(dm):
            for i in range(N):
                val = float(DataList[k][i])
                s = s + val
                if ind > 1:
                    s_sq = s_sq + val*val
                if ind > 2:
                    s_cube = s_cube + val*val*val
            N = float(N)
            mean = s/N
            if ind > 1:
                if binome:
                    if mean > 1:
                        bmean = mean/100.0
                        stdev = math.sqrt(bmean*(1-bmean))*100.0
                    else:
                        stdev = math.sqrt(mean*(1-mean))
                else:
                    var = s_sq/N - mean*mean
                    if var>0:
                        stdev = math.sqrt(var)
                    else:
                        stdev = 0
            else:
                stdev = None
            if ind > 2:
                if binome:
                    skew = 0
                else:
                    if stdev:
                        skew = (s_cube - 3.0*s*s_sq/N + 2.0*s*s*s/N/N)/stdev/stdev/stdev/N
                    else:
                        skew = 0
            else:
                skew = None
            if stdev:
                StdDev = math.sqrt(stdev*stdev + StdDev*StdDev)
            else:
                StdDev = 0
        return [mean,StdDev,skew]

    # Subtracting of patterns
    def __sub__(self,other):
        # set variables
        if self.getPatternName() != other.getPatternName():
            return None
        oWordList1 = self.getWordList()
        oWordList2 = other.getWordList()
        if not oWordList1 or not oWordList2:
            return None
        wl1 = oWordList1.getWordList(self.PatternType)
        wl2 = oWordList2.getWordList(other.PatternType)
        return self.getDistance(wl1,wl2,0,1,0)[0]

    # Summation of patterns
    def __add__(self,other):
        if (self.wlength != other.getWordLength()) or (self.PatternType != other.PatternType):
            return None
        oPattern = Pattern(self.wlength)
        oPattern.PatternType = self.PatternType
        for key in list(self.boxAGCT.keys()):
            oPattern.boxAGCT[key] = (self.boxAGCT[key]*self.getSeqLength() + other.boxAGCT[key]*other.getSeqLength())/(self.getSeqLength() + other.getSeqLength())
        oPattern.TotalWordNumber = self.TotalWordNumber + other.TotalWordNumber
        first_table = self.getNormalizationTable()
        second_table = other.getNormalizationTable()
        if first_table:
            first_words = first_table.getWords()
            second_words = second_table.getWords()
            for key in first_words:
                first_words[key] = (first_words[key]*self.getSeqLength() + second_words[key]*other.getSeqLength())/(self.getSeqLength() + other.getSeqLength())
            new_table = NormalizationTable(first_words)
            oPattern.setNormalizationTable(new_table)
        for key in list(self.words.keys()):
            oPattern.words[key][0] = self.words[key][0] + other.words[key][0]
        oPattern.setExpectation()
        return oPattern

###################################################################################################################
class WordList:
    def __init__(self,wlength):

        # ATTRIBUTES
        self.wlength = wlength
        self.wordList = []
        # Create library of pattern setting methods
        pmlib = PatternMethod()
        # Create list and dictionary of the words
        for i in range(4**self.wlength):
            word = pmlib.getNextWord(self.wlength)
            CurrWord = ['',0,0,None]
            CurrWord[0] = word
            CurrWord[1] = i
            self.wordList.append(CurrWord)

    # METHODS

    def getRange(self):
        if type(self.wordList[0][2]) == type([0,]):
            return len(self.wordList[0][2])
        else:
            return 1

    # n is int - index of values in the list
    def getWordList(self,pattern_type="n"):
        wl = []
        wl1 = []
        wl2 = []
        for item in self.wordList:
            if pattern_type == "d":
                wl1.append([item[0],item[1],item[2][0],0])
                wl2.append([item[0],item[1],item[2][1],0])
            else:
                wl1.append([item[0],item[1],item[2],0])
        wl.append(wl1)
        if len(wl2):
            wl.append(wl2)
        for l in wl:
            self.assignRanks(l)                
        return wl

    def assignRanks(self, WordList):
        WordList.sort(key=lambda ls: [ls[2],ls[0]], reverse=True)
        for j in range(len(WordList)):
            WordList[j][3] = j
        WordList.sort(key=lambda ls: ls[1], reverse=True)

    def getLength(self):
        return len(self.wordList)

    def getVariance(self):
        sum = 0
        sum_sq = 0
        N = len(self.wordList)
        for i in range(N):
            sum = sum + self.wordList[i][2]
            sum_sq = sum_sq + self.wordList[i][2] * self.wordList[i][2]
        val = (sum_sq - sum * sum/N)/(N - 1)
        return val

##############################################################################################################
class PatternMethod:
    def __init__(self, method="BOXAGCT"):
            # Method
        self.CurrMethod = method
            # Current word
        self.Word = ""
            # The letter of the word in focus
        self.Cursor = None
            # State of the word processing by any function
        self.flg_State = 0
            # Dictionary of methods of oligomer plot arrangment
        self.Methods = {}
            # Dictionary of letter functions 
        self.LetterFunk = {}
        
        self.setDictionary()

    # Methods
        # Block of functions for method BOXAGCT
    def A_BOXAGCT(self):
        self.changeLetters("G","C",self.Cursor)
        self.flg_State = 0       

    def G_BOXAGCT(self):
        self.changeLetters("A","C",self.Cursor)
        if self.Cursor >= 0:
            self.Cursor = self.Cursor - 1
        else:
            self.flg_State = 0
            
    def C_BOXAGCT(self):
        self.changeLetters("T","A",self.Cursor)
        if self.Cursor >= 0:
            self.flg_State = 0
        else:
            self.Cursor = self.Cursor - 1

    def T_BOXAGCT(self):
        self.changeLetters("C","A",self.Cursor)
        self.Cursor = self.Cursor - 1
 
    def changeLetters(self, a, b, position):
        if position >= 0:
            self.Word = self.Word[:position] + a + self.Word[(position+1):]
        else:
            if position == -1:
                self.Word = self.Word[:position] + b + self.Word[:(position+1)]
            else:
                self.Word = self.Word[:position] + b + self.Word[(position+1):]

    def setDictionary(self):
        # initialisation of dictionary of functions for method BOXAGCT
        self.LetterFunk["A"] = self.A_BOXAGCT
        self.LetterFunk["G"] = self.G_BOXAGCT
        self.LetterFunk["C"] = self.C_BOXAGCT
        self.LetterFunk["T"] = self.T_BOXAGCT
        self.Methods["BOXAGCT"] = self.LetterFunk

    def getNextWord(self, wlength):
        if self.Word == '':
            self.Word = "A"*wlength
        else:
            # create new word
            self.Cursor = wlength - 1
            self.flg_State = 1
            while self.flg_State == 1:
                self.Methods[self.CurrMethod][self.Word[self.Cursor]]()
        return self.Word

    def getWordList(self, wlength):
        if wlength != 0:
            self.Word = ''
            wordList = []
            for i in range(4**wlength):
                wordList.append(self.getNextWord(wlength))
            return wordList

###################################################################################################################
# words is a dictionary {'word':frequency}
class NormalizationTable:
    def __init__(self, words, wlength):
        
        # ATTRIBUTES
        self.words = words
        self.wlength = wlength

    # METHODS

    def getWordLikelihood(self, word):
        frame = len(list(self.words.keys())[0])
        pVal = self.words[word[:frame]]
        subword = word[1:frame]
        subset = self.setSubset(subword)
        start = 1
        stop = frame+1
        while stop <= len(word):
            pVal = pVal * subset[word[start:stop]]
            start = start + 1
            subword = word[start:stop]
            subset = self.setSubset(subword)
            stop = stop + 1
        return pVal

    def setSubset(self, subword):
        subset = {}
        wordlist = list(self.words.keys())
        TotalLikelihood = 0
        for word in wordlist:
            if subword == word[:len(subword)]:
                subset[word] = self.words[word]
                TotalLikelihood = TotalLikelihood + subset[word]
        SubSet = {}
        SubSet.update(subset)
        keys = list(SubSet.keys())
        for key in keys:
            if TotalLikelihood:
                SubSet[key] = SubSet[key]/TotalLikelihood
            else:
                SubSet[key] = 0
        return SubSet
    
    def getWords(self):
        words = {}
        words.update(self.words)
        return words
    
    def setWords(self,words):
        self.words = {}
        self.words.update(words)

    def getWordLength(self):
        return self.wlength

########################################################################
# Command line interface
class Interface:
    def __init__(self,options=None, completed_tasks=None, default=False, motif_setting="", max_sites_for_verification=3000, graph_title=""):
        self.oValidator = InterfaceValidator()
        self.motif_setting = motif_setting
        self.IO = seq_io.IO()
        self.graph_title = graph_title
        self.output = None
        self.max_sites_for_verification = max_sites_for_verification
        self.oCompletedTasks = completed_tasks
        self.scenaria = {"MGE":
            {
            "n0_4mer:D":{"mode":"sigma","condition":"bigger than","val1":2.0,"val2":0.0,"main":True},
            "n1_4mer:GRV/n1_4mer:RV":{"mode":"absolute","condition":"bigger than","val1":2.0,"val2":0.0,"main":False},
            "n0_4mer:PS":{"mode":"absolute","condition":"smaller than","val1":55.0,"val2":0.0,"main":False},
            }}   
        self.task_list = []
        self.options = {"-c":"MGE",         # name of the scenario
                   "-g":"",                 # GBK file
                   "-u":"yes",              # use BLASTN to filter rrn clusters
                   "-n":"no",               # use BLASTN to search tRNA on the borders of GI
                   "-l":"8000",             # sliding window length
                   "-b":"2000",             # sliding window big step
                   "-m":"500",              # sliding window medium step
                   "-s":"100",              # sliding window small step
                   "-r":"1000000",          # reference sequence step
                   "-e":"Contrasting",      # refinement No | Contrasting | Iteration | Contrasting/Iteration 
                   "-i":"input",            # input folder
                   "-o":"output",           # output folder
                   "-f":"no",               # save GI sequemces: no/fasta/gbk/gbk+fasta
                   "-w":"",                 # motif setting
                   "-v":"Yes",              # save SVG file
                   "-p":75,                 # promoter length
                   "-std":"off"             # strand to display
                }
        
        # Available tasks to perform
        self.task_categories = [ 
            {"title":"MOD","info":{"description":"Modified nucleotides","scenario":None,
                "k-pattern":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'bars','line_color':'','fill_color':'black'}}},
            {"title":"MGE","info":{"description":"Genomic islands","scenario":self.options['-c'],
                "k-pattern":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'box','line_color':'','fill_color':'pink'}}},
            {"title":"GC","info":{"description":"GC-content","scenario":None,
                "k-pattern":False,"modifiers":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"AT","info":{"description":"AT-content","scenario":None,
                "k-pattern":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"GCS","info":{"description":"GC-skew","scenario":None,
                "k-pattern":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"ATS","info":{"description":"AT-skew","scenario":None,
                "k-pattern":False,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"D","info":{"description":"k-mer pattern deviation","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"GD","info":{"description":"k-mer generalized pattern deviation","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"PS","info":{"description":"k-mer pattern skew","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"GPS","info":{"description":"k-mer generalized pattern skew","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"V","info":{"description":"k-mer variance","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            {"title":"GV","info":{"description":"k-mer generalized variance","scenario":None,
                "k-pattern":True,"modifiers":{"subtrahend":0,"divisor":0},'frames':[],'completed':False,
                "settings":{'style':'line','line_color':'blue','fill_color':''}}},
            ]
        
        # Tasks titled in this list cannot be removed from the list by the users
        # Use self.set_task(title, flg_set_as_immutable=True) to set this list
        self.immutable_tasks = []

        # Setting arguments
        if not default:
            if options:
                self.options.update(options)
            else:
                self.main_menu()

    # Execute selected program
    # GBK_files - list of input files in genbank or fasta formats
    # sites - list of locations of modified nucleotides/motifs [[start,end],...],...]
    # Each input sequence file is provided with a list of modified sites, len(GBK_files) == len(sites)
    def execute(self, GBK_files=[], top_indend=150, modified_sites=[], motif_setting="", verified=False):
        if GBK_files:
            self.options['-g'] = GBK_files if isinstance(GBK_files,list) else [GBK_files]
        # Check if the list of input sequence files corresponds to the list of modified files
        if len(self.options['-g']) != len(modified_sites):
            raise ValueError(f"Length of input sequence files {len(self.options['-g'])} must equal the list of modified sites {len(modified_sites)}!")
            
        options = self.oValidator.validate(self.options,list(self.scenaria.keys()))
        if options:
            self.options = options
        else:
            raise ValueError(f"Option settings did not pass validation:\n{options}")
            
        # Motif settings for the report
        if motif_setting:
            self.motif_setting = motif_setting
        
        # Conver task list to task dictionary
        d_tasks = dict(zip([t['title'].upper() for t in self.task_list],[t['info'] for t in self.task_list]))
        oMainWin = Main(task_list=d_tasks, options=self.options, completed_tasks=self.oCompletedTasks, 
            max_sites_for_verification=self.max_sites_for_verification, motif_setting=self.motif_setting, graph_title=self.graph_title)
        # Pass GBK file location
        self.output = oMainWin.process(self.options['-g'], immutable_tasks=self.immutable_tasks, top_indend=top_indend, 
            modified_sites=modified_sites, flg_exclude_strand=self.options['-std'], verified=verified)
        self.oCompletedTasks = oMainWin.getCompletedTasks()
        return oMainWin.get_svg(), oMainWin.get_stat_report()
        
    # Add a task to the task list by the task's title
    def set_task(self, title, pattern=None, flg_set_as_immutable=False):    # pattern is like n0_4mer
        #selected_task = list(filter(lambda t: t['title'] == title, self.task_categories))
        selected_task = [category for category in self.task_categories if category['title'] == title]
        if not selected_task:
            tools.alert(f"Task '{title}' does not exists!")
            return
        task = selected_task[0]
        # Check if any scenario is associated with this task
        if task['info']['scenario']:
            if task['info']['scenario'] not in self.scenaria:
                tools.alert(f"Scenario {task['info']['scenario']} is not recognized!")
                return
            # Replace name of scenario with task setting
            task['info']['scenario'] = self.scenaria[task['info']['scenario']]
        self.task_list.append(task)
        if self.task_list[-1]['info']['k-pattern'] != None:
            if pattern == 'None':
                tools.alert(f"Task {title} requires setting a k-mer pattern, which was not provided!")
                del self.task_list[-1]
                return
        if flg_set_as_immutable:
            # This task cannot be removed from the list by the users
            self.immutable_tasks.append(title + (f":{pattern}" if pattern else ""))
            
    def getCompletedTasks(self):
        return self.oCompletedTasks
    
    # show command prompt interface
    def edit_tasklist(self):
        response = ""
        while response != "Q":
            print()
            print("Edit task list")
            print()
            print("  A   Add a new task to the list;")
            print("  R   Remove a task from the list;")
            print("  Q   Quit;")
            try:
                response = input("?")
            except:
                continue
            response = str.upper(response)
            if response == "Q":
                return
            added_tasks = [t['title'] + (':'+t['info']['k-pattern'] if t['info']['k-pattern'] else '') for t in self.task_list]
            if self.immutable_tasks:
                #added_tasks = list(filter(lambda t: t['title'] not in self.immutable_tasks, added_tasks))
                added_tasks = [task for task in added_tasks if task['title'] not in self.immutable_tasks]
            if response == "A":
                self.task_list += self.add_task(added_tasks)
                return
            elif response == "R":
                while True:
                    item_to_remove = self.remove_task(added_tasks)
                    if item_to_remove == None:
                        continue
                    if not item_to_remove:
                        return
                    del added_tasks[item_to_remove]
                    del self.task_list[item_to_remove]
            else:
                continue

    # Add a new task to the list
    def add_task(self, used_tasks=[]):
        added_tasks = []
        response = ""
        while response not in ("S","Q"):
            locked_tasks = used_tasks + [t['title'] for t in added_tasks]
            available_tasks = []
            print()
            print("Add tasks and press 'S' to save or 'Q' to exit without saving:")
            print()
            for task in self.task_categories:
                if task['title'] not in locked_tasks and task['title'] not in self.immutable_tasks:
                    print(f"\t{task['title']}:\t{task['info']['description']};")
                    available_tasks.append(task['title'])
            print()
            if added_tasks:
                print(f"Currently selected tasks: {', '.join([t['title'] + (':' + t['info']['k-pattern'] if t['info']['k-pattern'] else '') for t in added_tasks])}")
                print()
            response = input("?").upper()
            if response == "Q":
                return []
            elif response == "S":
                return added_tasks
            elif response in available_tasks:
                #added_tasks.append(copy.deepcopy(list(filter(lambda t: t['title'] == response, self.task_categories))[0]))
                added_tasks.append(copy.deepcopy([task for task in self.task_categories if t['title'] == response][0]))
                if added_tasks[-1]['info']['k-pattern'] == True:
                    added_tasks[-1]['info']['k-pattern'] = self.set_pattern_type()
                    new_task_title = f"{added_tasks[-1]['title']}:{added_tasks[-1]['info']['k-pattern']}"
                    if len(added_tasks) > 1 and new_task_title in added_tasks[:-1]:
                        tools.alert(f"Task {new_task_title} has been added already!")
                        del added_tasks[-1]
            else:
                continue
    
    # Select pattern type
    def set_pattern_type(self):
        while True:
            print()
            response = input("Enter pattern type in format n#_##mer\n" + 
                "where ## is k of the k-mer and # - is normalization k-mer;\n" +
                "# < ##, for instance: n0_4mer:____")
    
            # 1. Check the format
            pattern = re.compile(r"^n(\d)_(\d{1,2})mer$")
            match = pattern.match(response)
        
            if not match:
                print("\nInvalid format! Input must follow the pattern 'n#_##mer', e.g., 'n0_4mer'.")
                continue
        
            # Extract # and ## from the input
            norm_kmer = int(match.group(1))
            kmer = int(match.group(2))
        
            # 2. Check that ## is in the range from 1 to 7
            if not (1 <= kmer <= 7):
                print("\nInvalid k-mer range! ## must be between 1 and 7.")
                continue
        
            # 3. Check that # is in the range from 0 to ##
            if not (0 <= norm_kmer < kmer):
                print(f"\nInvalid normalization k-mer! # must be less than ##, that is {kmer}.")
                continue
        
            return response 
    
    # Remove task from the list
    def remove_task(self, tasks):
        print()
        print("Remove task from the list")
        print()
        print("  0   escape;")
        for i in range(1,len(tasks),1):
            print(f"\t{i}\t{tasks[i]};")
        print()
        index = input("Select the task by its index: ")
        try:
            index = int(index)
        except:
            tools.alert("Select integer index of the task to remove!")
            return None
        if index >= len(tasks):
            tools.alert(f"Enter an integer from 0 to {len(tasks)-1}")
            return None
        return index
                        

    def save_scripts(self,path=None):
        if path==None:
            path = os.path.join("lib","scripts.txt")
        output = ""
        for script_name in self.scripts:
            output += script_name + "\n"
            for task in self.scripts[script_name]:
                output += task + "\n"
                for name in self.scripts[script_name][task]:
                    output += name + ":" + str(self.scripts[script_name][task][name]) + ","
                output = output[:-1] + "\n"
            output += "END\n"
        self.IO.save(output,path)

    def select_scenario(self):
        options = list(self.scripts.keys())
        while 5 > 0:
            print()
            print("Select scenario")
            print()
            print("  0   Quit")
            for i in range(len(options)):
                print("  " + str(i+1) + "   " + options[i])
            print()
            index = input("Select scenario by the index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(options):
                    self.options["-c"] = options[index-1]
                    self.task_list = {}
                    self.task_list.update(self.scripts[self.options["-c"]])
                    return
                elif index == 0:
                    return
                else:
                    print()
                    print("\tWrong scenario index")
                    print()
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(options)))
                print()
                continue
        
    def add_scenario(self):
        print()
        script_name = input("Enter name for the current scenario: ")
        print()
        if not script_name:
            return
        self.scripts[script_name] = {}
        self.scripts[script_name].update(self.task_list)
        self.save_scripts()
        
    def remove_scenario(self):
        options = list(self.scripts.keys())
        while 5 > 0:
            print()
            print("Remove scenario")
            print()
            print("  0   Quit")
            for i in range(len(options)):
                print("  " + str(i+1) + "   " + options[i])
            print()
            index = input("Select scenario by the index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(options):
                    del self.scripts[options[index-1]]
                elif index == 0:
                    self.save_scripts()
                    return
                else:
                    print()
                    print("\tWrong scenario index")
                    print()
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(options)))
                print()
                continue

    def select_task(self):
        index = None
        while index != 0:
            codes = []
            count = 1
            print()
            print("Select task category")
            print()
            print("  0   escape;")
            for task_code in list(self.task_categories.keys()):
                print("  " + str(count) + "   " + task_code + " (" + self.task_categories[task_code] + ");")
                codes.append(task_code)
                count += 1
            print()
            index = input("Select the task category by its index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(codes):
                    return codes[index-1]
                elif index != 0:
                    print()
                    print("\tWrong category index")
                    print()
                else:
                    break
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(codes)))
                print()
                continue
        return None

    def set_conditions(self,task,conditions):
        response = ""
        while response != "Q":
            print()
            print("Select condition for the task " + task)
            print()
            print("  G   Bigger than")
            print("  S   Smaller than")
            print("  B   Between")
            print("  M   " + conditions["mode"])
            print("  Q   quit;")
            try:
                response = input("?")
            except:
                continue
            response = str.upper(response)
            if response == "G":
                conditions["condition"] = "bigger than"
            elif response == "S":
                conditions["condition"] = "smaller than"
            elif response == "B":
                conditions["condition"] = "between"
            elif response == "M":
                if conditions["mode"] == "absolute":
                    conditions["mode"] = "sigma"
                elif conditions["mode"] == "sigma":
                    conditions["mode"] = "fraction"
                else:
                    conditions["mode"] = "absolute"
                continue
            elif response == "Q":
                return conditions
            else:
                continue
            
            question = conditions["condition"]
            while 5 > 0:
                print()
                try:
                    val = input(question + ": ")
                except:
                    print()
                    if question == "and":
                        print("Enter the upper limit of variation:")
                    else:
                        print("Enter the limit of variation:")
                    print()
                    continue
                try:
                    val = float(val)
                    if question == "and":
                        conditions["val2"] = val
                    else:
                        conditions["val1"] = val
                except:
                    print()
                    print("\tEnter a float point number")
                    print()
                    continue
                if conditions["condition"] != "between" or question == "and":
                    break
                question = "and"
            return conditions

    def main_menu(self):
        # Terminal interface
        response = ''
        while response != "Q":
            self.print_menu()
            try:
                response = input("?")
                response = str.upper(response)
            except:
                continue
            if response in ("Y","Q"):
                return self.options
            if response == "T":
                self.edit_tasklist()
            else:
                continue

    def print_menu(self):
            print("Genome pattern bar 2024/26/11")
            print()
            print("Settings for this run:\n")
            print("  T   Tasks to perform?\t: ")
            if len(self.task_list) < 2:
                print("      No parameters are selected for plotting.")
            else:
                for i in range(1,len(self.task_list),1):
                    print(f"      {i}: {self.task_list[i]['title']} - {self.task_list[i]['info']['description']} " + 
                        f"{self.task_list[i]['info']['k-pattern'] if self.task_list[i]['info']['k-pattern'] else ''};")
            print("  Q   to quit")
            print()
            print("Y to accept these settings, type the letter for one to change or Q to quit")

########################################################################
class PlotSettings:
    def __init__(self):
        pass
        
    def __call__(self, gbk_file="", default=False):
        oInterface = Interface(default=default)
        if gbk_file:
            oInterface.options['-g'] = gbk_file
        return oInterface.options

########################################################################       
# Validator
class InterfaceValidator:
    def __init__(self):
        self.IO = seq_io.IO()
        self.options = {"-c":"",# name of the scenario
                   "-g":"",     # GBK file
                   "-u":"",     # use BLASTN to filter rrn clusters
                   "-n":"",     # use BLASTN to search tRNA on the borders of GI
                   "-l":0,      # sliding window length
                   "-b":0,      # sliding window big step
                   "-m":0,      # sliding window medium step
                   "-s":0,      # sliding window small step
                   "-r":0,      # refeence window step
                   "-e":"",     # refeence window step
                   "-i":"",     # input folder
                   "-o":"",     # output folder
                   "-f":"",     # save GI sequemces: no/fasta/gbk/gbk+fasta
                   "-v":"",     # save SVG file
                   "-w":"",     # motif setting
                   "-p":75,     # promoter length
                   "-std":"off" # exclude reverse-complement strand
                }
        
    def validate(self, options, scenario):
        for key in options:
            if key not in self.options:
                tools.alert(f"Wrong option {key};")
                return 
            self.options[key] = options[key]
        if self.options['-g']:
            for fname in self.options['-g']:
                if not os.path.exists(fname):
                    tools.alert(f"File {self.options['-g']} does not exist!")
                    return
        if self.options["-c"] not in scenario:
            tools.alert(f"Wrong scenarium name {self.options['-c']}\nMust be in {'|'.join(scenario)};")
            return 
        if self.options["-u"] not in ("yes", "no"):
            tools.alert(f"Wrong option -u {self.options['-u']}\nMust be in yes|no;")
            return
        if self.options["-n"] not in ("yes", "no"):
            tools.alert(f"Wrong option -n {self.options['-n']}\nMust be in yes|no;")
            return
        try:
            self.options["-l"] = int(self.options["-l"])
        except:
            tools.alert("Sliding window length -l must be an integer;")
            return
        try:
            self.options["-b"] = int(self.options["-b"])
        except:
            tools.alert("Sliding window step -b must be an integer;")
            return
        try:
            self.options["-m"] = int(self.options["-m"])
        except:
            tools.alert("Sliding window step -m must be an integer;")
            return
        try:
            self.options["-s"] = int(self.options["-s"])
        except:
            tools.alert("Sliding window step -s must be an integer;")
            return
        try:
            self.options["-r"] = int(self.options["-r"])
        except:
            tools.alert("Reference window step -r must be an integer;")
            return
        if self.options["-s"] > self.options["-m"]:
            tools.alert("Small window step -s must be smaller or equal to the medium step -m;")
            return
        if self.options["-m"] > self.options["-b"]:
            tools.alert("Medium window step -m must be smaller or equal to the big step -b;")
            return
        if self.options["-b"] > self.options["-l"]:
            tools.alert("Big window step -b must be smaller or equal to the window length -l;")
            return
        if not os.path.exists(self.options["-i"]):
            tools.alert(f"Folder -i {self.options['-i']} does not exist;")
            return
        if not os.path.exists(self.options["-o"]):
            folder_name = self.IO.new_folder(self.options["-o"])
            if folder_name is None or not os.path.exists(folder_name):
                tools.alert(f"Folder -o {self.options['-o']} does not exist and cannot be created;")
                return
            self.options["-o"] = folder_name
        if self.options["-f"] not in ("no", "fasta", "gbk", "fasta+gbk"):
            tools.alert(f"Wrong option -f {self.options['-f']}\nMust be in {'|'.join(['no', 'fasta', 'gbk', 'gbk+fasta'])};")
            return
        if self.options["-v"].upper() not in ("YES", "NO"):
            tools.alert(f"Wrong option -v {self.options['-v']}\nMust be in yes|no;")
            return
        try:
            self.options["-p"] = int(self.options["-p"])
        except:
            tools.alert(f"Wrong option '-p' {self.options['-p']} for promoter length. Must be a positive integer!")
            return
        if self.options['-p'] < 0:
            tools.alert(f"Wrong option '-p' {self.options['-p']} for promoter length. Must be a positive integer!")
            return
        if self.options['-std'] not in ('off','leading','lagging'):
            self.options['-std'] = 'off'
        for path in self.options['-g']:
            if not os.path.exists(path):
                tools.alert(f"File {path} does not exists!")
                return
        return self.options
                
###############################################################################

if __name__ == "__main__":
    oInterface = Interface()
    
