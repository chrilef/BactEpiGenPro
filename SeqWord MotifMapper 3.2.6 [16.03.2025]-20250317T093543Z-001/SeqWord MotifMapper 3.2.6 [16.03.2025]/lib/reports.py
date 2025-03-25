import copy
import container, tools

########################################################################
class Table(container.Collection):
    def __init__(self, title=""):
        container.Collection.__init__(self,title)

########################################################################
class GeneIslandCollection(container.Collection):
    def __init__(self, title=""):
        container.Collection.__init__(self,title)

########################################################################
class GeneIsland:
    def __init__(self, start, end, title="", stat={}, seq=""):
        try:
            self.start = int(start)
            self.end = int(end)
        except:
            raise ValueError("Start and End values passed to GeneIsland object must be integers!")
        if not title:
            title = f"{start}..{end}"
        self.title = self.location = title
        self.genes = container.Collection()
        self.stat = {}
        if stat:
            self.stat = copy.deepcopy(stat)
        self.seq = ""
            
########################################################################
class Gene:
    def __init__(self, title):
        self.title = title
        self.locus_tag = ""
        self.name = ""
        self.start = 0
        self.end = 0
        self.product = ""

    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)  # Dynamically set the attribute value
            else:
                tools.alert(f"'Value {key}' is not a valid attribute of reports.Gene class!")

########################################################################
class StatPlotReport:
    def __init__(self, sample_size, title="", motif_setting="", verified=False):
        # Sample size
        self.sample_size = sample_size
        # Verification of locations of modified nucleotides
        self.verified = verified
        # Description of motifs used for modified site filtration
        self.motif_setting = motif_setting
        # Genome title
        self.title = title
        
        # Properties like 'GC', 'GCS', ...
        self.tasks = {}
        
        # Mobile genetic elements
        self.MGE = {'GI':GeneIslandCollection(self.title), 'statistics':StatUnit("MGE")}
        
        # coding/non-coding statistics
        self.genestat = {'coding':StatUnit("coding"), 'promoter':StatUnit("promoter"), 'noncoding':StatUnit("noncoding")}
        
        # contig statistics - list of StatUnit objects
        self.contigs = container.Collection("contigs")
        
    def get(self, motif_setting="<unidentified>"):
        return self.stat_report(motif_setting)
        
    def get_report(self, report="stat report", motif_setting="<unidentified>"):
        if report == "stat report":
            return self.stat_report(motifs)
        else:
            raise ValueError(f"Unknown report type {report}!")
            
    def stat_report(self, motif_setting="<unidentified>"):
        # Description of motifs used for modified site filtration
        if motif_setting != "<unidentified>" or not self.motif_setting:
            self.motif_setting = motif_setting
            
        report = [f"\n{'#' * 60}\nSTATISTICS: {self.title}, {self.sample_size} modified nucleotides\n\tfiltered by motifs {self.motif_setting}.\n",
            f"\t\tLocations of modified bases were {'' if self.verified else 'NOT '}VERIFIED"]
        
        # Add task report
        if self.tasks:
            report.append("\tGenome properties:")
            for task in self.tasks:
                for stat_method in self.tasks[task]:
                    for algorithm in self.tasks[task][stat_method]:
                        oStatUnit = self.tasks[task][stat_method][algorithm]
                        report.append(f"\t\t{task}\t{stat_method} {(algorithm)}\t{oStatUnit.correlation:.2f} " + 
                            f"[from {oStatUnit.confidence_interval[0]:.2f} to {oStatUnit.confidence_interval[1]:.2f}] with p-value = {oStatUnit.p_value:.3f}")
            report.append("")
            
        # Add MGE report
        if len(self.MGE['GI']):
            report.append("\tDistribution across MGE:")
            for oGI in self.MGE['GI']:
                report.append(f"\t\t[{oGI.start}..{oGI.end}];")
            report.append("")
            report.append(f"\tZ-score: {self.MGE['statistics'].Z_score:.2f} \u00B1 {self.MGE['statistics'].std_error:.2f} with p-value = {self.MGE['statistics'].p_value:.3f}\n" +
                f"\t\t\texpected / observed: {round(self.MGE['statistics'].expected_number)} / {round(self.MGE['statistics'].observed_number)}\n")
            
        # Add CDS/nonCDS distribution
        for key in ['coding', 'promoter', 'noncoding']:
            if self.genestat[key].title:
                report.append(f"\tDistribution across {key}:")
                report.append(f"\t\tZ-score: {self.genestat[key].Z_score:.2f} \u00B1 {self.genestat[key].std_error:.2f} with p-value = {self.genestat[key].p_value:.3f}\n" +
                    f"\t\t\texpected / observed: {round(self.genestat[key].expected_number)} / {self.genestat[key].observed_number}\n")
                
        # Add contig statistics
        if len(self.contigs):
            report.append("\tDistribution across contigs:")
            for contig in self.contigs:
                report.append(f"\t\t[{contig.start}..{contig.end}] " +
                    f"Z-score: {contig.Z_score:.2f} \u00B1 {contig.std_error:.2f} with p-value = {contig.p_value:.3f}\n" +
                    f"\t\t\texpected / observed: {round(contig.expected_number)} / {contig.observed_number}\n")
            if self.contigs.para['anova'] <= 0.05:
                report.append(f"\tThe hypothesis of biased distribution of modified sites across contigs is confimed with p-value {self.contigs.para['anova']};\n")
            else:
                report.append(f"\tThe hypothesis of biased distribution of modified sites across contigs is rejected with p-value {self.contigs.para['anova']};\n")
                
        return "\n".join(report)
        
    def set_MGE_stat(self, **kwargs):
        self.MGE['statistics'].set(**kwargs) 
        
    def set_MGE_loci(self, locations, genes, seq=""):
        if len(locations) != len(genes):
            raise ValueError("Length of MGE locations and genes must be the same!")
            
        for i in range(len(locations)):
            locus, k_mer_stat = locations[i]
            start, end = [int(v) for v in (locus.split("-") if locus.find("-") > -1 else locus.split(".."))]
            oGeneIsland = GeneIsland(start, end, locus, k_mer_stat, seq[start - 1:end] if end <= len(seq) else "")
            # Append gene objects to oGeneIsland object
            for cds in genes[i]:
                # sdc is a dict: {'locus_tag': 'SA150_1349', 'name': '.', 'product': 'Phage tail tape measure', 'start': 1484615, 'end': 1490815}
                # Must have either 'locus_tag' or 'start' and 'end'
                oGene = Gene(title = cds['locus_tag'] if 'locus_tag' in cds else f"{cds['start']}..{cds['end']}")
                oGene.set(**cds)
                oGeneIsland.genes.append(oGene)
            # Append oGeneIsland object to GeneIslandCollection
            self.MGE['GI'].append(oGeneIsland)
            
    def set_task(self, title, statistics):
        self.tasks[title] = {}
        for key in statistics:
            if key == 'correlation':
                correlation, conf_interval, p_value, algorithm = statistics[key]
                self.tasks[title][key] = {algorithm:StatUnit(f"Correlation {algorithm}")}
                self.tasks[title][key][algorithm].set(correlation = correlation, 
                    confidence_interval = conf_interval,
                    p_value = p_value,
                    algorithm = algorithm)
                    
    def set_CDS_nonCDS_stat(self, algorithm, statistics):   # algorithm = Z-score | LD
        for key in statistics:
            self.genestat[key].set(**dict(zip([algorithm, 'std_error', 'p_value', 'expected_number', 'observed_number'], statistics[key])))
        
    def set_contig_statistics(self, contigs, anova_p_value):
        self.contigs.para['anova'] = anova_p_value
        for contig in contigs:
            oContig = StatUnit(contig['title'])
            oContig.set(**dict(zip(['start', 'end', 'Z_score', 'std_error', 'p_value', 'expected_number', 'observed_number'],
                [contig['start'], contig['end']] + list(contig['Z-score']))))
            self.contigs.append(oContig)

########################################################################
class StatUnit:
    def __init__(self, title=""):
        self.title = title
        self.average = None
        self.standard_deviation = None
        self.variance = None
        self.std_error = None
        self.confidence_interval = []
        self.p_value = None
        self.adjasted_p_value = None
        self.anova_p_value = None
        self.Z_score = None
        self.LD = None
        self.expected_number = 0
        self.observed_number = 0
        self.lincage_disequilibrium = None
        self.correlation = {}
        self.algorithm = ""
        self.start = 0
        self.end = 0
        
    def set(self, **kwargs):
        for key, value in kwargs.items():
            key = key.replace("-","_")
            if hasattr(self, key):
                setattr(self, key, value)  # Dynamically set the attribute value
            else:
                raise ValueError(f"'{key}' is not a valid attribute of StatUnit")
