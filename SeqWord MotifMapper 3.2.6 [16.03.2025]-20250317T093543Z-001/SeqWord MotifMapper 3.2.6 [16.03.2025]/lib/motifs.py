import sys, os, string, re, copy
from functools import reduce
sys.path.append(os.path.join("..","lib"))
import blast, seq_io, tools, progressbar

##############################################################################################################
class Motif:
    def __init__(self,reference,binpath,tmppath,word="",motif="",motifs_or_sites='sites',modbase_location='0',reverse_modbase_location='0',
                score_cutoff=21,promoter_length=75,context_mismatches=0,motif_mismatch="No",filtered_loci=[],modified_or_unmodified_sites='M',
                strand=0, max_sites_for_verification=3000):
        self.oIO = seq_io.IO()
        self.flg_show_modified_sites = True if modified_or_unmodified_sites == "M" else False
        self.filtered_loci = filtered_loci
        self.oBlast = None
        self.gbk = None
        self.gff = []
        self.reference = reference
        self.binpath = binpath
        self.source_path = tmppath
        self.word = word
        self.motif=motif
        self.search_mode = motifs_or_sites
        self.modbase_location = str(modbase_location)
        self.reverse_modbase_location = str(reverse_modbase_location)
        self.score_cutoff = 21
        self.strand = strand
        self.max_sites_for_verification = int(max_sites_for_verification)
        try:
            self.score_cutoff = int(score_cutoff)
        except:
            path
        self.promoter_length = int(promoter_length)
        self.context_mismatches = context_mismatches
        self.motif_mismatch = motif_mismatch
        self.output_table = []              # Output table
        self.found_sites = []               # List of found modified entries
        self.unmodified_sites = []          # List of found unmodified entries
        self.overlapped_motifs = []         # List of overlapped motifs sharing the same modified site
        self.loci_in_refseq = []            # Canonical motifs found in the reference sequence
        self.loci_in_GFF = []               # Canonical motifs found in GFF and passed through filterring
        self.filtered_loci_in_refseq = []   # Canonical motifs of the reference sequence passed through filtering
        self.modified_loci_in_GFF = []      # Canonical motifs overlapping with not masked loci of the reference sequence and modified
        self.unmodified_loci_in_GFF = []    # Canonical motifs overlapping with not masked loci of the reference sequence and unmodified
        self.masked_loci_in_refseq = []     # Canonical motifs found in GFF overlapping masked loci of the reference sequence
        self.masked_loci_in_GFF = []        # Canonical motifs overlapping with not masked loci of the reference sequence
        self.novel_GFF_modpredictions = []  # Novel modified sites predicted in GFF, which were not present in refseq
        self.not_found_loci = []            # Modified motifs in GFF file, which were not found in the reference sequence
        self.num_instances = 0              # Number of canonical motifs found in GFF file, which match to the refseq and not masked
        self.num_found_sites = 0            # Number of found modified/unmodified canonical motifs
        self.num_expected_sites = 0         # Number of expected modified/unmodified sites
        self.num_found_motifs = 0           # Number of found motifs
        self.num_overlapped_motifs = 0      # Number of overlapped motifs sharing the same modified site
        self.statistics = []                # Programrun statistics
        self.table_titles = ['Motif','Strand','Location','Site','Match','Description','Site location relative to TSC','Annotation']
        self.success = False
        
    def __str__(self):
        return self.tostring()
    
    def execute(self, gff_data, gbk_data, echo=False):
        # Store the initial GFF entries
        self.success = False
        if isinstance(gbk_data, dict):
            gbk_data = gbk_data['object']
        # Check GFF data, which can be either a path, or the parced data
        if isinstance(gff_data, str):
            gff_path = gff_data
            if not os.path.exists(gff_path):
                return self.success
            gff = self.oIO.readGFF(path=gff_path, mode='dictionary')
            gff['input file'] = gff_path
        else:
            gff = gff_data
        # Record GFF data
        self.gff = gff if isinstance(gff, dict) else {'Heading':[],'Body':gff}
        # Check that the entry list is not empty
        if not self.gff['Body']:
            tools.alert("Filtered GFF records are empty!")
            return self.success
        
        # Check GBK data, which can be either a path, or the parces data and store locally at self.gbk
        if isinstance(gbk_data, str):
            self.gbk = self.oIO.read(gbk_infile,"genbank")
        else:
            self.gbk = gbk_data
        
        # Create output table
        self.output_table = self.search(generic_motif=self.word, 
            direct=self.modbase_location, reverse=self.reverse_modbase_location, 
            promoter=self.promoter_length)
        # Record result statistics
        self.output_table.insert(0,["\n" + f"Genome {self.reference}" "\n" + '\n'.join(self.statistics) + ";\n"])
        self.success = True
        return self.success
        '''
        try:
            # Create output table
            self.output_table = self.search(generic_motif=self.word, 
                direct=self.modbase_location, reverse=self.reverse_modbase_location, 
                promoter=self.promoter_length)
            # Record result statistics
            self.output_table.insert(0,[f"\nGenome {self.reference}\n{'\n'.join(self.statistics)};\n"])
            self.success = True
            return self.success
        except Exception as e:
            tools.alert(f"Error: {e}")
            return self.success
        '''

    def get_run_statistics(self, sites, flg_palindrom):
        # Find modified site location relatively to motif start location
        def _modified_site_location(motif_borders, site_location):
            start, end = [int(v) for v in (motif_borders.split("..") if motif_borders.find("..") > -1 else motif_borders.find("-"))]
            return int(site_location) - start + 1
            
        '''
        self.statistics = [
            "Number of motifs found in reference sequence: %d" % len(self.loci_in_refseq),
            "Number of filtered motifs in reference sequence after masking: %d" % len(self.filtered_loci_in_refseq),
            "Expected number of modified sites: %d" % self.num_expected_sites,
            "Number of sites found in GFF file: %d" % len(self.loci_in_GFF),
            "Number of masked sites taken from GFF file: %d" % len(self.masked_loci_in_GFF),
            "Number of sites prdecited in GFF file and not found: %d" % len(self.not_found_loci),
            "Number of sites found in GFF file, which were not present in the reference sequence: %d" % len(self.novel_GFF_modpredictions),
            "Total number of instances: %d" % self.num_instances
        ]
        '''
        self.statistics = []
        if self.flg_show_modified_sites:
            modified_point_statistics = tools.dereplicate_and_count([[int(site['Absolute point']) if
                site['Strand'] == site['Motif strand'] else -int(site['Absolute point']), 1] for site in sites])
            t = ""
        else:
            modified_point_statistics = tools.dereplicate_and_count([
                [_modified_site_location(motif_borders = site['Location'], site_location = site['Site']), 1] for 
                site in sites])
            t = "un"
        
        modified_locations = (
            "\t" + "; ".join(
                [f"{t}modified at location {item[0]} = {item[1]}" for item in sorted(modified_point_statistics, reverse=True)]
            )
        )

        self.statistics += [f"Motif {self.word} is modified at sites {', '.join([str(item[0]) for item in sorted(modified_point_statistics, reverse=True)])};",
            f"{self.num_found_motifs} fully {t}methylated motifs out of {len(self.filtered_loci_in_refseq)} found;" if self.search_mode == "motifs" 
                else f"{len(sites)} {t}modified sites were identified out of {self.num_expected_sites} found motifs;",
            modified_locations
        ]

    def modified_sites(self):
        if self.reverse_modbase_location != "0":
            return "%s and %s" % (self.modbase_location,self.reverse_modbase_location)
        return str(self.modbase_location)
        
    def tostring(self, text_to_insert=""):
        # Final modifications to the output
        def _modify_output(text_to_insert):
            # Create a copy of the table heading
            table_titles = copy.deepcopy(self.table_titles)
            # Add additional titles
            table_titles = table_titles[:6] + ["Motif strand", "Score", "Coverage", "Site location relative to TSC", "Annotation"]
            
            # Combine motif finding results with statplot data, if available
            output_table = [self.output_table[0], ["\n" + text_to_insert if isinstance(text_to_insert,str) else "\n".join(text_to_insert)], table_titles]
            output_table += [[self.output_table[i][key] for key in self.table_titles] for i in range(1,len(self.output_table),1)]

            # Remove motif strand column
            output_table[2] = output_table[2][:6] + output_table[2][7:]
            for i in range(3,len(output_table)):
                values = output_table[i][5].split("\t")
                output_table[i][5] = "\t".join([values[0]] + values[2:])
            
            return output_table            
            
        # If the process of finding motifs was not successful, return an empty string
        if not self.success:
            return ""
        # Final modification to the output table
        output_table = _modify_output(text_to_insert)
        # Return output as a tab separated table
        return "\n".join(map(lambda item: "\t".join(item), output_table))
    
    # Check the results of site verification
    def get_entries(self):
        def _select_location(location):
            '''
            If location is within coding sequence and in the range of other gene promoter, return within CDS location.
            '''
            if location == 'n/a':
                return location
            location = [int(v) for v in str(location).split(",")]
            return max(location)
            
        entries = []
        if self.flg_show_modified_sites:
            sites = self.found_sites
        else:
            sites = self.unmodified_sites
        for found_entry in sites:
            entry = copy.deepcopy(found_entry)
            if 'verified site' in entry:
                if isinstance(entry['verified site'],dict):
                    entry['site'] = int(entry['verified site']['Site'])
                    locations = [int(v) for v in entry['verified site']['Location'].split("..")]
                    entry['start'] = min(locations)
                    entry['end'] = max(locations)
                    entry['nucleotide'] = entry['verified site']['nucleotide']
                    entry['data']['context'] = entry['verified site']['context']
                    entry['data']['Annotation'] = entry['verified site']['Annotation']
                    entry['data']['Site location relative to TSC'] = _select_location(entry['verified site']['Site location relative to TSC'])
                    entry['verified site'] = True
            else:
                entry['verified site'] = False
            entries.append(entry)
        return entries
   
    def get_modbase_list(self):
        entries = self.get_entries()
        if not all([item['verified site'] for item in entries]) and self.max_sites_for_verification > 0 and len(entries) <= self.max_sites_for_verification:
            # Verify unverified sites
            entries = self.verify_sites(sites=entries, gff_data=self.gff)
        # in the format [[start,end,strand,{data}],...]
        entries = [[int(entry['start']), int(entry['end']), 
            1 if str(entry['strand']) in ('DIR','+','1') else -1, 
            {'Absolute location':entry['site'], 'Site':entry['site'], 'Annotation':entry['data']['Annotation'], 
                'Site location relative to TSC':entry['data']['Site location relative to TSC']}
                if 'data' in entry else {}] 
            for entry in entries]
        return entries
        
    def get_site_entries(self):
        return [site for site in self.output_table[1:]]

    def get_word(self):
        return self.word
        
    def get_motif(self):
        return self.motif
    
    def get_modbase_location(self):
        return self.modbase_location
    
    def get_rev_modbase_location(self):
        return self.reverse_modbase_location

    def get_modbase_dict(self):
        def get_color(is_matched):
            if is_matched == "Yes":
                return "blue"
            elif is_matched == "None":
                return "green"
            return "red"
        
        modbases = {}
        for i in range(1,len(self.output_table)):
            item = self.output_table[i]
            #moltype,direction,score,coverage = item['Description'].split("|")[0].split("\t")
            moltype,direction,score,coverage = self.get_entry_field(item, 'Description').split("|")[0].split("\t")
            modbases[i-1] = [self.get_entry_field(item,'Location'),
                            get_color(self.get_entry_field(item, 'Match')),"1" if self.get_entry_field(item,'Strand').upper() in ('DIR',"1","+") else "-1",
                            ["%s\t%s\t%s\t%s\t%s" % (self.get_entry_field(item,'Location'),moltype,score,coverage,self.get_entry_field(item,'Motif'))]
                ]
        return modbases

    def get_annotation(self,gene,strand=0):
        if gene == "non-coding":
            return gene
        description = "[%d..%d]" % (int(gene.location.start),int(gene.location.end))
        if 'locus_tag' in gene.qualifiers:
            description = "%s %s" % (gene.qualifiers['locus_tag'][0],description)
        if 'gene' in gene.qualifiers:
            description += "; %s" % gene.qualifiers['gene'][0]
        if 'product' in gene.qualifiers:
            description += "; %s" % gene.qualifiers['product'][0]
        return description
    
    def overlap_with_filtered_regions(self,location):
        for locus in self.filtered_loci:
            if location >= locus[0] and location < locus[1]:
                return True
        return False

    # Find genes overapping modified site location on the same strand
    def gene_overlap(self, cds,location,promoter_length,strand):
        # Genes located on the DNA strand opposit to the site strand are not considered
        if cds.strand != strand:
            return False
        # Genes are selected if the site is located in its body or in the specified promoter region.  
        if ((cds.strand==1 and int(cds.location.start)-promoter_length <= location <= int(cds.location.end)) or
            (cds.strand==-1 and int(cds.location.start) <= location <= int(cds.location.end)+promoter_length)):
            if cds.strand==1:
                cds.qualifiers['TSS'] = [location - int(cds.location.start)]
            else:
                cds.qualifiers['TSS'] = [int(cds.location.end) - location]
            if cds.qualifiers['TSS'][0] >= 0:
                return cds
            cds.qualifiers['promoter'] = ["Yes"]
            return cds
        return False    

    # Create a dictionary object containing all information about the modified site
    def get_site(self, entry, length, k, promoter, genes, complement_site, motif="", 
        strand_adjustment=True, flg_palindrom = False, relative_positions=True, echo=False):
                    
        cds = []

        # Read parameters
        strand = entry['strand'] if 'strand' in entry else entry['Strand']
        if 'Location' in entry:
            start,end = [int(v) for v in (entry['Location'].split("..") if entry['Location'].find("..") > -1 else entry['Location'].split("-"))]
        else:
            start = entry['start'] if 'start' in entry else entry['Start']
            end = entry['end'] if 'end' in entry else entry['End']
        site = entry['site'] if 'site' in entry else (entry['Site'] if 'Site' in entry else start)
        context = entry["data"]["context"] if "data" in entry else tools.get_locus_by_coordinates(seq=seq, start=start, end=end,
            flg_reverse = True if strand=="-" else False,
            flg_complement = True if strand=="-" else False).upper()
        coverage = entry["data"]["coverage"] if "data" in entry else "na"
        if 'Motif' in entry and not motif:
            motif = entry['Motif']
        match = entry['match'] if 'match' in entry else (entry['Match'] if 'Match' in entry else 'None')
        score = entry['score'] if 'score' in entry else "0"
        modtype = entry['modtype'] if 'modtype' in entry else "unmodified"
        
        motif_strand_location = -1 if complement_site else 1
        if str(strand) in ("-", "REV", "-1"):
            strand = -1
            site_strand_location = -1
        else:
            strand = 1
            site_strand_location = 1
        cds =["nd"]
        description = annotation = ""
        TSS = 'n/a'
        location = f"{start}..{end}"
        # Set motif coordinates containing modified nucleotide
        if motif:
            if relative_positions:
                site_location = int(start) + k - 1 if strand == 1 else int(end) - k + 1
            else:
                site_location = int(start) + k - 1
            # Strand adjusments
            if strand_adjustment:
                # Correct site location
                if site_strand_location == -1:
                    site_location -= 1
                # Correct site strand location
                if complement_site and site_strand_location == 1:
                    site_strand_location = -1
                elif not complement_site and site_strand_location == -1:
                    site_strand_location = 1
            # Modified base
            #base = motif[(site_location - start) if motif_strand_location == 1 else (end - site_location)].upper()
            #base = motif[site_location - start].upper()
            base = context[20].upper()
                
        elif context:
            # Bulk sites given without motifs for verification
            site_location = site
            motif = context[21-k:21-k+length]
            # Modified base
            base = context[20].upper()
            
        else:
            site_location = site
            # Modified base
            base = ""
        
        # Add information about gene, promoter or non-coding region
        if genes:
            cds = [gene for gene in genes if self.gene_overlap(cds=gene, location=int(site), promoter_length=promoter, strand=site_strand_location)]
            if not cds:
                cds =["non-coding"]
            else:
                TSS = ",".join([str(gene.qualifiers["TSS"][0]) for gene in cds])
            annotation = "|".join([self.get_annotation(gene,strand) for gene in cds])
        # Get complement sequence
        if complement_site and strand_adjustment:
            strand = int(-1 * strand)
        
        # Strand adjustment
        if flg_palindrom:
            motif_strand_location = int(-1 * motif_strand_location)
            
        if strand == -1 and (strand_adjustment or flg_palindrom):
            motif = tools.reverse_complement(motif, reverse=False)

        # Return dictionary
        oSite = dict(zip(self.table_titles,
                [motif,'DIR' if site_strand_location == 1 else 'REV',location,str(site_location),match,"%s\t%s\t%s\t%s" % 
                (modtype,'REV' if motif_strand_location == -1 else 'DIR',str(score),coverage),TSS,annotation]))
        # Add additional parameters
        oSite['Motif strand'] = 'REV' if motif_strand_location == -1 else 'DIR'
        oSite['Absolute point'] = k
        oSite['Base'] = base
        return oSite
                
    # Create a dictionary object containing all information about the unmodified site
    def get_unmodified_site(self, entry, promoter_length, genes, flg_palindrom = False, echo=False):
        motif = entry['word']
        point = int(entry['point'])
        start = int(entry['start'])
        end = int(entry['end'])
        site_location = int(entry['site'])
        base = entry['data']['context'][20].upper()
        score = 0
        location = f"{start}..{end}"
        site_strand_location = 1 if entry['strand'] == "+" else -1
        motif_strand_location = 1 if entry['motif_strand'] == "+" else -1
        cds =["nd"]
        description = annotation = ""
        TSS  = coverage = 'n/a'
        match = "None"
        modtype = "unmodified"
        if genes:
            cds = [gene for gene in genes if self.gene_overlap(cds=gene, location=int(site_location), promoter_length=promoter_length, strand=site_strand_location)]
            if not cds:
                cds =["non-coding"]
            else:
                TSS = ",".join([str(gene.qualifiers["TSS"][0]) for gene in cds])
            annotation = "|".join([self.get_annotation(gene,site_strand_location) for gene in cds])

        # Return dictionary
        oSite = dict(zip(self.table_titles,
                [motif,'DIR' if site_strand_location == 1 else 'REV', location, str(site_location), match, "%s\t%s\t%s\t%s" % 
                (modtype, 'REV' if motif_strand_location == -1 else 'DIR', str(score), coverage), TSS, annotation]))
        # Add additional parameters
        oSite['Motif strand'] = 'REV' if motif_strand_location == -1 else 'DIR'
        oSite['Absolute point'] = abs(point)
        oSite['Point'] = point
        oSite['Base'] = base
        return oSite
                
    def get_entry_field(self, item, field):
        return item[field] if field in item else item['verified site'][field]

    # Block of dereplicators
    '''
    Filter by exect nucleotide positions
    {'Motif': 'GTGCAC', 'Strand': 'REV', 'Location': '576132..576137', 'Site': '576134', 'Match': 'None', 
    'Description': 'unmodified\tDIR\t0\tna', 'Site location relative to TSC': 'n/a', 'Annotation': 'non-coding'}
    '''
    def dereplicate_sites(self, ls):
        def _get_location(d):
            return int(d['Site'] if 'Site' in d else (d['site'] if 'site' in d else d['start']))
        if len(ls) < 2:
            return ls
        ls.sort(key=lambda d: _get_location(d))
        for i in range(len(ls)-1,0,-1):
            if _get_location(ls[i-1]) == _get_location(ls[i]):
                del ls[i]
        return ls

    '''
    Dereplicate entries by locations at direct and reverse-complement strands
    {'genome': 'N.', 'method': 'kinModCall', 'modtype': 'm4C', 'start': '1633', 'end': '1633', 'score': '398', 'strand': '+', 'para': '.', 
    'data': {'coverage': '379', 'context': 'GACCGAAAAAGACCAGGGAGCACCTGAAAGCTTGGAGGACG', 'IPDRatio': '3.73', 'identificationQv': '446'}, 
    'nucleotide': 'C'}
    '''
    def dereplicate_entries(self, ls): 
        if len(ls) < 2:
            return ls
        ls.sort(key=lambda d: [int(d['start']),d['strand']])
        for i in range(len(ls)-1,0,-1):
            if (ls[i]["strand"]+str(ls[i]["start"]) == ls[i-1]["strand"]+str(ls[i-1]["start"]) or
                (i > 1 and ls[i]["strand"]+str(ls[i]["start"]) == ls[i-2]["strand"]+str(ls[i-2]["start"]))):
                del ls[i]
        return ls

    '''
    Dereplicate list of motifs by coordinates and strands
    '''
    def dereplicate_motif_list(self, ls):     # site = [start, end, strand]
        if len(ls) < 2:
            return ls
        ls.sort()
        for i in range(len(ls)-1, 0, -1):
            if ls[i] == ls[i-1]:
                del ls[i]
        return ls
    
    '''
    Dereplicate list of motifs by coordinates and strands
    '''
    def dereplicate_motifs(self, ls):     # site = [start, end, strand]
        if len(ls) < 2:
            return ls
        ls.sort()
        for i in range(len(ls)-1, 0, -1):
            if ls[i] == ls[i-1]:
                del ls[i]
        return ls
    
    def find_context(self,entry,seq,occupied_sites):
        context = entry['verified site']['context'] if 'verified site' in entry else entry["data"]["context"]
        word = entry['verified site']['Motif'] if 'verified site' in entry else entry['word']
        start = entry['verified site']['Site'] if 'verified site' in entry else int(entry['start'])
        strand = entry['strand']
        locations = []
        p = seq.find(context)
        if p > -1:
            locations.append(p)
        while p != -1:
            p = seq.find(context,p+1)
            if p > -1:
                locations.append(p)
        for i in range(len(locations)-1,-1,-1):
            p = locations[i]
            if (p >= len(context)+len(seq)//2 and ("%s%d..%d" % (strand,len(seq)-p-len(context)+1,len(seq)-p)) not in occupied_sites):
                locations[i] = [len(seq)-p-len(context)//2,"-",len(seq)-p-len(context)+1,len(seq)-p]
            elif ("%s%d..%d" % (strand,p+1,p+len(context))) not in occupied_sites:
                locations[i] = [p+1+len(context)//2,"+",p+1,p+len(context)]
            else:
                del locations[i]
                
        # Location was not found
        if not len(locations):
            return 0,strand
            
        # Several alternative locations were found
        if len(locations) > 1:
            # Sort alternative locations by closeness to the expected location
            locations.sort(key = lambda ls: abs(int(entry['start'])-ls[0]))
            
        start,strand,site_location,end = locations[0]
        occupied_sites.append(f"{strand}{start}..{end}")
        if int(start) > len(seq) // 2:
            start = len(seq) - int(end) + (1 if 'verified site' in entry else 21)
            #print("motifs:535", start, strand, entry)
        return start,strand

    # Verify base location using BLASTN
    def match(self, seq, entry, n, occupied_sites, dbname):        
        start,strand = self.find_context(entry,seq,occupied_sites)
        if start:
            return start,strand,True
        seq = seq[:len(seq)//2]
        context = entry["data"]["context"]
        context_file = "%s_record.fa" % dbname
        self.oIO.save(">record %d\n%s" % (n,context),os.path.join(self.source_path,context_file),"text")
        fasta_file = dbname+".fa"
        if not os.path.exists(os.path.join(self.source_path,fasta_file)):
            self.oIO.save(">genome\n%s" % seq,os.path.join(self.source_path,fasta_file),"text")
            self.oBlast.create_db(fasta_file,dbname)
        self.oBlast.execute(context_file,dbname)
        hsps = [hsp for hsp in self.oBlast.get_matches(len(context), self.context_mismatches) if f"{hsp.get_strand()}{hsp.sbjct_start}..{hsp.sbjct_end}" not in occupied_sites]

        if not hsps:
            return 0,0,False
        matches = [hsp for hsp in hsps if entry['word'].upper() in hsp.sbjct.upper()]

        if not matches:
            matches = hsps
        if len(matches) > 1:
            matches = [[len(context) - hsp.identities, hsp.measure_distance(entry['start'], len(context) // 2), hsp] for hsp in matches]
            matches.sort()
            matches = [ls[-1] for ls in matches]
        top_match = matches[0]
        strand = top_match.get_strand()
        start = int(top_match.sbjct_start-top_match.query_start+1+len(context)//2)
        occupied_sites.append("%s%d..%d" % (strand,top_match.sbjct_start,top_match.sbjct_end))
        return start,strand,top_match.sbjct.upper().find(entry['word'].upper()) > -1
        
    def search(self, generic_motif, direct="", reverse="", promoter=0, sequence="", context_length=41, echo=False):
        # Identify sought motifs in reference sequence
        
        def _set_loci(motif, wlength, seq):
            # Check motif information
            if not motif:
                motif = self.word       # motif like G[A,G,T]GC[A,C,T]C
            if not motif:
                return []

            tools.msg("\tSearching for instances of motifs...")
            # Identify locations of the given motif in the reference sequence
            p = 0
            while p < len(seq)-1:
                locus = re.search(motif,seq[p:])
                if not locus:
                    break
                p += locus.span()[0]
                if p < len(seq)/2:
                    self.loci_in_refseq.append([p+1,p+wlength,"+"])
                else:
                    self.loci_in_refseq.append([len(seq)-p-wlength+1,len(seq)-p,"-"])
                p += 1
                
            # Remove motifs found within filtered regions
            self.filtered_loci_in_refseq = self.get_available_loci()
            self.masked_loci_in_refseq = [locus for locus in self.loci_in_refseq if locus not in self.filtered_loci_in_refseq]
            
        def _is_palindrom(motif, points, wlength):
            # Check identity of the motif with its reverse complemet
            if motif.upper() != tools.reverse_complement(motif).upper():
                return False
            # Check symmetry of modified sites, but allow no sites or 1 site
            if len(points) > 1 and len(points) % 2 == 1:
                return False
            points = sorted([int(v) for v in points])
            # Check if site symmetry
            for i in range(len(points)//2):
                if points[i] + wlength + 1 != points[-i - 1]:
                    return False
            # The motif is palindromic
            return True

        def _fitting(entry):
            score = int(entry['score'])
            if self.score_cutoff <= score:
                return True
            return False
    
        def _masking(entry):
            if self.filtered_loci != [] and self.overlap_with_filtered_regions(entry["start"]):
                return True
                
        def _get_base(motif, point):
            if not point:
                return ""
            base = motif[abs(point) - 1].upper()
            return base if point >= 0 else tools.reverse_complement(base)
        
        # Reset all counted parameters
        self.reset()
        
        # Counter of rounds
        round_counter = 0
        
        # Take GFF records
        ls_sites = []
        if not self.gff or "Body" not in self.gff or not self.gff["Body"]:
            return ls_sites
        '''
        Processing locations of direct and reverse modified sites within the given motif
        For example: generic_motif,direct,reverse = GDGCHC, 4, 3
        They are converted to points [4, -3]
        '''
        points = []
        if direct:
            points += [int(v) for v in direct.split(",") if v.strip()]
        if reverse:
            points += [-int(v) for v in reverse.split(",") if v.strip()]
        
        # Dereplication is needed when both points are 0
        points = tools.dereplicate(points)
        
        # Bases expected at modified points
        expected_bases = dict(zip(points,[_get_base(motif=generic_motif, point=point) for point in points]))
            
        # Word length
        wlength = len(generic_motif)
        
        # Check if the motif is a palindrom, like GDGCHC
        flg_palindrom = _is_palindrom(motif=generic_motif, points=points, wlength=wlength)
        
        if flg_palindrom and len(points) == 1 and points[0]:
            # Add a symetric point
            v = int(points[0])
            sign = -(v // abs(v))
            additional_point =  sign * (wlength - abs(v) + 1)
            points.append(additional_point)
            points.sort()

        # Convert generic motif like GDGCHC to motif tamplete G[A,G,T]GC[A,C,T]C
        motif = tools.compile_motif(generic_motif)
        
        # Convert motif template into lists of all posible direct and reverse complement words
        dir_motifs = []
        if direct and direct != "0":
            dir_motifs = tools.inlist_motifs(motif)
        # Reverse motifs
        rev_motifs = []
        if reverse and reverse != "0":
            rev_motifs = tools.inlist_motifs(motif,True)
        
        '''
        Set gene and sequence information
        Sequence can be passed with 'sequence argument or taken from GBK file if provided
        '''
        genes = []
        dirseq = sequence
        revseq = ""
        if self.gbk:
            genes = [ft for ft in self.gbk.features if ft.type=="CDS"]
            if not dirseq:
                dirseq = str(self.gbk.seq)
                revseq = str(self.gbk.reverse_complement().seq)
        if not revseq:
            revseq = tools.reverse_complement(dirseq)
        
        if not dirseq:
            raise ValueError("Sequence information has not been passed!")
        # Referense sequence
        seq = dirseq + revseq
            
        # Identify all sought motif sequences self.loci_in_refseq:
        _set_loci(motif=motif, wlength=wlength, seq=seq)
        
        # Set BLASTN object
        self.oBlast = blast.BLAST("dna",self.binpath,self.source_path)
        # Generate temporary file name for BLASTN database
        random_file_name = tools.get_random_filename(self.source_path)
        
        # Identified sites to exclude duplications
        identified_sites = []   # Exact base positions
            
        # Loop across modified locations
        for k in points:
            point = k
            # For palindroms, reverse complement sites will be dublicated
            if flg_palindrom and k != 0 and round_counter >= len(points) // 2:
                continue

            occupied_sites = []
            #self.unmodified_loci_in_GFF.append([])
            ls_sites.append([])
            
            #### FILTERING
            searched_motifs = dir_motifs
            # rem - text heading for progress bar
            rem = f"{k}'s {motif[abs(k)-1]}" if k and motif else f"Verify {self.gff['Body']} sites "            # complement_site = False
            complement_site = False
            if k < 0:
                k += (wlength + 1)
                searched_motifs = rev_motifs
                complement_site = True
                rem = f"-{k}'s {tools.reverse_complement(motif)[abs(k)-1]}"

            # Identify modified site position assuming that the context sequence is 41 bp
            p = context_length // 2
            if k:
                p = context_length // 2 - k + 1
            
            # Progress bar heading
            message = f"{motif} at {rem}: " if motif else rem
            if len(message) > 50:
                message = "Run: "
            
            # Filtering GFF records that can be removed
            canonical_motifs = copy.deepcopy(self.gff["Body"])
            if searched_motifs:
                canonical_motifs = [entry for entry in canonical_motifs if entry["data"]["context"][p:p+wlength] in searched_motifs]
                
            # Dereplicate palindroms
            if flg_palindrom:
                canonical_motifs = self.dereplicate_entries(canonical_motifs)
            
            if not len(canonical_motifs):
                tools.alert(f"Motif {motif} was not found among modified sites!")
                break
            

            # Progress bar
            marker = ""
            bar = progressbar.indicator(len(canonical_motifs),message)
            for i in range(len(canonical_motifs)):
                '''
                entry = {'genome': '598_chr', 'method': 'kinModCall', 'modtype': 'm6A', 'start': '2399', 'end': '2399', 
                'score': '499', 'strand': '+', 'para': '.', 'data': {'coverage': '394', 'context': 'ATTATGTCTGACAATGAAGAAGACATTGATATCTTCTTTGC', 
                'IPDRatio': '4.94', 'identificationQv': '477'}, 'nucleotide': 'A'}
                '''               
                entry = canonical_motifs[i]
                # Set searched word
                entry['word'] = entry["data"]["context"][p:p+wlength]
                # Counting canonical motifs found in GFF
                self.loci_in_GFF.append(entry)
                if self.gbk:
                    # Verify location of modified site predicted in GFF by seqrching the context sequence in the provided reference sequence
                    location, strand, match = self.match(seq=seq, entry=entry, n=i, occupied_sites=occupied_sites, dbname=random_file_name)
                    # Processing motifs which were not found in the reference sequence
                    if not location or (self.motif_mismatch=="No" and not match):
                        entry['match'] = "No"
                        # # Generate site dictionary to save all relevant information and store it to the list of not found sites
                        self.not_found_loci.append(self.get_site(entry=entry, length=wlength, k=k, promoter=promoter, genes=genes, 
                            motif=entry['word'], complement_site=complement_site, flg_palindrom = flg_palindrom))
                        # Counting not found motifs
                        bar(i)
                        continue
                    else:
                        entry['match'] = "Yes"
                        
                    # Set exact site location and coordinates of the respective motif 
                    entry['strand'] = strand
                    entry['site'] = location
                    if strand == "+":
                        entry['start'] = location - k + 1
                        entry['end'] = location + wlength - k
                    else:
                        entry['start'] = location - wlength + k
                        entry['end'] = location + k - 1
                        
                        if not complement_site:
                            entry['word'] = "".join(reversed(list(entry['word'])))
                        
                    # Identify and set modified nucleotide
                    entry['nucleotide'] = seq[entry['site']-1].upper() if strand == "+" else tools.reverse_complement(seq[entry['site']-1]).upper()
                
                # Generate site dictionary to save all relevant information
                site = self.get_site(entry=entry, length=wlength, k=k, promoter=promoter, genes=genes, motif=entry['word'],
                        complement_site=complement_site, strand_adjustment=False, flg_palindrom = flg_palindrom)
                        
                # Add identified site if it was not added before; however, it allows dublication when whole motifs are searched
                if ((self.search_mode == "sites" and int(site['Site']) not in identified_sites) or self.search_mode == "motifs"):
                    # list of identified sites
                    ls_sites[-1].append(site)
                    # Store verified site location data to the current GFF entry
                    verified_site = copy.deepcopy(site)
                    verified_site['context'] = tools.get_locus_by_coordinates(seq=seq, start=int(site['Site'])-20, end=int(site['Site'])+20,
                        flg_reverse = True if verified_site['Strand']=="REV" else False,
                        flg_complement = True if verified_site['Strand']=="REV" else False).upper()
                    verified_site['nucleotide'] = verified_site['context'][20]

                    entry['verified site'] = verified_site
                    # Record current GFF entry
                    self.found_sites.append(entry)
                    # Additionally add the exact location of the verified site to avoid record dublication
                    if self.search_mode == "sites":
                        identified_sites.append(int(site['Site']))
            
                # Update the progress bar
                bar(i)
                
            # Stop progress bar
            bar.stop()
            round_counter += 1

            # Check correspondance of the found and the expected bases
            if point:
                for site in ls_sites[-1]:
                    if expected_bases[point] != site['Base']:
                        print(f"motifs:867. Found base {site['Base']} does not correspond to the expected base {expected_bases[point]} in motif {site['Motif']}.")
            
            # Record found modified sites 
            self.modified_loci_in_GFF.append([site for site in ls_sites[-1]])
            
            # Identify sites, which were predicted in GFF, but were abcent in the reference sequence
            all_sites = [f"{ls[2]}{ls[0]}..{ls[1]}" for ls in self.filtered_loci_in_refseq]
            modified_site_locations = [item["Strand"].replace("DIR", "+").replace("REV", "-") + item["Location"] for item in ls_sites[-1]]
            self.novel_GFF_modpredictions = [s for s in modified_site_locations if s not in all_sites]

            # Remove temporary files
            self.clean(self.source_path,random_file_name)
            
        
        # Calculate total number of identified modified sites including those, which were not found in the reference sequence
        self.num_instances = len(self.filtered_loci_in_refseq)+len(self.novel_GFF_modpredictions) # Filtered sites in refseq plus modifications newly predicted in GFF
        
        # Search for unmodified sites
        if not self.flg_show_modified_sites and points:
            tools.msg("Searching for unmodified sites...")
            
            if not flg_palindrom:
                # Count overlapped motifs like in sequence CRGKGA|TCMCYG for motif CRGKGATC for unmodified sites
                expected_sites = []
                for point in points:
                    expected_sites += [(int(start) + abs(point) - 1 if strand == "+" else int(end) - abs(point) + 1) for 
                        start, end, strand in self.filtered_loci_in_refseq]
                self.num_overlapped_motifs = len(expected_sites) - len(tools.dereplicate(copy.deepcopy(expected_sites)))

            # Collect all modified sites
            if len(ls_sites) > 1:
                found_modified_sites = list(reduce(lambda ls1,ls2: ls1 + ls2, ls_sites))
            else:
                found_modified_sites = ls_sites[0]
                
            # Identify unmodified sites
            if self.search_mode == "motifs" and len(points) > 1:
                # Search for motifs where all sites are unmodified
                unmodified_sites = self.get_unmodified_motifs(modified_sites=found_modified_sites, 
                    seq=seq, length=wlength, points=points, promoter=promoter, genes=genes, 
                    expected_bases=expected_bases, flg_palindrom = flg_palindrom)
            else:
                # Search for individual unmodified sites
                unmodified_sites = self.get_unmodified_sites(modified_sites=found_modified_sites, 
                    seq=seq, length=wlength, points=points, promoter=promoter, genes=genes, 
                    expected_bases=expected_bases, flg_palindrom = flg_palindrom)
            ls_sites = [unmodified_sites]
        
        # Verification of a bulk list of sites
        if not points:
            ls_sites.append(self.get_unmodified_sites([], seq, wlength, [0], promoter,genes))
        
        # Dereplicate list of motifs by coordinates and strand
        self.filtered_loci_in_refseq = self.dereplicate_motif_list(self.filtered_loci_in_refseq)
        
        if self.search_mode == "motifs" and len(points) > 1:
            # Searching for methylated motifs
            if flg_palindrom:
                self.num_expected_sites = len(self.filtered_loci_in_refseq) / len(points)
            else:
                self.num_expected_sites = len(self.filtered_loci_in_refseq)
            # Select fully modified motifs
            ls_sites = self.select_overlaps(ls_sites,points)
            # Dereplicate motifs
            ls_sites = self.dereplicate_sites(ls_sites)
            self.num_found_motifs = len(ls_sites)//len(points)
        else:
            # Set of methylated or unmethylated sites
            if len(ls_sites) > 1:
                ls_sites = reduce(lambda ls1,ls2: ls1+ls2, ls_sites)
            else:
                ls_sites = ls_sites[0]
                                
            if flg_palindrom:
                self.num_expected_sites = len(self.filtered_loci_in_refseq)
            else:
                # Count overlapped motifs like in sequence CRGKGA|TCMCYG for motif CRGKGATC for modified sites
                if not self.flg_show_modified_sites:
                    expected_sites = []
                    for point in points:
                        expected_sites += [(int(start) + abs(point) - 1 if strand == "+" else int(end) - abs(point) + 1) for 
                            start, end, strand in self.filtered_loci_in_refseq]
                    self.num_overlapped_motifs = len(expected_sites) - len(tools.dereplicate(copy.deepcopy(expected_sites)))
                
                self.num_expected_sites = len(self.filtered_loci_in_refseq) * len(points) - self.num_overlapped_motifs
            
        if self.flg_show_modified_sites and points:
            self.num_found_sites = (len(list(reduce(lambda ls1,ls2: ls1 + ls2, self.modified_loci_in_GFF))) if len(self.modified_loci_in_GFF) > 1 else
                (len(self.modified_loci_in_GFF[0]) if len(self.modified_loci_in_GFF) else 0))
            # Run statistics
            run_statistics = self.get_run_statistics(ls_sites, flg_palindrom)
        else:
            unmodified_sites = self.unmodified_loci_in_GFF
            unmodified_sites = self.dereplicate_sites(unmodified_sites)
            num_unmodified = len(unmodified_sites)
            self.num_found_sites = num_unmodified
            # Run statistics
            run_statistics = self.get_run_statistics(self.unmodified_loci_in_GFF, flg_palindrom)
        
        # Return outputs
        return ls_sites
    
    def get_unmodified_motifs(self, modified_sites ,seq, length, points, promoter, genes, expected_bases, flg_palindrom=False):
        modified_site_locations = [f"{item['strand']}:{item['verified site']['Location']}" for item in  self.found_sites]
        # Filter by strand
        if self.strand and not flg_palindrom:
            filtered_loci_in_refseq = [[start,end,strand] for start,end,strand in self.filtered_loci_in_refseq if
                (strand == "+" and self.strand == 1) or (strand == "-" and self.strand == -1)]
        else:
            filtered_loci_in_refseq = self.filtered_loci_in_refseq
        # Get unmodified site locations
        unmodified_sites = [ls for ls in filtered_loci_in_refseq if f"{ls[2]}:{ls[0]}..{ls[1]}" not in modified_site_locations]
        '''
        As an example, points in GDgcHC are [4, -3]
        '''
        sites = []
        for k in points: 
            for item in unmodified_sites:
                start, end, strand = item
                if flg_palindrom and strand == "-":
                    continue
                    
                motif_strand = strand
                site_strand = "+" if ((strand=="+" and k >= 0) or (strand=="-" and k < 0)) else "-"
                site_location = int(start) + abs(k) - 1
                entry = dict(zip(
                ['word', 'start', 'end', 'site', 'strand', 'score', 'modtype', 'match', 'data', 'site_strand', 'point', 'motif_strand'],
                ["", int(start), int(end), site_location, site_strand, 0, "unmodified", "None", {'context': 
                        tools.get_locus_by_coordinates(seq=seq, start = site_location - 20, end = site_location + 20, 
                        flg_complement = True if (site_strand == "-") else False, 
                        flg_reverse = True if (site_strand == "-") else False).upper(), 'coverage': 'na'}, 
                        site_strand, k, motif_strand]))
                              
                # Set word
                entry['word'] = tools.get_locus_by_coordinates(seq=seq, start=int(start), end=int(end), 
                    flg_complement = True if (site_strand == "-") else False, 
                    flg_reverse = True if (site_strand == "-") else False) 
                
                # Convert entry to site object
                site = self.get_unmodified_site(entry=entry, genes=genes, promoter_length=promoter, flg_palindrom=flg_palindrom)
                # Check whether unmodified base correspond to the expected base 
                if expected_bases[k] != site['Base']:
                    if not flg_palindrom:
                        tools.alert(f"Base identity control failed! In motif {site['Motif']} at location {entry['point']} " + 
                            f"base {expected_bases[int(entry['point'])]} is expected, but {site['Base']} is found.")
                    continue
                                            
                # Add site to the list
                sites.append(site)
        
        # Set run statistics
        self.get_run_statistics(sites=sites, flg_palindrom=flg_palindrom)
             
        self.unmodified_loci_in_GFF = self.dereplicate_sites([site for site in sites if site])
            
        # Stop program execution if the number of unmodified sites is too big.
        if self.max_sites_for_verification and len(self.unmodified_loci_in_GFF) > self.max_sites_for_verification:
            raise ValueError(f"Number of unmodified sites found, {len(self.unmodified_loci_in_GFF)}, is above the maximum number of sites to display, which is {self.max_sites_for_verification}!")
        
        oVerifiedSites = self.verify_sites(self.unmodified_loci_in_GFF)
        self.unmodified_sites += oVerifiedSites.get_entries()
        return copy.deepcopy(self.unmodified_loci_in_GFF)
    
    def get_unmodified_sites(self, modified_sites, seq, length, points, promoter, genes, expected_bases, flg_palindrom=False):
        '''
        As an example, points in GDgcHC are [4, -3]
        found sites = {'genome': '598_chr', 'method': 'kinModCall', 'modtype': 'm6A', 'start': 2397, 'end': 2408, 'score': '499', 
        'strand': '+', 'para': '.', 'data': {'coverage': '394', 'context': 'ATTATGTCTGACAATGAAGAAGACATTGATATCTTCTTTGC', 
        'IPDRatio': '4.94', 'identificationQv': '477'}, 'nucleotide': 'A', 'word': 'GAAGACATTGAT', 'match': 'Yes', 'site': 2399,
         
        'verified site': {'Motif': 'GAAGACATTGAT', 'Strand': 'DIR', 'Location': '2397..2408', 'Site': '2399', 'Match': 'Yes', 
        'Description': 'm6A\tDIR\t499\t394', 'Site location relative to TSC': '659', 'Annotation': "'K8B78_00010' [1740..2873]; 'dnaN'; 
        'DNA polymerase III subunit beta'", 'Motif strand': 'DIR', 'Absolute point': 3, 'Base': 'A', 
        'context': 'ATTATGTCTGACAATGAAGAAGACATTGATATCTTCTTTGC', 'nucleotide': 'A'}}
        '''
        
        # Convert found and verified modified site entries into a list of locations
        modified_site_locations = [int(item['verified site']['Site']) for item in  self.found_sites]
        
        # Create a list of all available locations found in the reference sequence
        available_sites = []
        
        # available_sites = [[site_location, point, start, end, site_strand, motif_strand],...]
        if flg_palindrom:
            for k in points:
                available_sites += [[int(start) + abs(k) - 1, k, int(start), int(end), "+" if k >= 0 else "-", strand] for 
                    start,end,strand in self.filtered_loci_in_refseq]
            # For palindroms, location on complement strand is dublication of the direct strand
            available_sites = [site for site in available_sites if site[-1] == "+"]
        else:
            # Filter by strand
            if self.strand:
                filtered_loci_in_refseq = [[start,end,strand] for start,end,strand in self.filtered_loci_in_refseq if
                    (strand == "+" and self.strand == 1) or (strand == "-" and self.strand == -1)]
            else:
                filtered_loci_in_refseq = self.filtered_loci_in_refseq
            for k in points:
                available_sites += [[
                    (int(start) + k - 1 if strand=="+" else int(end) - k + 1) if k > 0 else 
                        (int(end) + k + 1 if strand=="-" else int(start) + abs(k) - 1), 
                    k, int(start), int(end), 
                    "+" if ((strand=="+" and k >= 0) or (strand=="-" and k < 0)) else "-", strand] for 
                start,end,strand in filtered_loci_in_refseq]
                
        # Identify unmodified sites
        unmodified_loci_in_GFF = [ls for ls in available_sites if int(ls[0]) not in modified_site_locations]

        # Convert list of unmodified sites into a list of dictionary objects 
        unmodified_sites_in_GFF = [
            dict(zip(
                ['word', 'start', 'end', 'site', 'strand', 'score', 'modtype', 'match', 'data', 'point', 'motif_strand'],
                [tools.get_locus_by_coordinates(seq=seq, start=int(start), end=int(end)), 
                    int(start), int(end), site, site_strand, 0, "unmodified", "None", {'context': 
                        tools.get_locus_by_coordinates(seq=seq, start = int(site) - 20, end = int(site) + 20, 
                        flg_complement = True if (site_strand == "-") else False, 
                        flg_reverse = True if (site_strand == "-") else False).upper(), 
                        'coverage': 'na'}, point, motif_strand]
            )) for site,point,start,end,site_strand,motif_strand in unmodified_loci_in_GFF
        ]
        

        # Convert list of dictionary objects into a list of site objects
        sites = []
        for entry in unmodified_sites_in_GFF:
            site = self.get_unmodified_site(entry=entry, genes=genes, promoter_length=promoter, flg_palindrom=flg_palindrom)

            # Check whether unmodified base correspond to the expected base 
            if expected_bases[int(entry['point'])] != site['Base']:
                if not flg_palindrom:
                    tools.alert(f"Base identity control failed! In motif {site['Motif']} at location {entry['point']} " + 
                        f"base {expected_bases[int(entry['point'])]} is expected, but {site['Base']} is found.")
                continue
                
            sites.append(site)

        # Set run statistics
        self.get_run_statistics(sites=sites, flg_palindrom=flg_palindrom)
             
        # Remove empty site objects
        self.unmodified_loci_in_GFF = self.dereplicate_sites([site for site in sites if site])
            
        # Stop program execution if the number of unmodified sites is too big.
        if self.max_sites_for_verification and len(self.unmodified_loci_in_GFF) > self.max_sites_for_verification:
            raise ValueError(f"Number of unmodified sites found, {len(self.unmodified_loci_in_GFF)}, is above the maximum number of sites to display, which is {self.max_sites_for_verification}!")
        # Verify locations of unmodified sites
        oVerifiedSites = self.verify_sites(self.unmodified_loci_in_GFF)
        # If verified list is empty
        if not oVerifiedSites:
            return []
        # Create a list of entry objects of unmodified sites and return a copy of the list
        self.unmodified_sites += oVerifiedSites.get_entries()
        return copy.deepcopy(self.unmodified_loci_in_GFF)
    
    # This function is used only for unmodified palindroms
    # Set modifiednucleotide to site['Base']
    def reset_modified_base(self, site):
        # Check base
        start, end = [int(v) for v in (site['Location'].split("..") if site['Location'].find("..") > -1 else site['Location'].split("-"))]
        site_location = int(site['Site'])
        motif = site['Motif']
        if site_location > end:
            site['Base'] = ""
        else:
            site['Base'] = motif[site_location - start].upper()
        # Check strand assuming that palindrom motif can be only on DIR strand
        site['Motif strand'] = "DIR"
        site['Strand'] = "DIR" if site['Absolute point'] >= 0 else "REV"
    
    def get_available_loci(self):
        return [locus for locus in self.loci_in_refseq if not self.overlap_with_filtered_regions(int(locus[0]))]
        
    def get_all_motif_entries(self,seq):
        ls = []
        if not seq:
            return ls
        p = 0
        for w in self.loci_in_refseq:
            p = seq.find(w,p+1)
            if p < len(seq)/2:
                start = p+1
                end = p+len(w)
            else:
                start = len(seq)-p
                end = start + len(w) - 1
            if not self.overlap_with_filtered_regions(start):
                ls.append([start,end])
        return ls
    
    def get_description(self):
        return "%s,%s,%s" % (self.word,self.modbase_location,self.reverse_modbase_location)
    
    def select_overlaps(self, ls_items, points):
        def get_location(item):
            return item['Location'] if 'Location' in item else item['verified site']['Location']
            
        if not points:
            return []
            
        ls_items = sorted(reduce(lambda ls1, ls2: ls1 + ls2, copy.deepcopy(ls_items)), key=lambda item: [int(self.get_entry_field(item, 'Site')),self.get_entry_field(item, 'Location')])
        collection = [[copy.deepcopy(item), 1] for item in ls_items]
        
        for i in range(len(collection) - 1, 0, -1):
            if self.get_entry_field(collection[i - 1][0], 'Location') == self.get_entry_field(collection[i][0], 'Location'):
                # Check doublicates
                if self.get_entry_field(collection[i - 1][0], 'Site') != self.get_entry_field(collection[i][0], 'Site'):
                    collection[i - 1][1] += collection[i][1]
                del collection[i]
        
        collection = [ls for ls in collection if ls[1] >= len(points) and ls[1] % len(points) == 0]
        # Filter sites by locations 'start..end'
        locations = [self.get_entry_field(ls[0], 'Location') for ls in collection]
        ls_items = [item for item in ls_items if self.get_entry_field(item, 'Location') in locations]
        # Filter entries by site locations
        loci = [int(v) for v in [self.get_entry_field(ls[0], 'Site') for ls in collection]]
        if self.flg_show_modified_sites:
            entries = self.found_sites
        else:
            entries = self.unmodified_sites
        selected_entries = [entry for entry in entries if int(entry['site']) in loci]
        if self.flg_show_modified_sites:
            self.found_sites = selected_entries
        else:
            self.unmodified_sites = selected_entries
        
        return ls_items 
        
    def verify_sites(self, sites, gff_data=None):
        '''
        site = {'Motif': 'CTAG', 'Strand': 'DIR', 'Location': '1136743..1136746', 'Site': '1136743', 'Match': 'None', 
            'Description': 'unmodified\tDIR\t0\tna', 'Site location relative to TSC': '2067', 
            'Annotation': "'SIVHX_1221' [1134676..1136805]; 'Uncharacterized protein'"}
        gff_entry = {'genome': 'SVX1_Haloferax_genome', 'method': 'kinModCall', 'modtype': 'm4C', 
            'start': '33754', 'end': '33754', 'score': '38', 'strand': '+', 'para': '.', 
            'data': {'coverage': '52', 'context': 'GCCCCGCACGTTTCAACAACCTAGACAATTAATAGTGTTTA', 'IPDRatio': '1.85', 'identificationQv': '25'}, 'nucleotide': 'C'}
        '''
        # Utilities
        def _get_marker_string(entry):
            return f"{'DIR' if entry['strand'] in ('DIR','+','1') else 'REV'}:{entry['start']}..{entry['end']}"
            
        def _select_gff_records(ls_to_select, gff_data=[]): #ls_to_select = ['DIR:1..2000','REV:100..800',...]                
            return [entry for entry in gff_data if _get_marker_string(entry) in ls_to_select]
            
        def _get_strand(strand):
            return '+' if str(strand) in ('DIR','+','1') else '-'
            
        def _get_context(seq, site_location, strand):
            context = str(seq[int(site_location)-21:int(site_location)+20])
            #return context
            return context if strand == "+" else tools.reverse_complement(context)
        
        if gff_data:
            # Select records from GFF data
            gff_data = {'Heading':[], 'Body':_select_gff_records(ls_to_select = [f"{item['Strand']}:{item['Site']}..{item['Site']}" for 
                item in sites], gff_data = gff_data['Body'] if isinstance(gff_data,dict) else gff_data)}
        else:
            # Generate pseudo GFF records from site data
            seq = self.gbk['object'].seq if isinstance(self.gbk,dict) else self.gbk.seq
            gff_data = {'Heading':[],'Body':[]}
            gff_data['Body'] = [{'genome': '.', 'method': 'kinModCall', 'modtype': '.', 
                    'start': str(item['Site']), 'end': str(item['Site']), 'score': '0', 'strand': _get_strand(item['Strand']), 'para': '.', 
                    'data': {'coverage': '0', 'context': _get_context(seq = seq, site_location = int(item['Site']), strand = _get_strand(item['Strand'])), 
                        'IPDRatio': '0.0', 'identificationQv': '0'}, 
                        'nucleotide': seq[int(item['Site']) - 1] if _get_strand(item['Strand']) == "+" else tools.reverse_complement(seq[int(item['Site']) - 1])}
                for item in sites]
        
        if not len(gff_data['Body']):
            return []
            
        oMotifObj = Motif(reference=self.reference, binpath=self.binpath, tmppath=self.source_path,
            promoter_length=self.promoter_length, context_mismatches=0, motif_mismatch="No",
            filtered_loci=self.filtered_loci, strand=self.strand)
        success = oMotifObj.execute(gff_data = gff_data, gbk_data = self.gbk, echo=True)
        if success:
            return oMotifObj
        else:
            raise ValueError("Error with verification of selected loci!")
        
    def clean(self,path,generic_name):
        #fnames = list(filter(lambda fn: fn.find(generic_name) > -1, os.listdir(path)))
        fnames = [fn for fn in os.listdir(path) if fn.find(generic_name) > -1]
        for fname in fnames:
            os.remove(os.path.join(path,fname))
            
    def reset(self):
        self.output_table = []              # Output table
        self.loci_in_refseq = []            # Canonical motifs found in the reference sequence
        self.loci_in_GFF = []               # Canonical motifs found in GFF and passed through filterring
        self.filtered_loci_in_refseq = []   # Canonical motifs of the reference sequence passed through filtering
        self.modified_loci_in_GFF = []      # Canonical motifs overlapping with not masked loci of the reference sequence and modified
        self.unmodified_loci_in_GFF = []    # Canonical motifs overlapping with not masked loci of the reference sequence and unmodified
        self.masked_loci_in_refseq = []     # Canonical motifs found in GFF overlapping masked loci of the reference sequence
        self.masked_loci_in_GFF = []        # Canonical motifs overlapping with not masked loci of the reference sequence
        self.not_found_loci = []            # Modified motifs in GFF file, which were not found in the reference sequence
        self.num_instances = 0              # Number of canonical motifs found in GFF file - not found motifs
        self.num_found_sites = 0            # Number of found modified/unmodified canonical motifs
        self.num_expected_sites = 0         # Number of expected modified/unmodified sites
        self.num_found_motifs = 0           # Number of found motifs
        self.num_overlapped_motifs = 0      # Number of overlapped motifs sharing the same modified site
        self.oBlast = None
       
