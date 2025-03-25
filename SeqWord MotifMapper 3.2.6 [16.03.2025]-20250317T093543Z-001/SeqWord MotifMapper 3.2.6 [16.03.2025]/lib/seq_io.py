import os, sys, math, pickle, string, time, copy
import numpy as np
from datetime import datetime
from functools import reduce
import tools, container

##############################################################################################################
class IO:
    def __init__(self,path=""):
        self.oParser = Parser(path)
        
    def read(self,path,
                data_format="text",         # text | fasta | genbank | gff | binary
                inlist=False,               # used with 'text', split lines by \n
                separator="",               # used with 'text', split lines by \t
                strip_symbol="",            # used with 'text', strip given symbols
                reverse_complement=False,   # used with FASTA and GenBank to return reverse_complement sequences
                dbkey='$db$',               # used with binary, main key
                splkey='$suppl$',           # used with binary, supplementary key
                    ):
            if data_format.upper()=="TEXT":
                return self.open_text_file(path,inlist,separator,strip_symbol)
            if data_format.upper()=="BINARY":
                return self.open_binary_file(path,dbkey,splkey)
            if data_format.upper() in ("FASTA","FA"<"FST"):
                return self.read_seq(path,seq_format="fasta",reverse_complement=reverse_complement)
            if data_format.upper() in ("GENBANK","GBK","GB"):
                return self.read_seq(path,seq_format="genbank",reverse_complement=reverse_complement)
            if data_format.upper()=="GFF":
                return self.readGFF(path)
        
    def parse(self,path,
                data_format="text",         # text | fasta | genbank | gff | binary
                inlist=False,               # used with 'text', split lines by \n
                separator="",               # used with 'text', split lines by \t
                strip_symbol="",            # used with 'text', strip given symbols
                concatenate=False,          # used with FASTA and GenBank to return concatenated sequences
                reverse_complement=False,   # used with FASTA and GenBank to return reverse_complement sequences
                reverse_contigs=False,      # used with FASTA and GenBank to reverse order of contigs
                dbkey='$db$',               # used with binary, main key
                splkey='$suppl$',           # used with binary, supplementary key
                    ):
            if data_format.upper()=="TEXT":
                return self.open_text_file(self,path,inlist,separator,strip_symbol)
            if data_format.upper()=="BINARY":
                return self.open_binary_file(path,dbkey,splkey)
            if data_format.upper() in ("FASTA","FA"<"FST"):
                return self.parse_seq(path,seq_format="fasta",reverse_complement=reverse_complement,concatenate=concatenate,reverse_contigs=reverse_contigs)
            if data_format.upper() in ("GENBANK","GBK","GB"):
                return self.parse_seq(path,seq_format="genbank",reverse_complement=reverse_complement,concatenate=concatenate,reverse_contigs=reverse_contigs)
            if data_format.upper()=="GFF":
                return self.readGFF(path)
        
    #### Collection of save/open functions
    # fasta - FASTA formated collection of sequences    
    def save(self,data,path,
                data_format="text",     # text | binary
                dbkey='$db$',           # used with binary, main key
                splkey='$suppl$',       # used with binary, supplementary key
                suppl_data=None         # used with binary, supplementary data
                    ):
        if data_format=="binary":
            return self.save_binary_file(data,path,dbkey,suppl_data,splkey)
        with open(path, "w") as ofp:
            ofp.write(data)
            ofp.flush()
        return path
    
    def open_text_file(self,path,flg_inlist=False,sep="",strip_symbol=""):
        if not os.path.exists(path):
            return ""
        f = open(path)
        strText = f.read()
        f.close()
        strText = strText.replace("\"","")
        if flg_inlist:
            strText = strText.split("\n")
            if strip_symbol:
                #strText = list(map(lambda line: line.strip(strip_symbol), strText))
                strText = [line.strip(strip_symbol) for line in strText]
            if sep:
                #strText = list(map(lambda item: item.split(sep), strText))
                strText = [line.split(sep) for line in strText]
        return strText
    
    # Parse GFF file
    def readGFF(self, path, mode="object"): # mode = "object" | "dictionary"
        def get_entry(ls):
            # Convert GFF fields into dictionary fields
            entry = dict(zip(["genome","method","modtype","start","end","score","strand","para","data"],ls))
            data = [s.split("=") for s in entry["data"].split(";")]
            entry["data"] = dict(zip(map(lambda i: data[i][0], range(len(data))),map(lambda i: data[i][1], range(len(data)))))
            # Add modified nucleotide field
            context_sequence = entry['data']['context']
            modified_nucleotide = context_sequence[len(context_sequence)//2].upper()
            entry['nucleotide'] = modified_nucleotide
            return entry
        
        dGFF = {"Heading":[],"Body":[]}
        if not path or not os.path.exists(path):
            raise ValueError(f"")
            
        data = self.read(path,"text",inlist=True,separator="\t",strip_symbol=" ")
        for i in range(len(data)):
            if data[i][:2]=="##":
                dGFF["Heading"].append(data[i])
            elif len(data[i])==9:
                record = get_entry(data[i])
                dGFF["Body"].append(record)
            else:
                continue
        dGFF['Body'].sort(key=lambda d: [d['strand'],int(d['start'])])
        # Return GFF dictionary
        if mode == "dictionary":
            return dGFF
        # Create and return GFF object
        oGFF = GFF(path=path, heading=dGFF['Heading'], records=dGFF['Body'])
        return oGFF
    
    # Working with sequence files in fasta and genbank formats
    def parse_seq(self,path,seq_format="fasta",reverse_complement=False,concatenate=False,reverse_contigs=False):
        oParser = Parser()
        return oParser.parse(path,seq_format=seq_format,reverse_complement=reverse_complement,concatenate=concatenate,reverse_contigs=reverse_contigs)
    
    def read_seq(self,path,seq_format="fasta",reverse_complement=False):
        oParser = Parser()
        return oParser.read(path,seq_format=seq_format,reverse_complement=reverse_complement)
    
    # Working with binary files
    def openDBFile(self, fname, key='$db$', supplkey="$suppl$"):
        with open(fname, 'rb') as file:
            try:
                data = pickle.load(file)
                DB = data[key]
                supplementary = data[supplkey]
                return fname,DB,supplementary
            except:
                raise TypeError(f"File {fname} has wrong formatting or corrupted!")
    
    def saveDBFile(data, fname, supplementary=None, key='$db$', supplkey="$suppl$", flg_appendData=None):
        tmp_fname = os.path.join(os.path.dirname(fname),"~"+os.path.basename(fname))
        if os.path.exists(tmp_fname):
            os.remove(tmp_fname)
        # Open a file in binary write mode
        with open(tmp_fname, 'wb') as file:
            # Serialize and save the dictionary to the file
            pickle.dump({key:data,supplkey:supplementary}, file)
        if os.path.exists(fname):
            os.remove(fname)
        os.rename(tmp_fname,os.path.join(os.path.dirname(tmp_fname),os.path.basename(tmp_fname)[1:]))
        return fname
    
    # Create new folder
    def new_folder(self,folder_name):
        try:
            os.mkdir(folder_name)
            return folder_name
        except:
            return None
    
    # copy text files
    def copy(self,inpath,outpath,
                data_format="text",     # text | binary
                dbkey='$db$',           # used with binary, main key
                splkey='$suppl$',       # used with binary, supplementary key
                    ):
        if not os.path.exists(inpath):
            raise ValueError(f"Path {path} does not exist!")
        if data_format.upper()=="TEXT":
            try:
                self.save(read(inpath,data_format),outpath,data_format)
            except:
                raise TypeError(f"File {path} cannot be copied!")
        elif data_format.upper()=="BINARY":
            try:
                data,suppl_data = read(inpath,data_format,dbkey=dbkey,splkey=splkey)
                self.save(data,outpath,data_format,dbkey=dbkey,splkey=splkey,suppl_data=suppl_data)
            except:
                raise TypeError(f"File {path} cannot be copied!")
        return outpath
    

##############################################################################################################
# Collection of parsers
class Parser:
    def __init__(self, path="", seq_format=""): # seq_format = 'genbank','fasta'
        self.path = path
        self.seq_format = seq_format
    
    # Interface functions read and parse to accept user's requests and identify data type
    def parse(self,path="",seq_format="",reverse_complement=False,concatenate=False,reverse_contigs=False):  # Seq format can be predefined by users, or identified by file extension
        path,seq_format = self._check_path_format(path,seq_format)
        if seq_format.upper() == 'GENBANK':
            oObj = self._parse_genbank(path)
        elif seq_format.upper() == 'FASTA':
            oObj = self._parse_fasta(path)
        elif seq_format.upper() == 'GFF':
            return self._parse_gff(path)
        else:
            raise TypeError(f"Parsing of {seq_format} files is not supported!")
        if reverse_contigs:
            oObj = oObj.reverse()
        if reverse_complement:
            oObj = oObj.reverse_complement()
        if concatenate:
            oObj = oObj.get_concatenated()
        elif len(oObj)==0:
            raise ValueError(f"No sequences were found in file {path}")
        return oObj
                
    def read(self,path="",seq_format="",reverse_complement=False):  # Seq format can be predefined by users, or identified by file extension
        oObj = self.parse(path,seq_format=seq_format,reverse_complement=reverse_complement)
        if len(oObj) > 1:
            tools.msg(f"WARNING: file {path} contains more than one sequence!\nOnly first sequence is returned. Use 'parse' command to get all sequences of this file")
        return oObj[0]
            
    # Check input path and file format
    def _check_path_format(self,path="",seq_format=""):
        if not path:
            path = self.path
        if not seq_format:
            seq_format = self.seq_format
            if not seq_format:
                seq_format = self._define_seq_format(path)
        return path,seq_format
            
    # Seq formats are defined by file extensions
    def _define_seq_format(self,path):
        if path and path[path.rfind(".")+1:].upper() in ("GBK","GB","GBF"):
            return "GENBANK"
        if path and path[path.rfind(".")+1:].upper() in ("FASTA","FA","FST","FNA","FAA","FFN","FRN"):
            return "FASTA"
        return ""
        
    # Remove from text line Linux elements and quotes
    def format_line(self,line):
        return line.replace("\r","").replace("\"","'")
        
    # FASTA format parser
    def _parse_fasta(self,path):
        # Define the path to fasta file
        fasta_file = path
        
        # Initialize variables to store record information
        oFASTA = Fasta(os.path.basename(path[:path.rfind(".")]))
        # Open and read the GenBank file
        with open(fasta_file, "r") as file:
            lines = file.readlines()
            in_sequence = False
            
            lnum = 0
            while lnum < len(lines):
                line = self.format_line(lines[lnum])
                # Start capturing sequence data
                if line.startswith("ORIGIN"):
                    in_sequence = True
                    current_sequence = ""
                    
                # End sequence data capture
                elif in_sequence and line.startswith("//"):
                    in_sequence = False
                    # If there's a current sequence, add it to the list
                    if current_sequence:
                        oFASTA.append(current_sequence)
                
                # Capture sequence lines
                elif in_sequence:
                    current_sequence += line[10:].replace(" ", "").replace("\n", "").upper()
                lnum += 1
                        
        # Concatenate sequences with 50 N's in between
        oFASTA.concatenate()
        return oFASTA
    
    # GenBank format parser
    def _parse_genbank(self,path):
        # Define the path to GenBank file
        genbank_file = path
        
        # Create GBF as a collection of GBK records
        oGBF = GBF(os.path.basename(path[:path.rfind(".")]))
        # Open and read the GenBank file
        with open(genbank_file, "r") as file:
            # Initialize variables to store record information
            heading = []
            current_sequence = ""
            current_feature = None
            current_features = []
            lines = file.readlines()
            in_sequence = False
            in_features = False
            title = ""
            
            lnum = 0
            while lnum < len(lines):
                line = self.format_line(lines[lnum])
                # Capture heading lines
                if in_sequence == False and in_features == False:
                    heading.append(line.replace("\n",""))
                    if len(heading) == 1:
                        title = line[12:33].strip()
                        
                # Scip base count line
                if line.startswith("BASE COUNT"):
                    lnum += 1
                    continue
                
                # Start capturing sequence data
                if line.startswith("ORIGIN"):
                    in_sequence = True
                    in_features = False
                    current_sequence = ""
                    
                # End sequence data capture
                elif in_sequence and line.startswith("//"):
                    in_sequence = False
                    # Add new GBK record
                    oGBF.append(GBK(title,"\n".join(heading),SeqRecord(current_sequence,title),current_features))
                    current_features = []
                
                # Capture sequence lines
                elif in_sequence:
                    current_sequence += line[10:].replace(" ", "").replace("\n", "").upper()
                
                # Capture feature lines
                elif in_features==True:
                    if len(line) > 22 and line[21]=="/" and current_feature != None:
                        try:
                            key,value = line.strip().replace("\"","")[1:].split("=")
                        except:
                            lnum += 1
                            continue
                        # concatenate values in multiple lines
                        lnum += 1
                        line = lines[lnum]
                        while lnum < len(lines) and len(line) > 22 and line[5]==" " and line[21] != "/":
                            if key == "translation":
                                value += line.strip().replace("\"","")
                            else:
                                value += " " + line.strip().replace("\"","")
                            lnum += 1
                            line = lines[lnum]
                        current_feature.qualifiers[key] = [value.replace("\n","")]
                        continue
                    elif len(line) > 5 and line[5] != " ":
                        if current_feature:
                            current_features.append(current_feature)
                        feature_type = line[5:line.find(" ",5)]
                        strand = 1
                        line = line.replace("join(","").replace(")","")
                        if line.find("complement") > -1:
                            strand = -1
                            line = line.replace("complement(","")
                        exons = line[21:].strip().split(",")
                        current_feature = Feature(feature_type,strand,exons)
                        
                # Start capturing feature lines
                elif line.startswith("FEATURES"):
                    in_features = True
                
                # End feature data capture
                elif line.startswith("ORIGIN") or line.startswith("//"):
                    in_features = False
                    current_feature = None
                
                lnum += 1
        oGBF.concatenate()
        return oGBF
    
    # GFF file format parser
    def _parse_gff(self,path):
        # Define the path to GenBank file
        gff_file = path        
        # Create GFF as a collection of records representing each line in GFF file
        oGFF = GFF(os.path.basename(path[:path.rfind(".")]))
        # Open and read the GenBank file
        with open(gff_file, "r") as file:
            # Initialize variables to store record information
            heading = []
            in_features = False
            lines = file.readlines()
            lnum = 0
            while lnum < len(lines):
                line = self.format_line(lines[lnum])
                # Capture heading lines
                if line.startswith("##"):
                    heading.append(line)
                # Capture data lines as records
                elif len(line):
                    try:
                        entry = dict(zip(["genome","method","modtype","start","end","score","strand","para","data"],[s.strip() for s in line.split("\t")]))
                        data = list(map(lambda s: s.split("="), entry["data"].split(";")))
                        entry["data"] = dict(zip([data[i][0] for i in range(len(data))], [data[i][1] for i in range(len(data))]))
                        oGFF.append(Record(f"{entry['strand']}{entry['start']}..{entry['end']}",entry))
                    except:
                        raise ValueError(f"Line {line} cannot be parsed as a GFF record!")
                    
                else:
                    pass
                lnum += 1
        if heading:
            oGFF['Heading'] = "\n".join(heading)
        oGFF['Body'].sort(lambda d: [d.strand,int(d.start)])
        return oGFF
        
##############################################################################################################
class Fasta(container.Collection):
    def __init__(self,title="",spacer_length=50):
        container.Collection.__init__(self,title)
        self.description = self.title
        self.Seq = self.seq = None
        
    def __add__(self,other):
        if not isinstance(other,Fasta):
            raise TypeError("unsupported types of operands!")
        self.extend(other.get())
        self.concatenate()
        
    def __repr__(self):
        output = [f"{self.title} <{self.size} record(s)>"]
        for oGBK in self:
            output.append(f"\t>{oSeq.title}\n\t{str(oSeq)[:20]}...")
        return "\n".join(output)
        
    def base_count(self,nucleotides=[],concatenated=False):
        if nucleotides == []:
            nucleotides = ['a','c','g','t']
        if concatenated:
            return " ".join([f"{nuc} {self.Seq.base_count(nuc)}" for nuc in nucleotides])        
        return "\n".join([" ".join([f"{nuc} {oSeq.base_count(nuc)}" for nuc in nucleotides]) for oSeq in self])
                                    
    def get_concatenated(self):
        if self.Seq==None:
            self.concatenate()
        return self.Seq
            
    def concatenate(self):
        # Concatenate sequences with 50 N's in between
        if self.size > 1:
            sequences = [str(oSeq) for oSeq in self.get()]
            concatenated_sequence = self.spacer.join(sequences)
            self.Seq = self.seq = SeqRecord(concatenated_sequence,self.title)
        
        elif self.size:
            self.Seq = self.seq = self[0]
            
    def reverse_contigs(self):
        return self.reverse()
            
    def reverse_complement(self,reverse=False):
        # Create a copy
        seq_list = [oSeq.reverse_complement() for oSeq in self.get()]
        if reverse:
            seq_list.reverse()
        self.push(seq_list)
        self.concatenate()

    # return data as a formated text
    def format(self,concatenated_sequence=False):
        if concatenated_sequence:
            return self.Seq.format("fasta")
        return "\n".join([oSeq.format("fasta") for oSeq in self.get()])
    
    def copy(self,reverse_complement=False):
        oFasta = Fasta(self.title,len(self.spacer))
        oFasta.push(self.get())
        oFasta.concatenate()
        return oFasta
    
##############################################################################################################
class SeqRecord:
    def __init__(self,seq="",title="",moltype="DNA"):
        self.title = title
        self.description = self.title
        self.Seq = seq
        self.seq = self.Seq
        self.moltype = moltype
        
    def __len__(self):
        return len(self.Seq)
        
    def __str__(self):
        return self.Seq
        
    def __getitem__(self,index):
        if isinstance(index, slice):
            return SeqRecord(self.Seq[index.start:index.stop])
        return self.Seq[index]

    # Return GBK formated lines of sequence
    def _format_fasta(self,line_length = 100):
        fasta = [f">{self.title}"]
        fasta += [self.Seq[i:i + line_length].strip() if i <= len(self.Seq) - line_length else self.Seq[i:]
            for i in range(0, len(self.Seq), line_length)]
        return "\n".join(fasta)
        
    def _format_genbank(self,indend = 9,line_length = 66):
        seq = " ".join([self.Seq[i:i + 10] if i <= len(self.Seq) - 10 else self.Seq[i:] for i in range(0, len(self.Seq), 10)]).lower()
        seq = [seq[i:i + line_length].strip() if i <= len(seq) - line_length else seq[i:] for i in range(0, len(seq), line_length)]
        num = 1
        for i in range(len(seq)):
            length = len(seq[i].replace(" ",""))
            seq[i] = " "*(indend-len(str(num)))+str(num)+" " + seq[i]
            num += length
        return seq
        
    def format(self,seq_format="FASTA"):
        if seq_format.upper()=="FASTA":
            return self._format_fasta()
        elif seq_format.upper()=="GENBANK":
            return self._format_genbank()
            
    def base_count(self,nuc):
        return self.Seq.upper().count(nuc.upper())
            
    def reverse_complement(self,seq):
        # Create a translation table for nucleotides, including ambiguous ones.
        # This uses the str.maketrans() function to map each nucleotide to its complement.
        complement_table = str.maketrans(
            "ATGCRYWSKMBDHVNatgcrywskmbdhvn",
            "TACGYRWSMKVHDBNatgcyrwsmkvhdbn"
        )        
        # Translate the sequence using the translation table, then reverse it
        return seq.translate(complement_table)[::-1]    
    
    def translate(self,dna_seq="",strand=1):
        if not dna_seq and self.moltype=="DNA":
            dna_seq = self.Seq
        elif isinstance(dna_seq,SeqRecord) and dna_seq.moltype=="DNA":
            dna_seq = str(dna_seq.Seq)
        elif isinstance(dna_seq,str):
            pass
        else:
            raise TypeError(f"Object {dna_seq} cannot be translated!")
        if strand == -1:
            dna_seq = self.reverse_complement(dna_seq)
        # Bacterial codon table 11 (standard)
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        protein_seq = ''
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            # Translate the codon to an amino acid and append to the protein sequence
            protein_seq += codon_table.get(codon.upper(), 'X')  # 'X' for unknown codons
        return SeqRecord(protein_seq,self.title,"protein")
    
    def copy(self,reverse_complement=False):
        if reverse_complement:
            return SeqRecord(self.reverse_complement(self.Seq),self.title)
        return SeqRecord(self.Seq,self.title,self.moltype)
    
##############################################################################################################
class Location:
    def __init__(self,exons=["0..0"]):
        self.exons = list(map(lambda s: list(map(lambda v: int(v.replace(">","").replace("<","")), s.split(".."))), exons))
        #self.exons = [[int(v.replace(">", "").replace("<", "")) for v in s.split("..")] for s in exons]
        #self.exons = [ exons]
        self.start = self.exons[0][0]
        self.end = self.exons[-1][1]
        
    # Format GBK feature location
    def get_exons(self,start=0,complement=0):
        # 'start' will be added to start and locations
        # if 'complement' > 0, the valu means the total sequence length
        if not complement:
            return [f"{locus[0] + start}..{locus[1] + start}" for locus in self.exons]        
        return [f"{complement - locus[1] + start}..{complement - locus[0] + start}" for locus in self.exons]
            
    def format(self):
        exons = self.get_exons()
        if len(exons) > 1:
            return f"join({','.join(exons)})"
        return ','.join(exons)
        
    def copy(self,start=0,complement=0):
        return Location(self.get_exons(start=start,complement=complement))

##############################################################################################################
class GBF(container.Collection):
    def __init__(self,title="",spacer_length=50):
        container.Collection.__init__(self,title)
        self.description = self.title
        self.GBK = None
        self.spacer = "N" * spacer_length
        
    def __add__(self,other):
        if not isinstance(other,GBF):
            raise TypeError("unsupported types of operands!")
        self.extend(other.get())
        self.concatenate()
        
    def __repr__(self):
        output = [f"{self.title} <{self.size} record(s)>"]
        for oGBK in self:
            output.append(f"\t>{oGBK.title}\n\t{str(oGBK.Seq)[:20]}...")
        return "\n".join(output)
        
    def base_count(self,nucleotides=[],concatenated=False):
        if nucleotides == []:
            nucleotides = ['a','c','g','t']
        if concatenated:
            return " ".join([f"{nuc} {self.GBK.Seq.base_count(nuc)}" for nuc in nucleotides])
        return "\n".join([" ".join([f"{nuc} {oGBK.Seq.base_count(nuc)}" for nuc in nucleotides]) for oGBK in self])
                                    
    def concatenate(self):
        # Concatenate sequences with 50 N's separators
        if len(self) > 1:
            sequences = [str(oGBK.Seq) for oGBK in self]
            concatenated_sequence = self.spacer.join(sequences)
            features = []
            contig_start = 0
            contig_counter = 1
            for i in range(len(sequences)):
                sequence = sequences[i]
                contig_end = contig_start + len(sequence)
                # Add contig features
                features.append([Feature("contig",1,[f"{contig_start}..{contig_end}"])]+self[i].get_features(contig_start))
                features[-1][0].qualifiers["note"] = [f"Contig_{contig_counter}"]
                contig_start = contig_end + 1
                contig_counter += 1
            self.GBK = GBK(self.title,"",SeqRecord(concatenated_sequence,self.title),features)
        
        elif len(self) == 1:
            self.GBK = self[0]
            
    def get_concatenated(self):
        if self.GBK==None:
            self.concatenate()
        return self.GBK
            
    def reverse_complement(self,reverse=False):
        # Create a copy
        gbk_list = [oGBK.reverse_complement() for oGBK in self]
        if reverse:
            gbk_list.reverse()
        self.push(gbk_list)
        self.concatenate()

    def reverse_contigs(self):
        return self.reverse()
            
    # return data as a formated text
    def format(self,seq_format="GENBANK",concatenated_sequence=False):
        if concatenated_sequence:
            if seq_format.upper()=="GENBANK":
                return self.GBK._format_genbank()
            elif seq_format.upper()=="FASTA":
                return self.GBK._format_fasta()
            else:
                tools.msg("Wrong file format %s!" % seq_format)
                return ""
        if seq_format.upper()=="GENBANK":
            #return "\n".join(list(map(lambda oGBK: oGBK._format_genbank(), self.get())))
            return "\n".join([oGBK._format_genbank() for oGBK in self.get()])
        elif seq_format.upper()=="FASTA":
            #return "\n".join(list(map(lambda oGBK: oGBK._format_fasta(), self.get())))
            return "\n".join([oGBK._format_fasta() for oGBK in self.get()])
        else:
            tools.msg("Wrong file format %s!" % seq_format)
            return ""
    
    def copy(self,reverse_complement=False):
        oGBF = GBF(self.title,len(self.spacer))
        oGBF.push(self.get())
        oGBF.concatenate()
        return oGBF

##############################################################################################################
class GBK:
    def __init__(self,title="",heading="",oSeq=None,features=[]):
        self.title = title
        self.description = self.title
        self.heading = heading
        self.description = self.title
        self.accession = self.title
        self.Seq = self.seq = oSeq
        self.features = features
                
    def __len__(self):
        return len(self.Seq)

    def __repr__(self):
        return f"Title: {self.title}\nSeq. Length: {len(self.Seq)}\nNum. Features: {len(self.features)}"
        
    def __str__(self):
        return self.__repr__()
        
    def __getitem__(self,index):
        if isinstance(index, slice):
            oCopy = self.copy()
            # Set heading and title
            oCopy.heading = ""
            oCopy.title = self.title+"slice"
            oCopy.description = oCopy.title
            # Seq concatenated sequence
            oCopy.Seq += self.Seq[index.start:index.end]
            # Add features from the second GBK after an adjastment of their locations
            oCopy.features += [
                ft.copy(start=-index.start) 
                for ft in other.features 
                if ft.location.start >= index.start and ft.location.end <= index.stop
            ]
            return oCopy
        else:
            return self.Seq[index]
        
    def __add__(self,other):
        if not isinstance(other, GBK):
            raise TypeError("Both operands must be 'GBK' objects")
        # Create a copy
        oCopy = self.copy()
        # Set heading and title
        oCopy.heading = ""
        oCopy.title = self.title+"_concat"
        oCopy.description = oCopy.title
        # Seq concatenated sequence
        oCopy.Seq += other.Seq
        oCopy.seq = oCopy.Seq
        # Add features from the second GBK after an adjastment of their locations
        oCopy.features += [ft.copy(start=len(str(self.Seq))) for ft in other.features]
        return oCopy
    
    # Simmilar to addition, but separate sequences with N's
    def concatenate(self,other,separator_length=50):
        if not isinstance(other, GBK):
            raise TypeError("Both operands must be 'GBK' objects")
        # Create a copy
        oCopy = self.copy()
        # Set heading and title
        oCopy.heading = ""
        oCopy.title = self.title+"_concat"
        oCopy.description = oCopy.title
        # Seq concatenated sequence
        oCopy.Seq += "N"*separator_length + other.Seq
        oCopy.seq = oCopy.Seq
        # Add features from the second GBK after an adjastment of their locations
        oCopy.features += [ft.copy(start=len(str(self.Seq)))+separator_length for ft in other.features]
        return oCopy
    
    def rearrange(self,moltype="CDS",tag="",gene="",product="",strand=1):
        shift = 100
        # Create a copy
        oCopy = self.copy()
        # List of features
        #features = list(map(lambda ft: ft.copy(), oCopy.features))
        features = [ft.copy() for ft in oCopy.features]
        if features and moltype:
            features = [ft for ft in features if ft.type == moltype]        
        if features and tag:
            features = [ft for ft in features if 'locus_tag' in ft.qualifiers and ft.qualifiers['locus_tag'] and ft.qualifiers['locus_tag'][0] == tag]        
        if features and gene:
            features = [ft for ft in features if 'gene' in ft.qualifiers and ft.qualifiers['gene'] and ft.qualifiers['gene'][0] == gene]        
        if features and product:
            features = [ft for ft in features if 'product' in ft.qualifiers and ft.qualifiers['product'] and ft.qualifiers['product'][0] == product]
        if not len(feature):
            return None
            
        ft = features[0]
        reverse_complement = 0
        if strand != ft.strand:
            reverse_complement = len(self.Seq)
        if reverse_complement:
            if ft.location.end+shift >= len(self.Seq):
                shift = len(self.Seq)-ft.location.end-1
            oCopy = oCopy[ft.location.end+shift:]+oCopy[:ft.location.end+shift]
            return oCopy.reverse_complement()
        if ft.location.start-shift < 0:
            shift = ft.location.start
        return oCopy[ft.location.start-shift:]+oCopy[:ft.location.start-shift]
        
    def get_features(self,increment=0):
        return [ft.copy(increment) for ft in self.features]
        
    def reverse_complement(self):
        # Create a copy
        return self.copy(True)

    # return data as a formated text
    def format(self,seq_format="GENBANK"):
        if seq_format.upper()=="GENBANK":
            return self._format_genbank()
        elif seq_format.upper()=="FASTA":
            return self._format_fasta()
        else:
            tools.msg("Wrong file format %s!" % seq_format)
            return ""
            
    # Get features in fasta format
    def features2fasta(self,moltype="CDS",seqtype="PROT"):
        if moltype:
            features = [ft for ft in self.features if ft.type==moltype]
        else:
            features = self.features
        if moltype.upper() == "DNA":
            return "\n".join([f">{oFT.get_feature_title()}\n{str(self.Seq[oFT.start - 1:oFT.end])}" for oFT in self])
        elif moltype.upper() in ("PROT","PROTEIN","AMC"):
            fasta = []
            for oFT in features:
                if 'translation' in oFT.qualifiers:
                    fasta.append(f">{oFT.get_feature_title()}\n{oFT.qualifiers['translation'][0]}")
                else:
                    fasta.append(f">{oFT.get_feature_title()}\n{self.Seq.translate(str(self.Seq[oFT.start-1:oFT.end]),oFT.strand)}")
            return "\n".join(fasta)
        else:
            raise ValueError(f"Sequence type {moltype} is not supported!")
    
    def _format_genbank(self):
        # Get first GBK line
        gbk = [self._get_first_line()]
        # Format heading
        heading = self.heading.split("\n")
        if len(heading) > 1:
            gbk += heading[1:]
        # Format features and sequence
        gbk += (
            sum([ft.format() for ft in self.features], []) +
            ["BASE COUNT   " + " ".join([f"{nuc} {self.Seq.base_count(nuc)}" for nuc in ['a', 'c', 'g', 't']]), "ORIGIN"] +
            self.Seq.format("genbank") +
            ["//", ""]
        )            
        return "\n".join(gbk)
        
    def _format_fasta(self):
        return self.Seq.format("fasta")
        
    def _get_first_line(self):
        # Format title
        title = self.title
        if len(title) > 20:
            title = title[:17]+"..."
        # Get the current date
        current_date = datetime.now()        
        # Format the date as "UNK DD-MMM-YYYY"
        formatted_date = "UNK " + current_date.strftime("%d-%b-%Y").upper()
        # Format first line of GBK file
        return ("LOCUS       %s%s%d bp%sDNA              %s" %
                (title," "*(21-len(title)),len(str(self.Seq))," "*(11-len(self.Seq)),formatted_date))
                
    def copy(self, reverse_complement=False):
        if reverse_complement:
            return GBK(
                self.title,
                self.heading,
                self.seq.copy(reverse_complement=True),
                [ft.copy(complement=len(self.Seq)) for ft in self.features]
            )
        
        return GBK(
            self.title,
            self.heading,
            self.Seq.copy(),
            [ft.copy() for ft in self.features]
        )

##############################################################################################################
class Feature:
    def __init__(self,ftype="",strand=0,exons=["0..0"]):
        self.type = ftype
        self.strand = strand
        self.location = Location(exons)
        self.qualifiers = {}
    
    def __repr__(self):
        return f"Type: {self.type}\tLocation: {self.location.start}..{self.location.end}; strand: {self.strand}"
        
    def __str__(self):
        return self.__repr__()
        
    # Generate feature title
    def _get_feature_title(self):
        feature_title = self.type+":"
        if 'locus_tag' in self.qualifiers:
            feature_title += f" [{self.qualifiers['locus_tag'][0]}]"
        if 'gene' in self.qualifiers:
            feature_title += f" ({self.qualifiers['gene'][0]})"
        if 'product' in self.qualifiers:
            feature_title += f" {self.qualifiers['product'][0]}"
        return feature_title

    # Return GBK formated qualifier data
    def _format_qualifier(self,key):
        output = []
        if key not in list(self.qualifiers.keys()):
            return output
        for i in range(len(self.qualifiers[key])):
            output += self._format_feature_text(key,self.qualifiers[key][i])
        return output
        
    # Return GBK formated lines of feature data
    def _format_feature_text(self,key,value):
        field_length = 80
        indend = 21
        length = field_length-indend
        output = [" "*indend]
        line = ("/"+key+"=\""+value+"\"").split(" ")
        while line:
            while len(output[-1]) < field_length:
                text_element = line[0]
                if len(text_element) > length:
                    line[0] = line[0][field_length-len(output[-1])-1:]
                    text_element = text_element[:field_length-len(output[-1])]
                elif len(output[-1])+len(text_element) > field_length:
                    break
                else:
                    line = line[1:]
                output[-1] += text_element+" "
                if not line:
                    break
            output.append(" "*indend)
        output = [s for s in output if len(s.replace(" ","")) > 1]
        return output
    
    # Format GBK feature location
    def _format_location(self):
        # Check strand of the feature
        if self.strand == 1:
            return self.location.format()
        else:
            return f"complement({self.location.format()})"
            
    # Return GBK formated feature
    def format(self):
        # Format feature heading
        output = [" "*5 + self.type + " "*(16-len(self.type)) + self._format_location()]
        # Add qualifiers
        for key in list(self.qualifiers.keys()):
            output += self._format_qualifier(key)
        return output
        
    # Return a unique identifier
    def get_tag(self,index=0):
        if 'locus_tag' in self.qualifiers:
            return self.qualifiers['locus_tag'][0]
        return f"{self.type}_{index}"
        
    def copy(self,start=0,complement=0): 
        # 'start' will be added to start and locations
        # if 'complement' > 0, the valu means the total sequence length
        oCopy = Feature(self.type, self.strand, [exon for exon in self.location.get_exons(start=start, complement=complement)])
        oCopy.location = self.location.copy(start=start,complement=complement)
        oCopy.qualifiers = copy.deepcopy(self.qualifiers)
        return oCopy
        
##############################################################################################################
class GFF(container.Collection):
    def __init__(self, path="", heading="", records=[]):
        container.Collection.__init__(self, path)
        self.path = self.title
        self.heading = self.Heading = heading
        self.records = self.body = self.Body = container.Container()
        for record in records:
            self.records.append(GFF_record(**record))
        
    def __str__(self):
        return f"File: {self.path}\n{len(self.body)} records"
        
    def copy(self):
        oGffCopy = GFF(path=self.title)
        oGffCopy.heading = self.heading
        oGffCopy.body = self.body.copy()
        return oGffCopy
        
##############################################################################################################
class Record:
    def __init__(self, **attributes):

        # Parse attributes
        for key, value in list(attributes.items()):
            if isinstance(value, dict):
                setattr(self, key, Record(**value))
            else:
                setattr(self, key, value)
                    
    def __getitem__(self,attr):
        if isinstance(attr, str) and hasattr(self,attr):
            return eval(compile(f"self.{index_start}", "<string>", "eval"))
        raise ValueError(f"Object {type(self)} has no attribute {str(attr)}!")
    
    def __setitem__(self,attr,value):
        if isinstance(attr, str) and hasattr(self,attr):
            setattr(self, attr, value)
        raise ValueError(f"Object {type(self)} has no attribute {str(attr)}!")
        
##############################################################################################################

class GFF_record(Record):
    def __init__(self,  **attributes):
        Record.__init__(self, **attributes)

##############################################################################################################
if __name__ == "__main__":
    oSeqIO = IO()
    path = os.path.join("..","input","S.aureus_150.gbk")
    oGBK = oSeqIO.read(path,"genbank")
    oSeqIO.save(oGBK.format("genbank"),"S.aureus_150.gbk")
    oSeqIO.save(oGBK.format("fasta"),"S.aureus_150.fa")
    
