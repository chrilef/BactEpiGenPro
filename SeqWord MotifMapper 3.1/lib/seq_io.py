import os, sys, math, shelve, string, time, copy
import numpy as np
from datetime import datetime
from functools import reduce
import tools

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
        ofp = open(path, "w")
        ofp.write(data)
        ofp.flush()
        ofp.close()
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
                strText = list(map(lambda line: line.strip(strip_symbol), strText))
            if sep:
                strText = list(map(lambda item: item.split(sep), strText))
        return strText
    
    # Parse GFF file
    def readGFF(self,path):
        def get_entry(ls):
            entry = dict(zip(["genome","method","modtype","start","end","score","strand","para","data"],ls))
            data = list(map(lambda s: s.split("="), entry["data"].split(";")))
            entry["data"] = dict(zip(map(lambda i: data[i][0], range(len(data))),map(lambda i: data[i][1], range(len(data)))))
            return entry
        
        GFF = {"Heading":[],"Body":[]}
        if not os.path.exists(path):
            return GFF
        data = self.read(path,"text",inlist=True,separator="\t",strip_symbol=" ")
        for i in range(len(data)):
            if type(data[i])==type("s") and len(data[i]) > 2 and data[i][:2]=="##":
                GFF["Heading"].append(data[i])
            elif len(data[i])==9:
                GFF["Body"].append(get_entry(data[i]))
            else:
                continue
        GFF['Body'].sort(key=lambda d: [d['strand'],int(d['start'])])
        return GFF
    
    '''
    def readGFF(self,path):
        print("seq_io:105")
        oParser = Parser()
        return oParser.read(path,seq_format="GFF")
    '''
        
    # Working with sequence files in fasta and genbank formats
    def parse_seq(self,path,seq_format="fasta",reverse_complement=False,concatenate=False,reverse_contigs=False):
        oParser = Parser()
        return oParser.parse(path,seq_format=seq_format,reverse_complement=reverse_complement,concatenate=concatenate,reverse_contigs=reverse_contigs)
    
    def read_seq(self,path,seq_format="fasta",reverse_complement=False):
        oParser = Parser()
        return oParser.read(path,seq_format=seq_format,reverse_complement=reverse_complement)
    
    # Working with binary files
    def save_binary_file(self,data,fname,dbkey='data',suppl_data=None,splkey='$suppl$'):
        f = shelve.open(fname)
        f[dbkey] = data
        f[splkey] = suppl_data
        f.close()

    def open_binary_file(self,fname=None, dbkey='$db$', splkey='$suppl$'):
        if not fname or not os.path.exists(fname):
            return
        try:
            f = shelve.open(fname)
            self.oParser = {dbkey:f[dbkey],splkey:{}}
            if f.has_key(splkey):
                self.oParser[splkey].update(f[splkey])
            f.close()
            return fname,self.oParser[dbkey],self.oParser[splkey]
        except:
            return None

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
        if len(oObj) > 1:
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
        return line.replace("\r","").replace("\l","").replace("\"","'")
        
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
            print("parser:207",len(lines))
            while lnum < len(lines):
                line = self.format_line(lines[lnum])
                # Capture heading lines
                if line.startswith("##"):
                    heading.append(line)
                # Capture data lines as records
                elif len(line):
                    entry = dict(zip(["genome","method","modtype","start","end","score","strand","para","data"],list(map(lambda s: s.strip(),line.split("\t")))))
                    data = list(map(lambda s: s.split("="), entry["data"].split(";")))
                    entry["data"] = dict(zip(map(lambda i: data[i][0], range(len(data))),map(lambda i: data[i][1], range(len(data)))))
                    oGFF.append(Record(f"{entry['strand']}{entry['start']}..{entry['end']}",**entry))
                    '''
                    try:
                        entry = dict(zip(["genome","method","modtype","start","end","score","strand","para","data"],list(map(lambda s: s.strip(),line.split("\t")))))
                        data = list(map(lambda s: s.split("="), entry["data"].split(";")))
                        entry["data"] = dict(zip(map(lambda i: data[i][0], range(len(data))),map(lambda i: data[i][1], range(len(data)))))
                        oGFF.append(Record(f"{entry['strand']}{entry['start']}..{entry['end']}",entry))
                    except:
                        raise ValueError(f"Line {line} cannot be parsed as a GFF record!")
                    '''
                else:
                    pass
                lnum += 1
        if heading:
            oGFF['Heading'] = "\n".join(heading)
        oGFF['Body'].sort("lambda d: [d.strand],int(d.start)]")
        return oGFF
        
########################################################################
class ContainerElement:
    def __init__(self,Obj,index,link=None):
        if not hasattr(Obj,'title'):
            Obj.title = str(index)
        self.title = Obj.title
        self.index = index
        self.link = link
        self.Obj = Obj
        
    def link_function(self,command,args=[]):
        if command == "increment":
            self.index += 1
        elif command == "decrement":
            self.index -= 1
        else:
            pass
            
########################################################################
class Container:
    def __init__(self,title=""):
        self.title = title
        self.container = np.array([], dtype=object)
        self.titles = {}
        self.size = 0
    
    def __add__(self,other):
        if not isinstance(other,Collection):
            raise TypeError("unsupported types of operands!")
        return self.copy().extend(other.get())
        
    def __sub__(self,other):
        if not isinstance(other,Collection):
            raise TypeError("unsupported types of operands!")
        oCopy = self.copy()
        for key in other.get_titles():
            if key in list(self.titles.keys()):
                oCopy.__delitem__(key)
        return oCopy
        
    def __len__(self):
        return self.size
        
    def __contains__(self,key):
        if isinstance(key,str):
            if not self.has(key):
                return False 
        elif isinstance(key,int):
            if key > self.size:
                return False
        else:
            raise ValueError("Container must be refered either by text title or integer index!")
        return True
    
    def __iter__(self):
        if not len(self.container):
            return iter([])
        return iter(list(map(lambda obj: obj.Obj, self.container)))
    
    def __getitem__(self,index):
        if isinstance(index, slice):
            return self.get(index.start,index.stop)
        return self.get(index)
    
    def __setitem__(self,key,Obj):
        index,title = self._parse_key(key)
        if title=="~attribute":
            setattr(self, key, Obj)
            return
        link = None
        if index:
            link = self.container[index-1].Obj.link
        Obj = ContainerElement(Obj,index,link)
        self.titles[title] = Obj
        self.container[index] = Obj
    
    def __delitem__(self,key):
        index,title = self._parse_key(key)
        self.container = np.delete(self.container,index)
        del self.titles[title]
        self.size -= 1
        
    def _parse_key(self,key):
        if not key:
            key = 0
        if isinstance(key,str):
            if not self.has(key):
                if hasattr(self,key):
                    return key,"~attribute"
                raise ValueError(f"There is no object with title {key}!") 
            index = self.titles[key].index
            title = key
        elif isinstance(key,int):
            index = key
            title = self.container[key].title
        else:
            raise ValueError("Container must be refered either by text title or integer index!")
        return index,title
    
    def _get_title(self,obj):
        if hasattr(obj, 'title'):
            return obj.title
        return ""

    def _check_object_title(self,Obj):
        if not hasattr(Obj,'title') or not Obj.title:
            Obj.title = str(self.size)
        if self.has(Obj.title):
            raise KeyError(f"Title {Obj.title} is already occupied!")
        return True
            
    def clear(self):
        self.container = np.array([], dtype=object)
        self.titles = {}
        self.size = 0

    def has(self,title):
        return title in list(self.titles.keys())
    
    def index(self,title):
        if not self.has(title):
            return -1
        return self.titles[title].Obj.index
    
    def append(self,Obj):
        self._check_object_title(Obj)
        link = None
        if self.size:
            link = self.container[self.size-1].link
        Obj = ContainerElement(Obj,self.size,link)
        self.titles[Obj.title] = Obj
        self.container = np.append(self.container, Obj)
        self.size += 1
        
    def insert(self,key,Obj):
        index,title = self._parse_key(key)
        if abs(index) >= self.size:
            self.container.append(Obj)
            return
        self._check_object_title(Obj)
        Obj = ContainerElement(Obj,index,None)
        if index == 0:
            self.container[0].link = Obj.link_function
        else:
            Obj.link = self.container[index-1].link_function
            self.container[index].link = Obj.link_function
        self.container = np.insert(self.container, index, Obj)
        self.titles[Obj.title] = Obj
        self.size += 1            
        
    def extend(self,ls):
        for Obj in ls:
            if hasattr(Obj,'title') and self.has(Obj.title):
                continue
            self.append(Obj)
        
    def get_titles(self):
        return list(self.titles.keys())
    
    def get(self,start_key=None,stop_key=None):
        if start_key==None and stop_key==None:
            return list(map(lambda obj: obj.Obj.copy(), self.container))
        index_start,title_start = self._parse_key(start_key)
        if title_start=="~attribute":
            return eval(compile(f"self.{index_start}", "<string>", "eval"))
        if stop_key:
            index_stop,title_stop = self._parse_key(stop_key)
            indices = [index_start,index_stop]
            indices.sort()
            return list(map(lambda i: self.container[i].Obj, range(indices[0],indices[1]+1,1)))
        else:
            return self.container[index_start].Obj
            
    def dict(self):
        return self.titles
        
    def push(self,obj_ls):
        if self.size and obj_ls and type(self[0]) != type(obj_ls[0]):
            raise TypeError("Object in the collection and the list are of different types!")
        self.clear()
        for Obj in obj_ls:
            self.append(Obj)
            
    def sort(self,key="",reverse=False):
        obj_list = self.get()
        self.clear()
        # if key==None, sorting is escaped
        if key != None:
            if key=="":
                obj_list.sort(key=lambda Obj: Obj.title)
            else:
                try:
                    eval(compile(f"obj_list.sort(key={key})", "<string>", "eval"))
                except:
                    raise TypeError(f"command {key} cannot be used for sorting these objects!")
        if reverse:
            obj_list.reverse()
        for Obj in obj_list:
            self.append(Obj)
            
    def sorted(self,key="",reverse=False):
        return self.copy().sort(key,reverse)
        
    def reverse(self):
        return self.sorted(key=None,reverse=True)
            
    def copy(self):
        oNewContainer = Collection(self.title)
        for i in range(self.size):
            try:
                oNewContainer.append(self.container[i].Obj.copy())
            except:
                oNewContainer.append(copy.deepcopy(self.container[i].Obj))
        return oNewContainer

##############################################################################################################
class Collection(Container):
    def __init__(self,title=""):
        Container.__init__(self,title)

##############################################################################################################
class Fasta(Collection):
    def __init__(self,title="",spacer_length=50):
        Collection.__init__(self,title)
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
            return " ".join(list(map(lambda nuc: f"{nuc} {self.Seq.base_count(nuc)}", nucleotides)))
        return "\n".join(list(map(lambda oSeq: " ".join(list(map(lambda nuc: f"{nuc} {oSeq.base_count(nuc)}", nucleotides))), self)))
                                
    def get_concatenated(self):
        if self.Seq==None:
            self.concatenate()
        return self.Seq
            
    def concatenate(self):
        # Concatenate sequences with 50 N's in between
        if self.size > 1:
            sequences = list(map(lambda oSeq: str(oSeq), self))
            concatenated_sequence = self.spacer.join(sequences)
            self.Seq = self.seq = SeqRecord(concatenated_sequence,self.title)
        
        elif self.size:
            self.Seq = self.seq = self[0]
            
    def reverse_contigs(self):
        return self.reverse()
            
    def reverse_complement(self,reverse=False):
        # Create a copy
        seq_list = list(map(lambda oSeq: oSeq.reverse_complement(), self))
        if reverse:
            seq_list.reverse()
        self.push(seq_list)
        self.concatenate()

    # return data as a formated text
    def format(self,concatenated_sequence=False):
        if concatenated_sequence:
            return self.Seq.format("fasta")
        return "\n".join(list(map(lambda oSeq: oSeq.format("fasta"), self.get())))
    
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
        fasta += list(map(lambda i: self.Seq[i:i+line_length].strip() if i <= len(self.Seq)-line_length else self.Seq[i:], range(0,len(self.Seq),line_length)))
        return "\n".join(fasta)
        
    def _format_genbank(self,indend = 9,line_length = 66):
        seq = " ".join(list(map(lambda i: self.Seq[i:i+10] if i <= len(self.Seq)-10 else self.Seq[i:], range(0,len(self.Seq),10)))).lower()
        seq = list(map(lambda i: seq[i:i+line_length].strip() if i <= len(seq)-line_length else seq[i:], range(0,len(seq),line_length)))
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
        self.start = self.exons[0][0]
        self.end = self.exons[-1][1]
        
    # Format GBK feature location
    def get_exons(self,start=0,complement=0):
        # 'start' will be added to start and locations
        # if 'complement' > 0, the valu means the total sequence length
        if not complement:
            return list(map(lambda locus: f"{locus[0]+start}..{locus[1]+start}", self.exons))
        return list(map(lambda locus: f"{complement-locus[1]+start}..{complement-locus[0]+start}", self.exons))
        
    def format(self):
        exons = self.get_exons()
        if len(exons) > 1:
            return f"join({','.join(exons)})"
        return ','.join(exons)
        
    def copy(self,start=0,complement=0):
        return Location(self.get_exons(start=start,complement=complement))

##############################################################################################################
class GBF(Collection):
    def __init__(self,title="",spacer_length=50):
        Collection.__init__(self,title)
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
            return " ".join(list(map(lambda nuc: f"{nuc} {self.GBK.Seq.base_count(nuc)}", nucleotides)))
        return "\n".join(list(map(lambda GBK: " ".join(list(map(lambda nuc: f"{nuc} {oGBK.Seq.base_count(nuc)}", nucleotides))), self)))
                                
    def concatenate(self):
        # Concatenate sequences with 50 N's in between
        if self.size > 1:
            sequences = list(map(lambda oGBK: str(oGBK.Seq), self))
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
        
        elif self.size:
            self.GBK = self[0]
            
    def get_concatenated(self):
        if self.GBK==None:
            self.concatenate()
        return self.GBK
            
    def reverse_complement(self,reverse=False):
        # Create a copy
        gbk_list = list(map(lambda oGBK: oGBK.reverse_complement(), self))
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
            return "\n".join(list(map(lambda oGBK: oGBK._format_genbank(), self.get())))
        elif seq_format.upper()=="FASTA":
            return "\n".join(list(map(lambda oGBK: oGBK._format_fasta(), self.get())))
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
            oCopy.features += list(map(lambda ft: ft.copy(start=-index.start), 
                list(filter(lambda ft: ft.location.start >= index.start and ft.location.end <= index.stop, other.features))))
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
        oCopy.features += list(map(lambda ft: ft.copy(start=len(str(self.Seq))), other.features))
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
        oCopy.features += list(map(lambda ft: ft.copy(start=len(str(self.Seq)))+separator_length, other.features))
        return oCopy
    
    def rearrange(self,moltype="CDS",tag="",gene="",product="",strand=1):
        shift = 100
        # Create a copy
        oCopy = self.copy()
        # List of features
        features = list(map(lambda ft: ft.copy(), oCopy.features))
        if len(features) and moltype:
            features = list(filter(lambda ft: ft.type==moltype, features))
        if len(features) and tag:
            features = list(filter(lambda ft: 'locus_tag' in ft.qualifiers and ft.qualifiers['locus_tag'] and ft.qualifiers['locus_tag'][0]==tag, features))
        if len(features) and gene:
            features = list(filter(lambda ft: 'gene' in ft.qualifiers and ft.qualifiers['gene'] and ft.qualifiers['gene'][0]==gene, features))
        if len(features) and product:
            features = list(filter(lambda ft: 'product' in ft.qualifiers and ft.qualifiers['product'] and ft.qualifiers['product'][0]==product, features))
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
        return list(map(lambda ft: ft.copy(increment), self.features))
        
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
            features = list(filter(lambda oFT: oFT.type==moltype, self.features))
        else:
            features = self.features
        if moltype.upper() == "DNA":
            return "\n".join(list(map(lambda oFT: f">{oFT.get_feature_title()}\n{str(self.Seq[oFT.start-1:oFT.end])}")))
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
        gbk += (reduce(lambda ls1,ls2: ls1+ls2, map(lambda ft: ft.format(), self.features))+
            ["BASE COUNT   "+" ".join(list(map(lambda nuc: f"{nuc} {self.Seq.base_count(nuc)}", ['a','c','g','t']))),"ORIGIN"]+\
            self.Seq.format("genbank")+["//",""])
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
                
    def copy(self,reverse_complement=False):
        if reverse_complement:
            return GBK(self.title,self.heading,self.seq.copy(reverse_complement=True),list(map(lambda ft: ft.copy(complement=len(self.Seq)), self.features)))
        return GBK(self.title,self.heading,self.Seq.copy(),list(map(lambda ft: ft.copy(), self.features)))

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
        output = list(filter(lambda s: len(s.replace(" ","")) > 1, output))
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
        oCopy = Feature(self.type,self.strand,list(map(lambda exon: exon, self.location.get_exons(start=start,complement=complement))))
        oCopy.location = self.location.copy(start=start,complement=complement)
        oCopy.qualifiers = copy.deepcopy(self.qualifiers)
        return oCopy
        
##############################################################################################################
class GFF(Collection):
    def __init__(self,title=""):
        Collection.__init__(self,title)
        self.description = self.title
        self.body = self.Body = self.container
        self.heading = self.Heading = ""
        
    def copy(self):
        oGffCopy = GFF(self.title)
        oGffCopy.heading = self.heading
        for oRec in self:
            oGffCopy.append(oRec.copy())
        return oGffCopy
        
##############################################################################################################
class Record:
    def __init__(self,title="",**attributes):
        self.title = title
        for key, value in attributes.items():
            if isinstance(value, dict):
                setattr(self, key, Record("",**value))
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
class GffRecord(Record):
    def __init__(self,title="",**attributes):
        Record.__init__(self,title,**attributes)

##############################################################################################################
if __name__ == "__main__":
    oSeqIO = IO()
    path = os.path.join("..","input","S.aureus_150.gbk")
    oGBK = oSeqIO.read(path,"genbank")
    print(oGBK)
    #print(oGBK.features[1])
    #print(oGBK.features[1].qualifiers)
    #5/0
    oSeqIO.save(oGBK.format("genbank"),"S.aureus_150.gbk")
    oSeqIO.save(oGBK.format("fasta"),"S.aureus_150.fa")
    
