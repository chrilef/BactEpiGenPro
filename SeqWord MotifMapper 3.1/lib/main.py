import sys, os, string
sys.path.append("lib")
import seq_io, tools, motifs, atlas

###############################################################################
# Command line interface
class Interface:
    def __init__(self,options=None):
        self.reference_genes = {"gyra":1.91}
        self.oValidator = Validator()
        self.oIO = seq_io.IO()
        self.cwd = "."
        self.filtered_loci = []
        if __name__ == "__main__":
            self.cwd = ".."

        self.options = {
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-x":"",             # binpath 
               "-y":"",             # tmppath 
               "-i":"",             # input GFF file name
               "-g":"",             # Genome GBK file
               "-m":"",             # filter regions
               "-w":"",             # motif
               "-d":0,              # in strand modified nucleotide
               "-r":0,              # reverse strand modified nucleotide
               "-c":0,              # score cut-off
               "-n":0,              # context mismatches
               "-z":"Yes",          # allow motif mismatch
               "-p":0,              # promoter sequence length
               "-f":"M",            # Save graphs
                }
        if options:
            self.options.update(options)
            
            if not self.options['-x']:
                self.options['-x'] = os.path.join(self.cwd,"lib","bin")
            if not self.options['-y']:
                self.options['-y'] = os.path.join(self.options['-x'],"tmp")
            
            valid = self.oValidator.validate(self.options)
            if valid:
                self.execute(self.options['-i'])
            else:
                self.main_menu()
        else:
            self.main_menu()
            
    def main_menu(self):
        response = ''
        while response != "Q":
            print("\nSeqWord MotifMapper 3.1 2024/02/19")
            print()
            print("Settings for this run:\n")
            print("  I    Input GFF file\t\t\t: " + str(self.options["-i"]))
            print("  G    Genome GBK file\t\t\t: " + str(self.options["-g"]))
            print("  M    Filter regions\t\t\t: " + str(self.options["-m"]))
            print("  W    Motif word\t\t\t: " + str(self.options["-w"]))
            print("  S    Search for\t\t\t: " + str(self.options["-s"]))
            try:
                print("  D    Modified base location\t\t: %s" % str(self.options["-d"]))
            except:
                print("  D    Modified base location\t\t: 0")
                self.options['-d'] = "0"
            try:
                print("  R    Reversed modification\t\t: %s" % str(self.options["-r"]))
            except:
                print("  R    Reversed modification\t\t: 0")
                self.options['-r'] = "0"
            print("  F    Methylated/Unmethylated sites\t: " + self.options['-f'])
            try:
                print("  C    Cut-off score\t\t\t: %d" % int(self.options["-c"]))
            except:
                print("  C    Cut-off score\t\t: 0")
                self.options['-c'] = 0
            try:
                print("  P    Promoter sequence length\t\t: %d" % int(self.options["-p"]))
            except:
                print("  P    Promoter sequence length\t\t: 0")
                self.options['-p'] = 0
            try:
                print("  N    Context mismatches\t\t: %d" % int(self.options["-n"]))
            except:
                print("  N    Context mismatches\t\t: 0")
                self.options['-n'] = 0
            
            print("  Z    Allow motif mismatch\t\t: " + self.options['-z'])
            print("  L    Set last used options\t\t; ")
            print("  Q    to quit;")
            print()
            print("Y to accept these settings, type the letter of option to change setting, or Q to quit")
            print()
            response = input("?").upper()
            print()
            if response == "Y":
                if self.oValidator.validate(self.options):
                    self.execute(self.options['-i'])
                else:
                    tools.alert("Check run options!")
                continue
            elif response == "I":
                self.options['-i'] = ""
                generic_fname = input("Enter input GFF file name? ")
                if self.oValidator.validate_field(os.path.join(self.options['-u'],generic_fname),"-i"):
                    self.options['-i'] = generic_fname
                else:
                    tools.alert("Check if the file %s is in the folder '%s'" % (generic_fname,self.options['-u']))
                continue
            elif response == "G":
                self.options['-g'] = ""
                generic_fname = input("Enter input GBK file name? ")
                if self.oValidator.validate_field(os.path.join(self.options['-u'],generic_fname),"-g"):
                    self.options['-g'] = generic_fname
                else:
                    tools.alert("Check if the file %s is in the folder '%s'" % (generic_fname,self.options['-u']))
                continue
            elif response == "M":
                self.options['-m'] = ""
                generic_fname = input("Enter file name with region locations? ")
                if not generic_fname:
                    continue
                if self.oValidator.validate_field(os.path.join(self.options['-u'],generic_fname),"-m"):
                    self.options['-m'] = generic_fname
                else:
                    tools.alert("Check if the file %s is in the folder '%s'" % (generic_fname,self.options['-u']))
                continue
            elif response == "W":
                self.options['-w'] = ""
                word = input("Enter motif? ").upper()
                if not word:
                    continue
                elif self.oValidator.validate_word(word):
                    self.options['-w'] = word
                else:
                    tools.alert("Not appropriate nucleotide motif '%s'" % word)
                continue
            elif response == "D":
                v = input("Enter modified base location from 1 to %d? " % len(self.options['-w']))
                if self.oValidator.validate_location(v,"-d",self.options['-w']):
                    self.options['-d'] = str(v)
                else:
                    tools.alert("Modified base location must be an integer from 1 to %d" % len(self.options['-w']))
                continue
            elif response == "R":
                v = input("Enter reversed modified base location from 0 to %d? " % len(self.options['-w']))
                if self.oValidator.validate_location(v,"-r",self.options['-w']):
                    self.options['-r'] = str(v)
                else:
                    tools.alert("Reversed modified base location must be an integer from 0 to %d" % len(self.options['-w']))
                continue
            elif response == "C":
                v = input("Enter a positive number of score cut-off? ")
                if self.oValidator.validate_posnumber("-c",v):
                    self.options['-c'] = int(v)
                continue
            elif response == "N":
                v = input("Enter allowed number of context mesmatches? ")
                if self.oValidator.validate_posnumber("-n",v,0,5):
                    self.options['-n'] = int(v)
                continue
            elif response == "P":
                v = input("Enter a positive number of expected promoter sequence length? ")
                if self.oValidator.validate_posnumber("-p",v):
                    self.options['-p'] = int(v)
                continue
            elif response == "S":
                if self.options['-s'] == "motifs":
                    self.options['-s'] = "sites"
                else:
                    self.options['-s'] = "motifs"
            elif response == "F":
                if self.options['-f'] == "M":
                    self.options['-f'] = "U"
                else:
                    self.options['-f'] = "M"
            elif response == "Z":
                if self.options['-z'] == "Yes":
                    self.options['-z'] = "No"
                else:
                    self.options['-z'] = "Yes"
            elif response == "S":
                self._set_options()
                continue
            elif response == "L":
                self._set_options()
                continue
            elif response == "H":
                self.show_help()
                continue
            elif response == "Q":
                print()
                print("The program was terminated")
                exit()
            else:
                continue
            
    def _get_outfilename(self,gff_filename,oMotifObj):
        def suffix(modbase_location,rev_modbase_location):
            return "_%s,%s_%d" % (modbase_location,rev_modbase_location,int(self.options['-c']))
        def prefix():
            prx = ""
            if self.options['-s'] == "motifs":
                prx = "w"
            if self.options['-f'] == "U":
                prx = "u"+prx
            return prx
        return os.path.join(self.options['-o'],"%s%s_%s%s.txt" % (prefix(),gff_filename[:-4],oMotifObj.get_word(),
            suffix(str(oMotifObj.get_modbase_location()),str(oMotifObj.get_rev_modbase_location()))))
                
    def _parse_loci(self,fname):
        if not os.path.exists(fname):
            return []
        try:
            loci = self.oIO.read(fname,"text",inlist=tTrue)
        except:
            return []
        if not loci:
            return []
        loci = list(map(lambda item: item.strip(), loci))
        loci = list(filter(lambda item: item, loci))
        if not loci:
            return []
        success = False
        for symbol in ("-",".."," "):
            if loci[0].find(symbol) > 0:
                success = True
                break
        if not success:
            return []
        loci = list(map(lambda item: item.split(symbol), loci))
        loci = list(filter(lambda item: len(item)==2, loci))
        if not loci:
            return []
        try:
            loci = list(map(lambda item: [int(item[0]),int(item[1])], loci))
            list(map(lambda item: item.sort(), loci))
        except:
            return []
        return loci
    
    def _get_oMotifObj(self):
        find_methylated_sites = False
        if self.options['-f'] == "M":
            find_methylated_sites = True
        oMotifObj = motifs.Motif(reference=self.options['-g'],binpath=self.options['-x'],tmppath=self.options['-y'],
            word=self.options['-w'],search_mode=self.options['-s'],modbase_location=self.options['-d'],reverse_modbase_location=self.options['-r'],
            score_cutoff=self.options['-c'],promoter_length=self.options['-p'],context_mismatches=self.options['-n'],
            motif_mismatch=self.options['-z'],filtered_loci=self.filtered_loci,flg_show_methylated_sites=find_methylated_sites)
        return oMotifObj
    
    def _parse_motif_list(self):
        MotifObjs = []
        settings = list(filter(lambda ls: len(ls[0]), self.oIO.read(os.path.join(self.options['-u'],self.options['-w']),"text",inlist=True,separator=",",strip_symbol=" ")))
        for dataset in settings:
            try:
                self.options['-w'] = dataset[0]
                if len(dataset) > 1:
                    self.options['-d'] = dataset[1]
                else:
                    self.options['-d'] = "0"
                if len(dataset) > 2:
                    self.options['-r'] = dataset[2]
                else:
                    self.options['-r'] = "0"
                if len(dataset) > 3 and dataset[3]=="u":
                    self.options['-f'] = "U"
                else:
                    self.options['-f'] = "M"
                MotifObjs.append(self._get_oMotifObj())
            except:
                tools.alert("Motif object has not been created for %s!" % dataset)
                MotifObjs.append(None)
        return MotifObjs
            
    def _set_options(self):
        if os.path.exists(os.path.join("lib","info")):
            info = self.oIO.read(os.path.join("lib","info"),"text",inlist=True,separator="\t",strip_symbol=" ")
            info = list(filter(lambda item: len(item)==2, info))
            for key,value in info:
                 if key in self.options:
                    self.options[key] = value
                    
    def _save_options(self):
        self.oIO.save("\n".join(map(lambda item: "%s\t%s" % (item[0],str(item[1])), self.options.items())),os.path.join("lib","info"),"text")
        
    def show_help(self):
        pass

    # Execute selected program
    def execute(self,generic_fname,first_name=""):
        gbk_path = os.path.join(self.options['-u'],self.options['-g'])
        if not self.options['-g'] or not os.path.exists(gbk_path):
            gbk_path = None
        gff_files = list(filter(lambda f: f[-4:].upper()==".GFF", os.listdir(self.options['-u'])))
        if self.options['-i']:
            gff_files = list(filter(lambda fname: fname.find(self.options['-i']) > -1, gff_files))
        if self.options["-m"] and os.path.exists(os.path.join(self.options['-u'],self.options['-m'])):
            try:
                self.filtered_loci = list(map(lambda ls: [int(ls[0]),int(ls[1])], 
                    self.oIO.read(os.path.join(self.options['-u'],self.options['-m']),"text",inlist=True,separator="-",strip_symbol=" ")))
            except:
                tools.msg("Error when parsing %s file" % self.options['-m'])
                return
        if os.path.exists(os.path.join(self.options['-u'],self.options['-w'])):
            oMotifObjLs = self._parse_motif_list()
        else:
            oMotifObjLs = [self._get_oMotifObj()]
        for oMotifObj in oMotifObjLs:
            if oMotifObj == None:
                continue
            for gff_file in gff_files:
                tools.msg("Searching for motifs %s in %s" % (oMotifObj.get_description(),gff_file))
                oMotifObj.reset()
                success = oMotifObj.execute(os.path.join(self.options['-u'],gff_file),gbk_path)
                if success:
                    outfile = self._get_outfilename(gff_file,oMotifObj)
                    self.oIO.save(oMotifObj.tostring(),outfile,"text")
                    if self.options['-f'] in ("M","U"):
                        if self.options['-s'] == 'motifs':
                            msg_found = "%d/%d" % (oMotifObj.num_found_motifs,oMotifObj.num_expected_sites)
                        elif self.options['-f'] == "M":
                            msg_found = "%d/%d" % (oMotifObj.num_found_sites,oMotifObj.num_expected_sites)
                        elif self.options['-f'] == "U":
                            msg_found = "%d/%d" % (oMotifObj.num_found_sites,oMotifObj.num_expected_sites)
                        else:
                            msg_found = ""
                        if self.options['-g']:
                            oAtlas = atlas.Main(os.path.join(self.options['-u'],self.options['-g']),oMotifObj.get_word(),oMotifObj.get_modbase_dict(),self.filtered_loci)
                            svg = oAtlas.svg([oMotifObj.get_modbase_location(),oMotifObj.get_rev_modbase_location(),msg_found])
                            self.oIO.save(svg,outfile[:-4]+".svg","text")
                    self._save_options()
                else:
                    tools.alert("Error during program execution with the input file %s searching for motif %s!" % 
                        (gff_file,"%s_%s,%s" % (oMotifObj.get_word(),oMotifObj.get_modbase_location(),oMotifObj.get_rev_modbase_location())))
        
###############################################################################
# Validator
class Validator:
    def __init__(self):
        self.cwd = ""
        self.oIO = seq_io.IO()
        self.fields = ("-o","-u","-x","-y","-z","-i","-g","-m","-n","-w","-d","-r","-c","-p","-s")
        if __name__ == "__main__":
            self.cwd = ".."
        
    def validate(self,options,field=""):
        if not field:
            return self.validate_all(options)
        if field in ("-u","-o","-x","-y"):
            return self.validate_path(options[field],field)
        if field == "-w":
            return self.validate_word(options[field],options['-u'])
        if field in ("-i","-z"):
            return True
        if field == "-f":
            return options[field] in ("M","U")
        if field == "-g":
            return self.validate_input_files(os.path.join(options['-u'],options[field]),["GBK","GB"])
        if field == "-m":
            return self.validate_input_files(os.path.join(options['-u'],options[field]),[])
        if field in ('-d','-r'):
            return self.validate_field(options[field],field,options['-w'])
        if field in ("-c","-p","-n"):
            return self.validate_field(options[field],field)
        if field == "-s":
            return self.validate_field(options[field],field)
        return False
    
    def validate_field(self,v,field,para=0):
        if field == "-u":
            return self.validate_path(v,"-u")
        if field == "-o":
            return self.validate_path(v,"-o")
        if field == "-x":
            return self.validate_path(v,"-x")
        if field == "-y":
            return self.validate_path(v,"-y")
        if field == "-w":
            return self.validate_word(v)
        if field in ("-i","-z"):
            return True
        if field == "-f":
            return v in ("M","U")
        if field == "-g":
            return self.validate_input_files(v,["GBK","GB"])
        if field == "-m":
            return self.validate_input_files(v,[])
        if field in ('-d','-r'):
            return self.validate_location(v,field,para)
        if field in ("-c","-p"):
            return self.validate_posnumber(field,v)
        if field == "-n":
            return self.validate_posnumber(field,v,0,5)
        if field == "-s":
            if v not in ("sites","motifs"):
                tools.alert("Option -s must be either 'sites' or 'motifs'")
                return False
            return True
        return False
        
    def validate_all(self,options):
        for p in self.fields:
            valid = self.validate(options,p)
            if not valid:
                #tools.alert(f"Option {p} is not valid!")
                return False
        return True
    
    def validate_input_files(self,generic_fname,extensions,flg_allow_empty=False):
        if flg_allow_empty and not generic_fname:
            return True
        if os.path.exists(generic_fname) and (not extensions or generic_fname[generic_fname.rfind(".")+1:].upper() in extensions):
            return True
        tools.alert("File %s not found or its extension is not in %s" % (generic_fname,str(extensions)))
        return False
    
    def validate_path(self,path,field=""):
        if os.path.exists(os.path.join(self.cwd,path)):
            return True
        tools.alert("%s path %s does not exist!" % (field,path))
        return False
    
    def validate_word(self,word,input_folder="input"):
        if tools.is_dna_sequence(word):
            return True
        if os.path.exists(os.path.join(input_folder,word)):
            return True
        tools.alert("Motif must be either a DNA sequence or a path to the list of motifs in folder %s!" % input_folder)
        return False
        
    def validate_location(self,p,field,word):
        p = str(p)
        try:
            values = list(map(lambda v: int(v), p.split(",")))
        except:
            return False
        for v in values:
            if v < 0 or v > len(word):
                return False
        return True
    
    def validate_posnumber(self,field,n,lower_cutoff=0,upper_cutoff=0):
        try:
            n = int(n)
        except:
            tools.alert("An interger number is expected in field %s!" % field)
            return False
        if n < lower_cutoff or (upper_cutoff and n > upper_cutoff):
            if upper_cutoff:
                tools.alert("An interger number in range %d..%d is expected in field %s!" % (lower_cutoff,upper_cutoff,field))
            else:
                tools.alert("An interger number bigger than %d is expected in field %s!" % (lower_cutoff,field))
            return False
        return True
