import os, re
import seq_io, tools

# Get program run options with values by default
def get_options():
    oValidator = Validator()
    return oValidator.get_options()

# Get dictionary linking long and short argument names
def get_long_arguments():
    oValidator = Validator()
    return oValidator.get_long_arguments()

###############################################################################
# Menu
class Menu:
    def __init__(self, options, version="", date_of_creation=""):
        self.options = options
        self.oValidator = Validator()
        self.version = version
        self.date_of_creation = date_of_creation
        self.flg_show_all_menu = False
        self.completed = False
        # Input/Output object
        self.oIO = seq_io.IO()
        
    def show(self):
        while not self.completed:
            self.main_menu()
        return self.options

    def _set_options(self):
        if os.path.exists(os.path.join("lib","info")):
            info = self.oIO.read(os.path.join("lib","info"),"text",inlist=True,separator="\t",strip_symbol=" ")
            info = list(filter(lambda item: len(item)==2, info))
            for key,value in info:
                 if key in self.options:
                    self.options[key] = value
                    
    def main_menu(self):
        # SUBMENUES
        # General menu
        def general_menu(response):
            # Additional menu
            def general_additional_menu(response):
                # List of regions to filter
                # Sliding window length
                if response == "WL":
                    v = input(f"Enter sliding window length (positive integer from 5000 to 50000)? ")
                    if self.oValidator.validate_posnumber('-wl', v, lower_cutoff=5000, upper_cutoff=50000):
                        self.options['-wl'] = v
                        return True
                    else:
                        return False
                    
                # Sliding window step
                elif response == "WS":
                    v = input(f"Enter sliding window step (positive integer from 100 to 5000, but smaller or equal to sliding window length)? ")
                    if self.oValidator.validate_posnumber('-ws', v, lower_cutoff=100, upper_cutoff=int(self.options['-wl'])):
                        self.options['-ws'] = v
                        return True
                    else:
                        return False
                    
                elif response == "M":
                    self.options['-m'] = ""
                    generic_fname = input("Enter file name with filtered region locations? ")
                    if not generic_fname:
                        return False
                        
                    if self.oValidator.validate_field(generic_fname, "-m", os.path.join(self.options['-u'], self.options['-d'])):
                        self.options['-m'] = generic_fname
                        return True
                    else:
                        tools.alert(f"Check if the file {generic_fname} is in the folder '{os.path.join(self.options['-u'],self.options['-d'])}'")
                        return False

                # Promoter region length (additional option)
                elif response == "P":
                    v = input("Enter a positive number of expected promoter sequence length? ")
                    if self.oValidator.validate_posnumber("-p", v):
                        self.options['-p'] = int(v)
                        return True
                    return False
                
                # Set maximum sites for verification
                elif response == "R":
                    v = input("Enter a positive number of maximum number of modified sites for verification? ")
                    if self.oValidator.validate_posnumber("-r", v):
                        self.options['-r'] = int(v)
                        return True
                    else:
                        return False
                        
                # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
                elif response == "STD":
                    if self.options['-std'].upper() in ("LEADING", "+", "1", "DIR"):
                        self.options['-std'] = "lagging"
                    elif self.options['-std'].upper() in ("LAGGING", "-", "-1", "REV"):
                        self.options['-std'] = "off"
                    elif self.options['-std'].upper() in ("OFF", "0", "BOTH"):
                        self.options['-std'] = "leading"
                    else:
                        self.options['-std'] = "leading"
                    return True
            
                # Output graph format
                elif response == "OGF":
                    self.options['-ogf'] = self.select_output_graph_format(current_format=self.options['-ogf'])
                    return True
                    
                # Set DPI
                elif response == "DPI" and self.oValidator.is_raster_format(self.options['-ogf']):
                    v = input("Enter graph DPI value in range from 100 to 600: ")
                    try:
                        v = int(v)
                    except:
                        tools.alert("DPI must be an integer from 100 to 600!")
                        return False
                    if self.oValidator.validate_posnumber('-dpi', v,lower_cutoff=100, upper_cutoff=600):
                        self.options['-dpi'] = v
                        return True
                    return False
                
                # Context mismatches
                elif response == "N":
                    v = input("Enter allowed number of context mismatches? ")
                    if self.oValidator.validate_posnumber("-n", v, 0, 5):
                        self.options['-n'] = int(v)
                        return False
                    return True
                    
                # allow motif mismatch
                elif response == "Z":
                    self.options['-z'] = "No" if self.options['-z'] == "Yes" else "Yes"
                    return True
                    
                else:
                    return False
                
            # Run program execution
            if response == "Y":
                if self.oValidator.validate(self.options):
                    self.completed = True
                    return False
                else:
                    tools.alert("Check program run options!")
                    return True
                
            # Show additional options
            elif response == "~":
                if self.flg_show_all_menu:
                    self.flg_show_all_menu = False
                else:
                    self.flg_show_all_menu = True
                return True
            
            # MAIN MENU
            # Input GFF file
            elif response == "I":
                self.options['-i'] = ""
                generic_fname = input("Enter input GFF file name? ")
                if self.oValidator.validate_field(generic_fname, "-i", os.path.join(self.options['-u'],self.options['-d'])):
                    self.options['-i'] = generic_fname
                    return False
                else:
                    tools.alert(f"File {generic_fname} is not found in directory '{os.path.join(self.options['-u'],self.options['-d'])}'")
                    return True
                
            # Input GBK file
            elif response == "G":
                self.options['-g'] = ""
                generic_fname = input("Enter input GBK file name? ")
                if self.oValidator.validate_field(generic_fname, "-g", os.path.join(self.options['-u'], self.options['-d'])):
                    self.options['-g'] = generic_fname
                    return False
                else:
                    tools.alert(f"File {generic_fname} is not found in directory '{os.path.join(self.options['-u'],self.options['-d'])}'")
                    return True
                    
            # Subdirectory
            elif response == "D":
                subdirectory = input("Enter subdirectory name? ")
                if self.oValidator.validate_field(os.path.join(self.options['-u'], subdirectory), "-d"):
                    self.options['-d'] = subdirectory
                    return False
                else:
                    tools.alert(f"Subdirectory {subdirectory} is not found in directory '{os.path.join(self.options['-u'],self.options['-d'])}'")
                    return True
                
            # Generic file name
            elif response == "FT":
                fname = input("Enter generic output file name? ")
                self.options['-ft'] = self.oValidator.validate_file_name(fname)
                return True

            # Additional options of general menu
            if self.flg_show_all_menu:
                return general_additional_menu(response)
            else:
                return False
        
        # Circular map menu execution
        def circular_map_menu(response):
            # Additional options
            def circular_additional_menu(response):
                
                # CM Graph title
                if response == "CMT":
                    self.options['-cmt'] = input("Enter CM graph title? ")
                    return True
                return False

            
            # Set motif like GATC,1,-1
            if response == "W":
                self.options['-w'] = ""
                word = input("Enter motif like GATC,2,-2? ").upper()
                if not word:
                    return True
                elif self.oValidator.validate_word(word):
                    self.options['-w'] = word
                    return False
                else:
                    tools.alert(f"Not appropriate nucleotide motif '{word}'")
                    return True
            
            # Sites/Motifs 
            elif response == "S":
                self.options['-s'] = "sites" if self.options['-s'] == "motifs" else "motifs"
                return False
                
            # modified/Unmodified sites
            elif response == "F":
                self.options['-f'] = "U" if self.options['-f'] == "M" else "M"
                return False
            
            if self.flg_show_all_menu:
                return circular_additional_menu(response)
            else:
                return False
            
        # DOT-PLOT MENU
        def dotplot_menu(response):
            # dotplot additional menu
            def dotplot_additional_menu(response):
                # Score cut-off
                if response == "DPC":
                    v = input("Enter a positive number of score cut-off? ")
                    if self.oValidator.validate_posnumber("-dpc", v):
                        self.options['-dpc'] = int(v)
                        return True
                    else:
                        return False
                
                elif response == "DPW":
                    v = input("Enter a positive number of upper coverage limit? ")
                    if self.oValidator.validate_posnumber("-dpw", v):
                        self.options['-dpw'] = int(v)
                        return True
                    else:
                        return False
                    
                elif response == "DPS":
                    v = input("Enter a positive number of upper score limit? ")
                    if self.oValidator.validate_posnumber("-dps", v):
                        self.options['-dps'] = int(v)
                        return True
                    else:
                        return False
                        
                # DP Graph title
                elif response == "DPT":
                    self.options['-dpt'] = input("Enter DP graph title? ")
                    return True
                return False


            # Set dot-plot filter
            if response == "SFT":
                field, values = self.set_dotplot_filter()
                if field in ('-dpn', '-dpm', '-dpf'):
                    self.options[field] = values
                    if field == '-dpm' and values:
                        self.options['-dpn'] = ""
                    elif field == '-dpn' and values:
                        self.options['-dpm'] = ""
                    return False
                else:
                    return True
                
            # Additional options
            if self.flg_show_all_menu:
                return dotplot_additional_menu(response)
            else:
                return False
        
        # STATPLOT MENU
        def stat_panel_menu(response):
            def stat_panel_additional_menu(response):
                # Score cutoff
                if response == "C":
                    v = input("Enter a positive number of score cut-off? ")
                    if self.oValidator.validate_posnumber("-c", v):
                        self.options['-c'] = int(v)
                        return True
                    return False
                    
                # SP Graph title
                elif response == "SPT":
                    self.options['-sp'] = input("Enter SP graph title? ")
                    return True
                return False
                
            # Request filter setting
            if response == "SFT" and all([self.options[field].upper()[0] == "N" for field in ('-mm', '-dp')]):
                field, values = self.set_dotplot_filter()
                if field in ('-dpn', '-dpm', '-dpf'):
                    self.options[field] = values
                    if field == '-dpm' and values:
                        self.options['-dpn'] = ""
                    elif field == '-dpn' and values:
                        self.options['-dpm'] = ""
                    return True
                else:
                    return False

            # List of tasks
            elif response == "TSK":
                tasks = input(f"Enter comma-separated tasks: 'gc' for GC-content, 'gcs' for GC-skew and 'mge' for mobile genetic elements. ")
                if self.oValidator.validate_tasks(tasks):
                    self.options['-tsk'] = tasks
                    return True
                else:
                    return False
                
            # Additional options
            if self.flg_show_all_menu:
                return stat_panel_additional_menu(response)
            else:
                return False
        
        response = ''
        
        #### Main Cycle
        while response != "Q" and not self.completed:
            # Print menu
            self.print_main_menu()
            # Ask for user response
            print()
            response = input("?").upper()
            print()

            # GRAPH SELECTION
            # Show circular graph
            if response == "MM":
                self.options['-mm'] = "N" if self.options['-mm'][0].upper() == "Y" else "Y"
                continue
            # Show dot-plot
            elif response == "DP":
                self.options['-dp'] = "N" if self.options['-dp'][0].upper() == "Y" else "Y"
                continue
            # Show statistical panel
            elif response == "SP":
                self.options['-sp'] = "N" if self.options['-sp'][0].upper() == "Y" else "Y"
                continue
            
            # Execute general menu options
            if any([self.options[field].upper()[0] == "Y" for field in ('-mm', '-dp', '-sp')]):                
                if general_menu(response):
                    continue
            
            # Execute circular methylation map options
            if self.options['-mm'].upper()[0] == "Y":
                if circular_map_menu(response):
                    continue
            
            # Execute dotplot menu options
            if self.options['-dp'].upper()[0] == "Y":
                if dotplot_menu(response):
                    continue
            
            # Execute stat panel menu options
            # STAT PANEL MENU
            if self.options['-sp'].upper()[0] == "Y":
                if stat_panel_menu(response):
                    continue
                
            # Commands
            if response == "L":
                self._set_options()
            elif response == "H":
                self.show_help()
            elif response == "Q":
                tools.alert("The program is terminated")
                exit()
            else:
                continue
    
    def get_current_filter(self):
        filter_setting = ""
        if self.options['-dpn'].strip():
            filter_setting = f"Nucleotides: {self.options['-dpn']}; "
        if self.options['-dpm'].strip():
            filter_setting += f"Modification types: {self.options['-dpm']}; "
        if self.options['-dpf'].strip():
            filter_setting += f"Motifs: {self.options['-dpf']}; "
        return filter_setting
        
    def set_dotplot_filter(self):
        print()
        print("Filters\t\t\t\tCurrent setting:\n")
        print(f"\tDPN    Nucleotides filter\t: {self.options['-dpn']}")
        print(f"\tDPM    Methylation type filter\t: {self.options['-dpm']}")
        print(f"\tDPF    Motif filter\t\t: {self.options['-dpf']}")
        while True:
            response = input("Select filter or press Q to quit: ").upper()
            if response == "Q":
                return "", ""
            # Nucleotides to display
            elif response == "DPN":
                nucleotides = input("Enter modified nucleotide to display? ")
                if not self.oValidator.validate_field(nucleotides, "-dpn"):
                    tools.alert("Acceptable nucleotides are A,C,T,G")
                    continue
                return "-dpn", nucleotides
            # Modification types to display
            elif response == "DPM":
                mtypes = input("Enter modification types to display: m4C or m6A? ")
                if not self.oValidator.validate_field(mtypes, "-dpm"):
                    tools.alert("Acceptable modification types are m4C,m6A")
                    continue
                return "-dpm", mtypes
            # Filter motifs
            elif response == "DPF":
                motifs = input("Enter semicolon separated motifs to display or to filter out (-), like GATC,2,-2; -AGNCT,1,-1? ")
                if not self.oValidator.validate_field(motifs, "-dpf"):
                    tools.alert("Filter motif setting should be like GATC,2,-2; -AGNCT,1,-1")
                    continue
                return "-dpf", motifs
                
    def print_main_menu(self):
        # GRAPH SPECIFIC MENUS
        # Methylation map menu
        def print_MM_menu():
            # Additional option
            def print_additional_MM_menu():
                print(f"\tCMT  CM graph title\t\t: {self.options['-cmt']}")
                
            # Main options
            print("Methylation (circular) map settings:")
            print(f"\tW    Motif word\t\t\t: {self.options['-w']}")
            print(f"\tS    Search for\t\t\t: {self.options['-s']}")
            print(f"\tF    Modified/Unmodified\t: {self.options['-f']}")
    
            # Additional menu options
            if self.flg_show_all_menu:
                print_additional_MM_menu()
            
            # Separator
            print("=" * 41)
                
        def print_DP_menu():
            # Additional option
            def print_additional_DP_menu():
                
                print(f"\tDPC    Cut-off score\t\t: {self.options['-dpc']}")
                print(f"\tDPW    Maximum coverage\t\t: {self.options['-dpw']}")
                print(f"\tDPS    Maximum score\t\t: {self.options['-dps']}")
                print(f"\tDPT    DP graph title\t\t: {self.options['-dpt']}")
                
            # Main options
            print("Dot-plot graph settings:")
            print(f"\tSFT    Set \t\t\t: {self.get_current_filter()}")
            
            # Additional menu options
            if self.flg_show_all_menu:
                print_additional_DP_menu()
    
            # Separator
            print("=" * 41)
                
        def print_SP_menu():
            # Additional option            
            def print_additional_SP_menu():
                
                if all([self.options[field].upper()[0] == "N" for field in ('-mm', '-dp')]):
                    try:
                        print(f"\tC  Cut-off score\t\t: {int(self.options['-c'])}")
                    except:
                        print("\tC  Cut-off score\t\t: 0")
                        self.options['-c'] = 0
                
                print(f"\tSPT  SP graph title\t\t: {self.options['-spt']}")

            # Main options
            print("Statistical panel settings:")
            if all([self.options[field].upper()[0] == "N" for field in ('-mm', '-dp')]):
                # Set words or nucleotides for statistical evaluation
                print(f"\tSFT    Set filters\t\t: {self.get_current_filter()}")
            print(f"\tTSK    Genome properties\t: {self.options['-tsk']}")
        
            # Additional menu options
            if self.flg_show_all_menu:
                print_additional_SP_menu()

            # Separator
            print("=" * 41)
                
        #### MAIN MENU
        print(f"\nSeqWord Motif Mapper {self.version} {self.date_of_creation}")
        print()
        print("Settings for this run:\n")
        
        # General settings - printed if at least one graph is set for creation
        if any([self.options[field].upper()[0] == "Y" for field in ('-mm', '-dp', '-sp')]):
            print("General settings")
            print(f"  D    Subdirectory\t\t\t: {self.options['-d']}")
            print(f"  I    Input GFF file\t\t\t: {self.options['-i']}")
            print(f"  G    Genome GBK file\t\t\t: {self.options['-g']}")
            print(f"  FT   Generic file title\t\t: {self.options['-ft']}")
            # Additional menu options
            if self.flg_show_all_menu:
                print(f"  STD    Strand\t\t\t: {self.options['-std']}")
                print(f"  OGF  Output graph format\t\t: {self.options['-ogf']}")
                
                if self.oValidator.is_raster_format(self.options['-ogf'].strip().upper()):
                    print(f"  DPI  Graph dpi\t\t\t: {self.options['-dpi']}")
                
                try:
                    print(f"  WL   Sliding window length\t\t: {int(self.options['-wl'])}")
                    if int(self.options['-wl']) > 0:
                        try:
                            ws = int(self.options['-ws'])
                            if ws > int(self.options['-wl']):
                                self.options['-ws'] = self.options['-wl']
                        except:
                            self.options['-ws'] = self.options['-wl']
                        print(f"  WS   Sliding window step\t\t: {self.options['-ws']}")                        
                except:
                    print("  WL   Sliding window length\t: 0")
                    self.options['-wl'] = "0"

                try:
                    print(f"  P    Promoter sequence length\t\t: {int(self.options['-p'])}")
                except:
                    print("  P    Promoter sequence length\t\t: 0")
                    self.options['-p'] = 0
                    
                print(f"  M    Filter regions\t\t\t: {self.options['-m']}")
                
                print(f"  R    Max verification size\t\t: {self.options['-r']}")
                
                try:
                    print(f"  N    Context mismatches\t\t: {int(self.options['-n'])}")
                except:
                    print("  N    Context mismatches\t\t: 0")
                    self.options['-n'] = 0
                    
                print(f"  Z    Allow motif mismatch\t\t: {self.options['-z']}")

            # Separator
            print("=" * 41)                    
                       
        # Selecting types of output graphs
        print(f"  MM    Create methylation map\t\t: {self.options['-mm']}")
        if self.options['-mm'].upper()[0] == "Y":
            print_MM_menu()
                    
        # Dot-plot graph menu
        print(f"  DP    Create dot-plot graph\t\t: {self.options['-dp']}")
        if self.options["-dp"].upper()[0] == "Y":
            print_DP_menu()
        
        # Statistical panel
        print(f"  SP    Create statistical panel\t: {self.options['-sp']}")
        if self.options["-sp"].upper()[0] == "Y":
            print_SP_menu()
        
        # Services
        print("Services")
        print("  ~    show/hide additional menu options:")
        print("  L    set last used options\t\t; ")
        print("  Q    to quit\t\t\t\t;")
        print()
        if any([self.options[field].upper()[0] == "Y" for field in ('-mm', '-dp', '-sp')]):
            print("Y to accept these settings, type the letter of option to change setting, or Q to quit")
                
    # Select available output grap format
    def select_output_graph_format(self, current_format="HTML"):
        print("\nAvailable output graph formats:")
        output_graph_formats = list(self.oValidator.output_graph_formats.keys())
        for i in range(len(output_graph_formats)):
            graph_format = output_graph_formats[i]
            print(f"\t{i + 1}: {graph_format}")
        while True:
            graph_format = input(f"Enter graph format name, or its number, or press 'Q' to keep the current graph format {current_format}: ")
            if graph_format.strip().upper() == "Q":
                return current_format
            if graph_format.strip().upper() in output_graph_formats:
                return graph_format.strip().upper()
            try:
                i = int(graph_format) - 1
            except:
                tools.alert(f"Entered graph format {graph_format} is not recognized")
                continue
            if i < 0 or i >= len(output_graph_formats):
                tools.alert(f"Graph format number must be in range from 1 to {len(output_graph_formats)}!")
                continue
            return output_graph_formats[i]      
        
###############################################################################
# Validator
class Validator:
    def __init__(self):
        self.cwd = ""
        if __name__ == "__main__":
            self.cwd = ".."
        self.oIO = seq_io.IO()
        self.available_tasks = ("GC", "GCS", "MGE")
        self.available_strands = ("LEADING","LAGGING","OFF","-1","1","0","DIR","REV","BOTH")
        self.echo = True
        
        # Output graph formats
        self.output_graph_formats = {"HTML":False, "SVG":False, "EPS":False, "PDF":False, "JPG":True, "PNG":True, "TIF":True, "BMP":True}
        #self.output_graph_formats = {"HTML":False, "SVG":False, "PDF":False}

        # List of options with values by default            
        self.template_options = {
               # General settings
               "--general": {
                   "-d":"",             # subfolder
                   "-i":"",             # input GFF file name
                   "-g":"",             # genome GBK file
                   "-ft":"",            # Generic file title
                   "-ogf":"HTML",       # Output graph format
                   "-dpi":300,          # Raster dpi
                   "-m":"",             # filter regions
                   "-p":75,             # promoter sequence length
                   "-r":10000,          # maximum number of sites for verificatio, 0 - to skip checking
                   "-n":2,              # context mismatches
                   "-z":"Yes",          # allow motif mismatch
                   "-std":"off",        # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
                   "-u":"input",        # input folder
                   "-o":"output",       # output folder
                   "-tmp":"",           # tmp folder; if not set, default is ./lib/bin/tmp/
                   "-x":"",             # binpath; if not set, default is ./lib/bin/
                },
               # Circular map settings
               "--circular_map": {
                   "-mm":"N",           # Create dot-plot graph
                   "-w":"",             # motif
                   "-f":"M",            # modified/unmodified M | U
                   "-s":"sites",        # sites/motifs S | M
                   "-wl":8000,          # Sliding window length
                   "-ws":2000,          # Sliding window step
                   "-c":21,             # score cut-off
                   "-cmt":"",           # Circular map graph title
                },
               # Dot-plot settings
               "--dotplot": {
                   "-dp":"N",           # Create dot-plot graph
                   "-dpn":"A,C",        # comma separated nucleotides A,C,G,T to display
                   "-dpm":"",           # comma separated methylation types m6A,m4C to display
                   "-dpf":"",           # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
                   "-dpc":21,           # modification score cutoff
                   "-dpw":0,            # maximum coverage (X) value
                   "-dps":0,            # maximum score (Y) value
                   "-dpt":"",           # Dotplot graph title
               },
               # Statistics settings
               "--statplot": {
                   "-sp":"N",           # Create statistics panel
                   "-tsk":"gc, gcs",    # comma-separated tsks like gc, gcs, mge
                   "-spt":"",           # Statistical plot graph title
                },
        }
        
        # List of long names of arguments
        self.args = {
           "--input_GFF":"-i",              # input GFF file name
           "--input_GBK":"-g",              # genome GBK file
           "--project_directory":"-d",      # subfolder
           "--generic_file_name":"-ft",     # Generic graph title
           "--output_graph_format":"-ogf",  # Output graph format
           "--filter_file":"-m",            # filter regions
           "--promoter_length":"-p",        # promoter sequence length
           "--maximum_sites":"-r",          # maximum number of sites for verificatio, 0 - to skip checking
           "--blast_context_mismatch":"-n", # context mismatches
           "--blast_motif_mismatch":"-z",   # allow motif mismatch
           "--input_folder":"-u",           # input folder
           "--output_folder":"-o",          # output folder
           "--tmp_folder":"-tmp",           # output folder
           "--bin_folder":"-x",             # bin path 
           # Circular map settings
           "--circular_map":"-mm",          # Generate circular plot graph
           "--cmap_motif":"-w",            # motif like GATC,2,-2
           "--sites_or_motifs":"-s",        # Search for sites | motifs
           "--modified_or_unmodified":"-f", # modified/Unmetylated M | U
           "--window_length":"-wl",         # Sliding window length
           "--window_step":"-ws",           # Sliding window step
           "--cmap_score_cutoff":"-c",      # score cut-off
           "--cmap_graph_title":"-cmt",     # Circular map graph title
           # Dot-plot settings
           "--dotplot":"-dp",              # Generate dot-plot graph
           "--nucleotides":"-dpn",          # comma separated nucleotides A,C,G,T to display
           "--methylation_types":"-dpm",    # comma separated methylation types m6A,m4C to display
           "--dotplot_motifs":"-dpf",       # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
           "--dotplot_score_cutoff":"-dpc", # methylation score cutoff
           "--maximum_coverage":"-dpw",     # maximum coverage (X) value
           "--maximum_score":"-dps",        # maximum score (Y) value
           "--dp_graph_title":"-dpt",       # Dotplot graph title
           # Statistics settings
           "--statplot":"-sp",              # Create statistics panel
           "--tasks":"-tsk",                # comma-separated tsks like gc, gcs, mge
           "--strand":"-std",               # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
           "--sp_graph_title":"-spt",       # Statistical plot graph title
        }
    
    def get_long_arguments(self):
        return self.args

    def get_options(self):
        options = {}
        for key in self.template_options:
            options.update(self.template_options[key])
        return options
    
    def get_option_list(self, plots = {}):
        options = []
        for key in self.template_options:
            # Return options only for requested graphs
            if plots and key in self.args and plots[self.args[key]][0].upper() == "N":
                continue
            options += list(self.template_options[key].keys())
        return options
        
    def get_option_categories(self):
        return list(self.template_options.keys())
        
    def get_plot_options(self, options):
        return {"-mm":options["-mm"], 
            "-dp":options["-dp"], 
            "-sp":options["-sp"]}
            
    def is_raster_format(self, file_format):
        if file_format in self.output_graph_formats and self.output_graph_formats[file_format]:
            return True
        return False
    
    def validate(self, options, field="", echo=True):
        if not field:
            return self.validate_all(options, echo=echo)
        if field not in self.get_option_list():
            tools.alert(f"Unknown option {field}!")
            return False
        if field == "-d":
            return self.validate_field(options[field], field, para=options['-u'])
        if field == "-w":
            return self.validate_field(options[field], field, options=options)
        if field in ("-i", "-g", "-m", "-dpf"):
            return self.validate_field(options[field], field, para=os.path.join(options['-u'],options['-d']))
        return self.validate_field(options[field],field)
    
    def validate_field(self, v, field, para='', options={}):
        if field in ("-u", "-o", "-x", "-y", "-tmp"):
            return self.validate_path(v,field)
        if field == "-d":
            return self.validate_path(v,field,para)
        if field in ("-w", "-dpf"):
            return any([self.validate_word(v, options), 
                self.validate_input_files(para,[])])
        if field in ("-z", "-t", "-ft", "-cmt", "-dpt", "-spt"):
            return True
        if field == "-ogf":
            return v in self.output_graph_formats
        if field in ("-mm","-dp","-sp"):
            return self.validate_switches(v)
        if field == "-f":
            return v.upper()[0] in ("M","U")
        if field == "-i":
            return self.validate_input_files(os.path.join(para, v),["GFF","GF"])
        if field == "-g":
            return self.validate_input_files(os.path.join(para, v),["GBK","GB","GBF"])
        if field == "-m":
            return self.validate_input_files(os.path.join(para, v),[],True)
        if field in ("-r","-c","-p","-dpc","-dpw","-dps"):
            return self.validate_posnumber(field, v)
        if field == "-n":
            return self.validate_posnumber(field,v,0,5)
        if field == "-wl":
            return self.validate_posnumber(field,v,lower_cutoff=5000, upper_cutoff=50000)
        if field == "-ws":
            return self.validate_posnumber(field,v,lower_cutoff=100, upper_cutoff=5000)
        if field == "-dpi":
            return self.validate_posnumber(field,v,lower_cutoff=100, upper_cutoff=600)
        if field == "-s":
            if v not in ("sites","motifs"):
                if self.echo:
                    tools.alert("Option -s must be either 'sites' or 'motifs'")
                return False
            return True
        if field == "-dpn":
            return self.validate_nucleotides(v)
        if field == "-dpm":
            return self.validate_modification_types(v)
        if field == "-dpf":
            return self.validate_motifs(v)
        if field == "-tsk":
            return self.validate_tasks(v)
        if field == "-std":
            return self.validate_strand(v)
        return False
        
    def validate_all(self, options, echo=True):
        self.echo = echo
        plot_setting = self.get_plot_options(options = options)
        if all([field not in options or options[field].upper()[0] == "N" for field in ('-mm', '-dp', '-sp')]):
            tools.alert("Generation of at least one graph mast be requested!")
            return False
        for p in self.get_option_list(plots = plot_setting):
            if p not in ('-i','-m'):
                if options['-mm'].upper().startswith("N") and not p.upper().startswith("DP"):
                    continue
                if options['-dp'].upper().startswith("N") and p.upper().startswith("DP"):
                    continue
            valid = self.validate(options,p)
            if not valid:
                if self.echo:
                    tools.alert(f"Option {p} is not valid!")
                return False
        return True
    
    def validate_input_files(self, generic_fname, extensions, flg_allow_empty=False):
        if flg_allow_empty:
            return True
        if os.path.exists(generic_fname) and os.path.isfile(generic_fname) and (not extensions or generic_fname[generic_fname.rfind(".")+1:].upper() in extensions):
            return True
        if self.echo and extensions:
            tools.alert("File %s not found or its extension is not in %s" % (generic_fname,str(extensions)))
        return False
    
    def validate_file_name(self, fname, replacement="_"):
        for symbol in ("/","\\","\"","|","*","<",">",":","?"):
            fname = fname.replace(symbol,replacement)
        return fname
    
    def validate_path(self,path,field="",parent=""):
        path = os.path.join(self.cwd,parent,path)
        if os.path.exists(path):
            return True
        tools.alert("%s path %s does not exist!" % (field,path))
        return False
    
    def validate_word(self, word, options={}):
        if '-mm' in options and options['-mm'].upper()[0] == 'Y' and not word:
            if self.echo:
                tools.alert("Motif must not be empty!")
            return False
        return self.validate_motifs(word)

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
        lower_cutoff = int(lower_cutoff)
        upper_cutoff = int(upper_cutoff)
        if n < lower_cutoff or (upper_cutoff and n > upper_cutoff):
            if upper_cutoff:
                tools.alert("An interger number in range %d..%d is expected in field %s!" % (lower_cutoff,upper_cutoff,field))
            else:
                tools.alert("An interger number bigger than %d is expected in field %s!" % (lower_cutoff,field))
            return False
        return True
        
    def validate_nucleotides(self,nucleotides):
        if not nucleotides:
            return True
        return all([N.strip() in ['A','T','G','C'] for N in nucleotides.upper().split(",")])
        
    def validate_modification_types(self,mtypes):
        if not mtypes:
            return True
        return all([N.strip() in ['M4C','M6A'] for N in mtypes.upper().split(",")])
        
    def validate_tasks(self,tasks):
        if not tasks:
            return True
        if all([s.strip().upper() in self.available_tasks for s in tasks.split(',')]):
            return True
        else:
            tools.alert(f"Unknown task in the list {tasks.upper()}!\nAvailable tasks are: {self.available_tasks}")
            return False
            
    def validate_switches(self, switch):
        if switch and switch.upper()[0] in ("Y","N"):
            return True
        return False
        
    def validate_strand(self, strand):
        if strand and str(strand).upper() in self.available_strands:
            return True
        return False
        
    def validate_motifs(self,motifs):   # string of motifs is like "GATC,2,-2; -AGNCT,1,-1", where numbers identify modified nucleotides on both strands
        if not motifs:
            return True
        pattern = r"\b[AUGCTRIDHMNSVWYK]+\b"
        ls = [s.strip() for s in motifs.split(";")]
        for motif_set in ls:
            try:
                motif_set = [s.strip() for s in motif_set.split(",")]
                motif = motif_set[0][1:] if motif_set[0].startswith("-") else motif_set[0]
                positions = [abs(int(v)) for v in motif_set[1:]]
                if not re.fullmatch(pattern, motif, re.IGNORECASE):
                    if self.echo:
                        tools.alert(f"Wrong motif: {motif}. Motifs must include only nucleotide letter like AGNCT")
                    return False
                if any([v > len(motif) or v == 0 for v in positions]):
                    if self.echo:
                        tools.alert(f"Check modified nucleotide positions: {','.join([str(v) for v in positions])}.\n" + 
                        f"Positions must be positive or negative integers not bigger than the length of the motif: {len(motif)},\n" +
                        "and cannot be 0.")
                    return False
            except:
                if self.echo:
                    tools.alert(f"Wrong set of motifs: {motif_set}")
                return False
        return True
        
