import sys, os, string
path = os.path.join(os.getcwd(),"lib")
sys.path.append(path)
import main, tools, cui

version = "3.2.6"
date_of_creation = "03/03/2025"

def show_help():
    with open(os.path.join(path,"help")) as f:
        tools.msg(f.read())

###############################################################################
if __name__ == "__main__":
    options = cui.get_options()
    #### TEMP
    '''
    options = {
               # General settings
               "-i":"Nanohalo.gff",             # input GFF file name
               "-g":"Nanohalo.gbk",             # genome GBK file
               "-d":"DPANN",             # subfolder
               "-t":"",             # Generic graph title
               "-m":"",             # filter regions
			   "-gt":"",            # Graph title
			   "-ft":"",            # Generic file title
			   "-ogf":"SVG",       # Output graph format: SVG, HTML, PDF, EPS, JPG, JPEG, TIF, TIFF, PNG, BMP
               "-c":0,              # score cut-off
               "-n":2,              # context mismatches
               "-z":"Yes",          # allow motif mismatch
               "-p":120,             # promoter sequence length
               "-r":0,           # maximum number of sites for verificatio, 0 - to skip checking
               # Circular map settings
               "-mm":"Y",         # Generate circular plot graph
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-x":"",             # bin path 
               "-y":"",             # tmp path 
               "-w":"GDGCHC,4,-4",             # motif like GATC,2,-2
               "-s":"motifs",        # Search for sites | motifs
               "-f":"M",            # modified/unmodified M | U
               "-wl":8000,          # Sliding window length
               "-ws":2000,          # Sliding window step
               # Dot-plot part
               "-dp":"N",         # Generate dot-plot graph
               "-dpn":"C",          # comma separated nucleotides A,C,G,T to display
               "-dpm":"",           # comma separated methylation types m6A,m4C to display
               "-dpf":"-GDGCHC,4,-4",           # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
               "-dpc":21,           # methylation score cutoff
               "-dpw":0,            # maximum coverage (X) value
               "-dps":0,            # maximum score (Y) value
               # Statistics settings
               "-sp":"N",           # Create statistics panel
               "-tsk":"gc, gcs",        # comma-separated tsks like gc, gcs, mge
               "-std":"off",           # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
             }
    
    options = {
               # General settings
               "-d":"A.borkumensis",             # subfolder
               "-i":"FTB2.gff",             # input GFF file name
               "-g":"A.borkumensis_SK2.gbk",             # genome GBK file
               "-t":"",             # Generic graph title
               "-m":"",             # filter regions
			   "-ft":"my_output",            # Generic file title
			   "-ogf":"SVG",       # Output graph format: HTML, SVG, PDF
               "-n":2,              # context mismatches
               "-z":"Yes",          # allow motif mismatch
               "-p":120,             # promoter sequence length
               "-r":0,           # maximum number of sites for verificatio, 0 - to skip checking
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-tmp":"",       # tmp folder
               # Circular map settings
               "-mm":"N",         # Generate circular plot graph
               "-x":"",             # bin path 
               "-c":0,              # score cut-off
               #"-w":"AGGCCT,5,-5",             # motif like GATC,2,-2
               "-w":"GaTNNNNNGtGG,2,-3",             # motif like GATC,2,-2
               "-s":"sites",        # Search for sites | motifs
               "-f":"M",            # modified/unmodified M | U
               "-wl":8000,          # Sliding window length
               "-ws":2000,          # Sliding window step
               "-cmt":"A. borcumensis GaTNNNNNGtGG,2,-3",           # Circular map graph title
               # Dot-plot part
               "-dp":"Y",         # Generate dot-plot graph
               "-dpn":"C",          # comma separated nucleotides A,C,G,T to display
               "-dpm":"",           # comma separated methylation types m6A,m4C to display
               #"-dpf":"-GaTNNNNNGtGG,2,-3",           # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
               "-dpf":"-AgGCcT,5,-5",           # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
               "-dpc":21,           # modification score cutoff
               "-dpw":0,            # maximum coverage (X) value
               "-dps":0,            # maximum score (Y) value
               "-dpt":"C -AgGCcT,5,-5",           # Dotplot graph title
               # Statistics settings
               "-sp":"Y",           # Create statistics panel
               "-tsk":"gc, gcs, mge",        # comma-separated tsks like gc, gcs, mge
               "-std":"lagging",           # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
               "-spt":"",           # Statplot graph title
             }
    '''
    options = {
               # General settings
               "-d":"E.coli",             # subfolder
               "-i":"NC2.E.coli.gff",             # input GFF file name
               "-g":"E.coliBAA196.gbk",             # genome GBK file
               "-t":"",             # Generic graph title
               "-m":"",             # filter regions
			   "-ft":"NC2.E.coli",            # Generic file title
			   "-ogf":"SVG",       # Output graph format: HTML, SVG, PDF
               "-n":2,              # context mismatches
               "-z":"No",          # allow motif mismatch
               "-p":120,             # promoter sequence length
               "-r":0,           # maximum number of sites for verificatio, 0 - to skip checking
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-tmp":"",       # tmp folder
               # Circular map settings
               "-mm":"N",         # Generate circular plot graph
               "-x":"",             # bin path 
               "-c":0,              # score cut-off
               "-w":"GATC,2,-2",             # motif like GATC,2,-2
               #"-w":"cRGKGATC,1,6,-2",             # motif like GATC,2,-2
               #"-w":"CRGKGATCMCYG,1,6,-1,-6",             # motif like GATC,2,-2
               "-s":"sites",        # Search for sites | motifs
               "-f":"U",            # modified/unmodified M | U
               "-wl":8000,          # Sliding window length
               "-ws":2000,          # Sliding window step
               "-cmt":"",           # Circular map graph title
               # Dot-plot part
               "-dp":"Y",         # Generate dot-plot graph
               "-dpn":"A,C",          # comma separated nucleotides A,C,G,T to display
               "-dpm":"",           # comma separated methylation types m6A,m4C to display
               "-dpf":"CRGKGATC,1,6,-2",           # simicolone separated methylation motifs to include or exclude, like GATC,2,-2; -AGNCT,1,-1
               "-dpc":21,           # modification score cutoff
               "-dpw":0,            # maximum coverage (X) value
               "-dps":0,            # maximum score (Y) value
               "-dpt":"",           # Dotplot graph title
               # Statistics settings
               "-sp":"Y",           # Create statistics panel
               "-tsk":"gc, gcs, mge",        # comma-separated tsks like gc, gcs, mge
               "-std":"lagging",           # exclude strand (leading | lagging | off, off by default; or -1, 1, 0, or dir, rev, both)
               "-spt":"",           # Statplot graph title
             }
    
    args = sys.argv[1:]
    if args:
        long_arguments = cui.get_long_arguments()
        for i in range(0,len(args),2):
            key = args[i].replace("\"","").replace("'","")
            # Show help
            if key in ("-h", "-H", "--help"):
                show_help()
                exit()
            # Show version
            if key in ("-v", "-V", "--version"):
                tools.msg(f"Version {version} created on {date_of_creation}")
                exit()
            # Translate long argument to option key    
            if key in long_arguments:
                key = long_arguments[key]
            # Process wrong argument
            if key not in options:
                raise IOError("Unknown argument " + key + "!")
            # Set option value to the respective argument
            if i <= len(args)-2:
                options[key] = args[i+1]
    
    oMain = main.Interface(options, version, date_of_creation)
