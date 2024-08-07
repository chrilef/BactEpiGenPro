import sys, os, string
path = os.path.join(os.getcwd(),"lib")
sys.path.append(path)
import main

version = "3.0"

def show_help():
    f = open(os.path.join(path,"help"))
    print(f.read())
    f.close()

###############################################################################
if __name__ == "__main__":
    options = {
               "-u":"input",        # input folder
               "-o":"output",       # output folder
               "-x":"",             # binpath 
               "-y":"",             # tmppath 
               "-i":"",             # input GFF file name
               "-g":"",             # Genome GBK file
               "-m":"",             # filter regions
               "-w":"",             # motif
               "-s":"sites",        # Search for sites | motifs
               "-d":0,              # in strand modified nucleotide
               "-r":0,              # reverse strand modified nucleotide
               "-c":0,              # score cut-off
               "-n":2,              # context mismatches
               "-z":"Yes",          # allow motif mismatch
               "-p":75,             # promoter sequence length
               "-f":"M",            # Save graphs
             }

    arguments = sys.argv[1:]
    if arguments:
        for i in range(0,len(arguments),2):
            key = arguments[i].replace("\"","").replace("'","")
            if key in ("-h", "-H", "--help"):
                show_help()
                exit()
            if key in ("-v", "-V", "--version"):
                print(version)
                exit()
            if key not in options:
                raise IOError("Unknown argument " + key + "!")
            if i <= len(arguments)-2:
                options[key] = arguments[i+1]
    
    oMain = main.Interface(options)
