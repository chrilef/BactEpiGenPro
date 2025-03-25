import sys, os, string
path = os.path.join(os.getcwd(),"lib")
sys.path.append(path)
import main, tools, cui

version = "3.2.6"
date_of_creation = "16/03/2025"

def show_help():
    with open(os.path.join(path,"help")) as f:
        tools.msg(f.read())

###############################################################################
if __name__ == "__main__":
    options = cui.get_options()
    #### TEMP
    
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
   