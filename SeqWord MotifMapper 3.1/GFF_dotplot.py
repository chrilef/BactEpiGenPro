import sys, os, string, math

def svg(data,title,MaxX=0,MaxY=0,width=800,height=500):
    def round(v):
        return 10*math.ceil(v/10)
    def modification(l,mod):
        if mod=="modified_base":
            mod += "-%s" % l
        return mod

    opacity = 0.8
    font_size = 12
    left_border = 80
    top_border = 10
    bottom_border = 80
    X0 = left_border
    Y0 = height-bottom_border
    svg = ["<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height)]
    # Y axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />" % 
            (X0,top_border,X0,height-bottom_border+10))
    # X axis
    svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />" % 
            (10,Y0,width,Y0))
    #### Data processing
    dot_colors = {"m6A":"red","m4C":"darkgreen",
        "modified_base-G":"darkred","modified_base-C":"forestgreen","modified_base-A":"deeppink","modified_base-T":"blue"}
    x_values = list(map(lambda i: float(data[i]['data']['coverage']), range(len(data))))
    y_values = list(map(lambda i: float(data[i]['score']), range(len(data))))
    nucmod = list(map(lambda i: modification(data[i]['data']['nucleotide'],data[i]['modification']), range(len(data))))
    if not MaxX:
        MaxX = round(max(x_values))
    if not MaxY:
        MaxY = round(max(y_values))
    x_step = float(width-X0)/MaxX
    y_step = float(height-top_border-bottom_border)/MaxY
    for i in range(10):
        # Y axis
        yV = int(MaxY - i*MaxY/10)
        y = Y0-y_step*yV
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />" % 
                (X0-10,y,X0,y))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>" %
            (X0-10,y,font_size,"end",yV))
        # X axis
        xV = int(MaxX - i*MaxX/10)
        x = X0+x_step*xV
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />" % 
                (x,Y0+10,x,Y0))
        svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%d</text>" %
            (x,Y0+10+font_size,font_size,"middle",xV))
    
    for i in range(len(x_values)):
        x = X0+x_values[i]*x_step
        y = Y0-y_values[i]*y_step
        if nucmod[i] in ("m6A","m4C"):
            svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" stroke=\"black\" stroke-width=\"1\" fill=\"%s\" fill-opacity=\"%f\"/>" % (x,y,dot_colors[nucmod[i]],opacity))
        else:
            svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" fill=\"%s\" fill-opacity=\"%f\"/>" % (x,y,dot_colors[nucmod[i]],opacity))
    # Titles
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+60,top_border+font_size+10,font_size+10,"start",title))    
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+(width-X0)/2,height-bottom_border/2,font_size+10,"middle","Coverage"))    
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\" transform=\"rotate(-90 %d %d)\">%s</text>" %
        (X0/2,(height-top_border-bottom_border)/2,font_size+10,"middle",X0/2,(height-top_border-bottom_border)/2,"NucMod Score"))    
    # Legend
    y = height-font_size
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" stroke=\"black\" stroke-width=\"1\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+50,y,dot_colors["m6A"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+60,y+3,font_size+5,"start","- m6A;"))
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" stroke=\"black\" stroke-width=\"1\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+150,y,dot_colors["m4C"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+160,y+3,font_size+5,"start","- m4C;"))
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+250,y,dot_colors["modified_base-A"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+260,y+3,font_size+5,"start","- modA;"))
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+350,y,dot_colors["modified_base-T"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+360,y+3,font_size+5,"start","- modT;"))
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+450,y,dot_colors["modified_base-G"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+460,y+3,font_size+5,"start","- modG;"))
    svg.append("<circle cx=\"%f\" cy=\"%f\" r=\"5\" fill=\"%s\" fill-opacity=\"%f\"/>" % (X0+550,y,dot_colors["modified_base-C"],opacity))
    svg.append("<text x=\"%d\" y=\"%d\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"%d\" style=\"text-anchor:%s\">%s</text>" %
        (X0+560,y+3,font_size+5,"start","- modC;"))
    
    svg.append("</svg>")
    return "\n".join(svg)

def execute(input_file, output_file, title="", cutoff=0, nucleotides=[], mtypes=[], width=0, height=0):
    dots = []
    # Open input file
    with open(input_file) as f:
        # Parse GFF file
        lines = f.readlines()
        dots = []
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            # Create GFF records as dictionaries
            dots.append(dict(zip(["modification","start","end","score","strand","dot","data"],line.split("\t")[2:])))
            dots[-1]['data'] = dict(zip(list(map(lambda s: s.split("=")[0].strip(), dots[-1]['data'].split(";"))),
                list(map(lambda s: s.split("=")[1].strip(), dots[-1]['data'].split(";")))))
            dots[-1]['data']['nucleotide'] = dots[-1]['data']['context'][len(dots[-1]['data']['context'])//2]
            
        # Filter GFF records by modified nucleotides
        if nucleotides[0]:
            dots = list(filter(lambda dot: dot['data']['nucleotide'] in nucleotides, dots))
        
        # Filter GFF records by types of methylation
        if mtypes[0]:
            dots = list(filter(lambda dot: dot['modification'] in mtypes, dots))
        
        # Filter GFF records by score cutoff values
        if cutoff:
            dots = list(filter(lambda dot: int(dot['score']) >= cutoff, dots))
        
        # Set maximal X (coverage) and Y (score) values
        try:
            if not width:
                width = max(list(map(lambda i: int(dots[i]['data']['coverage']), range(len(dots)))))
                width = int(str(width)[0]+("0"*(len(str(width))-1)))+10**(len(str(width))-1)
            if not height:
                height = max(list(map(lambda i: int(dots[i]['score']), range(len(dots)))))
                height = int(str(height)[0]+("0"*(len(str(height))-1)))+10**(len(str(height))-1)
        except:
            raise IOError(f"File {input_file} is empty or corrupted!")
    
    # Set graph title
    if not title:
        title = os.path.basename(input_file[:-4])
        
    with open(output_file,"w") as outfile:
        outfile.write(svg(dots,title,width,height))
        
def show_help():
    help_string = '''
    
    Program GFF DotPlot 24.07.2024
    
    usage: python GFF_dotplot.py [-h] [-v] /path/my_file.gff [options]
    
    options:
        -o: output SVG files, empty by default to save the output into 
            the input folder under the input file name.
        -n: comma separated nucleotides A,C,G,T; empty by default.
        -m: comma separated methylation types m4C,m6A; empty by default.
        -t: graph title; empty by default.
        -c: score cut-off; 0 by default.
        -w: WIDTH - maximum coverage (X) value; 0 (AUTO) by default.
        -s: HEIGHT - maximum score (Y) value; 0 (AUTO) by default.
    '''
    print(help_string)
       
###############################################################################
if __name__ == "__main__":
    __version__ = "GFF_dotplot 24.07.2024"
    args = {
               "-o":"",             # output folder
               "-n":"",             # comma separated nucleotides A,C,G,T
               "-m":"",             # comma separated methylation types m4C,m6A
               "-t":"",             # Graph title
               "-c":0,              # score cut-off
               "-w":0,              # WIDTH - maximal coverage value
               "-s":0,              # HEIGHT - maximal score value
            }

    # Parsing arguments provided in the command line
    arguments = sys.argv[1:]
    if len(arguments) < 1:
        show_help()
        exit()

    # Check input GFF file
    if arguments[0] in ("-h", "-H", "--help"):
        show_help()
        exit()
    if arguments[0] in ("-v", "-V", "--version"):
        print(__version__)
        exit()
    INPUT_FILE = arguments[0]
    if not INPUT_FILE or not os.path.exists(INPUT_FILE):
        raise IOError(f"File {INPUT_FILE} does not exist!")
    if not INPUT_FILE.lower().endswith('.gff'):
        raise IOError(f"Input file {INPUT_FILE} is not a GFF file!")
    
    for i in range(1,len(arguments),2):
        key = arguments[i].replace("\"","").replace("'","")
        if key not in args:
            raise IOError("Unknown argument " + key + "!")
        if i <= len(arguments)-2:
            args[key] = arguments[i+1]
    
    # Check output SVG file
    if not args['-o']:
        OUTPUT_FILE = INPUT_FILE[:-4] + ".svg"
    else:
        OUTPUT_FILE = args['-o']
    if not OUTPUT_FILE.lower().endswith(".svg"):
        OUTPUT_FILE += ".svg"
    
    # Check other argument settings
    try:
        HEIGHT = int(args['-s'])
    except:
        raise IOError("'-s' must be 0 or a positive integer!")
    try:
        WIDTH = int(args['-w'])
    except:
        raise IOError("'-w' must be 0 or a positive integer!")
    try:
        SCORE_CUTOFF = int(args['-c'])
    except:
        raise IOError("'-c' must be 0 or a positive integer!")
    NUCLEOTIDES = list(map(lambda s: s.strip().upper(), args['-n'].split(',')))
    MTYPES = list(map(lambda s: s.strip(), args['-m'].split(',')))

    execute(input_file=INPUT_FILE, output_file=OUTPUT_FILE, title=args['-t'],
        cutoff=SCORE_CUTOFF, nucleotides=NUCLEOTIDES, mtypes=MTYPES, width=WIDTH, height=HEIGHT)
