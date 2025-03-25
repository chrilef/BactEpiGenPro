import sys, os, string, math, re
import GI_finder, seq_io, tools
from svg import SVG_dotplot as svg

def execute(input_file, output_file, title="", cutoff=0, nucleotides=[], mtypes=[], motifs=[], flg_exclude_strand="off",
            filtered_regions=[], width=0, height=0, plot_settings={}):
    nucleotides = list(map(lambda s: s.strip().upper(), nucleotides.split(',')))
    mtypes = list(map(lambda s: s.strip(), mtypes.split(',')))
    motifs_to_include, motifs_to_exclude = parse_motifs(motifs)
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
    
    # Filter GFF records by score cutoff values
    if cutoff:
        dots = list(filter(lambda dot: int(dot['score']) >= cutoff, dots))
        
    # Filter sites from excluded regions
    if filtered_regions:
        dots = filter_regions(sites=dots, regions=filtered_regions)
            
    # Filter GFF records by modified nucleotides
    if nucleotides and nucleotides[0]:
        dots = list(filter(lambda dot: dot['data']['nucleotide'] in nucleotides, dots))
        
    # Filter GFF records by types of methylation
    if mtypes and mtypes[0]:
        dots = list(filter(lambda dot: dot['modification'] in mtypes, dots))
        
    # Exclude methylation on the reverse-complement strand
    if flg_exclude_strand.upper() in ('LEADING', '+'):
        dots = [d for d in dots if d['strand'] == "+"]
    if flg_exclude_strand.upper() in ('LAGGING', '-'):
        dots = [d for d in dots if d['strand'] == "-"]

    # Filter GFF records by motifs
    if motifs_to_include:
        dots = filter_motifs(ls=dots, motifs=motifs_to_include, flg_include=True)
    if motifs_to_exclude:
        dots = filter_motifs(ls=dots, motifs=motifs_to_exclude, flg_include=False)
        
    # Set maximal X (coverage) and Y (score) values
    try:
        if not width:
            width = max(list(map(lambda i: int(dots[i]['data']['coverage']), range(len(dots)))))
            width = int(str(width)[0]+("0"*(len(str(width))-1)))+10**(len(str(width))-1)
        if not height:
            height = max(list(map(lambda i: int(dots[i]['score']), range(len(dots)))))
            height = int(str(height)[0]+("0"*(len(str(height))-1)))+10**(len(str(height))-1)
    except:
        tools.alert(f"File {input_file} is either empty, or corrupted, or all modified nucleotides were filtered out!")
        return 
    
    # Set graph title
    if not title:
        title = os.path.basename(input_file[:-4])
    
    # I/O module
    oIO = seq_io.IO()
    # Generate DotPlot SVG file
    svg_data = svg(data=dots, title=title, width=width, height=height)
    svg_data.set_svg()
    #Generate DensityPlot SVG
    if plot_settings:
        oDensityPlot = GI_finder.Interface(options=plot_settings)
        oDensityPlot.set_task("MOD", flg_set_as_immutable=True)
        oDensityPlot.set_task("GC")
        oDensityPlot.set_task("GCS")
        oDensityPlot.set_task("MGE")
        height = 0
        if svg_data:
            height = svg_data.height
            
        # Receive a list of modified nucleotides/motifs in the format [[start,end,strand],...]
        sites = [[int(site['start']), int(site['end']), 1 if site['strand'] == '+' else -1] for site in dots]
        # Pass input GBK file and locations of modified nucleotides
        PlotSVG = oDensityPlot.execute(GBK_files = plot_settings['-g'], top_indend = height, modified_sites=[sites])
        svg_data += PlotSVG
        #oIO.save(PlotSVG, output_file[:-4] + "_DensPlot.svg")
    
    oIO.save(svg_data.get(), output_file)
    
def parse_motifs(motifs):
    if not motifs:
        return [], []
    # Sort out the included and excluded motifs by sign '-' in front of the motifs, like -GATC,2,-2; AGCT,1,-1
    motifs_to_include = [motif.strip() for motif in motifs.split(";") if not motif.strip().startswith("-")]
    motifs_to_exclude = [motif[1:].strip() for motif in motifs.split(";") if motif.strip().startswith("-")]
    # Split each motif by commas to motifs and numbers of modified nucleotides
    motifs_to_include = [s.split(",") for s in motifs_to_include]
    motifs_to_exclude = [s.split(",") for s in motifs_to_exclude]
    # Convert modified nucleotide locations to integers
    motifs_to_include = [[ls[0].strip().upper()] + [int(v) for v in ls[1:]] for ls in motifs_to_include]
    motifs_to_exclude = [[ls[0].strip().upper()] + [int(v) for v in ls[1:]] for ls in motifs_to_exclude]    
    return motifs_to_include, motifs_to_exclude

# Exclude modified sites within given genomic regions
def filter_regions(sites, regions):
    """
    Filters sites to exclude those with start or end positions within specified regions.

    Args:
        sites (list of dict): List of site dictionaries, each containing 'start' and 'end' keys.
        regions (list of list): List of regions, where each region is a list [start, end].

    Returns:
        list of dict: Filtered list of sites excluding those within the regions.
    """
    def is_within_regions(start, end, region_bounds):
        for region_start, region_end in region_bounds:
            if region_start <= int(start) <= region_end or region_start <= int(end) <= region_end:
                return True
        return False

    # Convert region boundaries to a list of tuples for easier comparison
    region_bounds = [(region[0], region[1]) for region in regions]

    # Filter out sites that fall within the specified regions
    filtered_sites = [
        site for site in sites
        if not is_within_regions(site['start'], site['end'], region_bounds)
    ]

    return filtered_sites
    
def filter_motifs(ls, motifs, flg_include=True):
    """
    Filters a list of records based on motifs and modified nucleotide positions.

    Args:
        ls (list): List of dictionaries, each containing a 'data' key with a 'context' field.
        motifs (list): List of motifs. Each motif is a list with the motif string and nucleotide positions.
        flg_include (bool): If True, include matching records. If False, exclude them.

    Returns:
        list: Filtered list of records.
    """
    selected_records = []
    for motif_set in motifs:
        motif_template = motif_set[0].upper()  # Motif string
        for n in motif_set[1:]:
            # Validate modified nucleotide location
            if n == 0 or abs(n) > len(motif_template):
                raise ValueError("Modified nucleotide location cannot be 0 or bigger than the motif length")

            # Reverse complement if the position is negative
            if n > 0:
                motif = motif_template
            else:
                rev_motif = tools.reverse_complement(motif_template)
                # Exclud doubling the records for palindroms
                if rev_motif == motif:
                    continue
                motif = rev_motif
                n = abs(n)
                # Do not run palindromic motifs twice
                if motif in motifs and n in motif_set[1:]:
                    continue
            # Compile motif into regex
            reg_motif = re.compile(tools.compile_motif(motif))

            # Filter records based on match results
            if flg_include:
                # Add selected records
                selected_records += [rec for rec in ls if match(rec, reg_motif, n, motif_template_length=len(motif))]
            else:
                # Remove selected records
                ls = [rec for rec in ls if not match(rec, reg_motif, n, motif_template_length=len(motif))]
    
    if flg_include == False:
        return ls
    return selected_records

def match(rec, reg_motif, n, motif_template_length):
    """
    Checks if a record's context matches the motif aligned at position n.

    Args:
        rec (dict): A record containing a 'data' key with a 'context' field (41-character string).
        reg_motif (regex): Compiled regex motif to match against the context.
        n (int): Position in the motif to align with the central nucleotide.

    Returns:
        bool: True if the context matches the motif, False otherwise.
    """
    context = rec['data']['context']  # 41-character string
    central_nucleotide_index = 20    # Central nucleotide index in the context sequence

    # Calculate the start and end positions for matching the motif
    start = central_nucleotide_index - (n - 1)
    end = start + motif_template_length  # Use the regex pattern length for alignment

    # Ensure the alignment falls within the bounds of the context
    if start < 0 or end > len(context):
        return False

    # Extract the fragment of the context sequence
    context_fragment = context[start:end]

    # Match the fragment against the motif regex
    return bool(reg_motif.match(context_fragment))
                
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

##############################################################################
if __name__ == "__main__":
    __version__ = "GFF_dotplot 26.11.2024"
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

    execute(input_file=INPUT_FILE, output_file=OUTPUT_FILE, title=args['-t'],
        cutoff=SCORE_CUTOFF, nucleotides=NUCLEOTIDES, mtypes=MTYPES, width=WIDTH, height=HEIGHT)
