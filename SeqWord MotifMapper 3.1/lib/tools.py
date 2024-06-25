from functools import reduce
import os, string, random

nuc_table = {
    "A":"A",
    "T":"T",
    "G":"G",
    "C":"C",
    "U":"U",
    "R":"[A,G]",
    "Y":"[C,T]",
    "S":"[G,C]",
    "W":"[A,T]",
    "K":"[G,T]",
    "M":"[A,C]",
    "B":"[C,G,T]",
    "V":"[A,C,G]",
    "D":"[A,G,T]",
    "H":"[A,C,T]",
    "N":"[A,C,G,T]",
    }
nucleotides = [["A","T"],["T","A"],["G","C"],["C","G"],["R","Y"],["Y","R"],["S","S"],
        ["W","W"],["K","M"],["M","K"],["B","V"],["V","B"],["D","H"],["H","D"],["N","N"]]

rev_translation = [["A","T"],["T","A"],["G","C"],["C","G"],["R","Y"],["Y","R"],
                    ["K","M"],["M","K"],["B","V"],["V","B"],["D","H"],["H","D"]]

def msg(msg):
    print()
    print(msg)
    print()

def alert(msg):
    print()
    print(msg)
    print()

def format_number(num,dig,zoom=0):
    return str(int((10**(dig+zoom))*num)/float(10**dig))

def format_numeric_string(num):
    outnum = str(num)
    if len(outnum) <= 3:
        return outnum
    for i in range(len(outnum)-3,0,-3):
        outnum = outnum[:i]+","+outnum[i:]
    return outnum

def format_string(text,L=30,flg_fill=False):
    if len(text) < L:
        if flg_fill:
            text += (" "*(L-len(text)))
        return text
    else:
        return text[:L-3]+"..."
    
def get_random_filename(path='',length=10):
    letters = string.ascii_lowercase
    fname = "?"
    while fname:
        fname = "~"+"".join(random.choice(letters) for i in range(length))
        if not os.path.exists(os.path.join(path,fname)):
            break
    return fname

def dereplicate(ls):
    if len(ls) < 2:
        return
    ls.sort()
    for i in range(len(ls)-1,0,-1):
        if ls[i] == ls[i-1]:
            del ls[i]
    return ls

def get_modbase_description(ft):
    if 'description' not in ft.qualifiers:
        return {}
    data = list(map(lambda item: item.split("="), ft.qualifiers['description'][0].replace(" ","").split(";")))
    return dict(zip(list(map(lambda item: item[0], data)),list(map(lambda item: item[1], data))))

def compile_motif(word=""):
    if not word:
        return "".join(table.keys())
    letters = list(word.upper())
    return "".join(list(map(lambda l: nuc_table[l], letters)))

def inlist_motifs(motif,flg_reverse_complement = False):
    ls = []
    p = 0
    wlength = len(motif)
    while p < wlength:
        if motif[p] == "[":
            s = motif.find("]",p+1)
            if not ls:
                ls = motif[p+1:s].split(",")
            else:
                ls = reduce(lambda ls1,ls2: ls1+ls2, list(map(lambda l: list(map(lambda w: w+l, ls)), motif[p+1:s].split(","))))
            p = s+1
        else:
            if not ls:
                ls = [motif[p]]
            else:
                ls = list(map(lambda w: w+motif[p], ls))
            p += 1
    if flg_reverse_complement:
        ls = list(map(lambda w: reverse_complement(w), ls))
    return ls
def generate_motif(lngth,ls):
    if not ls:
        return ""
    if len(ls)==1:
        return ls[0]
    items = nuc_table.items()
    table = dict(zip(list(map(lambda item: item[1], items)),list(map(lambda item: item[0], items))))
    table = {
        "[A]":"A",
        "[T]":"T",
        "[G]":"G",
        "[C]":"C",
        "[U]":"U",
        "[A,G]":"R",
        "[C,T]":"Y",
        "[C,G]":"S",
        "[A,T]":"W",
        "[G,T]":"K",
        "[A,C]":"M",
        "[C,G,T]":"B",
        "[A,G,T]":"D",
        "[A,C,T]":"H",
        "[A,C,G]":"V",
        "[A,C,G,T]":"N"
        }
    motif = []
    for i in range(lngth):
        motif.append(table[str(dereplicate(list(map(lambda item: item[i], ls)))).replace("'","").replace(" ","")])
    return "".join(motif)

def reverse_complement(w):
    transtable = "".maketrans("".join(list(map(lambda ls: ls[0], rev_translation))),"".join(list(map(lambda ls: ls[1], rev_translation))))
    w = list(w.upper().translate(transtable))
    w.reverse()
    return "".join(w)

def is_dna_sequence(seq,ambiguous=True):
    letters = list(map(lambda ls: ls[0], nucleotides))
    if not ambiguous:
        letters = letters[:4]
    return len(list(filter(lambda val: val==False, list(map(lambda l: l.upper() in letters, list(seq))))))==0
    
def get_site_location_relative_to_start_codon(genes,site,strand,modbase_location,reverse_modbase_location):
    if not genes or genes==["non-coding"]:
        return "NA"
    locations = []
    start,end = list(map(lambda v: int(v), site.split("..")))
    for cds in genes:
        if cds.strand==strand:
            p = modbase_location
        elif reverse_modbase_location:
            p = end-reverse_modbase_location+1
        else:
            continue
        if strand == 1:
            locations.append(start-int(cds.location.start)+p-2)
        else:
            locations.append(int(cds.location.end)-end+p-2)
    if len(locations) > 1:
        locations.sort()
        if all(list(map(lambda v: v > 0, locations))):
            locations = [min(locations)]
        elif any(list(map(lambda v: v > 0, locations))):
            locations = [max(locations)]
    return ",".join(list(map(lambda v: str(v), locations)))

###############################################################################
if __name__ == "__main__":
    print(inlist_motifs("CTRGAWA"))
