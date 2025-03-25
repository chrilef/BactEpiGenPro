import sys, os, string, math
from svg import SVG_circular as SVG
from functools import reduce
sys.path.append(os.path.join("..","lib"))
import seq_io, tools

########################################################################
class Main:
    def __init__(self, seqfile, title="", graph_title="", graph_format="SVG", modbases={}, 
            loci=[], tmp_folder="tmp", window_length=5000, window_step=2000): # loci = [[1,2],[100,200],...)
        # ATTRIBUTES
        self.tmp_folder = tmp_folder
        self.SVG = None
        self.oIO = seq_io.IO()
        self.input_path = "input"
        self.output_path = "output"
        self.task_list = {}
        self.tasks = {"GC-content":{},"GC-skew":{}}
        self.window = int(window_length)
        self.step = int(window_step)
        
        self.sliding_windows = [] 
        self.fdr_cutoff_p = 1
        
        self.seqfile = seqfile
        self.title = title
        self.graph_title = graph_title
        self.graph_format = graph_format
        self.loci = loci
        self.marked_locus_color = "red"
        self.marked_locus_title = "filtered"
        
        #### TODO
        self.contigs = []
        self.genes = {}
        self.operons = modbases
        self.nrps = {}

    # METHODS
    
    def get_contigs(self,fname):
        contigs = []
        gbk = self.oIO.read(fname,"genbank")
        contigs = list(filter(lambda oFT: oFT.type=="contig", gbk.features))
        return list(map(lambda feature: "%d..%d" % (feature.location.start,feature.location.end), contigs))
    
    def svg(self,info=[]):
        return self.process(info)

    # process source data files
    def process(self,info=[]):
        # get file list
        if self.seqfile[-4:].upper() == ".GBK":
            seqlist = self.getSequenceFromGBK(self.seqfile)
            self.contigs = self.get_contigs(self.seqfile)
        elif self.seqfile[-5:].upper() == ".GBFF":
            seqlist = self.getSeqFromGBFF(self.seqfile)
        else:
            seqlist = self.getSeqFromFASTA(self.seqfile)
        if not seqlist:
            return
        acc = ""
        for seqname in seqlist:
            genes = {}
            if seqlist[seqname]['dataset'] and seqlist[seqname]['dataset']:
                genes.update({oFT.get_tag(): oFT for oFT in seqlist[seqname]['dataset']})
            
            seq = seqlist[seqname]["sequence"]
            acc = "Unknown"
            if seqlist[seqname]["accession"]:
                acc = seqlist[seqname]["accession"]
            if len(seq) == 0:
                continue
            seqlength = len(seq)
            self.SVG = SVG(seqname = (seqname if not self.graph_title else self.graph_title), seqlength = seqlength, task_number = len(self.tasks), 
                title = self.title, graph_title=self.graph_title, graph_format = self.graph_format, tmp_path = self.tmp_folder, info=info)
            start = 0
            stop = start+self.window
            while stop < seqlength:
                key = "%d-%d" % (start,stop)
                for task in self.tasks:
                    self.tasks[task][key] = self.get_value(task,seq[start:stop])
                start += self.step
                stop = start+self.window
            for task in self.tasks:
                task_description = None
                if task == "GC-content":
                    task_description = {'statistics': [50.0, 3.0], 'condition': ''}
                elif task == "GC-skew":
                    task_description = {'statistics': [1.0, .1], 'condition': ''}
                self.SVG.add_task(task,task_description,self.tasks[task])
            s = 0
            color = "brown"
            for i in range(len(self.contigs)):
                contig = self.contigs[i]
                if s < 3:
                    s += 1
                else:
                    s = 0
                lb,rb = list(map(lambda v: int(v),contig.split("..")))
                self.SVG.add_gene_as_bars(lb,rb,color,"[%d-%d]" % (lb,rb),"rev",s)
                color = "green"
            s = 0
            #### Loci
            for i in range(len(self.loci)):
                lb,rb = self.loci[i]
                self.SVG.add_gi(lb,rb,self.marked_locus_color,"[%d-%d]" % (lb,rb))
            #### ModNuc
            mygenes = list(self.operons.keys())
            mygenes.sort()
            locations = []
            for i in range(len(mygenes)):
                gene = mygenes[i]
                if type(self.operons[gene])==type([]):
                    strnd = self.operons[gene][2]
                    if strnd in ("1","dir"):
                        strnd = 'dir'
                        color = "blue"
                    else:
                        strnd = 'rev'
                        color = "red"
                    if self.operons[gene][1]:
                        color = self.operons[gene][1] 
                    lb,rb = list(map(lambda v: int(v),self.operons[gene][0].split("..")))
                    locations.append([lb,rb])
                    title = self.operons[gene][-1][0].replace("\t",";")
                else:
                    color = "blue"
                    strnd = "dir"
                    lb,rb = list(map(lambda v: int(v),self.operons[gene].split("-")))
                    locations.append([lb,rb])
                    title = self.operons[gene]
                if s < 2:
                    s += 1
                else:
                    s = 0
                self.SVG.add_gene_as_triangles(lb,rb,color,title,strnd,s)
            self.SVG.set_svg()
            #svg_code = self.SVG.get_svg()
            return self.SVG
    
    def get_statistics(self,data):
        data = data.values()
        avrg = float(reduce(lambda a,b: a+b,data))/len(data)
        data[0] = data[0]**2
        return [avrg,float(reduce(lambda a,b: a+b**2,data))/len(data)-avrg**2]
    
    def get_id(self,description):
        id = description[description.rfind(":olei")+1:]
        if id.find(" ") > -1:
            return ""
        return id
    
    def getSeqFromFASTA(self,fname,names=[]):
        oSeq = self.oIO.parse(fname,"fasta")
        return {oSeq.title:{"sequence":oSeq.seq,"dataset":None}}

    def getSequenceFromGBK(self,fname):
        oGBK = self.oIO.read(fname,"genbank")
        return {oGBK.title:{"sequence":str(oGBK.Seq),"accession":oGBK.title,"dataset":list(filter(lambda oFT: oFT.type=="CDS", oGBK.features))}}

    # check if the proposed fname already exists in the database
    # if so, change the name or return the proposed name    
    def getSeqFromGBFF(self,fname,names=[]):
        oGBF = self.oIO.parse(fname,"genbank")
        return {oGBK.title: str(oGBK.Seq) for oGBK in oGBF}
        
    def get_value(self,task,seq):
        if task == "GC-content":
            return 100.0*float(seq.upper().count("G")+seq.upper().count("C"))/len(seq)
        elif task == "GC-skew":
            c = seq.upper().count("C")
            g = seq.upper().count("G")
            if not c:
                return 0
            else:
                return 1.0 + 2.0*float(g-c)/(g+c)
        else:
            return 0


