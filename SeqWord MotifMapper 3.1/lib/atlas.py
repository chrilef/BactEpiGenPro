import sys, os, string, math
from functools import reduce
sys.path.append(os.path.join("..","lib"))
import seq_io, tools

class Main:
    def __init__(self,seqfile,title="",modbases={},loci=[]): # loci = [[1,2],[100,200],...)
        # ATTRIBUTES
        self.oIO = seq_io.IO()
        self.input_path = "input"
        self.output_path = "output"
        self.task_list = {}
        self.tasks = {"GC-content":{},"GC-skew":{}}
        self.window = 8000
        self.step = 2000
        
        self.seqfile = seqfile
        self.title = title
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
                #genes.update(seqlist[seqname]['dataset']['Gene map'])
                genes.update({oFT.get_tag(): oFT for oFT in seqlist[seqname]['dataset']})
            
            seq = seqlist[seqname]["sequence"]
            acc = "Unknown"
            if seqlist[seqname]["accession"]:
                acc = seqlist[seqname]["accession"]
            if len(seq) == 0:
                continue
            seqlength = len(seq)
            oSVG = SVG(seqname,seqlength,len(self.tasks),self.title,info)
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
                    #task_description = {'statistics': self.get_statistics(self.tasks[task]), 'condition': ''}
                    task_description = {'statistics': [1.0, .1], 'condition': ''}
                    #print self.tasks[task].keys()[0],self.tasks[task][self.tasks[task].keys()[0]]
                oSVG.add_task(task,task_description,self.tasks[task])
            s = 0
            color = "brown"
            for i in range(len(self.contigs)):
                contig = self.contigs[i]
                if s < 3:
                    s += 1
                else:
                    s = 0
                lb,rb = list(map(lambda v: int(v),contig.split("..")))
                oSVG.add_gene_as_bars(lb,rb,color,"[%d-%d]" % (lb,rb),"rev",s)
                color = "green"
            s = 0
            #### Loci
            for i in range(len(self.loci)):
                lb,rb = self.loci[i]
                oSVG.add_gi(lb,rb,self.marked_locus_color,"[%d-%d]" % (lb,rb))
            #### ModNuc
            mygenes = list(self.operons.keys())
            mygenes.sort()
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
                    title = self.operons[gene][-1][0].replace("\t",";")
                else:
                    color = "blue"
                    strnd = "dir"
                    lb,rb = list(map(lambda v: int(v),self.operons[gene].split("-")))
                    title = self.operons[gene]
                if s < 2:
                    s += 1
                else:
                    s = 0
                oSVG.add_gene_as_triangles(lb,rb,color,title,strnd,s)
            return oSVG.get_svg()
    
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

##############################################################################################################
class SVG:
    def __init__(self,seqname,seqlength,task_number,title="",info=[]):
        self.seqlength = seqlength
        self.seqname = seqname
        self.task_number = task_number
        self.title = title
        self.info = info
        self.flg_wrongLoci = False
        self.added_tasks = 0
        self.colors = ["black","blue","red","green","brown","orange","pink","magenta"]
        self.size = 550
        self.indend = 35
        self.r = (self.size - 2.0*self.indend)/2.0
        self.center = (self.size-2.0*self.indend)/2.0+self.indend
        self.svg_holder = ""
        self.task_lines = {"graphs":{},"circles":""}
        self.mge_boxes = []
        self.cycle_start = 2.0
        
    def _format_title(self):
        if not self.title or len(self.info) < 2:
            return self.title
        direct = []
        reverse = []
        if self.info[0] != "0":
            direct = list(map(lambda v: int(v)-1, self.info[0].split(",")))
        if self.info[1] != "0":
            reverse = list(map(lambda v: int(v)-1, self.info[1].split(",")))
        points = direct+reverse
        points.sort()
        title = self.title
        if points:
            if points[0] in direct:
                title = self.title[:points[0]] + ("<tspan fill='red' text-decoration='underline'>%s</tspan>" % self.title[points[0]])
            else:
                title = self.title[:points[0]] + ("<tspan fill='blue' font-style='italic'>%s</tspan>" % self.title[points[0]])
        for i in range(1,len(points),1):
            if points[i] in direct:
                title += self.title[points[i-1]+1:points[i]]+("<tspan fill='red' text-decoration='underline'>%s</tspan>" % self.title[points[i]])
            else:
                title += self.title[points[i-1]+1:points[i]]+("<tspan fill='blue' font-style='italic'>%s</tspan>" % self.title[points[i]])
        if points:
            title += self.title[points[-1]+1:]
        return title
        
        return list(map(lambda v: ("%s<tspan fill='red' text-decoration='underline'>%s</tspan>%s" % 
                                (self.title[:v-1],self.title[v-1],self.title[self.info[0]:])), points))
    def clear_tasks(self):
        self.task_lines = {"graphs":{},"circles":""}
        self.added_tasks = 0
       
    def add_task(self,task,task_description,windows):
        def get_start_location(s):
            try:
                val = int(s.split("-")[0])
            except:
                val = int(s.split("..")[0])
            return val

        if not task_description:
            task_description = {'statistics': [50.0, 5.0], 'condition': ''}
        mean,stdev = task_description["statistics"]
        band = self.r/2.0/(self.task_number+1)
        mid = self.r - self.r/6.0 - band/2.0 - self.added_tasks*band
        if not task_description["condition"]:
            self.task_lines["circles"] = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"white\" stroke=\"grey\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid))
        elif task_description["condition"] == "bigger than":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"aliceblue"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"white",band,mid,mean,stdev,float(task_description["val1"]))
        elif task_description["condition"] == "smaller than":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"white"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"aliceblue",band,mid,mean,stdev,float(task_description["val1"]))
        elif task_description["condition"] == "between":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"white"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"aliceblue",band,mid,mean,stdev,float(task_description["val1"]),float(task_description["val2"]))
        self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
            (self.center,self.center,mid-band/2.0))
        start = 0
        VAL = None
        wins = list(windows.keys())
        wins.sort(key=lambda s: get_start_location(s), reverse=True)
        self.task_lines['graphs'][task] = ""
        for win in wins:
            try:
                first,second = list(map(lambda v: float(v),win.split("-")))
            except:
                first,second = list(map(lambda v: float(v),win.split("..")))
            stop = (first+second)/2.0
            deviation = band*(windows[win]-mean)/6.0/stdev
            if VAL == None:
                VAL = deviation
                start = stop
                a = math.pi/self.cycle_start - 2.0*math.pi*start/self.seqlength
                x = self.center + (mid+VAL)*math.cos(a)
                y = self.center - (mid+VAL)*math.sin(a)
                self.task_lines['graphs'][task] += ("<path style=\"stroke:%s;stroke-width:1.0\" d=\"M%f,%f" % (self.colors[self.added_tasks],x,y))
                continue
            b = math.pi/self.cycle_start - 2.0*math.pi*stop/self.seqlength
            x = self.center + (mid+deviation)*math.cos(b)
            y = self.center - (mid+deviation)*math.sin(b)
            self.task_lines['graphs'][task] += ("L%f,%f" % (x,y))
            VAL = deviation
            start = stop
        self.task_lines['graphs'][task] += "Z\" />\n"
        self.task_lines['graphs'][task] += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (50,self.size+25*self.added_tasks,100,self.size+25*self.added_tasks,self.colors[self.added_tasks]))
        self.task_lines['graphs'][task] += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start;fill:black\">%s</text>\n" % 
            (125,self.size+3+25*self.added_tasks,task))
        self.added_tasks += 1
    
    def task_inner_circle(self,mode,fill_color,band,mid,mean,stdev,value1,value2=None):
        if mode=="absolute":
            r = mid+band*(value1-mean)/6.0/stdev
            if r > mid+band/2.0-1:
                r = mid+band/2.0-1
            #### TEMP
            r = 160
            line = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"none\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,r,fill_color))
        elif mode=="sigmas":
            r = mid+band*value1/6.0
            if r > mid+band/2.0-1:
                r = mid+band/2.0-1
            #### TEMP
            r = 160
            line = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"none\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,r,fill_color))
        elif mode=="fraction":
            pass
        if value2 != None:
            line += "\n" + self.task_inner_circle(self,mode,"white",band,mid,mean,stdev,value2)
        return line
        
    def add_gene_as_bars(self,lb,rb,color="",title="",strand="dir",i=0):
        #s_val = [15,20,25]
        s_val = [30,36,42,48]
        s = s_val[i]
        lb = int(lb)
        rb = int(rb)
        if strand == "rev":
            r1 = self.r-s
            r2 = r1-1
        else:
            r1 = self.r+s
            r2 = r1+1
        a = math.pi/2.0 - 2.0*math.pi*lb/self.seqlength
        x = self.center + r1*math.cos(a)
        y = self.center - r1*math.sin(a)
        if color:
            if color == "grey":
                self.flg_wrongLoci = True
            color = "fill=\"%s\" stroke-width=\"3\" stroke=\"%s\"" % (color,color)
        if title:
            self.mge_boxes.append("<a xlink:title=\"%s\">" % title)
        self.mge_boxes.append("<path %s d=\"M%f,%f" % (color,x,y))
        step = self.seqlength/360.0
        stop = lb + step
        b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
        while stop <= rb:
            x,y = self.get_coordinates(stop,r1)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop += step
        x,y = self.get_coordinates(stop,r2)
        self.mge_boxes[-1] += ("L%f,%f" % (x,y))
        stop -= step
        while stop >= lb:
            x,y = self.get_coordinates(stop,r2)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop -= step
        self.mge_boxes[-1] += "Z\" />\n"
        if title:
            self.mge_boxes.append("</a>")
        
    def add_gene_as_triangles(self,lb,rb,color="",title="",strand="dir",i=0):
        s_val = [15,20,25]
        s = s_val[i]
        lb = int(lb)
        rb = int(rb)
        mp = float(rb+lb)/2.0
        lb = int(math.floor(mp-float(self.seqlength)/720.0))
        rb = int(math.floor(mp+float(self.seqlength)/720.0))
        if strand == "rev":
            r1 = self.r-s
            r2 = r1+10
        else:
            r1 = self.r+s
            r2 = r1-10
        a = math.pi/2.0 - 2.0*math.pi*lb/self.seqlength
        x = self.center + r1*math.cos(a)
        y = self.center - r1*math.sin(a)
        if color:
            if color == "grey":
                self.flg_wrongLoci = True
            color = "fill=\"%s\" stroke=\"%s\"" % (color,color)
        if title:
            self.mge_boxes.append("<a xlink:title=\"%s\">" % title)
        self.mge_boxes.append("<path %s d=\"M%f,%f" % (color,x,y))
        step = self.seqlength/1800.0
        stop = lb + step
        b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
        while stop <= rb:
            x,y = self.get_coordinates(stop,r1)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop += step
        x,y = self.get_coordinates(float(rb+lb)/2.0,r2)
        self.mge_boxes[-1] += ("L%f,%f" % (x,y))
        self.mge_boxes[-1] += "Z\" />\n"
        if title:
            self.mge_boxes.append("</a>")
    
    def add_gi(self,lb,rb,color="",title=""):
        lb = int(lb)
        rb = int(rb)
        r1 = self.r-5
        r2 = self.r+5
        a = math.pi/2.0 - 2.0*math.pi*lb/self.seqlength
        x = self.center + r1*math.cos(a)
        y = self.center - r1*math.sin(a)
        if color:
            if color == "grey":
                self.flg_wrongLoci = True
            color = "fill=\"%s\"" % color
        if title:
            self.mge_boxes.append("<a href=\"\" xlink:title=\"%s\">" % title)
        self.mge_boxes.append("<path %s d=\"M%f,%f" % (color,x,y))
        step = self.seqlength/360.0
        stop = lb + step
        b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
        while stop <= rb:
            b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
            x = self.center + r1*math.cos(b)
            y = self.center - r1*math.sin(b)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop += step
        x = self.center + r2*math.cos(b)
        y = self.center - r2*math.sin(b)
        self.mge_boxes[-1] += ("L%f,%f" % (x,y))
        stop -= step
        while stop >= lb:
            b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
            x = self.center + r2*math.cos(b)
            y = self.center - r2*math.sin(b)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop -= step
        self.mge_boxes[-1] += "Z\" />\n"
        if title:
            self.mge_boxes.append("</a>")
        
    def get_coordinates(self,point,r):
        b = math.pi/2.0 - 2.0*math.pi*point/self.seqlength
        x = self.center + r*math.cos(b)
        y = self.center - r*math.sin(b)
        return x,y
    
    def get_svg(self):
        self.svg_holder = "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewbox=\"0 0 %d %d\">\n" % (self.size,self.size+25*(self.task_number+1))
        self.svg_holder += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,(self.size-2*self.indend)/2+self.indend,
            (self.size-2*self.indend)/2))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,self.indend-20,
            (self.size-2*self.indend)/2+self.indend,self.indend+20))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.size-self.indend-20,(self.size-2*self.indend)/2+self.indend,
            self.size-self.indend+20,(self.size-2*self.indend)/2+self.indend))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,self.size-self.indend-20,
            (self.size-2*self.indend)/2+self.indend,self.size-self.indend+20))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.indend-20,(self.size-2*self.indend)/2+self.indend,
            self.indend+20,(self.size-2*self.indend)/2+self.indend))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:middle\" font-size=\"15\">%s</text>" %
            (self.size/2,10,self.seqname))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s</text>" %
            ((self.size-2*self.indend)/2+self.indend+10,self.indend-10,tools.format_numeric_string(self.seqlength)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s Mbp</text>" %
            (self.size-self.indend+10,(self.size-2*self.indend)/2+self.indend-10,tools.format_number(self.seqlength/4.0,2,-6)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s Mbp</text>" %
            ((self.size-2*self.indend)/2+self.indend+10,self.size-self.indend+15,tools.format_number(self.seqlength/2.0,2,-6)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:end\">%s Mbp</text>" %
            (self.indend-10,(self.size-2*self.indend)/2+self.indend-10,tools.format_number(3.0*self.seqlength/4.0,2,-6)))
        dot = 100000
        while dot < self.seqlength:
            a = math.pi/2.0 - 2.0*math.pi*dot/self.seqlength
            if dot%1000000 == 0:
                s = 3.0
                fill_color = "black"
            else:
                s = 2.5
                fill_color = "grey"
            self.svg_holder += ("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"%s\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
                (self.center+self.r*math.cos(a),self.center-self.r*math.sin(a),s,fill_color))
            dot += 100000
        if self.task_lines["circles"]:
            self.svg_holder += "\n"+self.task_lines["circles"]+"\n"
        if self.task_lines['graphs']:
            self.svg_holder += "<g style=\"fill:none;stroke-linejoin:round;stroke-miterlimit:10\">\n"
            for task in self.task_lines['graphs']:
                self.svg_holder += self.task_lines['graphs'][task]
            self.svg_holder += "</g>\n"
        if self.mge_boxes:
            self.svg_holder += "<g style=\"fill:pink;stroke:black\">\n"
            for box in self.mge_boxes:
                self.svg_holder += box
            self.svg_holder += "</g>\n"
        if len(self.info) > 2:
            self.svg_holder += ("<text x=\"275\" y=\"265\" style=\"text-anchor:middle\" font-size=\"%d\" font-weight=\"bold\">%s</text>\n" % 
                (18,self.info[-1]))
        font_size = 25
        if len(self.title) > 20:
            font_size = 10
        if len(self.title) > 10:
            font_size = 18
        self.svg_holder += ("<text x=\"275\" y=\"295\" style=\"text-anchor:middle\" font-size=\"%d\" font-weight=\"bold\">%s</text>\n" % 
            (font_size,self._format_title()))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (50,self.size-25,100,self.size-25,"green"))
        self.svg_holder +=  ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start;fill:black\">%s</text>\n" % 
            (125,self.size-20,"Contigs"))

        if self.flg_wrongLoci:
            self.svg_holder += "<path  style=\"fill:grey;stroke:black\" d=\"M350,600L375,600L375,615L350,615Z\" />\n"
            self.svg_holder += "<text x=\"405\" y=\"615\" style=\"text-anchor:start;fill:black\">Falsely selected rrn operons;</text>\n"
        return self.svg_holder + "</svg>"

##############################################################################################################
