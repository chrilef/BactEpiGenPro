import sys, os, string, math, subprocess, re
import seq_io

##############################################################################################################
class container:
    def __init__(self):
        self.container = []
    
    def __len__(self):
        return len(self.container)
    
    def __getitem__(self,key):
        if type(key)==type(0):
            if len(self) <= key:
                return
            return self.container[key]
        elif type(key)==type(""):
            for record in self.container:
                if record.title == key:
                    return record
        return 
    
    def __contains__(self,key):
        for record in self.container:
            if record.title == key:
                return True
        return False
    
    def __iter__(self):
        if not self.container:
            return iter([])
        records = []
        for record in self.container:
            records.append(record.copy())
        return iter(records)

    # Container items have to have the method copy()
    def get(self):
        container = []
        for record in self.container:
            container.append(record.copy())
        return container

    def format_title(self,title,length,space): # space - number of characters per line
        words = ("%s (%d bp)" % (title,length)).split(" ")
        title = [""]
        for word in words:
            title[-1] += word+" "
            if len(title[-1]) >= space:
                title.append("")
        return title
    
    def sort_container(self,sort_fn=None,flg_reverse=False):
        if sort_fn == None:
            self.container.sort()
        self.container.sort(key=sort_fn,reverse=flg_reverse)
    
    def remove_all(self):
        self.container = []
                    
##############################################################################################################
class BLAST(container):
    ###################################################
    ####  program in blastn,blastp,bl2seq,bl2seqp
    ####  path - executables
    ####  ref - path, sequence, blast database or fasta
    ####        file to be converted to blastdb
    ####  query - path or sequence
    ###################################################
    def __init__(self,seqtype="dna",binpath="",source_path=""):
        self.seqtype = seqtype
        self.path = binpath
        self.source_path = source_path
        self.query = ""
        self.sbjct = ""
        self.cline = ""
        container.__init__(self)
    
    def _set_cline(self,program="blast"):
        if self.seqtype == "dna":
            algorithm = "blastn"
            flg_prot = "F"
            dbtype = "nucl"
        elif self.seqtype == "protein":
            algorithm = "blastp"
            flg_prot = "T"
            dbtype = "prot"
        else:
            print()
            print("Sequences type was not specified!")
            print()
            return
        if not self.query or not self.sbjct:
            print()
            print("Either query %s or subject %s are wrong!" % (self.query,self.sbjct))
        query = os.path.join(self.source_path,self.query)
        sbjct = os.path.join(self.path,self.sbjct)
        if sys.platform == "win32":
            if program == "formatdb":
                self.cline = f"{os.path.join(self.path, 'formatdb')} -i {query} -p {flg_prot} -n {sbjct}"
            elif program == "blast":
                default = "-G 6 -E 2 -F F -q -2 -r 1 -e 1.0"
                self.cline = f"{os.path.join(self.path, 'blastall')} -p {algorithm} -i {query} -d {sbjct} {default}"
        elif sys.platform == "linux2":
            if program == "formatdb":
                self.cline = f"{os.path.join(self.path, 'makeblastdb')} -in {query} -dbtype {dbtype} -out {self.sbjct}"
            elif program == "blast":
                default = "-gapopen 6 -gapextend 2 -evalue 1.0 -dust no -soft_masking false"
                self.cline = f"{os.path.join(self.path, algorithm)} -query {query} -db {sbjct} {default}"
    def execute(self,query="",sbjct="",print_output=False):
        self.remove_all()
        self.query = query
        self.sbjct = sbjct
        self._set_cline()
        if print_output:
            print("blast:116",self.cline)
        output = self.process(self.cline)
        if not output:
            return False
        if print_output:
            print("blast:120",output)
        oParser = parser("blast",output)
        self.container = oParser()
        return True
    
    def create_db(self,fasta_file,dbname):
        self.query = fasta_file
        self.sbjct = dbname
        self._set_cline("formatdb")
        self.process(self.cline)
    
    def process(self,cline,flg_wait=False):
        process = subprocess.Popen(cline,stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=(sys.platform!="win32"))
        if flg_wait:
            process.wait()
        process.stdin.close()
        output = process.stdout.read()
        process.stdout.close()
        return output
    
    def tostring(self,e_threshold=None):
        query_summary = header = output = ""
        sbjct_summary = {}
        for record in self:
            query_name = record.title
            header += "Query: " + query_name + "\n\n"
            output += "Query: " + query_name + "\n\n"
            for i in range(len(record.descriptions)):
                info = record.descriptions[i]
                if e_threshold != "None" and e_threshold != None and e_threshold < info.e:
                    continue
                alignment = record[i]
                header += "\tSubject: " + alignment.source + ("\t%f\t%f" % (info.score,info.e)) + "\n"
                if alignment.source not in sbjct_summary:
                    sbjct_summary[alignment.source] = []
                sbjct_summary[alignment.source].append(query_name)
                output += ("\tSubject: "+alignment.title+"\n\n")
                for hsp in alignment:
                    if hsp.positives == None:
                        output += ("\tScore = %f, E-value = %f,\n\tIdentities = %f\n\n" %
                            (hsp.score,hsp.expect,hsp.identities))
                    else:
                        output += ("\tScore = %f, E-value = %f,\n\tIdentities = %f, Positives = %f\n\n" %
                            (hsp.score,hsp.expect,hsp.identities,hsp.positives))
                    point = 0
                    q_gaps = s_gaps = 0
                    while point < hsp.alignment_length:
                        end = point+60
                        if end > hsp.alignment_length:
                            end = hsp.alignment_length
                        q_substr = hsp.query[point:end].upper()
                        s_substr = hsp.sbjct[point:end].upper()
                        qs = q_gaps
                        ss = s_gaps
                        q_gaps += q_substr.count("-")
                        s_gaps += s_substr.count("-")
                        output += "%g\t\t%s\t\t%g\n" % (hsp.query_start+point-qs,q_substr,hsp.query_start+end-1-q_gaps)
                        output += "%g\t\t%s\t\t%g\n\n" % (hsp.sbjct_start+point-ss,s_substr,hsp.sbjct_start+end-1-s_gaps)
                        point += 60
                output += "\n"
            header += "\n"
        header += "#"*60
        if sbjct_summary:
            query_summary += "Subject summary:\n"
            sbjct_summary = sbjct_summary.items()
            if len(sbjct_summary) > 1:
                sbjct_summary.sort(key=lambda ls: ls[1])
            for item in sbjct_summary:
                query_summary += "\t".join(["",str(len(item[1])),item[0]])+"\n"
        return "\n\n".join([query_summary,header,output])
    
    def svg(self,X=25,Y=25,width=900,height=200,flg_finish=True):
        svg = []
        title_width=200
        title_height = 50
        title_start = width-title_width+5
        if not len(self):
            return ""
        for record in self:
            # Find the longest sequence
            max_length = record.query_length
            for sbjct in record:
                if sbjct.sbjct_length > max_length:
                    max_length = sbjct.sbjct_length
            query_width = title_width+float(width-title_width)*record.query_length/max_length
            query_svg = record.svg(X,Y,self.query_genemap,query_width,title_height,title_width,title_start,False)
            for sbjct in record:
                alignment_summary = sbjct.summarize()
                svg.append("<text x=\"%d\" y=\"%d\">Summarized score = %d; Best expectation = %f</text>" % 
                    (X,Y,alignment_summary.score,alignment_summary.expect))
                Y += 20
                if self.program == "bl2seq" and self.seqtype in ("dna","codon") and len(record)==1:
                    sbjct_width = title_width+float(width-title_width)*sbjct.sbjct_length/max_length
                    span = height-2*title_height
                    sbjct_svg = sbjct.svg(X,Y+span,self.sbjct_genemap,sbjct_width,title_height,title_width,title_start,False)
                    hsp_svg = sbjct.svg_hsps(float(width-title_width)/max_length,X,Y+title_height,width-title_width,span-title_height,False)
                else:
                    pass
                Y += height
                svg.extend([query_svg,hsp_svg,sbjct_svg])
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,Y+600))
            svg.append("</svg>")
        return "\n".join(svg)

    def get_record(self,title):
        for record in self.container:
            if record.title == title:
                return record.copy()
        return None
        
    def get_top_alignment(self):
        if len(self.container) and len(self.container[0].container):
            if len(self.container) > 1:
                self.sort_container(sort_fn=lambda aln: aln.score, flg_reverse=True)
            record = self.container[0]
            if len(record.container[0]) > 1:
                record.container[0].sort_container(sort_fn=lambda obj: obj.identities, flg_reverse=True)
            hsps = list(map(lambda hsp: hsp.copy(), record.container[0]))
            return [record.container[0].e,record.container[0].score,record.container[0].title,hsps]
        return None
    
    def get_matches(self,query_length,mismatches=0):
        hsps = []
        for record in self.container:
            for alignment in record.container:
                hsps += list(map(lambda rec: rec.copy(), list(filter(lambda hsp: query_length-hsp.identities+hsp.query.count("-")+hsp.sbjct.count("-") <= int(mismatches), alignment.container))))
        return hsps
    
##############################################################################################################
class blast_record(container):
    def __init__(self,query_length,title=""):
        self.query_length = query_length
        self.title = title
        self.top_score = 0
        self.summerized_score = 0.0
        self.descriptions = []
        container.__init__(self)
        
    def copy(self):
        record = blast_record(self.query_length,self.title)
        for description in self.descriptions:
            record.add_description(description.title,description.score,description.e)
        for alignment in self.container:
            record.add_alignment(alignment.copy())
        return record
        
    def add_description(self,title,score,e):
        description = blast_description(title,score,e)
        self.descriptions.append(description)
    
    def add_alignment(self,alignment):
        self.container.append(alignment)
    
    def svg(self,X=5,Y=25,genemap=[],width=800,height=50,title_width=100,title_start=0,flg_finish=True):
        svg = []
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%f\" />" %
            (X,Y+height/2,X+width-title_width,Y+height/2,"red",3.0))
        # GENE MAP
        c = float(width-title_width)/self.query_length
        for gene in genemap:
            try:
                lb,rb = list(map(lambda s: int(s),gene.split("-")))
            except:
                lb,rb = list(map(lambda s: int(s),gene.split("..")))
            if lb < 0:
                lb = 0
            if rb > self.query_length:
                rb = self.query_length
            if (genemap[gene]['remark'].find("hypothetical") > -1 or
                genemap[gene]['name'].find("hypothetical") > -1 or
                genemap[gene]['remark'].find("unknown") > -1 or
                genemap[gene]['name'].find("unknown") > -1):
                color = "grey"
            else:
                color = "green"
            bar_height = height/5
            shift = height/7
            if genemap[gene]['direction'] == "rev":
                shift = height-bar_height-height/7
            svg.append("<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:%s;stroke:%s\" />" %
                (X+c*lb,Y+shift,c*(rb-lb),bar_height,color,"grey"))
        # TITLE
        title = self.format_title(self.title,self.query_length,title_width/5)
        for subtitle in title:
            svg.append("<text x=\"%d\" y=\"%d\">%s</text>" % (X+title_start,Y+5,subtitle))
            Y += 15
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def sort(self):
        if len(self.container) > 1:
            self.container.sort(key=lambda aln: [aln.e,-aln.score], reverse=True)

##############################################################################################################
class blast_description:
    def __init__(self,title,score,e):
        self.title = title
        self.score = score
        self.e = e
        
##############################################################################################################
class blast_alignment(container):
    def __init__(self,sbjct_length=0,title="",score=0,e=None,hsps=[]):
        self.sbjct_length = sbjct_length
        self.title = title
        self.source = self.parse_title()
        self.score = score
        self.e = e
        container.__init__(self)
        self.container.extend(hsps)

    def copy(self):
        hsps = []
        for hsp in self.container:
            hsps.append(hsp.copy())
        alignment = blast_alignment(self.sbjct_length,self.title,self.score,self.e,hsps)
        return alignment
    
    def __add__(self,other=None):
        new_alignment = self.copy()
        if not other:
            return new_alignment
        new_alignment.title += "; "+other.title
        new_alignment.source += "; "+other.source
        new_alignment.score += other.score
        new_alignment.e = min(new_alignment.e,other.e)
        new_alignment.hsps.extend(other.hsps)
        return new_alignment

    def summarize(self):
        if not len(self.container):
            return None
        summarized_alignment = self.container[0].copy()
        for i in range(1,len(self.container)):
            summarized_alignment += self.container[i]
        return summarized_alignment
    
    def get_score(self):
        if not self.container:
            return None
        score = 0
        for hsp in self.container:
            score += hsp.score
        return score
    
    def get_expect(self):
        if not self.container:
            return None
        if len(self.container) == 1:
            return self.container[0].expect
        expect = self.container[0].expect
        for i in range(len(self.container),1):
            if self.container[i].expect < expect:
                expect = self.container[i].expect
        return expect
    
    def svg(self,X=5,Y=25,genemap=[],width=800,height=50,title_width=100,title_start=0,flg_finish=True):
        svg = []
        svg.append("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"%f\" />" %
            (X,Y+height/2,X+width-title_width,Y+height/2,"red",3.0))
        # GENE MAP
        c = float(width-title_width)/self.sbjct_length
        for gene in genemap:
            try:
                lb,rb = list(map(lambda s: int(s),gene.split("-")))
            except:
                lb,rb = list(map(lambda s: int(s),gene.split("..")))
            if lb < 0:
                lb = 0
            if rb > self.sbjct_length:
                rb = self.sbjct_length
            if (genemap[gene]['remark'].find("hypothetical") > -1 or
                genemap[gene]['name'].find("hypothetical") > -1 or
                genemap[gene]['remark'].find("unknown") > -1 or
                genemap[gene]['name'].find("unknown") > -1):
                color = "grey"
            else:
                color = "green"
            bar_height = height/5
            shift = height/7
            if genemap[gene]['direction'] == "rev":
                shift = height-bar_height-height/7
            svg.append("<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:%s;stroke:%s\" />" %
                (X+c*lb,Y+shift,c*(rb-lb),bar_height,color,"grey"))
        # TITLE
        title = self.format_title(self.title,self.sbjct_length,title_width/5)
        for subtitle in title:
            svg.append("<text x=\"%d\" y=\"%d\">%s</text>" % (X+title_start,Y+5,subtitle))
            Y += 15
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def svg_hsps(self,c,X=5,Y=25,width=800,height=150,flg_finish=True):
        svg = []
        for hsp in self:
            if not hsp.strand or hsp.strand=="Plus/Plus":
                svg.append("<path d=\"M%f,%dL%f,%dL%f,%dL%f,%dZ\" style=\"fill:%s;stroke:%s;opacity:%f\" />" %
                    (X+c*(hsp.query_start-1),Y,X+c*(hsp.query_end-1),Y,
                    X+c*(hsp.sbjct_end-1),Y+height,X+c*(hsp.sbjct_start-1),Y+height,
                    "blue","grey",80.0))
            elif hsp.strand=="Plus/Minus":
                svg.append("<path d=\"M%f,%dL%f,%dL%f,%dL%f,%dZ\" style=\"fill:%s;stroke:%s;opacity:%f\" />" %
                    (X+c*(hsp.query_start-1),Y,X+c*(hsp.query_end-1),Y,
                    X+c*(hsp.sbjct_start-1),Y+height,X+c*(hsp.sbjct_end-1),Y+height,
                    "blue","grey",80.0))
        if flg_finish:
            svg.insert(0,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewbox=\"0 0 %d %d\">" % (width,height))
            svg.append("</svg>")
        return "\n".join(svg)
    
    def parse_title(self):
        if not self.title:
            return ""
        p = self.title.find("; from ")
        if p == -1:
            return self.title
        source = self.title.replace(", complete sequence.","")
        return source[p+7:]

##############################################################################################################
class blast_hsp:
    def __init__(self,score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits):
        self.score = score
        self.expect = e
        self.alignment_length = int(aln_length)
        self.identities = identities
        self.positives = positives
        self.gaps = gaps
        self.strand = strand
        self.query_start = qlb
        self.sbjct_start = slb
        self.query_end = qrb
        self.sbjct_end = srb
        self.query = query
        self.sbjct = sbjct
        self.match = hits
    
    def __str__(self):
        return "Query [%d..%d]; Sbjct [%d..%d]" % (self.query_start,self.query_end,self.sbjct_start,self.sbjct_end)
    
    def __repr__(self):
        return "Query [%d..%d]; Sbjct [%d..%d]" % (self.query_start,self.query_end,self.sbjct_start,self.sbjct_end)
    
    def measure_distance(self,location,shift=0,flg_abs=True):
        if self.strand == "Plus/Plus":
            dist = self.sbjct_start + int(shift) - int(location)
        else:
            dist = self.sbjct_end - int(shift) - int(location)
        if flg_abs:
            return abs(dist)
        return dist
    
    def get_strand(self):
        if self.strand == "Plus/Plus":
            return "+"
        return "-"
    
    def copy(self):
        hsp = blast_hsp(self.score,
                            self.expect,
                            self.alignment_length,
                            self.identities,
                            self.positives,
                            self.gaps,
                            self.strand,
                            self.query_start,
                            self.sbjct_start,
                            self.query_end,
                            self.sbjct_end,
                            self.query,
                            self.sbjct,
                            self.match)
        return hsp

    def __add__(self,other=None):
        new_hsp = self.copy()
        if not other:
            return new_hsp
        new_hsp.score += other.score
        new_hsp.expect = min(new_hsp.expect,other.expect)
        if other.identities:
            new_hsp.identities += other.identities
        if other.positives:
            new_hsp.positives += other.positives
        if other.gaps:
            new_hsp.gaps += other.gaps
        new_hsp.strand = ""
        new_hsp.query_start = min(new_hsp.query_start,other.query_start)
        new_hsp.query_stop = max(new_hsp.query_end,other.query_end)
        new_hsp.sbjct_start = min(new_hsp.sbjct_start,other.sbjct_start)
        new_hsp.sbjct_stop = max(new_hsp.sbjct_end,other.sbjct_end)
        new_hsp.query = ""
        new_hsp.sbjct = ""
        new_hsp.match = ""
        return new_hsp

##############################################################################################################
class parser:
    def __init__(self,program,raw_text):
        try:
            raw_text = raw_text.decode('utf-8')
        except:
            pass
        if sys.platform == "win32":
            self.oParser = win_parser(program,raw_text)
        elif sys.platform == "linux2":
            self.oParser = linux_parser(program,raw_text)
        
    def __call__(self):
        return self.oParser._parse()
     
##############################################################################################################
class sys_parser:
    def __init__(self,program,raw_text):
        self.program = program
        self.raw_text = raw_text
        self.dataset = []
        
    def _parse(self):
        self.raw_text = self.raw_text.replace("\r","")
        q = self.raw_text.find("Query= ")
        if q == -1:
            return
        while q != -1:
            t = self.raw_text.find("Query= ",q+1)
            title = self.raw_text[q+7:self.raw_text.find("\n",q)]
            query_length = int(self.raw_text[self.raw_text.find("         (",q)+10:self.raw_text.find(" letters)\n",q)].replace(",",""))
            record = blast_record(query_length,title)
            a = self.raw_text.find("Sequences producing significant alignments",q,t)
            if a == -1:
                q = self.raw_text.find("Query= ",q+1)
                continue
            else:
                a = self.raw_text.find("\n\n",a)+2
            b = self.raw_text.find(">",a)-2
            headers = self.raw_text[a:b].split("\n")
            for hit in headers:
                if not hit:
                    continue
                e,score,name = self._parse_spsa(hit)
                if e == None:
                    continue
                record.add_description(name,score,e)
                alignment = self._parse_alignment(name,self.raw_text[q:t])
                if not alignment:
                    alignment = blast_alignment()
                record.add_alignment(alignment)
            self.dataset.append(record)
            q = self.raw_text.find("Query= ",q+1)
        return self.dataset
            
    def _parse_spsa(self,hit):
        name = hit[:73]
        while name[-1] == " ":
            name = name[:-1]
        try:
            e = self._format_e(hit[73:])
            score = int(name[name.rfind(" ")+1:])
        except:
            return None,None,""
        name = name[:name.rfind(" ")].strip()
        if not name:
            name = "no name"
        return e,score,name
            
    def _parse_alignment(self,name,input):
        if len(name) > 20:
            name = name[:21]
        p = input.find(">"+name)
        if p == -1:
            return
        l = input.find("          Length = ",p)
        title = input[p+1:l]
        sbjct_length = int(input[l+19:input.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        d = input.find(">",p+1)
        if d == -1:
            d = input.find("  Database:",p+1)
        block = input[p:d]
        b = block.find(" Score = ")
        while b > -1:
            q = block.find("Query:",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:block.find(" Score = ",q)])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
            b = block.find(" Score = ",b+1)
        alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        return alignment
    
    def _parse_hsp_header(self,block):
        s = block.find(" Score = ")
        m = block.find(",   Method:")
        d = block.find(" Identities = ")
        p = block.find(" Positives = ")
        g = block.find("Gaps = ",d)
        t = block.find("Strand")
        score = float(block[s+9:block.find(" bits ")])
        if m > -1:
            e = self._format_e(block[block.find("Expect = ")+9:m])
        else:
            e = self._format_e(block[block.find("Expect = ")+9:block.find("\n",s)])
        aln_length = int(block[block.find("/",d)+1:block.find(" (",d)])
        identities = int(block[d+13:block.find("/",d)])
        positives = None
        if p > -1:
            positives = int(block[p+12:block.find("/",p)])
        gaps = 0
        if g > -1:
            gaps = int(block[g+7:block.find("/",g)])
        strand = ""
        if t > -1:
            strand = block[block.find("Plus",t):block.find("\n",t)]
            strand = strand.replace(" ","")
        return score,e,aln_length,identities,positives,gaps,strand
        
    def _parse_hsp_body(self,block):
        lines = block.split("\n")
        query = sbjct = hits = ""
        i = j = 0
        while i < len(lines):
            if i == 0:
                qlb = int(lines[i][6:lines[i].find(" ",8)])
                slb = int(lines[i+2][6:lines[i+2].find(" ",8)])
            if len(lines[i]) < 6 or lines[i][:6] != "Query:":
                i += 1
                continue
            indend = lines[i].rfind(" ",0,14)+1
            dedend = lines[i].find(" ",indend+1)
            query += lines[i][indend:dedend].upper()
            hits += lines[i+1][indend:dedend]
            sbjct += lines[i+2][indend:dedend].upper()
            j = i
            i += 3
        qrb = int(lines[j][lines[j].rfind(" "):])
        srb = int(lines[j+2][lines[j+2].rfind(" "):])
        return qlb,slb,qrb,srb,query,sbjct,hits
        
    def _format_e(self,e):  # e as string;
        if e.find(".") > -1:
            e = float(e)
        else:
            d = e.find("e-")
            try:
                x = int(e[:d])
            except:
                x = 1
            y = int(e[d+2:])
            e = x*10.0**(-y)
        return e

##############################################################################################################
class win_parser(sys_parser):
    def __init__(self,program,raw_text):
        sys_parser.__init__(self,program,raw_text)

##############################################################################################################
class linux_parser(sys_parser):
    def __init__(self,program,raw_text):
        sys_parser.__init__(self,program,raw_text)

    def _parse(self):
        self.raw_text = self.raw_text.replace("\r","")
        #### TEMP
        IO = seq_io.IO()
        IO.save(self.raw_text,"lintmp.out","text")
        if self.program in ("bl2seq","bl2seqp"):
            return self._parse_bl2seq()
        q = self.raw_text.find("Query= ")
        if q == -1:
            return
        while q != -1:
            t = self.raw_text.find("Query= ",q+1)
            l = self.raw_text.find("Length=",q+1)+7
            title = self.raw_text[q+7:self.raw_text.find("\n",q)]
            query_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
            record = blast_record(query_length,title)
            a = self.raw_text.find("Sequences producing significant alignments",q,t)
            if a == -1:
                q = self.raw_text.find("Query= ",q+1)
                continue
            else:
                a = self.raw_text.find("\n\n",a)+2
            b = self.raw_text.find(">",a)-2
            headers = self.raw_text[a:b].split("\n")
            for hit in headers:
                if not hit:
                    continue
                result = self._parse_spsa(hit)
                e,score,name = result
                if e == None:
                    continue
                record.add_description(name,score,e)
                alignment = self._parse_alignment(name,self.raw_text[q:t])
                if not alignment:
                    alignment = blast_alignment()
                record.add_alignment(alignment)
            self.dataset.append(record)
            q = self.raw_text.find("Query= ",q+1)
        return self.dataset
            
    def _parse_bl2seq(self):
        q = self.raw_text.find("Query= ")
        l = self.raw_text.find("Length=",q+1)+7
        title = self.raw_text[q+7:self.raw_text.find("\n",q)]
        query_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
        record = blast_record(query_length,title)
        alignment = self._parse_bl2seq_alignment()
        record.add_description(alignment.title,alignment.get_score(),alignment.get_expect())
        record.add_alignment(alignment)
        self.dataset.append(record)
        return self.dataset
    
    def _parse_spsa(self,hit):
        name = hit[2:78]
        while name[-1] == " ":
            name = name[:-1]
        try:
            e = self._format_e(hit[78:])
            score = int(name[name.rfind(" ")+1:])
        except:
            return None,None,""
        name = name[:name.rfind(" ")]
        while name[-1] == " ":
            name = name[:-1]
        return e,score,name
            
    def _parse_alignment(self,name,input):
        if len(name) > 20:
            name = name[:21]
        p = input.find("> "+name)
        if p == -1:
            return
        l = input.find("Length=",p)+7
        title = input[p+2:input.find("\n",p)]
        sbjct_length = int(input[l:input.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        d = input.find(">",p+1)
        if d == -1:
            d = input.find("  Database:",p+1)
        block = input[p:d]
        b = block.find(" Score = ")
        while b > -1:
            q = block.find("Query",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:block.find(" Score = ",b+1)])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
            b = block.find(" Score = ",b+1)
        alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        return alignment
            
    def _parse_bl2seq_alignment(self):
        p = self.raw_text.find("Subject=")+8
        title = self.raw_text[p:self.raw_text.find("\n",p)]
        l = self.raw_text.find("Length=",p)+7
        sbjct_length = int(self.raw_text[l:self.raw_text.find("\n",l)].replace(",",""))
        for symbol in ["\n","\r","\t","  "]:
            title = title.replace(symbol,"")
        hsps = []
        block = self.raw_text[p:]
        b = block.find(" Score = ")
        while b > -1:
            q = block.find("Query",b)
            score,e,aln_length,identities,positives,gaps,strand = self._parse_hsp_header(block[b:q])
            b = block.find(" Score = ",b+1)
            qlb,slb,qrb,srb,query,sbjct,hits = self._parse_hsp_body(block[q:b])
            hsps.append(blast_hsp(score,e,aln_length,identities,positives,gaps,strand,qlb,slb,qrb,srb,query,sbjct,hits))
        alignment = blast_alignment(sbjct_length,title,score,e,hsps)
        return alignment
