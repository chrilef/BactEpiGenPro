import os, math, time
import tools
import numpy as np

########################################################################
class SVG:
    def __init__(self, svg=[], width=1050, height=0, top_indend = 150, tmp_path="tmp", title=""):
        self.SVG = svg
        self.title = title
        self.width = width
        self.height = height
        self.top_indend = top_indend
        self.tmp_path = tmp_path
        # Output graph formats
        self.output_graph_formats = ["HTML", "SVG", "PDF", "EPS", "JPG", "PNG", "TIF", "BMP"]
        
    def __add__(self, oSVG):
        svg1 = self.inlist()
        svg2 = oSVG.inlist()
        width1,width2 = [obj.get_width() for obj in [self,oSVG]]
        height1,height2 = [obj.get_height() for obj in [self,oSVG]]
        svg1[0] = f"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 {max(width1,width2)} {max(height1,height2)}\">"
        if svg2[-1] != "</svg>":
            svg2.append("</svg>")
        svg = "\n".join(svg1[:-1]) + "\n" + "\n".join(svg2[1:])
        return SVG(svg=svg, width=max(width1,width2), top_indend=self.top_indend, title=self.title)
        
    def __str__(self):
        if self.SVG:
            return self.SVG.get()
        return ""
        
    def get(self):
        # Convert SVG to a list of lines
        svg_code_list = self.SVG if isinstance(self.SVG,list) else self.SVG.strip().split("\n")
        if not svg_code_list:
            return ""
        # Check the first line
        if not svg_code_list[0].startswith("<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\""):
            svg_code_list[0] = "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" " + svg_code_list[0][svg_code_list[0].find("viewBox="):]
        # Check the last line
        if not svg_code_list[-1].endswith("</svg>"):
            svg_code_list.append("</svg>")
        # Return SVG code as a text in multiple lines
        return "\n".join(svg_code_list)
        
    def tostring(self):
        return self.get()
        
    def inlist(self):
        return self.SVG.split("\n") if isinstance(self.SVG,str) else self.SVG
        
    def get_height(self):
        SVG = self.inlist() if isinstance(self.SVG,str) else self.SVG
        if self.SVG:
            return int(SVG[0][SVG[0].find("viewBox=\"") + len("viewBox=\""):-2].split(" ")[3])
        return 0
        
    def get_width(self):
        SVG = self.inlist() if isinstance(self.SVG,str) else self.SVG
        if SVG:
            return int(SVG[0][SVG[0].find("viewBox=\"") + len("viewBox=\""):-2].split(" ")[2])
        return 0
        
    def get_top_indend(self):
        return self.top_indend
        
    def save(self, output_file, output_format, dpi=300):
        # Check output file and folder availability
        if not (output_file and os.path.exists(os.path.dirname(output_file))):
            raise ValueError(f"Output file was not specified or output path {os.path.dirname(output_file)} does not exists!")
           
        # Generate and check SVG code
        svg_code = self.get()
        if not svg_code:
            raise ValueError(f"SVG code was not generated!")
        
        # Check output format
        output_format = output_format.strip().upper()
        if output_format not in self.output_graph_formats:
            raise ValueError(f"Unknown format {output_format} for graphic saving!")
        
        if output_format == "SVG":
            return self.save_svg(output_file, svg_code)
        else:
            success = self.convert_SVG(svg_code, output_file, output_format, dpi)
            # If graph was not save in the selected format, it is saving as HTML, then as SVG
            if not success:
                tools.msg("Attempting to save the file in HTML...")
                success = self.convert_SVG(svg_code, output_file, "HTML")
                if not success:
                    tools.msg("Attempting to save the file in SVG...")
                    return self.save_svg(output_file, svg_code)
            return success
                
    # Save SVG file
    def save_svg(self, output_file, svg_code=""):
        if not svg_code:
            svg_code = self.get()
        if not svg_code:
            raise ValueError(f"SVG code was not generated!")
            
        try:
            with open(output_file + ".svg", "w", encoding="utf-8") as f:
                f.write(svg_code)
        except Exception as e:
            tools.alert(f"Problem with saving SVG file! Error: {e}")
            return False
        return True
        
        
    # SVG conversion to alternative graphical formats
    ''' This is the main manager of the conversion functions 
        Output format can be SVG, HTML, PDF, EPS, JPG, JPEG, TIF, TIFF, PNG, BMP
    '''
    def convert_SVG(self, svg_code, output_file, output_format, dpi=300):
        tools.msg(f"Converting SVG to {output_format.upper()}...")
        if output_format.upper() == "HTML":
            return self.embed_svg_in_html(svg_code, output_file)
        elif output_format.upper() == "EPS":
            return self.convert_svg_to_eps(svg_code, output_file)
        elif output_format.upper() in ("PDF", "JPG", "JPEG", "PNG", "TIF", "TIFF", "BMP"):
            try:
                import cairosvg
            except ImportError:
                tools.alert(f"Conversion of the SVG code to {output_format.upper()} cannot be performed as the required module 'cairosvg' is not installed on this computer!")
                return False
            
            # Use high-quality raster conversion for raster formats
            if output_format.upper() in ("JPG", "JPEG", "PNG", "TIF", "TIFF", "BMP"):
                return self.convert_svg_to_high_quality_raster(svg_code, output_file, output_format, dpi)
            else:
                convertors = {
                    "PDF": self.convert_svg_to_pdf,
                    "EPS": self.convert_svg_to_eps
                }
                return convertors[output_format.upper()](svg_code, output_file)
            
    # SVG to HTML
    def embed_svg_in_html(self, svg_code, output_html):
        html_template = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{self.title}</title>
            <style>
                .svg-container {{
                    width: 50%;
                    height: {self.get_height() if self.get_height() else 500}px; /* Set desired height */
                    overflow: auto;
                    border: 1px solid black;
                }}
            </style>
        </head>
        <body>
            <div class="svg-container">
                {svg_code}
            </div>
        </body>
        </html>
        """
        try:
            with open(output_html + ".html", "w", encoding="utf-8") as f:
                f.write(html_template)
        except Exception as e:
            tools.alert(f"Problem with converting SVG to HTML! Error: {e}")
            return False
        return True
        
    def convert_svg_to_pdf(self, svg_code, output_pdf="output.pdf", dpi=None):
        """Convert SVG to PDF format."""
        try:
            import cairosvg
            cairosvg.svg2pdf(bytestring=svg_code.encode('utf-8'), write_to=output_pdf + ".pdf")
        except Exception as e:
            tools.alert(f"Problem with converting SVG to PDF! Error: {e}")
            return False
        return True
    
    import os

    def convert_svg_to_eps(self, svg_code, output_path):
        """Convert SVG to EPS (PostScript) format."""
        try:
            import cairosvg
            from PIL import Image
            from pdf2image import convert_from_path
            
            # Define tmp folder path
            os.makedirs(self.tmp_path, exist_ok=True)  # Ensure tmp folder exists
            
            # Step 1: Convert SVG to temporary PDF in tmp folder
            temp_pdf = os.path.join(self.tmp_path, f"{tools.random_id(6)}_temp.pdf")
            cairosvg.svg2pdf(bytestring=svg_code.encode('utf-8'), write_to=temp_pdf)
            
            # Check if temp PDF was created
            if not os.path.exists(temp_pdf):
                raise FileNotFoundError(f"Temporary PDF not created: {temp_pdf}")
            
            # Step 2: Convert PDF to EPS using pdf2image and Pillow
            temp_image = os.path.join(self.tmp_path, f"{tools.random_id(6)}_temp_image.png")
            pages = convert_from_path(temp_pdf, dpi=300, output_folder=self.tmp_path)
            pages[0].save(temp_image, "PNG")  # Save first page as image
            
            # Open the image and save as EPS
            img = Image.open(temp_image)
            img.save(output_path + ".eps" if not output_path.upper().endswith(".EPS") else output_path, "EPS")
            
            # Clean up temporary files
            os.remove(temp_pdf)
            os.remove(temp_image)
        except Exception as e:
            tools.alert(f"Problem with converting SVG to EPS! Error: {e}")
            return False
        return True        

    def convert_svg_to_high_quality_raster(self, svg_code, output_path, output_format, dpi=300):
        """Convert SVG to high-quality raster by first converting to EPS."""
        try:
            # Define tmp folder path
            os.makedirs(self.tmp_path, exist_ok=True)  # Ensure tmp folder exists
            
            # Convert SVG to EPS
            temp_eps = os.path.join(self.tmp_path, f"{tools.random_id(6)}_temp.eps")
            self.convert_svg_to_eps(svg_code, temp_eps)
            
            from PIL import Image
            
            # Open the EPS and convert to desired raster format
            img = Image.open(temp_eps)
            img.load(scale=2)  # Increase the scale for better quality
            
            # Convertors
            if output_format.upper() in ("JPEG", "JPG"):
                img = img.convert("RGB")  # JPEG does not support transparency
                img.save(output_path + ".jpg" if not output_path.upper().endswith(".JPG") else output_path, "JPEG", dpi=(dpi, dpi))
            elif output_format.upper() == "PNG":
                img.save(output_path + ".PNG" if not output_path.upper().endswith(".PNG") else output_path, "PNG", dpi=(dpi, dpi))
            elif output_format.upper() in ("TIF", "TIFF"):
                img.save(output_path + ".tif" if not output_path.upper().endswith(".TIF") else output_path, "TIFF", dpi=(dpi, dpi))
            elif output_format.upper() == "BMP":
                img = img.convert("RGB")  # BMP does not support transparency
                img.save(output_path + ".bmp" if not output_path.upper().endswith(".BMP") else output_path, "BMP")
            else:
                tools.alert(f"Unsupported raster format: {output_format}")
                return False

        except Exception as e:
            tools.alert(f"Problem with high-quality raster conversion! Error: {e}")
            return False
        return True
                
##############################################################################################################
class SVG_dotplot(SVG):
    def __init__(self, data={}, graph_title="", MinX=0, MaxX=0, MinY=0, MaxY=0, width=800, height=600, top_indend = 150, title="", legend="plot"):
        SVG.__init__(self, width=1000, top_indend = 150, height=500, title="")
        self.SVG = [f"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 {self.width} {self.height}\">"]
        self.graph_title = graph_title
        self.data = data
        self.MinX = float(MinX)
        self.MaxX = float(MaxX)
        self.MinY = float(MinY)
        self.MaxY = float(MaxY)
        self.radius = 5
        self.opacity = 0.8
        self.font_size = 12
        self.left_border = 80
        self.top_border = 10
        self.bottom_border = 80
        self.graph_border = 5
        self.X0 = self.left_border
        self.Y0 = self.height - self.bottom_border
        self.legend = legend    # plot | line | None
        self.dot_colors = {"m6A":"red","m4C":"darkgreen",
            "modified_base-G":"darkred","modified_base-C":"forestgreen","modified_base-A":"deeppink","modified_base-T":"blue"}
            
    # Prepare graph title
    def get_graph_title(self):
        if not self.graph_title:
            return ""
        if len(self.graph_title) <= 15:
            return self.graph_title 
        elif 15 < len(self.graph_title) <= 18:
            return f"{self.graph_title[:13]}..."
        else:
            return f"{self.graph_title[:6]}...{self.graph_title[-6:]}"
            
    def set_svg(self, data={}, MinX=0, MaxX=0, MinY=0, MaxY=0):
        def round(v):
            return 10*math.ceil(v/10)
            
        def modification(l,mod):
            if mod=="modified_base":
                mod += f"-{l}"
            return mod
            
        def _get_context(entry):
            return entry['data']['context'][14:20].lower() + entry['data']['context'][20].upper() + entry['data']['context'][21:27].lower()
        
        if data:
            self.data = data
        if MinX:
            self.MinX = float(MinX)
        if MaxY:
            self.MaxY = float(MaxY)
        if MinX:
            self.MinX = float(MinX)
        if MaxY:
            self.MaxY = float(MaxY)
            
        # Y axis
        self.SVG.append(f"<line x1=\"{self.X0}\" y1=\"{self.top_border}\" x2=\"{self.X0}\" y2=\"{self.height - self.bottom_border + 10}\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />")
        # X axis
        self.SVG.append(f"<line x1=\"10\" y1=\"{self.Y0}\" x2=\"{self.width}\" y2=\"{self.Y0}\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />")
        #### Data processing
        nucmod = list(map(lambda i: modification(self.data[i]['nucleotide'],self.data[i]['modtype']), range(len(self.data))))
        x_values = [float(v) for v in list(map(lambda i: float(self.data[i]['data']['coverage']), range(len(self.data))))]
        y_values = [float(v) for v in list(map(lambda i: float(self.data[i]['score']), range(len(self.data))))]
        # Add context titles to nodes
        titles = [f"{self.data[i]['strand']}{tools.format_numeric_string(self.data[i]['start'])}..{tools.format_numeric_string(self.data[i]['start'])} - [{_get_context(self.data[i])}]" \
            if self.data[i]['start'] != self.data[i]['end'] \
            else f"{self.data[i]['strand']}{tools.format_numeric_string(self.data[i]['start'])} - [{_get_context(self.data[i])}]" 
            for i in range(len(self.data))]
        modification_types = sorted(tools.dereplicate(nucmod))

        if not self.MinX:
            self.MinX = int(min(x_values)) - self.graph_border
        if not self.MaxX:
            self.MaxX = int(max(x_values))
            
        if not self.MinY:
            self.MinY = int(min(y_values)) - self.graph_border
        if not self.MaxY:
            self.MaxY = int(max(y_values))
        
        x_step = float(self.width-self.X0)/(self.MaxX - self.MinX + self.graph_border)
        y_step = float(self.height-self.top_border-self.bottom_border)/(self.MaxY - self.MinY + self.graph_border)
        for i in range(10):
            # Y axis
            yV = int(self.MaxY - self.MinY - i*self.MaxY/10)
            if yV < 0:
                continue
            y = self.Y0-y_step*yV
            self.SVG.append(f"<line x1=\"{self.X0 - 10}\" y1=\"{y}\" x2=\"{self.X0}\" y2=\"{y}\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />")
            self.SVG.append(f"<text x=\"{self.X0 - 10}\" y=\"{y}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size}\" style=\"text-anchor:end\">{yV + self.MinY + self.graph_border}</text>")
            # X axis
            xV = int(self.MaxX - self.MinX - i*self.MaxX/10)
            if xV < 0:
                continue
            x = self.X0+x_step*xV
            self.SVG.append(f"<line x1=\"{x}\" y1=\"{self.Y0 + 10}\" x2=\"{x}\" y2=\"{self.Y0}\" stroke=\"black\" stroke-width=\"1\" stroke-linejoin=\"round\" />")
            self.SVG.append(f"<text x=\"{x}\" y=\"{self.Y0 + 10 + self.font_size}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size}\" style=\"text-anchor:middle\">{xV + self.MinX + self.graph_border}</text>")
        
        # Add dots
        for i in range(len(x_values)):
            x = self.X0 + (x_values[i] - self.MinX) * x_step
            y = self.Y0 - (y_values[i] - self.MinY) * y_step
            if nucmod[i] in ("m6A","m4C"):
                self.SVG.append(f"<circle cx=\"{x}\" cy=\"{y}\" r=\"5\" stroke=\"black\" stroke-width=\"1\" fill=\"{self.dot_colors[nucmod[i]]}\" fill-opacity=\"{self.opacity}\">\n" +
                    f"<title>{titles[i]}</title></circle>")
            else:
                self.SVG.append(f"<circle cx=\"{x}\" cy=\"{y}\" r=\"{self.radius}\" fill=\"{self.dot_colors[nucmod[i]]}\" fill-opacity=\"{self.opacity}\">\n" +
                    f"<title>{titles[i]}</title></circle>")
        
        # Titles
        self.SVG.append(f"<text x=\"{self.X0 + 60}\" y=\"{self.top_border + self.font_size + 10}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 10}\" style=\"text-anchor:start\">{self.title}</text>")    
        self.SVG.append(f"<text x=\"{self.X0 + (self.width - self.X0)/2}\" y=\"{self.height - self.bottom_border/2}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 10}\" style=\"text-anchor:middle\">Coverage</text>")    
        self.SVG.append(f"<text x=\"{self.X0 / 2}\" y=\"{(self.height - self.top_border - self.bottom_border)/2}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 10}\" style=\"text-anchor:middle\" transform=\"rotate(-90 {self.X0 / 2} {(self.height - self.top_border - self.bottom_border)/2})\">NucMod Score</text>")  

        # Statistics plot:
        if self.legend == "plot":
            # Rectangle SVG element
            x_tab_1 = 55
            x_tab_2 = 60
            x_tab_3 = 155
            x_tab_4 = 120
            self.SVG.append("<g title=\"Legend plot\">")
            outline_attr = f"width=\"170\" height=\"{40 + 15 * len(modification_types)}\" stroke=\"#DBDBDB\" stroke-width=\"0.5\" fill=\"white\" fill-opacity=\"0.2\""
            self.SVG.append(f'<rect x="{self.left_border + 50}" y="{self.top_border}" {outline_attr} />')
            self.SVG.append(f"<text x=\"{self.left_border + x_tab_1}\" y=\"{self.top_border + 15}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">Title: {self.get_graph_title()}</text>")
            self.SVG.append(f"<text x=\"{self.left_border + x_tab_3}\" y=\"{self.top_border + 15}\" font-family=\"Times New Roman\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">{self.title}</text>")
            self.SVG.append(f"<text x=\"{self.left_border + x_tab_1}\" y=\"{self.top_border + 30}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">Sites in total:</text>")
            self.SVG.append(f"<text x=\"{self.left_border + x_tab_3}\" y=\"{self.top_border + 30}\" font-family=\"Times New Roman\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">{tools.format_numeric_string(len(self.data))}</text>")
            for i in range(len(modification_types)):
                modtype = modification_types[i]
                shift = 45 + i * 15
                mt = f"mod{modtype[-1]}" if modtype.startswith("modified_base-") else modtype
                count = len([s for s in nucmod if s == modtype])
                self.SVG.append(f"<text x=\"{self.left_border + x_tab_2}\" y=\"{self.top_border + shift}\" font-family=\"Times New Roman\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">{mt}:</text>")
                self.SVG.append(f"<circle cx=\"{self.left_border + x_tab_4}\" cy=\"{self.top_border + shift - self.radius}\" r=\"{self.radius}\" stroke=\"black\" stroke-width=\"1\" fill=\"{self.dot_colors[modification_types[i]]}\" fill-opacity=\"{self.opacity}\"/>")
                self.SVG.append(f"<text x=\"{self.left_border + x_tab_3}\" y=\"{self.top_border + shift}\" font-family=\"Times New Roman\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">{tools.format_numeric_string(count)};</text>")
                
            self.SVG.append("</g>")
        
        # Legend
        if self.legend == "line":
            self.SVG.append("<g title=\"Legend line\">")
            y = self.height - self.font_size - 10
            tab = 10
            for i in range(len(modification_types)):
                modtype = modification_types[i]
                shift = 50 + i * 100
                mt = f"mod{modtype[-1]}" if modtype.startswith("modified_base-") else modtype
                self.SVG.append(f"<circle cx=\"{self.X0 + shift}\" cy=\"{y}\" r=\"5\" stroke=\"black\" stroke-width=\"1\" fill=\"{self.dot_colors[modification_types[i]]}\" fill-opacity=\"{self.opacity}\"/>")
                self.SVG.append(f"<text x=\"{self.X0 + shift + tab}\" y=\"{y + 3}\" font-family=\"Times New Roman\" font-weight=\"bold\" font-size=\"{self.font_size + 5}\" style=\"text-anchor:start\">- {mt};</text>")
            self.SVG.append("</g>")

        self.SVG.append("</svg>")

##############################################################################################################
class SVG_linear(SVG):
    def __init__(self, seqname, seqlength, task_list, graph_title="", width=800, top_indend=150, title=""):
        SVG.__init__(self, width=width, top_indend=top_indend, title=title)
        self.seqlength = seqlength
        self.seqname = seqname
        self.task_list = task_list
        self.title = title
        self.graph_title = graph_title
        self.flg_wrongLoci = False
        self.added_tasks = 0
        self.colors = ["black","brown","blue","magenta","green","orange","pink"]
        self.main_font_size = 18
        self.large_font_size = 25
        self.small_font_size = 10
        self.width = width
        self.left_indend = 80
        self.right_indend = 5
        self.top_indend = top_indend
        self.title_height = 0
        if self.title:
            self.title_height = 20
        self.band_height = 80
        self.work_width = self.width-self.left_indend-self.right_indend
        self.work_height = self.top_indend + self.title_height + self.band_height*(len(self.task_list) - len([key for key in ("MGE","contigs") if key in self.task_list]))
        self.chromosome_width = 30
        self.legend_height = self.work_height+self.chromosome_width+150
        self.legend_step = 20
        legend_lines = len(self.task_list)
        if legend_lines < 5:
            legend_lines = 5
        self.height = self.work_height + 20
        self.right_field_width = 40
        self.length = self.work_width - self.right_field_width
        self.SVG = ""
        self.task_lines = {"graphs":{},"plots":[]}
        self.mge_boxes = []
        
    def clear_tasks(self):
        self.task_lines = {"graphs":{},"plots":[]}
        self.added_tasks = 0
        
    def X(self,pos):
        return self.left_indend + self.length*float(pos)/self.seqlength
        
    def Y(self,v,a=1.0,mode="",mean=0,stdev=0):
        if not mode or mode == "absolute":
            y = self.work_height - a*v
        elif mode == "relative":
            y = self.work_height - a*float(v-mean)/stdev
        if y > self.work_height:
            y = self.work_height
        elif y < self.top_indend:
            y = self.top_indend
        return y
        
    def add_task(self,task,task_description,windows, graphs={}, statistics={}, flg_StartEnd="entire", shift=0):
        if flg_StartEnd in ("end","entire"):
            self.added_tasks += 1

        Y = self.top_indend + self.added_tasks*self.band_height
        height = self.band_height
        self.task_lines['graphs'][task] = self.frame(loci=windows, title = task, 
            X = self.left_indend, Y = Y, width = self.length, height = self.band_height,
            flg_StartEnd = flg_StartEnd, histogram="line", graphs=graphs, statistics=statistics)
    
    # Box is created when flg_StartEnd is "start" or "entire"; when flg_StartEnd is None or "end", only path is returned;
    # histogram = line | bars
    def frame(self, loci, X, Y, width, height, title="", histogram="line", rigth_field_indend=50, right_field_top=10, flg_StartEnd=None, 
        length=0, graphs={}, statistics={}):
        # Read graphic settings
        line_color = graphs['line color'] if 'line color' in graphs else "#2560F2"
        background_color = graphs['background color'] if 'background color' in graphs else "white"
        outline_color = graphs['outline color'] if 'outline color' in graphs else "grey"
        central_line = graphs['central line'] if 'central line' in graphs else 'auto'
        maxVal = graphs['maxVal'] if 'maxVal' in graphs else 'auto'
        minVal = graphs['minVal'] if 'minVal' in graphs else 'auto'
        
        # Parse loci keys and values
        loci_parsed = []
        if isinstance(loci,dict):
            loci = list(loci.items())
        for key, value in loci:
            start, end = map(int, (key.strip().split('-') if key.find("-") > -1 else key.strip().split('..')))
            loci_parsed.append((start, end, value))
        loci_parsed.sort()  # Sort loci by start values
    
        # Calculate length if not provided
        if length == 0:
            length = max(end for _, end, _ in loci_parsed)
    
        # Determine maxVal and minVal
        if maxVal == "auto":
            maxVal = max(value for _, _, value in loci_parsed)
        if minVal == "auto":
            minVal = min(value for _, _, value in loci_parsed)
    
        # Validate maxVal and minVal
        if minVal >= maxVal:
            raise ValueError(f"minVal {minVal} must be less than maxVal {maxVal}")
    
        # Calculate average value for central_line if not provided
        if central_line == 'auto':
            central_line = sum(value for _, _, value in loci_parsed) / len(loci_parsed)
            
        if flg_StartEnd in ("start", "entire"):
            # Group element
            group = "<g>"
            if title:
                group += f"\n<title>{title}</title>"
            
            # Rectangle SVG element
            outline_attr = f'stroke="{outline_color}" stroke-width="0.5"' if outline_color else 'stroke="none"'
            rect = f'<rect x="{X}" y="{Y - height}" width="{width}" height="{height}" fill="{background_color}" {outline_attr} />\n'
            
            # Add min/max values
            rect += (f"<text x=\"{X + width + 1}\" y=\"{Y - height + right_field_top}\" font-family=\"Times New Roman\" " + 
                f"font-size=\"{self.small_font_size}\" style=\"text-anchor:start\">{maxVal:.2f}</text>\n")
            if minVal:
                rect += (f"<text x=\"{X + width + 1}\" y=\"{Y - 1}\" font-family=\"Times New Roman\" font-size=\"{self.small_font_size}\" " +
                f"style=\"text-anchor:start\">{minVal:.2f}</text>")
        
            # Central line SVG element
            y_central = Y - height * (central_line - minVal) / (maxVal - minVal)
            central_line_svg = f'<line x1="{X}" y1="{y_central}" x2="{X + width}" y2="{y_central}" stroke="grey" stroke-width="0.5" />'
            if (central_line - minVal) / (maxVal - minVal) > 0.2 and (maxVal - central_line) / (maxVal - minVal) > 0.2:
                # Add central line value
                central_line_svg += (f"\n<text x=\"{X + width + 1}\" y=\"{y_central}\" font-family=\"Times New Roman\" " + 
                    f"font-size=\"{self.small_font_size}\" style=\"text-anchor:start\">{central_line:.2f}</text>\n")
            
        if flg_StartEnd in ("end", "entire"):
            # Add task description and statistics
            shift = 15
            line = 0
            description = ""
            if title:
                description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top}\" font-family=\"Times New Roman\" " +
                    f"font-size=\"{self.main_font_size}\" style=\"text-anchor:start;font-weight:bold\">{title}</text>\n")
                line += 1

            # Read statistics
            correlation = LD = Z_score = None
            subtitles = [""]
            if statistics and all([key not in ("correlation","LD","Z-score") for key in list(statistics.keys())]):
                subtitles = list(statistics.keys())
            # Add subtitles and statistics data
            for subtitle in subtitles:
                stat = statistics[subtitle] if subtitle else statistics
                method = "Corr."
                if "correlation" in stat and stat['correlation']:
                    correlation, correlation_error, correlation_p_value = stat['correlation'][:3]
                    if len(stat['correlation']) > 3:
                        method = stat['correlation'][3]
                if "LD" in stat and stat['LD']:
                    LD, LD_error, LD_p_value = stat['LD'][:3]
                if "Z-score" in stat and stat['Z-score']:
                    Z_score, Z_error, Z_p_value = stat['Z-score'][:3]
                # Add subtitle
                if subtitle:
                    description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top  + line*shift}\" " + 
                        f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" " +
                        f"style=\"text-anchor:start;font-weight:bold\">{subtitle}</text>\n")
                    line += 1
                # Add correlation statistics
                if correlation != None:
                    if isinstance(correlation_error,list) or isinstance(correlation_error,tuple):
                        description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top + line*shift}\" " + 
                            f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" style=\"text-anchor:start\">" + 
                            f"{method}: {correlation:.2f} (ci: from {correlation_error[0]:.2f} to {correlation_error[1]:.2f});</text>\n")
                    elif correlation_error == None:
                        description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top + line*shift}\" " + 
                            f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" style=\"text-anchor:start\">" + 
                            f"{method}: {correlation:.2f}; p = {correlation_p_value:.3f};</text>\n")
                    else:
                        description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top + line*shift}\" " + 
                            f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" style=\"text-anchor:start\">" + 
                            f"{method} = {correlation:.2f} &#177; {correlation_error:.2f}; p = {correlation_p_value:.3f};</text>\n")
                    line += 1
                # Add LD statistics
                if LD != None:
                    description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top + line*shift}\" " + 
                        f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" style=\"text-anchor:start\">" + 
                        f"LD = {LD:.2f} &#177; {LD_error:.2f}; p = {LD_p_value:.3f};</text>\n")
                    line += 1
                    
                if Z_score != None:
                    description += (f"<text x=\"{X + width + rigth_field_indend}\" y=\"{Y - height + right_field_top + line*shift}\" " + 
                        f"font-family=\"Times New Roman\" font-size=\"{self.main_font_size - 2}\" style=\"text-anchor:start\">" + 
                        f"Z = {Z_score:.2f} &#177; {Z_error:.2f}; p = {Z_p_value:.3f};</text>\n")
                    line += 1
                    
                correlation = LD = Z_score = None
                        
        # Plot histogram
        path_data = []
        if histogram == "line":
            # Line histogram (path)
            for i in range(len(loci_parsed)):
                start, end, value = loci_parsed[i]
                x = X + width * ((start + end) / 2) / length
                y = Y - height * (value - minVal) / (maxVal - minVal)
                path_data.append((x, y))
        
            # Construct path element
            if len(path_data) > 1:
                path_d = f"M {path_data[0][0]},{path_data[0][1]}"
                for x, y in path_data[1:]:
                    path_d += f" L {x},{y}"
            else:
                path_d = ""
            histogram_svg = f'<path d="{path_d}" fill="none" stroke="{line_color}" stroke-width="1" />' if path_d else ""
        elif histogram == "bars":
            # Bar histogram
            for i in range(len(loci_parsed)):
                start, end, value = loci_parsed[i]
                x = X + width * ((start + end) / 2) / length
                y1 = Y - height * (central_line - minVal) / (maxVal - minVal)
                y2 = Y - height * (value - minVal) / (maxVal - minVal)
                bar_width = width * (end - start) / length
                path_data.append((x, y1, y2, bar_width, f"{value} [{start}..{end}]"))
    
            # Construct bar element
            if len(path_data) > 1:
                bars = ""
                for x, y1, y2, bar_width, title in path_data:
                    bars += (f"<line x1=\"{x}\" y1=\"{y1}\" x2=\"{x}\" y2=\"{y2}\" fill=\"none\" stroke=\"{line_color}\" stroke-width=\"{bar_width}\" stroke-linejoin=\"round\">\n" +
                        f"<title>{title}</title></line>")
            else:
                bars = ""
            histogram_svg = bars

        # Combine SVG elements based on flg_StartEnd
        if flg_StartEnd == "entire":
            svg_content = f"{group}\n{rect}\n{central_line_svg}\n{description}\n{histogram_svg}\n</g>\n"
        elif flg_StartEnd == "start":
            svg_content = f"{group}\n{rect}\n{central_line_svg}\n{histogram_svg}\n"
        elif flg_StartEnd == "end":
            svg_content = f"n{description}\n{histogram_svg}\n</g>\n"
        else:
            svg_content = f"{histogram_svg}\n"
            
        return svg_content
    
    def add_gi(self,lb,rb,color="",title=""):
        x1 = self.X(int(lb))
        x2 = self.X(int(rb))
        if x2-x1 < 1:
            x2 = x1+1
        if color:
            if color == "grey":
                self.flg_wrongLoci = True
            color = f"fill=\"{color}\""
        if not title:
            title = f"{lb}..{rb}"
        self.mge_boxes.append(
            f"<rect x=\"{x1}\" y=\"{self.work_height - (self.chromosome_width - 2) / 2}\" "
            f"width=\"{x2 - x1}\" height=\"{self.chromosome_width - 2}\" fill=\"pink\" stroke=\"none\">\n"
            f"<title>{title}</title></rect>\n"
            )
        
    def set_svg(self, verified=True):
        self.SVG = f"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 {self.width + 250} {self.height}\">\n"
        # White dot
        self.SVG += "<rect width=\"1\" height=\"1\" style=\"fill:wight;stroke-width:0;stroke-opacity:0.0\" />\n"
        # Draw cutoff plots
        if self.task_lines["plots"]:
            if len(self.task_lines["plots"]) > 1:
                self.task_lines["plots"].sort()
                self.task_lines["plots"].reverse()
            self.SVG += "\n".join([item[1] for item in self.task_lines["plots"]])+"\n"
            
        # Draw graphs
        if self.task_lines['graphs']:
            for task in self.task_lines['graphs']:
                frame = self.task_lines['graphs'][task]
                if not frame.endswith("</g>\n"):
                    frame += "</g>\n"
                self.SVG += frame
        
        # Add graph title
        if self.graph_title:
            self.SVG += (f"<text x=\"{self.left_indend}\" y=\"{self.top_indend - 10}\" style=\"text-anchor:start; \">{self.graph_title}</text>")
            
        # Draw verification
        if verified:
            self.SVG += (f"<text x=\"{self.width - self.right_indend}\" y=\"{self.top_indend - 10}\" style=\"text-anchor:end; fill:green;\">VERIFIED</text>")
        else:
            self.SVG += (f"<text x=\"{self.width - self.right_indend}\" y=\"{self.top_indend - 10}\" style=\"text-anchor:end; fill:red;\">NOT VERIFIED</text>")
                
        # Chromosome line
        self.SVG += (f"<line x1=\"{self.left_indend}\" y1=\"{self.work_height}\" x2=\"{self.left_indend + self.length}\" y2=\"{self.work_height}\" " +
            f"fill=\"none\" stroke=\"green\" stroke-width=\"{self.chromosome_width}\" stroke-linejoin=\"round\" />\n")
        # Contig lines
        line_width = 3
        line_shift = 4
        if 'contigs' in self.task_list and self.task_list['contigs']:
            for i in range(len(self.task_list['contigs']['contig borders'])):
                start = self.task_list['contigs']['contig borders'][i]['start']
                end = self.task_list['contigs']['contig borders'][i]['end']
                seqlength = self.task_list['contigs']['whole genome length']
                overall_P_values = self.task_list['contigs']['p-value']
                z_score, z_error, p_value, expectation, observation = self.task_list['contigs']['contig borders'][i]['Z-score']
                shift = line_shift if i % 2 == 0 else 2 * line_shift
                contig_X1 = self.left_indend + np.round(start * self.length / seqlength)
                contig_X2 = self.left_indend + np.round(end * self.length / seqlength)
                contig_Y = self.work_height+self.chromosome_width//2+shift
                color = "green"
                title = (self.task_list['contigs']['contig borders'][i]['title'] if self.task_list['contigs']['contig borders'][i]['title'] else
                    f"contig_{i+1} [{start}..{end}]")
                title += f"; Z={z_score:.1f}, p={p_value:.2f}"
                if p_value <= 0.05:
                    if z_score > 0:
                        color = "brown"
                    else:
                        color = "blue"
                self.SVG += (f"<line x1=\"{contig_X1}\" y1=\"{contig_Y}\" x2=\"{contig_X2}\" " + 
                    f"y2=\"{contig_Y}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"{line_width}\" stroke-linejoin=\"round\">\n" +
                    f"<title>{title}</title></line>\n")
                self.SVG += (f"<rect x=\"{contig_X1}\" y=\"{contig_Y + line_width}\" " + 
                    f"width=\"{contig_X2 - contig_X1}\" height=\"{20 - shift}\" fill=\"{color}\" fill-opacity=\"0.1\">\n" +
                    f"<title>{title}</title></rect>\n")
        else:
            self.SVG += (f"<line x1=\"{self.left_indend}\" y1=\"{self.work_height+self.chromosome_width//2+4}\" x2=\"{self.left_indend + self.length}\" " + 
                f"y2=\"{self.work_height+self.chromosome_width//2+line_shift}\" fill=\"none\" stroke=\"green\" stroke-width=\"{line_width}\" stroke-linejoin=\"round\" />\n")
        if "MOD" in self.task_list:
            task = self.task_list['MOD']
            statistics = {}
            if "MGE" in self.task_list and ('LD' in self.task_list["MGE"] or 'Z-score' in self.task_list["MGE"]):
                if 'LD' in self.task_list["MGE"]:
                    statistics["Distribution across MGE:"] = {"LD":self.task_list["MGE"]["LD"]}
                if 'Z-score' in self.task_list["MGE"]:
                    statistics["Distribution across MGE:"] = {"Z-score":self.task_list["MGE"]["Z-score"]}
            if "MOD" in self.task_list and ('LD' in self.task_list["MOD"] or 'Z-score' in self.task_list["MOD"]):
                if 'LD' in self.task_list["MOD"]:
                    key = "LD"
                else:
                    key = "Z-score"
                statistics["Distribution across CDS:"] = {key:self.task_list["MOD"][key]["coding"]}
                statistics["Distribution across non-CDS:"] = {key:self.task_list["MOD"][key]["noncoding"]}
                statistics["Distribution across promoters:"] = {key:self.task_list["MOD"][key]["promoter"]}
            if "Contigs" in self.task_list and 'LD' in self.task_list["MGE"]:
                statistics["Distribution accross contigs:"] = {"LD":self.task_list["MGE"]["LD"]}
            self.SVG += self.frame(loci=task['frames'], X=self.left_indend, Y=self.work_height - self.chromosome_width / 2 - 1,
                width=self.length, height=self.band_height, title=task['description'], histogram=task['settings']['style'], 
                flg_StartEnd='entire', graphs={}, statistics=statistics)
        
        # Draw notches and titles
        big_step = 1000000
        small_step = 100000
        if self.seqlength >= 10000000:
            small_step = 500000
        point = small_step
        while point <= self.seqlength:
            x = float(self.X(point))
            if point % big_step == 0:
                self.SVG += (f"<line x1=\"{x}\" y1=\"{self.work_height - self.chromosome_width / 2 - 10}\" x2=\"{x}\" " + 
                    f"y2=\"{self.work_height + self.chromosome_width / 2 + 10}\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" />\n")
                self.SVG += (f"<text x=\"{x}\" y=\"{self.work_height + self.chromosome_width / 2 + 50}\" font-family=\"Times New Roman\" " +
                    f"font-size=\"10\" style=\"text-anchor:middle\">{self.format_number(point, 1, -6)} Mbp</text>")
            else:
                self.SVG += (f"<circle cx=\"{x}\" cy=\"{self.work_height + self.chromosome_width / 2 + 30}\" r=\"2\" fill=\"grey\" stroke=\"black\" stroke-width=\"1.0\" />\n")
            point += small_step
            
        # Chromosome length
        self.SVG += (f"<text x=\"{self.left_indend + self.length + 50}\" y=\"{self.work_height + self.chromosome_width / 2 + 50}\" style=\"text-anchor:start\">{self.format_numeric_string(self.seqlength)} bp</text>")
        
        # Draw genomic islands
        if self.mge_boxes:
            self.SVG += "<g style=\"fill:pink;stroke:black\">\n"
            for box in self.mge_boxes:
                self.SVG += box
            self.SVG += "</g>\n"

        return self.SVG + "</svg>"
        
    def get_svg(self):
        return self.get()

    def format_number(self,num,dig,zoom=0):
        return int((10**(dig+zoom))*num)/float(10**dig)
    
    def format_numeric_string(self,num):
        outnum = str(num)
        if len(outnum) <= 3:
            return outnum
        for i in range(len(outnum)-3,0,-3):
            outnum = outnum[:i]+","+outnum[i:]
        return outnum

########################################################################
class SVG_circular(SVG):
    def __init__(self, seqname, seqlength, task_number, graph_title="", graph_format="SVG", width=0, height=0, top_indend=0, title="", tmp_path="tmp", info=[]):
        # Default width
        if not width:
            width = 550
        # Default height
        if not height:
            height = width + 25*(task_number+1)           
        SVG.__init__(self, width=width, top_indend=top_indend, height=height, title=title, tmp_path=tmp_path)
        self.seqlength = seqlength
        self.seqname = seqname
        self.graph_title = graph_title
        self.graph_format = graph_format
        self.task_number = task_number
        self.info = info
        self.flg_wrongLoci = False
        self.main_font_size = 18
        self.large_font_size = 25
        self.small_font_size = 10
        self.added_tasks = 0
        self.colors = ["black","blue","red","green","brown","orange","pink","magenta"]
        self.indend = 35
        self.r = (self.width - 2.0*self.indend)/2.0
        self.center = (self.width-2.0*self.indend)/2.0+self.indend
        self.task_lines = {"graphs":{},"circles":""}
        self.mge_boxes = []
        self.cycle_start = 2.0
        
    def _format_title(self):
        if not self.title or len(self.info) < 2 or self.graph_format.upper() not in ("SVG", "HTML"):
            return self.title
        direct = []
        reverse = []
        if self.info[0] != "":
            direct = list(map(lambda v: int(v)-1, self.info[0].split(",")))
        if self.info[1] != "":
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
        band = self.r/2.0/(self.task_number+1.5)
        mid = self.r - self.r/6.0 - band/2.0 - (self.added_tasks+1)*band
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
            (50,self.width+25*self.added_tasks,100,self.width+25*self.added_tasks,self.colors[self.added_tasks]))
        self.task_lines['graphs'][task] += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start;fill:black\">%s</text>\n" % 
            (125,self.width+3+25*self.added_tasks,task))
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
    
    def set_svg(self):
        self.SVG = f"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 {self.width} {self.height}\">\n"
        self.SVG += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
            ((self.width-2*self.indend)/2+self.indend,(self.width-2*self.indend)/2+self.indend,
            (self.width-2*self.indend)/2))
        self.SVG += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.width-2*self.indend)/2+self.indend,self.indend-20,
            (self.width-2*self.indend)/2+self.indend,self.indend+20))
        self.SVG += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.width-self.indend-20,(self.width-2*self.indend)/2+self.indend,
            self.width-self.indend+20,(self.width-2*self.indend)/2+self.indend))
        self.SVG += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.width-2*self.indend)/2+self.indend,self.width-self.indend-20,
            (self.width-2*self.indend)/2+self.indend,self.width-self.indend+20))
        self.SVG += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.indend-20,(self.width-2*self.indend)/2+self.indend,
            self.indend+20,(self.width-2*self.indend)/2+self.indend))
        # Add graph title
        self.SVG += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:end\" font-size=\"15\">%s</text>" %
            (self.width+self.indend,self.small_font_size,self.seqname))
        self.SVG += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s</text>" %
            ((self.width-2*self.indend)/2+self.indend+10,self.indend-10,tools.format_numeric_string(self.seqlength)))
        self.SVG += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s Mbp</text>" %
            (self.width-self.indend+10,(self.width-2*self.indend)/2+self.indend-10,tools.format_number(self.seqlength/4.0,2,-6)))
        self.SVG += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s Mbp</text>" %
            ((self.width-2*self.indend)/2+self.indend+10,self.width-self.indend+15,tools.format_number(self.seqlength/2.0,2,-6)))
        self.SVG += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:end\">%s Mbp</text>" %
            (self.indend-10,(self.width-2*self.indend)/2+self.indend-10,tools.format_number(3.0*self.seqlength/4.0,2,-6)))
        dot = 100000
        while dot < self.seqlength:
            a = math.pi/2.0 - 2.0*math.pi*dot/self.seqlength
            if dot%1000000 == 0:
                s = 3.0
                fill_color = "black"
            else:
                s = 2.5
                fill_color = "grey"
            self.SVG += ("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"%s\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
                (self.center+self.r*math.cos(a),self.center-self.r*math.sin(a),s,fill_color))
            dot += 100000
        if self.task_lines["circles"]:
            self.SVG += "\n"+self.task_lines["circles"]+"\n"
        if self.task_lines['graphs']:
            self.SVG += "<g style=\"fill:none;stroke-linejoin:round;stroke-miterlimit:10\">\n"
            for task in self.task_lines['graphs']:
                self.SVG += self.task_lines['graphs'][task]
            self.SVG += "</g>\n"
        if self.mge_boxes:
            self.SVG += "<g style=\"fill:pink;stroke:black\">\n"
            for box in self.mge_boxes:
                self.SVG += box
            self.SVG += "</g>\n"
        if len(self.info) > 2:
            self.SVG += ("<text x=\"275\" y=\"265\" style=\"text-anchor:middle\" font-size=\"%d\" font-weight=\"bold\">%s</text>\n" % 
                (self.main_font_size,self.info[-1]))
        font_size = self.large_font_size
        if len(self.title) > 20:
            font_size = self.small_font_size
        if len(self.title) > 10:
            font_size = self.main_font_size
        self.SVG += ("<text x=\"275\" y=\"295\" style=\"text-anchor:middle\" font-size=\"%d\" font-weight=\"bold\">%s</text>\n" % 
            (font_size,self._format_title()))
        self.SVG += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (50,self.width-25,100,self.width-25,"green"))
        self.SVG +=  ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start;fill:black\">%s</text>\n" % 
            (125,self.width-20,"Contigs"))

        if self.flg_wrongLoci:
            self.SVG += "<path  style=\"fill:grey;stroke:black\" d=\"M350,600L375,600L375,615L350,615Z\" />\n"
            self.SVG += f"<text x=\"405\" y=\"615\" style=\"text-anchor:start;fill:black;font-size:{self.main_font_size}\">Falsely selected rrn operons;</text>\n"
        self.SVG += "</svg>"
        return self.SVG
        
    def get_svg(self):
        return self.get()
