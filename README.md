# BactEpiGenPro - Analyse and Visualise Methylation Patterns in Bacterial Genomes.

--- Under Development --- 

This software aims to provide a user friendly platform to determine methylation patterns and their effects on gene expression in bacterial genomes. Here you can download and run these tools on your own machine or from our website at http://begp.bi.up.ac.za. Additionally, we provide other potentially useful scripts for your methylation analysis.

To visualise canonical methylation in your samples, you can run the SeqWord Motif Mapper program.

<pre>
Program SeqWord MotifMapper v3.0
Author: O. Reva (oleg.reva@up.ac.za)
Last time modified: 22 May, 2024
Usage: python run.py [-arguments]
python run.py -h / -H / --help - show this help
python run.py -v / -V / --version - show version
Arguments:
	   -u:       # <folder> input folder, 'input' by default;
	   -o:       # <folder> output folder, 'output' by default;
	   -y:       # <folder> tmp folder, 'tmp' by default;
	   -i:       # <file name> input GFF file in folder 'input'
	   -g:       # <file name> reference GBK file in folder 'input'
	   -m:       # <file name> file with masked regions in folder 'input'
	   -w:       # <word> motif, like 'CTAG'
	   -d:       # <INT,INT> comma separated modified nucleotides in forward stran
	   -r:       # <INT,INT> comma separated modified nucleotides in reverse stran
	   -c:       # <FLOAT> score cut-off, not set by default
	   -n:       # <INT> number context mismatches, 2 by default
	   -z:       # <Yes/No> allow motif mismatch, Yes by default
	   -p:       # <INT> length of upstream promoter sequence, 75 by default
	   -f:       # 'M' - show methylated motifs (default); 'U' - show unmenthylated motifs;
</pre>

Each parameter is explained below.

## 1. Dataset

You will need the following files:

- Annotated reference genome in GenBank Flat File Format (.gbk).

- Modified nucleotides in General Feature Format (.gff).

- Sequence regions to ignore in .gbk format (optional).

An example dataset, S.aureus_150.gbk and S.aureus_150.gff, can be found in the inputs folder. Additionally, you can use the figure below to guide your methylation analysis.<br>

![image](https://github.com/user-attachments/assets/0531cdb4-a1e7-4a28-be32-7a935c9b2b56)<br>

There are a few ways to obtain these files.

### GBK files:

- Find GenBank files from a public database such as the NCBI.
  
- Use a genome annotation tool such as PROKKA and specify .gbk as your preferred output format.
  
- Already have an annotated genome in FASTA format? Convert using fa2gbk.py.

### GFF files:

We highly recommend using PacBio SMRT sequencing for methylation analysis. The SMRT Link software can be used to identify modified bases (see ipdSummary.sh).

## 2. Parameters

### Basic

- Motif: A specific nucleotide sequence recognised by methyltransferases.
  
- Modified base locations: Comma separated integer list of modified base positions within the motif. This value can be specified for both the forward and reverse strands.
  
- Show methylated motifs: Show or hide methylated sites in the output.
  
- Searching mode: Defines the searching mode for methylated sites in a motif. Sites are reported if one or more nucleotides are methylated, while motifs are reported only if all specified nucleotide positions are 
  methylated.

### Advanced

- Filter regions: Sequence regions to ignore in .gbk format.
  
- Allow context mismatches: Toggle context mismatches on or off (default = on).
  
- Context mismatches: The maximum number of context sequence mismatches allowed per methylation site (default = 2).
  
- Cut-off score: Sets the strictness of the search. This is defined as the likelihood that a given nucleotide is methylated (default = unset). For example, a cut-off of 21 equates to p = 0.01.
  
- Promoter sequence length: The length of the upstream promoter sequence (default = 75 bp).

## 3. Outputs

You will receive two outputs:

1. Overview of genomic methylation sites.
   
<img src='https://github.com/user-attachments/assets/8f6f4f36-5480-40c9-8d1c-56254f623825' width="500" height="500"><br>

2. Text output including methylation motifs, locations, and gene annotations.<br>
   
<img src='https://github.com/user-attachments/assets/e12e9ced-a57a-4b81-a6e8-1f91118b8a12' ><br>

Additionally, you can plot sequencing depth against NucMod scores using GFF_dotplot.py, located in the SeqWord Motif Mapper folder.<br>

<pre>
Program GFF_dotplot.py
Author: O. Reva (oleg.reva@up.ac.za)
Last time modified: 28 July, 2024
Usage: python GFF_dotplot.py [-h] [-v] /path/my_file.gff [options]
options:
        -o: output SVG files, empty by default to save the output into
            the input folder under the input file name.
        -n: comma separated nucleotides A,C,G,T; empty by default.
        -m: comma separated methylation types m4C,m6A; empty by default.
        -t: graph title; empty by default.
        -c: score cut-off; 0 by default.
        -w: WIDTH - maximum coverage (X) value; 0 (AUTO) by default.
        -s: HEIGHT - maximum score (Y) value; 0 (AUTO) by default.
</pre><br>

<img src= 'https://github.com/user-attachments/assets/166b87f7-e252-4e73-8fdb-9ab349f52065'>

