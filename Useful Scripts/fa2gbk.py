# This script converts assembled genomes in FASTA format to GBK format.
# Before you continue, make sure to download python3
# Run the following command in your command line interface: pip3 install biopython
from Bio import SeqIO
#To be sure, use full file paths. 
inp_file="Your fasta file"  # your fasta file 
out_file="Your new gbk file" # the name of your new gbk file 
seq = list(SeqIO.parse(inp_file, "fasta"))
for i in seq:
    i.annotations['molecule_type'] = 'DNA'
count = SeqIO.write(seq, out_file, "genbank")
