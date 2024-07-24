#!/bin/bash
#This script creates ipdSummary.sh jobs for all samples in a given directory      
err_msg="Please enter inputs in the following order:\n
		-Path to input directory containing fasta and bam inputs. This can include multiple samples that share a common identifier 
		(please only use '.' and '_' for your file naming conventions)\n
		-Input Directory containing BAM and Fasta files \n
		(Subsets of a sample i.e. sample1_1 and sample1_2 will incremented as follows: ipdSummary.sample1.qsub & ipdSummary.sample1.2.qsub).\n
		-Path to output directory\"\n
        	

if [ $# -eq 0 ];then
echo -e $err_msg
exit 1
fi
for ((i=1; i<=2;i++)); do
arg="${!i}"
if [ -z "$arg" ] || [ "$arg" = "-h" ];then
echo -e $err_msg
exit 1
fi
done
echo "Run this script with argument -h for help"
script_dir=$(echo "$PWD")
script="ipdSummary_updated.sh"
for fasta in ${1}/*.fa;
do
base=$(basename "$fasta")
id=$(echo "$base" | sed 's/\.(fa\|fasta\|fna)$//'| cut -d '.' -f1 | cut -d '_' -f1 )
bam=$(find "${1}" -name "${id}[._]*.bam")
echo $bam
echo $fasta
output_file="${script_dir}/ipdSummary.${id}.qsub"
if [ -f "$output_file" ]; then
output_file="${script_dir}/ipdSummary.${id}.2.qsub"
id=${id}".2"
else
output_file="${script_dir}/ipdSummary.${id}.qsub"
fi
cp "$script" "$output_file"
sed -i -e "s|\$1|$fasta|g" \
       -e "s|\$2|$bam|g" \
       -e "s|\$3|$2|g" \
       -e "s|\$4|$id|g" \
       "$output_file"
done


                                  
