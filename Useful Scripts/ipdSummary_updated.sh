#!/bin/bash

#*** Enter your own server's job queuing parameters below:
#PBS -l nodes=1:ppn=8
#PBS -l walltime=30:00:00
#PBS -q normal
#PBS -o /nlustre/users/chrisl/msc/PacBio/stdout
#PBS -e /nlustre/users/chrisl/msc/PacBio/stderr
#PBS -N ipdSummary

#make sure to install these modules:
module load samtools-1.10
module load bowtie2-2.4.1
module load smrtlink_12.0.0.177059
module load bcftools-1.7
module load perl-5.32.0

# The next section's values are populated by call_ipdSummary.sh
inp_ref="$1"
inp_bam="$2"
base=$(basename "$inp_bam")
bam_name=$(echo "$base" | sed 's/\.bam//')
base=$(basename "$inp_ref")
ref_name=$(echo $base | sed 's/\.\(fa\|fasta\|fna\)$//')
cd "$3" #output directory		
id="$4" #new directory to store inputs and outputs
mkdir $id
cd $id
cp --no-clobber "$inp_ref" .
cp --no-clobber "$inp_bam" .

#use this command only if you experience chemistry or kinetic analysis related errors
ccs-kinetics-bystrandify \
    "${bam_name}.bam" \
    "ccs_kinetics_bystrandify.bystrand.${bam_name}.bam"
bam_name=ccs_kinetics_bystrandify.bystrand.${bam_name}
#*** Grab desired output directory (must exist)***
#make new directory for our results
#will be named after the sample id 

# Samtools indexing
samtools faidx "${ref_name}.fa"

# Create XML dataset based on 1 or multiple raw BAM files
dataset create --type SubreadSet "${bam_name}.xml" \
"${bam_name}.bam"

# Create MMI index of the reference genome
pbmm2 index ${ref_name}.fa \
"${ref_name}.mmi"

# Align and sort raw PacBio reads against the reference
pbmm2 align --sort ${ref_name}.mmi \
"${bam_name}.xml" \
"${bam_name}.alignment.bam"

# Identify modified nucleotides and store the report to GFF files
ipdSummary "${bam_name}.alignment.bam" \
    --reference "${ref_name}.fa" \
    --identify m6A,m4C \
    --gff "${bam_name}.dnamod.gff"

# Search for canonical motifs
motifMaker find --fasta ${inp_ref} \
    --gff ${bam_name}.dnamod.gff --minScore 20 \
    --output ${ref_name}.motifs.csv

#Create and save consensus sequence from PacBio read alignment
samtools index -b ${bam_name}.alignment.bam 
samtools mpileup -uf ${inp_ref} \
${bam_name}.alignment.bam | bcftools call -c | vcfutils.pl vcf2fq > ${bam_name}.fastq
