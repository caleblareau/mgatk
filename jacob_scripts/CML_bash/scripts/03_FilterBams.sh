#!/bin/bash

#$ -cwd
#$ -q short
#$ -N FilterBams
#$ -cwd
#$ -j y
#$ -e logs/FB.err
#$ -o logs/FB.out
#$ -t 1-2289
#$ -tc 100
#$ -l h_vmem=16g 

# Software
source /broad/software/scripts/useuse
reuse .samtools-1.3.1
reuse Java-1.8
picard="/seq/software/picard-public/current/picard.jar"

# Directories
out=/broad/landerlab/sankaran/Mito/16_CML/filtered

# Variables
bamlist="samplesForPileup.txt"
mtfasta="mtDNA.fasta"

# Get sample from bam list
bam=$(cat $bamlist | awk -v ln=$SGE_TASK_ID "NR==ln")
sample=$(basename $bam .bam)
echo $sample

# Keep only reads with < 4 mismatches, no indels, perfectly paired reads only, no duplicates
samtools view -H $bam > $out/$sample.header.sam
#samtools view $bam | grep -E "NM:i:(0|1|2|3)	" | grep -E "XO:i:0	" | cat $out/$sample.header.sam - | samtools view -f 0x2 -b - > $out/$sample.filtered.bam
samtools view $bam | grep -E "nM:i:(0|1|2|3|4)$" | grep -E "NH:i:1	" | cat $out/$sample.header.sam - | samtools view -b - > $out/$sample.filtered.bam && samtools index $out/$sample.filtered.bam

