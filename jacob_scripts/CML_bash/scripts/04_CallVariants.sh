#!/bin/bash

#$ -cwd
#$ -q gsa
#$ -N CallVariants
#$ -cwd
#$ -j y
#$ -e logs/CV.err
#$ -o logs/CV.out
#$ -l virtual_free=64g 

# Software
source /broad/software/scripts/useuse
reuse .samtools-1.3.1
reuse .bcftools-1.3.1 

# Directories
out=/broad/landerlab/sankaran/Mito/16_CML/filtered
calls=/broad/landerlab/sankaran/Mito/16_CML/calls

# Variables
bamlist="samplesForPileup.txt"

# Make list of QCd positions
cat $out/*.SNPs.Q0.QCpos.txt | sort | uniq > $out/QChetvars.txt
cat $out/QChetvars.txt | awk 'BEGIN {OFS=FS="\t"} {print "chrM", $1-1, $1}' > $out/QChetvars.bed

# Get sample from bam list
samtools mpileup -d8000 -r chrM -Q20 -t AD,ADF,ADR --skip-indels -b $bamlist -v > $calls/All.SNPs.vcf.gz
bcftools norm -m-both $calls/All.SNPs.vcf.gz > $calls/All.SNPs.MNP.vcf.gz
bcftools query -f '[%CHROM\t%POS\t%ALT\t%SAMPLE\tINFO/%AD\n]' $calls/All.SNPs.MNP.vcf.gz > $calls/All.SNPs.MNP.txt
cat $calls/All.SNPs.MNP.txt | sed -e 's/INFO\///g' -e 's/,/\t/g' | awk 'BEGIN {OFS=FS="\t"} {print $2, $3, $4, $6}' | gzip > $calls/All.SNPs.MNP.clean.txt



