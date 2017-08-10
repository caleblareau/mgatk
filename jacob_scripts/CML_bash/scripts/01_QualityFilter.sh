#!/bin/bash

#$ -cwd
#$ -q short
#$ -N QualityFilter
#$ -cwd
#$ -j y
#$ -e logs/QF.err
#$ -o logs/QF.out
#$ -t 1-2289
#$ -tc 100
#$ -l h_vmem=12g 

# Software
source /broad/software/scripts/useuse
reuse .samtools-1.3.1
reuse R-3.3
reuse .bcftools-1.3.1 

# Directories
out=/broad/landerlab/sankaran/Mito/16_CML/filtered

# Variables
bamlist="samplesForPileupFilt2.txt"
mtfasta="/broad/landerlab/sankaran/Mito/16_CML/scripts/mtDNA.fasta"

# Get sample from bam list
bam=$(cat $bamlist | awk -v ln=$SGE_TASK_ID "NR==ln")
sample=$(basename $bam .bam)
echo $sample

# Calculate depth for each unfiltered mitochondrial bam
samtools depth $bam | awk '{sum+=$3} END {print sum/16571}' > $out/$sample.SNPs.Q0.depth.txt

# Calculate BAQ for each unfiltered mitochondrial bam
samtools mpileup -d8000 -f $mtfasta -r chrM -Q0 -t AD,ADF,ADR --skip-indels -v $bam > $out/$sample.SNPs.Q0.BAQ.vcf.gz
tabix $out/$sample.SNPs.Q0.BAQ.vcf.gz
bcftools query -f '[%POS\tINFO/%I16\n]' $out/$sample.SNPs.Q0.BAQ.vcf.gz | sed -e 's/INFO\///g' -e 's/,/\t/g' > $out/$sample.SNPs.Q0.BAQ.txt

# Calculate BQ for each unfiltered mitochondrial bam
samtools mpileup -d8000 -f $mtfasta -r chrM -Q0 -B -t AD,ADF,ADR --skip-indels -v $bam > $out/$sample.SNPs.Q0.BQ.vcf.gz
tabix $out/$sample.SNPs.Q0.BQ.vcf.gz
bcftools query -f '[%POS\tINFO/%I16\n]' $out/$sample.SNPs.Q0.BQ.vcf.gz | sed -e 's/INFO\///g' -e 's/,/\t/g' > $out/$sample.SNPs.Q0.BQ.txt

# Call Rscript to plot BQ and BAQ and output high quality positions
Rscript MitoBaseQuality_RNA.R $out/$sample.SNPs.Q0



