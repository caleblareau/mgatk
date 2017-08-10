#!/bin/bash

#$ -cwd
#$ -q short
#$ -N Blacklist
#$ -cwd
#$ -j y
#$ -e logs/BL.err
#$ -o logs/BL.out
#$ -l virtual_free=32g 

# Software
source /broad/software/scripts/useuse
reuse R-3.3

# Directories
out=/broad/landerlab/sankaran/Mito/16_CML/filtered

# Variables
bamlist="samplesForPileupFilt2.txt"
outbamlist=$(basename $bamlist .txt)

# Combine all sample BAQ inputs
rm $out/All.SNPs.Q0.BAQ.txt
while read bam; do
	sample=$(basename $bam .bam)
	echo $sample
	depth=$(cut -f1 $out/$sample.SNPs.Q0.depth.txt)
	echo $depth
	if (($(echo "$depth > 100" | bc) == 1)); then
		echo $bam >> $outbamlist.Dfilt.txt
		cat $out/$sample.SNPs.Q0.BAQ.txt >> $out/All.SNPs.Q0.BAQ.txt
	fi
done < $bamlist

# Call Rscript to make blacklist
Rscript MitoBlacklist_RNA.R $out/All.SNPs.Q0



