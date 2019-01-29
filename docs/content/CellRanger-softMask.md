## Increasing mito reads / uniform coverage

_Caleb Lareau_


## Overview

We consistently observe low-coverage regions in certain mtDNA regions due to homology with 
nuclear DNA. In the setting of scATAC-seq, we expect to observe significantly more reads 
truly generating from the mtDNA than the corresponding nuclear segment. Thus, the idea
is to soft-mask regions of the nuclear genome that share consistent homology with the mtDNA,
which in effect forces observed reads to the mtDNA when aligned. 

## Generating soft-mask bed file

Several common reference genomes have been provided already, but a workflow for how to 
generate your own is discussed in the blacklist repo. 

## Updating CellRanger to have a soft-masked nuclear genome

One can follow this code using the reference data from CellRanger 1.0.0 to produce 
a soft-mask 

One exception is that I couldn't figure out how to make the `.flat` and `.gdx` files,
nor could I figure out where they were important in the 10X protocol, but the pipeline
evidently does need them. So, I just keep them. 

```
mv genome.fa old_genome.fa
mv genome.fa.flat old_genome.fa.flat
mv genome.fa.gdx old_genome.fa.gdx
maskFasta -fi old_genome.fa -bed mito_blacklist.bed  -fo genome.fa

# Cleanup old files
rm old_genome.fa 
rm genome.fa.*

# Make new reference genomes
rm hg19_10X.*
bsub -q big bowtie2-build genome.fa hg19_10X
bwa index genome.fa
samtools faidx genome.fa

# Put the flat files back; these won't be soft-masked, but they aren't being used by the aligner
mv old_genome.fa.flat genome.fa.flat
mv old_genome.fa.gdx genome.fa.gdx
```

<br><br>