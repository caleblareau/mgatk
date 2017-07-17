# Built-in annotation 

For several popular genome builds, we've included a variety of useful annotation files
to optimize the workflow for **proatac** users. Currently, the following reference 
genomes are supported--

```
hg19  hg38  mm9  mm10  hg19_mm9_c
```

## Mitochondrial Blacklist

Due to the well-documented affinity of the Tn5 transposase for nucleosome-free DNA, a large 
proportion of reads in ATAC, scATAC, and related assays come from mitochondrial DNA. To mitigate
the effect of the Tn5 mitochondrial affinity, we've provided a computational workflow
and annotation 

For the supported genomes in **proatac**, we've digested the mitochondrial genome into
20 base pair reads and mapped these against the nuclear DNA to determine areas that may 
be prone to false-positive signal induced by true mitochondrial reads. These blacklisted
regions are combined with the
[ENCODE Project's Blacklist Regions](https://sites.google.com/site/anshulkundaje/projects/blacklists).
Putative peaks overlapping with these loci are filtered out automatically as part of the
**proatac** workflow. A working repository to reproduce the synthetic mapping of the 
digested mitochondrial genome to the nuclear genome can be found in the
[repository here](https://github.com/buenrostrolab/mitoblacklist).

## TSS Annotation

Genomic loci associated with transcription start sites (TSS) are particularly useful for 
quantifying the efficacy of the ATAC-Seq protocol both for insert loci distributions as well
as a read enrichment. Once again, for the supported genomes, **proatac** automatically considers
a set of TSS loci and uses them for quality control and annotation purposes. A
[simple repository](https://github.com/buenrostrolab/tss-annotation) to show where
these files where coordinated from and how they were pre-processed is embedded. 

## Bedtools Genomes

For some operations, having the lengths of the chromsomes of the reference genomes is 
useful. For the supported genomes in **proatac**, we took these chromosome sizes from 
the [bedtools git page](https://github.com/arq5x/bedtools/tree/master/genomes), where
these two column files are considered "bedtools genomes."

## Missing a reference genome?

If there is a particular reference genome that you'd like to see incorporated into **proatac**,
please submit an [issue on GitHub](https://github.com/buenrostrolab/proatac/issues) or provide
a pull-request for each of the files listed above as well as some internal variables (see
the **proatac** python class for more information about what is required for a built-in
supported genome in this tool). 

<br><br>