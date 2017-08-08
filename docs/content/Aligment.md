# Mouse TCR Analysis Workflow 

This page contains an end-to-end workflow `(.fastq.gz -> mgatk)` for the 
analysis of scRNA-Seq data with `mgatk`. In brief, we show to download data
from EBI, align with `STAR` using recommended settings, and finally process
the mitochondrial genome data with `mgatk`.

### Download .fastq.gz data

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/006/ERR1146416/ERR1146416_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/006/ERR1146416/ERR1146416_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1146417/ERR1146417_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1146417/ERR1146417_2.fastq.gz
...
```

[Full download file](wget_mouseTCR.txt)

We then created a temporary file called `ERRs.txt`

```
ERR1146416
ERR1146417
...
```

### Alignment


```
sh runner_mouseTCR.sh ERRs.txt
```

where the contents of `runner_mouseTCR.sh` are as follows:

```
$ cat runner_mouseTCR.sh

#!/bin/bash

SRR_IDS=$(cat $1 |tr "\n" " ")

for SRR in $SRR_IDS
do
echo $SRR

STAR --genomeDir /data/aryee/pub/genomes/mm10/STAR --readFilesIn "../fastq/${SRR}_1.fastq.gz" "../fastq/${SRR}_2.fastq.gz" --readFilesCommand zcat --outFileNamePrefix ${SRR}
samtools view -H "${SRR}Aligned.out.sam" > "${SRR}.sam"
awk '$3 == "chrM" {print $0}' "${SRR}Aligned.out.sam" >> "${SRR}.sam"
samtools view -Sb "${SRR}.sam" | samtools sort > "${SRR}.mito.bam" && samtools index "${SRR}.mito.bam"
samtools view  "${SRR}.mito.bam" | wc -l > "${SRR}.mitoreads.txt"

rm "${SRR}Aligned.out.sam"
rm "${SRR}.sam"
rm "${SRR}Log.out"
rm "${SRR}SJ.out.tab"
rm "${SRR}Log.progress.out"

done

```

[Full alignment/processing file](runner_mouseTCR.sh)

**Note:** one will want to place the correct location of the `mm10` STAR genome build. Moreover,
one will want to verify that "chrM" is the correct mitochondrial chromosome (other builds
may have 'MT' or a similar chromosome). 

<br><br>