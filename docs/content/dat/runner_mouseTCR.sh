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
