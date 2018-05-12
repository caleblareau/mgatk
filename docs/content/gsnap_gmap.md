## Aligning with GMAP/GSNAP

0) Install GMAP/GSNAP
available at http://research-pub.gene.com/gmap/

Recall: GMAP is for RNA; GSNAP is for DNA

```
# Download
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-03-25.tar.gz
tar -zxvf  gmap-gsnap-2018-03-25.tar.gz 

# Install
./configure --prefix=`pwd`
make && make install

```

1) Build reference genome with the mitochondrial genome specified as circular
```
gmap_build -d hg38 hg38.fasta -c chrM -D /data/aryee/pub/genomes/gmap_gsnap
```

2) Sample alignment

```
sample=$1

gsnap --nthreads 4 --gunzip -D /data/aryee/pub/genomes/gmap_gsnap -d hg38 "fastq_trimmed/${sample}_1.trim.fastq.gz" "fastq_trimmed/${sample}_2.trim.fastq.gz" -A sam | samtools view -bS - | samtools sort -@ 4 - -o "${sample}.st.bam"

```


**Notes:**

```
XC: Indicates whether the alignment crosses over the origin of a
circular chromosome.  If so, the string XC:A:+ is printed.

samtools view  miseq_bulk1.gsnap.bam | grep XC:A:+
```
