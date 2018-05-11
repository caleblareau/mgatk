
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
gsnap --gunzip -D/data/aryee/pub/genomes/gmap_gsnap -d hg38 miseq/fastq/Bulk1_R1.fastq.gz miseq/fastq/Bulk1_R2.fastq.gz -A sam | samtools view -Sbh - | samtools sort -@ 4 - -o miseq_bulk1.gsnap.bam

```


**Notes:**

```
XC: Indicates whether the alignment crosses over the origin of a
circular chromosome.  If so, the string XC:A:+ is printed.

samtools view  miseq_bulk1.gsnap.bam | grep XC:A:+
```
