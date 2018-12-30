## Reconciling 10X scATAC with mgatk

_Caleb Lareau_

### Preprocessing with CellRanger
```
/apps/lab/aryee/cellranger-atac-1.0.0/cellranger-atac count --fastqs fastq --id test1 --sample Mix_Fix_1h --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0 --localcores 16
```

### Taking output for mgatk

```
bsub -q big-multi -n 16 mgatk bcall -i possorted_bam.bam -n CR_test1 -o CR_test1_mgatk -bt CB -b filtered_peak_bc_matrix/barcodes.tsv -c 16
```

### Comparison with bap

To compare with *bap*, which handles the reverse-complimenting differently, 
we have to take the reverse compliment of the barcodes passing the knee...

```
cat test1/outs/filtered_peak_bc_matrix/barcodes.tsv | cut -d "-" -f 1  | perl -ple 'y/ACGT/TGCA/ and $_ = reverse unless /^>/' > MixFix1h_bapFriendly_barcodes.tsv
```

Next, we can use these modified CellRanger barcodes to genotype the *bap* output

```
bsub -q big-multi -n 16 mgatk bcall -i Mix_Fix_1h_bap/mito/Mix_Fix_1h.mito.bam -n Bap_test1 -o Bap_test1_mgatk -bt XB -b MixFix1h_bapFriendly_barcodes.tsv -c 16
```

This should facilitate a direct comparison of the two pre-processing techniques w.r.t. mitochondrial coverage and genotypes.

<br><br>