## Reconciling 10X scATAC with mgatk

_Caleb Lareau_

### Preprocessing with CellRanger

Here's a basic call for preprocessing scATAC samples on Erisone using CellRanger--

```
/apps/lab/aryee/cellranger-atac-1.0.0/cellranger-atac count \
  --fastqs fastq --id test1 --sample sample_name \
  --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0 \
  --localcores 16
```

### Taking output for mgatk

From inside the `outs` in CellRanger output, we can launch a genotyping call.

```
mgatk bcall -i possorted_bam.bam\
  -n CR_test1 -o CR_test1_mgatk -c 16\
  -bt CB -b filtered_peak_bc_matrix/barcodes.tsv
```

<br><br>
