
## Standard usage
```
mgatk call -i humanbam -o out -n glio
```


## Using bcall

There are two options: 1) known barcodes to parse 2) unknown barcodes (discover based on # of mito reads)

**Option 1**
```
mgatk bcall -i barcode/test_barcode.bam -n bc1 -o bc1d -bt CB -b barcode/test_barcodes.txt -z
mgatk tenx -i barcode/test_barcode.bam -n bc1 -o bc1dmem -bt CB -b barcode/test_barcodes.txt -c 2

```

**Option 2**
```
mgatk bcall -i barcode/test_barcode.bam -n bc2 -o bc2d -bt DB -mb 200 -z
```

## Filtering out UMI barcodes

```
mgatk bcall -i barcode/test_2.umi.bam -bt CB -z -g GRCh37 -ub UB -n test2-umi -o test2_umi
```

## Deletions

Find them
```
mgatk-del-find -i pearsonbam/CACCACTAGGAGGCGA-1.qc.bam
```

count them
```
mgatk-del -i pearsonbam -z -lc 6073,5000 -rc 13095,5000
```
