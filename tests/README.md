
## Standard usage
```
mgatk call -i humanbam -o out -n glio
```

## Using the one
**This should give the same output as what's above
expect for some minor file names (output directory / base name)**

```
mgatk one -i humanbam/MGH60-P6-A11.mito.bam -o oneout
mgatk one -i humanbam/MGH60-P6-B01.mito.bam -o oneout
mgatk one -i humanbam/MGH97-P8-H02.mito.bam -o oneout
mgatk one -i humanbam/MGH97-P8-H03.mito.bam -o oneout

mgatk gather -i oneout -n glioOne
```

## Using bcall

There are two options: 1) known barcodes to parse 2) unknown barcodes (discover based on # of mito reads)

**Option 1**
```
mgatk bcall -i barcode/test_barcode.bam -n bc1 -o bc1d -bt DB -b barcode/test_barcodes.txt -z
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
mgatk-del-find -i humanbam/MGH97-P8-H02.mito.bam
```


