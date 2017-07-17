# Optimized adaptor trimming

Compared to most conventional tools, **proatac** offers two key features that optimize
the adaptor trimming component of the workflow. First, the trimmers are written in
efficient, parallelized C code for optimal performance. Secondly, the adaptor 
trimming makes no assumptions about the adaptor sequence _a priori_ compared to most
tools that require inputting the adaptor sequence ahead of time. These features ensure
an improved alignment rate and optimized computational time with little user overhead. The
current implementation of these features in **proatac** uses compiled linux and
mac binary files from the [PEAT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2)
software to perform adaptor trimming. 

## Bulk and C1-based single cell trimming.

Uses the naive trimming as described above without regard for the particular sequence. 

## Droplet-based barcode quantification and annotation

More to come but this is important. Need to split samples based on the particular barcode. 

<br><br>