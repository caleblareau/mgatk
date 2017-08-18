# User-defined mitochondrial genome

Given the `.fa` or `.fasta` file that was part of your reference genome for whichever
tool, one can extract only the mitochondrial genome  using the following command for
`chrM` (update your mitochondrial
genome name accordingly-- other examples may be `MT`). Here's what we did for the `mm10`
genome--

```
samtools faidx mm10.fa chrM  >> mm10.fasta
```

With the abridged `.fasta` file generated, one can supply it as the value to the 
`--mito-genome`

<br><br>