## Verifying your installation

After installing mgatk (`pip install mgatk`, or `pip install -e .` from a source
checkout), a few quick commands will confirm the CLI and its dependencies
(R + Bioconductor packages, snakemake) are wired up correctly. Run these
from the repository root:

```
mgatk --version
mgatk support
mgatk check -i tests/humanbam/MGH60-P6-A11.mito.bam -o /tmp/mgatk_check -n check -bt CB -g hg19
```

- `mgatk --version` confirms the package installed and the CLI entry point resolves.
- `mgatk support` lists the built-in reference genomes (rCRS, hg19, hg38, mm10, etc.) with no other dependencies required.
- `mgatk check` validates a real `.bam` file end-to-end against a reference genome without running the full pipeline; a passing run ends with `mgatk check passed!`.

To confirm the full pipeline (including duplicate removal and final `.rds`
generation with `R`), run a real `call` on the bundled test data:

```
mgatk call -i tests/humanbam -o /tmp/mgatk_call_test -n verify -g hg19
```

A successful run produces `A/C/G/T/coverage` matrices, a depth table, and an
`.rds`/`.signac.rds` file under `/tmp/mgatk_call_test/final/`.

## Running the automated test suite

```
pip install -e ".[test]"
pytest
```

The suite is split into two kinds of tests:

- **Unit tests** (`test_mgatkhelp.py`, `test_variant_calling.py`) exercise pure
  helper functions and processing logic directly; these always run and have
  no external dependencies beyond what `pip install` already provides.
- **End-to-end smoke tests** (in `test_cli.py` and `test_deletioncalling.py`)
  invoke the real `mgatk`/`mgatk-del`/`mgatk-del-find` CLIs against the bundled
  `.bam` files and check the actual output files. These require `snakemake`
  and the R/Bioconductor packages mgatk depends on (`data.table`,
  `SummarizedExperiment`, `GenomicRanges`, `Matrix`, `dplyr`); if any of those
  aren't available, the corresponding tests are skipped
  automatically (see the `requires_*` markers in `tests/conftest.py`) rather
  than failing.

## Manual usage examples

The commands below exercise each mode by hand against the bundled test data;
useful for exploring behavior or reproducing an issue outside of pytest.

### Standard usage
```
mgatk call -i humanbam -o out -n glio -g hg19
```
Here: we specify the reference genome b/c it was aligned to hg19 (16571 bp) instead of rCRS (16569), which is the default in `mgatk`

### Using bcall

There are two options: 1) known barcodes to parse 2) unknown barcodes (discover based on # of mito reads)

**Option 1**
```
mgatk bcall -i barcode/test_barcode.bam -n bc1 -o bc1d -bt CB -b barcode/test_barcodes.txt -z
mgatk tenx -i barcode/test_barcode.bam -n bc1 -o bc1dmem -bt CB -b barcode/test_barcodes.txt -c 2

```

**Option 2**
```
mgatk bcall -i barcode/test_barcode.bam -n bc2 -o bc2d -bt CB -mb 200 -z
```

### Filtering out UMI barcodes

```
mgatk bcall -i barcode/test_2.umi.bam -bt CB -z -g GRCh37 -ub UB -n test2-umi -o test2_umi
```

### Deletions

Find them
```
mgatk-del-find -i pearsonbam/CACCACTAGGAGGCGA-1.qc.bam
```

count them
```
mgatk-del -i pearsonbam -z -lc 6073,5000 -rc 13095,5000
```
