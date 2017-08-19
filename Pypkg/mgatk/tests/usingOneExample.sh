#!/bin/bash

# Standard usage
mgatk call -i humanbam -o out -n glio

# This should give the same output as what's above
# expect for some minor file names (output directory / base name)

mgatk one -i humanbam/MGH60-P6-A11.mito.bam -o oneout
mgatk one -i humanbam/MGH60-P6-B01.mito.bam -o oneout
mgatk one -i humanbam/MGH97-P8-H02.mito.bam -o oneout
mgatk one -i humanbam/MGH97-P8-H03.mito.bam -o oneout

mgatk gather -i oneout -n glioOne


