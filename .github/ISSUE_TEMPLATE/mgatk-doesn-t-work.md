---
name: mgatk doesn't work
about: The expected output files aren't there
title: ''
labels: bug
assignees: caleblareau

---

**Describe the bug**

A clear and concise description of what the bug is, including the command that you executed. 

**A summary of .log files**

Please post the output of `cat mgatk_output_folder/logs/*log`. If you check this, and it's not clear what your issue is, then I can try to help. 

**Post an ls -lRh of mgatk_output_folder**

This command is particularly useful if you run the tool in the `-z` mode. Let us see what files were successfully generated as it may inform where the bug occurred. 

**Describe the sequencing assay being analyzed**

Is the assay 10x 3' scRNA? or a nuclear prep? If so, there may be such low mitochondria that the data will be very difficult for `mgatk` to process. Smart-seq2 and/or mtscATAC-seq data are appropriate, but many others are not. Think carefully if your genomics library has a sufficient amount of reads (at a minimum 4-5% but ideally more than 20%). You can estimate this by looking at `samtools idxstats` of the .bam file that you are supplying. 

**Clarify if the execution successful on the test data provided in the repository**

This will help us narrow down whether it's a data-intrinsic problem or an installation problem.

**Additional context**

Add any other context about the problem here, including for example if you have run the tool successfully in other contexts.
