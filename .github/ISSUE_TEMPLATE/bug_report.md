---
name: Bug report
about: Report a crash, installation problem, or unexpected behavior
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what went wrong.

**Command run**
The exact `mgatk`/`mgatk-del`/`mgatk-del-find` command you ran, including all flags.

**Error output / traceback**
Paste the full error message or traceback here (use a code block).

**Expected behavior**
What you expected to happen instead.

**Environment**
 - `mgatk --version` output:
 - Install method: [pip / conda / from source]
 - OS: [e.g. Ubuntu 22.04, macOS 14]
 - Python version: [`python --version`]
 - Relevant tool versions if applicable: `R --version`, `snakemake --version`

**Sanity check**
Does `mgatk check -i <your.bam> -o <out> -n check -bt <tag> -g <genome>` pass? See [tests/README.md](../../tests/README.md#verifying-your-installation) for a few quick commands to confirm your install is healthy before digging further.

**Additional context**
Add any other context about the problem here.
