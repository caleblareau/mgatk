<p align="left">
  <img src="media/logo.png" width="50%"/>
</p>

[![Build Status](https://travis-ci.com/aryeelab/mgatk.svg?token=snx22Bgp4cRvvH32vAmH&branch=master)](https://travis-ci.com/aryeelab/mgatk)
[![PyPI version](https://badge.fury.io/py/mgatk.svg)](https://pypi.python.org/pypi/mgatk)
[![Snakemake](https://img.shields.io/badge/snakemake-≥4.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

The **mgatk** package was developed by Caleb Lareau.

[Source code is made freely available](http://github.com/caleblareau/mgatk)
and a packaged install version is provided through [PyPi](https://pypi.python.org/pypi/mgatk/).
<br><br>

## About
This repository houses the **mgatk** package, a python-based command line interface for
processing `.bam` files with mitochondrial reads and generating high-quality heteroplasmy 
estimation from sequencing data. This package places a special emphasis on mitochondrial
genotypes generated from single-cell genomics data, primarily `mscATAC-seq`, but is generally
applicable across other assays. 
<br><br>

## Installation

**Recommended:**
First, create a `python` virtual environment in some working directory to keep things tidy:

```
python3 -m venv venv3
source venv3/bin/activate
```

Next, install `mgatk` from [PyPi](https://pypi.org/project/mgatk/):

```
pip3 install mgatk
```

## Workflow Overview

A detailed description of the workflow including a description of various parameter
settings is discussed in depth in the [**mgatk** documentation](https://github.com/caleblareau/mgatk/wiki).
Below is a brief overview of our workflow

<br><br>

## Installation/Documentation/FAQ/.etc

Check out the [**mgatk** documentation](https://github.com/caleblareau/mgatk/wiki) for detailed
installation instructions, dependency configuration, and other information.
<br><br>

