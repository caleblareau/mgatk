<p align="left">
  <br><br><br>
  <img src="docs/content/media/logo.png" width="50%"/>
</p>

[![TeamCity (simple build status)](https://img.shields.io/teamcity/http/teamcity.jetbrains.com/s/bt345.svg)]()
[![Documentation Status](https://readthedocs.org/projects/mgatk/badge/?version=latest)](http://mgatk.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)


The **mgatk** package was developed by Caleb Lareau and Jacob Ulirsch and is maintained by
[Caleb](mailto:caleblareau@g.harvard.edu) in the
[Aryee Lab](https://aryeelab.org) and [Buenrostro Lab](https://buenrostrolab.com).
Source code is made [freely available](http://github.com/aryeelab/mgatk)
and a packaged install version is provided through [PyPi](https://pypi.python.org/pypi/mgatk/).
<br><br>

## About
We provide two interfaces for **mgatk**, including a python-based command line interface for
processing `.bam` files with mitochondrial reads and generating quality control, blacklist, 
and variant calls. These can be read into our `R` package that enables quality control, annotation, 
sample deconvolution, transcriptomic and epigenomic heritability analyses, as well as 
visualizations and diagnostics. 
<br><br>

## Installation

**Python-**
```
pip3 install mgatk
```

**R-**
```
devtools::install_github("aryeelab/mgatk", subdir="Rpkg/mgatk")
```

**or**
```
devtools::install_github("aryeelab/mgatk/Rpkg/mgatk")
```

## Install via github

**Python-**
```
pip3 install git+https://github.com/aryeelab/mgatk.git#subdirectory=Pypkg/mgatk
```

**R-**
```
devtools::install_github("aryeelab/mgatk", subdir="Rpkg/mgatk")
```

R package vignette build with [pkgdown](https://github.com/hadley/pkgdown).

## Workflow Overview

A detailed description of the workflow including a description of various parameter
settings is discussed in depth in the [**mgatk** documentation](http://mgatk.readthedocs.io).
Below is a brief overview of our workflow

<br><br>

## Installation/Documentation/FAQ/.etc

Check out the [**mgatk** documentation](http://mgatk.readthedocs.io) for detailed
installation instructions, dependency configuration, and other information.
<br><br>

### Questions/comments/feedback
are always welcomed. The easiest way for us to have correspondence (if appropriate/interesting
for the public) is through raising a [new issue](https://github.com/aryeelab/mgatk/issues/new)
on the GitHub source. For private concerns, email [Caleb](mailto:caleblareau@g.harvard.edu). 
<br><br><br>
