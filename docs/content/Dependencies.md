# Getting proatac running
**proatac** has a few dependencies that are listed below with relevant hyperlinks for 
installation instructions from the source. To quickly determine what may be lacking in
your system, try running **proatac** with the [default.yaml](yaml/CLmac.yaml) file
(more on that [here](#yaml)) using the `--check` flag. To do this, we'll first clone
the repository

```
git clone https://github.com/buenrostrolab/proatac.git
proatac yaml/default --check
```

If you get a message saying that the check was succesful, then you're most likely
ready to begin analyzing data. However, if you run into one or more error messages, 
you are likely missing the necessarily software. Make sure that

- [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and relevant index for analysis. 
- [java language](https://www3.ntu.edu.sg/home/ehchua/programming/howto/JDK_Howto.html)
- [macs2](https://github.com/taoliu/MACS)

We note that macs2 though also a PyPi package is only compatible with Python 2.7
whereas **proatac** is a Python 3 package. There's a good chance that macs2
is already living in your environment if you are reading this help page, which can
be tested using the following--

```
which macs2
```

and hopefully seeing a valid path. If not, one solution for macs2 install is to create
a separate python2 virtual environment using the following commands -- 

```
python2 -m venv venv2
source venv2/bin/active

pip install numpy
pip install wheel
pip install macs2
```

- [R language](https://www.r-project.org/) and package dependencies
(see [wiki/Rpackages](https://github.com/buenrostrolab/proatac/wiki/Rpackages) for more information). 

- [samtools](http://www.htslib.org/download/)
