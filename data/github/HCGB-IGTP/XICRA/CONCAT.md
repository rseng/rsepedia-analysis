# XICRA: Small RNAseq pipeline for paired-end reads.

## Table of Contents

- [Decription](#description)
  * [Installation](#installation)
  * [XICRA.stats](#xicrastats)
  * [Etymology](#etymology)
- [Supplementary information](#supplementary-information)
- [Documentation](#documentation)
- [License](#license)
- [Citation](#citation)

## Description

XICRA is a python pipeline developed in multiple separated modules that it is designed to take paired end fastq reads, trim adapters and low-quality base pairs positions, and merge reads (R1 & R2) that overlap. Using joined reads it describes all major RNA biotypes present in the samples including miRNA and isomiRs, tRNA fragments (tRFs) and piwi associated RNAs (piRNAs). 

This pipeline resulted from the observation that potential artifacts derived from sequencing errors and/or data processing could result in an overestimation of abundance and diversity of miRNA isoforms. Paired end sequencing improves isomiR calling in small RNA sequencing data. Internal variation isomiR calls are frequent artifacts in single read sequencing data. Internal sequence variant isomiRs in single read sequencing mode may be false positives. See additional detail in the original publication [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04128-1)

So far, XICRA produces a miRNA or tRNA analysis at the isomiR or tRF level using joined reads or single-end reads, multiple software at the user selection and following a standardization procedure. Results are generated for each sample analyzed and summarized for all samples in a single expression matrix. This information can be processed at the miRNA or isomiR level (single sequence) but also summarizing for each isomiR variant type. This information can be easily accessed using the accompanied R package 
[XICRA.stats](https://github.com/HCGB-IGTP/XICRA.stats). Although the pipeline is designed to take paired-end reads, it also accepts single-end reads. 

See additional details on the code [here](https://github.com/HCGB-IGTP/XICRA/tree/master/XICRA_pip). The workflow of the pipeline is described in the following image.
<img src="XICRA_pip/docs/source/images/workflow/XICRA_pipeline.png" alt="Workflow" width="950"/>

### Installation

XICRA will require python v3.6 and java (we tested in openjdk 14 2020-03-17).

The XICRA python pipeline is available in `pip` and also available using `conda`.   

XICRA depends on multiple third party software that we have listed [here](https://github.com/HCGB-IGTP/XICRA/blob/master/XICRA_pip/README.md#dependencies).

We encourage you to install XICRA and all dependencies using the `conda` environment we created and following these instructions. 

Unfortunately, a couple of executables are not available neither as a `conda` or `pip` packages. These packages are `miraligner` and `sRNAbench`. We have generated a `shell` script to retrieve and include within your `conda environment`.

To create a new conda environment, install third party software, install XICRA and missing dependencies, do as follows:

```sh
## clone repo
git clone https://github.com/HCGB-IGTP/XICRA.git

## move to folder
cd XICRA

## create conda environemt
conda env create -n XICRA -f XICRA_pip/devel/conda/requirements.txt

## activate
conda activate XICRA

## install latest python code
pip install XICRA

## install missing software
sh XICRA_pip/XICRA/config/software/installer.sh
```

To check everything is fine, try executing the `config` module:
```sh
XICRA config
```

On the other hand, if you might have already installed software and available within your path, you might only need to install it using the [XICRA pip](https://pypi.org/project/XICRA/) module. We encourage you to installed it within a python environment. See as an example the description [here](https://github.com/HCGB-IGTP/XICRA/blob/master/XICRA_pip/README.md#python-environment)
 
### xicrastats

We additionally provide a supplementary R package for parsing and plotting some XICRA results. See additional details [here](https://github.com/HCGB-IGTP/XICRA.stats).

Install it in `R` using:

```R
# Install XICRA.stats version from GitHub:
# install.packages("devtools")
devtools::install_github("HCGB-IGTP/XICRA.stats")
```


### Etymology
XICRA is the name in Catalan for a small cup of coffee or chocoate. Also, it was used as a small measure of milk, oil or wine (125 ml). See additional details here: https://ca.wikipedia.org/wiki/Xicra

_Una xicra és tassa petita, més aviat alta i estreta, emprada expressament per a prendre la xocolata desfeta o cafè. També és una unitat de mesura de volum per a líquids que es feia servir a Catalunya per a l'oli, vi, o llet._ 

## Supplementary Information
In this repository we provide supplementary information for the original paper describing the method. See additional details in folder BMC_bioinformatics_paper or [here](BMC_bioinformatics_paper/README.md)

## Documentation
For a full documentation and details visit Read the Docs site [here](https://xicra.readthedocs.io/). 

See a brief example on how to install and run XICRA [here](https://github.com/HCGB-IGTP/XICRA/tree/master/XICRA_pip#example)

## License 

MIT License

Copyright (c) 2020 HCGB-IGTP

See additional details [here](XICRA_pip/LICENSE)

Developed and maintained by Jose F. Sanchez-Herrero and Lauro Sumoy at HCGB-IGTP

http://www.germanstrias.org/technology-services/genomica-bioinformatica/

## Citation
Sanchez Herrero, J.F., Pluvinet, R., Luna de Haro, A. et al. Paired-end small RNA sequencing reveals a possible overestimation in the isomiR sequence repertoire previously reported from conventional single read data analysis. BMC Bioinformatics 22, 215 (2021). https://doi.org/10.1186/s12859-021-04128-1

# XICRA: Small RNAseq pipeline for paired-end reads

## Description

XICRA is a python pipeline developed in multiple separated modules that it is designed to take 
paired end fastq reads, trim adapters and low-quality base pairs positions, and merge reads (R1 & R2) 
that overlap. Using joined reads it describes all major RNA biotypes present in the samples including 
miRNA and isomiRs, tRNA fragments (tRFs) and piwi associated RNAs (piRNAs).

So far, XICRA produces a miRNA analysis at the isomiR level using joined reads, multiple software at the 
user selection and following a standardization procedure. Results are generated for each sample analyzed and 
summarized for all samples in a single expression matrix. This information can be processed at the miRNA or 
isomiR level (single sequence) but also summarizing for each isomiR variant type. This information can be 
easily accessed using the accompanied R package XICRA.stats. Although the pipeline is designed to take 
paired-end reads, it also accepts single-end reads.

## Installation

XICRA will require python v3.7 and java (we tested in openjdk 14 2020-03-17).

The XICRA python pipeline is available in `pip` and also available using `conda`.

XICRA depends on multiple third party software that we have listed below.

### Dependencies 

Python XICRA module will install itself along some python modules dependencies (pandas, multiqc, pybedtools, biopython etc.). 

But additionally, XICRA depends on third party software that we listed in the following [table](https://github.com/HCGB-IGTP/XICRA/blob/master/XICRA_pip/XICRA/config/software/soft_dependencies.csv).

### Conda environment

We encourage you to install XICRA and all dependencies using the `conda` environment we created and following these instructions. 

To create a new conda environment, install third party software, install XICRA and missing dependencies, do as follows: 

1) Get requirements file from XICRA git repo

```sh
wget https://raw.githubusercontent.com/HCGB-IGTP/XICRA/master/XICRA_pip/devel/conda/environment.yml
```

2) Create environment and install required packages using conda: 

```sh
conda env create -n XICRA -f environment.yml
```

3) Activate environment and install XICRA
```sh
## activate
conda activate XICRA

## install latest python code
pip install XICRA
```

4) Install missing software:  Unfortunately, a couple of executables are not available neither as a `conda` or `pip` packages. These packages are `miraligner` and `sRNAbench`. We have generated a `shell` script to retrieve and include within your `conda environment`.

```sh
## install missing software
sh XICRA_pip/XICRA/config/software/installer.sh
```

To check everything is fine, try executing the `config` module:
```sh
XICRA config
```

### Python environment

If you are not using a `conda` environment as you might have previously installed all dependencies, we encourage you to create a python environment containing all python modules required for XICRA. See as an example this code:

```sh
## create enviroment
python3 -m venv XICRA_env

## activate it
source XICRA_env/bin/activate

## install XICRA and dependencies
pip install XICRA

## execute XICRA
XICRA -h
```

## Documentation

See a full documentation, user guide and manual in [here](https://readthedocs.org/)

## Example
Here we include a brief example on how to use XICRA.

First, we create a python environment and will install XICRA and dependencies. See example details shown before.
Then, we can test XICRA by using an example of 100 miRNA simulated and provideded within the repository as an example of simulation.

```sh
## run XICRA example
ln -s ~/BMC_bioinformatics_paper/simulation/example/reads/

## prepare reads
XICRA prep --input reads/ --output_folder test_XICRA

## join reads
XICRA join --input test_XICRA --noTrim

## create miRNA analysis
XICRA miRNA --input test_XICRA --software miraligner sRNAbench

## explore results
ls test_XICRA/report/
```

## Documentation
For a full documentation and details visit Read the Docs site [here](https://xicra.readthedocs.io/). 

See a brief example on how to install and run XICRA [here](https://github.com/HCGB-IGTP/XICRA/tree/master/XICRA_pip#example)

## License 

MIT License

Copyright (c) 2020 HCGB-IGTP

See additional details [here](XICRA_pip/LICENSE)

Developed and maintained by Jose F. Sanchez-Herrero and Lauro Sumoy at HCGB-IGTP

http://www.germanstrias.org/technology-services/genomica-bioinformatica/

## Citation
Sanchez Herrero, J.F., Pluvinet, R., Luna de Haro, A. et al. Paired-end small RNA sequencing reveals a possible overestimation in the isomiR sequence repertoire previously reported from conventional single read data analysis. BMC Bioinformatics 22, 215 (2021). https://doi.org/10.1186/s12859-021-04128-1

## Authors
Antonio Luna de Haro (v0.1)
Jose F Sanchez-Herrero (v1.0)	
# XICRA documentation

This is the top level build directory for the XICRA documentation.
  
Documentation is written using [Sphinx](http://www.sphinx-doc.org/en/master/), a python documentation system built 
using [reStructuredText](http://docutils.sourceforge.net/rst.html). 

Find the documentation hosted in readthedocs.org: 


Building the documentation
--------------------------

To build the docs see deatils in file [documentation for developers](https://readthedocs.org/) [XICRA/docs/source/devel/documentation.rst]. <TODO>


Organization
------------

See details of the organization of this docs directory in docs/devel/documentation.rst file. <TODO># Instructions for developers

## Create conda env for developers
We have included files in the folder `devel/conda` for the build and configuration of the conda package.

To create a new conda environment, install third party software, install XICRA and missing dependencies for XICRA development, do as follows:

```sh
## clone repo
git clone https://github.com/HCGB-IGTP/XICRA.git

## move to folder XICRA_pip
cd XICRA/XICRA_pip

## create conda environemt
conda env create -n XICRA -f ./devel/conda/requirements.txt

## activate
conda activate XICRA

## install latest python code
pip install -r ./devel/pypi/requirements.txt
pip install -e .

## install missing software
sh ./XICRA/config/software/installer.sh
```

To check everything is fine, try executing the `config` module:
```sh
XICRA config
```


## Instruction for creating releases

One on hand, we can create a new `pip` package, also, we would create a `conda` release. Ideally, all would be concordant with Github code releases.

### Create pip package release

XICRA python code is distributed in pip (https://pypi.org/project/XICRA/). 

Here, we have generated a few scripts in `devel/pypi` including all commands for the generation of the pip package. First, we would need to clean previous builds, then create a new distrubution to finally upload it to pip website. 

Basically, in the main git folder, type:

#### clean distribution packages
```
sh devel/pypi/clean_devel.sh
```

This code basically removes old builds in folders: `build`, `dist` and `XICRA.egg-info`:


#### create distribution files
Now, we would create a new distribution, type either:

```
sh devel/pypi/create_distro.sh
```

or

```
python setup.py sdist bdist_wheel
```

Folders `build`, `dist` and `XICRA.egg-info` will be created.

#### Upload to pip

Basically, XICRA pip package is hosted under my user (I might change it in the future to a admin user).

To upload the builds generated, we first need to create a file named as `.pypirc` in the main git directory. Include the code:

```
$ nano .pypirc
[distutils] 
index-servers=pypi
[pypi] 
repository = https://upload.pypi.org/legacy/ 
username =jfsanchezherrero
```

Then, upload the build using the followind command:
```
sh devel/pypi/upload_pypi.sh
```

or

```
python -m twine upload dist/*
```

NOTE: You will need `jfsanchezherrero` password to successfully upload the code.

### references
See additional details on this topic in following linkgs:

https://dzone.com/articles/executable-package-pip-install
https://packaging.python.org/tutorials/packaging-projects/
# Supplementary Information
In this folder we provide additional information regarding the data included in the original paper.

The in-house dataset generated ([GSE155370](analysis_GSE155370/)) and the retrieved dataset from NCBI SRA ([GSE114923](analysis_GSE114923/))

We include the simulation command details and recipe we followed for the XICRA simulations (See details [here](simulation/README.md)).

We also provided a folder [code](code/) containin several scripts to replicate the plots included in the original paper.
# Recipe for XICRA simulations

To illustrate the potential of paired-end reads at the isomiR level analysis we have generated computer simulations to test the impact of technical errors from single end or paired-end reads. We followed the guidelines previously described for isomiR computer simulations by [Amsel et al. 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1772-z)

We created biological variation and technical variation using multiple high throughput sequencing profiles and evaluated the performance of the simulation using sensitivity and precision of the isomiRs detected under several circumstances. 

## Biological variation

We created artificial miRNA isoforms from Homo sapiens mature and hairpin sequences in miRBase using the bioinformatic scripts previously described at [isomiR-Benchmark](https://github.com/DanielAmsel/isomiR-Benchmark). We additionally added the canonical fasta sequence of each miRNA to the variant dataset generated. 

```sh
mkdir mirBase
cd mirBase

## download latest
wget -nd ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
wget -nd ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

# extract gunzip files
gunzip *

## get hsapiens mature and hairpin
grep 'sapiens' mature.fa | awk '{print $1}' | sed 's/>//' > mature_hsapiens.ids.txt
grep 'sapiens' hairpin.fa | awk '{print $1}' | sed 's/>//' > hairpin_hsapiens.ids.txt

## get hsapiens fasta files
seqtk subseq hairpin.fa hairpin_hsapiens.ids.txt > hairpin_hsapiens.fa
seqtk subseq mature.fa mature_hsapiens.ids.txt > mature_hsapiens.fa

## Simulate isomiRs and reads
cd ../
mkdir isomiR_simulations
cd isomiR_simulations

## simulate isomiR
perl ~/isomiR-Benchmark/create_isomiRs.pl ../mirBase/mature_hsapiens.fa ../mirBase/hairpin_hsapiens.fa

## Concat all sequences & change character
cd ..
cat isomiR_simulation/*fa | sed 's/|/::/' > all_seqs.fa

# add canonical & conver U->T & rename/add character
python rename_canonical.py mature_hsapiens.fa mature_hsapiens_renamed.fa
cat mature_hsapiens_renamed.fa >> all_seqs.fa


```

From the variant frequency table generated we select 100 miRNAs and discarded some random variants and generated random distribution frequencies of the variant types generated. To simplify the interpretation and evaluation, for each variant type, we selected a single isoform. For each variant type we included the corresponding frequencies generated to a total amount of 100 sequences for each mature miRNA. 

```sh
## get frequency
python get_freq.py all_seqs.fa fasta all_seqs.freqs.csv
```

## Technical simulation

For each biological dataset generated, we used ART (version Mount Rainier 2016–06-05 ) with Illumina HiSeq2500 and MiSeq-v1 sequencing system profiles using paired-end mode to simulate next generation sequencing (NGS) reads. We grouped all isoforms according to length and generated NGS simulation for each length subset to finally merged them all in a single file for each read accordingly. We used a 10x sequencing coverage for each input fasta sequence. 

As previously noted (Amsel, 2017), due to the nature of the ART simulation we had to parse and omit about half the total reads generated as they were reverse complemented. We only discarded reverse complemented R1 reads and its R2 counterpart accordingly. Due to the implementation based on frequencies that we did in the biological variation procedure, we made sure when applying the same coverage for each sequence, and discarding reverse complement, that the biological variation frequencies generated would be maintained in the NGS simulation. The observed range would vary from 5-500 counts for each single variant type simulated. 

We evaluated the performance of using PE reads for miRNA isomiR analysis using the NGS simulation datasets and the pipeline XICRA. For each dataset, we used the miRNA module using paired-end mode and single end mode for the R1 and R2 reads. For the paired-end mode we initially joined reads using two different join percentage difference cutoff (fastq-join parameter) to test the effect of using 100% perfect R1 and R2 reads or allowing the default difference (8%) along the minimum default overlap length cutoff (6 bp). For the single-end mode, we used the total R1 reads simulated or the total R2, reversed complemented using seqtk software, respectively. We generated a miRNA analysis at the isomiR level using the three different software available within XICRA: miraligner, sRNAbench and optimir.

All these steps mention here are implemented in simulation_sender.py

```sh
## simulate NGS reads using art
python simulation_sender.py --fasta all_seqs.fa --folder example --reads PE --seqSys HS25 --fcov 10 -t 2 --freqs all_seqs.freqs.csv -n 10 -r 10 --art_bin ./art_illumina -m 50 -s 5 --database db_mirbase
```

## Performance evaluation

Using the biological variation frequencies generated as true positives for each dataset, we evaluated the amount of isomiRs detected for each software and each type of read. For paired-end reads, we also used a different percentage difference cutoff. 

For each isomiR, the detected counts were classified as: True positives (TP) when observed counts matched the expected counts; false positives (FP) when observed counts exceeded the expected counts and were wrongly assigned; false negatives (FN) when observed counts did not get to the minimum expected counts. We calculated the sensitivity or recall as TP/(TP+FN) and the precision or specificity as TP/(TP+FP). We also reported True Negatives (TN) when expected counts were not observed and new generation isomiRs when new variants or miRNA appeared and were not expected. 
## NCBI Project
https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA473134

## PUBMED id
https://pubmed.ncbi.nlm.nih.gov/30006517/
doi: 10.1038/s41598-018-28854-4

## GEO
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114923
.. _contents:

********
Contents
********

Documentation
=============

.. only:: html

    :Version: |version|
    :Date: |today|

.. toctree::
   :maxdepth: 2

   user_guide/user_guide_index.rst
   devel/devel_index.rst
   info/info_index.rst
   glossary/glossary_index.rst
   bib/zbib_index.rst
   

Indices and tables
==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`:orphan:

.. title:: XICRA documentation 

XICRA Documentation
*******************

.. only:: html 

    :Version: |version|
    :Date: |today|
    
Introduction
============

``XICRA`` is a pipeline that integrates multiple bioinformatic tools for the analysis of
paired-end or single end reads from small RNA-seq data. It describes all major RNA biotypes 
present in the samples including miRNA and isomiRs, tRNA fragments (tRFs) and piwi associated 
RNAs (piRNAs). Results are generated for each sample and summarized for all samples in a single 
expression matrix.

The pipeline is written in Python_ with a modular architecture and based on open-source software 
and databases engines. The design of this bioinformatic tool allows miRNA analysis at different 
levels.

Multiple tasks are performed by several modules including:
 - preparation of raw data
 - quality analysis and trimming of the adapters
 - merge reads (R1 & R2) that overlap
 - map reads to reference genome and annotation 
 - quantification of RNA types
 - miRNA and isomiR quantification (including the variant type)
 - preparation of results for integrative visualization
 
The tool uses the miRTop database and its notation to quantify and report the miRNAs present in each sample. 
With the resulting matrix for each sample the analysis can be performed at the miRNA, isomiR, or variant type 
level.

The ``XICRA`` documentation includes: 

.. #
.. TODO
.. #

- A :doc:`User's Guide <user_guide/user_guide_index>` to get started.
- An example :doc:`Tutorial <tutorial/tutorial_index>`.
- A :doc:`Quick Start Guide <user_guide/quickStart>`.
- A list of :doc:`user_guide/faqs`
- Some :doc:`developer Guidelines <devel/devel_index>` to contribute to the project.
- Additional :doc:`information <info/info_index>` of the ``XICRA`` project.
- A list of :doc:`Glossary <glossary/glossary_index>` terms.
- A list of :doc:`Bibliography <bib/zbib_index>`


.. ##################
.. _pipeline-scheme:
.. ##################

Pipeline Scheme
===============

Here we show the scheme of the ``XICRA`` bioinformatic tool. It is divided in six main actions:
   
   #. **Preparation of the input data**: preparation of the fastq files from a sequencing run.
   
   #. **Quality analysis**: with quality check programs attending the input provided.
   
   #. **Adapters trimming**: for each read the adapter sequences are filtered out.
   
   #. **Read joining**: joins sequencing reads (paired-end).
   
   #. **Mapping reads and feature counts**: generates a RNA biotype analysis, quantifying each RNA type present in the samples.
   
   #. **miRNA analysis**: generates a miRNA analysis, with isomiR quantification and variant type information.
   
This information can be easily accessed using the accompanied R package XICRA.stats_. Although the pipeline is designed to take
paired-end reads, it also accepts single-end reads.


.. image:: images/workflow/XICRA_pipeline.png
   :width: 1500pt


Table of Contents
=================
.. toctree::
   :maxdepth: 4

   contents.rst
 
.. include:: links.inc
.. _glossary:

********
Glossary
********
.. #
.. TODO: add termns
.. #

.. glossary::

   miRNA
      Class of small non-coding RNAs (ncRNAs), have an average length of 21–23 nucleotides.
      They are endogenous regulatory molecules that modulate gene expression post-transcriptionally 
      by inducing target mRNA silencing and decay. Additional roles beyond negative modulation 
      of mRNA function have also been proposed.
      
   isomiR
      miRNAs can frequently appear in the form of multiple sequence variants or isoforms (termed isomiRs).
   
   Reads
      An inferred sequence of base pairs (or base pair probabilities) corresponding to all or part of a single DNA fragment. 
   
   mirBase
      A searchable database of published miRNA sequences and annotation.
      
   mirTop
      Command line tool to annotate with a standard naming miRNAs e isomiRs.
   
   Pypi
      Python package index (Pypi, https://pypi.org/) is a repository of software for the 
      Python programming language.
   
   reStructuredText (ReST)
      reStructuredText format (http://docutils.sourceforge.net/rst.html). An easy-to-read, what-you-see-is-what-you-get 
      plaintext markup syntax and parser system.
   
   Sphinx
      Python module for documentation. See details in http://www.sphinx-doc.org/en/master/
   

   
      
     ######################
Additional information
######################

.. only:: html

    :Version: |version|
    :Date: |today|
    

.. contents::

.. # 
.. TODO: complete this page
..

.. _background:

Project Status
==============

The project is available on `XICRA github`_. You can 
`report issues`_ or contribute to the project by `forking the project`_
and `creating pull requests`_. 
See additional details on XICRA :ref:`developer guidelines<developer>` and how to :ref:`work with Git<git-guidelines>`.


Background
==========

MicroRNAs (miRNAs), a class of small non-coding RNAs (ncRNAs), have an average length of 21–23 nucleotides (nt). 
They have been widely studied as endogenous regulatory molecules that modulate gene expression post-transcriptionally 
by inducing target mRNA silencing and decay. Additional roles beyond negative modulation of mRNA function have 
also been proposed. 

Most miRNA expression studies based on next generation sequencing (NGS) performed to this date have summarized all 
the reads mapping to a specific miRNA locus or miRNA sequence with or without mismatches and assign it to a single 
miRNA entity (a miRBase reference database entry). However, this type of analysis neglects the fact that not all reads 
are identical to the mature reference sequence in miRBase. Small RNA sequencing NGS methodology has revealed that 
miRNAs can frequently appear in the form of multiple sequence variants or isoforms (termed isomiRs).

``XICRA`` is a pipeline to analyze small RNA sequencing (small RNA-seq) data, which allows detecting and quantifying 
isomiR level variation in miRNA.


.. _license:

License
=======

This is brief description of the license terms for ``XICRA``. 


.. _credits:

Credits
=======

Give credit to who deserves


.. _citation:

Citation
========

Please cite the ``XICRA`` project, code or documentation using the following links:

.. github_stats:

Github statistics
=================

These are Github statistics generated for each release of the code.

.. include:: ./github_stats.txt
   :literal:

.. _whats_new:

What's new?
===========

See below information on differences among each release of the code generated.

.. include:: ./whats_new.txt
   :literal:
  
 
.. include:: ../links.inc
    .. ########################
.. _format-fastq-files:
.. ########################

Sample name details
===================

``XICRA`` accepts a wide range of files names. 

Format for fastq files can be: name.fastq.gz, name_1.fastq.gz, name_R2.fastq.gz, name_L001_R1.fastq.gz, etc.

Here we provide some guidelines of the format.


Length limitation
-----------------

There is a limitation for the sample ID ('name') of 10 characters.
** XICRA provides an option to rename samples if necessary: module prep option --rename **

Single end files
----------------

It is possible to provide NGS single-end files although some steps of the process could be accomplished using single-end files.
** Use option --single-end in the different XICRA modules **

name.fastq.gz

name.fastq

name.fq

Paired-end files
----------------

Paired-end files are full supported. The format for these files are:

name_1.fastq.gz or name_R1.fastq.gz

name_2.fastq.gz or name_R2.fastq.gz


Lane information:
-----------------

In some cases, files might contain lane information (*L00x* and/or *00x*).
``XICRA`` supports these names as long as follow these examples:

name_L00x_R1.fastq.gz   name_L00x_R2.fastq.gz

name_L00x_1.fastq.gz name_L00x_2.fastq.gz

name_L00x_R1_00x.fastq.gz  name_L00x_R2_00x.fastq.gz

** If you need to merge fastq files from different lanes, use option --merge within module prep **


Extensions:
-----------

name_L00x_R2.fastq   

name_L00x_R2.fq

name_L00x_R2.fastq.gz   

name_L00x_R2.fq.gz


********
Tutorial
********
.. #
.. TODO
.. #

.. only:: html

    :Version: |version|
    :Date: |today|

.. toctree::
    :maxdepth: 2.. ############
.. index-bib:
.. ############
 
Bibliography
^^^^^^^^^^^^

.. bibliography:: ./references.bib
	:cited:
	:style: unsrt.. title:: API guide

API Overview
------------
.. toctree::
   :hidden:
   
   XICRA.rst
   modules/modules_index.rst
   
.. 
.. TODO
..   scripts/scripts_index.rst
..   other_tools.rst
..   data.rst
..   third_party.rst
..
   

.. only:: html

    :Version: |version|
    :Date: |today|

Here we include and API_ documentation for ``XICRA`` in order to provide third parties to use the 
functionality of ``XICRA`` application.

This pipeline is composed of multiple modules and scripts that are separated in a main script:

- ``XICRA.py``

This main script integrates and connects all available modules and analysis arranged in 
several directories:

..
.. TODO: review this, other tools? data? also in contents
..

- ``XICRA/modules``
- ``XICRA/scripts``
- ``XICRA/data``
- ``XICRA/other_tools``
- ``XICRA/config``


Contents
^^^^^^^^

XX

.. toctree::
   :maxdepth: 2

   XICRA.rst
   modules/modules_index.rst
   
.. 
.. TODO
..   scripts/scripts_index.rst
..   config/config_index.rst
..   other_tools.rst
..   data.rst
..
   
 
.. ############ 
.. include:: ../links.inc
 

 prep
====

This module prepares fastq files from a sequencing run for later usage.

For each sample, a folder called as the sample will be generated conatining
the links to the raw data (or copied raw data files) that belong to the sample. 


Workflow
--------
TODO: build image


Functions
---------

.. automodule:: XICRA.modules.trimm
   :members:
 


.. include:: ../../links.incmiRNA
=====

..
.. TODO: link trimm module, link miRNA from user guide
..

This module analyses the joined or single reads, that have been 
previously trimmed of different samples. To do so, it uses three 
possible softwares, which can be run at the same time if desired. 

For each sample, an expression matrix using the miRTop nomenclature
is generated, containing the information of the counts at miRNA or
isomiR level; also describing the variant type of each isomiR.

In addition, an expression matrix comparing all the samples is created.


Workflow
--------
TODO: build image 


Functions
---------

.. automodule:: XICRA.modules.miRNA
   :members:
   
Other functions
---------------
This module calls the function:

.. currentmodule:: XICRA.scripts.generate_DE
.. autofunction:: generate_DE

.. include:: ../../links.incQC
===

This module calls fastqc_, a quality check program, to analyze the 
quality of each sample. 

For each sample, a folder called "fastqc" is generated with the 
quality analysis, reported in an html file.

In addition, by default, a multiQC_ report is generated for all the samples.


Workflow
--------
TODO: build image

Functions
---------

.. automodule:: XICRA.modules.qc
   :members:
   
Other modules
-------------
The ``QC`` module also uses the ``multiQC_report`` module:

.. automodule:: XICRA.scripts.multiQC_report
   :members:

.. include:: ../../links.inc.. _modules:

Modules
=======

.. only:: html

    :Version: |version|
    :Date: |today|

Developer guidelines for ``XICRA`` modules.

.. toctree::
   :maxdepth: 1

   prep.rst
   QC.rst
   trimm.rst
   join.rst
   biotype.rst
   miRNA.rst
biotype
=======
.. 
.. TODO: XICRA.scripts.RNAbiotype functions description
..


This module works in two steps. First, it maps the unjoined reads
to a given reference genome. Secondly, it generates a summary of the RNA
types found in each sample. 

Mapping
-------
The mapping step is performed with STAR_ software and the code is in a separate
script, ``mapReads.py``. To introduce the reference genome the are two options:

1. Introduce the path to the fasta sequence of the human genome, option -\ -fasta, and
   indicate the path to a reference genome annotation in GTF format,  -\ -annotation.
2. Introduce the STAR indexed genome file, option -\ -genomeDir. This directory, 
   called 'STAR_index' is created after the execution of STAR, thus, if it is 
   already generated it can be reused, check the STAR_ documentation.

The -\ -limitRAM parameter is also important. It indicates the maximum RAM (in bytes) that 
the computation will use to prevent the computer collapse. Note that, the STAR software requires high 
values of RAM  in order to do the mapping. Thus, a 'regular' laptop may not be able to perform 
this computation. We are currently working in the implementation of alternatives to this software 
that can be used with any computer. 

As a result, the mapping step generates a 'map' folder with a BAM file and a report from MultiQC
for each sample.

Feature counts
--------------
The second step performs the RNA biotype analysis  using the featureCounts_ software, 
code located in the script ``RNAbiotype.py``. After the computation for each sample, a final
MultiQC report is generated to compare the composition of each sample and, eventually, detect outliyers.


Workflow
--------
TODO: build image 


Functions
---------

.. automodule:: XICRA.modules.biotype
   :members:
   
Other functions
---------------
This module calls other scripts:

.. automodule:: XICRA.scripts.mapReads
   :members:

.. automodule:: XICRA.scripts.RNAbiotype
   :members:

.. include:: ../../links.inctrimm
=====

This module cuts the introduced adapter sequences of the input reads.

For each sample, a folder called "trimm" is generated with the trimmed
reads, called as the original ones +"_trimm".

In addition, a multiQC report is generated.


Workflow
--------
TODO: build image

Functions
---------

.. automodule:: XICRA.modules.trimm
   :members:
   
Other modules
-------------
The ``trimm`` module also uses the ``multiQC_report`` module:

.. automodule:: XICRA.scripts.multiQC_report
   :members:

.. include:: ../../links.inc.. _XICRA-main:

XICRA
""""""""""""""
.. automodule:: XICRA
    :members:

.. _developer:

Developer Guidelines
====================

.. only:: html

    :Version: |version|
    :Date: |today|

Guidelines for developing the ``XICRA`` project.

.. toctree::
   :maxdepth: 1
    
   ../api/api_index.rst
   developer_guidelines.rst
   ../data/data_index.rst
   documentation/documentation.rst
   working_with_git.rst.. _extern_progs:

extern_progs
============
.. automodule:: XICRA.config.extern_progs
    :members:
.. _config:

Config
======

.. only:: html

    :Version: |version|
    :Date: |today|

Developer guidelines for ``XICRA`` configuration scripts.

.. toctree::
   :maxdepth: 1
    
   extern_progs.rst
   set_config.rst
.. _set-configuration:

set_config
==========
.. automodule:: XICRA.config.set_config
    :members:
.. _scripts:

Scripts
=======

.. only:: html

    :Version: |version|
    :Date: |today|

Developer guidelines for ``XICRA`` scripts.

.. toctree::
   :maxdepth: 1

   RNAbiotype.rst
   fastqc_caller.rst
   functions.rst
   generate_DE.rst
   multiQC_report.rst
   reads2tabular.rst
   sampleParser.rst

.. _reads2tabular:

reads2tabular
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/reads2tabular.svg
    :align: center

.. automodule:: XICRA.scripts.reads2tabular
    :members:
    :undoc-members:


.. _sampleParser:

sampleParser
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/sampleParser.svg
    :align: center

.. automodule:: XICRA.scripts.sampleParser
    :members:
    :undoc-members:

.. _RNAbiotype:

RNAbiotype
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/RNAbiotype.svg
    :align: center

.. automodule:: XICRA.scripts.RNAbiotype
    :members:
    :undoc-members:

.. _functions:

functions
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/functions.svg
    :align: center

.. automodule:: XICRA.scripts.functions
    :members:
    :undoc-members:

.. _generate_DE:

generate_DE
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/generate_DE.svg
    :align: center

.. automodule:: XICRA.scripts.generate_DE
    :members:
    :undoc-members:

.. _fastqc_caller:

fastqc_caller
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/fastqc_caller.svg
    :align: center

.. automodule:: XICRA.scripts.fastqc_caller
    :members:
    :undoc-members:


.. _multiQC_report:

multiQC_report
==========================================
This script contains several functions. Here we show a graph representation of the different functions and relationships among them:

.. image:: ../../images/python_graph/multiQC_report.svg
    :align: center

.. automodule:: XICRA.scripts.multiQC_report
    :members:
    :undoc-members:

.. _users-guide-index:

############
User's Guide
############

.. only:: html

    :Version: |version|
    :Date: |today|

This User Guide provides information on ``XICRA`` usage and documentation. It is divided in several main chapters:

.. toctree::
    :maxdepth: 1

    installation/installing.rst
    quickStart.rst
    modules/user_guide_modules_index.rst
    ../tutorial/tutorial_index.rst
    report/report_index.rst
    software/software_index.rst
    info/info_index.rst
    faqs.rst.. ################################################
.. _faqs:

.. #
.. TODO: complete this page:
.. #
Frequently Asked Questions (FAQs)
*********************************

.. contents::

This is a collection of FAQs for ``XICRA`` tutorial and interpretation of results. 

.. seealso::

    Please read further information on each section following the links provided or
    the main :ref:`User Guide for XICRA<users-guide-index>`.

Installation
============

- What is required for the installation?

In order to correctly install ``XICRA``, it is necessary to have python3 and 
python development and virtual environment libraries installed. See additional details 
in section :ref:`System requirements<system-requirements>`.

- The common errors during the installation are:



- How do I know the conda environment is activated?

Once you execute the activation of the environment, either using the script ``conda activate XICRA`` or 
by executing ``source activate XICRA``, you should see a tag in the command-line, ``(XICRA)``, prompt as shown in 
the following image.

.. image:: ../images/FAQs_conda.jpeg
   :align: center

Once you deactivate the environment this tag should disappear. 

Read additional details in Conda_ official documentation website.

.. #
.. TODO: add other blocs as:
.. Quality control
.. ===============
.. #


.. ######################
.. include:: ../links.inc
.. ################################################
.. _quickStart:

Quick start
***********

This is quick start guide for impatient people. 

Installation
============
Build and activate a ``conda`` environment and get ``XICRA`` repository:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git
   
   ## create conda environemt
   conda env create -n XICRA -f XICRA_pip/devel/conda/environment.yml
   
   ## activate the environment 
   conda activate XICRA
   
   ## install latest python code
   pip install XICRA
   
   ## install missing software
   sh XICRA_pip/XICRA/config/software/installer.sh
   
Sometimes, the activation of the environment is done with ``source activate XICRA`` 
instead of ``conda activate XICRA``. 
Check everything is fine by executing the ``config`` module:

.. code-block:: sh

   XICRA config
   

Execute XICRA
=============
The functionallity of ``XICRA`` is divided in modules. Each of them in charge of 
different parts of the analysis. Here we show how to run each of the
modules:

   #. **Prepation of the data**: ``prep``, ``QC``, ``trimm``, ``join``.
   
      - ``XICRA prep --input input_folder --output output_folder``: preparation of the data.
      - ``XICRA QC --input input_folder``: quality analysis of the samples.
      - ``XICRA trimm --input input_folder --adapters_a XXXXX --adapters_A YYYYYYYY``: trimming the adapters. 
         * ``--adapters_a XXXXX``: sequence of the 3' adapter of R1
         * ``--adapters_A  YYYYYYYY``: sequence of the 3' adapters of R2.
      - ``XICRA join --input input_folder``: join paired end reads.

   #. **Quantification of RNA types**: ``biotype``
   
      - ``XICRA biotype --input input_folder --threads X --fasta file.fa --annotation file.gtf --limitRAM YYYY`` 
         * ``--threads X``: number of cores enabled to run the analysis.
         * ``--fasta file.fa``: fasta file with the reference genome to map reads.
         * ``--annotation file.gtf``: reference genome annotation in GTF format. 
         * ``--limitRAM YYYY``: given RAM in bytes to run the analysis. **Note that:** this process consumes high RAM values.

   #. **miRNA analysis**: ``miRNA``
   
      - ``XICRA miRNA --input input_folder --threads X --software analyisis_tool``
         * ``--threads X``: number of cores enabled to run the analysis.
         * ``--software analyisis_tool``: three options available, ``miraligner``, ``optimir`` and ``sRNAbench``.    

Take into account that these are basic examples, each of the modules can be run with other 
different parameters. To see all the possible parameters of each module run ``-h``. For 
example, for the ``trimm`` module:

.. code-block:: sh

   XICRA trimm -h
   

Test example
------------

Here we include a brief example on how to use ``XICRA``.

First, we create a ``conda`` environment and install ``XICRA`` and its dependencies. 
See example details shown before. Then, we can test ``XICRA`` by using an example 
of 100 miRNA simulated and provided within the repository as an example of simulation.

.. code-block:: sh

   ## run XICRA example
   ln -s ~/BMC_bioinformatics_paper/simulation/example/reads/

   ## prepare reads
   XICRA prep --input reads/ --output_folder test_XICRA

   ## join reads
   XICRA join --input test_XICRA --noTrim

   ## create miRNA analysis
   XICRA miRNA --input test_XICRA --software miraligner optimir

   ## explore results
   ls test_XICRA/report/
   
As a result, we will end up with a folder for each of the run analysis for every sample. 
Thus, in the ``miRNA`` folder, we will obtain the miRNA results for our
samples, performed with two different softwares ``miraligner`` and ``optimir``. 


User data
---------

In the presented example, nor ``QC``, neither ``trimm`` steps were necessary (that is why 
``--noTrim`` was added in the ``join`` command). However, with real data, running ``QC`` 
is highly recommended to check the quality of the samples (and filter outliers if necessary).
Running the ``trimm`` command to eliminate the adapters will be necesary for NGS data. 

On the other hand, the ``biotype`` was skipped as well. This step is only informative: it maps 
and annotates the reads and quantifies the RNA types present in each of the samples. 
**Note that:** this step is very time and RAM consuming.

The ``miRNA`` analysis can be performed whith any of the three available softwares ``miraligner``,
``optimir`` and ``sRNAbench`` (they can be combined as seen in the example). Unfortunately, the
donwloading of ``sRNAbench`` is no longer available, thus, only users with the software
already installed will be able to run the miRNA analysis with it. 

Finally, ``XICRA`` is also able to analyse single end reads, in this case ``--single_end`` should
be added in the commands (and no ``join`` step would be necessary). 

.. ###########
.. _deactivate-env:
.. ###########

Deactivate environment
----------------------

After finished the execution of any ``XICRA`` module or script, it is convenient 
to deactivate the environment. You can just close the terminal but it would be more appopiate
to conveniently deactivate the environment first.

To do so, execute one of the following commands:

.. code-block:: sh
  
   conda deactivate XICRA 
   
.. code-block:: sh
  
   source deactivate XICRA 
.. ########################
.. _userguide-info:
.. ########################

Input File Format Information
*****************************
      
.. contents::

   
.. include:: format_fastq_files.rst.. ########################
.. _format-fastq-files:
.. ########################

Sample name details
===================

``XICRA`` accepts a wide range of file names. The format for fastq files can be:

- name.fastq.gz
- name_1.fastq.gz, with '1' or '2' to specify the read
- name_R2.fastq.gz, 'R1' or 'R2' to specify the read
- name_L001_R1.fastq.gz, with the lane information as 'L00X' or '00X' after the name
- name_L001_R1_001.fastq.gz, with '00X' at the end of the file name. This naming is used 
  when the fastq of a sample had been cut in different files.
- name_L001_XYZ_R1_001.fastq.gz, there can be extra info for each file, XYZ.

The input file names should be structured considering the following aspects:

Length limitation
-----------------

There is a limitation for the sample ID ('name') of 25 characters.

``XICRA`` provides an option to rename samples if necessary in the module ``prep``, option **-**\ **-rename**.

Extensions
----------
The suported extensions are:

- name_L00x_R2.fastq   
- name_L00x_R2.fq
- name_L00x_R2.fastq.gz   
- name_L00x_R2.fq.gz

Single-end files
----------------

It is possible to provide NGS single-end files although some steps of the process could not be accomplished 
using single-end files. 

- name.fastq.gz
- name.fastq
- name.fq

Use option **-**\ **-single-end in the different XICRA modules.**

Paired-end files
----------------

Paired-end files are full supported. The format for these files are:

- name_1.fastq.gz or name_R1.fastq.gz
- name_2.fastq.gz or name_R2.fastq.gz

No parameter is needed in to specify this kind of files.

Lane information
----------------

Files might contain lane information (*L00x* and/or *00x*).
``XICRA`` supports these names as long as follow these examples:

- name_L00x_R1.fastq.gz, name_L00x_R2.fastq.gz
- name_L00x_1.fastq.gz, name_L00x_2.fastq.gz


Name extensions
---------------

It can also be the case that the reads of a sample are divided in different files. In those cases, the files
should contain a name final extension: 

- name1_L001_R1_001.fastq.gz, name1_L001_R2_001.fastq.gz
- name1_L001_R1_002.fastq.gz, name1_L001_R2_002.fastq.gz
- name1_L002_R1_001.fastq.gz, name1_L002_R2_001.fastq.gz
- name1_L002_R1_002.fastq.gz, name1_L002_R2_002.fastq.gz

Extra information
-----------------
In some cases, files might contain other extra information. In the following example, XYZ is the extra information:

- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz
- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz


Sample identification
=====================

``XICRA`` will store the names of all the input files. After that, it will identify the samples. 
It can be the case that more than one file belong to the same sample. In order to pass this information to 
``XICRA``, a combination of the following parameters may be needed depending on the characteristics of the
input file names:

Option: -\ -include_lane
------------------------

If you want to include lane tags (*L00X*, *00X*) into each  each sample name (differentiate samples considering the lane):
Use option **-**\ **-include_lane within each module** and the lane tag will also be used to identify samples.

However, if you want to consider as a single sample the different lanes, you need to merge the corresponding fastq files:
Use option **-**\ **-merge_Reads** within module ``prep``.

As an example, considering the input files:

- name1_L001_R1.fastq.gz, name1_L001_R2.fastq.gz
- name1_L002_R1.fastq.gz, name1_L002_R2.fastq.gz
- name1_L003_R1.fastq.gz, name1_L003_R2.fastq.gz
- name1_L004_R1.fastq.gz, name1_L004_R2.fastq.gz

   #. By adding the option ``--include_lane`` in **all modules**, ``XICRA`` will identify four samples:
   
      - Sample 1: name1_L001_R1, name1_L001_R2
      - Sample 2: name1_L002_R1, name1_L002_R2
      - Sample 3: name1_L003_R1, name1_L003_R2
      - Sample 4: name1_L004_R1, name1_L004_R2
      
      Remember to use option **-**\ **-include_lane within each module**.
      
   #. By adding the options ``--include_lane --merge_Reads`` **within module prep**, ``XICRA`` will only 
      identify one sample, merging all the corresponding files:
   
      - Sample 1: sample1_R1, sample1_R2
      
      

Option: -\ -include_all
-----------------------

In some cases, files might contain other extra information and it is necessary to use all the information of the 
file name to identify samples:

If that is the case use **-**\ **-include_all in al modules** .

If you want to merge fastq files that only differ in the final extension (_001, _002, ...):

Use options **-**\ **-merge_Reads** **-**\ **-include_all within module prep** and only 
**-**\ **-include_all in the rest of the modules**.

As an example, considering the input files:

- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz
- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz
- name1_L002_XYZ_R1_001.fastq.gz, name1_L002_XYZ_R2_001.fastq.gz
- name1_L002_XYZ_R1_002.fastq.gz, name1_L002_XYZ_R2_002.fastq.gz

   #. By adding the option ``--include_all`` in **all modules**, ``XICRA`` will identify four samples:
   
      - Sample 1: name1_L001_XYZ_R1_001, name1_L001_XYZ_R2_001
      - Sample 2: name1_L001_XYZ_R1_002, name1_L001_XYZ_R2_002
      - Sample 3: name1_L002_XYZ_R1_001, name1_L002_XYZ_R2_001
      - Sample 4: name1_L002_XYZ_R1_002, name1_L002_XYZ_R2_002
      
      Remember to use option **-**\ **-include_all within each module**.

   #. By adding the options ``--include_all --merge_Reads`` **within module prep**, ``XICRA`` will identify two samples:
   
      - Sample 1: name1_L001_XYZ_R1, name1_L001_XYZ_R2
      - Sample 2: name1_L002_XYZ_R1, name1_L002_XYZ_R2
      
      Remember to use option **-**\ **-include_all within each module**.
      




.. ############################
.. _prep-description:
.. ############################

prep
====

This module prepares fastq files from a sequencing run for later usage.

It initially checks the length of the name and advises the user to rename samples 
if exceeded a length (10 characters). Along ``XICRA`` there are a few string 
length limitations by different software that need to be sort out from the beginning 
of the process.

See additional details of this name format limitations under user-guide section: :ref:`format fastq files<format-fastq-files>`

Once sample names are checked this module allows the user to copy files into the project 
folder initiate or only link them using a symbolic link to avoid duplicated raw data. 


In addition, when multiples files have been generated for the same
sample e.g different lanes this module can merge them. It concatenates these files according the common
identifier and generates a unique file, one per paired-read if necessary, check :ref:`format fastq files<format-fastq-files>`.


If using the the ``--project`` option, this module would create a single folder for each
sample and a folder named as ``raw`` including the linked or copied fastq files. See additional
details of project folder organization :ref:`here<project-organization>`.

.. ##################
.. _run-prep:
.. ##################
How to run the prep module
--------------------------
Executing the following:

.. code-block:: sh

   XICRA prep -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA prep help

   :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA prep Input/Output

   :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
   :param --output_folder: Output folder. Name for the project folder.
   :param --single_end: Single end files. Default mode is paired-end. Default OFF.
   :param --batch: Provide this option if input is a file containing multiple paths instead a path.      
   :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
   :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
   :param --detached: Isolated mode. No project folder initiated for further steps. Default OFF.
   :param --include_lane: Include the lane tag (*L00X*) in the sample identification. See -\ -help_format for additional details. Default OFF.
   :param --include_all: Include all file name characters in the sample identification. See -\ -help_format for additional details. Default OFF.

   :type input: string
   :type output_folder: string
   :type in_sample: string 
   :type ex_sample: string

.. function:: Module XICRA prep options

   :param --threads: Number of CPUs to use. Default: 2.
   :param --copy_reads: Instead of generating symbolic links, copy files into output folder. Default OFF.
   :param --merge_Reads: Merge files corresponding to the same sample. Used in combination with -\ -include_lane and -\ -include_all will produce different results. Please check, -\ -help_format or https://xicra.readthedocs.io/en/latest/user_guide/info/info_index.html
   :param --rename: File containing original name and final name for each sample separated by comma. No need to provide a name for each pair if paired-end files. If provided with option '-\ -merge_Reads', the merged files would be renamed accordingly.
   
   :type threads: int
   :type rename: string
   
.. function:: Module XICRA prep additional information
  
   :param --help_format: Show additional help on name format for files.
   :param --help_project: Show additional help on the project scheme.
   :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/prep>`.


Output of prep for each sample
-------------------------------
After the ``prep`` module excution, a data folder will be generated. Inside it, for each sample, a new folder will be created 
(called as the sample). These sample folders will contain links to all the raw files that correspond to each sample (if -\ -copy_reads has
been selected, instead of links the copied files). 



.. include:: ../../links.inc.. ############################
.. _miRNA-description:

miRNA
=====

The ``XICRA`` pipeline provides this module, ``miRNA``, to generate the miRNA analysis.
MicroRNAs (miRNAs), a class of small non-coding RNAs, have an average length 
of 21–23 nucleotides (nt) and modulate gene expression post-transcriptionally. Most miRNA 
expression studies based on next generation sequencing (NGS) data, summarize all the reads 
mapping to a specific miRNA locus or miRNA sequence with or without mismatches and assign
it to a single miRNA entity, this is, a unique miRBase_ reference database entry (miRBase
database is a searchable database of published miRNA sequences and annotation).

However, this type of analysis neglects the fact that **not all reads are identical to the** 
**reference sequence in miRBase, which is called the canonical sequence**.
Small RNA sequencing NGS methodology has revealed that miRNAs can frequently appear in the 
form of multiple sequence variants or isoforms (termed isomiRs). Each isomiR can modulate 
gene expression post-transcriptionally. 


With the ``miRNA`` module we capture all the available information at different levels, 
the analysis can be done at miRNA or isomiR level. 

Its functionalyty is divided in three steps: 
   #. miRNA analysis
   #. Standarization of the results 
   #. Expression count matrix generation with miRTop_ (Command line tool to annotate with a standard naming miRNAs e isomiRs).


The analysis can be performed with three different softwares:

- Miraligner_: maps small RNA data to miRBase repository.
- Optimir_: algorithm for integrating available genome-wide genotype data into miRNA sequence alignment analysis.
- sRNAbench_: application for processing small-RNA data obtained from NGS platforms. 
  Unfortunately, the downloading of this tool is no longer available, thus, only users 
  with ``sRNAbench`` already installed will be able to run the ``XICRA`` analysis with it. 

According to our tests, published in this article_, ``miraligner`` is the software with 
the best performance for the miRNA analysis. 

How to run the miRNA module
---------------------------
Executing the following:

.. code-block:: sh

   XICRA miRNA -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA miRNA help
      
   :param -h --help: Show this help message and exit. 


.. function:: Module XICRA miRNA Input/Output

   :param --input: Folder containing a project or reads, according to the mode selected. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See --help_format for additional details. REQUIRED.
   :param --output_folder: Output folder.
   :param --single_end: Single end files. Default mode is paired-end. Default OFF.
   :param --batch: Provide this option if input is a file containing multiple paths instead a path.      
   :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
   :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
   :param --detached: Isolated mode. --input is a folder containing fastq reads. Provide a unique path o several using --batch option.
   :param --include_lane: Include the lane tag (*L00X*) in the sample name. See --help_format for additional details. Default OFF.
   :param --include_all: Include all characters as tag name before read pair, if any. See --help_format for additional details. Default OFF.
   :param --noTrimm: Use non-trimmed reads (or not containing '_trim' in the name).
   
   :type input: string
   :type output_folder: string
   :type in_sample: string 
   :type ex_sample: string

      
.. function:: Module XICRA miRNA options

   :param --threads: Number of CPUs to use. Default: 2. 
   :param --species: Species tag ID. Default: hsa (Homo sapiens).
   :param --database: Path to store miRNA annotation files downloaded: miRBase, miRCarta, etc.
   :param --miRNA_gff: miRBase hsa GFF file containing miRNA information.
   :param --hairpinFasta: miRNA hairpin fasta file.
   :param --matureFasta: miRNA mature fasta file.
   :param --miRBase_str: miRBase str information.
   
   :type threads: int 
   :type species: string 
   :type database: string
   :type miRNA_gff: string
   :type hairpinFasta: string
   :type matureFasta: string
   :type miRBase_str: string
   
   
.. function:: Module XICRA miRNA software

   :param --software: Software to analyze miRNAs, sRNAbench, optimir, miraligner. Provide several input if desired separated by a space. REQUIRED.
   
   :type software: string

.. function:: Module XICRA miRNA additional information
  
   :param --help_format: Show additional help on name format for files.
   :param --help_project: Show additional help on the project scheme.
   :param --help_miRNA: Show additional help on the miRNA paired-end reads process.
   :param --debug: Show additional message for debugging purposes.

- For further information of the module functionallity, check this :doc:`page <../../api/modules/miRNA>`.

Output of miRNA for each sample
-------------------------------
As the rest of the modules, the ``miRNA`` module will generate a folder in each of the 
sample directories called "miRNA". Inside this folder another two will be created for 
each of the softwares selected. For example, if we have executed ``--software optimir miraligner``
we will obtain four output folders:

- data/sampleName/miRNA/optimiR
- data/sampleName/miRNA/miraligner 
- data/sampleName/miRNA/optimiR_miRTop
- data/sampleName/miRNA/miraligner_miRTop

The first two folders will store the outputs of the corresponding softwares in their particular
format.

The folders ended in "_miRTop" will contain the results in the miRTop standarized format. 

Finally, the expression count matrix will be stored in .tsv format. Following the previous 
example, these files would be located in:

- data/sampleName/miRNA/optimiR_miRTop/counts/mirtop.tsv
- data/sampleName/miRNA/miraligner_miRTop/counts/mirtop.tsv


Expression count matrix for each sample
---------------------------------------

As a result for each sample (and software used) we will end up with a table like this, mirtop.tsv, 
called the expression count matrix. Here we can see an example of this matrix:

..
.. TODO: resize the table
..

.. csv-table::
   :file: mirtop_example.csv
   :header-rows: 1

- UID: unique identifier (UID) for each sequence defined by miRTop.
- Read: DNA sequence.
- miRNA: miRNA precursor, identifier defined by miRBase for each miRNA canonical sequence.
- Variant: variant type of each isomiR, 'NA' for the canonical sequence (checkout the miRTop `variant nomenclature`_).
- The following four columns indicate the amount of base pairs added or substracted, compared to the canonical sequence.
- SampleName: raw read count expression for this sample.

Output of miRNA, comparing samples
----------------------------------

On the other hand, as other modules, ``miRNA`` also builds an output to compare samples. In the folder
report/miRNA, three different files will be created for each software executed. For example, if we have run 
``--software miraligner``, we will obtain the following files: 

- **report/miRNA/miRNA_expression-miraligner_dup.csv**: Matrix with the number of reads of each UID of each sample that are duplicated. 
  Normally, they occur when some bases are added at the beginning and the end, so it cannot be differentiated if 
  they are 3p:+1;5p:+2 or 3p:+2;5p:+1. In those cases, they will both have the same UID. They are removed 
  (they typically have very few counts).
- **report/miRNA/miRNA_expression-miraligner.csv**: Final matrix (without the duplicated UIDs). Number of counts of each UID of each sample,
  to be further analyzed with ``R``. 
- **report/miRNA/miRNA_expression-miraligner_seq.csv**: table with the DNA sequence corresponding to each UID. 

The analysis of the matrix stored in miRNA_expression-miraligner.csv can be done at the isomiR level, differenciating by
UID, variant type or miRNA (just considering the miRNA identifier).  It can be done with the package XICRA.stats_.

.. include:: ../../links.inc************
Main modules
************

.. only:: html

    :Version: |version|
    :Date: |today|

There are multiple modules available that perform several :ref:`steps<pipeline-scheme>` and generate multiple results. 

``XICRA`` modules require several command-line arguments and options to run. There are a number of shared 
:ref:`arguments<shared-arguments>` among all modules and some others specific of each module and specified 
within each one.

.. toctree::
   :maxdepth: 1
   
   config.rst
   prep.rst
   QC.rst
   trimm.rst
   join.rst
   biotype.rst
   miRNA.rst
   
.. _shared-arguments:

Command-line shared arguments
*****************************

Here we include a brief description of the shared command-line arguments for some of ``XICRA`` modules.
   
   **Mode:**
   
   --project         Project mode. Requires as ``--input`` a folder containing an initialized ``BacterialTyper`` project [Default].
   
   --detached        Isolated mode. ``--input`` is a folder containining fastq reads. Provide a unique path o several using ``--batch`` option
   
   **Input/Output:**
   
   --input  string         Folder containing a project or reads, according to the mode selected. Files could be ``.fastq/.fq`` or ``fastq.gz/.fq.gz``. See ``--help_format`` for additional details.
   
   --single_end         Single end files [Default OFF]. Default mode is paired-end.
   
   --batch        Provide this option if input is a file containing multiple paths instead a path.
   
   --in_sample string         File containing a list of samples to include (one per line) from input folder(s) [Default OFF].
   
   --ex_sample string         File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].
   
   **Options:**
   
   --threads int        Number of CPUs to use [Default: 2]
   
   **Additional information:**
   
   --debug        Show additional message for debugging purposes.


.. _project-organization:

Details of the XICRA project folder
***********************************

TODO.. ############################
.. _qc-description:
.. ############################

QC
===

This module calls fastqc_, a quality check program, to analyze the 
quality of each sample. 

By default, creates a final MultiQC_ report with all the samples. It is useful 
to check if there are outliers among the input samples. It can be disabled using
the option -\ -skip_report. 

.. ##################
.. _run-qc:
.. ##################
How to run the QC module
--------------------------
Executing the following:

.. code-block:: sh

   XICRA QC -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA QC help

    :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA QC Input/Output

    :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
    :param --output_folder: Output folder. Name for the project folder.
    :param --batch: Provide this option if input is a file containing multiple paths instead a path.
    :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
    :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
    :param --detached: Isolated mode. No project folder initiated for further steps. Default OFF.
    :param --include_lane: Include the lane tag (*L00X*) in the sample identification. See -\ -help_format for additional details. Default OFF.
    :param --include_all: Include all file name characters in the sample identification. See -\ -help_format for additional details. Default OFF.      

    :type input: string
    :type output_folder: string
    :type in_sample: string 
    :type ex_sample: string

.. function:: Module XICRA QC options

    :param --skip_report: Do not report statistics using MultiQC report module. Default OFF.
    :param --threads: Number of CPUs to use. Default: 2.
   
    :type threads: int
    :type rename: string

.. function:: Module XICRA QC additional information
  
    :param --help_format: Show additional help on name format for files.
    :param --help_project: Show additional help on the project scheme.
    :param --help_multiqc: Show additional help on the multiQC module.
    :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/QC>`.


Output of QC
------------
After the ``QC`` module excution,in the data folder for each sample, a new directory will be created.
It will be called 'fastq' and will contain the the quality analysis.

In addition, if -\ -skip_report was not selected, a folder called 'report' will also be generated. In
this folder, the different reports that may be created after the execution of other modules (as ``trimm``)
will be stored. These reports will always show information of all the samples. In this case, a subfolder 
called 'fastqc' will be built with the MultiQC report in order to visualize the quality 
analysis of all samples.


.. include:: ../../links.inc.. ############################
.. _config-description:

config
======

The ``XICRA`` pipeline provides this module ``config`` to help the user
to check if the multiple dependencies and requirements are fulfilled.

We encourage ``XICRA`` users to run this module after the installation to
check whether the multiple requirements and dependencies are correctly installed.

.. seealso:: Additional information on ``XICRA`` configuration and requirements

   - :doc:`Installation <../installation/installing>` 
   
   - :doc:`Requirements <../installation/requirements>` 
   
   - :doc:`config module (API) <../../api/modules/config>`


How to run the config module
----------------------------

Once you have installed ``XICRA`` you should be able to run in the command line the pipeline.

If you type ``XICRA`` you should see a prompt with the different modules available. Following 
the pipeline name type the module of interest, in this case, ``config``. As an example:

.. code-block:: sh

   XICRA config -h

The different options and parameters for this module should appear in the command line prompt:


.. function:: Module XICRA config
   
   :param -h, --help: Show this help message and exit. 
   :param --debug: Show additional message for debugging purposes.

Basically, this ``XICRA config`` module allows the user to check if the requirements are fulfilled. 

.. code-block:: sh

   XICRA config

.. include:: ../../links.inc.. ############################
.. _biotype-description:
.. ############################

biotype
=======
This module generates a RNA biotype analysis. The aim of this computation 
is to check if there are samples with a very different configuration, outliers
that show different proportions of uniquely mapped, multimapped or no mapped
reads or with different quantities of miRNA, misc_RNA,... than the rest.
If this happens, it could be due to possible differences in sample manipulation, 
extraction, library preparation, sequencing, etc. Those samples should be 
excluded. 

The module is divided in two steps: 

1. Mapping the reads: performed with STAR_ software.

2. Feature counts: perforfed with featureCounts_.

The output of this module is ‘descriptive’, thus, we won’t need the output to continue 
with the analysis, it is an informative step to know the RNA types that are present
in each sample.  

**Note that**: the mapping process performed with STAR software requires high RAM values. 
We are working on the implementation of alternatives to be able to execute the ``biotype``
module in computers with less RAM capacity. 

.. ##################
.. _run-biotype:
.. ##################
How to run the biotype module
------------------------------
Executing the following:

.. code-block:: sh

   XICRA biotype -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA biotype help

    :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA QC Input/Output

    :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
    :param --output_folder: Output folder. Name for the project folder.
    :param --single_end: Single end files. Default OFF. Default mode is paired-end.
    :param --batch: Provide this option if input is a file containing multiple paths instead a path.
    :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
    :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
    :param --detached: Isolated mode. No project folder initiated for further steps. Default OFF.
    :param --include_lane: Include the lane tag (*L00X*) in the sample identification. See -\ -help_format for additional details. Default OFF.
    :param --include_all: Include all file name characters in the sample identification. See -\ -help_format for additional details. Default OFF.      

    :type input: string
    :type output_folder: string
    :type in_sample: string 
    :type ex_sample: string

.. function:: Module XICRA biotype options

    :param --threads: Number of CPUs to use. Default: 2.
    :param --annotation: Reference genome annotation in GTF format.
    :param --limitRAM: limitRAM parameter for STAR mapping. Default 20 Gbytes.
    :param --noTrim: Use non-trimmed reads [or not containing '_trim' in the name].
    :param --skip_report: Do not report statistics using MultiQC report module. Default OFF. See details in -\ -help_multiqc

   
    :type threads: int
    :type annotation: string
    :type limitRAM: int

.. function:: Module XICRA biotype parameters

    :param --no_multiMapping: Set NO to counting multimapping in the feature count.By default, multimapping reads are allowed. Default: False
    :param --stranded STRANDED: Select if reads are stranded [1], reverse stranded [2] or non-stranded [0], Default: 0.

.. function:: Module XICRA biotype reference genome

    :param --fasta: Reference genome to map reads.
    :param --genomeDir: STAR genomeDir for reference genome.


.. function:: Module XICRA biotype additional information
  
    :param --help_format: Show additional help on name format for files.
    :param --help_project: Show additional help on the project scheme.
    :param --help_RNAbiotype: Show additional help on the RNAbiotype paired-end reads process.
    :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/biotype>`.

Output of biotype
-----------------
Inside the data folder of each sample, a 'map' directory will be generated containing a
report from the mapping of MultiQC. After that, a final report will also be created in the
'report' folder with the featureCounts information of all samples. 

.. include:: ../../links.inc.. ############################
.. _trimm-description:
.. ############################

trimm
=====
This module trimms sequencing adapters that could be present in next generation sequencing 
files. Adapters have to be ligated to every single DNA molecule during library preparation. 
For Illumina short read sequencing, the corresponding protocols involve (in most cases) a 
fragmentation step, followed by the ligation of certain oligonucleotides to the 5’ and 
3’ ends. These 5' and 3' adapter sequences have important functions in Illumina sequencing, 
since they hold barcoding sequences, forward/reverse primers (for paired-end sequencing) 
and the important binding sequences for immobilizing the fragments to the flowcell and 
allowing bridge-amplification. 

.. image:: ../../images/modules/trimm/trimm_reads.jpeg
   :width: 300pt
   :align: center
   
In `Illumina sequencing`_, adapter sequences will only occur at the 3' end of the read and only 
if the read length is longer than the insert size:    

.. image:: ../../images/modules/trimm/trimm_length.jpeg
   :width: 300pt
   :align: center

However, in the case of miRNAs, there will be adapter contamination in the 3’ reads but also
in the 5', due to their short nature.  As the adapters are synthetic they need to be removed.
Note that, the 3’ adapter of the R2 will be the same as the 5’ of R1, and the 5' of R2 will
be the same as the 3' of R1.

The adapter sequences must be introduced by the user to be trimmed, using the parameters -\ -adapters_a and 
-\ -adapters_A, for the adapters attached to the 3' of read 1 and read 2 correspondingly, see this :ref:`section<run-trimm>`. 
The process is done by the software Cutadapt_. 

.. ##################
.. _run-trimm:
.. ##################
How to run the trimm module
---------------------------
Executing the following:

.. code-block:: sh

   XICRA trimm -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA trimm help

   :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA trimm Input/Output

   :param --input: Folder containing a project or reads, according to the mode selected. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See --help_format for additional details. REQUIRED.
   :param --output_folder: Output folder.
   :param --single_end: Single end files. Default mode is paired-end. Default OFF.
   :param --batch: Provide this option if input is a file containing multiple paths instead a path.      
   :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
   :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
   :param --detached: Isolated mode. --input is a folder containing fastq reads. Provide a unique path o several using --batch option.
   :param --include_lane: Include the lane tag (*L00X*) in the sample identification. See -\ -help_format for additional details. Default OFF.
   :param --include_all: IInclude all file name characters in the sample identification. See -\ -help_format for additional details. Default OFF.
   
   :type input: string
   :type output_folder: string
   :type in_sample: string 
   :type ex_sample: string

.. function:: Module XICRA trimm options

   :param --adapters_a: Sequence of an adapter ligated to the 3' end (of read 1). See --help_trimm_adapters for further information.
   :param --adapters_A: Sequence of an adapter ligated to the 3' read in pair (of read 2). See -\ -help_trimm_adapters for further information.
   :param --extra: provide extra options for cutadapt trimming process. See -\ -help_trimm_adapters for further information.
   :param --skip_report: Do not report statistics using MultiQC report module. [Default OFF]. See details in --help_multiqc
   :param --threads: Number of CPUs to use. Default: 2. 
   
   :type threads: int 
   :type adapters_a: string 
   :type adapters_A: string
   :type extra: string
   
.. function:: Module XICRA trimm additional information
  
   :param --help_format: Show additional help on name format for files.
   :param --help_project: Show additional help on the project scheme.
   :param --help_trimm_adapters: Show additional help of the trimm module.
   :param --help_multiqc: Show additional help on the multiQC module.
   :param --debug: Show additional message for debugging purposes.

- For further information of the module functionallity, check this :doc:`page <../../api/modules/trimm>`.

Output of trimm for each sample
-------------------------------
After the ``trimm`` module excution, for each sample, a new fastq file is generated in the in the "trimm" folder of the sample 
with the same name as the raw one + "_trimm". 

It generates by default a MultiQC report in the folder "report/trimm", this report may be useful to detect outliers or 
if there are similar amounts of trimmed sequences for all samples. 

.. include:: ../../links.inc.. ########################
.. _software-details:
.. ########################

.. #
.. TODO: complete this page
.. #

Software Information
********************
      
This is a general guide for the software that we employed along the ``XICRA`` pipeline.

.. contents::

.. ########################
.. _fastqc-description:
.. ########################

FastQC
======

FastQC :cite:`fastqc` is a quality control tool for high throughput sequence data. It aims to provide a simple 
way to do some quality control checks on raw sequence data coming from high throughput sequencing 
pipelines. It provides a modular set of analysis which you can use to give a quick impression of 
whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are

- Import of data from BAM, SAM or FastQ files (any variant)

- Providing a quick overview to tell you in which areas there may be problems

- Summary graphs and tables to quickly assess your data

- Export of results to an HTML based permanent report

- Offline operation to allow automated generation of reports without running the interactive application

Read further information about FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


.. ########################
.. _cutadapt-description:
.. ########################

Cutadapt
========


.. ########################
.. _STAR-description:
.. ########################

STAR
====


.. ########################
.. _featureCounts-description:
.. ########################

Feature counts
==============



.. ########################
.. _mirTop-description:
.. ########################

mirTop
======


.. ########################
.. _miraligner-description:
.. ########################

Miraligner
==========

.. ########################
.. _mintmap-description:
.. ########################

MINTmap
=======


.. ########################
.. _sRNAbench-description:
.. ########################

sRNAbench
=========


.. ########################
.. _optimir-description:
.. ########################

Optimir
=======
.. _python-requirements:

Python requirements details
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is the content of file :file:`XICRA/config/python/python_requirements.txt` containing
all python modules and versions required by ``XICRA``.

.. include:: ../../../../XICRA/config/python/python_requirements.txt
   :literal:.. ########################
.. _installing:
.. ########################

Installation
************

.. toctree::
   :hidden:
   
   requirements.rst
   python-environment.rst
      
.. contents::

This is an installation guide for the ``XICRA`` pipeline. 

``XICRA`` is a pipeline composed of multiple lines of code 
that calls and automatizes the analysis from multiple bioinformatic 
tools. Several dependencies (python, perl, additional software and their dependencies 
and third-party software) are required for the ``XICRA`` analysis. 

You first need to get ``XICRA`` (from different sources available): 
you can either install ``XICRA`` latest version using Conda_ and `Python pip`_ 
or if you want to contribute or see additional details of the project, it's 
recommended you :ref:`install the latest development version<install-from-source>`.

Once you get the code, before running ``XICRA`` you must make sure 
you have all the dependencies fulfilled from section 
:ref:`Requirements and dependencies<Requirements-dependencies>` using the 
``XICRA config`` module.


.. ########################
.. _Requirements-dependencies:
.. ########################

Requirements and dependencies
=============================

.. ########################
.. _system-requirements:
.. ########################

System requirements
-------------------

Take into account that you may need some system requirements to install ``XICRA`` already available
in your system such as ``python3``, ``python3-dev`` and ``python3-venv`` or ``build-essential`` libraries among others. 

Check and install them by typing:

.. code-block:: sh
   
   sudo apt install python3 python3-dev python-dev python3-venv python-pip
   sudo apt install build-essential libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev

``XICRA`` will require ``python`` v3.6 and ``java`` (we tested in ``openjdk`` 14 2020-03-17).

XICRA requirements
------------------

We include details on the different modules required by ``XICRA`` to successfully run here:

.. toctree:: 
   requirements.rst

.. ########################
.. _install-XICRA:
.. ########################

Installing XICRA
================

If you want to run a stable version of ``XICRA``, the installation can be done with 
``conda`` or ``pip`` packages following instructions in :ref:`Installing from conda<install-from-conda>`
and :ref:`Installing from pip<install-from-pip>`, correspondingly. We encourage you
to install ``XICRA`` and all dependencies using the ``conda`` environment.

..
.. Is it necessary the second part of the installation, From source??
..

On the other hand if you are interested in contributing to ``XICRA`` development, 
running the latest source code, or just like to build everything yourself, you
should follow the :ref:`Installing from source<install-from-source>` instructions. 

Additionally, there are a number of dependencies that might be necessary to 
install or check within your system in both cases. Choose the appopiate choice 
according to your intalling ``XICRA`` option selected.

.. ####################
.. _install-from-conda:
.. ####################

Installing from conda
---------------------

Unfortunately, a couple of executables are not available neither as a ``conda`` or 
``pip`` packages. These packages are ``miraligner`` and ``sRNAbench``. We have generated a
shell script to retrieve ``miraligner`` and include it within your ``conda`` environment.
However, ``sRNAbench`` is no longer available to be downloaded, thus, only the users
with the software already installed will be able to run the ``XICRA`` analysis with 
``sRNAbench``.
To create a new ``conda`` environment, install third party software, install ``XICRA``
and missing dependencies, do as follows:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git

   ## move to folder
   cd XICRA

   ## create conda environemt
   conda env create -n XICRA -f XICRA_pip/devel/conda/requirements.txt

   ## activate
   conda activate XICRA

   ## install latest python code
   pip install XICRA

   ## install missing software
   sh XICRA_pip/XICRA/config/software/installer.sh
   
To check everything is fine, try executing the ``config`` module:

.. code-block:: sh

   XICRA config
  
We additionally provide a supplementary R package for parsing and 
plotting some ``XICRA`` results, called XICRA.stats_. See additional details here.

Install it in ``R` using:

.. code-block:: sh

   ## Install XICRA.stats version from GitHub:
   ## install.packages("devtools")
   devtools::install_github("HCGB-IGTP/XICRA.stats")
  
.. ##################
.. _install-from-pip:
.. ##################

Installing from pip
-------------------

If you are not using a ``conda`` environment as you might have previously installed all 
dependencies, we encourage you to create a ``python`` environment containing all ``python`` 
modules required for ``XICRA``.  

.. code-block:: sh

   ## create enviroment
   python3 -m venv XICRA_env

   ## activate it
   source XICRA_env/bin/activate

   ## install XICRA and dependencies
   pip install XICRA

   ## execute XICRA
   XICRA -h
  
   
Follow additional `pip installing instructions`_ to learn about installing packages.

Also, as ``XICRA`` relies in multiple dependencies, external packages and third-party 
software, we encourage you to once you install ``XICRA`` using ``pip``. Check for 
dependencies using the ``XICRA config`` module. See
details in ``XICRA config`` module :ref:`section<config-description>`.


.. ##################
.. _install-from-source:
.. ##################

Installing from source
----------------------

Under some circumstancies (develop, bugs fixed, etc) you might be interested 
in obtaining the latest code version. Take into account, that you will need to install 
dependencies and fulfill requirements to have a working distribution. 

.. ##############
.. _get-git-code:
.. ##############

Get source code, for developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have included files in the folder ``devel/conda`` for the build and 
configuration of the ``conda`` package.
The ``XICRA`` project uses git_ as a version control system. To get the code, 
you can grab the latest version from the `XICRA github`_ website, and 
follow the :ref:`Install XICRA from source<install-XICRA-source>` instructions.

Using the command-line, check you have a working distribution of ``git`` by typing 
``git --help`` or install it by typing:

.. code-block:: sh

   sudo apt update
   sudo apt upgrade
   sudo apt install git

Once you have ``git`` installed, to create a new conda environment, install third party
software, install ``XICRA`` and missing dependencies for ``XICRA`` development, do 
as follows:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git

   ## move to folder XICRA_pip
   cd XICRA/XICRA_pip

   ## create conda environemt
   conda env create -n XICRA -f ./devel/conda/requirements.txt

   ## activate
   conda activate XICRA

   ## install latest python code
   pip install -r ./devel/pypi/requirements.txt
   pip install -e .

   ## install missing software
   sh ./XICRA/config/software/installer.sh
   

.. ###############################
.. _install-XICRA-source:
.. ###############################

Install XICRA from source
^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have ``XICRA`` source code available you only need to include
the ``XICRA`` folder and main script within your path. 

.. code-block:: sh

   export PYTHONPATH=$PYTHONPATH":"$PWD"/XICRA"

   export PATH=$PATH":"$PWD"/XICRA.py"

Take into account that before running ``XICRA`` you have to make sure you have all the 
dependencies fulfilled from section :ref:`Requirements and dependencies<Requirements-dependencies>`.
You can either install them yourself, use appropiate scripts for this purpose or use the ``XICRA config``
module to check, update and install all dependencies required.

.. ###########
.. include:: ../../links.inc.. ########################
.. _python-modules-required:
.. ########################

Python modules
--------------

There are several extra python module requirements that 
are needed and are summarized in the following table: 

.. csv-table::
   :header: "Module", "Version"
   :file: ../../../../XICRA/config/python/python_requirements.csv

These modules might have extra dependencies. Details of the list of 
all modules required are listed in :file:`XICRA/config/python/python_requirements.txt`. 
And accessible here:

.. toctree:: 
   python-requirements.rst

Although these dependencies will be fulfilled during the ``XICRA``
installation with ``pip``, you might be interested in installing them yourself. 

Using ``pip`` we can install them all at a glance. 

.. code-block:: sh

   pip install -r ./XICRA/config/python/python_requirements.txt

But again, following installation recommendations, we encourage you to create and install them 
within a virtual environment (See section: :ref:`Python environment<virtual-env-XICRA>` 
section for details).

You can test the presence of these ``python`` modules using the ``XICRA config`` module. 
Once you identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``PYTHONPATH`` or environment or you can execute the ``XICRA config`` 
with ``install`` option.

.. ######################
.. _soft-dependencies:
.. ######################

Software dependencies
---------------------

Also, several software packages are also required. They are listed in
:file:`XICRA/config/software/dependencies.csv`, which is shown below:

.. csv-table::
   :header-rows: 1 
   :file: ../../../../XICRA/config/software/soft_dependencies.csv

Most of the software are common software that any person doing bioinformatics should have, so
you might have already available within your system. However, installing ``XICRA`` using 
``conda`` all the dependecies will be correctky installed.   

You can test for any missing software dependencies using the ``XICRA config`` module. Once you 
identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``$PATH`` or you can execute the ``X config`` 
with ``install`` option.


.. #### Include links
.. include:: ../../links.inc