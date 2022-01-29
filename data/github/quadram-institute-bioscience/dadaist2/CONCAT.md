# dadaist2

[![Dadaist2 logo](docs/img/dadaist2.png)](https://quadram-institute-bioscience.github.io/dadaist2/)

[![release](https://img.shields.io/github/v/release/quadram-institute-bioscience/dadaist2?label=github%20release)](https://github.com/quadram-institute-bioscience/dadaist2/releases)
[![version](
https://anaconda.org/bioconda/dadaist2/badges/version.svg)](https://bioconda.github.io/recipes/dadaist2/README.html)
[![Conda](https://anaconda.org/bioconda/dadaist2/badges/downloads.svg)](https://bioconda.github.io/recipes/dadaist2/README.html)
[![Build Status](https://www.travis-ci.com/quadram-institute-bioscience/dadaist2.svg?branch=master)](https://www.travis-ci.com/quadram-institute-bioscience/dadaist2)
[![Dadaist-CI](https://github.com/quadram-institute-bioscience/dadaist2/actions/workflows/blank.yml/badge.svg)](https://github.com/quadram-institute-bioscience/dadaist2/actions/workflows/blank.yml)

Standalone wrapper for [DADA2](https://benjjneb.github.io/dada2/index.html) package, to quickly generate a feature table and a
set of representative sequences from a folder with Paired End Illumina reads.

# [Documentation and tutorials](https://quadram-institute-bioscience.github.io/dadaist2)

Please check the [online documentation](https://quadram-institute-bioscience.github.io/dadaist2) and tutorials
for installation and usage notes.

## Authors
* Andrea Telatin, Quadram Institute Bioscience, UK
* Rebecca Ansorge, Quadram Institute Bioscience, UK
* Giovanni Birolo, University of Turin, Italy
* **dadaist2** is the Perl wrapper
* R scripts from: [https://github.com/qiime2/q2-dada2](https://github.com/qiime2/q2-dada2)
# DADA2 validation

DADA2 provides a [tutorial](https://benjjneb.github.io/dada2/tutorial_1_8.html)
based on [Mothur SOP dataset](https://mothur.org/wiki/miseq_sop/).


### Get the reads

We will store the reads in a directory called `subsample`:
```
wget "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
unzip miseqsopdata.zip
mv MiSeq_SOP subsample/
```

### Run DADA2 natively
```bash
# Set DADAIST as the path to the repository of Dadaist2
mkdir -p dada-subsample
Rscript $DADAIST/test/miseq-sop-compare/dada2-sop.R subsample/ dada-subsample/
```

### Run DADA2 via Qiime2 2021.2
```
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path subsample-for-qiime/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 154 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-q 2 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza \
  --p-n-threads 32
```

### Run DADA2 via Dadaist2

The initial test has been done with version **0.8.0**.

```bash
# Default parameters (fastp)
dadaist2 -r $DADAIST/refs/silva_v138.fa.gz \
  -i subsample/ \
  -o dadaist_0.8.0/subsample-silva-default \
  -m metadata.tsv \
  -t 32

# Skip primer trimming and pool samples for DADA (as DADA2 workflow does)
dadaist2 -r $DADAIST/refs/silva_v138.fa.gz \
  -i subsample/ \
  -o dadaist_0.8.0/subsample-silva-notrim-pool \
  --no-trim --dada-pool \
  -m metadata.tsv \
  -t 32
```


## Comparison of the results

The three pipelines produced similar feature tables
and representative sequences. The FASTA file stats are:

```text
┌────────────────────┬──────┬──────────┬───────┬─────┬─────┬─────┬──── ───┬─────┬─────┐
│ File               │ #Seq │ Total bp │ Avg   │ N50 │ N75 │ N90 │ auN    │ Min │ Max │
├────────────────────┼──────┼──────────┼───────┼─────┼─────┼─────┼────────┼─────┼─────┤
│ dada2_relabel      │ 137  │ 34644    │ 252.9 │ 253 │ 253 │ 252 │ 9.2383 │ 251 │ 255 │
│ dadaist080_relabel │ 138  │ 34902    │ 252.9 │ 253 │ 253 │ 252 │ 9.1701 │ 251 │ 255 │
│ qiime_relabel      │ 136  │ 31671    │ 232.9 │ 233 │ 233 │ 232 │ 8.5710 │ 231 │ 235 │
└────────────────────┴──────┴──────────┴───────┴─────┴─────┴─────┴────────┴─────┴─────┘
```

Qiime2 produced slighly shorter sequences because of trimming.

The files are merged and clustered using the `compare.sh` script.

The output is:
```text
      1   dada
      2   dada;dadaist
      2   dadaist;qiime
     11   dadaist
     15   dada;qiime
    119   dada;dadaist;qiime
```

80% of the sequences are identical among the samples---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior, is needed a link to the samples causing the problem

**Expected behavior**
A clear and concise description of what you expected to happen.

 **Environment:**
 - OS: [e.g. Ubuntu 18.04]
- Dadaist2 version:
 
**Additional context**
It is useful to:
1) attach the log file 
2) if the problem is with dadaist2, run it appending `--debug 2>&1 | tee dadaist-debug.log` and attach the second log
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.
 
**Additional context**
Add any other context or screenshots about the feature request here.
# Using some components from Dadaist2 in custom pipelines

# External docs

This directory contains pod files for non Perl binaries---
sort: 5
permalink: /custom-workflows
---

# Advanced usage 

## Use in HPC clusters

The supported installation of _Dadaist2_ is via Bioconda:
* The package can be natively installed via `conda` in the user's home
* A singularity (or Docker) image can be generated as suggested in the [installation page](../installation).
  

## Dadaist2 as a module provider

Dadaist2 has been released as a set of wrappers to allow implementing some of them
in existing pipelines. 

A minimal example is provided in the `nextflow/simple.nf` Nextflow script 
([link](https://github.com/quadram-institute-bioscience/dadaist2/tree/master/nextflow)), where
we delegate to Nextflow the parallelisation of the input reads trimming with `cutadapt`.

### Example

Running a workflow using Dadaist tools:
```bash
nextflow run simple.nf  -with-singularity dadaist2.simg \
  --reads "data/16S/*_R{1,2}_001.fastq.gz" --ref "refs/SILVA_SSU_r138_2019.RData"  
```

The minimal workflow runs parallelising adaptor trimming and collecting the results
for Dadaist main module:
```text
N E X T F L O W  ~  version 19.10.0
Launching `simple.nf` [angry_celsius] - revision: 053df3283c
Example pipeline
 =======================================
 taxonomy db  : refs/SILVA_SSU_r138_2019.RData
 reads        : data/16S/*_R{1,2}_001.fastq.gz
 outdir       : dadaist
executor >  slurm (4)
[81/39e622] process > cutadapt (F99_S0_L001) [100%] 3 of 3, cached: 3 ✔
[19/e589db] process > dada (1)               [100%] 1 of 1, cached: 1 ✔
[5a/893a7d] process > phyloseq (1)           [100%] 1 of 1 ✔
[ae/bd2ecb] process > normalize-alpha (1)    [100%] 1 of 1 ✔
```

Note that `-with-singularity` makes use of our container that provides:

* All Dadaist wrappers 
* R with commonly used libraries such as DADA2, PhyloSeq, Microbiome, vegan...
* Tools like VSEARCH, cutadapt, fastp, seqfu, multiqc...---
sort: 2
permalink: /introduction
---

# Dadaist2 features

## Fast track to R and numerical ecology

Dadaist is a simple and modular toolkit to streamline the generation of
plots and analyses for microbiome studies, based on the popular DADA2
algorithm for reads denoising.

What makes Dadaist an interesting alternative to other suites is the focus on reproducible downstream analyses (thanks to the automatic generation of a PhyloSeq object and the preparation of files ready to be analysed with MicrobiomeAnalyst or Rhea)

* Generation of a [*PhyloSeq*](https://joey711.github.io/phyloseq/) object, for immediate usage in R
* Possibility to run in the pipeline a _custom R script_ that starts from the PhyloSeq object
* Generation of [*MicrobiomeAnalyst*](https://www.microbiomeanalyst.ca)-compatible files. MicrobiomeAnalyst provides a _web-interface_ to perform a broad range of visualizations and analyses.
* Generation of [*Rhea*](https://lagkouvardos.github.io/Rhea/)-compatible files. Rhea is a standardized set of scripts "_designed to help easy implementation by users_".
* Generation of [*MultiQC*](https://multiqc.info) ready report


<img src="img/flow.png">

## Mitigation of Cross-talk noise

Using [dadaist2-crosstalk]({{ '/pages/dadaist2-crosstalk.html' | relative_url }}) it is possible to reduce
the noise introduced by the spillover of reads from one sample to another sample in the same sequencing lane,
using the [_UNCROSS2_ algorithm](https://www.biorxiv.org/content/10.1101/400762v1.full).

## Long targets workflow

A custom DADA2 workflow that does not rely on read merging to identify molecular species
longer than the sequencing reads length. See some [ITS notes]({{ '/notes/2_ITS.html' | relative_url }}).

## Advanced logs and notifications

Dadaist2 is both a collection of tools (to create your own pipeline, for example using NextFlow) and a standalone 
pipeline designed to be easy to run from a local computer. 


![popup](img/popup.png)

* Colored terminal output to follow the progress of the pipeline
* Optional notification popups to follow the progress of the major steps while doing something else
* Regular text logs are also collected in an [easy to browse HTML report](example-log.html). 

---
sort: 3
permalink: /tutorial
---

# Dadaist2: a first tutorial

## Get ready

[Install Dadaist2](/installation) and activate the Miniconda environment (if needed).

For this tutorial we will analyze three small samples present in the repository, so
we will download the repository as well:

```bash
git clone https://github.com/quadram-institute-bioscience/dadaist2
cd dadaist2
```

Let's start checking the number of reads per sample. The reads are in `data/16S`:

```bash
seqfu count --basename data/16S/*.gz
```

This will tell the number of reads, checking that the forward (R1) and reverse (R2)
pair have the same amount of reads. This should be the output produced:

```
F99_S0_L001_R1_001.fastq.gz	4553	Paired
A01_S0_L001_R1_001.fastq.gz	6137	Paired
A02_S0_L001_R1_001.fastq.gz	5414	Paired
```

:warning: Sample names should not begin with numbers (eg: `1_R1.fastq.gz`). This will be
prohibited in a future release, and will trigger a warning in the current release.

## Download a reference database

Dadaist2 provides a convenient tool to download some pre-formatted reference databases.
To have a list of the available to download:
```bash
dadaist2-getdb --list
```

This should produce something like:
```text
dada2-hitdb: HITdb is a reference taxonomy for Human Intestinal 16S rRNA genes
dada2-rdp-species-16: RDP taxonomic training data formatted for DADA2 (RDP trainset 16/release 11.5)
dada2-rdp-train-16: RDP taxonomic training data formatted for DADA2 (RDP trainset 16/release 11.5)
dada2-silva-138: SILVA release 138
dada2-silva-species-138: SILVA release 138 (species)
decipher-gtdb95: GTDB
decipher-silva-138: SILVA release 138 (Decipher)
decipher-unite-2020: UNITE 2020 (Decipher)
testset: FASTQ input, Small 16S dataset to test the suite
```
Some reference files are for DADA2, others are for DECIPHER. Dadaist2 will automatically select the
correct classifier based on the extension of the reference database.

The keyword before the column is the dataset name, to download it we need to choose a destination directory,
we can for example make a `refs` subdirectory:

```bash
mkdir -p refs
dadaist2-getdb -d decipher-silva-138 -o ./refs
```

This will place `SILVA_SSU_r138_2019.RData` in the output directory.

## Generate a mapping file

A metadata file is not mandatory, but it's easy to generate one to be extended with more columns if needed.

```bash
dadaist2-metadata -i data/16S > metadata.tsv
```
This should generate a file called _metadata.tsv_ with the following content:
```text
#SampleID  Files
A01        A01_S0_L001_R1_001.fastq.gz,A01_S0_L001_R2_001.fastq.gz
A02        A02_S0_L001_R1_001.fastq.gz,A02_S0_L001_R2_001.fastq.gz
F99        F99_S0_L001_R1_001.fastq.gz,F99_S0_L001_R2_001.fastq.gz
```

## Run the analysis

Dadaist2 provides options to:
* select the QC strategy (fastp, cutadapt of seqfu)
* select the taxonomy classifier (DECIPHER or DADA2 naive classifier)
* adjust various steps via command line parameters


As a first run, we recommend using the default parameters:
```bash
dadaist2 -i data/16S/ -o example-output -d refs/SILVA_SSU_r138_2019.RData -t 8 -m metadata.tsv
```

Briefly:
* `-i` points to the input directory containing paired end reads (by default recognised by `_R1` and `_R2` tags, but this can be customised)
* `-o` is the output directory
* `-d` is the reference database in DADA2 or DECIPHER format (we downloaded a DECIPHER database)
* `-m` link to the metadata file (if not supplied a blank one will be generated and used)
* `-t` is the number of processing threads

## The output directory

Notable files:
* **rep-seqs.fasta** representative sequences (ASVs) in FASTA format
* **rep-seqs-tax.fasta** representative sequences (ASVs) in FASTA format, with taxonomy labels as comments
* **feature-table.tsv** table of raw counts (after cross-talk removal if specified)
* **taxonomy.tsv** a text file with the taxonomy of each ASV (used to add the labels to the _rep-seqs-tax.fasta_)
* copy of the **metadata.tsv** file

Subdirectories:
* **MicrobiomeAnalyst** a set of files formatted to be used with the online (also available offline as R package) software [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/upload/OtuUploadView.xhtml).
* **Rhea** a directory with files to be used with the [Rhea pipeline](https://lagkouvardos.github.io/Rhea/), as well as some pre-calculated outputs (Normalization and Alpha diversity are done by default, as they don't require knowledge about metadata categories)
* **R** a directory with the PhyloSeq object
---
sort: 1
permalink: /installation
---

# Installation

## Install via Miniconda

The easiest (and recommended) way to install **dadaist2** is from the BioConda repository.
This requires  _Miniconda_ installed ([how to install it](https://docs.conda.io/en/latest/miniconda.html)).
We will first install _mamba_, that makes the installation faster (naturally you can skip the mamba installation if
you already use it).

```bash
conda install -y -c conda-forge mamba
mamba install -y -c conda-forge -c bioconda dadaist2
```

If you want to keep dadaist2 and its dependencies in a separate environment (**recommended**):

```bash
conda install -y -c conda-forge mamba
mamba create -n dadaist -c conda-forge -c bioconda dadaist2
# Then type `conda activate dadaist` to use it
```

## Install via Conda environment files

```note
This is the recommended way for reproducible analyses
```

When installing a package via Miniconda, some of its dependecies might change
from one installation to the other. To ensure the highest level or reproducibility
we are now offering curated YAML files that can be used to install the stable versions
of dadaist.

A list of environment files are available
[in the **env** directory](https://github.com/quadram-institute-bioscience/dadaist2/tree/master/env),
or you can download the last as shown below.

```bash
# Change URL as appropriate selecting from the list in the link above
wget -O dadaist2.yaml "https://quadram-institute-bioscience.github.io/dadaist2/dadaist2-$(uname).yaml"
conda env create --file dadaist2.yaml -n dadaist2
```

## Developmental snapshot

We recomment to use Miniconda also to test the last developmental snapshot, as Miniconda
can create an environment with all the required dependencies, then the binaries from the
repository can be used instead:

```bash
mamba create -n dadaist-last -c conda-forge -c bioconda --only-deps dadaist2-full
git clone https://github.com/quadram-institute-bioscience/dadaist2
export PATH="$PWD"/dadaist2/bin:"$PATH"
```

---

## Docker image

Dadaist2 is available from [DockerHub](https://hub.docker.com/r/andreatelatin/dadaist2) and
the image can be pulled via:

```
sudo docker pull andreatelatin/dadaist2:last
```

## Advanced topics

### Singularity definition
To manually build an image with the latest version from Bioconda the following definition file can be saved as `dadaist-stable.def`:

```singularity
Bootstrap: docker
From: centos:centos7.6.1810

%environment
    source /opt/software/conda/bin/activate /opt/software/conda_env
    export PATH=/opt/software/dadaist2/bin:$PATH
    PERL5LIB=''
    LANG=C
         
%post
    yum -y install epel-release wget which nano curl zlib-devel git free
    yum -y groupinstall "Development Tools"
    mkdir -p /opt/software
    cd /opt/software    
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
    sh ./Miniconda3-latest-Linux-x86_64.sh -p /opt/software/conda -b
    /opt/software/conda/bin/conda install -y mamba
    /opt/software/conda/bin/mamba create -c conda-forge -c bioconda  -c aghozlane -p /opt/software/conda_env -y dadaist2-full cutadapt=3.3 qax  r-gunifrac 

%runscript
    exec dadaist2 "$@"
```

and the image built as:

```bash
sudo singularity build dadaist2.simg dadaist2-stable.def
```

### Developmental snapshot via Singularity

To get the latest code from the repository (and also some databases) here we share a `dadaist2-dev.def` file:

```singularity
Bootstrap: docker
From: centos:centos7.6.1810

%environment
    source /opt/software/conda/bin/activate /opt/software/conda_env
    export PATH=/opt/software/dadaist2/bin:$PATH
    PERL5LIB=''
    LANG=C
         
%post
    yum -y install epel-release wget which nano curl zlib-devel git free
    yum -y groupinstall "Development Tools"
    mkdir -p /dadaist_databases/
    mkdir -p /opt/software
    cd /opt/software  
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
    sh ./Miniconda3-latest-Linux-x86_64.sh -p /opt/software/conda -b
    /opt/software/conda/bin/conda config --add channels defaults
    /opt/software/conda/bin/conda config --add channels conda-forge
    /opt/software/conda/bin/conda config --add channels bioconda
    /opt/software/conda/bin/conda install -y mamba
    /opt/software/conda/bin/mamba create  -c aghozlane -p /opt/software/conda_env -y dadaist2-full cutadapt=3.3 qax  r-gunifrac 
    source /opt/software/conda/bin/activate /opt/software/conda_env
    cd /opt/software
    git clone https://github.com/quadram-institute-bioscience/dadaist2
    ./dadaist2/bin/dadaist2-getdb -d "dada2-unite" -o /dadaist_databases/
    ./dadaist2/bin/dadaist2-getdb -d "decipher-silva-138" -o /dadaist_databases/

%runscript
    exec dadaist2 "$@"
```
and the image built as:
```bash
sudo singularity build dadaist2-dev.simg dadaist2-dev.def
```

### Docker file 
To build an image with the latest version from Bioconda the following definition file can be saved as `Dockerfile`:
```
Bootstrap: docker
From: centos:centos7.6.1810

%environment
    source /opt/software/conda/bin/activate /opt/software/conda_env
    export PATH=/opt/software/dadaist2/bin:$PATH
    PERL5LIB=''
    LANG=C
         
%post
    yum -y install epel-release wget which nano curl zlib-devel git free
    yum -y groupinstall "Development Tools"
    mkdir -p /opt/software
    cd /opt/software    
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
    sh ./Miniconda3-latest-Linux-x86_64.sh -p /opt/software/conda -b
    /opt/software/conda/bin/conda install -y mamba
    /opt/software/conda/bin/mamba create -c conda-forge -c bioconda  -c aghozlane -p /opt/software/conda_env -y dadaist2-full cutadapt=3.3 qax  r-gunifrac 

%runscript
    exec dadaist2 "$@"
```

and the image built as:
```bash
sudo singularity build dadaist2.simg dadaist2-stable.def
```

### Developmental snapshot with Docker

To get the latest code from the repository (and also some databases) here we share a `dadaist2-dev.def` file:
```
Bootstrap: docker
From: centos:centos7.6.1810

%environment
    source /opt/software/conda/bin/activate /opt/software/conda_env
    export PATH=/opt/software/dadaist2/bin:$PATH
    PERL5LIB=''
    LANG=C
         
%post
    yum -y install epel-release wget which nano curl zlib-devel git free
    yum -y groupinstall "Development Tools"
    mkdir -p /dadaist_databases/
    mkdir -p /opt/software
    cd /opt/software  
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
    sh ./Miniconda3-latest-Linux-x86_64.sh -p /opt/software/conda -b
    /opt/software/conda/bin/conda config --add channels defaults
    /opt/software/conda/bin/conda config --add channels conda-forge
    /opt/software/conda/bin/conda config --add channels bioconda
    /opt/software/conda/bin/conda install -y mamba
    /opt/software/conda/bin/mamba create  -c aghozlane -p /opt/software/conda_env -y dadaist2-full cutadapt=3.3 qax  r-gunifrac 
    source /opt/software/conda/bin/activate /opt/software/conda_env
    cd /opt/software
    git clone https://github.com/quadram-institute-bioscience/dadaist2
    ./dadaist2/bin/dadaist2-getdb -d "dada2-unite" -o /dadaist_databases/
    ./dadaist2/bin/dadaist2-getdb -d "decipher-silva-138" -o /dadaist_databases/

%runscript
    exec dadaist2 "$@"
```
and the image built as:
```bash
sudo singularity build dadaist2-dev.simg dadaist2-dev.def
```

## Provided scripts

**dadaist2** will install the following programs:

* `dadaist2`, the main program
* `dadaist2-...`, a set of wrappers and tools all using the _dadaist2-_ prefix to make them easy to find (using [TAB](https://www.howtogeek.com/195207/use-tab-completion-to-type-commands-faster-on-any-operating-system/))
<a href="https://quadram-institute-bioscience.github.io/dadaist2" description="Dadaist documentation">
<img src="img/dadaist.png">
</a>
<br/>

# Dadaist2: highway to R

:package: See the **[repository](https://github.com/quadram-institute-bioscience/dadaist2)**

[![release](https://img.shields.io/github/v/release/quadram-institute-bioscience/dadaist2?label=github%20release)](https://github.com/quadram-institute-bioscience/dadaist2/releases)
[![version](https://img.shields.io/conda/v/bioconda/dadaist2?label=bioconda)](https://bioconda.github.io/recipes/dadaist2/README.html)
[![Conda](https://img.shields.io/conda/dn/bioconda/dadaist2)](https://bioconda.github.io/recipes/dadaist2/README.html)
[![Build Status](https://www.travis-ci.com/quadram-institute-bioscience/dadaist2.svg?branch=master)](https://www.travis-ci.com/quadram-institute-bioscience/dadaist2)
[![Dadaist-CI](https://github.com/quadram-institute-bioscience/dadaist2/actions/workflows/blank.yml/badge.svg)](https://github.com/quadram-institute-bioscience/dadaist2/actions/workflows/blank.yml)


Standalone wrapper for [DADA2](https://benjjneb.github.io/dada2/index.html) package, to quickly generate a feature table and a
set of representative sequences from a folder with Paired End Illumina reads.
*Dadaist2* is designed to simplify the stream of data from the read processing to the statistical analysis and plots.

Developeb by Andrea Telatin, Rebecca Ansorge and Giovanni Birolo. 

Dadaist2 is a highway to downstream analyses:
* Generation of a [*PhyloSeq*](https://joey711.github.io/phyloseq/) object, for immediate usage in R
* Possibility to run in the pipeline a _custom R script_ that starts from the PhyloSeq object
* Generation of [*MicrobiomeAnalyst*](https://www.microbiomeanalyst.ca)-compatible files. MicrobiomeAnalyst provides a web-interface to performgi a broad range of visualizations and analyses.
* Generation of [*Rhea*](https://lagkouvardos.github.io/Rhea/)-compatible files. Rhea is a standardized set of scripts "_designed to help easy implementation by users_".

In addition to this, Dadaist:
* Can automatically detect quality boundaries or trim the primers
* Has a custom mode for _variable length_ amplicons (i.e. ITS), to detect features longer than the sum of the paired-end reads.
* Ships an open source implementation of *[UNCROSS2](https://www.biorxiv.org/content/10.1101/400762v1.full)* by Robert Edgar.
* Has a modular design that allows recycling parts of it in custom workflows.
* Prepares [a MultiQC-enabled](https://quadram-institute-bioscience.github.io/dadaist2/mqc/) overview of the experiment
* Produces an easy to inspect [HTML execution log](https://quadram-institute-bioscience.github.io/dadaist2/mqc/log.html)

You can [read more](introduction) about the features.

## The workflow

<img src="img/scheme_small.png">

## Contents

{% include list.liquid all=true %}
---
sort: 4
permalink: /usage
---
# Other tutorials

## "Mothur SOP" tutorial

The "[MiSeq SOP](https://mothur.org/wiki/miseq_sop/)" tutorial by Pat Schloss
is a well-known dataset coming from the paper:

> Kozich J. J. _et al.: (2013): **Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform**. Applied and Environmental Microbiology [link](https://aem.asm.org/content/79/17/5112).

We will use this dataset to show how easy is to reproduce the main results with dadaist.

### Get the dataset

The dataset consits of 20 paired-end FASTQ files coming from the full dataset of the paper.

```bash
wget "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
unzip miseqsopdata.zip
```

We can use SeqFu to check the number of reads:
```bash
seqfu stats -n MiSeq_SOP/*R1*.fastq
```

This will produce the following table:
```
┌─────────────────────────────────────────┬───────┬──────────┬───────┬─────┬─────┬─────┬─────┬─────┐
│ File                                    │ #Seq  │ Total bp │ Avg   │ N50 │ N75 │ N90 │ Min │ Max │
├─────────────────────────────────────────┼───────┼──────────┼───────┼─────┼─────┼─────┼─────┼─────┤
│ MiSeq_SOP/F3D0_S188_L001_R1_001.fastq   │ 7793  │ 1955914  │ 251.0 │ 251 │ 251 │ 251 │ 249 │ 251 │
│ MiSeq_SOP/F3D141_S207_L001_R1_001.fastq │ 5958  │ 1495347  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D142_S208_L001_R1_001.fastq │ 3183  │ 798912   │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D143_S209_L001_R1_001.fastq │ 3178  │ 797624   │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D144_S210_L001_R1_001.fastq │ 4827  │ 1211542  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D145_S211_L001_R1_001.fastq │ 7377  │ 1851563  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D146_S212_L001_R1_001.fastq │ 5021  │ 1260202  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D147_S213_L001_R1_001.fastq │ 17070 │ 4284410  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D148_S214_L001_R1_001.fastq │ 12405 │ 3113525  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D149_S215_L001_R1_001.fastq │ 13083 │ 3283570  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D150_S216_L001_R1_001.fastq │ 5509  │ 1382628  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D1_S189_L001_R1_001.fastq   │ 5869  │ 1473024  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D2_S190_L001_R1_001.fastq   │ 19620 │ 4924442  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D3_S191_L001_R1_001.fastq   │ 6758  │ 1696205  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D5_S193_L001_R1_001.fastq   │ 4448  │ 1116401  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D6_S194_L001_R1_001.fastq   │ 7989  │ 2005171  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D7_S195_L001_R1_001.fastq   │ 5129  │ 1287342  │ 251.0 │ 251 │ 251 │ 251 │ 248 │ 251 │
│ MiSeq_SOP/F3D8_S196_L001_R1_001.fastq   │ 5294  │ 1328727  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/F3D9_S197_L001_R1_001.fastq   │ 7070  │ 1774487  │ 251.0 │ 251 │ 251 │ 251 │ 250 │ 251 │
│ MiSeq_SOP/Mock_S280_L001_R1_001.fastq   │ 4779  │ 1199027  │ 250.9 │ 251 │ 251 │ 250 │ 249 │ 251 │
└─────────────────────────────────────────┴───────┴──────────┴───────┴─────┴─────┴─────┴─────┴─────┘
```

### Get the reference database
If you already downloaded a suitable reference database, you can skip this step.
We will use SILVA 138 (that is available both for DECIPHER and DADA2 naive classifier).
In this example we'll use the "DADA2" version:

First, we can check which SILVA database are available for download:
```bash
dadaist2-getdb --list silva
```

Then we can download one from the list (or all of them). The following command will
download the 'dada2-silva-138' database in our home directory, under the `refs` subdirectory.

```bash
dadaist2-getdb -d dada2-silva-138 -o ~/refs/
```


### Refactoring the metadata file

The dataset comes with a tabular file with metadata called `MiSeq_SOP/mouse.time.design`.
We can rename the header and this will suffice to work with Dadaist2 (that expected
  a `#SampleID` column):

```bash
sed 's/group/#SampleID/g' MiSeq_SOP/mouse.time.design > metadata.tsv
```

### Running Dadaist2

The minimum parameter to run Dadaist are:
* `-i`: directory with the input reads, with `_R1` and `_R2` tags
* `-o`: output directory
Recommended parameters are:
* `-d`: to specify a reference database (will enable taxonomy annotation)
* `-m`: metadata file (if not provided a blank metadata file is generated)
* `-t`: number of threads
  
```bash
dadaist2  -i MiSeq_SOP/ -o dadaist2-sop -m metadata.tsv -d ~/refs/silva_nr_v138_train_set.fa.gz
```

### Exploring the dataset with Microbiome Analyst

The output folder contains a subdirectory called _MicrobiomeAnalyst_.
We can go to [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/upload/OtuUploadView.xhtml)
and upload our table, taxonomy and metadata:

<img src="img/ma-form.png" width="500" height="250"  alt="MicrobiomeAnalyst input form">

After deciding if and how to normalize or rescale our data (we used the defaults settings here),
you'll land on the tool box page:

<img src="img/ma-tools.png" width="570" height="300" alt="MicrobiomeAnalyst menu">

We can perform a PCA analysis and check if there is a separation between the "Early" and "Late"
samples, as shown in the original paper:

<img src="img/ma-pca.png" >


# Documentation

```note
This section is under costruction
```

In this website:
 * [manual pages]({{ 'pages' | relative_url }})
---
sort: 12
---
## dadaist2-phyloseqMake

**dadaist2-phyloseqMake** - Generate PhyloSeq object from the command line

## Author

Andrea Telatin and Rebecca Ansorge

## Parameters

- _-i_, _--input_ DIR

    Directory containing the `MicrobiomeAnalyst` folder generated by Dadaist2.

- _-o_, _--output_

    Output filename. If omitted, a 'phyloseq.rds' file will be placed in the input directory.

## Source code and documentation

The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 15
---
## dadaist2-crosstalk
**dadaist2-crosstalk** is an open source implementation of the
UNCROSS2 ([https://www.biorxiv.org/content/10.1101/400762v1.full.pdf](https://www.biorxiv.org/content/10.1101/400762v1.full.pdf))
algorithm by Robert Edgar.

## Usage

    usage: dadaist2-crosstalk [-h] -i INPUT [-o OUTPUT] [-v] [-d] [--version]
    
    Denoise Illumina cross-talk from OTU tables
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            OTU table filename
      -o OUTPUT, --output OUTPUT
                            Cleaned OTU table filename
      -v, --verbose         Print extra information
      -d, --debug           Print debug information
      --version             show program's version number and exit

## Rationale

Metabarcoding experiments usually result in highly multiplexed sequencing runs 
where each sample is identified by a molecular barcode that is also sequenced. 
De-multiplexing errors can cause assignment of reads to the wrong sample. 
While the effect is often negligible, it may affect analysis when a fraction 
of reads leak from a high abundance sample to a small or even zero abundance sample, 
which can affect diversity rate estimation. 

Removal of crosstalk between samples in feature tables has been already tackled 
by the UNCROSS2 algorithm, implemented in the closed source USEARCH software.

## Input and Output

Input is a feature table in TSV format having the OTU/ASV IDs in the first column
and each following column being the count of occurrences in a sample.

The output is a denoised table with the same format.
---
sort: 16
---
## dadaist2-mergeseqs
This tool merges the two paired end denoised sequences as they appear in 
a DADA2 feature table when asking DADA2 not to join the reads.

## Synopsis

     Combine pairs in DADA2 unmerged tables

     Usage: 
     dadaist2-mergeseqs [options] -i dada2.tsv 

     Options:
       -i, --input-file FILE      FASTA or FASTQ file
       -f, --fasta FILE           Write new sequences to FASTA
       -p, --pair-spacer STRING   Pairs separator [default: NNNNNNNNNN]
       -s, --strip STRING         Remove this string from sample names
       -n, --seq-name STRING      Sequence string name [default: MD5]
       -m, --max-mismatches INT   Maximum allowed mismatches [default: 0]
       --id STRING                Features column name [default: #OTU ID]
       --verbose                  Print verbose output
    

## Input

The TSV table produced by DADA2 (first column with the actual representative sequence, each following column
with the counts per sample). 
Here a truncated example:

    #OTU ID       A01_R1.fastq.gz A02_R1.fastq.gz F99_R1.fastq.gz
    GGAATTTTG[..]GGGCTTAACCTNNNNNNNNNNGCGCTTA [..]AAACAG  1263    1544    1341
    GGAATCTTC[..]TTGGCTCAACCNNNNNNNNNNAGCGCAGG[..]AAACAG  100     21      19
    GGAATTTTGG[..]CTTAACCTNNNNNNNNNNGCGCAGGCGG[..]AAACAG  490 24  296

A stretch of **Ns** separates the denoised `_R1` sequence and the denoise `_R2`.

## Output

The output is a similar table after joining the reads, when possible.
If using `--verbose` a summary will be printed to the standard error at
the end.

    Total:8;Split:8;Joined:8
---
sort: 3
---
## dadaist2-alpha
**dadaist2-normalize** - Normalize OTU table using the **Rhea** protocol.
The Rhea protocol ([https://lagkouvardos.github.io/Rhea/](https://lagkouvardos.github.io/Rhea/)) is a complete
set of scripts to analyse microbiome files. 

This wrapper is part of the _AutoRhea_ script bundled with _Dadaist2_. 
If used, please, cite the Rhea paper (see below).

## Authors
Andrea Telatin and Rebecca Ansorge

## Usage
    dadaist2-normalize [options] -i TABLE -o OUTDIR 

- _-i_, _--input-table_ FILE

    Input file, the **normalized** OTU table

- _-o_, _--output-outdir_ DIR

    Output directory

- _-e_, _--effective-richness_ FLOAT

    Effective richness (default: 0.0025)

- _-s_, _--standard-richness_ INT

    Standard richness (default: 1000)

## Output files
- _alpha-diversity.tab_

    Table with: Richness, Shannon.Index, Shannon.Effective, Simpson.Index, Simpson.Effective, and Evenness
    for each sample

## Citation
If you use **Rhea** in your work please cite/attribute the original publication:

     Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) 
     Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. 
     PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
    

## Source code and documentation
This wrapper is part of **Dadaist2** freely available at 
[https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 14
---
## dadaist2-taxplot
**dadaist2-taxaplot** - Automatically plot taxonomy barbplots
from a PhyloSeq object

## Authors
Andrea Telatin and Rebecca Ansorge

## Usage
    dadaist2-taxaplot [options] -i PHYLOSEQFILE -o OUTDIR 

- _-i_, _--input-file_ FILE

    Input file, the **normalized** OTU table

- _-o_, _--output-dir_ DIR

    Output directory

## Output files
- _abundance\_bar\_plots.pdf_

    Stacked barchart taxonomy plot.

- _bubble\_plots.pdf_

    Bubble plot (PDF) with the top taxa found in the analysed samples.

## Source code and documentation
This wrapper is part of **Dadaist2** freely available at 
[https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 13
---
## dadaist2-taxonomy-binning
**dadaist2-taxonomy-binning** - Normalize OTU table using the **Rhea** protocol.
The Rhea protocol ([https://lagkouvardos.github.io/Rhea/](https://lagkouvardos.github.io/Rhea/)) is a complete
set of scripts to analyse microbiome files. 

This wrapper is part of the _AutoRhea_ script bundled with _Dadaist2_. 
If used, please, cite the Rhea paper (see below).

## Authors
Andrea Telatin and Rebecca Ansorge

## Usage
    dadaist2-taxonomy-binning -i TABLE -o OUTDIR 

- _-i_, _--input-table_ FILE

    Input OTU table, the **normalized** and with **taxonomy column**.
    The default name is `OTUs_Table-norm-rel-tax.tab`.

- _-o_, _--output-outdir_ DIR

    Output directory. 

## Output files
- _0.Kingdom.all.tab_

    Relative abundances a the Kingdom level

- _1.Phyla.all.tab_

    Relative abundances a the Phylum level

- _2.Classes.all.tab_

    Relative abundances a the Class level

- _3.Orders.all.tab_

    Relative abundances a the Order level

- _4.Families.all.tab_

    Relative abundances a the Family level

- _5.Genera.all.tab_

    Relative abundances a the Genus level

- _tax.summary.all.tab_

    Summary table (all ranks)

- _taxonomic-overview.pdf_

    Stacked bar plots in PDF format

## Citation
If you use **Rhea** in your work please cite/attribute the original publication:

     Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) 
     Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. 
     PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
    

## Source code and documentation
This wrapper is part of **Dadaist2** freely available at 
[https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 11
---
## dadaist2-phyloseqCheck
**dadaist2-phyloseqCheck** - Check PhyloSeq object from the command line

## Author
Andrea Telatin and Rebecca Ansorge

## Parameters
- _-i_, _--input_ FILE

    Input file in in PhyloSeq object (R Object)

- _-j_, _--json_

    Output in JSON format

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 9
---
## dadaist2-metadata
**dadaist2-metadata** - create a sample sheet from a list of Paired End FASTQ files,
that can be used as a template to add further columns.
This is automatically called by `dadaist2`, but it can be used to generate a valid
templeate to be extended with more columns.

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Synopsis
    makeSampleSheet [options] -i INPUT_DIR

## Parameters
- _-i_, _--input-directory_ DIRECTORY

    Directory containing the paired end files in FASTQ format, gzipped or not.

- _-o_, _--output-file_ FILE

    Output file, if not specified will be printed to STDOUT.

- _-1_, _--for-tag_ (and _-2_, _--rev-tag_) TAG

    Identifier for forward and reverse reads (default: \_R1 and \_R2)

- _-s_, _id-separator_ STRING

    Sample name separator (default: \_)

- _-f_, _--field-separator_ CHAR

    Separator in the output file table (default: \\t)

- _-h_, _--header-first-col_ COLNAME

    dadaist2-metadata of the first column header (default: #SampleID)

- _--add-full-path_

    Add a colum with the absolute path of the sample Reads

- _--add-mock-column_ COLNAME

    Add an extra column named `COLNAME` having as value what is specified by
    `--mock-value`

- _---mock-value_ VALUE

    Default value used to fill an optional column (requires `--add-mock-column`). Default "sample".

- _--version_

    Print version and exit.

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 2
---
## dadaist2-addTaxToFasta
**dadaist2-addTaxToFasta** - Add taxonomy annotation to the FASTA file with
the representative sequences

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Usage
    dadaist2-assigntax -i FASTA -o DIR -r REFERENCE [-t THREADS]

- _-f_, _--fasta_ FASTA

    Input file in FASTA format (or in DADA2 table format)

- _-o_, _--output_ FASTA

    Output file in FASTA format. If not provided will be printed to the standard output.

- _-t_, _--taxonomy_ FILE

    "taxonomy.tsv" file as produced by `dadaist2-assigntax`

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 17
---

# dadaist2-rundada

```note
This is a new wrapper, introduced in 1.2.0, and experimental
```

Wrapper for DADA2 without any modification of the input reads.
The input can be supplied either:

* As a single directory containing the reads (`-i DIRECTORY`), or
* As two separate directory, one for the forward reads (`-f FOR_DIR`) and one for the reverse reads (`-r REV_DIR`).

The latter is used as a compatibility layer and will be used by _dadaist2_ itself to invoke the wrapper.

## Synopsis

```text
usage: dadaist2-rundada [-h] [-i INPUT_DIR] [-f FOR_DIR] [-r REV_DIR] -o OUTPUT_DIR [--tmp TMP] [--fortag FORTAG] [--revertag REVERTAG] [--sample-separator SAMPLE_SEPARATOR]
                        [--sample-extension SAMPLE_EXTENSION] [-q TRUNC_QUAL] [-j] [-p] [--trunc-len-1 TRUNC_LEN_1] [--trunc-len-2 TRUNC_LEN_2] [--trim-left-1 TRIM_LEFT_1] [--trim-left-2 TRIM_LEFT_2]
                        [--maxee-1 MAXEE_1] [--maxee-2 MAXEE_2] [--chimera {none,pooled,consensus}] [--min-parent-fold MIN_PARENT_FOLD] [--n-learn N_LEARN] [-t THREADS] [--keep-temp] [--save-rds]
                        [--save-plots] [--log LOG] [--copy] [--skip-checks] [--verbose]

Run DADA2

optional arguments:
  -h, --help            show this help message and exit

Main:
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Input directory with both R1 and R2
  -f FOR_DIR, --for-dir FOR_DIR
                        Input directory with R1 reads
  -r REV_DIR, --rev-dir REV_DIR
                        Input directory with R2 reads
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory
  --tmp TMP             Temporary directory

Input filtering:
  --fortag FORTAG       String defining a file as forward [default: _R1]
  --revertag REVERTAG   String defining a file as reverse [default: _R2]
  --sample-separator SAMPLE_SEPARATOR
                        String acting as samplename separator [default: _]
  --sample-extension SAMPLE_EXTENSION
                        String acting as samplename extension [default: .fastq.gz]

DADA2 parameters:
  -q TRUNC_QUAL, --trunc-qual TRUNC_QUAL
                        Truncate at the first occurrence of a base with Q lower [default: 8]
  -j, --join            Join without merging
  -p, --pool            Pool samples
  --trunc-len-1 TRUNC_LEN_1
                        Position at which to truncate forward reads [default: 0]
  --trunc-len-2 TRUNC_LEN_2
                        Position at which to truncate reverse reads [default: 0]
  --trim-left-1 TRIM_LEFT_1
                        Number of nucleotide to trim from the beginning of forward reads [default: 0]
  --trim-left-2 TRIM_LEFT_2
                        Number of nucleotide to trim from the beginning of reverse reads [default: 0]
  --maxee-1 MAXEE_1     Maximum expected errors in forward reads [default: 1.0]
  --maxee-2 MAXEE_2     Maximum expected errors in reverse reads [default: 1.00]
  --chimera {none,pooled,consensus}
                        Chimera handling can be none, pooled or consensus [default: pooled]
  --min-parent-fold MIN_PARENT_FOLD
                        Minimum abundance of parents of a chimeric sequence (>1.0) [default: 1.0]
  --n-learn N_LEARN     Number of reads to learn the model, 0 for all [default: 0]

Other parameters:
  -t THREADS, --threads THREADS
                        Number of threads
  --keep-temp           Keep temporary files
  --save-rds            Save RDS file with DADA2 output
  --save-plots          Save Quality plots of the input reads (PDF)
  --log LOG             Log file
  --copy                Copy input files instead of symbolic linking
  --skip-checks         Do not check installation of dependencies
  --verbose             Verbose mode
```

## Output files

The output directory will contain:

* dada2.stats (table with the statistics of reads loss)
* dada2.tsv (main feature table)
* dada2.rds (R object with the table. if `--save-rds` is specified)
* quality_R1.pdf, quality_R2.pdf (quality plots, if `--save-plots` is specified)
* dada2.execution.log, dada2.execution.txt (wrapper log files)
---
sort: 6
---
## dadaist2-exporter
**dadaist2-exporter** - tool to export dadaist2 output into MicrobiomeAnalyst
compatible format. _MicrobiomeAnalyst_ can be used as an **R** module or
via the user-friendly website [https://www.microbiomeanalyst.ca/](https://www.microbiomeanalyst.ca/) and
_Rhea_ ([https://lagkouvardos.github.io/Rhea/](https://lagkouvardos.github.io/Rhea/)).

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Synopsis
dadaist2-exporter \[options\] -i INPUT\_DIR

## Parameters
- _-i_, _--input-directory_ DIRECTORY

    Directory containing the paired end files in FASTQ format, gzipped or not.

- _-o_, _--output-directory_ DIRECTORY

    Output directory, by default will be a subdirectory called `MicrobiomeAnalyst`
    inside the input directory.

- _--skip-rhea_

    Do not create the **Rhea** subdirectory and its files.

- _--skip-ma_

    Do not create the **MicrobiomeAnalyst** subdirectory and its files.

- _--version_

    Print version and exit.

## Output
The output directory will contain:

- _metadata.csv_

    Metadata file to be used in the omonymous field.

- _table.csv_

    Feature table to be used in the 'OTU/ASV table' field.

- _taxonomy.csv_

    Taxonomy table to be used in the 'Taxonomy table' field.

- _seqs.fa_

    Not used in MicrobiomeAnalyst, but kept for reference.

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 5
---

# Help pages

This section contains the help pages of the tools provided by **dadaist2** in alphabetical order.
All the wrappers documentation is updated at each commit, while other tools are update at each release.

{% include list.liquid %}
---
sort: 8
---
## dadaist2-importq2
**dadaist2-importq2** - create a PhyloSeq object from a set of
Qiime2 artifacts.

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Synopsis
    dadaist2-importq2 [options] 

## Parameters
- _-t_, _--feature-table_ ARTIFACT

    The feature table (e.g. from DADA2)

- _-m_, _--metadata-file_ FILE

    The metadata file used by Qiime2

- _-e_ _--tree_ ARTIFACT

    Rooted tree artifact.

- _-x_ _--taxonomy_ ARTIFACT

    Taxonomy table artifact.

- _-r_, _--rep-seqs_ ARTIFACT

    Representative sequences (e.g. from DADA2)

- _-o_, _--output-phyloseq_ FILE

    The filename for the PhyloSeq object to produce (default: phyloseq.rds)

- _--version_

    Print version and exit.

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 7
---
## dadaist2-getdb
**dadaist2-getdb** - download reference databases for dadaist2

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## List available databases
    dadaist2-getdb --list [query]

If a `query` keyword is specified, only matching entries will be printed.

## Download one or more databases
    dadaist2-getdb -d DATABASE_NAME [-o OUTPUT_DIR]

    dadaist2-getdb -d DB1 -d DB2 -d DB3 [-o OUTPUT_DIR]

    dadaist2-getdb -q QUERY_STRING

- _-d_, _--database_ ID

    Identifier of the database to be downloaded (list available database and their
    identifier name using `dadaist2-getdb --list`). This parameter can be repeated
    multiple times to download multiple databases.

- _-q_, _--query_ STRING

    Download all databases matching the query string ('.' for all)

- _-o_, _--output-dir_ DIR

    Output directory, or the current working directory if not specified.

- _-t_, _--temp-dir_ DIR

    Temporary directory (default: `$TMPDIR`).

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 5
---
## dadaist2-dada2fasta
**dadaist2-dada2fasta** - a program to process the feature table generated by DADA2
(that uses the sequences as feature names) and saves it as feature table (using
progressive feature names, or the MD5 of the sequences).

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Synopsis
    dadaist2-dada2fasta  -i dada2table.tsv -o table.tsv -r repseqs.fasta

## Parameters
## Main Parameters

- _-i_, _--input_ FILE

    Output produced by DADA2 (feature table tsv).

- _-o_, _--output-table_ FILE

    Output feature table

- _-r_, _--rep-seqs_ FILE

    Fasta output with the representative sequences.

- _-s_, _--strip-pattern_ STR

    Remove from the sample names this string, usually found as
    filename suffix (default: \_R1.fastq.gz)

- _-p_, _--otu-prefix_ STR

    Prefix used for the represenative sequences, by default the
    MD5 of the sequence (example: ASV)

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains the full _DOCUMENTATION_ and we recommend 
checking for updates.
---
sort: 1
---
## dadaist2
**dadaist2** - a shell wrapper for DADA2, to detect representative sequences and
generate a feature table starting from Illumina Paired End reads. 
This is the main program of the _dadaist2 toolkit_ that includes several
wrappers and utilities to streamline the analysis
of metabarcoding reads from the Linux shell to R.

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Synopsis
    dadaist2 [options] -i INPUT_DIR -o OUTPUT_DIR

## Parameters
## Main Parameters

- _-i_, _--input-directory_ DIRECTORY

    Directory containing the paired end files in FASTQ format, gzipped or not.

- _-o_, _--output-directory_ DIRECTORY

    Output directory (will be created). 
    It is recommended recomment using a new directory for each run.

- _-d_, _--database_ DATABASE

    Reference database in gzipped FASTA format. 
    Optional (default: skip) but highly recommended.

- _-m_, _--metadata_ FILE

    Metadata file in TSV format, first column must match sample IDs. 
    If not supplied a template will be autogenerated using `dadaist2-metadata`.

- _-t_, _--threads_ INT

    Number of threads (default: 2)

- _--primers_ FOR:REV

    Strip primers with _cutadapt_, supply both sequences separated by a colon.

- _-j_, _--just-concat_

    Do not try merging paired end reads, just concatenate. 

- _--fastp_

    Perform the legacy "fastp" QC

- _--no-trim_

    Do not trim primers (using fastp). 
    Equivalent to `--trim-primer-for 0 --trim-primer-rev 0`.

- _--dada-pool_

    Pool samples in DADA2 analysis (experimental)

## Input reads

We recommend to prepare a polished directory of input reads having
filenames like Samplename\_R1.fastq.gz and Samplename\_R2.fastq.gz.

Filename starting by numbers are not accepted.

- _-1_, _--for-tag_ TAG

    Tag to recognize forward reads (default: \_R1)

- _-2_, _--rev-tag_ TAG

    Tag to recognize reverse reads (default: \_R2)

- _-s_, _--id-separator_ CHAR

    Character separating the sample name from the rest of the filename (default: \_)

## Metabarcoding processing

- _--trunc-len-1_ and _--trunc-len-2_ INT

    Position at which truncate reads (forward and reverse, respectively).

- _-q_, _--min-qual_ FLOAT

    Minimum average quality for DADA2 truncation (default: 28)

- _--no-trunc_

    Do not truncate reads at the end
     (required for non-overlapping amplicons, like ITS)

- _--maxee1_, and _--maxee2_ FLOAT

    Maximum Expected Errors in R1 and R2, respectively (default: 1.0 and 1.5)

- _-s1_, _--trim-primer-for_ INT

    Trim primer from R1 read specifying the number of bases. Similarly
    use `-s2` (`--trim-primer-rev`) to remove the front bases from the
    reverse pair (R2). Default: 20 bases each side.

- _--save-rds_

    Save a copy of the RDS file (default: off)

- _--max-loss_ FLOAT

    After DADA2 run, check the amount of reads globally remaining from input
    to non-chimeric, abort if the ratio is below threshold (default: 0.2)

## Other parameters

- _--crosstalk_ 

    Remove crosstalk using the UNCROSS2 algorithm 
    as described here [https://doi.org/10.1101/400762](https://doi.org/10.1101/400762).

- _-p_, _--prefix_ STRING

    Prefix for the output FASTA file, if "MD5" is specified, the sequence MD5 hash
    will be used instead. Default is "ASV".

- _-l_, _--log-file_ FILE

    Filename for the program log.

- _--tmp-dir_ DIR

    Where to place the temporary directory (default are system temp dir
    or `$TMPDIR`).

- _--skip-tree_

    Do not generate tree. _Experimental|Not recommended_.

- _--skip-plots_

    Do not generate quality plots.

- _--popup_

    Display popup notifications (tested on MacOS and Ubuntu)

- _--quiet_

    Reduce verbosity

- _--verbose_ and _--debug_

    Increase reported information

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains the full _DOCUMENTATION_ and we recommend 
checking for updates.

The paper describing Dadaist2 was published in:

Ansorge, R.; Birolo, G.; James, S.A.; Telatin, A. _Dadaist2: A Toolkit to Automate and Simplify Statistical Analysis and Plotting of Metabarcoding Experiments_.
Int. J. Mol. Sci. 2021, 22, 5309. [https://doi.org/10.3390/ijms22105309](https://doi.org/10.3390/ijms22105309)
---
sort: 10
---
## dadaist2-normalize
**dadaist2-normalize** - Normalize OTU table using the **Rhea** protocol.
The Rhea protocol ([https://lagkouvardos.github.io/Rhea/](https://lagkouvardos.github.io/Rhea/)) is a complete
set of scripts to analyse microbiome files.

This wrapper is part of the _AutoRhea_ script bundled with _Dadaist2_.
If used, please, cite the Rhea paper (see below).

## Authors
Andrea Telatin and Rebecca Ansorge

## Usage
    dadaist2-normalize [options] -i TABLE -o OUTDIR

- _-i_, _--input-table_ FILE

    Input file in in PhyloSeq object (R Object)

- _-o_, _--output-outdir_ DIR

    Output directory

- _-r_, _--random-subsampling_

    Use random subsampling (default: off)

- _-f_, _fixed-value_

    Normalized using a fixed value (default: minimum)

- _-c_, _--cutoff_ INT

    Normalization cutoff (if _--fixed-value_ is used)

- _-n_, _--n-labels_ INT

    Highlight the INT  most undersampled samples

## Citation
If you use **Rhea** in your work please cite/attribute the original publication:

    Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017)
    Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons.
    PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836

## Source code and documentation
This wrapper is part of **Dadaist2** freely available at
[https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 4
---
## dadaist2-assigntax
**dadaist2-assigntax** - Assign Taxonomy

## Author
Andrea Telatin <andrea.telatin@quadram.ac.uk>

## Usage
    dadaist2-assigntax [options] -i FASTA -o DIR -r REFERENCE 

- _-m_, _--method_

    Taxonomy assignment method, either DECIPHER or DADA2
    (default: DECIPHER)

- _-i_, _--input_ FASTA

    Input file in FASTA format (or in DADA2 table format)

- _-o_, _--outdir_ DIR

    Output directory, or the current working directory if not specified.

- _-f_, _--fasta_ FILE

    Save taxonomy assigned FASTA file.

- _-r_, _--reference_ FILE

    RData file with the training set in DECIPHER format.

- _-t_, _--threads_ INT

    Number of threads to use.

- _-u_, _--underscore-join_

    Join taxa names that have spaces with an underscore (default:
    use double quotes)

## Source code and documentation
The program is freely available at [https://quadram-institute-bioscience.github.io/dadaist2](https://quadram-institute-bioscience.github.io/dadaist2)
released under the MIT licence. The website contains further DOCUMENTATION.
---
sort: 13
---
## dadaist2-mqc-report

**dadaist2-mqc-report** generates a MultiQC-ready folder starting from the data
available in a _Dadaist2_ run where taxonomy was assigned.

An [example report](https://quadram-institute-bioscience.github.io/dadaist2/mqc/) is
available.

## Synopsis

```
usage: dadaist2-mqc-report [-h] -i INPUT_DIR [-t TOPTAXA] -o OUTPUT_DIR

Produce multiqc report

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input-dir INPUT_DIR
  -t TOPTAXA, --toptaxa TOPTAXA
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
```

## Description

Beside checking the [execution logs](https://quadram-institute-bioscience.github.io/dadaist2/mqc/log.html), it's useful to have an overview of the whole experiment before performing
accurate analyses.

### Content of the report

* DADA2 denoising statistics: how many input reads, filtered reads, denoised reads, merged reads and final reads after chimera removal
* [Octave plots](https://www.biorxiv.org/content/10.1101/389833v1) 
to evaluate the distribution of the abundance counts
* The most abundant FASTA sequences, both classified and unclassified
* Taxonomic barplots at different levels


## Usage

MultiQC is required, an it's automatically installed with the `dadaist2-full` package
but not with the thinner `dadaist2`. See installation.

1. Generate the report files: `dadaist2-mqc-report -i DadaistOutput -o DadaistOutput/qc/`
2. Generate the report: `multiqc -c DadaistOutput/qc/config.yaml DadaistOutput/qc/`
3. Open the report with a browser (will be `multiqc_report.html` in the output folder)
---
sort: 6
---

# Integration with R

Dadaist2 relies on automation of R processes at the core, so as a design choice the wrapper and launchers all collect data to be fed to R Scripts. 
In this way it’s impossible to distinguish the output produced natively by R or from Dadaist scripts. 
For example the `dadaist2-assigntax` program contains the logic and the input collection/validation, while it’s companion script `D2-dada-taxonomy.R` is executed to perform the Dada2/Decipher taxonomy annotation. 
Moreover, by design, we promote the installation via Miniconda to ensure the creation of an environment 
with a dedicated R ecosystem with version control of R and its dependencies. 
The importance of this decision is also to make the current and future updates on the pipeline -- 
that we are committed to produce as part of our ongoing longitudinal studies -- 
will remain independent of the actual implementation as long as the R part is on a separate R script.
The comparison we made on the primary analysis is very important, 
as the primary analysis is a mix of R (DADA2) and external tools 
(quality check, adaptor removal etc) and the code is more complex and benefits 
from a functional comparison. We initially omitted the comparison of secondary 
analyses as they only require R commands (and thus have a complete overlap 
between Dadaist and the same commands run in R) and for the wider variability 
in secondary analysis that the user is enabled by Dadaist.

## Comparison of an R-based workflow

We performed a manual analysis in R requiring the generation of a PhyloSeq object, followed
by the generation of a taxonomy plot (using PhyloSeq) and a bubble plot with custom
R commands.

The manual R workflow [is available here](plot.html), while the Dadaist2 script used is `dadaist2-taxplot`.

## Results

The two workflow produce identical PhyloSeq objects and identical plots.


<img src="../img/r-taxa.png">

Taxonomy barplots (manually generated from R, _left_, and automatically generated, _right_)


<img src="../img/r-bubble.png">

Taxonomy bubble plots (manually generated from R, _left_, and automatically generated, _right_)

---
sort: 2
---

# ITS analysis

Dadaist2 has been designed to fully support variable lenght amplicons, including fungal ITS.

This is possible with:
* a primer removal tool that will detect and discard concatamers (fu-primers)
* the possibility to skip pair-end merging (with `--just-concat`, or `-j` for short) and to re-join when possible using `dadaist2-mergeseqs`, that takes DADA2 feature table as input.
* support for taxonomy assignment in non contiguous sequences

## Why are ITS amplicons different?

They are not.

**16S rDNA** is the most common target for metabarcoding, and while there are several possible primer-pairs available, the commonly accepted 
_rule of the thumb_ is to choose a primer pair that will produce and amplicon shorter than the length of the two paired end sequences
(if adopting paired-end), to ensure that the two fragments will be merged. This is possible because 16S rDNA is highly conserved in length
and the majority of the variation observed in the commonly sequenced areas is mostly due to sequence alterations, not insertions and deletions.

**ITS amplicons** (and other targets) have a wider variability in length. This has a very annoying effect immediately when preparing the libraries:
the shorter targets will be amplified with higher efficiency, making the amplification bias worse than it already is with less variable targets.

In short: with "ITS protocol" here we refer to "amplicons highly variable in length".

## How to analyse ITS amplicons

*TLDR*: simply add the `-j` flag (or `--just-concat`) when running the analysis.

Dadaist2 provides support for the denoising mode of DADA2 where the two pairs are not merged. This mode is not supported - for example - in the
Qiime2 plugin.

DADA2 itself can either merge or [not](https://github.com/benjjneb/dada2/issues/279), while Dadaist2 will merge the representative sequences
that overlap, leaving unmerged those which don't.

## Increased sensitivity: a simulation

We downloaded the [UNITE database](https://unite.ut.ee/repository.php) (95,481 sequence) and performed
an _in silico_ PCR using ITS1 and ITS4 primers as follows:
```
seqkit amplicon -F CTTGGTCATTTAGAGGAAGTAA -R GCTGCGTTCTTCATCGATGC unite.fasta > unite-amplicons.fa
```

this generated  2,833 sequences (of which 2,629 unique), with the following metrics:

```text
seqfu stats --nice --basename unite.fa unite-amplicons.fa
┌─────────────────┬────────┬────────────┬───────┬─────┬─────┬─────┬────────┬─────┬──────┐
│ File            │ #Seq   │ Total bp   │ Avg   │ N50 │ N75 │ N90 │ auN    │ Min │ Max  │
├─────────────────┼────────┼────────────┼───────┼─────┼─────┼─────┼────────┼─────┼──────┤
│ unite           │ 95,481 │ 59,469,042 │ 622.8 │ 595 │ 522 │ 479 │ 73.827 │ 141 │ 7491 │
│ unite-amplicons │ 2,833  │ 842,403    │ 297.4 │ 295 │ 266 │ 235 │ 63.436 │ 171 │ 1159 │
└─────────────────┴────────┴────────────┴───────┴─────┴─────┴─────┴────────┴─────┴──────┘
```

As we can see, there are sequences larger than 600 bp, so without any overlap in MiSeq 2x300 sequencing.
Considering at least 20 bases of overlap, there are 53 unique sequences longer than 550 bases.

```bash
seqfu derep unite-amplicons.fa --min-length 580  | seqfu stats
```

### Results

The resulting file was re-classified using DECIPHER (via `dadaist2-assigntax`) to check how the sequences would be classified if provided in full (only unique classifiactions are kept here):
```text
Fungi Ascomycota Dothideomycetes Asterinales Asterinaceae Blastacervulus
Fungi Ascomycota Dothideomycetes Capnodiales Mycosphaerellaceae Scleroramularia
Fungi Ascomycota Eurotiomycetes Chaetothyriales Herpotrichiellaceae Cladophialophora
Fungi Ascomycota Eurotiomycetes Chaetothyriales Herpotrichiellaceae Exophiala
Fungi Ascomycota Eurotiomycetes Chaetothyriales unidentified_434 unidentified_304
Fungi Ascomycota Eurotiomycetes Eurotiales Trichocomaceae Rasamsonia
Fungi Ascomycota Eurotiomycetes Phaeomoniellales Phaeomoniellaceae Pseudophaeomoniella
Fungi Ascomycota Geoglossomycetes Geoglossales Geoglossaceae Geoglossum
Fungi Ascomycota Lecanoromycetes Lecanorales Cladoniaceae Cladonia
Fungi Ascomycota Leotiomycetes Helotiales Chrysodiscaceae Chrysodisca
Fungi Ascomycota Leotiomycetes Helotiales Dermateaceae Calloria
Fungi Ascomycota Leotiomycetes Helotiales Dermateaceae Parafabraea
Fungi Ascomycota Leotiomycetes Helotiales Dermateaceae Pseudofabraea
Fungi Ascomycota Leotiomycetes Helotiales Helotiaceae Hymenoscyphus
Fungi Ascomycota Leotiomycetes Helotiales Helotiaceae Phaeohelotium
Fungi Ascomycota Leotiomycetes Helotiales Helotiales_fam_Incertae_sedis Chalara
Fungi Ascomycota Leotiomycetes Helotiales Helotiales_fam_Incertae_sedis Vestigium
Fungi Ascomycota Leotiomycetes Helotiales Hyaloscyphaceae Capitotricha
Fungi Ascomycota Leotiomycetes Helotiales Hyaloscyphaceae Incrucipulum
Fungi Ascomycota Leotiomycetes Helotiales Hyaloscyphaceae Lachnum
Fungi Ascomycota Leotiomycetes Helotiales Hyaloscyphaceae unidentified_191
Fungi Ascomycota Leotiomycetes Helotiales Leotiaceae Pezoloma
Fungi Ascomycota Leotiomycetes Helotiales Myxotrichaceae Oidiodendron
Fungi Ascomycota Leotiomycetes Helotiales unidentified_8 unidentified_5
Fungi Ascomycota Orbiliomycetes Orbiliales Orbiliaceae Dactylella
Fungi Ascomycota Orbiliomycetes Orbiliales Orbiliaceae Hyalorbilia
Fungi Ascomycota Orbiliomycetes Orbiliales Orbiliaceae Orbilia
Fungi Ascomycota Orbiliomycetes Orbiliales Orbiliaceae unidentified_448
Fungi Ascomycota Pezizomycetes Pezizales Sarcosomataceae Donadinia
Fungi Ascomycota Saccharomycetes Saccharomycetales unidentified_72 unidentified_51
Fungi Ascomycota Sordariomycetes Chaetosphaeriales Chaetosphaeriaceae Chloridium
Fungi Ascomycota Sordariomycetes Xylariales Xylariaceae Annulohypoxylon
Fungi Ascomycota unidentified unidentified unidentified unidentified
Fungi Basidiomycota Agaricomycetes Hymenochaetales Hymenochaetaceae Coltricia
Fungi Basidiomycota Agaricomycetes Polyporales unidentified_544 unidentified_376
Fungi Basidiomycota Agaricomycetes Russulales Lachnocladiaceae Dichostereum
Fungi Basidiomycota Cystobasidiomycetes Erythrobasidiales Erythrobasidiaceae Bannoa
Fungi Mucoromycota Endogonomycetes GS22 unidentified_1771 unidentified_1228
```

Then the same fasta FILE has been processed to simulate a joining of paired end, so the first and last 300 bp
have been kept while replacing the middle part with Ns.

Taxonomy classification didnt change.

The files are in `./test/its` in the repository:
* `taxonomy.tsv` is the classification of unique amplicons > 550bp
* `join/taxonomy.tsv` is the classification of amplicons after joining the first and the last 300 bp
---
sort: 3
---

# Numerical ecology

The primary analysis of a metabarcoding experiment processes a set of FASTQ files (raw reads) to generate:
* a set of _representative sequences_ (either Amplicon Sequence Variants, ASVs, or Operational Taxonomic Units, OTUs)
* a _feature table_ (or _contingency table_): a matrix of counts of hits against each representative sequence per sample

Additionally:
* _taxonomic annotation_ of each representative sequence
* a _phylogenetic tree_ of the representative sequences

These files can be analysed using the principles of *numerical ecology*, to 

## MicrobiomeAnalyst

[MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/upload/OtuUploadView.xhtml)
is both an R module and a webserver to perform a range of explorative analyses and statistical tests,
like:
* Compositional profiling
* Comparative analysis
* Functional analysis
* Taxon Set Enrichment Analysis

A [nature protocol](https://www.nature.com/articles/s41596-019-0264-1) is available.

## Rhea
Dadaist2 implements the [Rhea](https://lagkouvardos.github.io/Rhea/) workflow to normalize the feature table,
analyse the alpha and beta diversity, generate taxonomy barplots.

> Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) _Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons_. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836

Dadaist2 produce a _Rhea_ subdirectory with the input files to follow the full Rhea protocol. In addition some steps (those not requiring assumptions on the experiment) are performed automatically:
* Normalization (this can be invoked independently via _dadaist2-normalize_)
* Alpha diversity (this can be invoked independently via _dadaist2-alpha_)

## PhyloSeq

[PhyloSeq](https://joey711.github.io/phyloseq/) is an R module that allows several
analyses of microbiome datasets.

Dadaist2 conveniently produces a phyloseq object that can be loaded with:

```r
ps <- loadRDS("phyloseq.rds")
```---
sort: 1
---
# General information

`dadaist2` is a wrapper designed to perform a basic [DADA2](https://benjjneb.github.io/dada2/index.html)
analysis from the command line, starting from a directory containing a set of
paired-end samples in FASTQ (eventually gzipped) format.

The wrapper is written in _Perl_, and will run _R_ scripts via `Rscript --vanilla`.

## Main dependencies

* **Perl**, with some modules: _FASTX::Reader_, _File::Temp_, and other standard modules (notably _Pod::Usage_ and _Digest::MD5_)
* **R**, and some libraries: _dada2_, _phyloseq_, and _DECIPHER_
  * DADA2 is used for denoising, feature table generation and representative sequences picking
  * DECIPHER (or alternatively DADA2) is used for taxonomy assignments
* **cutadapt** for primer removal
  * **fastp** (alternative to SeqFu QC, no trimming is performed)
* **clustalo** (for multiple sequence alignment)
* **fasttree** (to generate a tree of the representative sequences)


## Bibliography
* Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W Han, Amy Jo A Johnson, and Susan P Holmes. **Dada2: high-resolution sample inference from illumina amplicon data**. Nature methods, 13(7):581, 2016. [doi:10.1038/nmeth.3869](https://doi.org/doi:10.1038/nmeth.3869).
* Sievers F, Wilm A, Dineen D, et al. **Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega**. Molecular Systems Biology. 2011 Oct;7:539. [doi:10.1038/msb.2011.75](https://doi.org/doi:10.1038/msb.2011.75).
* Price MN, Dehal PS, Arkin AP. **FastTree 2--approximately maximum-likelihood trees for large alignments**. Plos one. 2010 Mar;5(3):e9490. [doi:10.1371/journal.pone.0009490](https://doi.org/doi:10.1371/journal.pone.0009490).
* McMurdie PJ, Holmes S. **phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data**. Plos one. 2013 ;8(4):e61217. [doi:10.1371/journal.pone.0061217](https://doi.org/doi:10.1371/journal.pone.0061217).
* Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) **Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons**. PeerJ 5:e2836 [doi:10.7717/peerj.2836](https://doi.org/10.7717/peerj.2836)
* Martin M. **Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads** [doi:10.14806/ej.17.1.200](http://journal.embnet.org/index.php/embnetjournal/article/view/200)
* 
The DADA2 wrapper was based on the R scripts from:
* [q2-dada2, a Qiime2 plugin](https://github.com/qiime2/q2-dada2)
---
sort: 5
---

# Continuous integration

Dadaist2 is a pipeline wrapping third party tools and scripts (most notably _DADA2_, _DECIPHER_, _Rhea_)
and some custom components (_e. g._ _crosstalk_, _octave plots_ ...) in a coherent framework. 

Most functions are provided as modules, and we are committed to update each module ensuring that the 
output remains reliable and in line with the original tool. To ensure that the reliability is preserved
at each release, we set up a set of tests:

* continuous integration with CircleCI tests, at each commit:
  * the generation of a feature table (and representative sequences) without a taxonomy assignment
  * the generation of a feature table (and representative sequences) *with* a taxonomy assignment
  * that the taxonomy assignment module (standalone) works correctly

## Function tests

In addition to the continuous integration, there is a more complete set of tests that is performed at each
release:

* QC
* DADA2 denoising
* Taxonomy (DADA and DECIPHER)

### QC test

The QC is available via `cutadapt` or via `fastp`. 
### DADA2 test

DADA2 provides [a tutorial](https://benjjneb.github.io/dada2/tutorial_1_8.html)
(wrote for version 1.8), based on the [MiSeq SOP](https://mothur.org/wiki/miseq_sop/) 
dataset (from [Mothur](https://mothur.org)).

In the repository we have a script (here 
[dada2-sop.R](https://github.com/quadram-institute-bioscience/dadaist2/blob/master/test/miseq-sop-compare/dada2-sop.R))
that uses the exact commands provided by the tutorial, which does not include any QC step.

### Taxonomy (DADA and DECIPHER)

## Pipeline tests

A pipeline require that the utilized components (see _Function Tests_) work together,
generating the expected output to be fed to the downstream steps.

At each release we test:

* MicrobiomeAnalyst output
* Rhea output


---


---
sort: 4
---

# Development notes

{% include list.liquid %}
---
sort: 4
---

# MultiQC report

Dadaist2 generates a MultiQC ready directory to generate a QC report using [MultiQC](https://multiqc.info),
saving the files in the `./multiqc` directory inside the output directory.

The report can be generated if MultiQC is installed, otherwise:
```
cd $OUTPUT/multiqc
multiqc -f .
```

## Notable sections

* The most abundant **representative sequences** (separated in unclassified and classified) are reported to be manually inspected, to allow to check for contaminants (among the unclassified) and to control the taxonomy classification
* **Taxonomy plots** are reported for a preliminary overview
* A correlation **matrix** allows to check if replicates are indeed closely correlated
* **Octave plots**, based on the distribution of counts per sample, allow to evaluate the level of _noisy sequences_
* Beta diversity **matrix**


## Screenshot

An [example](https://quadram-institute-bioscience.github.io/dadaist2/mqc/) is available on this website.

![popup](../img/multiqc.png)

