# Citation

This is currently unpublished work. If you choose to use it, please contact Rob Edwards (raedwards@gmail.com) or Scott Handley (shandley@wustl.edu) for an appropriate citation.
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/latest_release_date.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/platforms.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/license.svg)](https://anaconda.org/bioconda/hecatomb)
[![Documentation Status](https://readthedocs.org/projects/hecatomb/badge/?version=latest)](https://hecatomb.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/downloads.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

![](docs/img/hecatombLogo.png)

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

**For detailed pipeline overview, installation, usage and customisation instructions,
[please refer to the documentation hosted at Read the Docs](https://hecatomb.readthedocs.io).**

## Quick start guide

### Running on HPC

Hecatomb is powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and greatly benefits from the use of 
Snakemake profiles for HPC Clusters.
[More information and example for setting up Snakemake profiles for Hecatomb in the documentation](https://hecatomb.readthedocs.io/en/latest/profiles/).

### Install

```bash
# create conda env and install
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

# activate conda env
conda activate hecatomb

# check the installation
hecatomb -h

# download the databases - you only have to do this once
  # locally: uses 32 threads by default
hecatomb install

  # HPC: using a profile named 'slurm'
hecatomb install --profile slurm
```

### Run the test dataset

```bash
# locally: uses 32 threads and 64 GB RAM by default
hecatomb run --test

# HPC: using a profile named 'slurm'
hecatomb run --test --profile slurm
```

### Current limitations

Hecatomb is currently designed to only work with paired-end reads. 
We have considered making a branch for single-end reads, but that is not currently available.

When you specify a directory of reads with `--reads`, Hecatomb expects paired sequencing reads in the format 
sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

When you specify a TSV file with `--reads`, Hecatomb expects a 3-column tab separated file with the first column
specifying a sample name, and the other columns the relative or full paths to the forward and reverse read files. e.g.

```text
sample1    /path/to/reads/sample1.1.fastq.gz    /path/to/reads/sample1.2.fastq.gz
sample2    /path/to/reads/sample2.1.fastq.gz    /path/to/reads/sample2.2.fastq.gz
```

### Dependencies

The only dependency you need to get up and running with Hecatomb is [conda](https://docs.conda.io/en/latest/).
Hecatomb relies on [conda](https://docs.conda.io/en/latest/) (and [mamba](https://github.com/mamba-org/mamba))
to ensure portability and ease of installation of its dependencies.
All of Hecatomb's dependencies are installed during installation or runtime, so you don't have to worry about a thing!

### Links

[Hecatomb @ bio.tools](https://bio.tools/hecatomb)

[Hecatomb @ WorkflowHub](https://workflowhub.eu/workflows/235)

This is the current conda build recipe for Hecatomb, currently hosted at 
[https://anaconda.org/beardymcjohnface](https://anaconda.org/beardymcjohnface)

```bash
# from hecatomb/build/
conda build hecatomb/
```

## Dependencies

We make sure we're using Python 3 and we obviously need Snakemake.
We otherwise only need to specify mamba (for Snakemake to build conda envs for each rule), 
and Pysam and Plotly for a few Snakemake 'Run' rules (as you can't specify envs for these).

## Build

The dependencies are installed, and the repo is downloaded and unpacked to the conda env.
We copy the contents of `bin/` and `snakemake/` which are all that is required to run the pipeline.
The `test_data/` files are copied as they are required for running the Hecatomb tests.
We also copy `docs/` and `mkdocs.yaml` if the user wants to build a local copy of the documentation 
(they would need to separately install mkdocs).

## Requirements

Use this to build the conda env that would be generated from a bioconda install. 
Used for checking dependency versions etc for bioconda installation.
[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)                                                            
![GitHub language count](https://img.shields.io/github/languages/count/shandley/hecatomb)

# Snakemake Implementations of hecatomb

[Snakemake](https://snakemake.readthedocs.io/) is now the default implementation for Hecatomb.
Please refer to [the documentation](https://hecatomb.readthedocs.io/en/latest/) for installation and usage.
These scripts are here for historical reasons and are unlikely to be useful to most users.
By default, they are not installed with the conda build for Hecatomb.


## primer_stats.sh
A bash script that compiles statistics about the primers detected during step_01 of contaminant_removal.smk.

The script automatically calls primer_stats.R to do the final multi-join. 
Of note, you will need [tidyverse](https://www.tidyverse.org) installed for this part to work.

The script will output two *.tsv files useful for determining what primers were detected in each sample.
## Snakemake

The Hecatomb pipeline is powered by [Snakemake](https://snakemake.readthedocs.io/) 
and makes extensive use of [Conda](https://docs.conda.io/en/latest/).
It comes preconfigured and includes a launcher script to make running the pipeline simple.
We chose to use Snakemake over other workflow managers largely because of our extensive experience with it.
Snakemake does a lot of heavy lifting to make our lives easier, and the pipeline better.
It manages all the pipeline jobs, as well as generates benchmarks and reports for everything, allows the pipeline to 
be naturally reentrant, parallel, portable, robust, containerised, and many other buzzwords!

![](img/hecatombSnakemake.png)

## Overview

The pipeline consists of three main sections: Read preprocessing, Virus identification, and Assembly.

![](img/hecatombPipeline.png)

## Preprocessing

Heactomb is designed to perform rigorous quality control prior to assembly and taxonomic assignment. 
This rigor is justified primarily by the philosophy of [Garbage In, Garbage out](https://en.wikipedia.org/wiki/Garbage_in,_garbage_out). 
More specifically, the following issues are dealt with to ensure that only non-contaminat biological sequence is used for downstream analysis.

1. Non-biological sequence removal (primers, adapters)
2. Host sequence removal
3. Removal of redundant sequences (clustering)
	- Creation of sequence count table
	- Calculation of sequence properties (e.g. GC content, tetramer frequencies)

Hecatomb comes with a number of host genomes for use with this preprocessing, and bespoke host genomes can be added by the user.

## Virus identification

Hecatomb uses mmseqs2 to assign taxonomy to individual sequences. 
It does so through an iterative search approach, which is a large part of the heart of hecatomb.

The first search takes all query seqeunces (the seqtable from the preprocessing steps) and queries them against a virus protein target database. 
This initial search is a translated (nt-to-aa) search. 
All of the reads assigned a viral taxonomy in this first search need to be confirmed to be truly viral in origin. 
The issue here is that querying a target database consisting solely of viral proteins does not permit each sequence to see if it really originates from a different domain of life. 
So while a sequence may be classified as viral at this stage, once queried against a more comprehensive sequence database 
(consisting of bacteria, plants, vertebrates, fungi, etc.) those sequences may have higher statistically similar sequences in those domains of life.

It is formally possible to just do this secondary search, skipping the first step querying each sequence against a virus-only protien databse. 
However, this can be computationally very expensive. 
For example, if you have 1e6 reads as your query input (sequences in your seqtable) and you query these against the 
virus-only protien database with 1.8e6 entries that is over 25 orders of magnitude smaller than UniRef50 which has 49,410,134 entries when this was written. 
The primary search against a virus-only database serves as a 'net' to (relatively) quickly capture sequences of interest. 
In our experience, however, the secondary search against the larger database is still required as many of the sequences 
assigned a viral taxonomy when only querying the smaller virus-only sequence database will reveal themselves to be entirely something else upon this secondary search. 
Thus, the iterative approach will expedite sample processing by capturing only a small portion of the input sequences as 'potentially' viral and confirming them in the secondary search as 'likely' viral.

This process wherein a sequence is classified as viral in the first search against a virus-only database and then no 
longer classfied as viral after the second search against a more comprehensive database is where hecatomb gets it's namesake. 
Many investigators, particularly those lacking access to large computing infrastructure will rely on the results of the 
primary search which, due to the relatively small target database size, can be done on commodity hardware. 
Examining the results at this stage will regularly identify a large number of viruses. 
Basically anything with any level of sequence identity to a viral protein will be called viral even if the sequence 
originated from a bacteria, plant, fungi or other organism. 
These sequences are 'sacrificed' through the secondary search when many of the viral calls from the first search are lost.

## Assembly

The preprocessing rule also goes ahead and does assembly as contigs are an important prerequisite for many downstream analysis.
The assembly strategy consists of several steps intended to be resource efficient and to maximise representation of all present species.
First, individual [MEGAHIT](https://github.com/voutcn/megahit) assemblies are produced for each sample.
The reads are then mapped to the combined assemblies and any unmapped reads undergo another round of assembly.
All contigs are then combined and merged into a non-redundant set of contigs with [Flye](https://github.com/fenderglass/Flye).

The assembly contigs are directly annotated with MMSeqs.
The assembly is also subject to a pseudo consensus annotation approach whereby the SeqTable sequences are mapped and their
Taxonomic assignments in the BigTable are combined with the read mapping information.
We find this useful with investigating the origins of contigs of interest.
## Commands

* `hecatomb install` - Install the databases (you should only need to do this once)
* `hecatomb run` - Run the pipeline
* `hecatomb listHosts` - List the currently-available host genomes
* `hecatomb addHost` - Add your own host genome
* `hecatomb config` - Copy the default config file to the current directory (for use with `--configfile`)

## Input

Hecatomb is currently designed to only work with paired-end reads.
You can either specify a directory of reads, and Hecatomb will infer the sample names and forward/reverse files, or,
you can specify a TSV file to explicitly assign sample names and point to the corresponding read files.
In either case you just use `--reads` and Hecatomb will figure out if it's a file or directory.

When you specify a directory of reads, e.g. `hecatomb run --reads readDir/`, 
Hecatomb expects paired sequencing reads in the format sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

When you specify a TSV file, e.g. `hecatomb run --reads samples.tsv`, 
Hecatomb expects a 3-column tab separated file with the first column specifying the sample name, 
and the other columns the relative or full paths to the forward and reverse read files. e.g.

```text
sample1    /path/to/reads/sample1.1.fastq.gz    /path/to/reads/sample1.2.fastq.gz
sample2    /path/to/reads/sample2.1.fastq.gz    /path/to/reads/sample2.2.fastq.gz
```

## Read annotation + assembly

By default, Hecatomb will annotate your reads and perform an assembly.
If you have more than 32 threads available, you can increase the threads provided to the pipeline with `--threads`:

```bash
hecatomb run --reads fastq/ --threads 64
```

If you're running on a HPC cluster, you should first set up a 
[Snakemake Profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
[More info and example for Hecatomb here](configuration.md#profiles-for-hpc-clusters).
Then you would specify your profile name when running Hecatomb.
Assuming your profile is called `slurm`:

```bash
hecatomb run --reads fastq/ --profile slurm
```

Running Hecatomb on a HPC with a Snakemake profile is THE BEST WAY to run the pipeline.

## Read annotation only

To optionally skip generating an assembly when running Hecatomb, 
the command is exactly the same as above with the addition of the `--skipAssembly` flag:

```bash
hecatomb run --reads fastq/ --profile slurm --skipAssembly
```

## Quicker read annotation

The pipeline bottleneck is the MMSeqs searches.
Use the `--fast` flag to run Hecatomb with less sensitive settings for MMSeqs.
In limited testing, we find it performs almost as well but with considerable runtime improvements.

```bash
hecatomb run --reads fastq/ --profile slurm --fast
```

## Specifying a host genome

Hecatomb includes a thorough host read removal step which utilises a processed host genome.
You can specify a host, or add your own.

By default, Hecatomb will use the human genome.
If your sample is from a different source you will need to specify the host genome for your sample source.

To see what host genomes are available:

```bash
hecatomb listHosts
```

The following should be available by default: 
bat, mouse, camel, celegans, macaque, rat, dog, cat, tick, mosquito, cow, human

So if you are working with mouse samples you would run:

```bash
hecatomb run --reads fastq/ --host mouse
```

## Add your own host genome

If the genome for the host you're working with isn't included in the available hosts, or you have a reference genome
which you think is better, you can add it with `addHost`.
This script will mask viral-like regions from your genome and add it to your Hecatomb host database.

You will need to specify the host genome FASTA file, as well as a name for this host.
Assuming you want to add the llama genome and the FASTA genome file is called `llama.fasta`:

```bash
hecatomb addHost --host llama --hostfa llama.fasta
```

You will then be able to run Hecatomb with your new host genome:

```bash
hecatomb run --reads fastq/ --host llama
```

This section assumes you have finished [Part 1: Preparation](tutorialPt1.md) by downloading the `bigtable.tsv.gz` and `metadata.tsv.gz`, 
and have loaded them in Rstudio into dataframes called `data` and `meta`.

## Inspect the dataframes

Let's look at the bigtable file:

```R
# in R/Rstudio
View(data)
```

[![](img/tuteDataTable.png)](img/tuteDataTable.png)

And the metadata:

```R
View(meta)
```

![](img/tuteMetaTable.png)

For the metadata table I've made sure the column header for the sample IDs are the same in both tables.
There was more metadata associated with the original samples, but we've simplified things here.

## Merge in the metadata

Merging in your metadata is easy.

The merge function with perform an inner join by default, or you can specify outer, and left- and right-outer.
This shouldn't matter if you have metadata for all of your samples.

```R
# save the merged tables to a new dataframe
dataMeta = merge(data, meta, by='sampleID')

# If your metadata is incomplete and you want to keep all samples
dataMeta = merge(data, meta, by='sampleID', all=T)
```

## Preliminary bigtable plots

Let's look at the raw alignments.
First, extract the viral hits to a new data frame. 

```R
viruses = dataMeta %>% 
    filter(kingdom=="Viruses")
```

I like to plot the alignment length against identity, and facet by viral family.
We show the different alignment types by color, and we can scale the point size by the cluster number.
I've added in `alpha=0.1` to set it to 10% opacity and the points will overlap a lot at this scale.

```R
ggplot(viruses) + 
    geom_point(
        aes(x=alnlen,y=pident,color=alnType,size=count),
        alpha=0.1) + 
    facet_wrap(~family)
```

[![](img/tuteAlnPidFam.png)](img/tuteAlnPidFam.png)

We can immediately see that a handful of viral families make up a majority of the viral hits.
You can use these plots to help guide filtering strategies.
We can divide the alignments into 'quadrants' by adding alignment length and percent identity thresholds,
for instance alignment length of 150 and percent identity of 75. 

```R
ggplot(viruses) +
    geom_point(
        aes(x=alnlen,y=pident,color=alnType,size=count),
        alpha=0.1) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
```

[![](img/tuteAlnPidFamQuad.png)](img/tuteAlnPidFamQuad.png)

We can see that for Adenoviridae and Parvoviridae the majority of hits occupy the top two quadrants, 
and we can be reasonably confident about these alignments.
For Podoviridae and Circoviridae, the majority of hits occupy the bottom two quadrants.
This could indicate that the viruses are only distantly related to the reference genomes in these families.

## Challenges

**Can you plot the Bacterial hits faceted by phylum?**

Don't facet the bacterial hits by family, there will be way too many panels.

**Can you plot the Adenoviridae hits for each sample?**

You can use `%>% filter(...)` inside the `ggplot()` function.


## Filtering strategies

Hecatomb is not intended to be a black box of predetermined filtering cutoffs that returns an immutable table of hits.
Instead, it delivers as much information as possible to empower the user to decide which hits they want to keep and which hits to purge.
Let's take our raw viral hits dataframe `viruses` and filter them to only keep the ones we are confident about.

The evalue is one of the most common metrics to use for filtering alignments.
Let's see what hits would be removed if we used a fairly stringent cutoff of 1e-20.

You can return a filtered dataframe like so:

```R
virusesFiltered = viruses %>% 
    filter(evalue<1e-20)
```

But let's instead add a flag to the original dataframe for purging hits:

```R
# mutate() will add or modify columns, ifelse() will return a value base on a condition
viruses = viruses %>% 
    mutate(filter=ifelse(evalue<1e-20,'pass','filter'))
```

And plot:

```R
ggplot(viruses) +
    geom_point(
        aes(x=alnlen,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family)
```

[![](img/tuteVirEvalFilt.png)](img/tuteVirEvalFilt.png)

The orange (or red?) sequences are destined to be removed, while the blue-green sequences will be kept.
Some viral families will be removed altogether, which is probably a good thing if they only have low quality hits.

Going back to the quadrant concept, you might only want to keep sequences above a certain length and percent identity:

```R
# this will overwrite the flags with the new designations
viruses = viruses %>% 
    mutate(filter=ifelse(alnlen>150 & pident>75,'pass','filter'))
```

Plot with the vertical and horizontal lines to match our cutoffs:

```R
ggplot(viruses) +
    geom_point(
        aes(x=alnlen,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
```

[![](img/tuteVirLenFilt.png)](img/tuteVirLenFilt.png)

There are many alignment metrics included in the bigtable for you to choose from.

**Other suggestions**

In our previous plots, we were colouring alignments by alignment type ('aa' or 'nt' under 'alnType')
Hecatomb will only attempt to identify nucleotide hits for a sequence if there is no protein hit.
We generally expect more protein than nucleotide hits as intergenic space is generally small in viral genomes.
A high nucleotide:protein hit ratio could be the result of false-positive hits, 
and you may wish to ignore nucleotide hits in some circumstances.

Hecatomb will try to assign taxonomy using a Lowest Common Ancestor approach ('LCA' under 'taxMethod').
If the LCA calculation fails, then Hecatomb defaults to the top hit ('TopHit' under 'taxMethod').
You may or may not wish to include the top hit-based annotations, 
or ignore top hit annotations if there are no accompanying LCA-based annotations.

## Challenge

**Filter your raw viral hits to only keep protein hits with an evalue < 1e-10**

Save it to `virusesFiltered` and move on to [Part 3: visualising annotations](tutorialPt3.md).
## About Snakemake profiles

Snakemake profiles are a must-have for running Snakemake pipelines on HPC clusters.
While they can be a pain to set up, you only need to do this once and then life is easy.
**If you just want to get up and running with as little effort as possible, [see the tutorial below](profiles.md#profile-installation-examples).**

For more information, check the [Snakemake documentation on profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), 
or our recent [blog post on Snakemake profiles](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated).

## Profiles for Hecatomb

Hecatomb ships with an example profile for the Slurm workload manager in `hecatomb/snakemake/profile/example_slurm/`.
The example _profile_ `config.yaml` file (not to be confused with the _Hecatomb_ config file) contains all the Snakemake 
options for jobs that are submitted to the scheduler.
Hecatomb expects the following in the cluster command:

 - `resources.time` for time in minutes
 - `resources.mem_mb` for requested memory in Mb
 - `threads` for requested CPUs
 
We have tried to use what we believe is the most common nomenclature for these variables in Snakemake pipelines 
in the hopes that Hecatomb is compatible with existing Snakemake profiles and the available 
[Cookiecutter profiles for Snakemake](https://github.com/Snakemake-Profiles).

We recommend redirecting STDERR and STDOUT messages to log files using the Snakemake variables `{rule}` and `{jobid}`, 
for instance like this `--output=logs/{rule}/{jobid}.out`.
You should also prepend the scheduler command with a command to make the log directories in case they don't exists
(it can cause errors for some schedulers), in this example like so: `mkdir -p logs/{rule}/ && sbatch ...`. 
This will make troubleshooting easier for jobs that fail due to scheduler issues.

The example profile includes a 'watcher' script.
Snakemake won't always pick up when a scheduler prematurely terminates a job, which is why we need a watcher script.
This line in the config file tells Snakemake how to check on the status of a job:
`cluster-status: ~/.config/snakemake/slurm/slurm-status.py` (be sure to check and update your file path).
The `slurm-status.py` script will query the scheduler with the jobid and report back to Snakemake on the job's status.

## Profile installation examples

We'll walk through two ways to set up a profile for the Slurm workload manager.
If your HPC uses a different workload manager, the process of installing a profile will be similar but different.

## Copy an example profile

We have provided an example profile for the Slurm workload manager that should work for most HPCs using Slurm.
Snakemake will look for profiles in your home directory at:

```text
~/.config/snakemake/
```

First create a directory for your new profile, we'll call it 'slurm':

```bash
mkdir -p ~/.config/snakemake/slurm
```

Now copy the files for the example slurm profile 
(you can view them [here on GitHub](https://github.com/shandley/hecatomb/blob/main/snakemake/profile/example_slurm/)):

```bash
# go to your new profile directory
cd ~/.config/snakemake/slurm/
# copy the files (either from GitHub or from where you installed Hecatomb)
wget https://raw.githubusercontent.com/shandley/hecatomb/ec1c62cfaddf29ade68cf4f33f4991fa07f9e6e0/snakemake/profile/example_slurm/config.yaml
wget https://raw.githubusercontent.com/shandley/hecatomb/ec1c62cfaddf29ade68cf4f33f4991fa07f9e6e0/snakemake/profile/example_slurm/slurm-status.py
```

This example includes the necessary `config.yaml` file for profiles and a watcher script called `slurm-status.py`.
Make the watcher script executable:

```bash
chmod +x ~/.config/snakemake/slurm/slurm-status.py
```

Done!
You can now use this profile with hecatomb:

```bash
hecatomb run --test --fast --profile slurm
```

## Create a profile with cookiecutter

Cookiecutter is a nifty tool for creating projects using a template.
For Snakemake profiles, Cookiecutter takes away a lot of the manual configuration steps involved with setting up a profile for a specific scheduler.
There are currently [Cookiecutter templates for Slurm, SGE, PBS, and several other workload managers](https://github.com/Snakemake-Profiles),
and Hecatomb is intended to be compatible with these profiles.

We will walk through installing the [Slurm profile using Cookiecutter](https://github.com/Snakemake-Profiles/slurm).
To begin, create a new directory for your profile:

```bash
mkdir -p ~/.config/snakemake/slurm
```

Move to this directory and run the Cookiecutter command for this profile:

```bash
cd ~/.config/snakemake/slurm/
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

Follow the prompts and you're done!
For our system, we did not need to specify anything, but you may need to specify account information for billing etc.
There is more detail on the [GitHub page for this profile](https://github.com/Snakemake-Profiles/slurm).

Use the profile with Hecatomb:

```bash
hecatomb run --test --fast --profile slurm
```
![](img/hecatombLogo.png)

A hecatomb is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to 'sacrifice' false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

![](img/hecatombBanner.png)

Hecatomb was developed in response to the challenges associated with the detection of viral sequences in metagenomes. 
Virus detection or virome profiling is typically performed on samples containing extensive host nucleic acid (e.g. a 
tissue biopsy) or nucleic acid from a variety of other organisms such as bacteria, fungi and archaea from soil samples 
or bacteria, fungi, archaea and food from mammalian stool samples. 
All of these non-viral nucleic acid types constitute the background and present a variety of issues that must be 
evaluated in order to confidently detect viral sequences.

Detection and quantitative analysis of viral sequences from mixed metagenomic data has other unique challenges. 
Viruses have a propensity to rapidly acquire mutations and evolve to be highly diverged from reference database sequences.
Coding and non-coding regions of viral genomes have disparate challenges (e.g. long-terminal 5' and 3' repeats).
Viruses are further clssified as either single or double stranded, and either DNA or RNA.
False-positives are rampant in many studies due to significant sequence homology between viruses and other domains of life.

Hecatomb utilises a unique search strategy to address these issues and enable researchers to accurately identify and quantify 
viruses in their metagenome samples.
## System requirements

We have set some sensible defaults for Hecatomb on the following minimum recommended system requirements:

 - 32 CPUs (or 16 CPUs with hyperthreading)
 - 64 GB of RAM
 - approximately 55 GB HDD space (for the databases)
 - Additional HDD space for temporary and output files (highly dependent on input file sizes)

Larger datasets may require more RAM (and CPUs to speed things along).
As such it is highly recommended to run Hecatomb on a HPC cluster.

## Dependencies

Hecatomb was developed as a set of [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflows, which are 
controlled by a python launcher for your convenience.
While the pipeline utilises a number of other programs and tools to run, these are all managed by the pipeline using 
[conda](https://docs.conda.io/en/latest/) and [mamba](https://github.com/mamba-org/mamba).

Install Hecatomb via conda:

```bash
# This command will create a new conda env called 'hecatomb' and will install Hecatomb and all of it's dependencies.
conda create -n hecatomb -c conda-forge -c bioconda hecatomb
```

That's it!

```bash
# To use Hecatomb, activate your new conda env.
conda activate hecatomb

# Check that it's installed
hecatomb -h
```

## Customisation

If you're running Hecatomb on a cluster it is highly recommended to use Snakemake profiles 
([Set up a profile for Snakemake](profiles.md)).

You may also want to customise the available resources for your system, whether you're using a cluster or running 
locally ([Advanced configuration](configuration.md)).

## Download the databases

Before running Hecatomb for the first time you will need to download the databases.
You will only need to do this step once (unless we update the databases).
You can rerun this step as much as you like; the pipeline will only download any database files that are missing.

```bash
# Either, run locally
hecatomb install

# Or, for a HPC cluster using a Snakemake profile (replace 'slurm' with your profile name)
hecatomb install --profile slurm
```

## Run the test dataset

Hecatomb comes with a test dataset that you can run which will take a few hours to complete.
Use `--test` in place of specifying your read directory with `--reads`.

```bash
# run locally
hecatomb run --test

# run on cluster using a Snakemake profile
hecatomb run --test --profile slurm
```

## Build the docs

These document pages can be built from the repo like so.

```bash
# install mkdocs
pip install mkdocs

# cd to your install directory
cd /path/to/hecatomb

# build
mkdocs build
```
# Contig read-based annotations

Hecatomb will produce an assembly and annotate your reads.
It will also combine the two, mapping your annotated reads to the assembly and 
building a table to combine the mapped coordinates and the reads' annotations.

Load the table and have a look:

```R
contigTable = read.csv('contigSeqTable.tsv.gz',header=T,sep='\t')

View(contigTable)
```

[![](img/tuteCtgTbl.png)](img/tuteCtgTbl.png)

It contains the sequences from the seqtable with contig mapping coordinates, 
mapping quality, counts, and taxonomic annotations.
You could also join this with the bigtable to bring in everything for each sequence.
The premise behind this table was to enable dynamic taxonomic assignments for 
the assembly contigs.
This is still a work in progress,
but we will show some example plots to demonstrate the idea.

Let's look at what reads are mapping to contig_6.
We will separate out the reads by genus and color by family:

```R
ggplot(contigTable %>% 
        filter(contigID=='contig_6')
    ) +
    geom_point(aes(x=start,y=genus,color=family))
```

[![](img/tuteCtg6.png)](img/tuteCtg6.png)

There is a nice clear consensus that this is a Mastadenovirus.
The handful of other reads we could fairly confidently reassign as such if we wanted.
Given the clear consensus, it might be possible to assign it to a specific species:

```R
ggplot(contigTable %>% 
        filter(contigID=='contig_6')
    ) +
    geom_point(aes(x=start,y=species,color=family))
```

[![](img/tuteCtg6Sp.png)](img/tuteCtg6Sp.png)

Unfortunately its closest hit in the reference databases is an unclassified species.
However, its next closest hits are to Rhesus adenoviruses, which makes sense for this Macaque sample.

Lets look at contig_56:

```R
ggplot(contigTable %>% 
           filter(contigID=='contig_56')
) +
    geom_point(aes(x=start,y=family,color=kingdom))
```

[![](img/tuteCtg56.png)](img/tuteCtg56.png)

Most of the hits remain unclassified, and the classified hits are split between kingdoms.
Let's see what species the sequences are hitting:

```R
ggplot(contigTable %>% 
           filter(contigID=='contig_56')
) +
    geom_point(aes(x=start,y=species,color=kingdom))
```

[![](img/tuteCtg56Sp.png)](img/tuteCtg6Sp.png)

The annotates hits are split between _Bacteroides_ phage species, 
and a _Ktedonobacter_ bacteria species. 
The contig could represent a novel phage species, or a novel prophage in a bacteria genome.

We think this is a cool analysis with lots of potential, but we would love your feedback!
This tutorial will walk through the process of running Hecatomb and performing some preliminary plots and analyses in R/Rstudio.
This first section details the Hecatomb run.

**You won't need to actually run Hecatomb, you can simply [download the files](tutorialPt1.md#download-the-files) and continue with the rest of the tutorial.**

## System information

For this tutorial I'll be running Hecatomb on a 16-core/32-thread workstation with 64 GB of RAM running Ubuntu 18.04 LTS.
While this is fine for smaller datasets, larger datasets will require HPC resources.

## New install

```bash
# create new conda env and install hecatomb
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

# activate env
conda activate hecatomb

# install the database files
hecatomb install
```

## Run Hecatomb

We will run hecatomb on the test dataset, which is a subset of the dataset described in [Handley et al (2016)](https://doi.org/10.1016/j.chom.2016.02.010).
We use the fast MMSeqs settings with the default 32 threads and generate a summary report. 
This will give us an assembly and some read annotations.

```bash
Hecatomb run --test --fast --report
```

We should now have all the files we need!

## Download the files

Download the following (gzipped) files generated by Hecatomb:

- [bigtable.tsv.gz](https://cloudstor.aarnet.edu.au/plus/s/549SmspixHnryiK/download)
- [seqtable.fasta.gz](https://cloudstor.aarnet.edu.au/plus/s/pZ0GXoTYgPf9aF4/download)
- [taxonLevelCounts.tsv.gz](https://cloudstor.aarnet.edu.au/plus/s/wGfkmgmsZhGkUaf/download)
- [contigSeqTable.tsv.gz](https://cloudstor.aarnet.edu.au/plus/s/flovOcyWIc94RPx/download)
- [sankey.svg](https://cloudstor.aarnet.edu.au/plus/s/gwfa9xcA9m0aRmg/download)
- [krona.html](https://cloudstor.aarnet.edu.au/plus/s/hFo1Rnx8h3rTXSu/download)
- [assembly.fasta](https://cloudstor.aarnet.edu.au/plus/s/bmTo2jzwB65eRsr/download)

I also prepared a metadata file. 
The samples in the test dataset are the samples for the deceased Macaques from the original study.
The metadata contains the individuals' gender and the vaccine that was administered.

- [metadata.tsv.gz](https://cloudstor.aarnet.edu.au/plus/s/65xBlEe4TNxvOCp/download)

The `Sankey.svg` and `krona.html` are automatically-generated plots.
You can have a look at these files but don't need to do anything with them for the tutorial.
The `bigtable.tsv.gz` is the read-based annotations and will form the basis for most of the subsequent analysis.
The `contigSeqTable.tsv.gz` combines the read-based annotations with the contig mapping coordinates.
The `taxonLevelCounts.tsv.gz` file summarises the sequence counts for each sample, at every taxon level, 
and will make it easy to perform some preliminary comparisons.

## Optional: prefilter the bigtable.tsv

The dataset in this tutorial is small and manageable, but if you're working on a large dataset, 
the size of the bigtable can be too big for loading into R. 
If you run into this situation check out [prefiltering the bigtable](advanced.md#prefilter-the-bigtable).

## Setting up R/Rstudio

The remainder of the tutorial will be in [R](https://www.r-project.org/) and [Rstudio](https://www.rstudio.com/), 
and will use the packages [dplyr](https://dplyr.tidyverse.org/), [tidyr](https://tidyr.tidyverse.org/), and [ggplot2](https://ggplot2.tidyverse.org/).
The installation and packages for this tutorial have been tested with a fresh installation of R (4.1.2) and RStudio (2021.09.1+372) for Windows.

Assuming you have installed R and Rstudio for your operating system, you should be able to install the packages like so:

```R
# in R/Rstudio
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")
```

Alternatively, you can install the whole tidyverse with:

```R
install.packages("tidyverse")
```

We'll also need 'rstatix' for some of the statistical tests:

```R
install.packages("rstatix")
```

Once installed, load the packages:

```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
```

## Load the data

Make sure you set your working directory and have downloaded the files into this directory.
We will load the `bigtable.tsv.gz` into a dataframe called `data` with `load.csv()`.
The file is tab-separated and contains a header.

```R
data = read.csv('bigtable.tsv.gz',header=T,sep='\t')
```

Next, we will load the `metadata.tsv.gz` file into a dataframe called `meta` in the same way.

```R
meta = read.csv('metadata.tsv.gz',header=T,sep='\t')
```

We'll load `taxonLevelCounts.tsv.gz` while we're at it.

```R
taxonCounts = read.csv('taxonLevelCounts.tsv.gz',header=T,sep='\t')
```

Keep the `contigSeqTable.tsv.gz` file handy, we'll load that into R as well later on.

You can now move onto [Part 2: filtering annotations](tutorialPt2.md).
## Prefilter the bigtable

If your bigtable is too big, you can remove the stuff you don't want using the linux command line before loading it into R or Python.

**Remove low quality hits**

Use `awk` to apply evalue/bitscore/alignemnt lenght/etcetera cutoffs:

```bash
# copy the header for your new file
head -1 bigtable.tsv > newBigtable.tsv

# filter low quality hits, e.g. using an evalue cutoff of say 1e-20
awk -F '\t' '$7<1e-20' bigtable.tsv >> newBigtable.tsv
```

**Remove irrelevant hits**

If you only plan on analysing viruses or say protein hits, you can remove everything else with `awk`:

```bash
# copy the header
head -1 bigtable.tsv > newBigtable.tsv

# return only rows where "Viruses" is in the "kingdom" column
akw -F '\t' '$24=="Viruses"' bigtable.tsv >> newBigtable.tsv

# alternatively, only return the "aa" protein-based hits
awk -F '\t' '$5=="aa"' bigtable.tsv >> newBigtable.tsv
```

**Remove irrelevant columns**

If you've finished prefiltering, or you only plan on using say the evalues for filtering, 
you can remove all the other alignment metrics, this time using `cut`:

```bash
# remove all the alignment metrics except for evalue
cut --complement -f8-22 bigtable.tsv > newBigtable.tsv
```

## Retrieving sequences of interest

Upon analysing your bigtable, you may have a collection of interesting hits that you want to investigate further.
To retrieve the relevant sequences from your `seqtable.fasta` file you will only need a list of the sequence IDs--the first column in the bigtable.

For instance, get all seqIDs for the Viral family 'Flaviviridae' and save it to a file:

```R
# in R using tidyr/dplyr get the seq IDs
flaviviridaeSeqs = data %>% 
    filter(family=='Flaviviridae') %>% 
    pull(seqID)

# print to a file
lapply(flaviviridaeSeqs, write, "flavSeqIDs.list", append=TRUE, ncolumns=1)
```

Use Samtools to grab all these sequences from the seqtable.fasta file:

```shell
# in bash
samtools faidx seqtable.fasta -r flavSeqIDs.list > flaviviridaeSeqs.fasta
```
The pipeline outputs a number of files for further analysis and exploration, as well as to provide an overview of the 
read preprocessing and distribution.

## Report

`report.html`

This file is generated by Snakemake and outlines a lot of information relating to the Hecatomb run.
Under the Results tabs are summary files for things like reads for each sample following different preprocessing steps
as well as some summary plots.

## SeqTable

`hecatomb_out/RESULTS/seqtable.fasta`

The SeqTable is the primary output of the read preprocessing and serves as the input for Taxonomic assignment.
It is composed of all the representative sequences from the clustered reads for all samples.
Samples are clustered individually, and the seq IDs for this fasta file follows the format `>sampleID:count:seqNumber`.
Here, `count` is the number of reads in that cluster which is important for statistical exploration.
Sequences are numbered sequentially (`seqNumber`) to ensure unique IDs.

## BigTable

`hecatomb_out/RESULTS/bigtable.tsv`

The BigTable is the main output of Taxonomic assignment and can be directly imported into R or Python.
The BigTable combines the seqtable IDs with their sampleID, counts, normalised counts, alignment information, taxonomic assignments and Baltimore classification.
This file is big, hence the name, but is designed to make merging with sample metadata, plotting, and statistical interrogation as easy as possible.

The header looks like this:

```text
seqID  sampleID  count  CPM  alnType  targetID  evalue  pident  fident  nident  mismatches  qcov  tcov  qstart  qend  qlen  tstart  tend  tlen  alnlen  bits  targetName  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

## TaxonLevelCounts

`hecatomb_report/taxonLevelCounts.tsv`

This file is derived from the BigTable and summarises the total sequence counts, for each sample, at all taxonomic levels.
The TaxonLevelCounts combines the sampleID with the taxonomic level for which the counts refer, the full taxonomic path, 
the taxon name, and the total and normalised read counts.
The purpose of this file is to expedite statistical interrogation of your data.
For instance, if you wanted to compare the numbers of say Flaviviridae reads between two groups of samples, 
those counts have already been collected, and you can simply run your analysis and plotting on the relevant slice of the table.  

The file looks something like this:

```text
sampleID    taxonLevel  taxonPath                                   taxonName       count   CPM
sample1     Kingdom     k_Bacteria                                  Bacteria        3162    3178.818
sample1     phylum      K_Viruses,p_Phixviricota                    Phixviricota    1216    1222.467
sample1     class       K_Viruses,p_Uroviricota,c_Caudoviricetes    Caudoviricetes  1234    1240.564
etc.
```

## Assembly

`hecatomb_out/RESULTS/assembly.fasta`

These are the contigs generated UNLESS you run Hecatomb with the `--skipAssembly` flag.
The assembly is used for producing the ContigSeqTable and ContigKrona plots, as well as the direct contig annotations.

## CONTIG ANNOTATIONS

TODO

## ContigSeqTable

`hecatomb_out/RESULTS/contigSeqTable.tsv`

The ContigSeqTable combines the read mapping information for the assembly with the read-based taxonomic assignments.
This file is intended to assist the user in identifying and binning assembly contigs by applying a consensus approach to contig taxonomic assignment.
The file includes the positional mapping information and can also enable investigation of more complex features such as 
chimeric contigs, recombination or horizontal transfer events.

The header looks like this:

```text
contigID  seqID  start  stop  len  qual  count  CPM  alnType  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

## Sankey

`hecatomb_report/Sankey.svg`

The sankey diagram shows the fate of all reads throughout the preprocessing and read-based taxonomic assignment steps.
It serves to visualise the read filtering and distribution of taxonomic assignment methods and give you an overall impression of how well things went.
The sankey diagram produced for the test dataset is shown below. 
This dataset is relatively rich in viral sequences and yet the majority of reads are either filtered or non-viral (that we know of).

[![](img/Sankey.svg)](img/Sankey.svg)

## krona.html and contigKrona.html

`hecatomb_report/krona.html`

`hecatomb_report/contigKrona.html`

The Krona plots are to assist in visual exploration of the read annotations.
krona.html is derived from the bigtable and shows the raw distribution of taxon assignments.
contigKrona.html is derived from the contigSeqTable and includeds the taxon assignment method (either tophit or LCA).
The contigKrona plot helps to visualise the distributions of topHit versus LCA assigned reads as well as the 
distributions over contigs of the identified species, and the distribution of taxonomic assignments for each contig.
## Where can I get help or support?

[Open an issue on github](https://github.com/shandley/hecatomb/issues).
There is no question too stupid.
The pipeline is still in the early days of development, so we are expecting that there will still be bugs in the code.

## The pipeline died, what went wrong?

A great many things may have gone wrong!
You can simply try rerunning the pipeline, and it will pick up where it left off. 

We try to make the logging as thorough as possible and have meaningful error messages where we can.
Look through the terminal output produced by Snakemake for the failed job.
It will appear in red text for the terminal, or you can look through the snakemake log in `.snakemake/logs/`.
Each rule has its own directory in `hecatomb_out/STDERR/` for log files, 
but the Snakemake error code will tell you exactly what file to look at.

If you're running with a profile on a HPC, the job may have failed for reasons relating to the scheduler.
For instance, the memory or time may have exceeded what was requested.
Look in your scheduler logs to see if this is the case.
In the example profile, the slurm logs are saved in `logs/{rule}/{jobid}`.
The rule name and jobid are supplied in the Snakemake error messages, making it easy to find the relevant logs.

## I have too much data and am having memory problems

You can run your samples in batches, and even split samples into separate runs if you really need.
The `bigtable.tsv` files from different runs can be concatenated (just remember to remove the extra headers),
and your assemblies can be combined using Flye with the `--subassemblies` function.

If your sequencing depth is exceptionally high you could also subsample your FASTQ files 
(and run a couple at full coverage as a control).

Depending on which jobs are failing, you can reduce the cluster threshold (`CLUSTERID` in your `config.yaml`) to help 
reduce the size of your seqtable, and use more stringent e-value cutoffs for your MMSeqs searches 
(also in your `config.yaml`) to reduce the output file sizes.

## Hecatomb takes too long, what can I do?

First, try running Hecatomb with the `--fast` flag.
The MMSeqs steps are by far the most time consuming steps. 
The `--fast` flag will tell Hecatomb to use MMSeqs settings that are much much faster, but not quite as sensitive.
You should also configure your installation to utilise as many CPUs and as much memory as possible.
See [default resources config](https://hecatomb.readthedocs.io/en/latest/advanced/#default-resources) for more info.
If you don't need an assembly, you can skip those steps with `--skipAssembly` which will also save a bit of time.

## I've run the pipeline, now what?

Have a look at [the tutoria](#) which goes through some example plots and analyses.

This section assumes you have finished Tutorials [Part 1](tutorialPt1.md) and [Part 2](tutorialPt2.md)

## Unfiltered taxon counts

The `taxonLevelCounts.tsv` file is intended to make it quick and easy to compare virus hit counts across samples.
Let's take a look at it:

```R
View(taxonCounts)
```

[![](img/taxCountTable.png)](img/taxCountTable.png)

The file contains the counts and normalised counts for each taxonomic level, for each sample, 
based on the raw (unfiltered) hit counts.
If we wanted to compare the number of Adenoviridae hits between each sample we could pull out the Adenoviridae counts and plot them:

```R
# get adenoviridae counts
adenoCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Adenoviridae')

# plot
ggplot(adenoCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
```

[![](img/tuteAdenoBar.png)](img/tuteAdenoBar.png)

We can take this one step further and plot all the viral family normalised counts for each sample, in a stacked bar chart:

```R
# get all viral family counts
viralCounts = taxonCounts %>% 
    filter(taxonLevel=='family',grepl('k_Viruses',taxonPath))

# plot
ggplot(viralCounts) +
    geom_bar(aes(x=sampleID,y=count,fill=taxonName),position='stack',stat='identity') +
    coord_flip()
```

[![](img/tuteViralCounts.png)](img/tuteViralCounts.png)

## Generating taxon counts

The `taxonLevelCounts.tsv` file is convenient for comparing the raw counts,
but you will likely want to generate new counts from your _filtered_ hits.
Recreate the above plot from the filtered hits by first summing the counts
or normalisedCounts, e.g. at the family level:

```R
# Answer for "Challenge: Filter your raw viral hits to only keep protein hits with an evalue < 1e-10"
virusesFiltered = viruses %>% 
    filter(alnType=='aa',evalue<1e-10)

# collect the filtered counts
viralFiltCounts = virusesFiltered %>% 
    group_by(sampleID,family) %>% 
    summarise(n = sum(CPM))
```

Then plot again. 
This time we use `position='fill'` to make it look like 16s data, so we can confuse people.
We'll also add `theme_bw()` because grey is ugly:

```R
ggplot(viralFiltCounts) +
    geom_bar(aes(x=sampleID,y=n,fill=family),position='fill',stat='identity') +
    coord_flip() +
    theme_bw()
```

[![](img/tuteViralFiltCounts.png)](img/tuteViralFiltCounts.png)

These count tables we will use for plotting and some statistical comparisons.

# Challenge

**Make a stacked bar chart of the viral families for the Male and Female monkeys**

![](img/tuteGenderCounts.png)

# Visualising groups

We have a few viral families that are very prominent in our samples.
Let's see if there is a difference in viral loads according to our sample groups.
Collect sample counts for _Microviridae_.
Include the metadata group in `group_by()` so you can use it in the plot.

```R
microCounts = virusesFiltered %>% 
    group_by(family,sampleID,MacGuffinGroup) %>% 
    filter(family=='Microviridae') %>% 
    summarise(n = sum(CPM))
```

And plot. I like jitter plots but boxplots or violin plots might work better if you have hundreds of samples.

```R
ggplot(microCounts) +
    geom_jitter(aes(x=MacGuffinGroup,y=n),width = 0.1) +
    theme_bw()
```

![](img/tuteMicrovirJitter.png)

Let's do the same for _Podoviridae_.

```R
# collect counts
podoCounts = virusesFiltered %>% 
    group_by(family,sampleID,MacGuffinGroup) %>% 
    filter(family=='Podoviridae') %>% 
    summarise(n = sum(CPM))

# plot
ggplot(podoCounts) +
    geom_jitter(aes(x=MacGuffinGroup,y=n),width = 0.1) +
    theme_bw()
```

![](img/tutePodoJitter.png)

# Challenge

**Could gender be a good predictor of viral load for these families?**

While the MacGuffinGroup looks promising for _Podoviridae_, 
we'll need to [move on to Part 4: statistical tests](tutorialPt4.md) to find out for sure. 
A complete list of all database files is available in the Hecatomb `config.yaml` file.

## Contaminants

`databases/contaminants/`

This database consists of a collection of NEBNext and TruSeq adapters, primerB sequences, 
and a collection of cloning vectors that are common sources of false positives. 
These sequences are used in preprocessing to trim reads of non-biological contaminant sequences.

## Host genomes

`databases/host/`

Contamination of viral metagenomes by host DNA can be a significant burden and source of false positive in viral annotation.
We use a host reference genome to filter out any host reads prior to annotation.
However, host genomes typically contain large amounts of virus-like and virus-derived sequence.
This could lead to erroneous removal of true viral reads.
We therefore must process host reference genomes to mask viral-like sequence prior to using it for preprocessing.
Hecatomb comes with several processed host genomes and a tool for users to add their own genomes.
See [Adding your own host genome](usage.md#adding-your-own-host-genome) for more info.

## Amino acid databases

`databases/aa/`

**Primary AA**

The primary AA database is used to quickly identify any reads that are viral-like.
This database is all UniProt viral protein entires clusterd at 99% identity. 
This clustering greatly reduces the size of the database with minimal sacrifice to information content. 
For example, there were 4,635,541 protein entires for viruses in UniProt on October 29, 2020. 
When clustered at 99% identity only 1,840,337 representative sequences remained, 
this is about a 60% reduction in target query space.

**Secondary AA**

The secondary AA database is used to check if any viral-like reads (from the primary AA search) have a better non-viral match.
TODO ...

## Nucleotide databases

`databases/nt/`

**Primary NT**

The primary NT database is used to quickly identify any viral-like reads that were not identified in the primary AA search.
TODO ...

**Secondary NT**

The secondary NT database is used to check if the viral-like reads (from the primary NT search) have a better non-viral match.
TODO ...

## Taxonomy

`databases/tax/`

Hecatomb utilises the NCBI taxonomy database (taxdump) with [TaxonKit](https://github.com/shenwei356/taxonkit) 
for converting taxon IDs into complete lineages for the output bigtable.

## Tables

`databases/tables/`

These are a collection of supplementary tables.
Primarily, the Baltimore classifications are provided here to be merged with the bigtable annotations.
This section assumes you have completed [Tutorial Part3](tutorialPt3.md).

# Compare viral loads

In part 3, we compared the viral counts between the two sample groups for _Podoviridae_,
and it appeared as though group B had more viral sequence hits on average than group A.
We can compare the normalised counts for these two groups to see if they're significantly different.

**Student's T-test**

Let's check out the dataframe we made earlier that we'll be using for the test:

```R
View(podoCounts)
```

![](img/tutePodoCnts.png)

We'll use the base-r function `t.test()`, which takes two vectors--one with the 
group A counts and one with the group B counts.
We can use the `filter()` and `pull()` functions within the `t.test()` function like so:

```R
t.test(
    podoCounts %>% 
        filter(MacGuffinGroup=='A') %>% 
        pull(n),
    podoCounts %>% 
        filter(MacGuffinGroup=='B') %>% 
        pull(n),
    alternative='two.sided',
    paired=F,
    var.equal=T)    
```

```text
	Two Sample t-test

data:  podoCounts %>% filter(MacGuffinGroup == "A") %>% pull(n) and podoCounts %>% filter(MacGuffinGroup == "B") %>% pull(n)
t = -5.3033, df = 8, p-value = 0.0007255
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -615.5267 -242.4563
sample estimates:
mean of x mean of y 
 81.11595 510.10744 
```

**Note: for all the following plots, I've manually added the significance asterisks.**

![](img/tutePodoJitterTTest.png)

**Wilcoxon test**

You might prefer to perform a Wilcoxon test; the syntax is very similar to the t.test:

```R
wilcox.test(
    podoCounts %>% 
        filter(MacGuffinGroup=='A') %>% 
        pull(n),
    podoCounts %>% 
        filter(MacGuffinGroup=='B') %>% 
        pull(n),
    alternative='t',
    paired=F)
```

```text
	Wilcoxon rank sum exact test

data:  podoCounts %>% filter(MacGuffinGroup == "A") %>% pull(n) and podoCounts %>% filter(MacGuffinGroup == "B") %>% pull(n)
W = 0, p-value = 0.007937
alternative hypothesis: true location shift is not equal to 0
```

![](img/tutePodoJitterWilc.png)

**Dunn's test**

Let's use Dunn's test to check all the major families at the same time.
Dunn's is good for if you have three or more categories for a metadata field, such as our vaccine column.
First find out what the major families are by summing the hits for each family and sorting the table.

```R
# collect the family counts
viralFamCounts = virusesFiltered %>% 
    group_by(family) %>% 
    summarise(n=sum(CPM)) %>% 
    arrange(desc(n))

# update factor levels to the sorted order
viralFamCounts$family = factor(viralFamCounts$family,levels=viralFamCounts$family)
```

and plot:

```R
ggplot(viralFamCounts) +
    geom_bar(aes(x=family,y=n),stat='identity') +
    coord_flip()
```

![](img/tuteFamCnts.png)

Let's focus on _Siphoviridae_, _Adenoviridae_, _Podoviridae_, and _Microviridae_.
Collect summary counts for these families for each sample and include the metadata we want to use:

```R
viralMajorFamCounts = viruses %>% 
    filter(family %in% c('Siphoviridae','Adenoviridae','Podoviridae','Microviridae')) %>% 
    group_by(sampleID,family,vaccine) %>% 
    summarise(n=sum(CPM))
```

Now let's do the dunn's test for these families:

```R
viralMajorFamCounts %>% 
    group_by(family) %>% 
    dunn_test(n ~ vaccine,p.adjust.method='holm') %>% 
    add_significance()
```

```text
# A tibble: 12 × 10
   family       .y.   group1     group2        n1    n2 statistic      p  p.adj p.adj.signif
 * <chr>        <chr> <chr>      <chr>      <int> <int>     <dbl>  <dbl>  <dbl> <chr>       
 1 Adenoviridae n     Ad_alone   Ad_protein     3     2     0.467 0.641  1      ns          
 2 Adenoviridae n     Ad_alone   sham           3     4     0.438 0.661  1      ns          
 3 Adenoviridae n     Ad_protein sham           2     4    -0.105 0.916  1      ns          
 4 Microviridae n     Ad_alone   Ad_protein     4     2    -1.43  0.153  0.458  ns          
 5 Microviridae n     Ad_alone   sham           4     4    -0.584 0.559  0.681  ns          
 6 Microviridae n     Ad_protein sham           2     4     0.953 0.340  0.681  ns          
 7 Podoviridae  n     Ad_alone   Ad_protein     4     2     0.477 0.634  0.681  ns          
 8 Podoviridae  n     Ad_alone   sham           4     4     1.75  0.0798 0.240  ns          
 9 Podoviridae  n     Ad_protein sham           2     4     0.953 0.340  0.681  ns          
10 Siphoviridae n     Ad_alone   Ad_protein     4     2     2.48  0.0132 0.0395 *           
11 Siphoviridae n     Ad_alone   sham           4     4     1.40  0.161  0.322  ns          
12 Siphoviridae n     Ad_protein sham           2     4    -1.33  0.182  0.322  ns  
```

There's only one comparison that is significant.
Let's put it on a plot

```R
ggplot(viralMajorFamCounts) +
    geom_jitter(aes(x=vaccine,y=n)) +
    facet_wrap(~family)
```

![](img/tuteDunns.png)

# Compare presence/absence

You might not care about viral loads and instead are just interested in comparing the presence or absence of viruses.
For this you could use a Fisher's exact test.
To perform this test you need to assign a presence '1' or absence '0' for each viral family/genus/etc for each sample.
What number of hits you use for deciding if a virus is present is up to you.

I want to be sure about the alignments, so I'll apply some stringent filtering cutoffs.
Then I'll assign anything with _any_ hits as 'present' for that viral family.
Let's look at _Myoviridae_ ... for no particular reason.
 
```R
# apply a stringent filter
virusesStringent = viruses %>% 
    filter(evalue<1e-30,alnlen>150,pident>75,alnType=='aa')

# count hits for Myoviridae and score presence
# NOTE: the absent samples WONT be in the table yet, we have to add them in after
myovirPresAbs = virusesStringent %>% 
    filter(family=='Myoviridae') %>%
    group_by(sampleID) %>% 
    summarise(n=sum(CPM)) %>%
    mutate(present=ifelse(n>0,1,0))

# merge in the metadata
# Using all=T will do an outer join and bring in the absent samples
myovirPresAbs = merge(myovirPresAbs,meta,by='sampleID',all=T)

# convert na's to zeros
myovirPresAbs[is.na(myovirPresAbs)] = 0
```

To do the Fisher's exact test we need to specify a 2x2 grid;
The first column will be the number with _Myoviridae_ for each group.
The second column will be the numbers without for each group.

```R
# matrix rows
mtxGroupA = c(
    myovirPresAbs %>% 
        filter(MacGuffinGroup=='A',present==1) %>% 
        summarise(n=n()) %>% 
        pull(n),
    myovirPresAbs %>% 
        filter(MacGuffinGroup=='A',present==0) %>% 
        summarise(n=n()) %>% 
        pull(n))
mtxGroupB = c(
    myovirPresAbs %>% 
        filter(MacGuffinGroup=='B',present==1) %>% 
        summarise(n=n()) %>% 
        pull(n),
    myovirPresAbs %>% 
        filter(MacGuffinGroup=='B',present==0) %>% 
        summarise(n=n()) %>% 
        pull(n))

# create the 2x2 matrix
myovirFishMtx = matrix(c(mtxGroupA,mtxGroupB),nrow = 2)

# this bit is not necessary, but lets add row and col names to illustrate the matrix layout
colnames(myovirFishMtx) = c('GroupA','GroupB')
row.names(myovirFishMtx) = c('present','absent')

# view
View(myovirFishMtx)
```

![](img/tuteFishMtx.png)

```R
# Run Fisher's exact test
fisher.test(myovirFishMtx)
```

```text
	Fisher's Exact Test for Count Data

data:  myovirFishMtx
p-value = 0.04762
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.000000 0.975779
sample estimates:
odds ratio 
         0 
```

We don't have many samples, so our significance won't be great regardless.

In [Part 5](tutorialPt5.md) we will look at the contigs' read-based annotations.## Recommended customisation

If you're running Hecatomb on a HPC cluster, we absolutely recommend setting up a 
[Snakemake profile](profiles.md).

We also recommend reviewing the Snakemake `config.yaml` file in your Hecatomb installation directory.
The config file will be at `hecatomb/snakemake/config/config.yaml` and you can find your installation directory with:

```bash
which hecatomb
```

## Changing the Hecatomb configuration

The Hecatomb configuration file `hecatomb/snakemake/config/config.yaml` contains settings related to resources and 
cutoffs for various stages of the pipeline. The different config settings are outlined further on.
You can permanently change the behaviour of your Hecatomb installation by modifying the values in this config file.
Alternatively, you can specify your own config file. To do this, you first copy the system default config file like so:

```bash
hecatomb config
```

You can then edit your new `hecatomb.config.yaml` file to suit your needs and used it in a Hecatomb run like so:

```bash
hecatomb run --configfile hecatomb.config.yaml
```

## Database location

The databases are large (~55 GB) and if your Hecatomb installation is on a partition with limited on space,
you might want to specify a new location to house the database files. 
By default, this config setting is blank and the pipeline will use the install location (`hecatomb/databases/`).
You can specify the directory in the Hecatomb config file (`hecatomb/snakemake/config/config.yaml`) under `Databases: `, 
e.g:

```yaml
Databases: /scratch/HecatombDatabases
```

and run or rerun the installation 

```bash
hecatomb install
```

## Default resources

The Hecatomb config file contains some sensible defaults for resources.
While these should work for most datasets, they may fail for larger ones.
You may also have more CPUs etc at your disposal and want to minimise runtime of the pipeline.
Currently, the slowest steps are the MMSeqs searches; increasing the CPUs and RAM could significantly improve runtime.
The other settings (for Megahit and Minimap2, BBTools, and misc) will probably only show modest improvement.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Memory for MMSeqs in megabytes (e.g 64GB = 64000, recommend >= 64000)
MMSeqsMem: 64000
# Threads for MMSeqs (recommend >= 16)
MMSeqsCPU: 32
# Max runtime in minutes for MMSeqs (this is only enforced by the Snakemake profile)
MMSeqsTimeMin: 5760  # 4 days

# Memory for BBTools in megabytes (recommend >= 16000)
BBToolsMem: 16000
# CPUs for BBTools (recommend >= 8)
BBToolsCPU: 8

# Memory for Megahit/Flye in megabytes (recommend >= 32000)
MhitMem: 32000
# CPUs for Megahit/Flye in megabytes (recommend >= 16)
MhitCPU: 16

# Memory for slightly RAM-hungry jobs in megabytes (recommend >= 16000)
MiscMem: 16000
# CPUs for slightly RAM-hungry jobs (recommend >= 2)
MiscCPU: 2

# Default memory in megabytes (for use with --profile)
defaultMem: 2000
# Default time in minutes (for use with --profile)
defaultTime: 1440
# Default concurrent jobs (for use with --profile)
defaultJobs: 100
```

## Preprocessing settings

There are many filtering etc. cutoff values that are specified in the Hecatomb config file.
For instance `READ_MINLENGTH: ` specifies the minimum allowed read length after trimming.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Preprocessing
QSCORE: 15 # Read quality trimming score (rule remove_low_quality in 01_preprocessing.smk)
READ_MINLENGTH: 90 # Minimum read length during QC steps (rule remove_low_quality in 01_preprocessing.smk)
CONTIG_MINLENGTH: 1000 # Read minimum length (rule contig_reformating_and_stats in 01_preprocessing.smk)
ENTROPY: 0.5 # Read minimum entropy (rule remove_low_quality in 01_preprocessing.smk)
ENTROPYWINDOW: 25 # entropy window for low qual read filter

# CLUSTER READS TO SEQTABLE (MMSEQS EASY-LINCLUST)
 # -c = req coverage of target seq
 # --min-seq-id = req identity [0-1] of alignment
linclustParams:
 --kmer-per-seq-scale 0.3
 -c 0.7
 --cov-mode 1
 --min-seq-id 0.95
 --alignment-mode 3
```

There are additional settings further down in the config file for users that are familiar with MMSeqs, 
as well as some settings that you should not alter.

## Alignment filtering

Hecatomb has settings for filtering MMSeqs alignments at each stage of the search strategy.
By default, we use a lenient e-value cutoff to maximise the identification of viral sequences in the primary searches,
and a more stringent e-value cutoff for the multi-kingdom search.
You can lower the evalue cutoffs (`-e`) to improve runtime performance at the cost of lower recall.
The `--min-lenghth` should be the same or lower than the preprocessing cutoffs.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# ALIGNMENT FILTERING CUTOFFS
  # --min-length for AA should be equal or less than 1/3 of READ_MINLENGTH
  # --min-length for NT should be equal or less than READ_MINLENGTH
filtAAprimary:
 --min-length 30
 -e 1e-3
filtAAsecondary:
 --min-length 30
 -e 1e-5
filtNTprimary:
 --min-length 90
 -e 1e-3
filtNTsecondary:
 --min-length 90
 -e 1e-5
```

## Alignment settings

Hecatomb can perform MMSeqs alignments using either sensitive (default) or fast (`--fast`) parameters.
You can tweak the setting in the config file but you should consult the MMSeqs documentation before making any changes.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# PERFORMANCE SETTINGS - SEE MMSEQS DOCUMENTATION FOR DETAILS
# sensitive AA search
perfAA:
 --start-sens 1
 --sens-steps 3
 -s 7
 --lca-mode 2
 --shuffle 0
# fast AA search
perfAAfast:
 -s 4.0
 --lca-mode 1
 --shuffle 0
# sensitive NT search
perfNT:
 --start-sens 2
 -s 7
 --sens-steps 3
# fast NT search
perfNTfast:
 -s 4.0
```

## Additional Snakemake commands

As mentioned, Hecatomb is powered by Snakemake but runs via a launcher for your convenience.
The launcher--called with `hecatomb`--lets you specify the directory with your reads, host genome, where to save the results,
whether to do an assembly, and either specify the number of threads to use or a profile to use.
Snakemake itself has many command line options and the launcher can pass additional commands to Snakemake using the `--snake` option.

One such example is if you're not production ready you might wish to do a 'dry-run', where the run is simulated but no 
jobs are submitted, just to see if everything is configured correctly.
To do that, Snakemake needs the dry run flag (`--dry-run`, `--dryrun`, or `-n`).
In Hecatomb, you can pass this flag like so:

```bash
hecatomb run --reads fasq/ --profile slurm --snake=--dry-run
```

Hecatomb prints the Snakemake command to the terminal window before running and you should see these additional options 
added to the Snakemake command:

```bash
hecatomb run --test --snake=--dry-run
```
```text
Running Hecatomb
Running snakemake command:
snakemake -j 32 --use-conda --conda-frontend mamba --rerun-incomplete --printshellcmds \
  --nolock --conda-prefix /scratch/hecatomb/snakemake/workflow/conda \
  --dry-run -s /scratch/hecatomb/snakemake/workflow/Hecatomb.smk \
  -C Reads=/scratch/hecatomb/test_data Host=human Output=hecatomb_out Assembly=False
Building DAG of jobs...
```

Have a look at the full list of available Snakemake options with `snakemake --help`. 
The launcher will pass anything in `--snake=` verbatim, so use with care.

**NOTE: Don't use hecatomb's --snake arg to pass additional config settings with Snakemake's -C/--config arg**,
instead, follow the directions above [for changing the config of a Hecatomb run](configuration.md#changing-the-hecatomb-configuration).
Read counts for all samples following Step 05: removal of vector contaminants (PhiX + NCBI UniVecDB).
These are the hits for the Secondary viral NT search.
Note: total read counts, not rep seq counts.These are the hits for the Secondary multi-kingdom AA search.
Note: total read counts, not rep seq counts.These are the hits for the Primary viral AA search.
Note: total read counts, not rep seq counts.Read counts for all samples following Step 09: clustering similiar sequences.
Note: only using R1 reads.These are the hits for the Primary viral NT search.
Note: total read counts, not rep seq counts.This is a sankey diagram visualising the breakdown of read preprocessing and taxonomic assignments.Read counts for all samples following Step 03: removal of primer-free adapter (5' and 3').
Read counts for all samples following Step 02: removal of 3' read-through contamination.
Read counts for all samples following Step 04: removal of adapter-free primer (5' and 3').Table summarising the read counts at all taxon levels per sample.
Merge this table with your metadata and use for preliminary statistical analysis.
Read counts for all samples before ANY read preprocessing.
Read counts for all samples following Step 06: removal of low quality bases and short reads.
Read counts for all samples following Step 01: removal of 5' primer.
Read counts for all samples following Step 07: removal of host contamination.
Read counts for all samples following Step 08: removal of exact duplicates.
Note: we only use R1 reads from here.