# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# FAQs

- [Why do I get FTP errors when I download hundreds or thousands of genomes (e.g. bacteria or viruses) with biomartr?](https://github.com/HajkD/biomartr/issues/12)
    - > Short Answer: NCBI and Ensembl limit the server access for individual connections. Hence, when hundreds of retrieval queries are sent to download hundreds of genomes, it may happen that some of the queries get cancelled after a certain count. In that case, simply re-run the corresponding `meta.retrieval` function and it will pick up where it left off. 
- [How is the `biomartr` package different from the `BiomaRt` package that I have been using so far?](https://github.com/HajkD/biomartr/issues/11)
    - > Short Answer: `biomartr` is designed for the retrieval of big datasets such as multiple genomes, proteomes, annotation, etc whereas `biomaRt` is designed for the retrieval of small data such as sequences of individual genes 
- [How can I retrieve GO terms for human genes using ensembl gene IDs?](https://github.com/HajkD/biomartr/issues/5)
- [How can I download coding sequences or protein sequences for all sequenced species?](https://www.biostars.org/p/9202/#235962)
- [Where can I Download the human reference genome in fasta format?](https://www.biostars.org/p/1796/#236039)
biomartr
========

<!-- badges: start -->
[![](https://badges.ropensci.org/93_status.svg)](https://github.com/ropensci/onboarding/issues/93)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/biomartr.svg?branch=master)](https://travis-ci.org/ropensci/biomartr)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/biomartr)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/biomartr)](https://github.com/metacran/cranlogs.app)
[![Paper link](https://img.shields.io/badge/Published%20in-Bioinformatics-126888.svg)](https://academic.oup.com/bioinformatics/article/33/8/1216/2931816)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-biomartr/README.html)
<!-- badges: end -->
 
## Genomic Data Retrieval with R

### Motivation:

This package is born out of my own frustration to automate the genomic data retrieval process to create computationally reproducible scripts for large-scale genomics studies. Since I couldn't find easy-to-use and fully reproducible software libraries I sat down and tried to implement a framework that would enable anyone to automate and standardize the genomic data retrieval process. I hope that this package is useful to others as well and that it helps to promote reproducible research in genomics studies.

I happily welcome anyone who wishes to contribute to this project :) Just drop me an email.

Please find a detailed [documentation here](https://docs.ropensci.org/biomartr/articles/).


### Citation

__Please cite `biomartr` if it was helpful for your research. This will allow me to
continue maintaining this project in the future.__

> Drost HG, Paszkowski J. __Biomartr: genomic data retrieval with R__. *Bioinformatics* (2017) 33(8): 1216-1217. [doi:10.1093/bioinformatics/btw821](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btw821/2931816/Biomartr-genomic-data-retrieval-with-R).


### Short package description:

The vastly growing number of sequenced genomes allows us to perform a new type of biological research.
Using a comparative approach these genomes provide us with new insights on how biological information is encoded 
on the molecular level and how this information changes over evolutionary time.

The first step, however, of any genome based study is to retrieve genomes and their annotation from databases. To automate the
retrieval process of this information on a meta-genomic scale, the `biomartr` package provides interface functions for genomic sequence retrieval and functional annotation retrieval. The major aim of `biomartr` is to facilitate computational reproducibility and large-scale handling of genomic data for (meta-)genomic analyses.
In addition, `biomartr` aims to address the `genome version crisis`. With `biomartr` users can now control and be informed 
about the genome versions they retrieve automatically. Many large scale genomics studies lack this information
and thus, reproducibility and data interpretation become nearly impossible when documentation of genome version information
gets neglected.

In detail, `biomartr` automates genome, proteome, CDS, RNA, Repeats, GFF/GTF (annotation), genome assembly quality, and metagenome project data retrieval from the major biological databases such as

- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
- [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/)
- [ENSEMBL](https://www.ensembl.org/index.html)
- [ENSEMBLGENOMES](http://ensemblgenomes.org) (as of April 2019 - `ENSEMBL` and `ENSEMBLGENOMES` were joined - see [details here](http://www.ensembl.info/2019/03/08/joint-rest-server-for-ensembl-and-ensembl-genomes-in-ensembl-96/))
- [UniProt](http://www.uniprot.org)

Furthermore, an interface to the `Ensembl Biomart` database allows users to retrieve functional annotation for genomic loci using a novel and organism centric search strategy. In addition, users can [download entire databases](https://github.com/HajkD/biomartr/blob/master/vignettes/Database_Retrieval.Rmd) such as 

- `NCBI RefSeq` 
- `NCBI nr` 
- `NCBI nt`
- `NCBI Genbank`
- `ENSEMBL` 

with only one command.

### Similar Work

The main difference between the [BiomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html) package and the [biomartr](https://docs.ropensci.org/biomartr/) package is that `biomartr` extends the `functional annotation retrieval` procedure of `BiomaRt` and __in addition__ provides useful retrieval functions for genomes, proteomes, coding sequences, gff files, RNA sequences, Repeat Masker annotations files, and functions for the retrieval of entire databases such as `NCBI nr` etc.

Please consult the [Tutorials section](https://docs.ropensci.org/biomartr/#tutorials) for more details.

`In the context of functional annotation retrieval` the `biomartr` package allows users to screen available marts using only the scientific name of an organism of interest instead of first searching for marts and datasets which support a particular organism of interest (which is required when using the `BiomaRt` package). Furthermore, `biomartr` allows you to search for particular topics when searching for attributes and filters. I am aware that the similar naming of the packages is unfortunate, but it arose due to historical reasons (please find a detailed explanation here: https://github.com/ropensci/biomartr/blob/master/FAQs.md and here [#11](https://github.com/ropensci/biomartr/issues/11)).

I also dedicated [an entire vignette to compare](https://docs.ropensci.org/biomartr/articles/Functional_Annotation.html) the `BiomaRt` and `biomartr` package functionality in the context of `Functional Annotation` (where their functionality overlaps which comprises about only 20% of the overall functionality of the biomartr package).

### Feedback
>__I truly value your opinion and improvement suggestions. Hence, I would be extremely grateful if you could take this 1 minute and 3 question survey (https://goo.gl/forms/Qaoxxjb1EnNSLpM02) so that I can learn how to improve `biomartr` in the best possible way. Many many thanks in advance.__


## Installation

The `biomartr` package relies on some [Bioconductor](https://www.bioconductor.org/install/) tools and thus requires
installation of the following packages:

```r
# Install core Bioconductor packages
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
# Install package dependencies
BiocManager::install("Biostrings")
BiocManager::install("biomaRt")

```

Now users can install `biomartr` from CRAN:

```r
# install biomartr 0.9.2
install.packages("biomartr", dependencies = TRUE)
```

## Installation with Bioconda

With an activated Bioconda channel (see [2. Set up channels](http://bioconda.github.io/user/install.html#set-up-channels)), install with:

```
conda install r-biomartr
```

and update with:

```
conda update r-biomartr
```

or use the docker container:
```
docker pull quay.io/biocontainers/r-biomartr:<tag>
```
(check [r-biomartr/tags](https://quay.io/repository/biocontainers/r-biomartr?tab=tags) for valid values for <tag>)


## Example

### Collection Retrieval

The automated retrieval of collections (= Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats files)
will make sure that the genome file of an organism will match the CDS, proteome, RNA, GFF, etc file
and was generated using the same genome assembly version. One aspect of why genomics studies
fail in computational and biological reproducibility is that it is not clear whether CDS, proteome, RNA, GFF, etc files
used in a proposed analysis were generated using the same genome assembly file denoting the same genome assembly version.
To avoid this seemingly trivial mistake we encourage users to retrieve
genome file collections using the `biomartr` function `getCollection()`
and attach the corresponding output as Supplementary Data
to the respective genomics study to ensure computational and biological reproducibility.


```r
# download collection for Saccharomyces cerevisiae
biomartr::getCollection( db = "refseq", organism = "Saccharomyces cerevisiae")
```

Internally, the `getCollection()` function will now generate a folder named `refseq/Collection/Saccharomyces_cerevisiae`
and will store all genome and annotation files for `Saccharomyces cerevisiae` in the same folder.
In addition, the exact genoem and annotation version will be logged in the `doc` folder.

Internally, a text file named `doc_Saccharomyces_cerevisiae_db_refseq.txt` is generated. The information stored in this log file is structured as follows:

```
File Name: Saccharomyces_cerevisiae_assembly_stats_refseq.txt
Organism Name: Saccharomyces_cerevisiae
Database: NCBI refseq
URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt
Download_Date: Wed Jun 27 15:21:51 2018
refseq_category: reference genome
assembly_accession: GCF_000146045.2
bioproject: PRJNA128
biosample: NA
taxid: 559292
infraspecific_name: strain=S288C
version_status: latest
release_type: Major
genome_rep: Full
seq_rel_date: 2014-12-17
submitter: Saccharomyces Genome Database
```

In an ideal world this reference file could then be included as supplementary information in any
life science publication that relies on genomic information so that
reproducibility of experiments and analyses becomes achievable.


### Genome retrieval of hundreds of genomes using only one command

Download all mammalian vertebrate genomes from `NCBI RefSeq` via:

```r
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "genome")
```
All geneomes are stored in the folder named according to the kingdom.
In this case `vertebrate_mammalian`. Alternatively, users can specify
the `out.folder` argument to define a custom output folder path.

### Platforms

> Find `biomartr` also at [OmicTools](https://omictools.com/biomartr-tool).

### Frequently Asked Questions (FAQs)

Please find [all FAQs here](https://github.com/ropensci/biomartr/blob/master/FAQs.md).

### Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/biomartr/issues


## Tutorials

Getting Started with `biomartr`:

- [NCBI Database Retrieval](https://docs.ropensci.org/biomartr/articles/Database_Retrieval.html)
- [Genomic Sequence Retrieval](https://docs.ropensci.org/biomartr/articles/Sequence_Retrieval.html)
- [Meta-Genome Retrieval](https://docs.ropensci.org/biomartr/articles/MetaGenome_Retrieval.html)
- [Functional Annotation](https://docs.ropensci.org/biomartr/articles/Functional_Annotation.html)
- [BioMart Examples](https://docs.ropensci.org/biomartr/articles/BioMart_Examples.html)


Users can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r
# source the biomartr package
library(biomartr)

# look for all tutorials (vignettes) available in the biomartr package
# this will open your web browser
browseVignettes("biomartr")
```

## NEWS
The current status of the package as well as a detailed history of the functionality of each version of `biomartr` can be found in the [NEWS](https://docs.ropensci.org/biomartr/news/index.html) section.


## Install Developer Version
Some bug fixes or new functionality will not be available on CRAN yet, but in the developer version here on GitHub. To download and install the most recent version of `biomartr` run:

```r
# install the current version of biomartr on your system
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ropensci/biomartr")
```

## Genomic Data Retrieval

#### Meta-Genome Retrieval

* `meta.retrieval()` : Perform Meta-Genome Retieval from NCBI of species belonging to the same kingdom of life or to the same taxonomic subgroup
* `meta.retrieval.all()` : Perform Meta-Genome Retieval from NCBI of the entire kingdom of life
* `getMetaGenomes()` : Retrieve metagenomes from NCBI Genbank
* `getMetaGenomeAnnotations()` : Retrieve annotation *.gff files for metagenomes from NCBI Genbank
* `listMetaGenomes()` : List available metagenomes on NCBI Genbank
* `getMetaGenomeSummary()` : Helper function to retrieve the assembly_summary.txt file from NCBI genbank metagenomes
* `clean.retrieval()`: Format meta.retrieval output

#### Genome Retrieval

* `listGenomes()` : List all genomes available on NCBI and ENSEMBL servers
* `listKingdoms()` : list the number of available species per kingdom of life on NCBI and ENSEMBL servers
* `listGroups()` : list the number of available species per group on NCBI and ENSEMBL servers
* `getKingdoms()` : Retrieve available kingdoms of life
* `getGroups()` : Retrieve available groups for a kingdom of life
* `is.genome.available()` : Check Genome Availability  NCBI and ENSEMBL servers
* `getCollection()` : Retrieve a Collection: Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats
* `getGenome()` : Download a specific genome stored on NCBI and ENSEMBL servers
* `getGenomeSet()` : Genome Retrieval of multiple species
* `getProteome()` : Download a specific proteome stored on NCBI and ENSEMBL servers
* `getProteomeSet()` : Proteome Retrieval of multiple species
* `getCDS()` : Download a specific CDS file (genome) stored on NCBI and ENSEMBL servers
* `getCDSSet()` : CDS Retrieval of multiple species
* `getRNA()` : Download a specific RNA file stored on NCBI and ENSEMBL servers
* `getRNASet()` : RNA Retrieval of multiple species
* `getGFF()` : Genome Annotation Retrieval from NCBI (`*.gff`) and ENSEMBL (`*.gff3`) servers
* `getGTF()` : Genome Annotation Retrieval (`*.gtf`) from ENSEMBL servers
* `getRepeatMasker() :` Repeat Masker TE Annotation Retrieval
* `getAssemblyStats()` : Genome Assembly Stats Retrieval from NCBI
* `getKingdomAssemblySummary()` : Helper function to retrieve the assembly_summary.txt files from NCBI for all kingdoms
* `getMetaGenomeSummary()` : Helper function to retrieve the assembly_summary.txt files from NCBI genbank metagenomes
* `getSummaryFile()` : Helper function to retrieve the assembly_summary.txt file from NCBI for a specific kingdom
* `getENSEMBLInfo()` : Retrieve ENSEMBL info file
* `getGENOMEREPORT()` : Retrieve GENOME_REPORTS file from NCBI

#### Import Downloaded Files
* `read_genome()` : Import genomes as Biostrings or data.table object
* `read_proteome()` : Import proteome as Biostrings or data.table object
* `read_cds()` : Import CDS as Biostrings or data.table object
* `read_gff()` : Import GFF file
* `read_rna()` : Import RNA file
* `read_rm()` : Import Repeat Masker output file
* `read_assemblystats()` : Import Genome Assembly Stats File

#### Database Retrieval

* `listNCBIDatabases()` : Retrieve a List of Available NCBI Databases for Download
* `download.database()` : Download a NCBI database to your local hard drive
* `download.database.all()` : Download a complete NCBI Database such as e.g. `NCBI nr` to your local hard drive

### BioMart Queries

* `biomart()` : Main function to query the BioMart database
* `getMarts()` : Retrieve All Available BioMart Databases
* `getDatasets()` : Retrieve All Available Datasets for a BioMart Database
* `getAttributes()` : Retrieve All Available Attributes for a Specific Dataset
* `getFilters()` : Retrieve All Available Filters for a Specific Dataset
* `organismBM()` : Function for organism specific retrieval of available BioMart marts and datasets
* `organismAttributes()` : Function for organism specific retrieval of available BioMart attributes
* `organismFilters()` : Function for organism specific retrieval of available BioMart filters

### Performing Gene Ontology queries

#### Gene Ontology

* `getGO()` : Function to retrieve GO terms for a given set of genes


## Download Developer Version On Windows Systems

```r
# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/biomartr", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/biomartr", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("biomartr", lib.loc = "C:/Program Files/R/R-3.1.1/library")
```

### Troubleshooting on Windows Machines

- Install `biomartr` on a Win 8 laptop: [solution](https://github.com/HajkD/orthologr/issues/1) ( Thanks to Andres Romanowski )

# Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/biomartr/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.


# biomartr 1.0.2

### New Functions

- New function `check_annotation_biomartr()` helps to check whether downloaded GFF or GTF files are corrupt. Find more details [here](https://github.com/lawremi/rtracklayer/issues/15)

- new function `getCollectionSet()` allows users to retrieve a Collection: Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats of multiple species

Example:

```r
# define scientific names of species for which
# collections shall be retrieved
organism_list <- c("Arabidopsis thaliana", 
                   "Arabidopsis lyrata", 
                   "Capsella rubella")
# download the collection of Arabidopsis thaliana from refseq
# and store the corresponding genome file in '_ncbi_downloads/collection'
 getCollectionSet( db       = "refseq", 
             organism = organism_list, 
             path = "set_collections")
```

### New Features 

- the `getGFF()` function receives a new argument `remove_annotation_outliers` to enable users to remove corrupt lines from a GFF file
Example:

```r
Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana", remove_annotation_outliers = TRUE)
```

- the `getGFFSet()` function receives a new argument `remove_annotation_outliers` to enable users to remove corrupt lines from a GFF file

- the `getGTF()` function receives a new argument `remove_annotation_outliers` to enable users to remove corrupt lines from a GTF file

- adding a new message system to `biomartr::organismBM()`, `biomartr::organismAttributes()`, and `biomartr::organismFilters()` so that large API queries don't seem so unresponsive

- `getCollection()` receives new arguments `release`, `remove_annotation_outliers`, and `gunzip` that will now be passed on to downstream retrieval functions

- the `getGTF()`, `getGenome()` and `getGenomeSet()` functions receives a new argument `assembly_type = "toplevel"` to enable users to choose between toplevel and primary assembly when using ensembl database. Setting `assembly_type = "primary_assembly"` will save a lot a space on hard drives for people using large ensembl genomes.

- all `get*()` functions with `release` argument now check if the ENSEMBL release is >45 (Many thanks to @Roleren #31 #61) 

- in all `get*()` functions, the `readr::write_tsv(path = )` was exchanged to `readr::write_tsv(file = )`, since the `readr` package version > 1.4.0 is depreciating the `path` argument. 

- `tbl_df()` was deprecated in dplyr 1.0.0.
Please use `tibble::as_tibble()` instead. -> adjusted `organismBM()` accordingly

- `custom_download()`, `getGENOMEREPORT()`, and other download functions now have specified `withr::local_options(timeout = max(30000000, getOption("timeout")))` which extends the default 60sec timeout to 30000000sec



### Bug Fixes

- Fixing bug where genome availability check in `getCollection()` was only performed in `NCBI RefSeq` and not in other databases due to a constant used in `is.genome.available()` rather than a variable (Many thanks to Takahiro Yamada for catching the bug) #53

- fixing an issue that caused the `read_cds()` function to fail in `data.table` mode (Many thanks to Clement Kent) #57

- fixing an `SSL` bug that was found on `Ubuntu 20.04` systems #66 (Many thanks to Håkon Tjeldnes)

- fixing global variable issue that caused `clean.retrieval()` to fail when no documentation file was in a `meta.retrieval()` folder

- The NCBI recently started adding `NA` values as FTP file paths in their `species summary files` for species without reference genomes. As a result `meta.retrieval()` stopped working, because no FTP paths were found for some species. This issue was now fixed by adding the filter rule `!is.na(ftp_path)` into all `get*()` functions (Many thanks for making me aware of this issue Ashok Kumar Sharma #34 and Dominik Merges #72) 

- Fixing an issue in `custom_download()` where the `method` argument was causing issues when downloading from `https` directed `ftp` sites (Many thanks to @cmatKhan) #76 

- Fixing issue when trying to combine multiple summary-stats files where NA's were present in the list item that was passed along for combination in `meta.retrieval()` #73 (Many thanks to Dominik Merges)

[biomartr 0.9.2](https://github.com/ropensci/biomartr/releases/tag/v0.9.1)
- minor changes to comply with CRAN policy regarding Internet access failure 
-> Instead of using warnings or error messages, only gentle messages are allowed to be used



[biomartr 0.9.0](https://github.com/ropensci/biomartr/releases/tag/v0.9.0)
===========

__Please be aware that as of April 2019, ENSEMBLGENOMES
was retired ([see details here](http://www.ensembl.info/2019/03/08/joint-rest-server-for-ensembl-and-ensembl-genomes-in-ensembl-96/)). Hence, all `biomartr` functions were updated
and won't support data retrieval from `ENSEMBLGENOMES` servers anymore.__

### New Functions

- New function `clean.retrieval()` enables formatting and automatic unzipping of meta.retrieval output (find out more here: https://docs.ropensci.org/biomartr/articles/MetaGenome_Retrieval.html#un-zipping-downloaded-files)
- New function `getGenomeSet()` allows users to easily retrieve genomes of multiple specified species. 
In addition, the genome summary statistics for all retrieved species will be stored as well to provide
users with insights regarding the genome assembly quality of each species. This file can be used as Supplementary Information file
in publications to facilitate reproducible research.
- New function `getProteomeSet()` allows users to easily retrieve proteomes of multiple specified species
- New function `getCDSSet()` allows users to easily retrieve coding sequences of multiple specified species
- New function `getGFFSet()` allows users to easily retrieve GFF annotation files of multiple specified species
- New function `getRNASet()` allows users to easily retrieve RNA sequences of multiple specified species
- New function `summary_genome()` allows users to retrieve summary statistics for a genome assembly file to assess 
the influence of genome assembly qualities when performing comparative genomics tasks
- New function `summary_cds()` allows users to retrieve summary statistics for a coding sequence (CDS) file.
We noticed, that many CDS files stored in NCBI or ENSEMBL databases contain sequences that aren't divisible by 3 (division into codons).
This makes it difficult to divide CDS into codons for e.g. codon alignments or translation into protein sequences. In
addition, some CDS files contain a significant amount of sequences that do not start with AUG (start codon).
This function enables users to quantify how many of these sequences exist in a downloaded CDS file to process
these files according to the analyses at hand.


### New Features of Existing Functions 

- the default value of argument `reference` in `meta.retrieval()` changed from `reference = TRUE` to `reference = FALSE`.
This way all genomes (reference AND non-reference) genomes will be downloaded by default. This is what users seem to prefer.
- `getCollection()` now also retrieves `GTF` files when `db = 'ensembl'`
- `getAssemblyStats()` now also performs md5 checksum test
- all md5 checksum tests now retrieve the new md5checkfile format from NCBI RefSeq and Genbank
- `getGTF()`: users can now specify the NCBI Taxonomy ID or Accession ID in addition to the scientific name in argument 'organism' to retrieve genome assemblies 
- `getGFF()`: users can now specify the NCBI Taxonomy ID or Accession ID for ENSEMBL in addition to the scientific name in argument 'organism' to retrieve genome assemblies 
- `getMarts()` will now throw an error when BioMart servers cannot be reached (#36)
- `getGenome()` now also stores the genome summary statistics (see `?summary_genome()`) for the retrieved species in the `documentation` folder to provide
users with insights regarding the genome assembly quality
- In all get*() functions the default for argument `reference` is now set from `reference = TRUE` to `reference = FALSE` (= new default)
- all `get*()` functions now received a new argument `release` which allows users to retrieve
specific release versions of genomes, proteomes, etc from `ENSEMBL` and `ENSEMBLGENOMES`
- all `get*()` functions received two new arguments `clean_retrieval` and  `gunzip` which
allows users to upzip the downloaded files directly in the `get*()` function call and rename
the file for more convenient downstream analyses


[biomartr 0.8.0](https://github.com/ropensci/biomartr/releases/tag/v0.8.0)
===========

### New Functions

- new function `getCollection()` for retrieval of a collection: the genome sequence,
protein sequences, gff files, etc for a particular species

### New Functionality of Existing Functions 

- `getProteome()` can now retrieve proteomes from the [UniProt](http://www.uniprot.org/) database by specifying `getProteome(db = "uniprot")`.
An example can be found [here](https://github.com/ropensci/biomartr/blob/master/vignettes/Sequence_Retrieval.Rmd#example-retrieval-uniprot)

- `is.genome.available()` now prints out more useful interactive messages when searching for available organisms 

- `is.genome.available()` can now handle `taxids` and `assembly_accession ids` in addition to the scientific name when
specifying argument `organism`
An example can be found [here](https://github.com/ropensci/biomartr/blob/master/vignettes/Sequence_Retrieval.Rmd#example-ncbi-refseq)

- `is.genome.available()` can now check for organism availability in the UniProt database

- `getGenome()`: users can now specify the NCBI Taxonomy ID or Accession ID in addition to the scientific name in argument 'organism' to retrieve genome assemblies 

- `getProteome()`: users can now specify the NCBI Taxonomy ID or Accession ID in addition to the scientific name in argument 'organism' to retrieve proteomes 

- `getCDS()`: users can now specify the NCBI Taxonomy ID or Accession ID in addition to the scientific name in argument 'organism' to retrieve CDS  

- `getRNA()`: users can now specify the NCBI Taxonomy ID or Accession ID in addition to the scientific name in argument 'organism' to retrieve RNAs 

- `is.genome.available()`: argument order was changed from is.genome.available(organism, details, db) to is.genome.available(db, organism, details) to be logically more consistent
with all `get*()` functions
- `meta.retrieval` receives a new argument `restart_at_last` to indicate whether or not the download process when re-running the `meta.retrieval` function
shall pick up at the last species or whether it should crawl through all existing files to check the md5checksum
- `meta.retrieval` now generates an csv overview file in the `doc` folder which stores genome version, date, origin, etc information for
all downloaded organisms and can be directly used as Supplementary Data file in publications to increase computational and biological reproducibility of the genomics study
- `download.database.all()` can now skip already downloaded files and internally removes corrupted files with non-matching md5checksum. Re-downloading of currupted
files and be performed by simply re-running the `download.database.all()` command


[biomartr 0.7.0](https://github.com/ropensci/biomartr/releases/tag/v0.7.0)
===========

### Function changes

- the function `meta.retrieval()` will now pick up the download at the organism
where it left off and will report which species have already been retrieved 

- all `get*()` functions and the `meta.retrieval()` function receive a new argument `reference` which allows users to retrieve not-reference or not-representative genome versions when downloading from NCBI RefSeq or NCBI Genbank

- the argument order in `meta.retrieval()` changed from `meta.retrieval(kingdom, group, db, ...)` to `meta.retrieval(db,kingdom, group, ...)` to make the argument order more consistent with the `get*()` functions

- the argument order in `getGroups()` changed from `getGroups(kingdom, db)` to `getGroups(db, kingdom)` to make the argument order more consistent with the `get*()` and `meta.retrieval()` functions


### New Functions

- new internal functions `existingOrganisms()` and `existingOrganisms_ensembl()`
which check the organisms that have already been downloaded

biomartr 0.5.2
===========

### Bug fixes

- fixing bug (https://github.com/ropensci/biomartr/issues/6) that caused incorrect filtering condition when more than one entry for an organism is present in the assemblysummary.txt file at NCBI (Thanks to @kalmeshv)


biomartr 0.5.1
===========

### Bug fixes

- fixing a bug in `exists.ftp.file()` and `getENSEMBLGENOMES.Seq()` that caused bacterial genome, proteome, etc retrieval to fail due to the wrong construction of a query ftp request https://github.com/HajkD/biomartr/issues/7
(Many thanks to @dbsseven)

- fix a major bug in which organisms having no representative genome would generate NULL paths that subsequently crashed the `meta.retrieval()` function when it tried to print out the result paths.

### New Functions

- new function `getRepeatMasker()` for retrieval of Repeat Masker output files  

- new function `getGTF()` for genome annotation retrieval from `ensembl` and `ensemblgenomes` in `gtf` format (Thanks for suggesting it Ge Tan)

- new function `getRNA()` to perform RNA Sequence Retrieval from NCBI and ENSEMBL databases (Thanks for suggesting it @carlo-berg)

- new function `read_rna()` for importing Repeat Masker output files downloaded with `getRepeatMasker()`

- new function `read_rm()` for importing RNA downloaded with `getRNA()` as Biostrings or data.table object

- new helper function `custom_download()` that aims to make the download process more robust and stable
-> In detail, the download process is now adapting to the operating system, e.g. using either `curl` (macOS), `wget` (Linux), or `wininet` (Windows)




### Function changes

- function name `listDatabases()` has been renamed `listNCBIDatabases()`. In `biomartr` version 0.6.0 the function name `listDatabases()` will be depreciated

- `meta.retieval()` and `meta.retieval.all()` now allow the bulk retrieval of GTF files for `type = 'ensembl'` and `type = 'esnemblgenomes'` via `type = "gtf"`. See `getGTF()` for more details.

- `meta.retieval()` and `meta.retieval.all()` now allow the bulk retrieval of RNA files via `type = "rna"`. See `getRNA()` for more details.

- `meta.retieval()` and `meta.retieval.all()` now allow the bulk retrieval of Repeat Masker output files via `type = "rm"`. See `getRepeatMasker()` for more details.

- all `get*()` retrieval functions now skip the download of a particular file if it already exists in the specified file path

- `download.database()` and `download.database.all()` now internally perform md5 check sum checks to make sure that the file download was successful

- `download.database()` and `download.database.all()` now return the file paths of the downloaded file so that it is easier to use these
functions when constructing pipelines, e.g. `download.database() %>% ...` or `download.database.all() %>% ...`.

- `meta.retrieval()` and `meta.retrieval.all()` now return the file paths of the downloaded file so that it is easier to use these
functions when constructing pipelines, e.g. `meta.retrieval() %>% ...` or `meta.retrieval() %>% ...`.

- `getGenome()`, `getProteome()`, `getCDS()`, `getRNA()`, `getGFF()`, and `getAssemblyStats()` now internally perform md5 checksum tests
to make sure that files are retrieved intact.


biomartr 0.4.0
===========

### Bug fixes

- fixing a major bug https://github.com/HajkD/biomartr/issues/6 that caused that in all `get*()` (genome, proteome, gff, etc.) and `meta.retrieval*()` functions
 the meta retrieval process errored and terminated whenever NCBI or ENSEMBL didn't
store all types of sequences for a particular organism: genome, proteome, cds, etc. This has been fixed now and function calls
such as `meta.retrieval(kingdom = "bacteria", db = "genbank", type = "proteome")` should work properly now (Thanks to @ARamesh123 for making me aware if this bug). Hence, this bug affected all attempts to download all proteome sequences e.g. for bacteria and viruses, because NCBI does not store genome AND proteome information for all bacterial or viral species. 


### New Functions

- new function `getAssemblyStats()` allows users to retrieve the genome assembly stats file from NCBI RefSeq or Genbank, e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_stats.txt

- new function `read_assemblystats()` allows to import the genome assembly stats file from NCBI RefSeq or Genbank that was retrieved
using the `getAssemblyStats()` function

### Function changes

- `meta.retrieval()` and `meta.retrieval.all()` can now also download genome assembly stats for all selected species

- `meta.retrieval()` receives a new argument `group` that allows users to retrieve species belonging to a subgroup instead of the entire kingdom.
Available groups can be retrieved with `getGroups()`.

- functions `getSubgroups()` and `listSubgroups()` have been removed and their initial functionality
has been merged and integrated into `getGroups()` and `listGroups()`

- `listGroups()` receives a new argument `details` that allows users to retrieve the organism names that belong to the corresponding subgroups

- `getGroups()` is now based on `listGroups()`

- internal function `getGENOMESREPORT()` is now exported and available to the user

- all `organism*()` functions now also support Ensembl Plants, Ensembl Metazoa, Ensembl Protist, and Ensembl Fungi (Thanks for pointing out [Alex Gabel](https://github.com/AlexGa))

- `getMarts()` and `getDatasets()` now also support Ensembl Plants, Ensembl Metazoa, Ensembl Protist, and Ensembl Fungi (Thanks for pointing out [Alex Gabel](https://github.com/AlexGa))


### Vignette updates

- Vignette `Meta-Genome Retrieval` has more examples how to download genomes of species that belong to the same subgroup


biomartr 0.3.0
===========

### Bug fixes

- Fixing a bug https://github.com/HajkD/biomartr/issues/2 based on the [readr package](https://github.com/tidyverse/readr) that affected the `getSummaryFile()`, `getKingdomAssemblySummary()`, `getMetaGenomeSummary()`,
`getENSEMBL.Seq()` and `getENSEMBLGENOMES.Seq()` functions causing quoted lines in the `assembly_summary.txt` to be omitted when reading these files. This artefact caused that e.g. instead of information of 80,000 Bacteria genomes only 40,000 (which non-quotations) were read (Thanks to [Xin Wu](https://github.com/alartin)).


biomartr 0.2.1
===========

In this version of `biomartr` the `organism*()` functions were adapted to the new [ENSEMBL 87 release](http://www.ensembl.info/blog/2016/12/08/ensembl-87-has-been-released/)
in which organism name specification in the Biomart description column [was changed](https://github.com/HajkD/biomartr/issues/1)
from a scientific name convention to a mix of common name and scientific name convention.

- all `organism*()` functions have been adapted to the new ENSEMBL 87 release organism name notation that is used in the Biomart description

- fixing error handling bug that caused commands such as `download.database(db = "nr.27.tar.gz")` to not execute properly

biomartr 0.2.0
===========

In this version, `biomartr` was extended to now retrieve genome, proteome, CDS, GFF and meta-genome data
also from [ENSEMBL](http://www.ensembl.org/index.html) and [ENSEMLGENOMES](http://ensemblgenomes.org/).
Furthermore, all NCBI retrieval functions were updated to the new server folder structure standards of NCBI.


### New Functions

- new meta-retrieval function `meta.retrieval.all()` allows users to download all individual genomes of all kingdoms of life with one command

- new metagenome retrieval function `getMetaGenomes()` allows users to retrieve metagenome projects from NCBI Genbank

- new metagenome retrieval function `getMetaGenomeAnnotations()` allows users to retrieve annotation files for genomes belonging to a metagenome project stored at NCBI Genbank

- new retrieval function `getGFF()` allows users to retrieve annotation (*.gff) files for specific genomes from NCBI and ENSEMBL databases

- new import function `read_gff()` allowing users to import GFF files downloaded with `getGFF()`

- new internal functions to check for availability of ENSEMBL or ENSEMBLGENOMES databases

- new database retrieval function `download.database.all()` allows users to download entire NCBI databases with one command

- new function `listMetaGenomes()` allowing users to list available metagenomes on NCBI Genbank

- new external helper function `getSummaryFile()` to retrieve the assembly_summary.txt file from NCBI

- new external helper function `getKingdomAssemblySummary()` to retrieve the assembly_summary.txt files from NCBI for all kingdoms and combine them
into one big data.frame

- new function `listKingdoms()` allows users to list the number of available species per kingdom of life

- new function `listGroups()` allows users to list the number of available species per group

- new function `listSubgroups()` allows users to list the number of available species per subgroup

- new function `getGroups()` allows users to retrieve available groups for a kingdom of life

- new function `getSubgroups()` allows users to retrieve available subgroups for a kingdom of life

- new external helper function `getMetaGenomeSummary()` to retrieve the assembly_summary.txt files from NCBI genbank metagenomes

- new internal helper function `getENSEMBL.Seq()` acting as main interface function to communicate with the ENSEMBL database API for sequence retrieval

- new internal helper function `getENSEMBLGENOMES.Seq()` acting as main interface function to communicate with the ENSEMBL database API for sequence retrieval

- new internal helper function `getENSEMBL.Annotation()` acting as main interface function to communicate with the ENSEMBL database API for GFF retrieval

- new internal helper function `getENSEMBLGENOMES.Annotation()` acting as main interface function to communicate with the ENSEMBL database API for GFF retrieval

- new internal helper function `get.ensemblgenome.info()` to retrieve general organism information from ENSEMBLGENOMES 

- new internal helper function `get.ensembl.info()` to retrieve general organism information from ENSEMBL 

- new internal helper function `getGENOMEREPORT()` to retrieve the genome reports file from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt

- new internal helper function `connected.to.internet()` enabling internet connection check

### Function changes

- functions `getGenome()`, `getProteome()`, and `getCDS()` now can also in addition to NCBI retrieve genomes, proteomes or CDS from  [ENSEMBL](http://www.ensembl.org/index.html) and [ENSEMLGENOMES](http://ensemblgenomes.org/)

- the functions `getGenome()`, `getProteome()`, and `getCDS()` were completely re-written and now use the assembly_summary.txt files
provided by NCBI to retrieve the download path to the corresponding genome. Furthermore, these functions now lost the `kingdom` argument.
Users now only need to specify the organism name and not the kingdom anymore. Furthermore, all `get*` functions now
return the path to the downloaded genome so that this path can be used as input to all `read_*` functions.

- `download_databases()` has been renamed to `download.databases()` to be more consistent with other function notation

- the argument `db_format` was removed from `listDatabases()` and `download.database()` because it was misleading

- the command `listDatabases("all")` now returns all available NCBI databases that can be retrieved with `download.database()`

- `download.database()` now internally checks if input database specified by the user is actually available on NCBI servers

- the documentary file generated by `getGenome()`, `getProteome()`, and `getCDS()` is now extended to store more details about the downloaded genome

- argument `database` in `is.genome.available()` and `listGenomes()` has been renamed to `db` to be consistent with all other sequence retrieval functions

- `is.genome.available()` now also checks availability of organisms in ENSEMBL. See `db = "ensembl"`

- the argument `db_name` in `listDatabases()` has been renamed `db` to be more consistent with the notation in other functions

- the argument `name` in `download.database()` has been renamed `db` to be more consistent with the notation in other functions

- `getKingdoms()` now retrieves also kingdom information for ENSEMBL and ENSEMBLGENOMES

- `getKingdoms()` received new argument `db` to specify from which database (e.g. `refseq`, `genbank`, `ensembl` or `ensemblgenomes`) kingdom information shall be retrieved

- `getKingdoms(db = "refseq")` received one more member: `"viral"`, allowing the genome retrieval of all viruses

- argument `out.folder` in `meta.retrieval()` has been renamed to `path` to be more consistent with other retrieval functions

- all `read_*` functions now received a new argument `obj.type` allowing users to choose between storing input genomes as Biostrings object or data.table object

- all `read_*` functions now have `format = "fasta"` as default

- the `kingdom` argument in the `listGenomes()` function was renamed to `type`, now allowing users to specify not only specify kingdoms,
but also groups and subgroups. Use: `listGenomes(type = "kingdom")` or `listGenomes(type = "group")` or `listGenomes(type = "subgroup")`

- the `listGenomes()` function receives a new argument `subset` to specify a subset of the selected `type` argument. E.g. `subset = "Eukaryota"` when specifying
`type = "kingdom"`


### Vignette updates

- new Vignette `Meta-Genome Retrieval`
- Update examples and extend `Introduction` Vignette 
- Update examples and extend `Database Retrieval` Vignette 
- Update examples and extend `Sequence Retrieval` Vignette
- Update examples and extend `Functional Annotation` Vignette


biomartr 0.1.0
===========

- fixing a parsing error of the file `ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt`
The problem was that comment lines were introduced and columns couldn't be parsed correctly anymore. This caused that genomes, proteomes, and CDS files could not be downloaded properly. This has been fixed now.

- genomes, proteome, and CDS as well as meta-genomes can now be retrieved
from RefSeq and Genbank (not only RefSeq); only `getCDS()` does not have genebank access,
becasue genbank does not provide CDS sequences

- adding new function `meta.retrieval()` to mass retrieve genomes for entire kingdoms of life 

- fixed a major bug in `organismBM()` causing the function to fail. The failure of
this function affected all downstream `organism*()` functions. Bug is now fixed and everything
works properly

- updated Vignettes

biomartr 0.0.3
===========

- updating unit tests for new API

- fixing API problems that caused all BioMart related functions to fail

- fixing retrieval problems in `getCDS()`, `getProteome()`, and `getGenome()`

- the `listDatabases()` function now has a new option `db_name = "all"` allowing users to list all available databases stored on NCBI 

### Vignettes
- adding new vignette: Database Retrieval
- update the vignettes: Phylotranscriptomics, Sequence Retrieval, and Functional Annotation


biomartr 0.0.2
===========

### Vignettes
- adding vignettes: Introduction, Functional Annotation, Phylotranscriptomics, and Sequence Retrieval

biomartr 0.0.1
===========

Release Version
