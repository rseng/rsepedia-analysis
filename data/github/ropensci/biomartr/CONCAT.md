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
---
title: "Functional Annotation with biomartr"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Functional Annotation with biomartr}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```


## Functional Annotation Retrieval from `Ensembl Biomart`

> **_NOTE:_** To make sure that you have a sufficiently stable (internet) connection between R and the respective databases, please set the default `timeout` setting __on your local machine__ from 60sec to at least 30000sec before running any retrieval functions via:

```r
options(timeout = 30000)
```

### Getting Started 

The [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart) database enables users to retrieve a vast diversity of annotation data
for specific organisms. Initially, Steffen Durinck and Wolfgang Huber provided a powerful interface between
the R language and [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart) by implementing the R package [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html). 

The purpose of the `biomaRt` package was to mimic the ENSEMBL BioMart database structure to construct queries that can be sent to the Application Programming Interface (API) of BioMart. Although, this procedure was very useful in the past, it seems not intuitive from an organism centric point of view. Usually, users wish to download functional annotation for a particular organism of interest. However, the BioMart and thus the `biomaRt` package require that users already know in which `mart` and  `dataset` the organism of interest will be found which requires significant efforts of searching and screening. In addition, once the  `mart` and  `dataset` of a particular organism of interest were found and specified the user must again learn which `attribute` has to be specified to retrieve the functional annotation information of interest. 

The new functionality implemented in the `biomartr` package aims to overcome this 
search bottleneck by extending the functionality of the [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html) package. The new `biomartr` package introduces a more organism cantered annotation retrieval concept which does not require to screen for `marts`, `datasets`, and `attributes` beforehand. With `biomartr` users only need to specify the `scientific name` of the organism of interest to then retrieve available `marts`, `datasets`, and `attributes` for the corresponding organism of interest.   

This paradigm shift enables users to quickly construct queries to the BioMart database without having to learn the particular database structure and organization of BioMart. 

The following sections will introduce users to the functionality and data retrieval precedures of `biomartr` and will show how `biomartr`
extends the functionality of the initial [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html) package.

### The old `biomaRt` query methodology		

The best way to get started with the _old_ methodology presented by the established [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html) package is to understand the workflow of its data retrieval process. The query logic of the `biomaRt` package derives from the database organization of [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart) which stores a vast diversity of annotation data
for specific organisms. In detail, the [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart) database is organized into so called:		
 `marts`, `datasets`, and `attributes`. `Marts` denote a higher level category of functional annotation such as `SNP` (e.g. for functional annotation of particular single nucleotide polymorphisms (SNPs)) or `FUNCGEN` (e.g. for functional annotation of regulatory regions or relationsships of genes). `Datasets` denote the
particular species of interest for which functional annotation is available __within__ this specific `mart`. It can happen that
`datasets` (= particular species of interest) are available in one `mart` (= higher category of functional annotation) but not in an other `mart`.
For the actual retrieval of functional annotation information users must then specify the `type` of functional annotation information they 
wish to retrieve. These `types` are called `attributes` in the `biomaRt` notation.
 
Hence, when users wish to retrieve information for a specific organism of interest, they first need to specify a particular `mart` and `dataset` in which the information of the corresponding organism of interest can be found. Subsequently they can specify the `attributes` argument to retrieve a particular type of functional annotation (e.g. Gene Ontology terms).

The following section shall illustrate how `marts`, `datasets`, and `attributes` could be explored
using `biomaRt` before the `biomartr` package existed.
 		
The availability of `marts`, `datasets`, and `attributes` can be checked by the following functions:		

```{r,eval=FALSE}		
# install the biomaRt package		
# source("http://bioconductor.org/biocLite.R")		
# biocLite("biomaRt")		
# load biomaRt		
library(biomaRt)		
# look at top 10 databases		
head(biomaRt::listMarts(host = "https://www.ensembl.org"), 10)		
```	

Users will observe that several `marts` providing annotation for specific classes of organisms or groups of organisms are available.		

For our example, we will choose the `hsapiens_gene_ensembl` `mart` and list all available datasets that are element of this `mart`.		

```{r,eval=FALSE}		
head(biomaRt::listDatasets(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")), 10)		
```		

The `useMart()` function is a wrapper function provided by `biomaRt` to connect a selected BioMart database (`mart`) with a corresponding dataset stored within this `mart`.		
We select dataset `hsapiens_gene_ensembl` and now check for available attributes (annotation data) that can be accessed for `Homo sapiens` genes.	

```{r,eval=FALSE}		
head(biomaRt::listAttributes(biomaRt::useDataset(
                                         dataset = "hsapiens_gene_ensembl", 		
                                         mart    = useMart("ENSEMBL_MART_ENSEMBL",		
                                         host    = "https://www.ensembl.org"))), 10)		
```		
		
Please note the nested structure of this attribute query. For an attribute query procedure an additional wrapper function named `useDataset()` is needed in which `useMart()` and a corresponding dataset needs to be specified. The result is a table storing the name of available attributes for		
_Homo sapiens_ as well as a short description.		

Furthermore, users can retrieve all filters for _Homo sapiens_ that can be specified by the actual BioMart query process.		
		
```{r,eval=FALSE}		
 head(biomaRt::listFilters(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", 		
                                               mart    = useMart("ENSEMBL_MART_ENSEMBL",		
                                               host    = "https://www.ensembl.org"))), 10)		
```		
 		
After accumulating all this information, it is now possible to perform an actual BioMart query by using the `getBM()` function.		
 		
In this example we will retrieve attributes: `start_position`,`end_position` and `description`		
for the _Homo sapiens_ gene `"GUCA2A"`.		
 		
Since the input genes are `ensembl gene ids`, we need to specify the `filters` argument `filters = "hgnc_symbol"`.		
 		
```{r,eval=FALSE}		
 # 1) select a mart and data set		
 mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", 		
                    mart    = useMart("ENSEMBL_MART_ENSEMBL",		
                    host    = "https://www.ensembl.org"))		
 		
 # 2) run a biomart query using the getBM() function		
 # and specify the attributes and filter arguments		
 geneSet <- "GUCA2A"		
 		
 resultTable <- biomaRt::getBM(attributes = c("start_position","end_position","description"),		
                      filters    = "hgnc_symbol", 		
                      values     = geneSet, 		
                      mart       = mart)		
 		
 resultTable 		
```		
 		
When using `getBM()` users can pass all attributes retrieved by `listAttributes()` to the `attributes` argument of the `getBM()` function.		
  
  
## Extending `biomaRt` using the new query system of the `biomartr` package

### Getting Started with `biomartr`

This query methodology provided by `Ensembl Biomart` and the `biomaRt` package is a very well defined approach
for accurate annotation retrieval. Nevertheless, when learning this query methodology it (subjectively)
seems non-intuitive from the user perspective. Therefore, the `biomartr` package provides another
query methodology that aims to be more organism centric.

Taken together, the following workflow allows users to perform fast BioMart queries for 
attributes using the `biomart()` function implemented in this `biomartr` package:

1) get attributes, datasets, and marts via : `organismAttributes()`

2) choose available biological features (filters) via: `organismFilters()`

3) specify a set of query genes: e.g. retrieved with `getGenome()`, `getProteome()` or `getCDS()`

4) specify all arguments of the `biomart()` function using steps 1) - 3) and
perform a BioMart query

__Note that dataset names change very frequently due to the update of dataset versions.
So in case some query functions do not work properly, users should check with
`organismAttributes(update = TRUE)` whether or not their dataset name has been changed.
For example, `organismAttributes("Homo sapiens", topic = "id", update = TRUE)`
might reveal that the dataset `ENSEMBL_MART_ENSEMBL` has changed.__


## Retrieve marts, datasets, attributes, and filters with biomartr

### Retrieve Available Marts

The `getMarts()` function allows users to list all available databases that can be accessed through BioMart interfaces.

```{r,eval=FALSE}
# load the biomartr package
library(biomartr)

# list all available databases
biomartr::getMarts()
```

```
     mart                  version                       
   <chr>                 <chr>                         
 1 ENSEMBL_MART_ENSEMBL  Ensembl Genes 104             
 2 ENSEMBL_MART_MOUSE    Mouse strains 104             
 3 ENSEMBL_MART_SEQUENCE Sequence                      
 4 ENSEMBL_MART_ONTOLOGY Ontology                      
 5 ENSEMBL_MART_GENOMIC  Genomic features 104          
 6 ENSEMBL_MART_SNP      Ensembl Variation 104         
 7 ENSEMBL_MART_FUNCGEN  Ensembl Regulation 104        
 8 plants_mart           Ensembl Plants Genes 51       
 9 plants_variations     Ensembl Plants Variations 51  
10 fungi_mart            Ensembl Fungi Genes 51        
11 fungi_variations      Ensembl Fungi Variations 51   
12 protists_mart         Ensembl Protists Genes 51     
13 protists_variations   Ensembl Protists Variations 51
14 metazoa_mart          Ensembl Metazoa Genes 51      
15 metazoa_variations    Ensembl Metazoa Variations 51 
```

### Retrieve Available Datasets from a Specific Mart

Now users can select a specific database to list all available data sets that can be accessed through this database. In this example we choose
the `ENSEMBL_MART_ENSEMBL` database.

```{r,eval=FALSE}
head(biomartr::getDatasets(mart = "ENSEMBL_MART_ENSEMBL") , 5)
```

```
 dataset                 description                        version     
  <chr>                   <chr>                              <chr>       
1 fcatus_gene_ensembl     Cat genes (Felis_catus_9.0)        Felis_catus…
2 umaritimus_gene_ensembl Polar bear genes (UrsMar_1.0)      UrsMar_1.0  
3 ogarnettii_gene_ensembl Bushbaby genes (OtoGar3)           OtoGar3     
4 lcrocea_gene_ensembl    Large yellow croaker genes (L_cro… L_crocea_2.0
5 sformosus_gene_ensembl  Asian bonytongue genes (fSclFor1.… fSclFor1.1 
```

Now you can select the dataset `hsapiens_gene_ensembl` and list all available attributes that can be retrieved from this dataset.

```{r,eval=FALSE}
tail(biomartr::getDatasets(mart = "ENSEMBL_MART_ENSEMBL") , 38)
```

```
1 csabaeus_gene_ensembl    Vervet-AGM genes (ChlSab1.1)  ChlSab1.1      
 2 chircus_gene_ensembl     Goat genes (ARS1)             ARS1           
 3 mmulatta_gene_ensembl    Macaque genes (Mmul_10)       Mmul_10        
 4 mmonoceros_gene_ensembl  Narwhal genes (NGI_Narwhal_1) NGI_Narwhal_1  
 5 csemilaevis_gene_ensembl Tongue sole genes (Cse_v1.0)  Cse_v1.0       
 6 cpbellii_gene_ensembl    Painted turtle genes (Chryse… Chrysemys_pict…
 7 clanigera_gene_ensembl   Long-tailed chinchilla genes… ChiLan1.0      
 8 catys_gene_ensembl       Sooty mangabey genes (Caty_1… Caty_1.0       
 9 tguttata_gene_ensembl    Zebra finch genes (bTaeGut1_… bTaeGut1_v1.p  
10 nleucogenys_gene_ensembl Gibbon genes (Nleu_3.0)       Nleu_3.0       
# … with 28 more rows
```

### Retrieve Available Attributes from a Specific Dataset

Now that you have selected a database (`hsapiens_gene_ensembl`) and a dataset (`hsapiens_gene_ensembl`),
users can list all available attributes for this dataset using the `getAttributes()` function.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# list all available attributes for dataset: hsapiens_gene_ensembl
head( biomartr::getAttributes(mart    = "ENSEMBL_MART_ENSEMBL", 
                              dataset = "hsapiens_gene_ensembl"), 10 )
```

```
Starting retrieval of attribute information from mart ENSEMBL_MART_ENSEMBL and dataset hsapiens_gene_ensembl ...
                            name                  description
1                ensembl_gene_id               Gene stable ID
2        ensembl_gene_id_version       Gene stable ID version
3          ensembl_transcript_id         Transcript stable ID
4  ensembl_transcript_id_version Transcript stable ID version
5             ensembl_peptide_id            Protein stable ID
6     ensembl_peptide_id_version    Protein stable ID version
7                ensembl_exon_id               Exon stable ID
8                    description             Gene description
9                chromosome_name     Chromosome/scaffold name
10                start_position              Gene start (bp)
```

### Retrieve Available Filters from a Specific Dataset

Finally, the `getFilters()` function allows users to list available filters
for a specific dataset that can be used for a `biomart()` query.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# list all available filters for dataset: hsapiens_gene_ensembl
head( biomartr::getFilters(mart    = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "hsapiens_gene_ensembl"), 10 )
```

```
Starting retrieval of filters information from mart ENSEMBL_MART_ENSEMBL and dataset hsapiens_gene_ensembl ...
                 name                            description
1     chromosome_name               Chromosome/scaffold name
2               start                                  Start
3                 end                                    End
4          band_start                             Band Start
5            band_end                               Band End
6        marker_start                           Marker Start
7          marker_end                             Marker End
8       encode_region                          Encode region
9              strand                                 Strand
10 chromosomal_region e.g. 1:100:10000:-1, 1:100000:200000:1
```

## Organism Specific Retrieval of Information

In most use cases, users will work with a single or a set of model organisms. In this process they will mostly be
interested in specific annotations for this particular model organism. The `organismBM()`
function addresses this issue and provides users with an organism centric query to `marts` and `datasets`
which are available for a particular organism of interest.


__Note__ that when running the following functions for the first time, the data retrieval procedure will take some time, due to the remote access to BioMart. The corresponding result is then saved in a `*.txt` file named `_biomart/listDatasets.txt` within the `tempdir()` folder, allowing subsequent queries to be performed much faster.
The `tempdir()` folder, however, will be deleted after a new R session was established. In this case
the inital call of the subsequent functions again will take time to retrieve all organism specific data from the BioMart database.

This concept of locally storing all organism specific database linking information available in BioMart into
an internal file allows users to significantly speed up subsequent retrieval queries for that particular organism.


```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# retrieving all available datasets and biomart connections for
# a specific query organism (scientific name)
biomartr::organismBM(organism = "Homo sapiens")
```

```
 Starting retrieval of all available BioMart datasets for Homo sapiens ...
Datasets for the following marts will be retrieved:                                                         
                    mart                        version
1   ENSEMBL_MART_ENSEMBL              Ensembl Genes 104
2     ENSEMBL_MART_MOUSE              Mouse strains 104
3  ENSEMBL_MART_SEQUENCE                       Sequence
4  ENSEMBL_MART_ONTOLOGY                       Ontology
5   ENSEMBL_MART_GENOMIC           Genomic features 104
6       ENSEMBL_MART_SNP          Ensembl Variation 104
7   ENSEMBL_MART_FUNCGEN         Ensembl Regulation 104
8            plants_mart        Ensembl Plants Genes 51
9      plants_variations   Ensembl Plants Variations 51
10            fungi_mart         Ensembl Fungi Genes 51
11      fungi_variations    Ensembl Fungi Variations 51
12         protists_mart      Ensembl Protists Genes 51
13   protists_variations Ensembl Protists Variations 51
14          metazoa_mart       Ensembl Metazoa Genes 51
Processing mart ENSEMBL_MART_ENSEMBL ...
Processing mart ENSEMBL_MART_MOUSE ...
Processing mart ENSEMBL_MART_SEQUENCE ...
Processing mart ENSEMBL_MART_ONTOLOGY ...
Processing mart ENSEMBL_MART_GENOMIC ...
Processing mart ENSEMBL_MART_SNP ...
Processing mart ENSEMBL_MART_FUNCGEN ...
Processing mart plants_mart ...
Processing mart plants_variations ...
Processing mart fungi_mart ...
Processing mart fungi_variations ...
Processing mart protists_mart ...
Processing mart protists_variations ...
Processing mart metazoa_mart ...
# A tibble: 15 × 5                                                                                          
   organism_name description               mart      dataset      version
   <chr>         <chr>                     <chr>     <chr>        <chr>  
 1 hsapiens      Human genes (GRCh38.p13)  ENSEMBL_… hsapiens_ge… GRCh38…
 2 hsapiens      Human sequences (GRCh38.… ENSEMBL_… hsapiens_ge… GRCh38…
 3 hsapiens      encode                    ENSEMBL_… hsapiens_en… GRCh38…
 4 hsapiens      marker_feature_end        ENSEMBL_… hsapiens_ma… GRCh38…
 5 hsapiens      marker_feature            ENSEMBL_… hsapiens_ma… GRCh38…
 6 hsapiens      karyotype_end             ENSEMBL_… hsapiens_ka… GRCh38…
 7 hsapiens      karyotype_start           ENSEMBL_… hsapiens_ka… GRCh38…
 8 hsapiens      Human Somatic Short Vari… ENSEMBL_… hsapiens_sn… GRCh38…
 9 hsapiens      Human Structural Variant… ENSEMBL_… hsapiens_st… GRCh38…
10 hsapiens      Human Short Variants (SN… ENSEMBL_… hsapiens_snp GRCh38…
11 hsapiens      Human Somatic Structural… ENSEMBL_… hsapiens_st… GRCh38…
12 hsapiens      Human Regulatory Evidenc… ENSEMBL_… hsapiens_pe… GRCh38…
13 hsapiens      Human Regulatory Feature… ENSEMBL_… hsapiens_re… GRCh38…
14 hsapiens      Human Other Regulatory R… ENSEMBL_… hsapiens_ex… GRCh38…
15 hsapiens      Human miRNA Target Regio… ENSEMBL_… hsapiens_mi… GRCh38…
```

The result is a table storing all `marts` and `datasets` from which annotations can be retrieved
for _Homo sapiens_. Furthermore, a short description as well as the version of the data set
being accessed (very useful for publications) is returned.

Users will observe that 3 different `marts` provide 6 different `datasets` storing annotation information for
_Homo sapiens_.

> **_Please note__*, however, that scientific names of organisms must be written correctly! For ex. "Homo Sapiens" will be treated differently (not recognized) than "Homo sapiens" (recognized).__

Similar to the `biomaRt` package query methodology, users need to specify `attributes` and `filters` to be able to perform
accurate BioMart queries. Here the functions `organismAttributes()` and `organismFilters()` provide useful and intuitive
concepts to obtain this information.


```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# return available attributes for "Homo sapiens"
head(biomartr::organismAttributes("Homo sapiens"), 20)
```

```
1 ensembl_gene_id               Gene stable ID         hsapiens_ge… ENSEMBL_M…
 2 ensembl_gene_id_version       Gene stable ID version hsapiens_ge… ENSEMBL_M…
 3 ensembl_transcript_id         Transcript stable ID   hsapiens_ge… ENSEMBL_M…
 4 ensembl_transcript_id_version Transcript stable ID … hsapiens_ge… ENSEMBL_M…
 5 ensembl_peptide_id            Protein stable ID      hsapiens_ge… ENSEMBL_M…
 6 ensembl_peptide_id_version    Protein stable ID ver… hsapiens_ge… ENSEMBL_M…
 7 ensembl_exon_id               Exon stable ID         hsapiens_ge… ENSEMBL_M…
 8 description                   Gene description       hsapiens_ge… ENSEMBL_M…
 9 chromosome_name               Chromosome/scaffold n… hsapiens_ge… ENSEMBL_M…
10 start_position                Gene start (bp)        hsapiens_ge… ENSEMBL_M…
11 end_position                  Gene end (bp)          hsapiens_ge… ENSEMBL_M…
12 strand                        Strand                 hsapiens_ge… ENSEMBL_M…
13 band                          Karyotype band         hsapiens_ge… ENSEMBL_M…
14 transcript_start              Transcript start (bp)  hsapiens_ge… ENSEMBL_M…
15 transcript_end                Transcript end (bp)    hsapiens_ge… ENSEMBL_M…
16 transcription_start_site      Transcription start s… hsapiens_ge… ENSEMBL_M…
17 transcript_length             Transcript length (in… hsapiens_ge… ENSEMBL_M…
18 transcript_tsl                Transcript support le… hsapiens_ge… ENSEMBL_M…
19 transcript_gencode_basic      GENCODE basic annotat… hsapiens_ge… ENSEMBL_M…
20 transcript_appris             APPRIS annotation      hsapiens_ge… ENSEMBL_M…
```

Users will observe that the `organismAttributes()` function returns a data.frame storing attribute names, data sets, and marts which
are available for `Homo sapiens`. After the ENSEMBL release 87 the `ENSEMBL_MART_SEQUENCE` service provided
by Ensembl does not work properly and thus the `organismAttributes()` function prints out warning messages to make the user aware
when certain marts provided by Ensembl do not work properly, yet. 

An additional feature provided by `organismAttributes()` is the `topic` argument. The `topic` argument allows users to to search for specific attributes,  topics, or categories for faster filtering.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for attribute topic "id"
head(biomartr::organismAttributes("Homo sapiens", topic = "id"), 20)
```

```
  name                          description            dataset      mart      
   <chr>                         <chr>                  <chr>        <chr>     
 1 ensembl_gene_id               Gene stable ID         hsapiens_ge… ENSEMBL_M…
 2 ensembl_gene_id_version       Gene stable ID version hsapiens_ge… ENSEMBL_M…
 3 ensembl_transcript_id         Transcript stable ID   hsapiens_ge… ENSEMBL_M…
 4 ensembl_transcript_id_version Transcript stable ID … hsapiens_ge… ENSEMBL_M…
 5 ensembl_peptide_id            Protein stable ID      hsapiens_ge… ENSEMBL_M…
 6 ensembl_peptide_id_version    Protein stable ID ver… hsapiens_ge… ENSEMBL_M…
 7 ensembl_exon_id               Exon stable ID         hsapiens_ge… ENSEMBL_M…
 8 study_external_id             Study external refere… hsapiens_ge… ENSEMBL_M…
 9 go_id                         GO term accession      hsapiens_ge… ENSEMBL_M…
10 dbass3_id                     DataBase of Aberrant … hsapiens_ge… ENSEMBL_M…
11 dbass5_id                     DataBase of Aberrant … hsapiens_ge… ENSEMBL_M…
12 hgnc_id                       HGNC ID                hsapiens_ge… ENSEMBL_M…
13 protein_id                    INSDC protein ID       hsapiens_ge… ENSEMBL_M…
14 mim_morbid_description        MIM morbid description hsapiens_ge… ENSEMBL_M…
15 mim_morbid_accession          MIM morbid accession   hsapiens_ge… ENSEMBL_M…
16 mirbase_id                    miRBase ID             hsapiens_ge… ENSEMBL_M…
17 refseq_peptide                RefSeq peptide ID      hsapiens_ge… ENSEMBL_M…
18 refseq_peptide_predicted      RefSeq peptide predic… hsapiens_ge… ENSEMBL_M…
19 wikigene_id                   WikiGene ID            hsapiens_ge… ENSEMBL_M…
20 mobidblite                    MobiDBLite             hsapiens_ge… ENSEMBL_M…
```

Now, all `attribute names` having `id` as part of their `name` are being returned.

Another example is `topic = "homolog"`.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for attribute topic "homolog"
head(biomartr::organismAttributes("Homo sapiens", topic = "homolog"), 20)
```

```
  <chr>                                         <chr>           <chr>   <chr> 
 1 mspretus_homolog_ensembl_gene                 Algerian mouse… hsapie… ENSEM…
 2 mspretus_homolog_associated_gene_name         Algerian mouse… hsapie… ENSEM…
 3 mspretus_homolog_ensembl_peptide              Algerian mouse… hsapie… ENSEM…
 4 mspretus_homolog_chromosome                   Algerian mouse… hsapie… ENSEM…
 5 mspretus_homolog_chrom_start                  Algerian mouse… hsapie… ENSEM…
 6 mspretus_homolog_chrom_end                    Algerian mouse… hsapie… ENSEM…
 7 mspretus_homolog_canonical_transcript_protein Query protein … hsapie… ENSEM…
 8 mspretus_homolog_subtype                      Last common an… hsapie… ENSEM…
 9 mspretus_homolog_orthology_type               Algerian mouse… hsapie… ENSEM…
10 mspretus_homolog_perc_id                      %id. target Al… hsapie… ENSEM…
11 mspretus_homolog_perc_id_r1                   %id. query gen… hsapie… ENSEM…
12 mspretus_homolog_goc_score                    Algerian mouse… hsapie… ENSEM…
13 mspretus_homolog_wga_coverage                 Algerian mouse… hsapie… ENSEM…
14 mspretus_homolog_dn                           dN with Algeri… hsapie… ENSEM…
15 mspretus_homolog_ds                           dS with Algeri… hsapie… ENSEM…
16 mspretus_homolog_orthology_confidence         Algerian mouse… hsapie… ENSEM…
17 vpacos_homolog_ensembl_gene                   Alpaca gene st… hsapie… ENSEM…
18 vpacos_homolog_associated_gene_name           Alpaca gene na… hsapie… ENSEM…
19 vpacos_homolog_ensembl_peptide                Alpaca protein… hsapie… ENSEM…
20 vpacos_homolog_chromosome                     Alpaca chromos… hsapie… ENSEM…
```

Or `topic = "dn"` and `topic = "ds"` for `dn` and `ds` value retrieval.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for attribute topic "dn"
head(biomartr::organismAttributes("Homo sapiens", topic = "dn"))
```

```
  name                  description            dataset               mart      
  <chr>                 <chr>                  <chr>                 <chr>     
1 cdna_coding_start     cDNA coding start      hsapiens_gene_ensembl ENSEMBL_M…
2 cdna_coding_end       cDNA coding end        hsapiens_gene_ensembl ENSEMBL_M…
3 mspretus_homolog_dn   dN with Algerian mouse hsapiens_gene_ensembl ENSEMBL_M…
4 vpacos_homolog_dn     dN with Alpaca         hsapiens_gene_ensembl ENSEMBL_M…
5 pformosa_homolog_dn   dN with Amazon molly   hsapiens_gene_ensembl ENSEMBL_M…
6 cpalliatus_homolog_dn dN with Angola colobus hsapiens_gene_ensembl ENSEMBL_M…
```

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for attribute topic "ds"
head(biomartr::organismAttributes("Homo sapiens", topic = "ds"))
```

```
  name                description            dataset               mart        
  <chr>               <chr>                  <chr>                 <chr>       
1 ccds                CCDS ID                hsapiens_gene_ensembl ENSEMBL_MAR…
2 cds_length          CDS Length             hsapiens_gene_ensembl ENSEMBL_MAR…
3 cds_start           CDS start              hsapiens_gene_ensembl ENSEMBL_MAR…
4 cds_end             CDS end                hsapiens_gene_ensembl ENSEMBL_MAR…
5 mspretus_homolog_ds dS with Algerian mouse hsapiens_gene_ensembl ENSEMBL_MAR…
6 vpacos_homolog_ds   dS with Alpaca         hsapiens_gene_ensembl ENSEMBL_MAR…
```

Analogous to the `organismAttributes()` function, the `organismFilters()` function returns
all filters that are available for a query organism of interest.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# return available filters for "Homo sapiens"
head(biomartr::organismFilters("Homo sapiens"), 20)
```

```
   name                                description          dataset    mart    
   <chr>                               <chr>                <chr>      <chr>   
 1 chromosome_name                     Chromosome/scaffold… hsapiens_… ENSEMBL…
 2 start                               Start                hsapiens_… ENSEMBL…
 3 end                                 End                  hsapiens_… ENSEMBL…
 4 band_start                          Band Start           hsapiens_… ENSEMBL…
 5 band_end                            Band End             hsapiens_… ENSEMBL…
 6 marker_start                        Marker Start         hsapiens_… ENSEMBL…
 7 marker_end                          Marker End           hsapiens_… ENSEMBL…
 8 encode_region                       Encode region        hsapiens_… ENSEMBL…
 9 strand                              Strand               hsapiens_… ENSEMBL…
10 chromosomal_region                  e.g. 1:100:10000:-1… hsapiens_… ENSEMBL…
11 with_ccds                           With CCDS ID(s)      hsapiens_… ENSEMBL…
12 with_chembl                         With ChEMBL ID(s)    hsapiens_… ENSEMBL…
13 with_clone_based_ensembl_gene       With Clone-based (E… hsapiens_… ENSEMBL…
14 with_clone_based_ensembl_transcript With Clone-based (E… hsapiens_… ENSEMBL…
15 with_dbass3                         With DataBase of Ab… hsapiens_… ENSEMBL…
16 with_dbass5                         With DataBase of Ab… hsapiens_… ENSEMBL…
17 with_ens_hs_transcript              With Ensembl Human … hsapiens_… ENSEMBL…
18 with_ens_hs_translation             With Ensembl Human … hsapiens_… ENSEMBL…
19 with_entrezgene_trans_name          With EntrezGene tra… hsapiens_… ENSEMBL…
20 with_embl                           With European Nucle… hsapiens_… ENSEMBL…
```

The `organismFilters()` function also allows users to search for filters that correspond to
a specific topic or category.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for filter topic "id"
head(biomartr::organismFilters("Homo sapiens", topic = "id"), 20)
```

```
   name                          description               dataset     mart    
   <chr>                         <chr>                     <chr>       <chr>   
 1 with_protein_id               With INSDC protein ID ID… hsapiens_g… ENSEMBL…
 2 with_mim_morbid               With MIM morbid ID(s)     hsapiens_g… ENSEMBL…
 3 with_refseq_peptide           With RefSeq peptide ID(s) hsapiens_g… ENSEMBL…
 4 with_refseq_peptide_predicted With RefSeq peptide pred… hsapiens_g… ENSEMBL…
 5 ensembl_gene_id               Gene stable ID(s) [e.g. … hsapiens_g… ENSEMBL…
 6 ensembl_gene_id_version       Gene stable ID(s) with v… hsapiens_g… ENSEMBL…
 7 ensembl_transcript_id         Transcript stable ID(s) … hsapiens_g… ENSEMBL…
 8 ensembl_transcript_id_version Transcript stable ID(s) … hsapiens_g… ENSEMBL…
 9 ensembl_peptide_id            Protein stable ID(s) [e.… hsapiens_g… ENSEMBL…
10 ensembl_peptide_id_version    Protein stable ID(s) wit… hsapiens_g… ENSEMBL…
11 ensembl_exon_id               Exon ID(s) [e.g. ENSE000… hsapiens_g… ENSEMBL…
12 dbass3_id                     DataBase of Aberrant 3' … hsapiens_g… ENSEMBL…
13 dbass5_id                     DataBase of Aberrant 5' … hsapiens_g… ENSEMBL…
14 hgnc_id                       HGNC ID(s) [e.g. HGNC:10… hsapiens_g… ENSEMBL…
15 protein_id                    INSDC protein ID(s) [e.g… hsapiens_g… ENSEMBL…
16 mim_morbid_accession          MIM morbid accession(s) … hsapiens_g… ENSEMBL…
17 mirbase_id                    miRBase ID(s) [e.g. hsa-… hsapiens_g… ENSEMBL…
18 refseq_peptide                RefSeq peptide ID(s) [e.… hsapiens_g… ENSEMBL…
19 refseq_peptide_predicted      RefSeq peptide predicted… hsapiens_g… ENSEMBL…
20 wikigene_id                   WikiGene ID(s) [e.g. 1]   hsapiens_g… ENSEMBL…
```

## Construct BioMart queries with `biomartr`

The short introduction to the functionality of
`organismBM()`, `organismAttributes()`, and `organismFilters()`
will allow users to perform BioMart queries in a very intuitive 
organism centric way. The main function to perform BioMart queries
is `biomart()`.


For the following examples we will assume that we are interested in the annotation of specific genes from the _Homo sapiens_ proteome. We want to map the corresponding refseq gene id to a set of other gene ids used in other databases. For this purpose, first we need consult the `organismAttributes()` function.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)

head(biomartr::organismAttributes("Homo sapiens", topic = "id"))
```
```
 name                          description                  dataset    mart   
  <chr>                         <chr>                        <chr>      <chr>  
1 ensembl_gene_id               Gene stable ID               hsapiens_… ENSEMB…
2 ensembl_gene_id_version       Gene stable ID version       hsapiens_… ENSEMB…
3 ensembl_transcript_id         Transcript stable ID         hsapiens_… ENSEMB…
4 ensembl_transcript_id_version Transcript stable ID version hsapiens_… ENSEMB…
5 ensembl_peptide_id            Protein stable ID            hsapiens_… ENSEMB…
6 ensembl_peptide_id_version    Protein stable ID version    hsapiens_… ENSEMB…
```

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# retrieve the proteome of Homo sapiens from refseq
file_path <- biomartr::getProteome( db       = "refseq",
                                    organism = "Homo sapiens",
                                    path     = file.path("_ncbi_downloads","proteomes") )

Hsapiens_proteome <- biomartr::read_proteome(file_path, format = "fasta")

# remove splice variants from id
gene_set <- unlist(sapply(strsplit(Hsapiens_proteome@ranges@NAMES[1:5], ".",fixed = TRUE), function(x) x[1]))

result_BM <- biomartr::biomart( genes      = gene_set, # genes were retrieved using biomartr::getGenome()
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = c("ensembl_gene_id","ensembl_peptide_id"), # attributes were selected with biomartr::getAttributes()
                                filters    = "refseq_peptide") # specify what ID type was stored in the fasta file retrieved with biomartr::getGenome()

result_BM 
```

```
  refseq_peptide ensembl_gene_id ensembl_peptide_id
1      NP_000005 ENSG00000175899    ENSP00000323929
2      NP_000006 ENSG00000156006    ENSP00000286479
3      NP_000007 ENSG00000117054    ENSP00000359878
4      NP_000008 ENSG00000122971    ENSP00000242592
5      NP_000009 ENSG00000072778    ENSP00000349297
```


The `biomart()` function takes as arguments a set of genes (gene ids specified in the `filter` argument), the corresponding `mart` and `dataset`, as well as the `attributes` which shall be returned.

## Gene Ontology 

The `biomartr` package also enables a fast and intuitive retrieval of GO terms
and additional information via the `getGO()` function. Several databases can be selected
to retrieve GO annotation information for a set of query genes. So far, the `getGO()`
function allows GO information retrieval from the [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart) database. 

In this example we will retrieve GO information for a set of _Homo sapiens_ genes
stored as `hgnc_symbol`.

### GO Annotation Retrieval via BioMart

The `getGO()` function takes several arguments as input to retrieve GO information from BioMart. 
First, the scientific name of the `organism` of interest needs to be specified. Furthermore, a set of
`gene ids` as well as their corresponding `filter` notation (`GUCA2A` gene ids have `filter` notation `hgnc_symbol`; see `organismFilters()` for details)
need to be specified. The `database` argument then defines the database from which GO information shall be retrieved.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# search for GO terms of an example Homo sapiens gene
GO_tbl <- biomartr::getGO(organism = "Homo sapiens", 
                          genes    = "GUCA2A",
                          filters  = "hgnc_symbol")

GO_tbl
```


Hence, for each _gene id_ the resulting table stores all annotated GO terms found in [Ensembl Biomart](http://ensemblgenomes.org/info/access/biomart).

---
title: "Ensembl BioMart Examples"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Ensembl BioMart Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```

## Use Case #1: Functional Annotation of Genes Sharing a Common Evolutionary History

Evolutionary Transcriptomics aims to predict stages or periods of evolutionary conservation in
biological processes on the transcriptome level. However, finding genes sharing a [common evolutionary history](https://github.com/HajkD/myTAI/blob/master/vignettes/Enrichment.Rmd) could reveal how the the biological process might have evolved in the first place.

In this `Use Case` we will combine functional and biological annotation obtained with `biomartr` with enriched genes obtained with [PlotEnrichment()](https://github.com/HajkD/myTAI/blob/master/vignettes/Enrichment.Rmd). 

### Step 1

For the following example we will use the dataset an enrichment analyses found in [PlotEnrichment()](https://github.com/HajkD/myTAI/blob/master/vignettes/Enrichment.Rmd).

Install and load the [myTAI](https://github.com/HajkD/myTAI) package:

```{r, eval=FALSE}
# install myTAI
install.packages("myTAI")

# load myTAI
library(myTAI)
```

Download the `Phylostratigraphic Map` of _D. rerio_:

```r
# download the Phylostratigraphic Map of Danio rerio
# from Sestak and Domazet-Loso, 2015
```

The dataset comes from `supplementary file 3` of this publication: https://academic.oup.com/mbe/article/32/2/299/1058654#77837069

After downloading `supplementary file 3`, you will find the file “TableS3-2.xlsx” which can be used in the following biomartr functions.



Read the `*.xlsx` file storing the `Phylostratigraphic Map` of _D. rerio_ and format it for the use with `myTAI`:

```r
# install the readxl package
install.packages("readxl")

# load package readxl
library(readxl)

# read the excel file
DrerioPhyloMap.MBEa <- read_excel("TableS3-2.xlsx", sheet = 1, skip = 4)

# format Phylostratigraphic Map for use with myTAI
Drerio.PhyloMap <- DrerioPhyloMap.MBEa[ , 1:2]

# have a look at the final format
head(Drerio.PhyloMap)
```

```
  Phylostrata            ZFIN_ID
1           1 ZDB-GENE-000208-13
2           1 ZDB-GENE-000208-17
3           1 ZDB-GENE-000208-18
4           1 ZDB-GENE-000208-23
5           1  ZDB-GENE-000209-3
6           1  ZDB-GENE-000209-4
```

Now, `Drerio.PhyloMap` stores the `Phylostratigraphic Map` of _D. rerio_ which is used
as background set to perform enrichment analyses with `PlotEnrichment()` from `myTAI`.

### Enrichment Analyses

Now, the `PlotEnrichment()` function visualizes the over- and underrepresented `Phylostrata` of brain specific genes when compared with the total number of genes stored in the `Phylostratigraphic Map` of _D. rerio_.


```{r,eval=FALSE}
library(readxl)

# read expression data (organ specific genes) from Sestak and Domazet-Loso, 2015
Drerio.OrganSpecificExpression <- read_excel("TableS3-2.xlsx", sheet = 2, skip = 3)

# select only brain specific genes
Drerio.Brain.Genes <- unlist(unique(na.omit(Drerio.OrganSpecificExpression[ , "brain"])))

# visualize enriched Phylostrata of genes annotated as brain specific
PlotEnrichment(Drerio.PhyloMap,
               test.set     = Drerio.Brain.Genes,
               measure      = "foldchange",
               use.only.map = TRUE,
               legendName   = "PS")
```

Users will observe that for example brain genes deriving from PS5 are significantly enriched. 


Now we can select all brain genes originating in PS5 using the `SelectGeneSet()` function from `myTAI`. Please notice that `SelectGeneSet()` can only be used with phylostratigraphic maps only (`use.map.only = TRUE` argument) since myTAI version > 0.3.0.

```{r,eval=FALSE}
BrainGenes <- SelectGeneSet(ExpressionSet = Drerio.PhyloMap,
                            gene.set      = Drerio.Brain.Genes,
                            use.only.map  = TRUE)

# select only brain genes originating in PS5
BrainGenes.PS5 <- BrainGenes[which(BrainGenes[ , "Phylostrata"] == 5), ]

# look at the results
head(BrainGenes.PS5)
```

```
      Phylostrata           ZFIN_ID
14851           5 ZDB-GENE-000210-6
14852           5 ZDB-GENE-000210-7
14853           5 ZDB-GENE-000328-4
14856           5 ZDB-GENE-000411-1
14857           5 ZDB-GENE-000427-4
14860           5 ZDB-GENE-000526-1
```

Now users can perform the `biomart()` function to obtain the functional annotation of brain genes originating in PS5.

For this purpose, first we need to find the filter name of the corresponding gene ids such as `ZDB-GENE-000210-6`.

```{r, eval=FALSE}
# find filter for zfin.org ids
organismFilters("Danio rerio", topic = "zfin_id")
```

```
                            name                           description             dataset
52                  with_zfin_id                       with ZFIN ID(s) drerio_gene_ensembl
53  with_zfin_id_transcript_name          with ZFIN transcript name(s) drerio_gene_ensembl
103                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-060825-136] drerio_gene_ensembl
274                 with_zfin_id                       with ZFIN ID(s)    drerio_gene_vega
286                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-121214-212]    drerio_gene_vega
366                 with_zfin_id                       with ZFIN ID(s) drerio_gene_ensembl
367 with_zfin_id_transcript_name          with ZFIN transcript name(s) drerio_gene_ensembl
417                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-060825-136] drerio_gene_ensembl
588                 with_zfin_id                       with ZFIN ID(s)    drerio_gene_vega
600                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-121214-212]    drerio_gene_vega
680                 with_zfin_id                       with ZFIN ID(s) drerio_gene_ensembl
681 with_zfin_id_transcript_name          with ZFIN transcript name(s) drerio_gene_ensembl
731                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-060825-136] drerio_gene_ensembl
902                 with_zfin_id                       with ZFIN ID(s)    drerio_gene_vega
914                      zfin_id ZFIN ID(s) [e.g. ZDB-GENE-121214-212]    drerio_gene_vega
                    mart
52  ENSEMBL_MART_ENSEMBL
53  ENSEMBL_MART_ENSEMBL
103 ENSEMBL_MART_ENSEMBL
274 ENSEMBL_MART_ENSEMBL
286 ENSEMBL_MART_ENSEMBL
366 ENSEMBL_MART_ENSEMBL
367 ENSEMBL_MART_ENSEMBL
417 ENSEMBL_MART_ENSEMBL
588 ENSEMBL_MART_ENSEMBL
600 ENSEMBL_MART_ENSEMBL
680 ENSEMBL_MART_ENSEMBL
681 ENSEMBL_MART_ENSEMBL
731 ENSEMBL_MART_ENSEMBL
902 ENSEMBL_MART_ENSEMBL
914 ENSEMBL_MART_ENSEMBL
```

Now users can retrieve the corresponding GO attribute of _D. rerio_ with `organismAttributes`.

```{r,eval=FALSE}
# find go attribute term for D. rerio
organismAttributes("Danio rerio", topic = "go")
```

```
                                              name                             description
33                                           go_id                       GO Term Accession
36                                 go_linkage_type                   GO Term Evidence Code
38                            goslim_goa_accession                 GOSlim GOA Accession(s)
39                          goslim_goa_description                  GOSlim GOA Description
516                  ggorilla_homolog_ensembl_gene                 Gorilla Ensembl Gene ID
517  ggorilla_homolog_canomical_transcript_protein      Canonical Protein or Transcript ID
518               ggorilla_homolog_ensembl_peptide              Gorilla Ensembl Protein ID
519                    ggorilla_homolog_chromosome                 Gorilla Chromosome Name
520                   ggorilla_homolog_chrom_start           Gorilla Chromosome Start (bp)
521                     ggorilla_homolog_chrom_end             Gorilla Chromosome End (bp)
522                ggorilla_homolog_orthology_type                           Homology Type
523                       ggorilla_homolog_subtype                                Ancestor
524          ggorilla_homolog_orthology_confidence    Orthology confidence [0 low, 1 high]
525                       ggorilla_homolog_perc_id   % Identity with respect to query gene
526                    ggorilla_homolog_perc_id_r1 % Identity with respect to Gorilla gene
527                            ggorilla_homolog_dn                                      dN
528                            ggorilla_homolog_ds                                      dS
1240                                         go_id                                   GO ID
1241                                      quick_go                             Quick GO ID
1370                                         go_id                       GO Term Accession
1373                               go_linkage_type                   GO Term Evidence Code
1375                          goslim_goa_accession                 GOSlim GOA Accession(s)
1376                        goslim_goa_description                  GOSlim GOA Description
1853                 ggorilla_homolog_ensembl_gene                 Gorilla Ensembl Gene ID
1854 ggorilla_homolog_canomical_transcript_protein      Canonical Protein or Transcript ID
1855              ggorilla_homolog_ensembl_peptide              Gorilla Ensembl Protein ID
1856                   ggorilla_homolog_chromosome                 Gorilla Chromosome Name
1857                  ggorilla_homolog_chrom_start           Gorilla Chromosome Start (bp)
1858                    ggorilla_homolog_chrom_end             Gorilla Chromosome End (bp)
1859               ggorilla_homolog_orthology_type                           Homology Type
1860                      ggorilla_homolog_subtype                                Ancestor
1861         ggorilla_homolog_orthology_confidence    Orthology confidence [0 low, 1 high]
1862                      ggorilla_homolog_perc_id   % Identity with respect to query gene
1863                   ggorilla_homolog_perc_id_r1 % Identity with respect to Gorilla gene
1864                           ggorilla_homolog_dn                                      dN
1865                           ggorilla_homolog_ds                                      dS
2577                                         go_id                                   GO ID
2578                                      quick_go                             Quick GO ID
2707                                         go_id                       GO Term Accession
2710                               go_linkage_type                   GO Term Evidence Code
2712                          goslim_goa_accession                 GOSlim GOA Accession(s)
2713                        goslim_goa_description                  GOSlim GOA Description
3190                 ggorilla_homolog_ensembl_gene                 Gorilla Ensembl Gene ID
3191 ggorilla_homolog_canomical_transcript_protein      Canonical Protein or Transcript ID
3192              ggorilla_homolog_ensembl_peptide              Gorilla Ensembl Protein ID
3193                   ggorilla_homolog_chromosome                 Gorilla Chromosome Name
3194                  ggorilla_homolog_chrom_start           Gorilla Chromosome Start (bp)
3195                    ggorilla_homolog_chrom_end             Gorilla Chromosome End (bp)
3196               ggorilla_homolog_orthology_type                           Homology Type
3197                      ggorilla_homolog_subtype                                Ancestor
3198         ggorilla_homolog_orthology_confidence    Orthology confidence [0 low, 1 high]
3199                      ggorilla_homolog_perc_id   % Identity with respect to query gene
3200                   ggorilla_homolog_perc_id_r1 % Identity with respect to Gorilla gene
3201                           ggorilla_homolog_dn                                      dN
3202                           ggorilla_homolog_ds                                      dS
3914                                         go_id                                   GO ID
3915                                      quick_go                             Quick GO ID
                 dataset                 mart
33   drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
36   drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
38   drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
39   drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
516  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
517  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
518  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
519  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
520  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
521  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
522  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
523  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
524  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
525  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
526  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
527  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
528  drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1240    drerio_gene_vega ENSEMBL_MART_ENSEMBL
1241    drerio_gene_vega ENSEMBL_MART_ENSEMBL
1370 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1373 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1375 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1376 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1853 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1854 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1855 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1856 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1857 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1858 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1859 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1860 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1861 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1862 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1863 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1864 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
1865 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
2577    drerio_gene_vega ENSEMBL_MART_ENSEMBL
2578    drerio_gene_vega ENSEMBL_MART_ENSEMBL
2707 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
2710 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
2712 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
2713 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3190 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3191 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3192 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3193 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3194 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3195 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3196 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3197 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3198 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3199 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3200 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3201 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3202 drerio_gene_ensembl ENSEMBL_MART_ENSEMBL
3914    drerio_gene_vega ENSEMBL_MART_ENSEMBL
3915    drerio_gene_vega ENSEMBL_MART_ENSEMBL
```

Now users can specify the filter `zfin_id` and attribute `go_id` to retrieve the GO terms of corresponding gene ids (__Please note that this will take some time__).

```{r, eval=FALSE}
# retrieve GO terms of D. rerio brain genes originating in PS5
GO_tbl.BrainGenes <- biomart(genes      = unlist(BrainGenes.PS5[ , "ZFIN_ID"]),
                             mart       = "ENSEMBL_MART_ENSEMBL",
                             dataset    = "drerio_gene_ensembl",
                             attributes = "go_id",
                             filters    = "zfin_id")

head(GO_tbl.BrainGenes)
```

```
            zfin_id      go_id
1 ZDB-GENE-000210-6 GO:0060037
2 ZDB-GENE-000210-6 GO:0046983
3 ZDB-GENE-000210-7 GO:0046983
4 ZDB-GENE-000328-4 GO:0007275
5 ZDB-GENE-000328-4 GO:0007166
6 ZDB-GENE-000328-4 GO:0035567
```





---
title: "NCBI Database Retrieval"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NCBI Database Retrieval}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```

# Retrieve Sequence Databases from NCBI

> **_NOTE:_** To make sure that you have a sufficiently stable (internet) connection between R and the respective databases, please set the default `timeout` setting __on your local machine__ from 60sec to at least 300000sec before running any retrieval functions via:

```r
options(timeout = 300000)
```

## Getting Started
NCBI stores a variety of specialized database such as [Genbank, RefSeq, Taxonomy, SNP, etc.](https://www.ncbi.nlm.nih.gov/guide/data-software/) on their servers.
The `download.database()` and `download.database.all()` functions implemented in `biomartr` allows users to download these databases from NCBI.
This process might be very useful for downstream analyses such as sequence searches with e.g. BLAST. For this purpose see the R package [metablastr](https://github.com/HajkD/metablastr) which aims to seamlessly inegrate `biomartr` based genomic data retrieval with downsteam large-scale BLAST searches.

* [1. List available NCBI databases with `listNCBIDatabases()`](#ist-available-databases)
* [2. Download NCBI databases with `download.database.all()`](#download-ncbi-databases)
    - [2.1 Example NCBI nr retrieval](#example-ncbi-nr)
    - [2.2 Example NCBI nt retrieval](#example-ncbi-nt)
    - [2.3 Example NCBI RefSeq retrieval](#example-ncbi-refseq)
    - [2.4 Example Protein Database (PDB) retrieval](#example-pdb)
    - [2.5 Example NCBI Taxonomy retrieval](#example-ncbi-taxonomy)
    - [2.6 Example NCBI Swissprot retrieval](#example-ncbi-swissprot)
    - [2.7 Example NCBI CDD Delta retrieval](#example-ncbi-cdd-delta)


## List available databases

Before downloading specific databases from NCBI users might want to list available databases. Using the `listNCBIDatabases()` function users can retrieve a list of  available databases stored on NCBI.

```{r, eval=FALSE}
# retrieve a list of available sequence databases at NCBI
biomartr::listNCBIDatabases(db = "all")
```

```
[1] "16S_ribosomal_RNA"                     "16S_ribosomal_RNA-nucl-metadata"      
 [3] "18S_fungal_sequences"                  "18S_fungal_sequences-nucl-metadata"   
 [5] "28S_fungal_sequences"                  "28S_fungal_sequences-nucl-metadata"   
 [7] "Betacoronavirus"                       "Betacoronavirus-nucl-metadata"        
 [9] "blastdb-manifest"                      "blastdb-metadata-1-1"                 
[11] "blastdbv5"                             "cdd_delta"                            
[13] "cloud"                                 "env_nr"                               
[15] "env_nr-prot-metadata"                  "env_nt"                               
[17] "env_nt-nucl-metadata"                  "FASTA"                                
[19] "human_genome"                          "human_genome-nucl-metadata"           
[21] "ITS_eukaryote_sequences"               "ITS_eukaryote_sequences-nucl-metadata"
[23] "ITS_RefSeq_Fungi"                      "ITS_RefSeq_Fungi-nucl-metadata"       
[25] "landmark"                              "landmark-prot-metadata"               
[27] "LSU_eukaryote_rRNA"                    "LSU_eukaryote_rRNA-nucl-metadata"     
[29] "LSU_prokaryote_rRNA"                   "LSU_prokaryote_rRNA-nucl-metadata"    
[31] "mito"                                  "mito-nucl-metadata"                   
[33] "mouse_genome"                          "mouse_genome-nucl-metadata"           
[35] "nr"                                    "nr-prot-metadata"                     
[37] "nt"                                    "nt-nucl-metadata"                     
[39] "pataa"                                 "pataa-prot-metadata"                  
[41] "patnt"                                 "patnt-nucl-metadata"                  
[43] "pdbaa"                                 "pdbaa-prot-metadata"                  
[45] "pdbnt"                                 "pdbnt-nucl-metadata"                  
[47] "ref_euk_rep_genomes"                   "ref_euk_rep_genomes-nucl-metadata"    
[49] "ref_prok_rep_genomes"                  "ref_prok_rep_genomes-nucl-metadata"   
[51] "ref_viroids_rep_genomes"               "ref_viroids_rep_genomes-nucl-metadata"
[53] "ref_viruses_rep_genomes"               "ref_viruses_rep_genomes-nucl-metadata"
[55] "refseq_protein"                        "refseq_protein-prot-metadata"         
[57] "refseq_rna"                            "refseq_rna-nucl-metadata"             
[59] "refseq_select_prot"                    "refseq_select_prot-prot-metadata"     
[61] "refseq_select_rna"                     "refseq_select_rna-nucl-metadata"      
[63] "SSU_eukaryote_rRNA"                    "SSU_eukaryote_rRNA-nucl-metadata"     
[65] "swissprot"                             "swissprot-prot-metadata"              
[67] "taxdb"                                 "taxdb-metadata"                       
[69] "tsa_nr"                                "tsa_nr-prot-metadata"                 
[71] "tsa_nt"                                "tsa_nt-nucl-metadata"                 
[73] "v4"                                    "v5" 
```


However, in case users already know which database they would like to retrieve
they can filter for the exact files by specifying the NCBI database name. In the following example all sequence files that are part of the `NCBI nr` database shall be retrieved.


First, the `listNCBIDatabases(db = "nr")` allows to list all files corresponding to the `nr` database. 

```{r, eval=FALSE}
# show all NCBI nr files
biomartr::listNCBIDatabases(db = "nr")
```

```
[1] "nr.00.tar.gz"          "nr.01.tar.gz"         
 [3] "nr.02.tar.gz"          "nr.03.tar.gz"         
 [5] "nr.04.tar.gz"          "nr.05.tar.gz"         
 [7] "nr.06.tar.gz"          "nr.07.tar.gz"         
 [9] "nr.08.tar.gz"          "nr.09.tar.gz"         
[11] "nr.10.tar.gz"          "nr.11.tar.gz"         
[13] "nr.12.tar.gz"          "nr.13.tar.gz"         
[15] "nr.14.tar.gz"          "nr.15.tar.gz"         
[17] "nr.16.tar.gz"          "nr.17.tar.gz"         
[19] "nr.18.tar.gz"          "nr.19.tar.gz"         
[21] "nr.20.tar.gz"          "nr.21.tar.gz"         
[23] "nr.22.tar.gz"          "nr.23.tar.gz"         
[25] "nr.24.tar.gz"          "nr.25.tar.gz"         
[27] "nr.26.tar.gz"          "nr.27.tar.gz"         
[29] "nr.28.tar.gz"          "nr.29.tar.gz"         
[31] "nr.30.tar.gz"          "nr.31.tar.gz"         
[33] "nr.32.tar.gz"          "nr.33.tar.gz"         
[35] "nr.34.tar.gz"          "nr.35.tar.gz"         
[37] "nr.36.tar.gz"          "nr.37.tar.gz"         
[39] "nr.38.tar.gz"          "nr.39.tar.gz"         
[41] "nr.40.tar.gz"          "nr.41.tar.gz"         
[43] "nr.42.tar.gz"          "nr.43.tar.gz"         
[45] "nr.44.tar.gz"          "nr.45.tar.gz"         
[47] "nr.46.tar.gz"          "nr-prot-metadata.json"
[49] "nr.47.tar.gz"          "nr.48.tar.gz"         
[51] "nr.49.tar.gz"          "nr.50.tar.gz"         
[53] "nr.51.tar.gz"          "nr.52.tar.gz"         
[55] "nr.53.tar.gz"          "nr.54.tar.gz"         
[57] "nr.55.tar.gz"
```

The output illustrates that the `NCBI nr` database has been separated into several sub-data-packages.

Further examples are:


```{r, eval=FALSE}
# show all NCBI nt files
biomartr::listNCBIDatabases(db = "nt")
```

```
[1] "nt.00.tar.gz"          "nt.01.tar.gz"         
 [3] "nt.02.tar.gz"          "nt.03.tar.gz"         
 [5] "nt.04.tar.gz"          "nt.05.tar.gz"         
 [7] "nt.06.tar.gz"          "nt.07.tar.gz"         
 [9] "nt.08.tar.gz"          "nt.09.tar.gz"         
[11] "nt.10.tar.gz"          "nt.11.tar.gz"         
[13] "nt.12.tar.gz"          "nt.13.tar.gz"         
[15] "nt.14.tar.gz"          "nt.15.tar.gz"         
[17] "nt.16.tar.gz"          "nt.17.tar.gz"         
[19] "nt.18.tar.gz"          "nt.19.tar.gz"         
[21] "nt.20.tar.gz"          "nt.21.tar.gz"         
[23] "nt.22.tar.gz"          "nt.23.tar.gz"         
[25] "nt.24.tar.gz"          "nt.25.tar.gz"         
[27] "nt.26.tar.gz"          "nt.27.tar.gz"         
[29] "nt.28.tar.gz"          "nt.29.tar.gz"         
[31] "nt.30.tar.gz"          "nt.31.tar.gz"         
[33] "nt.32.tar.gz"          "nt.33.tar.gz"         
[35] "nt.34.tar.gz"          "nt.35.tar.gz"         
[37] "nt.36.tar.gz"          "nt.37.tar.gz"         
[39] "nt-nucl-metadata.json" "nt.38.tar.gz"         
[41] "nt.39.tar.gz"          "nt.40.tar.gz"         
[43] "nt.41.tar.gz"          "nt.42.tar.gz"         
[45] "nt.43.tar.gz"          "nt.44.tar.gz"         
[47] "nt.45.tar.gz"          "nt.46.tar.gz"         
[49] "nt.47.tar.gz"          "nt.48.tar.gz"
```


```{r, eval=FALSE}
# show all NCBI RefSeq (only proteomes)
head(biomartr::listNCBIDatabases(db = "refseq_protein"), 20)
```

```
[1] "refseq_protein.00.tar.gz" "refseq_protein.01.tar.gz"
 [3] "refseq_protein.02.tar.gz" "refseq_protein.03.tar.gz"
 [5] "refseq_protein.04.tar.gz" "refseq_protein.05.tar.gz"
 [7] "refseq_protein.06.tar.gz" "refseq_protein.07.tar.gz"
 [9] "refseq_protein.08.tar.gz" "refseq_protein.09.tar.gz"
[11] "refseq_protein.10.tar.gz" "refseq_protein.11.tar.gz"
[13] "refseq_protein.12.tar.gz" "refseq_protein.13.tar.gz"
[15] "refseq_protein.14.tar.gz" "refseq_protein.15.tar.gz"
[17] "refseq_protein.16.tar.gz" "refseq_protein.17.tar.gz"
[19] "refseq_protein.18.tar.gz" "refseq_protein.19.tar.gz"
```


```{r, eval=FALSE}
# show all NCBI RefSeq (only RNA)
biomartr::listNCBIDatabases(db = "refseq_rna")
```

```
[1] "refseq_rna.00.tar.gz"          "refseq_rna.01.tar.gz"         
 [3] "refseq_rna.02.tar.gz"          "refseq_rna.03.tar.gz"         
 [5] "refseq_rna.04.tar.gz"          "refseq_rna.05.tar.gz"         
 [7] "refseq_rna.06.tar.gz"          "refseq_rna.07.tar.gz"         
 [9] "refseq_rna.08.tar.gz"          "refseq_rna-nucl-metadata.json"
[11] "refseq_rna.09.tar.gz" 
```

```{r, eval=FALSE}
# show NCBI swissprot
biomartr::listNCBIDatabases(db = "swissprot")
```

```
[1] "swissprot.tar.gz"             "swissprot-prot-metadata.json"
```

```{r, eval=FALSE}
# show NCBI PDB
biomartr::listNCBIDatabases(db = "pdb")
```

```
 [1] "pdbaa.tar.gz"             "pdbnt.tar.gz"            
[3] "pdbaa-prot-metadata.json" "pdbnt-nucl-metadata.json"
```

```{r, eval=FALSE}
# show NCBI Human database
biomartr::listNCBIDatabases(db = "human")
```

```
1] "human_genome.00.tar.gz"          "human_genome.01.tar.gz"         
[3] "human_genome-nucl-metadata.json"
```


__Please not that all lookup and retrieval function will only work properly when a sufficient internet connection is provided.__

In a next step users can use the `listNCBIDatabases()` and `download.database.all()` functions to retrieve 
all files corresponding to a specific NCBI database.

## Download NCBI databases

Using the same search strategy by specifying the database name as described above, users can now download these databases using the `download.database.all()` function.

For downloading only single files users can type:

### Example NCBI nr

```{r, eval=FALSE}
# download the entire NCBI nr database
biomartr::download.database.all(db = "nr", path = "nr")
```
This command will download the pre-formatted (by makeblastdb formatted) database version is retrieved.

Using this command, all `NCBI nr` files are loaded into the `nr` folder (`path = "nr"`). For each data package, `biomartr` checks the `md5checksum` of the downloaded file and the file stored online to make sure that internet connection losses didn't currupt the file. In case you see a warning message notifying you about not-matching `md5checksum` values, please re-download the corresponding data package by re-running the `download.database.all()` command. From my own experience this can happen when server connections or internet connections are not very stable during the download process of large data chunks.


### Example NCBI nt

The same approach can be applied to all other databases mentioned above, e.g.:

```{r, eval=FALSE}
# download the entire NCBI nt database
biomartr::download.database.all(db = "nt", path = "nt")
```

### Example NCBI RefSeq

```{r, eval=FALSE}
# download the entire NCBI refseq (protein) database
biomartr::download.database.all(db = "refseq_protein", path = "refseq_protein")
```


### Example PDB

```{r, eval=FALSE}
# download the entire NCBI PDB database
biomartr::download.database.all(db = "pdb", path = "pdb")
```


### Example NCBI Taxonomy

Download NCBI Taxonomy via:

```{r, eval=FALSE}
# download the entire NCBI taxonomy database
biomartr::download.database.all(db = "taxdb", path = "taxdb")
```

```
Starting download of the files: taxdb.tar.gz, taxdb.btd, taxdb.bti ...
This download process may take a while due to the large size of the individual data chunks ...
Starting download process of file: taxdb.tar.gz ...
Checking md5 hash of file: taxdb.tar.gz ...
The md5 hash of file 'taxdb.tar.gz' matches!
File 'taxdb/taxdb.tar.gz has successfully been retrieved.
Starting download process of file: taxdb.btd ...
Checking md5 hash of file: taxdb.btd ...
The md5 hash of file 'taxdb.btd' matches!
File 'taxdb/taxdb.btd has successfully been retrieved.
Starting download process of file: taxdb.bti ...
Checking md5 hash of file: taxdb.bti ...
The md5 hash of file 'taxdb.bti' matches!
File 'taxdb/taxdb.bti has successfully been retrieved.
Download process is finished and files are stored in 'taxdb'.
```

### Example NCBI Swissprot

Download NCBI Swissprot via:

```{r, eval=FALSE}
# download the entire NCBI swissprot database
biomartr::download.database.all(db = "swissprot", path = "swissprot")
```

### Example NCBI CDD Delta

Download NCBI CDD Delta via:

```{r, eval=FALSE}
# download the entire NCBI CDD Delta database
biomartr::download.database.all(db = "cdd_delta", path = "cdd_delta")
```



For each data package, `biomartr` checks the `md5checksum` of the downloaded file and the file stored online to make sure that internet connection losses didn't currupt the file. In case you see a warning message notifying you about not-matching `md5checksum` values, please re-download the corresponding data package. From my own experience this can happen when server connections or internet connections are not very stable during the download process of large data chunks.

__Please notice that most of these databases are very large, so users should take of of providing a stable internet connection throughout the download process.__
---
title: "Meta-Genome Retrieval"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Meta-Genome Retrieval}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```

> **_NOTE:_** To make sure that you have a sufficiently stable (internet) connection between R and the respective databases, please set the default `timeout` setting __on your local machine__ from 60sec to at least 30000sec before running any retrieval functions via:

```r
options(timeout = 30000)
```


## Topics
- [1. Perform Meta-Genome Retrieval for Specific Kingdoms of Life](#perform-meta-genome-retrieval)
    - [1.1 Retrieve Genomic Sequences](#retrieve-genomic-sequences)
        - [1.1.1 Retrieval from `NCBI RefSeq`](#retrieval-from-ncbi-refseq)
        - [1.1.2 Retrieval from `NCBI Genbank`](#retrieval-from-ncbi-genbank)
        - [1.1.3 Retrieval from `ENSEMBL`](#retrieval-from-ensembl)
    - [1.2 Retrieve Protein Sequences](#retrieve-protein-sequences)
        - [1.2.1 Retrieval from `NCBI RefSeq`](#retrieval-from-ncbi-refseq-1)
        - [1.2.2 Retrieval from `NCBI Genbank`](#retrieval-from-ncbi-genbank-1)
        - [1.2.3 Retrieval from `ENSEMBL`](#retrieval-from-ensembl-1)
    - [1.3 Retrieve CDS Sequences](#retrieve-cds-sequences)
        - [1.3.1 Retrieval from `NCBI RefSeq`](#retrieval-from-ncbi-refseq-2)
        - [1.3.2 Retrieval from `NCBI Genbank`](#retrieval-from-ncbi-genbank-2)
        - [1.3.3 Retrieval from `ENSEMBL`](#retrieval-from-ensembl-2)
    - [1.4 Retrieve GFF files](#retrieve-gff-files)
    - [1.5 Retrieve GTF files](#retrieve-gtf-files)
    - [1.6 Retrieve RNA sequences](#retrieve-rna-sequences)
        - [1.6.1 Retrieval from `NCBI RefSeq`](#retrieval-from-ncbi-refseq-3)
        - [1.6.2 Retrieval from `NCBI Genbank`](#retrieval-from-ncbi-genbank-3)
        - [1.6.3 Retrieval from `ENSEMBL`](#retrieval-from-ensembl-3)
    - [1.7 Retrieve Repeat Masker Repeat Annotation File](#retrieve-repeat-masker-sequences)
        - [1.7.1 Retrieval from `NCBI RefSeq`](#retrieval-from-ncbi-refseq-4)
        - [1.7.2 Retrieval from `NCBI Genbank`](#retrieval-from-ncbi-genbank-4)
- [2. Retrieve groups or subgroups of species ](#retrieve-groups-or-subgroups-of-species)
    - [2.1 Example retrieval of all `Gammaproteobacteria` genomes from `NCBI RefSeq`](#example-retrieval-of-all-gammaproteobacteria-genomes-from-ncbi-refseq)
    - [2.2 Example retrieval of all `Adenoviridae` genomes from `NCBI RefSeq`](#example-retrieval-of-all-adenoviridae-genomes-from-ncbi-refseq)
- [3. Meta retrieval of genome assembly quality information](#meta-retrieval-of-genome-assembly-quality-information)
- [4. Retrieve data from metagenome projects such as `human gut metagenome` project from `NCBI Genbank`](#metagenome-project-retrieval-from-ncbi-genbank)
- [5. Retrieve Individual Genomes for all Species in the Tree of Life](#retrieve-individual-genomes-for-all-species-in-the-tree-of-life)
    - [5.1 Example: Genome Retrieval](#genome-retrieval)
    - [5.2 Example: Proteome Retrieval](#proteome-retrieval)

## Perform Meta-Genome Retrieval

The number of genome assemblies generated and stored in sequence databases is growing exponentially every year. With the availability of this growing amount of genomic data, meta-genomics studies become more and more popular. By using this bulk of genomes for comparing them to thousands of other genomes
new structural patterns and evolutionary insights can be obtained.
However, the first step in any meta-genomics study is the retrieval of the genomes, proteomes, coding sequences or annotation files that shall be compared
and investigated. For this purpose, the `meta.retrieval()` and `meta.retrieval.all()` functions allows users to perform straightforward meta-genome retrieval of hundreds of genomes, proteomes, CDS, etc in R. Finally, in addition to the retrieved sequence information the `meta.retrieval()` and `meta.retrieval.all()` functions will
generate a `summary file` which contains information about the genome version, genome status, submitter, etc for each organism to promote
computational and scientific reproducibility of the meta-genomics study at hand. This `summary file` can for example be attached as `Supplementary Data` to the respective study. 


#### Getting Started

The `meta.retrieval()` and `meta.retrieval.all()` functions aim to simplify the genome retrieval and computational reproducibility process for meta-genomics studies. Both functions allow users to either download genomes, proteomes, CDS, etc for species within a specific kingdom or subgroup of life (`meta.retrieval()`) or of all species of all kingdoms (`meta.retrieval.all()`). Before `biomartr` users had to write `shell` scripts to download respective genomic data. However, since many meta-genomics packages exist for the R programming language, I implemented this functionality for easy integration into existing R workflows and for easier reproducibility.

For example, the pipeline logic of the [magrittr](https://github.com/smbache/magrittr) package can be used with
`meta.retrieval()` and `meta.retrieval.all()` as follows.

```{r,eval=FALSE}
# download all vertebrate genomes, then apply ...
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "genome") %>% ...
```

Here `...` denotes any subsequent meta-genomics analysis. Hence, `meta.retrieval()` enables the pipeline methodology for meta-genomics. 

### Retrieve Genomic Sequences

To retrieve a list of all available kingdoms stored in the `NCBI RefSeq`,  `NCBI Genbank`, and `ENSEMBL` databases users can consult the `getKingdoms()` function which stores a list of all available kingdoms of life for the corresponding database. 

### Example `NCBI RefSeq`:

```{r, eval=FALSE}
getKingdoms(db = "refseq")
```

```
[1] "archaea"              "bacteria"             "fungi"                "invertebrate"        
[5] "plant"                "protozoa"             "vertebrate_mammalian" "vertebrate_other"    
[9] "viral"
```

### Example `NCBI Genbank`:

```{r, eval=FALSE}
getKingdoms(db = "genbank")
```

```
[1] "archaea"              "bacteria"             "fungi"               
[4] "invertebrate"         "plant"                "protozoa"            
[7] "vertebrate_mammalian" "vertebrate_other"
```

In these examples the difference betwenn `db = "refseq"` and `db = "genbank"` is that `db = "genbank"` does not store `viral` information.


### Example `ENSEMBL`

```{r, eval=FALSE}
getKingdoms(db = "ensembl")
```

```
[1] "EnsemblVertebrates"                                             
```

The `ENSEMBL` database does not differentiate between different kingdoms, but specialized on storing high-quality reference genomes of diverse biological disciplines.

#### Retrieval from `NCBI RefSeq`

Download all mammalian vertebrate genomes from `NCBI RefSeq`.

```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "genome", reference = FALSE)
```

The argument `kingdom` specifies the kingdom selected with `getKingdoms()` from which
genomes of organisms shall be retrieved. The `db` argument specifies the database
from which respective genomes shall be downloaded. The argument `type` specifies
that `genome assembly` files shall be retrieved. The argument `reference` indicates whether or not a genome shall be downloaded if it isn't marked in the database as either a `reference genome` or a `representative genome`. Options are:

- `reference = FALSE` (__Default__): all organisms (reference, representative, and non-representative genomes) are downloaded.
- `reference = TRUE`: organisms that are downloaded must be either a `reference` or `representative genome`. Thus, most genomes which are usually non-reference genomes will not be downloaded and the user will retrieve much less organisms than are stored in databases.

When running this command all geneomes are stored in a folder which is either named according to the kingdom
(in this case `vertebrate_mammalian`). Alternatively, users can specify
the `out.folder` argument to define a custom output folder path. 

Internally, in this example `meta.retrieval()` will generate a folder named `vertebrate_mammalian`
in which respective genomes will be stored. In addition, the `vertebrate_mammalian`
folder contains a folder named `documentation` which stores individual documentation
files for each individual organism and a `summary file` which stores documentation
for all retrieved organisms. This `summary file` can then be used as `Supplementary Data`
in studies to promote computational reproducibility.

An example documentation file of an individual organism looks like this:

```
File Name: Mus_musculus_genomic_genbank.gff.gz
Organism Name: Mus_musculus
Database: NCBI genbank
URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.7_GRCm38.p5/GCA_000001635.7_GRCm38.p5_genomic.gff.gz
Download_Date: Mon Nov 14 12:43:45 2016
refseq_category: reference genome
assembly_accession: GCA_000001635.7
bioproject: PRJNA20689
biosample: NA
taxid: 10090
infraspecific_name: NA
version_status: latest
release_type: Patch
genome_rep: Full
seq_rel_date: 2016-06-29
submitter: Genome Reference Consortium
```


An example `summary file` of all organism looks like this (here we use 105 Plant species as an example):

```
# A tibble: 105 x 16
   file_name     organism   url             database path  refseq_category
   <chr>         <chr>      <chr>           <chr>    <chr> <chr>          
 1 Aegilops_tau… Aegilops_… ftp://ftp.ncbi… refseq   Prot… representative…
 2 Amborella_tr… Amborella… ftp://ftp.ncbi… refseq   Prot… representative…
 3 Ananas_comos… Ananas_co… ftp://ftp.ncbi… refseq   Prot… representative…
 4 Arabidopsis_… Arabidops… ftp://ftp.ncbi… refseq   Prot… representative…
 5 Arabidopsis_… Arabidops… ftp://ftp.ncbi… refseq   Prot… reference geno…
 6 Arachis_dura… Arachis_d… ftp://ftp.ncbi… refseq   Prot… representative…
 7 Arachis_ipae… Arachis_i… ftp://ftp.ncbi… refseq   Prot… representative…
 8 Asparagus_of… Asparagus… ftp://ftp.ncbi… refseq   Prot… representative…
 9 Auxenochlore… Auxenochl… ftp://ftp.ncbi… refseq   Prot… representative…
10 Bathycoccus_… Bathycocc… ftp://ftp.ncbi… refseq   Prot… representative…
# ... with 95 more rows, and 10 more variables: assembly_accession <chr>,
#   bioproject <chr>, biosample <chr>, taxid <int>,
#   infraspecific_name <chr>, version_status <chr>, release_type <chr>,
#   genome_rep <chr>, seq_rel_date <date>, submitter <chr>
```

#### Restarting a corrupted download

Unfortunately, when downloading large amounts of genomes the NCBI RefSeq database
limits the number of queries from an individual IP address. This causes
that the download process might stop or time out at a particular step.
To overcome this limitation users can simply __re-run__ the `meta.retrieval()` command they previously executed 
and specify the argument `restart_at_last` which has the following two options:

- If `restart_at_last = TRUE` (__Default__) then `meta.retrieval()` will skip all organisms that are already present in the folder and will start downloading all remaining species (thus will pick up from where the initial download process stopped). However, this way `meta.retrieval()` will not be able to check whether already downloaded organism files are corrupted or not by checking the `md5 checksum` of the respective file. Thus, I recommend to download the last organism before `meta.retrieval()` stopped manually using `getGenome()` to make sure that the respective file is not corrupted.
- If `restart_at_last = FALSE` then `meta.retrieval()` will start from the beginning and crawl through already downloaded organism files and check whether already downloaded organism files are corrupted or not by checking the `md5 checksum` (this procedure takes longer than `restart_at_last = TRUE`). After checking existing files the function will start downloading all remaining organisms.

#### Un-zipping downloaded files

After downloading genomes users can format the output of `meta.retrieval()` by first un-zipping downloaded files and renaming them for more convenient downstream data analysis (e.g. from `Saccharomyces_cerevisiae_cds_from_genomic_refseq.fna.gz` to `Scerevisiae.fa`).

The easiest way to use `clean.retrieval()` in combination with `meta.retrieval()` is to use the pipe operator from the `magrittr` package:

```r
library(magrittr)
meta.retrieval(kingdom = "vertebrate_mammalian", 
               db = "refseq", 
               type = "genome") %>% 
    clean.retrieval()
```

In the first step, genome assembly files are downloaded with `meta.retrieval` and
subsequently (`%>%`) un-zipped and re-named using `clean.retrieval()`.



Example `Bacteria`
```{r,eval=FALSE}
# download all bacteria genomes
meta.retrieval(kingdom = "bacteria", db = "refseq", type = "genome", reference = FALSE)
```

Example `Viruses`
```{r,eval=FALSE}
# download all virus genomes
meta.retrieval(kingdom = "viral", db = "refseq", type = "genome", reference = FALSE)
```

Example `Archaea`
```{r,eval=FALSE}
# download all archaea genomes
meta.retrieval(kingdom = "archaea", db = "refseq", type = "genome", reference = FALSE)
```

Example `Fungi`
```{r,eval=FALSE}
# download all fungi genomes
meta.retrieval(kingdom = "fungi", db = "refseq", type = "genome", reference = FALSE)
```

Example `Plants`
```{r,eval=FALSE}
# download all plant genomes
meta.retrieval(kingdom = "plant", db = "refseq", type = "genome", reference = FALSE)
```

Example `Invertebrates`
```{r,eval=FALSE}
# download all invertebrate genomes
meta.retrieval(kingdom = "invertebrate", db = "refseq", type = "genome", reference = FALSE)
```

Example `Protozoa`
```{r,eval=FALSE}
# download all invertebrate genomes
meta.retrieval(kingdom = "protozoa", db = "refseq", type = "genome", reference = FALSE)
```

#### Retrieval from `NCBI Genbank`

Alternatively, download all mammalian vertebrate genomes from `NCBI Genbank`, e.g.

```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "genome", reference = FALSE)
```

Example `Bacteria`
```{r,eval=FALSE}
# download all bacteria genomes
meta.retrieval(kingdom = "bacteria", db = "genbank", type = "genome", reference = FALSE)
```

Example `Archaea`
```{r,eval=FALSE}
# download all archaea genomes
meta.retrieval(kingdom = "archaea", db = "genbank", type = "genome", reference = FALSE)
```

Example `Fungi`
```{r,eval=FALSE}
# download all fungi genomes
meta.retrieval(kingdom = "fungi", db = "genbank", type = "genome", reference = FALSE)
```

Example `Plants`
```{r,eval=FALSE}
# download all plant genomes
meta.retrieval(kingdom = "plant", db = "genbank", type = "genome", reference = FALSE)
```

Example `Invertebrates`
```{r,eval=FALSE}
# download all invertebrate genomes
meta.retrieval(kingdom = "invertebrate", db = "genbank", type = "genome", reference = FALSE)
```

Example `Protozoa`
```{r,eval=FALSE}
# download all invertebrate genomes
meta.retrieval(kingdom = "protozoa", db = "genbank", type = "genome", reference = FALSE)
```

#### Retrieval from `ENSEMBL`

```{r,eval=FALSE}
# download all genomes from ENSEMBL
meta.retrieval(kingdom = "Ensembl", db = "ensembl", type = "genome", reference = FALSE)
```

### Retrieve groups or subgroups of species 

In case users do not wish to retrieve genomes from an entire kingdom, but rather from
a group or subgoup (e.g. from species belonging to the `Gammaproteobacteria` class, a subgroup of the `bacteria` kingdom), they can use the following workflow.

#### Example retrieval of all `Gammaproteobacteria` genomes from `NCBI RefSeq`:

First, users can again consult the `getKingdoms()` function to retrieve kingdom information.

```{r, eval=FALSE}
getKingdoms(db = "refseq")
```

```
[1] "archaea"              "bacteria"             "fungi"                "invertebrate"        
[5] "plant"                "protozoa"             "vertebrate_mammalian" "vertebrate_other"    
[9] "viral"
```

In this example, we will choose the `bacteria` kingdom. Now, the `getGroups()` function allows users to obtain
available subgroups of the `bacteria` kingdom.

```{r, eval=FALSE}
getGroups(db = "refseq", kingdom = "bacteria")
```

```
 [1] "Acidithiobacillia"                     "Acidobacteriia"                       
 [3] "Actinobacteria"                        "Alphaproteobacteria"                  
 [5] "Aquificae"                             "Armatimonadetes"                      
 [7] "Bacteroidetes/Chlorobi group"          "Balneolia"                            
 [9] "Betaproteobacteria"                    "Blastocatellia"                       
[11] "Candidatus Kryptonia"                  "Chlamydiae"                           
[13] "Chloroflexi"                           "Cyanobacteria/Melainabacteria group"  
[15] "Deinococcus-Thermus"                   "delta/epsilon subdivisions"           
[17] "Endomicrobia"                          "Fibrobacteres"                        
[19] "Firmicutes"                            "Fusobacteriia"                        
[21] "Gammaproteobacteria"                   "Gemmatimonadetes"                     
[23] "Kiritimatiellaeota"                    "Nitrospira"                           
[25] "Planctomycetes"                        "Spirochaetia"                         
[27] "Synergistia"                           "Tenericutes"                          
[29] "Thermodesulfobacteria"                 "Thermotogae"                          
[31] "unclassified Acidobacteria"            "unclassified Bacteria (miscellaneous)"
[33] "unclassified Proteobacteria"           "Verrucomicrobia"                      
[35] "Zetaproteobacteria" 
```

Please note, that the `kingdom` argument specified in `getGroups()` needs to match with an available kingdom retrieved with `getKingdoms()`.
It is also important that in both cases: `getKingdoms()` and `getGroups()` the same database should be specified.


Now we choose the group `Gammaproteobacteria` and specify the `group` argument in the `meta.retrieval()` function.

```{r, eval=FALSE}
meta.retrieval(kingdom = "bacteria", group = "Gammaproteobacteria", db = "refseq", type = "genome", reference = FALSE)
```

Using this command, all bacterial (`kingdom = "bacteria"`) genomes (`type = "genome"`) that belong to the group `Gammaproteobacteria` (`group = "Gammaproteobacteria"`) will be retrieved from NCBI RefSeq (`db = "refseq"`).

Alternatively, `Gammaproteobacteria` genomes can be retrieved from NCBI Genbank by exchanging `db = "refseq"` to `db = "genbank"`.
If users wish to download proteome, CDS, or GFF files instead of genomes, they can specify the argument: `type = "proteome"`, `type = "cds"`, or `type = "gff"`.

#### Example retrieval of all `Adenoviridae` genomes from `NCBI RefSeq`:

Retrieve groups for viruses.
```{r, eval=FALSE}
getGroups(db = "refseq", kingdom = "viral")
```

```
 [1] "Adenoviridae"                                        "Alloherpesviridae"                                  
  [3] "Alphaflexiviridae"                                   "Alphatetraviridae"                                  
  [5] "Alvernaviridae"                                      "Amalgaviridae"                                      
  [7] "Ampullaviridae"                                      "Anelloviridae"                                      
  [9] "Apple fruit crinkle viroid"                          "Apple hammerhead viroid-like circular RNA"          
 [11] "Apscaviroid"                                         "Arenaviridae"                                       
 [13] "Arteriviridae"                                       "Ascoviridae"                                        
 [15] "Asfarviridae"                                        "Astroviridae"                                       
 [17] "Avsunviroid"                                         "Baculoviridae"                                      
 [19] "Barnaviridae"                                        "Benyviridae"                                        
 [21] "Betaflexiviridae"                                    "Bicaudaviridae"                                     
 [23] "Birnaviridae"                                        "Bornaviridae"                                       
 [25] "Bromoviridae"                                        "Bunyaviridae"                                       
 [27] "Caliciviridae"                                       "Carmotetraviridae"                                  
 [29] "Caulimoviridae"                                      "Cherry leaf scorch small circular viroid-like RNA 1"
 [31] "Cherry small circular viroid-like RNA"               "Chrysoviridae"                                      
 [33] "Circoviridae"                                        "Closteroviridae"                                    
 [35] "Cocadviroid"                                         "Coleviroid"                                         
 [37] "Coronaviridae"                                       "Corticoviridae"                                     
 [39] "Cystoviridae"                                        "Dicistroviridae"                                    
 [41] "Elaviroid"                                           "Endornaviridae"                                     
 [43] "Filoviridae"                                         "Flaviviridae"                                       
 [45] "Fusarividae"                                         "Fuselloviridae"                                     
 [47] "Gammaflexiviridae"                                   "Geminiviridae"                                      
 [49] "Genomoviridae"                                       "Globuloviridae"                                     
 [51] "Grapevine latent viroid"                             "Guttaviridae"                                       
 [53] "Hepadnaviridae"                                      "Hepeviridae"                                        
 [55] "Herpesviridae"                                       "Hostuviroid"                                        
 [57] "Hypoviridae"                                         "Hytrosaviridae"                                     
 [59] "Iflaviridae"                                         "Inoviridae"                                         
 [61] "Iridoviridae"                                        "Lavidaviridae"                                      
 [63] "Leviviridae"                                         "Lipothrixviridae"                                   
 [65] "Luteoviridae"                                        "Malacoherpesviridae"                                
 [67] "Marnaviridae"                                        "Marseilleviridae"                                   
 [69] "Megabirnaviridae"                                    "Mesoniviridae"                                      
 [71] "Microviridae"                                        "Mimiviridae"                                        
 [73] "Mulberry small circular viroid-like RNA 1"           "Mymonaviridae"                                      
 [75] "Myoviridae"                                          "Nanoviridae"                                        
 [77] "Narnaviridae"                                        "Nimaviridae"                                        
 [79] "Nodaviridae"                                         "Nudiviridae"                                        
 [81] "Nyamiviridae"                                        "Ophioviridae"                                       
 [83] "Orthomyxoviridae"                                    "Other"                                              
 [85] "Papillomaviridae"                                    "Paramyxoviridae"                                    
 [87] "Partitiviridae"                                      "Parvoviridae"                                       
 [89] "Pelamoviroid"                                        "Permutotetraviridae"                                
 [91] "Persimmon viroid"                                    "Phycodnaviridae"                                    
 [93] "Picobirnaviridae"                                    "Picornaviridae"                                     
 [95] "Plasmaviridae"                                       "Pneumoviridae"                                      
 [97] "Podoviridae"                                         "Polydnaviridae"                                     
 [99] "Polyomaviridae"                                      "Pospiviroid"                                        
[101] "Potyviridae"                                         "Poxviridae"                                         
[103] "Quadriviridae"                                       "Reoviridae"                                         
[105] "Retroviridae"                                        "Rhabdoviridae"                                      
[107] "Roniviridae"                                         "Rubber viroid India/2009"                           
[109] "Rudiviridae"                                         "Secoviridae"                                        
[111] "Siphoviridae"                                        "Sphaerolipoviridae"                                 
[113] "Sunviridae"                                          "Tectiviridae"                                       
[115] "Togaviridae"                                         "Tombusviridae"                                      
[117] "Totiviridae"                                         "Turriviridae"                                       
[119] "Tymoviridae"                                         "unclassified"                                       
[121] "unclassified Pospiviroidae"                          "Virgaviridae"
```

Now we can choose `Adenoviridae` as group argument for the `meta.retrieval()` function.

```{r, eval=FALSE}
meta.retrieval(kingdom = "viral", group = "Adenoviridae", db = "refseq", type = "genome", reference = FALSE)
```

Again, by exchanging `type = "genome"` by either `type = "proteome"`, `type = "cds"`, `type = "rna"`, `type = "assemblystats"`, or `type = "gff"`, users can retrieve proteome, CDS, RNA, genome assembly statistics or GFF files instead of genomes.


### Meta retrieval of genome assembly quality information

Although much effort is invested to increase the genome assembly quality
when new genomes are published or new versions are released, the influence
of genome assembly quality on downstream analyses cannot be neglected.
A rule of thumb is, that the larger the genome the more prone it is to
genome assembly errors and therefore, a reduction of assembly quality.

In [Veeckman et al., 2016](http://www.plantcell.org/content/28/8/1759.short) the authors conclude:

> As yet, no uniform metrics or standards are in place to estimate the completeness of a genome assembly or
> the annotated gene space, despite their importance for downstream analyses

In most metagenomics studies, however, the influence or bias of genome assembly quality on the 
outcome of the analysis (e.g. comparative genomics, annotation, etc.) is neglected. To better 
grasp the genome assembly quality, the NCBI databases store genome assembly statistics of some
species for which genome assemblies are available. An example assembly statistics report can
be found at: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_stats.txt. 

The `biomartr` package allows users to retrieve these genome assembly stats file in an automated way by specifying the argument `type = "assemblystats"` and `combine = TRUE`.
Please make sure that `combine = TRUE` when selecting `type = "assemblystats"`.

```{r, eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# retrieve genome assembly stats for all mammal genome assemblies
# and store these stats in a data.frame
mammals.gc <- meta.retrieval(kingdom = "vertebrate_mammalian", 
                             db      = "refseq", 
                             type    = "assemblystats", 
                             combine = TRUE)

mammals.gc
```

```
                    species total_length spanned_gaps unspanned_gaps region_count scaffold_count
                      <chr>        <int>        <int>          <int>        <int>          <int>
1  Ornithorhynchus anatinus   1995607322       243698            137            0         200283
2      Sarcophilus harrisii           NA       201317              0            0          35974
3      Dasypus novemcinctus           NA       268413              0            0          46559
4       Erinaceus europaeus           NA       219764              0            0           5803
5         Echinops telfairi           NA       269444              0            0           8402
6           Pteropus alecto   1985975446       104566              0            0          65598
7     Rousettus aegyptiacus   1910250568          559              0            0             NA
8        Callithrix jacchus           NA       184972           2242            0          16399
9  Cebus capucinus imitator           NA       133441              0            0           7156
10          Cercocebus atys           NA        65319              0            0          11433
# ... with 89 more rows, and 9 more variables: scaffold_N50 <int>, scaffold_L50 <int>,
#   scaffold_N75 <int>, scaffold_N90 <int>, contig_count <int>, contig_N50 <int>, total_gap_length <int>,
#   molecule_count <int>, top_level_count <int>
```

Analogously, this information can be retrieved for each kingdom other than `kingdom = "vertebrate_mammalian"`. Please consult `getKingdoms()` for available kingdoms.

### Metagenome project retrieval from NCBI Genbank

NCBI Genbank stores [metagenome projects](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/metagenomes/) in addition to species specific genome, proteome or CDS sequences. To retrieve these metagenomes users can perform the following combination of commands:

First, users can list the project names of available metagenomes by typing

```{r,eval=FALSE}
# list available metagenomes at NCBI Genbank
listMetaGenomes()
```

```
[1] "metagenome"                     "human gut metagenome"           "epibiont metagenome"           
 [4] "marine metagenome"              "soil metagenome"                "mine drainage metagenome"      
 [7] "mouse gut metagenome"           "marine sediment metagenome"     "termite gut metagenome"        
[10] "hot springs metagenome"         "human lung metagenome"          "fossil metagenome"             
[13] "freshwater metagenome"          "saltern metagenome"             "stromatolite metagenome"       
[16] "coral metagenome"               "mosquito metagenome"            "fish metagenome"               
[19] "bovine gut metagenome"          "chicken gut metagenome"         "wastewater metagenome"         
[22] "microbial mat metagenome"       "freshwater sediment metagenome" "human metagenome"              
[25] "hydrothermal vent metagenome"   "compost metagenome"             "wallaby gut metagenome"        
[28] "groundwater metagenome"         "gut metagenome"                 "sediment metagenome"           
[31] "ant fungus garden metagenome"   "food metagenome"                "hypersaline lake metagenome"   
[34] "hydrocarbon metagenome"         "activated sludge metagenome"    "viral metagenome"              
[37] "bioreactor metagenome"          "wasp metagenome"                "permafrost metagenome"         
[40] "sponge metagenome"              "aquatic metagenome"             "insect gut metagenome"         
[43] "activated carbon metagenome"    "anaerobic digester metagenome"  "rock metagenome"               
[46] "terrestrial metagenome"         "rock porewater metagenome"      "seawater metagenome"           
[49] "scorpion gut metagenome"        "soda lake metagenome"           "glacier metagenome"
```

Internally the `listMetaGenomes()` function downloads the assembly_summary.txt file from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/metagenomes/ to retrieve
available metagenome information. This procedure might take a few seconds during the first run of `listMetaGenomes()`. Subsequently, the assembly_summary.txt file will be stored in the `tempdir()` directory to
achieve a much faster access of this information during following uses of `listMetaGenomes()`.

In case users wish to retrieve detailed information about available metagenome projects they can specify the `details = TRUE` argument.

```{r,eval=FALSE}
# show all elements of the data.frame
options(tibble.print_max = Inf)
# detailed information on available metagenomes at NCBI Genbank
listMetaGenomes(details = TRUE)
```

```
# A tibble: 857 x 21
   assembly_accession bioproject    biosample     wgs_master refseq_category  taxid species_taxid
                <chr>      <chr>        <chr>          <chr>           <chr>  <int>         <int>
1     GCA_000206185.1 PRJNA32359 SAMN02954317 AAGA00000000.1              na 256318        256318
2     GCA_000206205.1 PRJNA32355 SAMN02954315 AAFZ00000000.1              na 256318        256318
3     GCA_000206225.1 PRJNA32357 SAMN02954316 AAFY00000000.1              na 256318        256318
4     GCA_000208265.2 PRJNA17779 SAMN02954240 AASZ00000000.1              na 256318        256318
5     GCA_000208285.1 PRJNA17657 SAMN02954268 AATO00000000.1              na 256318        256318
6     GCA_000208305.1 PRJNA17659 SAMN02954269 AATN00000000.1              na 256318        256318
7     GCA_000208325.1 PRJNA16729 SAMN02954263 AAQL00000000.1              na 256318        256318
8     GCA_000208345.1 PRJNA16729 SAMN02954262 AAQK00000000.1              na 256318        256318
9     GCA_000208365.1 PRJNA13699 SAMN02954283 AAFX00000000.1              na 256318        256318
10    GCA_900010595.1 PRJEB11544 SAMEA3639840 CZPY00000000.1              na 256318        256318
# ... with 847 more rows, and 14 more variables: organism_name <chr>, infraspecific_name <chr>,
#   isolate <chr>, version_status <chr>, assembly_level <chr>, release_type <chr>, genome_rep <chr>,
#   seq_rel_date <date>, asm_name <chr>, submitter <chr>, gbrs_paired_asm <chr>, paired_asm_comp <chr>,
#   ftp_path <chr>, excluded_from_refseq <chr>
```

Finally, users can retrieve available metagenomes using `getMetaGenomes()`. The `name`
argument receives the metagenome project name retrieved with `listMetaGenomes()`.
The `path` argument specifies the folder path in which corresponding genomes shall be stored.

```{r,eval=FALSE}
# retrieve all genomes belonging to the human gut metagenome project
getMetaGenomes(name = "human gut metagenome", path = file.path("_ncbi_downloads","human_gut"))
```

```
1] "The metagenome of 'human gut metagenome' has been downloaded to '_ncbi_downloads/human_gut'."
  [1] "_ncbi_downloads/human_gut/GCA_000205525.2_ASM20552v2_genomic.fna.gz"
  [2] "_ncbi_downloads/human_gut/GCA_000205765.1_ASM20576v1_genomic.fna.gz"
  [3] "_ncbi_downloads/human_gut/GCA_000205785.1_ASM20578v1_genomic.fna.gz"
  [4] "_ncbi_downloads/human_gut/GCA_000207925.1_ASM20792v1_genomic.fna.gz"
  [5] "_ncbi_downloads/human_gut/GCA_000207945.1_ASM20794v1_genomic.fna.gz"
  [6] "_ncbi_downloads/human_gut/GCA_000207965.1_ASM20796v1_genomic.fna.gz"
  [7] "_ncbi_downloads/human_gut/GCA_000207985.1_ASM20798v1_genomic.fna.gz"
  [8] "_ncbi_downloads/human_gut/GCA_000208005.1_ASM20800v1_genomic.fna.gz"
  [9] "_ncbi_downloads/human_gut/GCA_000208025.1_ASM20802v1_genomic.fna.gz"
 [10] "_ncbi_downloads/human_gut/GCA_000208045.1_ASM20804v1_genomic.fna.gz"
 [11] "_ncbi_downloads/human_gut/GCA_000208065.1_ASM20806v1_genomic.fna.gz"
 [12] "_ncbi_downloads/human_gut/GCA_000208085.1_ASM20808v1_genomic.fna.gz"
 [13] "_ncbi_downloads/human_gut/GCA_000208105.1_ASM20810v1_genomic.fna.gz"
 [14] "_ncbi_downloads/human_gut/GCA_000208125.1_ASM20812v1_genomic.fna.gz"
 [15] "_ncbi_downloads/human_gut/GCA_000208145.1_ASM20814v1_genomic.fna.gz"
 [16] "_ncbi_downloads/human_gut/GCA_000208165.1_ASM20816v1_genomic.fna.gz"
 ...
```

Internally, `getMetaGenomes()` creates a folder specified in the `path` argument.
Genomes associated with the metagenomes project specified in the `name` argument
will then be downloaded and stored in this folder. As return value `getMetaGenomes()`
returns the file paths to the genome files which can then be used as input to the `read*()` functions.

Alternatively or subsequent to the metagenome retrieval, users can retrieve annotation files of genomes belonging to a metagenome project
selected with `listMetaGenomes()` by using the `getMetaGenomeAnnotations()` function.


```{r,eval=FALSE}
# retrieve all genomes belonging to the human gut metagenome project
getMetaGenomeAnnotations(name = "human gut metagenome", path = file.path("_ncbi_downloads","human_gut","annotations"))
```

```
[1] "The annotations of metagenome 'human gut metagenome' have been downloaded and stored at '_ncbi_downloads/human_gut/annotations'."
  [1] "_ncbi_downloads/human_gut/annotations/GCA_000205525.2_ASM20552v2_genomic.gff.gz"
  [2] "_ncbi_downloads/human_gut/annotations/GCA_000205765.1_ASM20576v1_genomic.gff.gz"
  [3] "_ncbi_downloads/human_gut/annotations/GCA_000205785.1_ASM20578v1_genomic.gff.gz"
  [4] "_ncbi_downloads/human_gut/annotations/GCA_000207925.1_ASM20792v1_genomic.gff.gz"
  [5] "_ncbi_downloads/human_gut/annotations/GCA_000207945.1_ASM20794v1_genomic.gff.gz"
  [6] "_ncbi_downloads/human_gut/annotations/GCA_000207965.1_ASM20796v1_genomic.gff.gz"
  [7] "_ncbi_downloads/human_gut/annotations/GCA_000207985.1_ASM20798v1_genomic.gff.gz"
  [8] "_ncbi_downloads/human_gut/annotations/GCA_000208005.1_ASM20800v1_genomic.gff.gz"
  [9] "_ncbi_downloads/human_gut/annotations/GCA_000208025.1_ASM20802v1_genomic.gff.gz"
 [10] "_ncbi_downloads/human_gut/annotations/GCA_000208045.1_ASM20804v1_genomic.gff.gz"
 [11] "_ncbi_downloads/human_gut/annotations/GCA_000208065.1_ASM20806v1_genomic.gff.gz"
 [12] "_ncbi_downloads/human_gut/annotations/GCA_000208085.1_ASM20808v1_genomic.gff.gz"
 [13] "_ncbi_downloads/human_gut/annotations/GCA_000208105.1_ASM20810v1_genomic.gff.gz"
 [13] "_ncbi_downloads/human_gut/annotations/GCA_000208105.1_ASM20810v1_genomic.gff.gz"
 [14] "_ncbi_downloads/human_gut/annotations/GCA_000208125.1_ASM20812v1_genomic.gff.gz"
 [15] "_ncbi_downloads/human_gut/annotations/GCA_000208145.1_ASM20814v1_genomic.gff.gz"
 [16] "_ncbi_downloads/human_gut/annotations/GCA_000208165.1_ASM20816v1_genomic.gff.gz"
 ...
```

The file paths of the downloaded `*.gff` are retured by `getMetaGenomeAnnotations()` and can be used
as input for the `read.gff()` function in the [seqreadr](https://github.com/HajkD/seqreadr) package.


#### Retrieve Protein Sequences

Download all mammalian vertebrate proteomes.

#### Retrieval from `NCBI RefSeq`:
```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "proteome", reference = FALSE)
```


#### Retrieval from `NCBI Genbank`:
```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "proteome", reference = FALSE)
```


#### Retrieval from `ENSEMBL`:

```{r,eval=FALSE}
# download all Ensembl proteome sequneces
meta.retrieval(kingdom = "Ensembl", db = "ensembl", type = "proteome", reference = FALSE)
```

#### Retrieve CDS Sequences

Download all mammalian vertebrate CDS from RefSeq (Genbank does not store CDS data).

#### Retrieval from `NCBI RefSeq`:
```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "cds", reference = FALSE)
```

#### Retrieval from `NCBI Genbank`:
```{r,eval=FALSE}
# download all vertebrate genomes
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "cds", reference = FALSE)
```

#### Retrieval from `ENSEMBL`:

```{r,eval=FALSE}
# download all Ensembl CDS sequneces
meta.retrieval(kingdom = "Ensembl", db = "ensembl", type = "cds", reference = FALSE)
```

### Retrieve GFF files

Download all mammalian vertebrate gff files.

Example `NCBI RefSeq`:
```{r,eval=FALSE}
# download all vertebrate gff files
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "gff", reference = FALSE)
```


Example `NCBI Genbank`:
```{r,eval=FALSE}
# download all vertebrate gff files
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "gff", reference = FALSE)
```

### Retrieve GTF files

Download all mammalian vertebrate gtf files.

Example `ENSEMBL`:
```{r,eval=FALSE}
# download all vertebrate gff files
meta.retrieval(kingdom = "Ensembl", db = "ensembl", type = "gtf", reference = FALSE)
```


### Retrieve RNA sequences

Download all mammalian vertebrate RNA sequences from `NCBI RefSeq` and `NCBI Genbank`.

#### Retrieval from `NCBI RefSeq`:
```{r,eval=FALSE}
# download all vertebrate RNA sequneces
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "rna", reference = FALSE)
```

#### Retrieval from `NCBI Genbank`:
```{r,eval=FALSE}
# download all vertebrate RNA sequneces
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "rna", reference = FALSE)
```

#### Retrieval from `ENSEMBL`:

```{r,eval=FALSE}
# download all Ensembl RNA sequneces
meta.retrieval(kingdom = "Ensembl", db = "ensembl", type = "rna", reference = FALSE)
```

### Retrieve Repeat Masker Sequences

Download all mammalian vertebrate Repeat Masker Annotation files from `NCBI RefSeq` and `NCBI Genbank`.

#### Retrieval from `NCBI RefSeq`:
```{r,eval=FALSE}
# download all vertebrate RNA sequneces
meta.retrieval(kingdom = "vertebrate_mammalian", db = "refseq", type = "rm", reference = FALSE)
```

#### Retrieval from `NCBI Genbank`:
```{r,eval=FALSE}
# download all vertebrate RNA sequneces
meta.retrieval(kingdom = "vertebrate_mammalian", db = "genbank", type = "rm", reference = FALSE)
```

Users can obtain alternative kingdoms using `getKingdoms()`.

## Retrieve Individual Genomes for all Species in the Tree of Life

If users wish to download the all genomes, proteome, CDS, or gff files for all species
available in RefSeq or Genbank, they can use the `meta.retrieval.all()` function for this purpose.

### Genome Retrieval

Example `RefSeq`:
```{r,eval=FALSE}
# download all geneomes stored in RefSeq
meta.retrieval.all(db = "refseq", type = "genome", reference = FALSE)
```

Example `Genbank`:
```{r,eval=FALSE}
# download all geneomes stored in Genbank
meta.retrieval.all(db = "genbank", type = "genome", reference = FALSE)
```

### Proteome Retrieval

Example `RefSeq`:
```{r,eval=FALSE}
# download all proteome stored in RefSeq
meta.retrieval.all(db = "refseq", type = "proteome", reference = FALSE)
```

Example `Genbank`:
```{r,eval=FALSE}
# download all proteome stored in Genbank
meta.retrieval.all(db = "genbank", type = "proteome", reference = FALSE)
```

Again, by exchanging `type = "proteome"` by either 

- `type = "genome"`
- `type = "cds"`
- `type = "rna"`
- `type = "assemblystats"`
- `type = "gff"`

users can retrieve genome, CDS, RNA, genome assembly statistics or GFF files instead of proteomes.

---
title: "Sequence Retrieval"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Sequence Retrieval}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```


## Biological Sequence Retrieval 

The `biomartr` package allows users to retrieve biological sequences in a very simple and intuitive way.

Using `biomartr`, users can retrieve either genomes, proteomes, CDS, RNA, GFF, and genome assembly statistics data using the specialized functions:


> **_NOTE:_** To make sure that you have a sufficiently stable (internet) connection between R and the respective databases, please set the default `timeout` setting __on your local machine__ from 60sec to at least 30000sec before running any retrieval functions via:

```r
options(timeout = 30000)
```

## Topics

* [0. Getting Started with Sequence Retrieval](#getting-started-with-sequence-retrieval)
    - [0.1 Listing the total number of available genomes](#listing-the-total-number-of-available-genomes)
* [1. Genome retrieval with `getGenome()`](#genome-retrieval)
    - [A. list total number of available genomes in each database](#listing-the-total-number-of-available-genomes)
    - [B. dealing with kingdoms, groups, and subgroups ](#retrieving-kingdom-group-and-subgroup-information)
    - [1.1 from NCBI RefSeq](#example-ncbi-refseq)
    - [1.2 from NCBI Genbank](#example-ncbi-genbank)
    - [1.3 from ENSEMBL](#example-ensembl)
* [1.B. Multiple genome retrieval with `getGenomeSet()`](#genomeset-retrieval)
* [2. Proteome retrieval with `getProteome()`](#proteome-retrieval)
    - [2.1 from NCBI RefSeq](#example-ncbi-refseq-1)
    - [2.2 from NCBI Genbank](#example-ncbi-genbank-1)
    - [2.3 from ENSEMBL](#example-ensembl-1)
    - [2.5 from UniProt](#example-retrieval-uniprot)
* [2.B. Multiple proteome retrieval with `getProteomeSet()`](#proteomeset-retrieval)
* [3. Coding sequence retrieval with `getCDS()`](#cds-retrieval)
    - [3.1 from NCBI RefSeq](#example-ncbi-refseq-2)
    - [3.2 from NCBI Genbank](#example-ncbi-genbank-2)
    - [3.3 from ENSEMBL](#example-ensembl-2)
* [3.B. Multiple coding sequence retrieval with `getCDSSet()`](#cdsset-retrieval)
* [4. RNA retrieval with `getRNA()`](#rna-retrieval)
    - [4.1 from NCBI RefSeq](#example-ncbi-refseq-3)
    - [4.2 from NCBI Genbank](#example-ncbi-genbank-3)
    - [4.3 from ENSEMBL](#example-ensembl-3)
* [4.B. Multiple RNA retrieval with `getRNASet()`](#rnaset-retrieval)
* [5. GFF retrieval with `getGFF()`](#retrieve-the-annotation-file-of-a-particular-genome)
    - [5.1 from NCBI RefSeq](#example-ncbi-refseq-4)
    - [5.2 from NCBI Genbank](#example-ncbi-genbank-4)
    - [5.3 from ENSEMBL](#example-ensembl-4)
* [5.B. Multiple GFF retrieval with `getGFFSet()`](#gffset-retrieval)
* [6. GTF retrieval with `getGTF()`](#retrieve-the-annotation-file-of-a-particular-genome)
    - [6.1 from NCBI RefSeq](#example-ncbi-refseq-5)
    - [6.2 from NCBI Genbank](#example-ncbi-genbank-5)
    - [6.3 from ENSEMBL](#example-ensembl-5)
* [7. Repeat Masker Annotation file retrieval with `getRepeatMasker()`](#repeat-masker-retrieval)
    - [7.1 from NCBI RefSeq](#example-ncbi-refseq-6)
* [8. Genome Assembly Stats Retrieval with `getAssemblyStats()`](#genome-assembly-stats-retrieval)
    - [8.1 from NCBI RefSeq](#example-ncbi-refseq-7)
    - [8.2 from NCBI Genbank](#example-ncbi-genbank-6)
* [9. Collection (Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats) retrieval with `getCollection()`](#collection-retrieval)
    - [9.1 from NCBI RefSeq](#example-ncbi-refseq-8)
    - [9.2 from NCBI Genbank](#example-ncbi-genbank-7)
    - [9.3 from NCBI ENSEMBL](#example-ensembl-6)

## Getting Started with Sequence Retrieval

First users can check whether or not the genome, proteome, CDS, RNA, GFF, GTF, or genome assembly statistics of their interest is available for download.

Using the scientific name of the organism of interest, users can check whether the 
corresponding genome is available via the `is.genome.available()` function.

Please note that the first time you run this command it might take a while, because
during the initial execution of this function all necessary information is retrieved from NCBI
and then stored locally. All further runs are then much faster.

### Example `NCBI RefSeq` (?is.genome.available): 

```r
# checking whether or not the Homo sapiens 
# genome is avaialable for download
is.genome.available(db = "refseq", organism = "Homo sapiens")
```

```
A reference or representative genome assembly is available for 'Homo sapiens'.
[1] TRUE
```

In the case of the human genome, there are more than one entry in the NCBI RefSeq
database (see message). By specifying the argument 'details = TRUE' users can
retrieve all information for 'Homo sapiens' that is stored in NCBI RefSeq.

```r
# checking whether or not the Homo sapiens 
# genome is avaialable for download
is.genome.available(db = "refseq", organism = "Homo sapiens", details = TRUE)
```

```
  assembly_accessi… bioproject biosample wgs_master refseq_category taxid
  <chr>             <chr>      <chr>     <chr>      <chr>           <int>
1 GCF_000001405.38  PRJNA168   NA        NA         reference geno…  9606
# ... with 15 more variables: species_taxid <int>, organism_name <chr>,
#   infraspecific_name <chr>, isolate <chr>, version_status <chr>,
#   assembly_level <chr>, release_type <chr>, genome_rep <chr>,
#   seq_rel_date <date>, asm_name <chr>, submitter <chr>,
#   gbrs_paired_asm <chr>, paired_asm_comp <chr>, ftp_path <chr>,
#   excluded_from_refseq <chr>
```

Here, we find one possible versions of the human genome having the `assembly_accession`
ID `GCF_000001405.38`.

When retrieving a genome with e.g. `getGenome()` the `organism` argument can
also be specified using the `assembly_accession` ID instead of the scientific name of the 
organism. This is true for all `get*()` functions. Hence, instead of writing
`getGenome(db = "refseq", organism = "Homo sapiens")`, users can specify 
`getGenome(db = "refseq", organism = "GCF_000001405.38")`. 
This is particularly useful when more than one entry is available for one organism.

__Please note that the `assembly_accession` id might change internally in external databases such as NCBI RefSeq. Thus when
writing automated scripts using the `assembly_accession` id they might stop working due to the
internal change of ids at RefSeq. Recently, I had the case where the `assembly_accession` id was changed
at RefSeq from `GCF_000001405.37` to `GCF_000001405.38`. Thus, scripts based on screening for entries with `is.genome.available()`
stopped working, because the id `GCF_000001405.37` couldn't be found anymore.__

In some cases users will find several entries for the same scientific name.
This might be due to the fact that ecotypes, strains, or other types of sub-species
are available in the respective databases.

For example. when we search for the bacterium `Mycobacterium tuberculosis` in NCBI RefSeq we get 5377 hits.


```r
is.genome.available(organism = "Mycobacterium tuberculosis", db = "refseq", details = TRUE)
```

```
A tibble: 6,744 x 21
   assembly_accessi… bioproject  biosample  wgs_master  refseq_category
   <chr>             <chr>       <chr>      <chr>       <chr>          
 1 GCF_000729745.1   PRJNA224116 SAMN02899… JPFP000000… na             
 2 GCF_000729755.1   PRJNA224116 SAMN02899… JPFQ000000… na             
 3 GCF_000729765.1   PRJNA224116 SAMN02899… JPFR000000… na             
 4 GCF_000749605.1   PRJNA224116 SAMN02902… JQES000000… na             
 5 GCF_000749615.1   PRJNA224116 SAMN02902… JQER000000… na             
 6 GCF_000749625.1   PRJNA224116 SAMN02902… JQEQ000000… na             
 7 GCF_000749665.1   PRJNA224116 SAMN02902… JQEV000000… na             
 8 GCF_000749675.1   PRJNA224116 SAMN02902… JQET000000… na             
 9 GCF_000749685.1   PRJNA224116 SAMN02902… JQEN000000… na             
10 GCF_000749725.1   PRJNA224116 SAMN02902… JQEM000000… na             
 … with 6,734 more rows, and 16 more variables: taxid <int>,
   species_taxid <int>, organism_name <chr>,
   infraspecific_name <chr>, isolate <chr>, version_status <chr>,
   assembly_level <chr>, release_type <chr>, genome_rep <chr>,
   seq_rel_date <date>, asm_name <chr>, submitter <chr>,
   gbrs_paired_asm <chr>, paired_asm_comp <chr>, ftp_path <chr>,
   excluded_from_refseq <chr>
```

When looking at the names of these organisms we see that they consist of different strains of `Mycobacterium tuberculosis`.

```r
tail(is.genome.available(organism = "Mycobacterium tuberculosis", 
                         db = "refseq", 
                         details = TRUE)$organism_name)
```

```
[1] "Mycobacterium tuberculosis CDC1551"     
[2] "Mycobacterium tuberculosis H37Ra"       
[3] "Mycobacterium tuberculosis F11"         
[4] "Mycobacterium tuberculosis KZN 1435"    
[5] "Mycobacterium tuberculosis SUMu001"     
[6] "Mycobacterium tuberculosis str. Haarlem"
```

Now we can use the `assembly_accession` id to retrieve the `Mycobacterium tuberculosis`
strain we are interested in, e.g. `Mycobacterium tuberculosis CDC1551`.

```r
MtbCDC1551 <- getGenome(db = "refseq", 
                        organism = "GCF_000008585.1",  
                        path   = file.path("_ncbi_downloads","genomes"), 
                        reference = FALSE)
                        
MtbCDC1551_genome <- read_genome(MtbCDC1551)

MtbCDC1551_genome
```

```
A DNAStringSet instance of length 1
      width seq                             names               
[1] 4403837 TTGACCGATGACCC...AGGGAGATACGTCG NC_002755.2 Mycob...
```

#### Using the NCBI Taxonomy ID instead of the scientific name to screen for organism availability 

Instead of specifying the `scientific name` of the organism of interest users can specify the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) identifier (= `taxid`)
of the corresponding organism. For example, the `taxid` of `Homo sapiens` is `9606`.
Now users can specify `organism = "9606"` to retrieve entries for `Homo sapiens`:

```r
# checking availability for Homo sapiens using its taxid
is.genome.available(db = "refseq", organism = "9606", details = TRUE)
```

```
  assembly_accession bioproject biosample wgs_master refseq_category  taxid
  <chr>              <chr>      <chr>     <chr>      <chr>            <int>
1 GCF_000001405.38   PRJNA168   NA        NA         reference genome  9606
# ... with 15 more variables: species_taxid <int>, organism_name <chr>,
#   infraspecific_name <chr>, isolate <chr>, version_status <chr>,
#   assembly_level <chr>, release_type <chr>, genome_rep <chr>,
#   seq_rel_date <date>, asm_name <chr>, submitter <chr>,
#   gbrs_paired_asm <chr>, paired_asm_comp <chr>, ftp_path <chr>,
#   excluded_from_refseq <chr>
```

#### Using the accession ID instead of the scientific name or taxid to screen for organism availability

Finally, instead of specifying either the `scientific name` of the organism of interest nor the `taxid`,
users can specify the accession ID of the organism of interest.
In the following example we use the accession ID of `Homo sapiens` (= `GCF_000001405.38`):
```r
# checking availability for Homo sapiens using its taxid
is.genome.available(db = "refseq", organism = "GCF_000001405.38", details = TRUE)
```

```
assembly_accession bioproject biosample wgs_master refseq_category  taxid
  <chr>              <chr>      <chr>     <chr>      <chr>            <int>
1 GCF_000001405.38   PRJNA168   NA        NA         reference genome  9606
# ... with 15 more variables: species_taxid <int>, organism_name <chr>,
#   infraspecific_name <chr>, isolate <chr>, version_status <chr>,
#   assembly_level <chr>, release_type <chr>, genome_rep <chr>,
#   seq_rel_date <date>, asm_name <chr>, submitter <chr>,
#   gbrs_paired_asm <chr>, paired_asm_comp <chr>, ftp_path <chr>,
#   excluded_from_refseq <chr>
```

#### A small negative example

In some cases your organism of interest will not be available in `NCBI RefSeq`.
Here an example:

```r
# check genome availability for Candida glabrata
is.genome.available(db = "refseq", organism = "Candida glabrata")
```

```r
Unfortunatey, no entry for 'Candida glabrata' was found in the 'refseq' database. 
Please consider specifying 'db = genbank' or 'db = ensembl' 
to check whether 'Candida glabrata' is availabe in these databases.
[1] FALSE
```

When now checking the availability in `NCBI Genbank` we find that indeed
'Candida glabrata' is available:

```r
# check genome availability for Candida glabrata
is.genome.available(db = "genbank", organism = "Candida glabrata")
```

```
Only a non-reference genome assembly is available for 'Candida glabrata'. 
Please make sure to specify the argument 'reference = FALSE' when running any get*() function.
[1] TRUE
```

Although, an entry is available the `is.genome.available()` warns us that
only a non-reference genome is available for 'Candida glabrata'.
To then retrieve the e.g. genome etc. files users need to specify the
`reference = FALSE` argument in the `get*()` functions.

For example:

```r
# retrieve non-reference genome 
getGenome(db = "genbank", organism = "Candida glabrata", reference = FALSE)
```

```
Unfortunately no genome file could be found for organism 'Candida glabrata'. Thus, the download of this organism has been omitted. Have you tried to specify 'reference = FALSE' ?
[1] "Not available"
```

### Example `NCBI Genbank` (?is.genome.available): 

Test whether or not the genome of `Homo sapiens` is present at NCBI Genbank.

```r
# checking whether or not the Homo sapiens 
# genome is avaialable for download
is.genome.available(db = "genbank", organism = "Homo sapiens")
```

```
A reference or representative genome assembly is available for 'Homo sapiens'.                           
More than one entry was found for 'Homo sapiens'. Please consider to run the function 'is.genome.available()' and specify 'is.genome.available(organism = Homo sapiens, db = genbank, details = TRUE)'. This will allow you to select the 'assembly_accession' identifier that can then be specified in all get*() functions.
[1] TRUE
```

Now with `details = TRUE`.

```r
# checking whether or not the Homo sapiens 
# genome is avaialable for download
is.genome.available(db = "genbank", organism = "Homo sapiens", details = TRUE)
```


```
 A tibble: 1,041 × 23                                                                                   
   assembly_accession bioproject biosample  wgs_master refseq_category
   <chr>              <chr>      <chr>      <chr>      <chr>          
 1 GCA_000001405.28   PRJNA31257 NA         NA         reference geno…
 2 GCA_000002115.2    PRJNA1431  SAMN02981… AADD00000… na             
 3 GCA_000002125.2    PRJNA19621 SAMN02981… ABBA00000… na             
 4 GCA_000002135.3    PRJNA10793 NA         NA         na             
 5 GCA_000004845.2    PRJNA42199 SAMN00003… ADDF00000… na             
 6 GCA_000005465.1    PRJNA42201 NA         DAAB00000… na             
 7 GCA_000181135.1    PRJNA28335 SAMN00001… ABKV00000… na             
 8 GCA_000185165.1    PRJNA59877 SAMN02981… AEKP00000… na             
 9 GCA_000212995.1    PRJNA19621 SAMN02981… ABSL00000… na             
10 GCA_000252825.1    PRJNA19621 SAMN02981… ABBA00000… na             
# … with 1,031 more rows, and 18 more variables: taxid <int>,
#   species_taxid <int>, organism_name <chr>,
#   infraspecific_name <chr>, isolate <chr>, version_status <chr>,
#   assembly_level <chr>, release_type <chr>, genome_rep <chr>,
#   seq_rel_date <date>, asm_name <chr>, submitter <chr>,
#   gbrs_paired_asm <chr>, paired_asm_comp <chr>, ftp_path <chr>,
#   excluded_from_refseq <chr>, X22 <chr>, X23 <chr>
```

As you can see there are several versions of the `Homo sapiens` genome
available for download from NCBI Genbank. Using the `assembly_accession`
id will now allow to specify which version shall be retrieved.


### Using `is.genome.available()` with ENSEMBL

Users can also specify `db = "ensembl"` to retrieve available organisms provided by ENSEMBL.
Again, users might experience a delay in the execution of this function when running it for the first time.
This is due to the download of ENSEMBL information which is then stored internally to enable a much faster
execution of this function in following runs. The corresponding information files are stored at `file.path(tempdir(), "ensembl_summary.txt")`.


### Example `ENSEMBL` (?is.genome.available): 

```r
# cheking whether Homo sapiens is available in the ENSEMBL database
is.genome.available(db = "ensembl", organism = "Homo sapiens")
```

```
A reference or representative genome assembly is available for 'Homo sapiens'.
[1] TRUE
```

```r
# retrieve details for Homo sapiens
is.genome.available(db = "ensembl", organism = "Homo sapiens", details = TRUE)
```

```
  division           assembly accession  release name  taxon_id strain
  <chr>              <chr>    <chr>        <int> <chr>    <int> <chr> 
1 EnsemblVertebrates GRCh38   GCA_00000…     104 homo…     9606 NA    
# … with 3 more variables: display_name <chr>, common_name <chr>,
#   strain_collection <chr>
```

Again, users can either specify the `taxid` or `accession id` when searching for organism entries.

```r
# retrieve details for Homo sapiens using taxid
is.genome.available(db = "ensembl", organism = "9606", details = TRUE)
```

```
division           assembly accession  release name  taxon_id strain
  <chr>              <chr>    <chr>        <int> <chr>    <int> <chr> 
1 EnsemblVertebrates GRCh38   GCA_00000…     104 homo…     9606 NA    
# … with 3 more variables: display_name <chr>, common_name <chr>,
#   strain_collection <chr>
```

```r
# retrieve details for Homo sapiens using accession id
is.genome.available(organism = "GCA_000001405.28", db = "ensembl", details = TRUE)
```

```
 division           assembly accession  release name  taxon_id strain
  <chr>              <chr>    <chr>        <int> <chr>    <int> <chr> 
1 EnsemblVertebrates GRCh38   GCA_00000…     104 homo…     9606 NA    
# … with 3 more variables: display_name <chr>, common_name <chr>,
#   strain_collection <chr>
```

__Please note that the `accession` id can change internally at ENSEMBL. E.g. in a recent case
the `accession` id changed from `GCA_000001405.25` to `GCA_000001405.27`. Hence, please
be careful and take this issue into account when you build automated retrieval scripts that are based
on `accession` ids.__

### Example `UniProt` (?is.genome.available):

Users can also check the availability of proteomes in the [UniProt database](http://www.uniprot.org)
by specifying:


```r
# retrieve information from UniProt
is.genome.available(db = "uniprot", "Homo sapiens", details = FALSE)
```

```
A reference or representative genome assembly is available for 'Homo sapiens'.
More than one entry was found for 'Homo sapiens'. Please consider to run the function 'is.genome.available()' and specify 'is.genome.available(organism = Homo sapiens, db = uniprot, details = TRUE)'. This will allow you to select the 'assembly_accession' identifier that can then be specified in all get*() functions.
[1] TRUE
```

or with details:

```r
# retrieve information from UniProt
is.genome.available(db = "uniprot", "Homo sapiens", details = TRUE)
```

```
 A tibble: 29 × 16
   name          description        isReferenceProte… isRepresentativ…
   <chr>         <chr>              <lgl>             <lgl>           
 1 Homo sapiens… Homo sapiens (Hom… TRUE              TRUE            
 2 Human associ… NA                 FALSE             TRUE            
 3 Human respir… NA                 FALSE             FALSE           
 4 Human respir… NA                 FALSE             FALSE           
 5 Human respir… NA                 FALSE             FALSE           
 6 Human respir… NA                 FALSE             FALSE           
 7 Human respir… NA                 FALSE             FALSE           
 8 Human respir… NA                 FALSE             FALSE           
 9 Human respir… NA                 FALSE             FALSE           
10 Human respir… NA                 FALSE             FALSE           
# … with 19 more rows, and 12 more variables:
#   genomeAssembly <df[,4]>, dbReference <list>, component <list>,
#   reference <list>, annotationScore <df[,1]>, scores <list>,
#   upid <chr>, modified <dbl>, taxonomy <int>, source <chr>,
#   superregnum <chr>, strain <chr>
```

Users can also search available species at UniProt via `taxid` or `upid` id.

Here `9606` defines the taxonomy id for `Homo sapiens`.

```r
# retrieve information from UniProt
is.genome.available(db = "uniprot", "9606", details = TRUE)
```

```
 A tibble: 3 × 15
  name  description isReferenceProt… isRepresentativ… genomeAssembly$…
  <chr> <chr>       <lgl>            <lgl>            <chr>           
1 Homo… Homo sapie… TRUE             TRUE             Ensembl         
2 Homo… NA          FALSE            FALSE            ENA/EMBL        
3 Homo… NA          FALSE            FALSE            ENA/EMBL        
# … with 10 more variables: dbReference <list>, component <list>,
#   reference <list>, annotationScore <df[,1]>, scores <list>,
#   upid <chr>, modified <dbl>, taxonomy <int>, source <chr>,
#   superregnum <chr>
```

Here `UP000005640` defines the `upid` for `Homo sapiens`.

```r
# retrieve information from UniProt
is.genome.available(db = "uniprot", "UP000005640", details = TRUE)
```

```
 name  description isReferenceProt… isRepresentativ… genomeAssembly$…
  <chr> <chr>       <lgl>            <lgl>            <chr>           
1 Homo… Homo sapie… TRUE             TRUE             Ensembl         
# … with 10 more variables: dbReference <list>, component <list>,
#   reference <list>, annotationScore <df[,1]>, scores <list>,
#   upid <chr>, modified <dbl>, taxonomy <int>, source <chr>,
#   superregnum <chr>
```



In general, the argument `db` specifies from which database (`refseq`, `genbank`, `ensembl` or `uniprot`) organism
information shall be retrieved. Options are:

- `db = 'refseq'`
- `db = 'genbank'`
- `db = 'ensembl'`
- `db = 'uniprot'`


### Listing the total number of available genomes 

In some cases it might be useful to check how many genomes (in total) are available in the different databases.

Users can determine the total number of available genomes using the `listGenomes()` function.

Example `refseq`: 

```{r,eval=FALSE}
length(listGenomes(db = "refseq"))
```

```
[1] 24910
```

Example `genbank`: 

```{r,eval=FALSE}
length(listGenomes(db = "genbank"))
```

```
[1] 35298
```

Example `ensembl`: 

```{r,eval=FALSE}
length(listGenomes(db = "ensembl"))
```

```
[1] 310
```


Hence, currently 24910 genomes (including all kingdoms of life) are stored on `NCBI RefSeq` (as of 16/11/2021).

### Retrieving kingdom, group and subgroup information

Using this example users can retrieve the number of available species for each kingdom of life:

Example `refseq`:

```{r,eval=FALSE}
# the number of genomes available for each kingdom
listKingdoms(db = "refseq")
```

```
    Archaea  Bacteria Eukaryota   Viruses 
      330     11797      1002     10637 
```

Example `genbank`:

```{r,eval=FALSE}
# the number of genomes available for each kingdom
listKingdoms(db = "genbank")
```

```
  Archaea  Bacteria Eukaryota 
     1651     25904      7742
```

Example `ENSEMBL`:

```{r,eval=FALSE}
# the number of genomes available for each kingdom
listKingdoms(db = "ensembl")
```

```
EnsemblVertebrates 
               310
```

### Analogous computations can be performed for `groups` and `subgroups`

Unfortunately, `ENSEMBL` does not provide group or subgroup information.
Therefore, group and subgroup listings are limited to `refseq` and `genbank`.


Example `refseq`:

```{r,eval=FALSE}
# the number of genomes available for each group
listGroups(db = "refseq")
```

```
 Abditibacteriota 
                                               1 
                               Acidithiobacillia 
                                               8 
                                  Acidobacteriia 
                                              22 
                                Ackermannviridae 
                                              65 
                                  Actinobacteria 
                                            2640 
                                    Adenoviridae 
                                              64 
                                    Adomaviridae 
                                               2 
                                    Aliusviridae 
                                               2 
                               Alloherpesviridae 
                                              13 
                               Alphaflexiviridae 
                                              62 
                             Alphaproteobacteria 
                                            1762 
  ...
```

Example `genbank`:

```{r,eval=FALSE}
# the number of genomes available for each group
listGroups(db = "genbank")
```

```
                       Abditibacteriota                       Acidithiobacillia 
                                      1                                       8 
                         Acidobacteriia                          Actinobacteria 
                                     24                                    1596 
                    Alphaproteobacteria                              Amphibians 
                                   1604                                       6 
                          Apicomplexans                               Aquificae 
                                     47                                       7 
                           Archaeoglobi                         Armatimonadetes 
                                      5                                      38 
                            Ascomycetes                Bacteria candidate phyla 
                                    689                                    3532 
           Bacteroidetes/Chlorobi group                               Balneolia 
                                   2247                                      44 
                         Basidiomycetes                      Betaproteobacteria 
                                    204                                     751 
                                  Birds                          Blastocatellia 
                                     80                                       2 
                           Caldisericia                  candidate division WS1 
                                      1                                       1 
        candidate division Zixibacteria               Candidatus Aegiribacteria 
                                     17                                       1 
             Candidatus Aenigmarchaeota               Candidatus Bathyarchaeota 
                                     14                                      42 
               Candidatus Cloacimonetes               Candidatus Diapherotrites 
                                     88                                      11 
            Candidatus Fermentibacteria            Candidatus Geothermarchaeota 
                                      8                                       3 
           Candidatus Heimdallarchaeota              Candidatus Hydrogenedentes 
                                      4                                      10 
                Candidatus Korarchaeota                    Candidatus Kryptonia 
                                      1                                       4 
        Candidatus Lambdaproteobacteria              Candidatus Latescibacteria 
                                      6                                      13 
               Candidatus Lokiarchaeota               Candidatus Marinimicrobia 
                                      2                                      92 
               Candidatus Marsarchaeota                Candidatus Micrarchaeota 
                                     15                                      35 
                Candidatus Moduliflexus             Candidatus Muproteobacteria 
                                      1                                      14 
           Candidatus Nanohaloarchaeota                Candidatus Odinarchaeota 
                                     16                                       1 
                Candidatus Omnitrophica                Candidatus Pacearchaeota 
                                    126                                      41 
               Candidatus Parvarchaeota                Candidatus Tectomicrobia 
                                     11                                       6 
               Candidatus Thorarchaeota                 Candidatus Vecturithrix 
                                      7                                       1 
         Candidatus Verstraetearchaeota               Candidatus Woesearchaeota 
                                      5                                      36 
                             Chlamydiae                             Chloroflexi 
                                     43                                     403 
                         Chrysiogenetes                    Coprothermobacterota 
                                      1                                       1 
                          Crenarchaeota     Cyanobacteria/Melainabacteria group 
                                     54                                     172 
                        Deferribacteres                     Deinococcus-Thermus 
                                      8                                      27 
             delta/epsilon subdivisions                            Dictyoglomia 
                                    615                                       1 
                          Elusimicrobia                            Endomicrobia 
                                      1                                       2 
                  environmental samples                           Fibrobacteres 
                                      3                                      15 
                             Firmicutes                                  Fishes 
                                   2663                                     168 
                              Flatworms                           Fusobacteriia 
                                     35                                      29 
                    Gammaproteobacteria                        Gemmatimonadetes 
                                   2041                                      82 
                            Green Algae                            Hadesarchaea 
                                     30                                       5 
                           Halobacteria                              Holophagae 
                                    162                                       9 
                      Hydrogenophilalia                                 Insects 
                                     32                                     280 
                           Kinetoplasts                      Kiritimatiellaeota 
                                     33                                       3 
                            Land Plants                           Lentisphaerae 
                                    282                                      32 
                                Mammals                         Methanobacteria 
                                    145                                      31 
                           Methanococci                         Methanomicrobia 
                                      2                                      77 
                  Methanonatronarchaeia                             Methanopyri 
                                      2                                       1 
                          Nanoarchaeota nitrifying bacterium enrichment culture 
                                     11                                       1 
                            Nitrospinae                              Nitrospira 
                                     19                                      31 
                            Oligoflexia                                   Other 
                                     87                                      15 
                          Other Animals                             Other Fungi 
                                    125                                      58 
                           Other Plants                          Other Protists 
                                      1                                     128 
                         Planctomycetes                                Reptiles 
                                    171                                      21 
                           Rhodothermia                              Roundworms 
                                      4                                      91 
                           Solibacteres                            Spirochaetia 
                                      7                                     115 
                            Synergistia                             Tenericutes 
                                     43                                     141 
                         Thaumarchaeota                           Theionarchaea 
                                     92                                       2 
                            Thermococci                   Thermodesulfobacteria 
                                     21                                       7 
                         Thermoplasmata                             Thermotogae 
                                    135                                      23 
             unclassified Acidobacteria    unclassified Archaea (miscellaneous) 
                                     94                                      52 
  unclassified Bacteria (miscellaneous)                unclassified Caldiserica 
                                    197                                      11 
           unclassified Calditrichaeota              unclassified Elusimicrobia 
                                      2                                     101 
             unclassified Euryarchaeota                unclassified Nitrospirae 
                                    348                                      90 
            unclassified Proteobacteria            unclassified Rhodothermaeota 
                                    130                                       4 
              unclassified Spirochaetes              unclassified Synergistetes 
                                     37                                       3 
               unclassified Thermotogae                     uncultured archaeon 
                                      1                                       1 
            uncultured archaeon A07HB70             uncultured archaeon A07HN63 
                                      1                                       1 
            uncultured archaeon A07HR60             uncultured archaeon A07HR67 
                                      1                                       1 
                   uncultured bacterium                         Verrucomicrobia 
                                      1                                     296 
                     Zetaproteobacteria 
                                     41 
```

> **__Note__** that when running the `listGenomes()` function for the first time, it might take a while until the function returns any results, because necessary information need to be downloaded from NCBI and ENSEMBL databases. All subsequent executions of `listGenomes()` will then respond very fast, because they will access the corresponding files stored on your hard drive.

## Downloading Biological Sequences and Annotations

After checking for the availability of sequence information for an organism of interest,
the next step is to download the corresponding genome, proteome, CDS, or GFF file.
The following functions allow users to download proteomes, genomes, CDS and GFF files from several
database resources such as: [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/), [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/), [ENSEMBL](http://www.ensembl.org/index.html). When a corresponding proteome, genome, CDS or GFF file was
loaded to your hard-drive, a documentation `*.txt` file is generated storing `File Name`, `Organism`,
`Database`, `URL`, `DATE`, `assembly_accession`, `bioproject`, `biosample`, `taxid`, `version_status`, `release_type`,
`seq_rel_date` etc. information of the retrieved file. This way a better reproducibility of proteome, genome, CDS and GFF versions
used for subsequent data analyses can be achieved.


### Genome Retrieval

The easiest way to download a genome is to use the `getGenome()` function.

In this example we will download the genome of `Homo sapiens`.

The `getGenome()` function is an interface function to the [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/), [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/),
[ENSEMBL](http://www.ensembl.org/index.html) databases from
which corresponding genomes can be retrieved.

The `db` argument specifies from which database genome assemblies in `*.fasta` file format shall be retrieved.

Options are:

- `db = "refseq"` for retrieval from [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/)
- `db = "genbank"` for retrieval from [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/)
- `db = "ensembl"` for retrieval from [ENSEMBL](http://www.ensembl.org/index.html)


Furthermore, users need to specify the `scientific name`, the `taxid` (= [NCBI Taxnonomy](https://www.ncbi.nlm.nih.gov/taxonomy) identifier), or the `accession identifier` of the organism of interest for which a genome assembly shall be downloaded, e.g. `organism = "Homo sapiens"` or `organism = "9606"` or `organism = "GCF_000001405.37"`.
Finally, the `path` argument specifies the folder path in which the corresponding assembly shall be locally stored. In case users would like to store the genome file at a different location,
they can specify the `path = file.path("put","your","path","here")` argument (e.g. `file.path("_ncbi_downloads","genomes")`).

### Example NCBI RefSeq:

```{r,eval=FALSE}
# download the genome of Homo sapiens from refseq
# and store the corresponding genome file in '_ncbi_downloads/genomes'
HS.genome.refseq <- getGenome( db       = "refseq", 
                               organism = "Homo sapiens",
                               path     = file.path("_ncbi_downloads","genomes") )
```

In this example, `getGenome()` creates a directory named `'_ncbi_downloads/genomes'` into which the corresponding genome named `GCF_000001405.34_GRCh38.p8_genomic.fna.gz` is downloaded. The return value of `getGenome()` is the folder path to the downloaded genome file
that can then be used as input to the `read_genome()` function. The variable `HS.genome.refseq` stores the path to the downloaded genome. 

Users can also omit the `path` argument if they wish to store the genome in their current working directory. E.g.:

```{r,eval=FALSE}
# download the genome of Homo sapiens from refseq
# and store the corresponding genome file in '_ncbi_downloads/genomes'
HS.genome.refseq <- getGenome( db       = "refseq", 
                               organism = "Homo sapiens")
```


Subsequently, users can use the `read_genome()` function to import the genome into the R session. Users can choose to work with the genome sequence in R either as [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) object (`obj.type = "Biostrings"`) or [data.table](https://github.com/Rdatatable/data.table/wiki) object
(`obj.type = "data.table"`) by specifying the `obj.type` argument of the `read_genome()` function. 

```{r,eval=FALSE}
# import downloaded genome as Biostrings object
Human_Genome <- read_genome(file = HS.genome.refseq)

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_Genome
```

```
  A DNAStringSet instance of length 551
          width seq                                                names               
  [1] 248956422 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN NC_000001.11 Homo...
  [2]    175055 GAATTCAGCTGAGAAGAACAGGCA...TGTTTGTCAGTCACATAGAATTC NT_187361.1 Homo ...
  [3]     32032 AGGGGTCTGCTTAGAGAGGGTCTC...TGACTTACGTTGACGTGGAATTC NT_187362.1 Homo ...
  [4]    127682 GATCGAGACTATCCTGGCTAACAC...ATTGTCAATTGGGACCTTTGATC NT_187363.1 Homo ...
  [5]     66860 GAATTCATTCGATGACGATTCCAT...AAAAAACTCTCAGCCACGAATTC NT_187364.1 Homo ...
  ...       ... ...
[547]    170148 TTTCTTTCTTTTTTTTTTTTTTGT...GTCACAGGACTCATGGGGAATTC NT_187685.1 Homo ...
[548]    215732 TGTGGTGAGGACCCTTAAGATCTA...GTCACAGGACTCATGGGGAATTC NT_187686.1 Homo ...
[549]    170537 TCTACTCTCCCATGCTTGCCTCGG...GTCACAGGACTCATGGGGAATTC NT_187687.1 Homo ...
[550]    177381 GATCTATCTGTATCTCCACAGGTG...GTCACAGGACTCATGGGGAATTC NT_113949.2 Homo ...
[551]     16569 GATCACAGGTCTATCACCCTATTA...CCCTTAAATAAGACATCACGATG NC_012920.1 Homo ...
```

Internally, a text file named `doc_Homo_sapiens_db_refseq.txt` is generated. The information stored in this log file is structured as follows:

```
File Name: Homo_sapiens_genomic_refseq.fna.gz
Organism Name: Homo_sapiens
Database: NCBI refseq
URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/
GCF_000001405.35_GRCh38.p9/GCF_000001405.35_GRCh38.p9_genomic.fna.gz
Download_Date: Sat Oct 22 12:41:07 2016
refseq_category: reference genome
assembly_accession: GCF_000001405.35
bioproject: PRJNA168
biosample: NA
taxid: 9606
infraspecific_name: NA
version_status: latest
release_type: Patch
genome_rep: Full
seq_rel_date: 2016-09-26
submitter: Genome Reference Consortium
```

In addition, the genome summary statistics for the retrieved species is stored locally (`doc_Homo_sapiens_db_refseq_summary_statistics.tsv`) to provide users with insights regarding the genome assembly quality (see `?summary_genome()` for details).
This file can be used as `Supplementary Information` file in publications to facilitate reproducible research.
Most comparative genomics studies do not consider differences in genome assembly qualities when comparing the genomes of diverse species.
This way, they expose themselves to technical artifacts that might generate patterns mistaken to be of biological relevance whereas
in reality they just reflect the difference in genome assembly quality. Considering the quality of genome assemblies when downloading
the genomic sequences will help researchers to avoid these pitfalls. 

The summary statistics include:

- `genome_size_mbp`: Genome size in mega base pairs

- `n50_mbp`: The N50 contig size of the genome assembly in mega base pairs

- `n_seqs`: The number of chromosomes/scaffolds/contigs of the genome assembly file

- `n_nnn`: The absolute number of NNNs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_nnn`: The percentage (relative frequency) of NNNs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file

- `genome_entropy`: The Shannon Entropy of the genome assembly file (median entropy over all individual chromosome entropies)

- `n_gc`: The total number of GCs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_gc`: The (relative frequency) of GCs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file




In summary, the `getGenome()` and `read_genome()` functions allow users to retrieve genome assemblies by specifying
the scientific name of the organism of interest and allow them to import the retrieved genome assembly e.g. as `Biostrings` object.
Thus, users can then perform the `Biostrings notation` to work with downloaded genomes and can rely on the log
file generated by `getGenome()` to better document the source and version of genome assemblies used for subsequent studies.

Alternatively, users can perform the pipeline logic of the [magrittr](https://github.com/smbache/magrittr) package:

```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Human_Genome <- getGenome( db       = "refseq", 
                           organism = "Homo sapiens",
                           path     = file.path("_ncbi_downloads","genomes")) %>%
    read_genome()
```

```{r,eval=FALSE}
Human_Genome
```

```
  A DNAStringSet instance of length 551
          width seq                                                names               
  [1] 248956422 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN NC_000001.11 Homo...
  [2]    175055 GAATTCAGCTGAGAAGAACAGGCA...TGTTTGTCAGTCACATAGAATTC NT_187361.1 Homo ...
  [3]     32032 AGGGGTCTGCTTAGAGAGGGTCTC...TGACTTACGTTGACGTGGAATTC NT_187362.1 Homo ...
  [4]    127682 GATCGAGACTATCCTGGCTAACAC...ATTGTCAATTGGGACCTTTGATC NT_187363.1 Homo ...
  [5]     66860 GAATTCATTCGATGACGATTCCAT...AAAAAACTCTCAGCCACGAATTC NT_187364.1 Homo ...
  ...       ... ...
[547]    170148 TTTCTTTCTTTTTTTTTTTTTTGT...GTCACAGGACTCATGGGGAATTC NT_187685.1 Homo ...
[548]    215732 TGTGGTGAGGACCCTTAAGATCTA...GTCACAGGACTCATGGGGAATTC NT_187686.1 Homo ...
[549]    170537 TCTACTCTCCCATGCTTGCCTCGG...GTCACAGGACTCATGGGGAATTC NT_187687.1 Homo ...
[550]    177381 GATCTATCTGTATCTCCACAGGTG...GTCACAGGACTCATGGGGAATTC NT_113949.2 Homo ...
[551]     16569 GATCACAGGTCTATCACCCTATTA...CCCTTAAATAAGACATCACGATG NC_012920.1 Homo ...
```

#### Use `taxid` id for genome retrieval

Alternatively, instead of specifying the scientific name in the argument `organism`
users can specify the `taxonomy` id of the corresponding organism. Here, we specify
the taxonomy id `559292` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "refseq",
            organism = "559292") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```


```
A DNAStringSet instance of length 17
       width seq                                      names               
 [1]  230218 CCACACCACACCCACACAC...GGGTGTGGTGTGTGTGGG NC_001133.9 Sacch...
 [2]  813184 AAATAGCCCTCATGTACGT...TGTGGTGTGTGGGTGTGT NC_001134.8 Sacch...
 [3]  316620 CCCACACACCACACCCACA...GTGGGTGTGGTGTGTGTG NC_001135.5 Sacch...
 [4] 1531933 ACACCACACCCACACCACA...GTAGTAAGTAGCTTTTGG NC_001136.10 Sacc...
 [5]  576874 CGTCTCCTCCAAGCCCTGT...ATTTTCATTTTTTTTTTT NC_001137.3 Sacch...
 ...     ... ...
[13]  924431 CCACACACACACCACACCC...TGTGGTGTGTGTGTGGGG NC_001145.3 Sacch...
[14]  784333 CCGGCTTTCTGACCGAAAT...TGTGGGTGTGGTGTGGGT NC_001146.8 Sacch...
[15] 1091291 ACACCACACCCACACCACA...TGTGTGGGTGTGGTGTGT NC_001147.6 Sacch...
[16]  948066 AAATAGCCCTCATGTACGT...TTTAATTTCGGTCAGAAA NC_001148.4 Sacch...
[17]   85779 TTCATAATTAATTTTTTAT...TATAATATAATATCCATA NC_001224.1 Sacch...
```

#### Use `assembly_accession` id for genome retrieval

Alternatively, instead of specifying the scientific name or taxonomy in the argument `organism`
users can specify the `assembly_accession` id of the corresponding organism. Here, we specify
the `assembly_accession` id `GCF_000146045.2` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "refseq",
            organism = "GCF_000146045.2") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```

```
 A DNAStringSet instance of length 17
       width seq                                                    names               
 [1]  230218 CCACACCACACCCACACACCCACACA...GTGGTGTGGGTGTGGTGTGTGTGGG NC_001133.9 Sacch...
 [2]  813184 AAATAGCCCTCATGTACGTCTCCTCC...GTGTGGGTGTGGTGTGTGGGTGTGT NC_001134.8 Sacch...
 [3]  316620 CCCACACACCACACCCACACCACACC...GGGTGTGGTGGGTGTGGTGTGTGTG NC_001135.5 Sacch...
 [4] 1531933 ACACCACACCCACACCACACCCACAC...AATAAAGGTAGTAAGTAGCTTTTGG NC_001136.10 Sacc...
 [5]  576874 CGTCTCCTCCAAGCCCTGTTGTCTCT...GGGTTTCATTTTCATTTTTTTTTTT NC_001137.3 Sacch...
 ...     ... ...
[13]  924431 CCACACACACACCACACCCACACCAC...GTGTGGGTGTGGTGTGTGTGTGGGG NC_001145.3 Sacch...
[14]  784333 CCGGCTTTCTGACCGAAATTAAAAAA...GGGTGTGTGTGGGTGTGGTGTGGGT NC_001146.8 Sacch...
[15] 1091291 ACACCACACCCACACCACACCCACAC...GAGAGAGTGTGTGGGTGTGGTGTGT NC_001147.6 Sacch...
[16]  948066 AAATAGCCCTCATGTACGTCTCCTCC...TTTTTTTTTTAATTTCGGTCAGAAA NC_001148.4 Sacch...
[17]   85779 TTCATAATTAATTTTTTATATATATA...GCTTAATTATAATATAATATCCATA NC_001224.1 Sacch...
```

### Example NCBI Genbank:

Genome retrieval from `NCBI Genbank`.

```{r,eval=FALSE}
# download the genome of Homo sapiens from Genbank
# and store the corresponding genome file in '_ncbi_downloads/genomes'
HS.genome.genbank <- getGenome( db       = "genbank", 
                                organism = "Homo sapiens",
                                path     = file.path("_ncbi_downloads","genomes") )
```

```{r,eval=FALSE}
# import downloaded genome as Biostrings object
Human_Genome <- read_genome(file = HS.genome.genbank)
```

```{r,eval=FALSE}
# look at the Biostrings object
Human_Genome
```

```
  A DNAStringSet instance of length 551
          width seq                                                names               
  [1] 248956422 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN CM000663.2 Homo s...
  [2]    175055 GAATTCAGCTGAGAAGAACAGGCA...TGTTTGTCAGTCACATAGAATTC KI270706.1 Homo s...
  [3]     32032 AGGGGTCTGCTTAGAGAGGGTCTC...TGACTTACGTTGACGTGGAATTC KI270707.1 Homo s...
  [4]    127682 GATCGAGACTATCCTGGCTAACAC...ATTGTCAATTGGGACCTTTGATC KI270708.1 Homo s...
  [5]     66860 GAATTCATTCGATGACGATTCCAT...AAAAAACTCTCAGCCACGAATTC KI270709.1 Homo s...
  ...       ... ...
[547]    170148 TTTCTTTCTTTTTTTTTTTTTTGT...GTCACAGGACTCATGGGGAATTC KI270931.1 Homo s...
[548]    215732 TGTGGTGAGGACCCTTAAGATCTA...GTCACAGGACTCATGGGGAATTC KI270932.1 Homo s...
[549]    170537 TCTACTCTCCCATGCTTGCCTCGG...GTCACAGGACTCATGGGGAATTC KI270933.1 Homo s...
[550]    177381 GATCTATCTGTATCTCCACAGGTG...GTCACAGGACTCATGGGGAATTC GL000209.2 Homo s...
[551]     16569 GATCACAGGTCTATCACCCTATTA...CCCTTAAATAAGACATCACGATG J01415.2 Homo sap...
```


#### Use `taxonomy` id for genome retrieval

Alternatively, instead of specifying the scientific name in the argument `organism`
users can specify the `taxonomy` id of the corresponding organism. Here, we specify
the taxonomy id `559292` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "genbank",
            organism = "559292") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```

```
A DNAStringSet instance of length 16
       width seq                                names               
 [1]  230218 CCACACCACACCCACA...TGTGGTGTGTGTGGG BK006935.2 TPA_in...
 [2]  813184 AAATAGCCCTCATGTA...GGTGTGTGGGTGTGT BK006936.2 TPA_in...
 [3]  316620 CCCACACACCACACCC...GGTGTGGTGTGTGTG BK006937.2 TPA_in...
 [4] 1531933 ACACCACACCCACACC...GTAAGTAGCTTTTGG BK006938.2 TPA_in...
 [5]  576874 CGTCTCCTCCAAGCCC...TTCATTTTTTTTTTT BK006939.2 TPA_in...
 ...     ... ...
[12] 1078177 CACACACACACACCAC...ACATGAGGGCTATTT BK006945.2 TPA_in...
[13]  924431 CCACACACACACCACA...GGTGTGTGTGTGGGG BK006946.2 TPA_in...
[14]  784333 CCGGCTTTCTGACCGA...GGGTGTGGTGTGGGT BK006947.3 TPA_in...
[15] 1091291 ACACCACACCCACACC...GTGGGTGTGGTGTGT BK006948.2 TPA_in...
[16]  948066 AAATAGCCCTCATGTA...AATTTCGGTCAGAAA BK006949.2 TPA_in...
```

#### Use `assembly_accession` id for genome retrieval


Alternatively, instead of specifying the scientific name or taxonomy in the argument `organism`
users can specify the `assembly_accession` id of the corresponding organism. Here, we specify
the `assembly_accession` id `GCA_000146045.2` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "genbank",
            organism = "GCA_000146045.2") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```

```
A DNAStringSet instance of length 16
       width seq                                  names               
 [1]  230218 CCACACCACACCCACAC...GTGTGGTGTGTGTGGG BK006935.2 TPA_in...
 [2]  813184 AAATAGCCCTCATGTAC...TGGTGTGTGGGTGTGT BK006936.2 TPA_in...
 [3]  316620 CCCACACACCACACCCA...GGGTGTGGTGTGTGTG BK006937.2 TPA_in...
 [4] 1531933 ACACCACACCCACACCA...AGTAAGTAGCTTTTGG BK006938.2 TPA_in...
 [5]  576874 CGTCTCCTCCAAGCCCT...TTTCATTTTTTTTTTT BK006939.2 TPA_in...
 ...     ... ...
[12] 1078177 CACACACACACACCACC...TACATGAGGGCTATTT BK006945.2 TPA_in...
[13]  924431 CCACACACACACCACAC...TGGTGTGTGTGTGGGG BK006946.2 TPA_in...
[14]  784333 CCGGCTTTCTGACCGAA...TGGGTGTGGTGTGGGT BK006947.3 TPA_in...
[15] 1091291 ACACCACACCCACACCA...TGTGGGTGTGGTGTGT BK006948.2 TPA_in...
[16]  948066 AAATAGCCCTCATGTAC...TAATTTCGGTCAGAAA BK006949.2 TPA_in...
```

In addition, the genome summary statistics for the retrieved species is stored locally (`doc_saccharomyces_cerevisiae_db_genbank_summary_statistics.tsv`) to provide users with insights regarding the genome assembly quality (see `?summary_genome()` for details).
This file can be used as `Supplementary Information` file in publications to facilitate reproducible research.
Most comparative genomics studies do not consider differences in genome assembly qualities when comparing the genomes of diverse species.
This way, they expose themselves to technical artifacts that might generate patterns mistaken to be of biological relevance whereas
in reality they just reflect the difference in genome assembly quality. Considering the quality of genome assemblies when downloading
the genomic sequences will help researchers to avoid these pitfalls. 

The summary statistics include:

- `genome_size_mbp`: Genome size in mega base pairs

- `n50_mbp`: The N50 contig size of the genome assembly in mega base pairs

- `n_seqs`: The number of chromosomes/scaffolds/contigs of the genome assembly file

- `n_nnn`: The absolute number of NNNs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_nnn`: The percentage (relative frequency) of NNNs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file

- `genome_entropy`: The Shannon Entropy of the genome assembly file (median entropy over all individual chromosome entropies)

- `n_gc`: The total number of GCs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_gc`: The (relative frequency) of GCs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file

### Example ENSEMBL:

```{r,eval=FALSE}
# download the genome of Homo sapiens from ENSEMBL
# and store the corresponding genome file in '_ncbi_downloads/genomes'
HS.genome.ensembl <- getGenome( db       = "ensembl", 
                                organism = "Homo sapiens",
                                path     = file.path("_ncbi_downloads","genomes") ,
                                assembly_type = "primary_assembly")
```

```{r,eval=FALSE}
# import downloaded genome as Biostrings object
Human_Genome <- read_genome(file = HS.genome.ensembl)
```

```{r,eval=FALSE}
# look at the Biostrings object
Human_Genome
```

```
  A DNAStringSet instance of length 524
          width seq                                                names               
  [1] 248956422 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN 1 dna:chromosome ...
  [2] 242193529 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN 2 dna:chromosome ...
  [3] 198295559 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN 3 dna:chromosome ...
  [4] 190214555 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN 4 dna:chromosome ...
  [5] 181538259 NNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNN 5 dna:chromosome ...
  ...       ... ...
[520]       993 GCCCCACGTCCGGGAGGGAGGTGG...GAGGGAGGTGGGGGGTCAGCCCT KI270539.1 dna:sc...
[521]       990 TTTCATAGAGCATGTTTGAAACAC...TCAGAAACTTGTTGTGATGTGTG KI270385.1 dna:sc...
[522]       981 AGATTTCGTTGGAACGGGATAAAC...TCAGCTTTCAAACACTCTTTTTG KI270423.1 dna:sc...
[523]       971 ATTTGCGATGTGTGTTCTCAACTA...TTGGATAGCTTTGAAGTTTTCGT KI270392.1 dna:sc...
[524]       970 AAGTGGATATTTGGATAGCTTTGA...TCCTCAATAACAGAGTTGAACCT KI270394.1 dna:sc...
```
If you are using `db = "ensembl"` you can set `assembly_type = "primary_assembly"`. For the human genome the toplevel genome assembly size is > 70 GB uncompressed, while primary assembly is just a few GB. 
Users can also specify the `release` argument which denotes the database release version of ENSEMBL (`db = "ensembl"`) in case they wish to download archived vgenome versions. Default is release = NULL meaning that the most recent database version is used.

In addition, the genome summary statistics for the retrieved species is stored locally (`doc_homo_sapiens_db_ensembl_summary_statistics.tsv`) to provide users with insights regarding the genome assembly quality (see `?summary_genome()` for details).
This file can be used as `Supplementary Information` file in publications to facilitate reproducible research.
Most comparative genomics studies do not consider differences in genome assembly qualities when comparing the genomes of diverse species.
This way, they expose themselves to technical artifacts that might generate patterns mistaken to be of biological relevance whereas
in reality they just reflect the difference in genome assembly quality. Considering the quality of genome assemblies when downloading
the genomic sequences will help researchers to avoid these pitfalls. 

The summary statistics include:

- `genome_size_mbp`: Genome size in mega base pairs

- `n50_mbp`: The N50 contig size of the genome assembly in mega base pairs

- `n_seqs`: The number of chromosomes/scaffolds/contigs of the genome assembly file

- `n_nnn`: The absolute number of NNNs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_nnn`: The percentage (relative frequency) of NNNs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file

- `genome_entropy`: The Shannon Entropy of the genome assembly file (median entropy over all individual chromosome entropies)

- `n_gc`: The total number of GCs (over all chromosomes or scaffolds or contigs) in the genome assembly file

- `rel_gc`: The (relative frequency) of GCs (over all chromosomes or scaffolds or contigs) compared to the total number of nucleotides in the genome assembly file

#### Use `taxonomy` id for genome retrieval

Alternatively, instead of specifying the scientific name in the argument `organism`
users can specify the `taxonomy` id of the corresponding organism. Here, we specify
the taxonomy id `4932` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "ensembl",
            organism = "4932") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```

```
  A DNAStringSet instance of length 17
       width seq                            names               
 [1]  230218 CCACACCACACCCA...TGGTGTGTGTGGG I dna:chromosome ...
 [2]  813184 AAATAGCCCTCATG...TGTGTGGGTGTGT II dna:chromosome...
 [3]  316620 CCCACACACCACAC...TGTGGTGTGTGTG III dna:chromosom...
 [4] 1531933 ACACCACACCCACA...AAGTAGCTTTTGG IV dna:chromosome...
 [5]  576874 CGTCTCCTCCAAGC...CATTTTTTTTTTT V dna:chromosome ...
 ...     ... ...
[13]  924431 CCACACACACACCA...TGTGTGTGTGGGG XIII dna:chromoso...
[14]  784333 CCGGCTTTCTGACC...GTGTGGTGTGGGT XIV dna:chromosom...
[15] 1091291 ACACCACACCCACA...GGGTGTGGTGTGT XV dna:chromosome...
[16]  948066 AAATAGCCCTCATG...TTTCGGTCAGAAA XVI dna:chromosom...
[17]   85779 TTCATAATTAATTT...TATAATATCCATA Mito dna:chromoso...
```

#### Use `assembly_accession` id for genome retrieval


Alternatively, instead of specifying the scientific name or taxonomy in the argument `organism`
users can specify the `assembly_accession` id of the corresponding organism. Here, we specify
the `assembly_accession` id `GCA_000146045.2` which encodes the species `Saccharomyces cerevisiae`.


```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import genome as Biostrings object
Scerevisiae_Genome <- getGenome(
            db       = "ensembl",
            organism = "GCA_000146045.2") %>%
    read_genome()
```

```{r,eval=FALSE}
Scerevisiae_Genome
```

```
  A DNAStringSet instance of length 17
       width seq                                names               
 [1]  230218 CCACACCACACCCACA...TGTGGTGTGTGTGGG I dna:chromosome ...
 [2]  813184 AAATAGCCCTCATGTA...GGTGTGTGGGTGTGT II dna:chromosome...
 [3]  316620 CCCACACACCACACCC...GGTGTGGTGTGTGTG III dna:chromosom...
 [4] 1531933 ACACCACACCCACACC...GTAAGTAGCTTTTGG IV dna:chromosome...
 [5]  439888 CACACACACCACACCC...GTGTGGTGTGTGTGT IX dna:chromosome...
 ...     ... ...
[13] 1078177 CACACACACACACCAC...ACATGAGGGCTATTT XII dna:chromosom...
[14]  924431 CCACACACACACCACA...GGTGTGTGTGTGGGG XIII dna:chromoso...
[15]  784333 CCGGCTTTCTGACCGA...GGGTGTGGTGTGGGT XIV dna:chromosom...
[16] 1091291 ACACCACACCCACACC...GTGGGTGTGGTGTGT XV dna:chromosome...
[17]  948066 AAATAGCCCTCATGTA...AATTTCGGTCAGAAA XVI dna:chromosom...
```

### GenomeSet Retrieval

The `getGenomeSet()` function enables users to retrieve genome files for
multiple species. This is convenient when users wish to bulk download a particular
set of species. Internally, a folder named `set_genomes` is generated in which
genomes and genome info files are stored. In addition, users can specify the arguments: `clean_retrieval` and `gunzip` (both are `TRUE` by default) to clean file names and automatically unzip downloaded files. 

Example:

```r
# specify the species names
download_species <- c("Arabidopsis thaliana", 
                      "Arabidopsis lyrata", 
                      "Capsella rubella")
# retrieve these three species from NCBI RefSeq                       
biomartr::getGenomeSet("refseq", organisms = download_species, path = "set_genomes")
```

If the download process was interrupted, users can re-run the function
and it will only download missing genomes. In cases users wish
to download everything again and updating existing genomes, they may
specify the argument `update = TRUE`.

### Proteome Retrieval

The `getProteome()` function is an interface function to the [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/), [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/),
[ENSEMBL](http://www.ensembl.org/index.html), and [UniProt](http://uniprot.org/) databases from
which corresponding proteomes can be retrieved. It works analogous to `getGenome()`. 

The `db` argument specifies from which database proteomes in `*.fasta` file format shall be retrieved.

Options are:

- `db = "refseq"` for retrieval from [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/)
- `db = "genbank"` for retrieval from [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/)
- `db = "ensembl"` for retrieval from [ENSEMBL](http://www.ensembl.org/index.html)
- `db = "uniprot" ` for retrieval from [UniProt](http://uniprot.org/)


Furthermore, again users need to specify the scientific name of the organism of interest for which a proteomes shall be downloaded, e.g. `organism = "Homo sapiens"`.
Finally, the `path` argument specifies the folder path in which the corresponding proteome shall be locally stored. In case users would like to store the proteome file at a different location,
they can specify the `path = file.path("put","your","path","here")` argument.

### Example `NCBI RefSeq`:

```{r,eval=FALSE}
# download the proteome of Homo sapiens from refseq
# and store the corresponding proteome file in '_ncbi_downloads/proteomes'
HS.proteome.refseq <- getProteome( db       = "refseq", 
                                   organism = "Homo sapiens",
                                   path     = file.path("_ncbi_downloads","proteomes"))
```

In this example, `getProteome()` creates a directory named `'_ncbi_downloads/proteomes'` into which the corresponding genome named `GCF_000001405.34_GRCh38.p8_protein.faa.gz` is downloaded. The return value of `getProteome()` is the folder path to the downloaded proteome file
that can then be used as input to the `read_proteome()` function. The variable `HS.proteome.refseq` stores the path to the downloaded proteome. 
Subsequently, users can use the `read_proteome()` function to import the proteome into the R session. Users can choose to work with the proteome sequence in R either as [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) object (`obj.type = "Biostrings"`) or [data.table](https://github.com/Rdatatable/data.table/wiki) object
(`obj.type = "data.table"`) by specifying the `obj.type` argument of the `read_proteome()` function. 


```{r,eval=FALSE}
# import proteome as Biostrings object
Human_Proteome <- read_proteome(file = HS.proteome.refseq)
```

```{r,eval=FALSE}
Human_Proteome
```

```
A AAStringSet instance of length 113620
         width seq                          names               
     [1]  1474 MGKNKLLHPSLVL...YNAPCSKDLGNA NP_000005.2 alpha...
     [2]   290 MDIEAYFERIGYK...LVPKPGDGSLTI NP_000006.2 aryla...
     [3]   421 MAAGFGRCCRVLR...IVAREHIDKYKN NP_000007.1 mediu...
     [4]   412 MAAALLARASGPA...VIAGHLLRSYRS NP_000008.1 short...
     [5]   655 MQAARMAASLGRQ...RGGVVTSNPLGF NP_000009.1 very ...
     ...   ... ...
[113616]    98 MPLIYMNIMLAFT...LDYVHNLNLLQC YP_003024034.1 NA...
[113617]   459 MLKLIVPTIMLLP...SLNPDIITGFSS YP_003024035.1 NA...
[113618]   603 MTMHTTMTTLTLT...FFPLILTLLLIT YP_003024036.1 NA...
[113619]   174 MMYALFLLSVGLV...GVYIVIEIARGN YP_003024037.1 NA...
[113620]   380 MTPMRKTNPLMKL...ISLIENKMLKWA YP_003024038.1 cy...
```

Alternatively, users can perform the pipeline logic of the [magrittr](https://github.com/smbache/magrittr) package:

```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import proteome as Biostrings object
Human_Proteome <- getProteome( db       = "refseq", 
                               organism = "Homo sapiens",
                               path     = file.path("_ncbi_downloads","proteomes")) %>%
    read_proteome()
```

```{r,eval=FALSE}
Human_Proteome
```

```
 A AAStringSet instance of length 113620
         width seq                          names               
     [1]  1474 MGKNKLLHPSLVL...YNAPCSKDLGNA NP_000005.2 alpha...
     [2]   290 MDIEAYFERIGYK...LVPKPGDGSLTI NP_000006.2 aryla...
     [3]   421 MAAGFGRCCRVLR...IVAREHIDKYKN NP_000007.1 mediu...
     [4]   412 MAAALLARASGPA...VIAGHLLRSYRS NP_000008.1 short...
     [5]   655 MQAARMAASLGRQ...RGGVVTSNPLGF NP_000009.1 very ...
     ...   ... ...
[113616]    98 MPLIYMNIMLAFT...LDYVHNLNLLQC YP_003024034.1 NA...
[113617]   459 MLKLIVPTIMLLP...SLNPDIITGFSS YP_003024035.1 NA...
[113618]   603 MTMHTTMTTLTLT...FFPLILTLLLIT YP_003024036.1 NA...
[113619]   174 MMYALFLLSVGLV...GVYIVIEIARGN YP_003024037.1 NA...
[113620]   380 MTPMRKTNPLMKL...ISLIENKMLKWA YP_003024038.1 cy...
```

### Example `NCBI Genbank`:

```{r,eval=FALSE}
# download the proteome of Homo sapiens from genbank
# and store the corresponding proteome file in '_ncbi_downloads/proteomes'
HS.proteome.genbank <- getProteome( db       = "genbank", 
                                   organism = "Homo sapiens",
                                   path     = file.path("_ncbi_downloads","proteomes"))
```


```{r,eval=FALSE}
# import proteome as Biostrings object
Human_Proteome <- read_proteome(file = HS.proteome.genbank)
```

```{r,eval=FALSE}
Human_Proteome
```

```
A AAStringSet instance of length 13
     width seq                              names               
 [1]   318 MPMANLLLLIVPILI...VSMPITISSIPPQT AAB58943.1 NADH d...
 [2]   347 MNPLAQPVIYSTIFA...TLLLPISPFMLMIL AAB58944.1 NADH d...
 [3]   513 MFADRWLFSTNHKDI...PPYHTFEEPVYMKS AAB58945.1 cytoch...
 [4]   227 MAHAAQVGLQDATSP...IPLKIFEMGPVFTL AAB58946.1 cytoch...
 [5]    68 MPQLNTTVWPTMITP...WTKICSLHSLPPQS AAB58947.1 ATPase...
 ...   ... ...
 [9]    98 MPLIYMNIMLAFTIS...YGLDYVHNLNLLQC AAB58951.1 NADH d...
[10]   459 MLKLIVPTIMLLPLT...LLSLNPDIITGFSS AAB58952.1 NADH d...
[11]   603 MTMHTTMTTLTLTSL...SFFFPLILTLLLIT AAB58953.1 NADH d...
[12]   174 MMYALFLLSVGLVMG...FVGVYIVIEIARGN AAB58954.1 NADH d...
[13]   380 MTPMRKTNPLMKLIN...PTISLIENKMLKWA AAB58955.3 cytoch...
```

### Example `ENSEMBL`:

```{r,eval=FALSE}
# download the proteome of Homo sapiens from ENSEMBL
# and store the corresponding proteome file in '_ncbi_downloads/proteomes'
HS.proteome.ensembl <- getProteome( db       = "ensembl", 
                                   organism = "Homo sapiens",
                                   path     = file.path("_ncbi_downloads","proteomes"))
```

```{r,eval=FALSE}
# import proteome as Biostrings object
Human_Proteome <- read_proteome(file = HS.proteome.ensembl)
```

```{r,eval=FALSE}
Human_Proteome
```

```
  A AAStringSet instance of length 107844
         width seq                          names               
     [1]     3 PSY                          ENSP00000451515.1...
     [2]     4 TGGY                         ENSP00000452494.1...
     [3]     2 EI                           ENSP00000451042.1...
     [4]     4 GTGG                         ENSP00000487941.1...
     [5]     4 GTGG                         ENSP00000488240.1...
     ...   ... ...
[107840]   135 MLQKKIEEKDLKV...LNHICKVPLAIK ENSP00000495237.1...
[107841]   166 MEHAFTPLEPLLS...IIQEESLIYLLQ ENSP00000496198.1...
[107842]    42 MEVPTAYMISPKE...GLKSENTMLLRC ENSP00000495723.1...
[107843]   508 MPSMLERISKNLV...GTLSLLQQLAEA ENSP00000496548.1...
[107844]   508 MPSMLERISKNLV...GTLSLLQQLAEA ENSP00000494855.1...
```

Users can also specify the `release` argument which denotes the database release version of ENSEMBL (`db = "ensembl"`) in case they wish to download archived vgenome versions. Default is release = NULL meaning that the most recent database version is used.

### Example Retrieval `Uniprot`:

Another  way of retrieving proteome sequences is from the [UniProt](http://www.uniprot.org/) database.

```{r,eval=FALSE}
# download the proteome of Mus musculus from UniProt
# and store the corresponding proteome file in '_uniprot_downloads/proteomes'
Mm.proteome.uniprot<- getProteome( db       = "uniprot", 
                                   organism = "Mus musculus",
                                   path     = file.path("_uniprot_downloads","proteomes"))
```


```{r,eval=FALSE}
# import proteome as Biostrings object
Mouse_Proteome <- read_proteome(file = Mm.proteome.uniprot)
```

```{r,eval=FALSE}
Mouse_Proteome
```

```
A AAStringSet instance of length 84522
        width seq                                  names               
    [1]   781 MATQADLMELDMAMEPD...LPPGDSNQLAWFDTDL sp|Q02248|CTNB1_M...
    [2]  2531 MPRLLTPLLCLTLLPAL...PTTMPSQITHIPEAFK sp|Q01705|NOTC1_M...
    [3]   437 MLLLLARCFLVILASSL...LDSETMHPLGMAVKSS sp|Q62226|SHH_MOU...
    [4]   380 MKKPIGILSPGVALGTA...VKCKKCTEIVDQFVCK sp|P22725|WNT5A_M...
    [5]   387 MEESQSDISLELPLSQE...SRHKKTMVKKVGPDSD sp|P02340|P53_MOU...
    ...   ... ...
[84518]   459 MLKIILPSLMLLPLTWL...LILLTTSPKLITGLTM tr|A0A0F6PZ84|A0A...
[84519]   380 MTNMRKTHPLFKIINHS...LMPISGIIEDKMLKLY tr|U3TEV9|U3TEV9_...
[84520]   381 MTNMRKTHPLFKIINHS...MPISGIIEDKMLKLYP tr|A0A0F6PXN8|A0A...
[84521]   172 MNNYIFVLSSLFLVGCL...WSLFAGIFIIIEITRD tr|A0A0F6PXK9|A0A...
[84522]   381 MTNMRKTHPLFKIINHS...MPISGIIEDKMLKLYP tr|A0A0F6PYR5|A0A...
```

### ProteomeSet Retrieval

The `getProteomeSet()` function enables users to retrieve proteome files for
multiple species. This is convenient when users wish to bulk download a particular
set of species. Internally, a folder named `set_proteomes` is generated in which
proteomes and proteome info files are stored. In addition, users can specify the arguments: `clean_retrieval` and `gunzip` (both are `TRUE` by default) to clean file names and automatically unzip downloaded files. 

Example:

```r
# specify the species names
download_species <- c("Arabidopsis thaliana", 
                      "Arabidopsis lyrata", 
                      "Capsella rubella")
# retrieve these three species from NCBI RefSeq                       
biomartr::getProteomeSet("refseq", organisms = download_species, path = "set_proteomes")
```

If the download process was interrupted, users can re-run the function
and it will only download missing genomes. In cases users wish
to download everything again and updating existing genomes, they may
specify the argument `update = TRUE`.

### CDS Retrieval

The `getCDS()` function is an interface function to the [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/), [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/),
[ENSEMBL](http://www.ensembl.org/index.html) databases from
which corresponding CDS files can be retrieved. It works analogous to `getGenome()` and `getProteome()`. 

The `db` argument specifies from which database proteomes in `*.fasta` file format shall be retrieved.

Options are:

- `db = "refseq"` for retrieval from [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/)
- `db = "genbank"` for retrieval from [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/)
- `db = "ensembl"` for retrieval from [ENSEMBL](http://www.ensembl.org/index.html)


Furthermore, again users need to specify the scientific name of the organism of interest for which a proteomes shall be downloaded, e.g. `organism = "Homo sapiens"`.
Finally, the `path` argument specifies the folder path in which the corresponding CDS file shall be locally stored. In case users would like to store the CDS file at a different location,
they can specify the `path = file.path("put","your","path","here")` argument.

### Example `NCBI RefSeq`:

```{r,eval=FALSE}
# download the genome of Homo sapiens from refseq
# and store the corresponding genome CDS file in '_ncbi_downloads/CDS'
HS.cds.refseq <- getCDS( db       = "refseq", 
                         organism = "Homo sapiens",
                         path     = file.path("_ncbi_downloads","CDS"))
```


In this example, `getCDS()` creates a directory named `'_ncbi_downloads/CDS'` into which the corresponding genome named `Homo_sapiens_cds_from_genomic_refseq.fna.gz` is downloaded. The return value of `getCDS()` is the folder path to the downloaded genome file
that can then be used as input to the `read_cds()` function. The variable `HS.cds.refseq` stores the path to the downloaded CDS file. 
Subsequently, users can use the `read_cds()` function to import the genome into the R session. Users can choose to work with the genome sequence in R either as [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) object (`obj.type = "Biostrings"`) or [data.table](https://github.com/Rdatatable/data.table/wiki) object
(`obj.type = "data.table"`) by specifying the `obj.type` argument of the `read_cds()` function. 

```{r,eval=FALSE}
# import downloaded CDS as Biostrings object
Human_CDS <- read_cds(file     = HS.cds.refseq, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_CDS
```

```
 A BStringSet instance of length 114967
          width seq                                                names               
     [1]    918 ATGGTGACTGAATTCATTTTTCTG...CACATTCTAGTGTAAAGTTTTAG lcl|NC_000001.11_...
     [2]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [3]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [4]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [5]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     ...    ... ...
[114963]    297 ATGCCCCTCATTTACATAAATATT...ACCTAAACCTACTCCAATGCTAA lcl|NC_012920.1_c...
[114964]   1378 ATGCTAAAACTAATCGTCCCAACA...CATCATTACCGGGTTTTCCTCTT lcl|NC_012920.1_c...
[114965]   1812 ATAACCATGCACACTACTATAACC...TAACCCTACTCCTAATCACATAA lcl|NC_012920.1_c...
[114966]    525 ATGATGTATGCTTTGTTTCTGTTG...TTGAGATTGCTCGGGGGAATAGG lcl|NC_012920.1_c...
[114967]   1141 ATGACCCCAATACGCAAAACTAAC...AAACAAAATACTCAAATGGGCCT lcl|NC_012920.1_c...
```

Internally, a text file named `doc_Homo_sapiens_db_refseq.txt` is generated. The information stored in this log file is structured as follows:

```
File Name: Homo_sapiens_cds_from_genomic_refseq.fna.gz
Organism Name: Homo_sapiens
Database: NCBI refseq
URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/
GCF_000001405.35_GRCh38.p9/GCF_000001405.35_GRCh38.p9_cds_from_genomic.fna.gz
Download_Date: Sun Oct 23 17:19:05 2016
refseq_category: reference genome
assembly_accession: GCF_000001405.35
bioproject: PRJNA168
biosample: NA
taxid: 9606
infraspecific_name: NA
version_status: latest
release_type: Patch
genome_rep: Full
seq_rel_date: 2016-09-26
submitter: Genome Reference Consortium
```

In summary, the `getCDS()` and `read_cds()` functions allow users to retrieve CDS files by specifying
the scientific name of the organism of interest and allow them to import the retrieved CDS file e.g. as `Biostrings` object.
Thus, users can then perform the `Biostrings notation` to work with downloaded CDS and can rely on the log
file generated by `getCDS()` to better document the source and version of CDS used for subsequent studies.

Alternatively, users can perform the pipeline logic of the [magrittr](https://github.com/smbache/magrittr) package:

```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import CDS as Biostrings object
Human_CDS <- getCDS( db       = "refseq", 
                     organism = "Homo sapiens",
                     path     = file.path("_ncbi_downloads","CDS")) %>%
    read_cds(obj.type = "Biostrings")
```

```{r,eval=FALSE}
Human_CDS
```

```
 A BStringSet instance of length 114967
          width seq                                                names               
     [1]    918 ATGGTGACTGAATTCATTTTTCTG...CACATTCTAGTGTAAAGTTTTAG lcl|NC_000001.11_...
     [2]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [3]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [4]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     [5]    402 ATGAGTGACAGCATCAACTTCTCT...CAGGACCCAGGCACAGGCATTAG lcl|NC_000001.11_...
     ...    ... ...
[114963]    297 ATGCCCCTCATTTACATAAATATT...ACCTAAACCTACTCCAATGCTAA lcl|NC_012920.1_c...
[114964]   1378 ATGCTAAAACTAATCGTCCCAACA...CATCATTACCGGGTTTTCCTCTT lcl|NC_012920.1_c...
[114965]   1812 ATAACCATGCACACTACTATAACC...TAACCCTACTCCTAATCACATAA lcl|NC_012920.1_c...
[114966]    525 ATGATGTATGCTTTGTTTCTGTTG...TTGAGATTGCTCGGGGGAATAGG lcl|NC_012920.1_c...
[114967]   1141 ATGACCCCAATACGCAAAACTAAC...AAACAAAATACTCAAATGGGCCT lcl|NC_012920.1_c...
```

### Example `NCBI Genbank`:

```{r,eval=FALSE}
# download the genome of Homo sapiens from genbank
# and store the corresponding genome CDS file in '_ncbi_downloads/CDS'
HS.cds.genbank <- getCDS( db       = "genbank", 
                         organism = "Homo sapiens",
                         path     = file.path("_ncbi_downloads","CDS"))
```

```{r,eval=FALSE}
# import downloaded CDS as Biostrings object
Human_CDS <- read_cds(file     = HS.cds.genbank, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_CDS
```

```
  A BStringSet instance of length 13
     width seq                                                     names               
 [1]   956 ATACCCATGGCCAACCTCCTACTCCT...ATCTCCAGCATTCCCCCTCAAACCTA lcl|J01415.2_cds_...
 [2]  1042 ATTAATCCCCTGGCCCAACCCGTCAT...CTCCCCTTTTATACTAATAATCTTAT lcl|J01415.2_cds_...
 [3]  1542 ATGTTCGCCGACCGTTGACTATTCTC...AAGAACCCGTATACATAAAATCTAGA lcl|J01415.2_cds_...
 [4]   684 ATGGCACATGCAGCGCAAGTAGGTCT...AAATAGGGCCCGTATTTACCCTATAG lcl|J01415.2_cds_...
 [5]   207 ATGCCCCAACTAAATACTACCGTATG...TTCATTCATTGCCCCCACAATCCTAG lcl|J01415.2_cds_...
 ...   ... ...
 [9]   297 ATGCCCCTCATTTACATAAATATTAT...ATAACCTAAACCTACTCCAATGCTAA lcl|J01415.2_cds_...
[10]  1378 ATGCTAAAACTAATCGTCCCAACAAT...CGACATCATTACCGGGTTTTCCTCTT lcl|J01415.2_cds_...
[11]  1812 ATAACCATGCACACTACTATAACCAC...TCCTAACCCTACTCCTAATCACATAA lcl|J01415.2_cds_...
[12]   525 ATGATGTATGCTTTGTTTCTGTTGAG...TAATTGAGATTGCTCGGGGGAATAGG lcl|J01415.2_cds_...
[13]  1141 ATGACCCCAATACGCAAAACTAACCC...TGAAAACAAAATACTCAAATGGGCCT lcl|J01415.2_cds_...
```

### Example `ENSEMBL`:

```{r,eval=FALSE}
# download the genome of Homo sapiens from ensembl
# and store the corresponding genome CDS file in '_ncbi_downloads/CDS'
HS.cds.ensembl <- getCDS( db       = "ensembl", 
                         organism = "Homo sapiens",
                         path     = file.path("_ncbi_downloads","CDS"))
```

```{r,eval=FALSE}
# import downloaded CDS as Biostrings object
Human_CDS <- read_cds(file     = HS.cds.ensembl, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_CDS
```

```
  A BStringSet instance of length 102915
          width seq                                                names               
     [1]     13 ACTGGGGGATACG                                      ENST00000448914.1...
     [2]     12 GGGACAGGGGGC                                       ENST00000631435.1...
     [3]     12 GGGACAGGGGGC                                       ENST00000632684.1...
     [4]      9 CCTTCCTAC                                          ENST00000434970.2...
     [5]      8 GAAATAGT                                           ENST00000415118.1...
     ...    ... ...
[102911]   1665 ATGCTACTGCCACTGCTGCTGTCC...ATGCAGAAGTCAAGTTCCAATGA ENST00000436984.6...
[102912]   1920 ATGCTACTGCCACTGCTGCTGTCC...ATGCAGAAGTCAAGTTCCAATGA ENST00000439889.6...
[102913]   2094 ATGCTACTGCCACTGCTGCTGTCC...ATGCAGAAGTCAAGTTCCAATGA ENST00000339313.9...
[102914]    466 ATGCTACTGCCACTGCTGCTGTCC...AGCCCTGGACCTCTCTGTGCAGT ENST00000529627.1...
[102915]    559 ATGCGGAGATGCTACTGCCACTGC...CCTCACCTGCCATGTGGACTTCT ENST00000530476.1...
```

Users can also specify the `release` argument which denotes the database release version of ENSEMBL (`db = "ensembl"`) in case they wish to download archived vgenome versions. Default is release = NULL meaning that the most recent database version is used.

### CDSSet Retrieval

The `getCDSSet()` function enables users to retrieve CDS files for
multiple species. This is convenient when users wish to bulk download a particular
set of species. Internally, a folder named `set_cds` is generated in which
CDS and CDS info files are stored. In addition, users can specify the arguments: `clean_retrieval` and `gunzip` (both are `TRUE` by default) to clean file names and automatically unzip downloaded files. 

Example:

```r
# specify the species names
download_species <- c("Arabidopsis thaliana", 
                      "Arabidopsis lyrata", 
                      "Capsella rubella")
# retrieve these three species from NCBI RefSeq                       
biomartr::getCDSSet("refseq", organisms = download_species, path = "set_cds")
```

If the download process was interrupted, users can re-run the function
and it will only download missing genomes. In cases users wish
to download everything again and updating existing genomes, they may
specify the argument `update = TRUE`.

### RNA Retrieval

The `getRNA()` function is an interface function to the [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/), [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/),
[ENSEMBL](http://www.ensembl.org/index.html) databases from
which corresponding RNA files can be retrieved. It works analogous to `getGenome()`, `getProteome()`, and `getCDS()`. 

The `db` argument specifies from which database proteomes in `*.fasta` file format shall be retrieved.

Options are:

- `db = "refseq"` for retrieval from [NCBI RefSeq](http://www.ncbi.nlm.nih.gov/refseq/about/)
- `db = "genbank"` for retrieval from [NCBI Genbank](http://www.ncbi.nlm.nih.gov/genbank/about/)
- `db = "ensembl"` for retrieval from [ENSEMBL](http://www.ensembl.org/index.html)

Furthermore, again users need to specify the scientific name of the organism of interest for which a proteomes shall be downloaded, e.g. `organism = "Homo sapiens"`.
Finally, the `path` argument specifies the folder path in which the corresponding RNA file shall be locally stored. In case users would like to store the RNA file at a different location,
they can specify the `path = file.path("put","your","path","here")` argument.

### Example `NCBI RefSeq`:

```{r,eval=FALSE}
# download the RNA of Homo sapiens from refseq
# and store the corresponding RNA file in '_ncbi_downloads/RNA'
HS.rna.refseq <- getRNA( db       = "refseq", 
                         organism = "Homo sapiens",
                         path     = file.path("_ncbi_downloads","RNA"))
```


In this example, `getRNA()` creates a directory named `'_ncbi_downloads/RNA'` into which the corresponding RNA file named `Homo_sapiens_rna_from_genomic_refseq.fna.gz` is downloaded. The return value of `getRNA()` is the folder path to the downloaded genome file
that can then be used as input to the `read_rna()` function. The variable `HS.rna.refseq` stores the path to the downloaded RNA file. 
Subsequently, users can use the `read_cds()` function to import the genome into the R session. Users can choose to work with the genome sequence in R either as [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) object (`obj.type = "Biostrings"`) or [data.table](https://github.com/Rdatatable/data.table/wiki) object
(`obj.type = "data.table"`) by specifying the `obj.type` argument of the `read_rna()` function. 

```{r,eval=FALSE}
# import downloaded RNA as Biostrings object
Human_rna <- read_rna(file     = HS.rna.refseq, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_rna
```

```
   A BStringSet instance of length 164136
          width seq                                                                                       names               
     [1]   1652 CTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGT...CACAGCTAGAGATCCTTTATTAAAAGCACACTGTTGGTTTCTG lcl|NC_000001.11_...
     [2]   1769 TCCGGCAGAGCGGAAGCGGCGGCGGGAGCTTCCGGGAGGGCGG...ACCAACAGTGTGCTTTTAATAAAGGATCTCTAGCTGTGCAGGA lcl|NC_000001.11_...
     [3]     68 TGTGGGAGAGGAACATGGGCTCAGGACAGCGGGTGTCAGCTTGCCTGACCCCCATGTCGCCTCTGTAG                      lcl|NC_000001.11_...
     [4]     23 TGACCCCCATGTCGCCTCTGTAG                                                                   lcl|NC_000001.11_...
     [5]     23 GAGAGGAACATGGGCTCAGGACA                                                                   lcl|NC_000001.11_...
     ...    ... ...
[164132]     59 GAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCA                               lcl|NC_012920.1_t...
[164133]     71 ACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTA                   lcl|NC_012920.1_t...
[164134]     69 GTTCTTGTAGTTGAAATACAACGATGGTTTTTCATATCATTGGTCGTGGTTGTAGTCCGTGCGAGAATA                     lcl|NC_012920.1_t...
[164135]     66 GTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACA                        lcl|NC_012920.1_t...
[164136]     68 CAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGA                      lcl|NC_012920.1_t...
```

Internally, a text file named `doc_Homo_sapiens_db_refseq.txt` is generated. The information stored in this log file is structured as follows:

```
File Name: Homo_sapiens_rna_from_genomic_refseq.fna.gz
Organism Name: Homo_sapiens
Database: NCBI refseq
URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_rna_from_genomic.fna.gz
Download_Date: Wed Mar 15 16:46:45 2017
refseq_category: reference genome
assembly_accession: GCF_000001405.36
bioproject: PRJNA168
biosample: NA
taxid: 9606
infraspecific_name: NA
version_status: latest
release_type: Patch
genome_rep: Full
seq_rel_date: 2017-01-06
submitter: Genome Reference Consortium
```

In summary, the `getRNA()` and `read_rna()` functions allow users to retrieve RNA files by specifying
the scientific name of the organism of interest and allow them to import the retrieved RNA file e.g. as `Biostrings` object.
Thus, users can then perform the `Biostrings notation` to work with downloaded RNA and can rely on the log
file generated by `getRNA()` to better document the source and version of RNA used for subsequent studies.

Alternatively, users can perform the pipeline logic of the [magrittr](https://github.com/smbache/magrittr) package:

```{r,eval=FALSE}
# install.packages("magrittr")
library(magrittr)
# import RNA as Biostrings object
Human_rna <- getRNA( db       = "refseq", 
                     organism = "Homo sapiens",
                     path     = file.path("_ncbi_downloads","RNA")) %>%
    read_cds(obj.type = "Biostrings")
```

```{r,eval=FALSE}
Human_rna
```

```
   A BStringSet instance of length 164136
          width seq                                                                                       names               
     [1]   1652 CTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGT...CACAGCTAGAGATCCTTTATTAAAAGCACACTGTTGGTTTCTG lcl|NC_000001.11_...
     [2]   1769 TCCGGCAGAGCGGAAGCGGCGGCGGGAGCTTCCGGGAGGGCGG...ACCAACAGTGTGCTTTTAATAAAGGATCTCTAGCTGTGCAGGA lcl|NC_000001.11_...
     [3]     68 TGTGGGAGAGGAACATGGGCTCAGGACAGCGGGTGTCAGCTTGCCTGACCCCCATGTCGCCTCTGTAG                      lcl|NC_000001.11_...
     [4]     23 TGACCCCCATGTCGCCTCTGTAG                                                                   lcl|NC_000001.11_...
     [5]     23 GAGAGGAACATGGGCTCAGGACA                                                                   lcl|NC_000001.11_...
     ...    ... ...
[164132]     59 GAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCA                               lcl|NC_012920.1_t...
[164133]     71 ACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTA                   lcl|NC_012920.1_t...
[164134]     69 GTTCTTGTAGTTGAAATACAACGATGGTTTTTCATATCATTGGTCGTGGTTGTAGTCCGTGCGAGAATA                     lcl|NC_012920.1_t...
[164135]     66 GTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACA                        lcl|NC_012920.1_t...
[164136]     68 CAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGA                      lcl|NC_012920.1_t...
```

### Example `NCBI Genbank`:

```{r,eval=FALSE}
# download the RNA of Homo sapiens from genbank
# and store the corresponding genome RNA file in '_ncbi_downloads/RNA'
HS.rna.genbank <- getRNA( db       = "genbank", 
                         organism = "Homo sapiens",
                         path     = file.path("_ncbi_downloads","RNA"))
```

```{r,eval=FALSE}
# import downloaded RNA as Biostrings object
Human_rna <- read_cds(file     = HS.rna.genbank, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_rna
```

```
  A BStringSet instance of length 24
     width seq                                                                                            names               
 [1]    71 GTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACA                        lcl|J01415.2_trna...
 [2]   954 AATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACA...AAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAAC lcl|J01415.2_rrna...
 [3]    69 CAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGA                          lcl|J01415.2_trna...
 [4]  1559 GCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCT...ATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTT lcl|J01415.2_rrna...
 [5]    75 GTTAAGATGGCAGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTTAACA                    lcl|J01415.2_trna...
 ...   ... ...
[20]    59 GAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCA                                    lcl|J01415.2_trna...
[21]    71 ACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTA                        lcl|J01415.2_trna...
[22]    69 GTTCTTGTAGTTGAAATACAACGATGGTTTTTCATATCATTGGTCGTGGTTGTAGTCCGTGCGAGAATA                          lcl|J01415.2_trna...
[23]    66 GTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACA                             lcl|J01415.2_trna...
[24]    68 CAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGA                           lcl|J01415.2_trna...
```

### Example `ENSEMBL`:

```{r,eval=FALSE}
# download the RNA of Homo sapiens from ensembl
# and store the corresponding genome RNA file in '_ncbi_downloads/RNA'
HS.rna.ensembl <- getRNA( db       = "ensembl", 
                          organism = "Homo sapiens",
                          path     = file.path("_ncbi_downloads","RNA"))
```

```{r,eval=FALSE}
# import downloaded RNA as Biostrings object
Human_rna <- read_cds(file     = HS.rna.ensembl, 
                      obj.type = "Biostrings")

```

```{r,eval=FALSE}
# look at the Biostrings object
Human_rna
```

```
  A BStringSet instance of length 36701
         width seq                                                                                        names               
    [1]    104 GTGCTCACTTTGGCAACATACATACTAAAATTGGACGGATACAG...GCACAAGGATGACATGCAAATTCATGAAGCATTCCATATTTTT ENST00000516494.2...
    [2]    164 TTAACTACCTGACAGAGGAGATACTGTGATCATGAAAGTGGTTT...CATAATTTCTGGTGGTAGGGGACTGCGTTCATGTTCTCCCCTA ENST00000627793.1...
    [3]    114 ACACTGGTTTCTCTTCAGATCGAATAAATCTTTCGCCTTTTACT...AGTTATAAGCTAATTTTTTGTAAGCCTTGCCCTGGGGAGGCAG ENST00000629478.1...
    [4]     64 CAGTGTTACACCTCTTTTAGAATTTATCTATCAGGTTTTCCAGTGTTCACTGAAATTTGTCTCT                           ENST00000629187.1...
    [5]    107 GTGTTGGCCTGGGCAGCACGTATACTAAAGTTGGAATGACACAG...GCGCAAGGATGATGTGCAAATTCGTGACAAGTTCCATATTTTT ENST00000631612.1...
    ...    ... ...
[36697]     90 GGCTAGTCCAAATGTAGTGGTGTTCCAAACTAATTAATCACAACCAGTTACAGATTTCTTGTTTTCTTTTCCACTCACACTTAGCCTTAA ENST00000410951.1...
[36698]    109 GGCTGATCTGAAGATGATGAGTTATCTCAATTGATTGTTCAGCC...TCTATTCTTTCCTCTCTTCTCACTACTGCACTTGGCTAGGAAA ENST00000410462.2...
[36699]    109 GGTGGGTCCAAAGGTAGCAAGTTATCTCAATTGATCACAGTCAG...CATTCTATCACCCCTTCTCATTACTGTACTTGACTAGTCTTTT ENST00000364501.1...
[36700]    293 GGATATGAGGGCGATCTGGCAGGGACATCTGTCACCCCACTGAT...AAAATTAGCTGGGCATAGTGGCGTGCACCTGTCGTCCTAGCTA ENST00000365097.1...
[36701]    172 TTGAAGGCGTGGAGACTGAAGTCCTCTCTATATCCACAGAACAG...AAGAGGGCTGTTCAGTCTCCATGCCCTTCAATCCTTGGCTACT ENST00000615087.1...
```

Users can also specify the `release` argument which denotes the database release version of ENSEMBL (`db = "ensembl"`) in case they wish to download archived vgenome versions. Default is release = NULL meaning that the most recent database version is used.

### RNASet Retrieval

The `getRNASet()` function enables users to retrieve RNA files for
multiple species. This is convenient when users wish to bulk download a particular
set of species. Internally, a folder named `set_rna` is generated in which
RNA and RNA info files are stored. In addition, users can specify the arguments: `clean_retrieval` and `gunzip` (both are `TRUE` by default) to clean file names and automatically unzip downloaded files. 

Example:

```r
# specify the species names
download_species <- c("Arabidopsis thaliana", 
                      "Arabidopsis lyrata", 
                      "Capsella rubella")
# retrieve these three species from NCBI RefSeq                       
biomartr::getRNASet("refseq", organisms = download_species, path = "set_rna")
```

If the download process was interrupted, users can re-run the function
and it will only download missing genomes. In cases users wish
to download everything again and updating existing genomes, they may
specify the argument `update = TRUE`.

### Retrieve the annotation file of a particular genome

Finally, users can download the corresponding annotation `.gff` files for particular genomes of interest using the `getGFF()` or alternatively for `ensembl` the `getGTF()` function.

### Example `NCBI RefSeq`:

```{r,eval=FALSE}
# download the GFF file of Homo sapiens from refseq
# and store the corresponding file in '_ncbi_downloads/annotation'
HS.gff.refseq <- getGFF( db       = "refseq", 
                         organism = "Homo sapiens", 
                         path = file.path("_ncbi_downloads","annotation"))
```

After downloading the `.gff` file, users can import the `.gff` file with `read_gff()`.

```{r,eval=FALSE}
# import downloaded GFF file
Human_GFF <- read_gff(file = HS.gff.refseq)
```

```{r,eval=FALSE}
Human_GFF
```

```
          seqid     source       type start       end score strand phase
          <chr>      <chr>      <chr> <int>     <int> <dbl>  <chr> <dbl>
1  NC_000001.11     RefSeq     region     1 248956422     0      +     0
2  NC_000001.11 BestRefSeq       gene 11874     14409     0      +     0
3  NC_000001.11 BestRefSeq transcript 11874     14409     0      +     0
4  NC_000001.11 BestRefSeq       exon 11874     12227     0      +     0
5  NC_000001.11 BestRefSeq       exon 12613     12721     0      +     0
6  NC_000001.11 BestRefSeq       exon 13221     14409     0      +     0
7  NC_000001.11 BestRefSeq       gene 14362     29370     0      -     0
8  NC_000001.11 BestRefSeq transcript 14362     29370     0      -     0
9  NC_000001.11 BestRefSeq       exon 29321     29370     0      -     0
10 NC_000001.11 BestRefSeq       exon 24738     24891     0      -     0
```

#### Removing corrupt lines from downloaded GFF files

In some cases, `GFF` files stored at NCBI databases include corrupt lines that have more than 65000 characters. This [leads to problems](https://github.com/lawremi/rtracklayer/issues/15)
when trying to import such annotation files into R. To overcome this issue users can specify the `remove_annotation_outliers = TRUE`
argument to remove such outlier lines and overwrite the downloaded annotation file. This will make any downstream analysis with this annotation file much more reliable.

Example:

```r
Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana", remove_annotation_outliers = TRUE)
```

```
Starting GFF retrieval of 'Arabidopsis thaliana' from refseq ...


|============================================| 100%   52 MB
File _ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff.gz exists already. Thus, download has been skipped.
Importing '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff.gz' ...
|============================================| 100%  434 MB
Reading annotation file '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff.gz' and removing all outlier lines with number of characters greater 65000 ...
Overwriting '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff.gz' with removed outlier lines ...
Unzipping file _ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff.gz' ...
The new annotation file was created and has been stored at '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_refseq.gff'.
The outlier lines were stored in a separate file to explore at '/var/folders/3x/6bbw6ds1039gpwny1m0hn8r80000gp/T//RtmpuVjnsC/Arabidopsis_thaliana_genomic_refseq.gff.gz_outlier_lines.txt'.
```


### Example `NCBI Genbank`:

```{r,eval=FALSE}
# download the GFF file of Homo sapiens from genbank
# and store the corresponding file in '_ncbi_downloads/annotation'
HS.gff.genbank <- getGFF( db       = "genbank", 
                         organism = "Homo sapiens", 
                         path = file.path("_ncbi_downloads","annotation"))
```

After downloading the `.gff` file, users can import the `.gff` file with `read_gff()`.

```{r,eval=FALSE}
# import downloaded GFF file
Human_GFF <- read_gff(file = HS.gff.genbank)
```

```{r,eval=FALSE}
# show all elements of the data.frame
# options(tibble.print_max = Inf)
Human_GFF
```

```
        seqid  source       type     start       end score strand phase
        <chr>   <chr>      <chr>     <int>     <int> <dbl>  <chr> <dbl>
1  CM000663.2 Genbank     region         1 248956422     0      +     0
2  CM000663.2 Genbank centromere 122026460 125184587     0      +     0
3  KI270706.1 Genbank     region         1    175055     0      +     0
4  KI270707.1 Genbank     region         1     32032     0      +     0
5  KI270708.1 Genbank     region         1    127682     0      +     0
6  KI270709.1 Genbank     region         1     66860     0      +     0
7  KI270710.1 Genbank     region         1     40176     0      +     0
8  KI270711.1 Genbank     region         1     42210     0      +     0
9  KI270712.1 Genbank     region         1    176043     0      +     0
10 KI270713.1 Genbank     region         1     40745     0      +     0
```

#### Removing corrupt lines from downloaded GFF files

In some cases, `GFF` files stored at NCBI databases include corrupt lines that have more than 65000 characters. This [leads to problems](https://github.com/lawremi/rtracklayer/issues/15)
when trying to import such annotation files into R. To overcome this issue users can specify the `remove_annotation_outliers = TRUE`
argument to remove such outlier lines and overwrite the downloaded annotation file. This will make any downstream analysis with this annotation file much more reliable.

Example:

```r
Ath_path <- biomartr::getGFF(db = "genbank",
             organism = "Arabidopsis thaliana",
             remove_annotation_outliers = TRUE)
```

```
Starting GFF retrieval of 'Arabidopsis thaliana' from genbank ...

Completed!
Now continue with species download ...
GFF download of Arabidopsis_thaliana is completed!
Checking md5 hash of file: _ncbi_downloads/annotation/Arabidopsis_thaliana_md5checksums.txt ...
The md5 hash of file '_ncbi_downloads/annotation/Arabidopsis_thaliana_md5checksums.txt' matches!
Importing '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_genbank.gff.gz' ...
Reading annotation file '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_genbank.gff.gz' and removing all outlier lines with number of characters greater 65000 ...
Overwriting '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_genbank.gff.gz' with removed outlier lines ...
Unzipping file _ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_genbank.gff.gz' ...
The new annotation file was created and has been stored at '_ncbi_downloads/annotation/Arabidopsis_thaliana_genomic_genbank.gff'.
The outlier lines were stored in a separate file to explore at '/var/folders/3x/6bbw6ds1039gpwny1m0hn8r80000gp/T//RtmpuVjnsC/Arabidopsis_thaliana_genomic_genbank.gff.gz_outlier_lines.txt'.
```

### Example `ENSEMBL`:

```{r,eval=FALSE}
# download the GFF file of Homo sapiens from ENSEMBL
# and store the corresponding file in 'ensembl/annotation'
HS.gff.ensembl <- getGFF( db       = "ensembl", 
                         organism = "Homo sapiens", 
                         path = file.path("ensembl","annotation"))
```

After downloading the `.gff` file, users can import the `.gff` file with `read_gff()`.

```{r,eval=FALSE}
# import downloaded GFF file
Human_GFF <- read_gff(file = HS.gff.ensembl)
```

```{r,eval=FALSE}
# show all elements of the data.frame
# options(tibble.print_max = Inf)
Human_GFF
```

```
   seqid source              type start       end   score strand phase
   <int>  <chr>             <chr> <int>     <int>   <chr>  <chr> <dbl>
1      1 GRCh38        chromosome     1 248956422       .      .     0
2      1      . biological_region 10469     11240 1.3e+03      .     0
3      1      . biological_region 10650     10657   0.999      +     0
4      1      . biological_region 10655     10657   0.999      -     0
5      1      . biological_region 10678     10687   0.999      +     0
6      1      . biological_region 10681     10688   0.999      -     0
7      1      . biological_region 10707     10716   0.999      +     0
8      1      . biological_region 10708     10718   0.999      -     0
9      1      . biological_region 10735     10747   0.999      -     0
10     1      . biological_region 10737     10744   0.999      +     0
```

Users can also specify the `release` argument which denotes the database release version of ENSEMBL (`db = "ensembl"`) in case they wish to download archived vgenome versions. Default is release = NULL meaning that the most recent database version is used.


#### Removing corrupt lines from downloaded GFF files

In some cases, `GFF` files stored at NCBI databases include corrupt lines that have more than 65000 characters. This [leads to problems](https://github.com/lawremi/rtracklayer/issues/15)
when trying to import such annotation files into R. To overcome this issue users can specify the `remove_annotation_outliers = TRUE`
argument to remove such outlier lines and overwrite the downloaded annotation file. This will make any downstream analysis with this annotation file much more reliable.


Alternatively for `getGTF()`:

```{r,eval=FALSE}
# download the GTF file of Homo sapiens from ENSEMBL
# and store the corresponding file in 'ensembl/annotation'
HS.gtf.ensembl <- getGTF( db       = "ensembl", 
                         organism = "Homo sapiens", 
                         path = file.path("ensembl","annotation"),
                         assembly_type = "primary_assembly")
```

Taken together, `getGFF()` or `getGTF()` in combination with `getGenome()`, `getProteome()`, `getRNA()` and `getCDS()` allows users to retrieve the genome information together with the
corresponding `.gff` or `gtf` annotation file to make sure that they both have the same version and origin.


### GFFSet Retrieval

The `getGFFSet()` function enables users to retrieve GFF files for
multiple species. This is convenient when users wish to bulk download a particular
set of species. Internally, a folder named `set_gff` is generated in which
GFF and GFF info files are stored. In addition, users can specify the arguments: `clean_retrieval` and `gunzip` (both are `TRUE` by default) to clean file names and automatically unzip downloaded files. 

Example:

```r
# specify the species names
download_species <- c("Arabidopsis thaliana", 
                      "Arabidopsis lyrata", 
                      "Capsella rubella")
# retrieve these three species from NCBI RefSeq                       
biomartr::getGFFSet("refseq", organisms = download_species, path = "set_gff")
```

If the download process was interrupted, users can re-run the function
and it will only download missing genomes. In cases users wish
to download everything again and updating existing genomes, they may
specify the argument `update = TRUE`.

In some cases, `GFF` files stored at NCBI databases include corrupt lines that have more than 65000 characters. This [leads to problems](https://github.com/lawremi/rtracklayer/issues/15)
when trying to import such annotation files into R. To overcome this issue users can specify the `remove_annotation_outliers = TRUE`
argument to remove such outlier lines and overwrite the downloaded annotation file. This will make any downstream analysis with this annotation file much more reliable.


## Repeat Masker Retrieval

[Repeat Masker](http://www.repeatmasker.org) is a tool for screening DNA sequences for interspersed repeats and low complexity DNA sequences.
NCBI stores the `Repeat Masker` for sevel species in their databases and can be retrieved using `getRepeatMasker()` and imported
via `read_rm()`.

### Example `NCBI RefSeq`:
```{r,eval=FALSE}
# download repeat masker annotation file for Homo sapiens
Hsapiens_rm <- getRepeatMasker( db       = "refseq", 
                                 organism = "Homo sapiens", 
                                 path = file.path("refseq","TEannotation"))
```

Now users can import the TE annotation file using `read_rm()`.

```{r,eval=FALSE}
# import TE annotation file
Hsapiens_rm_import <- read_rm("refseq/TEannotation/Homo_sapiens_rm_refseq.out.gz")
# look at data
Hsapiens_rm_import
```



## Genome Assembly Stats Retrieval

By specifying the scientific name of an organism of interest the corresponding genome assembly stats file storing the assembly statistics of the organism of interest can be downloaded and stored locally. 

### Example `NCBI RefSeq`:
```{r,eval=FALSE}
# download genome assembly stats file for Homo sapiens
Hsapiens_stats <- getAssemblyStats( db       = "refseq", 
                                 organism = "Homo sapiens", 
                                 path = file.path("refseq","AssemblyStats"))
```

Now users can import the TE annotation file using `read_rm()`.

```{r,eval=FALSE}
# import TE annotation file
Hsapiens_stats_import <- read_assemblystats(Hsapiens_stats)
# look at data
Hsapiens_stats_import
```

```
 A tibble: 1,196 x 6
   unit_name molecule_name molecule_type sequence_type statistic         value
   <chr>     <chr>         <chr>         <chr>         <chr>             <int>
 1 all       all           all           all           total-length         NA
 2 all       all           all           all           spanned-gaps        661
 3 all       all           all           all           unspanned-gaps      349
 4 all       all           all           all           region-count        317
 5 all       all           all           all           scaffold-count      875
 6 all       all           all           all           scaffold-N50   59364414
 7 all       all           all           all           scaffold-L50         17
 8 all       all           all           all           scaffold-N75   31699399
 9 all       all           all           all           scaffold-N90    7557127
10 all       all           all           all           contig-count       1536
```

### Example `NCBI Genbank`:
```{r,eval=FALSE}
# download genome assembly stats file for Homo sapiens
Hsapiens_stats <- getAssemblyStats( db       = "genbank", 
                                 organism = "Homo sapiens", 
                                 path = file.path("genbank","AssemblyStats"))
```

Now users can import the TE annotation file using `read_rm()`.

```{r,eval=FALSE}
# import TE annotation file
Hsapiens_stats_import <- read_assemblystats(Hsapiens_stats)
# look at data
Hsapiens_stats_import
```


```
 A tibble: 1,196 x 6
   unit_name molecule_name molecule_type sequence_type statistic         value
   <chr>     <chr>         <chr>         <chr>         <chr>             <int>
 1 all       all           all           all           total-length         NA
 2 all       all           all           all           spanned-gaps        661
 3 all       all           all           all           unspanned-gaps      349
 4 all       all           all           all           region-count        317
 5 all       all           all           all           scaffold-count      875
 6 all       all           all           all           scaffold-N50   59364414
 7 all       all           all           all           scaffold-L50         17
 8 all       all           all           all           scaffold-N75   31699399
 9 all       all           all           all           scaffold-N90    7557127
10 all       all           all           all           contig-count       1536
```

## Collection Retrieval

The automated retrieval of collections (= Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats)
will make sure that the genome file of an organism will match the CDS, proteome, RNA, GFF, etc file
and was generated using the same genome assembly version. One aspect of why genomics studies
fail in computational and biological reproducibility is that it is not clear whether CDS, proteome, RNA, GFF, etc files
used in a proposed analysis were generated using the same genome assembly file denoting the same genome assembly version.
To avoid this seemingly trivial mistake we encourage users to retrieve
genome file collections using the `biomartr` function `getCollection()`
and attach the corresponding output as Supplementary Data
to the respective genomics study to ensure computational and biological reproducibility.


By specifying the scientific name of an organism of interest a collection consisting of the genome file, proteome file, CDS file, RNA file, GFF file, Repeat Masker file, AssemblyStats file of the organism of interest can be downloaded and stored locally.


### Example `NCBI RefSeq`:


```{r,eval=FALSE}
# download collection for Saccharomyces cerevisiae
getCollection( db = "refseq", 
               organism = "Saccharomyces cerevisiae", 
               path = file.path("refseq","Collections"))
```

Internally, the `getCollection()` function will now generate a folder named `refseq/Collection/Saccharomyces_erevisiae`
and will store all genome and annotation files for `Saccharomyces cerevisiae` in the same folder.
In addition, the exact genoem and annotation version will be logged in the `doc` folder.

```
Genome download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The genome of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_genomic_refseq.fna.gz'.
Starting proteome retrieval of 'Saccharomyces cerevisiae' from refseq ...


Proteome download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The proteome of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_protein_refseq.faa.gz' .
Starting CDS retrieval of 'Saccharomyces cerevisiae' from refseq ...


CDS download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The genomic CDS of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_cds_from_genomic_refseq.fna.gz' .
Starting GFF retrieval of 'Saccharomyces cerevisiae' from refseq ...


GFF download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The *.gff annotation file of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_genomic_refseq.gff.gz'.
Starting RNA retrieval of 'Saccharomyces cerevisiae' from refseq ...


RNA download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The genomic RNA of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_rna_from_genomic_refseq.fna.gz' .
Starting Repeat Masker retrieval of 'Saccharomyces cerevisiae' from refseq ...


RepeatMasker download is completed!
Checking md5 hash of file: refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt ...
The md5 hash of file 'refseq/Collections/Saccharomyces_cerevisiae_md5checksums.txt' matches!
The Repeat Masker output file of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_rm_refseq.out.gz'.
Starting assembly quality stats retrieval of 'Saccharomyces cerevisiae' from refseq ...

Genome assembly quality stats file download completed!
The assembly statistics file of 'Saccharomyces_cerevisiae' has been downloaded to 'refseq/Collections' and has been named 'Saccharomyces_cerevisiae_assembly_stats_refseq.txt'.
Collection retrieval finished successfully!


We retrieved the genome assembly and checked the annotation for 'Saccharomyces cerevisiae' (taxid: 559292, strain=S288C) from 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt' (accession: GCF_000146045.2, bioproject: PRJNA128, biosample: NA) using the biomartr R package (Drost and Paszkowski, 2017).
```

> Users can now simply attach the output folder as supplementary data in their study and state in the materials sections. This way, computational and biological reproducibility can be standardized and indeed will become trivial in the context of genome version and annotation version.



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMetaGenomeAnnotations.R
\name{getMetaGenomeAnnotations}
\alias{getMetaGenomeAnnotations}
\title{Retrieve annotation *.gff files for metagenomes from NCBI Genbank}
\usage{
getMetaGenomeAnnotations(
  name,
  path = file.path("_ncbi_downloads", "metagenome", "annotations")
)
}
\arguments{
\item{name}{metagenome name retrieved by \code{\link{listMetaGenomes}}.}

\item{path}{a character string specifying the location (a folder) 
in which the corresponding metagenome annotations shall be stored. 
Default is
\code{path} = \code{file.path("_ncbi_downloads","metagenome","annotations")}.}
}
\description{
Retrieve available annotation *.gff files for metagenomes 
from NCBI Genbank. NCBI Genbank allows users
to download entire metagenomes and their annotations of several metagenome 
projects. This function downloads available metagenomes that can then be 
downloaded via \code{\link{getMetaGenomes}}.
}
\examples{
\dontrun{
# Frist, retrieve a list of available metagenomes
listMetaGenomes()

# Now, retrieve the 'human gut metagenome'
getMetaGenomeAnnotations(name = "human gut metagenome")
}
}
\seealso{
\code{\link{getMetaGenomes}}, \code{\link{listMetaGenomes}}, 
\code{\link{getGFF}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDatasets.R
\name{getDatasets}
\alias{getDatasets}
\title{Retrieve All Available Datasets for a BioMart Database}
\usage{
getDatasets(mart)
}
\arguments{
\item{mart}{a character string specifying the database (mart) for 
which datasets shall be listed.}
}
\description{
This funcion queries the BioMart API and returns a table
storing all available datasets for a selected BioMart databases.
}
\examples{
\dontrun{
# search for available datasets
# getMarts()
# choose database: "ENSEMBL_MART_ENSEMBL"
head(getDatasets("ENSEMBL_MART_ENSEMBL"), 10)
}
}
\seealso{
\code{\link{getMarts}}, \code{\link{getAttributes}}, 
\code{\link{getFilters}}, \code{\link{organismBM}}, 
\code{\link{organismFilters}}, \code{\link{organismAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCDS.R
\name{getCDS}
\alias{getCDS}
\title{Coding Sequence Retrieval}
\usage{
getCDS(
  db = "refseq",
  organism,
  reference = FALSE,
  release = NULL,
  gunzip = FALSE,
  path = file.path("_ncbi_downloads", "CDS")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{there are three options to characterize an organism: 
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning that the most recent database version is used.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{path}{a character string specifying the location (a folder) 
in which the corresponding CDS file shall be stored. 
Default is \code{path} = \code{file.path("_ncbi_downloads","CDS")}.}
}
\value{
File path to downloaded CDS file.
}
\description{
Main retrieval function for coding sequences (CDS) 
of an organism of interest.
By specifying the scientific name of an organism of interest the 
corresponding fasta-file storing the CDS information for the organism 
of interest can be downloaded and stored locally. CDS files can be retrieved 
from several databases.
}
\examples{
\dontrun{
# download the genome of Arabidopsis thaliana from refseq
# and store the corresponding genome CDS file in '_ncbi_downloads/CDS'
file_path <- getCDS( db       = "refseq", 
             organism = "Arabidopsis thaliana", 
             path     = file.path("_ncbi_downloads","CDS"))

Ath_CDS <- read_cds(file_path, format = "fasta")

}
}
\seealso{
\code{\link{getGenome}}, \code{\link{getProteome}}, 
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{getRepeatMasker}}, 
\code{\link{getAssemblyStats}}, \code{\link{meta.retrieval}}, 
\code{\link{read_cds}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomartr-package.R
\docType{package}
\name{biomartr-package}
\alias{biomartr-package}
\alias{biomartr}
\title{Genomic Data Retrieval}
\description{
This package interacts with a suite of web Application 
Programming Interfaces and FTP sites to perform automated genomic data 
retieval and annotation information retrieval.
}
\section{About}{

To automate the retrieval process on a meta-genomic scale, this package
 provides useful interface functions for genomic sequence retrieval and 
 functional annotation retrieval. 
 The major aim of \code{biomartr} is to facilitate computational 
 reproducibility and large-scale handling of genomic data for 
 (meta-)genomic analyses. 
  
In detail, \code{biomartr} aims to provide users with an easy to use 
framework to obtain genome, proteome, CDS, GFF (annotation), genome 
assembly quality, and metagenome project data. Furthermore, an interface to 
the Ensembl Biomart database allows users to retrieve functional annotation 
for genomic loci.
Users can download entire databases
such as 
\itemize{
\item \code{NCBI RefSeq}
\item \code{NCBI nr}
\item \code{NCBI nt}
\item \code{NCBI Genbank}
\item \code{NCBI nt}
\item \code{Ensembl}
\item \code{Ensembl Genomes}
\item \code{UniProt}
}
}

\author{
Hajk-Georg Drost \email{hgd23@cam.ac.uk}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRNASet.R
\name{getRNASet}
\alias{getRNASet}
\title{RNA Retrieval of multiple species}
\usage{
getRNASet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = TRUE,
  gunzip = TRUE,
  update = FALSE,
  path = "set_RNAs"
)
}
\arguments{
\item{db}{a character string specifying the database from which the RNA
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the names of the organisms than shall be retrieved.
There are three available options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a RNA shall be downloaded if it isn't marked
in the database as either a reference RNA or a representative RNA}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{clean_retrieval}{logical value indicating whether or not downloaded files shall be renamed for more convenient downstream data analysis.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{update}{a logical value indicating whether or not files that were already downloaded and are still present in the 
output folder shall be updated and re-loaded (\code{update = TRUE} or whether the existing file shall be retained \code{update = FALSE} (Default)).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding RNAs shall be stored. Default is
\code{path} = \code{"set_RNAs"}.}
}
\value{
File path to downloaded RNAs.
}
\description{
Main RNA retrieval function for a set of organism of interest.
By specifying the scientific names of the organisms of interest the corresponding fasta-files storing the RNA of the organisms of interest
will be downloaded and stored locally. RNA files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory 'set_RNAs' to store
the RNAs of interest as fasta files for future processing.
In case the corresponding fasta file already exists within the
'set_RNAs' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
getRNASet("refseq", organisms = c("Arabidopsis thaliana", 
                                  "Arabidopsis lyrata", 
                                  "Capsella rubella"))
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getRNASet}}, \code{\link{getProteomeSet}}, \code{\link{getGFFSet}},
\code{\link{getCDS}}, \code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_annotation_biomartr.R
\name{check_annotation_biomartr}
\alias{check_annotation_biomartr}
\title{Check whether an annotation file contains outlier lines}
\usage{
check_annotation_biomartr(annotation_file, remove_annotation_outliers = FALSE)
}
\arguments{
\item{annotation_file}{a file path tp the annotation file.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input \code{annotation_file}? 
If yes, then the initial \code{annotation_file} will be overwritten and the removed outlier lines will be stored at \code{\link{tempdir}}
for further exploration.}
}
\description{
Some annotation files include lines with character lengths greater than 65000. This causes problems when trying to import such annotation files into R using \code{\link[rtracklayer]{import}}.
To overcome this issue, this function screens for such lines
in a given annotation file and removes these lines so that
\code{\link[rtracklayer]{import}} can handle the file.
}
\examples{
\dontrun{
# download an example annotation file from NCBI RefSeq
Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana")
# run annotation file check on the downloaded file
biomartr::check_annotation_biomartr(Ath_path)
# several outlier lines were detected, thus we re-run the
# function using 'remove_annotation_outliers = TRUE'
# to remove the outliers and overwrite the file
biomartr::check_annotation_biomartr(Ath_path, remove_annotation_outliers = TRUE)
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refseqOrganisms.R
\name{refseqOrganisms}
\alias{refseqOrganisms}
\title{Retrieve All Organism Names Stored on refseq}
\usage{
refseqOrganisms()
}
\description{
This function extracts all organism names (scientific names)
for which genomes, proteomes, and CDS files are stored on 
the NCBI refseq server.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProteome.R
\name{getProteome}
\alias{getProteome}
\title{Proteome Retrieval}
\usage{
getProteome(
  db = "refseq",
  organism,
  reference = TRUE,
  release = NULL,
  gunzip = FALSE,
  path = file.path("_ncbi_downloads", "proteomes")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
\item \code{db = "uniprot"}
}}

\item{organism}{there are three options to characterize an organism: 
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{path}{a character string specifying the location (a folder) in which 
the corresponding proteome shall be stored. Default is 
\code{path} = \code{file.path("_ncbi_downloads","proteomes")}.}
}
\value{
File path to downloaded proteome.
}
\description{
Main proteome retrieval function for an organism of interest.
By specifying the scientific name of an organism of interest the 
corresponding fasta-file storing the proteome of the organism of interest
can be downloaded and stored locally. Proteome files can be retrieved from 
several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/
 
 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory '_ncbi_downloads/proteomes' to store
the proteome of interest as fasta file for future processing.
}
\examples{
\dontrun{

# download the proteome of Arabidopsis thaliana from refseq
# and store the corresponding proteome file in '_ncbi_downloads/proteomes'
file_path <- getProteome( db       = "refseq", 
             organism = "Arabidopsis thaliana", 
             path     = file.path("_ncbi_downloads","proteomes") )

Ath_proteome <- read_proteome(file_path, format = "fasta")

# download the proteome of Arabidopsis thaliana from genbank
# and store the corresponding proteome file in '_ncbi_downloads/proteomes'
file_path <- getProteome( db       = "genbank", 
             organism = "Arabidopsis thaliana", 
             path     = file.path("_ncbi_downloads","proteomes") )

Ath_proteome <- read_proteome(file_path, format = "fasta")
}
}
\seealso{
\code{\link{getGenome}}, \code{\link{getCDS}}, \code{\link{getGFF}},
\code{\link{getRNA}}, \code{\link{getRepeatMasker}}, 
\code{\link{getAssemblyStats}}, \code{\link{meta.retrieval}}, 
\code{\link{read_proteome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_genome.R
\name{read_genome}
\alias{read_genome}
\title{Import Genome Assembly as Biostrings or data.table object}
\usage{
read_genome(file, format = "fasta", obj.type = "Biostrings", ...)
}
\arguments{
\item{file}{a character string specifying the path to the file 
storing the genome.}

\item{format}{a character string specifying the file format used to store 
the genome, e.g. \code{format = "fasta"} (default) or \code{format = "gbk"}.}

\item{obj.type}{a character string specifying the object stype in which 
the genomic sequence shall be represented. Either as 
\code{obj.type = "Biostrings"} (default) or 
as \code{obj.type = "data.table"}.}

\item{...}{additional arguments that are used by the 
\code{\link[seqinr]{read.fasta}} function.}
}
\value{
Either a \code{Biostrings} or \code{data.table} object.
}
\description{
This function reads an organism specific genome stored in a 
defined file format.
}
\details{
This function takes a string specifying the path to the genome file
of interest as first argument (e.g. the path returned by 
\code{\link{getGenome}}).
}
\seealso{
\code{\link{getGenome}}, \code{\link{read_proteome}}, 
\code{\link{read_cds}}, \code{\link{read_gff}}, \code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listDatabases.R, R/listNCBIDatabases.R
\name{listDatabases}
\alias{listDatabases}
\alias{listNCBIDatabases}
\title{Retrieve a List of Available NCBI Databases for Download}
\usage{
listDatabases(db = "nr", update = FALSE)

listNCBIDatabases(db = "nr", update = FALSE)
}
\arguments{
\item{db}{a character string specifying the name of the database that shall 
be searched for.}

\item{update}{a logical value specifying whether or not the local 
listDatabases.txt file shall be updated by remote access to NCBI.}
}
\description{
This function allows you to retrieve a list of database names 
and versions that can be downloaded from correspondning servers.

Database retrieval is crucial for most biological studies and analyses.
There is a vast diversity of databases that can be accessed remotely or that 
can be downloaded to your local machine. This function provides an interface 
to databases that can be downloaded from NCBI servers and lists all available
databases and their database version to be able to select an appropriate 
database for download with \code{\link{download.database}}.
}
\examples{
\dontrun{
# retrieve all versions of the NCBI 'nr' database that can be downloaded
listNCBIDatabases(db = "nr")

# analogous:
# listNCBIDatabases(db = "cdd")
# listNCBIDatabases(db = "nt")
# listNCBIDatabases(db = "gss")
# listNCBIDatabases(db = "refseq_protein")
}
}
\seealso{
\code{\link{download.database}}, \code{\link{download.database.all}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMetaGenomes.R
\name{getMetaGenomes}
\alias{getMetaGenomes}
\title{Retrieve metagenomes from NCBI Genbank}
\usage{
getMetaGenomes(name, path = file.path("_ncbi_downloads", "metagenome"))
}
\arguments{
\item{name}{metagenome name retrieved by \code{\link{listMetaGenomes}}.}

\item{path}{a character string specifying the location (a folder) in 
which the corresponding metagenome shall be stored. 
Default is \code{path} = \code{file.path("_ncbi_downloads","metagenome")}.}
}
\description{
Retrieve available metagenomes from NCBI Genbank. 
NCBI Genbank allows users to download entire metagenomes of several 
metagenome projects. This function downloads available metagenomes that can 
then be downloaded via \code{\link{getMetaGenomes}}.
}
\examples{
\dontrun{
# Frist, retrieve a list of available metagenomes
listMetaGenomes()

# Now, retrieve the 'human gut metagenome'
getMetaGenomes(name = "human gut metagenome")
}
}
\seealso{
\code{\link{getMetaGenomeAnnotations}}, 
\code{\link{listMetaGenomes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_cds.R
\name{read_cds}
\alias{read_cds}
\title{Import CDS as Biostrings or data.table object}
\usage{
read_cds(
  file,
  format = "fasta",
  obj.type = "Biostrings",
  delete_corrupt = FALSE,
  ...
)
}
\arguments{
\item{file}{a character string specifying the path to the file storing 
the CDS.}

\item{format}{a character string specifying the file format used to store 
the genome, e.g. \code{format = "fasta"} (default) or \code{format = "gbk"}.}

\item{obj.type}{a character string specifying the object stype in which the 
genomic sequence shall be represented. 
Either as \code{obj.type = "Biostrings"} (default) or as 
\code{obj.type = "data.table"}.}

\item{delete_corrupt}{a logical value specifying whether potential CDS 
sequences that cannot be divided by 3 shall be
be excluded from the the dataset. Default is \code{delete_corrupt = FALSE}.}

\item{...}{additional arguments that are used by
\code{\link[seqinr]{read.fasta}}.}
}
\value{
A data.table storing the gene id in the first column and the 
corresponding sequence as string in the second column.
}
\description{
This function reads an organism specific CDS stored in a 
defined file format.
}
\details{
The \code{read.cds} function takes a string specifying the path 
to the cds file of interest as first argument.

It is possible to read in different proteome file standards such as 
\emph{fasta} or \emph{genebank}.

CDS stored in fasta files can be downloaded from 
http://www.ensembl.org/info/data/ftp/index.html.
}
\seealso{
\code{\link{getCDS}}, \code{\link{read_genome}}, 
\code{\link{read_proteome}}, \code{\link{read_gff}}, \code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.retrieval.all.R
\name{meta.retrieval.all}
\alias{meta.retrieval.all}
\title{Perform Meta-Genome Retrieval of all organisms in all kingdoms of life}
\usage{
meta.retrieval.all(db = "refseq", type = "genome", reference = FALSE)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved: 
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"} 
\item \code{db = "emsembl"}
\item \code{db = "ensemblgenomes"}
}}

\item{type}{type of sequences that shall be retrieved. Options are:
\itemize{
 \item \code{type = "genome"} :
 for genome assembly retrieval; see also \code{\link{getGenome}}), 
 \item \code{type = "proteome"} :
 (for proteome retrieval; see also \code{\link{getProteome}}),
 \item \code{type = "cds"} :
 (for coding sequence retrieval; see also \code{\link{getCDS}}),
 \item \code{type = "gff"} :
(for annotation file retrieval in gff format; see also \code{\link{getGFF}}),
\item \code{type = "gtf"} :
(for annotation file retrieval in gtf format 
(only for ensembl and ensemblgenomes); see also \code{\link{getGTF}}),
 \item \code{type = "rna"} :
 (for RNA file retrieval in fasta format; see also \code{\link{getRNA}}),
 \item \code{type = "rm"} :
 (for Repeat Masker output file retrieval; see also 
 \code{\link{getRepeatMasker}}),
 \item \code{type = "assemblystats"} (for genome assembly quality stats 
 file retrieval; see also \code{\link{getAssemblyStats}}).
 }}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database 
as either a reference genome or a representative genome. Options are:
\itemize{
\item \code{reference = FALSE} (Default): all organisms (reference, representative, and non-representative genomes) are downloaded.
\item \code{reference = TRUE}: organisms that are downloaded must be either a reference or representative genome. Thus, most genomes which are usually non-reference genomes
will not be downloaded.
}}
}
\value{
a character vector storing the file paths of the retrieved files.
}
\description{
Download genomes, proteomes, cds, gff, rna, or assembly stats 
files of individual species of all kingdoms of life.
}
\details{
This function aims to perform bulk retrieval of all genomes 
of species for all kingdoms of life.
}
\examples{
\dontrun{
# download all genomes from refseq
meta.retrieval.all(db = "refseq", type = "genome")
# download all vertebrate genomes from genbank
meta.retrieval.all(db = "genbank", type = "genome")
# download all vertebrate genomes from ensemblgenomes
meta.retrieval.all(db = "genbank", type = "ensemblgenomes")
}
}
\seealso{
\code{\link{meta.retrieval}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getENSEMBLInfo.R
\name{getENSEMBLInfo}
\alias{getENSEMBLInfo}
\title{Retrieve ENSEMBL info file}
\usage{
getENSEMBLInfo()
}
\description{
Retrieve species and genome information from 
http://rest.ensembl.org/info/species?content-type=application/json/.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getReleases.R
\name{getReleases}
\alias{getReleases}
\title{Retrieve available database releases or versions 
of ENSEMBL}
\usage{
getReleases(db = "ensembl")
}
\arguments{
\item{db}{a character string specifying the database from which 
available resease versions shall be retrieved:
\itemize{
\item \code{db = "ensembl"}
}}
}
\description{
Retrieve available database releases or versions of 
ENSEMBL.
}
\examples{
\dontrun{ 
# retrieve available resease versions of ENSEMBL
getReleases("ensembl")
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_rm.R
\name{read_rm}
\alias{read_rm}
\title{Import Repeat Masker output file}
\usage{
read_rm(file)
}
\arguments{
\item{file}{a character string specifying the path to the file 
storing the Repeat Masker output (e.g. retrieved with \code{\link{getRepeatMasker}}).}
}
\description{
This function reads an organism specific 
Repeat Masker output file.
}
\details{
This function takes a string specifying the path to the 
Repeat Masker output file of interest as first argument.
}
\seealso{
\code{\link{getRepeatMasker}}, \code{\link{read_genome}}, 
\code{\link{read_proteome}}, \code{\link{read_gff}}, \code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGroups.R
\name{getGroups}
\alias{getGroups}
\title{Retrieve available groups for a kingdom of life (only available for NCBI RefSeq and NCBI Genbank)}
\usage{
getGroups(db = "refseq", kingdom)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved: 
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
}

Default is \code{db = "refseq"}.}

\item{kingdom}{a character string specifying for which kingdom of life 
groups shall be retrieved. See \code{\link{getKingdoms}} for details.}
}
\description{
A short list of available groups for a kingdom of life.
}
\examples{
# get possible kigdom names
getKingdoms(db = "refseq")
\dontrun{
# retrieve subgroups for vertebrate_mammalian available from refseq
getGroups(db = "refseq", kingdom = "vertebrate_mammalian")

# get possible kigdom names
getKingdoms(db = "genbank")
# retrieve subgroups for vertebrate_mammalian available from genbank
getGroups(db = "genbank", kingdom = "vertebrate_mammalian")
}
}
\seealso{
\code{\link{meta.retrieval}}, \code{\link{getGenome}}, 
\code{\link{getProteome}}, \code{\link{getCDS}}, \code{\link{getKingdoms}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRepeatMasker.R
\name{getRepeatMasker}
\alias{getRepeatMasker}
\title{Repeat Masker Retrieval}
\usage{
getRepeatMasker(
  db = "refseq",
  organism,
  reference = FALSE,
  path = file.path("_ncbi_downloads", "repeatmasker")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
}}

\item{organism}{a character string specifying the scientific name of the 
organism of interest, e.g. \code{organism = "Homo sapiens"}.}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{path}{a character string specifying the location (a folder) in which 
the corresponding file shall be stored. Default is 
\code{path} = \code{file.path("_ncbi_downloads","repeatmasker")}.}
}
\value{
File path to downloaded Repeat Masker output file.
}
\description{
Main Repeat Masker output retrieval function for an 
organism of interest.
By specifying the scientific name of an organism of interest the 
corresponding Repeat Masker file storing the genome of the organism of 
interest can be downloaded and stored locally. 
Repeat Masker files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
 
and creates a directory '_ncbi_downloads/repeatmasker' to store
the files of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/repeatmasker' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{

# download the Repeat Masker output file of Arabidopsis thaliana from refseq
# and store the corresponding genome file in '_ncbi_downloads/genomes'
file_path <- getRepeatMasker( db       = "refseq", 
             organism = "Arabidopsis thaliana", 
             path = file.path("_ncbi_downloads","repeatmasker"))

Ath_repeatmasker <- read_rm(file_path)


# download the Repeat Masker output file of Arabidopsis thaliana from genbank
# and store the corresponding genome file in '_ncbi_downloads/genomes'
file_path <- getRepeatMasker( db       = "genbank", 
             organism = "Arabidopsis thaliana", 
             path = file.path("_ncbi_downloads","repeatmasker"))

Ath_repeatmasker <- read_rm(file_path)
}

}
\seealso{
\code{\link{getProteome}}, \code{\link{getCDS}}, 
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}}, 
\code{\link{read_rm}}, \code{\link{getGenome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_assemblystats.R
\name{read_assemblystats}
\alias{read_assemblystats}
\title{Import Genome Assembly Stats File}
\usage{
read_assemblystats(file, type = "raw")
}
\arguments{
\item{file}{a character string specifying the path to the file storing 
the Genome Assembly Stats file.}

\item{type}{either \code{type = "raw"} to import the entire genome assembly 
stats file or \code{type = "stats"} to import overall statistics including 
all chromosomes, mitochondria and plastids.}
}
\description{
This function reads an organism specific Genome Assembly 
Stats file that was retrieved with \code{\link{getAssemblyStats}}.
}
\details{
This function takes a string specifying the path to the Genome 
Assembly Stats file of interest (e.g. the path returned by 
\code{\link{getAssemblyStats}}) and imports it.
}
\seealso{
\code{\link{getAssemblyStats}}, \code{\link{read_genome}}, 
\code{\link{read_proteome}}, \code{\link{read_cds}}, \code{\link{read_gff}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getKingdoms.R
\name{getKingdoms}
\alias{getKingdoms}
\title{Retrieve available kingdoms of life}
\usage{
getKingdoms(db = "refseq")
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved: \code{db = "refseq"}, \code{db = "genbank"}, 
\code{db = "ensembl"}, \code{db = "ensemblgenomes"}. 
Default is \code{db = "refseq"}.}
}
\description{
A short list of available kingdoms of life
}
\examples{
# retrieve kingdoms available from refseq
getKingdoms(db = "refseq")

# retrieve kingdoms available from genbank
getKingdoms(db = "genbank")
}
\seealso{
\code{\link{meta.retrieval}}, \code{\link{getGenome}}, 
\code{\link{getProteome}}, \code{\link{getCDS}}, \code{\link{getGroups}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCollectionSet.R
\name{getCollectionSet}
\alias{getCollectionSet}
\title{Retrieve a Collection: Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats of multiple species}
\usage{
getCollectionSet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = FALSE,
  gunzip = TRUE,
  update = FALSE,
  remove_annotation_outliers = TRUE,
  path = "set_collections"
)
}
\arguments{
\item{db}{a character string specifying the database from which the collection
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the scientific names of the organisms for
which collections shall be retrieved. There are three options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a collection shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning that the most recent database version is used.}

\item{clean_retrieval, }{a logical, default FALSE. Cleaning file names for more convenient
downstream processing.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{update}{a logical, default FALSE. The existing file will be retained if existing.
If TRUE, will download and overwrite the file.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input annotation_file?
If yes, then the initial annotation_file will be overwritten and the removed outlier lines
will be stored at \code{\link{tempdir}} for further exploration.}

\item{path}{a character string specifying the location (a folder) in which
the corresponding collection shall be stored. Default is
\code{path} = \code{file.path("_db_downloads","collections")}.}
}
\value{
File path to downloaded genome.
}
\description{
Main collection retrieval function for an organism of interest.
By specifying the scientific name of an organism of interest a collection consisting of
the genome file, proteome file, CDS file, RNA file, GFF file, Repeat Masker file, AssemblyStats
file of the organism of interest
can be downloaded and stored locally. Collections can be retrieved from
several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory '_ncbi_downloads/collection' to store
the genome of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/collection' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
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
}

}
\seealso{
\code{\link{getCollection}}, \code{\link{getGenomeSet}}, \code{\link{getProteomeSet}}, \code{\link{getCDSSet}},
\code{\link{getGenome}}, \code{\link{getProteome}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenome.R
\name{getGenome}
\alias{getGenome}
\title{Genome Retrieval}
\usage{
getGenome(
  db = "refseq",
  organism,
  reference = FALSE,
  release = NULL,
  gunzip = FALSE,
  path = file.path("_ncbi_downloads", "genomes"),
  assembly_type = "toplevel",
  kingdom_assembly_summary_file = NULL
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{there are three options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{a numeric, the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used. \code{release = 75} would for human would give the stable
GRCh37 release in ensembl. Value must be > 46, since ensembl did not structure their data
if the standard format before that.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{path}{a character string specifying the location (a folder) in which
the corresponding genome shall be stored. Default is
\code{path} = \code{file.path("_ncbi_downloads","genomes")}.}

\item{assembly_type}{a character string specifying from which assembly type the genome
shall be retrieved from (ensembl only, else this argument is ignored):
Default is
\code{assembly_type = "toplevel")}.
This will give you all multi-chromosomes (copies of the same chromosome with small variations).
As an example the toplevel fasta genome in human is over 70 GB uncompressed.
To get primary assembly with 1 chromosome variant per chromosome:
\code{assembly_type = "primary_assembly")}.
As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed:}
}
\value{
File path to downloaded genome.
}
\description{
Main genome retrieval function for an organism of interest.
By specifying the scientific name of an organism of interest the
corresponding fasta-file storing the genome of the organism of interest
can be downloaded and stored locally. Genome files can be retrieved from
several databases. In addition, the genome summary statistics for the
retrieved species is stored locally to provide users with
insights regarding the genome assembly quality (see \code{\link{summary_genome}} for details).
This is useful when comparing genomes with large difference in genome assembly qualities.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory '_ncbi_downloads/genomes' to store
the genome of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/genomes' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{

# download the genome of Arabidopsis thaliana from refseq
# and store the corresponding genome file in '_ncbi_downloads/genomes'
file_path <- getGenome( db       = "refseq",
             organism = "Arabidopsis thaliana",
             path = file.path("_ncbi_downloads","genomes"))

Ath_genome <- read_genome(file_path, format = "fasta")


# download the genome of Arabidopsis thaliana from genbank
# and store the corresponding genome file in '_ncbi_downloads/genomes'
file_path <- getGenome( db       = "genbank",
             organism = "Arabidopsis thaliana",
             path = file.path("_ncbi_downloads","genomes"))

Ath_genome <- read_genome(file_path, format = "fasta")
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getProteome}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{getRepeatMasker}},
\code{\link{getAssemblyStats}}, \code{\link{summary_genome}},
\code{\link{meta.retrieval}}, \code{\link{meta.retrieval.all}}, \code{\link{read_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGFFSet.R
\name{getGFFSet}
\alias{getGFFSet}
\title{GFF retrieval of multiple species}
\usage{
getGFFSet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = TRUE,
  gunzip = TRUE,
  remove_annotation_outliers = FALSE,
  update = FALSE,
  path = "set_GFF"
)
}
\arguments{
\item{db}{a character string specifying the database from which the GFF
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the names of the organisms than shall be retrieved.
There are three available options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a GFF shall be downloaded if it isn't marked
in the database as either a reference GFF or a representative GFF}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{clean_retrieval}{logical value indicating whether or not downloaded files shall be renamed for more convenient downstream data analysis.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input \code{annotation_file}? 
If yes, then the initial \code{annotation_file} will be overwritten and the removed outlier lines will be stored at \code{\link{tempdir}}
for further exploration.}

\item{update}{a logical value indicating whether or not files that were already downloaded and are still present in the 
output folder shall be updated and re-loaded (\code{update = TRUE} or whether the existing file shall be retained \code{update = FALSE} (Default)).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding CDSs shall be stored. Default is
\code{path} = \code{"set_CDS"}.}
}
\value{
File path to downloaded CDSs.
}
\description{
Main GFF retrieval function for a set of organism of interest.
By specifying the scientific names of the organisms of interest the corresponding fasta-files storing the GFF of the organisms of interest
will be downloaded and stored locally. GFF files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory 'set_CDSs' to store
the CDSs of interest as fasta files for future processing.
In case the corresponding fasta file already exists within the
'set_CDSs' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
getGFFSet("refseq", organisms = c("Arabidopsis thaliana", 
                                  "Arabidopsis lyrata", 
                                  "Capsella rubella"))
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getProteomeSet}}, \code{\link{getCDSSet}}, \code{\link{getRNASet}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_cds}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGFF.R
\name{getGFF}
\alias{getGFF}
\title{Genome Annotation Retrieval (GFF3)}
\usage{
getGFF(
  db = "refseq",
  organism,
  reference = FALSE,
  release = NULL,
  gunzip = FALSE,
  remove_annotation_outliers = FALSE,
  path = file.path("_ncbi_downloads", "annotation")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{a character string specifying the scientific name of the 
organism of interest, e.g. \code{organism = "Homo sapiens"}.}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input \code{annotation_file}? 
If yes, then the initial \code{annotation_file} will be overwritten and the removed outlier lines will be stored at \code{\link{tempdir}}
for further exploration.}

\item{path}{a character string specifying the location (a folder) in which 
the corresponding annotation file shall be stored. Default is 
\code{path = file.path("_ncbi_downloads","genomes")}.}
}
\value{
File path to downloaded annotation file.
}
\description{
Main retrieval function for GFF files of an 
organism of interest. By specifying the scientific name of an organism of 
interest the corresponding gff file storing the annotation  for the organism 
of interest can be downloaded and stored locally. GFF files can be retrieved 
from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
 
and creates a directory '_ncbi_downloads/annotation' to store
the genome of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/annotation' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
# download the annotation of Arabidopsis thaliana from refseq
# and store the corresponding genome file in '_ncbi_downloads/annotation'
getGFF( db       = "refseq", 
               organism = "Arabidopsis thaliana", 
               path = file.path("_ncbi_downloads","annotation"))


# download the genome of Arabidopsis thaliana from genbank
# and store the corresponding genome file in '_ncbi_downloads/annotation'
getGFF( db       = "genbank", 
               organism = "Arabidopsis thaliana", 
               path = file.path("_ncbi_downloads","annotation"))

}

}
\seealso{
\code{\link{getProteome}}, \code{\link{getCDS}}, 
\code{\link{getGenome}}, \code{\link{getRNA}}, \code{\link{getRepeatMasker}}, 
\code{\link{getAssemblyStats}}, \code{\link{meta.retrieval}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.genome.available.R
\name{is.genome.available}
\alias{is.genome.available}
\title{Check Genome Availability}
\usage{
is.genome.available(db = "refseq", organism, details = FALSE)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
\item \code{db = "uniprot"}
}}

\item{organism}{there are three options to characterize an organism: 
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{details}{a logical value specifying whether or not details on genome 
size, kingdom, etc. shall be printed to the console intead of a 
boolean value.}
}
\value{
a logical value specifing whether or not the genome of the input 
organism is available. In case \code{details} = \code{TRUE} only a character 
string specifying the genome details is being returned.
}
\description{
This function checks the availability of a given genome on the 
NBCI servers specified as scientific name.
}
\details{
Internally this function calls the \code{\link{listGenomes}} function to 
detect all available genomes and checks whether or not the specified organism
is available for download.
}
\examples{
\dontrun{
# checking whether the Homo sapiens genome is stored on NCBI
is.genome.available(organism = "Homo sapiens", db = "refseq")

# and printing details
is.genome.available(organism = "Homo sapiens", db = "refseq", details = TRUE)

# checking whether the Homo sapiens genome is stored on ENSEMBL
is.genome.available(organism = "Homo sapiens", db = "ensembl")

# and printing details
is.genome.available(organism = "Homo sapiens",
                    details = TRUE, 
                    db = "ensembl")
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getKingdomAssemblySummary.R
\name{getKingdomAssemblySummary}
\alias{getKingdomAssemblySummary}
\title{Retrieve and summarise the assembly_summary.txt files from 
NCBI for all kingdoms}
\usage{
getKingdomAssemblySummary(db)
}
\arguments{
\item{db}{database name. E.g. \code{refseq} or \code{genbank}.}
}
\description{
Retrieval function of the assembly_summary.txt file 
from NCBI for all kingdoms.
The assembly_summary.txt files store available species on NCBI.
}
\examples{
\dontrun{
test <- getKingdomAssemblySummary(db = "refseq")
test
}
}
\seealso{
\code{\link{getSummaryFile}}, \code{\link{getMetaGenomeSummary}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFilters.R
\name{getFilters}
\alias{getFilters}
\title{Retrieve All Available Filters for a Specific Dataset}
\usage{
getFilters(mart, dataset)
}
\arguments{
\item{mart}{a character string specifying the database (mart) for which 
datasets shall be listed.}

\item{dataset}{a character string specifying the dataset for which filters 
shall be listed.}
}
\description{
This funcion queries the BioMart API and returns a table
storing all available filters for a specific dataset.
}
\examples{
\dontrun{
# search for available datasets
# getMarts()
# choose database (mart): "ENSEMBL_MART_ENSEMBL"
# head(getDatasets(mart = "ENSEMBL_MART_ENSEMBL"), 10)
# choose dataset: "hsapiens_gene_ensembl"
head(getFilters(mart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl") , 5)
}
}
\seealso{
\code{\link{getMarts}}, \code{\link{getDatasets}}, 
\code{\link{getAttributes}}, \code{\link{organismBM}}, 
\code{\link{organismFilters}}, \code{\link{organismAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_proteome.R
\name{read_proteome}
\alias{read_proteome}
\title{Import Proteome as Biostrings or data.table object}
\usage{
read_proteome(file, format = "fasta", obj.type = "Biostrings", ...)
}
\arguments{
\item{file}{a character string specifying the path to the file storing 
the proteome.}

\item{format}{a character string specifying the file format used to store the 
genome, e.g. \code{format = "fasta"} (default) or \code{format = "gbk"}.}

\item{obj.type}{a character string specifying the object stype in which the 
genomic sequence shall be represented. 
Either as \code{obj.type = "Biostrings"} (default) or as 
\code{obj.type = "data.table"}.}

\item{...}{additional arguments that are used by
\code{\link[seqinr]{read.fasta}}.}
}
\value{
Either a \code{Biostrings} or \code{data.table} object.
}
\description{
This function reads an organism specific proteome stored in a 
defined file format.
}
\details{
This function takes a string specifying the path to the 
proteome file of interest as first argument.

It is possible to read in different proteome file standards such as 
\emph{fasta} or \emph{genebank}.
}
\seealso{
\code{\link{getProteome}}, \code{\link{read_genome}}, 
\code{\link{read_gff}}, \code{\link{read_cds}}, \code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAssemblyStats.R
\name{getAssemblyStats}
\alias{getAssemblyStats}
\title{Genome Assembly Stats Retrieval}
\usage{
getAssemblyStats(
  db = "refseq",
  organism,
  reference = FALSE,
  type = "download",
  path = file.path("_ncbi_downloads", "genomeassembly_stats")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{a character string specifying the scientific name of the 
organism of interest, e.g. \code{organism = "Homo sapiens"}.}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{type}{shall only the file be retrieved (default) 
\code{type = "download"} or should the corresponding file be downloaded and 
subsequently be imported \code{type = "import"}.}

\item{path}{a character string specifying the location (a folder) in
which the corresponding file shall be stored. Default is 
\code{path} = \code{file.path("_ncbi_downloads","genomeassembly_stats")}.}
}
\value{
File path to downloaded genome assembly stats file.
}
\description{
Main genome assembly stats retrieval function for an organism
of interest. By specifying the scientific name of an organism of interest the
corresponding  genome assembly stats file storing the assembly statistics of 
the organism of interest can be downloaded and stored locally. 
Genome assembly stats files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
 
to retrieve available scientific names of organisms and creates a directory
'_ncbi_downloads/genomeassembly_stats' to store
the Genome Assembly Stats of interest as text file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/genomeassembly_stats' folder and is
accessible within the workspace, no download process will be performed.

An example genome assembly stats file can be found here:
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/
GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_stats.txt.
}
\examples{
\dontrun{
# download the genome assembly stats file of Saccharomyces cerevisiae
# from NCBI RefSeq
# and store the corresponding genome file in 
# '_ncbi_downloads/genomeassembly_stats'
file_path <- getAssemblyStats( db = "refseq", 
                 organism = "Saccharomyces cerevisiae", 
                 path = file.path("_ncbi_downloads","genomeassembly_stats"))
# import the raw file as it is downloaded
Scerevisiae.stats <- read_assemblystats(file_path, type = "raw")

# download the genome assembly stats file of Saccharomyces cerevisiae
# from NCBI RefSeq 
# and import overall statistics of the genome assembly
Scerevisiae.stats.import <- getAssemblyStats( db = "refseq", 
                 organism = "Saccharomyces cerevisiae",
                 type = "import", 
                 path = file.path("_ncbi_downloads","genomeassembly_stats"))
}

}
\seealso{
\code{\link{getGenome}}, \code{\link{getProteome}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}}, 
\code{\link{read_assemblystats}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listMetaGenomes.R
\name{listMetaGenomes}
\alias{listMetaGenomes}
\title{List available metagenomes on NCBI Genbank}
\usage{
listMetaGenomes(details = FALSE)
}
\arguments{
\item{details}{a boolean value specifying whether only the scientific names 
of stored metagenomes shall be returned (\code{details = FALSE}) or all 
information such as "organism_name","bioproject", 
etc (\code{details = TRUE}).}
}
\description{
List available metagenomes on NCBI genbank. NCBI genbank 
allows users to download entire metagenomes of several metagenome projects. 
This function lists all available metagenomes that can then be downloaded via
\code{\link{getMetaGenomes}}.
}
\examples{
\dontrun{
# retrieve available metagenome projects at NCBI Genbank
listMetaGenomes()
# retrieve detailed information on available metagenome projects 
# at NCBI Genbank
listMetaGenomes(details = TRUE)
}
}
\seealso{
\code{\link{getMetaGenomes}}, \code{\link{getMetaGenomeSummary}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart.R
\name{biomart}
\alias{biomart}
\title{Main BioMart Query Function}
\usage{
biomart(genes, mart, dataset, attributes, filters, ...)
}
\arguments{
\item{genes}{a character vector storing the gene ids of a organisms
of interest to be queried against BioMart.}

\item{mart}{a character string specifying the mart to be used. Users
can obtain available marts using \code{\link{getMarts}}.}

\item{dataset}{a character string specifying the dataset within the mart to
be used, e.g. \code{dataset} = \code{"hsapiens_gene_ensembl"}.}

\item{attributes}{a character vector specifying the attributes that shall be
used, e.g. \code{attributes} = 
\code{c("start_position","end_position","description")}.}

\item{filters}{a character vector specifying the filter (query key) for the 
BioMart query, e.g. \code{filter} = \code{"ensembl_gene_id"}.}

\item{...}{additional parameters for the
\code{\link[biomaRt]{getBM}} function.}
}
\value{
A data.table storing the initial query gene vector in 
the first column, the output
gene vector in the second column, and all attributes in
the following columns.
}
\description{
This function takes a set of gene ids and the biomart
specifications and performs a biomart query for the given set of gene ids.
}
\details{
This function is the main query function of the biomartr package.

It enables to fastly access annotations of a given gene set based 
on the \pkg{biomaRt} package
implemented by Steffen Durinck et al.
}
\examples{
\dontrun{
# 1) select a mart
getMarts()

# we will select mart 'plants_mart' and search for available datasets
getDatasets(mart = "plants_mart")

# we choose dataset 'athaliana_eg_gene' and run biomart()
# using mart: 'plants_mart', dataset: "athaliana_eg_gene"
# attributes: c("start_position","end_position","description")
# for an example gene set of Arabidopsis thaliana:
# c("AT1G06090", "AT1G06100", "AT1G06110", "AT1G06120",
#    "AT1G06130", "AT1G06200")

biomart(genes      = c("AT1G06090", "AT1G06100",
                       "AT1G06110", "AT1G06120",
                       "AT1G06130", "AT1G06200"),
        mart       = "plants_mart",
        dataset    = "athaliana_eg_gene",
        attributes = c("start_position","end_position","description"),
        filters    = "ensembl_gene_id")
}

}
\seealso{
\code{\link{organismFilters}}, \code{\link{organismBM}},
\code{\link[biomaRt]{listAttributes}}, \code{\link[biomaRt]{getBM}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organismFilters.R
\name{organismFilters}
\alias{organismFilters}
\title{Retrieve Ensembl Biomart filters for a qyery organism}
\usage{
organismFilters(organism, update = FALSE, topic = NULL)
}
\arguments{
\item{organism}{a character string specifying the scientific name of a 
query organism.}

\item{update}{a logical value specifying whether or not the local 
listMart.txt, listDatasets.txt, and listFilters_organism.txt files shall be 
updated by remote access to BioMart.}

\item{topic}{a character string specifying a topic (category) of filters, 
e.g. \code{topic} = \code{"id"}.}
}
\value{
a data.frame storing corresponding filter names, description, 
datasets, and marts.
}
\description{
In addition to the \code{\link{organismBM}} and 
\code{\link{organismAttributes}} functions, this function
returns all available filters that can be accessed through different marts 
and datasets for a given query organism.
}
\details{
For a given query organism, this function retrieves all available 
filters that can be accessed through different marts and datasets.

Sometimes the same filter names correspond to different datasets and 
marts causing problems when using \code{\link{getMarts}}. 
The approach introduced by this function provides (again) a organism centric 
way of accessing organism specific filters.

The \code{topic} argument allows the user to search for specific filters 
topics/categories for faster selection.
}
\note{
When you run this function for the first time, the data retrieval procedure 
will take some time, due to the remote access to BioMart. The corresponding 
result is then saved in a *.txt file within the \code{\link{tempdir}}
directory  named "_biomart/listMarts.txt","_biomart/listDatasets.txt", and 
"_biomart/listFilters_organism.txt", allowing subsequent queries to perform 
much faster.
}
\examples{
\dontrun{ 
# search for filter topic "id" 
head(organismFilters("Homo sapiens", topic = "id"), 20)
}
}
\references{
\url{http://biomart.org/}

Mapping identifiers for the integration of genomic datasets with the
R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan
Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

BioMart and Bioconductor: a powerful link between biological databases and
microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean
Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21,
3439-3440 (2005).
}
\seealso{
\code{\link{organismBM}}, \code{\link{organismAttributes}}, 
\code{\link{getAttributes}}, \code{\link{getDatasets}}, 
\code{\link{getMarts}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getENSEMBL.Seq.R
\name{getENSEMBL.Seq}
\alias{getENSEMBL.Seq}
\title{Helper function for retrieving biological sequence files from ENSEMBL}
\usage{
getENSEMBL.Seq(
  organism,
  type = "dna",
  id.type = "toplevel",
  release = NULL,
  path
)
}
\arguments{
\item{organism}{scientific name of the organism of interest.}

\item{type}{biological sequence type.}

\item{id.type}{a character, default "toplevel". id type of assembly, either toplevel or primary_assembly usually.}

\item{release}{a numeric, the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used. \code{release = 75} would for human would give the stable
GRCh37 release in ensembl. Value must be > 46, since ensembl did not structure their data
if the standard format before that.}

\item{path}{location where file shall be stored.}
}
\value{
either a character path to downloaded file, or a logical FALSE, specifying failure.
}
\description{
This function downloads gff files of query
organisms from ENSEMBL.
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organismAttributes.R
\name{organismAttributes}
\alias{organismAttributes}
\title{Retrieve Ensembl Biomart attributes for a query organism}
\usage{
organismAttributes(organism, update = FALSE, topic = NULL)
}
\arguments{
\item{organism}{a character string specifying the scientific name 
of a query organism.}

\item{update}{a logical value specifying whether or not the local 
listMart.txt, listDatasets.txt, and listAttributes_organism.txt files shall 
be updated by remote access to BioMart.}

\item{topic}{a character string specifying a topic (category) of attributes, 
e.g. \code{topic} = \code{"id"}.}
}
\value{
a data.frame storing corresponding attribute names, description, 
datasets, and marts.
}
\description{
In addition to the \code{\link{organismBM}} function, 
this function returns all available attributes that can be accessed through 
different marts and datasets for a given query organism.
}
\details{
For a given query organism, this function retrieves all available attributes 
that can be accessed through different marts and datasets.

Sometimes the same attribute names correspond to different datasets and 
marts causing problems when using \code{\link{getMarts}}. The approach 
introduced by this function provides (again) a organism centric way of 
accessing organism specific attributes.

The \code{topic} argument allows the user to search for specific attribute 
topics/categories for faster filtering.
}
\note{
When you run this function for the first time, the data retrieval procedure 
will take some time, due to the remote access to BioMart. The corresponding 
result is then saved in a *.txt file within the \code{\link{tempdir}}
directory named "_biomart/listMarts.txt","_biomart/listDatasets.txt", and 
"_biomart/listAttributes_organism.txt", allowing subsequent queries to 
perform much faster.
}
\examples{
\dontrun{ 
# search for attribute topic id
head(organismAttributes("Homo sapiens", topic = "id"), 20)
}
}
\references{
\url{http://biomart.org/}

Mapping identifiers for the integration of genomic datasets with the
R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan
Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

BioMart and Bioconductor: a powerful link between biological databases and
microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean
Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21,
3439-3440 (2005).
}
\seealso{
\code{\link{organismFilters}}, \code{\link{organismBM}}, 
\code{\link{biomart}}, \code{\link[biomaRt]{listAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRNA.R
\name{getRNA}
\alias{getRNA}
\title{RNA Sequence Retrieval}
\usage{
getRNA(
  db = "refseq",
  organism,
  reference = FALSE,
  release = NULL,
  path = file.path("_ncbi_downloads", "RNA")
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved: 
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"} 
}}

\item{organism}{there are three options to characterize an organism: 
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{path}{a character string specifying the location (a folder) in which 
the corresponding
CDS file shall be stored. Default is 
\code{path} = \code{file.path("_ncbi_downloads","RNA")}.}
}
\value{
File path to downloaded RNA file.
}
\description{
Main retrieval function for RNA sequences of an organism 
of interest. By specifying the scientific name of an organism of interest the
corresponding fasta-file storing the RNA information for the organism 
of interest can be downloaded and stored locally. 
RNA files can be retrieved from several databases.
}
\examples{
\dontrun{
# download the RNA of Arabidopsis thaliana from refseq
# and store the corresponding RNA file in '_ncbi_downloads/RNA'
file_path <- getRNA( db       = "refseq", 
             organism = "Arabidopsis thaliana", 
             path     = file.path("_ncbi_downloads","RNA"))

Ath_RNA <- read_rna(file_path, format = "fasta")
}
}
\seealso{
\code{\link{getGenome}}, \code{\link{getProteome}}, 
\code{\link{getGTF}}, \code{\link{getGFF}}, \code{\link{getRepeatMasker}}, 
\code{\link{getAssemblyStats}}, \code{\link{meta.retrieval}}, 
\code{\link{read_cds}}, \code{\link{getCDS}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.retrieval.R
\name{clean.retrieval}
\alias{clean.retrieval}
\title{Format \code{\link{meta.retrieval}} output}
\usage{
clean.retrieval(x, gunzip = TRUE)
}
\arguments{
\item{x}{a vector containing file paths to the output files generated by \code{\link{meta.retrieval}}.}

\item{gunzip}{a logical value indicating whether or not files should only be renamed (\code{gunzip = FALSE}) or renamed AND unzipped (\code{gunzip}).}
}
\description{
Process the output of \code{\link{meta.retrieval}} by first
un-zipping downloaded files and renaming them for more convenient downstream data analysis.
}
\details{
The output of \code{\link{meta.retrieval}} usually contains compressed sequence files
and a naming convention based on the database the respective file was retrieved from (e.g. \code{Saccharomyces_cerevisiae_cds_from_genomic_refseq.fna.gz}). 
This function helps to format the \code{\link{meta.retrieval}} output files by
\itemize{
\item 1) Automatically uncompress all sequence files in the \code{meta.retrieval} output folder
\item 2) Automatically rename files from e.g. \code{Saccharomyces_cerevisiae_cds_from_genomic_refseq.fna.gz} to \code{Scerevisiae.fa}.
This allows more convenient downstream analyses and visualizations.
}
}
\examples{
\dontrun{
# The easiest way to use 'clean.retrieval()' in combination with
# 'meta.retrieval()' is to use the pipe operator from the 'magrittr' package
library(magrittr)
meta.retrieval(kingdom = "vertebrate_mammalian", 
               db = "refseq", 
               type = "genome") \%>\% clean.retrieval()
}
}
\seealso{
\code{\link{meta.retrieval}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCDSSet.R
\name{getCDSSet}
\alias{getCDSSet}
\title{CDS retrieval of multiple species}
\usage{
getCDSSet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = TRUE,
  gunzip = TRUE,
  update = FALSE,
  path = "set_CDS"
)
}
\arguments{
\item{db}{a character string specifying the database from which the CDS
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the names of the organisms than shall be retrieved.
There are three available options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a CDS shall be downloaded if it isn't marked
in the database as either a reference CDS or a representative CDS.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{clean_retrieval}{logical value indicating whether or not downloaded files shall be renamed for more convenient downstream data analysis.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{update}{a logical value indicating whether or not files that were already downloaded and are still present in the 
output folder shall be updated and re-loaded (\code{update = TRUE} or whether the existing file shall be retained \code{update = FALSE} (Default)).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding CDSs shall be stored. Default is
\code{path} = \code{"set_CDS"}.}
}
\value{
File path to downloaded CDSs.
}
\description{
Main CDS retrieval function for a set of organism of interest.
By specifying the scientific names of the organisms of interest the corresponding fasta-files storing the CDS of the organisms of interest
will be downloaded and stored locally. CDS files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory 'set_CDSs' to store
the CDSs of interest as fasta files for future processing.
In case the corresponding fasta file already exists within the
'set_CDSs' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
getCDSSet("refseq", organisms = c("Arabidopsis thaliana", 
                                  "Arabidopsis lyrata", 
                                  "Capsella rubella"))
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getProteomeSet}},
\code{\link{getRNASet}}, \code{\link{getGFFSet}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_cds}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_genome.R
\name{summary_genome}
\alias{summary_genome}
\title{Retrieve summary statistics for a genome assembly file}
\usage{
summary_genome(file, organism)
}
\arguments{
\item{file}{file path to a genome assembly file in \code{fasta} format.}

\item{organism}{character string specifying the organism at hand.}
}
\description{
A summary statistics of specific genome features is generated.
These statistics are useful to assess the genome quality of retrieved genome assemblies
when performing comparative genomics tasks. This way, users can assess whether or not
patterns found based on genome comparisons aren't just a technical artifact of
differences in genome assembly quality.
}
\details{
The summary statistics include: 
\itemize{
\item \code{genome_size_mbp}: Genome size in mega base pairs
\item \code{n50_mbp}: The N50 contig size of the genome assembly in mega base pairs
\item \code{n_seqs}: The number of chromosomes/scaffolds/contigs of the genome assembly file
\item \code{n_nnn}: The absolute number of NNNs (over all chromosomes or scaffolds or contigs) in the genome assembly file
\item \code{rel_nnn}: The percentage (relative frequency) of NNNs (over all chromosomes or scaffolds or contigs) compared to the total number of 
nucleotides in the genome assembly file
\item \code{genome_entropy}: The \code{Shannon Entropy} of the genome assembly file (median entropy over all individual chromosome entropies)
\item \code{n_gc}: The total number of GCs 
(over all chromosomes or scaffolds or contigs) in the genome assembly file
\item \code{rel_gc}: The (relative frequency) of GCs 
(over all chromosomes or scaffolds or contigs) compared to the total number of 
nucleotides in the genome assembly file
}
}
\examples{
\dontrun{
# retrieve genome from NCBI RefSeq
Sc <- biomartr::getGenome(db = "refseq", organism = "Saccharomyces cerevisiae")
# compute genome assembly summary statistics
Sc_genome_summary <- summary_genome(file = Sc, organism = "Saccharomyces cerevisiae")
# look at results
Sc_genome_summary
}
}
\seealso{
\code{\link{summary_cds}}, \code{\link{getCollection}}, \code{\link{getGenome}}, \code{\link{read_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.retrieval.R
\name{meta.retrieval}
\alias{meta.retrieval}
\title{Perform Meta-Genome Retrieval}
\usage{
meta.retrieval(
  db = "refseq",
  kingdom,
  group = NULL,
  type = "genome",
  restart_at_last = TRUE,
  reference = FALSE,
  combine = FALSE,
  path = NULL
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome 
shall be retrieved:

\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"} 
\item \code{db = "emsembl"}
}}

\item{kingdom}{a character string specifying the kingdom of the organisms 
of interest, e.g. 

\itemize{
\item For \code{NCBI RefSeq}:
\itemize{
\item \code{kingdom = "archaea"}
\item \code{kingdom = "bacteria"}
\item \code{kingdom = "fungi"}
\item \code{kingdom = "invertebrate"}
\item \code{kingdom = "plant"}
\item \code{kingdom = "protozoa"}
\item \code{kingdom = "viral"}
\item \code{kingdom = "vertebrate_mammalian"}
\item \code{kingdom = "vertebrate_other"}
}
\item For \code{NCBI Genbank}:
\itemize{
\item \code{kingdom = "archaea"}
\item \code{kingdom = "bacteria"}
\item \code{kingdom = "fungi"}
\item \code{kingdom = "invertebrate"}
\item \code{kingdom = "plant"}
\item \code{kingdom = "protozoa"}
\item \code{kingdom = "vertebrate_mammalian"}
\item \code{kingdom = "vertebrate_other"}
}
\item For \code{ENSEMBL}:
\itemize{
\item \code{kingdom = "Ensembl"}
}
}

Available kingdoms can be retrieved with \code{\link{getKingdoms}}.}

\item{group}{only species belonging to this subgroup will be downloaded. 
Groups can be retrieved with \code{\link{getGroups}}.}

\item{type}{type of sequences that shall be retrieved. Options are:

\itemize{
 \item \code{type = "genome"} :
 (for genome assembly retrieval; see also \code{\link{getGenome}}), 
 \item \code{type = "proteome"} :
 (for proteome retrieval; see also \code{\link{getProteome}}),
 \item \code{type = "cds"} :
 (for coding sequence retrieval; see also \code{\link{getCDS}}),
 \item \code{type = "gff"} :
(for annotation file retrieval in gff format; see also \code{\link{getGFF}}),
\item \code{type = "gtf"} :
(for annotation file retrieval in gtf format (only for ensembl and
 ensemblgenomes); see also \code{\link{getGTF}})
 \item \code{type = "rna"} :
 (for RNA file retrieval in fasta format; see also \code{\link{getRNA}}),
 \item \code{type = "rm"} :
 (for Repeat Masker output file retrieval; see also 
 \code{\link{getRepeatMasker}}),
 \item \code{type = "assemblystats"} :
 (for genome assembly quality stats file retrieval; 
 see also \code{\link{getAssemblyStats}}).
 }}

\item{restart_at_last}{a logical value indicating whether or not \code{meta.retrieval} should pick up at the last species when re-running the function.
\itemize{
\item If \code{restart_at_last = TRUE} (Default) then \code{meta.retrieval} will skip all organisms that are already present in the folder 
and will start downloading all remaining species. However, this way \code{meta.wretrieval} will not be able to check whether
already downloaded organism files are corrupted or not by checking the md5 checksum.
\item If \code{restart_at_last = FALSE} then \code{meta.retrieval} will start from the beginning and crawl through already downloaded
organism files and check whether already downloaded organism files are corrupted or not by checking the md5 checksum.
After checking existing files the function will start downloading all remaining organisms.
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked in the database 
as either a reference genome or a representative genome. Options are:
\itemize{
\item \code{reference = FALSE} (Default): all organisms (reference, representative, and non-representative genomes) are downloaded.
\item \code{reference = TRUE}: organisms that are downloaded must be either a reference or representative genome. Thus, most genomes which are usually non-reference genomes
will not be downloaded.
}}

\item{combine}{just in case \code{type = "assemblystats"} is specified, shall
assemby stats of individual species be imported and combined to a 
\code{\link{data.frame}}?}

\item{path}{path to the folder in which downloaded genomes shall be stored. 
By default the kingdom name is used to name the output folder.}
}
\value{
a character vector storing the file paths of the retrieved files.
}
\description{
Download genomes, proteomes, cds, gff, rna, or assembly stats 
files of all species within a kingdom of life. After downloading users
can unzip all files using \code{\link{clean.retrieval}}.
}
\details{
This function aims to perform bulk retrieval of the genomes, 
proteomes, cds, etc. of species that belong to the same kingdom of life or 
to the same subgroup.
}
\examples{
\dontrun{
# get all available kingdoms for refseq
getKingdoms(db = "refseq")
# download all vertebrate genomes from refseq
meta.retrieval(kingdom = "vertebrate_mammalian", 
               db = "refseq", 
               type = "genome")

# get all available kingdoms for genbank
getKingdoms(db = "genbank")
# download all vertebrate genomes from genbank
meta.retrieval(kingdom = "vertebrate_mammalian", 
               db = "genbank", 
               type = "genome")


# In case users do not wish to retrieve genomes from an entire kingdom, 
# but rather from a subgoup (e.g. from species belonging to the 
# Gammaproteobacteria class, a subgroup of the bacteria kingdom), 
# they can use the following workflow"
# First, users can again consult the getKingdoms() function to retrieve 
# kingdom information.
getKingdoms(db = "refseq")

# In this example, we will choose the bacteria kingdom. 
# Now, the getGroups() function allows users to obtain available 
# subgroups of the bacteria kingdom.
getGroups(db = "refseq", kingdom = "bacteria")

# Now we choose the group Gammaproteobacteria and specify 
# the group argument in the meta.retrieval() function
meta.retrieval(kingdom = "bacteria", 
   roup = "Gammaproteobacteria", 
   db = "refseq", 
   type = "genome")
}
}
\seealso{
\code{\link{meta.retrieval.all}}, \code{\link{getCollection}}, \code{\link{clean.retrieval}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAttributes.R
\name{getAttributes}
\alias{getAttributes}
\title{Retrieve All Available Attributes for a Specific Dataset}
\usage{
getAttributes(mart, dataset)
}
\arguments{
\item{mart}{a character string specifying the database (mart) 
for which datasets shall be listed.}

\item{dataset}{a character string specifying the dataset for which 
attributes shall be listed.}
}
\description{
This function queries the BioMart Interface and returns a table
storing all available attributes for a specific dataset.
}
\examples{
\dontrun{
# search for available datasets
getMarts()

# choose database (mart): ENSEMBL_MART_ENSEMBL
# and get a table of all available datasets from this BioMart database
head(getDatasets(mart = "ENSEMBL_MART_ENSEMBL"), 10)

# choose dataset: "hsapiens_gene_ensembl"
head(getAttributes(mart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl") , 5)
}
}
\seealso{
\code{\link{getMarts}}, \code{\link{getDatasets}},
\code{\link{getFilters}}, \code{\link{organismBM}}, 
\code{\link{organismFilters}}, \code{\link{organismAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGTF.R
\name{getGTF}
\alias{getGTF}
\title{Genome Annotation Retrieval (GTF)}
\usage{
getGTF(
  db = "ensembl",
  organism,
  remove_annotation_outliers = FALSE,
  path = file.path("ensembl", "annotation"),
  assembly_type = "toplevel",
  release = NULL
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome
shall be retrieved:
\itemize{
\item \code{db = "ensembl"}
}}

\item{organism}{a character string specifying the scientific name of the
organism of interest, e.g. \code{organism = "Homo sapiens"}.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input \code{annotation_file}?
If yes, then the initial \code{annotation_file} will be overwritten and the removed outlier lines will be stored at \code{\link{tempdir}}
for further exploration.}

\item{path}{a character string specifying the location (a folder) in which
the corresponding annotation file shall be stored. Default is
\code{path = file.path("ensembl","annotation")}.}

\item{assembly_type}{a character string specifying from which assembly type the genome
shall be retrieved from (ensembl only, else this argument is ignored):
Default is
\code{assembly_type = "toplevel")}.
This will give you all multi-chromosomes (copies of the same chromosome with small variations).
As an example the toplevel fasta genome in human is over 70 GB uncompressed.
To get primary assembly with 1 chromosome variant per chromosome:
\code{assembly_type = "primary_assembly")}.
As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed:}

\item{release}{a numeric, the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used. \code{release = 75} would for human would give the stable
GRCh37 release in ensembl. Value must be > 46, since ensembl did not structure their data
if the standard format before that.}
}
\value{
File path to downloaded annotation file.
}
\description{
Main retrieval function for GTF files of an
organism of interest. By specifying the scientific name of an organism of
interest the corresponding GTF file storing the annotation  for the organism
of interest can be downloaded and stored locally. GTF files can be retrieved
from several databases.
}
\details{
Internally this function loads the the overview.txt file from ENSEMBL:
and creates a directory 'ensembl/annotation' to store
the genome of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'ensembl/annotation' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
# download the annotation of Homo sapiens from ensembl
# and store the corresponding genome file in 'ensembl/annotation'
getGTF(db            = "ensembl",
       organism      = "Homo sapiens",
       path          = file.path("ensembl","annotation"))

getGTF(db            = "ensembl",
       organism      = "Homo sapiens",
       path          = file.path("ensembl","annotation"),
       assembly_type = "primary_assembly")

}

}
\seealso{
\code{\link{getProteome}}, \code{\link{getCDS}},
\code{\link{getGenome}}, \code{\link{getRNA}}, \code{\link{getRepeatMasker}},
\code{\link{getAssemblyStats}}, \code{\link{meta.retrieval}},
\code{\link{getGFF}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_gff.R
\name{read_gff}
\alias{read_gff}
\title{Import GFF File}
\usage{
read_gff(file)
}
\arguments{
\item{file}{a character string specifying the path to the file 
storing the CDS.}
}
\value{
Either a \code{Biostrings} or \code{data.table} object.
}
\description{
This function reads an organism specific CDS stored 
in a defined file format.
}
\details{
This function takes a string specifying the path to the GFF file
of interest (e.g. the path returned by \code{\link{getGFF}}).
}
\seealso{
\code{\link{getGenome}}, \code{\link{read_genome}}, 
\code{\link{read_proteome}}, \code{\link{read_cds}}, \code{\link{read_rna}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSummaryFile.R
\name{getSummaryFile}
\alias{getSummaryFile}
\title{Helper function to retrieve the assembly_summary.txt file from NCBI}
\usage{
getSummaryFile(db, kingdom)
}
\arguments{
\item{db}{database name. E.g. \code{refseq} or \code{genbank}.}

\item{kingdom}{kingdom for which assembly_summary.txt file shall be 
retrieved. See also \code{\link{getKingdoms}}.}
}
\description{
Retrieval function of the assembly_summary.txt file from NCBI.
}
\examples{
\dontrun{ 
test <- getSummaryFile("refseq","plant")
test
}
}
\seealso{
\code{\link{getKingdomAssemblySummary}}, 
\code{\link{getMetaGenomeSummary}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMetaGenomeSummary.R
\name{getMetaGenomeSummary}
\alias{getMetaGenomeSummary}
\title{Retrieve the assembly_summary.txt file from NCBI genbank metagenomes}
\usage{
getMetaGenomeSummary()
}
\description{
Retrieval function of the assembly_summary.txt file 
from NCBI genbank metagenomes.
This files stores all available metagenome projects on NCBI Genbank.
}
\examples{
\dontrun{
meta.summary <- getMetaGenomeSummary()
meta.summary
}
}
\seealso{
\code{\link{getKingdomAssemblySummary}}, 
\code{\link{getSummaryFile}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listKingdoms.R
\name{listKingdoms}
\alias{listKingdoms}
\title{List number of available genomes in each kingdom of life}
\usage{
listKingdoms(db = "refseq")
}
\arguments{
\item{db}{a character string specifying the database for which genome 
availability shall be checked, 
e.g. \code{db = "refseq"}, \code{db = "genbank"}, \code{db = "ensembl"}, 
\code{db = "ensemblgenomes"}.}
}
\description{
Users can retrieve the available number of sequenced 
genomes per kingdom.
}
\examples{
\dontrun{
# list number of available genomes in refseq for each kingdom of life
listKingdoms(db = "refseq")
# example for genbank
listKingdoms(db = "genbank")
# example for ensembl
listKingdoms(db = "ensembl")
# example for ensemblgenomes
listKingdoms(db = "ensemblgenomes")
}
}
\seealso{
\code{\link{listGenomes}}, \code{\link{is.genome.available}}, 
\code{\link{listGroups}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGENOMEREPORT.R
\name{getGENOMEREPORT}
\alias{getGENOMEREPORT}
\title{Retrieve NCBI GENOME_REPORTS file}
\usage{
getGENOMEREPORT()
}
\description{
Retrieves NCBI GENOME_REPORTS file from 
ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt.
}
\examples{
\dontrun{ 
report <- getGENOMEREPORT()
report
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getENSEMBLGENOMESInfo.R
\name{getENSEMBLGENOMESInfo}
\alias{getENSEMBLGENOMESInfo}
\title{Retrieve ENSEMBLGENOMES info file}
\usage{
getENSEMBLGENOMESInfo()
}
\description{
Retrieve species and genome information from 
http://rest.ensemblgenomes.org/info/species?content-type=application/json/.
}
\examples{
\dontrun{ 
info.file <- getENSEMBLGENOMESInfo()
info.file
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMarts.R
\name{getMarts}
\alias{getMarts}
\title{Retrieve information about available Ensembl Biomart databases}
\usage{
getMarts()
}
\description{
This funcion queries the Ensembl Biomart API and returns a table
storing information about all available Ensembl Biomart databases.
}
\examples{
\dontrun{
# get a table of all available databases from Ensembl Biomart
getMarts()
 }
}
\seealso{
\code{\link{getDatasets}}, \code{\link{getAttributes}}, 
\code{\link{getFilters}}, \code{\link{organismBM}}, 
\code{\link{organismFilters}}, \code{\link{organismAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_database.R
\name{download.database}
\alias{download.database}
\title{Download a NCBI Database to Your Local Hard Drive}
\usage{
download.database(db, path = "database")
}
\arguments{
\item{db}{a character string specifying the database that shall be downloaded
(selected from \code{\link{listDatabases}}).}

\item{path}{a character string specifying the location (a folder) in
which the corresponding database shall be stored.
Default is \code{path = "database"}.
In case this folder does not exist yet, it will be created.}
}
\value{
File path to the downloaded database file.
}
\description{
This function allows users to download a database selected by
\code{\link{listDatabases}} to their local hard drive.
}
\details{
This function downloads large databases to your hard drive.
For this purpose a folder
named \code{database} (default) is created and the correspondning
database then stored in this folder.
}
\examples{
\dontrun{
  # search for available NCBI nr databases
  listNCBIDatabases(db = "nr")
  # select NCBI nr version 27 =  "nr.27.tar.gz"
  # and download it to your hard drive
  # -> please note that large databases take some time for download!
  download.database(db = "nr.27.tar.gz")
}
}
\seealso{
\code{\link{download.database.all}}, \code{\link{listDatabases}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCollection.R
\name{getCollection}
\alias{getCollection}
\title{Retrieve a Collection: Genome, Proteome, CDS, RNA, GFF, Repeat Masker, AssemblyStats}
\usage{
getCollection(
  db = "refseq",
  organism,
  reference = TRUE,
  release = NULL,
  gunzip = FALSE,
  remove_annotation_outliers = FALSE,
  path = file.path("_db_downloads", "collections")
)
}
\arguments{
\item{db}{a character string specifying the database from which the collection 
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organism}{there are three options to characterize an organism: 
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a collection shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning that the most recent database version is used.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{remove_annotation_outliers}{shall outlier lines be removed from the input annotation_file?
If yes, then the initial annotation_file will be overwritten and the removed outlier lines 
will be stored at \code{\link{tempdir}} for further exploration.}

\item{path}{a character string specifying the location (a folder) in which 
the corresponding collection shall be stored. Default is 
\code{path} = \code{file.path("_db_downloads","collections")}.}
}
\value{
File path to downloaded genome.
}
\description{
Main collection retrieval function for an organism of interest.
By specifying the scientific name of an organism of interest a collection consisting of
the genome file, proteome file, CDS file, RNA file, GFF file, Repeat Masker file, AssemblyStats 
file of the organism of interest
can be downloaded and stored locally. Collections can be retrieved from 
several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
 
and creates a directory '_ncbi_downloads/collection' to store
the genome of interest as fasta file for future processing.
In case the corresponding fasta file already exists within the
'_ncbi_downloads/collection' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{

# download the collection of Arabidopsis thaliana from refseq
# and store the corresponding genome file in '_ncbi_downloads/collection'
 getCollection( db       = "refseq", 
             organism = "Arabidopsis thaliana", 
             path = file.path("_db_downloads","collections"))
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getProteomeSet}}, \code{\link{getCDSSet}}, 
\code{\link{getGenome}}, \code{\link{getProteome}}, \code{\link{getCDS}}, 
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}}, 
\code{\link{read_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listGenomes.R
\name{listGenomes}
\alias{listGenomes}
\title{List All Available Genomes either by kingdom, group, or subgroup}
\usage{
listGenomes(db = "refseq", type = "all", subset = NULL, details = FALSE)
}
\arguments{
\item{db}{a character string specifying the database for which genome 
availability shall be checked. Available options are:
\itemize{
\item \code{db = "refseq"} 
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{type}{a character string specifying a potential filter of available 
genomes. Available options are:
\itemize{
\item \code{type = "all"}
\item \code{type = "kingdom"} 
\item \code{type = "group"}
\item \code{type = "subgroup"}
}}

\item{subset}{a character string or character vector specifying a subset of 
\code{type}. E.g. if users are interested in retrieving all
\code{Eukaryota} species, they can specify: \code{type = "kingdom"} and 
\code{subset = "Eukaryota"}.}

\item{details}{a boolean value specifying whether only the scientific names 
of stored genomes shall be returned (details = FALSE) or all information such
as 
\itemize{
\item \code{organism_name}
\item \code{kingdoms}
\item \code{group}
\item \code{subgroup} 
\item \code{file_size_MB}, etc.
}}
}
\description{
This function retrieves the names of all genomes available on 
the NCBI ftp:// server and stores the results in a file named 'overview.txt' 
inside the directory _ncbi_downloads' that is built inside the workspace.
}
\details{
Internally this function loads the the overview.txt file from NCBI 
and creates a directory '_ncbi_downloads' in the \code{temdir()}
folder to store the overview.txt file for future processing. In case the 
overview.txt file already exists within the '_ncbi_downloads' folder and is 
accessible within the workspace, no download process will be performed again.
}
\note{
Please note that the ftp:// connection relies on the NCBI or ENSEMBL server 
and cannot be accurately accessed via a proxy.
}
\examples{
\dontrun{
# print details for refseq
listGenomes(db = "refseq") 
# print details for all plants in refseq
listGenomes(db = "refseq", type = "kingdom")
# print details for all plant groups in refseq
listGenomes(db = "refseq", type = "group")
# print details for all plant subgroups in refseq
listGenomes(db = "refseq", type = "subgroup")
}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listGroups.R
\name{listGroups}
\alias{listGroups}
\title{List number of available genomes in each taxonomic group}
\usage{
listGroups(db = "refseq", kingdom = "all", details = FALSE)
}
\arguments{
\item{db}{a character string specifying the database for which genome 
availability shall be checked. Available options are:
\itemize{
\item \code{db = "refseq"} 
\item \code{db = "genbank"}
}}

\item{kingdom}{a kingdom specification retrieved by 
\code{\link{getKingdoms}}.}

\item{details}{shall all species corresponding to the specified 
\code{kingdom} be returned? Default is \code{details = FALSE}.}
}
\description{
Users can retrieve the available number of sequenced 
genomes per group. Only available for \code{db = "refseq"} and 
\code{db = "genbank"}.
}
\examples{
\dontrun{
# example for refseq
listGroups(db = "refseq")
# example for genbank
listGroups(db = "genbank")
### in case groups should be specified by kingdom
# first, retrieve available kingdom names
listKingdoms()
# now we choose kingdom "bacteria"
listGroups(db = "refseq", kingdom = "bacteria")
# or
listGroups(db = "genbank", kingdom = "bacteria")
}
}
\seealso{
\code{\link{listGenomes}}, \code{\link{is.genome.available}}, 
\code{\link{listKingdoms}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGO.R
\name{getGO}
\alias{getGO}
\title{Gene Ontology Query}
\usage{
getGO(organism, genes, filters, ...)
}
\arguments{
\item{organism}{a character string specifying the scientific name 
of a query organism.}

\item{genes}{a character vector storing the gene ids of a organisms 
of interest to be queried against Ensembl Biomart.}

\item{filters}{a character vector specifying the filter (query key) for 
the Ensembl Biomart query, e.g. \code{filter} = \code{"ensembl_gene_id"}.}

\item{...}{additional parameters that can be passed to the 
\code{\link{biomart}} function.}
}
\description{
This function takes a gene id as character vector from a given 
query organism and returns the corresponding GO terms and additional GO 
information.
}
\details{
This function takes the scientific name of a query organism, a set of genes 
for which GO terms and additional information shall be retrieved, and a 
filter argument that specifies the attribute for the query genes.
}
\examples{
\dontrun{ 
GO_tbl <- getGO(organism = "Arabidopsis thaliana", 
                genes    = c("AT1G06090", "AT1G06100"),
                filters  = "ensembl_gene_id")

# look at the result
head(GO_tbl)
}
}
\seealso{
\code{\link{biomart}}, \code{\link{organismFilters}}, 
\code{\link{organismBM}}, \code{\link[biomaRt]{getBM}}, 
\code{\link{getMarts}}, \code{\link{getDatasets}}, 
\code{\link{getFilters}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organismBM.R
\name{organismBM}
\alias{organismBM}
\title{Retrieve Ensembl Biomart marts and datasets for a query organism}
\usage{
organismBM(organism = NULL, update = FALSE)
}
\arguments{
\item{organism}{a character string specifying the scientific name of a 
query organism. Default is \code{organism} = \code{NULL}. In this case all 
available biomart connections are returned.}

\item{update}{a logical value specifying whether or not the local 
listMart.txt and listDatasets.txt files shall be updated by remote access
 to BioMart.}
}
\description{
This function returns either all available biomart connections 
for all available organisms for which biomart access is possible, or 
(when specified) returns all organism specific biomart connections.
}
\details{
This function collects all available biomart connections and returns a table 
storing the organism for which biomart connections are available as well as 
the corresponding mart and database.
}
\note{
When you run this function for the first time, the data retrieval procedure 
will take some time, due to the remote access to BioMart. The corresponding 
result is then saved in a *.txt file named "_biomart/listDatasets.txt" in the
\code{\link{tempdir}} directory, allowing subsequent queries to perform
 much faster.
}
\examples{
\dontrun{ 
# returning all available biomart connections
head(organismBM(), 20)
# retrieving all available datasets and biomart connections for
# a specific query organism (scientific name)
organismBM(organism = "Homo sapiens")
# you can also update the downloaded version using 
# the "update = TRUE" argument
head(organismBM(update = TRUE), 20)
}
}
\references{
\url{http://biomart.org/}

Mapping identifiers for the integration of genomic datasets with the
R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan
Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

BioMart and Bioconductor: a powerful link between biological databases and
microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean
Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21,
3439-3440 (2005).
}
\seealso{
\code{\link{getMarts}}, \code{\link{getDatasets}}, 
\code{\link{biomart}}, \code{\link{organismFilters}}, 
\code{\link{organismAttributes}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProteomeSet.R
\name{getProteomeSet}
\alias{getProteomeSet}
\title{Proteome retrieval of multiple species}
\usage{
getProteomeSet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = TRUE,
  gunzip = TRUE,
  update = FALSE,
  path = "set_proteomes"
)
}
\arguments{
\item{db}{a character string specifying the database from which the proteome
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the names of the organisms than shall be retrieved.
There are three available options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a proteome shall be downloaded if it isn't marked
in the database as either a reference proteome or a representative proteome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{clean_retrieval}{logical value indicating whether or not downloaded files shall be renamed for more convenient downstream data analysis.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{update}{a logical value indicating whether or not files that were already downloaded and are still present in the 
output folder shall be updated and re-loaded (\code{update = TRUE} or whether the existing file shall be retained \code{update = FALSE} (Default)).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding proteomes shall be stored. Default is
\code{path} = \code{"set_proteomes"}.}
}
\value{
File path to downloaded proteomes.
}
\description{
Main proteome retrieval function for a set of organism of interest.
By specifying the scientific names of the organisms of interest the corresponding fasta-files storing the proteome of the organisms of interest
will be downloaded and stored locally. proteome files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory 'set_proteomes' to store
the proteomes of interest as fasta files for future processing.
In case the corresponding fasta file already exists within the
'set_proteomes' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
getProteomeSet("refseq", organisms = c("Arabidopsis thaliana", 
                                      "Arabidopsis lyrata", 
                                       "Capsella rubella"))
}

}
\seealso{
\code{\link{getGenomeSet}}, \code{\link{getCDSSet}},
\code{\link{getRNASet}}, \code{\link{getGFFSet}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_proteome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_database_all.R
\name{download.database.all}
\alias{download.database.all}
\title{Download all elements of an NCBI databse}
\usage{
download.database.all(db, path = NULL)
}
\arguments{
\item{db}{a character string specifying the database that shall be downloaded
(selected from \code{\link{listDatabases}}).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding
database shall be stored. In case this folder does not exist yet,
it will be created.}
}
\value{
A character vector storing the file paths of the downloaded databases.
}
\description{
The \code{\link{download.database}} functions allows users to
retrieve individual packages of a NCBI database. This function is designed to
retrieve the entire database selected by the users (hence all packages 
corresponding to this database).
}
\examples{
\dontrun{
# search for available NCBI databases
  listNCBIDatabases(db = "all")
# choose database NCBI nr and download compelete database
  download.database.all(db = "nr", path = "nr")
}
}
\seealso{
\code{\link{download.database}}, \code{\link{listNCBIDatabases}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_cds.R
\name{summary_cds}
\alias{summary_cds}
\title{Retrieve summary statistics for a coding sequence (CDS) file}
\usage{
summary_cds(file, organism)
}
\arguments{
\item{file}{file path to a CDS file in \code{fasta} format.}

\item{organism}{character string specifying the organism at hand.}
}
\description{
A summary statistics of specific CDS features is returned.
}
\details{
The summary statistics include:
\itemize{
\item \code{total_seqs}: 
\item \code{nnn_abs}: The total number of NNN's 
(over all chromosomes/scaffolds/contigs) in all coding sequences combined
\item \code{nnn_perc}: The percentage (relative frequency) of NNN's 
(over all chromosomes/scaffolds/contigs) compared to the total number of 
nucleotides of all coding sequences
}
}
\seealso{
\code{\link{getCollection}}, \code{\link{getCDS}}, \code{\link{read_cds}}, \code{\link{summary_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenomeSet.R
\name{getGenomeSet}
\alias{getGenomeSet}
\title{Genome Retrieval of multiple species}
\usage{
getGenomeSet(
  db = "refseq",
  organisms,
  reference = FALSE,
  release = NULL,
  clean_retrieval = TRUE,
  gunzip = TRUE,
  update = FALSE,
  path = "set_genomes",
  assembly_type = "toplevel"
)
}
\arguments{
\item{db}{a character string specifying the database from which the genome
shall be retrieved:
\itemize{
\item \code{db = "refseq"}
\item \code{db = "genbank"}
\item \code{db = "ensembl"}
}}

\item{organisms}{a character vector storing the names of the organisms than shall be retrieved.
There are three available options to characterize an organism:
\itemize{
\item by \code{scientific name}: e.g. \code{organism = "Homo sapiens"}
\item by \code{database specific accession identifier}: e.g. \code{organism = "GCF_000001405.37"} (= NCBI RefSeq identifier for \code{Homo sapiens})
\item by \code{taxonomic identifier from NCBI Taxonomy}: e.g. \code{organism = "9606"} (= taxid of \code{Homo sapiens})
}}

\item{reference}{a logical value indicating whether or not a genome shall be downloaded if it isn't marked
in the database as either a reference genome or a representative genome.}

\item{release}{the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used.}

\item{clean_retrieval}{logical value indicating whether or not downloaded files shall be renamed for more convenient downstream data analysis.}

\item{gunzip}{a logical value indicating whether or not files should be unzipped.}

\item{update}{a logical value indicating whether or not files that were already downloaded and are still present in the 
output folder shall be updated and re-loaded (\code{update = TRUE} or whether the existing file shall be retained \code{update = FALSE} (Default)).}

\item{path}{a character string specifying the location (a folder) in which
the corresponding genomes shall be stored. Default is
\code{path} = \code{"set_genomes"}.}

\item{assembly_type}{a character string specifying from which assembly type the genome
shall be retrieved from (ensembl only, else this argument is ignored):
Default is
\code{assembly_type = "toplevel")}.
This will give you all multi-chromosomes (copies of the same chromosome with small variations).
As an example the toplevel fasta genome in human is over 70 GB uncompressed.
To get primary assembly with 1 chromosome variant per chromosome:
\code{assembly_type = "primary_assembly")}.
As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed:}
}
\value{
File path to downloaded genomes.
}
\description{
Main genome retrieval function for a set of organism of interest.
By specifying the scientific names of the organisms of interest the corresponding fasta-files storing the genome of the organisms of interest
will be downloaded and stored locally. Genome files can be retrieved from several databases.
}
\details{
Internally this function loads the the overview.txt file from NCBI:

 refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

 genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/

and creates a directory 'set_genomes' to store
the genomes of interest as fasta files for future processing.
In case the corresponding fasta file already exists within the
'set_genomes' folder and is accessible within the workspace,
no download process will be performed.
}
\examples{
\dontrun{
getGenomeSet("refseq", organisms = c("Arabidopsis thaliana", 
                                     "Arabidopsis lyrata", 
                                     "Capsella rubella"))
}

}
\seealso{
\code{\link{getProteomeSet}}, \code{\link{getCDSSet}},
\code{\link{getRNASet}}, \code{\link{getGFFSet}}, \code{\link{getCDS}},
\code{\link{getGFF}}, \code{\link{getGTF}}, \code{\link{getRNA}}, \code{\link{meta.retrieval}},
\code{\link{read_genome}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_rna.R
\name{read_rna}
\alias{read_rna}
\title{Import RNA as Biostrings or data.table object}
\usage{
read_rna(file, format = "fasta", obj.type = "Biostrings", ...)
}
\arguments{
\item{file}{a character string specifying the path to the file 
storing the RNA.}

\item{format}{a character string specifying the file format used to store the
genome, e.g. \code{format = "fasta"} (default) or \code{format = "gbk"}.}

\item{obj.type}{a character string specifying the object stype in which the 
genomic sequence shall be represented. 
Either as \code{obj.type = "Biostrings"} (default) or 
as \code{obj.type = "data.table"}.}

\item{...}{additional arguments that are used by 
\code{\link[seqinr]{read.fasta}}.}
}
\value{
A data.table storing the gene id in the first column and the 
corresponding sequence as string in the second column.
}
\description{
This function reads an organism specific RNA stored in a 
defined file format.
}
\details{
This function takes a string specifying the path to the RNA file
of interest as first argument. It is possible to read in different proteome 
file standards such as \emph{fasta} or \emph{genebank}.
}
\seealso{
\code{\link{getRNA}}, \code{\link{read_genome}}, 
\code{\link{read_proteome}}, \code{\link{read_gff}}
}
\author{
Hajk-Georg Drost
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getENSEMBL.gtf.R
\name{getENSEMBL.gtf}
\alias{getENSEMBL.gtf}
\title{Helper function for retrieving gff files from ENSEMBL}
\usage{
getENSEMBL.gtf(
  organism,
  type = "dna",
  id.type = "toplevel",
  path,
  release = NULL
)
}
\arguments{
\item{organism}{scientific name of the organism of interest.}

\item{type}{biological sequence type.}

\item{id.type}{a character, default "toplevel". id type of assembly, either toplevel or primary_assembly usually.}

\item{path}{location where file shall be stored.}

\item{release}{a numeric, the database release version of ENSEMBL (\code{db = "ensembl"}). Default is \code{release = NULL} meaning
that the most recent database version is used. \code{release = 75} would for human would give the stable
GRCh37 release in ensembl. Value must be > 46, since ensembl did not structure their data
if the standard format before that.}
}
\value{
character filepath to download file, returns FALSE if failed.
}
\description{
This function downloads gff
files of query organisms from ENSEMBL.
}
\author{
Hajk-Georg Drost
}
