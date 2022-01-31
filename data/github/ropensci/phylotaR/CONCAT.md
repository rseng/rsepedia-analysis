# Automated Retrieval of Orthologous DNA Sequences from GenBank <img src="https://raw.githubusercontent.com/ropensci/phylotaR/master/logo.png" height="300" align="right"/>
[![Build Status](https://travis-ci.org/ropensci/phylotaR.svg?branch=master)](https://travis-ci.org/ropensci/phylotaR) [![Coverage Status](https://coveralls.io/repos/github/ropensci/phylotaR/badge.svg?branch=master)](https://coveralls.io/github/ropensci/phylotaR?branch=master) [![](https://badges.ropensci.org/187_status.svg)](https://github.com/ropensci/onboarding/issues/187) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/phylotaR)](https://CRAN.R-project.org/package=phylotaR)

R implementation of the PhyLoTa sequence cluster pipeline. For more information see the accompanying website. Tested and demonstrated on Unix and Windows. **Find out more by visiting the [phylotaR website](https://ropensci.github.io/phylotaR/).**

## Install
From CRAN:

```r
install.packages('phylotaR')
```

Or, download the development package from GitHub:

```r
devtools::install_github(repo='ropensci/phylotaR', build_vignettes=TRUE)
```

**Full functionality depends on a local copy of BLAST+ (>= 2.0.0)**. For details on downloading and compiling BLAST+ on your machine please visit the [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

## Pipeline

`phylotaR` runs the PhyLoTa pipeline in four automated stages: identify and retrieve taxonomic information on all descendent nodes of the taxonomic group of interest (`taxise`), download sequence data for every identified node (`download`), identify orthologous clusters using BLAST (`cluster`), and identify sister clusters for sets of clusters identified in the previous stage (`cluster^2`) After these stages are complete, `phylotaR` provides tools for exploring, identifying and exporting suitable clusters for subsequent analysis.

![phylotaR pipeline](https://raw.githubusercontent.com/ropensci/phylotaR/master/other/stages.png)

For more information on the pipeline and how it works see the publication, [phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R](https://doi.org/10.3390/life8020020).

## Running

At a minimum all a user need do is provide the taxonomic ID of their chosen taxonomic group of interest. For example, if you were interested in primates, you can visit the [NCBI taxonomy home page](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/) and search primates to look up their ID. After identifying the ID, the `phylotaR` pipeline can be run with the following script.

```r
library(phylotaR)
wd <- '[FILEPATH TO WORKING DIRECTORY]'
ncbi_dr <- '[FILEPATH TO COMPILED BLAST+ TOOLS]'
txid <- 9443  # primates ID
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
run(wd = wd)
```

The pipeline can be stopped and restarted at any point without loss of data. For more details on this script, how to change parameters, check the log and details of the pipeline, please check out the package vignette.

```r
library(phylotaR)
vignette("phylotaR")
```

## Timings

How long does it take for a phylotaR pipeline to complete? Below is a table listing the runtimes in minutes for different demonstration, taxonomic groups. 

Taxon|N. taxa|N. sequences|N. clusters|Taxise (mins.)|Download (mins.)|Cluster (mins.)|Cluster2 (mins.)|Total (mins.)|
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
Anisoptera|1175|11432|796|1.6|23|48|0.017|72|
Acipenseridae|51|2407|333|0.1|6.9|6.4|0.017|13|
Tinamiformes|25|251|98|0.067|2.4|0.18|0.017|2.7|
Aotus|13|1499|193|0.067|3.2|0.6|0|3.9|
Bromeliaceae|1171|9833|724|1.2|28|37|0.033|66|
Cycadidae|353|8331|540|0.32|19|18|0.033|37|
Eutardigrada|261|960|211|0.3|11|1.8|0.05|14|
Kazachstania|40|623|101|0.1|20|3|0.05|23|
Platyrrhini|212|12731|3112|0.35|51|6.9|1.2|60|

To run these same demonstrations see [demos/demo_run.R](https://github.com/ropensci/phylotaR/blob/master/demos/demo_run.R).

## License

MIT

## Version

Version 1.

## Authors

Dom Bennett (maintainer, R package dev), Hannes Hettling (workhouse code dev), Rutger Vos, Alexander Zizka and Alexandre Antonelli

## Reference

Bennett, D., Hettling, H., Silvestro, D., Zizka, A., Bacon, C., Faurby, S., … Antonelli, A. (2018). phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R. *Life*, **8**(2), 20. [DOI:10.3390/life8020020](https://doi.org/10.3390/life8020020)

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: Processing GenBank for molecular phylogenetics research. *Systematic Biology*, **57**(3), 335–346. [DOI:10.1080/10635150802158688](https://doi.org/10.1080/10635150802158688)

----

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# phylotaR 1.2.0

## `outsider` integration

* Now work with `outsider`, reduing user need to install BLAST+
* New parameter: `db_only`, avoid downloads by fetching exclusively from
`restez` database
* Excluding Transcriptome Shotgun Assembly sequences, making downloading faster

# phylotaR 1.1.1

## Minor fixes

* Validity check for `SeqRec`

# phylotaR 1.1.0

## restez integration

* Now works with [`restez`](https://ropensci.github.io/restez/articles/4_phylotar.html)
* Multiple IDs in `txid` argument
* Fancy new colours
* Easier process killing

# phylotaR 1.0.0
# Contributing

You are very welcome to help out in the development of phylotaR. If you have any ideas for future features than please add them to the [issues page](https://github.com/AntonelliLab/phylotaR/issues). If you have the impetus to add those features yourself, then please fork and send a pull request.

## Areas for possible contribution

### Alternatives to BLAST

Currently, phylotaR only makes use of BLAST. BLAST is good because it can work with any sequence type and is very sensitive. It is, however, slower than alternatives local alignment search tools. A set of generic input, run and output functions for working with any BLAST-alternative with an accompanying vignette would allow a user to any of the available alternatives. 

### BLAST API

A big issue with phylotaR is the need to install and run a local copy of BLAST. One possibility is add API functionality so that a user can send the jobs to the cloud instead of having to install their own version of BLAST.

### Inputting one's own sequences

A user may wish to make use of sequences they have generated themselves in conjunction with those available from GenBank. Currently there is no effort to allow a user to do this. It is a little tricky to do because of phylotaR's reliance on IDs.

### Multiple taxids

Many users wish to run a single phylotaR run for multiple, potentially paraphyletic, taxonomic IDs.

### A user determined taxonomy

Currrently, phylotaR depends on NCBI's taxonomy. In theory, a user could provide their own newick tree representing their preferred taxonomy instead.

### RefSeq

Identify sequences that are orthologous to a specified sequence.

## How to contribute

To contribute you will need a GitHub account and to have basic knowledge of the R language. You can then create a fork of the
repo in your own GitHub account and download the repository to your local machine. `devtools` is recommended.

```r
devtools::install_github('[your account]/phylotaR')
```

All new functions must be tested. For every new file in `R/`, a new test file must be created in `tests/testthat/`. To test the
package and make sure it meets CRAN guidelines use `devtools`. 

```r
devtools::test()
devtools::check_cran()
```

For help, refer to Hadley Wickham's book, [R packages](http://r-pkgs.had.co.nz/).

## Style guide

phylotaR is being developed for submission to ROpenSci. This means the package and its code should meet ROpenSci style and
standards. For example, function names should be all lowercase, separated by underscores and the last word should, ideally, be
a verb.

```
# e.g.
species_ids_retrieve()  # good
sppIDs()                # not great
sp.IDS_part2()          # really bad
sigNXTprt.p()           # awful
```

It is best to make functions small, with specific names. Feel free to break up code into multiple separate files (e.g. tools,
helper functions, stages ...). For more details and better explanations refer to the ROpenSci [guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md).
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->


# CONTRIBUTING #

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the phylotaR project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 

<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>

---
title: "phylotaR Tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{phylotaR Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

The first step to running a phylogenetic analysis is the identification of overlapping sequences. Often orthology is determined by pairing sequences whose gene names match (e.g. COI sequences with COI sequences, rbcl sequences with rbcl sequences). Problems can arise however if gene names differ between authors, if different gene sections are represented or if sequences are mislabelled. These issues can be especially problematic for large-scale analyses where individual errors cannot be detected.

[PhyLoTa](https://www.ncbi.nlm.nih.gov/pubmed/18570030) is a pipeline that uses an alignment search tool to identify orthologous sequences without the need for gene name matching. For a given parental taxonomic group, the pipeline will search through available sequences hosted on GenBank and identify orthologous sequence clusters. A user is then able to survey the identified clusters and select the ones which best suit their phylogenetic analysis needs, e.g. by selecting the clusters that maximise the number of taxonomic groups.

This R pacakge, `phylotaR`, is an R implementation of this pipeline. In this vignette we will demonstrate how to run PhyLoTa using a small taxonomic group. The pipeline is composed of four automated stages (taxise, download, cluster, cluster2) and a final user-performed stage of cluster selection.


# Installing NCBI BLAST+ Tools

The PhyLoTa pipeline uses BLAST to identify orthologous sequence clusters. In order to run phylotaR, a local copy of the BLAST software must be installed on your computer. **Installing the phylotaR package does not install BLAST, it must be installed separately**. To install BLAST+, please see the NCBI website's [installation instructions](https://www.ncbi.nlm.nih.gov/books/NBK279671/).

# Pipeline

## Setup

For demonstration purposes we will run the pipeline on a small taxonomic group. Because they are charismatic and relatively well-studied, we will select the Night Monkey genus, [Aotus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=9504). Now that we have decided on a taxonomic group we need to find out its unique taxonomic ID. This can be looked up by navigating to the [NCBI taxonomy webpage](https://www.ncbi.nlm.nih.gov/taxonomy) and searching 'Aotus'. Doing this, we can see that Aotus ID is **9504**. We will need this number for specifying the parameters in our pipeline. (Notice, that there is also a plant genus called Aotus.)

To begin a run, we will need to create a new folder that will contain all the output files generated by the `phylotaR` pipeline. Since we are running the analysis on the Aotus genus, let's call the folder `aotus/`. Now we have our working directory folder created, we can now open R and run the following code.

```{r setup, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
library(phylotaR)
wd <- '[YOUR PATH TO AOTUS FOLDER]'
ncbi_dr <- '[YOUR PATH TO NCBI BLAST TOOLS]'
txid <- 9504
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr, v = TRUE)
```

The above imports the `phylotaR` package and initiates a cache that will contain the pipeline parameters. For this tutorial we will keep the parameters as their default. See the function `parameters()` for a complete list and description of all the parameters and their default values. For more detailed information on the parameters please see the publication, [phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R](https://doi.org/10.3390/life8020020). `wd` must be a file path to the folder we called `aotus/`. `ncbi_dr` must be a file path to the folder containing all the NCBI BLAST+ tools -- see above 'Installing NCBI BLAST+ Tools'. Depending on your system and how you installed the tools, they may be in your system path in which case you can simply supply '.' to the `ncbi_dr` argument. On my computer I provide the path to the where the `blastn` executable is located, e.g. `/usr/local/ncbi/blast/bin/`. Running `setup()` will verify whether the BLAST tools are installed correctly.

## Running

After `setup()` has been run we can run the pipeline with the following command.

```{r running, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
run(wd = wd)
```

This will run all the automated stages of the pipeline: taxise, download, cluster and cluster2. The first of these stages looks up all the taxonomic information available on the descendants of the parent ID provided, `txid`. The second downloads representative sequences for all identified descendants. No additional arguments are required other than `wd` which specifies the working directory that contains the cache and all parameters as set up by `setup()`. In this folder you will also find a `log.txt` that reports detailed information on the progression of the pipeline as well as all the output files generated by each stage. Additionally, you will see session info and a blast version text files. These files, along with the log, can help debugging if any errors occur. The whole pipeline can complete in around 2 minutes for Aotus using default parameters. Aotus, however, is a genus of only 13 taxa, larger clades will take much longer particularly during the download stage.

## Restarting

The pipeline can be halted and restarted. The cache records all downloaded and generated data by the pipeline. If there is a system crash or the user wishes to halt the program, the pipeline can be restarted from the same point it stopped with the function `restart()`. Additionally, due to the potential random nature of the pipeline, a user may wish to re-run the pipeline from certain stages. This can be achived by first using `reset()` followed by `restart()`. For example, in the code below a completed pipeline is reset to 'cluster' and then restarted. After running these commands, the pipeline will run as if it has only just completed the download stage. Note, all resets and restarts are recorded in the log.

```{r restarting, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
reset(wd = wd, stage = 'cluster')
restart(wd = wd)
```

### Changing parameters

Paramaters can always be set by a user at the initiation of a folder with the `setup()` function. To change the parameter values after a folder has already been set up, a user can use `parameters_reset()`. For example, if the download stage is taking particularly long, the `btchsz` could be increased. This would raise the number of sequences downloaded per request. (Note, too high a `btchsz` may cause your NCBI Entrez access being limited.)

```{r parameters reset, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
# use ctrl+c or Esc to halt R
# increase the btchsz from the default to 300
parameters_reset(wd = wd, parameters = 'btchz', values = 300)
restart(wd = wd)
# ^ restart from whatever point it was halted
```

## Cluster selection

After a pipeline has completed, the identified clusters can be interrogated. We can generate a phylota object using `read_phylota()` but in the code below we will load a pre-existing phylota object from the package data. The phylota object contains cluster, sequence and taxonomic information on all the clusters. It has 6 data slots: cids, sids, txids, txdct, sqs, clstrs, prnt_id and prnt_nm. Each of these slots can be accessed with `@`, see ?\`Phylota-class\` for more information. The `phylotaR` package has a range of functions for probing clusters in a phylota object. For example, if we want to know how many different taxonomic groups are represented by each cluster we can use `get_ntaxa()`.

```{r selection1, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
library(phylotaR)
# pre-load already run aotus from package data
data('aotus')
all_clusters <- aotus
print(all_clusters)
# otherwise, run:
# all_clusters <- read_phylota(wd)
cids <- all_clusters@cids
n_taxa <- get_ntaxa(phylota = all_clusters, cid = cids)
```

We can then drop all the clusters with fewer than 6 taxa and create a new phylota object using the `drop_cls()` function. Let's then take a peak of the now smaller object's clusters using `summary()`.

```{r selection2, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
keep <- cids[n_taxa > 6]
selected <- drop_clstrs(phylota = all_clusters, cid = keep)
smmry <- summary(selected)
print(smmry)
```

This summary provides information on each cluster in the phylota object, such as median sequence length, MAD score (the mean alignment density, values closer to 1 indicate all the sequences are of a similar length), most common words in the sequence descriptions and feature names. Let's select the second ID in this table for further investigation. We can extract its cluster and sequences records in the following way.

```{r selection3, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
cid <- smmry[2, 'ID']
# get the cluster record
cluster_record <- selected@clstrs[[cid]]
# use the seq. IDs to get the sequence records
seq_records <- selected@sqs[cluster_record@sids]
# extract a single record
seq_record <- seq_records[[seq_records@ids[[1]]]]
summary(seq_record)
# get the sequence
seq <- rawToChar(seq_record@sq)
print(substr(x = seq, start = 1, stop = 80))
```

We could extract and write out each of the sequences for a cluster in the above manner. Handily, however, `phylotaR` comes with some tools to make outputting sequences easier. First because there are multiple sequences per taxon, we need to select a single representative sequence. We can do this with the `drop_by_rank()` function. With this function we choose a taxonomic rank at which we would like our sequences to be represented. The function then chooses the 'best' sequence representing each taxon for that rank using a range of criteria. With this new phylota object, we can then extract the scientific names and write out the sequences.

```{r selection4, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
# choose best sequence per species
reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1)
# get txids at the species level for each sequence
txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
# look up name for txids
scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# clean the names
scientific_names <- gsub('\\.', '', scientific_names)
scientific_names <- gsub('\\s+', '_', scientific_names)
print(scientific_names)
# look up sequence IDs for our chosen cluster
sids <- reduced@clstrs[[cid]]@sids
# write out
write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names,
          outfile = file.path(tempdir(), 'cytb.fasta'))
# ^ to avoid clutter, we're writing to a temporary folder
```

## Testing output

We can sanity check our cluster sequences by running a very quick phylogenetic analysis using mafft and raxml. The below code will use the cluster to generate an alignment and a tree through R. In order for the code to run, it requires the installation of mafft and raxml and, additionally, may require tweaking to work on your system.

```{r testing, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE}
library(ape)
system('mafft --auto cytb.fasta > alignment.fasta')
system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n aotus -s alignment.fasta'))
tree <- read.tree(file = 'RAxML_bestTree.aotus')
plot(tree, no.margin = TRUE, type = 'unrooted')
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{cache_setup}
\alias{cache_setup}
\title{Set-up a cache}
\usage{
cache_setup(ps, ovrwrt = FALSE)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}

\item{ovrwrt}{Overwrite existing cache? Default FALSE.}
}
\description{
Creates a cache of parameters in the wd.
}
\details{
Warning: overwriting with this function will delete the
existing cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{clade_select}}, \code{\link{clstr2_calc}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4-tools.R
\name{clstrs_merge}
\alias{clstrs_merge}
\title{Merge joined clusters}
\usage{
clstrs_merge(jnd_clstrs, txdct)
}
\arguments{
\item{jnd_clstrs}{List of joined clusters}

\item{txdct}{Taxonomic dictionary}
}
\value{
list of ClstrRecs
}
\description{
Takes a list of joined clusters and computes each
data slot to create a single merged cluster. txdct is required for
parent look-up.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-api.R
\name{batcher}
\alias{batcher}
\title{Download in batches}
\usage{
batcher(ids, func, ps, lvl = 0)
}
\arguments{
\item{ids}{Vector of record ids}

\item{func}{Downloader function}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
Vector of records

vector of rentrez function results
}
\description{
Run downloader function in batches for sequences or
taxonomic records
}
\seealso{
Other run-private: \code{\link{blast_clstr}},
  \code{\link{blast_filter}}, \code{\link{blast_setup}},
  \code{\link{blast_sqs}}, \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{Phylota-class}
\alias{Phylota-class}
\alias{Phylota-method}
\alias{as.character,Phylota-method}
\alias{show,Phylota-method}
\alias{print,Phylota-method}
\alias{str,Phylota-method}
\alias{summary,Phylota-method}
\alias{[[,Phylota,character-method}
\title{Phylota object}
\usage{
\S4method{as.character}{Phylota}(x)

\S4method{show}{Phylota}(object)

\S4method{print}{Phylota}(x)

\S4method{str}{Phylota}(object, max.level = 2L, ...)

\S4method{summary}{Phylota}(object)

\S4method{[[}{Phylota,character}(x, i)
}
\arguments{
\item{x}{\code{Phylota} object}

\item{object}{\code{Phylota} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}

\item{i}{Either sid or cid}
}
\description{
Phylota table contains all sequence, cluster and taxonomic
information from a phylotaR pipeline run.
}
\section{Slots}{

\describe{
\item{\code{cids}}{IDs of all clusters}

\item{\code{sids}}{IDs of all sequences}

\item{\code{txids}}{IDs of all taxa}

\item{\code{sqs}}{All sequence records as SeqArc}

\item{\code{clstrs}}{All cluster records as ClstrArc}

\item{\code{txdct}}{Taxonomic dictionary as TaxDict}

\item{\code{prnt_id}}{Parent taxonomic ID}

\item{\code{prnt_nm}}{Parent taxonomic name}
}}

\examples{
data('aotus')
# this is a Phylota object
# it contains cluster, sequence and taxonomic information from a phylotaR run
show(aotus)
# you can access its different data slots with @
aotus@cids   # cluster IDs
aotus@sids   # sequence IDs
aotus@txids  # taxonomic IDs
aotus@clstrs # clusters archive
aotus@sqs    # sequence archive
aotus@txdct  # taxonomic dictionary
# see all of the available slots
(slotNames(aotus))
# access different sequences and clusters with [[
(aotus[['0']])              # cluster record 0
(aotus[[aotus@sids[[1]]]])  # first sequence record
# get a summary of the whole object
(summary(aotus))
# the above generates a data.frame with information on each cluster:
# ID - unique id in the object
# Type - cluster type
# Seed - most connected sequence
# Parent - MRCA of all represented taxa
# N_taxa - number of NCBI recognised taxa
# N_seqs - number of sequences
# Med_sql - median sequence length
# MAD - Maximum alignment density, values close to 1 indicate all sequences in
#  the cluster have a similar length.
# Definition - most common words (and frequency) in sequence definitions
# Feature - most common feature name (and frequency)
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage1.R
\name{taxdict_gen}
\alias{taxdict_gen}
\title{Generate taxonomic dictionary}
\usage{
taxdict_gen(txids, recs, ps)
}
\arguments{
\item{txids}{Vector of taxonomic IDs}

\item{recs}{List of taxonomic records}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
TaxDict
}
\description{
Takes a vector of txids and a list
of taxonomic records and returns a taxonomic dictionary.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2-tools.R
\name{seqrec_augment}
\alias{seqrec_augment}
\title{Augment sequence records list}
\usage{
seqrec_augment(sqs, txdct)
}
\arguments{
\item{sqs}{List of SeqRecs}

\item{txdct}{Taxonomic Dictionary}
}
\value{
SeqArc
}
\description{
Add taxids to records and convert to archive.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotaR.R
\name{list_clstrrec_slots}
\alias{list_clstrrec_slots}
\title{List all ClstrRec slots}
\usage{
list_clstrrec_slots()
}
\value{
vector
}
\description{
Returns a vector of all available ClstrRec slots of type
character, integer and numeric.
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_tx_slot}
\alias{get_tx_slot}
\title{Get slot data for each taxon record}
\usage{
get_tx_slot(phylota, txid, slt_nm = list_taxrec_slots())
}
\arguments{
\item{phylota}{Phylota object}

\item{txid}{Taxonomic ID}

\item{slt_nm}{Slot name}
}
\value{
vector or list
}
\description{
Get slot data for taxa(s)
}
\examples{
data('aotus')
random_txid <- sample(aotus@txids, 1)
(get_tx_slot(phylota = aotus, txid = random_txid, slt_nm = 'scnm'))
# see list_taxrec_slots() for available slots
(list_taxrec_slots())
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstr_subtree}
\alias{clstr_subtree}
\title{Cluster all sequences descending from a txid}
\usage{
clstr_subtree(txid, sqs, txdct, dds, ps, lvl)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{sqs}{Sequence object of all downloaded sequences}

\item{txdct}{Taxonomic dictionary}

\item{dds}{Vector of direct descendants}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
ClstrArc
}
\description{
Identifies clusters from sequences associated with a
txid and all its descendants. Clusters returned by this function
will thus be of cl_type 'subtree'.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-api.R
\name{download_obj_check}
\alias{download_obj_check}
\title{Check an object returned from rentrez function}
\usage{
download_obj_check(obj)
}
\arguments{
\item{obj}{Object returned from rentrez function}
}
\value{
T/F
}
\description{
Returns T/F. Checks if object returned from rentrez
function is as expected.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{error}}, \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{blast_sqs}
\alias{blast_sqs}
\title{BLAST All vs All}
\usage{
blast_sqs(txid, typ, sqs, ps, lvl)
}
\arguments{
\item{txid}{Taxonomic node ID, numeric}

\item{typ}{Cluster type, 'direct' or 'subtree'}

\item{sqs}{Sequences}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
blast_res data.frame or NULL
}
\description{
Return BLAST results from BLASTing all vs all for
given sequences. Returns NULL if no BLAST results generated.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-api.R
\name{search_and_cache}
\alias{search_and_cache}
\title{Run rentrez function and cache results}
\usage{
search_and_cache(func, args, fnm, ps)
}
\arguments{
\item{func}{rentrez function}

\item{args}{rentrez function arguments, list}

\item{fnm}{rentrez function name}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
rentrez function results
}
\description{
Safely run a rentrez function. If the query fails, the
function will retry. All query results are cached. To remove cached
data use hard reset.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{run}
\alias{run}
\title{Run phylotaR pipeline}
\usage{
run(wd, nstages = 4)
}
\arguments{
\item{wd}{Working directory}

\item{nstages}{Number of total stages to run, max 4.}
}
\description{
Run the entire phylotaR pipeline. All generated files will be
stored in the wd. The process can be stopped at anytime and  restarted with
\code{restart}. \code{nstages} must be a numeric value representing the
number of stages that will be run. Stages are run in the following order:
1 - taxise, 2 - download, 3 - cluster and 4 - cluster2.

For example, specifying \code{nstages} = 3, will run taxise, download and
cluster. Stages can also be run individually, see linked functions below.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  # e.g. "/usr/local/ncbi/blast/bin/"
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  run(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{setup}},
  \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{blastcache_load}
\alias{blastcache_load}
\title{Load BLAST results from cache}
\usage{
blastcache_load(sids, wd)
}
\arguments{
\item{sids}{Sequence IDs}

\item{wd}{Working dir}
}
\value{
blast_res data.frame or NULL
}
\description{
Run to load cached BLAST results.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{is_txid_in_sq}
\alias{is_txid_in_sq}
\title{Is txid in sequence?}
\usage{
is_txid_in_sq(phylota, txid, sid)
}
\arguments{
\item{phylota}{Phylota}

\item{txid}{Taxonomic ID}

\item{sid}{Sequence ID}
}
\value{
boolean
}
\description{
Checks if given txid is represented by sequence by
looking at sequence source organism's lineage.
}
\examples{
data(tinamous)
sid <- tinamous@sids[[1]]
sq <- tinamous[[sid]]
txid <- sq@txid
# expect true
is_txid_in_sq(phylota = tinamous, txid = txid, sid = sid)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4-tools.R
\name{clstrs_join}
\alias{clstrs_join}
\title{Join clusters for merging}
\usage{
clstrs_join(blast_res, seed_ids, all_clstrs, ps)
}
\arguments{
\item{blast_res}{Seed sequence BLAST results}

\item{seed_ids}{Seed sequence IDs}

\item{all_clstrs}{List of all clusters}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
list of joined clusters
}
\description{
Uses seed sequence BLAST results and IDs to join clusters
identified as sisters into single clusters. Resulting object is of joined
clusters, merging is required to reformat the clusters for subsequent
analysis.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-txdct.R
\name{descendants_get}
\alias{descendants_get}
\title{Get descendants}
\usage{
descendants_get(id, txdct, direct = FALSE)
}
\arguments{
\item{id}{txid}

\item{txdct}{TaxDict}

\item{direct}{T/F, return only direct descendants?}
}
\value{
vector
}
\description{
Look-up either direct or all taxonomic descendants of
a node from taxonomic dictionary.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{download_obj_check}},
  \code{\link{error}}, \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-log.R
\name{warn}
\alias{warn}
\title{Write warning message to log}
\usage{
warn(ps, ...)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}

\item{...}{Message elements for concatenating}
}
\description{
Inform a user if a potential error has occurred in log.txt.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-sqs.R
\name{rawseqrec_breakdown}
\alias{rawseqrec_breakdown}
\title{Breakdown a sequence record into its features}
\usage{
rawseqrec_breakdown(record_parts, ps)
}
\arguments{
\item{record_parts}{list of record elements from a GenBank record}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
list of SeqRecs
}
\description{
Takes GenBank record's elements and returns a SeqRec. For
sequences with lots of features, the sequence is broken down into these
features provided they are of the right size. Sequences are either returned
as features or whole sequence records, never both.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{sids_load}
\alias{sids_load}
\title{Load sids from cache}
\usage{
sids_load(wd, txid)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Taxonomic ID, numeric}
}
\value{
vector of sids
}
\description{
Load sids downloaded by \code{sids_get} function.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_save}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{plot_phylota_pa}
\alias{plot_phylota_pa}
\title{Plot presence/absence matrix}
\usage{
plot_phylota_pa(phylota, cids, txids, cnms = cids, txnms = txids)
}
\arguments{
\item{phylota}{Phylota object}

\item{cids}{Vector of cluster IDs}

\item{txids}{Vector of taxonomic IDs}

\item{cnms}{Cluster names}

\item{txnms}{Taxonomic names}
}
\value{
geom_object
}
\description{
Plot presence/absence of taxa by each
cluster in phylota object.
}
\details{
Cluster names and taxonomic names can be given to the function, by
default IDs are used.
}
\examples{
library(phylotaR)
data(cycads)
# drop all but first ten
cycads <- drop_clstrs(cycads, cycads@cids[1:10])
# plot all
p <- plot_phylota_pa(phylota = cycads, cids = cycads@cids, txids = cycads@txids)
print(p)  # lots of information, difficult to interpret
# get genus-level taxonomic names
genus_txids <- get_txids(cycads, txids = cycads@txids, rnk = 'genus')
genus_txids <- unique(genus_txids)
# dropping missing
genus_txids <- genus_txids[genus_txids !=  '']
genus_nms <- get_tx_slot(cycads, genus_txids, slt_nm = 'scnm')
# make alphabetical for plotting
genus_nms <- sort(genus_nms, decreasing = TRUE)
# generate geom_object
p <- plot_phylota_pa(phylota = cycads, cids = cycads@cids, txids = genus_txids,
                     txnms = genus_nms)
# plot
print(p)  # easier to interpret
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{sids_save}
\alias{sids_save}
\title{Save sids to cache}
\usage{
sids_save(wd, txid, sids)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Taxonomic ID, numeric}

\item{sids}{sids}
}
\description{
Saves sids downloaded
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-sqs.R
\name{seqrec_gen}
\alias{seqrec_gen}
\title{Generate sequence record}
\usage{
seqrec_gen(accssn, nm, txid, sq, dfln, orgnsm, ml_typ, rec_typ, vrsn, age,
  lctn = NULL)
}
\arguments{
\item{accssn}{Accession ID}

\item{nm}{Sequence name}

\item{txid}{Taxonomic ID of source organism}

\item{sq}{Sequence}

\item{dfln}{Definition line}

\item{orgnsm}{Source organism name}

\item{ml_typ}{Molecule type}

\item{rec_typ}{Sequence record type}

\item{vrsn}{Accession version}

\item{age}{Number of days since upload}

\item{lctn}{Location numbers for features, e.g. '1..200'}
}
\value{
SeqRec
}
\description{
Creates an S4 SeqRec
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_get}},
  \code{\link{sids_check}}, \code{\link{sids_get}},
  \code{\link{sids_load}}, \code{\link{sids_save}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-txdct.R
\name{rank_get}
\alias{rank_get}
\title{Get rank}
\usage{
rank_get(txid, txdct)
}
\arguments{
\item{txid}{txid}

\item{txdct}{TaxDict}
}
\value{
character
}
\description{
Look-up taxonomic rank from dictionary.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3.R
\name{clstrs_calc}
\alias{clstrs_calc}
\title{Calculate clusters for all sequences in wd}
\usage{
clstrs_calc(txdct, ps)
}
\arguments{
\item{txdct}{Taxonomic dictionary}

\item{ps}{Parameters list, generated with parameters()}
}
\description{
Loop through downloaded sequences for each clade and
hierarchically find clusters using BLAST.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotaR.R
\name{parameters}
\alias{parameters}
\title{Default parameters}
\usage{
parameters(wd = ".", txid = character(), mkblstdb = "", blstn = "",
  v = FALSE, ncps = 1, mxnds = 1e+05, mdlthrs = 3000,
  mnsql = 250, mxsql = 2000, mxrtry = 100, mxsqs = 50000,
  mxevl = 1e-10, mncvrg = 51, btchsz = 100, db_only = FALSE,
  outsider = FALSE,
  srch_trm = "NOT predicted[TI] NOT \\"whole genome shotgun\\"[TI] NOT unverified[TI] NOT \\"synthetic construct\\"[Organism] NOT refseq[filter] NOT TSA[Keyword]",
  date = Sys.Date())
}
\arguments{
\item{wd}{The working directory where all output files are saved.}

\item{txid}{Taxonomic group of interest, allows vectors.}

\item{mkblstdb}{File path to makeblastdb}

\item{blstn}{File path to blastn}

\item{v}{Print progress statements to console?
Statements will always be printed to log.txt.}

\item{ncps}{The number of threads to use in the local-alignment search tool.}

\item{mxnds}{The maximum number of nodes descending from a taxonomic group.
If there are more than this number, nodes at the lower taxonomic level are
analysed.}

\item{mdlthrs}{'Model organism threshold'. Taxa with more sequences than this
number will be considered model organisms and a random mdlthrs subset of
their sequences will be downloaded.}

\item{mnsql}{The minimum length of sequence in nucleotide base pairs to
download.}

\item{mxsql}{The maximum length of sequence in nucleotide base pairs to
download. Any longer sequences will be ignored.}

\item{mxrtry}{The maximum number of attempts to make when downloading.}

\item{mxsqs}{The maximum number of sequences to BLAST in all-vs-all searches.
If there are more sequences for a node, BLAST is performed at the lower
taxonomic level.}

\item{mxevl}{The maximum E-value for a successful BLAST.}

\item{mncvrg}{The maximum percentile coverage defining an overlapping BLAST
hit. Sequences with BLAST matches with lower values are not considered
orthologous.}

\item{btchsz}{Batch size when querying NCBI}

\item{db_only}{Take sequences only from \code{restez} library? TRUE/FALSE.
If TRUE, internet is still required (for taxonomic look-up) and a
\code{restez} will need to be set up.}

\item{outsider}{Use \code{om..blast}? TRUE/FALSE. If TRUE, a module for
running \code{blastn} will be installed and all BLAST commands will be run
through it. \code{outsider} package is required.}

\item{srch_trm}{Sequence NCBI search term modifier. Use this parameter to
change the default search term options. Default: avoid predicted, WGS,
unverified, synthetic, RefSeq and Transcriptome Shotgun Assembly sequences.}

\item{date}{Date when pipeline was initiated}
}
\value{
list
}
\description{
Returns a parameter list with default parameter values.
}
\details{
This function is NOT used to change the parameters in a folder.
Use parameters_reset() instead. The purpose of this function is to describe
the paramaters and present their default values.
}
\concept{public-pipeline}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstrarc_join}
\alias{clstrarc_join}
\title{Join two cluster archive}
\usage{
clstrarc_join(clstrarc_1, clstrarc_2)
}
\arguments{
\item{clstrarc_1}{ClstrArc}

\item{clstrarc_2}{ClstrArc}
}
\value{
ClstrArc
}
\description{
Take two ClstrArc classes and join them into a single
ClstrArc.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline-tools.R
\name{stage_args_check}
\alias{stage_args_check}
\title{Check stage arguments}
\usage{
stage_args_check(to, frm)
}
\arguments{
\item{to}{ending stage}

\item{frm}{starting stage}
}
\value{
character, stage message
}
\description{
Ensures stage arguments are valid, raises an error if not.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3.R
\name{clusters_run}
\alias{clusters_run}
\title{Run the cluster stage}
\usage{
clusters_run(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Run the third stage of the phylotaR pipeline, cluster.
This stage hierarchically traverses the taxonomy identifying all
direct and subtree clusters from downloaded sequences. Any
taxonomic nodes too small for cluster identification are placed
into paraphyletic clusters.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run stages
  taxise_run(wd = wd)
  download_run(wd = wd)
  clusters_run(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-entrez.R
\name{txnds_count}
\alias{txnds_count}
\title{Count number of descending taxonomic nodes}
\usage{
txnds_count(txid, ps)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
integer
}
\description{
Searches NCBI taxonomy and returns number of descendants
taxonomic nodes (species, genera ...) of ID.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{reset}
\alias{reset}
\title{Reset a phylotaR pipeline run}
\usage{
reset(wd, stage, hard = FALSE)
}
\arguments{
\item{wd}{Working directory}

\item{stage}{Name of stage to which the pipeline will be reset}

\item{hard}{T/F, delete all cached data?}
}
\description{
Resets the pipeline to a specified stage.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run taxise
  taxise_run(wd = wd)
  # reset back to taxise as if it has not been run
  reset(wd = 'aotus', stage = 'taxise')
  # run taxise again ....
  taxise_run(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{restart}},
  \code{\link{run}}, \code{\link{setup}},
  \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstrrec_gen}
\alias{clstrrec_gen}
\title{Generate list of clusters}
\usage{
clstrrec_gen(clstr_list, txid, sqs, typ)
}
\arguments{
\item{clstr_list}{List of list of cluster descriptions}

\item{txid}{Taxonomic node ID}

\item{sqs}{Sequence records}

\item{typ}{Cluster type}
}
\value{
list of ClstrRecs
}
\description{
Takes a list of lists of cluster descriptions
and generates ClstrRecs.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{clstrs_save}
\alias{clstrs_save}
\title{Save clusters to cache}
\usage{
clstrs_save(wd, txid, clstrs)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Taxonomic ID, numeric}

\item{clstrs}{cluster list}
}
\description{
Saves clusters generated by \code{clstr_all} to cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{cmdln}},
  \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-txdct.R
\name{parent_get}
\alias{parent_get}
\title{Get taxonomic parent}
\usage{
parent_get(id, txdct)
}
\arguments{
\item{id}{txid(s)}

\item{txdct}{TaxDict}
}
\value{
Character
}
\description{
Look-up MRCA of taxonomic id(s) from taxonomic
dictionary
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_nsqs}
\alias{get_nsqs}
\title{Count number of sequences}
\usage{
get_nsqs(phylota, cid)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s)}
}
\value{
vector
}
\description{
Count the number of sequences in a cluster(s).
}
\examples{
data("cycads")
# count seqs for a random 10 clusters
random_cids <- sample(cycads@cids, 10)
nsqs <- get_nsqs(phylota = cycads, cid = random_cids)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-entrez.R
\name{sids_get}
\alias{sids_get}
\title{Return random set of sequence IDs}
\usage{
sids_get(txid, direct, ps, retmax = 100, hrdmx = 1e+05)
}
\arguments{
\item{txid}{NCBI taxon identifier}

\item{direct}{Node-level only or subtree as well? Default FALSE.}

\item{ps}{Parameters list, generated with parameters()}

\item{retmax}{Maximum number of sequences when querying model
organisms. The smaller the more random, the larger the faster.}

\item{hrdmx}{Absolute maximum number of sequence IDs to download
in a single query.}
}
\value{
vector of IDs
}
\description{
For a given txid return a random set of 
sequences associated.
}
\details{
For model organisms downloading all IDs can a take long
time or even cause an xml parsing error. For any search with more
than hrdmx sequences, this function we will run multiple small
searches downloading retmax seq IDs at a time with different
retstart values to generate a semi-random vector of sequence IDs.
For all other searches, all IDs will be retrieved. Note, it makes
no sense for mdlthrs in parameters to be greater than hrdmx in this
function.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_load}}, \code{\link{sids_save}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-sqs.R
\name{seqarc_gen}
\alias{seqarc_gen}
\title{Generate sequence archive}
\usage{
seqarc_gen(seqrecs)
}
\arguments{
\item{seqrecs}{List of SeqRecs}
}
\value{
SeqArc
}
\description{
Creates an S4 SeqArc from list of SeqRecs
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-sqs.R
\name{gb_extract}
\alias{gb_extract}
\title{Extract elements from a raw GenBank record}
\usage{
gb_extract(record)
}
\arguments{
\item{record}{raw GenBank text record}
}
\value{
list of GenBank elements
}
\description{
Returns a list of elements from a GenBank record such as
'organism', 'sequence' and features.
}
\details{
Uses restez extract functions. See restez package for more details.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4-tools.R
\name{clstrs_renumber}
\alias{clstrs_renumber}
\title{Renumber cluster IDs}
\usage{
clstrs_renumber(clstrrecs)
}
\arguments{
\item{clstrrecs}{List of clusters}
}
\value{
ClstrArc
}
\description{
Returns a ClstrArc with ID determined by the number
of sequences in each cluster.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4-tools.R
\name{seeds_blast}
\alias{seeds_blast}
\title{BLAST seed sequences}
\usage{
seeds_blast(sqs, ps)
}
\arguments{
\item{sqs}{All seed sequences to be BLASTed}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
blast res data.frame
}
\description{
Runs all-v-all blast for seed sequences.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seq_download}},
  \code{\link{seqarc_gen}}, \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{obj_load}
\alias{obj_load}
\title{Load a named object from the cache}
\usage{
obj_load(wd, nm)
}
\arguments{
\item{wd}{Working directory}

\item{nm}{Object name}
}
\value{
object, multiple formats possible
}
\description{
Loads an object from the cache as stored by
\code{obj_save}.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_save}}, \code{\link{outfmt_get}},
  \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{plot_phylota_treemap}
\alias{plot_phylota_treemap}
\title{Plot treemap of Phylota object}
\usage{
plot_phylota_treemap(phylota, cids = NULL, txids = NULL, cnms = cids,
  txnms = txids, with_labels = TRUE, area = c("ntx", "nsq", "ncl"),
  fill = c("NULL", "typ", "ntx", "nsq", "ncl"))
}
\arguments{
\item{phylota}{Phylota object}

\item{cids}{Cluster IDs}

\item{txids}{Taxonomic IDs}

\item{cnms}{Cluster names}

\item{txnms}{Taxonomic names}

\item{with_labels}{Show names per box?}

\item{area}{What determines the size per box?}

\item{fill}{What determines the coloured fill per box?}
}
\value{
geom_object
}
\description{
Treemaps show relative size with boxes. The user can
explore which taxa or clusters are most represented either by
sequence or cluster number. If cluster IDs are provided, the plot
is made for clusters. If taxonomic IDs are provided, the plot is
made for taxa.
}
\details{
The function can take a long time to run for large Phylota
objects over many taxonomic IDs because searches are made across
lineages. The idea of the function is to assess the data dominance
of specific clusters and taxa.
}
\examples{
data("tinamous")
# Plot clusters, size by n. sq, fill by n. tx
p <- plot_phylota_treemap(phylota = tinamous, cids = tinamous@cids,
                          area = 'nsq', fill = 'ntx')
print(p)
# Plot taxa, size by n. sq, fill by ncl
txids <- get_txids(tinamous, txids = tinamous@txids, rnk = 'genus')
txids <- txids[txids !=  '']
txids <- unique(txids)
txnms <- get_tx_slot(tinamous, txids, slt_nm = 'scnm')
p <- plot_phylota_treemap(phylota = tinamous, txids = txids, txnms = txnms,
                          area = 'nsq', fill = 'ncl')
print(p)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4.R
\name{clstr2_calc}
\alias{clstr2_calc}
\title{Cluster sets of clusters identified in cluster stage}
\usage{
clstr2_calc(ps)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}
}
\description{
Loads cluster sets from cache. Extracts seed sequences
and runs all-v-all BLAST of seeds to identify sister clusters.
Sisters are then merged. An object of all sequences and clusters
is then saved in cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-api.R
\name{safely_connect}
\alias{safely_connect}
\title{Safely run rentrez function}
\usage{
safely_connect(func, args, fnm, ps)
}
\arguments{
\item{func}{rentrez function}

\item{args}{rentrez function arguments, list}

\item{fnm}{rentrez function name}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
rentrez function results
}
\description{
Safely run a rentrez function. If the query fails,
the function will retry.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-blast.R
\name{blastdb_gen}
\alias{blastdb_gen}
\title{Generate a BLAST database}
\usage{
blastdb_gen(sqs, dbfl, ps)
}
\arguments{
\item{sqs}{Sequences}

\item{dbfl}{Outfile for database}

\item{ps}{Parameters list, generated with parameters()}
}
\description{
Generate BLAST database in wd for given sequences.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastn_run}},
  \code{\link{cache_rm}}, \code{\link{cache_setup}},
  \code{\link{clade_select}}, \code{\link{clstr2_calc}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage1-tools.R
\name{tax_download}
\alias{tax_download}
\title{Download taxonomic records}
\usage{
tax_download(ids, ps)
}
\arguments{
\item{ids}{Vector of taxonomic IDs}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
list of list
}
\description{
Downloads one batch of taxonomic
records.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{ncbicache_load}
\alias{ncbicache_load}
\title{Retrieve cached NCBI query}
\usage{
ncbicache_load(fnm, args, wd)
}
\arguments{
\item{fnm}{NCBI Entrez function name}

\item{args}{Args used for function}

\item{wd}{Working directory}
}
\value{
rentrez result
}
\description{
Run this function to load cached NCBI queries.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline-tools.R
\name{blast_setup}
\alias{blast_setup}
\title{Ensures NCBI BLAST tools are installed}
\usage{
blast_setup(d, v, wd, otsdr)
}
\arguments{
\item{d}{Directory to NCBI BLAST tools}

\item{v}{v, T/F}

\item{wd}{Working directory}

\item{otsdr}{Run through \code{outsider}?}
}
\value{
list
}
\description{
Ensures NCBI BLAST executables are installed on the system. Tests
version number of BLAST tools.
}
\details{
BLAST tools must be version >= 2.0
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_sqs}}, \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-txdct.R
\name{taxtree_gen}
\alias{taxtree_gen}
\title{Generate taxonomic tree}
\usage{
taxtree_gen(prinds, ids, root, ps)
}
\arguments{
\item{prinds}{Vector of integers indicating preceding node.}

\item{ids}{Vector of taxonomic IDs}

\item{root}{ID of root taxon}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
TreeMan

TreeMan class
}
\description{
Generate a taxonomic tree for
easy look up of taxonomic parents and descendants.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4.R
\name{clusters2_run}
\alias{clusters2_run}
\title{Run the cluster2 stage}
\usage{
clusters2_run(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Run the fourth stage of the phylotaR pipeline,
cluster2. Identify clusters at higher taxonomic levels by
merging sister clusters.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run stages
  taxise_run(wd = wd)
  download_run(wd = wd)
  clusters_run(wd = wd)
  clusters2_run(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{restart}
\alias{restart}
\title{Restart a phylotaR pipeline run}
\usage{
restart(wd, nstages = 4)
}
\arguments{
\item{wd}{Working directory}

\item{nstages}{Number of total stages to run, max 4.}
}
\description{
Restarts the running of a pipeline
as started with \code{run}.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # run and stop after 10 seconds
  R.utils::withTimeout(expr = {
    run(wd = wd)
  }, timeout = 10)
  # use ctrl+c or Esc to kill without a timelimit
  # and restart with ....
  restart(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{run}}, \code{\link{setup}},
  \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-entrez.R
\name{sqs_count}
\alias{sqs_count}
\title{Count number of sequences for txid}
\usage{
sqs_count(txid, ps, direct = FALSE)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{ps}{Parameters list, generated with parameters()}

\item{direct}{Node-level only or subtree as well? Default FALSE.}
}
\value{
integer
}
\description{
Return the number of sequences associated with a
taxonomic ID on NCBI GenBank.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotaR.R
\name{list_ncbi_ranks}
\alias{list_ncbi_ranks}
\title{List all NCBI Ranks}
\usage{
list_ncbi_ranks()
}
\value{
vector
}
\description{
Returns a vector of all
NCBI taxonomic ranks in descending order.
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{TaxDict-class}
\alias{TaxDict-class}
\alias{TaxDict-method}
\alias{as.character,TaxDict-method}
\alias{show,TaxDict-method}
\alias{print,TaxDict-method}
\alias{str,TaxDict-method}
\alias{summary,TaxDict-method}
\title{Taxonomic record dictionary}
\usage{
\S4method{as.character}{TaxDict}(x)

\S4method{show}{TaxDict}(object)

\S4method{print}{TaxDict}(x)

\S4method{str}{TaxDict}(object, max.level = 2L, ...)

\S4method{summary}{TaxDict}(object)
}
\arguments{
\item{x}{\code{TaxDict} object}

\item{object}{\code{TaxDict} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}
}
\description{
Taxonomic dictionary contains a taxonomic
tree and NCBI taxonomy data for all taxonomic IDs.
}
\section{Slots}{

\describe{
\item{\code{txids}}{Taxonomic IDs of taxon records}

\item{\code{recs}}{Environment of records}

\item{\code{prnt}}{Parent taxonomic ID}

\item{\code{txtr}}{Taxonomic tree}
}}

\examples{
data('aotus')
txdct <- aotus@txdct
# this is a TaxDict object
# it contains taxonomic information, including records and tree
show(txdct)
# you can access its different data slots with @
txdct@txids  # taxonomic IDs
txdct@recs   # taxonomic records environment
txdct@txtr   # taxonomic tree
txdct@prnt   # MRCA
# access any record through the records environment
txdct@recs[[txdct@txids[[1]]]]
# for interacting with the taxonomic tree, see the treeman package
summary(txdct@txtr)
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxRec-class}},
  \code{\link{clusters2_run}}, \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-drop.R
\name{drop_by_rank}
\alias{drop_by_rank}
\title{Reduce clusters to specific rank}
\usage{
drop_by_rank(phylota, rnk = "species", keep_higher = FALSE, n = 10,
  choose_by = c("pambgs", "age", "nncltds"), greatest = c(FALSE, FALSE,
  TRUE))
}
\arguments{
\item{phylota}{Phylota object}

\item{rnk}{Taxonomic rank}

\item{keep_higher}{Keep higher taxonomic ranks?}

\item{n}{Number of sequences per taxon}

\item{choose_by}{Vector of selection functions}

\item{greatest}{Greatest of lowest for each choose_by function}
}
\value{
phylota
}
\description{
Identifies higher level taxa for each sequence in clusters for
given rank. Selects representative sequences for each unique taxon using the
choose_by functions. By default, the function will choose the top ten
sequences by first sorting by those with fewest number of ambiguous
sequences, then by youngest, then by sequence length.
}
\examples{
data("dragonflies")
# For faster computations, let's only work with the 5 clusters.
dragonflies <- drop_clstrs(phylota = dragonflies, cid = dragonflies@cids[10:15])


# We can use drop_by_rank() to reduce to 10 sequences per genus for each cluster
(reduced_1 <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 10,
                           choose_by = c('pambgs', 'age', 'nncltds'),
                           greatest = c(FALSE, FALSE, TRUE)))

# We can specify what aspects of the sequences we would like to select per genus
# By default we select the sequences with fewest ambiguous nucleotides (e.g.
# we avoid Ns), the youngest age and then longest sequence.
# We can reverse the 'greatest' to get the opposite.
(reduced_2 <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 10,
                           choose_by = c('pambgs', 'age', 'nncltds'),
                           greatest = c(TRUE, TRUE, FALSE)))


# Leading to smaller sequnces ...
r1_sqlngth <- mean(get_sq_slot(phylota = reduced_1,
                                sid = reduced_1@sids, slt_nm = 'nncltds'))
r2_sqlngth <- mean(get_sq_slot(phylota = reduced_2,
                                sid = reduced_2@sids, slt_nm = 'nncltds'))
(r1_sqlngth > r2_sqlngth)
# ... with more ambigous characters ....
r1_pambgs <- mean(get_sq_slot(phylota = reduced_1, sid = reduced_1@sids,
                              slt_nm = 'pambgs'))
r2_pambgs <- mean(get_sq_slot(phylota = reduced_2, sid = reduced_2@sids,
                              slt_nm = 'pambgs'))
(r1_pambgs < r2_pambgs)
# .... and older ages (measured in days since being added to GenBank).
r1_age <- mean(get_sq_slot(phylota = reduced_1, sid = reduced_1@sids,
                           slt_nm = 'age'))
r2_age <- mean(get_sq_slot(phylota = reduced_2, sid = reduced_2@sids,
                           slt_nm = 'age'))
(r1_age < r2_age)


# Or... we can simply reduce the clusters to just one sequence per genus
(dragonflies <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 1))
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_clstrs}},
  \code{\link{drop_sqs}}, \code{\link{get_clstr_slot}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{ClstrRec-class}
\alias{ClstrRec-class}
\alias{ClstrRec-method}
\alias{as.character,ClstrRec-method}
\alias{show,ClstrRec-method}
\alias{print,ClstrRec-method}
\alias{str,ClstrRec-method}
\alias{summary,ClstrRec-method}
\title{Cluster record}
\usage{
\S4method{as.character}{ClstrRec}(x)

\S4method{show}{ClstrRec}(object)

\S4method{print}{ClstrRec}(x)

\S4method{str}{ClstrRec}(object, max.level = 2L, ...)

\S4method{summary}{ClstrRec}(object)
}
\arguments{
\item{x}{\code{ClstrRec} object}

\item{object}{\code{ClstrRec} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}
}
\description{
Cluster record contains all information on a cluster.
}
\section{Slots}{

\describe{
\item{\code{id}}{Cluster ID, integer}

\item{\code{sids}}{Sequence IDs}

\item{\code{nsqs}}{Number of sequences}

\item{\code{txids}}{Source txids for sequences}

\item{\code{ntx}}{Number of taxa}

\item{\code{typ}}{Cluster type: direct, subtree or merged}

\item{\code{seed}}{Seed sequence ID}

\item{\code{prnt}}{Parent taxonomic ID}
}}

\examples{
data('aotus')
clstrrec <- aotus@clstrs@clstrs[[1]]
# this is a ClstrRec object
# it contains cluster information
show(clstrrec)
# you can access its different data slots with @
clstrrec@id     # cluster id
clstrrec@sids   # sequence IDs
clstrrec@nsqs   # number of sequences
clstrrec@txids  # taxonomic IDs of sequences
clstrrec@ntx    # number unique taxonomic IDs
clstrrec@typ    # cluster type: merged, subtree, direct or paraphyly
clstrrec@prnt   # MRCA of all taxa
clstrrec@seed   # most inter-connected sequence
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{blastcache_save}
\alias{blastcache_save}
\title{Save BLAST results to cache}
\usage{
blastcache_save(sids, wd, obj)
}
\arguments{
\item{sids}{Sequence IDs}

\item{wd}{Working dir}

\item{obj}{BLAST result}
}
\description{
Run whenever local BLAST runs are made to save results
in cache in case the pipeline is run again.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{update_phylota}
\alias{update_phylota}
\title{Update slots}
\usage{
update_phylota(phylota)
}
\arguments{
\item{phylota}{Phylota}
}
\value{
Phylota
}
\description{
After change, run to update slots.
}
\seealso{
Other tools-private: \code{\link{mk_txid_in_sq_mtrx}},
  \code{\link{summary_phylota}}
}
\concept{tools-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{sqs_save}
\alias{sqs_save}
\title{Save sequences to cache}
\usage{
sqs_save(wd, txid, sqs)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Taxonomic ID, numeric}

\item{sqs}{Sequences}
}
\description{
Saves sequences downloaded
}
\details{
Used within the \code{dwnld} function. Saves
sequence data by txid in cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{ncbicache_save}
\alias{ncbicache_save}
\title{Save NCBI query result to cache}
\usage{
ncbicache_save(fnm, args, wd, obj)
}
\arguments{
\item{fnm}{NCBI Entrez function name}

\item{args}{Args used for function}

\item{wd}{Working directory}

\item{obj}{NCBI query result}
}
\description{
Run whenever NCBI queries are made to save results in
cache in case the pipeline is run again.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-blast.R
\name{outfmt_get}
\alias{outfmt_get}
\title{Determine 'outformat' format}
\usage{
outfmt_get(ps)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}
}
\value{
character
}
\description{
Depending on operating system, BLAST may or may not require ""
around \code{outfmt}. This function will run a micro BLAST analysis
to test. It will return the \code{outfmt} for use in \code{blastn_run}.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2.R
\name{clade_select}
\alias{clade_select}
\title{Get all node IDs that will be processed}
\usage{
clade_select(txdct, ps)
}
\arguments{
\item{txdct}{TxDct}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
vector of txids
}
\description{
All nodes with less than maximum number
of nodes and sequences.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clstr2_calc}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{progress_reset}
\alias{progress_reset}
\title{Reset progress}
\usage{
progress_reset(wd, stg)
}
\arguments{
\item{wd}{Working directory}

\item{stg}{Stage to which the pipeline will be reset}
}
\description{
Reset progress to an earlier completed stage.
}
\details{
For example, resetting the progress to 'download'
mark stages 'download', 'cluster' and 'cluster2' as un-run.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{progress_init}
\alias{progress_init}
\title{Initialise progress list in cache}
\usage{
progress_init(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Creates a progress list recording each stage run in
cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{summary_phylota}
\alias{summary_phylota}
\title{Summarise clusters in Phylota Table}
\usage{
summary_phylota(phylota)
}
\arguments{
\item{phylota}{Phylota object}
}
\description{
Generates a summary data.frame from all clusters in
Phylota object.
}
\seealso{
Other tools-private: \code{\link{mk_txid_in_sq_mtrx}},
  \code{\link{update_phylota}}
}
\concept{tools-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-drop.R
\name{drop_clstrs}
\alias{drop_clstrs}
\title{Drop cluster records from phylota object}
\usage{
drop_clstrs(phylota, cid)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s) to be kept}
}
\value{
phylota
}
\description{
Drops all clusters except those
identified by user.
}
\examples{
data("dragonflies")
# specify cids to *keep*
random_cids <- sample(dragonflies@cids, 100)
# drop an entire cluster
nbefore <- length(dragonflies@cids)
dragonflies <- drop_clstrs(phylota = dragonflies, cid = random_cids)
nafter <- length(dragonflies@cids)
# now there are only 100 clusters
(nafter < nbefore)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_sqs}}, \code{\link{get_clstr_slot}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{obj_check}
\alias{obj_check}
\title{Check if an object exists}
\usage{
obj_check(wd, nm)
}
\arguments{
\item{wd}{Working directory}

\item{nm}{Object name}
}
\value{
T/F
}
\description{
Check if an object exists in the cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_load}},
  \code{\link{obj_save}}, \code{\link{outfmt_get}},
  \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2-tools.R
\name{hierarchic_download}
\alias{hierarchic_download}
\title{Hierarchically get sequences for a txid}
\usage{
hierarchic_download(txid, txdct, ps, lvl = 0)
}
\arguments{
\item{txid}{Taxonomic node ID, numeric}

\item{txdct}{Taxonomic dictionary}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
Vector of SeqRecs
}
\description{
Looks up and downloads sequences for a taxonomic ID.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_ntaxa}
\alias{get_ntaxa}
\title{Count number of unique taxa}
\usage{
get_ntaxa(phylota, cid = NULL, sid = NULL, rnk = NULL,
  keep_higher = FALSE)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s)}

\item{sid}{Sequence ID(s)}

\item{rnk}{Taxonomic rank}

\item{keep_higher}{Keep higher taxonomic ranks?}
}
\value{
vector
}
\description{
Count the number of unique taxa represented by cluster(s) or
sequences in phylota table Use rnk to specify a taxonomic level to count. If
NULL counts will be made to the lowest level reported on NCBI.
}
\examples{
data('bromeliads')
# how many species are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'species'))
# how many genera are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'genus'))
# how many families are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'family'))
# use list_ncbi_ranks() to see available rank names
(list_ncbi_ranks())
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{sids_check}
\alias{sids_check}
\title{Check if sids exist}
\usage{
sids_check(wd, txid)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Taxonomic ID, numeric}
}
\value{
T/F
}
\description{
Check if sids are already downloaded for a txid.
}
\details{
#' @name sqs_load
#' @title Load sequences from cache
#' @description Load sequences downloaded by \code{dwnld} function.
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @family run-private
#' @return SeqArc
sqs_load <- function(wd, txid) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'sqs')
  if (!file.exists(d)) {
    stop('`sqs` not in cache. Have you run the download stage?')
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  if (!file.exists(fl)) {
    stop(paste0('[', txid, '] not in `sqs` of cache.'))
  }
  sqs <- try(readRDS(file = fl), silent = TRUE)
  if (inherits(sqs, 'try-error')) {
    file.remove(fl)
  }
  sqs
}
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_get}},
  \code{\link{sids_load}}, \code{\link{sids_save}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{read_phylota}
\alias{read_phylota}
\title{Generate a Phylota object in R}
\usage{
read_phylota(wd)
}
\arguments{
\item{wd}{Working directory}
}
\value{
Phylota
}
\description{
Creates a Phylota object containing information on
clusters, sequences and taxonomy from the working directory of a
completed pipeline.
}
\examples{
\dontrun{
  
  # Note, this example requires a wd with a completed phylotaR run
  phylota <- read_phylota(wd)
}
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline-tools.R
\name{parameters_setup}
\alias{parameters_setup}
\title{Set Up Parameters}
\usage{
parameters_setup(wd, ncbi_execs, overwrite = FALSE, ...)
}
\arguments{
\item{wd}{Working directory}

\item{ncbi_execs}{File directories for NCBI tools, see \code{blast_setup()}}

\item{overwrite}{Overwrite existing cache?}

\item{...}{Set parameters, see parameters()}
}
\description{
Initiates cache of parameters.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parent_get}}, \code{\link{progress_init}},
  \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-drop.R
\name{drop_sqs}
\alias{drop_sqs}
\title{Drop sequences in a cluster}
\usage{
drop_sqs(phylota, cid, sid)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID}

\item{sid}{Sequence ID(s) to be kept}
}
\value{
phylota
}
\description{
Drop all sequences in a cluster except those identified by user.
}
\examples{
data("dragonflies")
# drop random sequences from cluster 0
clstr <- dragonflies[['0']]
# specify the sids to *keep*
sids <- sample(clstr@sids, 100)
(dragonflies <- drop_sqs(phylota = dragonflies, cid = '0', sid = sids))
# Note, sequences dropped may be represented in other clusters
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{get_clstr_slot}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage1.R
\name{txids_get}
\alias{txids_get}
\title{Searches for descendant taxonomic IDs}
\usage{
txids_get(ps, retmax = 10000)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}

\item{retmax}{integer, maximum number of IDs to return per query}
}
\value{
Vector of txids

vector of ids
}
\description{
Searches NCBI taxonomy for all descendant taxonomic nodes.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{cache_rm}
\alias{cache_rm}
\title{Delete a cache}
\usage{
cache_rm(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Deletes a cache from a wd.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_setup}},
  \code{\link{clade_select}}, \code{\link{clstr2_calc}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_stage_times}
\alias{get_stage_times}
\title{Get run times for different stages}
\usage{
get_stage_times(wd)
}
\arguments{
\item{wd}{Working directory}
}
\value{
list of runtimes in minutes
}
\description{
Get slot data for taxa(s)
}
\examples{
\dontrun{
  
  # Note, this example requires a wd with a completed phylotaR run
  # return a named list of the time take in minutes for each stage
  get_stage_times(wd = wd)
}
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-log.R
\name{error}
\alias{error}
\title{Write error message to log}
\usage{
error(ps, ...)
}
\arguments{
\item{ps}{Parameters list, generated with parameters()}

\item{...}{Message elements for concatenating}
}
\description{
Inform a user if an error has occurred in log.txt,
halt pipeline.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_txids}
\alias{get_txids}
\title{Get taxonomic IDs by rank}
\usage{
get_txids(phylota, cid = NULL, sid = NULL, txids = NULL,
  rnk = NULL, keep_higher = FALSE)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID}

\item{sid}{Sequence ID(s)}

\item{txids}{Vector of txids}

\item{rnk}{Taxonomic rank}

\item{keep_higher}{Keep higher taxonomic IDs?}
}
\value{
vector
}
\description{
Return taxonomic IDs for a vector of sequence IDs or all
sequences in a cluster. User can specify what rank the IDs should be
returned. If NULL, the lowest level is returned.
}
\details{
txids can either be provided by user or they can be determined for
a vector of sids or for a cid. If keep_higher is TRUE, any sequence that has
a identity that is higher than the given rank will be returned. If FALSE,
these sequences will return ''.
}
\examples{
data('bromeliads')
# get all the genus IDs and names
genus_ids <- get_txids(phylota = bromeliads, txids = bromeliads@txids,
                       rnk = 'genus')
genus_ids <- unique(genus_ids)
# drop empty IDs -- this happens if a given lineage has no ID for specified rank
genus_ids <- genus_ids[genus_ids != '']
# get names
(get_tx_slot(phylota = bromeliads, txid = genus_ids, slt_nm = 'scnm'))
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstr_sqs}
\alias{clstr_sqs}
\title{Identify clusters from sequences}
\usage{
clstr_sqs(txid, sqs, ps, lvl, typ = c("direct", "subtree", "paraphyly"))
}
\arguments{
\item{txid}{Taxonomic ID}

\item{sqs}{Sequence object of sequences to be BLASTed}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}

\item{typ}{Direct, subtree or paraphyly?}
}
\description{
Given a sequence object, this function will generate
a list of cluster objects using BLAST
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotaR.R
\name{list_taxrec_slots}
\alias{list_taxrec_slots}
\title{List all TaxRec slots}
\usage{
list_taxrec_slots()
}
\value{
vector
}
\description{
Returns a vector of all available TaxRec slots of type
character, integer and numeric.
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{progress_read}
\alias{progress_read}
\title{Read the progress from cache}
\usage{
progress_read(wd)
}
\arguments{
\item{wd}{Working directory}
}
\value{
stage name, character, or NA is complete
}
\description{
Return the last completed stage using the cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{SeqRec-class}
\alias{SeqRec-class}
\alias{SeqRec-method}
\alias{as.character,SeqRec-method}
\alias{show,SeqRec-method}
\alias{print,SeqRec-method}
\alias{str,SeqRec-method}
\alias{summary,SeqRec-method}
\title{Sequence record}
\usage{
\S4method{as.character}{SeqRec}(x)

\S4method{show}{SeqRec}(object)

\S4method{print}{SeqRec}(x)

\S4method{str}{SeqRec}(object, max.level = 2L, ...)

\S4method{summary}{SeqRec}(object)
}
\arguments{
\item{x}{\code{SeqRec} object}

\item{object}{\code{SeqRec} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}
}
\description{
Sequence record contains sequence data.
}
\details{
Sequence is stored as raw. Use rawToChar().
}
\section{Slots}{

\describe{
\item{\code{id}}{Unique ID}

\item{\code{nm}}{Best-guess sequence name}

\item{\code{accssn}}{Accession}

\item{\code{vrsn}}{Accession version}

\item{\code{url}}{URL}

\item{\code{txid}}{Taxonomic ID of source taxon}

\item{\code{orgnsm}}{Scientific name of source taxon}

\item{\code{sq}}{Sequence}

\item{\code{dfln}}{Definition line}

\item{\code{ml_typ}}{Molecule type, e.g. DNA}

\item{\code{rec_typ}}{Record type: Whole or feature}

\item{\code{nncltds}}{Number of nucleotides}

\item{\code{nambgs}}{Number of ambiguous nucleotides}

\item{\code{pambgs}}{Proportion of ambiguous nucleotides}

\item{\code{gcr}}{GC ratio}

\item{\code{age}}{Number of days between sequence upload and running pipeline}
}}

\examples{
data('aotus')
seqrec <- aotus@sqs@sqs[[1]]
# this is a SeqRec object
# it contains sequence records
show(seqrec)
# you can access its different data slots with @
seqrec@id       # sequence ID, accession + feature location
seqrec@nm       # feature name, '' if none
seqrec@accssn   # accession
seqrec@vrsn     # accession version
seqrec@url      # NCBI GenBank URL
seqrec@txid     # Taxonomic ID
seqrec@orgnsm   # free-text organism name
seqrec@sq       # sequence, in raw format
seqrec@dfln     # sequence definition
seqrec@ml_typ   # molecule type
seqrec@rec_typ  # whole record or feature
seqrec@nncltds  # sequence length
seqrec@nambgs   # number of non-ATCGs
seqrec@pambgs   # proportion of non-ATCGs
seqrec@gcr      # GC-ratio
seqrec@age      # days since being added to GenBank
# get the sequence like so....
(rawToChar(seqrec@sq))
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{TaxDict-class}}, \code{\link{TaxRec-class}},
  \code{\link{clusters2_run}}, \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline-tools.R
\name{stages_run}
\alias{stages_run}
\title{Sequentially run each stage}
\usage{
stages_run(wd, to, frm, stgs_msg, rstrt = FALSE)
}
\arguments{
\item{wd}{Working directory}

\item{to}{Total number of stages to run}

\item{frm}{Starting stage to run from}

\item{stgs_msg}{Printout stage message for log}

\item{rstrt}{Restarting, T/F}
}
\description{
Runs stages from \code{frm} to \code{to}. Records stage progress
in cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-calc.R
\name{calc_wrdfrq}
\alias{calc_wrdfrq}
\title{Calculate word frequencies}
\usage{
calc_wrdfrq(phylota, cid, min_frq = 0.1, min_nchar = 1,
  type = c("dfln", "nm"), ignr_pttrn = "[^a-z0-9]")
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s)}

\item{min_frq}{Minimum frequency}

\item{min_nchar}{Minimum number of characters for a word}

\item{type}{Definitions (dfln) or features (nm)}

\item{ignr_pttrn}{Ignore pattern, REGEX for text to ignore.}
}
\value{
list
}
\description{
For all sequences in a cluster(s) calculate the frequency of
separate words in either the sequence definitions or the reported feature
name.
}
\details{
By default, anything that is not alphanumeric is  ignored. 'dfln'
and 'nm' match the slot names in a SeqRec, see list_seqrec_slots().
}
\examples{
data('dragonflies')
# work out what gene region the cluster is likely representing with word freqs.
random_cids <- sample(dragonflies@cids, 10)
# most frequent words in definition line
(calc_wrdfrq(phylota = dragonflies, cid = random_cids, type = 'dfln'))
# most frequent words in feature name
(calc_wrdfrq(phylota = dragonflies, cid = random_cids, type = 'nm'))
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{drop_by_rank}}, \code{\link{drop_clstrs}},
  \code{\link{drop_sqs}}, \code{\link{get_clstr_slot}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_clstr_slot}
\alias{get_clstr_slot}
\title{Get slot data for each cluster record}
\usage{
get_clstr_slot(phylota, cid, slt_nm = list_clstrrec_slots())
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID}

\item{slt_nm}{Slot name}
}
\value{
vector
}
\description{
Get slot data for cluster(s)
}
\examples{
data('aotus')
random_cid <- sample(aotus@cids, 1)
(get_clstr_slot(phylota = aotus, cid = random_cid, slt_nm = 'seed'))
# see list_clstrrec_slots() for available slots
(list_clstrrec_slots())
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-blast.R
\name{blastn_run}
\alias{blastn_run}
\title{Launch blastn}
\usage{
blastn_run(dbfl, outfl, ps)
}
\arguments{
\item{dbfl}{Database file}

\item{outfl}{Output file}

\item{ps}{Parameters list, generated with parameters()}
}
\description{
Use \code{blastn} to BLAST all-vs-all using a BLAST
database.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{cache_rm}}, \code{\link{cache_setup}},
  \code{\link{clade_select}}, \code{\link{clstr2_calc}},
  \code{\link{clstr_all}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstr_direct}
\alias{clstr_direct}
\title{Cluster sequences directly associated with txid}
\usage{
clstr_direct(txid, sqs, txdct, ps, lvl)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{sqs}{Sequence object of all downloaded sequences}

\item{txdct}{Taxonomic dictionary}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
ClstrArc
}
\description{
In GenBank certain sequences may only be associated
with a higher level taxon (e.g. genus, family ...). This function
generates clusters from these sequences, alone. This function
identifies such sequences in the sequence object and generates
a list of clusters of cl_type 'direct'.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-blast.R
\name{blast_filter}
\alias{blast_filter}
\title{Filter BLAST results}
\usage{
blast_filter(blast_res, ps, lvl = 3)
}
\arguments{
\item{blast_res}{BLAST results}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
data.frame blast res
}
\description{
Given a BLAST output, filters query-subject pairs
such that only HSPs with a coverage greater than \code{mncvrg}
(specified in the pipeline parameters) remain. Filters both:
query-subject and subject-query pairs, if one of the coverages is
insufficient. HSP coverage is obtained from the BLAST column
\code{qcovs}.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_setup}},
  \code{\link{blast_sqs}}, \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2-tools.R
\name{seqrec_get}
\alias{seqrec_get}
\title{seqrec_get}
\usage{
seqrec_get(txid, ps, direct = FALSE, lvl = 0)
}
\arguments{
\item{txid}{NCBI taxonomic ID}

\item{ps}{Parameters list, generated with parameters()}

\item{direct}{Node-level only or subtree as well? Default FALSE.}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
Vector of sequence records
}
\description{
Downloads sequences from GenBank in batches.
}
\details{
If a restez database is available and the number of sequences to
retrieve is less than 'btchsz', the function will look the sequences up
from the database rather than download.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{sids_check}}, \code{\link{sids_get}},
  \code{\link{sids_load}}, \code{\link{sids_save}},
  \code{\link{sqs_count}}, \code{\link{sqs_save}},
  \code{\link{stage_args_check}}, \code{\link{stages_run}},
  \code{\link{tax_download}}, \code{\link{taxdict_gen}},
  \code{\link{taxtree_gen}}, \code{\link{txids_get}},
  \code{\link{txnds_count}}, \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{blast_clstr}
\alias{blast_clstr}
\title{Cluster BLAST Results}
\usage{
blast_clstr(blast_res)
}
\arguments{
\item{blast_res}{BLAST results}
}
\value{
List of list

list of cluster descriptions
}
\description{
Find single-linkage clusters from BLAST results.
Identifies seed sequence.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_filter}}, \code{\link{blast_setup}},
  \code{\link{blast_sqs}}, \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{ClstrArc-class}
\alias{ClstrArc-class}
\alias{ClstrArc-method}
\alias{as.character,ClstrArc-method}
\alias{show,ClstrArc-method}
\alias{print,ClstrArc-method}
\alias{str,ClstrArc-method}
\alias{summary,ClstrArc-method}
\alias{[[,ClstrArc,character-method}
\alias{[,ClstrArc,character,missing,missing-method}
\title{Cluster record archive}
\usage{
\S4method{as.character}{ClstrArc}(x)

\S4method{show}{ClstrArc}(object)

\S4method{print}{ClstrArc}(x)

\S4method{str}{ClstrArc}(object, max.level = 2L, ...)

\S4method{summary}{ClstrArc}(object)

\S4method{[[}{ClstrArc,character}(x, i)

\S4method{[}{ClstrArc,character,missing,missing}(x, i, j, ...,
  drop = TRUE)
}
\arguments{
\item{x}{\code{ClstrArc} object}

\item{object}{\code{ClstrArc} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}

\item{i}{cid(s)}

\item{j}{Unused}

\item{drop}{Unused}
}
\description{
Multiple cluster records.
}
\section{Slots}{

\describe{
\item{\code{ids}}{Vector of cluster record IDs}

\item{\code{clstrs}}{List of ClstrArc named by ID}
}}

\examples{
data('aotus')
clstrarc <- aotus@clstrs
# this is a ClstrArc object
# it contains cluster records
show(clstrarc)
# you can access its different data slots with @
clstrarc@ids     # unique cluster ID
clstrarc@clstrs  # list of cluster records
# access cluster records [[
(clstrarc[[clstrarc@ids[[1]]]])  # first cluster record
# generate new cluster archives with [
(clstrarc[clstrarc@ids[1:10]])  # first 10 clusters
}
\seealso{
Other run-public: \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotaR.R
\name{list_seqrec_slots}
\alias{list_seqrec_slots}
\title{List all SeqRec slots}
\usage{
list_seqrec_slots()
}
\value{
vector
}
\description{
Returns a vector of all available SeqRec slots of type
character, integer and numeric.
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{parameters_load}
\alias{parameters_load}
\title{Load parameters from cache}
\usage{
parameters_load(wd)
}
\arguments{
\item{wd}{Working directory}
}
\value{
Parameters list
}
\description{
Parameters are held in cache, use this function to
load parameters set for a wd.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_setup}},
  \code{\link{parent_get}}, \code{\link{progress_init}},
  \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{write_sqs}
\alias{write_sqs}
\title{Write out sequences}
\usage{
write_sqs(phylota, outfile, sid, sq_nm = sid, width = 80)
}
\arguments{
\item{phylota}{Phylota}

\item{outfile}{Output file}

\item{sid}{Sequence ID(s)}

\item{sq_nm}{Sequence name(s)}

\item{width}{Maximum number of characters in a line, integer}
}
\description{
Write out sequences, as .fasta, for a given vector of IDs.
}
\details{
The user can control the output definition lines of the sequences using the
sq_nm. By default sequences IDs are used. Note, ensure the sq_nm are in the
same order as sid.
}
\examples{
data('aotus')
# get sequences for a cluster and write out
random_cid <- sample(aotus@cids, 1)
sids <- aotus[[random_cid]]@sids
write_sqs(phylota = aotus, outfile = file.path(tempdir(), 'test.fasta'),
          sq_nm = 'my_gene', sid = sids)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-entrez.R
\name{searchterm_gen}
\alias{searchterm_gen}
\title{Construct GenBank Search Term}
\usage{
searchterm_gen(txid, ps, direct = FALSE)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{ps}{Parameters list, generated with parameters()}

\item{direct}{Node-level only or subtree as well? Default FALSE.}
}
\value{
character, search term
}
\description{
Construct search term for searching GenBank's
nucleotide database. Limits the maximum size of sequences, avoids
whole genome shotguns, predicted, unverified and synthetic
sequences.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{seeds_blast}}, \code{\link{seq_download}},
  \code{\link{seqarc_gen}}, \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-system.R
\name{cmdln}
\alias{cmdln}
\title{Run a command via terminal/command prompt}
\usage{
cmdln(cmd, args, ps, lgfl = NULL)
}
\arguments{
\item{cmd}{Command to be run}

\item{args}{Vector of command arguments, each parameter and value
must be a separate element}

\item{ps}{Paramters}

\item{lgfl}{File to which stdout/err will be written}
}
\value{
status, integer or character
}
\description{
Provide the command and arguments as a vector.
Also can take a lgfl to which all stdout and stderr is written.
If lgfl is not provided, a list is returned of 'status', 'stdout'
and 'stderr'. Else only the status is returned - 1 success, 0
failed.
}
\details{
Note, stdout/err are returned as 'raw'. Use rawToChar() to
convert to characters.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{is_txid_in_clstr}
\alias{is_txid_in_clstr}
\title{Is txid in cluster?}
\usage{
is_txid_in_clstr(phylota, txid, cid)
}
\arguments{
\item{phylota}{Phylota}

\item{txid}{Taxonomic ID}

\item{cid}{Cluster ID}
}
\value{
boolean
}
\description{
Checks if given txid is represented by any of the
sequences of a cluster by searching through all the sequence search
organism lineages.
}
\examples{
data(tinamous)
cid <- tinamous@cids[[1]]
clstr <- tinamous[[cid]]
sq <- tinamous[[clstr@sids[[1]]]]
txid <- sq@txid
# expect true
is_txid_in_clstr(phylota = tinamous, txid = txid, cid = cid)
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_sq_slot}},
  \code{\link{get_stage_times}}, \code{\link{get_tx_slot}},
  \code{\link{get_txids}}, \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{parameters_reset}
\alias{parameters_reset}
\title{Change parameters in a working directory}
\usage{
parameters_reset(wd, parameters, values)
}
\arguments{
\item{wd}{Working directory}

\item{parameters}{Parameters to be changed, vector.}

\item{values}{New values for each parameter, vector.}
}
\description{
Reset parameters after running \code{setup()}.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # run
  # run(wd = wd) # not running in test
  # use ctrl+c or Esc to kill
  # change parameters, e.g. min and max sequence lengths
  parameters_reset(wd = 'aotus', parameters = c('mnsql', 'mxsql'),
                   values = c(300, 1500))
  # see ?parameters
  # restart
  restart(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2.R
\name{seq_download}
\alias{seq_download}
\title{Download sequences for txids}
\usage{
seq_download(txids, txdct, ps)
}
\arguments{
\item{txids}{Taxonomic node IDs, numeric vector}

\item{txdct}{Taxonomic dictionary}

\item{ps}{Parameters list, generated with parameters()}
}
\description{
Look up and download all sequences for given
taxonomic IDs.
}
\details{
Sequence downloads are cached.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seqarc_gen}}, \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-special.R
\name{mk_txid_in_sq_mtrx}
\alias{mk_txid_in_sq_mtrx}
\title{Return matrix of txid in sequence}
\usage{
mk_txid_in_sq_mtrx(phylota, txids, sids = phylota@sids)
}
\arguments{
\item{phylota}{Phylota}

\item{txids}{Taxonomic IDs}

\item{sids}{Sequence IDs}
}
\value{
matrix
}
\description{
Searches through lineages of sequences' source organisms to
determine whether each txid is represented by the sequence.
}
\seealso{
Other tools-private: \code{\link{summary_phylota}},
  \code{\link{update_phylota}}
}
\concept{tools-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstr_all}
\alias{clstr_all}
\title{Hierarchically cluster all sequences of a txid}
\usage{
clstr_all(txid, sqs, txdct, ps, lvl = 0)
}
\arguments{
\item{txid}{Taxonomic ID}

\item{sqs}{Sequence object of all downloaded sequences}

\item{txdct}{Taxonomic dictionary}

\item{ps}{Parameters list, generated with parameters()}

\item{lvl}{Integer, number of message indentations indicating code
depth.}
}
\value{
ClstrArc
}
\description{
Identifies all direct and subtree clusters for a taxonomic ID.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_direct}},
  \code{\link{clstr_sqs}}, \code{\link{clstr_subtree}},
  \code{\link{clstrarc_gen}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_sq_slot}
\alias{get_sq_slot}
\title{Get slot data for each sequence}
\usage{
get_sq_slot(phylota, cid = NULL, sid = NULL,
  slt_nm = list_seqrec_slots())
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID}

\item{sid}{Sequence ID(s)}

\item{slt_nm}{Slot name}
}
\value{
vector
}
\description{
Get slot data for either or sequences in a cluster of
a vector of sequence IDs. Use list_seqrec_slots() for a list of
available slots.
}
\examples{
data('aotus')
random_sid <- sample(aotus@sids, 1)
(get_sq_slot(phylota = aotus, sid = random_sid, slt_nm = 'dfln'))
# see list_seqrec_slots() for available slots
(list_seqrec_slots())
}
\seealso{
Other tools-public: \code{\link{calc_mad}},
  \code{\link{calc_wrdfrq}}, \code{\link{drop_by_rank}},
  \code{\link{drop_clstrs}}, \code{\link{drop_sqs}},
  \code{\link{get_clstr_slot}}, \code{\link{get_nsqs}},
  \code{\link{get_ntaxa}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{setup}
\alias{setup}
\title{Set-up parameters}
\usage{
setup(wd, txid, ncbi_dr = ".", v = FALSE, overwrite = FALSE,
  outsider = FALSE, ...)
}
\arguments{
\item{wd}{Working directory}

\item{txid}{Root taxonomic ID(s), vector or numeric}

\item{ncbi_dr}{Directory to NCBI BLAST tools, default '.'}

\item{v}{Verbose, T/F}

\item{overwrite}{Overwrite existing cache?}

\item{outsider}{Run through \code{outsider}? T/F}

\item{...}{Additional parameters}
}
\description{
Set up working directory with parameters.
}
\details{
See \code{\link{parameters}}() for a description of all parameters
and their defaults. You can change parameters after a folder has been set up
with \code{\link{parameters_reset}}().
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  # e.g. "/usr/local/ncbi/blast/bin/"
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # see ?parameters for all available parameter options
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-sqs.R
\name{seqrec_convert}
\alias{seqrec_convert}
\title{Convert raw Entrez gb text record to SeqRecs}
\usage{
seqrec_convert(raw_recs, ps)
}
\arguments{
\item{raw_recs}{Raw text records returned from Entrez fetch}

\item{ps}{Parameters list, generated with parameters()}
}
\value{
SeqRecs
}
\description{
Parses returned sequences features with Entrez, returns one or
more SeqRec objects for each raw record.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-calc.R
\name{calc_mad}
\alias{calc_mad}
\title{Calculate MAD score}
\usage{
calc_mad(phylota, cid)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s)}
}
\value{
vector
}
\description{
For all sequences in a cluster(s) the MAD score.
}
\details{
MAD is a measure of the deviation in sequence length of a cluster.
Values range from 0 to 1. Clusters with values close to 1 have sequences with
similar lengths.
}
\examples{
data("bromeliads")
random_cids <- sample(bromeliads@cids, 10)
(calc_mad(phylota = bromeliads, cid = random_cids))
}
\seealso{
Other tools-public: \code{\link{calc_wrdfrq}},
  \code{\link{drop_by_rank}}, \code{\link{drop_clstrs}},
  \code{\link{drop_sqs}}, \code{\link{get_clstr_slot}},
  \code{\link{get_nsqs}}, \code{\link{get_ntaxa}},
  \code{\link{get_sq_slot}}, \code{\link{get_stage_times}},
  \code{\link{get_tx_slot}}, \code{\link{get_txids}},
  \code{\link{is_txid_in_clstr}},
  \code{\link{is_txid_in_sq}},
  \code{\link{list_clstrrec_slots}},
  \code{\link{list_ncbi_ranks}},
  \code{\link{list_seqrec_slots}},
  \code{\link{list_taxrec_slots}},
  \code{\link{plot_phylota_pa}},
  \code{\link{plot_phylota_treemap}},
  \code{\link{read_phylota}}, \code{\link{write_sqs}}
}
\concept{tools-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage3-tools.R
\name{clstrarc_gen}
\alias{clstrarc_gen}
\title{Generate cluster archive container class}
\usage{
clstrarc_gen(clstrrecs)
}
\arguments{
\item{clstrrecs}{list of ClstrRecs}
}
\value{
ClstrArc
}
\description{
Takes a list of ClstrRecs, returns a ClstrArc.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_join}},
  \code{\link{clstrrec_gen}}, \code{\link{clstrs_calc}},
  \code{\link{clstrs_join}}, \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage2.R
\name{download_run}
\alias{download_run}
\title{Run download stage}
\usage{
download_run(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Run the second stage of phylotaR, download. This stage
downloads sequences for all nodes with sequence numbers less than
mxsqs. It hierarchically traverses the taxonomy for each node and
downloads direct and subtree sequences for all descendants.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run stages
  taxise_run(wd = wd)
  download_run(wd = wd)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{TaxRec-class}
\alias{TaxRec-class}
\alias{TaxRec-method}
\alias{as.character,TaxRec-method}
\alias{show,TaxRec-method}
\alias{print,TaxRec-method}
\alias{str,TaxRec-method}
\alias{summary,TaxRec-method}
\title{Taxonomic record}
\usage{
\S4method{as.character}{TaxRec}(x)

\S4method{show}{TaxRec}(object)

\S4method{print}{TaxRec}(x)

\S4method{str}{TaxRec}(object, max.level = 2L, ...)

\S4method{summary}{TaxRec}(object)
}
\arguments{
\item{x}{\code{TaxRec} object}

\item{object}{\code{TaxRec} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}
}
\description{
Taxonomic dictionary contains a taxonomic
tree and NCBI taxonomy data for all taxonomic IDs.
}
\section{Slots}{

\describe{
\item{\code{id}}{Taxonomic ID}

\item{\code{scnm}}{Scientific name}

\item{\code{cmnm}}{Common name}

\item{\code{rnk}}{Rank}

\item{\code{lng}}{Lineage}

\item{\code{prnt}}{Parent}
}}

\examples{
data('aotus')
taxrec <- aotus@txdct@recs[[aotus@txdct@txids[[1]]]]
# this is a TaxRec object
# it contains NCBI's taxonomic information for a single node
show(taxrec)
# you can access its different data slots with @
taxrec@id    # taxonomic ID
taxrec@scnm  # scientific name
taxrec@cmnm  # common name, '' if none
taxrec@rnk   # rank
taxrec@lng   # lineage information: list of IDs and ranks
taxrec@prnt  # parent ID
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{clusters2_run}}, \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage1.R
\name{taxise_run}
\alias{taxise_run}
\title{Run taxise stage}
\usage{
taxise_run(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Run the first stage of phylotaR, taxise. This looks up
all descendant taxonomic nodes for a given taxonomic ID. It then
looks up relevant taxonomic information and generates a taxonomic
dictionary for user interaction after phylotaR has completed.
}
\details{
Objects will be cached.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run stages
  taxise_run(wd = wd)
}
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqArc-class}},
  \code{\link{SeqRec-class}}, \code{\link{TaxDict-class}},
  \code{\link{TaxRec-class}}, \code{\link{clusters2_run}},
  \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{SeqArc-class}
\alias{SeqArc-class}
\alias{SeqArc-method}
\alias{as.character,SeqArc-method}
\alias{show,SeqArc-method}
\alias{print,SeqArc-method}
\alias{str,SeqArc-method}
\alias{summary,SeqArc-method}
\alias{[[,SeqArc,character-method}
\alias{[,SeqArc,character,missing,missing-method}
\title{Sequence record archive}
\usage{
\S4method{as.character}{SeqArc}(x)

\S4method{show}{SeqArc}(object)

\S4method{print}{SeqArc}(x)

\S4method{str}{SeqArc}(object, max.level = 2L, ...)

\S4method{summary}{SeqArc}(object)

\S4method{[[}{SeqArc,character}(x, i)

\S4method{[}{SeqArc,character,missing,missing}(x, i, j, ..., drop = TRUE)
}
\arguments{
\item{x}{\code{SeqArc} object}

\item{object}{\code{SeqArc} object}

\item{max.level}{Maximum level of nesting for str()}

\item{...}{Further arguments for str()}

\item{i}{sid(s)}

\item{j}{Unused}

\item{drop}{Unused}
}
\description{
Multiple sequence records containing sequence data.
}
\details{
Sequences are stored as raw. Use rawToChar().
}
\section{Slots}{

\describe{
\item{\code{ids}}{Vector of Sequence Record IDs}

\item{\code{nncltds}}{Vector of sequence lengths}

\item{\code{nambgs}}{Vector of number of ambiguous nucleotides}

\item{\code{txids}}{Vector source txid associated with each sequence}

\item{\code{sqs}}{List of SeqRecs named by ID}
}}

\examples{
data('aotus')
seqarc <- aotus@sqs
# this is a SeqArc object
# it contains sequence records
show(seqarc)
# you can access its different data slots with @
seqarc@ids     # sequence IDs defined as accession + feature position
seqarc@nncltds # number of nucleotides of all sequences  
seqarc@nambgs  # number of ambiguous nucleotides of all sequences
seqarc@txids   # all the taxonomic IDs for all sequences
seqarc@sqs     # list of all SeqRecs
# access sequence records [[
(seqarc[[seqarc@ids[[1]]]])  # first sequence record
# generate new sequence archives with [
(seqarc[seqarc@ids[1:10]])  # first 10 sequences
}
\seealso{
Other run-public: \code{\link{ClstrArc-class}},
  \code{\link{ClstrRec-class}},
  \code{\link{Phylota-class}}, \code{\link{SeqRec-class}},
  \code{\link{TaxDict-class}}, \code{\link{TaxRec-class}},
  \code{\link{clusters2_run}}, \code{\link{clusters_run}},
  \code{\link{parameters_reset}}, \code{\link{reset}},
  \code{\link{restart}}, \code{\link{run}},
  \code{\link{setup}}, \code{\link{taxise_run}}
}
\concept{run-public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{progress_save}
\alias{progress_save}
\title{Save current progress}
\usage{
progress_save(wd, stg)
}
\arguments{
\item{wd}{Working directory}

\item{stg}{Stage}
}
\description{
Stores the pipeline progress in the cache.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-log.R
\name{info}
\alias{info}
\title{Write info message to log}
\usage{
info(lvl, ps, ...)
}
\arguments{
\item{lvl}{Integer, number of message indentations indicating code
depth.}

\item{ps}{Parameters list, generated with parameters()}

\item{...}{Message elements for concatenating}
}
\description{
Inform a user via log.txt of pipeline progress.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{obj_save}},
  \code{\link{outfmt_get}}, \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools-cache.R
\name{obj_save}
\alias{obj_save}
\title{Save a named object in the cache}
\usage{
obj_save(wd, obj, nm)
}
\arguments{
\item{wd}{Working directory}

\item{obj}{Object}

\item{nm}{Object name}
}
\description{
Save an object in the cache that can be loaded by
\code{obj_load}.
}
\seealso{
Other run-private: \code{\link{batcher}},
  \code{\link{blast_clstr}}, \code{\link{blast_filter}},
  \code{\link{blast_setup}}, \code{\link{blast_sqs}},
  \code{\link{blastcache_load}},
  \code{\link{blastcache_save}}, \code{\link{blastdb_gen}},
  \code{\link{blastn_run}}, \code{\link{cache_rm}},
  \code{\link{cache_setup}}, \code{\link{clade_select}},
  \code{\link{clstr2_calc}}, \code{\link{clstr_all}},
  \code{\link{clstr_direct}}, \code{\link{clstr_sqs}},
  \code{\link{clstr_subtree}}, \code{\link{clstrarc_gen}},
  \code{\link{clstrarc_join}}, \code{\link{clstrrec_gen}},
  \code{\link{clstrs_calc}}, \code{\link{clstrs_join}},
  \code{\link{clstrs_merge}},
  \code{\link{clstrs_renumber}}, \code{\link{clstrs_save}},
  \code{\link{cmdln}}, \code{\link{descendants_get}},
  \code{\link{download_obj_check}}, \code{\link{error}},
  \code{\link{gb_extract}},
  \code{\link{hierarchic_download}}, \code{\link{info}},
  \code{\link{ncbicache_load}},
  \code{\link{ncbicache_save}}, \code{\link{obj_check}},
  \code{\link{obj_load}}, \code{\link{outfmt_get}},
  \code{\link{parameters_load}},
  \code{\link{parameters_setup}}, \code{\link{parent_get}},
  \code{\link{progress_init}}, \code{\link{progress_read}},
  \code{\link{progress_reset}},
  \code{\link{progress_save}}, \code{\link{rank_get}},
  \code{\link{rawseqrec_breakdown}},
  \code{\link{safely_connect}},
  \code{\link{search_and_cache}},
  \code{\link{searchterm_gen}}, \code{\link{seeds_blast}},
  \code{\link{seq_download}}, \code{\link{seqarc_gen}},
  \code{\link{seqrec_augment}},
  \code{\link{seqrec_convert}}, \code{\link{seqrec_gen}},
  \code{\link{seqrec_get}}, \code{\link{sids_check}},
  \code{\link{sids_get}}, \code{\link{sids_load}},
  \code{\link{sids_save}}, \code{\link{sqs_count}},
  \code{\link{sqs_save}}, \code{\link{stage_args_check}},
  \code{\link{stages_run}}, \code{\link{tax_download}},
  \code{\link{taxdict_gen}}, \code{\link{taxtree_gen}},
  \code{\link{txids_get}}, \code{\link{txnds_count}},
  \code{\link{warn}}
}
\concept{run-private}
