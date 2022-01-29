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

