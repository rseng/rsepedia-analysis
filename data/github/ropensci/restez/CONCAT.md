
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- devtools::rmarkdown::render("README.Rmd") -->

<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->

# Locally query GenBank <img src="https://raw.githubusercontent.com/ropensci/restez/master/logo.png" height="200" align="right"/>

[![Build
Status](https://travis-ci.org/ropensci/restez.svg?branch=master)](https://travis-ci.org/ropensci/restez)
[![Coverage
Status](https://coveralls.io/repos/github/ropensci/restez/badge.svg?branch=master)](https://coveralls.io/github/ropensci/restez?branch=master)
[![ROpenSci
status](https://badges.ropensci.org/232_status.svg)](https://github.com/ropensci/onboarding/issues/232)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/restez)](https://CRAN.R-project.org/package=restez)
[![DOI](https://zenodo.org/badge/129107980.svg)](https://zenodo.org/badge/latestdoi/129107980)
[![status](http://joss.theoj.org/papers/6eb3ba7dddbdab8788a430eb62fc3841/status.svg)](http://joss.theoj.org/papers/6eb3ba7dddbdab8788a430eb62fc3841)

> NOTE: `restez` is no longer available on CRAN due to the archiving of
> a key dependency. It can still be installed via GitHub. The issue is
> being dealt with and hopefully a new version of `restez` will be
> available on CRAN soon.

Download parts of [NCBI’s GenBank](https://www.ncbi.nlm.nih.gov/nuccore)
to a local folder and create a simple SQL-like database. Use ‘get’ tools
to query the database by accession IDs.
[rentrez](https://github.com/ropensci/rentrez) wrappers are available,
so that if sequences are not available locally they can be searched for
online through [Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25500/).

See the [detailed
tutorials](https://ropensci.github.io/restez/articles/restez.html) for
more information.

## Introduction

*Vous entrez, vous rentrez et, maintenant, vous …. restez\!*

Downloading sequences and sequence information from GenBank and related
NCBI taxonomic databases is often performed via the NCBI API, Entrez.
Entrez, however, has a limit on the number of requests and downloading
large amounts of sequence data in this way can be inefficient. For
programmatic situations where multiple Entrez calls are made,
downloading may take days, weeks or even months.

This package aims to make sequence retrieval more efficient by allowing
a user to download large sections of the GenBank database to their local
machine and query this local database either through package specific
functions or Entrez wrappers. This process is more efficient as GenBank
downloads are made via NCBI’s FTP using compressed sequence files. With
a good internet connection and a middle-of-the-road computer, a database
comprising 20 GB of sequence information can be generated in less than
10
minutes.

<img src="https://raw.githubusercontent.com/ropensci/restez/master/paper/outline.png" height="500" align="center"/>

## Installation

<!--
`restez` is available via CRAN and can be installed:


```r
install.packages("restez")
```
-->

The package can currently only be installed through GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/restez")
```

(It was previously available via CRAN but was archived due to a key
dependency [MonetDBLite](https://github.com/MonetDB/MonetDBLite-R) being
no longer available.)

## Quick Examples

> For more detailed information on the package’s functions and detailed
> guides on downloading, constructing and querying a database, see the
> [detailed
> tutorials](https://ropensci.github.io/restez/articles/restez.html).

### Setup

``` r
# Warning: running these examples may take a few minutes
library(restez)
#> -------------
#> restez v1.0.2
#> -------------
#> Remember to restez_path_set() and, then, restez_connect()
# choose a location to store GenBank files
restez_path_set(rstz_pth)
```

``` r
# Run the download function
db_download()
# after download, create the local database
db_create()
```

### Query

``` r
# connect, ensure safe disconnect after finishing
restez_connect()
#> Remember to run `restez_disconnect()`
# get a random accession ID from the database
id <- sample(list_db_ids(), 1)
#> Warning in list_db_ids(): Number of ids returned was limited to [100].
#> Set `n=NULL` to return all ids.
# you can extract:
# sequences
seq <- gb_sequence_get(id)[[1]]
str(seq)
#>  chr "ACTCTGACTTTTTACTGTATATAAAAACAGCTTTTTGGTTTATACTTGAATTCAGGAATAACCAAGCAGGTGTAAATATGCCAGCGCAAGAACAGCAAATTT"
# definitions
def <- gb_definition_get(id)[[1]]
print(def)
#> [1] "Unidentified RNA clone P10.7"
# organisms
org <- gb_organism_get(id)[[1]]
print(org)
#> [1] "unidentified"
# or whole records
rec <- gb_record_get(id)[[1]]
cat(rec)
#> LOCUS       AF040899                 102 bp    RNA     linear   UNA 06-MAR-1998
#> DEFINITION  Unidentified RNA clone P10.7.
#> ACCESSION   AF040899
#> VERSION     AF040899.1
#> KEYWORDS    .
#> SOURCE      unidentified
#>   ORGANISM  unidentified
#>             unclassified sequences.
#> REFERENCE   1  (bases 1 to 102)
#>   AUTHORS   Pan,W.S., Ji,X.Y., Wang,H.T. and Zhong,Y.S.
#>   TITLE     RNA from plasma of patient NO.10
#>   JOURNAL   Unpublished
#> REFERENCE   2  (bases 1 to 102)
#>   AUTHORS   Pan,W.S., Ji,X.Y., Wang,H.T. and Zhong,Y.S.
#>   TITLE     Direct Submission
#>   JOURNAL   Submitted (31-DEC-1997) Department of Applied Molecular Biology,
#>             Microbiology & Epidemiology Institution, 20 Dongdajie Street,
#>             Fengtai, Beijing 100071, China
#> FEATURES             Location/Qualifiers
#>      source          1..102
#>                      /organism="unidentified"
#>                      /mol_type="genomic RNA"
#>                      /db_xref="taxon:32644"
#>                      /clone="P10.7"
#>                      /note="from the plasma of patient no.10, a person infected
#>                      by an unknown hepatitis virus"
#> ORIGIN      
#>         1 actctgactt tttactgtat ataaaaacag ctttttggtt tatacttgaa ttcaggaata
#>        61 accaagcagg tgtaaatatg ccagcgcaag aacagcaaat tt
#> //
```

### Entrez wrappers

``` r
# use the entrez_* wrappers to access GB data
res <- entrez_fetch(db = 'nucleotide', id = id, rettype = 'fasta')
cat(res)
#> >AF040899.1 Unidentified RNA clone P10.7
#> ACTCTGACTTTTTACTGTATATAAAAACAGCTTTTTGGTTTATACTTGAATTCAGGAATAACCAAGCAGG
#> TGTAAATATGCCAGCGCAAGAACAGCAAATTT
# if the id is not in the local database
# these wrappers will search online via the rentrez package
res <- entrez_fetch(db = 'nucleotide', id = c('S71333.1', id),
                    rettype = 'fasta')
#> [1] id(s) are unavailable locally, searching online.
cat(res)
#> >AF040899.1 Unidentified RNA clone P10.7
#> ACTCTGACTTTTTACTGTATATAAAAACAGCTTTTTGGTTTATACTTGAATTCAGGAATAACCAAGCAGG
#> TGTAAATATGCCAGCGCAAGAACAGCAAATTT
#> 
#> >S71333.1 alpha 1,3 galactosyltransferase [New World monkeys, mermoset lymphoid cell line B95.8, mRNA Partial, 1131 nt]
#> ATGAATGTCAAAGGAAAAGTAATTCTGTCGATGCTGGTTGTCTCAACTGTGATTGTTGTGTTTTGGGAAT
#> ATATCAACAGCCCAGAAGGCTCTTTCTTGTGGATATATCACTCAAAGAACCCAGAAGTTGATGACAGCAG
#> TGCTCAGAAGGACTGGTGGTTTCCTGGCTGGTTTAACAATGGGATCCACAATTATCAACAAGAGGAAGAA
#> GACACAGACAAAGAAAAAGGAAGAGAGGAGGAACAAAAAAAGGAAGATGACACAACAGAGCTTCGGCTAT
#> GGGACTGGTTTAATCCAAAGAAACGCCCAGAGGTTATGACAGTGACCCAATGGAAGGCGCCGGTTGTGTG
#> GGAAGGCACTTACAACAAAGCCATCCTAGAAAATTATTATGCCAAACAGAAAATTACCGTGGGGTTGACG
#> GTTTTTGCTATTGGAAGATATATTGAGCATTACTTGGAGGAGTTCGTAACATCTGCTAATAGGTACTTCA
#> TGGTCGGCCACAAAGTCATATTTTATGTCATGGTGGATGATGTCTCCAAGGCGCCGTTTATAGAGCTGGG
#> TCCTCTGCGTTCCTTCAAAGTGTTTGAGGTCAAGCCAGAGAAGAGGTGGCAAGACATCAGCATGATGCGT
#> ATGAAGACCATCGGGGAGCACATCTTGGCCCACATCCAACACGAGGTTGACTTCCTCTTCTGCATGGATG
#> TGGACCAGGTCTTCCAAGACCATTTTGGGGTAGAGACCCTGGGCCAGTCGGTGGCTCAGCTACAGGCCTG
#> GTGGTACAAGGCAGATCCTGATGACTTTACCTATGAGAGGCGGAAAGAGTCGGCAGCATATATTCCATTT
#> GGCCAGGGGGATTTTTATTACCATGCAGCCATTTTTGGAGGAACACCGATTCAGGTTCTCAACATCACCC
#> AGGAGTGCTTTAAGGGAATCCTCCTGGACAAGAAAAATGACATAGAAGCCGAGTGGCATGATGAAAGCCA
#> CCTAAACAAGTATTTCCTTCTCAACAAACCCTCTAAAATCTTATCTCCAGAATACTGCTGGGATTATCAT
#> ATAGGCCTGCCTTCAGATATTAAAACTGTCAAGCTATCATGGCAAACAAAAGAGTATAATTTGGTTAGAA
#> AGAATGTCTGA
restez_disconnect()
```

## Contributing

Want to contribute? Check the [contributing
page](https://ropensci.github.io/restez/CONTRIBUTING.html).

## Version

Release version 1.

## Licence

MIT

## Citation

Bennett et al. (2018). restez: Create and Query a Local Copy of GenBank
in R. *Journal of Open Source Software*, 3(31), 1102.
<https://doi.org/10.21105/joss.01102>

## References

Benson, D. A., Karsch-Mizrachi, I., Clark, K., Lipman, D. J., Ostell,
J., & Sayers, E. W. (2012). GenBank. *Nucleic Acids Research*,
40(Database issue), D48–D53. <http://doi.org/10.1093/nar/gkr1202>

Winter DJ. (2017) rentrez: An R package for the NCBI eUtils API. *PeerJ
Preprints* 5:e3179v2 <https://doi.org/10.7287/peerj.preprints.3179v2>

## Maintainer

[Dom
Bennett](https://github.com/DomBennett)

-----

[![ropensci\_footer](http://ropensci.org/public_images/ropensci_footer.png)](http://ropensci.org)
# restez (in development)

## Bug fixes

* Check internet connection in mainland China

# restez 1.0.0

## Post-review version of `restez` released.

* Download and query parts of GenBank from within R
# Contributing

You are very welcome to help out in the development of restez. The NCBI databases and resources are vast, it is not possible
for a single person from a single discpline within the biological sciences to effectively realise all of what NCBI has to offer.
Restez needs your help!

If you have any ideas for future features than please add them to the [issues page](https://github.com/AntonelliLab/restez/issues).
If you have the guile, time and inspriation to add those features yourself, then please fork and send a pull request.

## Areas for possible contribution

### Protein database

Currently restez only downloads the nucleotide database (synonyms: GenBank, nuccore). All of the restez functions for downloading,
creating and querying the nucleotide database could easily be copied and reengineered for working with a protein database.

### Taxonomy

Likewise the taxonomic database could also be downloaded and integrated into the restez framework. The taxonomic database is,
however, stored different from the nucleotide given its unique nature. Better would be to make use of currently existing packages
(e.g taxizedb) and integrate them into restez. This may allow users to have simple entrez_search functions when looking up
sequence ids (e.g. retreiving all sequences associated with a particular taxonomic group.)

### Retmodes

Restez tries to recreate the output from rentrez wrappers as best it can. For any queries with retmodes that are not text-based
(e.g. xml, feature tables), however, restez must search online. This is because all NCBI FTP downloads are text based 'flatfiles'.
Currently, the only way to recreate non-text based data would be to convert the download flatfiles to other formats upon request.
This is far from an ideal scenario as there would be no guarrantee that a restez result would match a rentrez result for the same
query.

## How to contribute

To contribute you will need a GitHub account and to have basic knowledge of the R language. You can then create a fork of the
repo in your own GitHub account and download the repository to your local machine. `devtools` is recommended.

```r
devtools::install_github('[your account]/restez')
```

All new functions must be tested. For every new file in `R/`, a new test file must be created in `tests/testthat/`. To test the
package and make sure it meets CRAN guidelines use `devtools`. 

```r
devtools::test()
devtools::check_cran()
```

For help, refer to Hadley Wickham's book, [R packages](http://r-pkgs.had.co.nz/).

## Style guide

Restez is being developed for eventual submission to ROpenSci. This means the package and its code should meet ROpenSci style and
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
---
title: 'restez: Create and Query a Local Copy of GenBank in R'
tags:
  - GenBank
  - nucleotides
  - R
  - DNA
  - sequence
  - NCBI
  - rstats
authors:
 - name: Dominic J. Bennett
   orcid: 0000-0003-2722-1359
   affiliation: "1, 2"
 - name: Hannes Hettling
   affiliation: "3"
 - name: Daniele Silvestro
   orcid: 0000-0003-0100-0961
   affiliation: "1, 2"
 - name: Rutger Vos
   orcid: 0000-0001-9254-7318
   affiliation: "3"
 - name: Alexandre Antonelli
   orcid: 0000-0003-1842-9297
   affiliation: "1, 2, 4"
affiliations:
 - name: Gothenburg Global Biodiversity Centre, Box 461, SE-405 30 Gothenburg, Sweden
   index: 1
 - name: Department of Biological and Environmental Sciences, University of Gothenburg, Box 461, SE-405 30 Gothenburg, Sweden
   index: 2
 - name: Naturalis Biodiversity Center, P.O. Box 9517, 2300 RA Leiden, The Netherlands
   index: 3
 - name: Gothenburg Botanical Garden, SE 41319 Gothenburg, Sweden
   index: 4
date: 27 November 2018
bibliography: paper.bib
---

# Summary

Downloading sequences and sequence information from GenBank [@Benson2013] and related NCBI databases is often performed via the NCBI API, Entrez [@Ostell2002]. Entrez, however, has a limit on the number of requests, thus downloading large amounts of sequence data in this way can be inefficient. For situations where a large number of Entrez calls is made, downloading may take days, weeks or even months and could result in a user’s IP address being blacklisted from the NCBI services due to server overload. Additionally, Entrez limits the number of entries that can be retrieved at once, requiring a user to develop code for querying in batches.

The `restez` package [@restez_z] aims to make sequence retrieval more efficient by allowing a user to download the GenBank database, either in its entirety or in subsets, to their local machine and query this local database instead. This process is more time efficient as GenBank downloads are made via NCBI’s FTP server using compressed sequence files. With a good internet connection and a computer with currently standard capabilities, a database comprising 7 GB of sequence information (i.e. the total sequence data available for Rodentia as of 27 June 2018) can be generated in less than 10 minutes. (For an outline of the functions and structure of `restez`, see Figure 1.)

![The functions and file structure for downloading, setting up and querying a local copy of GenBank](https://raw.githubusercontent.com/ropensci/restez/master/paper/outline.png)

## Rentrez integration

`rentrez` [@Winter2017] is a popular R package for querying NCBI’s databases via Entrez in R. To maximize the compatibility of `restez`, we implemented wrapper functions with the same names and arguments as the `rentrez` equivalents. Whenever a wrapper function is called the local database copy is searched first. If IDs are missing in the local database a secondary call to Entrez is made via the internet. This allows for easy employment of `restez` in scripts and packages that are already using `rentrez`. At a minimum, a user currently using `rentrez` will only need to create a local subset of the GenBank database, call `restez` instead of `rentrez` and ensure the `restez` database is connected.

------

# Examples

## A small example

After a restez database has been set-up, we can retrieve all the sequences from an `rentrez::entrez_search()` with a single command.

```r
# Use rentrez to search for accession IDs of interest
# Sequences in fasta format can then be retrieved with entrez_fetch
res <- rentrez::entrez_fetch(db = 'nuccore', id = ids, rettype = 'fasta')
# ^ likely to raise an error if too many IDs
res <- restez::entrez_fetch(db = 'nuccore', id = ids, rettype = 'fasta')
# ^ not likely to raise an error
```

## A large example

`phylotaR` is an R package for  retrieving and identifying orthologous sequence clusters from GenBank as a first step in a phylogenetic analysis [@Bennett2018]. Because the package runs an automated pipeline, multiple queries to GenBank via Entrez are made using the `rentrez` package. As a result, for large taxonomic groups containing well-sequenced organisms the pipeline can take a long time to complete.

```r
library(phylotaR)
# run phylotaR pipeline for New World Monkeys
txid <- 9479  # taxonomic ID
setup(wd = 'nw_monkeys', txid = txid)
run(wd = wd)
# ^ takes around 40 minutes
```

We can download and create a local copy of the primates GenBank locally and re-run the above code with a library call to `restez` for speed-up gains and increased code reliability.

```r
# setup database
library(restez)
# Specify path to a local directory in which database will be stored
# Make sure you have sufficient disk space!
restez_path_set(filepath = 'restez_db')
db_download(db = 'nucleotide') # Interactively download GenBank data
db_create(db = 'nucleotide')
```
Now when re-running the first `phylotaR` code block with the inclusion of the `restez` package, the procedure completes approximately eight times faster.

```r
# run phylotaR again
library(phylotaR)
library(restez)
restez_path_set(filepath = 'restez_db')
txid <- 9479
setup(wd = 'nw_monkeys', txid = txid)
run(wd = wd)
# ^ takes around 5 minutes
```

**For more detailed and up-to-date examples and tutorials, see the `restez` GitHub page [@restez_gh].**

# Availability

`restez` is open source software made available under the MIT license. It can be installed through CRAN [@restez_cran], `install.package("restez")`, or from its GitHub source code repository using the `devtools` package, e.g. as follows: `devtools::install_github("ropensci/restez")`

# Funding

This package has been developed as part of the supersmartR project [@supersmartR] which has received funding through A.A. (from the Swedish Research Council [B0569601], the Swedish Foundation for Strategic Research, a Wallenberg Academy Fellowship, the Faculty of Sciences at the University of Gothenburg, the Wenner-Gren Foundations, and the David Rockefeller Center for Latin American Studies at Harvard University) and through D.S. (from the Swedish Research Council [2015-04748]).

# References
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

Please note that the restez project is released with a
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
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- devtools::rmarkdown::render("README.Rmd") -->
<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
# don't run callr
assignInNamespace(x = 'custom_download2', value = restez:::custom_download,
                  ns = 'restez')
assignInNamespace(x = 'gb_build2', value = restez:::gb_build,
                  ns = 'restez')
```

# Locally query GenBank <img src="https://raw.githubusercontent.com/ropensci/restez/master/logo.png" height="200" align="right"/>

[![Build Status](https://travis-ci.org/ropensci/restez.svg?branch=master)](https://travis-ci.org/ropensci/restez) [![Coverage Status](https://coveralls.io/repos/github/ropensci/restez/badge.svg?branch=master)](https://coveralls.io/github/ropensci/restez?branch=master) [![ROpenSci status](https://badges.ropensci.org/232_status.svg)](https://github.com/ropensci/onboarding/issues/232) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/restez)](https://CRAN.R-project.org/package=restez) [![DOI](https://zenodo.org/badge/129107980.svg)](https://zenodo.org/badge/latestdoi/129107980) [![status](http://joss.theoj.org/papers/6eb3ba7dddbdab8788a430eb62fc3841/status.svg)](http://joss.theoj.org/papers/6eb3ba7dddbdab8788a430eb62fc3841)

> NOTE: `restez` is no longer available on CRAN due to the archiving of a key dependency. It can still be installed via GitHub. The issue is being dealt with and hopefully a new version of `restez` will be available on CRAN soon.

Download parts of [NCBI's GenBank](https://www.ncbi.nlm.nih.gov/nuccore) to a local folder and create a simple SQL-like database. Use 'get' tools to query the database by accession IDs. [rentrez](https://github.com/ropensci/rentrez) wrappers are available, so that if sequences are not available locally they can be searched for online through [Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25500/).

See the [detailed  tutorials](https://ropensci.github.io/restez/articles/restez.html) for more information.

## Introduction

*Vous entrez, vous rentrez et, maintenant, vous .... restez!*

Downloading sequences and sequence information from GenBank and related NCBI taxonomic databases is often performed via the NCBI API, Entrez. Entrez, however, has a limit on the number of requests and downloading large amounts of sequence data in this way can be inefficient. For programmatic situations where multiple Entrez calls are made, downloading may take days, weeks or even months.

This package aims to make sequence retrieval more efficient by allowing a user to download large sections of the GenBank database to their local machine and query this local database either through package specific functions or Entrez wrappers. This process is more efficient as GenBank downloads are made via NCBI's FTP using compressed sequence files. With a good internet connection and a middle-of-the-road computer, a database comprising 20 GB of sequence information can be generated in less than 10 minutes.

<img src="https://raw.githubusercontent.com/ropensci/restez/master/paper/outline.png" height="500" align="center"/>

## Installation

<!--
`restez` is available via CRAN and can be installed:

```{r cran-installation, include=TRUE, echo=TRUE, eval=FALSE}
install.packages("restez")
```
-->

The package can currently only be installed through GitHub:

```{r gh-installation, include=TRUE, echo=TRUE, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/restez")
```

(It was previously available via CRAN but was archived due to a key dependency [MonetDBLite](https://github.com/MonetDB/MonetDBLite-R) being no longer available.)

## Quick Examples

> For more detailed information on the package's functions and detailed guides on downloading, constructing and querying a database, see the [detailed  tutorials](https://ropensci.github.io/restez/articles/restez.html).

### Setup

```{r presetup, include=FALSE}
rstz_pth <- tempdir()
restez::restez_path_set(filepath = rstz_pth)
restez::db_delete(everything = TRUE)
```
```{r restez-setup, echo=TRUE, eval=TRUE, results='hide'}
# Warning: running these examples may take a few minutes
library(restez)
# choose a location to store GenBank files
restez_path_set(rstz_pth)
```
``` {r reassigndownload, include=FALSE, eval=TRUE}
db_download <- function() {
  restez::db_download(preselection = '20')
}
```
```{r gb-download, echo=TRUE, eval=TRUE, results='hide'}
# Run the download function
db_download()
# after download, create the local database
db_create()
```

### Query

```{r query, echo=TRUE, eval=TRUE}
# connect, ensure safe disconnect after finishing
restez_connect()
# get a random accession ID from the database
id <- sample(list_db_ids(), 1)
# you can extract:
# sequences
seq <- gb_sequence_get(id)[[1]]
str(seq)
# definitions
def <- gb_definition_get(id)[[1]]
print(def)
# organisms
org <- gb_organism_get(id)[[1]]
print(org)
# or whole records
rec <- gb_record_get(id)[[1]]
cat(rec)
```

### Entrez wrappers

```{r entrez, echo=TRUE, eval=TRUE}
# use the entrez_* wrappers to access GB data
res <- entrez_fetch(db = 'nucleotide', id = id, rettype = 'fasta')
cat(res)
# if the id is not in the local database
# these wrappers will search online via the rentrez package
res <- entrez_fetch(db = 'nucleotide', id = c('S71333.1', id),
                    rettype = 'fasta')
cat(res)
restez_disconnect()
```

## Contributing

Want to contribute? Check the [contributing page](https://ropensci.github.io/restez/CONTRIBUTING.html).

## Version

Release version 1.

## Licence

MIT

## Citation

Bennett et al. (2018).  restez:  Create and Query a Local Copy of GenBank in R.
*Journal of Open Source Software*, 3(31), 1102. https://doi.org/10.21105/joss.01102

## References

Benson, D. A., Karsch-Mizrachi, I., Clark, K., Lipman, D. J., Ostell, J., &
Sayers, E. W. (2012). GenBank. *Nucleic Acids Research*, 40(Database issue),
D48–D53. http://doi.org/10.1093/nar/gkr1202

Winter DJ. (2017) rentrez: An R package for the NCBI eUtils API.
*PeerJ Preprints* 5:e3179v2 https://doi.org/10.7287/peerj.preprints.3179v2

## Maintainer

[Dom Bennett](https://github.com/DomBennett)

-----

[![ropensci_footer](http://ropensci.org/public_images/ropensci_footer.png)](http://ropensci.org)
---
title: "3. Advanced parsing of a GenBank record"
output: rmarkdown::html_vignette
---



In this tutorial we will demonstrate how the `gb_extract()` function works. `restez` downloads and stores all GenBank records in text format. Ordinarily, to be able to extract relevant bits of information from a text record in a systematic way we would need to make use of [regular expressions](https://en.wikipedia.org/wiki/Regular_expression). `gb_extract()` uses regular expressions, so we don't have to. Here's how it works.

## A GenBank record

restez comes with an example GenBank record, AY952423, which can be viewed [online](https://www.ncbi.nlm.nih.gov/nuccore/AY952423.1). We can retrieve the record with `rentrez` or load it using the `data()`. We can visualise the record in R as it appears online using `cat()` -- like print but with newline spaces parsed correctly.


```r
library(restez)
# record <- rentrez::entrez_fetch(db = 'nucleotide', id = 'AY952423', rettype = 'gb', retmode = 'text')
data(record)
cat(record)
#> LOCUS       AY952423                2623 bp    DNA     linear   PLN 17-APR-2005
#> DEFINITION  Livistona chinensis tRNA-Lys (trnK) gene, partial sequence; and
#>             matK gene, complete sequence; chloroplast.
#> ACCESSION   AY952423
#> VERSION     AY952423.1
#> KEYWORDS    .
#> SOURCE      chloroplast Livistona chinensis
#>   ORGANISM  Livistona chinensis
#>             Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
#>             Spermatophyta; Magnoliophyta; Liliopsida; Arecaceae; Coryphoideae;
#>             Livistoneae; Livistoninae; Livistona.
#> REFERENCE   1  (bases 1 to 2623)
#>   AUTHORS   Li,X.X. and Zhou,Z.K.
#>   TITLE     Monocotyledons phylogeny based on three genes (matK, rbcL and 18S
#>             rDNA) sequences
#>   JOURNAL   Unpublished
#> REFERENCE   2  (bases 1 to 2623)
#>   AUTHORS   Li,X.X. and Zhou,Z.K.
#>   TITLE     Direct Submission
#>   JOURNAL   Submitted (03-MAR-2005) Taxonomical and Ethnobotanical Department,
#>             Kunming Institute of Botany, The Chinese Academy of Sciences,
#>             Heilongtan, Kunming, Yunnan 650204, China
#> FEATURES             Location/Qualifiers
#>      source          1..2623
#>                      /organism="Livistona chinensis"
#>                      /organelle="plastid:chloroplast"
#>                      /mol_type="genomic DNA"
#>                      /db_xref="taxon:115492"
#>      gene            <1..>2623
#>                      /gene="trnK"
#>                      /note="tRNA-Lys"
#>      intron          <1..>2623
#>                      /gene="trnK"
#>      gene            813..2347
#>                      /gene="matK"
#>      misc_feature    813..2347
#>                      /gene="matK"
#>                      /note="similar to maturase K"
#> ORIGIN      
#>         1 attggggttg ctaactcaac ggtagagtac tcggctttta agtgcgacta tcatctttta
#>        61 cacatttgga tgaagtaagg aattcgtcca gactattggt agagtctata agaccacgac
#>       121 tgatcctgaa aggtaatgaa tggaaaaaat agcatgtcgt acgtaataca atgagaaact
#>       181 tgtaatttct tattgtaatt ttttaagtag aactttgagt ttatccttac tggatcatta
#>       241 caaaaatatt gtattttatt tttggaaggg gacgaaaaaa aggaaattcc caacatttat
#>       301 tgtttggtct aatgaataaa tggatagggg cctagggtag ggcccaattt ttgtaaaaca
#>       361 aaaagcaacg agcttatgtt cttaatttga ataattaccc gatctaatta gatgttaaaa
#>       421 ataaattagt gccagatgtg gtaaagggtt ctactgtaag tggacctttt tttttttttt
#>       481 ttatgaatcc tacctattat ctattatgga ttaaagatgg atgtgtataa gaagaagtat
#>       541 actgataaag agaatttttc caaagtcaaa agagcaatcg ggttgcaaaa ataaaggatt
#>       601 tttacctccg agaattataa attaattgga tcaaaaggag aggaaaaagt ctgtgattgg
#>       661 actccttcta tccgcgggta tgggtatata gtaggtatat atgtatattt gtatactata
#>       721 taaattacat gccctgttct gaccgtattg cactatgtat tatttgataa tccaagaaat
#>       781 gcctcctact tctggttcaa gtagaaatga aaatagaaga attacaagaa tatttagaaa
#>       841 aagatagatc tcggcaacaa cacttcctat acccactttt ctttcaggag tatatttatg
#>       901 cacttgctca tgattatggg tttaaagggt tcgatttttt acgaacctat ggaaattggg
#>       961 ggttatgata ataaatctag ttcagtactt gtaaaacatt taattactcg aatgtatcaa
#>      1021 cagaattatt tgatttattc tgttaatgaa tctaaccaaa atcgattgat tgagcataac
#>      1081 aattcttttt attctcaaat gatatctgaa gtttttgcga tcattgcaga aattccattc
#>      1141 tctcagcaat tactattttc tcttcgagga aaaaagaata ccaaaatctc agactttacg
#>      1201 atctattcat tcaatatttc cctttttaga agacaaatta tcacatttaa actatgtgtc
#>      1261 agatatatta ataccctatc ccatccattt ggaaatcttg gtgcaaattc ttcaatgctg
#>      1321 gatccaagat gtttcttctt tgcatttatt gcgattcttt ctccacgaac atcataatgg
#>      1381 gaatagtttt ttttttccaa agaaatcctt ttcaaaagaa aataaaagac tctttcgatt
#>      1441 cctatataat tcttatgtat ctgaatgtga atttgtctta gtgtttcttc gtaaacaatc
#>      1501 ctcttattta caatcaaaat cctatggaat ctttcttgag cgaacacatt tctatggaag
#>      1561 aatggaacat cttatagtag tgtgtcataa ttattgtcag aaggcctttt gggtcttcaa
#>      1621 ggatcctttt atgcattatg ttcgatatca aggaaaagca attctggcat caaaaggatc
#>      1681 ttatcttttg atgaagaaat ggagatgtca tcttgtcaat ttctggcaat attattttca
#>      1741 tttttgggct cagccttaca gaatttcaat aaaccaatta ggaaatcatt ccttctattt
#>      1801 tctcggttat ctttcaagtg tattaaaaaa tacttcgtct gtaaggaatc aaatgctaga
#>      1861 gaattccttt ttaatagata ctattactaa taaattggat accatagtcc cagttcttcc
#>      1921 tcttattgga tctttgtcta aagctaaatt ttgtaccgta tccgggcatc ctagtagtaa
#>      1981 gccaatctgg acggatttat cggattctga tattattgat agatttggtc ggatatgtag
#>      2041 aaatctttct cattattata gtggatcctc aaaaaaacag agcttatatc gaataaggta
#>      2101 tatacttcga ctttcttgtg ctagaacttt agctcgtaaa cataaaagta cagtacgtgc
#>      2161 ttttttgcaa agattaggtt cggaattatt agaagaattc tttacagaag aagaaggagt
#>      2221 tgtttttttg atttcccaaa agaacaaaac ctcttttcct ctctataggt cacatagaga
#>      2281 acgcatttgg tatttggata ttatccatat taatgaattg gtgaattcat ttatgatggg
#>      2341 gcgataagcc cctataaaat aagaaatata aattttttct aatgtctaat aaatagacga
#>      2401 caaattcatt aattttcatt ctgaaatgct catctagtag tgtagtgatt gaatcaactg
#>      2461 agtattcaaa atttttagac aaacttctag ggatagaagt ttgttttatc tgtatacata
#>      2521 ggtaaagtcg tgtgcaatga aaaatgcaag cacgatttgg ggagagataa ttttctctat
#>      2581 tgtaacaaat aaaaattatc tactccatcc gactagttaa tcg
#> //
```

## Extracting

We can extract different elements of the above record with `gb_extract()` ....


```r
# such as the LOCUS information ...
(gb_extract(record = record, what = 'locus'))
#>     accession        length           mol          type        domain          date 
#>    "AY952423"        "2623"         "DNA"      "linear"         "PLN" "17-APR-2005"
# the accession
(gb_extract(record = record, what = 'accession'))
#> [1] "AY952423"
# the accession + version
(gb_extract(record = record, what = 'version'))
#> [1] "AY952423.1"
# the organism name
(gb_extract(record = record, what = 'organism'))
#> [1] "Livistona chinensis"
# the sequence definition line
(gb_extract(record = record, what = 'definition'))
#> [1] "Livistona chinensis tRNA-Lys (trnK) gene, partial sequence; and matK gene, complete sequence; chloroplast"
# the keywords (this record doesn't have any ....)
(gb_extract(record = record, what = 'keywords'))
#> [1] ""
# even the features as a list object
features <- gb_extract(record = record, what = 'features')
print(features[[1]])
#> $type
#> [1] "source"
#> 
#> $location
#> [1] "1..2623"
#> 
#> $organism
#> [1] "Livistona chinensis"
#> 
#> $organelle
#> [1] "plastid:chloroplast"
#> 
#> $mol_type
#> [1] "genomic DNA"
#> 
#> $db_xref
#> [1] "taxon:115492gene            <1..>2623"
#> 
#> $gene
#> [1] "trnK"
#> 
#> $note
#> [1] "tRNA-Lysintron          <1..>2623"
# and of course the sequence itself
seq <- gb_extract(record = record, what = 'sequence')
str(seq)
#>  chr "ATTGGGGTTGCTAACTCAACGGTAGAGTACTCGGCTTTTAAGTGCGACTATCATCTTTTACACATTTGGATGAAGTAAGGAATTCGTCCAGACTATTGGTAGAGTCTATAA"| __truncated__
```

## From the database

You can try out the above functions yourself on any sequence record by downloading them through the [`rentrez` package](https://github.com/ropensci/rentrez) using `entrez_fetch(db='nucleotide', rettype='gb')`. Or why not test them out using any of the records from the rodents database?



```r
library(restez)
restez_path_set(rodents_path)
restez_connect()
#> Remember to run `restez_disconnect()`
(rand_id <- sample(suppressWarnings(list_db_ids()), 1))
#> [1] "AB008118"
record <- gb_record_get(rand_id)
(gb_extract(record = record, what = 'features'))
#> [[1]]
#> [[1]]$type
#> [1] "source"
#> 
#> [[1]]$location
#> [1] "1..714"
#> 
#> [[1]]$organism
#> [1] "Mus musculus"
#> 
#> [[1]]$mol_type
#> [1] "genomic DNA"
#> 
#> [[1]]$strain
#> [1] "129"
#> 
#> [[1]]$db_xref
#> [1] "taxon:10090"
#> 
#> 
#> [[2]]
#> [[2]]$type
#> [1] "gene"
#> 
#> [[2]]$location
#> [1] "complement(380..510)"
#> 
#> [[2]]$gene
#> [1] "Limk-2"
#> 
#> 
#> [[3]]
#> [[3]]$type
#> [1] "exon"
#> 
#> [[3]]$location
#> [1] "complement(380..510)"
#> 
#> [[3]]$gene
#> [1] "Limk-2"
#> 
#> [[3]]$note
#> [1] "1a"
restez_disconnect()
```

## Next up

**[Running phylotaR with restez](https://ropensci.github.io/restez/articles/4_phylotar.html)**
---
title: "5. Tips and Tricks"
output: rmarkdown::html_vignette
---



## Multiple restez paths

It is not advisable to download the entire GenBank database to your machine. Equally, it is best to limit the size of a database. Databases that are too large will be slow to query and are more likely to cause memory issues. For example, you may actually make a query that demands more memory than is available on your machine. One solution to instead set multiple `restez` paths on your machine.

You can either set up a path for different domains. Or you could download for a single set of domains and then create a database from the same downloaded files using the `alt_restez_path` argument. Do also make use of `restez_path_unset` to disconnect and unset the `restez` path.


```r
# a larger database from the same download files in rodents_path
db_create(alt_restez_path = rodents_path, max_length = 2000)
```


## Connecting and disconnecting

Always ensure you disconnect after connecting to a `restez` path. Not doing so may lead to some strange database errors such as 'seg faults' or you may even be prevented from connecting to a database again until you restart R. In scripts you should always place `restez_disconnect()` as the end of the script or when you have stopped making queries. If you are making queries from your own custom function you should use `on.exit`. This allows you to run 'clean up' code whenever a function exits, even if it errors.



```r
suppressMessages(library(restez))
random_definition <- function() {
  suppressMessages(restez_connect())
  on.exit(restez_disconnect())
  if (restez_ready()) {
    # deliberate mistake
    id <- sample(list_db_ids(n = NULL), 1)[[1]]
    return(gb_definition_get(id))
  }
}
restez_path_set(rodents_path)
(definition <- random_definition())
```

```
##                                                                                              KU614570 
## "Mus musculus clone PD151104P3E10 immunoglobulin heavy chain variable region (Igh) mRNA, partial cds"
```

```r
# not connected outside of function!
(restez_ready())
```

```
## [1] FALSE
```

## Which domain?

The `db_download` function lists the various possible GenBank domains that can be downloaded. You can work out which GenBank domain a sequence belongs to by its three letter code towards the end of its locus. For example, the top of the record for this sequence indicates it is in the rodent domain.

```
LOCUS       LT548182                 456 bp    DNA     linear   ROD 23-NOV-2016
DEFINITION  TPA_inf: Cavia porcellus GLNH gene for globin H.
ACCESSION   LT548182
VERSION     LT548182.1
```

## Database performance and behaviour

The `restez` package database is built with [`MonetDBlite`](https://github.com/hannesmuehleisen/MonetDBLite-R).
If you encounter any errors that include the phrase "Server says", then an issue is
likely to have occurred within the database. Please raise such issues with
[GitHub](https://github.com/ropensci/restez/issues). But keep the following
factors in mind:

* Is your request from the database likely to return an object too large for
your computer's RAM? If the size of database is 5GB then it is likely that
a request pulling all of the sequence data and information into an R session
will be around 5GB as well.
* Are you building and storing the database on a separate USB drive? It has
been noted that database behaviour can be unusual on separate USB drives. When
an issue, please provide information about your USB drive's format, size and USB
connections.
---
title: "Create and Query a Local Copy of GenBank in R"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{restez_tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# hide line spinner
assignInNamespace(x = 'custom_download2', value = restez:::custom_download,
                  ns = 'restez')
```

# 1. Download GenBank

In the first section of this introduction to `restez` we will explain how to download parts of GenBank on to your computer and from these downloads, construct a local database. GenBank can be downloaded from NCBI using the File Transfer Protocol (FTP), the online portal for these downloads can be viewed from your web browser,
[FTP NCBI GenBank](ftp://ftp.ncbi.nlm.nih.gov/genbank/). GenBank is hosted at this portal as 1000s of compressed sequence files called flatfiles. New versions of these flatfiles are released once a month as new sequence information is uploaded. For information on the latest GenBank release, see the [NCBI statistics page](https://www.ncbi.nlm.nih.gov/genbank/statistics/).

`restez` downloads a selection of these flatfiles and unpacks them into an SQL-like database (specifically, [MonetDB](https://www.monetdb.org/)) that can be queried using `restez`'s functions. Because there are potentially huge amounts of sequence information, all of which is unlikely to be of interest to a user, `restez` allows you to select the sets of flatfiles to download based on GenBank's domains. For the most part these sets of flatfile types are determined by taxonomy, for example, a user can download all primate or plant sequences.

Here, we will explain how to set up `restez`, download genbank and then create a database from the downloaded files.

## 1.1 Setting up

Before we can download anything, we need to tell `restez` where we would like to store the downloaded files and the eventual GenBank database. To do this we use the `restez_path_set` function. With this function a user can provide a system file path to a folder within which the function will create a new folder called 'restez'. In this example, we will create a `restez` folder in a temporary folder.

```{r presetup, include=FALSE}
rstz_pth <- tempdir()
restez::restez_path_set(filepath = rstz_pth)
restez::db_delete(everything = TRUE)
```
```{r setting up, echo=TRUE}
library(restez)
# rstz_pth <- tempdir()
restez_path_set(filepath = rstz_pth)
```

A user must always set the `restez` path in order to use the `restez` library and we are reminded to do this whenever we load the package. (We are also reminded to run `restez_connect` but we'll come to that a bit later.) Additionally, by running `restez_path_set` on a new folder messages are printed to console telling us that the path was created and that a downloads subfolder was created.

> **Note**, whenever you intend to use `restez` you will need to specify the `restez` path. If you are using it regularly, it might be a good idea to set it using your .Rprofile.

## 1.2 Download

Now we can download GenBank files using `db_download()`. This is an interactive function, it looks up the latest GenBank release, parses the release notes and prints to console the available sets of flatfiles that can be downloaded. For this example, we will select the smallest available domain which is 'Unannotated' and can be pre-specified with '20'. Normally, you would run this function without the `preselection` argument and just wait for the function to prompt you for a domain selection.

After launching, the download function will ask you whether it is OK to download and set up a database for the selected sequence files. If the selection is likely to be large you can always quit the process using either `Esc` or `Ctrl+c`. Please note that the diskspace estimates are preliminary. The package trys to be conservative in its estimates.

```{r download, echo=TRUE}
db_download(preselection = '20')
```

Now the download has begun, we just need to wait for it to finish. For the largest domains, this can take quite a while and may be prone to server errors. If you ever encounter a failed downloaded file it is possible to re-run the `db_download` function again and -- so long as the same selection is made -- the function will only download the missing files.

## 1.3 Create a database

After the download has completed, we need to use `db_create()` to create our local database. This looks in the downloads folder of our `restez` path, breaks these files up into separate GenBank records and adds them to the SQL-like database. Again, for very large collections of sequences this can take quite a while.

Before running the `db_create` function, however, we must first establish a connection to the database in the `restez` path using the `restez_connect` function. As with SQL databases we must always make sure to connect and then disconnect after we have interacted with a database.

```{r create database, echo=TRUE}
restez_connect()
db_create()
restez_disconnect()
```

`db_create()` allows a user to specify minimum and maximum sequence lengths. It's always a good idea to limit the number of sequences in a database so that look-up times are faster. If you know you are only interested in certain lengths of sequences it is a good idea to limit the sequences in the database at this stage. You can always run `db_create()` again to change the limits. You will simply need to delete the original database first with `db_delete`.

## 1.4 Checking the setup

After the download and the database creation steps are complete, we can confirm our setup using `restez_status`.

```{r confirm setup, echo=TRUE}
library(restez)
restez_path_set(rstz_pth)
restez_connect()
restez_status()
restez_disconnect()
```

The status function allows a user to always touch base with their `restez` set-up. It can be run at any point, before and/or after download or database creation. It is also designed to provide useful messages to a user for what they need to do next in order for them to make queries to their database. If you ever get confused, run `restez_status` to see what is happening.

Additionally, if you are developing your own functions that depend on `restez`, you can use `restez_ready` which will simply return TRUE or FALSE on whether the database can be queried.

```{r confirm setup2, echo=TRUE}
library(restez)
restez_path_set(rstz_pth)
restez_connect()
if (restez_ready()) {
  print('Database is ready to be queried!')
} else {
  print('Database cannot be queried :-(')
}
restez_disconnect()
```

# 2. Query

Once a restez database has been set up we can query the database using the `gb_*_get()` functions. These functions allow us to retrieve specific columns in the SQL-like database: 'sequence', 'definition', 'accession', 'version' and 'organism'. Also, they allow us to get the whole text formatted record and sequence data in fasta format. In this example, we can use `list_db_ids()` to identify accession numbers in the database. We could also use `entrez_search()` provided the database contains sequences of interest to us, see ['search and fetch'](https://ropensci.github.io/restez/articles/2_search_and_fetch.html).

> Always remember to disconnect after you have finished your querying!

```{r query, echo=TRUE}
library(restez)
restez_path_set(rstz_pth)
restez_connect()
ids <- suppressWarnings(list_db_ids(db = 'nucleotide', n = 100))
(id <- sample(ids, 1))
# sequence
str(gb_sequence_get(id))
# definition
(gb_definition_get(id))
# version
(gb_version_get(id))
# organism
(gb_organism_get(id))
# fasta
cat(gb_fasta_get(id))
# Note, for large databases these requests can take a long time to complete.
# Always try and limit the size of the database by specifying min and max
# sequence lengths with db_create()
# Also note, if an id is not present in the database nothing is returned
cat(gb_fasta_get(id = c(id, 'notanid')))
# Always remember to disconnect
restez_disconnect()
```

Additionally, for more flexibility and options for extracting sequence record information see [GenBank record parsing](https://ropensci.github.io/restez/articles/3_parsing.html).

# 3. Entrez

Entrez wrappers are part of the `restez` package. These allow a user to make use of the local GenBank using functions that were built for [`rentrez`](https://github.com/ropensci/rentrez). This minimises the amount of coding changes required for any Entrez dependent code.

>Currently, only `entrez_fetch()` is available with restez and only text formatted rettypes are allowed.

```{r entrez, echo=TRUE}
library(restez)
restez_path_set(rstz_pth)
restez_connect()
ids <- suppressWarnings(list_db_ids(db = 'nucleotide', n = 100))
(id <- sample(ids, 1))
# get fasta sequences with entrez
res <- restez::entrez_fetch(db = 'nucleotide', id = id, rettype = 'fasta')
cat(res)
# entrez_fetch will also search via Entrez for any ids not in the db
plant_sequence <- 'AY952423'
res <- restez::entrez_fetch(db = 'nucleotide', id = c(id, plant_sequence),
                            rettype = 'fasta')
cat(res)
restez_disconnect()
```

# 4. More information

For more information about `restez` see the other tutorials:

1. [Build a database of all rodents](https://ropensci.github.io/restez/articles/1_rodents.html)
2. [How to search for and fetch sequences](https://ropensci.github.io/restez/articles/2_search_and_fetch.html)
3. [Advanced parsing of a GenBank record](https://ropensci.github.io/restez/articles/3_parsing.html)
4. [Running phylotaR with restez](https://ropensci.github.io/restez/articles/4_phylotar.html)
5. [Tips and tricks](https://ropensci.github.io/restez/articles/5_tips_and_tricks.html)
---
title: "2. How to search for and fetch sequences"
date: "2020-01-07"
output: rmarkdown::html_vignette
---



A downside of the `restez` approach is we are unable to then search the local database for sequences of interest as the database cannot possibly contain sufficient metadata (most notably taxonomic information) to perform as rich a search as online via NCBI. As a result we must perform the sequence discovery process independently from `restez`. In this tutorial we will demonstrate how to perform online searches using the `rentrez` package and then use the search results to look up sequences in a `restez` database.

To run all the code in this tutorial you will need to have already set up the rodents database, see [Build a database of all rodents](https://ropensci.github.io/restez/articles/1_rodents.html).

## Search NCBI

Let's pretend we're interested in COI sequences for a group of rodents called the Sciuromorpha. We can create an NCBI search term; use `entrez_search` to perform the search; and then retrieve the accession IDs using `entrez_fetch` with the `acc` rettype. `entrez_fetch` returns a single text that will need to be split up by the newline character. Finally, we should drop any version number after the downloaded accessions. This is not strictly necessary, but it makes it easier to check our results from `restez` later.



```r
# 33553 - Sciuromorpha - squirrel-like things
search_term <- 'txid33553[Organism:exp] AND COI [GENE] AND 100:1000[SLEN]'
search_object <- rentrez::entrez_search(db = 'nucleotide', term = search_term,
                                        use_history = TRUE, retmax = 0)
accessions <- rentrez::entrez_fetch(db = 'nucleotide',
                                    web_history = search_object$web_history,
                                    rettype = 'acc')
accessions <- strsplit(x = accessions, split = '\\n')[[1]]
accessions <- sub(pattern = '\\.[0-9]+', replacement = '', x = accessions)
print(length(accessions))
#> [1] 462
print(accessions[1:10])
#>  [1] "MN326076" "MN326074" "MN326066" "MN326062" "MN326059" "MN326045" "HM380207" "KY033226" "KX859270" "KX859268"
```

## Retrieve sequences

To fetch the sequences from the rodents database, we can just use the `gb_fasta_get` function.



```r
library(restez)
restez_path_set(rodents_path)
restez_connect()
#> Remember to run `restez_disconnect()`
coi_sequences <- gb_fasta_get(id = accessions)
str(coi_sequences[[1]])
#>  chr ">FJ808614.1 Glis glis isolate A1g02_E11-F cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial\n"| __truncated__
# Are all accessions in results?
all(accessions %in% names(coi_sequences))
#> [1] TRUE
```

## Comparing to Entrez

Can we not just use `rentrez` to do the fetching as well? Yes, but `restez` can be a lot faster. NCBI limits the number of requests per user, often to as little as 100 items per request with varying time delays. Additionally for any programmatic retrieval of sequences using an online server can never be as reliable as a local copy.


```r
# time via restez
system.time(expr = {
  coi_sequences <- gb_fasta_get(id = accessions)
  })
#>    user  system elapsed 
#>   0.475   0.967   1.748
# time via Entrez
system.time(expr = {
  coi_sequences_p1 <- rentrez::entrez_fetch(db = 'nucleotide',
                                            id = accessions[1:100],
                                            rettype = 'fasta')
  coi_sequences_p2 <- rentrez::entrez_fetch(db = 'nucleotide',
                                            id = accessions[101:200],
                                            rettype = 'fasta')
  coi_sequences_p3 <- rentrez::entrez_fetch(db = 'nucleotide',
                                            id = accessions[201:300],
                                            rettype = 'fasta')
  coi_sequences_p4 <- rentrez::entrez_fetch(db = 'nucleotide',
                                            id = accessions[301:400],
                                            rettype = 'fasta')
  coi_sequences_p5 <- rentrez::entrez_fetch(db = 'nucleotide',
                                            id = accessions[401:456],
                                            rettype = 'fasta')
  })
#>    user  system elapsed 
#>   0.136   0.037   7.207
# always disconnect
restez_disconnect()
```
<!-- Below is no longer relevant now that the size of sequences in the db has been limited.
## Missing

A user should know that if an ID cannot be found in the local database no error or warning is raised. This is why it can be good practice to test whether all the provided IDs are in the returned named vector. In this example, we can see that not all the accession IDs that were provided are in the returned `coi_sequences`. Why is that?


```r
# Are all accessions in results?
all(accessions %in% names(coi_sequences))
#> [1] TRUE
# .... no
```

This is because `restez` only downloads GenBank sequences, and these 'missing' sequence IDs are RefSeq sequences. GenBank, however, acts as the source database for RefSeq -- which is the case for other NCBI databases too -- and all the missing RefSeq sequences can also be found in GenBank under a different ID.

We can look up alternative IDs and test whether they are in our `accesssions` vector, so:


```r
(accessions[!accessions %in% names(coi_sequences)])
#> character(0)
# NC* refers to RefSeq sequences and are not available through restez
# The sequence exists in GB under a different id which we can find like so
smmry <- rentrez::entrez_summary(db = 'nucleotide', id = 'NC_027278')
# This ID does exist in our results.
(smmry$assemblyacc %in% accessions)
#> [1] FALSE
```
-->
## Next up

**[Advanced parsing of a GenBank record](https://ropensci.github.io/restez/articles/3_parsing.html)**
---
title: "4. Running phylotaR with restez"
output: rmarkdown::html_vignette
---



In this tutorial we will showcase how a `restez` database can be used to speed up a [phylotaR](https://github.com/ropensci/phylotaR) run. `phylotaR` runs an automated pipeline for identifying ortholgous gene clusters as the first step in a phylogenetic analysis. A user provides a taxonomic identity and the pipeline downloads all relevant sequences and identifies clusters using a local-alignment search tool. For more information on `phylotaR` see its [published article](https://doi.org/10.3390/life8020020).

By using `restez` in conjunction with `phylotaR`, we will not only being saving time, but also improving the chances of a successful `phylotaR` run -- often NCBI Entrez limits the number of requests or even rejects requests from IP addresses that are making too many. Note, however, that the gains in using `restez` with `phylotaR` only make sense if you make use of the `restez` database multiple times or if you wish to radically increase the maximum number of sequences to download per taxon (by default it is only 3,000). Also, note that using a `restez` database does not currently eliminate the need for an internet connection. `phylotaR` still needs to look up taxonomic information and must also identify relevant sequence IDs using Entrez (this may change in the future as `restez` develops).

We will run a `phylotaR` run for the rodent subfamily, Dipodinae, the clade that contains jerboas. Because this clade is so small, it should not take too long to run it. We will, however, limit to showcasing only the download stage of the `phylotaR` pipeline, where `restez` is required.

## Install phylotaR

```r
# Currently only the development version of phylotaR can work with restez
# It must be installed from GitHub using devtools
devtools::install_github(repo = 'ropensci/phylotaR')
```

## Setup
Since we will be running the `phylotaR` pipeline for Dipodinae, we can use the rodents database we created [before](https://ropensci.github.io/restez/articles/1_rodents.html). We do not need to set-up the `phylotaR` pipeline any differently with a `restez` database, except to ensure we have set the `restez` path to the rodent database. (Connections to the `restez` database will be made automatically, `restez_connect` or `restez_disconnect` do not need to be run.)

```r
library(phylotaR)

# Restez (no need to call package)
restez::restez_path_set(filepath = rodents_path)

# Vars
wd <- 'dipodinae'
dir.create(wd)
txid <- 35737  # Dipodinae
mxsql <- 500
ncbi_dr <- '[PATH/TO/BLAST]' # e.g. '/usr/local/ncbi/blast/bin'

# setup
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mxsql = 500)
```

## Run


```r
# run just the first two stages for this demonstration
taxise_run(wd)
#> --------------------------------------------
#> Starting stage TAXISE: [2020-01-07 09:26:04]
#> --------------------------------------------
#> Searching taxonomic IDs ...
#> Downloading taxonomic records ...
#> . [1-28]
#> Generating taxonomic dictionary ...
#> ---------------------------------------------
#> Completed stage TAXISE: [2020-01-07 09:26:08]
#> ---------------------------------------------
download_run(wd)
#> ----------------------------------------------
#> Starting stage DOWNLOAD: [2020-01-07 09:26:08]
#> ----------------------------------------------
#> Identifying suitable clades ...
#> Identified [1] suitable clades.
#> Connecting to restez database ...
#> Downloading hierarchically ...
#> Working on parent [id 35737]: [1/1] ...
#> . + whole subtree ...
#> . . Getting [100 sqs] from restez database...
#> Successfully retrieved [100 sqs] in total.
#> -----------------------------------------------
#> Completed stage DOWNLOAD: [2020-01-07 09:26:12]
#> -----------------------------------------------
```

Sequences that cannot be found locally will be downloaded via Entrez. Sequences may not be found locally for three sets of reasons, 1. `phylotaR` may have identified non-GenBank sequences (e.g. RefSeq), 2. the local database may have size limits (ours has a min of 100 and a max 1000), and 3. the local database is out of date, new GenBank releases are only made once a month.

> Note: If the `phylotaR` messages do not explicit state that they are using a `restez` database, then, you may not have set the `restez` path.

Compared to running without a `restez` database, phylotaR download can run [**2x** faster](https://github.com/ropensci/restez/blob/master/other/phylotar_demo.R).

## Next up

**[Tips and tricks](https://ropensci.github.io/restez/articles/5_tips_and_tricks.html)**
---
title: "1. Build a database of all rodents"
date: "2020-01-07"
output: rmarkdown::html_vignette
---



In this first tutorial we are going to build a database for all rodents. The rodents are a good test case for playing with `restez` as they are a relatively small domain in GenBank but still have charismatic organisms that people are familiar enough with to understand. To keep things extra fast, we will also limit the number of sequences in the database by limiting the sequence sizes between 100 and 1000.

The database you build here will be used again in later tutorials and you may wish to experiment with it yourself. Therefore it is best to locate a suitable place in your harddrive where you would like to store it for later reference. In this tutorial and in others, we will always refer to the rodents' `restez` path with the variable `rodents_path`.

Setting up the rodents database will likely take a long time. The exact time will depend on your internet speeds and machine specs. For reference, this vigenette was written on a MacBook Air (2013) via WiFi with a download speed of 13 MBPS. With this setup, downloading the database took 26 minutes and creating the database took 59 minutes.

## Download


```r
library(restez)
# set the restez path to a memorable location
restez_path_set(rodents_path)
# download for domain 15
db_download(preselection = '15')
```

## Build

```r
library(restez)
restez_path_set(rodents_path)
db_create(min_length = 100, max_length = 1000)
```

## Check status

```r
library(restez)
#> -------------
#> restez v1.0.2
#> -------------
#> Remember to restez_path_set() and, then, restez_connect()
#> 
#> Attaching package: 'restez'
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     record
restez_path_set(rodents_path)
restez_connect()
#> Remember to run `restez_disconnect()`
restez_status()
#> Checking setup status at  ...
#> ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Restez path ...
#> ... Path '[RODENTS PATH]/restez'
#> ... Does path exist? 'Yes'
#> ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Download ...
#> ... Path '[RODENTS PATH]/restez/downloads'
#> ... Does path exist? 'Yes'
#> ... N. files 35
#> ... N. GBs 4.21
#> ... GenBank division selections 'Rodent'
#> ... GenBank Release 235
#> ... Last updated '2020-01-03 12:45:46'
#> ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Database ...
#> ... Path '[RODENTS PATH]/restez/sql_db'
#> ... Does path exist? 'Yes'
#> ... N. GBs 0.68
#> ... Is database connected? 'Yes'
#> ... Does the database have data? 'Yes'
#> ... Number of sequences 223200
#> ... Min. sequence length 100
#> ... Max. sequence length 1000
#> ... Last_updated '2020-01-03 13:39:35'
restez_disconnect()
```

## Next up

**[How to search for and fetch sequences](https://ropensci.github.io/restez/articles/2_search_and_fetch.html)**
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{restez_path_check}
\alias{restez_path_check}
\title{Check restez filepath}
\usage{
restez_path_check()
}
\description{
Raises error if restez path does
not exist.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_by_patterns}
\alias{extract_by_patterns}
\title{Extract by keyword}
\usage{
extract_by_patterns(record, start_pattern, end_pattern = "\\n")
}
\arguments{
\item{record}{GenBank record in text format, character}

\item{start_pattern}{REGEX pattern indicating the point to
start extraction, character}

\item{end_pattern}{REGEX pattern indicating the point to
stop extraction, character}
}
\value{
character or NULL
}
\description{
Search through GenBank record for a keyword and
return text up to the end_pattern.
}
\details{
The start_pattern should be any of the capitalized elements
in a GenBank record (e.g. LOCUS, DESCRIPTION, ACCESSION).
The end_pattern depends on how much of the selected element
a user wants returned. By default, the extraction will stop
at the next newline.
If keyword or end pattern not found, returns NULL.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_clean_sequence}
\alias{extract_clean_sequence}
\title{Extract clean sequence from sequence part}
\usage{
extract_clean_sequence(seqrecpart)
}
\arguments{
\item{seqrecpart}{Sequence part of a GenBank record, character}
}
\value{
character
}
\description{
Return clean sequence from seqrecpart of a GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-get-tools.R
\name{list_db_ids}
\alias{list_db_ids}
\title{List database IDs}
\usage{
list_db_ids(db = "nucleotide", n = 100)
}
\arguments{
\item{db}{character, database name}

\item{n}{Maximum number of IDs to return, if NULL returns all}
}
\value{
vector of characters
}
\description{
Return a vector of all IDs in
a database.
}
\details{
Warning: can return very large vectors
for large databases.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
# Warning: not recommended for real databases
#  with potentially millions of IDs
restez_connect()
all_ids <- list_db_ids()


# What shall we do with these IDs?
# ... how about make a mock fasta file
seqs <- gb_sequence_get(id = all_ids)
defs <- gb_definition_get(id = all_ids)
# paste together
fasta_seqs <- paste0('>', defs, '\\n', seqs)
fasta_file <- paste0(fasta_seqs, collapse = '\\n')
cat(fasta_file)


# delete after example
db_delete(everything = TRUE)
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_create}}, \code{\link{db_delete}},
  \code{\link{db_download}}, \code{\link{demo_db_create}},
  \code{\link{is_in_db}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_fasta_get}
\alias{gb_fasta_get}
\title{Get fasta from GenBank}
\usage{
gb_fasta_get(id, width = 70)
}
\arguments{
\item{id}{character, sequence accession ID(s)}

\item{width}{integer, maximum number of characters in a line}
}
\value{
named vector of fasta sequences, if no results found NULL
}
\description{
Get sequence and definition data in FASTA format. Equivalent to
\code{rettype='fasta'} in \code{\link[rentrez]{entrez_fetch}}.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(fasta <- gb_fasta_get(id = 'demo_1'))
(fastas <- gb_fasta_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_definition_get}},
  \code{\link{gb_organism_get}},
  \code{\link{gb_record_get}},
  \code{\link{gb_sequence_get}},
  \code{\link{gb_version_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomartr-tools.R
\name{custom_download}
\alias{custom_download}
\title{Helper function to perform customized downloads}
\usage{
custom_download(...)
}
\arguments{
\item{...}{additional arguments that shall be passed to
\code{\link[downloader]{download}}}
}
\description{
To achieve the most stable download experience,
ftp file downloads are customized for each operating system.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\author{
Hajk-Georg Drost
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{db_sqlngths_log}
\alias{db_sqlngths_log}
\title{Log the min and max sequence lengths}
\usage{
db_sqlngths_log(min_lngth, max_lngth)
}
\arguments{
\item{min_lngth}{Minimum length}

\item{max_lngth}{Maximum length}
}
\description{
Log the min and maximum sequence length used in the created db.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{db_sqlngths_get}
\alias{db_sqlngths_get}
\title{Return the minimum and maximum sequence lengths in db}
\usage{
db_sqlngths_get()
}
\value{
vector of integers
}
\description{
Returns the maximum and minimum sequence lengths as set by the
user upon db creation.
}
\details{
If no file found, returns empty character vector.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{restez_ready}
\alias{restez_ready}
\title{Is restez ready?}
\usage{
restez_ready()
}
\value{
Logical
}
\description{
Returns TRUE if a restez SQL database is available, connected
and has data. Use restez_status() for more information.
}
\examples{
library(restez)
fp <- tempdir()
restez_path_set(filepath = fp)
demo_db_create(n = 5)
restez_connect()
(restez_ready())
restez_disconnect()
db_delete(everything = TRUE)
(restez_ready())
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_disconnect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_set}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_version_get}
\alias{gb_version_get}
\title{Get version from GenBank}
\usage{
gb_version_get(id)
}
\arguments{
\item{id}{character, sequence accession ID(s)}
}
\value{
named vector of versions, if no results found NULL
}
\description{
Return the accession version
for an accession ID.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(ver <- gb_version_get(id = 'demo_1'))
(vers <- gb_version_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_definition_get}},
  \code{\link{gb_fasta_get}},
  \code{\link{gb_organism_get}},
  \code{\link{gb_record_get}},
  \code{\link{gb_sequence_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{restez_connect}
\alias{restez_connect}
\title{Connect to the restez database}
\usage{
restez_connect()
}
\description{
Sets a connection to the local database. If database
connection cannot be made, an error is returned.
}
\examples{
library(restez)
fp <- tempdir()
restez_path_set(filepath = fp)
demo_db_create(n = 5)
restez_connect()
restez_status()
restez_disconnect()
db_delete(everything = TRUE)
# Errors
# restez_status()
}
\seealso{
Other setup: \code{\link{restez_disconnect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_set}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_ready}}, \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_organism}
\alias{extract_organism}
\title{Extract organism}
\usage{
extract_organism(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return organism name from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{restez_path_set}
\alias{restez_path_set}
\title{Set restez path}
\usage{
restez_path_set(filepath)
}
\arguments{
\item{filepath}{character, valid filepath to the folder where the
database should be stored.}
}
\description{
Specify the filepath for the local GenBank database.
}
\details{
Adds 'restez_path' to options(). In this path
the folder 'restez' will be created and all downloaded and
database files will be stored there.
}
\examples{
\dontrun{
library(restez)
restez_path_set(filepath = 'path/to/where/you/want/files/to/download')
}
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_disconnect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_ready}}, \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_sql_query}
\alias{gb_sql_query}
\title{Query the GenBank SQL}
\usage{
gb_sql_query(nm, id)
}
\arguments{
\item{nm}{character, column name}

\item{id}{character, sequence accession ID(s)}
}
\value{
data.frame
}
\description{
Generic query function for retrieving
data from the SQL database for the get functions.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{gbrelease_get}
\alias{gbrelease_get}
\title{Get the GenBank release number in the restez path}
\usage{
gbrelease_get()
}
\value{
character
}
\description{
Returns the GenBank release number. Returns empty character
if none found.
}
\details{
If no file found, returns empty character vector.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_log}}, \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{gb_extract}
\alias{gb_extract}
\title{Extract elements of a GenBank record}
\usage{
gb_extract(record, what = c("accession", "version", "organism",
  "sequence", "definition", "locus", "features", "keywords"))
}
\arguments{
\item{record}{GenBank record in text format, character}

\item{what}{Which element to extract}
}
\value{
character or list of lists (what='features') or named character
vector (what='locus')
}
\description{
Return elements of GenBank record e.g. sequence, definition ...
}
\details{
This function uses a REGEX to extract particular elements of a
GenBank record. All of the what options return a single character with the
exception of 'locus' or 'keywords' that return character vectors and
'features' that returns a list of lists for all features.


The accuracy of these functions cannot be guaranteed due to the enormity of
the GenBank database. But the function is regularly tested on a range of
GenBank records.

Note: all non-latin1 characters are converted to '-'.
}
\examples{
library(restez)
data('record')
(gb_extract(record = record, what = 'locus'))
}
\concept{parse}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callr-tools.R
\name{gb_build2}
\alias{gb_build2}
\title{Callr version of gb_build()}
\usage{
gb_build2(dpth, seq_files, max_length, min_length)
}
\arguments{
\item{dpth}{Download path (where seq_files are stored)}

\item{seq_files}{.seq.tar seq file names}

\item{max_length}{Maximum sequence length.}

\item{min_length}{Minimum sequence length.}
}
\value{
Logical
}
\description{
Runs \code{\link{gb_build}} in callr.
This allows the user to kill the process. Additionally, the process will
print spinning dots to indicate it is still active.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build}},
  \code{\link{gb_df_create}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{restez_path_get}
\alias{restez_path_get}
\title{Get restez path}
\usage{
restez_path_get()
}
\value{
character
}
\description{
Return filepath to where the restez
database is stored.
}
\examples{
library(restez)
# set a restez path with a tempdir
restez_path_set(filepath = tempdir())
# check what the set path is
(restez_path_get())
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_disconnect}},
  \code{\link{restez_path_set}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_ready}}, \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-tools.R
\name{identify_downloadable_files}
\alias{identify_downloadable_files}
\title{Identify downloadable files}
\usage{
identify_downloadable_files()
}
\value{
data.frame
}
\description{
Searches through the release notes
for a GenBank release to find all listed .seq files.
Returns a data.frame for all .seq files and their
description.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}}, \code{\link{last_add_get}},
  \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callr-tools.R
\name{custom_download2}
\alias{custom_download2}
\title{Callr version of custom_download()}
\usage{
custom_download2(url, destfile)
}
\arguments{
\item{url}{URL of source file, character.}

\item{destfile}{filepath to where the file should be saved.}
}
\description{
Runs \code{\link{custom_download}} as an independent R process.
This allows the user to kill the process. Additionally, the process will
print spinning dots to indicate it is still active.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-setup-tools.R
\name{gb_df_create}
\alias{gb_df_create}
\title{Create GenBank data.frame}
\usage{
gb_df_create(accessions, versions, organisms, definitions, sequences,
  records)
}
\arguments{
\item{accessions}{character, vector of accessions}

\item{versions}{character, vector of accessions + versions}

\item{organisms}{character, vector of organism names}

\item{definitions}{character, vector of sequence definitions}

\item{sequences}{character, vector of sequences}

\item{records}{character, vector of GenBank records in text format}
}
\value{
data.frame
}
\description{
Make data.frame from columns vectors for
nucleotide entries. As part of gb_df_generate().
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-tools.R
\name{latest_genbank_release_notes}
\alias{latest_genbank_release_notes}
\title{Download the latest GenBank Release Notes}
\usage{
latest_genbank_release_notes()
}
\description{
Downloads the latest GenBank release notes to a user's restez
download path.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{sql_path_get}
\alias{sql_path_get}
\title{Get SQL path}
\usage{
sql_path_get()
}
\value{
character
}
\description{
Return path to where SQL database is stored.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{status_class}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_locus}
\alias{extract_locus}
\title{Extract locus}
\usage{
extract_locus(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
named character vector
}
\description{
Return locus information from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-tools.R
\name{testdatadir_get}
\alias{testdatadir_get}
\title{Get test data directory}
\usage{
testdatadir_get()
}
\description{
Get the folder containing test data.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/status-tools.R
\name{restez_status}
\alias{restez_status}
\title{Check restez status}
\usage{
restez_status(gb_check = FALSE)
}
\arguments{
\item{gb_check}{Check whether last download was from latest GenBank release?
Default FALSE.}
}
\value{
Status class
}
\description{
Report to console current setup status of restez.
}
\details{
Always remember to run \code{\link{restez_connect}} before running
this function. Set gb_check=TRUE to see if your downloads are up-to-date.
}
\examples{
library(restez)
fp <- tempdir()
restez_path_set(filepath = fp)
demo_db_create(n = 5)
restez_connect()
restez_status()
restez_disconnect()
db_delete(everything = TRUE)
# Errors:
# restez_status()
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_disconnect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_set}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_ready}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rentrez-wrappers.R
\name{entrez_fetch}
\alias{entrez_fetch}
\title{Entrez fetch}
\usage{
entrez_fetch(db, id = NULL, rettype, retmode = "", ...)
}
\arguments{
\item{db}{character, name of the database}

\item{id}{vector, unique ID(s) for record(s)}

\item{rettype}{character, data format}

\item{retmode}{character, data mode}

\item{...}{Arguments to be passed on to rentrez}
}
\value{
character string containing the file created
}
\description{
Wrapper for rentrez::entrez_fetch.
}
\details{
Attempts to first search local database with user-specified
parameters, if the record is missing in the database, the function then
calls rentrez::entrez_fetch to search GenBank remotely.

\code{rettype='fasta'} and \code{rettype='gb'} are respectively equivalent to 
\code{\link{gb_fasta_get}} and \code{\link{gb_record_get}}.
}
\note{
It is advisable to call restez and rentrez functions with '::' notation
rather than library() calls to avoid namespace issues. e.g.
restez::entrez_fetch().
}
\section{Supported return types and modes}{

XML retmode is not supported. Rettypes 'seqid', 'ft', 'acc' and 'uilist'
are also not supported.
}

\examples{
library(restez)
restez_path_set(tempdir())
demo_db_create(n = 5)
restez_connect()
# return fasta record
fasta_res <- entrez_fetch(db = 'nucleotide',
                          id = c('demo_1', 'demo_2'),
                          rettype = 'fasta')
cat(fasta_res)
# return whole GB record in text format
gb_res <- entrez_fetch(db = 'nucleotide',
                       id = c('demo_1', 'demo_2'),
                       rettype = 'gb')
cat(gb_res)
# NOT RUN
# whereas these request would go through rentrez
# fasta_res <- entrez_fetch(db = 'nucleotide',
#                           id = c('S71333', 'S71334'),
#                           rettype = 'fasta')
# gb_res <- entrez_fetch(db = 'nucleotide',
#                        id = c('S71333', 'S71334'),
#                        rettype = 'gb')

# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
\code{\link[rentrez]{entrez_fetch}}
}
\concept{entrez}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-tools.R
\name{mock_def}
\alias{mock_def}
\title{Mock def}
\usage{
mock_def(i)
}
\arguments{
\item{i}{integer, iterator}
}
\value{
character
}
\description{
Make a mock sequence definition.
Designed to be part of a loop.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_sequence_get}
\alias{gb_sequence_get}
\title{Get sequence from GenBank}
\usage{
gb_sequence_get(id)
}
\arguments{
\item{id}{character, sequence accession ID(s)}
}
\value{
named vector of sequences, if no results found NULL
}
\description{
Return the sequence(s) for a record(s)
from the accession ID(s).
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(seq <- gb_sequence_get(id = 'demo_1'))
(seqs <- gb_sequence_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_definition_get}},
  \code{\link{gb_fasta_get}},
  \code{\link{gb_organism_get}},
  \code{\link{gb_record_get}}, \code{\link{gb_version_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-tools.R
\name{char}
\alias{char}
\title{Print green}
\usage{
char(x)
}
\arguments{
\item{x}{Text to print, character}
}
\value{
coloured character encoding, character
}
\description{
Print to console green text to indicate a name/filepath/text
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{check_connection}},
  \code{\link{cleanup}}, \code{\link{connected}},
  \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{filename_log}
\alias{filename_log}
\title{Write filenames to log files}
\usage{
filename_log(fl, fp)
}
\arguments{
\item{fl}{file name, character}

\item{fp}{filepath to log file, character}
}
\description{
Record a filename in a log file along with GB release and time.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{flatfile_read}},
  \code{\link{gb_build2}}, \code{\link{gb_build}},
  \code{\link{gb_df_create}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{connected}
\alias{connected}
\title{Is restez connected?}
\usage{
connected()
}
\value{
Logical
}
\description{
Returns TRUE if a restez SQL database has been connected.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-get-tools.R
\name{is_in_db}
\alias{is_in_db}
\title{Is in db}
\usage{
is_in_db(id, db = "nucleotide")
}
\arguments{
\item{id}{character, sequence accession ID(s)}

\item{db}{character, database name}
}
\value{
named vector of booleans
}
\description{
Determine whether an id(s)
is/are present in a database.
}
\examples{
library(restez)
# set the restez path to a temporary dir
restez_path_set(filepath = tempdir())
# create demo database
demo_db_create(n = 5)
restez_connect()
# in the demo, IDs are 'demo_1', 'demo_2' ...
ids <- c('thisisnotanid', 'demo_1', 'demo_2')
(is_in_db(id = ids))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_create}}, \code{\link{db_delete}},
  \code{\link{db_download}}, \code{\link{demo_db_create}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-setup-tools.R
\name{gb_df_generate}
\alias{gb_df_generate}
\title{Generate GenBank records data.frame}
\usage{
gb_df_generate(records, min_length = 0, max_length = NULL)
}
\arguments{
\item{records}{character, vector of GenBank records in text format}

\item{min_length}{Minimum sequence length, default 0.}

\item{max_length}{Maximum sequence length, default NULL.}
}
\value{
data.frame
}
\description{
For a list of records, construct a data.frame
for insertion into SQL database.
}
\details{
The resulting data.frame has five columns: accession,
organism, raw_definition, raw_sequence, raw_record.
The prefix 'raw_' indicates the data has been converted to the
raw format, see ?charToRaw, in order to save on RAM.
The raw_record contains the entire GenBank record in text format.

Use max and min sequence lengths to minimise the size of the database.
All sequences have to be at least as long as min and less than or equal
in length to max, unless max is NULL in which there is no maximum length.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/status-tools.R
\name{print.status}
\alias{print.status}
\title{Print method for status class}
\usage{
\method{print}{status}(x)
}
\arguments{
\item{x}{Status object}
}
\description{
Prints to screen the three sections of the status class.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{quiet_connect}}, \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{last_entry_get}
\alias{last_entry_get}
\title{Return the last entry}
\usage{
last_entry_get(fp)
}
\arguments{
\item{fp}{Filepath, character}
}
\value{
vector
}
\description{
Return the last entry from a tab-delimited log file.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{slctn_log}
\alias{slctn_log}
\title{Log the GenBank selection made by a user}
\usage{
slctn_log(selection)
}
\arguments{
\item{selection}{selected GenBank sequences, named vector}
}
\description{
This function is called whenever a user makes a selection with
the \code{\link{db_download}}. It records GenBank numbers selections.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{sql_path_get}},
  \code{\link{status_class}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{dwnld_path_get}
\alias{dwnld_path_get}
\title{Get dwnld path}
\usage{
dwnld_path_get()
}
\value{
character
}
\description{
Return path to folder where raw .seq files
are stored.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-tools.R
\name{predict_datasizes}
\alias{predict_datasizes}
\title{Print file size predictions to screen}
\usage{
predict_datasizes(uncompressed_filesize)
}
\arguments{
\item{uncompressed_filesize}{GBs of the stated filesize, numeric}
}
\description{
Predicts the file sizes of the downloads and the database
from the GenBank filesize information. Conversion factors are based on
previous restez downloads.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{print.status}},
  \code{\link{quiet_connect}}, \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-tools.R
\name{cleanup}
\alias{cleanup}
\title{Clean up test data}
\usage{
cleanup()
}
\description{
Removes all temporary test data created.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{connected}},
  \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{last_add_get}
\alias{last_add_get}
\title{Return date and time of the last added sequence}
\usage{
last_add_get()
}
\value{
character
}
\description{
Return the date and time of the last added sequence as
determined using the 'add_log.tsv'.
}
\details{
If no file found, returns empty character vector.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-setup-tools.R
\name{gb_sql_add}
\alias{gb_sql_add}
\title{Add to GenBank SQL database}
\usage{
gb_sql_add(df)
}
\arguments{
\item{df}{Records data.frame}
}
\description{
Add records data.frame to SQL-like database.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{restez_disconnect}
\alias{restez_disconnect}
\title{Disconnect from restez database}
\usage{
restez_disconnect()
}
\description{
Safely disconnect from the restez connection
}
\examples{
library(restez)
fp <- tempdir()
restez_path_set(filepath = fp)
demo_db_create(n = 5)
restez_connect()
restez_status()
# always remember to disconnect from a database when you've finished
restez_disconnect()
db_delete(everything = TRUE)
# Errors
# restez_status()
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_set}},
  \code{\link{restez_path_unset}},
  \code{\link{restez_ready}}, \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-setup-tools.R
\name{db_delete}
\alias{db_delete}
\title{Delete database}
\usage{
db_delete(everything = FALSE)
}
\arguments{
\item{everything}{T/F, delete the whole restez folder as well?}
}
\description{
Delete the local SQL database and/or restez folder.
}
\details{
Any connected database will be automatically disconnected.
}
\examples{
library(restez)
fp <- tempdir()
restez_path_set(filepath = fp)
demo_db_create(n = 10)
db_delete(everything = FALSE)
# Will not run: gb_sequence_get(id = 'demo_1')
# only the SQL database is deleted
db_delete(everything = TRUE)
# Now returns NULL
(restez_path_get())
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_create}}, \code{\link{db_download}},
  \code{\link{demo_db_create}}, \code{\link{is_in_db}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-setup-tools.R
\name{flatfile_read}
\alias{flatfile_read}
\title{Read flatfile sequence records}
\usage{
flatfile_read(flpth)
}
\arguments{
\item{flpth}{Path to .seq file}
}
\value{
list of GenBank records in text format
}
\description{
Read records from a .seq file.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{gb_build2}}, \code{\link{gb_build}},
  \code{\link{gb_df_create}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_definition}
\alias{extract_definition}
\title{Extract definition}
\usage{
extract_definition(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return definition from GenBank record.
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{restez_rl}
\alias{restez_rl}
\title{Restez readline}
\usage{
restez_rl(prompt)
}
\arguments{
\item{prompt}{character, display text}
}
\value{
character
}
\description{
Wrapper for base readline.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-setup-tools.R
\name{gb_build}
\alias{gb_build}
\title{Read and add .seq files to database}
\usage{
gb_build(dpth, seq_files, max_length, min_length)
}
\arguments{
\item{dpth}{Download path (where seq_files are stored)}

\item{seq_files}{.seq.tar seq file names}

\item{max_length}{Maximum sequence length}

\item{min_length}{Minimum sequence length}
}
\value{
Logical
}
\description{
Given a list of seq_files, read and add the contents of the
files to a SQL-like database. If any errors during the process, FALSE is
returned.
}
\details{
This function will automatically connect to the restez database.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_df_create}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-setup-tools.R
\name{db_download}
\alias{db_download}
\title{Download database}
\usage{
db_download(db = "nucleotide", overwrite = FALSE,
  preselection = NULL)
}
\arguments{
\item{db}{Database type, only 'nucleotide' currently available.}

\item{overwrite}{T/F, overwrite pre-existing downloaded files?}

\item{preselection}{character of user input}
}
\value{
T/F, if all files download correctly, TRUE else FALSE.
}
\description{
Download .seq.tar files from the latest GenBank release. The
user interactively selects the parts of GenBank to download (e.g. primates,
plants, bacteria ...)
}
\details{
The downloaded files will appear in the restez filepath under downloads.
}
\examples{
\dontrun{
library(restez)
restez_path_set(filepath = 'path/for/downloads')
db_download()
}
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_create}}, \code{\link{db_delete}},
  \code{\link{demo_db_create}}, \code{\link{is_in_db}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_definition_get}
\alias{gb_definition_get}
\title{Get definition from GenBank}
\usage{
gb_definition_get(id)
}
\arguments{
\item{id}{character, sequence accession ID(s)}
}
\value{
named vector of definitions, if no results found NULL
}
\description{
Return the definition line
for an accession ID.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(def <- gb_definition_get(id = 'demo_1'))
(defs <- gb_definition_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_fasta_get}},
  \code{\link{gb_organism_get}},
  \code{\link{gb_record_get}},
  \code{\link{gb_sequence_get}},
  \code{\link{gb_version_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{seshinfo_log}
\alias{seshinfo_log}
\title{Log the system session information in restez path}
\usage{
seshinfo_log()
}
\description{
Records the session and system information to file.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{setup}}, \code{\link{slctn_get}},
  \code{\link{slctn_log}}, \code{\link{sql_path_get}},
  \code{\link{status_class}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_organism_get}
\alias{gb_organism_get}
\title{Get organism from GenBank}
\usage{
gb_organism_get(id)
}
\arguments{
\item{id}{character, sequence accession ID(s)}
}
\value{
named vector of definitions, if no results found NULL
}
\description{
Return the organism name
for an accession ID.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(org <- gb_organism_get(id = 'demo_1'))
(orgs <- gb_organism_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_definition_get}},
  \code{\link{gb_fasta_get}}, \code{\link{gb_record_get}},
  \code{\link{gb_sequence_get}},
  \code{\link{gb_version_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{last_dwnld_get}
\alias{last_dwnld_get}
\title{Return date and time of the last download}
\usage{
last_dwnld_get()
}
\value{
character
}
\description{
Return the date and time of the last download as determined
using the 'download_log.tsv'.
}
\details{
If no file found, returns empty character vector.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez-tools.R
\name{entrez_gb_get}
\alias{entrez_gb_get}
\title{Get Entrez GenBank record}
\usage{
entrez_gb_get(id, ...)
}
\arguments{
\item{id}{vector, unique ID(s) for record(s)}

\item{...}{arguments passed on to rentrez}
}
\value{
character string containing the file created
}
\description{
Return gb and gbwithparts format as expected from
an Entrez call. If not all IDs are returned, will
run rentrez::entrez_fetch.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-tools.R
\name{mock_org}
\alias{mock_org}
\title{Mock org}
\usage{
mock_org(i)
}
\arguments{
\item{i}{integer, iterator}
}
\value{
character
}
\description{
Make a mock sequence organism.
Designed to be part of a loop.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_rec}}, \code{\link{mock_seq}},
  \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{connection_get}
\alias{connection_get}
\title{Retrieve restez connection}
\usage{
connection_get()
}
\value{
connection
}
\description{
Safely acquire the restez connection. Raises error if no
connection set.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomartr-tools.R
\name{check_connection}
\alias{check_connection}
\title{Helper function to test if a stable internet connection
can be established.}
\usage{
check_connection()
}
\description{
All retrieval functions need a stable
internet connection to work properly. This internal function pings
the google homepage and throws an error if it cannot be reached.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{cleanup}}, \code{\link{connected}},
  \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\author{
Hajk-Georg Drost
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-setup-tools.R
\name{demo_db_create}
\alias{demo_db_create}
\title{Create demo database}
\usage{
demo_db_create(db_type = "nucleotide", n = 100)
}
\arguments{
\item{db_type}{character, database type}

\item{n}{integer, number of mock sequences}
}
\description{
Creates a local mock SQL database
from package test data for demonstration purposes.
No internet connection required.
}
\examples{
library(restez)
# set the restez path to a temporary dir
restez_path_set(filepath = tempdir())
# create demo database
demo_db_create(n = 5)
restez_connect()
# in the demo, IDs are 'demo_1', 'demo_2' ...
(gb_sequence_get(id = 'demo_1'))

# Delete a demo database after an example
db_delete(everything = TRUE)
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_create}}, \code{\link{db_delete}},
  \code{\link{db_download}}, \code{\link{is_in_db}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{has_data}
\alias{has_data}
\title{Does the connected database have data?}
\usage{
has_data()
}
\value{
Logical
}
\description{
Returns TRUE if a restez SQL database has data.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_keywords}
\alias{extract_keywords}
\title{Extract keywords}
\usage{
extract_keywords(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character vector
}
\description{
Return keywords as list from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-tools.R
\name{stat}
\alias{stat}
\title{Print blue}
\usage{
stat(...)
}
\arguments{
\item{...}{Any number of text arguments to print, character}
}
\value{
coloured character encoding, character
}
\description{
Print to console blue text to indicate a number/statistic.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/status-tools.R
\name{status_class}
\alias{status_class}
\title{Generate a list class for storing status information}
\usage{
status_class()
}
\value{
Status class
}
\description{
Creates a three-part list for holding information on the
status of the restez file path.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
\name{record}
\alias{record}
\docType{data}
\title{
Example GenBank record
}
\description{
Example GenBank record in text format for demonstration purposes.
}
\usage{data("record")}
\format{
  A large character object containing record information and DNA sequence.
}
\source{
\url{https://www.ncbi.nlm.nih.gov/nuccore/AY952423.1}
}
\references{
GenBank
}
\examples{
data(record)
cat(record)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez-tools.R
\name{entrez_fasta_get}
\alias{entrez_fasta_get}
\title{Get Entrez fasta}
\usage{
entrez_fasta_get(id, ...)
}
\arguments{
\item{id}{vector, unique ID(s) for record(s)}

\item{...}{arguments passed on to rentrez}
}
\value{
character string containing the file created
}
\description{
Return fasta format as expected from
an Entrez call. If not all IDs are returned, will
run rentrez::entrez_fetch.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-tools.R
\name{mock_rec}
\alias{mock_rec}
\title{Mock rec}
\usage{
mock_rec(i, definition = NULL, accession = NULL, version = NULL,
  organism = NULL, sequence = NULL)
}
\arguments{
\item{i}{integer, iterator}

\item{definition}{character}

\item{accession}{character}

\item{version}{character}

\item{organism}{character}

\item{sequence}{character}
}
\value{
character
}
\description{
Create a mock GenBank record for demo-ing and testing purposes.
Designed to be part of a loop. Accession, organism... etc. are optional
arguments.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_seq}},
  \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{gbrelease_check}
\alias{gbrelease_check}
\title{Check if the last GenBank release number is the latest}
\usage{
gbrelease_check()
}
\value{
logical
}
\description{
Returns TRUE if the GenBank release number is the most recent
GenBank release available.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}}, \code{\link{gbrelease_get}},
  \code{\link{gbrelease_log}}, \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-tools.R
\name{mock_gb_df_generate}
\alias{mock_gb_df_generate}
\title{Generate mock GenBank records data.frame}
\usage{
mock_gb_df_generate(n)
}
\arguments{
\item{n}{integer, number of entries}
}
\value{
data.frame
}
\description{
Make a mock nucleotide data.frame
for entry into a demonstration SQL database.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_inforecpart}
\alias{extract_inforecpart}
\title{Extract the information record part}
\usage{
extract_inforecpart(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return information part from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-tools.R
\name{latest_genbank_release}
\alias{latest_genbank_release}
\title{Retrieve latest GenBank release number}
\usage{
latest_genbank_release()
}
\value{
character
}
\description{
Downloads the latest GenBank release number and returns it.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_version}
\alias{extract_version}
\title{Extract version}
\usage{
extract_version(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return accession + version ID from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{slctn_get}
\alias{slctn_get}
\title{Retrieve GenBank selections made by user}
\usage{
slctn_get()
}
\value{
character vector
}
\description{
Returns the selections made by the user.
}
\details{
If no file found, returns empty character vector.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_log}}, \code{\link{sql_path_get}},
  \code{\link{status_class}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{dwnld_rcrd_log}
\alias{dwnld_rcrd_log}
\title{Log a downloaded file in the restez path}
\usage{
dwnld_rcrd_log(fl)
}
\arguments{
\item{fl}{file name, character}
}
\description{
This function is called whenever a file is successfully
downloaded. A row entry is added to the 'download_log.tsv' in the user's
restez path containing the file name, the GB release number and the time of
successfully download. The log is to help users keep track of when they
downloaded files and to determine if the downloaded files are out of date.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{readme_log}
\alias{readme_log}
\title{Create README in restez_path}
\usage{
readme_log()
}
\description{
Write notes for the curious sorts who peruse the restez_path.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_features}
\alias{extract_features}
\title{Extract features}
\usage{
extract_features(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
list of lists
}
\description{
Return feature table as list from GenBank record
}
\details{
If element is not found, empty list returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-tools.R
\name{setup}
\alias{setup}
\title{Set up test common test data}
\usage{
setup()
}
\description{
Creates temporary test folders.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{slctn_get}},
  \code{\link{slctn_log}}, \code{\link{sql_path_get}},
  \code{\link{status_class}}, \code{\link{stat}},
  \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_seqrecpart}
\alias{extract_seqrecpart}
\title{Extract the sequence record part}
\usage{
extract_seqrecpart(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return sequence part from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-tools.R
\name{mock_seq}
\alias{mock_seq}
\title{Mock seq}
\usage{
mock_seq(i, sqlngth = 10)
}
\arguments{
\item{i}{integer, iterator}

\item{sqlngth}{integer, sequence length}
}
\value{
character
}
\description{
Make a mock sequence. Designed to be part of a loop.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restez.R
\docType{package}
\name{restez}
\alias{restez}
\alias{restez-package}
\title{restez: Create and Query a Local Copy of GenBank in R}
\description{
The restez package comes with five families of functions:
setup, database, get, entrez and internal/private.
}
\section{Setup functions}{

These functions allow a user to set the filepath for where the GenBank files
should be stored, create connections and verify these settings.
}

\section{Database functions}{

These functions download specific parts of GenBank and create the local
SQL-like database.
}

\section{GenBank functions}{

These functions allow a user to query the local SQL-like database. A
user can use an NCBI accession ID to retrieve sequences or whole GenBank
records.
}

\section{Entrez functions}{

The entrez functions are wrappers to the \code{entrez_*} functions in the
rentrez package. e.g the restez's entrez_fetch will first try to search the
local database, if it fails it will then call rentrez's
\code{\link[rentrez]{entrez_fetch}} with the same arguments.
}

\section{Private/internal functions}{

These functions work behind the scenes to make everything work. If you're
curious you can read their documentation using the form
\code{?restez:::functionname}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{add_rcrd_log}
\alias{add_rcrd_log}
\title{Log files added to the SQL database in the restez path}
\usage{
add_rcrd_log(fl)
}
\arguments{
\item{fl}{filename, character}
}
\description{
This function is called whenever sequence files have been
successfully added to the nucleotide SQL database. Row entries are added to
'add_lot.tsv' in the user's restez path containing the filename, GB release
numbers and the time of successful adding.
The log is to help users keep track of when sequences have been added.
}
\seealso{
Other private: \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{dir_size}
\alias{dir_size}
\title{Calculate the size of a directory}
\usage{
dir_size(fp)
}
\arguments{
\item{fp}{File path, character}
}
\value{
numeric
}
\description{
Returns the size of directory in GB
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez-tools.R
\name{message_missing}
\alias{message_missing}
\title{Produce message of missing IDs}
\usage{
message_missing(n)
}
\arguments{
\item{n}{Number of missing IDs}
}
\description{
Sends message to console stating number of missing IDs.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_sequence}
\alias{extract_sequence}
\title{Extract sequence}
\usage{
extract_sequence(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return sequence from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-tools.R
\name{file_download}
\alias{file_download}
\title{Download a file}
\usage{
file_download(fl, overwrite = FALSE)
}
\arguments{
\item{fl}{character, base filename (e.g. gbpri9.seq) to be
downloaded}

\item{overwrite}{T/F}
}
\value{
T/F
}
\description{
Download a GenBank .seq.tar file. Check
the file has downloaded properly. If not, returns FALSE.
If overwrite is true, any previous file will be overwritten.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{filename_log}}, \code{\link{flatfile_read}},
  \code{\link{gb_build2}}, \code{\link{gb_build}},
  \code{\link{gb_df_create}}, \code{\link{gb_df_generate}},
  \code{\link{gb_sql_add}}, \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-tools.R
\name{extract_accession}
\alias{extract_accession}
\title{Extract accession}
\usage{
extract_accession(record)
}
\arguments{
\item{record}{GenBank record in text format, character}
}
\value{
character
}
\description{
Return accession ID from GenBank record
}
\details{
If element is not found, '' returned.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filepath-tools.R
\name{restez_path_unset}
\alias{restez_path_unset}
\title{Unset restez path}
\usage{
restez_path_unset()
}
\description{
Set the restez path to NULL
}
\details{
Any connected database will be automatically disconnected.
}
\seealso{
Other setup: \code{\link{restez_connect}},
  \code{\link{restez_disconnect}},
  \code{\link{restez_path_get}},
  \code{\link{restez_path_set}},
  \code{\link{restez_ready}}, \code{\link{restez_status}}
}
\concept{setup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection-tools.R
\name{quiet_connect}
\alias{quiet_connect}
\title{Quiely connect to the restez database}
\usage{
quiet_connect()
}
\description{
Quiet version of restez_connect for automatic connections.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-setup-tools.R
\name{db_create}
\alias{db_create}
\title{Create new NCBI database}
\usage{
db_create(db_type = "nucleotide", min_length = 0, max_length = NULL,
  alt_restez_path = NULL)
}
\arguments{
\item{db_type}{character, database type}

\item{min_length}{Minimum sequence length, default 0.}

\item{max_length}{Maximum sequence length, default NULL.}

\item{alt_restez_path}{Alternative restez path if you would like to use the
downloads from a different restez path.}
}
\description{
Create a new local SQL database from downloaded files.
Currently only GenBank/nucleotide/nuccore database is supported.
}
\details{
All .seq.gz files are added to the database. A user can specify sequence
limit sizes for those sequences to be added to the database -- smaller
databases are faster to search.

Alternatively, a user can use the \code{alt_restez_path} to add the files
from an alternative restez file path. For example, you may wish to have a
database of all environmental sequences but then an additional smaller one of
just the sequences with lengths below 100 bp. Instead of having to download
all environmental sequences twice, you can generate multiple restez databases
using the same downloaded files from a single restez path.

This function will not overwrite a pre-existing database. Old databases must
be deleted before a new one can be created. Use \code{\link{db_delete}} with
everything=FALSE to delete an SQL database.

Connections/disconnections to the database are made automatically.
}
\examples{
\dontrun{
library(restez)
restez_path_set(filepath = 'path/for/downloads/and/database')
db_download()
db_create()
}
}
\seealso{
Other database: \code{\link{count_db_ids}},
  \code{\link{db_delete}}, \code{\link{db_download}},
  \code{\link{demo_db_create}}, \code{\link{is_in_db}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-get-tools.R
\name{count_db_ids}
\alias{count_db_ids}
\title{Return the number of ids}
\usage{
count_db_ids(db = "nucleotide")
}
\arguments{
\item{db}{character, database name}
}
\value{
integer
}
\description{
Return the number of ids in a user's restez database.
}
\details{
Requires an open connection. If no connection or db 0 is returned.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(count_db_ids())

# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other database: \code{\link{db_create}},
  \code{\link{db_delete}}, \code{\link{db_download}},
  \code{\link{demo_db_create}}, \code{\link{is_in_db}},
  \code{\link{list_db_ids}}
}
\concept{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log-tools.R
\name{gbrelease_log}
\alias{gbrelease_log}
\title{Log the GenBank release number in the restez path}
\usage{
gbrelease_log(release)
}
\arguments{
\item{release}{GenBank release number, character}
}
\description{
This function is called whenever db_download is run. It logs the
GB release number in the 'gb_release.txt' in the user's restez path.
The log is to help users keep track of whether their database if out of date.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{cat_line}}, \code{\link{char}},
  \code{\link{check_connection}}, \code{\link{cleanup}},
  \code{\link{connected}}, \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gb-get-tools.R
\name{gb_record_get}
\alias{gb_record_get}
\title{Get record from GenBank}
\usage{
gb_record_get(id)
}
\arguments{
\item{id}{character, sequence accession ID(s)}
}
\value{
named vector of records, if no results found NULL
}
\description{
Return the entire GenBank record for an accession ID.
Equivalent to \code{rettype='gb'} in \code{\link[rentrez]{entrez_fetch}}.
}
\examples{
library(restez)
restez_path_set(filepath = tempdir())
demo_db_create(n = 5)
restez_connect()
(rec <- gb_record_get(id = 'demo_1'))
(recs <- gb_record_get(id = c('demo_1', 'demo_2')))


# delete demo after example
db_delete(everything = TRUE)
}
\seealso{
Other get: \code{\link{gb_definition_get}},
  \code{\link{gb_fasta_get}},
  \code{\link{gb_organism_get}},
  \code{\link{gb_sequence_get}},
  \code{\link{gb_version_get}}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-tools.R
\name{cat_line}
\alias{cat_line}
\title{Cat lines}
\usage{
cat_line(...)
}
\arguments{
\item{...}{Text to print, character}
}
\description{
Helper function for printing lines to console. Automatically
formats lines by adding newlines.
}
\seealso{
Other private: \code{\link{add_rcrd_log}},
  \code{\link{char}}, \code{\link{check_connection}},
  \code{\link{cleanup}}, \code{\link{connected}},
  \code{\link{connection_get}},
  \code{\link{custom_download2}},
  \code{\link{custom_download}},
  \code{\link{db_sqlngths_get}},
  \code{\link{db_sqlngths_log}}, \code{\link{dir_size}},
  \code{\link{dwnld_path_get}},
  \code{\link{dwnld_rcrd_log}},
  \code{\link{entrez_fasta_get}},
  \code{\link{entrez_gb_get}},
  \code{\link{extract_accession}},
  \code{\link{extract_by_patterns}},
  \code{\link{extract_clean_sequence}},
  \code{\link{extract_definition}},
  \code{\link{extract_features}},
  \code{\link{extract_inforecpart}},
  \code{\link{extract_keywords}},
  \code{\link{extract_locus}},
  \code{\link{extract_organism}},
  \code{\link{extract_seqrecpart}},
  \code{\link{extract_sequence}},
  \code{\link{extract_version}},
  \code{\link{file_download}}, \code{\link{filename_log}},
  \code{\link{flatfile_read}}, \code{\link{gb_build2}},
  \code{\link{gb_build}}, \code{\link{gb_df_create}},
  \code{\link{gb_df_generate}}, \code{\link{gb_sql_add}},
  \code{\link{gb_sql_query}},
  \code{\link{gbrelease_check}},
  \code{\link{gbrelease_get}}, \code{\link{gbrelease_log}},
  \code{\link{has_data}},
  \code{\link{identify_downloadable_files}},
  \code{\link{last_add_get}}, \code{\link{last_dwnld_get}},
  \code{\link{last_entry_get}},
  \code{\link{latest_genbank_release_notes}},
  \code{\link{latest_genbank_release}},
  \code{\link{message_missing}}, \code{\link{mock_def}},
  \code{\link{mock_gb_df_generate}},
  \code{\link{mock_org}}, \code{\link{mock_rec}},
  \code{\link{mock_seq}}, \code{\link{predict_datasizes}},
  \code{\link{print.status}}, \code{\link{quiet_connect}},
  \code{\link{readme_log}},
  \code{\link{restez_path_check}}, \code{\link{restez_rl}},
  \code{\link{seshinfo_log}}, \code{\link{setup}},
  \code{\link{slctn_get}}, \code{\link{slctn_log}},
  \code{\link{sql_path_get}}, \code{\link{status_class}},
  \code{\link{stat}}, \code{\link{testdatadir_get}}
}
\concept{private}
