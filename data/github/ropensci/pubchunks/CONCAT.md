

pubchunks
=========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/pubchunks)](https://cranchecks.info/pkgs/pubchunks)
[![R-check](https://github.com/ropensci/pubchunks/workflows/R-check/badge.svg)](https://github.com/ropensci/pubchunks/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/pubchunks/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/pubchunks)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/pubchunks)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/pubchunks)](https://cran.r-project.org/package=pubchunks)

## Get chunks of XML articles


## Package API

 - pub_tabularize
 - pub_guess_publisher
 - pub_sections
 - pub_chunks
 - pub_providers

The main workhorse function is `pub_chunks()`. It allows you to pull out sections of articles from many different publishers (see next section below) WITHOUT having to know how to parse/navigate XML. XML has a steep learning curve, and can require quite a bit of Googling to sort out how to get to different parts of an XML document. 

The other main function is `pub_tabularize()` - which takes the output of `pub_chunks()` and coerces into a data.frame for easier downstream processing.

## Supported publishers/sources

- eLife
- PLOS
- Entrez/Pubmed
- Elsevier
- Hindawi
- Pensoft
- PeerJ
- Copernicus
- Frontiers
- F1000 Research

If you know of other publishers or sources that provide XML let us know by [opening an issue](https://github.com/ropensci/pubchunks/issues).

We'll continue adding additional publishers.


## Installation

Stable version


```r
install.packages("pubchunks")
```

Development version from GitHub


```r
remotes::install_github("ropensci/pubchunks")
```

Load library


```r
library('pubchunks')
```

## Working with files


```r
x <- system.file("examples/10_1016_0021_8928_59_90156_x.xml", 
  package = "pubchunks")
```


```r
pub_chunks(x, "abstract")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: abstract
#>   showing up to first 5: 
#>    abstract (n=1): Abstract
#>                
#>                   This pa ...
pub_chunks(x, "title")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
pub_chunks(x, "authors")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: authors
#>   showing up to first 5: 
#>    authors (n=1): Chetaev, D.N
pub_chunks(x, c("title", "refs"))
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title, refs
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
#>    refs (n=6): Watson G.N.. 1949. Teoriia besselevykh funktsii. N
```

The output of `pub_chunks()` is a list with an S3 class `pub_chunks` to make 
internal work in the package easier. You can easily see the list structure 
by using `unclass()`.

## Working with the xml already in a string


```r
xml <- paste0(readLines(x), collapse = "")
pub_chunks(xml, "title")
#> <pub chunks>
#>   from: character
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
```

## Working with xml2 class object


```r
xml <- paste0(readLines(x), collapse = "")
xml <- xml2::read_xml(xml)
pub_chunks(xml, "title")
#> <pub chunks>
#>   from: xml_document
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
```

## Working with output of fulltext::ft_get()


```r
install.packages("fulltext")
```


```r
library("fulltext")
x <- fulltext::ft_get('10.1371/journal.pone.0086169')
pub_chunks(fulltext::ft_collect(x), sections="authors")
#> $plos
#> $plos$`10.1371/journal.pone.0086169`
#> <pub chunks>
#>   from: xml_document
#>   publisher/journal: plos/PLoS ONE
#>   sections: authors
#>   showing up to first 5: 
#>    authors (n=4): nested list
#> 
#> 
#> attr(,"ft_data")
#> [1] TRUE
```

## Coerce pub_chunks output into data.frame's


```r
x <- system.file("examples/elife_1.xml", package = "pubchunks")
res <- pub_chunks(x, c("doi", "title", "keywords"))
pub_tabularize(res)
#>                   doi                                          title
#> 1 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 2 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 3 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 4 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 5 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 6 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#>                       keywords .publisher
#> 1                     microRNA      elife
#> 2            nonsense mutation      elife
#> 3 nonsense-mediated mRNA decay      elife
#> 4                          APC      elife
#> 5             intron retention      elife
#> 6  premature termination codon      elife
```

## Get a random XML article


```r
library(rcrossref)
library(dplyr)

res <- cr_works(filter = list(
    full_text_type = "application/xml", 
    license_url="http://creativecommons.org/licenses/by/4.0/"))
links <- bind_rows(res$data$link) %>% filter(content.type == "application/xml")
download.file(links$URL[1], (i <- tempfile(fileext = ".xml")))
pub_chunks(i)
#> <pub chunks>
#>   from: file
#>   publisher/journal: unknown/NA
#>   sections: all
#>   showing up to first 5: 
#>    front (n=0): 
#>    body (n=0): 
#>    back (n=0): 
#>    title (n=0): 
#>    doi (n=0):
download.file(links$URL[13], (j <- tempfile(fileext = ".xml")))
pub_chunks(j)
#> <pub chunks>
#>   from: file
#>   publisher/journal: hindawi/BioMed Research International
#>   sections: all
#>   showing up to first 5: 
#>    front (n=2): nested list
#>    body (n=49): Oxidative stress and Reactive Oxygen Species (ROS)
#>    back (n=4): nested list
#>    title (n=1): Selected Enzyme Inhibitory Effects of Euphorbia ch ...
#>    doi (n=1): 10.1155/2018/1219367
download.file(links$URL[20], (k <- tempfile(fileext = ".xml")))
pub_chunks(k)
#> <pub chunks>
#>   from: file
#>   publisher/journal: hindawi/Case Reports in Pathology
#>   sections: all
#>   showing up to first 5: 
#>    front (n=2): nested list
#>    body (n=16): Bonnetti et al. first noted in 1992 an unusual cel
#>    back (n=3): nested list
#>    title (n=1): An Inguinal Perivascular Epithelioid Cell Tumor Me ...
#>    doi (n=1): 10.1155/2018/5749421
```




## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pubchunks/issues).
* License: MIT
* Get citation information for `pubchunks`: `citation(package = 'pubchunks')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
pubchunks 0.4.0
===============

### NEW FEATURES

* `pub_chunks()` gains new parameter `extract` to determine the final step of extracting text from XML format articles - use either `as.character()` or `xml2::xml_text()` (#11)

### MINOR IMPROVEMENTS

* fix to `sections="aff"` in `pub_chunks()` - keep all rows in merging data from different affilitation metadata parts to not lose any data (#7)


pubchunks 0.3.0
===============

### MINOR IMPROVEMENTS

* improvements to `pub_chunks()` with `sections="refs"` (fetching references from an article). all reference extraction was simply extracting the entire reference, which lost all white space. now there are custom reference parsers for the various types of reference formats. there still are outstanding issues, so do please point out where ref. extraction is not correct (#2)
* fix pub_guess_publisher examples error (#10)

pubchunks 0.2.2
===============

### MINOR IMPROVEMENTS

* fixed failing example (#9)

pubchunks 0.2.0
===============

### MINOR IMPROVEMENTS

* most section options in `pub_chunks()` now have defaults for extracting the section, and return NULL/empty list when not found (#3) (#4)
* improvements to `print.pub_chunks` so that the printed object contains more information (publisher/journal title) and more accurate ('character' used to include xml as character string and file paths, but are separated now). in addition, we state that the first 5 sections are printed so the user knows there could be more (#8)
* fix `pub_tabularize()` to accept list outputs from `pub_chunks()` (#5)

pubchunks 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There were no problems with the 1 reverse dependency.

--------

This version introduced a new parameter in a function and fixes an XML parsing issue.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/pubchunks/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/pubchunks.git`
* Make sure to track progress upstream (i.e., on our version of `pubchunks` at `ropensci/pubchunks`) by doing `git remote add upstream https://github.com/ropensci/pubchunks.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/pubchunks`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
