bold
====


[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran version](https://www.r-pkg.org/badges/version/bold)](https://cran.r-project.org/package=bold)
[![cran checks](https://cranchecks.info/badges/worst/bold)](https://cranchecks.info/pkgs/bold)
[![R-check](https://github.com/ropensci/bold/workflows/R-check/badge.svg)](https://github.com/ropensci/bold/actions/)
[![codecov.io](https://codecov.io/github/ropensci/bold/coverage.svg?branch=master)](https://codecov.io/github/ropensci/bold?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/bold)](https://github.com/r-hub/cranlogs.app)

`bold` accesses BOLD barcode data.

The Barcode of Life Data Systems (BOLD) is designed to support the generation and application of DNA barcode data. The platform consists of four main modules: a data portal, a database of barcode clusters, an educational portal, and a data collection workbench.

This package retrieves data from the BOLD database of barcode clusters, and allows for searching of over 1.7M public records using multiple search criteria including sequence data, specimen data, specimen *plus* sequence data, as well as trace files.

Documentation for the BOLD API: http://v4.boldsystems.org/index.php/api_home

See also the taxize book for more options for taxonomic workflows with BOLD: https://taxize.dev/

## Installation

__Installation instructions__

__Stable Version__


```r
install.packages("bold")
```

__Development Version__

Install `sangerseqR` first (used in function `bold::bold_trace()` only)

For R < 3.5


```r
source("http://bioconductor.org/biocLite.R")
biocLite("sangerseqR")
```

For R >= 3.5


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sangerseqR")
```

Then install `bold`


```r
remotes::install_github("ropensci/bold")
```


## Usage

```r
library("bold")
```


### Search for sequence data only

Default is to get a list back


```r
bold_seq(taxon='Coelioxys')[[1]]
#> $id
#> [1] "ABEE117-17"
#> 
#> $name
#> [1] "Coelioxys elongata"
#> 
#> $gene
#> [1] "ABEE117-17"
#> 
#> $sequence
#> [1] "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTATCATTATATACATATCATCCTTCCCCATCAGTTGATTTAGCAATTTTTTYTTTACATTTATCAGGAATTTYTTYTATTATCGGATCAATAAATTTTATTGTAACAATTTTAATAATAAAAAATTATTCAATAAATTATAATCAAATACCTTTATTTCCATGATCAATTTTAATTACTACAATTTTATTATTATTATCATTACCTGTATTAGCAGGAGCTATTACAATATTATTATTTGATCGTAATTTAAATTCATCATTTTTTGACCCAATAGGAGGAGGAGATCCTATTTTATATCAACATTTATTTTG------------------------------------"
```

You can optionally get back the `crul` response object


```r
res <- bold_seq(taxon='Coelioxys', response=TRUE)
res$response_headers
#> $status
#> [1] "HTTP/2 200 "
#> 
#> $server
#> [1] "nginx"
#> 
#> $date
#> [1] "Mon, 20 Apr 2020 16:11:50 GMT"
#> 
#> $`content-type`
#> [1] "application/x-download"
#> 
#> $`x-powered-by`
#> [1] "PHP/5.3.15"
#> 
#> $`content-disposition`
#> [1] "attachment; filename=fasta.fas"
#> 
#> $`x-frame-options`
#> [1] "SAMEORIGIN"
#> 
#> $`x-content-type-options`
#> [1] "nosniff"
#> 
#> $`x-xss-protection`
#> [1] "1; mode=block"
```

### Search for specimen data only

By default you download `tsv` format data, which is given back to you as a `data.frame`


```r
res <- bold_specimens(taxon='Osmia')
head(res[,1:8])
#>      processid   sampleid recordID catalognum   fieldnum
#> 1  BEECA373-06 05-NT-0373   514740            05-NT-0373
#> 2  BEECA607-06 05-NT-0607   516959            05-NT-0607
#> 3  BEECA963-07 01-OR-0790   554153            01-OR-0790
#> 4  BEECB358-07 04-WA-1076   596920 BBSL697174 04-WA-1076
#> 5  BEECB438-07 00-UT-1157   597000 BBSL432653 00-UT-1157
#> 6 BEECC1176-09 02-UT-2849  1060879 BBSL442586 02-UT-2849
#>                    institution_storing collection_code      bin_uri
#> 1   York University, Packer Collection              NA BOLD:AAI2013
#> 2   York University, Packer Collection              NA BOLD:AAC8510
#> 3   York University, Packer Collection              NA BOLD:ABZ3184
#> 4 Utah State University, Logan Bee Lab              NA BOLD:AAC5797
#> 5 Utah State University, Logan Bee Lab              NA BOLD:AAF2159
#> 6   York University, Packer Collection              NA BOLD:AAE5368
```

### Search for specimen plus sequence data

By default you download `tsv` format data, which is given back to you as a `data.frame`


```r
res <- bold_seqspec(taxon='Osmia', sepfasta=TRUE)
res$fasta[1:2]
#> $`BEECA373-06`
#> [1] "-ATTTTATATATAATTTTTGCTATATGATCAGGTATAATCGGATCAGCAATAAGAATTATTATTCGTATAGAATTAAGAATTCCTGGTTCATGAATTTCAAATGATCAAACTTATAACTCTTTAGTAACTGCTCATGCTTTTTTAATAATTTTTTTCTTAGTTATACCTTTTTTAATTGGAGGATTTGGAAATTGATTAATTCCTTTAATATTAGGAATCCCGGATATAGCTTTCCCTCGAATAAATAATATTAGATTTTGACTTTTACCCCCTTCATTAATATTATTACTTTTAAGAAATTTTATAAATCCAAGACCAGGTACTGGATGAACTGTTTATCCTCCTCTTTCTTCTCATTTATTTCATTCTTCTCCTTCAGTTGATATAGCCATTTTTTCTTTACATATTTCCGGTTTATCTTCTATTATAGGTTCGTTAAATTTTATTGTTACAATTATTATAATAAAAAATATTTCTTTAAAACATATCCAATTACCTTTATTTCCATGATCTGTTTTTATTACTACTATCTTATTACTTTTTTCTTTACCTGTTTTAGCAGGAGCTATTACTATATTATTATTTGATCGAAATTTTAATACTTCATTTTTTGATCCTACAGGAGGTGGAGATCCAATCCTTTATCAACATTTATTT"
#> 
#> $`BEECA607-06`
#> [1] "AATATTATATATAATTTTTGCTTTGTGATCTGGAATAATTGGTTCATCTATAAGAATTATTATTCGTATAGAATTAAGAATTCCTGGTTCATGAATTTCAAATGATCAAGTTTATAATTCATTAGTTACAGCTCATGCTTTTTTAATAATTTTTTTTTTAGTTATACCATTTATAATTGGAGGATTTGGAAATTGATTAGTTCCTTTAATATTAGGAATTCCTGATATAGCTTTTCCTCGAATAAATAATATTAGATTTTGATTATTACCACCATCATTAATACTTTTACTTTTAAGAAATTTTTTTAATCCAAGTTCAGGAACTGGATGAACTGTATATCCTCCTCTTTCATCATATTTATTTCATTCTTCACCTTCTGTTGATTTAGCTATTTTTTCTCTTCATATATCAGGTTTATCTTCTATTATAGGTTCATTAAACTTTATTGTAACTATTATTATAATAAAAAATATTTCTTTAAAGTATATTCAATTGCCATTATTTCCATGATCTGTTTTTATTACTACAATTCTTTTATTATTATCATTACCAGTTTTAGCAGGTGCTATTACTATATTATTATTTGATCGAAATTTTAATACTTCATTTTTTGATCCTACAGGAGGGGGAG--------------------------"
```

Or you can index to a specific sequence like


```r
res$fasta['GBAH0293-06']
#> $`GBAH0293-06`
#> [1] "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAATGTTAGGGATTCCAGATATAGCTTTTCCACGAATAAATAATATTAGATTTTGACTGTTACCTCCATCTTTAATATTATTACTTTTAAGAAATTTTTTAAATCCAAGTCCTGGAACAGGATGAACAGTTTATCCTCCTTTATCATCAAATTTATTTCATTCTTCTCCTTCAGTTGATTTAGCAATTTTTTCTTTACATATTTCAGGTTTATCTTCTATTATAGGTTCATTAAATTTTATTGTTACAATTATTATAATAAAAAATATTTCTTTAAAATATATTCAATTACCTTTATTTTCTTGATCTGTATTTATTACTACTATTCTTTTATTATTTTCTTTACCTGTATTAGCTGGAGCTATTACTATATTATTATTTGATCGAAATTTTAATACATCTTTTTTTGATCCAACAGGAGGGGGAGATCCAATTCTTTATCAACATTTATTTTGATTTTTTGGTCATCCTGAAGTTTATATTTTAATTTTACCTGGATTTGGATTAATTTCTCAAATTATTTCTAATGAAAGAGGAAAAAAAGAAACTTTTGGAAATATTGGTATAATTTATGCTATATTAAGAATTGGACTTTTAGGTTTTATTGTT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
```

### Get trace files

This function downloads files to your machine - it does not load them into your R session - but prints out where the files are for your information.


```r
x <- bold_trace(ids = 'ACRJP618-11', progress = FALSE)
read_trace(x$ab1)
#> Number of datapoints: 8877
#> Number of basecalls: 685
#> 
#> Primary Basecalls: NNNNNNNNNNNNNNNNNNGNNNTTGAGCAGGNATAGTAGGANCTTCTCTTAGTCTTATTATTCGAACAGAATTAGGAAATCCAGGATTTTTAATTGGAGATGATCAAATCTACAATACTATTGTTACGGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAATTGATTAGTTCCCCTTATACTAGGAGCCCCAGATATAGCTTTCCCTCGAATAAACAATATAAGTTTTTGGCTTCTTCCCCCTTCACTATTACTTTTAATTTCCAGAAGAATTGTTGAAAATGGAGCTGGAACTGGATGAACAGTTTATCCCCCACTGTCATCTAATATTGCCCATAGAGGTACATCAGTAGATTTAGCTATTTTTTCTTTACATTTAGCAGGTATTTCCTCTATTTTAGGAGCGATTAATTTTATTACTACAATTATTAATATACGAATTAACAGTATAAATTATGATCAAATACCACTATTTGTGTGATCAGTAGGAATTACTGCTTTACTCTTATTACTTTCTCTTCCAGTATTAGCAGGTGCTATCACTATATTATTAACGGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGGAGATCCAATTTTATATCAACATTTATTTTGATTTTTTGGACNTCNNNNAAGTTTAAN
#> 
#> Secondary Basecalls:
```

### Large data

Sometimes with `bold_seq()` you request a lot of data, which can cause problems due 
to BOLD's servers. 

An example is the taxonomic name _Arthropoda_. When you send a request like 
`bold_seq(taxon = "Arthropoda")` BOLD attempts to give you back sequences
for all records under _Arthropoda_. This, as you can imagine, is a lot of 
sequences. 



```r
library("taxize")
```

Using `taxize::downstream` get children of _Arthropoda_


```r
x <- downstream("Arthropoda", db = "ncbi", downto = "class")
#> ══  1 queries  ═══════════════
#> ✔  Found:  Arthropoda
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
nms <- x$Arthropoda$childtaxa_name
```

Optionally, check that the name exists in BOLD's data. Any that are not in 
BOLD will give back a row of NAs


```r
checks <- bold_tax_name(nms)
# all is good
checks[,1:5]
#>     taxid         taxon tax_rank tax_division parentid
#> 1   26059   Pycnogonida    class     Animalia       20
#> 2      63     Arachnida    class     Animalia       20
#> 3      74   Merostomata    class     Animalia       20
#> 4  493944     Pauropoda    class     Animalia       20
#> 5   80390      Symphyla    class     Animalia       20
#> 6      85     Diplopoda    class     Animalia       20
#> 7      75     Chilopoda    class     Animalia       20
#> 8      82       Insecta    class     Animalia       20
#> 9     372    Collembola    class     Animalia       20
#> 10 734357       Protura    class     Animalia       20
#> 11     84     Remipedia    class     Animalia       20
#> 12     73 Cephalocarida    class     Animalia       20
#> 13     68  Branchiopoda    class     Animalia       20
#> 14 765970   Hexanauplia    class     Animalia       20
#> 15     69  Malacostraca    class     Animalia       20
#> 16 889450 Ichthyostraca    class     Animalia       20
#> 17     NA          <NA>     <NA>         <NA>       NA
#> 18     80     Ostracoda    class     Animalia       20
```

Then pass those names to `bold_seq()`. You could pass all names in at once,
but we're trying to avoid the large data request problem here, so run each 
one separately with `lapply` or a for loop like request. 


```r
out <- lapply(nms, bold_seq)
```

## Citation

Get citation information for `bold` in R by running: `citation(package = 'bold')`

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/bold/issues)
* License: MIT
* Get citation information for `bold` in R doing `citation(package = 'bold')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
bold 1.2.0
==========

### MINOR IMPROVEMENTS

* vignettes fix (#77)

bold 1.1.0
==========

### MINOR IMPROVEMENTS

* fix a failing test (#73)


bold 1.0.0
==========

### MINOR IMPROVEMENTS

* change base url for all requests to https from http (#70)
* fixed a warning arising from use of `bold_seqspec()` - we now set the encoding to "UTF-8" before parsing the string to XML (#71)
* `bold_seqspec()` fix: capture "Fatal errors" returned by BOLD servers and pass that along to the user with advice (#66)
* add "Marker" and "Large requests" documentation sections to both `bold_seq()` and `bold_seqspec()`. the marker section details that the marker parameter doesn't actually filter results that you get - but you can filter them yourself. the large requests section gives some caveats associated with large data requests and outlines how to sort it out (#61)


bold 0.9.0
==========

### MINOR IMPROVEMENTS

* improved test coverage (#58)
* allow curl options to be passed into `bold_identify_parents()` (#64)
* fix instructions in README for package `sangerseqR` - instructions depend on which version of R is being used (#65) thanks @KevCaz

### BUG FIXES

* fixes in package for `_R_CHECK_LENGTH_1_LOGIC2_` (#57)
* `bold_identify()` fix: ampersands needed to be escaped (#62) thanks @devonorourke


bold 0.8.6
==========

### MINOR IMPROVEMENTS

* tests that make HTTP requests now use package `vcr` to cache responses, speeds up tests significantly, and no longer relies on an internet connection (#55) (#56)
* `bold_seq()`: sometimes on large requests, the BOLD servers time out, and give back partial output but don't indicate that there was an error. We catch this kind of error now, throw a message for the user, and the function gives back the partial output given by the server. Also added to the documentation for `bold_seq()` and in the README that if you run into this problem try to do many queries that will result in smaller set of results instead of one or fewer larger queries (#52) (#53)
* `bold_seq()`: remove return characters (`\r` and `\n`) from sequences (#54)


bold 0.8.0
==========

### MINOR IMPROVEMENTS

* link to taxize bookdown book in readme and vignette (#51)
* `bold_identify_parents()` gains many new parameters (`taxid`, `taxon`, `tax_rank`, `tax_division`, `parentid`, `parentname`, `taxonrep`, `specimenrecords`) to filter parents based on any of a number of fields - should solve problem where multiple parents found for a single taxon, often in different kingdoms (#50)
* add note in docs of `bold_identify()` that the function uses `lapply` internally, so queries with lots of sequences can take a long time

### BUG FIXES

* fix `bold_specimens()`: use `rawToChar()` on raw bytes instead of `parse()` from `crul` (#47)

bold 0.5.0
==========

### NEW FEATURES

* Now using BOLD's v4 API throughout the package. This was essentially
just a change of the BASE URL for each request (#30) 
* Now using `crul` for HTTP requests. Only really affects users in that
specifying curl options works slightly differenlty (#42) 

### BUG FIXES 

* `marker` parameter in `bold_seqspec` was and maybe still is not working, 
in the sense that using the parameter doesn't always limit results to the 
marker you specify. Not really fixed - watch out for it, and filter after you
get results back to get markers you want. (#25) 
* Fixed bug in `bold_identify_parents` - was failing when no match for a
parent name. (#41) thx @VascoElbrecht  
* `tsv` results were erroring in `bold_specimens` and other fxns (#46) - fixed
by switching to new BOLD v4 API (#30)

### MINOR IMPROVEMENTS

* Namespace calls to base pkgs for `stats` and `utils` - replaced 
`is` with `inherits` (#39) 



bold 0.4.0
==========

### NEW FEATURES

* New function `bold_identify_parents()` to add taxonomic information
to the output of `bold_identif()`. We take the taxon names from `bold_identify`
output, and use `bold_tax_name` to get the taxonomic ID, passing it to 
`bold_tax_id` to get the parent names, then attaches those to the input data. 
There are two options given what you put for the `wide` parameter. If `TRUE`
you get data.frames of the same dimensions with parent rank name and ID 
as new columns (for each name going up the hierarchy) - while if `FALSE` 
you get a long data.frame. thanks @dougwyu for inspiring this (#36) 

### MINOR IMPROVEMENTS

* replace `xml2::xml_find_one` with `xml2::xml_find_first` (#33)
* Fix description of `db` options in `bold_identify` man file - 
COX1 and COX1_SPECIES were switched (#37) thanks for pointing that out 
@dougwyu

### BUG FIXES

* Fix to `bold_tax_id` for when some elements returned from the BOLD 
API were empty/`NULL` (#32) thanks @fmichonneau !!


bold 0.3.5
==========

### MINOR IMPROVEMENTS

* Added more tests to the test suite (#28)

### BUG FIXES

* Fixed a bug in an internal data parser (#27)

bold 0.3.4
==========

### NEW FEATURES

* Added a code of conduct

### MINOR IMPROVEMENTS

* Switched to `xml2` from `XML` as the XML parser for this package (#26)
* Fixes to `bold_trace()` to create dir and tar file when it doesn't
already exist

### BUG FIXES

* Fixed odd problem where sometimes resulting data from HTTP request
was garbled on `content(x, "text")`, so now using `rawToChar(content(x))`,
which works (#24)


bold 0.3.0
==========

### MINOR IMPROVEMENTS

* Explicitly import non-base R functions (#22)
* Better package level manual file


bold 0.2.6
==========

### MINOR IMPROVEMENTS

* `sangerseqR` package now in Suggests for reading trace files, and is only used in `bold_trace()`
function.
* General code tidying, reduction of code duplication.
* `bold_trace()` gains two new parameters: `overwrite` to choose whether to overwrite an existing
file of the same name or not, `progress` to show a progress bar for downloading or not.
* `bold_trace()` gains a print method to show a tidy summary of the trace file downloaded.

### BUG FIXES

* Fixed similar bugs in `bold_tax_name()` (#17) and `bold_tax_id()` (#18) in which species that were missing from the BOLD database returned empty arrays but 200 status codes. Parsing those as failed attempts now. Also fixes problem in taxize in `bold_search()` that use these two functions.

bold 0.2.0
==========

### NEW FEATURES

* Package gains two new functions for working with the BOLD taxonomy APIs: `bold_tax_name()` and `bold_tax_id()`, which search for taxonomic data from BOLD using either names or BOLD identifiers, respectively. (#11)
* Two new packages in Imports: `jsonlite` and `reshape`.

### MINOR IMPROVEMENTS

* Added new taxonomy API functions to the vignette (#14)
* Added reference URLS to all function doc files to allow easy reference for the appropriate API docs.
* `callopts` parameter changed to `...` throughout the package, so that passing on options to `httr::GET` is done via named parameters, e.g., `config=verbose()`. (#13)
* Added examples of doing curl debugging throughout man pages.


bold 0.1.2
==========

### MINOR IMPROVEMENTS

* Improved the vignette (#8)
* Added small function to print helpful message when user inputs no parameters or zero length parameter values.

### BUG FIXES

* Fixed some broken tests with the new `httr` (v0.4) (#9), and added a few more tests (#7)


bold 0.1.0
==========

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.0.5 patched
* win-builder (devel and release)
* GitHub Actions (linux, macos, windows)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 downstream dependency.
(Summary at <https://github.com/ropensci/bold/blob/master/revdep/README.md>).
No problems were found.

-----

This version fixes the rmarkdown/markdown dependency issue for vignettes that Kurt emailed maintainers about.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/bold/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/bold.git`
* Make sure to track progress upstream (i.e., on our version of `bold` at `ropensci/bold`) by doing `git remote add upstream https://github.com/ropensci/bold.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/bold`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                  |
|:--------|:--------------------------------------|
|version  |R version 4.0.2 RC (2020-06-14 r78689) |
|os       |macOS Catalina 10.15.5                 |
|system   |x86_64, darwin17.0                     |
|ui       |X11                                    |
|language |(EN)                                   |
|collate  |en_US.UTF-8                            |
|ctype    |en_US.UTF-8                            |
|tz       |US/Pacific                             |
|date     |2020-06-17                             |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|bold    |1.0.0 |1.1.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*bold
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran version](https://www.r-pkg.org/badges/version/bold)](https://cran.r-project.org/package=bold)
[![cran checks](https://cranchecks.info/badges/worst/bold)](https://cranchecks.info/pkgs/bold)
[![R-check](https://github.com/ropensci/bold/workflows/R-check/badge.svg)](https://github.com/ropensci/bold/actions/)
[![codecov.io](https://codecov.io/github/ropensci/bold/coverage.svg?branch=master)](https://codecov.io/github/ropensci/bold?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/bold)](https://github.com/r-hub/cranlogs.app)

`bold` accesses BOLD barcode data.

The Barcode of Life Data Systems (BOLD) is designed to support the generation and application of DNA barcode data. The platform consists of four main modules: a data portal, a database of barcode clusters, an educational portal, and a data collection workbench.

This package retrieves data from the BOLD database of barcode clusters, and allows for searching of over 1.7M public records using multiple search criteria including sequence data, specimen data, specimen *plus* sequence data, as well as trace files.

Documentation for the BOLD API: http://v4.boldsystems.org/index.php/api_home

See also the taxize book for more options for taxonomic workflows with BOLD: https://taxize.dev/

## Installation

__Installation instructions__

__Stable Version__

```{r eval=FALSE}
install.packages("bold")
```

__Development Version__

Install `sangerseqR` first (used in function `bold::bold_trace()` only)

For R < 3.5

```{r eval=FALSE}
source("http://bioconductor.org/biocLite.R")
biocLite("sangerseqR")
```

For R >= 3.5

```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sangerseqR")
```

Then install `bold`

```{r eval=FALSE}
remotes::install_github("ropensci/bold")
```


## Usage
```{r}
library("bold")
```


### Search for sequence data only

By default you download `fasta` file, which is given back to you as a `data.frame`

```{r}
res <- bold_seq(taxon='Coelioxys')
head(res)
```

You can optionally get back the `crul` response object

```{r}
res <- bold_seq(taxon='Coelioxys', response=TRUE)
res$response_headers
```

### Search for specimen data only

By default you download `tsv` format data, which is given back to you as a `data.frame`

```{r}
res <- bold_specimens(taxon='Osmia')
head(res[,1:8])
```

### Search for specimen plus sequence data

By default you download `tsv` format data, which is given back to you as a `data.frame`

```{r}
res <- bold_seqspec(taxon='Osmia', sepfasta=TRUE)
res$fasta[1:2]
```

Or you can index to a specific sequence like

```{r}
res$fasta['GBAH0293-06']
```

### Get trace files

This function downloads files to your machine - it does not load them into your R session - but prints out where the files are for your information.

```{r}
x <- bold_trace(ids = 'ACRJP618-11', progress = FALSE)
read_trace(x$ab1)
```

### Large data

Sometimes with `bold_seq()` you request a lot of data, which can cause problems due 
to BOLD's servers. 

An example is the taxonomic name _Arthropoda_. When you send a request like 
`bold_seq(taxon = "Arthropoda")` BOLD attempts to give you back sequences
for all records under _Arthropoda_. This, as you can imagine, is a lot of 
sequences. 


```{r}
library("taxize")
```

Using `taxize::downstream` get children of _Arthropoda_

```{r}
x <- downstream("Arthropoda", db = "ncbi", downto = "class")
nms <- x$Arthropoda$childtaxa_name
```

Optionally, check that the name exists in BOLD's data. Any that are not in 
BOLD will give back a row of NAs

```{r}
checks <- bold_tax_name(nms)
# all is good
checks[,1:5]
```

Then pass those names to `bold_seq()`. You could pass all names in at once,
but we're trying to avoid the large data request problem here, so run each 
one separately with `lapply` or a for loop like request. 

```{r eval = FALSE}
out <- lapply(nms, bold_seq)
```

## Citation

Get citation information for `bold` in R by running: `citation(package = 'bold')`

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/bold/issues)
* License: MIT
* Get citation information for `bold` in R doing `citation(package = 'bold')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
---
title: "bold introduction"
author: "Scott Chamberlain"
date: "2020-04-20"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bold introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



`bold` is an R package to connect to BOLD Systems (https://www.boldsystems.org/) via their API. Functions in `bold` let you search for sequence data, specimen data, sequence + specimen data, and download raw trace files.

### bold info

+ BOLD home page: https://boldsystems.org/
+ BOLD API docs: https://v4.boldsystems.org/index.php/api_home

See also the taxize book for more options for taxonomic workflows with BOLD: https://taxize.dev/

### Using bold

**Install**

Install `bold` from CRAN



```r
install.packages("bold")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/bold")
```

Load the package


```r
library("bold")
```


### Search for taxonomic names via names

`bold_tax_name` searches for names with names.


```r
bold_tax_name(name = 'Diplura')
#>    taxid   taxon tax_rank tax_division parentid       parentname
#> 1 603673 Diplura    genus     Protista    53974 Scytosiphonaceae
#> 2 734358 Diplura    class     Animalia       20       Arthropoda
#>   specimenrecords taxonrep     representitive_image.image
#> 1               6     <NA>                           <NA>
#> 2             308  Diplura BSOIL/GBOL20120+1540220834.jpg
#>   representitive_image.apectratio   input
#> 1                              NA Diplura
#> 2                           0.841 Diplura
```


```r
bold_tax_name(name = c('Diplura', 'Osmia'))
#>    taxid   taxon tax_rank tax_division parentid       parentname
#> 1 603673 Diplura    genus     Protista    53974 Scytosiphonaceae
#> 2 734358 Diplura    class     Animalia       20       Arthropoda
#> 3   4940   Osmia    genus     Animalia     4962     Megachilinae
#>   specimenrecords taxonrep     representitive_image.image
#> 1               6     <NA>                           <NA>
#> 2             308  Diplura BSOIL/GBOL20120+1540220834.jpg
#> 3            3007    Osmia              BUSA/IMG_3061.jpg
#>   representitive_image.apectratio   input
#> 1                              NA Diplura
#> 2                           0.841 Diplura
#> 3                           1.486   Osmia
```


### Search for taxonomic names via BOLD identifiers

`bold_tax_id` searches for names with BOLD identifiers.


```r
bold_tax_id(id = 88899)
#>   input taxid   taxon tax_rank tax_division parentid parentname
#> 1 88899 88899 Momotus    genus     Animalia    88898  Momotidae
```


```r
bold_tax_id(id = c(88899, 125295))
#>    input  taxid      taxon tax_rank tax_division parentid  parentname
#> 1  88899  88899    Momotus    genus     Animalia    88898   Momotidae
#> 2 125295 125295 Helianthus    genus      Plantae   151101 Asteroideae
#>     taxonrep
#> 1       <NA>
#> 2 Helianthus
```


### Search for sequence data only

The BOLD sequence API gives back sequence data, with a bit of metadata.

The default is to get a list back


```r
bold_seq(taxon = 'Coelioxys')[1:2]
#> [[1]]
#> [[1]]$id
#> [1] "ABEE117-17"
#> 
#> [[1]]$name
#> [1] "Coelioxys elongata"
#> 
#> [[1]]$gene
#> [1] "ABEE117-17"
#> 
#> [[1]]$sequence
#> [1] "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTATCATTATATACATATCATCCTTCCCCATCAGTTGATTTAGCAATTTTTTYTTTACATTTATCAGGAATTTYTTYTATTATCGGATCAATAAATTTTATTGTAACAATTTTAATAATAAAAAATTATTCAATAAATTATAATCAAATACCTTTATTTCCATGATCAATTTTAATTACTACAATTTTATTATTATTATCATTACCTGTATTAGCAGGAGCTATTACAATATTATTATTTGATCGTAATTTAAATTCATCATTTTTTGACCCAATAGGAGGAGGAGATCCTATTTTATATCAACATTTATTTTG------------------------------------"
#> 
#> 
#> [[2]]
#> [[2]]$id
#> [1] "ABEE193-17"
#> 
#> [[2]]$name
#> [1] "Coelioxys rufescens"
#> 
#> [[2]]$gene
#> [1] "ABEE193-17"
#> 
#> [[2]]$sequence
#> [1] "---------------------------------------------------------------------------------------------------------------------------------------------AATAATGATATAATTTATAACTCTTTTATTACAGCTCATGCATTTTTAATAATTTTTTTTTTAGTTATACCTTTTTTAATTGGAGGATTTGGAAATTGATTAGCTCCTTTAATATTAGGAGCTCCAGATATAGCATTCCCTCGAATAAATAATATTAGATTTTGATTATTACCTCCTTCTTTATTAATATTATTAACTAGTAATTTAATTAATCCTAGACCAGGAACAGGATGAACAATTTATCCTCCTTTATCTTTATATAATTATCATCCTTCACCATCAGTAGATTTAGCAATTTTTTCTTTACATTTATCAGGAGTATCATCTATTATTGGTTCAATAAATTTTATTGTAACAATTTTATTAATAAAAAATTATTCAATAAATTATAATCAAATACCTTTATTCCCAKGATCAATTTTAATCACTACAATTTTATTATTATTATCTTTGCCTGTTTTAGCAGGAGCAATTACAATATTATTATTTGATCGAAATCTAAATTCATCCTTTTTTGACCCTTTAGGAGGAGGGGATCCAATTTTATACCAACATTTATTTTGATTTTTTGGACATCC---------------------"
```

You can optionally get back the `crul` response object


```r
res <- bold_seq(taxon = 'Coelioxys', response = TRUE)
res$response_headers
#> $status
#> [1] "HTTP/2 200 "
#> 
#> $server
#> [1] "nginx"
#> 
#> $date
#> [1] "Mon, 20 Apr 2020 16:23:04 GMT"
#> 
#> $`content-type`
#> [1] "application/x-download"
#> 
#> $`x-powered-by`
#> [1] "PHP/5.3.15"
#> 
#> $`content-disposition`
#> [1] "attachment; filename=fasta.fas"
#> 
#> $`x-frame-options`
#> [1] "SAMEORIGIN"
#> 
#> $`x-content-type-options`
#> [1] "nosniff"
#> 
#> $`x-xss-protection`
#> [1] "1; mode=block"
```

You can do geographic searches


```r
bold_seq(geo = "USA")
#> [[1]]
#> [[1]]$id
#> [1] "GBAN1777-08"
#> 
#> [[1]]$name
#> [1] "Macrobdella decora"
#> 
#> [[1]]$gene
#> [1] "GBAN1777-08"
#> 
#> [[1]]$sequence
#> [1] "---------------------------------ATTGGAATCTTGTATTTCTTATTAGGTACATGATCTGCTATAGTAGGGACCTCTATA---AGAATAATTATTCGAATTGAATTAGCTCAACCTGGGTCGTTTTTAGGAAAT---GATCAAATTTACAATACTATTGTTACTGCTCATGGATTAATTATAATTTTTTTTATAGTAATACCTATTTTAATTGGAGGGTTTGGTAATTGATTAATTCCGCTAATA---ATTGGTTCTCCTGATATAGCTTTTCCACGTCTTAATAATTTAAGATTTTGATTACTTCCGCCATCTTTAACTATACTTTTTTGTTCATCTATAGTCGAAAATGGAGTAGGTACTGGATGGACTATTTACCCTCCTTTAGCAGATAACATTGCTCATTCTGGACCTTCTGTAGATATA---GCAATTTTTTCACTTCATTTAGCTGGTGCTTCTTCTATTTTAGGTTCATTAAATTTTATTACTACTGTAGTTAATATACGATGACCAGGGATATCTATAGAGCGAATTCCTTTATTTATTTGATCCGTAATTATTACTACTGTATTGCTATTATTATCTTTACCAGTATTAGCAGCT---GCTATTTCAATATTATTAACAGATCGTAACTTAAATACTAGATTTTTTGACCCAATAGGAGGAGGGGATCCTATTTTATTCCAACATTTATTTTGATTTTTTGGCCACCCTGAAGTTTATATTTTAATTTTACCAGGATTTGGAGCTATTTCTCATGTAGTAAGTCATAACTCT---AAAAAATTAGAACCGTTTGGATCATTAGGGATATTATATGCAATAATTGGAATTGCAATTTTAGGTTTTATTGTTTGAGCACATCATATATTTACAGTAGGTCTTGATGTAGATACACGAGCTTATTTTACAGCAGCTACAATAGTTATTGCTGTTCCTACAGGAATTAAAGTATTTAGGTGATTG---GCAACT"
#> 
#> 
#> [[2]]
#> [[2]]$id
#> [1] "GBAN1780-08"
#> 
#> [[2]]$name
#> [1] "Haemopis terrestris"
#> 
#> [[2]]$gene
#> [1] "GBAN1780-08"
#> 
#> [[2]]$sequence
#> [1] "---------------------------------ATTGGAACWTTWTATTTTATTTTNGGNGCTTGATCTGCTATATTNGGGATCTCAATA---AGGAATATTATTCGAATTGAGCCATCTCAACCTGGGAGATTATTAGGAAAT---GATCAATTATATAATTCATTAGTAACAGCTCATGGATTAATTATAATTTTCTTTATGGTTATGCCTATTTTGATTGGTGGGTTTGGTAATTGATTACTACCTTTAATA---ATTGGAGCCCCTGATATAGCTTTTCCTCGATTAAATAATTTAAGTTTTTGATTATTACCACCTTCATTAATTATATTGTTAAGATCCTCTATTATTGAAAGAGGGGTAGGTACAGGTTGAACCTTATATCCTCCTTTAGCAGATAGATTATTTCATTCAGGTCCATCGGTAGATATA---GCTATTTTTTCATTACATATAGCTGGAGCATCATCTATTTTAGGCTCATTAAACTTTATTTCTACAATTATTAATATACGAATTAAAGGTATAAGATCTGATCGAGTACCTTTATTTGTATGATCAGTTGTTATTACAACAGTTCTGTTATTATTGTCTTTACCTGTTTTAGCTGCA---GCTATTACTATATTATTAACAGATCGTAATTTAAATACTACTTTTTTTGATCCTATAGGAGGTGGAGATCCAGTATTGTTTCAACACTTATTTTGATTTTTTGGTCATCCAGAAGTATATATTTTGATTTTACCAGGATTTGGAGCAATTTCTCATATTATTACAAATAATTCT---AAAAAATTGGAACCTTTTGGATCTCTTGGTATAATTTATGCTATAATTGGAATTGCAGTTTTAGGGTTTATTGTATGAGCCCATCATATATTTACTGTAGGATTAGATGTTGATACTCGAGCTTATTTTACAGCAGCTACTATAGTTATTGCTGTTCCTACTGGTATTAAAGTTTTTAGGTGATTA---GCAACA"
#> 
#> 
#> [[3]]
#> [[3]]$id
#> [1] "GBNM0293-06"
#> 
#> [[3]]$name
#> [1] "Steinernema carpocapsae"
#> 
#> [[3]]$gene
#> [1] "GBNM0293-06"
#> 
#> [[3]]$sequence
#> [1] "---------------------------------------------------------------------------------ACAAGATTATCTCTTATTATTCGTTTAGAGTTGGCTCAACCTGGTCTTCTTTTGGGTAAT---GGTCAATTATATAATTCTATTATTACTGCTCATGCTATTCTTATAATTTTTTTCATAGTTATACCTAGAATAATTGGTGGTTTTGGTAATTGAATATTACCTTTAATATTGGGGGCTCCTGATATAAGTTTTCCACGTTTGAATAATTTAAGTTTTTGATTGCTACCAACTGCTATATTTTTGATTTTAGATTCTTGTTTTGTTGACACTGGTTGTGGTACTAGTTGAACTGTTTATCCTCCTTTGAGG---ACTTTAGGTCACCCTGGYAGAAGTGTAGATTTAGCTATTTTTAGTCTTCATTGTGCAGGAATTAGCTCAATTTTAGGGGCTATTAATTTTATATGTACTACAAAAAATCTTCGTAGTAGTTCTATTTCTTTGGAACATATAAGACTTTTTGTTTGGGCTGTTTTTGTTACTGTTTTTTTATTAGTTTTATCTTTACCTGTTTTAGCTGGTGCTATTACTATGCTTTTAACAGACCGTAATTTAAATACTTCTTTTTTT------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
#> 
#> 
#> [[4]]
#> [[4]]$id
#> [1] "NEONV108-11"
#> 
#> [[4]]$name
#> [1] "Aedes thelcter"
#> 
#> [[4]]$gene
#> [1] "NEONV108-11"
#> 
#> [[4]]$sequence
#> [1] "AACTTTATACTTCATCTTCGGAGTTTGATCAGGAATAGTTGGTACATCATTAAGAATTTTAATTCGTGCTGAATTAAGTCAACCAGGTATATTTATTGGAAATGACCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTATTATAATTTTCTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGACTAGTTCCTCTAATATTAGGAGCCCCAGATATAGCTTTCCCTCGAATAAATAATATAAGTTTTTGAATACTACCTCCCTCATTAACTCTTCTACTTTCAAGTAGTATAGTAGAAAATGGATCAGGAACAGGATGAACAGTTTATCCACCTCTTTCATCTGGAACTGCTCATGCAGGAGCCTCTGTTGATTTAACTATTTTTTCTCTTCATTTAGCCGGAGTTTCATCAATTTTAGGGGCTGTAAATTTTATTACTACTGTAATTAATATACGATCTGCAGGAATTACTCTTGATCGACTACCTTTATTCGTTTGATCTGTAGTAATTACAGCTGTTTTATTACTTCTTTCACTTCCTGTATTAGCTGGAGCTATTACAATACTATTAACTGATCGAAATTTAAATACATCTTTCTTTGATCCAATTGGAGGAGGAGACCCAATTTTATACCAACATTTATTT"
#> 
#> 
#> [[5]]
#> [[5]]$id
#> [1] "NEONV109-11"
#> 
#> [[5]]$name
#> [1] "Aedes thelcter"
#> 
#> [[5]]$gene
#> [1] "NEONV109-11"
#> 
#> [[5]]$sequence
#> [1] "AACTTTATACTTCATCTTCGGAGTTTGATCAGGAATAGTTGGTACATCATTAAGAATTTTAATTCGTGCTGAATTAAGTCAACCAGGTATATTTATTGGAAATGACCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTATTATAATTTTCTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGACTAGTTCCTCTAATATTAGGAGCCCCAGATATAGCTTTCCCTCGAATAAATAATATAAGTTTTTGAATACTACCTCCCTCATTAACTCTTCTACTTTCAAGTAGTATAGTAGAAAATGGGTCAGGAACAGGATGAACAGTTTATCCACCTCTTTCATCTGGAACTGCTCATGCAGGAGCCTCTGTTGATTTAACTATTTTTTCTCTTCATTTAGCCGGAGTTTCATCAATTTTAGGGGCTGTAAATTTTATTACTACTGTAATTAATATACGATCTGCAGGAATTACTCTTGATCGACTACCTTTATTCGTTTGATCTGTAGTAATTACAGCTGTTTTATTACTTCTTTCACTTCCTGTATTAGCTGGAGCTATTACAATACTATTAACTGATCGAAATTTAAATACATCTTTCTTTGACCCAATTGGAGGGGGAGACCCAATTTTATACCAACATTTATTT"
```

And you can search by researcher name


```r
bold_seq(researchers = 'Thibaud Decaens')[[1]]
#> $id
#> [1] "COLNO015-09"
#> 
#> $name
#> [1] "Ptiliidae"
#> 
#> $gene
#> [1] "COLNO015-09"
#> 
#> $sequence
#> [1] "AACCTTGTATTTTATGTTCGGNGCTTGAGCTGGAATAGTCGGGACAAGTTTGAGTCTCCTTATCCGAACTGAACTCGGCACTCCAGGTTCACTAATTGGAGACGACCAAATCTACAACGTAATCGTAACAGCTCATGCTTTTGTGATGATTTTTTTTATGGTCATGCCGATTTTAATCGGAGGTTTCGGAAACTGACTTGTTCCCTTGATACTTGGAGCCCCTGATATAGCTTTCCCTCGAATGAACAACATAAGGTTCTGGCTCCTCCCCCCTTCTCTGACCCTACTTTTAATGAGAAGGATAGTAGAAAGCGGAGCAGGAACAGGGTGAACAGTTTATCCTCCCCTAGCTTCAAATATTGCTCATGGTGGAGCATCAGTTGATTTGGCGATTTTTAGGCTCCATTTAGCTGGAATCTCTTCCATTCTAGGAGCCGTAAATTTCATCACAACAATTCTCAACATACGAACTCCTCAAATAAGGTTTGATCAAATGCCATTGTTTGTTTGGGCTGTTGGAATTACAGCTCTCCTTCTTCTTCTTTCACTTCCNGTTTTAGCCGGAGCTATCACAATATTGCTAAC-------------------------------------------------------------------------"
```

by taxon IDs


```r
bold_seq(ids = c('ACRJP618-11', 'ACRJP619-11'))
#> [[1]]
#> [[1]]$id
#> [1] "ACRJP618-11"
#> 
#> [[1]]$name
#> [1] "Lepidoptera"
#> 
#> [[1]]$gene
#> [1] "ACRJP618-11"
#> 
#> [[1]]$sequence
#> [1] "------------------------TTGAGCAGGCATAGTAGGAACTTCTCTTAGTCTTATTATTCGAACAGAATTAGGAAATCCAGGATTTTTAATTGGAGATGATCAAATCTACAATACTATTGTTACGGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAATTGATTAGTTCCCCTTATACTAGGAGCCCCAGATATAGCTTTCCCTCGAATAAACAATATAAGTTTTTGGCTTCTTCCCCCTTCACTATTACTTTTAATTTCCAGAAGAATTGTTGAAAATGGAGCTGGAACTGGATGAACAGTTTATCCCCCACTGTCATCTAATATTGCCCATAGAGGTACATCAGTAGATTTAGCTATTTTTTCTTTACATTTAGCAGGTATTTCCTCTATTTTAGGAGCGATTAATTTTATTACTACAATTATTAATATACGAATTAACAGTATAAATTATGATCAAATACCACTATTTGTGTGATCAGTAGGAATTACTGCTTTACTCTTATTACTTTCTCTTCCAGTATTAGCAGGTGCTATCACTATATTATTAACGGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGGAGATCCAATTTTATATCAACATTTATTT"
#> 
#> 
#> [[2]]
#> [[2]]$id
#> [1] "ACRJP619-11"
#> 
#> [[2]]$name
#> [1] "Lepidoptera"
#> 
#> [[2]]$gene
#> [1] "ACRJP619-11"
#> 
#> [[2]]$sequence
#> [1] "AACTTTATATTTTATTTTTGGTATTTGAGCAGGCATAGTAGGAACTTCTCTTAGTCTTATTATTCGAACAGAATTAGGAAATCCAGGATTTTTAATTGGAGATGATCAAATCTACAATACTATTGTTACGGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAATTGATTAGTTCCCCTTATACTAGGAGCCCCAGATATAGCTTTCCCTCGAATAAACAATATAAGTTTTTGGCTTCTTCCCCCTTCACTATTACTTTTAATTTCCAGAAGAATTGTTGAAAATGGAGCTGGAACTGGATGAACAGTTTATCCCCCACTGTCATCTAATATTGCCCATAGAGGTACATCAGTAGATTTAGCTATTTTTTCTTTACATTTAGCAGGTATTTCCTCTATTTTAGGAGCGATTAATTTTATTACTACAATTATTAATATACGAATTAACAGTATAAATTATGATCAAATACCACTATTTGTGTGATCAGTAGGAATTACTGCTTTACTCTTATTACTTTCTCTTCCAGTATTAGCAGGTGCTATCACTATATTATTAACGGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGGAGATCCAATTTTATATCAACATTTATTT"
```

by container (containers include project codes and dataset codes)


```r
bold_seq(container = 'ACRJP')[[1]]
#> $id
#> [1] "ACRJP003-09"
#> 
#> $name
#> [1] "Lepidoptera"
#> 
#> $gene
#> [1] "ACRJP003-09"
#> 
#> $sequence
#> [1] "AACATTATATTTTATTTTTGGGATCTGATCTGGAATAGTAGGGACATCTTTAAGTATACTAATTCGAATAGAACTAGGAAATCCTGGATGTTTAATTGGGGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCCATTATAATTGGAGGTTTTGGCAATTGACTTGTACCATTAATATTAGGAGCCCCTGATATAGCATTTCCCCGAATAAATAATATAAGATTTTGACTTCTTCCCCCCTCATTAATTTTATTAATTTCAAGAAGAATTGTTGAAAATGGAGCAGGAACAGGATGAACAGTCTATCCTCCATTATCTTCTAATATTGCGCATAGAGGATCCTCTGTTGATTTAGCTATTTTCTCACTTCATTTAGCAGGAATTTCTTCTATTTTAGGAGCAATTAATTTTATTACAACTATTATTAATATACGAATAAATAATTTACTTTTTGACCAAATACCTCTATTTGTTTGAGCAGTAGGTATTACAGCTGTTCTTCTTTTATTATCATTACCAGTATTAGCAGGAGCAATTACCATACTATTAACAGATCGTAATTTAAATACTTCTTTCTTTGATCCTGCTGGAGGAGGAGATCCAATTTTATACCAACATTTATTT"
```

by bin (a bin is a _Barcode Index Number_)


```r
bold_seq(bin = 'BOLD:AAA5125')[[1]]
#> $id
#> [1] "ASARD6776-12"
#> 
#> $name
#> [1] "Eacles ormondei"
#> 
#> $gene
#> [1] "ASARD6776-12"
#> 
#> $sequence
#> [1] "AACTTTATATTTTATTTTTGGAATTTGAGCAGGTATAGTAGGAACTTCTTTAAGATTACTAATTCGAGCAGAATTAGGTACCCCCGGATCTTTAATTGGAGATGACCAAATTTATAATACCATTGTAACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGATTAGTACCCCTAATACTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGACTATTACCCCCATCTTTAACCCTTTTAATTTCTAGAAGAATTGTCGAAAATGGAGCTGGAACTGGATGAACAGTTTATCCCCCCCTTTCATCTAATATTGCTCATGGAGGCTCTTCTGTTGATTTAGCTATTTTTTCCCTTCATCTAGCTGGAATCTCATCAATTTTAGGAGCTATTAATTTTATCACAACAATCATTAATATACGACTAAATAATATAATATTTGACCAAATACCTTTATTTGTATGAGCTGTTGGTATTACAGCATTTCTTTTATTGTTATCTTTACCTGTACTAGCTGGAGCTATTACTATACTTTTAACAGATCGAAACTTAAATACATCATTTTTTGACCCAGCAGGAGGAGGAGATCCTATTCTCTATCAACATTTATTT"
```

And there are more ways to query, check out the docs for `?bold_seq`.


### Search for specimen data only

The BOLD specimen API doesn't give back sequences, only specimen data. By default you download `tsv` format data, which is given back to you as a `data.frame`


```r
res <- bold_specimens(taxon = 'Osmia')
head(res[,1:8])
#>      processid   sampleid recordID catalognum   fieldnum
#> 1  BEECA373-06 05-NT-0373   514740            05-NT-0373
#> 2  BEECA607-06 05-NT-0607   516959            05-NT-0607
#> 3  BEECA963-07 01-OR-0790   554153            01-OR-0790
#> 4  BEECB358-07 04-WA-1076   596920 BBSL697174 04-WA-1076
#> 5  BEECB438-07 00-UT-1157   597000 BBSL432653 00-UT-1157
#> 6 BEECC1176-09 02-UT-2849  1060879 BBSL442586 02-UT-2849
#>                    institution_storing collection_code      bin_uri
#> 1   York University, Packer Collection              NA BOLD:AAI2013
#> 2   York University, Packer Collection              NA BOLD:AAC8510
#> 3   York University, Packer Collection              NA BOLD:ABZ3184
#> 4 Utah State University, Logan Bee Lab              NA BOLD:AAC5797
#> 5 Utah State University, Logan Bee Lab              NA BOLD:AAF2159
#> 6   York University, Packer Collection              NA BOLD:AAE5368
```

You can optionally get back the data in `XML` format


```r
bold_specimens(taxon = 'Osmia', format = 'xml')
```


```r
<?xml version="1.0" encoding="UTF-8"?>
<bold_records  xsi:noNamespaceSchemaLocation="http://www.boldsystems.org/schemas/BOLDPublic_record.xsd"  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <record>
    <record_id>1470124</record_id>
    <processid>BOM1525-10</processid>
    <bin_uri>BOLD:AAN3337</bin_uri>
    <specimen_identifiers>
      <sampleid>DHB 1011</sampleid>
      <catalognum>DHB 1011</catalognum>
      <fieldnum>DHB1011</fieldnum>
      <institution_storing>Marjorie Barrick Museum</institution_storing>
    </specimen_identifiers>
    <taxonomy>
```

You can choose to get the `crul` response object back if you'd rather work with the raw data returned from the BOLD API.


```r
res <- bold_specimens(taxon = 'Osmia', format = 'xml', response = TRUE)
res$url
#> [1] "https://v4.boldsystems.org/index.php/API_Public/specimen?taxon=Osmia&format=xml"
res$status_code
#> [1] 200
res$response_headers
#> $status
#> [1] "HTTP/2 200 "
#> 
#> $server
#> [1] "nginx"
#> 
#> $date
#> [1] "Mon, 20 Apr 2020 16:23:41 GMT"
#> 
#> $`content-type`
#> [1] "application/x-download"
#> 
#> $`x-powered-by`
#> [1] "PHP/5.3.15"
#> 
#> $`content-disposition`
#> [1] "attachment; filename=bold_data.xml"
#> 
#> $`x-frame-options`
#> [1] "SAMEORIGIN"
#> 
#> $`x-content-type-options`
#> [1] "nosniff"
#> 
#> $`x-xss-protection`
#> [1] "1; mode=block"
```

### Search for specimen plus sequence data

The specimen/sequence combined API gives back specimen and sequence data. Like the specimen API, this one gives by default `tsv` format data, which is given back to you as a `data.frame`. Here, we're setting `sepfasta=TRUE` so that the sequence data is given back as a list, and taken out of the `data.frame` returned so the `data.frame` is more manageable.


```r
res <- bold_seqspec(taxon = 'Osmia', sepfasta = TRUE)
res$fasta[1:2]
#> $`BEECA373-06`
#> [1] "-ATTTTATATATAATTTTTGCTATATGATCAGGTATAATCGGATCAGCAATAAGAATTATTATTCGTATAGAATTAAGAATTCCTGGTTCATGAATTTCAAATGATCAAACTTATAACTCTTTAGTAACTGCTCATGCTTTTTTAATAATTTTTTTCTTAGTTATACCTTTTTTAATTGGAGGATTTGGAAATTGATTAATTCCTTTAATATTAGGAATCCCGGATATAGCTTTCCCTCGAATAAATAATATTAGATTTTGACTTTTACCCCCTTCATTAATATTATTACTTTTAAGAAATTTTATAAATCCAAGACCAGGTACTGGATGAACTGTTTATCCTCCTCTTTCTTCTCATTTATTTCATTCTTCTCCTTCAGTTGATATAGCCATTTTTTCTTTACATATTTCCGGTTTATCTTCTATTATAGGTTCGTTAAATTTTATTGTTACAATTATTATAATAAAAAATATTTCTTTAAAACATATCCAATTACCTTTATTTCCATGATCTGTTTTTATTACTACTATCTTATTACTTTTTTCTTTACCTGTTTTAGCAGGAGCTATTACTATATTATTATTTGATCGAAATTTTAATACTTCATTTTTTGATCCTACAGGAGGTGGAGATCCAATCCTTTATCAACATTTATTT"
#> 
#> $`BEECA607-06`
#> [1] "AATATTATATATAATTTTTGCTTTGTGATCTGGAATAATTGGTTCATCTATAAGAATTATTATTCGTATAGAATTAAGAATTCCTGGTTCATGAATTTCAAATGATCAAGTTTATAATTCATTAGTTACAGCTCATGCTTTTTTAATAATTTTTTTTTTAGTTATACCATTTATAATTGGAGGATTTGGAAATTGATTAGTTCCTTTAATATTAGGAATTCCTGATATAGCTTTTCCTCGAATAAATAATATTAGATTTTGATTATTACCACCATCATTAATACTTTTACTTTTAAGAAATTTTTTTAATCCAAGTTCAGGAACTGGATGAACTGTATATCCTCCTCTTTCATCATATTTATTTCATTCTTCACCTTCTGTTGATTTAGCTATTTTTTCTCTTCATATATCAGGTTTATCTTCTATTATAGGTTCATTAAACTTTATTGTAACTATTATTATAATAAAAAATATTTCTTTAAAGTATATTCAATTGCCATTATTTCCATGATCTGTTTTTATTACTACAATTCTTTTATTATTATCATTACCAGTTTTAGCAGGTGCTATTACTATATTATTATTTGATCGAAATTTTAATACTTCATTTTTTGATCCTACAGGAGGGGGAG--------------------------"
```

Or you can index to a specific sequence like


```r
res$fasta['GBAH0293-06']
#> $`GBAH0293-06`
#> [1] "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAATGTTAGGGATTCCAGATATAGCTTTTCCACGAATAAATAATATTAGATTTTGACTGTTACCTCCATCTTTAATATTATTACTTTTAAGAAATTTTTTAAATCCAAGTCCTGGAACAGGATGAACAGTTTATCCTCCTTTATCATCAAATTTATTTCATTCTTCTCCTTCAGTTGATTTAGCAATTTTTTCTTTACATATTTCAGGTTTATCTTCTATTATAGGTTCATTAAATTTTATTGTTACAATTATTATAATAAAAAATATTTCTTTAAAATATATTCAATTACCTTTATTTTCTTGATCTGTATTTATTACTACTATTCTTTTATTATTTTCTTTACCTGTATTAGCTGGAGCTATTACTATATTATTATTTGATCGAAATTTTAATACATCTTTTTTTGATCCAACAGGAGGGGGAGATCCAATTCTTTATCAACATTTATTTTGATTTTTTGGTCATCCTGAAGTTTATATTTTAATTTTACCTGGATTTGGATTAATTTCTCAAATTATTTCTAATGAAAGAGGAAAAAAAGAAACTTTTGGAAATATTGGTATAATTTATGCTATATTAAGAATTGGACTTTTAGGTTTTATTGTT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
```

### Get trace files

This function downloads files to your machine - it does not load them into your R session - but prints out where the files are for your information.


```r
bold_trace(taxon = 'Osmia', quiet = TRUE)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_tax_name.R
\name{bold_tax_name}
\alias{bold_tax_name}
\title{Search BOLD for taxonomy data by taxonomic name}
\usage{
bold_tax_name(name, fuzzy = FALSE, response = FALSE, ...)
}
\arguments{
\item{name}{(character) One or more scientific names. required.}

\item{fuzzy}{(logical) Whether to use fuzzy search or not (default: \code{FALSE})}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\description{
Search BOLD for taxonomy data by taxonomic name
}
\details{
The \code{dataTypes} parameter is not supported in this function.
If you want to use that parameter, get an ID from this function and pass
it into \code{bold_tax_id}, and then use the \code{dataTypes} parameter.
}
\examples{
\dontrun{
bold_tax_name(name='Diplura')
bold_tax_name(name='Osmia')
bold_tax_name(name=c('Diplura','Osmia'))
bold_tax_name(name=c("Apis","Puma concolor","Pinus concolor"))
bold_tax_name(name='Diplur', fuzzy=TRUE)
bold_tax_name(name='Osm', fuzzy=TRUE)

## get http response object only
bold_tax_name(name='Diplura', response=TRUE)
bold_tax_name(name=c('Diplura','Osmia'), response=TRUE)

## Names with no data in BOLD database
bold_tax_name("Nasiaeshna pentacantha")
bold_tax_name(name = "Cordulegaster erronea")
bold_tax_name(name = "Cordulegaster erronea", response=TRUE)

## curl debugging
bold_tax_name(name='Diplura', verbose = TRUE)
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=taxonomy
}
\seealso{
\code{\link[=bold_tax_id]{bold_tax_id()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_identify.R
\name{bold_identify}
\alias{bold_identify}
\title{Search for matches to sequences against the BOLD COI database.}
\usage{
bold_identify(sequences, db = "COX1", response = FALSE, ...)
}
\arguments{
\item{sequences}{(character) Returns all records containing matching marker
codes. Required. One or more. See Details.}

\item{db}{(character) The database to match against, one of COX1,
COX1_SPECIES, COX1_SPECIES_PUBLIC, OR COX1_L604bp. See Details for
more information.}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}, main
purpose being curl debugging

BOLD only allows one sequence per query. We internally \code{lapply}
over the input values given to the sequences` parameter to search
with one sequence per query. Remember this if you have a lot of sequences -
you are doing a separate query for each one, so it can take a long time -
if you run into errors let us know.}
}
\value{
A data.frame with details for each specimen matched. if a
failed request, returns \code{NULL}
}
\description{
Search for matches to sequences against the BOLD COI database.
}
\section{db parmeter options}{

\itemize{
\item COX1 Every COI barcode record on BOLD with a minimum sequence
length of 500bp (warning: unvalidated library and includes records without
species level identification). This includes many species represented by
only one or two specimens as well as all species with interim taxonomy. This
search only returns a list of the nearest matches and does not provide a
probability of placement to a taxon.
\item COX1_SPECIES Every COI barcode record with a species level
identification and a minimum sequence length of 500bp. This includes
many species represented by only one or two specimens as well as  all
species with interim taxonomy.
\item COX1_SPECIES_PUBLIC All published COI records from BOLD and GenBank
with a minimum sequence length of 500bp. This library is a collection of
records from the published projects section of BOLD.
\item OR COX1_L604bp Subset of the Species library with a minimum sequence
length of 640bp and containing both public and private records. This library
is intended for short sequence identification as it provides maximum overlap
with short reads from the barcode region of COI.
}
}

\section{Named outputs}{

To maintain names on the output list of data make sure to pass in a
named list to the \code{sequences} parameter. You can for example,
take a list of sequences, and use \code{\link[stats:setNames]{stats::setNames()}} to set names.
}

\examples{
\dontrun{
seq <- sequences$seq1
res <- bold_identify(sequences=seq)
head(res[[1]])
head(bold_identify(sequences=seq, db='COX1_SPECIES')[[1]])
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=idengine
}
\seealso{
\code{\link[=bold_identify_parents]{bold_identify_parents()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_identify_parents.R
\name{bold_identify_parents}
\alias{bold_identify_parents}
\title{Add taxonomic parent names to a data.frame}
\usage{
bold_identify_parents(
  x,
  wide = FALSE,
  taxid = NULL,
  taxon = NULL,
  tax_rank = NULL,
  tax_division = NULL,
  parentid = NULL,
  parentname = NULL,
  taxonrep = NULL,
  specimenrecords = NULL,
  ...
)
}
\arguments{
\item{x}{(data.frame/list) list of data.frames - the output from a call to
\code{\link[=bold_identify]{bold_identify()}}. or a single data.frame from the output from same.
required.}

\item{wide}{(logical) output in long or wide format. See Details.
Default: \code{FALSE}}

\item{taxid}{(character) A taxid name. Optional. See \code{Filtering} below.}

\item{taxon}{(character) A taxon name. Optional. See \code{Filtering} below.}

\item{tax_rank}{(character) A tax_rank name. Optional. See \code{Filtering}
below.}

\item{tax_division}{(character) A tax_division name. Optional. See
\code{Filtering} below.}

\item{parentid}{(character) A parentid name. Optional. See \code{Filtering}
below.}

\item{parentname}{(character) A parentname name. Optional. See \code{Filtering}
below.}

\item{taxonrep}{(character) A taxonrep name. Optional. See \code{Filtering}
below.}

\item{specimenrecords}{(character) A specimenrecords name. Optional.
See \code{Filtering} below.}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}, main
purpose being curl debugging}
}
\value{
a list of the same length as the input
}
\description{
Add taxonomic parent names to a data.frame
}
\details{
This function gets unique set of taxonomic names from the input
data.frame, then queries \code{\link[=bold_tax_name]{bold_tax_name()}} to get the
taxonomic ID, passing it to \code{\link[=bold_tax_id]{bold_tax_id()}} to get the parent
names, then attaches those to the input data.

Records in the input data that do not have matches for parent names
simply get NA values in the added columns.
}
\section{Filtering}{

The parameters \code{taxid}, \code{taxon}, \code{tax_rank}, \code{tax_division},
\code{parentid}, \code{parentname},\code{taxonrep}, and \code{specimenrecords} are not used
in the search sent to BOLD, but are used in filtering the data down
to a subset that is closer to the target you want. For all these
parameters, you can use regex strings since we use \code{\link[=grep]{grep()}} internally
to match. Filtering narrows down to the set that matches your query,
and removes the rest. The data.frame that we filter on with these
parameters internally is the result of a call to the \code{\link[=bold_tax_name]{bold_tax_name()}}
function.
}

\section{wide vs long format}{

When \code{wide = FALSE} you get many rows for each record. Essentially,
we \code{cbind} the taxonomic classification onto the one row from the
result of \code{\link[=bold_identify]{bold_identify()}}, giving as many rows as there are
taxa in the taxonomic classification.

When \code{wide = TRUE} you get one row for each record - thus the
dimensions of the input data stay the same. For this option, we take just
the rows for taxonomic ID and name for each taxon in the taxonomic
classification, and name the columns by the taxon rank, so you get
\code{phylum} and \code{phylum_id}, and so on.
}

\examples{
\dontrun{
df <- bold_identify(sequences = sequences$seq2)

# long format
out <- bold_identify_parents(df)
str(out)
head(out[[1]])

# wide format
out <- bold_identify_parents(df, wide = TRUE)
str(out)
head(out[[1]])

x <- bold_seq(taxon = "Satyrium")
out <- bold_identify(c(x[[1]]$sequence, x[[13]]$sequence))
res <- bold_identify_parents(out)
res

x <- bold_seq(taxon = 'Diplura')
out <- bold_identify(vapply(x, "[[", "", "sequence")[1:20])
res <- bold_identify_parents(out)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold-package.R
\docType{package}
\name{bold-package}
\alias{bold-package}
\alias{bold}
\title{bold}
\description{
bold: A programmatic interface to the Barcode of Life data
}
\section{About}{


This package gives you access to data from BOLD System
http://www.boldsystems.org/ via their API
(http://v4.boldsystems.org/index.php/api_home)
}

\section{Functions}{

\itemize{
\item \code{\link[=bold_specimens]{bold_specimens()}} - Search for specimen data
\item \code{\link[=bold_seq]{bold_seq()}} - Search for and retrieve sequences
\item \code{\link[=bold_seqspec]{bold_seqspec()}} - Get sequence and specimen data together
\item \code{\link[=bold_trace]{bold_trace()}} - Get trace files - saves to disk
\item \code{\link[=read_trace]{read_trace()}} - Read trace files into R
\item \code{\link[=bold_tax_name]{bold_tax_name()}} - Get taxonomic names via input names
\item \code{\link[=bold_tax_id]{bold_tax_id()}} - Get taxonomic names via BOLD identifiers
\item \code{\link[=bold_identify]{bold_identify()}} - Search for match given a COI sequence
}

Interestingly, they provide xml and tsv format data for the specimen data,
while  they provide fasta data format for the sequence data. So for the
specimen data  you can get back raw XML, or a data frame parsed from the
tsv data, while for sequence data you get back a list (b/c sequences are
quite long and would make a data frame unwieldy).
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_stats.R
\name{bold_stats}
\alias{bold_stats}
\title{Get BOLD stats}
\usage{
bold_stats(
  taxon = NULL,
  ids = NULL,
  bin = NULL,
  container = NULL,
  institutions = NULL,
  researchers = NULL,
  geo = NULL,
  dataType = "drill_down",
  response = FALSE,
  ...
)
}
\arguments{
\item{taxon}{(character) Returns all records containing matching taxa. Taxa
includes the ranks of  phylum, class, order, family, subfamily, genus,
and species.}

\item{ids}{(character) Returns all records containing matching IDs. IDs
include Sample IDs, Process IDs, Museum IDs and Field IDs.}

\item{bin}{(character) Returns all records contained in matching BINs. A
BIN is defined by a Barcode Index Number URI.}

\item{container}{(character) Returns all records contained in matching
projects or datasets. Containers include project codes and dataset codes}

\item{institutions}{(character) Returns all records stored in matching
institutions. Institutions are the Specimen Storing Site.}

\item{researchers}{(character) Returns all records containing matching
researcher names. Researchers include collectors and specimen identifiers.}

\item{geo}{(character) Returns all records collected in matching geographic
sites. Geographic sites includes countries and province/states.}

\item{dataType}{(character) one of "overview" or "drill_down" (default).
"drill_down": a detailed summary of information which provides record
counts by (BINs, Country, Storing Institution, Species). "overview":
the total counts of (BINs, Countries, Storing Institutions, Orders,
Families, Genus, Species)}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\description{
Get BOLD stats
}
\examples{
\dontrun{
x <- bold_stats(taxon='Osmia')
x$total_records
x$records_with_species_name
x$bins
x$countries
x$depositories
x$order
x$family
x$genus
x$species

# just get all counts
lapply(Filter(is.list, x), "[[", "count")

res <- bold_stats(taxon='Osmia', response=TRUE)
res$url
res$status_code
res$response_headers

# More than 1 can be given for all search parameters
bold_stats(taxon=c('Coelioxys','Osmia'))

## curl debugging
### These examples below take a long time, so you can set a timeout so that
### it stops by X sec
bold_stats(taxon='Osmia', verbose = TRUE)
# bold_stats(geo='Costa Rica', timeout_ms = 6)
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=webservices
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_seq.R
\name{bold_seq}
\alias{bold_seq}
\title{Search BOLD for sequences.}
\usage{
bold_seq(
  taxon = NULL,
  ids = NULL,
  bin = NULL,
  container = NULL,
  institutions = NULL,
  researchers = NULL,
  geo = NULL,
  marker = NULL,
  response = FALSE,
  ...
)
}
\arguments{
\item{taxon}{(character) Returns all records containing matching taxa. Taxa
includes the ranks of  phylum, class, order, family, subfamily, genus,
and species.}

\item{ids}{(character) Returns all records containing matching IDs. IDs
include Sample IDs, Process IDs, Museum IDs and Field IDs.}

\item{bin}{(character) Returns all records contained in matching BINs. A
BIN is defined by a Barcode Index Number URI.}

\item{container}{(character) Returns all records contained in matching
projects or datasets. Containers include project codes and dataset codes}

\item{institutions}{(character) Returns all records stored in matching
institutions. Institutions are the Specimen Storing Site.}

\item{researchers}{(character) Returns all records containing matching
researcher names. Researchers include collectors and specimen identifiers.}

\item{geo}{(character) Returns all records collected in matching geographic
sites. Geographic sites includes countries and province/states.}

\item{marker}{(character) Returns all records containing matching
marker codes.}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\value{
A list with each element of length 4 with slots for id, name,
gene, and sequence.
}
\description{
Get sequences for a taxonomic name, id, bin, container, institution,
researcher, geographic, place, or gene.
}
\section{Large requests}{

Some requests can lead to errors. These often have to do with requesting
data for a rank that is quite high in the tree, such as an Order,
for example, Coleoptera. If your request is taking a long time,
it's likely that something will go wrong on the BOLD server side,
or we'll not be able to parse the result here in R because
R can only process strings of a certain length. \code{bold}
users have reported errors in which the resulting response from
BOLD is so large that we could not parse it.

A good strategy for when you want data for a high rank is to
do many separate requests for lower ranks within your target
rank. You can do this manually, or use the function
\code{taxize::downstream} to get all the names of a lower
rank within a target rank. There's an example in the README
(https://docs.ropensci.org/bold/#large-data)
}

\section{If a request times out}{

This is likely because you're request was for a large number of
sequences and the BOLD service timed out. You still should get
some output, those sequences that were retrieved before the time
out happened. As above, see the README
(https://docs.ropensci.org/bold/#large-data) for an example of
dealing with large data problems with this function.
}

\section{Marker}{

Notes from BOLD on the \code{marker} param:
"All markers for a specimen matching the search string will be returned.
ie. A record with COI-5P and ITS will return sequence data for both
markers even if only COI-5P was specified."

You will likely end up with data with markers that you did not request -
just be sure to filter those out as needed.
}

\examples{
\dontrun{
res <- bold_seq(taxon='Coelioxys')
bold_seq(taxon='Aglae')
bold_seq(taxon=c('Coelioxys','Osmia'))
bold_seq(ids='ACRJP618-11')
bold_seq(ids=c('ACRJP618-11','ACRJP619-11'))
bold_seq(bin='BOLD:AAA5125')
bold_seq(container='ACRJP')
bold_seq(researchers='Thibaud Decaens')
bold_seq(geo='Ireland')
bold_seq(geo=c('Ireland','Denmark'))

# Return the http response object for detailed Curl call response details
res <- bold_seq(taxon='Coelioxys', response=TRUE)
res$url
res$status_code
res$response_headers

## curl debugging
### You can do many things, including get verbose output on the curl 
### call, and set a timeout
bold_seq(taxon='Coelioxys', verbose = TRUE)[1:2]
# bold_seqspec(taxon='Coelioxys', timeout_ms = 10)
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=webservices
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold-package.R
\docType{data}
\name{sequences}
\alias{sequences}
\title{List of 3 nucleotide sequences to use in examples for the
\code{\link[=bold_identify]{bold_identify()}} function}
\description{
List of 3 nucleotide sequences to use in examples for the
\code{\link[=bold_identify]{bold_identify()}} function
}
\details{
Each sequence is a character string, of lengths 410, 600, and 696.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_tax_id.R
\name{bold_tax_id}
\alias{bold_tax_id}
\title{Search BOLD for taxonomy data by BOLD ID.}
\usage{
bold_tax_id(
  id,
  dataTypes = "basic",
  includeTree = FALSE,
  response = FALSE,
  ...
)
}
\arguments{
\item{id}{(integer) One or more BOLD taxonomic identifiers. required.}

\item{dataTypes}{(character) Specifies the datatypes that will be
returned. 'all' returns all data. 'basic' returns basic taxon information.
'images' returns specimen images.}

\item{includeTree}{(logical) If \code{TRUE} (default: \code{FALSE}), returns a list
containing information for parent taxa as well as the specified taxon.}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\description{
Search BOLD for taxonomy data by BOLD ID.
}
\examples{
\dontrun{
bold_tax_id(id=88899)
bold_tax_id(id=88899, includeTree=TRUE)
bold_tax_id(id=88899, includeTree=TRUE, dataTypes = "stats")
bold_tax_id(id=c(88899,125295))

## dataTypes parameter
bold_tax_id(id=88899, dataTypes = "basic")
bold_tax_id(id=88899, dataTypes = "stats")
bold_tax_id(id=88899, dataTypes = "images")
bold_tax_id(id=88899, dataTypes = "geo")
bold_tax_id(id=88899, dataTypes = "sequencinglabs")
bold_tax_id(id=88899, dataTypes = "depository")
bold_tax_id(id=c(88899,125295), dataTypes = "geo")
bold_tax_id(id=c(88899,125295), dataTypes = "images")

## Passing in NA
bold_tax_id(id = NA)
bold_tax_id(id = c(88899,125295,NA))

## get http response object only
bold_tax_id(id=88899, response=TRUE)
bold_tax_id(id=c(88899,125295), response=TRUE)

## curl debugging
bold_tax_id(id=88899, verbose = TRUE)
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=taxonomy
}
\seealso{
\code{\link[=bold_tax_name]{bold_tax_name()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_trace.R
\name{bold_trace}
\alias{bold_trace}
\alias{read_trace}
\title{Get BOLD trace files}
\usage{
bold_trace(
  taxon = NULL,
  ids = NULL,
  bin = NULL,
  container = NULL,
  institutions = NULL,
  researchers = NULL,
  geo = NULL,
  marker = NULL,
  dest = NULL,
  overwrite = TRUE,
  progress = TRUE,
  ...
)

read_trace(x)
}
\arguments{
\item{taxon}{(character) Returns all records containing matching taxa. Taxa
includes the ranks of  phylum, class, order, family, subfamily, genus,
and species.}

\item{ids}{(character) Returns all records containing matching IDs. IDs
include Sample IDs, Process IDs, Museum IDs and Field IDs.}

\item{bin}{(character) Returns all records contained in matching BINs. A
BIN is defined by a Barcode Index Number URI.}

\item{container}{(character) Returns all records contained in matching
projects or datasets. Containers include project codes and dataset codes}

\item{institutions}{(character) Returns all records stored in matching
institutions. Institutions are the Specimen Storing Site.}

\item{researchers}{(character) Returns all records containing matching
researcher names. Researchers include collectors and specimen identifiers.}

\item{geo}{(character) Returns all records collected in matching geographic
sites. Geographic sites includes countries and province/states.}

\item{marker}{(character) Returns all records containing matching
marker codes.}

\item{dest}{(character) A directory to write the files to}

\item{overwrite}{(logical) Overwrite existing directory and file?}

\item{progress}{(logical) Print progress or not. NOT AVAILABLE FOR NOW.
HOPEFULLY WILL RETURN SOON.}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}

\item{x}{Object to print or read.}
}
\description{
Get BOLD trace files
}
\examples{
\dontrun{
# Use a specific destination directory
bold_trace(taxon='Bombus', geo='Alaska', dest="~/mytarfiles")

# Another example
# bold_trace(ids='ACRJP618-11', dest="~/mytarfiles")
# bold_trace(ids=c('ACRJP618-11','ACRJP619-11'), dest="~/mytarfiles")

# read file in
x <- bold_trace(ids=c('ACRJP618-11','ACRJP619-11'), dest="~/mytarfiles")
(res <- read_trace(x$ab1[2]))

# The progress dialog is pretty verbose, so quiet=TRUE is a nice touch,
# but not by default
# Beware, this one take a while
# x <- bold_trace(taxon='Osmia', quiet=TRUE)

if (requireNamespace("sangerseqR", quietly = TRUE)) {
 library("sangerseqR")
 primarySeq(res)
 secondarySeq(res)
 head(traceMatrix(res))
}
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=webservices
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_specimens.R
\name{bold_specimens}
\alias{bold_specimens}
\title{Search BOLD for specimens.}
\usage{
bold_specimens(
  taxon = NULL,
  ids = NULL,
  bin = NULL,
  container = NULL,
  institutions = NULL,
  researchers = NULL,
  geo = NULL,
  response = FALSE,
  format = "tsv",
  ...
)
}
\arguments{
\item{taxon}{(character) Returns all records containing matching taxa. Taxa
includes the ranks of  phylum, class, order, family, subfamily, genus,
and species.}

\item{ids}{(character) Returns all records containing matching IDs. IDs
include Sample IDs, Process IDs, Museum IDs and Field IDs.}

\item{bin}{(character) Returns all records contained in matching BINs. A
BIN is defined by a Barcode Index Number URI.}

\item{container}{(character) Returns all records contained in matching
projects or datasets. Containers include project codes and dataset codes}

\item{institutions}{(character) Returns all records stored in matching
institutions. Institutions are the Specimen Storing Site.}

\item{researchers}{(character) Returns all records containing matching
researcher names. Researchers include collectors and specimen identifiers.}

\item{geo}{(character) Returns all records collected in matching geographic
sites. Geographic sites includes countries and province/states.}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{format}{(character) One of xml, json, tsv (default). tsv format gives
back a data.frame object. xml gives back parsed XML as \code{xml_document}
object. 'json' (JavaScript Object Notation) and 'dwc' (Darwin Core Archive)
are supported in theory, but the JSON can be malformed, so we don't support
that here, and the DWC option actually returns TSV.}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\description{
Search BOLD for specimens.
}
\examples{
\dontrun{
bold_specimens(taxon='Osmia')
bold_specimens(taxon='Osmia', format='xml')
bold_specimens(taxon='Osmia', response=TRUE)
res <- bold_specimens(taxon='Osmia', format='xml', response=TRUE)
res$url
res$status_code
res$response_headers

# More than 1 can be given for all search parameters
bold_specimens(taxon=c('Coelioxys','Osmia'))

## curl debugging
### These examples below take a long time, so you can set a timeout so that 
### it stops by X sec
head(bold_specimens(taxon='Osmia', verbose = TRUE))
# head(bold_specimens(geo='Costa Rica', timeout_ms = 6))
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=webservices
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_seqspec.R
\name{bold_seqspec}
\alias{bold_seqspec}
\title{Get BOLD specimen + sequence data.}
\usage{
bold_seqspec(
  taxon = NULL,
  ids = NULL,
  bin = NULL,
  container = NULL,
  institutions = NULL,
  researchers = NULL,
  geo = NULL,
  marker = NULL,
  response = FALSE,
  format = "tsv",
  sepfasta = FALSE,
  ...
)
}
\arguments{
\item{taxon}{(character) Returns all records containing matching taxa. Taxa
includes the ranks of  phylum, class, order, family, subfamily, genus,
and species.}

\item{ids}{(character) Returns all records containing matching IDs. IDs
include Sample IDs, Process IDs, Museum IDs and Field IDs.}

\item{bin}{(character) Returns all records contained in matching BINs. A
BIN is defined by a Barcode Index Number URI.}

\item{container}{(character) Returns all records contained in matching
projects or datasets. Containers include project codes and dataset codes}

\item{institutions}{(character) Returns all records stored in matching
institutions. Institutions are the Specimen Storing Site.}

\item{researchers}{(character) Returns all records containing matching
researcher names. Researchers include collectors and specimen identifiers.}

\item{geo}{(character) Returns all records collected in matching geographic
sites. Geographic sites includes countries and province/states.}

\item{marker}{(character) Returns all records containing matching marker
codes. See Details.}

\item{response}{(logical) Note that response is the object that returns
from the Curl call, useful for debugging, and getting detailed info on
the API call.}

\item{format}{(character) One of xml or tsv (default). tsv format gives
back a data.frame object. xml gives back parsed xml as a}

\item{sepfasta}{(logical) If \code{TRUE}, the fasta data is separated into
a list with names matching the processid's from the data frame.
Default: \code{FALSE}}

\item{...}{Further args passed on to
\code{\link[crul:verb-GET]{crul::verb-GET}}, main purpose being curl
debugging}
}
\value{
Either a data.frame, parsed xml, a http response object, or a list
with length two (a data.frame w/o nucleotide data, and a list with
nucleotide data)
}
\description{
Get BOLD specimen + sequence data.
}
\section{Large requests}{

Some requests can lead to errors. These often have to do with requesting
data for a rank that is quite high in the tree, such as an Order,
for example, Coleoptera. If your request is taking a long time,
it's likely that something will go wrong on the BOLD server side,
or we'll not be able to parse the result here in R because
R can only process strings of a certain length. \code{bold}
users have reported errors in which the resulting response from
BOLD is so large that we could not parse it.

A good strategy for when you want data for a high rank is to
do many separate requests for lower ranks within your target
rank. You can do this manually, or use the function
\code{taxize::downstream} to get all the names of a lower
rank within a target rank. There's an example in the README
(https://docs.ropensci.org/bold/#large-data)
}

\section{If a request times out}{

This is likely because you're request was for a large number of
sequences and the BOLD service timed out. You still should get
some output, those sequences that were retrieved before the time
out happened. As above, see the README
(https://docs.ropensci.org/bold/#large-data) for an example of
dealing with large data problems with this function.
}

\section{Marker}{

Notes from BOLD on the \code{marker} param:
"All markers for a specimen matching the search string will be returned.
ie. A record with COI-5P and ITS will return sequence data for both
markers even if only COI-5P was specified."

You will likely end up with data with markers that you did not request -
just be sure to filter those out as needed.
}

\examples{
\dontrun{
bold_seqspec(taxon='Osmia')
bold_seqspec(taxon='Osmia', format='xml')
bold_seqspec(taxon='Osmia', response=TRUE)
res <- bold_seqspec(taxon='Osmia', sepfasta=TRUE)
res$fasta[1:2]
res$fasta['GBAH0293-06']

# records that match a marker name
res <- bold_seqspec(taxon="Melanogrammus aeglefinus", marker="COI-5P")

# records that match a geographic locality
res <- bold_seqspec(taxon="Melanogrammus aeglefinus", geo="Canada")

## curl debugging
### You can do many things, including get verbose output on the curl call,
### and set a timeout
head(bold_seqspec(taxon='Osmia', verbose = TRUE))
## timeout
# head(bold_seqspec(taxon='Osmia', timeout_ms = 1))
}
}
\references{
http://v4.boldsystems.org/index.php/resources/api?type=webservices
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_filter.R
\name{bold_filter}
\alias{bold_filter}
\title{Filter BOLD specimen + sequence data (output of bold_seqspec)}
\usage{
bold_filter(x, by, how = "max")
}
\arguments{
\item{x}{(data.frame) a data.frame, as returned from
\code{\link[=bold_seqspec]{bold_seqspec()}}. Note that some combinations of parameters
in \code{\link[=bold_seqspec]{bold_seqspec()}} don't return a data.frame. Stops with
error message if this is not a data.frame. Required.}

\item{by}{(character) the column by which to group. For example,
if you want the longest sequence for each unique species name, then
pass \strong{species_name}. If the column doesn't exist, error
with message saying so. Required.}

\item{how}{(character) one of "max" or "min", which get used as
\code{which.max} or \code{which.min} to get the longest or shortest
sequence, respectively. Note that we remove gap/alignment characters
(\code{-})}
}
\value{
a tibble/data.frame
}
\description{
Picks either shortest or longest sequences, for a given grouping variable
(e.g., species name)
}
\examples{
\dontrun{
res <- bold_seqspec(taxon='Osmia')
maxx <- bold_filter(res, by = "species_name")
minn <- bold_filter(res, by = "species_name", how = "min")

vapply(maxx$nucleotides, nchar, 1, USE.NAMES = FALSE)
vapply(minn$nucleotides, nchar, 1, USE.NAMES = FALSE)
}
}
