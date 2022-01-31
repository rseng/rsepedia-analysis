traits
=======



[![cran checks](https://cranchecks.info/badges/worst/traits)](https://cranchecks.info/pkgs/traits)
[![Build Status](https://travis-ci.org/ropensci/traits.svg?branch=master)](https://travis-ci.org/ropensci/traits)
[![codecov](https://codecov.io/gh/ropensci/traits/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/traits)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/traits)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/traits)](https://CRAN.R-project.org/package=traits)

R client for various sources of species trait data.

Docs: https://docs.ropensci.org/traits/

What is a trait? A "trait" for the purposes of this package is broadly defined as an aspect of a species that can be described or measured, such as physical traits (size, length, height, color), behavioral traits (running speed, etc.), and even variables that make up the niche of the species (e.g., habitat).

Included in `traits` with the associated function prefix or function name:

<table>
<colgroup>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
</colgroup>

<thead>
<tr>
  <th style="text-align:left;">Souce</th>
  <th style="text-align:left;">Function prefix</th>
  <th style="text-align:left;">Link</th>
</tr>
</thead>

<tbody>
<tr>
  <td style="text-align:left;">BETYdb</td>
  <td style="text-align:left;"><code>betydb_</code></td>
  <td style="text-align:left;">https://www.betydb.org/</td>
</tr>
<tr>
  <td style="text-align:left;">NCBI</td>
  <td style="text-align:left;"><code>ncbi_</code></td>
  <td style="text-align:left;">https://www.ncbi.nlm.nih.gov/</td>
</tr>
<tr>
  <td style="text-align:left;">Encylopedia of Life</td>
  <td style="text-align:left;"><code>traitbank_</code></td>
  <td style="text-align:left;">https://github.com/EOL/eol_website/blob/master/doc/api.md</td>
</tr>
<tr>
  <td style="text-align:left;">Birdlife International</td>
  <td style="text-align:left;"><code>birdlife_</code></td>
  <td style="text-align:left;">https://www.birdlife.org/</td>
</tr>
<tr>
  <td style="text-align:left;">LEDA Traitbase</td>
  <td style="text-align:left;"><code>leda_</code></td>
  <td style="text-align:left;"></td>
</tr>
<tr>
  <td style="text-align:left;">Zanne et al. plant dataset</td>
  <td style="text-align:left;"><code>tr_zanne</code></td>
  <td style="text-align:left;"></td>
</tr>
<tr>
  <td style="text-align:left;">Amniote life history dataset</td>
  <td style="text-align:left;"><code>tr_ernest</code></td>
  <td style="text-align:left;"></td>
</tr>
</tbody>
</table>


Talk to us on the issues page (https://github.com/ropensci/traits/issues) if you know of a source of traits data with an API, and we'll see about including it.

## Installation

Stable CRAN version


```r
install.packages("traits")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/traits")
```


```r
library("traits")
library("dplyr")
```

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Zachary Foster](https://github.com/zachary-foster)
* [Ignasi Bartomeus](https://github.com/ibartomeus)
* [David LeBauer](https://github.com/dlebauer)
* [David Harris](https://github.com/davharris)
* [Chris Black](https://github.com/infotroph)
* [Rupert Collins](https://github.com/boopsboops)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/traits/issues).
* License: MIT
* Get citation information for `traits` in R doing `citation(package = 'traits')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
traits 0.5.0
============

### DEFUNCT

* `tr_usda()` is defunct. The API is down for good (#122)
* all `coral_*` functions are defunct (#124)

### MINOR IMPROVEMENTS

* replace httr with crul (#89)
* fix BETYdb tess (#123)

### BUG FIXES

* `ncbi_searcher()` fix: we weren't including the NCBI Entrez API key even when it was found  (#120)


traits 0.4.2
============

### MINOR IMPROVEMENTS

* betydb gains alias for API versions (#114)

### BUG FIXES

* `taxa_search`: removed `traitbank` option because `traits::traitbank()` used internally has completely changed and it's no longer feasible to do a straight-forward taxon search for all traits in EOL's Traitbank  (#115)


traits 0.4.0
============

### NEW FEATURES

* New package author: Chris Black (@infotroph) (#106) 
* betydb functions now can do pagination (#94)
* betydb functions gain progress parameter to optionally suppress the progress bar (#113)
* EOL Traitbank completely changed their query interface - function no longer works as it did before. for now, you have to specify your own query that's rather complex, see docs for help. Later on we can try to simplify queries for users (#112)

### MINOR IMPROVEMENTS

* table in README for different sources and clarify what traits are (#110) (#111)
* fixed link to Birdlife (#108)

### BUG FIXES

* fix to `ncbi_searcher()` to prevent failures in some cases (#107) thanks @zachary-foster
* fix to `ncbi_byid()`: ten new fields added to the output (#101) (#102) thanks @boopsboops


traits 0.3.0
============

### DEFUNCT

* Four functions are now defunct - those involving getting data
on whether a species if native/invasive in a particular region. 
See `?traits-defunct` for more information. Deprecated functions:
`eol_invasive_()`, `fe_native()`, `g_invasive()`, `is_native()` (#72)

### NEW FEATURES

* Gains new function `tr_ernest` for a dataset of Amniote life history
data (#60)
* Gains new function `tr_usda` for the USDA plants database (#61)
* Gains new function `tr_zanne` for a dataset of plant growth data (#73)
* BetyDB functions gain automatic paging of large requests where API supports it, i.e. not in v0 (#94)

### MINOR IMPROVEMENTS

* Change Coral database base URL to https (#99)
* Now requiring `readr > 1.0` (#76)
* Changed `ncbi_*()` functions to give back `NA` types that match
data.frame column classes to make combining easier (#96)
* replace `xml2::xml_find_one` with `xml2::xml_find_first` throughout (#97)
* namespace all fxn calls for base pkgs, remove from Imports (#98)
* BetyDB cleanup (#25) (#77) (#82) (#88) 

### BUG FIXES

* Fixed `birdlife*` functions that needed to change URL structure 
due to changes in the Birdlife website (#100)
* Fixes to `traitbank()` (#79) (#80) thanks @dschlaep !
* `ncbi_*()` fxns now use https (#95)


traits 0.2.0
============

### DEPRECATED

* Marked four functions as deprecated - those involving getting data
on whether a species if native/invasive in a particular region. 
See `?traits-deprecated` for more information. Deprecated functions:
`eol_invasive_()`, `fe_native()`, `g_invasive()`, `is_native()` (#63)

### MINOR IMPROVEMENTS

* Standardized outputs of all data - all data.frame column names should 
be lowercase now (#47)
* With all `httr::content()` calls now explicitly setting encoding to 
`UTF-8`, and parsing to `text`, then manually parsing either JSON
or XML later (#65)
* Replaced `XML` with `xml2` for XML parsing (#67)

traits 0.1.2
============

### NEW FEATURES

* `ncbi_searcher()` gains new parameter `fuzzy` to toggle fuzzy taxonomic ID search or exact search. (#34) (thx @mpnelsen)

### MINOR IMPROVEMENTS

* Importing only functions (via `importFrom`) used across all imports now.
In addition, `importFrom` for all non-base R pkgs, including `methods`,
`stats` and `utils` packages (#36)
* Changed the `trait` parameter in `traitbank()` function to `pageid`, 
because EOL expects a page identifier, which is associated with a taxon, 
not a trait. The previous parameter name was very misleading.

traits 0.1.0
============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.2 patched
* ubuntu 14.04 (on travis-ci), R 4.0.2 
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

------

This version makes some functions defunct, and fixes some bugs.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/traits/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/traits.git`
* Make sure to track progress upstream (i.e., on our version of `traits` at `ropensci/traits`) by doing `git remote add upstream https://github.com/ropensci/traits.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/traits`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-08-25                                  |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|traits  |0.4.2 |0.5.0      |*  |
|crayon  |NA    |1.3.4.9000 |*  |

# Revdeps

## Failed to check (1)

|package   |version |error |warning |note |
|:---------|:-------|:-----|:-------|:----|
|metacoder |0.3.4   |1     |        |     |

*Wow, no problems at all. :)*# metacoder

<details>

* Version: 0.3.4
* Source code: https://github.com/cran/metacoder
* URL: https://grunwaldlab.github.io/metacoder_documentation/
* BugReports: https://github.com/grunwaldlab/metacoder/issues
* Date/Publication: 2020-04-29 19:40:03 UTC
* Number of recursive dependencies: 150

Run `revdep_details(,"metacoder")` for more info

</details>

## In both

*   checking whether package ‘metacoder’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/sckott/github/ropensci/traits/revdep/checks.noindex/metacoder/new/metacoder.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘metacoder’ ...
** package ‘metacoder’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/sckott/github/ropensci/traits/revdep/library.noindex/metacoder/Rcpp/include' -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c RcppExports.cpp -o RcppExports.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/sckott/github/ropensci/traits/revdep/library.noindex/metacoder/Rcpp/include' -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c repel_boxes.cpp -o repel_boxes.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -Wl,-rpath,/Library/Frameworks/R.framework/Resources/lib /Library/Frameworks/R.framework/Resources/lib/libc++abi.1.dylib -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o metacoder.so RcppExports.o repel_boxes.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
clang-7: error: no such file or directory: '/Library/Frameworks/R.framework/Resources/lib/libc++abi.1.dylib'
make: *** [metacoder.so] Error 1
ERROR: compilation failed for package ‘metacoder’
* removing ‘/Users/sckott/github/ropensci/traits/revdep/checks.noindex/metacoder/new/metacoder.Rcheck/metacoder’

```
### CRAN

```
* installing *source* package ‘metacoder’ ...
** package ‘metacoder’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/sckott/github/ropensci/traits/revdep/library.noindex/metacoder/Rcpp/include' -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c RcppExports.cpp -o RcppExports.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/sckott/github/ropensci/traits/revdep/library.noindex/metacoder/Rcpp/include' -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c repel_boxes.cpp -o repel_boxes.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -Wl,-rpath,/Library/Frameworks/R.framework/Resources/lib /Library/Frameworks/R.framework/Resources/lib/libc++abi.1.dylib -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o metacoder.so RcppExports.o repel_boxes.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
clang-7: error: no such file or directory: '/Library/Frameworks/R.framework/Resources/lib/libc++abi.1.dylib'
make: *** [metacoder.so] Error 1
ERROR: compilation failed for package ‘metacoder’
* removing ‘/Users/sckott/github/ropensci/traits/revdep/checks.noindex/metacoder/old/metacoder.Rcheck/metacoder’

```
traits
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![cran checks](https://cranchecks.info/badges/worst/traits)](https://cranchecks.info/pkgs/traits)
[![Build Status](https://travis-ci.org/ropensci/traits.svg?branch=master)](https://travis-ci.org/ropensci/traits)
[![codecov](https://codecov.io/gh/ropensci/traits/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/traits)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/traits)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/traits)](https://CRAN.R-project.org/package=traits)

R client for various sources of species trait data.

Docs: https://docs.ropensci.org/traits/

What is a trait? A "trait" for the purposes of this package is broadly defined as an aspect of a species that can be described or measured, such as physical traits (size, length, height, color), behavioral traits (running speed, etc.), and even variables that make up the niche of the species (e.g., habitat).

Included in `traits` with the associated function prefix or function name:

<table>
<colgroup>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
<col style="text-align:left;"/>
</colgroup>

<thead>
<tr>
  <th style="text-align:left;">Souce</th>
  <th style="text-align:left;">Function prefix</th>
  <th style="text-align:left;">Link</th>
</tr>
</thead>

<tbody>
<tr>
  <td style="text-align:left;">BETYdb</td>
  <td style="text-align:left;"><code>betydb_</code></td>
  <td style="text-align:left;">https://www.betydb.org/</td>
</tr>
<tr>
  <td style="text-align:left;">NCBI</td>
  <td style="text-align:left;"><code>ncbi_</code></td>
  <td style="text-align:left;">https://www.ncbi.nlm.nih.gov/</td>
</tr>
<tr>
  <td style="text-align:left;">Encylopedia of Life</td>
  <td style="text-align:left;"><code>traitbank_</code></td>
  <td style="text-align:left;">https://github.com/EOL/eol_website/blob/master/doc/api.md</td>
</tr>
<tr>
  <td style="text-align:left;">Birdlife International</td>
  <td style="text-align:left;"><code>birdlife_</code></td>
  <td style="text-align:left;">https://www.birdlife.org/</td>
</tr>
<tr>
  <td style="text-align:left;">LEDA Traitbase</td>
  <td style="text-align:left;"><code>leda_</code></td>
  <td style="text-align:left;"></td>
</tr>
<tr>
  <td style="text-align:left;">Zanne et al. plant dataset</td>
  <td style="text-align:left;"><code>tr_zanne</code></td>
  <td style="text-align:left;"></td>
</tr>
<tr>
  <td style="text-align:left;">Amniote life history dataset</td>
  <td style="text-align:left;"><code>tr_ernest</code></td>
  <td style="text-align:left;"></td>
</tr>
</tbody>
</table>


Talk to us on the issues page (https://github.com/ropensci/traits/issues) if you know of a source of traits data with an API, and we'll see about including it.

## Installation

Stable CRAN version

```{r eval=FALSE}
install.packages("traits")
```

Or development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/traits")
```

```{r}
library("traits")
library("dplyr")
```

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Zachary Foster](https://github.com/zachary-foster)
* [Ignasi Bartomeus](https://github.com/ibartomeus)
* [David LeBauer](https://github.com/dlebauer)
* [David Harris](https://github.com/davharris)
* [Chris Black](https://github.com/infotroph)
* [Rupert Collins](https://github.com/boopsboops)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/traits/issues).
* License: MIT
* Get citation information for `traits` in R doing `citation(package = 'traits')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: traits introduction
author: Scott Chamberlain
date: "2020-08-25"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{traits introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---





```r
library("traits")
```

## BetyDB

Get trait data for Willow (_Salix_ spp.)


```r
(salix <- betydb_search("Salix Vcmax"))
#> # A tibble: 14 x 36
#>    checked result_type    id citation_id site_id treatment_id sitename city 
#>      <int> <chr>       <int>       <int>   <int>        <int> <chr>    <chr>
#>  1       1 traits      39217         430     645         1342 ""       Saare
#>  2       1 traits      39218         430     645         1343 ""       Saare
#>  3       1 traits      39219         430     645         1344 ""       Saare
#>  4       1 traits      39220         430     645         1345 ""       Saare
#>  5       1 traits      25405          51      NA            1  <NA>    <NA> 
#>  6       1 traits      39213         430     645         1342 ""       Saare
#>  7       1 traits      39214         430     645         1343 ""       Saare
#>  8       1 traits      39215         430     645         1344 ""       Saare
#>  9       1 traits      39216         430     645         1345 ""       Saare
#> 10       1 traits      39221         430     645         1342 ""       Saare
#> 11       1 traits      39222         430     645         1343 ""       Saare
#> 12       1 traits      39223         430     645         1344 ""       Saare
#> 13       1 traits      39224         430     645         1345 ""       Saare
#> 14       1 traits      37519         381     602         1220  <NA>    <NA> 
#> # … with 28 more variables: lat <dbl>, lon <dbl>, scientificname <chr>,
#> #   commonname <chr>, genus <chr>, species_id <int>, cultivar_id <int>,
#> #   author <chr>, citation_year <int>, treatment <chr>, date <chr>, time <chr>,
#> #   raw_date <chr>, month <int>, year <int>, dateloc <chr>, trait <chr>,
#> #   trait_description <chr>, mean <dbl>, units <chr>, n <int>, statname <chr>,
#> #   stat <dbl>, notes <chr>, access_level <int>, cultivar <chr>, entity <lgl>,
#> #   method_name <lgl>
# equivalent:
# (out <- betydb_search("willow"))
```

Summarise data from the output `data.frame`


```r
library("dplyr")
salix %>%
  group_by(scientificname, trait) %>%
  mutate(.mean = as.numeric(mean)) %>%
  summarise(mean = round(mean(.mean, na.rm = TRUE), 2),
            min = round(min(.mean, na.rm = TRUE), 2),
            max = round(max(.mean, na.rm = TRUE), 2),
            n = length(n))
#> # A tibble: 4 x 6
#> # Groups:   scientificname [4]
#>   scientificname                  trait  mean   min   max     n
#>   <chr>                           <chr> <dbl> <dbl> <dbl> <int>
#> 1 Salix                           Vcmax  65    65    65       1
#> 2 Salix dasyclados                Vcmax  46.1  34.3  56.7     4
#> 3 Salix sachalinensis × miyabeana Vcmax  79.3  79.3  79.3     1
#> 4 Salix viminalis                 Vcmax  43.0  20.0  61.3     8
```

## NCBI sequence data

Get sequences by id


```r
ncbi_byid(ids = "360040093")
#>                  taxon
#> 1 Eristalis transversa
#>                                                                                                                                                                                       taxonomy
#> 1 Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Holometabola; Diptera; Brachycera; Muscomorpha; Syrphoidea; Syrphidae; Eristalinae; Eristalini; Eristalis
#>                                                                                                             gene_desc
#> 1 Eristalis transversa voucher CNC:Diptera:102013 cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial
#>       organelle     gi_no     acc_no keyword   specimen_voucher
#> 1 mitochondrion 360040093 JN991986.1 BARCODE CNC:Diptera:102013
#>               lat_lon
#> 1 38.4623 N 79.2417 W
#>                                                       country
#> 1 USA: Virginia, Reddish Knob Lookout, 14.5km W Briery Branch
#>                                                              paper_title
#> 1 The evolution of imperfect mimicry in hover flies (Diptera: Syrphidae)
#>       journal first_author uploaded_date length
#> 1 Unpublished   Penny,H.D.   03-NOV-2012    658
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             sequence
#> 1 tactttatattttgtatttggaacatgagcgggtatagtaggaacttcattaagaattttaattcgagctgaattaggtcatccaggtgcattaattggtgatgatcaaatttataatgttattgtaacagctcatgcttttgttataattttttttatagtaatacctattataattggaggatttggaaattgattagtaccacttatattaggagctccagatatagcattccctcgaataaataatataagtttctgattattacctccttctttaactctattattagtaagaagtatagtagaaaatggggctggaacaggatgaacagtttatcctccattatcaagtaatattgcacatggaggagcctcagttgatttagcaattttttcacttcacttatcaggaatatcatctattttaggtgcagtaaattttattacaacagttattaatatacgatcaacaggaattacttatgatcgtatacctttatttgtttgatctgttgctattacagctttattattattattatcattaccagtactagcaggagctattacaatattattaactgatcgaaatttaaatacatcattctttgatccagcaggaggaggagaccctatcctgtaccaacacttattc
```

Get sequences searching by taxonomic name


```r
out <- ncbi_searcher(taxa = "Umbra limi", seqrange = "1:2000")
#> using sleep: 0
#> ══  1 queries  ═══════════════
#> ✔  Found:  Umbra+limi
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
head(out)
#>        taxon length
#> 1 Umbra limi    298
#> 2 Umbra limi    761
#> 3 Umbra limi    765
#> 4 Umbra limi    764
#> 5 Umbra limi    743
#> 6 Umbra limi    758
#>                                                                                                gene_desc
#> 1 Umbra limi voucher Umbra_limi_GLF16S large subunit ribosomal RNA gene, partial sequence; mitochondrial
#> 2                                        Umbra limi voucher NXG2012264 rhodopsin (Rho) gene, partial cds
#> 3                                         Umbra limi voucher NXG201250 rhodopsin (Rho) gene, partial cds
#> 4                                        Umbra limi voucher NXG2012183 rhodopsin (Rho) gene, partial cds
#> 5                                         Umbra limi voucher NXG201252 rhodopsin (Rho) gene, partial cds
#> 6                                        Umbra limi voucher NXG2012231 rhodopsin (Rho) gene, partial cds
#>     acc_no      gi_no
#> 1 MT549086 1847697133
#> 2 KX146134 1049488959
#> 3 KX146015 1049488721
#> 4 KX145969 1049488629
#> 5 KX145777 1049488245
#> 6 KX145759 1049488209
```

## EOL's traitbank trait data


```r
traitbank(query = "MATCH (n:Trait) RETURN n LIMIT 1;")
#> $columns
#> [1] "n"
#> 
#> $data
#> $data[[1]]
#>   metadata.id metadata.labels    data.eol_pk data.object_page_id
#> 1    22529388           Trait R20-PK20910350            46581789
#>                                                    data.resource_pk
#> 1 ReverseOf_globi:assoc:7296029-FBC:FB:SpecCode:4755-ATE-EOL_V2:281
#>   data.scientific_name
#> 1              Plantae
#>                                                                                                                                                                                                                                                             data.source
#> 1 Froese, R. and D. Pauly. Editors. 2019. FishBase. World Wide Web electronic publication. www.fishbase.org, version (08/2019). Accessed at <https://github.com/globalbioticinteractions/fishbase/archive/6ebceaacea18c6ff6c247182f9af8ad6fc05cc82.zip> on 25 May 2020.
```

## Birdlife International

Habitat data


```r
birdlife_habitat(22721692)
#>         id Habitat (level 1)                  Habitat (level 2) Importance
#> 1 22721692            Forest           Subtropical/Tropical Dry      major
#> 2 22721692            Forest Subtropical/Tropical Moist Montane      major
#> 3 22721692            Forest                          Temperate   suitable
#> 4 22721692         Shrubland Subtropical/Tropical High Altitude   suitable
#>     Occurrence
#> 1     breeding
#> 2 non-breeding
#> 3     breeding
#> 4     breeding
```

Threats data


```r
birdlife_threats(22721692)
#>          id                                                  threat1
#> 1  22721692                                Agriculture & aquaculture
#> 2  22721692                                Agriculture & aquaculture
#> 3  22721692                                  Biological resource use
#> 4  22721692                          Climate change & severe weather
#> 5  22721692                          Climate change & severe weather
#> 6  22721692                          Climate change & severe weather
#> 7  22721692 Invasive and other problematic species, genes & diseases
#> 8  22721692 Invasive and other problematic species, genes & diseases
#> 9  22721692 Invasive and other problematic species, genes & diseases
#> 10 22721692 Invasive and other problematic species, genes & diseases
#> 11 22721692 Invasive and other problematic species, genes & diseases
#> 12 22721692 Invasive and other problematic species, genes & diseases
#> 13 22721692                             Natural system modifications
#> 14 22721692                             Natural system modifications
#> 15 22721692                     Residential & commercial development
#> 16 22721692                       Transportation & service corridors
#>                                                                                   threat2
#> 1                             Annual & perennial non-timber crops - Agro-industry farming
#> 2                              Annual & perennial non-timber crops - Small-holder farming
#> 3  Logging & wood harvesting - Unintentional effects: (subsistence/small scale) [harvest]
#> 4                                                                                Droughts
#> 5                                                           Habitat shifting & alteration
#> 6                                                                       Storms & flooding
#> 7                                 Invasive non-native/alien species/diseases - Sus scrofa
#> 8                        Invasive non-native/alien species/diseases - Unspecified species
#> 9                            Problematic native species/diseases - Dendroctonus frontalis
#> 10                           Problematic native species/diseases - Odocoileus virginianus
#> 11                              Problematic native species/diseases - Unspecified species
#> 12                    Problematic species/disease of unknown origin - Unspecified species
#> 13                                     Fire & fire suppression - Trend Unknown/Unrecorded
#> 14                                                          Other ecosystem modifications
#> 15                                                                  Housing & urban areas
#> 16                                                                      Roads & railroads
#>                                                                 stresses
#> 1                     Ecosystem degradation, Ecosystem conversion, Other
#> 2                     Ecosystem degradation, Ecosystem conversion, Other
#> 3                                                  Ecosystem degradation
#> 4                                                  Ecosystem degradation
#> 5                            Ecosystem degradation, Ecosystem conversion
#> 6                                                  Ecosystem degradation
#> 7                            Ecosystem degradation, Ecosystem conversion
#> 8                                                      Species mortality
#> 9                                                  Ecosystem degradation
#> 10                                                 Ecosystem degradation
#> 11                   Ecosystem degradation, Reduced reproductive success
#> 12                                                     Species mortality
#> 13                           Ecosystem degradation, Ecosystem conversion
#> 14                           Ecosystem degradation, Ecosystem conversion
#> 15 Ecosystem degradation, Ecosystem conversion, Species mortality, Other
#> 16                                                     Species mortality
#>                                                      timing
#> 1                                 Agriculture & aquaculture
#> 2                                 Agriculture & aquaculture
#> 3                                   Biological resource use
#> 4                           Climate change & severe weather
#> 5                           Climate change & severe weather
#> 6                           Climate change & severe weather
#> 7  Invasive and other problematic species, genes & diseases
#> 8  Invasive and other problematic species, genes & diseases
#> 9  Invasive and other problematic species, genes & diseases
#> 10 Invasive and other problematic species, genes & diseases
#> 11 Invasive and other problematic species, genes & diseases
#> 12 Invasive and other problematic species, genes & diseases
#> 13                             Natural system modifications
#> 14                             Natural system modifications
#> 15                     Residential & commercial development
#> 16                       Transportation & service corridors
#>                                                                                     scope
#> 1                             Annual & perennial non-timber crops - Agro-industry farming
#> 2                              Annual & perennial non-timber crops - Small-holder farming
#> 3  Logging & wood harvesting - Unintentional effects: (subsistence/small scale) [harvest]
#> 4                                                                                Droughts
#> 5                                                           Habitat shifting & alteration
#> 6                                                                       Storms & flooding
#> 7                                 Invasive non-native/alien species/diseases - Sus scrofa
#> 8                        Invasive non-native/alien species/diseases - Unspecified species
#> 9                            Problematic native species/diseases - Dendroctonus frontalis
#> 10                           Problematic native species/diseases - Odocoileus virginianus
#> 11                              Problematic native species/diseases - Unspecified species
#> 12                    Problematic species/disease of unknown origin - Unspecified species
#> 13                                     Fire & fire suppression - Trend Unknown/Unrecorded
#> 14                                                          Other ecosystem modifications
#> 15                                                                  Housing & urban areas
#> 16                                                                      Roads & railroads
#>    severity  impact
#> 1   Ongoing Ongoing
#> 2   Ongoing Ongoing
#> 3   Ongoing Ongoing
#> 4   Ongoing Ongoing
#> 5    Future  Future
#> 6   Ongoing Ongoing
#> 7   Ongoing Ongoing
#> 8   Ongoing Ongoing
#> 9   Ongoing Ongoing
#> 10  Ongoing Ongoing
#> 11  Ongoing Ongoing
#> 12  Ongoing Ongoing
#> 13  Ongoing Ongoing
#> 14   Future  Future
#> 15  Ongoing Ongoing
#> 16  Ongoing Ongoing
```
---
title: BETYdb Tutorial
author: Scott Chamberlain
date: "2020-08-25"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{BETYdb Tutorial}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



[BETYdb](https://www.betydb.org/) is the _Biofuel Ecophysiological Traits and Yields Database_. You can get many different types of data from this database, including trait data. 

Function setup: All functions are prefixed with `betydb_`. Plural function names like `betydb_traits()` accept parameters and always give back a data.frame, while singlur function names like `betydb_trait()` accept an ID and give back a list. 

The idea with the functions with plural names is to search for either traits, species, etc., and with the singular function names to get data by one or more IDs.

## Load traits


```r
library("traits")
```

## Traits

Get trait data for _Miscanthus giganteus_


```r
out <- betydb_search(query = "Switchgrass Yield")
```

Summarise data from the output `data.frame`


```r
suppressPackageStartupMessages(library("dplyr"))
out %>%
  group_by(id) %>%
  summarise(mean_result = mean(as.numeric(mean), na.rm = TRUE)) %>%
  arrange(desc(mean_result))
```

```
#> # A tibble: 449 x 2
#>       id mean_result
#>    <int>       <dbl>
#>  1  1666        27.4
#>  2 16845        27  
#>  3  1669        26.4
#>  4 16518        26  
#>  5  1663        25.4
#>  6 16742        25  
#>  7  1594        24.8
#>  8  1674        22.7
#>  9  1606        22.5
#> 10  1665        22.5
#> # … with 439 more rows
```

Single trait


```r
betydb_trait(id = 10)
```

```
#> # A tibble: 1 x 13
#>      id description units notes created_at updated_at name  max   min  
#>   <int> <chr>       <chr> <chr> <chr>      <chr>      <chr> <chr> <chr>
#> 1    10 Leaf Perce… perc… <NA>  <NA>       2011-06-0… leafN 10    0.02 
#> # … with 4 more variables: standard_name <chr>, standard_units <chr>,
#> #   label <chr>, type <chr>
```

## Species

Single species, _Acacia karroothorn_


```r
betydb_specie(id = 10)
```

```
#> # A tibble: 1 x 10
#>      id spcd  genus species scientificname commonname notes created_at
#>   <int> <chr> <chr> <chr>   <chr>          <chr>      <chr> <chr>     
#> 1    10 <NA>  Acac… karroo  Acacia karroo  karrootho… <NA>  <NA>      
#> # … with 2 more variables: updated_at <chr>, acceptedsymbol <chr>
```

## Citations

Get citatons searching for _Miscanthus_


```r
betydb_citation(10)
```

```
#> # A tibble: 1 x 13
#>      id author  year title journal   vol pg    url   pdf   created_at updated_at
#>   <int> <chr>  <int> <chr> <chr>   <int> <chr> <chr> <chr> <chr>      <chr>     
#> 1    10 Casler  2003 Cult… Crop S…    43 2226… http… http… <NA>       <NA>      
#> # … with 2 more variables: doi <chr>, user_id <chr>
```

## Sites

Single site


```r
betydb_site(id = 1)
```

```
#> # A tibble: 1 x 8
#>   city    state country notes sitename greenhouse geometry           time_zone  
#>   <chr>   <chr> <chr>   <chr> <chr>    <lgl>      <chr>              <chr>      
#> 1 Aliart… <NA>  GR      <NA>  Aliartos FALSE      POINT (23.17 38.3… Europe/Ath…
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eol_invasive.R
\name{eol_invasive_}
\alias{eol_invasive_}
\title{Search for presence of taxonomic names in EOL invasive species databases}
\usage{
eol_invasive_(...)
}
\description{
Search for presence of taxonomic names in EOL invasive species databases
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betydb.R
\name{betydb}
\alias{betydb}
\alias{betydb_record}
\alias{betydb_trait}
\alias{betydb_specie}
\alias{betydb_citation}
\alias{betydb_site}
\alias{betydb_experiment}
\title{Search for traits from BETYdb}
\usage{
betydb_record(
  id,
  table,
  api_version = NULL,
  betyurl = NULL,
  fmt = NULL,
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)

betydb_trait(
  id,
  genus = NULL,
  species = NULL,
  api_version = NULL,
  betyurl = NULL,
  fmt = "json",
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)

betydb_specie(
  id,
  genus = NULL,
  species = NULL,
  api_version = NULL,
  betyurl = NULL,
  fmt = "json",
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)

betydb_citation(
  id,
  genus = NULL,
  species = NULL,
  api_version = NULL,
  betyurl = NULL,
  fmt = "json",
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)

betydb_site(
  id,
  api_version = NULL,
  betyurl = NULL,
  fmt = "json",
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)

betydb_experiment(
  id,
  api_version = NULL,
  betyurl = NULL,
  fmt = "json",
  key = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{id}{(integer) One or more ids for a species, site, variable, etc.}

\item{table}{(character) Name of the database table with which this ID is associated.}

\item{api_version}{(character) Which version of the BETY API should we query? One of "v0" or "beta". Default is \code{options("betydb_api_version")} if set, otherwise  "v0".}

\item{betyurl}{(string) url to target instance of betydb. Default is \code{options("betydb_url")} if set, otherwise "https:/www.betydb.org/"}

\item{fmt}{(character) Format to return data in, one of json, xml, csv. Only json
currently supported.}

\item{key}{(character) An API key. Use this or user/pwd combo. Save in your
\code{.Rprofile} file as \code{options(betydb_key = "your40digitkey")}. Optional}

\item{user, pwd}{(character) A user name and password. Use a user/pwd combo or an API key.
Save in your \code{.Rprofile} file as \code{options(betydb_user = "yournamehere")} and \code{options(betydb_pwd = "yourpasswordhere")}. Optional}

\item{progress}{show progress bar? default: \code{TRUE}}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}. Optional}

\item{genus}{(character) A genus name. Optional}

\item{species}{(character) A specific epithet. Optional}
}
\description{
Search for traits from BETYdb

Get details about a single item from a table
}
\details{
BETYdb includes a primary home page (betydb.org) focused on bioenergy crops as well as a network of harmonized
databases that support and share data among more focused research programs.

For a list of publicly accessible instances of BETYdb and the urls that can be queried,
see \url{https://pecan.gitbooks.io/betydb-documentation/content/distributed_betydb.html}

This package queries plant traits, phenotypes, biomass yields, and ecosystem functions.
It does not currently interface with the workflow and provenance data that support PEcAn Project (pecanproject.org) and TERRA REF (terraref.org) software.

API documentation: \url{https://pecan.gitbooks.io/betydb-data-access/content/API.html}
API endpoints are here: \url{https://www.betydb.org/api/docs}
This package currently uses the original 'v0' API by default.
To use a newer version, set \code{api_version}.
Newer versions of the API will support database inserts.
}
\section{Authentication}{

Defers to use API key first since it's simpler, but if you don't have
an API key, you can supply a username and password.
}

\section{Functions}{

Singular functions like \code{betydb_trait} accept an id and additional parameters,
and return a list of variable outputs depending on the inputs.

However, plural functions like \code{betydb_traits} accept query parameters, but not
ids, and always return a single data.frame.

\code{betydb_search("Search terms", ...)} is a convenience wrapper that passes all further arguments to \code{\link{betydb_query}(table = "search", search = "Search terms", ...)}. See there for details on possible arguments.
}

\examples{
\dontrun{
# General Search
out <- betydb_search(query = "Switchgrass Yield")
library("dplyr")
out \%>\%
 group_by(id) \%>\%
 summarise(mean_result = mean(as.numeric(mean), na.rm = TRUE)) \%>\%
 arrange(desc(mean_result))
# Get by ID
## Traits
betydb_trait(id = 10)
## Species
betydb_specie(id = 1)
## Citations
betydb_citation(id = 1)
## Site information
betydb_site(id = 795)
}
}
\references{
API documentation \url{https://pecan.gitbooks.io/betydb-data-access/content/API.html} and
https://www.betydb.org/api/docs
}
\seealso{
\code{\link{betydb_query}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/birdlife.R
\name{birdlife_habitat}
\alias{birdlife_habitat}
\title{Get bird habitat information from BirdLife/IUCN}
\usage{
birdlife_habitat(id)
}
\arguments{
\item{id}{A single IUCN species ID}
}
\value{
a \code{data.frame} with level 1 and level 2 habitat classes, as
well as importance ratings and occurrence type (e.g. breeding or
non-breeding). The habitat classification scheme is described
at https://www.iucnredlist.org/resources/classification-schemes
}
\description{
Get bird habitat information from BirdLife/IUCN
}
\examples{
\dontrun{
# Setophaga chrysoparia
birdlife_habitat(22721692)
# Passer domesticus
birdlife_habitat(103818789)
}
}
\seealso{
Other birdlife: 
\code{\link{birdlife_threats}()}
}
\author{
David J. Harris \email{harry491@gmail.com}
}
\concept{birdlife}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leda.R
\name{leda}
\alias{leda}
\title{Access LEDA trait data}
\usage{
leda(trait = "age_first_flowering", ...)
}
\arguments{
\item{trait}{(character) Trait to get. See Details.}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Access LEDA trait data
}
\details{
For parameter \code{trait}, one of age_first_flowering, branching,
buds_seasonality, buds_vertical_dist, canopy_height, dispersal_type,
leaf_distribution, ldmc_geo, leaf_mass, leaf_size, morphology_disperal,
growth_form, life_span, releasing_height, seed_longevity,
seed_mass, seed_number, seed_shape, shoot_growth_form, snp, ssd, tv,
or clonal_growth_organs

The following are not supported as they are too much of a pain to parse:
buoyancy, seed_bank, sla_geo
}
\examples{
\dontrun{
# Age of first flowering
leda(trait = "age_first_flowering")

# Seed number
leda("seed_number")

# Releasing height
leda(trait = "releasing_height")

# Clonal growth organs
leda(trait = "clonal_growth_organs")

all <- c("age_first_flowering", "branching", "buds_seasonality",
  "buds_vertical_dist", "canopy_height",
  "dispersal_type", "leaf_distribution", "ldmc_geo", "leaf_mass",
  "leaf_size", "morphology_disperal", "growth_form", "life_span",
  "releasing_height", "seed_longevity", "seed_mass",
  "seed_number", "seed_shape", "shoot_growth_form",
  "snp", "ssd", "tv", "clonal_growth_organs")
out <- list()
for (i in seq_along(all)) {
  cat(all[i], sep="\n")
  out[[i]] <- leda(all[i])
}
sapply(out, NROW)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tr_zanne.R
\name{tr_zanne}
\alias{tr_zanne}
\title{Zanne et al. plant dataset}
\usage{
tr_zanne(read = TRUE, ...)
}
\arguments{
\item{read}{(logical) read in csv files. Default: \code{TRUE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
paths to the files (character) if \code{read=FALSE} or
a list of data.frame's if \code{read=TRUE}
}
\description{
Zanne et al. plant dataset
}
\details{
This data is a dataset stored on Dryad (doi: 10.5061/dryad.63q27).
When using this data, cite the paper:

Zanne AE, Tank DC, Cornwell WK, Eastman JM, Smith SA, FitzJohn RG,
McGlinn DJ, O'Meara BC, Moles AT, Reich PB, Royer DL, Soltis DE, Stevens PF,
Westoby M, Wright IJ, Aarssen L, Bertin RI, Calaminus A, Govaerts R,
Hemmings F, Leishman MR, Oleksyn J, Soltis PS, Swenson NG, Warman L,
Beaulieu JM, Ordonez A (2014) Three keys to the radiation of angiosperms
into freezing environments. Nature 506(7486): 89-92.
http://dx.doi.org/10.1038/nature12872

As well as the Dryad data package:

Zanne AE, Tank DC, Cornwell WK, Eastman JM, Smith SA, FitzJohn RG,
McGlinn DJ, O'Meara BC, Moles AT, Reich PB, Royer DL, Soltis DE, Stevens PF,
Westoby M, Wright IJ, Aarssen L, Bertin RI, Calaminus A, Govaerts R,
Hemmings F, Leishman MR, Oleksyn J, Soltis PS, Swenson NG, Warman L,
Beaulieu JM, Ordonez A (2013) Data from: Three keys to the radiation of
angiosperms into freezing environments. Dryad Digital Repository.
http://dx.doi.org/10.5061/dryad.63q27.2
}
\examples{
\dontrun{
res <- tr_zanne()
res$tax_lookup
res$woodiness
res$freezing
res$leaf_phenology
}
}
\references{
http://datadryad.org/resource/doi:10.5061/dryad.63q27
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_byname.R
\name{ncbi_byname}
\alias{ncbi_byname}
\title{Retrieve gene sequences from NCBI by taxon name and gene names.}
\usage{
ncbi_byname(
  taxa,
  gene = "COI",
  seqrange = "1:3000",
  getrelated = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{taxa}{(character) Scientific name to search for.}

\item{gene}{(character) Gene or genes (in a vector) to search for.
See examples.}

\item{seqrange}{(character) Sequence range, as e.g., \code{"1:1000"}. This is the range of
sequence lengths to search for. So \code{"1:1000"} means search for sequences from 1 to 1000
characters in length.}

\item{getrelated}{(logical) If \code{TRUE}, gets the longest sequences of a species
in the same genus as the one searched for. If \code{FALSE}, returns nothing if no match
found.}

\item{verbose}{(logical) If \code{TRUE} (default), informative messages printed.}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
data.frame
}
\description{
Retrieve gene sequences from NCBI by taxon name and gene names.
}
\details{
Removes predicted sequences so you don't have to remove them.
Predicted sequences are those with accession numbers that have "XM_" or
"XR_" prefixes. This function retrieves one sequences for each species,
picking the longest available for the given gene.
}
\examples{
\dontrun{
# A single species
ncbi_byname(taxa="Acipenser brevirostrum")

# Many species
species <- c("Colletes similis","Halictus ligatus","Perdita californica")
ncbi_byname(taxa=species, gene = c("coi", "co1"), seqrange = "1:2000")
}
}
\seealso{
\code{\link[=ncbi_searcher]{ncbi_searcher()}}, \code{\link[=ncbi_byid]{ncbi_byid()}}
}
\author{
Scott Chamberlain
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traits-package.r
\docType{package}
\name{traits-package}
\alias{traits-package}
\alias{traits}
\title{traits - Species trait data from around the web}
\description{
Currently included in \code{traits} with the associated function name or
function prefix:
}
\details{
\itemize{
\item BETYdb http://www.betydb.org - \code{betydb_}
\item National Center for Biotechnology Information - NCBI
http://www.ncbi.nlm.nih.gov/ - \code{ncbi_}
\item Encyclopedia of Life Traitbank - \code{traitbank_}
\item Birdlife International https://www.birdlife.org/ -
\code{birdlife_}
\item LEDA Traitbase http://www.leda-traitbase.org/LEDAportal/index.jsp -
\code{leda_}
\item Zanne et al. plant dataset - \code{\link[=tr_zanne]{tr_zanne()}}
\item Amniote life history dataset - \code{\link[=tr_ernest]{tr_ernest()}}
}

See also \link{traits-defunct}
}
\author{
Scott Chamberlain

Ignasi Bartomeus

Zachary Foster

David LeBauer

David Harris

Rupert Collins
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_native.R
\name{is_native}
\alias{is_native}
\title{Check if a species is native somewhere}
\usage{
is_native(...)
}
\description{
Check if a species is native somewhere
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_byid.R
\name{ncbi_byid}
\alias{ncbi_byid}
\title{Retrieve gene sequences from NCBI by accession number.}
\usage{
ncbi_byid(ids, format = NULL, verbose = TRUE)
}
\arguments{
\item{ids}{(character) GenBank ids to search for. One or more. Required.}

\item{format}{(character) Return type, e.g., \code{"fasta"}. NOW IGNORED.}

\item{verbose}{(logical) If \code{TRUE} (default), informative messages
printed.}
}
\value{
data.frame of the form:
\itemize{
\item taxon - taxonomic name (may include some junk, but hard to parse off)
\item taxonomy - organism lineage
\item gene_desc - gene description
\item organelle - if mitochondrial or chloroplast
\item gi_no - GI number
\item acc_no - accession number
\item keyword - if official DNA barcode
\item specimen_voucher - museum/lab accession number of vouchered material
\item lat_lon - longitude/latitude of specimen collection event
\item country - country/location of specimen collection event
\item paper_title - title of study
\item journal - journal study published in (if published)
\item first_author - first author of study
\item uploaded_date - date sequence was uploaded to GenBank
\item length - sequence length
\item sequence - sequence character string
}
}
\description{
Retrieve gene sequences from NCBI by accession number.
}
\details{
If bad ids are included with good ones, the bad ones are
silently dropped. If all ids are bad you'll get a stop with error message.
}
\examples{
\dontrun{
# A single gene
ncbi_byid(ids="360040093")

# Many genes (with different accession numbers)
ncbi_byid(ids=c("360040093","347448433"))
}
}
\seealso{
\code{\link[=ncbi_searcher]{ncbi_searcher()}}, ncbi_byname()]
}
\author{
Scott Chamberlain, Rupert Collins
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coral.R
\name{coral-defunct}
\alias{coral-defunct}
\alias{coral_taxa}
\alias{coral_traits}
\alias{coral_locations}
\alias{coral_methodologies}
\alias{coral_resources}
\alias{coral_species}
\title{Search for coral data on coraltraits.org}
\usage{
coral_taxa(...)

coral_traits(...)

coral_locations(...)

coral_methodologies(...)

coral_resources(...)

coral_species(...)
}
\arguments{
\item{...}{ignored}
}
\description{
DEFUNCT: As far as we can tell the coraltraits.org API is down
for good
}
\references{
https://coraltraits.org/
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_search.R
\name{taxa_search}
\alias{taxa_search}
\title{Search for traits by taxa names}
\usage{
taxa_search(x, db, ...)
}
\arguments{
\item{x}{(character) Taxonomic name(s) to search for}

\item{db}{(character) only 'ncbi' for now - other options
maybe in the future}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
A \code{data.frame}
}
\description{
Search for traits by taxa names
}
\examples{
\dontrun{
taxa_search("Poa annua", db = "ncbi")
}
}
\author{
Scott Chamberlain
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tr_usda.R
\name{tr_usda}
\alias{tr_usda}
\title{USDA plants data}
\usage{
tr_usda(...)
}
\arguments{
\item{...}{ignored}
}
\description{
DEFUNCT
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{traits_cache}
\alias{traits_cache}
\title{Caching}
\description{
Manage cached \code{traits} package files with \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/traits")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{traits_cache$cache_path_get()} get cache path
\item \code{traits_cache$cache_path_set()} set cache path
\item \code{traits_cache$list()} returns a character vector of full
path file names
\item \code{traits_cache$files()} returns file objects with metadata
\item \code{traits_cache$details()} returns files with details
\item \code{traits_cache$delete()} delete specific files
\item \code{traits_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
traits_cache

# list files in cache
traits_cache$list()

# delete certain database files
# traits_cache$delete("file path")
# traits_cache$list()

# delete all files in cache
# traits_cache$delete_all()
# traits_cache$list()

# set a different cache path from the default
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_searcher.R
\name{ncbi_searcher}
\alias{ncbi_searcher}
\title{Search for gene sequences available for taxa from NCBI.}
\usage{
ncbi_searcher(
  taxa = NULL,
  id = NULL,
  seqrange = "1:3000",
  getrelated = FALSE,
  fuzzy = FALSE,
  limit = 500,
  entrez_query = NULL,
  hypothetical = FALSE,
  verbose = TRUE,
  sleep = 0L
)
}
\arguments{
\item{taxa}{(character) Scientific name to search for.}

\item{id}{(\code{character}) Taxonomic id to search for. Not compatible with
argument \code{taxa}.}

\item{seqrange}{(character) Sequence range, as e.g., \code{"1:1000"}. This is the range of
sequence lengths to search for. So \code{"1:1000"} means search for sequences from 1 to 1000
characters in length.}

\item{getrelated}{(logical) If \code{TRUE}, gets the longest sequences of a species
in the same genus as the one searched for. If \code{FALSE}, returns nothing if no match
found.}

\item{fuzzy}{(logical) Whether to do fuzzy taxonomic ID search or exact
search. If \code{TRUE}, we use \code{xXarbitraryXx[porgn:__txid<ID>]},
but if \code{FALSE}, we use \code{txid<ID>}. Default: \code{FALSE}}

\item{limit}{(\code{numeric}) Number of sequences to search for and return.
Max of 10,000. If you search for 6000 records, and only 5000 are found,
you will of course only get 5000 back.}

\item{entrez_query}{(\code{character}; length 1) An Entrez-format query to
filter results with. This is useful to search for sequences with specific
characteristics. The format is the same as the one used to seach genbank.
(\url{https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options})}

\item{hypothetical}{(\code{logical}; length 1) If \code{FALSE}, an attempt
will be made to not return hypothetical or predicted sequences judging from
accession number prefixs (XM and XR). This can result in less than the
\code{limit} being returned even if there are more sequences available,
since this filtering is done after searching NCBI.}

\item{verbose}{(logical) If \code{TRUE} (default), informative messages printed.}

\item{sleep}{(integer) number of seconds to sleep before each HTTP request.
use if running to 429 Too Many Requests errors from NCBI. default: 0
(no sleep)}
}
\value{
\code{data.frame} of results if a single input is given. A list of
\code{data.frame}s if multiple inputs are given.
}
\description{
Search for gene sequences available for taxa from NCBI.
}
\section{Authentication}{

NCBI rate limits requests. If you set an API key you have a higher rate limit.
Set your API key like \code{Sys.setenv(ENTREZ_KEY="yourkey")} or you can use
\code{?rentrez::set_entrez_key}. set verbose curl output (\code{crul::set_verbose()}) to
make sure your api key is being sent in the requests
}

\examples{
\dontrun{
# A single species
out <- ncbi_searcher(taxa="Umbra limi", seqrange = "1:2000")
# Get the same species information using a taxonomy id
out <- ncbi_searcher(id = "75935", seqrange = "1:2000")
# If the taxon name is unique, using the taxon name and id are equivalent
all(ncbi_searcher(id = "75935") ==  ncbi_searcher(taxa="Umbra limi"))
# If the taxon name is not unique, use taxon id
#  "266948" is the uid for the butterfly genus, but there is also a genus
#  of orchids with the
#  same name
nrow(ncbi_searcher(id = "266948")) ==  nrow(ncbi_searcher(taxa="Satyrium"))
# get list of genes available, removing non-unique
unique(out$gene_desc)
# does the string 'RAG1' exist in any of the gene names
out[grep("RAG1", out$gene_desc, ignore.case=TRUE),]

# A single species without records in NCBI
out <- ncbi_searcher(taxa="Sequoia wellingtonia", seqrange="1:2000",
  getrelated=TRUE)

# Many species, can run in parallel or not using plyr
species <- c("Salvelinus alpinus","Ictalurus nebulosus","Carassius auratus")
out2 <- ncbi_searcher(taxa=species, seqrange = "1:2000")
lapply(out2, head)
library("plyr")
out2df <- ldply(out2) # make data.frame of all
unique(out2df$gene_desc) # get list of genes available, removing non-unique
out2df[grep("12S", out2df$gene_desc, ignore.case=TRUE), ]

# Using the getrelated and entrez_query options
ncbi_searcher(taxa = "Olpidiopsidales", limit = 5, getrelated = TRUE,
            entrez_query = "18S[title] AND 28S[title]")

# get refseqs
one <- ncbi_searcher(taxa = "Salmonella enterica",
  entrez_query="srcdb_refseq[PROP]")
two <- ncbi_searcher(taxa = "Salmonella enterica")
}
}
\seealso{
\code{\link{ncbi_byid}}, \code{\link{ncbi_byname}}
}
\author{
Scott Chamberlain, Zachary Foster \email{zacharyfoster1989@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betydb.R
\name{betydb_query}
\alias{betydb_query}
\alias{betydb_search}
\title{Query a BETY table}
\usage{
betydb_query(
  ...,
  table = "search",
  key = NULL,
  api_version = NULL,
  betyurl = NULL,
  user = NULL,
  pwd = NULL,
  progress = TRUE
)

betydb_search(
  query = "Maple SLA",
  ...,
  include_unchecked = NULL,
  progress = TRUE
)
}
\arguments{
\item{...}{(named character) Columns to query, as key="value" pairs. Note that betydb_query passes these along to BETY with no check whether the requested keys exist in the specified table.}

\item{table}{(character) The name of the database table to query, or "search" (the default) for the traits and yields view}

\item{key}{(character) An API key. Use this or user/pwd combo.
Save in your \code{.Rprofile} file as \code{options(betydb_key = "your40digitkey")}. Optional}

\item{api_version}{(character) Which version of the BETY API should we query? One of "v0" or "beta".
Default is \code{options("betydb_api_version")} if set, otherwise  "v0".}

\item{betyurl}{(string) url to target instance of betydb.
Default is \code{options("betydb_url")} if set, otherwise "https:/www.betydb.org/"}

\item{user, pwd}{(character) A user name and password. Use a user/pwd combo or an API key.
Save in your \code{.Rprofile} file as \code{options(betydb_user = "yournamehere")} and \code{options(betydb_pwd = "yourpasswordhere")}. Optional}

\item{progress}{show progress bar? default: \code{TRUE}}

\item{query}{(character) A string containing one or more words to be queried across all columns of the "search" table.}

\item{include_unchecked}{(logical) Include results that have not been quality checked? Applies only to tables with a "checked" column: "search", "traits", "yields". Default is to exclude unchecked values.}
}
\value{
A data.frame with attributes containing request metadata, or NULL if the query produced no results
}
\description{
Query a BETY table
}
\details{
Use betydb_query to retrieve records from a table that match on all the column filters specified in '...'.
If no filters are specified, retrieves the whole table. In API versions that support it (i.e. not in v0), filter strings beginning with "~" are treated as regular expressions.
}
\examples{
\dontrun{
# literal vs regular expression vs anchored regular expression:
betydb_query(units = "Mg", table = "variables")
# NULL
betydb_query(units = "Mg/ha", table = "variables") \%>\% select(name) \%>\% c()
# $name
# [1] "a_biomass"                  "root_live_biomass"
# [3] "leaf_dead_biomass_in_Mg_ha" "SDM"

betydb_query(genus = "Miscanthus", table = "species") \%>\% nrow()
# [1] 10
(betydb_query(genus = "~misc", table = "species", api_version = "beta")
 \%>\% select(genus)
 \%>\% unique() \%>\% c())
# $genus
# [1] "Platymiscium" "Miscanthus"   "Dermiscellum"

(betydb_query(genus = "~^misc", table = "species", api_version = "beta")
 \%>\% select(genus)
 \%>\% unique() \%>\% c())
# $genus
# [1] "Miscanthus"
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{traits-defunct}
\alias{traits-defunct}
\title{Defunct functions in traits}
\description{
These functions have been removed.
}
\details{
\itemize{
\item \code{eol_invasive_}: This function has moved to a new package.
See \code{originr::eol}
\item \code{fe_native}: This function has moved to a new package.
See \code{originr::flora_europaea}
\item \code{g_invasive}: This function has moved to a new package.
See \code{originr::gisd}
\item \code{is_native}: This function has moved to a new package. See
\code{originr::is_native}
\item \code{tr_usda}: the API behind this function is down for good
\item \code{coral_locations}: API down for good, as far as I can tell
\item \code{coral_methodologies}: API down for good, as far as I can tell
\item \code{coral_resources}: API down for good, as far as I can tell
\item \code{coral_species}: API down for good, as far as I can tell
\item \code{coral_taxa}: API down for good, as far as I can tell
\item \code{coral_traits}: API down for good, as far as I can tell
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traits-package.r
\docType{data}
\name{plantatt}
\alias{plantatt}
\title{PLANTATT plant traits dataset}
\description{
PLANTATT plant traits dataset
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fe_native.R
\name{fe_native}
\alias{fe_native}
\title{Check species status (native, exotic, ...) from Flora Europaea webpage}
\usage{
fe_native(...)
}
\description{
Check species status (native, exotic, ...) from Flora Europaea webpage
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g_invasive.R
\name{g_invasive}
\alias{g_invasive}
\title{Check invasive species status for a set of species from GISD database}
\usage{
g_invasive(...)
}
\description{
Check invasive species status for a set of species from GISD database
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/birdlife.R
\name{birdlife_threats}
\alias{birdlife_threats}
\title{Get bird threat information from BirdLife/IUCN}
\usage{
birdlife_threats(id)
}
\arguments{
\item{id}{A single IUCN species ID}
}
\value{
a \code{data.frame} with the species ID and two levels of threat
descriptions, plus stresses, timing, scope, severity, and impact associated
with each stressor.
}
\description{
Get bird threat information from BirdLife/IUCN
}
\examples{
\dontrun{
# Setophaga chrysoparia
birdlife_threats(22721692)
# Aburria aburri
birdlife_threats(22678440)
}
}
\seealso{
Other birdlife: 
\code{\link{birdlife_habitat}()}
}
\author{
David J. Harris \email{harry491@gmail.com}
}
\concept{birdlife}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tr_ernest.R
\name{tr_ernest}
\alias{tr_ernest}
\title{Amniote life history dataset}
\usage{
tr_ernest(read = TRUE, ...)
}
\arguments{
\item{read}{(logical) read in csv files. Default: \code{TRUE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
paths to the files (character) if \code{read=FALSE} or
a list of data.frame's if \code{read=TRUE}
}
\description{
Amniote life history dataset
}
\details{
When using this data, cite the paper:

Myhrvold, N. P., Baldridge, E., Chan, B., Sivam, D., Freeman, D. L. and
Ernest, S. K. M. (2015), An amniote life-history database to perform
comparative analyses with birds, mammals, and reptiles. Ecology, 96: 3109.
https://doi.org/10.1890/15-0846R.1

As well as the Dryad data package:

L. Freeman, Daniel; P. Myhrvold, Nathan; Chan, Benjamin; Sivam, Dhileep;
Ernest, S. K. Morgan; Baldridge, Elita (2016): Full Archive. figshare.
https://doi.org/10.6084/m9.figshare.3563457.v1
}
\examples{
\dontrun{
res <- tr_ernest()
res$data
res$references
res$sparse
res$range_count
}
}
\references{
https://doi.org/10.1890/15-0846R.1
https://doi.org/10.6084/m9.figshare.3563457.v1
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traitbank.R
\name{traitbank}
\alias{traitbank}
\title{Search for traits from EOL's Traitbank}
\usage{
traitbank(query, key = NULL, ...)
}
\arguments{
\item{query}{(character) a query to the EOL Cypher service that holds
Traitbank data. required. no default query given. see examples}

\item{key}{(character) EOL Cypher query API key. required, either
passed in or as an environment variable}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
a list
}
\description{
Search for traits from EOL's Traitbank
}
\details{
\code{traitbank} is an interface to the EOL Cypher query.
Note that the previous interface EOL had for Traits has been completely
replaced - thus, this function is completely different. You no longer
query by EOL page id, but using the query language for a database
called Neo4J. See the docs for help. Later we plan to make a more
user friendly interface to get Traitbank data that doesn't
require knowing the Neo4J query syntax
}
\section{Authentication}{

You'll need an EOL cypher key to use this function. Get one by signing
in to your EOL account https://eol.org/users/sign_in then head to
https://eol.org/services/authenticate to get a key. Store your key
in your .Renviron file or similar under the name "EOL_CYPHER_KEY",
and we will use that key in this function. Alternatively, you can
pass in your key to the key parameter, but we do not recommend
doing that as you risk accidentally committing your key to the
public web.
}

\examples{
\dontrun{
# traitbank_query function
traitbank(query = "MATCH (n:Trait) RETURN n LIMIT 1;")

# traitbank function
res <- traitbank(query = "MATCH (n:Trait) RETURN n LIMIT 2;")
res
}
}
\references{
https://github.com/EOL/eol_website/blob/master/doc/api.md
https://github.com/EOL/eol_website/blob/master/doc/query-examples.md
}
