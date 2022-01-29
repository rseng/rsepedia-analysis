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
