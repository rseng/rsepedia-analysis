citecorp
=========



[![cran checks](https://cranchecks.info/badges/worst/citecorp)](https://cranchecks.info/pkgs/citecorp)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/citecorp/workflows/R-check/badge.svg)](https://github.com/ropensci/citecorp/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/citecorp)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/citecorp)](https://cran.r-project.org/package=citecorp)

Client for the Open Citations Corpus http://opencitations.net/ (OCC)

OCC created their own identifiers called Open Citation Identifiers (oci), e.g., 

```
020010009033611182421271436182433010601-02001030701361924302723102137251614233701000005090307
```

You are probably not going to be using oci identifiers, but rather DOIs and/or PMIDs
and/or PMCIDs. See `?oc_lookup` for methods for cross-walking among identifier types.

If you'd like to use the OpenCitations Sparql endpoint yourself you can find that
at http://opencitations.net/sparql


## Install

CRAN version


```r
install.packages("citecorp")
```

Development version


```r
remotes::install_github("ropensci/citecorp")
```


```r
library("citecorp")
```

## Methods for converting IDs


```r
oc_doi2ids("10.1097/igc.0000000000000609")
#>                            doi                           paper      pmcid
#> 1 10.1097/igc.0000000000000609 https://w3id.org/oc/corpus/br/1 PMC4679344
#>       pmid
#> 1 26645990
oc_pmid2ids("26645990")
#>                            doi                           paper      pmcid
#> 1 10.1097/igc.0000000000000609 https://w3id.org/oc/corpus/br/1 PMC4679344
#>       pmid
#> 1 26645990
oc_pmcid2ids("PMC4679344")
#>                            doi                           paper      pmcid
#> 1 10.1097/igc.0000000000000609 https://w3id.org/oc/corpus/br/1 PMC4679344
#>       pmid
#> 1 26645990
```

You can pass in more than one identifer to each of the above functions:


```r
oc_doi2ids(oc_dois[1:6])
#>                                  doi                                 paper
#> 1               10.1128/jvi.00758-10 https://w3id.org/oc/corpus/br/5357460
#> 2 10.1111/j.2042-3306.1989.tb02167.x  https://w3id.org/oc/corpus/br/589891
#> 3       10.1097/rli.0b013e31821eea45 https://w3id.org/oc/corpus/br/3931705
#> 4           10.1177/0148607114529597 https://w3id.org/oc/corpus/br/5016780
#> 5            10.1111/1567-1364.12217 https://w3id.org/oc/corpus/br/3819297
#> 6      10.1016/s0168-9525(99)01798-9 https://w3id.org/oc/corpus/br/4606537
#>        pmcid     pmid
#> 1 PMC2953162 20702630
#> 2       <NA>  2670542
#> 3       <NA> 21577119
#> 4       <NA> 24711119
#> 5       <NA> 25263709
#> 6       <NA> 10461200
```

## COCI methods

OpenCitations Index of Crossref open DOI-to-DOI references

If you don't load `tibble` you get normal data.frame's


```r
library(tibble)
doi1 <- "10.1108/jd-12-2013-0166"
# references
oc_coci_refs(doi1)
#> # A tibble: 37 x 7
#>    journal_sc author_sc timespan citing    oci             cited        creation
#>  * <chr>      <chr>     <chr>    <chr>     <chr>           <chr>        <chr>   
#>  1 no         no        P9Y2M5D  10.1108/… 02001010008361… 10.1001/jam… 2015-03…
#>  2 no         no        P41Y8M   10.1108/… 02001010008361… 10.1002/asi… 2015-03…
#>  3 no         no        P25Y6M   10.1108/… 02001010008361… 10.1002/(si… 2015-03…
#>  4 no         no        P17Y2M   10.1108/… 02001010008361… 10.1007/bf0… 2015-03…
#>  5 no         no        P2Y2M3D  10.1108/… 02001010008361… 10.1007/s10… 2015-03…
#>  6 no         no        P5Y8M27D 10.1108/… 02001010008361… 10.1007/s11… 2015-03…
#>  7 no         no        P2Y3M    10.1108/… 02001010008361… 10.1016/j.w… 2015-03…
#>  8 no         no        P1Y10M   10.1108/… 02001010008361… 10.1016/j.w… 2015-03…
#>  9 no         no        P12Y     10.1108/… 02001010008361… 10.1023/a:1… 2015-03…
#> 10 no         no        P13Y10M  10.1108/… 02001010008361… 10.1038/350… 2015-03…
#> # … with 27 more rows
# citations
oc_coci_cites(doi1)
#> # A tibble: 23 x 7
#>    journal_sc author_sc timespan  citing     oci               cited    creation
#>  * <chr>      <chr>     <chr>     <chr>      <chr>             <chr>    <chr>   
#>  1 no         no        P3Y       10.1145/3… 0200101040536030… 10.1108… 2018    
#>  2 no         no        P2Y5M     10.1057/s… 0200100050736280… 10.1108… 2017-08 
#>  3 no         no        P4Y1M1D   10.3233/d… 0200302030336132… 10.1108… 2019-04…
#>  4 no         no        P4Y5M10D  10.3233/d… 0200302030336132… 10.1108… 2019-08…
#>  5 no         no        P1Y0M14D  10.3233/s… 0200302030336283… 10.1108… 2016-03…
#>  6 no         no        P3Y10M12D 10.3233/s… 0200302030336283… 10.1108… 2019-01…
#>  7 no         no        P3Y6M     10.1142/s… 0200101040236280… 10.1108… 2018-09 
#>  8 no         no        P2Y11M20D 10.7554/e… 0200705050436142… 10.1108… 2018-03…
#>  9 no         no        P0Y       10.3346/j… 0200303040636192… 10.1108… 2015    
#> 10 no         no        P3Y       10.1007/9… 0200100000736090… 10.1108… 2018    
#> # … with 13 more rows
# metadata
oc_coci_meta(doi1)
#> # A tibble: 1 x 13
#>   doi   reference issue source_id citation page  volume author citation_count
#> * <chr> <chr>     <chr> <chr>     <chr>    <chr> <chr>  <chr>  <chr>         
#> 1 10.1… 10.1001/… 2     issn:002… 10.1145… 253-… 71     Peron… 23            
#> # … with 4 more variables: year <chr>, source_title <chr>, title <chr>,
#> #   oa_link <chr>
```


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/citecorp/issues)
* License: MIT
* Get citation information for `citecorp` in R doing `citation(package = 'citecorp')`
* Please note that this project is released with a [Contributor Code of Conduct][coc].
By participating in this project you agree to abide by its terms.

[sparqldsl]: https://github.com/ropensci/sparqldsl
[coc]: https://github.com/ropensci/citecorp/blob/master/CODE_OF_CONDUCT.md


[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
citecorp 0.3.0
==============

### NEW FEATURES

* `oc_coci_refs()`, `oc_coci_cites()`, and `oc_coci_citation()` now support more than one DOI or COCI identifier passed in. `oc_coci_meta()` already supported more than one identifier; now all four `oc_coci*` methods support more than one identifier. The caveat is that the API route behindn `oc_coci_meta()` supports >1 identifier in a single http request - whereas the API routes behind the other functions only support one identifier per request, so they make an http request for each identifier passed in (#6)

### MINOR IMPROVEMENTS

* added tests for Cloudflare http status codes (#7)


citecorp 0.2.2
==============

### BUG FIXES

* better check if http requests will work before running examples that require http requests (#5)


citecorp 0.2.0
==============

### NEW FEATURES

* the functions `oc_doi2ids`, `oc_pmid2ids`, and `oc_pmcid2ids` (see `?oc_lookup`) now all accept more than 1 identifier. the underlying SPARQL queries were modified to make this possible. (#1)
* three new package datasets (oc_dois, oc_pmids, oc_pmcids) with vectors of identifiers of each of those types for testing and examples

### BUG FIXES

* fixed an error in parsing results in `oc_lookup` functions - this fix was later changed as part of another fix, but still appreciated (#2) (#3) thanks @Selbosh
* fix to parsing in `oc_lookup` functions: only retrieve data.frame columns if they exist - prevents error on retrieving a column that doesn't exist (#4) thanks @Selbosh


citecorp 0.1.0
==============

### NEW FEATURES

* Released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local OS X install, R 3.6.3 patched
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* There are no reverse dependencies.

---

This version adds more tests and changes some functions to support multiple input object identifiers.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/citecorp/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/crul.git`
* Make sure to track progress upstream (i.e., on our version of `crul` at `ropensci/citecorp`) by doing `git remote add upstream https://github.com/ropensci/citecorp.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/citecorp`

### rOpenSci Forum 

Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
