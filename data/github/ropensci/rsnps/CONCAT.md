
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsnps

[![cran
checks](https://cranchecks.info/badges/worst/rsnps)](https://cranchecks.info/pkgs/rsnps/)
[![R build
status](https://github.com/ropensci/rsnps/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rsnps/actions)
[![Build
status](https://ci.appveyor.com/api/projects/status/d2lv98726u6t9ut5/branch/master)](https://ci.appveyor.com/project/sckott/rsnps/branch/master/)
[![codecov.io](https://codecov.io/github/ropensci/rsnps/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rsnps?branch=master)
[![cran
version](https://www.r-pkg.org/badges/version/rsnps)](https://cran.r-project.org/package=rsnps)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/rsnps?color=E664A4)](https://github.com/r-hub/cranlogs.app)

This package gives you access to data from OpenSNP and NCBI’s dbSNP SNP
database.

## NOTE

`rsnps` used to be `ropensnp`

## Data sources

This set of functions/package accesses data from:

-   openSNP.org
    -   <https://opensnp.org>
    -   See documentation on the openSNP API
        <https://opensnp.org/faq#api>
    -   See blog post about their API
        <https://opensnp.wordpress.com/2012/01/18/some-progress-on-the-api-json-endpoints/>
    -   Relevant functions:
        -   `allgensnp()`, `allphenotypes()`, `annotations()`,
            `download_users()`, `fetch_genotypes()`, `genotypes()`,
            `phenotypes()`, `phenotypes_byid()`, `users()`
-   NCBI’s dbSNP SNP database
    -   See <https://www.ncbi.nlm.nih.gov/snp/> for more details
    -   Relevant function:
        -   `ncbi_snp_query()`

## Install

Install from CRAN

``` r
install.packages("rsnps")
```

Or dev version

``` r
install.packages("remotes")
remotes::install_github("ropensci/rsnps")
```

``` r
library("rsnps")
```

## Usage

### NCBI dbSNP data

``` r
snps <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
ncbi_snp_query(snps)
```

    #> # A tibble: 4 × 16
    #>   query        chromosome        bp class rsid   gene   alleles ancestral_allele
    #>   <chr>        <chr>          <dbl> <chr> <chr>  <chr>  <chr>   <chr>           
    #> 1 rs332        7          117559593 del   rs121… "CFTR… TTT, d… TTT             
    #> 2 rs420358     1           40341239 snv   rs420… ""     A,C,G,T A               
    #> 3 rs1837253    5          111066174 snv   rs183… ""     T,C     T               
    #> 4 rs1209415715 9           41782316 snv   rs120… ""     T,A,C   T               
    #> # … with 8 more variables: variation_allele <chr>, seqname <chr>, hgvs <chr>,
    #> #   assembly <chr>, ref_seq <chr>, minor <chr>, maf <dbl>,
    #> #   maf_population <list>

### openSNP data

`genotypes()` function

``` r
genotypes('rs9939609', userid='1,6,8', df=TRUE)
```

    #>    snp_name snp_chromosome snp_position                 user_name user_id
    #> 1 rs9939609             16     53786615 Bastian Greshake Tzovaras       1
    #> 2 rs9939609             16     53786615              Nash Parovoz       6
    #> 3 rs9939609             16     53786615         Samantha B. Clark       8
    #>   genotype_id genotype
    #> 1           9       AT
    #> 2           5       AT
    #> 3           2       TT

`phenotypes()` function

``` r
out <- phenotypes(userid=1)
out$phenotypes$`Hair Type`
```

    #> $phenotype_id
    #> [1] 16
    #> 
    #> $variation
    #> [1] "straight"

For more detail, see the [vignette: rsnps
tutorial](https://github.com/ropensci/rsnps/tree/master/vignettes).

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/rsnps/issues/).
-   License: MIT
-   Get citation information for `rsnsps` in R doing
    `citation(package = 'rsnps')`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rsnps 0.5.0
===========

### NEW FEATURES

* ncbi_snp_query(): enable allele frequency for different reference populations, ncbi_snp_query() outputs now a tibble (#97).

### MINOR IMPROVEMENTS

* ncbi_snp_query(): replace calls to RJSONIO with equivalent in jsonlite (#98).
* unit tests for ncbi_snp_query(): added a tolerance to any allele frequency checks.
* move vignette source to /vignettes and precompute using an R script.

### DOCUMENTATION FIXES

* Updated vignette.

rsnps 0.4.0
===========

### MAJOR IMPROVEMENTS

NCBI / dbSNP changed their API:

* Rewrote `ncbi_snp_query` to accommodate the new API (#86, #88). 
* Removed the functions `ncbi_snp_query2` an `ncbi_snp_summary`. 

### MINOR IMPROVEMENTS

* Reordered `ncbi_snp_query` dataframe output to have chromosome and bp beside each other (#70).
* Changed `ncbi_snp_query` parameter (`SNPs`) to lower case (`snps`). 

### DOCUMENTATION FIXES

* Restructured and fixed a typo in `README.Rmd` and added link to vignette (#63).
* Added info of two new maintainers to `DESCRIPTION`. 
* Added relevant API links to vignette. 

### BUG FIXES

* Fixed the test for `allphenotypes` function by making it less specific (#72). 


rsnps 0.3.0
===========

### DEPRECATED AND DEFUNCT

* `ld_search()` is now defunct. The Broad Institute has taken down the SNAP service behind the function. (#46) (#53) (#60)

### NEW FEATURES

* the three NCBI functions gain a new parameter `key` for passing in an NCBI Entrez API key. You can alternatively (and we encourage this) store your key as an environment variable and we'll use that instead. The key allows you to have higher rate limits than without a key (#58)
* gains new function `ncbi_snp_summary()` for summary data on a SNP (#31)

### MINOR IMPROVEMENTS

* http requests are now done using `crul` instead of `httr` (#44)
* now using markdown formatted documentation (#56)
* documented in `ncbi_snp_query()` that we can not change the assembly (#49)

### BUG FIXES

* fix to `ncbi_snp_query2()`: when many IDs passed in, we were failing with a "URI too long" message. We now check how many Ids are passed in and do a POST request as needed  (#39)
* fixed problem in `ncbi_snp_query()` where it wasn't pulling out correctly the gene name and BP position (#25)



rsnps 0.2.0
===========

### NEW FEATURES

* `LDSearch()` is now `ld_search()`, but `LDSearch()` still works until 
the next CRAN release when it will be defunct (#33)
* `NCBI_snp_query()` is now `ncbi_snp_query()`, but `NCBI_snp_query()` still 
works until the next CRAN release when it will be defunct (#33)
* `NCBI_snp_query2()` is now `ncbi_snp_query2()`, but `NCBI_snp_query2()` still 
works until the next CRAN release when it will be defunct (#33)

### MINOR IMPROVEMENTS

* Namespace all base R package function calls (#21)
* Improve `httr::content` call to parse to text, and `encoding = "UTF-8"` 
(#24)
* Added tests for `ld_search()` (#12)
* Added tests for `ncbi_snp_query()` and `ncbi_snp_query2()` (#13)
* Added ancestral allele output to `ncbi_snp_query()` (#23)

### BUG FIXES

* Fix to `fetch_genotypes()`, was failing sometimes when the commented
metadata lines at top varied in length (#22)
* Fix to `ld_search()` (#32)


rsnps 0.1.6
===========

### MINOR IMPROVEMENTS

* All examples now in `\dontrun`. (#11)
* Added additional tests for `LDSearch()` and `NCBI_snp_query()`.
* Added a vignette.

### BUG FIXES

* Bugs fixed in `LDSearch()`, which were actually bugs in `NCBI_snp_query()`. (#9)
* Bug fixed in `NCBI_snp_query()` as chromosome might also be "X". 

rsnps 0.1.0
===========

### NEW FEATURES 

* Bug fixes to all openSNP functions.

rsnps 0.0.5
===========

### NEW FEATURES 

* released to CRAN
## Test environments

* local OS X install, R 4.1.1
* Ubuntu Linux 20.04.1 LTS (on R-hub), R 4.1.2
* Fedora Linux (on R-hub) R-devel
* Windows (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs. 

There is one NOTE that is only found on Windows (Server 2022, R-devel 64-bit): 

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

## Reverse dependencies

* We have run R CMD check on the 1 downstream dependency
(<https://github.com/ropensci/rsnps/blob/master/revdep/README.md>).
No problems were found related to this package.

---

This version includes a new feature and two minor improvement of the function `ncbi_snp_query`,
and a change in the vignette (we now pre-compile the vignette to avoid long runtimes). 


Thanks!
Julia Gustavsen and Sina Rüeger
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

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the [Issues page](https://github.com/{owner}/{repo}/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `{repo}` at `{owner}/{repo}`) by doing `git remote add upstream https://github.com/{owner}/{repo}.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `{owner}/{repo}`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

### Thanks for contributing!
---
name: Bug report
about: Create a report to help us improve rsnps
title: ''
labels: 'Bug'
assignees: ''

---

**Please do not submit bug reports for issues with the OpenSNP or NCBI’s dbSNP server itself. rsnps is only an API client, it does not have control over the server hosting the data.**

**Describe the bug**
A clear and concise description of what the bug is.

**Code to reproduce**
A clear set of code illustrating the steps to reproduce the behaviour.

**Expected behaviour**
A clear and concise description of what you expected to happen.

**OS and R versions (please complete the following information):**
 - OS: [e.g., Windows10, macOS, Linux]
 - R Version [e.g., 3.6.2]

**Session Info: `devtools::session_info()` or `sessionInfo()`**

<details> <summary><strong>Session Info</strong></summary>
  
```r

```

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for rsnps
title: ''
labels: 'Feature'
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. 

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.2 (2019-12-12) |
|os       |macOS Mojave 10.14.6         |
|system   |x86_64, darwin15.6.0         |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_CA.UTF-8                  |
|ctype    |en_CA.UTF-8                  |
|tz       |Europe/Zurich                |
|date     |2020-06-21                   |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|rsnps   |0.3.0 |0.3.2.9121 |*  |
|plyr    |NA    |1.8.6      |*  |
|xml2    |NA    |1.3.2      |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

rsnps
=====

```{r, eval=TRUE, echo=FALSE}
knitr::opts_chunk$set(
  warning=FALSE,
  message=FALSE,
  comment="#>"
)
```

[![cran checks](https://cranchecks.info/badges/worst/rsnps)](https://cranchecks.info/pkgs/rsnps/)
[![R build status](https://github.com/ropensci/rsnps/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rsnps/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/d2lv98726u6t9ut5/branch/master)](https://ci.appveyor.com/project/sckott/rsnps/branch/master/)
[![codecov.io](https://codecov.io/github/ropensci/rsnps/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rsnps?branch=master)
[![cran version](https://www.r-pkg.org/badges/version/rsnps)](https://cran.r-project.org/package=rsnps)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rsnps?color=E664A4)](https://github.com/r-hub/cranlogs.app)


This package gives you access to data from OpenSNP and NCBI's dbSNP SNP database.

## NOTE

`rsnps` used to be `ropensnp`


## Data sources

This set of functions/package accesses data from:

+ openSNP.org
	+ <https://opensnp.org>
	+ See documentation on the openSNP API <https://opensnp.org/faq#api>
	+ See blog post about their API <https://opensnp.wordpress.com/2012/01/18/some-progress-on-the-api-json-endpoints/>
	+ Relevant functions:
		+ `allgensnp()`, `allphenotypes()`, `annotations()`, `download_users()`, 
		`fetch_genotypes()`, `genotypes()`, `phenotypes()`, `phenotypes_byid()`, `users()`

+ NCBI's dbSNP SNP database
	+ See <https://www.ncbi.nlm.nih.gov/snp/> for more details
	+ Relevant function:
		+ `ncbi_snp_query()`

## Install

Install from CRAN

```{r eval=FALSE}
install.packages("rsnps")
```

Or dev version

```{r eval=FALSE}
install.packages("remotes")
remotes::install_github("ropensci/rsnps")
```

```{r}
library("rsnps")
```

## Usage

### NCBI dbSNP data

```{r}
snps <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
ncbi_snp_query(snps)
```

### openSNP data

`genotypes()` function

```{r}
genotypes('rs9939609', userid='1,6,8', df=TRUE)
```

`phenotypes()` function

```{r}
out <- phenotypes(userid=1)
out$phenotypes$`Hair Type`
```

For more detail, see the [vignette: rsnps tutorial](https://github.com/ropensci/rsnps/tree/master/vignettes).

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rsnps/issues/).
* License: MIT
* Get citation information for `rsnsps` in R doing `citation(package = 'rsnps')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
output:
  pdf_document: default
  html_document: default
---
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rsnps tutorial}
%\VignetteEncoding{UTF-8}
-->



# rsnps tutorial

## Install and load library

When available on CRAN


```r
install.packages("rsnps")
```

Or get from Github


```r
install.packages("devtools")
devtools::install_github("ropensci/rsnps")
```



```r
library(rsnps)
```

## OpenSNP data

### All Genotypes

Get genotype data for all users at a particular SNP from [OpenSNP](https://opensnp.org):


```r
x <- allgensnp(snp='rs7412')
head(x)
#>     name chromosome position                name   id genotype_id local_genotype
#> 1 rs7412         19 44908822        R.M. Holston   22           8             CC
#> 2 rs7412         19 44908822 Charles G. Sullivan 5326        3834             CC
#> 3 rs7412         19 44908822   Glenn Allen Nolen   19           7             CC
#> 4 rs7412         19 44908822        Angel Harris  495         223             CC
#> 5 rs7412         19 44908822           Mom to AG  387         173             CC
#> 6 rs7412         19 44908822            kevinmcc  285         118             CC
```


### All Phenotypes

Get all phenotypes, their variations, and how many users have data available for a given phenotype

Get all data


```r
x <- allphenotypes(df = TRUE)
head(x)
#>   id characteristic known_variations number_of_users
#> 1  1      Eye color            Brown            1665
#> 2  1      Eye color      Brown-green            1665
#> 3  1      Eye color       Blue-green            1665
#> 4  1      Eye color        Blue-grey            1665
#> 5  1      Eye color            Green            1665
#> 6  1      Eye color             Blue            1665
```

Output a list, then call the characteristic of interest by 'id' or 'characteristic'


```r
datalist <- allphenotypes()
```

Get a list of all characteristics you can call


```r
names(datalist)[1:10]
#>  [1] "Eye color"                        "Lactose intolerance"              "Handedness"                      
#>  [4] "white skin"                       "Ability to find a bug in openSNP" "Beard Color"                     
#>  [7] "Hair Color"                       "Ability to Tan"                   "Height"                          
#> [10] "Hair Type"
```

Get data.frame for _ADHD_


```r
datalist[["ADHD"]]
#>    id characteristic                                                                                             known_variations
#> 1  29           ADHD                                                                                                        False
#> 2  29           ADHD                                                                                                         True
#> 3  29           ADHD                                                                               Undiagnosed, but probably true
#> 4  29           ADHD                                                                                                           No
#> 5  29           ADHD                                                                                                          Yes
#> 6  29           ADHD                                                                                                Not diagnosed
#> 7  29           ADHD                                                                  Diagnosed as not having but with some signs
#> 8  29           ADHD                                                                                                  Mthfr c677t
#> 9  29           ADHD                                                                                                    Rs1801260
#> 10 29           ADHD                                                                                                  Adult onset
#> 11 29           ADHD                                                                   Diagnosed as "other hyperkinetic disorder"
#> 12 29           ADHD                                                                                 Blonde, european, green eyes
#> 13 29           ADHD                                                                                                      Extreme
#> 14 29           ADHD Diagnosed as hyperactive type, though it is my belief that adhd is simply a normal trait such as eye color. 
#>    number_of_users
#> 1              325
#> 2              325
#> 3              325
#> 4              325
#> 5              325
#> 6              325
#> 7              325
#> 8              325
#> 9              325
#> 10             325
#> 11             325
#> 12             325
#> 13             325
#> 14             325
```

Get data.frame for _mouth size_ and _SAT Writing_


```r
datalist[c("mouth size","SAT Writing")]
#> $`mouth size`
#>    id characteristic     known_variations number_of_users
#> 1 120     mouth size               Medium             202
#> 2 120     mouth size                Small             202
#> 3 120     mouth size                Large             202
#> 4 120     mouth size Slightly wide mouth              202
#> 
#> $`SAT Writing`
#>    id characteristic                                        known_variations number_of_users
#> 1  41    SAT Writing                                                     750             110
#> 2  41    SAT Writing                                      Tested before 2005             110
#> 3  41    SAT Writing                                                     800             110
#> 4  41    SAT Writing                                     Country with no sat             110
#> 5  41    SAT Writing                                                     N/a             110
#> 6  41    SAT Writing                                 Never & have ba & above             110
#> 7  41    SAT Writing                                                     720             110
#> 8  41    SAT Writing                         Did well - don't remember score             110
#> 9  41    SAT Writing                                                     511             110
#> 10 41    SAT Writing                                                     760             110
#> 11 41    SAT Writing                                                     780             110
#> 12 41    SAT Writing                                                     700             110
#> 13 41    SAT Writing Not part of sat when i took test in august 1967 at uiuc             110
#> 14 41    SAT Writing                                 Not part of sat in 1961             110
#> 15 41    SAT Writing                                                     620             110
#> 16 41    SAT Writing                                                     560             110
```

### Annotations

Get just the metadata


```r
annotations(snp = 'rs7903146', output = 'metadata')
#>          .id        V1
#> 1       name rs7903146
#> 2 chromosome        10
#> 3   position 112998590
```

Just from PLOS journals


```r
annotations(snp = 'rs7903146', output = 'plos')[c(1:2),]
#>                   author
#> 1        Maggie C. Y. Ng
#> 2 André Gustavo P. Sousa
#>                                                                                                                                      title
#> 1 Meta-Analysis of Genome-Wide Association Studies in African Americans Provides Insights into the Genetic Architecture of Type 2 Diabetes
#> 2                                  Genetic Variants of Diabetes Risk and Incident Cardiovascular Events in Chronic Coronary Artery Disease
#>           publication_date number_of_readers                                          url                          doi
#> 1 2014-08-07T00:00:00.000Z             11650 https://doi.org/10.1371/journal.pgen.1004517 10.1371/journal.pgen.1004517
#> 2 2011-01-20T00:00:00.000Z              2482 https://doi.org/10.1371/journal.pone.0016341 10.1371/journal.pone.0016341
```

Just from SNPedia


```r
annotations(snp = 'rs7903146', output = 'snpedia')
#>                                               url                                                          summary
#> 1 http://www.snpedia.com/index.php/Rs7903146(C;C) Normal (lower) risk of Type 2 Diabetes and Gestational Diabetes.
#> 2 http://www.snpedia.com/index.php/Rs7903146(C;T)     1.4x increased risk for diabetes (and perhaps colon cancer).
#> 3 http://www.snpedia.com/index.php/Rs7903146(T;T)                            2x increased risk for Type-2 diabetes
```

Get all annotations


```r
annotations(snp = 'rs7903146', output = 'all')[1:5,]
#>        .id              author
#> 1 mendeley           T E Meyer
#> 2 mendeley      Camilla Cervin
#> 3 mendeley Nicholette D Palmer
#> 4 mendeley      Ashis K Mondal
#> 5 mendeley        Julian Munoz
#>                                                                                                                                title
#> 1                                                Diabetes genes and prostate cancer in the Atherosclerosis Risk in Communities study
#> 2                                                        Diabetes in Adults , Type 1 Diabetes , and Type 2 Diabetes GENETICS OF LADA
#> 3                                Association of TCF7L2 gene polymorphisms with reduced acute insulin response in Hispanic Americans.
#> 4                  Genotype and tissue-specific effects on alternative splicing of the transcription factor 7-like 2 gene in humans.
#> 5 Polymorphism in the transcription factor 7-like 2 (TCF7L2) gene is associated with reduced insulin secretion in nondiabetic women.
#>   publication_year number_of_readers open_access
#> 1             2010                 3        TRUE
#> 2             2008                 2       FALSE
#> 3             2008                 8       FALSE
#> 4             2010                13        TRUE
#> 5             2006                10        TRUE
#>                                                                                                                                      url
#> 1                              http://www.mendeley.com/research/diabetes-genes-prostate-cancer-atherosclerosis-risk-communities-study-4/
#> 2                                        http://www.mendeley.com/research/diabetes-adults-type-1-diabetes-type-2-diabetes-genetics-lada/
#> 3              http://www.mendeley.com/research/association-tcf7l2-gene-polymorphisms-reduced-acute-insulin-response-hispanic-americans/
#> 4        http://www.mendeley.com/research/genotype-tissuespecific-effects-alternative-splicing-transcription-factor-7like-2-gene-humans/
#> 5 http://www.mendeley.com/research/polymorphism-transcription-factor-7like-2-tcf7l2-gene-associated-reduced-insulin-secretion-nondiabet/
#>                                              doi publication_date summary first_author pubmed_link journal trait pvalue
#> 1 19/2/558 [pii]\\r10.1158/1055-9965.EPI-09-0902             <NA>    <NA>         <NA>        <NA>    <NA>  <NA>     NA
#> 2                         10.2337/db07-0299.Leif             <NA>    <NA>         <NA>        <NA>    <NA>  <NA>     NA
#> 3                           10.1210/jc.2007-1225             <NA>    <NA>         <NA>        <NA>    <NA>  <NA>     NA
#> 4                           10.1210/jc.2009-2064             <NA>    <NA>         <NA>        <NA>    <NA>  <NA>     NA
#> 5                              10.2337/db06-0574             <NA>    <NA>         <NA>        <NA>    <NA>  <NA>     NA
#>   pvalue_description confidence_interval
#> 1               <NA>                <NA>
#> 2               <NA>                <NA>
#> 3               <NA>                <NA>
#> 4               <NA>                <NA>
#> 5               <NA>                <NA>
```

### Download

Download genotype data for a user from 23andme or other repo. (not evaluated in this example)


```r
data <- users(df=TRUE)
head(data[[1]])
fetch_genotypes(url = data[[1]][1,"genotypes.download_url"], rows=15)
```

### Genotype user data

Genotype data for one or multiple users


```r
genotypes(snp='rs9939609', userid=1)
#> $snp
#> $snp$name
#> [1] "rs9939609"
#> 
#> $snp$chromosome
#> [1] "16"
#> 
#> $snp$position
#> [1] "53786615"
#> 
#> 
#> $user
#> $user$name
#> [1] "Bastian Greshake Tzovaras"
#> 
#> $user$id
#> [1] 1
#> 
#> $user$genotypes
#> $user$genotypes[[1]]
#> $user$genotypes[[1]]$genotype_id
#> [1] 9
#> 
#> $user$genotypes[[1]]$local_genotype
#> [1] "AT"
```



```r
genotypes('rs9939609', userid='1,6,8', df=TRUE)
#>    snp_name snp_chromosome snp_position                 user_name user_id genotype_id genotype
#> 1 rs9939609             16     53786615 Bastian Greshake Tzovaras       1           9       AT
#> 2 rs9939609             16     53786615              Nash Parovoz       6           5       AT
#> 3 rs9939609             16     53786615         Samantha B. Clark       8           2       TT
```



```r
genotypes('rs9939609', userid='1-2', df=FALSE)
#> [[1]]
#> [[1]]$snp
#> [[1]]$snp$name
#> [1] "rs9939609"
#> 
#> [[1]]$snp$chromosome
#> [1] "16"
#> 
#> [[1]]$snp$position
#> [1] "53786615"
#> 
#> 
#> [[1]]$user
#> [[1]]$user$name
#> [1] "Bastian Greshake Tzovaras"
#> 
#> [[1]]$user$id
#> [1] 1
#> 
#> [[1]]$user$genotypes
#> [[1]]$user$genotypes[[1]]
#> [[1]]$user$genotypes[[1]]$genotype_id
#> [1] 9
#> 
#> [[1]]$user$genotypes[[1]]$local_genotype
#> [1] "AT"
#> 
#> 
#> 
#> 
#> 
#> [[2]]
#> [[2]]$snp
#> [[2]]$snp$name
#> [1] "rs9939609"
#> 
#> [[2]]$snp$chromosome
#> [1] "16"
#> 
#> [[2]]$snp$position
#> [1] "53786615"
#> 
#> 
#> [[2]]$user
#> [[2]]$user$name
#> [1] "Senficon"
#> 
#> [[2]]$user$id
#> [1] 2
#> 
#> [[2]]$user$genotypes
#> list()
```

### Phenotype user data

Get phenotype data for one or multiple users



```r
phenotypes(userid=1)$phenotypes[1:3]
#> $`Caffeine dependence`
#> $`Caffeine dependence`$phenotype_id
#> [1] 538
#> 
#> $`Caffeine dependence`$variation
#> [1] "No"
#> 
#> 
#> $`hair on ear`
#> $`hair on ear`$phenotype_id
#> [1] 254
#> 
#> $`hair on ear`$variation
#> [1] "No"
#> 
#> 
#> $`Third Nipple`
#> $`Third Nipple`$phenotype_id
#> [1] 259
#> 
#> $`Third Nipple`$variation
#> [1] "None"
```


```r
phenotypes(userid='1,6,8', df=TRUE)[[1]][1:10,]
#>                                phenotype phenotypeID                                                           variation
#> 1                    Caffeine dependence         538                                                                  No
#> 2                            hair on ear         254                                                                  No
#> 3                           Third Nipple         259                                                                None
#> 4                             Alcoholism         485                                                                None
#> 5         Alcohol Consumption (per week)         484                                                                   0
#> 6  Allergy to artificial grape flavoring         352                                                                  No
#> 7                       inverted nipples         583                                                                None
#> 8    Do you prefer python, matlab, or R?         585                                                          Python & R
#> 9                      Political Compass         276 Economic Left/Right: -8.88  Social Libertarian/Authoritarian: -9.49
#> 10               Sweat eating spicy food         219                                                                 Yes
```



```r
out <- phenotypes(userid='1-8', df=TRUE)
lapply(out, head)
#> $`Bastian Greshake Tzovaras`
#>                               phenotype phenotypeID variation
#> 1                   Caffeine dependence         538        No
#> 2                           hair on ear         254        No
#> 3                          Third Nipple         259      None
#> 4                            Alcoholism         485      None
#> 5        Alcohol Consumption (per week)         484         0
#> 6 Allergy to artificial grape flavoring         352        No
#> 
#> $Senficon
#>   phenotype phenotypeID variation
#> 1   no data     no data   no data
#> 
#> $`no info on user_3`
#>   phenotype phenotypeID variation
#> 1   no data     no data   no data
#> 
#> $`no info on user_4`
#>   phenotype phenotypeID variation
#> 1   no data     no data   no data
#> 
#> $`no info on user_5`
#>   phenotype phenotypeID variation
#> 1   no data     no data   no data
#> 
#> $`Nash Parovoz`
#>                          phenotype phenotypeID        variation
#> 1         Y-DNA Haplogroup (ISOGG)         150        J-FGC5206
#> 2  The Dress: Perception of colour         338   White and gold
#> 3           Number of wisdom teeth          57                4
#> 4 Ability to find a bug in openSNP           5   extremely high
#> 5              Lactose intolerance           2 lactose-tolerant
#> 6                       white skin           4        Caucasian
#> 
#> $`no info on user_7`
#>   phenotype phenotypeID variation
#> 1   no data     no data   no data
#> 
#> $`Samantha B. Clark`
#>                             phenotype phenotypeID           variation
#> 1                            Gambling         539                  No
#> 2                 Caffeine dependence         538                  No
#> 3            Dietary supplements used         534                 b12
#> 4                                Diet         533 Vegan / plant-based
#> 5                   Tooth sensitivity         532         Sweet, cold
#> 6 OCD - Obsessive-Compulsive Disorder         555                  No
```

### All known variations

Get all known variations and all users sharing that phenotype for one phenotype(-ID).


```r
phenotypes_byid(phenotypeid=12, return_ = 'desc')
#> $id
#> [1] 12
#> 
#> $characteristic
#> [1] "Beard Color"
#> 
#> $description
#> [1] "coloration of facial hair"
```



```r
phenotypes_byid(phenotypeid=12, return_ = 'knownvars')
#> $known_variations
#> $known_variations[[1]]
#> [1] "Red"
#> 
#> $known_variations[[2]]
#> [1] "Blonde"
#> 
#> $known_variations[[3]]
#> [1] "Red-brown"
#> 
#> $known_variations[[4]]
#> [1] "Red-blonde-brown-black(in diferent parts i have different color,for example near the lips blond-red"
#> 
#> $known_variations[[5]]
#> [1] "No beard-female"
#> 
#> $known_variations[[6]]
#> [1] "Brown-black"
#> 
#> $known_variations[[7]]
#> [1] "Blonde-brown"
#> 
#> $known_variations[[8]]
#> [1] "Black"
#> 
#> $known_variations[[9]]
#> [1] "Dark brown with minor blondish-red"
#> 
#> $known_variations[[10]]
#> [1] "Brown-grey"
#> 
#> $known_variations[[11]]
#> [1] "Red-blonde-brown-black"
#> 
#> $known_variations[[12]]
#> [1] "Blond-brown"
#> 
#> $known_variations[[13]]
#> [1] "Brown, some red"
#> 
#> $known_variations[[14]]
#> [1] "Brown"
#> 
#> $known_variations[[15]]
#> [1] "Brown-gray"
#> 
#> $known_variations[[16]]
#> [1] "Never had a beard"
#> 
#> $known_variations[[17]]
#> [1] "I'm a woman"
#> 
#> $known_variations[[18]]
#> [1] "Black-brown-blonde"
#> 
#> $known_variations[[19]]
#> [1] "Was red-brown now mixed with gray,"
#> 
#> $known_variations[[20]]
#> [1] "Red-blonde-brown"
#> 
#> $known_variations[[21]]
#> [1] "Dark brown w/few blonde & red hairs"
#> 
#> $known_variations[[22]]
#> [1] "Dark blonde with red and light blonde on goatee area."
#> 
#> $known_variations[[23]]
#> [1] "Black with few red hairs"
#> 
#> $known_variations[[24]]
#> [1] "Black, graying"
#> 
#> $known_variations[[25]]
#> [1] "Red, moustache still is, beard mostly white"
#> 
#> $known_variations[[26]]
#> [1] "Blonde/brown-some black-and red on chin-all starting to gray"
#> 
#> $known_variations[[27]]
#> [1] "Dark brown"
#> 
#> $known_variations[[28]]
#> [1] "Every possible color, most hair shafts have more than one color at different points along the shaft"
#> 
#> $known_variations[[29]]
#> [1] "Black with few white hairs"
#> 
#> $known_variations[[30]]
#> [1] "Brown ginger"
#> 
#> $known_variations[[31]]
#> [1] "Dark blonde"
#> 
#> $known_variations[[32]]
#> [1] "Black - going white due to age"
#> 
#> $known_variations[[33]]
#> [1] "N/a"
```



```r
phenotypes_byid(phenotypeid=12, return_ = 'users')[1:10,]
#>    user_id                                                                                           variation
#> 1       22                                                                                                 Red
#> 2        1                                                                                              Blonde
#> 3       26                                                                                           red-brown
#> 4       10 Red-Blonde-Brown-Black(in diferent parts i have different color,for example near the lips blond-red
#> 5       14                                                                                     No beard-female
#> 6       42                                                                                         Brown-black
#> 7       45 Red-Blonde-Brown-Black(in diferent parts i have different color,for example near the lips blond-red
#> 8       16                                                                                        blonde-brown
#> 9        8                                                                                     No beard-female
#> 10     661                                                                                         Brown-black
```


## NCBI SNP data

### dbSNP

Query NCBI's [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) for information on a set of SNPs. 

An example with four markers, where one has been merged, and one has been withdrawn from NCBI.


```r
snps <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
(dbsnp_info <- ncbi_snp_query(snps))
#> # A tibble: 4 × 16
#>   query  chromosome     bp class rsid  gene  alleles ancestral_allele variation_allele seqname hgvs  assembly ref_seq minor    maf
#>   <chr>  <chr>       <dbl> <chr> <chr> <chr> <chr>   <chr>            <chr>            <chr>   <chr> <chr>    <chr>   <chr>  <dbl>
#> 1 rs332  7          1.18e8 del   rs12… "CFT… TTT, d… TTT              delTTT           NC_000… NC_0… GRCh38.… <NA>    <NA>  NA    
#> 2 rs420… 1          4.03e7 snv   rs42… ""    A,C,G,T A                C,G,T            NC_000… NC_0… GRCh38.… <NA>    <NA>  NA    
#> 3 rs183… 5          1.11e8 snv   rs18… ""    T,C     T                C                NC_000… NC_0… GRCh38.… T       C      0.726
#> 4 rs120… 9          4.18e7 snv   rs12… ""    T,A,C   T                A,C              NC_000… NC_0… GRCh38.… <NA>    <NA>  NA    
#> # … with 1 more variable: maf_population <list>
```

The maf column contains the minor allele frequency from the GnomAD database (if available). All population specific allele frequencies can be accessed through the column `maf_population` which returns a list.

```r
dbsnp_info$maf_population
#> [[1]]
#>   study ref_seq Minor MAF
#> 1            NA    NA  NA
#> 
#> [[2]]
#>             study ref_seq Minor       MAF
#> 1          ALSPAC       A     C 0.8227815
#> 2        Estonian       A     C 0.7895089
#> 3       GENOME_DK       A     C 0.8750000
#> 4            GoNL       A     C 0.8266533
#> 5          KOREAN       A     C 0.9658703
#> 6  NorthernSweden       A     C 0.8183333
#> 7          Qatari       A     C 0.8379630
#> 8        SGDP_PRJ       A     C 0.9175824
#> 9        Siberian       A     C 0.8333333
#> 10          TOMMO       A     C 0.9589499
#> 11         TOPMED       A     C 0.8765689
#> 12         TOPMED       A     C 0.8767313
#> 13        TWINSUK       A     C 0.8193096
#> 14     Vietnamese       A     C 0.9952830
#> 15  dbGaP_PopFreq       A     C 0.7991653
#> 16         KOREAN       A     G 0.0000000
#> 17         KOREAN       A     T 0.0000000
#> 18  dbGaP_PopFreq       A     T 0.0000000
#> 
#> [[3]]
#>             study ref_seq Minor       MAF
#> 1     1000Genomes       T     C 0.6178115
#> 2          ALSPAC       T     C 0.7477945
#> 3       Daghestan       T     C 0.6856128
#> 4        Estonian       T     C 0.7037946
#> 5       GENOME_DK       T     C 0.7250000
#> 6          GnomAD       T     C 0.7257767
#> 7            GoNL       T     C 0.7274549
#> 8   HGDP_Stanford       T     C 0.6602687
#> 9          HapMap       T     C 0.6054025
#> 10         KOREAN       T     C 0.3969283
#> 11        Korea1K       T     C 0.3733624
#> 12 NorthernSweden       T     C 0.6850000
#> 13     PAGE_STUDY       T     C 0.6673868
#> 14     PRJEB36033       T     C 1.0000000
#> 15     PRJEB37584       T     C 0.4141414
#> 16         Qatari       T     C 0.7824074
#> 17       SGDP_PRJ       T     C 0.7670940
#> 18       Siberian       T     C 0.7826087
#> 19          TOMMO       T     C 0.3405728
#> 20         TOPMED       T     C 0.7110490
#> 21         TOPMED       T     C 0.7196758
#> 22        TWINSUK       T     C 0.7437972
#> 23     Vietnamese       T     C 0.4074074
#> 24  dbGaP_PopFreq       T     C 0.7300999
#> 
#> [[4]]
#>           study ref_seq Minor MAF
#> 1 dbGaP_PopFreq       T     A   0
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-defunct.R
\name{NCBI_snp_query2}
\alias{NCBI_snp_query2}
\title{This function is defunct.}
\usage{
NCBI_snp_query2(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allphenotypes.R
\name{allphenotypes}
\alias{allphenotypes}
\title{Get all openSNP phenotypes, their variations, and how many users have data
available for a given phenotype.}
\usage{
allphenotypes(df = FALSE, ...)
}
\arguments{
\item{df}{Return a data.frame of all data. The column known_variations
can take multiple values, so the other columns id, characteristic, and
number_of_users are replicated in the data.frame. Default: \code{FALSE}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
data.frame of results, or list if \code{df=FALSE}
}
\description{
Either return data.frame with all results, or output a list, then call
the charicteristic by id (parameter = "id") or name (parameter =
"characteristic").
}
\examples{
\dontrun{
# Get all data
allphenotypes(df = TRUE)

# Output a list, then call the characterisitc of interest by 'id' or
# 'characteristic'
datalist <- allphenotypes()
names(datalist) # get list of all characteristics you can call
datalist[["ADHD"]] # get data.frame for 'ADHD'
datalist[c("mouth size","SAT Writing")] # get data.frame for 'ADHD'
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{swap}
\alias{swap}
\title{Swap Elements in a Vector}
\usage{
swap(vec, from, to = names(from), ...)
}
\arguments{
\item{vec}{A character vector, or vector coercable to character.}

\item{from}{A vector of elements to map from.}

\item{to}{A vector of elements to map to.}

\item{...}{Optional arguments passed to \code{\link[=match]{match()}}}
}
\description{
Swap Elements in a Vector
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_users.R
\docType{data}
\name{rsnpsCache}
\alias{rsnpsCache}
\title{rsnps environment}
\format{
An object of class \code{environment} of length 0.
}
\usage{
rsnpsCache
}
\description{
rsnps environment
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_users.R
\name{read_users}
\alias{read_users}
\title{Read in openSNP user files from local storage.}
\usage{
read_users(name = NULL, id = NULL, path = NULL, ...)
}
\arguments{
\item{name}{User name}

\item{id}{User id}

\item{path}{Path to file to read from.}

\item{...}{Parameters passed on to \code{\link[=read.table]{read.table()}}}
}
\value{
A data.frame.
}
\description{
Beware, these tables can be large. Check your RAM before executing. Or
possibly read in a subset of the data. This function reads in the
whole kitten kaboodle.
}
\details{
If you specify a name or id, this function reads environment variables
written in the function download_users, and then searches against those
variables for the path to the file saved. Alternatively, you can supply
the path.
}
\examples{
\dontrun{
# dat <- read_users(name = "kevinmcc")
# head(dat)
# dat <- read_users(id = 285)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDSearch.R
\name{ld_search}
\alias{ld_search}
\title{This function is defunct.}
\usage{
ld_search(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_snp_api.R
\name{get_frequency}
\alias{get_frequency}
\title{Internal function to get the frequency of the variants from
different studies.}
\usage{
get_frequency(Class, primary_info)
}
\arguments{
\item{Class}{What kind of variant is the rsid. Accepted options are "snv", "snp" and "delins".}

\item{primary_info}{refsnp entry read in JSON format}
}
\description{
Internal function to get the frequency of the variants from
different studies.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tryget}
\alias{tryget}
\title{Tryget}
\usage{
tryget(x)
}
\description{
Tryget
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgensnp.R
\name{allgensnp}
\alias{allgensnp}
\title{Get openSNP genotype data for all users at a particular snp.}
\usage{
allgensnp(snp = NA, ...)
}
\arguments{
\item{snp}{(character) A SNP name}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
data.frame of genotypes for all users at a certain SNP
}
\description{
Get openSNP genotype data for all users at a particular snp.
}
\examples{
\dontrun{
x <- allgensnp(snp = 'rs7412')
head(x)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_snp_api.R
\name{ncbi_snp_query}
\alias{ncbi_snp_query}
\title{Query NCBI's refSNP for information on a set of SNPs via the API}
\usage{
ncbi_snp_query(snps)
}
\arguments{
\item{snps}{(character) A vector of SNPs (rs numbers).}
}
\value{
A dataframe with columns:
\itemize{
\item query: The rs ID that was queried.
\item chromosome: The chromosome that the marker lies on.
\item bp: The chromosomal position, in base pairs, of the marker,
as aligned with the current genome used by dbSNP. we add 1 to the base
pair position in the BP column in the output data.frame to agree with
what the dbSNP website has.
\item rsid: Reference SNP cluster ID. If the rs ID queried
has been merged, the up-to-date name of the ID is returned here, and
a warning is issued.
\item class: The rsid's 'class'. See
\url{https://www.ncbi.nlm.nih.gov/projects/SNP/snp_legend.cgi?legend=snpClass}
for more details.
\item gene: If the rsid lies within a gene (either within the exon
or introns of a gene), the name of that gene is returned here; otherwise,
\code{NA}. Note that
the gene may not be returned if the rsid lies too far upstream or downstream
of the particular gene of interest.
\item alleles: The alleles associated with the SNP if it is a
SNV; otherwise, if it is an INDEL, microsatellite, or other kind of
polymorphism the relevant information will be available here.
\item minor: The allele for which the MAF is computed,
given it is an SNV; otherwise, \code{NA}.
\item maf: The minor allele frequency of the SNP, given it is an SNV.
This is drawn from the current global reference population used by NCBI (GnomAD).
\item ancestral_allele: allele as described in the current assembly
\item variation_allele: difference to the current assembly
\item seqname - Chromosome RefSeq reference.
\item hgvs -  full hgvs notation for variant
\item assembly - which assembly was used for the annotations
\item ref_seq - sequence in reference assembly
\item maf_population - dataframe of all minor allele frequencies reported, with columns study,
reference allele, alternative allele (minor) and minor allele frequency.
}
}
\description{
This function queries NCBI's refSNP for information related to the latest
dbSNP build and latest reference genome for information on the vector
of snps submitted.
}
\details{
This function currently pulling data for Assembly 38 - in particular
note that if you think the BP position is wrong, that you may be
hoping for the BP position for a different Assembly.

Note that you are limited in the to a max of one query per second
and concurrent queries are not allowed.
If users want to set curl options when querying for the SNPs they can do so by using
httr::set_config/httr::with_config
}
\examples{
\dontrun{
## an example with both merged SNPs, non-SNV SNPs, regular SNPs,
## SNPs not found, microsatellite
SNPs <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
ncbi_snp_query(SNPs)
# ncbi_snp_query("123456") ##invalid: must prefix with 'rs'
ncbi_snp_query("rs420358")
ncbi_snp_query("rs332") # warning that its merged into another, try that
ncbi_snp_query("rs121909001")
ncbi_snp_query("rs1837253")
ncbi_snp_query("rs1209415715")
ncbi_snp_query("rs111068718")
ncbi_snp_query(snps='rs9970807')

ncbi_snp_query("rs121909001")
ncbi_snp_query("rs121909001", verbose = TRUE)
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/projects/SNP/}

\url{https://pubmed.ncbi.nlm.nih.gov/31738401/} SPDI model
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_snp_api.R
\name{get_placements}
\alias{get_placements}
\title{Internal function to get the position, alleles, assembly, hgvs notation}
\usage{
get_placements(primary_info)
}
\arguments{
\item{primary_info}{refsnp entry read in JSON format}
}
\description{
Internal function to get the position, alleles, assembly, hgvs notation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/users.R
\name{users}
\alias{users}
\title{Get openSNP users.}
\usage{
users(df = FALSE, ...)
}
\arguments{
\item{df}{Return data.frame (\code{TRUE}) or not (\code{FALSE}). Default: \code{FALSE}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
List of openSNP users, their ID numbers, and XX if available.
}
\description{
Get openSNP users.
}
\examples{
\dontrun{
# just the list
data <- users(df = FALSE)
data

# get a data.frame of the users data
data <- users(df = TRUE)
data[[1]] # users with links to genome data
data[[2]] # users without links to genome data
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_users.R
\name{download_users}
\alias{download_users}
\title{Download openSNP user files.}
\usage{
download_users(name = NULL, id = NULL, dir = "~/", ...)
}
\arguments{
\item{name}{User name}

\item{id}{User id}

\item{dir}{Directory to save file to}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
File downloaded to directory you specify (or default), nothing
returned in R.
}
\description{
Download openSNP user files.
}
\examples{
\dontrun{
# Download a single user file, by id
download_users(id = 14)

# Download a single user file, by user name
download_users(name = 'kevinmcc')

# Download many user files
lapply(c(14,22), function(x) download_users(id=x))
read_users(id=14, nrows=5)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-defunct.R
\name{rsnps-defunct}
\alias{rsnps-defunct}
\title{Defunct functions in rsnps}
\description{
\itemize{
\item \code{LDSearch()}: Function name changed to \link{ld_search}
\item \code{ld_search()}: The Broad Institute took the service down, see
https://www.broadinstitute.org/snap/snap
\item \code{NCBI_snp_query()}: Function name changed to \link{ncbi_snp_query}
\item \code{NCBI_snp_query2()}: Function name changed to \link{ncbi_snp_query}
\item \code{ncbi_snp_summary()}: Function name changed to \link{ncbi_snp_query}
\item \code{ncbi_snp_query2()}: Function name changed to \link{ncbi_snp_query}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsnps-package.R
\docType{package}
\name{rsnps-package}
\alias{rsnps-package}
\alias{rsnps}
\title{Get SNP (Single-Nucleotide Polymorphism) Data on the Web}
\description{
This package gives you access to data from OpenSNP (https://opensnp.org)
via their API (https://opensnp.org/faq#api) and NCBI's dbSNP SNP database
(https://www.ncbi.nlm.nih.gov/snp).
}
\section{NCBI Authenication}{

This applies the function \code{\link[=ncbi_snp_query]{ncbi_snp_query()}}:

You can optionally use an API key, if you do it will
allow higher rate limits (more requests per time period)

If you don't have an NCBI API key, get one at
https://www.ncbi.nlm.nih.gov/account/

Create your key from your account. After generating your key
set an environment variable as \code{ENTREZ_KEY} in .Renviron.

\code{ENTREZ_KEY='youractualkeynotthisstring'}

You can optionally pass in your API key to the key parameter in NCBI
functions in this package. However, it's much better from a security
perspective to set an environment variable.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Kevin Ushey \email{kevinushey@gmail.com}

Hao Zhu \email{haozhu233@gmail.com}

Sina Rüeger \email{sina.rueeger@gmail.com}

Julia Gustavsen \email{j.gustavsen@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbi_snp_api.R
\name{get_gene_names}
\alias{get_gene_names}
\title{Internal function to get gene names.}
\usage{
get_gene_names(primary_info)
}
\arguments{
\item{primary_info}{refsnp entry read in JSON format}
}
\description{
If multiple gene names are encountered they are collapsed with a
"/".
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypes.R
\name{genotypes}
\alias{genotypes}
\title{Get openSNP genotype data for one or multiple users.}
\usage{
genotypes(snp = NA, userid = NA, df = FALSE, ...)
}
\arguments{
\item{snp}{SNP name.}

\item{userid}{ID of openSNP user.}

\item{df}{Return data.frame (\code{TRUE}) or not (\code{FALSE}). Default: \code{FALSE}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}]}
}
\value{
List (or data.frame) of genotypes for specified user(s) at a
certain SNP.
}
\description{
Get openSNP genotype data for one or multiple users.
}
\examples{
\dontrun{
genotypes(snp='rs9939609', userid=1)
genotypes('rs9939609', userid='1,6,8', df=TRUE)
genotypes('rs9939609', userid='1-2', df=FALSE)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{flip}
\alias{flip}
\title{Flip Genotypes}
\usage{
flip(SNP, sep = "", outSep = sep)
}
\arguments{
\item{SNP}{A vector of genotypes for a particular locus.}

\item{sep}{The separator between each allele.}

\item{outSep}{The output separator to use.}
}
\description{
Given a set of genotypes, flip them.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{annotations}
\alias{annotations}
\title{Get all openSNP phenotypes, their variations, and how many users have data
available for a given phenotype.}
\usage{
annotations(
  snp = NA,
  output = c("all", "plos", "mendeley", "snpedia", "metadata"),
  ...
)
}
\arguments{
\item{snp}{SNP name.}

\item{output}{Name the source or sources you want annotations from (options
are: 'plos', 'mendeley', 'snpedia', 'metadata'). 'metadata' gives the
metadata for the response.}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
data.frame of results
}
\description{
Either return data.frame with all results, or output a list, then call
the charicteristic by id (parameter = "id") or
name (parameter = "characteristic").
}
\examples{
\dontrun{
# Get all data
## get just the metadata
annotations(snp = 'rs7903146', output = 'metadata')

## just from plos
annotations(snp = 'rs7903146', output = 'plos') 

## just from snpedia
annotations(snp = 'rs7903146', output = 'snpedia')

## get all annotations
annotations(snp = 'rs7903146', output = 'all') 
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-defunct.R
\name{ncbi_snp_summary}
\alias{ncbi_snp_summary}
\alias{ncbi_snp_query2}
\alias{NCBI_snp_query}
\title{This function is defunct.}
\usage{
ncbi_snp_summary(...)

ncbi_snp_query2(...)

NCBI_snp_query(...)
}
\description{
This function is defunct.

This function is defunct.

This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_genotypes.R
\name{fetch_genotypes}
\alias{fetch_genotypes}
\title{Download openSNP genotype data for a user}
\usage{
fetch_genotypes(url, rows = 100, filepath = NULL, quiet = TRUE, ...)
}
\arguments{
\item{url}{(character) URL for the download. See example below of
function use.}

\item{rows}{(integer) Number of rows to read in. Useful for getting a
glimpse of  the data. Negative and other invalid values are ignored,
giving back all data.  Default: 100}

\item{filepath}{(character) If none is given the file is saved to a
temporary file,  which will be lost after your session is closed. Save
to a file if you want to  access it later.}

\item{quiet}{(logical) Should download progress be suppressed. Default:
\code{TRUE}}

\item{...}{Further args passed on to \code{\link[=download.file]{download.file()}}}
}
\value{
data.frame for a single user, with four columns:
\itemize{
\item rsid (character)
\item chromosome (integer)
\item position (integer)
\item genotype (character)
}
}
\description{
Download openSNP genotype data for a user
}
\details{
Beware, not setting the rows parameter means that you download
the entire file, which can be large (e.g., 15MB), and so take a while
to download depending on your connection speed. Therefore, rows is set to
10 by default to sort of protect the user.

Internally, we use \code{\link[=download.file]{download.file()}} to download each file, then
\code{\link[=read.table]{read.table()}} to read the file to a data.frame.
}
\examples{
\dontrun{
# get a data.frame of the users data
data <- users(df = TRUE)
head( data[[1]] ) # users with links to genome data
mydata <- fetch_genotypes(url = data[[1]][1,"genotypes.download_url"], 
  file="~/myfile.txt")

# see some data right away
mydata

# Or read in data later separately
read.table("~/myfile.txt", nrows=10)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{split_to_df}
\alias{split_to_df}
\title{Split a Vector of Strings Following a Regular Structure}
\usage{
split_to_df(x, sep, fixed = FALSE, perl = TRUE, useBytes = FALSE, names = NULL)
}
\arguments{
\item{x}{a vector of strings.}

\item{sep}{the delimiter / \link{regex} you wish to split your strings on.}

\item{fixed}{logical. If \code{TRUE}, we match \code{sep} exactly;
otherwise, we use regular expressions. Has priority over \code{perl}.}

\item{perl}{logical. Should perl-compatible regexps be used?}

\item{useBytes}{logical. If \code{TRUE}, matching is done byte-by-byte rather than
character-by-character.}

\item{names}{optional: a vector of names to pass to the returned \code{data.frame}.}
}
\description{
This function takes a vector of strings following a regular
structure, and converts that vector into a \code{data.frame}, split
on that delimiter. A nice wrapper to \code{\link[=strsplit]{strsplit()}}, essentially
\itemize{
\item the primary bonus is the automatic coersion to a \code{data.frame}.
}
}
\seealso{
\code{\link[=strsplit]{strsplit()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenotypes.R
\name{phenotypes}
\alias{phenotypes}
\title{Get openSNP phenotype data for one or multiple users.}
\usage{
phenotypes(userid = NA, df = FALSE, ...)
}
\arguments{
\item{userid}{ID of openSNP user.}

\item{df}{Return data.frame (\code{TRUE}) or not (\code{FALSE}). Default: \code{FALSE}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
List of phenotypes for specified user(s).
}
\description{
Get openSNP phenotype data for one or multiple users.
}
\examples{
\dontrun{
phenotypes(userid=1)
phenotypes(userid='1,6,8', df=TRUE)
phenotypes(userid='1-8', df=TRUE)

# coerce to data.frame
library(plyr)
df <- ldply(phenotypes(userid='1-8', df=TRUE))
head(df); tail(df)

# pass on curl options
phenotypes(1, verbose = TRUE)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes_byid}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-defunct.R
\name{LDSearch}
\alias{LDSearch}
\title{This function is defunct.}
\usage{
LDSearch(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenotypes_byid.R
\name{phenotypes_byid}
\alias{phenotypes_byid}
\title{Get all openSNP known variations and all users sharing that phenotype for
one phenotype(-ID).}
\usage{
phenotypes_byid(
  phenotypeid = NA,
  return_ = c("description", "knownvars", "users"),
  ...
)
}
\arguments{
\item{phenotypeid}{ID of openSNP phenotype.}

\item{return_}{Return data.frame (\code{TRUE}) or not (\code{FALSE}). Default: \code{FALSE}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
List of description of phenotype, list of known variants, or
data.frame of variants for each user with that phenotype.
}
\description{
Get all openSNP known variations and all users sharing that phenotype for
one phenotype(-ID).
}
\examples{
\dontrun{
phenotypes_byid(phenotypeid=12, return_ = 'desc')
phenotypes_byid(phenotypeid=12, return_ = 'knownvars')
phenotypes_byid(phenotypeid=12, return_ = 'users')

# pass on curl options
phenotypes_byid(phenotypeid=12, return_ = 'desc', verbose = TRUE)
}
}
\seealso{
Other opensnp-fxns: 
\code{\link{allgensnp}()},
\code{\link{allphenotypes}()},
\code{\link{annotations}()},
\code{\link{download_users}()},
\code{\link{fetch_genotypes}()},
\code{\link{genotypes}()},
\code{\link{phenotypes}()},
\code{\link{users}()}
}
\concept{opensnp-fxns}
