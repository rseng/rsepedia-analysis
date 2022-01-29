gistr
=======



[![cran checks](https://cranchecks.info/badges/worst/gistr)](https://cranchecks.info/pkgs/gistr)
[![R-check](https://github.com/ropensci/gistr/workflows/R-check/badge.svg)](https://github.com/ropensci/gistr/actions/)
[![codecov.io](https://codecov.io/github/ropensci/gistr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/gistr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/gistr)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/gistr)](https://cran.r-project.org/package=gistr)

`gistr` is a light interface to GitHub's gists for R.

Get started with the docs: https://docs.ropensci.org/gistr

## See also:

* [git2r](https://github.com/ropensci/git2r) an R client for the libgit2 C library by Stefan Widgren
* [gert](https://github.com/r-lib/gert) Simple git client for R by Jeroen Ooms
* [gistfo](https://github.com/MilesMcBain/gistfo) for turning your untitled RStudio tabs into gists!

## Install

Stable version from CRAN


```r
install.packages("gistr")
```

Or dev version from GitHub.


```r
remotes::install_github("ropensci/gistr")
```


```r
library("gistr")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/gistr/issues).
* License: MIT
* Get citation information for `gistr` in R doing `citation(package = 'gistr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
gistr in development
===============

### BUG FIXES

* fix `gist_auth()`: at some point `httr::oauth2.0_token` stopped returning the token in the `headers` slot; can't figure out when this change happened; fix is to get the token from a different place in the returned object; changes to `gist_auth()` to access that new path to the token (#83)


gistr 0.9.0
===============

### MINOR IMPROVEMENTS

* replace httr with crul for all but the oauth handling (#68)
* fix to internal fxn `stopstatus()` to handle correctly detecting scope header issues (#82)
* fixed old urls that have changed


gistr 0.5.0
===============

### MINOR IMPROVEMENTS

* vignette gains an example of round-tripping a data.frame to a gist then back from the gist to a data.frame (#78) (#79) thanks @jsta
* update package docs throughout to tell users to make sure to create a GitHub PAT (personal access token) with gist scope (#80)

### BUG FIXES

* fix to `gist_create()`: fail if both `files` and `code` params are `NULL` (the user should pass something in) (#72) thanks @maelle


gistr 0.4.2
===============

### NEW FEATURES

* `gist()` gains a parameter `revision` to request a specific revision of a gists. note that the returned brief print out of the gist in your console may not vary from revision to revision, but the underlying data has the correct data for the revision (#71)

### MINOR IMPROVEMENTS

* affecting all functions that create data `gist_create()`, `gist_create_git()`, `gist_create_obj()`, `update()`, `delete()`: GitHub for good reason gives a 404 when there are authentication issues. A common problem is that a user has incorrect or missing scopes. We now attempt to detect this scope problem specifically and throw a message when that happens (#70)
* toggle whether we index to a git path with `@` vs. `$` depending on `git2r` package version; for an upcoming version of `git2r` (#74)


gistr 0.4.0
===============

### MINOR IMPROVEMENTS

* Change all `dplyr::rbind_all` instances to `dplyr::bind_rows` (#69)

### BUG FIXES

* Fix to `gists()` internals for when `github.username` not set 
and user selects `what = "mineall"` - now stops with informative
message about setting `github.username` option (#66) (#67) thanks @sboysel


gistr 0.3.6
===============

### MINOR IMPROVEMENTS

* Added more tests for `as.gist()`

### BUG FIXES

* Fix to `as.gist.list()` method to not break sometimes when not all keys
returned in JSON content from github API (#63)
* Fix to `update()` to work correctly for deleting files. didn't previously
set `null`'s correctly (#64)

gistr 0.3.4
===============

### NEW FEATURES

* Gained new function `gist_create_git()` - creates gists using `git` 
instead of the GitHub Gists HTTP API. Uses the package `git2r` 
internally to do the `git` things. (#50) This function has been 
around a while, but not in the CRAN version, so a few other fixes
of note in case you're interested: (#56) (#57) (#58) (#59) (#61)

### MINOR IMPROVEMENTS

* Added new manual file `?create_gists` with details of the three different
ways to create gists, how they differ, and why there are three different
functions to create a gist. (57f13a711fb7a1514caee6a858d4cda31d614e6f)

### BUG FIXES

* Fix to `tabl()` to give back cleaner data output, returning main
metadata for each gist in a single data.frame, then forks and 
history in separate data.frame's if they exist. Makes for easier 
understanding and manipulation downstream. (#54)

gistr 0.3.0
===============

### NEW FEATURES

* Gained new function `gistr_save()` to save gist files to disk easily and optionally open them in your editor/R GUI (#47). In addition, files saved to a directory, with the dir named by the gist id (#49)
* `gist()` now accepts either a gist ID or full or partial URL's for a gist (#48)

### MINOR IMPROVEMENTS

* Can now optionally use `rmarkdown::render()` with `gist_create()` (#52)
* Explicitly import non-base R pkg functions, so importing from `utils`, `methods`, and `stats` (#53)
* Can now toggle use of `rmarkdown` package with a parameter in `gist_create()` (#52)
* Better error messages from the GitHub API (#42)

### BUG FIXES

* Fixed problem with `httr` `v1` where empty list not allowed to pass to 
the `query` parameter in `GET` (#51)

gistr 0.2.0
===============

### NEW FEATURES

* `gistr_create()` can now optionally include source file if `knit=TRUE` using the new
parameter `include_source` (#19)
* new function `gist_create_obj()` to create a gist directly from a R object, like
numeric, list, character, data.frame, matrix (#36) (#44)
* new function `gist_map()` to open a full page map in your default browser of a gist
after gist creation (#23)
* new function `tabl()` (weird function name to avoid the `table` function in base R).
This function goal is to make it easier to play with gist data. Data given back from the
GitHub API is great, but is in nested list format (after conversion from JSON) - so
is rather hard to manipulate. `tabl()` makes a data.frame from output of `gist()`,
`gists()`, `as.gist()`, and `commits()` (#25)

### MINOR IMPROVEMENTS

* `gistr_create()` works with `.Rnw` files, and example `.Rnw` file included in the package. (#20)
* Added ability in `gist_create()` to optionally include the source file passed into
the function call when `knit=TRUE` (#19)
* Added ability to inject imgur hooks into a knitted document so that images can be rendered in a gist automatically. The GitHub HTTP API doesn't allow binary uploads
(e.g., images), so the parameter `imgur_inject` uploads your images to imgur
and embeds links to the images in your document. (#33)
* Improved information on truncation. If you request a gist that is larger than 1MB,
the returned object says it's truncated. You can download the whole thing using
the `raw_url`, or for larger than 10 MB to the `git_pull_url`. (#26)

### BUG FIXES

* Fixed unicode problem on Windows (#37)
* Improved error catching (#28)
* `gist_create()` now works for an R script, didn't before (#29)

gistr 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.2 patched
* ubuntu 14.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

I have run R CMD check on the 9 downstream dependencies. Summary at <https://github.com/ropensci/gistr/tree/master/revdep>. None had problems.

---

This version fixes and improves some internals; no user facing changes.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/gistr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/gistr.git`
* Make sure to track progress upstream (i.e., on our version of `gistr` at `ropensci/gistr`) by doing `git remote add upstream https://github.com/ropensci/gistr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/gistr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and the below code block and proceed :) 

WHEN SHOWING EXAMPLES, DO NOT SHARE YOUR CREDENTIALS!
-->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
Stuff and things
------------------

## head of mtcars


```r
head(mtcars)
```

```
##                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
## Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
## Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
## Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
## Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
## Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
## Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
```

## head of iris


```r
head(iris)
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
## 5          5.0         3.6          1.4         0.2  setosa
## 6          5.4         3.9          1.7         0.4  setosa
```


## Scatter plot 


```r
plot(mpg ~ cyl, data=mtcars)
```

![plot of chunk unnamed-chunk-2](https://i.imgur.com/Ux99Wkj.png)

## Bar plot


```r
barplot(VADeaths)
```

![plot of chunk unnamed-chunk-3](https://i.imgur.com/IUwySie.png)

## Histogram


```r
hist(iris$Petal.Length)
```

![plot of chunk unnamed-chunk-4](https://i.imgur.com/6bPyEQT.png)
I shall generate random numbers.

    x <- rnorm(1000)

And I shall summarize them.

    summary(x)

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -2.83791 -0.65021  0.06139  0.01810  0.67508  3.02495
### Install and load alm


```r
install_github("ropensci/alm")
```


```r
library("alm")
```

### PLOS article data

The default in the `alm` package is for the PLOS ALM app. You do need to get an API key first here [http://alm.plos.org/](http://alm.plos.org/). You can pass in the `key` parameter or store in your `.Rprofile` file and pass in that way, or do `options(PlosApiKey = "yourkey")` and that will be stored for your current R session.


```r
alm_ids(doi='10.1371/journal.pone.0036240')
```

```
## $meta
##   total total_pages page error
## 1     1           1    1    NA
## 
## $data
##                       .id  pdf  html readers comments likes total
## 1               citeulike   NA    NA       5       NA    NA     5
## 2                crossref   NA    NA      NA       NA    NA     5
## 3                  nature   NA    NA      NA       NA    NA     1
## 4                  pubmed   NA    NA      NA       NA    NA     5
## 5                  scopus   NA    NA      NA       NA    NA     8
## 6                 counter 1152 16122      NA       NA    NA 17321
## 7        researchblogging   NA    NA      NA       NA    NA     1
## 8                     wos   NA    NA      NA       NA    NA     6
## 9                     pmc   64   194      NA       NA    NA   258
## 10               facebook   NA    NA      72       57    56   185
## 11               mendeley   NA    NA      70       NA    NA    70
## 12                twitter   NA    NA      NA      161    NA   161
## 13              wikipedia   NA    NA      NA       NA    NA     0
## 14          scienceseeker   NA    NA      NA       NA    NA     0
## 15         relativemetric   NA    NA      NA       NA    NA 43647
## 16                  f1000   NA    NA      NA       NA    NA     0
## 17               figshare    0     0      NA       NA     0     0
## 18              pmceurope   NA    NA      NA       NA    NA     5
## 19          pmceuropedata   NA    NA      NA       NA    NA     0
## 20            openedition   NA    NA      NA       NA    NA     0
## 21              wordpress   NA    NA      NA       NA    NA     0
## 22                 reddit   NA    NA      NA        0     0     0
## 23               datacite   NA    NA      NA       NA    NA     0
## 24             copernicus   NA    NA      NA       NA    NA     0
## 25        articlecoverage   NA    NA      NA        0    NA     0
## 26 articlecoveragecurated   NA    NA      NA        0    NA     0
## 27          plos_comments   NA    NA      NA        3    NA     4
```

### Crossref article data

You need to get a Crossref ALM API key first here [http://alm.labs.crossref.org/docs/Home](http://alm.labs.crossref.org/docs/Home), and pass in a different URL


```r
url <- "http://alm.labs.crossref.org/api/v5/articles"
alm_ids(doi='10.1371/journal.pone.0086859', url = url, key = getOption("crossrefalmkey"))
```

```
## $meta
##   total total_pages page error
## 1     1           1    1    NA
## 
## $data
##              .id pdf html readers comments likes total
## 1       crossref  NA   NA      NA       NA    NA     0
## 2       mendeley  NA   NA      NA       NA    NA     0
## 3       facebook  NA   NA      NA       NA    NA     0
## 4            pmc  NA   NA      NA       NA    NA     0
## 5      citeulike  NA   NA      NA       NA    NA     0
## 6         pubmed  NA   NA      NA       NA    NA     0
## 7      wordpress  NA   NA      NA       NA    NA     0
## 8         reddit  NA   NA      NA       NA    NA     0
## 9      wikipedia  NA   NA      NA       NA    NA     2
## 10      datacite  NA   NA      NA       NA    NA     0
## 11     pmceurope  NA   NA      NA       NA    NA     0
## 12 pmceuropedata  NA   NA      NA       NA    NA     0
```

### Public Knowledge Project (PKP) article data

You need to get a PKP ALM API key first here [http://pkp-alm.lib.sfu.ca/](http://pkp-alm.lib.sfu.ca/), and pass in a different URL


```r
url <- 'http://pkp-alm.lib.sfu.ca/api/v5/articles'
alm_ids(doi='10.3402/gha.v7.23554', url = url, key = getOption("pkpalmkey"))
```

```
## $meta
##   total total_pages page error
## 1     1           1    1    NA
## 
## $data
##                 .id pdf html readers comments likes total
## 1         citeulike  NA   NA       0       NA    NA     0
## 2            pubmed  NA   NA      NA       NA    NA     0
## 3         wikipedia  NA   NA      NA       NA    NA     0
## 4          mendeley  NA   NA       1       NA    NA     1
## 5          facebook  NA   NA       3        0     0     3
## 6            nature  NA   NA      NA       NA    NA     0
## 7  researchblogging  NA   NA      NA       NA    NA     0
## 8          crossref  NA   NA      NA       NA    NA     0
## 9     scienceseeker  NA   NA      NA       NA    NA     0
## 10        pmceurope  NA   NA      NA       NA    NA     0
## 11    pmceuropedata  NA   NA      NA       NA    NA     0
## 12      openedition  NA   NA      NA       NA    NA     0
## 13        wordpress  NA   NA      NA       NA    NA     0
## 14           reddit  NA   NA      NA       NA    NA     0
## 15       copernicus  NA   NA      NA       NA    NA     0
## 16           scopus  NA   NA      NA       NA    NA     0
## 17              pmc  NA   NA      NA       NA    NA     0
## 18   twitter_search  NA   NA      NA        0    NA     0
## 19            f1000  NA   NA      NA       NA    NA     0
```

<br>
__el fin!__


## Scatter plot 


```r
plot(mpg ~ cyl, data=mtcars)
```

![plot of chunk unnamed-chunk-2](https://i.imgur.com/UVoPCi3.png)

## Bar plot


```r
barplot(VADeaths)
```

![plot of chunk unnamed-chunk-3](https://i.imgur.com/13Re4qX.png)

## Histogram


```r
hist(iris$Petal.Length)
```

![plot of chunk unnamed-chunk-4](https://i.imgur.com/EItnHcZ.png)
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.5                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-07-28                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|gistr   |0.5.0 |0.9.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*