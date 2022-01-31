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

*Wow, no problems at all. :)**Wow, no problems at all. :)*gistr
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("gistr")
```

Or dev version from GitHub.

```{r eval=FALSE}
remotes::install_github("ropensci/gistr")
```

```{r}
library("gistr")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/gistr/issues).
* License: MIT
* Get citation information for `gistr` in R doing `citation(package = 'gistr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
```{r echo=FALSE}
library("knitr")
# upload all images to imgur.com
opts_knit$set(
  upload.fun = imgur_upload,
  base.url = NULL
)
```

```{r}
library("ggplot2")
```

## Scatter plot

```{r tidy=FALSE}
ggplot(mtcars, aes(cyl, hp)) + 
  geom_point() + 
  theme_grey(base_size = 18)
```

## Bar plot

```{r}
ggplot(iris, aes(Species, Sepal.Length)) + 
  stat_boxplot() +
  theme_grey(base_size = 18)
```
---
title: "Some plots and stuff"
author: "Foo Bar"
date: "June 30, 2015"
---

```{r echo=FALSE}
knitr::opts_knit$set(upload.fun = knitr::imgur_upload, base.url = NULL)
```

```{r}
library("ggplot2")
```

## Scatter plot

```{r tidy=FALSE}
ggplot(mtcars, aes(cyl, hp)) + 
  geom_point() + 
  theme_grey(base_size = 18)
```

## Bar plot

```{r}
ggplot(iris, aes(Species, Sepal.Length)) + 
  stat_boxplot() +
  theme_grey(base_size = 18)
```
## write table to file

```{r}
write.csv(mtcars, file = "mtcars.csv")
```
Stuff and things
------------------

## head of mtcars

```{r}
head(mtcars)
```

## head of iris

```{r}
head(iris)
```
```{r echo=FALSE}
library("knitr")
# upload all images to imgur.com
opts_knit$set(
  upload.fun = imgur_upload, 
  base.url = NULL
)
```

## Scatter plot 

```{r}
plot(mpg ~ cyl, data=mtcars)
```

## Bar plot

```{r}
barplot(VADeaths)
```

## Histogram

```{r}
hist(iris$Petal.Length)
```
```{r echo=FALSE}
knitr::opts_knit$set(upload.fun = knitr::imgur_upload, base.url = NULL)
```

## Scatter plot 

```{r}
plot(mpg ~ cyl, data=mtcars)
```

## Bar plot

```{r}
barplot(VADeaths)
```

## Histogram

```{r}
hist(iris$Petal.Length)
```
---
title: gistr introduction
author: Scott Chamberlain
date: "2020-07-28"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{gistr introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---





## Install

Stable version from CRAN



```r
install.packages("gistr")
```

Development version from GitHub



```r
remotes::install_github("ropensci/gistr")
```



```r
library("gistr")
```

## Authentication

There are two ways to authorise gistr to work with your GitHub account:

* Generate a personal access token (PAT) **with the gist scope selected** at [https://help.github.com/articles/creating-an-access-token-for-command-line-use](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token) and record it in the `GITHUB_PAT` environment variable.
  - To test out this approach, execute this in R: `Sys.setenv(GITHUB_PAT = "blahblahblah")`, where "blahblahblah" is the PAT you got from GitHub. Then take `gistr` out for a test drive.
  - If that works, you will probably want to define the GITHUB_PAT environment variable in a file such as `~/.bash_profile` or `~/.Renviron`.
* Interactively login into your GitHub account and authorise with OAuth.

Using the PAT is recommended.

Using the `gist_auth()` function you can authenticate separately first, or if you're not authenticated, this function will run internally with each function call. If you have a PAT, that will be used, if not, OAuth will be used.


```r
gist_auth()
```

## Workflow

In `gistr` you can use pipes, introduced perhaps first in R in the package `magrittr`, to pass outputs from one function to another. If you have used `dplyr` with pipes you can see the difference, and perhaps the utility, of this workflow over the traditional workflow in R. You can use a non-piping or a piping workflow with `gistr`. Examples below use a mix of both workflows. Here is an example of a piping workflow (with some explanation):


```r
file <- system.file("examples", "alm.md", package = "gistr")
gists(what = "minepublic")[[1]] %>% # List my public gists, and index to get just the 1st one
  add_files(file) %>% # Add a new file to that gist
  update() # update sends a PATCH command to the Gists API to add the file to your gist online
```

And a non-piping workflow that does the same exact thing:


```r
file <- system.file("examples", "alm.md", package = "gistr")
g <- gists(what = "minepublic")[[1]]
g <- add_files(g, file)
update(g)
```

Or you could string them all together in one line (but it's rather difficult to follow what's going on because you have to read from the inside out)


```r
file <- system.file("examples", "alm.md", package = "gistr")
update(add_files(gists(what = "minepublic")[[1]], file))
```

## Rate limit information


```r
rate_limit()
#> Rate limit: 5000
#> Remaining:  4998
#> Resets in:  59 minutes
```

## List gists

Limiting to a few results here to keep it brief


```r
gists(per_page = 2)
#> [[1]]
#> <gist>721a433293af4cb1fb0f66d7ccb37339
#>   URL: https://gist.github.com/721a433293af4cb1fb0f66d7ccb37339
#>   Description: Atividades 8o Ano
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:22Z / 2020-07-28T20:18:23Z
#>   Files: 8o_ano.md
#>   Truncated?: FALSE
#> 
#> [[2]]
#> <gist>d1064ad4f7fa9cd30b4409ade98fbd92
#>   URL: https://gist.github.com/d1064ad4f7fa9cd30b4409ade98fbd92
#>   Description: Averdonk Groessen
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:11Z / 2020-07-28T20:18:11Z
#>   Files: averdonk-groessen.markdown, index.pug, script.js, scripts, style.sass, styles
#>   Truncated?: FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
```

Since a certain date/time


```r
gists(since='2014-05-26T00:00:00Z', per_page = 2)
#> [[1]]
#> <gist>721a433293af4cb1fb0f66d7ccb37339
#>   URL: https://gist.github.com/721a433293af4cb1fb0f66d7ccb37339
#>   Description: Atividades 8o Ano
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:22Z / 2020-07-28T20:18:23Z
#>   Files: 8o_ano.md
#>   Truncated?: FALSE
#> 
#> [[2]]
#> <gist>d1064ad4f7fa9cd30b4409ade98fbd92
#>   URL: https://gist.github.com/d1064ad4f7fa9cd30b4409ade98fbd92
#>   Description: Averdonk Groessen
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:11Z / 2020-07-28T20:18:11Z
#>   Files: averdonk-groessen.markdown, index.pug, script.js, scripts, style.sass, styles
#>   Truncated?: FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
```

Request different types of gists, one of public, minepublic, mineall, or starred.


```r
gists('minepublic', per_page = 2)
#> [[1]]
#> <gist>792675323dea1961ce038b2e051d66d4
#>   URL: https://gist.github.com/792675323dea1961ce038b2e051d66d4
#>   Description: 
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:17:34Z / 2020-07-28T20:17:42Z
#>   Files: code.R
#>   Truncated?: FALSE
#> 
#> [[2]]
#> <gist>98badbaf7334a8a5707418543e178901
#>   URL: https://gist.github.com/98badbaf7334a8a5707418543e178901
#>   Description: a new cool gist
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:17:33Z / 2020-07-28T20:17:33Z
#>   Files: stuff.md
#>   Truncated?: FALSE
```


## List a single commit


```r
gist(id = 'f1403260eb92f5dfa7e1')
#> <gist>f1403260eb92f5dfa7e1
#>   URL: https://gist.github.com/f1403260eb92f5dfa7e1
#>   Description: Querying bitly from R 
#>   Public: TRUE
#>   Created/Edited: 2014-10-15T20:40:12Z / 2015-08-29T14:07:43Z
#>   Files: bitly_r.md
#>   Truncated?: FALSE
```

## Create gist

You can pass in files


```r
file <- system.file("examples", "stuff.md", package = "gistr")
gist_create(file, description='a new cool gist', browse = FALSE)
#> <gist>f7aa5e70e2fb40f4d92e972bcfa6224f
#>   URL: https://gist.github.com/f7aa5e70e2fb40f4d92e972bcfa6224f
#>   Description: a new cool gist
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:27Z / 2020-07-28T20:18:27Z
#>   Files: stuff.md
#>   Truncated?: FALSE
```

Or, wrap `gist_create()` around some code in your R session/IDE, with just the function name, and a `{'` at the start and a `}'` at the end.


```r
gist_create(code={'
x <- letters
numbers <- runif(8)
numbers

[1] 0.3229318 0.5933054 0.7778408 0.3898947 0.1309717 0.7501378 0.3206379 0.3379005
'})
```


```r
gist_create(code={'
x <- letters
numbers <- runif(8)
numbers

[1] 0.3229318 0.5933054 0.7778408 0.3898947 0.1309717 0.7501378 0.3206379 0.3379005
'}, browse=FALSE)
#> <gist>394103ff248f2be67bbccce863a89ca7
#>   URL: https://gist.github.com/394103ff248f2be67bbccce863a89ca7
#>   Description: 
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:29Z / 2020-07-28T20:18:29Z
#>   Files: code.R
#>   Truncated?: FALSE
```

### knit and create

You can also knit an input file before posting as a gist:


```r
file <- system.file("examples", "stuff.Rmd", package = "gistr")
gist_create(file, description='a new cool gist', knit=TRUE)
#> <gist>4162b9c53479fbc298db
#>   URL: https://gist.github.com/4162b9c53479fbc298db
#>   Description: a new cool gist
#>   Public: TRUE
#>   Created/Edited: 2014-10-27T16:07:31Z / 2014-10-27T16:07:31Z
#>   Files: stuff.md
#>   Truncated?: FALSE
```

Or code blocks before (note that code blocks without knitr block demarcations will result in unexecuted code):


```r
gist_create(code={'
x <- letters
(numbers <- runif(8))
'}, knit=TRUE)
#> <gist>ec45c396dee4aa492139
#>   URL: https://gist.github.com/ec45c396dee4aa492139
#>   Description:
#>   Public: TRUE
#>   Created/Edited: 2014-10-27T16:09:09Z / 2014-10-27T16:09:09Z
#>   Files: file81720d1ceff.md
#>   Truncated?: FALSE
```

## knit code from file path, code block, or gist file

knit a local file


```r
file <- system.file("examples", "stuff.Rmd", package = "gistr")
run(file, knitopts = list(quiet=TRUE)) %>% gist_create(browse = FALSE)
#> <gist>ceeddf11ebf775a1fe465f6791e15323
#>   URL: https://gist.github.com/ceeddf11ebf775a1fe465f6791e15323
#>   Description: 
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:31Z / 2020-07-28T20:18:31Z
#>   Files: stuff.md
#>   Truncated?: FALSE
```



knit a code block (knitr code block notation missing, do add that in) (result not shown)


```r
run({'
x <- letters
(numbers <- runif(8))
'}) %>% gist_create
```

knit a file from a gist, has to get file first (result not shown)


```r
gists('minepublic')[[1]] %>% run() %>% update()
```

## working with images

The GitHub API doesn't let you upload binary files (e.g., images) via their HTTP API, which we use in `gistr`. There is a workaround.

If you are using `.Rmd` or `.Rnw` files, you can set `imgur_inject = TRUE` in `gistr_create()` so that imgur knit options are injected at the top of your file so that images will be uploaded to imgur. Alternatively, you can do this yourself, setting knit options to use imgur.

A file already using imgur


```r
file <- system.file("examples", "plots_imgur.Rmd", package = "gistr")
gist_create(file, knit=TRUE)
#> <gist>1a6e7f7d6ddb739fce0b
#>   URL: https://gist.github.com/1a6e7f7d6ddb739fce0b
#>   Description:
#>   Public: TRUE
#>   Created/Edited: 2015-03-19T00:20:48Z / 2015-03-19T00:20:48Z
#>   Files: plots_imgur.md
```

A file _NOT_ already using imgur


```r
file <- system.file("examples", "plots.Rmd", package = "gistr")
gist_create(file, knit=TRUE, imgur_inject = TRUE)
#> <gist>ec9987ad245bbc668c72
#>   URL: https://gist.github.com/ec9987ad245bbc668c72
#>   Description:
#>   Public: TRUE
#>   Created/Edited: 2015-03-19T00:21:13Z / 2015-03-19T00:21:13Z
#>   Files: plots.md
#>   Truncated?: FALSE
```

## List commits on a gist


```r
gists()[[1]] %>% commits()
#> [[1]]
#> <commit>
#>   Version: b640145362bb2272c5f0a06ae2aee3241b8f6f59
#>   User: sckott
#>   Commited: 2020-07-28T20:18:28Z
#>   Commits [total, additions, deletions]: [5,5,0]
```

## Star a gist

Star


```r
gist('cbb0507082bb18ff7e4b') %>% star()
#> <gist>cbb0507082bb18ff7e4b
#>   URL: https://gist.github.com/cbb0507082bb18ff7e4b
#>   Description: This is my technical interview cheat sheet.  Feel free to fork it or do whatever you want with it.  PLEASE let me know if there are any errors or if anything crucial is missing.  I will add more links soon.
#>   Public: TRUE
#>   Created/Edited: 2014-05-02T19:43:13Z / 2018-04-16T21:11:53Z
#>   Files: The Technical Interview Cheat Sheet.md
#>   Truncated?: FALSE
```

Unstar


```r
gist('cbb0507082bb18ff7e4b') %>% unstar()
#> <gist>cbb0507082bb18ff7e4b
#>   URL: https://gist.github.com/cbb0507082bb18ff7e4b
#>   Description: This is my technical interview cheat sheet.  Feel free to fork it or do whatever you want with it.  PLEASE let me know if there are any errors or if anything crucial is missing.  I will add more links soon.
#>   Public: TRUE
#>   Created/Edited: 2014-05-02T19:43:13Z / 2018-04-16T21:27:36Z
#>   Files: The Technical Interview Cheat Sheet.md
#>   Truncated?: FALSE
```

## Edit a gist

Add files


```r
file <- system.file("examples", "alm.md", package = "gistr")
gists(what = "minepublic")[[1]] %>%
  add_files(file) %>%
  update()
#> <gist>394103ff248f2be67bbccce863a89ca7
#>   URL: https://gist.github.com/394103ff248f2be67bbccce863a89ca7
#>   Description: 
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:29Z / 2020-07-28T20:18:34Z
#>   Files: alm.md, code.R
#>   Truncated?: FALSE, FALSE
```

Delete files


```r
file <- system.file("examples", "alm.md", package = "gistr")
gists(what = "minepublic")[[1]] %>%
  delete_files(file) %>%
  update()
#> <gist>394103ff248f2be67bbccce863a89ca7
#>   URL: https://gist.github.com/394103ff248f2be67bbccce863a89ca7
#>   Description: 
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:29Z / 2020-07-28T20:18:36Z
#>   Files: code.R
#>   Truncated?: FALSE
```

## Open a gist in your default browser


```r
gists()[[1]] %>% browse()
```

> Opens the gist in your default browser

## Get embed script


```r
gists()[[1]] %>% embed()
#> [1] "<script src=\"https://gist.github.com/anupamgogoi-wso2/9230e4ade0074320ec95992ba2fbf518.js\"></script>"
```

## List forks

Returns a list of `gist` objects, just like `gists()`


```r
gist(id='1642874') %>% forks(per_page=2)
#> [[1]]
#> <gist>1642989
#>   URL: https://gist.github.com/1642989
#>   Description: Spline Transition
#>   Public: TRUE
#>   Created/Edited: 2012-01-19T21:45:20Z / 2019-10-23T20:09:07Z
#>   Files: 
#>   Truncated?: 
#> 
#> [[2]]
#> <gist>1643051
#>   URL: https://gist.github.com/1643051
#>   Description: Line Transition (Broken)
#>   Public: TRUE
#>   Created/Edited: 2012-01-19T21:51:30Z / 2019-10-23T20:08:44Z
#>   Files: 
#>   Truncated?:
```

## Fork a gist

Returns a `gist` object


```r
g <- gists()
(forked <- g[[ sample(seq_along(g), 1) ]] %>% fork())
#> <gist>72a9327aac0cab09900b04a027069ae5
#>   URL: https://gist.github.com/72a9327aac0cab09900b04a027069ae5
#>   Description: Rimworld output log published using HugsLib
#>   Public: TRUE
#>   Created/Edited: 2020-07-28T20:18:43Z / 2020-07-28T20:18:43Z
#>   Files: output_log.txt
#>   Truncated?: FALSE
```





## Example use cases

_Round-trip storage of a data frame_

Maybe you want to use a gist to share some data as an alternative to `dput`? We can do this by writing our `data.frame` to a temporary buffer and passing it to `gist_create`. We can read the data back from the gist by passing its content to `read.csv`.



```r
data(iris)

str <- ''
tc  <- textConnection('str', 'w', local = TRUE)
write.csv(iris, file = tc, row.names = FALSE)
close(tc)

content <- list(content=paste(as.character(str), collapse='\n'))

gistr::gist_create(code = {
  content$content
}, description = "using a gist as a data store", 
filename = "iris.csv")
#> <gist>c7dfe593f4944df4818df884689734f9
#>   URL: https://gist.github.com/c7dfe593f4944df4818df884689734f9
#>   Description: using a gist as a data store
#>   Public: TRUE
#>   Created/Edited: 2019-07-18T14:23:23Z / 2019-07-18T14:23:23Z
#>   Files: iris.csv
#>   Truncated?: FALSE

output <- read.csv(
  text = gist(gists(what = "minepublic", per_page = 1)[[1]]$id)$
    files$iris.csv$content)

identical(output, iris)
#> TRUE
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/star.R
\name{star}
\alias{star}
\alias{unstar}
\alias{star_check}
\title{Star a gist}
\usage{
star(gist, ...)

unstar(gist, ...)

star_check(gist, ...)
}
\arguments{
\item{gist}{A gist object or something that can be coerced to a gist object.}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
A message, and a gist object, the same one input to the function.
}
\description{
Star a gist
}
\examples{
\dontrun{
id <- '4ac33b9c00751fddc7f8'
gist(id) \%>\% star()
gist(id) \%>\% star_check()
gist(id) \%>\% unstar()
gist(id) \%>\% unstar() \%>\% star()
gist(id) \%>\% star_check()
gist(id) \%>\%
  star() \%>\%
  star_check()
  
# pass in a url
x <- "https://gist.github.com/expersso/4ac33b9c00751fddc7f8"
gist(x) \%>\% star
gist(x) \%>\% unstar
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse}
\alias{browse}
\title{Open a gist on GitHub}
\usage{
browse(gist, what = "html")
}
\arguments{
\item{gist}{A gist object or something that can be coerced to a gist object.}

\item{what}{One of html (default), json, forks, commits, or comments.}
}
\description{
Open a gist on GitHub
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run}
\alias{run}
\title{Run a .Rmd file}
\usage{
run(x, filename = "code.R", knitopts = list())
}
\arguments{
\item{x}{Input, one of: code wrapped in curly brackets and quotes, a file
path to an .Rmd file, or a gist.}

\item{filename}{Name of the file to create, only used if \code{code}
parameter is used. Default to \code{code.R}}

\item{knitopts}{(list) List of variables passed on to
\code{\link[knitr:knit]{knitr::knit()}}}
}
\value{
A path, unless a gist object is passed in, in which case a gist
object is returned.
}
\description{
Run a .Rmd file
}
\examples{
\dontrun{
# run a local file
file <- system.file("examples", "stuff.Rmd", package = "gistr")
run(file) \%>\% gist_create

# run code
run({'
```{r}
x <- letters
(numbers <- runif(8))
```
'}) \%>\% gist_create

# run a file from a gist, has to get file first
gists('minepublic')[[2]] \%>\% run() \%>\% update()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_create_obj.R
\name{gist_create_obj}
\alias{gist_create_obj}
\title{Create a gist from an R object}
\usage{
gist_create_obj(
  x = NULL,
  description = "",
  public = TRUE,
  browse = TRUE,
  pretty = TRUE,
  filename = "file.txt",
  ...
)
}
\arguments{
\item{x}{An R object, any of data.frame, matrix, list, character, numeric}

\item{description}{(character) Brief description of gist (optional)}

\item{public}{(logical) Whether gist is public (default: \code{TRUE})}

\item{browse}{(logical) To open newly create gist in default browser
(default: \code{TRUE})}

\item{pretty}{(logical) For data.frame and matrix objects, create
a markdown table. If FALSE, pushes up json. (default: \code{TRUE})}

\item{filename}{Name of the file to create. Default: \code{file.txt}}

\item{...}{Further args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
Create a gist from an R object
}
\details{
This function is specifically for going from R objects to a gist,
whereas \code{\link[=gist_create]{gist_create()}} is for going from files or executing code
}
\examples{
\dontrun{
## data.frame
### by default makes pretty table in markdown format
row.names(mtcars) <- NULL
gist_create_obj(mtcars)
gist_create_obj(iris)
### or just push up json
gist_create_obj(mtcars, pretty = FALSE)

## matrix
gist_create_obj(as.matrix(mtcars))
## list
gist_create_obj(apply(mtcars, 1, as.list))
## character
gist_create_obj("hello, world")
## numeric
gist_create_obj(runif(10))

## Assign a specific file name
gist_create_obj("
## header2

hey there!", filename = "my_markdown.md")
}
}
\seealso{
\code{\link[=gist_create]{gist_create()}}, \code{\link[=gist_create_git]{gist_create_git()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gistr-package.R
\name{create_gists}
\alias{create_gists}
\title{Create gists}
\description{
Creating gists in \code{gistr} can be done with any of
three functions:
\itemize{
\item \code{\link[=gist_create]{gist_create()}} - Create gists from files or code blocks,
using  the GitHub HTTP API. Because this function uses the GitHub HTTP API,
it does not work for binary files. However, you can get around this for
images by using knitr's hook to upload images to eg., imgur. In addition,
it's difficult to include artifacts from the knit-ing process.
\item \code{\link[=gist_create_git]{gist_create_git()}} - Create gists from files or code
blocks, using  git. Because this function uses git, you have more
flexibility than with the above function: you can include any binary files,
and can easily upload all artifacts.
\item \code{\link[=gist_create_obj]{gist_create_obj()}} - Create gists from R objects: data.frame, list,
character string, matrix, or numeric. Uses the GitHub HTTP API.
}

It may seem a bit odd to have three separate functions for creating gists.
\code{\link[=gist_create]{gist_create()}} was created first, and was out for a bit, so when
we had the idea to create gists via git (\code{\link[=gist_create_git]{gist_create_git()}}) and
from R objects (\code{\link[=gist_create_obj]{gist_create_obj()}}), it made sense to have a
different API for creating gists via the HTTP API, git, and from R objects.
We could have thrown everything into \code{\link[=gist_create]{gist_create()}}, but it would
have been a massive function, with far too many parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabl.R
\name{tabl}
\alias{tabl}
\alias{tabl_data}
\title{Make a table from gist or commit class or a list of either}
\usage{
tabl(x, ...)

tabl_data(x)
}
\arguments{
\item{x}{Either a gist or commit class object or a list of either}

\item{...}{Ignored}
}
\value{
A data.frame or list of data.frame's
}
\description{
Make a table from gist or commit class or a list of either
}
\details{
For commits we return a single data.frame. For gists, we always
return a list so that we are returning data consistently,
regardless of variable return data. So you can always index to the main
data.frame with gist metadata and file info by doing \code{result$data},
and likewise for  forks \code{result$forks} and history
\code{result$history}
}
\examples{
\dontrun{
# from a gist object
x <- as.gist('f1403260eb92f5dfa7e1')
res <- tabl(x)
res$data
res$forks
res$history

# from a list
ss <- gists('minepublic')
tabl(ss[1:3])
lapply(tabl(ss[1:3]), "[[", "data")
# index to data slots, but also make single data.frame
tabl_data(tabl(ss[1:3]))
## manipulate with dplyr
library("dplyr")
tabl_data(tabl(ss[1:30])) \%>\% 
  select(id, description, owner_login) \%>\% 
  filter(grepl("gist gist gist", description))

# commits
x <- gists()[[2]] \%>\% commits()
tabl(x[[1]])

## many
x <- sapply(gists(per_page = 100), commits)
tabl(x) \%>\%
  select(id, login, change_status.total, url) \%>\% 
  filter(change_status.total > 50)
  
# pass in a url
gist("https://gist.github.com/expersso/4ac33b9c00751fddc7f8") \%>\% tabl
## many
gg <- gists()
(urls <- vapply(gg, "[[", "", "html_url"))
lapply(urls[1:5], as.gist) \%>\% tabl()

# gist with forks and history
gist('1642874') \%>\% tabl

# gist with history, no forks
gist('c96d2e453c95d0166408') \%>\% tabl 
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist.R
\name{gist}
\alias{gist}
\alias{as.gist}
\title{Get a gist}
\usage{
gist(id, revision = NULL, ...)

as.gist(x)
}
\arguments{
\item{id}{(character) A gist id, or a gist URL}

\item{revision}{(character) A sha. optional}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}

\item{x}{Object to coerce. Can be an integer (gist id), string
(gist id), a gist, or an list that can be coerced to a gist.}
}
\description{
Get a gist
}
\details{
If a file is larger than ~1 MB, the content of the file given back
is truncated, so you won't get the entire contents. In the return S3 object
that's printed, we tell you at the bottom whether each file is truncated or
not. If a file is, simply get the \code{raw_url} URL for the file (see
example below), then retrieve from that. If the file is very big, you may
need to clone the file using git, etc.
}
\examples{
\dontrun{
gist('f1403260eb92f5dfa7e1')

as.gist('f1403260eb92f5dfa7e1')
as.gist(10)
as.gist(gist('f1403260eb92f5dfa7e1'))

# get a specific revision of a gist
id <- 'c1e2cb547d9f22bd314da50fe9c7b503'
gist(id, 'a5bc5c143beb697f23b2c320ff5a8dacf960b0f3')
gist(id, 'b70d94a8222a4326dff46fc85bc69d0179bd1da2')
gist(id, '648bb44ab9ae59d57b4ea5de7d85e24103717e8b')
gist(id, '0259b13c7653dc95e20193133bcf71811888cbe6')

# from a url, or partial url
x <- "https://gist.github.com/expersso/4ac33b9c00751fddc7f8"
x <- "gist.github.com/expersso/4ac33b9c00751fddc7f8"
x <- "gist.github.com/4ac33b9c00751fddc7f8"
x <- "expersso/4ac33b9c00751fddc7f8"
as.gist(x)

ids <- sapply(gists(), "[[", "id")
gist(ids[1])
gist(ids[2])
gist(ids[3])
gist(ids[4])

gist(ids[1]) \%>\% browse()

## If a gist file is > a certain size it is truncated
## in this case, we let you know in the return object that it is truncated
## e.g.
(bigfile <- gist(id = "b74b878fd7d9176a4c52"))
## then get the raw_url, and retrieve the file
url <- bigfile$files$`plossmall.json`$raw_url
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forks.R
\name{forks}
\alias{forks}
\title{List forks on a gist}
\usage{
forks(gist, page = NULL, per_page = 30, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{page}{(integer) Page number to return.}

\item{per_page}{(integer) Number of items to return per page. Default 30.
Max 100.}

\item{...}{Further named args to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
A list of gist class objects
}
\description{
List forks on a gist
}
\examples{
\dontrun{
gist(id='1642874') \%>\% forks(per_page=2)
gist(id = "8172796") \%>\% forks()

# pass in a url
gist("https://gist.github.com/expersso/4ac33b9c00751fddc7f8") \%>\% forks 
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/delete.R
\name{delete}
\alias{delete}
\title{Delete a gist}
\usage{
delete(gist, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Delete a gist
}
\examples{
\dontrun{
gists("minepublic")[[29]] \%>\% delete()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_limit.R
\name{rate_limit}
\alias{rate_limit}
\title{Get rate limit information}
\usage{
rate_limit(...)
}
\arguments{
\item{...}{Named args to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Get rate limit information
}
\examples{
\dontrun{
rate_limit()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forks.R
\name{fork}
\alias{fork}
\title{Fork a gist}
\usage{
fork(gist, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{...}{Further named args to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
A gist class object
}
\description{
Fork a gist
}
\examples{
\dontrun{
# fork a gist
w <- gists()[[1]] \%>\% fork()

# browse to newly forked gist
browse(w)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gists.R
\name{gists}
\alias{gists}
\title{List gists}
\usage{
gists(what = "public", since = NULL, page = NULL, per_page = 30, ...)
}
\arguments{
\item{what}{(character) What gists to return. One of public, minepublic,
mineall, or starred. If an id is given for a gist, this parameter is ignored.}

\item{since}{(character) A timestamp in ISO 8601 format:
YYYY-MM-DDTHH:MM:SSZ. Only gists updated at or after this time are returned.}

\item{page}{(integer) Page number to return.}

\item{per_page}{(integer) Number of items to return per page. Default 30.
Max 100.}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
List public gists, your own public gists, all your gists, by gist id, or
query by date.
}
\details{
When \code{what = "mineall"}, we use
\code{getOption("github.username")} internally to get your GitHub user name.
Make sure to set your GitHub user name
as an R option like \code{options(github.username = "foobar")} in your
\code{.Rprofile} file. If we can't find you're user name, we'll stop with an
error.
}
\examples{
\dontrun{
# Public gists
gists()
gists(per_page=2)
gists(page=3)
# Public gists created since X time
gists(since='2014-05-26T00:00:00Z')
# Your public gists
gists('minepublic')
gists('minepublic', per_page=2)
# Your private and public gists
gists('mineall')
# Your starred gists
gists('starred')
# pass in curl options
gists(per_page=1, verbose=TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commits.R
\name{commits}
\alias{commits}
\title{List gist commits}
\usage{
commits(gist, page = NULL, per_page = 30, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{page}{(integer) Page number to return.}

\item{per_page}{(integer) Number of items to return per page.
Default 30. Max 100.}

\item{...}{Further named args to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
List gist commits
}
\examples{
\dontrun{
gists()[[1]] \%>\% commits()
gist(id = '1f399774e9ecc9153a6f') \%>\% commits(per_page = 5)

# pass in a url
gist("https://gist.github.com/expersso/4ac33b9c00751fddc7f8") \%>\% commits
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_auth.R
\name{gist_auth}
\alias{gist_auth}
\title{Authorize with GitHub.}
\usage{
gist_auth(app = gistr_app, reauth = FALSE)
}
\arguments{
\item{app}{An \code{\link[httr:oauth_app]{httr::oauth_app()}} for GitHub. The default uses an
application \code{gistr_oauth} created by Scott Chamberlain.}

\item{reauth}{(logical) Force re-authorization?}
}
\value{
a named list, with a single slot for \code{Authorization}, with a single
element with the token - this is the expected format needed when passed
down to the http request
}
\description{
This function is run automatically to allow gistr to access your GitHub
account.
}
\details{
There are two ways to authorise gistr to work with your GitHub account:
\itemize{
\item Generate a personal access token with the gist scope selected, and set it
as the \code{GITHUB_PAT} environment variable per session using \code{Sys.setenv}
or across sessions by adding it to your \code{.Renviron} file or similar.
See
https://help.github.com/articles/creating-an-access-token-for-command-line-use
for help
\item Interactively login into your GitHub account and authorise with OAuth.
}

Using \code{GITHUB_PAT} is recommended.
}
\examples{
\dontrun{
gist_auth()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_files.R
\name{add_files}
\alias{add_files}
\alias{update_files}
\alias{delete_files}
\alias{rename_files}
\title{Add files to a gist object}
\usage{
add_files(gist, ...)

update_files(gist, ...)

delete_files(gist, ...)

rename_files(gist, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Add files to a gist object
}
\examples{
\dontrun{
add_files("~/stuff.Rmd")
# update_files()
# delete_files()
# rename_files()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_create_git.R
\name{gist_create_git}
\alias{gist_create_git}
\title{Create a gist via git instead of the GitHub Gists HTTP API}
\usage{
gist_create_git(
  files = NULL,
  description = "",
  public = TRUE,
  browse = TRUE,
  knit = FALSE,
  code = NULL,
  filename = "code.R",
  knitopts = list(),
  renderopts = list(),
  include_source = FALSE,
  artifacts = FALSE,
  imgur_inject = FALSE,
  git_method = "ssh",
  sleep = 1,
  ...
)
}
\arguments{
\item{files}{Files to upload. this or \code{code} param must be passed}

\item{description}{(character) Brief description of gist (optional)}

\item{public}{(logical) Whether gist is public (default: TRUE)}

\item{browse}{(logical) To open newly create gist in default browser
(default: TRUE)}

\item{knit}{(logical) Knit code before posting as a gist? If the file
has a \code{.Rmd}  or \code{.Rnw} extension, we run the file with
\code{\link[knitr]{knit}}, and if it has a \code{.R} extension, then
we use \code{\link[rmarkdown]{render}}}

\item{code}{Pass in any set of code. This can be a single R object, or
many lines of code wrapped in quotes, then curly brackets (see examples
below). this or \code{files} param must be passed}

\item{filename}{Name of the file to create, only used if \code{code}
parameter is used. Default to \code{code.R}}

\item{knitopts, renderopts}{(list) List of variables passed on to
\code{\link[knitr]{knit}}, or \code{\link[rmarkdown]{render}}}

\item{include_source}{(logical) Only applies if \code{knit=TRUE}. Include
source file in the gist in addition to the knitted output.}

\item{artifacts}{(logical/character) Include artifacts or not.
If \code{TRUE}, includes all artifacts. Or you can pass in a file extension
to only upload artifacts of certain file exensions. Default: \code{FALSE}}

\item{imgur_inject}{(logical) Inject \code{\link[knitr]{imgur_upload}}
into your \code{.Rmd} file to upload files to \url{https://imgur.com/}.
This will be ignored  if the file is a sweave/latex file because the
rendered pdf can't be uploaded anyway. Default: FALSE}

\item{git_method}{(character) One of ssh (default) or https. If a remote
already exists, we use that remote, and this parameter is ignored.}

\item{sleep}{(integer) Seconds to sleep after creating gist, but before
collecting metadata on the gist. If uploading a lot of stuff, you may want to
set this to a higher value, otherwise, you may not get accurate metadata for
your gist. You can of course always refresh afterwards by calling \code{gist}
with your gist id.}

\item{...}{Further args passed on to \code{\link[crul]{verb-POST}}}
}
\description{
Create a gist via git instead of the GitHub Gists HTTP API
}
\details{
Note that when \code{browse=TRUE} there is a slight delay in when
we open up the gist in your default browser and when the data will display
in the gist. We could have this function sleep a while and guess when it
will be ready, but instead we open your gist right after we're done sending
the data to GitHub. Make sure to refresh the page if you don't see your
content right away.

Likewise, the object that is returned from this function call may not have
the updated and correct file information. You can retrieve that easily by
calling \code{\link[=gist]{gist()}} with the gist id.

This function uses git instead of the HTTP API, and thus requires
the R package \code{git2r}. If you don't have \code{git2r} installed, and
try to use this function, it will stop and tell you to install \code{git2r}.

This function using git is better suited than \code{\link[=gist_create]{gist_create()}}
for use cases involving:
\itemize{
\item Big files - The GitHub API allows only files of up to 1 MB in size.
Using git we can get around that limit.
\item Binary files - Often artifacts created are binary files like
\code{.png}. The GitHub API doesn't allow transport of binary files, but
we can do that with git.
}

Another difference between this function and \code{\link[=gist_create]{gist_create()}} is
that this function can collect all artifacts coming out of a knit process.

If a gist is somehow deleted, or the remote changes, when you try to push
to the same gist again, everything should be fine. We now use
\code{tryCatch} on the  push attempt, and if it fails, we'll add a new
remote (which means a new gist), and push again.
}
\examples{
\dontrun{
# prepare a directory and a file
unlink("~/gitgist", recursive = TRUE)
dir.create("~/gitgist")
file <- system.file("examples", "stuff.md", package = "gistr")
writeLines(readLines(file), con = "~/gitgist/stuff.md")

# create a gist
gist_create_git(files = "~/gitgist/stuff.md")

## more than one file can be passed in
unlink("~/gitgist2", recursive = TRUE)
dir.create("~/gitgist2")
file.copy(file, "~/gitgist2/")
cat("hello world", file = "~/gitgist2/hello_world.md")
list.files("~/gitgist2")
gist_create_git(c("~/gitgist2/stuff.md", "~/gitgist2/hello_world.md"))

# Include all files in a directory
unlink("~/gitgist3", recursive = TRUE)
dir.create("~/gitgist3")
cat("foo bar", file="~/gitgist3/foobar.txt")
cat("hello", file="~/gitgist3/hello.txt")
list.files("~/gitgist3")
gist_create_git("~/gitgist3")

# binary files
png <- system.file("examples", "file.png", package = "gistr")
unlink("~/gitgist4", recursive = TRUE)
dir.create("~/gitgist4")
file.copy(png, "~/gitgist4/")
list.files("~/gitgist4")
gist_create_git(files = "~/gitgist4/file.png")

# knit files first, then push up
# note: by default we don't upload images, but you can do that, 
# see next example
rmd <- system.file("examples", "plots.Rmd", package = "gistr")
unlink("~/gitgist5", recursive = TRUE)
dir.create("~/gitgist5")
file.copy(rmd, "~/gitgist5/")
list.files("~/gitgist5")
gist_create_git("~/gitgist5/plots.Rmd", knit = TRUE)

# collect all/any artifacts from knitting process
arts <- system.file("examples", "artifacts_eg1.Rmd", package = "gistr")
unlink("~/gitgist6", recursive = TRUE)
dir.create("~/gitgist6")
file.copy(arts, "~/gitgist6/")
list.files("~/gitgist6")
gist_create_git("~/gitgist6/artifacts_eg1.Rmd", knit = TRUE, 
   artifacts = TRUE)

# from a code block
gist_create_git(code={'
x <- letters
numbers <- runif(8)
numbers

[1] 0.3229318 0.5933054 0.7778408 0.3898947 0.1309717 0.7501378 0.3206379 0.3379005
'}, filename="my_cool_code.R")

# Use https instead of ssh
png <- system.file("examples", "file.png", package = "gistr")
unlink("~/gitgist7", recursive = TRUE)
dir.create("~/gitgist7")
file.copy(png, "~/gitgist7/")
list.files("~/gitgist7")
gist_create_git(files = "~/gitgist7/file.png", git_method = "https")
}
}
\seealso{
\code{\link[=gist_create]{gist_create()}}, \code{\link[=gist_create_obj]{gist_create_obj()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gistr-package.R
\docType{package}
\name{gistr-package}
\alias{gistr-package}
\alias{gistr}
\title{R client for GitHub gists}
\description{
R client for GitHub gists.
}
\details{
gistr allows you to peform actions on gists, including listing, forking,
starring, creating, deleting, updating, etc.

There are two ways to authorise gistr to work with your GitHub account:
\itemize{
\item Generate a personal access token (PAT) at
\url{https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token}
and record it in the \code{GITHUB_PAT} envar.
\item Interactively login into your GitHub account and authorise with OAuth.
}

Using the \code{GITHUB_PAT} is recommended.
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Ramnath Vaidyanathan \email{ramnath.vaidya@gmail.com}

Karthik Ram \email{karthik.ram@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_save.R
\name{gist_save}
\alias{gist_save}
\alias{gist_open}
\title{Save gist files to disk}
\usage{
gist_save(gist, path = ".")

gist_open(x)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{path}{Root path to write to, a directory, not a file b/c a gist can
contain  many files. A folder is created with name of the gist id within
this root directory.  File names will be the same as given in the gist.}

\item{x}{An object of class \code{gist_files} (the output from
\code{\link[=gist_save]{gist_save()}}}
}
\value{
An object of class \code{gist_files}, S3 object containing file
paths
}
\description{
Save gist files to disk
}
\details{
\code{gist_save}: files are written into a new folder, named by the gist id,
e.g., \code{a65ac7e56b7b3f746913}

\code{gist_open}: opens files in your editor/R GUI. Internally, uses
\code{\link[=file.edit]{file.edit()}} to open files, using \code{getOption("editor")} to
open the files. If you're in R.app or RStudio, or other IDE's, files will
open in the IDE (I think).
}
\examples{
\dontrun{
gist("a65ac7e56b7b3f746913") \%>\% gist_save()
gist("a65ac7e56b7b3f746913") \%>\% gist_save() \%>\% gist_open()
gist("https://gist.github.com/expersso/4ac33b9c00751fddc7f8") \%>\%
  gist_save()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_create.R
\name{gist_create}
\alias{gist_create}
\title{Create a gist}
\usage{
gist_create(
  files = NULL,
  description = "",
  public = TRUE,
  browse = TRUE,
  code = NULL,
  filename = "code.R",
  knit = FALSE,
  knitopts = list(),
  renderopts = list(),
  include_source = FALSE,
  imgur_inject = FALSE,
  rmarkdown = FALSE,
  ...
)
}
\arguments{
\item{files}{Files to upload. this or \code{code} param must be passed}

\item{description}{(character) Brief description of gist (optional)}

\item{public}{(logical) Whether gist is public (default: TRUE)}

\item{browse}{(logical) To open newly create gist in default browser
(default: TRUE)}

\item{code}{Pass in any set of code. This can be a single R object, or
many lines of code wrapped in quotes, then curly brackets (see examples
below). this or \code{files} param must be passed}

\item{filename}{Name of the file to create, only used if \code{code}
parameter is used. Default to \code{code.R}}

\item{knit}{(logical) Knit code before posting as a gist? If the file
has a \code{.Rmd}  or \code{.Rnw} extension, we run the file with
\code{\link[knitr]{knit}}, and if it has a \code{.R} extension, then
we use \code{\link[rmarkdown]{render}}}

\item{knitopts, renderopts}{(list) List of variables passed on to
\code{\link[knitr]{knit}}, or \code{\link[rmarkdown]{render}}}

\item{include_source}{(logical) Only applies if \code{knit=TRUE}. Include
source file in the gist in addition to the knitted output.}

\item{imgur_inject}{(logical) Inject \code{\link[knitr]{imgur_upload}}
into your \code{.Rmd} file to upload files to \url{https://imgur.com/}.
This will be ignored  if the file is a sweave/latex file because the
rendered pdf can't be uploaded anyway. Default: FALSE}

\item{rmarkdown}{(logical) If \code{TRUE}, use \code{\link[rmarkdown:render]{rmarkdown::render()}} instead of
\code{\link[knitr:knit]{knitr::knit()}} to render the document.}

\item{...}{Further args passed on to \code{\link[crul]{verb-POST}}}
}
\description{
Create a gist
}
\examples{
\dontrun{
file <- tempfile()
cat("hello world", file = file)
gist_create(files=file, description='a new cool gist')

file1 <- tempfile()
file2 <- tempfile()
cat("foo bar", file = file1)
cat("foo bar", file = file2)
gist_create(files=c(file1, file2), description='spocc demo files')

# include any code by passing to the code parameter
gist_create(code={'
x <- letters
numbers <- runif(10)
numbers
'})

# Knit an .Rmd file before posting as a gist
file <- system.file("examples", "stuff.Rmd", package = "gistr")
gist_create(file, description='a new cool gist', knit=TRUE)

file <- system.file("examples", "plots.Rmd", package = "gistr")
gist_create(file, description='some plots', knit=TRUE)

# an .Rnw file
file <- system.file("examples", "rnw_example.Rnw", package = "gistr")
gist_create(file)
gist_create(file, knit=TRUE)

# Knit code input before posting as a gist
gist_create(code={'
```{r}
x <- letters
(numbers <- runif(8))
```
'}, knit=TRUE)

url <- "https://raw.githubusercontent.com/ropensci/geojsonio/master/inst/examples/zillow_or.geojson"
json <- crul::HttpClient$new(url)$get()$parse("UTF-8")
gist_create(code = json, filename = "zillow_or.geojson")

# Knit and include source file, so both files are in the gist
file <- system.file("examples", "stuff.Rmd", package = "gistr")
gist_create(file, knit=TRUE, include_source=TRUE)

gist_create(code={'
```{r}
x <- letters
(numbers <- runif(8))
```
'}, filename="code.Rmd", knit=TRUE, include_source=TRUE)

# Uploading images created during knit process
## using imgur - if you're file uses imgur or similar, you're good
file <- system.file("examples", "plots_imgur.Rmd", package = "gistr")
cat(readLines(file), sep = "\n") # peek at file
gist_create(file, knit=TRUE)
## if not, GitHub doesn't allow upload of binary files via the HTTP API 
## (which gistr uses) - so see gist_create_git(), which uses git
file <- system.file("examples", "plots.Rmd", package = "gistr")
gist_create(file, knit=TRUE, imgur_inject = TRUE)
## works with ggplot2 as well
file <- system.file("examples", "ggplot_imgur.Rmd", package = "gistr")
gist_create(file, knit=TRUE)

# Render `.R` files
file <- system.file("examples", "example1.R", package = "gistr")
cat(readLines(file), sep = "\n") # peek at file
gist_create(file, knit = TRUE)
gist_create(file, knit = TRUE, include_source = TRUE)
## many files
(file1 <- system.file("examples", "example1.R", package = "gistr"))
(file2 <- system.file("examples", "example2.R", package = "gistr"))
cat(readLines(file1), sep = "\n") # peek at file
cat(readLines(file2), sep = "\n") # peek at file
gist_create(files=list(file1, file2), knit = TRUE)
### three at once, some .R and some .Rmd
file3 <- system.file("examples", "plots_imgur.Rmd", package = "gistr")
gist_create(files=list(file1, file2, file3), knit = TRUE)
gist_create(files=list(file1, file2, file3), knit = TRUE, 
  include_source = TRUE)

# Use rmarkdown::render instead of knitr::knit
file <- system.file("examples", "rmarkdown_eg.Rmd", package = "gistr")
gist_create(file, knit = TRUE, rmarkdown = TRUE, imgur_inject = TRUE,
   renderopts = list(output_format = "md_document"))
}
}
\seealso{
\code{\link[=gist_create_obj]{gist_create_obj()}}, \code{\link[=gist_create_git]{gist_create_git()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gist_map.R
\name{gist_map}
\alias{gist_map}
\title{Opens a full screen map after uploading a geojson file}
\usage{
gist_map(x, browse = TRUE)
}
\arguments{
\item{x}{An object of class \code{gist} generated by
\code{\link[=gist_create]{gist_create()}} or \code{\link[=gist_create_obj]{gist_create_obj()}}}

\item{browse}{Default: \code{TRUE}. Set to \code{FALSE} if you don't want to
automatically browse to the URL.}
}
\description{
Takes a gist object and a input geojson file name and renders fullscreen map
}
\examples{
\dontrun{
file <- system.file("examples", "ecoengine_eg.geojson", package = "gistr")
gist_id <- gist_create(file, browse = FALSE)
gist_map(gist_id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{update}
\alias{update}
\title{Update/modify a gist}
\usage{
update(gist, description = gist$description, ...)
}
\arguments{
\item{gist}{A gist object or something coerceable to a gist}

\item{description}{(character) Brief description of gist (optional)}

\item{...}{Curl options passed on to \code{\link[crul]{verb-GET}}}
}
\value{
an object of class \code{gist}
}
\description{
Update/modify a gist
}
\examples{
\dontrun{
file1 <- system.file("examples", "alm.md", package = "gistr")
file2 <- system.file("examples", "zoo.json", package = "gistr")

# add new files
gists(what = "minepublic")[[3]] \%>\%
 add_files(file1, file2) \%>\%
 update()

# update existing files
### file name has to match to current name
gists(what = "minepublic")[[3]] \%>\%
 update_files(file1) \%>\%
 update()

# delete existing files
### again, file name has to match to current name
gists(what = "minepublic")[[3]] \%>\%
 delete_files(file1, file2) \%>\%
 update()

# rename existing files
# For some reason, this operation has to upload the content too
### first name is old file name with path (must match), and second is 
### new file name (w/o path)
## add first
gists(what = "minepublic")[[3]] \%>\%
 add_files(file1, file2) \%>\%
 update()
## then rename
gists(what = "minepublic")[[3]] \%>\%
 rename_files(list(file1, "newfile.md")) \%>\%
 update()
### you can pass in many renames
gists(what = "minepublic")[[3]] \%>\%
 rename_files(list(file1, "what.md"), list(file2, "new.json")) \%>\%
 update()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embed.R
\name{embed}
\alias{embed}
\title{Get embed script for a gist}
\usage{
embed(gist)
}
\arguments{
\item{gist}{A gist object or something that can be coerced to a gist object.}
}
\description{
Get embed script for a gist
}
\examples{
\dontrun{
gists()[[1]] \%>\% embed()

# pass in a url
gist("https://gist.github.com/expersso/4ac33b9c00751fddc7f8") \%>\% embed
}
}
