ots
===



[![R-check](https://github.com/ropensci/ots/workflows/R-check/badge.svg)](https://github.com/ropensci/ots/actions?query=workflow%3AR-check)

`ots` is an R client to retrieve data from various ocean time series datasets, including:

* BATS
* HOT
* CALCOFI
* LTER Kelp
* UOPG
* more to come...

Jump over to the issues page to suggest data sets to include or comment on ongoing data source integration progress.

What's the point of getting data from the web in R? This way we only have to solve the problem of how to efficiently get a dataset once, then you can benefit from that. In addition, this should allow you to get any changes to the dataset that appear, or corrections. Last, getting data programatically in R should get you one step closer to a reproducible workflow, one that makes science easier primarily for yourself, and for others using your work.

## Install


```r
install.packages("devtools")
devtools::install_github("ropensci/ots")
```


```r
library("ots")
```

## Easy integration with dplyr


```r
library('dplyr')
tbl_df(bats("zooplankton")$data) %>% 
  filter(sieve_size > 1000) %>% 
  group_by(cruise) %>% 
  summarise(mean_water_vol = mean(water_vol))
```

## BATS - Zooplankton dataset


```r
bats("zooplankton")
```

## BATS - Production dataset


```r
bats("production")
```

## HOT dataset


```r
hot()
```

## Channels Islands National Park kelp data


```r
kelp("benthic_cover")
```

## CALCOFI data


```r
calcofi('hydro_cast')
```

## UOPG data

Various datasets available through this source - in this example getting data from Biowatt, and getting the meteorology data. Note that we still need to fix the column names...


```r
(biowatt_met <- uopg(dataset = 'biowatt', type = "meteorology"))
```

## More coming...

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ots/issues).
* License: MIT
* Get citation information for `ots` in R doing `citation(package = 'ots')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
ots 0.1.0
===============

### NEW FEATURES

* released to CRAN
I have read and agree to the the CRAN policies at 
http://cran.r-project.org/web/packages/policies.html

R CMD CHECK passed on my local OS X install with R 3.2.2 and
R development version, Ubuntu running on Travis-CI, and Windows
R 3.2.2 and devel on Win-Builder.

This is a new submission.

Thanks! 
Scott Chamberlain
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

Please note that the `ots` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
