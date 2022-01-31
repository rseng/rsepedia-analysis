# git2rdata <img src="man/figures/logo.svg" align="right" alt="git2rdata logo" width="120">

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/git2rdata)](https://cran.r-project.org/package=git2rdata)
[![Rdoc](https://www.rdocumentation.org/badges/version/git2rdata)](https://www.rdocumentation.org/packages/git2rdata)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://badges.ropensci.org/263_status.svg)](https://github.com/ropensci/software-review/issues/263)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![DOI](https://zenodo.org/badge/147685405.svg)](https://zenodo.org/badge/latestdoi/147685405)
[![codecov](https://codecov.io/gh/ropensci/git2rdata/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/git2rdata)
![GitHub forks](https://img.shields.io/github/forks/ropensci/git2rdata.svg?style=social)
![GitHub stars](https://img.shields.io/github/stars/ropensci/git2rdata.svg?style=social)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ropensci/git2rdata.svg)
![GitHub repo size](https://img.shields.io/github/repo-size/ropensci/git2rdata.svg)
<!-- badges: end -->

<p style="display:none">Please visit the git2rdata website at https://ropensci.github.io/git2rdata/. The vignette code on the website link to a rendered version of the vignette. Functions have a link to their help file.</p>

## Rationale

The `git2rdata` package is an R package for writing and reading dataframes as plain text files. 
A metadata file stores important information.

1. Storing metadata allows to maintain the classes of variables. 
By default, `git2rdata` optimizes the data for file storage. 
The optimization is most effective on data containing factors. 
The optimization makes the data less human readable.
The user can turn this off when they prefer a human readable format over smaller files.
Details on the implementation are available in `vignette("plain_text", package = "git2rdata")`.
1. Storing metadata also allows smaller row based [diffs](https://en.wikipedia.org/wiki/Diff) between two consecutive [commits](https://en.wikipedia.org/wiki/Commit_(version_control)). 
This is a useful feature when storing data as plain text files under version control. 
Details on this part of the implementation are available in `vignette("version_control", package = "git2rdata")`. 
Although we envisioned `git2rdata` with a [git](https://git-scm.com/) workflow in mind, you can use it in combination with other version control systems like [subversion](https://subversion.apache.org/) or [mercurial](https://www.mercurial-scm.org/).
1. `git2rdata` is a useful tool in a reproducible and traceable workflow. 
`vignette("workflow", package = "git2rdata")` gives a toy example.
1. `vignette("efficiency", package = "git2rdata")` provides some insight into the efficiency of file storage, git repository size and speed for writing and reading.

## Why Use Git2rdata?

- You can store dataframes as plain text files.
- The dataframe you read identical information content as the one you wrote.
    - No changes in data type.
    - Factors keep their original levels, including their order.
    - Date and date-time format are unambiguous, documented in the metadata.
- The data and the metadata are in a standard and open format, making it readable by other software.
- `git2rdata` checks the data and metadata during the reading. 
`read_vc()` informs the user if there is tampering with the data or metadata.
- Git2rdata integrates with the [`git2r`](https://cran.r-project.org/package=git2r) package for working with git repository from R.
    - Another option is using git2rdata solely for writing to disk and handle the plain text files with your favourite version control system outside of R.
- The optimization reduces the required disk space by about 30% for both the working directory and the git history. 
- Reading data from a HDD is 30% faster than `read.table()`, writing to a HDD takes about 70% more time than `write.table()`.
- Git2rdata is useful as a tool in a reproducible and traceable workflow. 
See `vignette("workflow", package = "git2rdata")`.
- You can detect when a file was last modified in the git history. 
Use this to check whether an existing analysis is obsolete due to new data. 
This allows to not rerun up to date analyses, saving resources.

## Talk About `git2rdata` at useR!2019 in Toulouse, France

<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/sbRPmakBFqo" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Installation

Install from CRAN

```r
install.packages("git2rdata")
```

Install the development version from GitHub

```r
# installation requires the "remotes" package
# install.package("remotes")

# install with vignettes (recommended)
remotes::install_github(
  "ropensci/git2rdata", 
  build = TRUE, 
  dependencies = TRUE, 
  build_opts = c("--no-resave-data", "--no-manual")
)
# install without vignettes
remotes::install_github("ropensci/git2rdata"))
```

## Usage in Brief

The user stores dataframes with `write_vc()` and retrieves them with `read_vc()`. 
Both functions share the arguments `root` and `file`. 
`root` refers to a base location where to store the dataframe. 
It can either point to a local directory or a local git repository. 
`file` is the file name to use and can include a path relative to `root`. 
Make sure the relative path stays within `root`.

```r
# using a local directory
library(git2rdata)
root <- "~/myproject" 
write_vc(my_data, file = "rel_path/filename", root = root)
read_vc(file = "rel_path/filename", root = root)
root <- git2r::repository("~/my_git_repo") # git repository
```

More details on store dataframes as plain text files in `vignette("plain_text", package = "git2rdata")`.

```r
# using a git repository
library(git2rdata)
repo <- repository("~/my_git_repo")
pull(repo)
write_vc(my_data, file = "rel_path/filename", root = repo, stage = TRUE)
commit(repo, "My message")
push(repo)
read_vc(file = "rel_path/filename", root = repo)
```

Please read `vignette("version_control", package = "git2rdata")` for more details on using git2rdata in combination with version control.

## What Data Sizes Can Git2rdata Handle?

The recommendation for git repositories is to use files smaller than 100 MiB, a repository size less than 1 GiB and less than 25k files. 
The individual file size is the limiting factor. 
Storing the airbag dataset ([`DAAG::nassCDS`](https://cran.r-project.org/package=DAAG)) with `write_vc()` requires on average 68 (optimized) or 97 (verbose) byte per record. 
The file reaches the 100 MiB limit for this data after about 1.5 million (optimized) or 1 million (verbose) observations. 

Storing a 90% random subset of the airbag dataset requires 370 kiB (optimized) or 400 kiB (verbose) storage in the git history. 
Updating the dataset with other 90% random subsets requires on average 60 kiB (optimized) to 100 kiB (verbose) per commit. 
The git history reaches the limit of 1 GiB after 17k (optimized) to 10k (verbose) commits.

Your mileage might vary.

## Citation

Please use the output of `citation("git2rdata")`

## Folder Structure

- `R`: The source scripts of the [R](https://cran.r-project.org/) functions with documentation in [Roxygen](https://CRAN.R-project.org/package=roxygen2) format
- `man`: The help files in [Rd](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Rd-format) format
- `inst/efficiency`: pre-calculated data to speed up `vignette("efficiency", package = "git2rdata")`
- `testthat`: R scripts with unit tests using the [testthat](https://CRAN.R-project.org/package=testthat) framework
- `vignettes`: source code for the vignettes describing the package
- `man-roxygen`: templates for documentation in Roxygen format
- `pkgdown`: source files for the `git2rdata` [website](https://ropensci.github.io/git2rdata/)
- `.github`: guidelines and templates for contributors

```
git2rdata
├── .github 
├─┬ inst
│ └── efficiency
├── man 
├── man-roxygen 
├── pkgdown
├── R
├─┬ tests
│ └── testthat
└── vignettes
```

## Contributions

`git2rdata` welcomes contributions. 
Please read our [Contributing guidelines](https://github.com/ropensci/git2rdata/blob/master/.github/CONTRIBUTING.md) first. 
The `git2rdata` project has a [Contributor Code of Conduct](https://github.com/ropensci/git2rdata/blob/master/.github/CODE_OF_CONDUCT.md). 
By contributing to this project, you agree to abide by its terms.

[![rOpenSci footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# git2rdata 0.3.1

* Use `icuSetCollate()` to define a standardised sorting.

# git2rdata 0.3.0

## New features

* `write_vc()` gains an optional `split_by` argument.
  See `vignette("split_by")` for more details.
* `rename_variable()` efficiently renames variables in a stored `git2rdata`
  object.

## Bugfixes

* `read_vc()`, `is_git2rdata()` and `is_git2rmeta()` now yield a better message
  when both the data and metadata are missing.

# git2rdata 0.2.2

* Use the [checklist](https://packages.inbo.be/checklist/) package for CI.

# git2rdata 0.2.1

## Bugfixes

* Explicitly use the `stringsAsFactors` of `data.frame()` in the examples and 
  unit tests if the dataframe contains characters. 
  The upcoming change in default value of `stringsAsFactors` requires this 
  change. 
  See https://developer.r-project.org/Blog/public/2020/02/16/stringsasfactors/index.html

# git2rdata 0.2.0

## BREAKING FEATURES

* Calculation of data hash has changed (#53). 
  You must use `upgrade_data()` to read data stored by an older version.
* `is_git2rdata()` and `upgrade_data()` no longer not test equality in data
  hashes (but `read_vc()` still does).
* `write_vc()` and `read_vc()` fail when `file` is a location outside of `root`
  (#50).
* Reordering factor levels requires `strict = TRUE`.

## Bugfixes

* Linux and Windows machines now generated the same data hash (#49).

## NEW FEATURES

* Internal sorting uses the "C" locale, regardless of the current locale.
* `read_vc()` reads older stored in an older version (#44). 
  When the version is too old, it prompts to `upgrade_data()`.
* Improve `warnings()` and `error()` messages.
* Use vector version of logo.

# git2rdata 0.1

* Transfer to rOpenSci.
* Use new logo (@peterdesmet, #37).
* Add estimate of upper bound of the number of commits.

# git2rdata 0.0.5

* `upgrade_data()` uses the same order of the metadata as `write_vc()`.

# git2rdata 0.0.4

## BREAKING FEATURES

* `write_vc()` stores the `git2rdata` version number to the metadata.
  Use `upgrade_data()` to update existing data.

## NEW FEATURES

* `read_vc()` checks the meta data hash. A mismatch results in an error.
* The meta data gains a data hash.
  A mismatch throws a warning when reading the object.
  This tolerates updating the data by other software, while informing the user
  that such change occurred.
* `is_git2rmeta()` validates metadata.
* `list_data()` lists files with valid metadata. 
* `rm_data()` and `prune_meta()` remove files with valid metadata. 
  They don't touch `tsv` file without metadata or `yml` files not associated
  with `git2rdata`.
* Files with invalid metadata yield a warning with `list_data()`, `rm_data()`
  and `prune_meta()`.

## Bugfixes

* `write_vc()` and `relabel()` handle empty strings (`''`) in characters and 
  factors (#24).
* `read_vc()` no longer treats `#` as a comment character.
* `read_vc()` handles non ASCII characters on Windows.

## Other changes

* Use a faster algorithm to detect duplicates (suggestion by @brodieG). 
* Improve documentation.
* Fix typo's in documentation, vignettes and README.
* Add a rOpenSci review badge to the README.
* The README mentions on upper bound on the size of dataframes.
* Set lifecycle to "maturing" and repo status to "active".
* The functions handle `root` containing regex expressions.
* Rework `vignette("workflow", package = "git2rdata")`.
* Update timings in `vignette("efficiency", package = "git2rdata")`
* Minor tweaks in `vignette("plain_text", package = "git2rdata")`

# git2rdata 0.0.3

* Fix typo's in documentation, vignettes and README.

# git2rdata 0.0.2

## BREAKING CHANGES

* `meta()` appends the metadata as a list to the objects rather than in YAML
  format.
* `yaml::write_yaml()` writes the metadata list in YAML format.
* `write_vc()` now uses the 'strict' argument instead of 'override'.
* `rm_data()` removes the data files. Use `prune_meta()` to remove left-over
  metadata files (#9).

## NEW FEATURES

* Vignette on [efficiency](https://ropensci.github.io/git2rdata/articles/efficiency.html) added (#2).
* Three separate vignettes instead of one large vignette.
    * Focus on the [plain text format](https://ropensci.github.io/git2rdata/articles/plain_text.html).
    * Focus on [version control](https://ropensci.github.io/git2rdata/articles/version_control.html).
    * Focus on [workflows](https://ropensci.github.io/git2rdata/articles/workflow.html).
* S3 methods replace the old S4 methods (#8).
* Optimized factors use stable indices. Adding or removing levels result in
  smaller diffs (#13).
* Use `relabel()` to alter factor levels without changing their index (#13).
* `write.table()` stores the raw data instead of `readr::write_tsv()` (#7).
  This avoids the `readr` dependency.
* `write_vc()` and `read_vc()` use the current working directory as default root
  (#6, @florisvdh).
* The user can specify a string to code missing values (default = `NA`).
  This allows the storage of the character string `"NA"`.
* `write_vc()` returns a list of issues which potentially result in large diffs.
* `list_data()` returns a vector with dataframes in the repository.

## Other changes

* `write_vc()` allows to use a custom `NA` string.
* Each helpfile contains a working example (#11).
* README updated (#12).
    * Updated the rationale with links to the vignettes.
    * `git2rdata` has a hexagon sticker logo.
    * Add the [![DOI](https://zenodo.org/badge/147685405.svg)](https://zenodo.org/badge/latestdoi/147685405).
    * The installation instructions use `remotes` and build the vignettes.
* We removed `auto_commit()` because of limited extra functionality over
  `git2r::commit()`.

# git2rdata 0.0.1

## NEW FEATURES

* Use `readr` to write and read plain text files.
* Allow storage of strings with "NA" or special characters.
* Handle ordered factors.
* Stop handling complex numbers.
## Test environments
* local
    * ubuntu 18.04.5 LTS, R 4.0.3
* github actions
    * macOS-latest, release
    * windows-latest, release
    * ubuntu 20.04, devel
    * ubuntu 16.04, oldrel
    * checklist package: ubuntu 20.04.1, R 4.0.3
* r-hub
    * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
    * Ubuntu Linux 16.04 LTS, R-release, GCC
    * Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 note

r-hub gave a few false positive notes

* Windows Server 2008 R2 SP1, R-devel, 32/64 bit

```
Possibly mis-spelled words in DESCRIPTION:
  rdata (28:22, 31:33, 36:20, 40:48, 41:20, 43:24, 44:62, 45:62)
  workflow (41:37, 44:15, 44:36)
```

* Fedora Linux, R-devel, clang, gfortran

```
Possibly mis-spelled words in DESCRIPTION:
  rdata (28:22, 31:33, 36:20, 40:48, 41:20, 43:24, 44:62, 45:62)
```

Ubuntu Linux 16.04 LTS, R-release, GCC failed on r-hub because ICU is not
available on that build.

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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
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

Please note that the git2rdata project is released with a
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
---
name: Bug Report
about: Create a report to help us improve
---
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- Please paste your sessioninfo::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

```r
# replace this with your reproducible example
```

<details> <summary><strong>Session Info</strong></summary>
```r
# replace this by with the output from sessioninfo::session_info() or sessionInfo()
```
</details>
---
name: Feature Request
about: Suggest an idea for this project
---
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- Please describe the feature you propose as detailed as possible. 
If possible, add some code examples below indicating how the updated code should work.
Consider writing it as a <a href = "https://testthat.r-lib.org/">testthat</a> unit test-->

```r
# a silly feature request, fails under the current implementation
b <- 1 + 1
stopifnot(all.equal(b == 3))
```

```r
# testthat version of the silly feature request
expect_equal(1 + 1, 3)
```
---
title: "Optimizing Storage for Version Control"
author: "Thierry Onkelinx"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Optimizing Storage for Version Control}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(width = 83)
```

## Introduction

This vignette focuses on what `git2rdata` does to make storing dataframes under version control more efficient and convenient. 
`vignette("plain_text", package = "git2rdata")` describes all details on the actual file format. 
Hence we will not discuss the `optimize` and `na` arguments to the `write_vc()` function.

We will not illustrate the efficiency of `write_vc()` and `read_vc()`.
`vignette("efficiency", package = "git2rdata")` covers those topics.

## Setup

```{r initialise}
# Create a directory in tempdir
root <- tempfile(pattern = "git2r-")
dir.create(root)
# Create dummy data
set.seed(20190222)
x <- data.frame(
  x = sample(LETTERS),
  y = factor(
    sample(c("a", "b", NA), 26, replace = TRUE),
    levels = c("a", "b", "c")
  ),
  z = c(NA, 1:25),
  abc = c(rnorm(25), NA),
  def = sample(c(TRUE, FALSE, NA), 26, replace = TRUE),
  timestamp = seq(
    as.POSIXct("2018-01-01"),
    as.POSIXct("2019-01-01"),
    length = 26
  ),
  stringsAsFactors = FALSE
)
str(x)
```

## Assumptions

A critical assumption made by `git2rdata` is that the dataframe itself contains all information. 
Each row is an observation, each column is a variable. 
The dataframe has `colnames` but no `rownames`. 
This implies that two observations switching place does not alter the information content. 
Nor does switching two variables.

Version control systems like [git](https://git-scm.com/), [subversion](https://subversion.apache.org/) or [mercurial](https://www.mercurial-scm.org/) focus on accurately keeping track of _any_ change in the files. 
Two observations switching place in a plain text file _is_ a change, although the information content^[_sensu_ `git2rdata`] doesn't change.
`git2rdata` helps the user to prepare the plain text files in such a way that any change in the version history is an actual change in the information content.

## Sorting Observations

Version control systems often track changes in plain text files based on row based differences. 
In layman's terms they record lines removed from and inserted in the file at what location. 
Changing an existing line implies removing the old version and inserting the new one. 
The minimal example below illustrates this.

Original version

```
A,B
1,10
2,11
3,12
```

Altered version. 
The row containing `1, 10` moves to the last line. 
The row containing `3,12` changed to `3,0`.

```
A,B
2,11
3,0
1,10
```

Diff between original and altered version. Notice than we have a deletion of two lines and two insertions.

```diff
A,B
-1,10
2,11
-3,12
+3,0
+1,10
```

Ensuring that the observations are always sorted in the same way thus helps minimizing the diff. The sorted version of the same altered version looks like the example below. 

```
A,B
1,10
2,11
3,0
```

Diff between original and the sorted alternate version. Notice that all changes revert to actual changes in the information content. Another benefit is that changes are easily spotted in the diff. A deletion without insertion on the next line is a removed observation. An insertion without preceding deletion is a new observation. A deletion followed by an insertion is an updated observation.

```diff
A,B
1,10
2,11
-3,12
+3,0
```

This is where the `sorting` argument comes into play. 
If this argument is not provided when writing a file for the first time, it will yield a warning about the lack of sorting. 
`write_vc()` then writes the observations in their current order. 
New versions of the file will not apply any sorting either, leaving this burden to the user. 
The changed hash for the data file illustrates this in the example below. 
The metadata hash remains the same.

```{r row_order}
library(git2rdata)
write_vc(x, file = "row_order", root = root)
write_vc(x[sample(nrow(x)), ], file = "row_order", root = root)
```

`sorting` should contain a vector of variable names. 
The observations are automatically sorted along these variables. 
Now we get an error because the set of sorting variables has changed. 
The metadata stores the set of sorting variables. 
Changing the sorting can potentially lead to large diffs, which `git2rdata` tries to avoid as much as possible.

From this moment on we will store the output of `write_vc()` in an object reduce output.

```{r apply_sorting, error = TRUE}
fn <- write_vc(x, "row_order", root, sorting = "y")
```

Using `strict = FALSE` turns such errors into warnings and allows to update the file. Notice that we get a new warning: the variable we used for sorting resulted in ties, thus the order of the observations is not guaranteed to be stable. The solution is to use more or different variables. We'll need `strict = FALSE` again to override the change in sorting variables.

```{r update_sorting}
fn <- write_vc(x, "row_order", root, sorting = "y", strict = FALSE)
fn <- write_vc(x, "row_order", root, sorting = c("y", "x"), strict = FALSE)
```

Once we have defined the sorting, we may omit the `sorting` argument when writing new versions. 
`write_vc` uses the sorting as defined in the existing metadata.
It checks for potential ties.
Ties results in a warning.

```{r update_sorted}
print_file <- function(file, root, n = -1) {
  fn <- file.path(root, file)
  data <- readLines(fn, n = n)
  cat(data, sep = "\n")
}
print_file("row_order.yml", root, 7)
fn <- write_vc(x[sample(nrow(x)), ], "row_order", root)
fn <- write_vc(x[sample(nrow(x)), ], "row_order", root, sorting = c("y", "x"))
fn <- write_vc(x[sample(nrow(x), replace = TRUE), ], "row_order", root)
```

## Sorting Variables

The order of the variables (columns) has an even bigger impact on a row based diff. Let's revisit our minimal example. Suppose that we swap `A` and `B` from our [original example](#sorting-observations). The new data looks as below.

```
B,A
10,1
11,2
13,3
```

The resulting diff is maximal because every single row changed. 
Yet none of the information changed. 
Hence, maintaining column order is crucial when storing dataframes as plain text files under version control. 
The `vignette("efficiency", package = "git2rdata")` vignette illustrates this on a more realistic data set.

```diff
-A,B
+B,A
-1,10
+10,1
-2,11
+11,2
-3,13
+13,3
```

When `write_vc()` writes a dataframe for the first time, it stores the original order of the columns in the metadata.
From that moment on, `write_vc()` uses the order stored in the metadata. 
The example below writes the same data set twice. 
The second version contains identical information but randomizes the order of the observations and the columns. 
The sorting by the internals of `write_vc()` will undo this randomization, resulting in an unchanged file.

```{r variable_order}
write_vc(x, "column_order", root, sorting = c("x", "abc"))
print_file("column_order.tsv", root, n = 5)
write_vc(x[sample(nrow(x)), sample(ncol(x))], "column_order", root)
print_file("column_order.tsv", root, n = 5)
```

## Handling Factors Optimized

`vignette("plain_text", package = "git2rdata")` and `vignette("efficiency", package = "git2rdata")` illustrate how we can store a factor more efficiently when storing their index in the data file and the indices and labels in the metadata. 
We take this even a bit further: what happens if new data arrives and we need an extra factor level? 

```{r factor}
old <- data.frame(color = c("red", "blue"), stringsAsFactors = TRUE)
write_vc(old, "factor", root, sorting = "color")
print_file("factor.yml", root)
```

Let's add an observation with a new factor level. If we store the updated dataframe in a new file, we see that the indices are different. The factor level `"blue"` remains unchanged, but `"red"` becomes the third level and get index `3` instead of index `2`. This could lead to a large diff whereas the potential semantics (and thus the information content) are not changed.

```{r factor2}
updated <- data.frame(
  color = c("red", "green", "blue"), 
  stringsAsFactors = TRUE
)
write_vc(updated, "factor2", root, sorting = "color")
print_file("factor2.yml", root)
```

When we try to overwrite the original data with the updated version, we get an error because there is a change in factor levels and / or indices. In this specific case, we decided that the change is OK and force the writing by setting `strict = FALSE`. Notice that the original labels (`"blue"` and `"red"`) keep their index, the new level (`"green"`) gets the first available index number.

```{r factor_update, error = TRUE}
write_vc(updated, "factor", root)
fn <- write_vc(updated, "factor", root, strict = FALSE)
print_file("factor.yml", root)
```

The next example removes the `"blue"` level and switches the order of the remaining levels. 
Notice that the meta data retains the existing indices. 
The order of the labels and indices reflects their new ordering.

```{r factor_deleted}
deleted <- data.frame(
  color = factor(c("red", "green"), levels = c("red", "green"))
)
write_vc(deleted, "factor", root, sorting = "color", strict = FALSE)
print_file("factor.yml", root)
```

Changing a factor to an ordered factor or _vice versa_ will also keep existing level indices.

```{r factor_ordered}
ordered <- data.frame(
  color = factor(c("red", "green"), levels = c("red", "green"), ordered = TRUE)
)
write_vc(ordered, "factor", root, sorting = "color", strict = FALSE)
print_file("factor.yml", root)
```

## Relabelling a Factor

The example below will store a dataframe, relabel the factor levels and store it again using `write_vc()`. 
Notice the update of both the labels and the indices. 
Hence creating a large diff, where updating the labels would do.

```{r}
write_vc(old, "write_vc", root, sorting = "color")
print_file("write_vc.yml", root)
relabeled <- old
# translate the color names to Dutch
levels(relabeled$color) <- c("blauw", "rood")
write_vc(relabeled, "write_vc", root, strict = FALSE)
print_file("write_vc.yml", root)
```

We created `relabel()`, which changes the labels in the meta data while maintaining their indices. 
It takes three arguments: the name of the data file, the root and the change. 
`change` accepts two formats, a list or a dataframe. 
The name of the list must match with the variable name of a factor in the data. 
Each element of the list must be a named vector, the name being the existing label and the value the new label. 
The dataframe format requires a `factor`, `old` and `new` variable with one row for each change in label. 

```{r}
write_vc(old, "relabel", root, sorting = "color")
relabel("relabel", root, change = list(color = c(red = "rood", blue = "blauw")))
print_file("relabel.yml", root)
relabel(
  "relabel", root, 
  change = data.frame(
    factor = "color", old = "blauw", new = "blue", stringsAsFactors = TRUE
  )
)
print_file("relabel.yml", root)
```

A _caveat_: `relabel()` does not make sense when the data file uses verbose storage.
The verbose mode stores the factor labels and not their indices, in which case relabelling a label will always yield a large diff. 
Hence, `relabel()` requires the optimized storage. 
---
title: "Storing Large Dataframes"
author: "Thierry Onkelinx"
output: 
  rmarkdown::html_vignette:
        fig_caption: yes
vignette: >
  %\VignetteIndexEntry{Storing Large Dataframes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{git2r}
  %\VignetteDepends{microbenchmark}
  %\VignetteDepends{ggplot2}
---

```{r setup, include = FALSE}
library(knitr)
opts_chunk$set(
  fig.height = 4, fig.width = 6,
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)
inbo_colours <- c("#959B38", "#729BB7", "#E87837", "#BDDDD7", "#E4E517",
                  "#843860", "#C04384", "#C2C444", "#685457")
theme_inbo <- function(base_size = 12, base_family = "") {
  rect_bg <- "white"
  legend_bg <- "white"
  panel_bg <- "#F3F3F3"
  panel_grid <- "white"
  plot_bg <- "white"
  half_line <- base_size / 2
  theme(
    line = element_line(colour = "black", size = 0.5, linetype = 1,
                        lineend = "butt"),
    rect = element_rect(fill = rect_bg, colour = "black", size = 0.5,
                        linetype = 1),
    text = element_text(family = base_family, face = "plain",
                        colour = "#843860", size = base_size, hjust = 0.5,
                        vjust = 0.5, angle = 0, lineheight = 0.9,
                        margin = margin(), debug = FALSE),
    axis.line = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text = element_text(size = rel(0.8)),
    axis.text.x = element_text(margin = margin(t = 0.8 * half_line / 2),
                               vjust = 1),
    axis.text.x.top = NULL,
    axis.text.y = element_text(margin = margin(r = 0.8 * half_line / 2),
                               hjust = 1),
    axis.text.y.right = NULL,
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(colour = "black"),
    axis.title.x = element_text(
      margin = margin(t = 0.8 * half_line, b = 0.8 * half_line / 2)
    ),
    axis.title.x.top = NULL,
    axis.title.y = element_text(
      margin = margin(r = 0.8 * half_line, l = 0.8 * half_line / 2),
      angle = 90
    ),
    axis.title.y.right = NULL,
    legend.background = element_rect(colour = NA, fill = legend_bg),
    legend.key = element_rect(fill = panel_bg, colour = NA),
    legend.key.size = unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.margin = NULL,
    legend.spacing = unit(0.2, "cm"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.text = element_text(size = rel(0.8)),
    legend.text.align = NULL,
    legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0,
                                colour = "black"),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = margin(t = half_line, r = half_line, b = half_line,
                               l = half_line),
    legend.box.background = element_rect(colour = NA, fill = legend_bg),
    legend.box.spacing = unit(0.2, "cm"),
    panel.background = element_rect(fill = panel_bg, colour = NA),
    panel.border = element_blank(),
    panel.grid = element_line(colour = panel_grid),
    panel.grid.minor = element_line(colour = panel_grid, size = 0.25),
    panel.spacing = unit(half_line, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = FALSE,
    strip.background = element_rect(fill = "#8E9DA7", colour = NA),
    strip.text = element_text(size = rel(0.8), colour = "#F3F3F3"),
    strip.text.x = element_text(margin = margin(t = half_line, b = half_line)),
    strip.text.y = element_text(margin = margin(r = half_line, l = half_line),
                                angle = -90),
    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"),
    strip.placement = "outside",
    plot.background = element_rect(colour = NA, fill = plot_bg),
    plot.title = element_text(size = rel(1.2),
                              margin = margin(0, 0, half_line, 0)),
    plot.subtitle = element_text(size = rel(1),
                                 margin = margin(0, 0, half_line, 0)),
    plot.caption = element_text(size = rel(0.6),
                                margin = margin(0, 0, half_line, 0)),
    plot.margin = margin(t = half_line, r = half_line, b = half_line,
                         l = half_line),
    plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5),
    plot.tag.position = "topleft",
    complete = TRUE
  )
}
theme_set(theme_inbo())
update_geom_defaults("line", list(colour = "#356196"))
update_geom_defaults("hline", list(colour = "#356196"))
update_geom_defaults("boxplot", list(colour = "#356196"))
update_geom_defaults("smooth", list(colour = "#356196"))
```

## Introduction

Sometimes, a large dataframe has one or more variables with a small number of unique combinations.
E.g. a dataframe with one or more factor variables.
Storing the entire dataframe as a single text file requires storing lots of replicated data.
Each row stores the information for every variable, even if a subset of these variables remains constant over a subset of the data.

In such a case we can use the `split_by` argument of `write_vc()`.
This will store the large dataframe over a set of tab separated files.
One file for every combination of the variables defined by `split_by`.
Every partial data file holds the other variables for one combination of `split_by`.
We remove the `split_by` variables from the partial data files, reducing their size.
We add an `index.tsv` containing the combinations of the `split_by` variables and a unique hash for each combination.
This hash becomes the base name of the partial data files.

Splitting the dataframe into smaller files makes them easier to handle in version control system.
The total size depends on the amount of replication in the dataframe.
More on that in the next section.

## When to Split the Dataframe

Let's set the following variables:

-   $s$: the average number of bytes to store a single line of the `split_by` variables.

-   $r$: the average number of bytes to store a single line of the remaining variables.

-   $h_s$: the number of bytes to store the header of the `split_by` variables.

-   $h_r$: the number of bytes to store the header of the remaining variables.

-   $N$: the number of rows in the dataframe.

-   $N_s$: the number of unique combinations of the `split_by` variables.

Storing the dataframe with `write_vc()` without `split_by` requires $h_s + h_r + 1$ bytes for the header and $s + r + 1$ bytes for every observation.
The total number of bytes is $T_0 = h_s + h_r + 1 + N (s + r + 1)$.
Both $+ 1$ originate from the tab character to separate the `split_by` variables from the remaining variables.

Storing the dataframe with `write_vc()` with `split_by` requires an index file to store the combinations of the `split_by` variables.
It will use $h_s$ bytes for the header and $N_s s$ for the data.
The headers of the partial data files require $N_s h_r$ bytes ($N_s$ files and $h_r$ byte per file).
The data in the partial data files require $N r$ bytes.
The total number of bytes is $T_s = h_s + N_s s + N_s h_r + N r$.

We can look at the ratio of $T_s$ over $T_0$.

$$\frac{T_s}{T_0} = \frac{h_s + N_s s + N_s h_r + N r}{h_s + h_r + 1 + N (s + r + 1)}$$

Let's simplify the equation by assuming that we need an equal amount of character for the headers and the data ($h_s = s$ and $h_r = r$).

$$\frac{T_s}{T_0} = \frac{s + N_s s + N_s r + N r}{s + r + 1 + N (s + r + 1)}$$

$$\frac{T_s}{T_0} = \frac{s + N_s s + N_s r + N r}{s + r + 1 + N s + N r + N}$$

Let assume that $s = a r$ with $0 < a$ and $N_s = b N$ with $0 < b < 1$.

$$\frac{T_s}{T_0} = \frac{a r + N a b r + N b r + N r}{a r + r + 1 + N a r + N r + N}$$

$$\frac{T_s}{T_0} = \frac{(a + N a b + N b + N) r}{(N + 1) (a r + r + 1)}$$

$$\frac{T_s}{T_0} = \frac{a + N a b + N b + N}{(N + 1) (a + 1 + 1 / r)}$$ $$\frac{T_s}{T_0} = \frac{a + (a b + b + 1) N }{(N + 1) (a + 1 + 1 / r)}$$

When $N$ is large, we can state that $a \lll N$ and $N / (N + 1) \approx 1$.

$$\frac{T_s}{T_0} \approx \frac{a b + b + 1}{a + 1 + 1 / r}$$

```{r ratio, fig.cap = "Storage space required using `split_by` relative to storing a single file.", echo = FALSE}
combinations <- expand.grid(
  a = c(0.25, 0.5, 1, 2, 4),
  b = seq(0, 1, length = 41),
  r = c(10, 100, 1000)
)
combinations$ratio <- with(
  combinations,
  (a * b + b + 1) / (a + 1 + 1 / r)
)
ggplot(combinations, aes(x = b, y = ratio, colour = factor(a))) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line() +
  facet_wrap(~ paste("r =", r)) +
  scale_x_continuous(
    expression(b~{"="}~N[s]~{"/"}~N),
    labels = function(x) {
      paste0(100 * x, "%")
    }
  ) +
  scale_y_continuous(
    "Relative amount of disk space",
    labels = function(x) {
      paste0(100 * x, "%")
    }
  ) +
  scale_colour_manual(
    "a = s / r",
    values = inbo_colours,
    labels = c("1/4", "1/2", "1", "2", "4")
  )
```

The figure illustrates that using `split_by` is more efficient when the number of unique combinations ($N_s$) of the `split_by` variables is much smaller than the number of rows in the dataframe ($N$).
The efficiency also increases when the storage for a single combination of `split_by` variables ($s$) is larger than the storage needed for a single line of the remain variables ($r$).
The storage needed for a single line of the remain variables ($r$) doesn't influence the efficiency.

## Benchmarking

```{r load_data, echo = FALSE}
airbag <- readRDS(
  system.file("efficiency", "airbag.rds", package = "git2rdata")
)
```

```{r set_tmp_dir}
library(git2rdata)
root <- tempfile("git2rdata-split-by")
dir.create(root)
```

```{r get_write_timings, eval = system.file("split_by", "write_timings.rds", package = "git2rdata") == ""}
library(microbenchmark)
mb <- microbenchmark(
  part_1 = write_vc(airbag, "part_1", root, sorting = "X"),
  part_2 = write_vc(airbag, "part_2", root, sorting = "X", split_by = "airbag"),
  part_3 = write_vc(airbag, "part_3", root, sorting = "X", split_by = "abcat"),
  part_4 = write_vc(
    airbag, "part_4", root, sorting = "X", split_by = c("airbag", "sex")
  ),
  part_5 = write_vc(airbag, "part_5", root, sorting = "X", split_by = "dvcat"),
  part_6 = write_vc(
    airbag, "part_6", root, sorting = "X", split_by = "yearacc"
  ),
  part_15 = write_vc(
    airbag, "part_15", root, sorting = "X", split_by = c("dvcat", "abcat")
  ),
  part_45 = write_vc(
    airbag, "part_45", root, sorting = "X", split_by = "yearVeh"
  ),
  part_270 = write_vc(
    airbag, "part_270", root, sorting = "X", split_by = c("yearacc", "yearVeh")
  )
)
mb$time <- mb$time / 1e6
```

```{r store_write_timings, echo = FALSE}
if (system.file("split_by", "write_timings.rds", package = "git2rdata") == "") {
  dir.create(file.path("..", "inst", "split_by"), showWarnings = FALSE)
  saveRDS(mb, file.path("..", "inst", "split_by", "write_timings.rds"))
} else {
  mb <- readRDS(
    system.file("split_by", "write_timings.rds", package = "git2rdata")
  )
}
```

Splitting the dataframe over more than one file takes more time to write the data.
The log time seems to increase quadratic with log number of parts.

```{r plot_write_timings, echo = FALSE, fig.cap = "Boxplot of the write timings for different number of parts."}
mb$combinations <- as.integer(gsub("part_", "", levels(mb$expr)))[mb$expr]
ggplot(mb, aes(x = combinations, y = time)) +
  geom_boxplot(aes(group = combinations)) +
  scale_x_log10("Number of parts") +
  scale_y_log10("Time (in milliseconds)")
```

```{r get_read_timings, eval = system.file("split_by", "read_timings.rds", package = "git2rdata") == ""}
mb_r <- microbenchmark(
  part_1 = read_vc("part_1", root),
  part_2 = read_vc("part_2", root),
  part_3 = read_vc("part_3", root),
  part_4 = read_vc("part_4", root),
  part_5 = read_vc("part_5", root),
  part_6 = read_vc("part_6", root),
  part_15 = read_vc("part_15", root),
  part_45 = read_vc("part_45", root),
  part_270 = read_vc("part_270", root)
)
mb_r$time <- mb_r$time / 1e6
```

```{r store_read_timings, echo = FALSE}
if (system.file("split_by", "read_timings.rds", package = "git2rdata") == "") {
  saveRDS(mb_r, file.path("..", "inst", "split_by", "read_timings.rds"))
} else {
  mb_r <- readRDS(
    system.file("split_by", "read_timings.rds", package = "git2rdata")
  )
}
```

A small number of parts does not seem to affect the read timings much.
Above ten parts, the required time for reading seems to increase.
The log time seems to increase quadratic with log number of parts.

```{r plot_read_timings, echo = FALSE, fig.cap = "Boxplot of the read timings for the different number of parts."}
mb_r$combinations <- as.integer(gsub("part_", "", levels(mb_r$expr)))[mb_r$expr]
ggplot(mb_r, aes(x = combinations, y = time)) +
  geom_boxplot(aes(group = combinations)) +
  scale_x_log10("Number of parts") +
  scale_y_log10("Time (in milliseconds)")
```
---
title: "Efficiency Relative to Storage and Time"
author: "Thierry Onkelinx"
output: 
  rmarkdown::html_vignette:
        fig_caption: yes
vignette: >
  %\VignetteIndexEntry{Efficiency Relative to Storage and Time}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{git2r}
  %\VignetteDepends{microbenchmark}
  %\VignetteDepends{ggplot2}
---

```{r setup, include = FALSE}
library(knitr)
opts_chunk$set(
  fig.height = 4, fig.width = 6,
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)
inbo_colours <- c("#959B38", "#729BB7", "#E87837", "#BDDDD7", "#E4E517", 
                  "#843860", "#C04384", "#C2C444", "#685457")
theme_inbo <- function(base_size = 12, base_family = "") {
  rect.bg <- "white"
  legend.bg <- "white"
  panel.bg <- "#F3F3F3"
  panel.grid <- "white"
  plot.bg <- "white"
  half_line <- base_size / 2
  theme(
    line = element_line(colour = "black", size = 0.5, linetype = 1, 
                        lineend = "butt"),
    rect = element_rect(fill = rect.bg, colour = "black", size = 0.5, 
                        linetype = 1),
    text = element_text(family = base_family, face = "plain", 
                        colour = "#843860", size = base_size, hjust = 0.5, 
                        vjust = 0.5, angle = 0, lineheight = 0.9, 
                        margin = margin(), debug = FALSE),
    axis.line = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text = element_text(size = rel(0.8)),
    axis.text.x = element_text(margin = margin(t = 0.8 * half_line / 2), 
                               vjust = 1),
    axis.text.x.top = NULL,
    axis.text.y = element_text(margin = margin(r = 0.8 * half_line / 2), 
                               hjust = 1),
    axis.text.y.right = NULL,
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(colour = "black"),
    axis.title.x = element_text(
      margin = margin(t = 0.8 * half_line, b = 0.8 * half_line / 2)
    ),
    axis.title.x.top = NULL,
    axis.title.y = element_text(
      margin = margin(r = 0.8 * half_line, l = 0.8 * half_line / 2),
      angle = 90
    ),
    axis.title.y.right = NULL,
    legend.background = element_rect(colour = NA, fill = legend.bg),
    legend.key = element_rect(fill = panel.bg, colour = NA),
    legend.key.size = unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.margin = NULL,
    legend.spacing = unit(0.2, "cm"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.text = element_text(size = rel(0.8)),
    legend.text.align = NULL,
    legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0, 
                                colour = "black"),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = margin(t = half_line, r = half_line, b = half_line, 
                               l = half_line),
    legend.box.background = element_rect(colour = NA, fill = legend.bg),
    legend.box.spacing = unit(0.2, "cm"),
    panel.background = element_rect(fill = panel.bg, colour = NA),
    panel.border = element_blank(),
    panel.grid = element_line(colour = panel.grid),
    panel.grid.minor = element_line(colour = panel.grid, size = 0.25),
    panel.spacing = unit(half_line, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = FALSE,
    strip.background = element_rect(fill = "#8E9DA7", colour = NA),
    strip.text = element_text(size = rel(0.8), colour = "#F3F3F3"),
    strip.text.x = element_text(margin = margin(t = half_line, b = half_line)),
    strip.text.y = element_text(margin = margin(r = half_line, l = half_line),
                                angle = -90),
    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"),
    strip.placement = "outside",
    plot.background = element_rect(colour = NA, fill = plot.bg),
    plot.title = element_text(size = rel(1.2), 
                              margin = margin(0, 0, half_line, 0)),
    plot.subtitle = element_text(size = rel(1),
                                 margin = margin(0, 0, half_line, 0)),
    plot.caption = element_text(size = rel(0.6),
                                margin = margin(0, 0, half_line, 0)),
    plot.margin = margin(t = half_line, r = half_line, b = half_line, 
                         l = half_line),
    plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5),
    plot.tag.position = "topleft",
    complete = TRUE
  )
}
theme_set(theme_inbo())
update_geom_defaults("line", list(colour = "#356196"))
update_geom_defaults("hline", list(colour = "#356196"))
update_geom_defaults("boxplot", list(colour = "#356196"))
```

## Introduction

This vignette compares storage and retrieval of data by `git2rdata` with other standard R functionality. We consider `write.table()` and `read.table()` for data stored in a plain text format. `saveRDS()` and `readRDS()` use a compressed binary format.

To get some meaningful results, we will use the `nassCDS` dataset from the [DAAG](https://www.rdocumentation.org/packages/DAAG/versions/1.22/topics/nassCDS) package. 
We'll avoid the dependency on the package by directly downloading the data.

```{r download_data, eval = system.file("efficiency", "airbag.rds", package = "git2rdata") == ""}
airbag <- read.csv(
  "https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/nassCDS.csv"
)
airbag$dead <- airbag$dead == "dead"
airbag$airbag <- airbag$airbag == "airbag"
airbag$seatbelt <- airbag$seatbelt == "belted"
airbag$dvcat <- as.ordered(airbag$dvcat)
```

```{r load_data, echo = FALSE}
if (system.file("efficiency", "airbag.rds", package = "git2rdata") == "") {
  saveRDS(airbag, file.path("..", "inst", "efficiency", "airbag.rds"))
} else {
  airbag <- readRDS(
    system.file("efficiency", "airbag.rds", package = "git2rdata")
  )
}
```

```{r data_structure}
str(airbag)
```

## Data Storage 

### On a File System

We start by writing the dataset as is with `write.table()`, `saveRDS()`, `write_vc()` and `write_vc()` without storage optimization. Note that `write_vc()` uses optimization by default. Since `write_vc()` creates two files for each data set, we take their combined file size into account.

```{r set_tmp_dir}
library(git2rdata)
root <- tempfile("git2rdata-efficient")
dir.create(root)
```

```{r file_size}
write.table(airbag, file.path(root, "base_R.tsv"), sep = "\t")
base_size <- file.size(file.path(root, "base_R.tsv"))

saveRDS(airbag, file.path(root, "base_R.rds"))
rds_size <- file.size(file.path(root, "base_R.rds"))

fn <- write_vc(airbag, "airbag_optimize", root, sorting = "X")
optim_size <- sum(file.size(file.path(root, fn)))

fn <- write_vc(airbag, "airbag_verbose", root, sorting = "X", optimize = FALSE)
verbose_size <- sum(file.size(file.path(root, fn)))
```

Since the data is highly compressible, `saveRDS()` yields the smallest file at the cost of having a binary file format. Both `write_vc()` formats yield smaller files than `write.table()`. 
Partly because `write_vc()` doesn't store row names and doesn't use quotes unless needed. 
The difference between the optimized and verbose version of `write_vc()` is, in this case, solely due to the way `write_vc()` stores factors in the data (`tsv`) file. 
The optimized version stores the indices of the factor whereas the verbose version stores the levels. 
For example: `airbag$dvcat` has 5 levels with short labels (on average 5 character), storing the index requires 1 character. 
This results in more compact files.

```{r table_file_size, echo = FALSE}
kable(
  data.frame(
    method = c("saveRDS()", "write_vc(), optimized", "write_vc(), verbose", 
               "write.table()"),
    file_size = c(rds_size, optim_size, verbose_size, base_size) / 2 ^ 10,
    relative = c(rds_size, optim_size, verbose_size, base_size) / base_size
  ),
  caption = "Resulting file sizes (in kB) and file sizes relative to the size of write.table().",
  digits = 2
)
```

The reduction in file size when storing in factors depends on the length of the labels, the number of levels and the number of observations. 
The figure below illustrates the strong gain as soon as the level labels contain more than two characters. 
The gain is less pronounced when the factor has a large number of levels. 
The optimization fails in extreme cases with short factor labels and a high number of levels.

```{r factor_label_length, echo = FALSE, fig.cap = "Effect of the label length on the efficiency of storing factor optimized, assuming 1000 observations", warning = FALSE}
ratio <- function(label_length = 1:20, n_levels = 9, n_obs = 1000) {
  meta_length <- 63 + (5 + label_length + floor(log10(n_levels))) * n_levels
  optimized <- n_obs * mean(ceiling(log10(seq_len(n_levels) + 1)))
  verbose <- n_obs * label_length
  ifelse(
    62 ^ label_length >= n_levels,
    (optimized + meta_length) / (verbose + meta_length),
    NA
  )
}
lengths <- 1:50
f_ratio <- rbind(
  data.frame(
    label_length = lengths,
    levels = 1000,
    ratio = ratio(lengths, 1000, 1000)
  ),
  data.frame(
    label_length = lengths,
    levels = 100,
    ratio = ratio(lengths, 100, 1000)
  ),
  data.frame(
    label_length = lengths,
    levels = 10,
    ratio = ratio(lengths, 10, 1000)
  ),
  data.frame(
    label_length = lengths,
    levels = 3,
    ratio = ratio(lengths, 3, 1000)
  )
)
f_ratio$levels <- factor(f_ratio$levels, levels = c(1000, 100, 10, 3))
ggplot(f_ratio, aes(x = label_length, y = ratio, colour = levels)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line() +
  scale_x_continuous("label length (characters)") +
  scale_y_continuous("optimized bytes / verbose bytes", 
                     breaks = seq(0, 1.25, by = 0.25)) +
  scale_colour_manual("number of \nlevels", values = inbo_colours)
```

The effect of the number of observations is mainly due to the overhead of storing the metadata. The importance of this overhead increases when the number of observations is small.

```{r factor_observations, echo = FALSE, fig.cap = "Effect of the number of observations on the efficiency of storing factor optimized assuming labels with 10 characters"}
n_obs <- 10 ^ seq(log10(1), log10(10000), length = 41)
f_ratio <- rbind(
  data.frame(
    observations = n_obs,
    levels = 3,
    ratio = sapply(n_obs, ratio, label_length = 10, n_levels = 3)
  ),
  data.frame(
    observations = n_obs,
    levels = 10,
    ratio = sapply(n_obs, ratio, label_length = 10, n_levels = 10)
  ),
  data.frame(
    observations = n_obs,
    levels = 100,
    ratio = sapply(n_obs, ratio, label_length = 10, n_levels = 100)
  ),
  data.frame(
    observations = n_obs,
    levels = 1000,
    ratio = sapply(n_obs, ratio, label_length = 10, n_levels = 1000)
  )
)
f_ratio$levels <- factor(f_ratio$levels, levels = c(1000, 100, 10, 3))
ggplot(f_ratio, aes(x = observations, y = ratio, colour = levels)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line() +
  scale_x_log10() +
  scale_y_continuous("optimized bytes / verbose bytes", 
                     breaks = seq(0, 1.25, by = 0.25)) +
  scale_colour_manual("number of \nlevels", values = inbo_colours)
```

### In Git Repositories

Here we will simulate how much space the data requires to store the history in a git repository. 
We will create a git repository for each method and store different subsets of the same data. 
Each commit contains a new version of the data. Each version is a random sample containing 90% of the observations of the `airbag` data. 
Two consecutive versions of the subset will have about 90% of the observations in common.

After writing each version, we commit the file, perform garbage collection (`git gc`) on the git repository and then calculate the size of the git history (`git count-objects -v`). 

```{r git_size, eval = system.file("efficiency", "git_size.rds", package = "git2rdata") == ""}
library(git2r)
tmp_repo <- function() {
  root <- tempfile("git2rdata-efficient-git")
  dir.create(root)
  repo <- init(root)
  config(repo, user.name = "me", user.email = "me@me.com")
  return(repo)
}
commit_and_size <- function(repo, filename) {
  add(repo, filename)
  commit(repo, "test", session = TRUE)
  git_size <- system(
    sprintf("cd %s\ngit gc\ngit count-objects -v", dirname(repo$path)), 
    intern = TRUE
  )
  git_size <- git_size[grep("size-pack", git_size)]
  as.integer(gsub(".*: (.*)", "\\1", git_size))
}

repo_wt <- tmp_repo()
repo_wts <- tmp_repo()
repo_rds <- tmp_repo()
repo_wvco <- tmp_repo()
repo_wvcv <- tmp_repo()

repo_size <- replicate(
  100, {
    observed_subset <- rbinom(nrow(airbag), size = 1, prob = 0.9) == 1
    this <- airbag[
      sample(which(observed_subset)), 
      sample(ncol(airbag))
    ]
    this_sorted <- airbag[observed_subset, ]
    fn_wt <- file.path(workdir(repo_wt), "base_R.tsv")
    write.table(this, fn_wt, sep = "\t")
    fn_wts <- file.path(workdir(repo_wts), "base_R.tsv")
    write.table(this_sorted, fn_wts, sep = "\t")
    fn_rds <- file.path(workdir(repo_rds), "base_R.rds")
    saveRDS(this, fn_rds)
    fn_wvco <- write_vc(this, "airbag_optimize", repo_wvco, sorting = "X")
    fn_wvcv <- write_vc(
      this, "airbag_verbose", repo_wvcv, sorting = "X", optimize = FALSE
    )
    c(
      write.table = commit_and_size(repo_wt, fn_wt),
      write.table.sorted = commit_and_size(repo_wts, fn_wts),
      saveRDS = commit_and_size(repo_rds, fn_rds), 
      write_vc.optimized = commit_and_size(repo_wvco, fn_wvco), 
      write_vc.verbose = commit_and_size(repo_wvcv, fn_wvcv)
    )
})
```

```{r store_git_size, echo = FALSE}
if (system.file("efficiency", "git_size.rds", package = "git2rdata") == "") {
  saveRDS(repo_size, file.path("..", "inst", "efficiency", "git_size.rds"))
} else {
  repo_size <- readRDS(
    system.file("efficiency", "git_size.rds", package = "git2rdata")
  )
}
```

Each version of the data has on purpose a random order of observations and variables. This is what would happen in a worst case scenario as it would generate the largest possible diff. We also test `write.table()` with a stable ordering of the observations and variables. 

The randomised `write.table()` yields the largest git repository, converging to about `r sprintf("%.1f", repo_size["write.table", 100] / repo_size["write.table.sorted", 100])` times the size of a git repository based on the sorted `write.table()`. `saveRDS()` yields a `r sprintf("%.0f%%", 100 - 100 * repo_size["saveRDS", 100] / repo_size["write.table", 100])` reduction in repository size compared to the randomised `write.table()`, but still is `r sprintf("%.1f", repo_size["saveRDS", 100] / repo_size["write.table.sorted", 100])` times larger than the sorted `write.table()`. 
Note that the gain of storing binary files in a git repository is much smaller than the gain in individual file size because git compresses its history. 
The optimized `write_vc()` starts at `r sprintf("%.0f%%", 100 * repo_size["write_vc.optimized", 1] / repo_size["write.table.sorted", 1])` and converges toward `r sprintf("%.0f%%", 100 * repo_size["write_vc.optimized", 100] / repo_size["write.table.sorted", 100])`, the verbose version starts at `r sprintf("%.0f%%", 100 * repo_size["write_vc.verbose", 1] / repo_size["write.table.sorted", 1])` and converges towards `r sprintf("%.0f%%", 100 * repo_size["write_vc.verbose", 100] / repo_size["write.table.sorted", 100])`. 
Storage size is a lot smaller when using `write_vc()` with optimization. 
The verbose option of `write_vc()` has little the gain in storage size.
Another advantage is that `write_vc()` stores metadata.

```{r plot_git_size, echo = FALSE, fig.cap = "Size of the git history using the different storage methods."}
rs <- lapply(
  rownames(repo_size),
  function(x) {
    if (x == "saveRDS") {
      fun <- "saveRDS"
      optimized = "yes"
    } else if (x == "write_vc.optimized") {
      fun <- "write_vc"
      optimized = "yes"
    } else if (x == "write_vc.verbose") {
      fun <- "write_vc"
      optimized = "no"
    } else if (x == "write.table") {
      fun <- "write.table"
      optimized = "no"
    } else if (x == "write.table.sorted") {
      fun <- "write.table"
      optimized = "yes"
    }
    data.frame(commit = seq_along(repo_size[x, ]), size = repo_size[x, ], 
               rel_size = repo_size[x, ] / repo_size["write.table.sorted", ],
               fun = fun, optimized = optimized, stringsAsFactors = FALSE)
  }
)
rs <- do.call(rbind, rs)
rs$optimized <- factor(rs$optimized, levels = c("yes", "no"))
ggplot(rs, aes(x = commit, y = size / 2^10, colour = fun, linetype = optimized)) +
  geom_line() +
  scale_y_continuous("repo size (in MiB)") +
  scale_colour_manual("function", values = inbo_colours)
```

```{r plot_rel_git_size, echo = FALSE, fig.cap = "Relative size of the git repository when compared to write.table()."}
ggplot(rs, aes(x = commit, y = rel_size, colour = fun, linetype = optimized)) +
  geom_line() +
  scale_y_continuous("size relative to sorted write.table()", breaks = 0:10) + 
  scale_colour_manual("function", values = inbo_colours)
```

## Timings

The code below runs a microbenchmark on the four methods. A microbenchmark runs the code a hundred times and yields a distribution of timings for each expression.

### Writing Data

```{r get_file_timings, eval = system.file("efficiency", "file_timings.rds", package = "git2rdata") == ""}
library(microbenchmark)
mb <- microbenchmark(
  write.table = write.table(airbag, file.path(root, "base_R.tsv"), sep = "\t"),
  saveRDS = saveRDS(airbag, file.path(root, "base_R.rds")),
  write_vc.optim = write_vc(airbag, "airbag_optimize", root, sorting = "X"),
  write_vc.verbose = write_vc(airbag, "airbag_verbose", root, sorting = "X", 
                              optimize = FALSE)
)
mb$time <- mb$time / 1e6
```

```{r store_file_timings, echo = FALSE}
if (system.file("efficiency", "file_timings.rds", package = "git2rdata") == "") {
  saveRDS(mb, file.path("..", "inst", "efficiency", "file_timings.rds"))
} else {
  mb <- readRDS(
    system.file("efficiency", "file_timings.rds", package = "git2rdata")
  )
}
```

```{r median_write, echo = FALSE}
median_time <- aggregate(time ~ expr, data = mb, FUN = median)
write_ratio <- 100 * median_time$time / 
  median_time$time[median_time$expr == "write.table"]
names(write_ratio) <- median_time$expr
```


`write_vc()` takes `r paste(sprintf("%.0f%%", -100 + write_ratio[grep("write_vc", names(write_ratio))]), collapse = " to ")` more time than `write.table()` because it needs to prepare the metadata and sort the observations and variables. 
When overwriting existing data, `write_vc()` checks the new data against the existing metadata. 
`saveRDS()` requires `r sprintf("%.0f%%", write_ratio["saveRDS"])` of the time that `write.table()` needs.

```{r plot_file_timings, echo = FALSE, fig.cap = "Boxplot of the write timings for the different methods."}
mb$expr <- reorder(mb$expr, mb$time, FUN = median)
levels(mb$expr) <- gsub("write_vc\\.", "write_vc\n", levels(mb$expr))
ggplot(mb, aes(x = expr, y = time)) +
  geom_boxplot() +
  scale_y_continuous("Time (in milliseconds)", limits = c(0, NA)) +
  theme(axis.title.x = element_blank())
```

### Reading Data

```{r get_read_timings, eval = system.file("efficiency", "read_timings.rds", package = "git2rdata") == ""}
mb <- microbenchmark(
  read.table = read.table(file.path(root, "base_R.tsv"), header = TRUE, 
                          sep = "\t"),
  readRDS = readRDS(file.path(root, "base_R.rds")),
  read_vc.optim = read_vc("airbag_optimize", root),
  read_vc.verbose = read_vc("airbag_verbose", root)
)
mb$time <- mb$time / 1e6
```

```{r store_read_timings, echo = FALSE}
if (system.file("efficiency", "read_timings.rds", package = "git2rdata") == "") {
  saveRDS(mb, file.path("..", "inst", "efficiency", "read_timings.rds"))
} else {
  mb <- readRDS(
    system.file("efficiency", "read_timings.rds", package = "git2rdata")
  )
}
```

```{r median_read, echo = FALSE}
median_time <- aggregate(time ~ expr, data = mb, FUN = median)
read_ratio <- 100 * median_time$time / 
  median_time$time[median_time$expr == "read.table"]
names(read_ratio) <- median_time$expr
```

The timings on reading the data is another story. Reading the binary format takes about `r sprintf("%.0f%%", read_ratio["readRDS"])` of the time needed to read the standard plain text format using `read.table()`. `read_vc()` takes about `r sprintf("%.0f%%", read_ratio["read_vc.optim"])` (optimized) and `r sprintf("%.0f%%", read_ratio["read_vc.verbose"])` (verbose) of the time needed by `read.table()`, which at first seems strange because `read_vc()` calls `read.table()` to read the files and has some extra work to convert the variables to the correct data type. The main difference is that `read_vc()` knows the required data type _a priori_ and passes this info to `read.table()`. Otherwise, `read.table()` has to guess the correct data type from the file.

```{r plot_read_timings, echo = FALSE, fig.cap = "Boxplots for the read timings for the different methods."}
mb$expr <- factor(
  mb$expr, 
  levels = c("readRDS", "read.table", "read_vc.optim", "read_vc.verbose")
)
levels(mb$expr) <- gsub("read_vc\\.", "read_vc\n", levels(mb$expr))
ggplot(mb, aes(x = expr, y = time)) +
  geom_boxplot() +
  scale_y_continuous("Time (in milliseconds)", limits = c(0, NA)) +
  theme(axis.title.x = element_blank())
```
---
title: "Suggested Workflow for Storing a Variable Set of Dataframes under Version Control"
author: "Thierry Onkelinx"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Suggested Workflow for Storing a Variable Set of Dataframes under Version Control}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{git2r}
---

```{r setup, include = FALSE}
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(20120225)
```

## Introduction

This vignette describes a suggested workflow for storing a snapshot of dataframes as git2rdata objects under version control. The workflow comes in two flavours: 

  1. A single repository holding both the data and the analysis code. The single repository set-up is simple. A single reference (e.g. commit) points to both the data and the code. 
  1. One repository holding the data and a second repository holding the code. The data and the code have an independent history under a two repository set-up. Documenting the analysis requires one reference to each repository. Such a set-up is useful for repeating the same analysis (stable code) on updated data.

In this vignette we use a `git2r::repository()` object as the root. This adds git functionality to `write_vc()` and `read_vc()`, provided by the [`git2r`](https://cran.r-project.org/package=git2r) package. This allows to pull, stage, commit and push from within R.

Each commit in the data git repository describes a complete snapshot of the data at the time of the commit. 
The difference between two commits can consist of changes in existing git2rdata object (updated observations, new observations, deleted observations or updated metadata). 
Besides updating the existing git2rdata objects, we can also add new git2rdata objects or remove existing ones. 
We need to track such higher level addition and deletions as well.

We illustrate the workflow with a mock analysis on the `datasets::beaver1` and `datasets::beaver2` datasets.

## Setup

We start by initializing a git repository. `git2rdata` assumes that is already done. 
We'll use the `git2r` functions to do so. 
We start by creating a local bare repository. 
In practice we will use a remote on an external server (GitHub, Gitlab, Bitbucket, ...). 
The example below creates a local git repository with an upstream git repository. 
Any other workflow to create a similar structure is fine.

```{r initialize}
# initialize a bare git repo to be used as remote
remote <- tempfile("git2rdata-workflow-remote")
remote <- normalizePath(remote, winslash = "/")
dir.create(remote)
git2r::init(remote, bare = TRUE)

# initialize a local git repo
path <- tempfile("git2rdata-workflow")
path <- normalizePath(path, winslash = "/")
dir.create(path)
init_repo <- git2r::clone(remote, path, progress = FALSE)
git2r::config(init_repo, user.name = "me", user.email = "me@me.com")
# add an initial commit with .gitignore file
writeLines("*extra*", file.path(path, ".gitignore"))
git2r::add(init_repo, ".gitignore", force = TRUE)
git2r::commit(init_repo, message = "Initial commit")
# push initial commit to remote
git2r::push(init_repo, "origin", "refs/heads/master")
rm(init_repo)
```

## Structuring Git2rdata Objects Within a Project

`git2rdata` imposes a minimal structure. 
Both the `.tsv` and the `.yml` file need to be in the same folder. 
That's it. 
For the sake of simplicity, in this vignette we dump all git2rdata objects at the root of the repository. 

This might not be good idea for real project. 
We recommend to use at least a different directory tree for each import script. 
This directory can go into the root of a data repository.
It goes in the `data` directory in case of a data and code repository. 
Or the `inst` directory in case of an R package.

Your project might need a different directory structure. 
Feel free to choose the most relevant data structure for your project.

## Storing Dataframes _ad Hoc_ into a Git Repository

### First Commit

In the first commit we use `datasets::beaver1`.
We connect to the git repository using `repository()`.
Note that this assumes that `path` is an existing git repository.
Now we can write the dataset as a git2rdata object in the repository.
If the `root` argument of `write_vc()` is a `git_repository`, it gains two extra arguments: `stage` and `force`.
Setting `stage = TRUE`, will automatically stage the files written by `write_vc()`.

```{r store_data_1}
library(git2rdata)
repo <- repository(path)
fn <- write_vc(beaver1, "beaver", repo, sorting = "time", stage = TRUE)
```

We can use `status()` to check that `write_vc()` wrote and staged the required files.
Then we `commit()` the changes.

```{r avoid_subsecond_commit, echo = FALSE}
Sys.sleep(1.2)
```


```{r commit_data_1}
status(repo)
cm1 <- commit(repo, message = "First commit")
cat(cm1$message)
```

### Second Commit

The second commit adds `beaver2`.

```{r store_data_2}
fn <- write_vc(beaver2, "extra_beaver", repo, sorting = "time", stage = TRUE)
status(repo)
```

Notice that `extra_beaver` is not listed in the `status()`, although `write_vc()` wrote it to the repository. 
The reason is that we set a `.gitignore` which contains `"*extra*`, so git ignores any git2rdata object with a name containing "extra". 
We force git to stage it by setting `force = TRUE`.

```{r avoid_subsecond_commit2, echo = FALSE}
Sys.sleep(1.2)
```

```{r}
status(repo, ignored = TRUE)
fn <- write_vc(beaver2, "extra_beaver", repo, sorting = "time", stage = TRUE, 
               force = TRUE)
status(repo)
cm2 <- commit(repo, message = "Second commit")
```

### Third Commit

Now we decide that a single git2rdata object containing the data of both beavers is more relevant. 
We add an ID variable for each of the animals. 
This requires updating the `sorting` to avoid ties. 
And `strict = FALSE` to update the metadata. 
The "extra_beaver" git2rdata object is no longer needed so we remove it. 
We use `all = TRUE` to stage the removal of "extra_beaver" while committing the changes.

```{r avoid_subsecond_commit3, echo = FALSE}
Sys.sleep(1.2)
```

```{r store_data_3}
beaver1$beaver <- 1
beaver2$beaver <- 2
beaver <- rbind(beaver1, beaver2)
fn <- write_vc(beaver, "beaver", repo, sorting = c("beaver", "time"), 
               strict = FALSE, stage = TRUE)
file.remove(list.files(path, "extra", full.names = TRUE))
status(repo)
cm3 <- commit(repo, message = "Third commit", all = TRUE)
status(repo)
```

## Scripted Workflow for Storing Dataframes

We strongly recommend to add git2rdata object through an import script instead of adding them [_ad hoc_](#storing-dataframes-ad-hoc-into-a-git-repository). Store this script in the (analysis) repository. It documents the creation of the git2rdata objects. Rerun this script whenever updated data becomes available. 

Old versions of the import script and the associated git2rdata remain available through the version control history. Remove obsolete git2rdata objects from the import script. This keeps both the import script and the working directory tidy and minimal.

Basically, the import script should create all git2rdata objects within a given directory tree. 
This gives the advantage that we start the import script by clearing any existing git2rdata object in this directory. 
If the import script no longer creates a git2rdata object, it gets removed without the need to track what git2rdata objects existed in the previous version.

The brute force method of removing all files or all `.tsv` / `.yml` pairs is not a good idea. This removes the existing metadata which we need for efficient storage (see `vignette("efficiency", package = "git2rdata")`). A better solution is to use `rm_data()` on the directory at the start of the import script. This removes all `.tsv` files which have valid metadata. The existing metadata remains untouched at this point.

Then write all git2rdata objects and stage them. Unchanged objects will not lead to a diff, even if we first deleted and then recreated them. The script won't recreate the `.tsv` file of obsolete git2rdata objects. Use `prune_meta()` to remove any leftover metadata files.

Commit and push the changes at the end of the script.

Below is an example script recreating the "beaver" git2rdata object from the [third commit](#third-commit).

```{r eval = FALSE}
# load package
library(git2rdata)
# step 1: setup the repository and data path
repo <- repository(".")
data_path <- "data/beaver"
# step 1b: sync the repository with the remote
pull(repo = repo)
# step 2: remove all existing data files
rm_data(root = repo, path = data_path, stage = TRUE)

# step 3: write all relevant git2rdata objects to the data path
beaver1$beaver <- 1
beaver2$beaver <- 2
body_temp <- rbind(beaver1, beaver2)
fn <- write_vc(x = body_temp, file = file.path(data_path, "body_temperature"), 
               root = repo, sorting = c("beaver", "time"), stage = TRUE)

# step 4: remove any dangling metadata files
prune_meta(root = repo, path = data_path, stage = TRUE)

# step 5: commit the changes
cm <- commit(repo = repo, message = "import")
# step 5b: sync the repository with the remote
push(repo = repo)
```

## R Package Workflow for Storing Dataframes

We recommend a two repository set-up in case of recurring analyses. 
These are relative stable analyses which have to run with some frequency on updated data (e.g. once a month). 
That makes it worthwhile to convert the analyses into an R package. 
Split long scripts into a set of shorter functions which are much easier to document and maintain. 
An R package offers lots of [functionality](http://r-pkgs.had.co.nz/check.html) out of the box to check the quality of your code.

The example below converts the import script above into a function. 
We illustrate how you can use Roxygen2 (see `vignette("roxygen2", package = "roxygen2")`) tags to document the function and to list its dependencies.
Note that we added `session = TRUE` to `commit()`. 
This will append the `sessionInfo()` at the time of the commit to the commit message. 
Thus documenting all loaded R packages and their version. 
This documents to code used to create the git2rdata object since your analysis code resides in a dedicated package with its own version number. 
We strongly recommend to run the import from a fresh R session. 
Then the `sessionInfo()` at commit time contains those packages with are strictly required for the import.
Consider running the import from the command line. e.g. `Rscript -e 'mypackage::import_body_temp("path/to/root")'`.

```{r eval = FALSE}
#' Import the beaver body temperature data
#' @param path the root of the git repository
#' @importFrom git2rdata repository pull rm_data write_vc prune_meta commit push
#' @export
import_body_temp <- function(path) {
  # step 1: setup the repository and data path
  repo <- repository(path)
  data_path <- "data/beaver"
  # step 1b: sync the repository with the remote
  pull(repo = repo)
  # step 2: remove all existing data files
  rm_data(root = repo, path = data_path, stage = TRUE)
  
  # step 3: write all relevant git2rdata objects to the data path
  beaver1$beaver <- 1
  beaver2$beaver <- 2
  body_temp <- rbind(beaver1, beaver2)
  fn <- write_vc(x = body_temp, file = file.path(data_path, "body_temperature"), 
                 root = repo, sorting = c("beaver", "time"), stage = TRUE)
  
  # step 4: remove any dangling metadata files
  prune_meta(root = repo, path = data_path, stage = TRUE)
  
  # step 5: commit the changes
  cm <- commit(repo = repo, message = "import", session = TRUE)
  # step 5b: sync the repository with the remote
  push(repo = repo)
}
```

## Analysis Workflow with Reproducible Data

The example below is a small trivial example of a standardized analysis in which documents the source of the data by describing the name of the data, the repository URL and the commit. 
We can use this information when reporting the results. This makes the data underlying the results traceable.

```{r standardized_analysis}
analysis <- function(ds_name, repo) {
  ds <- read_vc(ds_name, repo)
  list(
    dataset = ds_name,
    repository = git2r::remote_url(repo),
    commit = recent_commit(ds_name, repo, data = TRUE),
    model = lm(temp ~ activ, data = ds)
  )
}
report <- function(x) {
  knitr::kable(
    coef(summary(x$model)),
    caption = sprintf("**dataset:** %s  \n**commit:** %s  \n**repository:** %s", 
                      x$dataset, x$commit$commit, x$repository)
  )
}
```

In this case we can run every analysis by looping over the list of datasets in the repository.

```{r run_current_analyses, results = "asis"}
repo <- repository(path)
current <- lapply(list_data(repo), analysis, repo = repo)
names(current) <- list_data(repo)
result <- lapply(current, report)
junk <- lapply(result, print)
```

The example below does the same thing for the first and second commit. 

```{r run_previous_analyses, results = "asis"}
# checkout first commit
git2r::checkout(cm1)
# do analysis
previous <- lapply(list_data(repo), analysis, repo = repo)
names(previous) <- list_data(repo)
result <- lapply(previous, report)
junk <- lapply(result, print)
# checkout second commit
git2r::checkout(cm2)
# do analysis
previous <- lapply(list_data(repo), analysis, repo = repo)
names(previous) <- list_data(repo)
result <- lapply(previous, report)
junk <- lapply(result, print)
```

If you inspect the reported results you'll notice that all the output (coefficients and commit hash) for "beaver" object is identical for the first and second commit. 
This makes sense since the "beaver" object didn't change during the second commit. 
The output for the current (third) commit is different because the dataset changed.

### Long running analysis

Imagine the case where an individual analysis takes a while to run. 
We store the most recent version of each analysis and add the information from `recent_commit()`. 
When preparing the analysis, you can run `recent_commit()` again on the dataset and compare the commit hash with that one of the available analysis. 
If the commit hashes match, then the data hasn't changed. 
Then there is no need to rerun the analysis^[assuming the code for running the analysis didn't change.], saving valuable computing resources and time. 
---
title: "Getting Started Storing Dataframes as Plain Text"
author: "Thierry Onkelinx"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started Storing Dataframes as Plain Text}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)
options(width = 83)
```

## Introduction

This vignette motivates why we wrote `git2rdata` and illustrates how you can use it to store dataframes as plain text files.

### Maintaining Variable Classes

R has different options to store dataframes as plain text files from R. 
Base R has `write.table()` and its companions like `write.csv()`. 
Some other options are `data.table::fwrite()`, `readr::write_delim()`, `readr::write_csv()` and `readr::write_tsv()`. 
Each of them writes a dataframe as a plain text file by converting all variables into characters. 
After reading the file, they revert this conversion. 
The distinction between `character` and `factor` gets lost in translation.
`read.table()` converts by default all strings to factors, `readr::read_csv()` keeps by default all strings as character. 
These functions cannot recover the factor levels.
These functions determine factor levels based on the observed levels in the plain text file. 
Hence factor levels without observations will disappear. 
The order of the factor levels is also determined by the available levels in the plain text file, which can be different from the original order.

The `write_vc()` and `read_vc()` functions from `git2rdata` keep track of the class of each variable and, in case of a factor, also of the factor levels and their order. Hence this function pair preserves the information content of the dataframe. The `vc` suffix stands for **v**ersion **c**ontrol as these functions use their full capacity in combination with a version control system.

## Efficiency Relative to Storage and Time

### Optimizing File Storage

Plain text files require more disk space than binary files. 
This is the price we have to pay for a readable file format. 
The default option of `write_vc()` is to create file as compact as possible. 
Since we use a tab delimited file format, we can omit quotes around character variables. 
This saves 2 bytes per row for each character variable. 
`write_vc` add quotes automatically in the exceptional cases when we needed them, e.g. to store a string that contains tab or newline characters.
We don't add quotes to row-variable combinations where we don't need them.

Since we store the class of each variable, we can further reduce the file size by following rules:

- Store a `logical` as 0 (FALSE), 1 (TRUE) or NA to the data.
- Store a `factor` as its indices in the data. 
Store the index, labels of levels and their order in the metadata.
- Store a `POSIXct` as a numeric to the data. 
Store the class and the origin in the metadata. 
Store and return timestamps as UTC.
- Store a `Date` as an integer to the data. 
Store the class and the origin in the metadata.

Storing the factors, POSIXct and Date as their index, makes them less user readable. The user can turn off this optimization when user readability is more important than file size.

### Optimized for Version Control

Another main goal of `git2rdata` is to optimise the storage of the plain text files under version control. `write_vc()` and `read_vc()` has methods for interacting with [git](https://git-scm.com/) repositories using the `git2r` framework. Users who want to use git without `git2r` or use a different version control system (e.g. [Subversion](https://subversion.apache.org/), [Mercurial](https://www.mercurial-scm.org/)), still can use `git2rdata` to write the files to disk and uses their preferred workflow on version control. 

Hence, `write_vc()` will always perform checks to look for changes which potentially lead to large diffs. More details on this in `vignette("version_control", package = "git2rdata")`. Some problems will always yield a warning. Other problems will yield an error by default. The user can turn these errors into warnings by setting the `strict = FALSE` argument.

As this vignette ignores the part on version control, we will always use `write_vc(strict = FALSE)` and hide the warnings to improve the readability.

## Basic Usage

Let's start by setting up the environment. We need a directory to store the data and a dataframe to store.

```{r}
# Create a directory in tempdir
path <- tempfile(pattern = "git2r-")
dir.create(path)
# Create dummy data
set.seed(20190222)
x <- data.frame(
  x = sample(LETTERS),
  y = factor(
    sample(c("a", "b", NA), 26, replace = TRUE),
    levels = c("a", "b", "c")
  ),
  z = c(NA, 1:25),
  abc = c(rnorm(25), NA),
  def = sample(c(TRUE, FALSE, NA), 26, replace = TRUE),
  timestamp = seq(
    as.POSIXct("2018-01-01"),
    as.POSIXct("2019-01-01"),
    length = 26
  ),
  stringsAsFactors = FALSE
)
str(x)
```

## Storing Optimized

Use `write_vc()` to store the dataframe. 
The `root` argument refers to the base directory where we store the data. 
The `file` argument becomes the base name of the files. 
The data file gets a `.tsv` extension, the metadata file a `.yml` extension. 
`file` can include a relative path starting from `root`.

```{r first_write}
library(git2rdata)
write_vc(x = x, file = "first_test", root = path, strict = FALSE)
```

`write_vc()` returns a vector of relative paths to the raw data and metadata files. 
The names of this vector contains the hashes of these files. 
We can have a look at both files. 
We'll display the first 10 rows of the raw data. 
Notice that the YAML format of the metadata has the benefit of being both human and machine readable. 

```{r manual_data}
print_file <- function(file, root, n = -1) {
  fn <- file.path(root, file)
  data <- readLines(fn, n = n)
  cat(data, sep = "\n")
}
print_file("first_test.tsv", path, 10)
print_file("first_test.yml", path)
```


## Storing Verbose

Adding `optimize = FALSE` to `write_vc()` will keep the raw data in a human readable format. The metadata file is slightly different. The most obvious is the `optimize: no` tag and the different hash. Another difference is the metadata for POSIXct and Date classes. They will no longer have an origin tag but a format tag.

```{r write_verbose}
write_vc(x = x, file = "verbose", root = path, optimize = FALSE, strict = FALSE)
```

```{r manual_verbose_data}
print_file("verbose.tsv", path, 10)
print_file("verbose.yml", path)
```

## Efficiency Relative to File Storage

Storing dataframes optimized or verbose has an impact on the required file size. 
The [efficiency](efficiency.html#on-a-file-system) vignette give a comparison.

## Reading Data

You retrieve the data with `read_vc()`. 
This function will reinstate the variables to their original state.

```{r first_read}
y <- read_vc(file = "first_test", root = path)
all.equal(x, y, check.attributes = FALSE)
y2 <- read_vc(file = "verbose", root = path)
all.equal(x, y2, check.attributes = FALSE)
```

`read_vc()` requires the meta data. 
It cannot handle dataframe not stored by `write_vc()`.

## Missing Values

`write_vc()` has an `na` argument which specifies the string which to use for missing values. 
Because we avoid using quotes, this string must be different from any character value in the data. 
This includes factor labels with verbose data storage. 
`write_vc()` checks this and will always return an error, even with `strict = FALSE`.

```{r echo = FALSE, results = "hide"}
stopifnot("X" %in% x$x, "b" %in% x$y)
```

```{r na_string, error = TRUE}
write_vc(x, "custom_na", path, strict = FALSE, na = "X", optimize = FALSE)
write_vc(x, "custom_na", path, strict = FALSE, na = "b", optimize = FALSE)
write_vc(x, "custom_na", path, strict = FALSE, na = "X")
write_vc(x, "custom_na", path, strict = FALSE, na = "b")
```

Please note that `write_vc()` uses the same  NA string for the entire dataset, thus for every variable. 

```{r manual_na_data}
print_file("custom_na.tsv", path, 10)
print_file("custom_na.yml", path, 4)
```

The default string for missing values is `"NA"`. We recommend to keep this default, as long as the dataset permits it. A first good alternative is an empty string (`""`). If that won't work either, you'll have to use your imagination. Try to keep it short, clear and robust^[robust in the sense that you won't need to change it later].

```{r empty_na}
write_vc(x, "custom_na", path, strict = FALSE, na = "")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.R
\name{prune_meta}
\alias{prune_meta}
\alias{prune_meta.git_repository}
\title{Prune Metadata Files}
\usage{
prune_meta(root = ".", path = NULL, recursive = TRUE, ...)

\method{prune_meta}{git_repository}(root, path = NULL, recursive = TRUE, ..., stage = FALSE)
}
\arguments{
\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{path}{the directory in which to clean all the data files. The directory
is relative to \code{root}.}

\item{recursive}{remove files in subdirectories too.}

\item{...}{parameters used in some methods}

\item{stage}{stage the changes after removing the files. Defaults to \code{FALSE}.}
}
\value{
returns invisibly a vector of removed files names. The paths are
relative to \code{root}.
}
\description{
Removes all \strong{valid} metadata (\code{.yml} files) from the \code{path} when they don't
have accompanying data (\code{.tsv} file). \strong{Invalid} metadata triggers a warning
without removing the metadata file.

Use this function with caution since it will remove all valid metadata files
without asking for confirmation. We strongly recommend to use this
function on files under version control. See
\code{vignette("workflow", package = "git2rdata")} for some examples on how to use
this.
}
\examples{
## on file system

# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# store a dataframe as git2rdata object. Capture the result to minimise
# screen output
junk <- write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
# write a standard tab separate file (non git2rdata object)
write.table(iris, file = file.path(root, "standard.tsv"), sep = "\t")
# write a YAML file
yml <- list(
  authors = list(
   "Research Institute for Nature and Forest" = list(
       href = "https://www.inbo.be/en")))
yaml::write_yaml(yml, file = file.path(root, "_pkgdown.yml"))

# list the git2rdata objects
list_data(root)
# list the files
list.files(root, recursive = TRUE)

# remove all .tsv files from valid git2rdata objects
rm_data(root, path = ".")
# check the removal of the .tsv file
list.files(root, recursive = TRUE)
list_data(root)

# remove dangling git2rdata metadata files
prune_meta(root, path = ".")
# check the removal of the metadata
list.files(root, recursive = TRUE)
list_data(root)


## on git repo

# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# store a dataframe
write_vc(iris[1:6, ], "iris", repo, sorting = "Sepal.Length", stage = TRUE)
# check that the dataframe is stored
status(repo)
list_data(repo)

# commit the current version and check the git repo
commit(repo, "add iris data", session = TRUE)
status(repo)

# remove the data files from the repo
rm_data(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# remove dangling metadata
prune_meta(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# clean up
junk <- file.remove(
  list.files(root, full.names = TRUE, recursive = TRUE), root)
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{read_vc}()},
\code{\link{relabel}()},
\code{\link{rename_variable}()},
\code{\link{rm_data}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relabel.R
\name{relabel}
\alias{relabel}
\title{Relabel Factor Levels by Updating the Metadata}
\usage{
relabel(file, root = ".", change)
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{change}{either a \code{list} or a \code{data.frame}. In case of a \code{list} is a
named \code{list} with named \code{vectors}. The names of list elements must match the
names of the variables. The names of the vector elements must match the
existing factor labels. The values represent the new factor labels. In case
of a \code{data.frame} it needs to have the variables \code{factor} (name of the
factor), \code{old} (the old) factor label and \code{new} (the new factor label).
\code{relabel()} ignores all other columns.}
}
\value{
invisible \code{NULL}.
}
\description{
Imagine the situation where we have a dataframe with a factor variable and we
have stored it with \code{write_vc(optimize = TRUE)}. The raw data file contains
the factor indices and the metadata contains the link between the factor
index and the corresponding label. See
\code{vignette("version_control", package = "git2rdata")}. In such a case,
relabelling a factor can be fast and lightweight by updating the metadata.
}
\examples{

# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a dataframe and store it as an optimized git2rdata object.
# Note that write_vc() uses optimization by default.
# Stage and commit the git2rdata object.
ds <- data.frame(
  a = c("a1", "a2"),
  b = c("b2", "b1"),
  stringsAsFactors = TRUE
)
junk <- write_vc(ds, "relabel", repo, sorting = "b", stage = TRUE)
cm <- commit(repo, "initial commit")
# check that the workspace is clean
status(repo)

# Define new labels as a list and apply them to the git2rdata object.
new_labels <- list(
  a = list(a2 = "a3")
)
relabel("relabel", repo, new_labels)
# check the changes
read_vc("relabel", repo)
# relabel() changed the metadata, not the raw data
status(repo)
git2r::add(repo, "relabel.*")
cm <- commit(repo, "relabel using a list")

# Define new labels as a dataframe and apply them to the git2rdata object
change <- data.frame(
  factor = c("a", "a", "b"),
  old = c("a3", "a1", "b2"),
  new = c("c2", "c1", "b3"),
  stringsAsFactors = TRUE
)
relabel("relabel", repo, change)
# check the changes
read_vc("relabel", repo)
# relabel() changed the metadata, not the raw data
status(repo)

# clean up
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{prune_meta}()},
\code{\link{read_vc}()},
\code{\link{rename_variable}()},
\code{\link{rm_data}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_vc.R
\name{read_vc}
\alias{read_vc}
\title{Read a Git2rdata Object from Disk}
\usage{
read_vc(file, root = ".")
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}
}
\value{
The \code{data.frame} with the file names and hashes as attributes.
}
\description{
\code{read_vc()} handles git2rdata objects stored by \code{write_vc()}. It reads and
verifies the metadata file (\code{.yml}). Then it reads and verifies the raw data.
The last step is back-transforming any transformation done by \code{meta()} to
return the \code{data.frame} as stored by \code{write_vc()}.

\code{read_vc()} is an S3 generic on \code{root} which currently handles \code{"character"}
(a path) and \code{"git-repository"} (from \code{git2r}). S3 methods for other version
control system could be added.
}
\examples{
## on file system

# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# write a dataframe to the directory
write_vc(iris[1:6, ], file = "iris", root = root, sorting = "Sepal.Length")
# check that a data file (.tsv) and a metadata file (.yml) exist.
list.files(root, recursive = TRUE)
# read the git2rdata object from the directory
read_vc("iris", root)

# store a new version with different observations but the same metadata
write_vc(iris[1:5, ], "iris", root)
list.files(root, recursive = TRUE)
# Removing a column requires version requires new metadata.
# Add strict = FALSE to override the existing metadata.
write_vc(
  iris[1:6, -2], "iris", root, sorting = "Sepal.Length", strict = FALSE
)
list.files(root, recursive = TRUE)
# storing the orignal version again requires another update of the metadata
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Width", strict = FALSE)
list.files(root, recursive = TRUE)
# optimize = FALSE stores the data more verbose. This requires larger files.
write_vc(
  iris[1:6, ], "iris2", root, sorting = "Sepal.Width", optimize = FALSE
)
list.files(root, recursive = TRUE)



## on git repo using a git2r::git-repository

# initialise a git repo using the git2r package
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# store a dataframe in git repo.
write_vc(iris[1:6, ], file = "iris", root = repo, sorting = "Sepal.Length")
# This git2rdata object is not staged by default.
status(repo)
# read a dataframe from a git repo
read_vc("iris", repo)

# store a new version in the git repo and stage it in one go
write_vc(iris[1:5, ], "iris", repo, stage = TRUE)
status(repo)

# store a verbose version in a different gir2data object
write_vc(
  iris[1:6, ], "iris2", repo, sorting = "Sepal.Width", optimize = FALSE
)
status(repo)

# clean up
junk <- file.remove(
  list.files(root, full.names = TRUE, recursive = TRUE), root)
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{prune_meta}()},
\code{\link{relabel}()},
\code{\link{rename_variable}()},
\code{\link{rm_data}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexport.R
\name{repository}
\alias{repository}
\title{Re-exported Function From \code{git2r}}
\description{
See \code{\link[git2r]{repository}} in \code{git2r}.
}
\seealso{
Other version_control: 
\code{\link{commit}()},
\code{\link{pull}()},
\code{\link{push}()},
\code{\link{recent_commit}()},
\code{\link{status}()}
}
\concept{version_control}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_variable.R
\name{rename_variable}
\alias{rename_variable}
\alias{rename_variable.character}
\alias{rename_variable.default}
\alias{rename_variable.git_repository}
\title{Rename a Variable}
\usage{
rename_variable(file, change, root = ".", ...)

\method{rename_variable}{character}(file, change, root = ".", ...)

\method{rename_variable}{default}(file, change, root, ...)

\method{rename_variable}{git_repository}(file, change, root, ..., stage = FALSE, force = FALSE)
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{change}{A named vector with the old names as values and the new names
as names.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{...}{parameters used in some methods}

\item{stage}{Logical value indicating whether to stage the changes after
writing the data. Defaults to \code{FALSE}.}

\item{force}{Add ignored files. Default is FALSE.}
}
\value{
invisible \code{NULL}.
}
\description{
The raw data file contains a header with the variable names.
The metadata list the variable names and their type.
Changing a variable name and overwriting the \code{git2rdata} object with result
in an error.
Because it will look like removing an existing variable and adding a new one.
Overwriting the object with \code{strict = FALSE} potentially changes the order of
the variables, leading to a large diff.
}
\details{
This function solves this by only updating the raw data header and the
metadata.
}
\examples{

# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# Create a dataframe and store it as an optimized git2rdata object.
# Note that write_vc() uses optimization by default.
# Stage and commit the git2rdata object.
ds <- data.frame(
  a = c("a1", "a2"),
  b = c("b2", "b1"),
  stringsAsFactors = TRUE
)
junk <- write_vc(ds, "rename", repo, sorting = "b", stage = TRUE)
cm <- commit(repo, "initial commit")
# check that the workspace is clean
status(repo)

# Define change.
change <- c(new_name = "a")
rename_variable(file = "rename", change = change, root = repo)
# check the changes
read_vc("rename", repo)
status(repo)

# clean up
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{prune_meta}()},
\code{\link{read_vc}()},
\code{\link{relabel}()},
\code{\link{rm_data}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.R
\name{rm_data}
\alias{rm_data}
\alias{rm_data.git_repository}
\title{Remove Data Files From Git2rdata Objects}
\usage{
rm_data(root = ".", path = NULL, recursive = TRUE, ...)

\method{rm_data}{git_repository}(
  root,
  path = NULL,
  recursive = TRUE,
  ...,
  stage = FALSE,
  type = c("unmodified", "modified", "ignored", "all")
)
}
\arguments{
\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{path}{the directory in which to clean all the data files. The directory
is relative to \code{root}.}

\item{recursive}{remove files in subdirectories too.}

\item{...}{parameters used in some methods}

\item{stage}{stage the changes after removing the files. Defaults to FALSE.}

\item{type}{Defines the classes of files to remove. \code{unmodified} are files in
the git history and unchanged since the last commit. \code{modified} are files in
the git history and changed since the last commit. \code{ignored} refers to file
listed in a \code{.gitignore} file. Selecting \code{modified} will remove both
\code{unmodified} and \code{modified} data files. Selecting \code{ìgnored} will remove
\code{unmodified}, \code{modified} and \code{ignored} data files. \code{all} refers to all
visible data files, including \code{untracked} files.}
}
\value{
returns invisibly a vector of removed files names. The paths are
relative to \code{root}.
}
\description{
Remove the data (\code{.tsv}) file from all valid git2rdata objects at the \code{path}.
The metadata remains untouched. A warning lists any git2rdata object with
\strong{invalid} metadata. The function keeps any \code{.tsv} file with
invalid metadata or from non-git2rdata objects.

Use this function with caution since it will remove all valid data files
without asking for confirmation. We strongly recommend to use this
function on files under version control. See
\code{vignette("workflow", package = "git2rdata")} for some examples on how to use
this.
}
\examples{
## on file system

# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# store a dataframe as git2rdata object. Capture the result to minimise
# screen output
junk <- write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
# write a standard tab separate file (non git2rdata object)
write.table(iris, file = file.path(root, "standard.tsv"), sep = "\t")
# write a YAML file
yml <- list(
  authors = list(
   "Research Institute for Nature and Forest" = list(
       href = "https://www.inbo.be/en")))
yaml::write_yaml(yml, file = file.path(root, "_pkgdown.yml"))

# list the git2rdata objects
list_data(root)
# list the files
list.files(root, recursive = TRUE)

# remove all .tsv files from valid git2rdata objects
rm_data(root, path = ".")
# check the removal of the .tsv file
list.files(root, recursive = TRUE)
list_data(root)

# remove dangling git2rdata metadata files
prune_meta(root, path = ".")
# check the removal of the metadata
list.files(root, recursive = TRUE)
list_data(root)


## on git repo

# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# store a dataframe
write_vc(iris[1:6, ], "iris", repo, sorting = "Sepal.Length", stage = TRUE)
# check that the dataframe is stored
status(repo)
list_data(repo)

# commit the current version and check the git repo
commit(repo, "add iris data", session = TRUE)
status(repo)

# remove the data files from the repo
rm_data(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# remove dangling metadata
prune_meta(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# clean up
junk <- file.remove(
  list.files(root, full.names = TRUE, recursive = TRUE), root)
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{prune_meta}()},
\code{\link{read_vc}()},
\code{\link{relabel}()},
\code{\link{rename_variable}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.R
\docType{methods}
\name{meta}
\alias{meta}
\alias{meta.character}
\alias{meta.factor}
\alias{meta.logical}
\alias{meta.POSIXct}
\alias{meta.Date}
\alias{meta.data.frame}
\title{Optimize an Object for Storage as Plain Text and Add Metadata}
\usage{
meta(x, ...)

\method{meta}{character}(x, na = "NA", ...)

\method{meta}{factor}(x, optimize = TRUE, na = "NA", index, strict = TRUE, ...)

\method{meta}{logical}(x, optimize = TRUE, ...)

\method{meta}{POSIXct}(x, optimize = TRUE, ...)

\method{meta}{Date}(x, optimize = TRUE, ...)

\method{meta}{data.frame}(
  x,
  optimize = TRUE,
  na = "NA",
  sorting,
  strict = TRUE,
  split_by = character(0),
  ...
)
}
\arguments{
\item{x}{the vector.}

\item{...}{further arguments to the methods.}

\item{na}{the string to use for missing values in the data.}

\item{optimize}{If \code{TRUE}, recode the data to get smaller text files.
If \code{FALSE}, \code{meta()} converts the data to character.
Defaults to \code{TRUE}.}

\item{index}{An optional named vector with existing factor indices.
The names must match the existing factor levels.
Unmatched levels from \code{x} will get new indices.}

\item{strict}{What to do when the metadata changes. \code{strict = FALSE}
overwrites the data and the metadata with a warning listing the changes,
\code{strict = TRUE} returns an error and leaves the data and metadata as is.
Defaults to \code{TRUE}.}

\item{sorting}{an optional vector of column names defining which columns to
use for sorting \code{x} and in what order to use them.
The default empty \code{sorting} yields a warning.
Add \code{sorting} to avoid this warning.
Strongly recommended in combination with version control.
See \code{vignette("efficiency", package = "git2rdata")} for an illustration of
the importance of sorting.}

\item{split_by}{An optional vector of variables name to split the text files.
This creates a separate file for every combination.
We prepend these variables to the vector of \code{sorting} variables.}
}
\value{
the optimized vector \code{x} with \code{meta} attribute.
}
\description{
Prepares a vector for storage. When relevant, \code{meta()} optimizes the object
for storage by changing the format to one which needs less characters. The
metadata stored in the \code{meta} attribute, contains all required information to
back-transform the optimized format into the original format.

In case of a data.frame, \code{meta()} applies itself to each of the columns. The
\code{meta} attribute becomes a named list containing the metadata for each column
plus an additional \code{..generic} element. \code{..generic} is a reserved name for
the metadata and not allowed as column name in a \code{data.frame}.

\code{\link{write_vc}} uses this function to prepare a dataframe for storage.
Existing metadata is passed through the optional \code{old} argument. This
argument intended for internal use.
}
\note{
The default order of factor levels depends on the current locale.
See \code{\link{Comparison}} for more details on that.
The same code on a different locale might result in a different sorting.
\code{meta()} ignores, with a warning, any change in the order of factor levels.
Add \code{strict = FALSE} to enforce the new order of factor levels.
}
\examples{
meta(c(NA, "'NA'", '"NA"', "abc\tdef", "abc\ndef"))
meta(1:3)
meta(seq(1, 3, length = 4))
meta(factor(c("b", NA, "NA"), levels = c("NA", "b", "c")))
meta(factor(c("b", NA, "a"), levels = c("a", "b", "c")), optimize = FALSE)
meta(factor(c("b", NA, "a"), levels = c("a", "b", "c"), ordered = TRUE))
meta(
  factor(c("b", NA, "a"), levels = c("a", "b", "c"), ordered = TRUE),
  optimize = FALSE
)
meta(c(FALSE, NA, TRUE))
meta(c(FALSE, NA, TRUE), optimize = FALSE)
meta(complex(real = c(1, NA, 2), imaginary = c(3, NA, -1)))
meta(as.POSIXct("2019-02-01 10:59:59", tz = "CET"))
meta(as.POSIXct("2019-02-01 10:59:59", tz = "CET"), optimize = FALSE)
meta(as.Date("2019-02-01"))
meta(as.Date("2019-02-01"), optimize = FALSE)
}
\seealso{
Other internal: 
\code{\link{is_git2rdata}()},
\code{\link{is_git2rmeta}()},
\code{\link{upgrade_data}()}
}
\concept{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexport.R
\name{status}
\alias{status}
\title{Re-exported Function From \code{git2r}}
\description{
See \code{\link[git2r]{status}} in \code{git2r}.
}
\seealso{
Other version_control: 
\code{\link{commit}()},
\code{\link{pull}()},
\code{\link{push}()},
\code{\link{recent_commit}()},
\code{\link{repository}()}
}
\concept{version_control}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_git2rdata.R
\name{is_git2rdata}
\alias{is_git2rdata}
\title{Check Whether a Git2rdata Object is Valid.}
\usage{
is_git2rdata(file, root = ".", message = c("none", "warning", "error"))
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{message}{a single value indicating the type of messages on top of the
logical value. \code{"none"}: no messages, \code{"warning"}: issue a warning in case of
an invalid metadata file. \code{"error"}: an invalid metadata file results in an
error. Defaults to \code{"none"}.}
}
\value{
A logical value. \code{TRUE} in case of a valid git2rdata object.
Otherwise \code{FALSE}.
}
\description{
A valid git2rdata object has valid metadata.
}
\examples{
# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# store a file
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
# check the stored file
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# Remove the metadata from the existing git2rdata object. Then it stops
# being a git2rdata object.
junk <- file.remove(file.path(root, "iris.yml"))
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# recreate the file and remove the data and keep the metadata. It stops being
# a git2rdata object, but the metadata remains valid.
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
junk <- file.remove(file.path(root, "iris.tsv"))
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# clean up
junk <- file.remove(list.files(root, full.names = TRUE), root)
}
\seealso{
Other internal: 
\code{\link{is_git2rmeta}()},
\code{\link{meta}()},
\code{\link{upgrade_data}()}
}
\concept{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexport.R
\name{pull}
\alias{pull}
\title{Re-exported Function From \code{git2r}}
\description{
See \code{\link[git2r]{pull}} in \code{git2r}.
}
\seealso{
Other version_control: 
\code{\link{commit}()},
\code{\link{push}()},
\code{\link{recent_commit}()},
\code{\link{repository}()},
\code{\link{status}()}
}
\concept{version_control}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexport.R
\name{commit}
\alias{commit}
\title{Re-exported Function From \code{git2r}}
\description{
See \code{\link[git2r]{commit}} in \code{git2r}.
}
\seealso{
Other version_control: 
\code{\link{pull}()},
\code{\link{push}()},
\code{\link{recent_commit}()},
\code{\link{repository}()},
\code{\link{status}()}
}
\concept{version_control}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_git2rmeta.R
\name{is_git2rmeta}
\alias{is_git2rmeta}
\title{Check Whether a Git2rdata Object Has Valid Metadata.}
\usage{
is_git2rmeta(file, root = ".", message = c("none", "warning", "error"))
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{message}{a single value indicating the type of messages on top of the
logical value. \code{"none"}: no messages, \code{"warning"}: issue a warning in case of
an invalid metadata file. \code{"error"}: an invalid metadata file results in an
error. Defaults to \code{"none"}.}
}
\value{
A logical value. \code{TRUE} in case of a valid metadata file. Otherwise
\code{FALSE}.
}
\description{
Valid metadata is a file with \code{.yml} extension. It has a top level item
\code{..generic}. This item contains \code{git2rdata} (the version number), \code{hash} (a
hash on the metadata) and \code{data_hash} (a hash on the data file). The version
number must be the current version.
}
\examples{
# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# store a file
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
# check the stored file
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# Remove the metadata from the existing git2rdata object. Then it stops
# being a git2rdata object.
junk <- file.remove(file.path(root, "iris.yml"))
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# recreate the file and remove the data and keep the metadata. It stops being
# a git2rdata object, but the metadata remains valid.
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
junk <- file.remove(file.path(root, "iris.tsv"))
is_git2rmeta("iris", root)
is_git2rdata("iris", root)

# clean up
junk <- file.remove(list.files(root, full.names = TRUE), root)
}
\seealso{
Other internal: 
\code{\link{is_git2rdata}()},
\code{\link{meta}()},
\code{\link{upgrade_data}()}
}
\concept{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_vc.R
\name{write_vc}
\alias{write_vc}
\alias{write_vc.character}
\alias{write_vc.git_repository}
\title{Store a Data.Frame as a Git2rdata Object on Disk}
\usage{
write_vc(
  x,
  file,
  root = ".",
  sorting,
  strict = TRUE,
  optimize = TRUE,
  na = "NA",
  ...,
  split_by
)

\method{write_vc}{character}(
  x,
  file,
  root = ".",
  sorting,
  strict = TRUE,
  optimize = TRUE,
  na = "NA",
  ...,
  split_by = character(0)
)

\method{write_vc}{git_repository}(
  x,
  file,
  root,
  sorting,
  strict = TRUE,
  optimize = TRUE,
  na = "NA",
  ...,
  stage = FALSE,
  force = FALSE
)
}
\arguments{
\item{x}{the \code{data.frame}.}

\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{sorting}{an optional vector of column names defining which columns to
use for sorting \code{x} and in what order to use them.
The default empty \code{sorting} yields a warning.
Add \code{sorting} to avoid this warning.
Strongly recommended in combination with version control.
See \code{vignette("efficiency", package = "git2rdata")} for an illustration of
the importance of sorting.}

\item{strict}{What to do when the metadata changes. \code{strict = FALSE}
overwrites the data and the metadata with a warning listing the changes,
\code{strict = TRUE} returns an error and leaves the data and metadata as is.
Defaults to \code{TRUE}.}

\item{optimize}{If \code{TRUE}, recode the data to get smaller text files.
If \code{FALSE}, \code{meta()} converts the data to character.
Defaults to \code{TRUE}.}

\item{na}{the string to use for missing values in the data.}

\item{...}{parameters used in some methods}

\item{split_by}{An optional vector of variables name to split the text files.
This creates a separate file for every combination.
We prepend these variables to the vector of \code{sorting} variables.}

\item{stage}{Logical value indicating whether to stage the changes after
writing the data. Defaults to \code{FALSE}.}

\item{force}{Add ignored files. Default is FALSE.}
}
\value{
a named vector with the file paths relative to \code{root}. The names
contain the hashes of the files.
}
\description{
A git2rdata object consists of two files.
The \code{".tsv"} file contains the raw data as a plain text tab separated file.
The \code{".yml"} contains the metadata on the columns in plain text YAML format.
See \code{vignette("plain text", package = "git2rdata")} for more details on the
implementation.
}
\note{
\code{..generic} is a reserved name for the metadata and is a forbidden
column name in a \code{data.frame}.
}
\examples{
## on file system

# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# write a dataframe to the directory
write_vc(iris[1:6, ], file = "iris", root = root, sorting = "Sepal.Length")
# check that a data file (.tsv) and a metadata file (.yml) exist.
list.files(root, recursive = TRUE)
# read the git2rdata object from the directory
read_vc("iris", root)

# store a new version with different observations but the same metadata
write_vc(iris[1:5, ], "iris", root)
list.files(root, recursive = TRUE)
# Removing a column requires version requires new metadata.
# Add strict = FALSE to override the existing metadata.
write_vc(
  iris[1:6, -2], "iris", root, sorting = "Sepal.Length", strict = FALSE
)
list.files(root, recursive = TRUE)
# storing the orignal version again requires another update of the metadata
write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Width", strict = FALSE)
list.files(root, recursive = TRUE)
# optimize = FALSE stores the data more verbose. This requires larger files.
write_vc(
  iris[1:6, ], "iris2", root, sorting = "Sepal.Width", optimize = FALSE
)
list.files(root, recursive = TRUE)



## on git repo using a git2r::git-repository

# initialise a git repo using the git2r package
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# store a dataframe in git repo.
write_vc(iris[1:6, ], file = "iris", root = repo, sorting = "Sepal.Length")
# This git2rdata object is not staged by default.
status(repo)
# read a dataframe from a git repo
read_vc("iris", repo)

# store a new version in the git repo and stage it in one go
write_vc(iris[1:5, ], "iris", repo, stage = TRUE)
status(repo)

# store a verbose version in a different gir2data object
write_vc(
  iris[1:6, ], "iris2", repo, sorting = "Sepal.Width", optimize = FALSE
)
status(repo)

# clean up
junk <- file.remove(
  list.files(root, full.names = TRUE, recursive = TRUE), root)
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{list_data}()},
\code{\link{prune_meta}()},
\code{\link{read_vc}()},
\code{\link{relabel}()},
\code{\link{rename_variable}()},
\code{\link{rm_data}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_data.R
\name{list_data}
\alias{list_data}
\title{List Available Git2rdata Files Containing Data}
\usage{
list_data(root = ".", path = ".", recursive = TRUE)
}
\arguments{
\item{root}{the \code{root} of the repository. Either a path or a \code{git-repository}}

\item{path}{relative \code{path} from the \code{root}. Defaults to the \code{root}}

\item{recursive}{logical.  Should the listing recurse into directories?}
}
\value{
A character vector of git2rdata object names, including their
relative path.
}
\description{
The function returns the names of all valid git2rdata objects. This implies
\code{.tsv} files with a matching \strong{valid} metadata file (\code{.yml}). \strong{Invalid}
metadata files result in a warning. The function ignores \strong{valid} metadata
files without matching raw data (\code{.tsv}).
}
\examples{
## on file system

# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# store a dataframe as git2rdata object. Capture the result to minimise
# screen output
junk <- write_vc(iris[1:6, ], "iris", root, sorting = "Sepal.Length")
# write a standard tab separate file (non git2rdata object)
write.table(iris, file = file.path(root, "standard.tsv"), sep = "\t")
# write a YAML file
yml <- list(
  authors = list(
   "Research Institute for Nature and Forest" = list(
       href = "https://www.inbo.be/en")))
yaml::write_yaml(yml, file = file.path(root, "_pkgdown.yml"))

# list the git2rdata objects
list_data(root)
# list the files
list.files(root, recursive = TRUE)

# remove all .tsv files from valid git2rdata objects
rm_data(root, path = ".")
# check the removal of the .tsv file
list.files(root, recursive = TRUE)
list_data(root)

# remove dangling git2rdata metadata files
prune_meta(root, path = ".")
# check the removal of the metadata
list.files(root, recursive = TRUE)
list_data(root)


## on git repo

# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo-")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# store a dataframe
write_vc(iris[1:6, ], "iris", repo, sorting = "Sepal.Length", stage = TRUE)
# check that the dataframe is stored
status(repo)
list_data(repo)

# commit the current version and check the git repo
commit(repo, "add iris data", session = TRUE)
status(repo)

# remove the data files from the repo
rm_data(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# remove dangling metadata
prune_meta(repo, path = ".")
# check the removal
list_data(repo)
status(repo)

# clean up
junk <- file.remove(
  list.files(root, full.names = TRUE, recursive = TRUE), root)
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other storage: 
\code{\link{prune_meta}()},
\code{\link{read_vc}()},
\code{\link{relabel}()},
\code{\link{rename_variable}()},
\code{\link{rm_data}()},
\code{\link{write_vc}()}
}
\concept{storage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upgrade_data.R
\name{upgrade_data}
\alias{upgrade_data}
\alias{upgrade_data.git_repository}
\title{Upgrade Files to the New Version}
\usage{
upgrade_data(file, root = ".", verbose, ..., path)

\method{upgrade_data}{git_repository}(
  file,
  root = ".",
  verbose = TRUE,
  ...,
  path,
  stage = FALSE,
  force = FALSE
)
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.
Defaults to the current working directory (\code{"."}).}

\item{verbose}{display a message with the update status. Defaults to \code{TRUE}.}

\item{...}{parameters used in some methods}

\item{path}{specify \code{path} instead of \code{file} to update all git2rdata objects
in this directory and it's subdirectories. \code{path} is relative to \code{root}. Use
\code{path = "."} to upgrade all git2rdata objects under \code{root}.}

\item{stage}{Logical value indicating whether to stage the changes after
writing the data. Defaults to \code{FALSE}.}

\item{force}{Add ignored files. Default is FALSE.}
}
\value{
the git2rdata object names.
}
\description{
Updates the data written by older versions to the current data format
standard. Works both on a single file and (recursively) on a path. The
\code{".yml"} file must contain a \code{"..generic"} element. \code{upgrade_data()} ignores
all other files.
}
\examples{
# create a directory
root <- tempfile("git2rdata-")
dir.create(root)

# write dataframes to the root
write_vc(iris[1:6, ], file = "iris", root = root, sorting = "Sepal.Length")
write_vc(iris[5:10, ], file = "subdir/iris", root = root,
         sorting = "Sepal.Length")
# upgrade a single git2rdata object
upgrade_data(file = "iris", root = root)
# use path = "." to upgrade all git2rdata objects under root
upgrade_data(path = ".", root = root)

# clean up
junk <- file.remove(list.files(root, full.names = TRUE), root)
}
\seealso{
Other internal: 
\code{\link{is_git2rdata}()},
\code{\link{is_git2rmeta}()},
\code{\link{meta}()}
}
\concept{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/git2rdata_package.R
\docType{package}
\name{git2rdata-package}
\alias{git2rdata}
\alias{git2rdata-package}
\title{git2rdata: Store and Retrieve Data.frames in a Git Repository}
\description{
The git2rdata package is an R package for writing and reading
    dataframes as plain text files.  A metadata file stores important
    information.  1) Storing metadata allows to maintain the classes of
    variables.  By default, git2rdata optimizes the data for file storage.
    The optimization is most effective on data containing factors.  The
    optimization makes the data less human readable.  The user can turn
    this off when they prefer a human readable format over smaller files.
    Details on the implementation are available in vignette("plain_text",
    package = "git2rdata").  2) Storing metadata also allows smaller row
    based diffs between two consecutive commits.  This is a useful feature
    when storing data as plain text files under version control.  Details
    on this part of the implementation are available in
    vignette("version_control", package = "git2rdata").  Although we
    envisioned git2rdata with a git workflow in mind, you can use it in
    combination with other version control systems like subversion or
    mercurial.  3) git2rdata is a useful tool in a reproducible and
    traceable workflow.  vignette("workflow", package = "git2rdata") gives
    a toy example.  4) vignette("efficiency", package = "git2rdata")
    provides some insight into the efficiency of file storage, git
    repository size and speed for writing and reading.  Please cite using
    <doi:10.5281/zenodo.1485309>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://ropensci.github.io/git2rdata/}
  \item Report bugs at \url{https://github.com/ropensci/git2rdata/issues}
}

}
\author{
\strong{Maintainer}: Thierry Onkelinx \email{thierry.onkelinx@inbo.be} (\href{https://orcid.org/0000-0001-8804-4216}{ORCID})

Other contributors:
\itemize{
  \item Floris Vanderhaeghe \email{floris.vanderhaeghe@inbo.be} (\href{https://orcid.org/0000-0002-6378-6229}{ORCID}) [contributor]
  \item Peter Desmet \email{peter.desmet@inbo.be} (\href{https://orcid.org/0000-0002-8442-8025}{ORCID}) [contributor]
  \item Els Lommelen \email{els.lommelen@inbo.be} (\href{https://orcid.org/0000-0002-3481-5684}{ORCID}) [contributor]
  \item Research Institute for Nature and Forest \email{info@inbo.be} [copyright holder, funder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexport.R
\name{push}
\alias{push}
\title{Re-exported Function From \code{git2r}}
\description{
See \code{\link[git2r]{push}} in  \code{git2r}.
}
\seealso{
Other version_control: 
\code{\link{commit}()},
\code{\link{pull}()},
\code{\link{recent_commit}()},
\code{\link{repository}()},
\code{\link{status}()}
}
\concept{version_control}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recent_commit.R
\name{recent_commit}
\alias{recent_commit}
\title{Retrieve the Most Recent File Change}
\usage{
recent_commit(file, root, data = FALSE)
}
\arguments{
\item{file}{the name of the git2rdata object. Git2rdata objects cannot
have dots in their name. The name may include a relative path. \code{file} is a
path relative to the \code{root}.
Note that \code{file} must point to a location within \code{root}.}

\item{root}{The root of a project. Can be a file path or a \code{git-repository}.}

\item{data}{does \code{file} refers to a data object (\code{TRUE}) or to a file
(\code{FALSE})?
Defaults to \code{FALSE}.}
}
\value{
a \code{data.frame} with \code{commit}, \code{author} and \code{when} for the most recent
commit that adds op updates the file.
}
\description{
Retrieve the most recent commit that added or updated a file or git2rdata
object. This does not imply that file still exists at the current HEAD as it
ignores the deletion of files.

Use this information to document the current version of file or git2rdata
object in an analysis. Since it refers to the most recent change of this
file, it remains unchanged by committing changes to other files. You can
also use it to track if data got updated, requiring an analysis to
be rerun. See \code{vignette("workflow", package = "git2rdata")}.
}
\examples{
# initialise a git repo using git2r
repo_path <- tempfile("git2rdata-repo")
dir.create(repo_path)
repo <- git2r::init(repo_path)
git2r::config(repo, user.name = "Alice", user.email = "alice@example.org")

# write and commit a first dataframe
# store the output of write_vc() minimize screen output
junk <- write_vc(iris[1:6, ], "iris", repo, sorting = "Sepal.Length",
                 stage = TRUE)
commit(repo, "important analysis", session = TRUE)
list.files(repo_path)
Sys.sleep(1.1) # required because git doesn't handle subsecond timings

# write and commit a second dataframe
junk <- write_vc(iris[7:12, ], "iris2", repo, sorting = "Sepal.Length",
                 stage = TRUE)
commit(repo, "important analysis", session = TRUE)
list.files(repo_path)
Sys.sleep(1.1) # required because git doesn't handle subsecond timings

# write and commit a new version of the first dataframe
junk <- write_vc(iris[7:12, ], "iris", repo, stage = TRUE)
list.files(repo_path)
commit(repo, "important analysis", session = TRUE)



# find out in which commit a file was last changed

# "iris.tsv" was last updated in the third commit
recent_commit("iris.tsv", repo)
# "iris.yml" was last updated in the first commit
recent_commit("iris.yml", repo)
# "iris2.yml" was last updated in the second commit
recent_commit("iris2.yml", repo)
# the git2rdata object "iris" was last updated in the third commit
recent_commit("iris", repo, data = TRUE)

# remove a dataframe and commit it to see what happens with deleted files
file.remove(file.path(repo_path, "iris.tsv"))
prune_meta(repo, ".")
commit(repo, message = "remove iris", all = TRUE, session = TRUE)
list.files(repo_path)

# still points to the third commit as this is the latest commit in which the
# data was present
recent_commit("iris", repo, data = TRUE)

#' clean up
junk <- file.remove(
  rev(list.files(repo_path, full.names = TRUE, recursive = TRUE,
                 include.dirs = TRUE, all.files = TRUE)),
  repo_path)
}
\seealso{
Other version_control: 
\code{\link{commit}()},
\code{\link{pull}()},
\code{\link{push}()},
\code{\link{repository}()},
\code{\link{status}()}
}
\concept{version_control}
