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
