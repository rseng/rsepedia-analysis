
# gittargets <img src='man/figures/logo-readme.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/486_status.svg)](https://github.com/ropensci/software-review/issues/486)
[![CRAN](https://www.r-pkg.org/badges/version/gittargets)](https://CRAN.R-project.org/package=gittargets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/gittargets/workflows/check/badge.svg)](https://github.com/ropensci/gittargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/gittargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/gittargets)
[![lint](https://github.com/ropensci/gittargets/workflows/lint/badge.svg)](https://github.com/ropensci/gittargets/actions?query=workflow%3Alint)

Pipelines with the [`targets`](https://docs.ropensci.org/targets/) R
package skip steps that are up to already date. Although this behavior
reduces the runtime of subsequent runs, it comes at the cost of
overwriting previous results. Ordinarily, you would need to rerun the
pipeline in order to recover any overwritten targets. However,
`gittargets` preserves historical output, creating version control
snapshots of data store. Each data snapshot remembers the
contemporaneous Git commit of the pipeline source code, so you can
recover the right data when you navigate the Git history. In other
words, `gittargets` makes it possible to switch commits or branches
without invalidating the pipeline. You can simply check out the
up-to-date targets from the past instead of taking the time to recompute
them from scratch.

## Prerequisites

1.  Familiarity with the [R programming
    language](https://www.r-project.org/), covered in [R for Data
    Science](https://r4ds.had.co.nz/).
2.  [Data science workflow management best
    practices](https://rstats.wtf/index.html).
3.  [Git](https://git-scm.com), covered in [Happy Git and GitHub for the
    useR](https://happygitwithr.com).
4.  [`targets`](https://docs.ropensci.org/targets/), which has resources
    on the [documentation
    website](https://docs.ropensci.org/targets/#how-to-get-started).
5.  Familiarity with the [`targets` data
    store](https://books.ropensci.org/targets/files.html#internal-data-files).

## Installation

The package is available to install from any of the following sources.

| Type        | Source   | Command                                                                     |
|-------------|----------|-----------------------------------------------------------------------------|
| Release     | CRAN     | `install.packages("gittargets")`                                            |
| Development | GitHub   | `remotes::install_github("ropensci/gittargets")`                            |
| Development | rOpenSci | `install.packages("gittargets", repos = "https://ropensci.r-universe.dev")` |

You will also need command line Git, available at
<https://git-scm.com/downloads>.[1] Please make sure Git is reachable
from your system path environment variables. To control which Git
executable `gittargets` uses, you may set the `TAR_GIT` environment
variable with `usethis::edit_r_environ()` or `Sys.setenv()`. You will
also need to configure your user name and user email at the global level
using the instructions at
<https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup>
(or `gert::git_config_global_set()`). Run `tar_git_ok()` to check
installation and configuration.

``` r
tar_git_ok()
#> ✓ Git binary: /path/to/git
#> ✓ Git config global user name: your_user_name
#> ✓ Git config global user email: your_email@example.com
#> [1] TRUE
```

There are also backend-specific installation requirements and
recommendations in the [package
vignettes](https://docs.ropensci.org/gittargets/articles/index.html).

## Motivation

Consider an example pipeline with source code in `_targets.R` and output
in the [data
store](https://books.ropensci.org/targets/files.html#internal-data-files).

``` r
# _targets.R
library(targets)
list(
  tar_target(data, airquality),
  tar_target(model, lm(Ozone ~ Wind, data = data)) # Regress on wind speed.
)
```

Suppose you run the pipeline and confirm that all targets are up to
date.

``` r
tar_make()
#> • start target data
#> • built target data
#> • start target model
#> • built target model
#> • end pipeline
```

``` r
tar_outdated()
#> character(0)
```

It is good practice to track the source code in a [version control
repository](https://git-scm.com) so you can revert to previous commits
or branches. However, the [data
store](https://books.ropensci.org/targets/files.html#internal-data-files)
is usually too large to keep in the same repository as the code, which
typically lives in a cloud platform like [GitHub](https://github.com)
where space and bandwidth are pricey. So when you check out an old
commit or branch, you revert the code, but not the data. In other words,
your targets are out of sync and out of date.

``` r
gert::git_branch_checkout(branch = "other-model")
```

``` r
# _targets.R
library(targets)
list(
  tar_target(data, airquality),
  tar_target(model, lm(Ozone ~ Temp, data = data)) # Regress on temperature.
)
```

``` r
tar_outdated()
#> [1] "model"
```

## Usage

With `gittargets`, you can keep your targets up to date even as you
check out code from different commits or branches. The specific steps
depend on the data backend you choose, and each supported backend has a
[package
vignette](https://docs.ropensci.org/gittargets/articles/index.html) with
a walkthrough. For example, the most important steps of the [Git data
backend](https://docs.ropensci.org/gittargets/articles/git.html) are as
follows.

1.  Create the source code and run the pipeline at least once so the
    [data
    store](https://books.ropensci.org/targets/files.html#internal-data-files)
    exists.
2.  `tar_git_init()`: initialize a [Git](https://git-scm.com)/[Git
    LFS](https://git-lfs.github.com) repository for the [data
    store](https://books.ropensci.org/targets/files.html#internal-data-files).
3.  Bring the pipeline up to date (e.g. with
    [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html))
    and commit any changes to the source code.
4.  `tar_git_snapshot()`: create a data snapshot for the current code
    commit.
5.  Develop the pipeline. Creating new code commits and code branches
    early and often, and create data snapshots at key strategic
    milestones.
6.  `tar_git_checkout()`: revert the data to the appropriate prior
    snapshot.

## Performance

`targets` generates a large amount of data in `_targets/objects/`, and
data snapshots and checkouts may take a long time. To work around
performance limitations, you may wish to only snapshot the data at the
most important milestones of your project. Please refer to the [package
vignettes](https://docs.ropensci.org/gittargets/articles/index.html) for
specific recommendations on optimizing performance.

## Future directions

The first data versioning system in `gittargets` uses
[Git](https://git-scm.com), which is designed for source code and may
not scale to enormous amounts of compressed data. Future releases of
`gittargets` may explore alternative data backends more powerful than
[Git LFS](https://git-lfs.github.com).

## Alternatives

Newer versions of the `targets` package (&gt;= 0.9.0) support continuous
data versioning through [Amazon Web Services](https://aws.amazon.com)
for [S3 buckets](https://aws.amazon.com/s3/) with [versioning
enabled](https://docs.aws.amazon.com/AmazonS3/latest/userguide/Versioning.html).
In this approach, `targets` tracks the version ID of each [AWS-backed
target](https://books.ropensci.org/targets/storage_amazon.html). That
way, when the metadata file reverts to a prior version, the pipeline
automatically uses prior versions of targets that were up to date at the
time the metadata was written. This approach has two distinct advantages
over `gittargets`:

1.  Cloud storage reduces the burden of local storage for large data
    pipelines.
2.  Target data is uploaded and tracked continuously, which means the
    user does not need to proactively take data snapshots.

However, not all users have access to [AWS](https://aws.amazon.com), not
everyone is able or willing to pay the monetary costs of cloud storage
for every single version of every single target, and uploads and
downloads to and from the cloud may bottleneck some pipelines.
`gittargets` fills this niche with a data versioning system that is

1.  Entirely local, and
2.  Entirely opt-in: users pick and choose when to register data
    snapshots, which consumes less storage than continuous snapshots or
    continuous cloud uploads to a versioned S3 bucket.

## Code of Conduct

Please note that the `gittargets` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Citation

``` r
citation("gittargets")
#> 
#> To cite gittargets in publications use:
#> 
#>   William Michael Landau (2021). gittargets: Version Control for the
#>   targets Package. https://docs.ropensci.org/gittargets/,
#>   https://github.com/ropensci/gittargets.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {gittargets: Version Control for the Targets Package},
#>     author = {William Michael Landau},
#>     note = {https://docs.ropensci.org/gittargets/, https://github.com/ropensci/gittargets},
#>     year = {2021},
#>   }
```

[1] `gert` does not have these requirements, but `gittargets` does not
exclusively rely on `gert` because `libgit2` does not automatically work
with [Git LFS](https://git-lfs.github.com).
# gittargets 0.0.1.9000

* Hard reset after checkout in `tar_git_checkout()` in order to recover potentially deleted files (#11).

# gittargets 0.0.1

* Join rOpenSci.
* Rewrite README to motivate the use case.
* Remove workflow diagram.
* Simplify snapshot model diagram.
* Fix the documentation of the `ref` argument of `tar_git_checkout()`.
* Add a section to the `git.Rmd` vignette on code merges.
* Allow command line Git tests to run locally on Windows.
* First version.
# Contributing

Development is a community effort, and we welcome participation.

## Code of Conduct

Please note that the `gittargets` project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Discussions

At <https://github.com/ropensci/gittargets/discussions>, you can post general questions, brainstorm ideas, and ask for help.

## Issues

<https://github.com/ropensci/gittargets/issues> is for bug reports, performance issues, maintenance tasks, and feature requests. When you post, please abide by the following guidelines.

* Before posting a new issue or discussion topic, please take a moment to search for existing similar threads in order to avoid duplication.
* For bug reports: if you can, please install the latest GitHub version of `gittargets` (i.e. `remotes::install_github("ropensci/gittargets")`) and verify that the issue still persists.
* Describe your issue in prose as clearly and concisely as possible.
* For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Development

External code contributions are extremely helpful in the right circumstances. Here are the recommended steps.

1. Prior to contribution, please propose your idea in a discussion topic or issue thread so you and the maintainer can define the intent and scope of your work.
2. [Fork the repository](https://help.github.com/articles/fork-a-repo/).
3. Follow the [GitHub flow](https://guides.github.com/introduction/flow/index.html) to create a new branch, add commits, and open a pull request.
4. Discuss your code with the maintainer in the pull request thread.
5. If everything looks good, the maintainer will merge your code into the project.

Please also follow these additional guidelines.

* Respect the architecture and reasoning of the package. Depending on the scope of your work, you may want to read the design documents (package vignettes).
* If possible, keep contributions small enough to easily review manually. It is okay to split up your work into multiple pull requests.
* Format your code according to the [tidyverse style guide](https://style.tidyverse.org/) and check your formatting with the `lint_package()` function from the [`lintr`](https://github.com/jimhester/lintr) package.
* For new features or functionality, add tests in `tests`. Tests that can be automated should go in `tests/testthat/`. Tests that cannot be automated should go in `tests/interactive/`. For features affecting performance, it is good practice to add profiling studies to `tests/performance/`.
* Check code coverage with `covr::package_coverage()`. Automated tests should cover all the new or changed functionality in your pull request.
* Run overall package checks with `devtools::check()` and `goodpractice::gp()`
* Describe your contribution in the project's [`NEWS.md`](https://github.com/ropensci/gittargets/blob/main/NEWS.md) file. Be sure to mention relevent GitHub issue numbers and your GitHub name as done in existing news entries.
* If you feel contribution is substantial enough for official author or contributor status, please add yourself to the `Authors@R` field of the [`DESCRIPTION`](https://github.com/ropensci/gittargets/blob/main/DESCRIPTION) file.
# Prework

* [ ] I understand and agree to the [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
* [ ] I have already submitted a [discussion topic](https://github.com/ropensci/gittargets/discussions) or [issue](https://github.com/ropensci/gittargets/issues) to discuss my idea with the maintainer.

# Related GitHub issues and pull requests

* Ref: #

# Summary

Please explain the purpose and scope of your contribution.
---
name: Maintenance
about: "Something in gittargets needs work: updates, documentation, etc. Not a bug, performance issue, or new feature."
title: ""
labels: "type: maintenance"
assignees: ""
---

## Prework

* [ ] Read and agree to the [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/gittargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/gittargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the issue.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.

## Reproducible example

* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).
---
name: Bug
about: Please do not submit a bug report unless your issue is a genuine bug in gittargets and not a known limitation, usage error, or issue from another package that gittargets depends on.
title: ""
labels: "type: bug"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/gittargets/blob/main/CONTRIBUTING.md).
* [ ] Confirm that your issue is a genuine bug in the `gittargets` package itself and not a user error, known limitation, or issue from another package that `gittargets` depends on. For example, if you get errors running `tar_make_clustermq()`, try isolating the problem in a reproducible example that runs `clustermq` and not `gittargets`. And for miscellaneous troubleshooting, please post to [discussions](https://github.com/ropensci/gittargets/discussions) instead of [issues](https://github.com/ropensci/gittargets/issues).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/gittargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the bug.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Expected result

What should have happened? Please be as specific as possible.

## Diagnostic information

* A [reproducible example](https://github.com/tidyverse/reprex).
* Session info, available through `sessionInfo()` or [`reprex(si = TRUE)`](https://github.com/tidyverse/reprex).
* A stack trace from `traceback()` or `rlang::trace_back()`.
* The [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity) of the GitHub commit of `gittargets` currently installed. `packageDescription("gittargets")$GithubSHA1` shows you this.
---
name: New feature
about: Suggest a new feature.
title: ""
labels: "type: new feature"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/gittargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/gittargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] New features take time and effort to create, and they take even more effort to maintain. So if the purpose of the feature is to resolve a struggle you are encountering personally, please consider first posting a "trouble" or "other" issue so we can discuss your use case and search for existing solutions first.
* [ ] Format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Proposal

Please describe the new feature. If applicable, write a minimal example in R code or pseudo-code to show input, usage, and desired output.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.
---
name: Performance
about: "Runtime, memory, or storage inefficiency"
title: ""
labels: "topic: performance"
assignees: wlandau

---

## Prework

* [ ] Read and agree to the [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/gittargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/gittargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the performance issue.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Benchmarks

How poorly does `gittargets` perform? To find out, we recommend you use the [`proffer`](https://github.com/r-prof/proffer) package and take screenshots of the results displayed in your browser.

```r
library(gittargets)
library(proffer)
px <- pprof({
  # All your gittargets code goes here.
})
```
The image of the guitar is in the public domain (no copyright, unlimited commercial use): https://openclipart.org/detail/10117/guitar-1
---
output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# gittargets <img src='man/figures/logo-readme.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/486_status.svg)](https://github.com/ropensci/software-review/issues/486)
[![CRAN](https://www.r-pkg.org/badges/version/gittargets)](https://CRAN.R-project.org/package=gittargets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/gittargets/workflows/check/badge.svg)](https://github.com/ropensci/gittargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/gittargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/gittargets)
[![lint](https://github.com/ropensci/gittargets/workflows/lint/badge.svg)](https://github.com/ropensci/gittargets/actions?query=workflow%3Alint)

Pipelines with the [`targets`](https://docs.ropensci.org/targets/) R package skip  steps that are up to already date. Although this behavior reduces the runtime of subsequent runs, it comes at the cost of overwriting previous results. Ordinarily, you would need to rerun the pipeline in order to recover any overwritten targets. However, `gittargets` preserves historical output, creating version control snapshots of data store. Each data snapshot remembers the contemporaneous Git commit of the pipeline source code, so you can recover the right data when you navigate the Git history. In other words, `gittargets` makes it possible to switch commits or branches without invalidating the pipeline. You can simply check out the up-to-date targets from the past instead of taking the time to recompute them from scratch.

## Prerequisites

1. Familiarity with the [R programming language](https://www.r-project.org/), covered in [R for Data Science](https://r4ds.had.co.nz/).
1. [Data science workflow management best practices](https://rstats.wtf/index.html).
1. [Git](https://git-scm.com), covered in [Happy Git and GitHub for the useR](https://happygitwithr.com).
1. [`targets`](https://docs.ropensci.org/targets/), which has resources on the [documentation website](https://docs.ropensci.org/targets/#how-to-get-started).
1. Familiarity with the [`targets` data store](https://books.ropensci.org/targets/files.html#internal-data-files).

## Installation

The package is available to install from any of the following sources.

Type | Source | Command
---|---|---
Release | CRAN | `install.packages("gittargets")`
Development | GitHub | `remotes::install_github("ropensci/gittargets")`
Development | rOpenSci | `install.packages("gittargets", repos = "https://ropensci.r-universe.dev")`

You will also need command line Git, available at <https://git-scm.com/downloads>.^[`gert` does not have these requirements, but `gittargets` does not exclusively rely on `gert` because `libgit2` does not automatically work with [Git LFS](https://git-lfs.github.com).] Please make sure Git is reachable from your system path environment variables. To control which Git executable `gittargets` uses, you may set the `TAR_GIT` environment variable with `usethis::edit_r_environ()` or `Sys.setenv()`. You will also need to configure your user name and user email at the global level using the instructions at <https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup> (or `gert::git_config_global_set()`). Run `tar_git_ok()` to check installation and configuration.

```{r, eval = FALSE}
tar_git_ok()
#> ✓ Git binary: /path/to/git
#> ✓ Git config global user name: your_user_name
#> ✓ Git config global user email: your_email@example.com
#> [1] TRUE
```

There are also backend-specific installation requirements and recommendations in the [package vignettes](https://docs.ropensci.org/gittargets/articles/index.html).

## Motivation

Consider an example pipeline with source code in `_targets.R` and output in the [data store](https://books.ropensci.org/targets/files.html#internal-data-files). 

```{r, eval = FALSE}
# _targets.R
library(targets)
list(
  tar_target(data, airquality),
  tar_target(model, lm(Ozone ~ Wind, data = data)) # Regress on wind speed.
)
```

Suppose you run the pipeline and confirm that all targets are up to date.

```{r, eval = FALSE}
tar_make()
#> • start target data
#> • built target data
#> • start target model
#> • built target model
#> • end pipeline
```

```{r, eval = FALSE}
tar_outdated()
#> character(0)
```

It is good practice to track the source code in a [version control repository](https://git-scm.com) so you can revert to previous commits or branches. However, the [data store](https://books.ropensci.org/targets/files.html#internal-data-files) is usually too large to keep in the same repository as the code, which typically lives in a cloud platform like [GitHub](https://github.com) where space and bandwidth are pricey. So when you check out an old commit or branch, you revert the code, but not the data. In other words, your targets are out of sync and out of date.

```{r, eval = FALSE}
gert::git_branch_checkout(branch = "other-model")
```

```{r, eval = FALSE}
# _targets.R
library(targets)
list(
  tar_target(data, airquality),
  tar_target(model, lm(Ozone ~ Temp, data = data)) # Regress on temperature.
)
```

```{r, eval = FALSE}
tar_outdated()
#> [1] "model"
```

## Usage

With `gittargets`, you can keep your targets up to date even as you check out code from different commits or branches. The specific steps depend on the data backend you choose, and each supported backend has a [package vignette](https://docs.ropensci.org/gittargets/articles/index.html) with a walkthrough. For example, the most important steps of the [Git data backend](https://docs.ropensci.org/gittargets/articles/git.html) are as follows.

1. Create the source code and run the pipeline at least once so the [data store](https://books.ropensci.org/targets/files.html#internal-data-files) exists.
1. `tar_git_init()`: initialize a [Git](https://git-scm.com)/[Git LFS](https://git-lfs.github.com) repository for the [data store](https://books.ropensci.org/targets/files.html#internal-data-files).
1. Bring the pipeline up to date (e.g. with [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html)) and commit any changes to the source code.
1. `tar_git_snapshot()`: create a data snapshot for the current code commit.
1. Develop the pipeline. Creating new code commits and code branches early and often, and create data snapshots at key strategic milestones.
1. `tar_git_checkout()`: revert the data to the appropriate prior snapshot.

## Performance

`targets` generates a large amount of data in `_targets/objects/`, and data snapshots and checkouts may take a long time. To work around performance limitations, you may wish to only snapshot the data at the most important milestones of your project. Please refer to the [package vignettes](https://docs.ropensci.org/gittargets/articles/index.html) for specific recommendations on optimizing performance. 

## Future directions

The first data versioning system in `gittargets` uses [Git](https://git-scm.com), which is designed for source code and may not scale to enormous amounts of compressed data. Future releases of `gittargets` may explore alternative data backends more powerful than [Git LFS](https://git-lfs.github.com).

## Alternatives

Newer versions of the `targets` package (>= 0.9.0) support continuous data versioning through [Amazon Web Services](https://aws.amazon.com) for [S3 buckets](https://aws.amazon.com/s3/) with [versioning enabled](https://docs.aws.amazon.com/AmazonS3/latest/userguide/Versioning.html). In this approach, `targets` tracks the version ID of each [AWS-backed target](https://books.ropensci.org/targets/storage_amazon.html). That way, when the metadata file reverts to a prior version, the pipeline automatically uses prior versions of targets that were up to date at the time the metadata was written. This approach has two distinct advantages over `gittargets`:

1. Cloud storage reduces the burden of local storage for large data pipelines.
2. Target data is uploaded and tracked continuously, which means the user does not need to proactively take data snapshots. 

However, not all users have access to [AWS](https://aws.amazon.com), not everyone is able or willing to pay the monetary costs of cloud storage for every single version of every single target, and uploads and downloads to and from the cloud may bottleneck some pipelines. `gittargets` fills this niche with a data versioning system that is 

1. Entirely local, and 
2. Entirely opt-in: users pick and choose when to register data snapshots, which consumes less storage than continuous snapshots or continuous cloud uploads to a versioned S3 bucket.

## Code of Conduct

Please note that the `gittargets` project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Citation

```{r, warning = FALSE}
citation("gittargets")
```
---
title: "Tutorial: Git data backend"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Tutorial: Git data backend}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(gert)
library(gittargets)
library(targets)
tmp <- tempfile()
dir.create(tmp)
knitr::opts_knit$set(root.dir = tmp)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("NOT_CRAN"), "true") &&
    tar_git_ok(verbose = FALSE)
)
```

This tutorial shows how to use `gittargets` with the Git-based data versioning backend. Before proceeding, please read the [README](https://github.com/ropensci/gittargets/blob/main/README.md) file or [documentation website front page](https://docs.ropensci.org/gittargets/) for an overview of the package.

## Installation

Please begin with the installation instructions on the [documentation website](https://docs.ropensci.org/gittargets/). In addition, if your `targets` pipeline generates large data files, consider installing [Git LFS](https://git-lfs.github.com). The Git data backend in `gittargets` automatically opts into [Git LFS](https://git-lfs.github.com), so you should not need to do any manual configuration to reap the performance benefits.

## Remotes

This backend is uses local [Git](https://git-scm.com) only. It is possible to push the snapshotted [data store](https://books.ropensci.org/targets/files.html#internal-data-files) to a cloud service like [GitHub](https://github.com), [GitLab](https://about.gitlab.com), or [Bitbucket](https://bitbucket.org), but this is the user's responsibility. Pipelines usually generate large data files, and [GitHub](https://github.com) and its peers have file size limitations. Also, `gittargets` automatically opts into [Git LFS](https://git-lfs.github.com) locally (unless `git_lfs` is `FALSE` in `tar_git_init()`), and [Git LFS](https://git-lfs.github.com) on the cloud is a paid service.

## Overall workflow

The most important steps of the Git data backend are as follows. The rest of this vignette walks through these steps in greater depth.

1. Create the source code and run the pipeline at least once so the [data store](https://books.ropensci.org/targets/files.html#internal-data-files) exists.
1. `tar_git_init()`: initialize a [Git](https://git-scm.com)/[Git-LFS](https://git-lfs.github.com) repository for the [data store](https://books.ropensci.org/targets/files.html#internal-data-files).
1. Bring the pipeline up to date (e.g. with [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html)) and commit any changes to the source code.
1. `tar_git_snapshot()`: create a data snapshot for the current code commit.
1. Develop the pipeline. Creating new code commits and code branches early and often, and create data snapshots at key strategic milestones.
1. `tar_git_checkout()`: revert the data to the appropriate prior snapshot.

## Create code

To begin development, we write `_targets.R` file for a [`targets`](https://docs.ropensci.org/targets/) pipeline. [`targets`](https://docs.ropensci.org/targets/) can handle large complex pipelines for [machine learning](https://github.com/wlandau/targets-keras), [Bayesian data analysis](https://github.com/wlandau/rmedicine2021-pipeline), and much more. However, this tutorial focuses on a much simpler pipeline for the sake of pedagogical simplicity.

```{r, eval = FALSE}
# _targets.R
library(targets)
list(
  tar_target(data, datasets::airquality),
  tar_target(result, summary(data))
)
```

```{r, echo = FALSE}
tar_script(
  list(
    tar_target(data, datasets::airquality),
    tar_target(result, summary(data))
  )
)
```

## Run pipeline

With our target script in hand, we run the pipeline.^[<https://books.ropensci.org/targets/hpc.html> describes heavy-duty alternatives to `tar_make()`.]

```{r}
tar_make()
```

We inspect the output with `tar_read()`.

```{r}
tar_read(result)
```

## Commit code

We usually iterate between writing code and running the pipeline until we have a decent set of results. After that, we commit the code to a [Git](https://git-scm.com/docs/git-pack-refs) repository, which may or may not live on [GitHub](https://github.com).^[Alternatives to GitHub include GitLab and Bitbucket.] [Happy Git with R](https://happygitwithr.com) is a great way to learn Git, and the [`gert`](https://docs.ropensci.org/gert/) package is a convenient way to interact with Git from R.

```{r, message = FALSE, output = FALSE, results = "hide"}
library(gert)
git_init()
git_add("_targets.R")
git_commit("Begin analyzing the airquality dataset")
git_branch_create("airquality")
```

## Snapshot data

Before we snapshot the data, we should check that the code is up to date in the Git repository and the targets are up to date in the pipeline. The `tar_git_status()` function is an easy way to do this.^[Helper functions `tar_git_status_code()`, `tar_git_status_targets()`, and `tar_git_status_data()` each generate a piece of the `tar_git_status()` output.]

```{r}
tar_git_status()
```

Our code and pipeline look ready for a data snapshot. First, we initialize the data repository with `tar_git_init()`. `tar_git_init()` writes a `.gitattributes` file in the data store that automatically opts into [Git LFS](https://git-lfs.github.com). If you have [Git LFS](https://git-lfs.github.com) but do not wish to use it, please remove the `.gitattributes` after calling `tar_git_init()`.

```{r}
tar_git_init()
```

Then, we create our first data commit with `tar_git_snapshot()`.^[Ordinarily, `tar_git_snapshot()` shows runs `tar_git_status()` and prompts the user to confirm the snapshot. But in this example, we skip this step.]

```{r, eval = FALSE}
tar_git_snapshot()
```

```{r, echo = FALSE}
tar_git_snapshot(status = FALSE)
```

### Snapshot model

![](./snapshot-model-git.png)

In the Git data backend, a data snapshot is a special kind of Git commit (gray boxes above). Each data commit is part of a data branch (vertical dashed lines above), and each data branch is specific to the current code commit (green and brown boxes above). In fact, each data branch name is of the form `"code=<SHA1>"`, where `<SHA1>` is the Git SHA1 hash of the corresponding code commit. You can always create a data snapshot, but it will supersede any prior data snapshot you already have for the current code commit. To revert to a previous data snapshots for a given code snapshot, you will need to manually enter the repository and check out the relevant data commit.

## Repeat

Development typically happens in cycles: develop the code, run the pipeline, commit the code, snapshot the data, and repeat. Not all code commits need a data snapshot, especially if the [`targets`](https://docs.ropensci.org/targets/) pipeline generates a lot of data. But even then, it is helpful to snapshot the data at key milestones, e.g. if an alternative research question comes up and it is desirable to create a new Git branch for the code. For example, suppose we wish to apply the same pipeline to a different dataset. The code changes:

```{r, eval = FALSE}
# _targets.R
library(targets)
list(
  tar_target(data, datasets::UKgas), # different dataset
  tar_target(result, summary(data))
)
```

```{r, echo = FALSE}
tar_script(
  list(
    tar_target(data, datasets::UKgas),
    tar_target(result, summary(data))
  )
)
```

We run the pipeline and inspect the new output.

```{r}
tar_make()
```

```{r}
tar_read(result)
```

We put the code in a new Git branch.

```{r}
git_branch_create("UKgas")
git_add("_targets.R")
git_commit("Switch to UKgas dataset")
```

Finally, we create a data snapshot for the new code commit.

```{r, eval = FALSE}
tar_git_snapshot()
```

```{r, echo = FALSE}
tar_git_snapshot(status = FALSE)
```

## View log

Now, suppose we want to switch the project back to the original dataset (`airquality`). To transition completely, we need to revert both the code and the data. If we only revert the code, then the data store will sill reflect the `UKgas` dataset, and none of our targets will be up to date. At this point, it is a good time to pause and check the `gittargets` log to see which code commits have available data snapshots.^[If you chose not to call `tar_git_snapshot()` for some code commits, then not all your code commits will have available data snapshots.]

```{r}
tar_git_log()
```

## Check out code

To check out the old `airquality` code, we use `gert::git_branch_checkout()`.

```{r, message = FALSE, output = FALSE, results = "hide"}
git_branch_checkout("airquality")
```

But because we did not revert the data, our results still reflect the `UKgas` dataset.

```{r}
tar_read(result)
```

Thus, all our targets are out of date.

```{r}
tar_outdated()
```

## Check out data

To bring our targets back up to date with the `airquality` data, instead of beginning a potentially long computation with `tar_make()`, we can check out the data snapshot that matches our current code commit.

```{r}
tar_git_checkout()
```

Now, our results reflect the `airquality` dataset we previously analyzed.

```{r}
tar_read(result)
```

And all our targets are up to date.

```{r}
tar_outdated()
```

## Merges

It is common to [merge](https://git-scm.com/docs/git-merge) code branches into one another. When this happens, a new merge commit is created in the code repository, and the data repository remains unchanged. In fact, the only change is that the code repository is now at a new commit that has no data snapshot yet. If you wish, simply create a new data snapshot with `tar_git_snapshot()`. If the code commit immediately prior had an up-to-date data snapshot of its own, then the new snapshot for the merge commit should cost little storage or runtime.

## Custom data files

Only files inside the `targets` data store are tracked in a `gittargets` data snapshot. If your pipeline requires custom external files, you may put them in a folder called `_targets/user/`. That way, `gittargets` will automatically put them under data version control in the next snapshot.

## Performance

If your `targets` pipeline generates large data files, consider installing [Git LFS](https://git-lfs.github.com). Once you install [Git LFS](https://git-lfs.github.com), it should just work on your project right out of the box because `tar_git_init()` writes the following to `_targets/.gitattributes`:

```{sh, eval = FALSE}
objects filter=lfs diff=lfs merge=lfs -text
```

In addition, every data snapshot with `tar_git_snapshot()` creates a new Git branch. With thousands of commits and thus thousands of branches, performance may suffer unless you ensure `pack_refs` is `TRUE` in the function call (default).^[Alternatively, you can call `tar_git_snapshot(pack_refs = FALSE)` and then afterwards run `git pack-refs --all`](https://git-scm.com/docs/git-pack-refs) in the command line with your working directory inside `_targets/`.]
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_status_code.R
\name{tar_git_status_code}
\alias{tar_git_status_code}
\title{Status of the code repository (Git)}
\usage{
tar_git_status_code(code = getwd())
}
\arguments{
\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}
}
\value{
If the code repository exists, the return value is the data frame
produced by \code{gert::git_status(repo = code)}. If the code has no Git
repository, then the return value is \code{NULL}.
}
\description{
Show the Git status of the code repository.
}
\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's file space.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
list.files("_targets", all.files = TRUE)
gert::git_init()
tar_git_init()
tar_git_status_code()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_status_targets.R
\name{tar_git_status_targets}
\alias{tar_git_status_targets}
\title{Status of the targets (Git)}
\usage{
tar_git_status_targets(
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store"),
  reporter = targets::tar_config_get("reporter_outdated"),
  envir = parent.frame(),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter)
)
}
\arguments{
\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[targets:tar_script]{tar_script()}},
\code{\link[targets:tar_config_get]{tar_config_get()}}, and \code{\link[targets:tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[targets:tar_config_get]{tar_config_get()}} and \code{\link[targets:tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[targets:tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}
}
\value{
A \code{tibble} with the names of outdated targets.
}
\description{
Show which targets are outdated.
}
\details{
This function has prettier output than \code{targets::tar_outdated()},
and it mainly serves \code{\link[=tar_git_status]{tar_git_status()}}.
}
\examples{
targets::tar_dir({ # Containing code does not modify the user's file space.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
list.files("_targets", all.files = TRUE)
tar_git_status_targets()
})
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_snapshot.R
\name{tar_git_snapshot_menu}
\alias{tar_git_snapshot_menu}
\title{Data snapshot menu (Git)}
\usage{
tar_git_snapshot_menu(
  commit,
  message,
  code,
  script,
  store,
  stash_gitignore,
  reporter,
  envir,
  callr_function,
  callr_arguments
)
}
\arguments{
\item{commit}{Character of length 1, Git SHA1 hash of the code commit
that will correspond to the data snapshot (if created).}

\item{message}{Optional Git commit message of the data snapshot.
If \code{NULL}, then the message is the Git commit message of the
matching code commit.}

\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[targets:tar_script]{tar_script()}},
\code{\link[targets:tar_config_get]{tar_config_get()}}, and \code{\link[targets:tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{stash_gitignore}{Logical of length 1, whether to temporarily
stash the \code{.gitignore} file of the data store. See the
"Stashing .gitignore" section for details.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[targets:tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}
}
\value{
Integer of length 1: \code{2L} if the user agrees to snapshot,
\code{1L} if the user declines.
}
\description{
Check the project status and show an interactive menu
for \code{\link[=tar_git_snapshot]{tar_git_snapshot()}}.
}
\examples{
# See the examples of tar_git_snapshot().
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_ok.R
\name{tar_git_ok}
\alias{tar_git_ok}
\title{Check Git}
\usage{
tar_git_ok(verbose = TRUE)
}
\arguments{
\item{verbose}{Whether to print messages to the console.}
}
\value{
Logical of length 1, whether Git is installed and configured
correctly.
}
\description{
Check if Git is installed and if \code{user.name} and \code{user.email}
are configured globally.
}
\details{
You can install Git from \url{https://git-scm.com/downloads/}
and configure your identity using the instructions at
\url{https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup}.
You may find it convenient to run \code{gert::git_config_global()}
with \code{name} equal to \code{user.name} and \code{user.email}.
}
\examples{
tar_git_ok()
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_snapshot.R
\name{tar_git_snapshot}
\alias{tar_git_snapshot}
\title{Snapshot the data repository (Git).}
\usage{
tar_git_snapshot(
  message = NULL,
  ref = "HEAD",
  code = getwd(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store"),
  stash_gitignore = TRUE,
  reporter = targets::tar_config_get("reporter_outdated"),
  envir = parent.frame(),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  status = interactive(),
  force = FALSE,
  pack_refs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{message}{Optional Git commit message of the data snapshot.
If \code{NULL}, then the message is the Git commit message of the
matching code commit.}

\item{ref}{Character of length 1, reference
(branch name, Git SHA1 hash, etc.) of the code commit
that will map to the new data snapshot. Defaults to the commit
checked out right now.}

\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[targets:tar_script]{tar_script()}},
\code{\link[targets:tar_config_get]{tar_config_get()}}, and \code{\link[targets:tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{stash_gitignore}{Logical of length 1, whether to temporarily
stash the \code{.gitignore} file of the data store. See the
"Stashing .gitignore" section for details.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[targets:tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{status}{Logical of length 1, whether to print the project status
with \code{\link[=tar_git_status]{tar_git_status()}} and ask whether a snapshot should be created.}

\item{force}{Logical of length 1. Force checkout the data branch
of an existing data snapshot of the current code commit?}

\item{pack_refs}{Logical of length 1, whether to run \verb{git pack-refs --all}
in the data store after taking the snapshot. Packing references
improves efficiency when the number of snapshots is large.
Learn more at \url{https://git-scm.com/docs/git-pack-refs}.}

\item{verbose}{Logical of length 1, whether to print R console messages
confirming that a snapshot was created.}
}
\description{
Snapshot the Git data repository of a \code{targets} project.
}
\details{
A Git-backed \code{gittargets} data snapshot is a special kind of
Git commit. Every data commit is part of a branch specific to
the current code commit.
That way, when you switch branches or commits in the code,
\code{tar_git_checkout()} checks out the latest data snapshot
that matches the code in your workspace.
That way, your targets can stay up to date even as you
transition among multiple branches.
}
\section{Stashing .gitignore}{

The \code{targets} package writes a \code{.gitignore} file to new data stores
in order to prevent accidental commits to the code Git repository.
Unfortunately, for \code{gittargets}, this automatic \code{.gitignore} file
interferes with proper data versioning. So by default, \code{gittargets}
temporarily stashes it to a hidden file called \code{.gittargets_gitignore}
inside the data store. If your R program crashes while the stash
is active, you can simply move it manually back to \code{.gitignore}
or run \code{tar_git_status_data()} to restore the stash automatically
if no \code{.gitignore} already exists.
}

\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's filespace.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
gert::git_init()
gert::git_add("_targets.R")
gert::git_commit("First commit")
tar_git_init()
tar_git_snapshot(status = FALSE)
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_init.R
\name{tar_git_init}
\alias{tar_git_init}
\title{Initialize a data repository (Git).}
\usage{
tar_git_init(
  store = targets::tar_config_get("store"),
  stash_gitignore = TRUE,
  git_lfs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{stash_gitignore}{Logical of length 1, whether to temporarily
stash the \code{.gitignore} file of the data store. See the
"Stashing .gitignore" section for details.}

\item{git_lfs}{Logical, whether to automatically opt into Git LFS to track
large files in \verb{_targets/objects} more efficiently. If \code{TRUE}
and Git LFS is installed, it should work automatically. If \code{FALSE},
you can always opt in later by running \verb{git lfs track objects}
inside the data store.}

\item{verbose}{Logical of length 1, whether to print messages to the
R console.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Initialize a Git repository for a \code{targets} data store.
}
\details{
\code{tar_git_init()} also writes a \code{.gitattributes} file to the
store to automatically track target output date with \code{git-lfs}
if it is installed on your system.
}
\section{Stashing .gitignore}{

The \code{targets} package writes a \code{.gitignore} file to new data stores
in order to prevent accidental commits to the code Git repository.
Unfortunately, for \code{gittargets}, this automatic \code{.gitignore} file
interferes with proper data versioning. So by default, \code{gittargets}
temporarily stashes it to a hidden file called \code{.gittargets_gitignore}
inside the data store. If your R program crashes while the stash
is active, you can simply move it manually back to \code{.gitignore}
or run \code{tar_git_status_data()} to restore the stash automatically
if no \code{.gitignore} already exists.
}

\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's file space.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
tar_git_init()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_package.R
\docType{package}
\name{gittargets-package}
\alias{gittargets-package}
\title{targets: Dynamic Function-Oriented Make-Like Declarative Pipelines for R}
\description{
Pipelines with the \code{targets} R package skip  steps that
are up to already date. Although this behavior reduces the runtime
of subsequent runs, it comes at the cost of overwriting previous
results. So if the pipeline source code is under version control,
and if you revert to a previous commit or branch,
the data will no longer be up to date with the code you
just checked out. Ordinarily, you would need to rerun the
pipeline in order to recover the targets you had before.
However, \code{gittargets} preserves historical output,
creating version control snapshots of data store.
Each data snapshot remembers the contemporaneous Git commit
of the pipeline source code, so you can recover the right
data when you navigate the Git history. In other words,
\code{gittargets} makes it possible to switch commits or branches
without invalidating the pipeline. You can simply check out
the up-to-date targets from the past instead of taking the
time to recompute them from scratch.
}
\concept{help}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_checkout.R
\name{tar_git_checkout}
\alias{tar_git_checkout}
\title{Check out a snapshot of the data (Git)}
\usage{
tar_git_checkout(
  ref = "HEAD",
  code = getwd(),
  store = targets::tar_config_get("store"),
  force = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{ref}{Character of length 1. SHA1 hash, branch name,
or other reference in the code repository
that points to a code commit. (You can also identify the code
commit by supplying a data branch of the form \verb{code=<SHA1>}.)
Defaults to \code{"HEAD"}, which points to the currently
checked out code commit.

Once the desired code commit is identified,
\code{tar_git_snapshot()} checks out the latest corresponding data snapshot.
There may be earlier data snapshots corresponding to this code commit,
but \code{tar_git_snapshot()} only checks out the latest one.
To check out an earlier superseded data snapshot,
you will need to manually use command line Git in the data repository.

If \code{tar_git_snapshot()} cannot find a data snapshot for the
desired code commit, then it will throw an error.
For a list of commits in the current code branch
that have available data snapshots, see the \code{commit_code}
column of the output of \code{\link[=tar_git_log]{tar_git_log()}}.}

\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{force}{ignore conflicts and overwrite modified files}

\item{verbose}{Logical of length 1, whether to print R console messages
confirming that a snapshot was created.}
}
\value{
Nothing (invisibly).
}
\description{
Check out a snapshot of the data associated with
a particular code commit (default: \code{HEAD}).
}
\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's filespace.
# Work on an initial branch.
targets::tar_script(tar_target(data, "old_data"))
targets::tar_make()
targets::tar_read(data) # "old_data"
gert::git_init()
gert::git_add("_targets.R")
gert::git_commit("First commit")
gert::git_branch_create("old_branch")
tar_git_init()
# Work on a new branch.
tar_git_snapshot(status = FALSE, verbose = FALSE)
targets::tar_script(tar_target(data, "new_data"))
targets::tar_make()
targets::tar_read(data) # "new_data"
gert::git_branch_create("new_branch")
gert::git_add("_targets.R")
gert::git_commit("Second commit")
tar_git_snapshot(status = FALSE, verbose = FALSE)
# Go back to the old branch.
gert::git_branch_checkout("old_branch")
# The target is out of date because we only reverted the code.
targets::tar_outdated()
# But tar_git_checkout() lets us restore the old version of the data!
tar_git_checkout()
targets::tar_read(data) # "old_data"
# Now, the target is up to date! And we did not even have to rerun it!
targets::tar_outdated()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_log.R
\name{tar_git_log}
\alias{tar_git_log}
\title{Data snapshots of a code branch (Git)}
\usage{
tar_git_log(
  code = getwd(),
  store = targets::tar_config_get("store"),
  branch = gert::git_branch(repo = code),
  max = 100
)
}
\arguments{
\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{branch}{Character of length 1, name of the code repository branch
to query. Defaults to the currently checked-out code branch.}

\item{max}{Positive numeric of length 1, maximum number of code commits
to inspect for the given branch.}
}
\value{
A data frame of information about
data snapshots and code commits.
}
\description{
Show all the data snapshots of a code branch.
}
\details{
By design, \code{tar_git_log()} only queries a single
code branch at a time. This allows \code{tar_git_log()}
to report more detailed information about the snapshots
of the given code branch.
To query all data snapshots over all branches, simply run
\code{gert::git_branch_list(local = TRUE, repo = "_targets")}.
The valid snapshots show \code{"code=<SHA1>"} in the \code{name} column,
where \verb{<SHA1>} is the Git commit hash of the code commit
corresponding to the data snapshot.
}
\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's filespace.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
gert::git_init()
gert::git_add("_targets.R")
gert::git_commit("First commit")
tar_git_init()
tar_git_snapshot(status = FALSE, verbose = FALSE)
tar_git_log()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_status.R
\name{tar_git_status}
\alias{tar_git_status}
\title{Status of the project (Git)}
\usage{
tar_git_status(
  code = getwd(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store"),
  stash_gitignore = TRUE,
  reporter = targets::tar_config_get("reporter_outdated"),
  envir = parent.frame(),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter)
)
}
\arguments{
\item{code}{Character of length 1, directory path to the code repository,
usually the root of the \code{targets} project.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[targets:tar_script]{tar_script()}},
\code{\link[targets:tar_config_get]{tar_config_get()}}, and \code{\link[targets:tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{stash_gitignore}{Logical of length 1, whether to temporarily
stash the \code{.gitignore} file of the data store. See the
"Stashing .gitignore" section for details.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[targets:tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}
}
\value{
\code{NULL} (invisibly). Status information is printed
to the R console.
}
\description{
Print the status of the code repository,
the data repository, and the targets.
}
\section{Stashing .gitignore}{

The \code{targets} package writes a \code{.gitignore} file to new data stores
in order to prevent accidental commits to the code Git repository.
Unfortunately, for \code{gittargets}, this automatic \code{.gitignore} file
interferes with proper data versioning. So by default, \code{gittargets}
temporarily stashes it to a hidden file called \code{.gittargets_gitignore}
inside the data store. If your R program crashes while the stash
is active, you can simply move it manually back to \code{.gitignore}
or run \code{tar_git_status_data()} to restore the stash automatically
if no \code{.gitignore} already exists.
}

\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's files pace.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
list.files("_targets", all.files = TRUE)
gert::git_init()
tar_git_init()
tar_git_status()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_data}()},
\code{\link{tar_git_status_targets}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_git_status_data.R
\name{tar_git_status_data}
\alias{tar_git_status_data}
\title{Status of the data repository (Git)}
\usage{
tar_git_status_data(
  store = targets::tar_config_get("store"),
  stash_gitignore = TRUE
)
}
\arguments{
\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[targets:tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{stash_gitignore}{Logical of length 1, whether to temporarily
stash the \code{.gitignore} file of the data store. See the
"Stashing .gitignore" section for details.}
}
\value{
If the data repository exists, the return value is the data frame
produced by \code{gert::git_status(repo = store)}. If the data store has no Git
repository, then the return value is \code{NULL}.
}
\description{
Show the Git status of the data repository.
}
\section{Stashing .gitignore}{

The \code{targets} package writes a \code{.gitignore} file to new data stores
in order to prevent accidental commits to the code Git repository.
Unfortunately, for \code{gittargets}, this automatic \code{.gitignore} file
interferes with proper data versioning. So by default, \code{gittargets}
temporarily stashes it to a hidden file called \code{.gittargets_gitignore}
inside the data store. If your R program crashes while the stash
is active, you can simply move it manually back to \code{.gitignore}
or run \code{tar_git_status_data()} to restore the stash automatically
if no \code{.gitignore} already exists.
}

\examples{
if (Sys.getenv("TAR_EXAMPLES") == "true" && tar_git_ok(verbose = FALSE)) {
targets::tar_dir({ # Containing code does not modify the user's file space.
targets::tar_script(tar_target(data, 1))
targets::tar_make()
list.files("_targets", all.files = TRUE)
gert::git_init()
tar_git_init()
tar_git_status_data()
})
}
}
\seealso{
Other git: 
\code{\link{tar_git_checkout}()},
\code{\link{tar_git_init}()},
\code{\link{tar_git_log}()},
\code{\link{tar_git_ok}()},
\code{\link{tar_git_snapshot}()},
\code{\link{tar_git_status_code}()},
\code{\link{tar_git_status_targets}()},
\code{\link{tar_git_status}()}
}
\concept{git}
