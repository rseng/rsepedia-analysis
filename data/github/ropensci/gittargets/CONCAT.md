
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
