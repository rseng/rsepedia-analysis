
# stantargets <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/430_status.svg)](https://github.com/ropensci/software-review/issues/430)
[![joss](https://joss.theoj.org/papers/10.21105/joss.03193/status.svg)](https://doi.org/10.21105/joss.03193)
[![zenodo](https://zenodo.org/badge/315447649.svg)](https://zenodo.org/badge/latestdoi/315447649)
[![R
Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
<!--
[![cran](http://www.r-pkg.org/badges/version/stantargets)](https://cran.r-project.org/package=stantargets)
-->
[![active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/stantargets/workflows/check/badge.svg)](https://github.com/ropensci/stantargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/stantargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/stantargets)
[![lint](https://github.com/ropensci/stantargets/workflows/lint/badge.svg)](https://github.com/ropensci/stantargets/actions?query=workflow%3Alint)

Bayesian data analysis usually incurs long runtimes and cumbersome
custom code, and the process of prototyping and deploying custom
[Stan](https://mc-stan.org) models can become a daunting software
engineering challenge. To ease this burden, the `stantargets` R package
creates [Stan](https://mc-stan.org) pipelines that are concise,
efficient, scalable, and tailored to the needs of Bayesian
statisticians. Leveraging
[`targets`](https://docs.ropensci.org/targets/), `stantargets` pipelines
automatically parallelize the computation and skip expensive steps when
the results are already up to date. Minimal custom user-side code is
required, and there is no need to manually configure branching, so
`stantargets` is easier to use than
[`targets`](https://docs.ropensci.org/targets/) and
[`CmdStanR`](https://mc-stan.org/cmdstanr/) directly. `stantargets` can
access all of [`cmdstanr`](https://github.com/stan-dev/cmdstanr)’s major
algorithms (MCMC, variational Bayes, and optimization) and it supports
both single-fit workflows and multi-rep simulation studies.

## Prerequisites

1.  The [prerequisites of the `targets` R
    package](https://docs.ropensci.org/targets/index.html#prerequisites).
2.  Basic familiarity with
    [`targets`](https://docs.ropensci.org/targets/): watch minutes 7
    through 40 of [this video](https://youtu.be/Gqn7Xn4d5NI?t=439), then
    read [this
    chapter](https://books.ropensci.org/targets/walkthrough.html) of the
    [user manual](https://books.ropensci.org/targets/).
3.  Familiarity with Bayesian Statistics and
    [Stan](https://mc-stan.org/). Prior knowledge of
    [`cmdstanr`](https://mc-stan.org/cmdstanr/) helps.

## How to get started

Read the `stantargets`
[introduction](https://docs.ropensci.org/stantargets/articles/introduction.html)
and
[simulation](https://docs.ropensci.org/stantargets/articles/simulation.html)
vignettes, and use <https://docs.ropensci.org/stantargets/> as a
reference while constructing your own workflows. Visit
<https://github.com/wlandau/stantargets-example-validation> for an
example project based on the [simulation
vignette](https://docs.ropensci.org/stantargets/articles/simulation.html).
The example has an [RStudio Cloud
workspace](https://rstudio.cloud/project/2466069) which allows you to
run the project in a web browser.

## Example projects

| Description                                                                                                        | Link                                                |
| ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------- |
| Validating a minimal Stan model                                                                                    | <https://github.com/wlandau/targets-stan>           |
| Using Target Markdown and `stantargets` to validate a Bayesian longitudinal model for clinical trial data analysis | <https://github.com/wlandau/rmedicine2021-pipeline> |

## Installation

Install the GitHub development version to access the latest features and
patches.

``` r
remotes::install_github("ropensci/stantargets")
```

The [CmdStan](https://github.com/stan-dev/cmdstan) command line
interface is also required.

``` r
cmdstanr::install_cmdstan()
```

If you have problems installing
[CmdStan](https://github.com/stan-dev/cmdstan), please consult the
[installation guide of
`cmdstanr`](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) and the
[installation guide of
CmdStan](https://mc-stan.org/docs/2_26/cmdstan-guide/cmdstan-installation.html).
Alternatively, the [Stan discourse](https://discourse.mc-stan.org) is a
friendly place to ask Stan experts for help.

## Usage

First, write a [`_targets.R`
file](https://books.ropensci.org/targets/walkthrough.html) that loads
your packages, defines a function to generate
[Stan](https://mc-stan.org/) data, and lists a pipeline of targets. The
target list can call target factories like
[`tar_stan_mcmc()`](https://docs.ropensci.org/stantargets/reference/tar_stan_mcmc.html)
as well as ordinary targets with
[`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html).
The following minimal example is simple enough to contain entirely
within the `_targets.R` file, but for larger projects, you may wish to
store functions in separate files as in the
[`targets-stan`](https://github.com/wlandau/targets-stan) example.

``` r
# _targets.R
library(targets)
library(stantargets)

generate_data <- function() {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_stan_mcmc(
    name = example,
    stan_files = "x.stan",
    data = generate_data()
  )
)
```

Run
[`tar_visnetwork()`](https://docs.ropensci.org/targets/reference/tar_visnetwork.html)
to check `_targets.R` for correctness, then call
[`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html)
to run the pipeline. Access the results using
[`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html),
e.g. `tar_read(example_summary_x)`. Visit the [introductory
vignette](https://docs.ropensci.org/stantargets/articles/introduction.html)
to read more about this example.

## How it works behind the scenes

`stantargets` supports specialized [target
factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories)
that create ensembles of [target
objects](https://docs.ropensci.org/targets/reference/tar_target.html)
for [`cmdstanr`](https://github.com/stan-dev/cmdstanr) workflows. These
[target
factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories)
abstract away the details of
[`targets`](https://docs.ropensci.org/targets/) and
[`cmdstanr`](https://github.com/stan-dev/cmdstanr) and make both
packages easier to use. For details, please read the [introductory
vignette](https://docs.ropensci.org/stantargets/articles/introduction.html).

## Help

If you have trouble using `stantargets`, you can ask for help in the
[GitHub discussions
forum](https://github.com/ropensci/stantargets/discussions/categories/help).
Because the purpose of `stantargets` is to combine
[`targets`](https://docs.ropensci.org/targets/) and
[`cmdstanr`](https://github.com/stan-dev/cmdstanr), your issue may have
something to do with one of the latter two packages, a [dependency of
`targets`](https://github.com/ropensci/targets/blob/4e3ef2a6c986f558a25e544416f480fc01236b6b/DESCRIPTION#L49-L88),
or [Stan](https://mc-stan.org) itself. When you troubleshoot, peel back
as many layers as possible to isolate the problem. For example, if the
issue comes from [`cmdstanr`](https://github.com/stan-dev/cmdstanr),
create a [reproducible example](https://reprex.tidyverse.org) that
directly invokes [`cmdstanr`](https://github.com/stan-dev/cmdstanr)
without invoking `stantargets`. The GitHub discussion and issue forums
of those packages, as well as the [Stan
discourse](https://discourse.mc-stan.org), are great resources.

## Participation

Development is a community effort, and we welcome discussion and
contribution. Please note that this package is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Citation

``` r
citation("stantargets")
#> 
#> To cite stantargets in publications use:
#> 
#>   Landau, W. M., (2021). The stantargets R package: a workflow
#>   framework for efficient reproducible Stan-powered Bayesian data
#>   analysis pipelines. Journal of Open Source Software, 6(60), 3193,
#>   https://doi.org/10.21105/joss.03193
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {The stantargets {R} package: a workflow framework for efficient reproducible {S}tan-powered {B}ayesian data analysis pipelines},
#>     author = {William Michael Landau},
#>     journal = {Journal of Open Source Software},
#>     year = {2021},
#>     volume = {6},
#>     number = {60},
#>     pages = {3193},
#>     url = {https://doi.org/10.21105/joss.03193},
#>   }
```
# stantargets 0.0.3.9001


# stantargets 0.0.3

* Update docs to changes in `cmdstanr`, `posterior`, and `targets`.

# stantargets 0.0.2

* Reference JOSS paper.

# stantargets 0.0.1

* Skip tests if CmdStan is not installed (@sakrejda).
* Use custom `generate_data()` function in the docs, as opposed to `tar_stan_example_data()` directly (@sakrejda).
* Add the `pedantic` argument for compilation (@sakrejda).
* Reduce dependencies on some `rlang` functions like `sym()` (@sakrejda).
* Change `trn()` to `if_any()` (@mattwarkentin, @sakrejda, @tjmahr).
* Add @sakrejda and @mattwarkentin as reviewers in the `DESCRIPTION`.
* Talk about the R package and system dependencies of `stantargets` in the README (@mattwarkentin).
* Throw an error earlier if the Stan file does not exist (@sakrejda, @mattwarkentin)
* Use `@format` `roxygen2` tag for data generation (@mattwarkentin).
* Use `@family` go cross-reference functions (@mattwarkentin).
* Elaborate on the roles and return values of specific targets generated by target factories (@mattwarkentin).
* Undergo rOpenSci peer review and transition to rOpenSci.
* Link to an example project.

# stantargets 0.0.0.9002

* Return the executable file after the Stan source file in model compilation targets.
* Replace the `log` argument with `stdout` and `stderr` (#23).
* Switch meaning of `%||%` and `%|||%` to conform to historical precedent.

# stantargets 0.0.0.9001

* Join on data to summary output using `.join_data` in the Stan data (#18).
* Pre-compile models for testing and add an environment variable to skip tests that always force recompilation (#19).
* Load packages for any target computing summaries.

# stantargets 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributing

Development is a community effort, and we welcome participation.

## Code of Conduct

By participating in this project, you agree to abide by the [code of conduct](https://ropensci.org/code-of-conduct/).

## Issues

Anyone can start or contribute to an [issue](https://github.com/ropensci/stantargets/issues) or [discussion thread](https://github.com/ropensci/stantargets/discussions). Issues are mainly for bug reports and package maintenance, and discussions are for usage help and brainstorming. Please respect the following guidelines.

* Before posting a new issue or discussion, please take a moment to search for existing threads in order to avoid duplication.
* For bug reports: if you can, please install the latest GitHub version of `stantargets` (i.e. `remotes::install_github("ropensci/stantargets")`) and verify that the issue still persists.
* Describe your issue in prose as clearly and concisely as possible.
* For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Development

External code contributions are extremely helpful in the right circumstances. Here are the recommended steps.

1. Prior to contribution, please propose your idea in a [new issue thread](https://github.com/ropensci/stantargets/issues) so you and the maintainer can define the intent and scope of your work.
2. [Fork the repository](https://help.github.com/articles/fork-a-repo/).
3. Follow the [GitHub flow](https://guides.github.com/introduction/flow/index.html) to create a new branch, add commits, and open a pull request.
4. Discuss your code with the maintainer in the pull request thread.
5. If everything looks good, the maintainer will merge your code into the project.

Please also follow these additional guidelines.

* Respect the architecture and reasoning of the package. Depending on the scope of your work, you may want to read the design documents (package vignettes).
* If possible, keep contributions small enough to easily review manually. It is okay to split up your work into multiple pull requests.
* Format your code according to the [tidyverse style guide](https://style.tidyverse.org/) and check your formatting with the `lint_package()` function from the [`lintr`](https://github.com/jimhester/lintr) package.
* Check code coverage with `covr::package_coverage()`. Automated tests should cover all the new or changed functionality in your pull request.
* Run overall package checks with `devtools::check()` and `goodpractice::gp()`
* Describe your contribution in the project's [`NEWS.md`](https://github.com/ropensci/stantargets/blob/main/NEWS.md) file. Be sure to mention relevant GitHub issue numbers and your GitHub name as done in existing news entries.
* If you feel contribution is substantial enough for official author or contributor status, please add yourself to the `Authors@R` field of the [`DESCRIPTION`](https://github.com/ropensci/stantargets/blob/main/DESCRIPTION) file.
# Prework

* [ ] I understand and agree to this repository's [code of conduct](https://ropensci.org/code-of-conduct/).
* [ ] I understand and agree to this repository's [contributing guidelines](https://github.com/ropensci/stantargets/blob/main/CONTRIBUTING.md).
* [ ] I have already submitted an issue to the [issue tracker](http://github.com/ropensci/stantargets/issues) to discuss my idea with the maintainer.

# Related GitHub issues and pull requests

* Ref: #

# Summary

Please explain the purpose and scope of your contribution.
---
name: Maintenance
about: "Something in targets needs work: updates, documentation, etc. Not a bug, performance issue, or new feature."
title: ""
labels: "type: new maintenance"
assignees: ""
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/targets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
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
about: Something is wrong with stantargets.
title: ""
labels: "type: bug"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/stantargets/blob/main/CONTRIBUTING.md).
* [ ] Confirm that your issue is most likely a genuine bug in `stantargets` and not a known limitation, usage error, or bug in another package that `stantargets` depends on.
* [ ] If there is [already a relevant issue](https://github.com/ropensci/stantargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the bug.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Expected result

What should have happened? Please be as specific as possible.

## Diagnostic information

* A [reproducible example](https://github.com/tidyverse/reprex).
* Session info, available through `sessionInfo()` or [`reprex(si = TRUE)`](https://github.com/tidyverse/reprex).
* A stack trace from `traceback()` or `rlang::trace_back()`.
* The [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity) of the GitHub commit of `stantargets` currently installed. `packageDescription("stantargets")$GithubSHA1` shows you this.
---
name: New feature
about: Suggest a new feature.
title: ""
labels: "type: new feature"
assignees: wlandau
---

## Prework

* [ ] I understand and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and  [contributing guidelines](https://github.com/ropensci/stantargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/stantargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] New features take time and effort to create, and they take even more effort to maintain. So if the purpose of the feature is to resolve a struggle you are encountering personally, please consider first posting a "trouble" or "other" issue so we can discuss your use case and search for existing solutions first.

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

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/stantargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/stantargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
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

How poorly does `stantargets` perform? To find out, we recommend you use the [`proffer`](https://github.com/r-prof/proffer) package and take screenshots of the results displayed in your browser.

```r
library(stantargets)
library(proffer)
px <- pprof({
  # All your stantargets code goes here.
})
```
---
title: 'The stantargets R package: a workflow framework for efficient reproducible Stan-powered Bayesian data analysis pipelines'
tags:
- R
- reproducibility
- high-performance computing
- pipeline
- workflow
- Make
- Bayesian
- Stan
date: "2020"
output: pdf_document
authors:
- name: William Michael Landau
  orcid: 0000-0003-1878-3253
  email: will.landau@gmail.com
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Eli Lilly and Company
  index: 1
---

# Statement of need

Real-world Bayesian data analysis is a thorough process of investigation and experimentation. Statisticians iteratively refine and compare multiple models to improve inference and understand model behavior [@bayesworkflow]. Intermediate results inform downstream model-building decisions, and the final models do not always agree with the initial choices. A typical step of this empirical development process requires a computational method such as Markov chain Monte Carlo to approximate the posterior distribution of the model parameters given the data [@bda3]. Even with fast and flexible probabilistic programming languages like Stan [@stan], this computation can take several minutes or hours to complete, and successive iterations become expensive enough to obstruct the research.

# Summary

The [`stantargets`](https://github.com/ropensci/stantargets) R package [@stantargets] reduces the practical burdens of developing and maintaining Bayesian data analysis workflows with Stan. It expresses the models, datasets, and inferential results as interdependent components of a formal pipeline, tracks these components for changes, and automatically reruns the affected steps in response to these changes, optionally with distributed computing on a cluster. If a step is already up to date with its upstream dependencies, [`stantargets`](https://github.com/ropensci/stantargets) automatically skips it, potentially saving hours of runtime. When the whole pipeline is up to date, the user has tangible evidence that the output matches the underlying code and data, which affirms reproducibility.

The [`stantargets`](https://github.com/ropensci/stantargets) package is an extension of [`cmdstanr`](https://github.com/stan-dev/cmdstanr) [@cmdstanr], a lightweight interface to Stan, and [`targets`](https://github.com/ropensci/targets) [@targets], a general-purpose pipeline toolkit for reproducible research and high-performance computing. [`stantargets`](https://github.com/ropensci/stantargets) builds [`targets`](https://github.com/ropensci/targets)-powered pipelines specifically tailored to Bayesian statistics with [`cmdstanr`](https://github.com/stan-dev/cmdstanr), from single-run workflows to large-scale simulation studies. Using domain knowledge, [`stantargets`](https://github.com/ropensci/stantargets) abstracts away burdensome low-level configuration details and streamlines pipeline construction, freeing Bayesian statisticians to focus less on software development and more on model development. [`stantargets`](https://github.com/ropensci/stantargets) is part of the [R Targetopia](https://wlandau.github.io/targetopia) [@targetopia], an emerging ecosystem of R packages to democratize reproducible analysis pipelines across multiple domains of Statistics and data science.

# References
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

# stantargets <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/430_status.svg)](https://github.com/ropensci/software-review/issues/430)
[![joss](https://joss.theoj.org/papers/10.21105/joss.03193/status.svg)](https://doi.org/10.21105/joss.03193)
[![zenodo](https://zenodo.org/badge/315447649.svg)](https://zenodo.org/badge/latestdoi/315447649)
[![R Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
<!--
[![cran](http://www.r-pkg.org/badges/version/stantargets)](https://cran.r-project.org/package=stantargets)
-->
[![active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/stantargets/workflows/check/badge.svg)](https://github.com/ropensci/stantargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/stantargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/stantargets)
[![lint](https://github.com/ropensci/stantargets/workflows/lint/badge.svg)](https://github.com/ropensci/stantargets/actions?query=workflow%3Alint)


Bayesian data analysis usually incurs long runtimes and cumbersome custom code, and the process of prototyping and deploying custom [Stan](https://mc-stan.org) models can become a daunting software engineering challenge. To ease this burden, the `stantargets` R package creates [Stan](https://mc-stan.org) pipelines that are concise, efficient, scalable, and tailored to the needs of Bayesian statisticians. Leveraging [`targets`](https://docs.ropensci.org/targets/), `stantargets` pipelines automatically parallelize the computation and skip expensive steps when the results are already up to date. Minimal custom user-side code is required, and there is no need to manually configure branching, so `stantargets` is easier to use than [`targets`](https://docs.ropensci.org/targets/) and [`CmdStanR`](https://mc-stan.org/cmdstanr/) directly. `stantargets` can access all of [`cmdstanr`](https://github.com/stan-dev/cmdstanr)'s major algorithms (MCMC, variational Bayes, and optimization) and it supports both single-fit workflows and multi-rep simulation studies.

## Prerequisites

1. The [prerequisites of the `targets` R package](https://docs.ropensci.org/targets/index.html#prerequisites).
1. Basic familiarity with [`targets`](https://docs.ropensci.org/targets/): watch minutes 7 through 40 of [this video](https://youtu.be/Gqn7Xn4d5NI?t=439), then read [this chapter](https://books.ropensci.org/targets/walkthrough.html) of the [user manual](https://books.ropensci.org/targets/).
1. Familiarity with Bayesian Statistics and [Stan](https://mc-stan.org/). Prior knowledge of [`cmdstanr`](https://mc-stan.org/cmdstanr/) helps.

## How to get started

Read the `stantargets` [introduction](https://docs.ropensci.org/stantargets/articles/introduction.html) and [simulation](https://docs.ropensci.org/stantargets/articles/simulation.html) vignettes, and use <https://docs.ropensci.org/stantargets/> as a reference while constructing your own workflows. Visit <https://github.com/wlandau/stantargets-example-validation> for an example project based on the [simulation vignette](https://docs.ropensci.org/stantargets/articles/simulation.html). The example has an [RStudio Cloud workspace](https://rstudio.cloud/project/2466069) which allows you to run the project in a web browser.

## Example projects

Description | Link
---|---
Validating a minimal Stan model | <https://github.com/wlandau/targets-stan>
Using Target Markdown and `stantargets` to validate a Bayesian longitudinal model for clinical trial data analysis | <https://github.com/wlandau/rmedicine2021-pipeline>

## Installation

Install the GitHub development version to access the latest features and patches.

```{r, eval = FALSE}
remotes::install_github("ropensci/stantargets")
```

The [CmdStan](https://github.com/stan-dev/cmdstan) command line interface is also required.

```{r, eval = FALSE}
cmdstanr::install_cmdstan()
```

If you have problems installing [CmdStan](https://github.com/stan-dev/cmdstan), please consult the [installation guide of `cmdstanr`](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) and the [installation guide of CmdStan](https://mc-stan.org/docs/2_26/cmdstan-guide/cmdstan-installation.html). Alternatively, the [Stan discourse](https://discourse.mc-stan.org) is a friendly place to ask Stan experts for help.

## Usage

First, write a [`_targets.R` file](https://books.ropensci.org/targets/walkthrough.html) that loads your packages, defines a function to generate [Stan](https://mc-stan.org/) data, and lists a pipeline of targets. The target list can call target factories like [`tar_stan_mcmc()`](https://docs.ropensci.org/stantargets/reference/tar_stan_mcmc.html) as well as ordinary targets with [`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html). The following minimal example is simple enough to contain entirely within the `_targets.R` file, but for larger projects, you may wish to store functions in separate files as in the [`targets-stan`](https://github.com/wlandau/targets-stan) example.

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

generate_data <- function() {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_stan_mcmc(
    name = example,
    stan_files = "x.stan",
    data = generate_data()
  )
)
```

Run [`tar_visnetwork()`](https://docs.ropensci.org/targets/reference/tar_visnetwork.html) to check `_targets.R` for correctness, then call [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html) to run the pipeline. Access the results using [`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html), e.g. `tar_read(example_summary_x)`. Visit the [introductory vignette](https://docs.ropensci.org/stantargets/articles/introduction.html) to read more about this example.

## How it works behind the scenes

`stantargets` supports specialized [target factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories) that create ensembles of [target objects](https://docs.ropensci.org/targets/reference/tar_target.html) for [`cmdstanr`](https://github.com/stan-dev/cmdstanr) workflows. These [target factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories) abstract away the details of [`targets`](https://docs.ropensci.org/targets/) and [`cmdstanr`](https://github.com/stan-dev/cmdstanr) and make both packages easier to use. For details, please read the [introductory vignette](https://docs.ropensci.org/stantargets/articles/introduction.html).

## Help

If you have trouble using `stantargets`, you can ask for help in the [GitHub discussions forum](https://github.com/ropensci/stantargets/discussions/categories/help). Because the purpose of `stantargets` is to combine [`targets`](https://docs.ropensci.org/targets/) and [`cmdstanr`](https://github.com/stan-dev/cmdstanr), your issue may have something to do with one of the latter two packages, a [dependency of `targets`](https://github.com/ropensci/targets/blob/4e3ef2a6c986f558a25e544416f480fc01236b6b/DESCRIPTION#L49-L88), or [Stan](https://mc-stan.org) itself. When you troubleshoot, peel back as many layers as possible to isolate the problem. For example, if the issue comes from [`cmdstanr`](https://github.com/stan-dev/cmdstanr), create a [reproducible example](https://reprex.tidyverse.org) that directly invokes [`cmdstanr`](https://github.com/stan-dev/cmdstanr) without invoking `stantargets`. The GitHub discussion and issue forums of those packages, as well as the [Stan discourse](https://discourse.mc-stan.org), are great resources.

## Participation

Development is a community effort, and we welcome discussion and contribution. Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Citation

```{r, warning = FALSE}
citation("stantargets")
```
---
title: "Introduction to stantargets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to stantargets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
# With the root.dir option below,
# this vignette runs the R code in a temporary directory
# so new files are written to temporary storage
# and not the user's file space.
knitr::opts_knit$set(root.dir = fs::dir_create(tempfile()))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
if (identical(Sys.getenv("NOT_CRAN", unset = "false"), "false")) {
  knitr::opts_chunk$set(eval = FALSE)
}
library(cmdstanr)
library(dplyr)
library(targets)
library(stantargets)
if (identical(Sys.getenv("IN_PKGDOWN"), "true")) {
  cmdstanr::install_cmdstan()
}
```

The `stantargets` package makes it easy to run a single Stan model and keep track of the results. [`cmdstanr`](https://github.com/stan-dev/cmdstanr) fits the models, and [`targets`](https://docs.ropensci.org/targets/) manages the workflow and helps avoid unnecessary computation.

First, write a Stan model file.

```{r}
lines <- "data {
  int <lower = 1> n;
  vector[n] x;
  vector[n] y;
  real true_beta;
}
parameters {
  real beta;
}
model {
  y ~ normal(x * beta, 1);
  beta ~ normal(0, 1);
}"
writeLines(lines, "x.stan")
```

A typical workflow proceeds as follows:

1. Prepare a list of input data to Stan, including vector elements `x` and `y`.
1. Fit the Stan model using the list of input data.
1. Use the fitted model object to compute posterior summaries and convergence diagnostics.
1. Use the fitted model object to extract posterior draws of parameters and store them in a tidy data frame.
1. Use the fitted model to compute Hamiltonian Monte Carlo (HMC) diagnostics.

`stantargets` expresses this workflow using the [`tar_stan_mcmc()`](https://docs.ropensci.org/stantargets/reference/tar_stan_mcmc.html) function. To use it in a [`targets`](https://docs.ropensci.org/targets/) pipeline, invoke it from the `_targets.R` script of the project. 

```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function(n = 10) {
    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, x * true_beta, 1)
    list(n = n, x = x, y = y, true_beta = true_beta)
  }
  list(
    tar_stan_mcmc(
      example,
      "x.stan",
      generate_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

# The _targets.R file ends with a list of target objects
# produced by stantargets::tar_stan_mcmc(), targets::tar_target(), or similar.
list(
  tar_stan_mcmc(
    example,
    "x.stan",
    generate_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
```

Above, `tar_stan_mcmc(example, ...)` only *defines* the pipeline. It does not actually run Stan, it declares the targets that will eventually run Stan. Run `tar_manifest()` to show specific details about the targets.

```{r}
tar_manifest()
```

Each target listed above is responsible for a piece of the workflow.

* `example_file_x`: Reproducibly track changes to the Stan model file.
* `example_data`: Run the code you supplied to the `data` argument of `tar_stan_mcmc()` and return a dataset compatible with Stan.
* `example_mcmc_x`: Run the MCMC and return an object of class `CmdStanMCMC`.
* `example_draws_X`: Return a friendly `tibble` of the posterior draws from `example`. Uses the `$draws()` method. Suppress with `draws = FALSE` in `tar_stan_mcmc()`.
* `example_summaries_x`: Return a friendly `tibble` of the posterior summaries from `example`. Uses the `$summary()` method. Suppress with `summary = FALSE` in `tar_stan_mcmc()`.
* `example_diagnostics_x`: Return a friendly `tibble` of the sampler diagnostics from `example`. Uses the `$sampler_diagnostics()` method. Suppress with `diagnostics = FALSE` in `tar_stan_mcmc()`.

The suffix `_x` comes from the base name of the model file, in this case `x.stan`. If you supply multiple model files to the `stan_files` argument, all the models share the same dataset, and the suffixes distinguish among the various targets.

[`targets`](https://docs.ropensci.org/targets/) produces a graph to show the dependency relationships among the targets. Below, the MCMC depends on the model file and the data, and the draws, summary, and diagnostics depend on the MCMC.

```{r, output = FALSE, message = FALSE}
tar_visnetwork(targets_only = TRUE)
```

Run the computation with `tar_make()`.

```{r, output = FALSE}
tar_make()
```

The output lives in a special folder called `_targets/` and you can retrieve it with functions `tar_load()` and `tar_read()` (from [`targets`](https://docs.ropensci.org/targets/)).

```{r}
tar_read(example_summary_x)
```

At this point, all our results are up to date because their dependencies did not change.

```{r}
tar_make()
```

But if we change the underlying code or data, some of the targets will no longer be valid, and they will rerun during the next `tar_make()`. Below, we change the Stan model file, so the MCMC reruns while the data is skipped. This behavior saves time and enhances reproducibility.

```{r}
write(" ", file = "x.stan", append = TRUE)
```

```{r}
tar_outdated()
```

```{r}
tar_visnetwork(targets_only = TRUE)
```

```{r, output = FALSE}
tar_make()
```

At this point, we can add more targets and custom functions for additional post-processing.

```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function(n = 10) {
    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, x * true_beta, 1)
    list(n = n, x = x, y = y, true_beta = true_beta)
  }
  list(
    tar_stan_mcmc(
      example,
      "x.stan",
      generate_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    ),
    tar_stan_summary(
      custom_summary,
      fit = example_mcmc_x,
      summaries = list(~posterior::quantile2(.x, probs = c(0.25, 0.75)))
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_stan_mcmc(
    example,
    "x.stan",
    generate_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_stan_summary(
    custom_summary,
    fit = example_mcmc_x,
    summaries = list(~posterior::quantile2(.x, probs = c(0.25, 0.75)))
  )
)
```

In the graph, our new `custom_summary` target should be connected to the upstream `example` target, and only `custom_summary` should be out of date.

```{r}
tar_visnetwork(targets_only = TRUE)
```

In the next `tar_make()`, we skip the expensive MCMC and just run the custom summary.

```{r, output = FALSE, warning = FALSE}
tar_make()
```

```{r}
tar_read(custom_summary)
```

## Multiple models

`tar_stan_mcmc()` and related functions allow you to supply multiple models to `stan_files`. If you do, each model will run on the same dataset. Consider a new model `y.stan`.

```{r}
lines <- "data {
  int <lower = 1> n;
  vector[n] x;
  vector[n] y;
  real true_beta;
}
parameters {
  real beta;
}
model {
  y ~ normal(x * x * beta, 1); // Regress on x^2 instead of x.
  beta ~ normal(0, 1);
}"
writeLines(lines, "y.stan")
```

To include this `y.stan`, we add it to the `stan_files` argument of `tar_stan_mcmc()`.

```{r, echo = FALSE}
library(targets)
file.copy("x.stan", "y.stan")
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function(n = 10) {
    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, x * true_beta, 1)
    list(n = n, x = x, y = y, true_beta = true_beta)
  }
  list(
    tar_stan_mcmc(
      example,
      c("x.stan", "y.stan"), # another model
      generate_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    ),
    tar_stan_summary(
      custom_summary,
      fit = example_mcmc_x,
      summaries = list(~posterior::quantile2(.x, probs = c(0.25, 0.75)))
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_stan_mcmc(
    example,
    c("x.stan", "y.stan"), # another model
    generate_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_stan_summary(
    custom_summary,
    fit = example_mcmc_x,
    summaries = list(~posterior::quantile2(.x, probs = c(0.25, 0.75)))
  )
)
```

In the graph below, notice how the `*_x` targets and `*_y` targets are both connected to `example_data` upstream.

```{r}
tar_visnetwork(targets_only = TRUE)
```

## Generated quantities

It is possible to use the `CmdStanMCMC` object from one run to simulate generated quantities downstream. For example, the `tar_stan_gq_rep_summaries()` function takes a single `CmdStanMCMC` object, produces multiple replications of generated quantities from multiple models, and aggregates the summaries from each. The following pipeline uses this technique to repeatedly draw from the posterior predictive distribution.

```{r}
lines <- "data {
  int <lower = 1> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real beta;
}
model {
  y ~ normal(x * beta, 1);
  beta ~ normal(0, 1);
}
generated quantities {
  real y_rep[n] = normal_rng(x * beta, 1); // posterior predictive draws
}"
writeLines(lines, "gen.stan")
```

```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function(n = 10) {
    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, x * true_beta, 1)
    list(n = n, x = x, y = y, true_beta = true_beta)
  }
  list(
    tar_stan_mcmc(
      example,
      "x.stan",
      generate_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    ),
    tar_stan_gq_rep_summary(
      postpred,
      stan_files = "gen.stan",
      fitted_params = example_mcmc_x, # one CmdStanFit object
      data = generate_data(), # multiple simulated datasets
      batches = 2, # 2 dynamic branches
      reps = 5, # 5 replications per branch
      quiet = TRUE,
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_stan_mcmc(
    example,
    "x.stan",
    generate_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_stan_gq_rep_summary(
    postpred,
    stan_files = "gen.stan",
    fitted_params = example_mcmc_x, # one CmdStanFit object
    data = generate_data(), # Function runs once per rep.
    batches = 2, # 2 dynamic branches
    reps = 5, # 5 replications per branch
    quiet = TRUE,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
```

Since we have defined many objects in the pipeline, it is extra important to check the graph to be sure everything is connected.

```{r}
tar_visnetwork(targets_only = TRUE)
```

Then, we run the computation. The original MCMC is already up to date, so we only run the targets relevant to the generated quantities.

```{r}
tar_make()
```

Finally, we read the summaries of posterior predictive samples.

```{r}
tar_read(postpred)
```

## More information

For more on [`targets`](https://docs.ropensci.org/targets/), please visit the reference website <https://docs.ropensci.org/targets/> or the user manual <https://books.ropensci.org/targets/>. The manual walks though advanced features of `targets` such as [high-performance computing](https://books.ropensci.org/targets/hpc.html) and [cloud storage support](https://books.ropensci.org/targets/cloud.html).
---
title: "Bayesian simulation pipelines with stantargets"
output: rmarkdown::html_vignette
bibliography: simulation.bib
vignette: >
  %\VignetteIndexEntry{Bayesian simulation pipelines with stantargets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
# With the root.dir option below,
# this vignette runs the R code in a temporary directory
# so new files are written to temporary storage
# and not the user's file space.
knitr::opts_knit$set(root.dir = fs::dir_create(tempfile()))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
if (identical(Sys.getenv("NOT_CRAN", unset = "false"), "false")) {
  knitr::opts_chunk$set(eval = FALSE)
}
library(cmdstanr)
library(dplyr)
library(targets)
library(stantargets)
if (identical(Sys.getenv("IN_PKGDOWN"), "true")) {
  cmdstanr::install_cmdstan()
}
```

## Background

The [introductory vignette](https://docs.ropensci.org/stantargets/articles/introduction.html) vignette caters to Bayesian data analysis workflows with few datasets to analyze. However, it is sometimes desirable to run one or more Bayesian models repeatedly across multiple simulated datasets. Examples:

1. Validate the implementation of a Bayesian model using simulation.
2. Simulate a randomized controlled experiment to explore frequentist properties such as power and Type I error.

This vignette focuses on (1).

## Example project

Visit <https://github.com/wlandau/stantargets-example-validation> for an example project based on this vignette. The example has an [RStudio Cloud workspace](https://rstudio.cloud/project/2466069) which allows you to run the project in a web browser.

## Simulation-based validation study

This particular example uses the concept of calibration that Bob Carpenter [explains here](https://statmodeling.stat.columbia.edu/2017/04/12/bayesian-posteriors-calibrated/) [@carpenter2017]. The goal is to simulate multiple datasets from the model below, analyze each dataset, and assess how often the estimated posterior intervals cover the true parameters from the prior predictive simulations. If coverage is no systematically different from nominal, this is evidence that the model was implemented correctly. The quantile method by @cook2006 generalizes this concept, and simulation-based calibration [@talts2020] generalizes further. The interval-based technique featured in this vignette is not as robust as SBC, but it may be more expedient for large models because it does not require visual inspection of multiple histograms.

```{r}
lines <- "data {
  int <lower = 1> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] beta;
}
model {
  y ~ normal(beta[1] + x * beta[2], 1);
  beta ~ normal(0, 1);
}"
writeLines(lines, "model.stan")
```

Next, we define a pipeline to simulate multiple datasets and fit each dataset with the model. In our data-generating function, we put the true parameter values of each simulation in a special `.join_data` list. `stantargets` will automatically join the elements of `.join_data` to the correspondingly named variables in the summary output. This will make it super easy to check how often our posterior intervals capture the truth. As for scale, generate 10 datasets (5 batches with 2 replications each) and run the model on each of the 10 datasets.^[Internally, each batch is a [dynamic branch target](https://books.ropensci.org/targets/dynamic.html), and the number of replications determines the amount of work done within a branch. In the general case, [batching](https://books.ropensci.org/targets/dynamic.html#batching) is a way to find the right compromise between target-specific overhead and the horizontal scale of the pipeline.] By default, each of the 10 model runs computes 4 MCMC chains with 2000 MCMC iterations each (including burn-in) and you can adjust with the `chains`, `iter_sampling`, and `iter_warmup` arguments of `tar_stan_mcmc_rep_summary()`.

```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  # Use computer memory more sparingly:
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  simulate_data <- function(n = 10L) {
    beta <- rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- rnorm(n, beta[1] + x * beta[2], 1)
    list(
      n = n,
      x = x,
      y = y,
      .join_data = list(beta = beta)
    )
  }
  list(
    tar_stan_mcmc_rep_summary(
      model,
      "model.stan",
      simulate_data(),
      batches = 5, # Number of branch targets.
      reps = 2, # Number of model reps per branch target.
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      ),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)
options(crayon.enabled = FALSE)
# Use computer memory more sparingly:
tar_option_set(memory = "transient", garbage_collection = TRUE)

simulate_data <- function(n = 10L) {
  beta <- rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- rnorm(n, beta[1] + x * beta[2], 1)
  list(
    n = n,
    x = x,
    y = y,
    .join_data = list(beta = beta)
  )
}

list(
  tar_stan_mcmc_rep_summary(
    model,
    "model.stan",
    simulate_data(), # Runs once per rep.
    batches = 5, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    ),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
```

We now have a pipeline that runs the model 10 times: 5 batches (branch targets) with 2 replications per batch.

```{r}
tar_visnetwork()
```

Run the computation with `tar_make()`

```{r, output = FALSE, warning = FALSE}
tar_make()
```

The result is an aggregated data frame of summary statistics, where the `.rep` column distinguishes among individual replicates. We have the posterior intervals for `beta` in columns `q2.5` and `q97.5`. And thanks to `.join_data` in `simulate_data()`, there is a special `.join_data` column in the output to indicate the true value of each parameter from the simulation.

```{r}
tar_load(model)
model
```

Now, let's assess how often the estimated 95% posterior intervals capture the true values of `beta`. If the model is implemented correctly, the coverage value below should be close to 95%. (Ordinarily, we would [increase the number of batches and reps per batch](https://books.ropensci.org/targets/dynamic.html#batching) and [run batches in parallel computing](https://books.ropensci.org/targets/hpc.html).)

```{r}
library(dplyr)
model %>%
  group_by(variable) %>%
  summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
```
For maximum reproducibility, we should express the coverage assessment as a custom function and a target in the pipeline.

```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  simulate_data <- function(n = 10L) {
    beta <- rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- rnorm(n, beta[1] + x * beta[2], 1)
    list(
      n = n,
      x = x,
      y = y,
      .join_data = list(beta = beta)
    )
  }
  list(
    tar_stan_mcmc_rep_summary(
      model,
      "model.stan",
      simulate_data(),
      batches = 5, # Number of branch targets.
      reps = 2, # Number of model reps per branch target.
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      ),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    ),
    tar_target(
      coverage,
      model %>%
        group_by(variable) %>%
        summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

simulate_data <- function(n = 10L) {
  beta <- rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- rnorm(n, beta[1] + x * beta[2], 1)
  list(
    n = n,
    x = x,
    y = y,
    .join_data = list(beta = beta)
  )
}

list(
  tar_stan_mcmc_rep_summary(
    model,
    "model.stan",
    simulate_data(),
    batches = 5, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    ),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_target(
    coverage,
    model %>%
      group_by(variable) %>%
      summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
  )
)
```

The new `coverage` target should the only outdated target, and it should be connected to the upstream `model` target.

```{r}
tar_visnetwork()
```

When we run the pipeline, only the coverage assessment should run. That way, we skip all the expensive computation of simulating datasets and running MCMC multiple times.

```{r, output = FALSE, warning = FALSE}
tar_make()
```

```{r}
tar_read(coverage)
```

## Multiple models

`tar_stan_rep_mcmc_summary()` and similar functions allow you to supply multiple Stan models. If you do, each model will share the the same collection of datasets. Suppose we have a new model, `model2.stan`.

```{r}
lines <- "data {
  int <lower = 1> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] beta;
}
model {
  y ~ normal(beta[1] + x * x * beta[2], 1); // Regress on x^2 instead of x.
  beta ~ normal(0, 1);
}"
writeLines(lines, "model2.stan")
```

To set up the simulation workflow to run on both models, we add `model2.stan` to the `stan_files` argument of `tar_stan_rep_mcmc_summary()`. And in the coverage summary below, we group by `.name` to compute a coverage statistic for each model.


```{r, echo = FALSE}
library(targets)
tar_script({
  library(stantargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  simulate_data <- function(n = 10L) {
    beta <- rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- rnorm(n, beta[1] + x * beta[2], 1)
    list(
      n = n,
      x = x,
      y = y,
      .join_data = list(beta = beta)
    )
  }
  list(
    tar_stan_mcmc_rep_summary(
      model,
      c("model.stan", "model2.stan"), # another model
      simulate_data(),
      batches = 5,
      reps = 2,
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      ),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    ),
    tar_target(
      coverage,
      model %>%
        group_by(.name, variable) %>%
        summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(stantargets)

simulate_data <- function(n = 10L) {
  beta <- rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- rnorm(n, beta[1] + x * beta[2], 1)
  list(
    n = n,
    x = x,
    y = y,
    .join_data = list(beta = beta)
  )
}

list(
  tar_stan_mcmc_rep_summary(
    model,
    c("model.stan", "model2.stan"), # another model
    simulate_data(),
    batches = 5,
    reps = 2,
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    ),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_target(
    coverage,
    model %>%
      group_by(.name, variable) %>%
      summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
  )
)
```

In the graph below, notice how targets `model_model` and `model_model2` are both connected to `model_data` upstream. Downstream, `model` is equivalent to `dplyr::bind_rows(model_model, model_model2)`, and it will have special columns `.name` and `.file` to distinguish among all the models.

```{r}
tar_visnetwork()
```

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_package.R
\docType{package}
\name{stantargets-package}
\alias{stantargets-package}
\title{targets: Targets Archetypes for Stan}
\description{
Bayesian data analysis usually incurs long runtimes
and cumbersome custom code. A pipeline toolkit tailored to
Bayesian statisticians, the \code{stantargets} R package leverages
\code{targets} and \code{cmdstanr} to ease these burdens.
\code{stantargets} makes it super easy to set up scalable
Stan pipelines that automatically parallelize the computation
and skip expensive steps when the results are already up to date.
Minimal custom code is required, and there is no need to manually
configure branching, so usage is much easier than \code{targets} alone.
\code{stantargets} can access all of \code{cmdstanr}'s major algorithms
(MCMC, variational Bayes, and optimization) and it supports
both single-fit workflows and multi-rep simulation studies.
}
\seealso{
\url{https://docs.ropensci.org/stantargets/}, \code{\link[=tar_stan_mcmc]{tar_stan_mcmc()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle_rep_summary.R
\name{tar_stan_mle_rep_summary}
\alias{tar_stan_mle_rep_summary}
\title{Multiple optimization runs per model with summaries}
\usage{
tar_stan_mle_rep_summary(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  init_alpha = NULL,
  iter = NULL,
  tol_obj = NULL,
  tol_rel_obj = NULL,
  tol_grad = NULL,
  tol_rel_grad = NULL,
  tol_param = NULL,
  history_size = NULL,
  sig_figs = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mle_rep_summaries()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
The specific target objects returned by
\code{tar_stan_mle_rep_summary(name = x, , stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run maximum likelihood
once per dataset.
Each dynamic branch returns a tidy data frames of maximum likelihood
estimates corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of maximum likelihood estimates.
}
}
\description{
\code{tar_stan_mle_rep_summaries()} creates targets
to run maximum likelihood multiple times per model and
save the MLEs in a long-form summary-like data frame.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$optimize()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mle_rep_summary(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other optimization: 
\code{\link{tar_stan_mle_rep_draws}()},
\code{\link{tar_stan_mle}()}
}
\concept{optimization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc_rep.R
\name{tar_stan_mcmc_rep}
\alias{tar_stan_mcmc_rep}
\title{Multiple MCMCs per model with tidy output}
\usage{
tar_stan_mcmc_rep(
  name,
  stan_files,
  data = list(),
  output_type = c("summary", "draws", "diagnostics"),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  chain_ids = seq_len(chains),
  threads_per_chain = NULL,
  iter_warmup = NULL,
  iter_sampling = NULL,
  save_warmup = FALSE,
  thin = NULL,
  max_treedepth = NULL,
  adapt_engaged = TRUE,
  adapt_delta = NULL,
  step_size = NULL,
  metric = NULL,
  metric_file = NULL,
  inv_metric = NULL,
  init_buffer = NULL,
  term_buffer = NULL,
  window = NULL,
  fixed_param = FALSE,
  sig_figs = NULL,
  validate_csv = TRUE,
  show_messages = TRUE,
  data_copy = character(0),
  variables = NULL,
  inc_warmup = FALSE,
  summaries = NULL,
  summary_args = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{Code to generate a single replication of a simulated dataset.
The workflow simulates multiple datasets, and each
model runs on each dataset. To join data on to the model
summaries, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{inc_warmup}{(logical) Should warmup draws be included? Defaults to
\code{FALSE}. Ignored except when used with \link[cmdstanr]{CmdStanMCMC} objects.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Developers can consult the design specification at
\url{https://books.ropensci.org/targets-design/}
to learn about the structure and composition of target objects.
}
\description{
Internal function for replicated MCMC.
Users should not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_output.R
\name{tar_stan_output}
\alias{tar_stan_output}
\title{Post-process Stan output}
\usage{
tar_stan_output(
  fit,
  output_type,
  summaries,
  summary_args,
  variables,
  inc_warmup,
  data,
  data_copy
)
}
\arguments{
\item{fit}{A Stan fit object.}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{variables}{(character vector) The variables to include.}

\item{inc_warmup}{Logical, whether to include the warmup draws.}

\item{data}{List, Stan dataset.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}
}
\value{
A data frame of user-friendly Stan output.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle_rep.R
\name{tar_stan_mle_rep_run}
\alias{tar_stan_mle_rep_run}
\title{Run a Stan model and return only the summaries.}
\usage{
tar_stan_mle_rep_run(
  stan_file,
  stan_name,
  stan_path,
  data,
  output_type,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  algorithm,
  init_alpha,
  iter,
  sig_figs,
  tol_obj,
  tol_rel_obj,
  tol_grad,
  tol_rel_grad,
  tol_param,
  history_size,
  data_copy,
  variables,
  summaries,
  summary_args
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}
}
\value{
A data frame of posterior summaries.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc_rep.R
\name{tar_stan_mcmc_rep_run}
\alias{tar_stan_mcmc_rep_run}
\title{Run a Stan model and return only the summaries.}
\usage{
tar_stan_mcmc_rep_run(
  stan_file,
  stan_name,
  stan_path,
  data,
  output_type,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  chains,
  parallel_chains,
  chain_ids,
  threads_per_chain,
  iter_warmup,
  iter_sampling,
  save_warmup,
  thin,
  max_treedepth,
  adapt_engaged,
  adapt_delta,
  step_size,
  metric,
  metric_file,
  inv_metric,
  init_buffer,
  term_buffer,
  window,
  fixed_param,
  sig_figs,
  validate_csv,
  show_messages,
  data_copy,
  inc_warmup,
  variables,
  summaries,
  summary_args
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{stan_name}{Friendly suffix of the Stan model target.}

\item{stan_path}{Original path to the input Stan file.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{inc_warmup}{(logical) Should warmup draws be included? Defaults to
\code{FALSE}. Ignored except when used with \link[cmdstanr]{CmdStanMCMC} objects.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}
}
\value{
A data frame of posterior summaries.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_example_data.R
\name{tar_stan_example_data}
\alias{tar_stan_example_data}
\title{Simulate example data for \code{\link[=tar_stan_example_file]{tar_stan_example_file()}}.}
\format{
A list with the following elements:
\itemize{
\item \code{n}: integer, number of data points.
\item \code{x}: numeric, covariate vector.
\item \code{y}: numeric, response variable.
\item \code{true_beta}: numeric of length 1, value of the regression
coefficient \code{beta} used in simulation.
\item \code{.join_data}: a list of simulated values to be appended
to as a \code{.join_data} column in the output of
targets generated by functions such as
\code{\link[=tar_stan_mcmc_rep_summary]{tar_stan_mcmc_rep_summary()}}. Contains the
regression coefficient \code{beta} (numeric of length 1)
and prior predictive data \code{y} (numeric vector).
}
}
\usage{
tar_stan_example_data(n = 10L)
}
\arguments{
\item{n}{Integer of length 1, number of data points.}
}
\value{
List, dataset compatible with the model file from
\code{\link[=tar_stan_example_file]{tar_stan_example_file()}}.
}
\description{
An example dataset compatible with the model file
from \code{\link[=tar_stan_example_file]{tar_stan_example_file()}}.
}
\details{
The \code{tar_stan_example_data()} function draws a Stan
dataset from the prior predictive distribution of the
model from \code{\link[=tar_stan_example_file]{tar_stan_example_file()}}. First, the
regression coefficient \code{beta} is drawn from its standard
normal prior, and the covariate \code{x} is computed.
Then, conditional on the \code{beta} draws and the covariate,
the response vector \code{y} is drawn from its
Normal(\code{x * beta}, 1) likelihood.
}
\examples{
tar_stan_example_data()
}
\seealso{
Other examples: 
\code{\link{tar_stan_example_file}()}
}
\concept{examples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb_rep.R
\name{tar_stan_vb_rep_run}
\alias{tar_stan_vb_rep_run}
\title{Run a Stan model and return only the summaries.}
\usage{
tar_stan_vb_rep_run(
  stan_file,
  stan_name,
  stan_path,
  data,
  output_type,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  algorithm,
  iter,
  grad_samples,
  elbo_samples,
  eta,
  adapt_engaged,
  adapt_iter,
  tol_rel_obj,
  eval_elbo,
  output_samples,
  sig_figs,
  data_copy,
  variables,
  summaries,
  summary_args
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{stan_name}{Friendly suffix of the Stan model target.}

\item{stan_path}{Original path to the input Stan file.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation?}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}
}
\value{
A data frame of posterior summaries.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_compile.R
\name{tar_stan_compile}
\alias{tar_stan_compile}
\title{Stan model compilation}
\usage{
tar_stan_compile(
  name,
  stan_file,
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_compile()} returns a target object to compile a Stan file.
The return value of this target is a character vector
containing the Stan model source file and compiled
executable file. A change in either file
will cause the target to rerun in the next run of the pipeline.
See the "Target objects" section for background.
}
\description{
\code{tar_stan_compile()} creates a target
to compile a Stan model and return the
original Stan model file. Does not compile the model
if the compilation is already up to date.
}
\details{
Most of the arguments are passed to the
\verb{$compile()} method of the \code{CmdStanModel} class.
For details, visit \url{https://mc-stan.org/cmdstanr/reference/}.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(tar_stan_compile(compiled_model, path))
})
targets::tar_make()
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc_rep_summary.R
\name{tar_stan_mcmc_rep_summary}
\alias{tar_stan_mcmc_rep_summary}
\title{Multiple MCMCs per model with summaries}
\usage{
tar_stan_mcmc_rep_summary(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  chain_ids = seq_len(chains),
  threads_per_chain = NULL,
  iter_warmup = NULL,
  iter_sampling = NULL,
  save_warmup = FALSE,
  thin = NULL,
  max_treedepth = NULL,
  adapt_engaged = TRUE,
  adapt_delta = NULL,
  step_size = NULL,
  metric = NULL,
  metric_file = NULL,
  inv_metric = NULL,
  init_buffer = NULL,
  term_buffer = NULL,
  window = NULL,
  fixed_param = FALSE,
  sig_figs = NULL,
  validate_csv = TRUE,
  show_messages = TRUE,
  data_copy = character(0),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{Code to generate a single replication of a simulated dataset.
The workflow simulates multiple datasets, and each
model runs on each dataset. To join data on to the model
summaries, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mcmc_rep_summary()} returns a list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mcmc_rep_summary(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to the
model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run MCMC once per dataset.
Each dynamic branch returns a tidy data frames of summaries.
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of summaries.
}
}
\description{
Targets to run MCMC multiple times and
save only the summary output from each run.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$sample()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc_rep_summary(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other MCMC: 
\code{\link{tar_stan_mcmc_rep_diagnostics}()},
\code{\link{tar_stan_mcmc_rep_draws}()},
\code{\link{tar_stan_mcmc}()}
}
\concept{MCMC}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq_rep_summary.R
\name{tar_stan_gq_rep_summary}
\alias{tar_stan_gq_rep_summary}
\title{Multiple runs of generated quantities per model with summaries}
\usage{
tar_stan_gq_rep_summary(
  name,
  stan_files,
  data = list(),
  fitted_params,
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  output_dir = NULL,
  sig_figs = NULL,
  parallel_chains = getOption("mc.cores", 1),
  threads_per_chain = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{fitted_params}{(multiple options) The parameter draws to use. One of
the following:
\itemize{
\item A \link[cmdstanr]{CmdStanMCMC} or \link[cmdstanr]{CmdStanVB} fitted model object.
\item A \link[posterior:draws_array]{posterior::draws_array} (for MCMC) or \link[posterior:draws_matrix]{posterior::draws_matrix} (for
VB) object returned by CmdStanR's \code{\link[cmdstanr:fit-method-draws]{$draws()}} method.
\item A character vector of paths to CmdStan CSV output files.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_gq_rep_summaries()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_gq_rep_summary(name = x, stan_files = "y.stan")}
returns a list of target objects:
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to the
model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run generated quantities
once per dataset.
Each dynamic branch returns a tidy data frames of summaries
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of summaries.
}
}
\description{
\code{tar_stan_gq_rep_summaries()} creates targets
to run generated quantities multiple times and
save only the summaries from each run.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$generate_quantities()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc(
    your_model,
    stan_files = c(x = path),
    data = tar_stan_example_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_stan_gq_rep_summary(
    generated_quantities,
    stan_files = path,
    data = tar_stan_example_data(),
    fitted_params = your_model_mcmc_x,
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other generated quantities: 
\code{\link{tar_stan_gq_rep_draws}()},
\code{\link{tar_stan_gq}()}
}
\concept{generated quantities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb_rep.R
\name{tar_stan_vb_rep}
\alias{tar_stan_vb_rep}
\title{Multiple iterations per model of variational Bayes with tidy output}
\usage{
tar_stan_vb_rep(
  name,
  stan_files,
  data = list(),
  output_type = c("summary", "draws"),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  iter = NULL,
  grad_samples = NULL,
  elbo_samples = NULL,
  eta = NULL,
  adapt_engaged = NULL,
  adapt_iter = NULL,
  tol_rel_obj = NULL,
  eval_elbo = NULL,
  output_samples = NULL,
  sig_figs = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation?}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Developers can consult the design specification at
\url{https://books.ropensci.org/targets-design/}
to learn about the structure and composition of target objects.
}
\description{
Internal function. Users should not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_compile.R
\name{tar_stan_compile_run}
\alias{tar_stan_compile_run}
\title{Compile a Stan model and return the model file.}
\usage{
tar_stan_compile_run(
  stan_file,
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}
}
\value{
Character of length 1, the value of \code{stan_file}.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq.R
\name{tar_stan_gq_run}
\alias{tar_stan_gq_run}
\title{Compile and run a Stan model and return the \code{CmdStanFit} object.}
\usage{
tar_stan_gq_run(
  stan_file,
  data,
  fitted_params,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  output_dir,
  sig_figs,
  parallel_chains,
  threads_per_chain,
  variables
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{fitted_params}{(multiple options) The parameter draws to use. One of
the following:
\itemize{
\item A \link[cmdstanr]{CmdStanMCMC} or \link[cmdstanr]{CmdStanVB} fitted model object.
\item A \link[posterior:draws_array]{posterior::draws_array} (for MCMC) or \link[posterior:draws_matrix]{posterior::draws_matrix} (for
VB) object returned by CmdStanR's \code{\link[cmdstanr:fit-method-draws]{$draws()}} method.
\item A character vector of paths to CmdStan CSV output files.
}}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}
}
\value{
A \code{CmdStanFit} object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb.R
\name{tar_stan_vb}
\alias{tar_stan_vb}
\title{One variational Bayes run per model with multiple outputs}
\usage{
tar_stan_vb(
  name,
  stan_files,
  data = list(),
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  iter = NULL,
  grad_samples = NULL,
  elbo_samples = NULL,
  eta = NULL,
  adapt_engaged = NULL,
  adapt_iter = NULL,
  tol_rel_obj = NULL,
  eval_elbo = NULL,
  output_samples = NULL,
  sig_figs = NULL,
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  draws = TRUE,
  summary = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of Stan model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{fs::path_ext_remove(basename(stan_files))} will be used
as target name suffixes. If \code{stan_files} is a named vector,
the suffixed will come from \code{names(stan_files)}.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{variables}{(character vector) The variables to include.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{draws}{Logical, whether to create a target for posterior draws.
Saves \code{posterior::as_draws_df(fit$draws())} to a compressed \code{tibble}.
Convenient, but duplicates storage.}

\item{summary}{Logical, whether to create a target for
\code{fit$summary()}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_vb()} returns a list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_vb(name = x, stan_files = "y.stan", ...)} are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: run the R expression in the \code{data} argument to produce
a Stan dataset for the model. Returns a Stan data list.
\item \code{x_vb_y}: run variational Bayes on the model and the dataset.
Returns a \code{cmdstanr} \code{CmdStanVB} object with all the results.
\item \code{x_draws_y}: extract draws from \code{x_vb_y}.
Omitted if \code{draws = FALSE}.
Returns a tidy data frame of draws.
\item \code{x_summary_y}: extract compact summaries from \code{x_vb_y}.
Returns a tidy data frame of summaries.
Omitted if \code{summary = FALSE}.
}
}
\description{
Targets to run a Stan model once with
variational Bayes and save multiple outputs.
}
\details{
Most of the arguments are passed to the \verb{$compile()},
\verb{$variational()}, and \verb{$summary()} methods of the \code{CmdStanModel} class.
If you previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_vb(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    variables = "beta",
    summaries = list(~quantile(.x, probs = c(0.25, 0.75))),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other variational Bayes: 
\code{\link{tar_stan_vb_rep_draws}()},
\code{\link{tar_stan_vb_rep_summary}()}
}
\concept{variational Bayes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle.R
\name{tar_stan_mle}
\alias{tar_stan_mle}
\title{One optimization run per model with multiple outputs}
\usage{
tar_stan_mle(
  name,
  stan_files,
  data = list(),
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  init_alpha = NULL,
  iter = NULL,
  tol_obj = NULL,
  tol_rel_obj = NULL,
  tol_grad = NULL,
  tol_rel_grad = NULL,
  tol_param = NULL,
  history_size = NULL,
  sig_figs = NULL,
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  draws = TRUE,
  summary = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of Stan model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{fs::path_ext_remove(basename(stan_files))} will be used
as target name suffixes. If \code{stan_files} is a named vector,
the suffixed will come from \code{names(stan_files)}.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{variables}{(character vector) The variables to include.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{draws}{Logical, whether to create a target for posterior draws.
Saves \code{posterior::as_draws_df(fit$draws())} to a compressed \code{tibble}.
Convenient, but duplicates storage.}

\item{summary}{Logical, whether to create a target for
\code{fit$summary()}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mle()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mle(name = x, stan_files = "y.stan", ...)}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: run the R expression in the \code{data} argument to produce
a Stan dataset for the model. Returns a Stan data list.
\item \code{x_mle_y}: run generated quantities on the model and the dataset.
Returns a \code{cmdstanr} \code{CmdStanGQ} object with all the results.
\item \code{x_draws_y}: extract maximum likelihood estimates from \code{x_mle_y}
in draws format.
Omitted if \code{draws = FALSE}.
Returns a wide data frame of MLEs.
\item \code{x_summary_y}: extract MLEs from from \code{x_mle_y} in summary format.
Returns a long data frame of MLEs.
Omitted if \code{summary = FALSE}.
}
}
\description{
\code{tar_stan_mle()} creates targets to optimize a Stan model once
per model and separately save draws-like output and summary-like output.
}
\details{
Most of the arguments are passed to the \verb{$compile()},
\verb{$optimize()}, and \verb{$summary()} methods of the \code{CmdStanModel} class.
If you previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mle(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other optimization: 
\code{\link{tar_stan_mle_rep_draws}()},
\code{\link{tar_stan_mle_rep_summary}()}
}
\concept{optimization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq_rep.R
\name{tar_stan_gq_rep_run}
\alias{tar_stan_gq_rep_run}
\title{Run a Stan model and return only the summaries.}
\usage{
tar_stan_gq_rep_run(
  stan_file,
  stan_name,
  stan_path,
  data,
  output_type,
  fitted_params,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  output_dir,
  sig_figs,
  parallel_chains,
  threads_per_chain,
  data_copy,
  variables,
  summaries,
  summary_args
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{stan_name}{Friendly suffix of the Stan model target.}

\item{stan_path}{Original path to the input Stan file.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{fitted_params}{(multiple options) The parameter draws to use. One of
the following:
\itemize{
\item A \link[cmdstanr]{CmdStanMCMC} or \link[cmdstanr]{CmdStanVB} fitted model object.
\item A \link[posterior:draws_array]{posterior::draws_array} (for MCMC) or \link[posterior:draws_matrix]{posterior::draws_matrix} (for
VB) object returned by CmdStanR's \code{\link[cmdstanr:fit-method-draws]{$draws()}} method.
\item A character vector of paths to CmdStan CSV output files.
}}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}
}
\value{
A data frame of posterior summaries.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_summary.R
\name{tar_stan_summary}
\alias{tar_stan_summary}
\title{One summary of a \code{CmdStanFit} object}
\usage{
tar_stan_summary(
  name,
  fit,
  data = NULL,
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  format = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{fit}{Symbol, name of a \code{CmdStanFit} object or an upstream target
that returns a \code{CmdStanFit} object.}

\item{data}{Code to generate the \code{data} for the Stan model.}

\item{variables}{(character vector) The variables to include.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_summary()} returns target object to
summarize a \code{CmdStanFit} object. The return value of the target
is a tidy data frame of summaries returned by the \verb{$summary()}
method of the \code{CmdStanFit} object.
See the "Target objects" section for background.
}
\description{
Create a target to run the \verb{$summary()}
method of a \code{CmdStanFit} object.
}
\details{
\code{\link[=tar_stan_mcmc]{tar_stan_mcmc()}} etc. with \code{summary = TRUE} already gives you a
target with output from the \verb{$summary()} method.
Use \code{tar_stan_summary()} to create additional specialized summaries.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
# First, write your Stan model file, e.g. model.stan.
# Then in _targets.R, write a pipeline like this:
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
# Running inside a temporary directory to avoid
# modifying the user's file space. The file "model.stan"
# created below lives in a temporary directory.
# This satisfies CRAN policies.
tar_stan_example_file("model.stan")
targets::tar_script({
library(stantargets)
list(
  # Run a model and produce default summaries.
  tar_stan_mcmc(
    your_model,
    stan_files = "model.stan",
    data = tar_stan_example_data()
  ),
  # Produce a more specialized summary
  tar_stan_summary(
    your_summary,
    fit = your_model_mcmc_model,
    data = your_model_data_model,
    variables = "beta",
    summaries = list(~quantile(.x, probs = c(0.25, 0.75)))
  )
)}, ask = FALSE)
targets::tar_make()
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq_rep_draws.R
\name{tar_stan_gq_rep_draws}
\alias{tar_stan_gq_rep_draws}
\title{Multiple runs of generated quantities per model with draws}
\usage{
tar_stan_gq_rep_draws(
  name,
  stan_files,
  data = list(),
  fitted_params,
  batches = 1L,
  reps = 1L,
  combine = FALSE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  output_dir = NULL,
  sig_figs = NULL,
  parallel_chains = getOption("mc.cores", 1),
  threads_per_chain = NULL,
  variables = NULL,
  data_copy = character(0),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
  garbage_collection = TRUE,
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{fitted_params}{(multiple options) The parameter draws to use. One of
the following:
\itemize{
\item A \link[cmdstanr]{CmdStanMCMC} or \link[cmdstanr]{CmdStanVB} fitted model object.
\item A \link[posterior:draws_array]{posterior::draws_array} (for MCMC) or \link[posterior:draws_matrix]{posterior::draws_matrix} (for
VB) object returned by CmdStanR's \code{\link[cmdstanr:fit-method-draws]{$draws()}} method.
\item A character vector of paths to CmdStan CSV output files.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_gq_rep_draws()} returns a list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_gq_rep_draws(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run generated quantities
once per dataset.
Each dynamic branch returns a tidy data frames of draws
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of draws.
}
}
\description{
\code{tar_stan_gq_rep_draws()} creates targets
to run generated quantities multiple times and
save only the draws from each run.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$sample()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc(
    your_model,
    stan_files = c(x = path),
    data = tar_stan_example_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile(),
    refresh = 0
  ),
  tar_stan_gq_rep_draws(
    generated_quantities,
    stan_files = path,
    data = tar_stan_example_data(),
    fitted_params = your_model_mcmc_x,
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other generated quantities: 
\code{\link{tar_stan_gq_rep_summary}()},
\code{\link{tar_stan_gq}()}
}
\concept{generated quantities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc.R
\name{tar_stan_mcmc}
\alias{tar_stan_mcmc}
\title{One MCMC per model with multiple outputs}
\usage{
tar_stan_mcmc(
  name,
  stan_files,
  data = list(),
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  chain_ids = seq_len(chains),
  threads_per_chain = NULL,
  iter_warmup = NULL,
  iter_sampling = NULL,
  save_warmup = FALSE,
  thin = NULL,
  max_treedepth = NULL,
  adapt_engaged = TRUE,
  adapt_delta = NULL,
  step_size = NULL,
  metric = NULL,
  metric_file = NULL,
  inv_metric = NULL,
  init_buffer = NULL,
  term_buffer = NULL,
  window = NULL,
  fixed_param = FALSE,
  sig_figs = NULL,
  validate_csv = TRUE,
  show_messages = TRUE,
  variables = NULL,
  inc_warmup = FALSE,
  summaries = list(),
  summary_args = list(),
  draws = TRUE,
  diagnostics = TRUE,
  summary = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of Stan model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{fs::path_ext_remove(basename(stan_files))} will be used
as target name suffixes. If \code{stan_files} is a named vector,
the suffixed will come from \code{names(stan_files)}.}

\item{data}{Code to generate the \code{data} for the Stan model.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{inc_warmup}{(logical) Should warmup draws be included? Defaults to
\code{FALSE}. Ignored except when used with \link[cmdstanr]{CmdStanMCMC} objects.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{draws}{Logical, whether to create a target for posterior draws.
Saves \code{posterior::as_draws_df(fit$draws())} to a compressed \code{tibble}.
Convenient, but duplicates storage.}

\item{diagnostics}{Logical, whether to create a target for
\code{posterior::as_draws_df(fit$sampler_diagnostics())}.
Saves \code{posterior::as_draws_df(fit$draws())} to a compressed \code{tibble}.
Convenient, but duplicates storage.}

\item{summary}{Logical, whether to create a target for
\code{fit$summary()}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the non-data-frame
targets such as the Stan data and any CmdStanFit objects.
Please choose an all=purpose
format such as \code{"qs"} or \code{"aws_qs"} rather than a file format like
\code{"file"} or a data frame format like \code{"parquet"}. For more on storage
formats, see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mcmc()} returns a list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mcmc(name = x, stan_files = "y.stan", ...)}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to the
model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: run the R expression in the \code{data} argument to produce
a Stan dataset for the model. Returns a Stan data list.
\item \code{x_mcmc_y}: run MCMC on the model and the dataset.
Returns a \code{cmdstanr} \code{CmdStanMCMC} object with all the results.
\item \code{x_draws_y}: extract draws from \code{x_mcmc_y}.
Omitted if \code{draws = FALSE}.
Returns a tidy data frame of draws.
\item \code{x_summary_y}: extract compact summaries from \code{x_mcmc_y}.
Returns a tidy data frame of summaries.
Omitted if \code{summary = FALSE}.
\item \code{x_diagnostics}: extract HMC diagnostics from \code{x_mcmc_y}.
Returns a tidy data frame of HMC diagnostics.
Omitted if \code{diagnostics = FALSE}.
}
}
\description{
\code{tar_stan_mcmc()} creates targets to run one MCMC
per model and separately save summaries draws, and diagnostics.
}
\details{
Most of the arguments are passed to the \verb{$compile()},
\verb{$sample()}, and \verb{$summary()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    variables = "beta",
    summaries = list(~quantile(.x, probs = c(0.25, 0.75))),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other MCMC: 
\code{\link{tar_stan_mcmc_rep_diagnostics}()},
\code{\link{tar_stan_mcmc_rep_draws}()},
\code{\link{tar_stan_mcmc_rep_summary}()}
}
\concept{MCMC}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb_rep_draws.R
\name{tar_stan_vb_rep_draws}
\alias{tar_stan_vb_rep_draws}
\title{Multiple variational Bayes runs per model with draws}
\usage{
tar_stan_vb_rep_draws(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = FALSE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  iter = NULL,
  grad_samples = NULL,
  elbo_samples = NULL,
  eta = NULL,
  adapt_engaged = NULL,
  adapt_iter = NULL,
  tol_rel_obj = NULL,
  eval_elbo = NULL,
  output_samples = NULL,
  sig_figs = NULL,
  data_copy = character(0),
  variables = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
  garbage_collection = TRUE,
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation?}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_vb_rep_summaries()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_vb_rep_draws(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run variational Bayes
once per dataset.
Each dynamic branch returns a tidy data frames of draws
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of draws.
}
}
\description{
\code{tar_stan_vb_rep_draws()} creates targets to run
variational Bayes multiple times per model and
save only the draws from each run.
}
\details{
Draws could take up a lot of storage. If storage becomes
excessive, please consider thinning the draws or using
\code{tar_stan_vb_rep_summaries()} instead.

Most of the arguments are passed to the \verb{$compile()}
and \verb{$variational()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_vb_rep_draws(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other variational Bayes: 
\code{\link{tar_stan_vb_rep_summary}()},
\code{\link{tar_stan_vb}()}
}
\concept{variational Bayes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_example_file.R
\name{tar_stan_example_file}
\alias{tar_stan_example_file}
\title{Write an example Stan model file.}
\usage{
tar_stan_example_file(path = tempfile(pattern = "", fileext = ".stan"))
}
\arguments{
\item{path}{Character of length 1, file path to write the model file.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Overwrites the file at \code{path} with a built-in example
Stan model file.
}
\examples{
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
writeLines(readLines(path))
}
\seealso{
Other examples: 
\code{\link{tar_stan_example_data}()}
}
\concept{examples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc_rep_diagnostics.R
\name{tar_stan_mcmc_rep_diagnostics}
\alias{tar_stan_mcmc_rep_diagnostics}
\title{Multiple MCMCs per model with sampler diagnostics}
\usage{
tar_stan_mcmc_rep_diagnostics(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = FALSE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  chain_ids = seq_len(chains),
  threads_per_chain = NULL,
  iter_warmup = NULL,
  iter_sampling = NULL,
  save_warmup = FALSE,
  thin = NULL,
  max_treedepth = NULL,
  adapt_engaged = TRUE,
  adapt_delta = NULL,
  step_size = NULL,
  metric = NULL,
  metric_file = NULL,
  inv_metric = NULL,
  init_buffer = NULL,
  term_buffer = NULL,
  window = NULL,
  fixed_param = FALSE,
  sig_figs = NULL,
  validate_csv = TRUE,
  show_messages = TRUE,
  inc_warmup = FALSE,
  data_copy = character(0),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
  garbage_collection = TRUE,
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{Code to generate a single replication of a simulated dataset.
The workflow simulates multiple datasets, and each
model runs on each dataset. To join data on to the model
summaries, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{inc_warmup}{(logical) Should warmup draws be included? Defaults to
\code{FALSE}. Ignored except when used with \link[cmdstanr]{CmdStanMCMC} objects.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mcmc_rep_diagnostics()} returns
a list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mcmc_rep_diagnostics(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to the
model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run MCMC once per dataset.
Each dynamic branch returns a tidy data frames of HMC diagnostics
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of HMC diagnostics.
}
}
\description{
\code{tar_stan_mcmc_rep_diagnostics()} creates targets
to run MCMC multiple times per model and save only the sampler
diagnostics from each run.
}
\details{
Saved diagnostics could get quite large in storage,
so please use thinning if necessary.

Most of the arguments are passed to the \verb{$compile()}
and \verb{$generate_quantities()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc_rep_diagnostics(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other MCMC: 
\code{\link{tar_stan_mcmc_rep_draws}()},
\code{\link{tar_stan_mcmc_rep_summary}()},
\code{\link{tar_stan_mcmc}()}
}
\concept{MCMC}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq.R
\name{tar_stan_gq}
\alias{tar_stan_gq}
\title{Generated quantities on an existing CmdStanFit object}
\usage{
tar_stan_gq(
  name,
  stan_files,
  data = list(),
  fitted_params,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  output_dir = NULL,
  sig_figs = NULL,
  parallel_chains = getOption("mc.cores", 1),
  threads_per_chain = NULL,
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  draws = TRUE,
  summary = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of Stan model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{fs::path_ext_remove(basename(stan_files))} will be used
as target name suffixes. If \code{stan_files} is a named vector,
the suffixed will come from \code{names(stan_files)}.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{fitted_params}{Symbol, name of a \code{CmdStanFit} object
computed in a previous target: for example, the
\verb{*_mcmc_*} target from \code{\link[=tar_stan_mcmc]{tar_stan_mcmc()}}. Must be a subclass
that \verb{$generate_quantities()} can accept as \code{fitted_params}.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{variables}{(character vector) The variables to include.}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{draws}{Logical, whether to create a target for posterior draws.
Saves \code{posterior::as_draws_df(fit$draws())} to a compressed \code{tibble}.
Convenient, but duplicates storage.}

\item{summary}{Logical, whether to create a target for
\code{fit$summary()}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_gq()} returns list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_gq(name = x, stan_files = "y.stan", ...)}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to the model
file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: run the R expression in the \code{data} argument to produce
a Stan dataset for the model. Returns a Stan data list.
\item \code{x_gq_y}: run generated quantities on the model and the dataset.
Returns a \code{cmdstanr} \code{CmdStanGQ} object with all the results.
\item \code{x_draws_y}: extract draws from \code{x_gq_y}.
Omitted if \code{draws = FALSE}.
Returns a tidy data frame of draws.
\item \code{x_summary_y}: extract compact summaries from \code{x_gq_y}.
Returns a tidy data frame of summaries.
Omitted if \code{summary = FALSE}.
}
}
\description{
\code{tar_stan_gq()} creates targets to run
the generated quantities of a Stan model and save
draws and summaries separately.
}
\details{
Most of the arguments are passed to the \verb{$compile()},
\verb{$generate_quantities()}, and \verb{$summary()} methods
of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc(
    your_model,
    stan_files = c(x = path),
    data = tar_stan_example_data(),
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  ),
  tar_stan_gq(
    custom_gq,
    stan_files = path, # Can be a different model.
    fitted_params = your_model_mcmc_x,
    data = your_model_data, # Can be a different dataset.
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other generated quantities: 
\code{\link{tar_stan_gq_rep_draws}()},
\code{\link{tar_stan_gq_rep_summary}()}
}
\concept{generated quantities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc.R
\name{tar_stan_mcmc_run}
\alias{tar_stan_mcmc_run}
\title{Compile and run a Stan model and return the \code{CmdStanFit} object.}
\usage{
tar_stan_mcmc_run(
  stan_file,
  data,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  chains,
  parallel_chains,
  chain_ids,
  threads_per_chain,
  iter_warmup,
  iter_sampling,
  save_warmup,
  thin,
  max_treedepth,
  adapt_engaged,
  adapt_delta,
  step_size,
  metric,
  metric_file,
  inv_metric,
  init_buffer,
  term_buffer,
  window,
  fixed_param,
  sig_figs,
  validate_csv,
  show_messages,
  variables,
  inc_warmup
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{Character of length 1. If \code{"original"}, then
\code{cmdstan} will compile the source file right before running
it (or skip compilation if the binary is up to date). This
assumes the worker has access to the file. If the worker
is running on a remote computer that does not have access
to the model file, set to \code{"copy"} instead. \code{compile = "copy"}
means the pipeline will read the lines of the original Stan model file
and send them to the worker. The worker writes the lines
to a local copy and compiles the model from there, so it
no longer needs access to the original Stan model file on your
local machine. However, as a result, the Stan model re-compiles
every time the main target reruns.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}
}
\value{
A \code{CmdStanFit} object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle_rep_draws.R
\name{tar_stan_mle_rep_draws}
\alias{tar_stan_mle_rep_draws}
\title{Multiple optimization runs per model with draws}
\usage{
tar_stan_mle_rep_draws(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  init_alpha = NULL,
  iter = NULL,
  sig_figs = NULL,
  tol_obj = NULL,
  tol_rel_obj = NULL,
  tol_grad = NULL,
  tol_rel_grad = NULL,
  tol_param = NULL,
  history_size = NULL,
  data_copy = character(0),
  variables = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mle_rep_draws()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mcmc_rep_draws(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run maximum likelihood
once per dataset.
Each dynamic branch returns a tidy data frames of maximum likelihood
estimates corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of maximum likelihood estimates.
}
}
\description{
\code{tar_stan_mle_rep_draws()} creates targets
to run maximum likelihood multiple times per model and
save the MLEs in a wide-form draws-like data frame.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$optimize()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mle_rep_draws(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other optimization: 
\code{\link{tar_stan_mle_rep_summary}()},
\code{\link{tar_stan_mle}()}
}
\concept{optimization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle.R
\name{tar_stan_mle_run}
\alias{tar_stan_mle_run}
\title{Compile and run a Stan model and return a \code{CmdStanMLE} object.}
\usage{
tar_stan_mle_run(
  stan_file,
  data,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  algorithm,
  init_alpha,
  iter,
  sig_figs,
  tol_obj,
  tol_rel_obj,
  tol_grad,
  tol_rel_grad,
  tol_param,
  history_size,
  variables
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{Character of length 1. If \code{"original"}, then
\code{cmdstan} will compile the source file right before running
it (or skip compilation if the binary is up to date). This
assumes the worker has access to the file. If the worker
is running on a remote computer that does not have access
to the model file, set to \code{"copy"} instead. \code{compile = "copy"}
means the pipeline will read the lines of the original Stan model file
and send them to the worker. The worker writes the lines
to a local copy and compiles the model from there, so it
no longer needs access to the original Stan model file on your
local machine. However, as a result, the Stan model re-compiles
every time the main target reruns.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}
}
\value{
A \code{CmdStanFit} object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mle_rep.R
\name{tar_stan_mle_rep}
\alias{tar_stan_mle_rep}
\title{Multiple optimization runs per model with tidy output}
\usage{
tar_stan_mle_rep(
  name,
  stan_files,
  data = list(),
  output_type = c("summary", "draws"),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  init_alpha = NULL,
  iter = NULL,
  tol_obj = NULL,
  tol_rel_obj = NULL,
  tol_grad = NULL,
  tol_rel_grad = NULL,
  tol_param = NULL,
  history_size = NULL,
  sig_figs = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The optimization algorithm. One of \code{"lbfgs"},
\code{"bfgs"}, or \code{"newton"}. The control parameters below are only available
for \code{"lbfgs"} and \verb{"bfgs}. For their default values and more details see
the CmdStan User's Guide. The default values can also be obtained by
running \code{cmdstanr_example(method="optimize")$metadata()}.}

\item{init_alpha}{(positive real) The initial step size parameter.}

\item{iter}{(positive integer) The maximum number of iterations.}

\item{tol_obj}{(positive real) Convergence tolerance on changes in objective function value.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on relative changes in objective function value.}

\item{tol_grad}{(positive real) Convergence tolerance on the norm of the gradient.}

\item{tol_rel_grad}{(positive real) Convergence tolerance on the relative norm of the gradient.}

\item{tol_param}{(positive real) Convergence tolerance on changes in parameter value.}

\item{history_size}{(positive integer) The size of the history used when
approximating the Hessian. Only available for L-BFGS.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Developers can consult the design specification at
\url{https://books.ropensci.org/targets-design/}
to learn about the structure and composition of target objects.
}
\description{
Internal function. Users should not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb.R
\name{tar_stan_vb_run}
\alias{tar_stan_vb_run}
\title{Compile and run a Stan model and return a \code{CmdStanVB} object.}
\usage{
tar_stan_vb_run(
  stan_file,
  data,
  compile,
  quiet,
  stdout,
  stderr,
  dir,
  pedantic,
  include_paths,
  cpp_options,
  stanc_options,
  force_recompile,
  seed,
  refresh,
  init,
  save_latent_dynamics,
  output_dir,
  algorithm,
  iter,
  grad_samples,
  elbo_samples,
  eta,
  adapt_engaged,
  adapt_iter,
  tol_rel_obj,
  eval_elbo,
  output_samples,
  sig_figs,
  variables
)
}
\arguments{
\item{stan_file}{(string) The path to a \code{.stan} file containing a Stan
program. The helper function \code{\link[cmdstanr:write_stan_file]{write_stan_file()}} is provided for cases when
it is more convenient to specify the Stan program as a string.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{compile}{Character of length 1. If \code{"original"}, then
\code{cmdstan} will compile the source file right before running
it (or skip compilation if the binary is up to date). This
assumes the worker has access to the file. If the worker
is running on a remote computer that does not have access
to the model file, set to \code{"copy"} instead. \code{compile = "copy"}
means the pipeline will read the lines of the original Stan model file
and send them to the worker. The worker writes the lines
to a local copy and compiles the model from there, so it
no longer needs access to the original Stan model file on your
local machine. However, as a result, the Stan model re-compiles
every time the main target reruns.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}
}
\value{
A \code{CmdStanFit} object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_vb_rep_summary.R
\name{tar_stan_vb_rep_summary}
\alias{tar_stan_vb_rep_summary}
\title{Multiple iterations per model of variational Bayes with summaries}
\usage{
tar_stan_vb_rep_summary(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  algorithm = NULL,
  iter = NULL,
  grad_samples = NULL,
  elbo_samples = NULL,
  eta = NULL,
  adapt_engaged = NULL,
  adapt_iter = NULL,
  tol_rel_obj = NULL,
  eval_elbo = NULL,
  output_samples = NULL,
  sig_figs = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = list(),
  summary_args = list(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{algorithm}{(string) The algorithm. Either \code{"meanfield"} or
\code{"fullrank"}.}

\item{iter}{(positive integer) The \emph{maximum} number of iterations.}

\item{grad_samples}{(positive integer) The number of samples for Monte Carlo
estimate of gradients.}

\item{elbo_samples}{(positive integer) The number of samples for Monte Carlo
estimate of ELBO (objective function).}

\item{eta}{(positive real) The step size weighting parameter for adaptive
step size sequence.}

\item{adapt_engaged}{(logical) Do warmup adaptation?}

\item{adapt_iter}{(positive integer) The \emph{maximum} number of adaptation
iterations.}

\item{tol_rel_obj}{(positive real) Convergence tolerance on the relative norm
of the objective.}

\item{eval_elbo}{(positive integer) Evaluate ELBO every Nth iteration.}

\item{output_samples}{(positive integer) Number of approximate posterior
samples to draw and save.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_vb_rep_summaries()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_vb_rep_summary(name = x,  stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with paths to
the model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run variational Bayes
once per dataset.
Each dynamic branch returns a tidy data frames of summaries
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of summaries.
}
}
\description{
\code{tar_stan_vb_rep_summaries()} creates targets
to run variational Bayes multiple times and
save only the summary output from each run.
}
\details{
Most of the arguments are passed to the \verb{$compile()}
and \verb{$variational()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_vb_rep_summary(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other variational Bayes: 
\code{\link{tar_stan_vb_rep_draws}()},
\code{\link{tar_stan_vb}()}
}
\concept{variational Bayes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_mcmc_rep_draws.R
\name{tar_stan_mcmc_rep_draws}
\alias{tar_stan_mcmc_rep_draws}
\title{Multiple MCMC runs per model with draws}
\usage{
tar_stan_mcmc_rep_draws(
  name,
  stan_files,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = FALSE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  refresh = NULL,
  init = NULL,
  save_latent_dynamics = FALSE,
  output_dir = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  chain_ids = seq_len(chains),
  threads_per_chain = NULL,
  iter_warmup = NULL,
  iter_sampling = NULL,
  save_warmup = FALSE,
  thin = NULL,
  max_treedepth = NULL,
  adapt_engaged = TRUE,
  adapt_delta = NULL,
  step_size = NULL,
  metric = NULL,
  metric_file = NULL,
  inv_metric = NULL,
  init_buffer = NULL,
  term_buffer = NULL,
  window = NULL,
  fixed_param = FALSE,
  sig_figs = NULL,
  validate_csv = TRUE,
  show_messages = TRUE,
  inc_warmup = FALSE,
  variables = NULL,
  data_copy = character(0),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
  garbage_collection = TRUE,
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{Code to generate a single replication of a simulated dataset.
The workflow simulates multiple datasets, and each
model runs on each dataset. To join data on to the model
summaries, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{refresh}{(non-negative integer) The number of iterations between
printed screen updates. If \code{refresh = 0}, only error messages will be
printed.}

\item{init}{(multiple options) The initialization method to use for the
variables declared in the parameters block of the Stan program:
\itemize{
\item A real number \code{x>0}. This initializes \emph{all} parameters randomly between
\verb{[-x,x]} (on the \emph{unconstrained} parameter space);
\item The number \code{0}. This initializes \emph{all} parameters to \code{0};
\item A character vector of paths (one per chain) to JSON or Rdump files
containing initial values for all or some parameters. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} to write \R objects to JSON files compatible with
CmdStan.
\item A list of lists containing initial values for all or some parameters. For
MCMC the list should contain a sublist for each chain. For optimization and
variational inference there should be just one sublist. The sublists should
have named elements corresponding to the parameters for which you are
specifying initial values. See \strong{Examples}.
\item A function that returns a single list with names corresponding to the
parameters for which you are specifying initial values. The function can
take no arguments or a single argument \code{chain_id}. For MCMC, if the function
has argument \code{chain_id} it will be supplied with the chain id (from 1 to
number of chains) when called to generate the initial values. See
\strong{Examples}.
}}

\item{save_latent_dynamics}{(logical) Should auxiliary diagnostic information
about the latent dynamics be written to temporary diagnostic CSV files?
This argument replaces CmdStan's \code{diagnostic_file} argument and the content
written to CSV is controlled by the user's CmdStan installation and not
CmdStanR (for some algorithms no content may be written). The default
is \code{FALSE}, which is appropriate for almost every use case. To save the
temporary files created when \code{save_latent_dynamics=TRUE} see the
\code{\link[cmdstanr:fit-method-save_output_files]{$save_latent_dynamics_files()}}
method.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{chains}{(positive integer) The number of Markov chains to run. The
default is 4.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{chain_ids}{(integer vector) A vector of chain IDs. Must contain
\code{chains} unique positive integers. If not set, the default chain IDs are
used (integers starting from \code{1}).}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{iter_warmup}{(positive integer) The number of warmup iterations to run
per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_warmup}.}

\item{iter_sampling}{(positive integer) The number of post-warmup iterations
to run per chain. Note: in the CmdStan User's Guide this is referred to as
\code{num_samples}.}

\item{save_warmup}{(logical) Should warmup iterations be saved? The default
is \code{FALSE}. If \code{save_warmup=TRUE} then you can use
\link[cmdstanr:fit-method-draws]{$draws(inc_warmup=TRUE)} to include warmup when
accessing the draws.}

\item{thin}{(positive integer) The period between saved samples. This should
typically be left at its default (no thinning) unless memory is a problem.}

\item{max_treedepth}{(positive integer) The maximum allowed tree depth for
the NUTS engine. See the \emph{Tree Depth} section of the CmdStan User's Guide
for more details.}

\item{adapt_engaged}{(logical) Do warmup adaptation? The default is \code{TRUE}.
If a precomputed inverse metric is specified via the \code{inv_metric} argument
(or \code{metric_file}) then, if \code{adapt_engaged=TRUE}, Stan will use the
provided inverse metric just as an initial guess during adaptation. To turn
off adaptation when using a precomputed inverse metric set
\code{adapt_engaged=FALSE}.}

\item{adapt_delta}{(real in \verb{(0,1)}) The adaptation target acceptance
statistic.}

\item{step_size}{(positive real) The \emph{initial} step size for the discrete
approximation to continuous Hamiltonian dynamics. This is further tuned
during warmup.}

\item{metric}{(string) One of \code{"diag_e"}, \code{"dense_e"}, or \code{"unit_e"},
specifying the geometry of the base manifold. See the \emph{Euclidean Metric}
section of the CmdStan User's Guide for more details. To specify a
precomputed (inverse) metric, see the \code{inv_metric} argument below.}

\item{metric_file}{(character vector) The paths to JSON or
Rdump files (one per chain) compatible with CmdStan that contain
precomputed inverse metrics. The \code{metric_file} argument is inherited from
CmdStan but is confusing in that the entry in JSON or Rdump file(s) must be
named \code{inv_metric}, referring to the \emph{inverse} metric. We recommend instead
using CmdStanR's \code{inv_metric} argument (see below) to specify an inverse
metric directly using a vector or matrix from your \R session.}

\item{inv_metric}{(vector, matrix) A vector (if \code{metric='diag_e'}) or a
matrix (if \code{metric='dense_e'}) for initializing the inverse metric. This
can be used as an alternative to the \code{metric_file} argument. A vector is
interpreted as a diagonal metric. The inverse metric is usually set to an
estimate of the posterior covariance. See the \code{adapt_engaged} argument
above for details about (and control over) how specifying a precomputed
inverse metric interacts with adaptation.}

\item{init_buffer}{(nonnegative integer) Width of initial fast timestep
adaptation interval during warmup.}

\item{term_buffer}{(nonnegative integer) Width of final fast timestep
adaptation interval during warmup.}

\item{window}{(nonnegative integer) Initial width of slow timestep/metric
adaptation interval.}

\item{fixed_param}{(logical) When \code{TRUE}, call CmdStan with argument
\code{"algorithm=fixed_param"}. The default is \code{FALSE}. The fixed parameter
sampler generates a new sample without changing the current state of the
Markov chain; only generated quantities may change. This can be useful
when, for example, trying to generate pseudo-data using the generated
quantities block. If the parameters block is empty then using
\code{fixed_param=TRUE} is mandatory. When \code{fixed_param=TRUE} the \code{chains} and
\code{parallel_chains} arguments will be set to \code{1}.}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{validate_csv}{(logical) When \code{TRUE} (the default), validate the
sampling results in the csv files. Disable if you wish to manually read in
the sampling results and validate them yourself, for example using
\code{\link[cmdstanr:read_cmdstan_csv]{read_cmdstan_csv()}}.}

\item{show_messages}{(logical) When \code{TRUE} (the default), prints all
informational messages, for example rejection of the current proposal.
Disable if you wish silence these messages, but this is not recommended
unless you are very sure that the model is correct up to numerical error.
If the messages are silenced then the \verb{$output()} method of the resulting
fit object can be used to display all the silenced messages.}

\item{inc_warmup}{(logical) Should warmup draws be included? Defaults to
\code{FALSE}. Ignored except when used with \link[cmdstanr]{CmdStanMCMC} objects.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
\code{tar_stan_mcmc_rep_draws()} returns a
list of target objects. See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{stan_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_stan_mcmc_rep_draws(name = x, stan_files = "y.stan")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the Stan model file. Returns
a character vector with the paths to the
model file and compiled executable.
\item \code{x_lines_y}: read the Stan model file for safe transport to
parallel workers. Omitted if \code{compile = "original"}.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple datasets
by repeatedly running the R expression in the \code{data} argument.
Each dynamic branch returns a batch of Stan data lists that \code{x_y}
supplies to the model.
\item \code{x_y}: dynamic branching target to run MCMC once per dataset.
Each dynamic branch returns a tidy data frames of draws
corresponding to a batch of Stan data from \code{x_data}.
\item \code{x}: combine all branches of \code{x_y} into a single non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame of draws.
}
}
\description{
\code{tar_stan_mcmc_rep_draws()} creates targets
to run MCMC multiple times per model and
save only the draws from each run.
}
\details{
Draws could take up a lot of storage. If storage becomes
excessive, please consider thinning the draws or using
\code{tar_stan_mcmc_rep_summaries()} instead.

Most of the arguments are passed to the \verb{$compile()}
and \verb{$sample()} methods of the \code{CmdStanModel} class. If you
previously compiled the model in an upstream \code{\link[=tar_stan_compile]{tar_stan_compile()}}
target, then the model should not recompile.
}
\section{Target objects}{

Most \code{stantargets} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (Sys.getenv("TAR_LONG_EXAMPLES") == "true") {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(stantargets)
# Do not use temporary storage for stan files in real projects
# or else your targets will always rerun.
path <- tempfile(pattern = "", fileext = ".stan")
tar_stan_example_file(path = path)
list(
  tar_stan_mcmc_rep_draws(
    your_model,
    stan_files = path,
    data = tar_stan_example_data(),
    batches = 2,
    reps = 2,
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
\seealso{
Other MCMC: 
\code{\link{tar_stan_mcmc_rep_diagnostics}()},
\code{\link{tar_stan_mcmc_rep_summary}()},
\code{\link{tar_stan_mcmc}()}
}
\concept{MCMC}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_summary.R
\name{tar_stan_summary_join_data}
\alias{tar_stan_summary_join_data}
\title{Join some Stan data to summary output}
\usage{
tar_stan_summary_join_data(summaries, data)
}
\arguments{
\item{summaries}{A data frame of Stan posterior summaries.}

\item{data}{Code to generate the \code{data} for the Stan model.}
}
\value{
A data frame of user-friendly Stan output.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_stan_gq_rep.R
\name{tar_stan_gq_rep}
\alias{tar_stan_gq_rep}
\title{Multiple runs of generated quantities per model with
tidy output}
\usage{
tar_stan_gq_rep(
  name,
  stan_files,
  data = quote(list()),
  fitted_params,
  output_type = c("summary", "draws"),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  compile = c("original", "copy"),
  quiet = TRUE,
  stdout = NULL,
  stderr = NULL,
  dir = NULL,
  pedantic = FALSE,
  include_paths = NULL,
  cpp_options = list(),
  stanc_options = list(),
  force_recompile = FALSE,
  seed = NULL,
  output_dir = NULL,
  sig_figs = NULL,
  parallel_chains = getOption("mc.cores", 1),
  threads_per_chain = NULL,
  data_copy = character(0),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name for the collection of targets.
Serves as a prefix for target names.}

\item{stan_files}{Character vector of paths to known existing Stan model
files created before running the pipeline.}

\item{data}{(multiple options) The data to use for the variables specified in
the data block of the Stan program. One of the following:
\itemize{
\item A named list of \R objects with the names corresponding to variables
declared in the data block of the Stan program. Internally this list is then
written to JSON for CmdStan using \code{\link[cmdstanr:write_stan_json]{write_stan_json()}}. See
\code{\link[cmdstanr:write_stan_json]{write_stan_json()}} for details on the conversions performed on \R objects
before they are passed to Stan.
\item A path to a data file compatible with CmdStan (JSON or \R dump). See the
appendices in the CmdStan manual for details on using these formats.
\item \code{NULL} or an empty list if the Stan program has no data block.
}}

\item{fitted_params}{(multiple options) The parameter draws to use. One of
the following:
\itemize{
\item A \link[cmdstanr]{CmdStanMCMC} or \link[cmdstanr]{CmdStanVB} fitted model object.
\item A \link[posterior:draws_array]{posterior::draws_array} (for MCMC) or \link[posterior:draws_matrix]{posterior::draws_matrix} (for
VB) object returned by CmdStanR's \code{\link[cmdstanr:fit-method-draws]{$draws()}} method.
\item A character vector of paths to CmdStan CSV output files.
}}

\item{output_type}{Type of output to create, either \code{"summaries"},
\code{"draws"}, or \code{"diagnostics"}.}

\item{batches}{Number of batches. Each batch is a sequence
of branch targets containing multiple reps. Each rep
generates a dataset and runs the model on it.}

\item{reps}{Number of replications per batch.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{compile}{(logical) Do compilation? The default is \code{TRUE}. If \code{FALSE}
compilation can be done later via the \code{\link[cmdstanr:model-method-compile]{$compile()}}
method.}

\item{quiet}{(logical) Should the verbose output from CmdStan during
compilation be suppressed? The default is \code{TRUE}, but if you encounter an
error we recommend trying again with \code{quiet=FALSE} to see more of the
output.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{dir}{(string) The path to the directory in which to store the CmdStan
executable (or \code{.hpp} file if using \verb{$save_hpp_file()}). The default is the
same location as the Stan program.}

\item{pedantic}{(logical) Should pedantic mode be turned on? The default is
\code{FALSE}. Pedantic mode attempts to warn you about potential issues in your
Stan program beyond syntax errors. For details see the \href{https://mc-stan.org/docs/reference-manual/pedantic-mode.html}{\emph{Pedantic mode} chapter} in
the Stan Reference Manual. \strong{Note:} to do a pedantic check for a model
that is already compiled use the
\code{\link[cmdstanr:model-method-check_syntax]{$check_syntax()}} method instead.}

\item{include_paths}{(character vector) Paths to directories where Stan
should look for files specified in \verb{#include} directives in the Stan
program.}

\item{cpp_options}{(list) Any makefile options to be used when compiling the
model (\code{STAN_THREADS}, \code{STAN_MPI}, \code{STAN_OPENCL}, etc.). Anything you would
otherwise write in the \code{make/local} file.}

\item{stanc_options}{(list) Any Stan-to-C++ transpiler options to be used
when compiling the model. See the \strong{Examples} section below as well as the
\code{stanc} chapter of the CmdStan Guide for more details on available options:
https://mc-stan.org/docs/cmdstan-guide/stanc.html.}

\item{force_recompile}{(logical) Should the model be recompiled even if was
not modified since last compiled. The default is \code{FALSE}.}

\item{seed}{(positive integer(s)) A seed for the (P)RNG to pass to CmdStan.
In the case of multi-chain sampling the single \code{seed} will automatically be
augmented by the the run (chain) ID so that each chain uses a different
seed. The exception is the transformed data block, which defaults to
using same seed for all chains so that the same data is generated for all
chains if RNG functions are used. The only time \code{seed} should be specified
as a vector (one element per chain) is if RNG functions are used in
transformed data and the goal is to generate \emph{different} data for each
chain.}

\item{output_dir}{(string) A path to a directory where CmdStan should write
its output CSV files. For interactive use this can typically be left at
\code{NULL} (temporary directory) since CmdStanR makes the CmdStan output
(posterior draws and diagnostics) available in \R via methods of the fitted
model objects. The behavior of \code{output_dir} is as follows:
\itemize{
\item If \code{NULL} (the default), then the CSV files are written to a temporary
directory and only saved permanently if the user calls one of the \verb{$save_*}
methods of the fitted model object (e.g.,
\code{\link[cmdstanr:fit-method-save_output_files]{$save_output_files()}}). These temporary
files are removed when the fitted model object is
\link[base:gc]{garbage collected} (manually or automatically).
\item If a path, then the files are created in \code{output_dir} with names
corresponding to the defaults used by \verb{$save_output_files()}.
}}

\item{sig_figs}{(positive integer) The number of significant figures used
when storing the output values. By default, CmdStan represent the output
values with 6 significant figures. The upper limit for \code{sig_figs} is 18.
Increasing this value will result in larger output CSV files and thus an
increased usage of disk space.}

\item{parallel_chains}{(positive integer) The \emph{maximum} number of MCMC chains
to run in parallel. If \code{parallel_chains} is not specified then the default
is to look for the option \code{"mc.cores"}, which can be set for an entire \R
session by \code{options(mc.cores=value)}. If the \code{"mc.cores"} option has not
been set then the default is \code{1}.}

\item{threads_per_chain}{(positive integer) If the model was
\link[cmdstanr:model-method-compile]{compiled} with threading support, the number of
threads to use in parallelized sections \emph{within} an MCMC chain (e.g., when
using the Stan functions \code{reduce_sum()} or \code{map_rect()}). This is in
contrast with \code{parallel_chains}, which specifies the number of chains to
run in parallel. The actual number of CPU cores used use is
\code{parallel_chains*threads_per_chain}. For an example of using threading see
the Stan case study \href{https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html}{Reduce Sum: A Minimal Example}.}

\item{data_copy}{Character vector of names of scalars in \code{data}.
These values will be inserted as columns in the output data frame
for each rep. To join more than just scalars, include a \code{.join_data}
element of your Stan data list with names and dimensions corresponding
to those of the model. For details, read
\url{https://docs.ropensci.org/stantargets/articles/simulation.html}.}

\item{variables}{(character vector) Optionally, the names of the variables
(parameters, transformed parameters, and generated quantities) to read in.
\itemize{
\item If \code{NULL} (the default) then all variables are included.
\item If an empty string (\code{variables=""}) then none are included.
\item For non-scalar variables all elements or specific elements can be selected:
\itemize{
\item \code{variables = "theta"} selects all elements of \code{theta};
\item \code{variables = c("theta[1]", "theta[3]")} selects only the 1st and 3rd elements.
}
}}

\item{summaries}{Optional list of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{Optional list of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frame
of posterior summaries. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{format_df}{Character of length 1, storage format of the data frame
targets such as posterior draws. We recommend efficient data frame formats
such as \code{"feather"} or \code{"aws_parquet"}. For more on storage formats,
see the help file of \code{targets::tar_target()}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Developers can consult the design specification at
\url{https://books.ropensci.org/targets-design/}
to learn about the structure and composition of target objects.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
