
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
