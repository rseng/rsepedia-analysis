
# jagstargets <img src='man/figures/logo.png' align="right" height="139"/>

[![JOSS](https://joss.theoj.org/papers/759f48d9ae7bc57e318e2d0ecc00569e/status.svg)](https://joss.theoj.org/papers/759f48d9ae7bc57e318e2d0ecc00569e)
[![ropensci](https://badges.ropensci.org/425_status.svg)](https://github.com/ropensci/software-review/issues/425)
[![DOI](https://zenodo.org/badge/321076424.svg)](https://zenodo.org/badge/latestdoi/321076424)
[![R
Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![cran](https://www.r-pkg.org/badges/version/jagstargets)](https://cran.r-project.org/package=jagstargets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/jagstargets/workflows/check/badge.svg)](https://github.com/ropensci/jagstargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/jagstargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/gittargets)
[![lint](https://github.com/ropensci/jagstargets/workflows/lint/badge.svg)](https://github.com/ropensci/jagstargets/actions?query=workflow%3Alint)

Bayesian data analysis usually incurs long runtimes and cumbersome
custom code, and the process of prototyping and deploying custom
[JAGS](https://mcmc-jags.sourceforge.io) models can become a daunting
software engineering challenge. To ease this burden, the `jagstargets` R
package creates [JAGS](https://mcmc-jags.sourceforge.io) pipelines that
are concise, efficient, scalable, and tailored to the needs of Bayesian
statisticians. Leveraging
[`targets`](https://docs.ropensci.org/targets/), `jagstargets` pipelines
automatically parallelize the computation and skip expensive steps when
the results are already up to date. Minimal custom user-side code is
required, and there is no need to manually configure branching, so
`jagstargets` is easier to use than
[`targets`](https://docs.ropensci.org/targets/) and
[`R2jags`](https://CRAN.R-project.org/package=R2jags) directly.

## Prerequisites

1.  The [prerequisites of the `targets` R
    package](https://docs.ropensci.org/targets/#prerequisites).
2.  Basic familiarity with
    [`targets`](https://docs.ropensci.org/targets/): watch minutes 7
    through 40 of [this video](https://youtu.be/Gqn7Xn4d5NI?t=439), then
    read [this
    chapter](https://books.ropensci.org/targets/walkthrough.html) of the
    [user manual](https://books.ropensci.org/targets/).
3.  Familiarity with Bayesian Statistics and
    [JAGS](https://mcmc-jags.sourceforge.io/). Prior knowledge of
    [`rjags`](https://cran.r-project.org/package=rjags) or
    [`R2jags`](https://cran.r-project.org/package=R2jags) helps.

## How to get started

Read the `jagstargets` [introductory
vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html),
and then use <https://docs.ropensci.org/jagstargets/> as a reference
while constructing your own workflows. If you need to analyze large
collections of simulated datasets, please consult the [simulation
vignette](https://docs.ropensci.org/jagstargets/articles/simulation.html).

## Installation

`jagstargets` requires the user to install
[JAGS](https://mcmc-jags.sourceforge.io/),
[`rjags`](https://CRAN.R-project.org/package=rjags), and
[`R2jags`](https://CRAN.R-project.org/package=R2jags) beforehand. You
can install JAGS from <https://mcmc-jags.sourceforge.io/>, and you can
install the rest from CRAN.

``` r
install.packages(c("rjags", "R2jags"))
```

Then, install the latest release from CRAN.

``` r
install.packages("jagstargets")
```

Alternatively, install the GitHub development version to access the
latest features and patches.

``` r
install.packages("remotes")
remotes::install_github("ropensci/jagstargets")
```

## Usage

Begin with one or more models: for example, the simple regression model
below with response variable *y* and covariate *x*.

<center>
<img src="./man/figures/model.gif">
</center>

Next, write a JAGS model file for each model like the `model.jags` file
below.

``` jags
model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * beta, 1)
  }
  beta ~ dnorm(0, 1)
}
```

To begin a reproducible analysis pipeline with this model, write a
[`_targets.R` file](https://books.ropensci.org/targets/walkthrough.html)
that loads your packages, defines a function to generate JAGS data, and
lists a pipeline of targets. The target list can call target factories
like
[`tar_jags()`](https://docs.ropensci.org/jagstargets/reference/tar_jags.html)
as well as ordinary targets with
[`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html).
The following minimal example is simple enough to contain entirely
within the `_targets.R` file, but for larger projects, you may wish to
store functions in separate files as in the
[`targets-stan`](https://github.com/wlandau/targets-stan) example.

``` r
# _targets.R
library(targets)
library(jagstargets)

generate_data <- function() {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  out <- list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_jags(
    example,
    jags_files = "model.jags", # You provide this file.
    parameters.to.save = "beta",
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
e.g. `tar_read(tar_read(example_summary_x)`. Visit the [introductory
vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html)
to read more about this example.

## How the package works

`jagstargets` supports specialized [target
factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories)
that create ensembles of [target
objects](https://docs.ropensci.org/targets/reference/tar_target.html)
for [`R2jags`](https://CRAN.R-project.org/package=R2jags) workflows.
These [target
factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories)
abstract away the details of
[`targets`](https://docs.ropensci.org/targets/) and
[`R2jags`](https://CRAN.R-project.org/package=R2jags) and make both
packages easier to use. For details, please read the [introductory
vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html).

## Help

If you have trouble using `jagstargets`, you can ask for help in the
[GitHub discussions
forum](https://github.com/ropensci/jagstargets/discussions/categories/help).
Because the purpose of `jagstargets` is to combine
[`targets`](https://docs.ropensci.org/targets/) and
[`R2jags`](https://CRAN.R-project.org/package=R2jags), your issue may
have something to do with one of the latter two packages, a [dependency
of
`targets`](https://github.com/ropensci/targets/blob/4e3ef2a6c986f558a25e544416f480fc01236b6b/DESCRIPTION#L49-L88),
or [`R2jags`](https://CRAN.R-project.org/package=R2jags) itself. When
you troubleshoot, peel back as many layers as possible to isolate the
problem. For example, if the issue comes from
[`R2jags`](https://CRAN.R-project.org/package=R2jags), create a
[reproducible example](https://reprex.tidyverse.org) that directly
invokes [`R2jags`](https://CRAN.R-project.org/package=R2jags) without
invoking `jagstargets`. The GitHub discussion and issue forums of those
packages are great resources.

## Participation

Development is a community effort, and we welcome discussion and
contribution. By participating in this project, you agree to abide by
the [code of conduct](https://ropensci.org/code-of-conduct/) and the
[contributing
guide](https://github.com/ropensci/jagstargets/blob/main/CONTRIBUTING.md).

## Citation

``` r
citation("jagstargets")
#> 
#> To cite jagstargets in publications use:
#> 
#>   Landau, W. M., (2021). The jagstargets R package: a reproducible
#>   workflow framework for Bayesian data analysis with JAGS. Journal of
#>   Open Source Software, 6(68), 3877,
#>   https://doi.org/10.21105/joss.03877
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {The jagstargets R package: a reproducible workflow framework for Bayesian data analysis with JAGS},
#>     author = {William Michael Landau},
#>     journal = {Journal of Open Source Software},
#>     year = {2021},
#>     volume = {6},
#>     number = {68},
#>     pages = {3877},
#>     url = {https://doi.org/10.21105/joss.03877},
#>   }
```
# jagstargets 1.0.1

* Adjust logo size for README.

# jagstargets 1.0.1

* Reference JOSS paper.

# jagstargets 1.0.0

* Add Zenodo and badge.
* Fix bibliography capitalization.

# jagstargets 0.0.1

* Complete rOpenSci review.
* Allow installation/checks to pass without `rjags` or `R2jags` installed (#18, @jeroen).

# jagstargets 0.0.0.9002

* Replace the `log` argument with `stdout` and `stderr`.
* Switch meaning of `%||%` and `%|||%` to conform to historical precedent.

# jagstargets 0.0.0.9001

* Use `.join_data` list element instead of arguments `copy_data` and `omit_data`.

# jagstargets 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributing

Development is a community effort, and we welcome participation.

## Code of Conduct

By participating in this project, you agree to abide by the [code of conduct](https://ropensci.org/code-of-conduct/).

## Discussions

At <https://github.com/ropensci/jagstargets/discussions>, you can post general questions, brainstorm ideas, and ask for help.

## Issues

<https://github.com/ropensci/jagstargets/issues> is for bug reports, performance issues, maintenance tasks, and feature requests. When you post, please abide by the following guidelines.

* Before posting a new issue, please take a moment to search for existing similar issues in order to avoid duplication.
* For bug reports: if you can, please install the latest GitHub version of `jagstargets` (i.e. `remotes::install_github("ropensci/jagstargets")`) and verify that the issue still persists.
* Describe your issue in prose as clearly and concisely as possible.
* For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Development

External code contributions are extremely helpful in the right circumstances. Here are the recommended steps.

1. Prior to contribution, please propose your idea in a [new issue thread](https://github.com/ropensci/jagstargets/issues) so you and the maintainer can define the intent and scope of your work.
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
* Describe your contribution in the project's [`NEWS.md`](https://github.com/ropensci/jagstargets/blob/main/NEWS.md) file. Be sure to mention relevant GitHub issue numbers and your GitHub name as done in existing news entries.
* If you feel contribution is substantial enough for official author or contributor status, please add yourself to the `Authors@R` field of the [`DESCRIPTION`](https://github.com/ropensci/jagstargets/blob/main/DESCRIPTION) file.
# Prework

* [ ] I understand and agree to this repository's [code of conduct](https://ropensci.org/code-of-conduct/).
* [ ] I understand and agree to this repository's [contributing guidelines](https://github.com/ropensci/jagstargets/blob/main/CONTRIBUTING.md).
* [ ] I have already submitted an issue to the [issue tracker](http://github.com/ropensci/jagstargets/issues) to discuss my idea with the maintainer.

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
about: Something is wrong with jagstargets.
title: ""
labels: "type: bug"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/jagstargets/blob/main/CONTRIBUTING.md).
* [ ] Confirm that your issue is most likely a genuine bug in `jagstargets` and not a known limitation or usage error.
* [ ] If there is [already a relevant issue](https://github.com/ropensci/jagstargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
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
* The [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity) of the GitHub commit of `jagstargets` currently installed. `packageDescription("jagstargets")$GithubSHA1` shows you this.
---
name: New feature
about: Suggest a new feature.
title: ""
labels: "type: new feature"
assignees: wlandau
---

## Prework

* [ ] I understand and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and  [contributing guidelines](https://github.com/ropensci/jagstargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/jagstargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
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

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com//jagstargets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com//jagstargets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com//targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the performance issue.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com//targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Benchmarks

How poorly does `jagstargets` perform? To find out, we recommend you use the [`proffer`](https://github.com/r-prof/proffer) package and take screenshots of the results displayed in your browser.

```r
library(jagstargets)
library(proffer)
px <- pprof({
  # All your jagstargets code goes here.
})
```
---
title: 'The jagstargets R package: a reproducible workflow framework for Bayesian data analysis with JAGS'
tags:
- R
- reproducibility
- high-performance computing
- pipeline
- workflow
- Make
- Bayesian
- JAGS
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

Researchers who perform Bayesian statistics regularly experiment with models to refine, compare, and understand them [@bayesworkflow]. Each fitted model informs subsequent model-building choices, and even for experienced practitioners, the investigation often leads to final models that differ from the ones originally proposed. Fitting a model usually means applying Markov chain Monte Carlo (MCMC) or a similar method to approximate the full joint posterior distribution of the parameters [@bda3]. Flexible probabilistic programming languages such as JAGS have made model specification quick and straightforward [@jags], but computation time is still a bottleneck. A workflow can take several minutes or hours to run, and in subsequent iterations, researchers struggle to keep the results up to date with frequent changes to the models, code, and data.

# Summary

The [`jagstargets`](https://docs.ropensci.org/jagstargets/) R package [@jagstargets] is a workflow toolkit for Bayesian data analysis with JAGS. It helps users express a Bayesian statistical modeling exercise as a formal pipeline with dedicated steps for data generation, analysis, and summarization. Pipelines can be customized for popular use cases: for example, the analysis of a single dataset using multiple alternative models, or a large simulation study to test that a model is implemented correctly [@carpenter2017]. These pipelines, which can be visualized and executed using [`targets`](https://docs.ropensci.org/targets/) [@targets], support key features that increase the efficiency and reproducibility of Bayesian workflows. The steps of a pipeline are automatically orchestrated using optional distributed computing, and up-to-date tasks are automatically skipped if the upstream code and data did not change since the last run. Thus, researchers can quickly iterate on Bayesian workflows while maintaining agreement between the results and the underlying code, models, and datasets.

The [`jagstargets`](https://docs.ropensci.org/jagstargets/) package is a bridge between packages [`R2jags`](https://github.com/suyusung/R2jags) [@R2jags] and [`targets`](https://docs.ropensci.org/targets/). [`targets`](https://docs.ropensci.org/targets/) is a general-purpose pipeline toolkit for reproducible research and high-performance computing. It can be used on its own for Bayesian data analysis, but it does not natively support  domain-specific functionality for Bayesian workflows or JAGS. [`jagstargets`](https://docs.ropensci.org/jagstargets/), on the other hand, leverages the existing capabilities of [`targets`](https://docs.ropensci.org/targets/) to more easily and concisely express computational pipelines for Bayesian data analysis, from single analyses or large-scale simulation studies. It abstracts away laborious error-prone tasks such as invoking JAGS, tidying posterior samples, and batching thousands of MCMC runs to improve computational efficiency, all of which reduce the user-side burden of writing custom R code for a [`targets`](https://docs.ropensci.org/targets/) pipeline. [`jagstargets`](https://docs.ropensci.org/targets/) is similar to the [`stantargets`](https://docs.ropensci.org/stantargets/) R package [@stantargets], the latter of which streamlines pipeline construction of Bayesian data analysis pipelines with Stan [@stan].

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

# jagstargets <img src='man/figures/logo.png' align="right" height="139"/>

[![JOSS](https://joss.theoj.org/papers/759f48d9ae7bc57e318e2d0ecc00569e/status.svg)](https://joss.theoj.org/papers/759f48d9ae7bc57e318e2d0ecc00569e)
[![ropensci](https://badges.ropensci.org/425_status.svg)](https://github.com/ropensci/software-review/issues/425)
[![DOI](https://zenodo.org/badge/321076424.svg)](https://zenodo.org/badge/latestdoi/321076424)
[![R Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![cran](https://www.r-pkg.org/badges/version/jagstargets)](https://cran.r-project.org/package=jagstargets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/jagstargets/workflows/check/badge.svg)](https://github.com/ropensci/jagstargets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/jagstargets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/gittargets)
[![lint](https://github.com/ropensci/jagstargets/workflows/lint/badge.svg)](https://github.com/ropensci/jagstargets/actions?query=workflow%3Alint)

Bayesian data analysis usually incurs long runtimes and cumbersome custom code, and the process of prototyping and deploying custom [JAGS](https://mcmc-jags.sourceforge.io) models can become a daunting software engineering challenge. To ease this burden, the `jagstargets` R package creates [JAGS](https://mcmc-jags.sourceforge.io) pipelines that are concise, efficient, scalable, and tailored to the needs of Bayesian statisticians. Leveraging [`targets`](https://docs.ropensci.org/targets/), `jagstargets` pipelines automatically parallelize the computation and skip expensive steps when the results are already up to date. Minimal custom user-side code is required, and there is no need to manually configure branching, so `jagstargets` is easier to use than [`targets`](https://docs.ropensci.org/targets/) and [`R2jags`](https://CRAN.R-project.org/package=R2jags) directly.

## Prerequisites

1. The [prerequisites of the `targets` R package](https://docs.ropensci.org/targets/#prerequisites).
1. Basic familiarity with [`targets`](https://docs.ropensci.org/targets/): watch minutes 7 through 40 of [this video](https://youtu.be/Gqn7Xn4d5NI?t=439), then read [this chapter](https://books.ropensci.org/targets/walkthrough.html) of the [user manual](https://books.ropensci.org/targets/).
1. Familiarity with Bayesian Statistics and [JAGS](https://mcmc-jags.sourceforge.io/). Prior knowledge of [`rjags`](https://cran.r-project.org/package=rjags) or [`R2jags`](https://cran.r-project.org/package=R2jags) helps.

## How to get started

Read the `jagstargets` [introductory vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html), and then use <https://docs.ropensci.org/jagstargets/> as a reference while constructing your own workflows. If you need to analyze large collections of simulated datasets, please consult the [simulation vignette](https://docs.ropensci.org/jagstargets/articles/simulation.html).

## Installation

`jagstargets` requires the user to install [JAGS](https://mcmc-jags.sourceforge.io/), [`rjags`](https://CRAN.R-project.org/package=rjags), and [`R2jags`](https://CRAN.R-project.org/package=R2jags) beforehand. You can install JAGS from <https://mcmc-jags.sourceforge.io/>, and you can install the rest from CRAN.

```{r, eval = FALSE}
install.packages(c("rjags", "R2jags"))
```

Then, install the latest release from CRAN.

```{r, eval = FALSE}
install.packages("jagstargets")
```

Alternatively, install the GitHub development version to access the latest features and patches.

```{r, eval = FALSE}
install.packages("remotes")
remotes::install_github("ropensci/jagstargets")
```

## Usage

Begin with one or more models: for example, the simple regression model below with response variable $y$ and covariate $x$.

<center>
<img src="./man/figures/model.gif">
</center>

Next, write a JAGS model file for each model like the `model.jags` file below.

```jags
model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * beta, 1)
  }
  beta ~ dnorm(0, 1)
}
```

To begin a reproducible analysis pipeline with this model, write a [`_targets.R` file](https://books.ropensci.org/targets/walkthrough.html) that loads your packages, defines a function to generate JAGS data, and lists a pipeline of targets. The target list can call target factories like [`tar_jags()`](https://docs.ropensci.org/jagstargets/reference/tar_jags.html) as well as ordinary targets with [`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html). The following minimal example is simple enough to contain entirely within the `_targets.R` file, but for larger projects, you may wish to store functions in separate files as in the [`targets-stan`](https://github.com/wlandau/targets-stan) example.

```{r, eval = FALSE}
# _targets.R
library(targets)
library(jagstargets)

generate_data <- function() {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  out <- list(n = n, x = x, y = y, true_beta = true_beta)
}

list(
  tar_jags(
    example,
    jags_files = "model.jags", # You provide this file.
    parameters.to.save = "beta",
    data = generate_data()
  )
)
```

Run [`tar_visnetwork()`](https://docs.ropensci.org/targets/reference/tar_visnetwork.html) to check `_targets.R` for correctness, then call [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html) to run the pipeline. Access the results using [`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html), e.g. `tar_read(tar_read(example_summary_x)`. Visit the [introductory vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html) to read more about this example.

## How the package works

`jagstargets` supports specialized [target factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories) that create ensembles of [target objects](https://docs.ropensci.org/targets/reference/tar_target.html) for [`R2jags`](https://CRAN.R-project.org/package=R2jags) workflows. These [target factories](https://ropensci.org/blog/2021/02/03/targets/#target-factories) abstract away the details of [`targets`](https://docs.ropensci.org/targets/) and [`R2jags`](https://CRAN.R-project.org/package=R2jags) and make both packages easier to use. For details, please read the [introductory vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html).

## Help

If you have trouble using `jagstargets`, you can ask for help in the [GitHub discussions forum](https://github.com/ropensci/jagstargets/discussions/categories/help). Because the purpose of `jagstargets` is to combine [`targets`](https://docs.ropensci.org/targets/) and [`R2jags`](https://CRAN.R-project.org/package=R2jags), your issue may have something to do with one of the latter two packages, a [dependency of `targets`](https://github.com/ropensci/targets/blob/4e3ef2a6c986f558a25e544416f480fc01236b6b/DESCRIPTION#L49-L88), or [`R2jags`](https://CRAN.R-project.org/package=R2jags) itself. When you troubleshoot, peel back as many layers as possible to isolate the problem. For example, if the issue comes from [`R2jags`](https://CRAN.R-project.org/package=R2jags), create a [reproducible example](https://reprex.tidyverse.org) that directly invokes [`R2jags`](https://CRAN.R-project.org/package=R2jags) without invoking `jagstargets`. The GitHub discussion and issue forums of those packages are great resources.

## Participation

Development is a community effort, and we welcome discussion and contribution. By participating in this project, you agree to abide by the [code of conduct](https://ropensci.org/code-of-conduct/) and the [contributing guide](https://github.com/ropensci/jagstargets/blob/main/CONTRIBUTING.md).

## Citation

```{r, warning = FALSE}
citation("jagstargets")
```
---
title: "Introduction to jagstargets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to jagstargets}
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
skip <- identical(Sys.getenv("NOT_CRAN", unset = "false"), "false") ||
  !requireNamespace("rjags", quietly = TRUE) ||
  !requireNamespace("R2jags", quietly = TRUE)
if (skip) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(R2jags)
}
library(dplyr)
library(targets)
library(jagstargets)
```

The `jagstargets` package makes it easy to run a single jags model and keep track of the results. [`R2jags`](https://github.com/suyusung/R2jags) fits the models, and [`targets`](https://docs.ropensci.org/targets/) manages the workflow and helps avoid unnecessary computation.

Consider the simple regression model below with response variable `y` and covariate `x`.


$$
\begin{aligned}
y_i &\stackrel{\text{iid}}{\sim} \text{Normal}(x_i \beta, 1) \\
\beta &\sim \text{Normal}(0, 1)
\end{aligned}
$$

We write this model in the JAGS model file below.

```{r}
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * beta, 1)
  }
  beta ~ dnorm(0, 1)
}"
writeLines(lines, "x.jags")
```

A typical workflow proceeds as follows:

1. Prepare a list of input data to JAGS, including vector elements `x` and `y`.
1. Fit the JAGS model using the list of input data.
1. Use the fitted model object to compute posterior summaries and convergence diagnostics.
1. Use the fitted model object to extract posterior draws of parameters and store them in a tidy data frame.
1. If there are other models to compare, use the fitted model object to compute the deviance information criterion (DIC).

`jagstargets` encapsulates this workflow with the [`tar_jags()`](https://docs.ropensci.org/jagstargets/reference/tar_jags.html) function. To use it in a [`targets`](https://docs.ropensci.org/targets/) pipeline, invoke it from the `_targets.R` script of the project.

```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  list(
    tar_jags(
      example,
      jags_files = "x.jags",
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(jagstargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  out <- list(n = n, x = x, y = y)
}

# The _targets.R file ends with a list of target objects
# produced by jagstargets::tar_jags(), targets::tar_target(), or similar.
list(
  tar_jags(
    example,
    jags_files = "x.jags",
    parameters.to.save = "beta",
    data = generate_data()
  )
)
```

[`tar_jags()`](https://docs.ropensci.org/jagstargets/reference/tar_jags.html) only *defines* the pipeline. It does not actually run JAGS, it declares the targets that will eventually run JAGS. The specific targets are as follows. Run `tar_manifest()` to show specific details about the targets declared.

```{r}
tar_manifest()
```

Each target is responsible for a piece of the workflow.

* `example_file_x`: Reproducibly track changes to the jags model file.
* `example_data`: Run the code you supplied to the `data` argument of `tar_jags()` and return a dataset compatible with JAGS.
* `example_mcmc_x`: Run the MCMC and return an object of class `rjags` from [`R2jags`](https://github.com/suyusung/R2jags).
* `example_draws_x`: Return a friendly `tibble` of the posterior draws from `example`.
* `example_summaries_x`: Return a friendly `tibble` of the posterior summaries from `example`. Uses [`posterior::summarize_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
* `example_dic_x`: Return a friendly `tibble` with each model's DIC and penalty.

The suffix `_x` comes from the base name of the model file, in this case `x.jags`. If you supply multiple model files to the `jags_files` argument, all the models share the same dataset, and the suffixes distinguish among the various targets.

The targets depend on one another: for example, `example_mcmc_x` takes `example_data` as input. [`targets`](https://docs.ropensci.org/targets/) can visualize the dependency relationships in a dependency graph, which is helpful for understanding the pipeline and troubleshooting issues.

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

But if we change the underlying code or data, some of the targets will no longer be valid, and they will rerun during the next `tar_make()`. Below, we change the jags model file, so the MCMC reruns while the data is skipped. This behavior saves time and enhances reproducibility.

```{r}
write(" ", file = "x.jags", append = TRUE)
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

At this point, we can add more targets and custom functions for additional post-processing. See below for a custom summary target (which is equivalent to customizing the `summaries` argument of `tar_jags()`.)

```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  list(
    tar_jags(
      example,
      jags_files = "x.jags",
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
    ),
    tar_target(
      custom_summary,
      posterior::summarize_draws(
        dplyr::select(example_draws_x, -.draw),
        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
      )
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(jagstargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  out <- list(n = n, x = x, y = y)
}

list(
  tar_jags(
    example,
    jags_files = "x.jags",
    parameters.to.save = "beta",
    data = generate_data()
  ),
  tar_target(
    custom_summary,
    posterior::summarize_draws(
      dplyr::select(example_draws_x, -.draw),
      ~posterior::quantile2(.x, probs = c(0.25, 0.75))
    )
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

`tar_jags()` and related functions allow you to supply multiple models to `jags_files`. If you do, each model will run on the same dataset. Consider a new model, `y.jags`.

```{r}
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * x[i] * beta, 1) # Regress on x^2 instead of x.
  }
  beta ~ dnorm(0, 1)
}"
writeLines(lines, "y.jags")
```

Below, we add `y.jags` to the `jags_files` argument of `tar_jags()`.

```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
tar_script({
  library(targets)
  library(jagstargets)
  list(
    tar_jags(
      example,
      jags_files = c("x.jags", "y.jags"),
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
    ),
    tar_target(
      custom_summary,
      posterior::summarize_draws(
        dplyr::select(example_draws_x, -.draw),
        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
      )
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(jagstargets)

generate_data <- function(n = 10) {
  true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, x * true_beta, 1)
  out <- list(n = n, x = x, y = y)
}

list(
  tar_jags(
    example,
    jags_files = c("x.jags", "y.jags"),
    parameters.to.save = "beta",
    data = generate_data()
  ),
  tar_target(
    custom_summary,
    posterior::summarize_draws(
      dplyr::select(example_draws_x, -.draw),
      ~posterior::quantile2(.x, probs = c(0.25, 0.75))
    )
  )
)
```

In the graph below, notice how the `*_x` targets and `*_y` targets are both connected to `example_data` upstream.

```{r}
tar_visnetwork(targets_only = TRUE)
```

## More information

For more on [`targets`](https://docs.ropensci.org/targets/), please visit the reference website <https://docs.ropensci.org/targets/> or the user manual <https://books.ropensci.org/targets/>. The manual walks though advanced features of `targets` such as [high-performance computing](https://books.ropensci.org/targets/hpc.html) and [cloud storage support](https://books.ropensci.org/targets/cloud.html).
---
title: "Bayesian simulation pipelines with jagstargets"
output: rmarkdown::html_vignette
bibliography: simulation.bib
vignette: >
  %\VignetteIndexEntry{Bayesian simulation pipelines with jagstargets}
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
skip <- identical(Sys.getenv("NOT_CRAN", unset = "false"), "false") ||
  !requireNamespace("rjags", quietly = TRUE) ||
  !requireNamespace("R2jags", quietly = TRUE)
if (skip) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(R2jags)
}
library(dplyr)
library(targets)
library(jagstargets)
```

The [introductory vignette](https://docs.ropensci.org/jagstargets/articles/introduction.html) vignette caters to Bayesian data analysis workflows with few datasets to analyze. However, it is sometimes desirable to run one or more Bayesian models repeatedly across many simulated datasets. Examples:

1. Validate the implementation of a Bayesian model, using simulation to determine how reliably the model estimates the parameters under known data-generating scenarios.
2. Simulate a randomized controlled experiment to explore frequentist properties such as power and Type I error.

This vignette focuses on (1). The goal of this particular example to simulate multiple datasets from the model below, analyze each dataset, and assess how often the estimated posterior intervals cover the true parameters from the prior predictive simulations. The quantile method by @cook2006 generalizes this concept, and simulation-based calibration [@talts2020] generalizes further. The interval-based technique featured in this vignette is not as robust as SBC, but it may be more expedient for large models because it does not require visual inspection of multiple histograms.

Consider a simple regression model with a continuous response `y` with a covariate `x`.

$$
\begin{aligned}
y_i &\stackrel{\text{iid}}{\sim} \text{Normal}(\beta_1 + x_i \beta_2, 1) \\
\beta_1, \beta_2 &\stackrel{\text{iid}}{\sim} \text{Normal}(0, 1)
\end{aligned}
$$

We write this model in a JAGS model file.

```{r}
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(beta[1] + x[i] * beta[2], 1)
  }
  for (i in 1:2) {
    beta[i] ~ dnorm(0, 1)
  }
}"
writeLines(lines, "model.jags")
```

Next, we define a pipeline to simulate multiple datasets and fit each dataset with the model. In our data-generating function, we put the true parameter values of each simulation in a special `.join_data` list. `jagstargets` will automatically join the elements of `.join_data` to the correspondingly named variables in the summary output. This will make it super easy to check how often our posterior intervals capture the truth. As for scale, generate 20 datasets (5 batches with 4 replications each) and run the model on each of the 20 datasets.^[Internally, each batch is a [dynamic branch target](https://books.ropensci.org/targets/dynamic.html), and the number of replications determines the amount of work done within a branch. In the general case, [batching](https://books.ropensci.org/targets/dynamic.html#batching) is a way to find the right compromise between target-specific overhead and the horizontal scale of the pipeline.] By default, each of the 20 model runs computes 3 MCMC chains with 2000 MCMC iterations each (including burn-in) and you can adjust with the `n.chains` and `n.iter` arguments of `tar_jags_rep_summary()`.

```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function (n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      "model.jags",
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5, # Number of branch targets.
      reps = 4, # Number of model reps per branch target.
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
library(targets)
library(jagstargets)
options(crayon.enabled = FALSE)
# Use computer memory more sparingly:
tar_option_set(memory = "transient", garbage_collection = TRUE)

generate_data <- function (n = 10L) {
  beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
  # Elements of .join_data get joined on to the .join_data column
  # in the summary output next to the model parameters
  # with the same names.
  .join_data <- list(beta = beta)
  list(n = n, x = x, y = y, .join_data = .join_data)
}

list(
  tar_jags_rep_summary(
    model,
    "model.jags",
    data = generate_data(),
    parameters.to.save = "beta",
    batches = 5, # Number of branch targets.
    reps = 4, # Number of model reps per branch target.
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    )
  )
)
```

We now have a pipeline that runs the model 10 times: 5 batches (branch targets) with 4 replications per batch.

```{r}
tar_visnetwork()
```

Run the computation with `tar_make()`

```{r, output = FALSE, warning = FALSE}
tar_make()
```

The result is an aggregated data frame of summary statistics, where the `.rep` column distinguishes among individual replicates. We have the posterior intervals for `beta` in columns `q2.5` and `q97.5`. And thanks to the `.join_data` list we included in `generate_data()`, our output has a `.join_data` column with the true values of the parameters in our simulations.

```{r}
tar_load(model)
model
```

Now, let's assess how often the estimated 95% posterior intervals capture the true values of `beta`. If the model is implemented correctly, the coverage value below should be close to 95%. (Ordinarily, we would [increase the number of batches and reps per batch](https://books.ropensci.org/targets/dynamic.html#batching) and [run batches in parallel computing](https://books.ropensci.org/targets/hpc.html).)

```{r}
library(dplyr)
model %>%
  group_by(variable) %>%
  dplyr::summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
```

For maximum reproducibility, we should express the coverage assessment as a custom function and a target in the pipeline.

```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  generate_data <- function (n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    # Elements of .join_data get joined on to the .join_data column
    # in the summary output next to the model parameters
    # with the same names.
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      "model.jags",
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5, # Number of branch targets.
      reps = 4, # Number of model reps per branch target.
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    ),
    tar_target(
      coverage,
      model %>%
        group_by(variable) %>%
        summarize(
          coverage = mean(q2.5 < .join_data & .join_data < q97.5),
          .groups = "drop"
        )
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
# packages needed to define the pipeline:
library(targets)
library(jagstargets)

tar_option_set(
  packages = "dplyr", # packages needed to run the pipeline
  memory = "transient", # memory efficiency
  garbage_collection = TRUE # memory efficiency
)

generate_data <- function (n = 10L) {
  beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
  # Elements of .join_data get joined on to the .join_data column
  # in the summary output next to the model parameters
  # with the same names.
  .join_data <- list(beta = beta)
  list(n = n, x = x, y = y, .join_data = .join_data)
}

list(
  tar_jags_rep_summary(
    model,
    "model.jags",
    data = generate_data(),
    parameters.to.save = "beta",
    batches = 5, # Number of branch targets.
    reps = 4, # Number of model reps per branch target.
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    )
  ),
  tar_target(
    coverage,
    model %>%
      group_by(variable) %>%
      summarize(
        coverage = mean(q2.5 < .join_data & .join_data < q97.5),
        .groups = "drop"
      )
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

`tar_jags_rep_mcmc_summary()` and similar functions allow you to supply multiple jags models. If you do, each model will share the the same collection of datasets. Below, we add a new `model2.jags` file to the `jags_files` argument of `tar_jags_rep_mcmc_summary()`. In the coverage summary below, we group by `.name` to compute a coverage statistic for each model.

```{r}
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(beta[1] + x[i] * x[i] * beta[2], 1) # Regress on x^2 instead of x.
  }
  for (i in 1:2) {
    beta[i] ~ dnorm(0, 1)
  }
}"
writeLines(lines, "model2.jags")
```


```{r, echo = FALSE}
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  generate_data <- function (n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    # Elements of .join_data get joined on to the .join_data column
    # in the summary output next to the model parameters
    # with the same names.
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      c("model.jags", "model2.jags"), # another model
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5,
      reps = 4,
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    ),
    tar_target(
      coverage,
      model %>%
        group_by(.name) %>%
        summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
    )
  )
})
```

```{r, eval = FALSE}
# _targets.R
# packages needed to define the pipeline:
library(targets)
library(jagstargets)

tar_option_set(
  packages = "dplyr", # packages needed to run the pipeline
  memory = "transient", # memory efficiency
  garbage_collection = TRUE # memory efficiency
)

generate_data <- function (n = 10L) {
  beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
  x <- seq(from = -1, to = 1, length.out = n)
  y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
  # Elements of .join_data get joined on to the .join_data column
  # in the summary output next to the model parameters
  # with the same names.
  .join_data <- list(beta = beta)
  list(n = n, x = x, y = y, .join_data = .join_data)
}

list(
  tar_jags_rep_summary(
    model,
    c("model.jags", "model2.jags"), # another model
    data = generate_data(),
    parameters.to.save = "beta",
    batches = 5,
    reps = 4,
    variables = "beta",
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.975))
    )
  ),
  tar_target(
    coverage,
    model %>%
      group_by(.name) %>%
      summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
  )
)
```

In the graph below, notice how targets `model_model1` and `model_model2` are both connected to `model_data` upstream. Downstream, `model` is equivalent to `dplyr::bind_rows(model_model1, model_model2)`, and it will have special columns `.name` and `.file` to distinguish among all the models.

```{r}
tar_visnetwork()
```

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_rep_dic.R
\name{tar_jags_rep_dic}
\alias{tar_jags_rep_dic}
\title{Tidy DIC output from multiple MCMCs per model}
\usage{
tar_jags_rep_dic(
  name,
  jags_files,
  parameters.to.save,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  n.cluster = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = as.integer(n.iter/2),
  n.thin = 1,
  jags.module = c("glm", "dic"),
  inits = NULL,
  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister"),
  jags.seed = 123,
  stdout = NULL,
  stderr = NULL,
  progress.bar = "text",
  refresh = 0,
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

\item{jags_files}{Character vector of JAGS model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{tools::file_path_sans_ext(basename(jags_files))} will be used
as target name suffixes. If \code{jags_files} is a named vector,
the suffixed will come from \code{names(jags_files)}.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{Code to generate the \code{data} list for the JAGS model.
Optionally include a \code{.join_data} element to join parts of the data
to correspondingly named parameters in the summary output.
See the vignettes for details.}

\item{batches}{Number of batches. Each batch runs a model \code{reps} times.}

\item{reps}{Number of replications per batch. Ideally, each rep
should produce its own random dataset using the code
supplied to \code{data}.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frames
of posterior summaries and other data frames returned by targets.
We recommend efficient data frame formats
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
\code{tar_jags_rep_dic()} returns list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{jags_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_jags_rep_dic(name = x, jags_files = "y.jags")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the JAGS model file. Returns
a character vector of length 1 with the path to the JAGS
model file.
\item \code{x_lines_y}: read the contents of the JAGS model file
for safe transport to parallel workers.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple JAGS
datasets from the R expression in the \code{data} argument.
Each dynamic branch returns a batch of JAGS data lists.
\item \code{x_y}: run JAGS on each dataset from \code{x_data}.
Each dynamic branch returns a tidy data frame of DIC
results for each batch of data.
\item \code{x}: combine all the batches from \code{x_y} into a non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame with all DIC info
from all the branches of \code{x_y}.
}
}
\description{
Run multiple MCMCs on simulated datasets
and return DIC and the effective number of parameters
for each run.
}
\details{
The MCMC targets use \code{R2jags::jags()} if \code{n.cluster} is \code{1} and
\code{R2jags::jags.parallel()} otherwise. Most arguments to \code{tar_jags()}
are forwarded to these functions.
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
if (requireNamespace("R2jags", quietly = TRUE)) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(jagstargets)
# Do not use a temp file for a real project
# or else your targets will always rerun.
tmp <- tempfile(pattern = "", fileext = ".jags")
tar_jags_example_file(tmp)
list(
  tar_jags_rep_dic(
    your_model,
    jags_files = tmp,
    data = tar_jags_example_data(),
    parameters.to.save = "beta",
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags.R
\name{tar_jags}
\alias{tar_jags}
\title{One MCMC per model with multiple outputs}
\usage{
tar_jags(
  name,
  jags_files,
  parameters.to.save,
  data = list(),
  summaries = list(),
  summary_args = list(),
  n.cluster = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = as.integer(n.iter/2),
  n.thin = 1,
  jags.module = c("glm", "dic"),
  inits = NULL,
  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister"),
  jags.seed = 1,
  stdout = NULL,
  stderr = NULL,
  progress.bar = "text",
  refresh = 0,
  draws = TRUE,
  summary = TRUE,
  dic = TRUE,
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

\item{jags_files}{Character vector of JAGS model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{tools::file_path_sans_ext(basename(jags_files))} will be used
as target name suffixes. If \code{jags_files} is a named vector,
the suffixed will come from \code{names(jags_files)}.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{Code to generate the \code{data} list for the JAGS model.
Optionally include a \code{.join_data} element to join parts of the data
to correspondingly named parameters in the summary output.
See the vignettes for details.}

\item{summaries}{List of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{List of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{draws}{Logical, whether to create a target for posterior draws.
Saves draws as a compressed \code{posterior::as_draws_df()} \code{tibble}.
Convenient, but duplicates storage.}

\item{summary}{Logical, whether to create a target to store a small
data frame of posterior summary statistics and convergence diagnostics.}

\item{dic}{Logical, whether to create a target with deviance
information criterion (DIC) results.}

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
targets such as the JAGS data and any JAGS fit objects.
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
\code{tar_jags()} returns list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{jags_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_jags(name = x, jags_files = "y.jags", ...)} returns a list
of \code{targets::tar_target()} objects:
\itemize{
\item \code{x_file_y}: reproducibly track the JAGS model file. Returns
a character vector of length 1 with the path to the JAGS
model file.
\item \code{x_lines_y}: read the contents of the JAGS model file
for safe transport to parallel workers.
Returns a character vector of lines in the model file.
\item \code{x_data}: run the R expression in the \code{data} argument to produce
a JAGS dataset for the model. Returns a JAGS data list.
\item \code{x_mcmc_y}: run MCMC on the model and dataset.
Returns an \code{rjags} object from \code{R2jags} with all the MCMC results.
\item \code{x_draws_y}: extract posterior samples from \code{x_mcmc_y}.
Returns a tidy data frame of MCMC draws. Omitted if \code{draws = FALSE}.
\item \code{x_summary_y}: extract posterior summaries from \code{x_mcmc_y}.
Returns a tidy data frame of MCMC draws.
Omitted if \code{summary = FALSE}.
\item \code{x_dic}: extract deviance information criterion (DIC) info
from \code{x_mcmc_y}. Returns a tidy data frame of DIC info.
Omitted if \code{dic = FALSE}.
}
}
\description{
Targets to run a JAGS model once with MCMC
and save multiple outputs.
}
\details{
The MCMC targets use \code{R2jags::jags()} if \code{n.cluster} is \code{1} and
\code{R2jags::jags.parallel()} otherwise. Most arguments to \code{tar_jags()}
are forwarded to these functions.
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
if (requireNamespace("R2jags", quietly = TRUE)) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(jagstargets)
# Do not use a temp file for a real project
# or else your targets will always rerun.
tmp <- tempfile(pattern = "", fileext = ".jags")
tar_jags_example_file(tmp)
list(
  tar_jags(
    your_model,
    jags_files = tmp,
    data = tar_jags_example_data(),
    parameters.to.save = "beta",
    stdout = R.utils::nullfile(),
    stderr = R.utils::nullfile()
  )
)
}, ask = FALSE)
targets::tar_make()
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_rep.R
\name{tar_jags_rep}
\alias{tar_jags_rep}
\title{Tidy output from multiple MCMCs per model.}
\usage{
tar_jags_rep(
  name,
  jags_files,
  parameters.to.save,
  data = quote(list()),
  batches = 1L,
  reps = 1L,
  output = c("summary", "draws", "dic"),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  combine = TRUE,
  n.cluster = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = as.integer(n.iter/2),
  n.thin = 1,
  jags.module = c("glm", "dic"),
  inits = NULL,
  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister"),
  jags.seed = 1,
  stdout = NULL,
  stderr = NULL,
  progress.bar = "text",
  refresh = 0,
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

\item{jags_files}{Character vector of JAGS model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{tools::file_path_sans_ext(basename(jags_files))} will be used
as target name suffixes. If \code{jags_files} is a named vector,
the suffixed will come from \code{names(jags_files)}.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{Code to generate the \code{data} list for the JAGS model.
Optionally include a \code{.join_data} element to join parts of the data
to correspondingly named parameters in the summary output.
See the vignettes for details.}

\item{batches}{Number of batches. Each batch runs a model \code{reps} times.}

\item{reps}{Number of replications per batch. Ideally, each rep
should produce its own random dataset using the code
supplied to \code{data}.}

\item{output}{Character of length 1 denoting the type of output \code{tibble}
to return: \code{"draws"} for MCMC samples (which could take up a lot of space)
\code{"summary"} for lightweight posterior summary statistics,
and \code{"dic"} for the overall deviance information criterion
and effective number of parameters}

\item{variables}{Character vector of model parameter names.
The output posterior summaries are restricted to these variables.}

\item{summaries}{List of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{List of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frames
of posterior summaries and other data frames returned by targets.
We recommend efficient data frame formats
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
\code{tar_jags_rep(name = x, jags_files = "y.jags")}
returns a list of \code{targets::tar_target()} objects:
\itemize{
\item \code{x_file_y}: reproducibly track the jags model file.
\item \code{x_lines_y}: contents of the jags model file.
\item \code{x_data}: dynamic branching target with simulated datasets.
\item \code{x_y}: dynamic branching target with tidy data frames of MCMC summaries.
\item \code{x}: combine all the model-specific summary targets into
a single data frame with columns to distinguish among the models.
Suppressed if \code{combine} is \code{FALSE}.
}
}
\description{
Internal function. Users should not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_df.R
\name{tar_jags_df}
\alias{tar_jags_df}
\title{Select a strategic piece of \code{R2jags} output.}
\usage{
tar_jags_df(
  fit,
  data,
  output = c("draws", "summary", "dic"),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL
)
}
\arguments{
\item{fit}{\code{R2jags} object.}

\item{data}{A list, the original JAGS dataset.}

\item{output}{Character of length 1 denoting the type of output \code{tibble}
to return: \code{"draws"} for MCMC samples (which could take up a lot of space)
\code{"summary"} for lightweight posterior summary statistics,
and \code{"dic"} for the overall deviance information criterion
and effective number of parameters}

\item{variables}{Character vector of model parameter names.
The output posterior summaries are restricted to these variables.}

\item{summaries}{List of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{List of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}
}
\value{
A data frame of \code{R2jags} output. Depends on the \code{output} argument.
}
\description{
Not a user-side function. Do not call directly.
Exported for infrastructure reasons only.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_package.R
\docType{package}
\name{jagstargets-package}
\alias{jagstargets-package}
\alias{jagstargets}
\title{jagstargets: Targets for JAGS Workflows}
\description{
Bayesian data analysis usually incurs long runtimes
and cumbersome custom code. A pipeline toolkit tailored to
Bayesian statisticians, the \code{jagstargets} R package leverages
\code{targets} and \code{R2jags} to ease this burden.
\code{jagstargets} makes it super easy to set up scalable
JAGS pipelines that automatically parallelize the computation
and skip expensive steps when the results are already up to date.
Minimal custom code is required, and there is no need to manually
configure branching, so usage is much easier than \code{targets} alone.
}
\seealso{
\url{https://docs.ropensci.org/jagstargets/}, \code{\link[=tar_jags]{tar_jags()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_example_data.R
\name{tar_jags_example_data}
\alias{tar_jags_example_data}
\title{Simulate example JAGS data.}
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
\code{\link[=tar_jags_rep_summary]{tar_jags_rep_summary()}}. Contains the
regression coefficient \code{beta} (numeric of length 1)
and prior predictive data \code{y} (numeric vector).
}
}
\usage{
tar_jags_example_data(n = 10L)
}
\arguments{
\item{n}{Integer of length 1, number of data points.}
}
\value{
List, dataset compatible with the model file from
\code{\link[=tar_jags_example_file]{tar_jags_example_file()}}. The output has a \code{.join_data}
element so the true value of \code{beta} from the simulation
is automatically appended to the \code{beta} rows of the
summary output.
}
\description{
An example dataset compatible with the model file
from \code{\link[=tar_jags_example_file]{tar_jags_example_file()}}. The output has a \code{.join_data}
element so the true value of \code{beta} from the simulation
is automatically appended to the \code{beta} rows of the
summary output.
}
\details{
The \code{tar_jags_example_data()} function draws a JAGS
dataset from the prior predictive distribution of the
model from \code{\link[=tar_jags_example_file]{tar_jags_example_file()}}. First, the
regression coefficient \code{beta} is drawn from its standard
normal prior, and the covariate \code{x} is computed.
Then, conditional on the \code{beta} draws and the covariate,
the response vector \code{y} is drawn from its
Normal(\code{x * beta}, 1) likelihood.
}
\examples{
tar_jags_example_data()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_rep_summary.R
\name{tar_jags_rep_summary}
\alias{tar_jags_rep_summary}
\title{Tidy posterior summaries from multiple MCMCs per model}
\usage{
tar_jags_rep_summary(
  name,
  jags_files,
  parameters.to.save,
  data = list(),
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  batches = 1L,
  reps = 1L,
  combine = TRUE,
  n.cluster = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = as.integer(n.iter/2),
  n.thin = 1,
  jags.module = c("glm", "dic"),
  inits = NULL,
  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister"),
  jags.seed = 123,
  stdout = NULL,
  stderr = NULL,
  progress.bar = "text",
  refresh = 0,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
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

\item{jags_files}{Character vector of JAGS model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{tools::file_path_sans_ext(basename(jags_files))} will be used
as target name suffixes. If \code{jags_files} is a named vector,
the suffixed will come from \code{names(jags_files)}.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{Code to generate the \code{data} list for the JAGS model.
Optionally include a \code{.join_data} element to join parts of the data
to correspondingly named parameters in the summary output.
See the vignettes for details.}

\item{variables}{Character vector of model parameter names.
The output posterior summaries are restricted to these variables.}

\item{summaries}{List of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{List of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{batches}{Number of batches. Each batch runs a model \code{reps} times.}

\item{reps}{Number of replications per batch. Ideally, each rep
should produce its own random dataset using the code
supplied to \code{data}.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frames
of posterior summaries and other data frames returned by targets.
We recommend efficient data frame formats
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
\code{tar_jags_rep_summary()} returns list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{jags_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_jags_rep_dic(name = x, jags_files = "y.jags")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the JAGS model file. Returns
a character vector of length 1 with the path to the JAGS
model file.
\item \code{x_lines_y}: read the contents of the JAGS model file
for safe transport to parallel workers.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple JAGS
datasets from the R expression in the \code{data} argument.
Each dynamic branch returns a batch of JAGS data lists.
\item \code{x_y}: run JAGS on each dataset from \code{x_data}.
Each dynamic branch returns a tidy data frame of summaries
for each batch of data.
\item \code{x}: combine all the batches from \code{x_y} into a non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame with all summaries
from all the branches of \code{x_y}.
}
}
\description{
Run multiple MCMCs on simulated datasets
and return posterior summaries and the effective number of parameters
for each run.
}
\details{
The MCMC targets use \code{R2jags::jags()} if \code{n.cluster} is \code{1} and
\code{R2jags::jags.parallel()} otherwise. Most arguments to \code{tar_jags()}
are forwarded to these functions.
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
if (requireNamespace("R2jags", quietly = TRUE)) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(jagstargets)
# Do not use a temp file for a real project
# or else your targets will always rerun.
tmp <- tempfile(pattern = "", fileext = ".jags")
tar_jags_example_file(tmp)
list(
  tar_jags_rep_summary(
    your_model,
    jags_files = tmp,
    data = tar_jags_example_data(),
    parameters.to.save = "beta",
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_rep_draws.R
\name{tar_jags_rep_draws}
\alias{tar_jags_rep_draws}
\title{Tidy posterior draws from multiple MCMCs per model}
\usage{
tar_jags_rep_draws(
  name,
  jags_files,
  parameters.to.save,
  data = list(),
  batches = 1L,
  reps = 1L,
  combine = FALSE,
  n.cluster = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = as.integer(n.iter/2),
  n.thin = 1,
  jags.module = c("glm", "dic"),
  inits = NULL,
  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister"),
  jags.seed = 123,
  stdout = NULL,
  stderr = NULL,
  progress.bar = "text",
  refresh = 0,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = "qs",
  format_df = "fst_tbl",
  error = targets::tar_option_get("error"),
  memory = "transient",
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

\item{jags_files}{Character vector of JAGS model files. If you
supply multiple files, each model will run on the one shared dataset
generated by the code in \code{data}. If you supply an unnamed vector,
\code{tools::file_path_sans_ext(basename(jags_files))} will be used
as target name suffixes. If \code{jags_files} is a named vector,
the suffixed will come from \code{names(jags_files)}.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{Code to generate the \code{data} list for the JAGS model.
Optionally include a \code{.join_data} element to join parts of the data
to correspondingly named parameters in the summary output.
See the vignettes for details.}

\item{batches}{Number of batches. Each batch runs a model \code{reps} times.}

\item{reps}{Number of replications per batch. Ideally, each rep
should produce its own random dataset using the code
supplied to \code{data}.}

\item{combine}{Logical, whether to create a target to
combine all the model results
into a single data frame downstream. Convenient, but
duplicates data.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the data frames
of posterior summaries and other data frames returned by targets.
We recommend efficient data frame formats
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
\code{tar_jags_rep_draws()} returns list of target objects.
See the "Target objects" section for
background.
The target names use the \code{name} argument as a prefix, and the individual
elements of \code{jags_files} appear in the suffixes where applicable.
As an example, the specific target objects returned by
\code{tar_jags_rep_dic(name = x, jags_files = "y.jags")}
are as follows.
\itemize{
\item \code{x_file_y}: reproducibly track the JAGS model file. Returns
a character vector of length 1 with the path to the JAGS
model file.
\item \code{x_lines_y}: read the contents of the JAGS model file
for safe transport to parallel workers.
Returns a character vector of lines in the model file.
\item \code{x_data}: use dynamic branching to generate multiple JAGS
datasets from the R expression in the \code{data} argument.
Each dynamic branch returns a batch of JAGS data lists.
\item \code{x_y}: run JAGS on each dataset from \code{x_data}.
Each dynamic branch returns a tidy data frame of draws
for each batch of data.
\item \code{x}: combine all the batches from \code{x_y} into a non-dynamic target.
Suppressed if \code{combine} is \code{FALSE}.
Returns a long tidy data frame with all draws
from all the branches of \code{x_y}.
}
}
\description{
Run multiple MCMCs on simulated datasets
and return posterior samples and the effective number of parameters
for each run.
}
\details{
The MCMC targets use \code{R2jags::jags()} if \code{n.cluster} is \code{1} and
\code{R2jags::jags.parallel()} otherwise. Most arguments to \code{tar_jags()}
are forwarded to these functions.
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
if (requireNamespace("R2jags", quietly = TRUE)) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
library(jagstargets)
# Do not use a temp file for a real project
# or else your targets will always rerun.
tmp <- tempfile(pattern = "", fileext = ".jags")
tar_jags_example_file(tmp)
list(
  tar_jags_rep_draws(
    your_model,
    jags_files = tmp,
    data = tar_jags_example_data(),
    parameters.to.save = "beta",
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags.R
\name{tar_jags_run}
\alias{tar_jags_run}
\title{Run a JAGS model and return the whole output object.}
\usage{
tar_jags_run(
  jags_lines,
  parameters.to.save,
  data,
  inits,
  n.cluster,
  n.chains,
  n.iter,
  n.burnin,
  n.thin,
  jags.module,
  RNGname,
  jags.seed,
  stdout,
  stderr,
  progress.bar,
  refresh
)
}
\arguments{
\item{jags_lines}{Character vector of lines from a JAGS model file.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}
}
\value{
An \code{R2jags} output object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_rep.R
\name{tar_jags_rep_run}
\alias{tar_jags_rep_run}
\title{Run a batch of iterations for a jags model
and return only strategic output.}
\usage{
tar_jags_rep_run(
  jags_lines,
  jags_name,
  jags_file,
  parameters.to.save,
  data,
  variables = NULL,
  summaries = NULL,
  summary_args = NULL,
  reps,
  output,
  n.cluster = n.cluster,
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  jags.module = jags.module,
  inits = inits,
  RNGname = RNGname,
  jags.seed = jags.seed,
  stdout = stdout,
  stderr = stderr,
  progress.bar = progress.bar,
  refresh = refresh
)
}
\arguments{
\item{jags_lines}{Character vector of lines from a JAGS file.}

\item{jags_name}{Friendly suffix of the jags model target.}

\item{jags_file}{Original path to the input jags file.}

\item{parameters.to.save}{Model parameters to save, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{data}{A list, the original JAGS dataset.}

\item{variables}{Character vector of model parameter names.
The output posterior summaries are restricted to these variables.}

\item{summaries}{List of summary functions passed to \code{...} in
\code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{summary_args}{List of summary function arguments passed to
\code{.args} in \code{posterior::summarize_draws()} through \verb{$summary()}
on the \code{CmdStanFit} object.}

\item{output}{Character of length 1 denoting the type of output \code{tibble}
to return: \code{"draws"} for MCMC samples (which could take up a lot of space)
\code{"summary"} for lightweight posterior summary statistics,
and \code{"dic"} for the overall deviance information criterion
and effective number of parameters}

\item{n.cluster}{Number of parallel processes, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.chains}{Number of MCMC chains, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.iter}{Number if iterations (including warmup), passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.burnin}{Number of warmup iterations, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{n.thin}{Thinning interval, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.module}{Character vector of JAGS modules to load, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{inits}{Initial values of model parameters, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{RNGname}{Choice of random number generator, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{jags.seed}{Seeds to apply to JAGS, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{stdout}{Character of length 1, file path to write the stdout stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stdout.
Does not apply to messages, warnings, or errors.}

\item{stderr}{Character of length 1, file path to write the stderr stream
of the model when it runs. Set to \code{NULL} to print to the console.
Set to \code{R.utils::nullfile()} to suppress stderr.
Does not apply to messages, warnings, or errors.}

\item{progress.bar}{Type of progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}

\item{refresh}{Frequency for refreshing the progress bar, passed to
\code{R2jags::jags()} or \code{R2jags::jags.parallel()}.
See the argument documentation of the
\code{R2jags::jags()} and \code{R2jags::jags.parallel()} help files for details.}
}
\value{
A data frame of posterior summaries.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_example_file.R
\name{tar_jags_example_file}
\alias{tar_jags_example_file}
\title{Write an example JAGS model file.}
\usage{
tar_jags_example_file(path = tempfile(pattern = "", fileext = ".jags"))
}
\arguments{
\item{path}{Character of length 1, file path to write the model file.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Overwrites the file at \code{path} with a built-in example
JAGS model file.
}
\examples{
path <- tempfile(pattern = "", fileext = ".jags")
tar_jags_example_file(path = path)
writeLines(readLines(path))
}
