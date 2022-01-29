
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
