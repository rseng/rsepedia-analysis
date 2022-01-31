
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medoutcon`

<!-- badges: start -->

[![R-CMD-check](https://github.com/nhejazi/medoutcon/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/medoutcon/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/medoutcon/master.svg)](https://codecov.io/github/nhejazi/medoutcon?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5809519.svg)](https://doi.org/10.5281/zenodo.5809519)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03979/status.svg)](https://doi.org/10.21105/joss.03979)
<!-- badges: end -->

> Efficient Causal Mediation Analysis for the Natural and Interventional
> Effects

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](https://kararudolph.github.io/)

-----

## What’s `medoutcon`?

The `medoutcon` R package provides facilities for efficient estimation
of path-specific (in)direct effects that measure the impact of a
treatment variable \(A\) on an outcome variable \(Y\), through a direct
path (through \(A\) only) and an indirect path (through a set of
mediators \(M\) only). In the presence of an intermediate
<b>med</b>iator-<b>out</b>come <b>con</b>founder \(Z\), itself affected
by the treatment \(A\), these correspond to the *interventional*
(in)direct effects described by Dı́az et al. (2020), though similar (yet
less general) effect definitions and/or estimation strategies have
appeared in VanderWeele, Vansteelandt, and Robins (2014), Rudolph et al.
(2017), Zheng and van der Laan (2017), and Benkeser and Ran (2021). When
no intermediate confounders are present, these effect definitions
simplify to the well-studied *natural* (in)direct effects, and our
estimators are analogs of those formulated by Zheng and van der Laan
(2012). Both an efficient one-step bias-corrected estimator with
cross-fitting (Pfanzagl and Wefelmeyer 1985; Zheng and van der Laan
2011; Chernozhukov et al. 2018) and a cross-validated targeted minimum
loss estimator (TMLE) (van der Laan and Rose 2011; Zheng and van der
Laan 2011) are made available. `medoutcon` integrates with the [`sl3` R
package](https://github.com/tlverse/sl3) (Coyle et al. 2021) to leverage
statistical machine learning in the estimation procedure.

-----

## Installation

Install the most recent *stable release* from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/medoutcon")
```

-----

## Example

To illustrate how `medoutcon` may be used to estimate stochastic
interventional (in)direct effects of the exposure (`A`) on the outcome
(`Y`) in the presence of mediator(s) (`M`) and a mediator-outcome
confounder (`Z`), consider the following example:

``` r
library(data.table)
library(tidyverse)
#> ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
#> ✔ ggplot2 3.3.5     ✔ purrr   0.3.4
#> ✔ tibble  3.1.6     ✔ dplyr   1.0.7
#> ✔ tidyr   1.1.4     ✔ stringr 1.4.0
#> ✔ readr   2.1.1     ✔ forcats 0.5.1
#> ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::between()   masks data.table::between()
#> ✖ dplyr::filter()    masks stats::filter()
#> ✖ dplyr::first()     masks data.table::first()
#> ✖ dplyr::lag()       masks stats::lag()
#> ✖ dplyr::last()      masks data.table::last()
#> ✖ purrr::transpose() masks data.table::transpose()
library(medoutcon)
#> medoutcon v0.1.6: Efficient Natural and Interventional Causal Mediation Analysis
set.seed(1584)

# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  ## baseline covariates
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure
  a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))

  ## mediator -- could be multivariate
  m <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w[, -3] + a - z)))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + m)))

  ## construct output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

# set seed and simulate example data
example_data <- make_example_data()
w_names <- str_subset(colnames(example_data), "W")
m_names <- str_subset(colnames(example_data), "M")

# quick look at the data
head(example_data)
#>    W_1 W_2 W_3 A Z M Y
#> 1:   1   0   1 0 0 0 1
#> 2:   0   1   0 0 0 1 0
#> 3:   1   1   1 1 0 1 1
#> 4:   0   1   1 0 0 1 0
#> 5:   0   0   0 0 0 1 1
#> 6:   1   0   1 1 0 1 0

# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   effect = "direct",
                   estimator = "onestep")
os_de
#> Interventional Direct Effect
#> Estimator: onestep
#> Estimate: -0.075
#> Std. Error: 0.056
#> 95% CI: [-0.186, 0.035]

# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     effect = "direct",
                     estimator = "tmle")
tmle_de
#> Interventional Direct Effect
#> Estimator: tmle
#> Estimate: -0.078
#> Std. Error: 0.059
#> 95% CI: [-0.193, 0.037]
```

For details on how to use data adaptive regression (machine learning)
techniques in the estimation of nuisance parameters, consider consulting
the vignette that accompanies the package.

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/medoutcon/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/medoutcon/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `medoutcon` R package, please cite the following:

``` 
    @article{diaz2020nonparametric,
      title={Non-parametric efficient causal mediation with intermediate
        confounders},
      author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S and Rudolph, Kara E
        and {van der Laan}, Mark J},
      year={2020},
      url = {https://arxiv.org/abs/1912.09936},
      doi = {10.1093/biomet/asaa085},
      journal={Biometrika},
      volume = {108},
      number = {3},
      pages = {627--641},
      publisher={Oxford University Press}
    }

    @article{hejazi2022medoutcon-joss,
      author = {Hejazi, Nima S and Rudolph, Kara E and D{\'\i}az,
        Iv{\'a}n},
      title = {{medoutcon}: Nonparametric efficient causal mediation
        analysis with machine learning in {R}},
      year = {2022},
      doi = {10.21105/joss.03979},
      url = {https://doi.org/10.21105/joss.03979},
      journal = {Journal of Open Source Software},
      publisher = {The Open Journal}
    }

    @software{hejazi2022medoutcon-rpkg,
      author={Hejazi, Nima S and D{\'\i}az, Iv{\'a}n and Rudolph, Kara E},
      title = {{medoutcon}: Efficient natural and interventional causal
        mediation analysis},
      year  = {2022},
      doi = {10.5281/zenodo.5809519},
      url = {https://github.com/nhejazi/medoutcon},
      note = {R package version 0.1.6}
    }
```

-----

## License

© 2020-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2020-2022 Nima S. Hejazi
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

<div id="refs" class="references">

<div id="ref-benkeser2020nonparametric">

Benkeser, David, and Jialu Ran. 2021. “Nonparametric Inference for
Interventional Effects with Multiple Mediators.” *Journal of Causal
Inference*. <https://doi.org/10.1515/jci-2020-0018>.

</div>

<div id="ref-chernozhukov2018double">

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

</div>

<div id="ref-coyle-gh-sl3">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2021. *`sl3`: Modern Machine Learning Pipelines for Super
Learning* (version 1.4.4). <https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-diaz2020nonparametric">

Dı́az, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2020. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” *Biometrika* 108 (3): 627–41.
<https://doi.org/10.1093/biomet/asaa085>.

</div>

<div id="ref-pfanzagl1985contributions">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88. <https://doi.org/10.1007/978-1-4612-5769-1>.

</div>

<div id="ref-rudolph2017robust">

Rudolph, Kara E, Oleg Sofrygin, Wenjing Zheng, and Mark J van der Laan.
2017. “Robust and Flexible Estimation of Stochastic Mediation Effects: A
Proposed Method and Example in a Randomized Trial Setting.”
*Epidemiologic Methods* 7 (1). <https://doi.org/10.1515/em-2017-0007>.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vanderweele2014effect">

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300.
<https://doi.org/10.1097/ede.0000000000000034>.

</div>

<div id="ref-zheng2011cross">

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning: Causal Inference
for Observational and Experimental Data*, 459–74. Springer.
<https://doi.org/10.1007/978-1-4419-9782-1_27>.

</div>

<div id="ref-zheng2012targeted">

———. 2012. “Targeted Maximum Likelihood Estimation of Natural Direct
Effects.” *International Journal of Biostatistics* 8 (1).
<https://doi.org/10.2202/1557-4679.1361>.

</div>

<div id="ref-zheng2017longitudinal">

———. 2017. “Longitudinal Mediation Analysis with Time-Varying Mediators
and Exposures, with Application to Survival Outcomes.” *Journal of
Causal Inference* 5 (2). <https://doi.org/10.1515/jci-2016-0006>.

</div>

</div>
# medoutcon 0.1.5

* For user clarity, the name of the argument for providing externally computed
  observation-level weights has changed (from `ext_weights`) to `svy_weights`.
* Support for the natural direct and indirect effects has been added, requiring
  the addition of the new internal argument `effect_type` across functions for
  estimation, including `cv_eif()`, `est_onestep()`, and `est_tml()`. When
  `Z = NULL` is set in `medoutcon()`, a natural effect estimate corresponding to
  the argument `effect` is returned instead of an interventional effect.
* The `summary()` and `print()` methods have been updated to allow handling of
  natural effects and counterfactual means under arbitrary contrasts.

# medoutcon 0.1.0

* An initial public release of this package, version 0.1.0, which includes
  support for external observation-level weights.
# Contributing to `medoutcon` development

We, the authors of the `medoutcon` R package, use the same guide as is used for
contributing to the development of the popular `tidyverse` ecosystem of R
packages. This document is simply a formal re-statement of that fact.

The goal of this guide is to help you get up and contributing to `medoutcon` as
quickly as possible. The guide is divided into two main pieces:

* Filing a bug report or feature request in an issue.
* Suggesting a change via a pull request.

## Issues

When filing an issue, the most important thing is to include a minimal
reproducible example so that we can quickly verify the problem, and then figure
out how to fix it. There are three things you need to include to make your
example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.

2.  The easiest way to include **data** is to use `dput()` to generate the R
    code to recreate it.

3.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:

    * make sure you've used spaces and your variable names are concise, but
      informative

    * use comments to indicate where your problem lies

    * do your best to remove everything that is not related to the problem.
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a
fresh R session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the output
of `sessionInfo()`.)

## Pull requests

To contribute a change to `medoutcon`, you follow these steps:

1. Create a branch in git and make your changes.
2. Push branch to GitHub and issue pull request (PR).
3. Discuss the pull request.
4. Iterate until either we accept the PR or decide that it's not a good fit for
   `medoutcon`.

Each of these steps are described in more detail below. This might feel
overwhelming the first time you get set up, but it gets easier with practice.

If you're not familiar with git or GitHub, please start by reading
<http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivates the
   need for change. Please describe the problem your PR addresses and show
   how your pull request solves it as concisely as possible.

   Also include this motivation in `NEWS` so that when a new release of
   `medoutcon` comes out it's easy for users to see what's changed. Add your
   item at the top of the file and use markdown for formatting. The
   news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.

3.  __Use `medoutcon` coding style__. To do so, please follow the [official
    `tidyverse` style guide](http://style.tidyverse.org). Maintaining a
    consistent style across the whole code base makes it much easier to jump
    into the code. If you're modifying existing `medoutcon` code that doesn't
    follow the style guide, a separate pull request to fix the style would be
    greatly appreciated.

4.  If you're adding new parameters or a new function, you'll also need
    to document them with [`roxygen2`](https://github.com/klutometis/roxygen).
    Make sure to re-run `devtools::document()` on the code before submitting.

This seems like a lot of work but don't worry if your pull request isn't
perfect. It's a learning process. A pull request is a process, and unless
you've submitted a few in the past it's unlikely that your pull request will be
accepted as is. Please don't submit pull requests that change existing
behaviour. Instead, think about how you can add a new feature in a minimally
invasive way.
---
title: "`medoutcon`: Nonparametric efficient causal mediation analysis with machine learning in `R`"
tags:
  - causal inference
  - machine learning
  - semiparametric estimation
  - mediation analysis
  - natural direct effect
  - interventional direct effect
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1
  - name: Kara E. Rudolph
    orcid: 0000-0002-9417-7960
    affiliation: 2
  - name: Iván Díaz
    orcid: 0000-0001-9056-2047
    affiliation: 1
affiliations:
  - name: Division of Biostatistics, Department of Population Health Sciences, Weill Cornell Medicine, USA
    index: 1
  - name: Department of Epidemiology, Mailman School of Public Health, Columbia University, USA
    index: 2
date: 30 December 2021
bibliography: ../inst/REFERENCES.bib
---

# Summary

Science is most often concerned with questions of _mechanism_. In myriad
applications, only the portion of the causal effect of an exposure on an outcome
through a particular pathway under study is of interest. The study of such
path-specific, or mediation, effects has a rich history, first undertaken
scientifically by @wright1921correlation and @wright1934method. Today, the study
of such effects has attracted a great deal of attention in statistics and causal
inference, inspired by applications in disciplines ranging from epidemiology and
vaccinology to psychology and economics. Examples include understanding the
biological mechanisms by which vaccines causally alter infection risk
[@hejazi2020efficient; @benkeser2021inference], assessing the effect of novel
pharmacological therapies on substance abuse disorder relapse
[@hejazi2021nonparametric; @rudolph2020explaining], and evaluating the effects
of housing vouchers on adolescent development [@rudolph2021helped]. The
`medoutcon` `R` package provides researchers in each of these disciplines, and
in others, with the tools necessary to implement statistically efficient
estimators of the _interventional_ direct and indirect effects
[@diaz2020nonparametric] (for brevity, henceforth, (in)direct effects),
a recently formulated set of causal effects robust to the presence of
confounding of the mediator-outcome relationship by the exposure. In cases where
such confounding is a nonissue, the interventional (in)direct effects
[@vanderweele2014effect] reduce to the well-studied _natural_ (in)direct effects
[@robins1992identifiability; @pearl2001direct], for which `medoutcon` provides
efficient estimators similar to those of @zheng2012targeted. By readily
incorporating the use of machine learning in the estimation of nuisance
parameters (through integration with the `sl3` `R` package [@coyle-gh-sl3] of
the `tlverse` ecosystem [@vdl2022targeted]), `medoutcon` incorporates
state-of-the-art non/semi-parametric estimation techniques, facilitating their
adoption in a vast array of settings.

# Statement of Need

While there is demonstrable interest in causal mediation analysis in a large
variety of disciplines, thoughtfully implementing data analysis strategies based
on recent developments in this area is challenging. Contributions in the causal
inference and statistics literature largely fall into two key areas. Broadly,
the study of identification outlines novel causal effect parameters with
properties desirable in real-world settings (e.g., the interventional effects,
which can be learned under mediator-outcome confounding) and untestable
assumptions under which a statistical functional corresponds to a causal
estimand of interest. A complementary line of study develops non/semi-parametric
efficiency theory for the statistical functionals outlined in the causal
identification literature, allowing for their robust estimation with modern
techniques from machine learning. Neither concerns itself with opening the door
to applying these estimators in real-world data analyses. Moreover, the
implementation of open source software for efficient estimators of causal
effects is complex -- for such a task, the data scientist must be
knowledgeable of causal inference, semiparametric statistical theory, machine
learning, and the intersection of these disciplines, and that is to forego
mention of research software engineering best practices, including, for example,
unit/regression testing and automated continuous integration. The `medoutcon`
`R` package is a free, open source implementation of non/semi-parametric
efficient estimators of the natural and interventional (in)direct effects,
providing data scientists in research and in industry with access to
state-of-the-art statistical methodology for causal mediation analysis. Its
estimators have been interrogated in simulation studies and applied in
real-world data analyses. To the best of our knowledge, no other `R` package
provides similarly convenient access to multiply robust, non/semi-parametric
efficient estimators of causal mediation effects with a flexible interface to
accommodate machine learning of nuisance parameters.

# Natural and Interventional Causal Mediation Effects

To evaluate the causal effects of an exposure on an outcome through mediating
pathways, let's consider a dataset of $n$ units, where the observed data on
a single unit is assumed to have been generated by a nonparametric structural
equation model (NPSEM) [@pearl2009causality]:
\begin{align*}
  W &= f_W(U_W); A = f_A(W, U_A); Z=f_Z(W, A, U_Z);\\
  M &= f_M(W, A, Z, U_M); Y = f_Y(W, A, Z, M, U_Y),
\end{align*}
where $W$ are baseline (pre-exposure) covariates, $A \in \{0,1\}$ is the
(binary) exposure of interest, $Z$ is an intermediate confounder of the
mediator-outcome relationship and is affected by exposure $A$, $M$ represents
mediating variables, and $Y$ is the outcome. This NPSEM admits an equivalent
representation as a directed acyclic graph (or DAG), in which each variable is
a node and dependencies are represented by directed paths between the nodes. The
natural (in)direct effects cannot generally be identified (i.e., learned from
the observed data) in the presence of intermediate confounding, so, for now, we
make the simplifying assumption that the intermediate variable $Z$ is absent. In
this simple case, the population average treatment effect (ATE) -- that is, the
total effect of $A$ on $Y$, comparing two exposure contrasts $\{a', a^{\star}\}$
-- may be decomposed into the natural direct effect (NDE) and the natural
indirect effect (NIE) as
\begin{equation*}
  \mathbb{E}[Y(a') - Y(a^{\star})] =
    \underbrace{\mathbb{E}[Y(a', M(a')) - Y(a',
      M(a^{\star}))]}_{\text{Indirect effect (through $M$)}} +
    \underbrace{\mathbb{E}[Y(a', M(a^{\star})) - Y(a^{\star},
      M(a^{\star}))]}_{\text{Direct effect (not through $M$)}},
\end{equation*}
where the _counterfactual_ variables $Y(\cdot)$ are _potential outcomes_
[@imbens2015causal; @hernan2021causal] -- that is, $Y(a')$ is the value that the
outcome would take when the exposure is set to level $a'$, possibly contrary to
fact. Similarly, $M(a^{\star})$ is the value that the mediators would take when
the exposure is set to level $a^{\star}$, as the result of an intervention, for
example. The NIE captures the effect of the exposure $A$ on $Y$ through the
mediating variables $M$ while the NDE captures the effect of $A$ on $Y$ through
all other pathways. @robins1992identifiability and @pearl2001direct
independently studied this decomposition within the potential outcomes and NPSEM
frameworks, respectively. In both cases, the NDE and NIE are derived from the
ATE by introducing a decomposition term that deterministically sets the values
of the exposure and mediators to differing values by the application of _static_
interventions. As regards estimation, @tchetgen2012semiparametric and
@zheng2012targeted outlined non/semi-parametric efficiency theory for developing
estimators of the NDE and NIE and proposed efficient estimators of these causal
quantities.

The presence of intermediate confounders $Z$ often cannot be ruled out in
real-world data analysis scenarios. Such post-exposure variables, which are
affected by $A$ and affect both $M$ and $Y$, complicate efforts to disentangle
the effect of $A$ on $Y$ through paths involving $M$ and other paths.
Recognizing the limitations of the natural effects in these settings,
@didelez2006direct, @petersen2006estimation, @vanderweele2014effect, and
@rudolph2017robust, among others, contributed to the development of the
interventional (in)direct effects. Unlike the decomposition strategy that
delineates the NDE and NIE, these effects require a more sophisticated approach
to identification, relying upon _stochastic_ interventions on the mediator(s),
which require random draws from the mediator's post-intervention distribution
rather than the setting of fixed counterfactual values. Specifically, for the
two exposure contrasts $\{a', a^{\star}\}$, the effect of $A$ on $Y$ can be
defined as the difference in expected outcome in the hypothetical worlds in
which $(A,M) = (a', G_{a'})$ versus $(A,M) = (a^{\star}, G_{a^{\star}})$. Here,
$G_a$ denotes a random draw from the conditional distribution of $M_a$
conditional on $W$, as defined by a stochastic intervention. The direct and
indirect effects are defined as follows
\begin{equation*}
\mathbb{E}[Y(a', G_{a'}) - Y(a^{\star}, G_{a^{\star}})] =
  \underbrace{\mathbb{E}[Y(a', G_{a'}) - Y(a',
    G_{a^{\star}})]}_{\text{Indirect effect (through $M$)}} +
  \underbrace{\mathbb{E}[Y(a', G_{a^{\star}}) - Y(a^{\star},
      G_{a^{\star}})]}_{\text{Direct effect (not through $M$)}}.
\end{equation*}
Like the NDE, this interventional direct effect measures the effects through all
paths avoiding the mediating variables. Analogous to the NIE, the interventional
indirect effect measures the effect through paths involving the mediators. Note,
however, that natural and interventional mediation effects have different
interpretations. That is, the interventional indirect effect measures the effect
of fixing the exposure at $a'$ while setting the mediator to a random draw
$G_{a^{\star}}$ (i.e., under an intervention setting the exposure to
$a^{\star}$) versus a random draw $G_{a'}$ (i.e., after setting the exposure to
$a'$), given covariates $W$. Intuitively, the interventional effects remain
identifiable under intermediate confounding since the stochastic intervention on
the mediators breaks the relationship between $Z$ and $M$. Prior to the work of
@diaz2020nonparametric, and contemporaneous developments by
@benkeser2020nonparametric, non/semi-parametric efficiency theory for the
interventional (in)direct effects was unavailable. Recently, a novel family of
interventional effects, accommodating flexible stochastic interventions on the
exposure [@hejazi2021nonparametric], have been formulated as well.

# `medoutcon`'s Scope

Development of the `medoutcon` package began as a software accompaniment to the
theoretical developments of @diaz2020nonparametric. Where the investigations
of these authors outlined efficient estimators of the interventional (in)direct
effects, `medoutcon` implements these efficient estimators. Implemented in the
`R` language and environment for statistical computing [@R], `medoutcon` aims to
provide a simple application user interface (API) for convenience in a variety
of data analytic applications. Specifically, `medoutcon` -- via a single,
user-facing eponymous function `medoutcon()` -- provides access to both one-step
and targeted minimum loss (TML) estimators of these causal (in)direct effects.
State-of-the-art machine learning algorithms, including ensemble modeling
[@vdl2007super], may readily be used for the estimation of relevant nuisance
parameters, through a design that tightly couples `medoutcon` with the `sl3` `R`
package [@coyle-gh-sl3]. Cross-fitting is automatically incorporated, via the
`origami` `R` package [@coyle2018origami; @coyle-cran-origami], in computing the
efficient estimators, allowing for some common but restrictive theoretical
regularity conditions to be relaxed [@bickel1993efficient; @zheng2011cross;
@chernozhukov2017double].

Beyond implementing the interventional (in)direct effects, `medoutcon`
additionally allows for the natural (in)direct effects to be estimated when
intermediate confounders are omitted from the call to the `medoutcon()` function
(i.e., by setting `Z = NULL`). This feature is based on a correspondence between
the identifying statistical functionals of the natural and interventional
(in)direct effects in the absence of intermediate confounding. In this
simplified case, the efficient estimators of the interventional (in)direct
effects formulated by @diaz2020nonparametric are analogous to the efficient
estimators of the natural (in)direct effects formulated by @zheng2012targeted.
By supporting this case, `medoutcon` serves as a one-stop tool for estimating
these scientifically informative causal mediation effects, allowing for
practicing data scientists and applied statisticians to deploy cutting-edge
estimators of the natural and interventional (in)direct effects through
a unified API.

# Availability

The `medoutcon` package is publicly available [via
GitHub](https://github.com/nhejazi/medoutcon), with plans for submission to the
Comprehensive `R` Archive Network, pending the inclusion of its dependencies
(`sl3`, in particular) in that repository. Use of the `medoutcon` package has
been extensively documented in the package's `README`, a vignette, and its
[documentation website](https://code.nimahejazi.org/medoutcon). Ongoing
development of the package incorporates research and data science software
engineering best practices, including a suite of unit tests and automated
continuous integration checking. `medoutcon` has and will continue to be used in
the teaching of conference workshops on modern causal mediation analysis (e.g.,
see [recent materials from SER
2021](https://code.nimahejazi.org/ser2021_mediation_workshop/)).

# Acknowledgments

NSH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

# References

---
output:
  rmarkdown::github_document
bibliography: "inst/REFERENCES.bib"
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# R/`medoutcon`

<!-- badges: start -->
[![R-CMD-check](https://github.com/nhejazi/medoutcon/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/medoutcon/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/nhejazi/medoutcon/master.svg)](https://codecov.io/github/nhejazi/medoutcon?branch=master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5809519.svg)](https://doi.org/10.5281/zenodo.5809519)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03979/status.svg)](https://doi.org/10.21105/joss.03979)
<!-- badges: end -->

> Efficient Causal Mediation Analysis for the Natural and Interventional Effects

__Authors:__ [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](https://kararudolph.github.io/)

---

## What's `medoutcon`?

The `medoutcon` R package provides facilities for efficient estimation of
path-specific (in)direct effects that measure the impact of a treatment variable
$A$ on an outcome variable $Y$, through a direct path (through $A$ only) and an
indirect path (through a set of mediators $M$ only). In the presence of an
intermediate <b>med</b>iator-<b>out</b>come <b>con</b>founder $Z$, itself
affected by the treatment $A$, these correspond to the _interventional_
(in)direct effects described by @diaz2020nonparametric, though similar (yet less
general) effect definitions and/or estimation strategies have appeared in
@vanderweele2014effect, @rudolph2017robust, @zheng2017longitudinal, and
@benkeser2020nonparametric. When no intermediate confounders are present, these
effect definitions simplify to the well-studied _natural_ (in)direct effects,
and our estimators are analogs of those formulated by @zheng2012targeted.  Both
an efficient one-step bias-corrected estimator with cross-fitting
[@pfanzagl1985contributions; @zheng2011cross; @chernozhukov2018double] and a
cross-validated targeted minimum loss estimator (TMLE) [@vdl2011targeted;
@zheng2011cross] are made available. `medoutcon` integrates with the [`sl3` R
package](https://github.com/tlverse/sl3) [@coyle-gh-sl3] to leverage statistical
machine learning in the estimation procedure.

---

## Installation

Install the most recent _stable release_ from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

```{r gh-master-installation, eval=FALSE}
remotes::install_github("nhejazi/medoutcon")
```

---

## Example

To illustrate how `medoutcon` may be used to estimate stochastic interventional
(in)direct effects of the exposure (`A`) on the outcome (`Y`) in the presence of
mediator(s) (`M`) and a mediator-outcome confounder (`Z`), consider the
following example:

```{r example, warning=FALSE}
library(data.table)
library(tidyverse)
library(medoutcon)
set.seed(1584)

# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  ## baseline covariates
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure
  a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))

  ## mediator -- could be multivariate
  m <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w[, -3] + a - z)))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + m)))

  ## construct output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

# set seed and simulate example data
example_data <- make_example_data()
w_names <- str_subset(colnames(example_data), "W")
m_names <- str_subset(colnames(example_data), "M")

# quick look at the data
head(example_data)

# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   effect = "direct",
                   estimator = "onestep")
os_de

# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     effect = "direct",
                     estimator = "tmle")
tmle_de
```

For details on how to use data adaptive regression (machine learning) techniques
in the estimation of nuisance parameters, consider consulting the vignette that
accompanies the package.

---

## Issues

If you encounter any bugs or have any specific feature requests, please [file an
issue](https://github.com/nhejazi/medoutcon/issues).

---

## Contributions

Contributions are very welcome. Interested contributors should consult our
[contribution
guidelines](https://github.com/nhejazi/medoutcon/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

---

## Citation

After using the `medoutcon` R package, please cite the following:

        @article{diaz2020nonparametric,
          title={Non-parametric efficient causal mediation with intermediate
            confounders},
          author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S and Rudolph, Kara E
            and {van der Laan}, Mark J},
          year={2020},
          url = {https://arxiv.org/abs/1912.09936},
          doi = {10.1093/biomet/asaa085},
          journal={Biometrika},
          volume = {108},
          number = {3},
          pages = {627--641},
          publisher={Oxford University Press}
        }

        @article{hejazi2022medoutcon-joss,
          author = {Hejazi, Nima S and Rudolph, Kara E and D{\'\i}az,
            Iv{\'a}n},
          title = {{medoutcon}: Nonparametric efficient causal mediation
            analysis with machine learning in {R}},
          year = {2022},
          doi = {10.21105/joss.03979},
          url = {https://doi.org/10.21105/joss.03979},
          journal = {Journal of Open Source Software},
          publisher = {The Open Journal}
        }

        @software{hejazi2022medoutcon-rpkg,
          author={Hejazi, Nima S and D{\'\i}az, Iv{\'a}n and Rudolph, Kara E},
          title = {{medoutcon}: Efficient natural and interventional causal
            mediation analysis},
          year  = {2022},
          doi = {10.5281/zenodo.5809519},
          url = {https://github.com/nhejazi/medoutcon},
          note = {R package version 0.1.6}
        }

---

## License

&copy; 2020-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license. See below
for details:
```
MIT License

Copyright (c) 2020-2022 Nima S. Hejazi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## References

---
title: "Efficient causal mediation analysis with the natural and interventional effects"
author: "[Nima Hejazi](https://nimahejazi.org), [Iván
  Díaz](https://www.idiaz.xyz/), and [Kara
  Rudolph](https://kararudolph.github.io/)"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: ../inst/REFERENCES.bib
vignette: >
  %\VignetteIndexEntry{Efficient causal mediation analysis with the natural and interventional effects}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Background and Motivations

An exposure of interest often affects an outcome directly, or indirectly by the
mediation of some intermediate variables. Identifying and quantifying the
mechanisms underlying causal effects is an increasingly popular endeavor in
public health, medicine, and the social sciences, as knowledge of such
mechanisms can improve understanding of both _why and how_ treatments can be
effective. Such mechanistic knowledge may be arguably even more important in
cases where treatments result in unanticipated ineffective or even harmful
effects.

Traditional techniques for mediation analysis fare poorly in the face of
intermediate confounding. Classical parameters like the natural (in)direct
effects face a lack of identifiability in cases where mediator-outcome (i.e.,
intermediate) confounders affected by exposure complicate the relationship
between the exposure, mediators, and outcome. @diaz2020nonparametric provide a
theoretical and computational study of the properties of newly developed
interventional (in)direct effect estimands within the non-parametric statistical
model. Among their contributions, @diaz2020nonparametric

- derive the efficient influence function (EIF), an key object in
  semiparametric efficiency theory;
- use the EIF to develop two asymptotically optimal, non-parametric estimators,
  each of which is capable of leveraging machine learning for the estimation of
  nuisance parameters; and
- present theoretical conditions under which their proposed estimators are
  consistent, multiply robust, and efficient.

# Problem Setup and Notation

The problem addressed by the work of @diaz2020nonparametric may be represented
by the following nonparametric structural equation model (NPSEM):
\begin{align*}
  W &= f_W(U_W); A = f_A(W, U_A); Z=f_Z(W, A, U_Z);\\ \nonumber
  M &= f_M(W, A, Z, U_M); Y = f_Y(W, A, Z, M, U_Y).
\end{align*}
In the NPSEM, $W$ denotes a vector of observed pre-treatment covariates, $A$
denotes a categorical treatment variable, $Z$ denotes an intermediate confounder
affected by treatment, $M$ denotes a (possibly multivariate) mediator, and $Y$
denotes a continuous or binary outcome.  The vector of exogenous factors
$U=(U_W,U_A,U_Z,U_M,U_Y)$, and the functions $f$, are assumed deterministic but
unknown.  Importantly, the NPSEM encodes a time-ordering between these variables
and allows the evaluation of counterfactual quantities defined by intervening on
a set of nodes of the NPSEM.  The observed data unit can be represented by the
random variable $O = (W, A, Z, M, Y)$; we consider access to $O_1, \ldots, O_n$,
a sample of $n$ i.i.d. observations of $O$.

@diaz2020nonparametric additionally define the following parameterizations,
familiarity with which will be useful for using the [`medoutcon` `R`
package](https://github.com/nhejazi/medoutcon). In particular, these authors
define $g(a \mid w)$ as the probability mass function of $A = a$ conditional on
$W = w$ and use $h(a \mid m, w)$ to denote the probability mass function of $A
= a$ conditional on $(M, W) = (m, w)$. Further, @diaz2020nonparametric use
$b(a, z, m, w)$ to denote the outcome regression function $\mathbb{E}(Y \mid A
= a, Z = z, M = m, W = w)$, as well as $q(z \mid a,w)$ and $r(z \mid a, m, w)$
to denote the corresponding conditional densities of $Z$.

# Interventional (In)Direct Effects

@diaz2020nonparametric define the _total effect_ of $A$ on $Y$ in terms of a
contrast between two user-supplied values $a', a^{\star} \in \mathcal{A}$.
Examination of the NPSEM reveals that there are four paths involved in this
effect, namely $A \rightarrow Y$, $A \rightarrow M \rightarrow Y$, $A
\rightarrow Z \rightarrow Y$, and $A \rightarrow Z \rightarrow M \rightarrow Y$.
Mediation analysis has classically considered the _natural direct effect_ (NDE)
and the _natural indirect effect_ (NIE), which are defined as
$\mathbb{E}_c(Y_{a', M_{a^{\star}}} - Y_{a^{\star}, M_{a^{\star}}})$ and
$\mathbb{E}_c(Y_{a',M_{a'}} - Y_{a',M_{a^{\star}}})$, respectively. The natural
direct effect measures the effect through paths _not_ involving the mediator
($A \rightarrow Y$ and $A \rightarrow Z \rightarrow Y$), whereas the natural
indirect effect measures the effect through paths involving the mediator
($A \rightarrow M \rightarrow Y$ and $A \rightarrow Z \rightarrow M \rightarrow
Y$). As the sum of the natural direct and indirect effects equals the average
treatment effect $\mathbb{E}_c(Y_1-Y_0)$, this effect decomposition is
appealing. Unfortunately, the natural direct and indirect effects are
not generally identified in the presence of an intermediate confounder affected
by treatment.

To circumvent this issue, @diaz2020nonparametric define the direct and indirect
effects using stochastic interventions on the mediator, following a strategy
previously outlined by @vanderweele2014effect and @rudolph2017robust, among
others. Let $G_a$ denote a random draw from the conditional distribution of
$M_a$ conditional on $W$. Consider the effect of $A$ on $Y$ defined as the
difference in expected outcome in hypothetical worlds in which $(A,M) = (a',
G_{a'})$ versus $(A,M) = (a^{\star}, G_{a^{\star}})$ with probability one, which
may be decomposed into direct and indirect effects as follows
\begin{equation*}
\mathbb{E}_c(Y_{a', G_{a'}} - Y_{a^{\star}, G_{a^{\star}}}) =
  \underbrace{\mathbb{E}_c(Y_{a', G_{a'}} - Y_{a',
    G_{a^{\star}}})}_{\text{Indirect effect (through $M$)}} +
  \underbrace{\mathbb{E}_c(Y_{a', G_{a^{\star}}} - Y_{a^{\star},
      G_{a^{\star}}})}_{\text{Direct effect (not through $M$)}}.
\end{equation*}
Like the natural direct effect, this interventional direct effect measures the
effects through paths not involving the mediator. Likewise, the interventional
indirect effect measures the effect through paths involving the mediator. Note,
however, that natural and interventional mediation effects have different
interpretations. That is, the interventional indirect effect measures the effect
of fixing the exposure at $a'$ while setting the mediator to a random draw
$G_{a^{\star}}$ from those with exposure $a'$ versus a random draw $G_{a'}$ from
those with exposure $a^{\star}$, given covariates $W$. As is clear from the
effect decomposition, the term $\theta_c = \mathbb{E}_c(Y_{a', G_{a^{\star}}})$
is required for estimation of both the interventional direct and indirect
effects; thus, @diaz2020nonparametric focus on estimation of this quantity.
Importantly, it has been shown that $\theta_c$ is identified by the statistical
functional
\begin{equation*}
  \theta = \int b(a', z, m, w) q(z \mid a', w) p(m \mid a^{\star}, w)
    p(w) d\nu(w,z,m)
\end{equation*}
under a set of standard identifiability conditions [@vanderweele2014effect],
which are further reviewed in @diaz2020nonparametric.

# Efficient Estimation

@diaz2020nonparametric define two efficient estimators of their interventional
(in)direct effects. These are based on the one-step estimation and targeted
minimum loss (TML) estimation frameworks, respectively. Briefly, both estimation
strategies proceed in two stages, starting by first constructing initial
estimates of the nuisance parameters present in the EIF, then proceeding to
apply distinct bias-correction strategies in their second stages. Both
estimation strategies require an assumption about the behavior of initial
estimators of the nuisance parameters (specifically, that these lie in a Donsker
class); however, the need for such an assumption may be avoided by making use of
cross-validation in the fitting fo initial estimators. The `medoutcon` `R`
package requires the use of cross-validation in the construction of these
initial estimates, resulting in cross-fitted one-step and and cross-validated
TML estimators [@klaassen1987consistent; @zheng2011cross;
@chernozhukov2018double].

The one-step estimator $\hat{\theta}_{\text{os}}$ is constructed by adding the
empirical mean of the EIF (evaluated at initial estimates of the nuisance
parameters) to the substitution estimator. By constrast, the TML estimator
$\hat{\theta}_{\text{tmle}}$ updates the components of the substitution
estimator via logistic tilting models formulated to ensure that relevant score
equations appearing in the EIF are (approximately) solved.  While the estimators
are asymptotically equivalent, TML estimators have been shown to exhibit
superior finite-sample performance, making them potentially more reliable than
one-step estimators. For the exact form of the EIF as well as those of the
one-step and TML estimators, consult @diaz2020nonparametric.

# Data Analysis Example

## Setting up the data example

Now, we'll take a look at how to estimate the interventional direct and indirect
effects using a simulated data example. @diaz2020nonparametric illustrate the
use of their estimators of these effects in an application in which they seek to
elucidate the mechanisms behind the unintended harmful effects that a housing
intervention had on adolescent girls' risk behavior.

First, let's load a few required packages and set a seed for our simulation.

```{r setup, message=FALSE, warning=FALSE}
library(data.table)
library(medoutcon)
library(sl3)
set.seed(75681)
n_obs <- 500
```

Next, we'll generate a very simple simulated dataset. The function
`make_example_data`, defined below, generates three binary baseline covariates
$W = (W_1, W_2, W_3)$, a binary exposure variable $A$, a single binary mediateor
$M$ an intermediate confounder $Z$ that affects the mediator $M$ and is itself
affected by the exposure $A$, and, finally, a binary outcome $Y$ that is a
function of $(W, A, Z, M)$.

```{r make_example_data, message=FALSE, warning=FALSE}
# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  ## baseline covariates
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure
  a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))

  ## mediator -- could be multivariate
  m <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w[, -3] + a - z)))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + m)))

  ## construct output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

# set seed and simulate example data
example_data <- make_example_data(n_obs)
w_names <- stringr::str_subset(colnames(example_data), "W")
m_names <- stringr::str_subset(colnames(example_data), "M")
```

Now, let's take a quick look at our simulated data:

```{r example_data, message=FALSE, warning=FALSE}
# quick look at the data
head(example_data)
```

As noted above, all covariates in our dataset are binary; however, note that
this need not be the case for using our methodology --- in particular, the only
current limitation is that the intermediate confounder $Z$ must be binary when
using our implemented TML estimator of the (in)direct effects.

Using this dataset, we'll proceed to estimate the interventional (in)direct
effects. In order to do so, we'll need to estimate several nuisance parameters,
including the exposure mechanism $g(A \mid W)$, a re-parameterized exposure
mechanism that conditions on the mediators $h(A \mid M, W)$, the outcome
mechanism $b(Y \mid M, Z, A, W)$, and two variants of the intermediate
confounding mechanism $q(Z \mid A, W)$ and $r(Z \mid M, A, W)$. In order to
estimate each of these nuisance parameters flexibly, we'll rely on data adaptive
regression strategies in order to avoid the potential for (parametric) model
misspecification.

<!--
Note that there are two additional nuisance parameters that must also be
estimated ($u$ and $v$), which are themselves functions of the other nuisance
parameters.  We recommend estimating these via the highly adaptive lasso, which
is the...
-->

## Ensemble learning of nuisance functions

As we'd like to rely on flexible, data adaptive regression strategies for
estimating each of the nuisance parameters $(g, h, b, q, r)$, we require a
method for choosing among or combining the wide variety of available regression
strategies. For this, we recommend the use of the Super Learner algorithm for
ensemble machine learning [@vdl2007super].  The recently developed [`sl3` R
package](https://tlverse.org/sl3) [@coyle2020sl3] provides a unified interface
for deploying a wide variety of machine learning algorithms (simply called
_learners_ in the `sl3` nomenclature) as well as for constructing Super Learner
ensemble models of such learners. For a complete guide on using the `sl3` R
package, consider consulting https://tlverse.org/sl3, or https://tlverse.org
(and https://github.com/tlverse) for the `tlverse` ecosystem, of which `sl3` is
an integral part.

To construct an ensemble learner using a handful of popular machine learning
algorithms, we'll first instantiate variants of learners from the appropriate
classes for each algorithm, and then create a Super Learner ensemble via the
`Lrnr_sl` class. Below, we demonstrate the construction of an ensemble learner
based on a modeling library including an intercept model, a main-terms GLM,
$\ell_1$-penalized Lasso regression, an elastic net regression that equally
weights the $\ell_1$ and $\ell_2$ penalties, random forests (`ranger`), and the
highly adaptive lasso (HAL):

```{r make_sl, message=FALSE, warning=FALSE}
# instantiate learners
mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new(family = binomial())
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, family = "binomial", nfolds = 3)
enet_lrnr <- Lrnr_glmnet$new(alpha = 0.5, family = "binomial", nfolds = 3)
rf_lrnr <- Lrnr_ranger$new(num.trees = 200)

# for HAL, use linear probability formulation, with bounding in unit interval
hal_gaussian_lrnr <- Lrnr_hal9001$new(
  family = "gaussian",
  fit_control = list(
    max_degree = 3,
    n_folds = 3,
    use_min = TRUE,
    type.measure = "mse"
  )
)
bound_lrnr <- Lrnr_bound$new(bound = 1e-6)
hal_bounded_lrnr <- Pipeline$new(hal_gaussian_lrnr, bound_lrnr)

# create learner library and instantiate super learner ensemble
lrnr_lib <- Stack$new(mean_lrnr, fglm_lrnr, enet_lrnr, lasso_lrnr,
                      rf_lrnr, hal_bounded_lrnr)
sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())
```

While we recommend the use of a Super Learner ensemble model like the one
constructed above in practice, such a library will be too computationally
intensive for our examples. To reduce computation time, we construct a simpler
library, using only a subset of the above learning algorithms:

```{r make_simple_sl, message=FALSE, warning=FALSE}
# create simpler learner library and instantiate super learner ensemble
lrnr_lib <- Stack$new(mean_lrnr, fglm_lrnr, lasso_lrnr, rf_lrnr)
sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())
```

Having set up our ensemble learner, we're now ready to estimate each of the
interventional effects using the efficient estimators exposed in the `medoutcon`
package.

## Estimating the direct effect

We're now ready to estimate the interventional direct effect. This direct effect
is computed as a contrast between the interventions $(a' = 1, a^{\star} = 0)$
and $(a' = 0, a^{\star} = 0)$. In particular, our efficient estimators of the
interventional direct effect proceed by constructing estimators
$\hat{\theta}(a' = 1, a^{\star} = 0)$ and $\hat{\theta}(a' = 0, a^{\star} = 0)$.
Then, an efficient estimator of the direct effect is available by application
of the delta method, that is, $\hat{\theta}^{\text{DE}} =
\hat{\theta}(a' = 1, a^{\star} = 0) - \hat{\theta}(a' = 0, a^{\star} = 0)$.
Applying the same principle to the EIF estimates, one can derive variance
estimates and construct asymptotically correct Wald-style confidence intervals
for $\hat{\theta}^{\text{DE}}$.

The `medoutcon` package makes the estimation task quite simple, as only a single
call to the eponymous `medoutcon` function is required. As demonstrated below,
we need only feed in each component of the observed data $O = (W, A, Z, M, Y)$
(of which $W$ and $M$ can be multivariate), specify the effect type, and the
estimator. Additionally, for each nuisance parameter we may specify a separate
regression function --- in the examples below, we use the simpler Super Learner
ensemble constructed above for fitting each nuisance function, but this need not
be the case (i.e., different estimators may be used for each nuisance function).

First, we examine the one-step estimator of the interventional direct effect.
Recall that the one-step estimator is constructed by adding the mean of the EIF
(evaluated at initial estimates of the nuisance parameters) to the substitution
estimator. As noted above, this is done separately for each of the two contrasts
$(a' = 0, a^{\star} = 0)$ and $(a' = 1, a^{\star} = 0)$. Thus, the one-step
estimator of this direct effect is constructed by application of the delta
method to each of the one-step estimators (and EIFs) for these contrasts.

```{r de_os, message=FALSE, warning=FALSE}
# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   g_learners = sl_lrnr,
                   h_learners = sl_lrnr,
                   b_learners = sl_lrnr,
                   q_learners = sl_lrnr,
                   r_learners = sl_lrnr,
                   effect = "direct",
                   estimator = "onestep",
                   estimator_args = list(cv_folds = 2))
summary(os_de)
```

From the output of the summary method, we note that the one-step estimate of
the interventional direct effect $\hat{\theta}_{\text{os}}^{\text{DE}}$ is
`r round(os_de$theta, 3)`, with 95% confidence interval
[`r round(unname(confint(os_de, 0.95)[1]), 3)`,
`r round(unname(confint(os_de, 0.95)[3]), 3)`].

Next, let's compare the one-step estimate to the TML estimate. Analogous to the
case of the one-step estimator, the TML estimator can be evaluated via a single
call to the `medoutcon` function:

```{r de_tmle, message=FALSE, warning=FALSE}
# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     g_learners = sl_lrnr,
                     h_learners = sl_lrnr,
                     b_learners = sl_lrnr,
                     q_learners = sl_lrnr,
                     r_learners = sl_lrnr,
                     effect = "direct",
                     estimator = "tmle",
                     estimator_args = list(cv_folds = 2, max_iter = 5))
summary(tmle_de)
```

From the output of the summary method, we note that the TML estimate of the
interventional direct effect $\hat{\theta}_{\text{tmle}}^{\text{DE}}$ is
`r round(tmle_de$theta, 3)`, with 95% confidence interval
[`r round(unname(confint(tmle_de, 0.95)[1]), 3)`,
`r round(unname(confint(tmle_de, 0.95)[3]), 3)`]. Here, we recall that the TML
estimator generally exhibits better finite-sample performance than the one-step
estimator [@vdl2011targeted; @vdl2018targeted], so the TML estimate is likely to
be more reliable in our modest sample size of $n =$ `r n_obs`.

## Estimating the indirect effect

Estimation of the interventional indirect effect proceeds similarly to the
strategy discussed above for the corresponding direct effect. An efficient
estimator can be computed as a contrast between the interventions $(a' = 1,
a^{\star} = 0)$ and $(a' = 1, a^{\star} = 1)$. Specifically, our efficient
estimators of the interventional indirect effect proceed by constructing
estimators $\hat{\theta}(a' = 1, a^{\star} = 0)$ and $\hat{\theta}(a' = 1,
a^{\star} = 1)$.  Then, application of the delta method yields an efficient
estimator of the indirect effect, that is, $\hat{\theta}^{\text{IE}} =
\hat{\theta}(a' = 1, a^{\star} = 0) - \hat{\theta}(a' = 1, a^{\star} = 1)$. The
same principle may be applied to the EIF estimates to derive variance estimates
and construct asymptotically correct Wald-style confidence intervals for
$\hat{\theta}^{\text{IE}}$.

Now, we examine the one-step estimator of the interventional indirect effect.
The one-step estimator is constructed by adding the mean of the EIF
(evaluated at initial estimates of the nuisance parameters) to the substitution
estimator. As noted above, this is done separately for each of the two contrasts
$(a' = 1, a^{\star} = 1)$ and $(a' = 1, a^{\star} = 0)$. Thus, the one-step
estimator of this indirect effect is constructed by application of the delta
method to each of the one-step estimators (and EIFs) for the contrasts.

```{r ie_os, message=FALSE, warning=FALSE}
# compute one-step estimate of the interventional indirect effect
os_ie <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   g_learners = sl_lrnr,
                   h_learners = sl_lrnr,
                   b_learners = sl_lrnr,
                   q_learners = sl_lrnr,
                   r_learners = sl_lrnr,
                   effect = "indirect",
                   estimator = "onestep")
summary(os_ie)
```

From the output of the summary method, we note that the one-step estimate of
the interventional indirect effect $\hat{\theta}_{\text{os}}^{\text{IE}}$ is
`r round(os_ie$theta, 3)`, with 95% confidence interval
[`r round(unname(confint(os_ie, 0.95)[1]), 3)`,
`r round(unname(confint(os_ie, 0.95)[3]), 3)`].

As before, let's compare the one-step estimate to the TML estimate. Analogous
to the case of the one-step estimator, the TML estimator can be evaluated via a
single call to the `medoutcon` function, as demonstrated below

```{r ie_tmle, message=FALSE, warning=FALSE}
# compute targeted minimum loss estimate of the interventional indirect effect
tmle_ie <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     g_learners = sl_lrnr,
                     h_learners = sl_lrnr,
                     b_learners = sl_lrnr,
                     q_learners = sl_lrnr,
                     r_learners = sl_lrnr,
                     effect = "indirect",
                     estimator = "tmle")
summary(tmle_ie)
```

From the output of the summary method, we note that the TML estimate of the
interventional indirect effect $\hat{\theta}_{\text{tmle}}^{\text{IE}}$ is
`r round(tmle_ie$theta, 3)`, with 95% confidence interval
[`r round(unname(confint(tmle_ie, 0.95)[1]), 3)`,
`r round(unname(confint(tmle_ie, 0.95)[3]), 3)`]. As before, the TML estimator
provides better finite-sample performance than the one-step estimator, so it may
be preferred in this example.

## References

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{cv_eif}
\alias{cv_eif}
\title{EIF for stochastic interventional (in)direct effects}
\usage{
cv_eif(
  fold,
  data_in,
  contrast,
  g_learners,
  h_learners,
  b_learners,
  q_learners,
  r_learners,
  u_learners,
  v_learners,
  effect_type = c("interventional", "natural"),
  w_names,
  m_names
)
}
\arguments{
\item{fold}{Object specifying cross-validation folds as generated by a call
to \code{\link[origami]{make_folds}}.}

\item{data_in}{A \code{data.table} containing the observed data with columns
are in the order specified by the NPSEM (Y, M, Z, A, W), with column names
set appropriately based on the input data. Such a structure is merely a
convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{NULL} has
no effect, as the value of the argument \code{effect} is instead used to
define the contrasts. To override \code{effect}, provide a \code{numeric}
double vector, giving the values of a' and a*, e.g., \code{c(0, 1)}.}

\item{g_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the propensity score.}

\item{h_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a parameterization of the
propensity score that conditions on the mediators.}

\item{b_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the outcome regression.}

\item{q_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the treatment and potential
baseline covariates.}

\item{r_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the mediators, the treatment,
and potential baseline confounders.}

\item{u_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{v_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{effect_type}{A \code{character} indicating whether components of the
interventional or natural (in)direct effects are to be estimated. In the
case of the natural (in)direct effects, estimation of several nuisance
parameters is unnecessary.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by \code{\link{medoutcon}}.}
}
\description{
EIF for stochastic interventional (in)direct effects
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{fit_nuisance_v}
\alias{fit_nuisance_v}
\title{Fit pseudo-outcome regression conditioning on treatment and baseline}
\usage{
fit_nuisance_v(
  train_data,
  valid_data,
  contrast,
  learners,
  b_out,
  q_out,
  m_names,
  w_names
)
}
\arguments{
\item{train_data}{A \code{data.table} containing observed data, with columns
in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
appropriately based on the original input data. Such a structure is merely
a convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{valid_data}{A holdout data set, with columns exactly matching those
appearing in the preceding argument \code{data}, to be used for estimation
via cross-fitting. Not optional for this nuisance parameter.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{c(0, 1)}
assumes a binary intervention node \code{A}.}

\item{learners}{\code{\link[sl3]{Stack}}, or other learner class (inheriting
from \code{\link[sl3]{Lrnr_base}}), containing a set of learners from
\pkg{sl3}, to be used in fitting a model for this nuisance parameter.}

\item{b_out}{Output from the internal function for fitting the outcome
regression \code{\link{fit_out_mech}}.}

\item{q_out}{Output from the internal function for fitting the mechanism of
the intermediate confounder while conditioning on the mediators, i.e.,
\code{\link{fit_moc_mech}}, setting \code{type = "q"}.}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by a call to the wrapper function \code{\link{medoutcon}}.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}
}
\description{
Fit pseudo-outcome regression conditioning on treatment and baseline
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{fit_moc_mech}
\alias{fit_moc_mech}
\title{Fit intermediate confounding mechanism with(out) conditioning on mediators}
\usage{
fit_moc_mech(
  train_data,
  valid_data = NULL,
  contrast,
  learners,
  m_names,
  w_names,
  type = c("q", "r")
)
}
\arguments{
\item{train_data}{A \code{data.table} containing observed data, with columns
in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
appropriately based on the original input data. Such a structure is merely
a convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{valid_data}{A holdout data set, with columns exactly matching those
appearing in the preceding argument \code{data}, to be used for estimation
via cross-fitting. Optional, defaulting to \code{NULL}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{c(0, 1)}
assumes a binary intervention node \code{A}.}

\item{learners}{\code{\link[sl3]{Stack}}, or other learner class (inheriting
from \code{\link[sl3]{Lrnr_base}}), containing a set of learners from
\pkg{sl3}, to be used in fitting a model for the intermediate confounding
mechanism, i.e., q = E[z|a',W] and r = E[z|a',m,w]).}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by a call to the wrapper function \code{\link{medoutcon}}.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}

\item{type}{A \code{character} vector indicating whether to condition on the
mediators (M) or not. Specifically, this is an option for specifying one of
two types of nuisance regressions: "r" is defined as the component that
conditions on the mediators (i.e., r = E[z|a',m,w]) while "q" is defined as
the component that does not (i.e., q = E[z|a',w]).}
}
\description{
Fit intermediate confounding mechanism with(out) conditioning on mediators
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{est_onestep}
\alias{est_onestep}
\title{One-step estimator for stochastic interventional (in)direct effects}
\usage{
est_onestep(
  data,
  contrast,
  g_learners,
  h_learners,
  b_learners,
  q_learners,
  r_learners,
  u_learners,
  v_learners,
  w_names,
  m_names,
  y_bounds,
  effect_type = c("interventional", "natural"),
  svy_weights = NULL,
  cv_folds = 5L
)
}
\arguments{
\item{data}{A \code{data.table} containing the observed data, with columns
in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
appropriately based on the original input data. Such a structure is merely
a convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{NULL} has
no effect, as the value of the argument \code{effect} is instead used to
define the contrasts. To override \code{effect}, provide a \code{numeric}
double vector, giving the values of a' and a*, e.g., \code{c(0, 1)}.}

\item{g_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the propensity score.}

\item{h_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a parameterization of the
propensity score that conditions on the mediators.}

\item{b_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the outcome regression.}

\item{q_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the treatment and potential
baseline covariates.}

\item{r_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the mediators, the treatment,
and potential baseline confounders.}

\item{u_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{v_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by \code{\link{medoutcon}}.}

\item{y_bounds}{A \code{numeric} double indicating the minimum and maximum
observed values of the outcome variable Y prior to its being re-scaled to
the unit interval.}

\item{effect_type}{A \code{character} indicating whether components of the
interventional or natural (in)direct effects are to be estimated. In the
case of the natural (in)direct effects, estimation of several nuisance
parameters is unnecessary.}

\item{svy_weights}{A \code{numeric} vector of observation-level weights that
have been computed externally. Such weights are used in the construction of
a re-weighted estimator.}

\item{cv_folds}{A \code{numeric} integer specifying the number of folds to
be created for cross-validation. Use of cross-validation allows for entropy
conditions on the one-step estimator to be relaxed. For compatibility with
\code{\link[origami]{make_folds}}, this value specified must be greater
than or equal to 2; the default is to create 5 folds.}
}
\description{
One-step estimator for stochastic interventional (in)direct effects
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{est_tml}
\alias{est_tml}
\title{TML estimator for stochastic interventional (in)direct effects}
\usage{
est_tml(
  data,
  contrast,
  g_learners,
  h_learners,
  b_learners,
  q_learners,
  r_learners,
  u_learners,
  v_learners,
  w_names,
  m_names,
  y_bounds,
  effect_type = c("interventional", "natural"),
  svy_weights = NULL,
  cv_folds = 5L,
  max_iter = 5L,
  tiltmod_tol = 10
)
}
\arguments{
\item{data}{A \code{data.table} containing the observed data, with columns
in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
appropriately based on the original input data. Such a structure is merely
a convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{NULL} has
no effect, as the value of the argument \code{effect} is instead used to
define the contrasts. To override \code{effect}, provide a \code{numeric}
double vector, giving the values of a' and a*, e.g., \code{c(0, 1)}.}

\item{g_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the propensity score.}

\item{h_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a parameterization of the
propensity score that conditions on the mediators.}

\item{b_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the outcome regression.}

\item{q_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the treatment and potential
baseline covariates.}

\item{r_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the mediators, the treatment,
and potential baseline confounders.}

\item{u_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{v_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by \code{\link{medoutcon}}.}

\item{y_bounds}{A \code{numeric} double indicating the minimum and maximum
observed values of the outcome variable Y prior to its being re-scaled to
the unit interval.}

\item{effect_type}{A \code{character} indicating whether components of the
interventional or natural (in)direct effects are to be estimated. In the
case of the natural (in)direct effects, estimation of several nuisance
parameters is unnecessary.}

\item{svy_weights}{A \code{numeric} vector of observation-level weights that
have been computed externally. Such weights are used in the construction of
a re-weighted estimator.}

\item{cv_folds}{A \code{numeric} value specifying the number of folds to be
created for cross-validation. Use of cross-validation allows for entropy
conditions on the TML estimator to be relaxed. Note: for compatibility with
\code{\link[origami]{make_folds}}, this value  must be greater than or
equal to 2; the default is to create 10 folds.}

\item{max_iter}{A \code{numeric} integer giving the maximum number of steps
to be taken for the iterative procedure to construct a TML estimator.}

\item{tiltmod_tol}{A \code{numeric} indicating the maximum step size to be
taken when performing TMLE updates based on logistic tilting models. When
the step size of a given update exceeds this value, the update is avoided.}
}
\description{
TML estimator for stochastic interventional (in)direct effects
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{scale_to_unit}
\alias{scale_to_unit}
\title{Scale values to the unit interval}
\usage{
scale_to_unit(vals)
}
\arguments{
\item{vals}{A \code{numeric} vector of values to be scaled into the closed
interval [0, 1].}
}
\description{
Scale values to the unit interval
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print.medoutcon}
\alias{print.medoutcon}
\title{Print method for interventional mediation effect estimate objects}
\usage{
\method{print}{medoutcon}(x, ...)
}
\arguments{
\item{x}{An object of class \code{medoutcon}.}

\item{...}{Other options (not currently used).}
}
\description{
The \code{print} method for objects of class \code{medoutcon}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bound_precision}
\alias{bound_precision}
\title{Bounding to numerical precision}
\usage{
bound_precision(vals, tol = 1e-06)
}
\arguments{
\item{vals}{A \code{numeric} vector of values in the interval [0, 1].}

\item{tol}{A \code{numeric} indicating the tolerance limit to which extreme
values should be truncated. Realizations of \code{val} less than \code{tol}
are truncated to \code{tol} while those greater than (1 - \code{tol}) are
truncated to (1 - \code{tol}).}
}
\description{
Bounds extreme values to a specified tolerance level, for use with sensitive
quantities that must be transformed, e.g., via \code{\link[stats]{qlogis}}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{fit_treat_mech}
\alias{fit_treat_mech}
\title{Fit propensity scores for treatment contrasts}
\usage{
fit_treat_mech(
  train_data,
  valid_data = NULL,
  contrast,
  learners,
  m_names,
  w_names,
  type = c("g", "h")
)
}
\arguments{
\item{train_data}{A \code{data.table} containing the observed data; columns
are in the order specified by the NPSEM (Y, M, Z, A, W), with column names
set appropriately based on the input data. Such a structure is merely a
convenience utility to passing data around to the various core estimation
routines and is automatically generated \code{\link{medoutcon}}.}

\item{valid_data}{A holdout data set, with columns exactly matching those
appearing in the preceding argument \code{train_data}, to be used for
estimation via cross-fitting. Optional, defaulting to \code{NULL}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{c(0, 1)}
assumes a binary intervention node \code{A}.}

\item{learners}{\code{\link[sl3]{Stack}}, or other learner class (inheriting
from \code{\link[sl3]{Lrnr_base}}), containing a set of learners from
\pkg{sl3}, to be used in fitting a propensity score models, i.e., g :=
P(A = 1 | W) and h := P(A = 1 | M, W).}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by \code{\link{medoutcon}}.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}

\item{type}{A \code{character} indicating which of the treatment mechanism
variants to estimate. Option \code{"g"} corresponds to the propensity score
g(A|W) while option \code{"h"} conditions on the mediators h(A|M,W).}
}
\description{
Fit propensity scores for treatment contrasts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medoutcon.R
\name{medoutcon}
\alias{medoutcon}
\title{Efficient estimation of stochastic interventional (in)direct effects}
\usage{
medoutcon(
  W,
  A,
  Z,
  M,
  Y,
  obs_weights = rep(1, length(Y)),
  svy_weights = NULL,
  effect = c("direct", "indirect"),
  contrast = NULL,
  g_learners = sl3::Lrnr_glm_fast$new(),
  h_learners = sl3::Lrnr_glm_fast$new(),
  b_learners = sl3::Lrnr_glm_fast$new(),
  q_learners = sl3::Lrnr_glm_fast$new(),
  r_learners = sl3::Lrnr_glm_fast$new(),
  u_learners = sl3::Lrnr_hal9001$new(),
  v_learners = sl3::Lrnr_hal9001$new(),
  estimator = c("tmle", "onestep"),
  estimator_args = list(cv_folds = 5L, max_iter = 5L, tiltmod_tol = 10)
)
}
\arguments{
\item{W}{A \code{matrix}, \code{data.frame}, or similar object corresponding
to a set of baseline covariates.}

\item{A}{A \code{numeric} vector corresponding to a treatment variable. The
parameter of interest is defined as a location shift of this quantity.}

\item{Z}{A \code{numeric} vector corresponding to an intermediate confounder
affected by treatment (on the causal pathway between the intervention A,
mediators M, and outcome Y, but unaffected itself by the mediators).}

\item{M}{A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
similar corresponding to a set of mediators (on the causal pathway between
the intervention A and the outcome Y).}

\item{Y}{A \code{numeric} vector corresponding to an outcome variable.}

\item{obs_weights}{A \code{numeric} vector of observation-level weights.
The default is to give all observations equal weighting.}

\item{svy_weights}{A \code{numeric} vector of observation-level weights that
have been computed externally, such as survey sampling weights. Such
weights are used in the construction of re-weighted efficient estimators.}

\item{effect}{A \code{character} indicating whether to compute the direct
or the indirect effect as discussed in <https://arxiv.org/abs/1912.09936>.
This is ignored when the argument \code{contrast} is provided. By default,
the direct effect is estimated.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default value of \code{NULL} has
no effect, as the value of the argument \code{effect} is instead used to
define the contrasts. To override \code{effect}, provide a \code{numeric}
double vector, giving the values of a' and a*, e.g., \code{c(0, 1)}.}

\item{g_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the propensity score.}

\item{h_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a parameterization of the
propensity score that conditions on the mediators.}

\item{b_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for the outcome regression.}

\item{q_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the treatment and potential
baseline covariates.}

\item{r_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a model for a nuisance regression of
the intermediate confounder, conditioning on the mediators, the treatment,
and potential baseline confounders.}

\item{u_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{v_learners}{A \code{\link[sl3]{Stack}} object, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
for in the efficient influence function.}

\item{estimator}{The desired estimator of the direct or indirect effect (or
contrast-specific parameter) to be computed. Both an efficient one-step
estimator using cross-fitting and a cross-validated targeted minimum loss
estimator (TMLE) are available. The default is the TML estimator.}

\item{estimator_args}{A \code{list} of extra arguments to be passed (via
\code{...}) to the function call for the specified estimator. The default
is chosen so as to allow the number of folds used in computing the one-step
or TML estimators to be easily adjusted. In the case of the TML estimator,
the number of update (fluctuation) iterations is limited, and a tolerance
is included for the updates introduced by the tilting (fluctuation) models.}
}
\description{
Efficient estimation of stochastic interventional (in)direct effects
}
\examples{
# here, we show one-step and TML estimates of the interventional direct
# effect; the indirect effect can be evaluated by a straightforward change
# to the penultimate argument. the natural direct and indirect effects can
# be evaluated by omitting the argument Z (inappropriate in this example).
# create data: covariates W, exposure A, post-exposure-confounder Z,
#              mediator M, outcome Y
n_obs <- 200
w_1 <- rbinom(n_obs, 1, prob = 0.6)
w_2 <- rbinom(n_obs, 1, prob = 0.3)
w <- as.data.frame(cbind(w_1, w_2))
a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))
z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))
m_1 <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w + a - z)))
m_2 <- rbinom(n_obs, 1, plogis(rowSums(w - a - z)))
m <- as.data.frame(cbind(m_1, m_2))
y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + rowSums(m))))

# one-step estimate of the interventional direct effect
os_de <- medoutcon(
  W = w, A = a, Z = z, M = m, Y = y,
  effect = "direct",
  estimator = "onestep"
)

# TML estimate of the interventional direct effect
# NOTE: improved variance estimate and de-biasing from targeting procedure
tmle_de <- medoutcon(
  W = w, A = a, Z = z, M = m, Y = y,
  effect = "direct",
  estimator = "tmle"
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{summary.medoutcon}
\alias{summary.medoutcon}
\title{Summary for interventional mediation effect estimate objects}
\usage{
\method{summary}{medoutcon}(object, ..., ci_level = 0.95)
}
\arguments{
\item{object}{An object of class \code{medoutcon}, as produced by invoking
\code{\link{medoutcon}}.}

\item{...}{Other arguments. Not currently used.}

\item{ci_level}{A \code{numeric} indicating the level of the confidence
interval to be computed.}
}
\description{
Print a convenient summary for objects of \code{S3} class \code{medoutcon}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{scale_from_unit}
\alias{scale_from_unit}
\title{Scale values from the unit interval to their original scale}
\usage{
scale_from_unit(scaled_vals, max_orig, min_orig)
}
\arguments{
\item{scaled_vals}{A \code{numeric} vector of values scaled to lie in the
closed interval [0, 1] by use of \code{\link{scale_to_unit}}.}

\item{max_orig}{The maximum of the values on the original scale.}

\item{min_orig}{The minimum of the values on the original scale.}
}
\description{
Scale values from the unit interval to their original scale
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bound_propensity}
\alias{bound_propensity}
\title{Bounding propensity scores}
\usage{
bound_propensity(vals, bounds = c(0.001, 0.999))
}
\arguments{
\item{vals}{A \code{numeric} vector of values in the interval [0, 1].}

\item{bounds}{A \code{numeric} vector containing two values, the first being
the minimum allowable value and the second being the maximum allowable for
values appearing in the vector \code{vals} (the previous argument).}
}
\description{
Bounds estimated propensity score values to be within a specified range.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{fit_nuisance_u}
\alias{fit_nuisance_u}
\title{Fit pseudo-outcome regression conditioning on mediator-outcome confounder}
\usage{
fit_nuisance_u(
  train_data,
  valid_data,
  learners,
  b_out,
  q_out,
  r_out,
  g_out,
  h_out,
  w_names
)
}
\arguments{
\item{train_data}{A \code{data.table} containing observed data, with columns
in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
appropriately based on the original input data. Such a structure is merely
a convenience utility to passing data around to the various core estimation
routines and is automatically generated by \code{\link{medoutcon}}.}

\item{valid_data}{A holdout data set, with columns exactly matching those
appearing in the preceding argument \code{data}, to be used for estimation
via cross-fitting. NOT optional for this nuisance parameter.}

\item{learners}{\code{\link[sl3]{Stack}}, or other learner class (inheriting
from \code{\link[sl3]{Lrnr_base}}), containing a set of learners from
\pkg{sl3}, to be used in fitting a model for this nuisance parameter.}

\item{b_out}{Output from the internal function for fitting the outcome
regression \code{\link{fit_out_mech}}.}

\item{q_out}{Output from the internal function for fitting the mechanism of
the intermediate confounder while conditioning on mediators, i.e.,
\code{\link{fit_moc_mech}}, setting \code{type = "q"}.}

\item{r_out}{Output from the internal function for fitting the mechanism of
the intermediate confounder without conditioning on mediators, i.e.,
\code{\link{fit_moc_mech}}, setting \code{type = "r"}.}

\item{g_out}{Output from the internal function for fitting the treatment
mechanism without conditioning on mediators \code{\link{fit_treat_mech}}.}

\item{h_out}{Output from the internal function for fitting the treatment
mechanism conditioning on the mediators \code{\link{fit_treat_mech}}.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}
}
\description{
Fit pseudo-outcome regression conditioning on mediator-outcome confounder
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mechanisms.R
\name{fit_out_mech}
\alias{fit_out_mech}
\title{Fit outcome regression}
\usage{
fit_out_mech(
  train_data,
  valid_data = NULL,
  contrast,
  learners,
  m_names,
  w_names
)
}
\arguments{
\item{train_data}{A \code{data.table} containing the observed data, with
columns in the order specified by the NPSEM (Y, M, Z, A, W), with column
names set based on the original input data. Such a structure is merely a
convenience utility to passing data around to the various core estimation
routines and is automatically generated \code{\link{medoutcon}}.}

\item{valid_data}{A holdout data set, with columns exactly matching those
appearing in the preceding argument \code{data}, to be used for estimation
via cross-fitting. Optional, defaulting to \code{NULL}.}

\item{contrast}{A \code{numeric} double indicating the two values of the
intervention \code{A} to be compared. The default of \code{c(0, 1)} assumes
a binary intervention node \code{A}.}

\item{learners}{\code{\link[sl3]{Stack}}, or other learner class (inheriting
from \code{\link[sl3]{Lrnr_base}}), containing a set of learners from
\pkg{sl3}, to be used in fitting the outcome regression, i.e., b(A,Z,M,W).}

\item{m_names}{A \code{character} vector of the names of the columns that
correspond to mediators (M). The input for this argument is automatically
generated by \code{\link{medoutcon}}.}

\item{w_names}{A \code{character} vector of the names of the columns that
correspond to baseline covariates (W). The input for this argument is
automatically generated by \code{\link{medoutcon}}.}
}
\description{
Fit outcome regression
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{confint.medoutcon}
\alias{confint.medoutcon}
\title{Confidence intervals for interventional mediation effect estimates}
\usage{
\method{confint}{medoutcon}(object, parm = seq_len(object$theta), level = 0.95, ...)
}
\arguments{
\item{object}{An object of class \code{medoutcon}, as produced by invoking
\code{\link{medoutcon}}, for which a confidence interval is to be computed.}

\item{parm}{A \code{numeric} vector indicating indices of \code{object$est}
for which to return confidence intervals.}

\item{level}{A \code{numeric} indicating the level of the confidence
interval to be computed.}

\item{...}{Other arguments. Not currently used.}
}
\description{
Compute confidence intervals for objects of class \code{medoutcon}, which
contain estimates produced by \code{\link{medoutcon}}.
}
