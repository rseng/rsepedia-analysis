
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

