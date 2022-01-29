# mcboost

<!-- badges: start -->
[![tic](https://github.com/mlr-org/mcboost/workflows/tic/badge.svg?branch=main)](https://github.com/mlr-org/mcboost/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN Status](https://www.r-pkg.org/badges/version-ago/mcboost)](https://cran.r-project.org/package=mcboost)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03453/status.svg)](https://doi.org/10.21105/joss.03453)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Mattermost](https://img.shields.io/badge/chat-mattermost-orange.svg)](https://lmmisld-lmu-stats-slds.srv.mwn.de/mlr_invite/)
<!-- badges: end -->

## What does it do?

**mcboost** implements Multi-Calibration Boosting ([Hebert-Johnson et al., 2018](https://proceedings.mlr.press/v80/hebert-johnson18a.html); [Kim et al., 2019](https://arxiv.org/pdf/1805.12317.pdf)) for the multi-calibration of a machine learning model's prediction. Multi-Calibration works best in scenarios where the underlying data & labels are unbiased but a bias is introduced within the algorithm's fitting procedure. This is often the case, e.g. when an algorithm fits a majority population while ignoring or under-fitting minority populations.

For more information and example, see the package's [website](https://mlr-org.github.io/mcboost/).

More details with respect to usage and the procedures can be found in the package vignettes.

## Installation

The current version can be downloaded from CRAN using:

```r
install.packages("mcboost")
```

You can install the development version of mcboost from **Github** with:

```r
remotes::install_github("mlr-org/mcboost")
```

## Usage

Post-processing with `mcboost` needs three components. We start with an initial prediction model (1) and an auditing algorithm (2) that may be customized by the user. The auditing algorithm then runs Multi-Calibration-Boosting on a labeled auditing dataset (3). The resulting model can be used for obtaining multi-calibrated predictions.

<p align="center">
  <img src="https://github.com/mlr-org/mcboost/raw/main/paper/MCBoost.png" />
</p>

## Example

In this simple example, our goal is to improve calibration
for an `initial predictor`, e.g. a ML algorithm trained on
an initial task.
Internally, `mcboost` often makes use of `mlr3` and learners that come with `mlr3learners`.


``` r
library(mcboost)
library(mlr3)
```

First we set up an example dataset.

```r
  #  Example Data: Sonar Task
  tsk = tsk("sonar")
  tid = sample(tsk$row_ids, 100) # 100 rows for training
  train_data = tsk$data(cols = tsk$feature_names, rows = tid)
  train_labels = tsk$data(cols = tsk$target_names, rows = tid)[[1]]
```

To provide an example, we assume that we have already a learner `l` which we train below.
We can now wrap this initial learner's predict function for use with `mcboost`, since `mcboost` expects the initial model to be specified as a `function` with `data` as input.

```r
  l = lrn("classif.rpart")
  l$train(tsk$clone()$filter(tid))

  init_predictor = function(data) {
    # Get response prediction from Learner
    p = l$predict_newdata(data)$response
    # One-hot encode and take first column
    one_hot(p)
  }
```

We can now run Multi-Calibration Boosting by instantiating the object and calling the `multicalibrate` method.
Note, that typically, we would use Multi-Calibration on a separate validation set!
We furthermore select the auditor model, a `SubpopAuditorFitter`,
in our case a `Decision Tree`:

```r
  mc = MCBoost$new(
    init_predictor = init_predictor,
    auditor_fitter = "TreeAuditorFitter")
  mc$multicalibrate(train_data, train_labels)
```

Lastly, we predict on new data.

```r
tstid = setdiff(tsk$row_ids, tid) # held-out data
test_data = tsk$data(cols = tsk$feature_names, rows = tstid)
mc$predict_probs(test_data)
```

### Multi-Calibration

While `mcboost` in its defaults implements Multi-Accuracy ([Kim et al., 2019](https://arxiv.org/pdf/1805.12317.pdf)),
it can also multi-calibrate predictors ([Hebert-Johnson et al., 2018](http://proceedings.mlr.press/v80/hebert-johnson18a.html)).
In order to achieve this, we have to set the following hyperparameters:

```r
  mc = MCBoost$new(
    init_predictor = init_predictor,
    auditor_fitter = "TreeAuditorFitter",
    num_buckets = 10,
    multiplicative = FALSE
  )
```

## MCBoost as a PipeOp

`mcboost` can also be used within a `mlr3pipeline` in order to use at the full end-to-end pipeline (in the form of a `GraphLearner`).

```r
  library(mlr3)
  library(mlr3pipelines)
  gr = ppl_mcboost(lrn("classif.rpart"))
  tsk = tsk("sonar")
  tid = sample(1:208, 108)
  gr$train(tsk$clone()$filter(tid))
  gr$predict(tsk$clone()$filter(setdiff(1:208, tid)))
```



## Further Examples

The `mcboost` vignettes [**Basics and Extensions**](https://mlr-org.github.io/mcboost/articles/mcboost_basics_extensions.html) and [**Health Survey Example**](https://mlr-org.github.io/mcboost/articles/mcboost_example.html) demonstrate a lot of interesting showcases for applying `mcboost`.


## Contributing

This R package is licensed under the LGPL-3.
If you encounter problems using this software (lack of documentation, misleading or wrong documentation, unexpected behaviour, bugs, …) or just want to suggest features, please open an issue in the issue tracker.
Pull requests are welcome and will be included at the discretion of the maintainers.

As this project is developed with [mlr3's](https://github.com/mlr-org/mlr3/) style guide in mind, the following resources can be helpful
to individuals wishing to contribute: Please consult the [wiki](https://github.com/mlr-org/mlr3/wiki/) for a [style guide](https://github.com/mlr-org/mlr3/wiki/Style-Guide), a [roxygen guide](https://github.com/mlr-org/mlr3/wiki/Roxygen-Guide) and a [pull request guide](https://github.com/mlr-org/mlr3/wiki/PR-Guidelines).

### Code of Conduct

Please note that the mcboost project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## Citing mcboost

If you use `mcboost`, please cite our package as well as the two papers it is based on:

```
  @article{pfisterer2021,
    author = {Pfisterer, Florian and Kern, Christoph and Dandl, Susanne and Sun, Matthew and 
    Kim, Michael P. and Bischl, Bernd},
    title = {mcboost: Multi-Calibration Boosting for R},
    journal = {Journal of Open Source Software},
    doi = {10.21105/joss.03453},
    url = {https://doi.org/10.21105/joss.03453},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {64},
    pages = {3453}
  }
  # Multi-Calibration
  @inproceedings{hebert-johnson2018,
    title = {Multicalibration: Calibration for the ({C}omputationally-Identifiable) Masses},
    author = {Hebert-Johnson, Ursula and Kim, Michael P. and Reingold, Omer and Rothblum, Guy},
    booktitle = {Proceedings of the 35th International Conference on Machine Learning},
    pages = {1939--1948},
    year = {2018},
    editor = {Jennifer Dy and Andreas Krause},
    volume = {80},
    series = {Proceedings of Machine Learning Research},
    address = {Stockholmsmässan, Stockholm Sweden},
    publisher = {PMLR}
  }
  # Multi-Accuracy
  @inproceedings{kim2019,
    author = {Kim, Michael P. and Ghorbani, Amirata and Zou, James},
    title = {Multiaccuracy: Black-Box Post-Processing for Fairness in Classification},
    year = {2019},
    isbn = {9781450363242},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    url = {https://doi.org/10.1145/3306618.3314287},
    doi = {10.1145/3306618.3314287},
    booktitle = {Proceedings of the 2019 AAAI/ACM Conference on AI, Ethics, and Society},
    pages = {247--254},
    location = {Honolulu, HI, USA},
    series = {AIES '19}
  }
```
# mcboost (development version)

# mcboost 0.4.0
* [Experimental] mcboost now has experimental support for *survival* tasks.
  See `MCBoostSurv` and the corresponding vignette "MCBoostSurv - Basics" for more information.
* We have published an article about mcboost in the Journal of Open Source Software: "https://joss.theoj.org/papers/10.21105/joss.03453". See citation("mcboost") for the citation info.


# mcboost 0.3.3
* Auditors can now also update weights if correlations are negative by switching the sign of the update direction as intended in the paper.
* Minor adaptions to improve stability of unit tests

# mcboost 0.3.2
* Minor adpations to improve stability of unit tests

# mcboost 0.3.1

* Fixed a bug for additive weight updates, were updates went
  in the wrong direction.
* Added new parameter `eval_fulldata` that allows to compute
  auditor effect across the full sample (as opposed to the bucket).

# mcboost 0.3.0

* First CRAN-ready version of the package.
* Added a `NEWS.md` file to track changes to the package.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
## Reason for resubmission

- Update citation file
- Minor bug fixes and refactoring
- New functionality: mcboost for survival

## R CMD check

- No NOTEs, WARNINGs or ERRORs

## R-HUB

- All r-hub checks pass without NOTEs, WARNINGs or ERRORs.
- PREPERROR likely due to https://github.com/r-hub/rhub/issues/448. 
  All checks show "Status: success"### Florian Pfisterer

Florian Pfisterer implemented and extended the main part of the package, heavily influcenced by
an unpublished python code-base written in large parts by Matthew. Florian furthermore worked on the interaction with
**mlr3** by integrating **mlr3** learners as Auditing Mechanism as well as exporting functionality to integrate
**mcboost** as a `PipeOp` into **mlr3pipelines**.

### Christoph Kern

Christoph Kern prepared and contributed to the vignettes and co-authored the summary paper. Christoph contributed (very) moderately to the python code underlying this package and helped conceptionally in transitioning from the python code to the R implementation.

### Susanne Dandl

Susanne Dandl reviewed the R package, she provided advice on extensions, extended and improved vignettes
and worked towards thorough unit testing of the different methods.

### Matthew Sun

Matthew Sun wrote the initial Python implementation of MCBoost, with feedback and oversight provided by Michael.
His version guided large parts of mcboost's current design and architecture.

### Michael P. Kim

Michael P. Kim is a coauthor of the research papers that introduced Multi-Calibration.
Michael oversaw the development of the initial python implementation of MCBoost
and provided additional advice and directions in the development of this R package.

### Bernd Bischl

Oversaw the package development and provided feedback with respect to API design, implementation details and methodology.
---
title: 'mcboost: Multi-Calibration Boosting for R'
tags:
  - R
  - Multi-Calibration
  - Multi-Accuracy
  - Boosting
  - Post-Processing
  - Fair ML
authors:
  - name: Florian Pfisterer^[Corresponding author]
    orcid: 0000-0001-8867-762X
    affiliation: 1
  - name: Christoph Kern
    orcid: 0000-0001-7363-4299
    affiliation: 2
  - name: Susanne Dandl
    orcid: 0000-0003-4324-4163
    affiliation: 1
  - name: Matthew Sun
    affiliation: 3
  - name: Michael P. Kim
    affiliation: 4
  - name: Bernd Bischl
    orcid: 0000-0001-6002-6980
    affiliation: 1
affiliations:
 - name: Ludwig Maximilian University of Munich
   index: 1
 - name: University of Mannheim
   index: 2
 - name: Princeton University
   index: 3
 - name: UC Berkeley
   index: 4
date: 01 June 2021
bibliography: paper.bib
---

# Summary

Given the increasing usage of automated prediction systems in the context of high-stakes decisions, a growing body of research focuses on methods for detecting and mitigating biases in algorithmic decision-making.
One important framework to audit for and mitigate biases in predictions is that of Multi-Calibration, introduced by @hebert-johnson2018.
The underlying fairness notion, Multi-Calibration, promotes the idea of multi-group fairness and requires calibrated predictions not only for marginal populations, but also for subpopulations that may be defined by complex intersections of many attributes.
A simpler variant of Multi-Calibration, referred to as Multi-Accuracy, requires unbiased predictions for large collections of subpopulations.
@hebert-johnson2018 proposed a boosting-style algorithm for learning multi-calibrated predictors.
@kim2019 demonstrated how to turn this algorithm into a post-processing strategy to achieve multi-accuracy, demonstrating empirical effectiveness across various domains.
This package provides a stable implementation of the multi-calibration algorithm, called MCBoost.
In contrast to other Fair ML approaches, MCBoost does not harm the overall utility of a prediction model, but rather aims at improving calibration and accuracy for large sets of subpopulations post-training.
MCBoost comes with strong theoretical guarantees, which have been explored formally in @hebert-johnson2018, @kim2019, @dwork-rankings, @dwork-oi and @kimkern2021.

`mcboost` implements Multi-Calibration Boosting for R.
`mcboost` is model agnostic and allows the user to post-process any supervised machine learning model.
It accepts initial models that fit binary outcomes or continuous outcomes with predictions that are in (or scaled to) the range [0, 1].
For convenience and ease of use, `mcboost` tightly integrates with the **mlr3** [@mlr3] machine learning eco-system in R by allowing to calibrate regression or classification models fitted either within or outside of mlr3.
Post-processing with `mcboost` starts with an initial prediction model that is passed on to an auditing algorithm that runs Multi-Calibration-Boosting on a labeled auditing dataset (Fig. 1). The resulting model can be used for obtaining multi-calibrated predictions.
`mcboost` includes two pre-defined learners for auditing (ridge regression and decision trees), and allows to easily adjust the learner and its parameters for Multi-Calibration Boosting.
Users may also specify a fixed set of subgroups, instead of a learner, on which predictions should be audited.
Furthermore, `mcboost` includes utilities to guard against overfitting to the auditing dataset during post-processing.

![Fig 1. Conceptual illustration of Multi-Calibration Boosting with `mcboost`.\label{fig:overview}](MCBoost.png)

# Statement of need

Given the ubiquitous use of machine learning models in crucial areas and growing concerns of biased predictions for minority subpopulations, Multi-Calibration Boosting should be widely accessible in the form of a free and open-source software package.
Prior to the development of `mcboost`, Multi-Calibration Boosting has not been released as a software package for R.

The results in @kim2019 highlight that MCBoost can improve classification accuracy for subpopulations in various settings, including gender detection with image data, income classification with survey data and disease prediction using biomedical data.
@Barda2020bias show that post-processing for Multi-Calibration can greatly improve calibration metrics of two medical risk assessment models when evaluated in subpopulations defined by intersections of age, sex, ethnicity, socioeconomic status and immigration history.
@Barda2020covid demonstrate that Multi-Calibration can also be used to adjust an initial classifier for a new task. They re-calibrate a baseline model for predicting the risk of severe respiratory infection with data on COVID-19 fatality rates in subpopulations, resulting in an accurate and calibrated COVID-19 mortality prediction model.


We hope that `mcboost` lets Multi-Calibration Boosting be utilized by a wide community of developers and data scientists to audit and post-process prediction models, and helps to promote fairness in machine learning and statistical estimation applications.

# Acknowledgements

We thank Matthew Sun for developing an initial Python implementation of MCBoost.
This work has been partially supported by the German Federal Ministry of Education and Research (BMBF) under Grant No. 01IS18036A. The authors of this work take full responsibilities for its content.

# References
