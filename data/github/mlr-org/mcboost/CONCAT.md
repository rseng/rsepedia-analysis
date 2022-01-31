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
---
title: "MCBoostSurv - Basics"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MCBoostSurv - Basics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```


```{r setup}
library("mcboost")
library("mlr3")
library("mlr3proba")
library("mlr3pipelines")
library("mlr3learners")
library("tidyverse")
set.seed(27099)
```


## Minimal Example: McBoostSurv

To show the basic functionality of `MCBoostSurv`, we provide a minimal example on 
the standard survival data set rats. After loading and pre-processing the data, we train
a `mlr3learner` on the training data. We instantiate a `MCBoostSurv` instance
with the default parameters. Then, we run the `$multicalibrate()` method on our data to start multi-calibration in survival analysis. With `$predict_probs()`, we can get
multicalibrated predictions. 

```{r}

#prepare task 
task = tsk("rats")
prep_pipe = po("encode", param_vals = list(method="one-hot")) 
prep = prep_pipe$train(list(task))[[1]]

#split data
train = prep$clone()$filter(1:199)
val = prep$clone()$filter(200:250)
test = prep$clone()$filter(256:300)

# get trained survival model 
baseline = lrn("surv.ranger")$train(train)

# initialize mcboost
mc_surv = MCBoostSurv$new(init_predictor = baseline)

# multicalibrate model 
mc_surv$multicalibrate(data = val$data(cols = val$feature_names), 
                       labels = val$data(cols = val$target_names))

# get new predictions
mc_surv$predict_probs(test$data(cols = test$feature_names))

```
## What does mcboost do?

Internally mcboostsurv runs the following procedure `max_iter` times (similar ro `mcboost`, just for distributions over time):

1. Predict on X using the model from the previous iteration, `init_predictor` in the first iteration.
1. Compute the residuals `res = y - y_hat` for all time points
1. Split predictions into `num_buckets` according to `y_hat` and time.
1. Fit the auditor (`auditor_fitter`) (here called`c(x)`) on the data in each bucket with target variable `r`.
1. Compute `misscal = mean(c(x) * res(x))`
1. if `misscal > alpha`:
    For the bucket with highest `misscal`, update the model using the prediction `c(x)`.
    else:
    Stop the procedure



## Multicalibrate model trained on PBC data

Based on this, we can now show multicalibration on a data set with two sensitive attributes (age and gender). Again, we load and pre-process the data.

### Load Dataset
```{r}
library(survival)
data_pbc = pbc %>%
    mutate(status = if_else(status == 2, 1, 0)
    ) %>%
    select(-id) %>%
    drop_na()

task_pbc = TaskSurv$new("pbc", backend = as_data_backend(data_pbc), 
                        time = "time", event = "status")


#Create data split

train_test = rsmp("holdout", ratio = 0.8)$instantiate(task_pbc)
train_g = train_test$train_set(1)
test_ids = train_test$test_set(1)
train_val = rsmp("holdout", ratio = 0.75)$instantiate(task_pbc$clone()$filter(train_g))
train_ids = train_val$train_set(1)
val_ids = train_val$test_set(1)

# Train distributional survival model 

xgb_distr = as_learner(ppl("distrcompositor",
               learner = as_learner(prep_pipe %>>% lrn("surv.xgboost"))))

xgb_distr$train(task_pbc$clone()$filter(train_ids))



```

### Mutlicalibrate survival model with validation data

```{r}
# initialize mcboost
 mcboost_learner = as_learner(
   prep_pipe %>>% ppl_mcboostsurv(
     learner = as_learner(prep_pipe %>>% xgb_distr), 
     param_vals = list(
       alpha = 1e-6,
       eta = 0.2,
       time_buckets = 2,
       num_buckets = 1 )
  )
)

# multicalibrate model
mcboost_learner$train(task_pbc$clone()$filter(val_ids))

# get new predictions
test_task = task_pbc$clone()$filter(test_ids)
pred_pbc_mc = mcboost_learner$predict(task_pbc$clone()$filter(test_ids))

pred_pbc_xgb = xgb_distr$predict(task_pbc$clone()$filter(test_ids))

```

### Development of IBS in the defined subgroups
```{r}

pred_pbc_xgb$score(msr("surv.graf"))
pred_pbc_mc$score(msr("surv.graf"))
```



---
title: "MCBoost - Basics and Extensions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MCBoost - Basics and Extensions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```


```{r setup}
library("mcboost")
library("mlr3")
set.seed(83007)
```


## Example 0: Multi-Accuracy in 6 lines of code

As a brief introduction we show how to use **mcboost** in only 6 lines of code.
For our example, we use the data from the *sonar* binary classification task.
We instantiate a `MCBoost` instance by specifying a `auditor_fitter`.
This `auditor_fitter` defines the splits into groups in each boosting iteration
based on the obtained residuals.
In this example, we choose a `Tree` based model.
Afterwards, we run the `$multicalibrate()` method on our data to start multi-calibration.
We only use the first 200 samples of the *sonar* data set to train our multi-calibrated model.


```{r}
tsk = tsk("sonar")
d = tsk$data(cols = tsk$feature_names)
l = tsk$data(cols = tsk$target_names)[[1]]
mc = MCBoost$new(auditor_fitter = "TreeAuditorFitter")
mc$multicalibrate(d[1:200,], l[1:200])
```

After the calibration, we use the model to predict on the left-out data (8 observations).

```{r}
mc$predict_probs(d[201:208,])
```


## What does mcboost do?

Internally mcboost runs the following procedure `max_iter` times:

1. Predict on X using the model from the previous iteration, `init_predictor` in the first iteration.
1. Compute the residuals `res = y - y_hat`
1. Split predictions into `num_buckets` according to `y_hat`.
1. Fit the auditor (`auditor_fitter`) (here called`c(x)`) on the data in each bucket with target variable `r`.
1. Compute `misscal = mean(c(x) * res(x))`
1. if `misscal > alpha`:
    For the bucket with highest `misscal`, update the model using the prediction `c(x)`.
    else:
    Stop the procedure

A lot more details can be found either in the code, or in the corresponding publications.


## Example 1: Multi-Accuracy Boosting on the Adult Dataset

First we download the data and create an `mlr3` classification task:

```{r}
library(data.table)
adult_train = fread(
  "https://raw.githubusercontent.com/Yorko/mlcourse.ai/master/data/adult_train.csv",
  stringsAsFactors = TRUE
)
adult_train$Country = NULL
adult_train$fnlwgt = NULL
train_tsk = TaskClassif$new("adult_train", adult_train, target = "Target")
```

We removed the features `Country` and `fnlwgt` since we expect them to have no predictive power.
`fnlwgt` means final weight and aims to allocate similar weights to people with similar demographic characteristics,
while `Country` has 42 distinct levels but 89 \% of the observations are from the United States.

### 1.1 Preprocessing

Then we do basic preprocessing:

  * Collapse rarest factors according to their prevalence
  * Drop missing factor levels
  * One-hot encode categorical variables
  * Impute NA's using a histogram approach

```{r}
library(mlr3pipelines)
pipe = po("collapsefactors", no_collapse_above_prevalence = 0.0006) %>>%
  po("fixfactors") %>>%
  po("encode") %>>%
  po("imputehist")
prep_task = pipe$train(train_tsk)[[1]]
```

In order to simulate settings where a sensitive feature is not available,
we remove the (dummy encoded) feature `Race` from the training task.

```{r}
prep_task$set_col_roles(c("Race.Amer.Indian.Eskimo", "Race.Asian.Pac.Islander", "Race.Black", "Race.Other", "Race.White"), remove_from = "feature")
```

Now we fit a `random forest`.

```{r}
library(mlr3learners)
l = lrn("classif.ranger", num.trees = 10L, predict_type = "prob")
l$train(prep_task)
```

### 1.2 MCBoost

A simple way to use the predictions from any `model` in **mcboost** is to wrap the predict
function and provide it as an initial predictor. This can be done from any model / any library.
Note, that we have to make sure, that our `init_predictor` returns a numeric vector of predictions.

```{r}
init_predictor = function(data) {
  l$predict_newdata(data)$prob[, 2]
}
```

As **mcboost** requires the data to be provided in `X, y` format (a `data.table` or `data.frame` of features and a
vector of labels), we create those two objects.

```{r}
data = prep_task$data(cols = prep_task$feature_names)
labels = 1 - one_hot(prep_task$data(cols = prep_task$target_names)[[1]])
```

We use a ridge regularized linear regression model as the auditor.

```{r}
mc = MCBoost$new(auditor_fitter = "RidgeAuditorFitter", init_predictor = init_predictor)
mc$multicalibrate(data, labels)
```

The `print` method additionally lists the average auditor values in the different buckets in each iteration:

```{r}
mc
```

### 1.3 Evaluation on Test Data

```{r}
adult_test = fread(
  "https://raw.githubusercontent.com/Yorko/mlcourse.ai/master/data/adult_test.csv",
  stringsAsFactors = TRUE
)
adult_test$Country = NULL
adult_test$fnlwgt = NULL

# The first row seems to have an error
adult_test = adult_test[Target != "",]
adult_test$Target = droplevels(adult_test$Target)

# Note, that we have to convert columns from numeric to integer here:
sdc = train_tsk$feature_types[type == "integer", id]
adult_test[, (sdc) := lapply(.SD, as.integer), .SDcols = sdc]

test_tsk = TaskClassif$new("adult_test", adult_test, target = "Target")
prep_test = pipe$predict(test_tsk)[[1]]
```

Now, we can again extract `X, y`.

```{r}
test_data = prep_test$data(cols = prep_test$feature_names)
test_labels = 1 - one_hot(prep_test$data(cols = prep_test$target_names)[[1]])
```

and **predict**.

```{r}
prs = mc$predict_probs(test_data)
```

The accuracy of the multi-calibrated model

```{r}
mean(round(prs) == test_labels)
```

is similar to the non-calibrated model.

```{r}
mean(round(init_predictor(test_data)) == test_labels)
```

But if we have a look at the bias for the different subpopulations of feature `Race`,
we can see that the predictions got more calibrated.
Note that we did not explicitly give neither the initial model
nor the auditor access to the feature `Race`.

```{r}
# Get bias per subgroup for multi-calibrated predictor
adult_test$biasmc = (prs - test_labels)
adult_test[, .(abs(mean(biasmc)), .N), by = .(Race)]
# Get bias per subgroup for initial predictor
adult_test$biasinit = (init_predictor(test_data) - test_labels)
adult_test[, .(abs(mean(biasinit)), .N), by = .(Race)]
```

### 1.4 The Auditor Effect

We can also obtain the auditor effect after multicalibration.
This indicates "how much" each observation has been affected by multi-calibration (on average across iterations).

```{r}
ae = mc$auditor_effect(test_data)
hist(ae)
```

We can see that there are a few instances with more pronounced effects, while most have actually only a low effect.

In order to get more insights, we compute quantiles of the
less and more effected population (median as cut-point) and analyze differences.

```{r}
effect = apply(test_data[ae >= median(ae[ae > 0]),], 2, quantile)
no_effect  = apply(test_data[ae < median(ae[ae>0]),], 2, quantile)
difference = apply((effect-no_effect), 2, mean)
difference[difference > 0.1]
```

There seems to be a difference in some variables like `Education` and `Marital_Status`.

We can further analyze the individuals:

```{r}
test_data[ae >= median(ae[ae>0]), names(which(difference > 0.1)), with = FALSE]
```

### Predicting using only the first 'n' iterations

Multi-calibration is an iterative procedure.
The `t` parameter can be used to predict using only the first `t` iterations.
This then predicts using only the first `t` iterations of the multi-calibration procedure.

```{r}
prs = mc$predict_probs(test_data, t = 3L)
```


## Example 2: MCBoost with non-mlr3 models: GLM

`mcboost` does not require your model to be a `mlr3` model.
As an input, `mcboost` expects a function `init_predictor` that takes as input `data` and returns a prediction.


```{r}
tsk = tsk("sonar")
data = tsk$data()[, Class := as.integer(Class) - 1L]
mod = glm(data = data, formula = Class ~ .)
```

The `init_predictor` could then use the `glm` model:

```{r}
init_predictor = function(data) {
  predict(mod, data)
}
```

... and we can calibrate this predictor.

```{r}
d = data[, -1]
l = data$Class
mc = MCBoost$new(init_predictor = init_predictor)
mc$multicalibrate(d[1:200,], l[1:200])
mc$predict_probs(d[201:208,])
```


## Example 3: Avoiding Overfitting in MCBoost

Very often `MCBoost`'s calibration is very aggressive and tends to overfit.
This section tries to introduce a method to regularize against this overfitting.

### 3.1 CVLearner

In this section we use a
`Cross-Validated` learner that predicts on held-out data during the training phase. This idea is based on Wolpert (1992)'s Stacked Generalization.
Other, simpler methods include choosing a smaller step size `eta` or reducing the number of `iters`.

```{r}
tsk = tsk("sonar")
```

As an `init_predictor` we again use a `ranger` model from mlr3 and
construct an init predictor using the convenience function provided by `mcboost`.

```{r}
learner = lrn("classif.ranger", predict_type = "prob")
learner$train(tsk)
init_predictor = mlr3_init_predictor(learner)
```

... and we can calibrate this predictor.
This time, we use a `CVTreeAuditorFitter` instead of a `TreeAuditorFitter`. This allows us to avoid
overfitting similar to a technique coined `stacked generalization` first described by Wolpert in 1992.
Note, that this can sometimes take a little longer since each learner is cross-validated using `3` folds (default).

```{r}
d = data[, -1]
l = data$Class
mc = MCBoost$new(init_predictor = init_predictor, auditor_fitter=CVTreeAuditorFitter$new(), max_iter = 2L)
mc$multicalibrate(d[1:200,], l[1:200])
mc$predict_probs(d[201:208,])
```

### 3.2 Data Splitting

We can also use a fresh chunk of the validation data in each iteration. `mcboost` implements two strategies, `"bootstrap"` and `"split"`. While `"split"` simply splits up the data,  `"bootstrap"` draws a new bootstrap sample of the data in each iteration.

```{r}
tsk = tsk("sonar")
```

Again, we use a `ranger` mlr3 model as our initial predictor:

```{r}
learner = lrn("classif.ranger", predict_type = "prob")
learner$train(tsk)
init_predictor = mlr3_init_predictor(learner)
```

and we can now calibrate:

```{r}
d = data[, -1]
l = data$Class
mc = MCBoost$new(
  init_predictor = init_predictor,
  auditor_fitter= TreeAuditorFitter$new(),
  iter_sampling = "bootstrap"
)
mc$multicalibrate(d[1:200,], l[1:200])
mc$predict_probs(d[201:208,])
```


## Example 4: Adjusting the SubPop Fitter

For this example, we use the *sonar* dataset once again:

```{r}
tsk = tsk("sonar")
data = tsk$data(cols = tsk$feature_names)
labels = tsk$data(cols = tsk$target_names)[[1]]
```

### 4.1 LearnerAuditorFitter

The Subpop-fitter can be easily adjusted by constructing it from a `LearnerAuditorFitter`.
This allows for using any **mlr3** learner.
See [here](https://mlr3extralearners.mlr-org.com/articles/learners/list_learners.html) for a list of available learners.

```{r}
rf = LearnerAuditorFitter$new(lrn("regr.rpart", minsplit = 10L))
mc = MCBoost$new(auditor_fitter = rf)
mc$multicalibrate(data, labels)
```

The `TreeAuditorFitter` and `RidgeAuditorFitter` are two instantiations of this Fitter with pre-defined learners. By providing their character strings the fitter could be automatically constructed.

### 4.2 SubpopAuditorFitter & SubgroupAuditorFitter

In some occasions, instead of using a `Learner`, we might want to use a fixed set of subgroups.
Those can either be defined from the data itself or provided from the outside.

**Splitting via the dataset**

In order to split the data into groups according to a set of columns, we use a `SubpopAuditorFitter`
together with a list of `subpops`. Those define the group splits to multi-calibrate on.
These splits can be either a `character` string, referencing a binary variable in the data
or a `function` that, when evaluated on the data, returns a binary vector.

In order to showcase both options, we add a binary variable to our `data`:

```{r}
data[, Bin := sample(c(1, 0), nrow(data), replace = TRUE)]
```

```{r}
rf = SubpopAuditorFitter$new(list(
  "Bin",
  function(data) {data[["V1"]] > 0.2},
  function(data) {data[["V1"]] > 0.2 | data[["V3"]] < 0.29}
))
```

```{r}
mc = MCBoost$new(auditor_fitter = rf)
mc$multicalibrate(data, labels)
```

And we can again apply it to predict on new data:

```{r}
mc$predict_probs(data)
```

**Manually defined masks**

If we want to add the splitting from the outside, by supplying binary masks for the
rows of the data, we can provide manually defined masks.
Note, that the masks have to correspond with the number of rows in the dataset.

```{r}
rf = SubgroupAuditorFitter$new(list(
  rep(c(0, 1), 104),
  rep(c(1, 1, 1, 0), 52)
))
```

```{r}
mc = MCBoost$new(auditor_fitter = rf)
mc$multicalibrate(data, labels)
```

During prediction, we now have to supply a set of masks for the prediction data.

```{r}
predict_masks = list(
  rep(c(0, 1), 52),
  rep(c(1, 1, 1, 0), 26)
)
```

```{r}
mc$predict_probs(data[1:104,], subgroup_masks = predict_masks)
```


## Example 5: Multi-Calibrating data with missing values using a pipeline

When data has missing values or other non-standard columns, we often have to pre-process data in order
to be able to fit models.
Those preprocessing steps can be embedded into the `SubPopFitter` by using a **mlr3pipelines** Pipeline.
The following code shows a brief example:

```{r}
tsk = tsk("penguins")
# first we convert to a binary task
row_ids = tsk$data(cols = c("species", "..row_id"))[species %in% c("Adelie", "Gentoo")][["..row_id"]]
tsk$filter(row_ids)$droplevels()
tsk
```

```{r}
library("mlr3pipelines")
library("mlr3learners")

# Convert task to X,y
X = tsk$data(cols = tsk$feature_names)
y = tsk$data(cols = tsk$target_names)

# Our inital model is a pipeline that imputes missings and encodes categoricals
init_model = as_learner(po("encode") %>>% po("imputehist") %>>%
  lrn("classif.glmnet", predict_type = "prob"))
# And we fit it on a subset of the data in order to simulate a poorly performing model.
init_model$train(tsk$clone()$filter(row_ids[c(1:9, 160:170)]))
init_model$predict(tsk)$score()

# We define a pipeline that imputes missings and encodes categoricals
auditor = as_learner(po("encode") %>>% po("imputehist") %>>% lrn("regr.rpart"))

mc = MCBoost$new(auditor_fitter = auditor, init_predictor = init_model)
mc$multicalibrate(X, y)
```

and we can observe where it improved:

```{r}
mc
```


## Example 6: Multi-Calibration Regression

We abuse the `Communities & Crime` dataset in order to showcase how `mcboost` can be used in a regression setting.

First we download the data and create an `mlr3` regression task:

```{r}
library(data.table)
library(mlr3oml)
oml = OMLData$new(42730)
data = oml$data

tsk = TaskRegr$new("communities_crime", data, target = "ViolentCrimesPerPop")
```

Currently, **mcboost** only allows to work with targets between 0 and 1.
Luckily, our target variable's values are already in that range, but
if they were not, we could simply scale them to [0;1] before our analysis.

```{r}
summary(data$ViolentCrimesPerPop)
```

We again split our task into **train** and **test**.
We do this in `mlr3` by simply setting some (here 500) row roles to `"validation"`.

```{r}
tsk$set_row_roles(sample(tsk$row_roles$use, 500), "validation")
```

### 6.1 Preprocessing

Then we do basic preprocessing, since we do not have any categorical
variables, we only impute NA's using a histogram approach.

```{r}
library(mlr3pipelines)
pipe =  po("imputehist")
prep_task = pipe$train(list(tsk))[[1]]

prep_task$set_col_roles(c("racepctblack", "racePctWhite", "racePctAsian", "racePctHisp", "community"), remove_from = "feature")
```

Now we fit our first `Learner`: A `random forest`.

```{r}
library(mlr3learners)
l = lrn("regr.ranger", num.trees = 10L)
l$train(prep_task)
```

### 6.2 MCBoost

A simple way to use the predictions from any `Model` in **mcboost** is to wrap the predict
function and provide it as an initial predictor. This can be done from any model / any library.
Note, that we have to make sure, that our `init_predictor` returns a numeric vector of predictions.

```{r}
init_predictor = function(data) {
  l$predict_newdata(data)$response
}
```

As **mcboost** requires the data to be provided in `X, y` format (a `data.table` or `data.frame` of features and a
vector of labels), we create those two objects.

```{r}
data = prep_task$data(cols = prep_task$feature_names)
labels = prep_task$data(cols = prep_task$target_names)[[1]]
```

```{r}
mc = MCBoost$new(auditor_fitter = "RidgeAuditorFitter", init_predictor = init_predictor, eta = 0.1)
mc$multicalibrate(data, labels)
```

### 6.3 Evaluation on Test Data

We first create the test task by setting the `validation` rows to `use`, and then
use our preprocessing `pipe's`  predict function to also impute missing values
for the validation data. Then we again extract features `X` and target `y`.

```{r}
test_task = tsk$clone()
test_task$row_roles$use = test_task$row_roles$validation
test_task = pipe$predict(list(test_task))[[1]]
test_data = test_task$data(cols = tsk$feature_names)
test_labels = test_task$data(cols = tsk$target_names)[[1]]
```

and **predict**.

```{r}
prs = mc$predict_probs(test_data)
```

Now we can compute the MSE of the multi-calibrated model

```{r}
mean((prs - test_labels)^2)
```

and compare to the non-calibrated version:

```{r}
mean((init_predictor(test_data) - test_labels)^2)
```

But looking at sub-populations we can see that the predictions got
more calibrated.
Since we cannot show all subpopulations we only show the MSE for the feature `racepctblack`.

```{r}
test_data$se_mcboost = (prs - test_labels)^2
test_data$se_init = (init_predictor(test_data) - test_labels)^2

test_data[, .(mcboost = mean(se_mcboost), initial = mean(se_init), .N), by = .(racepctblack > 0.5)]
```
---
title: "MCBoost - Health Survey Example"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MCBoost - Health Survey Example}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{}
---

```{r, echo = FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

```{r setup, message = FALSE}
library(tidyverse)
library(PracTools)
library(ranger)
library(neuralnet)
library(formattable)
library(mlr3)
library(mlr3learners)
library(mcboost)
```

## Data and Setup

This vignette presents two typical use cases of MCBoost with data from a health survey. The goal is to post-process two initial prediction models for multi-accuracy using different flavors of MCBoost, and to eventually compare the naive and post-processed predictors overall and for subpopulations. The first scenario starts with a neural net and, as an example, evaluates the initial and post-processed predictors with a focus on subgroup accuracy after running MCBoost. The second scenario uses a random forest and evaluates the initial and post-processed predictors with respect to subgroup calibration.

We use data derived from the National Health Interview Survey (NHIS 2003), which includes demographic and health-related variables for 21,588 individuals. This data can directly be included from the `PracTools` package.

```{r, eval = TRUE}
data(nhis.large)
```
We can obtain more information using: 

```{r, eval = FALSE}
?nhis.large
```


In the following, our outcome of interest is whether an individual is covered by any type of health insurance (`notcov`, 1 = not covered, 0 = covered). We additionally prepare two sets of variables:

- Predictor variables (age, parents in household, education, income, employment status, physical or other limitations)
- Subpopulation variables (sex, hispanic ethnicity, race)

The second set of variables will not be used for training the initial prediction models, but will be our focus when it comes to evaluating prediction performance for subgroups.

Before we training an initial model, we preprocess the data:

 * We encode categorical features as `factor`.
 * We explicitly assign `NA`s in categorical features to a dedicated factor level
 * We drop `NA`s in the outcome variable `notcov`
 * We encode `notcov` as a factor variable instead of a dummy variable (1 = `notcov`, 0 = `cov`)
 * We create a new feature `inv_wt` as the inverse of survey weights `svwyt`

```{r}
categorical <- c("age.grp", "parents", "educ", "inc.grp", "doing.lw",
  "limited", "sex", "hisp", "race")

nhis <- nhis.large %>%
  mutate_at(categorical, as.factor) %>%
  mutate_at(categorical, fct_explicit_na) %>%
  drop_na(notcov) %>%
  select(all_of(categorical), notcov, svywt, ID)

nhis$notcov <- factor(ifelse(nhis$notcov == 1, "notcov", "cov"))

nhis_enc <- data.frame(model.matrix(notcov ~ ., data = nhis)[,-1])
nhis_enc$notcov <- nhis$notcov
nhis_enc$sex <- nhis$sex
nhis_enc$hisp <- nhis$hisp
nhis_enc$race <- nhis$race
nhis_enc$inv_wt <- (1 / nhis$svywt)
```

The pre-processed NHIS data will be split into three datasets:

- A training set `train` for training the initial prediction models (55 \% of data)
- An auditing set `post` for post-processing the initial models with MCBoost (20 \%)
- A test set `test`for model evaluation (25 \%)

To increase the difficulty of the prediction task, we sample from the NHIS data such that the prevalence of demographic subgroups in the test data differs from their prevalence in the training and auditing data. This is achieved by employing weighted sampling from NHIS (variable `inv_wt` from above).

```{r}
set.seed(2953)

test <- nhis_enc %>% slice_sample(prop = 0.25, weight_by = inv_wt)

nontest_g <- nhis_enc %>% anti_join(test, by = "ID")

train_g <- nontest_g %>% slice_sample(prop = 0.75)

post <- nontest_g %>% anti_join(train_g, by = "ID") %>% select(-ID, -svywt, -inv_wt, -c(sex:race))

train <- train_g %>% select(-ID, -svywt, -inv_wt, -c(sex:race), -c(sex2:race3))
```

As a result, non-hispanic white individuals (`hisp2`) are overrepresented and hispanic individuals are underrepresented in both the training and auditing set, compared to their prevalence in the test set.

```{r}
train_g %>% summarise_at(vars(sex2:race3), mean)
# hispanic individuals
1 - sum(train_g %>% summarise_at(vars(hisp2:hisp4), mean))
```

```{r}
post %>% summarise_at(vars(sex2:race3), mean)
# hispanic individuals
1 - sum(post %>% summarise_at(vars(hisp2:hisp4), mean))
```

```{r}
test %>% summarise_at(vars(sex2:race3), mean)
# hispanic individuals
1 - sum(test %>% summarise_at(vars(hisp2:hisp4), mean))
```

## Scenario 1: Improve Subgroup Accuracy

We train an initial model for predicting healthcare coverage with the training set. Here, we use a neural network with one hidden layer, rather naively with little tweaking.

```{r, message = FALSE}
nnet <- neuralnet(notcov ~ .,
  hidden = 5,
  linear.output = FALSE,
  err.fct = 'ce',
  threshold = 0.5,
  lifesign = 'full',
  data = train
)
```

### MCBoost Auditing

We prepare a function that allows us to pass the predictions of the model to MCBoost for post-processing.

```{r}
init_nnet = function(data) {
  predict(nnet, data)[, 2]
}
```

To showcase different use cases of MCBoost, we prepare two post-processing data sets based on the auditing set. The first set includes only the predictor variables that were used by the initial models, whereas the second set will allow post-processing based on our demographic subgroups of interest (sex, hispanic ethnicity, race).

```{r}
d1 <- select(post, -c(notcov, sex2:race3))
d2 <- select(post, -notcov)
l <- 1 - one_hot(post$notcov)
```

We initialize two custom auditors for MCBoost: Ridge regression with a small penalty on model complexity, and a `SubpopAuditorFitter` with a fixed set of subpopulations.

```{r}
ridge = LearnerAuditorFitter$new(lrn("regr.glmnet", alpha = 0, lambda = 2 / nrow(post)))

pops = SubpopAuditorFitter$new(list("sex2", "hisp2", "hisp3", "hisp4", "race2", "race3"))
```

The ridge regression will only be given access to the initial predictor variables when post-processing the neural net predictions with the auditing data. In contrast, we guide the subpop-fitter to audit the initial predictions explicitly on the outlined subpopulations (sex, hispanic ethnicity, race). In summary, we have:

- `nnet`: Initial neural net
- `nnet_mc_ridge`: Neural net, post-processed with ridge regression and the initial set of predictor variables
- `nnet_mc_subpop`: Neural net, post-processed with a fixed set of subpopulations

```{r}
nnet_mc_ridge = MCBoost$new(init_predictor = init_nnet,
                            auditor_fitter = ridge,
                            multiplicative = TRUE,
                            partition = TRUE,
                            max_iter = 15)
nnet_mc_ridge$multicalibrate(d1, l)

nnet_mc_subpop = MCBoost$new(init_predictor = init_nnet,
                             auditor_fitter = pops,
                             partition = TRUE,
                             max_iter = 15)
nnet_mc_subpop$multicalibrate(d2, l)
```

### Model Evaluation

Next, we use the initial and post-processed models to predict the outcome in the test data. We compute predicted probabilities and class predictions.

```{r}
test$nnet <- predict(nnet, newdata = test)[, 2]
test$nnet_mc_ridge <- nnet_mc_ridge$predict_probs(test)
test$nnet_mc_subpop <- nnet_mc_subpop$predict_probs(test)

test$c_nnet <- round(test$nnet)
test$c_nnet_mc_ridge <- round(test$nnet_mc_ridge)
test$c_nnet_mc_subpop <- round(test$nnet_mc_subpop)
test$label <- 1 - one_hot(test$notcov)
```

Here we compare the overall accuracy of the initial and post-processed models. Overall, we observe little differences in performance.

```{r}
mean(test$c_nnet == test$label)
mean(test$c_nnet_mc_ridge == test$label)
mean(test$c_nnet_mc_subpop == test$label)
```

However, we might be concerned with model performance for smaller subpopulations. In the following, we focus on subgroups defined by 2-way conjunctions of sex, hispanic ethnicity, and race.

```{r, warning = FALSE}
test <- test %>%
  group_by(sex, hisp) %>%
  mutate(sex_hisp = cur_group_id()) %>%
  group_by(sex, race) %>%
  mutate(sex_race = cur_group_id()) %>%
  group_by(hisp, race) %>%
  mutate(hisp_race = cur_group_id()) %>%
  ungroup()

grouping_vars <- c("sex", "hisp", "race", "sex_hisp", "sex_race", "hisp_race")

eval <- map(grouping_vars, group_by_at, .tbl = test) %>%
  map(summarise,
      'accuracy_nnet' = mean(c_nnet == label),
      'accuracy_nnet_mc_ridge' = mean(c_nnet_mc_ridge == label),
      'accuracy_nnet_mc_subpop' = mean(c_nnet_mc_subpop == label),
      'size' = n()) %>%
  bind_rows()
```

We evaluate classification accuracy on these subpopulations, and order the results according to the size of the selected subgroups (`size`). Subgroup accuracy varies between methods, with MCBoost-Ridge (`nnet_mc_ridge`) and MCBoost-Subpop (`nnet_mc_subpop`) stabilizing subgroup performance when compared to the initial model, respectively.

```{r}
eval %>%
  arrange(desc(size)) %>%
  select(size, accuracy_nnet:accuracy_nnet_mc_subpop) %>%
  round(., digits = 3) %>%
  formattable(., lapply(1:nrow(eval), function(row) {
  area(row, col = 2:4) ~ color_tile("transparent", "lightgreen")
    }))
```

## Scenario 2: Improve Subgroup Calibration

In this scenario, we use a random forest with the default settings of the ranger package as the initial predictor.

```{r}
rf <- ranger(notcov ~ ., data = train, probability = TRUE)
```

### MCBoost Auditing

We again prepare a function to pass the predictions to MCBoost for post-processing.

```{r}
init_rf = function(data) {
  predict(rf, data)$prediction[, 2]
}
```

We use two custom auditors for MCBoost, i.e., ridge and lasso regression with different penalties on model complexity.

```{r}
ridge = LearnerAuditorFitter$new(lrn("regr.glmnet", alpha = 0, lambda = 2 / nrow(post)))

lasso = LearnerAuditorFitter$new(lrn("regr.glmnet", alpha = 1, lambda = 40 / nrow(post)))
```

The ridge regression will only be given access to the initial predictor variables when post-processing the random forest predictions. In contrast, we allow the lasso regression to audit the initial predictions both with the initial predictors and the subpopulations (sex, hispanic ethnicity, race). In summary, we have:

- `rf`: Initial random forest
- `rf_mc_ridge`: Random forest, post-processed with ridge regression and the initial set of predictor variables
- `rf_mc_lasso`: Random forest, post-processed with lasso regression and the extended set of predictors

```{r}
rf_mc_ridge = MCBoost$new(init_predictor = init_rf,
                          auditor_fitter = ridge,
                          multiplicative = TRUE,
                          partition = TRUE,
                          max_iter = 15)
rf_mc_ridge$multicalibrate(d1, l)

rf_mc_lasso = MCBoost$new(init_predictor = init_rf,
                          auditor_fitter = lasso,
                          multiplicative = TRUE,
                          partition = TRUE,
                          max_iter = 15)
rf_mc_lasso$multicalibrate(d2, l)
```

### Model Evaluation

We again compute predicted probabilities and class predictions using the initial and post-processed models.

```{r}
test$rf <- predict(rf, test)$prediction[, 2]
test$rf_mc_ridge <- rf_mc_ridge$predict_probs(test)
test$rf_mc_lasso <- rf_mc_lasso$predict_probs(test)

test$c_rf <- round(test$rf)
test$c_rf_mc_ridge <- round(test$rf_mc_ridge)
test$c_rf_mc_lasso <- round(test$rf_mc_lasso)
```

Here we compare the overall accuracy of the initial and post-processed models. As before, we observe small differences in overall performance.

```{r}
mean(test$c_rf == test$label)
mean(test$c_rf_mc_ridge == test$label)
mean(test$c_rf_mc_lasso == test$label)
```

However, we might be concerned with calibration in subpopulations. In the following we focus on subgroups defined by 2-way conjunctions of sex, hispanic ethnicity, and race.

```{r}
eval <- map(grouping_vars, group_by_at, .tbl = test) %>%
  map(summarise,
      'bias_rf' = abs(mean(rf) - mean(label))*100,
      'bias_rf_mc_ridge' = abs(mean(rf_mc_ridge) - mean(label))*100,
      'bias_rf_mc_lasso' = abs(mean(rf_mc_lasso) - mean(label))*100,
      'size' = n()) %>%
  bind_rows()
```

This evaluation focuses on the difference between the average predicted risk of healthcare non-coverage and the observed proportion of non-coverage in the test data for subgroups. Considering the MCBoost-Ridge (`rf_mc_ridge`) and MCBoost-Lasso (`rf_mc_lasso`) results, post-processing with MCBoost reduces bias for many subpopulations.

```{r}
eval %>%
  arrange(desc(size)) %>%
  select(size, bias_rf:bias_rf_mc_lasso) %>%
  round(., digits = 3) %>%
  formattable(., lapply(1:nrow(eval), function(row) {
  area(row, col = 2:4) ~ color_tile("lightgreen", "transparent")
    }))
```
---
title: "MCBoost step by step"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MCBoost step by step}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(mcboost)
```

```{r}
library(data.table)
n = 10000L
x = rnorm(n, 5, 2.5)
s = sample(1:2, n, replace = TRUE)
itcpt = c(5.2, 1.1)
betas = c(-.12, .7)
y = x * betas[s] + itcpt[s] + rnorm(n, -3, .4) + 1
dt = data.table(x = x, s = s, y = y)
dt[, yprob := 1 / (1 + exp(-(y - mean(y)))) + rnorm(n, 0, abs(0.1*(s-2)))]
dt[, train:=FALSE][1:(ceiling(n/2)), train := TRUE]
dt[, y := as.integer(runif(n) > yprob)]
```

```{r}
library(ggplot2)
ggplot(dt) + geom_point(aes(x=x,y=yprob,color=factor(s)))
dt[, mean(y), by = s]
```

```{r}
mod = glm(y ~ x, data = dt[train == TRUE,], family = binomial())
dt[, yh := predict(mod, dt, type = "response")]
dt[, .(mean((yh > 0.5) == y), .N), by = .(s, train)]
```

```{r}
ggplot(dt) + geom_point(aes(x=yprob,y=yh,color=factor(s)))
dt[, mean(y), by = s]
```


```{r}
init_predictor = function(data) {
  predict(mod, data, type = "response")
}
mc = MCBoost$new(
  auditor_fitter = "TreeAuditorFitter",
  init_predictor = init_predictor,
  max_iter = 10L,
  eta = .2,
  multiplicative = TRUE
)
mc$multicalibrate(dt[, c("s", "x"), with = FALSE], dt$y)
dt[, yh_mc := mc$predict_probs(dt)]
ggplot(dt) + geom_point(aes(x=yprob,y=yh_mc,color=factor(s)))
dt[, mean(y), by = s]
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuditorFitters.R
\name{AuditorFitter}
\alias{AuditorFitter}
\title{AuditorFitter Abstract Base Class}
\value{
\code{list} with items\cr
\itemize{
\item \code{corr}: pseudo-correlation between residuals and learner prediction.
\item \code{l}: the trained learner.
}
}
\description{
Defines an \code{AuditorFitter} abstract base class.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{AuditorFitter$new()}}
\item \href{#method-fit_to_resid}{\code{AuditorFitter$fit_to_resid()}}
\item \href{#method-fit}{\code{AuditorFitter$fit()}}
\item \href{#method-clone}{\code{AuditorFitter$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a \code{\link{AuditorFitter}}.
This is an abstract base class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AuditorFitter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit_to_resid"></a>}}
\if{latex}{\out{\hypertarget{method-fit_to_resid}{}}}
\subsection{Method \code{fit_to_resid()}}{
Fit to residuals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AuditorFitter$fit_to_resid(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit (mostly used internally, use \code{fit_to_resid}).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AuditorFitter$fit(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCBoostSurv.R
\name{MCBoostSurv}
\alias{MCBoostSurv}
\title{Multi-Calibration Boosting}
\description{
Implements Multi-Calibration Boosting by Hebert-Johnson et al. (2018) and
Multi-Accuracy Boosting by Kim et al. (2019) for the multi-calibration of a
machine learning model's prediction for survival models.
Multi-Calibration works best in scenarios where the underlying data & labels are unbiased
but a bias is introduced within the algorithm's fitting procedure. This is often the case,
e.g. when an algorithm fits a majority population while ignoring or under-fitting minority
populations.\cr
Expects initial models that predict probobilities (between 0 and 1) for different time points.
The method defaults to \verb{Multi-Accuracy Boosting} as described in Kim et al. (2019).
In order to obtain behaviour as described in Hebert-Johnson et al. (2018) set
\code{multiplicative=FALSE} and \code{num_buckets} to 10.
\itemize{
For additional details, please refer to the relevant publications:
\item{Hebert-Johnson et al., 2018. Multicalibration: Calibration for the (Computationally-Identifiable) Masses.
Proceedings of the 35th International Conference on Machine Learning, PMLR 80:1939-1948.
https://proceedings.mlr.press/v80/hebert-johnson18a.html.}{}
\item{Kim et al., 2019. Multiaccuracy: Black-Box Post-Processing for Fairness in Classification.
Proceedings of the 2019 AAAI/ACM Conference on AI, Ethics, and Society (AIES '19).
Association for Computing Machinery, New York, NY, USA, 247–254.
https://dl.acm.org/doi/10.1145/3306618.3314287}{}
}
}
\section{Super class}{
\code{\link[mcboost:MCBoost]{mcboost::MCBoost}} -> \code{MCBoostSurv}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{max_iter}}{\code{\link{integer}} \cr
The maximum number of iterations of the multi-calibration/multi-accuracy method.}

\item{\code{alpha}}{\code{\link{numeric}} \cr
Accuracy parameter that determines the stopping condition.}

\item{\code{eta}}{\code{\link{numeric}} \cr
Parameter for multiplicative weight update (step size).}

\item{\code{num_buckets}}{\code{\link{integer}} \cr
The number of buckets to split into in addition to using the whole sample.}

\item{\code{bucket_strategy}}{\code{\link{character}} \cr
Currently only supports "simple", even split along probabilities.
Only relevant for \code{num_buckets} > 1.}

\item{\code{rebucket}}{\code{\link{logical}} \cr
Should buckets be re-calculated at each iteration?}

\item{\code{eval_fulldata}}{\code{\link{logical}} \cr
Should auditor be evaluated on the full data?}

\item{\code{partition}}{\code{\link{logical}} \cr
True/False flag for whether to split up predictions by their "partition"
(e.g., predictions less than 0.5 and predictions greater than 0.5).}

\item{\code{multiplicative}}{\code{\link{logical}} \cr
Specifies the strategy for updating the weights (multiplicative weight vs additive).}

\item{\code{iter_sampling}}{\code{\link{character}} \cr
Specifies the strategy to sample the validation data for each iteration.}

\item{\code{auditor_fitter}}{\code{\link{AuditorFitter}} \cr
Specifies the type of model used to fit the residuals.}

\item{\code{predictor}}{\code{\link{function}} \cr
Initial predictor function.}

\item{\code{iter_models}}{\code{\link{list}} \cr
Cumulative list of fitted models.}

\item{\code{iter_partitions}}{\code{\link{list}} \cr
Cumulative list of data partitions for models.}

\item{\code{iter_corr}}{\code{\link{list}} \cr
Auditor correlation in each iteration.}

\item{\code{auditor_effects}}{\code{\link{list}} \cr
Auditor effect in each iteration.}

\item{\code{time_points}}{\code{\link{integer}} \cr
Times included in the prediction (columnames)}

\item{\code{time_buckets}}{\code{\link{integer}} \cr
The number of buckets to split the time points (columns) of the prediction.}

\item{\code{bucket_strategies}}{\code{\link{character}} \cr
Possible bucket_strategies in McBoostSurv.
Only relevant for \code{time_buckets} > 1.
\code{even_splits}: split buckets evenly
\code{quantiles}: split buckets by quantiles}

\item{\code{bucket_aggregation}}{\code{\link{function}} \cr
If not NULL, predictions are not selected by time/probability,
but by time/individual. Individuals are selected by aggregated value per
individual (e.g. mean).
Only relevant for \code{time_buckets} > 1.}

\item{\code{max_time_quantile}}{\code{\link{double}} \cr
Time quantile which should be evaluated and multicalibrated.
Similar to a 75\%-Integrated Brier Score.}

\item{\code{time_points_eval}}{\code{\link{integer}} \cr
Vector of time_points that should be evaluated.}

\item{\code{loss}}{\code{\link{character}} \cr
Loss function which is optimized during boosting.
\code{censored_brier}: censored version of the integrated brier score
\code{brier}: uncensored version of the integrated brier score
\code{censored_brier_proper}: proper version of the censored version of the integrated brier score
For more details, we are referring to https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.graf.html.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{MCBoostSurv$new()}}
\item \href{#method-clone}{\code{MCBoostSurv$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="MCBoost" data-id="auditor_effect">}\href{../../mcboost/html/MCBoost.html#method-auditor_effect}{\code{mcboost::MCBoost$auditor_effect()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="MCBoost" data-id="multicalibrate">}\href{../../mcboost/html/MCBoost.html#method-multicalibrate}{\code{mcboost::MCBoost$multicalibrate()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="MCBoost" data-id="predict_probs">}\href{../../mcboost/html/MCBoost.html#method-predict_probs}{\code{mcboost::MCBoost$predict_probs()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="MCBoost" data-id="print">}\href{../../mcboost/html/MCBoost.html#method-print}{\code{mcboost::MCBoost$print()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a multi-calibration instance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoostSurv$new(
  max_iter = 25,
  alpha = 1e-04,
  eta = 0.1,
  num_buckets = 1,
  partition = ifelse(num_buckets > 1, TRUE, FALSE),
  time_buckets = 2L,
  max_time_quantile = 1,
  bucket_strategy = "even_splits",
  bucket_aggregation = NULL,
  rebucket = FALSE,
  eval_fulldata = FALSE,
  multiplicative = TRUE,
  auditor_fitter = "RidgeAuditorFitter",
  subpops = NULL,
  default_model_class = LearnerSurvKaplan,
  init_predictor = NULL,
  loss = "censored_brier",
  iter_sampling = "none"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{max_iter}}{\code{\link{integer}} \cr
The maximum number of iterations of the multi-calibration/multi-accuracy method.
Default \code{5L}.}

\item{\code{alpha}}{\code{\link{numeric}} \cr
Accuracy parameter that determines the stopping condition. Default \code{1e-4}.}

\item{\code{eta}}{\code{\link{numeric}} \cr
Parameter for multiplicative weight update (step size). Default \code{1.0}.}

\item{\code{num_buckets}}{\code{\link{integer}} \cr
The number of buckets to split into in addition to using the whole sample. Default \code{2L}.}

\item{\code{partition}}{\code{\link{logical}} \cr
True/False flag for whether to split up predictions by their "partition"
(e.g., predictions less than 0.5 and predictions greater than 0.5).
Defaults to \code{TRUE} (multi-accuracy boosting).}

\item{\code{time_buckets}}{\code{\link{integer}} \cr
The number of buckets to split the time points (columns) of the prediction.}

\item{\code{max_time_quantile}}{\code{\link{double}} \cr
Time quantile which should be evaluated and multicalibrated.
Can be used to perform multi-calibration only up to the \code{max_time_quantile} percent of timepoints.
Initialized to \code{1}.}

\item{\code{bucket_strategy}}{\code{\link{character}} \cr
Bucketstragy for bucketing.
\code{even_splits}: split buckets evenly
\code{quantiles}: split buckets by quantiles}

\item{\code{bucket_aggregation}}{\code{\link{function}} \cr
If not NULL, predictions are not selected by time/probability,
but by time/individual. Individuals are selected by aggregated value per
individual (e.g. mean).
Only relevant for \code{time_buckets} > 1.}

\item{\code{rebucket}}{\code{\link{logical}} \cr
Should buckets be re-done at each iteration? Default \code{FALSE}.}

\item{\code{eval_fulldata}}{\code{\link{logical}} \cr
Should the auditor be evaluated on the full data or on the respective bucket for determining
the stopping criterion? Default \code{FALSE}, auditor is only evaluated on the bucket.
This setting keeps the implementation closer to the Algorithm proposed in the corresponding
multi-accuracy paper (Kim et al., 2019) where auditor effects are computed across the full
sample (i.e. eval_fulldata = TRUE).}

\item{\code{multiplicative}}{\code{\link{logical}} \cr
Specifies the strategy for updating the weights (multiplicative weight vs additive).
Defaults to \code{TRUE} (multi-accuracy boosting). Set to \code{FALSE} for multi-calibration.}

\item{\code{auditor_fitter}}{\code{\link{AuditorFitter}}|\code{\link{character}}|\code{\link[mlr3:Learner]{mlr3::Learner}} \cr
Specifies the type of model used to fit the
residuals. The default is \code{\link{RidgeAuditorFitter}}.
Can be a \code{character}, the name of a \code{\link{AuditorFitter}}, a \code{\link[mlr3:Learner]{mlr3::Learner}} that is then
auto-converted into a \code{\link{LearnerAuditorFitter}} or a custom \code{\link{AuditorFitter}}.}

\item{\code{subpops}}{\code{\link{list}} \cr
Specifies a collection of characteristic attributes
and the values they take to define subpopulations
e.g. list(age = c('20-29','30-39','40+'), nJobs = c(0,1,2,'3+'), ,..).}

\item{\code{default_model_class}}{\code{Predictor} \cr
The class of the model that should be used as the init predictor model if
\code{init_predictor} is not specified. Defaults to \code{ConstantPredictor} which
predicts a constant value.}

\item{\code{init_predictor}}{\code{\link{function}}|\code{\link[mlr3:Learner]{mlr3::Learner}} \cr
The initial predictor function to use (i.e., if the user has a pretrained model).
If a \code{mlr3} \code{Learner} is passed, it will be autoconverted using \code{mlr3_init_predictor}.
This requires the \code{\link[mlr3:Learner]{mlr3::Learner}} to be trained.}

\item{\code{loss}}{\code{\link{character}} \cr
#' Loss function which is optimized during boosting.
\code{censored_brier}: censored version of the integrated brier score
\code{brier}: uncensored version of the integrated brier score
\code{censored_brier_proper}: proper version of the censored version of the integrated brier score
For more details, we are referring to https://mlr3proba.mlr-org.com/reference/mlr_measures_surv.graf.html.}

\item{\code{iter_sampling}}{\code{\link{character}} \cr
How to sample the validation data for each iteration?
Can be \code{bootstrap}, \code{split} or \code{none}.\cr
"split" splits the data into \code{max_iter} parts and validates on each sample in each iteration.\cr
"bootstrap" uses a new bootstrap sample in each iteration.\cr
"none" uses the same dataset in each iteration.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoostSurv$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PipelineMCBoost.R
\name{ppl_mcboostsurv}
\alias{ppl_mcboostsurv}
\title{Multi-calibration pipeline (for survival models)}
\usage{
ppl_mcboostsurv(learner = lrn("surv.kaplan"), param_vals = list())
}
\arguments{
\item{learner}{(mlr3)\code{\link[mlr3:Learner]{mlr3::Learner}}\cr
Initial learner.
Defaults to \code{lrn("surv.kaplan")}.
Note: An initial predictor can also be supplied via the \code{init_predictor} parameter.
The learner is internally wrapped into a \code{PipeOpLearnerCV}
with \code{resampling.method = "insample"} as a default.
All parameters can be adjusted through the resulting Graph's \code{param_set}.}

\item{param_vals}{\code{list} \cr
List of parameter values passed on to \code{MCBoostSurv$new}}
}
\value{
(mlr3pipelines) \code{\link{Graph}}
}
\description{
Wraps MCBoostSurv in a Pipeline to be used with \code{mlr3pipelines}.
For now this assumes training on the same dataset that is later used
for multi-calibration.
}
\examples{
library("mlr3pipelines")
gr = ppl_mcboostsurv()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCBoost.R
\name{MCBoost}
\alias{MCBoost}
\title{Multi-Calibration Boosting}
\description{
Implements Multi-Calibration Boosting by Hebert-Johnson et al. (2018) and
Multi-Accuracy Boosting by Kim et al. (2019) for the multi-calibration of a
machine learning model's prediction.
Multi-Calibration works best in scenarios where the underlying data & labels are unbiased
but a bias is introduced within the algorithm's fitting procedure. This is often the case,
e.g. when an algorithm fits a majority population while ignoring or under-fitting minority
populations.\cr
Expects initial models that fit binary outcomes or continuous outcomes with
predictions that are in (or scaled to) the 0-1 range.
The method defaults to \verb{Multi-Accuracy Boosting} as described in Kim et al. (2019).
In order to obtain behaviour as described in Hebert-Johnson et al. (2018) set
\code{multiplicative=FALSE} and \code{num_buckets} to 10.
\itemize{
For additional details, please refer to the relevant publications:
\item{Hebert-Johnson et al., 2018. Multicalibration: Calibration for the (Computationally-Identifiable) Masses.
Proceedings of the 35th International Conference on Machine Learning, PMLR 80:1939-1948.
https://proceedings.mlr.press/v80/hebert-johnson18a.html.}{}
\item{Kim et al., 2019. Multiaccuracy: Black-Box Post-Processing for Fairness in Classification.
Proceedings of the 2019 AAAI/ACM Conference on AI, Ethics, and Society (AIES '19).
Association for Computing Machinery, New York, NY, USA, 247–254.
https://dl.acm.org/doi/10.1145/3306618.3314287}{}
}
}
\examples{
# See vignette for more examples.
# Instantiate the object
\dontrun{
mc = MCBoost$new()
# Run multi-calibration on training dataset.
mc$multicalibrate(iris[1:100, 1:4], factor(sample(c("A", "B"), 100, TRUE)))
# Predict on test set
mc$predict_probs(iris[101:150, 1:4])
# Get auditor effect
mc$auditor_effect(iris[101:150, 1:4])
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{max_iter}}{\code{\link{integer}} \cr
The maximum number of iterations of the multi-calibration/multi-accuracy method.}

\item{\code{alpha}}{\code{\link{numeric}} \cr
Accuracy parameter that determines the stopping condition.}

\item{\code{eta}}{\code{\link{numeric}} \cr
Parameter for multiplicative weight update (step size).}

\item{\code{num_buckets}}{\code{\link{integer}} \cr
The number of buckets to split into in addition to using the whole sample.}

\item{\code{bucket_strategy}}{\code{\link{character}} \cr
Currently only supports "simple", even split along probabilities.
Only relevant for \code{num_buckets} > 1.}

\item{\code{rebucket}}{\code{\link{logical}} \cr
Should buckets be re-calculated at each iteration?}

\item{\code{eval_fulldata}}{\code{\link{logical}} \cr
Should auditor be evaluated on the full data?}

\item{\code{partition}}{\code{\link{logical}} \cr
True/False flag for whether to split up predictions by their "partition"
(e.g., predictions less than 0.5 and predictions greater than 0.5).}

\item{\code{multiplicative}}{\code{\link{logical}} \cr
Specifies the strategy for updating the weights (multiplicative weight vs additive).}

\item{\code{iter_sampling}}{\code{\link{character}} \cr
Specifies the strategy to sample the validation data for each iteration.}

\item{\code{auditor_fitter}}{\code{\link{AuditorFitter}} \cr
Specifies the type of model used to fit the residuals.}

\item{\code{predictor}}{\code{\link{function}} \cr
Initial predictor function.}

\item{\code{iter_models}}{\code{\link{list}} \cr
Cumulative list of fitted models.}

\item{\code{iter_partitions}}{\code{\link{list}} \cr
Cumulative list of data partitions for models.}

\item{\code{iter_corr}}{\code{\link{list}} \cr
Auditor correlation in each iteration.}

\item{\code{auditor_effects}}{\code{\link{list}} \cr
Auditor effect in each iteration.}

\item{\code{bucket_strategies}}{\code{\link{character}} \cr
Possible bucket_strategies in McBoostSurv.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{MCBoost$new()}}
\item \href{#method-multicalibrate}{\code{MCBoost$multicalibrate()}}
\item \href{#method-predict_probs}{\code{MCBoost$predict_probs()}}
\item \href{#method-auditor_effect}{\code{MCBoost$auditor_effect()}}
\item \href{#method-print}{\code{MCBoost$print()}}
\item \href{#method-clone}{\code{MCBoost$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a multi-calibration instance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$new(
  max_iter = 5,
  alpha = 1e-04,
  eta = 1,
  num_buckets = 2,
  partition = ifelse(num_buckets > 1, TRUE, FALSE),
  bucket_strategy = "simple",
  rebucket = FALSE,
  eval_fulldata = FALSE,
  multiplicative = TRUE,
  auditor_fitter = NULL,
  subpops = NULL,
  default_model_class = ConstantPredictor,
  init_predictor = NULL,
  iter_sampling = "none"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{max_iter}}{\code{\link{integer}} \cr
The maximum number of iterations of the multi-calibration/multi-accuracy method.
Default \code{5L}.}

\item{\code{alpha}}{\code{\link{numeric}} \cr
Accuracy parameter that determines the stopping condition. Default \code{1e-4}.}

\item{\code{eta}}{\code{\link{numeric}} \cr
Parameter for multiplicative weight update (step size). Default \code{1.0}.}

\item{\code{num_buckets}}{\code{\link{integer}} \cr
The number of buckets to split into in addition to using the whole sample. Default \code{2L}.}

\item{\code{partition}}{\code{\link{logical}} \cr
True/False flag for whether to split up predictions by their "partition"
(e.g., predictions less than 0.5 and predictions greater than 0.5).
Defaults to \code{TRUE} (multi-accuracy boosting).}

\item{\code{bucket_strategy}}{\code{\link{character}} \cr
Currently only supports "simple", even split along probabilities.
Only taken into account for \code{num_buckets} > 1.}

\item{\code{rebucket}}{\code{\link{logical}} \cr
Should buckets be re-done at each iteration? Default \code{FALSE}.}

\item{\code{eval_fulldata}}{\code{\link{logical}} \cr
Should the auditor be evaluated on the full data or on the respective bucket for determining
the stopping criterion? Default \code{FALSE}, auditor is only evaluated on the bucket.
This setting keeps the implementation closer to the Algorithm proposed in the corresponding
multi-accuracy paper (Kim et al., 2019) where auditor effects are computed across the full
sample (i.e. eval_fulldata = TRUE).}

\item{\code{multiplicative}}{\code{\link{logical}} \cr
Specifies the strategy for updating the weights (multiplicative weight vs additive).
Defaults to \code{TRUE} (multi-accuracy boosting). Set to \code{FALSE} for multi-calibration.}

\item{\code{auditor_fitter}}{\code{\link{AuditorFitter}}|\code{\link{character}}|\code{\link[mlr3:Learner]{mlr3::Learner}} \cr
Specifies the type of model used to fit the
residuals. The default is \code{\link{RidgeAuditorFitter}}.
Can be a \code{character}, the name of a \code{\link{AuditorFitter}}, a \code{\link[mlr3:Learner]{mlr3::Learner}} that is then
auto-converted into a \code{\link{LearnerAuditorFitter}} or a custom \code{\link{AuditorFitter}}.}

\item{\code{subpops}}{\code{\link{list}} \cr
Specifies a collection of characteristic attributes
and the values they take to define subpopulations
e.g. list(age = c('20-29','30-39','40+'), nJobs = c(0,1,2,'3+'), ,..).}

\item{\code{default_model_class}}{\code{Predictor} \cr
The class of the model that should be used as the init predictor model if
\code{init_predictor} is not specified. Defaults to \code{ConstantPredictor} which
predicts a constant value.}

\item{\code{init_predictor}}{\code{\link{function}}|\code{\link[mlr3:Learner]{mlr3::Learner}} \cr
The initial predictor function to use (i.e., if the user has a pretrained model).
If a \code{mlr3} \code{Learner} is passed, it will be autoconverted using \code{mlr3_init_predictor}.
This requires the \code{\link[mlr3:Learner]{mlr3::Learner}} to be trained.}

\item{\code{iter_sampling}}{\code{\link{character}} \cr
How to sample the validation data for each iteration?
Can be \code{bootstrap}, \code{split} or \code{none}.\cr
"split" splits the data into \code{max_iter} parts and validates on each sample in each iteration.\cr
"bootstrap" uses a new bootstrap sample in each iteration.\cr
"none" uses the same dataset in each iteration.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-multicalibrate"></a>}}
\if{latex}{\out{\hypertarget{method-multicalibrate}{}}}
\subsection{Method \code{multicalibrate()}}{
Run multi-calibration.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$multicalibrate(data, labels, predictor_args = NULL, audit = FALSE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{labels}}{\code{\link{numeric}}\cr
One-hot encoded labels (of same length as data).}

\item{\code{predictor_args}}{\code{\link{any}} \cr
Arguments passed on to \code{init_predictor}. Defaults to \code{NULL}.}

\item{\code{audit}}{\code{\link{logical}} \cr
Perform auditing? Initialized to \code{TRUE}.}

\item{\code{...}}{\code{\link{any}} \cr
Params passed on to other methods.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{NULL}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-predict_probs"></a>}}
\if{latex}{\out{\hypertarget{method-predict_probs}{}}}
\subsection{Method \code{predict_probs()}}{
Predict a dataset with multi-calibrated predictions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$predict_probs(x, t = Inf, predictor_args = NULL, audit = FALSE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\code{\link{data.table}} \cr
Prediction data.}

\item{\code{t}}{\code{\link{integer}} \cr
Number of multi-calibration steps to predict. Default: \code{Inf} (all).}

\item{\code{predictor_args}}{\code{\link{any}} \cr
Arguments passed on to \code{init_predictor}. Defaults to \code{NULL}.}

\item{\code{audit}}{\code{\link{logical}} \cr
Should audit weights be stored? Default \code{FALSE}.}

\item{\code{...}}{\code{\link{any}} \cr
Params passed on to the residual prediction model's predict method.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{\link{numeric}}\cr
Numeric vector of multi-calibrated predictions.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-auditor_effect"></a>}}
\if{latex}{\out{\hypertarget{method-auditor_effect}{}}}
\subsection{Method \code{auditor_effect()}}{
Compute the auditor effect for each instance which are the cumulative
absolute predictions of the auditor. It indicates "how much"
each observation was affected by multi-calibration on average across iterations.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$auditor_effect(
  x,
  aggregate = TRUE,
  t = Inf,
  predictor_args = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\code{\link{data.table}} \cr
Prediction data.}

\item{\code{aggregate}}{\code{\link{logical}} \cr
Should the auditor effect be aggregated across iterations? Defaults to \code{TRUE}.}

\item{\code{t}}{\code{\link{integer}} \cr
Number of multi-calibration steps to predict. Defaults to \code{Inf} (all).}

\item{\code{predictor_args}}{\code{\link{any}} \cr
Arguments passed on to \code{init_predictor}. Defaults to \code{NULL}.}

\item{\code{...}}{\code{\link{any}} \cr
Params passed on to the residual prediction model's predict method.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{\link{numeric}} \cr
Numeric vector of auditor effects for each row in \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Prints information about multi-calibration.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{\code{any}\cr
Not used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCBoost$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{package}
\name{mcboost-package}
\alias{mcboost}
\alias{mcboost-package}
\title{mcboost: Multi-Calibration Boosting}
\description{
Implements 'Multi-Calibration Boosting' (2018) <https://proceedings.mlr.press/v80/hebert-johnson18a.html> and 'Multi-Accuracy Boosting' (2019) <arXiv:1805.12317> for the multi-calibration of a machine learning model's prediction. 'MCBoost' updates predictions for sub-groups in an iterative fashion in order to mitigate biases like poor calibration or large accuracy differences across subgroups. Multi-Calibration works best in scenarios where the underlying data & labels are unbiased, but resulting models are. This is often the case, e.g. when an algorithm fits a majority population while ignoring or under-fitting minority populations.
}
\references{
Kim et al., 2019: Multiaccuracy: Black-Box Post-Processing for Fairness in Classification.
Hebert-Johnson et al., 2018: Multicalibration: Calibration for the ({C}omputationally-Identifiable) Masses.
Pfisterer F, Kern C, Dandl S, Sun M, Kim M, Bischl B (2021).
\dQuote{mcboost: Multi-Calibration Boosting for R.}
\emph{Journal of Open Source Software}, \bold{6}(64), 3453.
\doi{10.21105/joss.03453}, \url{https://joss.theoj.org/papers/10.21105/joss.03453}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/mlr-org/mcboost}
  \item Report bugs at \url{https://github.com/mlr-org/mcboost/issues}
}

}
\author{
\strong{Maintainer}: Florian Pfisterer \email{pfistererf@googlemail.com} (\href{https://orcid.org/0000-0001-8867-762X}{ORCID})

Other contributors:
\itemize{
  \item Susanne Dandl \email{susanne.dandl@stat.uni-muenchen.de} (\href{https://orcid.org/0000-0003-4324-4163}{ORCID}) [contributor]
  \item Christoph Kern \email{c.kern@uni-mannheim.de} (\href{https://orcid.org/0000-0001-7363-4299}{ORCID}) [contributor]
  \item Carolin Becker [contributor]
  \item Bernd Bischl \email{bernd_bischl@gmx.net} (\href{https://orcid.org/0000-0001-6002-6980}{ORCID}) [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PipelineMCBoost.R
\name{ppl_mcboost}
\alias{ppl_mcboost}
\title{Multi-calibration pipeline}
\usage{
ppl_mcboost(learner = lrn("classif.featureless"), param_vals = list())
}
\arguments{
\item{learner}{(mlr3)\code{\link[mlr3:Learner]{mlr3::Learner}}\cr
Initial learner. Internally wrapped into a \code{PipeOpLearnerCV}
with \code{resampling.method = "insample"} as a default.
All parameters can be adjusted through the resulting Graph's \code{param_set}.
Defaults to \code{lrn("classif.featureless")}.
Note: An initial predictor can also be supplied via the \code{init_predictor} parameter.}

\item{param_vals}{\code{list} \cr
List of parameter values passed on to \code{MCBoost$new}.}
}
\value{
(mlr3pipelines) \code{\link{Graph}}
}
\description{
Wraps MCBoost in a Pipeline to be used with \code{mlr3pipelines}.
For now this assumes training on the same dataset that is later used
for multi-calibration.
}
\examples{
  \dontrun{
  library("mlr3pipelines")
  gr = ppl_mcboost()
  }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PipeOpLearnerPred.R, R/PipeOpMCBoost.R
\name{mlr_pipeops_mcboost}
\alias{mlr_pipeops_mcboost}
\alias{PipeOpLearnerPred}
\alias{PipeOpMCBoost}
\title{Multi-Calibrate a Learner's Prediction}
\format{
\code{\link{R6Class}} inheriting from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.

\code{\link{R6Class}} inheriting from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}
\description{


Post-process a learner prediction using multi-calibration.
For more details, please refer to \url{https://arxiv.org/pdf/1805.12317.pdf} (Kim et al. 2018)
or the help for \code{\link{MCBoost}}.
If no \code{init_predictor} is provided, the preceding learner's predictions
corresponding to the \code{prediction} slot are used as an initial predictor for \code{MCBoost}.
}
\section{Construction}{
\preformatted{PipeOpLearnerPred$new(learner, id = NULL, param_vals = list())

* `learner` :: [`Learner`][mlr3::Learner] \\cr
  [`Learner`][mlr3::Learner] to  prediction, or a string identifying a
  [`Learner`][mlr3::Learner] in the [`mlr3::mlr_learners`] [`Dictionary`][mlr3misc::Dictionary].
* `id` :: `character(1)`
  Identifier of the resulting object, internally defaulting to the `id` of the [`Learner`][mlr3::Learner] being wrapped.
* `param_vals` :: named `list`\\cr
  List of hyperparameter settings, overwriting the hyperparameter settings that would otherwise be set during construction. Default `list()`.


[mlr3::Learner]: R:mlr3::Learner
[mlr3::Learner]: R:mlr3::Learner
[mlr3::Learner]: R:mlr3::Learner
[`mlr3::mlr_learners`]: R:\%60mlr3::mlr_learners\%60
[mlr3misc::Dictionary]: R:mlr3misc::Dictionary
[mlr3::Learner]: R:mlr3::Learner
}

\preformatted{PipeOpMCBoost$new(id = "mcboost", param_vals = list())
}
\itemize{
\item \code{id} :: \code{character(1)}
Identifier of the resulting  object, default \code{"threshold"}.
\item \code{param_vals} :: named \code{list}\cr
List of hyperparameter settings, overwriting the hyperparameter settings that would otherwise be set during construction.
See \code{MCBoost} for a comprehensive description of all hyperparameters.
}
}

\section{Input and Output Channels}{

\code{\link{PipeOpLearnerPred}} has one input channel named \code{"input"}, taking a \code{\link[mlr3:Task]{Task}} specific to the \code{\link[mlr3:Learner]{Learner}}
type given to \code{learner} during construction; both during training and prediction.

\code{\link{PipeOpLearnerPred}} has one output channel named \code{"output"}, producing a \code{\link[mlr3:Task]{Task}} specific to the \code{\link[mlr3:Learner]{Learner}}
type given to \code{learner} during construction; both during training and prediction.


During training, the input and output are \code{"data"} and \code{"prediction"}, two \code{\link[mlr3:TaskClassif]{TaskClassif}}.
A \code{\link[mlr3:PredictionClassif]{PredictionClassif}} is required as input and returned as output during prediction.
}

\section{State}{



The \verb{$state} is a \code{MCBoost} Object as obtained from \code{MCBoost$new()}.
}

\section{Parameters}{

The \verb{$state} is set to the \verb{$state} slot of the \code{\link[mlr3:Learner]{Learner}} object, together with the \verb{$state} elements inherited from
\code{\link[mlr3pipelines:PipeOpTaskPreproc]{mlr3pipelines::PipeOpTaskPreproc}}. It is a named \code{list} with the inherited members, as well as:
\itemize{
\item \code{model} :: \code{any}\cr
Model created by the \code{\link[mlr3:Learner]{Learner}}'s \verb{$.train()} function.
\item \code{train_log} :: \code{\link{data.table}} with columns \code{class} (\code{character}), \code{msg} (\code{character})\cr
Errors logged during training.
\item \code{train_time} :: \code{numeric(1)}\cr
Training time, in seconds.
\item \code{predict_log} :: \code{NULL} | \code{\link{data.table}} with columns \code{class} (\code{character}), \code{msg} (\code{character})\cr
Errors logged during prediction.
\item \code{predict_time} :: \code{NULL} | \code{numeric(1)}
Prediction time, in seconds.
}


\itemize{
\item \code{max_iter} :: \code{integer}\cr
A integer specifying the number of multi-calibration rounds. Defaults to 5.
}
}

\section{Fields}{

Fields inherited from \code{\link{PipeOp}}, as well as:
\itemize{
\item \code{learner} :: \code{\link[mlr3:Learner]{Learner}}\cr
\code{\link[mlr3:Learner]{Learner}} that is being wrapped. Read-only.
\item \code{learner_model} :: \code{\link[mlr3:Learner]{Learner}}\cr
\code{\link[mlr3:Learner]{Learner}} that is being wrapped. This learner contains the model if the \code{PipeOp} is trained. Read-only.
}


Only fields inherited from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}

\section{Methods}{

Methods inherited from \code{\link[mlr3pipelines:PipeOpTaskPreproc]{mlr3pipelines::PipeOpTaskPreproc}}/\code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.


Only methods inherited from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}

\examples{
\dontrun{
gr = gunion(list(
  "data" = po("nop"),
  "prediction" = po("learner_cv", lrn("classif.rpart"))
)) \%>>\%
  PipeOpMCBoost$new()
tsk = tsk("sonar")
tid = sample(1:208, 108)
gr$train(tsk$clone()$filter(tid))
gr$predict(tsk$clone()$filter(setdiff(1:208, tid)))
}
}
\seealso{
https://mlr3book.mlr-org.com/list-pipeops.html

https://mlr3book.mlr-org.com/list-pipeops.html

Other PipeOps: 
\code{\link{mlr_pipeops_mcboostsurv}}

Other PipeOps: 
\code{\link{mlr_pipeops_mcboostsurv}}
}
\concept{PipeOps}
\section{Super classes}{
\code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}} -> \code{\link[mlr3pipelines:PipeOpTaskPreproc]{mlr3pipelines::PipeOpTaskPreproc}} -> \code{PipeOpLearnerPred}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{learner}}{The wrapped learner.}

\item{\code{learner_model}}{The wrapped learner's model(s).}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PipeOpLearnerPred$new()}}
\item \href{#method-clone}{\code{PipeOpLearnerPred$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="predict">}\href{../../mlr3pipelines/html/PipeOp.html#method-predict}{\code{mlr3pipelines::PipeOp$predict()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="print">}\href{../../mlr3pipelines/html/PipeOp.html#method-print}{\code{mlr3pipelines::PipeOp$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="train">}\href{../../mlr3pipelines/html/PipeOp.html#method-train}{\code{mlr3pipelines::PipeOp$train()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a Learner Predictor PipeOp. Can be used to wrap trained or untrainted
mlr3 learners.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpLearnerPred$new(learner, id = NULL, param_vals = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{learner}}{\code{\link{Learner}}\cr
The learner that should be wrapped.}

\item{\code{id}}{\code{\link{character}} \cr
The \code{PipeOp}'s id. Defaults to "mcboost".}

\item{\code{param_vals}}{\code{\link{list}} \cr
List of hyperparameters for the \code{PipeOp}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpLearnerPred$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
\section{Super class}{
\code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}} -> \code{PipeOpMCBoost}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{predict_type}}{Predict type of the PipeOp.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PipeOpMCBoost$new()}}
\item \href{#method-clone}{\code{PipeOpMCBoost$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="predict">}\href{../../mlr3pipelines/html/PipeOp.html#method-predict}{\code{mlr3pipelines::PipeOp$predict()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="print">}\href{../../mlr3pipelines/html/PipeOp.html#method-print}{\code{mlr3pipelines::PipeOp$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="train">}\href{../../mlr3pipelines/html/PipeOp.html#method-train}{\code{mlr3pipelines::PipeOp$train()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a Multi-Calibration PipeOp.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpMCBoost$new(id = "mcboost", param_vals = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{\code{\link{character}} \cr
The \code{PipeOp}'s id. Defaults to "mcboost".}

\item{\code{param_vals}}{\code{\link{list}} \cr
List of hyperparameters for the \code{PipeOp}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpMCBoost$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PipeOpMCBoostSurv.R
\name{mlr_pipeops_mcboostsurv}
\alias{mlr_pipeops_mcboostsurv}
\alias{PipeOpMCBoostSurv}
\title{Multi-Calibrate a Learner's Prediction (Survival Model)}
\format{
\code{\link{R6Class}} inheriting from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}
\description{
Post-process a survival learner prediction using multi-calibration.
For more details, please refer to \url{https://arxiv.org/pdf/1805.12317.pdf} (Kim et al. 2018)
or the help for \code{\link{MCBoostSurv}}.
If no \code{init_predictor} is provided, the preceding learner's predictions
corresponding to the \code{prediction} slot are used as an initial predictor for \code{MCBoostSurv}.
}
\section{Construction}{
\preformatted{PipeOpMCBoostSurv$new(id = "mcboostsurv", param_vals = list())
}
\itemize{
\item \code{id} :: \code{character(1)}
Identifier of the resulting  object, default \code{"threshold"}.
\item \code{param_vals} :: named \code{list}\cr
List of hyperparameter settings, overwriting the hyperparameter settings that would otherwise be set during construction.
See \code{MCBoostSurv} for a comprehensive description of all hyperparameters.
}
}

\section{Input and Output Channels}{

During training, the input and output are \code{"data"} and \code{"prediction"}, two \code{\link[mlr3proba:TaskSurv]{TaskSurv}}.
A \code{\link[mlr3proba:PredictionSurv]{PredictionSurv}} is required as input and returned as output during prediction.
}

\section{State}{

The \verb{$state} is a \code{MCBoostSurv} Object as obtained from \code{MCBoostSurv$new()}.
}

\section{Parameters}{

\itemize{
\item \code{max_iter} :: \code{integer}\cr
A integer specifying the number of multi-calibration rounds. Defaults to 5.
}
}

\section{Fields}{

Only fields inherited from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}

\section{Methods}{

Only methods inherited from \code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}}.
}

\examples{
library(mlr3)
library(mlr3pipelines)
# Attention: gunion inputs have to be in the correct order for now.
\dontrun{
gr = gunion(list(
  "data" = po("nop"),
  "prediction" = po("learner_pred", lrn("surv.ranger"))
)) \%>>\%
  PipeOpMCBoostSurv$new()
tsk = tsk("rats")
tid = sample(1:300, 100)
gr$train(tsk$clone()$filter(tid))
gr$predict(tsk$clone()$filter(setdiff(1:300, tid)))
}
}
\seealso{
https://mlr3book.mlr-org.com/list-pipeops.html

Other PipeOps: 
\code{\link{mlr_pipeops_mcboost}}
}
\concept{PipeOps}
\section{Super class}{
\code{\link[mlr3pipelines:PipeOp]{mlr3pipelines::PipeOp}} -> \code{PipeOpMCBoostSurv}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{predict_type}}{Predict type of the PipeOp.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PipeOpMCBoostSurv$new()}}
\item \href{#method-clone}{\code{PipeOpMCBoostSurv$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="predict">}\href{../../mlr3pipelines/html/PipeOp.html#method-predict}{\code{mlr3pipelines::PipeOp$predict()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="print">}\href{../../mlr3pipelines/html/PipeOp.html#method-print}{\code{mlr3pipelines::PipeOp$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mlr3pipelines" data-topic="PipeOp" data-id="train">}\href{../../mlr3pipelines/html/PipeOp.html#method-train}{\code{mlr3pipelines::PipeOp$train()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a Multi-Calibration PipeOp (for Survival).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpMCBoostSurv$new(id = "mcboostsurv", param_vals = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{\code{\link{character}} \cr
The \code{PipeOp}'s id. Defaults to "mcboostsurv".}

\item{\code{param_vals}}{\code{\link{list}} \cr
List of hyperparameters for the \code{PipeOp}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PipeOpMCBoostSurv$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuditorFitters.R
\name{SubpopAuditorFitter}
\alias{SubpopAuditorFitter}
\title{Static AuditorFitter based on Subpopulations}
\value{
\code{\link{AuditorFitter}}\cr

\code{list} with items\cr
\itemize{
\item \code{corr}: pseudo-correlation between residuals and learner prediction.
\item \code{l}: the trained learner.
}
}
\description{
Used to assess multi-calibration based on a list of
binary valued columns: \code{subpops} passed during initialization.
}
\examples{
  \dontrun{
  library("data.table")
  data = data.table(
    "AGE_NA" = c(0, 0, 0, 0, 0),
    "AGE_0_10" =  c(1, 1, 0, 0, 0),
    "AGE_11_20" = c(0, 0, 1, 0, 0),
    "AGE_21_31" = c(0, 0, 0, 1, 1),
    "X1" = runif(5),
    "X2" = runif(5)
  )
  label = c(1,0,0,1,1)
  pops = list("AGE_NA", "AGE_0_10", "AGE_11_20", "AGE_21_31", function(x) {x[["X1" > 0.5]]})
  sf = SubpopAuditorFitter$new(subpops = pops)
  sf$fit(data, label - 0.5)
  }
}
\seealso{
Other AuditorFitter: 
\code{\link{CVLearnerAuditorFitter}},
\code{\link{LearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}}
}
\concept{AuditorFitter}
\section{Super class}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{SubpopAuditorFitter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{subpops}}{\code{\link{list}} \cr
List of subpopulation indicators.
Initialize a SubpopAuditorFitter}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{SubpopAuditorFitter$new()}}
\item \href{#method-fit}{\code{SubpopAuditorFitter$fit()}}
\item \href{#method-clone}{\code{SubpopAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initializes a \code{\link{SubpopAuditorFitter}} that
assesses multi-calibration within each group defined
by the \verb{subpops'. Names in }subpops` must correspond to
columns in the data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubpopAuditorFitter$new(subpops)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{subpops}}{\code{\link{list}} \cr
Specifies a collection of characteristic attributes
and the values they take to define subpopulations
e.g. list(age = c('20-29','30-39','40+'), nJobs = c(0,1,2,'3+'), ,..).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the learner and compute correlation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubpopAuditorFitter$fit(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubpopAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuditorFitters.R
\name{SubgroupAuditorFitter}
\alias{SubgroupAuditorFitter}
\title{Static AuditorFitter based on Subgroups}
\value{
\code{\link{AuditorFitter}}\cr

\code{list} with items\cr
\itemize{
\item \code{corr}: pseudo-correlation between residuals and learner prediction.
\item \code{l}: the trained learner.
}
}
\description{
Used to assess multi-calibration based on a list of
binary \code{subgroup_masks} passed during initialization.
}
\examples{
 \dontrun{
 library("data.table")
 data = data.table(
   "AGE_0_10" =  c(1, 1, 0, 0, 0),
   "AGE_11_20" = c(0, 0, 1, 0, 0),
   "AGE_21_31" = c(0, 0, 0, 1, 1),
   "X1" = runif(5),
   "X2" = runif(5)
 )
 label = c(1,0,0,1,1)
 masks = list(
   "M1" = c(1L, 0L, 1L, 1L, 0L),
   "M2" = c(1L, 0L, 0L, 0L, 1L)
 )
 sg = SubgroupAuditorFitter$new(masks)
 }
}
\seealso{
Other AuditorFitter: 
\code{\link{CVLearnerAuditorFitter}},
\code{\link{LearnerAuditorFitter}},
\code{\link{SubpopAuditorFitter}}
}
\concept{AuditorFitter}
\section{Super class}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{SubgroupAuditorFitter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{subgroup_masks}}{\code{\link{list}} \cr
List of subgroup masks.
Initialize a SubgroupAuditorFitter}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{SubgroupAuditorFitter$new()}}
\item \href{#method-fit}{\code{SubgroupAuditorFitter$fit()}}
\item \href{#method-clone}{\code{SubgroupAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initializes a \code{\link{SubgroupAuditorFitter}} that
assesses multi-calibration within each group defined
by the `subpops'.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubgroupAuditorFitter$new(subgroup_masks)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{subgroup_masks}}{\code{\link{list}} \cr
List of subgroup masks. Subgroup masks are list(s) of integer masks,
each with the same length as data to be fitted on.
They allow defining subgroups of the data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the learner and compute correlation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubgroupAuditorFitter$fit(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SubgroupAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{mlr3_init_predictor}
\alias{mlr3_init_predictor}
\title{Create an initial predictor function from a trained mlr3 learner}
\usage{
mlr3_init_predictor(learner)
}
\arguments{
\item{learner}{\code{\link[mlr3:Learner]{mlr3::Learner}}
A trained learner used for initialization.}
}
\value{
\code{\link{function}}
}
\description{
Create an initial predictor function from a trained mlr3 learner
}
\examples{
 \dontrun{
 library("mlr3")
 l = lrn("classif.featureless")$train(tsk("sonar"))
 mlr3_init_predictor(l)
 }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{one_hot}
\alias{one_hot}
\title{One-hot encode a factor variable}
\usage{
one_hot(labels)
}
\arguments{
\item{labels}{\code{\link{factor}}\cr
Factor to encode.}
}
\value{
\code{\link{integer}}\cr
Integer vector of encoded labels.
}
\description{
One-hot encode a factor variable
}
\examples{
 \dontrun{
 one_hot(factor(c("a", "b", "a")))
 }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuditorFitters.R
\name{CVLearnerAuditorFitter}
\alias{CVLearnerAuditorFitter}
\alias{CVTreeAuditorFitter}
\alias{CVRidgeAuditorFitter}
\title{Cross-validated AuditorFitter from a Learner}
\value{
\code{\link{AuditorFitter}}\cr

\code{list} with items\cr
\itemize{
\item \code{corr}: pseudo-correlation between residuals and learner prediction.
\item \code{l}: the trained learner.
}
}
\description{
CVLearnerAuditorFitter returns the cross-validated predictions
instead of the in-sample predictions.

Available data is cut into complementary subsets (folds).
For each subset out-of-sample predictions are received by training a model
on all other subsets and predicting afterwards on the left-out subset.
}
\section{Functions}{
\itemize{
\item \code{CVTreeAuditorFitter}: Cross-Validated auditor based on rpart

\item \code{CVRidgeAuditorFitter}: Cross-Validated auditor based on glmnet
}}

\seealso{
Other AuditorFitter: 
\code{\link{LearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}

Other AuditorFitter: 
\code{\link{LearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}

Other AuditorFitter: 
\code{\link{LearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}
}
\concept{AuditorFitter}
\section{Super class}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{CVLearnerAuditorFitter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{learner}}{\code{CVLearnerPredictor}\cr
Learner used for fitting residuals.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{CVLearnerAuditorFitter$new()}}
\item \href{#method-fit}{\code{CVLearnerAuditorFitter$fit()}}
\item \href{#method-clone}{\code{CVLearnerAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define a \code{CVAuditorFitter} from a learner.
Available instantiations:\cr \code{\link{CVTreeAuditorFitter}} (rpart) and
\code{\link{CVRidgeAuditorFitter}} (glmnet).
See \code{\link[mlr3pipelines:mlr_pipeops_learner_cv]{mlr3pipelines::PipeOpLearnerCV}} for more information on
cross-validated learners.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVLearnerAuditorFitter$new(learner, folds = 3L)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{learner}}{\code{\link[mlr3:Learner]{mlr3::Learner}}\cr
Regression Learner to use.}

\item{\code{folds}}{\code{\link{integer}}\cr
Number of folds to use for PipeOpLearnerCV. Defaults to 3.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the cross-validated learner and compute correlation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVLearnerAuditorFitter$fit(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVLearnerAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
\section{Super classes}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{\link[mcboost:CVLearnerAuditorFitter]{mcboost::CVLearnerAuditorFitter}} -> \code{CVTreeAuditorFitter}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{CVTreeAuditorFitter$new()}}
\item \href{#method-clone}{\code{CVTreeAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="CVLearnerAuditorFitter" data-id="fit">}\href{../../mcboost/html/CVLearnerAuditorFitter.html#method-fit}{\code{mcboost::CVLearnerAuditorFitter$fit()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define a cross-validated AuditorFitter from a rpart learner
See \code{\link[mlr3pipelines:mlr_pipeops_learner_cv]{mlr3pipelines::PipeOpLearnerCV}} for more information on
cross-validated learners.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVTreeAuditorFitter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVTreeAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
\section{Super classes}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{\link[mcboost:CVLearnerAuditorFitter]{mcboost::CVLearnerAuditorFitter}} -> \code{CVRidgeAuditorFitter}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{CVRidgeAuditorFitter$new()}}
\item \href{#method-clone}{\code{CVRidgeAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="CVLearnerAuditorFitter" data-id="fit">}\href{../../mcboost/html/CVLearnerAuditorFitter.html#method-fit}{\code{mcboost::CVLearnerAuditorFitter$fit()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define a cross-validated AuditorFitter from a glmnet learner.
See \code{\link[mlr3pipelines:mlr_pipeops_learner_cv]{mlr3pipelines::PipeOpLearnerCV}} for more information on
cross-validated learners.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVRidgeAuditorFitter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CVRidgeAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{make_survival_curve}
\alias{make_survival_curve}
\title{Make every row monotonically decreasing in order to obtain the survival property.
Additionally, many predicitions need 1 as a first value and 0 as a last value.
(e.g. \code{PredictionSurv} needs this attribute.)}
\usage{
make_survival_curve(prediction)
}
\arguments{
\item{prediction}{\code{\link{data.table}}
Data.table with predictions. Every row is survival probability for the corresponding time.
Every column corresponds to  a specific time point.}
}
\value{
\code{\link{data.table}}
}
\description{
Make every row monotonically decreasing in order to obtain the survival property.
Additionally, many predicitions need 1 as a first value and 0 as a last value.
(e.g. \code{PredictionSurv} needs this attribute.)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuditorFitters.R
\name{LearnerAuditorFitter}
\alias{LearnerAuditorFitter}
\alias{TreeAuditorFitter}
\alias{RidgeAuditorFitter}
\title{Create an AuditorFitter from a Learner}
\value{
\code{\link{AuditorFitter}}\cr

\code{list} with items\cr
\itemize{
\item \code{corr}: pseudo-correlation between residuals and learner prediction.
\item \code{l}: the trained learner.
}
}
\description{
Instantiates an AuditorFitter that trains a \code{\link[mlr3:Learner]{mlr3::Learner}}
on the data.
}
\section{Functions}{
\itemize{
\item \code{TreeAuditorFitter}: Learner auditor based on rpart

\item \code{RidgeAuditorFitter}: Learner auditor based on glmnet
}}

\seealso{
Other AuditorFitter: 
\code{\link{CVLearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}

Other AuditorFitter: 
\code{\link{CVLearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}

Other AuditorFitter: 
\code{\link{CVLearnerAuditorFitter}},
\code{\link{SubgroupAuditorFitter}},
\code{\link{SubpopAuditorFitter}}
}
\concept{AuditorFitter}
\section{Super class}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{LearnerAuditorFitter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{learner}}{\code{LearnerPredictor}\cr
Learner used for fitting residuals.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{LearnerAuditorFitter$new()}}
\item \href{#method-fit}{\code{LearnerAuditorFitter$fit()}}
\item \href{#method-clone}{\code{LearnerAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define an \code{AuditorFitter} from a Learner.
Available instantiations:\cr \code{\link{TreeAuditorFitter}} (rpart) and
\code{\link{RidgeAuditorFitter}} (glmnet).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LearnerAuditorFitter$new(learner)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{learner}}{\code{\link[mlr3:Learner]{mlr3::Learner}}\cr
Regression learner to use.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the learner and compute correlation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LearnerAuditorFitter$fit(data, resid, mask)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{\code{\link{data.table}}\cr
Features.}

\item{\code{resid}}{\code{\link{numeric}}\cr
Residuals (of same length as data).}

\item{\code{mask}}{\code{\link{integer}}\cr
Mask applied to the data. Only used for \code{SubgroupAuditorFitter}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LearnerAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
\section{Super classes}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{\link[mcboost:LearnerAuditorFitter]{mcboost::LearnerAuditorFitter}} -> \code{TreeAuditorFitter}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{TreeAuditorFitter$new()}}
\item \href{#method-clone}{\code{TreeAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="LearnerAuditorFitter" data-id="fit">}\href{../../mcboost/html/LearnerAuditorFitter.html#method-fit}{\code{mcboost::LearnerAuditorFitter$fit()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define a AuditorFitter from a rpart learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TreeAuditorFitter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TreeAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
\section{Super classes}{
\code{\link[mcboost:AuditorFitter]{mcboost::AuditorFitter}} -> \code{\link[mcboost:LearnerAuditorFitter]{mcboost::LearnerAuditorFitter}} -> \code{RidgeAuditorFitter}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RidgeAuditorFitter$new()}}
\item \href{#method-clone}{\code{RidgeAuditorFitter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="AuditorFitter" data-id="fit_to_resid">}\href{../../mcboost/html/AuditorFitter.html#method-fit_to_resid}{\code{mcboost::AuditorFitter$fit_to_resid()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="mcboost" data-topic="LearnerAuditorFitter" data-id="fit">}\href{../../mcboost/html/LearnerAuditorFitter.html#method-fit}{\code{mcboost::LearnerAuditorFitter$fit()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Define a AuditorFitter from a glmnet learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RidgeAuditorFitter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RidgeAuditorFitter$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
