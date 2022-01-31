# emba

<!-- badges: start -->
[![R build status](https://github.com/bblodfon/emba/workflows/R-CMD-check/badge.svg)](https://github.com/bblodfon/emba/actions)
[![codecov](https://codecov.io/gh/bblodfon/emba/branch/master/graph/badge.svg)](https://codecov.io/gh/bblodfon/emba)
[![CRAN status](https://www.r-pkg.org/badges/version/emba)](https://cran.r-project.org/package=emba)
[![Downloads](https://cranlogs.r-pkg.org/badges/emba)](https://cran.r-project.org/package=emba)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02583/status.svg)](https://doi.org/10.21105/joss.02583)
<!-- badges: end -->

Analysis and visualization of an ensemble of boolean models for biomarker discovery in cancer cell networks.

The package allows to easily load the simulation data results of the [DrugLogics](https://github.com/druglogics) software pipeline that is used to predict synergistic drug combinations in cancer cell lines.
It has generic functions that can be used to split a boolean model dataset to model groups with regards to the models predictive performance (number of *true positive* predictions/*Matthews correlation coefficient* score) or synergy prediction based on a given set of *gold standard* synergies and find the average activity difference per network node between all model group pairs.
Thus, given user-specific thresholds, important nodes (*biomarkers*) can be accessed in the sense that they make the models predict specific synergies (*synergy biomarkers*) or have better performance in general (*performance biomarkers*).

Lastly, if the boolean models have a [specific equation form](https://druglogics.github.io/druglogics-doc/gitsbe-description.html#default-equation) and differ only in their link operator, *link operator* biomarkers can also be found.

## Install

CRAN version:
```
install.packages("emba")
```

Development version:
```
remotes::install_github("bblodfon/emba")
```

## Usage

Check the [Get Started guide](https://bblodfon.github.io/emba/articles/emba.html).

For an earlier example usage of this package (version `0.1.1`), see this [analysis](https://druglogics.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/) performed on multiple boolean model datasets.

## Cite

- Formatted citation:

Zobolas et al., (2020). emba: R package for analysis and visualization of biomarkers in boolean model ensembles. Journal of Open Source Software, 5(53), 2583, https://doi.org/10.21105/joss.02583

- BibTeX citation:
```
@article{Zobolas2020,
  doi = {10.21105/joss.02583},
  url = {https://doi.org/10.21105/joss.02583},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {53},
  pages = {2583},
  author = {John Zobolas and Martin Kuiper and Åsmund Flobak},
  title = {emba: R package for analysis and visualization of biomarkers in boolean model ensembles},
  journal = {Journal of Open Source Software}
}
```

## Code of Conduct

Please note that the emba project is released with a [Contributor Code of Conduct](https://bblodfon.github.io/emba/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
# emba 0.1.8

- add option `all.ss` in function `get_stable_state_from_models_dir()` to include all stable states in the returned `tibble` object
- fix bug in `get_stable_state_from_models_dir()` (return empty tibble when no models had stable state)

# emba 0.1.7

- Input functions that read model directories with `.gitsbe` files, now disregard other kind of files that might be inside these directories.
- add minimum package dependencies in `DESCRIPTION` file
- add JOSS paper

# emba 0.1.6

- Fixed test for the `update_biomarker_files` function (writes to `tmpdir()` instead of the user's library directory)

# emba 0.1.5

- Finally added tests to the package! **Coverage is now 97%**.
- Used the [pkgdown](https://github.com/r-lib/pkgdown/) package to create static html documentation for emba. [Check it here](https://bblodfon.github.io/emba/index.html)!
- Change MCC calculation to return 0 when undefined/`NaN` MCC scores were produced (which is the correct limiting value - see [Chicco at al. (2020)](https://doi.org/10.1186/s12864-019-6413-7)). Thus, the previous versions handling of `NaN` MCC scores, is now deprecated.
- Add the `penalty` parameter to account for the difference in model group size when calculating the average activity or link operator data differences. This minimizes the bias in the returned biomarkers.
    - For the implementation check the function `emba::get_vector_diff` and the corresponding [StackOverflow question](https://math.stackexchange.com/questions/3547139/formula-for-weighted-average-difference).
    - To get the same results as with previous versions of this library, use `penalty=0` in the general `emba::biomarker_*` functions (though the results will probably be very biased and that's why the default value for the `penalty` is now **0.1**).
- Changed documentation to specify that the `models.stable.state` parameter used in various functions can take any values in the [0,1] interval and not just 0 (*inactive*) and 1 (*active*).
- The following functions do not take the redundant parameter `models` anymore:
  - `emba::get_avg_link_operator_diff_mat_based_on_tp_predictions`
  - `emba::get_avg_activity_diff_mat_based_on_tp_predictions`
  - `emba::get_avg_activity_diff_based_on_tp_predictions`
- Refactor several of the functions that load the results from the DrugLogics pipeline:
  - If a model has less or more than 1 stable state, it's discarded and a message is printed.
  - Return value is now a `data.frame` object instead of a `matrix`.
  - The models names do not have the annoying `.gitsbe` extension anymore.
  - These changes affect the following functions: `emba::get_link_operators_from_models_dir`, `emba::get_stable_state_from_models_dir` and `emba::get_model_names`.
- The general functions `emba::biomarker_mcc_analysis` and `emba::biomarker_tp_analysis` do not use the `calculate.subsets.stats` input option anymore.
The `emba::biomarker_synergy_analysis` continues to do so and now also calculates and returns all possible synergy set and subset pairs that miss just one of the model predicted synergies (`emba::get_synergy_comparison_sets`).
- Various small bug fixes and other code refactoring :)

# emba 0.1.4

- `get_synergy_scores` now supports reading both *ensemble-wise* and *model-wise* synergies files
- add `calculate.subsets.stats` option to the general analysis functions (`biomarker_*`) that decides if the powerset of the observed synergies and the number of models predicting each subset is going to be calculated. 
The default value is set to `FALSE` to save computation time :)

# emba 0.1.3

- add function `get_synergy_scores`
- fixed test that used a randomly generated matrix

# emba 0.1.2

- add function `get_avg_link_operator_diff_based_on_synergy_set_cmp`
- add function `get_avg_link_operator_diff_based_on_specific_synergy_prediction`
- add function `filter_network` - to use for visualizing induced subgraphs
- update dependencies (set `usefun` min version to 0.4.3)

# emba 0.1.1

- Optimized `count_models_that_predict_synergies` function and added tests for it. For a benchmark see
relative [Stack Overflow thread](https://stackoverflow.com/questions/58380043/optimize-r-code-for-row-operations-on-ternary-data-frame)

# emba 0.1.0

- Added a `NEWS.md` file to track changes to the package
- Transferred functions from separate R scripts to the package
- Finished code refactoring and splitting to different modules
- Finished writing documentation for all functions
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
reported to the community leaders responsible for enforcement at this [e-mail](bblodfonbblodfon@gmail.com). 
All complaints will be reviewed and investigated promptly and fairly.

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

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 2.0, available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at https://www.contributor-covenant.org/faq.
Translations are available at https://www.contributor-covenant.org/translations.
---
title: 'emba: R package for analysis and visualization of biomarkers in boolean model ensembles'
tags:
  - R
  - boolean networks
  - logical modeling
  - biomarkers
  - mechanistic models
  - drug synergies
  - anti-cancer drug combinations
  - druglogics
authors:
  - name: John Zobolas
    orcid: 0000-0002-3609-8674
    affiliation: 1, 2
  - name: Martin Kuiper
    orcid: 0000-0002-1171-9876
    affiliation: 1
  - name: Åsmund Flobak
    orcid: 0000-0002-3357-425X
    affiliation: 2, 3
affiliations:
 - name: Department of Biology, Norwegian University of Science and Technology (NTNU), Trondheim, Norway
   index: 1
 - name: Department of Clinical and Molecular Medicine, Norwegian University of Science and Technology (NTNU), Trondheim, Norway
   index: 2
 - name: The Cancer Clinic, St. Olav’s Hospital, Trondheim, Norway
   index: 3
date: 31 July 2020
bibliography: paper.bib
---

# Introduction

Computational modeling of cellular systems has been one of the most powerful tools used to build interpretable knowledge of biological processes and help identify molecular mechanisms that drive diseases such as cancer [@Aldridge2006].
In particular, the use of logical modeling has proven to be a substantially useful approach, since it allows the easy construction, simulation and analysis of predictive models, capable of providing a qualitative and insightful view on the extremely complex landscape of biological systems [@Morris2010; @Wang2012; @Abou-Jaoude2016].
Such mechanistic models, with the systematic integration of prior knowledge and experimental data, have been extensively used to better understand what drives deregulation of signal transduction, the outcome of which is the manifestation of diseases [@Traynard2017].
Furthermore, their explanatory power has been used to provide insights into a drug’s mode of action, investigate the mechanisms of resistance to drugs [@Eduati2017] and suggest new therapeutic combination candidates, among others [@Flobak2015].

One of the major challenges in systems medicine, has been the identification of scientifically validated, predictive biomarkers that correlate with patient response to given therapies.
The analysis of biological predictive markers of pharmacologic response can not only further our understanding of the systemic processes involved in diseases but can also help to classify patients into groups with similar responses to specific therapeutic interventions, advancing personalized medicine [@Senft2017].
In addition, the identification of biomarkers in tumor cells (e.g. mutations) has enabled the discovery of drug targets which are utilized in combinatorial molecular-targeted therapies - a strategy which aims to treat specific patient subgroups and has shown larger overall survival rates and reduced side-effects than monotherapy [@Al-Lazikani2012].
Despite the huge advancements towards drug combination therapy, genetic heterogeneity, drug resistance and drug combination synergy mechanisms still pose fundamental challenges to clinicians, modelers and lab researchers.

To help bridge the model simulation results with the (clinical) laboratory observations, several optimization methods have been used, such as model calibration, parameter estimation and sensitivity analysis.
These methods also allow us to determine which model parameters have the biggest influence in the overall behaviour of the system [@Aldridge2006].
For example, in @Frohlich2018, a computational framework that allowed for the efficient parameterization and contextualization of a large-scale cancer signaling network, was used to predict combination treatment outcome from single drug data.
This model was calibrated to fit and accurately describe specific cell-line experimental data, while enabling the identification of biomarkers of drug sensitivity as well as molecular mechanisms that affect drug resistance.
Furthermore, in @Dorier2016, a network optimization approach which topologically parameterized boolean models according to a genetic algorithm was used, in order to best match the experimentally observed behaviour.
This method resulted in an ensemble of boolean models which can be used to simulate response under drug perturbations in order to assess the underlying mechanisms and to generate new testable hypotheses.
Such an aggregation of best-fit models (wisdom of the crowds) has been shown to be quite robust and effective for model prediction performance [@Marbach2012].

# Statement of need

There is a plethora of software tools devoted to the qualitative modeling and analysis of biological networks.
The Consortium for the development of Logical Models and Tools (CoLoMoTo) is a community effort which aims to standardize the representation of logical networks and provide a common repository of methods and tools to analyze these networks [@Naldi2015].
Furthermore, to facilitate the access to several software logical modeling tools and enable reproducible computational workflows, the CoLoMoTo Interactive Notebook was introduced as a unified computational framework [@Naldi2018a].
The incorporated tools are accessed via a common programming interface (though originally implemented in different programming languages e.g. Java, Python, C\texttt{++} and R) and offer a collection of features like accessing online model repositories [@Helikar2012], model editing [@Naldi2018b], dynamical analysis (finding attractors, stochastic simulations, reachability properties, model-checking techniques) [@Mussel2010; @Klarner2016; @Stoll2017; @Pauleve2017; @Naldi2018] and model parameterization/optimization to fit perturbation signaling data [@Terfve2012; @Gjerga2020].
Despite the diverse and multi-purpose logical modeling tools that exist, there is still a lack of data analysis-oriented software that assists with the discovery of predictive biomarkers in ensembles of parameterized boolean networks that have been subject to drug combination perturbations.

The `emba` R package aims to fill that gap and provide a first implementation of such a novel software.
Initially, it was designed as a complementary software tool, to help the analysis of the parameterized boolean model ensembles which were produced by modules from the DrugLogics NTNU software pipeline (see respective documentation [@dl-doc]).
Later, we generalized most of the functions in the package and modularized them to package-essential (that form the core of the `emba` package) and various general-purpose yet useful functions (that are now part of the dependency package `usefun` [@R-usefun]).

# Summary

The main functionality of the `emba` R package is to find *performance* and *synergy* biomarkers.
Performance biomarkers are nodes in the input boolean networks whose activity state and/or model parameterization affects the predictive performance of those models.
The prediction performance can be assessed via the number of true positive predictions or the Matthews correlation coefficient score which is more robust to class imbalances [@Chicco2020].
On the other hand, synergy biomarkers are nodes which provide hints for the mechanisms behind the complex process of synergy manifestation in drug combination datasets.

For more information, see our “Get started guide” and the reference manual in the package website [@emba-site].
Several analyses using the `emba` R package are available in a separate repository [@gma-git].
Future developments will include the implementation of a method for the identification of *topology* biomarkers, where we will be able to assess which interactions in the network are important for the manifestation of synergies in specific cell-contexts.

# Acknowledgements

This work was supported by ERACoSysMed grant *COLOSYS* (JZ, MK) and The NTNU Strategic Research Area *NTNU Health* (AF).

# References
---
title: "emba"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{emba}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Intro

The *emba* R package name stands for *Ensemble (Boolean) Model Biomarker Analysis*.
It's main purpose is to be used on a dataset consisted of an **ensemble of boolean models**.
These models are usually (but not necessarily) different versions of the same initial model, parameterized in different ways (e.g. some boolean operators in the model equations have changed from *OR* to *AND* or vice-versa).
A prerequisite for using this package, is that this model dataset must be tested in-silico (using some computational way) against a list of drug combinations, in order to assess which drugs combinations behave *synergistically* for which models.
An example software that generates such boolean model ensembles and performs a comprehensive drug response analysis on them is the DrugLogics NTNU software pipeline (see respective [documentation](https://druglogics.github.io/druglogics-doc/)).

Given a list of *gold-standard* (lab-observed/verified) **synergies** ^[Note that the assessment of these synergies based on experimental data (usually High-Throughput Screening data) is an analysis on its own], this package enables the easy grouping of the models into different classes based on a specific performance metric evaluation.
This model classification enables the discovery and visualization of **biomarkers** - nodes whose *activity* and/or boolean model *parameterization* might affect either the prediction performance of those models or the manifestation of the predicted synergies.

In the next sections we will describe the main inputs and outputs of the *general analysis* functions (which group a lot of functionality into one) and provide some insights on the implementation behind.
Biomarkers will be assessed and visualized using a test dataset generated from the DrugLogics software mentioned above.

The complementary R package [usefun](https://github.com/bblodfon/usefun) has various helpful functions that are used both inside the *emba* package and during the analysis below.

For further analyses using this package on boolean model ensemble datasets see this [GitHub repository](https://github.com/druglogics/gitsbe-model-analysis/).
See also [an example](https://druglogics.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/) that demonstrates all the intermediate steps included in the *general analysis* functions as well as other miscellaneous usages that this guide does not cover.
Lastly, you might also want to check a [nice presentation](https://bblodfon.github.io/r-pres/digital_life_2019.html#12) I made for a conference about this package.

# Setup

```{r setup, message=FALSE}
# libraries
library(emba)
library(usefun)
library(dplyr)
library(knitr)
library(Ckmeans.1d.dp)

# wrapper to change printing to invisible
pr = function(x) invisible(x)
```

# Input

## Test dataset

The test dataset we will use has $7500$ boolean models with $139$ nodes each. 
It helps to think of each boolean model as a **network of nodes** where the edges represent either activation or inhibition of the corresponding target and the nodes activity can be either active (1) or inactive (0).

The models have been assessed for synergy against a total of $153$ drug combinations.

```{r input-1}
data.list = readRDS(url("https://github.com/bblodfon/emba/blob/main/vignettes/data.rds?raw=true"))

model.predictions = data.list$model.predictions
models.stable.state = data.list$models.stable.state
models.link.operator = data.list$models.equations
observed.synergies = data.list$observed.synergies

# (x,y) coordinates for visualization
nice.layout = data.list$nice.layout
# model network as an igraph object
net = data.list$net

# drug combinations
drug.combos = colnames(model.predictions)

# change model names (shorter names for readability)
model.names = paste0("model", 1:7500)
rownames(model.predictions) = model.names
rownames(models.stable.state) = model.names
rownames(models.link.operator) = model.names
```

## Model Predictions

This data represents the results of **in-silico testing the boolean models against a drug combination dataset**.
More specifically, the model predictions is a `data.frame` whose values (corresponding to a specific **model-drug combination element**) can be one of the following:

- 0 (no synergy predicted)
- 1 (synergy was predicted)
- `NA` (in case the model couldn't be assessed for synergy, e.g. there were no stable states in either the drug combination perturbed model or in any of the two single-drug perturbed models).

```{r input-2}
model.predictions[1:5, 77:84] %>% kable(caption = "Model predictions example")
```

## Model stable states

Each **model** must have a *stable state configuration* where the **nodes** have fixed to either 0 (inactive state) or 1 (active state).
In other words, a **fixpoint attractor**.
Of course, if a model has multiple attractors or other methods are used to derive a solution to the system of boolean equations that is the model itself, then *continuous* activity state values (in the $[0,1]$ interval) are also supported.

```{r input-3}
models.stable.state[1:5, 5:11] %>% kable(caption = "Model stable states example")
```

## Model link operators

This is a non-essential input for the functions we will use, but we include it here since the test dataset supports it.
It is a way to represent the **structure (parameterization)** of the boolean models in the dataset.

If each boolean model is a list of boolean equations of the form:

`T = (A1 OR A2 OR ...) AND NOT (I1 OR I2 OR ...)`

, where the `A` and `I` nodes are the activating and inhibiting regulators respectively of the **target node** `T` and the `AND NOT` is the **link (balance) operator**, we can specify a `data.frame` object whose values (corresponding to a specific **model-target node element**) can be one of the following:

- 0 (`AND NOT` link operator)
- 1 (`OR NOT` link operator) 
- 0.5 (if the target node does not have *both* activating and inhibiting regulators and thus the corresponding boolean equation has **no link operator**)

```{r input-4}
models.link.operator[1:5, 1:10] %>% kable(caption = "Models link operator example")
```

Note that in the test dataset, the nodes (columns of the `models.link.operator` object) who didn't have a link operator are pruned.

## Observed (GS) synergies

A list of *gold standard (GS)* drug combinations which have been termed as **synergistic** via experimental and/or other computational methods.
These drug combinations must be a subset of the ones tested in the models (the column names of the `model.predictions` data).

```{r input-5, results='asis'}
usefun::pretty_print_vector_values(observed.synergies, vector.values.str = "observed synergies")
```

# Performance biomarkers

*Performance biomarkers* are nodes in our studied networks (boolean models) whose activity state and/or boolean model parameterization (link operator) affects the prediction performance of those models.
These nodes can be thus used as indicators of either *activity* or *structural* changes that have a positive effect on the prediction performance of our models.

The model performance can be assessed via various ways. 
In this package we offer two ways to group the models to different classification categories: either based on the 
**number of true positive (TP) predictions** or on the **Matthews correlation coefficient (MCC) score** with respect to the drug combination dataset tested for synergy.
The function `emba::biomarker_tp_analysis()` is used for the former classification and the function `emba::biomarker_mcc_analysis()` for the latter.
Note that it's generally better to use the MCC classification, since it's a more robust performance evaluation metric compared to the number of TP predictions, since it takes into account all of the four confusion matrix values.

When the models have been grouped to different classification categories, their nodes activity or boolean model parameterization can be summarised in each group and compared to the others, obtaining thus the expected biomarkers using the methodology described below.

## TP-based analysis

We use the `emba::biomarker_tp_analysis()` function with the specified inputs:
```{r tp-analysis-1}
tp.analysis.res = emba::biomarker_tp_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  penalty = 0.1,
  threshold = 0.55)
```

The `penalty` term is used to reduce the bias when model groups have different *sizes*.
For example, if I were to compare the average activity of nodes between two groups of models, with respective group sizes 5 and 1000, then the result would be heavily biased towards the group with the larger size, making thus the quality of the results coming out of this comparison questionable.
As such, with `penalty` values closer to 0, more bias is introduced and we expect more biomarkers to be found.
The default value of $0.1$ is a good rule-of-thumb choice for minimizing such biases.
See more info on `emba::get_vector_diff()`.

---

As a first result, we get the predicted synergies - i.e. the drug combinations that are a subset of the observed ones and were **predicted by at least one** of the models in the dataset:
```{r tp-analysis-2, results = 'asis'}
usefun::pretty_print_vector_values(tp.analysis.res$predicted.synergies, vector.values.str = "predicted synergies")
```

The percentage of true positive predicted synergies is thus `r round(100*length(tp.analysis.res$predicted.synergies)/length(observed.synergies), digits = 1)`%.
Such a low number might be a sign that the models quality is poor (need for a different parameterization) or other reasons like incorrect assessment of the gold standard synergies, etc.

The next informative barplot shows the distribution of models according to their true positive predictions:
```{r tp-analysis-3, fig.align='center', fig.width=7, fig.height=5.6}
pr(emba::make_barplot_on_models_stats(table(tp.analysis.res$models.synergies.tp), 
  title = "True Positive Synergy Predictions",
  xlab = "Number of maximum correctly predicted synergies",
  ylab = "Number of models"))
```

- The **maximum number of predicted synergies by any individual model** is 3
- There are only 2 models in total that could predict these 3 synergies
- Almost half of the models make no true positive predictions
- This model classification is largely *skewed*

Next result we get is the **average activity differences** per network node for all group classifications:
```{r tp-analysis-4}
tp.analysis.res$diff.state.tp.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","STAT3","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Activity Difference Matrix")
```

- Rows represent the different classification group matchings, e.g. (1,2) means the models that predicted 1 TP synergy vs the models that predicted 2 TP synergies.
- All values are in the $[-1,1]$ interval. 
The more negative the activity difference value, the more *inhibited* the node is in the better performance models (e.g. `STAT3` node).
The more positive the activity difference value, the more *active* the node is in the better performance models (e.g. `CEBPA` node).
- Based on a **user-given** `threshold` level, a node is declared as an **activity biomarker** if it's highest absolute value surpasses that threshold (see `emba::get_biomarkers()` for more info).

In our case, `threshold = 0.55` and thus `CEBPA` and `PSEN1` are returned as active biomarkers:
```{r tp-analysis-5, results='asis'}
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.active,
  vector.values.str = "active biomarkers")
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.inhibited,
  vector.values.str = "inhibited biomarkers")
```

With the models initial network as an [igraph](https://igraph.org/r/) object (see `emba::construct_network()` on how to create such a `net` object), we can visualize every row of the above matrix as follows:
```{r tp-analysis-6, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::plot_avg_state_diff_graph(net, tp.analysis.res$diff.state.tp.mat["(2,3)",], 
  layout = nice.layout, title = "Bad models (2 TP) vs Good models (3 TP)"))
```

Note that with less `penalty`, more bias would be introduced and thus more biomarkers would be found (even for a higher chosen `threshold`):
```{r tp-analysis-7, results='asis'}
tp.analysis.res.biased = emba::biomarker_tp_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  penalty = 0,
  threshold = 0.7)

usefun::pretty_print_vector_values(tp.analysis.res.biased$biomarkers.tp.active,
  vector.values.str = "active biomarkers")

usefun::pretty_print_vector_values(tp.analysis.res.biased$biomarkers.tp.inhibited,
  vector.values.str = "inhibited biomarkers")
```

Last result we get is the **average link operator differences** per network node (whose boolean equation had a link operator) for all group classifications:
```{r tp-analysis-8}
tp.analysis.res$diff.link.tp.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","STAT3","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Link Operator Difference Matrix")
```

- Rows again represent the different classification group matchings, e.g. (1,2) means the models that predicted 1 TP synergy vs the models that predicted 2 TP synergies.
- All values are in the $[-1,1]$ interval.
A value closer to $-1$ means that on average, the node's boolean equation has the *AND NOT* link operator in the better performance models (e.g. `STAT3` node).
A value closer to $1$ means that on average, the node's boolean equation has mostly the *OR NOT* link operator in the better performance models (e.g. `CEBPA` node).
- Based on the given `threshold` level, a node is declared as a **link operator biomarker** if it's highest absolute value surpasses that threshold (see `emba::get_biomarkers()` for more info).

In our case, `threshold = 0.55` and thus `CEBPA` is returned as an `OR` link operator biomarker:
```{r tp-analysis-9, results='asis'}
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.or,
  vector.values.str = "'OR' biomarkers")
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.and,
  vector.values.str = "'AND' biomarkers")
```

We can also visualize every row of the average link operator differences matrix as follows:
```{r tp-analysis-10, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::plot_avg_link_operator_diff_graph(net, tp.analysis.res$diff.link.tp.mat["(2,3)",], 
  layout = nice.layout, title = "Bad models (2 TP) vs Good models (3 TP)"))
```

Interpreting the result regarding the `CEBPA` biomarker, we look back at its boolean equation and we see that the higher performance models must have the `OR NOT` link operator in order for `CEBPA` to be in an active (ON) state (an `AND NOT` results mostly on an inhibited state for `CEBPA`):

`CEBPA = (GSK3B OR MAP2K1 OR MEK1/2) OR NOT CTNNB1`

## MCC-based analysis

We use the `emba::biomarker_mcc_analysis()` function with the specified inputs:
```{r mcc-analysis-1}
mcc.analysis.res = emba::biomarker_mcc_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  threshold = 0.65,
  num.of.mcc.classes = 4,
  penalty = 0.2)
```

- The `penalty` term is used to reduce the bias when model groups have different *sizes* (default value is $0.1$).
See more info about this on the TP-based analysis above and on the documentation of the function `emba::get_vector_diff()`.
- We can choose the number of model groups to be created (`num.of.mcc.classes` parameter, with default value $5$).
Internally, the function `Ckmeans.1d.dp()` is used to perform an optimal univariate *k-means* clustering on the models MCC scores, i.e. **it groups the models to different MCC classes**, with higher classes having higher MCC scores (corresponding thus to better performance models).

---

First result is the predicted synergies, which are the same as the ones found with the TP-based analysis (the model predictions did not change).
As such, the drug combinations which were **predicted by at least one of the models** in the dataset are:
```{r mcc-analysis-2, results = 'asis'}
usefun::pretty_print_vector_values(mcc.analysis.res$predicted.synergies, vector.values.str = "predicted synergies")
```

We can get a first idea of the range and distribution of the models MCC scores with the next barplot:
```{r mcc-analysis-3, fig.align='center', fig.width=7, fig.height=5.8}
pr(emba::make_barplot_on_models_stats(table(mcc.analysis.res$models.mcc), 
  title = "MCC scores", xlab = "MCC value", 
  ylab = "Number of models", cont.values = TRUE))
```

- There are no relatively bad models (MCC values close to $-1$)
- Most of the models perform a little better than random prediction ($MCC > 0$)

We can also plot the **MCC-model histogram**, which in addition shows the estimated *density* (how many models) and *width* (MCC range) of each MCC class:
```{r mcc-analysis-4, fig.align='center', fig.width=7, fig.height=5.5}
models.mcc = mcc.analysis.res$models.mcc
num.of.mcc.classes = 4

res = Ckmeans.1d.dp(x = models.mcc, k = num.of.mcc.classes)
models.cluster.ids = res$cluster

pr(emba::plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes))
```

Next result we get is the **average activity differences** per network node for all group classifications:
```{r mcc-analysis-5}
mcc.analysis.res$diff.state.mcc.mat %>%
  as.data.frame() %>%
  select(c("AKT","PPM1A","PTEN","PSEN1","PTK2","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Activity Difference Matrix")
```

- Rows represent the different classification group matchings, e.g. (1,2) means the models that were in the 1st MCC class vs the models that were in the 2nd MCC class.
- All values are in the $[-1,1]$ interval. 
The more negative the activity difference value, the more inhibited the node is in the better performance models (e.g. `PTK2`). 
The more positive the activity difference value, the more active the node is in the better performance models (e.g. `PPM1A,PTEN`).
- Based on a **user-given** `threshold` level, a node is declared as an **activity biomarker** if it’s highest absolute value surpasses that threshold (see `emba::get_biomarkers()` for more info).

In our case, `threshold = 0.65` and thus `PTEN` and `PPM1A` are returned as active biomarkers and `PTK2` as an inhibited biomarker:
```{r mcc-analysis-6, results='asis'}
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.active,
  vector.values.str = "active biomarkers")
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.inhibited,
  vector.values.str = "inhibited biomarkers")
```

Note that looking at the respective boolean equations:

- `PPM1A = PTEN`
- `PTK2 = not PTEN`

we conclude that **the only activity biomarker of interest** is `PTEN` as it's the only regulator whose state directly influences the `PPM1A` and `PTK2` nodes.

With the models initial network as an [igraph](https://igraph.org/r/) object (see `emba::construct_network()` on how to create such a `net` object), we can visualize every row of the above matrix as follows:
```{r mcc-analysis-7, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::plot_avg_state_diff_graph(net, mcc.analysis.res$diff.state.mcc.mat["(1,4)",], 
  layout = nice.layout, title = "Bad models (MCC Class 1) vs Good models (MCC Class 4)"))
```

Last result we get is the **average link operator differences** per network node (whose boolean equation had a link operator) for all group classifications:
```{r mcc-analysis-8}
mcc.analysis.res$diff.link.mcc.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","CEBPA","STAT3","JAK1")) %>% # show only part of the matrix
  kable(caption = "Average Link Operator Difference Matrix")
```

- Rows again represent the different classification group matchings, e.g. (1,2) means the models that were in the 1st MCC class vs the models that were in the 2nd MCC class.
- All values are in the $[-1,1]$ interval.
A value closer to $-1$ means that on average, the node's boolean equation has the *AND NOT* link operator in the better performance models (e.g. `STAT3` node).
A value closer to $1$ means that on average, the node's boolean equation has mostly the *OR NOT* link operator in the better performance models (e.g. `PTEN` node).
- Based on the given `threshold` level, a node is declared as a **link operator biomarker** if it's highest absolute value surpasses that threshold (see `emba::get_biomarkers()` for more info).

In our case, `threshold = 0.65` and thus `PTEN` is returned as an `OR` link operator biomarker:
```{r mcc-analysis-9, results='asis'}
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.or,
  vector.values.str = "'OR' biomarkers")
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.and,
  vector.values.str = "'AND' biomarkers")
```

We can also visualize every row of the average link operator differences matrix as in the following example:
```{r mcc-analysis-10, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::plot_avg_link_operator_diff_graph(net, mcc.analysis.res$diff.link.mcc.mat["(1,4)",], 
  layout = nice.layout, title = "Bad models (MCC Class 1) vs Good models (MCC Class 4)"))
```

## Comparing the 2 methods

Overall, we note that **using the more robust MCC score** to classify the models according to their prediction performance on the drug combination dataset they were tested on, **produces more reliable biomarkers** compared to using the simple number of true positive predictions.
In addition, the biomarker results are different between the 2 methods, e.g. the TP-analysis revealed `CEBPA` as an active state performance biomarker whereas the MCC-based analysis showed `PTEN` to be so.

# Synergy biomarkers

*Synergy* biomarkers are nodes in our studied networks (boolean models) whose activity state and/or boolean model parameterization (link operator) affects the manifestation of synergies.
These nodes can be thus used as indicators of either *activity* or *structural* changes that make the models **predict specific drug combinations as synergistic**.

The core idea behind the implementation is that the models are now classified to groups based on whether they **predict or not each one of the predicted synergies** (which for the test dataset are the same 5 as found with the previous analyses).
Thus, by comparing the average node activity or boolean model parameterization from the group that predicted a drug combination as a **synergy** vs the group that predicted it to be an **antagonism**, we can derive biomarkers for that drug combination.

The function used to perform such an analysis is the `emba::biomarker_tp_analysis()`:

```{r synergy-analysis-1}
synergy.analysis.res = emba::biomarker_synergy_analysis(
  model.predictions,
  models.stable.state,
  models.link.operator,
  observed.synergies,
  threshold = 0.5,
  calculate.subsets.stats = TRUE,
  penalty = 0.1)
```

---

Now in addition to the predicted synergies set, we get all the **subsets** for which a model predicted all drug combinations in that subset as *synergistic*.
We can visualize this result with the `emba::make_barplot_on_synergy_subset_stats()` function:
```{r synergy-analysis-2, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::make_barplot_on_synergy_subset_stats(
  synergy.analysis.res$synergy.subset.stats,
  threshold.for.subset.removal = 1, 
  bottom.mar = 9))
```

- Almost half of the models ($3478$) predict none of the *gold standard* synergies
- The `PI-D1` synergy is predicted by almost all the rest of the models
- There were only $2$ models that predicted $3$ gold standard synergies, i.e. the set `BI-PD,PD-PI,PI-D1` which is the maximum number of predicted synergies by an individual model

Next result is the matrix of *activity state differences* vectors, one for each of the predicted synergies:
```{r synergy-analysis-3, fig.align='center', fig.width=7, fig.height=5}
synergy.analysis.res$diff.state.synergies.mat[,1:8] %>% 
  kable(caption = "Average State Differences per Synergy Predicted", digits = 3)
```

- Every row refers to a different predicted synergy in the dataset. 
The columns are the network nodes.
- Every row is the result of **comparing two model groups**: the *synergistic group* (models that predicted the row-annotated drug combination as **synergistic**) vs the *antagonistic group* (models that predicted the row-annotated drug combination as **antagonistic**).
- All values are in the $[−1,1]$ interval.
The more negative (positive) the activity difference value, the more inhibited (active) the node is in the synergistic model group.
- Based on a **user-given** `threshold` level, a node for a specific synergy is declared as an **activity biomarker** if it’s highest absolute value surpasses that threshold (see `emba::get_biomarkers()` for more info).

Every row of the above matrix can be also network-plotted.
We show for example the average state difference graph for the `PI-D1` synergy:
```{r synergy-analysis-4, fig.align='center', fig.width=7, fig.height=5.5}
pr(emba::plot_avg_state_diff_graph(net,
  synergy.analysis.res$diff.state.synergies.mat["PI-D1",],
  layout = nice.layout, title = "Prediction of PI-D1 (Good Models: synergy, Bad Models: antagonism)"))
```

Given the **user-defined** `threshold` ($0.5$) we also get as a result the activity state biomarkers:
```{r synergy-analysis-5, fig.align='center', fig.width=7, fig.height=5}
# prune nodes (columns) that were not found as biomarkers for any predicted synergy
biomarker.act.mat = usefun::prune_columns_from_df(
  df = synergy.analysis.res$activity.biomarkers, value = 0)

biomarker.act.mat[, 4:12] %>% # show only part of the matrix
  kable(caption = "Activity State Biomarkers Per Synergy Predicted")
```

- In the above matrix, $1$ means *active state biomarker*, $-1$ *inhibited state biomarker* and $0$ means not a biomarker.
- Using the following code, you can filter and derive your own activity biomarkers matrix:
```{r, eval = FALSE}
# define your own threshold
my.thres = 0.76
activity.biomarkers.new = as.data.frame(apply(
  synergy.analysis.res$diff.state.synergies.mat, c(1,2), 
  usefun::get_ternary_class_id, my.thres))
```

Note that there were predicted synergies (rows in the above matrix), for which we couldn't find activity biomarkers (row was all zeros).
This is justifiable, since **the number of models in the synergistic and antagonistic model groups can be fairly unbalanced** and the `penalty` term is used to correct this bias.
For example, comparing the models that predict `AK-PD` as synergistic vs those that predict it as antagonistic, we have:

```{r synergy-analysis-6, results='asis'}
drug.comb = "AK-PD"

syn.models.num = sum(model.predictions[, drug.comb] == 1 & !is.na(model.predictions[, drug.comb]))
ant.models.num  = sum(model.predictions[, drug.comb] == 0 & !is.na(model.predictions[, drug.comb]))

usefun::pretty_print_string(paste0("Number of models (AK-PD): #Synergistic: ", syn.models.num, ", #Antagonistic: ", ant.models.num))
```

Lastly, the `synergy.analysis.res$diff.link.synergies.mat` result is a matrix that contains the **average link operator differences** per network node (whose boolean equation had a link operator) when comparing the synergistic vs antagonistic model groups for each predicted synergy.
The corresponding **link operator biomarkers** (based on the given `threshold`) are given in the `synergy.analysis.res$link.operator.biomarkers` output.

# R Session Info

```{r r-session-info, comment=""}
sessionInfo()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_fitness_from_models_dir}
\alias{get_fitness_from_models_dir}
\title{Load the models fitness scores}
\usage{
get_fitness_from_models_dir(models.dir)
}
\arguments{
\item{models.dir}{string. A dir with \emph{.gitsbe} files/models}
}
\value{
a numeric vector with elements the fitness scores and the names of the
models included in the \emph{names} attribute.
}
\description{
Use this function to merge the fitness scores from all models into a single
vector (the fitness score is a value between 0 and 1 and denotes how close
was the model fitted to one or more training data observations). Each model's
fitness value is loaded from the respective \emph{.gitsbe} file that can be
found inside the given \code{models.dir} directory (other kind of files are
discarded).
}
\examples{

models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
models.fitness = get_fitness_from_models_dir(models.dir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_mat_based_on_mcc_clustering}
\alias{get_avg_activity_diff_mat_based_on_mcc_clustering}
\title{Get average activity difference matrix based on MCC clustering}
\usage{
get_avg_activity_diff_mat_based_on_mcc_clustering(
  models.mcc,
  models.stable.state,
  num.of.mcc.classes,
  penalty = 0
)
}
\arguments{
\item{models.mcc}{a numeric vector of Matthews Correlation Coefficient (MCC)
scores, one for each model. The \emph{names} attribute holds the models' names.
Can be the result of using the function \code{\link{calculate_models_mcc}}.}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{num.of.mcc.classes}{numeric. A positive integer larger than 2 that
signifies the number of mcc classes (groups) that we should split the models
MCC values.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of
average node activity state differences} between two groups of models where
the classification was based on the models' MCC values.
Rows represent the different classification group matchings, e.g. (1,2) means the
models that belonged to the 1st group of MCC values vs the models that
belonged to the 2nd group. The columns represent the network's node names.
Values are in the [-1,1] interval.
}
\description{
This function splits the Matthews correlation coefficient (MCC) scores
of the models to specific groups using the \pkg{Ckmeans.1d.dp}
package (groups are denoted by ids, e.g. 1,2,3, etc.
where a larger id corresponds to a group of models with higher MCC scores)
and for each pairwise
combination of group id matchings (e.g. (0,1), (1,3), etc.), it uses the
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}}
function, comparing thus all groups of models that belong to different
MCC classes while taking into account the given \code{penalty} factor and the
number of models in each respective model MCC group.
}
\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{get_x_axis_values}
\alias{get_x_axis_values}
\title{Get the refined x-axis values}
\usage{
get_x_axis_values(models.stats, there.is.one.NaN.category, cont.values)
}
\arguments{
\item{models.stats}{table object, the result of using \link[base]{table} on
a (numeric) vector. Usually it represents some models statistics summary -
counts for each TP prediction value for example.}

\item{there.is.one.NaN.category}{logical. Is there one \emph{NaN} category?
(check is done before on the \emph{names} attribute of the \code{models.stats})}

\item{cont.values}{logical. If TRUE, the values of the x-axis will be trimmed
to 3 digits after the decimal point. Otherwise, they will be returned as they
are.}
}
\description{
This function returns the x-axis values that are going to be used by
\link[emba]{make_barplot_on_models_stats} to render the bar plot.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{print_biomarkers_per_predicted_synergy}
\alias{print_biomarkers_per_predicted_synergy}
\title{Print biomarkers for each predicted synergy}
\usage{
print_biomarkers_per_predicted_synergy(
  biomarkers.dir,
  predicted.synergies,
  html.output = TRUE
)
}
\arguments{
\item{biomarkers.dir}{string. It specifies the full path name (without the
ending character \emph{/}) of the directory which holds the biomarker files
for each drug combination in the \code{predicted.synergies}.
The biomarker files must be formatted as: \emph{\%drug.comb\%_biomarkers_active} or
\emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is an element of
the \code{predicted.synergies} vector. If the files are not properly formatted
or don't even exist, zero biomarkers are reported.}

\item{predicted.synergies}{a character vector of the synergies (drug
combination names) that were predicted by \strong{at least one} of the models
in the dataset.}

\item{html.output}{logical. If TRUE, it makes the printed output nice for
an HTML document. Default value: TRUE.}
}
\description{
Print biomarkers for each predicted synergy
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_observed_synergies_per_cell_line}
\alias{get_observed_synergies_per_cell_line}
\title{Get observed synergies per cell line}
\usage{
get_observed_synergies_per_cell_line(cell.line.dirs, drug.combos)
}
\arguments{
\item{cell.line.dirs}{a character vector of the cell line directories, in the
form of \emph{\{path\}/cell_line_name}. The cell line name directory
should be different for each element of the vector as we use it to fill in the
\code{rownames} of the result \code{data.frame} object. Inside each cell line directory
we read the observed synergies from a file called \emph{observed_synergies}
(if it exists and is non-empty). This file has the names of the observed
drug combinations, one in each line.}

\item{drug.combos}{a character vector with elements the names of all the drug combinations
that were tested in the analysis.}
}
\value{
a data.frame, whose columns represent the drug combinations tested
and the rows the cell lines. Possible values for each \emph{cell line-drug combination}
element are either \emph{1} (an observed synergy) or \emph{0} (non-observed synergy).
}
\description{
Use this function to get the observed synergies from the respective
files inside the given list of cell line directories.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}
\alias{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}
\title{Get average activity difference matrix based on specific synergy prediction}
\usage{
get_avg_activity_diff_mat_based_on_specific_synergy_prediction(
  model.predictions,
  models.stable.state,
  predicted.synergies,
  penalty = 0
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{predicted.synergies}{a character vector of the synergies (drug
combination names) that were predicted by \strong{at least one} of the models
in the dataset. It must be a subset of the column names (the drug combinations)
of the \code{model.predictions} object.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of
average node activity state differences} between two groups of models where
the classification for each individual row was based on the prediction or not
of a specific synergistic drug combination.
The row names are the predicted synergies, one per row, while the columns
represent the network's node names. Values are in the [-1,1] interval.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}}
function on a vector of drug combinations that were observed as synergistic
(e.g. by experiments) but also found as such by at least one of the models
(these drug combinations are the \code{predicted.synergies}).
}
\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_based_on_tp_predictions}
\alias{get_avg_activity_diff_based_on_tp_predictions}
\title{Get the average activity difference based on the number of true positives}
\usage{
get_avg_activity_diff_based_on_tp_predictions(
  models.synergies.tp,
  models.stable.state,
  num.low,
  num.high,
  penalty = 0
)
}
\arguments{
\item{models.synergies.tp}{an integer vector of TP values. The \emph{names}
attribute holds the models' names and must be a subset of the row names
of the \code{models.stable.state} parameter. Consider using the function
\code{\link{calculate_models_synergies_tp}}.}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{num.low}{integer. The number of true positives representing the 'bad'
model class.}

\item{num.high}{integer. The number of true positives representing the 'good'
model class. This number has to be strictly higher than \code{num.low}.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a numeric vector with values in the [-1,1] interval (minimum and maximum
possible average difference) and with the \emph{names} attribute representing the name
of the nodes.
}
\description{
This function splits the models to 'good' and 'bad' based on the number of true
positive predictions: \emph{num.high} TPs (good) vs \emph{num.low} TPs (bad).
Then, for each network node, it finds the node's average activity in each of
the two classes (a value in the [0,1] interval) and then subtracts the
'bad' average activity value from the good' one, taking into account the
given \code{penalty} factor and the number of models in each respective
model group.
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node is more \strong{inhibited} in the 'good' models compared to the
'bad' ones while a value closer to 1 means that the node is more \strong{activated}
in the 'good' models. A value closer to 0 indicates that the activity of that
node is \strong{not so much different} between the 'good' and 'bad' models and
so it won't not be a node of interest when searching for indicators of better
performance (higher number of true positives) in the good models.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{biomarker_mcc_analysis}
\alias{biomarker_mcc_analysis}
\title{Biomarker analysis based on MCC model classification}
\usage{
biomarker_mcc_analysis(
  model.predictions,
  models.stable.state,
  models.link.operator = NULL,
  observed.synergies,
  threshold,
  num.of.mcc.classes = 5,
  penalty = 0.1
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models).}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.
Note that the rows (models) have to be in the same order as in the \code{model.predictions}
parameter.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names (same order as in the \code{model.predictions}
parameter) whereas the column names specify
the network nodes (gene, proteins, etc.). Possible values for each
\emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
\emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
by both activating and inhibiting regulators (no link operator). Default value:
NULL (no analysis on the models parameterization regarding the mutation of the
boolean equation link operator will be done).}

\item{observed.synergies}{a character vector with elements the names of the
drug combinations that were found as synergistic. This should be a subset of
the tested drug combinations, that is the column names of the \code{model.predictions}
parameter.}

\item{threshold}{numeric. A number in the [0,1] interval, above which (or
below its negative value) a biomarker will be registered in the returned result.
Values closer to 1 translate to a more strict threshold and thus less
biomarkers are found.}

\item{num.of.mcc.classes}{numeric. A positive integer larger than 2 that
signifies the number of mcc classes (groups) that we should split the models
MCC values. Default value: 5.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.1.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a list with various elements:
\itemize{
  \item \code{predicted.synergies}: a character vector of the synergies (drug
  combination names) that were predicted by \strong{at least one} of the models
  in the dataset.
  \item \code{models.mcc}: a numeric vector of MCC scores, one for each model.
  Values are in the [-1,1] interval.
  \item \code{diff.state.mcc.mat}: a matrix whose rows are \strong{vectors of
  average node activity state differences} between two groups of models where
  the classification was based on the \emph{MCC score} of each model and was
  found using an optimal univariate k-means clustering method
  (\code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}).
  Rows represent the different classification group matchings, e.g. (1,2)
  means the models that were classified into the first MCC class vs the models
  that were classified in the 2nd class (higher is better). The columns
  represent the network's node names. Values are in the [-1,1] interval.
  \item \code{biomarkers.mcc.active}: a character vector whose elements are
  the names of the \emph{active state} biomarkers. These nodes appear more
  active in the better performance models.
  \item \code{biomarkers.mcc.inhibited}: a character vector whose elements are
  the names of the \emph{inhibited state} biomarkers. These nodes appear more
  inhibited in the better performance models.
  \item \code{diff.link.mcc.mat}: a matrix whose rows are \strong{vectors of
  average node link operator differences} between two groups of models where
  the classification was based on the \emph{MCC score} of each model and was
  found using an optimal univariate k-means clustering method
  (\code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}).
  Rows represent the different classification group matchings, e.g. (1,2)
  means the models that were classified into the first MCC class vs the models
  that were classified in the 2nd class (higher is better).
  The columns represent the network's node names. Values are in the [-1,1] interval.
  \item \code{biomarkers.mcc.or}: a character vector whose elements are
  the names of the \emph{OR} link operator biomarkers. These nodes have
  mostly the \emph{OR} link operator in their respective boolean equations
  in the better performance models.
  \item \code{biomarkers.mcc.and}: a character vector whose elements are
  the names of the \emph{AND} link operator biomarkers. These nodes have
  mostly the \emph{AND} link operator in their respective boolean equations
  in the better performance models.
}
}
\description{
Use this function to perform a full biomarker analysis on an ensemble boolean model
dataset where the model classification is based on the \emph{Matthews correlation
coefficient score (MCC)}. This analysis enables the discovery of \emph{performance
biomarkers}, nodes whose activity and/or boolean model parameterization (link
operator) affects the prediction performance of the models (as measured by
the MCC score).
}
\seealso{
Other general analysis functions: 
\code{\link{biomarker_synergy_analysis}()},
\code{\link{biomarker_tp_analysis}()}
}
\concept{general analysis functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_mat_based_on_tp_predictions}
\alias{get_avg_activity_diff_mat_based_on_tp_predictions}
\title{Get average activity difference matrix based on the number of true positives}
\usage{
get_avg_activity_diff_mat_based_on_tp_predictions(
  models.synergies.tp,
  models.stable.state,
  penalty = 0
)
}
\arguments{
\item{models.synergies.tp}{an integer vector of TP values. The \emph{names}
attribute must hold the models' names. Consider using the function
\code{\link{calculate_models_synergies_tp}}.}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names (same order as in the \code{models.synergies.tp}
parameter) whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of
average node activity state differences} between two groups of models where
the classification was based on the number of true positive predictions.
Rows represent the different classification group matchings, e.g. (1,2) means the
models that predicted 1 TP synergy vs the models that predicted 2 TP
synergies and the columns represent the network's node names.
Values are in the [-1,1] interval.
}
\description{
This function finds all the TP values of the models given (e.g. 0,1,2,3) and
generates every pairwise combination (e.g. the group matchings: (0,1), (1,3),
etc.). Then, it uses the \code{\link{get_avg_activity_diff_based_on_tp_predictions}}
function on each generated classification group matching, comparing thus all
groups of models with different true positive (TP) values, while taking into
account the given \code{penalty} factor and the number of models in each
respective model group.
}
\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{count_models_that_predict_synergies}
\alias{count_models_that_predict_synergies}
\title{Count models that predict a set of synergies}
\usage{
count_models_that_predict_synergies(drug.comb.vec, model.predictions)
}
\arguments{
\item{drug.comb.vec}{a character vector. Elements are (synergistic) drug
combinations, each one being a string in the form \emph{A-B} - no spaces
between the names and the hyphen '-')}

\item{model.predictions}{a \code{data.frame} object with rows as the models
and columns the drug combinations tested. Possible values for each
\emph{model-drug combination element} are either \emph{0} (no synergy
predicted), \emph{1} (synergy was predicted) or \emph{NA}}
}
\value{
the number of models that predict the given drug combination set
(have a value of 1 in the respective columns of the \code{model.predictions}
data.frame). If the given set is empty, we return the number of models that
predicted no synergies at all (after the \emph{NA} values are discarded, the
number of rows in the \code{model.predictions} data.frame that have only zero
values)
}
\description{
Use this function to find the number of models that predict a given set of
drug combinations (usually the ones found as synergies).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_mcc_classes_hist}
\alias{plot_mcc_classes_hist}
\title{Plot histogram of the MCC classes}
\usage{
plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes)
}
\arguments{
\item{models.mcc}{a numeric vector of Matthews
Correlation Coefficient (MCC) scores, one for each model.
The \emph{names} attribute may hold the models' names (but it is not required).}

\item{models.cluster.ids}{a numeric vector of cluster ids assigned to each
model. It can be the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
with input the models' MCC values (\code{models.mcc}) and the number of clusters
(\code{num.of.mcc.classes}).}

\item{num.of.mcc.classes}{numeric. A positive integer (>2) that signifies the
number of mcc classes (groups) that we should split the models MCC values.}
}
\description{
This function is a wrapper of the \code{\link[Ckmeans.1d.dp]{ahist}} function
for plotting nicely the distribution of the MCC models' values.
}
\examples{
models.mcc = c(-0.04, -0.17, 0.15, -0.24, -0.02 , 0.27, -0.42 , 0.38)
models.cluster.ids = c(2,2,3,1,2,3,1,3)
num.of.mcc.classes = 3
plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{biomarker_synergy_analysis}
\alias{biomarker_synergy_analysis}
\title{Biomarker analysis per synergy predicted}
\usage{
biomarker_synergy_analysis(
  model.predictions,
  models.stable.state,
  models.link.operator = NULL,
  observed.synergies,
  threshold,
  calculate.subsets.stats = FALSE,
  penalty = 0.1
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models).}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network
nodes (gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.
Note that the rows (models) have to be in the same order as in the \code{model.predictions}
parameter.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names (same order as in the \code{model.predictions}
parameter) whereas the column names specify
the network nodes (gene, proteins, etc.). Possible values for each
\emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
\emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
by both activating and inhibiting regulators (no link operator). Default value:
NULL (no analysis on the models parameterization regarding the mutation of the
boolean equation link operator will be done).}

\item{observed.synergies}{a character vector with elements the names of the
drug combinations that were found as synergistic. This should be a subset of
the tested drug combinations, that is the column names of the \code{model.predictions}
parameter.}

\item{threshold}{numeric. A number in the [0,1] interval, above which (or
below its negative value) a biomarker will be registered in the returned result.
Values closer to 1 translate to a more strict threshold and thus less
biomarkers are found.}

\item{calculate.subsets.stats}{logical. If \emph{TRUE}, then the results will
include a vector of integers, representing the number of models that predicted
every subset of the given \code{observed.synergies} (where at least one model
predicts every synergy in the subset). The default value is \emph{FALSE}, since
the powerset of the predicted \code{observed.synergies} can be very large to compute.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.1.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a list with various elements:
\itemize{
  \item \code{predicted.synergies}: a character vector of the synergies (drug
  combination names) that were predicted by \strong{at least one} of the models
  in the dataset.
  \item \code{synergy.subset.stats}: an integer vector with elements the number
  of models the predicted each \strong{observed synergy subset} if the
  \emph{calculate.subsets.stats} option is enabled.
  \item \code{synergy.comparison.sets}: a \code{data.frame} with pairs of
  \emph{(set, subset)} for each model-predicted synergy where each respective
  subset misses just one synergy from the larger set (present only if the
  \emph{calculate.subsets.stats} option is enabled). Can be used to refine
  the synergy biomarkers by comparing any two synergy sets with the functions
  \code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}} or
  \code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}}.
  \item \code{diff.state.synergies.mat}: a matrix whose rows are
  \strong{vectors of average node activity state differences} between two
  groups of models where the classification for each individual row was based
  on the prediction or not of a specific synergistic drug combination. The
  row names are the predicted synergies, one per row, while the columns
  represent the network's node names. Values are in the [-1,1] interval.
  \item \code{activity.biomarkers}: a \code{data.frame} object with rows
  the \code{predicted synergies} and columns the nodes (column names of the
  \code{models.stable.states} matrix). Possible values for each
  \emph{synergy-node} element are either \emph{1} (\emph{active state}
  biomarker), \emph{-1} (\emph{inhibited state} biomarker) or \emph{0} (not
  a biomarker) for the given \code{threshold} value.
  \item \code{diff.link.synergies.mat}: a matrix whose rows are
  \strong{vectors of average node link operator differences} between two
  groups of models where the classification for each individual row was
  based on the prediction or not of a specific synergistic drug combination.
  The row names are the predicted synergies, one per row, while the columns
  represent the network's node names. Values are in the [-1,1] interval.
  \item \code{link.operator.biomarkers}: a \code{data.frame} object with rows
  the \code{predicted synergies} and columns the nodes (column names of the
  \code{models.link.operator} matrix). Possible values for each
  \emph{synergy-node} element are either \emph{1} (\emph{OR} link operator
  biomarker), \emph{-1} (\emph{AND} link operator biomarker) or \emph{0} (not
  a biomarker) for the given \code{threshold} value.
}
}
\description{
Use this function to discover \emph{synergy biomarkers}, i.e. nodes whose
activity and/or boolean equation parameterization (link operator) affect the
manifestation of synergies in the boolean models. Models are classified to groups based on
whether they predict or not each of the predicted synergies.
}
\seealso{
Other general analysis functions: 
\code{\link{biomarker_mcc_analysis}()},
\code{\link{biomarker_tp_analysis}()}
}
\concept{general analysis functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_edges_from_topology_file}
\alias{get_edges_from_topology_file}
\title{Get the edges from a specified topology}
\usage{
get_edges_from_topology_file(topology.file)
}
\arguments{
\item{topology.file}{string. The name of the .sif file (can be a full path
name).}
}
\value{
a matrix with as many rows as in the .sif topology file (each row is
an edge) and 4 columns defining the source and target node name, the
regulation (activation or inhibition) and the color (green or red) of the
signed interaction.
}
\description{
Use this function to read a topology .sif file (either space or tab-delimited)
and get a matrix of network edges specifying the source and target name, the
regulation effect (activation or inhibition) and the color (green or red) of
each interaction.
}
\examples{
topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
edges = get_edges_from_topology_file(topology.file)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emba.R
\docType{package}
\name{emba}
\alias{emba}
\title{emba}
\description{
Analysis and visualization of an ensemble of
boolean models for biomarker discovery in cancer cell networks.
}
\details{
For a complete list of functions, use \code{library(help = "emba")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_model_names}
\alias{get_model_names}
\title{Get the model names}
\usage{
get_model_names(models.dir)
}
\arguments{
\item{models.dir}{string. A directory with \emph{.gitsbe} files/models
(non-gitsbe files are disregarded).}
}
\value{
a character vector of the model names, corresponding to the names
of the \emph{.gitsbe} files (the extension is pruned).
}
\description{
Get the model names
}
\examples{

models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
models = get_model_names(models.dir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_models_synergies_tp}
\alias{calculate_models_synergies_tp}
\title{Count the predictions of the observed synergies per model (TP)}
\usage{
calculate_models_synergies_tp(observed.model.predictions)
}
\arguments{
\item{observed.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{synergistic}
(\emph{positive results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}
}
\value{
an integer vector with elements the number of true positive predictions
per model. The model names are given in the \emph{names} attribute (same order
as in the \emph{rownames} attribute of the observed.model.predictions
\code{data.frame}).
}
\description{
Since the given \code{observed.model.predictions} data.frame has only the
positive results, this function returns the total number of 1's in each row.
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_mcc}()},
\code{\link{calculate_models_mcc}()},
\code{\link{calculate_models_synergies_fn}()},
\code{\link{calculate_models_synergies_fp}()},
\code{\link{calculate_models_synergies_tn}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{update_biomarker_files}
\alias{update_biomarker_files}
\title{Update biomarker files for a specific synergy}
\usage{
update_biomarker_files(
  biomarkers.dir,
  drug.comb,
  biomarkers.active.new,
  biomarkers.inhibited.new,
  method = "replace"
)
}
\arguments{
\item{biomarkers.dir}{string. It specifies the full path name of the
directory (without the ending character \emph{/}) which holds the biomarker
files for the synergistic drug combination
specified in the parameter \code{drug.comb}. The biomarker files must be
formatted as: \emph{\%drug.comb\%_biomarkers_active} or
\emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is the value
of the \code{drug.comb} parameter.}

\item{drug.comb}{string. The drug combination (e.g. "A-B") that will be used
to identify the related biomarker files.}

\item{biomarkers.active.new}{a numeric vector whose \emph{names} attribute
includes the node names of the (newly found) \emph{active biomarkers} for the specified
synergy. The values of the vector are the average activity difference of
each node, derived from a comparison between 2 different groups of models.}

\item{biomarkers.inhibited.new}{a numeric vector whose \emph{names} attribute
includes the node names of the (newly found) \emph{inhibited biomarkers} for the specified
synergy. The values of the vector are the average activity difference of
each node, derived from a comparison between 2 different groups of models.}

\item{method}{string. It specifies the method to use to update the biomarker
files when there are \emph{common nodes} between the 'old' and 'new' biomarkers:
\enumerate{
  \item \code{replace}(DEFAULT): we discard the 'old' biomarkers and keep
  only the 'new' ones
  \item \code{prune.to.common}: we keep only the common biomarkers
  \item \code{extend}: we add to the 'old' set of biomarkers the extra ones
  from the 'new' set that are not non-common to the 'old' ones, extending
  thus the 'old' biomarker set
}}
}
\description{
This function gets the (previously-found or 'old') synergy biomarkers from their
respective files and if any of these files are empty (no 'old' biomarkers
found) or non-existent, the 'new' biomarkers (given as input vector parameters) are
automatically saved. When the 'new' biomarkers \strong{share common nodes}
with the 'old' biomarkers, there exist 3
possible ways to combine the results, given by the \code{method} parameter.
If \strong{no common nodes} exist, no matter the \code{method} selected,
the 'new' biomarkers are added to the 'old' ones.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_link_operator_diff_based_on_synergy_set_cmp}
\alias{get_avg_link_operator_diff_based_on_synergy_set_cmp}
\title{Get the average link operator difference based on the comparison of two synergy sets}
\usage{
get_avg_link_operator_diff_based_on_synergy_set_cmp(
  synergy.set.str,
  synergy.subset.str,
  model.predictions,
  models.link.operator,
  penalty = 0
)
}
\arguments{
\item{synergy.set.str}{a string of drug combinations, comma-separated. The
number of the specified combinations must be larger than the ones defined
in the \code{synergy.subset.str} parameter. They also must be included in the
tested drug combinations, i.e. the columns of the \code{model.predictions}
parameter.}

\item{synergy.subset.str}{a string of drug combinations, comma-separated.
There must be at least one combination defined and all of them should also
be included in the \code{synergy.set.str} parameter.}

\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
(\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a numeric vector with values in the [-1,1] interval (minimum and
maximum possible average difference) and with the names attribute
representing the name of the nodes.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}}
which splits the models to 'good' and 'bad' based on the predictions
of two different synergy sets, one of them being a subset of the other.
The 'good' models are those that predict the \code{synergy.set.str}
(e.g. "A-B,A-C,B-C") while the 'bad' models are those that predict the
\code{synergy.subset.str} (e.g. "A-B,B-C"). Then, for each network node,
the function finds the node's average link operator value in each of the two classes
(a value in the [0,1] interval, 0 being \emph{AND NOT} and 1 being \emph{OR NOT})
and then subtracts the bad class average link operator value from the good one,
taking into account the given \code{penalty} factor and the number of models
in the 'good' and 'bad' class respectively.
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node's boolean equation has the \strong{AND NOT} link operator in the
models that predicted the extra synergy(-ies) that are included in the
\code{synergy.set.str} but not in the \code{synergy.subset.str}, whereas a
value closer to 1 means that the node's boolean equation has mostly the
\strong{OR NOT} link operator in these models.
These nodes are potential \strong{link operator biomarkers} because the
structure of their respective boolean equations (denoted by their link operator)
can influence the prediction performance of a model and make it predict the
extra synergy(-ies). A value closer to 0 indicates that the link operator in
the node's boolean equation is \strong{not so much different} between the
models that predicted the synergy set and those that predicted it's subset,
so it won't not be a node of interest when searching for potential link operator
biomarkers for the extra synergy(-ies).
A value exactly equal to 0 can also mean that this node didn't not have a link operator
in its boolean equation, again making it a non-important indicator of difference
in model performance.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_avg_link_operator_diff_graphs}
\alias{plot_avg_link_operator_diff_graphs}
\title{Plot the graphs from an average link operator differences matrix}
\usage{
plot_avg_link_operator_diff_graphs(net, diff.mat, layout = NULL)
}
\arguments{
\item{net}{igraph graph object}

\item{diff.mat}{a matrix whose rows are \strong{vectors of average node link
operator differences} between two groups of models based on some kind of
classification (e.g. number of TP predictions) and whose names are set in the \code{rownames}
attribute of the matrix (usually denoting the different classification
groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
that predicted 2 TP synergies, if the classification is done by number of TP
predictions). The columns represent the network's node names. Could be the
result of using the function \code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}}.}

\item{layout}{a (nx2) numeric matrix of x-y coordinates (2 columns) for each
of the nodes (n) in the \code{net} igraph object. If NULL, we use the default
layout provided by \code{\link[igraph]{layout_nicely}}.}
}
\description{
This function presents a convenient way to use the
\code{\link{plot_avg_link_operator_diff_graph}} function multiple times.
}
\seealso{
Other network plotting functions: 
\code{\link{plot_avg_link_operator_diff_graph}()},
\code{\link{plot_avg_state_diff_graph_vis}()},
\code{\link{plot_avg_state_diff_graphs}()},
\code{\link{plot_avg_state_diff_graph}()}
}
\concept{network plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_model_predictions}
\alias{get_model_predictions}
\title{Load the models predictions data}
\usage{
get_model_predictions(model.predictions.file)
}
\arguments{
\item{model.predictions.file}{a tab-delimited file (for the specific format
check the example below)}
}
\value{
a \code{data.frame} object with rows the models and columns the
drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)
}
\description{
Use this function to read a file that has the model predictions data
and output it to a \code{data.frame} object.
}
\examples{

model.predictions.file = system.file("extdata", "model_predictions",
  package = "emba", mustWork = TRUE)
model.predictions = get_model_predictions(model.predictions.file)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_mcc}
\alias{calculate_mcc}
\title{Calculate Matthews correlation coefficient vector}
\usage{
calculate_mcc(tp, tn, fp, fn)
}
\arguments{
\item{tp}{numeric vector of TPs}

\item{tn}{numeric vector of TNs}

\item{fp}{numeric vector of FPs}

\item{fn}{numeric vector of FNs}
}
\value{
a numeric vector of MCC values, each value being in the [-1,1]
interval. If any of the four sums of the MCC formula are zero, then we return
an MCC score of zero, which can be shown to be the correct limiting value (model
is no better than a random predictor, see Chicco et al. (2020),
\doi{10.1186/s12864-019-6413-7}).
}
\description{
Use this function to calculate the MCC scores given vectors of \emph{TP} (true
positives), \emph{FP} (false positives), \emph{TN} (true negatives) and \emph{FN}
(false negatives) values.
Note that the input vectors have to be of the same size and have one-to-one value
correspondence for the output MCC vector to make sense.
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_models_mcc}()},
\code{\link{calculate_models_synergies_fn}()},
\code{\link{calculate_models_synergies_fp}()},
\code{\link{calculate_models_synergies_tn}()},
\code{\link{calculate_models_synergies_tp}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_neighbors}
\alias{get_neighbors}
\title{Get neighbor nodes}
\usage{
get_neighbors(net, nodes)
}
\arguments{
\item{net}{igraph object}

\item{nodes}{character vector of node names}
}
\value{
a character vector of all the unique neighbors of the given
\code{nodes} in the \code{net} graph.
}
\description{
Given an igraph network object and vector of node names, this function
returns the set of unique neighbor nodes considering both ingoing and outgoing
edges (the closed neighbourhood node set).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_link_operator_diff_mat_based_on_tp_predictions}
\alias{get_avg_link_operator_diff_mat_based_on_tp_predictions}
\title{Get average link operator difference matrix based on the number of true positives}
\usage{
get_avg_link_operator_diff_mat_based_on_tp_predictions(
  models.synergies.tp,
  models.link.operator,
  penalty = 0
)
}
\arguments{
\item{models.synergies.tp}{an integer vector of TP values. The \emph{names}
attribute must hold the models' names. Consider using the function
\code{\link{calculate_models_synergies_tp}}.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names (same order as in the \code{models.synergies.tp}
parameter) whereas the column names specify the network nodes (gene, proteins, etc.).
Possible values for each \emph{model-node
element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
(\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of average node link operator
differences} between two groups of models based on some kind of classification
(e.g. number of TP predictions) and whose names are set in the \code{rownames}
attribute of the matrix (usually denoting the different classification
groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
that predicted 2 TP synergies, if the classification is done by number of TP
predictions). The columns represent the network's node names. Values are in
the [-1,1] interval.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}}
function with the parameter \code{models.link.operator} as input in the place of
\code{models.stable.state}, since the two matrices representing the two inputs
have the same data format (rows represent models, columns represent nodes,
and each value is a number in the [0,1] interval).
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node's boolean equation has the \strong{AND NOT} link operator in the
'good' models compared to the 'bad' ones while a value closer to 1 means that
the node's boolean equation has mostly the \strong{OR NOT} link operator
in the 'good' models. A value closer to 0 indicates that the link operator in
the node's boolean equation is \strong{not so much different} between the
'good' and 'bad' models and so it won't not be a node of interest when
searching for indicators of better performance (higher number of true positives)
in the parameterization of the good models (the boolean equations). A value
exactly equal to 0 can also mean that this node didn't not have a link operator
in its boolean equation, again making it a non-important indicator of difference
in model performance.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_models_synergies_fp}
\alias{calculate_models_synergies_fp}
\title{Count the predictions of the non-synergistic drug combinations per model (FP)}
\usage{
calculate_models_synergies_fp(unobserved.model.predictions)
}
\arguments{
\item{unobserved.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{non-synergistic}
(\emph{negative results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}
}
\value{
an integer vector with elements the number of false positive predictions
per model. The model names are given in the \emph{names} attribute (same order
as in the \emph{rownames} attribute of the unobserved.model.predictions
\code{data.frame}).
}
\description{
Since the given \code{unobserved.model.predictions} data.frame has only the
negative results, this function returns the total number of 1's in each row.
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_mcc}()},
\code{\link{calculate_models_mcc}()},
\code{\link{calculate_models_synergies_fn}()},
\code{\link{calculate_models_synergies_tn}()},
\code{\link{calculate_models_synergies_tp}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_models_based_on_mcc_class_id}
\alias{get_models_based_on_mcc_class_id}
\title{Get models based on the MCC class id}
\usage{
get_models_based_on_mcc_class_id(class.id, models.cluster.ids, models.mcc)
}
\arguments{
\item{class.id}{an integer specifying the class id.}

\item{models.cluster.ids}{a numeric vector of cluster ids assigned to each
model. It is the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
with input the sorted vector of the models' MCC values.}

\item{models.mcc}{a numeric sorted vector of Matthews Correlation Coefficient (MCC)
scores, one for each model.
The \emph{names} attribute holds the models' names.}
}
\value{
a character vector of model names
}
\description{
This helper function finds all the models that belong to a specific MCC
cluster, i.e. their MCC values belong to the same cluster id.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{make_barplot_on_synergy_subset_stats}
\alias{make_barplot_on_synergy_subset_stats}
\title{Bar plot of observed synergy subsets}
\usage{
make_barplot_on_synergy_subset_stats(
  synergy.subset.stats,
  threshold.for.subset.removal,
  bottom.margin,
  ylim.add = 0,
  cell.line = NULL
)
}
\arguments{
\item{synergy.subset.stats}{integer vector with values the amount of models
that predicted each synergy subset, defined as a comma-separated string of
drug combinations in the \emph{names} attribute of the vector}

\item{threshold.for.subset.removal}{integer. Use it to discard elements of
the \code{synergy.subset.stats} vector that are strictly less than the
specified threshold}

\item{bottom.margin}{integer used to vertically fit in the names of the drug
combinations in the x-axis (specified in inches). The best \code{bottom.margin}
value depends on the \emph{maximum size} of a synergy subset as defined in the
\code{names} attribute of the \code{synergy.subset.stats}.
Some rules of thumb are:
size = 1 => bottom.margin = 4,
size = 2 => bottom.margin = 6,
size = 3 => bottom.margin = 9,
size = 4 => bottom.margin = 12, etc.}

\item{ylim.add}{integer. Signifies the height to add to the upper
\code{ylim} parameter on the barplot, in addition to the maximum bar height
across the whole plot. Default value is 0.}

\item{cell.line}{string. The name of the cell line to be used in the title
of the produced plot. Default value: NULL (the cell line name will not be
added to the title).}
}
\description{
Use this function to easily make a barplot that shows the amount of models
that predicted each synergy subset out of the set of all observed synergies.
}
\examples{
synergy.subset.stats = c(1,4,3,2)
names(synergy.subset.stats) = c("A-B", "B-C", "C-A", "C-D")
make_barplot_on_synergy_subset_stats(synergy.subset.stats,
threshold.for.subset.removal = 0, bottom.margin = 4, ylim.add = 0.5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_observed_model_predictions}
\alias{get_observed_model_predictions}
\title{Subset the model predictions to the (true) observed synergies}
\usage{
get_observed_model_predictions(model.predictions, observed.synergies)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{observed.synergies}{a character vector with elements the names of the
drug combinations that were found as synergistic}
}
\value{
a \code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{synergistic}
(\emph{positive results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)
}
\description{
Subset the model predictions to the (true) observed synergies
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{get_node_colors}
\alias{get_node_colors}
\title{Get the node colors}
\usage{
get_node_colors(net, diff, col)
}
\arguments{
\item{net}{an igraph graph object with the node names defined in \code{V(net)$name}}

\item{diff}{numeric vector. Every value is in the [-1,1] interval and
represents the average activity difference of each node. The node names have
to be specified in the \emph{names} attribute of the given \code{diff} vector
and have to be the same as in \code{V(net)$name}.}

\item{col}{a character vector of colors to do the color interpolation in the
[-1,1] interval. Usually a two-element vector specifying the colors matching
the start and end of the interval (-1 and 1 respectively) or a three-element
vector specifying the colors matching the values -1, 0 and 1 (can be more of
course, you get the idea).}
}
\value{
a character vector of hex color codes where the \emph{names} attribute
corresponds to the nodes of the given igraph object. Will be used to fill in
the \code{V(net)$color} property of the \code{net} object. If there are nodes
that are part of the network object \code{net} but not present in the \code{diff}
vector, then a \emph{NA} value will be given for the color of these nodes.
}
\description{
This function splits the [-1,1] interval into \strong{2000} smaller
ones and matches each value of the \code{diff} vector to a specific hex color
code, using a spline interpolation of the colors as defined in the \code{col}
parameter.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_unobserved_model_predictions}
\alias{get_unobserved_model_predictions}
\title{Subset the model predictions to the (false) non-observed synergies}
\usage{
get_unobserved_model_predictions(model.predictions, observed.synergies)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{observed.synergies}{a character vector with elements the names of the
drug combinations that were found as synergistic}
}
\value{
a \code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{non-synergistic}
(\emph{negative results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)
}
\description{
Subset the model predictions to the (false) non-observed synergies
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_node_names}
\alias{get_node_names}
\title{Get the node names}
\usage{
get_node_names(models.dir)
}
\arguments{
\item{models.dir}{string. A directory with at least one \emph{.gitsbe} file/model.}
}
\value{
a character vector of the node names (protein and/or gene names)
}
\description{
This function uses the first \emph{.gitsbe} file that it finds inside the given
directory to output a vector of the network node names (which should be the
same for every model)
}
\examples{

models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
nodes = get_node_names(models.dir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_alt_drugname}
\alias{get_alt_drugname}
\title{Get alternative drug combination name}
\usage{
get_alt_drugname(drug.comb)
}
\arguments{
\item{drug.comb}{a string in the form \emph{drugname.1-drugname.2} (no
spaces between the names and the hyphen '-')}
}
\value{
the alternative, yet equivalent drug combination
}
\description{
Use this function on a string \emph{A-B} that represents a drug combination,
to get the reverse combination name - \emph{B-A} - for testing/checking data.
}
\examples{
drug.comb = "A-B"
alt.drug.comb = get_alt_drugname(drug.comb)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{is_comb_element_of}
\alias{is_comb_element_of}
\title{Is drug combination element of given vector?}
\usage{
is_comb_element_of(drug.comb, comb.vector)
}
\arguments{
\item{drug.comb}{a string in the form \emph{A-B} (no spaces between the names
and the hyphen '-')}

\item{comb.vector}{a character vector of drug combinations, each one in the
form \emph{drugname.1-drugname.2}}
}
\value{
logical, depending if the drug combination is element of the given
vector or not.
}
\description{
Use this function to determine if a drug combination is part of a vector of
other drug combinations. We take care only of pair-wise drug combinations and
an internal check is done for alternative drug names, e.g. we check if
\emph{A-B} combination is included, but also for \emph{B-A}.
}
\examples{
# TRUE
is_comb_element_of("A-B", c("E-F", "A-B"))
is_comb_element_of("B-A", c("E-F", "A-B"))

# FALSE
is_comb_element_of("A-B", c("E-F", "A-D"))
is_comb_element_of("A-B", c())

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_avg_state_diff_graphs}
\alias{plot_avg_state_diff_graphs}
\title{Plot the graphs from an average state differences matrix}
\usage{
plot_avg_state_diff_graphs(net, diff.mat, layout = NULL)
}
\arguments{
\item{net}{igraph graph object}

\item{diff.mat}{a matrix whose rows are \strong{vectors of average node activity
state differences} between two groups of models based on some kind of classification
(e.g. number of TP predictions) and whose names are set in the \code{rownames}
attribute of the matrix (usually denoting the different classification
groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
that predicted 2 TP synergies, if the classification is done by number of TP
predictions). The columns represent the network's node names. Could be the
result of using the function \code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}}.}

\item{layout}{a (nx2) numeric matrix of x-y coordinates (2 columns) for each
of the nodes (n) in the \code{net} igraph object. If NULL, we use the default
layout provided by \code{\link[igraph]{layout_nicely}}.}
}
\description{
This function presents a convenient way to use the function
\code{\link{plot_avg_state_diff_graph}} multiple times.
}
\seealso{
Other network plotting functions: 
\code{\link{plot_avg_link_operator_diff_graphs}()},
\code{\link{plot_avg_link_operator_diff_graph}()},
\code{\link{plot_avg_state_diff_graph_vis}()},
\code{\link{plot_avg_state_diff_graph}()}
}
\concept{network plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_based_on_mcc_clustering}
\alias{get_avg_activity_diff_based_on_mcc_clustering}
\title{Get the average activity difference based on MCC clustering}
\usage{
get_avg_activity_diff_based_on_mcc_clustering(
  models.mcc,
  models.stable.state,
  mcc.class.ids,
  models.cluster.ids,
  class.id.low,
  class.id.high,
  penalty = 0
)
}
\arguments{
\item{models.mcc}{a numeric vector of Matthews Correlation Coefficient (MCC)
scores, one for each model. The \emph{names} attribute holds the models' names.
Can be the result of using the function \code{\link{calculate_models_mcc}}.}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{mcc.class.ids}{a numeric vector of group/class ids starting from 1,
e.g. \code{c(1,2,3)} (3 MCC classes).}

\item{models.cluster.ids}{a numeric vector of cluster ids assigned to each
model. It is the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
with input the vector of the models' MCC values.}

\item{class.id.low}{integer. This number specifies the MCC class id of the
'bad' models.}

\item{class.id.high}{integer. This number specifies the MCC class id of the
'good' models and needs to be strictly higher than \code{class.id.low}.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a numeric vector with values in the [-1,1] interval (minimum and maximum
possible average difference) and with the \emph{names} attribute representing the name
of the nodes.
}
\description{
This function splits the models to 'good' and 'bad' based on an MCC value
clustering method: \emph{class.id.high} denotes the group id with the higher MCC
values (good model group) vs \emph{class.id.low} which denotes the group id with
the lower MCC values (bad model group). Then, for each network node, the function
finds the node's average activity in each of the two classes (a value in
the [0,1] interval) and then subtracts the bad class average activity value from
the good one, taking into account the given \code{penalty} factor and the
number of models in each respective model group.
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node is more \strong{inhibited} in the 'good' models compared to the
'bad' ones while a value closer to 1 means that the node is more \strong{activated}
in the 'good' models. A value closer to 0 indicates that the activity of that
node is \strong{not so much different} between the 'good' and 'bad' models and
so it won't not be a node of interest when searching for indicators of better
performance (higher MCC score/class) in the good models.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{get_synergy_biomarkers_per_cell_line}
\alias{get_synergy_biomarkers_per_cell_line}
\title{Get synergy biomarkers per cell line}
\usage{
get_synergy_biomarkers_per_cell_line(biomarkers.dirs)
}
\arguments{
\item{biomarkers.dirs}{a character vector of the biomarker directories, in the
form of \emph{\{path\}/cell_line_name/\{dir\}}. The cell line name directory
should be different for each element of the vector as we use it to fill in the
\code{rownames} of each cell line-specific \code{data.frame} object.
Inside each \emph{\{dir\}} (the directory name does not matter, but 'biomarkers'
is a good choice), we read the synergy biomarkers from a file (if it
exists and is non-empty) with the name \emph{biomarkers_per_synergy}. This file
has as first row the node names (columns) while every next row starts with the row name
(drug combination name) followed by a series of numbers from the ternary set
\{1,-1,0\}, denoting thus which nodes where found as active biomarkers for that
synergy, inhibited or not at all as biomarkers.}
}
\value{
a list of cell line-specific data frames (each element
from the list takes its name from the respective cell line).
Each cell-line specific \code{data.frame} object has as rows the
\strong{true positive predicted synergies} for that particular cell line
and columns the network nodes (should be the same for all cell lines).
Possible values for each \emph{synergy-node}
element in each cell line-specific \code{data.frame} are either \emph{1}
(\emph{active state} biomarker), \emph{-1}
(\emph{inhibited state} biomarker) or \emph{0} (not a biomarker).
}
\description{
Use this function to get the synergy biomarkers for each cell line.
The biomarkers must be stored in a single file inside each given cell line-specific
directory.
}
\examples{
dir = system.file("extdata", "AGS", "bio", package = "emba", mustWork = TRUE)
res_list = get_synergy_biomarkers_per_cell_line(biomarkers.dirs = c(dir))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_observed_synergies}
\alias{get_observed_synergies}
\title{Load the observed synergies data}
\usage{
get_observed_synergies(file, drug.combinations.tested = NULL)
}
\arguments{
\item{file}{string. The name of the file, can be a full path. See example
below for the format of an observed synergies file.}

\item{drug.combinations.tested}{a character vector with drug combinations
as elements. Default value: NULL.}
}
\value{
a character vector with elements the names of the drug combinations
that were found as synergistic
}
\description{
Use this function to read a file that has the observed synergies data and
output it to a character vector. If \code{drug.combinations.tested}
is NULL (the default), no data validation is done, otherwise we check that
the observed synergies are indeed a subset of the tested drug combinations.
}
\examples{
observed.synergies.file = system.file("extdata", "observed_synergies",
  package = "emba", mustWork = TRUE)
observed.synergies = get_observed_synergies(observed.synergies.file)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{add_numbers_above_the_bars}
\alias{add_numbers_above_the_bars}
\title{Add numbers horizontally above the bars of a barplot}
\usage{
add_numbers_above_the_bars(stats, bp, color)
}
\arguments{
\item{stats}{a numeric vector}

\item{bp}{the result of \strong{\code{barplot}} command, usually a numeric
vector or matrix}

\item{color}{string. The color for the numbers}
}
\description{
Add numbers horizontally above the bars of a barplot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}
\alias{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}
\title{Get average link operator difference matrix based on specific synergy prediction}
\usage{
get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction(
  model.predictions,
  models.link.operator,
  predicted.synergies,
  penalty = 0
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names (same order as in the
\code{model.predictions} parameter) whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
(\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).}

\item{predicted.synergies}{a character vector of the synergies (drug
combination names) that were predicted by \strong{at least one} of the models
in the dataset. It must be a subset of the column names (the drug combinations)
of the \code{model.predictions} object.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of average node link operator
differences} between two groups of models where the classification for each
individual row was based on the prediction or not of a specific synergistic
drug combination.
The row names are the predicted synergies, one per row, while the columns
represent the network's node names. Values are in the [-1,1] interval.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}}
function with the parameter \code{models.link.operator} as input in the place of
\code{models.stable.state}, since the two matrices representing the two inputs
have the same data format (rows represent models, columns represent nodes,
and each value is a number in the [0,1] interval).
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node's boolean equation has the \strong{AND NOT} link operator in the
models that predicted the specified synergy while a value closer to 1 means that
the node's boolean equation has mostly the \strong{OR NOT} link operator
in these models. A value closer to 0 indicates that the link operator in
the node's boolean equation is \strong{not so much different} between the
models that predicted the synergy and those that did not and so it won't not
be a node of interest when searching for \emph{synergy biomarkers} - nodes
whose parameterization (value of the link operator) affects the manifestation
of synergy. A value exactly equal to 0 can also mean that this node didn't
not have a link operator in its boolean equation (making it thus a non-important
node with regard to the parameterization).
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_models_synergies_tn}
\alias{calculate_models_synergies_tn}
\title{Count the non-synergies of the non-synergistic drug combinations per model (TN)}
\usage{
calculate_models_synergies_tn(unobserved.model.predictions)
}
\arguments{
\item{unobserved.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{non-synergistic}
(\emph{negative results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}
}
\value{
an integer vector with elements the number of true negative predictions
per model. The model names are given in the \emph{names} attribute (same order
as in the \emph{rownames} attribute of the unobserved.model.predictions
\code{data.frame}).
}
\description{
Since the given \code{unobserved.model.predictions} data.frame has only the
negative results, this function returns the total number of 0's \emph{and}
NA's in each row.
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_mcc}()},
\code{\link{calculate_models_mcc}()},
\code{\link{calculate_models_synergies_fn}()},
\code{\link{calculate_models_synergies_fp}()},
\code{\link{calculate_models_synergies_tp}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{get_biomarkers}
\alias{get_biomarkers}
\title{Get total biomarkers from average data differences matrix}
\usage{
get_biomarkers(diff.mat, threshold)
}
\arguments{
\item{diff.mat}{a matrix whose rows are vectors of average node data
differences between two groups of models based on some kind of classification
(e.g. number of TP predictions) and whose names are set in the \code{rownames}
attribute of the matrix (usually denoting the different classification
groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
that predicted 2 TP synergies, if the classification is done by number of TP
predictions). The columns represent the network's node names.}

\item{threshold}{numeric. A number in the [0,1] interval, above which (or
below its negative value) a biomarker will be registered in the returned result.
Values closer to 1 translate to a more strict threshold and thus less
biomarkers are found.}
}
\value{
a list with two elements:
 \itemize{
   \item \code{biomarkers.pos}: a character vector that includes the node
   names of the \emph{positive} biomarkers
   \item \code{biomarkers.neg}: a character vector that includes the node
   names of the \emph{negative} biomarkers
}
}
\description{
Use this function to find all biomarkers across multiple
performance classification group matchings based on a given threshold between
0 and 1.
}
\section{Details}{

This function uses the \code{\link{get_biomarkers_per_type}} function
to get the biomarkers (nodes) of both types (positive and negative) from the
average data differences matrix. The logic behind the biomarker selection is
that if there is at least one value in a column of the \code{diff.mat} matrix
that surpasses the threshold given, then the corresponding node (name of the
column) is returned as a biomarker.
This means that for a single node, if at least one value that represents an average data
difference (for example, the average activity state difference) between any
of the given classification group comparisons is above the given threshold
(or below the negative symmetric threshold), then a \emph{positive}
(\emph{negative}) biomarker is reported.

In the case of a node which is found to surpass the
significance threshold level given \emph{both negatively and positively},
we will keep it as a biomarker
in the category which corresponds to the \strong{comparison of the highest
classification groups}. For example, if the data comes from a model performance
classification based on the MCC score and in the comparison of the MCC classes
(1,3) the node of interest had an average difference of \emph{-0.89} (a negative
biomarker) while for the comparison of the (3,4) MCC classes it had a value
of \emph{0.91} (a positive biomarker), then we will keep that node \emph{only as a
positive biomarker}. The logic behind this is that
the 'higher' performance-wise are the classification groups that we compare,
the more sure we are that the average data difference corresponds to a
\emph{better indicator} for the type of the biomarker found.
}

\seealso{
Other biomarker functions: 
\code{\link{get_biomarkers_per_type}()}
}
\concept{biomarker functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{validate_observed_synergies_data}
\alias{validate_observed_synergies_data}
\title{Validate observed synergies data}
\usage{
validate_observed_synergies_data(observed.synergies, drug.combinations.tested)
}
\arguments{
\item{observed.synergies}{a non-empty character vector of drug combinations}

\item{drug.combinations.tested}{a non-empty character vector of drug combinations}
}
\value{
NULL if no errors found, otherwise stops execution.
}
\description{
This function checks that the observed synergies are part (a subset) of the
tested drug combinations
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_stable_state_from_models_dir}
\alias{get_stable_state_from_models_dir}
\title{Load the models stable state data}
\usage{
get_stable_state_from_models_dir(models.dir, all.ss = FALSE)
}
\arguments{
\item{models.dir}{string. A directory with \emph{.gitsbe} files/models.
\strong{Do not} include the ending path character in the string (\emph{/}).
Only files that include the string \emph{gitsbe} are parsed.}

\item{all.ss}{logical. Should all stable states be included in the returned
object? Default value is \emph{FALSE} (only the 1 stable state models are included).}
}
\value{
The format of the returned object depends on the \code{all.ss} value.
If:
\itemize{
  \item \code{all.ss} is \emph{FALSE} (default): a \code{data.frame} (nxm)
  with n models and m nodes. The row names
  specify the models names (taken from the file names without the \emph{gitsbe}
  extension) whereas the column names specify the name of the
  network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
  element} are either \emph{0} (inactive node) or \emph{1} (active node). If a
  \emph{.gitsbe} file/model has zero (0) or more than 1 stable states, a diagnostic
  message is printed and the corresponding model is discarded, i.e. it will not
  be included in the returned \code{data.frame} object.
  \item \code{all.ss} is \emph{TRUE}: a \code{tibble} object where each row
  stores a separate stable state and the columns correspond to network nodes
  (as before) with an extra last column that has the name of the model that
  produced that stable state. As such, models that have multiple stable states
  will occupy several rows with the last column having the same name/model.
  Models with no stable states are discarded.
}
}
\description{
Use this function to merge the stable states from all boolean models into a single
\code{data.frame} object. The models stable states are loaded from \emph{.gitsbe}
files that can be found inside the given \code{models.dir} directory.
}
\examples{

models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
models.stable.state = get_stable_state_from_models_dir(models.dir)
models.stable.state = get_stable_state_from_models_dir(models.dir, all.ss = TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_synergy_comparison_sets}
\alias{get_synergy_comparison_sets}
\title{Get synergy comparison sets}
\usage{
get_synergy_comparison_sets(synergy.subset.stats)
}
\arguments{
\item{synergy.subset.stats}{integer vector with values the amount of models
that predicted each synergy subset, defined as a comma-separated string of
drug combinations in the \emph{names} attribute of the vector. It can be the
result of using the function \code{\link{get_synergy_subset_stats}}.}
}
\value{
\code{data.frame} object with 3 columns. For each row, the 1st column defines a
\emph{single synergy} of interest (e.g. drug combination "A-B"), the 2nd a
\emph{synergy set} that includes the single one (e.g. the set "F-G,A-B,C-D")
and the 3rd the \emph{synergy subset} of the \emph{set} that does not include
the single synergy of the first column (e.g. "F-G,C-D").
}
\description{
This helper function identifies pairs of (\emph{set}, \emph{subset}) for each
synergy (implicitly given through the \code{synergy.subset.stats} object) where
each respective \emph{subset} misses just one synergy from the larger \emph{set}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_models_synergies_fn}
\alias{calculate_models_synergies_fn}
\title{Count the non-synergies of the observed synergies per model (FN)}
\usage{
calculate_models_synergies_fn(observed.model.predictions)
}
\arguments{
\item{observed.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found/observed as \strong{synergistic}
(\emph{negative results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}
}
\value{
an integer vector with elements the number of false negative predictions
per model. The model names are given in the \emph{names} attribute (same order
as in the \emph{rownames} attribute of the observed.model.predictions
\code{data.frame}).
}
\description{
Since the given \code{observed.model.predictions} data.frame has only the
positive results, this function returns the total number of 0's \emph{and}
NA's in each row.
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_mcc}()},
\code{\link{calculate_models_mcc}()},
\code{\link{calculate_models_synergies_fp}()},
\code{\link{calculate_models_synergies_tn}()},
\code{\link{calculate_models_synergies_tp}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_synergy_scores}
\alias{get_synergy_scores}
\title{Get synergy scores from file}
\usage{
get_synergy_scores(file_name, file_type = "ensemblewise")
}
\arguments{
\item{file_name}{string. The name of the file, can be a full path.}

\item{file_type}{string. The type of input file, can be either \emph{ensemblewise}
(default option) or \emph{modelwise}.}
}
\value{
a \code{tibble} containing a representation of the data in the file.
}
\description{
Use this function to read \href{https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output}{Drabme's
output files} that have synergy scores for a list of tested drug perturbations.
}
\details{
Two types of files can be read: the \emph{model-wise} and the \emph{ensemble-wise}
synergies:
\itemize{
  \item The \emph{model-wise} synergies data is structured as a 3-column
table, first being the names of the tested drug combinations, second the
number of models that predicted that combination as synergistic and the third
the number of models that predicted the drug combination as non-synergistic.
  \item The \emph{ensemble-wise} synergies data is structured as a 2-column table,
the first being the names of drug combinations and the second the ensemble-wise
synergy scores for those perturbations.
 }

Note that no matter the file type, the first line of the file is always skipped
and the columns must be \emph{tab-separated}.
See example below for the format of such files generated by the
\href{https://druglogics.github.io/druglogics-doc/drabme.html}{Drabme} module.
}
\examples{

ensemblewise_synergies_file = system.file("extdata", "ensemblewise_synergies",
  package = "emba", mustWork = TRUE)
modelwise_synergies_file = system.file("extdata", "modelwise_synergies",
  package = "emba", mustWork = TRUE)
data_ensemble = get_synergy_scores(ensemblewise_synergies_file) # file_type = "ensemblewise"
data_modelwise = get_synergy_scores(modelwise_synergies_file, file_type = "modelwise")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_based_on_synergy_set_cmp}
\alias{get_avg_activity_diff_based_on_synergy_set_cmp}
\title{Get the average activity difference based on the comparison of two synergy sets}
\usage{
get_avg_activity_diff_based_on_synergy_set_cmp(
  synergy.set.str,
  synergy.subset.str,
  model.predictions,
  models.stable.state,
  penalty = 0
)
}
\arguments{
\item{synergy.set.str}{a string of drug combinations, comma-separated. The
number of the specified combinations must be larger than the ones defined
in the \code{synergy.subset.str} parameter. They also must be included in the
tested drug combinations, i.e. the columns of the \code{model.predictions}
parameter.}

\item{synergy.subset.str}{a string of drug combinations, comma-separated.
There must be at least one combination defined and all of them should also
be included in the \code{synergy.set.str} parameter.}

\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a numeric vector with values in the [-1,1] interval (minimum and
maximum possible average difference) and with the names attribute
representing the name of the nodes.
}
\description{
This function splits the models to 'good' and 'bad' based on the predictions
of two different synergy sets, one of them being a subset of the other.
The 'good' models are those that predict the \code{synergy.set.str}
(e.g. "A-B,A-C,B-C") while the 'bad' models are those that predict the
\code{synergy.subset.str} (e.g. "A-B,B-C"). Then, for each network node,
the function finds the node's average activity in each of the two classes
(a value in the [0,1] interval) and then subtracts the bad class average
activity value from the good one, taking into account the given \code{penalty}
factor and the number of models in the 'good' and 'bad' class respectively.
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node is more \strong{inhibited} in the models that predicted the extra
synergy(-ies) that are included in the \code{synergy.set.str} but not in the
\code{synergy.subset.str}, whereas a value closer to 1 means that the node is
more \strong{activated} in these models. These nodes are \strong{potential
biomarkers} because their activity state can influence the prediction
performance of a model and make it predict the extra synergy(-ies).
A value closer to 0 indicates that the activity of that
node is \strong{not so much different} between the models that predicted the
synergy set and those that predicted it's subset, so it won't be a node
of interest when searching for potential biomarkers for the extra synergy(-ies).
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_avg_state_diff_graph_vis}
\alias{plot_avg_state_diff_graph_vis}
\title{Plot the graph of average state differences (visNetwork)}
\usage{
plot_avg_state_diff_graph_vis(net, diff, nodes.size = 20, title)
}
\arguments{
\item{net}{igraph graph object (to be translated to a \code{visNetwork} object)}

\item{diff}{numeric vector. Every value must be in the [-1,1] interval and
represents the average activity difference of each node. The node names have
to be specified in the \emph{names} attribute of the given vector. For example,
\code{diff} could be the result of using the function
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}}}

\item{nodes.size}{an integer specifying the size of the nodes. Default value: 20.}

\item{title}{string. The title of the \code{visNetwork} plot.}
}
\description{
This function uses the \code{\link[visNetwork]{visNetwork}} package to plot a
network of nodes. The nodes are positioned by default in a hierarchical layout
and their colors are derived using the \code{diff} values and the
\code{\link{get_node_colors}} function. The color of each node indicates how
much more inhibited or active that node is, when comparing the average model
classified in the 'good' category vs the average 'bad' one.
}
\examples{
topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
net = construct_network(topology.file)
diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51,-0.86,-0.80)
names(diff) = c("A","C","B","D","W","I","E","J","F","K")
plot_avg_state_diff_graph_vis(net, diff, title = "TEST")

}
\seealso{
Other network plotting functions: 
\code{\link{plot_avg_link_operator_diff_graphs}()},
\code{\link{plot_avg_link_operator_diff_graph}()},
\code{\link{plot_avg_state_diff_graphs}()},
\code{\link{plot_avg_state_diff_graph}()}
}
\concept{network plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_link_operator_diff_mat_based_on_mcc_clustering}
\alias{get_avg_link_operator_diff_mat_based_on_mcc_clustering}
\title{Get average link operator difference matrix based on MCC clustering}
\usage{
get_avg_link_operator_diff_mat_based_on_mcc_clustering(
  models.mcc,
  models.link.operator,
  num.of.mcc.classes,
  penalty = 0
)
}
\arguments{
\item{models.mcc}{a numeric vector of Matthews Correlation Coefficient (MCC)
scores, one for each model. The \emph{names} attribute holds the models' names.
Can be the result of using the function \code{\link{calculate_models_mcc}}.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names (same order as in the
\code{models.mcc} parameter) whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
(\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).}

\item{num.of.mcc.classes}{numeric. A positive integer larger than 2 that
signifies the number of mcc classes (groups) that we should split the models
MCC values.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a matrix whose rows are \strong{vectors of average node link operator
differences} between two groups of models where
the classification was based on the models' MCC values.
Rows represent the different classification group matchings, e.g. (1,2) means the
models that belonged to the 1st group of MCC values vs the models that
belonged to the 2nd group. The columns represent the network's node names.
Values are in the [-1,1] interval.
}
\description{
This function uses the \code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}}
function with the parameter \code{models.link.operator} as input in the place of
\code{models.stable.state}, since the two matrices representing the two inputs
have the same data format (rows represent models, columns represent nodes,
and each value is a number in the [0,1] interval).
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node's boolean equation has the \strong{AND NOT} link operator in the
'good' models compared to the 'bad' ones while a value closer to 1 means that
the node's boolean equation has mostly the \strong{OR NOT} link operator
in the 'good' models. A value closer to 0 indicates that the link operator in
the node's boolean equation is \strong{not so much different} between the
'good' and 'bad' models and so it won't not be a node of interest when
searching for indicators of better performance (higher average MCC value)
in the parameterization of the good models (the boolean equations). A value
exactly equal to 0 can also mean that this node didn't not have a link operator
in its boolean equation, again making it a non-important indicator of difference
in model performance.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{filter_network}
\alias{filter_network}
\title{Filter the network's vertices}
\usage{
filter_network(net, nodes, level)
}
\arguments{
\item{net}{an igraph object.}

\item{nodes}{character vector of node names. It must be a subset of the nodes
of the \emph{net} object.}

\item{level}{integer. Can be only 0, 1 or 2 and specifies the neighbourhood
depth of the result graph.}
}
\value{
an induced subgraph of the \code{net} igraph object.
}
\description{
Produce an \strong{induced subgraph} of the given \emph{net} igraph object.
How many vertices/nodes will be kept in the result graph object is determined
by the initial nodes given and the level provided. A level equal to 0 corresponds
to a subgraph with only the given nodes, a level equal to 1 to a subgraph with
the nodes + their neighbors (the closed neighbourhood set where every node is within
1 edge distnace from the given ones) and a level equal to 2 to a subgraph
with the nodes + their neighbors + the nodes neighbor neighbors!
(so the neighbourhood of the neighbourhood or every node is within 2 edges
distance from the given ones).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{get_synergy_biomarkers_from_dir}
\alias{get_synergy_biomarkers_from_dir}
\title{Get synergy biomarkers from dir}
\usage{
get_synergy_biomarkers_from_dir(
  predicted.synergies,
  biomarkers.dir,
  models.dir = NULL,
  node.names = NULL
)
}
\arguments{
\item{predicted.synergies}{a character vector of the synergies (drug
combination names) that were predicted by \strong{at least one} of the models
in the dataset.}

\item{biomarkers.dir}{string. It specifies the full path name of the
directory which holds the biomarker files (without the ending character
\emph{/}). The biomarker files must be formatted as:
\emph{\%drug.comb\%_biomarkers_active} or
\emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is an element
of the \code{predicted.synergies} vector.}

\item{models.dir}{string. A directory with \emph{.gitsbe} files/models. It's
needed in order to call \code{\link{get_node_names}}.}

\item{node.names}{a character vector which has the names of the nodes. If it's
not NULL, then it will be used instead of the \code{models.dir} parameter.
The \code{node.names} should include all the nodes that are reported as
biomarkers in the biomarker files inside the \code{biomarkers.dir} directory.
Note that the biomarker nodes in the files will be included in the returned
\code{data.frame} object no matter the \code{node.names} specified.
Default value: NULL.}
}
\value{
a data.frame, whose columns represent the network nodes and the
rows the predicted synergies. Possible values for each \emph{synergy-node}
element are either \emph{1} (\emph{active state} biomarker), \emph{-1}
(\emph{inhibited state} biomarker) or \emph{0} (not a biomarker or the node
is not at all present in the network or the drug combination is not a
synergistic one).
}
\description{
This function reads the synergy biomarker files inside the given directory and merges
the results into a \code{data.frame} which it returns. This functions should
be used when the synergy biomarker results are in separate files inside the
directory given (see \code{biomarkers.dir} parameter).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{get_perf_biomarkers_per_cell_line}
\alias{get_perf_biomarkers_per_cell_line}
\title{Get performance biomarkers per cell line}
\usage{
get_perf_biomarkers_per_cell_line(biomarkers.dirs, node.names)
}
\arguments{
\item{biomarkers.dirs}{a character vector of the biomarker directories, in the
form of \emph{\{path\}/cell_line_name/\{dir\}}. The cell line name directory
should be different for each element of the vector as we use it to fill in the
\code{rownames} of the result \code{data.frame} object. Inside each \emph{\{dir\}}
(the directory name does not matter, but 'biomarkers' is a good choice),
we read the biomarkers from two files (if they exist and are non-empty):
\emph{biomarkers_active} and \emph{biomarkers_inhibited}, which have the
active and inhibited performance biomarkers for each cell line (these files
have a list of node names/biomarkers, one in each line).}

\item{node.names}{a character vector of the node names used in the analysis.
The biomarker names taken from the files inside the given directories must be
a subset of this vector.}
}
\value{
a data.frame, whose columns represent the network nodes and the
rows the cell lines. Possible values for each \emph{cell line-node}
element are either \emph{1} (\emph{active state} biomarker), \emph{-1}
(\emph{inhibited state} biomarker) or \emph{0} (not a biomarker).
}
\description{
Use this function to get the performance biomarkers from the respective
files inside the given list of directories.
}
\examples{
dir = system.file("extdata", "AGS", "bio", package = "emba", mustWork = TRUE)
res = get_perf_biomarkers_per_cell_line(biomarkers.dirs = c(dir),
  node.names = paste0("x", 1:20))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_avg_activity_diff_based_on_specific_synergy_prediction}
\alias{get_avg_activity_diff_based_on_specific_synergy_prediction}
\title{Get average activity difference based on specific synergy prediction}
\usage{
get_avg_activity_diff_based_on_specific_synergy_prediction(
  model.predictions,
  models.stable.state,
  drug.comb,
  penalty = 0
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names whereas the column names specify the network nodes
(gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.}

\item{drug.comb}{string. The drug combination which will be used to split
the models. It must be included in the column names of the \code{model.predictions}
object.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a numeric vector with values in the [-1,1] interval (minimum and maximum
possible average difference) and with the \emph{names} attribute representing
the name of the nodes.
}
\description{
Given a specific drug combination, this function splits the models to
good (those that predicted that particular combination, i.e. found it
as synergistic - a value of \emph{1} in the \code{model.predictions}) and
bad (those that found it as non-synergistic - a value of \emph{0} in the
\code{model.predictions}). The models whose predicted value for that synergy is marked as
\emph{NA} are excluded from the analysis. Then, for each network node, the
function finds the node's average activity in each of the two model groups (a
value in the [0,1] interval) and then subtracts the bad group's average
activity value from the good one, taking into account the given \code{penalty}
factor and the number of models in each respective model group.
}
\section{Details}{

So, if a node has a value close to -1 it means that on average,
this node is more \strong{inhibited} in the models that predicted the specific
drug combination given, whereas a value closer to 1 means that the node is more
\strong{activated} in these models. A value closer to 0 indicates that the activity of that
node is \strong{not so much different} between the models that predicted the synergy and
those that did not and so it won't not be a node of interest when searching
for \emph{synergy biomarkers} - nodes whose activity is important for the
manifestation of the synergy.
}

\seealso{
\code{\link{get_vector_diff}}

Other average data difference functions: 
\code{\link{get_avg_activity_diff_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_activity_diff_based_on_tp_predictions}()},
\code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}()},
\code{\link{get_avg_link_operator_diff_based_on_synergy_set_cmp}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_mcc_clustering}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction}()},
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}()}
}
\concept{average data difference functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{assign_link_operator_value_to_equation}
\alias{assign_link_operator_value_to_equation}
\title{Assign link operator value to boolean equation}
\usage{
assign_link_operator_value_to_equation(equation)
}
\arguments{
\item{equation}{string. The boolean equation in the form
\eqn{Target *= (Activator or Activator or...)
and not (Inhibitor or Inhibitor or...)"}}
}
\value{
\strong{1} if the \code{equation} has the '\emph{or not}' link operator,
\strong{0} if the \code{equation} has the '\emph{and not}' link operator and
\strong{NA} if it has neither.
}
\description{
Assign link operator value to boolean equation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_avg_link_operator_diff_graph}
\alias{plot_avg_link_operator_diff_graph}
\title{Plot the graph of average link operator differences (igraph)}
\usage{
plot_avg_link_operator_diff_graph(net, diff, layout = NULL, title)
}
\arguments{
\item{net}{igraph graph object}

\item{diff}{numeric vector. Every value is in the [-1,1] interval and
represents the average link operator value difference of each node. The node
names have to be specified in the \emph{names} attribute of the given vector.
For example, \code{diff} could be the result of using the function
\code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}} and
getting one vector row from the output matrix.
A value closer to -1 means that the 'good' models have more of the \emph{AND NOT}
link operator in their respective boolean equations while a value closer to 1
means that the 'good' models have more of the \emph{OR NOT} link operator.}

\item{layout}{a (nx2) numeric matrix of x-y coordinates (2 columns) for each
of the nodes (n) in the \code{net} igraph object. If NULL, we use the default
layout provided by \code{\link[igraph]{layout_nicely}}.}

\item{title}{string. The title of the igraph plot}
}
\description{
This function uses the \code{\link[igraph]{plot.igraph}} package to plot a network
of nodes. The nodes are positioned according to the specified coordinates
given by the \code{layout} parameter and the colors are derived using the
\code{diff} values and the \code{\link{get_node_colors}} function. The color
of each node indicates if the node's boolean function has on average the
\emph{AND NOT} or the \emph{OR NOT} link operator when comparing the average
model classified in the 'good' category vs the average bad' one. A non-colored
node (white) will indicate nodes that do not have the link operator in their
respective boolean equation (where they function as the target).
}
\examples{
topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
net = construct_network(topology.file)
diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51)
names(diff) = c("A","C","B","D","W","I","E","J")
plot_avg_link_operator_diff_graph(net, diff, title = "TEST")

}
\seealso{
\code{\link{get_node_colors}}

Other network plotting functions: 
\code{\link{plot_avg_link_operator_diff_graphs}()},
\code{\link{plot_avg_state_diff_graph_vis}()},
\code{\link{plot_avg_state_diff_graphs}()},
\code{\link{plot_avg_state_diff_graph}()}
}
\concept{network plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_synergy_subset_stats}
\alias{get_synergy_subset_stats}
\title{Find the number of predictive models for every synergy subset}
\usage{
get_synergy_subset_stats(model.predictions, synergies)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models).}

\item{synergies}{a character vector with elements the synergistic drug
combinations. Note that these synergies should be a subset of the column
names of the \code{model.predictions} data.frame.}
}
\value{
an integer vector with elements the number of models the predicted
each synergy subset. The \emph{names} attribute has the names of each
synergistic drug combination subset, which are the drug combinations comma
separated (e.g. 'A-B,C-D').
}
\description{
Use this function to find for each possible subset of drug combinations out
of a given list of synergies, the number of models that predicted it given
the models' predictions. So, if for example the set of synergies is this one:
\{'A-B','C-D','E-F'\}, we want to know how many models predicted none of them,
just the single subsets (e.g. the \{'A-B'\}),
the two-element subsets (e.g. the \{'A-B','C-D'\}) and all 3 of them.
}
\section{Details}{

Note that if the \code{synergies} vector has more than 10-15 elements, then
this function might take long time to execute even with an optimal
implementation of \code{\link{count_models_that_predict_synergies}}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{make_barplot_on_models_stats}
\alias{make_barplot_on_models_stats}
\title{Bar plot of model stats}
\usage{
make_barplot_on_models_stats(
  models.stats,
  cell.line = NULL,
  title,
  xlab,
  ylab,
  cont.values = FALSE,
  threshold = 0,
  ylim.add = 0
)
}
\arguments{
\item{models.stats}{table object, the result of using \link[base]{table} on
a (numeric) vector. Usually it represents some models statistics summary -
counts for each TP prediction value for example.}

\item{cell.line}{string. The name of the cell line to be used in the title
of the produced plot. Default value: NULL (the cell line name will not be
added to the title)}

\item{title}{string. The title of the plot}

\item{xlab}{string. The title of the x-axis}

\item{ylab}{string. The title of the y-axis}

\item{cont.values}{logical. If TRUE, the values of the x-axis will be trimmed
to 3 digits after the decimal point. Default value: FALSE.}

\item{threshold}{integer. Values from the \code{model.stats} that are \emph{less
or equal} to the threshold will be pruned. Use it when there too many
categories and the figure appears too dense. Default value: 0}

\item{ylim.add}{integer. Signifies the height to add to the upper
\code{ylim} parameter on the barplot, in addition to the maximum bar height
across the whole plot. Default value is 0.}
}
\description{
Use this function to produce a bar plot when the input is the result of using
the \link[base]{table} function to a numeric vector
}
\examples{
x = c(rep(1,100), rep(2,423), rep(3,231), rep(NaN,531))
make_barplot_on_models_stats(models.stats = table(x, useNA = "ifany"),
title = "True Positives Distribution across models",
xlab = "Number of TP values", ylab = "Number of models")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff.R
\name{get_vector_diff}
\alias{get_vector_diff}
\title{Calculate difference vector with penalty term}
\usage{
get_vector_diff(vec1, vec2, m1 = 1, m2 = 1, penalty = 0)
}
\arguments{
\item{vec1}{numeric vector}

\item{vec2}{numeric vector}

\item{m1}{integer > 0}

\item{m2}{integer > 0}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty (\code{m1,m2} don't matter) and a value of 1 is the strickest possible
penalty. Default value is 0.}
}
\value{
the vector of differences between the two given vectors based on the
formula: \deqn{(vec1 - vec2) * w}, where \eqn{w = (min(m1,m2)/max(m1,m2))^p}
and \eqn{p = penalty}.

See also related \href{https://math.stackexchange.com/questions/3547139/formula-for-weighted-average-difference}{StackOverflow question}.
If \code{vec1} has \code{names}, the returned vector will have the same names
attribute as \code{vec1}.
}
\description{
This function calculates the difference between two given numeric vectors while
adding a penalty term (weight) to account for the number of models/instances that
each vector's values were calculated from. Thus, if the models/instances are
disproportionate and a penalty is included, the difference vector's values will
be changed accordingly to reflect that.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{print_model_and_drug_stats}
\alias{print_model_and_drug_stats}
\title{Print model and drug statistics}
\usage{
print_model_and_drug_stats(drug.combs, models, nodes, html.output)
}
\arguments{
\item{drug.combs}{integer. Number of drug combinations tested}

\item{models}{integer. Number of models tested}

\item{nodes}{integer. Number of network nodes}

\item{html.output}{logical. If TRUE, the printed output will look nice in an
HTML document}
}
\description{
Use this function to pretty print in an R notebook useful statistics for the
ensemble model analysis: how many drug combinations were tested by each model,
the number of models used and how many nodes each boolean network model had.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{construct_network}
\alias{construct_network}
\title{Construct igraph network graph}
\usage{
construct_network(topology.file, models.dir = NULL)
}
\arguments{
\item{topology.file}{string. The name of the .sif file (can be a full path
name).}

\item{models.dir}{string. A dir with \emph{.gitsbe} files/models. Default
value: NULL. If specified, it is used for the validation of the node names.}
}
\value{
an igraph graph object representing the network as defined in the
topology file
}
\description{
Use this function to create an igraph graph object based on the topology .sif
file given. It automatically sets various visualization graph properties and
checks if the node names from the topology file are the same as in the models
inside the given \code{models.dir} (if not NULL).
}
\examples{
topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
net = construct_network(topology.file)

}
\seealso{
\code{\link[igraph]{graph_from_data_frame}},
\code{\link{get_edges_from_topology_file}},
\code{\link{get_node_names}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{calculate_models_mcc}
\alias{calculate_models_mcc}
\title{Calculate the Matthews correlation coefficient for each model}
\usage{
calculate_models_mcc(
  observed.model.predictions,
  unobserved.model.predictions,
  number.of.drug.comb.tested
)
}
\arguments{
\item{observed.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found as \strong{synergistic}
(\emph{positive results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{unobserved.model.predictions}{\code{data.frame} object with rows the models
and columns the drug combinations that were found as \strong{non-synergistic}
(\emph{negative results}). Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models)}

\item{number.of.drug.comb.tested}{numeric. The total number of drug
combinations tested, which should be equal to the sum of the columns of the
\code{observed.model.predictions} and the \code{unobserved.model.predictions}.}
}
\value{
a numeric vector of MCC values, each value being in the [-1,1]
interval. The \emph{names} attribute holds the models' names
if applicable (i.e. the input \code{data.frames} have \emph{rownames}).
}
\description{
Calculate the Matthews correlation coefficient for each model
}
\seealso{
Other confusion matrix calculation functions: 
\code{\link{calculate_mcc}()},
\code{\link{calculate_models_synergies_fn}()},
\code{\link{calculate_models_synergies_fp}()},
\code{\link{calculate_models_synergies_tn}()},
\code{\link{calculate_models_synergies_tp}()}
}
\concept{confusion matrix calculation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_avg_state_diff_graph}
\alias{plot_avg_state_diff_graph}
\title{Plot the graph of average state differences (igraph)}
\usage{
plot_avg_state_diff_graph(net, diff, layout = NULL, title)
}
\arguments{
\item{net}{igraph graph object}

\item{diff}{numeric vector. Every value is in the [-1,1] interval and
represents the average activity difference of each node. The node names have
to be specified in the \emph{names} attribute of the given vector. For example,
\code{diff} could be the result of using the function
\code{\link{get_avg_activity_diff_based_on_tp_predictions}}.}

\item{layout}{a (nx2) numeric matrix of x-y coordinates (2 columns) for each
of the nodes (n) in the \code{net} igraph object. If NULL, we use the default
layout provided by \code{\link[igraph]{layout_nicely}}.}

\item{title}{string. The title of the igraph plot}
}
\description{
This function uses the \code{\link[igraph]{plot.igraph}} package to plot a network
of nodes. The nodes are positioned according to the specified coordinates
given by the \code{layout} parameter and the colors are derived using the
\code{diff} values and the \code{\link{get_node_colors}} function. The color
of each node indicates how much more inhibited or active that node is, when
comparing the average model classified in the 'good' category vs the average
'bad' one.
}
\examples{
topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
net = construct_network(topology.file)
diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51,-0.86,-0.80)
names(diff) = c("A","C","B","D","W","I","E","J","F","K")
plot_avg_state_diff_graph(net, diff, title = "TEST")

}
\seealso{
\code{\link{get_node_colors}}

Other network plotting functions: 
\code{\link{plot_avg_link_operator_diff_graphs}()},
\code{\link{plot_avg_link_operator_diff_graph}()},
\code{\link{plot_avg_state_diff_graph_vis}()},
\code{\link{plot_avg_state_diff_graphs}()}
}
\concept{network plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomarkers.R
\name{get_biomarkers_per_type}
\alias{get_biomarkers_per_type}
\title{Get biomarkers from average data differences matrix (per type)}
\usage{
get_biomarkers_per_type(diff.mat, threshold, type)
}
\arguments{
\item{diff.mat}{a matrix whose rows are vectors of average node data
differences between two groups of models based on some kind of classification
(e.g. number of TP predictions) and whose names are set in the \code{rownames}
attribute of the matrix (usually denoting the different classification
groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
that predicted 2 TP synergies, if the classification is done by number of TP
predictions). The columns represent the network's node names.}

\item{threshold}{numeric. A number in the [0,1] interval, above which (or
below its negative value) a biomarker will be registered in the returned result.
Values closer to 1 translate to a more strict threshold and thus less
biomarkers are found.}

\item{type}{character. Accepted values are \emph{positive} or \emph{negative}.}
}
\value{
a character vector that includes the node names that were found
either as \emph{positive} or \emph{negative}.
}
\description{
Use this function to find either positive or negative biomarkers across multiple
performance classification group matchings based on a given threshold between
0 and 1.
}
\details{
The logic behind the biomarker selection is that if there is at least one value
in a column of the \code{diff.mat} matrix that surpasses the threshold given, then the
corresponding node (name of the column) is return as a biomarker. This means
that for a single node, if at least one value that represents an average data
difference (for example, the average activity state difference) between any
of the given classification group comparisons is above the given threshold (or
below the negative symmetric threshold), then a \emph{positive} (\emph{negative})
biomarker is reported.
}
\seealso{
Other biomarker functions: 
\code{\link{get_biomarkers}()}
}
\concept{biomarker functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{get_link_operators_from_models_dir}
\alias{get_link_operators_from_models_dir}
\title{Load the models boolean equation link operator data}
\usage{
get_link_operators_from_models_dir(
  models.dir,
  remove.equations.without.link.operator = TRUE
)
}
\arguments{
\item{models.dir}{string. A directory path with \emph{.gitsbe} files/models.
\strong{Do not} include the ending path character in the string (\emph{/}).}

\item{remove.equations.without.link.operator}{logical. Should we keep the
nodes (columns in the returned matrix) which do not have both type of
regulators (so no link operator)? Default value: TRUE (remove these nodes).}
}
\value{
a \code{data.frame} (nxm) with n models and m nodes. The row names
specify the models' names whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
element} are either \emph{0} (\strong{and not} link operator), \emph{1}
(\strong{or not} link operator) or \emph{0.5} if the node is not targeted by
both activating and inhibiting regulators (no link operator).
}
\description{
Use this function to merge the link operator data used in the boolean equations
of the models into a single \code{data.frame} object. Every boolean model is defined by a series
of boolean equations in the form \eqn{Target *= (Activator or Activator or...)
and not (Inhibitor or Inhibitor or...)"}. The \strong{link operator} can be
either \emph{and not}, \emph{or not} or non-existent if the target has only
activating regulators or only inhibiting ones (the \emph{not} remains in the
latter case). The models are loaded from \emph{.gitsbe} files (and only these)
that can be found inside the given \code{models.dir} directory.
}
\examples{

models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
models.link.operator = get_link_operators_from_models_dir(models.dir)
models.link.operator.with.extra.nodes =
  get_link_operators_from_models_dir(models.dir, FALSE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{biomarker_tp_analysis}
\alias{biomarker_tp_analysis}
\title{Biomarker analysis based on TP model classification}
\usage{
biomarker_tp_analysis(
  model.predictions,
  models.stable.state,
  models.link.operator = NULL,
  observed.synergies,
  threshold,
  penalty = 0.1
)
}
\arguments{
\item{model.predictions}{a \code{data.frame} object with rows the models and
columns the drug combinations. Possible values for each \emph{model-drug combination
element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
predicted) or \emph{NA} (couldn't find stable states in either the drug
combination inhibited model or in any of the two single-drug inhibited models).}

\item{models.stable.state}{a \code{data.frame} (nxm) with n models and m nodes.
The row names specify the models' names whereas the column names specify the
network nodes (gene, proteins, etc.). Possible values for each \emph{model-node element}
can be between \emph{0} (inactive node) and \emph{1} (active node) inclusive.
Note that the rows (models) have to be in the same order as in the \code{model.predictions}
parameter.}

\item{models.link.operator}{a \code{data.frame} (nxm) with n models and m nodes. The row
names specify the models' names (same order as in the \code{model.predictions}
parameter) whereas the column names specify
the network nodes (gene, proteins, etc.). Possible values for each
\emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
\emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
by both activating and inhibiting regulators (no link operator). Default value:
NULL (no analysis on the models parameterization regarding the mutation of the
boolean equation link operator will be done).}

\item{observed.synergies}{a character vector with elements the names of the
drug combinations that were found as synergistic. This should be a subset of
the tested drug combinations, that is the column names of the \code{model.predictions}
parameter.}

\item{threshold}{numeric. A number in the [0,1] interval, above which (or
below its negative value) a biomarker will be registered in the returned result.
Values closer to 1 translate to a more strict threshold and thus less
biomarkers are found.}

\item{penalty}{value between 0 and 1 (inclusive). A value of 0 means no
penalty and a value of 1 is the strickest possible penalty. Default value is 0.1.
This penalty is used as part of a weighted term to the difference in a value of
interest (e.g. activity or link operator difference) between two group of
models, to account for the difference in the number of models from each
respective model group.}
}
\value{
a list with various elements:
\itemize{
  \item \code{predicted.synergies}: a character vector of the synergies (drug
  combination names) that were predicted by \strong{at least one} of the models
  in the dataset.
  \item \code{models.synergies.tp}: an integer vector of true positive (TP)
  values, one for each model.
  \item \code{diff.tp.mat}: a matrix whose rows are \strong{vectors of
  average node activity state differences} between two groups of models where
  the classification was based on the number of true positive predictions.
  Rows represent the different classification group matchings, e.g. (1,2) means the
  models that predicted 1 TP synergy vs the models that predicted 2 TP
  synergies and the columns represent the network's node names.
  Values are in the [-1,1] interval.
  \item \code{biomarkers.tp.active}: a character vector whose elements are
  the names of the \emph{active state} biomarkers. These nodes appear as more
  active in the better performance models.
  \item \code{biomarkers.tp.inhibited}: a character vector whose elements are
  the names of the \emph{inhibited state} biomarkers. These nodes appear as more
  inhibited in the better performance models.
  \item \code{diff.link.tp.mat}: a matrix whose rows are \strong{vectors of
  average node link operator differences} between two groups of models where
  the classification was based on the number of true positive predictions.
  Rows represent the different classification group matchings, e.g. (1,2) means the
  models that predicted 1 TP synergy vs the models that predicted 2 TP
  synergies and the columns represent the network's node names.
  Values are in the [-1,1] interval.
  \item \code{biomarkers.tp.or}: a character vector whose elements are
  the names of the \emph{OR} link operator biomarkers. These nodes have
  mostly the \emph{OR} link operator in their respective boolean equations
  in the better performance models.
  \item \code{biomarkers.tp.and}: a character vector whose elements are
  the names of the \emph{AND} link operator biomarkers. These nodes have
  mostly the \emph{AND} link operator in their respective boolean equations
  in the better performance models.
}
}
\description{
Use this function to perform a full biomarker analysis on an ensemble boolean model
dataset where the model classification is based on the number of \emph{true
positive} (TP) predictions. This analysis enables the discovery of \emph{performance
biomarkers}, nodes whose activity and/or boolean model parameterization (link
operator) affects the prediction performance of the models (as measured by
the number of TPs).
}
\seealso{
Other general analysis functions: 
\code{\link{biomarker_mcc_analysis}()},
\code{\link{biomarker_synergy_analysis}()}
}
\concept{general analysis functions}
