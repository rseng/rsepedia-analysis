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
