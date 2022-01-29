
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mikropml <a href='http://www.schlosslab.org/mikropml/'><img src='man/figures/logo.png' align="right" height="120" /></a>

> meek-ROPE em el

User-Friendly R Package for Supervised Machine Learning Pipelines

<!-- badges: start -->

[![check](https://github.com/SchlossLab/mikropml/workflows/check/badge.svg)](https://github.com/SchlossLab/mikropml/actions?query=workflow%3Acheck+branch%3Amain)
[![codecov](https://codecov.io/gh/SchlossLab/mikropml/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SchlossLab/mikropml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/SchlossLab/mikropml/blob/main/LICENSE.md)
[![CRAN](https://img.shields.io/cran/v/mikropml?color=blue&label=CRAN&logo=R)](https://CRAN.R-project.org/package=mikropml)
[![Conda](https://img.shields.io/conda/vn/conda-forge/r-mikropml)](https://anaconda.org/conda-forge/r-mikropml)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mikropml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03073/status.svg)](https://doi.org/10.21105/joss.03073)
<!-- badges: end -->

An interface to build machine learning models for classification and
regression problems. `mikropml` implements the ML pipeline described by
[Topçuoğlu *et al.* (2020)](https://doi.org/doi:10.1128/mBio.00434-20)
with reasonable default options for data preprocessing, hyperparameter
tuning, cross-validation, testing, model evaluation, and interpretation
steps. See the [website](http://www.schlosslab.org/mikropml/) for more
information, documentation, and examples.

## Installation

You can install the latest release from
[CRAN](https://cran.r-project.org/package=mikropml):

``` r
install.packages('mikropml')
```

or the development version from
[GitHub](https://github.com/SchlossLab/mikRopML):

``` r
# install.packages("devtools")
devtools::install_github("SchlossLab/mikropml")
```

or install from a terminal using
[conda](https://docs.conda.io/projects/conda/en/latest/index.html):

``` bash
conda install -c conda-forge r-mikropml
```

### Dependencies

  - Imports: caret, dplyr, e1071, glmnet, kernlab, MLmetrics,
    randomForest, rlang, rpart, stats, utils, xgboost
  - Suggests: doFuture, foreach, future, future.apply, ggplot2, knitr,
    progress, progressr, purrr, rmarkdown, testthat, tidyr

## Usage

Check out the [introductory
vignette](http://www.schlosslab.org/mikropml/articles/introduction.html)
for a quick start tutorial. For a more in-depth discussion, read [all
the vignettes](http://www.schlosslab.org/mikropml/articles/index.html)
and/or take a look at the [reference
documentation](http://www.schlosslab.org/mikropml/reference/index.html).

You can watch the Riffomonas Project series of [video
tutorials](https://www.youtube.com/playlist?list=PLmNrK_nkqBpKpzb9-vI4V7SdXC-jXEcmg)
covering mikropml and other skills related to machine learning.

We also provide an [example Snakemake
workflow](https://github.com/SchlossLab/mikropml-snakemake-workflow) for
running `mikropml` on an HPC.

## Help & Contributing

If you come across a bug, [open an
issue](https://github.com/SchlossLab/mikropml/issues) and include a
[minimal reproducible example](https://www.tidyverse.org/help/).

If you’d like to contribute, see our guidelines
[here](http://www.schlosslab.org/mikropml/CONTRIBUTING.html).

## Code of Conduct

Please note that the mikropml project is released with a [Contributor
Code of
Conduct](http://www.schlosslab.org/mikropml/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

## License

The mikropml package is licensed under [the MIT
license](https://github.com/SchlossLab/mikropml/blob/main/LICENSE.md).
Text and images included in this repository, including the mikropml
logo, are licensed under the [CC BY 4.0
license](https://creativecommons.org/licenses/by/4.0/).

## Citation

``` 

To cite mikRopML in publications, use:

  Topçuoğlu et al., (2021). mikropml: User-Friendly R Package for
  Supervised Machine Learning Pipelines. Journal of Open Source
  Software, 6(61), 3073, https://doi.org/10.21105/joss.03073

A BibTeX entry for LaTeX users is

  @Article{,
    title = {{mikropml}: User-Friendly R Package for Supervised Machine Learning Pipelines},
    author = {Begüm D. Topçuoğlu and Zena Lapp and Kelly L. Sovacool and Evan Snitkin and Jenna Wiens and Patrick D. Schloss},
    journal = {Journal of Open Source Software},
    year = {2021},
    month = {May},
    volume = {6},
    number = {61},
    pages = {3073},
    doi = {10.21105/joss.03073},
    url = {https://joss.theoj.org/papers/10.21105/joss.03073},
  }
```

## Why the name?

The word “mikrop” (pronounced “meek-ROPE”) is Turkish for “microbe”.
This package was originally implemented as a machine learning pipeline
for microbiome-based classification problems (see [Topçuoğlu *et al.*
2020](https://doi.org/10.1128/mBio.00434-20)). We realized that these
methods are applicable in many other fields too, but stuck with the name
because we like it\!
# development version 1.3.0

- Allow `kfold >= length(groups)` (#285, @kelly-sovacool).
    - When using the groups parameter, groups are kept together in cross-validation partitions when `kfold` <= the number of groups in the training set. Previously, an error was thrown if this condition was not met. Now, if there are not enough groups in the training set for groups to be kept together during CV, groups are allowed to be split up across CV partitions. 
- Report p-values for permutation feature importance (#288, @kelly-sovacool)

# mikropml 1.2.0

- New parameter `cross_val` added to `run_ml()` allows users to define their own custom cross-validation scheme (#278, @kelly-sovacool).
    - Also added a new parameter `calculate_performance`, which controls whether performance metrics are calculated (default: `TRUE`). Users may wish to skip performance calculations when training models with no cross-validation.
- New parameter `group_partitions` added to `run_ml()` allows users to control which groups should go to which partition of the train/test split (#281, @kelly-sovacool).
- Modified the `training_frac` parameter in `run_ml()` (#281, @kelly-sovacool).
    - By default, `training_frac` is a fraction between 0 and 1 that specifies how much of the dataset should be used in the training fraction of the train/test split.
    - Users can instead give `training_frac` a vector of indices that correspond to which rows of the dataset should go in the training fraction of the train/test split. This gives users direct control over exactly which observations are in the training fraction if desired.

# mikropml 1.1.1

- Fixed bugs related to grouping correlated features (#276, @kelly-sovacool).
    - Also, `group_correlated_features()` is now a user-facing function.

# mikropml 1.1.0

- New correlation method option for feature importance (#267, @courtneyarmour).
    - The default is still "spearman", and now you can use other methods supported by `stats::cor` with the `corr_method` parameter: `get_feature_importance(corr_method = "pearson")`
- There are now [video tutorials](https://www.youtube.com/playlist?list=PLmNrK_nkqBpKpzb9-vI4V7SdXC-jXEcmg) covering mikropml and other skills related to machine learning, created by @pschloss (#270).
- Fixed a bug where `preprocess_data()` converted the outcome column to a character vector (#273, @kelly-sovacool, @ecmaggioncalda).

# mikropml 1.0.0

- mikropml now has a logo created by @NLesniak!
- Made documentation improvements (#238, #231 @kelly-sovacool; #256 @BTopcuoglu).
- New option in `preprocess_data()`: `prefilter_threshold` (#240, @kelly-sovacool, @courtneyarmour).
    - Remove any features that appear in N=`prefilter_threshold` or fewer rows in the data.
    - Created function `remove_singleton_columns()` called by `preprocess_data()` to carry this out.
- New option in `get_feature_importance()`: `groups` (#246, @kelly-sovacool).
    - Provide custom groups of features to permute together during permutation importance.
    - `groups` is `NULL` by default; in this case, correlated features above `corr_thresh` are grouped together.
- `preprocess_data()` now replaces spaces in the outcome column with underscores (#247, @kelly-sovacool, @JonnyTran).
- Clarify in the intro vignette that we do not support multi-label outcomes. (#254, @zenalapp)
- Optional progress bar for `preprocess_data()` and `get_feature_importance()` using [the progressr package](https://github.com/HenrikBengtsson/progressr) (#257, @kelly-sovacool, @JonnyTran, @FedericoComoglio).
- The mikropml paper is soon to be published in [JOSS](https://joss.theoj.org/papers/10.21105/joss.03073)!

# mikropml 0.0.2

- Fixed a test failure on Solaris.
- Fixed multiple test failures with R 3.6.2 due to `stringsAsFactors` behavior.
- Made minor documentation improvements.
- Moved `rpart` from Suggests to Imports for consistency with other packages used during model training.

# mikropml 0.0.1

This is the first release version of mikropml! 🎉

- Added a `NEWS.md` file to track changes to the package.
- Major new functions:
    - `run_ml()`
    - `preprocess_data()`
    - `plot_model_performance()`
    - `plot_hp_performance()`
- Support for ML methods in `run_ml()`:
    - `glmnet`: logistic and linear regression
    - `rf`: random forest
    - `rpart2`: decision trees
    - `svmRadial`: support vector machines
    - `xgbTree`: gradient-boosted trees
- New vignettes:
    - [Introduction](http://www.schlosslab.org/mikropml/articles/introduction.html)
    - [Preprocess data](http://www.schlosslab.org/mikropml/articles/preprocess.html)
    - [Hyperparameter tuning](http://www.schlosslab.org/mikropml/articles/tuning.html)
    - [Parallel processing](http://www.schlosslab.org/mikropml/articles/parallel.html)
    - [The mikropml paper](http://www.schlosslab.org/mikropml/articles/paper.html)
## Test environments

* GitHub Actions (ubuntu-16.04): devel, release, oldrel
* GitHub Actions (windows): release
* GitHub Actions (macOS): release
* win-builder: devel

## R CMD check results

0 errors | 0 warnings | 1 note

```
  URL: https://anaconda.org/conda-forge/r-mikropml
    From: inst/doc/paper.html
          README.md
    Status: 400
    Message: Bad Request
```

I believe this is a spurious note as the URL works in my local browser.

## revdepcheck results

No reverse dependencies found.
# Getting help with mikropml

Thanks for using mikropml!
Before filing an issue, there are a few places to explore and pieces to put together to make the process as smooth as possible.

## Make a reprex

Start by making a minimal **repr**oducible **ex**ample using the  [reprex](https://reprex.tidyverse.org/) package. 
If you haven't heard of or used reprex before, you're in for a treat! 
Seriously, reprex will make all of your R-question-asking endeavors easier (which is a pretty insane ROI for the five to ten minutes it'll take you to learn what it's all about). 
For additional reprex pointers, check out the [Get help!](https://www.tidyverse.org/help/) section of the tidyverse site.

## Where to ask?

Armed with your reprex, the next step is to figure out [where to ask](https://www.tidyverse.org/help/#where-to-ask). 

*   If it's a question: start with [community.rstudio.com](https://community.rstudio.com/), and/or StackOverflow. There are more people there to answer questions.  

*   If it's a bug: you're in the right place, [file an issue](https://github.com/SchlossLab/mikropml/issues/new).  
  
*   If you're not sure: let the community help you figure it out! 
    If your problem _is_ a bug or a feature request, you can easily return here and report it. 

Before opening a new issue, be sure to [search issues and pull requests](https://github.com/SchlossLab/mikropml/issues) to make sure the bug hasn't been reported and/or already fixed in the development version. 
By default, the search will be pre-populated with `is:issue is:open`. 
You can [edit the qualifiers](https://help.github.com/articles/searching-issues-and-pull-requests/)  (e.g. `is:pr`, `is:closed`) as needed. 
For example, you'd simply remove `is:open` to search _all_ issues in the repo, open or closed.

## What happens next?

To be as efficient as possible, development of tidyverse packages tends to be very bursty, so you shouldn't worry if you don't get an immediate response.
Typically we don't look at a repo until a sufficient quantity of issues accumulates, then there’s a burst of intense activity as we focus our efforts. 
That makes development more efficient because it avoids expensive context switching between problems, at the cost of taking longer to get back to you. 
This process makes a good reprex particularly important because it might be multiple months between your initial report and when we start working on it. 
If we can’t reproduce the bug, we can’t fix it!
## Issues
- Resolves # .

## Change(s) made
-
-

## Checklist

(~Strikethrough~ any points that are not applicable.)

- [ ] Write unit tests for any new functionality or bug fixes.
- [ ] Update docs if there are any API changes:
  - [ ] roxygen comments
  - [ ] vignettes
- [ ] Update `NEWS.md` if this includes any user-facing changes. 
- [ ] The check workflow succeeds on your most recent commit. **This is always required before the PR can be merged.**
# Contributor Covenant Code of Conduct

This document was adapted from the [Tidyverse Code of Conduct](https://tidyverse.tidyverse.org/CODE_OF_CONDUCT.html).

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
# Contributing to mikropml

This document was adapted from the [Tidyverse Contributing guide](https://tidyverse.tidyverse.org/CONTRIBUTING.html).

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("SchlossLab/mikropml", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the mikropml project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
---

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

Brief description of the problem

```r
# insert reprex here
```
---
title: "mikropml: User-Friendly R Package for Supervised Machine Learning Pipelines"
output: 
  rmarkdown::html_vignette:
    keep_md: true
tags:
  - R
  - machine learning
  - regression
  - classification
  - decision trees
  - random forest
  - xgboost
  - support vector machines
  - microbiology
author: Begüm D. Topçuoğlu, Zena Lapp, Kelly L. Sovacool, Evan Snitkin, Jenna Wiens, Patrick D. Schloss
authors:
  - name: Begüm D. Topçuoğlu^[co-first author]
    orcid: 0000-0003-3140-537X
    affiliation: "3, 4"
  - name: Zena Lapp^[co-first author]
    orcid: 0000-0003-4674-2176
    affiliation: 1
  - name: Kelly L. Sovacool^[co-first author]
    orcid: 0000-0003-3283-829X
    affiliation: 1
  - name: Evan Snitkin
    orcid: 0000-0001-8409-278X
    affiliation: "3, 5"
  - name: Jenna Wiens
    orcid: 0000-0002-1057-7722
    affiliation: 2
  - name: Patrick D. Schloss^[corresponding author]
    orcid: 0000-0002-6935-4275
    affiliation: 3
affiliations:
  - name: Department of Computational Medicine & Bioinformatics, University of Michigan
    index: 1
  - name: Department of Electrical Engineering & Computer Science, University of Michigan
    index: 2
  - name: Department of Microbiology & Immunology, University of Michigan
    index: 3
  - name: Exploratory Science Center, Merck & Co., Inc., Cambridge, Massachusetts, USA.
    index: 4
  - name: Department of Internal Medicine/Division of Infectious Diseases, University of Michigan
    index: 5
date: 2020
bibliography: paper.bib
vignette: >
  %\VignetteIndexEntry{mikropml paper}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




# Summary

![](mikropml-logo.png){ width=120px }

Machine learning (ML) for classification and prediction based on a set of
features is used to make decisions in healthcare, economics, criminal justice
and more. However, implementing an ML pipeline including preprocessing, model
selection, and evaluation can be time-consuming, confusing, and difficult. Here,
we present [`mikropml`](http://www.schlosslab.org/mikropml/) (prononced
"meek-ROPE em el"), an easy-to-use R package that implements ML pipelines using
regression, support vector machines, decision trees, random forest, or
gradient-boosted trees. The package is available on
[GitHub](https://github.com/SchlossLab/mikropml/),
[CRAN](https://cran.r-project.org/package=mikropml), and
[conda](https://anaconda.org/conda-forge/r-mikropml).

# Statement of need

Most applications of machine learning (ML) require reproducible steps for data
pre-processing, cross-validation, testing, model evaluation, and often
interpretation of why the model makes particular predictions. Performing these
steps is important, as failure to implement them can result in incorrect and
misleading results [@teschendorff_avoiding_2019; @wiens_no_2019].

Supervised ML is widely used to recognize patterns in large datasets and to make
predictions about outcomes of interest. Several packages including `caret`
[@kuhn_building_2008] and `tidymodels` [@kuhn_tidymodels_2020] in R,
`scikitlearn` [@pedregosa_scikit-learn_2011] in Python, and the H2O `autoML`
platform [@h2o_platform] allow scientists to train ML models with a variety of
algorithms. While these packages provide the tools necessary for each ML step,
they do not implement a complete ML pipeline according to good practices in the
literature. This makes it difficult for practitioners new to ML to easily begin
to perform ML analyses.

To enable a broader range of researchers to apply ML to their problem domains,
we created [`mikropml`](https://github.com/SchlossLab/mikropml/), an easy-to-use
R package [@r_core_team_r_2020] that implements the ML pipeline created by
Topçuoğlu _et al._ [@topcuoglu_framework_2020] in a single function that returns a trained model,
model performance metrics and feature importance. `mikropml` leverages
the `caret` package to support several ML algorithms: linear regression,
logistic regression, support vector machines with a radial basis kernel,
decision trees, random forest, and gradient boosted trees. It incorporates good
practices in ML training, testing, and model evaluation
[@topcuoglu_framework_2020;@teschendorff_avoiding_2019]. Furthermore, it
provides data preprocessing steps based on the FIDDLE (FlexIble Data-Driven
pipeLinE) framework outlined in Tang _et al._ [@tang_democratizing_2020] and
post-training permutation importance steps to estimate the importance of each
feature in the models trained [@breiman_random_2001; @fisher_all_2018].

`mikropml` can be used as a starting point in the application of ML to datasets
from many different fields. It has already been applied to microbiome data to
categorize patients with colorectal cancer [@topcuoglu_framework_2020], to
identify differences in genomic and clinical features associated with bacterial
infections [@lapp_machine_2020], and to predict gender-based biases in academic
publishing [@hagan_women_2020].

# mikropml package

The `mikropml` package includes functionality to preprocess the data, train ML
models, evaluate model performance, and quantify feature importance (Figure 1).
We also provide
[vignettes](http://www.schlosslab.org/mikropml/articles/index.html) and an
[example Snakemake
workflow](https://github.com/SchlossLab/mikropml-snakemake-workflow)
[@koster_snakemakescalable_2012] to showcase how to run an ideal ML pipeline
with multiple different train/test data splits. The results can be visualized
using helper functions that use `ggplot2` [@wickham_ggplot2_2016].

While mikropml allows users to get started quickly and facilitates
reproducibility, it is not a replacement for understanding the ML workflow which
is still necessary when interpreting results [@pollard_turning_2019]. To
facilitate understanding and enable one to tailor the code to their application,
we have heavily commented the code and have provided supporting documentation
which can be read [online](http://www.schlosslab.org/mikropml/).

## Preprocessing data

We provide the function `preprocess_data()` to preprocess features using several
different functions from the `caret` package. `preprocess_data()` takes
continuous and categorical data, re-factors categorical data into binary
features, and provides options to normalize continuous data, remove features
with near-zero variance, and keep only one instance of perfectly correlated
features. We set the default options based on those implemented in FIDDLE
[@tang_democratizing_2020]. More details on how to use `preprocess_data()` can
be found in the accompanying
[vignette](http://www.schlosslab.org/mikropml/articles/preprocess.html).

## Running ML

The main function in mikropml, `run_ml()`, minimally takes in the model choice
and a data frame with an outcome column and feature columns. For model choice,
`mikropml` currently supports logistic and linear regression [`glmnet`:
@friedman_regularization_2010], support vector machines with a radial basis
kernel [`kernlab`: @karatzoglou_kernlab_2004], decision trees [`rpart`:
@therneau_rpart_2019], random forest [`randomForest`: @liaw_classication_2002],
and gradient-boosted trees [`xgboost`:  @chen_xgboost_2020]. `run_ml()` randomly
splits the data into train and test sets while maintaining the distribution of
the outcomes found in the full dataset. It also provides the option to split the
data into train and test sets based on categorical variables (e.g. batch,
geographic location, etc.). `mikropml` uses the `caret` package
[@kuhn_building_2008] to train and evaluate the models, and optionally
quantifies feature importance. The output includes the best model built based on
tuning hyperparameters in an internal and repeated cross-validation step, model
evaluation metrics, and optional feature importances. Feature importances are
calculated using a permutation test, which breaks the relationship between the
feature and the true outcome in the test data, and measures the change in model
performance. This provides an intuitive metric of how individual features
influence model performance and is comparable across model types, which is
particularly useful for model interpretation [@topcuoglu_framework_2020]. Our
[introductory
vignette](http://www.schlosslab.org/mikropml/articles/introduction.html)
contains a comprehensive tutorial on how to use `run_ml()`.

![mikropml pipeline](mikRopML-pipeline.png){width=100%}

## Ideal workflow for running mikropml with many different train/test splits

To investigate the variation in model performance depending on the train and
test set used [@topcuoglu_framework_2020; @lapp_machine_2020], we provide
examples of how to `run_ml()` many times with different train/test splits and
how to get summary information about model performance on [a local
computer](http://www.schlosslab.org/mikropml/articles/parallel.html) or on a
high-performance computing cluster using a [Snakemake
workflow](https://github.com/SchlossLab/mikropml-snakemake-workflow).

## Tuning & visualization

One particularly important aspect of ML is hyperparameter tuning. We provide a
reasonable range of default hyperparameters for each model type. However
practitioners should explore whether that range is appropriate for their data,
or if they should customize the hyperparameter range. Therefore, we provide a
function `plot_hp_performance()` to plot the cross-validation performance metric
of a single model or models built using different train/test splits. This helps
evaluate if the hyperparameter range is being searched exhaustively and allows
the user to pick the ideal set. We also provide summary plots of test
performance metrics for the many train/test splits with different models using
`plot_model_performance()`. Examples are described in the accompanying [vignette
on hyperparameter
tuning](http://www.schlosslab.org/mikropml/articles/tuning.html).

## Dependencies

mikropml is written in R [@r_core_team_r_2020] and depends on several packages:
`dplyr` [@wickham_dplyr_2020], `rlang` [@henry_rlang_2020] and `caret`
[@kuhn_building_2008]. The ML algorithms supported by `mikropml` require:
`glmnet` [@friedman_regularization_2010], `e1071` [@meyer_e1071_2020], and
`MLmetrics` [@yan_mlmetrics_2016] for logistic regression, `rpart2`
[@therneau_rpart_2019] for decision trees, `randomForest`
[@liaw_classication_2002] for random forest, `xgboost` [@chen_xgboost_2020] for
xgboost, and `kernlab` [@karatzoglou_kernlab_2004] for support vector machines.
We also allow for parallelization of cross-validation and other steps using the
`foreach`, `doFuture`, `future.apply`, and `future` packages
[@bengtsson_futureapply_2020]. Finally, we use `ggplot2` for plotting
[@wickham_ggplot2_2016].

# Acknowledgments

We thank members of the Schloss Lab who participated in code clubs related to
the initial development of the pipeline, made documentation improvements, and
provided general feedback.
We also thank Nick Lesniak for designing the mikropml logo.

We thank the US Research Software Sustainability Institute (NSF #1743188) for
providing training to KLS at the Winter School in Research Software Engineering.

# Funding

Salary support for PDS came from NIH grant 1R01CA215574. KLS received support
from the NIH Training Program in Bioinformatics (T32 GM070449). ZL received
support from the National Science Foundation Graduate Research Fellowship
Program under Grant No. DGE 1256260. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the authors and do not
necessarily reflect the views of the National Science Foundation.

# Author contributions

BDT, ZL, and KLS contributed equally. Author order among the co-first authors
was determined by time since joining the project.

BDT, ZL, and KLS conceptualized the study and wrote the code. KLS structured the
code in R package form. BDT, ZL, JW, and PDS developed methodology. PDS, ES, and
JW supervised the project. BDT, ZL, and KLS wrote the original draft. All
authors reviewed and edited the manuscript.

# Conflicts of interest

None.

# References
We thank the reviewers for their thoughtful feedback. We have responded to all reviewer comments below. Feel free to let us know if you have additional questions.

## @JonnyTran reviews

> When testing `preprocess_data()` and `run_ml()` on a larger dataset with 1917 rows × 14426 columns, it would take an indeterminately long time. It should have a progress bar indicating the iteration #, and the estimated duration for each iteration, similar to Python's `tqdm`.

We opened and resolved issue [#249](https://github.com/SchlossLab/mikropml/issues/249) related to this suggestion.
We have implemented this feature for `preprocess_data()` and the feature importance part of `run_ml()`.
Unfortunately, we cannot do it for model training as we use `caret::train()`, which does not support having a progress bar.

> Additionally, it can also provide live updates of the loss and other metrics as the model is being trained, to help the user stop the run if it is not satisfactory.

Another limitation of using `caret::train()` to train the model is that there is no option to display updates on relevant metrics as the model is being trained.
As much of our package is based on using this function, we unfortunately cannot easily implement this suggestion.
(See issue [#250](https://github.com/SchlossLab/mikropml/issues/250))

> I also encountered a bug when running `run_ml()` on a dataset with class name containing spaces in the label column. I can mitigate the bug by replacing all `" "`'s with `"_"`'s.

Thank you for bringing this to our attention. We fixed this issue ([#244](https://github.com/SchlossLab/mikropml/issues/244)).

> Currently, it doesn't seem to be able to support multi-label classifications.

Unfortunately `caret` does not support this functionality so we believe it's outside the scope of this package (see issue [#252](https://github.com/SchlossLab/mikropml/issues/252)).

> The `run_ml()` function currently implements 5 off-the-shelve ML algorithms, while providing 12 other parameters for training criteria, hyperparameters, feature importance, etc.. If in the future it would support more algorithms, custom metrics, or training parameters, I'd imagine there'll be limitations imposed by the function arguments. I'd suggest the function to take in 3 objects, e.g. `run_ml(dataset, model, metrics, [args])`, where a metrics object can allow the user select standard metrics or define their own metric functions given the model output and true labels.

This is a very good point that we have responded to in issue [#251](https://github.com/SchlossLab/mikropml/issues/251).
We chose to structure `run_ml` in a way that is easier for our users, rather than extremely flexible.
If users want a more customizable experience, especially related to what models they want to use, we expect them to use other packages not developed for machine learning novices.
Furthermore, `run_ml()` already takes these objects:

1. `dataset`: The input dataset.
1.  `method`: The ML model to be used. While we only officially support 5 models (as these are the ones we tested), all of the models supported by `caret` (https://topepo.github.io/caret/available-models.html) should work in our package. If `caret` supports additional models in the future, these should also work in `mikropml`. We realize that the model options are not as generalizable as e.g. PyTorch, since users must choose from options supported by `caret`. However, our code heavily relies on `caret` to perform the underlying model training. Additionally, as `mikropml` is oriented toward beginner practitioners, we believe that it does not need to provide the option to include custom models.
1. `perf_metric_function` and `perf_metric_name`: The performance metric to be used. We chose sensible defaults, but the user can provide their own performance metrics if they would like.
1. `hyperparameters`: The values of hyperparameters in the model that the user would like to tune.

## @FedericoComoglio reviews

> I agree with @JonnyTran 's suggestion to display time to completion estimates and relevant metrics (e.g. loss or chosen metrics for current iteration) during training. The average user will find this helpful. This would also be the case when parallelizing `run_ml()` execution.

See response above.

> In a number of places across the vignette, "machine learning" (e.g. "run machine learning" in the "Preprocessing data" vignette) is used in place of a likely more appropriate "model training"/"train model"

Thank you for this suggestion. We have changed the wording in several places (issue [#253](https://github.com/SchlossLab/mikropml/issues/253) opened and resolved).

> For training random forest models, it is unfortunate that the number of trees can not be tuned. However, this is a consequence of design choices that pertain to the caret package and hence not directly related to this work.

We wholeheartedly agree.
