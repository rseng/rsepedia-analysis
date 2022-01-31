
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
---
output:
    github_document:
      html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

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

An interface to build machine learning models for classification and regression
problems. `mikropml` implements the ML pipeline described by [Topçuoğlu _et al._
(2020)](https://doi.org/doi:10.1128/mBio.00434-20) with reasonable default
options for data preprocessing, hyperparameter tuning, cross-validation,
testing, model evaluation, and interpretation steps. See the
[website](http://www.schlosslab.org/mikropml/) for more information,
documentation, and examples.

## Installation

You can install the latest release from
[CRAN](https://cran.r-project.org/package=mikropml):
```{r install_cran, eval = FALSE}
install.packages('mikropml')
```

or the development version from 
[GitHub](https://github.com/SchlossLab/mikRopML):

```{r install_github, eval = FALSE}
# install.packages("devtools")
devtools::install_github("SchlossLab/mikropml")
```

or install from a terminal using
[conda](https://docs.conda.io/projects/conda/en/latest/index.html):

```{bash conda, eval = FALSE}
conda install -c conda-forge r-mikropml
```


### Dependencies

```{r deps, echo = FALSE, message = FALSE, warning = FALSE}
library(dplyr)
description <- utils::packageDescription('mikropml', 
                                         fields = c('Imports', 'Suggests'))
deps <- lapply(names(description), 
               function (x) {
                 paste0('- ', x, ': ', 
                        description[[x]] %>% 
                          gsub("\n", " ", .))}
               ) %>% 
  unlist() %>% 
  paste(., collapse = '\n')
```

`r deps`

## Usage

Check out the [introductory
vignette](http://www.schlosslab.org/mikropml/articles/introduction.html) for a
quick start tutorial. For a more in-depth discussion, read [all the
vignettes](http://www.schlosslab.org/mikropml/articles/index.html) and/or take a
look at the [reference
documentation](http://www.schlosslab.org/mikropml/reference/index.html). 

You can watch the Riffomonas Project series of 
[video tutorials](https://www.youtube.com/playlist?list=PLmNrK_nkqBpKpzb9-vI4V7SdXC-jXEcmg) 
covering mikropml and other skills related to machine learning.

We also
provide an [example Snakemake
workflow](https://github.com/SchlossLab/mikropml-snakemake-workflow) for running
`mikropml` on an HPC.

## Help & Contributing

If you come across a bug, [open an
issue](https://github.com/SchlossLab/mikropml/issues) and include a [minimal
reproducible example](https://www.tidyverse.org/help/).

If you'd like to contribute, see our guidelines
[here](http://www.schlosslab.org/mikropml/CONTRIBUTING.html).

## Code of Conduct

Please note that the mikropml project is released with a [Contributor Code of
Conduct](http://www.schlosslab.org/mikropml/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

## License

The mikropml package is licensed under 
[the MIT license](https://github.com/SchlossLab/mikropml/blob/main/LICENSE.md).
Text and images included in this repository, including the mikropml logo, 
are licensed under the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).

## Citation

```{r cite, echo = FALSE, comment = ''}
citation(package = 'mikropml')
```

## Why the name?

The word "mikrop" (pronounced "meek-ROPE") is Turkish for "microbe". This
package was originally implemented as a machine learning pipeline for
microbiome-based classification problems (see [Topçuoğlu _et al._
2020](https://doi.org/10.1128/mBio.00434-20)). We realized that these methods
are applicable in many other fields too, but stuck with the name because we like
it!
---
title: "Parallel processing"
author: "Kelly L. Sovacool"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Parallel processing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
#NOT_CRAN <- TRUE
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

```{r setup}
library(mikropml)
library(dplyr)
```

## Speed up single runs

By default, `preprocess_data()` and `run_ml()` use only one process in series.
If you'd like to parallelize various steps of the pipeline to make them run
faster, install `foreach`, `future`, `future.apply`, and `doFuture`. Then,
register a future plan prior to calling `preprocess_data()` and `run_ml()`:

```{r register, eval = FALSE}
doFuture::registerDoFuture()
future::plan(future::multicore, workers = 2)
```

Above, we used the `multicore` plan to split the work across 2 cores. See the
[`future`
documentation](https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html)
for more about picking the best plan for your use case. Notably, `multicore`
does not work inside RStudio or on Windows; you will need to use `multisession`
instead in those cases.

After registering a future plan, you can call `preprocess_data()` and `run_ml()`
as usual, and they will run certain tasks in parallel.

```{r run_single}
otu_data_preproc <- preprocess_data(otu_mini_bin, 'dx')$dat_transformed
result1 <- run_ml(otu_data_preproc, 'glmnet')
```

## Call `run_ml()` multiple times in parallel in R

You can use functions from the `future.apply` package to call `run_ml()`
multiple times in parallel with different parameters. You will first need to run
`future::plan()` as above if you haven't already. Then, call `run_ml()` with
multiple seeds using `future_lapply()`:

```{r multi_seeds}
# NOTE: use more seeds for real-world data
results_multi <- future.apply::future_lapply(seq(100, 102), function(seed) {
  run_ml(otu_data_preproc, 'glmnet', seed = seed)
  }, future.seed = TRUE)
```

Each call to `run_ml()` with a different seed uses a different random split of
the data into training and testing sets. Since we are using seeds, we must set
`future.seed` to `TRUE` (see the [`future.apply`
documentation](https://cran.r-project.org/web/packages/future.apply/future.apply.pdf)
and [this blog
post](https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/)
for details on parallel-safe random seeds). This example uses only a few seeds
for speed and simplicity, but for real data we recommend using many more seeds
to get a better estimate of model performance.

In these examples, we used functions from the `future.apply` package to
`run_ml()` in parallel, but you can accomplish the same thing with parallel
versions of the `purrr::map()` functions using the `furrr` package (e.g.
`furrr::future_map_dfr()`).

Extract the performance results and combine into one dataframe for all seeds:

```{r bind_multi_seeds}
perf_df <- future.apply::future_lapply(results_multi, 
                                       function(result) {
                                         result[['performance']] %>% 
                                           select(cv_metric_AUC, AUC, method)
                                         },
                                       future.seed = TRUE) %>% 
  dplyr::bind_rows()
perf_df
```

### Multiple ML methods

You may also wish to compare performance for different ML methods. `mapply()`
can iterate over multiple lists or vectors, and `future_mapply()` works the same
way:

```{r multi_methods_seeds}
# NOTE: use more seeds for real-world data
param_grid <- expand.grid(seeds = seq(100, 102),
                          methods = c('glmnet', 'rf'))
results_mtx <- future.apply::future_mapply(
    function(seed, method) {
      run_ml(otu_data_preproc, method, seed = seed)
      },
    param_grid$seeds,
    param_grid$methods %>% as.character(),
    future.seed = TRUE
  )
```

Extract and combine the performance results for all seeds and methods:

```{r bind_multi_methods}
perf_df2 <- lapply(results_mtx['performance',], 
                   function(x) {
                     x %>% select(cv_metric_AUC, AUC, method)
                   }) %>% 
  dplyr::bind_rows()
perf_df2
```

Visualize the performance results (`ggplot2` is required):

```{r plot_perf}
perf_boxplot <- plot_model_performance(perf_df2)
perf_boxplot
```

`plot_model_performance()` returns a ggplot2 object. 
You can add layers to customize the plot:

```{r customize_plot}
perf_boxplot +
   theme_classic() +
   scale_color_brewer(palette = "Dark2") +
   coord_flip()
```

You can also create your own plots however you like using the performance
results.

## Live progress updates

`preprocess_data()` and `get_feature_importance()` support reporting live
progress updates using the `progressr` package. The format is up to you, but we
recommend using a progress bar like this:

```{r progress, eval = FALSE}
# optionally, specify the progress bar format with the `progress` package.
progressr::handlers(progressr::handler_progress(
    format = ":message :bar :percent | elapsed: :elapsed | eta: :eta",
    clear = FALSE,
    show_after = 0))
# tell progressr to always report progress in any functions that use it.
# set this to FALSE to turn it back off again.
progressr::handlers(global = TRUE)

# run your code and watch the live progress updates.
dat <- preprocess_data(otu_mini_bin, 'dx')$dat_transformed
#> Using 'dx' as the outcome column.
#> preprocessing ========================>-------  78% | elapsed:  1s | eta:  0s
results <- run_ml(dat, "glmnet", kfold = 2, cv_times = 2,
                  find_feature_importance = TRUE)
#> Using 'dx' as the outcome column.
#> Training the model...
#> Training complete.
#> Feature importance =========================== 100% | elapsed: 37s | eta:  0s
```

Note that some future backends support "near-live" progress updates, meaning the
progress may not be reported immediately when parallel processing with futures.
Read more on that [in the `progressr`
vignette](https://progressr.futureverse.org/articles/progressr-intro.html#near-live-versus-buffered-progress-updates-with-futures).
For more on `progressr` and how to customize the format of progress updates, see
the [`progressr` docs](https://progressr.futureverse.org/).

## Parallelizing with Snakemake

When parallelizing multiple calls to `run_ml()` in R as in the examples above,
all of the results objects are held in memory. This isn't a big deal for a small
dataset run with only a few seeds. However, for large datasets run in parallel
with, say, 100 seeds (recommended), you may run into problems trying to store
all of those objects in memory at once. One solution is to write the results
files of each `run_ml()` call, then concatenate them at the end. We show one way
to accomplish this with Snakemake in [an example Snakemake workflow
here](https://github.com/SchlossLab/mikropml-snakemake-workflow).
---
title: "Introduction to mikropml"
author: "Zena Lapp"
output: rmarkdown::html_vignette
bibliography: paper.bib
vignette: >
  %\VignetteIndexEntry{Introduction to mikropml}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The goal of `mikropml` is to make supervised machine learning (ML) easy for you
to run while implementing good practices for machine learning pipelines. All you
need to run the ML pipeline is one function: `run_ml()`. We've selected sensible default
arguments related to good practices [@topcuoglu_framework_2020;
@tang_democratizing_2020], but we allow you to change those arguments to tailor
`run_ml()` to the needs of your data.

This document takes you through all of the `run_ml()` inputs, both required and
optional, as well as the outputs.

In summary, you provide:

- A dataset with an outcome column and feature columns (rows are samples; unfortunately we do not support multi-label classification)
- Model choice (i.e. method)

And the function outputs:

- The trained model
- Model performance metrics
- (Optional) feature importance metrics

# It's running so slow!

Since I assume a lot of you won't read this entire vignette, I'm going to say
this at the beginning. If the `run_ml()` function is running super slow, you
should consider parallelizing. See `vignette("parallel")` for examples.

# Understanding the inputs

## The input data

The input data to `run_ml()` is a dataframe where each row is a sample or
observation. One column (assumed to be the first) is the outcome of interest,
and all of the other columns are the features. We package `otu_mini_bin` as a
small example dataset with `mikropml`.

```{r}
#install.packages("devtools")
#devtools::install_github("SchlossLab/mikropml")
library(mikropml)
head(otu_mini_bin)
```

Here, `dx` is the outcome column (normal or cancer), and there are 10 features
(`Otu00001` through `Otu00010`). Because there are only 2 outcomes, we will be
performing binary classification in the majority of the examples below. At the
bottom, we will also briefly provide examples of multi-class and continuous
outcomes. As you'll see, you run them in the same way as for binary
classification!

The feature columns are the amount of each [Operational Taxonomic Unit
(OTU)](https://en.wikipedia.org/wiki/Operational_taxonomic_unit) in microbiome
samples from patients with cancer and without cancer. The goal is to predict
`dx`, which stands for diagnosis. This diagnosis can be cancer or not based on
an individual's microbiome. No need to understand exactly what that means, but
if you're interested you can read more about it from the original paper
[@topcuoglu_framework_2020].

For real machine learning applications you'll need to use more features, but for
the purposes of this vignette we'll stick with this example dataset so
everything runs faster.

## The methods we support

All of the methods we use are supported by a great ML wrapper package
[`caret`](https://topepo.github.io/caret/), which we use to train our machine
learning models.

The methods we have tested (and their backend packages) are:

- Logistic/multiclass/linear regression (`"glmnet"`)
- Random forest (`"rf"`)
- Decision tree (`"rpart2"`)
- Support vector machine with a radial basis kernel (`"svmRadial"`)
- xgboost (`"xgbTree"`)

For documentation on these methods, as well as many others, you can look at the
[available models](https://topepo.github.io/caret/available-models.html) (or see
[here](https://topepo.github.io/caret/train-models-by-tag.html) for a list by
tag). While we have not vetted the other models used by `caret`, our function is
general enough that others might work. While we can't promise that we can help
with other models, feel free to [open an issue on
GitHub](https://github.com/SchlossLab/mikropml/issues) if you have questions
about other models and we _might_ be able to help.

We will first focus on `glmnet`, which is our default implementation of
L2-regularized logistic regression. Then we will cover a few other examples
towards the end.

# Before running ML

Before you execute `run_ml()`, you should consider preprocessing your data,
either on your own or with the `preprocess_data()` function. You can learn more
about this in the preprocessing vignette: `vignette("preprocess")`.

# The simplest way to `run_ml()`

As mentioned above, the minimal input is your dataset (`dataset`) and the
machine learning model you want to use (`method`).

You may also want to provide:

- The outcome column name. By default `run_ml()` will pick the first column, but it's best practice to specify the column name explicitly.
- A seed so that the results will be reproducible, and so that you get the same results as those you see here (i.e have the same train/test split).

Say we want to use logistic regression, then the method we will use is `glmnet`.
To do so, run the ML pipeline with:

```{r, eval = FALSE}
results <- run_ml(otu_mini_bin,
                  'glmnet',
                  outcome_colname = 'dx',
                  seed = 2019)
```
```{r, echo = FALSE}
# reduce vignette runtime by using precomputed results
results <- otu_mini_bin_results_glmnet
```

You'll notice a few things:

1. It takes a little while to run. This is because of some of the parameters we use.
1. There is a message stating that 'dx' is being used as the outcome column. This is what we want, but it's a nice sanity check!
1. There was a warning. Don't worry about this warning right now - it just means that some of the hyperparameters aren't a good fit - but if you're interested in learning more, see `vignette("tuning")`.

Now, let's dig into the output a bit.
The results is a list of 4 things:

```{r}
names(results)
```

`trained_model` is the trained model from `caret`.
There is a bunch of info in this that we won't get into,
because you can learn more from the `caret::train()` documentation.

```{r}
names(results$trained_model)
```

`test_data` is the partition of the dataset that was used for testing. In
machine learning, it's always important to have a held-out test dataset that is
not used in the training stage. In this pipeline we do that using `run_ml()`
where we split your data into training and testing sets. The training data are
used to build the model (e.g. tune hyperparameters, learn the data) and the test
data are used to evaluate how well the model performs.

```{r}
head(results$test_data)
```

`performance` is a dataframe of (mainly) performance metrics (1 column for
cross-validation performance metric, several for test performance metrics, and 2
columns at the end with ML method and seed):

```{r}
results$performance
```

When using logistic regression for binary classification, area under the
receiver-operator characteristic curve (AUC) is a useful metric to evaluate
model performance. Because of that, it's the default that we use for `mikropml`.
However, it is crucial to evaluate your model performance using multiple
metrics. Below you can find more information about other performance metrics and
how to use them in our package.

`cv_metric_AUC` is the AUC for the cross-validation folds for the training data.
This gives us a sense of how well the model performs on the training data.

Most of the other columns are performance metrics for the test data — the data
that wasn't used to build the model. Here, you can see that the AUC for the test
data is not much above 0.5, suggesting that this model does not predict much
better than chance, and that the model is overfit because the cross-validation
AUC (`cv_metric_AUC`, measured during training) is much higher than the testing
AUC. This isn't too surprising since we're using so few features with this
example dataset, so don't be discouraged. The default option also provides a
number of other performance metrics that you might be interested in, including
area under the precision-recall curve (prAUC).

The last columns of `results$performance` are the method and seed (if you set
one) to help with combining results from multiple runs (see
`vignette("parallel")`).

`feature_importance` has information about feature importance values if
`find_feature_importance = TRUE` (the default is `FALSE`). Since we used the
defaults, there's nothing here:

```{r}
results$feature_importance
```

# Customizing parameters

There are a few arguments that allow you to change how you execute `run_ml()`.
We've chosen reasonable defaults for you, but we encourage you to change these
if you think something else would be better for your data.

## Changing `kfold`, `cv_times`, and `training_frac`

- `kfold`: The number of folds to run for cross-validation (default: 5).
- `cv_times`: The number of times to run repeated cross-validation (default: 100).
- `training_frac`: The fraction of data for the training set (default: 0.8). The rest of the data is used for testing.

Here's an example where we change some of the default parameters:

```{r}
results_custom <- run_ml(otu_mini_bin,
                         'glmnet',
                         kfold = 2,
                         cv_times = 5,
                         training_frac = 0.5,
                         seed = 2019)
```

You might have noticed that this one ran faster — that's because we reduced `kfold` and `cv_times`.
This is okay for testing things out and may even be necessary for smaller datasets. 
But in general it may be better to have larger numbers for these parameters; 
we think the defaults are a good starting point [@topcuoglu_framework_2020].

### Custom training indices

When `training_frac` is a fraction between 0 and 1, a random sample of
observations in the dataset are chosen for the training set to satisfy the
`training_frac`. However, in some cases you might wish to control exactly which
observations are in the training set. You can instead assign `training_frac` a
vector of indices that correspond to which rows of the dataset should go in the
training set (all remaining sequences will go in the testing set).

```{r custom_train_indices, warning=FALSE}
n_obs <- otu_mini_bin %>% nrow()
training_size <- 0.8 * n_obs
training_rows <- sample(n_obs, training_size)
results_custom_train <- run_ml(otu_mini_bin,
                               'glmnet',
                               kfold = 2,
                               cv_times = 5,
                               training_frac = training_rows,
                               seed = 2019
                               )
```


## Changing the performance metric

There are two arguments that allow you to change what performance metric to use
for model evaluation, and what performance metrics to calculate using the test
data.

`perf_metric_function` is the function used to calculate the performance
metrics.

The default for classification is `caret::multiClassSummary()` and the default
for regression is `caret::defaultSummary()`. We'd suggest not changing this
unless you really know what you're doing.

`perf_metric_name` is the column name from the output of `perf_metric_function`.
We chose reasonable defaults (AUC for binary, logLoss for multiclass, and RMSE
for continuous), but the default functions calculate a bunch of different
performance metrics, so you can choose a different one if you'd like.

The default performance metrics available for classification are:

```{r, echo=FALSE}
# TODO: can we get these programmatically somehow instead of hard-coding them?
c("logLoss", "AUC", "prAUC", "Accuracy", "Kappa", "Mean_F1", "Mean_Sensitivity", "Mean_Specificity", "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", "Mean_Precision", "Mean_Recall", "Mean_Detection_Rate", "Mean_Balanced_Accuracy")
```

The default performance metrics available for regression are:

```{r, echo=FALSE}
c("RMSE", "Rsquared", "MAE")
```

Here's an example using prAUC instead of AUC:

```{r}
results_pr <- run_ml(otu_mini_bin, 
                     'glmnet', 
                     cv_times = 5, 
                     perf_metric_name = 'prAUC', 
                     seed = 2019)
```

You'll see that the cross-validation metric is prAUC, instead of the default AUC:

```{r}
results_pr$performance
```

## Using groups

The optional `groups` is a vector of groups to keep together when splitting the
data into train and test sets and for cross-validation. Sometimes it's
important to split up the data based on a grouping instead of just randomly.
This allows you to control for similarities within groups that you don't want to
skew your predictions (i.e. batch effects). For example, with biological data
you may have samples collected from multiple hospitals, and you might like to
keep observations from the same hospital in the same partition.

Here's an example where we split the data into train/test sets based on groups:

```{r custom_groups, warning=FALSE}
# make random groups
set.seed(2019)
grps <- sample(LETTERS[1:8], nrow(otu_mini_bin), replace=TRUE)
results_grp <- run_ml(otu_mini_bin, 
                      'glmnet', 
                      cv_times = 2, 
                      training_frac = 0.8, 
                      groups = grps, 
                      seed = 2019)
```

The one difference here is `run_ml()` will report how much of the data is in the
training set if you run the above code chunk. This can be a little finicky
depending on how many samples and groups you have. This is because it won't be
exactly what you specify with `training_frac`, since you have to include all of
one group in either the training set _or_ the test set.

### Controling how groups are assigned to partitions

When you use the `groups` parameter as above, by default `run_ml()` will assume
that you want all of the observations from each group to be placed in the same
partition of the train/test split. This makes sense when you want to use groups
to control for batch effects. However, in some cases you might prefer to control
exactly which groups end up in which partition, and you might even be okay with
some observations from the same group being assigned to different partitions.

For example, say you want groups A and B to be used for training, C and D for
testing, and you don't have a preference for what happens to the other groups.
You can give the `group_partitions` parameter a named list to specify which
groups should go in the training set and which should go in the testing set.

```{r group_partitions, warning=FALSE}
results_grp_part <- run_ml(otu_mini_bin, 
                      'glmnet', 
                      cv_times = 2, 
                      training_frac = 0.8, 
                      groups = grps, 
                      group_partitions = list(train = c('A', 'B'),
                                              test = c('C', 'D')
                                              ),
                      seed = 2019)
```

In the above case, all observations from A & B will be used for training, all
from C & D will be used for testing, and the remaining groups will be randomly
assigned to one or the other to satisfy the `training_frac` as closely as
possible.

In another scenario, maybe you want only groups A through F to be used for
training, but you also want to allow other observations not selected for
training from A through F to be used for testing:

```{r only_group_A_train, warning = FALSE}
results_grp_trainA <- run_ml(otu_mini_bin, 
                      'glmnet', 
                      cv_times = 2, 
                      kfold = 2,
                      training_frac = 0.5, 
                      groups = grps, 
                      group_partitions = list(train = c("A", "B", "C", "D", "E", "F"),
                                              test = c("A", "B", "C", "D", "E", "F", "G", "H")
                                              ),
                      seed = 2019)
```

If you need even more control than this, take a look at 
[setting custom training indices](#custom-training-indices). 
You might also prefer to  provide your own train control scheme with the
`cross_val` parameter in `run_ml()`.

# Finding feature importance

To find which features are contributing to predictive power, you can use
`find_feature_importance = TRUE`. How we use permutation importance to determine
feature importance is described in [@topcuoglu_framework_2020]. Briefly, it
permutes each of the features individually (or correlated ones together) and
evaluates how much the performance metric decreases. The more performance
decreases when the feature is randomly shuffled, the more important that feature
is. The default is `FALSE` because it takes a while to run and is only useful if
you want to know what features are important in predicting your outcome.

Let's look at some feature importance results:

```{r, eval = FALSE}
results_imp <- run_ml(otu_mini_bin,
  "rf",
  outcome_colname = "dx",
  find_feature_importance = TRUE,
  seed = 2019
)
```
```{r, echo = FALSE}
results_imp <- otu_mini_bin_results_rf
```

Now, we can check out the feature importances:

```{r}
results_imp$feature_importance
```

There are several columns:

1. `perf_metric`: The performance value of the permuted feature.
1. `perf_metric_diff`: The difference between the performance for the actual and permuted data (i.e. test performance minus permuted performance). Features with a larger `perf_metric_diff` are more important.
1. `pvalue`: the probability of obtaining the actual performance value under the null hypothesis.
1. `names`: The feature that was permuted.
1. `method`: The ML method used.
1. `perf_metric_name`: The peformance metric used.
1. `seed`: The seed (if set).

As you can see here, the differences are negligible (close to zero), which makes
sense since our model isn't great. If you're interested in feature importance,
it's especially useful to run multiple different train/test splits, as shown in
our 
[example snakemake workflow](https://github.com/SchlossLab/mikropml-snakemake-workflow/).

You can also choose to permute correlated features together using `corr_thresh`
(default: 1). Any features that are above the correlation threshold are permuted
together; i.e. perfectly correlated features are permuted together when using
the default value.

```{r}
results_imp_corr <- run_ml(otu_mini_bin,
                           'glmnet',
                           cv_times = 5,
                           find_feature_importance = TRUE,
                           corr_thresh = 0.2,
                           seed = 2019)
results_imp_corr$feature_importance
```

You can see which features were permuted together in the `names` column. Here
all 3 features were permuted together (which doesn't really make sense, but it's
just an example).

If you previously executed `run_ml()` without feature importance but now wish
to find feature importance after the fact, see the example code in the
`get_feature_importance()` documentation.

`get_feature_importance()` can show a live progress bar, see
`vignette("parallel")` for examples.

# Tuning hyperparameters (using the `hyperparameter` argument)

This is important, so we have a whole vignette about them. The bottom line is we
provide default hyperparameters that you can start with, but it's important to
tune your hyperparameters. For more information about what the default
hyperparameters are, and how to tune hyperparameters, see `vignette("tuning")`.

# Other models

Here are examples of how to train and evaluate other models.
The output for all of them is very similar, so we won't go into those details.

## Random forest

```{r, eval = FALSE}
results_rf <- run_ml(otu_mini_bin,
                     'rf',
                     cv_times = 5,
                     seed = 2019)
```

You can also change the number of trees to use for random forest 
(`ntree`; default: 1000). 
This can't be tuned using `rf` package implementation of random forest. 
Please refer to `caret` documentation if you are interested in 
other packages with random forest implementations.

```{r, eval = FALSE}
results_rf_nt <- run_ml(otu_mini_bin,
                        'rf',
                        cv_times = 5,
                        ntree = 10,
                        seed = 2019)
```

## Decision tree

```{r, eval = FALSE}
results_dt <- run_ml(otu_mini_bin,
                     'rpart2',
                     cv_times = 5,
                     seed = 2019)
```

## SVM

```{r, eval = FALSE}
results_svm <- run_ml(otu_mini_bin,
                      'svmRadial',
                      cv_times = 5,
                      seed = 2019)
```

If you get a message "maximum number of iterations reached", see [this issue](https://github.com/topepo/caret/issues/425) in caret.

# Other data

## Multiclass data

We provide `otu_mini_multi` with a multiclass outcome (three or more outcomes): 

```{r}
otu_mini_multi %>% dplyr::pull('dx') %>% unique()
```


Here's an example of running multiclass data:

```{r, eval = FALSE}
results_multi <- run_ml(otu_mini_multi,
                        outcome_colname = "dx",
                        seed = 2019
)
```
```{r, echo = FALSE}
results_multi <- otu_mini_multi_results_glmnet
```

The performance metrics are slightly different,
but the format of everything else is the same:

```{r}
results_multi$performance
```

## Continuous data

And here's an example for running continuous data, 
where the outcome column is numerical:

```{r, eval = FALSE}
results_cont <- run_ml(otu_mini_bin[, 2:11],
                       'glmnet',
                       outcome_colname = 'Otu00001',
                       seed = 2019)
```
```{r, echo = FALSE}
results_cont <- otu_mini_cont_results_glmnet
```

Again, the performance metrics are slightly different,
but the format of the rest is the same:

```{r}
results_cont$performance
```

# References
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

```{r setup, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
```{r render, eval = FALSE, echo = FALSE}
rmarkdown::render(here::here('vignettes','paper.Rmd'))
pkgdown::build_article('paper')
```

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
---
title: "Hyperparameter tuning"
author: "Begüm D. Topçuoğlu"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Hyperparameter tuning}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

One particularly important aspect of machine learning (ML) is hyperparameter tuning. 
A hyperparameter is a parameter that is set before the ML training begins. 
These parameters are tunable and they effect how well the model trains.
We must do a grid search for many hyperparameter possibilities and exhaust our search to pick the ideal value for the model and dataset. 
In this package, we do this during the cross-validation step.

Let's start with an example ML run. 
The input data to `run_ml()` is a dataframe where each row is a sample or observation.
One column (assumed to be the first) is the outcome of interest,
and all of the other columns are the features.
We package `otu_mini_bin` as a small example dataset with `mikropml`.

```{r}
#install.packages("devtools")
#devtools::install_github("SchlossLab/mikropml")
library(mikropml)
head(otu_mini_bin)
```

Before we train and evaluate a ML model, we can preprocess the data. 
You can learn more about this in the preprocessing vignette: `vignette("preprocess")`.

```{r}
preproc <- preprocess_data(dataset = otu_mini_bin,
                           outcome_colname = 'dx')
dat <- preproc$dat_transformed
```

We'll use `dat` for the following examples.

## The simplest way to `run_ml()`

As mentioned above, the minimal input is your dataset (`dataset`) and the machine learning model you want to use (`method`).

When we `run_ml()`, by default we do a 100 times repeated, 5-fold cross-validation, 
where we evaluate the hyperparameters in these 500 total iterations.

Say we want to run L2 regularized logistic regression. We do this with:

```{r, warning = FALSE}
results <- run_ml(dat,
                  'glmnet',
                  outcome_colname = 'dx',
                  cv_times = 100, 
                  seed = 2019)
```

You'll probably get a warning when you run this because the dataset is very small. If you want to learn more about that, check out the introductory vignette about training and evaluating a ML model: `vignette("introduction")`.

By default, `run_ml()` selects hyperparameters depending on the dataset and method used.

```{r}
results$trained_model
```

As you can see, the `alpha` hyperparameter is set to 0, which specifies L2 regularization.
`glmnet` gives us the option to run both L1 and L2 regularization.
If we change `alpha` to 1, we would run L1-regularized logistic regression. 
You can also tune `alpha` by specifying a variety of values between 0 and 1. 
When you use a value that is between 0 and 1, you are running elastic net.
The default hyperparameter `lambda` which adjusts the L2 regularization penalty is a range of values between 10^-4 to 10. 

When we look at the 100 repeated cross-validation performance metrics such as 
`AUC`, `Accuracy`, `prAUC` for each tested `lambda` value, 
we see that some are not appropriate for this dataset and some do better than others. 

```{r}
results$trained_model$results
```

## Customizing hyperparameters

In this example, we want to change the `lambda` values to provide a better range to test in the cross-validation step. 
We don't want to use the defaults but provide our own named list with new values. 

For example:

```{r}
new_hp <- list(alpha = 1, 
               lambda = c(0.00001, 0.0001, 0.001, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.1))
new_hp
```

Now let's run L2 logistic regression with the new `lambda` values:

```{r, warning = FALSE}
results <- run_ml(dat,
                  'glmnet',
                  outcome_colname = 'dx',
                  cv_times = 100,
                  hyperparameters = new_hp,
                  seed = 2019
                  )
results$trained_model
```

This time, we cover a larger and different range of `lambda` settings in cross-validation. 

How do we know which `lambda` value is the best one? 
To answer that, we need to run the ML pipeline on multiple data splits 
and look at the mean cross-validation performance of each `lambda` across those modeling experiments. 
We describe how to run the pipeline with multiple data splits in `vignette("parallel")`.

Here we train the model with the new `lambda` range we defined above. 
We run it 3 times each with a different seed, which will result in different
splits of the data into training and testing sets.
We can then use `plot_hp_performance` to see which `lambda` gives us the largest mean AUC value across modeling experiments. 


```{r, warning = FALSE}
results <- lapply(seq(100, 102), function(seed) {
   run_ml(dat, "glmnet", seed = seed, hyperparameters = new_hp)
 })
models <- lapply(results, function(x) x$trained_model)
hp_metrics <- combine_hp_performance(models)
plot_hp_performance(hp_metrics$dat, lambda, AUC)
```

As you can see, we get a mean maxima at `0.03` which is the best `lambda` value 
for this dataset when we run 3 data splits. 
The fact that we are seeing this maxima in the middle of our range and not at the edges, 
shows that we are providing a large enough range to exhaust our `lambda` search 
as we build the model.
We recommend the user to use this plot to make sure the best hyperparameter 
is not on the edges of the provided list. 
For a better understanding of the global maxima, 
it would be better to run more data splits by using more seeds. 
We picked 3 seeds to keep the runtime down for this vignette,
but for real-world data we recommend using many more seeds.

## Hyperparameter options

You can see which default hyperparameters would be used for your dataset with `get_hyperparams_list()`. 
Here are a few examples with built-in datasets we provide:

```{r}
get_hyperparams_list(otu_mini_bin, 'glmnet')
get_hyperparams_list(otu_mini_bin, 'rf')
get_hyperparams_list(otu_small, 'rf')
```

Here are the hyperparameters that are tuned for each of the modeling methods.
The output for all of them is very similar, so we won't go into those details.

### Regression

As mentioned above, `glmnet` uses the `alpha` parameter and `lambda` hyperparameter.
`alpha` of `0` is for L2 regularization (ridge).
`alpha` of `1` is for L1 regularization (lasso).
`alpha` in between is elastic net. You can also tune `alpha` like you would any other hyperparameter. 

Please refer to original `glmnet` documentation for more information:
https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

The default hyperparameters chosen by `run_ml()` are fixed for `glmnet`.

```{r, echo = FALSE}
mikropml:::set_hparams_glmnet()
```

### Random forest

When we run `rf` we are using the the `randomForest` package implementation.
We are tuning the `mtry` hyperparameter. 
This is the number of features that are randomly collected to be sampled at each tree node. 
This number needs to be less than the number of features in the dataset.
Please refer to the original documentation for more information:
https://cran.r-project.org/web/packages/randomForest/randomForest.pdf

By default, we take the square root of number of features in the dataset 
and we provide a range that is `[sqrt_features / 2, sqrt_features, sqrt_features * 2]`. 

For example if the number of features is 1000:

```{r, echo = FALSE}
mikropml:::set_hparams_rf(1000)
```

Similar to `glmnet` method, we can provide our own `mtry` range.

### Decision tree

When we run `rpart2`, we are running the `rpart` package implementation of decision tree.
We are tuning the `maxdepth` hyperparameter. 
This is the maximum depth of any node of the final tree.
Please refer to the original documentation for more information on maxdepth:
https://cran.r-project.org/web/packages/rpart/rpart.pdf

By default, we provide a range that is less than the number of features in the dataset. 

For example if we have 1000 features:

```{r, echo = FALSE}
mikropml:::set_hparams_rpart2(1000)
```

or 10 features:

```{r, echo = FALSE}
mikropml:::set_hparams_rpart2(10)
```


### SVM with radial basis kernel

When we run the `svmRadial` method, we are tuning the `C` and `sigma` hyperparameters. 
`sigma` defines how far the influence of a single training example reaches and `C` behaves as a regularization parameter.
Please refer to this great `sklearn` resource for more information on these hyperparameters:
https://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html

By default, we provide 2 separate range of values for the two hyperparameters. 

```{r, echo = FALSE}
mikropml:::set_hparams_svmRadial()
```

### XGBoost

When we run the `xgbTree` method, we are tuning the 
`nrounds`, `gamma`, `eta` `max_depth`, `colsample_bytree`, `min_child_weight` and `subsample` hyperparameters. 

You can read more about these hyperparameters here:
https://xgboost.readthedocs.io/en/latest/parameter.html

By default, we set the `nrounds`, `gamma`, `colsample_bytree` and `min_child_weight` 
to fixed values and we provide a range of values for `eta`, `max_depth` and `subsample`. 
All of these can be changed and optimized by the user by supplying a custom
named list of hyperparameters to `run_ml()`.

```{r, echo = FALSE}
mikropml:::set_hparams_xgbTree(1000)
```
---
title: "Preprocessing data"
author: "Zena Lapp"
output: rmarkdown::html_vignette
bibliography: paper.bib
vignette: >
  %\VignetteIndexEntry{Preprocessing data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Before training a model, it's often necessary and prudent to preprocess
your input data. We provide a function (`preprocess_data()`) to preprocess input
data. The defaults we chose are based on best practices used in
[FIDDLE](https://gitlab.eecs.umich.edu/mld3/FIDDLE/-/tree/master/)
[@tang_democratizing_2020]. Feel free to check out FIDDLE for more information
about data preprocessing!

`preprocess_data()` takes an input dataset where the rows are the samples and
the columns are the outcome variable and features. We preprocess the data as
follows:

- Remove missing outcome values.
- Convert any spaces in outcome names to underscores (`_`).
- Leave binary features as-is (except that categorical variables are converted to 0 and 1, and binary variables with missing features are split into two rows - see below for more details).
- Normalize continuous features using `caret::preProcess()` based on the method provided.
- Convert categorical features with more than 2 categories to 0 and 1 in multiple columns (one for each category, so each category has it's own column).
- Replace missing categorical data with 0.
- Impute missing continuous values with the median of the feature. 
- By default, remove all features with near-zero variance (option to also remove only features with zero variance).
- By default, collapse correlated features.

# It's running so slow!

Since I assume a lot of you won't read this entire vignette, I'm going to say
this at the beginning. If the `preprocess_data()` function is running super
slow, you should consider parallelizing it so it goes faster!
`preprocess_data()` also can report live progress updates. See
`vignette("parallel")` for details.


# Examples

We're going to start off simple and get more complicated, but if you want the
whole shebang at once, just scroll to the bottom.

First, we have to load `mikropml`:

```{r setup}
library(mikropml)
```

## Binary data

Let's start with only binary variables:

```{r}
# raw binary dataset
bin_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("no", "yes", "no"),
  var2 = c(0, 1, 1),
  var3 = factor(c("a","a","b"))
)
bin_df
```

In addition to the dataframe itself, you have to provide the name of the outcome column to `preprocess_data()`. Here's what the preprocessed data looks like:

```{r}
# preprocess raw binary data
preprocess_data(dataset = bin_df, outcome_colname = "outcome")
```

The output is a list: `dat_transformed` which has the transformed data, 
`grp_feats` which is a list of grouped features, and `removed_feats` which is a 
list of featuures that were removed. Here, `grp_feats` is `NULL` because there 
are no perfectly correlated features (e.g. `c(0,1,0)` and `c(0,1,0)`, or 
`c(0,1,0)` and `c(1,0,1)` - see below for more details). 

The first column (`var1`) in `dat_transformed` is a character and is changed to 
`var1_yes` that has zeros (no) and ones (yes). The values in the second column 
(`var2`) stay the same because it's already binary, but the name changes to 
`var2_1`. The third column (`var3`) is a factor and is also changed to binary 
where b is 1 and a is 0, as denoted by the new column name `var3_b`. 

## Categorical data

On to non-binary categorical data:

```{r}
# raw categorical dataset
cat_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c('a','b','c')
  )
cat_df
```

```{r}
# preprocess raw categorical data
preprocess_data(dataset = cat_df, outcome_colname = "outcome")
```

As you can see, this variable was split into 3 different columns - one for each
type (a, b, and c). And again, `grp_feats` is `NULL`.

## Continuous data

Now, looking at continuous variables:

```{r}
# raw continuous dataset
cont_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c(1,2,3)
  )
cont_df
```

```{r}
# preprocess raw continuous data
preprocess_data(dataset = cont_df, outcome_colname = "outcome")
```

Wow! Why did the numbers change? This is because the default is to normalize the 
data using `"center"` and `"scale"`. While this is often best practice, you may 
not want to normalize the data, or you may want to normalize the data in a 
different way. If you don't want to normalize the data, you can use 
`method=NULL`:

```{r, eval = FALSE}
# preprocess raw continuous data, no normalization
preprocess_data(dataset = cont_df, outcome_colname = "outcome", method = NULL)
```

You can also normalize the data in different ways. You can choose any method 
supported by the `method` argument of `caret::preProcess()` (see the 
`caret::preProcess()` docs for details). Note that these methods are only 
applied to continuous variables. 

Another feature of `preprocess_data()` is that if you provide continuous 
variables as characters, they will be converted to numeric:

```{r}
# raw continuous dataset as characters
cont_char_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("1","2","3")
  )
cont_char_df
```

```{r, eval = FALSE}
# preprocess raw continuous character data as numeric
preprocess_data(dataset = cont_char_df, outcome_colname = "outcome")
```

If you don't want this to happen, and you want character data to remain 
character data even if it can be converted to numeric, you can use 
`to_numeric=FALSE` and they will be kept as categorical:

```{r}
# preprocess raw continuous character data as characters
preprocess_data(dataset = cont_char_df, outcome_colname = "outcome", to_numeric = FALSE)
```

As you can see from this output, in this case the features are treated as groups 
rather than numbers (e.g. they are not normalized). 

## Collapse perfectly correlated features

By default, `preprocess_data()` collapses features that are perfectly positively 
or negatively correlated. This is because having multiple copies of those 
features does not add information to machine learning, and it makes `run_ml` 
faster.

```{r}
# raw correlated dataset
corr_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("no", "yes", "no"),
  var2 = c(0, 1, 0),
  var3 = c(1,0,1)
)
corr_df
```

```{r}
# preprocess raw correlated dataset
preprocess_data(dataset = corr_df, outcome_colname = "outcome")
```

As you can see, we end up with only one variable, as all 3 are grouped together.
Also, the second element in the list is no longer `NULL`. Instead, it tells you
that `grp1` contains `var1`, `var2`, and `var3`.

If you want to group positively correlated features, but not negatively
correlated features (e.g. for interpretability, or another downstream
application), you can do that by using `group_neg_corr=FALSE`:

```{r}
# preprocess raw correlated dataset; don't group negatively correlated features
preprocess_data(dataset = corr_df, outcome_colname = "outcome", group_neg_corr = FALSE)
```

Here, `var3` is kept on it's own because it's negatively correlated with `var1`
and `var2`. You can also choose to keep all features separate, even if they are
perfectly correlated, by using `collapse_corr_feats=FALSE`:

```{r}
# preprocess raw correlated dataset; don't group negatively correlated features
preprocess_data(dataset = corr_df, outcome_colname = "outcome", collapse_corr_feats = FALSE)
```

In this case, `grp_feats` will always be `NULL`.

## Data with near-zero variance

What if we have variables that are all zero, or all "no"? Those ones won't
contribute any information, so we remove them:

```{r}
# raw dataset with non-variable features
nonvar_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("no", "yes", "no"),
  var2 = c(0, 1, 1),
  var3 = c("no","no","no"),
  var4 = c(0,0,0),
  var5 = c(12,12,12)
)
nonvar_df
```

Here, `var3`, `var4`, and `var5` all have no variability, so these variables are
removed during preprocessing:

```{r}
# remove features with near-zero variance
preprocess_data(dataset = nonvar_df, outcome_colname = "outcome")
```

You can read the `caret::preProcess()` documentation for more information. By
default, we remove features with "near-zero variance" (`remove_var='nzv'`). This
uses the default arguments from `caret::nearZeroVar()`. However, particularly
with smaller datasets, you might not want to remove features with near-zero
variance. If you want to remove only features with zero variance, you can use
`remove_var='zv'`:

```{r}
# remove features with zero variance
preprocess_data(dataset = nonvar_df, outcome_colname = "outcome", remove_var = 'zv')
```

If you want to include all features, you can use the argument `remove_zv=NULL`.
For this to work, you cannot collapse correlated features (otherwise it errors
out because of the underlying `caret` function we use).

```{r}
# don't remove features with near-zero or zero variance
preprocess_data(dataset = nonvar_df, outcome_colname = "outcome", remove_var = NULL, collapse_corr_feats = FALSE)
```

If you want to be more nuanced in how you remove near-zero variance features
(e.g. change the default 10% cutoff for the percentage of distinct values out of
the total number of samples), you can use the `caret::preProcess()` function
after running `preprocess_data` with `remove_var=NULL` (see the
`caret::nearZeroVar()` function for more information).

## Missing data

`preprocess_data()` also deals with missing data. It:

- Removes missing outcome variables.
- Maintains zero variability in a feature if it already has no variability (i.e. the feature is removed if removing features with near-zero variance).
- Replaces missing binary and categorical variables with zero (after splitting into multiple columns).
- Replaces missing continuous data with the median value of that feature.

If you'd like to deal with missing data in a different way, please do that prior
to inputting the data to `preprocess_data()`.

### Remove missing outcome variables

```{r}
# raw dataset with missing outcome value
miss_oc_df <- data.frame(
  outcome = c("normal", "normal", "cancer",NA),
  var1 = c("no", "yes", "no","no"),
  var2 = c(0, 1, 1,1)
)
miss_oc_df
```

```{r}
# preprocess raw dataset with missing outcome value
preprocess_data(dataset = miss_oc_df, outcome_colname = "outcome")
```

### Maintain zero variability in a feature if it already has no variability

```{r}
# raw dataset with missing value in non-variable feature
miss_nonvar_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("no", "yes", "no"),
  var2 = c(NA, 1, 1)
)
miss_nonvar_df
```

```{r}
# preprocess raw dataset with missing value in non-variable feature
preprocess_data(dataset = miss_nonvar_df, outcome_colname = "outcome")
```

Here, the non-variable feature with missing data is removed because we removed
features with near-zero variance. If we maintained that feature, it'd be all
ones:

```{r}
# preprocess raw dataset with missing value in non-variable feature
preprocess_data(dataset = miss_nonvar_df, outcome_colname = "outcome", remove_var = NULL, collapse_corr_feats = FALSE)
```

### Replace missing binary and categorical variables with zero

```{r}
# raw dataset with missing value in categorical feature
miss_cat_df <- data.frame(
  outcome = c("normal", "normal", "cancer"),
  var1 = c("no", "yes", NA),
  var2 = c(NA, 1, 0)
)
miss_cat_df
```

```{r}
# preprocess raw dataset with missing value in non-variable feature
preprocess_data(dataset = miss_cat_df, outcome_colname = "outcome")
```

Here each binary variable is split into two, and the missing value is considered
zero for both of them.

### Replace missing continuous data with the median value of that feature

```{r}
# raw dataset with missing value in continuous feature
miss_cont_df <- data.frame(
  outcome = c("normal", "normal", "cancer","normal"),
  var1 = c(1,2,2,NA),
  var2 = c(1,2,3,NA)
)
miss_cont_df
```

Here we're not normalizing continuous features so it's easier to see what's
going on (i.e. the median value is used):

```{r}
# preprocess raw dataset with missing value in continuous feature
preprocess_data(dataset = miss_cont_df, outcome_colname = "outcome", method = NULL)
```

## Putting it all together

Here's some more complicated example raw data that puts everything we discussed together:

```{r}
test_df <- data.frame(
  outcome = c("normal", "normal", "cancer", NA),
  var1 = 1:4,
  var2 = c("a", "b", "c", "d"),
  var3 = c("no", "yes", "no", "no"),
  var4 = c(0, 1, 0, 0),
  var5 = c(0, 0, 0, 0),
  var6 = c("no", "no", "no", "no"),
  var7 = c(1, 1, 0, 0),
  var8 = c(5, 6, NA, 7),
  var9 = c(NA, "x", "y", "z"),
  var10 = c(1, 0, NA, NA),
  var11 = c(1, 1, NA, NA),
  var12 = c("1", "2", "3", "4")
)
test_df
```

Let's throw this into the preprocessing function with the default values:

```{r}
preprocess_data(dataset = test_df, outcome_colname = "outcome")
```

As you can see, we got several messages:

- One of the samples (row 4) was removed because the outcome value was missing.
- One of the variables in a feature with no variation had a missing value that was replaced with the the non-varying value (`var11`).
- Four categorical missing values were replaced with zero (`var9`). 
There are 4 missing rather than just 1 (like in the raw data) because we split the categorical variable into 4 different columns first.
- One missing continuous value was imputed using the median value of that feature (`var8`).

Additionally, you can see that the continuous variables were normalized, the
categorical variables were all changed to binary, and several features were
grouped together. The variables in each group can be found in `grp_feats`.

## Next step: train and evaluate your model!

After you preprocess your data (either using `preprocess_data()` or by
preprocessing the data on your own),
you're ready to train and evaluate machine learning models! 
Please see `run_ml()` information about training models. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{get_outcome_type}
\alias{get_outcome_type}
\title{Get outcome type.}
\usage{
get_outcome_type(outcomes_vec)
}
\arguments{
\item{outcomes_vec}{Vector of outcomes.}
}
\value{
Outcome type (continuous, binary, or multiclass).
}
\description{
If the outcome is numeric, the type is continuous.
Otherwise, the outcome type is binary if there are only two outcomes or
multiclass if there are more than two outcomes.
}
\examples{
get_outcome_type(c(1, 2, 1))
get_outcome_type(c("a", "b", "b"))
get_outcome_type(c("a", "b", "c"))
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_hp_performance}
\alias{plot_hp_performance}
\title{Plot hyperparameter performance metrics}
\usage{
plot_hp_performance(dat, param_col, metric_col)
}
\arguments{
\item{dat}{dataframe of hyperparameters and performance metric (e.g. from \code{get_hp_performance()} or \code{combine_hp_performance()})}

\item{param_col}{hyperparameter to be plotted. must be a column in \code{dat}.}

\item{metric_col}{performance metric. must be a column in \code{dat}.}
}
\value{
ggplot of hyperparameter performance.
}
\description{
Plot hyperparameter performance metrics
}
\examples{
# plot for a single `run_ml()` call
hp_metrics <- get_hp_performance(otu_mini_bin_results_glmnet$trained_model)
hp_metrics
plot_hp_performance(hp_metrics$dat, lambda, AUC)
\dontrun{
# plot for multiple `run_ml()` calls
results <- lapply(seq(100, 102), function(seed) {
  run_ml(otu_small, "glmnet", seed = seed)
})
models <- lapply(results, function(x) x$trained_model)
hp_metrics <- combine_hp_performance(models)
plot_hp_performance(hp_metrics$dat, lambda, AUC)
}
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}

Kelly Sovacool \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mikropml.R
\docType{package}
\name{mikropml}
\alias{mikropml}
\title{mikropml: User-Friendly R Package for Robust Machine Learning Pipelines}
\description{
\code{mikropml} implements supervised machine learning pipelines using regression,
support vector machines, decision trees, random forest, or gradient-boosted trees.
The main functions are \code{preprocess_data()} to process your data prior to
running machine learning, and \code{run_ml()} to run machine learning.
}
\section{Authors}{

\itemize{
\item Begüm D. Topçuoğlu (\href{https://orcid.org/0000-0003-3140-537X}{ORCID})
\item Zena Lapp (\href{https://orcid.org/0000-0003-4674-2176}{ORCID})
\item Kelly L. Sovacool (\href{https://orcid.org/0000-0003-3283-829X}{ORCID})
\item Evan Snitkin (\href{https://orcid.org/0000-0001-8409-278X}{ORCID})
\item Jenna Wiens (\href{https://orcid.org/0000-0002-1057-7722}{ORCID})
\item Patrick D. Schloss (\href{https://orcid.org/0000-0002-6935-4275}{ORCID})
}
}

\section{See vignettes}{

\itemize{
\item \href{http://www.schlosslab.org/mikropml/articles/introduction.html}{Introduction}
\item \href{http://www.schlosslab.org/mikropml/articles/preprocess.html}{Preprocessing data}
\item \href{http://www.schlosslab.org/mikropml/articles/tuning.html}{Hyperparameter tuning}
\item \href{http://www.schlosslab.org/mikropml/articles/parallel.html}{Parallel processing}
\item \href{http://www.schlosslab.org/mikropml/articles/paper.html}{The mikropml paper}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin_results_svmRadial}
\alias{otu_mini_bin_results_svmRadial}
\title{Results from running the pipline with svmRadial on \code{otu_mini_bin}}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_bin_results_svmRadial
}
\description{
Results from running the pipline with svmRadial on \code{otu_mini_bin}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{replace_spaces}
\alias{replace_spaces}
\title{Replace spaces in all elements of a character vector with underscores}
\usage{
replace_spaces(x, new_char = "_")
}
\arguments{
\item{x}{a character vector}

\item{new_char}{the character to replace spaces (default: \verb{_})}
}
\value{
character vector with all spaces replaced with \code{new_char}
}
\description{
Replace spaces in all elements of a character vector with underscores
}
\examples{
dat <- data.frame(
  dx = c("outcome 1", "outcome 2", "outcome 1"),
  a = 1:3, b = c(5, 7, 1)
)
dat$dx <- replace_spaces(dat$dx)
dat
}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_model_performance}
\alias{plot_model_performance}
\title{Plot performance metrics for multiple ML runs with different parameters}
\usage{
plot_model_performance(performance_df)
}
\arguments{
\item{performance_df}{dataframe of performance results from multiple calls to \code{run_ml()}}
}
\value{
A ggplot2 plot of performance.
}
\description{
ggplot2 is required to use this function.
}
\examples{
\dontrun{
# call `run_ml()` multiple times with different seeds
results_lst <- lapply(seq(100, 104), function(seed) {
  run_ml(otu_small, "glmnet", seed = seed)
})
# extract and combine the performance results
perf_df <- lapply(results_lst, function(result) {
  result[["performance"]]
}) \%>\%
  dplyr::bind_rows()
# plot the performance results
p <- plot_model_performance(perf_df)


# call `run_ml()` with different ML methods
param_grid <- expand.grid(
  seeds = seq(100, 104),
  methods = c("glmnet", "rf")
)
results_mtx <- mapply(
  function(seed, method) {
    run_ml(otu_mini_bin, method, seed = seed, kfold = 2)
  },
  param_grid$seeds, param_grid$methods
)
# extract and combine the performance results
perf_df2 <- dplyr::bind_rows(results_mtx["performance", ])
# plot the performance results
p <- plot_model_performance(perf_df2)

# you can continue adding layers to customize the plot
p +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  coord_flip()
}
}
\author{
Begüm Topçuoglu, \email{topcuoglu.begum@gmail.com}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_ml.R
\name{run_ml}
\alias{run_ml}
\title{Run the machine learning pipeline}
\usage{
run_ml(
  dataset,
  method,
  outcome_colname = NULL,
  hyperparameters = NULL,
  find_feature_importance = FALSE,
  calculate_performance = TRUE,
  kfold = 5,
  cv_times = 100,
  cross_val = NULL,
  training_frac = 0.8,
  perf_metric_function = NULL,
  perf_metric_name = NULL,
  groups = NULL,
  group_partitions = NULL,
  corr_thresh = 1,
  ntree = 1000,
  seed = NA
)
}
\arguments{
\item{dataset}{Dataframe with an outcome variable and other columns as features.}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{hyperparameters}{Dataframe of hyperparameters
(default \code{NULL}; sensible defaults will be chosen automatically).}

\item{find_feature_importance}{Run permutation importance (default: \code{FALSE}).
\code{TRUE} is recommended if you would like to identify features important for
predicting your outcome, but it is resource-intensive.}

\item{calculate_performance}{Whether to calculate performance metrics (default: \code{TRUE}).
You might choose to skip this if you do not perform cross-validation during model training.}

\item{kfold}{Fold number for k-fold cross-validation (default: \code{5}).}

\item{cv_times}{Number of cross-validation partitions to create (default: \code{100}).}

\item{cross_val}{a custom cross-validation scheme from \code{caret::trainControl()}
(default: \code{NULL}, uses \code{kfold} cross validation repeated \code{cv_times}).
\code{kfold} and \code{cv_times} are ignored if the user provides a custom cross-validation scheme.
See the \code{caret::trainControl()} docs for information on how to use it.}

\item{training_frac}{Fraction of data for training set (default: \code{0.8}). Rows
from the dataset will be randomly selected for the training set, and all
remaining rows will be used in the testing set. Alternatively, if you
provide a vector of integers, these will be used as the row indices for the
training set. All remaining rows will be used in the testing set.}

\item{perf_metric_function}{Function to calculate the performance metric to
be used for cross-validation and test performance. Some functions are
provided by caret (see \code{\link[caret:postResample]{caret::defaultSummary()}}).
Defaults: binary classification = \code{twoClassSummary},
multi-class classification = \code{multiClassSummary},
regression = \code{defaultSummary}.}

\item{perf_metric_name}{The column name from the output of the function
provided to perf_metric_function that is to be used as the performance metric.
Defaults: binary classification = \code{"ROC"},
multi-class classification = \code{"logLoss"},
regression = \code{"RMSE"}.}

\item{groups}{Vector of groups to keep together when splitting the data into
train and test sets. If the number of groups in the training set is larger
than \code{kfold}, the groups will also be kept together for cross-validation.
Length matches the number of rows in the dataset (default: \code{NULL}).}

\item{group_partitions}{Specify how to assign \code{groups} to the training and
testing partitions (default: \code{NULL}). If \code{groups} specifies that some
samples belong to group \code{"A"} and some belong to group \code{"B"}, then setting
\code{group_partitions = list(train = c("A", "B"), test = c("B"))} will result
in all samples from group \code{"A"} being placed in the training set, some
samples from \code{"B"} also in the training set, and the remaining samples from
\code{"B"} in the testing set. The partition sizes will be as close to
\code{training_frac} as possible. If the number of groups in the training set is
larger than \code{kfold}, the groups will also be kept together for
cross-validation.}

\item{corr_thresh}{For feature importance, group correlations
above or equal to \code{corr_thresh} (range \code{0} to \code{1}; default: \code{1}).}

\item{ntree}{For random forest, how many trees to use (default: \code{1000}).
Note that caret doesn't allow this parameter to be tuned.}

\item{seed}{Random seed (default: \code{NA}).
Your results will only be reproducible if you set a seed.}
}
\value{
Named list with results:
\itemize{
\item \code{trained_model}: Output of \code{\link[caret:train]{caret::train()}}, including the best model.
\item \code{test_data}: Part of the data that was used for testing.
\item \code{performance}: Dataframe of performance metrics. The first column is the cross-validation performance metric, and the last two columns are the ML method used and the seed (if one was set), respectively. All other columns are performance metrics calculated on the test data. This contains only one row, so you can easily combine performance dataframes from multiple calls to \code{run_ml()} (see \code{vignette("parallel")}).
\item \code{feature_importance}: If feature importances were calculated, a dataframe where each row is a feature or correlated group. The columns are the performance metric of the permuted data, the difference between the true performance metric and the performance metric of the permuted data (true - permuted), the feature name, the ML method, the performance metric name, and the seed (if provided). For AUC and RMSE, the higher perf_metric_diff is, the more important that feature is for predicting the outcome. For log loss, the lower perf_metric_diff is, the more important that feature is for predicting the outcome.
}
}
\description{
This function runs machine learning (ML), evaluates the best model,
and optionally calculates feature importance using the framework
outlined in Topçuoğlu \emph{et al.} 2020 (\doi{10.1128/mBio.00434-20}).
Required inputs are a dataframe with an outcome variable and other columns
as features, as well as the ML method.
See \code{vignette('introduction')} for more details.
}
\section{More details}{


For more details, please see \href{http://www.schlosslab.org/mikropml/articles/}{the vignettes}.
}

\examples{
\dontrun{

# regression
run_ml(otu_small, "glmnet",
  seed = 2019
)

# random forest w/ feature importance
run_ml(otu_small, "rf",
  outcome_colname = "dx",
  find_feature_importance = TRUE
)

# custom cross validation & hyperparameters
run_ml(otu_mini_bin[, 2:11],
  "glmnet",
  outcome_colname = "Otu00001",
  seed = 2019,
  hyperparameters = list(lambda = c(1e-04), alpha = 0),
  cross_val = caret::trainControl(method = "none"),
  calculate_performance = FALSE
)
}
}
\author{
Begüm Topçuoğlu, \email{topcuoglu.begum@gmail.com}

Zena Lapp, \email{zenalapp@umich.edu}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_val.R
\name{define_cv}
\alias{define_cv}
\title{Define cross-validation scheme and training parameters}
\usage{
define_cv(
  train_data,
  outcome_colname,
  hyperparams_list,
  perf_metric_function,
  class_probs,
  kfold = 5,
  cv_times = 100,
  groups = NULL,
  group_partitions = NULL
)
}
\arguments{
\item{train_data}{Dataframe for training model.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{hyperparams_list}{Named list of lists of hyperparameters.}

\item{perf_metric_function}{Function to calculate the performance metric to
be used for cross-validation and test performance. Some functions are
provided by caret (see \code{\link[caret:postResample]{caret::defaultSummary()}}).
Defaults: binary classification = \code{twoClassSummary},
multi-class classification = \code{multiClassSummary},
regression = \code{defaultSummary}.}

\item{class_probs}{Whether to use class probabilities (TRUE for categorical outcomes, FALSE for numeric outcomes).}

\item{kfold}{Fold number for k-fold cross-validation (default: \code{5}).}

\item{cv_times}{Number of cross-validation partitions to create (default: \code{100}).}

\item{groups}{Vector of groups to keep together when splitting the data into
train and test sets. If the number of groups in the training set is larger
than \code{kfold}, the groups will also be kept together for cross-validation.
Length matches the number of rows in the dataset (default: \code{NULL}).}

\item{group_partitions}{Specify how to assign \code{groups} to the training and
testing partitions (default: \code{NULL}). If \code{groups} specifies that some
samples belong to group \code{"A"} and some belong to group \code{"B"}, then setting
\code{group_partitions = list(train = c("A", "B"), test = c("B"))} will result
in all samples from group \code{"A"} being placed in the training set, some
samples from \code{"B"} also in the training set, and the remaining samples from
\code{"B"} in the testing set. The partition sizes will be as close to
\code{training_frac} as possible. If the number of groups in the training set is
larger than \code{kfold}, the groups will also be kept together for
cross-validation.}
}
\value{
Caret object for trainControl that controls cross-validation
}
\description{
Define cross-validation scheme and training parameters
}
\examples{
training_inds <- get_partition_indices(otu_small \%>\% dplyr::pull("dx"),
  training_frac = 0.8,
  groups = NULL
)
train_data <- otu_small[training_inds, ]
test_data <- otu_small[-training_inds, ]
cv <- define_cv(train_data,
  outcome_colname = "dx",
  hyperparams_list = get_hyperparams_list(otu_small, "glmnet"),
  perf_metric_function = caret::multiClassSummary,
  class_probs = TRUE,
  kfold = 5
)
}
\author{
Begüm Topçuoğlu, \email{topcuoglu.begum@gmail.com}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{get_hp_performance}
\alias{get_hp_performance}
\title{Get hyperparameter performance metrics}
\usage{
get_hp_performance(trained_model)
}
\arguments{
\item{trained_model}{trained model (e.g. from \code{run_ml()})}
}
\value{
Named list:
\itemize{
\item \code{dat}: Dataframe of performance metric for each group of hyperparameters.
\item \code{params}: Hyperparameters tuned.
\item \code{metric}: Performance metric used.
}
}
\description{
Get hyperparameter performance metrics
}
\examples{
get_hp_performance(otu_mini_bin_results_glmnet$trained_model)
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}

Kelly Sovacool \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_small}
\alias{otu_small}
\title{Small OTU abundance dataset}
\format{
A data frame with 60 rows and 61 variables.
The \code{dx} column is the diagnosis: healthy or cancerous (colorectal).
All other columns are OTU relative abundances.
}
\usage{
otu_small
}
\description{
A dataset containing relatives abundances of 60 OTUs for 60 human stool samples.
This is a subset of the data provided in \code{extdata/otu_large.csv}, which was
used in \href{https://journals.asm.org/doi/10.1128/mbio.00434-20}{Topçuoğlu \emph{et al.} 2020}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_cont_results_glmnet}
\alias{otu_mini_cont_results_glmnet}
\title{Results from running the pipeline with glmnet on \code{otu_mini_bin} with \code{Otu00001}
as the outcome}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_cont_results_glmnet
}
\description{
Results from running the pipeline with glmnet on \code{otu_mini_bin} with \code{Otu00001}
as the outcome
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_cv}
\alias{otu_mini_cv}
\title{Cross validation on \code{train_data_mini} with grouped features.}
\format{
An object of class \code{list} of length 27.
}
\usage{
otu_mini_cv
}
\description{
Cross validation on \code{train_data_mini} with grouped features.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_multi_results_glmnet}
\alias{otu_mini_multi_results_glmnet}
\title{Results from running the pipeline with glmnet on \code{otu_mini_multi} for
multiclass outcomes}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_multi_results_glmnet
}
\description{
Results from running the pipeline with glmnet on \code{otu_mini_multi} for
multiclass outcomes
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin_results_xgbTree}
\alias{otu_mini_bin_results_xgbTree}
\title{Results from running the pipline with xbgTree on \code{otu_mini_bin}}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_bin_results_xgbTree
}
\description{
Results from running the pipline with xbgTree on \code{otu_mini_bin}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{calc_perf_metrics}
\alias{calc_perf_metrics}
\title{Get performance metrics for test data}
\usage{
calc_perf_metrics(
  test_data,
  trained_model,
  outcome_colname,
  perf_metric_function,
  class_probs
)
}
\arguments{
\item{test_data}{Held out test data: dataframe of outcome and features.}

\item{trained_model}{Trained model from \code{\link[caret:train]{caret::train()}}.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{perf_metric_function}{Function to calculate the performance metric to
be used for cross-validation and test performance. Some functions are
provided by caret (see \code{\link[caret:postResample]{caret::defaultSummary()}}).
Defaults: binary classification = \code{twoClassSummary},
multi-class classification = \code{multiClassSummary},
regression = \code{defaultSummary}.}

\item{class_probs}{Whether to use class probabilities (TRUE for categorical outcomes, FALSE for numeric outcomes).}
}
\value{
Dataframe of performance metrics.
}
\description{
Get performance metrics for test data
}
\examples{
\dontrun{
results <- run_ml(otu_small, "glmnet", kfold = 2, cv_times = 2)
calc_perf_metrics(results$test_data,
  results$trained_model,
  "dx",
  multiClassSummary,
  class_probs = TRUE
)
}
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperparameters.R
\name{get_tuning_grid}
\alias{get_tuning_grid}
\title{Generate the tuning grid for tuning hyperparameters}
\usage{
get_tuning_grid(hyperparams_list, method)
}
\arguments{
\item{hyperparams_list}{Named list of lists of hyperparameters.}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}
}
\value{
The tuning grid.
}
\description{
Generate the tuning grid for tuning hyperparameters
}
\examples{
ml_method <- "glmnet"
hparams_list <- get_hyperparams_list(otu_small, ml_method)
get_tuning_grid(hparams_list, ml_method)
}
\author{
Begüm Topçuoğlu, \email{topcuoglu.begum@gmail.com}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_multi}
\alias{otu_mini_multi}
\title{Mini OTU abundance dataset with 3 categorical variables}
\format{
A data frame
The \code{dx} column is the colorectal cancer diagnosis: adenoma, carcinoma, normal.
All other columns are OTU relative abundances.
}
\usage{
otu_mini_multi
}
\description{
A dataset containing relatives abundances of OTUs for human stool samples
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin}
\alias{otu_mini_bin}
\title{Mini OTU abundance dataset}
\format{
A data frame
The \code{dx} column is the diagnosis: healthy or cancerous (colorectal).
All other columns are OTU relative abundances.
}
\usage{
otu_mini_bin
}
\description{
A dataset containing relatives abundances of OTUs for human stool samples
with a binary outcome, \code{dx}.
This is a subset of \code{otu_small}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{tidy_perf_data}
\alias{tidy_perf_data}
\title{Tidy the performance dataframe}
\usage{
tidy_perf_data(performance_df)
}
\arguments{
\item{performance_df}{dataframe of performance results from multiple calls to \code{run_ml()}}
}
\value{
Tidy dataframe with model performance metrics.
}
\description{
Used by \code{plot_model_performance()}.
}
\examples{
\dontrun{
# call `run_ml()` multiple times with different seeds
results_lst <- lapply(seq(100, 104), function(seed) {
  run_ml(otu_small, "glmnet", seed = seed)
})
# extract and combine the performance results
perf_df <- lapply(results_lst, function(result) {
  result[["performance"]]
}) \%>\%
  dplyr::bind_rows()
# make it pretty!
tidy_perf_data(perf_df)
}
}
\author{
Begüm Topçuoglu, \email{topcuoglu.begum@gmail.com}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin_results_rpart2}
\alias{otu_mini_bin_results_rpart2}
\title{Results from running the pipline with rpart2 on \code{otu_mini_bin}}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_bin_results_rpart2
}
\description{
Results from running the pipline with rpart2 on \code{otu_mini_bin}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperparameters.R
\name{get_hyperparams_list}
\alias{get_hyperparams_list}
\title{Set hyperparameters based on ML method and dataset characteristics}
\usage{
get_hyperparams_list(dataset, method)
}
\arguments{
\item{dataset}{Dataframe with an outcome variable and other columns as features.}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}
}
\value{
Named list of hyperparameters.
}
\description{
For more details see the vignette on \href{http://www.schlosslab.org/mikropml/articles/tuning.html}{hyperparameter tuning}.
}
\examples{
get_hyperparams_list(otu_mini_bin, "rf")
get_hyperparams_list(otu_small, "rf")
get_hyperparams_list(otu_mini_bin, "rpart2")
get_hyperparams_list(otu_small, "rpart2")
}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin_results_glmnet}
\alias{otu_mini_bin_results_glmnet}
\title{Results from running the pipline with L2 logistic regression on \code{otu_mini_bin} with feature importance and grouping}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_bin_results_glmnet
}
\description{
Results from running the pipline with L2 logistic regression on \code{otu_mini_bin} with feature importance and grouping
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_feats.R
\name{group_correlated_features}
\alias{group_correlated_features}
\title{Group correlated features}
\usage{
group_correlated_features(
  features,
  corr_thresh = 1,
  group_neg_corr = TRUE,
  corr_method = "spearman"
)
}
\arguments{
\item{features}{a dataframe with each column as a feature for ML}

\item{corr_thresh}{For feature importance, group correlations
above or equal to \code{corr_thresh} (range \code{0} to \code{1}; default: \code{1}).}

\item{group_neg_corr}{Whether to group negatively correlated features
together (e.g. c(0,1) and c(1,0)).}

\item{corr_method}{correlation method. options or the same as those supported
by \code{stats::cor}: spearman, pearson, kendall. (default: spearman)}
}
\value{
vector where each element is a group of correlated features
separated by pipes (\code{|})
}
\description{
Group correlated features
}
\examples{
features <- data.frame(
  a = 1:3, b = 2:4, c = c(1, 0, 1),
  d = (5:7), e = c(5, 1, 4), f = c(-1, 0, -1)
)
group_correlated_features(features)
}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_model.R
\name{train_model}
\alias{train_model}
\title{Train model using \code{\link[caret:train]{caret::train()}}.}
\usage{
train_model(
  model_formula,
  train_data,
  method,
  cv,
  perf_metric_name,
  tune_grid,
  ntree
)
}
\arguments{
\item{model_formula}{Model formula, typically created with \code{stats::as.formula()}.}

\item{train_data}{Training data. Expected to be a subset of the full dataset.}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}

\item{cv}{Cross-validation caret scheme from \code{define_cv()}.}

\item{perf_metric_name}{The column name from the output of the function
provided to perf_metric_function that is to be used as the performance metric.
Defaults: binary classification = \code{"ROC"},
multi-class classification = \code{"logLoss"},
regression = \code{"RMSE"}.}

\item{tune_grid}{Tuning grid from \code{get_tuning_grid()}.}

\item{ntree}{For random forest, how many trees to use (default: \code{1000}).
Note that caret doesn't allow this parameter to be tuned.}
}
\value{
Trained model from \code{\link[caret:train]{caret::train()}}.
}
\description{
Train model using \code{\link[caret:train]{caret::train()}}.
}
\examples{
\dontrun{
training_data <- otu_mini_bin_results_glmnet$trained_model$trainingData \%>\%
  dplyr::rename(dx = .outcome)
method <- "rf"
hyperparameters <- get_hyperparams_list(otu_mini_bin, method)
cross_val <- define_cv(training_data,
  "dx",
  hyperparameters,
  perf_metric_function = caret::multiClassSummary,
  class_probs = TRUE,
  cv_times = 2
)
tune_grid <- get_tuning_grid(hyperparameters, method)

rf_model <- train_model(
  stats::as.formula(paste("dx", "~ .")),
  training_data,
  method,
  cross_val,
  "AUC",
  tune_grid,
  1000
)
rf_model$results \%>\% dplyr::select(mtry, AUC, prAUC)
}
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition.R
\name{get_partition_indices}
\alias{get_partition_indices}
\title{Select indices to partition the data into training & testing sets.}
\usage{
get_partition_indices(
  outcomes,
  training_frac = 0.8,
  groups = NULL,
  group_partitions = NULL
)
}
\arguments{
\item{outcomes}{vector of outcomes}

\item{training_frac}{Fraction of data for training set (default: \code{0.8}). Rows
from the dataset will be randomly selected for the training set, and all
remaining rows will be used in the testing set. Alternatively, if you
provide a vector of integers, these will be used as the row indices for the
training set. All remaining rows will be used in the testing set.}

\item{groups}{Vector of groups to keep together when splitting the data into
train and test sets. If the number of groups in the training set is larger
than \code{kfold}, the groups will also be kept together for cross-validation.
Length matches the number of rows in the dataset (default: \code{NULL}).}

\item{group_partitions}{Specify how to assign \code{groups} to the training and
testing partitions (default: \code{NULL}). If \code{groups} specifies that some
samples belong to group \code{"A"} and some belong to group \code{"B"}, then setting
\code{group_partitions = list(train = c("A", "B"), test = c("B"))} will result
in all samples from group \code{"A"} being placed in the training set, some
samples from \code{"B"} also in the training set, and the remaining samples from
\code{"B"} in the testing set. The partition sizes will be as close to
\code{training_frac} as possible. If the number of groups in the training set is
larger than \code{kfold}, the groups will also be kept together for
cross-validation.}
}
\value{
Vector of row indices for the training set.
}
\description{
Use this function to get the row indices for the training set.
}
\details{
If \code{groups} is \code{NULL}, uses \link[caret]{createDataPartition}.
Otherwise, uses \code{create_grouped_data_partition()}.

Set the seed prior to calling this function if you would like your data
partitions to be reproducible (recommended).
}
\examples{
training_inds <- get_partition_indices(otu_mini_bin$dx)
train_data <- otu_mini_bin[training_inds, ]
test_data <- otu_mini_bin[-training_inds, ]
}
\author{
Kelly Sovacool, {sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{.data}
\alias{contr.ltfr}
\alias{!!}
\alias{:=}
\title{dplyr pipe}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{caret}{\code{\link[caret:dummyVars]{contr.ltfr}}}

  \item{dplyr}{\code{\link[dplyr:reexports]{\%>\%}}}

  \item{rlang}{\code{\link[rlang:nse-force]{!!}}, \code{\link[rlang:tidyeval-data]{.data}}, \code{\link[rlang:nse-force]{:=}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_data.R
\name{get_caret_processed_df}
\alias{get_caret_processed_df}
\title{Get preprocessed dataframe for continuous variables}
\usage{
get_caret_processed_df(features, method)
}
\arguments{
\item{features}{Dataframe of features for machine learning}

\item{method}{Methods to preprocess the data, described in
\code{\link[caret:preProcess]{caret::preProcess()}} (default: \code{c("center","scale")}, use \code{NULL} for
no normalization).}
}
\value{
Named list:
\itemize{
\item \code{processed}: Dataframe of processed features.
\item \code{removed}: Names of any features removed during preprocessing.
}
}
\description{
Get preprocessed dataframe for continuous variables
}
\examples{
get_caret_processed_df(mikropml::otu_small[, 2:ncol(otu_small)], c("center", "scale"))
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_bin_results_rf}
\alias{otu_mini_bin_results_rf}
\title{Results from running the pipline with random forest on \code{otu_mini_bin}}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_bin_results_rf
}
\description{
Results from running the pipline with random forest on \code{otu_mini_bin}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_data.R
\name{remove_singleton_columns}
\alias{remove_singleton_columns}
\title{Remove columns appearing in only \code{threshold} row(s) or fewer.}
\usage{
remove_singleton_columns(dat, threshold = 1)
}
\arguments{
\item{dat}{dataframe}

\item{threshold}{Number of rows. If a column only has non-zero & non-NA values
in \code{threshold} row(s) or fewer, it will be removed.}
}
\value{
dataframe without singleton columns
}
\description{
Removes columns which only have non-zero & non-NA values in \code{threshold} row(s) or fewer.
}
\examples{
remove_singleton_columns(data.frame(a = 1:3, b = c(0, 1, 0), c = 4:6))
remove_singleton_columns(data.frame(a = 1:3, b = c(0, 1, 0), c = 4:6), threshold = 0)
remove_singleton_columns(data.frame(a = 1:3, b = c(0, 1, NA), c = 4:6))
remove_singleton_columns(data.frame(a = 1:3, b = c(1, 1, 1), c = 4:6))
}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}

Courtney Armour
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{randomize_feature_order}
\alias{randomize_feature_order}
\title{Randomize feature order to eliminate any position-dependent effects}
\usage{
randomize_feature_order(dataset, outcome_colname)
}
\arguments{
\item{dataset}{Dataframe with an outcome variable and other columns as features.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}
}
\value{
Dataset with feature order randomized.
}
\description{
Randomize feature order to eliminate any position-dependent effects
}
\examples{
dat <- data.frame(
  outcome = c("1", "2", "3"),
  a = 4:6, b = 7:9, c = 10:12, d = 13:15
)
randomize_feature_order(dat, "outcome")
}
\author{
Nick Lesniak, \email{nlesniak@umich.edu}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_multi_group}
\alias{otu_mini_multi_group}
\title{Groups for otu_mini_multi}
\format{
An object of class \code{character} of length 490.
}
\usage{
otu_mini_multi_group
}
\description{
Groups for otu_mini_multi
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_importance.R
\name{get_feature_importance}
\alias{get_feature_importance}
\title{Get feature importance using the permutation method}
\usage{
get_feature_importance(
  trained_model,
  train_data,
  test_data,
  outcome_colname,
  perf_metric_function,
  perf_metric_name,
  class_probs,
  method,
  seed = NA,
  corr_thresh = 1,
  groups = NULL,
  nperms = 100,
  corr_method = "spearman"
)
}
\arguments{
\item{trained_model}{Trained model from \code{\link[caret:train]{caret::train()}}.}

\item{train_data}{Training data: dataframe of outcome and features.}

\item{test_data}{Held out test data: dataframe of outcome and features.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{perf_metric_function}{Function to calculate the performance metric to
be used for cross-validation and test performance. Some functions are
provided by caret (see \code{\link[caret:postResample]{caret::defaultSummary()}}).
Defaults: binary classification = \code{twoClassSummary},
multi-class classification = \code{multiClassSummary},
regression = \code{defaultSummary}.}

\item{perf_metric_name}{The column name from the output of the function
provided to perf_metric_function that is to be used as the performance metric.
Defaults: binary classification = \code{"ROC"},
multi-class classification = \code{"logLoss"},
regression = \code{"RMSE"}.}

\item{class_probs}{Whether to use class probabilities (TRUE for categorical outcomes, FALSE for numeric outcomes).}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}

\item{seed}{Random seed (default: \code{NA}).
Your results will only be reproducible if you set a seed.}

\item{corr_thresh}{For feature importance, group correlations
above or equal to \code{corr_thresh} (range \code{0} to \code{1}; default: \code{1}).}

\item{groups}{Vector of feature names to group together during permutation.
Each element should be a string with feature names separated by a pipe
character (\code{|}). If this is \code{NULL} (default), correlated features will be
grouped together based on \code{corr_thresh}.}

\item{nperms}{number of permutations to perform (default: \code{100}).}

\item{corr_method}{correlation method. options or the same as those supported
by \code{stats::cor}: spearman, pearson, kendall. (default: spearman)}
}
\value{
Data frame with performance metrics for when each feature (or group
of correlated features; \code{names}) is permuted (\code{perf_metric}), differences
between the actual test performance metric on and the permuted performance
metric (\code{perf_metric_diff}; test minus permuted performance), and the
p-value (\code{pvalue}: the probability of obtaining the actual performance
value under the null hypothesis). Features with a larger \code{perf_metric_diff}
are more important. The performance metric name (\code{perf_metric_name}) and
seed (\code{seed}) are also returned.
}
\description{
Calculates feature importance using a trained model and test data. Requires
the \code{future.apply} package.
}
\details{
For permutation tests, the p-value is the number of permutation statistics
that are greater than the test statistic, divided by the number of
permutations. In our case, the permutation statistic is the model performance
(e.g. AUROC) after randomizing the order of observations for one feature, and
the test statistic is the actual performance on the test data. By default we
perform 100 permutations per feature; increasing this will increase the
precision of estimating the null distribution, but also increases runtime.
The p-value represents the probability of obtaining the actual performance in
the event that the null hypothesis is true, where the null hypothesis is that
the feature is not important for model performance.
}
\examples{
\dontrun{
results <- run_ml(otu_small, "glmnet", kfold = 2, cv_times = 2)
names(results$trained_model$trainingData)[1] <- "dx"
get_feature_importance(results$trained_model,
  results$trained_model$trainingData, results$test_data,
  "dx",
  multiClassSummary, "AUC",
  class_probs = TRUE, method = "glmnet"
)

# optionally, you can group features together with a custom grouping
get_feature_importance(results$trained_model,
  results$trained_model$trainingData, results$test_data,
  "dx",
  multiClassSummary, "AUC",
  class_probs = TRUE, method = "glmnet",
  groups = c(
    "Otu00007", "Otu00008", "Otu00009", "Otu00011", "Otu00012",
    "Otu00015", "Otu00016", "Otu00018", "Otu00019", "Otu00020", "Otu00022",
    "Otu00023", "Otu00025", "Otu00028", "Otu00029", "Otu00030", "Otu00035",
    "Otu00036", "Otu00037", "Otu00038", "Otu00039", "Otu00040", "Otu00047",
    "Otu00050", "Otu00052", "Otu00054", "Otu00055", "Otu00056", "Otu00060",
    "Otu00003|Otu00002|Otu00005|Otu00024|Otu00032|Otu00041|Otu00053",
    "Otu00014|Otu00021|Otu00017|Otu00031|Otu00057",
    "Otu00013|Otu00006", "Otu00026|Otu00001|Otu00034|Otu00048",
    "Otu00033|Otu00010",
    "Otu00042|Otu00004", "Otu00043|Otu00027|Otu00049", "Otu00051|Otu00045",
    "Otu00058|Otu00044", "Otu00059|Otu00046"
  )
)

# the function can show a progress bar if you have the progressr package installed
## optionally, specify the progress bar format

progressr::handlers(progressr::handler_progress(
  format = ":message :bar :percent | elapsed: :elapsed | eta: :eta",
  clear = FALSE,
  show_after = 0
))
## tell progressr to always report progress
progressr::handlers(global = TRUE)
## run the function and watch the live progress udpates
feat_imp <- get_feature_importance(results$trained_model,
  results$trained_model$trainingData, results$test_data,
  "dx",
  multiClassSummary, "AUC",
  class_probs = TRUE, method = "glmnet"
)

# you can specify any correlation method supported by `stats::cor`:
feat_imp <- get_feature_importance(results$trained_model,
  results$trained_model$trainingData, results$test_data,
  "dx",
  multiClassSummary, "AUC",
  class_probs = TRUE, method = "glmnet",
  corr_method = "pearson"
)
}

}
\author{
Begüm Topçuoğlu, \email{topcuoglu.begum@gmail.com}

Zena Lapp, \email{zenalapp@umich.edu}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_data.R
\name{preprocess_data}
\alias{preprocess_data}
\title{Preprocess data prior to running machine learning}
\usage{
preprocess_data(
  dataset,
  outcome_colname,
  method = c("center", "scale"),
  remove_var = "nzv",
  collapse_corr_feats = TRUE,
  to_numeric = TRUE,
  group_neg_corr = TRUE,
  prefilter_threshold = 1
)
}
\arguments{
\item{dataset}{Dataframe with an outcome variable and other columns as features.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{method}{Methods to preprocess the data, described in
\code{\link[caret:preProcess]{caret::preProcess()}} (default: \code{c("center","scale")}, use \code{NULL} for
no normalization).}

\item{remove_var}{Whether to remove variables with near-zero variance
(\code{'nzv'}; default), zero variance (\code{'zv'}), or none (\code{NULL}).}

\item{collapse_corr_feats}{Whether to keep only one of perfectly correlated
features.}

\item{to_numeric}{Whether to change features to numeric where possible.}

\item{group_neg_corr}{Whether to group negatively correlated features
together (e.g. c(0,1) and c(1,0)).}

\item{prefilter_threshold}{Remove features which only have non-zero & non-NA
values N rows or fewer (default: 1). Set this to -1 to keep all columns at
this step. This step will also be skipped if \code{to_numeric} is set to
\code{FALSE}.}
}
\value{
Named list including:
\itemize{
\item \code{dat_transformed}: Preprocessed data.
\item \code{grp_feats}: If features were grouped together, a named list of the features corresponding to each group.
\item \code{removed_feats}: Any features that were removed during preprocessing (e.g. because there was zero variance or near-zero variance for those features).
}

If the \code{progressr} package is installed, a progress bar with time elapsed
and estimated time to completion can be displayed.
}
\description{
Function to preprocess your data for input into \code{\link[=run_ml]{run_ml()}}.
}
\section{More details}{


See the \href{http://www.schlosslab.org/mikropml/articles/preprocess.html}{preprocessing vignette}
for more details.

Note that if any values in \code{outcome_colname} contain spaces, they will be
converted to underscores for compatibility with \code{caret}.
}

\examples{
preprocess_data(mikropml::otu_small, "dx")

# the function can show a progress bar if you have the progressr package installed
## optionally, specify the progress bar format
progressr::handlers(progressr::handler_progress(
  format = ":message :bar :percent | elapsed: :elapsed | eta: :eta",
  clear = FALSE,
  show_after = 0
))
## tell progressor to always report progress
\dontrun{
progressr::handlers(global = TRUE)
## run the function and watch the live progress udpates
dat_preproc <- preprocess_data(mikropml::otu_small, "dx")
}
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}

Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{get_performance_tbl}
\alias{get_performance_tbl}
\title{Get model performance metrics as a one-row tibble}
\usage{
get_performance_tbl(
  trained_model,
  test_data,
  outcome_colname,
  perf_metric_function,
  perf_metric_name,
  class_probs,
  method,
  seed = NA
)
}
\arguments{
\item{trained_model}{Trained model from \code{\link[caret:train]{caret::train()}}.}

\item{test_data}{Held out test data: dataframe of outcome and features.}

\item{outcome_colname}{Column name as a string of the outcome variable
(default \code{NULL}; the first column will be chosen automatically).}

\item{perf_metric_function}{Function to calculate the performance metric to
be used for cross-validation and test performance. Some functions are
provided by caret (see \code{\link[caret:postResample]{caret::defaultSummary()}}).
Defaults: binary classification = \code{twoClassSummary},
multi-class classification = \code{multiClassSummary},
regression = \code{defaultSummary}.}

\item{perf_metric_name}{The column name from the output of the function
provided to perf_metric_function that is to be used as the performance metric.
Defaults: binary classification = \code{"ROC"},
multi-class classification = \code{"logLoss"},
regression = \code{"RMSE"}.}

\item{class_probs}{Whether to use class probabilities (TRUE for categorical outcomes, FALSE for numeric outcomes).}

\item{method}{ML method.
Options: \code{c("glmnet", "rf", "rpart2", "svmRadial", "xgbTree")}.
\itemize{
\item glmnet: linear, logistic, or multiclass regression
\item rf: random forest
\item rpart2: decision tree
\item svmRadial: support vector machine
\item xgbTree: xgboost
}}

\item{seed}{Random seed (default: \code{NA}).
Your results will only be reproducible if you set a seed.}
}
\value{
A one-row tibble with columns \code{cv_auroc}, column for each of the performance metrics for the test data \code{method}, and \code{seed}.
}
\description{
Get model performance metrics as a one-row tibble
}
\examples{
\dontrun{
results <- run_ml(otu_small, "glmnet", kfold = 2, cv_times = 2)
names(results$trained_model$trainingData)[1] <- "dx"
get_performance_tbl(results$trained_model, results$test_data,
  "dx",
  multiClassSummary, "AUC",
  class_probs = TRUE,
  method = "glmnet"
)
}

}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}

Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_val.R
\name{keep_groups_in_cv_partitions}
\alias{keep_groups_in_cv_partitions}
\title{Whether groups can be kept together in partitions during cross-validation}
\usage{
keep_groups_in_cv_partitions(groups, group_partitions, kfold)
}
\arguments{
\item{groups}{Vector of groups to keep together when splitting the data into
train and test sets. If the number of groups in the training set is larger
than \code{kfold}, the groups will also be kept together for cross-validation.
Length matches the number of rows in the dataset (default: \code{NULL}).}

\item{group_partitions}{Specify how to assign \code{groups} to the training and
testing partitions (default: \code{NULL}). If \code{groups} specifies that some
samples belong to group \code{"A"} and some belong to group \code{"B"}, then setting
\code{group_partitions = list(train = c("A", "B"), test = c("B"))} will result
in all samples from group \code{"A"} being placed in the training set, some
samples from \code{"B"} also in the training set, and the remaining samples from
\code{"B"} in the testing set. The partition sizes will be as close to
\code{training_frac} as possible. If the number of groups in the training set is
larger than \code{kfold}, the groups will also be kept together for
cross-validation.}

\item{kfold}{Fold number for k-fold cross-validation (default: \code{5}).}
}
\value{
\code{TRUE} if possible, \code{FALSE} otherwise
}
\description{
Whether groups can be kept together in partitions during cross-validation
}
\author{
Kelly Sovacool, \email{sovacool@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{otu_mini_cont_results_nocv}
\alias{otu_mini_cont_results_nocv}
\title{Results from running the pipeline with glmnet on \code{otu_mini_bin} with \code{Otu00001}
as the outcome column,
using a custom train control scheme that does not perform cross-validation}
\format{
An object of class \code{list} of length 4.
}
\usage{
otu_mini_cont_results_nocv
}
\description{
Results from running the pipeline with glmnet on \code{otu_mini_bin} with \code{Otu00001}
as the outcome column,
using a custom train control scheme that does not perform cross-validation
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{get_perf_metric_name}
\alias{get_perf_metric_name}
\title{Get default performance metric name}
\usage{
get_perf_metric_name(outcome_type)
}
\arguments{
\item{outcome_type}{Type of outcome (one of: \code{"continuous"},\code{"binary"},\code{"multiclass"}).}
}
\value{
Performance metric name.
}
\description{
Get default performance metric name for cross-validation.
}
\examples{
get_perf_metric_name("continuous")
get_perf_metric_name("binary")
get_perf_metric_name("multiclass")
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{combine_hp_performance}
\alias{combine_hp_performance}
\title{Combine hyperparameter performance metrics for multiple train/test splits}
\usage{
combine_hp_performance(trained_model_lst)
}
\arguments{
\item{trained_model_lst}{List of trained models.}
}
\value{
Named list:
\itemize{
\item \code{dat}: Dataframe of performance metric for each group of hyperparameters
\item \code{params}: Hyperparameters tuned.
\item \code{Metric}: Performance metric used.
}
}
\description{
Combine hyperparameter performance metrics for multiple train/test splits generated by, for instance, \href{http://www.schlosslab.org/mikropml/articles/parallel.html}{looping in R} or using a \href{https://github.com/SchlossLab/mikropml-snakemake-workflow}{snakemake workflow} on a high-performance computer.
}
\examples{
\dontrun{
results <- lapply(seq(100, 102), function(seed) {
  run_ml(otu_small, "glmnet", seed = seed, cv_times = 2, kfold = 2)
})
models <- lapply(results, function(x) x$trained_model)
combine_hp_performance(models)
}
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_metrics.R
\name{get_perf_metric_fn}
\alias{get_perf_metric_fn}
\title{Get default performance metric function}
\usage{
get_perf_metric_fn(outcome_type)
}
\arguments{
\item{outcome_type}{Type of outcome (one of: \code{"continuous"},\code{"binary"},\code{"multiclass"}).}
}
\value{
Performance metric function.
}
\description{
Get default performance metric function
}
\examples{
get_perf_metric_fn("continuous")
get_perf_metric_fn("binary")
get_perf_metric_fn("multiclass")
}
\author{
Zena Lapp, \email{zenalapp@umich.edu}
}
