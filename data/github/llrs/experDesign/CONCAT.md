
<!-- README.md is generated from README.Rmd. Please edit that file -->

# experDesign

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/experDesign)](https://CRAN.R-project.org/package=experDesign)
[![R build
status](https://github.com/llrs/experDesign/workflows/R-CMD-check/badge.svg)](https://github.com/llrs/experDesign/actions?workflow=R-CMD-check)
[![Coverage
status](https://codecov.io/gh/llrs/experDesign/branch/master/graph/badge.svg)](https://codecov.io/github/llrs/experDesign?branch=master)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03358/status.svg)](https://doi.org/10.21105/joss.03358)

<!-- badges: end -->

The goal of experDesign is to help you decide which samples go in which
batch, reducing the potential batch bias before performing an
experiment. It provides three main functions :

-   `design()`: Randomize the samples according to their variables.
-   `replicates()`: Selects some samples for replicates and randomizes
    the samples.
-   `spatial()`: Randomize the samples on a spatial grid.

## Installation

To install the latest version on
[CRAN](https://CRAN.R-project.org/package=experDesign) use:

``` r
install.packages("experDesign")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("llrs/experDesign")
```

## Example

We can use the survey dataset for the examples:

``` r
library("experDesign")
data(survey, package = "MASS") 
head(survey)
#>      Sex Wr.Hnd NW.Hnd W.Hnd    Fold Pulse    Clap Exer Smoke Height      M.I
#> 1 Female   18.5   18.0 Right  R on L    92    Left Some Never 173.00   Metric
#> 2   Male   19.5   20.5  Left  R on L   104    Left None Regul 177.80 Imperial
#> 3   Male   18.0   13.3 Right  L on R    87 Neither None Occas     NA     <NA>
#> 4   Male   18.8   18.9 Right  R on L    NA Neither None Never 160.00   Metric
#> 5   Male   20.0   20.0 Right Neither    35   Right Some Never 165.00   Metric
#> 6 Female   18.0   17.7 Right  L on R    64   Right Some Never 172.72 Imperial
#>      Age
#> 1 18.250
#> 2 17.583
#> 3 16.917
#> 4 20.333
#> 5 23.667
#> 6 21.000
```

The dataset has numeric, categorical values and some `NA`’s value.

# Picking samples for each batch

Imagine that we can only work in groups of 70, and we want to randomize
by Sex, Smoke, Age, and by writing hand.  
There are 1.6543999^{61} combinations some of them would be have in a
single experiment all the right handed students. We could measure all
these combinations but we can try to find an optimum value.

``` r
# To reduce the variables used:
omit <- c("Wr.Hnd", "NW.Hnd", "Fold", "Pulse", "Clap", "Exer", "Height", "M.I")
(keep <- colnames(survey)[!colnames(survey) %in% omit])
#> [1] "Sex"   "W.Hnd" "Smoke" "Age"
head(survey[, keep])
#>      Sex W.Hnd Smoke    Age
#> 1 Female Right Never 18.250
#> 2   Male  Left Regul 17.583
#> 3   Male Right Occas 16.917
#> 4   Male Right Never 20.333
#> 5   Male Right Never 23.667
#> 6 Female Right Never 21.000

# Looking for groups at most of 70 samples.
index <- design(pheno = survey, size_subset = 70, omit = omit)
index
#> $SubSet1
#>  [1]   8  10  11  13  16  18  21  25  26  32  38  42  43  44  50  54  55  60  65
#> [20]  67  68  72  74  75  76  79  88  89  92  97 109 124 140 143 145 147 150 155
#> [39] 163 164 173 176 181 182 183 185 187 194 200 206 207 208 209 210 213 217 218
#> [58] 228 231 236
#> 
#> $SubSet2
#>  [1]   2   3   7   9  12  15  20  23  28  36  37  39  45  49  51  52  57  59  61
#> [20]  63  64  66  69  70  71  82  84  85  86  99 100 104 112 114 115 119 126 135
#> [39] 137 146 153 159 160 166 169 171 174 175 180 186 191 195 204 211 215 219 221
#> [58] 227 232
#> 
#> $SubSet3
#>  [1]   4   6  14  24  29  34  35  40  46  47  53  56  58  73  78  80  81  87  95
#> [20]  96  98 101 102 103 105 107 108 117 118 120 121 122 125 129 131 132 133 134
#> [39] 136 139 142 148 152 154 156 157 162 172 189 190 199 205 220 226 229 230 233
#> [58] 234 235
#> 
#> $SubSet4
#>  [1]   1   5  17  19  22  27  30  31  33  41  48  62  77  83  90  91  93  94 106
#> [20] 110 111 113 116 123 127 128 130 138 141 144 149 151 158 161 165 167 168 170
#> [39] 177 178 179 184 188 192 193 196 197 198 201 202 203 212 214 216 222 223 224
#> [58] 225 237
```

We can transform then into a vector to append to the file or to pass to
the lab mate with:

``` r
head(batch_names(index))
#> [1] "SubSet4" "SubSet2" "SubSet2" "SubSet3" "SubSet4" "SubSet3"
```

# Previous work

The CRAN task View of [Experimental
Design](https://CRAN.R-project.org/view=ExperimentalDesign) includes
many packages relevant for designing an experiment before collecting
data, but none of them provides how to manage them once the samples are
already collected.

Two packages allow to distribute the samples on batches:

-   The
    [OSAT](https://bioconductor.org/packages/release/bioc/html/OSAT.html)
    package handles categorical variables but not numeric data. It
    doesn’t work with our data.

-   The [minDiff](https://github.com/m-Py/minDiff) package reported in
    [Stats.SE](https://stats.stackexchange.com/a/326015/105234), handles
    both numeric and categorical data. But it can only optimize for two
    nominal criteria. It doesn’t work for our data.

-   The [Omixer](https://bioconductor.org/packages/Omixer/) package
    handles both numeric and categorical data (converting categorical
    variables to numeric). But both the same way either Pearson’s
    Chi-squared Test if there are few samples or Kendall’s correlation.
    It does allow to protect some spots from being used.

If you are still designing the experiment and do not have collected any
data [DeclareDesign](https://cran.r-project.org/package=DeclareDesign)
might be relevant for you.

Question in
[Bioinformatics.SE](https://bioinformatics.stackexchange.com/q/4765/48)
I made before developing the package.

# Other

Please note that this project is released with a [Contributor Code of
Conduct](https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/).
By participating in this project you agree to abide by its terms.
# experDesign 0.1.1

* Added thesis advisors on the description.

* Update documentation

# experDesign 0.1.0

* Added reference to a new package Omixer on the README. 

* Fixed batches sizes errors.

* Speed increase (5x) on design, spatial and replicates.

* Update Code of Conduct.

* Add online documentation url.

# experDesign 0.0.4

* Remove BiocStyle dependency.
* Gain the ability to name the subsets.
* Add examples to all functions.
* Add function to consider spatial distribution on plates/machines.

# experDesign 0.0.2

* CRAN release

# experDesign 0.0.1

# experDesign 0.0.900

# experDesign 0.0.100

# experDesign 0.0.0.9000

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
reported to the community leaders responsible for enforcement at package's 
maintainer email. All complaints will be reviewed and investigated promptly and 
fairly.

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
This resubmission fixes some bugs and speed increase on some functions

## Test environments

* local R installation, R 3.6.3 and 4.0.1
* ubuntu 16.04 (on travis-ci), old-release, release, devel
* win-builder, old-release, release, devel
* Github Actions (windows), release
* Github Actions (macOS), release
* Github Actions Ubuntu 20.04, release and devel

## R CMD check results

0 errors | 0 warnings | 0 note

## Dependencies

No dependencies
Tests and Coverage
================
28 noviembre, 2018 17:05:00

-   [Coverage](#coverage)
-   [Unit Tests](#unit-tests)

This output is created by [covrpage](https://github.com/yonicd/covrpage).

Coverage
--------

Coverage summary is created using the [covr](https://github.com/r-lib/covr) package.

| Object                                             | Coverage (%) |
|:---------------------------------------------------|:------------:|
| experDesign                                        |     85.49    |
| [R/reporting.R](../R/reporting.R)                  |     0.00     |
| [R/evaluate\_category.R](../R/evaluate_category.R) |     36.84    |
| [R/utils.R](../R/utils.R)                          |     83.33    |
| [R/indexing.R](../R/indexing.R)                    |     89.29    |
| [R/QC.R](../R/QC.R)                                |     92.86    |
| [R/designer.R](../R/designer.R)                    |     96.67    |
| [R/entropy.R](../R/entropy.R)                      |    100.00    |
| [R/evaluate\_num.R](../R/evaluate_num.R)           |    100.00    |
| [R/evaluate.R](../R/evaluate.R)                    |    100.00    |

<br>

Unit Tests
----------

Unit Test summary is created using the [testthat](https://github.com/r-lib/testthat) package.

| file                                                       |    n|   time|  error|  failed|  skipped|  warning|
|:-----------------------------------------------------------|----:|------:|------:|-------:|--------:|--------:|
| [test\_batch-names.R](testthat/test_batch-names.R)         |    1|  0.002|      0|       0|        0|        0|
| [test\_create-subset.R](testthat/test_create-subset.R)     |    5|  0.005|      0|       0|        0|        0|
| [test\_design.R](testthat/test_design.R)                   |    2|  0.172|      0|       0|        0|        0|
| [test\_entropy.R](testthat/test_entropy.R)                 |    3|  0.003|      0|       0|        0|        0|
| [test\_evaluate-helper.R](testthat/test_evaluate-helper.R) |    2|  0.002|      0|       0|        0|        0|
| [test\_evaluate-index.R](testthat/test_evaluate-index.R)   |    1|  0.007|      0|       0|        0|        0|
| [test\_evaluate-mad.R](testthat/test_evaluate-mad.R)       |    1|  0.001|      0|       0|        0|        0|
| [test\_evaluate-mean.R](testthat/test_evaluate-mean.R)     |    1|  0.002|      0|       0|        0|        0|
| [test\_evaluate-na.R](testthat/test_evaluate-na.R)         |    1|  0.001|      0|       0|        0|        0|
| [test\_evaluate-orig.R](testthat/test_evaluate-orig.R)     |    6|  0.006|      0|       0|        0|        0|
| [test\_evaluate-sd.R](testthat/test_evaluate-sd.R)         |    1|  0.002|      0|       0|        0|        0|
| [test-extreme\_cases.R](testthat/test-extreme_cases.R)     |    4|  3.182|      0|       0|        0|        0|
| [test\_insert.R](testthat/test_insert.R)                   |    3|  0.003|      0|       0|        0|        0|
| [test-qc.R](testthat/test-qc.R)                            |    3|  0.003|      0|       0|        0|        0|
| [test\_replicates.R](testthat/test_replicates.R)           |    1|  0.223|      0|       0|        0|        0|
| [test\_simplify2matrix.R](testthat/test_simplify2matrix.R) |    1|  0.001|      0|       0|        0|        0|

<details closed> <summary> Show Detailed Test Results </summary>

| file                                                            | context             | test                 | status |    n|   time|
|:----------------------------------------------------------------|:--------------------|:---------------------|:-------|----:|------:|
| [test\_batch-names.R](testthat/test_batch-names.R#L6)           | batch\_names        | works                | PASS   |    1|  0.002|
| [test\_create-subset.R](testthat/test_create-subset.R#L5)       | create\_subset      | works                | PASS   |    2|  0.002|
| [test\_create-subset.R](testthat/test_create-subset.R#L10)      | create\_subset      | works well           | PASS   |    1|  0.001|
| [test\_create-subset.R](testthat/test_create-subset.R#L16)      | create\_subset      | use\_index           | PASS   |    2|  0.002|
| [test\_design.R](testthat/test_design.R#L8)                     | design              | works                | PASS   |    2|  0.172|
| [test\_entropy.R](testthat/test_entropy.R#L5)                   | entropy             | Extremes             | PASS   |    2|  0.002|
| [test\_entropy.R](testthat/test_entropy.R#L12)                  | entropy             | Ignores NA           | PASS   |    1|  0.001|
| [test\_evaluate-helper.R](testthat/test_evaluate-helper.R#L7)   | evaluate\_helper    | works                | PASS   |    1|  0.001|
| [test\_evaluate-helper.R](testthat/test_evaluate-helper.R#L15)  | evaluate\_helper    | mean                 | PASS   |    1|  0.001|
| [test\_evaluate-index.R](testthat/test_evaluate-index.R#L9_L13) | evaluate\_index     | works                | PASS   |    1|  0.007|
| [test\_evaluate-mad.R](testthat/test_evaluate-mad.R#L9)         | evaluate\_mad       | works                | PASS   |    1|  0.001|
| [test\_evaluate-mean.R](testthat/test_evaluate-mean.R#L9)       | evaluate\_mean      | works                | PASS   |    1|  0.002|
| [test\_evaluate-na.R](testthat/test_evaluate-na.R#L10)          | evaluate\_na        | works                | PASS   |    1|  0.001|
| [test\_evaluate-orig.R](testthat/test_evaluate-orig.R#L8)       | evaluate\_orig      | works                | PASS   |    6|  0.006|
| [test\_evaluate-sd.R](testthat/test_evaluate-sd.R#L9)           | evaluate\_sd        | works                | PASS   |    1|  0.002|
| [test-extreme\_cases.R](testthat/test-extreme_cases.R#L7_L10)   | test-extreme\_cases | extreme\_cases works | PASS   |    3|  3.167|
| [test-extreme\_cases.R](testthat/test-extreme_cases.R#L19)      | test-extreme\_cases | check\_index works   | PASS   |    1|  0.015|
| [test\_insert.R](testthat/test_insert.R#L8)                     | insert              | works                | PASS   |    3|  0.003|
| [test-qc.R](testthat/test-qc.R#L6)                              | qcSubset            | all                  | PASS   |    1|  0.001|
| [test-qc.R](testthat/test-qc.R#L12)                             | qcSubset            | by batch             | PASS   |    2|  0.002|
| [test\_replicates.R](testthat/test_replicates.R#L7)             | replicates          | works                | PASS   |    1|  0.223|
| [test\_simplify2matrix.R](testthat/test_simplify2matrix.R#L5)   | simplify2matrix     | works                | PASS   |    1|  0.001|

</details>

<details> <summary> Session Info </summary>

| Field    | Value                        |
|:---------|:-----------------------------|
| Version  | R version 3.5.1 (2018-07-02) |
| Platform | i686-pc-linux-gnu (32-bit)   |
| Running  | Ubuntu 16.04.5 LTS           |
| Language | en\_US                       |
| Timezone | Europe/Madrid                |

| Package  | Version |
|:---------|:--------|
| testthat | 2.0.1   |
| covr     | 3.2.1   |
| covrpage | 0.0.62  |

</details>

<!--- Final Status : pass --->
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# experDesign

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/experDesign)](https://CRAN.R-project.org/package=experDesign)
[![R build status](https://github.com/llrs/experDesign/workflows/R-CMD-check/badge.svg)](https://github.com/llrs/experDesign/actions?workflow=R-CMD-check)
[![Coverage status](https://codecov.io/gh/llrs/experDesign/branch/master/graph/badge.svg)](https://codecov.io/github/llrs/experDesign?branch=master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/142569201.svg)](https://zenodo.org/badge/latestdoi/142569201)


<!-- badges: end -->

The goal of experDesign is to help you decide which samples go in which batch, 
reducing the potential batch bias before performing an experiment.
It provides three main functions :

* `design()`: Randomize the samples according to their variables.
* `replicates()`: Selects some samples for replicates and randomizes the samples.
* `spatial()`: Randomize the samples on a spatial grid.

## Installation

To install the latest version on [CRAN](https://CRAN.R-project.org/package=experDesign) use:

```r
install.packages("experDesign")
```
You can install the development version from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("llrs/experDesign")
```

## Example

We can use the survey dataset for the examples:

```{r show}
library("experDesign")
data(survey, package = "MASS") 
head(survey)
```

The dataset has numeric, categorical values and some `NA`'s value.

# Picking samples for each batch

Imagine that we can only work in groups of 70, and we want to randomize by Sex, 
Smoke, Age, and by writing hand.  
There are `r choose(237, 70)` combinations some of them would be have in a 
single experiment all the right handed students. We could measure all these combinations
but we can try to find an optimum value.

```{r design, fig.show='hold'}
# To reduce the variables used:
omit <- c("Wr.Hnd", "NW.Hnd", "Fold", "Pulse", "Clap", "Exer", "Height", "M.I")
(keep <- colnames(survey)[!colnames(survey) %in% omit])
head(survey[, keep])

# Looking for groups at most of 70 samples.
index <- design(pheno = survey, size_subset = 70, omit = omit)
index
```

We can transform then into a vector to append to the file or to pass to the lab mate with:

```{r batch_names}
head(batch_names(index))
```



# Previous work

The CRAN task View of [Experimental Design](https://CRAN.R-project.org/view=ExperimentalDesign) includes many packages relevant for designing an experiment before collecting data, but none of them provides how to manage them once the samples are already collected.

Two packages allow to distribute the samples on batches:

- The [OSAT](https://bioconductor.org/packages/release/bioc/html/OSAT.html) package handles categorical 
variables but not numeric data. It doesn't work with our data.

 - The [minDiff](https://github.com/m-Py/minDiff) package reported in [Stats.SE](https://stats.stackexchange.com/a/326015/105234), handles both 
numeric and categorical data. But it can only optimize for two nominal criteria.
It doesn't work for our data.

 - The [Omixer](https://bioconductor.org/packages/Omixer/) package handles both 
numeric and categorical data (converting categorical variables to numeric). But both the same way either Pearson's Chi-squared Test if there are few samples or Kendall's correlation. It does allow to protect some spots from being used.

If you are still designing the experiment and do not have collected any data [DeclareDesign](https://cran.r-project.org/package=DeclareDesign) might be relevant for you.

Question in [Bioinformatics.SE](https://bioinformatics.stackexchange.com/q/4765/48) I made before developing the package.

# Other

Please note that this project is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/).
By participating in this project you agree to abide by its terms.
---
title: "exprDesign"
author: "Lluís Revilla Sancho"
date: "`r Sys.Date()`"
output:
  html_document:
    fig_caption: true
    code_folding: show
    self_contained: yes
    toc: true
    toc_float: true
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{exprDesign}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r knitsetup, message=FALSE, warning=FALSE, include=FALSE}
knitr::opts_knit$set(root.dir = ".")
knitr::opts_chunk$set(collapse = TRUE, warning = TRUE)
set.seed(445)
library("experDesign")
```

# Introduction

This package was developed to help prepare some samples to be send to a 
facility. It assumes that you have collected the samples and the information 
but you still need to do the experiment in several batches due to technical or 
practical limitations. The question that tries to answer is: 

> Which samples go with each batch?

Of all the possible combinations of samples, it looks for the combination which 
minimizes the differences between each subgroup according to the following rules:

 - If the variable is categorical it tires to randomize the variable across the subgroups.
 - If the variable is numeric it tries to distribute evenly following the original distribution of values.
 - If there are `NA` (not available values) it looks to distribute them randomly.

Even with this measures you might end up with some batch effect due to:
 - Confounding variables not provided for their randomization on the batches
  Sometimes due to being unknown, impossible to measure.
 - Lack of [replicates](https://en.wikipedia.org/wiki/Replication_(statistics))(samples with the same conditions)
  If you can't provide new replicates, aim to provide more technical replicates.
  Technical replicates mean reanalyzing the same sample twice or more, the more samples with technical replicates the more accurate your measures will be and easier to avoid or detect batch effects.
 - Processing
  If there is a change on the methodology, or you pause and resume later the sample collection there might be changes on the outcome due to external factors. 

# Previous work

Before building this package I would like to give credit to those that made 
also efforts in this direction:


The CRAN task View of [Experimental Design](https://CRAN.R-project.org/view=ExperimentalDesign) includes many packages relevant for designing an experiment before collecting data, but none of them provides how to manage them once the samples are already collected.

Two packages allow to distribute the samples on batches:

- The [OSAT](https://bioconductor.org/packages/release/bioc/html/OSAT.html) package handles categorical 
variables but not numeric data. It doesn't work with our data.

 - The [minDiff](https://github.com/m-Py/minDiff) package reported in [Stats.SE](https://stats.stackexchange.com/a/326015/105234), handles both 
numeric and categorical data. But it can only optimize for two nominal criteria.
It doesn't work for our data.

 - The [Omixer](https://bioconductor.org/packages/Omixer/) package handles both 
numeric and categorical data (converting categorical variables to numeric). But both the same way either Pearson's Chi-squared Test if there are few samples or Kendall's correlation. It does allow to protect some spots from being used.

If you are still designing the experiment and do not have collected any data [DeclareDesign](https://cran.r-project.org/package=DeclareDesign) might be relevant for you.

Question in [Bioinformatics.SE](https://bioinformatics.stackexchange.com/q/4765/48) I made before developing the package.


# Design of Experiment {#DoE}

Imagine you have some samples already collected and you want to distributed them in batches:
```{r experDesign_setup}
library("experDesign")
metadata <- expand.grid(height = seq(60, 80, 5), 
                        weight = seq(100, 300, 50),
                        sex = c("Male","Female"))
head(metadata, 15)
```

If you block incorrectly and end up with a group in a single batch we will end up with batch effect.
In order to avoid this `design()` helps you assign each sample to a batch (in this case each batch has 24 samples at most). First we can explore the number of samples and the number of batches:

```{r size}
size_data <- nrow(metadata)
size_batch <- 24
(batches <- optimum_batches(size_data, size_batch))
# So now the best number of samples for each batch is less than the available
(size <- optimum_subset(size_data, batches))
# The distribution of samples per batch
sizes_batches(size_data, size, batches)
```

Note that instead of using a whole batch and then leave a single sample on the third distributes all the samples in the three batches that will be needed.

# Randomization

We can directly look for the distribution of the samples given our max number of samples per batch:

```{r design}
d <- design(metadata, size_batch)
# It is a list but we can convert it to a vector with:
batch_names(d)
```

Naively one would either fill some batches fully or distribute them not evenly 
(the  first 17 packages together, the next 17 and so on). This solution ensures 
that the data is randomized. For more random distribution you can increase the number of iterations performed to calculate this distribution.

# Randomization and replicates

If you need space for replicates to control for batch effect you can use:

```{r replicates}
r <- replicates(metadata, size_batch, 5)
lengths(r)
r
```

Which seeks as controls the most diverse values and adds them to the samples 
distribution. Note that if the sample is already present on that batch is not added again, that's why the number of samples per batch is different from the design without replicates.

# Layout

Lastly, we can see how these samples would be distributed in a layout of 6x4:

```{r spatial}
s <- spatial(r, metadata, rows = LETTERS[1:6], columns = 1:4)
head(s)
```

# Report for easy on field usage

We can add the batches to the initial data with `inspect()`:

```{r report}
report <- inspect(r, metadata)
report2 <- inspect(s, report, index_name = "position")
head(report2)
```

And now we can see the batch and position of each sample


# Unbalanced setting

In the previous case the data was mostly balanced (check it out in the `orig` object) 
but let's create an unbalanced dataset to check it.

```{r unbalanced}
n <- 99
samples <- 100
unbalanced <- data.frame(Classroom = rep(c("A", "B"), each = samples/2),
                         Sex = c(rep("M", n), rep("F", samples-n)),
                         Age = rnorm(samples, mean = 25, sd = 3))
table(unbalanced)[, , 1:5]
```

In this dataset the classroom a single classroom has all the females (`r 50 -n`).

```{r unbalanced_design}
i <- design(unbalanced, 15)

# Mean entropy en each subset
rowMeans(evaluate_index(i, unbalanced)["entropy", , ])
# Original entropy on the dataset
evaluate_orig(unbalanced)["entropy", ]
# Dispersion of the entropy
apply(evaluate_index(i, unbalanced)["entropy", , ], 1, sd)
```

We can see that in this simple case where a single variable has all the other cases we approximately reached the same entropy levels. 

# Quality check

If you need a subset with the samples that are more diverse you can use the 
following function:

```{r QC}
data(survey, package = "MASS") 
head(survey)
samples <- extreme_cases(survey, size = 10)
survey[samples, ]
```

You can also test a given index with `check_index()`:

```{r check_index}
check_index(unbalanced, i)
```

Each row has information about how accurate is a given variable to the samples available (on this case `unbalanced`). 
Some variables are distributed more randomly than others on this index.

If we are not satisfied we could `design()` a new index increasing the iterations to obtain a potentially better distribution. 
If you want a good stratified randomization you should increase the iterations used 10 fold. 

# SessionInfo

```{r sessioninfo}
sessionInfo()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{qcSubset}
\alias{qcSubset}
\title{Random subset}
\usage{
qcSubset(index, size, each = FALSE)
}
\arguments{
\item{index}{A list of indices indicating which samples go to which subset.}

\item{size}{The number of samples that should be taken.}

\item{each}{A logical value if the subset should be taken from all the
samples or for each batch.}
}
\description{
Select randomly some samples from an index
}
\examples{
set.seed(50)
index <- create_subset(100, 50, 2)
QC_samples <- qcSubset(index, 10)
QC_samplesBatch <- qcSubset(index, 10, TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{extreme_cases}
\alias{extreme_cases}
\title{Select the subset of extreme cases to evaluation}
\usage{
extreme_cases(pheno, size, omit = NULL, iterations = 500)
}
\arguments{
\item{pheno}{Data.frame with the sample information.}

\item{size}{The number of samples to subset.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}

\item{iterations}{Numeric value of iterations that will be performed.}
}
\value{
A vector with the number of the rows that are selected.
}
\description{
Subset some samples that are mostly different.
}
\examples{
metadata <- expand.grid(height = seq(60, 80, 5), weight = seq(100, 300, 50),
 sex = c("Male","Female"))
sel <- extreme_cases(metadata, 10)
# We can see that it selected both Female and Males and wide range of height
# and weight:
metadata[sel, ]
}
\seealso{
\code{\link[=optimum]{optimum()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporting.R
\name{inspect}
\alias{inspect}
\title{Inspect the index}
\usage{
inspect(i, pheno, omit = NULL, index_name = "batch")
}
\arguments{
\item{i}{List of indices of samples per batch}

\item{pheno}{Data.frame with the sample information.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}

\item{index_name}{Column name of the index of the resulting data.frame.}
}
\value{
The data.frame with a new column batch with the name of the batch the sample goes to.
}
\description{
Given the index and the data of the samples append the batch assignment
}
\examples{
data(survey, package = "MASS")
columns <- c("Sex", "Age", "Smoke")
index <- design(pheno = survey[, columns], size_subset = 70,
                iterations = 10)
batches <- inspect(index, survey[, columns])
head(batches)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_num.R
\name{evaluate_mad}
\alias{evaluate_mad}
\title{Evaluate median absolute deviation}
\usage{
evaluate_mad(i, pheno)
}
\arguments{
\item{i}{List of indices}

\item{pheno}{Data.frame with information about the samples}
}
\value{
A vector with the mean difference between the median absolute deviation
of each group and the original mad.
}
\description{
Looks for the median absolute deviation values in each subgroup.
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
# Note that categorical columns will be omitted:
evaluate_mad(index, survey[, c("Sex", "Smoke", "Age")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}

Other functions to evaluate numbers: 
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_sd}()}
}
\concept{functions to evaluate numbers}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indexing.R
\name{create_subset}
\alias{create_subset}
\title{Create index of subsets of a data}
\usage{
create_subset(size_data, size_subset = NULL, n = NULL, name = "SubSet")
}
\arguments{
\item{size_data}{A numeric value of the amount of samples to distribute.}

\item{size_subset}{A numeric value with the amount of samples per batch.}

\item{n}{A numeric value with the number of batches.}

\item{name}{A character used to name the subsets, either a single one or a
vector the same size as \code{n}.}
}
\value{
A random list of indices of the samples.
}
\description{
Index of the samples grouped by batches.
}
\examples{
index <- create_subset(100, 50, 2)
}
\seealso{
\code{\link[=batch_names]{batch_names()}}, \code{\link[=use_index]{use_index()}} if you already
have a factor to be used as index.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{spatial}
\alias{spatial}
\title{Distribute the sample on the plate}
\usage{
spatial(
  index,
  pheno,
  omit = NULL,
  remove_positions = NULL,
  rows = LETTERS[1:5],
  columns = 1:10,
  iterations = 500
)
}
\arguments{
\item{index}{A list with the samples on each subgroup, as provided from
\code{design()} or \code{replicates()}.}

\item{pheno}{Data.frame with the sample information.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}

\item{remove_positions}{Character, name of positions.}

\item{rows}{Character, name of the rows to be used.}

\item{columns}{Character, name of the rows to be used.}

\item{iterations}{Numeric value of iterations that will be performed.}
}
\value{
The indices of which samples go with which batch.
}
\description{
This function assumes that to process the batch the samples are distributes in
a plate. Sometimes you know in advance the
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 25)
index2 <- spatial(index, survey[, c("Sex", "Smoke", "Age")], iterations = 25)
head(index2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experDesign-package.R
\docType{package}
\name{experDesign-package}
\alias{experDesign-package}
\alias{experDesign}
\title{experDesign: Expert experiment design in batches}
\description{
Enables easy distribution of samples per batch avoiding batch and
confounding effects by randomization of the variables in each batch.
}
\details{
The most important function is \code{\link[=design]{design()}}, which distributes
samples in batches according to the information provided.

To help in the bench there is the \code{\link[=inspect]{inspect()}} function that appends
the group to the data provided.
}
\author{
Lluís Revilla
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_num.R
\name{evaluate_sd}
\alias{evaluate_sd}
\title{Evaluates the mean of the numeric values}
\usage{
evaluate_sd(i, pheno)
}
\arguments{
\item{i}{List of indices}

\item{pheno}{Data.frame with the samples}
}
\value{
A matrix with the standard deviation value for each column for each
subset
}
\description{
Looks for the standard deviation of the numeric values
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
# Note that categorical columns will be omitted:
evaluate_sd(index, survey[, c("Sex", "Smoke", "Age")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()}

Other functions to evaluate numbers: 
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()}
}
\concept{functions to evaluate numbers}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_category.R
\name{evaluate_entropy}
\alias{evaluate_entropy}
\title{Evaluate entropy}
\usage{
evaluate_entropy(i, pheno)
}
\arguments{
\item{i}{list of numeric indices of the data.frame}

\item{pheno}{Data.frame with information about the samples}
}
\value{
Value to minimize
}
\description{
Looks if the nominal or character columns are equally distributed according
to the entropy and taking into account the independence between batches.
If any column is different in each row it is assumed to be the sample names
and thus omitted.
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
# Note that numeric columns will be omitted:
evaluate_entropy(index, survey[, c("Sex", "Smoke", "Age")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}

Other functions to evaluate categories: 
\code{\link{evaluate_independence}()},
\code{\link{evaluate_na}()}
}
\concept{functions to evaluate categories}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy.R
\name{evaluate_na}
\alias{evaluate_na}
\title{Evaluate the dispersion of NAs}
\usage{
evaluate_na(i, pheno)
}
\arguments{
\item{i}{list of numeric indices of the data.frame}

\item{pheno}{Data.frame}
}
\value{
The optimum value to reduce
}
\description{
Looks how are \code{NA} distributed in each subset
}
\examples{
samples <- 10
m <- matrix(rnorm(samples), nrow = samples)
m[sample(seq_len(samples), size = 5), ] <- NA # Some NA
i <- create_subset(samples, 3, 4) # random subsets
evaluate_na(i, m)
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}

Other functions to evaluate categories: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()}

Other functions to evaluate numbers: 
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_sd}()}
}
\concept{functions to evaluate categories}
\concept{functions to evaluate numbers}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_category.R
\name{evaluate_independence}
\alias{evaluate_independence}
\title{Compare independence by chisq.test}
\usage{
evaluate_independence(i, pheno)
}
\arguments{
\item{i}{Index of subsets.}

\item{pheno}{A data.frame with the information about the samples.}
}
\value{
Returns a vector with the p-values of the chisq.test between the
category and the subset.
}
\description{
Looks the independence between the categories and the batches.
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
# Note that numeric columns will be omitted:
evaluate_independence(index, survey[, c("Sex", "Smoke", "Age")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}

Other functions to evaluate categories: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_na}()}
}
\concept{functions to evaluate categories}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indexing.R
\name{batch_names}
\alias{batch_names}
\title{Name the batch}
\usage{
batch_names(i)
}
\arguments{
\item{i}{A list of numeric indices.}
}
\value{
A character vector with the names of the batch for each the index.
}
\description{
Given an index return the name of the batches the samples are in
}
\examples{
index <- create_subset(100, 50, 2)
batch <- batch_names(index)
head(batch)
}
\seealso{
\code{\link[=create_subset]{create_subset()}}, for the inverse look at
\code{\link[=use_index]{use_index()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/designer.R
\name{replicates}
\alias{replicates}
\title{Design a batch experiment with experimental controls}
\usage{
replicates(pheno, size_subset, controls, omit = NULL, iterations = 500)
}
\arguments{
\item{pheno}{Data.frame with the sample information.}

\item{size_subset}{Numeric value of the number of sample per batch.}

\item{controls}{The numeric value of the amount of technical controls per
batch.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}

\item{iterations}{Numeric value of iterations that will be performed.}
}
\value{
A index with some samples duplicated in the batches
}
\description{
To ensure that the batches are comparable some samples are processed in each
batch. This function allows to take into account that effect.
It uses the most different samples as controls as defined with \code{\link[=extreme_cases]{extreme_cases()}}.
}
\examples{
samples <- data.frame(L = letters[1:25], Age = rnorm(25))
index <- replicates(samples, 5, controls = 2, iterations = 10)
head(index)
}
\seealso{
\code{\link[=design]{design()}}, \code{\link[=extreme_cases]{extreme_cases()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{check_index}
\alias{check_index}
\title{Check index distribution on batches}
\usage{
check_index(pheno, index, omit = NULL)
}
\arguments{
\item{pheno}{Data.frame with the sample information.}

\item{index}{A list of indices indicating which samples go to which subset.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}
}
\value{
A matrix with the differences with the original data.
}
\description{
Report the statistics for each subset and variable compared to the original.
}
\details{
The closer the values are to 0, the less difference is with the original
distribution, so it is a better randomization.
}
\examples{
index <- create_subset(50, 24)
metadata <- expand.grid(height = seq(60, 80, 5), weight = seq(100, 300, 50),
                        sex = c("Male","Female"))
check_index(metadata, index)
}
\seealso{
Functions that create an index \code{\link[=design]{design()}}, \code{\link[=replicates]{replicates()}},
\code{\link[=spatial]{spatial()}}. See also \code{\link[=create_subset]{create_subset()}} for a random index.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimum.R
\name{optimum}
\alias{optimum}
\alias{optimum_batches}
\alias{optimum_subset}
\alias{sizes_batches}
\title{Optimum values for batches}
\usage{
optimum_batches(size_data, size_subset)

optimum_subset(size_data, batches)

sizes_batches(size_data, size_subset, batches)
}
\arguments{
\item{size_data}{A numeric value of the number of samples to use.}

\item{size_subset}{Numeric value of the number of sample per batch.}

\item{batches}{A numeric value of the number of batches.}
}
\value{
\describe{
\item{\code{optimum_batches}}{A numeric value with the number of batches to
use.}
\item{\code{optimum_subset}}{A numeric value with the maximum number of samples
per batch of the data.}
\item{\code{sizes_batches}}{A numeric vector with the number of samples in each
batch.}
}
}
\description{
Calculates the optimum values for number of batches or size of the batches.
If you need to do several batches it can be better to distribute it evenly
and add replicates.
}
\examples{
size_data <- 50
size_batch <- 24
(batches <- optimum_batches(size_data, size_batch))
# So now the best number of samples for each batch is less than the available
(size <- optimum_subset(size_data, batches))
# The distribution of samples per batch
sizes_batches(size_data, size, batches)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_num.R
\name{evaluate_mean}
\alias{evaluate_mean}
\title{Evaluates the mean of the numeric values}
\usage{
evaluate_mean(i, pheno)
}
\arguments{
\item{i}{List of indices}

\item{pheno}{Data.frame with information about the samples}
}
\value{
A matrix with the mean value for each column for each subset
}
\description{
Looks for the mean of the numeric values
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
# Note that categorical columns will be omitted:
evaluate_mean(index, survey[, c("Sex", "Smoke", "Age")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}

Other functions to evaluate numbers: 
\code{\link{evaluate_mad}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_sd}()}
}
\concept{functions to evaluate numbers}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate.R
\name{evaluate_orig}
\alias{evaluate_orig}
\title{Evaluate each variable provided}
\usage{
evaluate_orig(pheno)
}
\arguments{
\item{pheno}{Data.frame with information about the samples}
}
\value{
A matrix with the mean, standard deviation, MAD values of the
numeric variables, the entropy of the categorical, and the amount of
\code{NA} per variable.
}
\description{
Measure some summary statistics of the whole cohort of samples
}
\examples{
data(survey, package = "MASS")
evaluate_orig(survey[, c("Sex", "Age", "Smoke")])
}
\seealso{
Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_index}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_sd}()}
}
\concept{functions to evaluate samples}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indexing.R
\name{use_index}
\alias{use_index}
\title{Convert a factor to index}
\usage{
use_index(x)
}
\arguments{
\item{x}{A character or a factor to be used as index}
}
\description{
Convert a given factor to an accepted index
}
\examples{
plates <- c("P1", "P2", "P1", "P2", "P2", "P3", "P1", "P3", "P1", "P1")
use_index(plates)
}
\seealso{
You can use \code{\link[=evaluate_index]{evaluate_index()}} to evaluate how good an
index is. For the inverse look at  \code{\link[=batch_names]{batch_names()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy.R
\name{entropy}
\alias{entropy}
\title{Calculates the entropy}
\usage{
entropy(x)
}
\arguments{
\item{x}{A character or vector with two or more categories}
}
\value{
The numeric value of the Shannon entropy scaled between 0 and 1.
}
\description{
Calculates the entropy of a category. It uses the amount of categories to
scale between 0 and 1.
}
\note{
It omits the \code{NA} if present.
}
\examples{
entropy(c("H", "T", "H", "T"))
entropy(c("H", "T", "H", "T", "H", "H", "H"))
entropy(c("H", "T", "H", "T", "H", "H", NA))
entropy(c("H", "T", "H", "T", "H", "H"))
entropy(c("H", "H", "H", "H", "H", "H", NA))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporting.R
\name{distribution}
\alias{distribution}
\title{Distribution by batch}
\usage{
distribution(report, column)
}
\arguments{
\item{report}{A data.frame which must contain a batch column. Which can be
obtained with \code{\link[=inspect]{inspect()}}.}

\item{column}{The name of the column one wants to inspect.}
}
\value{
\code{TRUE} if the values are maximal distributed, otherwise \code{FALSE}.
}
\description{
Checks if all the values are maximally distributed in the several batches.
Aimed for categorical variables.
}
\examples{
data(survey, package = "MASS")
columns <- c("Sex", "Age", "Smoke")
index <- design(pheno = survey[, columns], size_subset = 70,
                iterations = 10)
batches <- inspect(index, survey[, columns])
distribution(batches, "Sex")
distribution(batches, "Smoke")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/designer.R
\name{design}
\alias{design}
\title{Design a batch experiment}
\usage{
design(pheno, size_subset, omit = NULL, iterations = 500, name = "SubSet")
}
\arguments{
\item{pheno}{Data.frame with the sample information.}

\item{size_subset}{Numeric value of the number of sample per batch.}

\item{omit}{Name of the columns of the \code{pheno} that will be omitted.}

\item{iterations}{Numeric value of iterations that will be performed.}

\item{name}{A character used to name the subsets, either a single one or a
vector the same size as \code{n}.}
}
\value{
The indices of which samples go with which batch.
}
\description{
Given some samples it distribute them in several batches, trying to have
equal number of samples per batch. It can handle both numeric and
categorical data.
}
\examples{
data(survey, package = "MASS")
index <- design(survey[, c("Sex", "Smoke", "Age")], size_subset = 50,
                iterations = 50)
index
}
\seealso{
The \verb{evaluate_*} functions and \code{\link[=create_subset]{create_subset()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate.R
\name{evaluate_index}
\alias{evaluate_index}
\title{Evaluates a data.frame}
\usage{
evaluate_index(i, pheno)
}
\arguments{
\item{i}{Index}

\item{pheno}{Data.frame with information about the samples}
}
\value{
An array of three dimensions with the mean, standard deviation
(\code{\link[=sd]{sd()}}), and median absolute deviation (\code{\link[=mad]{mad()}}) of the numeric variables, the
entropy of the categorical and the number of \code{NA} by each subgroup.
}
\description{
Measures several indicators per group
}
\examples{
data(survey, package = "MASS")
index <- create_subset(nrow(survey), 50, 5)
ev_index <- evaluate_index(index, survey[, c("Sex", "Smoke")])
ev_index["entropy", , ]
}
\seealso{
If you have already an index you can use \code{\link[=use_index]{use_index()}}.

Other functions to evaluate samples: 
\code{\link{evaluate_entropy}()},
\code{\link{evaluate_independence}()},
\code{\link{evaluate_mad}()},
\code{\link{evaluate_mean}()},
\code{\link{evaluate_na}()},
\code{\link{evaluate_orig}()},
\code{\link{evaluate_sd}()}
}
\concept{functions to evaluate samples}
