
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
