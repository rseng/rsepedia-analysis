
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stabm

[![R-CMD-check](https://github.com/bommert/stabm/workflows/R-CMD-check/badge.svg)](https://github.com/bommert/stabm/actions)
[![CRAN
Status](https://www.r-pkg.org/badges/version-ago/stabm)](https://cran.r-project.org/package=stabm)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03010/status.svg)](https://doi.org/10.21105/joss.03010)

`stabm` provides an implementation of many measures which assess the
stability of feature selection. The following stability measures are
currently included:

``` r
stabm::listStabilityMeasures()
#>                           Name Corrected Adjusted Minimum Maximum
#> 1               stabilityDavis     FALSE    FALSE       0       1
#> 2                stabilityDice     FALSE    FALSE       0       1
#> 3             stabilityHamming     FALSE    FALSE       0       1
#> 4   stabilityIntersectionCount      TRUE     TRUE    <NA>       1
#> 5  stabilityIntersectionGreedy      TRUE     TRUE    <NA>       1
#> 6     stabilityIntersectionMBM      TRUE     TRUE    <NA>       1
#> 7    stabilityIntersectionMean      TRUE     TRUE    <NA>       1
#> 8             stabilityJaccard     FALSE    FALSE       0       1
#> 9               stabilityKappa      TRUE    FALSE      -1       1
#> 10         stabilityLustgarten      TRUE    FALSE      -1       1
#> 11           stabilityNogueira      TRUE    FALSE      -1       1
#> 12         stabilityNovovicova     FALSE    FALSE       0       1
#> 13             stabilityOchiai     FALSE    FALSE       0       1
#> 14                stabilityPhi      TRUE    FALSE      -1       1
#> 15           stabilitySechidis     FALSE     TRUE    <NA>      NA
#> 16              stabilitySomol      TRUE    FALSE       0       1
#> 17         stabilityUnadjusted      TRUE    FALSE      -1       1
#> 18               stabilityWald      TRUE    FALSE     1-p       1
#> 19                 stabilityYu      TRUE     TRUE    <NA>       1
#> 20           stabilityZucknick     FALSE     TRUE       0       1
```

## Installation

You can install the released version of stabm from
[CRAN](https://cran.r-project.org/package=stabm) with:

``` r
install.packages("stabm")
```

For the development version, use
[devtools](https://cran.r-project.org/package=devtools):

``` r
devtools::install_github("bommert/stabm")
```

## Contributions

This R package is licensed under the
[LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.en.html). If you
encounter problems using this software (lack of documentation,
misleading or wrong documentation, unexpected behaviour, bugs, …) or
just want to suggest features, please open an issue in the [issue
tracker](https://github.com/bommert/stabm/issues). Pull requests are
welcome and will be included at the discretion of the author.

## Code of Conduct

Please note that the `stabm` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Citation

If you use stabm, please cite our [JOSS
article](https://doi.org/10.21105/joss.03010):

    @Article{stabm,
      title = {{stabm}: Stability Measures for Feature Selection},
      author = {Andrea Bommert and Michel Lang},
      journal = {Journal of Open Source Software},
      year = {2021},
      doi = {10.21105/joss.03010},
      publisher = {The Open Journal},
      volume = {6},
      number = {59},
      pages = {3010},
    }
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
---
title: 'stabm: Stability Measures for Feature Selection'
tags:
  - R
  - feature selection stability
  - stability measures
  - similarity measures
authors:
  - name: Andrea Bommert
    orcid: 0000-0002-1005-9351
    affiliation: 1
  - name: Michel Lang
    orcid: 0000-0001-9754-0393
    affiliation: 1
affiliations:
 - name: Faculty of Statistics, TU Dortmund University, 44221 Dortmund, Germany
   index: 1
date: 05 February 2021
bibliography: paper.bib
---

# Summary
The R [@R] package *stabm* provides functionality for quantifying the similarity of two or more sets.
For example, consider the two sets $\{A, B, C, D\}$ and $\{A, B, C, E\}$.
Intuitively, these sets are quite similar because their overlap is large compared to the cardinality of the two sets.
The R package *stabm* implements functions to express the similarity of sets by a real valued score.
Quantifying the similarity of sets is useful for comparing sets of selected features.
But also for many other tasks like similarity analyses of gene sets or text corpora, the R package *stabm* can be employed.

In the context of feature selection, the similarity of sets of selected features is assessed in order to determine the stability of a feature selection algorithm.
The stability of a feature selection algorithm is defined as the robustness of the set of selected features towards different data sets from the same data generating distribution [@kalousis2007stability].
For stability assessment, either *m* data sets from the same data generating process are available or *m* data sets are created from one data set.
The latter is often achieved with subsampling or random perturbations [@awada2012review].
Then, the feature selection algorithm of interest is applied to each of the *m* data sets, resulting in *m* feature sets.
To quantify the stability of the feature selection algorithm, the similarity of the *m* sets is calculated.
In the context of feature selection stability, set similarity measures are called stability measures.

The R package *stabm* provides an open-source implementation of the 20 stability measures displayed in the table below.
Argument checks are performed with checkmate [@lang2017checkmate] to provide helpful error messages.
It is publicly available on CRAN and on Github and it has only a few dependencies.

|Name | Reference|
|-----|----------|
|stabilityDavis | @davis2006reliable|
|stabilityDice | @dice1945measures|
|stabilityHamming | @dunne2002solutions|
|stabilityIntersectionCount | @bommert2020adjusted|
|stabilityIntersectionGreedy | @bommert2020adjusted|
|stabilityIntersectionMBM | @bommert2020adjusted|
|stabilityIntersectionMean | @bommert2020adjusted|
|stabilityJaccard | @jaccard1901etude|
|stabilityKappa | @carletta1996assessing|
|stabilityLustgarten | @lustgarten2009measuring|
|stabilityNogueira | @nogueira2018stability|
|stabilityNovovicova | @novovicova2009new|
|stabilityOchiai | @ochiai1957zoogeographic|
|stabilityPhi | @nogueira2016measuring|
|stabilitySechidis | @sechidis2020stability|
|stabilitySomol | @somol2008evaluating|
|stabilityUnadjusted | @bommert2020adjusted|
|stabilityWald | @wald2013stability|
|stabilityYu | @yu2012stable|
|stabilityZucknick | @zucknick2008comparing|

# Statement of Need
The R package *stabm* provides an implementation of many stability measures.
For theoretical and empirical comparative studies of the stability measures implemented in *stabm*, we refer to @bommert2017multicriteria, @bommert2020adjusted, @bommert2020integration, and @nogueira2018stability.
It has been demonstrated that considering the feature selection stability when fitting a predictive model often is beneficial for obtaining models with high predictive accuracy [@bommert2017multicriteria; @bommert2020integration; @schirra2016selection].
The stability measures implemented in the R package *stabm* have been employed in @bommert2017multicriteria, @bommert2020benchmark, @bommert2020adjusted, and @bommert2020integration.

# Related Software
A subset of the implemented stability measures is also available in other R or Python packages.
The R package *sets* [@meyer2009sets] and the Python package *scikit-learn* [@pedregosa2011scikit] provide an implementation of the Jaccard index [@jaccard1901etude] to assess the similarity of two sets.
The Python package *GSimPy* [@zhang2020gsimpy] implements the Jaccard index, the Dice index [@dice1945measures], and the Ochiai index [@ochiai1957zoogeographic].
The source code for the publication @nogueira2018stability provides an implementation of their stability measure in R, Python, and Matlab.

# Acknowledgements

This work was supported by the German Research Foundation (DFG), Project RA 870/7-1, and Collaborative Research Center SFB 876, A3.

# References
