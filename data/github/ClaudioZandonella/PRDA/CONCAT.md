Todo Moving Branch

- [ ] change branch in badges readme
- [ ] change branch in codecov code in travis.yml
- [ ] change branch in the link contributing "https://github.com/ClaudioZandonella/PRDA/blob/develop/CONTRIBUTING.md""
# Todo rOpenSci Packages

- [x] You can choose to use = over <- as long you are consistent with one choice within your package. 
- [x] README should include: https://devguide.ropensci.org/building.html#readme
- [x] We recommend not creating README.md directly, but from a README.Rmd file (an R Markdown file) if you have any demonstration code. (usethis::use_readme_rmd())
- [x] Add #' @noRd to internal functions.
- [x] Only use package startup messages when necessary (function masking for instance). Avoid package startup messages like “This is foobar 2.4-0” or citation guidance because they can be annoying to the user. Rely on documentation for such guidance.
- [x] https://devguide.ropensci.org/building.html#website
- [x] Use Imports instead of Depends for packages providing functions
- [x] use the goodpractice package (goodpractice::gp()) as a guide to improve your package, since most exceptions to it will need to be justified
- [x] Test coverage below 75% will likely require additional tests or explanation before being sent for review. Once you’ve set up CI, use your package’s code coverage report (cf this section of our book) to identify untested lines, and to add further tests. usethis::use_coverage()
- [x] Both test status and code coverage should be reported via badges in your package README.
- [x] R packages should have CI for all platforms when they contain Compiled code
- [x] https://devguide.ropensci.org/ci.html continuos integration
- [x] We urge package maintainers to make sure they are receiving GitHub notifications, as well as making sure emails from rOpenSci staff and CRAN maintainers are not going to their spam box.
- [x] If you would like your package to also be submitted to Journal of Open-Source Software (JOSS), it should include a paper.md file describing the package. More detail on JOSS’s requirements can be found at their website. https://joss.theoj.org/about#author_guidelines
- [x] use repostatus.org badges (which we recommend) https://www.repostatus.org/
- [ ] name function object_action()

send

- [ ] Next, open a new issue in the software review repository and fill out the template. https://github.com/ropensci/software-review/issues/new


other:

- [x] problema del welch's t test e correzione d di cohen
- [x] change name repo and package
- [x] Add more description on hypothetical effect size in vignette overwie and readme
- [ ] Update Zenodo
- [x] Updarte URL in DESCRIPTION
- [ ] Coverage for multiple CI (https://devguide.ropensci.org/ci.html#coverage)
- [ ] Problem with .covrignore function RcppArmadillo
- [x] Fail CI on windows

# Changed name

- [x] url in DESCRIPTION
- [x] url in readme
- [ ] zenodo
- [x] url in paper
- [x] url in vignette PRDA

Joss:

- [ ] https://joss.readthedocs.io/en/latest/submitting.html


After submission:

- [ ] Once you have submitted a package and it has passed editor checks, add a peer-review badge via... (https://devguide.ropensci.org/building.html#readme)
- [ ] After a package is accepted but before transfer, the rOpenSci footer should be added to the bottom of the README file with the following markdown line... (https://devguide.ropensci.org/building.html#readme)
- [ ] After transfer to our GitHub organization, rOpenSci Code of Conduct will apply to your project. Please add this text to the README (https://devguide.ropensci.org/collaboration.html#friendlyfiles)

## Todo list

- [ ] Calcolare Cohen's d for Welch test.
https://www.datanovia.com/en/lessons/t-test-effect-size-using-cohens-d-measure/ riporta:
$$d = \frac{M_x - M_y}{\sqrt{\frac{(Var_1 + Var_2)}{2}}}$$
Non capisco se Var si intende la varianza campionaria o la stima della varianza.
Vedere anche https://stats.stackexchange.com/questions/210352/do-cohens-d-and-hedges-g-apply-to-the-welch-t-test
https://arxiv.org/pdf/1901.09581.pdf

- [ ] Controllare anche come è calcolato il Cohen's d nel caso del paired t-test.
https://www.datanovia.com/en/lessons/t-test-effect-size-using-cohens-d-measure/ riporta:
$$d = \frac{nx-2}{nx-1.25} \times \frac{mx}{sd(x)}$$

- [ ] Come calcolare il valore critico di Cohen's d dato il valre della statistica test.
https://www.bwgriffin.com/gsu/courses/edur9131/content/Effect_Sizes_pdf5.pdf riporta delle formule da controllare.
https://www.frontiersin.org/articles/10.3389/fpsyg.2013.00863/full altre formule (meglio). Ci vanno oppure no le correzioni?
  - One.sampled t-test $d = t / \sqrt{n}$
  - Paired t-test $d = \frac{t}{\sqrt{n}}+ \mu$
  - Two sampled t.test $d =  t \sqrt{\frac{n_x+n_y}{n_x n_y}} + \mu$
  - Welch t-test $d =  t \sqrt{\frac{2}{n_x n_y}\frac{n_x * var_y+n_y*var_x}{var_x + var_y}} + \mu$ (noi facciamo in realtà il smpling con entrampe le pop variance = 1 quindi dovrebbe essere uguale al t-test?)
  
  
- [ ] Troncare le distribuzioni basta eliminare i valori non voluti e continuare il sampling?

- [ ] In prospective quando utilizzo una distribuzione per definire gli effetti considero raggiunto il livello di potenza considerando la media delle potenze. Potrebbe essere meglio considerare la mediana o la media va bene?






<!-- README.md is generated from README.Rmd. Please edit that file -->

# PRDA: Prospective and Retrospective Design Analysis

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
status](https://www.r-pkg.org/badges/version/PRDA)](https://CRAN.R-project.org/package=PRDA)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ClaudioZandonella/PRDA?branch=master&svg=true)](https://ci.appveyor.com/project/ClaudioZandonella/PRDA/branch/master)
[![Travis build
status](https://travis-ci.org/ClaudioZandonella/PRDA.svg?branch=master)](https://travis-ci.org/ClaudioZandonella/PRDA)
[![Codecov test
coverage](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master/graph/badge.svg)](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02810/status.svg)](https://doi.org/10.21105/joss.02810)
[![DOI](https://zenodo.org/badge/212573857.svg)](https://zenodo.org/badge/latestdoi/212573857)

<hr>

<!-- badges: end -->

{PRDA} allows performing a prospective or retrospective design analysis
to evaluate inferential risks (i.e., power, Type M error, and Type S
error) in a study considering Pearson’s correlation between two
variables or mean comparisons (one-sample, paired, two-sample, and
Welch’s *t*-test).

For an introduction to design analysis and a general overview of the
package see `vignette("PRDA")`. Examples for retrospective design
analysis and prospective design analysis are provided in
`vignette("retrospective")` and `vignette("prospective")` respectively.

All the documentation is available at
<https://claudiozandonella.github.io/PRDA/>.

## Installation

You can install the released version of PRDA from
[CRAN](https://CRAN.R-project.org/package=PRDA) with:

``` r
install.packages("PRDA")
```

And the development version from
[GitHub](https://github.com/ClaudioZandonella/PRDA/tree/master) with:

``` r
# install.packages("devtools")
devtools::install_github("ClaudioZandonella/PRDA",
                         build_vignettes = TRUE)
```

## The Package

{PRDA} package can be used for Pearson’s correlation between two
variables or mean comparisons (i.e., one-sample, paired, two-sample, and
Welch’s t-test) considering an hypothetical value of *ρ* or Cohen’s *d*
respectively. See `vignette("retrospective")` and
`vignette("prospective")` to know how to set function arguments for the
different effect types.

### Functions

In {PRDA} there are two main functions `retrospective()` and
`prospective()`.

#### • `retrospective()`

Given the hypothetical population effect size and the study sample size,
the function `retrospective()` performs a retrospective design analysis.
According to the defined alternative hypothesis and the significance
level, the inferential risks (i.e., Power level, Type M error, and Type
S error) are computed together with the critical effect value (i.e., the
minimum absolute effect size value that would result significant).

Consider a study that evaluated the correlation between two variables
with a sample of 30 subjects. Suppose that according to the literature
the hypothesized effect is *ρ* = .25. To evaluate the inferential risks
related to the study we use the function `retrospective()`.

``` r
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 30, 
              test_method = "pearson")
#> 
#>  Design Analysis
#> 
#> Hypothesized effect:  rho = 0.25 
#> 
#> Study characteristics:
#>    test_method   sample_n1   sample_n2   alternative   sig_level   df
#>    pearson       30          NULL        two_sided     0.05        28
#> 
#> Inferential risks:
#>    power   typeM   typeS
#>    0.27    1.826   0.003
#> 
#> Critical value(s): rho  =  ± 0.361
```

In this case, the statistical power is almost 30% and the associated
Type M error and Type S error are respectively around 1.80 and 0.003.
That means, statistical significant results are on average an
overestimation of 80% of the hypothesized population effect and there is
a .3% probability of obtaining a statistically significant result in the
opposite direction.

To know more about function arguments and further examples see the
function documentation `?retrospective` and `vignette("retrospective")`.

#### • `prospective()`

Given the hypothetical population effect size and the required power
level, the function `prospective()` performs a prospective design
analysis. According to the defined alternative hypothesis and the
significance level, the required sample size is computed together with
the associated Type M error, Type S error, and the critical effect value
(i.e., the minimum absolute effect size value that would result
significant).

Consider a study that will evaluate the correlation between two
variables. Knowing from the literature that we expect an effect size of
*ρ* = .25, the function `prospective()` can be used to compute the
required sample size to obtain a power of 80%.

``` r
prospective(effect_size = .25, power = .80, test_method = "pearson",
            display_message = FALSE)
#> 
#>  Design Analysis
#> 
#> Hypothesized effect:  rho = 0.25 
#> 
#> Study characteristics:
#>    test_method   sample_n1   sample_n2   alternative   sig_level   df 
#>    pearson       122         NULL        two_sided     0.05        120
#> 
#> Inferential risks:
#>    power   typeM   typeS
#>    0.797   1.119   0    
#> 
#> Critical value(s): rho  =  ± 0.178
```

The required sample size is \(n=122\), the associated Type M error is
around 1.10 and the Type S error is approximately 0.

To know more about function arguments and further examples see the
function documentation `?prospective` and `vignette("prospective")`.

### Hypothetical effect size

The hypothetical population effect size can be defined as a single value
according to previous results in the literature or experts indications.
Alternatively, {PRDA} allows users to specify a distribution of
plausible values to account for their uncertainty about the hypothetical
population effect size. To know how to specify the hypothetical effect
size according to a distribution and an example of application see
`vignette("retrospective")`.

## Contributing to PRDA

The PRDA package is still in the early stages of its life. Thus, surely
there are many bugs to fix and features to propose. Anyone is welcome to
contribute to the PRDA package.

Please note that this project is released under a [Contributor Code of
Conduct](https://www.contributor-covenant.org/). By contributing to this
project, you agree to abide by its terms.

#### Bugs and New Features

To propose a new feature or to report a bug, please open an issue on
[GitHub](https://github.com/ClaudioZandonella/PRDA/issues). See
[Community
guidelines](https://github.com/ClaudioZandonella/PRDA/blob/master/CONTRIBUTING.md).

#### Future Plans

  - Improve compute time by parallelizing the code
  - Implement design analysis in the case of linear regression models

## Citation

To cite {PRDA} in publications use:

Zandonella Callegher, C., Pastore, M., Andreella, A., Vesely, A.,
Toffalini, E., Bertoldo, G., & Altoè G. (2020). PRDA: Prospective and
Retrospective Design Analysis (Version 1.0.0). Zenodo.
<https://doi.org/10.5281/zenodo.4044214>

A BibTeX entry for LaTeX users is

    @Misc{,
        author       = {Zandonella Callegher, Claudio and Pastore, Massimiliano and Andreella, Angela and 
                        Vesely, Anna and Toffalini, Enrico and Bertoldo, Giulia and Altoè, Gianmarco},
        title        = {PRDA: Prospective and Retrospective Design 
                       Analysis},
        year         = 2020,
        publisher    = {Zenodo},
        version      = {1.0.0},
        doi          = {10.5281/zenodo.4044214},
        url          = {https://doi.org/10.5281/zenodo.4044214}
      }

## References

<div id="refs" class="references">

<div id="ref-altoeEnhancingStatisticalInference2020">

Altoè, Gianmarco, Giulia Bertoldo, Claudio Zandonella Callegher, Enrico
Toffalini, Antonio Calcagnì, Livio Finos, and Massimiliano Pastore.
2020. “Enhancing Statistical Inference in Psychological Research via
Prospective and Retrospective Design Analysis.” *Frontiers in
Psychology* 10. <https://doi.org/10.3389/fpsyg.2019.02893>.

</div>

<div id="ref-bertoldoDesigningStudiesEvaluating2020">

Bertoldo, Giulia, Claudio Zandonella Callegher, and Gianmarco Altoè.
2020. “Designing Studies and Evaluating Research Results: Type M and
Type S Errors for Pearson Correlation Coefficient.” Preprint. PsyArXiv.
<https://doi.org/10.31234/osf.io/q9f86>.

</div>

<div id="ref-gelmanPowerCalculationsAssessing2014">

Gelman, Andrew, and John Carlin. 2014. “Beyond Power Calculations:
Assessing Type S (Sign) and Type M (Magnitude) Errors.” *Perspectives on
Psychological Science* 9 (6): 641–51.
<https://doi.org/10.1177/1745691614551642>.

</div>

</div>

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

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
* Focusing on what is best not just for us as individuals, but for the
overall community

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

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
[INSERT CONTACT METHOD].
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

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

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
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contributing to PRDA

The PRDA package is still in the early stages of its life. Thus, surely there are many bugs to fix and features to propose. Anyone is welcome to contribute to the PRDA package.

Please note that this project is released under a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.

## How to report or fix a Bug

In the case of unexpected behaviours or errors, open an issue. Try to describe with adequate details the problem and, if possible, provide a minimal reproducible example. This will be very useful to solve the problem.

Alternatively, if you want to directly propose a solution, open a pull request. This should concern, however, only small fixes (e.s. documentation or typos). In the case of solutions involving a large amount of changes, it is better to open an issue to discuss and plan them before. 

## How to propose new features

The PRDA package is still in the early stages of its life, so if you want to propose a new feature first open an issue to discuss and plan them. See the Future Plan list below for possible new features.

### Future Plans

- Improve compute time by parallelizing the code
- Implement design analysis in the case of linear regression models
## Resubmission

This is a resubmission. In this version I have:

* Replaced \dontrun{} with \donttest{} in the functions examples that are executable in more than 5 sec.

* Removed the seed argument in the functions to avoid modifying the .GlobalEnv, in line with CRAN policies.


## Test environments
* Mac OS X 10.13 (on travis-ci):
  - release R 4.0.3
  - devel R 4.1.0 (2020-11-27 r79522)
* Ubuntu 16.04 (on travis-ci):
  - release R 4.0.2
  - devel R 4.1.0 (2020-11-26 r79513)
* win-builder:
  - release R 4.0.3
  - devel R 4.1.0 (2020-11-24 r79490)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Claudio Zandonella Callegher <claudiozandonella@gmail.com>'
> 
> New submission

This is my first submission.

> Possibly mis-spelled words in DESCRIPTION:
>  Alto� (43:33)
>  Bertoldo (43:57)
>  Gelman (41:5)
>  al (43:42, 43:69)
>  et (43:39, 43:66)

The encoding is UTF-8. In the first case the surname is Altoè and is correctly displayed in all vignettes and documentation. In the other cases names and words are correctly spelled.


---
title: 'PRDA: An R package for Prospective and Retrospective Design Analysis'
tags:
  - R
  - design analysis
  - power analysis
  - Type M error
  - Type S error
  - replicabiliyt
authors:
  - name: Claudio Zandonella Callegher
    orcid: 0000-0001-7721-6318
    affiliation: 1
  - name: Giulia Bertoldo
    orcid: 0000-0002-6960-3980
    affiliation: 1
  - name: Enrico Toffalini
    orcid: 0000-0002-1404-5133
    affiliation: 3
  - name: Anna Vesely
    orcid: 0000-0001-6696-2390
    affiliation: 2
  - name: Angela Andreella
    orcid: 0000-0002-1141-3041
    affiliation: 2
  - name: Massimiliano Pastore
    orcid: 0000-0002-7922-6365
    affiliation: 1
  - name: Gianmarco Altoè
    orcid: 0000-0003-1154-9528
    affiliation: 1
affiliations:
 - name: Department of Developmental Psychology and Socialisation, University of Padova, Padova, Italy
   index: 1
 - name: Department of Statistical Sciences, University of Padova, Padova, Italy
   index: 2
 - name: Department of General Psychology, University of Padova, Padova, Italy
   index: 3
date: "10 December, 2020"
bibliography: paper_JOSS.bib
editor_options: 
  chunk_output_type: console
---




# Summary

*Design Analysis* was introduced by @gelmanPowerCalculationsAssessing2014 as an extension of Power Analysis. Traditional power analysis has a narrow focus on statistical significance. Design analysis, instead, evaluates together with power levels also other inferential risks (i.e., Type M error and Type S error), to assess estimates uncertainty under hypothetical replications of a study.

Given an hypothetical value of effect size and study characteristics (i.e., sample size, statistical test directionality, significance level),
*Type M error* (Magnitude, also known as Exaggeration Ratio) indicates the factor by which a statistically significant effect is on average exaggerated. *Type S error* (Sign), instead, indicates the probability of finding a statistically significant result in the opposite direction to the hypothetical effect.

Although Type M error and Type S error depend directly on power level, they underline valuable information regarding estimates uncertainty that would otherwise be overlooked. This enhances researchers awareness about the inferential risks related to their studies and helps them in the interpretation of their results. However, design analysis is rarely applied in real research settings also for the lack of dedicated software.

To know more about design analysis consider @gelmanPowerCalculationsAssessing2014 and @luNoteTypeErrors2018. While, for an introduction to design analysis with examples in psychology see @altoeEnhancingStatisticalInference2020 and  @bertoldoDesigningStudiesEvaluating2020.


# Statement of need 

`PRDA` is an R package performing prospective or retrospective design analysis to evaluate inferential risks (i.e., power, Type M error, and Type S error) in a study considering Pearson's correlation between two variables or mean comparisons (one-sample, paired, two-sample, and Welch's *t*-test). *Prospective Design Analysis* is performed in the planning stage of a study to define the required sample size to obtain a given level of power. *Retrospective Design Analysis*, instead, is performed when the data have already been collected to evaluate the inferential risks associated with the study.

Another recent R package, `retrodesign` [@timmRetrodesignToolsType2019], allows conducting retrospective design analysis considering estimate of the unstandardized effect size (i.e., regression coefficient or mean difference) and standard error of the estimate. `PRDA` package, instead, considers standardized effect size (i.e., Pearson correlation coefficient or Cohen's *d*) and study sample size. These are more commonly used in research fields such as Psychology or Social Science, and therefore are implemented in `PRDA` to facilitate researchers' reasoning about design analysis. `PRDA`, additionally, offers the possibility to conduct a prospective design analysis and to account for the uncertainty about the hypothetical value of effect size. In fact, hypothetical effect size can be defined as a single value according to previous results in the literature or experts indications, or by specifying a distribution of plausible values.

The package is available from GitHub (https://github.com/ClaudioZandonella/PRDA) and CRAN (https://CRAN.R-project.org/package=PRDA). Documentation about the package is available at https://claudiozandonella.github.io/PRDA/.

# Examples

Imagine a study evaluating the relation a given personality trait (e.g., introversion) and math performance. Suppose that 20 participants were included in the study and results indicated a statistically significant correlation (e.g, $r = .55, p = .012$). The magnitude of the estimated correlation, however, is beyond what could be considered plausible in this field. 

## Retrospective design analysis

Suppose previous results in the literature indicate correlations in this area are more likely to be around $\rho = .25$. To evaluate the inferential risks associated with the study design, we can use the function `retrospective()`.


```r
library(PRDA)

set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 20, test_method = "pearson")
```

```
## 
## 	Design Analysis
## 
## Hypothesized effect:  rho = 0.25 
## 
## Study characteristics:
##    test_method   sample_n1   sample_n2   alternative   sig_level   df
##    pearson       20          NULL        two_sided     0.05        18
## 
## Inferential risks:
##    power   typeM   typeS
##    0.185   2.161   0.008
## 
## Critical value(s): rho  =  ± 0.444
```

In the output, we have the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. We obtained a statistical power of almost 20% that is associated with a Type M error of around 2.2 and a Type S error of 0.01. That means, statistical significant results are on average an overestimation of 120% of the hypothesized population effect and there is a 1% probability of obtaining a statistically significant result in the opposite direction. To know more about function arguments and examples see the function documentation and vignette.

### Effect size distribution

Alternatively, if no precise information about hypothetical effect size is available, researchers could specify a distribution of values  to account for their uncertainty. For example, they might define a normal distribution with mean of .25 and standard deviation of .1, truncated between .10 and 40.


```r
retrospective(effect_size = function(n) rnorm(n, .25, .1), sample_n1 = 20,
              test_method = "pearson", tl = .1, tu = .4, B = 1e3, 
              display_message = FALSE)
```

```
## Truncation could require long computational time
```

```
## 
## 	Design Analysis
## 
## Hypothesized effect:  rho ~ rnorm(n, 0.25, 0.1) [tl =  0.1 ; tu = 0.4 ]
##    n_effect   Min.    1st Qu.   Median   Mean    3rd Qu.   Max.
##    1000       0.101   0.197     0.25     0.252   0.308     0.4 
## 
## Study characteristics:
##    test_method   sample_n1   sample_n2   alternative   sig_level   df
##    pearson       20          NULL        two_sided     0.05        18
## 
## Inferential risks:
##         Min.    1st Qu.   Median   Mean       3rd Qu.   Max. 
## power   0.055   0.133     0.1880   0.203727   0.26600   0.449
## typeM   1.407   1.785     2.1645   2.347745   2.70075   5.263
## typeS   0.000   0.000     0.0060   0.017573   0.02300   0.246
## 
## Critical value(s): rho  =  ± 0.444
```

Consequently this time we obtained a distribution of values for power, Type M error, and Type S error. Summary information are provided in the output.

## Prospective design analysis

Given the previous results, researchers might consider planning a replication study to obtain more reliable results. The function `prospective()` can be used to compute the sample size needed to obtain a given level of power (e.g., power = 80%).


```r
prospective(effect_size = .25, power = .8, test_method = "pearson",
            display_message = FALSE)
```

```
## 
## 	Design Analysis
## 
## Hypothesized effect:  rho = 0.25 
## 
## Study characteristics:
##    test_method   sample_n1   sample_n2   alternative   sig_level   df 
##    pearson       122         NULL        two_sided     0.05        120
## 
## Inferential risks:
##    power   typeM   typeS
##    0.796   1.12    0    
## 
## Critical value(s): rho  =  ± 0.178
```

In the output, we have again the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. To obtain a power of around 80% the required sample size is $n = 122$, the associated Type M error is around 1.10 and the Type S error is approximately 0. To know more about function arguments and examples see the function documentation and vignette.

In `PRDA` there are no implemented functions to obtain graphical representations of the results. However, it is easy to access all the results and use them to create the plots according to your own needs and preferences. See vignettes for an example.




# References



