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



---
output: github_document
editor_options: 
  chunk_output_type: console
bibliography: vignettes/PRDA.bib
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)

library(PRDA)
```

# PRDA: Prospective and Retrospective Design Analysis

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN status](https://www.r-pkg.org/badges/version/PRDA)](https://CRAN.R-project.org/package=PRDA)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ClaudioZandonella/PRDA?branch=master&svg=true)](https://ci.appveyor.com/project/ClaudioZandonella/PRDA/branch/master)
[![Travis build status](https://travis-ci.org/ClaudioZandonella/PRDA.svg?branch=master)](https://travis-ci.org/ClaudioZandonella/PRDA)
[![Codecov test coverage](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master/graph/badge.svg)](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master)
[![DOI](https://zenodo.org/badge/212573857.svg)](https://zenodo.org/badge/latestdoi/212573857)
<hr>
<!-- badges: end -->

{PRDA} allows performing a prospective or retrospective design analysis to evaluate inferential risks (i.e., power, Type M error, and Type S error) in a study considering Pearson's correlation between two variables or mean comparisons (one-sample, paired, two-sample, and Welch's *t*-test). 

For an introduction to design analysis and a general overview of the package see `vignette("PRDA")`.
Examples for retrospective design analysis and prospective design analysis are provided in `vignette("retrospective")` and `vignette("prospective")` respectively.

All the documentation is available at https://claudiozandonella.github.io/PRDA/.

## Installation

You can install the released version of PRDA from [CRAN](https://CRAN.R-project.org/package=PRDA) with:

``` r
install.packages("PRDA")
```

And the development version from [GitHub](https://github.com/ClaudioZandonella/PRDA/tree/master) with:

``` r
# install.packages("devtools")
devtools::install_github("ClaudioZandonella/PRDA",
                         build_vignettes = TRUE)
```


## The Package

{PRDA} package can be used for Pearson's correlation between two variables or mean comparisons (i.e., one-sample, paired, two-sample, and Welch's t-test) considering an hypothetical value of *&rho;* or Cohen's *d* respectively. See `vignette("retrospective")` and `vignette("prospective")` to know how to set function arguments for the different effect types. 


### Functions

In {PRDA} there are two main functions `retrospective()` and `prospective()`.

#### &#8226; `retrospective()`

Given the hypothetical population effect size and the study sample size, the function `retrospective()` performs a retrospective design analysis. According to the defined alternative hypothesis and the significance level, the inferential risks (i.e., Power level, Type M error, and Type S error) are computed together with the critical effect value (i.e., the minimum absolute effect size value that would result significant).

Consider a study that evaluated the correlation between two variables with a sample of 30 subjects. Suppose that according to the literature the hypothesized effect is *&rho;* = .25. To evaluate the inferential risks related to the study we use the function `retrospective()`.

```{r retrospective,}
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 30, 
              test_method = "pearson")
```

In this case, the statistical power is almost 30% and the associated Type M error and Type S error are respectively around 1.80 and 0.003. That means, statistical significant results are on average an overestimation of 80% of the hypothesized population effect and there is a .3% probability of obtaining a statistically significant result in the opposite direction.

To know more about function arguments and further examples see the function documentation `?retrospective` and  `vignette("retrospective")`.

#### &#8226; `prospective()`

Given the hypothetical population effect size and the required power level, the function `prospective()` performs a prospective design analysis. According to the defined alternative hypothesis and the significance level, the required sample size is computed together with the associated Type M error, Type S error, and the critical effect value (i.e., the minimum absolute effect size value that would result significant).  

Consider a study that will evaluate the correlation between two variables. Knowing from the literature that we expect an effect size of *&rho;* = .25, the function `prospective()` can be used to compute the required sample size to obtain a power of 80%. 
```{r prospective}
prospective(effect_size = .25, power = .80, test_method = "pearson",
            display_message = FALSE)
```

The required sample size is $n=122$, the associated Type M error is around 1.10 and the Type S error is approximately 0. 

To know more about function arguments and further examples see the function documentation `?prospective` and `vignette("prospective")`.

### Hypothetical effect size

The hypothetical population effect size can be defined as a single value according to previous results in the literature or experts indications. Alternatively, {PRDA} allows users to specify a distribution of plausible values to account for their uncertainty about the hypothetical population effect size.  To know how to specify the hypothetical effect size according to a distribution and an example of application see `vignette("retrospective")`.

## Contributing to PRDA

The PRDA package is still in the early stages of its life. Thus, surely there are many bugs to fix and features to propose. Anyone is welcome to contribute to the PRDA package.

Please note that this project is released under a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.

#### Bugs and New Features

To propose a new feature or to report a bug, please open an issue on [GitHub](https://github.com/ClaudioZandonella/PRDA/issues). See [Community guidelines](https://github.com/ClaudioZandonella/PRDA/blob/master/CONTRIBUTING.md).

#### Future Plans

- Improve compute time by parallelizing the code
- Implement design analysis in the case of linear regression models

## Citation

To cite {PRDA} in publications use:

Zandonella Callegher, C., Pastore, M., Andreella, A., Vesely, A., Toffalini, E., Bertoldo, G., & Altoè G. (2020). PRDA: Prospective and Retrospective Design Analysis (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.4044214

A BibTeX entry for LaTeX users is
```{}
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
```

---
nocite: | 
  @altoeEnhancingStatisticalInference2020, @bertoldoDesigningStudiesEvaluating2020, @gelmanPowerCalculationsAssessing2014
...

## References

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
date: "`r format(Sys.time(), '%d %B, %Y')`"
bibliography: paper_JOSS.bib
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


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

```{r retro}
library(PRDA)

set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 20, test_method = "pearson")
```

In the output, we have the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. We obtained a statistical power of almost 20% that is associated with a Type M error of around 2.2 and a Type S error of 0.01. That means, statistical significant results are on average an overestimation of 120% of the hypothesized population effect and there is a 1% probability of obtaining a statistically significant result in the opposite direction. To know more about function arguments and examples see the function documentation and vignette.

### Effect size distribution

Alternatively, if no precise information about hypothetical effect size is available, researchers could specify a distribution of values  to account for their uncertainty. For example, they might define a normal distribution with mean of .25 and standard deviation of .1, truncated between .10 and 40.

```{r retro_dist}
retrospective(effect_size = function(n) rnorm(n, .25, .1), sample_n1 = 20,
              test_method = "pearson", tl = .1, tu = .4, B = 1e3, 
              display_message = FALSE)
```

Consequently this time we obtained a distribution of values for power, Type M error, and Type S error. Summary information are provided in the output.

## Prospective design analysis

Given the previous results, researchers might consider planning a replication study to obtain more reliable results. The function `prospective()` can be used to compute the sample size needed to obtain a given level of power (e.g., power = 80%).

```{r pro}
prospective(effect_size = .25, power = .8, test_method = "pearson",
            display_message = FALSE)
```

In the output, we have again the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. To obtain a power of around 80% the required sample size is $n = 122$, the associated Type M error is around 1.10 and the Type S error is approximately 0. To know more about function arguments and examples see the function documentation and vignette.

In `PRDA` there are no implemented functions to obtain graphical representations of the results. However, it is easy to access all the results and use them to create the plots according to your own needs and preferences. See vignettes for an example.




# References



---
title: "PRDA: Prospective and Retrospective Design Analysis"
description: "This vignette offers a general overview of the PRDA package, introducing the concept of design analysis and providing information about the main functions in the package."
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PRDA}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
bibliography: PRDA.bib
editor_options: 
  chunk_output_type: console
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(PRDA)
```


{PRDA} allows performing a prospective or retrospective design analysis to evaluate inferential risks (i.e., power, Type M error, and Type S error) in a study considering Pearson's correlation between two variables or mean comparisons (one-sample, paired, two-sample, and Welch's *t*-test). 

## Introduction to design analysis

The term *Design Analysis* was introduced by @gelmanPowerCalculationsAssessing2014 as a broader definition of Power Analysis. Traditional power analysis has a narrow focus on statistical significance. Design analysis, instead, evaluates together with power levels also other inferential risks (i.e., Type M error and Type S error), to assess estimates uncertainty under hypothetical replications of a study.

Given an *hypothetical effect size* and  the study characteristics (i.e., sample size, statistical test directionality, $\alpha$ level), design analysis evaluates:

- **Power**: the probability of the test rejecting the null hypothesis.
- **Type M error** (Magnitude): the factor by which a statistically significant effect is on average exaggerated, also known as *Exaggeration Ratio*.
- **Type S error** (Sign): the probability of finding a statistically significant result in the opposite direction to the hypothetical effect.

Moreover, @gelmanPowerCalculationsAssessing2014 distinguished between two types of design analysis according to the study phase:

- **Prospective Design Analysis**: if the analysis is performed in the planning stage of a study to define the sample size needed to obtain a required level of power.
- **Retrospective Design Analysis**: if the analysis is performed in a later stage when the data have already been collected. This is still useful to evaluate the inferential risks associated with the study.

It is important to do not mistake a retrospective design analysis for post-hoc power analysis. The former defines the hypothetical effect size according to previous results in the literature or experts indications, whereas the latter defines the hypothetical effect size based on the same study results and it is a widely-deprecated practice [@goodmanUsePredictedConfidence1994; @lenthStatisticalPowerCalculations12007; @Senn1304].

### Enhancing researchers awareness

Although Type M error and Type S error depend directly on power level, they underline valuable information regarding estimates uncertainty that would otherwise be overlooked. This enhances researchers awareness about the inferential risks related to their studies and helps them in the interpretation of the results.

Despite the lower chances, a statistically significant result could be obtained even in an underpowered study (e.g., power = 20%). This might seem a promising finding, and researchers might think that getting a statistically significant result in an underpowered study means the results must be reliable. Therefore, they would probably be even more confident in the interpretation of their results.

However, in this scenario statistically significant results are almost certain to be an overestimation of the population effect. As pointed out by @gelman_hill_vehtari_2020 "a key risk for a low-power study is not so much that it has a small chance of succeeding, but rather that an apparent success merely masks a larger failure" (p.292). This is also referred as the "Winner's curse", indicating that the apparent win in terms of a statistically significant result is an actual loss as the obtained estimate is inflated.

For example, in a study considering a two-sample *t*-test with 30 participants per group, if the hypothetical population effect size is small (e.g., Cohen's *d* of .25) the actual power is only 16%.  The associated Type M error is around 2.60 and  the Type S error is 0.01. That means, statistical significant results are on average an overestimation of 160% of the hypothesized population effect and there is a 1% probability of obtaining a statistically significant result in the opposite direction. 

In this scenario, knowing the type M and S errors, researchers would be much more cautious in interpreting the results and might consider carrying out a replication study to obtain more reliable results.

### More on design analysis

To know more about design analysis consider @gelmanPowerCalculationsAssessing2014. While, for an introduction to design analysis considering examples in psychology see @altoeEnhancingStatisticalInference2020 and  @bertoldoDesigningStudiesEvaluating2020.

## The package

Given a plausible value of effect size, {PRDA} performs a prospective or retrospective design analysis to evaluate the inferential risks (i.e., power, Type M error, and Type S error) related to the study design.

{PRDA} package can be used for Pearson's correlation between two variables or mean comparisons (i.e., one-sample, paired, two-sample, and Welch's t-test) considering an hypothetical value of $\rho$ or Cohen's *d* respectively. See [`vignette("retrospective")`](retrospective.html) and [`vignette("prospective")`](prospective.html) to know how to set function arguments for the different effect types. 

### Install

You can install the released version of PRDA from [CRAN](https://CRAN.R-project.org/package=PRDA) with:

``` r
install.packages("PRDA")
```

And the development version from [GitHub](https://github.com/ClaudioZandonella/PRDA/tree/master) with:

```{r setup, eval = F}
# If devtools is not installed yet: 
# install.packages( "devtools" )  
devtools::install_github("CaludioZandonella/PRDA",
                         build_vignettes = TRUE)
library(PRDA)
```

### Functions

In {PRDA} there are two main functions:

- **`retrospective()`**.
Given the hypothetical population effect size and the study sample size, the function `retrospective()` performs a retrospective design analysis. According to the defined alternative hypothesis and the significance level, the inferential risks (i.e., Power level, Type M error, and Type S error) are computed together with the critical effect value (i.e., the minimum absolute effect size value that would result significant). To know more about function arguments and examples see the function documentation `?retrospective` and  [`vignette("retrospective")`](retrospective.html).

- **`prospective()`**.
Given the hypothetical population effect size and the required power level, the function `prospective()` performs a prospective design analysis. According to the defined alternative hypothesis and the significance level, the required sample size is computed together with the associated Type M error, Type S error, and the critical effect value (i.e., the minimum absolute effect size value that would result significant).  To know more about function arguments and examples see the function documentation `?prospective` and [`vignette("prospective")`](prospective.html).

### Hypothetical effect size

The hypothetical population effect size can be defined as a single value according to previous results in the literature or experts indications. Alternatively, {PRDA} allows users to specify a distribution of plausible values to account for their uncertainty about the hypothetical population effect size.  To know how to specify the hypothetical effect size according to a distribution and an example of application see [`vignette("retrospective")`](retrospective.html).



## Case study

@eisenbergerDoesRejectionHurt2003 claimed that social and physical pain seem to share similar neural underpinnings. Their experiment included 13 participants, and they found a statistically significant correlation between perceived distress due to social exclusion and activity in the brain area associated with physical pain. However, the magnitude of the estimated correlation ($r = .88$) is beyond what could be considered plausible. In this field correlations are likely to be around $\rho = .25$ [for a complete discussion see @bertoldoDesigningStudiesEvaluating2020].

### Retrospective design analysis

The function `retrospective()` can be used to evaluate the inferential risks associated with the study.

```{r, retro}
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 13, test_method = "pearson")
```

In the output, we have the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. We obtained a statistical power of almost 13% that is associated with a Type M error of around 2.6 and a Type S error of 0.03. That means, statistical significant results are on average an overestimation of 160% of the hypothesized population effect and there is a 3% probability of obtaining a statistically significant result in the opposite direction.

To know more about function arguments and examples see the function documentation `?retrospective` and  [`vignette("retrospective")`](retrospective.html).

### Prospective design analysis

Considering the previous results, researchers might consider planning a replication study to obtain more reliable results. The function `prospective()` can be used to compute the sample size needed to obtain a required level of power (e.g., power = 80%).

```{r, pro}
prospective(effect_size = .25, power = .8, test_method = "pearson", 
            display_message = FALSE)
```

In the output, we have again the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. To obtain a power of around 80% the required sample size is $n = 126$, the associated Type M error is around 1.10 and the Type S error is approximately 0.

To know more about function arguments and examples see the function documentation `?prospective` and [`vignette("prospective")`](prospective.html).

## References



---
title: "Prospective Design Analysis"
description: "This vignette shows how to conduct a prospective design analysis using the function `prospective()`."
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{prospective}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
bibliography: PRDA.bib
editor_options: 
  chunk_output_type: console
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(PRDA)
```

Given the hypothetical population effect size and the required power level, the function `prospective()` performs a prospective design analysis. Prospective design analysis allows to define the sample size needed to obtain a required level of power computing also the associated inferential risks. Function arguments are:

```{r, eval=FALSE, echo = T}
retrospective(effect_size, power, ratio_n = 1,
              test_method = c("pearson", "two_sample", "welch",
                                        "paired", "one_sample"),
              alternative = c("two_sided","less","greater"),
              sig_level = .05, ratio_sd = 1, B = 1e4,
              tl = -Inf, tu = Inf, B_effect = 1e3,
              sample_range = c(2, 1000), tol = .01,
              display_message = TRUE)
```

Complete arguments description is provided in the function documentation `?prospective`. In the following sections, instead, different examples are presented. For further details about design analysis see @altoeEnhancingStatisticalInference2020 and @bertoldoDesigningStudiesEvaluating2020.


## Prospective Design Analysis for Correlation

To conduct a prospective design analysis considering a correlation between two variables, we need to specify `test_method = "pearson"` (default option). Note that, only Pearson's correlation is available, while the Kendall's $\tau$ and Spearman's $\rho$ are not implemented.

### Example 1: Pearson's correlation

Consider a study that will evaluate the correlation between two variables. Knowing from the literature that we expect an effect size of $\rho = .25$, which is the required sample size to obtain a power of 60%? We can use the function `prospective()` setting the argument `test_method = "pearson"`. 

```{r, example1}
set.seed(2020) # set seed to make results reproducible

prospective(effect_size = .25, power = .60, test_method = "pearson",
            display_message = TRUE)
```

The default option `display_message = TRUE` prints the different steps to find the required sample size and the progress bar. Note that, however, the progress bar is available only when `effect_size` is defined as a function. In the output, we have the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. To obtain a power of around 60% the required sample size is $n = 76$, the associated Type M error is almost 1.30 and the Type S error is approximately 0. Finally, the critical values (i.e., the minimum absolute effect size value that would result significant) are $\rho = \pm.226$. Note that correlation tests were conducted considering `"two_sided"` alternative hypothesis and a significance level of .05 (the default settings).


## Prospective Design Analysis for Means Comparison

To conduct a retrospective design analysis considering means comparisons, we need to specify the appropriate *t*-test (i.e., One-sample, Paired, Two-sample, or Welch's *t*-test) using the argument `test_method`. Arguments specifications for the different *t*-tests are presented in the following Table.


| Test                | `test_method` |     Other required arguments    |
| -------------------:|:-------------:|:-------------------------------:|
| One-sample *t*-test | `one_sample`  | `ratio_n = NULL`                |
| Paired *t*-test     | `paired`      | `ratio_n = 1`                   |
| Two-sample *t*-test | `two_sample`  | `ratio_n`                       |
| Welch's *t*-test    | `welch`       | `ratio_n` and `ratio_sd`        |

### Example 2: Paired *t*-test

Imagine that we are planning a study where the same group is measured twice (e.g., pre- and post-test). Knowing from the literature that we expect an effect size of $d = .35$, which is the required sample size to obtain a power of 80%? We can use the function `prospective()` specifying the corresponding arguments. We use the option `test_method = one_sample` for paired *t*-test  and set `ratio_n = 1`.

```{r, example2}
prospective(effect_size = .35, power = .8, test_method = "paired",
            ratio_n = 1, display_message = FALSE)
```

To obtain a power of 80%, the required sample size is $n=66$, the associated Type M error is around 1.10 and the Type S error is approximately 0. Finally, the critical values (i.e., the minimum absolute effect size value that would result significant) are $d  =  \pm 0.246$.

### Example 3: Two-sample *t*-test

Imagine now the case where two groups (e.g., treatment and control group) will be compared. However, we know in advance that the sample size in the two groups will be different (e.g., the number of participants in the treatment group could be limited due to strict selecting criteria). We can define the ratio between the sample size in the first group and in the second group using the `ratio_n` argument. Again, we hypothesize an effect size of $d = .35$, but this time we specify a one-sided alternative hypothesis and a significance level of .10. We can do that using respectively the arguments `alternative` and `sig_level`.

```{r, example3}
prospective(effect_size = .35, power = .80, ratio_n = .5, 
            test_method = "two_sample", alternative = "great", sig_level = .10, 
            display_message = FALSE)
```

The option `test_method = "two_sample"` is used to consider a two-sample *t*-test. To obtain a power of 80%, we would need at least 55 participants in the first group and 110 participants in the second group. The associated Type M error is almost 1.20 and the Type S error is approximately 0. Finally, the critical value is $d = .213$.


### Example 4: Welch's *t*-test

Consider again the previous example, but this time we do not assume homogeneity of variance between the two groups. We suppose, instead, that the ratio between the standard deviation of the first group and of the second group is 1.5. In this case the appropriate test is the Welch's *t*-test. We set the option `test_method = "welch"` and specify the argument `ratio_sd`.


```{r, example4}
prospective(effect_size = .35, power = .80, ratio_n = .5, test_method = "welch",
            ratio_sd = 1.5, alternative = "great", sig_level = .10, 
            display_message = FALSE)
```

Now, to obtain a power of 80%, we would need at least 63 participants in the first group and 126 participants in the second group. The associated Type M error is almost 1.20 and the Type S error is approximately 0. Finally, the critical value is $d = .212$. Results are really close to the previous ones.


## Population effect size distribution

Defining the hypothetical population effect size as a single value could be limiting.  Instead, researchers may prefer to use a probability distribution representing their uncertainty regarding the hypothetical population effect. Note that this could be interpreted as a prior distribution of the population effect in a Bayesian framework. 

To define the hypothetical population effect size (`effect_size`) according to a probability distribution, it is necessary to specify a function that allows sampling values from a given distribution. The function has to be defined as `function(n) my_function(n, ...)`, with only one single argument `n` representing the number of samples (e.g., `function(n) rnorm(n, mean = 0, sd = 1)`). See [`vignette("retrospective")`](retrospective.html) for further details.

### Example 5: Effect size distribution

Consider the same scenario as in the correlation example (Example 1). This time we define the hypothesized effect size according to a normal distribution with mean .30 and standard deviation .10. Moreover, to avoid unreasonable values we truncate the distribution between .15 and .45.

```{r, example5}
prospective(effect_size = function(n) rnorm(n, .3, .1), power = .60, 
            test_method = "pearson", tl = .15, tu = .45, B_effect = 500, 
            B = 500, display_message = FALSE)
```

Note that we adjusted `B_effect` and `B` to find a good trade-off between computational times and results accuracy. Differently from previous outputs, we have now a summary for the sampled effects distribution and for the inferential risks.

## Graphical representation

Currently there are no personalized plot functions in {PRDA}. However, it is easy to access all the results and use them to create the plots according to your needs. 

The function `prospective()` returns a list with class `"design_analysis"` that contains:

- `design_analysis` - a character string indicating the type of design analysis (prospective or retrospective).
- `call_arguments` - a list with all the arguments passed to the function. 
- `effect_info` - a list with all the information regarding the considered hypothetical population effect size. In particular, in `effect_samples` we find the vector with the sampled effects (or unique value in the case of a single value).
- `test_info` - a list with all the information regarding the test performed.
- `prospective_res` - a data frame with the results of the design analysis (i.e., `power`, `typeM`, and `typeS`).

Output complete description is provided in the function help page `?prospective`.

```{r, data_plot}
da_fit <- prospective(effect_size = function(n) rnorm(n, .3, .1), power = .60,
                      test_method = "pearson", tl = .15, tu = .45, 
                      B_effect = 500, B = 500, display_message = FALSE)

str(da_fit, max.level = 1)
```


Similarly to the examples provided in [`vignette("retrospective")`](retrospective.html), results can be used to create the plots according to your needs. See  [`vignette("retrospective")`](retrospective.html) for further details.

---
nocite: | 
  @gelmanPowerCalculationsAssessing2014
...

### References
---
title: "Retrospective Design Analysis"
description: "This vignette shows how to conduct a retrospective design analysis using the function `retrospective()`."
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{retrospective}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
bibliography: PRDA.bib
editor_options: 
  chunk_output_type: console
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
```

```{r, message=FALSE}
library(tidyverse)
library(ggplot2)
library(PRDA)
```

Given the hypothetical population effect size and the study sample size, the function `retrospective()` performs a retrospective design analysis. Retrospective design analysis allows to evaluate the inferential risks associated with a study design when the data have already been collected.  Function arguments are:

```{r, eval=FALSE, echo = T}
retrospective(effect_size, sample_n1, sample_n2 = NULL,
              test_method = c("pearson", "two_sample", "welch",
                              "paired", "one_sample")
              alternative = c("two_sided","less","greater"),
              sig_level = .05, ratio_sd = 1, B = 1e4, 
              tl = -Inf, tu = Inf, B_effect = 1e3,
              display_message = TRUE)
```

Complete arguments description is provided in the function documentation `?retrospective`. In the following sections, instead, different examples are presented. For further details about design analysis see @altoeEnhancingStatisticalInference2020 and @bertoldoDesigningStudiesEvaluating2020.


## Retrospective Design Analysis for Correlation

To conduct a retrospective design analysis considering a correlation between two variables, we need to specify `test_method = "pearson"` (default options). Note that, only Pearson's correlation is available, while the Kendall's $\tau$ and Spearman's $\rho$ are not implemented.

### Example 1: Pearson's correlation

Consider a study that evaluates the correlation between two variables with a sample of 30 subjects. Suppose that according to the literature the hypothesized effect is $\rho = .25$. We can evaluate the inferential risks related to the present study design setting the argument `test_method = "pearson"`.

```{r, example1}
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 30, test_method = "pearson")
```

In the output, we have the summary information about the hypothesized population effect, the study characteristics, and the inferential risks. We obtained a statistical power of almost 30% that is associated with a Type M error of around 1.80 and a Type S error of 0.003. That means, statistical significant results are on average an overestimation of 80% of the hypothesized population effect and there is a .3% probability of obtaining a statistically significant result in the opposite direction. Finally, the critical values (i.e., the minimum absolute effect size value that would result significant) are $\rho = \pm.36$. Note that it is not necessary to specify `sample_2` argument in the case of correlation and tests were conducted considering `"two_sided"` alternative hypothesis with a significance level of .05 (the default settings).


## Retrospective Design Analysis for Means Comparison

To conduct a retrospective design analysis considering means comparisons, we need to specify the appropriate *t*-test (i.e., One-sample, Paired, Two-sample, or Welch's *t*-test) using the argument `test_method`. Arguments specifications for the different *t*-tests are presented in the following Table.


| Test                | `test_method` |     Other required arguments    |
| -------------------:|:-------------:|:-------------------------------:|
| One-sample *t*-test | `one_sample`  | `sample_n2 = NULL`              |
| Paired *t*-test     | `paired`      | `sample_n2` equal to `sample_n1`|
| Two-sample *t*-test | `two_sample`  | `sample_n2`                     |
| Welch's *t*-test    | `welch`       | `sample_n2` and `ratio_sd`      |

### Example 2: Paired *t*-test

Consider a study where the same group ($n = 25$) was measured twice (e.g., pre- and post-test). Knowing from the literature that we would expect an effect size of $d = .35$, which are the inferential risks related to the present study design? We can use the function `retrospective()` specifying the corresponding arguments. We use the option `test_method = one_sample` for paired *t*-test  and set the same value for `sample_n1` and `sample_n2`.

```{r, example2}
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 25,
              test_method = "paired")
```

In this case, we obtained a statistical power of almost 40% that is associated with a Type M error of around 1.60 and the Type S error is approximately 0. That means, statistical significant results are on average an overestimation of 60% of the hypothesized population effect. Finally, the critical values (i.e., the minimum absolute effect size value that would result significant) are $d  =  \pm 0.413$.

### Example 3: Two-sample *t*-test

Consider a case where two groups (e.g., treatment and control group) with respectively 25 and 35 subjects were compared. Again, we hypothesize an effect size of $d = .35$, but this time we specify a one-sided alternative hypothesis and a significance level of .10. We can do that using respectively the arguments `alternative` and `sig_level`.

```{r, example3}
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 35,
              test_method = "two_sample", alternative = "great", 
              sig_level = .10, B = 1e5)
```

The option `test_method = "two_sample"` is used to consider a two-sample *t*-test and `B = 1e5` increases the number of replication for more accurate results. We obtained a statistical power of around 50% that is associated with a Type M error of almost 1.60. Note that in the case of one-sided tests, the type S error will be always 0 or 1 depending on whether the hypothesized effect is coherent with the alternative hypothesis. Finally, the critical value is $d = .34$.

### Example 4: Welch's *t*-test

Consider again the previous example, but this time we do not assume homogeneity of variance between the two groups. We suppose, instead, that the ratio between the standard deviation of the first group and of the second group is 1.5.  In this case the appropriate test is the Welch's *t*-test. We set the option `test_method = "welch"` and specify the argument `ratio_sd`.

```{r, example4}
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 35,
              test_method = "welch", ratio_sd = 1.5, alternative = "great", 
              sig_level = .10, B = 1e5)
```

We obtained a statistical power of 50%, which is associated with a Type M error of almost 1.65, and the critical value is $d = .35$. Results are really close to the previous ones with a slightly higher Type M error.


## Population effect size distribution

Defining the hypothetical population effect size as a single value could be limiting. Instead, researchers may prefer to use a probability distribution representing their uncertainty regarding the hypothetical population effect. Note that this could be interpreted as a prior distribution of the population effect in a Bayesian framework. 

To define the hypothetical population effect size (`effect_size`) according to a probability distribution, it is necessary to specify a function that allows sampling values from a given distribution. The function has to be defined as `function(n) my_function(n, ...)`, with only one single argument `n` representing the number of samples. For example, `function(n) rnorm(n, mean = 0, sd = 1)` would allow to sample from a normal distribution with mean 0 and standard deviation 1; or `function(n) sample(c(.1,.3,.5), n, replace = TRUE)` would allow to sample form a set of three equally plausible values. This allows users to define hypothetical effect size distribution according to their needs.

Argument `B_effect` defines the number of sampled effects. Increase the number to obtain more accurate results, although this will require more computational time (default is `B_effect = 1000`). To avoid long computational times when using a function to define the hypothetical population effect size, we suggest adjusting `B` (i.e., the number of simulation per each effect).

Optional arguments `tl` and `tu` allow truncating the sampling distribution defining the lower truncation point and upper truncation point respectively. Specifying truncation points is recommended as it allows avoiding unreasonable results in the case of unbounded distributions (i.e., too large effects, effects close to zero, or effects in the opposite direction from the expected ones). Note that if `effect_type = "correlation"`, distribution is automatically truncated between -1 and 1.


### Example 5: Effect size distribution

Consider the same scenario as in the correlation example (Example 1). This time we define the hypothesized effect size according to a normal distribution with mean .30 and standard deviation .10. Moreover, to avoid unreasonable values we truncate the distribution between .15 and .45. The argument `display_message = TRUE` (default) allows to print the progress bar. Note that the progress bar is available only when `effect_size` is defined as a function.

```{r, example5}
retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 30,
              test_method = "pearson", tl = .15, tu = .45, B_effect = 1e3, 
              B = 1e3, display_message = TRUE)
```

We adjusted `B_effect` and `B` to find a good trade-off between computational times and results accuracy. Differently from previous outputs, we have now a summary for the sampled effects distribution and for the inferential risks.

## Graphical representation

Currently there are no personalized plot functions in {PRDA}. However, it is easy to access all the results and use them to create the plots according to your needs. Here we present an example using {tidyverse} and {ggplot2}. 

The function `retrospective()` returns a list with class `"design_analysis"` that contains:

- `design_analysis` - a character string indicating the type of design analysis (prospective or retrospective).
- `call_arguments` - a list with all the arguments passed to the function. 
- `effect_info` - a list with all the information regarding the considered hypothetical population effect size. In particular, in `effect_samples` we find the vector with the sampled effects (or unique value in the case of a single value).
- `test_info` - a list with all the information regarding the test performed.
- `retrospective_res` - a data frame with the results of the design analysis (i.e., `power`, `typeM`, and `typeS`).

Output complete description is provided in the function help page `?retrospective`.

```{r, data_plot}
da_fit <- retrospective(effect_size = function(n) rnorm(n, .3, .1), 
                        sample_n1 = 30, test_method = "pearson",
                        tl = .15, tu = .45, B_effect = 1e3, B = 1e3, 
                        display_message = FALSE)

str(da_fit, max.level = 1)
```

Note that the inferential risks associated to the *i* value of the vector `effect_samples` are reported in the *i* row of `retrospective_res` dataframe. Thus, we can simply add `effect_samples` as a new column of the dataframe.

```{r}
data_plot <- da_fit$retrospective_res %>%
  mutate(effect = da_fit$effect_info$effect_samples)
```

Plotting the distribution of sampled effects, we can evaluate whether they accurately represent the intended distribution. If not, we could increase the number of sampled effects (`B_effect`).

```{r, fig.dim= c(4, 3), dev='png'}
ggplot(data_plot)+
  geom_histogram(aes(effect, y = ..density..),
                 col = "black", fill = "#00BFC4", alpha = .8,
                 breaks=seq(.15,.45,.02))+
  scale_x_continuous(breaks = seq(.1,.5,.05), limits = c(.1,.5))+
  labs(x = "Sampled Effects",
       y = "Density")+
  theme_bw()
```

We can also plot the distributions of Power, Type M error, and Type S error.

```{r, fig.dim = c(7.23, 2.5), dev='png'}
data_plot %>%
  pivot_longer(cols = c("power", "typeM", "typeS"), 
               names_to = "Criteria", values_to = "Value") %>%
  mutate(Criteria = recode(Criteria, power = "Power", typeM = "Type M",  typeS = "Type S")) %>%
  ggplot(aes(x = Value, y = ..density.., fill = Criteria)) +
  geom_histogram(col = "black", alpha = .7, bins = 15) + 
  facet_wrap(.~ Criteria, scales = "free") +
  labs(y = "Density") +
  theme_bw() +
  theme(legend.position = "none") 
 
  
```

---
nocite: | 
  @gelmanPowerCalculationsAssessing2014
...

### References

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retrospective.R
\name{retrospective}
\alias{retrospective}
\title{Retrospective Design Analysis}
\usage{
retrospective(
  effect_size,
  sample_n1,
  sample_n2 = NULL,
  test_method = c("pearson", "two_sample", "welch", "paired", "one_sample"),
  alternative = c("two_sided", "less", "greater"),
  sig_level = 0.05,
  ratio_sd = 1,
  B = 10000,
  tl = -Inf,
  tu = Inf,
  B_effect = 1000,
  display_message = TRUE
)
}
\arguments{
\item{effect_size}{a numeric value or function (see Details) indicating the
hypothetical population effect size.}

\item{sample_n1}{a numeric value indicating the sample size of the first
group.}

\item{sample_n2}{a numeric value indicating the sample size of the second
group. This argument is required when \code{test_method} is set to
\code{"two_sample"} or \code{"welch"}. In the case of \code{test_method =
"paired"}, set \code{sample_n2} equal to \code{sample_n1}. Whereas in the
case of \code{test_method = "one_sample"}, set \code{sample_n2} to
\code{NULL}. This argument is ignored for \code{test_method = "pearson"}.
See Test methods section in Details.}

\item{test_method}{a character string specifying the test type, must be one of
\code{"pearson"} (default, Pearson's correlation), \code{"two_sample"}
(independent two-sample \emph{t}-test), \code{"welch"} (Welch's
\emph{t}-test), \code{"paired"} (dependent \emph{t}-test for paired
samples), or \code{"one_sample"} (one-sample \emph{t}-test). You can specify
just the initial letters.}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of \code{"two_sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter.}

\item{sig_level}{a numeric value indicating the significance level on which
the alternative hypothesis is evaluated.}

\item{ratio_sd}{a numeric value indicating the ratio between the standard
deviation in the first group and in the second group. This argument is
needed in the case of Welch's \emph{t}-test.}

\item{B}{a numeric  value indicating the number of iterations. Increase the
number of iterations to obtain more stable results.}

\item{tl}{optional value indicating the lower truncation point if
\code{effect_size} is defined as a function.}

\item{tu}{optional value indicating the upper truncation point if
\code{effect_size} is defined as a function.}

\item{B_effect}{a numeric  value indicating the number of sampled effects
if \code{effect_size} is defined as a function. Increase the number to
obtain more stable results.}

\item{display_message}{a logical variable indicating whether to display or not
the progress bar. Not that this applies only when \code{effect_size} is
defined as a function.}
}
\value{
A list with class "design_analysis" containing the following
 components:
   \item{design_analysis}{a character string indicating the type of design
   analysis: "retrospective".}
   \item{call_arguments}{a list with all the arguments passed to the
   function and the raw function call.}
   \item{effect_info}{a list with all the information regarding the
   considered hypothetical population effect size. The list includes:
   \code{effect_type} indicating the type of effect; \code{effect_function}
   indicating the function from which effect are sampled or the string
   "single_value" if a single value was provided; \code{effect_summary}
   summary of the sampled effects; \code{effect_samples} vector with the
   sampled effects (or unique value in the case of a single value). if
   relevant \code{tl} and \code{tu} specifying the lower upper truncation
   point respectively.}
   \item{test_info}{a list with all the information regarding the test
   performed. The list includes: \code{test_method} character sting
   indicating the test method (i.e., "pearson", "one_sample", "paired",
   "two_sample", or "welch"); sample size (\code{sample_n1} and if relevant
   \code{sample_n2}), alternative hypothesis (\code{alternative}),
   significance level (\code{sig_level})  and  degrees of freedom (\code{df})
   of the statistical test; \code{critical_effect} the minimum absolute
   effect value that would result significant. Note that
   \code{critical_effect} in the case of \code{alternative = "two_sided"} is
   the absolute value and both positive and negative values should be
   considered.}
   \item{retrospective_res}{a data frame with the results of the design
   analysis. Columns names are \code{power}, \code{typeM}, and \code{typeS}.}
}
\description{
Given the hypothetical population effect size and the study sample size, the
function \code{retrospective()} performs a retrospective design analysis for
Pearson's correlation test between two variables or \emph{t}-test comparing
group means (Cohen's \emph{d}). According to the defined alternative
hypothesis and the significance level, inferential risks (i.e., Power level,
Type M error, and Type S error) are computed together with the critical
effect value (i.e., the minimum absolute effect size value that would result
significant).
}
\details{
Conduct a retrospective design analysis to evaluate inferential risks
 according to study design. A general overview is provided in the
 \code{vignette("retrospective")}.

 \strong{Population effect size}

 The hypothetical population effect size (\code{effect_size}) can be set to a
 single value or a function that allows sampling values from a given
 distribution. The function has to be defined as \code{function(n)
 my_function(n, ...)}, with only one single argument \code{n} representing
 the number of sampled values (e.g., \code{function(n) rnorm(n, mean = 0, sd
 = 1)}; \code{function(n) sample(c(.1,.3,.5), n, replace = TRUE)}). This
 allows users to define hypothetical effect size distribution according to
 their needs.

 Argument \code{B_effect} allows defining the number of sampled effects.
 Users can access sampled effects in the \code{effect_info} list included in
 the output to evaluate if the sample is representative of their
 specification. Increase the number to obtain more accurate results but it
 will require more computational time (default is 1000). To avoid long
 computational times, we suggest adjusting \code{B} when using a function to
 define the hypothetical population effect size.

 Optional arguments \code{tl} and \code{tu} allow truncating the sampling
 distribution specifying the lower truncation point and upper  truncation
 point respectively. Note that if \code{effect_type = "correlation"},
 distribution is automatically truncated between -1 and 1.

 \strong{Test methods}

 The function \code{retrospective()} performs a retrospective design analysis
 considering correlations between two variables or comparisons between group
 means.

 In the case of a correlation, only Pearson's correlation between two
 variables is available, whereas Kendall's \emph{tau} and Spearman's
 \emph{rho} are not implemented. The \code{test_method} argument has to be
 set to \code{"pearson"} (default) and the \code{effect_size} argument is
 used to define the hypothetical population effect size in terms of Pearson's
 correlation coefficient (\eqn{\rho}). The \code{sample_n2} argument is
 ignored.

 In the case of a comparison between group means, the \code{effect_size}
 argument is used to define the hypothetical population effect size in terms
 of Cohen's \emph{d} and the available \emph{t}-tests are selected specifying
 the argument \code{test_method}. For independent two-sample \emph{t}-test,
 use \code{"two_sample"} and indicate the sample size of the second group
 (\code{sample_n2}). For Welch's \emph{t}-test, use \code{"welch"} and
 indicate and indicate the sample size of the second group (\code{sample_n2})
 and the ratio between the standard deviation in the first group and in the
 second group (\code{ratio_sd}). For dependent \emph{t}-test for paired
 samples, use \code{"paired"} (\code{sample_n1} and \code{sample_n2} have to
 be equal). For one-sample \emph{t}-test, use \code{"one_sample"}
 (\code{sample_n2} has to be \code{NULL}).

 \strong{Study design}

 Study design can be further defined according to statistical test
 directionality and required \eqn{\alpha}-level using the arguments
 \code{alternative} and \code{sig_level} respectively.
}
\examples{

# Pearson's correlation
retrospective(effect_size = .3, sample_n1 = 25, test_method = "pearson")

# Two-sample t-test
retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 35,
              test_method = "two_sample")
# Welch t-test
retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 35,
              test_method = "welch", ratio_sd = 1.5)
# Paired t-test
retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 25,
              test_method = "paired")
# One-sample t-test
retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = NULL,
              test_method = "one_sample")




\donttest{
# Define effect_size using functions (long computational times)
# Remember to adjust B
retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 25,
              test_method = "pearson", tl = .15, B = 1e3)
retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 25,
              test_method = "one_sample", tl = .2, tu = .4, B = 1e3)
}

}
\references{
Altoè, G., Bertoldo, G., Zandonella Callegher, C., Toffalini, E.,
 Calcagnì, A., Finos, L., & Pastore, M. (2020). Enhancing Statistical
 Inference in Psychological Research via Prospective and Retrospective Design
 Analysis. Frontiers in Psychology, 10.
 \url{https://doi.org/10.3389/fpsyg.2019.02893}

 Bertoldo, G., Altoè, G., & Zandonella Callegher, C. (2020).
 Designing Studies and Evaluating Research Results: Type M and Type S Errors
 for Pearson Correlation Coefficient. Retrieved from
 \url{https://psyarxiv.com/q9f86/}

 Gelman, A., & Carlin, J. (2014). Beyond Power Calculations: Assessing Type S
 (Sign) and Type M (Magnitude) Errors. Perspectives on Psychological Science,
 9(6), 641–651. \url{https://doi.org/10.1177/1745691614551642}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRDA-package.R
\docType{package}
\name{PRDA}
\alias{PRDA}
\title{PRDA: Prospective and Retrospective Design Analysis.}
\description{
Given an hypothetical value of effect size, {PRDA} performs a prospective
or retrospective design analysis to evaluate the inferential risks (i.e.,
power, Type M error, and Type S error) related to the study design. See
\code{vignette("PRDA")} for a brief introduction to \emph{Design
Analysis}.
}
\details{
PRDA package can be used for Pearson's correlation between two variables
or mean comparisons (i.e., one-sample, paired, two-sample, and Welch's
t-test) considering an hypothetical value of \eqn{\rho} or Cohen's \emph{d}
respectively. See \code{vignette("retrospective")} for more details.
}
\section{Functions}{

In {PRDA} there are two main functions:
\itemize{
\item{\strong{\code{retrospective()}}}. Given the hypothetical population
effect size and the study sample size, the function \code{retrospective()}
performs a retrospective design analysis. According to the defined
alternative hypothesis and the significance level, the inferential risks
(i.e., Power level, Type M error, and Type S error) are computed together
with the critical effect value (i.e., the minimum absolute effect size value
that would result significant). To know more about function arguments and
examples see the function documentation
\code{\link[PRDA:retrospective]{?retrospective}} and
\code{vignette("retrospective")}.

\item{\strong{\code{prospective()}}}. Given the hypothetical population
effect size and the required power level, the function \code{prospective()}
performs a prospective design analysis. According to the defined alternative
hypothesis and the significance level, the required sample size is computed
together with the associated Type M error, Type S error, and the critical
effect value (i.e., the minimum absolute effect size value that would
result significant).  To know more about function arguments and examples see
the function documentation \code{\link[PRDA:prospective]{?prospective}}
and \code{vignette("prospective")}.
}
}

\section{Hypothetical Effect Size}{

The hypothetical population effect size can be defined as a single value
according to previous results in the literature or experts indications.
Alternatively, {PRDA} allows users to specify a distribution of plausible
values to account for their uncertainty about the hypothetical population
effect size.  To know how to specify the hypothetical effect size according
to a distribution and an example of application see
\code{vignette("retrospective")}.
}

\references{
Altoè, G., Bertoldo, G., Zandonella Callegher, C., Toffalini, E.,
 Calcagnì, A., Finos, L., & Pastore, M. (2020). Enhancing Statistical
 Inference in Psychological Research via Prospective and Retrospective Design
 Analysis. Frontiers in Psychology, 10.
 \url{https://doi.org/10.3389/fpsyg.2019.02893}

 Bertoldo, G., Altoè, G., & Zandonella Callegher, C. (2020, June 15).
 Designing Studies and Evaluating Research Results: Type M and Type S Errors
 for Pearson Correlation Coefficient. Retrieved from
 \url{https://psyarxiv.com/q9f86/}

 Gelman, A., & Carlin, J. (2014). Beyond Power Calculations: Assessing Type S
 (Sign) and Type M (Magnitude) Errors. Perspectives on Psychological Science,
 9(6), 641–651. \url{https://doi.org/10.1177/1745691614551642}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prospective.R
\name{prospective}
\alias{prospective}
\title{Prospective Design Analysis}
\usage{
prospective(
  effect_size,
  power,
  ratio_n = 1,
  test_method = c("pearson", "two_sample", "welch", "paired", "one_sample"),
  alternative = c("two_sided", "less", "greater"),
  sig_level = 0.05,
  ratio_sd = 1,
  B = 10000,
  tl = -Inf,
  tu = Inf,
  B_effect = 1000,
  sample_range = c(2, 1000),
  eval_power = c("median", "mean"),
  tol = 0.01,
  display_message = TRUE
)
}
\arguments{
\item{effect_size}{a numeric value or function (see Details) indicating the
hypothetical population effect size.}

\item{power}{a numeric value indicating the required power level.}

\item{ratio_n}{a numeric value indicating the ratio between the sample size in
the first group and in the second group. This argument is required when
\code{test_method} is set to \code{"two_sample"} or \code{"welch"}. In the
case of \code{test_method = "paired"}, set \code{ratio_n} to 1. Whereas in
the case of \code{test_method = "one_sample"}, set \code{ratio_n} to
\code{NULL}. This argument is ignored for \code{test_method = "pearson"}.
See Test methods section in Details.}

\item{test_method}{a character string specifying the test type, must be one of
\code{"pearson"} (default, Pearson's correlation), \code{"two_sample"}
(independent two-sample \emph{t}-test), \code{"welch"} (Welch's
\emph{t}-test), \code{"paired"} (dependent \emph{t}-test for paired
samples), or \code{"one_sample"} (one-sample \emph{t}-test). You can specify
just the initial letters.}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of "two_sided" (default), "greater" or "less". You can specify
just the initial letter.}

\item{sig_level}{a numeric value indicating the significance level on which
the alternative hypothesis is evaluated.}

\item{ratio_sd}{a numeric value indicating the ratio between the standard
deviation in the first group and in the second group. This argument is
required only in the case of Welch's \emph{t}-test.}

\item{B}{a numeric  value indicating the number of iterations. Increase the
number of iterations to obtain more stable results.}

\item{tl}{optional value indicating the lower truncation point if
\code{effect_size} is defined as a function.}

\item{tu}{optional value indicating the upper truncation point if
\code{effect_size} is defined as a function.}

\item{B_effect}{a numeric  value indicating the number of sampled effects
if \code{effect_size} is defined as a function. Increase the number to
obtain more stable results.}

\item{sample_range}{a length-2 numeric vector indicating the minimum and
maximum sample size of the first group (\code{sample_n1}).}

\item{eval_power}{a character string specifying the function used to summarize
the resulting distribution of power values. Must be one of "median"
(default) or "mean". You can specify just the initial letters. See Details.}

\item{tol}{a numeric value indicating the tolerance of required power level.}

\item{display_message}{a logical variable indicating whether to display or not
the information about computational steps and the progress bar. Not that the
progress bar is available only when \code{effect_size} is defined as a
function.}
}
\value{
A list with class "design_analysis" containing the following
 components:
   \item{design_analysis}{a character string indicating the type of design
   analysis: "prospective".}
   \item{call_arguments}{a list with all the arguments passed to the
   function and the raw function call.}
   \item{effect_info}{a list with all the information regarding the
   considered hypothetical population effect size. The list includes:
   \code{effect_type} indicating the type of effect; \code{effect_function}
   indicating the function from which effect are sampled or the string
   "single_value" if a single value was provided; \code{effect_summary}
   summary of the sampled effects; \code{effect_samples} vector with the
   sampled effects (or unique value in the case of a single value); if
   relevant \code{tl} and \code{tu} specifying the lower upper truncation
   point respectively.}
   \item{test_info}{a list with all the information regarding the test
   performed. The list includes: \code{test_method} character sting
   indicating the test method (i.e., "pearson", "one_sample", "paired",
   "two_sample", or "welch"); the required sample size (\code{sample_n1} and
   if relevant \code{sample_n2}), the alternative hypothesis
   (\code{alternative}), significance level (\code{sig_level})  and  degrees
   of freedom (\code{df}) of the statistical test; \code{critical_effect} the
   minimum absolute effect value that would result significant. Note that
   \code{critical_effect} in the case of \code{alternative = "two_sided"} is
   the absolute value and both positive and negative values should be
   considered.}
   \item{prospective_res}{a data frame with the results of the design
   analysis. Columns names are \code{power}, \code{typeM}, and \code{typeS}.}
}
\description{
Given the hypothetical population effect size and the required power level,
the function \code{prospective()} performs a prospective design analysis for
Pearson's correlation test between two variables or \emph{t}-test comparing
group means (Cohen's \emph{d}). According to the defined alternative
hypothesis and the significance level, the required sample size is computed
together with the associated Type M error, Type S error, and the critical
effect value (i.e., the minimum absolute effect size value that would
result significant).
}
\details{
Conduct a prospective design analysis to define the required sample
  size and the associated inferential risks according to study design. A
  general overview is provided in the \code{vignette("prospective")}.

  \strong{Population effect size}

  The hypothetical population effect size (\code{effect_size}) can be set to
  a single value or a function that allows sampling values from a given
  distribution. The function has to be defined as \code{function(n)
  my_function(n, ...)}, with only one single argument \code{n} representing
  the number of sampled values (e.g., \code{function(n) rnorm(n, mean = 0, sd
  = 1)}; \code{function(n) sample(c(.1,.3,.5), n, replace = TRUE)}). This
  allows users to define hypothetical effect size distribution according to
  their needs.

  Argument \code{B_effect} allows defining the number of sampled effects.
  Users can access sampled effects in the \code{effect_info} list included in
  the output to evaluate if the sample is representative of their
  specification. Increase the number to obtain more accurate results but it
  will require more computational time (default is 1000). To avoid long
  computational times, we suggest adjusting \code{B} when using a function to
  define the hypothetical population effect size.

  Optional arguments \code{tl} and \code{tu} allow truncating the sampling
  distribution specifying the lower truncation point and upper truncation
  point respectively. Note that if \code{effect_type = "correlation"},
  distribution is automatically truncated between -1 and 1.

  When a distribution of effects is specified, a corresponding distribution
  of power values is obtained as result. To evaluate whether the required
  level of power is obtained, user can decide between the median or the mean
  value as a summary of the distribution using the argument
  \code{eval_power}. They answer two different questions. Which is the
  required sample size to obtain 50% of the time a power equal or greater
  than the required level (median)?; Which is the required sample size to
  obtain on average a power equal or greater than the required level (mean)?

  \strong{Test methods}

  The function \code{retrospective()} performs a retrospective design
  analysis considering correlations between two variables or comparisons
  between group means.

  In the case of a correlation, only Pearson's correlation between two
  variables is available, whereas Kendall's \emph{tau} and Spearman's
  \emph{rho} are not implemented. The \code{test_method} argument has to be
  set to \code{"pearson"} (default) and the \code{effect_size} argument is
  used to define the hypothetical population effect size in terms of
  Pearson's correlation coefficient (\eqn{\rho}). The \code{ratio_n}
  argument is ignored.

  In the case of a comparison between group means, the \code{effect_size}
  argument is used to define the hypothetical population effect size in terms
  of Cohen's \emph{d} and the available \emph{t}-tests are selected
  specifying the argument \code{test_method}. For independent two-sample
  \emph{t}-test, use \code{"two_sample"} and indicate the ratio between the
  sample size of the first group and the second group (\code{ratio_n}). For
  Welch's \emph{t}-test, use \code{"welch"} and indicate the ratio between
  the sample size of the first group and the second group (\code{ratio_n})
  and the ratio between the standard deviation in the first group and in the
  second group (\code{ratio_sd}). For dependent \emph{t}-test for paired
  samples, use \code{"paired"} (\code{ratio_n} has to be 1). For one-sample
  \emph{t}-test, use \code{"one_sample"} (\code{ratio_n} has to be
  \code{NULL}).

  \strong{Study design}

  Study design can be further defined according to statistical test
  directionality and required \eqn{\alpha}-level using the arguments
  \code{alternative} and \code{sig_level} respectively.
}
\examples{

# Pearson's correlation
prospective(effect_size = .3, power = .8, test_method = "pearson", B = 1e3)

# Two-sample t-test
prospective(effect_size = .3, power = .8, ratio_n = 1.5,
            test_method = "two_sample", B = 1e3)
# Welch t-test
prospective(effect_size = .3, power = .8, ratio_n = 2,
            test_method = "welch", ratio_sd = 1.5, B = 1e3)
# Paired t-test
prospective(effect_size = .3, power = .8, ratio_n = 1,
            test_method = "paired", B = 1e3)
# One-sample t-test
prospective(effect_size = .3, power = .8, ratio_n = NULL,
            test_method = "one_sample", B = 1e3)



\donttest{
# Define effect_size using functions (long computational time)
prospective(effect_size = function(n) rnorm(n, .3, .1), power = .8,
            test_method = "pearson", B_effect = 500, B = 500, tl = .15)
prospective(effect_size = function(n) rnorm(n, .3, .1), power = .8,
            test_method = "two_sample", ratio_n = 1, B_effect = 500, B = 500,
            tl = .2, tu = .4)
}

}
\references{
Altoè, G., Bertoldo, G., Zandonella Callegher, C., Toffalini, E.,
 Calcagnì, A., Finos, L., & Pastore, M. (2020). Enhancing Statistical
 Inference in Psychological Research via Prospective and Retrospective Design
 Analysis. Frontiers in Psychology, 10.
 \url{https://doi.org/10.3389/fpsyg.2019.02893}

 Bertoldo, G., Altoè, G., & Zandonella Callegher, C. (2020).
 Designing Studies and Evaluating Research Results: Type M and Type S Errors
 for Pearson Correlation Coefficient. Retrieved from
 \url{https://psyarxiv.com/q9f86/}

 Gelman, A., & Carlin, J. (2014). Beyond Power Calculations: Assessing Type S
 (Sign) and Type M (Magnitude) Errors. Perspectives on Psychological Science,
 9(6), 641–651. \url{https://doi.org/10.1177/1745691614551642}
}
