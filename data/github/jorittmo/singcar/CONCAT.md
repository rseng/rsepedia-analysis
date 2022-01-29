
# singcar

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/singcar)](https://CRAN.R-project.org/package=singcar)
[![CRAN downloads total](http://cranlogs.r-pkg.org/badges/grand-total/singcar)](https://CRAN.R-project.org/package=singcar)
[![CRAN downloads month](https://cranlogs.r-pkg.org/badges/singcar)](https://CRAN.R-project.org/package=singcar)
[![Travis build status](https://travis-ci.com/jorittmo/singcar.svg?branch=master)](https://travis-ci.com/github/jorittmo/singcar)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

The aim of the R package `singcar` is to provide and encourage usage of
appropriate statistical methods for comparing a case against a control sample.
For instance, they may commonly be done in a neuropsychological context, in
which an individual has incurred a specific brain injury and we wish to test
whether this damage has led to an impairment of some cognitive function and
whether two different functions are dissociable. For many functions there is
normed data available which the patient can be compared against directly.
However, when this is not possible a control sample estimating the population,
against which we wish to compare the patient, must be used. Both frequentist and
Bayesian methods have been developed to do this, first and foremost by John
Crawford and Paul Garthwaite (Crawford et al., 2011; Crawford & Garthwaite,
2002, 2007, 2005; Crawford & Howell, 1998). It is these methods that `singcar`
implements. Power calculators for these tests are also provided. Although the
canonical applications for these tests are in Cognitive Neuropsychology or
Clinical Neuropsychology, they are potentially applicable to any circumstance in
which a measure taken from a single individual is to be compared against data
from a normative sample (i.e. a control group). It should be noted that these
statistical methods could also be applied as a general method of outlier
detection in small samples. 

To cite this package you can use:

* Rittmo, J, Ö., McIntosh, R, D. (2021). singcar: Comparing single cases to small samples in R. Journal of Open Source Software, 6(68), 3887, https://doi.org/10.21105/joss.03887


## Installation

`singcar` is now available on CRAN! To get this stable version run:

```R
install.packages("singcar")
library("singcar")
```

You can install the unstable(!) developmental version of `singcar` by running the following:

```R
install.packages("devtools")
library("devtools")
install_github("jorittmo/singcar")
library("singcar")
```

## Issues

Please report any bugs, issues or feature requests in the 
[issue tracker](https://github.com/jorittmo/singcar/issues).

## Contributions 

Contributions are highly appreciated and 
[contribution guidelines](CONTRIBUTING.md) prior to submitting a pull request.

## Example

The package comes with the dataset `size_weight_illusion`, a neuropsychological
dataset from an investigation of the size-weight illusion in DF, a patient with
visual form agnosia following following bilateral lesions to the lateral
occipital complex (Hassan et al., 2020). It was investigated whether
DF experienced visual size-weight illusion to the same extent as controls (n = 28)
and whether visual and kinaesthetic size-weight illusion could be dissociable.
Below follows examples of how to analyse this dataset using the tests provided
in `singcar`.

### Testing for a deficit

If we want to assess whether DF has an impairment compared to controls on
visual size-weight illusion we can test this using a modified two-sample t-test,
called TD (test of deficit: Crawford & Howell, 1998).

``` r
library(singcar)
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

TD(case = DF_V_SWI, controls = CON_V_SWI, conf_int = TRUE)

##  	Crawford-Howell (1998) t-test

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## t = -1.7243, df = 27, p-value = 0.04804
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Standardised case score (Z-CC), 95% CI [-2.34, -1.15]       Proportion below case (%), 95% CI [0.95, 12.47] 
##                                             -1.754857                                              4.804003 

```
This can similarly be tested with a Bayesian version of the same test,
yielding approximately (since this test is based on MCMC methods) the same output (Crawford & Garthwaite, 2007).

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient 
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

BTD(case = DF_V_SWI, controls = CON_V_SWI)

##  	Bayesian Test of deficit by Crawford and Garthwaite (2007)

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## est. z = -1.7388, df = 27, p-value = 0.04803
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Std. case score (Z-CC), 95% credible interval [-2.34, -1.15]    Proportion below case (%), 95% credible interval [0.96, 12.44] 
##                                                -1.754857                                                          4.802900

```
If the control sample for a study is not appropriately matched to the case on
variables such as e.g. age or education level it is appropriate to use tests
that account for this by allowing for the inclusion of covariates. Including
theoretically sound covariates is often a good idea. To do this
Crawford et al. (2011) extended their Bayesian verison of the TD. This test
assess the patient on the task of interest by essentially comparing him/her to
the controls with the same score on the covariate.

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting the coviariate below
DF_age <- size_weight_illusion[size_weight_illusion$PPT == "DF", "YRS"] # Patient
CON_age <- size_weight_illusion[size_weight_illusion$PPT != "DF", "YRS"] # Controls

## BTD_cov(case_task = DF_V_SWI, case_covar = DF_age, control_task = CON_V_SWI, control_covar = CON_age)

##  	Bayesian Test of deficit with Covariates

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## est. z = -1.6828, p-value = 0.05386
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Std. case difference (Z-CCC), 95% credible interval [-2.31, -1.10]     Proportion below case (%), 95% credible interval [1.05, 13.52] 
##                                                          -1.749556                                                           5.386088 

```
### Testing for a dissociation

If we want to assess whether DF has a dissociation between two functions we can
use a modified paired samples t-test to assess the size of the difference
between the case scores from the two tasks to the distribution of differences
between the tasks in the controls. This can however only be done directly using
the t-distribution if the tasks are measured on the same scale and is called the
unstandardised difference test (UDT: Crawford & Garthwaite, 2005). In the
`size_weight_illusion` dataset it is possible to use this test to whether
patient DF exhibits a dissociation between visual size-weight illusion and
kinaesthetic size-weight illusion because the visual and kinaesthetic conditions
are parallel versions of the same task, with different sensory cues. This would
be done as shown below:

``` r
library(singcar)
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinaesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

UDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI)

##  	Unstandardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## t = -0.6667, df = 27, p-value = 0.5106
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                       Standardised case score, task A (Z-CC)                       Standardised case score, task B (Z-CC) 
##                                                  -0.13647439                                                  -0.07931545 
## Standardised task discrepancy (Z-DCC), 95% CI [-1.53, -0.59]              Proportion below case (%), 95% CI [6.35, 27.68] 
##                                                  -1.06478887                                                  25.53097678  

```

Most often this is not possible because we wish to estimate abnormality of 
discrepancy on tasks that are not comparable. So otherwise, that
is if the scores must be standardised to be comparable, a statistic that
approximates the t-distribution has been developed and should be used (the
revised standardised difference test RSDT: Crawford & Garthwaite, 2005). 
The visual and kinaesthetic size-weight illusion will be used for illustrative
purposes here as well:

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinaesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

RSDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI)

##  	Revised Standardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## approx. abs. t = 1.015, df = 27, p-value = 0.3191
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                         Case score on task A as standard (z) score                         Case score on task B as standard (z) score 
##                                                        -1.7548574                                                         -0.7836956 
##  Std. effect size (Z-DCC) for task diff. between case and controls Proportion of control population with more extreme task difference 
##                                                         -1.0647889                                                         15.9560625 


```
A Bayesian version of this test was also developed (Crawford & Garthwaite, 2005),
however, unlike `TD` and `BTD` the `RSDT` and `BSDT` (Bayesian standardised
difference test) differ somewhat and `BSDT` has been shown to keep a better
control of Type I errors if a patient exhibits extreme deficits on both tasks of
interest. Therefore the `BSDT` is recommended above `RSDT`. The usage of the two
R functions is very similar. Since the `BSDT` is based on MCMC methods it can be
quite computationally intensive, depending on the number of iterations you choose.


``` r

# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinaesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

BSDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI, iter = 10^6)

##  	Bayesian Standardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## est. z = -1.0736, p-value = 0.3067
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                                                              Case score on task A as standard (z) score 
##                                                                                              -1.7548574 
##                                                              Case score on task B as standard (z) score 
##                                                                                              -0.7836956 
## Std. effect size (Z-DCC) for task diff. between case and controls, 95% credible interval [-1.71, -0.45] 
##                                                                                              -1.0647889 
## Proportion of control population with more extreme task difference, 95% credible interval [4.32, 32.47] 
##                                                                                              15.3357743 


```

Just as for `BTD` a version of `BSDT` allowing for
covariates has been developed. This test assess the patient on the discrepancy
between the tasks of interest by essentially comparing him/her to the controls
with the same score on the covariate. 

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

# Extracting the coviariate below
DF_age <- size_weight_illusion[size_weight_illusion$PPT == "DF", "YRS"] # Patient
CON_age <- size_weight_illusion[size_weight_illusion$PPT != "DF", "YRS"] # Controls

BSDT_cov(case_tasks = c(DF_V_SWI, DF_K_SWI ), case_covar = DF_age,
         control_tasks = cbind(CON_V_SWI, CON_K_SWI), control_covar = CON_age, iter = 10^5)
         
##  	Bayesian Standardised Difference Test with Covariates

## data:  Case score Y1: 0.03, Case score Y2: 0.10, Controls score Y1: 0.16, Controls score Y2: 0.18
## ave. z = -1.0229, p-value = 0.3303
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                                                               Case score on task X as standard (z) score 
##                                                                                                -1.754857 
##                                                               Case score on task Y as standard (z) score 
##                                                                                                -0.783696 
## Std. effect size (Z-DCCC) for task diff. between case and controls, 95% credible interval [-1.67, -0.40] 
##                                                                                                -1.064152 
##  Proportion of control population with more extreme task difference, 95% credible interval [4.77, 34.54] 
##                                                                                                16.510000       

```
All of the functions above can also take summary (mean, sd, control sample size) data as input.

### Power calculators

A further capacity of `singcar` is that it can be used to calculate power for
for these single case-control comparisons. Calculations for all Bayesian tests
and `RSDT` are simulation based and (especially the tests with covariates) can
be computationally intense. Calculators for `TD` and `UDT` (unstandardised
difference test) are exact (their power functions have been derived
analytically) and can both be used to find a specific sample size given a
desired power. For the other calculators all parameters must be given. Means and
standard deviations for the control population are at default set to 0 and 1
meaning that the case value will be interpreted as differences from the mean in
standard deviations, these parameter values can be changed as you like. Examples
are given below:

``` r
TD_power(case = -2, power = 0.8, mean = 0, sd = 1, alternative = "two.sided")

## [1] "Power(0.44280) will not increase more than 0.5% for any additional participant over n = 16"
##    n     power
## 1 16 0.4428042

TD_power(case = 70, sample_size = 10, mean = 100, sd = 15, alternative = "less", alpha = 0.1)

## [1] 0.7039033

RSDT_power(case_a = 70, case_b = 20, mean_a = 100, mean_b = 25, sd_a = 15, sd_b = 10, sample_size = 10) 

## [1] 0.5689


# Takes long time to compute
BTD_cov_power(case = -2, case_cov = 0, control_task = c(0, 1),
              control_covar = c(0, 1), cor_mat = diag(2), sample_size = 10)

# [1] 0.511

```

# References

Crawford, J., & Garthwaite, P. (2002). Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. Neuropsychologia, 40(8), 1196-1208. https://doi.org/10.1016/S0028-3932(01)00224-X

Crawford, J., & Garthwaite, P. (2007). Comparison of a single case to a control or normative sample in neuropsychology: Development of a Bayesian approach. Cognitive Neuropsychology, 24(4), 343-372. https://doi.org/10.1080/02643290701290146

Crawford, J., & Garthwaite, P. (2005). Testing for Suspected Impairments and Dissociations in Single-Case Studies in Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and Revised Tests for Dissociations. Neuropsychology, 19(3), 318-331. https://doi.org/10.1037/0894-4105.19.3.318

Crawford, J., Garthwaite, P., & Ryan, K. (2011). Comparing a single case to a control sample: Testing for neuropsychological deficits and dissociations in the presence of covariates. Cortex, 47(10), 1166-1178. https://doi.org/10.1016/j.cortex.2011.02.017

Crawford, J., & Howell, D. (1998). Comparing an Individual's Test Score Against Norms Derived from Small Samples. The Clinical Neuropsychologist, 12(4), 482-486. https://doi.org/10.1076/clin.12.4.482.7241

Hassan, E. K., Sedda, A., Buckingham, G., & McIntosh, R. D. (2020). The size-weight illusion in visual form agnosic patient DF. Neurocase, 1-8. https://doi.org/10.1080/13554794.2020.1800748
# singcar 0.1.0

* Released to CRAN 
* Updates for acceptance:
  - Added appropriate references and updated the description.
  - Removed a function that could be used to modify the global seed.
  - Updated documentation slightly.
  - Decreased runtime on examples.
  - Updated automatic tests.

# singcar 0.1.1

* Added a `NEWS.md` file to track changes to the package.
* Removed test dependent on RNG which failed irregularly.

# singcar 0.1.2

* Print output is now more compact.
* Confusing statistics without interpretative value have been removed from the output.
* `type = cairo` removed from vignette since Cairo-based devices are not always available.
* Fixed bug that made power calculators for tests of deficits behave unexpectedly.
* Exchanged URLs for DOI-links in documentation references.

# singcar 0.1.3

* Changed default correlation matrices for Bayesian power calculators.
* Vignette now updated with extensive description of the methods used.
* Edited print output.
* Source code has been more rigorously commented. 

# singcar 0.1.4

* Added paper for Journal of Open Source Software
* Added contribution guidelines and Code of Conduct


# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual
identity and orientation.

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

* The use of sexualized language or imagery, and sexual attention or advances of
  any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email address,
  without their explicit permission
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
version 2.1, available at
[https://www.contributor-covenant.org/version/2/1/code_of_conduct.html][v2.1].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available at
[https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.1]: https://www.contributor-covenant.org/version/2/1/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contributing to `singcar`

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to `singcar`!

`singcar` is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/jorittmo/singcar
[docs]: https://cran.r-project.org/web/packages/singcar/singcar.pdf
[issues]: https://github.com/jorittmo/singcar/issues
[new_issue]: https://github.com/jorittmo/singcar/issues/new
[email]: mailto:j.rittmo@gmail.com

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think `singcar` is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using `singcar` for a paper you are writing? Consider citing it, e.g. by using `utils::citation("singcar")`.

### Ask a question️

Using `singcar` and got stuck? Browse the [documentation][docs] or [vignette](https://cran.r-project.org/web/packages/singcar/vignettes/singcar_vignette.html) to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. I will try my best to address it. 

Want to ask a question in private? Contact the package maintainer by [email][email].

### Feature suggestions

Have an idea for a new `singcar` feature? Take a look at the [documentation][docs] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue].

### Report a bug or issue

Report bugs or issues as a new [new issue][new_issue] in the issue tracker, so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your session info (e.g. using `utils::sessionInfo()`).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed (but minimal) steps to reproduce the bug, preferably with a reproducible example using [reprex](https://reprex.tidyverse.org/). 

### Improve the documentation

Noticed a typo in the vignette? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code

Care to fix bugs or implement new functionality for `singcar`? Have a look at the [issue tracker][issues], leave a comment on things you are
working on and see the pull request guidelines below.

# Pull request guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

### Reference

These contribution guidelines have been adopted and adapted from [here](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c).
## Resubmission
This is a resubmission. In this version I have:

* Updated the vignette with extensive description of the methods used.
* Changed default correlation matrices for Bayesian power calculators.
* Slightly edited print output.
* Commented source code more rigorously. 

## Test environments
- local Windows 10 R installation, R 4.0.4
- ubuntu 16.04 (on travis-ci), R 4.0.4
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
0 errors | 0 warnings | 0 notes 

## Downstream dependencies
There are currently no downstream dependencies for this package.
---
title: '`singcar`: Comparing single cases to small samples in `R`'
tags:
  - R
  - Case control comparison
  - Neuropsychology
  - Brain lesion studies
  - Small samples
authors:
  - name: Jonathan Ö. Rittmo
    orcid: 0000-0001-5075-0166
    affiliation: 1
  - name: Robert D. McIntosh
    orcid: 0000-0002-7615-6699
    affiliation: 1
affiliations:
 - name: Human Cognitive Neuroscience, Psychology, University of Edinburgh, UK
   index: 1
date: 11 October 2021
bibliography: paper.bib

---

# Summary

Case-control comparisons are a class of statistical tests allowing researchers
to compare single cases to populations estimated from a sample. Such tests have
wide potential utility, but historically have been applied mostly in the fields
of cognitive and clinical neuropsychology, to infer whether individuals have
suffered significant cognitive changes as the consequence of a brain lesion. One
may wish to estimate whether that individual has abnormally low performance on
some cognitive ability, or if one cognitive ability is abnormally discrepant
with respect to another cognitive ability. John Crawford, Paul Garthwaite and
colleagues have developed several related methods to statistically test for
abnormality on a single variate and abnormality of the difference between two
variates when a single case is compared to a small sample, while controlling the
Type I error rate [e.g.,@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
This paper presents the `R` package `singcar` in which they are implemented. Due to recent discussion
on the fundamental power limits of these tests [@mcintoshPowerCalculationsSinglecase2020] the package also includes 
associated power calculators. 


# Statement of need

There are many reasons why researchers and clinicians might want to look at
single cases instead of at the average of some group. In certain fields, such as
neuropsychology, this need arises because the pattern of naturally-occurring
brain damage will be unique in each individual case. From a theoretical
perspective, this means that a single patient might be the only available source
of data for a given phenomenon. 
From a practical, clinical perspective,
diagnosis and description of the pattern of cognitive impairment is done at the
individual level. Individual brain-damaged patients are thus often compared to
the healthy population to assess changes in cognitive functioning. If we want to
assess the patient score on some variate Y, for which we do not know the
population parameters, these must be estimated from a sample. Thus, the
single-case of interest is compared to a control sample. There are many other
areas where the application of such methods could also be useful, for example
studies of uncommon human expertise.

As it represents the canonical field for the application of these methods, the
nomenclature of neuropsychology is adopted. An abnormally low score on
a single variate is referred to as a *deficit*, an important concept
for clinical and basic neuropsychology alike. For the latter area another
concept is also considered to be of cardinal importance: the ability to test for
an abnormally large discrepancy between two variates. This is referred to as a
*dissociation*, which is taken to provide evidence for some degree of
functional independence between two cognitive abilities. By charting
dissociations, a cognitive architecture of the mind can be theorized
[@shalliceNeuropsychologyMentalStructure1988].

During the last 20 years, a class of related methods have been developed for
case-control comparisons, allowing researchers to estimate abnormality and test
for deficits and dissociations in the single case, while controlling the Type I
error rate. These tests have been developed mainly by John Crawford and Paul
Garthwaite [e.g.,@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
John Crawford has provided free software packages to perform these tests, making
them available at https://homepages.abdn.ac.uk/j.crawford/pages/dept/psychom.htm. However,
these are available only as standalone compiled computer programs for Windows
operating systems. Many of these programs require manual input of summary
statistics, and output a static text file and for thorough documentation one must
consult the original publications. 

Our aim is to encourage and simplify usage of these methods by implementing them
in the package `singcar` for the `R` environment
[@rcoreteamLanguageEnvironmentStatistical2020], bringing them together in a
fully documented package with open source code that works across platforms.
Further advantages of `singcar` include an API that has more modifiable test
parameters. It is also possible to automate these tests if multiple analyses need to
be run for the purposes of data analysis or simulation studies. 
The development of Crawford and Garthwaite's methods has been focused around limiting Type I
errors, but to emphasise the importance of considering Type II errors we also provide power
calculators for each test function. Our hope in doing so is to increase awareness of power
for this methodology as well as to aid in the planning and design of experiments 
[@mcintoshPowerCalculationsSinglecase2020].

Note that the `R` package `singlecase` [@matthieusinglecase2008] contains some overlapping
functionality with `singcar`, but it has not been maintained since 2008 and lacks
core functionality such as tests allowing for the inclusion of covariates, and
power calculators. A recent study by @mitchell2020peripheral investigating 
peripheral reaching in Alzheimer’s disease and mild cognitive impairment exemplifies
the uses of these novel functionalities in `singcar`. 

# Functionality

`singcar` contains seven functions to estimate a
case's abnormality compared to a normal population estimated from a small
sample, three of them with regards to a single variate and four with regards to
the discrepancy between two variates. Both frequentist and Bayesian methods are
provided, all developed originally by Crawford and colleagues
[@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
Of special note for psychological research are the methods allowing the
inclusion of covariates [@crawfordComparingSingleCase2011] using Bayesian
regression techniques. These methods make matching the control sample to the
case less cumbersome. For rationale as well as mathematical and contextual
background of the methods consult the package vignette.

# References
