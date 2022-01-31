
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
---
author: 
  - name: Jonathan Ö. Rittmo
    affiliation: University of Edinburgh
    address: |
      | Human Cognitive Neuroscience
      | Psychology
    email: j.rittmo@gmail.com
  - name: Robert D. McIntosh
    affiliation: University of Edinburgh
    address: |
      | Human Cognitive Neuroscience
      | Psychology
    email: r.d.mcintosh@ed.ac.uk
title: "The R package singcar: Comparing Single Cases to Small Samples"
  # formatted: "The \\proglang{R} package \\pkg{singcar}: Comparing Single Cases to Small Samples"
  # plain:     "The R package singcar: Comparing Single Cases to Small Samples"
  # short:     "\\pkg{singcar}: Comparing Single Cases to Small Samples"
abstract: >
    Comparison of single cases to populations estimated from a sample has many
    potential applications. Historically, one major application is in the field of
    neuropsychology, where single-case statistics may be used to infer whether an
    individual has suffered a significant cognitive deficit as the consequence of a
    brain lesion. One may wish to estimate whether that individual has abnormally
    low performance on some cognitive ability, or if one cognitive ability is
    abnormally discrepant with respect to another cognitive ability. Several
    statistical methods have been developed to test for abnormality on a single
    variate and abnormality of the difference between two variates when a single
    case is compared to a sample, without losing control of the Type I error rate.
    This paper describes some of the main methods and presents a package in which
    they are implemented.
keywords: "single case comparisons, small samples, R"
  # formatted: [single case comparisons, small samples, "\\proglang{R}"]
  # plain:     [single case comparisons, small samples, R]
bibliography: ref.bib
link-citations: true
preamble: >
  \usepackage{amsmath}
  \usepackage{amsfonts}
  \usepackage{longtable}
  \usepackage{caption}
  \usepackage{pbox}
  \usepackage{booktabs}
# \usepackage{cleveref}
# \renewcommand{\eqref}{\Cref}
# \Crefformat{equation}{#2#1#3}
output:
    bookdown::html_document2:
      base_format: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The R package singcar: Comparing Single Cases to Small Samples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, setup, include=FALSE}
options(prompt = 'R> ', continue = '+ ')

library(singcar)
```

# Introduction

There are many reasons why researchers and clinicians might want to look at
single cases instead of at the average of some group. In certain fields, such as
neuropsychology, this need arises because the pattern of naturally-occurring
brain damage will be unique in each individual case. From a theoretical
perspective, this means that a single patient might be the only available source
of data for a given phenomenon. From a practical, clinical perspective,
diagnosis and description of the pattern of cognitive impairment is done at the
individual level. Individual brain-damaged patients are thus often compared to
the healthy population to assess changes in cognitive functioning. If we want to
assess the patient score on some variate Y, for which we do not know the
population parameters, these must be estimated from a sample. Thus, the
single-case of interest is compared to a control sample. There are many other
areas where the application of such methods could also be useful: for example
studies of uncommon human expertise, targeted quality checks in industries with
limited production output, or animal studies where rearing a large experimental
group might be infeasible.

As it represents the canonical field for the application of these methods, the
nomenclature of neuropsychology will here be adopted. An abnormally low score on
a single variate will be referred to as a *deficit*, an important concept
for clinical and basic neuropsychology alike. For the latter area another
concept is also considered to be of cardinal importance: the ability to test for
*dissociation*, which is taken to provide evidence for some degree of
functional independence between two cognitive abilities. By charting
dissociations, a cognitive architecture of the mind can be theorized
[@shalliceNeuropsychologyMentalStructure1988].

During the last 20 years methods have been developed to estimate and test for
both deficits and dissociations in the single case, while controlling the Type I
error rate, mainly by John Crawford and Paul Garthwaite
[e.g.,@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
They are available as standalone computer programs, only taking summary data as
input, at https://homepages.abdn.ac.uk/j.crawford/pages/dept/psychom.htm.
But the majority of these methods have not yet been implemented in any standard
statistical environment. By doing so in the package `singcar` for the `R`
environment [@rcoreteamLanguageEnvironmentStatistical2020], our aim
is to encourage and simplify their usage. Further, limiting Type II errors has
not received as much attention as limiting Type I errors in single-case
methodology [@mcintoshPowerCalculationsSinglecase2020]. By including novel tools
for power calculations, our hope is to increase awareness of the inherently low
power in this field as well as to aid in the planning and design of experiments.

The `R` package `singcar` contains seven functions to estimate a case's
abnormality compared to a normal population estimated from a small sample, three of them
with regards to a single variate and four with regards to the discrepancy between two variates.
Both frequentist and Bayesian methods are provided, all developed originally by
Crawford and colleagues [@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
Of special note for psychological research are the
methods allowing the inclusion of covariates [@crawfordComparingSingleCase2011]
using Bayesian regression techniques (Section \@ref(cov)). These methods make matching the control
sample to the case less cumbersome.

In Section \@ref(section2) the implemented methods are described in detail and in
Section \@ref(section3) the package is described and usage exemplified with a
dataset from a recent neuropsychological single case study.


# Comparing a single case to small samples {#section2}

## Background

In the neuropsychological application of the methods to be presented,
the variates of interest will be scores obtained by participants on tasks
assessing relevant cognitive functions. The variates will
thus often be referred to as task scores. Historically, it was not uncommon for
researchers to compare single cases against a control
population estimated from a sample by evaluating the case score as a
\(Z\) score from the estimated distribution. The \(p\) value associated with this
\(Z\) score would then be treated as an estimate of the case's abnormality. A
similar logic was sometimes applied for estimating the abnormality of the
discrepancy between two tasks, using a method developed by @payneStatisticsInvestigationIndividual1957:
researchers calculated a \(Z\) score based on the case's difference
between two standardised variates divided by the standard deviation of the
difference.

The \(Z\) score approach is of course problematic because it treats parameter estimations from
a restricted control sample as if they were population parameters. This could
be appropriate with very large samples, but
control samples in neuropsychology are often small (sometimes < 10), and it is
well known that the sampling distribution of the sample variance is right skewed
for small sample sizes. Hence, underestimation of variance would be more
probable than overestimation, inflating obtained \(Z\) scores and thus
Type I errors if the Z scores were used for hypothesis testing [@crawfordComparingIndividualTest1998].

The following sections describe methods that have been devised to
allow for the evaluation of the abnormality of a case against a restricted
control sample, whilst retaining appropriate control over the Type I error rate.
These methods include frequentist approaches (Section \@ref(sec22)) and Bayesian
approaches (Section \@ref(sec23)).


## Frequentist approaches {#sec22}

An appropriate method for comparing a single observation to the mean of a
sample was proposed within the biological sciences by @sokalBiometryPrinciplesPractice1981 (p.
227). It was popularized within neuropsychology by @crawfordComparingIndividualTest1998,
where its common application was as a (one-tailed) test of deficit (TD) to
determine whether the score of a single case with brain damage was abnormally
low (assuming that a low score indicates poorer performance) with respect to the
mean score of a control sample. The \(t\) distribution is used to account for
the underestimation of the sample variance. The basic approach is a modified two
samples \(t\) test where the case simply is treated as a sample of size 1. The
degrees of freedom for this distribution is \(n + 1 - 2 = n - 1\).

\begin{equation}
t_{n-1} = \frac{y^* - \overline{y}}{s \sqrt{\frac{n + 1}{n}}}
(\#eq:TD)
\label{eq:TD}
\end{equation}

Where \(y^*\) is the case score, \(\overline{y}, \ s\) the sample mean and sample
standard deviation respectively and \(n\) the size of the control sample. This
method does not allow estimation of a notional patient population. However,
since the question posed when using this test is about
the probability of the case being part of the control population, the potential alternative
affinity of the case is of no concern [@crawfordComparingIndividualTest1998]. This
test of deficit provides transparent control of the Type I error rate [@crawfordComparingSingleCase2009, @crawfordInferentialMethodsComparing2004, @crawfordSinglecaseResearchNeuropsychology2012].
Moreover, the \(p\) value obtained constitutes an unbiased point estimate of the
abnormality of the case, as demonstrated by @crawfordMethodsTestingDeficit2006. The
simple proof for this is given below.

The proportion of controls that would score lower on a random variable \(Y\) than
the case \(y^*\) is
\[
\mathbb{P}[Y< y^*]
\]
Subtracting \(\overline{y}\) from both sides of the inequality
and dividing by \(s \sqrt{\frac{n + 1}{n}}\)
\[
\mathbb{P}[Y< y^*] =\mathbb{P}\left[\frac{Y-\overline{y}}{s \sqrt{\frac{n + 1}{n}}} < \frac{y^* - \overline{y}}{s \sqrt{\frac{n + 1}{n}}}\right]
\]
The quantity to the left of the inequality, i.e., \(\frac{y-\overline{y}}{s \sqrt{\frac{n + 1}{n}}}\)
is \(t\) distributed with \(n-1\) degrees of freedom. Hence,
\[
\mathbb{P}[Y< y^*] =\mathbb{P}\left[t_{n-1} < \frac{y^* - \overline{y}}{s \sqrt{\frac{n + 1}{n}}}\right]
\]
The quantity to the right of the inequality is the test statistic from
Equation \@ref(eq:TD), hence \(\mathbb{P}[y< y^*]\) is the same as the \(p\) value obtained
from the test of deficit. This fact also makes the construction of confidence
intervals of this abnormality possible.

@crawfordPointIntervalEstimates2010 pointed out that, although the \(Z\) score
is not appropriate for estimating the abnormality of a case when the control
sample is small, it does provide a standardised effect size measure of
abnormality, similar to Cohen's d [@cohenStatisticalPowerAnalysis1988]. Insensitivity to sample size is
indeed a requirement for an effect size index and the proposed quantity was:
\begin{equation}
Z_{CC} = \frac{y^* - \overline{y}}{s}
(\#eq:zcc)
\label{eq:zcc}
\end{equation}
This is simply a \(Z\) score calculated from sample estimates. The
subscript (CC) indicates that it relates to a "case-controls" comparison. The
construction of confidence intervals for the point estimate of abnormality using
the \(Z_{CC}\) index of effect size was described by @crawfordInvestigationSingleCase2002.

Put \(p\) as the percentage of the population that would fall below a
case score, then the \(1-\frac{\alpha}{2}\) confidence interval for \(p\) is
constructed as follows: Let \(Z_{CC} = \frac{y^*-\overline{y}}{s}\) and \(y^* \neq \overline{y}\)
then \(Z_{CC}\) comes from a non-central \(t\) distribution on \(n-1\)
degrees of freedom. By deploying a search algorithm we find the value \(\delta_U\)
for the non-centrality parameter (NCP) of a non-central \(t\) distribution on \(n-1\)
degrees of freedom such that the \(100\frac{\alpha}{2}\) percentile equates
\(Z_{CC}\sqrt{n}\) and similarly we find the NCP \(\delta_L\) of a non-central
\(t\) distribution such that its \(100(1-\frac{\alpha}{2})\) percentile equates
\(Z_{CC}\sqrt{n}\). The upper and lower boundaries for \(Z_{CC}\) are then given by:
\[
Z_{CC_U} = \frac{\delta_U}{n}, \ \ Z_{CC_L} = \frac{\delta_L}{n}
\]
and the boundaries for \(p\) by:
\[
p_U = \Phi\left(\frac{\delta_U}{n}\right), \ p_L = \Phi\left(\frac{\delta_L}{n}\right)
\]
Where \(\Phi\) is the CDF of the standard normal distribution. Note that the above
equation is appropriate when the case score falls on the left side of the
distribution. If it were to fall on the right side, then the upper and lower
boundaries would be given by:
\[
p_U = 1 - \Phi\left(\frac{\delta_L}{n}\right), \ p_L = 1 - \Phi\left(\frac{\delta_U}{n}\right)
\]

As has been shown, estimating abnormality on a single variate is a simple
matter. For estimating abnormality of discrepancy, when the normally distributed
variates \(Y_1\) and \(Y_2\) can be compared without standardisation, the approach
is equally simple and in fact identical to the test of deficit, Equation \@ref(eq:TD),
except that the statistic is based on difference scores, taking the correlation
between the variates into account (\(\rho_{12}\)).
\begin{equation}
t_{n-1} = \frac{(y^*_1 - \overline{y}_1) - (y^* _2 - \overline{y}_2) }{ \sqrt{(s^2_1 +s^2_2 -2s_1 s_2 \rho_{12})(\frac{n+1}{n})}}
(\#eq:UDT)
\label{eq:UDT}
\end{equation}
Construction of confidence intervals for abnormality estimates based on this
unstandardised difference test (UDT) is done in an analogous way to that of the
confidence intervals of the TD [@crawfordTestingSuspectedImpairments2005].
However, the test is somewhat limited in its usefulness, because it is only
applicable if the two variates are measured on equivalent scales. Far more
common, at least in neuropsychology, is the need to assess discrepancies between
differently-scaled variates, which require standardisation to be comparable.

By standardising the variates (not taking sample size into account)
we get an effect size measure similar to \(Z_{CC}\), 
Equation \@ref(eq:zcc) [@crawfordPointIntervalEstimates2010]:
\begin{equation}
Z_{DCC}  = \frac{z^*_1 - z^*_2}{\sqrt{2-2\rho_{12}}}
(\#eq:PJ)
\label{eq:PJ}
\end{equation}
Where \(z_1^*\) and \(z_2^*\) are the standardised case scores on some task A
and some task B and the subscript (DCC) indicates "discrepancy-case-controls".
This quantity was proposed by @payneStatisticsInvestigationIndividual1957
as a significance test for discrepancies, but of course it is not appropriate
for such a purpose if small samples are used: we cannot use the \(Z\)
distribution for estimating the abnormality of the discrepancy for the same
reason we cannot use \(Z_{CC}\) to test for a deficit.
<!-- But also because we cannot use  \(Z\) scores to represent the case's scores -->
<!-- on each variate since the size of these will be inflated as well. -->

When the case scores on variates \(Y_1\) and \(Y_2\) are estimated from small
samples and they need to be standardised, it means that instead of estimating
the discrepancy between two normally distributed variates we need to estimate
the discrepancy between two \(t\) distributed variates.
Because linear combinations of correlated \(t\) distributed random variates are not
themselves \(t\) distributed this problem is non-trivial. The distribution of such
difference scores was examined by @garthwaiteDistributionDifferenceTwo2004. They used
asymptotic expansion to find functions of the sample correlation that, if used
as a denominator to the difference between the variates, would yield an
asymptotically \(t\) distributed quantity. They found:

\begin{equation}
\psi=\frac{\frac{(y^*_1-\overline{y}_1)}{s_{1}}-\frac{(y^*_2-\overline{y}_2)}{s_{2}}}{
\sqrt{
(\frac{n+1}{n})
\left( (2-2 \rho)+
\frac{2(1-\rho^{2})}{n-1}+
\frac{(5+c^{2})(1-\rho^{2})}{2(n-1)^{2}}+
\frac{\rho(1+c^{2})(1-\rho^{2})}{2(n-1)^{2}}\right)
}}
\end{equation}

Where \(\rho\) is the sample correlation between the tasks and \(c\) the critical
two-tailed \(t\) value with \(n-1\) degrees of freedom. @garthwaiteDistributionDifferenceTwo2004
demonstrated that \(\mathbb{P}[ \psi > c] \approx \mathbb{P}[t >c]\). To obtain a
precise probability for \(\psi\), one solves for \(\psi = c\), which gives a quantity
not dependent on a pre-specified critical value.
\(\psi = c\) is a quadratic equation in \(c^2\), choosing the positive root of which
yields:

\begin{align}
\begin{split}
c & = \sqrt{\frac{ -b + \sqrt{b^2 - 4ad}}{2a}}, \  \text{where} \\
a & = (1+r)(1-r^2), \\
b & =  (1-r)[4(n-1)^2+4(1+r)(n-1)+(1+r)(5+r)], \\
d & =  - 2\left[\frac{y^*_{1} - \overline{y}_1}{s_1}-\frac{y^*_2 -\overline{y}_2}{s_2}\right]^2\left(\frac{n(n-1)^2}{n+1}\right)
\end{split}
(\#eq:RSDT)
\label{eq:RSDT}
\end{align}

Where \(p = \mathbb{P}[t_{n-1}>c]\) is the estimate of abnormality, and is used for
significance testing. It should be noted that because \(c\) is quadratic and we
choose the positive root, the resultant statistic cannot be negative. 
If it is required that the test statistic expresses the direction of a negative
effect the correct sign must be imposed.

The quantity in Equation \@ref(eq:RSDT) is referred to as the "revised
standardised difference test"
\citep[RSDT:][]{crawfordTestingSuspectedImpairments2005}. The name stems from it
being preceded by a test similar to that of the unstandardised difference test
(Equation \@ref(eq:UDT)) but with standardised normal variates, developed by
@crawfordPayneJonesRevisited1998.
@crawfordTestingSuspectedImpairments2005 show with Monte Carlo
simulations that the RSDT is superior in controlling Type I errors compared to
their earlier $t$ test and the $Z$ score method of
@payneStatisticsInvestigationIndividual1957.

Even for very small sample sizes of \(n=5\), RSDT was shown to barely exceed the
specified 5\% error rate. However, if a case lacks a true discrepancy
between the variates but exhibits extreme scores on both of them (in the same direction), the error rate
of the RSDT increases steeply. @crawfordComparisonSingleCase2007 
showed that the RSDT starts to lose control of the
error for task scores being more extreme than two standard deviations away from
the mean on both variates, without them exhibiting a discrepancy. For task
scores at 8 standard deviations from the mean, the Type I error rate of the RSDT
was inflated to nearly 35\%.

Another major drawback of the RSDT is that it has proved difficult to construct
confidence limits on the point estimate of abnormality due to \(\psi\) only being
approximately \(t\) distributed. To remedy this and to be able to provide an exact
rather than just approximate estimate of abnormality Crawford and colleagues started
looking into Bayesian methodology.

## Bayesian approaches {#sec23}

In the frequentist framework, parameters are treated as fixed attributes of a
population and their estimations are thought to converge to the true value
across a series of trials. In the Bayesian framework, by contrast, parameters
are treated as random variables with associated probability distributions. To
estimate a parameter distribution, we can use any prior knowledge of that
parameter to assign probabilities to possible values of the parameter, forming
what is known as a prior distribution, or simply a prior. If no prior
information is available we may want to use a non-informative prior, the most
simple of which would assign equal probabilities to all possible parameter
values. The prior is updated when new information is obtained to form what is
called a posterior distribution. The posterior probability of a hypothesis
(i.e., a specified value of the parameter) is calculated by using Bayes theorem.
If we disregard the marginal probability of the data in Bayes theorem it can be
rewritten as:
\begin{equation*}
posterior \ \propto \ likelihood \times prior
\end{equation*}

What this is saying is that the posterior density of a hypothesis is
proportional (\(\propto\)) to the likelihood of the data under that hypothesis times the
prior probability of the hypothesis.

In many cases it is not feasible to calculate the posterior analytically and
instead Monte Carlo methods are often used. These methods are
mathematical algorithms that solve numerical problems by repeated sampling from
distributions. The algorithms used differ depending on the problem at hand, but
in general they are all building on rules of drawing random numbers based on the
likelihood and the prior, and observing the distribution formed after a large
number of iterations. The peak of such a distribution is often taken as an
estimation of the parameter of interest and when using non-informative priors
this often also corresponds to the maximum likelihood estimate. Hence, using a
non-informative prior often yields estimations with frequentist properties, a
feature that is useful if hypothesis testing is desired.

This was a requirement when @crawfordComparisonSingleCase2007 and
@crawfordComparingSingleCase2011 developed Bayesian versions of the tests
described in Section \@ref(sec22), now implemented in `singcar`. The
following sections describe the procedural details of these methods which are
all based on Monte Carlo simulations. The steps presented closely follow those
outlined in the original papers with just a few slight changes of notation.


### The Bayesian test of deficit {#BTD}

The Bayesian test of deficit (BTD) allows us to obtain an estimate of \(p\) (and accompanying credible
intervals), which is the proportion of controls that would obtain a value more
extreme than the case on the variate of interest.

Assume a sample of \(n\) controls on which we measure
some value \(y\) that is normally distributed with unknown mean \(\mu\) and unknown
variance \(\sigma^2\). Let \(\overline{y}\) and \(s^2\) denote the sample mean and
sample variance respectively and \(y^*\) the case score.
Because \(\mu\) and \(\sigma^2\) are unknown, the prior distribution of \(\mu\) is
conditioned on the prior distribution of \(\sigma^2\).
A non-informative prior \(\mu | \sigma^2 \sim \mathcal{N}(0, \ \infty)\) is
assumed. For the posterior we have that
the marginal distribution of \(\frac{(n-1)s^2}{\sigma^2} \sim \chi^2_{n-1}\) and
the posterior distribution for
\(\mu|\sigma^2 \sim \mathcal{N}(\overline{y}, \sigma^2/n)\), see e.g., @gelmanBayesianDataAnalysis2013 (p. 45 and 65).
To estimate \(p\), the following steps are iterated:

1. Let \(\psi\) be a random draw from a \(\chi^2\)-distribution on \(n-1 \ df\).
   Then let \(\hat{\sigma}^2_{(i)} = \frac{(n-1)s^2}{\psi}\) be the estimation of
   \(\sigma^2\) for this iteration, hence the subscript \((i)\).
2. Let \(z\) be a random draw from a standard normal distribution.
   Then let \(\hat{\mu}_{(i)}=\overline{y}+z\sqrt{(\hat{\sigma}_{(i)}^2/n)}\) be
   the estimate of \(\mu\) for this iteration.

3.  With estimates of \(\mu\) and \(\sigma\), \(p\) is calculated conditional
    on these estimates being the "correct" \(\mu\) and \(\sigma\). Let
    \(z^*_{(i)}= \frac{y^* - \hat{\mu}_{(i)}}{\sqrt{\hat{\sigma}_{(i)}^2}}\) be the standardised
    case score then \(\hat{p}_{(i)} =\Phi\left(z^*_{(i)}\right)\) or \(\hat{p}_{(i)} = 1-\Phi\left(z^*_{(i)}\right)\)
    depending on alternative hypothesis, is the estimate of \(p\) for this
    iteration. That is the probability of drawing a
    value more extreme than \(z^*_{(i)}\) from a standard normal distribution.

Repeating these steps a large number of times will yield a distribution
of \(\hat{p}\), the mean of which is taken as the point estimate of \(p\) and used
for hypothesis testing. Multiplied by 100 it gives a point estimate of the percentage
of the control population expected to exhibit a more extreme score than the case. If repeated
e.g., 1000 times, the 25th smallest and 25th largest \(\hat{p}_{(i)}\) would give the lower and
upper boundaries of the 95\% credible interval for \(p\). Similarly,
the 25th smallest and 25th largest values of \(z^*_{(i)}\) would give the lower and
upper boundaries of the 95\% credible interval for \(Z_{CC}\) the point
estimate of which is that of Equation \@ref(eq:zcc). @crawfordComparisonSingleCase2007
show that this method yields converging results to that of the frequentist test of deficit, Equation \@ref(eq:TD).

### The Bayesian standardised difference test {#BSDT}

The Bayesian standardised difference test (BSDT) follow a similar procedure to
that of BTD. However, we now assume a sample of \(n\) controls from which we
obtain values on the variates \(Y_1\) and \(Y_2\) following a bivariate normal
distribution and representing task A and task B. Let \(\overline{y}_1\) and
\(\overline{y}_2\) denote the sample means and
\begin{equation*}
\pmb{A}=\begin{bmatrix}
s^2_{1} & s_{12} \\
s_{12} & s^2_{2} \end{bmatrix}
\end{equation*}
the sample variance-covariance matrix, where \(s_{12}\) is the sample covariance,
then \(\pmb{S} =\pmb{A}(n-1)\) is the sums of squares and cross products (SSCP) matrix. It
is assumed that the observations come from a bivariate normal distribution with
unknown mean \(\pmb{\mu}\) and unknown variance \(\Sigma\):

\begin{equation*}
\pmb{\mu} = \begin{pmatrix}
\mu_1 \\
\mu_2 \end{pmatrix} \ \text{and} \ \Sigma=\begin{bmatrix}
\sigma^2_{1} & \sigma_{12} \\
\sigma_{12} & \sigma^2_{2} \end{bmatrix}
\end{equation*}

Let the case scores be denoted \(y_1^*\) and \(y_2^*\). Just as for the frequentist
dissociation tests we want to estimate the proportion \(p\) of the control
population that would exhibit a greater difference \(Y_1-Y_2\) than the case's
\(y_1^*-y_2^*\).

The multivariate generalization of the \(\chi^2\)-distribution is the Wishart
distribution. That is, a Wishart distribution of 1 dimension is a
\(\chi^2\)-distribution on \(n\) degrees of freedom. One can view the Wishart
distribution parameterised with \(n\) degrees of freedom and some
variance-covariance matrix as being a distribution of SSCP-matrices related to
that variance-covariance matrix. Similarly, one can view the inverse Wishart
distribution parameterised with the same degrees of freedom and a SSCP-matrix as
being a distribution of variance-covariance matrices related to that
SSCP-matrix.

The non-informative prior for \(f(\pmb\mu, \Sigma^{-1})\) chosen in
@crawfordComparisonSingleCase2007 was \(f(\pmb\mu, \Sigma^{-1}) \propto
|\Sigma|\) in favour of the perhaps more commonly used \(f(\pmb\mu, \Sigma^{-1})
\propto |\Sigma|^{(k+1)/2}\) ($k =$ number of variates) because, for unstandardised differences, it was shown
to yield identical interval estimates to the frequentist UDT. The posterior
marginal distribution of \(\Sigma^{-1}\) for this choice of prior is a Wishart
distribution with \(n\) degrees of freedom and scale matrix \(\pmb{S}^{-1}\).
The conditional distribution of \(\pmb\mu|\Sigma\) is then a multivariate normal
distribution with mean \([\overline{y}_1 \ \overline{y}_2]\) and variance
\(\Sigma/n\), see e.g., @gelmanBayesianDataAnalysis2013 (p. 72-73).
@crawfordComparisonSingleCase2007 refer to this as the "standard theory"
prior.

Using the standard theory prior produces good frequentist properties for
\(\sigma_1\) and \(\sigma_2\), but @bergerObjectivePriorsBivariate2008
have shown that the convergence to frequentist estimates is less good for
\(\rho\) and differences in means ($\mu_1 - \mu_2$). Instead, they recommend
that \(f(\pmb\mu, \Sigma^{-1}) \propto \frac{1}{\sigma_1\sigma_2(1-\rho^2)}\)
should be used as a general purpose prior for bivariate normal data.
To construct posteriors with this prior rejection sampling is used. Random
observations are drawn from an inverse Wishart on \(n-1\)
degrees of freedom and only a subset of the draws are accepted.
@crawfordComparingSingleCase2011 noticed however that this prior gave rise to too narrow
credible intervals and, with simulations, showed that if the sample size was
treated as \(n-1\) for estimations of \(\Sigma\), so that the intervals were
somewhat conservative, the frequentist properties of their estimations were
improved. @crawfordComparingSingleCase2011 recommend this and refer to it as the
"calibrated" prior. The procedure for sampling from the posterior using this
prior is given below:


1. Since the unbiased estimate of \(\Sigma\) is \(\pmb{S}/(n-1)\),
   we put \(\pmb{S}^*= (n-2)\pmb{S}/(n-1)\) to retain this property
   when reducing the sample size.

2. Generate an observation from an inverse Wishart distribution
   on \(n-2\) degrees of freedom and with scale matrix \(\pmb{S}^*\).
   Denote this:
   \[
   \hat{\Sigma} = \begin{bmatrix}
   \hat{\sigma}^2_{1} & \hat{\sigma}_{12} \\
   \hat{\sigma}_{12} & \hat{\sigma}^2_{2} \end{bmatrix}
   \]
   and put
   \[
   \hat{\rho}= \frac{\hat{\sigma}_{12}}{\sqrt{\hat{\sigma}^2_{1}\hat{\sigma}^2_{2}}}
   \]

3. Generate a value \(u\) from a uniform distribution such that \(u \sim U(0, 1)\).
   If \(u^2 \leq 1-\hat{\rho}^2\) we accept \(\hat\Sigma\) as the estimation of
   \(\Sigma\) for this iteration and denote it
   \begin{equation*}
   \hat{\Sigma}_{(i)}=\begin{bmatrix}
   \hat{\sigma}^2_{1(i)} & \hat{\sigma}_{12(i)} \\
   \hat{\sigma}_{12(i)} & \hat{\sigma}^2_{2(i)} \end{bmatrix}
   \end{equation*}
   otherwise
   we iterate the procedure until we have an accepted \(\hat{\Sigma}\).


Bayesians might argue that good frequentist properties are not necessary when
conducting Bayesian analyses. If so, one can use the "standard theory" prior
by generating a random observation from an inverse Wishart distribution on \(n\)
degrees of freedom and scale matrix \(\pmb{S}\) and denote it
\(\hat{\Sigma}_{(i)}\). With an estimate of \(\Sigma\) (regardless of the prior
chosen), follow the steps below to obtain an estimate of \(p\), that is the percentage
of the control population expected to show a more extreme score than the case.


1.  Let \(z_{r1}\) and \(z_{r2}\) be two random draws from a standard normal distribution.
    Perform Cholesky decomposition on \(\hat{\Sigma}_{(i)}\), that is finding the lower triangular matrix \(\pmb{T}\)
    such that \(\pmb{T}\pmb{T'}=\hat{\Sigma}_{(i)}\). Then
    \begin{equation*}
    \pmb{\hat{\mu}}_{(i)} = \begin{pmatrix}
    \hat{\mu}_{1(i)} \\
    \hat{\mu}_{2(i)} \end{pmatrix} = \begin{pmatrix}
    \overline{y}_1 \\
    \overline{y}_2 \end{pmatrix}+ \pmb{T} \begin{pmatrix}
    z_{r1} \\
    z_{r2} \end{pmatrix} / \sqrt{n}
    \end{equation*}
    is the estimation of \(\pmb{\mu}\) for this iteration.

2.  With estimations of \(\pmb{\mu}\) and \(\Sigma\) we can calculate \(p\), given that
    they are the the correct values of \(\pmb{\mu}\) and \(\Sigma\). If an unstandardised
    test is required put:
    \begin{equation*}
    z_{(i)}^* = \frac{(y_1^* - \hat{\mu}_{1(i)}) - (y^*_2 - \hat{\mu}_{2(i)})}
    {\sqrt{\hat{\sigma}^2_{1(i)}+\hat{\sigma}^2_{2(i)}-2\hat{\sigma}_{12(i)}}}
    \end{equation*}
    If a standardised test is required\footnote{In the typical application a standardised test will most often be
      required because different cognitive functions frequently are measured on different operationalised scales.}, put:
    \begin{equation*}
    z_{1(i)} = \frac{y_1^* - \hat{\mu}_{1(i)}}{\sqrt{\hat{\sigma}^2_{1(i)}}}, \ z_{2(i)} = \frac{y_2^* -
    \hat{\mu}_{2(i)}}{\sqrt{\hat{\sigma}^2_{2(i)}}}, \ \hat{\rho}_{(i)} = \frac{\hat{\sigma}_{12(i)}}{\sqrt{\hat{\sigma}_{1(i)}\hat{\sigma}_{2(i)}}}
    \end{equation*}
    and
    \begin{equation*}
    z^*_{(i)} = \frac{z_{1(i)} - z_{2(i)}}{\sqrt{2-2\hat{\rho}_{(i)}}}
    \end{equation*}

3.  Let \(\hat{p}_{(i)}\) be the tail area of a standard normal distribution
    less or greater than \(z^*_{(i)}\) (depending on alternative hypothesis).
    \(\hat{p}_{(i)}\) is then the estimate of \(p\) for this iteration.


Repeating these
steps a large number of times will yield a distribution of \(\hat{p}\), the
mean of which is taken as the point estimate of \(p\) and used for hypothesis testing.
Multiplied by 100 it gives the point estimate of the percentage of controls expected to exhibit
a more extreme task difference than the case. If repeated e.g., 1000 times, the
25th smallest and 25th largest \(\hat{p}_i\) is the lower and upper boundaries of the 95\%
credible interval for \(p\). Similarly,
the 25th smallest and 25th largest values of \(z^*_{(i)}\) gives the lower and
upper boundaries of the 95\% credible interval for \(Z_{DCC}\), the point
estimate of which is that of Equation \@ref(eq:PJ).

### Bayesian tests allowing for covariates {#cov}

@crawfordComparingSingleCase2011 extended the previously described Bayesian tests by
using Bayesian regression techniques that allowed the case's abnormality to be assessed
in the presence of covariates. These tests thus allow you to compare
the case's score on the task of interest conditioned on the results of the
controls having the same score as the case on the covariate(s). If a case has 15
years of education, his/her score on the task would be compared to the controls
with equal length of education.

In addition to improving the precision of these tests, by accounting for
spurious noise, this also facilitates the collection of larger control samples,
because it reduces the need to very closely match control samples to a single
case. Moreover, it means that the same control sample may be used for the
evaluation of multiple single-cases, since the matching on relevant covariates
of interest can be achieved statistically. One degree of freedom is lost for
each covariate included, so as a rule of thumb they should be included only if
they have a sample correlation with the variate(s) of interest stronger than
\(0.3\) [@crawfordComparingSingleCase2011]. Of course, the control sample should also ideally
bracket the case on the covariates, in order to avoid extrapolating "out
of sample". As for the earlier tests, the assumption of the variates of interest
are still that they follow a normal or bivariate normal distribution. However,
no assumptions are made about the distribution of the covariates.

Suppose we have data from a sample of size \(n\), with scores on \(m\) covariates
and \(k = 1, 2\) variates of interest, depending on whether we are testing for a
deficit or a discrepancy. We denote the covariates
\(\boldsymbol{X} = (X_1, \dots, X_m)\).

From these values we wish to estimate
\[
\boldsymbol{B} = \begin{bmatrix}
\boldsymbol{\beta}_1 & \cdots & \boldsymbol{\beta}_k
\end{bmatrix}
\]
Where \(\boldsymbol{\beta}_i\) is a vector of length \(m+1\) containing regression
coefficients for each covariate on the \(i\)th variate of interest, the first element
in \(\boldsymbol{\beta}_i\) being the intercept. We also
wish to estimate
\[
\Sigma \ | \ \boldsymbol{X} =
\begin{bmatrix}
\sigma^2_1 \ | \ \boldsymbol{X} & \rho\sigma_1\sigma_2 \ | \ \boldsymbol{X} \\
\rho\sigma_1\sigma_2 \ | \ \boldsymbol{X} & \sigma^2_2 \ | \ \boldsymbol{X}
\end{bmatrix}
\]
That is, the covariance matrix of the variates of interest conditioned upon the
covariates, for task \(Y_1\) and task \(Y_2\) i.e., when \(k=2\). If we wish to test for
a deficit, that is \(k = 1\) then \(\Sigma \ | \ \boldsymbol{X}\) is a \(1\times1\)
matrix containing the conditional variance of the variate of interest.
We assume \(\Sigma\) not to vary with the covariates.

The procedure will be outlined in the general case, that is for \(k = 1\)
or \(k = 2\). The main difference between testing for abnormality on a single
variate and testing for abnormal discrepancy between two variates is the recommended
specification of the prior.

Let \(\boldsymbol{X}\) be the \(n \times m+1\) design matrix on which we regress
\(\boldsymbol{Y}\), the \(n \times k\) response matrix.
\[
\boldsymbol{X} =
\begin{bmatrix}
1 & x_{11} & \cdots & x_{1m} \\
\vdots & \vdots & \ddots & \vdots \\
1 & x_{n1} & \cdots & x_{nm}
\end{bmatrix},\ \
\boldsymbol{Y} =
\begin{bmatrix}
 y_{11} & \cdots & y_{1k} \\
 \vdots & \ddots & \vdots \\
 y_{n1} & \cdots & y_{nk}
\end{bmatrix}
\]
The data estimates of \(\boldsymbol{B}\) and \(\Sigma\) are then
\[
\boldsymbol{B}^* = (\boldsymbol{X}'\boldsymbol{X})^{-1}\boldsymbol{X}'\boldsymbol{Y} \ \text{and} \
\Sigma^* = \frac{1}{n-m-1}(\boldsymbol{Y}-\boldsymbol{X}\boldsymbol{B}^*)'(\boldsymbol{Y}-\boldsymbol{X}\boldsymbol{B}^*)
\]
For the "standard theory" prior the posterior of \(\Sigma\) is an inverse
Wishart distribution with \(df = n - m\) when \(k=2\) and \(df=n-m-1\) when \(k
= 1\) [@crawfordComparingSingleCase2011,
@tiaoBayesianEstimationMultivariate1964] and scale matrix \((n-m-1)\Sigma^*\).
For each iteration we generate an estimate of \(\Sigma\) from this distribution.
For the "calibrated" prior (only used for \(k=2\)) the steps described in
\@ref(BSDT) are followed and we use \(\pmb{S}^* = (n-m-2)\pmb{S}/(n-m-1)\) as the
scale matrix when we sample from an inverse Wishart distribution on \((n-m-2)\)
degrees of freedom, where \(\pmb{S} = (n-m-1)\Sigma^*\). The generated
observation from either method is denoted \(\hat{\Sigma}_{(i)}\).
We have

\begin{equation}
\hat{\Sigma}_{(i)}=[\hat{s}^2_{(i)}] \ \text{when} \ k=1 \ \text{and} \
\hat{\Sigma}_{(i)} =
\begin{bmatrix}
\hat{s}^2_{1(i)} & \hat{s}_{12(i)}  \\
\hat{s}_{12(i)}  & \hat{s}^2_{2(i)}
\end{bmatrix} \ \text{when} \ k=2 \
(\#eq:sigma)
\label{eq:sigma}
\end{equation}

\(\boldsymbol{B}^*\) will be an \((m+1) \times k\) matrix. We turn this into a \(k(m + 1) \times 1\)
vector by concatenating the colums of regression coefficients in \(\boldsymbol{B}^*\), such
that \(\boldsymbol{B}^*_{\text{vec}}=(\boldsymbol{\beta}^{*'}_1, ...,\boldsymbol{\beta}^{*'}_k)'\).
Then this vector is the data estimate for
\(\boldsymbol{B}_{\text{vec}}=(\boldsymbol{\beta}'_1, ...,\boldsymbol{\beta}'_k)'\).
Take the kronecker product \(\boldsymbol{\hat\Sigma}_{(i)} \otimes (\boldsymbol{X'X})^{-1}\) and
denote this \(\boldsymbol{\Lambda}_{(i)}\). Then the posterior of \(\boldsymbol{B}_{\text{vec}}|\boldsymbol{\hat\Sigma}_{(i)}\)
is a multivariate normal distribution with mean vector \(\boldsymbol{B}^*_{\text{vec}}\)
and variance-covariance matrix \(\boldsymbol{\Lambda}_{(i)}\). Draw a random value
from this distribution such that
\(\hat{\boldsymbol{B}}^*_{\text{vec(i)}} =(\hat{\boldsymbol{\beta}}'_{1(i)}, ...,\hat{\boldsymbol{\beta}}'_{k(i)})'\)
and \(\hat{\boldsymbol{B}}_{(i)} = (\hat{\boldsymbol{\beta}}_{1(i)}, ...,\hat{\boldsymbol{\beta}}_{k(i)})\),
where each \(\hat{\boldsymbol{\beta}}_{j(i)}\) is a vector of length \(m+1\) and \(\hat{\boldsymbol{B}}_{(i)}\)
an \(m+1 \times k\) matrix.

Let \(\boldsymbol{x}^*\) be a vector of the case's values on the covariates.
The conditional values for the case on the tasks of interest
is then
\[
\hat{\boldsymbol{\mu}}_{(i)} = \hat{\boldsymbol{B}}_{(i)}\boldsymbol{x}^*
\]

We now estimate the effect size of the case using the conditional
estimations derived above. Depending on whether we are testing for a
deficit or a discrepancy (i.e., whether \(k=1\) or \(k=2\)) the calculations
of the effect sizes differ. @crawfordComparingSingleCase2011 call these effect
sizes \(Z_{CCC}\) and \(Z_{DCCC}\) for a deficit and a discrepancy respectively.
They are similar to \(Z_{CC}\) (Equation \@ref(eq:zcc)) and \(Z_{DCC}\) (Equation \@ref(eq:PJ)), however,
the extra \(C\) in the subscript indicates that they are conditional on
covariates. Denote \(y^*_1\) and \(y^*_2\) as the case's scores on the two variates
of interest. Given that we want to estimate deficiency of a case on \(Y_1\), we
calculate
\begin{equation}
\hat{Z}_{CCC(i)} = \frac{y^*_1-\hat{\mu}_{(i)}}{\hat{s}^2_{(i)}}
(\#eq:zccc)
\label{eq:zccc}
\end{equation}
For \(Z_{DCCC}\) we have two conditional means which will be denoted
\(\hat{\mu}_{1(i)}\) and \(\hat{\mu}_{2(i)}\) for \(Y_1\) and \(Y_2\)
respectively. We then calculate
\begin{equation}
\hat{Z}_{DCCC(i)} = \frac{\frac{y^*_1-\hat{\mu}_{1(i)}}{\hat{s}_{1(i)}}-\frac{y^*_2-\hat{\mu}_{2(i)}}
{\hat{s}_{2(i)}}}{\sqrt{2-2\hat{\rho}_{12(i)}}}
(\#eq:zdccc)
\label{eq:zdccc}
\end{equation}
Where \(\hat{s}_{1(i)}\) and \(\hat{s}_{2(i)}\) are conditional
standard deviations obtained from \(\hat{\Sigma}_{(i)}\) in Equation \@ref(eq:sigma).
The conditional correlation between \(Y_1\) and \(Y_2\) in the denominator is given by
\[\hat{\rho}_{12(i)} = \frac{\hat{s}_{12(i)}}{\sqrt{\hat{s}^2_{1(i)}\hat{s}^2_{2(i)}}}\]

We then find the tail area under the standard normal distribution that is less or greater than
\(\hat{Z}_{CCC(i)}\) or \(\hat{Z}_{DCCC(i)}\) depending on the problem at hand and the alternative
hypothesis specified. Denote the value obtained \(\hat{p}_{(i)}\) which then
is an estimate of \(p\). If testing for a deficit we would have:
\[
\hat{p}_{(i)}= \Phi(\hat{Z}_{CCC(i)})
\]

Iterating these steps a large number of times will yield a distribution of
\(\hat{p}\), the mean of which is taken as the point estimate of \(p\). This
estimate is used for significance testing and if multiplied by 100 it will
give an estimate of the percentage of controls expected to exhibit a more
extreme score than the case. If repeated e.g., 1000 times, the
25th smallest and 25th largest \(\hat{p}_{(i)}\) is the lower and upper boundaries of
the 95\% credible interval for \(p\).

To obtain a point estimate of \(Z_{CCC}\) and \(Z_{DCCC}\) we use Equation \@ref(eq:zccc)
and \@ref(eq:zdccc), but use the conditional means, standard deviations
and correlation calculated directly from the control sample.
The \(1-\alpha\) credible intervals for these
effect sizes are given by the \(\alpha/2\) and \(1-\alpha/2\) quantiles
of the sampled distribution of \(Z_{CCC}\) or \(Z_{DCCC}\). That is, just as
for \(p\), if repeated e.g., 1000 times, the
25th smallest and 25th largest value of \(\hat{Z}_{CCC(i)}\) or \(\hat{Z}_{DCCC(i)}\)
is the lower and upper boundaries of
the 95\% credible interval for \(Z_{CCC}\) or \(Z_{DCCC}\).


# The `singcar` package {#section3}

## Functions for testing abnormality

| Name                                                  | Function    | Source                                      | Reference              |
| -----------                                           | ----------- | -----------                                 | -----------            |
| Test of deficit                                       | `TD()`      | @crawfordComparingIndividualTest1998        | Equation \@ref(eq:TD)  |
| Bayesian test of deficit                              | `BTD()`     | @crawfordComparisonSingleCase2007           | Section \@ref(BTD)     |
| Bayesian test of deficit with covariates              | `BTD_cov()` | @crawfordComparingSingleCase2011            | Section \@ref(cov)     |
| Unstandardised difference test                        | `UDT()`     | @crawfordTestingSuspectedImpairments2005    | Equation \@ref(eq:UDT) |
| Revised standardised difference test                  | `RSDT()`    | @crawfordTestingSuspectedImpairments2005    | Equation \@ref(eq:RSDT)|
| Bayesian standardised difference test                 | `BSDT()`    | @crawfordComparisonSingleCase2007           | Section \@ref(BSDT)    |
| Bayesian standardised difference test with covariates | `BSDT_cov()`| @crawfordComparingSingleCase2011            | Section \@ref(cov)     |

Table: Main functions for significance testing of abnormality when using small control samples and principal relevant reference in prior literature.


To facilitate meta-analyses of studies using case-control comparisons,
all functions in the table above can take both summary statistics as
input as well as raw data. The output will be a list set to class `"htest"`,
for which the generic function `print` have a method. To use the package,
begin by installing and loading it:

```{r loadpack, eval=FALSE}
install.packages("singcar")
library("singcar")
```

## Example from neuropsychology

The package comes with the dataset `size_weight_illusion`, a
neuropsychological dataset from an investigation of the size-weight illusion in
DF, a patient with visual form agnosia following bilateral lesions to the
lateral occipital complex [@hassanSizeweightIllusionVisual2020]. The size-weight illusion
is a perceptual phenomenon in which smaller objects are perceived as heavier
during lifting than larger objects of equal weight [@buckinghamGettingGripHeaviness2014].
The illusion implies that sensory cues about object size affect the perception
of weight. It has been suggested that patient DF does not experience this
illusion in the same way as the healthy population when only visual cues about
object size are available [@dijkermanVisuomotorPerformancePatient2004a]. In contrast, when kinaesthetic (tactile) cues are
provided it is suggested that DF's experience of the illusion is unaffected
by her brain damage. In other words, patient DF was expected to have a deficit
in the visual size-weight illusion and exhibit an abnormally large discrepancy
(dissociation) between visual and kinaesthetic size-weight illusions. The
dataset consists of data from patient DF and 28 control participants, with the
variables sex, age and visual as well as kinaesthetic size-weight illusion. The
measure of the size-weight illusion is a scaled measure expressing the number of
grams weight difference perceived per cubic cm of volume change. Below follows
examples of how to analyse this dataset using the tests provided in
`singcar`.

```{r datahead}
head(size_weight_illusion)
```


### Testing for a deficit

The simplest way to test for an abnormality on a single variate is
to use the frequentist test of deficit. Start by extracting patient (patient DF
is the first observation) and control data from the relevant variate, in this
case the visual size-weight illusion:
```{r VSWI}
PAT_VSWI <- size_weight_illusion[1, "V_SWI"]
CON_VSWI <- size_weight_illusion[-1, "V_SWI"] 
```
Using the function `TD()` we then apply the formula in Equation \@ref(eq:TD).
The argument `conf_int_spec` specifies how fine grained the search
algorithm for the confidence interval should be. The arguments `sd`
and `sample_size` can be given if the test should be based on summary
statistics rather than raw data, the `controls` argument should then be the
mean of the control sample.
```{r TD}
TD(case = PAT_VSWI,
   controls = CON_VSWI,
   sd = NULL,
   sample_size = NULL,
   alternative = "less",
   conf_int = TRUE,
   conf_level = 0.95,
   conf_int_spec = 0.01,
   na.rm = FALSE)
```
This can similarly be tested with the Bayesian analogue which has a very similar syntax.
This test yields an output that converges to that of TD as the argument for the
number of iterations (`iter`) increase. The degrees of freedom shown in the
output below is the degrees of freedom for the \(\chi^2\) distribution from which
we sample \(\psi\), as described in Section \@ref(BTD).
```{r BTD}
set.seed(42)
BTD(case = PAT_VSWI,
    controls = CON_VSWI,
    sd = NULL,
    sample_size = NULL,
    alternative = "less",
    int_level = 0.95,
    iter = 10000,
    na.rm = FALSE)
```

If the control sample for a study is not appropriately matched to the case on
variables such as, for example, age or education level it is appropriate to use tests
that account for this by allowing for the inclusion of covariates. Including
theoretically sound covariates is often a good idea. @crawfordComparingSingleCase2011
recommends however to only include a covariate if it correlates \(\geq 0.3\) with
the variate of interest, because of the loss of degrees of freedom.

The function `BTD_cov()` allows for the inclusion of covariates and therefore to
assess the patient on the task of interest by essentially comparing him/her to
the controls with the same score on the covariate. Even though the correlation
between age and visual size-weight illusion is \(< 0.3\) it is included here as
a coviariate for demonstrative purposes. Start again by extracting the
scores on the covariate for the patient and for the control participants.
```{r covariate}
PAT_age <- size_weight_illusion[1, "YRS"] 
CON_age <- size_weight_illusion[-1, "YRS"]
```
Since `BTD_cov()` is somewhat computationally intense, the number of
iterations has been reduced in the example compared to `BTD()`. For actual analysis the number
of iterations should be based on required precision. It should be noted that
there is no restriction on the number of covariates used as long as \(n > m+1\). If more than one
covariate is used the case scores should be given as a vector of values
and the control scores should be given as a data frame or matrix with the
same number of columns as the number of values in the covariate vector of the case.

If summary statistics are used instead of raw data, the argument `use_sumstats`
must be set to `TRUE` and the correlation matrix for the covariates
and variate of interest must be given as well as the sample size. In addition,
the `control_covar` argument must be supplied as an \(m \times 2\)
matrix or data frame giving the mean of the covariate(s) in the first column
and the standard deviation in the second. The degrees of freedom shown in the
output below is the degrees of freedom for the inverse Wishart distribution from which
we sample \(\psi\), as described in Section \@ref(cov).
```{r BTD-cov}
BTD_cov(case_task = PAT_VSWI,
        case_covar = PAT_age,
        control_task = CON_VSWI,
        control_covar = CON_age,
        alternative = "less",
        int_level = 0.95,
        iter = 1000,
        use_sumstats = FALSE,
        cor_mat = NULL,
        sample_size = NULL)
```

### Testing for a dissociation

For assessing abnormal discrepancy between two variates the simplest
function to use is the unstandardised difference test, Equation \@ref(eq:UDT).
This test is in `singcar`called by the function `UDT()`. However, one
should use this only if the variates are known to come from equivalent distributions.
Otherwise, tests that can evaluate standardised scores
without inflating Type I errors should be used. In the frequentist framework the
appropriate test for this is the RSDT, Equation \@ref(eq:RSDT). In this example we wish
to estimate and test the abnormality of the discrepancy between visual
and kinaesthetic size-weight illusion in patient DF. That is, we want to compare the difference
between the variates exhibited by the patient and the distribution of differences
in the healthy control sample. Again, start by extracting
the patient and control scores for the second variate of interest.
```{r KSWI}
PAT_KSWI <- size_weight_illusion[1, "K_SWI"] 
CON_KSWI <- size_weight_illusion[-1, "K_SWI"] 
```
Using the function `RSDT()` we then apply the formula in
Equation \@ref(eq:RSDT). This test is most often used two-sided due to the fact the the
sign of the discrepancy solely depends on the order of the input.
This function does, however, not provide any confidence intervals. The syntax of
`UDT()` is very similar to that of `RSDT()`, the main difference being
options for confidence intervals as shown for `TD()`. If summary statistics
are used then the additional argument `r_ab`, which is the sample
correlation, must be set as well as the sample size, standard deviation and mean
for both variates.
```{r RSDT}
RSDT(case_a = PAT_VSWI,
     case_b = PAT_KSWI,
     controls_a = CON_VSWI,
     controls_b = CON_KSWI,
     sd_a = NULL,
     sd_b = NULL,
     sample_size = NULL,
     r_ab = NULL,
     alternative = "two.sided",
     na.rm = FALSE)
```

The Bayesian analogue of this test is recommended over the RSDT
 because it keeps a better control
over Type I errors when the case exhibits extreme deficits on both variates but
no discrepancy between them [@crawfordComparisonSingleCase2007].
The syntax of the two functions is similar, but the `BSDT()` comes with more
optional arguments. For example, one has the option of applying this test
without standardising the variates by setting the argument `unstandardised`
to `TRUE`. The output then converges to that of the frequentist UDT.
Furthermore, one can choose between priors. Setting the
argument `calibrated` to `FALSE` specifies the use of the "standard
theory" prior. If left to the default (`TRUE`) an accept-reject algorithm
is deployed for each simulation of \(\Sigma\), as described in Section
\@ref(BSDT). This default prior has been shown to have better frequentist
properties for estimating \(\rho\) and differences between means in a bivariate
normal distributions [@bergerObjectivePriorsBivariate2008] and is therefore
the default and recommended choice.
```{r BSDT}
BSDT(case_a = PAT_VSWI,
     case_b = PAT_KSWI,
     controls_a = CON_VSWI,
     controls_b = CON_KSWI,
     sd_a = NULL,
     sd_b = NULL,
     sample_size = NULL,
     r_ab = NULL,
     alternative = "two.sided",
     int_level = 0.95,
     iter = 10000,
     unstandardised = FALSE,
     calibrated = TRUE,
     na.rm = FALSE)
```

If analysing discrepancies between variates in the presence of covariates the syntax is slightly
different, requiring that one specifies the case's scores on the variates of
interest as a vector and the control scores as a data frame or matrix.
One has the option to choose between the "calibrated" and "standard theory"
prior, just as for `BSDT()`, where `calibrated = TRUE` is the
recommended and default behaviour. If using summary statistics
`use_sumstats` must be set to `TRUE` and the summary input should be
supplied to the arguments `control_tasks`/`control_covar` as data
frames or matrices with the means of each variable represented by the first
column and the standard deviations by the second. The case's scores should
be supplied as vectors.

```{r BSDT_cov}
BSDT_cov(case_tasks = c(PAT_VSWI, PAT_KSWI),
         case_covar = PAT_age,
         control_tasks = cbind(CON_VSWI, CON_KSWI),
         control_covar = CON_age,
         alternative = "two.sided",
         int_level = 0.95,
         calibrated = TRUE,
         iter = 1000,
         use_sumstats = FALSE,
         cor_mat = NULL,
         sample_size = NULL)
```

### Power calculators

A further capacity of `singcar`is that it can be used to calculate
statistical power of the tests. The notion of power when comparing cases to
control samples have been somewhat overlooked for this class of statistical tests.
In recent work [@mcintoshPowerCalculationsSinglecase2020], we argued that,
even though power is inherently limited in this paradigm, a priori calculations are
still useful for study design and interpretation in neuropsychological and other applications.
Calculating power for the test of deficit is similar to calculating power for
any \(t\) test and can be done analytically.
\begin{equation} power = 1 - \beta =
T_{n-1}\left(t_{\alpha, \ n-1} \Bigg\rvert \frac{x^* - \overline{x}}{\sigma
\sqrt{\frac{n+1}{n}}}\right) 
(\#eq:TDpower)
\label{eq:TDpower} 
\end{equation}
Where \(T_{n-1}(.\rvert \theta)\) is the cumulative distribution function for the
non-central \(t\) distribution with \(n-1\) degrees of freedom and non-centrality
parameter \(\frac{y^* - \overline{y}}{\sigma \sqrt{\frac{n+1}{n}}}\) (i.e., TD, Equation
\@ref(eq:TD) and \(t_{\alpha, \ n-1}\) is the \(\alpha\) quantile of the \emph{central}
\(t\) distribution on \(n-1\) degrees of freedom (note that this is for a one-sided
test). For the unstandardised difference test power is calculated in an
analogous way by putting Equation \@ref(eq:UDT) as the non-centrality parameter. Deriving
power for the other functions in an analytic manner is however not possible (the
RSDT is only approximately \(t\) distributed) and a Monte Carlo approach has been
used for these tests. To call any power calculator in the package one simply
uses the function names with `_power` added as a suffix.

So, for example, to calculate power for the test of deficit we call `TD_power()`.
The expected case score and either sample size or desired power must
be supplied. The mean and standard deviation of the control sample
can also be specified with the arguments `mean` and `sd`.
If not, they take the default values of 0 and 1 respectively so
that the case score is interpreted as distance from the mean
in standard deviations. A conventional \(\alpha\)-level of
\(0.05\) is assumed if nothing else is supplied. The alternative
hypothesis can also be specified by the argument `alternative`:
specify `"less"` (default) or `"greater"` for a one-tailed test, specify 
`"two.sided"` for a two-tailed test.
```{r TD-power}
TD_power(case = 70,
         mean = 100,
         sd = 15,
         sample_size = 16,
         power = NULL,
         alternative = "less",
         alpha = 0.05,
         spec = 0.005)
```

`TD_power()` can also be used to calculate required sample size for a
desired level of power. For example, if we specify a desired power level of 0.6,
leave `sample_size` to its default and let the rest of the arguments
be as in the previous example, we see from the output of the function that
power will not increase more than 0.5\% for any additional participant after a sample
size of 15. That is, the algorithm stops searching
when this level of specificity has been reached and we are nearing the
asymptotic maximum power for this effect size. We can increase the specificity
by lowering the `spec` argument.
```{r TDpowersize}
TD_power(case = 70,
         mean = 100,
         sd = 15,
         sample_size = NULL,
         power = 0.6,
         alternative = "less",
         alpha = 0.05,
         spec = 0.005)
```

Power calculators for the Bayesian tests of deficit cannot calculate
required sample size. This is because they rely on simulation methods
to estimate approximate power and deploying a search algorithm to find the required sample
size for a given level of power would be computationally too intense. The syntax
is otherwise relatively similar to that of `TD_power()`. For `BTD_power()`
we have the two extra arguments `nsim` and `iter`, indicating the number
of simulations used in the power function and by `BTD()`, respectively.
```{r BTDpower}
BTD_power(case = 70,
         mean = 100,
         sd = 15,
         sample_size = 15,
         alternative = "less",
         alpha = 0.05,
         nsim = 1000,
         iter = 1000)
```

The only difference in syntax of `BTD_cov_power()` is due to the inclusion of covariates.
The variate of interest must be specified as a vector where the first element
gives the mean and the second the standard deviation in the argument `control_task`.
The covariates can be specified similarly or as an \(m \times 2\) matrix where the first
column gives the means of each covariate and the second column gives the standard
deviations. The correlation matrix of the variates must be given as well. In the example
below, power is evaluated for a test taking two covariates into account, both with a mean
of 0 and a standard deviation of 1. The correlation is specified as a \(3\times 3\)
matrix with pairwise correlations of $0.3$. The default settings include only one
covariate having a $0.3$ correlation with the variate of interest.
This function is computationally intense and hence, the number
of simulations has, for the example below, been decreased.
```{r BTDcovpower}
covars <- matrix(c(0, 1,
                   0, 1), ncol = 2, byrow = TRUE)
BTD_cov_power(case = -2,
              case_cov = c(0.2, -0.6),
              control_task = c(0, 1),
              control_covar = covars,
              cor_mat = diag(3) + 0.3 - diag(c(0.3, 0.3, 0.3)),
              sample_size = 15,
              alternative = "less",
              alpha = 0.05,
              nsim = 100,
              iter = 100)
```

For the difference tests one must supply the expected case scores on both
variates as well as sample size. The means and standard deviations of the
control sample can also be specified. If unspecified, they take on the default values of
0 and 1 respectively, so that the expected case scores are interpreted as
distances from the means in standard deviations. `RSDT_power()`,
`BSDT_power()` and `UDT_power()` additionally require an estimate of
the sample correlation between the variates of interest, `r_ab`. If this is
not specified a correlation of 0.5 is assumed by default. For
`BSDT_cov_power()` the correlation matrix between the variates of interest
and the covariates must instead be supplied (i.e., at least a \(3\times3\) matrix
where the first correlation is the correlation between the variates of
interest).

The alternative hypothesis is by default assumed to be
`"two.sided"` since the direction of the effect is dependent on the
order of the inputs, but can be specified to be `"less"` or
`"greater"` as well. The syntax is similar for all three functions but with small
differences. For `UDT_power()` one can request required sample size for a
desired power, as for `TD_power()`. Calculators for the Bayesian tests
have the extra argument `calibrated` as to be able to specify the prior.
`BSDT_cov_power()` requires input in the same format as `BTD_cov_power()`
for both `control_tasks` and `control_covar`. The two examples
below demonstrate usage for `RSDT_power()` and `BSDT_cov_power()`.
```{r RSDTpower}
RSDT_power(case_a = 70,
           case_b = 55,
           mean_a = 100,
           mean_b = 50,
           sd_a = 15,
           sd_b = 10,
           r_ab = 0.5,
           sample_size = 15,
           alternative = "two.sided",
           alpha = 0.05,
           nsim = 1000)
```

```{r BSDTcovpower}
cor_mat <- matrix(c(1,   0.5, 0.6,
                    0.5,   1, 0.3,
                    0.6, 0.3,   1), ncol = 3, byrow = TRUE)
BSDT_cov_power(case_tasks = c(70, 55),
               case_cov = 65,
               control_tasks = matrix(c(100, 15,
                                        50, 10), ncol = 2, byrow = TRUE),
               control_covar = c(50, 25),
               cor_mat = cor_mat,
               sample_size = 15,
               alternative = "two.sided",
               alpha = 0.05,
               nsim = 100,
               iter = 100,
               calibrated = TRUE)
```

# Summary

The `singcar`package for`R` has been outlined and
its main functionalities described. The package consists of methods for
estimating abnormality of a single case when compared to a population estimated
by a small sample. Methods for estimating abnormality on a single variate have been
described as well as methods for estimating abnormality of the difference
between two variates. Historically, the use of these tests has mainly been
within the field of neuropsychology, but their potential applicability is far
wider, extending to other areas of psychology and perhaps even to completely new
fields.

In neuropsychology the methods developed by John Crawford and Paul Garthwaite
are frequently used, especially the test of deficit, Equation \@ref(eq:TD), and
the revised standardised difference test, Equation \@ref(eq:RSDT). However, the
Bayesian standardised difference test, which outperforms the RSDT regarding
control of Type I errors, has not gained as much traction. Neither have the more
flexible Bayesian methods allowing for the inclusion of covariates.
It is hoped that by providing them in a documented package for a
popular language such as`R`, they will receive further uptake.

The methods described have been developed for keeping transparent control over
Type I errors, but power calculators have been implemented in the
package as well. Consideration of power can assist researchers in study design
and in setting realistic expectations for what these types of statistical
hypothesis tests can achieve [@mcintoshPowerCalculationsSinglecase2020].
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BTD.R
\name{BTD}
\alias{BTD}
\title{Bayesian Test of Deficit}
\usage{
BTD(
  case,
  controls,
  sd = NULL,
  sample_size = NULL,
  alternative = c("less", "greater", "two.sided"),
  int_level = 0.95,
  iter = 10000,
  na.rm = FALSE
)
}
\arguments{
\item{case}{Case observation, can only be a single value.}

\item{controls}{Numeric vector of observations from the control sample. If
single value, treated as mean.}

\item{sd}{If input of controls is single value, the standard
deviation of the sample must be given as well.}

\item{sample_size}{If input of controls is single value, the size of the
sample must be given as well.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"less"} (default), \code{"greater"} or
\code{"two.sided"}. You can specify just the initial letter.}

\item{int_level}{Level of confidence for credible intervals, defaults to 95\%.}

\item{iter}{Number of iterations. Set to higher for more accuracy, set to
lower for faster calculations.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
\tabular{llll}{
\code{statistic}   \tab the mean z-value over \code{iter} number of
iterations \cr\cr  \code{parameter} \tab the degrees of freedom used to
specify the posterior distribution. \cr\cr \code{p.value}    \tab the mean p-value
for all simulated Z-scores.\cr\cr \code{estimate}    \tab estimated standardised difference
(Z-CC) and point estimate of p-value. \cr\cr \code{null.value}   \tab the
value of the difference under the null hypothesis.\cr\cr \code{interval}
\tab named numerical vector containing credibility level and intervals for
both Z-CC and estimated proportion. \cr\cr \code{desc}     \tab named
numerical containing descriptive statistics: mean and standard deviations of
controls as well as sample size. \cr\cr \code{alternative}     \tab a
character string describing the alternative hypothesis.\cr\cr \code{method}
\tab a character string indicating what type of test was performed.\cr\cr
\code{data.name} \tab a character string giving the name(s) of the data as
well as summaries. }
}
\description{
Takes a single observation and compares it to a distribution estimated by a
control sample using Bayesian methodology. Calculates standardised difference
between the case score and the mean of the controls and proportions falling
above or below the case score, as well as associated credible intervals. This
approach was developed by Crawford and Garthwaite (2007) but converge to the
results of \code{\link{TD}()}, which is faster. Returns the point estimate of
the standardised difference between the case score and the mean of the
controls and the point estimate of the p-value (i.e. the percentage of the
population that would be expected to obtain a lower or higher score,
depending on the alternative hypothesis). This test is based on random number
generation which means that results may vary between runs. This is by design
and the reason for not using \code{set.seed()} to reproduce results inside
the function is to emphasise the randomness of the test. To get more accurate
and stable results please increase the number of iterations by increasing
\code{iter} whenever feasible.
}
\examples{
BTD(case = -2, controls = 0, sd = 1, sample_size = 20, iter = 1000)

BTD(case = size_weight_illusion[1, "V_SWI"],
    controls = size_weight_illusion[-1, "V_SWI"], alternative = "l", iter = 1000)

}
\references{
Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a
control or normative sample in neuropsychology: Development of a Bayesian
approach. \emph{Cognitive Neuropsychology, 24}(4), 343-372.
\doi{10.1080/02643290701290146}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TD.R
\name{TD}
\alias{TD}
\title{Test of Deficit}
\usage{
TD(
  case,
  controls,
  sd = NULL,
  sample_size = NULL,
  alternative = c("less", "greater", "two.sided"),
  conf_int = TRUE,
  conf_level = 0.95,
  conf_int_spec = 0.01,
  na.rm = FALSE
)
}
\arguments{
\item{case}{Case observation, can only be a single value.}

\item{controls}{Numeric vector of observations from the control sample. If
single value, treated as mean.}

\item{sd}{If input of controls is single value, the standard
deviation of the sample must be given as well.}

\item{sample_size}{If input of controls is single value, the size of the
sample must be gven as well.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"less"} (default), \code{"greater"} or
\code{"two.sided"}. You can specify just the initial letter.}

\item{conf_int}{Initiates a search algorithm for finding confidence
intervals. Defaults to \code{TRUE}, set to \code{FALSE} for faster
calculation (e.g. for simulations).}

\item{conf_level}{Level of confidence for intervals, defaults to 95\%.}

\item{conf_int_spec}{The size of iterative steps for calculating confidence
intervals. Smaller values gives more precise intervals but takes longer to
calculate. Defaults to a specificity of 0.01.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list of class \code{"htest"} containing the following components:
\tabular{llll}{
\code{statistic}   \tab the value of the t-statistic.\cr\cr  \code{parameter}
\tab the degrees of freedom for the t-statistic.\cr\cr \code{p.value}    \tab
the p-value for the test.\cr\cr \code{estimate}    \tab estimated
standardised difference (Z-CC) and point estimate of p-value. \cr\cr
\code{null.value}   \tab the value of the difference under the null
hypothesis.\cr\cr \code{interval}     \tab named numerical vector containing
level of confidence and confidence intervals for both Z-CC and p-value. \cr\cr
\code{desc}     \tab named numerical containing descriptive statistics: mean
and standard deviations of controls as well as sample size and standard error
used in the t-formula. \cr\cr \code{alternative}     \tab a character string
describing the alternative hypothesis.\cr\cr \code{method} \tab a character
string indicating what type of t-test was performed.\cr\cr \code{data.name}
\tab a character string giving the name(s) of the data as well as
summaries. }
}
\description{
Crawford and Howell's (1998) modified t-test. Takes a single observation and
compares it to a distribution estimated by a control sample. Calculates
standardised difference between the case score and the mean of the controls
and proportions falling above or below the case score, as well as associated
confidence intervals.
}
\details{
Returns the point estimate of the standardised difference
between the case score and the mean of the controls and the point estimate
of the p-value (i.e. the percentage of the population that would be
expected to obtain a lower or higher score, depending on the alternative
hypothesis).
}
\section{Note of caution}{

Calculating the confidence intervals relies on finding non-centrality
parameters for non-central t-distributions. Depending on the degrees of
freedom, the confidence level and the effect size exact accuracy from the
\code{stats::qt()} function used can not be guaranteed. However, the
approximations should be good enough for most cases.
See \url{https://stat.ethz.ch/pipermail/r-help/2008-June/164843.html}.
}

\examples{
TD(case = -2, controls = 0, sd = 1, sample_size = 20)

TD(case = size_weight_illusion[1, "V_SWI"],
   controls = size_weight_illusion[-1, "V_SWI"], alternative = "l")

}
\references{
Crawford, J. R., & Howell, D. C. (1998). Comparing an Individual's Test Score
Against Norms Derived from Small Samples. \emph{The Clinical Neuropsychologist,
12}(4), 482 - 486. \doi{10.1076/clin.12.4.482.7241}

Crawford, J. R., & Garthwaite, P. H. (2002). Investigation of the single case
in neuropsychology: Confidence limits on the abnormality of test scores and
test score differences. \emph{Neuropsychologia, 40}(8), 1196-1208.
\doi{10.1016/S0028-3932(01)00224-X}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissoc_power.R
\name{RSDT_power}
\alias{RSDT_power}
\title{Power calculator for RSDT}
\usage{
RSDT_power(
  case_a,
  case_b,
  mean_a = 0,
  mean_b = 0,
  sd_a = 1,
  sd_b = 1,
  r_ab = 0.5,
  sample_size,
  alternative = c("two.sided", "greater", "less"),
  alpha = 0.05,
  nsim = 10000
)
}
\arguments{
\item{case_a}{A single value from the expected case observation on task A.}

\item{case_b}{A single value from the expected case observation on task B.}

\item{mean_a}{The expected mean from the control sample on task A. Defaults
to 0.}

\item{mean_b}{The expected mean from the control sample on task B. Defaults
to 0.}

\item{sd_a}{The expected standard deviation from the control sample on task
A. Defaults to 1.}

\item{sd_b}{The expected standard deviation from the control sample on task
B. Defaults to 1.}

\item{r_ab}{The expected correlation between the tasks. Defaults to 0.5}

\item{sample_size}{The size of the control sample, vary this parameter to see
how the sample size affects power.}

\item{alternative}{The alternative hypothesis. A string of either "two.sided"
(default) or "one.sided".}

\item{alpha}{The specified Type I error rate. This can also be varied, with
effects on power. Defaults to 0.05.}

\item{nsim}{The number of simulations to run. Higher number gives better
accuracy, but low numbers such as 10000 or even 1000 are usually sufficient
for the purposes of this calculator.}
}
\value{
Returns a single value approximating the power of the test for the
  given parameters.
}
\description{
Calculates approximate power, given sample size, using Monte Carlo
simulation, for specified case scores, means and standard deviations for the
control sample. The means and standard deviations defaults to 0 and 1
respectively, so if no other values are given the case scores are interpreted
as deviations from the mean in standard deviations. Hence, the effect size of
the dissociation (Z-DCC) would in that case be the difference between the two
case scores.
}
\examples{
RSDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
           sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20, nsim = 1000)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UDT.R
\name{UDT}
\alias{UDT}
\title{Unstandardised Difference Test}
\usage{
UDT(
  case_a,
  case_b,
  controls_a,
  controls_b,
  sd_a = NULL,
  sd_b = NULL,
  sample_size = NULL,
  r_ab = NULL,
  alternative = c("two.sided", "greater", "less"),
  conf_int = TRUE,
  conf_level = 0.95,
  conf_int_spec = 0.01,
  na.rm = FALSE
)
}
\arguments{
\item{case_a}{Case's score on task A.}

\item{case_b}{Case's score on task B.}

\item{controls_a}{Controls' scores on task A. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task A while mean and SD for task B.}

\item{controls_b}{Controls' scores on task B. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task B while mean and SD for task A.}

\item{sd_a}{If single value for task A is given as input you must
supply the standard deviation of the sample.}

\item{sd_b}{If single value for task B is given as input you must
supply the standard deviation of the sample.}

\item{sample_size}{If A or B is given as mean and SD you must supply the
sample size. If controls_a is given as vector and controls_b as mean and
SD, sample_size must equal the number of observations in controls_a.}

\item{r_ab}{If A and/or B is given as mean and SD you must supply the
correlation between the tasks.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter. Since the direction
of the expected effect depends on which task is set as A and which is set
as B, be very careful if changing this parameter.}

\item{conf_int}{Initiates a search algorithm for finding confidence
intervals. Defaults to \code{TRUE}, set to \code{FALSE} for faster
calculation (e.g. for simulations).}

\item{conf_level}{Level of confidence for intervals, defaults to 95\%.}

\item{conf_int_spec}{The size of iterative steps for calculating confidence
intervals. Smaller values gives more precise intervals but takes longer to
calculate. Defaults to a specificity of 0.01.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab the t-statistic. \cr\cr
  \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
  \code{p.value}    \tab the p-value of the test.\cr\cr \code{estimate} \tab
  unstandardised case scores, task difference and pont estimate of proportion
  control population expected to above or below the observed task difference.
  \cr\cr \code{control.desc}   \tab named numerical with descriptive
  statistics of the control samples. \cr\cr \code{null.value}   \tab the
  value of the difference under the null hypothesis.\cr\cr
  \code{alternative}     \tab a character string describing the alternative
  hypothesis.\cr\cr \code{method} \tab a character string indicating what
  type of test was performed.\cr\cr \code{data.name} \tab a character string
  giving the name(s) of the data}
}
\description{
A test on the discrepancy between two tasks in a single case, by comparison
to the mean of discrepancies of the same two tasks in a control sample. Use
\emph{only} when the two tasks are measured on the same scale with the same
underlying distribution because no standardisation is performed on task
scores. As a rule-of-thumb, the UDT may be applicable to pairs of tasks for
which it would be sensible to perform a paired t-test within the control
group. Calculates however a standardised effect size in the same manner as
\code{\link{RSDT}()}. This is original behaviour from Crawford and Garthwaite
(2005) but might not be appropriate. So use this standardised effect size
with caution. Calculates a standardised effect size of task discrepancy as
well as a point estimate of the proportion of the control population that
would be expected to show a more extreme discrepancy and respective
confidence intervals.
}
\details{
Running  \code{UDT} is equivalent to running \code{TD} on discrepancy scores
making it possible to run unstandardised tests with covariates by applying
\code{BTD_cov} to discrepancy scores.
}
\examples{
UDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
sd_b = 1, sample_size = 20, r_ab = 0.68)

UDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
 controls_a = size_weight_illusion[-1, "V_SWI"], controls_b = size_weight_illusion[-1, "K_SWI"])

}
\references{
Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
Suspected Impairments and Dissociations in Single-Case Studies in
Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
\doi{10.1037/0894-4105.19.3.318}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSDT.R
\name{BSDT}
\alias{BSDT}
\title{Bayesian Standardised Difference Test}
\usage{
BSDT(
  case_a,
  case_b,
  controls_a,
  controls_b,
  sd_a = NULL,
  sd_b = NULL,
  sample_size = NULL,
  r_ab = NULL,
  alternative = c("two.sided", "greater", "less"),
  int_level = 0.95,
  iter = 10000,
  unstandardised = FALSE,
  calibrated = TRUE,
  na.rm = FALSE
)
}
\arguments{
\item{case_a}{Case's score on task A.}

\item{case_b}{Case's score on task B.}

\item{controls_a}{Controls' scores on task A. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task A while mean and SD for task B.}

\item{controls_b}{Controls' scores on task A. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task B while mean and SD for task A.}

\item{sd_a}{If single value for task A is given as input you must
supply the standard deviation of the sample.}

\item{sd_b}{If single value for task B is given as input you must
supply the standard deviation of the sample.}

\item{sample_size}{If A or B is given as mean and SD you must supply the
sample size. If controls_a is given as vector and controls_b as mean and
SD, sample_size must equal the number of observations in controls_a.}

\item{r_ab}{If A or B is given as mean and SD you must supply the
correlation between the tasks.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter. Since the direction
of the expected effect depends on which task is set as A and which is set
as B, be very careful if changing this parameter.}

\item{int_level}{Level of confidence for credible intervals, defaults to 95\%.}

\item{iter}{Number of iterations, defaults to 10000. Greater number gives better
estimation but takes longer to calculate.}

\item{unstandardised}{Estimate z-value based on standardised or
unstandardised task scores. Set to \code{TRUE} only if tasks are measured on the
same scale with the same underlying distribution.}

\item{calibrated}{\code{TRUE} is default. Whether or not to use the standard theory (Jeffreys) prior
distribution (if set to \code{FALSE}) or a calibrated prior examined by
Berger and Sun (2008). The sample estimation of the covariance matrix is
based on the sample size being n - 1 when the calibrated prior is used. See
Crawford et al. (2011) for further information. Calibrated prior is
recommended.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab the mean z-value over \code{iter}
  number of iterations. \cr\cr \code{parameter} \tab the degrees of freedom
  used to specify the posterior distribution. \cr\cr \code{p.value}    \tab
  the mean p-value over \code{iter} number of iterations. \cr\cr
  \code{estimate} \tab case scores expressed as z-scores on task A and B.
  Standardised effect size (Z-DCC) of task difference between case and
  controls and point estimate of the proportion of the control population
  estimated to show a more extreme task difference. \cr\cr  \code{null.value}
  \tab the value of the difference under the null hypothesis.\cr\cr
  \code{alternative}     \tab a character string describing the alternative
  hypothesis.\cr\cr \code{method} \tab a character string indicating what
  type of test was performed.\cr\cr \code{data.name} \tab a character string
  giving the name(s) of the data}
}
\description{
A test on the discrepancy between two tasks in a single case, by comparison
to the discrepancy of means in the same two tasks in a control sample. Can
take both tasks measured on the same scale with the same underlying
distribution or tasks measured on different scales by setting
\code{unstandardised} to \code{TRUE} or \code{FALSE} (default). Calculates a
standardised effects size of task discrepancy as well as a point estimate of
the proportion of the control population that would be expected to show a
more extreme discrepancy as well as relevant credible intervals. This test
is based on random number generation which means that results may vary
between runs. This is by design and the reason for not using \code{set.seed()}
to reproduce results inside the function is to emphasise the randomness of
the test. To get more accurate and stable results please increase the number
of iterations by increasing \code{iter} whenever feasible. Developed by
Crawford and Garthwaite (2007).
}
\details{
Uses random generation of inverse wishart distributions from the
CholWishart package (Geoffrey Thompson, 2019).
}
\examples{
BSDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
sd_b = 1, sample_size = 20, r_ab = 0.68, iter = 100)

BSDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
 controls_a = size_weight_illusion[-1, "V_SWI"],
 controls_b = size_weight_illusion[-1, "K_SWI"], iter = 100)

}
\references{
Berger, J. O., & Sun, D. (2008). Objective Priors for the Bivariate Normal
Model. \emph{The Annals of Statistics, 36}(2), 963-982. JSTOR.

Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a
control or normative sample in neuropsychology: Development of a Bayesian
approach. \emph{Cognitive Neuropsychology, 24}(4), 343-372.
\doi{10.1080/02643290701290146}

Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing a single
case to a control sample: Testing for neuropsychological deficits and
dissociations in the presence of covariates. \emph{Cortex, 47}(10),
1166-1178. \doi{10.1016/j.cortex.2011.02.017}

Geoffrey Thompson (2019). CholWishart: Cholesky Decomposition of the Wishart
Distribution. R package version 1.1.0.
\url{https://CRAN.R-project.org/package=CholWishart}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSDT.R
\name{RSDT}
\alias{RSDT}
\title{Revised Standardised Difference Test}
\usage{
RSDT(
  case_a,
  case_b,
  controls_a,
  controls_b,
  sd_a = NULL,
  sd_b = NULL,
  sample_size = NULL,
  r_ab = NULL,
  alternative = c("two.sided", "greater", "less"),
  na.rm = FALSE
)
}
\arguments{
\item{case_a}{Case's score on task A.}

\item{case_b}{Case's score on task B.}

\item{controls_a}{Controls' scores on task A. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task A while mean and SD for task B.}

\item{controls_b}{Controls' scores on task B. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task B while mean and SD for task A.}

\item{sd_a}{If single value for task A is given as input you must
supply the standard deviation of the sample.}

\item{sd_b}{If single value for task B is given as input you must
supply the standard deviation of the sample.}

\item{sample_size}{If A or B is given as mean and SD you must supply the
sample size. If controls_a is given as vector and controls_b as mean and
SD, sample_size must equal the number of observations in controls_a.}

\item{r_ab}{If A or B is given as mean and SD you must supply the
correlation between the tasks.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter. Since the direction
of the expected effect depends on which task is set as A and which is set
as B, be very careful if changing this parameter.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab Returns the value of a approximate
  t-statistic, however, because of the underlying equation, it cannot be
  negative. See effect direction from Z-DCC. \cr\cr \code{parameter} \tab the
  degrees of freedom for the t-statistic.\cr\cr \code{p.value}    \tab the
  p-value for the test.\cr\cr \code{estimate} \tab case scores expressed as
  z-scores on task A and Y. Standardised effect size (Z-DCC) of task
  difference between case and controls and point estimate of the proportion
  of the control population estimated to show a more extreme task
  discrepancy. \cr\cr \code{sample.size}   \tab the size of the control
  sample\cr\cr \code{null.value}   \tab the value of the discrepancy under
  the null hypothesis.\cr\cr  \code{alternative}     \tab a character string
  describing the alternative hypothesis.\cr\cr \code{method} \tab a character
  string indicating what type of test was performed.\cr\cr \code{data.name}
  \tab a character string giving the name(s) of the data}
}
\description{
A test on the discrepancy between two tasks in a single case, by comparison
to the discrepancy of means in the same two tasks in a control sample.
Standardises task scores as well as task discrepancy, so the tasks do not
need to be measured on the same scale. Calculates a standardised effect size
(Z-DCC) of task discrepancy as well as a point estimate of the proportion of
the control population that would be expected to show a more extreme
discrepancy. Developed by Crawford and Garthwaite (2005).
}
\examples{
RSDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
sd_b = 1, sample_size = 20, r_ab = 0.68)

RSDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
 controls_a = size_weight_illusion[-1, "V_SWI"], controls_b = size_weight_illusion[-1, "K_SWI"])

}
\references{
Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
Suspected Impairments and Dissociations in Single-Case Studies in
Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
\doi{10.1037/0894-4105.19.3.318}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{size_weight_illusion}
\alias{size_weight_illusion}
\title{Data from one patient and 28 controls on the size-weight illusion}
\format{
A data frame with 29 rows and 6 variables:
\describe{
  \item{GROUP}{factor with patient (SC) or control group (HC)}
  \item{PPT}{participant identifier}
  \item{SEX}{gender of partcipants}
  \item{YRS}{age of participants}
  \item{V_SWI}{SWI measure from the visual task}
  \item{K_SWI}{SWI measure from the kinaesthetic task}
}
}
\source{
\url{https://osf.io/3s2fp/?view_only=50c8af0b39ee436b85d292b0a701cc3b}
}
\usage{
size_weight_illusion
}
\description{
A dataset containing data from 28 healthy controls and one patient, DF, with
visual form agnosia (inability to perceive the form of objects) from
bilateral lesions to the lateral occipital complex. The size-weight illusion
occurs when a person underestimates the weight of a larger item when compared
to a smaller of equal weight (Charpentier, 1891). From these data, one can
assess the magnitude of the illusion for patient DF by comparison to
age-matched controls under visual and kinaesthetic cue conditions. The
measure of the size-weight illusion is a scaled measure expressing the
number of grams weight difference perceived per cubic cm of volume change
(Hassan et al, 2020).
}
\references{
Hassan, E. K., Sedda, A., Buckingham, G., & McIntosh, R. D. (2020). The
size-weight illusion in visual form agnosic patient DF. Neurocase, 1-8.
https://doi.org/10.1080/13554794.2020.1800748
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deficit_power.R
\name{BTD_power}
\alias{BTD_power}
\title{Power calculator for BTD}
\usage{
BTD_power(
  case,
  mean = 0,
  sd = 1,
  sample_size,
  alternative = c("less", "greater", "two.sided"),
  alpha = 0.05,
  nsim = 1000,
  iter = 1000
)
}
\arguments{
\item{case}{A single value from the expected case observation.}

\item{mean}{The expected mean of the control sample.}

\item{sd}{The expected standard deviation of the control sample.}

\item{sample_size}{The size of the control sample, vary this parameter to see
how the sample size affects power.}

\item{alternative}{The alternative hypothesis. A string of either "less" (default),
"greater" or "two.sided".}

\item{alpha}{The specified Type I error rate. This can also be varied, with
effects on power.}

\item{nsim}{The number of simulations for the power calculation. Defaults to
1000 due to BTD already being computationally intense.}

\item{iter}{The number of simulations used by the BTD. Defaults to 1000.}
}
\value{
Returns a single value approximating the power of the test for the
  given parameters.
}
\description{
Calculates approximate power, given sample size, using Monte Carlo simulation for the
Bayesian test of deficit for a specified case score, mean and standard
deviation for the control sample. The mean and standard deviation defaults to
0 and 1, so if no other values are given the case score is interpreted as
deviation from the mean in standard deviations.
}
\examples{
BTD_power(case = -2, mean = 0, sd = 1, sample_size = 20)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissoc_power.R
\name{BSDT_cov_power}
\alias{BSDT_cov_power}
\title{Power calculator for BSDT_cov}
\usage{
BSDT_cov_power(
  case_tasks,
  case_cov,
  control_tasks = matrix(c(0, 0, 1, 1), ncol = 2),
  control_covar = c(0, 1),
  cor_mat = diag(3) + 0.3 - diag(c(0.3, 0.3, 0.3)),
  sample_size,
  alternative = c("two.sided", "greater", "less"),
  alpha = 0.05,
  nsim = 1000,
  iter = 1000,
  calibrated = TRUE
)
}
\arguments{
\item{case_tasks}{A vector of length 2. The expected case scores from the
tasks of interest.}

\item{case_cov}{A vector containing the expected case scores on all
covariates included.}

\item{control_tasks}{A 2x2 matrix or dataframe containing the expected means
(first column) and standard deviations (second column). Defaults to two
variables with means 0 and sd = 1.}

\item{control_covar}{A px2 matrix or dataframe containing the expected means
(first column) and standard deviations (second column), p being the number
of covariates. Defaults to one covariate with mean 0 and sd = 1.}

\item{cor_mat}{A correlation matrix containing the correlations of the tasks
of interest and the coviariate(s). The first two variables are treated as
the tasks of interest. Defaults pairwise correlations between the variates of 0.3.}

\item{sample_size}{Single value giving the size of the control sample for which you wish
to calculate power.}

\item{alternative}{The alternative hypothesis. A string of either "less",
"greater" or "two.sided" (default).}

\item{alpha}{The specified Type I error rate, default is 0.05. This can be
varied, with effects on power.}

\item{nsim}{The number of simulations for the power calculation. Defaults to
1000 due to BSDT already being computationally intense. Increase for better
accuracy.}

\item{iter}{The number of simulations used by the BSDT_cov, defaults to 1000.
Increase for better accuracy.}

\item{calibrated}{Whether or not to use the standard theory (Jeffreys) prior
distribution (if set to \code{FALSE}) or a calibrated prior. See Crawford
et al. (2011) for further information. Calibrated prior is recommended.}
}
\value{
Returns a single value approximating the power of the test for the
  given parameters.
}
\description{
Computationally intense. Lower \code{iter} and/or \code{nsim} for faster but
less precise calculations. Calculates approximate power, given sample size,
using Monte Carlo simulation for BSDT with covariates
for specified (expected) case score, means and standard deviations for the
control sample on the task of interest and included covariates. The number of
covariates defaults to 1, means and standard deviations for the tasks and
covariate default to 0 and 1, so if no other values are given the case scores
is interpreted as deviation from the mean in standard deviations for both tasks
and covariates.
}
\examples{
BSDT_cov_power(c(-2, 0), case_cov = c(0, 0, 0),
control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
sample_size = 10, cor_mat = diag(5), iter = 20, nsim = 20)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissoc_power.R
\name{UDT_power}
\alias{UDT_power}
\title{Power calculator for UDT}
\usage{
UDT_power(
  case_a,
  case_b,
  mean_a = 0,
  mean_b = 0,
  sd_a = 1,
  sd_b = 1,
  r_ab = 0.5,
  sample_size = NULL,
  power = NULL,
  alternative = c("two.sided", "greater", "less"),
  alpha = 0.05,
  spec = 0.005
)
}
\arguments{
\item{case_a}{A single value from the expected case observation on task A.}

\item{case_b}{A single value from the expected case observation on task B.}

\item{mean_a}{The expected mean from the control sample on task A. Defaults
to 0.}

\item{mean_b}{The expected mean from the control sample on task B. Defaults
to 0.}

\item{sd_a}{The expected standard deviation from the control sample on task
A. Defaults to 1.}

\item{sd_b}{The expected standard deviation from the control sample on task
B. Defaults to 1.}

\item{r_ab}{The expected correlation between the tasks. Defaults to 0.5}

\item{sample_size}{The size of the control sample, vary this parameter to see
how the sample size affects power. One of sample size or power must be
specified, not both.}

\item{power}{A single value between 0 and 1 specifying desired power for
calculating necessary sample size. One of sample size or power must be
specified, not both.}

\item{alternative}{The alternative hypothesis. A string of either "two.sided"
(default) or "one.sided".}

\item{alpha}{The specified Type I error rate. This can also be varied, with
effects on power. Defaults to 0.05.}

\item{spec}{A single value between 0 and 1. If desired power is given as
input the function will utilise a search algorithm to find the sample size
needed to reach the desired power. However, if the power specified is
greater than what is actually possible to achieve the algorithm could
search forever. Hence, when power does not increase substantially for any
additional participant in the sample, the algorithm stops. By default the
algorithm stops when power does not increase more than 0.5% for any added
participant, but by varying \code{spec}, this specificity can be changed.}
}
\value{
Either a single value of the exact power, if sample size is given. Or
  a dataframe consisting of both the sample size and the exact power such
  size would yield.
}
\description{
Calculates exact power given sample size or sample size given power, using
analytical methods for the frequentist test of deficit for a specified case
scores, means and standard deviations for the control sample. The means and
standard deviations defaults to 0 and 1 respectively, so if no other values
are given, the case scores are interpreted as deviations from the mean in
standard deviations. The returned value will approximate the power for the
given parameters.
}
\examples{
UDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
          sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20)
UDT_power(case_a = -3, case_b = -1, power = 0.8)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BTD_cov.R
\name{BTD_cov}
\alias{BTD_cov}
\title{Bayesian Test of Deficit with Covariates}
\usage{
BTD_cov(
  case_task,
  case_covar,
  control_task,
  control_covar,
  alternative = c("less", "two.sided", "greater"),
  int_level = 0.95,
  iter = 10000,
  use_sumstats = FALSE,
  cor_mat = NULL,
  sample_size = NULL
)
}
\arguments{
\item{case_task}{The case score from the task of interest. Must be a single
value.}

\item{case_covar}{A vector containing the case scores on all covariates
included. Can be of any length except 0, in that case use
\code{\link{BTD}}.}

\item{control_task}{A vector containing the scores from the controls on the
task of interest. Or a vector of length 2 containing the mean and standard
deviation of the task. In that order.}

\item{control_covar}{A vector, matrix or dataframe containing the control
scores on the covariates included. If matrix or dataframe each column
represents a covariate. Or a matrix or dataframe containing summary
statistics where the first column represents the means for each covariate
and the second column represents the standard deviation.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter.}

\item{int_level}{The probability level on the Bayesian credible intervals, defaults to 95\%.}

\item{iter}{Number of iterations to be performed. Greater number gives better
estimation but takes longer to calculate. Defaults to 10000.}

\item{use_sumstats}{If set to \code{TRUE}, \code{control_tasks} and
\code{control_covar} are treated as matrices with summary statistics. Where
the first column represents the means for each variable and the second
column represents the standard deviation.}

\item{cor_mat}{A correlation matrix of all variables included. NOTE: the
first variable should be the task of interest.}

\item{sample_size}{An integer specifying the sample size of the controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab the average z-value over
  \code{iter} number of iterations. \cr\cr \code{parameter} \tab the degrees
  of freedom used to specify the posterior distribution. \cr\cr
  \code{p.value}    \tab the average p-value over \code{iter} number of
  iterations. \cr\cr \code{estimate} \tab case scores expressed as z-scores
  on task X and Y. Standardised effect size (Z-CCC) of task difference
  between case and controls and point estimate of the proportion of the
  control population estimated to show a more extreme task difference. \cr\cr
  \code{null.value}   \tab the value of the difference between tasks under
  the null hypothesis.\cr\cr \code{interval} \tab named numerical vector
  containing level of confidence and confidence intervals for both effect
  size and p-value.\cr\cr \code{desc}     \tab data frame containing means
  and standard deviations for controls as well as case scores. \cr\cr
  \code{cor.mat} \tab matrix giving the correlations between the task of
  interest and the covariates included. \cr\cr \code{sample.size} \tab number
  of controls..\cr\cr \code{alternative}     \tab a character string
  describing the alternative hypothesis.\cr\cr \code{method} \tab a character
  string indicating what type of test was performed.\cr\cr \code{data.name}
  \tab a character string giving the name(s) of the data}
}
\description{
Takes a single observation and compares it to a distribution estimated by a
control sample, while controlling for the effect of covariates, using
Bayesian methodology. This test is used when assessing a case conditioned on
some other variable, for example, assessing abnormality when controlling for
years of education or sex. Under the null hypothesis the case is an
observation from the distribution of scores from the task of interest coming
from observations having the same score as the case on the covariate(s).
Returns a significance test, point and interval estimates of difference
between the case and the mean of the controls as well as point and interval
estimates of abnormality, i.e. an estimation of the proportion of controls
that would exhibit a more extreme conditioned score. This test is based on
random number generation which means that results may vary between runs. This
is by design and the reason for not using \code{set.seed()} to reproduce
results inside the function is to emphasise the randomness of the test. To
get more accurate and stable results please increase the number of iterations
by increasing \code{iter} whenever feasible. Developed by Crawford,
Garthwaite and Ryan (2011).
}
\details{
Uses random generation of inverse wishart distributions from the
CholWishart package (Geoffrey Thompson, 2019).
}
\examples{

BTD_cov(case_task = size_weight_illusion[1, "V_SWI"],
         case_covar = size_weight_illusion[1, "YRS"],
         control_task = size_weight_illusion[-1, "V_SWI"],
         control_covar = size_weight_illusion[-1, "YRS"], iter = 100)

}
\references{
Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing
a single case to a control sample: Testing for neuropsychological deficits
and dissociations in the presence of covariates. \emph{Cortex, 47}(10),
1166-1178. \doi{10.1016/j.cortex.2011.02.017}

Geoffrey Thompson (2019). CholWishart: Cholesky Decomposition of the Wishart
Distribution. R package version 1.1.0.
\url{https://CRAN.R-project.org/package=CholWishart}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singcar.R
\docType{package}
\name{singcar}
\alias{singcar}
\title{singcar: Comparing Single Cases to Small Samples}
\description{
The aim of \pkg{singcar} is to provide and encourage usage of appropriate
statistical methods for comparing a case against a control sample. For
instance, they may commonly be done in a neuropsychological context, in which
an individual has incurred a specific brain injury and we wish to test
whether this damage has led to an impairment of some cognitive function and
whether two different functions are dissociable. For many cognitive functions
there is normed data available which the patient can be compared against
directly. However, when this is not possible a control sample estimating the
population, against which we wish to compare the patient, must be used. Both
frequentist and Bayesian methods have been developed to do this, first and
foremost by John Crawford and Paul Garthwaite (Crawford et al., 2011;
Crawford & Garthwaite, 2002, 2007, 2005; Crawford & Howell, 1998). It is
these methods that \pkg{singcar} implements. Power calculators for each respective
test are also provided. Although the canonical applications for these tests
are in Cognitive Neuropsychology or Clinical Neuropsychology, they are
potentially applicable to any circumstance in which a measure taken from a
single individual is to be compared against data from a normative sample
(i.e. a control group). It should be noted that these statistical methods
could also be applied as a general method of outlier detection in small
samples.
}
\section{singcar functions}{

\code{\link{TD}()}

\code{\link{BTD}()}

\code{\link{BTD_cov}()}

\code{\link{RSDT}()}

\code{\link{UDT}()}

\code{\link{BSDT}()}

\code{\link{BSDT_cov}()}

\code{\link{TD_power}()}

\code{\link{BTD_power}()}

\code{\link{BTD_cov_power}()}

\code{\link{RSDT_power}()}

\code{\link{UDT_power}()}

\code{\link{BSDT_power}()}

\code{\link{BSDT_cov_power}()}
}

\references{
Crawford, J., & Garthwaite, P. (2002). Investigation of the single case in
neuropsychology: Confidence limits on the abnormality of test scores and test
score differences. Neuropsychologia, 40(8), 1196-1208.
https://doi.org/10.1016/S0028-3932(01)00224-X

Crawford, J., & Garthwaite, P. (2007). Comparison of a single case to a
control or normative sample in neuropsychology: Development of a Bayesian
approach. Cognitive Neuropsychology, 24(4), 343-372.
https://doi.org/10.1080/02643290701290146

Crawford, J., & Garthwaite, P. (2005). Testing for Suspected Impairments and
Dissociations in Single-Case Studies in Neuropsychology: Evaluation of
Alternatives Using Monte Carlo Simulations and Revised Tests for
Dissociations. Neuropsychology, 19(3), 318-331.
https://doi.org/10.1037/0894-4105.19.3.318

Crawford, J., Garthwaite, P., & Ryan, K. (2011). Comparing a single case to a
control sample: Testing for neuropsychological deficits and dissociations in
the presence of covariates. Cortex, 47(10), 1166-1178.
https://doi.org/10.1016/j.cortex.2011.02.017

Crawford, J., & Howell, D. (1998). Comparing an Individual's Test Score
Against Norms Derived from Small Samples. The Clinical Neuropsychologist,
12(4), 482-486. https://doi.org/10.1076/clin.12.4.482.7241
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissoc_power.R
\name{BSDT_power}
\alias{BSDT_power}
\title{Power calculator for BSDT}
\usage{
BSDT_power(
  case_a,
  case_b,
  mean_a = 0,
  mean_b = 0,
  sd_a = 1,
  sd_b = 1,
  r_ab = 0.5,
  sample_size,
  alternative = c("two.sided", "greater", "less"),
  alpha = 0.05,
  nsim = 1000,
  iter = 1000,
  calibrated = TRUE
)
}
\arguments{
\item{case_a}{A single value from the expected case observation on task A.}

\item{case_b}{A single value from the expected case observation on task B.}

\item{mean_a}{The expected mean from the control sample on task A. Defaults
to 0.}

\item{mean_b}{The expected mean from the control sample on task B. Defaults
to 0.}

\item{sd_a}{The expected standard deviation from the control sample on task
A. Defaults to 1.}

\item{sd_b}{The expected standard deviation from the control sample on task
B. Defaults to 1.}

\item{r_ab}{The expected correlation between the tasks. Defaults to 0.5}

\item{sample_size}{The size of the control sample, vary this parameter to see
how the sample size affects power.}

\item{alternative}{The alternative hypothesis. A string of either "two.sided"
(default) or "one.sided".}

\item{alpha}{The specified Type I error rate. This can be varied, with
effects on power. Defaults to 0.05.}

\item{nsim}{The number of simulations to run. Higher number gives better
accuracy, but low numbers such as 10000 or even 1000 are usually sufficient
for the purposes of this calculator. Defaults to 1000 due to the
computationally intense \code{BSTD}.}

\item{iter}{The number simulations used by \code{BSTD}. Defaults to 1000.}

\item{calibrated}{Whether or not to use the standard theory (Jeffreys) prior
distribution (if set to \code{FALSE}) or a calibrated prior. See Crawford
et al. (2011) for further information. Calibrated prior is recommended.}
}
\value{
Returns a single value approximating the power of the test for the
  given parameters.
}
\description{
Calculates approximate power, given sample size, using Monte Carlo
simulation, for specified case scores, means and standard deviations for the
control sample. The means and standard deviations default to 0 and 1
respectively, so if no other values are given the case scores are interpreted
as deviations from the mean in standard deviations. Hence, the effect size of
the dissociation (Z-DCC) would in that case be the difference between the two
case scores. Is computationally heavy and might therefore take a few seconds.
}
\examples{
BSDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
           sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20, nsim = 100, iter = 100)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deficit_power.R
\name{TD_power}
\alias{TD_power}
\title{Power calculator for TD}
\usage{
TD_power(
  case,
  mean = 0,
  sd = 1,
  sample_size = NULL,
  power = NULL,
  alternative = c("less", "greater", "two.sided"),
  alpha = 0.05,
  spec = 0.005
)
}
\arguments{
\item{case}{A single value from the expected case observation.}

\item{mean}{The expected mean of the control sample.}

\item{sd}{The expected standard deviation of the control sample.}

\item{sample_size}{The size of the control sample, vary this parameter to see
how the sample size affects power. One of sample size or power must be
specified, not both.}

\item{power}{A single value between 0 and 1 specifying desired power for
calculating necessary sample size. One of sample size or power must be
specified, not both.}

\item{alternative}{The alternative hypothesis. A string of either "less" (default),
"greater" or "two.sided".}

\item{alpha}{The specified Type I error rate. This can also be varied, with
effects on power.}

\item{spec}{A single value between 0 and 1. If desired power is given as
input the function will utilise a search algorithm to find the sample size
needed to reach the desired power. However, if the power specified is
greater than what is actually possible to achieve the algorithm could
search forever. Hence, when power does not increase substantially for
any additional participant in the sample, the algorithm stops.
By default the algorithm stops when power does not increase more
than 0.5\% for any added participant, but by varying \code{spec},
this specificity can be changed.}
}
\value{
Either a single value of the exact power, if sample size is given. Or
  a dataframe consisting of both the sample size and the exact power such
  size would yield.
}
\description{
Calculates exact power given sample size or sample size given power, using
analytical methods for the frequentist test of deficit for a specified case
score and mean and standard deviation for the control sample. The mean and
standard deviation defaults to 0 and 1, so if no other values are given the
case score is interpreted as deviation from the mean in standard deviations.
}
\examples{
TD_power(case = -2, mean = 0, sd = 1, sample_size = 20)
TD_power(case = -2, mean = 0, sd = 1, power = 0.8)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSDT_cov.R
\name{BSDT_cov}
\alias{BSDT_cov}
\title{Bayesian Standardised Difference Test with Covariates}
\usage{
BSDT_cov(
  case_tasks,
  case_covar,
  control_tasks,
  control_covar,
  alternative = c("two.sided", "greater", "less"),
  int_level = 0.95,
  calibrated = TRUE,
  iter = 10000,
  use_sumstats = FALSE,
  cor_mat = NULL,
  sample_size = NULL
)
}
\arguments{
\item{case_tasks}{A vector of length 2. The case scores from the two tasks.}

\item{case_covar}{A vector containing the case scores on all covariates
included.}

\item{control_tasks}{A matrix or dataframe with 2 columns and n rows
containing the control scores for the two tasks. Or if \code{use_sumstats}
is set to \code{TRUE} a 2x2 matrix or dataframe
containing summary statistics where the first column represents the means
for each task and the second column represents the standard deviation.}

\item{control_covar}{A matrix or dataframe containing the control scores on
the covariates included. Or if \code{use_sumstats}
is set to \code{TRUE} a matrix or dataframe containing summary
statistics where the first column represents the means for each covariate
and the second column represents the standard deviation.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter. Since the direction
of the expected effect depends on which task is set as A and which is set
as B, be very careful if changing this parameter.}

\item{int_level}{The probability level on the Bayesian credible intervals, defaults to 95\%.}

\item{calibrated}{Whether or not to use the standard theory (Jeffreys) prior
distribution (if set to \code{FALSE}) or a calibrated prior examined by
Berger and Sun (2008). The sample estimation of the covariance matrix is based on
the sample size being n - 1 when the calibrated prior is used. See Crawford
et al. (2011) for further information. Calibrated prior is recommended.}

\item{iter}{Number of iterations to be performed. Greater number gives better
estimation but takes longer to calculate. Defaults to 10000.}

\item{use_sumstats}{If set to \code{TRUE}, \code{control_tasks} and
\code{control_covar} are treated as matrices with summary statistics. Where
the first column represents the means for each variable and the second
column represents the standard deviation.}

\item{cor_mat}{A correlation matrix of all variables included. NOTE: the two
first variables should be the tasks of interest. Only needed if \code{use_sumstats}
is set to \code{TRUE}.}

\item{sample_size}{An integer specifying the sample size of the controls.
Only needed if \code{use_sumstats}
is set to \code{TRUE}.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab the average z-value over
  \code{iter} number of iterations. \cr\cr \code{parameter} \tab the degrees
  of freedom used to specify the posterior distribution. \cr\cr
  \code{p.value}    \tab the average p-value over \code{iter} number of
  iterations. \cr\cr \code{estimate} \tab case scores expressed as z-scores
  on task A and B. Standardised effect size (Z-DCCC) of task difference
  between case and controls and point estimate of the proportion of the
  control population estimated to show a more extreme task difference. \cr\cr
  \code{null.value}   \tab the value of the difference between tasks under
  the null hypothesis.\cr\cr \code{interval} \tab named numerical vector
  containing level of confidence and confidence intervals for both effect
  size and p-value.\cr\cr \code{desc}     \tab data frame containing means
  and standard deviations for controls as well as case scores. \cr\cr
  \code{cor.mat} \tab matrix giving the correlations between the tasks of
  interest and the covariates included. \cr\cr \code{sample.size}     \tab
  number of controls. \cr\cr \code{alternative} \tab a character string
  describing the alternative hypothesis.\cr\cr \code{method} \tab a character
  string indicating what type of test was performed.\cr\cr \code{data.name}
  \tab a character string giving the name(s) of the data}
}
\description{
Takes two single observations from a case on two variables (A and B) and
compares their standardised discrepancy to the discrepancies of the variables
in a control sample, while controlling for the effects of covariates, using
Bayesian methodology. This test is used when assessing a case conditioned on
some other variable, for example, assessing abnormality of discrepancy when
controlling for years of education or sex. Under the null hypothesis the case
is an observation from the distribution of discrepancies between the tasks of
interest coming from observations having the same score as the case on the
covariate(s). Returns a significance test, point and interval estimates of
difference between the case and the mean of the controls as well as point and
interval estimates of abnormality, i.e. an estimation of the proportion of
controls that would exhibit a more extreme conditioned score. This test is
based on random number generation which means that results may vary between
runs. This is by design and the reason for not using \code{set.seed()} to
reproduce results inside the function is to emphasise the randomness of the
test. To get more accurate and stable results please increase the number of
iterations by increasing \code{iter} whenever feasible. Developed by
Crawford, Garthwaite and Ryan (2011).
}
\details{
Uses random generation of inverse wishart distributions from the
CholWishart package (Geoffrey Thompson, 2019).
}
\examples{
BSDT_cov(case_tasks = c(size_weight_illusion[1, "V_SWI"],
                        size_weight_illusion[1, "K_SWI"]),
         case_covar = size_weight_illusion[1, "YRS"],
         control_tasks = cbind(size_weight_illusion[-1, "V_SWI"],
                               size_weight_illusion[-1, "K_SWI"]),
         control_covar = size_weight_illusion[-1, "YRS"], iter = 100)

}
\references{
Berger, J. O., & Sun, D. (2008). Objective Priors for the Bivariate Normal
Model. \emph{The Annals of Statistics, 36}(2), 963-982. JSTOR.

Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing a single
case to a control sample: Testing for neuropsychological deficits and
dissociations in the presence of covariates. \emph{Cortex, 47}(10),
1166-1178. \doi{10.1016/j.cortex.2011.02.017}

#' Geoffrey Thompson (2019). CholWishart: Cholesky Decomposition of the Wishart
Distribution. R package version 1.1.0.
\url{https://CRAN.R-project.org/package=CholWishart}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deficit_power.R
\name{BTD_cov_power}
\alias{BTD_cov_power}
\title{Power calculator for BTD_cov}
\usage{
BTD_cov_power(
  case,
  case_cov,
  control_task = c(0, 1),
  control_covar = c(0, 1),
  cor_mat = diag(2) + 0.3 - diag(c(0.3, 0.3)),
  sample_size,
  alternative = c("less", "greater", "two.sided"),
  alpha = 0.05,
  nsim = 1000,
  iter = 1000
)
}
\arguments{
\item{case}{A single value from the expected case observation on the task of
interest.}

\item{case_cov}{A vector of expected case observations from covariates of
interest.}

\item{control_task}{A vector of length 2 containing the expected mean and standard
deviation of the task of interest. In that order.}

\item{control_covar}{A matrix with 2 columns containing expected means (in the 1st
column) and standard deviations (in the 2nd column) of the included
covariates.}

\item{cor_mat}{A correlation matrix containing the correlations of the
task of interest and the coviariate(s). The first variable is treated as
the task of interest. Defaults to a correlation of 0.3 between the covariate
and the variate of interest.}

\item{sample_size}{Single value of the size of the sample for which you wish
to calculate power.}

\item{alternative}{The alternative hypothesis. A string of either "less" (default),
"greater" or "two.sided".}

\item{alpha}{The specified Type I error rate. This can also be varied, with
effects on power.}

\item{nsim}{The number of simulations for the power calculation. Defaults to
1000 due to BTD_cov already being computationally intense.}

\item{iter}{The number of simulations used by the BTD_cov. Defaults to 1000.}
}
\value{
Returns a single value approximating the power of the test for the
  given parameters.
}
\description{
Computationally intense. Lower \code{iter} and/or \code{nsim} for less exact
but faster calculations. Calculates approximate power, given sample size,
using Monte Carlo simulation for the Bayesian test of deficit with covariates
for specified (expected) case score, means and standard deviations for the
control sample on the task of interest and included covariates. The number of
covariates defaults to 1, means and standard deviations for the task and
covariate defaults to 0 and 1, so if no other values are given the case score
is interpreted as deviation from the mean in standard deviations for both task
and covariate.
}
\examples{
cor_mat = matrix(c(1, 0.2, 0.3, 0.2, 1, 0.4, 0.3, 0.4, 1), ncol = 3)

BTD_cov_power(case = -2, case_cov = c(105, 30), control_task = c(0, 1),
control_covar = matrix(c(100, 40, 15, 10), ncol = 2), sample_size = 15,
cor_mat = cor_mat, iter = 20, nsim = 20)
}
