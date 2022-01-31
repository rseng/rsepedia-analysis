
<!-- README.md is generated from README.Rmd. Please edit that file -->

# presize <img src='man/figures/logo.png' align="right" height="200">

[![](https://www.r-pkg.org/badges/version/presize?color=green)](https://cran.r-project.org/package=presize)
[![](https://img.shields.io/badge/dev%20version-0.2.4.9001-blue.svg)](https://github.com/CTU-Bern/presize)
[![](https://img.shields.io/badge/shiny%20app-0.2.3-silver.svg)](https://shiny.ctu.unibe.ch/presize)
[![Actions
Status](https://github.com/CTU-Bern/presize/workflows/R-CMD-fullcheck/badge.svg)](https://github.com/CTU-Bern/presize/actions)
[![Codecov test
coverage](https://codecov.io/gh/CTU-Bern/presize/branch/master/graph/badge.svg)](https://codecov.io/gh/CTU-Bern/presize?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03118/status.svg)](https://doi.org/10.21105/joss.03118)

[Bland (2009)](https://www.bmj.com/content/339/bmj.b3985) recommended to
base study sizes on the width of the confidence interval rather the
power of a statistical test. The goal of `presize` is to provide
functions for such precision based sample size calculations. For a given
sample size, the functions will return the precision (width of the
confidence interval), and vice versa.

## Installation

`presize` can be installed from CRAN in the usual manner:

``` r
install.packages("presize")
```

You can install the development version of `presize` with:

``` r
install.packages('presize', repos = 'https://ctu-bern.r-universe.dev')
```

## Overview

presize provides functions for

| Measure                               | Function         | Methods available                                                                                                        |
| ------------------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------ |
| **Descriptive measures**              |                  |                                                                                                                          |
| Mean                                  | `prec_mean`      |                                                                                                                          |
| Proportion                            | `prec_prop`      | Wilson, Agresti-Coull, exact, Wald (see Brown, Cai, and DasGupta 2001)                                                   |
| Rate                                  | `prec_rate`      | Score, variance stabilizing, exact, Wald (see Barker 2002)                                                               |
| **Absolute differences**              |                  |                                                                                                                          |
| Mean difference                       | `prec_meandiff`  |                                                                                                                          |
| Risk difference                       | `prec_riskdiff`  | Newcombe (Newcombe 1998), Miettinen-Nurminen (Miettinen and Nurminen 1985), Agresti-Caffo (Agresti and Caffo 2000), Wald |
| **Relative differences**              |                  |                                                                                                                          |
| Odds ratio                            | `prec_or`        | Gart, Wolff, independence smoothed logit (see Fagerland, Lydersen, and Laake 2015)                                       |
| Risk ratio                            | `prec_riskratio` | Koopman (Koopman 1984), Katz (Katz et al. 1978)                                                                          |
| Rate ratio                            | `prec_rateratio` | Rothman (Rothman and Greenland 2018)                                                                                     |
| **Correlation measures**              |                  |                                                                                                                          |
| Correlation coefficient               | `prec_cor`       | Pearson, Kendall, Spearman (see Bonnett and Wright 2000)                                                                 |
| Intraclass correlation                | `prec_icc`       | Bonnett (2002)                                                                                                           |
| Limit of agreement                    | `prec_lim_agree` | Bland and Altman (1986)                                                                                                  |
| Cohen’s kappa                         | `prec_kappa`     | Rotondi and Donner (2012)                                                                                                |
| **Diagnostic measures**               |                  |                                                                                                                          |
| Sensitivity<sup>1</sup>               | `prec_sens`      | As per `prec_prop`                                                                                                       |
| Specificity<sup>1</sup>               | `prec_spec`      | As per `prec_prop`                                                                                                       |
| Area under the curve                  | `prec_auc`       | Hanley and McNeil (1982)                                                                                                 |
| Negative likelilood ratio<sup>2</sup> | `preg_neg_lr`    | Simel, Samsa, and Matchar (1991)                                                                                         |
| Positive likelilood ratio<sup>2</sup> | `preg_pos_lr`    | Simel, Samsa, and Matchar (1991)                                                                                         |
| Generic likelilood ratio              | `preg_lr`        | Simel, Samsa, and Matchar (1991)                                                                                         |

<sup>1</sup> Simple wrappers for `prec_prop`.

<sup>2</sup> Wrappers for `prec_lr` with values provided via sens and
spec

## Example

Suppose we want to estimate the proportion of hospital admissions with
diabetes. Diabetes has a prevalence of approximately 10% (Emerging Risk
Factors Collaboration et al. (2010)). We assume a slightly higher
proportion of diabetics, 15%, as diabetes is a risk factor for a wide
range of conditions. We want to estimate the prevalence of diabetes to
within 5% (plus/minus 2.5%). With `presize`, this is simple. We use the
`prec_prop` (precision of a proportion) function and pass our 15% and 5%
as arguments `p` and `conf.width`:

``` r
library(presize) # load the package
prec_prop(p = 0.15, conf.width = 0.05)
#> Warning in prec_prop(p = 0.15, conf.width = 0.05): more than one method was
#> chosen, 'wilson' will be used
#> 
#>      sample size for a proportion with Wilson confidence interval. 
#> 
#>      p      padj        n conf.width conf.level       lwr       upr
#> 1 0.15 0.1517077 783.4897       0.05       0.95 0.1267077 0.1767077
#> 
#> NOTE: padj is the adjusted proportion, from which the ci is calculated.
```

In the n column, we see that we would need to ask 784 (rounding 783.5
up) patients to achieve the desired CI width. Disappointingly, we also
know that we only have funds to collect the data from 600 patients. We
wonder if 600 patients would yield sufficient precision - we could also
accept a CI width of 6% (plus/minus 3%). In such a case, we can pass the
arguments `p` and `n`.

``` r
prec_prop(p = 0.15, n = 600)
#> Warning in prec_prop(p = 0.15, n = 600): more than one method was chosen,
#> 'wilson' will be used
#> 
#>      precision for a proportion with Wilson confidence interval. 
#> 
#>      p      padj   n conf.width conf.level       lwr       upr
#> 1 0.15 0.1522266 600 0.05713404       0.95 0.1236596 0.1807936
#> 
#> NOTE: padj is the adjusted proportion, from which the ci is calculated.
```

Now we see that with 600 patients, the CI would have a width of 5.7%. We
are happy with this and continue planning our study with those values.
All of the functions listed in Table 1 can be used similarly.

We can also look at a range of scenarios simulatenously by passing a
vector to one of the arguments, which could be used to create something
analogous to a power curve:

``` r
prec_prop(p = 0.15, n = seq(600, 800, 50))
#> Warning in prec_prop(p = 0.15, n = seq(600, 800, 50)): more than one method was
#> chosen, 'wilson' will be used
#> 
#>      precision for a proportion with Wilson confidence interval. 
#> 
#>      p      padj   n conf.width conf.level       lwr       upr
#> 1 0.15 0.1522266 600 0.05713404       0.95 0.1236596 0.1807936
#> 2 0.15 0.1520563 650 0.05489329       0.95 0.1246097 0.1795030
#> 3 0.15 0.1519102 700 0.05289705       0.95 0.1254617 0.1783588
#> 4 0.15 0.1517835 750 0.05110386       0.95 0.1262316 0.1773355
#> 5 0.15 0.1516726 800 0.04948148       0.95 0.1269319 0.1764133
#> 
#> NOTE: padj is the adjusted proportion, from which the ci is calculated.
```

## Shiny app

An online interactive version of the package is available
[here](https://shiny.ctu.unibe.ch/presize). The app can also be launched
locally via `launch_presize_app()` in RStudio.

![](man/figures/app.png)<!-- -->

## Getting help

The package website, including more details on the functions, can be
found [here](https://ctu-bern.github.io/presize/).

If you have a question, feel free to make a thread on the
[discussion](https://github.com/CTU-Bern/presize/discussions) page.

If you encounter a bug, please create an
[issue](https://github.com/CTU-Bern/presize/issues).

## Contributing

Contributions to `presize` are welcome. If you have ideas, open an
[issue](https://github.com/CTU-Bern/presize/issues) or a [discussion
thread](https://github.com/CTU-Bern/presize/discussions) on GitHub.

If you want to contribute code, please feel free to fork the repository,
make your changes and make a pull request to have them integrated into
the package. New functionality should have accompanying tests and pass
continuous integration. See also the [contributing
guidelines](https://github.com/CTU-Bern/presize/blob/master/CONTRIBUTING.md).

## Funding

`presize` was largely developed at CTU Bern, with collaboration from CTU
Basel. Funding was provided by the Swiss Clinical Trial Organisation.

![](man/figures/SCTO_Platforms.png)<!-- -->

<!-- ![](man/fig/scto_ctu_member_cmyk.jpg) -->

## Citation [![DOI](https://joss.theoj.org/papers/10.21105/joss.03118/status.svg)](https://doi.org/10.21105/joss.03118)

If you use `presize`, please cite it in your publication as:  
Haynes et al., (2021). presize: An R-package for precision-based sample
size calculation in clinical research. Journal of Open Source Software,
6(60), 3118, <https://doi.org/10.21105/joss.03118>

### Acknowledgements

The package logo was created with
[`ggplot2`](https://ggplot2.tidyverse.org/) and
[`hexSticker`](https://github.com/GuangchuangYu/hexSticker) with icons
from [Font Awesome](https://fontawesome.com/) (via the [emojifont
package](https://github.com/GuangchuangYu/emojifont)).

## References

<div id="refs" class="references">

<div id="ref-ac2000">

Agresti, A, and B Caffo. 2000. “Simple and Effective Confidence
Intervals for Proportions and Differences of Proportions Result from
Adding Two Successes and Two Failures.” *The Americal Statistician* 54
(4): 280–88. <https://doi.org/10.2307/2685779>.

</div>

<div id="ref-barker2002">

Barker, L. 2002. “A Comparison of Nine Confidence Intervals for a
Poisson Parameter When the Expected Number of Events Is ≤ 5.” *The
Americal Statistician* 56 (2): 85–89.
<https://doi.org/10.1198/000313002317572736>.

</div>

<div id="ref-ba1986">

Bland, J M, and D G Altman. 1986. “Statistical Methods for Assessing
Agreement Between Two Methods of Clinical Measurement.” *Lancet*
i(8476): 307–10. <https://doi.org/10.1016/S0140-6736(86)90837-8>.

</div>

<div id="ref-bonnett2002">

Bonnett, D G. 2002. “Sample Size Requirements for Estimating Intraclass
Correlations with Desired Precision.” *Statistics in Medicine* 21:
1331–5. <https://doi.org/10.1002/sim.1108>.

</div>

<div id="ref-bw2000">

Bonnett, D G, and T A Wright. 2000. “Sample Size Requirements for
Estimating Pearson, Kendall and Spearman Correlations.” *Psychometrika*
65: 23–28. <https://doi.org/10.1007/BF02294183>.

</div>

<div id="ref-brown2001">

Brown, L D, T T Cai, and A DasGupta. 2001. “Interval Estimation for a
Binomial Proportion.” *Statistical Science* 16 (2): 101–17.
<https://doi.org/10.1214/ss/1009213286>.

</div>

<div id="ref-diab">

Emerging Risk Factors Collaboration, N Sarwar, P Gao, S R Seshasai, R
Gobin, S Kaptoge, E Di Angelantonio, et al. 2010. “Diabetes Mellitus,
Fasting Blood Glucose Concentration, and Risk of Vascular Disease: A
Collaborative Meta-Analysis of 102 Prospective Studies.” *Lancet* 375
(9733): 2215–22. <https://doi.org/10.1016/S0140-6736(10)60484-9>.

</div>

<div id="ref-fll2015">

Fagerland, M W, S Lydersen, and P Laake. 2015. “Recommended Confidence
Intervals for Two Independent Binomial Proportions.” *Statistical
Methods in Medical Research* 24 (2): 224–54.
<https://doi.org/10.1177/0962280211415469>.

</div>

<div id="ref-hm1982">

Hanley, J A, and B J McNeil. 1982. “The Meaning and Use of the Area
Under a Receiver Operating Characteristic (Roc) Curve.” *Radiology* 148:
29–36. <https://doi.org/10.1148/radiology.143.1.7063747>.

</div>

<div id="ref-kbap1978">

Katz, D, J Baptista, S P Azen, and M C Pike. 1978. “Obtaining Confidence
Intervals for the Risk Ratio in Cohort Studies.” *Biometrics* 34:
469–74. <https://doi.org/10.2307/2530610>.

</div>

<div id="ref-koopman1984">

Koopman, P A R. 1984. “Confidence Intervals for the Ratio of Two
Binomial Proportions.” *Biometrics* 40: 513–17.
<https://doi.org/10.2307/2531551>.

</div>

<div id="ref-mn1985">

Miettinen, O, and M Nurminen. 1985. “Comparative Analysis of Two Rates.”
*Statistics in Medicine* 4: 213–26.
<https://doi.org/10.1002/sim.4780040211>.

</div>

<div id="ref-newcombe1998">

Newcombe, R G. 1998. “Interval Estimation for the Difference Between
Independent Proportions: Comparison of Eleven Methods.” *Statistics in
Medicine* 17: 873–90.
[https://doi.org/10.1002/(sici)1097-0258(19980430)17:8\<873::aid-sim779\>3.0.co;2-i](https://doi.org/10.1002/\(sici\)1097-0258\(19980430\)17:8%3C873::aid-sim779%3E3.0.co;2-i).

</div>

<div id="ref-rg2018">

Rothman, K J, and S Greenland. 2018. “Planning Study Size Based on
Precision Rather Than Power.” *Epidemiology* 29: 599–603.
<https://doi.org/10.1097/EDE.0000000000000876>.

</div>

<div id="ref-rd2012">

Rotondi, M A, and A Donner. 2012. “A Confidence Interval Approach to
Sample Size Estimation for Interobserver Agreement Studies with Multiple
Raters and Outcomes.” *Journal of Clinical Epidemiology* 65: 778–84.
[https://doi.org/10.1016/j.jclinepi.2011.10.019](https://doi.org/10.1016/j.jclinepi.2011.10.019%20).

</div>

<div id="ref-simel1991">

Simel, D L, G P Samsa, and D B Matchar. 1991. “Likelihood Ratios with
Confidence: Sample Size Estimation for Diagnostic Test Studies.”
*Journal of Clinical Epidemiology* 44 (8): 763–70.
<https://doi.org/10.1016/0895-4356(91)90128-v>.

</div>

</div>
presize 0.3.0
-----------------------------------------

* BREAKING CHANGE - prec_mean argument changed from 'mu' to 'mean'.


presize 0.2.4.9001
-----------------------------------------

* logo
* installation via universe


presize 0.2.4.9000
-----------------------------------------

* addition of citation info
* update SCTO figures

presize 0.2.3
-----------------------------------------

* minor formatting of URLs and references to pass CRAN checks

presize 0.2.2
-----------------------------------------

* bug fix in shiny app for sens and spec
* more examples, better described
* add extra input validation to all functions

presize 0.2.1
-----------------------------------------

* add reset buttons to shiny app
* correct copy/paste error on ICC page

presize 0.2.0
-----------------------------------------

* various changes to ensure that all functions support vectors for different scenarios


presize 0.1.4
-----------------------------------------

* minor changes to readme and shinyapp (corrections of references/typos)

presize 0.1.3
-----------------------------------------

* more minor changes for CRAN 

presize 0.1.2
-----------------------------------------

* minor changes requested by CRAN 

presize 0.1.1
-----------------------------------------

* first version for CRAN
* minor clarifications to options/descriptions
* minor changes to shiny app

presize 0.1.0
-----------------------------------------

* initial 'final' version

presize 0.0.1.9007
-----------------------------------------

* addition of wrappers for `prec_lr` (`prec_pos_lr`, `prec_neg_lr`) to simplify positive/negative LRs


presize 0.0.1.9006
-----------------------------------------

* addition of method for likelihood ratios `prec_lr`

presize 0.0.1.9005
-----------------------------------------

* addition of function for Cohen's kappa

* update shiny app to include kappa

* POSSIBLE BREAKING CHANGE: arguments in `prec_rateratio` renamed from `*_exp` and `*_control` to `*1` and `*2` for consistency with other functions

presize 0.0.1.9004
-----------------------------------------

* addition of shiny app and pkgdown

presize 0.0.1.9003
-----------------------------------------

* `prec_sens` and `prec_spec` allow prev and conf.width

* multiple notes allowed in print method

* add contributing guidelines

* `prec_sens_spec` removed. Confidence intervals and sample sizes are quite different to other methods. A two step approach using `prec_sens` and `prec_spec` is instead recommended


presize 0.0.1.9002
-----------------------------------------

* addition of various tests

* addition of rate ratio method
# Contributing to presize

This guide gives some points to contributing to `presize`. Contributions do not
have to entail literally modifying/adding code, bug reports are equally
important.

* Bug reports should be made via issues

* Modifications or additions to code can be made via pull requests (PR)

## Issues

Be concise, but explain what the problem is. Adding code and output is useful,
if possible.
Consider commenting the code to improve clarity.


## Pull requests

To contribute modifications or additions to code, create a PR. The
steps are broadly as follows.

1. Fork `presize` to your GitHub profile.
1. Clone it to your local machine.
1. Make your changes and commit them to your local git.
1. Push your changes to GitHub.
1. Create a PR
1. Discuss with reviewer about any requested changes until the PR can be
closed.

Please describe the motivation behind your PR and provide unit tests via
testthat. If your PR addresses an issue, add that info in the description
e.g. closes #{issue-number}

---
title: '`presize`: An R-package for precision-based sample size calculation in clinical research'
tags:
  - R software
  - clinical trials
  - sample size calculation
authors:
 - name: Alan G Haynes
   orcid: 0000-0003-1374-081X
   affiliation: "1, 2"
 - name: Armando Lenz
   orcid: 0000-0001-5888-0846
   affiliation: "1, 2"
 - name: Odile Stalder
   orcid: 0000-0002-5563-2975
   affiliation: "1, 2"
 - name: Andreas Limacher
   orcid: 0000-0002-9094-9476
   affiliation: "1, 2"
affiliations:
 - name: CTU Bern, University of Bern, Bern, Switzerland
   index: 1
 - name: Statistics and Methodology Platform of the Swiss Clinical Trial Organisation (SCTO), Bern, Switzerland
   index: 2
bibliography: paper.bib
date: 2021-02-11
---

# Background

Sample size calculation is a crucial step for planning a clinical study. A study 
that is too small leads to inconclusive results; a study that is too large is a 
waste of resources. Either case might be unethical. 
Furthermore, @bland2009 called for a focus on the width of confidence intervals 
rather than the power of the test in sample size calculations. Indeed, many 
research projects aim to estimate a quantity rather than test a hypothesis,
sample size calculation approaches for which are largely missing from other software
packages. There are many software packages for hypothesis-based sample size calculation, 
such as [Stata](https://www.stata.com/), [PASS](https://www.ncss.com/software/pass/),
[G*Power](https://www.gpower.hhu.de), including many R packages, such as 
[`pwr`](https://CRAN.R-project.org/package=pwr) [@pwr] and
[`TrialSize`](https://CRAN.R-project.org/package=TrialSize) 
[[@trialsize]; see other packages detailed on the 
[CRAN Clinical Trials taskview](https://cran.r-project.org/web/views/ClinicalTrials.html)].

# Statement of need

To the best of our knowledge, only Stata provides precision-based approaches, and 
only then for a small number of statistics.
We have, therefore, developed an R package for precision-based sample size calculation, 
`presize`, which can be used within the R-environment or a shiny application.

# Development

`presize` is programmed in the R programming language [@cran], and offers sample size
calculation for estimation-based research. We implemented the most common
measures used in descriptive research, including descriptive, absolute and relative 
differences, correlation and diagnostic measures. 
There are two approaches for each measure. Firstly, based on a given sample size,
e.g. for a retrospective data analysis, the precision of an expected measure can
be calculated. Precision is expressed as the confidence interval around the
measure. The level of confidence can be specified; the usual 95%-confidence
interval (CI) is the default. Secondly, based on a given precision (i.e. CI), the 
sample size can be calculated. This is mainly of use in the planning of prospective 
studies, where the aim is to estimate a measure of interest with enough confidence. 
The available functions in the package require common input arguments; either the 
sample size or the width of the CI, and the level of the CI. Depending on the 
function, further specific input arguments are required such as the expected area 
under the curve (AUC) for test accuracy.


For ease-of-use, we have also implemented a Shiny application which can be used
[online](https://ctu-bern.shinyapps.io/presize) or from within the R-environment.
Values of parameters can be entered using rulers and numeric fields. Based on
the given parameters, the application will either display the sample size for a
given precision, or vice-versa, the precision for a given sample size. Moreover,
the application also displays the corresponding R-code, which allows copying the
command into the R-environment for further exploration as well as reproducibility.


# Usage

`presize` is available on [CRAN](https://CRAN.R-project.org/package=presize) or 
[GitHub](https://github.com/CTU-Bern/presize) and can be installed and loaded into 
the R session using

```r
# installation:
# install.packages("presize") # CRAN
# remotes::install_github('ctu-bern/presize') # development version on GitHub
library(presize)
```
As a brief example, suppose we want to estimate the proportion of hospital admissions 
with diabetes. Diabetes has a prevalence of approximately 10% [@diab]. We assume a 
slightly higher proportion of diabetics, 15%, as diabetes is a risk factor for a 
wide range of conditions. We want to estimate the prevalence of diabetes to within 
5% (plus/minus 2.5%). With `presize`, this is simple. We use the `prec_prop` 
(precision of a proportion) function and pass our 15% and 5% as arguments `p` 
and `conf.width`:

```r
prec_prop(p = 0.15, conf.width = 0.05)

     sample size for a proportion with Wilson confidence interval. 

     p      padj        n conf.width conf.level       lwr       upr
1 0.15 0.1517077 783.4897       0.05       0.95 0.1267077 0.1767077
```

In the `n` column, we see that we would need to ask 784 (rounding 783.5 up) 
patients to achieve the desired CI width. It is also possible to calculate the 
CI width with a given number of participants, for the case that we know roughly 
how many participants we could include: 

```r
prec_prop(p = 0.15, n = 600)

     sample size for a proportion with Wilson confidence interval. 

     p      padj   n conf.width conf.level       lwr       upr
1 0.15 0.1522266 600 0.05713404       0.95 0.1236596 0.1807936
```

# Discussion

We have developed a comprehensive and easy-to-use software package for precision-based 
sample size calculation. As far as we know, `presize` is the first package that comprises 
the most common summary measures used in estimation-based clinical research.
A limitation of the package is that it does not allow calculating the probability 
of a CI, i.e. the probability that a future confidence interval 
will have at least the desired precision. The functions currently return the average 
CI width. In practice, 50% of trials will yield narrower CIs and 
50% will yield wider CIs due to sampling variation. Providing a method to specify 
the coverage probability is one possible avenue for further development.

We often observe in our consulting activity that researchers try to implement a 
hypothesis-based approach into a project that is in fact purely descriptive. Reasons 
for this might be a lack of methodological understanding, but also a lack of appropriate 
tools. To conclude, we believe that our software package will facilitate the 
appropriate use of estimation-based sample size calculation in descriptive research projects.

# Acknowledgements
The development of `presize` was enabled through financial support of the [Swiss 
Clinical Trial Organisation (SCTO)](https://www.scto.ch/en) as part of its [Statistics and Methodology Platform](https://www.scto.ch/en/network/scto-platforms/statistics-and-methodology.html).
We also wish to thank statisticians of [CTU Bern](https://www.ctu.unibe.ch/) and 
[CTU Basel](https://www.unispital-basel.ch/ueber-uns/das-universitaetsspital/leitung/direktion/klinische-forschung/) 
for suggestions and testing. We also thank [Tom Kelly](https://github.com/TomKellyGenetics) 
and [Bruce Mecum](https://github.com/amoeba) for reviewing this paper and the package, 
and [Mark Jensen](https://github.com/majensen) for providing editorial support.

# References

<!-- Thanks for submitting a pull request to presize! -->

**Summary**

<!-- Describe the motivation behind your PR - refer to issues where relevant -->

A couple of reminders:
- [ ] unit tests written
- [ ] continuous integration (github actions) passing?




---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Anything else...
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---


---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
output: github_document
bibliography: paper/paper.bib
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# presize <img src='man/figures/logo.png' align="right" height="200">


[![](https://www.r-pkg.org/badges/version/presize?color=green)](https://cran.r-project.org/package=presize) 
`r badger::badge_custom("dev version", as.character(packageVersion("presize")), "blue", "https://github.com/CTU-Bern/presize")`
`r badger::badge_custom("shiny app", readLines("https://shiny.ctu.unibe.ch/version/presize"), "silver", "https://shiny.ctu.unibe.ch/presize")`
[![Actions Status](https://github.com/CTU-Bern/presize/workflows/R-CMD-fullcheck/badge.svg)](https://github.com/CTU-Bern/presize/actions)
[![Codecov test coverage](https://codecov.io/gh/CTU-Bern/presize/branch/master/graph/badge.svg)](https://codecov.io/gh/CTU-Bern/presize?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03118/status.svg)](https://doi.org/10.21105/joss.03118)

[Bland (2009)](https://www.bmj.com/content/339/bmj.b3985) recommended to
base study sizes on the width of the confidence interval rather the power of 
a statistical test. The goal of `presize` is to provide functions for such 
precision based sample size calculations. For a given sample size, the 
functions will return the precision (width of the confidence interval), and 
vice versa.

## Installation

`presize` can be installed from CRAN in the usual manner:

```{r cran-installation, eval = FALSE}
install.packages("presize")
```


You can install the development version of `presize` with:

```{r gh-installation, eval = FALSE}
install.packages('presize', repos = 'https://ctu-bern.r-universe.dev')
```

## Overview
presize provides functions for

Measure | Function | Methods available 
-------- | ---------- | --------
**Descriptive measures** | |
Mean | `prec_mean` |
Proportion | `prec_prop` | Wilson, Agresti-Coull, exact, Wald [see @brown2001]
Rate | `prec_rate` | Score, variance stabilizing, exact, Wald [see @barker2002]
**Absolute differences** | |
Mean difference | `prec_meandiff` |
Risk difference | `prec_riskdiff` | Newcombe [@newcombe1998], Miettinen-Nurminen [@mn1985], Agresti-Caffo [@ac2000], Wald
**Relative differences** | |
Odds ratio | `prec_or` | Gart, Wolff, independence smoothed logit [see @fll2015]
Risk ratio | `prec_riskratio` | Koopman [@koopman1984], Katz [@kbap1978]
Rate ratio | `prec_rateratio` | Rothman [@rg2018]
**Correlation measures** | |
Correlation coefficient | `prec_cor` | Pearson, Kendall, Spearman [see @bw2000]
Intraclass correlation | `prec_icc` | @bonnett2002
Limit of agreement | `prec_lim_agree` | @ba1986
Cohen's kappa | `prec_kappa` | @rd2012
**Diagnostic measures** | |
Sensitivity<sup>1</sup> | `prec_sens` | As per `prec_prop`
Specificity<sup>1</sup> | `prec_spec` | As per `prec_prop`
Area under the curve | `prec_auc` | @hm1982
Negative likelilood ratio<sup>2</sup> | `preg_neg_lr` | @simel1991
Positive likelilood ratio<sup>2</sup> | `preg_pos_lr` | @simel1991
Generic likelilood ratio | `preg_lr` | @simel1991

<sup>1</sup> Simple wrappers for `prec_prop`.

<sup>2</sup> Wrappers for `prec_lr` with values provided via sens and spec

## Example

Suppose we want to estimate the proportion of hospital admissions with diabetes. 
Diabetes has a prevalence of approximately 10% (@diab). We assume a 
slightly higher proportion of diabetics, 
15%, as diabetes is a risk factor for a wide range of conditions. We want to 
estimate the prevalence of diabetes to within 5% (plus/minus 2.5%). With `presize`,
this is simple. We use the `prec_prop` (precision of a proportion) function and pass 
our 15% and 5% as arguments `p` and `conf.width`:

```{r}
library(presize) # load the package
prec_prop(p = 0.15, conf.width = 0.05)
```


In the n column, we see that we would need to ask 784 (rounding 783.5 up) patients to achieve the desired CI width. 
Disappointingly, we also know that we only have funds to collect the data from 
600 patients. 
We wonder if 600 patients would yield sufficient precision - we could 
also accept a CI width of 6% (plus/minus 3%).
In such a case, we can pass the arguments `p` and `n`.

```{r}
prec_prop(p = 0.15, n = 600)
```

Now we see that with 600 patients, the CI would have a width of 
5.7%. We are happy with this and continue planning our study with those values. 
All of the functions listed in Table 1 can be used similarly.

We can also look at a range of scenarios simulatenously by passing a vector to 
one of the arguments, which could be used to create something analogous to a 
power curve: 

```{r}
prec_prop(p = 0.15, n = seq(600, 800, 50))
```


## Shiny app

An online interactive version of the package is available [here](https://shiny.ctu.unibe.ch/presize). The app can also be launched locally via `launch_presize_app()` in RStudio.

```{r, echo=FALSE}
knitr::include_graphics("man/figures/app.png")
```

## Getting help

The package website, including more details on the functions, can be found [here](https://ctu-bern.github.io/presize/).

If you have a question, feel free to make a thread on the [discussion](https://github.com/CTU-Bern/presize/discussions) page.

If you encounter a bug, please create an [issue](https://github.com/CTU-Bern/presize/issues).

## Contributing

Contributions to `presize` are welcome. If you have ideas, open an [issue](https://github.com/CTU-Bern/presize/issues) or a [discussion thread](https://github.com/CTU-Bern/presize/discussions) on GitHub. 

If you want to contribute code, please feel free to fork the repository, make your changes and make a pull request to have them integrated into the package.  New functionality should have accompanying tests and pass continuous integration. See also the [contributing guidelines](https://github.com/CTU-Bern/presize/blob/master/CONTRIBUTING.md).

## Funding

`presize` was largely developed at CTU Bern, with collaboration from CTU Basel. Funding was provided by the Swiss Clinical Trial Organisation.

```{r, echo=FALSE}
knitr::include_graphics("man/figures/SCTO_Platforms.png")
```

<!-- ![](man/fig/scto_ctu_member_cmyk.jpg) -->

## Citation [![DOI](https://joss.theoj.org/papers/10.21105/joss.03118/status.svg)](https://doi.org/10.21105/joss.03118)

If you use `presize`, please cite it in your publication as:  
Haynes et al., (2021). presize: An R-package for precision-based sample size calculation in clinical research. Journal of Open Source Software, 6(60), 3118, https://doi.org/10.21105/joss.03118


### Acknowledgements

The package logo was created with [`ggplot2`](https://ggplot2.tidyverse.org/) and [`hexSticker`](https://github.com/GuangchuangYu/hexSticker) with icons from [Font Awesome](https://fontawesome.com/) (via the [emojifont package](https://github.com/GuangchuangYu/emojifont)).

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differences.R
\name{prec_riskratio}
\alias{prec_riskratio}
\title{Sample size or precision for risk ratio}
\usage{
prec_riskratio(
  p1,
  p2,
  n1 = NULL,
  r = 1,
  conf.width = NULL,
  conf.level = 0.95,
  method = c("koopman", "katz"),
  ...
)
}
\arguments{
\item{p1}{risk among exposed.}

\item{p2}{risk among unexposed.}

\item{n1}{number of patients in exposed group.}

\item{r}{allocation ratio (relative size of unexposed and exposed cohort
(\code{n2} / \code{n1})).}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{method}{Exactly one of \code{koopman} (\emph{default}), \code{katz}.
Methods can be abbreviated.}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}
}
\description{
\code{prec_riskratio} returns the risk ratio and the sample size or the
precision for the provided proportions.
}
\details{
Exactly one of the parameters \code{n1} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

Koopman (\code{koopman}) provides an asymptotic score confidence interval
that is always consistent with Pearsons chi-squared test. It is the
recommended interval (Fagerland et al.).

Katz (\code{katz}) use a logarithmic transformation to calculate the
confidence interval. The CI cannot be computed if one of the proportions is
zero. If both proportions are 1, the estimate of the standard error becomes
zero, resulting in a CI of [1, 1].

\code{\link[stats]{uniroot}} is used to solve n for the katz, and koopman
method.
}
\examples{
# Validate function with example in Fagerland et al. (2015), Table 5.
prec_riskratio(p1 = 7/34, p2 = 1/34, n1 = 34, r = 1, met = "katz")
# 7 (0.91 to 54)
prec_riskratio(p1 = 7/34, p2 = 1/34, n1 = 34, r = 1, met = "koopman")
# 7 (1.21 to 43)

# Validate the Koopman method with example in Koopman (1984)
prec_riskratio(p1 = 36/40, p2 = 16/80, n1 = 40, r = 2, met = "koopman")
# 4.5 (2.94 to 7.15)
}
\references{
Fagerland MW, Lydersen S, and Laake P (2015). \emph{Recommended confidence
intervals for two independent binomial proportions}, Statistical methods in
medical research 24(2):224-254.

Katz D, Baptista J, Azen SP, and Pike MC (1978) \emph{Obtaining Confidence
Intervals for the Risk Ratio in Cohort Studies}, Biometrics 34:469-474.

Koopman PAR (1984) \emph{Confidence Intervals for the Ratio of Two Binomial
Proportions}, Biometrics 40:513-517.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_measures.R
\name{prec_kappa}
\alias{prec_kappa}
\title{Sample size or precision for Cohen's kappa}
\usage{
prec_kappa(
  kappa,
  n = NULL,
  raters = 2,
  n_category = 2,
  props,
  conf.width = NULL,
  conf.level = 0.95
)
}
\arguments{
\item{kappa}{expected value of Cohen's kappa.}

\item{n}{sample size.}

\item{raters}{number of raters (maximum of 6).}

\item{n_category}{number of categories of outcomes (maximum of 5).}

\item{props}{expected proportions of each outcome (should have length
\code{n_category}).}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_kappa} returns the sample size or the precision for the provided Cohen's kappa coefficient.
}
\details{
This function wraps the \code{FixedN} and \code{CI} functions in the
\code{kappaSize} package.
The \code{FixedN} functions in \code{kappaSize} return a one sided confidence
interval. The values that are passed to \code{kappaSize} ensure that two-sided
confidence intervals are returned, although we assume that confidence intervals
are symmetrical.
}
\examples{
# precision based on sample size
#   two categories with proportions of 30 and 70\\%, four raters
prec_kappa(kappa = .5, n = 200, raters = 4, n_category = 2, props = c(.3,.7))
# sample size to get a given precision
prec_kappa(kappa = .5, conf.width = .15, raters = 4, n_category = 2,
           props = c(.3,.7))

# as above, but with two scenarios for kappa
prec_kappa(kappa = c(.5, .75), conf.width = .15, raters = 4, n_category = 2,
           props = c(.3,.7))
prec_kappa(kappa = c(.5, .75), conf.width = c(.15, 0.3), raters = 4,
           n_category = 2, props = c(.3,.7))

}
\seealso{
\code{\link[kappaSize]{FixedNBinary}},
\code{\link[kappaSize]{FixedN3Cats}},
\code{\link[kappaSize]{CIBinary}},
\code{\link[kappaSize]{CI3Cats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/presize-package.R
\docType{package}
\name{presize-package}
\alias{presize}
\alias{presize-package}
\title{presize: Precision Based Sample Size Calculation}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Bland (2009) <doi:10.1136/bmj.b3985> recommended to
    base study sizes on the width of the confidence interval rather the power of 
    a statistical test. The goal of 'presize' is to provide functions for such 
    precision based sample size calculations. For a given sample size, the 
    functions will return the precision (width of the confidence interval), and 
    vice versa.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/CTU-Bern/presize}
  \item \url{https://ctu-bern.github.io/presize/}
  \item Report bugs at \url{https://github.com/CTU-Bern/presize/issues}
}

}
\author{
\strong{Maintainer}: Alan G. Haynes \email{alan.haynes@ctu.unibe.ch}

Authors:
\itemize{
  \item Armando Lenz \email{armando.lenz@ctu.unibe.ch}
  \item Andreas Limacher \email{andreas.limacher@ctu.unibe.ch}
}

Other contributors:
\itemize{
  \item Odile Stalder \email{odile.stalder@ctu.unibe.ch} [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differences.R
\name{prec_rateratio}
\alias{prec_rateratio}
\title{Sample size or precision for a rate ratio}
\usage{
prec_rateratio(
  n1 = NULL,
  rate1 = NULL,
  rate2 = 2 * rate1,
  prec.level = NULL,
  r = 1,
  conf.level = 0.95
)
}
\arguments{
\item{n1}{number of patients in exposed group.}

\item{rate1}{event rate in the exposed group.}

\item{rate2}{event rate in the unexposed group.}

\item{prec.level}{ratio of the upper limit over the lower limit of the
rate ratio confidence interval.}

\item{r}{allocation ratio (relative size of unexposed and exposed cohort
(\code{n2} / \code{n1})).}

\item{conf.level}{confidence level.}
}
\description{
\code{prec_rateratio} returns the sample size or the precision for the
provided proportions.
}
\details{
Exactly one of the parameters  \code{n1} or \code{conf.width} must be passed as
NULL, and that parameter is determined from the other. Event rates in the two
groups should also be provided (\code{rate1, rate2}). If only
\code{rate1} is provided, \code{rate2} is assumed to be 2 times
\code{rate1}.
}
\examples{
# 20 participants, a rate of 50\%  against a rate of 300\\%
prec_rateratio(20, .5, 3)
# sample size required to attain a CI whose upper limit is not more than 3.81 larger
#  than the lower limit
prec_rateratio(rate1 = .5, rate2 = 3, prec.level = 3.81)
}
\references{
Rothman KJ, Greenland S (2018). \emph{Planning Study Size Based on
  Precision Rather Than Power}. Epidemiology, 29:599-603.
  \doi{10.1097/EDE.0000000000000876}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_measures.R
\name{prec_icc}
\alias{prec_icc}
\title{Sample size or precision for an intraclass correlation}
\usage{
prec_icc(rho, k, n = NULL, conf.width = NULL, conf.level = 0.95)
}
\arguments{
\item{rho}{desired intraclass correlation.}

\item{k}{number of observations per n (subject).}

\item{n}{number of subjects.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_icc} returns the sample size or the precision for the given
intraclass correlation.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the others.

Sample size or precision is calculated according to formula 3 in Bonett
(2002), which is an approximation. Whether ICC is calculated for a one-way or
a two-way ANOVA does not matter in the approximation. As suggested by the
author, \eqn{5*rho} is added to n, if \eqn{k = 2} and \eqn{rho \ge 7}.

n is rounded up to the next whole number using \code{ceiling}.
}
\examples{
# Bonett (2002) gives an example using 4 raters, with an ICC of 0.85 and want
# a confidence width of 0.2. Bonett calculated that a sample size of 19.2 was
# required. This can be done via
prec_icc(0.85, 4, conf.width = 0.2)
# note that \code{presamp} rounds up to the nearist integer.

# Bonett then goes on to estimate the width given the sample size, finding a
# value 'close to 0.2':
prec_icc(0.85, 4, 20)
}
\references{
Bonett DG (2002). \emph{Sample size requirements for estimating
  intraclass correlations with desired precision}. Statistics in Medicine,
  21:1331-1335. \doi{10.1002/sim.1108}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differences.R
\name{prec_riskdiff}
\alias{prec_riskdiff}
\title{Sample size or precision for risk difference}
\usage{
prec_riskdiff(
  p1,
  p2,
  n1 = NULL,
  conf.width = NULL,
  r = 1,
  conf.level = 0.95,
  method = c("newcombe", "mn", "ac", "wald"),
  ...
)
}
\arguments{
\item{p1}{risk among exposed.}

\item{p2}{risk among unexposed.}

\item{n1}{number of patients in exposed group.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{r}{allocation ratio (relative size of exposed and unexposed cohort
(\code{n1} / \code{n2})).}

\item{conf.level}{confidence level.}

\item{method}{Exactly one of \code{newcombe} (\emph{default}), \code{mn}
(Miettinen-Nurminen), \code{ac} (Agresti-Caffo), \code{wald}. Methods can
be abbreviated.}

\item{...}{other options to uniroot (e.g. \code{tol})}
}
\description{
\code{prec_riskdiff} returns the risk difference and the sample size or the
precision for the provided proportions.
}
\details{
Exactly one of the parameters \code{n1} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

Newcombe (\code{newcombe}) proposed a confidence interval based on the wilson
score method for the single proportion (see \link{prec_prop}). The confidence
interval without continuity correction is implemented from equation 10 in
Newcombe (1998).

Miettinen-Nurminen (\code{mn}) provide a closed from equation for the
restricted maximum likelihood estimate . The implementation is based on
code provided by Yongyi Min on
\url{http://users.stat.ufl.edu/~aa/cda/R/two-sample/R2/index.html}.

Agresti-Caffo (\code{ac}) confidence interval is based on the Wald confidence
interval, adding 1 success to each cell of the 2 x 2 table (see Agresti and
Caffo 2000).

\code{\link[stats]{uniroot}} is used to solve n for the newcombe, ac, and mn
method.
}
\examples{
# proportions of 40 and 30\\%, 50 participants, how wide is the CI?
prec_riskdiff(p1 = .4, p2 = .3, n1 = 50)
# proportions of 40 and 30\\%, 50 participants, how many participants for a CI 0.2 wide?
prec_riskdiff(p1 = .4, p2 = .3, conf.width = .2)

# Validate Newcombe (1998)
prec_riskdiff(p1 = 56/70, p2 = 48/80, n1 = 70, r = 70/80, met = "newcombe")  # Table IIa
prec_riskdiff(p1 = 10/10, p2 = 0/10, n1 = 10, met = "newcombe")  # Table IIh

# multiple scenarios
prec_riskdiff(p1 = c(56/70, 9/10, 6/7, 5/56),
              p2 = c(48/80, 3/10, 2/7, 0/29),
              n1 = c(70, 10, 7, 56),
              r = c(70/80, 1, 1, 56/29),
              method = "wald")

}
\references{
Agresti A (2003) \emph{Categorical Data Analysis}, Second Edition, Wiley
Series in Probability and Statistics,
\doi{10.1002/0471249688}.

Agresti A and Caffo B (2000) \emph{Simple and Effective Confidence Intervals
for Proportions and Differences of Proportions Result from Adding Two
Successes and Two Failures}, The American Statistician, 54(4):280-288.

Miettinen O and Nurminen M (1985) \emph{Comparative analysis of two rates},
Statistics in Medicine, 4:213-226.

Newcombe RG (1998) \emph{Interval estimation for the difference between
independent proportions: comparison of eleven methods}, Statistics in
Medicine, 17:873-890.

Fagerland MW, Lydersen S, and Laake P (2015). \emph{Recommended confidence
intervals for two independent binomial proportions}, Statistical methods in
medical research 24(2):224-254.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.R
\name{prec_lr}
\alias{prec_lr}
\alias{prec_pos_lr}
\alias{prec_neg_lr}
\title{Sample size or precision for likelihood ratios}
\usage{
prec_lr(prev, p1, p2, n = NULL, conf.width = NULL, conf.level = 0.95, ...)

prec_pos_lr(
  prev,
  sens,
  spec,
  n = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  ...
)

prec_neg_lr(
  prev,
  sens,
  spec,
  n = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  ...
)
}
\arguments{
\item{prev}{disease/case prevalence in the study group.}

\item{p1}{proportion of positives in group 1 (e.g. sensitivity).}

\item{p2}{proportion of positives in group 2 (e.g. 1 - specificity).}

\item{n}{total group size.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level (defaults to 0.95).}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}

\item{sens}{sensitivity.}

\item{spec}{specificity.}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
These functions calculate the precision or sample size for likelihood ratios (LRs).
\code{prec_lr} is a generalized method for that can be used for positive and
negative LRs as well as conditional LRs.

\code{prec_pos_lr} is a wrapper to \code{prec_lr} to ease
calculations for positive likelihood ratios by allowing sensitivity and
specificity to be given explicitly.

\code{prec_neg_lr} is a wrapper to \code{prec_lr} to ease
calculations for negative likelihood ratios by allowing sensitivity and
specificity to be given explicitly.
}
\details{
These functions implement formula 10 from Simel et al 1991.
\code{prec_lr} is a generalized function allowing for many scenarios, while
\code{prec_pos_lr} and \code{prec_neg_lr} are specific to positive and
negative likelihood ratios in the 2*2 setting (e.g. disease status and test
positive/negative).

For the positive likelihood ratio (LR+), in a 2x2 style experiment, \code{p1}
should be sensitivity, \code{p2} should be 1-specificity. Alternatively, use
\code{prec_pos_lr}.

For the negative likelihood ratio (LR-), in a 2x2 style experiment, \code{p1}
should be 1-sensitivity, \code{p2} should be specificity. Alternatively, use
\code{prec_neg_lr}.

For conditional likelihood ratios with 3x2 tables, such as positive or
negative tests against inconclusive ones (yields), \code{p1} would be the
proportion of positive or negative tests in the diseased group and \code{p2}
would be the proportion of positive or negative tests in the non-diseased group.
}
\section{Functions}{
\itemize{
\item \code{prec_pos_lr}: "Positive likelihood ratio"

\item \code{prec_neg_lr}: "Negative likelihood ratio"
}}

\examples{
# equal numbers of diseased/non-diseased, 80\% sens, 73\% spec, 74 participants total
prec_lr(.5, .8, .27, 74)

# Simel et al 1991, problem 1 - LR+ CI width from N
# Sensitivity of a new test is at least 80\%, specificity is 73\% and the LR+
# is 2.96 (= 0.8/(1-0.73)). We have as many diseased as not diseased
# (n1 = n2, n = 2*n1 = 146.8, prevalence = .5)
prec_lr(prev = .5, p1 = .8, p2 = 1-.73, n = 146.8)
prec_pos_lr(prev = .5, sens = .8, spec = .73, n = 146.8)

# problem 1 of Simel et al actually derives n1 rather than the width of the
# confidence interval (ie N from CI width). If we know that the lower limit
# of the CI should be 2.0, the confidence interval width is approximately
# exp(2*(log(2.96) - log(2))) = 2.19 (approximate because the CI Of the LR
# is only symetrical on the log(LR) scale), which we can put in conf.width
prec_lr(prev = .5, p1 = .8, p2 = 1-.73, conf.width = 2.2)
# same, but using the wrapper to specify sens and spec
prec_pos_lr(prev = .5, sens = .8, spec = .73, conf.width = 2.2)

# Simel et al 1991, problem 2 - LR- CI width from N
# p1 = 1 - sens = .1, p2 = spec = .5
# n1 = n2, n = 160, prev = .5
prec_lr(prev = .5, p1 = .1, p2 = .5, n = 160)
# same, but using the wrapper to specify sens and spec
prec_neg_lr(prev = .5, sens = .9, spec = .5, n = 160)

}
\references{
Simel, DL, Samsa, GP and Matchar, DB (1991) \emph{Likelihood ratios with confidence: Sample size estimation for diagnostic test studies.} J Clin Epidemiol 44(8), 763-770
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differences.R
\name{prec_or}
\alias{prec_or}
\title{Sample size or precision for an odds ratio}
\usage{
prec_or(
  p1,
  p2,
  n1 = NULL,
  r = 1,
  conf.width = NULL,
  conf.level = 0.95,
  method = c("gart", "woolf", "indip_smooth"),
  ...
)
}
\arguments{
\item{p1}{risk among exposed.}

\item{p2}{risk among unexposed.}

\item{n1}{number of patients in exposed group.}

\item{r}{allocation ratio (relative size of unexposed and exposed cohort
(\code{n2} / \code{n1})).}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{method}{Exactly one of \code{indip_smooth} (\emph{default}),
\code{gart}, or \code{woolf}. Methods can be abbreviated.}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_or} returns the sample size or the precision for the
provided proportions.
}
\details{
Exactly one of the parameters \code{n1} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

Woolf (\code{woolf}), Gart (\code{gart}), and Independence-smoothed logit
(\code{indip_smooth}) belong to a general family of adjusted confidence
intervals, adding 0 (woolf) to each cell, 0.5 (gart) to each cell, or an
adjustment for each cell based on observed data (independence-smoothed). In
gart and indip_smooth, estimate of the CI is not possible if \eqn{p1 = 0}, in
which case the OR becomes 0, but the lower level of the CI is > 0. Further,
if \eqn{p1 = 1} and \eqn{p2 < 1}, or if \eqn{p1 > 0} and \eqn{p2 = 0}, the OR
becomes \eqn{\infty}, but the upper limit of the CI is finite. For the
approximate intervals, \code{gart} and \code{indip_smooth} are the
recommended intervals (Fagerland et al. 2011).

\code{\link[stats]{uniroot}} is used to solve n for the woolf, gart, and
indip_smooth method.
}
\examples{
# 10\\% events in one group, 15\\% in the other, 200 participants total
#  (= 100 in each group), estimate confidence interval width
prec_or(p1 = .1, p2 = .15, n1 = 200/2)
# formula by Gart
prec_or(p1 = .1, p2 = .15, n1 = 200/2, method = "gart")
# formula by Woolf
prec_or(p1 = .1, p2 = .15, n1 = 200/2, method = "woolf")

# 10\\% odds in one group, 15\\% in the other, desired CI width of 0.1,
#  estimate N
prec_or(p1 = .1, p2 = .15, conf.width = .1)
# formula by Gart
prec_or(p1 = .1, p2 = .15, conf.width = .1, method = "gart")
# formula by Woolf
prec_or(p1 = .1, p2 = .15, conf.width = .1, method = "woolf")

}
\references{
Fagerland MW, Lydersen S, Laake P (2015). \emph{Recommended
confidence intervals for two independent binomial proportions}. Statistical
Methods in Medical Research, 24(2):224-254.
\doi{10.1177/0962280211415469}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differences.R
\name{prec_meandiff}
\alias{prec_meandiff}
\title{Sample size or precision for a mean difference}
\usage{
prec_meandiff(
  delta,
  sd1,
  sd2 = sd1,
  n1 = NULL,
  r = 1,
  conf.width = NULL,
  conf.level = 0.95,
  variance = c("equal", "unequal"),
  ...
)
}
\arguments{
\item{delta}{difference in means between the two groups.}

\item{sd1}{standard deviation in group 1.}

\item{sd2}{standard deviation in group 2.}

\item{n1}{number of patients in group 1.}

\item{r}{allocation ratio (relative size of group 2 and group 1 (n2 / n1)).}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{variance}{\code{equal} (\emph{default}) or \code{unequal} variance.}

\item{...}{other options to uniroot (e.g. \code{tol})}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_meandiff} returns the sample size or the precision for the
provided mean difference and standard deviations.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.
}
\examples{
# mean difference of 5, SD of 2.5, CI width with 20 participants assuming equal variances
prec_meandiff(delta = 5, sd1 = 2.5, n1 = 20, var = "equal")
# mean difference of 5, SD of 2.5, number of participants for a CI width of 3,
#  assuming equal variances
prec_meandiff(delta = 5, sd1 = 2.5, conf.width = 3, var = "equal")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptive_stats.R
\name{prec_prop}
\alias{prec_prop}
\title{Sample size or precision for a proportion}
\usage{
prec_prop(
  p,
  n = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  method = c("wilson", "agresti-coull", "exact", "wald"),
  ...
)
}
\arguments{
\item{p}{proportion.}

\item{n}{number of observations.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{method}{The method to use to calculate precision. Exactly one method
may be provided. Methods can be abbreviated.}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements. In the wilson and
  agresti-coull formula, the p from which the confidence interval is
  calculated is adjusted by a term (i.e. \eqn{p + term \pm ci}). This
  adjusted p is returned in \code{padj}.
}
\description{
\code{prec_prop} returns the sample size or the precision for the provided
proportion.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

The wilson, agresti-coull, exact, and wald method are implemented. The
wilson method is suggested for small \code{n} (< 40), and the agresti-coull method
is suggested for larger \code{n} (see reference). The wald method is not suggested,
but provided due to its widely distributed use.

\code{\link[stats]{uniroot}} is used to solve \code{n} for the agresti-coull,
wilson, and exact methods. Agresti-coull can be abbreviated by ac.
}
\examples{
# CI width for 15\\% with 50 participants
prec_prop(0.15, n = 50)
# number of participants for 15\\% with a CI width of 0.2
prec_prop(0.15, conf.width = 0.2)
# confidence interval width for a range of scenarios between 10 and 90\\% with
#  100 participants via the wilson method
prec_prop(p = 1:9 / 10, n = 100, method = "wilson")
# number of participants for a range of scenarios between 10 and 90\\% with
#  a CI of 0.192 via the wilson method
prec_prop(p = 1:9 / 10, conf.width = .192, method = "wilson")
}
\references{
Brown LD, Cai TT, DasGupta A (2001) \emph{Interval Estimation for
  a Binomial Proportion}, Statistical Science, 16:2, 101-117,
  \doi{10.1214/ss/1009213286}
}
\seealso{
\code{\link[stats]{binom.test}}, \code{\link[binom]{binom.confint}}
  in package \pkg{binom}, and \code{\link[Hmisc]{binconf}} in package
  \pkg{Hmisc}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_measures.R
\name{prec_lim_agree}
\alias{prec_lim_agree}
\title{Sample size or precision for limit of agreement on Bland-Altman plots}
\usage{
prec_lim_agree(n = NULL, conf.width = NULL, conf.level = 0.95)
}
\arguments{
\item{n}{sample size.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_lim_agree} returns the sample size or the precision for the limit
of agreement, i.e. the confidence interval around the limit of agreement,
expressed in SD-units. It is an approximation based on the Normal distribution,
instead of a Student t distribution.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

The sample size and precision are calculated according to formulae in Bland &
Altman (1986). The CI width is a simple function of the sample size only.
}
\examples{
# calculate confidence interval width, given N
prec_lim_agree(200)
# calculate N given, confidence interval width
prec_lim_agree(conf.width = .1)
}
\references{
Bland & Altman (1986) \emph{Statistical methods for assessing agreement
between two methods of clinical measurement} Lancet i(8476):307-310
\doi{10.1016/S0140-6736(86)90837-8}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyApp.R
\name{launch_presize_app}
\alias{launch_presize_app}
\title{Presize shiny app}
\usage{
launch_presize_app()
}
\description{
Besides the programmatic approach to using presize, we also supply a shiny app,
enabling point-and-click interaction with the program. The app will open in a
new window. Select the appropriate method from the menu on the left and enter
the relevant parameters indicated in the panel on the right. The output is then
displayed lower down the page.
}
\details{
The main disadvantage to the app is that it only allows a single scenario at
a time.

The app is also available at \href{https://shiny.ctu.unibe.ch/presize/}{https://shiny.ctu.unibe.ch/presize/}.

\if{html}{\figure{app.png}{options: width="100\%"}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_measures.R
\name{prec_cor}
\alias{prec_cor}
\title{Sample size or precision for correlation coefficient}
\usage{
prec_cor(
  r,
  n = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  method = c("pearson", "kendall", "spearman"),
  ...
)
}
\arguments{
\item{r}{desired correlation coefficient.}

\item{n}{sample size.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{method}{Exactly one of \code{pearson} (\emph{default}), \code{kendall},
or \code{spearman}. Methods can be abbreviated.}

\item{...}{other options to uniroot (e.g. \code{tol})}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_cor} returns the sample size or the precision for the given
pearson, spearman, or kendall correlation coefficient.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

Sample size or precision is calculated according to formula 2 in Bonett and
Wright (2000). The use of pearson is only recommended, if \eqn{n \ge 25}. The
pearson correlation coefficient assumes bivariate normality. If the
assumption of bivariate normality cannot be met, spearman or kendall should
be considered.

n is rounded up to the next whole number using \code{ceiling}.

\code{\link[stats]{uniroot}} is used to solve n.
}
\examples{
# calculate confidence interval width...
# Pearson correlation coefficient
prec_cor(r = 0.5, n = 100)
# Kendall rank correlation coefficient (tau)
prec_cor(r = 0.5, n = 100, method = "kendall")
# Spearman's rank correlation coefficient
prec_cor(r = 0.5, n = 100, method = "spearman")
# calculate N required for a given confidence interval width...
# Pearson correlation coefficient
prec_cor(r = 0.5, conf.width = .15)
# Kendall rank correlation coefficient (tau)
prec_cor(r = 0.5, conf.width = .15, method = "kendall")
# Spearman's rank correlation coefficient
prec_cor(r = 0.5, conf.width = .15, method = "spearman")
}
\references{
Bonett DG, and Wright TA (2000) \emph{Sample size requirements
  for estimating Pearson, Kendall and Spearman correlations} Psychometrika
  65:23-28. \doi{10.1007/BF02294183}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptive_stats.R
\name{prec_mean}
\alias{prec_mean}
\title{Sample size or precision for a mean}
\usage{
prec_mean(
  mean,
  sd,
  n = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  ...,
  mu = NULL
)
}
\arguments{
\item{mean}{mean.}

\item{sd}{standard deviation.}

\item{n}{number of observations.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}

\item{mu}{deprecated argument}
}
\value{
Object of class "presize", a list with \code{mean} mean, \code{sd} standard deviation, \code{n} sample size,
\code{conf.width} precision (the width of the confidence interval),
\code{lwr} lower bound of confidence interval, \code{upr} upper bound of confidence interval,
 augmented with method and note elements.
}
\description{
\code{prec_mean} returns the sample size or the precision for the provided
mean and standard deviation.
}
\details{
Exactly one of the parameters \code{n} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

The precision is defined as the full width of the confidence interval. The
confidence interval calculated as \eqn{t(n - 1) * sd / sqrt(n)}, with t(n-1)
from the t-distribution with n-1 degrees of freedom.

\code{\link[stats]{uniroot}} is used to solve \code{n}.
}
\examples{
# mean of 5, SD of 2.5, whats the confidence interval width with 20 participants?
prec_mean(mean = 5, sd = 2.5, n = 20)
# mean of 5, SD of 2.5, how many participants for CI width of 2.34?
prec_mean(mean = 5, sd = 2.5, conf.width = 2.34)  # approximately the inverse of above
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.R
\name{prec_sens}
\alias{prec_sens}
\alias{prec_spec}
\title{Sample size and precision of sensitivity and specificity}
\usage{
prec_sens(
  sens,
  n = NULL,
  ntot = NULL,
  prev = NULL,
  conf.width = NULL,
  round = "ceiling",
  ...
)

prec_spec(
  spec,
  n = NULL,
  ntot = NULL,
  prev = NULL,
  conf.width = NULL,
  round = "ceiling",
  ...
)
}
\arguments{
\item{sens, spec}{proportions.}

\item{n}{number of observations.}

\item{ntot}{total sample size.}

\item{prev}{prevalence of cases/disease (i.e. proportion of \code{ntot} with
the disease).}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{round}{string, round calculated \code{n} up (\code{ceiling}) or down
(\code{floor}).}

\item{...}{options passed to prec_prop (e.g. method,
conf.width, conf.level).}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
Because sensitivity and specificity are simple proportions, these functions
act as wrappers for \code{prec_prop}.
}
\details{
If \code{ntot} and \code{prev} are given, they are used to calculate
  \code{n}.
}
\note{
Calculated \code{n} can take on non-integer numbers, but
  \code{prec_prop} requires integers, so the calculated \code{n} is rounded
  according to the approach indicated in \code{round}.
}
\examples{
  # confidence interval width with n
  prec_sens(.6, 50)
  # confidence interval width with ntot and prevalence (assuming 50\% prev)
  prec_sens(.6, ntot = 100, prev = .5)
  # sample size with confidence interval width
  prec_sens(.6, conf.width = 0.262)
}
\seealso{
\code{prec_prop}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptive_stats.R
\name{prec_rate}
\alias{prec_rate}
\title{Sample size or precision for a rate}
\usage{
prec_rate(
  r,
  x = NULL,
  conf.width = NULL,
  conf.level = 0.95,
  method = c("score", "vs", "exact", "wald"),
  ...
)
}
\arguments{
\item{r}{rate or rate ratio.}

\item{x}{number of events.}

\item{conf.width}{precision (the full width of the confidence interval).
Should not exceed 5 times \code{r}.}

\item{conf.level}{confidence level.}

\item{method}{The method to use to calculate precision. Exactly one method
may be provided. Methods can be abbreviated.}

\item{...}{other arguments to uniroot (e.g. \code{tol}).}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
\code{prec_rate} returns the sample size or the precision for the provided
rate.
}
\details{
Exactly one of the parameters \code{r} or \code{conf.width} must be passed as NULL,
and that parameter is determined from the other.

The \code{score}, variance stabilizing (\code{vs}), \code{exact}, and
\code{wald} method are implemented to calculate the rate and the precision.
For few events \code{x} (<5), the exact method is recommended.

If more than one method is specified or the method is miss-specified, the
'score' method will be used.

\code{\link[stats]{uniroot}} is used to solve n for the score and
exact method.
}
\examples{
# confidence interval width for a rate of 2.5 events per unit and 20 events,
#  using the score method
prec_rate(2.5, x = 20, met = "score")
# number of events to yield a CI width of 2.243 for a rate of 2.5 events per
#  unit and 20 events, using the score method
prec_rate(2.5, conf.width = 2.243, met = "score")
# confidence interval width for a rate of 2.5 events per unit and 20 events,
#  using the exact method
prec_rate(2.5, x = 20, met = "exact")
# vs and wald have the same conf.width, but different lwr and upr
prec_rate(2.5, x = 20, met = "vs")
prec_rate(2.5, x = 20, met = "wald")
}
\references{
Barker, L. (2002) \emph{A Comparison of Nine Confidence Intervals
for a Poisson Parameter When the Expected Number of Events is \eqn{\le} 5},
The American Statistician, 56:2, 85-89,
\doi{10.1198/000313002317572736}
}
\seealso{
\code{\link[stats]{poisson.test}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.R
\name{prec_auc}
\alias{prec_auc}
\title{Sample size or precision for AUC}
\usage{
prec_auc(auc, prev, n = NULL, conf.width = NULL, conf.level = 0.95, ...)
}
\arguments{
\item{auc}{AUC value.}

\item{prev}{prevalence.}

\item{n}{number of observations.}

\item{conf.width}{precision (the full width of the confidence interval).}

\item{conf.level}{confidence level.}

\item{...}{other arguments to \code{optimize}.}
}
\value{
Object of class "presize", a list of arguments (including the
  computed one) augmented with method and note elements.
}
\description{
Calculate the sample size from AUC, prevalence and confidence interval width
or the expected confidence interval width from AUC, prevalence and sample
size, following Hanley and McNeil (1982).
}
\details{
Sample size is derived by optimizing the difference between the difference
between the lower and upper limits of the confidence interval and
\code{conf.width}.
}
\examples{
# confidence interval width
N <- 500
prev <- .1
auc <- .65
(prec <- prec_auc(auc, prev, n = N))
cwidth <- prec$conf.width
# sample size
prec_auc(auc, prev, conf.width = cwidth)
}
\references{
Hanley, JA and McNeil, BJ (1982) \emph{The Meaning and Use of the Area under a Receiver Operating Characteristic (ROC) Curve.} Radiology 148, 29-36
}
