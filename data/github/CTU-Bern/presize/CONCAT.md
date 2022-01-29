
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
