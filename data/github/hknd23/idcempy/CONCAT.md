---
title: 'IDCeMPy: Python Package for Inflated Discrete Choice Models'

authors:
- affiliation: 1
  name: Nguyen K. Huynh
  orcid: 0000-0002-6234-7232
- affiliation: 2
  name: Sergio Béjar
  orcid: 0000-0002-9352-3892
- affiliation: 1
  name: Vineeta Yadav
- affiliation: 1
  name: Bumba Mukherjee
date: "21 June 2021"
output:
  html_document:
    df_print: paged
  pdf_document: default
bibliography: paper.bib

tags:
- Python, Inflated Ordered Probit Models, Generalized Inflated MNL Models

affiliations:
- index: 1
  name: Dept. of Political Science, Pennsylvania State University
- index: 2
  name: Dept. of Political Science, San Jose State University
  
---
# Summary

Scholars and data scientists often use discrete choice models to evaluate ordered dependent variables using the ordered probit model and unordered polytomous outcome measures via the multinomial logit (MNL) estimator [@greene2002nlogit; @JSSv074i10; @richards2018new]. These models, however, cannot account for the possibility that in many ordered and unordered polytomous choice outcomes, a disproportionate share of observations — stemming from two distinct data generating processes (DGPs) — fall into a single category which is thus “inflated.” For instance, ordered outcome measures of self-reported smoking behavior that range from 0 for “no smoking” to 3 for “smoking 20 cigarettes or more daily” contain excessive observations in the zero (no smoking) category that includes individuals who never smoke cigarettes and those who smoked previously but temporarily stop smoking because of an increase in cigarette costs [@harris2007zero; @greene2015inflated]. The “indifference” middle-category in ordered measures of immigration attitudes is inflated since it includes respondents who are genuinely indifferent about immigration and those who select “indifference” because of social desirability reasons [@bagozzi2012mixture; @brown2020modelling]. The baseline category of unordered polytomous variables of presidential vote choice is also often inflated as it includes non-voters who abstain from voting owing to temporary factors and routine non-voters who are disengaged from the political process [@campbell2008religion; @bagozzi2017distinguishing].  Inflated discrete choice models have been developed to address such category inflation in ordered and unordered polytomous outcome variables as failing to do so leads to model misspecification and incorrect inferences [@harris2007zero; @bagozzi2012mixture; @brown2020modelling].

`IDCeMPy` is an open-source Python package that enables researchers to fit three distinct sets of discrete choice models used by data scientists, economists, engineers, political scientists, and public health researchers: the Zero-Inflated Ordered Probit (ZiOP) model without and with correlated errors (ZiOPC model), Middle-Inflated Ordered Probit (MiOP) model without and with correlated errors (MiOPC), and Generalized-Inflated Multinomial Logit (GiMNL) models. Functions that fit the ZiOP(C) model in `IDCeMPy` evaluate zero-inflated ordered dependent variables that result from two DGPs, while functions that fit the MiOP(C) models account for inflated middle-category ordered outcomes that emerge from distinct DGPs. The functions in `IDCeMPy` that fit GiMNL models account for the large share and heterogeneous mixture of observations in the baseline and other lower outcome categories in unordered polytomous dependent variables. The primary location for the description of the functions that fit the models listed above is available at the [IDCeMPy package’s documentation website](https://idcempy.readthedocs.io/en/latest/).

# State of the Field

Software packages and code are available for estimating standard (non-inflated) discrete choice models. In the R environment, the `MASS` [@venables2002random] and `micEcon` [@henningsen2014micecon] packages fit binary and discrete choice models. The `Rchoice` [@JSSv074i10] package allows researchers to estimate binary and ordered probit and logit models as well as the Poisson model by employing various optimization routines. The proprietary LIMDEP package NLOGIT [@greene2002nlogit] fits conventional binary and ordered discrete choice models but is neither open-sourced nor freely available. The R `mlogit` [@croissant2012estimation] and `mnlogit` [@hasan2014fast] packages provide tools for working with conventional MNL models, while `gmnl` [@sarrias2017multinomial] and `PReMiuM` [@liverani2015premium] estimate MNL models that incorporate unit-specific heterogeneity. There are proprietary LIMDEP software and R code — but not an R package — that fit few inflated ordered probit and MNL models [@harris2007zero; @bagozzi2012mixture; @bagozzi2017distinguishing]. Outside R, the Python `biogeme` [@bierlaire2016pythonbiogeme] package fits mixed logit and MNL models. Further, @dale2021estimation’s ZiOP STATA command (but not package) fits the Zero-Inflated Ordered Probit without correlated errors. @xia2019gidm’s `gidm` STATA command fits discrete choice models without correlated errors for inflated zero and other lower-category discrete outcomes. 

The R or LIMDEP software, along with the STATA commands listed above, are undoubtedly helpful. However, to our knowledge, there are no R or Python packages to fit a variety of statistical models that account for the excessive (i.e., “inflated”) share of observations in the baseline, and other higher categories of ordered and unordered polytomous dependent variables, which are commonly analyzed across the natural and social sciences. As discussed below, our Python package `IDCeMPy` thus fills an important lacuna by providing an array of functions that fit a substantial range of inflated discrete choice models applicable across various disciplines.

# Statement of Need 

Although our `IDCeMPy` package also fits standard discrete choice models, what makes it unique is that unlike existing software, it offers functions to fit and assess the performance of both Zero-Inflated and Middle-Inflated Ordered Probit (OP) models without and with correlated errors as well as a set of Generalized-Inflated MNL models. The models included in `IDCeMPy` account for the excessive proportion of observations in any given ordered or unordered outcome category by combining a single binary probit or logit split-stage equation with either an ordered probit outcome stage (for the Zero and Middle-Inflated OP models) or an MNL outcome-stage equation. Users can treat the error terms from the two equations in the Zero and Middle-Inflated OP models as independent or correlated in the package’s estimation routines. `IDCeMPy` also provides functions to assess each included model’s goodness-of-fit via the AIC statistics, extract the covariates’ marginal effects from each model, and conduct Vuong tests for comparing the performance between the standard and inflated discrete choice models. 

The functions in `IDCeMPy` use quasi-Newton optimization methods such as the Broyden-Fletcher-Goldfarb-Shanno algorithm for Maximum-Likelihood-Estimation (MLE), which facilitates convergence and estimation speed. Another feature is that the coefficients, standard errors, and confidence intervals obtained for each model estimated in `IDCeMPy` are in `pandas.DataFrame` [@mckinney-proc-scipy-2010] format and are stored as class attribute `.coefs`. This allows for easy export to CSV or Excel, which makes it easier for users to perform diagnostic tests and extract marginal effects. `IDCeMPy` is thus essential as it provides a much-needed unified software package to fit statistical models to account for category inflation in several ordered and unordered outcome variables used across fields as diverse as economics, engineering, marketing, political science, public health, sociology, and transportation research. Users can employ the wide range of statistical models in `IDCeMPy` to assess:

- Zero-inflation in self-reported smoking behavior [@harris2007zero], demand for health treatment [@greene2015inflated], and accident injury-severity [@fountas2018analysis].

- Middle-category inflation in ordered measures of monetary policy [@brown2020modelling] and European Union (EU) membership attitudes [@elgun2007exposure].

- Inflated unordered polytomous outcomes such as transportation choice, environmental policy and consumer demand [@richards2018new], and Presidential vote choice [@campbell2008religion].

# Functionality and Applications

`IDCeMPy` contains the functions listed below to estimate via MLE the following inflated discrete choice models listed earlier:

* `opmod`; `iopmod`; `iopcmod`: Fits the ordered probit model, the Zero-Inflated (ZIOP) and Middle-Inflated ordered probit (MIOP) models without correlated errors, and the ZIOPC and MIOPC models that incorporate correlated errors.

* `opresults`; `iopresults`; `iopcresults`: Presents covariate estimates, Variance-Covariance (VCV) matrix, Log-Likelihood, and AIC statistics of the object models.

* `iopfit`; `iopcfit`: Computes fitted probabilities from each estimated model’s objects.

* `vuong_opiop`; `vuong_opiopc`: Calculates Vuong test statistic for comparing the performance of the OP with the ZiOP(C) and MiOP(C) models.

* `split_effects`; `ordered_effects`: Estimates marginal effects of covariates in the split-stage and outcome-stage respectively. 

* `mnlmod`; `gimnlmod`: Fits MNL model and Generalized-Inflated MNL models.

* `mnlresults`; `gimnlresults`; `vuong_gimnl`: Presents covariate estimates, VCV matrix, Log-Likelihood, and AIC statistics of `mnlmod`; `gimnlmod`. Vuong test statistic for comparing MNL to GIMNL models obtained from `vuong_gimnl`. 


Details about the functionality summarized above are available at the [package’s documentation website](https://idcempy.readthedocs.io/en/latest/), which is open-source and hosted by [ReadTheDocs](https://readthedocs.org/). The features of the functions in `IDCeMPy` that fit the 

(i)	ZiOP(C) models are presented using the ordered self-reported tobacco consumption dependent variable from the [2018 National
Youth Tobacco Dataset](https://www.cdc.gov/tobacco/data_statistics/surveys/nyts/index.htm) 

(ii)	MiOP(C) models are illustrated using the ordered EU support outcome variable from @elgun2007exposure 

(iii)	GiMNL models are evaluated using the unordered polytomous Presidential vote choice dependent variable from @campbell2008religion 


# Availability and Installation

`IDCeMPy` is open-source software made available under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0). It can be installed from [PyPI](https://pypi.org/project/idcempy/) or from its [GitHub repository](https://github.com/hknd23/idcempy). 

# References
# IDCeMPy: Python Package for Inflated Discrete Choice Models

*Nguyen K. Huynh, Sergio Bejar, Vineeta Yadav, Bumba Mukherjee*

<!-- badges: start -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03322/status.svg)](https://doi.org/10.21105/joss.03322)
[![PyPI](https://img.shields.io/pypi/v/idcempy?color=brightgreen)](https://pypi.org/project/idcempy/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/idcempy)](https://pypi.org/project/idcempy/)
[![Downloads](https://pepy.tech/badge/idcempy)](https://pepy.tech/project/idcempy)
[![Project Status: Active â€“ The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**IDCeMPy** is a Python package that provides functions to fit and assess the performance of the following distinct
sets of “inflated” discrete choice models.

* Fit the Zero-Inflated Ordered Probit (ZIOP) model without and with correlated errors (ZIOPC
model) to evaluate zero-inflated ordered choice outcomes that result from a dual data generating
process (d.g.p.).
* Fit the Middle-Inflated Ordered Probit (MIOP) model without and with correlated errors (MIOPC) to account for the inflated middle-category in ordered choice measures related to a dual d.g.p.
* Fit Generalized Inflated Multinomial Logit (GIMNL) models account for the predominant and heterogeneous share of observations in the baseline or any lower category in unordered polytomous choice outcomes.
* Compute AIC and Log-likelihood statistics and the Vuong Test statistic to assess the performance of each inflated discrete choice model in the package.

**IDCeMPy** uses Newton numerical optimization methods to estimate the inflated discrete choice models listed above via Maximum Likelihood Estimation (MLE).  
**IDCeMPY** is compatible with [Python](https://python.org) 3.7+

## Why **IDCeMPy**?

An excessive (“inflated”) share of observations—stemming from two distinct d.g.p’s—fall into a single choice category in many ordered and unordered polytomous outcome variables. Standard Ordered Probit and Multinomial Logit models cannot account for such category inflation which leads to biased inferences. Examples include,

*	The inflated zero-category of "no smoking" in ordered measures of self-reported smoking behavior is generated from nonsmokers who never smoke cigarettes and those who smoked previously but temporarily stopped smoking because of high cigarette prices.

*	The inflated "indifference" middle-category in ordered measures of immigration attitudes includes respondents truly indifferent to immigration and those that choose indifference for social desirability reasons.  

*	The inflated baseline or other lower outcome categories of unordered polytomous outcome measures of vote choice include nonvoters who temporarily abstain from voting and routine nonvoters who always abstain.

**IDCeMPy** includes the ZIOP(C) models for evaluating zero-inflated ordered choice outcomes that result from a dual d.g.p, the MIOP(C) models that address inflated middle-category ordered outcome measures arising from distinct d.g.p’s, and GIMNL models that account for inflated baseline or other categories for unordered polytomous outcomes.

Each inflated discrete choice model in this package addresses category inflation in one’s discrete outcome—unordered or unordered polytomous—of interest by jointly estimating a binary split-stage equation and an ordered or multinomial discrete choice outcome equation.   

## Functions in the **IDCeMPy** Package

| Function         | Description                                                                                                          |
| ---------------- | -------------------------------------------------------------------------------------------------------------------- |
| `opmod`; `iopmod`; `iopcmod` |Fits the ordered probit model, the Zero-Inflated (ZIOP) & Middle-Inflated ordered probit (MIOP) models without correlated errors, and the ZIOPC & MIOPC models that incorporate correlated errors. |
|`opresults`; `iopresults`; `iopcresults`| Presents covariate estimates, Variance-Covariance (VCV) matrix, and goodness-of-fit statistics (Log-Likelihood and AIC) of `opmod`, `iopmod`, `iopcmod`.|
| `iopfit`; `iopcfit`| Computes fitted probabilities from each estimated model's object.|
| `vuong_opiop`; `vuong_opiopc` | Calculates Vuong test statistic for comparing the OP model's performance to ZiOP(C) and MiOP(C) models.|
|`split_effects`; `ordered_effects`| Estimates marginal effects of covariates from the split and outcome-stage respectively.|
|`mnlmod`; `gimnlmod`| Fits MNL model and Generalized-Inflated MNL models.|
|`mnlresults`; `gimnlresults`; `vuong_gimnl`| Presents covariate estimates, VCV matrix, and goodness-of-fit statistics of `mnlmod`, `gimnlmod`. Vuong test statistic for comparing MNL to GIMNL models obtained from `vuong_gimnl`|

## Dependencies
- scipy
- numpy
- pandas

## Installation

From [PyPi](https://pypi.org/project/idcempy/):

```sh
pip install idcempy
```

From [GitHub](https://github.com/hknd23/idcempy/)

```sh
git clone https://github.com/hknd23/idcempy.git
cd idcempy
python setup.py install
```
On [readthedocs](https://idcempy.readthedocs.io/en/latest/) you will find the installation guide, a complete overview of each feature included in **IDCeMPy**, and example scripts of all the models.

## Using the Package

### Example 1: Zero-inflated Ordered Probit Model with Correlated Errors (ZIOPC)
We illustrate how **IDCeMPy** can be used to estimate the OP and ZIOP(C) models for zero-inflated ordered outcome variables by using the CDC's 2018 [National Youth Tobacco Dataset](https://www.cdc.gov/tobacco/data_statistics/surveys/nyts/index.htm). The self-reported ordered tobacco consumption outcome variable in this data ranges from 0 for "no smoking" to 4 for "15 or more cigarettes". The zero "no smoking" category contains excessive observations that include permanent nonsmokers who never smoke (non-inflated cases) and transient nonsmokers (inflated cases) who temporarily stopped smoking because of high cigarette prices.   

**IDCeMPy** allows users to fit the ordered probit (OP) and Zero-inflated Ordered Probit (ZIOP) model without and with correlated errors (ZIOPC). The application of the OP model (available from `opmod`) and ZIOP model without correlated errors (see `iopmod`) to the CDC's 2018 Tobacco Consumption data is provided in the package's documentation. We fit the Zero-Inflated Ordered Probit Model with correlated errors to this data below.

First, install pandas and matplotlib to import and visualize data (if the packages are not already installed):

```sh
pip install pandas
pip install matplotlib
pip install urllib
```

Then, import `IDCeMPy`, required packages, and dataset.

```python
from idcempy import zmiopc
import pandas as pd
import matplotlib.pyplot as plot
import urllib
url = 'https://github.com/hknd23/idcempy/raw/main/data/tobacco_cons.csv'
data = pd.read_csv(url)
```
Users can define the lists with the names of the covariates to include in the ZIOPC model's split-stage (**Z**), the OP outcome-stage (**X**) as well as the zero-inflated ordered outcome variable (**Y**).

```python
X = ['age', 'grade', 'gender_dum']
Y = ['cig_count']
Z = ['gender_dum']
```

The default value of the starting parameters is set to .01. Users can, however, define an array of starting parameters before estimating the `ziopc` model and add it as an argument in the `iopcmod` function. 

The following line of code creates a ziopc regression object model.

```python
ziopc_tob = zmiopc.iopcmod('ziopc', data, X, Y, Z, method='bfgs',
                    weights=1, offsetx=0, offsetz=0)
```
Users can estimate the ZIOP model without correlated errors by using `zmiopc.iopmod` and the parameter 'ziop'. Please note that the models with correlated errors estimated with `zmiopc.iopcmod` have substantially higher run-time than `zmiopc.iopmod`. The above model takes roughly 1 hour and 8 minutes (on Windows 10, Intel Core i7-2600, 16GB RAM). 

The results from the ZIOPC model for this application are stored in a class (`ZiopcModel`) with the following attributes:

* *coefs*: Model coefficients and standard errors
* *llik*: Log-likelihood
* *AIC*: Akaike information criterion
* *vcov*: Variance-covariance matrix

We can generate the covariate estimates, standard errors, *p* value and *t* statistics in the ZIOPC case by typing:

```python
print(ziopc_tob.coefs)
```

```python
                          Coef        SE     tscore             p       2.5%      97.5%

Probit Split-stage
-----------------------
intercept              9.538072  3.470689   2.748178  5.992748e-03   2.735521  16.340623
gender_dum            -9.165963  3.420056  -2.680062  7.360844e-03 -15.869273  -2.462654

OP Outcome-stage
-----------------------
age                   -0.028606  0.008883  -3.220369  1.280255e-03  -0.046016  -0.011196
grade                  0.177541  0.010165  17.465452  0.000000e+00   0.157617   0.197465
gender_dum             0.602136  0.053084  11.343020  0.000000e+00   0.498091   0.706182
cut1                   1.696160  0.044726  37.923584  0.000000e+00   1.608497   1.783822
cut2                  -0.758095  0.033462 -22.655678  0.000000e+00  -0.823679  -0.692510
cut3                  -1.812077  0.060133 -30.134441  0.000000e+00  -1.929938  -1.694217
cut4                  -0.705836  0.041432 -17.036110  0.000000e+00  -0.787043  -0.624630
rho                   -0.415770  0.074105  -5.610526  2.017123e-08  -0.561017  -0.270524
```

The Akaike Information Criterion (AIC) statistics for the ZIOPC model is given by,

```python
print(ziopc_tob.AIC)
```

```python
16061.716497590078
```
The AIC of the OP and ZIOP models reported in the documentation is 8837.44 and 10138.32, respectively.

`split_effects` creates a dataframe that provides values to illustrate via boxplots (with 95% Confidence Intervals) the marginal effect of the ZIOP(C) model's split-stage covariates on the first difference in the predicted probability that the zero-category observations are non-inflated. In the tobacco consumption example,`split_effects` provides and illustrates via boxplots (with 95% CIs) the first difference in the predicted probability of zero-category observations being permanent nonsmokers (non-inflated cases) when the dummy split-stage covariate 'gender_dum' changes from 0 (female) to 1 (male).

```python
ziopcgender_split = zmiopc.split_effects(ziopc_tob, 1)
ziopcgender_split.plot.box(grid='False')
plot.show()
```

<p align="center">
   <img src="https://github.com/hknd23/idcempy/raw/main/graphics/ziopc_split_gender.png" width="500" height="300" />
   <br>
   <em>Fig. 1: Marginal Effect of Gender on Probability of Permanent Nonsmoker</em>
</p>

`ordered_effects`creates a dataframe that provides values to illustrate the marginal effect of the ZIOP(C) model's outcome-stage covariates on the first difference in the predicted probability of each ordered outcome category conditional on the zero-category observations being non-inflated. In the example below, `ordered_effects`provides and illustrate via boxplots (with 95% CIs) the first difference in the predicted probability (with 95% CIs) of each 0 to 4 ordered category of the tobacco consumption outcome when the dummy outcome-stage covariate 'gender_dum' changes from 0 to 1, conditional on zero-category observations being non-inflated.   

```python
ziopcgender_ordered = zmiopc.ordered_effects(ziopc_tob, 2)
ziopcgender_ordered.plot.box(grid='False')
plot.show()
```

<p align="center">
   <img src="https://github.com/hknd23/idcempy/raw/main/graphics/ziopc_ordered_gender_0214.png" width="500" height="300" />
   <br>
   <em>Fig. 2: Marginal Effect of Gender on Self-Reported Tobacco Consumption</em>
</p>

Module `zmiopc` also provides the function `vuong_opiopc` that employs the Vuong test statistic to compare the performance of the standard OP model (also available through `opmod`) versus the ZIOPC model and also the OP versus ZIOP model. The Vuong statistics from comparing the OP and the ZIOPC model is given by,

```python
op_tob = zmiopc.opmod(data, X, Y)
zmiopc.vuong_opiopc(op_tob, ziopc_tob)
```

```python
6.576246015382724
```
The Vuong test statistic favors the OP over both the ZIOPC model and ZIOP model (see documentation).

### Example 2: Middle-inflated Ordered Probit Models with Correlated Errors (MIOPC)
We next illustrate how **IDCeMPy** can be employed to fit the OP and MIOP(C) models for inflated middle-category ordered outcome variables. This is done by using Elgün and Tillman's ([2007](https://journals.sagepub.com/doi/10.1177/1065912907305684)) survey-response data in which the ordered outcome measure of support for the European Union (EU) by Europeans is given by 1 for “a bad thing,” 2 for “neither good nor bad,” and 3 for “a good thing.” The middle (neither good nor bad) category in this ordered measure contains excessive observations that include informed respondents who opt for this category based on their knowledge about the EU and uninformed respondents who choose this category to save face.

**IDCeMPy** allows users to fit the OP and Middle-inflated Ordered Probit (MIOP) model without and with correlated errors (MIOPC). The application of the OP model from `opmod` and MIOP model without correlated errors from `iopmod` to the EU support data is provided in the package's documentation. Users can estimate the MIOP model without correlated errors by simply substituting 'miop' for 'miopc'.

We turn to fit the Middle-Inflated Ordered Probit Model with correlated errors (MIOPC) to the aforementioned data. To this end, first load the dataset.

```python
url = 'https://github.com/hknd23/idcempy/raw/main/data/EUKnowledge.dta'
data = pd.read_stata(url)
```

Users can define the lists with names of the covariates they would like to include in the MIOPC model's split-stage (**Z**) and the second-stage (**X**) as well as the name of the ordered "middle-inflated" outcome variable (**Y**).

```python
Y = ["EU_support_ET"]
X = ['Xenophobia', 'discuss_politics']
Z = ['discuss_politics', 'EU_Know_obj']
```

Run the model and print the results (this MiOPC specification takes roughly 31 minutes to finish):

```python
miopc_EU = zmiopc.iopcmod('miopc', data, X, Y, Z)
```

```python
print(miopc_EU.coefs)

                              Coef    SE  tscore     p   2.5%  97.5%
Probit Split-stage
---------------------------
int                         -0.129 0.021  -6.188 0.000 -0.170 -0.088
discuss_politics             0.192 0.026   7.459 0.000  0.142  0.243
EU_Know_obj                  0.194 0.027   7.154 0.000  0.141  0.248

OP Outcome-stage
---------------------------
Xenophobia                  -0.591 0.045 -13.136 0.000 -0.679 -0.502
discuss_politics            -0.029 0.021  -1.398 0.162 -0.070  0.012
cut1                        -1.370 0.044 -30.948 0.000 -1.456 -1.283
cut2                        -0.322 0.103  -3.123 0.002 -0.524 -0.120
rho                         -0.707 0.106  -6.694 0.000 -0.914 -0.500
```

The AIC statistic for the MIOPC model is obtained from,

```python
print(miopc_EU.AIC)
```

```python
21669.96812802041
```

The AIC statistics for the MIOP model is 21729.39 and the OP model is 22100.90 (see documentation). 

In this EU support example, the `split_effects` dataframe provides and illustrates via boxplots (with 95% CI) the first difference in the predicted probability of middle-category observations being informed respondents (non-inflated cases) when the split-stage covariate 'EU_know_obj' is increased by one standard deviation from its mean value (for continuous variables, the "=0" and "=1" box plots represents the mean and one standard deviation above mean value, respectively).

<p align="center">
   <img src="https://github.com/hknd23/idcempy/raw/main/graphics/MiOPC_Split_EU_Know_0214.png" width="500" height="300" />
   <br>
   <em>Fig. 3: Marginal Effect of EU Knowledge on Probability of Informed Respondents</em>
</p>

`ordered_effects()` calculates and illustrates via boxplots (with 95% CI) the first difference in predicted probabilities of each ordered outcome category of "EU Support" when the outcome-stage Xenophobia covariate is increased by 1 standard deviation from its mean value, conditional on middle-category observations being informed respondents.

```python
xeno = zmiopc.ordered_effects(miopc_EU, 2)
xeno.plot.box(grid='False')
```

<p align="center">
   <img src="https://github.com/hknd23/idcempy/raw/main/graphics/MiOPC_EU_Xenophobia_0214.png" width="500" height="300" />
   <br>
   <em>Fig. 4: Marginal Effect of Xenophobia on EU Support</em>
</p>

Users can call the function `vuong_opiopc` to employ the Vuong test stastic to compare the OP model to the MIOPC model and also the OP to the MIOP model. The Vuong test statistics from comparing the OP to the MIOPC model is,  

```python
op_EU = zmiopc.opmod(DAT, X, Y)
zmiopc.vuong_opiopc(op_EU, miopc_EU)
```

```python
-10.435718518003675
```
The Vuong test statistic thus favors the MIOPC over the OP model, and also the MIOP over the OP model (see documentation).

### Example 3: Generalized Inflated Multinomial Logit Models (GIMNL)
**IDCeMPy** also includes functions to fit the GIMNL and standard MNL models. The Generalized Inflated Multinomial Logit Models account for the inflated and thus heterogeneous share of observations that can exist in the baseline or any other category of unordered polytomous outcome variables. To save space, we focus on just presenting the "Baseline" Inflated MNL (i.e., BIMNL) model that addresses excessive observations in the baseline category of unordered outcome measures. We fit this BIMNL model to the 2004 Presidential vote choice data from [Campbell and Monson (2008)](https://academic.oup.com/poq/article-abstract/72/3/399/1836972). The 0,1,2 unordered-polytomous Presidential vote choice dependent variable in their data includes the following options: abstained (their MNL baseline category), Bush, or Kerry. The inflated baseline category incorporates excessive observations of abstained nonvoters who did not vote in the said elections due to temporary factors and routine nonvoters who never vote.   

Users can fit the standard MNL model(available from `mnlmod`) to the Campbell and Monson (2008) data, which is described in the documentation. To illustrate how users can fit the BIMNL model to this data, however, we begin by importing the `gimnl` module.

```python
from idcempy import gimnl
url= 'https://github.com/hknd23/idcempy/raw/main/data/replicationdata.dta'
data= pd.read_stata(url)
```
Define the unordered vote choice outcome variable in the BIMNL as **Y**, whose unordered categories are given by 0,1,2. Denote the covariates in this model's logit split-stage as **Z** and **X** for the MNL-outcome stage for each unordered category 1 and 2.  

```python
x = ['educ', 'party7', 'agegroup2']
z = ['educ', 'agegroup2']
y = ['vote_turn']
```

```python
reference = [0, 1, 2]
inflatecat = "baseline"
```
The argument `inflatecat` can be used to specify any unordered category as the inflated category in their unordered-polytomous outcome measure. Further, from the argument `reference`, users can select which category of the unordered outcome variable is the baseline ("reference") category by placing it first. Since the baseline ("0") category in the Presidential vote choice outcome measure is inflated, the following code fits the BIMNL Model,

```python
gimnl_2004vote = gimnl.gimnlmod(data, x, y, z, reference, inflatecat)
```

Print the estimates:

```python
                       Coef    SE  tscore     p    2.5%  97.5%
Logit Split-stage
----------------------
intercept            -4.935 2.777  -1.777 0.076 -10.379  0.508
educ                  1.886 0.293   6.441 0.000   1.312  2.460
agegroup2             1.295 0.768   1.685 0.092  -0.211  2.800

MNL Outcome Category 1
---------------------
intercept            -4.180 1.636  -2.556 0.011  -7.387 -0.974
educ                  0.334 0.185   1.803 0.071  -0.029  0.697
party7                0.454 0.057   7.994 0.000   0.343  0.566
agegroup2             0.954 0.248   3.842 0.000   0.467  1.441

MNL Outcome Category 2
----------------------
intercept             0.900 1.564   0.576 0.565  -2.166  3.966
educ                  0.157 0.203   0.772 0.440  -0.241  0.554
party7               -0.577 0.058  -9.928 0.000  -0.691 -0.463
agegroup2             0.916 0.235   3.905 0.000   0.456  1.376
```
The AIC statistic for the BIMNL model is given by,

```python
print(gimnl_2004vote.AIC)
```
```
1656.8324085039708
```

The AIC for the standard MNL model (see documentation) is 1657.19. The Vuong statistic for comparing the MNL to the BIMNL model in this case is, 

```python
mnl_2004vote = gimnl.mnlmod(data, x, y, reference)
gimnl.vuong_gimnl(mnl_2004vote, gimnl_2004vote)
```

```python
-1.2835338187781173
```

Users can employ the argument `inflatecat` to specify any unordered category as the inflated category (dictated by the distribution) in their unordered-polytomous outcome measure. If a higher category (say 1 or 2) is inflated in the 0,1,2 unordered-polytomous outcome measure, then users can specify `reference` and `inflatecat` as follows,
```python
gimnl.gimnlmod(data, x, y, z, reference, inflatecat = "second")
gimnl.gimnlmod(data, x, y, z, reference, inflatecat = "third")
```
## Contributions

The authors welcome and encourage new contributors to help test `IDCeMPy` and add new functionality. Issues can be raised by any users for questions and bug reports. For further details, see [Guidelines for Contributors](https://github.com/hknd23/idcempy/blob/main/CONTRIBUTING.md).
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
nkh8@psu.edu.
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
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to IDCeMPy

Welcome to `IDCeMPy` package, and thank you for considering contributing to this Python
package. Please note that you can contribute in many ways, including fixing a bug, requesting a new
feature, code patches, editing the documentation, and patch reviews. We appreciate your help in
improving any feature included in this package.

This contributing guide contains the necessary information you need to know to help build and
improve `IDCeMPy`. Please read and follow these guidelines as it will make both communication
and the contribution process easy and effective. We will reciprocate by fixing bugs, evaluating changes, and helping you finalize your pull requests.

## Code of Conduct

We treat the open-source community with respect and hold ourselves as well as other
contributors to high standards of communication. By contributing to this project, you agree to
uphold our following [Code of Conduct](https://github.com/hknd23/idcempy/blob/main/CODE_OF_CONDUCT.md):

- Focusing on what is best for the community.
- Respect differing viewpoints and accept constructive criticisms.
- Avoid conduct that could reasonably be considered inappropriate in a professional setting.

Our responsibilities as maintainers include the right to remove, edit, or reject comments,
commits, code, issues, and other contributions that are not aligned to this [Code of Conduct](https://github.com/hknd23/idcempy/blob/main/CODE_OF_CONDUCT.md).

## Bug Reports

We use GitHub issues to track public bugs. Note that a bug can be reported by simply opening a
new issue. An issue template for reporting bugs is available. Before openning a new issue,
please check whether a similar issue has already been reported.
A useful Bug Report is one that includes:

- A summary or background.
- Steps to reproduce. Please be as specific as possible in this regard.
- Give sample code and screenshots if you can. This makes it easier to understand, track and correct the
main issues that you are raising.
- Briefly explain what you expected would happen and what actually transpired (contrary
to your expectations).
- Finally, please note the steps you may have undertaken to address the bug.

## Feature Requests

We welcome any requests to add new features. The feature request issue template is
available to facilitate the process. That said, please ensure that your request to add a new
feature, for example, fits the objectives and scope of the `IDCeMPy`.
Please provide as much detail as possible and ensure that the feature 
that you intend to add is compatible with other features included in this package. 

## Issues

Contributors should use Issues to report problems with the library, request a new feature, or to
discuss potential changes before a PR is created. When a contributor creates a new Issue, a
template will be loaded that will guide him or her to collect and provide the necessary
information that we (the maintainers) need to investigate. If the contributor finds an Issue that
addresses the problem they are having, then he or she should add their own reproduction
information to the existing issue rather than creating a new one.

## Pull Requests

Pull requests to our libraries in `IDCeMPy` are more than welcome. Pull requests (PRs) are the
best way to propose changes to the codebase. In general, we use [GitHub flow](https://guides.github.com/introduction/flow/index.html) for PRs to:

- [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) the repository to your
own Github account and create your branch from master.

- If you have added code that should be tested, please add tests.
- Commit changes to the branch.
- If you have changed APIs, please update the documentation.
- Follow any formatting and testing guidelines specific to this repo.
- Add unit or integration tests for fixed or changed functionality (if a test suite already
exists). Ensure that the test suite passes.
- Push changes to your fork.
- Open a PR in our repository and follow the PR template so that we can review ANY
changes.

**Please ask first** if you'd like to embark in any significant pull request. 

The precess below details the steps that you should follow if you'd like your work considered for inclusion in `IDCeMPy`:

1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/<your-username>/<repo-name>
   # Navigate to the newly cloned directory
   cd <repo-name>
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/<upstream-owner>/<repo-name>
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout <dev-branch>
   git pull upstream <dev-branch>
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit changes in logical chunks. Please adhere to these [git commit
   message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html) or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits before making them public.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream <dev-branch>
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description.

**IMPORTANT**: By submitting a patch, you agree to allow the project owner to
license your work under the same license as that used by the project.
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---
**Describe the bug**
A clear and concise description of what the bug is.

**Expected behavior**
A clear and concise description of what you expected to happen. 

**Actual behavior**
A clear and concise description of the output that you received. 

**To Reproduce**
Steps to reproduce the behavior:
1. Import '...'
2. Run function '....'
3. See error

Ideally, the steps to reproduce the bugs should reproduce exactly the actual behavior. If applicable, include URL to data source. 

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS, Windows]
 - IDE [e.g. PyCharm, IDLE]
 - Version [e.g. 3.7]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: feature request
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
Contact
=======

For more information, please contact nkh8@psu.edu. 


Package Documentation: IDCeMPy
==============================

Guide
^^^^^

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   idcempy_tutorial 
   api
   contact

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
***************
IDCeMPy Package
***************

Description
===========
`IDCeMPy` is a Python package that provides functions to fit and assess the performance of the following distinct
sets of “inflated” discrete choice models.

* Fit the Zero-Inflated Ordered Probit (ZIOP) model without and with correlated errors (ZIOPC model) to evaluate zero-inflated ordered choice outcomes that results from a dual data generating process (d.g.p.).

* Fit the Middle-Inflated Ordered Probit (MIOP) model without and with correlated errors (MIOPC) to account for the inflated middle-category in ordered choice measures that relates to a dual d.g.p.

* Fit Generalized Inflated Multinomial Logit (GIMNL) models that account for the preponderant and heterogeneous share of observations in the baseline or any lower category in unordered polytomous choice outcomes.

* Compute AIC and Log-likelihood statistics and the Vuong Test statistic to assess the performance of each inflated discrete choice model in the package.

`IDCeMPy` uses Newton numerical optimization methods to estimate the models listed above via Maximum Likelihood Estimation (MLE).

When should you use `IDCeMPy`?
==============================

An excessive (“inflated”) share of observations—stemming from two distinct d.g.p’s—fall into a single choice category in many ordered and unordered polytomous outcome variables. Standard Ordered Probit and Multinomial Logit models cannot account for such category inflation which leads to biased inferences. Examples include,

* The inflated zero-category of "no smoking" in ordered measures of self-reported smoking behavior is generated from nonsmokers who never smoke cigarettes and those who smoked previously but temporarily stopped smoking because of high cigarette prices.

* The inflated "indifference" middle-category in ordered measures of immigration attitudes includes respondents truly indifferent to immigration and those that choose indifference for social desirability reasons.

* The inflated baseline or other lower outcome categories of unordered polytomous outcome measures of vote choice include nonvoters who temporarily abstain from voting and routine nonvoters who always abstain.

`IDCeMPy` includes the ZIOP(C) models for evaluating zero-inflated ordered choice outcomes that results from a dual d.g.p, the MIOP(C) models that address inflated middle-category ordered outcome measures arising from distinct d.g.p’s, and GIMNL models that account for inflated baseline or other categories for unordered polytomous outcomes.

Each inflated discrete choice model in this package addresses category inflation in one’s discrete outcome—unordered or unordered polytomous—of interest by jointly estimating a binary split-stage equation and an ordered or multinomial discrete choice outcome equation.

Installation
=============
The package can be installed in two different ways.

1. From `PyPi <https://pypi.org/>`__:

.. testcode::

  # Import the package

  pip install idcempy

2. From its `GitHub Repository <https://github.com/hknd23/idcempy/>`__:

.. testcode::

  # Import the package

  git clone https://github.com/hknd23/idcempy.git
  cd idcempy
  python setup.py install

Examples
========

Zero-inflated Ordered Probit (ZiOP) Model without Correlated Errors
--------------------------------------------------------------------
The `iopcmod` function estimates regression objects for "zero-inflated" and "middle-inflated" ordered probit models without correlated errors.  Below you will find instructions to estimate a ZiOP model.

We first import the required libraries, set up the package and import the dataset:

.. testcode::

  # Import the necessary libraries and package

  import numpy as np
  import pandas as pd
  import urllib
  from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset.

  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'

  # Read the dataset
  data=pd.read_csv(url)

Our data is now a `pandas` DataFrame, and we can proceed to estimate the ZiOP model as follows.

.. testcode::

  # First, you should define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `Data.Frame) used in ordered probit stage.
  # Z = Column names of covariates (from `Data.Frame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `Data.Frame`).

  X = ['age', 'grade', 'gender_dum']
  Z = ['gender_dum']
  Y = ['cig_count']

The package sets a default start value of .01 for all parameters.  Users can modify it by creating an array with their desired values, define such array as `pstart` and add it to as an argument in the model function.  

:func:`zmiopc.iopmod` estimates the ZiOP model and returns :class:`zmiopc.IopModel`.

.. testcode::

   # Model estimation:
  ziop_tob= zmiopc.iopmod('ziop', data, X, Y, Z, method='bfgs', weights= 1,offsetx= 0, offsetz=0)

Results from the model:

The following message will appear when the model has converged:

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 5060.160903
         Iterations: 79
         Function evaluations: 1000
         Gradient evaluations: 100

Object :class:`zmiopc.IopModel` stores model results and goodness-of-fit tests in its attributes 'coefs', 'AIC', 'llik', and 'vcov'.

Use the following line of code to see the estimates of coefficients:

.. testcode::

   print(ziop_tob.coefs)

.. testoutput::

                            Coef        SE      tscore        p           2.5%      97.5%
   cut1                   1.693797  0.054383  31.145912  0.000000e+00   1.587207   1.800387
   cut2                  -0.757830  0.032290 -23.469359  0.000000e+00  -0.821119  -0.694542
   cut3                  -1.804483  0.071237 -25.330846  0.000000e+00  -1.944107  -1.664860
   cut4                  -0.691907  0.052484 -13.183210  0.000000e+00  -0.794775  -0.589038
   Inflation: int         4.161455  3.864721   1.076780  2.815784e-01  -3.413398  11.736309
   Inflation: gender_dum -3.462848  3.857160  -0.897772  3.693074e-01 -11.022881   4.097185
   Ordered: age          -0.029139  0.013290  -2.192508  2.834282e-02  -0.055187  -0.003090
   Ordered: grade         0.177897  0.012133  14.661952  0.000000e+00   0.154116   0.201678
   Ordered: gender_dum    0.206509  0.034914   5.914823  3.322323e-09   0.138078   0.274940

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).  You can obtain them via the following commands:

.. testcode::

  print(ziop_tob.llik)
  print(ziop_tob.AIC)
  print(ziop_tob.vcov)

An example for the AIC:

.. testcode::

  print(ziop_tob.AIC)

.. testoutput::

  10138.321806674261

You can also extract predicted probabilities from the model:
:func:`zmiopc.iopfit` returns :class:`zmiopc.FittedVals` containing fitted probablities.

.. testcode::

  fittedziop = ziopc.iopfit(ziop_tob)
  print(fittedziopc.responsefull)

.. testoutput::

  array[[0.8822262  0.06879832 0.01455244 0.0242539  0.01016914]
 [0.84619828 0.08041296 0.01916279 0.03549797 0.01872801]
 [0.93105632 0.04349743 0.00831396 0.0127043  0.004428  ]
 ...
 [0.73347708 0.1291157  0.03295816 0.06500889 0.03944016]
 [0.87603805 0.06808193 0.01543795 0.02735256 0.01308951]
 [0.82681957 0.08778215 0.02153509 0.04095753 0.02290566]]

You can compute changes in predicted probabilities when the value of a variable changes.
This allows you to illustrate how changes in the split-probit covariates affect the probabilities of
being in one population versus another. The example below illustrates the marginal effects of the variable
'gender_dum' on the outcome variable in the ZiOP model estimated above.

.. testcode::

    ziopcgender = idcempy.split_effects(ziop_tob, 1, nsims = 10000)

The returned dataframe contains predicted probabilities when 'gender_dum' equals 0, and when 'gender_dum' equals 1.

You can also calculate the change in predicted probabilities of the outcome variable when the value of a covarariate changes, and plot those values.

.. testcode::

    gender = zmiopc.ordered_effects(ziop_tob, 2, nsims = 10000)
    gender.plot.box(grid='False')

Zero-inflated Ordered Probit (ZiOPC) with Correlated Errors
-----------------------------------------------------------
The package also includes the function `iopcmod` which fits "zero-inflated" ordered probit models (ZiOPC) under the assumption that the two errors are correlated with each other (i.e. correlated errors).

We first import the required libraries, set up the package and import the dataset:

.. testcode::
  # Import the necessary libraries and IDCeMPy.

  import numpy as np
  import pandas as pd
  import urllib
  from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset.

  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'

  # Read the imported dataset.
  data=pd.read_stata(url)

.. testcode::

  # First, you should define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `Data.Frame) used in ordered probit stage.
  # Z = Column names of covariates (from `Data.Frame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `Data.Frame`).

  X = ['age', 'grade', 'gender_dum']
  Z = ['gender_dum']
  Y = ['cig_count']

Our data is now a `pandas` DataFrame, and we can proceed to estimate the ZiOP model as follows.

.. testcode::

    ziopc_tob = zmiopc.iopcmod('ziopc', data, X, Y, Z, method='bfgs', weights=1, offsetx=0, offsetz=0)

The package sets a default start value of .01 for all parameters.  Users can modify it by creating an array with their desired values, define such array as `pstart` and add it to as an argument in the model function.

The results are stored in the attributes of :class:`zmiopc.IopCModel`.

.. testoutput::

         Current function value: 5060.051910
         Iterations: 119
         Function evaluations: 1562
         Gradient evaluations: 142

The following line of code prints the results

.. testcode::

    print(ziopc_tob.coefs)

.. testoutput::

                            Coef        SE     tscore             p       2.5%      97.5%
   cut1                   1.696160  0.044726  37.923584  0.000000e+00   1.608497   1.783822
   cut2                  -0.758095  0.033462 -22.655678  0.000000e+00  -0.823679  -0.692510
   cut3                  -1.812077  0.060133 -30.134441  0.000000e+00  -1.929938  -1.694217
   cut4                  -0.705836  0.041432 -17.036110  0.000000e+00  -0.787043  -0.624630
   Inflation: int         9.538072  3.470689   2.748178  5.992748e-03   2.735521  16.340623
   Inflation: gender_dum -9.165963  3.420056  -2.680062  7.360844e-03 -15.869273  -2.462654
   Ordered: age          -0.028606  0.008883  -3.220369  1.280255e-03  -0.046016  -0.011196
   Ordered: grade         0.177541  0.010165  17.465452  0.000000e+00   0.157617   0.197465
   Ordered: gender_dum    0.602136  0.053084  11.343020  0.000000e+00   0.498091   0.706182
   rho                   -0.415770  0.074105  -5.610526  2.017123e-08  -0.561017  -0.270524

To print the estimates of the log-likelihood, AIC, and Variance-Covariance matrix, you should type:

.. testcode::

  print(ziopc_tob.llik)
  print(ziopc_tob.AIC)
  print(ziopc_tob.vcov)

The AIC of the ziopc_tob model, for example, is:

.. testoutput::

  10140.103819465658

The predicted probabilities from the `ziopc_tob` model can ve obtained as follows.

:func:`zmiopc.iopcfit` returns :class:`zmiopc.FittedVals` containing fitted probablities.

.. testcode::

  fittedziopc = zmiopc.iopcfit(ziopc_tob)
  print(fittedziopc.responsefull)

.. testoutput::

  array[[0.88223509 0.06878162 0.01445941 0.0241296  0.01039428]
 [0.84550989 0.08074461 0.01940226 0.03589458 0.01844865]
 [0.93110954 0.04346074 0.00825639 0.01264189 0.00453143]
 ...
 [0.73401588 0.12891071 0.03267436 0.06438928 0.04000977]
 [0.87523652 0.06888286 0.01564958 0.0275354  0.01269564]
 [0.82678185 0.0875059  0.02171135 0.04135142 0.02264948]]

You can compute changes in predicted probabilities when the value of a variable changes.
This allows you to illustrate how the changes in the split-probit covariates affect the probabilities of being in one population versus another. The example below illustrates the marginal effects of the variable 'gender_dum' on the outcome variable in the ZiOPC model estimated in ths documentation.

.. testcode::

    ziopcgender = idcempy.split_effects(ziopc_tob, 1, nsims = 10000)

You can calculate the change in predicted probabilities of the outcome variable when the value of a covarariate changes.
In addition, a you can obtain a box plot that displays the change in predicted probabilities of the outcome variable in the ZiOPC model.

.. testcode::

   # Calculate change in predicted probabilities
   gender = zmiopc.ordered_effects(ziopc_tob, 1, nsims = 10000)

   # Box-plot of precicted probablilites
   gender.plot.box(grid='False')

Middle-inflated Ordered Probit (MiOP) without Correlated Errors
---------------------------------------------------------------
If your ordered outcome variable is inflated in the middle category, you should estimate a  Middle-inflated Ordered Probit (MiOP) model.

The following example uses data from Elgun and Tilam (`2007 <https://journals.sagepub.com/doi/10.1177/1065912907305684>`_).

We begin by loading the required libraries and IDCeMPy

.. testcode::

  # Import the necessary libraries and IDCeMPy.

  import numpy as np
  import pandas as pd
  import urllib
  from idcempy import zmiopc

Next, we load the dataset.
.. testcode::

    # Import and read the dataset
    url = 'https://github.com/hknd23/zmiopc/blob/main/data/'
    data2 = pd_read.stata(url)

We then define the lists with the names of the variables used in the model
.. testcode::

  # First, you should define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `Data.Frame) used in ordered probit stage.
  # Z = Column names of covariates (from `Data.Frame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `Data.Frame`).

  X = ['Xenophobia', 'discuss_politics']
  Z = ['discuss_politics', 'EU_Know_ob']
  Y = ['EU_support_ET']

Your data is now ready, and you can begin the estimation process.

:func:`zmiopc.iopmod` estimates the MiOP model and returns :class:`zmiopc.IopModel`.

.. testcode::

  # Model estimation:
  miop_EU = zmiopc.iopmod('miop', data, X, Y, Z, method='bfgs', weights= 1,offsetx= 0, offsetz=0)

The following message will appear when the model finishes converging.

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 10857.695490
         Iterations: 37
         Function evaluations: 488
         Gradient evaluations: 61  # See estimates:

Print the results.

.. testcode::

   print(miop_EU.coefs)

.. testoutput::

                                 Coef        SE       tscore         p         2.5%     97.5%
   cut1                        -1.159621  0.049373 -23.487133  0.000000e+00 -1.256392 -1.062851
   cut2                        -0.352743  0.093084  -3.789492  1.509555e-04 -0.535188 -0.170297
   Inflation: int              -0.236710  0.079449  -2.979386  2.888270e-03 -0.392431 -0.080989
   Inflation: discuss_politics  0.190595  0.035918   5.306454  1.117784e-07  0.120197  0.260993
   Inflation: EU_Know_obj       0.199574  0.020308   9.827158  0.000000e+00  0.159770  0.239379
   Ordered: Xenophobia         -0.663551  0.044657 -14.858898  0.000000e+00 -0.751079 -0.576024
   Ordered: discuss_politics    0.023784  0.029365   0.809964  4.179609e-01 -0.033770  0.081339

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).

.. testcode::

   # Print estimates of LL, AIC and VCOV
   print(miop_EU.llik)
   print(miop_EU.AIC)
   print(miop_EU.vcov)

For example, the AIC in this case is:

.. testcode::

   print(miop_EU.AIC)

.. testoutput::

   21729.390980849777

To estimate the predicted probabilities:

.. testcode::

   fittedmiop = zmiopc.iopcfit(miop_EU)
   print(fittedziopc.responsefull)

The package also allows you to simulates data from MiOP model results and compute changes in predicted probabilities when the value of a variable changes.
This allows you to illustrate how the changes in the split-probit covariates affect the probablilities of being in one population versus another.

.. testcode::

    miopxeno = idcempy.split_effects(miop_EU, 1, nsims = 10000)

To plot the predicted probabilities.

.. testcode::

     miopxeno.plot.box(grid='False')


You can calculate the change in predicted probabilities of the outcome variable when the value of a covarariate changes. The box plots below display the change in predicted probabilities of the outcome variable in the MiOP model estimated above when Xenophobia increases one standard deviation from its mean value.

.. testcode::

    xeno = zmiopc.ordered_effects(miop_EU, 2, nsims = 10000)
    xeno.plot.box(grid='False')

Middle-inflated Ordered Probit (MiOPC) Model with Correlated Errors
--------------------------------------------------------------------
You can estimate a Middle-inflated Ordered Probit (MiOPC) with correlated errors as follows.

We begin by loading the required libraries and IDCeMPy

.. testcode::

  # Import the necessary libraries and IDCeMPy.

  import numpy as np
  import pandas as pd
  import urllib
  from idcempy import zmiopc

Next, we load the dataset.

.. testcode::

    # Import and read the dataset
    url = 'https://github.com/hknd23/zmiopc/blob/main/data/'
    data2 = pd_read.stata(url)

We then define the lists with the names of the variables used in the model

.. testcode::

   # First, you should define a list of variable names of X, Z, and Y.
   # X = Column names of covariates (from `Data.Frame) used in ordered probit stage.
   # Z = Column names of covariates (from `Data.Frame`) used in split-population stage.
   # Y = Column name of ordinal outcome variable (from `Data.Frame`).

   X = ['Xenophobia', 'discuss_politics']
   Z = ['discuss_politics', EU_Know_ob]
   Y = ['EU_support_ET']

The model can be estimated as follows.

:func:`zmiopc.iopcmod` estimates the MiOPC model and returns :class:`zmiopc.IopcModel`.

.. testcode::

   # Model estimation:
   miopc_EU = zmiopc.iopcmod('miopc', pstartziop, data, X, Y, Z, method='bfgs', weights= 1,offsetx= 0, offsetz=0)

Now print(miopc_EU.coefs).

.. testoutput::

                                 Coef  SE     tscore  p     2.5%  97.5%
   cut1                        -1.370 0.044 -30.948 0.000 -1.456 -1.283
   cut2                        -0.322 0.103  -3.123 0.002 -0.524 -0.120
   Inflation: int              -0.129 0.021  -6.188 0.000 -0.170 -0.088
   Inflation: discuss_politics  0.192 0.026   7.459 0.000  0.142  0.243
   Inflation: EU_Know_obj       0.194 0.027   7.154 0.000  0.141  0.248
   Ordered: Xenophobia         -0.591 0.045 -13.136 0.000 -0.679 -0.502
   Ordered: discuss_politics   -0.029 0.021  -1.398 0.162 -0.070  0.012
   rho                         -0.707 0.106  -6.694 0.000 -0.914 -0.500

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).  You can obtain them via the following commands:

.. testcode::

   print(miopc_EU.llik)
   print(miopc_EU.AIC)
   print(miopc_EU.vcov)

To estimate the predicted probabilities:

.. testcode::

   fittedmiopc = zmiopc.iopcfit(miopc_EU)
   print(fittedziopc.responsefull)

The following line of code allows you to compute changes in predicted probabilities when the value of a variable changes.
This allows you to illustrate how the changes in the split-probit covariates affect the probablilities of being in one population versus another.

.. testcode::

    miopcxeno = idcempy.split_effects(miopc_EU, 1, nsims = 10000)

A box plot can illustrate the change in predicted probabilities.

.. testcode::

     miopcxeno.plot.box(grid='False')


To calculate the change in predicted probabilities of the outcome variable when the value of a covarariate changes. The box plots below display the change in predicted probabilities of the outcome variable in the MiOPC model estimated above when Xenophobia increases one standard deviation from its mean value.

.. testcode::

    xeno = zmiopc.ordered_effects(miopc_EU, 2, nsims = 10000)
    xeno.plot.box(grid='False')


The Standard Ordered Probit (OP) model
--------------------------------------

The package also includes a function that estimates a standard Ordered Probit (OP) model.
The OP model does not account for the "zero inflation", so it does not have a split-probit stage.

We first import the required libraries, set up the package and import the dataset:

.. testcode::

   # Import the necessary libraries and package

   import numpy as np
   import pandas as pd
   import urllib
   from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset.

  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'

  # Read the dataset
  data=pd.read_csv(url)

.. testcode::

     # Define a list of variable names X,Y:
     # X = Column names of covariates (from `Data.Frame`) in the OP equation
     # Y = Column name of outcome variable (from `Data.Frame`).

     X = ['age', 'grade', 'gender_dum']
     Y = ['cig_count']

Your data is not ready for estimation.

.. testcode::

  # Starting parameters for optimization:
  pstartop = np.array([.01, .01, .01, .01, .01, .01, .01])

  # Model estimation:
  op_tob = zmiopc.opmod(pstartop, data, X, Y, method='bfgs', weights=1, offsetx=0)

  # See estimates:
  print(ziop_tob.coefs)

Results from the model:

The following message will appear when the model has converged:

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 4411.710049
         Iterations: 10
         Function evaluations: 976
         Gradient evaluations: 121

:class:`zmiopc.OpModel` stores results from model estimation and other information in its attributes.
The following line of code to see the estimates of coefficients:

.. testcode::

   print(op_tob.coefs)

.. testoutput::

                Coef        SE     tscore         p      2.5%     97.5%
   cut1        1.696175  0.047320  35.844532  0.000000  1.603427  1.788922
   cut2       -0.705037  0.031650 -22.276182  0.000000 -0.767071 -0.643004
   cut3       -2.304405  0.121410 -18.980329  0.000000 -2.542369 -2.066441
   cut4        2.197381  0.235338   9.337141  0.000000  1.736119  2.658643
   age        -0.070615  0.007581  -9.314701  0.000000 -0.085474 -0.055756
   grade       0.233741  0.010336  22.614440  0.000000  0.213483  0.254000
   gender_dum  0.020245  0.032263   0.627501  0.530331 -0.042991  0.083482

Log-likelihood, AIC, and Variance-Covariance matrix can be extracted with:

.. testcode::

  print(op_tob.llik)
  print(op_tob.AIC)
  print(op_tob.vcov)

The Vuong Test
--------------

Harris and Zhao (`2007 <https://doi.org/10.1016/j.jeconom.2007.01.002>`__) suggest that a variant of the Vuong (`1989 <https://www.jstor.org/stable/1912557>`__) Test (with a v statistic) can be used to compare the performance of the ZiOP versus the standard Ordered Probit (OP) model using :func:`zmiopc.vuong_opiop`.

The formula to estimate the Vuong test is: 

.. math::

    v = \frac{\sqrt{N}(\frac{1}{N}\sum_{i}^{N}m_{i})}{\sqrt{\frac{1}{N}\sum_{i}^{N}(m_{i}-\bar{m})^{2}}}

where v < -1.96 favors the more general (ZiOP/ZiOPC) model, -1.96 < v < 1.96 lends no support to either model, and v > 1.96 supports the simpler (OP) model.

The OP and ZiOP models must have the same number of observations, and the OP must have the same number of covariates as ZiOP's OP stage. The statistic reveals that the OP model is preferred over the ZiOP model.

.. testcode::

   zmiopc.vuong_opiop(op_tob, ziop_tob)

.. testoutput::

   6.624742132792222

The Vuong test can also be implemented to compare the ZiOPC, MiOP and MiOPC models and the OP model.

Generalized Inflated Multinomial Logit (GiMNL) Model
----------------------------------------------------

The IDCeMPy package also includes a function that estimates General "inflated" Multinomial Logit models (GiMNL). GiMNL models minimize issues present when unordered polytomous outcome variables have an excessive share and heterogeneous pool of observations in the lower category.
Failing to account for such inflation could lead to inaccurate inferences.

To estimate the GiMNL model, we first import the library and the dataset introduced above.

.. testcode::

   from idcempy import gimnl
   url= 'https://github.com/hknd23/zmiopc/raw/main/data/replicationdata.dta'
   data= pd.read_stata(url)

We the define the list of covariates in the split-stage (z), the second-stage (x) and the outcome variable (y).

.. testcode::

   # Define the variable list for x, y and z.
   # x = Column names of covariates (from `Data.Frame`) in the outcome-stage.
   # z = Column names of covariates (from `Data.Frame`) in the split-stage.
   # y = Column names of outcome variable (from `Data.Frame`).

   x = ['educ', 'party7', 'agegroup2']
   z = ['educ', 'agegroup2']
   y = ['vote_turn']

Users can employ the argument `inflatecat` to specify any unordered category as the inflated category (dictated by the distribution) in their unordered-polytomous outcome measure. If a higher category (say 1) is inflated in a 0,1,2 unordered outcome measure.

We first need to specify the order of the outcome variable. Then, you need to define which category is "inflated."

.. testcode::

   order = [0, 1, 2]
   inflatecat = "baseline"

Further, employing the argument `reference`, users can select which category of the unordered outcome variable is the baseline ("reference") category by placing it first. Since the baseline ("0") category in the Presidential vote choice outcome measure is inflated, the following code fits the BIMNL Model.

.. testcode::

   gimnl_2004vote = gimnl.gimnlmod(data, x, y, z, order, inflatecat)

The following line of code prints the coefficients of the covariates.

.. testcode::

   print(gimnl_2004vote.coefs)

.. testoutput::

                          Coef   SE    tscore   p    2.5%   97.5%
   Inflation: int       -4.935 2.777  -1.777 0.076 -10.379  0.508
   Inflation: educ       1.886 0.293   6.441 0.000   1.312  2.460
   Inflation: agegroup2  1.295 0.768   1.685 0.092  -0.211  2.800
   1: int               -4.180 1.636  -2.556 0.011  -7.387 -0.974
   1: educ               0.334 0.185   1.803 0.071  -0.029  0.697
   1: party7             0.454 0.057   7.994 0.000   0.343  0.566
   1: agegroup2          0.954 0.248   3.842 0.000   0.467  1.441
   2: int                0.900 1.564   0.576 0.565  -2.166  3.966
   2: educ               0.157 0.203   0.772 0.440  -0.241  0.554
   2: party7            -0.577 0.058  -9.928 0.000  -0.691 -0.463
   2: agegroup2          0.916 0.235   3.905 0.000   0.456  1.376

The results from the model are stored in a :class:`gimnlModel` with the following attributes:

- coefs: Model coefficients and standard errors
- llik: Log-likelihood
- AIC: Akaike information criterion
- vcov: Variance-covariance matrix

You can, for example, print the AIC as follows.

.. testcode::

    print(gimnl_2004vote.AIC)
    
.. testoutput::

    1656.8324085039708

Using the function :py:func:`gimnl.mnlmod`, users can fit a standard Multinomial Logit Model (MNL) by specifying the list of **X**, **Y**, and baseline (using `reference`).

.. testcode::

   mnl_2004vote = gimnl.mnlmod(data, x, y, z, order)
   print(mnl_2004vote.coefs)

.. testoutput::

     Coef        SE  tscore     p   2.5%  97.5%
  1: int       -4.914 0.164 -29.980 0.000 -5.235 -4.593
  1: educ       0.455 0.043  10.542 0.000  0.371  0.540
  1: party7     0.462 0.083   5.571 0.000  0.300  0.625
  1: agegroup2  0.951 0.029  32.769 0.000  0.894  1.008
  2: int        0.172 0.082   2.092 0.036  0.011  0.334
  2: educ       0.282 0.031   9.011 0.000  0.221  0.343
  2: party7    -0.567 0.085  -6.641 0.000 -0.734 -0.399
  2: agegroup2  0.899 0.138   6.514 0.000  0.629  1.170

Similar to the GiMNL model, the AIC for the MNL model can also be given by:

.. testcode::

    print(mnl_2004vote.AIC)

.. testoutput::

    1657.192925769978

***************
IDCeMPy Package
***************

Description
===========
`IDCeMPy` is a Python package that provides functions to fit and assess the performance of the following distinct
sets of “inflated” discrete choice models:

* Fit the Zero-Inflated Ordered Probit (ZiOP) model without and with correlated errors (ZiOPC model) to evaluate zero-inflated ordered choice outcomes that results from a dual data generating process (d.g.p.).

* Fit the Middle-Inflated Ordered Probit (MiOP) model without and with correlated errors (MiOPC) to account for the inflated middle-category in ordered choice measures that relates to a dual d.g.p.

* Fit Generalized Inflated Multinomial Logit (GiMNL) models that account for the preponderant and heterogeneous share of observations in the baseline or any lower category in unordered polytomous choice outcomes.

* Compute AIC and Log-likelihood statistics and the Vuong Test statistic to assess the performance of each inflated discrete choice model in the package.

`IDCeMPy` uses Newton numerical optimization methods to estimate the models listed above via Maximum Likelihood Estimation (MLE).

When Should You use `IDCeMPy`?
==============================

An excessive (“inflated”) share of observations—stemming from two distinct d.g.p’s—fall into a single choice category in many ordered and unordered polytomous outcome variables.
Standard Ordered Probit and Multinomial Logit models cannot account for such category inflation which leads to biased inferences. Examples for such d.g.p’s include:

* The inflated zero-category of "no smoking" in ordered measures of self-reported smoking behavior is generated from nonsmokers who never smoke cigarettes and those who smoked previously but temporarily stopped smoking because of high cigarette prices.

* The inflated "indifference" middle-category in ordered measures of immigration attitudes includes respondents truly indifferent to immigration and those that choose indifference for social desirability reasons.

* The inflated baseline or other lower outcome categories of unordered polytomous outcome measures of vote choice include nonvoters who temporarily abstain from voting and routine nonvoters who always abstain.

`IDCeMPy` includes the ZIOP(C) models for evaluating zero-inflated ordered choice outcomes that results from a dual d.g.p, the MIOP(C) models that address inflated middle-category ordered outcome measures arising from distinct d.g.p’s, and GIMNL models that account for inflated baseline or other categories for unordered polytomous outcomes.

Each inflated discrete choice model in this package addresses category inflation in one’s discrete outcome—unordered or unordered polytomous—of interest by jointly estimating a binary split-stage equation and an ordered or multinomial discrete choice outcome equation.

Installation
============
The package can be installed in two different ways:

1. From `PyPi <https://pypi.org/project/idcempy/>`__:

.. testcode::

  $  pip install idcempy

2. From its `GitHub Repository <https://github.com/hknd23/idcempy/>`__:

.. testcode::

  $  git clone https://github.com/hknd23/idcempy.git
  $  cd idcempy
  $  python setup.py install

Examples
========

The examples below demonstrate how to use the `IDCeMPy` package to estimate the inflated discrete choice models ZiOP(C), MiOP(C), and GiMNL.
The example code files, with rough calculation of model run time, are available in the `/examples <https://github.com/hknd23/idcempy/tree/main/examples>`__ directory.
For each model example below, the run time is available as reference point. The specification used to record the times is Intel Core i7-2600 (3.40GHz Quad core), 16GB RAM.
Please note that for models in the `zmiopc` module, the run-time for models with correlated errors estimated with :func:`zmiopc.iopcmod` is substantially higher
than their without correlated errors counterparts using :func:`zmiopc.iopmod`. Other factors affecting run-time are the number of observations and the number of covariates.

The examples use the pandas, urllib, and matplotlib packages for importing and visualizing data:

.. testcode::

  $  pip install pandas
  $  pip install matplotlib
  $  pip install urllib

Zero-inflated Ordered Probit (ZiOP) Model without Correlated Errors
-------------------------------------------------------------------
The :func:`zmiopc.iopmod` function estimates regression objects for "zero-inflated" and "middle-inflated" ordered probit models without correlated errors.
This section provides instruction to estimate the ZiOP model using the self-reported smoking behavior as empirical example.

We first import the required libraries, set up the package and import the dataset:

.. testcode::

  # Import the necessary libraries and package
  import pandas as pd
  import urllib
  import matplotlib.pyplot as plot
  from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset as a pandas.DataFrame
  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'
  data = pd.read_csv(url)

The data is now a `pandas` DataFrame, and we can proceed to estimate the ZiOP model as follows.

.. testcode::

  # First, define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `DataFrame`) used in ordered probit stage.
  # Z = Column names of covariates (from `DataFrame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `DataFrame`).

  X = ['age', 'grade', 'gender_dum']
  Z = ['gender_dum']
  Y = ['cig_count']

The package sets a default start value of .01 for all parameters.
Users can specify their own starting parameters by creating a list or numpy.array with their desired values.

:func:`zmiopc.iopmod` estimates the ZiOP model and returns :class:`zmiopc.IopModel`.

.. testcode::

   # Model estimation:
   ziop_tob= zmiopc.iopmod('ziop', data, X, Y, Z, method = 'bfgs', weights = 1, offsetx = 0, offsetz = 0)

   # 'ziop' = model to be estimated. In this case 'ziop'
   # data = name of Pandas DataFrame
   # X = variables in the ordered probit stage.
   # Y = dependent variable.
   # Z = variables in the inflation stage.
   # method = method for optimization.  By default set to 'bfgs'
   # weights = weights.
   # offsetx = offset of X.  By Default is zero.
   # offsetz = offset of z


Results from the model:

The following message will appear when the model has converged:

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 5060.160903
         Iterations: 79
         Function evaluations: 1000
         Gradient evaluations: 100

The run-time for this model is 80.006 seconds (N= 9624).
Object :class:`zmiopc.IopModel` stores model results and goodness-of-fit tests in its attributes 'coefs', 'AIC', 'llik', and 'vcov'.

The following line of code prints the estimates of coefficients:

.. testcode::

   print(ziop_tob.coefs)

.. testoutput::

                            Coef        SE      tscore        p           2.5%      97.5%
   cut1                   1.693797  0.054383  31.145912  0.000000e+00   1.587207   1.800387
   cut2                  -0.757830  0.032290 -23.469359  0.000000e+00  -0.821119  -0.694542
   cut3                  -1.804483  0.071237 -25.330846  0.000000e+00  -1.944107  -1.664860
   cut4                  -0.691907  0.052484 -13.183210  0.000000e+00  -0.794775  -0.589038
   Inflation: int         4.161455  3.864721   1.076780  2.815784e-01  -3.413398  11.736309
   Inflation: gender_dum -3.462848  3.857160  -0.897772  3.693074e-01 -11.022881   4.097185
   Ordered: age          -0.029139  0.013290  -2.192508  2.834282e-02  -0.055187  -0.003090
   Ordered: grade         0.177897  0.012133  14.661952  0.000000e+00   0.154116   0.201678
   Ordered: gender_dum    0.206509  0.034914   5.914823  3.322323e-09   0.138078   0.274940

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object :class:`zmiopc.IopModel` also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).
They can be obtained via the following:

.. testcode::

  print(ziop_tob.llik)
  print(ziop_tob.AIC)
  print(ziop_tob.vcov)

An example for the AIC:

.. testcode::

  print(ziop_tob.AIC)

.. testoutput::

  10138.321806674261

The following funtion extracts predicted probabilities from the model:
:func:`zmiopc.iopfit` returns :class:`zmiopc.FittedVals` containing fitted probablities.

.. testcode::

  fittedziop = ziopc.iopfit(ziop_tob)

  # Print the predicted probabilities
  print(fittedziopc.responsefull)

.. testoutput::

  array[[0.8822262  0.06879832 0.01455244 0.0242539  0.01016914]
 [0.84619828 0.08041296 0.01916279 0.03549797 0.01872801]
 [0.93105632 0.04349743 0.00831396 0.0127043  0.004428  ]
 ...
 [0.73347708 0.1291157  0.03295816 0.06500889 0.03944016]
 [0.87603805 0.06808193 0.01543795 0.02735256 0.01308951]
 [0.82681957 0.08778215 0.02153509 0.04095753 0.02290566]]

:func:`zmiopc.split_effects` and :func:`zmiopc.ordered_effects` compute changes in predicted probabilities when the value of a variable changes in the Inflation or Ordered stages, respectively.

:func:`zmiopc.split_effects` computes how changes in the split-probit covariates affect the probabilities of
being in one population versus another. The example below illustrates the marginal effects of the variable
'gender_dum' on the outcome variable in the ZiOP model estimated above.

.. testcode::

    ziopcgender = zmiopc.split_effects(ziop_tob, 1, nsims = 10000)

The returned dataframe contains predicted probabilities when 'gender_dum' equals 0, and when 'gender_dum' equals 1.

Likewise, :func:`zmiopc.ordered_effects` can also calculate the change in predicted probabilities in each of the ordered outcomes in the ordered-probit stage when the value of a covarariate changes.
Results from :func:`zmiopc.split_effects` and :func:`zmiopc.ordered_effects` can be illustrated using `matplotlib` box plots:

.. testcode::

    gender = zmiopc.ordered_effects(ziop_tob, 2, nsims = 10000)

    # The box plot from the results:
    gender.plot.box(grid='False')

Zero-inflated Ordered Probit (ZiOPC) with Correlated Errors
-----------------------------------------------------------

The package also includes :func:`zmiopc.iopcmod` which fits "zero-inflated" ordered probit models (ZiOPC) under the assumption that the two errors are correlated with each other (i.e. correlated errors).

We first import the required libraries, set up the package and import the dataset:

.. testcode::

  # Import the necessary libraries and IDCeMPy.
  import pandas as pd
  import urllib
  import matplotlib.pyplot as plot
  from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset.
  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'

  # Define a `Pandas` DataFrame.
  data = pd.read_stata(url)

.. testcode::

  # First, define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `DataFrame`) used in ordered probit stage.
  # Z = Column names of covariates (from `DataFrame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `DataFrame`).

  X = ['age', 'grade', 'gender_dum']
  Z = ['gender_dum']
  Y = ['cig_count']

:func:`zmiopc.iopcmod` estimates the ZiOPC model using the keyword `'ziopc'` in the first argument:

.. testcode::

   ziopc_tob = zmiopc.iopcmod('ziopc', data, X, Y, Z, method = 'bfgs', weights = 1, offsetx = 0, offsetz = 0)

   # 'ziopc' = model to be estimated. In this case 'ziopc'
   # data = name of Pandas DataFrame
   # X = variables in the ordered probit stage.
   # Y = dependent variable.
   # Z = variables in the inflation stage.
   # method = method for optimization.  By default set to 'bfgs'
   # weights = weights.
   # offsetx = offset of X.  By Default is zero.
   # offsetz = offset of z

The run-time for this ZiOPC model is 4261.707 seconds. The results are stored in the attributes of :class:`zmiopc.IopCModel`.

.. testoutput::

         Current function value: 5060.051910
         Iterations: 119
         Function evaluations: 1562
         Gradient evaluations: 142

The following line of code prints the results:

.. testcode::

    print(ziopc_tob.coefs)

.. testoutput::

                            Coef        SE     tscore             p       2.5%      97.5%
   cut1                   1.696160  0.044726  37.923584  0.000000e+00   1.608497   1.783822
   cut2                  -0.758095  0.033462 -22.655678  0.000000e+00  -0.823679  -0.692510
   cut3                  -1.812077  0.060133 -30.134441  0.000000e+00  -1.929938  -1.694217
   cut4                  -0.705836  0.041432 -17.036110  0.000000e+00  -0.787043  -0.624630
   Inflation: int         9.538072  3.470689   2.748178  5.992748e-03   2.735521  16.340623
   Inflation: gender_dum -9.165963  3.420056  -2.680062  7.360844e-03 -15.869273  -2.462654
   Ordered: age          -0.028606  0.008883  -3.220369  1.280255e-03  -0.046016  -0.011196
   Ordered: grade         0.177541  0.010165  17.465452  0.000000e+00   0.157617   0.197465
   Ordered: gender_dum    0.602136  0.053084  11.343020  0.000000e+00   0.498091   0.706182
   rho                   -0.415770  0.074105  -5.610526  2.017123e-08  -0.561017  -0.270524

To print the estimates of the log-likelihood, AIC, and Variance-Covariance matrix:

.. testcode::

  # Print Log-Likelihood
  print(ziopc_tob.llik)

  # Print AIC
  print(ziopc_tob.AIC)

  # Print VCOV matrix
  print(ziopc_tob.vcov)

The AIC of the ziopc_tob model, for example, is:

.. testoutput::

  10140.103819465658

The predicted probabilities from the `ziopc_tob` model can be obtained with :func:`zmiopc.iopcfit` as follows.

.. testcode::

  # Define the model for which you want to estimate the predicted probabilities
  fittedziopc = zmiopc.iopcfit(ziopc_tob)

  # Print predicted probabilities
  print(fittedziopc.responsefull)

.. testoutput::

 array[[0.88223509 0.06878162 0.01445941 0.0241296  0.01039428]
 [0.84550989 0.08074461 0.01940226 0.03589458 0.01844865]
 [0.93110954 0.04346074 0.00825639 0.01264189 0.00453143]
 ...
 [0.73401588 0.12891071 0.03267436 0.06438928 0.04000977]
 [0.87523652 0.06888286 0.01564958 0.0275354  0.01269564]
 [0.82678185 0.0875059  0.02171135 0.04135142 0.02264948]]

Similar to the ZiOP model, :func:`zmiopc.split_effects` and :func:`zmiopc.ordered_effects` can also compute changes in predicted probabilities for the ZiOPC model.

.. testcode::

  ziopcgender = zmiopc.split_effects(ziopc_tob, 1, nsims = 10000)

.. testcode::

  # Calculate change in predicted probabilities
  gender = zmiopc.ordered_effects(ziopc_tob, 1, nsims = 10000)

  # Box-plot of precicted probabilities
  gender.plot.box(grid='False')

Middle-inflated Ordered Probit (MiOP) without Correlated Errors
---------------------------------------------------------------

A Middle-inflated Ordered Probit (MiOP) model should be estimated when the ordered outcome variable is inflated in the middle category.

The following example uses 2004 presidential vote data from Elgun and Tilam (`2007 <https://journals.sagepub.com/doi/10.1177/1065912907305684>`_).

We begin by loading the required libraries and IDCeMPy:

.. testcode::

  # Import the necessary libraries and IDCeMPy.
  import pandas as pd
  import urllib
  import matplotlib.pyplot as plot
  from idcempy import zmiopc

Next, we load the dataset:

.. testcode::

  # Import and read the dataset
  url = 'https://github.com/hknd23/idcempy/raw/main/data/EUKnowledge.dta'

  # Define a `Pandas` DataFrame
  data = pd_read.stata(url)

We then define the lists with the names of the variables used in the model

.. testcode::

  # First, define a list of variable names of X, Z, and Y.
  # X = Column names of covariates (from `DataFrame`) used in ordered probit stage.
  # Z = Column names of covariates (from `DataFrame`) used in split-population stage.
  # Y = Column name of ordinal outcome variable (from `DataFrame`).

  X = ['Xenophobia', 'discuss_politics']
  Z = ['discuss_politics', 'EU_Know_ob']
  Y = ['EU_support_ET']

After importing the dataset and specifying the list of variables from it, the MiOP model is estimated with the following step:

.. testcode::

 # Model estimation:
 miop_EU = zmiopc.iopmod('miop', data, X, Y, Z, method = 'bfgs', weights = 1,offsetx = 0, offsetz = 0)

 # 'miop' = Type of model to be estimated. In this case 'miop'
 # data = name of Pandas DataFrame
 # X = variables in the ordered probit stage.
 # Y = dependent variable.
 # Z = variables in the inflation stage.
 # method = method for optimization.  By default set to 'bfgs'
 # weights = weights.
 # offsetx = offset of X.  By Default is zero.
 # offsetz = offset of z

The following message will appear when the model finishes converging:

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 10857.695490
         Iterations: 37
         Function evaluations: 488
         Gradient evaluations: 61

The run-time for the model is: 18.886 seconds (N= 11887). Print the results of the model with:

.. testcode::

   print(miop_EU.coefs)

.. testoutput::

                                 Coef        SE       tscore         p         2.5%     97.5%
   cut1                        -1.159621  0.049373 -23.487133  0.000000e+00 -1.256392 -1.062851
   cut2                        -0.352743  0.093084  -3.789492  1.509555e-04 -0.535188 -0.170297
   Inflation: int              -0.236710  0.079449  -2.979386  2.888270e-03 -0.392431 -0.080989
   Inflation: discuss_politics  0.190595  0.035918   5.306454  1.117784e-07  0.120197  0.260993
   Inflation: EU_Know_obj       0.199574  0.020308   9.827158  0.000000e+00  0.159770  0.239379
   Ordered: Xenophobia         -0.663551  0.044657 -14.858898  0.000000e+00 -0.751079 -0.576024
   Ordered: discuss_politics    0.023784  0.029365   0.809964  4.179609e-01 -0.033770  0.081339

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object :class:`zmiopc.IopModel` also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).

.. testcode::

   # Print estimates of LL, AIC and VCOV

   # Print Log-Likelihood
   print(miop_EU.llik)

   # Print AIC
   print(miop_EU.AIC)

   # Print VCOV
   print(miop_EU.vcov)


:func:`zmiopc.iopfit` calculates the predicted probabilities for the MiOP model:

.. testcode::

   # Define the model for which you want to estimate the predicted probabilities
   fittedmiop = zmiopc.iopfit(miop_EU)

   # Print predicted probabilities
   print(fittedmiop.responsefull)

The MiOP model can also work with :func:`zmiopc.split_effects` and :func:`zmiopc.ordered_effects` to compute changes in predicted probabilities when the value of a variable changes:

.. testcode::

    # Define model from which predicted probabilities will be estimated and the number of simulations.
    miopxeno = zmiopc.split_effects(miop_EU, 1, nsims = 10000)

To plot the predicted probabilities:

.. testcode::

     # Get box plot of predicted probabilities
     miopxeno.plot.box(grid='False')

.. testcode::

    # Define model from which predicted probabilities will be estimated and the number of simulations.
    xeno = zmiopc.ordered_effects(miop_EU, 2, nsims = 10000)

    # Get box plot of predicted probabilities
    xeno.plot.box(grid='False')

Middle-inflated Ordered Probit (MiOPC) Model with Correlated Errors
-------------------------------------------------------------------

The steps to estimate the Middle-inflated Ordered Probit (MiOPC) with correlated errors is as follows:

First is importing the data and libraries:

.. testcode::

  # Import the necessary libraries and IDCeMPy.
  import pandas as pd
  import urllib
  import matplotlib.pyplot as plot
  from idcempy import zmiopc

Next, we load the dataset:

.. testcode::

  # Import and read the dataset
  url = 'https://github.com/hknd23/idcempy/raw/main/data/EUKnowledge.dta'

  # Define a `Pandas` DataFrame
  data = pd_read.stata(url)

We then define the lists with the names of the variables used in the model:

.. testcode::

   # First, define a list of variable names of X, Z, and Y.
   # X = Column names of covariates (from `DataFrame`) used in ordered probit stage.
   # Z = Column names of covariates (from `DataFrame`) used in split-population stage.
   # Y = Column name of ordinal outcome variable (from `DataFrame`).

   X = ['Xenophobia', 'discuss_politics']
   Z = ['discuss_politics', EU_Know_ob]
   Y = ['EU_support_ET']

The model can be estimated as follows:

.. testcode::

   # Model estimation
   miopc_EU = zmiopc.iopcmod('miopc', data, X, Y, Z, method = 'bfgs', weights = 1,offsetx = 0, offsetz =0 )

   # 'miopc' = Type of model to be estimated. In this case 'miopc'
   # data = name of Pandas DataFrame
   # X = variables in the ordered probit stage.
   # Y = dependent variable.
   # Z = variables in the inflation stage.
   # method = method for optimization.  By default set to 'BFGS'
   # weights = weights.
   # offsetx = offset of X.  By Default is zero.
   # offsetz = offset of z

The run-time for the model is: 1929.000 seconds (N= 11887). Print model coefficients:

.. testcode::

   print(miopc_EU.coefs)

.. testoutput::

                                 Coef  SE     tscore  p     2.5%  97.5%
   cut1                        -1.370 0.044 -30.948 0.000 -1.456 -1.283
   cut2                        -0.322 0.103  -3.123 0.002 -0.524 -0.120
   Inflation: int              -0.129 0.021  -6.188 0.000 -0.170 -0.088
   Inflation: discuss_politics  0.192 0.026   7.459 0.000  0.142  0.243
   Inflation: EU_Know_obj       0.194 0.027   7.154 0.000  0.141  0.248
   Ordered: Xenophobia         -0.591 0.045 -13.136 0.000 -0.679 -0.502
   Ordered: discuss_politics   -0.029 0.021  -1.398 0.162 -0.070  0.012
   rho                         -0.707 0.106  -6.694 0.000 -0.914 -0.500

In addition to coefficient estimates, the table also presents the standard errors, and confidence intervals.

The model object :class:`zmiopc.IopCModel` also stores three different diagnostic tests: (1) Log-likelihood, (2) Akaike Information Criteria (AIC), and Variance-Covariance Matrix (VCM).
They can be obtained via the following:

.. testcode::

   # Print Log-Likelihood
   print(miopc_EU.llik)

   # Print AIC
   print(miopc_EU.AIC)

   # Print VCCOV matrix
   rint(miopc_EU.vcov)

To calculate the predicted probabilities:

.. testcode::

   # Define model to fit
   fittedmiopc = zmiopc.iopcfit(miopc_EU)

   # Print predicted probabilities
   print(fittedziopc.responsefull)

The following line of code computes changes in predicted probabilities when the value of a chosen variable in the split stage changes:

.. testcode::

   # Define model from which effects will be estimated and number of simulations
   miopcxeno = zmiopc.split_effects(miopc_EU, 1, nsims = 10000)

A box plot can illustrate the change in predicted probabilities:

.. testcode::

    # Get box plot of predicted probabilities
    miopcxeno.plot.box(grid='False')


To calculate the change in predicted probabilities of the outcome variable in the outcome-stage when the value of a covarariate changes.
The box plots below display the change in predicted probabilities of the outcome variable in the MiOPC model estimated above when Xenophobia increases one standard deviation from its mean value.

.. testcode::

    # Define model from which effects will be estimated and number of simulations
    xeno = zmiopc.ordered_effects(miopc_EU, 2, nsims = 10000)

    # Get box plot of predicted probabilities
    xeno.plot.box(grid='False')


The Standard Ordered Probit (OP) model
--------------------------------------

The package also includes :func:`zmiopc.opmod` that estimates a standard Ordered Probit (OP) model.
The OP model does not account for "zero inflation" or "middle inflation," so it does not have a split-probit stage.

First, import the required libraries and data:

.. testcode::

  # Import the necessary libraries and package
  import pandas as pd
  import urllib
  from idcempy import zmiopc

  # Import the "Youth Tobacco Consumption" dataset.
  url='https://github.com/hknd23/zmiopc/blob/main/data/tobacco_cons.csv'

  # Define a `Pandas` DataFrame
  data = pd.read_csv(url)

The list of variable names for the Independent and Dependent variables needs to be specified:

.. testcode::

  # Define a list of variable names (strings) X,Y:
  # X = Column names of covariates (from `DataFrame`) in the OP equation
  # Y = Column name of outcome variable (from `DataFrame`).

  X = ['age', 'grade', 'gender_dum']
  Y = ['cig_count']

After importing the data and specifying the model, the following code fits the OP model:

.. testcode::

  # Model estimation:
  op_tob = zmiopc.opmod(data, X, Y, method = 'bfgs', weights = 1, offsetx  =0)

  # data = name of pandas DataFrame
  # X = variables in the ordered probit stage.
  # Y = dependent variable.
  # method = method for optimization.  By default set to 'bfgs'
  # weights = weights.
  # offsetx = offset of X.  By Default is zero.
  # offsetz = offset of z


The following message will appear when the model has converged:

.. testoutput::

         Warning: Desired error not necessarily achieved due to precision loss.
         Current function value: 4411.710049
         Iterations: 10
         Function evaluations: 976
         Gradient evaluations: 121

The model's run-time is 37.694 seconds (N= 9624). :class:`zmiopc.OpModel` stores results from model estimation and other information in its attributes.
The following line of code to see the estimates of coefficients:

.. testcode::

   # Print coefficients of the models
   print(op_tob.coefs)

.. testoutput::

                Coef        SE     tscore         p      2.5%     97.5%
   cut1        1.696175  0.047320  35.844532  0.000000  1.603427  1.788922
   cut2       -0.705037  0.031650 -22.276182  0.000000 -0.767071 -0.643004
   cut3       -2.304405  0.121410 -18.980329  0.000000 -2.542369 -2.066441
   cut4        2.197381  0.235338   9.337141  0.000000  1.736119  2.658643
   age        -0.070615  0.007581  -9.314701  0.000000 -0.085474 -0.055756
   grade       0.233741  0.010336  22.614440  0.000000  0.213483  0.254000
   gender_dum  0.020245  0.032263   0.627501  0.530331 -0.042991  0.083482

Log-likelihood, AIC, and Variance-Covariance matrix can be extracted with:

.. testcode::

  # Print Log-Likelihood
  print(op_tob.llik)

  # Print AIC
  print(op_tob.AIC)

  # Print VCOV matrix
  print(op_tob.vcov)

The Vuong Test
--------------

Harris and Zhao (`2007 <https://doi.org/10.1016/j.jeconom.2007.01.002>`__) suggest that a variant of the Vuong (`1989 <https://www.jstor.org/stable/1912557>`__)
Test (with a v statistic) can be used to compare the performance of the ZiOP versus the standard Ordered Probit (OP) model. The Vuong's test formula is:

.. math::

    v = \frac{\sqrt{N}(\frac{1}{N}\sum_{i}^{N}m_{i})}{\sqrt{\frac{1}{N}\sum_{i}^{N}(m_{i}-\bar{m})^{2}}}

where v < -1.96 favors the more general (ZiOP/ZiOPC) model, -1.96 < v < 1.96 lends no support to either model, and v > 1.96 supports the simpler (OP) model.

The OP and ZiOP models must have the same number of observations, and the OP must have the same number of covariates as ZiOP's OP stage.
The statistic below reveals that the OP model is preferred over the ZiOP model.

.. testcode::

   # Estimate Vuong test.  OP model first, ZIOP model specified next in this case
   zmiopc.vuong_opiop(op_tob, ziop_tob)

.. testoutput::

   6.624742132792222

The Vuong test can also be implemented to compare the ZiOPC, MiOP and MiOPC models with the OP model.

Generalized Inflated Multinomial Logit (GiMNL) Model
----------------------------------------------------

The :py:mod:`gimnl` module provides :func:`gimnl.gimnlmod` to estimate the General "inflated" Multinomial Logit models (GiMNL) with three outcomes in the dependent variable.
The GiMNL model minimize issues present when unordered polytomous outcome variables have an excessive share and heterogeneous pool of observations in the lower category.

Similar to the models in the :py:mod:`zmiopc` module, the first step is to import the libraries and 2004 presidential vote choice dataset.

.. testcode::

  # Import the module
  import pandas as pd
  import urllib
  from idcempy import gimnl

  # Load the dataset
  url= 'https://github.com/hknd23/zmiopc/raw/main/data/replicationdata.dta'

  # Define a `Pandas` DataFrame
  data = pd.read_stata(url)

We the define the list of covariates in the split-stage (z), the multinomial logit-stage (x) and the outcome variable (y).
The values of the dependent variable must be represented numerically as "0", "1", and "2" to represent each category.
To specify the baseline/reference category, users provide a three-element list for the `reference` argument (e.g [0,1,2]).
The first element of the list is the baseline/reference category.

.. testcode::

   # x = Column names of covariates (from `DataFrame`) in the outcome-stage.
   # z = Column names of covariates (from `DataFrame`) in the split-stage.
   # y = Column names of outcome variable (from `DataFrame`).

   x = ['educ', 'party7', 'agegroup2']
   z = ['educ', 'agegroup2']
   y = ['vote_turn']

The flexibility of :func:`gimnl.gimnlmod` allows users to customize the baseline and inflated categories.
Users can employ the argument `inflatecat` with `'baseline'`, `'second'`, or `'third'` to specify any unordered category as the inflated category (dictated by the distribution) in their unordered-polytomous outcome measure.
If `'baseline'` is selected, the first element (baseline/reference category) in `reference` is the inflated outcome.
Likewise, if `'second'` or `'third'` is selection, the second or third element will be the inflated outcome. The following code specifies the outcome '0' (Abstain) as both the baseline and inflated category.

.. testcode::

   # Define order of variables
   order = [0, 1, 2]

   # Define "inflation" category
   inflatecat = "baseline"

.. testcode::

   # Estimate the model
   gimnl_2004vote = gimnl.gimnlmod(data, x, y, z, method = 'bfgs', order, inflatecat)

   # data = name of pandas DataFrame.
   # x = variables in the ordered stage.
   # y = dependent variable.
   # z = variables in the inflation stage.
   # method = optimization method.  Default is 'bfgs'
   # order = order of variables.
   # inflatecat = inflated category.

The following line of code prints the coefficients of the covariates:

.. testcode::

   # Print coefficients
   print(gimnl_2004vote.coefs)

.. testoutput::

                          Coef   SE    tscore   p    2.5%   97.5%
   Inflation: int       -4.935 2.777  -1.777 0.076 -10.379  0.508
   Inflation: educ       1.886 0.293   6.441 0.000   1.312  2.460
   Inflation: agegroup2  1.295 0.768   1.685 0.092  -0.211  2.800
   1: int               -4.180 1.636  -2.556 0.011  -7.387 -0.974
   1: educ               0.334 0.185   1.803 0.071  -0.029  0.697
   1: party7             0.454 0.057   7.994 0.000   0.343  0.566
   1: agegroup2          0.954 0.248   3.842 0.000   0.467  1.441
   2: int                0.900 1.564   0.576 0.565  -2.166  3.966
   2: educ               0.157 0.203   0.772 0.440  -0.241  0.554
   2: party7            -0.577 0.058  -9.928 0.000  -0.691 -0.463
   2: agegroup2          0.916 0.235   3.905 0.000   0.456  1.376

The model's run-time is 16.646 seconds (N= 1341). The results from the model are stored in a :class:`gimnlModel` with the following attributes:

- coefs: Model coefficients and standard errors.
- llik: Log-likelihood.
- AIC: Akaike information criterion.
- vcov: Variance-covariance matrix.

For example, AIC can be printed as follows.

.. testcode::

  # Print Log_Likelihood
  print(gimnl_2004vote.llik)

  # Print AIC
  print(gimnl_2004vote.AIC)

  # Print VCOV matrix
  print(gimnl_2004vote.vcov)

Users can fit a standard three-category Multinomial Logit Model (MNL) by specifying the list of **x**, **y**, and baseline (using `reference`).

.. testcode::

   #Estimate the model
   mnl_2004vote = gimnl.mnlmod(data, x, y, method = 'bfgs')

   # data = name of Pandas DataFrame.
   # x = variables in MNL stage.
   # y = dependent variable
   # method = optimization method. Default is 'bfgs'

   # Print the coefficients
   print(mnl_2004vote.coefs)

.. testoutput::

     Coef        SE  tscore     p   2.5%  97.5%
  1: int       -4.914 0.164 -29.980 0.000 -5.235 -4.593
  1: educ       0.455 0.043  10.542 0.000  0.371  0.540
  1: party7     0.462 0.083   5.571 0.000  0.300  0.625
  1: agegroup2  0.951 0.029  32.769 0.000  0.894  1.008
  2: int        0.172 0.082   2.092 0.036  0.011  0.334
  2: educ       0.282 0.031   9.011 0.000  0.221  0.343
  2: party7    -0.567 0.085  -6.641 0.000 -0.734 -0.399
  2: agegroup2  0.899 0.138   6.514 0.000  0.629  1.170

The MNL model's run-time is 8.276 seconds. Similar to the GiMNL model, the AIC for the MNL model can also be given by:

.. testcode::

  # Print Log-Likelihood
  print(mnl_2004vote.AIC)

  # Print AIC
  print(mnl_2004vote.AIC)

  # Print VCOV matrix
  print(mnl_2004vote.vcov)

Contributions
=============

The authors welcome and encourage new contributors to help test `IDCeMPy` and add new functionality.
You can find detailed instructions on "how to contribute" to `IDCeMPy` `here <https://github.com/hknd23/idcempy/blob/main/CONTRIBUTING.md>`_.
.. toctree::
   :maxdepth: 3

API Reference
=============

Classes and functions in idcempy package. For tutorial on how to use the package to estimate the ZiOP(C) and MiOP(C) models, see :doc:`zmiopc_tutorial`.

For tutorial on the GiMNL module, see and :doc:`gimnl_tutorial`.


The zmiopc module
------------------

.. automodule:: zmiopc
   :members:
   :member-order: bysource

The gimnl module
------------------

.. automodule:: gimnl
  :members:
  :member-order: bysource
