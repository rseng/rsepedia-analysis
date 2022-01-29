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
