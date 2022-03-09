
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtplay\_1.6.3 (2021-11-01)

o resolved the issue occurred when fixing the item guessing parameters
to a specific (e.g., 0.1) in the `est_irt()` function.

o updated the `est_irt()` function to estimate the population latent
ability distribution only when all item parameters are fixed in a test
using the fixed item parameter calibration (FIPC).

o fixed the `est_irt()` function so that the log-likelihood, AIC, and
BIC can be computed based on both the fixed- and freely estimated items
when the FIPC is implemented. In the previous version, those valused
were computed based on only freely estimated items.

o added a new argument of ‘item.id’ in the ‘est\_irt()’ and
‘est\_item()’ functions where a user can provide item IDs.

o added a new ‘rdif()’ function which computes RDIF statistics (Lim,
Choe, & Han, under review; Lim, Choe, Han, Lee, & Hong, 2021) for
analyzing DIF.

o updated ‘plot.test.info()’ function so that (a) multiple item
information functions can be displayed in one panel by setting ‘overlap
= TRUE’ and (b) a plot of conditional standard error of estimation at a
test level can be shown by setting ‘csee = TRUE’.

o updated ‘est\_score()’ function to make it return NA values for
examinees who have all missing responses.

o fixed an error of ‘est\_score()’ function which occurs when an
examinee has missing data for all polytomous items or dichotomous items
(thanks to Craig Wells).

o fixed a few minor issues of ‘irtfit()’ function (thanks to Dimitrios
Zacharatos).

o Updated ‘bring.flexmirt()’ function to read the empirical histogram of
population distribution from “-prm.txt” file.

o Updated ‘run\_flexmirt()’ function to run flexMIRT in which version is
\>= 3.6.

o Updated ‘est\_item()’ function to produce a variance-covariance matrix
for item parameter estimates.

o Added ‘vcov()’ method for ‘est\_item’ function.

# irtplay\_1.6.2 (2020-12-14)

o Added ‘coef’, ‘logLik’ methods for ‘est\_item’ function, and added
‘coef’, ‘logLik’, and ‘vcov’ methods for ‘est\_irt’ function.

o Updated ‘est\_irt’ and ‘est\_item’ functions so that the argument of
‘verbose’ can suppress all messages of the parameter estimation
progress.

o Updated ‘bring.flexmirt’ function by including a new argument of
“rePrm.gpc”. In the previous version, the nominal model parameters in
the flexMIRT parameter output file are reparameterized into the (G)PCM
parameters when (G)PCM is fit to data by setting ‘rePrm = TRUE’. In the
new version, however, the nominal model parameters are reparameterized
into the (G)PCM slope/difficulty parameters only when ‘rePrm.gpc =
TRUE’.

# irtplay\_1.6.1 (2020-08-13)

o Included a new function of ‘post\_den’ to compute updated prior
(a.k.a. posterior) densities of the latent ability distribution given a
prior ability distribution, item parameters, and item response data.

o Updated ‘est\_score’ function so that the standard errors of ability
estimates can be computed using eigther the observed item information
function or the expected item information (a.k.a. Fisher information)
function when MLE, MLE with fence (MLEF), or MAP scoring method is used.

o Updated ‘est\_irt’ function to compute the loglikelihood-based fit
statistics of Akaike information criterion (AIC) and Bayesian
information criterion (BIC).

o Updated ‘est\_irt’ function to tally the number of freely estimated
parameters taking the mean and variance parameters of the latent ability
distribution into consideration when ‘fipc = TRUE’.

o Updated ‘est\_irt’ function to suppress printing the observed data
log-likelihood after each EM cycle using the argument of ‘verbose’.

o Fixed an error of the ‘est\_irt’ function when only dichotomous items
are used with ‘fipc = TRUE’. In that condition, an error message of
“subscript out of bounds” was returned in the previous version. No
error message is shown in the updated version. (thanks to Ahmet GUVEN)

o Fixed the ‘lwrc’ function so that it can return the probability
results even when only a single theta value is used.

# irtplay\_1.6.0 (2020-07-14)

The package has been updated significantly in this verstion. In this
version, I have:

o Updated ‘est\_score’ function to estimate ability parameters much
faster than the previous version of the function.

o Updated ‘est\_irt’ and ‘est\_item’ functions to estimate item
parameters much faster than the previous version of the functions.

o Updated ‘test.info’ function to compute items infomation and test
information much faster than the previous version of the function.

o Added an option to use a prior distribution of the item difficulty (or
threshold) parameters in ‘est\_irt’, ‘est\_item’, and ‘llike\_item’
functions.

o Solved unstable item parameter estimation of ‘est\_irt’ and
‘est\_item’ functions which occured when the scaling factor of ‘D’
is other than 1.0 and ‘use.aprior = TRUE’.

o Fixed an error which occured in the function ‘est\_irt’ when the data
set contains missing values and ‘fix.a.1pl = FALSE’.

# irtplay\_1.5.1 (2020-06-16)

o Included ‘summary’ method to summarize the IRT calibration results
from ‘est\_irt’ or ‘est\_item’ objects.

o Included a new function of ‘getirt’ to extract various estimates
results from ‘est\_irt’ or ‘est\_item’ objects.

o Fixed an error which happens when “DRM” is specified in the model name
in the function ‘est\_irt’.

o Included total computation time in the function ‘est\_irt’.

# irtplay\_1.5.0 (2020-04-12)

o Changed the title of ‘irtplay’ package to “Unidimensional Item
Response Theory Modeling”.

o Included a new function of ‘est\_irt’ to fit unidimensional IRT models
to mixture of dichotomous and polytomous item data using the marginal
maximum likelihood estimation with expectation-maximization (MMLE-EM;
Bock & Aitkin, 1981) algorithm.

o Included the fixed item parameter calibration (FIPC; Kim, 2006)
approach, which is one of useful online calibration methods, in the
function ‘est\_irt’.

o Updated the documentation to explain how to implement the new function
‘est\_irt’.

o Included well-known LSAT6 dichotomous response data set from Thissen
(1982).

o Fixed a problem of inaccurate item parameter estimation in the
function ‘est\_item’ when a prior distribution of the slope parameter is
used with a scaling factor other than D = 1.

o Updated the function ‘bring.flexmirt’ to read the item parameters of
the generalized partial credit model when the number of score categories
are two.

o Updated the function ‘est\_score’ to find a smart starting value when
MLE is used. More specifically, the smart starting value is a theta
value where the log-likelihood is the maximum at the highest peak.

# irtplay\_1.4.1 (2020-02-21)

o Included the function ‘run\_flexmirt’ to implement flexMIRT software
(Cai, 2017) through R.

o Applied a prior distribution to the slope parameters of the IRT 1PL
model when the slope parameters are constrained to be equal in the
function of ‘est\_item’.

o Fixed a problem of using staring values to estimate item parameters in
the function of ‘est\_item’.

# irtplay\_1.4.0 (2020-01-23)

o Fixed a non-convergence problem of the maximum likelihood estimation
with fences (MLEF) in the function of ‘est\_score’.

o Updated the description and introduction of the package.

o Updated the documentation to explain how to implement the function
“est\_item” in more detail.

o Updated the README.md file to explain how to implement the function
“est\_item” in more detail.

# irtplay\_1.3.0 (2019-11-17)

o Included the function ‘llike\_score’ to compute the loglikelihood
function of ability for an examinee.

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_item’ to exclude items that contains no
item response data during the item parameter estimation.

o Updated the function ‘est\_item’ to count the number of item responses
for each item used to estimate the item parameters.

o Updated the function ‘est\_score’ to find better starting values when
MLE is used.

o Updated the function ‘est\_score’ to address NaNs of gradient values
and NaNs of hessian values when MLE, MLEF, or MAP is used.

o Fixed a problem of the function ‘est\_score’, which returned an error
message when a vector of an examinee’s response data was used in the
argument of ‘x’.

# irtplay\_1.2.0 (2019-10-16)

o Fixed a problem of the function ‘est\_score’, which returned an error
message when only one dichotomous item or one polytomous item was
included in the item meta data set.

o Fixed a problem of the function ‘est\_item’, which returned an error
message when the inverse of hessian matrix is not obtainable.

o Included the ‘maximum likelihood estimation with fences scoring method
(Han, 2016) in the function ’est\_score’.

o Included the ‘inverse test characteristic curve (TCC)’ scoring method
(e.g., Stocking, 1996) in the function ‘est\_score’.

o Included the function ‘llike\_item’ to compute the loglikelihood
values of items.

# irtplay\_1.1.0 (2019-09-15)

o For the function ‘est\_item’, default parameters of a-parameter prior
distribution were revised

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_score’ to estimate an ability in a brute
force way when MLE or MAP fails to find the solution.

o Updated the function ‘irtfit’ to compute the likelihood ratio
chi-square fit statistic (G2; Mckinley & Mills, 1985).

# irtplay\_1.0.0 (2019-08-21)

o initial release on CRAN
irtplay
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

The goal of `irtplay` is to fit unidimensional item response theory
(IRT) models to mixture of dichotomous and polytomous data, calibrate
online item parameters (i.e., pretest and operational items), estimate
examinees abilities, and examine the IRT model-data fit on item-level in
different ways as well as provide useful functions related to
unidimensional IRT.

For the item parameter estimation, the marginal maximum likelihood
estimation with expectation-maximization (MMLE-EM) algorithm (Bock &
Aitkin, 1981) is used. For the online calibration, the fixed item
parameter calibration (FIPC) method (e.g., Kim, 2006) and the fixed
ability parameter calibration (FAPC) method (Ban, Hanson, Wang, Yi, &
Harris, 2001; stocking, 1988), often called Stocking’s Method A, are
provided. For the ability estimation, several popular scoring methods
(e.g., MLE, EAP, and MAP) are implemented. In terms of assessing the IRT
model-data fit, one of distinguished features of this package is that it
gives not only item fit statistics (e.g., chi-square fit statistic (X2;
e.g., Bock, 1960; Yen, 1981), likelihood ratio chi-square fit statistic
(G2; McKinley & Mills, 1985), infit and outfit statistics (Ames et al.,
2015), and S-X2 (Orlando & Thissen, 2000, 2003)) but also graphical
displays to look at residuals between the observed data and model-based
predictions (Hambleton, Swaminathan, & Rogers, 1991).

In addition, there are many useful functions such as analyzing
differential item functioning (DIF), computing asymptotic
variance-covariance matrices of item parameter estimates, importing item
and/or ability parameters from popular IRT software, running flexMIRT
(Cai, 2017) through R, generating simulated data, computing the
conditional distribution of observed scores using the Lord-Wingersky
recursion formula, computing item and test information functions,
computing item and test characteristic curve functions, and plotting
item and test characteristic curves and item and test information
functions.

## Installation

You can install the released version of irtplay from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("irtplay")
```

Or you can install the latest development version from Github using the
*devtools* package:

``` r
install.packages("devtools")
devtools::install_github("hwangQ/irtplay")
```

## 1\. Online item calibration with the fixed item parameter calibration (FIPC) method (e.g., Kim, 2006)

The fixed item parameter calibration (FIPC) is one of useful online item
calibration methods for computerized adaptive testing (CAT) to put the
parameter estimates of pretest items on the same scale of operational
item parameter estimates without post hoc linking/scaling (Ban, Hanson,
Wang, Yi, & Harris, 2001; Chen & Wang, 2016). In FIPC, the operational
item parameters are fixed to estimate the characteristic of the
underlying latent variable prior distribution when calibrating the
pretest items. More specifically, the underlying latent variable prior
distribution of the operational items is estimated during the
calibration of the pretest items to put the item parameters of the
pretest items on the scale of the operational item parameters (Kim,
2006). In `irtplay` package, FIPC is implemented with two main steps:

1.  Prepare a response data set and the item metadata of the fixed (or
    operational) items.
2.  Implement FIPC to estimate the item parameters of pretest items
    using the `est_irt()` function.

### (1) Preparing a data set

To run the `est_irt()` function, it requires two data sets:

1.  Item metadata set (i.e., model, score category, and item parameters.
    see the desciption of the argument `x` in the function `est_irt`).
2.  Examinees’ response data set for the items. It should be a matrix
    format where a row and column indicate the examinees and the items,
    respectively. The order of the columns in the response data set must
    be exactly the same as the order of rows of the item metadata.

### (2) Estimating the pretest item parameters

When FIPC is implemented in `est_irt()` function, the pretest item
parameters are estimated by fixing the operational item parameters. To
estimate the item parameters, you need to provide the item metadata in
the argument `x` and the response data in the argument `data`.

It is worthwhile to explain about how to prepare the item metadata set
in the argument `x`. A specific form of a data frame should be used for
the argument `x`. The first column should have item IDs, the second
column should contain the number of score categories of the items, and
the third column should include IRT models. The available IRT models are
“1PLM”, “2PLM”, “3PLM”, and “DRM” for dichotomous items, and “GRM” and
“GPCM” for polytomous items. Note that “DRM” covers all dichotomous
IRT models (i.e, “1PLM”, “2PLM”, and “3PLM”) and “GRM” and “GPCM”
represent the graded response model and (generalized) partial credit
model, respectively. From the fourth column, item parameters should be
included. For dichotomous items, the fourth, fifth, and sixth columns
represent the item discrimination (or slope), item difficulty, and item
guessing parameters, respectively. When “1PLM” or “2PLM” is specified
for any items in the third column, NAs should be inserted for the item
guessing parameters. For polytomous items, the item discrimination (or
slope) parameters should be contained in the fourth column and the item
threshold (or step) parameters should be included from the fifth to the
last columns. When the number of categories differs between items, the
empty cells of item parameters should be filled with NAs. See `est_irt`
for more details about the item metadata.

Also, you should specify in the argument `fipc = TRUE` and a specific
FIPC method in the argument `fipc.method`. Finally, you should provide a
vector of the location of the items to be fixed in the argument
`fix.loc`. For more details about implementing FIPC, see the description
of the `est_irt()` function.

When implementing FIPC, you can estimate both the emprical histogram and
the scale of latent variable prior distribution by setting `EmpHist =
TRUE`. If `EmpHist = FALSE`, the normal prior distribution is used
during the item parameter estimation and the scale of the normal prior
distribution is updated during the EM cycle.

If necessary, you need to specify whether prior distributions of item
slope and guessing parameters (only for the IRT 3PL model) are used in
the arguments of `use.aprior` and `use.gprior`, respectively. If you
decide to use the prior distributions, you should specify what
distributions will be used for the prior distributions in the arguments
of `aprior` and `gprior`, respectively. Currently three probability
distributions of Beta, Log-normal, and Normal distributions are
available.

In addition, if the response data include missing values, you must
indicate the missing value in argument `missing`.

Once the `est_irt()` function has been implemented, you’ll get a list of
several internal objects such as the item parameter estimates, standard
error of the parameter estimates.

## 2\. Online item calibration with the fixed ability parameter calibration method (e.g., Stocking, 1988)

In CAT, the fixed ability parameter calibration (FAPC), often called
Stocking’s Method A, is the relatively simplest and most straightforward
online calibration method, which is the maximum likelihood estimation of
the item parameters given the proficiency estimates. In CAT, FAPC can be
used to put the parameter estimates of pretest items on the same scale
of operational item parameter estimates and recalibrate the operational
items to evaluate the parameter drifts of the operational items (Chen &
Wang, 2016; Stocking, 1988). Also, FAPC is known to result in accurate,
unbiased item parameters calibration when items are randomly rather than
adaptively administered to examinees, which occurs most commonly with
pretest items (Ban et al., 2001; Chen & Wang, 2016). Using `irtplay`
package, the FAPC is implemented to calibrate the items with two main
steps:

1.  Prepare a data set for the calibration of item parameters (i.e.,
    item response data and ability estimates).
2.  Implement the FAPC to estimate the item parameters using the
    `est_item()` function.

### (1) Preparing a data set

To run the `est_item()` function, it requires two data sets:

1.  Examinees’ ability (or proficiency) estimates. It should be in the
    format of a numeric vector.
2.  Examinees’ response data set for the items. It should be in the
    format of matrix where a row and column indicate the examinees and
    the items, respectively. The order of the examinees in the response
    data set must be exactly the same as that of the examinees’ ability
    estimates.

### (2) Estimating the pretest item parameters

The `est_item()` function estimates the pretest item parameters given
the proficiency estimates. To estimate the item parameters, you need to
provide the response data in the argument `data` and the ability
estimates in the argument `score`.

Also, you should provide a string vector of the IRT models in the
argument `model` to indicate what IRT model is used to calibrate each
item. Available IRT models are “1PLM”, “2PLM”, “3PLM”, and “DRM” for
dichotomous items, and “GRM” and “GPCM” for polytomous items. “GRM” and
“GPCM” represent the graded response model and (generalized) partial
credit model, respectively. Note that “DRM” is considered as “3PLM” in
this function. If a single character of the IRT model is specified, that
model will be recycled across all items.

The `est_item()` function requires a vector of the number of score
categories for the items in the argument `cats`. For example, a
dichotomous item has two score categories. If a single numeric value is
specified, that value will be recycled across all items. If NULL and all
items are binary items (i.e., dichotomous items), it assumes that all
items have two score categories.

If necessary, you need to specify whether prior distributions of item
slope and guessing parameters (only for the IRT 3PL model) are used in
the arguments of `use.aprior` and `use.gprior`, respectively. If you
decide to use the prior distributions, you should specify what
distributions will be used for the prior distributions in the arguments
of `aprior` and `gprior`, respectively. Currently three probability
distributions of Beta, Log-normal, and Normal distributions are
available.

In addition, if the response data include missing values, you must
indicate the missing value in argument `missing`.

Once the `est_item` function has been implemented, you’ll get a list of
several internal objects such as the item parameter estimates, standard
error of the parameter estimates.

## 3\. The process of evaluating the IRT model-data fit

One way to assess goodness of IRT model-data fit is through an item fit
analysis by examining the traditional item fit statistics and looking at
the discrepancy between the observed data and model-based predictions.
Using `irtplay` package, the traditional approach of evaluating the IRT
model-data fit on item-level can be implemented with three main steps:

1.  Prepare a data set for the IRT item fit analysis (i.e., item
    metadata, ability estimates, and response data).
2.  Obtain the IRT fit statistics such as the X2, G2, infit, and outfit
    statistics using the `irtfit()` function.
3.  Based on the results of IRT model fit analysis (i.e., an object of
    class `irtfit`) obtained in step 2, draw the IRT residual plots
    (i.e., raw residual and standardized residual plots) using `plot`
    method.

### (1) Preparing a data set

Before conducting the IRT model fit analysis, it is necessary to prepare
a data set. To run the `irtfit()` function, it requires three data sets:

1.  Item metadata including the item ID, number of score categories, IRT
    models, and item parameters. The item metadata should be in the
    format of data frame. You can prepare the data either by using the
    `shape_df()` function or by creating a data frame of the item
    metadata by yourself. If you have output files of item parameter
    estimates obtained from one of the IRT software such as BILOG-MG 3,
    PARSCALE 4, flexMIRT, and mirt (R package), the item metadata can be
    easily obtained using the functions of `bring.bilog()`,
    `bring.parscale()`, `bring.flexmirt()`, `bring.mirt()`. See the
    functions of `irtfit()`, `test.info()`, or `simdat()` for more
    details about the item metadata format.
2.  Examinees’ ability (or proficiency) estimates. It should be in the
    format of a numeric vector.
3.  Examinees’ response data set for the items. It should be in the
    format of matrix where a row and column indicate the examinees and
    the items, respectively. The order of the examinees in the response
    data set must be exactly the same as that of the examinees’ ability
    estimates. The order of the items in the response data set must be
    exactly the same as that of the items in the item metadata.

### (2) Computing the IRT model-data fit statistics

The `irtfit()` function computes the traditional IRT item fit statistics
such as X2, G2, infit, and outfit statistics. To calculate the X2 and G2
statistics, two methods are available to divide the ability scale into
several groups. The two methods are “equal.width” for dividing the scale
by an equal length of the interval and “equal.freq” for dividing the
scale by an equal frequency of examinees. Also, you need to specify the
location of ability point at each group (or interval) where the expected
probabilities of score categories are calculated from the IRT models.
Available locations are “average” for computing the expected probability
at the average point of examinees’ ability estimates in each group and
“middle” for computing the expected probability at the midpoint of
each group.

To use the `irtfit()` function, you need to insert the item metadata in
the argument `x`, the ability estimates in the argument `score`, and the
response data in the argument `data`. If you want to divide the ability
scale into other than ten groups, you need to specify the number of
groups in the argument `n.width`. In addition, if the response data
include missing values, you must indicate the missing value in argument
`missing`.

Once the `irtfit()` function has been implemented, you’ll get the fit
statistic results and the contingency tables for every item used to
calculate the X2 and G2 fit statistics.

### (3) Drawing the IRT residual plots

Using the saved object of class `irtfit`, you can use the `plot` method
to evaluate the IRT raw residual and standardized residual plots.

Because the `plot` method can draw the residual plots for an item at a
time, you have to indicate which item will be examined. For this, you
can specify an integer value, which is the location of the studied item,
in the argument `item.loc`.

In terms of the raw residual plot, the argument `ci.method` is used to
select a method to estimate the confidence intervals among four methods.
Those methods are “wald” for the Wald interval, which is based on the
normal approximation (Laplace, 1812), “cp” for Clopper-Pearson interval
(Clopper & Pearson, 1934), “wilson” for Wilson score interval (Wilson,
1927), and “wilson.cr” for Wilson score interval with continuity
correction (Newcombe, 1998).

## 3\. Examples of implementing online calibration and evaluating the IRT model-data fit

``` r
library("irtplay")

##----------------------------------------------------------------------------
# 1. The example code below shows how to prepare
# the data sets and how to implement the fixed
# item parameter calibration (FIPC):
##----------------------------------------------------------------------------

## Step 1: prepare a data set In this example, we
## generated examinees' true proficiency
## parameters and simulated the item response
## data using the function 'simdat'.

## import the '-prm.txt' output file from
## flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt",
  package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df

# generate 1,000 examinees' latent abilities from
# N(0.4, 1.3)
set.seed(20)
score <- rnorm(1000, mean = 0.4, sd = 1.3)

# simulate the response data
sim.dat <- simdat(x = x, theta = score, D = 1)

## Step 2: Estimate the item parameters fit the
## 3PL model to all dichotmous items, fit the GRM
## model to all polytomous data, fix the five 3PL
## items (1st - 5th items) and three GRM items
## (53th to 55th items) also, estimate the
## empirical histogram of latent variable
fix.loc <- c(1:5, 53:55)
(mod.fix1 <- est_irt(x = x, data = sim.dat, D = 1,
  use.gprior = TRUE, gprior = list(dist = "beta",
    params = c(5, 16)), EmpHist = TRUE, Etol = 0.001,
  fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc,
  verbose = FALSE))
#> 
#> Call:
#> est_irt(x = x, data = sim.dat, D = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>     params = c(5, 16)), EmpHist = TRUE, Etol = 0.001, fipc = TRUE, 
#>     fipc.method = "MEM", fix.loc = fix.loc, verbose = FALSE)
#> 
#> Item parameter estimation using MMLE-EM. 
#> 36 E-step cycles were completed using 49 quadrature points.
#> First-order test: Convergence criteria are satisfied.
#> Second-order test: Solution is a possible local maximum.
#> Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Log-likelihood: -24845.46
summary(mod.fix1)
#> 
#> Call:
#> est_irt(x = x, data = sim.dat, D = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>     params = c(5, 16)), EmpHist = TRUE, Etol = 0.001, fipc = TRUE, 
#>     fipc.method = "MEM", fix.loc = fix.loc, verbose = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items: 55
#>  Number of Cases: 1000
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 147
#>  Number of fixed items: 8
#>  Number of E-step cycles completed: 36
#>  Maximum parameter change: 0.0009203804
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 2.93
#>  Standard error computation: 1.16
#>  Total computation: 4.33
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 49690.92
#>  Akaike Information Criterion (AIC): 49984.92
#>  Bayesian Information Criterion (BIC): 50706.36
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   3PLM   0.76    NA   1.46    NA   0.26    NA     NA    NA
#> 2    CMC2     2   3PLM   1.92    NA  -1.05    NA   0.18    NA     NA    NA
#> 3    CMC3     2   3PLM   0.93    NA   0.39    NA   0.10    NA     NA    NA
#> 4    CMC4     2   3PLM   1.05    NA  -0.41    NA   0.20    NA     NA    NA
#> 5    CMC5     2   3PLM   0.87    NA  -0.12    NA   0.16    NA     NA    NA
#> 6    CMC6     2   3PLM   1.47  0.15   0.61  0.09   0.07  0.03     NA    NA
#> 7    CMC7     2   3PLM   1.45  0.25   1.23  0.14   0.24  0.04     NA    NA
#> 8    CMC8     2   3PLM   0.80  0.11   0.82  0.21   0.12  0.05     NA    NA
#> 9    CMC9     2   3PLM   0.81  0.13   0.63  0.28   0.20  0.07     NA    NA
#> 10  CMC10     2   3PLM   1.55  0.21   0.16  0.14   0.18  0.05     NA    NA
#> 11  CMC11     2   3PLM   0.99  0.17  -0.01  0.33   0.32  0.08     NA    NA
#> 12  CMC12     2   3PLM   0.86  0.13   1.30  0.18   0.11  0.04     NA    NA
#> 13  CMC13     2   3PLM   1.48  0.26   1.61  0.12   0.18  0.03     NA    NA
#> 14  CMC14     2   3PLM   1.53  0.21   0.25  0.16   0.27  0.05     NA    NA
#> 15  CMC15     2   3PLM   1.53  0.18  -0.11  0.13   0.14  0.05     NA    NA
#> 16  CMC16     2   3PLM   2.16  0.22   0.02  0.07   0.08  0.03     NA    NA
#> 17  CMC17     2   3PLM   1.39  0.19   0.03  0.17   0.20  0.06     NA    NA
#> 18  CMC18     2   3PLM   1.36  0.27   1.34  0.16   0.27  0.04     NA    NA
#> 19  CMC19     2   3PLM   2.48  0.37  -0.94  0.12   0.19  0.06     NA    NA
#> 20  CMC20     2   3PLM   1.80  0.37  -1.21  0.26   0.40  0.10     NA    NA
#> 21  CMC21     2   3PLM   1.76  0.22  -0.98  0.17   0.21  0.07     NA    NA
#> 22  CMC22     2   3PLM   0.94  0.13  -0.51  0.27   0.19  0.08     NA    NA
#> 23  CMC23     2   3PLM   0.83  0.10  -0.37  0.23   0.13  0.06     NA    NA
#> 24  CMC24     2   3PLM   0.98  0.21   1.86  0.20   0.22  0.04     NA    NA
#> 25  CMC25     2   3PLM   0.63  0.09  -2.01  0.47   0.21  0.09     NA    NA
#> 26  CMC26     2   3PLM   1.13  0.14  -1.68  0.28   0.22  0.09     NA    NA
#> 27  CMC27     2   3PLM   1.19  0.14   0.01  0.16   0.14  0.05     NA    NA
#> 28  CMC28     2   3PLM   2.23  0.26  -0.13  0.09   0.15  0.04     NA    NA
#> 29  CMC29     2   3PLM   1.31  0.16  -1.32  0.22   0.19  0.08     NA    NA
#> 30  CMC30     2   3PLM   1.63  0.30   1.03  0.15   0.37  0.04     NA    NA
#> 31  CMC31     2   3PLM   1.03  0.15   0.93  0.17   0.15  0.05     NA    NA
#> 32  CMC32     2   3PLM   1.55  0.21  -0.75  0.20   0.26  0.08     NA    NA
#> 33  CMC33     2   3PLM   1.24  0.19  -1.09  0.30   0.31  0.10     NA    NA
#> 34  CMC34     2   3PLM   1.34  0.16   0.31  0.15   0.17  0.05     NA    NA
#> 35  CMC35     2   3PLM   1.24  0.15  -0.36  0.20   0.19  0.07     NA    NA
#> 36  CMC36     2   3PLM   1.06  0.17   1.05  0.17   0.15  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.11  0.26  -0.29  0.11   0.16  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.57  0.11  -0.30  0.55   0.26  0.10     NA    NA
#> 39   CFR1     5    GRM   2.09  0.13  -1.81  0.10  -1.14  0.07  -0.68  0.06
#> 40   CFR2     5    GRM   1.38  0.08  -0.70  0.08  -0.08  0.07   0.48  0.06
#> 41   AMC1     2   3PLM   1.25  0.18   0.62  0.16   0.18  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.79  0.22  -1.61  0.18   0.17  0.07     NA    NA
#> 43   AMC3     2   3PLM   1.37  0.17   0.64  0.12   0.12  0.04     NA    NA
#> 44   AMC4     2   3PLM   0.94  0.11  -0.22  0.23   0.16  0.06     NA    NA
#> 45   AMC5     2   3PLM   1.11  0.33   2.83  0.26   0.21  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.22  0.37   1.70  0.09   0.19  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.16  0.13   0.02  0.14   0.10  0.04     NA    NA
#> 48   AMC8     2   3PLM   1.31  0.16   0.33  0.15   0.18  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.22  0.13   0.30  0.12   0.09  0.04     NA    NA
#> 50  AMC10     2   3PLM   1.83  0.28   1.48  0.09   0.15  0.03     NA    NA
#> 51  AMC11     2   3PLM   1.68  0.22  -1.08  0.17   0.19  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.91  0.13  -0.82  0.35   0.26  0.09     NA    NA
#> 53   AFR1     5    GRM   1.14    NA  -0.37    NA   0.22    NA   0.85    NA
#> 54   AFR2     5    GRM   1.23    NA  -2.08    NA  -1.35    NA  -0.71    NA
#> 55   AFR3     5    GRM   0.88    NA  -0.76    NA  -0.01    NA   0.67    NA
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.24  0.05
#> 40   1.05  0.07
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.38    NA
#> 54  -0.12    NA
#> 55   1.25    NA
#>  Group Parameters: 
#>              mu  sigma2  sigma
#> estimates  0.40    1.88   1.37
#> se         0.04    0.08   0.03
# plot the estimated empirical histogram of
# latent variable prior distribution
(emphist <- getirt(mod.fix1, what = "weights"))
#>    theta       weight
#> 1  -6.00 2.301252e-10
#> 2  -5.75 1.434595e-09
#> 3  -5.50 8.649281e-09
#> 4  -5.25 5.019868e-08
#> 5  -5.00 2.782912e-07
#> 6  -4.75 1.456751e-06
#> 7  -4.50 7.085630e-06
#> 8  -4.25 3.134808e-05
#> 9  -4.00 1.227012e-04
#> 10 -3.75 4.100603e-04
#> 11 -3.50 1.119742e-03
#> 12 -3.25 2.386354e-03
#> 13 -3.00 3.889321e-03
#> 14 -2.75 5.167720e-03
#> 15 -2.50 6.767224e-03
#> 16 -2.25 1.053445e-02
#> 17 -2.00 1.742567e-02
#> 18 -1.75 2.236173e-02
#> 19 -1.50 2.498735e-02
#> 20 -1.25 3.374293e-02
#> 21 -1.00 4.224442e-02
#> 22 -0.75 4.948939e-02
#> 23 -0.50 7.710828e-02
#> 24 -0.25 8.024288e-02
#> 25  0.00 4.752956e-02
#> 26  0.25 5.636581e-02
#> 27  0.50 8.674316e-02
#> 28  0.75 6.936382e-02
#> 29  1.00 6.035964e-02
#> 30  1.25 6.234864e-02
#> 31  1.50 5.768393e-02
#> 32  1.75 4.254007e-02
#> 33  2.00 3.148380e-02
#> 34  2.25 3.111071e-02
#> 35  2.50 2.935626e-02
#> 36  2.75 1.950494e-02
#> 37  3.00 1.019680e-02
#> 38  3.25 5.154557e-03
#> 39  3.50 2.816466e-03
#> 40  3.75 1.753617e-03
#> 41  4.00 1.282170e-03
#> 42  4.25 1.097172e-03
#> 43  4.50 1.050444e-03
#> 44  4.75 1.046478e-03
#> 45  5.00 1.004637e-03
#> 46  5.25 8.721994e-04
#> 47  5.50 6.557355e-04
#> 48  5.75 4.168380e-04
#> 49  6.00 2.220842e-04
plot(emphist$weight ~ emphist$theta, xlab = "Theta",
  ylab = "Density", type = "h")
```

<img src="man/figures/README-example-1.png" width="70%" height="50%" />

``` r
##----------------------------------------------------------------------------
# 2. The example code below shows how to prepare
# the data sets and how to estimate the item
# parameters using the FAPC:
##----------------------------------------------------------------------------

## Step 1: prepare a data set In this example, we
## generated examinees' true proficiency
## parameters and simulated the item response
## data using the function 'simdat'. Because the
## true proficiency parameters are not known in
## reality, the true proficiencies would be
## replaced with the proficiency estimates for
## the calibration.

# import the '-prm.txt' output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt",
  package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df

# modify the item metadata so that some items
# follow 1PLM, 2PLM and GPCM
x[c(1:3, 5), 3] <- "1PLM"
x[c(1:3, 5), 4] <- 1
x[c(1:3, 5), 6] <- 0
x[c(4, 8:12), 3] <- "2PLM"
x[c(4, 8:12), 6] <- 0
x[54:55, 3] <- "GPCM"

# generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean = 0, sd = 1)

# simulate the response data
data <- simdat(x = x, theta = score, D = 1)

## Step 2: Estimate the item parameters 1) item
## parameter estimation: constrain the slope
## parameters of the 1PLM to be equal
(mod1 <- est_item(x, data, score, D = 1, fix.a.1pl = FALSE,
  use.gprior = TRUE, gprior = list(dist = "beta",
    params = c(5, 17)), use.startval = FALSE))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished.
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         17)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15830.66
summary(mod1)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         17)), use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 162
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.96
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31661.31
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.02  0.06   1.60  0.13     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.02    NA  -1.06  0.12     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.02    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.02    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.88  0.17   1.03  0.23   0.13  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.35  0.34   1.31  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.36  0.31   0.15  0.24   0.24  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.53  0.27   0.01  0.17   0.20  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.15  -0.41  0.22   0.16  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.27  0.38   1.42  0.20   0.22  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.25  0.32  -1.11  0.14   0.17  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.74  0.22   0.18  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.21  -1.25  0.23   0.20  0.08     NA    NA
#> 22  CMC22     2   3PLM   0.92  0.16  -0.55  0.28   0.19  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.10  0.22  -0.12  0.27   0.22  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.21  0.34   1.43  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.51  0.40   0.21  0.09     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.16  0.35   0.19  0.08     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.09  0.17   0.14  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.19  0.31  -0.17  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.48  0.54  -0.81  0.20   0.38  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.88  0.45   0.69  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.70  0.16   1.00  0.32   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.73  0.30  -0.78  0.21   0.26  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.45  0.28   0.19  0.08     NA    NA
#> 34  CMC34     2   3PLM   1.04  0.17   0.21  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.36  0.19  -0.45  0.17   0.16  0.06     NA    NA
#> 36  CMC36     2   3PLM   0.88  0.17   0.97  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.25  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.17  -0.31  0.32   0.20  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.81  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.70  0.25  -1.59  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.30  0.25   0.68  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.94  0.17  -0.18  0.26   0.18  0.07     NA    NA
#> 45   AMC5     2   3PLM   1.69  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.83  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.69  0.41   0.37  0.18   0.25  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.29   0.39  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.55  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.48  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.73  0.23  -1.02  0.15   0.16  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.95  0.20  -0.83  0.38   0.24  0.10     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02
# 2) item parameter estimation: fix the slope
# parameters of the 1PLM to 1
(mod2 <- est_item(x, data, score, D = 1, fix.a.1pl = TRUE,
  a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta",
    params = c(5, 17)), use.startval = FALSE))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished.
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 17)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15830.7
summary(mod2)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 17)), use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 161
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.92
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31661.4
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.00    NA   1.62  0.12     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.00    NA  -1.08  0.11     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.00    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.00    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.88  0.17   1.03  0.23   0.13  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.35  0.34   1.31  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.36  0.31   0.15  0.24   0.24  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.53  0.27   0.01  0.17   0.20  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.15  -0.41  0.22   0.16  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.27  0.38   1.42  0.20   0.22  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.25  0.32  -1.11  0.14   0.17  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.74  0.22   0.18  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.21  -1.25  0.23   0.20  0.08     NA    NA
#> 22  CMC22     2   3PLM   0.92  0.16  -0.55  0.28   0.19  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.10  0.22  -0.12  0.27   0.22  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.21  0.34   1.43  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.51  0.40   0.21  0.09     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.16  0.35   0.19  0.08     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.09  0.17   0.14  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.19  0.31  -0.17  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.48  0.54  -0.81  0.20   0.38  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.88  0.45   0.69  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.70  0.16   1.00  0.32   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.73  0.30  -0.78  0.21   0.26  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.45  0.28   0.19  0.08     NA    NA
#> 34  CMC34     2   3PLM   1.04  0.17   0.21  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.36  0.19  -0.45  0.17   0.16  0.06     NA    NA
#> 36  CMC36     2   3PLM   0.88  0.17   0.97  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.25  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.17  -0.31  0.32   0.20  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.81  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.70  0.25  -1.59  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.30  0.25   0.68  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.94  0.17  -0.18  0.26   0.18  0.07     NA    NA
#> 45   AMC5     2   3PLM   1.69  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.83  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.69  0.41   0.37  0.18   0.25  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.29   0.39  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.55  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.48  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.73  0.23  -1.02  0.15   0.16  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.95  0.20  -0.83  0.38   0.24  0.10     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02
# 3) item parameter estimation: fix the guessing
# parameters of the 3PLM to 0.2
(mod3 <- est_item(x, data, score, D = 1, fix.a.1pl = TRUE,
  fix.g = TRUE, a.val.1pl = 1, g.val = 0.2, use.startval = FALSE))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished.
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     fix.g = TRUE, a.val.1pl = 1, g.val = 0.2, use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15916.26
summary(mod3)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     fix.g = TRUE, a.val.1pl = 1, g.val = 0.2, use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 121
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.72
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31832.52
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.00    NA   1.62  0.12     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.00    NA  -1.08  0.11     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.00    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.00    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   2.25  0.29   0.83  0.08   0.20    NA     NA    NA
#> 7    CMC7     2   3PLM   1.00  0.17   1.22  0.19   0.20    NA     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.53  0.28   1.37  0.15   0.20    NA     NA    NA
#> 14  CMC14     2   3PLM   1.26  0.18   0.05  0.10   0.20    NA     NA    NA
#> 15  CMC15     2   3PLM   1.52  0.20   0.00  0.09   0.20    NA     NA    NA
#> 16  CMC16     2   3PLM   2.36  0.27   0.18  0.07   0.20    NA     NA    NA
#> 17  CMC17     2   3PLM   1.06  0.15  -0.30  0.12   0.20    NA     NA    NA
#> 18  CMC18     2   3PLM   1.16  0.23   1.38  0.19   0.20    NA     NA    NA
#> 19  CMC19     2   3PLM   2.30  0.30  -1.07  0.10   0.20    NA     NA    NA
#> 20  CMC20     2   3PLM   1.48  0.22  -1.71  0.19   0.20    NA     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.19  -1.25  0.16   0.20    NA     NA    NA
#> 22  CMC22     2   3PLM   0.93  0.14  -0.52  0.15   0.20    NA     NA    NA
#> 23  CMC23     2   3PLM   1.08  0.16  -0.17  0.12   0.20    NA     NA    NA
#> 24  CMC24     2   3PLM   1.13  0.23   1.40  0.19   0.20    NA     NA    NA
#> 25  CMC25     2   3PLM   0.82  0.15  -1.55  0.28   0.20    NA     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.14  0.31   0.20    NA     NA    NA
#> 27  CMC27     2   3PLM   1.27  0.17   0.22  0.10   0.20    NA     NA    NA
#> 28  CMC28     2   3PLM   2.22  0.27  -0.15  0.07   0.20    NA     NA    NA
#> 29  CMC29     2   3PLM   1.84  0.28  -1.19  0.14   0.20    NA     NA    NA
#> 30  CMC30     2   3PLM   1.19  0.20   0.35  0.11   0.20    NA     NA    NA
#> 31  CMC31     2   3PLM   0.76  0.15   1.15  0.23   0.20    NA     NA    NA
#> 32  CMC32     2   3PLM   1.62  0.22  -0.90  0.12   0.20    NA     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.16  -1.43  0.21   0.20    NA     NA    NA
#> 34  CMC34     2   3PLM   1.11  0.16   0.33  0.12   0.20    NA     NA    NA
#> 35  CMC35     2   3PLM   1.42  0.18  -0.36  0.10   0.20    NA     NA    NA
#> 36  CMC36     2   3PLM   1.00  0.17   1.15  0.18   0.20    NA     NA    NA
#> 37  CMC37     2   3PLM   2.30  0.27  -0.15  0.07   0.20    NA     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.14  -0.33  0.15   0.20    NA     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.45  0.22   0.57  0.10   0.20    NA     NA    NA
#> 42   AMC2     2   3PLM   1.72  0.25  -1.57  0.16   0.20    NA     NA    NA
#> 43   AMC3     2   3PLM   1.44  0.21   0.77  0.11   0.20    NA     NA    NA
#> 44   AMC4     2   3PLM   0.96  0.15  -0.13  0.13   0.20    NA     NA    NA
#> 45   AMC5     2   3PLM   1.83  0.53   2.10  0.25   0.20    NA     NA    NA
#> 46   AMC6     2   3PLM   3.32  0.68   1.50  0.09   0.20    NA     NA    NA
#> 47   AMC7     2   3PLM   1.48  0.21   0.25  0.09   0.20    NA     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.23   0.39  0.09   0.20    NA     NA    NA
#> 49   AMC9     2   3PLM   1.72  0.23   0.59  0.09   0.20    NA     NA    NA
#> 50  AMC10     2   3PLM   3.21  0.62   1.40  0.09   0.20    NA     NA    NA
#> 51  AMC11     2   3PLM   1.78  0.22  -0.96  0.11   0.20    NA     NA    NA
#> 52  AMC12     2   3PLM   0.91  0.15  -0.96  0.19   0.20    NA     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02
##----------------------------------------------------------------------------
# 3. The example code below shows how to prepare
# the data sets and how to conduct the IRT
# model-data fit analysis:
##----------------------------------------------------------------------------

## Step 1: prepare a data set for IRT In this
## example, we use the simulated mixed-item
## format CAT Data But, only items that have item
## responses more than 1,000 are assessed.

# find the location of items that have more than
# 1,000 item responses
over1000 <- which(colSums(simCAT_MX$res.dat, na.rm = TRUE) >
  1000)

# (1) item metadata
x <- simCAT_MX$item.prm[over1000, ]
dim(x)
#> [1] 113   7
print(x[1:10, ])
#>     id cats model     par.1      par.2 par.3 par.4
#> 2   V2    2  2PLM 0.9152754  1.3843593    NA    NA
#> 3   V3    2  2PLM 1.3454796 -1.2554919    NA    NA
#> 5   V5    2  2PLM 1.0862914  1.7114409    NA    NA
#> 6   V6    2  2PLM 1.1311496 -0.6029080    NA    NA
#> 7   V7    2  2PLM 1.2012407 -0.4721664    NA    NA
#> 8   V8    2  2PLM 1.3244155 -0.6353713    NA    NA
#> 10 V10    2  2PLM 1.2487125  0.1381082    NA    NA
#> 11 V11    2  2PLM 1.4413208  1.2276303    NA    NA
#> 12 V12    2  2PLM 1.2077273 -0.8017795    NA    NA
#> 13 V13    2  2PLM 1.1715456 -1.0803926    NA    NA
# (2) examinee's ability estimates
score <- simCAT_MX$score
length(score)
#> [1] 30000
print(score[1:100])
#>   [1] -0.30311440 -0.67224807 -0.73474583  1.76935738 -0.91017203 -0.28448278
#>   [7]  0.81656431 -1.66434615  0.59312008 -0.35182937  0.23129679 -0.93107524
#>  [13] -0.29971993 -0.32700449 -0.22271651  1.48912121 -0.92927809  0.43453041
#>  [19] -0.01795450 -0.28365286  0.01115173 -0.76101441  0.12144273  0.83096135
#>  [25]  1.96600585 -0.83510402 -0.40268865 -0.05605526  0.72398446 -0.16026059
#>  [31] -1.09011778  1.22126764 -0.13340360 -1.28230720 -1.05581980  0.83484173
#>  [37] -0.52136360 -0.66913590 -1.08580804  1.73214834  0.56950387  0.48016332
#>  [43] -0.03472720 -2.17577824  0.44127032  0.98913071  1.43861714 -1.08133809
#>  [49] -0.69016072  0.19325797  0.89998383  1.25383167 -1.09600809  0.50519143
#>  [55] -0.51707395 -0.39474484 -0.45031102  1.85675021  1.50768131  1.06011811
#>  [61] -0.41064797  1.10960278 -0.68853387 -0.59397660 -0.65326436  0.29147751
#>  [67] -1.86787473  1.04838050 -1.14582092  1.07395234 -0.03828693  0.08445559
#>  [73]  0.34582524  0.72300905  0.84448992 -1.86488055  0.77121937  1.66573208
#>  [79]  0.10311673 -0.50768866 -1.60992457 -0.23074682  0.16162326  0.26091160
#>  [85]  0.60682182  0.65415304 -0.69923141  1.07545766  0.24060267 -0.93542383
#>  [91]  1.24988766 -0.01826940  1.27403936  0.10985621 -1.19092047  0.79614598
#>  [97]  0.62302338 -0.89455596 -0.03472720  0.20250837
# (3) response data
data <- simCAT_MX$res.dat[, over1000]
dim(data)
#> [1] 30000   113
print(data[1:20, 1:6])
#>       Item.dc.2 Item.dc.3 Item.dc.5 Item.dc.6 Item.dc.7 Item.dc.8
#>  [1,]        NA        NA        NA        NA         0         1
#>  [2,]        NA        NA        NA        NA         0         1
#>  [3,]        NA        NA        NA        NA         1         1
#>  [4,]        NA        NA         0        NA        NA        NA
#>  [5,]        NA         1        NA         0         1         1
#>  [6,]        NA         0        NA         0         1         0
#>  [7,]        NA        NA        NA        NA        NA        NA
#>  [8,]        NA         1        NA         1         1         0
#>  [9,]        NA        NA         0        NA        NA        NA
#> [10,]        NA         0        NA         1         1         0
#> [11,]        NA        NA         0        NA        NA        NA
#> [12,]        NA         1        NA         0         1         1
#> [13,]        NA         0        NA         1         1         0
#> [14,]        NA        NA        NA        NA         1        NA
#> [15,]        NA        NA        NA         1         1         1
#> [16,]         1        NA         0        NA        NA        NA
#> [17,]        NA         0        NA         0         1         0
#> [18,]        NA        NA        NA        NA        NA        NA
#> [19,]        NA         0        NA         0         1         1
#> [20,]        NA        NA        NA        NA        NA        NA
## Step 2: Compute the IRT mode-data fit
## statistics (1) the use of 'equal.width'
fit1 <- irtfit(x = x, score = score, data = data, group.method = "equal.width",
  n.width = 11, loc.theta = "average", range.score = c(-4,
    4), D = 1, alpha = 0.05, missing = NA, overSR = 2.5)

# what kinds of internal objects does the results
# have?
names(fit1)
#> [1] "fit_stat"            "contingency.fitstat" "contingency.plot"   
#> [4] "item_df"             "individual.info"     "ancillary"          
#> [7] "call"
# show the results of the fit statistics
fit1$fit_stat[1:10, ]
#>     id      X2      G2 df.X2 df.G2 crit.value.X2 crit.value.G2 p.value.X2
#> 1   V2  75.070  75.209     8    10         15.51         18.31          0
#> 2   V3 186.880 168.082     8    10         15.51         18.31          0
#> 3   V5 151.329 139.213     8    10         15.51         18.31          0
#> 4   V6 178.409 157.911     8    10         15.51         18.31          0
#> 5   V7 185.438 170.360     9    11         16.92         19.68          0
#> 6   V8 209.653 193.001     8    10         15.51         18.31          0
#> 7  V10 267.444 239.563     9    11         16.92         19.68          0
#> 8  V11 148.896 133.209     7     9         14.07         16.92          0
#> 9  V12 139.295 125.647     9    11         16.92         19.68          0
#> 10 V13 128.422 117.439     9    11         16.92         19.68          0
#>    p.value.G2 outfit infit     N overSR.prop
#> 1           0  1.018 1.016  2018       0.364
#> 2           0  1.124 1.090 11041       0.636
#> 3           0  1.133 1.111  5181       0.727
#> 4           0  1.056 1.045 13599       0.545
#> 5           0  1.078 1.059 18293       0.455
#> 6           0  1.098 1.075 16163       0.636
#> 7           0  1.097 1.073 19702       0.727
#> 8           0  1.129 1.083 13885       0.455
#> 9           0  1.065 1.051 12118       0.636
#> 10          0  1.075 1.059 10719       0.545
# show the contingency tables for the first item
# (dichotomous)
fit1$contingency.fitstat[[1]]
#>      N freq.0 freq.1 obs.prop.0 obs.prop.1 exp.prob.0 exp.prob.1 raw_resid.0
#> 1    8      5      3  0.6250000  0.3750000  0.7627914  0.2372086 -0.13779141
#> 2   14      8      6  0.5714286  0.4285714  0.7121079  0.2878921 -0.14067932
#> 3   60     34     26  0.5666667  0.4333333  0.6708959  0.3291041 -0.10422928
#> 4  185     99     86  0.5351351  0.4648649  0.6230537  0.3769463 -0.08791853
#> 5  240    115    125  0.4791667  0.5208333  0.5765337  0.4234663 -0.09736699
#> 6  349    145    204  0.4154728  0.5845272  0.5301760  0.4698240 -0.11470327
#> 7  325    114    211  0.3507692  0.6492308  0.4784096  0.5215904 -0.12764036
#> 8  246     82    164  0.3333333  0.6666667  0.4419993  0.5580007 -0.10866594
#> 9  377    139    238  0.3687003  0.6312997  0.4086532  0.5913468 -0.03995295
#> 10 214     78    136  0.3644860  0.6355140  0.3394647  0.6605353  0.02502128
#>    raw_resid.1
#> 1   0.13779141
#> 2   0.14067932
#> 3   0.10422928
#> 4   0.08791853
#> 5   0.09736699
#> 6   0.11470327
#> 7   0.12764036
#> 8   0.10866594
#> 9   0.03995295
#> 10 -0.02502128
# (2) the use of 'equal.freq'
fit2 <- irtfit(x = x, score = score, data = data, group.method = "equal.freq",
  n.width = 11, loc.theta = "average", range.score = c(-4,
    4), D = 1, alpha = 0.05, missing = NA)

# show the results of the fit statistics
fit2$fit_stat[1:10, ]
#>     id      X2      G2 df.X2 df.G2 crit.value.X2 crit.value.G2 p.value.X2
#> 1   V2  77.967  78.144     9    11         16.92         19.68          0
#> 2   V3 202.035 181.832     9    11         16.92         19.68          0
#> 3   V5 146.383 135.908     9    11         16.92         19.68          0
#> 4   V6 140.038 133.287     9    11         16.92         19.68          0
#> 5   V7 188.814 177.526     9    11         16.92         19.68          0
#> 6   V8 211.279 196.328     9    11         16.92         19.68          0
#> 7  V10 259.669 239.292     9    11         16.92         19.68          0
#> 8  V11 166.427 150.419     9    11         16.92         19.68          0
#> 9  V12 145.789 134.690     9    11         16.92         19.68          0
#> 10 V13 141.283 132.270     9    11         16.92         19.68          0
#>    p.value.G2 outfit infit     N overSR.prop
#> 1           0  1.018 1.016  2018       0.727
#> 2           0  1.124 1.090 11041       0.636
#> 3           0  1.133 1.111  5181       0.727
#> 4           0  1.056 1.045 13599       0.545
#> 5           0  1.078 1.059 18293       0.455
#> 6           0  1.098 1.075 16163       0.545
#> 7           0  1.097 1.073 19702       0.636
#> 8           0  1.129 1.083 13885       0.636
#> 9           0  1.065 1.051 12118       0.364
#> 10          0  1.075 1.059 10719       0.636
# show the contingency table for the fourth item
# (polytomous)
fit2$contingency.fitstat[[4]]
#>       N freq.0 freq.1 obs.prop.0 obs.prop.1 exp.prob.0 exp.prob.1  raw_resid.0
#> 1  1241    967    274  0.7792103  0.2207897  0.8038510  0.1961490 -0.024640641
#> 2  1243    879    364  0.7071601  0.2928399  0.7161793  0.2838207 -0.009019180
#> 3  1243    784    459  0.6307321  0.3692679  0.6575849  0.3424151 -0.026852795
#> 4  1219    747    472  0.6127974  0.3872026  0.6049393  0.3950607  0.007858099
#> 5  1236    705    531  0.5703883  0.4296117  0.5613454  0.4386546  0.009042942
#> 6  1243    677    566  0.5446500  0.4553500  0.5279560  0.4720440  0.016694048
#> 7  1270    662    608  0.5212598  0.4787402  0.4925592  0.5074408  0.028700633
#> 8  1230    616    614  0.5008130  0.4991870  0.4491759  0.5508241  0.051637085
#> 9  1207    553    654  0.4581607  0.5418393  0.4027790  0.5972210  0.055381721
#> 10 1233    494    739  0.4006488  0.5993512  0.3509261  0.6490739  0.049722759
#> 11 1234    465    769  0.3768233  0.6231767  0.2630181  0.7369819  0.113805214
#>     raw_resid.1
#> 1   0.024640641
#> 2   0.009019180
#> 3   0.026852795
#> 4  -0.007858099
#> 5  -0.009042942
#> 6  -0.016694048
#> 7  -0.028700633
#> 8  -0.051637085
#> 9  -0.055381721
#> 10 -0.049722759
#> 11 -0.113805214
## Step 3: Draw the IRT residual plots 1. the
## dichotomous item (1) both raw and standardized
## residual plots using the object 'fit1'
plot(x = fit1, item.loc = 1, type = "both", ci.method = "wald",
  ylim.sr.adjust = TRUE)
```

<img src="man/figures/README-example-2.png" width="70%" height="50%" />

    #>                               theta   N freq.0 freq.1 obs.prop.0 obs.prop.1
    #> [-0.1218815,0.08512996] -0.02529272   3      3      0  1.0000000  0.0000000
    #> (0.08512996,0.2921415]   0.18431014   5      2      3  0.4000000  0.6000000
    #> (0.2921415,0.499153]     0.39488272  14      8      6  0.5714286  0.4285714
    #> (0.499153,0.7061645]     0.60618911  60     34     26  0.5666667  0.4333333
    #> (0.7061645,0.913176]     0.83531169 185     99     86  0.5351351  0.4648649
    #> (0.913176,1.120187]      1.04723712 240    115    125  0.4791667  0.5208333
    #> (1.120187,1.327199]      1.25232143 349    145    204  0.4154728  0.5845272
    #> (1.327199,1.53421]       1.47877397 325    114    211  0.3507692  0.6492308
    #> (1.53421,1.741222]       1.63898436 246     82    164  0.3333333  0.6666667
    #> (1.741222,1.948233]      1.78810197 377    139    238  0.3687003  0.6312997
    #> (1.948233,2.155245]      2.11166019 214     78    136  0.3644860  0.6355140
    #>                         exp.prob.0 exp.prob.1 raw_resid.0 raw_resid.1
    #> [-0.1218815,0.08512996]  0.7841844  0.2158156  0.21581559 -0.21581559
    #> (0.08512996,0.2921415]   0.7499556  0.2500444 -0.34995561  0.34995561
    #> (0.2921415,0.499153]     0.7121079  0.2878921 -0.14067932  0.14067932
    #> (0.499153,0.7061645]     0.6708959  0.3291041 -0.10422928  0.10422928
    #> (0.7061645,0.913176]     0.6230537  0.3769463 -0.08791853  0.08791853
    #> (0.913176,1.120187]      0.5765337  0.4234663 -0.09736699  0.09736699
    #> (1.120187,1.327199]      0.5301760  0.4698240 -0.11470327  0.11470327
    #> (1.327199,1.53421]       0.4784096  0.5215904 -0.12764036  0.12764036
    #> (1.53421,1.741222]       0.4419993  0.5580007 -0.10866594  0.10866594
    #> (1.741222,1.948233]      0.4086532  0.5913468 -0.03995295  0.03995295
    #> (1.948233,2.155245]      0.3394647  0.6605353  0.02502128 -0.02502128
    #>                               se.0       se.1 std_resid.0 std_resid.1
    #> [-0.1218815,0.08512996] 0.23751437 0.23751437   0.9086423  -0.9086423
    #> (0.08512996,0.2921415]  0.19366063 0.19366063  -1.8070560   1.8070560
    #> (0.2921415,0.499153]    0.12101070 0.12101070  -1.1625362   1.1625362
    #> (0.499153,0.7061645]    0.06066226 0.06066226  -1.7181899   1.7181899
    #> (0.7061645,0.913176]    0.03563007 0.03563007  -2.4675377   2.4675377
    #> (0.913176,1.120187]     0.03189453 0.03189453  -3.0527806   3.0527806
    #> (1.120187,1.327199]     0.02671560 0.02671560  -4.2934941   4.2934941
    #> (1.327199,1.53421]      0.02770914 0.02770914  -4.6064351   4.6064351
    #> (1.53421,1.741222]      0.03166362 0.03166362  -3.4318859   3.4318859
    #> (1.741222,1.948233]     0.02531791 0.02531791  -1.5780508   1.5780508
    #> (1.948233,2.155245]     0.03236968 0.03236968   0.7729850  -0.7729850
    #>                         raw.resid.0 raw.resid.1     se.0.1     se.1.1
    #> [-0.1218815,0.08512996]  0.21581559 -0.21581559 0.23751437 0.23751437
    #> (0.08512996,0.2921415]  -0.34995561  0.34995561 0.19366063 0.19366063
    #> (0.2921415,0.499153]    -0.14067932  0.14067932 0.12101070 0.12101070
    #> (0.499153,0.7061645]    -0.10422928  0.10422928 0.06066226 0.06066226
    #> (0.7061645,0.913176]    -0.08791853  0.08791853 0.03563007 0.03563007
    #> (0.913176,1.120187]     -0.09736699  0.09736699 0.03189453 0.03189453
    #> (1.120187,1.327199]     -0.11470327  0.11470327 0.02671560 0.02671560
    #> (1.327199,1.53421]      -0.12764036  0.12764036 0.02770914 0.02770914
    #> (1.53421,1.741222]      -0.10866594  0.10866594 0.03166362 0.03166362
    #> (1.741222,1.948233]     -0.03995295  0.03995295 0.02531791 0.02531791
    #> (1.948233,2.155245]      0.02502128 -0.02502128 0.03236968 0.03236968
    #>                         std.resid.0 std.resid.1
    #> [-0.1218815,0.08512996]   0.9086423  -0.9086423
    #> (0.08512996,0.2921415]   -1.8070560   1.8070560
    #> (0.2921415,0.499153]     -1.1625362   1.1625362
    #> (0.499153,0.7061645]     -1.7181899   1.7181899
    #> (0.7061645,0.913176]     -2.4675377   2.4675377
    #> (0.913176,1.120187]      -3.0527806   3.0527806
    #> (1.120187,1.327199]      -4.2934941   4.2934941
    #> (1.327199,1.53421]       -4.6064351   4.6064351
    #> (1.53421,1.741222]       -3.4318859   3.4318859
    #> (1.741222,1.948233]      -1.5780508   1.5780508
    #> (1.948233,2.155245]       0.7729850  -0.7729850
    # (2) the raw residual plots using the object
    # 'fit1'
    plot(x = fit1, item.loc = 1, type = "icc", ci.method = "wald",
      ylim.sr.adjust = TRUE)

<img src="man/figures/README-example-3.png" width="70%" height="50%" />

    #>                               theta   N freq.0 freq.1 obs.prop.0 obs.prop.1
    #> [-0.1218815,0.08512996] -0.02529272   3      3      0  1.0000000  0.0000000
    #> (0.08512996,0.2921415]   0.18431014   5      2      3  0.4000000  0.6000000
    #> (0.2921415,0.499153]     0.39488272  14      8      6  0.5714286  0.4285714
    #> (0.499153,0.7061645]     0.60618911  60     34     26  0.5666667  0.4333333
    #> (0.7061645,0.913176]     0.83531169 185     99     86  0.5351351  0.4648649
    #> (0.913176,1.120187]      1.04723712 240    115    125  0.4791667  0.5208333
    #> (1.120187,1.327199]      1.25232143 349    145    204  0.4154728  0.5845272
    #> (1.327199,1.53421]       1.47877397 325    114    211  0.3507692  0.6492308
    #> (1.53421,1.741222]       1.63898436 246     82    164  0.3333333  0.6666667
    #> (1.741222,1.948233]      1.78810197 377    139    238  0.3687003  0.6312997
    #> (1.948233,2.155245]      2.11166019 214     78    136  0.3644860  0.6355140
    #>                         exp.prob.0 exp.prob.1 raw_resid.0 raw_resid.1
    #> [-0.1218815,0.08512996]  0.7841844  0.2158156  0.21581559 -0.21581559
    #> (0.08512996,0.2921415]   0.7499556  0.2500444 -0.34995561  0.34995561
    #> (0.2921415,0.499153]     0.7121079  0.2878921 -0.14067932  0.14067932
    #> (0.499153,0.7061645]     0.6708959  0.3291041 -0.10422928  0.10422928
    #> (0.7061645,0.913176]     0.6230537  0.3769463 -0.08791853  0.08791853
    #> (0.913176,1.120187]      0.5765337  0.4234663 -0.09736699  0.09736699
    #> (1.120187,1.327199]      0.5301760  0.4698240 -0.11470327  0.11470327
    #> (1.327199,1.53421]       0.4784096  0.5215904 -0.12764036  0.12764036
    #> (1.53421,1.741222]       0.4419993  0.5580007 -0.10866594  0.10866594
    #> (1.741222,1.948233]      0.4086532  0.5913468 -0.03995295  0.03995295
    #> (1.948233,2.155245]      0.3394647  0.6605353  0.02502128 -0.02502128
    #>                               se.0       se.1 std_resid.0 std_resid.1
    #> [-0.1218815,0.08512996] 0.23751437 0.23751437   0.9086423  -0.9086423
    #> (0.08512996,0.2921415]  0.19366063 0.19366063  -1.8070560   1.8070560
    #> (0.2921415,0.499153]    0.12101070 0.12101070  -1.1625362   1.1625362
    #> (0.499153,0.7061645]    0.06066226 0.06066226  -1.7181899   1.7181899
    #> (0.7061645,0.913176]    0.03563007 0.03563007  -2.4675377   2.4675377
    #> (0.913176,1.120187]     0.03189453 0.03189453  -3.0527806   3.0527806
    #> (1.120187,1.327199]     0.02671560 0.02671560  -4.2934941   4.2934941
    #> (1.327199,1.53421]      0.02770914 0.02770914  -4.6064351   4.6064351
    #> (1.53421,1.741222]      0.03166362 0.03166362  -3.4318859   3.4318859
    #> (1.741222,1.948233]     0.02531791 0.02531791  -1.5780508   1.5780508
    #> (1.948233,2.155245]     0.03236968 0.03236968   0.7729850  -0.7729850
    #>                         raw.resid.0 raw.resid.1     se.0.1     se.1.1
    #> [-0.1218815,0.08512996]  0.21581559 -0.21581559 0.23751437 0.23751437
    #> (0.08512996,0.2921415]  -0.34995561  0.34995561 0.19366063 0.19366063
    #> (0.2921415,0.499153]    -0.14067932  0.14067932 0.12101070 0.12101070
    #> (0.499153,0.7061645]    -0.10422928  0.10422928 0.06066226 0.06066226
    #> (0.7061645,0.913176]    -0.08791853  0.08791853 0.03563007 0.03563007
    #> (0.913176,1.120187]     -0.09736699  0.09736699 0.03189453 0.03189453
    #> (1.120187,1.327199]     -0.11470327  0.11470327 0.02671560 0.02671560
    #> (1.327199,1.53421]      -0.12764036  0.12764036 0.02770914 0.02770914
    #> (1.53421,1.741222]      -0.10866594  0.10866594 0.03166362 0.03166362
    #> (1.741222,1.948233]     -0.03995295  0.03995295 0.02531791 0.02531791
    #> (1.948233,2.155245]      0.02502128 -0.02502128 0.03236968 0.03236968
    #>                         std.resid.0 std.resid.1
    #> [-0.1218815,0.08512996]   0.9086423  -0.9086423
    #> (0.08512996,0.2921415]   -1.8070560   1.8070560
    #> (0.2921415,0.499153]     -1.1625362   1.1625362
    #> (0.499153,0.7061645]     -1.7181899   1.7181899
    #> (0.7061645,0.913176]     -2.4675377   2.4675377
    #> (0.913176,1.120187]      -3.0527806   3.0527806
    #> (1.120187,1.327199]      -4.2934941   4.2934941
    #> (1.327199,1.53421]       -4.6064351   4.6064351
    #> (1.53421,1.741222]       -3.4318859   3.4318859
    #> (1.741222,1.948233]      -1.5780508   1.5780508
    #> (1.948233,2.155245]       0.7729850  -0.7729850
    # (3) the standardized residual plots using the
    # object 'fit1'
    plot(x = fit1, item.loc = 113, type = "sr", ci.method = "wald",
      ylim.sr.adjust = TRUE)

<img src="man/figures/README-example-4.png" width="70%" height="50%" />

    #>                           theta   N freq.0 freq.1 freq.2 freq.3 obs.prop.0
    #> [0.3564295,0.5199582] 0.3564295   1      1      0      0      0 1.00000000
    #> (0.5199582,0.6834869] 0.6081321   5      3      2      0      0 0.60000000
    #> (0.6834869,0.8470155] 0.7400138  15      5     10      0      0 0.33333333
    #> (0.8470155,1.010544]  0.8866202  55      5     15     34      1 0.09090909
    #> (1.010544,1.174073]   1.0821064 133      6     40     53     34 0.04511278
    #> (1.174073,1.337602]   1.2832293 260      8     37    153     62 0.03076923
    #> (1.337602,1.50113]    1.4747336  98      0     23     57     18 0.00000000
    #> (1.50113,1.664659]    1.5311735 306      0      7     85    214 0.00000000
    #> (1.664659,1.828188]   1.7632607 418      0      0    145    273 0.00000000
    #> (1.828188,1.991716]   1.8577191  69      0      0      0     69 0.00000000
    #> (1.991716,2.155245]   2.1021956 263      0      0      0    263 0.00000000
    #>                       obs.prop.1 obs.prop.2 obs.prop.3  exp.prob.0 exp.prob.1
    #> [0.3564295,0.5199582] 0.00000000  0.0000000 0.00000000 0.196511833 0.31299790
    #> (0.5199582,0.6834869] 0.40000000  0.0000000 0.00000000 0.130431315 0.27025065
    #> (0.6834869,0.8470155] 0.66666667  0.0000000 0.00000000 0.102631449 0.24407213
    #> (0.8470155,1.010544]  0.27272727  0.6181818 0.01818182 0.077129446 0.21379299
    #> (1.010544,1.174073]   0.30075188  0.3984962 0.25563910 0.051169395 0.17398130
    #> (1.174073,1.337602]   0.14230769  0.5884615 0.23846154 0.032494074 0.13632441
    #> (1.337602,1.50113]    0.23469388  0.5816327 0.18367347 0.020534959 0.10523862
    #> (1.50113,1.664659]    0.02287582  0.2777778 0.69934641 0.017857642 0.09707782
    #> (1.664659,1.828188]   0.00000000  0.3468900 0.65311005 0.009866868 0.06836048
    #> (1.828188,1.991716]   0.00000000  0.0000000 1.00000000 0.007690015 0.05880602
    #> (1.991716,2.155245]   0.00000000  0.0000000 1.00000000 0.003962699 0.03912356
    #>                       exp.prob.2 exp.prob.3  raw_resid.0  raw_resid.1
    #> [0.3564295,0.5199582]  0.3472692  0.1432210  0.803488167 -0.312997903
    #> (0.5199582,0.6834869]  0.3900531  0.2092649  0.469568685  0.129749354
    #> (0.6834869,0.8470155]  0.4043226  0.2489738  0.230701885  0.422594537
    #> (0.8470155,1.010544]   0.4127992  0.2962784  0.013779645  0.058934282
    #> (1.010544,1.174073]    0.4120662  0.3627831 -0.006056613  0.126770584
    #> (1.174073,1.337602]    0.3983962  0.4327853 -0.001724843  0.005983284
    #> (1.337602,1.50113]     0.3756891  0.4985373 -0.020534959  0.129455259
    #> (1.50113,1.664659]     0.3676106  0.5174539 -0.017857642 -0.074201998
    #> (1.664659,1.828188]    0.3299158  0.5918569 -0.009866868 -0.068360484
    #> (1.828188,1.991716]    0.3132481  0.6202558 -0.007690015 -0.058806016
    #> (1.991716,2.155245]    0.2690653  0.6878484 -0.003962699 -0.039123555
    #>                       raw_resid.2 raw_resid.3        se.0       se.1       se.2
    #> [0.3564295,0.5199582] -0.34726922 -0.14322104 0.397359953 0.46371351 0.47610220
    #> (0.5199582,0.6834869] -0.39005312 -0.20926492 0.150611412 0.19860274 0.21813376
    #> (0.6834869,0.8470155] -0.40432265 -0.24897378 0.078357401 0.11090564 0.12671381
    #> (0.8470155,1.010544]   0.20538261 -0.27809653 0.035974863 0.05528201 0.06638675
    #> (1.010544,1.174073]   -0.01356995 -0.10714402 0.019106171 0.03287157 0.04267975
    #> (1.174073,1.337602]    0.19006531 -0.19432375 0.010996190 0.02128019 0.03036171
    #> (1.337602,1.50113]     0.20594357 -0.31486387 0.014326112 0.03099761 0.04892172
    #> (1.50113,1.664659]    -0.08983281  0.18189246 0.007570744 0.01692484 0.02756294
    #> (1.664659,1.828188]    0.01697418  0.06125317 0.004834464 0.01234350 0.02299737
    #> (1.828188,1.991716]   -0.31324812  0.37974416 0.010516294 0.02832213 0.05583668
    #> (1.991716,2.155245]   -0.26906531  0.31215156 0.003873963 0.01195570 0.02734578
    #>                             se.3 std_resid.0 std_resid.1 std_resid.2
    #> [0.3564295,0.5199582] 0.35029812   2.0220663  -0.6749812  -0.7294006
    #> (0.5199582,0.6834869] 0.18191928   3.1177497   0.6533110  -1.7881373
    #> (0.6834869,0.8470155] 0.11165000   2.9442258   3.8103971  -3.1908333
    #> (0.8470155,1.010544]  0.06156999   0.3830354   1.0660662   3.0937290
    #> (1.010544,1.174073]   0.04169091  -0.3169978   3.8565422  -0.3179482
    #> (1.174073,1.337602]   0.03072722  -0.1568583   0.2811669   6.2600333
    #> (1.337602,1.50113]    0.05050741  -1.4333937   4.1762987   4.2096551
    #> (1.50113,1.664659]    0.02856568  -2.3587698  -4.3842081  -3.2591881
    #> (1.664659,1.828188]   0.02403956  -2.0409436  -5.5381760   0.7380923
    #> (1.828188,1.991716]   0.05842604  -0.7312476  -2.0763275  -5.6100775
    #> (1.991716,2.155245]   0.02857270  -1.0229057  -3.2723764  -9.8393733
    #>                       std_resid.3  raw.resid.0  raw.resid.1 raw.resid.2
    #> [0.3564295,0.5199582]  -0.4088547  0.803488167 -0.312997903 -0.34726922
    #> (0.5199582,0.6834869]  -1.1503175  0.469568685  0.129749354 -0.39005312
    #> (0.6834869,0.8470155]  -2.2299488  0.230701885  0.422594537 -0.40432265
    #> (0.8470155,1.010544]   -4.5167547  0.013779645  0.058934282  0.20538261
    #> (1.010544,1.174073]    -2.5699613 -0.006056613  0.126770584 -0.01356995
    #> (1.174073,1.337602]    -6.3241558 -0.001724843  0.005983284  0.19006531
    #> (1.337602,1.50113]     -6.2340132 -0.020534959  0.129455259  0.20594357
    #> (1.50113,1.664659]      6.3675177 -0.017857642 -0.074201998 -0.08983281
    #> (1.664659,1.828188]     2.5480159 -0.009866868 -0.068360484  0.01697418
    #> (1.828188,1.991716]     6.4995706 -0.007690015 -0.058806016 -0.31324812
    #> (1.991716,2.155245]    10.9248190 -0.003962699 -0.039123555 -0.26906531
    #>                       raw.resid.3      se.0.1     se.1.1     se.2.1     se.3.1
    #> [0.3564295,0.5199582] -0.14322104 0.397359953 0.46371351 0.47610220 0.35029812
    #> (0.5199582,0.6834869] -0.20926492 0.150611412 0.19860274 0.21813376 0.18191928
    #> (0.6834869,0.8470155] -0.24897378 0.078357401 0.11090564 0.12671381 0.11165000
    #> (0.8470155,1.010544]  -0.27809653 0.035974863 0.05528201 0.06638675 0.06156999
    #> (1.010544,1.174073]   -0.10714402 0.019106171 0.03287157 0.04267975 0.04169091
    #> (1.174073,1.337602]   -0.19432375 0.010996190 0.02128019 0.03036171 0.03072722
    #> (1.337602,1.50113]    -0.31486387 0.014326112 0.03099761 0.04892172 0.05050741
    #> (1.50113,1.664659]     0.18189246 0.007570744 0.01692484 0.02756294 0.02856568
    #> (1.664659,1.828188]    0.06125317 0.004834464 0.01234350 0.02299737 0.02403956
    #> (1.828188,1.991716]    0.37974416 0.010516294 0.02832213 0.05583668 0.05842604
    #> (1.991716,2.155245]    0.31215156 0.003873963 0.01195570 0.02734578 0.02857270
    #>                       std.resid.0 std.resid.1 std.resid.2 std.resid.3
    #> [0.3564295,0.5199582]   2.0220663  -0.6749812  -0.7294006  -0.4088547
    #> (0.5199582,0.6834869]   3.1177497   0.6533110  -1.7881373  -1.1503175
    #> (0.6834869,0.8470155]   2.9442258   3.8103971  -3.1908333  -2.2299488
    #> (0.8470155,1.010544]    0.3830354   1.0660662   3.0937290  -4.5167547
    #> (1.010544,1.174073]    -0.3169978   3.8565422  -0.3179482  -2.5699613
    #> (1.174073,1.337602]    -0.1568583   0.2811669   6.2600333  -6.3241558
    #> (1.337602,1.50113]     -1.4333937   4.1762987   4.2096551  -6.2340132
    #> (1.50113,1.664659]     -2.3587698  -4.3842081  -3.2591881   6.3675177
    #> (1.664659,1.828188]    -2.0409436  -5.5381760   0.7380923   2.5480159
    #> (1.828188,1.991716]    -0.7312476  -2.0763275  -5.6100775   6.4995706
    #> (1.991716,2.155245]    -1.0229057  -3.2723764  -9.8393733  10.9248190
    # 2. the polytomous item (1) both raw and
    # standardized residual plots using the object
    # 'fit1'
    plot(x = fit1, item.loc = 113, type = "both", ci.method = "wald",
      ylim.sr.adjust = TRUE)

<img src="man/figures/README-example-5.png" width="70%" height="50%" />

    #>                           theta   N freq.0 freq.1 freq.2 freq.3 obs.prop.0
    #> [0.3564295,0.5199582] 0.3564295   1      1      0      0      0 1.00000000
    #> (0.5199582,0.6834869] 0.6081321   5      3      2      0      0 0.60000000
    #> (0.6834869,0.8470155] 0.7400138  15      5     10      0      0 0.33333333
    #> (0.8470155,1.010544]  0.8866202  55      5     15     34      1 0.09090909
    #> (1.010544,1.174073]   1.0821064 133      6     40     53     34 0.04511278
    #> (1.174073,1.337602]   1.2832293 260      8     37    153     62 0.03076923
    #> (1.337602,1.50113]    1.4747336  98      0     23     57     18 0.00000000
    #> (1.50113,1.664659]    1.5311735 306      0      7     85    214 0.00000000
    #> (1.664659,1.828188]   1.7632607 418      0      0    145    273 0.00000000
    #> (1.828188,1.991716]   1.8577191  69      0      0      0     69 0.00000000
    #> (1.991716,2.155245]   2.1021956 263      0      0      0    263 0.00000000
    #>                       obs.prop.1 obs.prop.2 obs.prop.3  exp.prob.0 exp.prob.1
    #> [0.3564295,0.5199582] 0.00000000  0.0000000 0.00000000 0.196511833 0.31299790
    #> (0.5199582,0.6834869] 0.40000000  0.0000000 0.00000000 0.130431315 0.27025065
    #> (0.6834869,0.8470155] 0.66666667  0.0000000 0.00000000 0.102631449 0.24407213
    #> (0.8470155,1.010544]  0.27272727  0.6181818 0.01818182 0.077129446 0.21379299
    #> (1.010544,1.174073]   0.30075188  0.3984962 0.25563910 0.051169395 0.17398130
    #> (1.174073,1.337602]   0.14230769  0.5884615 0.23846154 0.032494074 0.13632441
    #> (1.337602,1.50113]    0.23469388  0.5816327 0.18367347 0.020534959 0.10523862
    #> (1.50113,1.664659]    0.02287582  0.2777778 0.69934641 0.017857642 0.09707782
    #> (1.664659,1.828188]   0.00000000  0.3468900 0.65311005 0.009866868 0.06836048
    #> (1.828188,1.991716]   0.00000000  0.0000000 1.00000000 0.007690015 0.05880602
    #> (1.991716,2.155245]   0.00000000  0.0000000 1.00000000 0.003962699 0.03912356
    #>                       exp.prob.2 exp.prob.3  raw_resid.0  raw_resid.1
    #> [0.3564295,0.5199582]  0.3472692  0.1432210  0.803488167 -0.312997903
    #> (0.5199582,0.6834869]  0.3900531  0.2092649  0.469568685  0.129749354
    #> (0.6834869,0.8470155]  0.4043226  0.2489738  0.230701885  0.422594537
    #> (0.8470155,1.010544]   0.4127992  0.2962784  0.013779645  0.058934282
    #> (1.010544,1.174073]    0.4120662  0.3627831 -0.006056613  0.126770584
    #> (1.174073,1.337602]    0.3983962  0.4327853 -0.001724843  0.005983284
    #> (1.337602,1.50113]     0.3756891  0.4985373 -0.020534959  0.129455259
    #> (1.50113,1.664659]     0.3676106  0.5174539 -0.017857642 -0.074201998
    #> (1.664659,1.828188]    0.3299158  0.5918569 -0.009866868 -0.068360484
    #> (1.828188,1.991716]    0.3132481  0.6202558 -0.007690015 -0.058806016
    #> (1.991716,2.155245]    0.2690653  0.6878484 -0.003962699 -0.039123555
    #>                       raw_resid.2 raw_resid.3        se.0       se.1       se.2
    #> [0.3564295,0.5199582] -0.34726922 -0.14322104 0.397359953 0.46371351 0.47610220
    #> (0.5199582,0.6834869] -0.39005312 -0.20926492 0.150611412 0.19860274 0.21813376
    #> (0.6834869,0.8470155] -0.40432265 -0.24897378 0.078357401 0.11090564 0.12671381
    #> (0.8470155,1.010544]   0.20538261 -0.27809653 0.035974863 0.05528201 0.06638675
    #> (1.010544,1.174073]   -0.01356995 -0.10714402 0.019106171 0.03287157 0.04267975
    #> (1.174073,1.337602]    0.19006531 -0.19432375 0.010996190 0.02128019 0.03036171
    #> (1.337602,1.50113]     0.20594357 -0.31486387 0.014326112 0.03099761 0.04892172
    #> (1.50113,1.664659]    -0.08983281  0.18189246 0.007570744 0.01692484 0.02756294
    #> (1.664659,1.828188]    0.01697418  0.06125317 0.004834464 0.01234350 0.02299737
    #> (1.828188,1.991716]   -0.31324812  0.37974416 0.010516294 0.02832213 0.05583668
    #> (1.991716,2.155245]   -0.26906531  0.31215156 0.003873963 0.01195570 0.02734578
    #>                             se.3 std_resid.0 std_resid.1 std_resid.2
    #> [0.3564295,0.5199582] 0.35029812   2.0220663  -0.6749812  -0.7294006
    #> (0.5199582,0.6834869] 0.18191928   3.1177497   0.6533110  -1.7881373
    #> (0.6834869,0.8470155] 0.11165000   2.9442258   3.8103971  -3.1908333
    #> (0.8470155,1.010544]  0.06156999   0.3830354   1.0660662   3.0937290
    #> (1.010544,1.174073]   0.04169091  -0.3169978   3.8565422  -0.3179482
    #> (1.174073,1.337602]   0.03072722  -0.1568583   0.2811669   6.2600333
    #> (1.337602,1.50113]    0.05050741  -1.4333937   4.1762987   4.2096551
    #> (1.50113,1.664659]    0.02856568  -2.3587698  -4.3842081  -3.2591881
    #> (1.664659,1.828188]   0.02403956  -2.0409436  -5.5381760   0.7380923
    #> (1.828188,1.991716]   0.05842604  -0.7312476  -2.0763275  -5.6100775
    #> (1.991716,2.155245]   0.02857270  -1.0229057  -3.2723764  -9.8393733
    #>                       std_resid.3  raw.resid.0  raw.resid.1 raw.resid.2
    #> [0.3564295,0.5199582]  -0.4088547  0.803488167 -0.312997903 -0.34726922
    #> (0.5199582,0.6834869]  -1.1503175  0.469568685  0.129749354 -0.39005312
    #> (0.6834869,0.8470155]  -2.2299488  0.230701885  0.422594537 -0.40432265
    #> (0.8470155,1.010544]   -4.5167547  0.013779645  0.058934282  0.20538261
    #> (1.010544,1.174073]    -2.5699613 -0.006056613  0.126770584 -0.01356995
    #> (1.174073,1.337602]    -6.3241558 -0.001724843  0.005983284  0.19006531
    #> (1.337602,1.50113]     -6.2340132 -0.020534959  0.129455259  0.20594357
    #> (1.50113,1.664659]      6.3675177 -0.017857642 -0.074201998 -0.08983281
    #> (1.664659,1.828188]     2.5480159 -0.009866868 -0.068360484  0.01697418
    #> (1.828188,1.991716]     6.4995706 -0.007690015 -0.058806016 -0.31324812
    #> (1.991716,2.155245]    10.9248190 -0.003962699 -0.039123555 -0.26906531
    #>                       raw.resid.3      se.0.1     se.1.1     se.2.1     se.3.1
    #> [0.3564295,0.5199582] -0.14322104 0.397359953 0.46371351 0.47610220 0.35029812
    #> (0.5199582,0.6834869] -0.20926492 0.150611412 0.19860274 0.21813376 0.18191928
    #> (0.6834869,0.8470155] -0.24897378 0.078357401 0.11090564 0.12671381 0.11165000
    #> (0.8470155,1.010544]  -0.27809653 0.035974863 0.05528201 0.06638675 0.06156999
    #> (1.010544,1.174073]   -0.10714402 0.019106171 0.03287157 0.04267975 0.04169091
    #> (1.174073,1.337602]   -0.19432375 0.010996190 0.02128019 0.03036171 0.03072722
    #> (1.337602,1.50113]    -0.31486387 0.014326112 0.03099761 0.04892172 0.05050741
    #> (1.50113,1.664659]     0.18189246 0.007570744 0.01692484 0.02756294 0.02856568
    #> (1.664659,1.828188]    0.06125317 0.004834464 0.01234350 0.02299737 0.02403956
    #> (1.828188,1.991716]    0.37974416 0.010516294 0.02832213 0.05583668 0.05842604
    #> (1.991716,2.155245]    0.31215156 0.003873963 0.01195570 0.02734578 0.02857270
    #>                       std.resid.0 std.resid.1 std.resid.2 std.resid.3
    #> [0.3564295,0.5199582]   2.0220663  -0.6749812  -0.7294006  -0.4088547
    #> (0.5199582,0.6834869]   3.1177497   0.6533110  -1.7881373  -1.1503175
    #> (0.6834869,0.8470155]   2.9442258   3.8103971  -3.1908333  -2.2299488
    #> (0.8470155,1.010544]    0.3830354   1.0660662   3.0937290  -4.5167547
    #> (1.010544,1.174073]    -0.3169978   3.8565422  -0.3179482  -2.5699613
    #> (1.174073,1.337602]    -0.1568583   0.2811669   6.2600333  -6.3241558
    #> (1.337602,1.50113]     -1.4333937   4.1762987   4.2096551  -6.2340132
    #> (1.50113,1.664659]     -2.3587698  -4.3842081  -3.2591881   6.3675177
    #> (1.664659,1.828188]    -2.0409436  -5.5381760   0.7380923   2.5480159
    #> (1.828188,1.991716]    -0.7312476  -2.0763275  -5.6100775   6.4995706
    #> (1.991716,2.155245]    -1.0229057  -3.2723764  -9.8393733  10.9248190
    # (2) the raw residual plots using the object
    # 'fit1'
    plot(x = fit1, item.loc = 113, type = "icc", ci.method = "wald",
      layout.col = 2, ylim.sr.adjust = TRUE)

<img src="man/figures/README-example-6.png" width="70%" height="50%" />

    #>                           theta   N freq.0 freq.1 freq.2 freq.3 obs.prop.0
    #> [0.3564295,0.5199582] 0.3564295   1      1      0      0      0 1.00000000
    #> (0.5199582,0.6834869] 0.6081321   5      3      2      0      0 0.60000000
    #> (0.6834869,0.8470155] 0.7400138  15      5     10      0      0 0.33333333
    #> (0.8470155,1.010544]  0.8866202  55      5     15     34      1 0.09090909
    #> (1.010544,1.174073]   1.0821064 133      6     40     53     34 0.04511278
    #> (1.174073,1.337602]   1.2832293 260      8     37    153     62 0.03076923
    #> (1.337602,1.50113]    1.4747336  98      0     23     57     18 0.00000000
    #> (1.50113,1.664659]    1.5311735 306      0      7     85    214 0.00000000
    #> (1.664659,1.828188]   1.7632607 418      0      0    145    273 0.00000000
    #> (1.828188,1.991716]   1.8577191  69      0      0      0     69 0.00000000
    #> (1.991716,2.155245]   2.1021956 263      0      0      0    263 0.00000000
    #>                       obs.prop.1 obs.prop.2 obs.prop.3  exp.prob.0 exp.prob.1
    #> [0.3564295,0.5199582] 0.00000000  0.0000000 0.00000000 0.196511833 0.31299790
    #> (0.5199582,0.6834869] 0.40000000  0.0000000 0.00000000 0.130431315 0.27025065
    #> (0.6834869,0.8470155] 0.66666667  0.0000000 0.00000000 0.102631449 0.24407213
    #> (0.8470155,1.010544]  0.27272727  0.6181818 0.01818182 0.077129446 0.21379299
    #> (1.010544,1.174073]   0.30075188  0.3984962 0.25563910 0.051169395 0.17398130
    #> (1.174073,1.337602]   0.14230769  0.5884615 0.23846154 0.032494074 0.13632441
    #> (1.337602,1.50113]    0.23469388  0.5816327 0.18367347 0.020534959 0.10523862
    #> (1.50113,1.664659]    0.02287582  0.2777778 0.69934641 0.017857642 0.09707782
    #> (1.664659,1.828188]   0.00000000  0.3468900 0.65311005 0.009866868 0.06836048
    #> (1.828188,1.991716]   0.00000000  0.0000000 1.00000000 0.007690015 0.05880602
    #> (1.991716,2.155245]   0.00000000  0.0000000 1.00000000 0.003962699 0.03912356
    #>                       exp.prob.2 exp.prob.3  raw_resid.0  raw_resid.1
    #> [0.3564295,0.5199582]  0.3472692  0.1432210  0.803488167 -0.312997903
    #> (0.5199582,0.6834869]  0.3900531  0.2092649  0.469568685  0.129749354
    #> (0.6834869,0.8470155]  0.4043226  0.2489738  0.230701885  0.422594537
    #> (0.8470155,1.010544]   0.4127992  0.2962784  0.013779645  0.058934282
    #> (1.010544,1.174073]    0.4120662  0.3627831 -0.006056613  0.126770584
    #> (1.174073,1.337602]    0.3983962  0.4327853 -0.001724843  0.005983284
    #> (1.337602,1.50113]     0.3756891  0.4985373 -0.020534959  0.129455259
    #> (1.50113,1.664659]     0.3676106  0.5174539 -0.017857642 -0.074201998
    #> (1.664659,1.828188]    0.3299158  0.5918569 -0.009866868 -0.068360484
    #> (1.828188,1.991716]    0.3132481  0.6202558 -0.007690015 -0.058806016
    #> (1.991716,2.155245]    0.2690653  0.6878484 -0.003962699 -0.039123555
    #>                       raw_resid.2 raw_resid.3        se.0       se.1       se.2
    #> [0.3564295,0.5199582] -0.34726922 -0.14322104 0.397359953 0.46371351 0.47610220
    #> (0.5199582,0.6834869] -0.39005312 -0.20926492 0.150611412 0.19860274 0.21813376
    #> (0.6834869,0.8470155] -0.40432265 -0.24897378 0.078357401 0.11090564 0.12671381
    #> (0.8470155,1.010544]   0.20538261 -0.27809653 0.035974863 0.05528201 0.06638675
    #> (1.010544,1.174073]   -0.01356995 -0.10714402 0.019106171 0.03287157 0.04267975
    #> (1.174073,1.337602]    0.19006531 -0.19432375 0.010996190 0.02128019 0.03036171
    #> (1.337602,1.50113]     0.20594357 -0.31486387 0.014326112 0.03099761 0.04892172
    #> (1.50113,1.664659]    -0.08983281  0.18189246 0.007570744 0.01692484 0.02756294
    #> (1.664659,1.828188]    0.01697418  0.06125317 0.004834464 0.01234350 0.02299737
    #> (1.828188,1.991716]   -0.31324812  0.37974416 0.010516294 0.02832213 0.05583668
    #> (1.991716,2.155245]   -0.26906531  0.31215156 0.003873963 0.01195570 0.02734578
    #>                             se.3 std_resid.0 std_resid.1 std_resid.2
    #> [0.3564295,0.5199582] 0.35029812   2.0220663  -0.6749812  -0.7294006
    #> (0.5199582,0.6834869] 0.18191928   3.1177497   0.6533110  -1.7881373
    #> (0.6834869,0.8470155] 0.11165000   2.9442258   3.8103971  -3.1908333
    #> (0.8470155,1.010544]  0.06156999   0.3830354   1.0660662   3.0937290
    #> (1.010544,1.174073]   0.04169091  -0.3169978   3.8565422  -0.3179482
    #> (1.174073,1.337602]   0.03072722  -0.1568583   0.2811669   6.2600333
    #> (1.337602,1.50113]    0.05050741  -1.4333937   4.1762987   4.2096551
    #> (1.50113,1.664659]    0.02856568  -2.3587698  -4.3842081  -3.2591881
    #> (1.664659,1.828188]   0.02403956  -2.0409436  -5.5381760   0.7380923
    #> (1.828188,1.991716]   0.05842604  -0.7312476  -2.0763275  -5.6100775
    #> (1.991716,2.155245]   0.02857270  -1.0229057  -3.2723764  -9.8393733
    #>                       std_resid.3  raw.resid.0  raw.resid.1 raw.resid.2
    #> [0.3564295,0.5199582]  -0.4088547  0.803488167 -0.312997903 -0.34726922
    #> (0.5199582,0.6834869]  -1.1503175  0.469568685  0.129749354 -0.39005312
    #> (0.6834869,0.8470155]  -2.2299488  0.230701885  0.422594537 -0.40432265
    #> (0.8470155,1.010544]   -4.5167547  0.013779645  0.058934282  0.20538261
    #> (1.010544,1.174073]    -2.5699613 -0.006056613  0.126770584 -0.01356995
    #> (1.174073,1.337602]    -6.3241558 -0.001724843  0.005983284  0.19006531
    #> (1.337602,1.50113]     -6.2340132 -0.020534959  0.129455259  0.20594357
    #> (1.50113,1.664659]      6.3675177 -0.017857642 -0.074201998 -0.08983281
    #> (1.664659,1.828188]     2.5480159 -0.009866868 -0.068360484  0.01697418
    #> (1.828188,1.991716]     6.4995706 -0.007690015 -0.058806016 -0.31324812
    #> (1.991716,2.155245]    10.9248190 -0.003962699 -0.039123555 -0.26906531
    #>                       raw.resid.3      se.0.1     se.1.1     se.2.1     se.3.1
    #> [0.3564295,0.5199582] -0.14322104 0.397359953 0.46371351 0.47610220 0.35029812
    #> (0.5199582,0.6834869] -0.20926492 0.150611412 0.19860274 0.21813376 0.18191928
    #> (0.6834869,0.8470155] -0.24897378 0.078357401 0.11090564 0.12671381 0.11165000
    #> (0.8470155,1.010544]  -0.27809653 0.035974863 0.05528201 0.06638675 0.06156999
    #> (1.010544,1.174073]   -0.10714402 0.019106171 0.03287157 0.04267975 0.04169091
    #> (1.174073,1.337602]   -0.19432375 0.010996190 0.02128019 0.03036171 0.03072722
    #> (1.337602,1.50113]    -0.31486387 0.014326112 0.03099761 0.04892172 0.05050741
    #> (1.50113,1.664659]     0.18189246 0.007570744 0.01692484 0.02756294 0.02856568
    #> (1.664659,1.828188]    0.06125317 0.004834464 0.01234350 0.02299737 0.02403956
    #> (1.828188,1.991716]    0.37974416 0.010516294 0.02832213 0.05583668 0.05842604
    #> (1.991716,2.155245]    0.31215156 0.003873963 0.01195570 0.02734578 0.02857270
    #>                       std.resid.0 std.resid.1 std.resid.2 std.resid.3
    #> [0.3564295,0.5199582]   2.0220663  -0.6749812  -0.7294006  -0.4088547
    #> (0.5199582,0.6834869]   3.1177497   0.6533110  -1.7881373  -1.1503175
    #> (0.6834869,0.8470155]   2.9442258   3.8103971  -3.1908333  -2.2299488
    #> (0.8470155,1.010544]    0.3830354   1.0660662   3.0937290  -4.5167547
    #> (1.010544,1.174073]    -0.3169978   3.8565422  -0.3179482  -2.5699613
    #> (1.174073,1.337602]    -0.1568583   0.2811669   6.2600333  -6.3241558
    #> (1.337602,1.50113]     -1.4333937   4.1762987   4.2096551  -6.2340132
    #> (1.50113,1.664659]     -2.3587698  -4.3842081  -3.2591881   6.3675177
    #> (1.664659,1.828188]    -2.0409436  -5.5381760   0.7380923   2.5480159
    #> (1.828188,1.991716]    -0.7312476  -2.0763275  -5.6100775   6.4995706
    #> (1.991716,2.155245]    -1.0229057  -3.2723764  -9.8393733  10.9248190
    # (3) the standardized residual plots using the
    # object 'fit1'
    plot(x = fit1, item.loc = 113, type = "sr", ci.method = "wald",
      layout.col = 4, ylim.sr.adjust = TRUE)

<img src="man/figures/README-example-7.png" width="70%" height="50%" />

    #>                           theta   N freq.0 freq.1 freq.2 freq.3 obs.prop.0
    #> [0.3564295,0.5199582] 0.3564295   1      1      0      0      0 1.00000000
    #> (0.5199582,0.6834869] 0.6081321   5      3      2      0      0 0.60000000
    #> (0.6834869,0.8470155] 0.7400138  15      5     10      0      0 0.33333333
    #> (0.8470155,1.010544]  0.8866202  55      5     15     34      1 0.09090909
    #> (1.010544,1.174073]   1.0821064 133      6     40     53     34 0.04511278
    #> (1.174073,1.337602]   1.2832293 260      8     37    153     62 0.03076923
    #> (1.337602,1.50113]    1.4747336  98      0     23     57     18 0.00000000
    #> (1.50113,1.664659]    1.5311735 306      0      7     85    214 0.00000000
    #> (1.664659,1.828188]   1.7632607 418      0      0    145    273 0.00000000
    #> (1.828188,1.991716]   1.8577191  69      0      0      0     69 0.00000000
    #> (1.991716,2.155245]   2.1021956 263      0      0      0    263 0.00000000
    #>                       obs.prop.1 obs.prop.2 obs.prop.3  exp.prob.0 exp.prob.1
    #> [0.3564295,0.5199582] 0.00000000  0.0000000 0.00000000 0.196511833 0.31299790
    #> (0.5199582,0.6834869] 0.40000000  0.0000000 0.00000000 0.130431315 0.27025065
    #> (0.6834869,0.8470155] 0.66666667  0.0000000 0.00000000 0.102631449 0.24407213
    #> (0.8470155,1.010544]  0.27272727  0.6181818 0.01818182 0.077129446 0.21379299
    #> (1.010544,1.174073]   0.30075188  0.3984962 0.25563910 0.051169395 0.17398130
    #> (1.174073,1.337602]   0.14230769  0.5884615 0.23846154 0.032494074 0.13632441
    #> (1.337602,1.50113]    0.23469388  0.5816327 0.18367347 0.020534959 0.10523862
    #> (1.50113,1.664659]    0.02287582  0.2777778 0.69934641 0.017857642 0.09707782
    #> (1.664659,1.828188]   0.00000000  0.3468900 0.65311005 0.009866868 0.06836048
    #> (1.828188,1.991716]   0.00000000  0.0000000 1.00000000 0.007690015 0.05880602
    #> (1.991716,2.155245]   0.00000000  0.0000000 1.00000000 0.003962699 0.03912356
    #>                       exp.prob.2 exp.prob.3  raw_resid.0  raw_resid.1
    #> [0.3564295,0.5199582]  0.3472692  0.1432210  0.803488167 -0.312997903
    #> (0.5199582,0.6834869]  0.3900531  0.2092649  0.469568685  0.129749354
    #> (0.6834869,0.8470155]  0.4043226  0.2489738  0.230701885  0.422594537
    #> (0.8470155,1.010544]   0.4127992  0.2962784  0.013779645  0.058934282
    #> (1.010544,1.174073]    0.4120662  0.3627831 -0.006056613  0.126770584
    #> (1.174073,1.337602]    0.3983962  0.4327853 -0.001724843  0.005983284
    #> (1.337602,1.50113]     0.3756891  0.4985373 -0.020534959  0.129455259
    #> (1.50113,1.664659]     0.3676106  0.5174539 -0.017857642 -0.074201998
    #> (1.664659,1.828188]    0.3299158  0.5918569 -0.009866868 -0.068360484
    #> (1.828188,1.991716]    0.3132481  0.6202558 -0.007690015 -0.058806016
    #> (1.991716,2.155245]    0.2690653  0.6878484 -0.003962699 -0.039123555
    #>                       raw_resid.2 raw_resid.3        se.0       se.1       se.2
    #> [0.3564295,0.5199582] -0.34726922 -0.14322104 0.397359953 0.46371351 0.47610220
    #> (0.5199582,0.6834869] -0.39005312 -0.20926492 0.150611412 0.19860274 0.21813376
    #> (0.6834869,0.8470155] -0.40432265 -0.24897378 0.078357401 0.11090564 0.12671381
    #> (0.8470155,1.010544]   0.20538261 -0.27809653 0.035974863 0.05528201 0.06638675
    #> (1.010544,1.174073]   -0.01356995 -0.10714402 0.019106171 0.03287157 0.04267975
    #> (1.174073,1.337602]    0.19006531 -0.19432375 0.010996190 0.02128019 0.03036171
    #> (1.337602,1.50113]     0.20594357 -0.31486387 0.014326112 0.03099761 0.04892172
    #> (1.50113,1.664659]    -0.08983281  0.18189246 0.007570744 0.01692484 0.02756294
    #> (1.664659,1.828188]    0.01697418  0.06125317 0.004834464 0.01234350 0.02299737
    #> (1.828188,1.991716]   -0.31324812  0.37974416 0.010516294 0.02832213 0.05583668
    #> (1.991716,2.155245]   -0.26906531  0.31215156 0.003873963 0.01195570 0.02734578
    #>                             se.3 std_resid.0 std_resid.1 std_resid.2
    #> [0.3564295,0.5199582] 0.35029812   2.0220663  -0.6749812  -0.7294006
    #> (0.5199582,0.6834869] 0.18191928   3.1177497   0.6533110  -1.7881373
    #> (0.6834869,0.8470155] 0.11165000   2.9442258   3.8103971  -3.1908333
    #> (0.8470155,1.010544]  0.06156999   0.3830354   1.0660662   3.0937290
    #> (1.010544,1.174073]   0.04169091  -0.3169978   3.8565422  -0.3179482
    #> (1.174073,1.337602]   0.03072722  -0.1568583   0.2811669   6.2600333
    #> (1.337602,1.50113]    0.05050741  -1.4333937   4.1762987   4.2096551
    #> (1.50113,1.664659]    0.02856568  -2.3587698  -4.3842081  -3.2591881
    #> (1.664659,1.828188]   0.02403956  -2.0409436  -5.5381760   0.7380923
    #> (1.828188,1.991716]   0.05842604  -0.7312476  -2.0763275  -5.6100775
    #> (1.991716,2.155245]   0.02857270  -1.0229057  -3.2723764  -9.8393733
    #>                       std_resid.3  raw.resid.0  raw.resid.1 raw.resid.2
    #> [0.3564295,0.5199582]  -0.4088547  0.803488167 -0.312997903 -0.34726922
    #> (0.5199582,0.6834869]  -1.1503175  0.469568685  0.129749354 -0.39005312
    #> (0.6834869,0.8470155]  -2.2299488  0.230701885  0.422594537 -0.40432265
    #> (0.8470155,1.010544]   -4.5167547  0.013779645  0.058934282  0.20538261
    #> (1.010544,1.174073]    -2.5699613 -0.006056613  0.126770584 -0.01356995
    #> (1.174073,1.337602]    -6.3241558 -0.001724843  0.005983284  0.19006531
    #> (1.337602,1.50113]     -6.2340132 -0.020534959  0.129455259  0.20594357
    #> (1.50113,1.664659]      6.3675177 -0.017857642 -0.074201998 -0.08983281
    #> (1.664659,1.828188]     2.5480159 -0.009866868 -0.068360484  0.01697418
    #> (1.828188,1.991716]     6.4995706 -0.007690015 -0.058806016 -0.31324812
    #> (1.991716,2.155245]    10.9248190 -0.003962699 -0.039123555 -0.26906531
    #>                       raw.resid.3      se.0.1     se.1.1     se.2.1     se.3.1
    #> [0.3564295,0.5199582] -0.14322104 0.397359953 0.46371351 0.47610220 0.35029812
    #> (0.5199582,0.6834869] -0.20926492 0.150611412 0.19860274 0.21813376 0.18191928
    #> (0.6834869,0.8470155] -0.24897378 0.078357401 0.11090564 0.12671381 0.11165000
    #> (0.8470155,1.010544]  -0.27809653 0.035974863 0.05528201 0.06638675 0.06156999
    #> (1.010544,1.174073]   -0.10714402 0.019106171 0.03287157 0.04267975 0.04169091
    #> (1.174073,1.337602]   -0.19432375 0.010996190 0.02128019 0.03036171 0.03072722
    #> (1.337602,1.50113]    -0.31486387 0.014326112 0.03099761 0.04892172 0.05050741
    #> (1.50113,1.664659]     0.18189246 0.007570744 0.01692484 0.02756294 0.02856568
    #> (1.664659,1.828188]    0.06125317 0.004834464 0.01234350 0.02299737 0.02403956
    #> (1.828188,1.991716]    0.37974416 0.010516294 0.02832213 0.05583668 0.05842604
    #> (1.991716,2.155245]    0.31215156 0.003873963 0.01195570 0.02734578 0.02857270
    #>                       std.resid.0 std.resid.1 std.resid.2 std.resid.3
    #> [0.3564295,0.5199582]   2.0220663  -0.6749812  -0.7294006  -0.4088547
    #> (0.5199582,0.6834869]   3.1177497   0.6533110  -1.7881373  -1.1503175
    #> (0.6834869,0.8470155]   2.9442258   3.8103971  -3.1908333  -2.2299488
    #> (0.8470155,1.010544]    0.3830354   1.0660662   3.0937290  -4.5167547
    #> (1.010544,1.174073]    -0.3169978   3.8565422  -0.3179482  -2.5699613
    #> (1.174073,1.337602]    -0.1568583   0.2811669   6.2600333  -6.3241558
    #> (1.337602,1.50113]     -1.4333937   4.1762987   4.2096551  -6.2340132
    #> (1.50113,1.664659]     -2.3587698  -4.3842081  -3.2591881   6.3675177
    #> (1.664659,1.828188]    -2.0409436  -5.5381760   0.7380923   2.5480159
    #> (1.828188,1.991716]    -0.7312476  -2.0763275  -5.6100775   6.4995706
    #> (1.991716,2.155245]    -1.0229057  -3.2723764  -9.8393733  10.9248190
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_info.R
\name{test.info}
\alias{test.info}
\alias{test.info.default}
\alias{test.info.est_item}
\alias{test.info.est_irt}
\title{Item and Test Information Function}
\usage{
test.info(x, ...)

\method{test.info}{default}(x, theta, D = 1, ...)

\method{test.info}{est_item}(x, theta, ...)

\method{test.info}{est_irt}(x, theta, ...)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}. See below for details.}

\item{...}{Further arguments passed to or from other methods.}

\item{theta}{A vector of theta values where item and test information values are computed.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}
}
\value{
This function returns an object of class \code{\link{test.info}}. This object contains item and test information values
given the specified theta values.
}
\description{
This function computes both item and test information functions (Hambleton et al., 1991) given a set of theta values.
}
\details{
A specific form of a data frame should be used for the argument \code{x}. The first column should have item IDs,
the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
In the \pkg{irtplay} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as 
the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item. 
Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because 
the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories, 
four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
  ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
  ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
  ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
}
And an example of a data frame for a mixed-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
  ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
  ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
  ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
  ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
  ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
  ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
}
See \code{IRT Models} section in the page of \code{\link{irtplay-package}} for more details about the IRT models used in the \pkg{irtplay} package. 
An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Default method to compute item and test information functions for a data frame \code{x} containing the item metadata.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.

\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.
}}

\examples{
## example 1.
## using the function "shape_df" to create a data frame of test metadata
# create a list containing the dichotomous item parameters
par.dc <- list(a=c(1.1, 1.2, 0.9, 1.8, 1.4),
               b=c(0.1, -1.6, -0.2, 1.0, 1.2),
               g=rep(0.2, 5))

# create a list containing the polytomous item parameters
par.py <- list(a=c(1.4, 0.6),
               d=list(c(0.0, -1.9, 1.2), c(0.4, -1.1, 1.5, 0.2)))

# create a numeric vector of score categories for the items
cats <- c(2, 4, 2, 2, 5, 2, 2)

# create a character vector of IRT models for the items
model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")

# create an item metadata set
test <- shape_df(par.dc=par.dc, par.py=par.py, cats=cats, model=model) # create a data frame

# set theta values
theta <- seq(-2, 2, 0.1)

# compute item and test information values given the theta values
test.info(x=test, theta=theta, D=1)


## example 2.
## using a "-prm.txt" file obtained from a flexMIRT
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# set theta values
theta <- seq(-2, 2, 0.1)

# compute item and test information values given the theta values
test.info(x=test_flex, theta=theta, D=1)

}
\references{
Hambleton, R. K., & Swaminathan, H., & Rogers, H. J. (1991) \emph{Fundamentals of item response theory}.
Newbury Park, CA: Sage.
}
\seealso{
\code{\link{plot.test.info}}, \code{\link{shape_df}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_traceline.R
\name{plot.traceline}
\alias{plot.traceline}
\title{Plot ICC and TCC}
\usage{
\method{plot}{traceline}(
  x,
  item.loc = NULL,
  score.curve = FALSE,
  overlap = FALSE,
  layout.col = 2,
  xlab.text,
  ylab.text,
  main.text,
  lab.size = 15,
  main.size = 15,
  axis.size = 15,
  line.color,
  line.size = 1,
  strip.size = 12,
  ...
)
}
\arguments{
\item{x}{An object of class \code{\link{traceline}}.}

\item{item.loc}{A numeric value indicating that the \emph{n}th item (or the location of item) is plotted.
If NULL, the TCC based on a total test form is drawn. Default is NULL.}

\item{score.curve}{Logical value. If TRUE, item score curve (i.e., a weighted sum of item category probabilities over the item scores) is plotted
in a panel. Otherwise, ICCs for all score categories are plotted in separate panels. For a dichotomous item, the item score curve is the same as
the ICC of score category 1. Ignored when \code{item.loc = NULL}. Default is FALSE.}

\item{overlap}{Logical value indicating whether multiple item score curves are plotted in one panel. 
If FALSE, the multiple item score curves are displayed with multiple panels, one for each.}

\item{layout.col}{An integer value indicating the number of columns in the panel when displaying ICCs for an item or 
when displaying multiple item scores with multiple panels.}

\item{xlab.text, ylab.text}{A title for the x and y axes.}

\item{main.text}{An overall title for the plot.}

\item{lab.size}{The size of xlab and ylab. Default is 15.}

\item{main.size}{The size of \code{main.text}. Default is 15.}

\item{axis.size}{The size of labels along the x and y axes. Default is 15.}

\item{line.color}{A character string specifying the color for a line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/} for more details
about colors used in ggplot2.}

\item{line.size}{The size of lines. Default is 1.}

\item{strip.size}{The size of facet labels when ICCs for an item are plotted.}

\item{...}{Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.}
}
\description{
This function plots item or test characteristic curve using the ggplot2 package. The item characteristic
(or category) curve (ICC) or item score curve is drawn for an individual item. The test characteristic curve (TCC) is drawn
based on a total test form.
}
\details{
All of the plots are drawn using the ggplot2 package.
If \code{item.loc = NULL}, the TCC based on the total test form is plotted. In the argument \code{item.loc},
a vector of positive integer values should be specified to indicate the \emph{n}th items among the total test form. For example,
if there are ten items in the test form and the score curves of the 1st, 2nd, and 3rd items should be plotted, then \code{item.loc = 1:3}.
}
\examples{
## example
## using a "-prm.txt" file obtained from a flexMIRT
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# set theta values
theta <- seq(-3, 3, 0.1)

# compute the item category probabilities and item/test
# characteristic functions given the theta values
x <- traceline(x=test_flex, theta, D=1)

# plot TCC based on the total test form
plot(x, item.loc=NULL)

# plot ICCs for the first item (dichotomous item)
plot(x, item.loc=1, score.curve=FALSE, layout.col=2)

# plot item score curve for the first item (dichotomous item)
plot(x, item.loc=1, score.curve=TRUE)

# plot item score curves for the first six dichotomous items
# with multiple panels 
plot(x, item.loc=1:6, score.curve=TRUE, overlap=FALSE)

# plot item score curve for the first six dichotomous items
# in one panel
plot(x, item.loc=1:6, score.curve=TRUE, overlap=TRUE)

# plot ICCs for the last item (polytomous item)
plot(x, item.loc=55, score.curve=FALSE, layout.col=2)

# plot item score curve for the last item (polytomous item)
plot(x, item.loc=55, score.curve=TRUE)

# plot item score curves for the last three polytomous items
# with multiple panels
plot(x, item.loc=53:55, score.curve=TRUE, overlap=FALSE)

# plot item score curves for the last three poltyomous items
# in one panel
plot(x, item.loc=53:55, score.curve=TRUE, overlap=TRUE)

}
\seealso{
\code{\link{traceline}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimCAT_MX.R
\docType{data}
\name{simCAT_MX}
\alias{simCAT_MX}
\title{Simulated mixed-item format CAT Data}
\format{
This data includes a list of length three. The first internal object is a data.frame of
the item pool consisting of 200 dichotomous items and 30 polytomous items. The dichotomous items were
calibrated with the IRT 3PL model and the polytomous items were calibrated with the generalized
partial credit model. All polytomous items have three score categories (i.e., 0, 1, 2). The second
internal object is the response data set including a sparse response data set of 30,000 examinees
for the items in the item pool. The third internal object is the examinee's ability estimates
for 30,000 examinees.
}
\usage{
simCAT_MX
}
\description{
This data set contains an item pool information, response data, and examinee's ability estimates.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSAT6.R
\docType{data}
\name{LSAT6}
\alias{LSAT6}
\title{LSAT6 data}
\format{
This data contains 1,000 dichotomous response patterns of five items obtained from
the Law School Admissions Test, section 6.
}
\usage{
LSAT6
}
\description{
Well-known LSAT6 dichotomous response data set from Thissen (1982).
}
\references{
Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
\emph{Psychometrika, 47}, 175-186.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sx2_fit.R
\name{sx2_fit}
\alias{sx2_fit}
\alias{sx2_fit.default}
\alias{sx2_fit.est_item}
\alias{sx2_fit.est_irt}
\title{S-X2 fit statistic}
\usage{
sx2_fit(x, ...)

\method{sx2_fit}{default}(
  x,
  data,
  D = 1,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  ...
)

\method{sx2_fit}{est_item}(
  x,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  ...
)

\method{sx2_fit}{est_irt}(
  x,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  ...
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object
of class \code{\link{est_item}} obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}}
obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}}
for more details about the item metadata. The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}.}

\item{...}{Further arguments passed to or from other methods.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{alpha}{A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test for \eqn{S-X^{2}} fit statistic. Default is .05.}

\item{min.collapse}{An integer value to indicate the minimum frequency of cells to be collapsed. Default is 1. See below for details.}

\item{norm.prior}{A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
c(0,1).}

\item{nquad}{An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 30.}

\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
using the function \code{\link{gen.weight}}. If missing, default values are used (see the arguments of \code{norm.prior} and \code{nquad}).}
}
\value{
This function returns a list. Within a list, several internal objects are contained such as:
\item{fit_stat}{A data frame containing the results of \eqn{S-X^{2}} fit statistics for all items.}
\item{item_df}{The item metadata specified in the argument \code{x}.}
\item{exp_freq}{A list containing the collapsed expected frequency tables for all items.}
\item{obs_freq}{A list containing the collapsed observed frequency tables for all items.}
\item{exp_prob}{A list containing the collapsed expected probability tables for all items.}
\item{obs_prop}{A list containing the collapsed observed proportion tables for all items.}
}
\description{
This function computes \eqn{S-X^{2}} (Orlando & Thissen, 2000, 2003) item fit statistic.
}
\details{
Often, very small expected frequencies in the contingency tables used to compute \eqn{\chi^{2}} fit statistics could
compromise the accuracy of the \eqn{\chi^{2}} approximation for their distribution (Orlando & Thissen, 2000).
To avoid this problem, Orlando and Thissen (2000) used an algorithm of collapsing adjacent test score groups to maintain
a minimum expected category frequency of 1. However, if Orlando and Thissen's cell collapsing approach is applied to polytomous data,
too much information would be lost (Kang & Chen, 2008). Thus, Kang and Chen (2008) collapsed adjacent cells of item score categories
for a specific score group to ensure a minimum expected category frequency of 1. The same collapsing strategies were applied
in the function \code{\link{sx2_fit}}. If a minimum expected category frequency needs to be set to different number, you can specify
the minimum value in the argument \code{min.collapse}.

Note that if "DRM" is specified for an item in the item metadata set, the item is considered as "3PLM" to compute degree of freedom of
the \eqn{S-X^{2}} fit statistic.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Default method to compute \eqn{S-X^{2}} fit statistics for a data frame \code{x} containing the item metadata.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.

\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.
}}

\examples{
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the first twenty dichotomous items and last polytomous item
# assuming that the test consists of twenty-one items
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:20, 55), ]

# generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

\donttest{
# compute fit statistics
fit <- sx2_fit(x=x, data=data, nquad=30)

# fit statistics
fit$fit_stat
}

}
\references{
Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item fit index for polytomous IRT models.
\emph{Journal of Educational Measurement, 45}(4), 391-406.

Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for dichotomous item response theory models.
\emph{Applied Psychological Measurement, 24}(1), 50-64.

Orlando, M., & Thissen, D. (2003). Further investigation of the performance of S-X2: An item fit index for use with
dichotomous item response theory models. \emph{Applied Psychological Measurement, 27}(4), 289-298.
}
\seealso{
\code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtmodel.R
\name{plm}
\alias{plm}
\title{Polytomous Response Model Probabilities (GRM and GPCM)}
\usage{
plm(theta, a, d, D = 1, pmodel = c("GRM", "GPCM"))
}
\arguments{
\item{theta}{A vector of ability values.}

\item{a}{A numeric value of item discrimination (or slope) parameter.}

\item{d}{A vector of item difficulty (or threshold) parameters.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function  (if set to 1.7).
Default is 1.}

\item{pmodel}{A character string indicating the polytomous model being used. Available models are "GRM" for
the the graded response model and "GPCM" for the (generalized) partial credit model.}
}
\value{
This function returns a vector or matrix. When a matrix is returned, rows indicate theta values and columns represent
categories of an item.
}
\description{
This function computes the probability of selecting a specific category for an item
for a given set of theta values using the graded response model, partial credit model, and generalized
partial credit model.
}
\details{
When the category probabilities are computed for an item with the partial credit model, \code{a = 1} for that item.
When \code{pmodel = "GPCM"}, \code{d} should include the item difficulty (or threshold) parameters. In the \pkg{irtplay} package, 
the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as the item location (or overall difficulty) 
parameter subtracted by the threshold parameter for unique score categories of the item. Note that when an GPCM item has \emph{K} 
unique score categories, \emph{K-1} item difficulty parameters are necessary because the item difficulty parameter for the first category 
boundary is always 0. For example, if an GPCM item has five score categories, four item difficulty parameters should be specified. 
For more details about the parameterization of the (generalized) partial credit model, See \code{IRT Models} section 
in the page of \code{\link{irtplay-package}}.
}
\examples{
## Category probabilities for an item with four categories
## using a generalized partial credit model
plm(theta=c(-0.2, 0, 0.5), a=1.4, d=c(-0.2, 0, 0.5), D=1, pmodel='GPCM')

## Category probabilities for an item with five categories
## using a graded response model
plm(theta=c(-0.2, 0, 0.5), a=1.2, d=c(-0.4, -0.2, 0.4, 1.5), D=1, pmodel='GRM')

}
\seealso{
\code{\link{drm}}, \code{\link{irtfit}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtmodel.R
\name{drm}
\alias{drm}
\title{Dichotomous Response Model Probabilities}
\usage{
drm(theta, a, b, g = NULL, D = 1)
}
\arguments{
\item{theta}{A vector of ability values.}

\item{a}{A vector of item discrimination (or slope) parameters.}

\item{b}{A vector of item difficulty (or threshold) parameters.}

\item{g}{A vector of item guessing parameters.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}
}
\value{
This function returns a vector or matrix. When a matrix is returned, rows indicate theta values and columns represent items.
}
\description{
This function computes the probability of correct answers for one or more items for a given set of theta values
using the IRT 1PL, 2PL, and 3PL models.
}
\details{
\code{g} does not need to be specified when the response probabilities of the 1PL and 2PL models are computed.
}
\examples{
## when vectors are used for both theta values and item parameters (3PLM)
drm(c(-0.1, 0.0, 1.5), a=c(1, 2), b=c(0, 1), g=c(0.2, 0.1), D=1)

## when vectors are only used for item parameters (2PLM)
drm(0.0, a=c(1, 2), b=c(0, 1), D=1)

## when vectors are only used for theta values (3PLM)
drm(c(-0.1, 0.0, 1.5), a=1, b=1, g=0.2, D=1)

}
\seealso{
\code{\link{plm}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary}
\alias{summary}
\alias{summary.est_irt}
\alias{summary.est_item}
\title{Summary of item calibration}
\usage{
summary(object, ...)

\method{summary}{est_irt}(object, ...)

\method{summary}{est_item}(object, ...)
}
\arguments{
\item{object}{An object of class \code{\link{est_irt}} or \code{\link{est_item}}.}

\item{...}{Further arguments passed to or from other methods.}
}
\description{
This function summarizes the IRT calibration results of \code{\link{est_irt}}
or \code{\link{est_item}} object.
}
\section{Methods (by class)}{
\itemize{
\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.
}}

\examples{
\donttest{
# fit the 1PL model to LSAT6 data and constrain the slope parameters to be equal
fit.1pl <- est_irt(data=LSAT6, D=1, model="1PLM", cats=2, fix.a.1pl=FALSE)

# summary of the estimation
summary(fit.1pl)
}

}
\seealso{
\code{\link{est_irt}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{bind.fill}
\alias{bind.fill}
\title{Bind Fill}
\usage{
bind.fill(List, type = c("rbind", "cbind"))
}
\arguments{
\item{List}{A list containing different length of numeric vectors}

\item{type}{A character string specifying whether rbind is used or cbind is used.}
}
\value{
A matrix.
}
\description{
This function creates a cbind matrix or rbind matrix using a list containing different length
of numeric vectors.
}
\examples{
# sample list
score_list <- list(item1=c(0:3), item2=c(0:2), item3=c(0:5), item3=c(0:4))

# examples
# 1) create a rbind with the sample score list
bind.fill(score_list, type="rbind")

# 2) create a cbind with the sample score list
bind.fill(score_list, type="cbind")

}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdif.R
\name{rdif}
\alias{rdif}
\title{IRT residual-based DIF (RDIF) detection framework}
\usage{
rdif(
  x,
  data,
  score = NULL,
  group,
  focal.name,
  D = 1,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("rdif_rs", "rdif_r", "rdif_s"),
  max.iter = 10,
  min.resp = NULL,
  method = "MLE",
  range = c(-4, 4),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}. See \code{\link{est_irt}}, \code{\link{irtfit}}, 
\code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{score}{A vector of examinees' ability estimates. If the abilities are not provided, \code{\link{rdif}} function estimates the abilities before 
computing RDIF statistics. See \code{\link{est_score}} for more details about scoring methods. Default is NULL.}

\item{group}{A numeric or character vector indicating group membership of examinees. The length of vector should the same with the number of rows 
in the response data matrix.}

\item{focal.name}{A single numeric or character indicating the level of group which corresponds to the focal group. 
For example, if \code{group = c(0, 1, 0, 1, 1)} and '1' indicates the focal group, then \code{focal.name = 1}.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{alpha}{A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test using the RDIF fit statistics.
Default is .05.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}

\item{purify}{A logical value indicating whether a purification process will be implemented or not. Default is FALSE.}

\item{purify.by}{A character string specifying a RDIF statistic with which the purification is implemented. Available statistics 
are "rdif_rs" for \eqn{RDIF_{RS}}, "rdif_r" for \eqn{RDIF_{R}}, and "rdif_s" for \eqn{RDIF_{S}}.}

\item{max.iter}{An positive integer value specifying the maximum number of iterations for the purification process. Default is 10.}

\item{min.resp}{An positive integer value specifying the minimum number of item responses for an examinee when scoring is conducted. 
Default is NULL. See details below for more information.}

\item{method}{A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation, 
"MAP" for the maximum a posteriori estimation, and "EAP" for the expected a posteriori estimation. Default method is "MLE".}

\item{range}{A numeric vector of two components to restrict the range of ability scale for the MLE. Default is c(-4, 4).}

\item{norm.prior}{A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
c(0,1). Ignored if \code{method} is "MLE".}

\item{nquad}{An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
Ignored if \code{method} is "MLE" or "MAP".}

\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP", default values are used (see the arguments
of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "MLE" or "MAP".}

\item{ncore}{The number of logical CPU cores to use. Default is 1. See \code{\link{est_score}} for details.}

\item{verbose}{A logical value. If TRUE, the progress messages of purification procedure are suppressed. Default is TRUE.}
}
\value{
This function returns a list of four internal objects. The four objects are: 
\item{no_purify}{A list of several sub-objects containing the results of DIF analysis without a purification procedure. The sub-objects are: 
    \describe{
      \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column 
       indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized, \eqn{RDIF_{S}}, 
       \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of 
       the focal group, sample size of the reference group, and total sample size, respectively. Note that \eqn{RDIF_{RS}} does not have its standardized 
       value because it is a \eqn{\chi^{2}} statistic.} 
      \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID, 
       mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, and
       covariance of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, respectively.}
      \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector 
       means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
      \item{score}{A vector of ability estimates used to compute the RDIF statistics.}
   }
}
\item{purify}{A logical value indicating whether the purification process was used.} 
\item{with_purify}{A list of several sub-objects containing the results of DIF analysis with a purification procedure. The sub-objects are:
    \describe{
      \item{purify.by}{A character string indicating which RDIF statistic is used for the purification. "rdif_r", "rdif_s", and "rdif_rs" refers to
       \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
      \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column 
       indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized, \eqn{RDIF_{S}}, 
       \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of 
       the focal group, sample size of the reference group, total sample size, and \emph{n}th iteration where the RDIF statistics were computed, 
       respectively.}
      \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID, 
       mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, covariance 
       of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, and \emph{n}th iteration where the RDIF statistics were computed, respectively.}
      \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector 
       means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
      \item{n.iter}{A total number of iterations repleated for the purification.}
      \item{score}{A vector of final purified ability estimates used to compute the RDIF statistics.}
      \item{complete}{A logical value indicating whether the purification process was completed. If FALSE, it means that the purification process 
       reached the maximum iteration number but it was not complete.}
    }
}
\item{alpha}{A significance \eqn{\alpha}-level used to compute the p-values of RDIF statistics.}
}
\description{
This function computes three RDIF statistics (Lim, Choe, & Han, Under review; Lim, Choe, Han, Lee, & Hong, 2021), 
which are \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, for each item. \eqn{RDIF_{R}} primarily 
captures the typical contrast in raw residual pattern between two groups caused by uniform DIF whereas 
\eqn{RDIF_{S}} primarily captures the typical contrast in squared residual pattern between two groups caused
by nonuniform DIF. \eqn{RDIF_{RS}} can reasonably capture both types of DIF.
}
\details{
The RDIF framework (Lim et al., Under review; Lim et al., 2021) consists of three IRT residual-based statistics: \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, 
and \eqn{RDIF_{RS}}. Under the null hypothesis that a test contains no DIF items, \eqn{RDIF_{R}} and \eqn{RDIF_{S}} follow 
normal distributions asymptotically. \eqn{RDIF_{RS}} is a based on a bivariate normal distribution of \eqn{RDIF_{R}} and 
\eqn{RDIF_{S}} statistics. Under the null hypothesis of no DIF items, it follows a \eqn{\chi^{2}} distribution asymptotically 
with 2 degrees of freedom. See Lim et al. (2021) for more details about RDIF framework. 

The \code{\link{rdif}} function computes all three RDIF statistics of \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}. The current 
version of \code{\link{rdif}} function only supports dichotomous item response data. To compute the three statistics, the \code{\link{rdif}} function 
requires (1) item parameter estimates obtained from aggregate data regardless of group membership, (2) examinees' ability estimates 
(e.g., MLE), and (3) examinees' item response data. Note that the ability estimates need to be computed using the aggregate data-based 
item parameter estimates. The item parameter estimates should be provided in the \code{x} argument, the ability estimates should 
be provided in the \code{score} argument, and the response data should be provided in the \code{data} argument. When the abilities 
are not given in the \code{score} argument (i.e., \code{score = NULL}), the \code{\link{rdif}} function estimates examinees' abilities 
automatically using the scoring method specified in the \code{method} argument (e.g., \code{method = "MLE"}). 

The \code{group} argument accepts a vector of either two distinct numeric or character variables. Between two distinct variable, one is to 
represent the reference group and another one is to represent the focal group. The length of the vector should be the same with the number
of rows in the response data and each value in the vector should indicate each examinee of the response data. Once the \code{gruop} is 
specified, a single numeric or character value needs to be provided in the \code{focal.name} argument to define which group variable in 
the \code{group} argument represents the focal group.   

As other DIF detection approaches, an iterative purification process can be implemented for the RDIF framework. 
When \code{purify = TRUE}, the purification process is implemented based on one of RDIF statistics specified in the \code{purify.by} 
argument (e.g, \code{purify.by="rdif_rs"}). At each iterative purification, examinees' latent abilities are computed using purified items and 
scoring method specified in the \code{method} argument. The iterative purification process stops when no further DIF items are found or 
the process reaches a predetermined limit of iteration, which can be specified in the \code{max.iter} argument. See Lim et al. (2021) 
for more details about the purification procedure.

Scoring with a few items entails large standard errors which in turn could compromise DIF detection with RDIF framework. 
The \code{min.resp} argument can be used to avoid using scores with large standard errors when computing the RDIF statistics, espeically 
during the purification process. For example, if \code{min.resp} is not NULL (e.g., \code{min.resp=5}), item responses of examinees 
whose tally of item responses are less than the specified minimum number are treated as missing values (i.e., NA). Accordingly, 
their ability estimates become missing values and are not used for computing the RDIF statistics. If \code{min.resp=NULL}, 
an examinee's score will be computed as long as there exists, at least, 1 item response for the examinee.
}
\examples{
\donttest{
# call library
library("dplyr")

## Uniform DIF detection 
###############################################
# (1) manipulate true uniform DIF data
###############################################
# import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select 36 of 3PLM items which are non-DIF items 
par_nstd <- 
  bring.flexmirt(file=flex_sam, "par")$Group1$full_df \%>\% 
  dplyr::filter(.data$model == "3PLM") \%>\% 
  dplyr::filter(dplyr::row_number() \%in\% 1:36) \%>\% 
  dplyr::select(1:6)
par_nstd$id <- paste0("nondif", 1:36)

# generate four new items to inject uniform DIF
difpar_ref <- 
  shape_df(par.dc=list(a=c(0.8, 1.5, 0.8, 1.5), b=c(0.0, 0.0, -0.5, -0.5), g=0.15), 
           item.id=paste0("dif", 1:4), cats=2, model="3PLM")

# manipulate uniform DIF on the four new items by adding constants to b-parameters 
# for the focal group
difpar_foc <- 
  difpar_ref \%>\% 
  dplyr::mutate_at(.vars="par.2", .funs=function(x) x + rep(0.7, 4))

# combine the 4 DIF and 36 non-DIF items for both reference and focal groups
# thus, the first four items have uniform DIF 
par_ref <- rbind(difpar_ref, par_nstd)
par_foc <- rbind(difpar_foc, par_nstd)

# generate the true thetas
set.seed(123)
theta_ref <- rnorm(500, 0.0, 1.0)
theta_foc <- rnorm(500, 0.0, 1.0)

# generate the response data
resp_ref <- simdat(par_ref, theta=theta_ref, D=1)
resp_foc <- simdat(par_foc, theta=theta_foc, D=1)
data <- rbind(resp_ref, resp_foc)

###############################################
# (2) estimate the item and ability parameters 
#     using the aggregate data
###############################################
# estimate the item parameters 
est_mod <- est_irt(data=data, D=1, model="3PLM")
est_par <- est_mod$par.est

# estimate the ability parameters using MLE
score <- est_score(x=est_par, data=data, method="MLE")$est.theta

###############################################
# (3) conduct DIF analysis
###############################################
# create a vector of group membership indicators 
# where '1' indicates the focal group 
group <- c(rep(0, 500), rep(1, 500))

# (a)-1 compute RDIF statistics by providing scores, 
#       and without a purification 
dif_nopuri_1 <- rdif(x=est_par, data=data, score=score, 
                     group=group, focal.name=1, D=1, alpha=0.05)
print(dif_nopuri_1)

# (a)-2 compute RDIF statistics by not providing scores 
#       and without a purification 
dif_nopuri_2 <- rdif(x=est_par, data=data, score=NULL, 
                     group=group, focal.name=1, D=1, alpha=0.05, 
                     method="MLE")
print(dif_nopuri_2)

# (b)-1 compute RDIF statistics with a purification 
#       based on \eqn{RDIF_{R}}
dif_puri_r <- rdif(x=est_par, data=data, score=score, 
                   group=group, focal.name=1, D=1, alpha=0.05, 
                   purify=TRUE, purify.by="rdif_r")
print(dif_puri_r)

# (b)-2 compute RDIF statistics with a purification 
#       based on \eqn{RDIF_{S}}
dif_puri_s <- rdif(x=est_par, data=data, score=score, 
                   group=group, focal.name=1, D=1, alpha=0.05, 
                   purify=TRUE, purify.by="rdif_s")
print(dif_puri_s)

# (b)-3 compute RDIF statistics with a purification 
#       based on \eqn{RDIF_{RS}}
dif_puri_rs <- rdif(x=est_par, data=data, score=score, 
                    group=group, focal.name=1, D=1, alpha=0.05, 
                    purify=TRUE, purify.by="rdif_rs")
print(dif_puri_rs)
}

}
\references{
Lim, H., Choe, E. M., & Han, K. T. (Under review). A residual-based differential item functioning detection framework in 
item response theory. \emph{Journal of Educational Measurement}. 

Lim, H., Choe, E. M., Han, K. T., Lee, S., & Hong, M. (2021, June). \emph{IRT residual approach 
to detecting DIF.} Paper presented at the Annual Meeting of the National Council on Measurement 
in Education. Online.
}
\seealso{
\code{\link{est_item}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, 
\code{\link{gen.weight}}, \code{\link{est_score}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_score.R
\name{est_score}
\alias{est_score}
\alias{est_score.default}
\alias{est_score.est_irt}
\title{Estimate examinees' ability (proficiency) parameters}
\usage{
est_score(x, ...)

\method{est_score}{default}(
  x,
  data,
  D = 1,
  method = "MLE",
  range = c(-4, 4),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  fence.a = 3,
  fence.b = NULL,
  se = TRUE,
  obs.info = TRUE,
  constant = 0.1,
  constraint = FALSE,
  range.tcc = c(-7, 7),
  missing = NA,
  ncore = 1,
  ...
)

\method{est_score}{est_irt}(
  x,
  method = "MLE",
  range = c(-4, 4),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  fence.a = 3,
  fence.b = NULL,
  se = TRUE,
  obs.info = TRUE,
  constant = 0.1,
  constraint = FALSE,
  range.tcc = c(-7, 7),
  missing = NA,
  ncore = 1,
  ...
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...) or an object of
class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{test.info}},
or \code{\link{simdat}} for more details about the item metadata. This data frame can be easily obtained using the function \code{\link{shape_df}}.}

\item{...}{additional arguments to pass to \code{parallel::makeCluster}.}

\item{data}{A matrix or vector containing examinees' response data for the items in the argument \code{x}. When a matrix is used, a row and column indicate
the examinees and items, respectively. When a vector is used, it should contains the item response data for an examinee.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{method}{A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation,
"MLEF" for the maximum likelihood estimation with fences, "MAP" for the maximum a posteriori estimation,
"EAP" for the expected a posteriori estimation, "EAP.SUM" for the expected a posteriori summed scoring, and "INV.TCC" for the inverse TCC scoring.
Default method is "MLE".}

\item{range}{A numeric vector of two components to restrict the range of ability scale for the MLE. Default is c(-4, 4).}

\item{norm.prior}{A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
c(0,1). Ignored if \code{method} is "MLE", "MLEF", or "INV.TCC".}

\item{nquad}{An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
Ignored if \code{method} is "MLE", "MLEF", "MAP", or "INV.TCC".}

\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP" or "EAP.SUM", default values are used (see the arguments
of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "MLE", "MLEF", "MAP", or "INV.TCC".}

\item{fence.a}{A numeric value specifying the item slope parameter (i.e., \emph{a}-parameter) for the two imaginary items in MLEF. See below for details.
Default is 3.0.}

\item{fence.b}{A numeric vector of two components specifying the lower and upper fences of item difficulty parameters (i.e., \emph{b}-parameters)
for the two imaginary items, respectively, in MLEF. When \code{fence.b = NULL}, the lower and upper fences of item difficulty parameters were
automatically set. See below for details. Default is NULL.}

\item{se}{A logical value. If TRUE, the standard errors of ability estimates are computed. However, if \code{method} is "EAP.SUM" or "INV.TCC", the standard
errors are always returned. Default is TRUE.}

\item{obs.info}{A logical value. If TRUE, the observed item information functions are used to compute the standard errors of ability estimates when "MLE", "MLEF", 
or "MAP" is specified in \code{method}. If FALSE, the expected item information (a.k.a. Fisher information) functions are used to compute the standard errors. 
Note that under the 1PL and 2PL models, the observed item information function is exactly equal to the expected item information function. Default is TRUE.}

\item{constant}{A numeric value used to adjust zero and perfect raw sum scores, or the raw sum score equal to the sum of item guessing parameters,
if necessary, to find estimable solutions for those raw sum scores when \code{method = "INV.TCC"}. The zero raw score is forced to become the score of "zero raw score + constant"
and the perfect raw score is forced to become the score of "perfect raw score - constant". If the 3PLM items are included in the item metadata,
the raw sum score equal to the sum of item guessing parameters is forced to become the score of "the raw sum score + constant". Default is .1.}

\item{constraint}{A logical value indicating whether the ability estimates will be restricted within a specific ability range
specified in the argument \code{range.tcc} when \code{method = "INV.TCC"}. If \code{constraint = TRUE}, all ability estimates less than the first value in the vector specified in
the argument \code{range.tcc} are transformed to the first value and all ability estimates greater than the second value in the vector specified in
the argument \code{range.tcc} are transformed to the second value. Also, when \code{constraint = TRUE} and the 3PLM items are contained
in the item metadata, linear interpolation method is used to find the ability estimates for the raw sum scores less than the sum of item guessing
parameters. When \code{constraint = FALSE} and the 3PLM items are contained in the item metadata, linear extrapolation method is used to find
the ability estimates for the raw sum scores less than the sum of item guessing parameters. See below for details. Default is FALSE.}

\item{range.tcc}{A numeric vector of two components to be used as the lower and upper bounds of ability estimates when \code{method = "INV.TCC"} and
\code{constraint = TRUE}. Default is c(-7, 7).}

\item{missing}{A value indicating missing values in the response data set. Default is NA. See below for details.}

\item{ncore}{The number of logical CPU cores to use. Default is 1. See below for details.}
}
\value{
A list including a vector of the ability estimates and a vector of the standard errors of ability estimates. When \code{method} is
"EAP.SUM" or "INV.TCC", raw sum scores of examinees and a table with the possible raw sum scores and corresponding ability estimates are returned as well.
}
\description{
This function estimates examinees' latent ability parameters. Available scoring methods are maximum likelihood estimation (MLE),
maximum likelihood estimation with fences (MLEF; Han, 2016), maximum a posteriori estimation (MAP; Hambleton et al., 1991),
expected a posteriori estimation (EAP; Bock & Mislevy, 1982), EAP summed scoring (Thissen et al., 1995; Thissen & Orlando, 2001),
and inverse test characteristic curve (TCC) scoring (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).
}
\details{
For MAP scoring method, only the normal prior distribution is available for the population distribution.

When there are missing data in the response data set, the missing value must be specified in \code{missing}. The missing data are taken into account
when either of MLE, MLEF, MAP, and EAP is used. However, there must be no missing data in the response data set when "EAP.SUM" or "INV.TCC" is used.
One of possible ways to use "EAP.SUM" or "INV.TCC" method when missing values exist is to remove rows with any missing values.

In the maximum likelihood estimation with fences (MLEF; Han, 2016), two 2PLM imaginary items are necessary. The first imaginary item serves as the lower
fence and its difficulty parameter (i.e., \emph{b}-parameters) should be lower than any difficulty parameter values in the test form. Likewise, the second
imaginary item serves as the upper fence and its difficulty parameter should be greater than any difficulty parameter values in the test form. Also, the two
imaginary items should have a very high item slope parameter (i.e., \emph{a}-parameter) value. See Han (2016) for more details.

When \code{fence.b = NULL} in MLEF, the function automatically sets the lower and upper fences of item difficulty parameters using two steps. More specifically,
in the first step, the lower fence of the item difficulty parameter is set to the greatest integer value less than the minimum of item difficulty parameters
in the item metadata and the upper fence of the item difficulty parameter is set to the smallest integer value greater than the maximum of item difficulty
parameters in the item metadata. Then, in the second step, if the lower fence set in the first step is greater than -3.5, the lower fence is constrained to -3.5
and if the upper fence set in the first step is less than 3.5, the upper fence is constrained to 3.5. Otherwise, the fence values of item difficulty parameters
set in the first step are used.

When "INV.TCC" method is used employing the IRT 3-parameter logistic model (3PLM) in a test, ability estimates for the raw sum scores less than the sum of item
guessing parameters are not attainable. In this case, either of linear interpolation and linear extrapolation can be applied. Note that
if \code{constraint = TRUE}, linear interpolation method is used. Otherwise, linear extrapolation method is used. Let \eqn{\theta_{min}} and
\eqn{\theta_{max}} be the minimum and maximum ability estimates and \eqn{\theta_{X}} be the ability estimate for the smallest raw score, X, greater than or equal
to the sum of item guessing parameters. When linear interpolation method is used, a linear line is constructed between two points of (x=\eqn{\theta_{min}}, y=0) and
(x=\eqn{\theta_{X}}, y=X). Because \code{constraint = TRUE}, \eqn{\theta_{min}} is the first value in the argument \code{range.tcc}.
When linear extrapolation method is used, a linear line is constructed using two points of (x=\eqn{\theta_{X}}, y=X) and
(x=\eqn{\theta_{max}}, y=maximum raw score). Then, ability estimates for the raw sum scores between zero and the smallest raw score greater than or equal
to the sum of item guessing parameters are found using the constructed linear line. When it comes to the scoring method of "INV.TCC", the standard errors of ability
estimates are computed using an approach suggested by Lim, Davey, and Wells (2020).

To speed up the ability estimation for MLE, MLEF, MAP, and EAP methods, this function applies a parallel process using multiple logical CPU cores.
You can set the number of logical CPU cores by specifying a positive integer value in the argument \code{ncore}. Default value is 1.

Note that the standard errors of ability estimates are computed using observed information functions for MLE, MLEF, and MAP methods.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Default method to estimate examinees' latent ability parameters using a data frame \code{x} containing the item metadata.

\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.
}}

\examples{
## the use of a "-prm.txt" file obtained from a flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# generate examinees abilities
set.seed(12)
theta <- rnorm(10)

# simulate the item response data
data <- simdat(x, theta, D=1)

\donttest{
# estimate the abilities using MLE
est_score(x, data, D=1, method="MLE", range=c(-4, 4), se=TRUE, ncore=2)

# estimate the abilities using MLEF with default fences of item difficulty parameters
est_score(x, data, D=1, method="MLEF", fence.a=3.0, fence.b=NULL, se=TRUE, ncore=2)

# estimate the abilities using MLEF with different fences of item difficulty parameters
est_score(x, data, D=1, method="MLEF", fence.a=3.0, fence.b=c(-5, 5), se=TRUE, ncore=2)

# estimate the abilities using MAP
est_score(x, data, D=1, method="MAP", norm.prior=c(0, 1), nquad=30, se=TRUE, ncore=2)

# estimate the abilities using EAP
est_score(x, data, D=1, method="EAP", norm.prior=c(0, 1), nquad=30, se=TRUE, ncore=2)

# estimate the abilities using EAP summed scoring
est_score(x, data, D=1, method="EAP.SUM", norm.prior=c(0, 1), nquad=30)

# estimate the abilities using inverse TCC scoring
est_score(x, data, D=1, method="INV.TCC", constant=0.1, constraint=TRUE, range.tcc=c(-7, 7))

}

}
\references{
Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Psychometrika, 35}, 179-198.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.

Han, K. T. (2016). Maximum likelihood score estimation method with fences for short-length tests and computerized adaptive tests.
\emph{Applied psychological measurement, 40}(4), 289-301.

Kolen, M. J. & Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
Springer

Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
\emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.

Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical approach to evaluate the performance of MST.
\emph{Journal of Educational Measurement}. DOI: 10.1111/jedm.12276.

Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
\emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.

Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
\emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
Measurement, 19}(1), 39-49.
}
\seealso{
\code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_item.R
\name{est_item}
\alias{est_item}
\title{Fixed ability parameter calibration}
\usage{
est_item(
  x = NULL,
  data,
  score,
  D = 1,
  model = NULL,
  cats = NULL,
  item.id = NULL,
  fix.a.1pl = FALSE,
  fix.a.gpcm = FALSE,
  fix.g = FALSE,
  a.val.1pl = 1,
  a.val.gpcm = 1,
  g.val = 0.2,
  use.aprior = FALSE,
  use.bprior = FALSE,
  use.gprior = TRUE,
  aprior = list(dist = "lnorm", params = c(0, 0.5)),
  bprior = list(dist = "norm", params = c(0, 1)),
  gprior = list(dist = "beta", params = c(5, 17)),
  missing = NA,
  use.startval = FALSE,
  control = list(eval.max = 500, iter.max = 500),
  verbose = TRUE
)
}
\arguments{
\item{x}{A data frame containing the item metadata. This metadata is necessary to obtain the information of
each item (i.e., number of score categories and IRT model) to be calibrated. You can easily create an empty
item metadata using the function \code{\link{shape_df}}. When \code{use.startval = TRUE}, the item parameters
specified in the item metadata are used as the starting values for the item parameter estimation.
If \code{x = NULL}, the arguments of \code{model} and \code{cats} must be specified. See \code{\link{irtfit}},
\code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{score}{A vector of examinees' ability estimates. Length of the vector must be the same as the number of rows in the
response data set.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{model}{A vector of character strings indicating what IRT model is used to calibrate each item. Available IRT models are
"1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function.
If a single character of the IRT model is specified, that model will be recycled across all items. This information is only required
when \code{x = NULL}.}

\item{cats}{A numeric vector specifying the number of score categories for each item. For example, a dichotomous
item has two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items
are binary items (i.e., dichotomous items), it assumes that all items have two score categories. This information is only required
when \code{x = NULL}.}

\item{item.id}{A character vector of item IDs. If NULL, the item IDs are generated automatically. Default is NULL.}

\item{fix.a.1pl}{A logical value. If TRUE, the slope parameters of the 1PLM items are fixed to a specific value specified in the argument
\code{a.val.1pl}. Otherwise, the slope parameters of all 1PLM items are constrained to be equal and estimated. Default is FALSE.}

\item{fix.a.gpcm}{A logical value. If TRUE, the GPCM items are calibrated with the partial credit model and the slope parameters of
the GPCM items are fixed to a specific value specified in the argument \code{a.val.gpcm}. Otherwise, the slope parameter of each GPCM item
is estimated. Default is FALSE.}

\item{fix.g}{A logical value. If TRUE, the guessing parameters of the 3PLM items are fixed to a specific value specified in the argument
\code{g.val}. Otherwise, the guessing parameter of each 3PLM item is estimated. Default is FALSE.}

\item{a.val.1pl}{A numeric value. This value is used to fixed the slope parameters of the 1PLM items.}

\item{a.val.gpcm}{A numeric value. This value is used to fixed the slope parameters of the GPCM items.}

\item{g.val}{A numeric value. This value is used to fixed the guessing parameters of the 3PLM items.}

\item{use.aprior}{A logical value. If TRUE, a prior distribution for the slope parameters is used for the parameter calibration
across all items. Default is FALSE.}

\item{use.bprior}{A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used for the parameter calibration
across all items. Default is FALSE.}

\item{use.gprior}{A logical value. If TRUE, a prior distribution for the guessing parameters is used for the parameter calibration
across all 3PLM items. Default is TRUE.}

\item{aprior}{A list containing the information of the prior distribution for item slope parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{bprior}{A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{gprior}{A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}

\item{use.startval}{A logical value. If TRUE, the item parameters provided in the item metadata (i.e., the argument \code{x}) are used as
the starting values for the item parameter estimation. Otherwise, internal starting values of this function are used. Default is FALSE.}

\item{control}{A list of control parameters to be passed to the optimization function of \code{nlminb()} in the \pkg{stats} package. The control parameters
set the conditions of the item parameter estimation process such as the maximum number of iterations. See \code{nlminb()} in the \pkg{stats} package for details.}

\item{verbose}{A logical value. If FALSE, all progress messages are suppressed. Default is TRUE.}
}
\value{
This function returns an object of class \code{\link{est_item}}. Within this object, several internal objects are contained such as:
\item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
\item{par.est}{A data frame containing the item parameter estimates.}
\item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
observed information functions.}
\item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
when interpreting the variance-covariance matrix of item parameter estimates.}
\item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.} 
\item{loglikelihood}{A sum of the log-likelihood values of the complete data set across all estimated items.}
\item{data}{A data frame of the examinees' response data set.}
\item{score}{A vector of the examinees' ability values used as the fixed effects.}
\item{scale.D}{A scaling factor in IRT models.}
\item{convergence}{A string indicating the convergence status of the item parameter estimation.}
\item{nitem}{A total number of items included in the response data.}
\item{deleted.item}{The items which have no item response data. Those items are excluded from the item parameter estimation.}
\item{npar.est}{A total number of the estimated parameters.}
\item{n.response}{An integer vector indicating the number of item responses for each item used to estimate the item parameters.}
\item{TotalTime}{Time (in seconds) spent for total compuatation.}

The internal objects can be easily extracted using the function \code{\link{getirt}}.
}
\description{
This function performs the fixed ability parameter calibration (FAPC), often called
Method A, which is the maximum likelihood estimation of item parameters given the ability
estimates (Baker & Kim, 2004; Ban, Hanson, Wang, Yi, & Harris, 2001; Stocking, 1988). Also, this could be 
considered as a special type of the joint maximum likelihood estimation where only one cycle of parameter 
estimation is implemented given the ability estimates (Birnbaum, 1968). FAPC is one of potentially useful 
online item calibration methods for computerized adaptive testing (CAT) to put the parameter estimates of 
pretest items on the same scale of operational item parameter estimates and recalibrate the operational 
items to evaluate the parameter drifts of the operational items (Chen & Wang, 2016; Stocking, 1988).
}
\details{
In most cases, the function \code{\link{est_item}} will return successfully converged item parameter estimates using
the default internal starting values. However, if there is a convergence problem in the calibration, one possible solution is using
different starting values. When the item parameter values are specified in the item metadata (i.e., the argument \code{x}), those values
can be used as the starting values for the item parameter calibration by setting \code{use.startval = TRUE}.
}
\examples{
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# modify the item metadata so that some items follow 1PLM, 2PLM and GPCM
x[c(1:3, 5), 3] <- "1PLM"
x[c(1:3, 5), 4] <- 1
x[c(1:3, 5), 6] <- 0
x[c(4, 8:12), 3] <- "2PLM"
x[c(4, 8:12), 6] <- 0
x[54:55, 3] <- "GPCM"

# generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

\donttest{
# 1) item parameter estimation: constrain the slope parameters of the 1PLM to be equal
(mod1 <- est_item(x, data, score, D=1, fix.a.1pl=FALSE, use.gprior=TRUE,
                  gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
summary(mod1)

# extract the item parameter estimates
getirt(mod1, what="par.est")

# 2) item parameter estimation: fix the slope parameters of the 1PLM to 1
(mod2 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, a.val.1pl=1, use.gprior=TRUE,
                  gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
summary(mod2)

# extract the standard error estimates
getirt(mod2, what="se.est")

# 3) item parameter estimation: fix the guessing parameters of the 3PLM to 0.2
(mod3 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, fix.g=TRUE, a.val.1pl=1, g.val=.2,
                  use.startval=FALSE))
summary(mod3)

# extract both item parameter and standard error estimates
getirt(mod2, what="estimates")

}

}
\references{
Baker, F. B., & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques.} CRC Press.

Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.

Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
\emph{Statistical theories of mental test scores} (pp. 397-479). Reading, MA: Addison-Wesley.

Chen, P., & Wang, C. (2016). A new online calibration method for multidimensional computerized adaptive testing.
\emph{Psychometrika, 81}(3), 674-701.

Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.
}
\seealso{
\code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{sx2_fit}},
\code{\link{traceline.est_item}}, \code{\link{getirt}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimCAT_DC.R
\docType{data}
\name{simCAT_DC}
\alias{simCAT_DC}
\title{Simulated single-item format CAT Data}
\format{
This data includes a list of length three. The first internal object is a data.frame of
the item pool consisting of 100 dichotomous items. The item parameters of the first 90 items were generated with
the IRT 2PL model and calibrated with the same model. However, the item parameters of the last 10 items were
generated with the IRT 3PL model but calibrated with the IRT 2PL model. The second internal object is the response
data set including a sparse response data set of 10,000 examinees for the items in the item pool.
The third internal object is the examinee's ability estimates for 10,000 examinees.
}
\usage{
simCAT_DC
}
\description{
This data set contains an item pool information, response data, and examinee's ability estimates.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglike_score.R
\name{llike_score}
\alias{llike_score}
\title{Loglikelihood of ability}
\usage{
llike_score(
  x,
  data,
  theta,
  D = 1,
  method = "MLE",
  norm.prior = c(0, 1),
  fence.a = 3,
  fence.b = NULL,
  missing = NA
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}.}

\item{data}{A matrix or vector containing examinees' response data for the items in the argument \code{x}. When a matrix is used, a row and column indicate
the examinees and items, respectively. When a vector is used, it should contains the item response data for an examinee.}

\item{theta}{A numeric vector of abilities of which loglikelihood values are computed.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{method}{A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation,
"MLEF" for the maximum likelihood estimation with fences, "MAP" for the maximum a posteriori estimation. Default method is "MLE".}

\item{norm.prior}{A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
c(0,1). Ignored if \code{method} is "MLE" or "MLEF".}

\item{fence.a}{A numeric value specifying the item slope parameter (i.e., \emph{a}-parameter) for the two imaginary items in MLEF. See below for details.
Default is 3.0.}

\item{fence.b}{A numeric vector of two components specifying the lower and upper fences of item difficulty parameters (i.e., \emph{b}-parameters)
for the two imaginary items, respectively, in MLEF. When \code{fence.b = NULL}, the lower and upper fences of item difficulty parameters were
automatically set. See below for details. Default is NULL.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}
}
\value{
A data frame of loglikelihood values. A row indicates the ability value where the loglikelihood is computed and
a column represents a response pattern, respectively.
}
\description{
This function computes the loglikelihood of abilities for examinees given the item parameters and response data.
}
\details{
The loglikelihood function of ability for an examinee can be computed given the item parameters and the examinee's response data for the items.
For example, if you want to examine the loglikelihood functions of abilities for two examinees given the same test items specified in the argument \code{x},
then you should provide the item response data matrix with two rows in the argument \code{data} and a vector of ability points where the loglikelihood values
need to be computed in the argument \code{theta}. Or if you want to examine the loglikelihood function of ability for an examinee given the test items
specified in the argument \code{x}, then you should provide the item response data matrix with one row (or a vector of item response data) in the argument
\code{data} and a vector of ability points where the loglikelihood values need to be computed in the argument \code{theta}.
}
\examples{
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# generate examinees' abilities from N(0, 1)
set.seed(10)
score <- rnorm(5, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

# set the ability values where the loglikelihood values are computed
theta <- seq(-3, 3, 0.5)

# compute the loglikelihood values (When MLE method is used)
llike_score(x=x, data=data, theta=theta, D=1, method="MLE")

}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape_df.R
\name{shape_df}
\alias{shape_df}
\title{Create a data frame of item metadata}
\usage{
shape_df(
  par.dc = list(a = NULL, b = NULL, g = NULL),
  par.py = list(a = NULL, d = NULL),
  item.id = NULL,
  cats,
  model,
  empty.par = FALSE
)
}
\arguments{
\item{par.dc}{A list containing three vectors of dichotomous item parameters. Namely, the item discrimination (a), item difficulty (b),
and item guessing parameters.}

\item{par.py}{A list containing a vector of polytomous item discrimination (or slope) parameters and a list of polytomous item threshold
(or step) parameters. In the list, the argument \code{a} should have a vector of slope parameters and the argument \code{d} should include
a list of threshold (or step) parameters. See below for more details.}

\item{item.id}{A character vector of item IDs. If NULL, an ID is automatically given to each item.}

\item{cats}{A vector containing the number of score categories for items.}

\item{model}{A character vector of IRT models corresponding to items. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for
dichotomous items, and "GRM" and "GPCM" for polytomous items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and
"3PLM") and "GRM" and "GPCM" represent the graded response model and (generalized) partial credit model, respectively.}

\item{empty.par}{A logical value to create an empty item meta. If TRUE, the number of score categories and corresponding IRT models should be specified
in the arguments of \code{cats} and \code{model}, respectively. In the empty item meta, the item slope parameter has a fixed value of 1, the item difficulty
(or threshold) parameter has a fixed value of 0, and the item guessing parameter has a fixed value of .2. Default is FALSE.}
}
\value{
This function returns a data frame.
}
\description{
This function creates a data frame which includes item meta (e.g., item parameter, categories, models ...) to be
used for the IRT model-data fit analysis as well as other analyses.
}
\details{
For any item where "1PLM" or "2PLM" is specified in \code{model}, the item guessing parameter will be NA. If \code{model} is
a vector of \eqn{length = 1}, the specified model is replicated across all items. As in the function \code{\link{simdat}}, it is important
to clearly specify \code{cats} according to the order of items in the test form when a data frame for a mixed-format test needs to be created.
See \code{\link{simdat}} for more details about how to specify \code{cats}.

When specifying item parameters in \code{par.dc} and \code{par.dc}, keep the order of item parameter types. For example,
in the list of \code{par.dc}, the order of items parameters should be the slope, the difficulty, and the guessing parameters.

When specifying item parameters in \code{par.dc}, note that in the list of the threshold (or step) parameters, each vector should contain
the threshold (or step) parameters for each item. When an item follows the (generalized) partial credit model, the item step parameters
are the overall item difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for each category. Thus, the number
of step parameters for item with m categories is m-1 because a step parameter for the first category does not affect the category probabilities.
}
\examples{
## a mixed-item format test form
## with five dichotomous and two polytomous items
# create a list containing the dichotomous item parameters
par.dc <- list(a=c(1.1, 1.2, 0.9, 1.8, 1.4),
               b=c(0.1, -1.6, -0.2, 1.0, 1.2),
               g=rep(0.2, 5))

# create a list containing the polytomous item parameters
par.py <- list(a=c(1.4, 0.6),
               d=list(c(0.0, -1.9, 1.2), c(0.4, -1.1, 1.5, 0.2)))

# create a numeric vector of score categories for the items
cats <- c(2, 4, 2, 2, 5, 2, 2)

# create a character vector of IRT models for the items
model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")

# create an item meta set
shape_df(par.dc=par.dc, par.py=par.py, cats=cats, model=model)

## an empty item meta with five dichotomous and two polytomous items
# create a numeric vector of score categories for the items
cats <- c(2, 4, 3, 2, 5, 2, 2)

# create a character vector of IRT models for the items
model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")

# create an empty item meta set
shape_df(cats=cats, model=model, empty.par=TRUE)

## an item meta for a single-item format test form with five dichotomous
shape_df(par.dc=par.dc, cats=rep(2, 5), model="DRM")


}
\seealso{
\code{\link{test.info}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_weight.R
\name{gen.weight}
\alias{gen.weight}
\title{Generate Weights}
\usage{
gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1, l = -4, u = 4, theta)
}
\arguments{
\item{n}{An integer identifying the number of theta (or node) values for which weights are generated. Default is 41.}

\item{dist}{A character string specifying a probability distribution from which the weights are generated. Available distributions are
"norm" for a normal distribution, "unif" for a uniform distribution, and "emp" for an empirical distribution.
When \code{dist = "norm"}, either \code{n} or \code{theta} can be specified, when \code{dist = "unif"},
only \code{n} can be used, and when \code{dist = "emp"}, only \code{theta} can be used.}

\item{mu, sigma}{A mean and standard deviation of a normal distribution.}

\item{l, u}{Lower and upper limits of a uniform distribution.}

\item{theta}{A vector of empirical theta (or node) values for which weights are generated.}
}
\value{
This function returns a data frame with two columns, where the first column has theta values (nodes) and the second column provides weights.
}
\description{
This function generates a set of weights based on a set of theta values to be used in the functions \code{\link{est_score}}
and \code{\link{sx2_fit}}.
}
\details{
When the argument \code{theta} is missing, \emph{n} weights can be generated from either the normal distribution or the uniform distribution.
Note that if \code{dist = "norm"}, gaussian quadrature points and weights from the normal distribution are generated. See
\code{gauss.quad.prob()} in the \pkg{statmod} package for more details.

When the argument \code{theta} is not missing, the weights corresponding to the provided theta values are generated. Specifically, if
\code{dist = "norm"}, normalized weights from the normal distribution are returned. If \code{dist = "emp"}, every specified theta value has the equal
values of normalized weights.
}
\examples{
## example 1
## generate 41 gaussian quadrature points and weights of normal distribution
gen.weight(n=41, dist = "norm", mu = 0, sigma = 1)

## example 2
## generate 41 theta values and weights from the uniform normal distribution,
## given the mininum value of -4 and the maximum value of 4
gen.weight(n=41, dist = "unif", l = -4, u = 4)

## example 3
## generate the normalized weights from the standardized normal distribution,
## given a set of theta values
theta <- seq(-4, 4, by=0.1)
gen.weight(dist = "norm", mu = 0, sigma = 1, theta = theta)

## example 4
## generate the same values of normalized weights for the theta values that are
## randomly sampled from the standardized normal distribution
theta <- rnorm(100)
gen.weight(dist = "emp", theta = theta)

}
\seealso{
\code{\link{est_score}}, \code{\link{sx2_fit}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lwrc.R
\name{lwrc}
\alias{lwrc}
\title{Lord-Wingersky Recursion Formula}
\usage{
lwrc(x = NULL, theta, prob = NULL, cats, D = 1)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
See \code{\link{irtfit}}, \code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}. If \code{prob = NULL}, this data frame is
used in the recursion formula. See below for details.}

\item{theta}{A vector of theta values where the conditional distribution of observed scores are computed.
The theta values are only required when a data frame is specified in the argument \code{x}.}

\item{prob}{A matrix containing the probability of answering each category of an item. Each row indicates an item and
each column represents each category of the item. When the number of categories differs between items, the empty cells
should be filled with zeros or NA values. If \code{x = NULL}, this probability matrix is used in the recursion Formula.}

\item{cats}{A numeric vector specifying the number of categories for each item. For example, a dichotomous
item has two categories. This information is only required when a probability matrix is specified in the argument
\code{prob}.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function
(if set to 1.7). Default is 1.}
}
\value{
When the \code{prob} argument is provided, this function returns a vector of the probabilities of obtaining every 
observed score on a test. When the \code{x} argument is specified, the function returns a matrix of conditional probabilities 
across all possible observed scores and theta values.
}
\description{
This function computes the conditional distributions of number-correct (or observed) scores
given probabilities of category responses to items or given a set of theta values using Lord and
Wingersky recursion formula (1984).
}
\details{
The Lord and Wingersky recursive algorithm is an efficient way of calculating the compound probabilities
of any number-correct scores on a test based on IRT models. This algorithm is particularly useful when computing
the IRT model-based observed score distribution for a test.

To compute the conditional distributions of observed scores, either the item metadata set specified in \code{x} or
the probability matrix specified in \code{prob} can be used.
}
\examples{
## example 1: when a matrix of probabilities is used as a data set
## this is an example from Kolen and Brennan (2004, p. 183)
# create a matrix of probabilities of getting correct and incorrect answers for three items
probs <- matrix(c(.74, .73, .82, .26, .27, .18), nrow=3, ncol=2, byrow = FALSE)

# create a vector of score categories for the three items
cats <- c(2,2,2)

# compute the conditional distributions of observed scores
lwrc(prob=probs, cats=cats)

## example 2: when a matrix of probabilities is used as a data set
## with a mixed-format test
# category probabilities for a dichotomous item
p1 <- c(0.2, 0.8, 0, 0, 0)
# category probabilities for a dichotomous item
p2 <- c(0.4, 0.6, NA, NA, NA)
# category probabilities for a polytomous item with five categories
p3 <- c(0.1, 0.2, 0.2, 0.4, 0.1)
# category probabilities for a polytomous item with three categories
p4 <- c(0.5, 0.3, 0.2, NA, NA)

# rbind the probability vectors
p <- rbind(p1, p2, p3, p4)

# create a vector of score categories for the four items
cats <- c(2, 2, 5, 3)

# compute the conditional distributions of observed scores
lwrc(prob=p, cats=cats)

## example 3: when a data frame for the item metadata is used instead of a probabiliy matrix
## with a mixed-format test
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# compute the conditional distributions of observed scores
lwrc(x=x, theta=seq(-1, 1, 0.1), D=1)

}
\references{
Kolen, M. J. & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
Springer.

Lord, F. & Wingersky, M. (1984). Comparison of IRT true score and equipercentile observed score equatings.
\emph{Applied Psychological Measurement, 8}(4), 453-461.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_info.R
\name{plot.test.info}
\alias{plot.test.info}
\title{Plot Item and Test Information Functions}
\usage{
\method{plot}{test.info}(
  x,
  item.loc = NULL,
  overlap = FALSE,
  csee = FALSE,
  xlab.text,
  ylab.text,
  main.text,
  lab.size = 15,
  main.size = 15,
  axis.size = 15,
  line.color,
  line.size = 1,
  layout.col = 4,
  strip.size = 12,
  ...
)
}
\arguments{
\item{x}{x An object of class \code{\link{test.info}}.}

\item{item.loc}{A vector of numeric values indicating that the item information functions of the \emph{n}th items
(or the location of items in a test form) are plotted. If NULL, the test information function for the total test form is drawn.
Default is NULL.}

\item{overlap}{Logical value indicating whether multiple item information functions are plotted in one panel. 
If FALSE, multiple item information functions are displayed in multiple panels, one for each.}

\item{csee}{Logical value indicating whether the function displays the conditional standard error of estimation (CSEE) at a test level. 
If FALSE, item/test information function is plotted. Note that the CSEE plot is displayed only at a test level.}

\item{xlab.text, ylab.text}{A title for the x and y axes.}

\item{main.text}{An overall title for the plot.}

\item{lab.size}{The size of xlab and ylab. Default is 15.}

\item{main.size}{The size of \code{main.text}. Default is 15.}

\item{axis.size}{The size of labels along the x and y axes. Default is 15.}

\item{line.color}{A character string specifying a color for the line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}
for more details about colors used in ggplot2.}

\item{line.size}{The size of lines. Default is 1.}

\item{layout.col}{An integer value indicating the number of columns in the panel when displaying the item information functions of
the multiple items. Default is 4.}

\item{strip.size}{The size of facet labels when the item information functions of the multiple items are drawn.}

\item{...}{Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.}
}
\description{
This function plots item or test information function given a specified theta values. In addition, 
this function displays conditional standard errors at a test level.
}
\details{
All of the plots are drawn using the ggplot2 package.
The object of class \code{\link{test.info}} can be obtained from the function \code{\link{test.info}}.
}
\examples{
## the use of a "-prm.txt" file obtained from a flexMIRT
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# set theta values
theta <- seq(-4, 4, 0.1)

# compute item and test information values given the theta values
x <- test.info(x=test_flex, theta=theta, D=1)

# draw a plot of the test information function
plot(x)

# draw a plot of the item information function for the second item
plot(x, item.loc=2)

# draw a plot of multiple item information functions across the multiple panels
plot(x, item.loc=1:8, overlap=FALSE)

# draw a plot of multiple item information functions across in one panel
plot(x, item.loc=1:8, overlap=TRUE)

# draw a plot of conditional standard error at a test level
plot(x, csee=TRUE)

}
\seealso{
\code{\link{test.info}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covirt.R
\name{covirt}
\alias{covirt}
\title{Asymptotic variance-covariance matrices of item parameter estimates}
\usage{
covirt(
  x,
  D = 1,
  nstd = 1000,
  pcm.loc = NULL,
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{nstd}{An integer value or a vector of integer values indicating a sample size. When a vector is specified, length of the vector must be
the same as the number of test items in the argument \code{x}. Default is 1,000. See below for details.}

\item{pcm.loc}{A vector of integer values indicating the locations of partial credit model (PCM) items. For the PCM items,
the variance-covariance matrices are computed only for the item category difficulty parameters. Default is NULL. See below for details.}

\item{norm.prior}{A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution.
Default is c(0,1).}

\item{nquad}{An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.}

\item{weights}{A two-column matrix or data frame containing the theta values (in the first column) and the weights (in the second column)
for the prior distribution. The weights and theta values can be easily obtained using the function \code{\link{gen.weight}}.
If NULL, default values are used for the prior distribution (see the arguments of \code{norm.prior} and \code{nquad}). Default is NULL.}
}
\value{
A list of two internal objects. The first internal object contains a list of the variance-covariance matrices of item parameter estimates.
The second internal object contains a list of the standard errors of item parameter estimates.
}
\description{
This function calculates the analytical asymptotic variance-covariance matrices (e.g., Li & Lissitz, 2004; Thissen & Wainer, 1982)
of item parameter estimates for dichotomous and polytomous IRT Models without examinee's responses to test items, 
given a set of item parameter estimates and sample size. The square roots of variance terms in the matrices can be used as the asymptotic 
standard errors of maximum likelihood item parameter estimates.
}
\details{
The standard errors obtained from the analytical approach are likely to represent lower bounds for the actual standard errors (Thissen & Wainer, 1982).
Therefore, they may be useful for assessing the degree of precision of a set of item parameter estimates when the corresponding standard errors of 
the estimates are not presented in literature or research reports.

Sometimes item parameters need to be estimated using different sample size. If the item parameters in the argument \code{x} were
calibrated with different number of examinees, a vector of different sample sizes should be specified in the argument \code{nstd}. Suppose
that you want to compute the variance-covariance matrices of five IRT 3PLM items and the five items were calibrated with 500, 600, 1,000, 2,000,
and 700 examinees, respectively. Then, \code{nstd = c(500, 600, 1000, 2000, 700)} must be specified.

Because you can specify only "GPCM" for both the partial credit model (PCM) or the generalized partial credit model (GPCM) in the item metadata,
you must indicate which items are the PCM items through the argument \code{pcm.loc}. This is because the item category difficulty parameters are estimated
from the PCM, meaning that the variance-covariance of item parameter estimates must be computed for the item category difficulty parameters. Suppose
that you want to compute the variance-covariance matrices of five polytomous items and the last two items were calibrated with the PCM. Then,
\code{pcm.loc = c(4, 5)} must be specified.
}
\examples{
## the use of a "-prm.txt" file obtained sfrom a flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the first two dichotomous items and last polytomous item
x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df[c(1:2, 55), ]

# compute the var-covariance matrices with sample size of 2,000
covirt(x, D=1, nstd=2000, norm.prior=c(0, 1), nquad=40)

}
\references{
Li, Y. & Lissitz, R. (2004). Applications of the analytically derived asymptotic standard errors of item response theory
item parameter estimates. \emph{Journal of educational measurement, 41}(2), 85-117.

Thissen, D. & Wainer, H. (1982). Weighted likelihood estimation of ability in item response theory. 
\emph{Psychometrika, 54}(3), 427-450.
}
\seealso{
\code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtfit.R
\name{irtfit}
\alias{irtfit}
\alias{irtfit.default}
\alias{irtfit.est_item}
\alias{irtfit.est_irt}
\title{Traditional IRT item fit statistics}
\usage{
irtfit(x, ...)

\method{irtfit}{default}(
  x,
  score,
  data,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  D = 1,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  ...
)

\method{irtfit}{est_item}(
  x,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  ...
)

\method{irtfit}{est_irt}(
  x,
  score,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  ...
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}. See below for details.}

\item{...}{Further arguments passed to or from other methods.}

\item{score}{A vector of examinees' ability estimates.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{group.method}{A character string indicating how to group examinees along the ability scale for computing the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
Available methods are "equal.width" for grouping examinees by dividing the ability scale into intervals of equal width and "equal.freq"
for grouping examinees by dividing the ability scale into intervals with equal frequencies of examinees. However, "equal.freq" does not
always guarantee exactly the same frequency of examinees for all groups. Default is "equal.width". To divide the ability scale, the range
of ability scale and the number of divided groups must be specified in the arguments of \code{range.score} and \code{n.width}, respectively.
See below for details.}

\item{n.width}{An integer value to specify the number of divided groups along the ability scale. Default is 10. See below for details.}

\item{loc.theta}{A character string to indicate the location of ability point at each group (or interval) where the expected probabilities
of score categories are calculated using the IRT models. Available locations are "average" for computing the expected probability
at the average point of examinees' ability estimates in each group and "middle" for computing the expected probability at the midpoint of each group.
Default is "average".}

\item{range.score}{A vector of two numeric values to restrict the range of ability scale. All ability estimates less than
the first value are transformed to the first value. All ability estimates greater than the second value are transformed to the second value.
If NULL, the minimum and maximum values of ability estimates in the argument \code{score} is used as the range of ability scale. Note that
selection of grouping method in the argument \code{group.method} has nothing to do with the range of ability scale. Default is NULL.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{alpha}{A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test for the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
Default is .05.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}

\item{overSR}{A numeric value to specify a criterion to find ability groups (or intervals) which have standardized residuals
greater than the specified value. Default is 2.}

\item{min.collapse}{An integer value to indicate the minimum frequency of cells to be collapsed when computing the \eqn{\chi^{2}} and \eqn{G^{2}}
fit statistics. Neighboring interval groups will be collapsed to avoid expected interval frequencies less than the specified minimum cell frequency.
Default is 1.}
}
\value{
This function returns an object of class \code{\link{irtfit}}. Within this object, several internal objects are contained such as:
\item{fit_stat}{A data frame containing the results of three IRT fit statistics (i.e., \eqn{\chi^{2}} and \eqn{G^{2}}, infit, outfit statistics) across
all evaluated items. In the data frame, the columns indicate item's ID, \eqn{\chi^{2}} fit statistic, \eqn{G^{2}} fit statistic, degrees of freedom for the \eqn{\chi^{2}},
degrees of freedom for the \eqn{G^{2}}, critical value for the \eqn{\chi^{2}}, critical value for the \eqn{G^{2}}, p-value for the \eqn{\chi^{2}},
p-value for the \eqn{G^{2}}, outfit statistic, infit statistic, the number of examinees used to compute the five fit statistics, and the proportion of
ability groups (or intervals), before collapsing the cells, that have standardized residuals greater than the specified criterion in the argument \code{overSR},
respectively.}
\item{contingency.fitstat}{A list of contingency tables used to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics for all items.
Note that the collapsing cell strategy is implemented to these contingency tables.}
\item{contingency.plot}{A list of contingency tables used to draw a raw and standardized residual plots (Hambleton et al., 1991) in the function of
\code{\link{plot.irtfit}}. Note that the collapsing cell strategy is \emph{not} implemented to these contingency tables.}
\item{individual.info}{A list of data frames including individual residual and variance values. Those information are used to compute
infit and outfit statistics.}
\item{item_df}{The item metadata specified in the argument \code{x}.}
\item{ancillary}{A list of ancillary information used in the item fit analysis.}
}
\description{
This function computes traditional IRT item fit statistics (i.e., \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
loglikelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley & Mills, 1985), and infit and outfit statistics (Ames et al., 2015)) and returns
contingency tables to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics. Note that caution is needed in interpreting the infit and
outfit statistics for non-Rasch models. The saved object of this function, especially the object of contingency tables,
is used in the function of \code{\link{plot.irtfit}} to draw a raw and standardized residual plots (Hambleton et al., 1991).
}
\details{
A specific form of a data frame should be used for the argument \code{x}. The first column should have item IDs,
the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
In the \pkg{irtplay} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as 
the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item. 
Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because 
the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories, 
four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
  ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
  ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
  ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
}
And an example of a data frame for a mixed-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
  ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
  ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
  ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
  ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
  ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
  ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
}
See \code{IRT Models} section in the page of \code{\link{irtplay-package}} for more details about the IRT models used in the \pkg{irtplay} package. 
An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.

To calculate the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics, two methods are used in the argument \code{group.method} to divide the ability scale
into several groups. If \code{group.method = "equal.width"}, the examinees are grouped based on equal length of intervals.
If \code{group.method = "equal.freq"}, the examinees are grouped so that all groups have equal frequencies. However, the grouping method
of "equal.freq" does guarantee that every group has the exactly same frequency of examinees. This is because the examinees are divided by
the same size of quantile.

When dividing the ability scale into intervals to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics, the intervals should be wide enough not to include
too small number of examinees. On the other hand, the interval should be narrow enough to include homogeneous examinees in terms of ability
(Hambleton et al, 1991). Thus, if you want to divide the ability scale into other than ten groups, you need to specify the number of groups
in the argument \code{n.width}. Yen (1981) fixed the number of groups to 10, whereas Bock (1960) allowed for any number of groups.

Regarding degrees of freedom (\emph{df}), the \eqn{\chi^{2}} is assumed to be distributed approximately as a chi-square with \emph{df} equal to
the number of groups less the number of the IRT model parameters (Ames et al., 2015) whereas the \eqn{G^{2}} is assumed to be distributed approximately
as a chi-square with \emph{df} equal to the number of groups (Ames et al., 2015; Muraki & Bock, 2003)

Note that if "DRM" is specified for an item in the item metadata set, the item is considered as "3PLM" to compute degrees of freedom of
the \eqn{\chi^{2}} fit statistic.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Default method to compute the traditional IRT item fit statistics for a data frame \code{x} containing the item metadata.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.

\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.
}}

\examples{
\donttest{
## example 1
## use the simulated CAT data
# find the location of items that have more than 10,000 responses
over10000 <- which(colSums(simCAT_MX$res.dat, na.rm=TRUE) > 10000)

# select the items that have more than 10,000 responses
x <- simCAT_MX$item.prm[over10000, ]

# select the response data for the items
data <- simCAT_MX$res.dat[, over10000]

# select the examinees' abilities
score <- simCAT_MX$score

# compute fit statistics
fit1 <- irtfit(x=x, score=score, data=data, group.method="equal.width",
               n.width=10, loc.theta="average", range.score=NULL, D=1, alpha=0.05,
               missing=NA, overSR=2)

# fit statistics
fit1$fit_stat

# contingency tables
fit1$contingency.fitstat


## example 2
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the first two dichotomous items and last polytomous item
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]

# generate examinees' abilities from N(0, 1)
set.seed(10)
score <- rnorm(1000, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

# compute fit statistics
fit2 <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
               n.width=11, loc.theta="average", range.score=c(-4, 4), D=1, alpha=0.05)

# fit statistics
fit2$fit_stat

# contingency tables
fit2$contingency.fitstat

# residual plots for the first item (dichotomous item)
plot(x=fit2, item.loc=1, type = "both", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)

# residual plots for the third item (polytomous item)
plot(x=fit2, item.loc=3, type = "both", ci.method = "wald", show.table=FALSE, ylim.sr.adjust=TRUE)

}

}
\references{
Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional Module on Item-Fit Statistics for Item Response Theory Models.
\emph{Educational Measurement: Issues and Practice, 34}(3), 39-48.

Bock, R.D. (1960), \emph{Methods and applications of optimal scaling}. Chapel Hill, NC: L.L. Thurstone Psychometric Laboratory.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.

McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
\emph{Applied Psychological Measurement, 9}, 49-57.

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com

Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric procedure for assessing goodness-of-fit in
item response theory. \emph{Applied Measurement in Education, 21}(1), 22-40.

Yen, W. M. (1981). Using simulation results to choose a latent trait model. \emph{Applied Psychological Measurement, 5}, 245-262.
}
\seealso{
\code{\link{plot.irtfit}}, \code{\link{shape_df}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_irt.R
\name{est_irt}
\alias{est_irt}
\title{Item parameter estimation using MMLE-EM algorithm}
\usage{
est_irt(
  x = NULL,
  data,
  D = 1,
  model = NULL,
  cats = NULL,
  item.id = NULL,
  fix.a.1pl = FALSE,
  fix.a.gpcm = FALSE,
  fix.g = FALSE,
  a.val.1pl = 1,
  a.val.gpcm = 1,
  g.val = 0.2,
  use.aprior = FALSE,
  use.bprior = FALSE,
  use.gprior = TRUE,
  aprior = list(dist = "lnorm", params = c(0, 0.5)),
  bprior = list(dist = "norm", params = c(0, 1)),
  gprior = list(dist = "beta", params = c(5, 16)),
  missing = NA,
  Quadrature = c(49, 6),
  weights = NULL,
  group.mean = 0,
  group.var = 1,
  EmpHist = FALSE,
  use.startval = FALSE,
  Etol = 1e-04,
  MaxE = 500,
  control = list(iter.max = 200),
  fipc = FALSE,
  fipc.method = "MEM",
  fix.loc = NULL,
  verbose = TRUE
)
}
\arguments{
\item{x}{A data frame containing the item metadata. This metadata is necessary to obtain the information of
each item (i.e., number of score categories and IRT model) to be calibrated. You can easily create an empty
item metadata using the function \code{\link{shape_df}}. When \code{use.startval = TRUE}, the item parameters
specified in the item metadata are used as the starting values for the item parameter estimation.
If \code{x = NULL}, the arguments of \code{model} and \code{cats} must be specified. Note that when \code{fipc = TRUE}
to implement the FIPC method, the item metadata of a test form must be provided in the argument \code{x}.
See below for details.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{model}{A vector of character strings indicating what IRT model is used to calibrate each item. Available IRT models are
"1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function.
If a single character of the IRT model is specified, that model will be recycled across all items. This information is only required
when \code{x = NULL} and \code{fipc = FALSE}.}

\item{cats}{A numeric vector specifying the number of score categories for each item. For example, a dichotomous
item has two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items
are binary items (i.e., dichotomous items), it assumes that all items have two score categories. This information is only required
when \code{x = NULL} and \code{fipc = FALSE}.}

\item{item.id}{A character vector of item IDs. If NULL, the item IDs are generated automatically. When \code{fipc = TRUE} and the Item IDs
are given by the \code{item.id} argument, the Item IDs in the \code{x} argument are overridden. Default is NULL.}

\item{fix.a.1pl}{A logical value. If TRUE, the slope parameters of the 1PLM items are fixed to a specific value specified in the argument
\code{a.val.1pl}. Otherwise, the slope parameters of all 1PLM items are constrained to be equal and estimated. Default is FALSE.}

\item{fix.a.gpcm}{A logical value. If TRUE, the GPCM items are calibrated with the partial credit model and the slope parameters of
the GPCM items are fixed to a specific value specified in the argument \code{a.val.gpcm}. Otherwise, the slope parameter of each GPCM item
is estimated. Default is FALSE.}

\item{fix.g}{A logical value. If TRUE, the guessing parameters of the 3PLM items are fixed to a specific value specified in the argument
\code{g.val}. Otherwise, the guessing parameter of each 3PLM item is estimated. Default is FALSE.}

\item{a.val.1pl}{A numeric value. This value is used to fixed the slope parameters of the 1PLM items.}

\item{a.val.gpcm}{A numeric value. This value is used to fixed the slope parameters of the GPCM items.}

\item{g.val}{A numeric value. This value is used to fixed the guessing parameters of the 3PLM items.}

\item{use.aprior}{A logical value. If TRUE, a prior distribution for the slope parameters is used for the parameter calibration
across all items. Default is FALSE.}

\item{use.bprior}{A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used for the parameter calibration
across all items. Default is FALSE.}

\item{use.gprior}{A logical value. If TRUE, a prior distribution for the guessing parameters is used for the parameter calibration
across all 3PLM items. Default is TRUE.}

\item{aprior}{A list containing the information of the prior distribution for item slope parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{bprior}{A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{gprior}{A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
and \code{dnorm()} in the \pkg{stats} package for more details.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}

\item{Quadrature}{A numeric vector of two components specifying the number of quadrature points (in the first component) and
the symmetric minimum and maximum values of these points (in the second component). For example, a vector of c(49, 6) indicates 49 rectangular
quadrature points over -6 and 6. The quadrature points are used in the E step of the EM algorithm. Default is c(49, 6).}

\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
using the function \code{\link{gen.weight}}. If NULL, a normal prior density is used based on the information provided in the arguments
of \code{Quadrature}, \code{group.mean}, and \code{group.var}). Default is NULL.}

\item{group.mean}{A numeric value to set the mean of latent variable prior distribution. Default is 0. This value is fixed to remove
the indeterminancy of item parameter scale when calibrating items. However, the scale of prior distribution is updated when FIPC is implemented.}

\item{group.var}{A positive numeric value to set the variance of latent variable prior distribution. Default is 1. This value is fixed to remove
the indeterminancy of item parameter scale when calibrating items. However, the scale of prior distribution is updated when FIPC is implemented.}

\item{EmpHist}{A logical value. If TRUE, the empirical histogram of the latent variable prior distribution is simultaneously estimated with
the item parameters using Woods's (2007) approach. The item parameters are calibrated against the estimated empirical histogram prior distribution.
See below for details.}

\item{use.startval}{A logical value. If TRUE, the item parameters provided in the item metadata (i.e., the argument \code{x}) are used as
the starting values for the item parameter estimation. Otherwise, internal starting values of this function are used. Default is FALSE.}

\item{Etol}{A positive numeric value. This value sets the convergence criterion for E steps of the EM algorithm. Default is 1e-4.}

\item{MaxE}{A positive integer value. This value determines the maximum number of the E steps in the EM algorithm. Default is 500.}

\item{control}{A list of control parameters to be passed to the optimization function of \code{nlminb()} in the \pkg{stats} package. The control parameters
set the conditions of M steps of the EM algorithm. For example, the maximum number of iterations in each of the iterative M steps can
be set by \code{control = list(iter.max=200)}. Default maximum number of iterations in each M step is 200. See \code{nlminb()} in the \pkg{stats} package
for other control parameters.}

\item{fipc}{A logical value. If TRUE, FIPC is implemented for item parameter estimation. See below for details.}

\item{fipc.method}{A character string specifying the FIPC method. Available methods include "OEM" for "No Prior Weights Updating and One EM Cycle
(NWU-OEM; Wainer & Mislevy, 1990)" and "MEM" for "Multiple Prior Weights Updating and Multiple EM Cycles (MWU-MEM; Kim, 2006)."
When \code{fipc.method = "OEM"}, the maximum number of the E steps of the EM algorithm is set to 1 no matter what number is specified
in the argument \code{MaxE}.}

\item{fix.loc}{A vector of positive integer values specifying the location of the items to be fixed in the item metadata (i.e., \code{x})
when the FIPC is implemented. For example, suppose that five items located in the 1st, 2nd, 4th, 7th, and 9th rows of the item metadata \code{x}
should be fixed. Then \code{fix.loc = c(1, 2, 4, 7, 9)}.}

\item{verbose}{A logical value. If FALSE, all progress messages including the process information on the EM algorithm are suppressed.
Default is TRUE.}
}
\value{
This function returns an object of class \code{\link{est_irt}}. Within this object, several internal objects are contained such as:
\item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
\item{par.est}{A data frame containing the item parameter estimates.}
\item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
observed information functions. The standard errors are estimated using the cross-production approximation method (Meilijson, 1989).}
\item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
when interpreting the variance-covariance matrix of item parameter estimates.}
\item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
\item{loglikelihood}{A sum of the log-likelihood values of the observed data set (marginal log-likelihood) across all estimated items.}
\item{aic}{A model fit statistic of Akaike information criterion based on the loglikelihood.}
\item{bic}{A model fit statistic of Bayesian information criterion based on the loglikelihood.}
\item{group.par}{A data frame containing the mean, variance, and standard deviation of latent variable prior distribution.}
\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the (updated) latent variable prior distribution.}
\item{posterior.dist}{A matrix of normalized posterior densities for all the response patterns at each of the quadrature points.
The row and column indicate the response pattern and the quadrature point, respectively.}
\item{data}{A data.frame of the examinees' response data set.}
\item{scale.D}{A scaling factor in IRT models.}
\item{ncase}{A total number of response patterns.}
\item{nitem}{A total number of items included in the response data.}
\item{Etol}{A convergence criteria for E steps of the EM algorithm.}
\item{MaxE}{The maximum number of E steps in the EM algorithm.}
\item{aprior}{A list containing the information of the prior distribution for item slope parameters.}
\item{gprior}{A list containing the information of the prior distribution for item guessing parameters.}
\item{npar.est}{A total number of the estimated parameters.}
\item{niter}{The number of EM cycles completed.}
\item{maxpar.diff}{A maximum item parameter change when the EM cycles were completed.}
\item{EMtime}{Time (in seconds) spent for the EM cycles.}
\item{SEtime}{Time (in seconds) spent for computing the standard errors of the item parameter estimates.}
\item{TotalTime}{Time (in seconds) spent for total compuatation.}
\item{test.1}{Status of the first-order test to report if the gradients has vanished sufficiently for the solution to be stable.}
\item{test.2}{Status of the second-order test to report if the information matrix is positive definite, which is a prerequisite
for the solution to be a possible maximum.}
\item{var.note}{A note to report if the variance-covariance matrix of item parameter estimates is obtainable from the information matrix.}
\item{fipc}{A logical value to indicate if FIPC was used.}
\item{fipc.method}{A method used for the FIPC.}
\item{fix.loc}{A vector of integer values specifying the location of the fixed items when the FIPC was implemented.}

The internal objects can be easily extracted using the function \code{\link{getirt}}.
}
\description{
This function fits unidimensional item response (IRT) models to a mixture of dichotomous and polytomous data using
marginal maximum likelihood estimation with expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin, 1981). This function also
implements the fixed item parameter calibration (FIPC; Kim, 2006). As Method A (Stocking, 1988), FIPC is one of useful online item
calibration methods for computerized adaptive testing (CAT) to put the parameter estimates of pretest items on the same scale of
operational item parameter estimates (Ban, Hanson, Wang, Yi, & Harris, 2001). For dichotomous items, IRT one-, two-, and three-parameter
logistic models are available. For polytomous items, the graded response model (GRM) and the (generalized) partial credit model (GPCM)
are available.
}
\details{
A specific form of a data frame should be used for the argument \code{x}. The first column should have item IDs,
the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
In the \pkg{irtplay} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as 
the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item. 
Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because 
the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories, 
four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
  ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
  ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
  ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
}
And an example of a data frame for a mixed-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
  ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
  ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
  ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
  ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
  ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
  ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
}
See \code{IRT Models} section in the page of \code{\link{irtplay-package}} for more details about the IRT models used in the \pkg{irtplay} package. 
An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.

To fit the IRT models to data, the IRT model and the number of score category information for the estimated items must be provided as well as
the item response data. There are two way to provide the IRT model and score category information. The first way is to provide the item metadata
to the argument \code{x}. As explained above, the item metadata can be easily created by the function \code{\link{shape_df}}. The second way is
specify the IRT models and the score category information into the arguments of \code{model} and \code{cats}. Thus, if \code{x=NULL}, the specified
information in \code{model} and \code{cats} are used.

To implement FIPC, however, the item metadata must be provided in the argument \code{x}. This is because the item parameters of the fixed items
in the item metadata are used to estimate the characteristic of the underlying latent variable prior distribution when calibrating the rest of freely estimated items.
More specifically, the underlying latent variable prior distribution of the fixed items is estimated during the calibration of the freely estimated items
to put the item parameters of the freely estimated items on the scale of the fixed item parameters (Kim, 2006).

In terms of approaches for FIPC, Kim (2006) described five different methods. Among them, two methods are available in the
function \code{\link{est_irt}}. The first method is "NWU-OEM" where uses just one E step in the EM algorithm, involving data from only the fixed items, and
just one M step, involving data from only non-fixed items. This method is suggested by Wainer and Mislevy (1990) in the context of online calibration. This method
can be implemented by setting \code{fipc.method = "OEM"}. The second method is "MWU-MEM" which iteratively updates the latent variable prior distribution and
finds the parameter estimates of the non-fixed items. In this method, the same procedure of NWU-OEM method is applied to the first EM cycle. From the second
EM cycle, both the parameters of non-fixed items and the weights of the prior distribution are concurrently updated. This method can be implemented by 
setting \code{fipc.method = "MEM"}. See Kim (2006) for more details.

When \code{EmpHist = TRUE}, the empirical histogram of latent variable prior distribution is simultaneously estimated with the item parameters. If \code{fipc = TRUE}
given \code{EmpHist = TRUE}, the scale parameters (e.g., mean and variance) of the empirical prior distribution are estimated as well. If \code{fipc = FALSE} given
\code{EmpHist = TRUE}, the scale parameters of the empirical prior distribution are fixed to the values specified in the arguments of \code{group.mean} and \code{group.var}.
When \code{EmpHist = FALSE}, the normal prior distribution is used during the item parameter estimation. If \code{fipc = TRUE} given \code{EmpHist = FALSE},
the scale parameters of the normal prior distribution are estimated as well as the item parameters. If \code{fipc = FALSE} given \code{EmpHist = FALSE},
the scale parameters of the normal prior distribution are fixed to the values specified in the arguments of \code{group.mean} and \code{group.var}.
}
\examples{
\donttest{

##------------------------------------------------------------------------------
# 1. item parameter estimation for the dichotomous item data (LSAT6)
##------------------------------------------------------------------------------
# fit the 1PL model to LSAT6 data and constrain the slope parameters to be equal
(mod.1pl.c <- est_irt(data=LSAT6, D=1, model="1PLM", cats=2, fix.a.1pl=FALSE))

# summary of the estimation
summary(mod.1pl.c)

# extract the item parameter estimates
getirt(mod.1pl.c, what="par.est")

# extract the standard error estimates
getirt(mod.1pl.c, what="se.est")

# fit the 1PL model to LSAT6 data and fix the slope parameters to 1.0
(mod.1pl.f <- est_irt(data=LSAT6, D=1, model="1PLM", cats=2, fix.a.1pl=TRUE, a.val.1pl=1))

# summary of the estimation
summary(mod.1pl.f)

# fit the 2PL model to LSAT6 data
(mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2))

# summary of the estimation
summary(mod.2pl)

# assess the fit of the 2PL model to the LSAT5 data using S-X2 fit statistic
(sx2fit.2pl <- sx2_fit(x=mod.2pl))

# compute the item and test information at several theta points
theta <- seq(-4, 4, 0.1)
(info.2pl <- test.info(x=mod.2pl, theta=theta))

# draw the test characteristic curve plot
(trace.2pl <- traceline(x=mod.2pl, theta=theta))
plot(trace.2pl)

# draw the item characteristic curve for the 1st item
plot(trace.2pl, item.loc=1)

# fit the 2PL model to LSAT6 data and
# estimate the empirical histogram of latent variable prior distribution
# also use a less stringent convergence criterion for E-step
(mod.2pl.hist <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2, EmpHist=TRUE, Etol=0.001))
(emphist <- getirt(mod.2pl.hist, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# fit the 3PL model to LSAT6 data and use the Beta prior distribution for
# the guessing parameters
(mod.3pl <- est_irt(data=LSAT6, D=1, model="3PLM", cats=2, use.gprior=TRUE,
                    gprior=list(dist="beta", params=c(5, 16))))

# summary of the estimation
summary(mod.3pl)

# fit the 3PL model to LSAT6 data, but fix the guessing parameters to be 0.2
(mod.3pl.f <- est_irt(data=LSAT6, D=1, model="3PLM", cats=2, fix.g=TRUE, g.val=0.2))

# summary of the estimation
summary(mod.3pl.f)

# fit the differnt dichotomous models to each item of LSAT6 data
# fit the constrained 1PL model to the 1st, 2nd, and 3rd items, fit the 2PL model to
# the 4th item, and fit the 3PL model to the 5th item with the Beta prior of
# the guessing parameter
(mod.drm.mix <- est_irt(data=LSAT6, D=1, model=c("1PLM", "1PLM", "1PLM", "2PLM", "3PLM"),
                        cats=2, fix.a.1pl=FALSE, use.gprior=TRUE,
                        gprior=list(dist="beta", params=c(5, 16))))
# summary of the estimation
summary(mod.drm.mix)

##------------------------------------------------------------------------------
# 2. item parameter estimation for the mixed-item format data (simulation data)
##------------------------------------------------------------------------------
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# modify the item metadata so that the 39th and 40th items follow GPCM
x[39:40, 3] <- "GPCM"

# generate 1,000 examinees' latent abilities from N(0, 1)
set.seed(37)
score1 <- rnorm(1000, mean=0, sd=1)

# simulate the response data
sim.dat1 <- simdat(x=x, theta=score1, D=1)

# fit the 3PL model to all dichotomous items, fit the GPCM model to 39th and 40th items,
# and fit the GRM model to the 53th, 54th, 55th items.
# use the beta prior distribution for the guessing parameters, use the log-normal
# prior distribution for the slope parameters, and use the normal prior distribution
# for the difficulty (or threshold) parameters.
# also, specify the argument 'x' to provide the IRT model and score category information
# for items
item.meta <- shape_df(item.id=x$id, cats=x$cats, model=x$model, empty.par=TRUE)
(mod.mix1 <- est_irt(x=item.meta, data=sim.dat1, D=1, use.aprior=TRUE, use.bprior=TRUE,
                     use.gprior=TRUE,
                     aprior=list(dist="lnorm", params=c(0.0, 0.5)),
                     bprior=list(dist="norm", params=c(0.0, 2.0)),
                     gprior=list(dist="beta", params=c(5, 16))))

# summary of the estimation
summary(mod.mix1)

# estimate examinees' latent scores given the item parameter estimates using the MLE
(score.mle <- est_score(x=mod.mix1, method = "MLE", range = c(-4, 4), ncore=2))

# compute the traditional fit statistics
(fit.mix1 <- irtfit(x=mod.mix1, score=score.mle$est.theta, group.method="equal.width",
                    n.width=10, loc.theta="middle"))

# residual plots for the first item (dichotomous item)
plot(x=fit.mix1, item.loc=1, type = "both", ci.method = "wald",
     show.table=TRUE, ylim.sr.adjust=TRUE)

# residual plots for the last item (polytomous item)
plot(x=fit.mix1, item.loc=55, type = "both", ci.method = "wald",
     show.table=FALSE, ylim.sr.adjust=TRUE)

# fit the 2PL model to all dichotomous items, fit the GPCM model to 39th and 40th items,
# and fit the GRM model to the 53th, 54th, 55th items.
# also, specify the arguments of 'model' and 'cats' to provide the IRT model and
# score category information for items
(mod.mix2 <- est_irt(data=sim.dat1, D=1,
                     model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
                     cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3))))

# summary of the estimation
summary(mod.mix2)

# fit the 2PL model to all dichotomous items, fit the GPCM model to 39th and 40th items,
# fit the GRM model to the 53th, 54th, 55th items, and estimate the empirical histogram
# of latent variable prior distribution.
# also, specify the arguments of 'model' and 'cats' to provide the IRT model and
# score category information for items
(mod.mix3 <- est_irt(data=sim.dat1, D=1,
                     model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
                     cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)), EmpHist=TRUE))
(emphist <- getirt(mod.mix3, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# fit the 2PL model to all dichotomous items,
# fit the PCM model to 39th and 40th items by fixing the slope parameters to 1,
# and fit the GRM model to the 53th, 54th, 55th items.
# also, specify the arguments of 'model' and 'cats' to provide the IRT model and
# score category information for items
(mod.mix4 <- est_irt(data=sim.dat1, D=1,
                     model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
                     cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)),
                     fix.a.gpcm=TRUE, a.val.gpcm=1))

# summary of the estimation
summary(mod.mix4)

##------------------------------------------------------------------------------
# 3. fixed item parameter calibration (FIPC) for the mixed-item format data
#    (simulation data)
##------------------------------------------------------------------------------
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# generate 1,000 examinees' latent abilities from N(0.4, 1.3)
set.seed(20)
score2 <- rnorm(1000, mean=0.4, sd=1.3)

# simulate the response data
sim.dat2 <- simdat(x=x, theta=score2, D=1)

# fit the 3PL model to all dichotomous items, fit the GRM model to all polytomous data,
# fix the five 3PL items (1st - 5th items) and three GRM items (53rd to 55th items)
# also, estimate the empirical histogram of latent variable
# use the MEM method.
fix.loc <- c(1:5, 53:55)
(mod.fix1 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
                     gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE,
                     Etol=1e-3, fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
(prior.par <- mod.fix1$group.par)
(emphist <- getirt(mod.fix1, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# summary of the estimation
summary(mod.fix1)

# fit the 3PL model to all dichotomous items, fit the GRM model to all polytomous data,
# fix the five 3PL items (1st - 5th items) and three GRM items (53rd to 55th items)
# at this moment, do estimate the empirical histogram of latent variable.
# instead, estimate the scale of normal prior distribution of latent variable
# use the MEM method.
fix.loc <- c(1:5, 53:55)
(mod.fix2 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
                     gprior=list(dist="beta", params=c(5, 16)), EmpHist=FALSE,
                     Etol=1e-3, fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
(prior.par <- mod.fix2$group.par)
(emphist <- getirt(mod.fix2, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# fit the 3PL model to all dichotomous items, fit the GRM model to all polytomous data,
# at this moment fix only the five 3PL items (1st - 5th items)
# and estimate the empirical histogram of latent variable.
# use the OEM method. Thus, only 1 EM cycle is used.
fix.loc <- c(1:5)
(mod.fix3 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
                     gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE,
                     Etol=1e-3, fipc=TRUE, fipc.method="OEM", fix.loc=fix.loc))
(prior.par <- mod.fix3$group.par)
(emphist <- getirt(mod.fix3, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# summary of the estimation
summary(mod.fix3)

# fit the 3PL model to all dichotomous items, fit the GRM model to all polytomous data,
# at this moment fix all 55 items and estimate only the latent ability distribution 
# using the MEM method.
fix.loc <- c(1:55)
(mod.fix4 <- est_irt(x=x, data=sim.dat2, D=1, EmpHist=TRUE,
                     Etol=1e-3, fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
(prior.par <- mod.fix4$group.par)
(emphist <- getirt(mod.fix4, what="weights"))
plot(emphist$weight ~ emphist$theta, type="h")

# summary of the estimation
summary(mod.fix4)

}

}
\references{
Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.

Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm.
\emph{Psychometrika, 46}, 443-459.

Kim, S. (2006). A comparative study of IRT fixed parameter calibration methods.
\emph{Journal of Educational Measurement, 43}(4), 355-381.

Meilijson, I. (1989). A fast improvement to the EM algorithm on its own terms.
\emph{Journal of the Royal Statistical Society: Series B (Methodological), 51}, 127-138.

Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.

Wainer, H., & Mislevy, R. J. (1990). Item response theory, item calibration, and proficiency estimation. In H. Wainer (Ed.),
\emph{Computer adaptive testing: A primer} (Chap. 4, pp.65-102). Hillsdale, NJ: Lawrence Erlbaum.

Woods, C. M. (2007). Empirical histograms in item response theory with ordinal data. \emph{Educational and Psychological Measurement, 67}(1), 73-87.
}
\seealso{
\code{\link{est_item}}, \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{sx2_fit}},
\code{\link{traceline.est_item}}, \code{\link{getirt}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trace.R
\name{traceline}
\alias{traceline}
\alias{traceline.default}
\alias{traceline.est_item}
\alias{traceline.est_irt}
\title{Compute Item/Test Characteristic Functions}
\usage{
traceline(x, ...)

\method{traceline}{default}(x, theta, D = 1, ...)

\method{traceline}{est_item}(x, theta, ...)

\method{traceline}{est_irt}(x, theta, ...)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object
of class \code{\link{est_item}} obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}}
obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}}
for more details about the item metadata. The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}.}

\item{...}{Further arguments passed to or from other methods.}

\item{theta}{A vector of theta values.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}
}
\value{
This function returns an object of class \code{\link{traceline}}. This object contains a list containing
the item category probabilities, item characteristic function, and test characteristic function.
}
\description{
This function computes the item category probabilities, item characteristic function, and
test characteristic function given a set of theta values. The returned object of this function can be used
to draw the item or test characteristic curve using the function \code{\link{plot.traceline}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Default method to compute the item category probabilities, item characteristic function, and
test characteristic function for a data frame \code{x} containing the item metadata.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.

\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.
}}

\examples{
## example
## using a "-prm.txt" file obtained from a flexMIRT
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# set theta values
theta <- seq(-3, 3, 0.5)

# compute the item category probabilities and item/test
# characteristic functions given the theta values
traceline(x=test_flex, theta, D=1)

}
\seealso{
\code{\link{plot.traceline}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bring_output.R
\name{bring.flexmirt}
\alias{bring.flexmirt}
\alias{bring.bilog}
\alias{bring.parscale}
\alias{bring.mirt}
\title{Import Item and Ability Parameters from IRT Software}
\usage{
bring.flexmirt(
  file,
  type = c("par", "sco"),
  rePrm = TRUE,
  rePrm.gpc = TRUE,
  n.factor = 1
)

bring.bilog(file, type = c("par", "sco"))

bring.parscale(file, type = c("par", "sco"))

bring.mirt(x)

bring.bilog(file, type = c("par", "sco"))

bring.parscale(file, type = c("par", "sco"))

bring.mirt(x)
}
\arguments{
\item{file}{A file name (including a directory) containing the item or ability parameters.}

\item{type}{A character string indicating a type of output file. Available types are "par" for a file
containing item parameter estimates and "sco" for a file containing ability parameter estimates.}

\item{rePrm}{A logical value. If TRUE and when the IRT dichotomous model (e.g., 3PLM) or GRM is fit to data, 
the item intercept and logit of item guessing parameters are reparameterized into the item difficulty 
and item guessing parameters, respectively. Default is TRUE.}

\item{rePrm.gpc}{A logical value. If TRUE and when (G)PCM is fit to data, the nominal model
parameters in the flexMIRT parameter output file are reparameterized into the (G)PCM slope/difficulty parameters. 
Default is TRUE.}

\item{n.factor}{A numeric value indicating the number of estimated factors. This argument should be specified
when \code{type = "sco"}. Default is 1.}

\item{x}{An output object obtained from the function \code{\link[mirt]{mirt}}.}
}
\value{
These functions return a list including several objects. Only for the output of flexMIRT, the results of
multiple group analysis can be returned. In that case, each element of the list contains the estimation results for
each group.
}
\description{
These functions import item and/or ability parameters from BILOG-MG 3, PARSCALE 4, flexMIRT, and
mirt (R package).
}
\details{
The \code{\link{bring.flexmirt}} was written by modifying the function \code{read.flexmirt}
(Pritikin, 2018). The functions \code{\link{bring.bilog}} and \code{\link{bring.parscale}}
were written by modifying the functions \code{read.bilog} and \code{read.parscale}
(Weeks, 2017), respectively.

The file extensions for item parameter and ability files, respectively, are: ".par" and ".sco"
for BILOG-MG and PARSCALE, and "-prm.txt" and "-sco.txt" for flexMIRT. For mirt, the name of the output
object is specified by the user.

Although \code{\link{bring.flexmirt}} is able to extract multidimensional item and ability parameter estimates,
this package only deals with unidimensional IRT methods.

For polytomous item parameters, \code{\link{bring.flexmirt}} and \code{\link{bring.mirt}} are able to import
the item parameters of the graded response model and the (generalized) partial credit model.
}
\note{
Regarding the item parameter files for any IRT software, only the internal object "full_df" in the returned list is
necessary for the IRT linking. The object "full_df" is a data frame containing the item metadata
in a test form (e.g., item parameters, number of categories, models). See \code{\link{test.info}}
or \code{\link{simdat}} for more details about the item metadata.

Also, when item parameters are estimated using the partial credit or the generalized partial credit model,
item step parameters are returned in the object "full_df". Item step parameters are the overall item difficulty (or location)
parameter subtracted by the difficulty (or threshold) parameter for each category. See \code{\link{irtfit}} for more details
about the parameterization of the (generalized) partial credit model.
}
\section{Sample Output Files of IRT software}{


To illustrate how to import the item parameter estimate files of PARSCALE 4 and flexMIRT
using \code{\link{bring.parscale}} and \code{\link{bring.flexmirt}}, two item parameter
estimate output files are included in this package.

Among the two output files, one of them is from PARSCALE 4 with a file extension of ".PAR"
(i.e., "parscale_sample.PAR") and another one is from flexMIRT
with a file extension of "-prm.txt" (i.e., "flexmirt_sample-prm.txt").

For the two item parameter estimate output files, both are mixed-format tests with 55 items
consisting of fifty dichotomous items following the IRT 3PL model and five polytomous items with five
categories following the graded response model. The examples below show how to import those output files.
}

\examples{
## example 1
# import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item meta data
bring.flexmirt(file=flex_sam, "par")$Group1$full_df

## example 2
## import the ".par" output file from PARSCALE
pscale_sam <- system.file("extdata", "parscale_sample.PAR", package = "irtplay")

# read item parameters and transform them to item meta data
bring.parscale(file=pscale_sam, "par")$full_df

}
\references{
Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
Chapel Hill, NC: Vector Psychometric Group.

Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment.
\emph{Journal of Statistical Software, 48}(6), 1-29.

Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests Using IRT-Based Methods.
\emph{Journal of Statistical Software, 35}(12), 1-33. URL http://www.jstatsoft.org/v35/i12/.

Pritikin, J. (2018). \emph{rpf: Response Probability Functions}. R package version 0.59.
https://CRAN.R-project.org/package=rpf

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com

Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). BILOG-MG 3: Multiple-group
IRT analysis and test maintenance for binary items [Computer Program]. Chicago, IL: Scientific
Software International. URL http://www.ssicentral.com
}
\seealso{
\code{\link{irtfit}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdat.R
\name{simdat}
\alias{simdat}
\title{Simulated Response Data}
\usage{
simdat(
  x = NULL,
  theta,
  a.dc,
  b.dc,
  g.dc = NULL,
  a.py,
  d.py,
  cats,
  pmodel,
  D = 1
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...). This data frame
can be easily obtained using the function \code{\link{shape_df}}. See below for details.}

\item{theta}{A vector of theta values.}

\item{a.dc}{A vector of item discrimination (or slope) parameters for dichotomous IRT models.}

\item{b.dc}{A vector of item difficulty (or threshold) parameters for dichotomous IRT models.}

\item{g.dc}{A vector of item guessing parameters for dichotomous IRT models.}

\item{a.py}{A vector of item discrimination (or slope) parameters for polytomous IRT models.}

\item{d.py}{A list containing vectors of item threshold (or step) parameters for polytomous IRT models.}

\item{cats}{A vector containing the number of score categories for items.}

\item{pmodel}{A vector of character strings specifying the polytomous model with which response data are simulated.
For each polytomous model, "GRM" for the graded response model or "GPCM" for the (generalized) partial credit model can be
specified.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}
}
\value{
This function returns a vector or a matrix. When a matrix is returned, rows indicate theta values and columns represent items.
}
\description{
This function generates a simulated response data for a single- or a mixed-format test forms. For dichotomous
item response data, the IRT 1PL, 2PL, and 3PL models are available. For polytomous item response data, the graded response model,
the partial credit model, and the generalized partial credit model are available.
}
\details{
There are two ways of generating the simulated response data.
The first way is by using the argument \code{x} to read in a data frame of item metadata. In the data frame, the first column should have item IDs,
the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
In the \pkg{irtplay} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as 
the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item. 
Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because 
the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories, 
four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
  ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
  ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
  ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
}
And an example of a data frame for a mixed-format test is as follows:
\tabular{lrlrrrrr}{
  ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
  ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
  ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
  ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
  ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
  ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
  ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
  ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
}
See \code{IRT Models} section in the page of \code{\link{irtplay-package}} for more details about the IRT models used in the \pkg{irtplay} package. 
An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.

The second way is by directly specifying item parameters for each item for which response data should be simulated
(i.e., without using a data frame, as shown in the examples that follow). In addition to item parameters,
\code{theta}, \code{cats}, \code{pmodel}, and  \code{D} should be specified as well. \code{g.dc} does not need to be specified when only
the 1PL and 2PL models are used for dichotomous item response data. For dichotomous items, 2s should be specified in \code{cats}.
For polytomous items, the number of unique score categories should be specified in \code{cats}. When a response data set is generated with
a mixed-format test, it is important to clearly specify \code{cats} according to the order of items in the test form. Suppose that the response
data of ten examinees are simulated with five items, including three dichotomous items and two polytomous items with three categories.
Also, suppose that the second and the forth items are the polytomous items. Then, \code{cats = c(2, 3, 2, 3, 2)} should be used.
Additionally, among those two polytomous items, if the first and second item response data are simulated from the graded response model
and generalized partial credit model, respectively, then \code{pmodel = c('GRM', 'GPCM')}.
}
\examples{
## example 1.
## simulates response data with a mixed-format test.
## for the first two polytomous items, the generalized partial credit model is used
## for the last polytomous item, the graded response model is used
# 100 examinees are sampled
theta <- rnorm(100)

# set item parameters for three dichotomous items with the 3PL model
a.dc <- c(1, 1.2, 1.3); b.dc <- c(-1, 0, 1); g.dc <- rep(0.2, 3)

# set item parameters for three polytomous item parameters
# note that 4, 4, and 5 categories are used for polytomous items
a.py <- c(1.3, 1.2, 1.7)
d.py <- list(c(-1.2, -0.3, 0.4), c(-0.2, 0.5, 1.6), c(-1.7, 0.2, 1.1, 2.0))

# create a numeric vector of score categoires for both dichotomous and polytomous item data
# this score category vector is used to specify the location of the polytomous items
cats <- c(2, 2, 4, 4, 5, 2)

# create a character vector of the IRT model for the polytomous items
pmodel <- c('GPCM', 'GPCM', 'GRM')

# simulate the response data
simdat(theta=theta, a.dc=a.dc, b.dc=b.dc, g.dc=NULL,
       a.py=a.py, d.py=d.py, cats=cats, pmodel=pmodel, D=1)


## example 2.
## simulates response data with a sigle-format test with the 2PL model.
# create a numeric vector of score categoires for the three 2PL model items
cats <- rep(2, 3)

# simulate the response data
simdat(theta=theta, a.dc=a.dc, b.dc=b.dc, cats=cats, D=1)

## example 3.
## the use of a "-prm.txt" file obtained from a flexMIRT
# import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# read item parameters and transform them to item metadata
test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df

# simulate the response data
simdat(x=test_flex, theta=theta, D=1) # use a data.farame of item meta information

}
\seealso{
\code{\link{drm}}, \code{\link{plm}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_flexmirt.R
\name{run_flexmirt}
\alias{run_flexmirt}
\title{Run flexMIRT through R}
\usage{
run_flexmirt(file.syntax, dir.flex = NULL, show.output.on.console = FALSE, ...)
}
\arguments{
\item{file.syntax}{A single string or vector containing the file path(s) of a flexmirt syntax file(s) to be run.
An example is “C:/Users/Data/irtmodel.flexmirt".}

\item{dir.flex}{A path of directory where flexMIRT is installed. The path may include a folder name with "flexMIRT" 
(e.g, flexMIRT3, flexMIRT 3.6). If NULL, a path where flexMIRT is installed will be searched in "C:/Program Files" and 
it will be used as a default path (e.g., "C:/Program Files/flexMIRT3", "C:/Program Files/flexMIRT 3.6").}

\item{show.output.on.console}{A logical value to indicate whether to capture the output of the command and show it on the R console.
Default is FALSE. See \code{\link[base]{system}}.}

\item{...}{Further arguments passed from the function \code{\link[base]{system}}.}
}
\value{
output files of flexMIRT
}
\description{
This function implements flexMIRT (Cai, 2017) to run a model specified in the syntax file of
flexMIRT (i.e., *.flexmirt) through R. To run this function, flexMIRT software must be installed in advance.
This function will be useful especially when conducting a simulation study using flexMIRT.
}
\details{
When a path of directory where flexMIRT (with a version < 3.6) is installed is provided 
in the argument \code{dir.flex}, the directory must include following six file of
\itemize{
  \item WinFlexMIRT.exe
  \item FlexMIRT_x64.exe
  \item FlexMIRT_x86.exe
  \item vpg.dll
  \item vpg.licensing.client.dll
  \item vpg.licensing.dll
}
When a path of directory where flexMIRT (with a version >= 3.6) is installed is provided 
in the argument \code{dir.flex}, the directory must include following six files of
\itemize{
  \item WinFlexMIRT.exe
  \item vpg.dll
  \item vpg.licensing.client.dll
  \item vpg.licensing.dll
  \item VPGLicenseClientNet.dll
}
and an additional directory of "Resources" that contains two files which are
\itemize{
  \item flexMIRT_x64_AVX.exe
  \item flexMIRT_x86_AVX.exe
}
}
\examples{

# Emxaples below will run when flexMIRT software is installed
# in a default path of "C:/Program Files/flexMIRT3".
# Otherwise provide a path where flexMIRT software is installed
# in the argument 'dir.flex'.

\dontrun{
# (1) run a single syntax file
# import an example of flexMIRT syntax file to run the item parameter estimation of IRT 3PL model
file.syntax <- system.file("extdata", "2PLM_example.flexmirt", package = "irtplay")

# run flexMIRT to estimate the item parameters of IRT 3PL model
run_flexmirt(file.syntax=file.syntax, dir.flex=NULL, show.output=TRUE)

# check the output file
out.file <- system.file("extdata", "2PLM_example-prm.txt", package = "irtplay")
bring.flexmirt(out.file, type="par")

# (2) run multiple syntax files
# import two examples of flexMIRT syntax files
file.syntax1 <- system.file("extdata", "2PLM_example.flexmirt", package = "irtplay")
file.syntax2 <- system.file("extdata", "3PLM_example.flexmirt", package = "irtplay")

# run flexMIRT to estimate the item parameters
run_flexmirt(file.syntax=c(file.syntax1, file.syntax2), dir.flex=NULL, show.output=FALSE)

# check the output file
out.file1 <- system.file("extdata", "2PLM_example-prm.txt", package = "irtplay")
out.file2 <- system.file("extdata", "3PLM_example-prm.txt", package = "irtplay")
bring.flexmirt(out.file1, type="par")
bring.flexmirt(out.file2, type="par")
}


}
\references{
Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
Chapel Hill, NC: Vector Psychometric Group.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_irtfit.R
\name{plot.irtfit}
\alias{plot.irtfit}
\title{Draw raw and standardized residual plots}
\usage{
\method{plot}{irtfit}(
  x,
  item.loc = NULL,
  type = "both",
  ci.method = c("wald", "cp", "wilson", "wilson.cr"),
  show.table = TRUE,
  layout.col = 2,
  xlab.text,
  ylab.text,
  main.text,
  lab.size = 15,
  main.size = 15,
  axis.size = 15,
  line.size = 1,
  point.size = 2.5,
  strip.size = 12,
  ylim.icc = c(0, 1),
  ylim.sr.adjust = FALSE,
  ylim.sr = c(-4, 4),
  ...
)
}
\arguments{
\item{x}{An object of class \code{\link{irtfit}}.}

\item{item.loc}{An integer value indicating that the \emph{n}th item (or the location of the item) is plotted. See below for
details.}

\item{type}{A character string indicating what type of residual plot is returned. Available options
are "icc" for the raw residual plot, "sr" for the standardized residual plot, and "both" for both of them.
Default is "both".}

\item{ci.method}{A character string indicating what method is used to estimate the confidence interval for the raw residual plot.
Available options are "wald" for Wald method, "cp" for Clopper-Pearson interval, "wilson" for Wilson score interval, and
"wilson.cr" for Wilson score interval with continuity correction. Default is "wald". See below for details.}

\item{show.table}{A logical value. If TRUE, a contingency table containing the information used to draw the residual
plots for the studied item is returned. This contingency table is the same as one contained in the internal object of \code{contingency.plot}
in the object of class \code{\link{irtfit}}. Default is TRUE.}

\item{layout.col}{An integer value indicating the number of columns in the panel when a polytomous item is used.
Default is 2.}

\item{xlab.text}{A title for the x axis. If missing, the default string is used.}

\item{ylab.text}{A title for the y axis. If \code{type = "both"}, two character strings can be
specified for the raw residual and standardized residual plots, respectively. If missing,
the default strings are used.}

\item{main.text}{An overall title for the plot. If \code{type = "both"}, two character strings
can be specified for the raw residual and standardized residual plots, respectively. If missing,
the default strings are used.}

\item{lab.size}{The size of xlab and ylab. Default is 15.}

\item{main.size}{The size of \code{main.text}. Default is 15.}

\item{axis.size}{The size of labels along the x and y axes. Default is 15.}

\item{line.size}{The size of lines. Default is 1.}

\item{point.size}{The size of points. Default is 2.5.}

\item{strip.size}{The size of facet labels. Default is 12.}

\item{ylim.icc}{A vector of two numeric values specifying the range of y axis for the raw residual plot. Default is c(0, 1).}

\item{ylim.sr.adjust}{A logical value. If TRUE, the range of y axis for the standardized residual plot is adjusted for each item.
If FALSE, the range of y axis for the standardized residual plot is fixed to the values specified in the argument \code{ylim.sr}.}

\item{ylim.sr}{A vector of two numeric values specifying the range of y axis for the standardized residual plot.
Default is c(-4, 4).}

\item{...}{Further arguments passed from the function \code{ggplot()} in the \pkg{ggplot2} package.}
}
\description{
This function provides graphical displays to look at residuals between the observed data
and model-based predictions (Hambleton, Swaminathan, & Rogers, 1991). This function gives two residual plots for
each score category of an item: (a) the raw residual plot and (b) the standardized residual plot. Note that
for dichotomous items the residual plots are drawn only for the score category of 1.
}
\details{
All of the plots are drawn using the ggplot2 package.

Once the results of the IRT model fit analysis are obtained from the function \code{\link{irtfit}},
an object of class \code{\link{irtfit}} can be used to draw the IRT raw residual and standardized residual plots. Especially, the information
contained in an internal object of \code{contingency.plot} are mainly used to draw the residual plots.

Because the residual plots are drawn for an item at a time, you have to indicate which item will be evaluated. For this,
you should specify an integer value, which is the location of the studied item, in the argument \code{item.loc}.
For example, if you want to draw the residual plots for the third item, then \code{item.loc = 3}.

In terms of the raw residual plot, the argument \code{ci.method} is used to select a method to estimate the confidence intervals
among four methods. Those methods are "wald" for the Wald interval, which is based on the normal approximation (Laplace, 1812),
"cp" for Clopper-Pearson interval (Clopper & Pearson, 1934), "wilson" for Wilson score interval (Wilson, 1927), and
"wilson.cr" for Wilson score interval with continuity correction (Newcombe, 1998).
See \url{https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval} for more details about
the binomial proportion confidence intervals. Note that the width of confidence interval is determined by the \eqn{\alpha}-level
specified in the argument \code{alpha} of the function \code{\link{irtfit}}.

Regarding the standardized residual plot, any standardized residuals greater than the specified criterion value
in the argument {\code{overSR}} of the function \code{\link{irtfit}} are displayed as triangles. Otherwise,
they are displayed as circles.
}
\examples{
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the first two dichotomous items and last polytomous item
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]

# generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(1000, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

\donttest{
# compute fit statistics
fit <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
               n.width=11, loc.theta="average", range.score=c(-4, 4), D=1, alpha=0.05, overSR=1.5)

# residual plots for the first item (dichotomous item)
plot(x=fit, item.loc=1, type = "both", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)

# residual plots for the third item (polytomous item)
plot(x=fit, item.loc=3, type = "both", ci.method = "wald", show.table=FALSE, ylim.sr.adjust=TRUE)

# raw residual plot for the third item (polytomous item)
plot(x=fit, item.loc=3, type = "icc", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)

# standardized residual plot for the third item (polytomous item)
plot(x=fit, item.loc=3, type = "sr", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)
}

}
\references{
Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial.
\emph{Biometrika, 26}(4), 404-413.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.

Laplace, P. S. (1820).\emph{Theorie analytique des probabilites} (in French). Courcier.

Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of seven methods.
\emph{Statistics in medicine, 17}(8), 857-872.

Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
\emph{Journal of the American Statistical Association, 22}(158), 209-212.
}
\seealso{
\code{\link{irtfit}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglike_item.R
\name{llike_item}
\alias{llike_item}
\title{Loglikelihood of Items}
\usage{
llike_item(
  x,
  data,
  score,
  D = 1,
  use.aprior = FALSE,
  use.bprior = FALSE,
  use.gprior = FALSE,
  aprior = list(dist = "lnorm", params = c(0, 0.5)),
  bprior = list(dist = "norm", params = c(0, 1)),
  gprior = list(dist = "beta", params = c(5, 17)),
  missing = NA
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
See \code{\link{irtfit}}, \code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}. If \code{prob = NULL}, this data frame is
used in the recursion formula. See below for details.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{score}{A vector of examinees' ability estimates. Length of the vector must be the same as the number of rows in the
response data set.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{use.aprior}{A logical value. If TRUE, a prior distribution for the slope parameters is used when computing the loglikelihood values
across all items. Default is FALSE.}

\item{use.bprior}{A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used when computing the loglikelihood values
across all items. Default is FALSE.}

\item{use.gprior}{A logical value. If TRUE, a prior distribution for the guessing parameters is used when computing the loglikelihood values
across all 3PLM items. Default is TRUE.}

\item{aprior}{A list containing the information of the prior distribution for item slope parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()}, and \code{dnorm()}
in the \pkg{stats} package for more details.}

\item{bprior}{A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()}, and \code{dnorm()}
in the \pkg{stats} package for more details.}

\item{gprior}{A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()}, and \code{dnorm()}
in the \pkg{stats} package for more details.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}
}
\value{
A vector of loglikelihood values. Each element represents a sum of loglikeihoods across all ability values for each item.
}
\description{
This function computes the loglikelihoods of individual items given the item parameters, ability values, and response data.
}
\examples{
## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the first two dichotomous items and last polytomous item
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]

# generate examinees' abilities from N(0, 1)
set.seed(10)
score <- rnorm(10, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

# compute the loglikelihood values (no priors are used)
llike_item(x, data, score, D=1, use.aprior=FALSE, use.gprior=FALSE)

}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_flexmirt.R
\name{write.flexmirt}
\alias{write.flexmirt}
\title{Write a "-prm.txt" file for flexMIRT}
\usage{
write.flexmirt(x, file = NULL, norm.pop = c(0, 1), rePrm = TRUE)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}.}

\item{file}{The destination file name.}

\item{norm.pop}{A numeric vector of two components specifying a mean and standard deviation of the normal
population distribution. Default is c(0,1).}

\item{rePrm}{A logical value indicating whether the item parameters in the item metadata
are the reparameterized item parameters. If TRUE, the item intercepts and logits of item guessing parameters
should be included in the item metadata. If FALSE, the item difficulty and item guessing parameters
should be included in the item metadata.}
}
\value{
A "-prm.txt" file.
}
\description{
This function writes an output file of "-prm.txt" for flexMIRT (Cai, 2017). The current version of this function
can be used only for the unidimensional IRT models.
}
\examples{
## use the simulated CAT data
# extract the item metadata
x <- simCAT_MX$item.prm

# set a name of "-prm.txt" file
temp_prm <- file.path(tempdir(), "temp-prm.txt")

# write out the "-prm.txt" file
write.flexmirt(x, file=temp_prm, norm.pop=c(0, 1), rePrm=FALSE)

}
\references{
Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
Chapel Hill, NC: Vector Psychometric Group.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_den.R
\name{post_den}
\alias{post_den}
\title{Updated prior (a.k.a. posterior) latent ability distribution}
\usage{
post_den(
  x,
  data,
  D = 1,
  Quadrature = c(49, 6),
  weights = NULL,
  group.mean = 0,
  group.var = 1,
  missing = NA
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...). 
See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}}  for more details about the item metadata.
This data frame can be easily obtained using the function \code{\link{shape_df}}.}

\item{data}{A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
the examinees and items, respectively.}

\item{D}{A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
Default is 1.}

\item{Quadrature}{A numeric vector of two components specifying the number of quadrature points (in the first component) and
the symmetric minimum and maximum values of these points (in the second component). For example, a vector of c(49, 6) indicates 49 rectangular
quadrature points over -6 and 6. Default is c(49, 6).}

\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
using the function \code{\link{gen.weight}}. If NULL, a normal prior density is used based on the information provided in the arguments
of \code{Quadrature}, \code{group.mean}, and \code{group.var}). Default is NULL.}

\item{group.mean}{A numeric value to set the mean of latent variable prior distribution. Default is 0.}

\item{group.var}{A positive numeric value to set the variance of latent variable prior distribution. Default is 1.}

\item{missing}{A value indicating missing values in the response data set. Default is NA.}
}
\value{
This function returns a list containing two internal objects. The first internal object is a data frame with two columns, 
where the first column has theta values (nodes) and the second column provides the weights of the posterior latent ability distribution. 
The second internal object is a data frame containing the mean, variance, and standard deviation of the distribution.
}
\description{
This function computes updated prior (a.k.a. posterior) densities of the latent ability distribution given 
a prior ability distribution, item parameters, and item response data.
}
\examples{

\donttest{
# fit the 2PL model to LSAT6 data
(mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2))

# extract the item parameter estimates
(x <- getirt(x=mod.2pl, what="par.est"))

# update the standard normal prior deisnty of the ability distribution
# using the estimated item parameters
(upd_prior <- post_den(x=x, data=LSAT6, D=1, group.mean=0, group.var=1))
}

}
\seealso{
\code{\link{shape_df}}, \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtplay-package.R
\docType{package}
\name{irtplay-package}
\alias{irtplay-package}
\title{irtplay: Unidimensional Item Response Theory Modeling}
\description{
Fit unidimensional item response theory (IRT) models to a mixture of dichotomous and polytomous data,
calibrate online item parameters (i.e., pretest and operational items), estimate examinees' abilities, 
and examine the IRT model-data fit on item-level in different ways as well as provide useful functions 
related to unidimensional IRT.

For the item parameter estimation, the marginal maximum likelihood estimation via expectation-maximization (MMLE-EM) algorithm
(Bock & Aitkin, 1981) is used. For the online calibration, the fixed item parameter calibration (FIPC) method (Kim, 2006) and 
the fixed ability parameter calibration (FAPC) method, (Ban, Hanson, Wang, Yi, & Harris, 2001; stocking, 1988),
often called Stocking's Method A, are provided. For the ability estimation, several popular scoring methods (e.g., MLE, EAP, and MAP) 
are implemented. In terms of assessing the IRT model-data fit, one of distinguished features of this package is that it gives 
not only item fit statistics (e.g., \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981), likelihood ratio \eqn{\chi^{2}} 
fit statistic (\eqn{G^{2}}; McKinley & Mills, 1985), infit and outfit statistics (Ames et al., 2015), and \eqn{S-X^{2}} 
(Orlando & Thissen, 2000, 2003)) but also graphical displays to look at residuals between the observed data and model-based 
predictions (Hambleton, Swaminathan, & Rogers, 1991).

In addition, there are many useful functions such as analyzing differential item functioning (DIF), computing asymptotic 
variance-covariance matrices of item parameter estimates,importing item and/or ability parameters from popular IRT software, 
running flexMIRT (Cai, 2017) through R, generating simulated data, computing the conditional distribution of observed scores 
using the Lord-Wingersky recursion formula, computing the loglikelihood of individual items, computing the loglikelihood 
of abilities, computing item and test information functions, computing item and test characteristic curve functions, and 
plotting item and test characteristic curves and item and test information functions.

\tabular{ll}{ Package: \tab irtplay\cr Version: \tab 1.6.3\cr Date: \tab
2021-11-04\cr Depends: \tab R (>= 3.6)\cr License: \tab GPL (>= 2)\cr }
}
\details{
Following five sections describe a) how to implement the online item calibration using FIPC, a) how to implement the online item
calibration using Method A, b) the process of evaluating the IRT model-data fit, c) two examples for the online calibration and
evaluating the IRT model-data fit, and d) IRT Models used in \pkg{irtplay} package.
}
\section{Online item calibration with the fixed item parameter calibration method (e.g., Kim, 2006)}{


The fixed item parameter calibration (FIPC) is one of useful online item calibration methods for computerized adaptive testing (CAT)
to put the parameter estimates of pretest items on the same scale of operational item parameter estimates without post hoc
linking/scaling (Ban, Hanson, Wang, Yi, & Harris, 2001; Chen & Wang, 2016). In FIPC, the operational item parameters are fixed to
estimate the characteristic of the underlying latent variable prior distribution when calibrating the pretest items. More specifically,
the underlying latent variable prior distribution of the operational items is estimated during the calibration of the pretest
items to put the item parameters of the pretest items on the scale of the operational item parameters (Kim, 2006). In the \pkg{irtplay}
package, FIPC is implemented with two main steps:

\enumerate{
  \item Prepare a response data set and the item metadata of the fixed (or operational) items.
  \item Implement FIPC to estimate the item parameters of pretest items using the \code{\link{est_irt}} function.
}

\describe{
  \item{1. Preparing a data set}{
  To run the \code{\link{est_irt}} function, it requires two data sets:

    \enumerate{
      \item Item metadata set (i.e., model, score category, and item parameters. see the desciption of the argument \code{x} in the function \code{\link{est_irt}}).
      \item Examinees' response data set for the items. It should be a matrix format where a row and column indicate the examinees and the items, respectively.
      The order of the columns in the response data set must be exactly the same as the order of rows of the item metadata.
    }
  }

  \item{2. Estimating the pretest item parameters}{
  When FIPC is implemented in \code{\link{est_irt}} function, the pretest item parameters are estimated by fixing the operational item parameters. To estimate the item
  parameters, you need to provide the item metadata in the argument \code{x} and the response data in the argument \code{data}.

  It is worthwhile to explain about how to prepare the item metadata set in the argument \code{x}. A specific form of a data frame should be used for
  the argument \code{x}. The first column should have item IDs, the second column should contain the number of score categories of the items, and the third
  column should include IRT models. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous
  items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded response model and
  (generalized) partial credit model, respectively. From the fourth column, item parameters should be included. For dichotomous items, the fourth, fifth,
  and sixth columns represent the item discrimination (or slope), item difficulty, and item guessing parameters, respectively. When "1PLM" or "2PLM" is
  specified for any items in the third column, NAs should be inserted for the item guessing parameters. For polytomous items, the item discrimination (or slope)
  parameters should be contained in the fourth column and the item threshold (or step) parameters should be included from the fifth to the last columns.
  When the number of categories differs between items, the empty cells of item parameters should be filled with NAs. See `est_irt` for more details about
  the item metadata.

  Also, you should specify in the argument \code{fipc = TRUE} and a specific FIPC method in the argument \code{fipc.method}. Finally, you should provide
  a vector of the location of the items to be fixed in the argument \code{fix.loc}. For more details about implementing FIPC, see the
  description of the function \code{\link{est_irt}}.

  When implementing FIPC, you can estimate both the emprical histogram and the scale of latent variable prior distribution by setting \code{EmpHist = TRUE}.
  If \code{EmpHist = FALSE}, the normal prior distribution is used during the item parameter estimation and the scale of the normal prior distribution is
  updated during the EM cycle.

  The \code{\link{est_item}} function requires a vector of the number of score categories for the items in the argument \code{cats}. For example, a dichotomous item has
  two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items are binary items
  (i.e., dichotomous items), it assumes that all items have two score categories.

  If necessary, you need to specify whether prior distributions of item slope and guessing parameters (only for the IRT 3PL model) are used in the arguments of
  \code{use.aprior} and \code{use.gprior}, respectively. If you decide to use the prior distributions, you should specify what distributions will be used for the prior
  distributions in the arguments of \code{aprior} and \code{gprior}, respectively. Currently three probability distributions of Beta, Log-normal, and Normal
  distributions are available.

  In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.

  Once the \code{\link{est_irt}} function has been implemented, you'll get a list of several internal objects such as the item parameter estimates,
  standard error of the parameter estimates.
  }
}
}

\section{Online item calibration with the fixed ability parameter calibration method (e.g., Stocking, 1988)}{

In CAT, the fixed ability parameter calibration (FAPC), often called Stocking's Method A, is the relatively simplest 
and most straightforward online calibration method, which is the maximum likelihood estimation of the item parameters 
given the proficiency estimates. In CAT, FAPC can be used to put the parameter estimates of pretest items on 
the same scale of operational item parameter estimates and recalibrate the operational items to evaluate the parameter 
drifts of the operational items (Chen & Wang, 2016; Stocking, 1988). Also, FAPC is known to result in accurate, unbiased 
item parameters calibration when items are randomly rather than adaptively administered to examinees, which occurs most 
commonly with pretest items (Ban, Hanson, Wang, Yi, & Harris, 2001; Chen & Wang, 2016). Using \pkg{irtplay} package, 
the FAPC is implemented to calibrate the items with two main steps:

\enumerate{
  \item Prepare a data set for the calibration of item parameters (i.e., item response data and ability estimates).
  \item Implement the FAPC to estimate the item parameters using the \code{\link{est_item}} function.
}

\describe{
  \item{1. Preparing a data set}{
  To run the \code{\link{est_item}} function, it requires two data sets:

    \enumerate{
      \item Examinees' ability (or proficiency) estimates. It should be in the format of a numeric vector.
      \item response data set for the items. It should be in the format of matrix where a row and column indicate
    the examinees and the items, respectively. The order of the examinees in the response data set must be exactly the same as that of the examinees' ability estimates.
    }
  }

  \item{2. Estimating the pretest item parameters}{
  The \code{\link{est_item}} function estimates the pretest item parameters given the proficiency estimates. To estimate the item parameters,
  you need to provide the response data in the argument \code{data} and the ability estimates in the argument \code{score}.

  Also, you should provide a string vector of the IRT models in the argument \code{model} to indicate what IRT model is used to calibrate each item.
  Available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent
  the graded response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function. If a single
  character of the IRT model is specified, that model will be recycled across all items.

  The \code{\link{est_item}} function requires a vector of the number of score categories for the items in the argument \code{cats}. For example, a dichotomous item has
  two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items are binary items
  (i.e., dichotomous items), it assumes that all items have two score categories.

  If necessary, you need to specify whether prior distributions of item slope and guessing parameters (only for the IRT 3PL model) are used in the arguments of
  \code{use.aprior} and \code{use.gprior}, respectively. If you decide to use the prior distributions, you should specify what distributions will be used for the prior
  distributions in the arguments of \code{aprior} and \code{gprior}, respectively. Currently three probability distributions of Beta, Log-normal, and Normal
  distributions are available.

  In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.

  Once the \code{\link{est_item}} function has been implemented, you'll get a list of several internal objects such as the item parameter estimates,
  standard error of the parameter estimates.
  }
}
}

\section{The process of evaluating the IRT model-data fit}{

One way to assess goodness of IRT model-data fit is through an item fit analysis by examining the traditional item fit statistics
and looking at the discrepancy between the observed data and model-based predictions. Using \pkg{irtplay} package, the traditional approach
of evaluating the IRT model-data fit on item-level can be implemented with three main steps:

\enumerate{
  \item Prepare a data set for the IRT item fit analysis (i.e., item metadata, ability estimates, and response data).
  \item Obtain the IRT fit statistics such as \eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit statistics using the function \code{\link{irtfit}}.
  \item Based on the results of IRT model fit analysis (i.e., an object of class \code{\link{irtfit}}) obtained in step 2,
draw the IRT residual plots (i.e., raw residual and standardized residual plots) using the function \code{\link{plot.irtfit}}.
}

\describe{
  \item{1. Preparing a data set}{
  Before conducting the IRT model fit analysis, it is necessary to prepare a data set. To run the function \code{\link{irtfit}}, it requires
  three data sets:

   \enumerate{
   \item Item metadata including the item ID, number of score categories, IRT models, and item parameters. The item metadata should be in the format of
   data frame. You can prepare the data either by using the function \code{\link{shape_df}} or by creating a data frame of the item metadata by yourself.
   If you have output files of item parameter estimates obtained from one of the IRT software such as BILOG-MG 3, PARSCALE 4, flexMIRT, and mirt (R package),
   the item metadata can be easily obtained using the functions of \code{\link{bring.bilog}}, \code{\link{bring.parscale}}, \code{\link{bring.flexmirt}},
   and \code{\link{bring.mirt}}. See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata format.
   \item Examinees' ability (or proficiency) estimates. It should be in the format of a numeric vector.
   \item Examinees' response data set for the items. It should be in the format of matrix where a row and column indicate the examinees and the items,
   respectively. The order of the examinees in the response data set must be exactly the same as that of the examinees' ability estimates. The order of the items
   in the response data set must be exactly the same as that of the items in the item metadata.
   }

 }

  \item{2. Computing the IRT model-data fit statistics}{
  The function \code{\link{irtfit}} computes the traditional IRT item fit statistics such as \eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit statistics.
  To calculate the \eqn{\chi^{2}} and \eqn{G^{2}} statistics, two methods are available to divide the ability scale into several groups. The two methods are "equal.width"
  for dividing the scale by an equal length of the interval and "equal.freq" for dividing the scale by an equal frequency of examinees. Also, you need to
  specify the location of ability point at each group (or interval) where the expected probabilities of score categories are calculated from the IRT models.
  Available locations are "average" for computing the expected probability at the average point of examinees' ability estimates in each group and "middle" for
  computing the expected probability at the midpoint of each group.

  To use the function \code{\link{irtfit}}, you need to insert the item metadata in the argument \code{x}, the ability estimates in the argument \code{score},
  and the response data in the argument \code{data}. If you want to divide the ability scale into other than ten groups, you need to specify the number of groups
  in the argument \code{n.width}. In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.

  Once the function \code{\link{irtfit}} has been implemented, you'll get the fit statistic results and the contingency tables for every item used
  to calculate the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
 }

  \item{3. Drawing the IRT residual plots}{
  Using the saved object of class \code{\link{irtfit}}, you can use the \code{\link{plot}} method to evaluate the IRT raw residual and standardized residual plots.

  Because the \code{\link{plot}} method can draw the residual plots for an item at a time, you have to indicate which item will be examined. For this,
  you can specify an integer value, which is the location of the studied item, in the argument \code{item.loc}.

  In terms of the raw residual plot, the argument \code{ci.method} is used to select a method to estimate the confidence intervals among four methods.
  Those methods are "wald" for the Wald interval, which is based on the normal approximation (Laplace, 1812), "cp" for Clopper-Pearson interval
  (Clopper & Pearson, 1934), "wilson" for Wilson score interval (Wilson, 1927), and "wilson.cr" for Wilson score interval with continuity correction
  (Newcombe, 1998).
 }
}
}

\section{Three examples of R script}{


The example code below shows how to implement the online calibration and how to evaluate the IRT model-data fit:\preformatted{
##---------------------------------------------------------------
# Attach the packages
library(irtplay)

##----------------------------------------------------------------------------
# 1. The example code below shows how to prepare the data sets and how to
#    implement the fixed item parameter calibration (FIPC):
##----------------------------------------------------------------------------

## Step 1: prepare a data set
## In this example, we generated examinees' true proficiency parameters and simulated
## the item response data using the function "simdat".

## import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# generate 1,000 examinees' latent abilities from N(0.4, 1.3)
set.seed(20)
score <- rnorm(1000, mean=0.4, sd=1.3)

# simulate the response data
sim.dat <- simdat(x=x, theta=score, D=1)

## Step 2: Estimate the item parameters
# fit the 3PL model to all dichotmous items, fit the GRM model to all polytomous data,
# fix the five 3PL items (1st - 5th items) and three GRM items (53th to 55th items)
# also, estimate the empirical histogram of latent variable
fix.loc <- c(1:5, 53:55)
(mod.fix1 <- est_irt(x=x, data=sim.dat, D=1, use.gprior=TRUE,
                    gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE, Etol=1e-3,
                    fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
summary(mod.fix1)

# plot the estimated empirical histogram of latent variable prior distribution
(emphist <- getirt(mod.fix1, what="weights"))
plot(emphist$weight ~ emphist$theta, xlab="Theta", ylab="Density")


##----------------------------------------------------------------------------
# 2. The example code below shows how to prepare the data sets and how to estimate
#    the item parameters using Method A:
##----------------------------------------------------------------------------

## Step 1: prepare a data set
## In this example, we generated examinees' true proficiency parameters and simulated
## the item response data using the function "simdat". Because, the true
## proficiency parameters are not known in reality, however, the true proficiencies
## would be replaced with the proficiency estimates for the calibration.

# import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")

# select the item metadata
x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df

# modify the item metadata so that some items follow 1PLM, 2PLM and GPCM
x[c(1:3, 5), 3] <- "1PLM"
x[c(1:3, 5), 4] <- 1
x[c(1:3, 5), 6] <- 0
x[c(4, 8:12), 3] <- "2PLM"
x[c(4, 8:12), 6] <- 0
x[54:55, 3] <- "GPCM"

# generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean=0, sd=1)

# simulate the response data
data <- simdat(x=x, theta=score, D=1)

## Step 2: Estimate the item parameters
# 1) item parameter estimation: constrain the slope parameters of the 1PLM to be equal
(mod1 <- est_item(x, data, score, D=1, fix.a.1pl=FALSE, use.gprior=TRUE,
                 gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
summary(mod1)

# 2) item parameter estimation: fix the slope parameters of the 1PLM to 1
(mod2 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, a.val.1pl=1, use.gprior=TRUE,
                 gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
summary(mod2)

# 3) item parameter estimation: fix the guessing parameters of the 3PLM to 0.2
(mod3 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, fix.g=TRUE, a.val.1pl=1, g.val=.2,
                 use.startval=FALSE))
summary(mod3)


##----------------------------------------------------------------------------
# 3. The example code below shows how to prepare the data sets and how to conduct
#    the IRT model-data fit analysis:
##----------------------------------------------------------------------------

## Step 1: prepare a data set for IRT
## In this example, we use the simulated mixed-item format of CAT Data
## But, only items that have examinees' responses more than 1,000 are assessed.

# find the location of items that have more than 1,000 item responses
over1000 <- which(colSums(simCAT_MX$res.dat, na.rm=TRUE) > 1000)

# (1) item metadata
x <- simCAT_MX$item.prm[over1000, ]

# (2) examinee's ability estimates
score <- simCAT_MX$score

# (3) response data
data <- simCAT_MX$res.dat[, over1000]

## Step 2: Compute the IRT mode-data fit statistics
# (1) the use of "equal.width"
fit1 <- irtfit(x=x, score=score, data=data, group.method="equal.width",
               n.width=10, loc.theta="average", range.score=NULL, D=1,
               alpha=0.05, missing=NA)

# what kinds of internal objects does the results have?
names(fit1)

# show the results of the fit statistics
fit1$fit_stat[1:10, ]

# show the contingency tables for the first item (dichotomous item)
fit1$contingency.fitstat[[1]]

# (2) the use of "equal.freq"
fit2 <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
               n.width=10, loc.theta="average", range.score=NULL, D=1,
               alpha=0.05, missing=NA)

# show the results of the fit statistics
fit2$fit_stat[1:10, ]

# show the contingency table for the fourth item (polytomous item)
fit2$contingency.fitstat[[4]]

## Step 3: Draw the IRT residual plots
# 1. for the dichotomous item
# (1) both raw and standardized residual plots using the object "fit1"
plot(x=fit1, item.loc=1, type = "both", ci.method = "wald",
     ylim.sr.adjust=TRUE)

# (2) the raw residual plots using the object "fit1"
plot(x=fit1, item.loc=1, type = "icc", ci.method = "wald",
     ylim.sr.adjust=TRUE)

# (3) the standardized residual plots using the object "fit1"
plot(x=fit1, item.loc=113, type = "sr", ci.method = "wald",
     ylim.sr.adjust=TRUE)

# 2. for the polytomous item
# (1) both raw and standardized residual plots using the object "fit1"
plot(x=fit1, item.loc=113, type = "both", ci.method = "wald",
     ylim.sr.adjust=TRUE)

# (2) the raw residual plots using the object "fit1"
plot(x=fit1, item.loc=113, type = "icc", ci.method = "wald",
     layout.col=2, ylim.sr.adjust=TRUE)

# (3) the standardized residual plots using the object "fit1"
plot(x=fit1, item.loc=113, type = "sr", ci.method = "wald",
     layout.col=4, ylim.sr.adjust=TRUE)
}
}

\section{IRT Models}{


In the \pkg{irtplay} package, both dichotomous and polytomous IRT models are available.
For dichotomous items, IRT one-, two-, and three-parameter logistic models (1PLM, 2PLM, and 3PLM) are used. 
For polytomous items, the graded response model (GRM) and the (generalized) partial credit model (GPCM) are used.
Note that the item discrimination (or slope) parameters should be fixed to 1 when the partial credit model is fit to data. 

In the following, let \eqn{Y} be the response of an examinee with latent ability \eqn{\theta} on an item and suppose that there 
are \eqn{K} unique score categories for each polytomous item. 

\describe{
  \item{IRT 1-3PL models}{
    For the IRT 1-3PL models, the probability that an examinee with \eqn{\theta} provides a correct answer for an item is given by,
     \deqn{P(Y = 1|\theta) = g + \frac{(1 - g)}{1 + exp(-Da(\theta - b))},}
    where \eqn{a} is the item discrimination (or slope) parameter, \eqn{b} represents the item difficulty parameter, 
    \eqn{g} refers to the item guessing parameter. \eqn{D} is a scaling factor in IRT models to make the logistic function 
    as close as possible to the normal ogive function when \eqn{D = 1.702}. When the 1PLM is used, \eqn{a} is either fixed to a constant 
    value (e.g., \eqn{a=1}) or constrained to have the same value across all 1PLM item data. When the IRT 1PLM or 2PLM is fit to data, 
    \eqn{g = 0} is set to 0.
  }
  \item{GRM}{
    For the GRM, the probability that an examinee with latent ability \eqn{\theta} responds to score category \eqn{k} (\eqn{k=0,1,...,K}) 
    of an item is a given by,
    \deqn{P(Y = k | \theta) = P^{*}(Y \ge k | \theta) - P^{*}(Y \ge k + 1 | \theta),}
    \deqn{P^{*}(Y \ge k | \theta) = \frac{1}{1 + exp(-Da(\theta - b_{k}))}, and}
    \deqn{P^{*}(Y \ge k + 1 | \theta) = \frac{1}{1 + exp(-Da(\theta - b_{k+1}))}, }
    
    where \eqn{P^{*}(Y \ge k | \theta} refers to the category boundary (threshold) function for score category \eqn{k} of an item
    and its formula is analogous to that of 2PLM. \eqn{b_{k}} is the difficulty (or threshold) parameter for category boundary 
    \eqn{k} of an item. Note that \eqn{P(Y = 0 | \theta) = 1 - P^{*}(Y \ge 1 | \theta)}
    and \eqn{P(Y = K-1 | \theta) = P^{*}(Y \ge K-1 | \theta)}. 
  }
  \item{GPCM}{
    For the GPCM, the probability that an examinee with latent ability \eqn{\theta} responds to score category \eqn{k} (\eqn{k=0,1,...,K}) 
    of an item is a given by,
     \deqn{P(Y = k | \theta) = \frac{exp(\sum_{v=0}^{k}{Da(\theta - b_{v})})}{\sum_{h=0}^{K-1}{exp(\sum_{v=0}^{h}{Da(\theta - b_{v})})}},}      
    where \eqn{b_{v}} is the difficulty parameter for category boundary \eqn{v} of an item. In other contexts, the difficulty parameter \eqn{b_{v}} 
    can also be parameterized as \eqn{b_{v} = \beta - \tau_{v}}, where \eqn{\beta} refers to the location (or overall difficulty) parameter 
    and \eqn{\tau_{jv}} represents a threshold parameter for score category \eqn{v} of an item. In the \pkg{irtplay} package, \eqn{K-1} difficulty 
    parameters are necessary when an item has \eqn{K} unique score categories because \eqn{b_{0}=0}. When a partial credit model is fit to data, \eqn{a} 
    is fixed to 1.
   }

}
}

\references{
Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional Module on Item-Fit Statistics for Item Response Theory Models.
\emph{Educational Measurement: Issues and Practice, 34}(3), 39-48.

Baker, F. B., & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques.} CRC Press.

Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.

Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
\emph{Statistical theories of mental test scores} (pp. 397-479). Reading, MA: Addison-Wesley.

Bock, R.D. (1960), \emph{Methods and applications of optimal scaling}. Chapel Hill, NC: L.L. Thurstone Psychometric Laboratory.

Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm.
\emph{Psychometrika, 46}, 443-459.

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Psychometrika, 35}, 179-198.

Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
Chapel Hill, NC: Vector Psychometric Group.

Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment.
\emph{Journal of Statistical Software, 48}(6), 1-29.

Chen, P., & Wang, C. (2016). A new online calibration method for multidimensional computerized adaptive testing.
\emph{Psychometrika, 81}(3), 674-701.

Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial.
\emph{Biometrika, 26}(4), 404-413.

Hambleton, R. K., & Swaminathan, H., & Rogers, H. J. (1991) \emph{Fundamentals of item response theory}.
Newbury Park, CA: Sage.

Han, K. T. (2016). Maximum likelihood score estimation method with fences for short-length tests and computerized adaptive tests.
\emph{Applied psychological measurement, 40}(4), 289-301.

Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item fit index for polytomous IRT models.
\emph{Journal of Educational Measurement, 45}(4), 391-406.

Kim, S. (2006). A comparative study of IRT fixed parameter calibration methods.
\emph{Journal of Educational Measurement, 43}(4), 355-381.

Kolen, M. J. & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
Springer.

Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
\emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.

Laplace, P. S. (1820).\emph{Theorie analytique des probabilites} (in French). Courcier.

Li, Y. & Lissitz, R. (2004). Applications of the analytically derived asymptotic standard errors of item response theory
item parameter estimates. \emph{Journal of educational measurement, 41}(2), 85-117.

Lim, H., Choe, E. M., & Han, K. T. (Under review). A residual-based differential item functioning detection framework in 
item response theory. \emph{Journal of Educational Measurement}. 

Lim, H., Choe, E. M., Han, K. T., Lee, S., & Hong, M. (2021, June). \emph{IRT residual approach 
to detecting DIF.} Paper presented at the Annual Meeting of the National Council on Measurement 
in Education. Online. 

Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical approach to evaluate the performance of MST.
\emph{Journal of Educational Measurement}. DOI: 10.1111/jedm.12276.

Lord, F. & Wingersky, M. (1984). Comparison of IRT true score and equipercentile observed score equatings.
\emph{Applied Psychological Measurement, 8}(4), 453-461.

McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
\emph{Applied Psychological Measurement, 9}, 49-57.

Meilijson, I. (1989). A fast improvement to the EM algorithm on its own terms.
\emph{Journal of the Royal Statistical Society: Series B (Methodological), 51}, 127-138.

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com

Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of seven methods.
\emph{Statistics in medicine, 17}(8), 857-872.

Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for dichotomous item response theory models.
\emph{Applied Psychological Measurement, 24}(1), 50-64.

Orlando, M., & Thissen, D. (2003). Further investigation of the performance of S-X2: An item fit index for use with
dichotomous item response theory models. \emph{Applied Psychological Measurement, 27}(4), 289-298.

Pritikin, J. (2018). \emph{rpf: Response Probability Functions}. R package version 0.59.
https://CRAN.R-project.org/package=rpf.

Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
\emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.

Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.

Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
\emph{Psychometrika, 47}, 175-186.

Thissen, D. & Wainer, H. (1982). Weighted likelihood estimation of ability in item response theory. 
\emph{Psychometrika, 54}(3), 427-450.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
Measurement, 19}(1), 39-49.

Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
\emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.

Wainer, H., & Mislevy, R. J. (1990). Item response theory, item calibration, and proficiency estimation. In H. Wainer (Ed.),
\emph{Computer adaptive testing: A primer} (Chap. 4, pp.65-102). Hillsdale, NJ: Lawrence Erlbaum.

Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests Using IRT-Based Methods.
\emph{Journal of Statistical Software, 35}(12), 1-33. URL http://www.jstatsoft.org/v35/i12/.

Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric procedure for assessing goodness-of-fit in
item response theory. \emph{Applied Measurement in Education, 21}(1), 22-40.

Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
\emph{Journal of the American Statistical Association, 22}(158), 209-212.

Woods, C. M. (2007). Empirical histograms in item response theory with ordinal data. \emph{Educational and Psychological Measurement, 67}(1), 73-87.

Yen, W. M. (1981). Using simulation results to choose a latent trait model. \emph{Applied Psychological Measurement, 5}, 245-262.

Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). BILOG-MG 3: Multiple-group
IRT analysis and test maintenance for binary items [Computer Program]. Chicago, IL: Scientific
Software International. URL http://www.ssicentral.com
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getirt.R
\name{getirt}
\alias{getirt}
\alias{getirt.est_irt}
\alias{getirt.est_item}
\title{Extract various elements from 'est_irt' or 'est_item' objects}
\usage{
getirt(x, ...)

\method{getirt}{est_irt}(x, what, ...)

\method{getirt}{est_item}(x, what, ...)
}
\arguments{
\item{x}{An object of class \code{\link{est_irt}} or \code{\link{est_item}}.}

\item{...}{Further arguments passed to or from other methods.}

\item{what}{A character string indicating what to extract.}
}
\description{
This function extracts various internal objects from an object of class \code{\link{est_irt}}
or \code{\link{est_item}}
.
}
\details{
Objects which can be extracted from the object of class \code{\link{est_irt}} include:

\describe{
\item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
\item{par.est}{A data frame containing the item parameter estimates.}
\item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
observed information functions. The standard errors are estimated using the cross-production approximation method (Meilijson, 1989).}
\item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
when interpreting the variance-covariance matrix of item parameter estimates.}
\item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
\item{loglikelihood}{A sum of the log-likelihood values of the observed data set (marginal log-likelihood) across all estimated items.}
\item{aic}{A model fit statistic of Akaike information criterion based on the loglikelihood.}
\item{bic}{A model fit statistic of Bayesian information criterion based on the loglikelihood.}
\item{group.par}{A data frame containing the mean, variance, and standard deviation of latent variable prior distribution.}
\item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
(in the second column) of the (updated) latent variable prior distribution.}
\item{posterior.dist}{A matrix of normalized posterior densities for all the response patterns at each of the quadrature points.
The row and column indicate the response pattern and the quadrature point, respectively.}
\item{data}{A data frame of the examinees' response data set.}
\item{scale.D}{A scaling factor in IRT models.}
\item{ncase}{A total number of response patterns.}
\item{nitem}{A total number of items included in the response data.}
\item{Etol}{A convergence criteria for E steps of the EM algorithm.}
\item{MaxE}{The maximum number of E steps in the EM algorithm.}
\item{aprior}{A list containing the information of the prior distribution for item slope parameters.}
\item{bprior}{A list containing the information of the prior distribution for item difficulty (or threshold) parameters.}
\item{gprior}{A list containing the information of the prior distribution for item guessing parameters.}
\item{npar.est}{A total number of the estimated parameters.}
\item{niter}{The number of EM cycles completed.}
\item{maxpar.diff}{A maximum item parameter change when the EM cycles were completed.}
\item{EMtime}{Time (in seconds) spent for the EM cycles.}
\item{SEtime}{Time (in seconds) spent for computing the standard errors of the item parameter estimates.}
\item{TotalTime}{Time (in seconds) spent for total compuatation.}
\item{test.1}{Status of the first-order test to report if the gradients has vanished sufficiently for the solution to be stable.}
\item{test.2}{Status of the second-order test to report if the information matrix is positive definite, which is a prerequisite
for the solution to be a possible maximum.}
\item{var.note}{A note to report if the variance-covariance matrix of item parameter estimates is obtainable from the information matrix.}
\item{fipc}{A logical value to indicate if FIPC was used.}
\item{fipc.method}{A method used for the FIPC.}
\item{fix.loc}{A vector of integer values specifying the location of the fixed items when the FIPC was implemented.}
}

Objects which can be extracted from the object of class \code{\link{est_item}} include:

\describe{
\item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
\item{par.est}{A data frame containing the item parameter estimates.}
\item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
observed information functions.}
\item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
when interpreting the variance-covariance matrix of item parameter estimates.}
\item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
\item{loglikelihood}{A sum of the log-likelihood values of the complete data set across all estimated items.}
\item{data}{A data frame of the examinees' response data set.}
\item{score}{A vector of the examinees' ability values used as the fixed effects.}
\item{scale.D}{A scaling factor in IRT models.}
\item{convergence}{A string indicating the convergence status of the item parameter estimation.}
\item{nitem}{A total number of items included in the response data.}
\item{deleted.item}{The items which have no item response data. Those items are excluded from the item parameter estimation.}
\item{npar.est}{A total number of the estimated parameters.}
\item{n.response}{An integer vector indicating the number of item responses for each item used to estimate the item parameters.}
\item{TotalTime}{Time (in seconds) spent for total compuatation.}
}

See \code{\link{est_irt}} and \code{\link{est_item}} for more details.
}
\section{Methods (by class)}{
\itemize{
\item \code{est_irt}: An object created by the function \code{\link{est_irt}}.

\item \code{est_item}: An object created by the function \code{\link{est_item}}.
}}

\examples{
\donttest{
# fit the 2PL model to LSAT6 data
mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2)

# extract the item parameter estimates
(est.par <- getirt(mod.2pl, what="par.est"))

# extract the standard error estimates
(est.se <- getirt(mod.2pl, what="se.est"))

# extract the variance-covariance matrix of item parameter estimates
(cov.mat <- getirt(mod.2pl, what="covariance"))
}

}
\seealso{
\code{\link{est_irt}}, \code{\link{est_item}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
