Installation
------------

This is an R package for estimation and prediction of the cosinor model for periodic data. It allows for covariate effects and provides tools for inference on the mean shift, amplitude, and acrophase. Check out the shiny app that illustrates the model here:

<http://sachsmc.shinyapps.io/cosinor-shinyapp/>

The package is on CRAN and can be installed as follows

``` r
install.packages("cosinor")
```

To install from github, you must first have the `devtools` package installed. Then run this command to install the latest version:

``` r
devtools::install_github("cosinor", "sachsmc")
```

Model details
-------------

For outcome $Y$, time variable $t$, and fixed period $D$ it is of interest to fit the periodic model

$$
Y(t) = \alpha + \beta_1 * \cos\{2 * \pi * t / D - \beta_2\} + \varepsilon
$$

where $\alpha$ is the intercept, $\beta_1$ is the amplitude, and $\beta_2$ is the acrophase (also called phase-shift). The $\varepsilon$ is an error term with mean 0.

This model transforms so that it can be fit using a simple linear model. Let $r = \cos\{2 * \pi * t / D\}$ and let $s = \sin\{2 * \pi * t / D\}$. Then we have

$$
Y(t) = \alpha + \gamma_1 * r + \gamma_2 * s + \varepsilon.
$$

The original coefficients can be recovered as follows:

$$
\beta_1 = \sqrt{\gamma_1^2 + \gamma_2^2}
$$

and

$$
\beta_2 = tan^{-1}(\gamma_2 / \gamma_1).
$$

In the package, $(\alpha, \gamma_1, \gamma_2)$ is estimated using `lm`, and inference on the transformed coefficients is obtained using the delta method.

Example usage
-------------

Load the package into your library:

``` r
library(cosinor)
```

The package comes with an example dataset called `vitamind` to illustrate the usage. Fit a basic cosinor model with covariates:

``` r
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind, period = 12)
```

The time variable is indicated by the `time` function in the formula. By default, the period length is 12. Covariates can be included directly in the formula to allow for mean shifts. To allow for differences in amplitude and acrophase by covariate values, include the covariate wrapped in the `amp.acro` function in the formula. Let's summarize the model.

``` r
summary(fit)
```

    ## Raw model coefficients:
    ##             estimate standard.error lower.CI upper.CI p.value
    ## (Intercept)  29.6898         0.4654  28.7776  30.6020  0.0000
    ## X             1.9019         0.8041   0.3258   3.4779  0.0180
    ## rrr           0.9308         0.6357  -0.3151   2.1767  0.1431
    ## sss           6.2010         0.6805   4.8673   7.5347  0.0000
    ## X:rrr         5.5795         1.1386   3.3479   7.8111  0.0000
    ## X:sss        -1.3825         1.1364  -3.6097   0.8447  0.2237
    ## 
    ## ***********************
    ## 
    ## Transformed coefficients:
    ##             estimate standard.error lower.CI upper.CI p.value
    ## (Intercept)  29.6898         0.4654  28.7776  30.6020   0.000
    ## [X = 1]       1.9019         0.8041   0.3258   3.4779   0.018
    ## amp           6.2705         0.6799   4.9378   7.6031   0.000
    ## [X = 1]:amp   8.0995         0.9095   6.3170   9.8820   0.000
    ## acr           1.4218         0.1015   1.2229   1.6207   0.000
    ## [X = 1]:acr   0.6372         0.1167   0.4084   0.8659   0.000

This prints the raw coefficients, and the transformed coefficients. The transformed coefficients display the amplitude and acrophase for the covariate = 1 group and the covariate = 0 group. Beware when using continuous covariates!

Testing is easy to do, tell the function which covariate to test, and whether to test the amplitude or acrophase. The global test is useful when you have multiple covariates in the model.

``` r
test_cosinor(fit, "X", param = "amp")
```

    ## Global test: 
    ## Statistic: 
    ## [1] 2.59
    ## 
    ## 
    ##  P-value: 
    ## [1] 0.1072
    ## 
    ##  Individual tests: 
    ## Statistic: 
    ## [1] 1.61
    ## 
    ## 
    ##  P-value: 
    ## [1] 0.1072
    ## 
    ##  Estimate and confidence interval[1] "1.83 (-0.4 to 4.05)"

The predict function allows you to estimate the mean value of your outcome for individuals. This is useful for adjusting for the seasonal effects on some measurement. Let's compare the raw values to the seasonally adjusted values of the vitamin D dataset.

``` r
summary(vitamind$Y)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   13.27   25.59   29.79   30.09   34.85   51.30

``` r
summary(predict(fit))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   16.10   25.98   29.26   29.69   33.59   46.48

Plotting the fitted curves is also easy to do. Currently only `ggplot2` is supported.

``` r
library(ggplot2)
ggplot.cosinor.lm(fit, x_str = "X")
```

![](README_files/figure-markdown_github+tex_math_dollars/unnamed-chunk-8-1.png)
---
title: "Cosinor"
author: "Michael Sachs"
date: "July 23, 2014"
output: 
    md_document:
        variant: markdown_github+tex_math_dollars
---


## Installation

This is an R package for estimation and prediction of the cosinor model for periodic data. It allows for covariate effects and provides tools for inference on the mean shift, amplitude, and acrophase. Check out the shiny app that illustrates the model here:

http://sachsmc.shinyapps.io/cosinor-shinyapp/

The package is on CRAN and can be installed as follows

```{r eval = FALSE}
install.packages("cosinor")
```

To install from github, you must first have the `devtools` package installed. Then
run this command to install the latest version:

```{r eval = FALSE}
devtools::install_github("cosinor", "sachsmc")
```

## Model details

For outcome $Y$, time variable $t$, and fixed period $D$ it is of interest to fit the periodic model

$$
Y(t) = \alpha + \beta_1 * \cos\{2 * \pi * t / D - \beta_2\} + \varepsilon
$$

where $\alpha$ is the intercept, $\beta_1$ is the amplitude, and $\beta_2$ is the acrophase (also called phase-shift). The $\varepsilon$ is an error term with mean 0. 

This model transforms so that it can be fit using a simple linear model. Let $r = \cos\{2 * \pi * t / D\}$ and let $s = \sin\{2 * \pi * t / D\}$. Then we have 

$$
Y(t) = \alpha + \gamma_1 * r + \gamma_2 * s + \varepsilon.
$$

The original coefficients can be recovered as follows:

$$
\beta_1 = \sqrt{\gamma_1^2 + \gamma_2^2}
$$

and

$$
\beta_2 = tan^{-1}(\gamma_2 / \gamma_1).
$$

In the package, $(\alpha, \gamma_1, \gamma_2)$ is estimated using  `lm`, and inference on the transformed coefficients is obtained using the delta method. 


## Example usage

Load the package into your library:

```{r}
library(cosinor)
```

The package comes with an example dataset called `vitamind` to illustrate the usage. Fit a basic cosinor model with covariates:

```{r}
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind, period = 12)
```

The time variable is indicated by the `time` function in the formula. By default, the period length is 12. Covariates can be included directly in the formula to allow for mean shifts. To allow for differences in amplitude and acrophase by covariate values, include the covariate wrapped in the `amp.acro` function in the formula. Let's summarize the model.

```{r}
summary(fit)
```

This prints the raw coefficients, and the transformed coefficients. The transformed coefficients display the amplitude and acrophase for the covariate = 1 group and the covariate = 0 group. Beware when using continuous covariates!

Testing is easy to do, tell the function which covariate to test, and whether to test the amplitude or acrophase. The global test is useful when you have multiple covariates in the model.

```{r}
test_cosinor(fit, "X", param = "amp")
```

The predict function allows you to estimate the mean value of your outcome for individuals. This is useful for adjusting for the seasonal effects on some measurement. Let's compare the raw values to the seasonally adjusted values of the vitamin D dataset.

```{r}
summary(vitamind$Y)
summary(predict(fit))
```

Plotting the fitted curves is also easy to do. Currently only `ggplot2` is supported.

```{r}
library(ggplot2)
ggplot.cosinor.lm(fit, x_str = "X")
```


<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{Usage}
-->


This is an R package for estimation and prediction of the cosinor model for periodic data. It allows for covariate effects and provides tools for inference on the mean shift, amplitude, and acrophase. 

The package can be installed from CRAN:

```{r eval = FALSE}
install.packages("cosinor")
```


To install from github, you must first have the `devtools` package installed. Then run this command to install the latest version:

```{r eval = FALSE}
devtools::install_github("cosinor", "sachsmc")
```

Load the package into your library:

```{r}
library(cosinor)
```

The package comes with an example dataset called `vitamind` to illustrate the usage. Fit a basic cosinor model with covariates:

```{r}
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind, period = 12)
```

The time variable is indicated by the `time` function in the formula. By default, the period length is 12. Covariates can be included directly in the formula to allow for mean shifts. To allow for differences in amplitude and acrophase by covariate values, include the covariate wrapped in the `amp.acro` function in the formula. Let's summarize the model.

```{r}
summary(fit)
```

This prints the raw coefficients, and the transformed coefficients. The transformed coefficients display the amplitude and acrophase for the covariate = 1 group and the covariate = 0 group. Beware when using continuous covariates!

Testing is easy to do, tell the function which covariate to test, and whether to test the amplitude or acrophase. The global test is useful when you have multiple covariates in the model.

```{r}
test_cosinor(fit, "X", param = "amp")
```

The predict function allows you to estimate the mean value of your outcome for individuals. This is useful for adjusting for the seasonal effects on some measurement. Let's compare the raw values to the seasonally adjusted values of the vitamin D dataset.

```{r}
summary(vitamind$Y)
summary(predict(fit))
```

Plotting the fitted curves is also easy to do. Currently only `ggplot2` is supported.

```{r}
library(ggplot2)
ggplot.cosinor.lm(fit, x_str = "X")
```

The package comes with an interactive web app that can be used to analyze your own data. If you have your data loaded in R as a `data.frame`, run the command 

```{r eval = FALSE}
cosinor_analyzer(vitamind)
```

to run the Shiny application. 

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{predict.cosinor.lm}
\alias{predict.cosinor.lm}
\title{Predict from a cosinor model}
\usage{
\method{predict}{cosinor.lm}(object, newdata, ...)
}
\arguments{
\item{object}{An object of class \code{cosinor.lm}}

\item{newdata}{Optional new data}

\item{...}{other arguments}
}
\description{
Given a time variable and optional covariates, generate predicted values from
a cosinor fit. Default prediction is the mean value, optionally can predict
at a given month
}
\examples{
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
predict(fit)
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{cosinor.lm}
\alias{cosinor.lm}
\title{Fit cosinor model}
\usage{
cosinor.lm(formula, period = 12, data, na.action = na.omit)
}
\arguments{
\item{formula}{Forumla specifying the model. Indicate the time variable with
\code{time()} and covariate effects on the amplitude and acrophase with
\code{amp.acro()}. See details for more information.}

\item{period}{Length of time for a complete period of the sine curve.}

\item{data}{Data frame where variable can be found}

\item{na.action}{What to do with missing data}
}
\description{
Given an outcome and time variable, fit the cosinor model with optional
covariate effects.
}
\details{
This defines special functions that are used in the formula to
  indicate the time variable and which covariates effect the amplitude. To
  indicate the time variable wrap the name of it in the function
  \code{time()}. To indicate a variable which affects the
  acrophase/amplitude, wrap the name in \code{amp.acro()}. This will then do
  all the tranformations for you. See examples for usage.
}
\examples{
cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
}
\references{
Tong, YL. Parameter Estimation in Studying Circadian Rhythms, Biometrics (1976). 32(1):85--94.
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{test_cosinor}
\alias{test_cosinor}
\title{Test for differences in a cosinor model}
\usage{
test_cosinor(object, x_str, param = "amp")
}
\arguments{
\item{object}{An object of class \code{cosinor.lm}}

\item{x_str}{Character naming the covariate whose amplitude/acrophase will be tested}

\item{param}{Character string naming the parameter to test, either "amp" for
  amplitude or "acr" for acrophase}
}
\description{
Given a time variable and optional covariates, generate inference a cosinor
fit. For the covariate named (or vector of covariates), this function
performs a Wald test comparing the group with covariates equal to 1 to the
group with covariates equal to 0. This may not be the desired result for
continuous covariates.
}
\examples{
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
test_cosinor(fit, "X", "amp")
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{get_varnames}
\alias{get_varnames}
\title{Extract variable names from terms object, handling specials}
\usage{
get_varnames(Terms)
}
\arguments{
\item{Terms}{a terms object}
}
\description{
Extract variable names from terms object, handling specials
}
\keyword{Internal}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{update_covnames}
\alias{update_covnames}
\title{Replace covariate names with descriptive text}
\usage{
update_covnames(names)
}
\arguments{
\item{names}{Coefficient names to update}
}
\description{
Replace covariate names with descriptive text
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{print.cosinor.lm}
\alias{print.cosinor.lm}
\title{Print cosinor model}
\usage{
\method{print}{cosinor.lm}(x, ...)
}
\arguments{
\item{x}{cosinor.lm object}

\item{...}{passed to summary}
}
\description{
Given an outcome and time variable, fit the cosinor model with optional covariate effects.
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{summary.cosinor.lm}
\alias{summary.cosinor.lm}
\title{Summarize a cosinor model}
\usage{
\method{summary}{cosinor.lm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{cosinor.lm}}

\item{...}{Currently unusued}
}
\description{
Given a time variable and optional covariates, generate inference a cosinor
fit. Gives estimates, confidence intervals, and tests for the raw parameters,
and for the mean, amplitude, and acrophase parameters. If the model includes
covariates, the function returns the estimates of the mean, amplitude, and
acrophase for the group with covariates equal to 1 and equal to 0. This may
not be the desired result for continuous covariates.
}
\examples{
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
summary(fit)
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{simulate_cosinor}
\alias{simulate_cosinor}
\title{Simulate data from a cosinor model}
\usage{
simulate_cosinor(n, beta.mean = 2, beta.amp = 0, beta.acro = 0)
}
\arguments{
\item{n}{Sample size}

\item{beta.mean}{Effect on the mean (intercept)}

\item{beta.amp}{Effect on the amplitude}

\item{beta.acro}{Effect on the acrophase}
}
\description{
This function simulates data from a cosinor model with a single covariate,
where the time scale is month, and optionally
allows for single covariate effects on the mean,
amplitude, and acrophase.
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\docType{data}
\name{vitamind}
\alias{vitamind}
\title{Vitamin D}
\format{A data frame with 3 variables: \code{X}, \code{Y},
  \code{time}.}
\usage{
vitamind
}
\description{
Simulated data set to illustrate the cosinor model. \code{Y}
is an outcome variable that varies of time \code{time} according
to a cosine curve. The binary covariate \code{X} is associated with the
mean and amplitude of the cosine curve.
}
\keyword{datasets}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{ggplot.cosinor.lm}
\alias{ggplot.cosinor.lm}
\title{Plot a cosinor model}
\usage{
ggplot.cosinor.lm(object, x_str = NULL)
}
\arguments{
\item{object}{An object of class \code{cosinor.lm}}

\item{x_str}{Character vector naming the covariate(s) to be plotted. May be NULL to plot overall curve}
}
\description{
Given a cosinor.lm model fit, generate a plot of the data with the fitted values.
Optionally allows for plotting by covariate levels 0 and 1.
}
\examples{
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
ggplot.cosinor.lm(fit, "X")
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{cosinor.lm.default}
\alias{cosinor.lm.default}
\title{Fit cosinor model}
\usage{
cosinor.lm.default(formula, ...)
}
\arguments{
\item{formula}{Forumla specifying the model. Indicate the time variable with \code{time()} and covariate effects on the
amplitude and acrophase with \code{amp.acro()}. See details.}

\item{...}{other arguments}
}
\description{
Given an outcome and time variable, fit the cosinor model with optional covariate effects.
}
\details{
This defines special functions that are used in the formula to indicate the time variable
and which covariates effect the amplitude. To indicate the time variable wrap the name of it in the function
\code{time()}. To indicate a variable which affects the acrophase/amplitude, wrap the name in
\code{amp.acro()}. This will then do all the tranformations for you. See examples for usage.
}
\examples{
cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{print.test_cosinor}
\alias{print.test_cosinor}
\title{Print results of test of cosinor model}
\usage{
\method{print}{test_cosinor}(x, ...)
}
\arguments{
\item{x}{test_cosinor object}

\item{...}{Arguments passed to \code{print}}
}
\description{
Print results of test of cosinor model
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{print.test}
\alias{print.test}
\title{Print test of model}
\usage{
\method{print}{test}(x)
}
\arguments{
\item{x}{test object}
}
\description{
Print test of model
}
\keyword{Internal}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{print.summary.cosinor.lm}
\alias{print.summary.cosinor.lm}
\title{Print the summary of a cosinor model}
\usage{
\method{print}{summary.cosinor.lm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{summary.cosinor.lm}}

\item{...}{Currently unusued}
}
\description{
Print the summary of a cosinor model
}
\examples{
fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
summary(fit)
}

% Generated by roxygen2 (4.0.1): do not edit by hand
\name{cosinor_analyzer}
\alias{cosinor_analyzer}
\title{Shiny application to demonstrate cosinor fit}
\usage{
cosinor_analyzer(data = vitamind)
}
\arguments{
\item{data}{Data frame to analyze}
}
\description{
Given a dataset, specify the outcome, time variable, and optional covariates. The app will then perform a cosinor analysis and plot the results.
}
\examples{
\dontrun{
library(shiny)
cosinor_analyzer(vitamind)
}
}

