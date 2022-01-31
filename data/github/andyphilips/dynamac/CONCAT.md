<img src="https://andyphilips.github.io/dynamac/img/logo.png" alt="drawing" width="200"/>

## Description
`dynamac` performs dynamic simulation and testing for single-equation ARDL models in R and Stata. Link to the webpage can be found [here](https://andyphilips.github.io/dynamac/).

## Installation (R)
[![Build Status](https://travis-ci.com/andyphilips/dynamac.svg?branch=master)](https://travis-ci.com/andyphilips/dynamac) [![Build status](https://ci.appveyor.com/api/projects/status/o8h5gdh5cuah359y?svg=true)](https://ci.appveyor.com/project/andyphilips/dynamac) ![](https://www.r-pkg.org/badges/version/dynamac) [![codecov](https://codecov.io/gh/andyphilips/dynamac/branch/master/graph/badge.svg)](https://codecov.io/gh/andyphilips/dynamac) [![metacran downloads](https://cranlogs.r-pkg.org/badges/dynamac)](https://cran.r-project.org/package=dynamac) [![doi](https://img.shields.io/badge/DOI:-10.32614%2FRJ--2018--076-blueviolet.svg)](https://doi.org/10.32614/RJ-2018-076) [![DOI](https://zenodo.org/badge/129036070.svg)](https://zenodo.org/badge/latestdoi/129036070)

`dynamac` in R can easily be installed from its CRAN repository:
```
install.packages("dynamac")
library(dynamac)
```

Alternatively, you can use the `devtools` package to load directly from GitHub:
```
library(devtools)
devtools::install_github("andyphilips/dynamac")
library(dynamac)
```
You should now be able to use the package.

## Installation (Stata)
To install `dynamac` in Stata type one of the following into the console (note that using `net install` calls directly from Github and will ensure the most up-to-date version):
```
net sj 18-4 st0545
```
```
findit dynamac
```
```
cap ado uninstall dynamac
cap ado uninstall dynardl
cap ado uninstall pssbounds
net install dynamac, from(https://github.com/andyphilips/dynamac/raw/gh-pages/Stata/src/)
```

Alternatively, you download the [raw files](https://github.com/andyphilips/dynamac/raw/gh-pages/Stata/src/) and either call directly to the .ado files or place them in your "ado/plus/" folder.

## I've found a bug!
Great! Please post it on the [issues page](https://github.com/andyphilips/dynamac/issues), or email the authors.
--- 
title: 'Exploring meaningful visual effects and quantities of interest from dynamic models through dynamac'  
authors:
  - name: Soren Jordan 
    affiliation: 1 
    orcid: 0000-0003-4201-1085 
  - name: Andrew Q. Philips 
    affiliation: 2 
affiliations: 
  - index: 1 
    name: Assistant Professor, Department of Political Science, Auburn University 
  - index: 2 
    name: Assistant Professor, Department of Political Science, University of Colorado Boulder 
date: 8 October 2020 
tags: 
  - time series
  - methodology
  - long-run effects
  - simulations
  - quantities of interest
bibliography: jordan-philips.bib 
---

# Summary 

dynamac [-@dynamacCRAN] implements a framework of estimating and interpreting autoregressive distributed lag (ARDL) models in R. dynamac uses stochastic simulation techniques [@jordan2018cointegration;@jordan2018dynamic] to easily recover traditional quantities of interest, such as short- and long-run effects, even from complex dynamic specifications, including models with cointegration. These simulation techniques bring a wider set of inferences from complex models to users in fields as diverse as environmental science [@jcp] and economics [@te]. Users can make dynamic inferences from their models, including measures of uncertainty, without the need for any formulae or algebraic solutions. Although other packages may estimate ARDL models (e.g., dynlm [@dynlmCRAN], which only provides model estimates; ardl [@ardlCRAN], which uses dynlm for estimation but provides support for the bounds test for cointegration and automated lag selection; dlagM [@dlagMCRAN], which allows for forecasting), dynamac is unique in that it provides post-estimation diagnostics for autocorrelation in the residuals and is designed specifically for providing inferences through counterfactual simulation.

The procedure first estimates the parameters using an ARDL model:
\begin{align}\label{eq:ardlgeneral}
\begin{split}
	y_t,\Delta y_t =  &\alpha_0 + \delta T + \sum_{p=1}^P \phi_p y_{t-p} + \sum_{l_1 = 0}^{L_1} \theta_{1l}x_{1,t-l_1} + \cdots + \sum_{l_k = 0}^{L_k} \theta_{kl}x_{t-l_k} + \\
	& \sum_{m=1}^M \alpha_m \Delta y_{t-m} + \sum_{q_1=0}^{Q_1} \beta_{1q_{1}} \Delta x_{1,t-q{_1}} + \cdots + \sum_{q_k=0}^{Q_k} \beta_{kq_k} \Delta x_{k,t-q{_k}} + \epsilon_t
\end{split}
\end{align} The dependent variable (appearing in level form or first-differences) is a function of a constant, a deterministic linear trend $T$, up to $P$ and $L$ lags of the dependent and independent variables, up to $M$ and $Q$ lags of the first-differenced dependent and independent variables, and  error term $\epsilon_t$.  Users should use a combination of theory, information criteria, unit-root, cointegration and residual-based tests to arrive at a more restrictive specification of Equation \ref{eq:ardlgeneral}.

Next, through `dynardl()`, the estimated parameters are used to draw $s$ simulations from a multivariate normal distribution, with mean $\hat{\beta}$ and variance obtained from the estimated variance-covariance matrix. Stable predicted values for $y$ using the $s$ simulations of $\hat{\beta}$ are then created, using starting values of the independent variables (usually means for continuous variables or modes for categorical ones). Users can choose from either expected values or---by incorporating fundamental model uncertainty---predicted values. Next, at a defined time period $t$, one of the independent variables is "shocked" by an amount determined by the researcher, and the effect on the dependent variable is observed across the resulting time periods graphically through `dynardl.simulation.plot()`.

This functionality simplifies the production and interpretation of dynamic models to a broader class of users. Basic interpretation of dynamic models usually involves calculating the short-run effect of a one-unit increase in an independent variable, as well as the long-run effect. The former is easy to estimate and observe (represented by a $\hat\beta$ coefficient in an ARDL model), but the latter---and its associated confidence intervals---are harder to obtain, since they involve a non-linear combination of parameter estimates. Inferences can also be obscured due to the complexity of the model (like including multiple lags and/or first-differences of dependent and independent variables). For short-run effects, increased model complexity encourages a focus on "boring" effects (the one-unit effect reported by the $\hat\beta$ in the regression output); for long-run effects, greater model complexity results in increasingly intractable closed-form estimates of the effects. Stochastic simulations can solve both of these problems by allowing for short-run dynamics beyond the traditional one-unit effect, and also by removing the need for analytic calculations when discussing long-run effects. 

# Plotting functionality
dynamac implements six new visualizations. These can be presented one at a time using the command `dynardl.simulation.plot()`, or all six using `dynardl.all.plots()`. The resulting plots are:

* The response path of the level of the dependent variable, given a change (a "shock", such as a one standard deviation increase or decrease) in a single independent variable at one point in time.
* Deviations between the predicted and average values of $y_t$: $\hat{y}_t - \bar{y}$. While this plot follows a response path identical to the plot above, if a shock dissipates and the series reverts to its mean, this is easier to observe if we do not have to subtract the stable starting value from the plot. If the series reverts to a value other than its mean, we would be especially interested in this new long-run mean in response to the shock. This is analogous to an impulse response function.
* The period-to-period response in $y_t$ by differencing the response path. This allows us to visualize successive movements over time.
* The size of the shock itself over time. This allows us to visualize dynamic persistence, so we can make statements like "how many time points until the shock is at 50 percent of its original magnitude?" or "how many time points until the shock is effectively zero?" This information is poorly represented by either the short-run or long-run effects as traditionally calculated.
* The cumulative response in $y_t$ to a shock in $x_t$. This helps visualize what the "final" effect of a shock is on the dependent variable. A variable might first cause a positive response, but ultimately may be overwhelmed by a negative movement. Unlike other plots, this considers the cumulative histories in each individual simulation when plotting the response in $y$ over time, making the confidence intervals more conservative. This is analogous to the traditional long-run effect of an ARDL model.
* Total movement that occurs in $y_t$ given a change in $x_t$, which is analogous to the "absolute" effect of a shock to $x_t$ on $y_t$. In this setup, each of the period-over-period changes sums towards the "total" movement in $y_t$, so that all of the movement in $y$ is attributed as an effect to the independent variable.

![Six quantities of interest from the ARDL equation $\Delta y_t = -0.8 y_{t-1} -2 \Delta x_t  + x_{t-1} + u_t$.](allplots-revised.pdf)

Figure 1 plots these visualizations from a hypothetical ARDL model. Starting with the top-left plot and moving clockwise, users can easily understand how a dependent variable responds over time to a shock in an independent variable, how it changes away from its average pre-shock value, how the change itself changes over time, the cumulative absolute changes in the dependent variable, the cumulative nature of the changes in the dependent variable, and how the shock decays over time. All six graphical quantities of interest provide a richer set of information about our inferences. These require no additional mathematical burden on the user, regardless of the underlying model complexity. 

# Conclusion

Dynamic models are complicated since they have effects that do not lend themselves to straightforward interpretation. Stochastic simulations can move our interpretations forward, even in the face of increased model complexity. As we have shown, there are six quantities of interest that are of interest to users, since they provide us with a richer understanding of the dynamic process taking place. dynamac now includes the ability to obtain all of these visualizations described above in a single function: ``dynardl.all.plots()``. This will help users expand the different types of quantities of interest they can use to assess statistical and substantive significance.

# References

---
title: "An Introduction to dynamac: Dynamic Inferences (and Cointegration Testing) from Autoregressive Distributed Lag Models"
author: "Soren Jordan and Andrew Q. Philips"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: timeseriesrj.bib
vignette: >
  %\VignetteIndexEntry{dynamac Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Autoregressive distributed lag (ARDL) models are an integral part of estimating scientific processes over time. However, as we extend their usefulness by adding richness in dynamic specifications (through multiple lags of variables, either in levels or differences, or lags of the dependent variable), we begin to challenge our ability to draw meaningful inferences from coefficients alone. Variables may appear in multiple time periods, in multiple forms, and might even be filtered through lagged values of the dependent variable. Coefficients tell us about the immediate effect of some variable but have little to say about the long-run effect. 

There is a better solution. Instead of performing intense operations to develop a closed-form or algebraic solution, we can rely on the power of computing to *simulate* the over-time response in the dependent variable of some model, given a corresponding change in an $x$ variable. We often call these "changes" in $x$ variables, and the associated response in the dependent variable $y$, a "counterfactual response" (a simulated response to a "shock" we control).

dynamac helps simulate these counterfactuals. More generally, though, it is built to make using and drawing inferences from single-equation ARDL models as easy as possible. Moreover, it helps users implement the useful cointegration test from Pearson, Shin, and Smith [-@pesaran2001bounds]: the ARDL-bounds testing procedure.

We illustrate the usefulness of these functions through example. After a brief discussion of the ARDL model in the general sense, we estimate a collection of these models and demonstrate both the challenges of the ARDL model in the abstract and the solutions of dynamac in particular. 

## ARDL models generally

The ARDL model has a general form where $y$, modeled in levels or differences, is a function of itself (in lagged levels or differences), up to $k$ variables $x$, either in contemporaneous (same period, or appearing at time $t$) levels, lagged levels, contemporaneous differences, or lagged differences. Conventionally, the number of lags of the dependent variable in levels is counted by $p$, the number of lags of the dependent variable in differences is counted by $m$, the number of lags of the independent variables in levels is counted by $l$, and the number of lags of the independent variables in differences is counted by $q$ (note: $l$ and $q$ can begin at 0, i.e, contemporaneous effects appearing at time $t$). 

The number of time periods of each variable can, theoretically, be quite large. Or, put differently, $p$, $m$, $l$, and $q$, especially $l$ and $q$, might be hard to account for. Analysts without strong theory about the nature of the effects often begin with restricting all but the contemporaneous and first lag of each series, or an ARDL model of the nature 

$$y_t =  \alpha_0 + \phi_1 y_{t-1} + \theta_{1,0} x_{1,t} + \theta_{1,1} x_{1,t-1} + \cdots + \theta_{k,0}x_{k,t} + \theta_{k,1}x_{k,t-1}+ \beta *T  + \epsilon_t$$
where $\alpha_0$ is a constant and $\beta *T$ is a trend term. (The exact nature of this equation depends on whether $y$ is stationary or integrated, as well as if differences or lagged differences are entered into the model. But we will get to this below.)

It's useful at this point to stop and think: if I have multiple lags of a variable $x$, and they are filtered through a lagged dependent variable $y_{t-1}$, it might be difficult to get a sense of the "total" effect of $x$ on $y$. This becomes more and more difficult as $x$ increases in lagged levels $l$ *and potentially also* lagged differences $q$. Our coefficient estimates, while estimated, are no longer as useful in direct interpretation. Put differently, the very flexibility of the ARDL model also undoes its usefulness in interpretation! So we might seek an alternative way of interpreting these models. dynamac provides a unified way of estimating ARDL models and interpreting their effects. It also provides a way of implementing a popular test for cointegration. We uncover these methods below. 

## Estimating an ARDL model: best practices and dynamac

dynamac includes two datasets. We will use the Wright [-@wrightpolitical] dataset on income inequality. We can look at the first few rows of this dataset

```{r}
library(dynamac)
head(ineq)
```

`concern` is public concern about income inequality; `incshare10` is the income share of the top ten percent; `urate` is the unemployment rate. Wright [-@wrightpolitical] argues that concern about income inequality grows as inequality itself worsens and economic conditions deteriorate. A simple model, then, is that concern is a function of past values of itself, the level of income share held by the top ten percent, *changing* levels of income share held by the top ten percent, as well *changing* levels of unemployment in the short term. 

$$\Delta Concern_t = \alpha_0 + \phi_1 Concern_{t-1} + \theta_1 Income Top 10_{t-1} + \beta_1 \Delta Income Top 10_t + \beta_2 \Delta Unemployment_t + \epsilon_t$$

where the residuals are white noise. Let's develop the corresponding ARDL model using dynamac.

### Understanding our time series

Step 1 in any time series analysis is a visual inspection of the series coupled with formal tests of stationarity: whether the series has a constant mean, variance, and covariance over time (so that it reverts back to mean level), or if the series instead violates any of these three conditions. We advocate for using the urca package for this [@pfaff2016package], which includes a variety of tools and tests. Simply plotting the series reveals the following:

```{r}
library(urca)
ts.plot(ineq$concern)
ts.plot(ineq$incshare10)
ts.plot(ineq$urate)
```

None of the three series looks especially well-behaved. `concern` rose quickly and has been moving sluggishly since. `incshare10` has only grown over time, which cannot be mean-reverting (like a stationary series). `urate` experiences steep peaks and valleys with significant interludes in between. All three series look integrated. But we should be more discerning by using empirical tests, also included in urca. So we can execute the Augmented Dickey-Fuller (ADF) test, Phillips-Perron test, DF-GLS test, and KPSS test on each series. On `concern` this looks like

```{r}
summary(ur.df(ineq$concern, type = c("none"), lags = 1))
summary(ur.pp(ineq$concern, type = c("Z-tau"), model = c("constant"), use.lag = 1))
summary(ur.ers(ineq$concern, type = c("DF-GLS"), model = c("constant"), lag.max = 1)) 
summary(ur.kpss(ineq$concern, type = c("mu"), use.lag = 1))
```

The ADF test, Phillips-Perron test, and DF-GLS test all have null hypotheses of a unit root, all of which we fail to reject. The KPSS test has a null hypothesis of stationarity which we *do* reject. We have compelling evidence, then, of integration (I[1]) in the series `concern`. We can check whether differencing is enough to make the series stationary by executing the same tests

```{r eval = FALSE}
summary(ur.df(diff(ineq$concern), type = c("none"), lags = 1))
summary(ur.pp(diff(ineq$concern), type = c("Z-tau"), model = c("constant"), use.lag = 1))
summary(ur.ers(diff(ineq$concern), type = c("DF-GLS"), model = c("constant"), lag.max = 1)) 
summary(ur.kpss(diff(ineq$concern), type = c("mu"), use.lag = 1))
```

Each of the above tests, when run, indicated that the differenced series of `concern` is stationary. Having gathered the basic information about the nature of the history of our variables, we might be itching to estimate some preliminary models. To this point, R hasn't offered much in the way of improving beyond the basic `lm()` framework for using regression to estimate time series models in a consistent syntax. We elaborate on this problem and illustrate our solution, `dynardl`.

### Estimating ARDL models with `dynardl` 

ARDL models are flexible, but their flexibility often results in variables of different lengths due to differencing and lagging. For instance, consider our simple model from above where `incshare10` appears as a first lag and a first difference, `urate` appears as a first difference, there is a lagged dependent variable, and `concern` is the dependent variable in differences. We can introduce a lag through the built-in `lshift` function in dynamac. The syntax is just `lshift(variable, num.lags)`, where `num.lags` is the number of periods for the variable to be lagged. We can also difference a series through `dshift`. The syntax is just `dshift(variable)` for a first difference. For instance, 

```{r}
head(ineq$incshare10)
head(lshift(ineq$incshare10, 1))
head(dshift(ineq$incshare10))
```

So the syntax for the simple model described would be easy to write. But notice the problem

```{r, eval = TRUE, error = TRUE}
summary(lm(diff(ineq$concern) ~ lshift(ineq$concern, 1) + lshift(ineq$incshare10, 1) + dshift(ineq$incshare10) + dshift(ineq$urate)))
```

Introducing the lag changed the variable lengths, and we probably don't want to have to think about time series operation in our model specification, anyway. Other software has unified this operation by introducing a common syntax, like d. for differencing and l. for lagging, or even l.d. for lag differences (what program could that be?). In R, though, it remains more of a nuisance than it should be. The `dynardl` function helps to remedy this challenge by encouraging the user only to think about model structure in exactly the way outlined above: which variables $x$ get lags $l$ and lagged differences $q$, etc. The function expects a basic formula of all the variables in the model, the data set (`data`), if there is one, then lagged levels (`lags`) and lagged differences (`lagdiffs`) as a list and differences (`diffs`) and contemporaneous levels (`levels`) as a vector. The sole exception is if the user wants to run the model in differences, s/he needs to specify `ec = TRUE` (in an effort to force us to think critically about the error-correction form of the model). For instance, our above model now becomes (ignoring the `simulate` option for now)

```{r, eval = FALSE}
dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        ec = TRUE, simulate = FALSE)
```

Purely hypothetically, if we wanted more lagged levels of `concern`, we would just change `1` to `c(1:5)` for lags at $t-1$ to $t-5$, or any other number of lags. If we wanted to include the first-difference of `urate` lagged at periods $t-1$ and $t-2$, we would now run

```{r, eval = FALSE}
dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("urate" = c(1, 2)),
        ec = TRUE, simulate = FALSE)
```

If the variable can appear in multiple periods (`lags` and `lagdiffs`), it must be specified as a list. If it can only appear contemporaneously (`levels` and `diffs`), it must be specified as a vector. (The alternative was to specify `levels` through a `lag` at time 0: this did not seem practical.) We can also add or suppress a constant from the model by specifying `constant = TRUE/FALSE`, and we can do the same with a linear trend term with `trend = TRUE/FALSE` (by default, constants are included and trends are not). 

The option `ec` specifies the nature of the dependent variable. If it is possibly error-correcting (`ec = TRUE`), it is run in differences, or period-over-period changes. If it is not (`ec = FALSE`), the dependent variable is run in levels.

At this point, `dynardl` is just a unifying way of thinking about time series models for estimation. Yet it is going to offer two other important advantages: interpreting the effects of our variables, and executing a powerful test for cointegration. We will start with the latter. Remember testing for I(1) processes earlier: each test indicated that `concern` was I(1). Running the same tests for each of `incshare10` and `urate` suggests that all three series are I(1). This might cause us concern about cointegration: the special relationship between I(1) series where the series are in a long-run relationship, even if they move apart in the short term.  Cointegrating series are difficult to *uncover*. Traditional methods, like the Engle and Granger [-@engle1987co] two-step method or likelihood-based approaches [@johansen1995likelihood] too often conclude cointegration when it does not exist, at least in smaller sample sizes [@philips2016have]. An alternative test [@pesaran2001bounds], which we refer to as the ARDL-bounds procedure, performs much better in small samples ($t < 80$), but it is more cumbersome to implement. The package dynamac is meant to resolve that deficiency by implementing critical value testing procedure for the user. Following Philips [-@philips2016have], it requires estimating the relationship in error-correction form, obtaining statistics of interest, and then testing them via the function `pssbounds`. 

### Cointegration testing using the ARDL-bounds procedure and `pssbounds`

The ARDL-bounds procedure begins with two preliminaries. First, we must ensure that the regressors are not of order I(2) or more and that any seasonal components have been removed. We demonstrated this above when we found that the first-difference of `incshare10` and `urate` were both stationary. In addition, there were no seasonal components, looking at the series (although more formal diagnostics are probably warranted). Second, we must ensure that the dependent variable *is* integrated of order I(1). And again, above, we demonstrated that `incshare10` is integrated.

The next step in ARDL-bounds is to estimate the error-correction form of the model. Two points are key: the independent variables that are potentially I(1) must be entered in levels, and the resulting residuals of the error-correction ARDL model must be white noise (random with no residual autocorrelation). Here, that means that we run our dependent variable in differences, we include the lagged levels of the dependent variable, we include the levels of the potentially cointegrating variable, `incshare10`, and the short-run effects of `urate` through differences. Returning to our model above, and using `dynardl`, this is now straightforward: 

```{r}
res1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        ec = TRUE, simulate = FALSE)
summary(res1)
```

Next we need to ensure that the residuals from this model are white noise. To help with this, we also introduce `dynardl.auto.correlated`. This expects a `dynardl` model and implements a few white-noise-residual tests with readable output that reminds of the null hypotheses. Here, we just run

```{r}
dynardl.auto.correlated(res1)
```

At a reasonable level of significance ($p < 0.10$), the BG test indicates that we reject the null hypothesis of no autocorrelation in the residuals of the model in `res`. Philips [-@philips2016have] indicates that the next step is to add lagged first-differences to our model. To add a lagged difference of $y$, we would run 

```{r}
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = FALSE)
summary(res2)
```

and then again test for autocorrelation.
 
```{r}
dynardl.auto.correlated(res2)
length(res2$model$residuals)
```

We can now be much more confident that there is no more autocorrelation in the residuals. At this point, we can execute the ARDL-bounds test procedure. We need a few pieces of information: the length of the time series, the $t$-statistic on the lagged dependent variable in the ARDL error-correction model, the "case" of the regression (a combination of whether the intercept, if any, and trend, if any, were restricted [@pesaran2001bounds]), and the number of regressors $k$ appearing in first lags (NOT including the lagged dependent variable). We also need the number of observations in the model time series and the $F$-statistic on the restriction that the first lags of each of the variables in the model (including the lagged dependent variable) are jointly equal to zero.

From our regression, we know the $t$-statistic on the lagged dependent variable is -3.684 (just looking at the output). Additionally, we estimated a model with an unrestricted intercept and no trend, which happens to be case 3 (which we would know by referencing Pesaran, Shin, and Smith [-@pesaran2001bounds]). There are $k = 1$ variables in first lags (`incshare10`, not including the LDV), and the number of `obs` in the model is 47 periods. We now only need the test of the restriction that the first lags are equal to zero. We can calculate this by hand through coefficient testing. If we have all of the coefficients, we just need to compare them to the set that are lagged levels:

```{r}
coef(res2$model)
# The second coefficient is the LDV, and the last coefficient is also a first lag. So they're both restricted

B <- coef(res2$model)
V <- vcov(res2$model)

# tag restrictions on LDV and l.1.incshare10
R <- matrix(c(0, 1, 0, 0, 0, 0, 
			  0, 0, 0, 0, 0, 1), byrow = T, nrow = 2)
k.plus1 <- sum(R)

# Restriction is that it is equal to 0
q <- 0
fstat <- (1/k.plus1)*t(R%*%B-q)%*%solve(R%*%V%*%t(R))%*%(R%*%B-q)	
fstat
```

To test for cointegration, we need special critical values for the F-test statistic. For this, we're ready for `pssbounds`

```{r}
pssbounds(obs = 47, fstat = 12.20351, tstat = -3.684, case = 3, k = 1)
```

Finally, a payoff. The $t$-statistic and $F$-statistic are situated between special sample critical values. Depending on the level of significance that we preselected, these values offer a number of different conclusions. If the $F$-statistic exceeds the upper I(1) critical value, we may conclude cointegration. Thus, we can be confident that we have a well-specified model. If the F-statistic falls below the I(0) critical value, Pesaran, Shin, and Smith [-@pesaran2001bounds] note that this indicates that all regressors are in fact stationary. Thus, cointegration cannot exist. If the F-statistic falls between the lower I(0) and upper I(1) critical values, the results are inconclusive. This means that cointegration may exist, but further testing and re-specification of the model is needed. Users that end up with these results are encouraged to see Philips [-@philips2016have], who provides a strategy for dealing with this.

In our model here, since the value of the $F$-statistic exceeds the critical value at the upper I(1) bound of the test at the 5% level, we may conclude that Income Top 10 (the variable in the model in levels) and the dependent variable, Concern, are in a cointegrating relationship. This is furthered by the $t$-statistic of -3.684 falling below the 5% critical value for I(1). Thus, we can be confident both in our model specification and the unique, long-run relationship that exists between the two variables.

Calculating all of those values by hand was unnecessarily involved, especially since we know most of the values are saved post-estimation anyway. We're further motivated not to have to reference the tables in Pesaran, Shin, and Smith [-@pesaran2001bounds] for the case of the regression every time we run `pssbounds`. Thus, our preferred way of estimating the ARDL-bounds test is just by passing the error-correction model to `pssbounds`. In other words, the test is equivalent by just running

```{r}
pssbounds(res2)
```

Pesaran, Shin, and Smith [-@pesaran2001bounds] also provide for testing that the intercept is also restricted, which they refer to as case 2. If we wanted to add this restriction to the test, we would use `restriction = TRUE` in `pssbounds`.

```{r}
pssbounds(res2, restriction = TRUE)
```

Of course, users are encouraged to read Pesaran, Shin, and Smith [-@pesaran2001bounds] and Philips [-@philips2016have] for the implications of this test. We merely note that cases 2 and 4 accessible through `restriction = TRUE`. We also note that, for most users most of the time, case 3 will be the most common case (an unrestricted intercept and no trend).

Any `dynardl` object run in an error-correction format is available for `pssbounds` post-estimation. This, of course, is not meant to imply that all `dynardl` models are meant to be tested for cointegration. Stationary dependent variables cannot be cointegrated with other variables. As such, we would want to run those models in levels, with `EC = FALSE`. Like always, we point users to the importance of pre-regression testing of the nature of our variables, most specifically whether or not they are integrated. `dynardl` will attempt to be as helpful as possible, but in the end, it is up to the user to know which model is appropriate and to ensure that the model is specified correctly.

We had another motivation. `dynardl` + `pssbounds` are a powerful combination for estimating and testing for cointegration in a unified framework. But once we have tested for cointegration, we still want inferences about the nature of the effect of some $x$ variable on the dependent variable. And these inferences become much more difficult to obtain as our models increase in complexity. For instance, the very lagged difference of $y$ that we introduced to help purge autocorrelation from the model we just estimated made interpreting the coefficients from that model much less straightforward. In the next section, we elaborate on the ability of `dynardl` to provide these inferences through the process we outlined earlier: counterfactual simulation of a response to a shock in some $x$ variable. These inferences do not require anything beyond what we have already covered: a `dynardl` model and a theoretically interesting $x$ variable.

### Counterfactual simulation using `dynardl`

ARDL models are useful because they are flexible, but their flexibility undermines our ability to make sense of the coefficients estimated in any given model. An alternative approach to interpretation is to use the coefficients from an estimated model to simulate meaningful responses in the dependent variable to counterfactual changes in an independent variable $x$ (that we control), allowing the change to filter through the various forms of the $x$ variable in the model, as well as different forms of the $y$ variable (like differences and lagged levels) that might be included. 

dynamac handles this simulated interpretation through a function we have already seen: `dynardl`. All we need to do is specify that `simulate = TRUE`, and we will produce simulated responses: we can observe how a change to a variable in a `dynardl` model produces a corresponding change in the dependent variable. Other arguments are required, but only one has no default: we need to know which $x$ variable to simulate a shock to (the `shockvar`). This "shock" means that, at a time specified by the user, the value of an $x$ variable will move to some level. If the variable is in levels or lagged levels, this means that its new value becomes the pre-shock average plus whatever the shock value is. If the variable is in differences or lagged differences, the shock lasts for one period (as a permanent change in a differenced variable would imply that it is changing every period!).

`dynardl` has reasonable defaults for all of the other required parameters: simulations default to a `range` of 20 periods, lines are drawn at roughly `sig` = 95% confidence, the `shockvar` is shocked by a standard deviation (the `shockval`), we use 1,000 `sims` to simulate the average response, we don't force the other $x$ variables to any values (we allow them to take their means, except for differenced variables, which we set to be zero: assuming that, period-over-period, there is no change), we allow for 20 periods of `burnin`, and we create predicted values of the dependent variable, rather than expected values. All of these options can be controlled by the user, but we'll return to that in a moment.

The simulation process is fairly straightforward. `dynardl` uses the coefficients from the model. It draws a number (specifically, `sims`) of values of the coefficients $\hat\beta$ from the estimated regression model. The distribution is assumed to be multivariate normal with mean of $\hat\beta$ and variance-covariance of the estimated model. Uncertainty is incorporated by simulating $\hat\sigma$ $^{2}$ as a scaled draw from the chi-squared distribution. These fit together by using values of the independent variables (usually their means: see the preceding paragraph) $X$, multiplying by $\hat\beta$ to get predicted values of the dependent variable $y$, and then using $\hat \sigma$ $^{2}$ to introduce uncertainty in the predicted values. If you want to exert more direct control over the values of any $x$ variables used, you can `forceset` them to any other value you like. This is not limited to means or integers; if you have any substantively interesting value you'd like to hold a variable at, you are free to specify whichever value you like. 

But what we're really interested in is the effect of some variable $x$. In the world of `dynardl`, this is our `shockvar`. We specify a variable from the model, tell `dynardl` how much we want it to theoretically change by (in `shockval`, defaulting to a standard deviation of the `shockvar`) and when we want it to change (`time`), and then observe the change in $y$.

Let's go back to our earlier model. Remember we ran

```{r, eval = FALSE}
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = FALSE)
summary(res2)
```

Here, we set `simulate = FALSE` because we were just illustrating estimation without simulations (given that simulations take a few seconds to estimate, we may only want to simulate changes when we have a final model, free of autocorrelation and the other things we tested for). Now, we'll turn `simulate = TRUE` and specify an $x$ variable to observe the response of. We'll observe the changes resulting from `incshare10`, given its lagged level demonstrated a significant effect. In other words, our `shockvar = "incshare10"`. By default, `dynardl` will shock it by its standard deviation, but we can observe any change we like with `shockval`. As with any other time in which stochastic values are involved, we should set a seed to ensure our results are directly replicable.

```{r}
set.seed(020990)
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE,
        shockvar = "incshare10")
summary(res2$model)
```

It looks somewhat goofy in the vignette, as `dynardl` displays a progress bar to inform you of how many simulations have been completed. But the payoff is in `res2$simulation`. We get to observe the effect of a standard-deviation shock in `incshare10` in the dependent variable. This response is best visualized in a plot.  To see this plot, the model can be saved to an object, and plots can be created with `dynardl.simulation.plot`.

```{r}
dynardl.simulation.plot(res2, type = "area", response = "levels")
```

Here, the shock has an immediate effect that does not dissipate for a long time. So long, in fact, that we might want to lengthen the `range` of the simulation. 

```{r}
set.seed(020990)
res3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10")
dynardl.simulation.plot(res3, type = "area", response = "levels")
```

If instead of an area plot we desired a spike plot:

```{r}
dynardl.simulation.plot(res3, type = "spike", response = "levels")
```

In `res3` there are actually three sets of output. The first is the model, the second is the information for `pssbounds`, and the last is for plotting. We mention this so that users can use their usual TeX shortcuts for creating tables, customizing plots, or testing using pssbounds later.

```{r}
res3$model
res3$pssbounds
res3$simulation
```

Only models where `ec = TRUE` will have `$pssbounds` information, as only those models can possibly be testing a cointegrating relationship. We suspect that no one will need them by hand, as you can just pass the whole object via `pssbounds(res3)` to get the same results. 

Other options might be of interest. In order to smooth our confidence intervals (note: this does NOT make us more confident), we can increase the number of simulations. Additionally, the plotting function allows the user to produce plot in grayscale (for publications). Consider the following (notice the `sims` argument of `dynardl` and the new arguments under `dynardl.simulation.plot`:

```{r}
set.seed(020990)
res4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30, sims = 10000,
        shockvar = "incshare10")
dynardl.simulation.plot(res4, type = "spike", response = "levels", bw = TRUE)
```

The full extent of these options is addressed in the documentation.

There is some question as to whether or not quantities of interest from these types of stochastic simulations are best summarized by the means or the medians of the resulting distributions. In most cases, the difference is likely to be minimal. In cases of skew, though, the median might serve significantly better than the mean (described in Rainey [-@raineytransform]). Here, we can do that by setting `qoi = median`. 

```{r}
set.seed(020990)
res5 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10", qoi = "median")
dynardl.simulation.plot(res5, type = "area", response = "levels")
```

There are a variety of quantities that are plottable from the simulations, outside of just the response of the dependent variable. `dynardl.simulation.plot` includes options for plotting the `levels` of the dependent variable, the levels but demeaned (`levels.from.mean`) of the dependent variable, and the period-over-period `diffs` of the changes in the simulated values. You can get a sense of how the shock to the independent variable is decaying through time by observing the differences in each period as an absolute value (how much is the dependent variable adjusting) through `shock.effect.decay.` You can also see the `cumulative.diffs` of those changes and the absolute value of the changes `cumulative.abs.diffs`. For the final two options, `fullsims = TRUE` must be specified in the `dynardl` simulation, which reports the full raw simulated values as a part of the `dynardl` object. These values are used to draw confidence intervals over the simulated changes. 

In addition, `dynardl.simulation.plot` will expect a value for when it should regard the changes in the dependent variable as noise, rather than *real* changes. The default `tol` is to disregard changes that are less than 0.1% of the mean of the dependent variable. Alternatively, we can specify a `last.period` in which we believe the dependent variable is responding. `dynardl.simulation.plot` also allows you to pass whatever normal plotting arguments you would use in the normal `...` way.

```{r}
set.seed(020990)
res6 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10", fullsims = TRUE)

par(mfrow = c(2, 3))
dynardl.simulation.plot(res6, type = "area", response = "levels")
dynardl.simulation.plot(res6, type = "area", response = "levels.from.mean")
dynardl.simulation.plot(res6, type = "area", response = "diffs")
dynardl.simulation.plot(res6, type = "area", response = "shock.effect.decay")
dynardl.simulation.plot(res6, type = "area", response = "cumulative.diffs", axes = F)
dynardl.simulation.plot(res6, type = "area", response = "cumulative.abs.diffs")
```

If we don't want to see the equilibriating behavior before the shock, we can draw the plot starting in a different period through `start.period.`

```{r}
dynardl.simulation.plot(res6, response = "shock.effect.decay", start.period = 9)
```

If you're in an exploratory phase of model building and prefer not to have to run the plots separately, you are free to combine all of them using `dynardl.all.plots`. 

```{r}
dynardl.all.plots(res6)
```

## Being smart(er than `dynardl`) about data and modeling

dynamac is meant to be a unifying package for estimating and interpreting times series ARDL models. It is not, however, a data manager. It assumes that your data are ordered, that the progression between time series makes sense (i.e. there is a consistent unit of time separating the ordered observations), that there are not problematic missing values, and that all the other usual pre-estimation data testing has been performed by the user. Users should know and explore their datasets well before passing those data to any statistical software.

Nor will dynamac stop you from running "bad idea" models. Users should be careful about the order of integration in their variables, whether seasonal unit roots are present, if variables have nuanced lagged-difference structures, and so on. We offer functions (like `dynardl.auto.correlated`) to help users make these decisions on their path to a final model. But care must be taken at every step.

## Bibliography


---
title: "An Introduction to dynamac: Dynamic Inferences (and Cointegration Testing) from Autoregressive Distributed Lag Models"
author: "Soren Jordan and Andrew Q. Philips"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: timeseriesrj.bib
vignette: >
  %\VignetteIndexEntry{dynamac Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Autoregressive distributed lag (ARDL) models are an integral part of estimating scientific processes over time. However, as we extend their usefulness by adding richness in dynamic specifications (through multiple lags of variables, either in levels or differences, or lags of the dependent variable), we begin to challenge our ability to draw meaningful inferences from coefficients alone. Variables may appear in multiple time periods, in multiple forms, and might even be filtered through lagged values of the dependent variable. Coefficients tell us about the immediate effect of some variable but have little to say about the long-run effect. 

There is a better solution. Instead of performing intense operations to develop a closed-form or algebraic solution, we can rely on the power of computing to *simulate* the over-time response in the dependent variable of some model, given a corresponding change in an $x$ variable. We often call these "changes" in $x$ variables, and the associated response in the dependent variable $y$, a "counterfactual response" (a simulated response to a "shock" we control).

dynamac helps simulate these counterfactuals. More generally, though, it is built to make using and drawing inferences from single-equation ARDL models as easy as possible. Moreover, it helps users implement the useful cointegration test from Pearson, Shin, and Smith [-@pesaran2001bounds]: the ARDL-bounds testing procedure.

We illustrate the usefulness of these functions through example. After a brief discussion of the ARDL model in the general sense, we estimate a collection of these models and demonstrate both the challenges of the ARDL model in the abstract and the solutions of dynamac in particular. 

## ARDL models generally

The ARDL model has a general form where $y$, modeled in levels or differences, is a function of itself (in lagged levels or differences), up to $k$ variables $x$, either in contemporaneous (same period, or appearing at time $t$) levels, lagged levels, contemporaneous differences, or lagged differences. Conventionally, the number of lags of the dependent variable in levels is counted by $p$, the number of lags of the dependent variable in differences is counted by $m$, the number of lags of the independent variables in levels is counted by $l$, and the number of lags of the independent variables in differences is counted by $q$ (note: $l$ and $q$ can begin at 0, i.e, contemporaneous effects appearing at time $t$). 

The number of time periods of each variable can, theoretically, be quite large. Or, put differently, $p$, $m$, $l$, and $q$, especially $l$ and $q$, might be hard to account for. Analysts without strong theory about the nature of the effects often begin with restricting all but the contemporaneous and first lag of each series, or an ARDL model of the nature 

$$y_t =  \alpha_0 + \phi_1 y_{t-1} + \theta_{1,0} x_{1,t} + \theta_{1,1} x_{1,t-1} + \cdots + \theta_{k,0}x_{k,t} + \theta_{k,1}x_{k,t-1}+ \beta *T  + \epsilon_t$$
where $\alpha_0$ is a constant and $\beta *T$ is a trend term. (The exact nature of this equation depends on whether $y$ is stationary or integrated, as well as if differences or lagged differences are entered into the model. But we will get to this below.)

It's useful at this point to stop and think: if I have multiple lags of a variable $x$, and they are filtered through a lagged dependent variable $y_{t-1}$, it might be difficult to get a sense of the "total" effect of $x$ on $y$. This becomes more and more difficult as $x$ increases in lagged levels $l$ *and potentially also* lagged differences $q$. Our coefficient estimates, while estimated, are no longer as useful in direct interpretation. Put differently, the very flexibility of the ARDL model also undoes its usefulness in interpretation! So we might seek an alternative way of interpreting these models. dynamac provides a unified way of estimating ARDL models and interpreting their effects. It also provides a way of implementing a popular test for cointegration. We uncover these methods below. 

## Estimating an ARDL model: best practices and dynamac

dynamac includes two datasets. We will use the Wright [-@wrightpolitical] dataset on income inequality. We can look at the first few rows of this dataset

```{r}
library(dynamac)
head(ineq)
```

`concern` is public concern about income inequality; `incshare10` is the income share of the top ten percent; `urate` is the unemployment rate. Wright [-@wrightpolitical] argues that concern about income inequality grows as inequality itself worsens and economic conditions deteriorate. A simple model, then, is that concern is a function of past values of itself, the level of income share held by the top ten percent, *changing* levels of income share held by the top ten percent, as well *changing* levels of unemployment in the short term. 

$$\Delta Concern_t = \alpha_0 + \phi_1 Concern_{t-1} + \theta_1 Income Top 10_{t-1} + \beta_1 \Delta Income Top 10_t + \beta_2 \Delta Unemployment_t + \epsilon_t$$

where the residuals are white noise. Let's develop the corresponding ARDL model using dynamac.

### Understanding our time series

Step 1 in any time series analysis is a visual inspection of the series coupled with formal tests of stationarity: whether the series has a constant mean, variance, and covariance over time (so that it reverts back to mean level), or if the series instead violates any of these three conditions. We advocate for using the urca package for this [@pfaff2016package], which includes a variety of tools and tests. Simply plotting the series reveals the following:

```{r}
library(urca)
ts.plot(ineq$concern)
ts.plot(ineq$incshare10)
ts.plot(ineq$urate)
```

None of the three series looks especially well-behaved. `concern` rose quickly and has been moving sluggishly since. `incshare10` has only grown over time, which cannot be mean-reverting (like a stationary series). `urate` experiences steep peaks and valleys with significant interludes in between. All three series look integrated. But we should be more discerning by using empirical tests, also included in urca. So we can execute the Augmented Dickey-Fuller (ADF) test, Phillips-Perron test, DF-GLS test, and KPSS test on each series. On `concern` this looks like

```{r}
summary(ur.df(ineq$concern, type = c("none"), lags = 1))
summary(ur.pp(ineq$concern, type = c("Z-tau"), model = c("constant"), use.lag = 1))
summary(ur.ers(ineq$concern, type = c("DF-GLS"), model = c("constant"), lag.max = 1)) 
summary(ur.kpss(ineq$concern, type = c("mu"), use.lag = 1))
```

The ADF test, Phillips-Perron test, and DF-GLS test all have null hypotheses of a unit root, all of which we fail to reject. The KPSS test has a null hypothesis of stationarity which we *do* reject. We have compelling evidence, then, of integration (I[1]) in the series `concern`. We can check whether differencing is enough to make the series stationary by executing the same tests

```{r eval = FALSE}
summary(ur.df(diff(ineq$concern), type = c("none"), lags = 1))
summary(ur.pp(diff(ineq$concern), type = c("Z-tau"), model = c("constant"), use.lag = 1))
summary(ur.ers(diff(ineq$concern), type = c("DF-GLS"), model = c("constant"), lag.max = 1)) 
summary(ur.kpss(diff(ineq$concern), type = c("mu"), use.lag = 1))
```

Each of the above tests, when run, indicated that the differenced series of `concern` is stationary. Having gathered the basic information about the nature of the history of our variables, we might be itching to estimate some preliminary models. To this point, R hasn't offered much in the way of improving beyond the basic `lm()` framework for using regression to estimate time series models in a consistent syntax. We elaborate on this problem and illustrate our solution, `dynardl`.

### Estimating ARDL models with `dynardl` 

ARDL models are flexible, but their flexibility often results in variables of different lengths due to differencing and lagging. For instance, consider our simple model from above where `incshare10` appears as a first lag and a first difference, `urate` appears as a first difference, there is a lagged dependent variable, and `concern` is the dependent variable in differences. We can introduce a lag through the built-in `lshift` function in dynamac. The syntax is just `lshift(variable, num.lags)`, where `num.lags` is the number of periods for the variable to be lagged. We can also difference a series through `dshift`. The syntax is just `dshift(variable)` for a first difference. For instance, 

```{r}
head(ineq$incshare10)
head(lshift(ineq$incshare10, 1))
head(dshift(ineq$incshare10))
```

So the syntax for the simple model described would be easy to write. But notice the problem

```{r, eval = TRUE, error = TRUE}
summary(lm(diff(ineq$concern) ~ lshift(ineq$concern, 1) + lshift(ineq$incshare10, 1) + dshift(ineq$incshare10) + dshift(ineq$urate)))
```

Introducing the lag changed the variable lengths, and we probably don't want to have to think about time series operation in our model specification, anyway. Other software has unified this operation by introducing a common syntax, like d. for differencing and l. for lagging, or even l.d. for lag differences (what program could that be?). In R, though, it remains more of a nuisance than it should be. The `dynardl` function helps to remedy this challenge by encouraging the user only to think about model structure in exactly the way outlined above: which variables $x$ get lags $l$ and lagged differences $q$, etc. The function expects a basic formula of all the variables in the model, the data set (`data`), if there is one, then lagged levels (`lags`) and lagged differences (`lagdiffs`) as a list and differences (`diffs`) and contemporaneous levels (`levels`) as a vector. The sole exception is if the user wants to run the model in differences, s/he needs to specify `ec = TRUE` (in an effort to force us to think critically about the error-correction form of the model). For instance, our above model now becomes (ignoring the `simulate` option for now)

```{r, eval = FALSE}
dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        ec = TRUE, simulate = FALSE)
```

Purely hypothetically, if we wanted more lagged levels of `concern`, we would just change `1` to `c(1:5)` for lags at $t-1$ to $t-5$, or any other number of lags. If we wanted to include the first-difference of `urate` lagged at periods $t-1$ and $t-2$, we would now run

```{r, eval = FALSE}
dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("urate" = c(1, 2)),
        ec = TRUE, simulate = FALSE)
```

If the variable can appear in multiple periods (`lags` and `lagdiffs`), it must be specified as a list. If it can only appear contemporaneously (`levels` and `diffs`), it must be specified as a vector. (The alternative was to specify `levels` through a `lag` at time 0: this did not seem practical.) We can also add or suppress a constant from the model by specifying `constant = TRUE/FALSE`, and we can do the same with a linear trend term with `trend = TRUE/FALSE` (by default, constants are included and trends are not). 

The option `ec` specifies the nature of the dependent variable. If it is possibly error-correcting (`ec = TRUE`), it is run in differences, or period-over-period changes. If it is not (`ec = FALSE`), the dependent variable is run in levels.

At this point, `dynardl` is just a unifying way of thinking about time series models for estimation. Yet it is going to offer two other important advantages: interpreting the effects of our variables, and executing a powerful test for cointegration. We will start with the latter. Remember testing for I(1) processes earlier: each test indicated that `concern` was I(1). Running the same tests for each of `incshare10` and `urate` suggests that all three series are I(1). This might cause us concern about cointegration: the special relationship between I(1) series where the series are in a long-run relationship, even if they move apart in the short term.  Cointegrating series are difficult to *uncover*. Traditional methods, like the Engle and Granger [-@engle1987co] two-step method or likelihood-based approaches [@johansen1995likelihood] too often conclude cointegration when it does not exist, at least in smaller sample sizes [@philips2016have]. An alternative test [@pesaran2001bounds], which we refer to as the ARDL-bounds procedure, performs much better in small samples ($t < 80$), but it is more cumbersome to implement. The package dynamac is meant to resolve that deficiency by implementing critical value testing procedure for the user. Following Philips [-@philips2016have], it requires estimating the relationship in error-correction form, obtaining statistics of interest, and then testing them via the function `pssbounds`. 

### Cointegration testing using the ARDL-bounds procedure and `pssbounds`

The ARDL-bounds procedure begins with two preliminaries. First, we must ensure that the regressors are not of order I(2) or more and that any seasonal components have been removed. We demonstrated this above when we found that the first-difference of `incshare10` and `urate` were both stationary. In addition, there were no seasonal components, looking at the series (although more formal diagnostics are probably warranted). Second, we must ensure that the dependent variable *is* integrated of order I(1). And again, above, we demonstrated that `incshare10` is integrated.

The next step in ARDL-bounds is to estimate the error-correction form of the model. Two points are key: the independent variables that are potentially I(1) must be entered in levels, and the resulting residuals of the error-correction ARDL model must be white noise (random with no residual autocorrelation). Here, that means that we run our dependent variable in differences, we include the lagged levels of the dependent variable, we include the levels of the potentially cointegrating variable, `incshare10`, and the short-run effects of `urate` through differences. Returning to our model above, and using `dynardl`, this is now straightforward: 

```{r}
res1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        ec = TRUE, simulate = FALSE)
summary(res1)
```

Next we need to ensure that the residuals from this model are white noise. To help with this, we also introduce `dynardl.auto.correlated`. This expects a `dynardl` model and implements a few white-noise-residual tests with readable output that reminds of the null hypotheses. Here, we just run

```{r}
dynardl.auto.correlated(res1)
```

At a reasonable level of significance ($p < 0.10$), the BG test indicates that we reject the null hypothesis of no autocorrelation in the residuals of the model in `res`. Philips [-@philips2016have] indicates that the next step is to add lagged first-differences to our model. To add a lagged difference of $y$, we would run 

```{r}
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = FALSE)
summary(res2)
```

and then again test for autocorrelation.
 
```{r}
dynardl.auto.correlated(res2)
length(res2$model$residuals)
```

We can now be much more confident that there is no more autocorrelation in the residuals. At this point, we can execute the ARDL-bounds test procedure. We need a few pieces of information: the length of the time series, the $t$-statistic on the lagged dependent variable in the ARDL error-correction model, the "case" of the regression (a combination of whether the intercept, if any, and trend, if any, were restricted [@pesaran2001bounds]), and the number of regressors $k$ appearing in first lags (NOT including the lagged dependent variable). We also need the number of observations in the model time series and the $F$-statistic on the restriction that the first lags of each of the variables in the model (including the lagged dependent variable) are jointly equal to zero.

From our regression, we know the $t$-statistic on the lagged dependent variable is -3.684 (just looking at the output). Additionally, we estimated a model with an unrestricted intercept and no trend, which happens to be case 3 (which we would know by referencing Pesaran, Shin, and Smith [-@pesaran2001bounds]). There are $k = 1$ variables in first lags (`incshare10`, not including the LDV), and the number of `obs` in the model is 47 periods. We now only need the test of the restriction that the first lags are equal to zero. We can calculate this by hand through coefficient testing. If we have all of the coefficients, we just need to compare them to the set that are lagged levels:

```{r}
coef(res2$model)
# The second coefficient is the LDV, and the last coefficient is also a first lag. So they're both restricted

B <- coef(res2$model)
V <- vcov(res2$model)

# tag restrictions on LDV and l.1.incshare10
R <- matrix(c(0, 1, 0, 0, 0, 0, 
			  0, 0, 0, 0, 0, 1), byrow = T, nrow = 2)
k.plus1 <- sum(R)

# Restriction is that it is equal to 0
q <- 0
fstat <- (1/k.plus1)*t(R%*%B-q)%*%solve(R%*%V%*%t(R))%*%(R%*%B-q)	
fstat
```

To test for cointegration, we need special critical values for the F-test statistic. For this, we're ready for `pssbounds`

```{r}
pssbounds(obs = 47, fstat = 12.20351, tstat = -3.684, case = 3, k = 1)
```

Finally, a payoff. The $t$-statistic and $F$-statistic are situated between special sample critical values. Depending on the level of significance that we preselected, these values offer a number of different conclusions. If the $F$-statistic exceeds the upper I(1) critical value, we may conclude cointegration. Thus, we can be confident that we have a well-specified model. If the F-statistic falls below the I(0) critical value, Pesaran, Shin, and Smith [-@pesaran2001bounds] note that this indicates that all regressors are in fact stationary. Thus, cointegration cannot exist. If the F-statistic falls between the lower I(0) and upper I(1) critical values, the results are inconclusive. This means that cointegration may exist, but further testing and re-specification of the model is needed. Users that end up with these results are encouraged to see Philips [-@philips2016have], who provides a strategy for dealing with this.

In our model here, since the value of the $F$-statistic exceeds the critical value at the upper I(1) bound of the test at the 5% level, we may conclude that Income Top 10 (the variable in the model in levels) and the dependent variable, Concern, are in a cointegrating relationship. This is furthered by the $t$-statistic of -3.684 falling below the 5% critical value for I(1). Thus, we can be confident both in our model specification and the unique, long-run relationship that exists between the two variables.

Calculating all of those values by hand was unnecessarily involved, especially since we know most of the values are saved post-estimation anyway. We're further motivated not to have to reference the tables in Pesaran, Shin, and Smith [-@pesaran2001bounds] for the case of the regression every time we run `pssbounds`. Thus, our preferred way of estimating the ARDL-bounds test is just by passing the error-correction model to `pssbounds`. In other words, the test is equivalent by just running

```{r}
pssbounds(res2)
```

Pesaran, Shin, and Smith [-@pesaran2001bounds] also provide for testing that the intercept is also restricted, which they refer to as case 2. If we wanted to add this restriction to the test, we would use `restriction = TRUE` in `pssbounds`.

```{r}
pssbounds(res2, restriction = TRUE)
```

Of course, users are encouraged to read Pesaran, Shin, and Smith [-@pesaran2001bounds] and Philips [-@philips2016have] for the implications of this test. We merely note that cases 2 and 4 accessible through `restriction = TRUE`. We also note that, for most users most of the time, case 3 will be the most common case (an unrestricted intercept and no trend).

Any `dynardl` object run in an error-correction format is available for `pssbounds` post-estimation. This, of course, is not meant to imply that all `dynardl` models are meant to be tested for cointegration. Stationary dependent variables cannot be cointegrated with other variables. As such, we would want to run those models in levels, with `EC = FALSE`. Like always, we point users to the importance of pre-regression testing of the nature of our variables, most specifically whether or not they are integrated. `dynardl` will attempt to be as helpful as possible, but in the end, it is up to the user to know which model is appropriate and to ensure that the model is specified correctly.

We had another motivation. `dynardl` + `pssbounds` are a powerful combination for estimating and testing for cointegration in a unified framework. But once we have tested for cointegration, we still want inferences about the nature of the effect of some $x$ variable on the dependent variable. And these inferences become much more difficult to obtain as our models increase in complexity. For instance, the very lagged difference of $y$ that we introduced to help purge autocorrelation from the model we just estimated made interpreting the coefficients from that model much less straightforward. In the next section, we elaborate on the ability of `dynardl` to provide these inferences through the process we outlined earlier: counterfactual simulation of a response to a shock in some $x$ variable. These inferences do not require anything beyond what we have already covered: a `dynardl` model and a theoretically interesting $x$ variable.

### Counterfactual simulation using `dynardl`

ARDL models are useful because they are flexible, but their flexibility undermines our ability to make sense of the coefficients estimated in any given model. An alternative approach to interpretation is to use the coefficients from an estimated model to simulate meaningful responses in the dependent variable to counterfactual changes in an independent variable $x$ (that we control), allowing the change to filter through the various forms of the $x$ variable in the model, as well as different forms of the $y$ variable (like differences and lagged levels) that might be included. 

dynamac handles this simulated interpretation through a function we have already seen: `dynardl`. All we need to do is specify that `simulate = TRUE`, and we will produce simulated responses: we can observe how a change to a variable in a `dynardl` model produces a corresponding change in the dependent variable. Other arguments are required, but only one has no default: we need to know which $x$ variable to simulate a shock to (the `shockvar`). This "shock" means that, at a time specified by the user, the value of an $x$ variable will move to some level. If the variable is in levels or lagged levels, this means that its new value becomes the pre-shock average plus whatever the shock value is. If the variable is in differences or lagged differences, the shock lasts for one period (as a permanent change in a differenced variable would imply that it is changing every period!).

`dynardl` has reasonable defaults for all of the other required parameters: simulations default to a `range` of 20 periods, lines are drawn at roughly `sig` = 95% confidence, the `shockvar` is shocked by a standard deviation (the `shockval`), we use 1,000 `sims` to simulate the average response, we don't force the other $x$ variables to any values (we allow them to take their means, except for differenced variables, which we set to be zero: assuming that, period-over-period, there is no change), we allow for 20 periods of `burnin`, and we create predicted values of the dependent variable, rather than expected values. All of these options can be controlled by the user, but we'll return to that in a moment.

The simulation process is fairly straightforward. `dynardl` uses the coefficients from the model. It draws a number (specifically, `sims`) of values of the coefficients $\hat\beta$ from the estimated regression model. The distribution is assumed to be multivariate normal with mean of $\hat\beta$ and variance-covariance of the estimated model. Uncertainty is incorporated by simulating $\hat\sigma$ $^{2}$ as a scaled draw from the chi-squared distribution. These fit together by using values of the independent variables (usually their means: see the preceding paragraph) $X$, multiplying by $\hat\beta$ to get predicted values of the dependent variable $y$, and then using $\hat \sigma$ $^{2}$ to introduce uncertainty in the predicted values. If you want to exert more direct control over the values of any $x$ variables used, you can `forceset` them to any other value you like. This is not limited to means or integers; if you have any substantively interesting value you'd like to hold a variable at, you are free to specify whichever value you like. 

But what we're really interested in is the effect of some variable $x$. In the world of `dynardl`, this is our `shockvar`. We specify a variable from the model, tell `dynardl` how much we want it to theoretically change by (in `shockval`, defaulting to a standard deviation of the `shockvar`) and when we want it to change (`time`), and then observe the change in $y$.

Let's go back to our earlier model. Remember we ran

```{r, eval = FALSE}
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = FALSE)
summary(res2)
```

Here, we set `simulate = FALSE` because we were just illustrating estimation without simulations (given that simulations take a few seconds to estimate, we may only want to simulate changes when we have a final model, free of autocorrelation and the other things we tested for). Now, we'll turn `simulate = TRUE` and specify an $x$ variable to observe the response of. We'll observe the changes resulting from `incshare10`, given its lagged level demonstrated a significant effect. In other words, our `shockvar = "incshare10"`. By default, `dynardl` will shock it by its standard deviation, but we can observe any change we like with `shockval`. As with any other time in which stochastic values are involved, we should set a seed to ensure our results are directly replicable.

```{r}
set.seed(020990)
res2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE,
        shockvar = "incshare10")
summary(res2$model)
```

It looks somewhat goofy in the vignette, as `dynardl` displays a progress bar to inform you of how many simulations have been completed. But the payoff is in `res2$simulation`. We get to observe the effect of a standard-deviation shock in `incshare10` in the dependent variable. This response is best visualized in a plot.  To see this plot, the model can be saved to an object, and plots can be created with `dynardl.simulation.plot`.

```{r}
dynardl.simulation.plot(res2, type = "area", response = "levels")
```

Here, the shock has an immediate effect that does not dissipate for a long time. So long, in fact, that we might want to lengthen the `range` of the simulation. 

```{r}
set.seed(020990)
res3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10")
dynardl.simulation.plot(res3, type = "area", response = "levels")
```

If instead of an area plot we desired a spike plot:

```{r}
dynardl.simulation.plot(res3, type = "spike", response = "levels")
```

In `res3` there are actually three sets of output. The first is the model, the second is the information for `pssbounds`, and the last is for plotting. We mention this so that users can use their usual TeX shortcuts for creating tables, customizing plots, or testing using pssbounds later.

```{r}
res3$model
res3$pssbounds
res3$simulation
```

Only models where `ec = TRUE` will have `$pssbounds` information, as only those models can possibly be testing a cointegrating relationship. We suspect that no one will need them by hand, as you can just pass the whole object via `pssbounds(res3)` to get the same results. 

Other options might be of interest. In order to smooth our confidence intervals (note: this does NOT make us more confident), we can increase the number of simulations. Additionally, the plotting function allows the user to produce plot in grayscale (for publications). Consider the following (notice the `sims` argument of `dynardl` and the new arguments under `dynardl.simulation.plot`:

```{r}
set.seed(020990)
res4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30, sims = 10000,
        shockvar = "incshare10")
dynardl.simulation.plot(res4, type = "spike", response = "levels", bw = TRUE)
```

The full extent of these options is addressed in the documentation.

There is some question as to whether or not quantities of interest from these types of stochastic simulations are best summarized by the means or the medians of the resulting distributions. In most cases, the difference is likely to be minimal. In cases of skew, though, the median might serve significantly better than the mean (described in Rainey [-@raineytransform]). Here, we can do that by setting `qoi = median`. 

```{r}
set.seed(020990)
res5 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10", qoi = "median")
dynardl.simulation.plot(res5, type = "area", response = "levels")
```

There are a variety of quantities that are plottable from the simulations, outside of just the response of the dependent variable. `dynardl.simulation.plot` includes options for plotting the `levels` of the dependent variable, the levels but demeaned (`levels.from.mean`) of the dependent variable, and the period-over-period `diffs` of the changes in the simulated values. You can get a sense of how the shock to the independent variable is decaying through time by observing the differences in each period as an absolute value (how much is the dependent variable adjusting) through `shock.effect.decay.` You can also see the `cumulative.diffs` of those changes and the absolute value of the changes `cumulative.abs.diffs`. For the final two options, `fullsims = TRUE` must be specified in the `dynardl` simulation, which reports the full raw simulated values as a part of the `dynardl` object. These values are used to draw confidence intervals over the simulated changes. 

In addition, `dynardl.simulation.plot` will expect a value for when it should regard the changes in the dependent variable as noise, rather than *real* changes. The default `tol` is to disregard changes that are less than 0.1% of the mean of the dependent variable. Alternatively, we can specify a `last.period` in which we believe the dependent variable is responding. `dynardl.simulation.plot` also allows you to pass whatever normal plotting arguments you would use in the normal `...` way.

```{r}
set.seed(020990)
res6 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = TRUE, range = 30,
        shockvar = "incshare10", fullsims = TRUE)

par(mfrow = c(2, 3))
dynardl.simulation.plot(res6, type = "area", response = "levels")
dynardl.simulation.plot(res6, type = "area", response = "levels.from.mean")
dynardl.simulation.plot(res6, type = "area", response = "diffs")
dynardl.simulation.plot(res6, type = "area", response = "shock.effect.decay")
dynardl.simulation.plot(res6, type = "area", response = "cumulative.diffs", axes = F)
dynardl.simulation.plot(res6, type = "area", response = "cumulative.abs.diffs")
```

If we don't want to see the equilibriating behavior before the shock, we can draw the plot starting in a different period through `start.period.`

```{r}
dynardl.simulation.plot(res6, response = "shock.effect.decay", start.period = 9)
```

If you're in an exploratory phase of model building and prefer not to have to run the plots separately, you are free to combine all of them using `dynardl.all.plots`. 

```{r}
dynardl.all.plots(res6)
```

## Being smart(er than `dynardl`) about data and modeling

dynamac is meant to be a unifying package for estimating and interpreting times series ARDL models. It is not, however, a data manager. It assumes that your data are ordered, that the progression between time series makes sense (i.e. there is a consistent unit of time separating the ordered observations), that there are not problematic missing values, and that all the other usual pre-estimation data testing has been performed by the user. Users should know and explore their datasets well before passing those data to any statistical software.

Nor will dynamac stop you from running "bad idea" models. Users should be careful about the order of integration in their variables, whether seasonal unit roots are present, if variables have nuanced lagged-difference structures, and so on. We offer functions (like `dynardl.auto.correlated`) to help users make these decisions on their path to a final model. But care must be taken at every step.

## Bibliography


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\docType{data}
\name{ineq}
\alias{ineq}
\title{Data on public concern about economic inequality}
\format{A data frame with 49 rows and 9 variables:
\describe{
  \item{year}{Year}
  \item{mood}{Public mood liberalism}
  \item{urate}{Unemployment rate}
  \item{concern}{Concern about economic inequality}
  \item{demcontrol}{Democratic control of congress}
  \item{incshare10}{Proportion of income of top 10 percent}
  \item{csentiment}{Consumer sentiment}
  \item{incshare01}{Proportion of income of top 1 percent}
}}
\source{
\url{http://dx.doi.org/10.7910/DVN/UYUU9G}
}
\usage{
data(ineq)
}
\description{
A dataset from: Wright, Graham. 2017. "The political
implications of American concerns about economic inequality."
Political Behavior 40(2): 321-346.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dshift}
\alias{dshift}
\title{Take first difference of a series}
\usage{
dshift(x)
}
\arguments{
\item{x}{a series to be differenced}
}
\value{
the differenced series
}
\description{
Take first difference of a series
}
\details{
\code{dshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
}
\examples{
x.var <- seq(0, 50, 5)
d.x.var <- dshift(x.var)
head(x.var)
head(d.x.var)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dynamac-deprecated}
\alias{dynamac-deprecated}
\alias{area.simulation.plot}
\alias{spike.simulation.plot}
\title{Deprecated functions in package \pkg{dynamac}}
\usage{
area.simulation.plot(x, response = "levels", bw = FALSE)

spike.simulation.plot(x, response = "levels", bw = FALSE)
}
\arguments{
\item{x}{a dynardl model with a simulation to be plotted}

\item{response}{whether the plot of the response should be shown in levels of the dependent variable (\code{levels}) 
or in changes from the mean of the dependent variable (\code{mean.changes}). The default is \code{levels}}

\item{bw}{should the colors be in black and white (for publication)? The default is \code{FALSE}}

\item{x}{a dynardl model with a simulation to be plotted}

\item{response}{whether the plot of the response should be shown in levels of the dependent variable (\code{levels}) 
or in changes from the mean of the dependent variable (\code{mean.changes}). The default is \code{levels}}

\item{bw}{should the colors be in black and white (for publication)? The default is \code{FALSE}}
}
\value{
an area plot

a spike plot
}
\description{
The functions listed below are deprecated and will be defunct
in the near future. When possible, functions with similar functionality are mentioned.
Help pages for deprecated functions are available at \code{help-("deprecated")}.
}
\section{\code{area.simulation.plot}}{

For \code{area.simulation.plot}, use \code{\link{dynardl.simulation.plot}}.
}

\section{\code{spike.simulation.plot}}{

For \code{spike.simulation.plot}, use \code{\link{dynardl.simulation.plot}}.
}

\seealso{
\code{link{dynamac-deprecated}}

\code{link{dynamac-deprecated}}
}
\author{
Soren Jordan and Andrew Q. Philips

Soren Jordan and Andrew Q. Philips
}
\keyword{internal}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{summary.dynardl}
\alias{summary.dynardl}
\title{Enable summary calls to \code{\link{dynardl}} model objects}
\usage{
\method{summary}{dynardl}(object, ...)
}
\arguments{
\item{object}{a \code{dynardl} model}

\item{...}{additional arguments in the generic summary call}
}
\value{
A summary of the fitted ARDL model.
}
\description{
Enable summary calls to \code{\link{dynardl}} model objects
}
\details{
\code{dynardl}, by default, stores regression results in \code{foo$model}. This calls those results directly with \code{summary}
}
\examples{
# Using the ineq data from dynamac
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = FALSE)
summary(ardl.model)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\docType{data}
\name{supreme.sup}
\alias{supreme.sup}
\title{Data on US Supreme Court Approval}
\format{A data frame with 42 rows and 9 variables:
\describe{
  \item{dcalc}{Supreme Court support}
  \item{l_dcalc}{Lagged Supreme Court spport}
  \item{iddiv}{Ideological divergence}
  \item{mooddev}{Mean deviation of Mood}
  \item{dirdev}{Mean deviation of percent liberal decisions}
  \item{sg}{Rulings against Solicitor General's amicus briefs}
  \item{laws}{Laws declared unconstitutional}
  \item{presapp}{Approval of president}
  \item{congapp}{Approval of Congress}
}}
\source{
\url{http://dx.doi.org/10.2307/2669280}
}
\usage{
data(supreme.sup)
}
\description{
A dataset from: Durr, Robert H., Andrew D. Martin, and Christina Wolbrecht. 2000. "Ideological divergence and public support for the Supreme Court." American Journal of Policial Science 44(4): 768-776.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{lshift}
\alias{lshift}
\title{Take lag transformation of a series}
\usage{
lshift(x, l)
}
\arguments{
\item{x}{a series to be lagged}

\item{l}{the number of lags}
}
\value{
the lagged series
}
\description{
Take lag transformation of a series
}
\details{
\code{lshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
}
\examples{
x.var <- runif(50)
l.1.x.var <- lshift(x.var, 1)
l.2.x.var <- lshift(x.var, 2)
head(x.var)
head(l.1.x.var)
head(l.2.x.var)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dynardl.auto.correlated}
\alias{dynardl.auto.correlated}
\title{Run a variety of autocorrelation tests on the residuals from a \code{\link{dynardl}} model}
\usage{
dynardl.auto.correlated(x, bg.type = "Chisq", digits = 3,
  order = NULL, object.out = FALSE)
}
\arguments{
\item{x}{a \code{dynardl} model}

\item{bg.type}{a character string for the type of Breusch-Godfrey test to run. The default is \code{Chisq}: the Chisq test statistic. The other option is \code{F}: the F-test statistic}

\item{digits}{the number of digits to round to when showing output. The default is \code{3}}

\item{order}{the maximum order of serial autocorrelation to test when executing the Breusch-Godfrey test}

\item{object.out}{if \code{TRUE}, and \code{dynardl.auto.correlated} is assigned to an object, the AIC, BIC, and results will be stored for the user's convenience}
}
\value{
The results of autocorrelation tests
}
\description{
Run a variety of autocorrelation tests on the residuals from a \code{\link{dynardl}} model
}
\details{
This is a simple and convenient way to test whether the residuals from the \code{dynardl} model are white noise. As an aside, this is also why \code{dynardl} has a \code{simulate = FALSE} argument: users can ensure the model has white noise residuals before estimating a potentially time-intensive simulation. The output also reminds the user of the null hypotheses for the autocorrelation tests
}
\examples{
# Using the ineq data from dynamac
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = FALSE)
dynardl.auto.correlated(ardl.model)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dynardl.all.plots}
\alias{dynardl.all.plots}
\title{Combine all of the potential plots of a simulated response in a \code{\link{dynardl}} model}
\usage{
dynardl.all.plots(x, type = "area", bw = FALSE, last.period = NULL,
  start.period = 1, tol = (abs(x$model$ymean) * 0.01),
  abs.errors = "none", ylim = NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
\item{x}{a \code{dynardl} model with a simulation to be plotted. Since all plots includes absolute cumulative differences, \code{fullsims} must be \code{TRUE} in the \code{dynardl} simulation}

\item{type}{whether the plot should be an area plot (\code{area}) or a spike plot (\code{spike})}

\item{bw}{should the colors be in black and white (for publication)? The default is \code{FALSE}}

\item{last.period}{when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify a specific period in which to stop calculating absolute cumulative differences. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}}

\item{start.period}{which period of the simulation to begin the plot with. You can view the equilibriating behavior of the dependent variable, or you can skip forward in time (maybe to just before the shock). The default is \code{1} (the first period of the simulation)}

\item{tol}{when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify the minimum amount of movement required to qualify as a non-noise change over time periods (for calculating absolute cumulative differences). The default is 0.1 percent of the mean of the dependent variable. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}}

\item{abs.errors}{when calculating confidence for the absolute cumulative effect, should differences accumulate in each time time period (\code{cumulate}, which could be explosive if the error in the model is large), should differences be observed at each time (\code{within.period}, which will have smaller values in equilibrium than when changing), or should only the values be plotted (\code{none})}

\item{ylim}{a user-defined y-limit to be used instead of the default (for instance, for shared axes. Use caution, as it will be passed to all plots)}

\item{xlab}{a user-defined x-label to be used instead of the default (use caution, as it will be passed to all plots)}

\item{ylab}{a user-defined y-label to be used instead of the default (use caution, as it will be passed to all plots)}

\item{...}{other arguments to be passed to the call to plot. Use caution, as they will be passed to all plots}
}
\value{
a 2 x 3 grid of the plots of the simulated dynardl model effects plots
}
\description{
Combine all of the potential plots of a simulated response in a \code{\link{dynardl}} model
}
\details{
When running \code{dynardl}, \code{simulate} must be \code{TRUE} so that there is a simulation to plot. Also, \code{fullsims} must be \code{TRUE} as the plot will contain absolute cumulative differences. See \code{\link{dynardl.simulation.plot}} for arguments to the individual plotting types
}
\examples{
# Using the ineq data in dynamac
# Shocking Income Top 10
set.seed(1)
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = TRUE, range = 30,
       shockvar = "incshare10", fullsims = TRUE)

# Shows all of the potential responses
dynardl.all.plots(ardl.model)	
# Same plot, but with spikeplot
dynardl.all.plots(ardl.model, type = "spike")  
# Grayscale plots
dynardl.all.plots(ardl.model, bw = TRUE)	 
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dynardl.simulation.plot}
\alias{dynardl.simulation.plot}
\title{Create a plot of a simulated response in a \code{\link{dynardl}} model}
\usage{
dynardl.simulation.plot(x, type = "area", response = "levels",
  bw = FALSE, last.period = NULL, tol = (abs(x$model$ymean) * 0.01),
  start.period = 1, abs.errors = "none", ylim = NULL, ylab = NULL,
  xlab = NULL, ...)
}
\arguments{
\item{x}{a \code{dynardl} model with a simulation to be plotted}

\item{type}{whether the plot should be an area plot (\code{area}) or a spike plot (\code{spike})}

\item{response}{whether the plot of the response should be shown in levels of the dependent variable (\code{levels}), levels from the mean of the dependent variable (\code{levels.from.mean}), period-over-period changes in the dependent variable (\code{diffs}), the absolute value of the (decreasing) change in the dependent variable  in each time period due to the shock (\code{shock.effect.decay}), the sum of the period-over-period changes (\code{cumulative.diffs}), or the absolute value of the cumulative differences (where negative effects are treated as positive) (\code{cumulative.abs.diffs}). The default is \code{levels}}

\item{bw}{should the colors be in black and white (for publication)? The default is \code{FALSE}}

\item{last.period}{when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify a specific period in which to stop calculating absolute cumulative differences. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}}

\item{tol}{when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify the minimum amount of movement required to qualify as a non-noise change over time periods (for calculating absolute cumulative differences). The default is 0.1 percent of the mean of the dependent variable. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}}

\item{start.period}{which period of the simulation to begin the plot with. You can view the equilibriating behavior of the dependent variable, or you can skip forward in time (maybe to just before the shock). The default is \code{1} (the first period of the simulation)}

\item{abs.errors}{when calculating confidence for the absolute cumulative effect, should differences accumulate in each time time period (\code{cumulate}, which could be explosive if the error in the model is large), should differences be observed at each time (\code{within.period}, which will have smaller values in equilibrium than when changing), or should only the values be plotted (\code{none}). The default is \code{none}}

\item{ylim}{a user-defined y-limit to be used instead of the default (for instance, for shared axes)}

\item{ylab}{a user-defined y-label to be used instead of the default}

\item{xlab}{a user-defined x-label to be used instead of the default}

\item{...}{other arguments to be passed to the call to plot}
}
\value{
a plot of the simulated dynardl model
}
\description{
Create a plot of a simulated response in a \code{\link{dynardl}} model
}
\details{
When running \code{dynardl}, \code{simulate} must be \code{TRUE} so that there is a simulation to plot. For types \code{cumulative.diffs} and \code{cumulative.abs.diffs}, \code{fullsims} must be \code{TRUE} in the \code{dynardl} simulation
}
\examples{
# Using the ineq data in dynamac
# Shocking Income Top 10
set.seed(1)
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = TRUE, range = 30,
       shockvar = "incshare10", fullsims = TRUE)

# Shows absolute levels
dynardl.simulation.plot(ardl.model)	
# Shows changes from mean level
dynardl.simulation.plot(ardl.model, response = "levels.from.mean")  
# Same plot, but with spikeplot
dynardl.simulation.plot(ardl.model, type = "spike", response = "levels.from.mean")  
# Grayscale plots
dynardl.simulation.plot(ardl.model, bw = TRUE)	 
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{dynardl}
\alias{dynardl}
\title{Estimate and simulate ARDL model}
\usage{
dynardl(formula, data = list(), lags = list(), diffs = c(),
  lagdiffs = list(), levels = c(), ec = FALSE, trend = FALSE,
  constant = TRUE, modelout = FALSE, noLDV = FALSE,
  simulate = FALSE, shockvar = list(),
  shockval = sd(data[[shockvar]], na.rm = T), time = 10,
  qoi = "mean", forceset = NULL, range = 20, burnin = 20,
  sims = 1000, sig = 95, expectedval = FALSE, fullsims = FALSE)
}
\arguments{
\item{formula}{a symbolic description of the model to be estimated. ARDL
models are estimated using linear regression}

\item{data}{an optional data frame or list containing the the variables in the model}

\item{lags}{a list of variables and their corresponding lags to be estimated}

\item{diffs}{a vector of variables to be differenced. Only first differences are supported}

\item{lagdiffs}{a list of variables to be included in lagged differences}

\item{levels}{a vector of variables to be included in levels}

\item{ec}{estimate model in error-correction form, (i.e., \code{y} appears in first-differences). By default, \code{ec} is set to \code{FALSE}, meaning \code{y} will appear in levels.}

\item{trend}{include a linear time trend. The default is \code{FALSE}}

\item{constant}{include a constant. The default is \code{TRUE}}

\item{modelout}{print the regression estimates in the console}

\item{noLDV}{do not add a lagged dependent variable (LDV) to ARDL models when omitted in formula (special thanks to Hannes Datta). This is not recommended}

\item{simulate}{simulate the reponse. Otherwise, just the regression model will be estimated. If \code{simulate = FALSE}, options \code{shockvar}, \code{shockval},  \code{time}, \code{qoi}, \code{forceset}, \code{range}, \code{burnin}, \code{sims}, \code{sig}, \code{expectedval}, and \code{fullsims} are ignored. The default is \code{FALSE} so that users can build models without needing to simulate the results each time. When \code{simulate = TRUE}, users are highly encouraged to set a seed before simulation, as with any stochastic exercise}

\item{shockvar}{the variable to be shocked in the counterfactual simulation. There is no default}

\item{shockval}{the amount by which the \code{shockvar} should be shocked. The default is one standard deviation of the shocked variable}

\item{time}{the time period in the simulation for the variable to be shocked}

\item{qoi}{summarize the response of the dependent variable with the mean or the median. Although the default is \code{mean}, if there is underlying skew in the distribution, it might be better summarized by \code{median}}

\item{forceset}{by default, in the simulations, variables in levels will be set to their means; variables in differences will be set to 0. Alternatively, users can set any variable in the model to a different value using a list in \code{forceset}. These values can be any user-defined value, including means, medians, percentiles, or other values of interest}

\item{range}{the range of the simulation to be conducted}

\item{burnin}{the number of time periods to disregard before recording the values. These do not include the \code{range}; in other words, they take place before the \code{range} specified above. Users can increase the number of \code{burnin} periods, but probably should not decrease them. The default is 20}

\item{sims}{the number of simulations to use in creating the quantities of interest (the response of the dependent variable). The default is 1000}

\item{sig}{the significance level (1 - \code{p}) that the user wants for the simulations. The default level is 95\% significance (\code{sig = 95})}

\item{expectedval}{if this is \code{TRUE}, the simulation will record the expected values of across the \code{sims} by averaging errors. The default is \code{FALSE}, since expected values do not account for stochastic error present in the model itself}

\item{fullsims}{whether all of the raw simulations should be stored in the model object. These are required for some of the more advanced plotting functions, especially those that use the simulations to derive confidence intervals about the size of the period-over-period differences. The default is \code{FALSE}}
}
\value{
\code{dynardl} should always return an estimated model. It may or may not be simulated, according to the user. But the relevant regression output, model residuals (which can be tested for autocorrelation), and simulated response (if created) are stored in a list if the model is assigned to an object
}
\description{
Estimate autoregressive distributed lag models and simulate interesting values (if desired)
}
\details{
Estimate an auto-regressive distributed lag model. Moreover, enable a graphical interpretation of the results (through \code{\link{dynardl.simulation.plot}}) by simulating the response of the dependent variable to shocks in one of the regressors, and enable the Pesaran, Shin, and Smith (2001) test for cointegration for error-correction models (through \code{\link{pssbounds}})
}
\examples{
# Using the inequality data from dynamac
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       ec = TRUE, simulate = FALSE)
summary(ardl.model)

# Adding a lagged difference of the dependent variable
ardl.model.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = FALSE)
summary(ardl.model.2)

# Does not work: levels and diffs must appear as a vector
# ardl.model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#      lags = list("concern" = 1, "incshare10" = 1),
#       levels = list("urate" = 1),
#       diffs = list("incshare10" = 1, "urate" = 1), 
#       lagdiffs = list("concern" = 1),
#       ec = TRUE, simulate = FALSE)

ardl.model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
       lags = list("concern" = 1, "incshare10" = 1),
       levels = c("urate"),
       diffs = c("incshare10", "urate"), 
       lagdiffs = list("concern" = 1),
       ec = TRUE, simulate = FALSE)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{ardl}
\keyword{estimation}
\keyword{simulation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{ldshift}
\alias{ldshift}
\title{Take the lagged first difference of a series}
\usage{
ldshift(x, l)
}
\arguments{
\item{x}{a series to be differenced}

\item{l}{the number of lags}
}
\value{
the lagged differenced series
}
\description{
Take the lagged first difference of a series
}
\details{
\code{ldshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
}
\examples{
x.var <- runif(50)
ld.1.x.var <- ldshift(x.var, 1)
ld.2.x.var <- ldshift(x.var, 2)
head(x.var)
head(ld.1.x.var)
head(ld.2.x.var)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\name{pssbounds}
\alias{pssbounds}
\title{Perform Pesaran, Shin, and Smith (2001) cointegration test}
\usage{
pssbounds(data = list(), obs = NULL, fstat = NULL, tstat = NULL,
  case = NULL, k = NULL, restriction = FALSE, digits = 3,
  object.out = FALSE)
}
\arguments{
\item{data}{an optional \code{\link{dynardl}} model. This option is highly recommended. Users are welcome to supply their own case, k regressors, t-statistic, F-statistic, and observations, but it is easier to have the model determine these quantities. If a \code{\link{dynardl}} model is supplied, user-supplied arguments are ignored}

\item{obs}{number of observations}

\item{fstat}{F-statistic of the joint test that variables in first lags are equal to zero: the specific restriction tested 
is \code{l.y + l.1.x1 + l.1.x2 + ... + l.1.xk = 0}, except in cases II and IV (see \code{restriction} and \code{case})}

\item{tstat}{t-statistic of the lagged dependent variable}

\item{case}{The case of the test, as per Pesaran, Shin, and Smith (2001). Case I: no intercept or trend; case II: restricted intercept, no trend; case III: unrestricted intercept with no trend; case IV: unrestricted intercept and restricted trend; case V: unrestricted intercept and trend. Case III is most frequently specified}

\item{k}{number of regressors appearing in levels in the estimated model, not including the lagged dependent variable}

\item{restriction}{if you design to test case II or IV of pssbounds, where it is assumed that the constant (case 2) or trend (case 4) are restricted in the resulting F-test, indicate that restriction = \code{TRUE}. If restriction = \code{TRUE} and there is no trend in the regression (trend = \code{FALSE} in \code{\link{dynardl}}), the F-test will include the constant in addition to the lagged dependent variable and lagged regressors in order to test for cointegration under the assumption of a restricted constant (see Pesaran, Shin and Smith [2001], case II). If restriction = \code{TRUE} and there is a trend in the regression (trend = \code{TRUE} in \code{\link{dynardl}}), the F-test will include the trend term in addition to the lagged dependent variable and lagged regressors in order to test for cointegration under the assumption of a restricted trend (see Pesaran, Shin and Smith [2001], case IV). If you are estimating the regular unrestricted ECM (this is more common), restriction = \code{FALSE}. The default is \code{FALSE}}

\item{digits}{the number of digits to round to when showing output. The default is \code{3}}

\item{object.out}{if \code{TRUE}, and \code{pssbounds} is assigned to an object, the test quantities will be stored for the user's convenience}
}
\description{
Perform Pesaran, Shin, and Smith (2001) cointegration test
}
\details{
pssbounds performs post-estimation cointegration testing using the bounds testing procedure from Pesaran, Shin, and Smith (2001). Since test statistics vary based on the number of \code{k} regressors, length of the series, these are required, in addition to F- and t-statistics
}
\examples{
# Using the ineq data from dynamac
# We can get all the values by hand
ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
        lags = list("concern" = 1, "incshare10" = 1),
        diffs = c("incshare10", "urate"), 
        lagdiffs = list("concern" = 1),
        ec = TRUE, simulate = FALSE)
summary(ardl.model)
pssbounds(obs = 47, fstat = 7.01578, tstat = -3.223, case = 3, k = 1)

# Or just pass a dynardl model.
pssbounds(ardl.model)
}
\author{
Soren Jordan and Andrew Q. Philips
}
\keyword{cointegration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamac.R
\docType{data}
\name{france.data}
\alias{france.data}
\title{Data on French Energy Consumption and GDP}
\format{A data frame with 53 rows and 4 variables:
\describe{
  \item{country}{Country}
  \item{year}{Year}
  \item{lnGDP_cons2010USD}{ln(GDP), constant 2010 US dollars}
  \item{lnenergy}{ln(energy consumption), millions tons oil equivalent}
}}
\usage{
data(france.data)
}
\description{
Data on GDP are from World Bank World Development Indicators. Data
on energy consumption are from the PB Statistical Review of World
Energy (June 2018).
}
\keyword{datasets}
