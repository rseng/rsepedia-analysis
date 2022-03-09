# bistablehistory 1.0.0
## First CRAN Release
* Initial CRAN Release

# bistablehistory 1.1.0
## Improvements
* Custom prior values for history parameters, intercept terms, history effect, and fixed effects.
* Simplified Stan code.
* predict() computes values from history, reducing fit object size.
* predict() returns a vector of length that matches original table.

## Bug Fixes
* Change to difference in history values instead of the weighted mean.
* Use for scale instead of rate for prediction for Gamma family.
* Spelling in documentation.
* Additional tests.
# Cumulative History Analysis For Bistable Perception Time Series

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/299245172.svg)](https://zenodo.org/badge/latestdoi/299245172)
[![CRAN status](https://www.r-pkg.org/badges/version/bistablehistory)](https://cran.r-project.org/package=bistablehistory)
<!-- badges: end -->

A package to compute a cumulative history for time-series of perceptual
dominance in bistable displays.

Estimates cumulative history, an estimate of accumulating
adaptation/prediction error for the dominant percept, for time-series
for continuously viewed bistable perceptual rivalry displays. Computes
cumulative history via a homogeneous first order differential process.
I.e., it assumes exponential growth/decay of the history as a function
of time and perceptually dominant state. Supports Gamma, log normal, and
normal distribution families.

For details on rationale please refer to
([Pastukhov & Braun, 2011](https://doi.org/10.1167/11.10.12)).

## Installation

For current stable version use

```r
install.packages("bistablehistory")
```

The master branch is the development version. To install it please use

```r
library("devtools")
install_github("alexander-pastukhov/bistablehistory", dependencies = TRUE)
```

### Note

This package uses [Stan](https://mc-stan.org), a "state-of-the-art
platform for statistical modeling and high-performance statistical
computation". Therefore, it depends on the package
[rstantools](https://cran.r-project.org/package=rstantools), which in
turn depends on the [rstan](https://cran.r-project.org/package=rstan)
package, which uses the [V8 JavaScript library](https://v8.dev), through
the [V8 R package](https://cran.r-project.org/package=V8).

Therefore, you will need to install the V8 JavaScript library on your
system, and it is recommended that you also install the V8 R package
beforehand. For detailed instructions, please see
https://github.com/jeroen/v8.

You will also need the [R package
curl](https://cran.r-project.org/package=curl), which depends on
`libcurl-*` in various operating systems. Please see the documentation
at https://cran.r-project.org/package=curl.

## Usage

The main function is `fit_cumhist` that takes a data frame with
time-series as the first argument. Minimally, you need to specify `state`
--- string with the column name that encodes perceptually dominant state
--- and either `duration` (column name with duration of individual
dominance phases) or `onset` (column name with onset times of individual
dominance phases). Thus, for a simplest case of a single subject and
single run/block measurement with all defaults (gamma distribution,
fitted cumulative history time constant but fixed mixed state value and
history mixing proportion) the call would be

```r
library(bistablehistory)
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                         state = "State",
                         duration = "Duration")
```

or, equivalently

```r
library(bistablehistory)
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                         state = "State",
                         onset = "Time")
```

Now you can look at the fitted value for history time constant via

```r
history_tau(gamma_fit)
```

and main effect of history for both parameters of gamma distribution

```r
coef(gamma_fit)
```

For further details please see vignettes on package usage ([Usage
examples](https://CRAN.R-project.org/package=bistablehistory/vignettes/usage-examples.html)
and [Cumulative
history](https://CRAN.R-project.org/package=bistablehistory/vignettes/cumulative-history.html))
and on an example of writing Stan code directly ([Writing Stan
code](https://CRAN.R-project.org/package=bistablehistory/vignettes/writing-stan-code.html)).
---
title: "Writing Stan code"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Writing Stan code}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The package allows for only limited models as, e.g., neither random slopes, nor interaction effects are allowed. Imposing this restriction was a design decision, as it would require duplicating functionality of general purposes packages. Instead, the package itself provides some basic fitting that should be sufficient for most simple cases. However, below you will find example of how to incorporate cumulative history into a model written in Stan. This way, you can achieve maximal flexibility but still save time by reusing the code.

## Stan model

This is a complete Stan code for a model with log-normal distribution for multiple runs from a single experimental session of a single participant. The history time-constant `tau` is fitted, whereas constants are used for other cumulative history parameters.

```{stan output.var="example_model", eval=FALSE}
data{
    // --- Complete time-series ---
    int<lower=1> rowsN;     // Number of rows in the COMPLETE multi-timeseries table including mixed phase.
    real duration[rowsN];   // Duration of a dominance/transition phase
    int istate[rowsN];      // Index of a dominance istate, 1 and 2 code for two competing clear states, 3 - transition/mixed.
    int is_used[rowsN];     // Whether history value must used to predict duration or ignored
                            // (mixed phases, warm-up period, last, etc.)
    int run_start[rowsN];   // 1 marks a beginning of the new time-series (run/block/etc.)
    real session_tmean[rowsN]; // Mean dominance phase duration for both CLEAR percepts. Used to scale time-constant.
    
    // --- A shorter clear-states only time-series ---
    int clearN;                  // Number of rows in the clear-states only time-series
    real clear_duration[clearN]; // Duration for clear percepts only.
    
    // --- Cumulative history parameters
    real<lower=0, upper=1> history_starting_values[2]; // Starting values for cumulative history at the beginning of the run
    real<lower=0, upper=1> mixed_state;                // Mixed state signal strength
}
parameters {
    real<lower=0> tau; // history time-constant
    
    // linear model for mu
    real a;
    real bH;
    
    // variance
    real<lower=0> sigma;
}
transformed parameters{
    vector[clearN] mu; // vector of computed mu for each clear percept
  
    {
        // temporary variables
        real current_history[2]; // current computed history
        real tau_H;              // tau in the units of time
        real dH;                 // computed history difference
        int iC = 1;              // Index of clear percepts used for fitting

        // matrix with signal levels
        matrix[2, 3] level = [[1, 0, mixed_state], 
                              [0, 1, mixed_state]];

        for(iT in 1:rowsN){
            // new time-series, recompute absolute tau and reset history state
            if (run_start[iT]){
                // reset history
                current_history = history_starting_values;

                // Recompute tau in units of time. 
                // This is relevant only for multiple sessions / participants.
                // However, we left this code for generality.
                tau_H = session_tmean[iT] * tau;
            }

            // for valid percepts, we use history to compute mu
            if (is_used[iT] == 1){
                // history difference
                dH = current_history[3-istate[iT]] - current_history[istate[iT]];

                // linear model for mu
                mu[iC] = a + bH * dH;
                iC += 1;
            }

            // computing history for the NEXT episode
            // see vignette on cumulative history
            for(iState in 1:2){
                current_history[iState] = level[iState, istate[iT]] + 
                  (current_history[iState] - level[iState, istate[iT]]) * exp(-duration[iT] / tau_H);
            }
        }
    }
}
model{
  // sampling individual parameters
  tau ~ lognormal(log(1), 0.75);
  a ~ normal(log(3), 5);
  bH ~ normal(0, 1);
  sigma ~ exponential(1);
  
  // sampling data using computed mu and sampled sigma
  clear_duration ~ lognormal(exp(mu), sigma);
}
```


## Data preparation
The `data` section defines model inputs. Hopefully, the comments make understanding it fairly straightforward. However, it has several features that although are not needed for the limited single session / single session make it easier to generalized the code for more complicated cases. 

For example, not all dominance phases are used for fitting. Specifically, all mixed perception phases, first dominance phase for each percept (not enough time to form reliably history) and last dominance phase (curtailed by the end of the block) are excluded. Valid dominance phases are marked in `is_used` vector. Their total number is stored in `clearN` variable and the actual dominance durations in `clear_duration`. The latter is not strictly necessary but allows us to avoid a loop and vectorize the sampling statement `clear_duration ~ lognormal(exp(mu), sigma);`.

In addition, `session_tmean` is a vector rather than a scalar. This is not necessary for a single session example here but we opted to use as it will better generalize for more complicated cases.

bistability package provides a service function `preprocess_data()` that simplifies the process of preparing the data. However, you need to perform the last step, forming a list of inputs for Stan sampling, yourself.
```{r eval=FALSE}
# function that checks data for internal consistency and returns a preprocessed table
df <- bistablehistory::preprocess_data(br_single_subject, 
                                       state="State",
                                       duration="Duration",
                                       run="Block")

# data for Stan model
stan_data <- list(
  # complete time-series
  rowsN = nrow(df),
  duration = df$duration,
  istate = df$istate,
  is_used = df$is_used,
  run_start = df$run_start,
  session_tmean = df$session_tmean,
  
  # only valid clear percepts
  clearN = sum(df$is_used),
  clear_duration = df$duration[df$is_used == 1],
  
  # history parameters, all fixed to default values
  history_starting_values = c(0, 0),
  mixed_state = 0.5
)
```

## Using the model
You can use this model either with `rstan` or `cmdstanr` packages. Below is in an example using `cmdstanr`, assuming that model file is called `example.stan`. 
```{r eval=FALSE}
# compile the model
model <- cmdstanr::cmdstan_model("example.stan")

# sample model
fit <- model$sample(data=stan_data, chains=1)

# extract posterior samples for tau parameter
tau <- fit$draws(variables = "tau")
```

---
title: "Usage examples"
output: rmarkdown::html_vignette
bibliography: usage-examples.bib  
vignette: >
  %\VignetteIndexEntry{Usage examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>"
)
```

## Minimal example

The main function is `fit_cumhist()` that takes a data frame with time-series as a first argument. In addition, you need to specify the name of the column that codes the perceptual state (`state` argument) and a column that holds either dominance phase duration (`duration`) or its onset (`onset`). The code below fits data using Gamma distribution (default family) for a single run of a single participant. By default, the function fits cumulative history time constant but uses default fixed mixed state value (`mixed_state = 0.5`) and initial history values (`history_init = 0`).

```{r minimal example duration, warning = FALSE, message = FALSE}
library(bistablehistory)

data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock,
                         state="State",
                         duration="Duration",
                         refresh=0)
```

Alternatively, you specify _onset_ of individual dominance phases that will be used to compute their duration.
```{r eval=FALSE}
gamma_fit <- fit_cumhist(br_singleblock,
                        state="State",
                        onset="Time")
```

You can look at the fitted value for history time constant using `history_tau()`
```{r}
history_tau(gamma_fit)
```


and main effect of history for both parameters of gamma distribution
```{r}
historyef(gamma_fit)
```

The following model is fitted for the example above, see also [companion vignette](cumulative-history.html) for details on cumulative history computation.
$$Duration[i] \sim  Gamma(shape[i], rate[i]) \\
log(shape[i]) = \alpha^{shape} + \beta^{shape}_H \cdot \Delta h[i] \\
log(rate[i]) = \alpha^{rate} + \beta^{rate}_H \cdot \Delta h[i] \\
\Delta h[i] = \text{cumulative_history}(\tau, \text{history_init})\\
\alpha^{shape}, \alpha^{rate} \sim Normal(log(3), 5) \\
\beta^{shape}_H, \beta^{rate}_H \sim Normal(0, 1) \\
\tau \sim Normal(log(1), 0.15)$$

## Passing Stan control parameters
You can pass Stan control parameters via `control` argument, e.g.,
```{r eval=FALSE}
gamma_fit <- fit_cumhist(br_singleblock,
                        state="State",
                        duration="Duration",
                        control=list(max_treedepth = 15,
                                     adapt_delta = 0.99))
```

See Stan documentation for details [@Carpenter2017].

## Run
By default, `fit_cumhist()` function assumes that the time-series represent a single run, so that history states are initialized only once at the very beginning. You can use `run` argument to pass the name of a column that specifies individual runs. In this case, history is initialized at the beginning of every run to avoid spill-over effects.

```{r eval=FALSE}
gamma_fit <- fit_cumhist(br_single_subject,
                        state="State",
                        onset="Time",
                        run="Block")
```

## Experimental session
Experimental session specifies which time-series were measured together and is used to compute an average dominance phase duration that, in turn, is used when computing cumulative history: $\tau_H = \tau \cdot <D>$, where $\tau$ is normalized time constant and $<D>$ is the mean dominance phase duration. This can be used to account for changes in overall alternation rate between different sessions (days), as, for example, participants new to the stimuli tend to "speed up" over the course of days [@Suzuki2007]. If you _do not_ specify `session` parameter then a single mean dominance phase duration is computed for all runs of a single subject.

## Random effect
The `random_effect` argument allows you to specify a name of the column that codes for a random effect, e.g., participant identity, bistable display (if different displays were used for a single participant), etc. If specified, it is used to fit a hierarchical model with random slopes for the _history effect_ ($\beta_H$). Note that we if random _independent_ intercepts are used as prior research suggest large differences in overall alternation rate between participants [@Brascamp2019].

Here, is the R code that specifies participants as random effect
```{r eval=FALSE}
gamma_fit <-  fit_cumhist(kde_two_observers,
                          state="State",
                          duration="Duration",
                          random_effect="Observer",
                          run="Block")
```

And here is the corresponding model, specified for the shape parameter only as identical formulas are used for the rate parameter as well. Here, $R_i$ codes for a random effect level (participant identity) and a non-centered parametrization is used for the pooled random slopes.

$$Duration[i] \sim  Gamma(shape[i], rate[i]) \\
log(shape[i]) = \alpha[R_i] + \beta_H[R_i] \cdot \Delta h[i] \\
\Delta H[i] = \text{cumulative_history}(\tau, \text{history_init})\\
\alpha[R_i] \sim Normal(log(3), 5) \\
\beta_H[R_i] = \beta^{pop}_H + \beta^{z}_H[R_i] \cdot \sigma^{pop}_H\\
\beta^{pop}_H \sim Normal(0, 1) \\
\beta^{z}_H[R_i] \sim Normal(0, 1) \\
\sigma^{pop}_H \sim Exponential(1) \\
\tau \sim Normal(log(1), 0.15)$$

Identical approach is take for $\tau$, if `tau=' "1|random"'` was specified and same holds for `mixed_state=' "1|random"'` argument, see below.

## Fixed effects
`fit_cumhist()` functions allows you to specify multiple fixed effect terms as a vector of strings. The implementation is restricted to:

* Only continuous (metric) independent variables should be used.
* A single value is fitted for each main effect, irrespective of whether a random effect was specified.
* You cannot specify an interaction either between fixed effects or between a fixed effect and cumulative history variable.

Although this limits usability of the fixed effects, these restrictions allowed for both a simpler model specification and a simpler underlying code. If you do require more complex models, please refer to  [companion vignette](writing-stan-code.html) that provides an example on writing model using Stan directly.

You can specify custom priors (a mean and a standard deviation of a prior normal distribution) via `history_effect_prior` and `fixed_effects_priors` arguments. The former accepts a vector with mean and standard deviation, whereas the latter takes a named list in format \code{list("<fixed parameter name>"=c(<mean>, <std>))}.

Once fitted, you can use `fixef()` function to extract a posterior distribution or its summary for each effect.

## Cumulative history parameters
`fit_cumhist()` function takes three parameters for cumulative history computation (see also [companion vignette](cumulative-history.html)):

* `tau` : a _normalized_ time constant in units of mean dominance phase duration.
* `mixed_state` : value used for mixed/transition state phases, defaults to `0.5`.
* `history_init` : an initial value for cumulative history at the onset of each run. Defaults to `0`.

Note that although `history_init` accepts only fixed values either a single value used for both states or a vector of two. In contrast, both fixed and fitted values can be used for the other three parameters. Here are possible function argument values

* a single positive number for `tau` or single number within [0, 1] range for `mixed_state`. In this case, the value is used directly for the cumulative history computation, which is default option for  `mixed_state`.
* `NULL` : a single value is fitted and used for all participants and runs. This is a default for `tau`.
* `'random'` : an _independent_ tau is fitted for each random cluster (participant, displays, etc.). `random_effect` argument must be specified.
* `'1|random'` : values for individual random cluster are sampled from a fitted population distribution (_pooled_ values). `random_effect` argument must be specified.

You can specify custom priors for each cumulative history parameter via `history_priors` argument by specifying mean and standard deviation of a prior normal distribution. The `history_priors` argument must be a named list, \code{list("<parameter name>"=c(<mean>, <std>))}, e.g., `history_priors = list("tau"=c(1, 0.15))`.

Once fitted, you can use `history_tau()` and `history_mixed_state()`functions to obtain a posterior distribution or its summary for each parameter.


## Distribution family
`fit_cumhist` currently supports three distributions: `'gamma'`, `'lognormal'`, and `'normal'`.


### Gamma
$$Duration[i] \sim Gamma(shape[i], rate[i])$$
For Gamma distribution independent linear models with a log link function are fitted for both shape and rate parameter. Priors for intercepts for both parameters are $\alpha ~ Normal(log(3), 5)$.

### Log-normal
$$Duration[i] \sim LogNormal(\mu[i], \sigma)$$
The $\mu$ parameter is computed via a linear model with a log link function. Priors for the intercept are $\alpha ~ Normal(log(3), 5)$. Prior for $\sigma$ was $\sigma \sim Exponential(1)$.

### Normal
$$Duration[i] \sim Normal(\mu[i], \sigma)$$
The $\mu$ parameter is computed via a linear model. Priors for the intercept are $\alpha ~ Normal(3, 5)$. Prior for $\sigma$ was $\sigma \sim Exponential(1)$.

## Model comparison
Models fits can be compared via information criteria. Specifically, the log likelihood is stored in a `log_lik` parameter that can be directly using `loo::extract_log_lik()` function (see package [@@loo]) or used to compute either a leave-one-out cross-validation (via `loo()` convenience function) or WAIC (via `waic()`). These are information criteria that can be used for model comparison the same way as Akaike (AIC), Bayesian (BIC), or deviance (DIC) information criteria. The latter can also be computed from log likelihood, however, WAIC and LOOCV are both preferred for multi-level models, see [@Vehtari2017]. The model comparison itself can be performed via `loo::loo_compare()` function of the `loo` package.

## Computing and using cumulative history
If you are interested in the cumulative history itself, you can extract from the fitted object via `extract_history()` function
```{r eval=FALSE}
H <- extract_history(gam_fit)
```

Alternatively, you can skip fitting and compute history directly using predefined values via `compute_history()`.
```{r eval=FALSE}
df <- compute_history(br_singleblock,
                      state="State",
                      duration="Duration", 
                      tau=1,
                      mixed_state=0.5,
                      history_init=0)
```

## References

<div id="refs"></div>

---
title: "Cumulative History"
author: "Alexander (Sasha) Pastukhov"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: cumulative-history.bib  
vignette: >
  %\VignetteIndexEntry{Cumulative History}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  comment = "#>"
)

library(dplyr)
library(ggplot2)
```

Certain stimuli, such as a Necker cube depicted below, are compatible with several, typically two, comparably likely perceptual interpretations. In case of the Necker cube, one can perceive as "upwards" or "downwards" and during continuous viewing the perception alternates between these alternatives (these alternations are schematically depicted on the right).

```{r out.width="80%", fig.align='center'}
knitr::include_graphics("nc-stimulus.png")
```

A distribution of these so-called dominance phases typically has a right-skewed distribution that is frequently fitted using Gamma distribution.

```{r out.width="80%", fig.align='center'}
knitr::include_graphics("gamma-distribution.png")
```

However, the individual dominance phases show a subtle but consistent serial dependence, see, for example, [@VanEe2009]. This serial dependence is thought to reflect accumulation of slow adaptation [@Pastukhov2013] or prediction error [@Weilnhammer2017] for the dominant percept. This slow accumulation process can be described via a homogeneous first order process. A brief summary on the formula and examples of accumulation over time for different initial values, signal strength, and time constants are presented below. For further details, please refer to [@PastukhovBraun2011].

---

The cumulative history for a perceptual state is computed via a homogeneous first order process (for details on first order linear differential equations please refer to chapter 2 in @Wilson1999):
$$\tag{1}\frac{dh_i}{dt} = \frac{1}{\tau} (-h_i + S_i(t))$$

where $\tau$ is the time constant, $h_i$ is cumulative history and $S_i(t)$ is current signal level for for the i<sup>th</sup> perceptual state, so that
$$\tag{2}
S(t) = \begin{cases}
  1 & \text{if state $i$ is dominant}\\
  0 & \text{if state $i$ is suppressed}\\
  S_{mixed} & \text{if it is a mixed/transition phase and $0 ≥ S_{mixed} ≥1 $}
\end{cases}$$

where $S_{mixed}$ corresponds to the `mixed_state` parameter that can be either specified (the `fit_cumhist()` function uses a default of `0.5`) or fitted. The general solution for the equation (see Theorem 1 in @Wilson1999, chapter 2, page 15) is

$$\tag{3}h_i(t) = A e^{-t/\tau} + \frac{1}{\tau} \int_{0}^{t} e^{-(t'-t)/\tau} S(t') dt'$$

where $A$ is chosen to satisfy the initial condition. Assuming a constant signal $S$, we obtain
$$\tag{4}h_i(t) = A e^{-t/\tau} + S_i \cdot (1 - e^{-t/\tau})$$

For the cumulative history, we are interested in $h_i(t + \Delta t)$: a change following a dominance phase that starts at time $t$, ends at time $t + \Delta t$, and has a constant signal strength $S_i$. Assuming that a dominance phase starts at $t=0$ and substituting into the equation 4, $h_i(0) = A$. In other words, constant $A$ is equal to cumulative history state before the dominance phase onset and, therefore, the relative signal strength during the dominance phase is determined by the difference between the signal strength and the initial cumulative history value: $S_i - h_i(0)$. Thus
$$\tag{5} h_i(\Delta t) = h_i(0) + (S - h_i(0)) \cdot (1 - e^{-\Delta t/\tau})$$
$$\tag{6} h_i(\Delta t) = S + (h_i(0) - S)  \cdot e^{-\Delta t/\tau}$$


The figure below shows accumulation over time for three different initial values ($x(0)$), signal strength ($S$), and and time constants ($tau$). Note that the package allows to either specify and fit both the time constant (argument `tau` in `fit_cumhist()` function) and the initial history value at the block (`history_init` argument).

```{r  fig.width=8, fig.height=4, out.width="80%", fig.align='center'}
t <- seq(0, 20, length.out=100)
df <- 
  data.frame(h0 = c(0, 0.7, 0.9), S = c(1, 0, 0.5), tau = c(1, 4, 2)) %>%
  group_by(h0, S, tau) %>%
  summarise(h0 = h0[1],
            S = S[1],
            tau = tau[1],
            t = t,
            h = S + (h0 - S) * exp(-t / tau),
            Parameters = sprintf("h(0) = %.1f, S = %.1f, tau = %d", h0[1], S[1], tau[1]),
            .groups="keep")

ggplot(df, aes(x=t, y=h, color=Parameters)) + 
  geom_line() + 
  xlab("dt [s]") +
  ylab("h(t + Δt)") + 
  ylab("h(t + \u0394t)") + 
  ylim(0, 1) + 
  theme(legend.position="top")
```

As for a bistable case there are two history states (one for each perceptual state), we compute a history as a difference of cumulative histories
$$\tag{7}\Delta h(t, \tau) = h_{suppressed}(t, \tau) - h_{dominant}(t, \tau) $$
where $h_{dominant}$ and $h_{suppressed}$ are history states for the currently dominant and suppressed states, respectively. _E.g._, if a left eye dominates during following phase, $h_{dominant} = h_{left}$ and $h_{suppressed} = h_{right}$ and vice versa.

## References

<div id="refs"></div>

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{kde_two_observers}
\alias{kde_two_observers}
\title{Multirun data for two participants, kinetic-depth effect display}
\format{
A data frame with 1186 rows and 5 variables:
\describe{
\item{Observer}{Participant ID}
\item{Block}{Run / block index}
\item{State}{Factor variable for state with levels \code{-1} and \code{1} coding two clear perceptual states and \code{-2} the mixed / transition phase}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
kde_two_observers
}
\description{
Multirun data for two participants, kinetic-depth effect display
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_fixed_history_parameter.R
\name{check_fixed_history_parameter}
\alias{check_fixed_history_parameter}
\title{Evaluates values for a fixed history parameter}
\usage{
check_fixed_history_parameter(param_name, param_value, randomN, upperLimit)
}
\arguments{
\item{param_name}{Name of the parameter.}

\item{param_value}{A single value or \code{randomN} numeric values.}

\item{randomN}{Number of levels for the random variable.}

\item{upperLimit}{Upper limit for a valid \code{param_value}.}
}
\value{
A numeric vector \code{randomN} long.
}
\description{
Expects either a single value within a valid range or
\code{randomN} values.
}
\examples{
check_fixed_history_parameter("tau", 1, 10, Inf)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_R2.R
\name{bayes_R2}
\alias{bayes_R2}
\alias{bayes_R2.cumhist}
\title{Computes R-squared using Bayesian R-squared approach.}
\usage{
\method{bayes_R2}{cumhist}(object, summary = TRUE, probs = c(0.055, 0.945), ...)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}

\item{...}{Unused.}
}
\value{
vector of values or a data.frame with summary
}
\description{
For detail refer to:
Andrew Gelman, Ben Goodrich, Jonah Gabry, and Aki Vehtari (2018).
R-squared for Bayesian regression models. The American Statistician
\doi{10.1080/00031305.2018.1549100} and
\url{https://avehtari.github.io/bayes_R2/bayes_R2.html}
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
bayes_R2(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bistablehistory-package.R
\docType{package}
\name{bistablehistory-package}
\alias{bistablehistory-package}
\alias{bistablehistory}
\title{Cumulative History Analysis for Bistable Perception Time Series}
\description{
Estimates cumulative history for time-series for continuously
viewed bistable perceptual rivalry displays. Computes cumulative history
via a homogeneous first order differential process. I.e., it assumes
exponential growth/decay of the history as a function time and perceptually
dominant state, Pastukhov & Braun (2011) \doi{10.1167/11.10.12}.
Supports Gamma, log normal, and normal distribution families.
Provides a method to compute history directly and example of using the
computation on a custom Stan code.
}
\references{
Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
}
\seealso{
\code{vignette("cumulative-history", package = "bistablehistory")}
\code{vignette("usage-examples", package = "bistablehistory")}
\code{vignette("writing-stan-code", package = "bistablehistory")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/historyef.R
\name{historyef}
\alias{historyef}
\title{Extract the history-effects estimates}
\usage{
historyef(object, summary = TRUE, probs = c(0.055, 0.945))
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}
}
\value{
data.frame with values or summary
}
\description{
Extracts models population-level coefficients history-specific terms
for every modeled distribution parameter.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
historyef(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_history_parameter.R
\name{extract_history_parameter}
\alias{extract_history_parameter}
\title{Extracts a history parameter as a matrix}
\usage{
extract_history_parameter(
  object,
  param_name,
  samplesN = NULL,
  link_function = NULL
)
}
\arguments{
\item{object}{A \link[=cumhist-class]{cumhist} object}

\item{param_name}{String, a name of the parameter}

\item{samplesN}{Number of samples, if NULL is computed from rstan (but it is cheaper to do this once).}

\item{link_function}{A link function to use (exp or inv.logit) or \code{NULL} for identity.}
}
\value{
Matrix with \code{samplesN} rows and randomN
(found in \code{object$data$randomN}) columns
}
\description{
Extracts a history parameter as a matrix with
\code{samplesN} rows and randomN (found in \code{object$data$randomN})
columns.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
extract_history_parameter(br_fit, "tau", link_function = exp)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_term.R
\name{extract_term_to_matrix}
\alias{extract_term_to_matrix}
\title{Extracts a term with one column per fixed or random-level into a matrix}
\usage{
extract_term_to_matrix(object, term)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{term}{String, term name}
}
\value{
Matrix
}
\description{
Extracts a 3D array for a term with  sample, linear-model,
random/fixed-effect order and returns a matrix with samples as rows
and columns in order 1) all random/fixed effects for lm1, 2) all
random/fixed effects for lm2, etc.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
a <- extract_term_to_matrix(br_fit, "a")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{br_singleblock}
\alias{br_singleblock}
\title{Single run for binocular rivalry stimulus}
\format{
A data frame with 76 rows and 6 variables:
\describe{
\item{Observer}{Participant ID, all rows contain \emph{"ap"}}
\item{Group}{Display, all rows contain \emph{"BR"}}
\item{Block}{Run / block index, all rows contain \emph{1}}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Index of a perceptually dominant state, \emph{1}, \emph{2} - clear perceptual state, \emph{3} mixed / transition phase}
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
br_singleblock
}
\description{
A single subject / single run dataset for binocular rivalry.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loo.R
\name{loo.cumhist}
\alias{loo.cumhist}
\title{Computes an efficient approximate leave-one-out
cross-validation via loo library. It can be used
for a model comparison via loo::loo_compare() function.}
\usage{
\method{loo}{cumhist}(x, ...)
}
\arguments{
\item{x}{A \link[=cumhist-class]{cumhist} object}

\item{...}{unused}
}
\value{
A named list, see \code{\link[loo:loo]{loo::loo()}} for details.
}
\description{
Computes an efficient approximate leave-one-out
cross-validation via loo library. It can be used
for a model comparison via loo::loo_compare() function.
}
\examples{
data(br_singleblock)
\donttest{
gamma_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
loo_gamma <- loo(gamma_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{br}
\alias{br}
\title{Binocular rivalry data}
\format{
A data frame with 3769 rows and 6 variables:
\describe{
\item{Observer}{Participant ID.}
\item{Display}{Display, all rows contain \code{"BR"}}
\item{Block}{Run / block index.}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
br
}
\description{
Dataset on binocular rivalry for eight participants.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumhist-class.R
\docType{class}
\name{cumhist-class}
\alias{cumhist-class}
\alias{cumhist}
\title{Class \code{cumhist}.}
\description{
Cumulative history model fitted to time-series data.
}
\details{
See \code{methods(class = "cumhist")} for an overview of available methods.
}
\section{Slots}{

\describe{
\item{\code{family}}{A \code{string} with distribution family.}

\item{\code{data}}{A \code{list} with preprocessed data.}

\item{\code{stanfit}}{a \code{\link[rstan:stanfit-class]{stanfit}} object.}
}}

\seealso{
\code{\link{fit_cumhist}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/history_parameter.R
\name{history_parameter}
\alias{history_parameter}
\title{Extract values of used or fitted history parameter}
\usage{
history_parameter(
  object,
  param,
  summary = TRUE,
  probs = c(0.055, 0.945),
  includePopulationLevel = TRUE
)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{param}{Parameter name: \code{"tau"} or \code{"mixed_state"}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}

\item{includePopulationLevel}{Logical, for pooled random effect only. Whether to include
population mean as a separate \code{"_population"} level, default to \code{TRUE}.}
}
\value{
A vector, if summary was not requested. Or a tibble with a summary or if a fixed value was used.
}
\description{
Extract values of used or fitted history parameter
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
history_parameter(br_fit, "tau")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_history_parameter.R
\name{evaluate_history_init}
\alias{evaluate_history_init}
\title{Evaluates validity of initial history values.}
\usage{
evaluate_history_init(history_init)
}
\arguments{
\item{history_init}{Either a single value or a pair of values within 0..1 range.}
}
\value{
A vector of two values
}
\description{
Checks number and range of values. If a scalar is passed, uses same value
for both states.
}
\examples{
evaluate_history_init(0.5)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{kde}
\alias{kde}
\title{Kinetic-depth effect data}
\format{
A data frame with 38698 rows and 6 variables:
\describe{
\item{Observer}{Participant ID.}
\item{Display}{Display, all rows contain \code{"KD"}}
\item{Block}{Run / block index.}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
kde
}
\description{
Dataset on kinetic-depth effect for eleven participants.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{br_single_subject}
\alias{br_single_subject}
\title{Single experimental session for binocular rivalry stimulus}
\format{
A data frame with 76 rows and 6 variables:
\describe{
\item{Observer}{Participant ID, all rows contain \emph{"ap"}}
\item{Display}{Display, all rows contain \emph{"BR"}}
\item{Block}{Run / block index}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Index of a perceptually dominant state, \emph{1}, \emph{2} - clear perceptual state, \emph{3} mixed / transition phase}
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
br_single_subject
}
\description{
A single subject / multiple runs dataset for binocular rivalry.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/history_parameter.R
\name{history_mixed_state}
\alias{history_mixed_state}
\title{Extract values of used or fitted history parameter mixed_state}
\usage{
history_mixed_state(
  object,
  summary = TRUE,
  probs = c(0.055, 0.945),
  includePopulationLevel = TRUE
)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}

\item{includePopulationLevel}{Logical, for pooled random effect only. Whether to include
population mean as a separate \code{"_population"} level, default to \code{TRUE}.}
}
\value{
A single value, if fixed value was used. A vector or a tibble, depending on the
option used (single intercept, independent or random intercepts), and whether summary was
requested.
}
\description{
A short-cut for \code{history_parameter(object, "mixed_state", ...)}.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
history_tau(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_data.R
\name{preprocess_data}
\alias{preprocess_data}
\title{Preprocesses time-series data for fitting}
\usage{
preprocess_data(
  data,
  state,
  duration = NULL,
  onset = NULL,
  random_effect = NULL,
  session = NULL,
  run = NULL
)
}
\arguments{
\item{data}{A table with one or many time-series.}

\item{state}{String, the name of the column that specifies
perceptual state. The column type should be a factor with
two or three levels (the third level is assumed to correspond to a
transition/mixed phase) or should be convertible to a two level
factor (as it would be impossible to infer the identity of transition/
mixed phase).}

\item{duration}{String, name of the column with duration of individual
perceptual dominance phases. Optional, you can specify \code{onset}
instead.}

\item{onset}{String, name of the column with onsets of the perceptual
dominance states. Optional, used to compute duration of the dominance
phases, if these are not provided explicitly via \code{duration}
parameter.}

\item{random_effect}{String, name of the column that identifies random effect,
e.g. individual participants, stimuli for a single participant, etc.
If omitted, no random effect is assumed. If specified and
there is more than one level (participant, stimulus, etc.), it is used
in a hierarchical model.}

\item{session}{String, name of the column that identifies unique
experimental session for which a mean dominance phase duration will
be computed (see \code{norm_tau} parameter). Code assumes that session
IDs are different within a participant but can be the same between them.
If omitted, a single mean dominance duration based on the entire time series
is used.}

\item{run}{String, name of the column that identifies unique runs/blocks.
If omitted, the data is assumed to belong to a single time series. Code
assumes that run IDs are different within an experimental session but
can be the same between the session. E.g. session A, runs 1, 2, 3.. and
session B, runs 1, 2, 3 but not session A, runs 1, 2, 1.}
}
\value{
A tibble with columns
\itemize{
\item \code{state}
\item \code{duration}
\item \code{random}
\item \code{irandom} - integer, index of \code{random} values,
\item \code{session}
\item \code{run}
\item \code{session_tmean} - numeric, mean duration of clear percepts for every combination of \code{random} and \code{session}.
\item \code{is_used} - integer, whether computed history value needs to be used for linear model fitting.
\item \code{run_start} - integer, 1 for the first row of the run time-series.
}
}
\description{
Performs sanity checks (e.g., whether \code{data} can be used as a data.frame),
computes duration of dominance phases (if necessary), assumes a single entry for
any missing \code{session}, \code{run}, \code{random_effect}.
}
\examples{
df <- preprocess_data(br_singleblock, state="State", duration="Duration")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_normal_prior.R
\name{check_normal_prior}
\alias{check_normal_prior}
\title{Checks for validity of values for use as normal distribution parameters.}
\usage{
check_normal_prior(values, parameter)
}
\arguments{
\item{values}{Parameters for normal distribution.}

\item{parameter}{Name of the parameter for which the prior is defined.}
}
\value{
Logical TRUE, if none of the tests fail
}
\description{
Should a pair of numeric values, second value should be non-zero.
Stops execution with an error.
}
\examples{
check_normal_prior(c(0, 1), "tau")
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fast_history_compute}
\alias{fast_history_compute}
\title{Computes cumulative history}
\usage{
fast_history_compute(df, normalized_tau, mixed_state, history_init)
}
\arguments{
\item{df}{DataFrame with \code{"state"} (integer, 1 and 2 clear state, 3 - mixed state), \code{"duration"} (double),
\code{"irandom"} (integer, 1-based index of a random cluster), \code{"run_start"} (integer, 1 for the first entry of
the run, 0 otherwise), \code{"session_tmean"} (double)}

\item{normalized_tau}{DoubleVector A normalized tau value for each random cluster / individual. Thus, its length must be
equal to the number of unique indexes in \code{df["irandom"]}.}

\item{mixed_state}{DoubleVector A values used for the mixed state for each random cluster / individual.
Thus, its length must be equal to the number of unique indexes in \code{df["irandom"]}.}

\item{history_init}{DoubleVector, size 2. Initial values of history for a run.}
}
\value{
NumericMatrix, size \code{df.nrows()} × 2. Computed history values for each state.
}
\description{
Computes cumulative history based on common \code{history} values and
\code{normalized_tau} and \code{mixed_state} that are defined for each
random cluster / individual.
}
\examples{
df <- preprocess_data(br_singleblock, state="State", duration="Duration")
fast_history_compute(df, 1, 0.5, c(0, 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_history.R
\name{compute_history}
\alias{compute_history}
\title{Computes cumulative history for the time-series}
\usage{
compute_history(
  data,
  state,
  duration = NULL,
  onset = NULL,
  random_effect = NULL,
  session = NULL,
  run = NULL,
  tau = 1,
  mixed_state = 0.5,
  history_init = 0
)
}
\arguments{
\item{data}{A table with time-series.}

\item{state}{String, the name of the column that specifies
perceptual state. The column type should be a factor with
two or three levels (the third level is assumed to correspond to a
transition/mixed phase) or should be convertible to a two level
factor (as it would be impossible to infer the identity of transition/
mixed phase).}

\item{duration}{String, name of the column with duration of individual
perceptual dominance phases. Optional, you can specify \code{onset}
instead.}

\item{onset}{String, name of the column with onsets of the perceptual
dominance states. Optional, used to compute duration of the dominance
phases, if these are not provided explicitly via \code{duration}
parameter.}

\item{random_effect}{String, name of the column that identifies random effect,
e.g. individual participants, stimuli for a single participant, etc.
If omitted, no random effect is assumed. If specified and
there is more than one level (participant, stimulus, etc.), it is used
in a hierarchical model.}

\item{session}{String, name of the column that identifies unique
experimental session for which a mean dominance phase duration will
be computed (see \code{norm_tau} parameter). Code assumes that session
IDs are different within a participant but can be the same between them.
If omitted, a single mean dominance duration based on the entire time series
is used.}

\item{run}{String, name of the column that identifies unique runs/blocks.
If omitted, the data is assumed to belong to a single time series. Code
assumes that run IDs are different within an experimental session but
can be the same between the session. E.g. session A, runs 1, 2, 3.. and
session B, runs 1, 2, 3 but not session A, runs 1, 2, 1.}

\item{tau}{Time constant of exponential growth/decay
normalized to the mean duration of clear percepts within each \code{session}.
Can be 1) a single positive number (>0) that is used for all participants and runs,
2) \code{NULL} (default) -  a \emph{single} value will be fitted for all participants and runs,
3) \code{"random"} - an independent tau is fitted for each random cluster,
4) \code{"1|random"}- a tau for a random cluster
is sampled from a population distribution, i.e., pooled parameter values via
a multilevel model.}

\item{mixed_state}{Specifies an activation level during
transition/mixed phases (state #3, see \code{state}). Either a single
number (range 0..1) that will be used as a fixed level or a vector
of two numbers \code{c(mu, kappa)} that specifies, correspondingly, mean
(range 0..1) and precision (>0) of beta proportion distribution, it
should be sampled from. Defaults to a fixed value of \code{0.5}.}

\item{history_init}{Initial value for cumulative history computation. Either
a numeric scalar in 0..1 range or a vector of two numbers in 0..1 range.
In the latter case, two histories will start at different levels.}
}
\value{
A matrix \code{nrow(data)} × 2 with computed history values
}
\description{
Computes cumulative history for each state in the time-series.
}
\examples{
df <- compute_history(br_singleblock, state = "State",
                      duration = "Duration", tau = 1,
                      mixed_state = 0.5, history_init = 0)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{predict_samples}
\alias{predict_samples}
\title{Computes prediction for a each sample.}
\usage{
predict_samples(
  family,
  fixedN,
  randomN,
  lmN,
  istate,
  duration,
  is_used,
  run_start,
  session_tmean,
  irandom,
  fixed,
  tau_ind,
  mixed_state_ind,
  history_init,
  a,
  bH,
  bF,
  sigma
)
}
\arguments{
\item{family}{int, distribution family: gamma (1), lognormal(2), or
normal (3).}

\item{fixedN}{int, number of fixed parameters (>= 0).}

\item{randomN}{int, number of random factors (>= 1).}

\item{lmN}{int, number of linear models (>= 1).}

\item{istate}{IntegerVector, zero-based perceptual state 0 or 1,
2 is mixed state.}

\item{duration}{DoubleVector, duration of a dominance phase.}

\item{is_used}{IntegerVector, whether dominance phase is used for
prediction (1) or not (0).}

\item{run_start}{IntegerVector, 1 whenever a new run starts.}

\item{session_tmean}{DoubleVector, average dominance phase duration.}

\item{irandom}{IntegerVector, zero-based index of a random effect.}

\item{fixed}{NumericMatrix, matrix with fixed effect values.}

\item{tau_ind}{NumericMatrix, matrix with samples of tau for each
random level.}

\item{mixed_state_ind}{NumericMatrix, matrix with samples of
mixed_state for each random level.}

\item{history_init}{DoubleVector, Initial values of history for a run}

\item{a}{NumericMatrix, matrix with samples of
a (intercept) for each random level.}

\item{bH}{NumericMatrix, matrix with sample of
bH for each linear model and random level.}

\item{bF}{NumericMatrix, matrix with sample of
bF for each linear model and fixed factor.}

\item{sigma}{DoubleVector, samples of sigma.}
}
\value{
NumericMatrix with predicted durations for each sample.
}
\description{
Computing prediction for each sample,
recomputing cumulative history and uses
fitted parameter values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.cumhist}
\alias{coef.cumhist}
\title{Extract Model Coefficients}
\usage{
\method{coef}{cumhist}(object, summary = TRUE, probs = c(0.055, 0.945), ...)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}

\item{...}{Unused.}
}
\value{
data.frame with values or summary
}
\description{
Extracts models population-level coefficients history-specific terms and
fixed-effect terms for every modeled distribution parameter.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock,
                      state = "State",
                      duration = "Duration",
                      fixed_effects = "Time")
coef(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{nc}
\alias{nc}
\title{Necker cube data}
\format{
A data frame with 3464 rows and 6 variables:
\describe{
\item{Observer}{Participant ID.}
\item{Display}{Display, all rows contain \code{"NC"}}
\item{Block}{Run / block index.}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\source{
\doi{10.1167/11.10.12}
}
\usage{
nc
}
\description{
Dataset on Necker cube for five participants.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.cumhist}
\alias{summary.cumhist}
\title{Summary for a cumhist object}
\usage{
\method{summary}{cumhist}(object, ...)
}
\arguments{
\item{object}{A \link[=cumhist-class]{cumhist} object}

\item{...}{Unused}
}
\value{
Nothing, console output only.
}
\description{
Summary for a cumhist object
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
summary(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_history.R
\name{extract_history}
\alias{extract_history}
\title{Computes history for a fitted model}
\usage{
extract_history(object)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}
}
\value{
A matrix of cumulative history values for each state
}
\description{
Computes history for a fitted model, uses only mean values
for each history parameter. Uses values for each random cluster,
if \code{"random"} or \code{"1|random"} parametrisation was used.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
extract_history(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/history_parameter.R
\name{history_tau}
\alias{history_tau}
\title{Extract values of used or fitted history parameter tau}
\usage{
history_tau(
  object,
  summary = TRUE,
  probs = c(0.055, 0.945),
  includePopulationLevel = TRUE
)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}

\item{includePopulationLevel}{Logical, for pooled random effect only. Whether to include
population mean as a separate \code{"_population"} level, default to \code{TRUE}.}
}
\value{
A single value, if fixed value was used. A vector or a tibble, depending on the
option used (single intercept, independent or random intercepts), and whether summary was
requested.
}
\description{
A short-cut for \code{history_parameter(object, "tau", ...)}.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
history_tau(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{br_contrast}
\alias{br_contrast}
\title{Binocular rivalry, variable contrast}
\format{
A data frame with 4616 rows and 6 variables:
\describe{
\item{Observer}{Participant ID.}
\item{Block}{Run / block index.}
\item{Contrast}{Contrast on scale from 0 to 1.}
\item{Time}{Time relative to the run onset in \emph{seconds}}
\item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
\item{Duration}{Duration of a dominance phase in \emph{seconds}. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
}
}
\usage{
br_contrast
}
\description{
Dataset on binocular rivalry with variable but equal
contrast for six participants.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waic.R
\name{waic.cumhist}
\alias{waic.cumhist}
\title{Computes widely applicable information criterion
(WAIC).}
\usage{
\method{waic}{cumhist}(x, ...)
}
\arguments{
\item{x}{A \link[=cumhist-class]{cumhist} object.}

\item{...}{Additional arguments (unused)}
}
\value{
A named list, see \code{\link[loo:waic]{loo::waic()}} for details.
}
\description{
Computes widely applicable information criterion
via \link[=loo-package]{loo} library. It can be used for a model comparison via
\link[loo:loo_compare]{loo::loo_compare()} function.
}
\examples{
 \donttest{
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
waic_gamma <- waic(gamma_fit)
normal_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration", family="normal")
waic_normal <- waic(normal_fit)
loo::loo_compare(waic_gamma, waic_normal)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_cumhist.R
\name{fit_cumhist}
\alias{fit_cumhist}
\title{Fits cumulative history for bistable perceptual rivalry displays.}
\usage{
fit_cumhist(
  data,
  state,
  duration = NULL,
  onset = NULL,
  random_effect = NULL,
  session = NULL,
  run = NULL,
  fixed_effects = NULL,
  tau = NULL,
  mixed_state = 0.5,
  history_init = 0,
  family = "gamma",
  history_priors = NULL,
  intercept_priors = NULL,
  history_effect_prior = NULL,
  fixed_effects_priors = NULL,
  chains = 1,
  cores = NULL,
  ...
)
}
\arguments{
\item{data}{A table with time-series.}

\item{state}{String, the name of the column that specifies
perceptual state. The column type should be a factor with
two or three levels (the third level is assumed to correspond to a
transition/mixed phase) or should be convertible to a two level
factor (as it would be impossible to infer the identity of transition/
mixed phase).}

\item{duration}{String, name of the column with duration of individual
perceptual dominance phases. Optional, you can specify \code{onset}
instead.}

\item{onset}{String, name of the column with onsets of the perceptual
dominance states. Optional, used to compute duration of the dominance
phases, if these are not provided explicitly via \code{duration}
parameter.}

\item{random_effect}{String, name of the column that identifies random
effect, e.g. individual participants, stimuli for a single participant,
etc. If omitted, no random effect is assumed. If specified and
there is more than one level (participant, stimulus, etc.), it is used
in a hierarchical model.}

\item{session}{String, name of the column that identifies unique
experimental session for which a mean dominance phase duration will
be computed (see \code{norm_tau} parameter). Code assumes that session
IDs are different within a participant but can be the same between them.
If omitted, a single mean dominance duration based on the entire time series
is used.}

\item{run}{String, name of the column that identifies unique runs/blocks.
If omitted, the data is assumed to belong to a single time series. Code
assumes that run IDs are different within an experimental session but
can be the same between the session. E.g. session A, runs 1, 2, 3.. and
session B, runs 1, 2, 3 but not session A, runs 1, 2, 1.}

\item{fixed_effects}{String or vector of strings. Name of column(s)
with values to be used for fitting an additional fixed effect(s). E.g.,
contrast in binocular rivalry, rotation speed for kinetic-depth effect,
etc.}

\item{tau}{Time constant of exponential growth/decay
normalized to the mean duration of clear percepts within each \code{session}.
Can be 1) a single positive number (>0) that is used for all participants and runs,
2) \code{NULL} (default) -  a \emph{single} value will be fitted for all participants and runs,
3) \code{"random"} - an independent tau is fitted for each random cluster,
4) \code{"1|random"}- a tau for a random cluster
is sampled from a population distribution, i.e., pooled parameter values via
a multilevel model.}

\item{mixed_state}{Specifies an activation level during
transition/mixed phases (state #3, see \code{state}). Either a single
number (range 0..1) that will be used as a fixed level or a vector
of two numbers \code{c(mu, kappa)} that specifies, correspondingly, mean
(range 0..1) and precision (>0) of beta proportion distribution, it
should be sampled from. Defaults to a fixed value of \code{0.5}.}

\item{history_init}{Initial value for cumulative history computation. Either
a numeric scalar in 0..1 range or a vector of two numbers in 0..1 range.
In the latter case, two histories will start at different levels.}

\item{family}{String, distribution used to fit duration of perceptual dominance
phases. Options include \code{"gamma"} (default), \code{"lognormal"}, and \code{"normal"}.}

\item{history_priors}{Named list of optional priors for population-level cumulative history
parameters. Must follow the format \code{list("tau"=c(1, 0.15))} with values coding mean
and standard deviation of the normal distribution.}

\item{intercept_priors}{A vector of optional priors for population-level intercept
parameter. Should be \code{c(<shape-mean>, <shape-sd>, <scale-mean>, <scale-sd>)}
format for Gamma family, \code{c(<mean>, <sd>)} for normal and lognormal families.
The values code mean and standard deviation of the normal distribution.}

\item{history_effect_prior}{A vector of options priors for population-level slope
of history effect. The values code mean and standard deviation of the normal distribution.
Defaults to mu=0, sigma=1.}

\item{fixed_effects_priors}{A named list of optional priors for fixed effects. Must
follow the format \code{list("<name-of-variable>"=c(<mu>, <sigma>))}, where \code{<mu>} and
\code{<sigma>} are mean and standard deviation of a normal distribution. Defaults to mu=0,
sigma=1.}

\item{chains}{Number of chains for sampling.}

\item{cores}{Number of CPU cores to use for sampling. If omitted, All cores are used.}

\item{...}{Additional arguments passed to \link[rstan:stanmodel-method-sampling]{rstan::sampling()} function.}
}
\value{
An object of class \link[=cumhist-class]{cumhist}
}
\description{
Fits a generalized linear model using cumulative history and
specified fixed effects.
}
\examples{
\donttest{
data(br_singleblock)
gamma_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixef.R
\name{fixef}
\alias{fixef}
\title{Extract the fixed-effects estimates}
\usage{
fixef(object, summary = TRUE, probs = c(0.055, 0.945))
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to 89\% credible interval.}
}
\value{
\code{tibble} with values or summary, \code{NULL} if not fixed effects were used.
}
\description{
Extracts models fixed-effect terms for every modeled distribution parameter.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock,
                      state = "State",
                      duration = "Duration",
                      fixed_effects = "Time")
fixef(br_fit)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.cumhist}
\alias{print.cumhist}
\title{Prints out cumhist object}
\usage{
\method{print}{cumhist}(x, ...)
}
\arguments{
\item{x}{A \link[=cumhist-class]{cumhist} object}

\item{...}{Unused}
}
\value{
Nothing, console output only.
}
\description{
Prints out cumhist object
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration", fixed_effects="Time")
br_fit
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.cumhist}
\alias{predict.cumhist}
\title{Computes predicted dominance phase durations using posterior predictive distribution.}
\usage{
\method{predict}{cumhist}(object, summary = TRUE, probs = NULL, full_length = TRUE, ...)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{summary}{Whether summary statistics should be returned instead of
raw sample values. Defaults to \code{TRUE}}

\item{probs}{The percentiles used to compute summary, defaults to NULL (no CI).}

\item{full_length}{Only for \code{summary = TRUE}, whether the summary table should
include rows with no predictions. I.e., rows with mixed phases, first/last dominance
phase in the run, etc. See \code{\link[=preprocess_data]{preprocess_data()}}. Defaults to \code{TRUE}.}

\item{...}{Unused}
}
\value{
If \code{summary=FALSE}, a numeric matrix iterationsN x clearN.
If \code{summary=TRUE} but \code{probs=NULL} a vector of mean predicted durations.
If \code{summary=TRUE} and \code{probs} is not \code{NULL}, a data.frame
with a column \emph{"Predicted"} (mean) and a column for each specified quantile.
}
\description{
Computes predicted dominance phase durations using fitted model. Returns predicted
values only for the dominance phases that were marked for use. I.e., excluding first
and last dominance phases, mixed phases, etc. See \code{\link[=preprocess_data]{preprocess_data()}}.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
predict(br_fit)

# full posterior prediction samples
predictions_samples <- predict(br_fit, summary=FALSE)
}
}
\seealso{
\code{\link{fit_cumhist}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_history_parameter.R
\name{evaluate_history_option}
\alias{evaluate_history_option}
\title{Evaluates whether and how to fit a cumulative history parameter.}
\usage{
evaluate_history_option(param_name, param_value, randomN, upperLimit)
}
\arguments{
\item{param_name}{Name of the parameter.}

\item{param_value}{Value from the \code{\link{fit_cumhist}} function call.}

\item{randomN}{Number of levels for the random variable.}

\item{upperLimit}{Upper limit for a valid \code{param_value}.}
}
\value{
a list with \code{<param_name>_option} and \code{fixed_<param_name>}.
}
\description{
Evaluation is based on the \code{param_value}.
\enumerate{
\item A single positive number (>0) that is used for all participants and runs.
\item \code{NULL} (default) -  a \emph{single} value will be fitted for all participants
and runs, also applied if \code{randomN == 1}.
\item \code{"random"} - an independent value is fitted for each random cluster.
\item \code{"1|random"}- a value for a random cluster is sampled from a population
distribution, i.e., pooled parameter values via a multilevel model.
}
}
\examples{
evaluate_history_option("tau", 1, 1, Inf)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_term.R
\name{extract_replicate_term_to_matrix}
\alias{extract_replicate_term_to_matrix}
\title{Extract a term and replicates it randomN times for each linear model}
\usage{
extract_replicate_term_to_matrix(object, term)
}
\arguments{
\item{object}{An object of class \link[=cumhist-class]{cumhist}}

\item{term}{String, term name}
}
\value{
Matrix
}
\description{
Extract a term and replicates it randomN times for each linear model.
Used for population mean or variance terms.
}
\examples{
\donttest{
br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
bH_mu <- extract_replicate_term_to_matrix(br_fit, "bH_mu")
}
}
