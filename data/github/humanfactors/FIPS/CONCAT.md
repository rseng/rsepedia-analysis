[![Travis build status](https://travis-ci.com/humanfactors/FIPS.svg?branch=master)](https://travis-ci.com/humanfactors/FIPS)
[![codecov](https://codecov.io/gh/humanfactors/FIPS/branch/master/graph/badge.svg)](https://codecov.io/gh/humanfactors/FIPS)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02340/status.svg)](https://doi.org/10.21105/joss.02340)

# Fatigue Impairment Prediction Suite (FIPS)

<img align="right" src="https://github.com/humanfactors/FIPS/blob/master/inst/logo/FIPS_logo.png?raw=true" alt="FIPSLOGO" width="200"/> 

> If you are measure sleep behaviour or want to predict fatigue, this package probably can help.

FIPS is an R package that provides researchers and practitioners with a comprehensive set of functions for applying and simulating from bio-mathematical models (BMMs) of fatigue.

FIPS includes a set of functions for transforming sleep and actigraphy data to the data frame structure required for executing BMM simulations (called a [`FIPS_df`](https://humanfactors.github.io/FIPS/reference/FIPS_df.html)). Importantly, FIPS includes a set of functions for simulating from and interpreting several forms of BMM, including the [Unified Model](https://www.sciencedirect.com/science/article/pii/S0022519313001811) and [Three Process Model](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0108679). All models are extendable and include customisable parameters. The core features of FIPS leverage R's flexible `S3` class system, making extensions straightforward.

## Installation
We have no plans for a CRAN submission. To install the latest version of FIPS:

```r
# install.packages('remotes') # if remotes not installed
remotes::install_github("humanfactors/FIPS")
```

# Core Features Example

Detailed information regarding the FIPS data formats can be found in the ["FIPS Simulation Walkthrough Vignette"](https://humanfactors.github.io/FIPS/articles/FIPS-simulation-walkthrough.html).

If you require a more manual user-interface to generate the FIPS-df, we have a web application that may help you generate the input CSV: https://github.com/Akhil-Eaga/FIPS-Data-Entry-System.


**Step 1:** Prior to simulation, FIPS requires sleep history data to be in a special format, called a [`FIPS_df`](https://humanfactors.github.io/FIPS/reference/FIPS_df.html) which contains all the information required for modelling (e.g., time awake, time asleep). This can be created with [`parse_sleepwake_sequence`](https://humanfactors.github.io/FIPS/reference/parse_sleepwake_sequence.html) or [`parse_sleeptimes`](https://humanfactors.github.io/FIPS/reference/parse_sleeptimes.html).

```r
my_FIPS_dataframe = FIPS::parse_sleepwake_sequence(
  seq = unit_sequence,  # A binary (1,0) vector 
  epoch = 5,            # Epoch in minutes of vector
  series.start = as.POSIXct("2020-05-21 08:00:00"))
```

**Step 2:** To run a model simulation, you use [`FIPS_simulate`](https://humanfactors.github.io/FIPS/reference/FIPS_simulate.html), which returns a `FIPS_simulation` with all model predictions/forecasts generated in the corresponding columns. Note that the formula argument is optional, and sensible defaults will be used if omitted.

```r
# Run a simulation with the three process model
TPM.simulation.results = FIPS::FIPS_simulate(
  FIPS_df = my_FIPS_dataframe,   # A FIPS_df
  modeltype = "TPM",             # Three Process Model
  pvec = TPM_make_pvec()         # Default parameter vector
  # formula = alertness ~ s + c + u + w  # An optional formula for output
)

# Run a simulation with the unified model
unified.simulation.results = FIPS::FIPS_simulate(
FIPS_df = my_FIPS_dataframe,   # A FIPS_df
modeltype = "unified",         # Unified model
pvec = unified_make_pvec()     # Default parameter vector
)  
```

```
$ print(TPM.simulation.results)
> ---------
> Model Type: TPM 
> Epoch Value: 5 minutes 
> Simulation duration: 60 hours 
> Time points: 721 
> Parameters used (pvec input): ...[Suppressed for README]...
> For descriptions of these parameters, inspect:  help(FIPS::TPM_make_pvec) 
> ---------
> # A tibble: 721 x 10
>    datetime             time wake_status sim_hours     s     c     w        u   KSS alertness
>    <dttm>              <dbl> <lgl>           <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl>     <dbl>
>  1 2020-05-21 08:00:00  8    TRUE           0       7.96 -1.67 -5.72 -0.00274 10.3       6.28
>  2 2020-05-21 08:05:00  8.08 TRUE           0.0833  7.94 -1.63 -5.04 -0.00549  9.84      6.31
>  3 2020-05-21 08:10:00  8.17 TRUE           0.167   7.93 -1.59 -4.45 -0.00919  9.47      6.33
>  4 2020-05-21 08:15:00  8.25 TRUE           0.25    7.91 -1.55 -3.92 -0.0138   9.14      6.35
>  5 2020-05-21 08:20:00  8.33 TRUE           0.333   7.89 -1.50 -3.46 -0.0194   8.85      6.37
>  6 2020-05-21 08:25:00  8.42 TRUE           0.417   7.88 -1.46 -3.05 -0.0258   8.59      6.39
>  7 2020-05-21 08:30:00  8.5  TRUE           0.5     7.86 -1.42 -2.69 -0.0332   8.36      6.41
>  8 2020-05-21 08:35:00  8.58 TRUE           0.583   7.85 -1.37 -2.37 -0.0415   8.16      6.43
>  9 2020-05-21 08:40:00  8.67 TRUE           0.667   7.83 -1.32 -2.09 -0.0506   7.98      6.46
> 10 2020-05-21 08:45:00  8.75 TRUE           0.75    7.81 -1.28 -1.84 -0.0606   7.82      6.48
> # ... with 711 more rows
```

```
$ print(unified.simulation.results)
> ---------
> Model Type: unified 
> Epoch Value: 5 minutes 
> Simulation duration: 60 hours 
> Time points: 721 
> Parameters used (pvec input): ...[Suppressed for README]...
> For descriptions of these parameters, inspect:  help(FIPS::unified_make_pvec) 
> ---------
> # A tibble: 721 x 9
>    datetime             time wake_status sim_hours      s      l     c     w fatigue
>    <dttm>              <dbl> <lgl>           <dbl>  <dbl>  <dbl> <dbl> <dbl>   <dbl>
>  1 2020-05-21 08:00:00  8    TRUE           0      0      0      0.335 1.14     1.38
>  2 2020-05-21 08:05:00  8.08 TRUE           0.0833 0.0502 0.0206 0.320 1.10     1.37
>  3 2020-05-21 08:10:00  8.17 TRUE           0.167  0.100  0.0412 0.305 1.06     1.36
>  4 2020-05-21 08:15:00  8.25 TRUE           0.25   0.150  0.0618 0.291 1.02     1.35
>  5 2020-05-21 08:20:00  8.33 TRUE           0.333  0.200  0.0824 0.276 0.978    1.34
>  6 2020-05-21 08:25:00  8.42 TRUE           0.417  0.250  0.103  0.261 0.941    1.33
>  7 2020-05-21 08:30:00  8.5  TRUE           0.5    0.300  0.123  0.247 0.906    1.32
>  8 2020-05-21 08:35:00  8.58 TRUE           0.583  0.349  0.144  0.232 0.872    1.31
>  9 2020-05-21 08:40:00  8.67 TRUE           0.667  0.399  0.164  0.218 0.839    1.30
> 10 2020-05-21 08:45:00  8.75 TRUE           0.75   0.448  0.185  0.204 0.807    1.29
> # ... with 711 more rows
```

**Step 3:** You now can access printing, summary and plot methods for the FIPS_simulation object. A [detailed tutorial of the plotting functionality](https://humanfactors.github.io/FIPS/articles/plotting.html) for FIPS is provided in the vignettes.

```r
plot(TPM.simulation.results)
summary(TPM.simulation.results)
print(TPM.simulation.results)
```

## What are BMMs?

BMMs are a class of biological phenomenological models which are used to predict the neuro-behavioural outcomes of fatigue (e.g., alertness, performance) using sleep-wake history. There are several different BMM implementations, but most have their roots in Borbély's (1982) two process model which stipulates that sleepiness/performance impairment functions in response to two processes: a circadian process and a homeostatic process. BMMs enable hypothesis testing of the latent factors underlying the relationships between sleep, fatigue, and human performance. For example, they enable researchers to estimate the relative contributions of homeostatic processes on fatigue, relative to endogenous circadian processes. These models are also frequently applied by defence and industrial sectors to support system safety as part of broader fatigue management strategies. FIPS is the first open-source BMM framework enabling practitioners to inspect, validate, and ideally extend BMMs. 

## Critical Model Information

Critical functions, arguments and outputs for each model are summarised in table below:

| Model                        | Three Process           | Unified                      |
|------------------------------|-------------------------|------------------------------|
| `FIPS_Simulate` argument     | `"TPM"`                 | `"unified"`                  |
| Parameter vector function    | `TPM_make_pvec()`       | `unified_make_pvec`          |
| Default formula              | `alertness ~ s + c + u` | `fatigue ~ s + I(κappa) * c` |
| Time varying parameter names | s, c, u, w              | s, l, c, w                   |
| Default prediction outputs   | KSS, alertness          | Fatigue (PVT lapses)         |

# Contributing and Support

We welcome contributions great or small from the community. It would incredibly useful to receive feedback via Github Issues for anything regarding the package, including: installation issues, bugs or unexpected behaviour, usability, feature requests or inquiries, or even something you don't understand in the tutorials about this class of models more generally. Please file a Github issue for any general support queries too.

# Terms for Academic Usage
In addition to the rights stipulated in the GNU Affero GPL-3, we request that all academic work leveraging FIPS provide a direct citation to the software package.

Wilson et al., (2020). FIPS: An R Package for Biomathematical Modelling of Human Fatigue Related Impairment. _Journal of Open Source Software_, 5(51), 2340, https://doi.org/10.21105/joss.02340

```tex
@article{Wilson2020,
  doi = {10.21105/joss.02340},
  url = {https://doi.org/10.21105/joss.02340},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {51},
  pages = {2340},
  author = {Michael David Wilson and Luke Strickland and Timothy Ballard},
  title = {FIPS: An R Package for Biomathematical Modelling of Human Fatigue Related Impairment},
  journal = {Journal of Open Source Software}
}
```
# Contribution Guidelines

Thank you for showing interest in FIPS! We are excited to know you are interested in being apart of this community.

We welcome contributions great or small from the community. Please do note that given the project is still in infancy,
we would appreciate if you could contact us prior to submitting a large pull request. This is to ensure that we aren't
doubling up on development efforts, and also we would love to hear from you!

Currently we are a young project, maintained by a small team (_n_ = 3). Our primary focus over 2020 is expanding the
FIPS user base. It would _incredibly_ useful to receive feedback via Github Issues for anything regarding the package,
including: installation issues, bugs or unexpected behaviour, usability, feature requests or inquiries, or even
something you don't understand in the tutorials about this class of models more generally. See specific details for each
of these requests below:

**Reporting a Bug:** If something has gone wrong with FIPS, would are grateful for your time to report the issue. Please
use the provided "bug report template" when filing an issue. Please try to provide a minimal reproducible example (e.g.,
via the suggested `reprex` package) where possible.

**General BMM Questions:** If you have any general questions about the biomathematical models implemented in this
package feel free to open a issue. We understand BMMs are a relatively novel class of models, and real research and
practice circumstances often demand solutions beyond existing papers. Please do read the package documentation first,
but it is useful for us to know if more information is required.

**Data Importing:** A particular area where contributions would be appreciated is sleep data formats. If you have a
relatively standard sleep data in a format that you are unable to import into FIPS, please do file an issue. While we
cannot promise the resources to ensure your exact format is supported in FIPS, we will endeavor to do so, and warmly
welcome collaborations aimed at expanding the supported data types.

**Feature Requests:** We also welcome feature requests, but would kindly ask for citations/references to any requests
for technical implementations.

## Pull Requests & Code Conventions

As mentioned, we would appreciate if you could contact us prior to submitting a large pull request. In particular, if
you are interested in adding support for another form of BMM, please do reach out as this would prompt us to develop
guidelines for adding models.

A brief list of coding conventions, which are somewhat atypical:

- Hard wrap source files to a 120 line margin (excluding comments or strings senstive to newlines)
- Markdown files should be either unwrapped or 'filled' to 100 or 120 line margin in Emacs.
- We use `<-` for function assignment, and `=` for variable assignment.
- Underscores (i.e., snake-case) should separate words in object names. Period separation is currently used for legacy supported functions.
- When naming objects, abbreviations and parameter names should use capitals as a first priority (e.g., `TPM_make_pvec` not `tpm_make_pvec`; but `unified_make_pvec` is fine).

Development should be conducted on a separate branch. Tests are run using the `testthat` package automatically with Rstudio default.
This can be invoked by pressing `ctrl + shift + T` in Windows or Linux. `roxygen2` is used for all documentation.

Thankyou!---
title: 'FIPS: An R Package for Biomathematical Modelling of Human Fatigue Related Impairment'
tags:
  - R
  - psychobiology
  - human factors
  - fatigue
  - biomathematical modelling
  - dynamic models
authors:
  - name: Michael David Wilson
    orcid: 0000-0003-4143-7308
    affiliation: 1
  - name: Luke Strickland
    orcid: 0000-0002-6071-6022
    affiliation: 1
  - name: Timothy Ballard
    orcid: 0000-0001-8875-4541
    affiliation: 2
affiliations:
 - name: Future of Work Institute, Curtin University
   index: 1
 - name: School of Psychology, University of Queensland
   index: 2
date: May 01 2020
bibliography: sleep.bib
---

# Summary

In many workplace contexts, accurate predictions of a human's fatigue
state can drastically improve system safety. Biomathematical models of
fatigue (BMMs) are a family of dynamic phenomenological models that
predict the neurobehavioural outcomes of fatigue (e.g., sleepiness,
performance impairment) based on sleep/wake history [@Dawson2017].
However, to-date there are no open source implementations of BMMs, and
this presents a significant barrier to their broadscale adoption by
researchers and industry practitioners.

`FIPS` is an open source R package [@R] to facilitate BMM research and
simulation. FIPS has implementations of several published
bio-mathematical models and includes functions for easily manipulating
sleep history data into the required data structures. FIPS also includes
default plot and summary methods to aid model interpretation. Model
objects follow tidy data conventions [@wickham2014tidy], enabling FIPS to be
integrated into existing research workflows of R users.

# Background on Biomathematical Models

Borbély's [-@Borbely1982] seminal two-process BMM specifies that
subjective fatigue is modulated by the additive interaction of two
biological processes: the homeostatic and the circadian. The homeostatic
process, denoted by *S*, is responsible for the increase in fatigue
during wake and the recovery from fatigue during sleep. Fluctuations in
process *S* are described by exponential functions with fixed lower and
upper asymptotes. The endogenous circadian process, denoted by *C*,
reflects the effect of the body clock on sleep propensity. The dynamics
of these processes are driven by a set of governing parameters (e.g.,
the phase of the circadian process). Figure 1 below shows the additive
effects of varying governing parameters of S and C. This model has
formed the basis of many other models that predict neurobehavioural
performance and fatigue based on sleep history 
[@Akerstedt2008; @peng18_improved; @ramakrishnan16_unified_model; @hursh_fatigue_2004].

![A parameter sensitivity plot of the Three Process Model. The _x_ axis
represents a 24 hour day, with the dark gray plot regions indicating sleep and
the light gray indicating wake. The top panel shows the homeostatic process with
five variations of the $\tau_{d}$ parameter and the centre panel shows the
circadian process with five variations of the $\varphi_{phase}$ parameter. The
bottom panel shows the multiplicative combinations of all unique S and C
processes from the previous panels. The plot was produced with functionality in
`FIPS`.](FIPS_Fatigue_2PM.png)

BMMs have a rich history of application in laboratory sleep deprivation studies
where they are used to understand the latent factors underlying human fatigue.
An important aim of these studies is to identify the governing parameter values
which provide the best account for the data at hand [@Reifman2004]. Further,
propriety BMM implementations are frequently applied in safety-critical
industries (e.g., aviation, military, mining operations) to support system
safety as part of broader fatigue management strategies 
[e.g., aiding in rostering decisions; see @Dawson2017; @dawson_modelling_2011].

Unfortunately, the broader adoption of BMMs by the cognitive and behavioral
sciences has been constrained by several factors. Firstly, BMM researchers have
typically only provided the mathematical derivations of their models (i.e.,
formulae) and not their computational implementations. This is a barrier to
reproducibility [@wilson2019all] because implementing BMMs and the required data
structures from the ground up requires substantial expertise and time
investment. Prior to FIPS, the only available BMM implementations were contained
within closed-source commercial software [e.g., SAFTE-FAST; @hursh_fatigue_2004].
Even for researchers and practitioners able to afford licenses, these tools
prohibit users from inspecting, modifying, independently evaluating, or
extending the code and contained models.

# Package Motivation and Features

`FIPS` aims to make BMM approaches accessible to a wider community and
assist researchers in conducting robust and reproducible biomathematical
modelling of fatigue. The package includes:

-   Functions to transform common sleep data formats into the
    longitudinal data structure required for conducting BMM simulation
    and estimation.

-   Well documented implementations of three forms of BMM: The Unified
    Model [@ramakrishnan16_unified_model] and Two- and Three-Process Models
    [@Akerstedt2008; @Reifman2004; @Borbely1982].

-   A function for plotting BMM outputs, including with observed data
    points. The visualisations are publication-ready, but flexibly
    adjusted via the `ggplot2` package.

The package also contains two vignettes: a walk-through of a sleep simulation
scenario which includes generating, transforming and analyzing the data;
and a detailed tutorial in plotting model outputs.

# FIPS Interface and Data Structures

Conducting a BMM simulation in FIPS requires users to generate a `FIPS_df`, a
tidy data frame containing a time series (based on sleep history data) of all
variables required to conduct BMM research. FIPS supports two sleep data
formats, each format is associated with a corresponding function that automatically
performs all required transformations to the `FIPS_df` format:

-   The `parse_sleeptimes` function transforms a data frame containing three
    vectors: sleep onset times, sleep end times (i.e., awakening), and the sleep
    episode sequence identifier. This format is human readable and well suited
    for individuals who are manually entering sleep history data (e.g., from a
    paper sleep diary).[^1]	

-   The `parse_sleepwake_sequence` function transforms a bit vector
    representing sleep (0) and wake (1) statuses, with each bit
    representing an equal epoch (e.g., 1 minute). While not very human
    readable, this data format is commonly output by actigraphy devices
    and their corresponding sleep/wake status algorithms (e.g.,
    Cole--Kripke). Consequently, this format is often supported by other
    proprietary BMM software.

[^1]: It should be noted that the `parse_sleeptimes` function may also be useful to users wishing to convert a set of human-readable sleep times (e.g., manually entered from a sleep diary) to a bit vector (i.e., a list of `1` or `0` values). This is because the `FIPS_df` contains a column representing the series as a bit vector (see `FIPS_df$wake_status_int`).

The resulting `FIPS_df` dataframe can then be sent to the simulation
dispatch function `FIPS_simulate` to execute a specific model. This
function requires: a `FIPS_df`, a model string (e.g., "unified"), and an associated
parameter vector for the selected model (`pvec`). Documentation is
provided for customizing each parameter in the `pvec`, with citations
for the default values. The returned `FIPS_simulation` data frame has
added columns for each time-step of the series, including model
predictions (e.g., alertness and sleepiness in the case of the three
process model), as well as time-varying model processes (e.g., the
circadian process, *c*, and homeostatic process, *s*). `FIPS_simulation`
predictions can be plotted by calling plot or `FIPS_plot` on the object,
and simulation and model configurations can be reviewed via `summary`
or `print`.

# Research & Future Development

FIPS is being actively developed as part of a broader research project
examining cognitive fatigue prediction in safety-critical workplace
environments. Given the extensibility of FIPS's implementation, we hope
other BMM researchers may consider collaborating to implement other BMMs
in the framework.

# Acknowledgments

The development of this package was supported by the Future of Work Institute at Curtin University and the Maritime Division of the Australian Defence Science and Technology Group (DSTG).

# License

This project is licensed under the "GNU Affero General Public License"
version 3 - see the LICENSE file for details

# References

---
name: Bug report
about: Create a bug report to improve the package
title: "[BUG]"
labels: Bug
assignees: ''

---

**Bug description and expected behaviour**
_Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, use the question template._

---

**How to reproduce:**
_Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>._


```r
# insert reprex here
```

---

**System Information**

```
# insert output of sessionInfo() here
```
---
title: "Plotting"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Plotting}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, include=FALSE}
library(FIPS)
library(dplyr)
library(tibble)
library(ggplot2)
library(lubridate)
```

## Introduction

In this vignette we briefly demonstrate how to use the FIPS
plotting function. Below we create sleeptimes and simulate an example 
sleep scenario to plot. This example includes initial sleep deprivation.

(Note: you should read `vignette("FIPS-simulation-walkthrough","FIPS")` first).

```{r}
# Simulation start date time (i.e., when you want first predictions to begin)
simulation.start = lubridate::ymd_hms('2018-05-01 11:00:00', tz = "Australia/Perth")
# Simulation end date time (i.e., when you want predictions to end)
simulation.end = lubridate::ymd_hms('2018-05-09 21:00:00', tz = "Australia/Perth")
# Sleep times, 5 hours per night, starting 23:00
example.sleeptimes <- tibble::tibble(
  sleep.start = seq(
    from = lubridate::ymd_hms('2018-05-03 02:00:00', tz = "Australia/Perth"), 
    to = lubridate::ymd_hms('2018-05-09 16:00:00', tz = "Australia/Perth"),
    by = '24 hours'),
  sleep.end = sleep.start + lubridate::dhours(3),
  sleep.id = rank(sleep.start))

simulated.dataframe = parse_sleeptimes(
  sleeptimes = example.sleeptimes,
  series.start = simulation.start,
  series.end = simulation.end,
  sleep.start.col = "sleep.start",
  sleep.end.col = "sleep.end",
  sleep.id.col = "sleep.id",
  roundvalue = 5)

unified.simulation.results = FIPS_simulate(
  FIPS_df = simulated.dataframe, # The FIPS_df
  modeltype = "unified",         # three process model
  pvec = unified_make_pvec()      # paramater vector
  )

```

FIPS contains a generic `plot` method for simulation results. By default this will
plot fatigue over time. Sleep is indicated by the grey rectangles, and the start
of each new day is indicated by the dashed lines.

```{r fig.width=7}
plot(unified.simulation.results)
```

Often it will be helpful to restrict the plotting between certain datetimes,
and plot additional data such as the model's time-varying process estimates (e.g., c, s):

```{r fig.width=7}
plot(unified.simulation.results, from='2018-05-03 13:30:00',
                     to = "2018-05-07 18:30:00", plot_stat=c("fatigue", "s", "c", "l"))
```

Observed data can be plotted against the sleep predictions. Below demonstrates
by sampling some data points from the previously generated data frame,
then adding some random noise to make some 'pretend' observed data.


```{r fig.width=7}

pretend_ratings <- unified.simulation.results %>% dplyr::filter(wake_status) 

pretend_ratings <- pretend_ratings[sample(1:length(pretend_ratings$datetime), 100),
                   c("datetime", "sim_hours", "fatigue")]

pretend_ratings$fatigue <- pretend_ratings$fatigue + rnorm(length(pretend_ratings$fatigue),0,0.7)

plot(unified.simulation.results, from='2018-05-03 13:30:00',
                     to = "2018-05-07 18:30:00", plot_stat=c("fatigue", "s", "c", "l"),
                     fatigue_CIs = FALSE, observed=pretend_ratings, observed_y="fatigue")

```


Below demonstrates an example plotting the three process model. Results are
similar except the output is alertness, rather than fatigue. This is
automatically handled for you if you are using the generic `plot` function.

```{r fig.width=7}
TPM.simulation.results = FIPS_simulate(
  FIPS_df = simulated.dataframe, # The FIPS_df
  modeltype = "TPM",         # three process model
  pvec = TPM_make_pvec()      # paramater vector
  )

plot(TPM.simulation.results)
```

It is important to note that the generic `plot()` function is simply calling `FIPS_plot`
under the hood. Users wishing full control may find the `FIPS_plot` function more intuitive to work with.

```{r, fig.width=7}
FIPS_plot(TPM.simulation.results, 
          plot_stat = c("alertness", "c", "u", "s")) +
  ggplot2::ggtitle("Three Process Model Simulation") +
  ggplot2::scale_color_brewer(palette="Set1")

```

---
title: "FIPS Simulation Walkthrough"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{FIPS Simulation Walkthrough}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

This vignette demonstrates how to generate a sleep times dataframe and then proceed to process this for modeling with the Three Process Model of Fatigue. The code here demonstrates the usage from a higher-level user perspective. 

Note that in future releases, we may seek to create functions that explicitly generate these sleep times for you. Furthermore, it is possible to simply create the sleeptimes in Excel or another program and transform to the format shown here. Therefore, the sleep generation features are presented here for practical purposes.

```{r setup, include=FALSE}
library(FIPS)
library(dplyr)
library(tibble)
library(ggplot2)
library(lubridate)
library(colorspace)
library(tidyr)
```

## Data Preparation

### Sleep Times Format

The dataframe below shows a prototypical FIPS 'sleep times' dataframe. This form of dataframe is intended for individuals who are manually inputting sleep history (e.g., from pencil forms to a spreadsheet). Below, we generate this data structure ourselves for convenience. Specifically, we want to simulate a scenario where an individual obtains exactly 7.5 hours of sleep per night for 6 nights. There are likely multiple ways to achieve this, but the generation function below allows sufficient flexibility to have continuously offset sleeptimes (e.g., by changing `by = '24 hours'` to another value).

The sleeptimes should correspond to only one participant (multiple people are not currently supported in FIPS without map functions). The default column names for this dataframe are `sleep.id`, `sleep.start`, and `sleep.end`. It is recommended that you explicitly specify timezones in all functions to avoid any silent errors relating to time zone specifications.

```{r}
example.sleeptimes <- tibble::tibble(
  sleep.start = seq(
    from = lubridate::ymd_hms('2018-05-01 23:00:00', tz = "Australia/Perth"), 
    to = lubridate::ymd_hms('2018-05-07 17:00:00', tz = "Australia/Perth"),
    by = '24 hours'),
  sleep.end = sleep.start + lubridate::dhours(7.5),
  sleep.id = rank(sleep.start))

print(example.sleeptimes)
```

Prior to actually conducting the simulation, you _must_ convert the sleep times to the continuous "FIPS_df" format. This format is a continuous time series style dataframe that contains calculated variables to be interpreted by the FIPS model functions. This dataframe contains the following headings: `datetime, sleep.id, wake_status, wake_status_int, change_point, switch_direction, status_duration, total_prev, time, day, sim_hours`. Information regarding these will be presented in the section below, but first let's quickly run through how to generate the FIPS dataframe from sleep times.

The `parse_sleeptimes` function from FIPS will takes in 'sleep times' generated previously, and note the arguments below.

- `sleeptimes` — The sleep times dataframe generated previously
- `series.start` — This is the start datetime of the entire simulation series.
- `series.end` — This will be the datetime of the entire simulation series
- `sleep.start.col` — This is the name of your `sleep.start` column in sleep times if changed from default above
- `sleep.end.col` — This is the name of your `sleep.end` column in sleep times if changed from default above
- `sleep.id.col` — This is the name of your `sleep.id` column in sleep times if changed from default above
- `roundvalue` — The epoch of the series (i.e., rounding value). Please read further information below!

*The setting of `roundvalue` is critically important.* It determines the epoch or spacing between observations in your dataset (in minutes). At a value of `1`, the simulation is updated every 1 minute. Consequently, this will increase the size (in rows) of your dataset by a factor of 5 (relative to `5` minutes). In most cases, 5, 10 or even 15 minutes should be sufficient.

Moreover, note that all sleep observations will have their datetime rounded to this value. For example, with a `roundvalue = 5` the datetime `2018-05-07 21:02:02` would be rounded to `2018-05-07 21:00:00`. For this reason, it is ideal if your `series.start` and `series.end` are rounded to the same epoch value as you request. Seconds should never be included in your datetime (biomathematical models are not sensitive to this time resolution anyway).

```{r}
# Simulation start date time (i.e., when you want first predictions to begin)
simulation.start = lubridate::ymd_hms('2018-05-01 07:00:00', tz = "Australia/Perth")

# Simulation end date time (i.e., when you want predictions to end)
# In this case it ends 
simulation.end = lubridate::ymd_hms('2018-05-07 21:00:00', tz = "Australia/Perth")

# The Continuous FIPS_df dataframe format
# This creates the format ready for simulation
simulated.dataframe = parse_sleeptimes(
  sleeptimes = example.sleeptimes,
  series.start = simulation.start,
  series.end = simulation.end,
  sleep.start.col = "sleep.start",
  sleep.end.col = "sleep.end",
  sleep.id.col = "sleep.id",
  roundvalue = 5)

print(simulated.dataframe)
```
### Bitvector Sequence Format

It is common to represent sleep history information as a bitvector (i.e., a sequence of 1's and 0's). In bitvectors, sleep and wake times are represented by 1's or 0's, with each bit representing an equal time duration (e.g., 1 minute in that status). The bitvector sequence also must be relative to a *start datetime*. This format is commonly outputted from actigraphy devices (a form of wearable sleep tracker) and corresponding sleep detection algorithms (e.g., Cole-Kripke). Other proprietary software packages (e.g., SAFTE-FAST) require imported data to be in bit vector form (though the exact required formats do vary).

`FIPS` expects bitvectors to repesent **sleep as 0** and **wake as 1**. The `parse_sleepwake_sequence` function can transform a bit vector sequence to a compliant `FIPS_df`. Below, an example of this data and steps to transform are provided, however, note that we do not use this dataframe again within this vignette.

```{r}
# Simulation start date time (i.e., when you want first predictions to begin)
simulation.start = lubridate::ymd_hms('2018-05-01 10:00:00', tz = "Australia/Perth")

# Example bitvector sequence. This typically would be imported directly via a textfile.
# Here we generate, though typically this would be returned by a ReadLines/read.delim/read.table
bv.sleep.sequence = rep(rep(c(1,0), 6), sample(20:40, 12))

bv.sim.dataframe = parse_sleepwake_sequence(
  seq = bv.sleep.sequence,
  series.start = simulation.start,
  epoch = 15)

print(bv.sim.dataframe)
```

## Modelling and the Simulation Dataframe

Now that you have generated the FIPS format, you can now apply the FIPS simulation functions to this to actually run BMM predictions on the series. In the example below, we will run a Three Process Model simulation over the `FIPS_df` series we just created. To do this, we will use the `FIPS_simulation` function, which takes in three arguments: a `FIPS_df` object, a specification of a `modeltype` (see help for model types currently implemented), and a `pvec` which is a vector of parameters for the model.

The parameter vectors provided by default in the package are those reported by [Ingre et al. (2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0108679). Please see the help files for citations and further information. These defaults are customisable, and if you are using Rstudio or Emacs ESS, you will get full autocompletion with descriptions of each parameter (See `help("TPM_make_pvec", "FIPS")`).

We will use the default parameter vector for the Three process model created by `FIPS::TPM_make_pvec()`, which takes on the following values:

```{r results="asis", echo=FALSE}

pdf <- data.frame("Parameter" = names(TPM_make_pvec()),
                  "Value" = unname(TPM_make_pvec()))

knitr::kable(list(pdf[1:7,], pdf[8:15,]), row.names=FALSE, digits=3)
```

Calling `FIPS_simulate()` will the produces a `FIPS_df` with model predictions and predicted time-varying process variables (e.g., `s`, `c`). The returned `FIPS_df` will also now inherit the `FIPS_simulation` class.  A `FIPS_simulation` object has attributes containing the parameter vector, the `modeltype` string, and the `pvec` used, and several other values used internally. A custom print function of the object will reveal all this information. Note that this print function does mask some of the columns for ease of reading.

```{r}
# Run a simulation with the three process model
TPM.simulation.results = FIPS_simulate(
  FIPS_df = simulated.dataframe, # The FIPS_df
  modeltype = "TPM",             # three process model
  pvec = TPM_make_pvec()       # parameter vector
  )

TPM.simulation.results
```

- `datetime` = vector of datetime stamps separated at equidistance intervals.
- `sleep.id` = a supplementary variable indicating sleep episode identifier.
- `wake_status` =  Awake (`T`) or asleep (`F`) at that epoch interval/epoch
- `wake_status_int` = Awake (1) or asleep (0) at that epoch interval/epoch
- `change_point` = Whether the individual changed wake status at that interval/epoch.
- `switch_direction` = Whether switch was to sleep or to wake
- `status_duration` = How long individual has been in status at that current time point
- `total_prev` = If a switch has occured, how long were that in the previous status.
- `time` = time of day in decimal hour
- `day` =  days into simulation
- `sim_hours` = hours simulation has run for total
- `s` — Estimate of S (homeostatic) process at that time
- `l` — Estimate of l (lower asmytope) process at that time
- `c` — Estimate of C (24-hour circadian) process at that time
- `w` — Estimate of W (sleep intertia) process at that time
- `u` — Estimate of U (12-hour circadian) process at that time
- `alertness` — Currently just `s + c + u`
- `KSS` — This is equal to `10.6 + -0.6 * (alertness)`

# Plotting
### Plot 1: Parameter Plots

FIPS provides a default plot method to visualize model predictions and time-varying process estimates over time, discussed in detail
[in the plotting vignette](plotting.html). These plots aid with debugging and understanding how the different model parameters contribute to predictions. 

```{r fig.width=7}
# Plot the whole time series
plot(TPM.simulation.results, plot_stat=c("alertness", "s", "c"))

#Narrow in on the first 90 hours
plot(TPM.simulation.results, plot_stat=c("alertness", "s", "c"),
      from= '2018-05-01 23:00:00',
      to= as_datetime('2018-05-01 23:00:00') + hours(90))

```

### Plot 2: Heatmap Prediction Plots

There are many other ways that visualizations could help inform theory and practice.
The below heat plot is an example of how dangerous times in a mission, according to the TPM predictions, 
could be visualized. In the example below, the regions in deep orange indicate a greatly increased risk of fatigue. However, the fact that the increased fatigue occurs at the start of the mission indicates the initalisation parameters (e.g., `S0`) are to blame, so this isn't a cause for significant concern.

```{r fig.width=5, fig.height=5}

if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package \"colorspace\" needed for this example to work. Please install it.",
      call. = FALSE)
  }

plot3_w = TPM.simulation.results %>%
  # Change epoch to 30 minutes for this plot (optional)
  mutate(time = strftime(datetime, format = "%H:%M")) %>% 
  separate(time, into = c("plot_hour", "plot_minute"), ":") %>%
  mutate(halfhour = if_else(plot_minute < 30, 0, 0.5) + as.numeric(plot_hour)) %>% 
  group_by(day, plot_hour, halfhour) %>% 
  summarise(KSS = mean(KSS)) %>% 
  filter(day > 1) %>% 
  # Start Plot
  ggplot(aes(day, halfhour, fill = KSS)) +
  geom_tile(aes(fill = KSS), color = "white") +
  scale_x_continuous(breaks = seq(0,7,1)) +
  scale_y_continuous(breaks = seq(0,23,1)) +
  labs(y = "Hour of Day", x = "Day of Mission") +
  colorspace::scale_fill_continuous_divergingx(name = "KSS", palette = "Zissou 1", limits = c(1,9), mid = 5) +
  ggtitle(label = "Three Process Model (KSS)", subtitle = "Heatmap of TPM predicted performance by 0.5hr") +
  theme(axis.title = element_text(size = 12)) 

plot3_w

```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_formula_interface.R
\name{get_bmm_model_frame}
\alias{get_bmm_model_frame}
\title{get_bmm_model_frame}
\usage{
get_bmm_model_frame(.FIPS_sim, model_formula)
}
\arguments{
\item{.FIPS_sim}{A FIPS_Simulation object}

\item{model_formula}{A formula describing how the time-varying processes predictors should be calculated for the predicted output.}
}
\value{
Minimal dataframe with only requird columns for computation
}
\description{
get_bmm_model_frame
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_make_pvec}
\alias{TPM_make_pvec}
\title{Make Three Process Model (TPM) Parameter Vector}
\usage{
TPM_make_pvec(
  la = 2.4,
  ha = 14.3,
  d = -0.0353,
  g = log((14.3 - 14)/(14.3 - 7.96))/8,
  bl = 12.2,
  Cm = 0,
  Ca = 2.5,
  p = 16.8,
  Um = -0.5,
  Ua = 0.5,
  Wc = -5.72,
  Wd = -1.51,
  S0 = 7.96,
  KSS_intercept = 10.6,
  KSS_beta = -0.6
)
}
\arguments{
\item{la}{Low asymptote. Minimum alertness allowed by homeostatic process (S)}

\item{ha}{High asymptote. Maximum alertness allowed by homeostatic process (S)}

\item{d}{Alertness decay rate. Rate at which alterness decays when awake (homestatic process)}

\item{g}{Fatigue recovery rate per unit time. Rate at which alertness recovers when asleep}

\item{bl}{Alertness level that breaks Sprime function. The alertness level at which low pressure sleep kicks in}

\item{Cm}{Mesor of C process (average level) for 24-hour circadian (C) process}

\item{Ca}{Amplitude of C process (extent to which peaks deviate from average level) for 24-hour circadian (C) process}

\item{p}{Default C process phase (i.e., peak) for 24-hour circadian (C) process}

\item{Um}{Mesor of U process (average level) for 12-hour circadian (U) process (dip in afternoon)}

\item{Ua}{Amplitude of U process process (extent to which peaks deviate from average level) for 12-hour circadian (U) process (dip in afternoon)}

\item{Wc}{Initial reduction in alertness for sleep inertia (W) process}

\item{Wd}{Recovery rate for sleep inertia (W) process}

\item{S0}{Initial value of homeostatic process (S)}

\item{KSS_intercept}{KSS transformation intercept}

\item{KSS_beta}{KSS transformation beta}
}
\value{
parameter vector
}
\description{
Creates and checks a TPM parameter vector. No arguments returns default settings from
Ingre, M., Van Leeuwen, W., Klemets, T., Ullvetter, C., Hough, S., Kecklund, G., Karlsson, D., & Åkerstedt, T. (2014). Validating and Extending the Three Process Model of Alertness in Airline Operations. \emph{PLoS ONE}, \emph{9}(10), e108679. \url{https://doi.org/10.1371/journal.pone.0108679}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Sp2fun}
\alias{TPM_Sp2fun}
\title{TPM_Sp2fun}
\usage{
TPM_Sp2fun(ha, g, bl, tas, breaktime)
}
\arguments{
\item{ha}{High asymptote (typically = 14.3)}

\item{g}{Rate of recovery (typically about -0.3813)}

\item{bl}{Break point (typically 12.2)}

\item{tas}{Time asleep (hours)}

\item{breaktime}{Break time (time from sleep until breakpoint is reached)}
}
\value{
S
}
\description{
Calculates S during lower pressure component of sleep after break point is reached
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_breaktimefun}
\alias{TPM_breaktimefun}
\title{Break time function (i.e., bt function)}
\usage{
TPM_breaktimefun(ha, g, bl, ss)
}
\arguments{
\item{ha}{High asymptote (typically = 14.3)}

\item{g}{Rate of recovery (typically about -0.3813)}

\item{bl}{Break point (typically 12.2)}

\item{ss}{S upon falling asleep}
}
\value{
breaktime
}
\description{
Break time function (i.e., bt function)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Ufun}
\alias{TPM_Ufun}
\title{TPM U Function (12-hour circadian process)}
\usage{
TPM_Ufun(Um, Ua, p, tod)
}
\arguments{
\item{Um}{Mesor of U process (typically = -0.5)}

\item{Ua}{Amplitude of U process (typically = 0.5)}

\item{p}{Default C phase (i.e., time of peak typically 16.8)}

\item{tod}{Time of day (in decimal hours)}
}
\value{
U
}
\description{
Calculates 12-hour circadian process.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_simulation_dispatch}
\alias{unified_simulation_dispatch}
\title{Unified Simulation Dispatcher}
\usage{
unified_simulation_dispatch(dat, pvec, model_formula)
}
\arguments{
\item{dat}{input dataframe (ensure this is a FIPS_df)}

\item{pvec}{a vector of default parameters, see [unified_pvec]}

\item{model_formula}{A formula expression object of desire model caluclation}
}
\description{
Constructor/dispatcher for Unified model simulations.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_simulate}
\alias{TPM_simulate}
\title{Simulate: Three Process Model}
\usage{
TPM_simulate(pvec, dat)
}
\arguments{
\item{pvec}{a vector of default parameters}

\item{dat}{input dataframe (ensure this is in FIPS format)}
}
\value{
dataframe with simulated values - where fatigue middle is estimate (if no error terms)
}
\description{
Simulates three process model over specified period.
Default parameters (la through S0) are constants used in previous applications of the model.
}
\details{
Access the modelled tibble directly by calling the object.
}
\section{Parameters}{


S function (homeostatic process)
la = low asymptote (default = 2.4)
ha = high asymptote (default = 14.3)
d = rate of decay in alertness when awake (default = -0.0353)
g = rate of recovery in alertness when asleep (default = log((ha-14)/(ha-7.96))/8)
bl = break level in alertness, time at which low pressure sleep kicks in (default = 12.2)

C function (24-circadian process)
Cm = average level of process (i.e., mesor; default = 0)
Ca = amplitude of process (default = 2.5)
p = phase or time at which process reaches its peak (default = 16.8)

U function (12-hour circadian process)
Um = average level of process (i.e., mesor; default = -0.5)
Ua = amplitude of process (default = 0.5)

W function (sleep intertia process)
Wc = Initial reduction in alertness at time of waking (default = -5.72)
Wd = Rate of recovery of alertness (default = -1.51)

Regression equation for converting alertness to KSS fatigue ratings
a = intercept (default = 10.6)
b = coefficient (default = -0.6)
}

\seealso{
TPM_make_pvec
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_simulation.R
\name{summary.FIPS_simulation}
\alias{summary.FIPS_simulation}
\title{summary.FIPS_simulation}
\usage{
\method{summary}{FIPS_simulation}(x)
}
\description{
summary.FIPS_simulation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Lfun}
\alias{unified_Lfun}
\title{Sleep Debt Penalty (L) Function during Wake
Calculates L / Sleep Debt Process during wake}
\usage{
unified_Lfun(l_at_wake, taw, tau_la, U0)
}
\arguments{
\item{l_at_wake}{Lower asymptote at wake onset}

\item{taw}{Time awake}

\item{tau_la}{Rate of change in lower asymptote}

\item{U0}{Upper asymptote}
}
\description{
Sleep Debt Penalty (L) Function during Wake
Calculates L / Sleep Debt Process during wake
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sleeptimesformat.R
\name{generate_presleep_times}
\alias{generate_presleep_times}
\title{Fill pre-observation wake times}
\usage{
generate_presleep_times(simulationstart, firstsleep, expand_by = 5)
}
\arguments{
\item{simulationstart}{start of simulation}

\item{firstsleep}{first sleep in the sleep dataframe}

\item{expand_by}{expand}
}
\value{
returns expanded tibble containing sleep.id = NA (due to waking) and wake_status = T
}
\description{
The first sleep is unlikely to also be the start of the mission simulation
Thus, this function fills the start of the tibble with the all times between
The mission start time and the first instance of sleep, intervaled by X minutes
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Sfun}
\alias{unified_Sfun}
\title{S Wake (Unified)
Calculates S during wake}
\usage{
unified_Sfun(s_at_wake, taw, tau_w, U0)
}
\arguments{
\item{s_at_wake}{S upon waking}

\item{taw}{Time awake}

\item{tau_w}{Controls rate of rise in S during wake}

\item{U0}{Upper asymptote}
}
\description{
S Wake (Unified)
Calculates S during wake
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Cfun}
\alias{TPM_Cfun}
\title{TPM C Function (24-hour circadian process)}
\usage{
TPM_Cfun(Cm, Ca, p, tod)
}
\arguments{
\item{Cm}{Mesor of C process (typically = 0)}

\item{Ca}{Amplitude of C process (typically = 2.5)}

\item{p}{Default C phase (i.e., time of peak typically 16.8)}

\item{tod}{Time of day (in decimal hours)}
}
\value{
C
}
\description{
Calculates 24-hour circadian process.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_formula_interface.R
\name{process_bmm_formula}
\alias{process_bmm_formula}
\title{process_bmm_formula}
\usage{
process_bmm_formula(.FIPS_sim, model_formula, pvec)
}
\arguments{
\item{.FIPS_sim}{A FIPS_Simulation object}

\item{model_formula}{A formula describing how the time-varying processes predictors should be calculated for the predicted output.}

\item{pvec}{A required pvec argument for the .FIPS_sim}
}
\value{
Minimal dataframe with only requird columns for computation
}
\description{
The pvec is required to ensure it is contained in the environment for expression evaluation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_plot.R
\name{plot.FIPS_simulation}
\alias{plot.FIPS_simulation}
\title{plot.FIPS_simulation}
\usage{
\method{plot}{FIPS_simulation}(
  x,
  from = NULL,
  to = NULL,
  plot_stat = NULL,
  fatigue_CIs = FALSE,
  observed = NULL,
  observed_y = NULL
)
}
\arguments{
\item{x}{A valid .FIPS_simulation series that has been simulated}

\item{from}{The starting datetime to be plotted}

\item{to}{The ending datetime to be plotted}

\item{plot_stat}{Which variables to plot}

\item{fatigue_CIs}{A logical indicating whether uncertainty intervals on fatigue should be plotted}

\item{observed}{A data frame with any observed sleepiness ratings or objective indicators to plot against predictions}

\item{observed_y}{The name of the observed sleepiness ratings in the observed data frame}
}
\description{
S3 plot method for FIPS_simulation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_simulation.R
\name{is_FIPS_simulation}
\alias{is_FIPS_simulation}
\title{Test if the object is a simmed FIPS_df}
\usage{
is_FIPS_simulation(x)
}
\arguments{
\item{x}{An object}
}
\value{
\code{TRUE} if the object inherits from the \code{inherits(x, "FIPS_df") & attr(x, "simmed") }.
}
\description{
This function returns \code{TRUE} for FIPS_df if a simulation has been run on it,
and \code{FALSE} for all other objects, including regular data frames.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sleeptimesformat.R
\name{round_times}
\alias{round_times}
\title{Round times by column}
\usage{
round_times(.data, colname, round_by = 5)
}
\arguments{
\item{.data}{The sleeptimes dataframe}

\item{colname}{the column required to be rounded}

\item{round_by}{Amount (in minutes) to round sleep times to}
}
\value{
The sleep dataframe with all sleep.start and sleep.end rounded to X minute interval
}
\description{
Round times by column
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\docType{data}
\name{pvec.threeprocess}
\alias{pvec.threeprocess}
\title{pvec.threeprocess
Here for compatability - will likely remove soon}
\format{
An object of class \code{numeric} of length 15.
}
\usage{
pvec.threeprocess
}
\description{
pvec.threeprocess
Here for compatability - will likely remove soon
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Lpfun}
\alias{unified_Lpfun}
\title{Sleep Debt Penalty (L) Function during Sleep
calculates L during sleep}
\usage{
unified_Lpfun(ls, tas, tau_la, U0)
}
\arguments{
\item{ls}{Lower asymptote upon falling asleep}

\item{tas}{Time asleep}

\item{tau_la}{Rate of change in lower asymptote}

\item{U0}{Upper asymptote}
}
\description{
Sleep Debt Penalty (L) Function during Sleep
calculates L during sleep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_singlevectorformat.R
\name{parse_sleepwake_sequence}
\alias{parse_sleepwake_sequence}
\title{Parse Sleep Wake Binary/Integer Sequence}
\usage{
parse_sleepwake_sequence(seq, epoch, series.start)
}
\arguments{
\item{seq}{A sequence of \code{0} (sleep) and \code{1} (wake) integers indicating sleep/wake status at that moment.}

\item{epoch}{Integer expressing length of each observations in the series (minutes).}

\item{series.start}{A POSIXct object indicating the start datetime of the simulation (i.e., pre-first sleep waking duration)}
}
\value{
\code{FIPS_df} formatted dataframe
}
\description{
It is common to have sleep wake history information in the form of a binary sequence.
This is the format used by SAFTE-FAST and other proprietary software.
Further, this format is often easily exported by actigraphy measurement software.
}
\examples{
 
start_date = as.POSIXct("2018-05-01 10:00:00")
bitvector_sequence = rep(rep(c(1,0), 6), sample(20:40, 12))
FIPSdf_from_bitvec = parse_sleepwake_sequence(
 seq = bitvector_sequence,
 series.start = start_date,
 epoch = 15)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Cfun}
\alias{unified_Cfun}
\title{Unified Circadian Process (C)
calculates C (circadian process)}
\usage{
unified_Cfun(tod, phi, tau = 24, A = 1)
}
\arguments{
\item{tod}{Time of day (in decimal hours)}

\item{phi}{Phase at beginning of the simulation (I think this should be 0 if t = tod)}

\item{tau}{Period of C process}

\item{A}{Amplitute of process}
}
\description{
Unified Circadian Process (C)
calculates C (circadian process)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Wfun}
\alias{unified_Wfun}
\title{Sleep inertia function (direct 3PM import)
Caclulates effect of sleep intertia on alterness}
\usage{
unified_Wfun(taw, wc, wd)
}
\arguments{
\item{taw}{Time awake}

\item{wc}{Extent of alertness reduction at time of waking (typically = -5.72)}

\item{wd}{Rate of recovery of alterness (typically = -1.51)}
}
\description{
Sleep inertia function (direct 3PM import)
Caclulates effect of sleep intertia on alterness
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sleeptimesformat.R
\name{parse_sleeptimes}
\alias{parse_sleeptimes}
\alias{sleeptimes_to_FIPSdf}
\title{Parse Sleep Times to FIPS_df}
\usage{
parse_sleeptimes(
  sleeptimes,
  series.start,
  series.end,
  roundvalue = 5,
  sleep.start.col,
  sleep.end.col,
  sleep.id.col
)

sleeptimes_to_FIPSdf(
  sleeptimes,
  series.start,
  series.end,
  roundvalue = 5,
  sleep.start.col,
  sleep.end.col,
  sleep.id.col
)
}
\arguments{
\item{sleeptimes}{A dataframe in the sleep time format (see help for more info)}

\item{series.start}{A POSIXct object indicating the start datetime of the simulation (i.e., pre-first sleep waking duration)}

\item{series.end}{A POSIXct object indicating the end datetime of the simulation}

\item{roundvalue}{An whole numeric (integer) value to round the sleep times to in minutes (\verb{default = 5 minutes}). Second precision not supported.}

\item{sleep.start.col}{\link{string} The column in the dataframe containing the sleep start times}

\item{sleep.end.col}{\link{string} The column name in the dataframe containing the sleep end times}

\item{sleep.id.col}{\link{string} A column name specifying the sleep id sequence (i.e., \code{1:n()})}
}
\value{
FIPS_df
}
\description{
This function parses a standardised sleeptime dataframe into the full FIPS format, ready for simulation and modelling.
The sleeptime format requires a sleep.id column (vector), a series of sleep times, and a series of corresponding wake times.
This format is the simplest to work with for human-readable or human-generated dataframes. See \link{parse_sleepwake_sequence} for
binary input methods.
}
\details{
It is crucial that that following conditions are met for all arguments:
\itemize{
\item Ensure that all specified datetimes for all datetime arguments are in an identical timezone.
\item Ensure that the minimum sleep start time is >= series.start
\item Ensure that the maximum wake time (sleep end) is <= series.end
\item Ensure that each sleep start is < the corresponding sleep.end
}
}
\examples{

 my_sleeptimes = tibble::tribble(
   ~sleep.id,          ~sleep.start,            ~sleep.end,
   1L, "2018-05-21 01:00:00", "2018-05-21 07:00:00",
   2L, "2018-05-21 23:00:00", "2018-05-22 04:00:00",
   3L, "2018-05-23 01:00:00", "2018-05-23 09:00:00") \%>\%
   dplyr::mutate(
     sleep.start = lubridate::ymd_hms(sleep.start),
     sleep.end = lubridate::ymd_hms(sleep.end))

 my_simstart = lubridate::ymd_hms('2018-05-20 22:00:00')
 my_simend   = lubridate::ymd_hms('2018-05-23 10:00:00')

 my_FIPS_df = parse_sleeptimes(
   sleeptimes = my_sleeptimes,
   series.start = my_simstart,
   series.end = my_simend,
   sleep.start.col = "sleep.start",
   sleep.end.col = "sleep.end",
   sleep.id.col = "sleep.id",
   roundvalue = 5)

}
\seealso{
For binary input parsing see: \link{parse_sleepwake_sequence}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_simulation.R
\name{print.FIPS_simulation}
\alias{print.FIPS_simulation}
\title{print.FIPS_simulation}
\usage{
\method{print}{FIPS_simulation}(x)
}
\description{
print.FIPS_simulation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_simulate.R
\name{FIPS_simulate}
\alias{FIPS_simulate}
\title{FIPS Simulation dispatcher}
\usage{
FIPS_simulate(FIPS_df, modeltype = NULL, pvec, model_formula = NULL)
}
\arguments{
\item{FIPS_df}{A valid FIPS_df series that has not already been modelled}

\item{modeltype}{String: either \code{"TPM"} (Three Process Model) or \code{"unified"}.}

\item{pvec}{Parameter vector (named list), see default pvecs for guidance.}

\item{model_formula}{An optional formula describing how the time-varying processes predictors should be calculated for the predicted output. See details.}
}
\value{
a FIPS_simulation object
}
\description{
\code{FIPS_simulate} is used to apply a particular BMM simulation to a \code{FIPS_df}.
It will dispatch to the selected model simulation type and then return the df with
the fitted model columns added on as a \code{FIPS_simulation} object.
}
\details{
If the formula argument is omitted, then default prediction values will be returned:
KSS and alertness for TPM, and fatigue (i.e., lapses) for unified.
\subsection{Formula Argument}{
\itemize{
\item The formula argument takes the form of an R formula expression (e.g., \code{y ~ c + s}).
\item The term of the left-hand side (e.g., \code{y}) defines the variable name
\item The term(s) on the right-hand side define how the variable is computed
\item All variables must be defined in your enviornment or in the FIPS_simulation output (e.g., \verb{s, c, l, w} for "TPM").
\item Parameter vector and other variable arguments must be placed in \code{I(expression)} format. For example,
\code{fatigue ~ s + c + I(pvec[["KSS_beta"]])}.
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_plot.R
\name{FIPS_plot}
\alias{FIPS_plot}
\title{FIPS Time Series Plot}
\usage{
FIPS_plot(
  dats,
  from = NULL,
  to = NULL,
  plot_stat = NULL,
  fatigue_CIs = FALSE,
  observed = NULL,
  observed_y = NULL
)
}
\arguments{
\item{dats}{A FIPS_simulation object (i.e., FIPS_df with simulation results)}

\item{from}{The starting datetime to be plotted}

\item{to}{The ending datetime to be plotted}

\item{plot_stat}{Which variables to plot}

\item{fatigue_CIs}{A logical indicating whether uncertainty intervals on fatigue should be plotted}

\item{observed}{A data frame with any observed sleepiness ratings or objective indicators to plot against predictions}

\item{observed_y}{The name of the observed sleepiness ratings in the observed data frame}
}
\value{
A ggplot2 object displaying fatigue and other requested processes over time
}
\description{
FIPS Time Series Plot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Spfun}
\alias{TPM_Spfun}
\title{TPM Sleep S Function}
\usage{
TPM_Spfun(ha, g, bl, ss, tas)
}
\arguments{
\item{ha}{High asymptote (typically = 14.3)}

\item{g}{Rate of recovery (typically about -0.3813)}

\item{bl}{Break point of sleep recovery (typically 12.2)}

\item{ss}{S at falling asleep}

\item{tas}{Time asleep (hours)}
}
\value{
S
}
\description{
Calculates S during sleep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sleeptimesformat.R
\name{expand_sleep_series}
\alias{expand_sleep_series}
\title{Expand Sleep Times to full vector}
\usage{
expand_sleep_series(.data, expand_by = 5)
}
\arguments{
\item{.data}{A sleeptimes dataframe}

\item{expand_by}{Amount (in minutes) to expand sleep times by}
}
\value{
Sleeptimedataframe with single columns vector for datetime and wake status
}
\description{
Turns the paired sleeptimes into a long single vectored datetime sequence
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_Spfun}
\alias{unified_Spfun}
\title{S Sleep (Unified)
Calculates S during sleep}
\usage{
unified_Spfun(ss, tas, tau_s, U0, tau_la, ls)
}
\arguments{
\item{ss}{S upon falling asleep}

\item{tas}{Time asleep}

\item{tau_s}{Controls rate of decay in S during sleep}

\item{U0}{Upper asymptote}

\item{tau_la}{Rate of change in lower asymptote}

\item{ls}{Lower asymptote at sleep onset}
}
\description{
S Sleep (Unified)
Calculates S during sleep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_check_pvec}
\alias{unified_check_pvec}
\title{Check Unified Parameter Vector}
\usage{
unified_check_pvec(pvec)
}
\arguments{
\item{pvec}{The pvec to check contains all required three process model parameters}
}
\value{
logical
}
\description{
Checks the pvec contains required parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_df.R
\name{is_FIPS_df}
\alias{is_FIPS_df}
\title{Test if the object is a \link{FIPS_df}}
\usage{
is_FIPS_df(x)
}
\arguments{
\item{x}{An object}
}
\value{
\code{TRUE} if the object inherits from the \code{inherits(x, "FIPS_df")}.
}
\description{
This function returns \code{TRUE} for FIPS_df,
and \code{FALSE} for all other objects, including regular data frames.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS-package.R
\docType{package}
\name{FIPS-package}
\alias{FIPS}
\alias{FIPS-package}
\title{FIPS: The Fatigue Impairment Prediction Suite}
\description{
The "Fatigue Impairment Prediction Suite" (FIPS) is currently under development and implemented in the R programming language.
FIPS provides researchers and practitioners comprehensive set of functions for applying bio-mathematical models (BMMs) of fatigue.
FIPS is under active development and implemented in the R programming language.
FIPS provides a set of well-documented functions for transforming sleep and actigraphy data to formats required for applying BMMs,
as well as a set of functions for simulating and interpreting BMMs with several kinds of models and customisable parameter settings.
}
\details{
BMMs are a class of biological phenomenological models which are used to predict the neuro-behavioural
outcomes of fatigue (e.g., alertness, performance) using sleep-wake history.
These models are frequently applied by defence and industrial sectors to support system
safety as part of broader fatigue management strategies.
FIPS is the first open-source BMM framework enabling practitioners to inspect, validate, and ideally extend BMMs.
Although there are several different implementations of BMMs, most have their roots in Borbély's (1982)
 two process model which models sleepiness/performance impairment as the sum of two processes: a circadian process and a homeostatic process.
}
\author{
\strong{Maintainer}: Michael David Wilson \email{michael.d.wilson@curtin.edu.au}

Authors:
\itemize{
  \item Luke Strickland \email{luke.strickland@curtin.edu.au}
  \item Timothy Ballard \email{timothy.ballard@uq.edu.au}
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Sfun}
\alias{TPM_Sfun}
\title{TPM Wake S Function}
\usage{
TPM_Sfun(la, d, sw, taw)
}
\arguments{
\item{la}{Lower asymptote (typically = 2.4)}

\item{d}{Decay in alertness (typically = -0.0353)}

\item{sw}{S upon waking}

\item{taw}{Time awake (hours)}
}
\value{
S
}
\description{
Calculates S process during wake
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIPS_df.R
\name{FIPS_df}
\alias{FIPS_df}
\title{The FIPS_df}
\usage{
FIPS_df(.data)
}
\arguments{
\item{.data}{A dataframe that matches the required FIPS_df structure}
}
\description{
All models implemented in FIPS are implemented to be run on a \code{FIPS_df} object ---
a dataframe containing a time series of all variables required to run \link{FIPS_simulate} to
generate a \code{FIPS_simulation} object (a subclass of \code{FIPS_df}).
}
\section{Specification}{

The specification for a \code{FIPS_df} is a dataframe object that contains the following
variables (columns). The only input \emph{required} to generate this is a series of sleep/wake times.
The FIPS_df is a tibble with the following variables (columns):
\itemize{
\item \code{datetime} = vector of datetime stamps separated at equidistance intervals.
\item \code{sleep.id} = a supplementary variable indicating sleep episode identifier.
\item \code{wake_status} =  Awake (\code{T}) or asleep (\code{F}) at that epoch interval/epoch
\item \code{wake_status_int} = Awake (1) or asleep (0) at that epoch interval/epoch
\item \code{change_point} = Whether the individual changed wake status at that interval/epoch.
\item \code{switch_direction} = Whether switch was to sleep or to wake
\item \code{status_duration} = How long individual has been in status at that current time point
\item \code{total_prev} = If a switch has occured, how long were that in the previous status.
\item \code{time} = time of day in decimal hour
\item \code{day} =  days into simulation
\item \code{sim_hours} = hours simulation has run for total
}

Note that is theoretically possible to generate a dataframe yourself and apply the \code{FIPS:::as_FIPS_df()} method.
Be cautious however that as of FIPS 0.1.0 there is no validator at object instantiation time. This is because the validation occurs
in internal functions specified in the \link{parse_sleeptimes}. If the API for creating your own \code{FIPS_df}
objects is opened to the user end this will be addressed.
}

\seealso{
See FIPS_simulation (internal) for additional columns added after a simulation is run.

Also see \link{parse_sleeptimes} for sleep times converter.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sleeptimesformat.R
\name{generate_postwake_times}
\alias{generate_postwake_times}
\title{Fill post-observation wake times}
\usage{
generate_postwake_times(simulationend, lastwake, expand_by = 5)
}
\arguments{
\item{simulationend}{start of simulation}

\item{lastwake}{first sleep in the sleep dataframe}

\item{expand_by}{expand value}
}
\value{
returns expanded tibble containing sleep.id = NA (due to waking) and wake_status = T
}
\description{
The last wake moment is unlikely to also be the end of the series.
This function fills constructs a tibble with the all times between
the final wake episode and the end of the series, intervaled by `expand_by` minutes
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_get_KSS_vector}
\alias{TPM_get_KSS_vector}
\title{TPM_add_KSS}
\usage{
TPM_get_KSS_vector(.FIPS_sim)
}
\arguments{
\item{.FIPS_sim}{a FIPS_simulation on the TPM}
}
\description{
Adds the KSS to a TPM FIPS_simulation. Will return error if FIPS_simulation not a TPM.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_make_pvec}
\alias{unified_make_pvec}
\title{Make Model Default (pvec) Parameters}
\usage{
unified_make_pvec(
  U0 = 24.12,
  L0 = 0,
  S0 = 0,
  phi = 2.02,
  kappa = 4.13,
  tau_s = 1,
  tau_w = 40,
  tau_la = 4.06 * 24,
  sigma = 1,
  wc = 1.14,
  wd = -0.46
)
}
\arguments{
\item{U0}{Upper asymptote (defaults to = 24.12)}

\item{L0}{Lower asymptote(defaults to = 0,   # (0.88 * 3 - 2)*1.74)}

\item{S0}{Initial starting point of S process (defaults to = 0,   # 1.11 + (1.74-1.11)*0.64)}

\item{phi}{Phase at beginning of the simulation (defaults to = 2.02)}

\item{kappa}{Relative influence of C process (defaults to = 4.13)}

\item{tau_s}{Controls rate of decay in S during sleep (defaults to = 1)}

\item{tau_w}{Controls rate of rise in S during wake (defaults to = 40)}

\item{tau_la}{Rate of change in lower asymptote (defaults to = 4.06*24)
I don't think we have any particular reason to claim sigma is Bayesian error.}

\item{sigma}{Bayesian error - ignore unless you have error calculations (defaults to = 1)}

\item{wc}{Sleep inertia: extent of alertness reduction at time of waking (typically = -5.72) (defaults to = 1.14)}

\item{wd}{Sleep inertia: exponential recovery of alertness (typically = -1.51) (defaults to = -0.4)}
}
\description{
The default unified model parameters from:
Ramakrishnan, S., Wesensten, N. J., Balkin, T. J., \& Reifman, J.
(2016). A Unified Model of Performance: Validation of its Predictions
across Different Sleep/Wake Schedules. \emph{Sleep}, \emph{39}(1),
249--262. \url{https://doi.org/10.5665/sleep.5358}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\docType{data}
\name{unified_pvec}
\alias{unified_pvec}
\title{Unified Model Default (pvec) Parameters}
\format{
An object of class \code{numeric} of length 11.
}
\usage{
unified_pvec
}
\description{
The default unified model parameters from:
}
\details{
Ramakrishnan, S., Wesensten, N. J., Balkin, T. J., \& Reifman, J.
(2016). A Unified Model of Performance: Validation of its Predictions
across Different Sleep/Wake Schedules. \emph{Sleep}, \emph{39}(1),
249--262. \url{https://doi.org/10.5665/sleep.5358}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Wfun}
\alias{TPM_Wfun}
\title{Sleep Inertia (W) Function}
\usage{
TPM_Wfun(Wc, Wd, taw)
}
\arguments{
\item{Wc}{Extent of alertness reduction at time of waking (typically = -5.72)}

\item{Wd}{Exponential recovery of alertness (typically = -1.51)}

\item{taw}{Time awake (hours)}
}
\value{
W
}
\description{
Calculates effect of sleep inertia on alertness
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_unifiedfatiguemodel.R
\name{unified_simulate}
\alias{unified_simulate}
\title{Simulate: Unified Model}
\usage{
unified_simulate(pvec, dat)
}
\arguments{
\item{pvec}{a vector of default parameters, see [unified_pvec]}

\item{dat}{input dataframe (ensure this is a FIPS_df)}
}
\value{
simulated dataset complete
}
\description{
Runs a full simulation of the 'Unified Model'.
}
\section{References}{


Rajdev, P., Thorsley, D., Rajaraman, S., Rupp, T. L., Wesensten, N. J.,
Balkin, T. J., \& Reifman, J. (2013). A unified mathematical model to
quantify performance impairment for both chronic sleep restriction and
total sleep deprivation. \emph{Journal of Theoretical Biology},
\emph{331}, 66--77. \url{https://doi.org/10.1016/j.jtbi.2013.04.013}

Ramakrishnan, S., Wesensten, N. J., Balkin, T. J., \& Reifman, J.
(2016). A Unified Model of Performance: Validation of its Predictions
across Different Sleep/Wake Schedules. \emph{Sleep}, \emph{39}(1),
249--262. \url{https://doi.org/10.5665/sleep.5358}
}

\seealso{
unified_make_pvec
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_check_pvec}
\alias{TPM_check_pvec}
\title{Check Three Process Parameter Vector}
\usage{
TPM_check_pvec(pvec)
}
\arguments{
\item{pvec}{The pvec to check contains all required three process model parameters}
}
\value{
logical
}
\description{
Checks the pvec contains required parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_simulation_dispatch}
\alias{TPM_simulation_dispatch}
\title{TPM Simulation Dispatcher}
\usage{
TPM_simulation_dispatch(dat, pvec, model_formula)
}
\arguments{
\item{dat}{input dataframe (ensure this is in FIPS format)}

\item{pvec}{a vector of default parameters}

\item{model_formula}{A formula expression object of desired model caluclation}
}
\description{
Constructor/dispatcher for TPM model simulations.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_threeprocessmodel.R
\name{TPM_Sp1fun}
\alias{TPM_Sp1fun}
\title{TPM_Sp1fun}
\usage{
TPM_Sp1fun(ha, g, bl, ss, tas)
}
\arguments{
\item{ha}{High asymptote (typically = 14.3)}

\item{g}{Rate of recovery (typically about -0.3813)}

\item{bl}{Break point (typically 12.2)}

\item{ss}{S upon falling asleep}

\item{tas}{Time asleep (hours)}
}
\value{
S
}
\description{
Calculates S during high pressure component of sleep before break point is reached
}
