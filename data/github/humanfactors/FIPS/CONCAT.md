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
