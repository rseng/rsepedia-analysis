---
title: 'mcbette: model comparison using babette'
tags:
  - R
  - phylogenetics
  - model comparison
  - nested sampling
  - babette
  - BEAST2
authors:
  - name: Richèl J.C. Bilderbeek
    orcid: 0000-0003-1107-7049
    affiliation: 1
affiliations:
 - name: Theoretical & Evolutionary Community Ecology, TRES, GELIFES, University of Groningen
   index: 1
date: 31 July 2020
bibliography: paper.bib
---

![](man/figures/mcbette_logo_50.png)

> The ``mcbette`` logo.

# Abstract

One can generate a phylogeny from a DNA alignment and a model of evolution.
Selecting an evolutionary model is non-trivial, as there are many.
``mcbette`` is an R package that determines the model that has most
evidence for having generated the alignment, from a set of models.
In this way, the model that is 'simple enough, but not simpler' can
be used to generate a phylogeny.

# Statement of need 

![](man/figures/combined.png)

> Constructing a species phylogeny (at the right) 
> from a DNA alignment (at the left)
> using an evolutionary model (the arrow). 
> ``mcbette`` allows for selecting an evolutionary model from a set of models.

``mcbette`` is an R package to do model comparison between
a set of evolutionary models on a DNA alignment, which allows
to select that model that is closest to the process consistent with
the DNA alignment and species tree.

Unlike other methods, ``mcbette`` can both be installed
and run from an R script, allowing one to run many 
analyses using different models, examine the results directly
from R and integrate ``mcbette`` into an existing R pipeline.

# Getting started

``mcbette`` is aimed at being used by anyone interested in phylogenetics
and assumes some basic knowledge about the field.
The BEAST book [@Drummond:2015] serves as an excellent starting point
about the field of phylogenetics, 
where the ``mcbette`` README and vignette show a simpled worked-out example.
The evolutionary models are those allowed by the ``babette`` 
R package [@Bilderbeek:2018], which consist of (among others) 
a site model, clock model and tree model (see 'Supported models' below for an overview). 
``babette`` is an R package to work with the phylogenetic 
tool BEAST2 [@Bouckaert:2019]. 
Additionally, ``mcbette`` uses the novel 'NS'
'BEAST2' package [@Russel:2019] to do the actual model comparison.

To see a demo of ``mcbette``, see the vignette:

```
vignette(topic = "demo", package = "mcbette")
```

# Quirks

``mcbette`` has two quirks. 
First, ``mcbette`` only works under Linux and Mac, because BEAST2 packages only 
work under Linux and Mac (that is, without using a GUI).
Second, ``mcbette`` uses the ``rJava``
package, because BEAST2 is written in Java. 
Getting ``rJava`` properly installed is the hardest part
to get ``mcbette`` working.

# Supported models

At the time of writing, these are the BEAST2 models that ``babette`` supports:

 * 1 site model: gamma site model
 * 4 nucleotide substitution models: JC (after Jukes and Cantor), 
   HKY (after Hasegawa, Kishino and Yano), TN (after Tamura and Nei), 
   generalized time-reversible model
 * 2 clock models: strict, relaxed log-normal
 * 5 tree models: birth-death, coalescent Bayesian skyline, 
   coalescent constant-population, coalescent exponential-population, Yule

To see these:

```
vignette(topic = "inference_models", package = "beautier")
```

# References


<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcbette

<!-- badges: start -->

[![peer-review](https://badges.ropensci.org/360_status.svg)](https://github.com/ropensci/software-review/issues/360)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02762/status.svg)](https://doi.org/10.21105/joss.02762)

| Branch    | [GitHub Actions](https://github.com/ropensci/mcbette/actions)                                      | [![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.org)                                                 | [![AppVeyor logo](man/figures/AppVeyor.png)](https://www.appveyor.com)                                                                                                    | [![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)                                                                                 |
| --------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `master`  | ![R-CMD-check](https://github.com/ropensci/mcbette/workflows/R-CMD-check/badge.svg?branch=master)  | [![Build Status](https://travis-ci.org/ropensci/mcbette.svg?branch=master)](https://travis-ci.org/ropensci/mcbette)  | [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/mcbette?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/mcbette)  | [![codecov.io](https://codecov.io/github/ropensci/mcbette/coverage.svg?branch=master)](https://codecov.io/github/ropensci/mcbette?branch=master)   |
| `develop` | ![R-CMD-check](https://github.com/ropensci/mcbette/workflows/R-CMD-check/badge.svg?branch=develop) | [![Build Status](https://travis-ci.org/ropensci/mcbette.svg?branch=develop)](https://travis-ci.org/ropensci/mcbette) | [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/mcbette?branch=develop&svg=true)](https://ci.appveyor.com/project/ropensci/mcbette) | [![codecov.io](https://codecov.io/github/ropensci/mcbette/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/mcbette?branch=develop) |

<!-- badges: end -->

![](man/figures/mcbette_logo.png)

`mcbette` allows for doing a Model Comparison using `babette` (hence the
name), that is, given a (say, DNA) alignment, it compares multiple
phylogenetic inference models to find the model that is likeliest to
generate that alignment. With this, one can find the phylogenetic
inference model that is simple enough, but not too simple.

To do so, `mcbette` uses `babette` \[Bilderbeek and Etienne, 2018\] with
the addition of using the BEAST2 \[Bouckaert et al., 2019\] nested
sampling package as described in \[Russell et al., 2019\].

## Installation

:warning: `mcbette` only works on Linux and Mac.

`mcbette` depends on the
[rJava](https://cran.r-project.org/package=rJava) and
[Rmpfr](https://cran.r-project.org/package=Rmpfr) packages.

On Linux, to install these, do (as root):

    apt install r-cran-rjava libmpfr-dev

After this, installing `mcbette` is easy:

``` r
install.packages("mcbette")
beastier::install_beast2()
mauricer::install_beast2_pkg("NS")
```

## Example

``` r
library(mcbette)
```

Suppose we have a DNA alignment, obtained by sequencing multiple
species:

``` r
fasta_filename <- beautier::get_beautier_path("anthus_aco_sub.fas")
ape::image.DNAbin(ape::read.FASTA(fasta_filename))
```

<img src="man/figures/README-example_alignment-1.png" width="100%" />

Note that this alignment holds too little information to really base a
publishable research on.

To create a posterior distribution of phylogenies from this alignment,
one needs to specify an inference model. An inference model is (among
others) a combination of a site model, clock model and tree model. See
the ‘Available models’ section to see the available models.

In this example we let two inference models compete.

Here is the default inference model:

``` r
inference_model_1 <- beautier::create_ns_inference_model()
message(
  paste(
    inference_model_1$site_model$name,
    inference_model_1$clock_model$name,
    inference_model_1$tree_prior$name
  )
)
#> JC69 strict yule
```

The JC69 site model assumes that the four DNA nucleotides are equally
likely to mutate from/to. The strict clock model assumes that mutation
rates of all species are equal. The Yule tree model assumes that new
species form at a constant rate for an extinction rate of zero.

The competing model has a different site, clock and tree model:

``` r
inference_model_2 <- inference_model_1
inference_model_2$site_model <- beautier::create_hky_site_model()
inference_model_2$clock_model <- beautier::create_rln_clock_model()
inference_model_2$tree_prior <- beautier::create_bd_tree_prior()
```

The HKY site model assumes that DNA substitution rates differ between
transitions (purine-to-purine or pyrimidine-to-pyrimidine) and
translations (purine-to-pyrimidine or the other way around). The relaxed
log-normal clock model assumes that mutation rates of all species are
differ, where all these rates together follow a log-normal distribution.
The birth-death tree model assumes that new species form and go extinct
at a constant rate.

`mcbette` shows the evidence (also called marginal likelihood) for each
inference model, which is the likelihood that a model has generated the
data.

``` r
if (can_run_mcbette()) {
  marg_liks <- est_marg_liks(
    fasta_filename = fasta_filename,
    inference_models = list(inference_model_1, inference_model_2)
  )
}
```

Here we display the marginal likelihoods as a table:

``` r
if (can_run_mcbette()) {
  knitr::kable(marg_liks)
}
```

The most important result are the model weights. When a model’s weight
is very close to one, one would prefer to use that inference model in
doing a Bayesian inference. If these model weights are rather similar,
one could argue to use either model.

Here we display the marginal likelihoods as a barplot:

``` r
if (can_run_mcbette()) {
  plot_marg_liks(marg_liks)
}
```

## Available models

The available site models:

``` r
beautier::get_site_model_names()
#> [1] "JC69" "HKY"  "TN93" "GTR"
```

The available clock models:

``` r
beautier::get_clock_model_names()
#> [1] "relaxed_log_normal" "strict"
```

The available tree models:

``` r
beautier::get_tree_prior_names()
#> [1] "birth_death"                    "coalescent_bayesian_skyline"   
#> [3] "coalescent_constant_population" "coalescent_exp_population"     
#> [5] "yule"
```

## Documentation

  - The mcbette vignette
  - [rOpenSci blog post: Call BEAST2 for Bayesian evolutionary analysis
    from R](https://ropensci.org/blog/2020/01/28/babette/)
  - [rOpenSci blog post: Selecting the Best Phylogenetic Evolutionary
    Model](https://ropensci.org/blog/2020/12/01/mcbette-selecting-the-best-inference-model/)

## FAQ

### Under which platforms does `mcbette` work?

`mcbette` only works on Linux and Mac, because BEAST2 package management
only works on those platforms.

### How do I let mcbette compare all models?

First, this is impossible, as there are infinitely many inference models
possible.

``` r
inference_models <- list()
i <- 1
for (site_model in beautier::create_site_models()) {
  for (clock_model in beautier::create_clock_models()) {
    for (tree_prior in beautier::create_tree_priors()) {
      inference_models[[i]] <- beautier::create_ns_inference_model(
        site_model = site_model,
        clock_model = clock_model,
        tree_prior = tree_prior
      )
      i <- i + 1
    }
  }
}
```

Now, `inference_models` holds a list of inference models, to be used
with `mcbette::est_marg_liks`.

### `Error: dir.exists(examples_folder) is not TRUE`

Currently, this line gives a suboptimal error message:

    beastier::get_beast2_example_filename("Primates.nex")

The error message is:

    Error: dir.exists(examples_folder) is not TRUE

The error message it should display is:

    Error: BEAST2 examples folder not found at path '[path]'.
    Maybe BEAST2 is not installed? 
    Tip: run 'beastier::install_beast2()'

This will be added in a future version of `beastier`.

## Using an existing BEAST2 installation

When BEAST2 is already installed, yet at a non-default location, one can
use the `beast2_bin_path` argument in `create_mcbette_beast2_options`.

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## References

  - Bilderbeek, Richel JC, and Rampal S. Etienne. “babette: BEAUti 2,
    BEAST 2 and Tracer for R.” Methods in Ecology and Evolution (2018).
    <https://doi.org/10.1111/2041-210X.13032>
  - Bouckaert R., Vaughan T.G., Barido-Sottani J., Duchêne S., Fourment
    M., Gavryushkina A., et al. (2019) BEAST 2.5: An advanced software
    platform for Bayesian evolutionary analysis. PLoS computational
    biology, 15(4), e1006650.
  - Russel, Patricio Maturana, et al. “Model selection and parameter
    inference in phylogenetics using nested sampling.” Systematic
    biology 68.2 (2019): 219-233.

# News

Newest versions at top.

## `mcbette` 1.14 (2021-05-14)

### NEW FEATURES

 * Use GitHub Actions as a continuous integration service

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.13 (2020-12-05)

### NEW FEATURES

 * Add `plot_marg_liks`

### MINOR IMPROVEMENTS

 * Add `check_marg_liks`

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.12 (2020-11-11)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Process CRAN feedback, thanks Gregor Seyer
 * Use `message` instead of `cat`/`print`
 * Fix URL in DESCRIPTION

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.11 (2020-10-16)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Depends on CRAN version of beautier (v2.4)

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.10 (2020-10-15)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Processed feedback rOpenSci
 * Processed feedback JOSS

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.9 (2020-10-10)

### NEW FEATURES

 * Transferred repository ownership to rOpenSci 

### MINOR IMPROVEMENTS

 * Processed feedback rOpenSci
 * Processed feedback JOSS

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.8.4 (2020-05-21)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Use CRAN versions of babette packages

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.8.3 (2020-02-24)

### NEW FEATURES

 * Add `mcbette_self_test` to self-test `mcbette`
 * Add `mcbette_report` to print all information needed for a bug report
 * Add `get_mcbette_state` and `set_mcbette_state` to store and restore
   the `mcbette` state. `check_mcbette_state` verifies if `mcbette` is in
   one of these three states: (1) both BEAST2 and the BEAST NS are installed,
   (2) only BEAST2 is installed (3) neither BEAST2 and the BEAST NS are 
   installed

### MINOR IMPROVEMENTS

 * Separated the package `mcbette` from a script running it.
   The latter can be found at 
   [https://github.com/ropensci/mcbette_run](https://github.com/ropensci/mcbette_run)
 * Add `README.Rmd`, build `README.md` from it
 * Add AppVeyor build, even though the core feature of `mcbette` 
   requires Linux or MacOS

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.8.2 (2020-01-15)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Added `can_run_mcbette` as a shorthand for three if-statements
 * Added `create_ns_inference_model` for convenience, as it is the same
   as `create_inference_model(mcmc = create_ns_mcmc(...), ...)`. This
   function will be moved to `beautier` in a future version
 * Vignette uses BEAST2 example alignment, instead of a simulated alignment
 * Prepare for rOpenSci review
 * Added more and better examples to the documentation

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.8.1 (2020-01-14)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Use CRAN version of `babette`

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.8 (2020-01-06)

### NEW FEATURES

 * Follow interface of beautier v2.3

### MINOR IMPROVEMENTS

 * Use CRAN versions of `beautier`, `beastier`, `tracerer`, `mauricer`
 * Processed all @lintr-bot's comments

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.7 (2019-09-10)

### NEW FEATURES

 * Higher number of particles do result in a better estimation
 * Renamed `est_marg_liks_from_models` to `est_marg_liks`,
   removed the old `est_marg_liks`

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None


## `mcbette` 1.6 (2019-08-27)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Creates BEAST2 temporary files in same folder as the BEAST2 working
   folder. This allows `mcbette` to run on the Groninger Peregrine
   computer cluster

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.5 (2019-08-18)

### NEW FEATURES

 * Show effective sample size in marginal likelihood estimation

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None


## `mcbette` 1.4 (2019-08-15)

### NEW FEATURES

 * Removed duplicate `epsilon` argument from `est_marg_liks_from_model`:
   use the `epsilon` supplied in the inference models' nested sampling
   MCMC

### MINOR IMPROVEMENTS

 * Better error message when using a CBS site model and too few taxa

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.3 (2019-08-14)

### NEW FEATURES

 * Added `est_marg_lik` to estimate a single marginal likelihood
 * Disallow failed marginal likelihood estimations in `est_marg_liks`:
   the resulting table will never contain `NA`s

### MINOR IMPROVEMENTS

 * Better error message when using a CBS site model and too few taxa
 * Builds on Bionic

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## `mcbette` 1.2

An untagged release

## `mcbette` 1.1 (2019-05-29)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Travis builds on three operating systems

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None
# Contributing

Awesome that you are reading this.

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage must be 100%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

Emails with code will not be accepted.

## Submitting bugs

Awesome. These are your options:

 * Add an Issue, with [a reprex](https://community.rstudio.com/t/faq-whats-a-reproducible-example-reprex-and-how-do-i-do-one/5219)
 * Submit a Pull Request, where the test is added to the `tests/testthat` folder.
   Fork from the `develop` branch
 * Send @richelbilderbeek an email (@richelbilderbeek will make an Issue of it)

Pull Requests should follow the same guidelines as 'Submitting code'.

## Branching policy

 * The `master` branch should always build successfully
 * The `development` branch is for developers

Hi @richelbilderbeek,

With this Pull Request I'd would like to [add reason].

Sure, I've read [CONTRIBUTING.md](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md) :+1:

Cheers, [your name]

---
name: Custom issue
about: Anything else
title: ''
labels: ''
assignees: ''

---


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
Script to reproduce the behavior:

```r
# Your R script here, without this comment :-)
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment:**
Show the results of running the following script:

```r
library(mcbette)
mcbette::mcbette_report()
```

**Additional context**
Add any other context about the problem here.
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
|site_model_name |clock_model_name |tree_prior_name | marg_log_lik| marg_log_lik_sd|    weight|
|:---------------|:----------------|:---------------|------------:|---------------:|---------:|
|JC69            |strict           |yule            |    -179.9266|        2.427167| 0.2823000|
|HKY             |strict           |yule            |    -182.2376|        1.992356| 0.0279950|
|TN93            |strict           |yule            |    -179.0459|        2.446301| 0.6810869|
|GTR             |strict           |yule            |    -183.4157|        2.636430| 0.0086181|
# `figures`

 * Bnilsson as Lady Macbeth from 
   https://commons.wikimedia.org/wiki/File:Bnilssonlmacbeth.jpg,
   which is in the public domain
 * Crown from [pngimg.com](http://pngimg.com/download/23834),
   licensed under Creative Commons 4.0 BY-NC
 * Arrow from http://www.clker.com/clipart-transparent-arrow-1.html

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# mcbette

<!-- badges: start -->
[![peer-review](https://badges.ropensci.org/360_status.svg)](https://github.com/ropensci/software-review/issues/360)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02762/status.svg)](https://doi.org/10.21105/joss.02762)

Branch   |[GitHub Actions](https://github.com/ropensci/mcbette/actions)                                     |[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.org)                                                 |[![AppVeyor logo](man/figures/AppVeyor.png)](https://www.appveyor.com)                                                                                                            |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|--------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/mcbette/workflows/R-CMD-check/badge.svg?branch=master) |[![Build Status](https://travis-ci.org/ropensci/mcbette.svg?branch=master)](https://travis-ci.org/ropensci/mcbette)  |[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/mcbette?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/mcbette)  |[![codecov.io](https://codecov.io/github/ropensci/mcbette/coverage.svg?branch=master)](https://codecov.io/github/ropensci/mcbette?branch=master)
`develop`|![R-CMD-check](https://github.com/ropensci/mcbette/workflows/R-CMD-check/badge.svg?branch=develop)|[![Build Status](https://travis-ci.org/ropensci/mcbette.svg?branch=develop)](https://travis-ci.org/ropensci/mcbette) |[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/mcbette?branch=develop&svg=true)](https://ci.appveyor.com/project/ropensci/mcbette) |[![codecov.io](https://codecov.io/github/ropensci/mcbette/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/mcbette?branch=develop)

<!-- badges: end -->


![](man/figures/mcbette_logo.png)

`mcbette` allows for doing a Model Comparison using `babette` (hence the name),
that is, given a (say, DNA) alignment, it compares multiple phylogenetic inference 
models to find the model that is likeliest to generate that alignment.
With this, one can find the phylogenetic inference model that is
simple enough, but not too simple.

To do so, `mcbette` uses `babette` [Bilderbeek and Etienne, 2018] with the addition
of using the BEAST2 [Bouckaert et al., 2019] nested sampling package
as described in [Russell et al., 2019].

## Installation

:warning: `mcbette` only works on Linux and Mac.

`mcbette` depends on the [rJava](https://cran.r-project.org/package=rJava) 
and [Rmpfr](https://cran.r-project.org/package=Rmpfr) packages.

On Linux, to install these, do (as root):

```
apt install r-cran-rjava libmpfr-dev
```

After this, installing `mcbette` is easy:

```r
install.packages("mcbette")
beastier::install_beast2()
mauricer::install_beast2_pkg("NS")
```

## Example

```{r load_mcbette}
library(mcbette)
```

Suppose we have a DNA alignment, obtained by sequencing multiple
species:

```{r example_alignment}
fasta_filename <- beautier::get_beautier_path("anthus_aco_sub.fas")
ape::image.DNAbin(ape::read.FASTA(fasta_filename))
```

Note that this alignment holds too little information to really base
a publishable research on.

To create a posterior distribution of phylogenies from this alignment,
one needs to specify an inference model. An inference model 
is (among others) a combination of a site model, clock model and tree 
model. See the 'Available models' section to see the available models.

In this example we let two inference models compete.

Here is the default inference model:

```{r example_inference_model_1}
inference_model_1 <- beautier::create_ns_inference_model()
message(
  paste(
    inference_model_1$site_model$name,
    inference_model_1$clock_model$name,
    inference_model_1$tree_prior$name
  )
)
```

The JC69 site model assumes that the four DNA nucleotides are equally
likely to mutate from/to. The strict clock model assumes that mutation
rates of all species are equal. The Yule tree model assumes that new
species form at a constant rate for an extinction rate of zero.

The competing model has a different site, clock and tree model:

```{r example_inference_model_2}
inference_model_2 <- inference_model_1
inference_model_2$site_model <- beautier::create_hky_site_model()
inference_model_2$clock_model <- beautier::create_rln_clock_model()
inference_model_2$tree_prior <- beautier::create_bd_tree_prior()
```

The HKY site model assumes that DNA substitution rates differ between
transitions (purine-to-purine or pyrimidine-to-pyrimidine) and 
translations (purine-to-pyrimidine or the other way around).
The relaxed log-normal clock model assumes that mutation
rates of all species are differ, where all these rates together follow
a log-normal distribution. 
The birth-death tree model assumes that new
species form and go extinct at a constant rate.

`mcbette` shows the evidence (also called marginal likelihood)
for each inference model, which is the likelihood that a model
has generated the data.

```{r example_est_marg_liks, cache=TRUE}
if (can_run_mcbette()) {
  marg_liks <- est_marg_liks(
    fasta_filename = fasta_filename,
    inference_models = list(inference_model_1, inference_model_2)
  )
}
```

Here we display the marginal likelihoods as a table:

```{r marg_liks_as_table}
if (can_run_mcbette()) {
  knitr::kable(marg_liks)
}
```

The most important result are the model weights.
When a model's weight is very close to one, 
one would prefer to use that inference model in doing a Bayesian inference.
If these model weights are rather similar, one could argue to use either
model. 

Here we display the marginal likelihoods as a barplot:

```{r plot_marg_liks}
if (can_run_mcbette()) {
  plot_marg_liks(marg_liks)
}
```

## Available models

The available site models:

```{r}
beautier::get_site_model_names()
```

The available clock models:

```{r}
beautier::get_clock_model_names()
```

The available tree models:

```{r}
beautier::get_tree_prior_names()
```

## Documentation

 * The mcbette vignette
 * [rOpenSci blog post: Call BEAST2 for Bayesian evolutionary analysis from R](https://ropensci.org/blog/2020/01/28/babette/)
 * [rOpenSci blog post: Selecting the Best Phylogenetic Evolutionary Model](https://ropensci.org/blog/2020/12/01/mcbette-selecting-the-best-inference-model/)

## FAQ

### Under which platforms does `mcbette` work?

`mcbette` only works on Linux and Mac, because BEAST2 package management
only works on those platforms.

### How do I let mcbette compare all models?

First, this is impossible, as there are infinitely many inference
models possible.

```{r}
inference_models <- list()
i <- 1
for (site_model in beautier::create_site_models()) {
  for (clock_model in beautier::create_clock_models()) {
    for (tree_prior in beautier::create_tree_priors()) {
      inference_models[[i]] <- beautier::create_ns_inference_model(
        site_model = site_model,
        clock_model = clock_model,
        tree_prior = tree_prior
      )
      i <- i + 1
    }
  }
}
```

Now, `inference_models` holds a list of inference models, to be used
with `mcbette::est_marg_liks`.

### `Error: dir.exists(examples_folder) is not TRUE`

Currently, this line gives a suboptimal error message:

```
beastier::get_beast2_example_filename("Primates.nex")
```

The error message is:

```
Error: dir.exists(examples_folder) is not TRUE
```

The error message it should display is:

```
Error: BEAST2 examples folder not found at path '[path]'.
Maybe BEAST2 is not installed? 
Tip: run 'beastier::install_beast2()'
```

This will be added in a future version of `beastier`.

## Using an existing BEAST2 installation

When BEAST2 is already installed, yet at a non-default location,
one can use the `beast2_bin_path` argument in
`create_mcbette_beast2_options`.

## Code of conduct

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

## References

 * Bilderbeek, Richel JC, and Rampal S. Etienne. "babette: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018). https://doi.org/10.1111/2041-210X.13032
 * Bouckaert R., Vaughan T.G., Barido-Sottani J., Duchêne S., Fourment M., Gavryushkina A., et al. (2019) BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. PLoS computational biology, 15(4), e1006650.
 * Russel, Patricio Maturana, et al. "Model selection and parameter inference in phylogenetics using nested sampling." Systematic biology 68.2 (2019): 219-233.

---
title: "Demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

![](mcbette_logo.png)

Given a DNA alignment, one wonders which phylogenetic inference model
fits that alignment best. `mcbette` ('Model Comparison using babette') 
can give the answer.

In this example, we use a 'BEAST2' example alignment and compare the fit
of two inference models on that alignment. We'll interpret the finding
in the end, concluding which inference model to use.

## Getting started

First, load `mcbette`:

```{r}
library(mcbette)
```

To use `mcbette`, BEAST2 and the BEAST2 `NS` package must be installed:

```{r}
if (rappdirs::app_dir()$os == "win") {
  message("'mcbette' can only run on Linux and MacOS")
} else if (!beastier::is_beast2_installed()) {
  message(
    "BEAST2 must be installed. ",
    "Tip: use 'beastier::install_beast2()'"
  )
} else if (!mauricer::is_beast2_ns_pkg_installed()) {
  message(
    "The BEAST2 'NS' package must be installed. ",
    "Tip: use 'mauricer::install_beast2_pkg(\"NS\")'"
  )
}
```

If you've just gotten a message that you need to install either
BEAST2 or the BEAST2 NS package, do so. The rest of this
vignette will be empty.

## Method

To run `mcbette`, we need a FASTA file with a DNA alignment in it.
We use one that is present in the `mcbette` package:

```{r}
fasta_filename <- system.file("extdata", "primates.fas", package = "mcbette")
```

Now we have the alignment saved in our FASTA file, we can display it:

```{r}
alignment <- ape::read.FASTA(fasta_filename)
image(alignment)
```

`mcbette` allows one to select the evolutionary model that has the
heighest evidence (aka marginal likelihood) for having generated that 
alignment. For more information how to set up an inference model, 
see the 'Inference models' vignette of the `beautier` package:

```
vignette("beautier", "inference_models")
```

In this example, we compare two evolutionary models on the alignment shown
above. Because the best evolutionary model will likely be used
in Bayesian inference, in this example we will use 'evolutionary model'
and 'inference model' interchangably.

One of the inference models is the default `babette` inference
model. 

```{r}
inference_model_1 <- beautier::create_ns_inference_model()
inference_model_1$site_model$name
```


only differing in their nucleotide substitution model:

  1. JC69: all nucleotides share the same mutation rate. 
     For example, the
     mutation rate from adenine to cytosine is the same as the
     mutation rate from guanine to thymine
  2. GTR: there is a different mutation rate from each nucleotide to each
     other nucleotide

```{r}
if (can_run_mcbette()) {
  # Create the two inference models
  inference_model_1 <- beautier::create_ns_inference_model(
    site_model = beautier::create_jc69_site_model()
  )
  inference_model_2 <- beautier::create_ns_inference_model(
    site_model = beautier::create_gtr_site_model()
  )
  # Shorten the run, by doing a short (dirty, unreliable) MCMC
  inference_model_1$mcmc <- beautier::create_test_ns_mcmc()
  inference_model_2$mcmc <- beautier::create_test_ns_mcmc()

  # Combine the two inference models
  inference_models <- c(list(inference_model_1), list(inference_model_2))

  # Compare the the two inference models
  marg_liks <- est_marg_liks(
    fasta_filename = fasta_filename,
    inference_models = inference_models
  )
  knitr::kable(marg_liks)
}
```

The results are interpreted by `interpret_marg_lik_estimates` as follows:

```{r}
if (can_run_mcbette()) {
  interpret_marg_lik_estimates(marg_liks)
}
```
---
title: "Demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

![](mcbette_logo.png)

Given a DNA alignment, one wonders which phylogenetic inference model
fits that alignment best. `mcbette` ('Model Comparison using babette') 
can give the answer.

In this example, we use a 'BEAST2' example alignment and compare the fit
of two inference models on that alignment. We'll interpret the finding
in the end, concluding which inference model to use.

## Getting started

First, load `mcbette`:

```{r}
library(mcbette)
```

To use `mcbette`, BEAST2 and the BEAST2 `NS` package must be installed:

```{r}
if (rappdirs::app_dir()$os == "win") {
  message("'mcbette' can only run on Linux and MacOS")
} else if (!beastier::is_beast2_installed()) {
  message(
    "BEAST2 must be installed. ",
    "Tip: use 'beastier::install_beast2()'"
  )
} else if (!mauricer::is_beast2_ns_pkg_installed()) {
  message(
    "The BEAST2 'NS' package must be installed. ",
    "Tip: use 'mauricer::install_beast2_pkg(\"NS\")'"
  )
}
```

If you've just gotten a message that you need to install either
BEAST2 or the BEAST2 NS package, do so. The rest of this
vignette will be empty.

## Method

To run `mcbette`, we need a FASTA file with a DNA alignment in it.
We use one that is present in the `mcbette` package:

```{r}
fasta_filename <- system.file("extdata", "primates.fas", package = "mcbette")
```

Now we have the alignment saved in our FASTA file, we can display it:

```{r}
alignment <- ape::read.FASTA(fasta_filename)
image(alignment)
```

`mcbette` allows one to select the evolutionary model that has the
heighest evidence (aka marginal likelihood) for having generated that 
alignment. For more information how to set up an inference model, 
see the 'Inference models' vignette of the `beautier` package:

```
vignette("beautier", "inference_models")
```

In this example, we compare two evolutionary models on the alignment shown
above. Because the best evolutionary model will likely be used
in Bayesian inference, in this example we will use 'evolutionary model'
and 'inference model' interchangably.

One of the inference models is the default `babette` inference
model. 

```{r}
inference_model_1 <- beautier::create_ns_inference_model()
inference_model_1$site_model$name
```


only differing in their nucleotide substitution model:

  1. JC69: all nucleotides share the same mutation rate. 
     For example, the
     mutation rate from adenine to cytosine is the same as the
     mutation rate from guanine to thymine
  2. GTR: there is a different mutation rate from each nucleotide to each
     other nucleotide

```{r}
if (can_run_mcbette()) {
  # Create the two inference models
  inference_model_1 <- beautier::create_ns_inference_model(
    site_model = beautier::create_jc69_site_model()
  )
  inference_model_2 <- beautier::create_ns_inference_model(
    site_model = beautier::create_gtr_site_model()
  )
  # Shorten the run, by doing a short (dirty, unreliable) MCMC
  inference_model_1$mcmc <- beautier::create_test_ns_mcmc()
  inference_model_2$mcmc <- beautier::create_test_ns_mcmc()
  
  
  # Combine the two inference models
  inference_models <- c(list(inference_model_1), list(inference_model_2))

  # Compare the the two inference models
  marg_liks <- est_marg_liks(
    fasta_filename = fasta_filename,
    inference_models = inference_models
  )
  knitr::kable(marg_liks)
}
```

The results are interpreted by `interpret_marg_lik_estimates` as follows:

```{r}
if (can_run_mcbette()) {
  interpret_marg_lik_estimates(marg_liks)
}
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_marg_liks.R
\name{is_marg_liks}
\alias{is_marg_liks}
\title{Determine if the \code{marg_liks} is valid}
\usage{
is_marg_liks(marg_liks, verbose = FALSE)
}
\arguments{
\item{marg_liks}{a table of (estimated) marginal likelihoods,
as, for example, created by \link{est_marg_liks}.
This \link{data.frame} has the following columns:
\itemize{
  \item \code{site_model_name}: name of the site model,
    must be an element of \link[beautier]{get_site_model_names}
  \item \code{clock_model_name}: name of the clock model,
    must be an element of \link[beautier]{get_clock_model_names}
  \item \code{tree_prior_name}: name of the tree prior,
    must be an element of \link[beautier]{get_tree_prior_names}
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
Use \link{get_test_marg_liks} to get a test \code{marg_liks}.
Use \link{is_marg_liks} to determine if a \code{marg_liks} is valid.
Use \link{check_marg_liks} to check that a \code{marg_liks} is valid.}

\item{verbose}{if TRUE show debug output}
}
\value{
TRUE if the argument is a valid \code{marg_liks},
  FALSE otherwise
}
\description{
Determine if the \code{marg_liks} is valid
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_marg_liks.R
\name{plot_marg_liks}
\alias{plot_marg_liks}
\title{Plot the \code{marg_liks}}
\usage{
plot_marg_liks(marg_liks)
}
\arguments{
\item{marg_liks}{a table of (estimated) marginal likelihoods,
as, for example, created by \link{est_marg_liks}.
This \link{data.frame} has the following columns:
\itemize{
  \item \code{site_model_name}: name of the site model,
    must be an element of \link[beautier]{get_site_model_names}
  \item \code{clock_model_name}: name of the clock model,
    must be an element of \link[beautier]{get_clock_model_names}
  \item \code{tree_prior_name}: name of the tree prior,
    must be an element of \link[beautier]{get_tree_prior_names}
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
Use \link{get_test_marg_liks} to get a test \code{marg_liks}.
Use \link{is_marg_liks} to determine if a \code{marg_liks} is valid.
Use \link{check_marg_liks} to check that a \code{marg_liks} is valid.}
}
\value{
a \link[ggplot2]{ggplot}
}
\description{
Plot the \code{marg_liks}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcbette.R
\docType{package}
\name{mcbette}
\alias{mcbette}
\title{mcbette: Model Comparison Using Babette}
\description{
'mcbette' does a model comparing using \link[babette]{babette},
where the models are Bayesian phylogenetic models,
as created by \link[beautier]{create_inference_model}.
}
\details{
The main function is \link{est_marg_liks},
which estimate the marginal likelihoods (aka evidence)
for one or more inference models, based on a single alignment.
Also, the marginal likelihoods are compared, resulting in a
relative weight for each model, where a relative weight of a model
close to \code{1.0} means that that model is way likelier than
the others.

In the process, multiple (temporary) files are created (where
\code{[x]} denotes the index in a list)

\itemize{
  \item \code{beast2_optionses[x]$input_filename}
    path to the the BEAST2 XML input file
  \item \code{beast2_optionses[x]$output_state_filename}
    path to the BEAST2 XML state file
  \item \code{inference_models[x]$mcmc$tracelog$filename}
    path to the BEAST2 trace file with parameter estimates
  \item \code{inference_models[x]$mcmc$treelog$filename}
    path to the BEAST2 \code{trees} file with the posterior trees
  \item \code{inference_models[x]$mcmc$screenlog$filename}
    path to the BEAST2 screen output file
}

These file can be deleted manually by \link[babette]{bbt_delete_temp_files},
else these will be deleted automatically by the operating system.
}
\examples{
if (can_run_mcbette()) {

  # An example FASTA file
  fasta_filename <- system.file("extdata", "simple.fas", package = "mcbette")

  inference_model_1 <- beautier::create_ns_inference_model(
    site_model = beautier::create_jc69_site_model()
  )
  inference_model_2 <- beautier::create_ns_inference_model(
    site_model = beautier::create_gtr_site_model()
  )

  # Shorten the run, by doing a short (dirty, unreliable) MCMC
  inference_model_1$mcmc <- beautier::create_test_ns_mcmc()
  inference_model_2$mcmc <- beautier::create_test_ns_mcmc()

  inference_models <- c(list(inference_model_1), list(inference_model_2))

  # Estimate the marginal log-likelihoods of the two models
  marg_liks <- est_marg_liks(
    fasta_filename = fasta_filename,
    inference_models = inference_models
  )

  # Interpret the results
  interpret_marg_lik_estimates(marg_liks)
}
}
\seealso{
Use \link{can_run_mcbette} to see if 'mcbette' can run.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_params_doc.R
\name{default_params_doc}
\alias{default_params_doc}
\title{Documentation of general function arguments.
This function does nothing.
It is intended to inherit function argument documentation.}
\usage{
default_params_doc(
  beast2_bin_path,
  beast2_folder,
  beast2_working_dir,
  beast2_options,
  beast2_optionses,
  clock_model,
  clock_models,
  epsilon,
  fasta_filename,
  inference_model,
  inference_models,
  marg_liks,
  mcbette_state,
  mcmc,
  os,
  rng_seed,
  site_model,
  site_models,
  tree_prior,
  tree_priors,
  verbose
)
}
\arguments{
\item{beast2_bin_path}{path to the the BEAST2 binary file}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{beast2_working_dir}{folder in which BEAST2 will run and
produce intermediate files.
By default, this is a temporary folder}

\item{beast2_options}{a \code{beast2_options} structure,
as can be created by \link[beastier]{create_mcbette_beast2_options}.}

\item{beast2_optionses}{list of one or more \code{beast2_options}
structures,
as can be created by \link[beastier]{create_mcbette_beast2_options}.
Use of reduplicated plural to achieve difference with
\code{beast2_options}}

\item{clock_model}{a clock model,
as can be created by \link[beautier]{create_clock_model}}

\item{clock_models}{a list of one or more clock models,
as can be created by \link[beautier]{create_clock_models}}

\item{epsilon}{measure of relative accuracy.
Smaller values result in longer, more precise estimations}

\item{fasta_filename}{name of the FASTA file}

\item{inference_model}{an inference model,
as can be created by \link[beautier]{create_inference_model}}

\item{inference_models}{a list of one or more inference models,
as can be created by \link[beautier]{create_inference_model}}

\item{marg_liks}{a table of (estimated) marginal likelihoods,
as, for example, created by \link{est_marg_liks}.
This \link{data.frame} has the following columns:
\itemize{
  \item \code{site_model_name}: name of the site model,
    must be an element of \link[beautier]{get_site_model_names}
  \item \code{clock_model_name}: name of the clock model,
    must be an element of \link[beautier]{get_clock_model_names}
  \item \code{tree_prior_name}: name of the tree prior,
    must be an element of \link[beautier]{get_tree_prior_names}
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
Use \link{get_test_marg_liks} to get a test \code{marg_liks}.
Use \link{is_marg_liks} to determine if a \code{marg_liks} is valid.
Use \link{check_marg_liks} to check that a \code{marg_liks} is valid.}

\item{mcbette_state}{the \link{mcbette} state,
which is a \link{list} with the following elements:
\itemize{
  \item{
    beast2_installed
      \link{TRUE} if BEAST2 is installed,
      \link{FALSE} otherwise
  }
  \item{
     ns_installed
      \link{NA} if BEAST2 is not installed.
      \link{TRUE} if the BEAST2 NS package is installed
      \link{FALSE} if the BEAST2 NS package is not installed
  }
}}

\item{mcmc}{an MCMC for the Nested Sampling run,
as can be created by \link[beautier]{create_mcmc_nested_sampling}}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}

\item{rng_seed}{a random number generator seed used for the BEAST2
inference}

\item{site_model}{a site model,
as can be created by \link[beautier]{create_site_model}}

\item{site_models}{a list of one or more site models,
as can be created by \link[beautier]{create_site_models}}

\item{tree_prior}{a tree prior,
as can be created by \link[beautier]{create_tree_prior}}

\item{tree_priors}{a list of one or more tree priors,
as can be created by \link[beautier]{create_tree_priors}}

\item{verbose}{if TRUE show debug output}
}
\description{
Documentation of general function arguments.
This function does nothing.
It is intended to inherit function argument documentation.
}
\note{
This is an internal function, so it should be marked with
  \code{@noRd}. This is not done, as this will disallow all
  functions to find the documentation parameters
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_marg_liks.R
\name{check_marg_liks}
\alias{check_marg_liks}
\title{Check if the \code{marg_liks} are of the same type as returned
by \link{est_marg_liks}.}
\usage{
check_marg_liks(marg_liks)
}
\arguments{
\item{marg_liks}{a table of (estimated) marginal likelihoods,
as, for example, created by \link{est_marg_liks}.
This \link{data.frame} has the following columns:
\itemize{
  \item \code{site_model_name}: name of the site model,
    must be an element of \link[beautier]{get_site_model_names}
  \item \code{clock_model_name}: name of the clock model,
    must be an element of \link[beautier]{get_clock_model_names}
  \item \code{tree_prior_name}: name of the tree prior,
    must be an element of \link[beautier]{get_tree_prior_names}
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
Use \link{get_test_marg_liks} to get a test \code{marg_liks}.
Use \link{is_marg_liks} to determine if a \code{marg_liks} is valid.
Use \link{check_marg_liks} to check that a \code{marg_liks} is valid.}
}
\description{
\link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_test_marg_liks.R
\name{get_test_marg_liks}
\alias{get_test_marg_liks}
\title{Get testing \code{marg_liks}}
\usage{
get_test_marg_liks()
}
\description{
Get testing \code{marg_liks}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_weights.R
\name{calc_weights}
\alias{calc_weights}
\title{Calculate the weights for each marginal likelihood}
\usage{
calc_weights(marg_liks)
}
\arguments{
\item{marg_liks}{(non-log) marginal likelihood estimates}
}
\value{
the weight of each marginal likelihood estimate,
which will sum up to 1.0
}
\description{
Calculate the weights for each marginal likelihood
}
\examples{
# Evidences (aka marginal likelihoods) can be very small
evidences <- c(0.0001, 0.0002, 0.0003, 0.0004)

# Sum will be 1.0
calc_weights(evidences)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_ns_pkg.R
\name{check_beast2_ns_pkg}
\alias{check_beast2_ns_pkg}
\title{Checks if the BEAST2 'NS' package is installed.}
\usage{
check_beast2_ns_pkg(beast2_bin_path = beastier::get_default_beast2_bin_path())
}
\arguments{
\item{beast2_bin_path}{path to the the BEAST2 binary file}
}
\description{
Checks if the BEAST2 'NS' package is installed.
Will \link{stop} if not
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/can_run_mcbette.R
\name{can_run_mcbette}
\alias{can_run_mcbette}
\title{Can 'mcbette' run?}
\usage{
can_run_mcbette(beast2_folder = beastier::get_default_beast2_folder())
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}
}
\description{
Can 'mcbette' run?
Will return \link{TRUE} if:
\itemize{
  \item (1) Running on Linux or MacOS
  \item (2) BEAST2 is installed
  \item (3) The BEAST2 NS package is installed
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcbette_self_test.R
\name{mcbette_self_test}
\alias{mcbette_self_test}
\title{Performs a minimal \link{mcbette} run}
\usage{
mcbette_self_test(beast2_folder = beastier::get_default_beast2_folder())
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}
}
\description{
Performs a minimal \link{mcbette} run
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mcbette_state.R
\name{get_mcbette_state}
\alias{get_mcbette_state}
\title{Get the current state of \link{mcbette}}
\usage{
get_mcbette_state(beast2_folder = beastier::get_default_beast2_folder())
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}
}
\value{
a \link{list} with the following elements:
\itemize{
  \item{
    beast2_installed
      \link{TRUE} if BEAST2 is installed,
      \link{FALSE} otherwise
  }
  \item{
     ns_installed
      \link{TRUE} if the BEAST2 NS package is installed
      \link{FALSE} if the BEAST2 or the BEAST2 NS package is not installed
  }
}
}
\description{
Get the current state of \link{mcbette}
}
\examples{
get_mcbette_state()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_mcbette_state.R
\name{set_mcbette_state}
\alias{set_mcbette_state}
\title{Set the \link{mcbette} state.}
\usage{
set_mcbette_state(
  mcbette_state,
  beast2_folder = beastier::get_default_beast2_folder(),
  verbose = FALSE
)
}
\arguments{
\item{mcbette_state}{the \link{mcbette} state,
which is a \link{list} with the following elements:
\itemize{
  \item{
    beast2_installed
      \link{TRUE} if BEAST2 is installed,
      \link{FALSE} otherwise
  }
  \item{
     ns_installed
      \link{NA} if BEAST2 is not installed.
      \link{TRUE} if the BEAST2 NS package is installed
      \link{FALSE} if the BEAST2 NS package is not installed
  }
}}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{verbose}{if TRUE show debug output}
}
\description{
Set the \link{mcbette} state to having BEAST2 installed with
or without installing the BEAST2 NS package.
}
\note{
In newer versions of BEAST2, BEAST2 comes pre-installed with the
BEAST2 NS package. For such a version, one cannot install BEAST2
without NS. A warning will be issues if one intends to only install
BEAST2 (i.e. without the BEAST2 NS package) and gets the BEAST2
NS package installed as a side effect as well.

Also, installing or uninstalling a BEAST2 package from a BEAST2
installation will affect all installations.
}
\examples{
mcbette_state <- get_mcbette_state()
mcbette_state$beast2_installed <- TRUE
mcbette_state$ns_installed <- TRUE
\donttest{
  set_mcbette_state(mcbette_state)
}
}
\seealso{
\itemize{
  \item Use \link{get_mcbette_state} to
    get the current \link{mcbette} state
  \item Use \link{check_mcbette_state} to
    check the current \link{mcbette} state
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcbette_report.R
\name{mcbette_report}
\alias{mcbette_report}
\title{Create a \link{mcbette} report,
to be used when reporting bugs}
\usage{
mcbette_report(beast2_folder = beastier::get_default_beast2_folder())
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}
}
\value{
nothing. It is intended that the output (not
the return value) is copy-pasted from screen.
}
\description{
Create a \link{mcbette} report,
to be used when reporting bugs
}
\examples{
mcbette_report()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_marg_lik_estimates.R
\name{interpret_marg_lik_estimates}
\alias{interpret_marg_lik_estimates}
\title{Interpret the marginal likelihood estimates}
\usage{
interpret_marg_lik_estimates(marg_liks)
}
\arguments{
\item{marg_liks}{a table of (estimated) marginal likelihoods,
as, for example, created by \link{est_marg_liks}.
This \link{data.frame} has the following columns:
\itemize{
  \item \code{site_model_name}: name of the site model,
    must be an element of \link[beautier]{get_site_model_names}
  \item \code{clock_model_name}: name of the clock model,
    must be an element of \link[beautier]{get_clock_model_names}
  \item \code{tree_prior_name}: name of the tree prior,
    must be an element of \link[beautier]{get_tree_prior_names}
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
Use \link{get_test_marg_liks} to get a test \code{marg_liks}.
Use \link{is_marg_liks} to determine if a \code{marg_liks} is valid.
Use \link{check_marg_liks} to check that a \code{marg_liks} is valid.}
}
\description{
Interpret the marginal likelihood estimates
as created by \link{est_marg_liks}.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_marg_liks.R
\name{est_marg_liks}
\alias{est_marg_liks}
\title{Estimate the marginal likelihoods for one or more inference models}
\usage{
est_marg_liks(
  fasta_filename,
  inference_models = list(beautier::create_inference_model(mcmc =
    beautier::create_ns_mcmc())),
  beast2_optionses = rep(list(beastier::create_mcbette_beast2_options()), times =
    length(inference_models)),
  verbose = FALSE,
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{fasta_filename}{name of the FASTA file}

\item{inference_models}{a list of one or more inference models,
as can be created by \link[beautier]{create_inference_model}}

\item{beast2_optionses}{list of one or more \code{beast2_options}
structures,
as can be created by \link[beastier]{create_mcbette_beast2_options}.
Use of reduplicated plural to achieve difference with
\code{beast2_options}}

\item{verbose}{if TRUE show debug output}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
a \link{data.frame} showing the estimated marginal likelihoods
(and its estimated error) per combination of models. Columns are:
\itemize{
  \item \code{site_model_name}: name of the site model
  \item \code{clock_model_name}: name of the clock model
  \item \code{tree_prior_name}: name of the tree prior
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{weight}: relative model weight, a value from 1.0 (all
    evidence is in favor of this model combination) to 0.0 (no
    evidence in favor of this model combination)
  \item \code{ess}: effective sample size of the marginal likelihood
    estimation
}
}
\description{
Estimate the marginal likelihoods (aka evidence)
for one or more inference models, based on a single alignment.
Also, the marginal likelihoods are compared, resulting in a
relative weight for each model, where a relative weight of a model
close to \code{1.0} means that that model is way likelier than
the others.
}
\details{
In the process, multiple (temporary) files are created (where
\code{[x]} denotes the index in a list)

\itemize{
  \item \code{beast2_optionses[x]$input_filename}
    path to the the BEAST2 XML input file
  \item \code{beast2_optionses[x]$output_state_filename}
    path to the BEAST2 XML state file
  \item \code{inference_models[x]$mcmc$tracelog$filename}
    path to the BEAST2 trace file with parameter estimates
  \item \code{inference_models[x]$mcmc$treelog$filename}
    path to the BEAST2 \code{trees} file with the posterior trees
  \item \code{inference_models[x]$mcmc$screenlog$filename}
    path to the BEAST2 screen output file
}

These file can be deleted manually by \link[babette]{bbt_delete_temp_files},
else these will be deleted automatically by the operating system.
}
\examples{
if (can_run_mcbette()) {

  # Use an example FASTA file
  fasta_filename <- system.file("extdata", "simple.fas", package = "mcbette")

  # Create two inference models
  inference_model_1 <- beautier::create_ns_inference_model(
    site_model = beautier::create_jc69_site_model()
  )
  inference_model_2 <- beautier::create_ns_inference_model(
    site_model = beautier::create_hky_site_model()
  )

  # Shorten the run, by doing a short (dirty, unreliable) MCMC
  inference_model_1$mcmc <- beautier::create_test_ns_mcmc()
  inference_model_2$mcmc <- beautier::create_test_ns_mcmc()

  # Combine the inference models
  inference_models <- list(inference_model_1, inference_model_2)

  # Create the BEAST2 options, that will write the output
  # to different (temporary) filanems
  beast2_options_1 <- beastier::create_mcbette_beast2_options()
  beast2_options_2 <- beastier::create_mcbette_beast2_options()

  # Combine the two BEAST2 options sets,
  # use reduplicated plural
  beast2_optionses <- list(beast2_options_1, beast2_options_2)

  # Compare the models
  marg_liks <- est_marg_liks(
    fasta_filename,
    inference_models = inference_models,
    beast2_optionses = beast2_optionses
  )

  # Interpret the results
  interpret_marg_lik_estimates(marg_liks)
}
}
\seealso{
\itemize{
  \item \link{can_run_mcbette}: see if 'mcbette' can run
  \item \link{est_marg_liks}: estimate multiple marginal likelihood of a
    single inference mode
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpret_bayes_factor.R
\name{interpret_bayes_factor}
\alias{interpret_bayes_factor}
\title{Interpret a Bayes factor}
\usage{
interpret_bayes_factor(bayes_factor)
}
\arguments{
\item{bayes_factor}{Bayes factor to be interpreted}
}
\value{
a string with the interpretation in English
}
\description{
Interpret a Bayes factor, using the interpretation from [1].
}
\details{
\itemize{
  \item [1] H. Jeffreys (1961). The Theory of Probability (3rd ed.).
    Oxford. p. 432
}
}
\examples{
interpret_bayes_factor(0.5)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_marg_lik.R
\name{est_marg_lik}
\alias{est_marg_lik}
\title{Estimate the marginal likelihood for an inference model.}
\usage{
est_marg_lik(
  fasta_filename,
  inference_model = beautier::create_ns_inference_model(),
  beast2_options = beastier::create_mcbette_beast2_options(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{fasta_filename}{name of the FASTA file}

\item{inference_model}{an inference model,
as can be created by \link[beautier]{create_inference_model}}

\item{beast2_options}{a \code{beast2_options} structure,
as can be created by \link[beastier]{create_mcbette_beast2_options}.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
a \link{list} showing the estimated marginal likelihoods
(and its estimated error), its items are::
\itemize{
  \item \code{marg_log_lik}: estimated marginal (natural) log likelihood
  \item \code{marg_log_lik_sd}: estimated error of \code{marg_log_lik}
  \item \code{esses} the Effective Sample Size
}
}
\description{
Estimate the marginal likelihood for an inference model.
}
\examples{
if (can_run_mcbette()) {

  # An example FASTA file
  fasta_filename <- system.file("extdata", "simple.fas", package = "mcbette")

  # A testing inference model with inaccurate (thus fast) marginal
  # likelihood estimation
  inference_model <- beautier::create_ns_inference_model()

  # Shorten the run, by doing a short (dirty, unreliable) MCMC
  inference_model$mcmc <- beautier::create_test_ns_mcmc()

  # Setup the options for BEAST2 to be able to call BEAST2 packages
  beast2_options <- beastier::create_mcbette_beast2_options()

  # Estimate the marginal likelihood
  est_marg_lik(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\seealso{
\itemize{
  \item \link{can_run_mcbette}: see if 'mcbette' can run
  \item \link{est_marg_liks}: estimate multiple marginal likelihoods
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mcbette_state.R
\name{check_mcbette_state}
\alias{check_mcbette_state}
\title{Check if the \code{mcbette_state} is valid.}
\usage{
check_mcbette_state(mcbette_state)
}
\arguments{
\item{mcbette_state}{the \link{mcbette} state,
which is a \link{list} with the following elements:
\itemize{
  \item{
    beast2_installed
      \link{TRUE} if BEAST2 is installed,
      \link{FALSE} otherwise
  }
  \item{
     ns_installed
      \link{NA} if BEAST2 is not installed.
      \link{TRUE} if the BEAST2 NS package is installed
      \link{FALSE} if the BEAST2 NS package is not installed
  }
}}
}
\description{
Check if the \code{mcbette_state} is valid.
Will \link{stop} otherwise.
}
\author{
Richèl J.C. Bilderbeek
}
