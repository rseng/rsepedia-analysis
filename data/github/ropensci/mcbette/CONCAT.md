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

