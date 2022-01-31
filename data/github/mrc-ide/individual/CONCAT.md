# Individual <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->
[![R build status](https://github.com/mrc-ide/individual/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/individual/actions)
[![codecov.io](https://codecov.io/github/mrc-ide/individual/coverage.svg)](https://codecov.io/github/mrc-ide/individual)
[![CRAN](https://www.r-pkg.org/badges/version/individual)](https://cran.r-project.org/package=individual)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03539/status.svg)](https://doi.org/10.21105/joss.03539)
<!-- badges: end -->

An R package for specifying and simulating individual-based models.

This package is designed to:

  1. encourage clear and testable components for defining your individual-based 
models, and
  2. provide memory efficient, fast code for executing your model

## Installation

The package can be installed from github using the "remotes" library

```R
library('remotes')
install_github('mrc-ide/individual')
```

Alternatively you can install individual directly from CRAN, but be aware that
the CRAN version may not be the most recent version of the package:

```R
install.packages("individual")
```

For development it is most convenient to run the code from source. You can
install the dependencies in RStudio by opening the project and selecting "Build" > "Install and Restart"

Command line users can execute:

```R
library('remotes')
install_deps('.', dependencies = TRUE)
```

Docker users can build a minimal image with

```bash
docker build . -f docker/Dockerfile -t [your image name]
```

Or if you would like devtools and documentation tools you can run

```bash
docker build . -f docker/Dockerfile.dev -t [your image name]
```

## Usage

We recommend first reading `vignette("Tutorial")` which describes
how to simulate a simple SIR model in "individual", and later `vignette("API")`
which describes in detail how to use the data structures in "individual" to
build more complicated models. If you are running into performance issues,
learn more about how to speed up your model in `vignette("Performance")`.

## Statement of need

Individual-based models are important tools for infectious disease epidemiology,
but practical use requires an implementation that is both comprehensible so that
code may be maintained and adapted, and fast. "individual" is an R package which
provides users a set of primitive classes using the [R6](https://github.com/r-lib/R6)
class system that define elements common to many tasks in infectious disease
modeling. Using R6 classes helps ensure that methods invoked on objects are
appropriate for that object type, aiding in testing and maintenance of models
programmed using "individual". Computation is carried out in C++ using 
[Rcpp](https://github.com/RcppCore/Rcpp) to link to R, helping achieve good
performance for even complex models.

"individual" provides a unique method to specify individual-based models compared
to other agent/individual-based modeling libraries, where users specify a type
for agents, which are subsequently stored in an array or other data structure.
In "individual", users instead instantiate a object for each variable which
describes some aspect of state, using the appropriate R6 class. Finding subsets
of individuals with particular combinations of state variables for further
computation can be efficiently accomplished with set operations, using a custom
bitset class implemented in C++. Additionally, the software makes no assumptions
on the types of models that may be simulated (*e.g.* mass action, network),
and updates are performed on a discrete time step.

We hope our software is useful to infectious disease modellers, ecologists, and
others who are interested in individual-based modeling in R.

## Contributing

Thank you! Please refer to the vignette on `vignette("Contributing")` for info on how to
contribute :)

## Alternatives

### Non R Software

 - [Repast](https://github.com/repast)
 - [Mesa](https://github.com/projectmesa/mesa)
 - [NetLogo](https://ccl.northwestern.edu/netlogo/)
 - [Agents.jl](https://github.com/JuliaDynamics/Agents.jl)

### Non R Software for Epi

 - [EpiFire](https://github.com/tjhladish/EpiFire)
 - [SimpactCyan](https://github.com/j0r1/simpactcyan)
 - [Numerus Model Builder](https://www.numerusinc.com)
 - NOVA
 - [EMOD](https://www.idmod.org/tools#emod)
 - [Pathogen.jl](https://github.com/jangevaare/Pathogen.jl), a package for individual-based simulation of common compartmental models.

### General R Packages

 - [nlrx](https://github.com/ropensci/nlrx), [RNetLogo](https://github.com/r-forge/rnetlogo), [NetLogoR](https://github.com/PredictiveEcology/NetLogoR) are NetLogo interfaces
 - [RRepast](https://github.com/antonio-pgarcia/RRepast) is a repast interface
 - [simecol](https://github.com/r-forge/simecol), provides classes and methods to enhance reproducibility of ecological models.

### R based DES

 - [simmeR](https://github.com/r-simmer/simmer)
 - [SpaDES](https://github.com/PredictiveEcology/SpaDES)

### R based IBMs

 - [IBMPopSim](https://github.com/DaphneGiorgi/IBMPopSim)
 - [ibm](https://github.com/roliveros-ramos/ibm)
 - [ibmcraftr](https://github.com/SaiTheinThanTun/ibmcraftr)

### R based Epi

 - [EpiModel](https://github.com/statnet/EpiModel)
 - [SimInf](https://github.com/stewid/SimInf)
 - [hybridModels](https://github.com/fernandosm/hybridModels)
 - epinet
 - [EpiDynamics](https://github.com/oswaldosantos/EpiDynamics)
 - [nosoi](https://github.com/slequime/nosoi), generate synthetic data for phylogenetic analysis
 - [EpiILMCT](https://github.com/waleedalmutiry/EpiILMCT)
 - [EpiILM](https://github.com/waleedalmutiry/EpiILM)
 - [SPARSEMODr](https://github.com/cran/SPARSEMODr)
# individual 0.1.8

  * [fix bug](https://github.com/mrc-ide/individual/pull/163) with C++ event
  listeners

# individual 0.1.7

  * Check for bad bitset input when queueing updates to CategoricalVariable objects 
  [PR here](https://github.com/mrc-ide/individual/pull/145)
  * Add "bench" package to suggests and include benchmarking scripts of major 
  functionality in `tests/performance` [PR here](https://github.com/mrc-ide/individual/pull/151)
  * Update to latest version of "testthat" package so that C++ tests of `IterableBitset`
  object can be run without giving LTO check errors (see `src/test-bitset.cpp`)
  * Add CITATION file
  * Add method `Bitset$clear` to zero out all set bits in bitsets [PR here](https://github.com/mrc-ide/individual/pull/157)
  * Fix bug ([issue here](https://github.com/mrc-ide/individual/issues/152)) where `DoubleVariable` and `IntegerVariable` updates could change size of the variable object [PR here](https://github.com/mrc-ide/individual/pull/156)
  
# individual 0.1.6

  * Added a `NEWS.md` file to track changes to the package.
  * Add Mac OS files to .gitignore
  * Update pkgdown reference organization.
  * Update [R-CMD-check workflow](https://github.com/r-lib/actions/tree/master/examples#standard-ci-workflow).
  * `Event.h` now defines class methods outside of the class definition for 
  easier readability, and add documentation.
  * `TargetedEvent$schedule` now dispatches to different C++ functions in `event.cpp`
  and `Event.h` depending on if input index is a bitset or vector (previous 
  behavior used bitset's $to_vector$ method in R to pass a vector).
  * `test-event.R` now only contains tests for `Event` class, new test file
  `test-targetedevent.R` contains a much updated suite of tests for the
  `TargetedEvent` class.
  * Fix bug where `CategoricalVariable` could be queued updates for indices in
  a vector that were outside the range of the population.
  * Update `Bitset$not` to operate in place. inplace = FALSE will be deprecated
    in 0.2.0
  * Rename the IterableBitset ~ operator to !

# individual 0.1.5

  * Added package logo.
  * Update DESCRIPTION and remove "reshape2" from suggested packages.
  * If given a `Bitset` for argument `index`, `queue_update` methods for 
  `IntegerVariable` and `DoubleVariable` pass the bitset directly to the C++ 
  functions `integer_variable_queue_update_bitset` and `double_variable_queue_update_bitset`
  rather than converting to vector and using vector methods.
  * `CategoricalVariable.h`, `IntegerVariable.h`, and `DoubleVariable.h` now define
  class methods outside of the class definition for easier readability, and add
  documentation.
  * `CategoricalVariable`, `IntegerVariable`, and `DoubleVariable` classes define
  a virtual destructor with default implementation.
  * `get_index_of_set` and `get_size_of_set_vector` methods for `IntegerVariable`
  now pass arguments by reference.
  * `get_values` method for `IntegerVariable` and `DoubleVariable` corrected to
  return value rather than reference.
  * add overload for `get_values` for `IntegerVariable` and `DoubleVariable` to
  accept `std::vector<size_t>` as argument rather than converting to bitset.
  * add function `bitset_to_vector_internal` to `IterableBitset.h`.
  * split `testthat/test/test-variables.R` into `testthat/test/test-categoricalvariable.R`,
  `testthat/test/test-integervariable.R`, and `testthat/test/test-doublevariable.R`
  * remove unnecessary `#include` statements from header files.
  * remove unnecessary comparisons for `size_t` types.
  
---
title: 'individual: An R package for individual-based epidemiological models'
tags:
  - R
  - epidemiology
  - individual based
  - agent based
  - infectious disease
  - simulation
  - stochastic
authors:
  - name: Giovanni D. Charles
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Sean L. Wu
    orcid: 0000-0002-5781-9493
    affiliation: 2
affiliations:
 - name:  MRC Centre for Global Infectious Disease Analysis, Abdul Latif Jameel Institute for Disease and Emergency Analytics (J-IDEA), Imperial College London, London, UK.
   index: 1
 - name: Division of Epidemiology and Biostatistics, School of Public Health, University of California, Berkeley, CA 94720, USA
   index: 2
date: 13 August 2017
bibliography: paper.bib
---

# Summary

`individual` is an R package which provides users a set of useful primitive elements
for specifying individual-based models (IBMs), also called agent-based models
(ABMs), with special attention to models for infectious disease epidemiology. 
Users build models by specifying variables for each characteristic describing individuals 
in the simulated population using data structures from the package.
`individual` provides efficient methods for finding
subsets of individuals based on these variables, or cohorts. Cohorts can then
be targeted for variable updates or scheduled for events.
Variable updates queued during a time step are executed at the end of a discrete time step,
and the code places no restrictions on how individuals are allowed to interact.
These data structures are designed to provide an intuitive way for users to turn their conceptual
model of a system into executable code, which is fast and memory efficient.

# Statement of need

Complex stochastic models are crucial for many tasks 
in infectious disease epidemiology.
Such models can formalize theory, generate synthetic data, evaluate counterfactual 
scenarios, forecast trends, and be used for statistical inference [@Ganyani:2021]. IBMs are a way to 
design disaggregated simulation models, usually contrasted
with mathematical models, which may model a density or concentration
of individuals, or otherwise lump individuals with similar attributes together in some
way [@Shalizi:2006]. For modeling finite numbers of individuals with significant
between-individual heterogeneity and complex dynamics, IBMs are a natural modeling
choice when a representation using mathematical models would be cumbersome
or impossible [@Willem:2017]. Even if an aggregated representation were feasible, there are many 
reasons why an individual-based representation is to be preferred. Synthetic data
may need to produce individual level outcomes, which aggregated models by their very 
nature are unable to provide [@Tracy:2018]. Other complexities, such as when events occur after
a random delay whose distribution differs from a Markovian
one, mean even aggregated models will need to store individual completion times,
necessitating more complex simulation algorithms and data structures; in such
cases it is often more straightforward to adopt an individual-based representation
from the start.

For practical use, individual-based models 
need to balance comprehensibility and speed. A fast model whose code is only
understood by the author can be difficult to use as a basis for scientific
exploration, which necessarily requires the development of multiple models to
test hypotheses or explore sensitivity to certain assumptions. On the
other hand a clear yet slow model can be practically unusable for tasks such as
uncertainty quantification or statistical inference. `individual`
provides a toolkit for users to write models that is general enough
to cover nearly all models of practical interest using simple, standardized code which is
fast enough to be useful for computationally taxing applications.

# State of the field

There are many software libraries for epidemiological simulation,
both in R and other programming languages. However, based on our review of
existing software, no other library exists in
the R language which provides users with a set of primitive elements for defining 
epidemiological models without imposing strong restrictions upon the type of model
that may be simulated (e.g.; compartmental, network, etc.), or limiting users to particular
mathematical forms for model dynamics.

### General R packages

Generic individual-based simulation packages in R include
IBMPopSim [@Giorgi:2020], ibm [@Oliveros:2016] and ibmcraftr [@Tun:2016].
IBMPopSim provides sophisticated simulation algorithms, but requires users to input C++ code
as a string which is then compiled, making it difficult to interface with the
existing R ecosystem.

### Epidemiological R packages

EpiModel [@Jenness:2018] allows the simulation of highly detailed discrete time
models on networks, relying on the statnet [@Handcock:statnet] project for
classes and algorithms. However due to its focus on directly transmitted
diseases, `individual` may be more applicable to other epidemiological situations 
such as vector-borne diseases. In addition it does not offer an interface for
compiled code.

hybridModels [@Fernando:2020], similarly provides tools for simulating epidemics
on dynamic networks. However, it is fully implemented in R, limiting the scope for
scale and optimisation.

Other packages in R are more specialised or restrict the model's transmission dynamics 
to specific mathematical forms (e.g.; mass action). These include SimInf [@Bauer:2016], 
nosoi [@Lequime:2020], SPARSEMODr [@Mihaljevic:2021], EpiILMCT [@Almutiry:2020] and
EpiILM [@Warriyar:2020].

# Design principles

Because in many epidemiological models the most important representation of state
is a finite set of mutually exclusive values, such as the Susceptible, Infectious, Recovered
classes from the well-known SIR model [@Allen:2017], `individual` uses a bitset to store these data.
At the R level users can call set operations (union, intersection,
complement, symmetric difference, set difference) which are implemented as bitwise
operations in the C++ source. This lets users write clear, highly efficient
code for updating their model, fully in R. 

In contrast to other individual-based modeling software, where users focus on
defining a type for simulated individuals,
in `individual` users instead define variables, one for each characteristic
of the simulated population.
Individual agents are defined by their their position in each bitset giving 
membership in a variable, or element in a vector of integers or floats.
This design is similar to a component system, a design pattern to help
decouple complicated types [@Nystrom:2014]. 
Because of this disaggregated representation of state, performing
operations to find and schedule cohorts of individuals benefits from fast bitwise operators.
This state representation is (to our knowledge), novel for epidemiological simulation. 
While @Rizzi:2018 proposed using a bitset to represent the state of each
simulated individual, the population was still stored as types in an array.

`individual` uses `Rcpp` [@Rcpp] to link to C++ source code, 
which underlies the data structures exposed to the user. 
The API for `individual` uses `R6` [@R6] classes at the R level
which users call to create, update, and query variables.
`individual` also provides a C++ header-only interface which advanced users
can link to from their R package.
Users can then write their own C++ code or benefit from other packages with
a compiled interface,
significantly enhancing the extensibility of `individual`'s API, and
documentation on interacting with `individual`'s C++ API is available in
the package [documentation](https://mrc-ide.github.io/individual/articles/Performance.html).

After a user has specified all the variables in their model, dynamics
are specified by processes which run each time step, and events which can be
scheduled to update specific cohorts in the future. The simulation loop then
executes processes, fires events and updates state on each discrete time step.

![A flow diagram for the simulation loop](sim_loop.png)

# Licensing and availability

`individual` is licensed under the MIT License, with all
source code stored at [GitHub](https://github.com/mrc-ide/individual).
Requests, suggestions, and bug reports are encouraged via
filing an [issue](https://github.com/mrc-ide/individual/issues).
A general guide on how to contribute to `individual` is available at
the [package's website](https://mrc-ide.github.io/individual/articles/Contributing.html).
The automated test coverage can be found at
[codecov.io](https://app.codecov.io/gh/mrc-ide/individual/). Example code can
be found in the
[tutorial section](https://mrc-ide.github.io/individual/articles/Tutorial.html)
of the package documentation.


# Acknowledgements

We would like to thank Dr. Pete Winskill and Dr. Oliver Watson for their
encouragement, testing and contributions to the repository. And Dr. Richard
Fitzjohn for his early technical feedback and expert advice on R package
development.

# References
---
title: "Contributing"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Contributing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Thank you for taking the time to contribute to Individual. We're grateful that you would take some time to make Individual-Based Modeling easier in R. Pull requests are very welcome.

## Issues

For major changes, please open an issue first to give other developers a heads up on what you would like to add or change. That way we can reduce any redundant efforts and give any resources if needed.

For bug reports please include:

 * R version
 * OS
 * Steps to recreate
 * Expected behaviour
 * Actual behaviour

## Git

We use Git on this project. Which means we use `master`, `dev`, `feat/*`, `bug/*`, `hotfix/*` and `release/*` branches. Please refer to [this post](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more information of each type of branch. NOTE: `bug/*` branches are `feat/*` branches which fix a bug.

Practically speaking, *all* new code contributions should be feature branches. You should branch off of the `dev` branch into one called `feat/[your feature name here]`. When we consider pull requests from forked repositories to the mrc-ide/individual repository, we will expect this convention.

We periodically merge `dev` into `master` for small release updates. These releases will appear on the [GitHub releases page](https://github.com/mrc-ide/individual/releases). Please use [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) as it helps us version Individual properly. 

Large releases will have a `release/*` branch to manage the transition. No new features will be merged into `release/*` branches. Only documentation and bug fixes will be considered.

## Code organisation

*R/simulation.R* - Contains the main entry point and configuration for models

*R/variables.R* - Defines classes for the available variables

*R/events.R* - Defines classes for the available events

*src/* - The C++ side of the R interface and tests

*inst/include/Variable.h* - The implementations of Variables

*inst/include/Event.h* - The implementations of Events

*tests/* - are divided into unit, integration and performance tests. Integration tests are
strongly recommended for large functions and unit tests for everything else

## Pull Requests

Here's a checklist for a successful PR:

 - Read your own diff
 - Describe your PR
 - Write any particular notes of interest for the reviewer
 - Check that your code passes all CI checks
 - Check that your code is mergeable
 
 These are the things we check for:
 
 - Do I understand the code?
 - Does the code look like it would work?
 - Does it work when run locally?
 - Is it tested enough?
 - Is it documented enough?

Our review process is based off of [RESIDE's PR review process](https://reside-ic.github.io/articles/pull-requests/)

## Microbenchmarks

We use [google benchmark](https://github.com/google/benchmark) for our
microbenchmarks. You can compile and run the benchmarks like this:

```
cd tests/performance
g++ *_benchmark.cpp -std=c++14 -lbenchmark -lpthread -o benchmark.out
./benchmark.out
```

## Wishlist

 * 90% test coverage
 * More Variables
 * Speed optimisations (tests TBC)
 * CRAN
 * Anything on the github issue board---
title: "Tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, echo=FALSE}
library(individual)
```

```{r conditional_block, eval=!pkgdown::in_pkgdown(),echo=F}
knitr::asis_output(
"## Table of Contents {#toc}

  1. [Introduction](#intro)
  2. [Specification](#spec)
  3. [Processes](#proc)
  4. [Events](#event)
  5. [Rendering](#render)
  5. [Simulation](#sim)"
)
```

## Introduction {#intro}

The [Susceptible-Infectious-Recovered (SIR) model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model) is the "hello, world!" model for infectious disease simulations, and here we describe how to build it in "individual". This tutorial will illustrate the use of events, processes, and rendering output.

### Specification {#spec}

To start, we should define some constants. The epidemic will be simulated in a population of 1000, where 5 persons are initially infectious, whose indices are randomly sampled. The effective contact rate $\beta$ will be a function of the deterministic $R_{0}$ and recovery rate $\gamma$. We also specify `dt`, which is the size of the time step $(\Delta t)$. Because individual's time steps are all of unit length, we scale transition probabilities by `dt` to create models with different sized steps, interpreting the discrete time model as a discretization of a continuous time model. If the maximum time is `tmax` then the overall number of time steps is `tmax/dt`.

```{r}
N <- 1e3
I0 <- 5
S0 <- N - I0
dt <- 0.1
tmax <- 100
steps <- tmax/dt
gamma <- 1/10
R0 <- 2.5
beta <- R0 * gamma
health_states <- c("S","I","R")
health_states_t0 <- rep("S",N)
health_states_t0[sample.int(n = N,size = I0)] <- "I"
```

Next, we will define the `individual::CategoricalVariable` which should store the model's "state".

```{r}
health <- CategoricalVariable$new(categories = health_states,initial_values = health_states_t0)
```

### Processes {#proc}

In order to model infection, we need a process. This is a function that takes only a single argument, `t`, for the current time step (unused here, but can model time-dependent processes, such as seasonality or school holiday). Within the function, we get the current number of infectious individuals, then calculate the *per-capita* force of infection on each susceptible person, $\lambda = \beta I/N$. Next we get a `individual::Bitset` containing those susceptible individuals and use the `sample` method to randomly select those who will be infected on this time step. The probability is given by $1 - e^{-\lambda\Delta t}$. This is the same as the CDF of an exponential random variate so we use `stats::pexp` to compute that quantity. Finally, we queue a state update for those individuals who were sampled.

```{r}
infection_process <- function(t){
  I <- health$get_size_of("I")
  foi <- beta * I/N
  S <- health$get_index_of("S")
  S$sample(rate = pexp(q = foi * dt))
  health$queue_update(value = "I",index = S)
}
```

### Events {#event}

Now we need to model recovery. For geometrically distributed infectious periods, we could use another process that randomly samples some individuals each time step to recover, but we'll use a `individual::TargetedEvent` to illustrate their use. The recovery event is quite simple, and the "listener" which is added is a function that is called when the event triggers, taking `target`, a `individual::Bitset` of scheduled individuals, as its second argument. Those individuals are scheduled for a state update within the listener function body. 

```{r}
recovery_event <- TargetedEvent$new(population_size = N)
recovery_event$add_listener(function(t, target) {
  health$queue_update("R", target)
})
```

Finally, we need to define a recovery process that queues future recovery events. We first get `individual::Bitset` objects of those currently infectious individuals and those who have already been scheduled for a recovery. Then, using bitwise operations, we get the intersection of already infectious persons with persons who have not been scheduled, precisely those who need to have recovery times sampled and recoveries scheduled. We sample those times from `stats::rgeom`, where the probability for recovery is $1-e^{-\gamma \Delta t}$. Note that we add one to the resulting vector, because by default R uses a "number of failures" parameterization rather than "number of trials", meaning it would be possible for an individual to be infectious for 0 time steps without the correction. Finally we schedule the recovery using the recovery event object.

We note at this point would be possible to queue the recovery event at the same time the infection state update was made, but we separate event and process for illustration of how the package works.

```{r}
recovery_process <- function(t){
  I <- health$get_index_of("I")
  already_scheduled <- recovery_event$get_scheduled()
  I$and(already_scheduled$not(inplace = TRUE))
  rec_times <- rgeom(n = I$size(),prob = pexp(q = gamma * dt)) + 1
  recovery_event$schedule(target = I,delay = rec_times)
}
```

### Rendering {#render}

The last thing to do before simulating the model is rendering output to plot. We use a `individual::Render` object which stores output from the model, which is tracked each time step. To do so requires another process, for which we use the "prefab" `individual::categorical_count_renderer_process`.

```{r}
health_render <- Render$new(timesteps = steps)
health_render_process <- categorical_count_renderer_process(
  renderer = health_render,
  variable = health,
  categories = health_states
)
```

### Simulation {#sim}

Finally, the simulation can be run by passing objects to the `individual::simulation_loop` function:

```{r}
simulation_loop(
  variables = list(health),
  events = list(recovery_event),
  processes = list(infection_process,recovery_process,health_render_process),
  timesteps = steps
)
```

We can easily plot the results by accessing the renderer.

```{r,dpi=300, fig.width=6, fig.height=4}
states <- health_render$to_dataframe()
health_cols <-  c("royalblue3","firebrick3","darkorchid3")
matplot(
  x = states[[1]]*dt, y = states[-1],
  type="l",lwd=2,lty = 1,col = adjustcolor(col = health_cols, alpha.f = 0.85),
  xlab = "Time",ylab = "Count"
)
legend(
  x = "topright",pch = rep(16,3),
  col = health_cols,bg = "transparent",
  legend = health_states, cex = 1.5
)
```---
title: "Performance"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Performance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, echo=FALSE}
library(individual)
```

```{r conditional_block, eval=!pkgdown::in_pkgdown(),echo=F}
knitr::asis_output(
"## Table of Contents {#toc}

  1. [Introduction](#intro)
  2. [Bitset](#bitset)
  3. [Prefabs](#prefab)
  4. [C++ Prefabs](#cpp_prefab)"
)
```

## Introduction {#intro}

Individual is designed for running big individual-based models. But if you find your model taking too long or consuming all of your memory, here are some things you can try to speed up your simulation.

### Bitset {#bitset}

The [bitset](https://en.wikipedia.org/wiki/Bit_array) data structure is used to record the presence or absence of an element in a finite set. `individual::Bitset` implements this data structure and is able to preform set operations extremely quickly using [bitwise operations](https://en.wikipedia.org/wiki/Bitwise_operation). Taking advantage of these operations can lead to very fast processes, when you need to find some particular subset of individuals.

Let's take a look at the recovery process in `vignette("Tutorial")`. A crucial operation here is to get all infectious individuals who are not already scheduled for recovery. The object `I` is a bitset containing all those individuals currently infectious, and `already_scheduled` is another bitset containing those individuals scheduled for a recovery event. Using `already_scheduled$not()` returns a new bitset of those individuals *not* in the set of already scheduled persons. This is passed to the `I$and()`, which modifies `I` in-place so that the result is the intersection of currently infectious persons and persons who have not yet been scheduled for a recovery event, which is precisely the set of people we want.

```{r,eval=FALSE}
recovery_process <- function(t){
  I <- health$get_index_of("I")
  already_scheduled <- recovery_event$get_scheduled()
  I$and(already_scheduled$not())
  rec_times <- rgeom(n = I$size(),prob = pexp(q = gamma * dt)) + 1
  recovery_event$schedule(target = I,delay = rec_times)
}
```

Bitsets can also be efficiently sampled using `Bitset$sample()`. This is used in the infection process of `vignette("Tutorial")`. Once the per-capita force of infection (probability of moving from S to I during this time step) is calculated, the bitset `S` is sampled with that probability which modifies it in-place. The number of elements remaining after being sampled is binomially distributed. The argument `rate` can also be specified as a vector of probabilities, one for each element in the bitset.

```{r,eval=FALSE}
infection_process <- function(t){
  I <- health$get_size_of("I")
  foi <- beta * I/N
  S <- health$get_index_of("S")
  S$sample(rate = pexp(q = foi * dt))
  health$queue_update(value = "I",index = S)
}
```

When creating a new Bitset, a user must specify the maximum size of the bitset. This is the maximum number of positive integers which the bitset can store. For example, if calling `Bitset$new(size = 100)`, the resulting object is able to store the presence or absence of integers between 1 and 100 (inclusive). Attempting to insert or remove elements outside of this range will result in an error.

Bitsets offer methods to preform unions (`Bitset$or()`), intersections (`Bitset$and()`), symmetric set difference (also known as exclusive or, `Bitset$xor()`), and set difference (`Bitset$set_difference()`) with other bitsets. These methods modify the bitset in-place. The method `Bitset$not()` gives the complement of a bitset, and returns a new `individual::Bitset` object, leaving the original bitset intact. Because these set operations use bitwise operations directly rather than more expensive relational operators, computations with bitsets are extremely fast. Taking advantage of bitset operations can help make processes in "individual" much faster.

This can be seen when implementing a common pattern in epidemiological models: sampling success or failure for a bitset of individuals, and then generating two bitsets to hold individuals sampled one way or the other. A first method might use `individual::filter_bitset`.

```{r,eval=FALSE}
n <- 1e4
bset <- Bitset$new(n)$insert(1:n)
probs <- runif(n)

keep <- probs >= 0.5

stay <- filter_bitset(bitset = bset,other = which(keep))
leave <- filter_bitset(bitset = bset,other = which(!keep))
```

This pattern is almost always slower than using the sample method with a set difference:

```{r,eval=FALSE}
stay <- bset$copy()
stay$sample(rate = probs)

leave <- bset$copy()$set_difference(stay)
```

In both instances the original bitset object `bset` is not modified. The latter pattern can be made even faster if the original may be modified by directly taking the set difference with it. For models with large population sizes, the speed differences can be substantial. 

Because a bitset stores integers in some finite set, it can be returned as an integer vector by using `Bitset$to_vector()`. However, this is a slow and expensive operation, as data must be copied into a new vector which is returned to R. If your model's dynamics require the frequent returning of integer vectors, an `individual::IntegerVariable` object will be more appropriate. However, for most discrete variables, and especially those which mirror compartments in mathematical models, bitset operations and `individual::CategoricalVariable` (which uses bitsets internally) should be preferred.

### Prefabs {#prefab}

Every time your processes ask for a variable, there is an overhead associated with moving simulation data into R, potentially incurring expensive copying of data.

Because many epidemiological models have similar state transitions, we've included several "prefab" processes and event listeners implemented in C++ which provide significant speed improvements and can be used out of the box in models. The functions return pointers which can be passed to the process list of `individual::simulate_loop` or event listeners just like closures in R. The processes available are:

  * `individual::bernoulli_process`: moves individuals from one categorical variable state to another at a constant probability
  * `individual::multi_probability_bernoulli_process`: moves individuals from one categorical variable state to another at a 
  individual level probability specified by a `individual::DoubleVariable` object
  * `individual::fixed_probability_multinomial_process`: moves individuals from one categorical variable state to a set of possible destination
  values with constant probability to leave and multinomially distributed choice of destination state.
  * `individual::multi_probability_multinomial_process`: moves individuals from one categorical variable state to a set of possible destination
  values with individual level probability to leave specified by a `individual::DoubleVariable` object and multinomially distributed choice of destination state.
  * `individual::infection_age_process`: Simulates infection for age-structured models, where individuals come into contact at a rate given by a mixing (contact) matrix.
  
Prefabs for event listeners and renderers:

  * `individual::update_category_listener`: event listener for `individual::TargetedEvent` objects which updates the categorical variable state when it fires.
  * `individual::reschedule_listener`: event listener for `individual::TargetedEvent` objects which schedules some new followup event when it fires.
  * `individual::categorical_count_renderer_process`: used for `individual::Render` objects that counts the size of each state in a categorical variable.

### C++ Prefabs {#cpp_prefab}

Unfortunately, we don't have a prefab for every situation. Please feel free to write one of your own!

These are the basic steps to add C++ processes to your R package:

1. Run `usethis::use_rcpp` to set your package up for C++ development.
2. Add `individual` to the `LinkingTo` section of your package DESCRIPTION.
3. If you package is named `mypackage`, create a header file containing `#include<individual.h>` in any of these locations:
    ```{cpp, eval=FALSE}
    src/mypackage_types.h
    src/mypackage_types.hpp
    inst/include/mypackage_types.h
    inst/include/mypackage_types.hpp
    ```
    Then this header file will be automatically included in `RcppExports.cpp`. For more information, see section "2.5 Types in Generated Code" in the [Rcpp Attributes vignette](https://CRAN.R-project.org/package=Rcpp/vignettes/Rcpp-attributes.pdf).
  
4. Create a file `src/Makecars` containing the line `CXX_STD = CXX14`. Because `individual` uses C++14 features, when compiling your package against it you must let the compiler know it should use the C++14 standard, otherwise it will not be able to compile. 
5. Write your process!

Processes in C++ are of type `process_t`, defined in `inst/include/common_types.h`. Types for listeners for `individual::Event` and `individual::TargetedEvent` are `listener_t` and `targeted_listener_t`, defined in `inst/include/Event.h`. Below is how the C++ implementation of
`multi_probability_bernoulli_process` is coded.

Note that the return type is a `Rcpp::XPtr` [(see here)](https://dirk.eddelbuettel.com/code/rcpp/html/classRcpp_1_1XPtr.html) to a `process_t`, which is implemented as a `std::function` [(see here)](http://www.cplusplus.com/reference/functional/function/) object, a C++ class that can
hold any callable type. The `Rcpp::XPtr` is initialized with a pointer to a `process_t` object, which itself holds a C++ [lambda function](https://en.cppreference.com/w/cpp/language/lambda), basically a function closure.

The lambda function captures the input arguments by value, and takes a single argument when called `t`, giving the current time step (just like process functions in R). Sampling those individuals to change state is implemented with the C++ API for these objects.

```{Rcpp,eval=FALSE}
#include <individual.h>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::XPtr<process_t> multi_probability_bernoulli_process_cpp(
    Rcpp::XPtr<CategoricalVariable> variable,
    const std::string from,
    const std::string to,
    const Rcpp::XPtr<DoubleVariable> rate_variable
){

    // make pointer to lambda function and return XPtr to R
    return Rcpp::XPtr<process_t>(
        new process_t([variable,rate_variable,from,to](size_t t){      

            // sample leavers with their unique prob
            individual_index_t leaving_individuals(variable->get_index_of(std::vector<std::string>{from}));
            std::vector<double> rate_vector = rate_variable->get_values(leaving_individuals);
            bitset_sample_multi_internal(leaving_individuals, rate_vector.begin(), rate_vector.end());

            variable->queue_update(to, leaving_individuals);

        }),
        true
    ); 
};

```

The exported function can be used normally in R when creating the list of processes:

```{r,eval=FALSE}
processes <- list(
  multi_probability_bernoulli_process_cpp(state, "I", "R", prob),
  other_process_1,
  other_process_2
)
```

That's everything you need to scale your models up to millions of individuals!
---
title: "API"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, echo=FALSE}
library(individual)
```

```{r conditional_block, eval=!pkgdown::in_pkgdown(),echo=F}
knitr::asis_output(
"## Table of Contents {#toc}

  1. [Introduction](#intro)
  2. [Variables](#var)
      1. [Categorical Variables](#cat_var)
      2. [Integer Variables](#int_var)
      3. [Double Variables](#dbl_var)
  3. [Processes](#proc)
      1. [Render](#rend)
  4. [Events](#event)
      1. [Targeted Events](#t_event)
  5. [Simulate](#sim)"
)
```
 
## Introduction {#intro}

This package defines a set of useful primitives for structuring and simulating individual based mechanistic models, with special emphasis on infectious disease modelling. *Structuring* means that the package lets you specify the state space of the model and transition rules for how state changes over time, and *simulating* means that the package defines a principled way to put state and transitions together to update state over a time step. Transitions are any functions that change state, and in "individual" come in two forms: [processes](#proc), which occur each time step and [events](#event), which are scheduled. The simulation updates over a discrete time step of unit length.

This vignette provides a description of the various primitive elements provided by the "individual" package for modelling and explains how they work together. Please see the `vignette("Tutorial")` for a practical example.

## Variables {#var}

In "individual", variables are how one defines state, and represent any attribute of an individual. While many variables will be dynamically updated throughout a simulation, they don't have to be. Specifying a baseline characteristic for each individual that doesn't change over a simulation will still be specified using a variable object.

There are 3 types of variable objects: `individual::CategoricalVariable` can represent a set of mutually exclusive categories, `individual::IntegerVariable` for representing integers, and `individual::DoubleVariable` for continuous values (real numbers).

### Categorical Variable {#cat_var}

Most epidemiological models will *require* use of categorical variables to define state, for example, the Susceptible, Infectious, and Recovered classes in the classic [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model). In general, for compartmental models, each set of compartments which describes an aspect of an individual should be mapped to a single `individual::CategoricalVariable`. For example, a model which classifies individuals as in an SIR state _and_ in a high or low risk groups could be represented with two categorical variable objects

There can be an unlimited number of categorical variable objects, but every individual can only take on a single value for each of them at any time. The allowable (categorical) state space for each individual is a contingency table where the margins are given by the values for each `individual::CategoricalVariable`.

`individual::CategoricalVariable` objects internally store state as a hash table of category values (with string keys), each of which contains a `individual::Bitset` recording the individuals in that state, making operations on `individual::CategoricalVariable` objects, or chains of operations extremely fast, and therefore should be the preferred variable type for discrete variables. However, cases in which discrete variables have a large number of values or when many values have few or no individuals in that category should use the `individual::IntegerVariable` instead.

### Integer Variable {#int_var}

An `individual::IntegerVariable` should be used for discrete variables that may either be technically unbounded, or whose set of possible values is so large that a `individual::CategoricalVariable` is impractical. `individual::IntegerVariable` objects internally stores state as a vector of integers of length equal to the number of individuals in the simulation, to be contrasted with the hash table and bitsets used in `individual::CategoricalVariable`. Examples of this variable type might include household for a model that simulated transmission of disease between a community of households, or age group for an age-structured model.

### Double Variable {#dbl_var}

A `individual::DoubleVariable` can be used for continuous variables that are represented as double precision floating-point numbers (hence, "double"). Such variables are internally stored as a vector of doubles of length equal to the number of individuals in the simulation. Examples of this variable type might include levels of immune response or concentration of some environmental pollutant in the blood.

## Processes {#proc}

In "individual", processes are functions which are called on each time step with a single argument, `t`, the current time step of the simulation. Processes are how most of the dynamics of a model are specified, and are implemented as [closures](http://adv-r.had.co.nz/Functional-programming.html#closures), which are functions that are able to access data not included in the list of function arguments.

An illustrative example can be found in the provided `individual::bernoulli_process`, which moves individuals from one value (source) of a `individual::CategoricalVariable` to another (destination) at a constant probability (such that the dwell time in the `from` state follows a [geometric distribution](https://en.wikipedia.org/wiki/Geometric_distribution)). 

```{r, eval=FALSE}
bernoulli_process <- function(variable, from, to, rate) {
  function(t) {
    variable$queue_update(
      to,
      variable$get_index_of(from)$sample(rate)
    )
  }
}
```

We can test empirically that dwell times do indeed follow a geometric distribution (and demonstrate uses of `individual::IntegerVariable` objects) by modifying slightly the function to take `int_variable`, an object of class `individual::IntegerVariable`. Each time the returned function closure is called it randomly samples a subset of individuals to move to the destination state, and records the time step at which those individuals move in the `individual::IntegerVariable`. We run the simple simulation loop until there are no more individuals left in the source state, and compare the dwell times in the source state to the results of `stats::rgeom`, and see that they are from the same distribution, using `stats::chisq.test` for comparison of discrete distributions.

```{r}
bernoulli_process_time <- function(variable, int_variable, from, to, rate) {
  function(t) {
    to_move <- variable$get_index_of(from)$sample(rate)
    int_variable$queue_update(values = t, index = to_move)
    variable$queue_update(to, to_move)
  }
}

n <- 5e4
state <- CategoricalVariable$new(c('S', 'I'), rep('S', n))
time <- IntegerVariable$new(initial_values = rep(0, n))

proc <- bernoulli_process_time(variable = state,int_variable = time,from = 'S',to = 'I',rate = 0.1)

t <- 0
while (state$get_size_of('S') > 0) {
  proc(t = t)
  state$.update()
  time$.update()
  t <- t + 1
}

times <- time$get_values()
chisq.test(times, rgeom(n,prob = 0.1),simulate.p.value = T)
```
By using `individual::IntegerVariable` and `individual::DoubleVariable` objects to modify the probability with which individuals move between states, as well as counters for recording previous state transitions and times, complex dynamics can be specified from simple building blocks. Several "prefab" functions are provided to make function closures corresponding to state transition processes common in epidemiological models.

### Render {#rend}

During each time step we often want to record (render) certain output from the simulation. Because rendering output occurs each time step, they are called as other processes during the simulation loop (this is another way we could have implemented the time counters in the previous example).

In "individual", rendering output is handled by `individual::Render` class objects, which build up the recorded data and can return a `base::data.frame` object for plotting and analysis. A `individual::Render` object is then stored in a function closure that is called each time step like other processes, and into which data may be recorded. Please note that the function closure must also store any variable objects whose values should be tracked, such as the prefab rendering process `individual::categorical_count_renderer_process`, taking a `individual::Render` object as the first argument and a `individual::CategoricalVariable` object as the second.

```{r,eval=FALSE}
categorical_count_renderer_process <- function(renderer, variable, categories) {
  function(t) {
    for (c in categories) {
      renderer$render(paste0(c, '_count'), variable$get_size_of(c), t)
    }
  }
}
```

## Events {#event}

While technically nothing more than variables and processes is required to simulate complex models in "individual", keeping track of auxiliary state, such as counters which track the amount of time until something happens (e.g., a non-geometric waiting time distribution) can quickly add unnecessary frustration and complexity. The `individual::Event` class is a way to model events that don't occur every time step, but are instead scheduled to fire at some future time after a delay. The scheduled events can be canceled, if something occurs first that should preempt them.

Events can be scheduled by processes; when an event fires, the event listener is called. The event listener is an arbitrary function that is called with the current time step `t` as its sole argument. Listeners can execute any state changes, and can themselves schedule future events (not necessarily of the same type). Listeners are implemented as function closures, just like [processes](#proc). Events can have any number of listeners attached to them which are called when the event fires.

### Targeted Events {#t_event}

Unlike `individual::Event` objects, whose listeners are called with a single argument `t`, when the listeners of `individual::TargetedEvent` objects are called, two arguments are passed. The first is the current time step `t`, and the second is a `individual::Bitset` object giving the indies of individuals which will be affected by the event. Because it stores these objects internally to schedule future updates, initialization of `individual::TargetedEvent` objects requires the total population size as an argument. When scheduling a `individual::TargetedEvent`, a user provides a vector or `individual::Bitset` of indicies for those individuals which will have the event scheduled. In addition, a user provides a delay which the event will occur after. The delay can either be a single value, so all individuals will experience the event after that delay, or a vector of values equal to the number of individuals who were scheduled. This allows, for example individual heterogeneity in distribution of waiting times for some event.

Much like regular [Events](#event), previously scheduled targeted events can be canceled, requiring either a `individual::Bitset` or vector of indices for individuals whose events should be canceled.

The code to define events is typically quite brief, although typically another process must be defined that queues them. A recovery event in an SIR model for example, might look like this, where `state` is a `individual::CategoricalVariable`:

```{r}
recovery_event <- TargetedEvent$new(population_size = 100)
recovery_event$add_listener(function(timestep, target) {
  state$queue_update('R', target)
})
```

## Simulate {#sim}

Using a combination of [variables](#var), [processes](#proc), and [events](#events), highly complex individual based models can be specified, and simulated with the function `individual::simulation_loop`, which defines the order in which state updates are preformed. 

It is possible to create conflicts between different processes, for example, if a single individual is subject to a death process and infection process, then both death and infection state changes could be scheduled for that individual during a time step. When processes or events create conflicting state updates, the *last* update will take precedence. Another way to make sure conflicts do not exist is to ensure only a single process is called for each state which samples each individual's update among all possible updates which could happen for those individuals, so that conflicting updates cannot be scheduled for individuals in a state.

```{r, echo=FALSE, fig.cap="Simulation loop flowchart", out.width = "100%", fig.align='center'}
knitr::include_graphics("sim_loop.png")
```

The equivalent code used in the simulation loop is below.

```{r, class.source="bg-primary", class.output="bg-primary",eval=FALSE}
simulation_loop <- function(
  variables = list(),
  events = list(),
  processes = list(),
  timesteps
  ) {
  if (timesteps <= 0) {
    stop('End timestep must be > 0')
  }
  for (t in seq_len(timesteps)) {
    for (process in processes) {
      execute_any_process(process, t)
    }
    for (event in events) {
      event$.process()
    }
    for (variable in variables) {
      variable$.update()
    }
    for (event in events) {
      event$.tick()
    }
  }
}
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{multi_probability_multinomial_process}
\alias{multi_probability_multinomial_process}
\title{Overdispersed multinomial process}
\usage{
multi_probability_multinomial_process(
  variable,
  source_state,
  destination_states,
  rate_variable,
  destination_probabilities
)
}
\arguments{
\item{variable}{a \code{\link{CategoricalVariable}} object.}

\item{source_state}{a string representing the source state.}

\item{destination_states}{a vector of strings representing the destination states.}

\item{rate_variable}{\code{\link{DoubleVariable}} giving individual probability of each individual in source state to leave}

\item{destination_probabilities}{probability vector of destination states.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Simulates a two-stage process where all individuals
in a given \code{source_state} sample whether to leave or not with a
individual probability specified by the \code{\link{DoubleVariable}}
object \code{rate_variable}; those who leave go to one of the \code{destination_states} with
probabilities contained in the vector \code{destination_probabilities}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{reschedule_listener}
\alias{reschedule_listener}
\title{Reschedule listener}
\usage{
reschedule_listener(event, delay)
}
\arguments{
\item{event}{a \code{\link[individual]{TargetedEvent}}.}

\item{delay}{the delay until the follow-up event.}
}
\description{
Schedules a follow-up event as the result of an event
firing.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integer_variable.R
\name{IntegerVariable}
\alias{IntegerVariable}
\title{IntegerVariable Class}
\description{
Represents a integer valued variable for an individual.
This class is similar to \code{\link[individual]{CategoricalVariable}},
but can be used for variables with unbounded ranges, or other situations where part
of an individual's state is better represented by an integer, such as
household or age bin.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{IntegerVariable$new()}}
\item \href{#method-get_values}{\code{IntegerVariable$get_values()}}
\item \href{#method-get_index_of}{\code{IntegerVariable$get_index_of()}}
\item \href{#method-get_size_of}{\code{IntegerVariable$get_size_of()}}
\item \href{#method-queue_update}{\code{IntegerVariable$queue_update()}}
\item \href{#method-.update}{\code{IntegerVariable$.update()}}
\item \href{#method-clone}{\code{IntegerVariable$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new IntegerVariable.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$new(initial_values)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{initial_values}}{a vector of the initial values for each individual}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_values"></a>}}
\if{latex}{\out{\hypertarget{method-get_values}{}}}
\subsection{Method \code{get_values()}}{
Get the variable values.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$get_values(index = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{index}}{optionally return a subset of the variable vector. If
\code{NULL}, return all values; if passed a \code{\link[individual]{Bitset}}
or integer vector, return values of those individuals.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_index_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_index_of}{}}}
\subsection{Method \code{get_index_of()}}{
Return a \code{\link[individual]{Bitset}} for individuals with some subset of values.
Either search for indices corresponding to values in \code{set}, or
for indices corresponding to values in range \eqn{[a,b]}. Either \code{set}
or \code{a} and \code{b} must be provided as arguments.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$get_index_of(set = NULL, a = NULL, b = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{set}}{a vector of values (providing \code{set} means \code{a,b} are ignored)}

\item{\code{a}}{lower bound}

\item{\code{b}}{upper bound}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_size_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_size_of}{}}}
\subsection{Method \code{get_size_of()}}{
Return the number of individuals with some subset of values.
Either search for indices corresponding to values in \code{set}, or
for indices corresponding to values in range \eqn{[a,b]}. Either \code{set}
or \code{a} and \code{b} must be provided as arguments.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$get_size_of(set = NULL, a = NULL, b = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{set}}{a vector of values (providing \code{set} means \code{a,b} are ignored)}

\item{\code{a}}{lower bound}

\item{\code{b}}{upper bound}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-queue_update"></a>}}
\if{latex}{\out{\hypertarget{method-queue_update}{}}}
\subsection{Method \code{queue_update()}}{
Queue an update for a variable. There are 4 types of variable update:

\enumerate{
 \item{Subset update: }{The argument \code{index} represents a subset of the variable to
update. The argument \code{values} should be a vector whose length matches the size of \code{index},
which represents the new values for that subset.}
 \item{Subset fill: }{The argument \code{index} represents a subset of the variable to
update. The argument \code{values} should be a single number, which fills the specified subset.}
 \item{Variable reset: }{The index vector is set to \code{NULL} and the argument \code{values}
replaces all of the current values in the simulation. \code{values} should be a vector
whose length should match the size of the population, which fills all the variable values in
the population}
 \item{Variable fill: }{The index vector is set to \code{NULL} and the argument \code{values}
should be a single number, which fills all of the variable values in 
the population.}
}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$queue_update(values, index = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{values}}{a vector or scalar of values to assign at the index}

\item{\code{index}}{is the index at which to apply the change, use \code{NULL} for the
fill options. If using indices, this may be either a vector of integers or
a \code{\link[individual]{Bitset}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.update"></a>}}
\if{latex}{\out{\hypertarget{method-.update}{}}}
\subsection{Method \code{.update()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$.update()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IntegerVariable$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{multi_probability_bernoulli_process}
\alias{multi_probability_bernoulli_process}
\title{Overdispersed Bernoulli process}
\usage{
multi_probability_bernoulli_process(variable, from, to, rate_variable)
}
\arguments{
\item{variable}{a \code{\link{CategoricalVariable}} object.}

\item{from}{a string representing the source state.}

\item{to}{a string representing the destination state.}

\item{rate_variable}{\code{\link{DoubleVariable}} giving individual probability of each individual in source state to leave.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Simulates a Bernoulli process where all individuals
in a given source state \code{from} sample whether or not 
to transition to destination state \code{to} with a
individual probability specified by the \code{\link{DoubleVariable}}
object \code{rate_variable}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{update_category_listener}
\alias{update_category_listener}
\title{Update category listener}
\usage{
update_category_listener(variable, to)
}
\arguments{
\item{variable}{a \code{\link[individual]{CategoricalVariable}} object.}

\item{to}{a string representing the destination category.}
}
\description{
Updates the category of a sub-population as the result of an
event firing, to be used in the \code{\link[individual]{TargetedEvent}}
class.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{bernoulli_process}
\alias{bernoulli_process}
\title{Bernoulli process}
\usage{
bernoulli_process(variable, from, to, rate)
}
\arguments{
\item{variable}{a categorical variable.}

\item{from}{a string representing the source category.}

\item{to}{a string representing the destination category.}

\item{rate}{the probability to move individuals between categories.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Simulate a process where individuals in a given \code{from} state
advance to the \code{to} state each time step with probability \code{rate}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{categorical_count_renderer_process}
\alias{categorical_count_renderer_process}
\title{Render Categories}
\usage{
categorical_count_renderer_process(renderer, variable, categories)
}
\arguments{
\item{renderer}{a \code{\link[individual]{Render}} object.}

\item{variable}{a \code{\link[individual]{CategoricalVariable}} object.}

\item{categories}{a character vector of categories to render.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Renders the number of individuals in each category.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event.R
\name{Event}
\alias{Event}
\title{Event Class}
\description{
Describes a general event in the simulation.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Event$new()}}
\item \href{#method-add_listener}{\code{Event$add_listener()}}
\item \href{#method-schedule}{\code{Event$schedule()}}
\item \href{#method-clear_schedule}{\code{Event$clear_schedule()}}
\item \href{#method-.tick}{\code{Event$.tick()}}
\item \href{#method-.process}{\code{Event$.process()}}
\item \href{#method-.process_listener}{\code{Event$.process_listener()}}
\item \href{#method-.process_listener_cpp}{\code{Event$.process_listener_cpp()}}
\item \href{#method-clone}{\code{Event$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialise an Event.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_listener"></a>}}
\if{latex}{\out{\hypertarget{method-add_listener}{}}}
\subsection{Method \code{add_listener()}}{
Add an event listener.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$add_listener(listener)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{listener}}{the function to be executed on the event, which takes a single
argument giving the time step when this event is triggered.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-schedule"></a>}}
\if{latex}{\out{\hypertarget{method-schedule}{}}}
\subsection{Method \code{schedule()}}{
Schedule this event to occur in the future.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$schedule(delay)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{delay}}{the number of time steps to wait before triggering the event,
can be a scalar or a vector of values for events that should be triggered
multiple times.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clear_schedule"></a>}}
\if{latex}{\out{\hypertarget{method-clear_schedule}{}}}
\subsection{Method \code{clear_schedule()}}{
Stop a future event from triggering.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$clear_schedule()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.tick"></a>}}
\if{latex}{\out{\hypertarget{method-.tick}{}}}
\subsection{Method \code{.tick()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$.tick()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.process"></a>}}
\if{latex}{\out{\hypertarget{method-.process}{}}}
\subsection{Method \code{.process()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$.process()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.process_listener"></a>}}
\if{latex}{\out{\hypertarget{method-.process_listener}{}}}
\subsection{Method \code{.process_listener()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$.process_listener(listener)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.process_listener_cpp"></a>}}
\if{latex}{\out{\hypertarget{method-.process_listener_cpp}{}}}
\subsection{Method \code{.process_listener_cpp()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$.process_listener_cpp(listener)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Event$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{infection_age_process}
\alias{infection_age_process}
\title{Infection process for age-structured models}
\usage{
infection_age_process(
  state,
  susceptible,
  exposed,
  infectious,
  age,
  age_bins,
  p,
  dt,
  mixing
)
}
\arguments{
\item{state}{a \code{\link{CategoricalVariable}} object.}

\item{susceptible}{a string representing the susceptible state (usually "S").}

\item{exposed}{a string representing the state new infections go to (usually "E" or "I").}

\item{infectious}{a string representing the infected and infectious  state (usually "I").}

\item{age}{a \code{\link{IntegerVariable}} giving the age of each individual.}

\item{age_bins}{the total number of age bins (groups).}

\item{p}{the probability of infection given a contact.}

\item{dt}{the size of the time step (in units relative to the contact rates in \code{mixing}).}

\item{mixing}{a mixing (contact) matrix between age groups.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Simulates infection for age-structured models, where
individuals contact each other at a rate given by some mixing (contact) matrix.
The force of infection on susceptibles in a given age class is computed as:
\deqn{\lambda_{i} = p \sum\limits_{j} C_{i,j} \left( \frac{I_{j}}{N_{j}} \right)  }
Where \eqn{C} is the matrix of contact rates, \eqn{p} is the probability of infection
per contact. The per-capita probability of infection for susceptible individuals is then:
\deqn{1 - e^{-\lambda_{i} \Delta t}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitset.R
\name{Bitset}
\alias{Bitset}
\title{A Bitset Class}
\description{
This is a data structure that compactly stores the presence of 
integers in some finite set (\code{max_size}), and can 
efficiently perform set operations (union, intersection, complement, symmetric
difference, set difference). 
WARNING: All operations are in-place so please use \code{$copy}
if you would like to perform an operation without destroying your current bitset.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{.bitset}}{a pointer to the underlying IterableBitset.}

\item{\code{max_size}}{the maximum size of the bitset.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Bitset$new()}}
\item \href{#method-insert}{\code{Bitset$insert()}}
\item \href{#method-remove}{\code{Bitset$remove()}}
\item \href{#method-clear}{\code{Bitset$clear()}}
\item \href{#method-size}{\code{Bitset$size()}}
\item \href{#method-or}{\code{Bitset$or()}}
\item \href{#method-and}{\code{Bitset$and()}}
\item \href{#method-not}{\code{Bitset$not()}}
\item \href{#method-xor}{\code{Bitset$xor()}}
\item \href{#method-set_difference}{\code{Bitset$set_difference()}}
\item \href{#method-sample}{\code{Bitset$sample()}}
\item \href{#method-choose}{\code{Bitset$choose()}}
\item \href{#method-copy}{\code{Bitset$copy()}}
\item \href{#method-to_vector}{\code{Bitset$to_vector()}}
\item \href{#method-clone}{\code{Bitset$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
create a bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$new(size, from = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{size}}{the size of the bitset.}

\item{\code{from}}{pointer to an existing IterableBitset to use; if \code{NULL}
make empty bitset, otherwise copy existing bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-insert"></a>}}
\if{latex}{\out{\hypertarget{method-insert}{}}}
\subsection{Method \code{insert()}}{
insert into the bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$insert(v)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v}}{an integer vector of elements to insert.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove"></a>}}
\if{latex}{\out{\hypertarget{method-remove}{}}}
\subsection{Method \code{remove()}}{
remove from the bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$remove(v)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v}}{an integer vector of elements (not indices) to remove.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clear"></a>}}
\if{latex}{\out{\hypertarget{method-clear}{}}}
\subsection{Method \code{clear()}}{
clear the bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$clear(v)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-size"></a>}}
\if{latex}{\out{\hypertarget{method-size}{}}}
\subsection{Method \code{size()}}{
get the number of elements in the set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$size()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-or"></a>}}
\if{latex}{\out{\hypertarget{method-or}{}}}
\subsection{Method \code{or()}}{
to "bitwise or" or union two bitsets.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$or(other)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{other}}{the other bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-and"></a>}}
\if{latex}{\out{\hypertarget{method-and}{}}}
\subsection{Method \code{and()}}{
to "bitwise and" or intersect two bitsets.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$and(other)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{other}}{the other bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-not"></a>}}
\if{latex}{\out{\hypertarget{method-not}{}}}
\subsection{Method \code{not()}}{
to "bitwise not" or complement a bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$not(inplace)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{inplace}}{whether to overwrite the current bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-xor"></a>}}
\if{latex}{\out{\hypertarget{method-xor}{}}}
\subsection{Method \code{xor()}}{
to "bitwise xor" or get the symmetric difference of two bitset
(keep elements in either bitset but not in their intersection).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$xor(other)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{other}}{the other bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_difference"></a>}}
\if{latex}{\out{\hypertarget{method-set_difference}{}}}
\subsection{Method \code{set_difference()}}{
Take the set difference of this bitset with another
(keep elements of this bitset which are not in \code{other}).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$set_difference(other)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{other}}{the other bitset.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-sample"></a>}}
\if{latex}{\out{\hypertarget{method-sample}{}}}
\subsection{Method \code{sample()}}{
sample a bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$sample(rate)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{rate}}{the success probability for keeping each element, can be
a single value for all elements or a vector of unique
probabilities for keeping each element.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-choose"></a>}}
\if{latex}{\out{\hypertarget{method-choose}{}}}
\subsection{Method \code{choose()}}{
choose k random items in the bitset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$choose(k)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{k}}{the number of items in the bitset to keep. The selection of
these k items from N total items in the bitset is random, and
k should be chosen such that \eqn{0 \le k \le N}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-copy"></a>}}
\if{latex}{\out{\hypertarget{method-copy}{}}}
\subsection{Method \code{copy()}}{
returns a copy the bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$copy()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_vector"></a>}}
\if{latex}{\out{\hypertarget{method-to_vector}{}}}
\subsection{Method \code{to_vector()}}{
return an integer vector of the elements
stored in this bitset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$to_vector()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Bitset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorical_variable.R
\name{CategoricalVariable}
\alias{CategoricalVariable}
\title{CategoricalVariable Class}
\description{
Represents a categorical variable for an individual.
This class should be used for discrete variables taking values in 
a finite set, such as infection, health, or behavioral state. It should
be used in preference to \code{\link[individual]{IntegerVariable}}
if possible because certain operations will be faster.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{CategoricalVariable$new()}}
\item \href{#method-get_index_of}{\code{CategoricalVariable$get_index_of()}}
\item \href{#method-get_size_of}{\code{CategoricalVariable$get_size_of()}}
\item \href{#method-get_categories}{\code{CategoricalVariable$get_categories()}}
\item \href{#method-queue_update}{\code{CategoricalVariable$queue_update()}}
\item \href{#method-.update}{\code{CategoricalVariable$.update()}}
\item \href{#method-clone}{\code{CategoricalVariable$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new CategoricalVariable
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$new(categories, initial_values)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{categories}}{a character vector of possible values}

\item{\code{initial_values}}{a character vector of the initial value for each
individual}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_index_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_index_of}{}}}
\subsection{Method \code{get_index_of()}}{
return a \code{\link[individual]{Bitset}} for individuals with the given \code{values}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$get_index_of(values)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{values}}{the values to filter}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_size_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_size_of}{}}}
\subsection{Method \code{get_size_of()}}{
return the number of individuals with the given \code{values}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$get_size_of(values)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{values}}{the values to filter}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_categories"></a>}}
\if{latex}{\out{\hypertarget{method-get_categories}{}}}
\subsection{Method \code{get_categories()}}{
return a character vector of possible values.
Note that the order of the returned vector may not be the same order
that was given when the variable was intitialized, due to the underlying
unordered storage type.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$get_categories()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-queue_update"></a>}}
\if{latex}{\out{\hypertarget{method-queue_update}{}}}
\subsection{Method \code{queue_update()}}{
queue an update for this variable
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$queue_update(value, index)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{the new value}

\item{\code{index}}{the indices of individuals whose value will be updated
to the one specified in \code{value}. This may be either a vector of integers or
a \code{\link[individual]{Bitset}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.update"></a>}}
\if{latex}{\out{\hypertarget{method-.update}{}}}
\subsection{Method \code{.update()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$.update()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CategoricalVariable$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/targeted_event.R
\name{TargetedEvent}
\alias{TargetedEvent}
\title{TargetedEvent Class}
\description{
Describes a targeted event in the simulation.
This is useful for events which are triggered for a sub-population.
}
\section{Super class}{
\code{\link[individual:Event]{individual::Event}} -> \code{TargetedEvent}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{TargetedEvent$new()}}
\item \href{#method-schedule}{\code{TargetedEvent$schedule()}}
\item \href{#method-get_scheduled}{\code{TargetedEvent$get_scheduled()}}
\item \href{#method-clear_schedule}{\code{TargetedEvent$clear_schedule()}}
\item \href{#method-.process_listener}{\code{TargetedEvent$.process_listener()}}
\item \href{#method-.process_listener_cpp}{\code{TargetedEvent$.process_listener_cpp()}}
\item \href{#method-clone}{\code{TargetedEvent$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="individual" data-topic="Event" data-id=".process">}\href{../../individual/html/Event.html#method-.process}{\code{individual::Event$.process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="individual" data-topic="Event" data-id=".tick">}\href{../../individual/html/Event.html#method-.tick}{\code{individual::Event$.tick()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="individual" data-topic="Event" data-id="add_listener">}\href{../../individual/html/Event.html#method-add_listener}{\code{individual::Event$add_listener()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialise a TargetedEvent.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$new(population_size)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{population_size}}{the size of the population.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-schedule"></a>}}
\if{latex}{\out{\hypertarget{method-schedule}{}}}
\subsection{Method \code{schedule()}}{
Schedule this event to occur in the future.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$schedule(target, delay)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target}}{the individuals to pass to the listener, this may be 
either a vector of integers or a \code{\link[individual]{Bitset}}.}

\item{\code{delay}}{the number of time steps to wait before triggering the event,
can be a scalar in which case all targeted individuals are scheduled for
for the same delay or a vector of values giving the delay for that
individual.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_scheduled"></a>}}
\if{latex}{\out{\hypertarget{method-get_scheduled}{}}}
\subsection{Method \code{get_scheduled()}}{
Get the individuals who are scheduled as a \code{\link[individual]{Bitset}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$get_scheduled()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clear_schedule"></a>}}
\if{latex}{\out{\hypertarget{method-clear_schedule}{}}}
\subsection{Method \code{clear_schedule()}}{
Stop a future event from triggering for a subset of individuals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$clear_schedule(target)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target}}{the individuals to clear, this may be either a vector of integers or
a \code{\link[individual]{Bitset}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.process_listener"></a>}}
\if{latex}{\out{\hypertarget{method-.process_listener}{}}}
\subsection{Method \code{.process_listener()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$.process_listener(listener)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.process_listener_cpp"></a>}}
\if{latex}{\out{\hypertarget{method-.process_listener_cpp}{}}}
\subsection{Method \code{.process_listener_cpp()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$.process_listener_cpp(listener)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TargetedEvent$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefab.R
\name{fixed_probability_multinomial_process}
\alias{fixed_probability_multinomial_process}
\title{Multinomial process}
\usage{
fixed_probability_multinomial_process(
  variable,
  source_state,
  destination_states,
  rate,
  destination_probabilities
)
}
\arguments{
\item{variable}{a \code{\link{CategoricalVariable}} object.}

\item{source_state}{a string representing the source state.}

\item{destination_states}{a vector of strings representing the destination states.}

\item{rate}{probability of individuals in source state to leave.}

\item{destination_probabilities}{probability vector of destination states.}
}
\value{
a function which can be passed as a process to \code{\link{simulation_loop}}.
}
\description{
Simulates a two-stage process where all individuals
in a given \code{source_state} sample whether to leave or not with probability
\code{rate}; those who leave go to one of the \code{destination_states} with
probabilities contained in the vector \code{destination_probabilities}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitset.R
\name{filter_bitset}
\alias{filter_bitset}
\title{Filter a bitset}
\usage{
filter_bitset(bitset, other)
}
\arguments{
\item{bitset}{the \code{\link{Bitset}} to filter}

\item{other}{the values to keep (may be a vector of intergers or another \code{\link{Bitset}})}
}
\description{
This non-modifying function returns a new \code{\link{Bitset}}
object of the same maximum size as the original but which only contains
those values at the indices specified by the argument \code{other}.
Indices in \code{other} may be specified either as a vector of integers or as
another bitset. Please note that filtering by another bitset is not a
"bitwise and" intersection, and will have the same behavior as providing
an equivalent vector of integer indices.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulation_loop}
\alias{simulation_loop}
\title{A premade simulation loop}
\usage{
simulation_loop(
  variables = list(),
  events = list(),
  processes = list(),
  timesteps
)
}
\arguments{
\item{variables}{a list of Variables}

\item{events}{a list of Events}

\item{processes}{a list of processes to execute on each timestep}

\item{timesteps}{the number of timesteps to simulate}
}
\description{
Run a simulation where event listeners take precedence 
over processes for state changes.
}
\examples{
population <- 4
timesteps <- 5
state <- CategoricalVariable$new(c('S', 'I', 'R'), rep('S', population))
renderer <- Render$new(timesteps)

transition <- function(from, to, rate) {
  return(function(t) {
    from_state <- state$get_index_of(from)
    state$queue_update(
      to,
      from_state$sample(rate)
    )
  })
}

processes <- list(
  transition('S', 'I', .2),
  transition('I', 'R', .1),
  transition('R', 'S', .05),
  categorical_count_renderer_process(renderer, state, c('S', 'I', 'R'))
)

simulation_loop(variables=list(state), processes=processes, timesteps=timesteps)
renderer$to_dataframe()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render.R
\name{Render}
\alias{Render}
\title{Render}
\description{
Class to render output for the simulation.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Render$new()}}
\item \href{#method-set_default}{\code{Render$set_default()}}
\item \href{#method-render}{\code{Render$render()}}
\item \href{#method-to_dataframe}{\code{Render$to_dataframe()}}
\item \href{#method-clone}{\code{Render$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialise a renderer for the simulation, creates the default state
renderers.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Render$new(timesteps)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{timesteps}}{number of timesteps in the simulation.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_default"></a>}}
\if{latex}{\out{\hypertarget{method-set_default}{}}}
\subsection{Method \code{set_default()}}{
Set a default value for a rendered output
renderers.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Render$set_default(name, value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{the variable to set a default for.}

\item{\code{value}}{the default value to set for a variable.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Update the render with new simulation data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Render$render(name, value, timestep)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{the variable to render.}

\item{\code{value}}{the value to store for the variable.}

\item{\code{timestep}}{the time-step of the data point.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_dataframe"></a>}}
\if{latex}{\out{\hypertarget{method-to_dataframe}{}}}
\subsection{Method \code{to_dataframe()}}{
Return the render as a \code{\link[base]{data.frame}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Render$to_dataframe()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Render$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/double_variable.R
\name{DoubleVariable}
\alias{DoubleVariable}
\title{DoubleVariable Class}
\description{
Represents a continuous variable for an individual.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{DoubleVariable$new()}}
\item \href{#method-get_values}{\code{DoubleVariable$get_values()}}
\item \href{#method-get_index_of}{\code{DoubleVariable$get_index_of()}}
\item \href{#method-get_size_of}{\code{DoubleVariable$get_size_of()}}
\item \href{#method-queue_update}{\code{DoubleVariable$queue_update()}}
\item \href{#method-.update}{\code{DoubleVariable$.update()}}
\item \href{#method-clone}{\code{DoubleVariable$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new DoubleVariable.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$new(initial_values)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{initial_values}}{a numeric vector of the initial value for each
individual.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_values"></a>}}
\if{latex}{\out{\hypertarget{method-get_values}{}}}
\subsection{Method \code{get_values()}}{
get the variable values.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$get_values(index = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{index}}{optionally return a subset of the variable vector. If
\code{NULL}, return all values; if passed a \code{\link[individual]{Bitset}}
or integer vector, return values of those individuals.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_index_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_index_of}{}}}
\subsection{Method \code{get_index_of()}}{
return a \code{\link[individual]{Bitset}} giving individuals 
whose value lies in an interval \eqn{[a,b]}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$get_index_of(a, b)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{a}}{lower bound}

\item{\code{b}}{upper bound}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_size_of"></a>}}
\if{latex}{\out{\hypertarget{method-get_size_of}{}}}
\subsection{Method \code{get_size_of()}}{
return the number of individuals whose value lies in an interval
Count individuals whose value lies in an interval \eqn{[a,b]}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$get_size_of(a, b)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{a}}{lower bound}

\item{\code{b}}{upper bound}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-queue_update"></a>}}
\if{latex}{\out{\hypertarget{method-queue_update}{}}}
\subsection{Method \code{queue_update()}}{
Queue an update for a variable. There are 4 types of variable update:
\enumerate{
 \item{Subset update: }{The argument \code{index} represents a subset of the variable to
update. The argument \code{values} should be a vector whose length matches the size of \code{index},
which represents the new values for that subset.}
 \item{Subset fill: }{The argument \code{index} represents a subset of the variable to
update. The argument \code{values} should be a single number, which fills the specified subset.}
 \item{Variable reset: }{The index vector is set to \code{NULL} and the argument \code{values}
replaces all of the current values in the simulation. \code{values} should be a vector
whose length should match the size of the population, which fills all the variable values in
the population}
 \item{Variable fill: }{The index vector is set to \code{NULL} and the argument \code{values}
should be a single number, which fills all of the variable values in 
the population.}
}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$queue_update(values, index = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{values}}{a vector or scalar of values to assign at the index.}

\item{\code{index}}{is the index at which to apply the change, use \code{NULL} for the
fill options. If using indices, this may be either a vector of integers or
a \code{\link[individual]{Bitset}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-.update"></a>}}
\if{latex}{\out{\hypertarget{method-.update}{}}}
\subsection{Method \code{.update()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$.update()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DoubleVariable$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
