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
