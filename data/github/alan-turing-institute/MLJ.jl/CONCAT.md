# Road map

February 2020; updated, May 2021

Please visit [contributing guidelines](CONTRIBUTING.md) if interested
in contributing to MLJ.

### Guiding goals

-   **Usability, interoperability, extensibility, reproducibility,**
	and **code transparency**.

-   Offer state-of-art tools for model **composition** and model
	**optimization** (hyper-parameter tuning)

-   Avoid common **pain-points** of other frameworks:

	-   identifying all models that solve a given task

	-   routine operations requiring a lot of code

	-   passage from data source to algorithm-specific data format

	-   probabilistic predictions: inconsistent representations, lack
		of options for performance evaluation

-   Add some focus to julia machine learning software development more
	generally

### Priorities

Priorities are somewhat fluid, depending on funding offers and
available talent. Rough priorities for the core development team at
present are marked with **†** below. However, we are always keen to
review external contributions in any area.

## Future enhancements

The following road map is more big-picture; see also [this
list](https://github.com/alan-turing-institute/MLJ.jl/issues/673).


### Adding models

- [ ] **Integrate deep learning** using [Flux.jl](https://github.com/FluxML/Flux.jl.git) deep learning.  [Done](https://github.com/FluxML/MLJFlux.jl) but can
  improve experience by:

  - [x] finishing iterative model wrapper [#139](https://github.com/alan-turing-institute/MLJ.jl/issues/139)

  - [ ] improving performance by implementing data front-end after (see [MLJBase
  #501](https://github.com/JuliaAI/MLJBase.jl/pull/501)) but see also [this relevant discussion](https://github.com/FluxML/MLJFlux.jl/issues/97).


-  [ ] Probabilistic programming:
   [Turing.jl](https://github.com/TuringLang/Turing.jl),
   [Gen](https://github.com/probcomp/Gen),
   [Soss.jl](https://github.com/cscherrer/Soss.jl.git)
   [#157](https://github.com/alan-turing-institute/MLJ.jl/issues/157)
   [discourse
   thread](https://discourse.julialang.org/t/ppl-connection-to-mlj-jl/28736)
   [done](https://github.com/tlienart/SossMLJ.jl) but experimental and
   requires:

   - [ ] extension of probabilistic scoring functions to
	 "distributions" that can only be sampled.

-   [ ] Feature engineering (python featuretools?, recursive feature
	elimination?)
	[#426](https://github.com/alan-turing-institute/MLJ.jl/issues/426) [MLJModels #314](https://github.com/JuliaAI/MLJModels.jl/issues/314)


### Enhancing core functionality

-   [x] Iterative model control [#139](https://github.com/alan-turing-institute/MLJ.jl/issues/139). [Done](https://github.com/JuliaAI/MLJIteration.jl)

-   [ ] **†** Add more tuning
	strategies. See [here](https://github.com/JuliaAI/MLJTuning.jl#what-is-provided-here)
	for complete
	wish-list. Particular focus on:

	- [x] random search
	([#37](https://github.com/alan-turing-institute/MLJ.jl/issues/37))
	(done)

	- [x] Latin hypercube
	  [done](https://github.com/JuliaAI/MLJTuning.jl/pull/96)

	- [ ] Bayesian methods, starting with Gaussian Process methods a
	  la PyMC3. Some preliminary research done .

	- [ ] POC for AD-powered gradient descent [#74](https://github.com/alan-turing-institute/MLJ.jl/issues/74)

	- [ ] Tuning with adaptive resource allocation, as in
	  Hyperband. This might be implemented elegantly with the help of
	  the recent `IterativeModel` wrapper, which applies, in
	  particular to `TunedModel` instances [see
	  here](https://alan-turing-institute.github.io/MLJ.jl/dev/controlling_iterative_models/#Using-training-losses,-and-controlling-model-tuning).

	- [ ] Genetic algorithms
[#38](https://github.com/alan-turing-institute/MLJ.jl/issues/38)

	- [ ] Particle Swarm Optization (current WIP, GSoC project @lhnguyen-vn)

	- [ ] tuning strategies for non-Cartesian spaces of models [MLJTuning
	#18](https://github.com/JuliaAI/MLJTuning.jl/issues/18), architecture search, and other AutoML workflows

- [ ]  Systematic benchmarking, probably modeled on
	[MLaut](https://arxiv.org/abs/1901.03678) [#69](https://github.com/alan-turing-institute/MLJ.jl/issues/74)

- [ ]   Give `EnsembleModel` more extendible API and extend beyond bagging
	(boosting, etc) and migrate to separate repository?
	[#363](https://github.com/alan-turing-institute/MLJ.jl/issues/363)

- [ ]  **†** Enhance complex model compostition:

	- [x] Introduce a canned
	stacking model wrapper ([POC](https://alan-turing-institute.github.io/DataScienceTutorials.jl/getting-started/stacking/)). WIP @olivierlabayle

	- [ ] Get rid of macros for creating pipelines and possibly
	implement target transforms as wrapper ([MLJBase
	#594](https://github.com/alan-turing-institute/MLJ.jl/issues/594))
	WIP @CameronBieganek and @ablaom


### Broadening scope

- [ ] Integrate causal and counterfactual methods for, example,
  applications to FAIRness; see [this
  proposal](https://julialang.org/jsoc/gsoc/MLJ/#causal_and_counterfactual_methods_for_fairness_in_machine_learning)

- [ ] Explore possibility of closer integration of Interpretable
  Machine Learning approaches, such as Shapley values and lime; see
  [Shapley.jl](https://gitlab.com/ExpandingMan/Shapley.jl),
  [ShapML.jl](https://github.com/nredell/ShapML.jl),
  [ShapleyValues.jl](https://github.com/slundberg/ShapleyValues.jl),
  [Shapley.jl (older)](https://github.com/frycast/Shapley.jl) and
  [this
  proposal](https://julialang.org/jsoc/gsoc/MLJ/#interpretable_machine_learning_in_julia)

- [ ]  Spin-off a stand-alone measures (loss functions) package
	(currently
	[here](https://github.com/JuliaAI/MLJBase.jl/tree/master/src/measures)). Introduce
	measures for multi-targets [MLJBase
	#502](https://github.com/JuliaAI/MLJBase.jl/issues/502).

- [ ] Add sparse data support and better support for NLP models; we
	could use [NaiveBayes.jl](https://github.com/dfdx/NaiveBayes.jl)
	as a POC (currently wrapped only for dense input) but the API
	needs finalizing first
	{#731](https://github.com/alan-turing-institute/MLJ.jl/issues/731). Probably
	need a new SparseTables.jl package.

- [x] POC for implementation of time series models classification
	[#303](https://github.com/alan-turing-institute/MLJ.jl/issues/303),
	[ScientificTypesBase #14](https://github.com/JuliaAI/ScientificTypesBase.jl/issues/14) POC is [here](https://github.com/JuliaAI/TimeSeriesClassification.jl)

- [ ] POC for time series forecasting, along lines of sktime; probably needs [MLJBase
	#502](https://github.com/JuliaAI/MLJBase.jl/issues/502)
	first, and someone to finish [PR on time series
	CV](https://github.com/JuliaAI/MLJBase.jl/pull/331). See also [this proposal](https://julialang.org/jsoc/gsoc/MLJ/#time_series_forecasting_at_scale_-_speed_up_via_julia)

- [ ]   Add tools or separate repository for visualization in MLJ.

	- [x] Extend visualization of tuning plots beyond two-parameters
	[#85](https://github.com/alan-turing-institute/MLJ.jl/issues/85)
	(closed).
	[#416](https://github.com/alan-turing-institute/MLJ.jl/issues/416)
	[Done](https://github.com/JuliaAI/MLJTuning.jl/pull/121) but might be worth adding alternatives suggested in issue.

	- [ ] visualizing decision boundaries ? [#342](https://github.com/alan-turing-institute/MLJ.jl/issues/342)

	- [ ] provide visualizations that MLR3 provides via [mlr3viz](https://github.com/mlr-org/mlr3viz)

- [ ] Extend API to accomodate outlier detection, as provided by [OutlierDetection.jl](https://github.com/davnn/OutlierDetection.jl) [#780](https://github.com/alan-turing-institute/MLJ.jl/issues/780) WIP @davn and @ablaom

- [ ] Add more pre-processing tools:

  - [x] missing value imputation using Gaussina Mixture Model. Done,
	via addition of BetaML model, `MissingImputator`.

  - [ ] improve `autotype` method (from ScientificTypes), perhaps by
	training on large collection of datasets with manually labelled
	scitype schema.
	
- [ ] Add integration with [MLFlow](https://julialang.org/jsoc/gsoc/MLJ/#mlj_and_mlflow_integration); see [this proposal](https://julialang.org/jsoc/gsoc/MLJ/#mlj_and_mlflow_integration)

- [ ] Extend integration with [OpenML](https://www.openml.org) WIP @darenasc


### Scalability

- [ ]   Roll out data front-ends for all models after  [MLJBase
  #501](https://github.com/JuliaAI/MLJBase.jl/pull/501)
  is merged.

- [ ]  Online learning support and distributed data
	[#60](https://github.com/alan-turing-institute/MLJ.jl/issues/60)

- [ ]  DAG scheduling for learning network training
	[#72](https://github.com/alan-turing-institute/MLJ.jl/issues/72)
	(multithreading first?)

- [ ]  Automated estimates of cpu/memory requirements
	[#71](https://github.com/alan-turing-institute/MLJ.jl/issues/71)

- [x] Add multithreading to tuning [MLJTuning
  #15](https://github.com/JuliaAI/MLJTuning.jl/issues/15)
  [Done](https://github.com/JuliaAI/MLJTuning.jl/pull/42).
The MLJ.jl package is licensed under the MIT "Expat" License:

Copyright (c) 2020: Edoardo Barp, Anthony Blaom, Gergö Bohner, Valentin Churavy, Harvey Devereux, Thibaut Lienart,
Franz J Király, Mohammed Nook, Annika Stechemesser, Sebastian Vollmer; Mike Innes in partnership with Julia Computing.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 
# Code Organization

![](material/MLJ_stack.svg)

*Dependency chart for MLJ repositories. Repositories with dashed
connections do not currently exist but are planned/proposed.*

Repositories of some possible interest outside of MLJ, or beyond
its conventional use, are marked with a ⟂ symbol:

* [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) is the
  general user's point-of-entry for choosing, loading, composing,
  evaluating and tuning machine learning models. It pulls in most code
  from other repositories described below.  MLJ also hosts the [MLJ
  manual](src/docs) which documents functionality across the
  repositories, with the exception of ScientificTypesBase, and
  MLJScientific types which host their own documentation. (The MLJ
  manual and MLJTutorials do provide overviews of scientific types.)

* [MLJModelInterface.jl](https://github.com/JuliaAI/MLJModelInterface.jl)
  is a lightweight package imported by packages implementing MLJ's
  interface for their machine learning models. It's only dependencies
  are ScientificTypesBase.jl (which depends only on the standard
  library module `Random`) and
  [StatisticalTraits.jl](https://github.com/JuliaAI/StatisticalTraits.jl)
  (which depends only on ScientificTypesBase.jl).

* (⟂)
  [MLJBase.jl](https://github.com/JuliaAI/MLJBase.jl) is
  a large repository with two main purposes: (i) to give "dummy"
  methods defined in MLJModelInterface their intended functionality
  (which depends on third party packages, such as
  [Tables.jl](https://github.com/JuliaData/Tables.jl),
  [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
  and
  [CategoricalArrays.jl](https://github.com/JuliaData/CategoricalArrays.jl));
  and (ii) provide functionality essential to the MLJ user that has
  not been relegated to its own "satellite" repository for some
  reason. See the [MLJBase.jl
  readme](https://github.com/JuliaAI/MLJBase.jl) for a
  detailed description of MLJBase's contents.

* [MLJModels.jl](https://github.com/JuliaAI/MLJModels.jl)
  hosts the *MLJ model registry*, which contains metadata on all the
  models the MLJ user can search and load from MLJ. Moreover, it
  provides the functionality for **loading model code** from MLJ on
  demand. Finally, it furnishes some commonly used transformers for
  data pre-processing, such as `ContinuousEncoder` and `Standardizer`.

* [MLJTuning.jl](https://github.com/JuliaAI/MLJTuning.jl)
  provides MLJ's `TunedModel` wrapper for hyper-parameter
  optimization, including the extendable API for tuning strategies,
  and selected in-house implementations, such as `Grid` and
  `RandomSearch`.
  
* [MLJEnsembles.jl](https://github.com/JuliaAI/MLJEnsembles.jl)
  provides MLJ's `EnsembleModel` wrapper, for creating homogenous
  model ensembles.
  
* [MLJIteration.jl](https://github.com/JuliaAI/MLJIteration.jl)
  provides the `IteratedModel` wrapper for controlling iterative
  models (snapshots, early stopping criteria, etc)
  
* [MLJSerialization.jl](https://github.com/JuliaAI/MLJSerialization.jl)
  provides functionality for saving MLJ machines to file
  
* [MLJOpenML.jl](https://github.com/JuliaAI/MLJOpenML.jl) provides
  integration with the [OpenML](https://www.openml.org) data science
  exchange platform
  
* (⟂)
  [MLJLinearModels.jl](https://github.com/JuliaAI/MLJLinearModels.jl)
  is an experimental package for a wide range of julia-native penalized linear models
  such as Lasso, Elastic-Net, Robust regression, LAD regression,
  etc. 

* [MLJFlux.jl](https://github.com/FluxML/MLJFlux.jl) an experimental
  package for gradient-descent models, such as traditional
  neural-networks, built with
  [Flux.jl](https://github.com/FluxML/Flux.jl), in MLJ.
  
* (⟂)
  [ScientificTypesBase.jl](https://github.com/JuliaAI/ScientificTypesBase.jl)
  is an ultra lightweight package providing "scientific" types,
  such as `Continuous`, `OrderedFactor`, `Image` and `Table`. It's
  purpose is to formalize conventions around the scientific
  interpretation of ordinary machine types, such as `Float32` and
  `DataFrame`.
  
* (⟂)
  [ScientificTypes.jl](https://github.com/JuliaAI/ScientificTypes.jl)
  articulates MLJ's own convention for the scientific interpretation of
  data.
  
* (⟂)
  [StatisticalTraits.jl](https://github.com/JuliaAI/StatisticalTraits.jl)
  An ultra lightweight package defining fall-back implementations for
  a collection of traits possessed by statistical objects.

* (⟂)
  [DataScienceTutorials](https://github.com/alan-turing-institute/DataScienceTutorials.jl)
  collects tutorials on how to use MLJ, which are deployed
  [here](https://alan-turing-institute.github.io/DataScienceTutorials.jl/)
<div align="center">
    <img src="material/MLJLogo2.svg" alt="MLJ" width="200">
</div>

<h2 align="center">A Machine Learning Framework for Julia
<p align="center">
  <a href="https://github.com/alan-turing-institute/MLJ.jl/actions">
    <img src="https://github.com/alan-turing-institute/MLJ.jl/workflows/CI/badge.svg"
         alt="Build Status">
  </a>
  <a href="https://alan-turing-institute.github.io/MLJ.jl/dev/">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg"
         alt="Documentation">
  </a>
  <a href="https://mybinder.org/v2/gh/alan-turing-institute/MLJ.jl/master?filepath=binder%2FMLJ_demo.ipynb">
  <img src="https://mybinder.org/badge_logo.svg"
       alt="Binder">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yelllow"
       alt="bibtex">
  </a>
  <a href="BIBLIOGRAPHY.md">
    <img src="https://img.shields.io/badge/cite-BibTeX-blue"
       alt="bibtex">
  </a>

</p>
</h2>


MLJ (Machine Learning in Julia) is a toolbox written in Julia
providing a common interface and meta-algorithms for selecting,
tuning, evaluating, composing and comparing over [160 machine learning
models](https://alan-turing-institute.github.io/MLJ.jl/dev/list_of_supported_models/)
written in Julia and other languages.

**New to MLJ?** Start [here](https://alan-turing-institute.github.io/MLJ.jl/dev/).

**Integrating an existing machine learning model into the MLJ
framework?** Start [here](https://alan-turing-institute.github.io/MLJ.jl/dev/quick_start_guide_to_adding_models/).

MLJ was initially created as a Tools, Practices and Systems project at
the [Alan Turing Institute](https://www.turing.ac.uk/)
in 2019. Current funding is provided by a [New Zealand Strategic
Science Investment
Fund](https://www.mbie.govt.nz/science-and-technology/science-and-innovation/funding-information-and-opportunities/investment-funds/strategic-science-investment-fund/ssif-funded-programmes/university-of-auckland/)
awarded to the University of Auckland.

MLJ been developed with the support of the following organizations:

<div align="center">
    <img src="material/Turing_logo.png" width = 100/>
    <img src="material/UoA_logo.png" width = 100/>
    <img src="material/IQVIA_logo.png" width = 100/>
    <img src="material/warwick.png" width = 100/>
    <img src="material/julia.png" width = 100/>
</div>


### The MLJ Universe

The functionality of MLJ is distributed over a number of repositories
illustrated in the dependency chart below. These repositories live at
the [JuliaAI](https://github.com/JuliaAI) umbrella organization.

<div align="center">
    <img src="material/MLJ_stack.svg" alt="Dependency Chart">
</div>

*Dependency chart for MLJ repositories. Repositories with dashed
connections do not currently exist but are planned/proposed.*

<br>
<p align="center">
<a href="CONTRIBUTING.md">Contributing</a> &nbsp;•&nbsp; 
<a href="ORGANIZATION.md">Code Organization</a> &nbsp;•&nbsp;
<a href="ROADMAP.md">Road Map</a> 
</br>

#### Contributors

*Core design*: A. Blaom, F. Kiraly, S. Vollmer

*Lead contributor*: A. Blaom

*Active maintainers*: A. Blaom, S. Okon, T. Lienart, D. Aluthge



# Community code of conduct

Members of the MLJ community should not infringe on the basic needs of
fellow members for respect, equality, honesty, inclusion, integrity,
and safety.

Formally, the [NumFOCUS Code of
Conduct](https://numfocus.org/code-of-conduct) governs all [MLJ
community activity](#mlj-community-participation), as defined
below. Violations and other complaints are not handled by NumFocus
but should be reported to a member of the [Steering
Committee](#steering-committee), whose contact details appear
below. The Committee will respond to reports on a case-by-case basis.

## MLJ community activity

"MLJ community activity" includes participation in any of the
following forums:

- the GitHub repository
  [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) (including
  all issues, discussions, and pull requests)

- all GitHub repositories in the [JuliaAI](https://github.com/JuliaAI)
  organization (including all issues, discussions and pull requests)

- all [Julia Discourse](https://discourse.julialang.org/about) topic
  threads related to any software provided by any of the above
  repositories

- all conversations on Julia [Slack](https://slack.com/intl/en-nz/)
  channels related to any software provided by any of the above
  mentioned repositories


## Steering Committee

- Sebastian Vollmer: sjvollmer@gmail.com
- Mark Gahegan: m.gahegan@auckland.ac.nz
- Anthony Blaom: anthony.blaom@gmail.com
## Contributing to the MLJ machine learning project

Contributions to MLJ are most welcome. Queries can be made through
issues or the Julia [slack
channel](https://julialang.org/slack/), #MLJ. 

- [Road map](ROADMAP.md)

- [Code organization](ORGANIZATION.md)


### Conventions

We follow
[this](https://nvie.com/posts/a-successful-git-branching-model/) git
work-flow and, in particular, ask that **all pull requests be made to
the`dev` branch** of the appropriate repo, and not to `master`. This
includes changes to documentation. All pull requests onto `master`
come from `dev` and generally precede a tagged release.

Contributors are kindly requested to adhere to the
[Blue](https://github.com/invenia/BlueStyle) style guide, with line
widths capped at 80 characters.


### Very brief design overview

MLJ has a basement level *model* interface, which must be implemented
for each new learning algorithm. Formally, each model is a `mutable
struct` storing hyperparameters and the implementer defines
model-dispatched `fit` and `predict` methods; for details, see
[here](docs/src/adding_models_for_general_use.md). The general user
interacts using *machines* which bind models with data and have an
internal state reflecting the outcomes of applying `fit!` and
`predict` methods on them. The model interface is pure "functional";
the machine interface more "object-oriented".

A generalization of machine, called a *nodal* machine, is a key
element of *learning networks* which combine several models together,
and form the basis for specifying new composite model types. See
[here](https://alan-turing-institute.github.io/MLJ.jl/dev/composing_models/)
for more on these.

MLJ code is now spread over [multiple repositories](ORGANIZATION.md).


# Citing MLJ

An overview of MLJ design:


[![DOI](https://joss.theoj.org/papers/10.21105/joss.02704/status.svg)](https://doi.org/10.21105/joss.02704)

```bibtex
@article{Blaom2020,
  doi = {10.21105/joss.02704},
  url = {https://doi.org/10.21105/joss.02704},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {55},
  pages = {2704},
  author = {Anthony D. Blaom and Franz Kiraly and Thibaut Lienart and Yiannis Simillides and Diego Arenas and Sebastian J. Vollmer},
  title = {{MLJ}: A Julia package for composable machine learning},
  journal = {Journal of Open Source Software}
}
```

An in-depth view of MLJ's model composition design:

[![arXiv](https://img.shields.io/badge/arXiv-2012.15505-<COLOR>.svg)](https://arxiv.org/abs/2012.15505)

```bibtex
@misc{blaom2020flexible,
      title={Flexible model composition in machine learning and its implementation in {MLJ}}, 
      author={Anthony D. Blaom and Sebastian J. Vollmer},
      year={2020},
      eprint={2012.15505},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```
---
title: 'MLJ: A Julia package for composable machine learning'
tags:
  - Julia
  - Machine Learning
  - model composition
  - stacking
  - ensembling
  - hyper-parameter tuning
authors:
  - name: Anthony D. Blaom
    orcid: 0000-0001-6689-886X
    affiliation: "1, 2, 3"
  - name: Franz Kiraly
    orcid: 0000-0002-9254-793X
    affiliation: "3, 4"
  - name: Thibaut Lienart
    orcid: 0000-0003-0872-7098
    affiliation: 3
  - name: Yiannis Simillides
    orcid: 0000-0002-0287-8699
    affiliation: 7
  - name: Diego Arenas
    orcid: 0000-0001-7829-6102
    affiliation: 6
  - name: Sebastian J. Vollmer
    orcid: 0000-0002-9025-0753
    affiliation: "3, 5"
affiliations:
 - name: University of Auckland, New Zealand
   index: 1
 - name: New Zealand eScience Infrastructure, New Zealand
   index: 2
 - name: Alan Turing Institute, London, United Kingdom
   index: 3
 - name: University College London, United Kingdom
   index: 4
 - name: University of Warwick, United Kingdom
   index: 5
 - name: University of St Andrews, St Andrews, United Kingdom
   index: 6
 - name: Imperial College London, United Kingdom
   index: 7
date: XX April 2020
bibliography: paper.bib
---

# Introduction

Statistical modeling, and the building of complex modeling pipelines,
is a cornerstone of modern data science. Most experienced
data scientists rely on high-level open source modeling toolboxes -
such as sckit-learn [@Pedregosa2001; @Buitinck2013] (Python); Weka
[@Holmes1994] (Java); mlr [@BischlEtal2016] and caret [@Kuhn2008]
(R) - for quick blueprinting, testing, and creation of
deployment-ready models.  They do this by providing a common interface
to atomic components, from an ever-growing model zoo, and by providing
the means to incorporate these into complex work-flows. Practitioners
are wanting to build increasingly sophisticated composite models, as
exemplified in the strategies of top contestants in machine learning
competitions such as Kaggle.

MLJ (Machine Learning in Julia) [@MLJdocs] is a toolbox written in
Julia that provides a common interface and meta-algorithms for
selecting, tuning, evaluating, composing and comparing machine model
implementations written in Julia and other languages. More broadly,
the MLJ project hopes to bring cohesion and focus to a number of
emerging and existing, but previously disconnected, machine learning
algorithms and tools of high quality, written in Julia. A welcome
corollary of this activity will be increased cohesion and synergy
within the talent-rich communities developing these tools.

In addition to other novelties outlined below, MLJ aims to provide
first-in-class model composition capabilities. Guiding goals of the
MLJ project have been usability, interoperability, extensibility, code
transparency, and reproducibility.


## Why Julia?

Nowadays, even technically competent users of scientific software will
prototype solutions using a high-level language such as python, R, or
MATLAB. However, to achieve satisfactory performance, such code
typically wraps performance critical algorithms written in a second
low-level language, such as C or FORTRAN. Through its use of an
extensible, hierarchical system of abstract types, just-in-time
compilation, and by replacing object-orientation with multiple
dispatch, Julia solves the ubiquitous "two language problem"
[@BezansonEtal2017]. With less technical programming knowledge,
experts in a domain of application can get under the hood of machine
learning software to broaden its applicability, and innovation can be
accelerated through a dramatically reduced software development cycle.

As an example of the productivity boost provided by the
single-language paradigm, we cite the DifferentialEquations.jl package
[@RackauckasNie2017], which, in a few short years of development
by a small team of domain experts, became the best package in its
class [@Rackauckas2017].

Another major advantage of a single-language solution is the ability
to automatically differentiate (AD) functions from their code
representations. The Flux.jl package [@Innes2018], for example,
already makes use of AD to allow unparalleled flexibility in neural
network design.

As a new language, Julia is high-performance computing-ready, and its
superlative meta-programming features allow developers to create
domain-specific syntax for user interaction.


## Novelties

**Composability.** In line with current trends in "auto-ML",
MLJ's design is largely predicated on the importance of model
composability. Composite models share all the behavior of regular
models, constructed using a new flexible "learning networks"
syntax. Unlike the toolboxes cited above, MLJ's composition syntax is
flexible enough to define stacked models, with out-of-sample
predictions for the base learners, as well as more routine linear
pipelines, which can include target transformations that are
learned. As in mlr, hyper-parameter tuning is implemented as a model
wrapper.

**A unified approach to probabilistic predictions.** In MLJ,
probabilistic prediction is treated as a first class feature,
leveraging Julia's type system. In particular, unnecessary
case-distinctions, and ambiguous conventions regarding the
representation of probabilities, are avoided.

**Scientific types** To help users focus less on data representation
(e.g., `Float32`, `DataFrame`) and more on the intended *purpose* or
*interpretation* of data, MLJ articulates model data requirements
using *scientific types* [@ScientificTypesBase], such as "continuous",
"ordered factor" or "table".

**Connecting models directly to arbitrary data containers**. A
user can connect models directly to tabular data in a manifold of
in-memory and out-of-memory formats by using a universal table
interface provided by the Tables.jl package [@Quinn].

**Finding the right model.** A model registry gives the user access to
model metadata without the need to actually load code defining the
model implementation. This metadata includes the model's data
requirements, for example, as well as a load path to enable MLJ to
locate the model interface code. Users can readily match models to
machine learning tasks, facilitating searches for an optimal model, a
search that can be readily automated.

**Tracking classes of categorical variables.** Finally, with
the help of scientific types and the CategoricalArrays.jl package
[@CategoricalArrays], users are guided to create safe
representations of categorical data, in which the complete pool of
possible classes is embedded in the data representation, and
classifiers preserve this information when making predictions. This
avoids a pain-point familiar in frameworks that simply recast
categorical data using integers: evaluating a classifier on the test
target, only to find the test data includes classes not seen in the
training data. Preservation of the original labels for these classes
also facilitates exploratory data analysis and interpretability.


# Scientific types

A scientific type is an ordinary Julia type (generally without
instances) reserved for indicating how some data should be
interpreted. Some of these types are shown in \autoref{fig1}.

<!-- ![Part of the scientific type hierarchy.\label{fig1}](scitypesII.svg) -->
![Part of the scientific type hierarchy.\label{fig1}](scitypesII.png)

To the scientific types, MLJ adds a specific *convention* specifying a
scientific type for every Julia object. The convention is expressed
through a single method `scitype`. So, for example, `scitype(x)`
returns `Continuous` whenever the type of `x` is a subtype of Julia's
`AbstractFloat` type, as in `scitype(3.14) == Continuous`. A tabular
data structure satisfying the Tables.jl interface, will always have
type `Table{K}`, where the type parameter `K` is the union of all
column scientific types.  A `coerce` method recasts machine types to
have the desired scientific type (interpretation), and a `schema`
method summarizes the machine and scientific types of tabular data.

Since scientific types are also Julia types, Julia's advanced type
system means scientific types can be organized in a type hierarchy.
It is straightforward to check the compatibility of data with a
model's scientific requirements and methods can be dispatched on
scientific type just as they would on ordinary types.

# Flexible and compact work-flows for performance evaluation and tuning

To evaluate the performance of some `model` object (specifying
the hyper-parameters of some supervised learning algorithm) using some
specified `resampling` strategy, and measured against some
battery of performance `measures`,  one runs:


```julia
evaluate(model, X, y, 
         resampling=CV(nfolds=6), 
		 measures=[L2HingeLoss(), BrierScore()])
```

which has (truncated) output

 `measure` | `measurement` | `per_fold` 
-------------|-----------------|-------------
L2HingeLoss  | 1.4             | [0.485, 1.58, 2.06, 1.09, 2.18, 1.03]
BrierScore{UnivariateFinite} | -0.702 | [-0.242, -0.788, -1.03, -0.545, -1.09, -0.514] 

As in mlr, hyper-parameter optimization is realized as a model
wrapper, which transforms a base model into a "self-tuning" version of
that model. That is, tuning is is abstractly specified before being
executed. This allows tuning to be integrated into work-flows
(learning networks) in multiple ways. A well-documented tuning
interface [@MLJTuning] allows developers to easily extend available
hyper-parameter tuning strategies.

We now give an example of syntax for wrapping a model called
`forest_model` in a random search tuning strategy, using
cross-validation, and optimizing the mean square loss. The `model` in
this case is a composite model with an ordinary hyper-parameter called
`bagging_fraction` and a *nested* hyper-parameter `atom.n_subfeatures`
(where `atom` is another model). The first two lines of code define
ranges for these parameters.

```julia
r1 = range(forest_model, :(atom.n_subfeatures), lower=1, upper=9)
r2 = range(forest_model, :bagging_fraction, lower=0.4, upper=1.0)
self_tuning_forest_model = TunedModel(model=forest_model,
                                      tuning=RandomSearch(),
                                      resampling=CV(nfolds=6),
                                      range=[r1, r2],
                                      measure=LPDistLoss(2),
                                      n=25)
```

In this random search example, default priors are assigned to each
hyper-parameter, but options exist to customize these. Both resampling
and tuning have options for parallelization; Julia has first class
support for both distributed and multi-threaded parallelism.

# A unified approach to probabilistic predictions and their evaluation

MLJ puts probabilistic models and deterministic models on equal
footing. Unlike most most frameworks, a supervised model is either
*probabilistic* - meaning it's `predict` method returns a distribution
object - *or* it is *deterministic* - meaning it returns
objects of the same scientific type as the training observations. To
use a probabilistic model to make deterministic predictions one can
wrap the model in a pipeline with an appropriate post-processing
function, or use additional `predict_mean`, `predict_median`,
`predict_mode` methods to deal with the common use-cases.

A "distribution" object returned by a probabilistic predictor is one
that can be sampled (using Julia's `rand` method) and queried
for properties. Where possible the object is in fact a
`Distribution` object from the Distributions.jl package
[@LinEtal2020], for which an additional `pdf` method for
evaluating the distribution's probability density or mass function
will be implemented, in addition to `mode`, `mean`
and `median` methods (allowing MLJ's fallbacks for
`predict_mean`, etc, to work).

One important distribution *not* provided by Distributions.jl is a
distribution for finite sample spaces with *labeled* elements (called
`UnivariateFinite`) which additionally tracks all possible classes of
the categorical variable it is modeling, and not just those observed
in training data.

By predicting distributions, instead of raw probabilities or
parameters, MLJ avoids a common pain point, namely deciding and
agreeing upon a convention about how these should be represented:
Should a binary classifier predict one probability or two? Are we
using the standard deviation or the variance here? What's the protocol
for deciding the order of (unordered) classes? How should multi-target
predictions be combined?, etc.

A case-in-point concerns performance measures (metrics) for
probabilistic models, such as cross-entropy and Brier loss. All
built-in probabilistic measures provided by MLJ are passed a
distribution in their prediction slot.

For an overview on probabilistic supervised learning we refer to
[@Gressmann2018].


# Model interfaces

In MLJ a *model* is just a struct storing the hyper-parameters
associated with some learning algorithm suggested by the struct name
(e.g., `DecisionTreeClassifier`) and that is all.  MLJ provides a
basic *model interface*, to be implemented by new machine learning
models, which is functional in style, for simplicity and maximal
flexibility. In addition to a `fit` and optional `update` method, one
implements one or more operations, such as `predict`, `transform` and
`inverse_transform`, acting on the learned parameters returned by
`fit`.

The optional `update` method allows one to avoid unnecessary
repetition of code execution (warm restart). The three main use-cases
are:

- **Iterative models.** If the only change to a random forest model is
  an increase in the number of trees by ten, for example, then not all
  trees need to be retrained; only ten new trees need to be trained.

- **Data preprocessing.** Avoid overheads associated with data
  preprocessing, such as coercion of data into an algorithm-specific
  type.

- **Smart training of composite models.** When tuning a simple
  transformer-predictor pipeline model using a holdout set, for
  example, it is unnecessary to retrain the transformer if only the
  predictor hyper-parameters change. MLJ implements "smart" retraining
  of composite models like this by defining appropriate `update`
  methods.

In the future MLJ will add an `update_data` method to support
models that can carry out on-line learning.

Presently, the general MLJ user is encouraged to interact through a
*machine interface* which sits on top of the model
interface. This makes some work-flows more convenient but, more
significantly, introduces a syntax which is more natural in the
context of model composition (see below). A *machine* is a
mutable struct that binds a model to data at construction, as in
`mach = machine(model, data)`, and which stores learned
parameters after the user calls `fit!(mach, rows=...)`. To
retrain with new hyper-parameters, the user can mutate `model`
and repeat the `fit!` call.

The operations `predict`, `transform`, etc are overloaded for
machines, which is how the user typically uses them, as in the call
`predict(mach, Xnew)`.


# Flexible model composition

Several limitations surrounding model composition are increasingly
evident to users of the dominant machine learning software
platforms. The basic model composition interfaces provided by the
toolboxes mentioned in the Introduction all share one or more of the
following shortcomings, which do not exist in MLJ:

- Composite models do not inherit all the behavior of ordinary
  models.

- Composition is limited to linear (non-branching) pipelines.

- Supervised components in a linear pipeline can only occur at the
  end of the pipeline.

- Only static (unlearned) target transformations/inverse
  transformations are supported.

- Hyper-parameters in homogeneous model ensembles cannot be coupled.

- Model stacking, with out-of-sample predictions for base learners,
  cannot be implemented.

- Hyper-parameters and/or learned parameters of component models are
  not easily inspected or manipulated (in tuning algorithms, for
  example)
  
- Composite models cannot implement multiple operations, for example,
  both a `predict` and `transform` method (as in clustering models) or
  both a `transform` and `inverse_transform` method.


We now sketch MLJ's composition API, referring the reader to
[@Blaom_I] for technical details, and to the MLJ documentation
[@MLJdocs; @MLJtutorials] for examples that will clarify how the
composition syntax works in practice.

Note that MLJ also provides "canned" model composition for common use
cases, such as non-branching pipelines and homogeneous ensembles,
which are not discussed further here.

Specifying a new composite model type is in two steps, *prototyping*
and *export*.

## Prototyping

In prototyping the user defines a so-called *learning network*, by
effectively writing down the same code she would use if composing the
models "by hand". She does this using the machine syntax, with which
she will already be familiar, from the basic `fit!`/`predict`
work-flow for single models. There is no need for the user to provide
production training data in this process. A dummy data set suffices,
for the purposes of testing the learning network as it is built.

<!-- ![Specifying prediction and training flows in a simple learning network. The network shown combines a ridge regressor with a learned target transformation (Box Cox).\label{fig2}](target_transformerVERTICAL.svg) -->
![Specifying prediction and training flows in a simple learning network. The network shown combines a ridge regressor with a learned target transformation (Box Cox).\label{fig2}](target_transformerVERTICAL.png)

The upper panel of Figure \autoref{fig2} illustrates a simple learning
network in which a continuous target `y` is "normalized" using a
learned Box Cox transformation, producing `z`, while PCA dimension
reduction is applied to some features `X`, to obtain `Xr`. A Ridge
regressor, trained using data from `Xr` and `z`, is then applied to
`Xr` to make a target prediction `ẑ`. To obtain a final prediction
`ŷ`, we apply the *inverse* of the Box Cox transform, learned
previously, to `ẑ`.

The lower "training" panel of the figure shows the three machines
which will store the parameters learned in training - the Box Cox
exponent and shift (`machine1`), the PCA projection
(`machine2`) and the ridge model coefficients and intercept
(`machine3`). The diagram additionally indicates where machines
should look for training data, and where to access model
hyper-parameters (stored in `box_cox`, `PCA` and
`ridge_regressor`).

The only syntactic difference between composing "by hand" and building
a learning network is that the training data must be wrapped in
"source nodes" (which can be empty if testing is not required) and the
`fit!` calls can be omitted, as training is now lazy. Each data
"variable" in the manual work-flow is now a node of a directed acyclic
graph encoding the composite model architecture. Nodes are callable,
with a node call triggering lazy evaluation of the `predict`,
`transform` and other operations in the network. Instead of
calling `fit!` on every machine, a single call to `fit!`
on a *node* triggers training of all machines needed to call
that node, in appropriate order. As mentioned earlier, training such a
node is "smart" in the sense that hyper-parameter changes to a model
only trigger retraining of necessary machines. So, for example, there
is no need to retrain the Box Cox transformer in the preceding example
if only the ridge regressor hyper-parameters have changed.

The syntax, then, for specifying the learning network shown
\autoref{fig2} looks like this:

```julia
X = source(X_dummy)        # or just source()
y = source(y_dummy)        # or just source()

machine1 = machine(box_cox, y)
z = transform(machine1, y)

machine2 = machine(PCA, X)
Xr = transform(machine2, X)

machine3 = machine(ridge_regressor, Xr, z)
ẑ = predict(machine3, Xr)

ŷ = inverse_transform(machine1, ẑ)

fit!(ŷ)  # to test training on the dummy data
ŷ()      # to test prediction on the dummy data
```

Note that the machine syntax is a mechanism allowing for multiple
nodes to point to the same learned parameters of a model, as in the
learned target transformation/inverse transformation above. They also
allow multiple nodes to share the same model (hyper-parameters) as in
homogeneous ensembles. And different nodes can be accessed during
training and "prediction" modes of operation, as in stacking.

## Export

In the second step of model composition, the learning network is
"exported" as a new stand-alone composite model type, with the
component models appearing in the learning network becoming default
values for corresponding hyper-parameters of the composite. This new
type (which is unattached to any particular data) can be instantiated
and used just like any other MLJ model (tuned, evaluated, etc). Under
the hood, training such a model builds a learning network, so that
training is "smart". Defining a new composite model type requires
generating and evaluating code, but this is readily implemented using
Julia's meta-programming tools, i.e., executed by the user with a
simple macro call.

# Future directions

There are plans to: (i) grow the number of models; (ii) enhance core
functionality, particularly around hyper-parameter optimization
[@MLJTuning]; and (iii) broaden scope, particularly around
probabilistic programming models, time series, sparse data and natural
language processing. A more comprehensive road map is linked from the
MLJ repository [@MLJ].

# Acknowledgments

We acknowledge valuable conversations with Avik Sengupta, Mike Innes,
mlr author Bernd Bischl, and IQVIA's Yaqub Alwan and Gwyn Jones. Seed
funding for the MLJ project has been provided by the Alan Turing
Institute's Tools, Practices and Systems programme, with special thanks
to Dr James Hethering, its former Programme Director, and Katrina
Payne. Mathematics for Real-World Systems Centre for Doctoral Training
at the University of Warwick provided funding for students exploring
the Julia ML ecosystem, who created an initial proof-of-concept.

**Code contributors.** D. Aluthge, D. Arenas, E. Barp, C. Bieganek,
A. Blaom, G. Bohner, M. K. Borregaard, D. Buchaca, V. Churavy,
H. Devereux, M. Giordano, J. Hoffimann, T. Lienart, M. Nook,
Z. Nugent, S. Okon, P. Oleśkiewicz, J. Samaroo, A. Shridar,
Y. Simillides, A. Stechemesser, S. Vollmer

# References
> 
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
<!--
A clear and concise description of what the bug is.
-->

**To Reproduce**
<!--
Add a Minimal, Complete, and Verifiable example (for more details, see e.g. 
https://stackoverflow.com/help/mcve

If the code is too long, feel free to put it in a public gist and link
it in the issue: https://gist.github.com
-->

```julia
[Past your code here.]
```

**Expected behavior**
<!--
A clear and concise description of what you expected to happen.
-->

**Additional context**
<!--
Add any other context about the problem here.
-->

**Versions**
<details>

<!--
Please run the following snippet and paste the output here.
using MLJ
...
-->

</details>

<!-- Thanks for contributing! -->
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
```@meta
EditURL = "<unknown>/../../../../MLJ/examples/lightning_tour/lightning_tour.jl"
```

# Lightning tour of MLJ

*For a more elementary introduction to MLJ, see [Getting
Started](https://alan-turing-institute.github.io/MLJ.jl/dev/getting_started/).*

**Note.** Be sure this file has not been separated from the
accompanying Project.toml and Manifest.toml files, which should not
should be altered unless you know what you are doing. Using them,
the following code block instantiates a julia environment with a tested
bundle of packages known to work with the rest of the script:

````@example lightning_tour
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
````

Assuming Julia 1.7

In MLJ a *model* is just a container for hyper-parameters, and that's
all. Here we will apply several kinds of model composition before
binding the resulting "meta-model" to data in a *machine* for
evaluation, using cross-validation.

Loading and instantiating a gradient tree-boosting model:

````@example lightning_tour
using MLJ
MLJ.color_off()

Booster = @load EvoTreeRegressor # loads code defining a model type
booster = Booster(max_depth=2)   # specify hyper-parameter at construction
````

````@example lightning_tour
booster.nrounds=50               # or mutate post facto
booster
````

This model is an example of an iterative model. As is stands, the
number of iterations `nrounds` is fixed.

### Composition 1: Wrapping the model to make it "self-iterating"

Let's create a new model that automatically learns the number of iterations,
using the `NumberSinceBest(3)` criterion, as applied to an
out-of-sample `l1` loss:

````@example lightning_tour
using MLJIteration
iterated_booster = IteratedModel(model=booster,
                                 resampling=Holdout(fraction_train=0.8),
                                 controls=[Step(2), NumberSinceBest(3), NumberLimit(300)],
                                 measure=l1,
                                 retrain=true)
````

### Composition 2: Preprocess the input features

Combining the model with categorical feature encoding:

````@example lightning_tour
pipe = ContinuousEncoder |> iterated_booster
````

### Composition 3: Wrapping the model to make it "self-tuning"

First, we define a hyper-parameter range for optimization of a
(nested) hyper-parameter:

````@example lightning_tour
max_depth_range = range(pipe,
                        :(deterministic_iterated_model.model.max_depth),
                        lower = 1,
                        upper = 10)
````

Now we can wrap the pipeline model in an optimization strategy to make
it "self-tuning":

````@example lightning_tour
self_tuning_pipe = TunedModel(model=pipe,
                              tuning=RandomSearch(),
                              ranges = max_depth_range,
                              resampling=CV(nfolds=3, rng=456),
                              measure=l1,
                              acceleration=CPUThreads(),
                              n=50)
````

### Binding to data and evaluating performance

Loading a selection of features and labels from the Ames
House Price dataset:

````@example lightning_tour
X, y = @load_reduced_ames;
nothing #hide
````

Binding the "self-tuning" pipeline model to data in a *machine* (which
will additionally store *learned* parameters):

````@example lightning_tour
mach = machine(self_tuning_pipe, X, y)
````

Evaluating the "self-tuning" pipeline model's performance using 5-fold
cross-validation (implies multiple layers of nested resampling):

````@example lightning_tour
evaluate!(mach,
          measures=[l1, l2],
          resampling=CV(nfolds=5, rng=123),
          acceleration=CPUThreads())
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

# Target Transformations

Some supervised models work best if the target variable has been
standardized, i.e., rescaled to have zero mean and unit variance.
Such a target transformation is learned from values of the training
target variable. In particular, one generally learns a different
transformation when training on a proper subset of the training
data. Good data hygiene prescribes that a new transformation should
be computed each time the supervised model is trained on new
data - for example in cross-validation.

Additionally, one generally wants to *inverse* transform the
predictions of the supervised model, so that final target predictions
are on the original scale.

All these concerns are addressed by wrapping the supervised model
using `TransformedTargetModel`:

```@setup 123
using MLJ
MLJ.color_off()
```

```@example 123
Ridge = @load RidgeRegressor pkg=MLJLinearModels verbosity=0
ridge = Ridge()
ridge2 = TransformedTargetModel(ridge, target=Standardizer())
```
Note that the all the original hyper-parameters, as well as those of
the `Standardizer`, are accessible as nested hyper-parameters of the
wrapped model, which can be trained or evaluated like any other:

```@example 123
X, y = make_regression(rng=1234)
y = 10^6*y
mach = machine(ridge2, X, y)
fit!(mach, rows=1:60, verbosity=0)
predict(mach, rows=61:62)
```

Training and predicting using `ridge2` as above means:

1. Standardizing the target `y` using the first 60 rows to get a new target `z`

2. Training the original `ridge` model using the first 60 rows of `X` and `z`

3. Calling `predict` on the machine trained in Step 2 on rows `61:62` of `X`

4. Applying the inverse scaling learned in Step 1 to those predictions (to get the final output shown above)

Since both `ridge` and `ridge2` return predictions on the original
scale, we can meaningfully compare the corresponding mean absolute
errors and see that the wrapped model appears to be better:

```@example 123
evaluate(ridge, X, y, measure=mae)
```

```@example 123
evaluate(ridge2, X, y, measure=mae)
```

Ordinary functions can also be used in target transformations but an
inverse must be explicitly specified:

```@example 123
ridge3 = TransformedTargetModel(ridge, target=y->log.(y), inverse=z->exp.(z))
X, y = @load_boston
evaluate(ridge3, X, y, measure=mae)
```

Without the log transform (ie, using `ridge`) we get the poorer
`mae` of 3.9.

```@docs
TransformedTargetModel
```
# MLJ Cheatsheet


## Starting an interactive MLJ session

```@repl cheat
using MLJ
MLJ_VERSION # version of MLJ for this cheatsheet
```

## Model search and code loading

`info("PCA")` retrieves registry metadata for the model called "PCA"

`info("RidgeRegressor", pkg="MultivariateStats")` retrieves metadata
for "RidgeRegresssor", which is provided by multiple packages

`models()` lists metadata of every registered model.

`models("Tree")` lists models with "Tree" in the model or package name.

`models(x -> x.is_supervised && x.is_pure_julia)` lists all supervised models written in pure julia.

`models(matching(X))` lists all unsupervised models compatible with input `X`.

`models(matching(X, y))` lists all supervised models compatible with input/target `X/y`.

With additional conditions:

```julia
models() do model
    matching(model, X, y) &&
    model.prediction_type == :probabilistic &&
        model.is_pure_julia
end
```

`Tree = @load DecisionTreeClassifier pkg=DecisionTree` imports "DecisionTreeClassifier" type and binds it to `Tree`
`tree = Tree()` to instantiate a `Tree`. 

`tree2  = Tree(max_depth=2)` instantiates a tree with different hyperparameter

`Ridge = @load RidgeRegressor pkg=MultivariateStats` imports a type for a model provided by multiple packages

For interactive loading instead use `@iload`


## Scitypes and coercion

`scitype(x)` is the scientific type of `x`. For example `scitype(2.4) == Continuous`

![scitypes_small.png](img/scitypes_small.png)

type                                       | scitype
-------------------------------------------|----------------------------------
`AbstractFloat`                            | `Continuous`
`Integer`                                  | `Count`
`CategoricalValue` and `CategoricalString` | `Multiclass` or `OrderedFactor`
`AbstractString`                           | `Textual`

*Figure and Table for common scalar scitypes*

Use `schema(X)` to get the column scitypes of a table `X`

`coerce(y, Multiclass)` attempts coercion of all elements of `y` into scitype `Multiclass`

`coerce(X, :x1 => Continuous, :x2 => OrderedFactor)` to coerce columns `:x1` and `:x2` of table `X`.

`coerce(X, Count => Continuous)` to coerce all columns with `Count` scitype to `Continuous`.


## Ingesting data

Split the table `channing` into target `y` (the `:Exit` column) and
features `X` (everything else), after a seeded row shuffling:

```julia
using RDatasets
channing = dataset("boot", "channing")
y, X =  unpack(channing, ==(:Exit); rng=123)
```

Same as above but exclude `:Time` column from `X`:

```julia
using RDatasets
channing = dataset("boot", "channing")
y, X =  unpack(channing,
               ==(:Exit),            # y is the :Exit column
               !=(:Time);            # X is the rest, except :Time
               rng=123)
```

Splitting row indices into train/validation/test, with seeded shuffling:

`train, valid, test = partition(eachindex(y), 0.7, 0.2, rng=1234)` for 70:20:10 ratio

For a stratified split:

`train, test = partition(eachindex(y), 0.8, stratify=y)`

Split a table or matrix `X`, instead of indices:

`Xtrain, Xvalid, Xtest = partition(X, 0.5, 0.3, rng=123)` 

Getting data from [OpenML](https://www.openml.org):

`table = OpenML.load(91)`

Creating synthetic classification data:

`X, y = make_blobs(100, 2)` (also: `make_moons`, `make_circles`)

Creating synthetic regression data:

`X, y = make_regression(100, 2)`

## Machine construction

Supervised case:

`model = KNNRegressor(K=1)` and `mach = machine(model, X, y)`

Unsupervised case:

`model = OneHotEncoder()` and `mach = machine(model, X)`

## Fitting

`fit!(mach, rows=1:100, verbosity=1, force=false)` (defaults shown)


## Prediction

Supervised case: `predict(mach, Xnew)` or `predict(mach, rows=1:100)`

Similarly, for probabilistic models: `predict_mode`, `predict_mean` and `predict_median`.

Unsupervised case: `transform(mach, rows=1:100)` or `inverse_transform(mach, rows)`, etc.


## Inspecting objects

`@more` gets detail on last object in REPL

`params(model)` gets nested-tuple of all hyperparameters, even nested ones

`info(ConstantRegressor())`, `info("PCA")`, `info("RidgeRegressor",
pkg="MultivariateStats")` gets all properties (aka traits) of registered models

`info(rms)` gets all properties of a performance measure

`schema(X)` get column names, types and scitypes, and nrows, of a table `X`

`scitype(X)` gets scientific type of `X`

`fitted_params(mach)` gets learned parameters of fitted machine

`report(mach)` gets other training results (e.g. feature rankings)


## Saving and retrieving machines

`MLJ.save("trained_for_five_days.jlso", mach)` to save machine `mach`

`predict_only_mach = machine("trained_for_five_days.jlso")` to deserialize.


## Performance estimation

`evaluate(model, X, y, resampling=CV(), measure=rms, operation=predict, weights=..., verbosity=1)`

`evaluate!(mach, resampling=Holdout(), measure=[rms, mav], operation=predict, weights=..., verbosity=1)`

`evaluate!(mach, resampling=[(fold1, fold2), (fold2, fold1)], measure=rms)`

## Resampling strategies (`resampling=...`)

`Holdout(fraction_train=0.7, rng=1234)` for simple holdout

`CV(nfolds=6, rng=1234)` for cross-validation

`StratifiedCV(nfolds=6, rng=1234)` for stratified cross-validation

`TimeSeriesSV(nfolds=4)` for time-series cross-validation

or a list of pairs of row indices:

`[(train1, eval1), (train2, eval2), ... (traink, evalk)]`

## Tuning

### Tuning model wrapper

`tuned_model = TunedModel(model=…, tuning=RandomSearch(), resampling=Holdout(), measure=…, operation=predict, range=…)`

### Ranges for tuning (`range=...`)

If `r = range(KNNRegressor(), :K, lower=1, upper = 20, scale=:log)`

then `Grid()` search uses `iterator(r, 6) == [1, 2, 3, 6, 11, 20]`.

`lower=-Inf` and `upper=Inf` are allowed.

Non-numeric ranges: `r = range(model, :parameter, values=…)`

Nested ranges: Use dot syntax, as in `r = range(EnsembleModel(atom=tree), :(atom.max_depth), ...)`

Can specify multiple ranges, as in `range=[r1, r2, r3]`. For more range options do `?Grid` or `?RandomSearch`


### Tuning strategies

`RandomSearch(rng=1234)` for basic random search

`Grid(resolution=10)` or `Grid(goal=50)` for basic grid search

Also available: `LatinHyperCube`, `Explicit` (built-in), `MLJTreeParzenTuning`, `ParticleSwarm`, `AdaptiveParticleSwarm` (3rd-party packages)


#### Learning curves

For generating plot of performance against parameter specified by `range`:

`curve = learning_curve(mach, resolution=30, resampling=Holdout(), measure=…, operation=predict, range=…, n=1)`

`curve = learning_curve(model, X, y, resolution=30, resampling=Holdout(), measure=…, operation=predict, range=…, n=1)`

If using Plots.jl:

`plot(curve.parameter_values, curve.measurements, xlab=curve.parameter_name, xscale=curve.parameter_scale)`


## Controlling iterative models

Requires: `using MLJIteration`

`iterated_model = IteratedModel(model=…, resampling=Holdout(), measure=…, controls=…, retrain=false)`


### Controls

Increment training: `Step(n=1)`

Stopping: `TimeLimit(t=0.5)` (in hours), `NumberLimit(n=100)`, `NumberSinceBest(n=6)`, `NotANumber()`, `Threshold(value=0.0)`, `GL(alpha=2.0)`, `PQ(alpha=0.75, k=5)`, `Patience(n=5)`

Logging: `Info(f=identity)`, `Warn(f="")`, `Error(predicate, f="")`

Callbacks: `Callback(f=mach->nothing)`, `WithNumberDo(f=n->@info(n))`, `WithIterationsDo(f=i->@info("num iterations: $i"))`, `WithLossDo(f=x->@info("loss: $x"))`, `WithTrainingLossesDo(f=v->@info(v))`

Snapshots: `Save(filename="machine.jlso")`

Wraps: `MLJIteration.skip(control, predicate=1)`, `IterationControl.with_state_do(control)`


## Performance measures (metrics)

Do `measures()` to get full list.

`info(rms)` to list properties (aka traits) of the `rms` measure


## Transformers

Built-ins include: `Standardizer`, `OneHotEncoder`, `UnivariateBoxCoxTransformer`, `FeatureSelector`, `FillImputer`, `UnivariateDiscretizer`, `ContinuousEncoder`, `UnivariateTimeTypeToContinuous`

Externals include: `PCA` (in MultivariateStats), `KMeans`, `KMedoids` (in Clustering).

`models(m -> !m.is_supervised)` to get full list


## Ensemble model wrapper

`EnsembleModel(atom=…, weights=Float64[], bagging_fraction=0.8, rng=GLOBAL_RNG, n=100, parallel=true, out_of_bag_measure=[])`


## Target transformation wrapper

`TransformedTargetModel(model=ConstantClassifier(), target=Standardizer())`

## Pipelines

`pipe = (X -> coerce(X, :height=>Continuous)) |> OneHotEncoder |> KNNRegressor(K=3)` 

Unsupervised:

`pipe = Standardizer |> OneHotEncoder`

Concatenation:

`pipe1 |> pipe2` or `model |> pipe` or `pipe |> model`, etc


## Define a supervised learning network:

`Xs = source(X)`
`ys = source(y)`

... define further nodal machines and nodes ...

`yhat = predict(knn_machine, W, ys)` (final node)


## Exporting a learning network as stand-alone model:

Supervised, with final node `yhat` returning point-predictions:

```julia
@from_network machine(Deterministic(), Xs, ys; predict=yhat) begin
    mutable struct Composite
	    reducer=network_pca
		regressor=network_knn
    end
```

Here `network_pca` and `network_knn` are models appearing in the
learning network.

Supervised, with `yhat` final node returning probabilistic predictions:

```julia
@from_network machine(Probabilistic(), Xs, ys; predict=yhat) begin
    mutable struct Composite
        reducer=network_pca
        classifier=network_tree
    end
```

Unsupervised, with final node `Xout`:

```julia
@from_network machine(Unsupervised(), Xs; transform=Xout) begin
    mutable struct Composite
	    reducer1=network_pca
		reducer2=clusterer
    end
end
```UnivariateTimeTypeToContinuous
# More on Probabilistic Predictors

Although one can call `predict_mode` on a probabilistic binary
classifier to get deterministic predictions, a more flexible strategy
is to wrap the model using `BinaryThresholdPredictor`, as this allows
the user to specify the threshold probability for predicting a
positive class. This wrapping converts a probablistic classifer into a
deterministic one.

The positive class is always the second class returned when calling
`levels` on the training target `y`.

```@docs
MLJModels.BinaryThresholdPredictor
```
# [List of Supported Models](@id model_list)

MLJ provides access to to a wide variety of machine learning models.
We are always looking for
[help](https://github.com/alan-turing-institute/MLJ.jl/blob/master/CONTRIBUTING.md)
adding new models or testing existing ones.  Currently available
models are listed below; for the most up-to-date list, run `using MLJ;
models()`.

* *experimental*: indicates the package is fairly new and/or is under
  active development; you can help by testing these packages and
  making them more robust,
* *medium*: indicates the package is fairly mature but may benefit
  from optimisations and/or extra features; you can help by suggesting
  either,
* *high*: indicates the package is very mature and functionalities are
  expected to have been fairly optimised and tested.

| Package | Models | Maturity | Note
| ------- | ------ | -------- | ----
[BetaML.jl](https://github.com/sylvaticus/BetaML.jl) | DecisionTreeClassifier, DecisionTreeRegressor, KernelPerceptronClassifier, PegasosClassifier, PerceptronClassifier, RandomForestClassifier | medium |
[Clustering.jl](https://github.com/JuliaStats/Clustering.jl) | KMeans, KMedoids | high | †
[DecisionTree.jl](https://github.com/bensadeghi/DecisionTree.jl) | DecisionTreeClassifier, DecisionTreeRegressor, AdaBoostStumpClassifier, RandomForestClassifier, RandomForestRegressor | high | 
[EvoTrees.jl](https://github.com/Evovest/EvoTrees.jl) | EvoTreeRegressor, EvoTreeClassifier, EvoTreeCount, EvoTreeGaussian | medium | gradient boosting models
[GLM.jl](https://github.com/JuliaStats/GLM.jl) | LinearRegressor, LinearBinaryClassifier, LinearCountRegressor | medium | †
[LIBSVM.jl](https://github.com/mpastell/LIBSVM.jl) | LinearSVC, SVC, NuSVC, NuSVR, EpsilonSVR, OneClassSVM | high | also via ScikitLearn.jl
[LightGBM.jl](https://github.com/IQVIA-ML/LightGBM.jl) | LightGBMClassifier, LightGBMRegressor | high | 
[MLJFlux.jl](https://github.com/FluxML/MLJFlux.jl) | NeuralNetworkRegressor, NeuralNetworkClassifier, MultitargetNeuralNetworkRegressor, ImageClassifier | experimental |
[MLJLinearModels.jl](https://github.com/JuliaAI/MLJLinearModels.jl) | LinearRegressor, RidgeRegressor, LassoRegressor, ElasticNetRegressor, QuantileRegressor, HuberRegressor, RobustRegressor, LADRegressor, LogisticClassifier, MultinomialClassifier | experimental |
[MLJModels.jl](https://github.com/JuliaAI/MLJModels.jl) (built-in) | StaticTransformer, FeatureSelector, FillImputer, UnivariateStandardizer, Standardizer, UnivariateBoxCoxTransformer, OneHotEncoder, ContinuousEncoder, ConstantRegressor, ConstantClassifier, BinaryThreshholdPredictor | medium |
[MLJText.jl](https://github.com/JuliaAI/MLJText.jl) | TfidfTransformer, BM25Transformer, BagOfWordsTransformer | low |
[MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl) | LinearRegressor, MultitargetLinearRegressor, RidgeRegressor, MultitargetRidgeRegressor, PCA, KernelPCA, ICA, LDA, BayesianLDA, SubspaceLDA, BayesianSubspaceLDA, FactorAnalysis, PPCA | high | 
[NaiveBayes.jl](https://github.com/dfdx/NaiveBayes.jl) | GaussianNBClassifier, MultinomialNBClassifier, HybridNBClassifier | experimental |
[NearestNeighborModels.jl](https://github.com/JuliaAI/NearestNeighborModels.jl) | KNNClassifier, KNNRegressor, MultitargetKNNClassifier, MultitargetKNNRegressor | high |
[OutlierDetectionNeighbors.jl](https://github.com/OutlierDetectionJL/OutlierDetectionNeighbors.jl) | ABODDetector, COFDetector, DNNDetector, KNNDetector, LOFDetector | medium | 
[OutlierDetectionNetworks.jl](https://github.com/OutlierDetectionJL/OutlierDetectionNetworks.jl) | AEDetector, DSADDetector, ESADDetector | medium | 
[OutlierDetectionPython.jl](https://github.com/OutlierDetectionJL/OutlierDetectionPython.jl) | ABODDetector, CBLOFDetector, COFDetector, COPODDetector, HBOSDetector, IForestDetector, KNNDetector, LMDDDetector, LOCIDetector, LODADetector, LOFDetector, MCDDetector, OCSVMDetector, PCADetector, RODDetector, SODDetector, SOSDetector | high | 
[ParallelKMeans.jl](https://github.com/PyDataBlog/ParallelKMeans.jl) | KMeans | experimental |
[PartialLeastSquaresRegressor.jl](https://github.com/lalvim/PartialLeastSquaresRegressor.jl) | PLSRegressor, KPLSRegressor | experimental |
[ScikitLearn.jl](https://github.com/cstjean/ScikitLearn.jl) | ARDRegressor, AdaBoostClassifier, AdaBoostRegressor, AffinityPropagation, AgglomerativeClustering, BaggingClassifier, BaggingRegressor, BayesianLDA, BayesianQDA, BayesianRidgeRegressor, BernoulliNBClassifier, Birch, ComplementNBClassifier, DBSCAN, DummyClassifier, DummyRegressor, ElasticNetCVRegressor, ElasticNetRegressor, ExtraTreesClassifier, ExtraTreesRegressor, FeatureAgglomeration, GaussianNBClassifier, GaussianProcessClassifier, GaussianProcessRegressor, GradientBoostingClassifier, GradientBoostingRegressor, HuberRegressor, KMeans, KNeighborsClassifier, KNeighborsRegressor, LarsCVRegressor, LarsRegressor, LassoCVRegressor, LassoLarsCVRegressor, LassoLarsICRegressor, LassoLarsRegressor, LassoRegressor, LinearRegressor, LogisticCVClassifier, LogisticClassifier, MeanShift, MiniBatchKMeans, MultiTaskElasticNetCVRegressor, MultiTaskElasticNetRegressor, MultiTaskLassoCVRegressor, MultiTaskLassoRegressor, MultinomialNBClassifier, OPTICS, OrthogonalMatchingPursuitCVRegressor, OrthogonalMatchingPursuitRegressor, PassiveAggressiveClassifier, PassiveAggressiveRegressor, PerceptronClassifier, ProbabilisticSGDClassifier, RANSACRegressor, RandomForestClassifier, RandomForestRegressor, RidgeCVClassifier, RidgeCVRegressor, RidgeClassifier, RidgeRegressor, SGDClassifier, SGDRegressor, SVMClassifier, SVMLClassifier, SVMLRegressor, SVMNuClassifier, SVMNuRegressor, SVMRegressor, SpectralClustering, TheilSenRegressor | high | †
[TSVD.jl](https://github.com/JuliaLinearAlgebra/TSVD.jl) | TSVDTransformer | high | 
[XGBoost.jl](https://github.com/dmlc/XGBoost.jl) | XGBoostRegressor, XGBoostClassifier, XGBoostCount | high |

**Note** (†): Some models are missing and assistance is welcome to
complete the interface. Post a message on the Julia #mlj Slack channel
if you would like to help, thanks!
```@raw html
<script async defer src="https://buttons.github.io/buttons.js"></script>

<div style="font-size:1.25em;font-weight:bold;">
  <a href="about_mlj"
    style="color: #389826;">About</a>           &nbsp;|&nbsp;
  <a href="https://alan-turing-institute.github.io/MLJ.jl/dev/about_mlj/#Installation" 
    style="color: #389826;">Install</a>         &nbsp;|&nbsp;
  <a href="https://alan-turing-institute.github.io/MLJ.jl/dev/about_mlj/#Learning-to-use-MLJ"
    style="color: #389826;">Learn</a>    &nbsp;|&nbsp;
  <a href="mlj_cheatsheet" style="color: #9558B2;">Cheatsheet</a>       &nbsp;|&nbsp;
  <a href="common_mlj_workflows" style="color: #9558B2;">Workflows</a>  &nbsp;|&nbsp;
  <a href="https://github.com/alan-turing-institute/MLJ.jl/" style="color: #9558B2;">For Developers</a> &nbsp;|&nbsp;
  <a href="https://mybinder.org/v2/gh/alan-turing-institute/MLJ.jl/master?filepath=binder%2FMLJ_demo.ipynb" style="color: #9558B2;">Live Demo</a> &nbsp;|&nbsp;
  <a href="third_party_packages" style="color: #9558B2;">3rd Party Packages</a>
</div>

<span style="color: #9558B2;font-size:4.5em;">
MLJ</span>
<br>
<span style="color: #9558B2;font-size:2.25em;font-style:italic;">
A Machine Learning Framework for Julia</span>
```

To support MLJ development, please cite these works or star the repo:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02704/status.svg)](https://doi.org/10.21105/joss.02704)  [![arXiv](https://img.shields.io/badge/arXiv-2012.15505-<COLOR>.svg)](https://arxiv.org/abs/2012.15505)

```@raw html
<a class="github-button" 
  href="https://github.com/alan-turing-institute/MLJ.jl" 
  data-icon="octicon-star" 
  data-size="large" 
  data-show-count="true" 
  aria-label="Star alan-turing-institute/MLJ.jl on GitHub">
  Star</a>
```

### Basics
[Getting Started](@ref) | 
[Working with Categorical Data](@ref) | 
[Common MLJ Workflows](@ref) |
[Machines](@ref) |
[MLJ Cheatsheet](@ref) 

### Data
[Working with Categorical Data](@ref) | 
[Preparing Data](@ref) |
[Generating Synthetic Data](@ref) |
[OpenML Integration](@ref)

### Models
[Model Search](@ref model_search) |
[Loading Model Code](@ref) |
[Transformers and Other Unsupervised Models](@ref) |
[More on Probabilistic Predictors](@ref) |
[Composing Models](@ref) |
[Simple User Defined Models](@ref) |
[List of Supported Models](@ref model_list) |
[Third Party Packages](@ref) 

### Meta-algorithms
[Evaluating Model Performance](@ref) |
[Tuning Models](@ref) |
[Controlling Iterative Models](@ref) |
[Learning Curves](@ref)

### Composition
[Composing Models](@ref) |
[Linear Pipelines](@ref) |
[Target Transformations](@ref) |
[Homogeneous Ensembles](@ref) |
[Model Stacking](@ref) |

### Customization and Extension
[Simple User Defined Models](@ref) |
[Quick-Start Guide to Adding Models](@ref) |
[Adding Models for General Use](@ref) |
[Composing Models](@ref) |
[Internals](@ref internals_section) |
[Modifying Behavior](@ref)

### Miscellaneous
[Weights](@ref) |
[Acceleration and Parallelism](@ref) |
[Performance Measures](@ref) 



# Model Stacking

In a model stack, as introduced by [Wolpert
(1992)](https://www.sciencedirect.com/science/article/abs/pii/S0893608005800231),
an adjucating model learns the best way to combine the predictions of
multiple base models. In MLJ, such models are constructed using the
`Stack` constructor. To learn more about stacking and to see how to
construct a stack "by hand" using [Learning Networks](@ref), see [this Data Science in Julia
tutorial](https://juliaai.github.io/DataScienceTutorials.jl/getting-started/stacking/)

```@docs
MLJBase.Stack
```
# Composing Models

Three common ways of combining multiple models together
have out-of-the-box implementations in MLJ:

- [Linear Pipelines](@ref) - for unbranching chains that take the
  output of one model (e.g., dimension reduction, such as `PCA`) and
  make it the input of the next model in the chain (e.g., a
  classification model, such as `EvoTreeClassifier`). To include
  transformations of the target variable in a supervised pipeline
  model, see [Target Transformations](@ref).
  
- [Homogeneous Ensembles](@ref) - for blending the predictions of
  multiple supervised models all of the same type, but which receive different
  views of the training data to reduce overall variance. The technique
  is known as observation
  [bagging](https://en.wikipedia.org/wiki/Bootstrap_aggregating). Bagging
  decision trees, like a `DecisionTreeClassifier`, gives what is known
  as a *random forest*, although MLJ also provides several canned
  random forest models.
  
- [Model Stacking](@ref) - for combining the predictions of a smaller
  number of models of possibly *different* type, with the help of an
  adjudicating model.
  
We note that composite models share all of the functionality of
ordinary models. Their main novelty is that they include other models
as hyper-parameters.

Finally, MLJ provides a powerful way to combine machine models in
flexible *learning networks*. By wrapping training data in *source
nodes* before calling functions like `machine`, `predict` and
`transform`, a complicated user workflow which already combines
multiple models is transformed into a blueprint for a new stand-alone
composite model type. For example, MLJ's `Stack` model is implemented
using a learning network. The remainder of this page is devoted to
explaining this advanced feature.

## Learning Networks

Below is a practical guide to the MLJ implementantion of learning
networks, which have been described more abstractly in the article:

[Anthony D. Blaom and Sebastian J. Voller (2020): Flexible model
composition in machine learning and its implementation in MLJ.
Preprint, arXiv:2012.15505](https://arxiv.org/abs/2012.15505)

We assuming familiarity with the basics outlined in [Getting
Started](index.md). The syntax for building a learning network is
essentially an extension of the basic syntax but with data containers
replaced with nodes of a graph.

It is important to distinguish between *learning networks* and the
composite MLJ model types they are used to define.

A *learning network* is a directed acyclic graph whose nodes are
objects that can be called to obtained data, either for training a
machine, or for using as input to an *operation*. An operation is
either:

- *static*, that is, an ordinary function, such as such as `+`, `log` or `vcat`; or

- *dynamic*, that is, an operation such as `predict` or `transform`
  which is dispatched on both data *and* a training outcome attached
  to some machine.

Since the result of calling a node depends on the outcome of
training events (and may involve lazy evaluation) one may think of a
node as "dynamic" data, as opposed to the "static" data appearing in
an ordinary MLJ workflow.

Different operations can dispatch on the same machine (i.e., can
access a common set of learned parameters) and different machines
can point to the same model (allowing for hyperparameter coupling).

By contrast, an *exported learning network* is a
learning network exported as a stand-alone, re-usable `Model` object,
to which all the MLJ `Model` meta-algorithms can be applied
(ensembling, systematic tuning, etc).

By specifying data at the source nodes of a learning network, one can
use and test the learning network as it is defined, which is also a
good way to understand how learning networks work under the hood. This
data, if specified, is ignored in the export process, for the exported
composite model, like any other model, is not associated with any data
until wrapped in a machine.

In MLJ learning networks treat the flow of information during training
and prediction/transforming separately.


### Building a simple learning network

The diagram below depicts a learning network which standardizes the
input data `X`, learns an optimal Box-Cox transformation for the
target `y`, predicts new target values using ridge regression, and
then inverse-transforms those predictions to restore them to the
original scale. (This represents a model we could alternatively build
using the `TransformedTargetModel` wrapper and a pipeline.) Here we
have only dynamic operations, labelled blue. The machines are in
green. Notice that two operations both use `box`, which stores the
learned Box-Cox parameters. The lower "Training" panel
indicates which nodes are used to train each machine, and what model
each machine is associated with.

![](img/target_transformer.png)

Looking ahead, we note that the new composite model type we will
create later will be assigned a single hyperparameter `regressor`, and the
learning network model `RidgeRegressor(lambda=0.1)` will become this
parameter's default value. Since model hyperparameters are mutable,
this regressor can be changed to a different one (e.g.,
`HuberRegressor()`).

For testing purposes, we'll use a small synthetic data set:

```@setup 42
using MLJ
MLJ.color_off()
const KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels
```

```@example 42
using Statistics
import DataFrames

x1 = rand(300)
x2 = rand(300)
x3 = rand(300)
y = exp.(x1 - x2 -2x3 + 0.1*rand(300))
X = DataFrames.DataFrame(x1=x1, x2=x2, x3=x3)

train, test  = partition(eachindex(y), 0.8); # hide
```
Step one is to wrap the data in *source nodes*:

```@example 42
Xs = source(X)
ys = source(y)
```

*Note.* One can omit the specification of data at the source nodes (by
writing instead `Xs = source()` and `ys = source()`) and
still export the resulting network as a stand-alone model using the
@from_network macro described later; see the example under [Static
operations on nodes](@ref). However, one will be unable to fit
or call network nodes, as illustrated below.

The contents of a source node can be recovered by simply calling the
node with no arguments:

```@example 42
ys()[1:2]
```

We label the nodes that we will define according to their outputs in
the diagram. Notice that the nodes `z` and `yhat` use the same
machine, namely `box`, for different operations.

To construct the `W` node we first need to define the machine `stand`
that it will use to transform inputs.

```@example 42
stand_model = Standardizer()
stand = machine(stand_model, Xs)
```
Because `Xs` is a node, instead of concrete data, we can call
`transform` on the machine without first training it, and the result
is the new node `W`, instead of concrete transformed data:

```@example 42
W = transform(stand, Xs)
```

To get actual transformed data we *call* the node appropriately, which
will require we first train the node. Training a node, rather than a
machine, triggers training of *all* necessary machines in the network.


```@example 42
fit!(W, rows=train)
W()           # transform all data
W(rows=test ) # transform only test data
W(X[3:4,:])   # transform any data, new or old
```

If you like, you can think of `W` (and the other nodes we will define)
as "dynamic data": `W` is *data*, in the sense that it an be called
("indexed") on rows, but *dynamic*, in the sense the result depends on
the outcome of training events.

The other nodes of our network are defined similarly:

```@example 42
RidgeRegressor = @load RidgeRegressor pkg=MultivariateStats
box_model = UnivariateBoxCoxTransformer()  # for making data look normally-distributed
box = machine(box_model, ys)
z = transform(box, ys)

ridge_model = RidgeRegressor(lambda=0.1)
ridge =machine(ridge_model, W, z)
zhat = predict(ridge, W)

yhat = inverse_transform(box, zhat); 
```
We are ready to train and evaluate the completed network. Notice that
the standardizer, `stand`, is *not* retrained, as MLJ remembers that
it was trained earlier:


```@example 42
fit!(yhat, rows=train);
rms(y[test], yhat(rows=test)) # evaluate
```
We can change a hyperparameters and retrain:

```@example 42
ridge_model.lambda = 0.01
fit!(yhat, rows=train); 
```
And re-evaluate:

```@example 42
rms(y[test], yhat(rows=test))
```

> **Notable feature.** The machine, `ridge::Machine{RidgeRegressor}`, is retrained, because its underlying model has been mutated. However, since the outcome of this training has no effect on the training inputs of the machines `stand` and `box`, these transformers are left untouched. (During construction, each node and machine in a learning network determines and records all machines on which it depends.) This behavior, which extends to exported learning networks, means we can tune our wrapped regressor (using a holdout set) without re-computing transformations each time a `ridge_model` hyperparameter is changed.


### Learning network machines

As we show shortly, a learning network needs to be "exported" to create a
new stand-alone model type. Instances of that type can be bound with
data in a machine, which can then be evaluated, for example. Somewhat
paradoxically, one can wrap a learning network in a certain kind of
machine, called a *learning network machine*, **before** exporting it,
and in fact, the export process actually requires us to do so. Since a
composite model type does not yet exist, one constructs the machine
using a "surrogate" model, whose name indicates the ultimate model
supertype (`Deterministic`, `Probabilistic`, `Unsupervised` or
`Static`). This surrogate model has no fields.

Continuing with the example above:

```@example 42
surrogate = Deterministic()
mach = machine(surrogate, Xs, ys; predict=yhat);
```

Notice that a key-word argument declares which node is for making
predictions, and the arguments `Xs` and `ys` declare which source
nodes receive the input and target data. With `mach` constructed in
this way, the code

```@setup 42
fit!(mach, verbosity=0)
predict(mach, X[test,:]);
```

```julia
fit!(mach)
predict(mach, X[test,:]);
```

is equivalent to

```@setup 42
fit!(yhat, verbosity=0)
yhat(X[test,:]);
```

```julia
fit!(yhat)
yhat(X[test,:]);
```

Like ordinary machine, once can call `report(mach)` and
`fitted_params(mach)`. While it's main purpose is for export (see
below), this machine can even be evaluated:

```@example 42
evaluate!(mach, resampling=CV(nfolds=3), measure=LPLoss(p=2))
```

For more on constructing learning network machines, see
[`machine`](@ref).

A learning network machine can also include additional internal state
in it's report (and so in the report of the corresponding exported
model). See [Exposing internal state of a learning network](@ref) for
this advanced feature.


## Exporting a learning network as a stand-alone model

Having satisfied that our learning network works on the synthetic
data, we are ready to export it as a stand-alone model.


### Method I: The @from_network macro

Having defined a learning network machine, `mach`, as above, the
following code defines a new model subtype `WrappedRegressor <:
Supervised` with a single field `regressor`:

```julia
@from_network mach begin
	mutable struct WrappedRegressor
		regressor=ridge_model
	end
end
```

Note the declaration of the default value `ridge_model`, *which must
refer to an actual model appearing in the learning network*. It can be
typed, as in the alternative declaration below, which also declares
some traits for the type (as shown by `info(WrappedRegressor)`; see
also [Trait declarations](@ref)).

```julia
@from_network mach begin
	mutable struct WrappedRegressor
		regressor::Deterministic=ridge_model
	end
	input_scitype = Table(Continuous,Finite)
	target_scitype = AbstractVector{<:Continuous}
end

```

We can now create an instance of this type and apply the
meta-algorithms that apply to any MLJ model:

```julia
julia> composite = WrappedRegressor()
WrappedRegressor(
	regressor = RidgeRegressor(
			lambda = 0.01))

X, y = @load_boston;
evaluate(composite, X, y, resampling=CV(), measure=l2, verbosity=0)
```

Since our new type is mutable, we can swap the `RidgeRegressor` out
for any other regressor:

```
KNNRegressor = @load KNNRegressor
composite.regressor = KNNRegressor(K=7)
julia> composite
WrappedRegressor(regressor = KNNRegressor(K = 7,
										  algorithm = :kdtree,
										  metric = Distances.Euclidean(0.0),
										  leafsize = 10,
										  reorder = true,
										  weights = :uniform,),) @ 2…63
```


### Method II: Finer control (advanced)

This section describes an advanced feature that can be skipped on a
first reading.

In Method I above, only models appearing in the network will appear as
hyperparameters of the exported composite model. There is a second
more flexible method for exporting the network, which allows finer
control over the exported `Model` struct, and which also avoids
macros. The two steps required are:

- Define a new `mutable struct` model type.

- Wrap the learning network code in a model `fit` method.

Let's start with an elementary illustration in the learning network we
just exported using Method I.

The `mutable struct` definition looks like this:

```@example 42
mutable struct WrappedRegressor2 <: DeterministicComposite
	regressor
end

# keyword constructor
WrappedRegressor2(; regressor=RidgeRegressor()) = WrappedRegressor2(regressor)
nothing #hide
```

The other supertype options are `ProbabilisticComposite`,
`IntervalComposite`, `UnsupervisedComposite` and `StaticComposite`.

We now simply cut and paste the code defining the learning network
into a model `fit` method (as opposed to a machine `fit!` method):


```@example 42
function MLJ.fit(model::WrappedRegressor2, verbosity::Integer, X, y)
	Xs = source(X)
	ys = source(y)

	stand_model = Standardizer()
	stand = machine(stand_model, Xs)
	W = transform(stand, Xs)

	box_model = UnivariateBoxCoxTransformer()
	box = machine(box_model, ys)
	z = transform(box, ys)

	ridge_model = model.regressor        ###
	ridge =machine(ridge_model, W, z)
	zhat = predict(ridge, W)

	yhat = inverse_transform(box, zhat)

	mach = machine(Deterministic(), Xs, ys; predict=yhat)
	return!(mach, model, verbosity)
end
```

This completes the export process.

Notes:

- The line marked `###`, where the new exported model's hyperparameter
  `regressor` is spliced into the network, is the only modification to
  the previous code.

- After defining the network there is the additional step of
  constructing and fitting a learning network machine (see above).

- The last call in the function `return!(mach, model, verbosity)`
  calls `fit!` on the learning network machine `mach` and splits it
  into various pieces, as required by the MLJ model interface. See
  also the [`return!`](@ref) doc-string.

- **Important note** An MLJ `fit` method is not allowed to mutate its
  `model` argument.

> **What's going on here?** MLJ's machine interface is built atop a more primitive *[model](simple_user_defined_models.md)* interface, implemented for each algorithm. Each supervised model type (eg, `RidgeRegressor`) requires model `fit` and `predict` methods, which are called by the corresponding *machine* `fit!` and `predict` methods. We don't need to define a  model `predict` method here because MLJ provides a fallback which simply calls the `predict` on the learning network machine created in the `fit` method.

#### A composite model coupling component model hyper-parameters

We now give a more complicated example of a composite model which
exposes some parameters used in the network that are not simply
component models. The model combines a clustering model (e.g.,
`KMeans()`) for dimension reduction with ridge regression, but has the
following "coupling" of the hyper parameters: The ridge regularization
depends on the number of clusters used (with less regularization for a
greater number of clusters) and a user-specified "coupling"
coefficient `K`.

```@example 42
RidgeRegressor = @load RidgeRegressor pkg=MLJLinearModels

mutable struct MyComposite <: DeterministicComposite
	clusterer     # the clustering model (e.g., KMeans())
	ridge_solver  # a ridge regression parameter we want to expose
	K::Float64    # a "coupling" coefficient
end

function MLJ.fit(composite::Composite, verbosity, X, y)

	Xs = source(X)
	ys = source(y)

	clusterer = composite.clusterer
	k = clusterer.k

	clustererM = machine(clusterer, Xs)
	Xsmall = transform(clustererM, Xs)

	# the coupling: ridge regularization depends on number of
	# clusters (and the specified coefficient `K`):
	lambda = exp(-composite.K/clusterer.k)

	ridge = RidgeRegressor(lambda=lambda, solver=composite.ridge_solver)
	ridgeM = machine(ridge, Xsmall, ys)

	yhat = predict(ridgeM, Xsmall)

	mach = machine(Deterministic(), Xs, ys; predict=yhat)
	return!(mach, composite, verbosity)

end

kmeans = (@load KMeans pkg=Clustering)()
my_composite = MyComposite(kmeans, nothing, 0.5)
```

```@example 42
evaluate(my_composite, X, y, measure=MeanAbsoluteError(), verbosity=0)
```

## Static operations on nodes

Continuing to view nodes as "dynamic data", we can, in addition to
applying "dynamic" operations like `predict` and `transform` to nodes,
overload ordinary "static" (unlearned) operations as well. These
operations can be ordinary functions (with possibly multiple
arguments) or they could be functions *with parameters*, such as "take
a weighted average of two nodes", where the weights are
parameters. Here we address the simpler case of ordinary functions. For
the parametric case, see "Static transformers" in [Transformers and Other Unsupervised Models](@ref)

Let us first give a demonstration of operations that work
out-of-the-box. These include:

- addition and scalar multiplication

- `exp`, `log`, `vcat`, `hcat`

- tabularization (`MLJ.table`) and matrixification (`MLJ.matrix`)

As a demonstration of some of these, consider the learning network
below that: (i) One-hot encodes the input table `X`; (ii) Log
transforms the continuous target `y`; (iii) Fits specified K-nearest
neighbour and ridge regressor models to the data; (iv) Computes an
average of the individual model predictions; and (v) Inverse
transforms (exponentiates) the blended predictions.

Note, in particular, the lines defining `zhat` and `yhat`, which
combine several static node operations.

```@example 42
RidgeRegressor = @load RidgeRegressor pkg=MultivariateStats
KNNRegressor = @load KNNRegressor

Xs = source()
ys = source()

hot = machine(OneHotEncoder(), Xs)

# W, z, zhat and yhat are nodes in the network:

W = transform(hot, Xs) # one-hot encode the input
z = log(ys)            # transform the target

model1 = RidgeRegressor(lambda=0.1)
model2 = KNNRegressor(K=7)

mach1 = machine(model1, W, z)
mach2 = machine(model2, W, z)

# average the predictions of the KNN and ridge models:
zhat = 0.5*predict(mach1, W) + 0.5*predict(mach2, W)

# inverse the target transformation
yhat = exp(zhat)
```

Exporting this learning network as a stand-alone model:

```julia
@from_network machine(Deterministic(), Xs, ys; predict=yhat) begin
	mutable struct DoubleRegressor
		regressor1=model1
		regressor2=model2
	end
end
```

To deal with operations on nodes not supported out-of-the box, one
can use the `@node` macro. Supposing, in the preceding example, we
wanted the geometric mean rather than arithmetic mean. Then, the
definition of `zhat` above can be replaced with

```julia
yhat1 = predict(mach1, W)
yhat2 = predict(mach2, W)
gmean(y1, y2) = sqrt.(y1.*y2)
zhat = @node gmean(yhat1, yhat2)
```

There is also a `node` function, which would achieve the same in this way:

```julia
zhat = node((y1, y2)->sqrt.(y1.*y2), predict(mach1, W), predict(mach2, W))
```

### More `node` examples

Here are some examples taken from MLJ source
(at work in the example above) for overloading common operations for nodes:

```julia
Base.log(v::Vector{<:Number}) = log.(v)
Base.log(X::AbstractNode) = node(log, X)

import Base.+
+(y1::AbstractNode, y2::AbstractNode) = node(+, y1, y2)
+(y1, y2::AbstractNode) = node(+, y1, y2)
+(y1::AbstractNode, y2) = node(+, y1, y2)
```

Here `AbstractNode` is the common super-type of `Node` and `Source`.

And a final example, using the `@node` macro to row-shuffle a table:

```julia
using Random
X = (x1 = [1, 2, 3, 4, 5],
	 x2 = [:one, :two, :three, :four, :five])
rows(X) = 1:nrows(X)

Xs = source(X)
rs  = @node rows(Xs)
W = @node selectrows(Xs, @node shuffle(rs))

julia> W()
(x1 = [5, 1, 3, 2, 4],
 x2 = Symbol[:five, :one, :three, :two, :four],)

```

## Exposing internal state of a learning network 

This section describes an advanced feature. 

Suppose you have a learning network that you would like to export as
a new stand-alone model type `MyModel`. Having bound `MyModel` to some
data in a machine `mach`, you would like to arrange that
`report(mach)` will record some additional information about the
internal state of the learning network that was built internally, when
you called `fit!(mach)`. This is possible by specifying the relevant
node or nodes when constructing the associated learning network
machine (see [Learning network machines](@ref) above) by including
the `report` keyword argument, as in 

```julia
mach = machine(Probabilistic(), Xs, ys, predict=yhat, report=(mean=N1, stderr=N2))
```

Here, `yhat`, `N1` and `N2` are nodes in the learning network and
`mean` and `stderr` are the desired key names for the machine's
report. After training this machine, `report(mach).mean` will return
the value of `N1()` when the underlying learning network was trained,
while `report(mach).stderr` will return the value of `N2()`.


Note that as `N1` and `N2` are called with no arguments, they do not
see "production" data, which is a point of difference with the
`predict` node `yhat`, which is called on the production data `Xnew`
on a call such as `predict(mach, Xnew)`. However, this also means the
nodes can have multiple origin nodes (query `origins` for
details). This is indeed the case in the following dummy example,
recording a training error in the composite model report:

```julia
using MLJ

import MLJModelInterface

struct MyModel <: ProbabilisticComposite
    model
end

function MLJModelInterface.fit(composite::MyModel, verbosity, X, y)

    Xs = source(X)
    ys = source(y)

    mach = machine(composite.model, Xs, ys)
    yhat = predict(mach, Xs)
    e = @node auc(yhat, ys)   # <------  node whose state we wish to export

    network_mach = machine(Probabilistic(),
                           Xs,
                           ys,
                           predict=yhat,
                           report=(training_error=e,))  # <------ how we export additional node(s)

    return!(network_mach, composite, verbosity)
end

X, y = make_moons()
composite = MyModel(ConstantClassifier())
mach = machine(composite, X, y) |> fit!
err = report(mach).training_error    # <------ accesssing the node state

yhat = predict(mach, rows=:);
@assert err ≈ auc(yhat, y) # true
```

## The learning network API

Two new julia types are part of learning networks: `Source` and `Node`.

Formally, a learning network defines *two* labeled directed acyclic
graphs (DAG's) whose nodes are `Node` or `Source` objects, and whose
labels are `Machine` objects. We obtain the first DAG from directed
edges of the form $N1 -> N2$ whenever $N1$ is an *argument* of $N2$
(see below). Only this DAG is relevant when calling a node, as
discussed in examples above and below. To form the second DAG
(relevant when calling or calling `fit!` on a node) one adds edges for
which $N1$ is *training argument* of the the machine which labels
$N1$. We call the second, larger DAG, the *completed learning network*
(but note only edges of the smaller network are explicitly drawn in
diagrams, for simplicity).


### Source nodes

Only source nodes reference concrete data. A `Source` object has a
single field, `data`.

```@docs
source(X)
rebind!
sources
origins
```

### Nodes

The key components of a `Node` are:

- An *operation*, which will either be *static* (a fixed function) or
  *dynamic* (such as `predict` or `transform`, dispatched on a machine).

- A *machine* on which to dispatch the operation (void if the
  operation is static). The training arguments of the machine are
  generally other nodes.

- Upstream connections to other nodes (including source nodes)
  specified by *arguments* (one for each argument of the operation).

```@docs
node
```

```@docs
@node
```

```@docs
@from_network
```

```@docs
return!
```

See more on fitting nodes at [`fit!`](@ref) and [`fit_only!`](@ref).
# Working with Categorical Data

## Scientific types for discrete data

Recall that models articulate their data requirements using scientific
types (see [Getting Started](@ref) or the [ScientificTypes.jl
documentation](https://JuliaAI.github.io/ScientificTypes.jl/dev/)). There
are three scientific types discrete data can have: `Count`,
`OrderedFactor` and `Multiclass`.


### Count data

In MLJ you cannot use integers to represent (finite) categorical
data. Integers are reserved for discrete data you want interpreted as
`Count <: Infinite`:

```@example hut
using MLJ # hide
scitype([1, 4, 5, 6])
```

The `Count` scientific type includes things like the number of phone
calls, or city populations, and other "frequency" data of a generally
unbounded nature.

That said, you may have data that is theoretically `Count`, but which
you coerce to `OrderedFactor` to enable the use of more models,
trusting to your knowledge of how those models work to inform an
appropriate interpretation.


### OrderedFactor and Multiclass data

Other integer data, such as the number of an animal's legs, or number
of rooms of homes, are generally coerced to `OrderedFactor <:
Finite`. The other categorical scientific type is `Multiclass <:
Finite`, which is for *unordered* categorical data. Coercing data to
one of these two forms is discussed under [ Detecting and coercing
improperly represented categorical data](@ref) below.


### Binary data

There is no separate scientific type for binary data. Binary data is
either `OrderedFactor{2}` if ordered, and `Multiclass{2}` otherwise.
Data with type `OrderedFactor{2}` is considered to have an intrinsic
"positive" class, e.g., the outcome of a medical test, and the
"pass/fail" outcome of an exam. MLJ measures, such as `true_positive`
assume the *second* class in the ordering is the "positive"
class. Inspecting and changing order is discussed in the next section.

If data has type `Bool` it is considered `Count` data (as `Bool <:
Integer`) and generally users will want to coerce such data to
`Multiclass` or `OrderedFactor`.


## Detecting and coercing improperly represented categorical data

One inspects the scientific type of data using `scitype` as shown
above. To inspect all column scientific types in a table
simultaneously, use `schema`. (The `scitype(X)` of a table `X`
contains a condensed form of this information used in type dispatch;
see
[here](https://github.com/JuliaAI/ScientificTypesBase.jl#more-on-the-table-type).)

```@example hut
import DataFrames.DataFrame
X = DataFrame(
                 name       = ["Siri", "Robo", "Alexa", "Cortana"],
                 gender     = ["male", "male", "Female", "female"],
                 likes_soup = [true, false, false, true],
                 height     = [152, missing, 148, 163],
                 rating     = [2, 5, 2, 1],
                 outcome    = ["rejected", "accepted", "accepted", "rejected"])
schema(X)
```

Coercing a single column:

```@example hut
X.outcome = coerce(X.outcome, OrderedFactor)
```

The *machine* type of the result is a `CategoricalArray`. For more on
this type see [Under the hood:
CategoricalValue and CategoricalArray](@ref) below.

Inspecting the order of the levels:

```@example hut
levels(X.outcome)
```

Since we wish to regard "accepted" as the positive class, it should
appear second, which we correct with the `levels!` function:

```@example hut
levels!(X.outcome, ["rejected", "accepted"])
levels(X.outcome)
```

!!! warning "Changing levels of categorical data"

    The order of levels should generally be changed
    early in your data science work-flow and then not again. Similar
    remarks apply to *adding* levels (which is possible; see the
    [CategorialArrays.jl documentation](https://juliadata.github.io/CategoricalArrays.jl/stable/)). MLJ supervised and unsupervised models assume levels
    and their order do not change.

Coercing all remaining types simultaneously:

```@example hut
Xnew = coerce(X, :gender     => Multiclass,
                 :likes_soup => OrderedFactor,
                 :height     => Continuous,
                 :rating     => OrderedFactor)
schema(Xnew)
```

For `DataFrame`s there is also in-place coercion, using `coerce!`.


## Tracking all levels

The key property of vectors of scientific type `OrderedFactor` and
 `Multiclass` is that the pool of all levels is not lost when
separating out one or more elements:

```@example hut
v = Xnew.rating
```

```@example hut
levels(v)
```

```@example hut
levels(v[1:2])
```

```@example hut
levels(v[2])
```
By tracking all classes in this way, MLJ
avoids common pain points around categorical data, such as evaluating
models on an evaluation set, only to crash your code because classes appear
there which were not seen during training.

By drawing test, validation and training data from a common data
structure (as described in [Getting Started](@ref), for example) one
ensures that all possible classes of categorical variables are tracked
at all times. However, this does not mitigate problems with new
*production* data, if categorical features there are missing classes
or contain previously unseen classes.


## New or missing levels in production data

!!! warning 

    Unpredictable behaviour may result whenever `Finite` categorical data presents in a production set with different classes (levels) from those presented in training

Consider, for example, the following naive workflow:

```@example hut
# train a one-hot encoder on some data:
x = coerce(["black", "white", "white", "black"], Multiclass)
X = DataFrame(x=x)

model = OneHotEncoder()
mach = machine(model, X) |> fit!

# one-hot encode new data with missing classes:
xproduction = coerce(["white", "white"], Multiclass)
Xproduction = DataFrame(x=xproduction)
Xproduction == X[2:3,:]
```

So far, so good. But the following operation throws an error:

```julia
julia> transform(mach, Xproduction) == transform(mach, X[2:3,:])
ERROR: Found category level mismatch in feature `x`. Consider using `levels!` to ensure fitted and transforming features have the same category levels.
```

The problem here is that `levels(X.x)` and `levels(Xproduction.x)` are different:

```@example hut
levels(X.x)
```

```@example hut
levels(Xproduction.x)
```

This could be anticipated by the fact that the training and production
data have different schema:

```@example hut
schema(X)
```

```@example hut
schema(Xproduction)
```

One fix is to manually correct the levels of the production data:

```@example hut
levels!(Xproduction.x, levels(x))
transform(mach, Xproduction) == transform(mach, X[2:3,:])
```

Another solution is to pack all production data with dummy rows based
on the training data (subsequently dropped) to ensure there are no
missing classes. Currently MLJ contains no general tooling to check
and fix categorical levels in production data (although one can check
that training data and production data have the same schema, to ensure
the *number* of classes in categorical data is consistent).


## Extracting an integer representation of Finite data

Occasionally, you may really want an integer representation of data
that currently has scitype `Finite`. For example, you are developer
wrapping an algorithm from an external package for use in MLJ, and
that algorithm uses integer representations. Use the `int` method for
this purpose, and use `decoder` to construct decoders for reversing
the transformation:

```@example hut
v = coerce(["one", "two", "three", "one"], OrderedFactor);
levels!(v, ["one", "two", "three"]);
v_int = int(v)
```

```@example hut
d = decoder(v); # or decoder(v[1])
d.(v_int)
```

## Under the hood: CategoricalValue and CategoricalArray

In MLJ the objects with `OrderedFactor` or `Multiclass` scientific
type have machine type `CategoricalValue`, from the
[CategoricalArrays.jl]
(https://juliadata.github.io/CategoricalArrays.jl/stable/) package.
In some sense `CategoricalValue`s are an implementation detail users
can ignore for the most part, as shown above. However, you may want
some basic understanding of these types, and those implementing MLJ's
model interface for new algorithms will have to understand them. For
the complete API, see the CategoricalArrays.jl
[documentation](https://juliadata.github.io/CategoricalArrays.jl/stable/). Here
are the basics:


To construct an `OrderedFactor` or `Multiclass` vector directly from
raw labels, one uses `categorical`:

```@example hut
using CategoricalArrays # hide
v = categorical(['A', 'B', 'A', 'A', 'C'])
typeof(v)
```

(Equivalent to the idiomatically MLJ `v = coerce(['A', 'B', 'A', 'A',
'C']), Multiclass)`.)

```@example hut
scitype(v)
```

```@example hut
v = categorical(['A', 'B', 'A', 'A', 'C'], ordered=true, compress=true)
```

```@example hut
scitype(v)
```

When you index a `CategoricalVector` you don't get a raw label, but
instead an instance of `CategoricalValue`. As explained above, this
value knows the complete pool of levels from vector from which it
came. Use `get(val)` to extract the raw label from a value `val`.

Despite the distinction that exists between a value (element) and a
label, the two are the same, from the point of `==` and `in`:

```@julia
v[1] == 'A' # true
'A' in v    # true
```


## Probabilistic predictions of categorical data

Recall from [Getting Started](@ref) that probabilistic classfiers
ordinarily predict `UnivariateFinite` distributions, not raw
probabilities (which are instead accessed using the `pdf` method.)
Here's how to construct such a distribution yourself:

```@example hut
v = coerce(["yes", "no", "yes", "yes", "maybe"], Multiclass)
d = UnivariateFinite([v[1], v[2]], [0.9, 0.1])
```

Or, equivalently,

```@example hut
d = UnivariateFinite(["no", "yes"], [0.9, 0.1], pool=v)
```

This distribution tracks *all* levels, not just the ones to which you
have assigned probabilities:

```@example hut
pdf(d, "maybe")
```

However, `pdf(d, "dunno")` will throw an error.

You can declare `pool=missing`, but then `"maybe"` will not be tracked:

```@example hut
d = UnivariateFinite(["no", "yes"], [0.9, 0.1], pool=missing)
levels(d)
```

To construct a whole *vector* of `UnivariateFinite` distributions,
simply give the constructor a matrix of probabilities:

```@example hut
yes_probs = rand(5)
probs = hcat(1 .- yes_probs, yes_probs)
d_vec = UnivariateFinite(["no", "yes"], probs, pool=v)
```

Or, equivalently:

```@julia
d_vec = UnivariateFinite(["no", "yes"], yes_probs, augment=true, pool=v)
```

For more options, see [`UnivariateFinite`](@ref).
# OpenML Integration

The [OpenML platform](https://www.openml.org) provides an integration
platform for carrying out and comparing machine learning solutions
across a broad collection of public datasets and software platforms.

Integration with OpenML API is presently limited to querying and
downloading datasets.

Documentation is [here](https://juliaai.github.io/OpenML.jl/stable/).

# Known Issues

Routine issues are posted
[here](https://github.com/alan-turing-institute/MLJ.jl/issues). Below
are some longer term issues and limitations.

#### ScikitLearn/MKL issue

For users of Mac OS using Julia 1.3 or higher, using ScikitLearn
models can lead to unexpected MKL errors due to an issue not related
to MLJ. See
[this Julia Discourse discussion](https://discourse.julialang.org/t/julia-1-3-1-4-on-macos-and-intel-mkl-error/36469/2) 
and
[this issue](https://github.com/JuliaPackaging/BinaryBuilder.jl/issues/700)
for context. 

A temporary workaround for this issue is to force the installation of
an older version of the `OpenSpecFun_jll` library. To install an
appropriate version, activate your MLJ environment and run

```julia
  using Pkg;
  Pkg.add(PackageSpec(url="https://github.com/tlienart/OpenSpecFun_jll.jl"))
```

#### Serialization for composite models with component models with custom serialization

See
[here](https://github.com/alan-turing-institute/MLJ.jl/issues/678). Workaround:
Instead of `XGBoost` models (the chief known case) use models from the
pure Julia package `EvoTrees`.

# Index of Methods

```@index
```

# Performance Measures

In MLJ loss functions, scoring rules, sensitivities, and so on, are
collectively referred to as *measures*. These include re-exported loss
functions from the
[LossFunctions.jl](https://github.com/JuliaML/LossFunctions.jl)
library, overloaded to behave the same way as the built-in measures.

To see list all measures, run `measures()`.  Further measures for
probabilistic predictors, such as proper scoring rules, and for
constructing multi-target product measures, are planned.  If you'd like
to see measure added to MLJ, post a comment
[here](https://github.com/JuliaAI/MLJBase.jl/issues/299).g

*Note for developers:* The measures interface and the built-in
measures described here are defined in MLJBase, but will ultimately live
in a separate package.


## Using built-in measures

These measures all have the common calling syntax

```julia
measure(ŷ, y)
```

or

```julia
measure(ŷ, y, w)
```

where `y` iterates over observations of some target variable, and `ŷ`
iterates over predictions (`Distribution` or `Sampler` objects in the
probabilistic case). Here `w` is an optional vector of sample weights,
or a dictionary of class weights, when these are supported by the
measure.

```@repl losses_and_scores
using MLJ
y = [1, 2, 3, 4];
ŷ = [2, 3, 3, 3];
w = [1, 2, 2, 1];
rms(ŷ, y) # reports an aggregrate loss
l2(ŷ, y, w) # reports per observation losses
y = coerce(["male", "female", "female"], Multiclass)
d = UnivariateFinite(["male", "female"], [0.55, 0.45], pool=y);
ŷ = [d, d, d];
log_loss(ŷ, y)
```

The measures `rms`, `l2` and `log_loss` illustrated here are actually
        instances of measure *types*. For, example, `l2 = LPLoss(p=2)` and
`log_loss = LogLoss() = LogLoss(tol=eps())`. Common aliases are
provided:

```@repl losses_and_scores
cross_entropy
```

## Traits and custom measures

Notice that `l1` reports per-sample evaluations, while `rms`
only reports an aggregated result. This and other behavior can be
gleaned from measure *traits* which are summarized by the `info`
method:

```@repl losses_and_scores
info(l1)
```

Query the doc-string for a measure using the name of its type:

```@repl losses_and_scores
rms
@doc RootMeanSquaredError # same as `?RootMeanSqauredError
```

Use `measures()` to list all measures, and `measures(conditions...)` to
search for measures with given traits (as you would [query
models](model_search.md)). The trait `instances` list the actual
callable instances of a given measure type (typically aliases for the
default instance).

```@docs
measures(conditions...)
```

A user-defined measure in MLJ can be passed to the `evaluate!`
method, and elsewhere in MLJ, provided it is a function or callable
object conforming to the above syntactic conventions. By default, a
custom measure is understood to:

- be a loss function (rather than a score)

- report an aggregated value (rather than per-sample evaluations)

- be feature-independent

To override this behaviour one simply overloads the appropriate trait,
as shown in the following examples:

```@repl losses_and_scores
y = [1, 2, 3, 4];
ŷ = [2, 3, 3, 3];
w = [1, 2, 2, 1];
my_loss(ŷ, y) = maximum((ŷ - y).^2);
my_loss(ŷ, y)
my_per_sample_loss(ŷ, y) = abs.(ŷ - y);
MLJ.reports_each_observation(::typeof(my_per_sample_loss)) = true;
my_per_sample_loss(ŷ, y)
my_weighted_score(ŷ, y) = 1/mean(abs.(ŷ - y));
my_weighted_score(ŷ, y, w) = 1/mean(abs.((ŷ - y).^w));
MLJ.supports_weights(::typeof(my_weighted_score)) = true;
MLJ.orientation(::typeof(my_weighted_score)) = :score;
my_weighted_score(ŷ, y)
X = (x=rand(4), penalty=[1, 2, 3, 4]);
my_feature_dependent_loss(ŷ, X, y) = sum(abs.(ŷ - y) .* X.penalty)/sum(X.penalty);
MLJ.is_feature_dependent(::typeof(my_feature_dependent_loss)) = true
my_feature_dependent_loss(ŷ, X, y)
```

The possible signatures for custom measures are: `measure(ŷ, y)`,
`measure(ŷ, y, w)`, `measure(ŷ, X, y)` and `measure(ŷ, X, y, w)`, each
measure implementing one non-weighted version, and possibly a second
weighted version.

*Implementation detail:* Internally, every measure is evaluated using
the syntax

```julia
MLJ.value(measure, ŷ, X, y, w)
```
and the traits determine what can be ignored and how `measure` is actually called. If `w=nothing` then the non-weighted form of `measure` is
dispatched.

## Using measures from LossFunctions.jl

The [LossFunctions.jl](https://github.com/JuliaML/LossFunctions.jl)
package includes "distance loss" functions for `Continuous` targets,
and "marginal loss" functions for `Finite{2}` (binary) targets. While the
LossFunctions.jl interface differs from the present one (for, example
binary observations must be +1 or -1), MLJ has overloaded instances
of the LossFunctions.jl types to behave the same as the built-in
types.

Note that the "distance losses" in the package apply to deterministic
predictions, while the "marginal losses" apply to probabilistic
predictions.


## List of measures

All measures listed below have a doc-string associated with the measure's
*type*. So, for example, do `?LPLoss` not `?l2`.

```@setup losses_and_scores
using DataFrames
```

```@example losses_and_scores
ms = measures()
types = map(ms) do m
    m.name
end
instance = map(ms) do m m.instances end
table = (type=types, instances=instance)
DataFrame(table)
```


## Other performance related tools

In MLJ one computes a confusion matrix by calling an instance of the
`ConfusionMatrix` measure type on the data:

```@docs
ConfusionMatrix
```

```@docs
roc_curve
```
# Learning Curves

A *learning curve* in MLJ is a plot of some performance estimate, as a
function of some model hyperparameter. This can be useful when tuning
a single model hyperparameter, or when deciding how many iterations
are required for some iterative model. The `learning_curve` method
does not actually generate a plot, but generates the data needed to do
so.

To generate learning curves you can bind data to a model by
instantiating a machine. You can choose to supply all available data,
as performance estimates are computed using a resampling strategy,
defaulting to `Holdout(fraction_train=0.7)`.

```@example hooking
using MLJ
X, y = @load_boston;

atom = (@load RidgeRegressor pkg=MLJLinearModels)()
ensemble = EnsembleModel(model=atom, n=1000)
mach = machine(ensemble, X, y)

r_lambda = range(ensemble, :(model.lambda), lower=1e-1, upper=100, scale=:log10)
curve = MLJ.learning_curve(mach;
                           range=r_lambda,
                           resampling=CV(nfolds=3),
                           measure=MeanAbsoluteError())
```
```julia
using Plots
plot(curve.parameter_values,
     curve.measurements,
     xlab=curve.parameter_name,
     xscale=curve.parameter_scale,
     ylab = "CV estimate of RMS error")
```

![](img/learning_curve42.png)

In the case the `range` hyperparameter is the number of iterations in
some iterative model, `learning_curve` will not restart the training
from scratch for each new value, unless a non-holdout `resampling`
strategy is specified (and provided the model implements an
appropriate `update` method). To obtain multiple curves (that are
distinct) you will need to pass the name of the model random number
generator, `rng_name`, and specify the random number generators to be
used using `rngs=...` (an integer automatically generates the number
specified):

```@example hooking
atom.lambda= 7.3
r_n = range(ensemble, :n, lower=1, upper=50)
curves = MLJ.learning_curve(mach;
                            range=r_n,
                            measure=MeanAbsoluteError(),
                            verbosity=0,
                            rng_name=:rng,
                            rngs=4)
```

```julia
plot(curves.parameter_values,
     curves.measurements,
     xlab=curves.parameter_name,
     ylab="Holdout estimate of RMS error")
```

![](img/learning_curve_n.png)


## API reference

```@docs
MLJTuning.learning_curve
```
# Transformers and Other Unsupervised Models

Several unsupervised models used to perform common transformations,
such as one-hot encoding, are available in MLJ out-of-the-box. These
are detailed in [Built-in transformers](@ref) below.

A transformer is *static* if it has no learned parameters. While such
a transformer is tantamount to an ordinary function, realizing it as
an MLJ static transformer (subtype of `Static <: Unsupervised`) can be
useful, especially if the function depends on parameters the user
would like to manipulate (which become *hyper-parameters* of the
model). The necessary syntax for defining your own static transformers
is described in [Static transformers](@ref) below.

Some unsupervised models, such as clustering algorithms, have a
`predict` method in addition to a `transform` method. We give an
example of this in [Transformers that also predict](@ref)

Finally we note that models that fit a distribution, or more generally
a sampler object, to some data, which are sometimes viewed as
unsupervised, are treated in MLJ as *supervised* models. See [Models
that learn a probability distribution](@ref) for an example.


## Built-in transformers

```@docs
MLJModels.Standardizer
MLJModels.OneHotEncoder
MLJModels.ContinuousEncoder
MLJModels.FillImputer
MLJModels.FeatureSelector
MLJModels.UnivariateBoxCoxTransformer
MLJModels.UnivariateDiscretizer
MLJModels.UnivariateTimeTypeToContinuous
```


## Static transformers

The main use-case for static transformers is for insertion into
[Linear Pipelines](@ref) or other exported learning networks (see [Composing
Models](@ref)). If a static transformer has no hyper-parameters, it is
tantamount to an ordinary function. An ordinary function can be
inserted directly into a pipeline; the situation for learning
networks is only slightly more complicated; see [Static operations on nodes](@ref).

The following example defines a new model type `Averager` to perform
the weighted average of two vectors (target predictions, for
example). We suppose the weighting is normalized, and therefore
controlled by a single hyper-parameter, `mix`.

```@setup boots
using MLJ
```

```@example boots
mutable struct Averager <: Static
    mix::Float64
end

MLJ.transform(a::Averager, _, y1, y2) = (1 - a.mix)*y1 + a.mix*y2
```

*Important.* Note the sub-typing `<: Static`.

Such static transformers with (unlearned) parameters can have
arbitrarily many inputs, but only one output. In the single input case
an `inverse_transform` can also be defined. Since they have no real
learned parameters, you bind a static transformer to a machine without
specifying training arguments.

```@example boots
mach = machine(Averager(0.5)) |> fit!
transform(mach, [1, 2, 3], [3, 2, 1])
```

Let's see how we can include our `Averager` in a learning network (see
[Composing Models](@ref)) to mix the predictions of two regressors,
with one-hot encoding of the inputs:

```@example boots
X = source()
y = source() 

ridge = (@load RidgeRegressor pkg=MultivariateStats)()
knn = (@load KNNRegressor)()
averager = Averager(0.5)

hotM = machine(OneHotEncoder(), X)
W = transform(hotM, X) # one-hot encode the input

ridgeM = machine(ridge, W, y)
y1 = predict(ridgeM, W)

knnM = machine(knn, W, y)
y2 = predict(knnM, W)

averagerM= machine(averager)
yhat = transform(averagerM, y1, y2)
```

Now we export to obtain a `Deterministic` composite model and then
instantiate composite model

```julia
learning_mach = machine(Deterministic(), X, y; predict=yhat)
Machine{DeterministicSurrogate} @772 trained 0 times.
  args:
    1:	Source @415 ⏎ `Unknown`
    2:	Source @389 ⏎ `Unknown`


@from_network learning_mach struct DoubleRegressor
       regressor1=ridge
       regressor2=knn
       averager=averager
       end

composite = DoubleRegressor()
julia> composite = DoubleRegressor()
DoubleRegressor(
    regressor1 = RidgeRegressor(
            lambda = 1.0),
    regressor2 = KNNRegressor(
            K = 5,
            algorithm = :kdtree,
            metric = Distances.Euclidean(0.0),
            leafsize = 10,
            reorder = true,
            weights = :uniform),
    averager = Averager(
            mix = 0.5)) @301

```

which can be can be evaluated like any other model:

```julia
composite.averager.mix = 0.25 # adjust mix from default of 0.5
julia> evaluate(composite, (@load_reduced_ames)..., measure=rms)
Evaluating over 6 folds: 100%[=========================] Time: 0:00:00
┌───────────┬───────────────┬────────────────────────────────────────────────────────┐
│ _.measure │ _.measurement │ _.per_fold                                             │
├───────────┼───────────────┼────────────────────────────────────────────────────────┤
│ rms       │ 26800.0       │ [21400.0, 23700.0, 26800.0, 25900.0, 30800.0, 30700.0] │
└───────────┴───────────────┴────────────────────────────────────────────────────────┘
_.per_observation = [missing]
_.fitted_params_per_fold = [ … ]
_.report_per_fold = [ … ]
```


## Transformers that also predict

Some clustering algorithms learn to label data by identifying a
collection of "centroids" in the training data. Any new input
observation is labeled with the cluster to which it is closest (this
is the output of `predict`) while the vector of all distances from the
centroids defines a lower-dimensional representation of the
observation (the output of `transform`). In the following example a
K-means clustering algorithm assigns one of three labels 1, 2, 3 to
the input features of the iris data set and compares them with the
actual species recorded in the target (not seen by the algorithm).

```julia
import Random.seed!
seed!(123)

X, y = @load_iris;
KMeans = @load KMeans pkg=ParallelKMeans
kmeans = KMeans()
mach = machine(kmeans, X) |> fit!

# transforming:
Xsmall = transform(mach);
selectrows(Xsmall, 1:4) |> pretty
julia> selectrows(Xsmall, 1:4) |> pretty
┌─────────────────────┬────────────────────┬────────────────────┐
│ x1                  │ x2                 │ x3                 │
│ Float64             │ Float64            │ Float64            │
│ Continuous          │ Continuous         │ Continuous         │
├─────────────────────┼────────────────────┼────────────────────┤
│ 0.0215920000000267  │ 25.314260355029603 │ 11.645232464391299 │
│ 0.19199200000001326 │ 25.882721893491123 │ 11.489658693899486 │
│ 0.1699920000000077  │ 27.58656804733728  │ 12.674412792260142 │
│ 0.26919199999998966 │ 26.28656804733727  │ 11.64392098898145  │
└─────────────────────┴────────────────────┴────────────────────┘

# predicting:
yhat = predict(mach);
compare = zip(yhat, y) |> collect;
compare[1:8]
8-element Array{Tuple{CategoricalValue{Int64,UInt32},CategoricalString{UInt32}},1}:
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")
 (1, "setosa")

compare[51:58]
8-element Array{Tuple{CategoricalValue{Int64,UInt32},CategoricalString{UInt32}},1}:
 (2, "versicolor")
 (3, "versicolor")
 (2, "versicolor")
 (3, "versicolor")
 (3, "versicolor")
 (3, "versicolor")
 (3, "versicolor")
 (3, "versicolor")

compare[101:108]
8-element Array{Tuple{CategoricalValue{Int64,UInt32},CategoricalString{UInt32}},1}:
 (2, "virginica")
 (3, "virginica")
 (2, "virginica")
 (2, "virginica")
 (2, "virginica")
 (2, "virginica")
 (3, "virginica")
 (2, "virginica")
```
# Frequently Asked Questions

## Julia already has a great machine learning toolbox, ScitkitLearn.jl. Why MLJ?

An alternative machine learning toolbox for Julia users is
[ScikitLearn.jl](https://github.com/cstjean/ScikitLearn.jl). Initially
intended as a Julia wrapper for the popular python library
[scikit-learn](https://scikit-learn.org/stable/), ML algorithms
written in Julia can also implement the ScikitLearn.jl
API. Meta-algorithms (systematic tuning, pipelining, etc) remain
python wrapped code, however.

While ScikitLearn.jl provides the Julia user with access to a mature
and large library of machine learning models, the scikit-learn API on
which it is modeled, dating back to 2007, is not likely to
evolve significantly in the future. MLJ enjoys (or will enjoy) several
features that should make it an attractive alternative in the longer
term:

- **One language.** ScikitLearn.jl wraps Python code, which in turn
  wraps C code for performance-critical routines. A Julia machine
  learning algorithm that implements the MLJ model interface is 100%
  Julia. Writing code in Julia is almost as fast as Python and
  well-written Julia code runs almost as fast as C. Additionally, a
  single language design provides superior interoperability. For
  example, one can implement: (i) gradient-descent tuning of
  hyperparameters, using automatic differentiation libraries such as
  Flux.jl; and (ii) GPU performance boosts without major code
  refactoring, using CuArrays.jl.

- **Registry for model metadata.** In ScikitLearn.jl the list of
  available models, as well as model metadata (whether a model handles
  categorical inputs, whether is can make probabilistic predictions,
  etc) must be gleaned from documentation. In MLJ, this information is
  more structured and is accessible to MLJ via a searchable model
  registry (without the models needing to be loaded).

- **Flexible API for model composition.** Pipelines in scikit-learn
  are more of an afterthought than an integral part of the original
  design. By contrast, MLJ's user-interaction API was predicated on
  the requirements of a flexible "learning network" API, one that
  allows models to be connected in essentially arbitrary ways
  (such as Wolpert model stacks). Networks can be built
  and tested in stages before being exported as first-class
  stand-alone models. Networks feature "smart" training (only
  necessary components are retrained after parameter changes) and will
  eventually be trainable using a DAG scheduler.

- **Clean probabilistic API.** The scikit-learn API does not specify a
  universal standard for the form of probabilistic predictions. By
  fixing a probabilistic API along the lines of the
  [skpro](https://github.com/alan-turing-institute/skpro) project, MLJ
  aims to improve support for Bayesian statistics and probabilistic
  graphical models.

- **Universal adoption of categorical data types.** Python's
  scientific array library NumPy has no dedicated data type for
  representing categorical data (i.e., no type that tracks the pool of
  *all* possible classes). Generally scikit-learn models deal with
  this by requiring data to be relabeled as integers. However, the
  naive user trains a model on relabeled categorical data only to
  discover that evaluation on a test set crashes their code because a
  categorical feature takes on a value not observed in training. MLJ
  mitigates such issues by insisting on the use of categorical data
  types, and by insisting that MLJ model implementations preserve the
  class pools. If, for example, a training target contains classes in
  the pool that do not actually appear in the training set, a
  probabilistic prediction will nevertheless predict a distribution
  whose support includes the missing class, but which is appropriately
  weighted with probability zero.

Finally, we note that a large number of ScikitLearn.jl models are now
wrapped for use in MLJ.
# Tuning Models

MLJ provides several built-in and third-party options for optimizing a
model's hyper-parameters.  The quick-reference table below omits some
advanced key-word options.

tuning strategy | notes |package to import | package providing core algorithm
----------------|-------|------------------|----------------------------------
[`Grid`](@ref)`(goal=nothing, resolution=10)` | shuffled by default; `goal` is upper bound for number of grid points | MLJ.jl or MLJTuning.jl | [MLJTuning.jl](https://github.com/FluxML/model-zoo)
[`RandomSearch`](@ref)`(rng=GLOBAL_RNG)` | with customizable priors |MLJ.jl or MLJTuning.jl   | [MLJTuning.jl](https://github.com/FluxML/model-zoo)
[`LatinHypercube`](@ref)`(rng=GLOBAL_RNG)` | with discrete parameter support | MLJ.jl or MLJTuning.jl | [LatinHypercubeSampling](https://github.com/MrUrq/LatinHypercubeSampling.jl)
`MLJTreeParzenTuning()` | See this [example](https://github.com/IQVIA-ML/TreeParzen.jl/blob/master/docs/examples/simple_mlj_demo/simple_mlj_demo.md) for usage | TreeParzen.jl | [TreeParzen.jl](https://github.com/IQVIA-ML/TreeParzen.jl) (port to Julia of [hyperopt](http://hyperopt.github.io/hyperopt/))
`ParticleSwarm(n_particles=3, rng=GLOBAL_RNG)` | Standard Kennedy-Eberhart algorithm, plus discrete parameter support | MLJParticleSwarmOptimization.jl | [MLJParticleSwarmOptimization.jl](https://github.com/JuliaAI/MLJParticleSwarmOptimization.jl/)
`AdaptiveParticleSwarm(n_particles=3, rng=GLOBAL_RNG)` | Zhan et al. variant with automated swarm coefficient updates, plus discrete parameter support | MLJParticleSwarmOptimization.jl | [MLJParticleSwarmOptimization.jl](https://github.com/JuliaAI/MLJParticleSwarmOptimization.jl/)
`Explicit()` | For an [explicit list](@ref explicit) of models of varying type | MLJ.jl or MLJTuning.jl | [MLJTuning.jl](https://github.com/FluxML/model-zoo)

Below we illustrate hyperparameter optimization using the
[`Grid`](@ref), [`RandomSearch`](@ref), [`LatinHypercube`](@ref) and
`Explicit` tuning strategies.

## Overview

In MLJ model tuning is implemented as a model wrapper. After wrapping
a model in a tuning strategy and binding the wrapped model to data in
a machine called `mach`, calling `fit!(mach)` instigates a search for
optimal model hyperparameters, within a specified `range`, and then
uses all supplied data to train the best model. To predict using that
model, one then calls `predict(mach, Xnew)`. In this way the wrapped
model may be viewed as a "self-tuning" version of the unwrapped
model. That is, wrapping the model simply transforms certain
hyper-parameters into *learned* parameters.

A corollary of the tuning-as-wrapper approach is that evaluating the
performance of a `TunedModel` instance, using [`evaluate!`](@ref)
implies nested resampling. This approach is inspired by
[MLR](https://mlr.mlr-org.com/articles/tutorial/nested_resampling.html). See
also [below](@ref explicit).

In MLJ, tuning is an *iterative* procedure, with an iteration parameter
`n`, the total number of model instances to be evaluated.
Accordingly, tuning can be controlled using MLJ's
[`IteratedModel`](@ref MLJIteration.IteratedModel) wrapper. After
familiarizing oneself with the `TunedModel` wrapper described below,
see [Controlling model tuning](@ref) for more on this advanced
feature.

For a more in-depth overview of tuning in MLJ, or for implementation
details, see the [MLJTuning
documentation](https://github.com/JuliaAI/MLJTuning.jl). For a
complete list of options see the [`TunedModel`](@ref) doc-string
below.


## Tuning a single hyperparameter using a grid search (regression example)

```@example goof
using MLJ
X = MLJ.table(rand(100, 10));
y = 2X.x1 - X.x2 + 0.05*rand(100);
Tree = @load DecisionTreeRegressor pkg=DecisionTree verbosity=0;
tree = Tree()
```

Let's tune `min_purity_increase` in the model above, using a
grid-search. To do so we will use the simplest `range` object, a
one-dimensional range object constructed using the `range` method:

```@example goof
r = range(tree, :min_purity_increase, lower=0.001, upper=1.0, scale=:log);
self_tuning_tree = TunedModel(model=tree,
							  resampling=CV(nfolds=3),
							  tuning=Grid(resolution=10),
							  range=r,
							  measure=rms);
```

Incidentally, a grid is generated internally "over the range" by calling the
`iterator` method with an appropriate resolution:

```@example goof
iterator(r, 5)
```

Non-numeric hyperparameters are handled a little differently:

```@example goof
selector = FeatureSelector();
r2 = range(selector, :features, values = [[:x1,], [:x1, :x2]]);
iterator(r2)
```

Unbounded ranges are also permitted. See the [`range`](@ref) and
[`iterator`](@ref) docstrings below for details, and the
[`sampler`](@ref) docstring for generating random samples from
one-dimensional ranges (used internally by the [`RandomSearch`](@ref)
strategy).

Returning to the wrapped tree model:

```@example goof
mach = machine(self_tuning_tree, X, y);
fit!(mach, verbosity=0)
```

We can inspect the detailed results of the grid search with
`report(mach)` or just retrieve the optimal model, as here:

```@example goof
fitted_params(mach).best_model
```

For more detailed information, we can look at `report(mach)`, for example:

```@example goof
entry = report(mach).best_history_entry
```

Predicting on new input observations using the optimal model, *trained
on all the data* bound to `mach`:

```@example goof
Xnew  = MLJ.table(rand(3, 10));
predict(mach, Xnew)
```

Or predicting on some subset of the observations bound to `mach`:

```@example goof
test = 1:3
predict(mach, rows=test)
```

For tuning using only a subset `train` of all observation indices,
specify `rows=train` in the above `fit!` call. In that case the above
`predict` calls would be based on training the optimal model on all
`train` rows.


## A probabilistic classifier example

Tuning a classifier is not essentially different from tuning a
regressor. A common gotcha however is to overlook the distinction
between supervised models that make point predictions (subtypes of
`Deterministic`) and those that make probabilistic predictions
(subtypes of `Probabilistic`). The `DecisionTreeRegressor` model in
the preceding illustration was deterministic, so in this example will
consider a probabilistic classifier:

```@example goof
info("KNNClassifier").prediction_type
```

```@example goof
X, y = @load_iris
KNN = @load KNNClassifier verbosity=0
knn = KNN()
```

We'll tune the hyperparameter `K` in the model above, using a
grid-search once more:

```@example goof
K_range = range(knn, :K, lower=5, upper=20);
```

Since the model is probabilistic, we can choose either: (i) a
probabilistic measure, such as `brier_loss`; or (ii) use a
deterministic measure, such as `misclassification_rate` (which means
`predict_mean` is called instead of `predict` under the hood).

**Case (i) - probabilistic measure**:

```@example goof
self_tuning_knn = TunedModel(model=knn,
							 resampling = CV(nfolds=4, rng=1234),
							 tuning = Grid(resolution=5),
							 range = K_range,
							 measure=BrierLoss());

mach = machine(self_tuning_knn, X, y);
fit!(mach, verbosity=0);
```

**Case (ii) - deterministic measure**:

```@example goof
self_tuning_knn = TunedModel(model=knn,
							 resampling = CV(nfolds=4, rng=1234),
							 tuning = Grid(resolution=5),
							 range = K_range,
							 measure=MisclassificationRate())

mach = machine(self_tuning_knn, X, y);
fit!(mach, verbosity=0);
```

Let's inspect the best model and corresponding evaluation of the
metric in case (ii):

```@example goof
entry = report(mach).best_history_entry
```

```@example goof
entry.model.K
```

Recall that fitting `mach` also retrains the optimal
model on all available data. The following is therefore an optimal
model prediction based on all available data:

```@example goof
predict(mach, rows=148:150)
```

### Specifying a custom measure

Users may specify a custom loss or scoring function.  Suppose, for
example, we define a new scoring function `custom_accuracy` by

```@example goof
custom_accuracy(y,yhat) = mean(y .== yhat);
```

In tuning, scores are maximised, while losses are minimised. By
default, a custom measure is assumed to be a loss rather than a score,
so we must also declare

```@example goof
MLJ.orientation(::typeof(custom_accuracy)) = :score
```

For full details on constructing custom measures, see [Traits and custom
measures](@ref).


```@example goof
self_tuning_knn = TunedModel(model=knn,
							 resampling = CV(nfolds=4),
							 tuning = Grid(resolution=5),
							 range = K_range,
							 measure = [custom_accuracy, MulticlassFScore()],
							 operation = predict_mode);

mach = machine(self_tuning_knn, X, y)
fit!(mach, verbosity=0)
entry = report(mach).best_history_entry
```

```@example goof
entry.model.K
```

## Tuning multiple nested hyperparameters

The `forest` model below has another model, namely a
`DecisionTreeRegressor`, as a hyperparameter:

```@example goof
tree = Tree() # defined above
forest = EnsembleModel(model=tree)
```

Ranges for nested hyperparameters are specified using dot syntax. In
this case we will specify a `goal` for the total number of grid
points:

```@example goof
r1 = range(forest, :(model.n_subfeatures), lower=1, upper=9);
r2 = range(forest, :bagging_fraction, lower=0.4, upper=1.0);
self_tuning_forest = TunedModel(model=forest,
									  tuning=Grid(goal=30),
									  resampling=CV(nfolds=6),
									  range=[r1, r2],
									  measure=rms);

X = MLJ.table(rand(100, 10));
y = 2X.x1 - X.x2 + 0.05*rand(100);

mach = machine(self_tuning_forest, X, y);
fit!(mach, verbosity=0);
```

We can plot the grid search results:

```julia
using Plots
plot(mach)
```

![](img/tuning_plot.png)

Instead of specifying a `goal`, we can declare a global `resolution`,
which is overriden for a particular parameter by pairing it's range
with the resolution desired. In the next example, the default
`resolution=100` is applied to the `r2` field, but a resolution of `3`
is applied to the `r1` field. Additionally, we ask that the grid
points be randomly traversed, and the the total number of evaluations
be limited to 25.

```@example goof
tuning = Grid(resolution=100, shuffle=true, rng=1234)
self_tuning_forest = TunedModel(model=forest,
									  tuning=tuning,
									  resampling=CV(nfolds=6),
									  range=[(r1, 3), r2],
									  measure=rms,
									  n=25);
fit!(machine(self_tuning_forest, X, y), verbosity=0);
```

For more options for a grid search, see [`Grid`](@ref) below.

## Tuning using a random search

Let's attempt to tune the same hyperparameters using a `RandomSearch`
tuning strategy. By default, bounded numeric ranges like `r1` and `r2`
are sampled uniformly (before rounding, in the case of the integer
range `r1`). Positive unbounded ranges are sampled using a Gamma
distribution by default, and all others using a (truncated) normal
distribution.

```@example goof
self_tuning_forest = TunedModel(model=forest,
									  tuning=RandomSearch(),
									  resampling=CV(nfolds=6),
									  range=[r1, r2],
									  measure=rms,
									  n=25);
X = MLJ.table(rand(100, 10));
y = 2X.x1 - X.x2 + 0.05*rand(100);
mach = machine(self_tuning_forest, X, y);
fit!(mach, verbosity=0)
```

```julia
using Plots
plot(mach)
```

![](img/random_search_tuning_plot.png)

The prior distributions used for sampling each hyperparameter can be
customized, as can the global fallbacks. See the
[`RandomSearch`](@ref) doc-string below for details.


## Tuning using Latin hypercube sampling

One can also tune the hyperparameters using the `LatinHypercube`
tuning stragegy.  This method uses a genetic based optimization
algorithm based on the inverse of the Audze-Eglais function, using the
library
[`LatinHypercubeSampling.jl`](https://github.com/MrUrq/LatinHypercubeSampling.jl).

We'll work with the data `X`, `y` and ranges `r1` and `r2` defined
above and instatiate a Latin hypercube resampling strategy:

```@example goof
latin = LatinHypercube(gens=2, popsize=120)
```

Here `gens` is the number of generations to run the optimisation for
and `popsize` is the population size in the genetic algorithm. For
more on these and other `LatinHypercube` parameters, refer to the
[LatinHypercubeSampling.jl](https://github.com/MrUrq/LatinHypercubeSampling.jl)
documentation. Pay attention that `gens` and `popsize` are not to be
confused with the iteration parameter `n` in the construction of a
corresponding `TunedModel` instance, which specifies the total number
of models to be evaluated, independent of the tuning strategy.

For this illustration we'll add a third, nominal,  hyper-parameter:

```@example goof
r3 = range(forest, :(model.post_prune), values=[true, false]);
self_tuning_forest = TunedModel(model=forest,
									  tuning=latin,
									  resampling=CV(nfolds=6),
									  range=[r1, r2, r3],
									  measure=rms,
									  n=25);
mach = machine(self_tuning_forest, X, y);
fit!(mach, verbosity=0)
```

```julia
using Plots
plot(mach)
```
![](img/latin_hypercube_tuning_plot.png)


## [Comparing models of different type and nested cross-validation](@id explicit)

Instead of mutating hyperparameters of a fixed model, one can instead
optimise over an explicit list of models, whose types are allowed to
vary. As with other tuning strategies, evaluating the resulting
`TunedModel` itself implies nested resampling (e.g., nested
cross-validation) which we now examine in a bit more detail.

```@example goof
tree = (@load DecisionTreeClassifier pkg=DecisionTree verbosity=0)()
knn = (@load KNNClassifier pkg=NearestNeighborModels verbosity=0)()
models = [tree, knn]
nothing # hide
```

The following model is equivalent to best in `models` by using 3-fold
cross-validation:

```@example goof
multi_model = TunedModel(models=models,
						 resampling=CV(nfolds=3),
						 measure=log_loss,
						 check_measure=false)
nothing # hide
```

Note that there is no need to specify a `tuning` strategy or `range`
but we do specify `models` (plural) instead of `model`. Evaluating
`multi_model` implies nested cross-validation (each model gets
evaluated 2 x 3 times):

```@example goof
X, y = make_blobs()

e = evaluate(multi_model, X, y,
			 resampling=CV(nfolds=2),
			 measure=log_loss,
			 verbosity=6)
```

Now, for example, we can get the best model for the first fold out of
the two folds:

```@example goof
e.report_per_fold[1].best_model
```

And the losses in the outer loop (these still have to be matched to
the best performing model):

```@example goof
e.per_fold
```

It is also possible to get the results for the nested evaluations.
For example, for the first fold of the outer loop and the second model:

```@example goof
e.report_per_fold[2].history[1]
```

## API

```@docs
MLJBase.range
MLJBase.iterator
MLJBase.sampler
Distributions.fit(::Type{D}, ::MLJBase.NumericRange) where D<:Distributions.Distribution
MLJTuning.TunedModel
MLJTuning.Grid
MLJTuning.RandomSearch
MLJTuning.LatinHypercube
```
# Preparing Data

## Splitting data

MLJ has two tools for splitting data. To split data *vertically* (that
is, to split by observations) use [`partition`](@ref). This is commonly applied to a
vector of observation *indices*, but can also be applied to datasets
themselves, provided they are vectors, matrices or tables.

To split tabular data *horizontally* (i.e., break up a table based on
feature names) use [`unpack`](@ref).

```@docs
MLJBase.partition
MLJBase.unpack
```

## Bridging the gap between data type and model requirements

As outlined in [Getting Started](@ref), it is important that the
[scientific type](https://github.com/JuliaAI/ScientificTypesBase.jl) of
data matches the requirements of the model of interest. For example,
while the majority of supervised learning models require input
features to be `Continuous`, newcomers to MLJ are sometimes
surprised at the disappointing results of [model queries](@ref
model_search) such as this one:

```@setup poot
using MLJ
```
```@example poot
X = (height   = [185, 153, 163, 114, 180],
     time     = [2.3, 4.5, 4.2, 1.8, 7.1],
     mark     = ["D", "A", "C", "B", "A"],
     admitted = ["yes", "no", missing, "yes"]);
y = [12.4, 12.5, 12.0, 31.9, 43.0]
models(matching(X, y))
```

Or are unsure about the source of the following warning:

```julia
Tree = @load DecisionTreeRegressor pkg=DecisionTree verbosity=0
tree = Tree();
julia> machine(tree, X, y)

julia> machine(tree, X, y)
┌ Warning: The scitype of `X`, in `machine(model, X, ...)` is incompatible with `model=DecisionTreeRegressor @378`:                                                                
│ scitype(X) = Table{Union{AbstractVector{Continuous}, AbstractVector{Count}, AbstractVector{Textual}, AbstractVector{Union{Missing, Textual}}}}
│ input_scitype(model) = Table{var"#s46"} where var"#s46"<:Union{AbstractVector{var"#s9"} where var"#s9"<:Continuous, AbstractVector{var"#s9"} where var"#s9"<:Count, AbstractVector{var"#s9"} where var"#s9"<:OrderedFactor}.
└ @ MLJBase ~/Dropbox/Julia7/MLJ/MLJBase/src/machines.jl:103
Machine{DecisionTreeRegressor,…} @198 trained 0 times; caches data
  args: 
    1:  Source @628 ⏎ `Table{Union{AbstractVector{Continuous}, AbstractVector{Count}, AbstractVector{Textual}, AbstractVector{Union{Missing, Textual}}}}`
    2:  Source @544 ⏎ `AbstractVector{Continuous}`
```

The meaning of the warning is:

- The input `X` is a table with column scitypes `Continuous`, `Count`, and `Textual` and `Union{Missing, Textual}`, which can also see by inspecting the schema:
	
```@example poot
schema(X)
```

- The model requires a table whose column element scitypes subtype `Continuous`, an incompatibility.

### Common data preprocessing workflows

There are two tools for addressing data-model type mismatches like the
above, with links to further documentation given below:

**Scientific type coercion:** We coerce machine types to obtain the
intended scientific interpretation. If `height` in the above example
is intended to be `Continuous`, `mark` is supposed to be
`OrderedFactor`, and `admitted` a (binary) `Multiclass`, then we can
do
        
  
```@example poot
X_coerced = coerce(X, :height=>Continuous, :mark=>OrderedFactor, :admitted=>Multiclass);
schema(X_coerced)
```

**Data transformations:** We carry out conventional data
transformations, such as missing value imputation and feature
encoding:
  
```@example poot
imputer = FillImputer()
mach = machine(imputer, X_coerced) |> fit!
X_imputed = transform(mach, X_coerced);
schema(X_imputed)
```

```@example poot
encoder = ContinuousEncoder()
mach = machine(encoder, X_imputed) |> fit!
X_encoded = transform(mach, X_imputed)
```

```@example poot
schema(X_encoded)
```

Such transformations can also be combined in a pipeline; see [Linear
Pipelines](@ref).


## Scientific type coercion

Scientific type coercion is documented in detail at
[ScientificTypesBase.jl](https://github.com/JuliaAI/ScientificTypesBase.jl). See
also the tutorial at the [this MLJ
Workshop](https://github.com/ablaom/MachineLearningInJulia2020)
(specifically,
[here](https://github.com/ablaom/MachineLearningInJulia2020/blob/master/tutorials.md#fixing-scientific-types-in-tabular-data))
and [this Data Science in Julia
tutorial](https://alan-turing-institute.github.io/DataScienceTutorials.jl/data/scitype/).

Also relevant is the section, [Working with Categorical Data](@ref).


## Data transformation

MLJ's Built-in transformers are documented at [Transformers and Other Unsupervised Models](@ref). The most relevant in the present context
  are: [`ContinuousEncoder`](@ref), [`OneHotEncoder`](@ref),
  [`FeatureSelector`](@ref) and [`FillImputer`](@ref). A Gaussian
  mixture models imputer is provided by BetaML, which can be loaded
       with

```julia
MissingImputator = @load MissingImputator pkg=BetaML
```

[This MLJ
Workshop](https://github.com/ablaom/MachineLearningInJulia2020), and the "End-to-end
examples" in [Data Science in Julia
tutorials](https://alan-turing-institute.github.io/DataScienceTutorials.jl/)
give further illustrations of data preprocessing in MLJ.
# Generating Synthetic Data

MLJ has a set of functions - `make_blobs`, `make_circles`,
`make_moons` and `make_regression` (closely resembling functions in
[scikit-learn](https://scikit-learn.org/stable/datasets/index.html#generated-datasets)
of the same name) - for generating synthetic data sets. These are
useful for testing machine learning models (e.g., testing user-defined
composite models; see [Composing Models](@ref))


##  Generating Gaussian blobs

```@docs
make_blobs
```


```@example foggy
using MLJ, DataFrames
X, y = make_blobs(100, 3; centers=2, cluster_std=[1.0, 3.0])
dfBlobs = DataFrame(X)
dfBlobs.y = y
first(dfBlobs, 3)
```

```julia
using VegaLite
dfBlobs |> @vlplot(:point, x=:x1, y=:x2, color = :"y:n") 
```




![svg](img/output_4_0.svg)




```julia
dfBlobs |> @vlplot(:point, x=:x1, y=:x3, color = :"y:n") 
```




![svg](img/output_5_0.svg)



##  Generating concentric circles

```@docs
make_circles
```


```@example foggy
using MLJ, DataFrames
X, y = make_circles(100; noise=0.05, factor=0.3)
dfCircles = DataFrame(X)
dfCircles.y = y
first(dfCircles, 3)
```



```julia
using VegaLite
dfCircles |> @vlplot(:circle, x=:x1, y=:x2, color = :"y:n") 
```




![svg](img/output_8_0.svg)



##  Sampling from two interleaved half-circles

```@docs
make_moons
```


```@example foggy
using MLJ, DataFrames
X, y = make_moons(100; noise=0.05)
dfHalfCircles = DataFrame(X)
dfHalfCircles.y = y
first(dfHalfCircles, 3)
```




```julia
using VegaLite
dfHalfCircles |> @vlplot(:circle, x=:x1, y=:x2, color = :"y:n") 
```




![svg](img/output_11_0.svg)



## Regression data generated from noisy linear models

```@docs
make_regression
```


```@example foggy
using MLJ, DataFrames
X, y = make_regression(100, 5; noise=0.5, sparse=0.2, outliers=0.1)
dfRegression = DataFrame(X)
dfRegression.y = y
first(dfRegression, 3)
```
	
# Evaluating Model Performance

MLJ allows quick evaluation of a supervised model's performance
against a battery of selected losses or scores.
For more on available performance measures, see
[Performance Measures](performance_measures.md).

In addition to hold-out and cross-validation, the user can specify
their own list of train/test pairs of row indices for resampling, or
define their own re-usable resampling strategies.

For simultaneously evaluating *multiple* models and/or data
sets, see [Benchmarking](benchmarking.md).

## Evaluating against a single measure

```@setup evaluation_of_supervised_models
using MLJ
MLJ.color_off()
```

```@repl evaluation_of_supervised_models
using MLJ
X = (a=rand(12), b=rand(12), c=rand(12));
y = X.a + 2X.b + 0.05*rand(12);
model = (@load RidgeRegressor pkg=MultivariateStats verbosity=0)()
cv=CV(nfolds=3)
evaluate(model, X, y, resampling=cv, measure=l2, verbosity=0)
```

Alternatively, instead of applying `evaluate` to a model + data, one
may call `evaluate!` on an existing machine wrapping the model in
data:

```@repl evaluation_of_supervised_models
mach = machine(model, X, y)
evaluate!(mach, resampling=cv, measure=l2, verbosity=0)
```

(The latter call is a mutating call as the learned parameters stored in the
machine potentially change. )

## Multiple measures

```@repl evaluation_of_supervised_models
evaluate!(mach,
          resampling=cv,
          measure=[l1, rms, rmslp1], verbosity=0)
```

## Custom measures and weighted measures

```@repl evaluation_of_supervised_models
my_loss(yhat, y) = maximum((yhat - y).^2);

my_per_observation_loss(yhat, y) = abs.(yhat - y);
MLJ.reports_each_observation(::typeof(my_per_observation_loss)) = true;

my_weighted_score(yhat, y) = 1/mean(abs.(yhat - y));
my_weighted_score(yhat, y, w) = 1/mean(abs.((yhat - y).^w));
MLJ.supports_weights(::typeof(my_weighted_score)) = true;
MLJ.orientation(::typeof(my_weighted_score)) = :score;

holdout = Holdout(fraction_train=0.8)
weights = [1, 1, 2, 1, 1, 2, 3, 1, 1, 2, 3, 1];
evaluate!(mach,
          resampling=CV(nfolds=3),
          measure=[my_loss, my_per_observation_loss, my_weighted_score, l1],
          weights=weights, verbosity=0)
```

## User-specified train/test sets

Users can either provide their own list of train/test pairs of row indices for resampling, as in this example:

```@repl evaluation_of_supervised_models
fold1 = 1:6; fold2 = 7:12;
evaluate!(mach,
          resampling = [(fold1, fold2), (fold2, fold1)],
          measure=[l1, l2], verbosity=0)
```

Or define their own re-usable `ResamplingStrategy` objects, - see
[Custom resampling strategies](@ref) below.


## Built-in resampling strategies


```@docs
MLJBase.Holdout
```

```@docs
MLJBase.CV
```

```@docs
MLJBase.StratifiedCV
```

```@docs
MLJBase.TimeSeriesCV
```

## Custom resampling strategies

To define your own resampling strategy, make relevant parameters of
your strategy the fields of a new type `MyResamplingStrategy <:
MLJ.ResamplingStrategy`, and implement one of the following methods:

```julia
MLJ.train_test_pairs(my_strategy::MyResamplingStrategy, rows)
MLJ.train_test_pairs(my_strategy::MyResamplingStrategy, rows, y)
MLJ.train_test_pairs(my_strategy::MyResamplingStrategy, rows, X, y)
```

Each method takes a vector of indices `rows` and return a
vector `[(t1, e1), (t2, e2), ... (tk, ek)]` of train/test pairs of row
indices selected from `rows`. Here `X`, `y` are the input and target
data (ignored in simple strategies, such as `Holdout` and `CV`).

Here is the code for the `Holdout` strategy as an example:

```julia
struct Holdout <: ResamplingStrategy
    fraction_train::Float64
    shuffle::Bool
    rng::Union{Int,AbstractRNG}

    function Holdout(fraction_train, shuffle, rng)
        0 < fraction_train < 1 ||
            error("`fraction_train` must be between 0 and 1.")
        return new(fraction_train, shuffle, rng)
    end
end

# Keyword Constructor
function Holdout(; fraction_train::Float64=0.7, shuffle=nothing, rng=nothing)
    if rng isa Integer
        rng = MersenneTwister(rng)
    end
    if shuffle === nothing
        shuffle = ifelse(rng===nothing, false, true)
    end
    if rng === nothing
        rng = Random.GLOBAL_RNG
    end
    return Holdout(fraction_train, shuffle, rng)
end

function train_test_pairs(holdout::Holdout, rows)
    train, test = partition(rows, holdout.fraction_train,
                          shuffle=holdout.shuffle, rng=holdout.rng)
    return [(train, test),]
end
```

## API

```@docs
MLJBase.evaluate!
MLJBase.evaluate
MLJBase.PerformanceEvaluation
```
# Glossary

Note: This glossary includes some detail intended mainly for MLJ developers.

## Basics

### hyper-parameters

Parameters on which some learning algorithm depends, specified before
the algorithm is applied, and where learning is interpreted in the
broadest sense. For example, PCA feature reduction is a
"preprocessing" transformation "learning" a projection from training
data, governed by a dimension hyperparameter. Hyper-Parameters in our
sense may specify configuration (eg, number of parallel processes)
even when this does not effect the end-product of learning. (But we
exclude verbosity level.)

### model (object of abstract type `Model`)

Object collecting together hyperameters of a single algorithm.  Models
are classified either as *supervised* or *unsupervised* models (eg,
"transformers"), with corresponding subtypes `Supervised <: Model` and
`Unsupervised <: Model`.


### fit-result (type generally defined outside of MLJ)

Also known as "learned" or "fitted" parameters, these are "weights",
"coefficients", or similar paramaters learned by an algorithm, after
adopting the prescribed hyper-parameters. For example, decision trees
of a random forest, the coefficients and intercept of a linear model,
or the rotation and projection matrices of PCA reduction scheme.


### operation

Data-manipulating operations (methods) parameterized by some
fit-result. For supervised learners, the `predict`, `predict_mean`,
`predict_median`, or `predict_mode` methods; for transformers, the
`transform` or `inverse_transform` method. An operation may also
refer to an ordinary data-manipulating method that does *not* depend
on a fit-result (e.g., a broadcasted logarithm) which is then called
*static* operation for clarity. An operation that is not static is
*dynamic*.

### machine (object of type `Machine`)

An object consisting of:

(1) A model

(2) A fit-result (undefined until training)

(3) *Training arguments* (one for each data argument of the model's
associated `fit` method). A training argument is data used for
training (subsampled by specifying `rows=...` in `fit!`) but also in
evaluation (subsampled by specifying `rows=...` in `predict`,
`predict_mean`, etc). Generally, there are two training arguments for
supervised models, and just one for unsuperivsed models. Each argument
is either a `Source` node, wrapping concrete data supplied to the
`machine` constructor, or a `Node`, in the case of a learning network
(see below). Both kinds of nodes can be *called* with an optional
`rows=...` keyword argument to (lazily) return concrete data.

In addition, machines store "report" metadata, for recording
algorithm-specific statistics of training (eg, internal estimate of
generalization error, feature importances); and they cache information
allowing the fit-result to be updated without repeating unnecessary
information.

Machines are trained by calls to a `fit!` method which may be
passed an optional argument specifying the rows of data to be used in
training.

For more, see the [Machines](@ref) section.


## Learning Networks and Composite Models

*Note:* Multiple machines in a learning network may share the same
model, and multiple learning nodes may share the same machine.

### source node (object of type `Source`)

A container for training data and point of entry for new data in a
learning network (see below).


###  node (object of type `Node`)

Essentially a machine (whose arguments are possibly other nodes)
wrapped in an associated operation (e.g., `predict` or
`inverse_transform`). It consists primarily of:

1. An operation, static or dynamic.
1. A machine, or `nothing` if the operation is static.
1. Upstream connections to other nodes, specified by a list of
   *arguments* (one for each argument of the operation). These are the
   arguments on which the operation "acts" when the node `N` is
   called, as in `N()`.



### learning network

An acyclic directed graph implicit in the connections of a collection
of source(s) and nodes. 


### wrapper

Any model with one or more other models as hyper-parameters.


### composite model

Any wrapper, or any learning network, "exported" as a model (see
[Composing Models](composing_models.md)).
# Common MLJ Workflows

## Data ingestion

```@setup workflows
# to avoid RDatasets as a doc dependency:
using MLJ; color_off()
import DataFrames
channing = (Sex = rand(["Male","Female"], 462),
            Entry = rand(Int, 462),
            Exit = rand(Int, 462),
            Time = rand(Int, 462),
            Cens = rand(Int, 462)) |> DataFrames.DataFrame
coerce!(channing, :Sex => Multiclass)
```

```julia
import RDatasets
channing = RDatasets.dataset("boot", "channing")

julia> first(channing, 4)
4×5 DataFrame
 Row │ Sex   Entry  Exit   Time   Cens
     │ Cat…  Int32  Int32  Int32  Int32
─────┼──────────────────────────────────
   1 │ Male    782    909    127      1
   2 │ Male   1020   1128    108      1
   3 │ Male    856    969    113      1
   4 │ Male    915    957     42      1
```

Inspecting metadata, including column scientific types:

```@example workflows
schema(channing)
```

Horizontally splitting data and shuffling rows.

Here `y` is the `:Exit` column and `X` everything else:

```@example workflows
y, X =  unpack(channing, ==(:Exit), rng=123);
nothing # hide
```

Here `y` is the `:Exit` column and `X` everything else except `:Time`:

```@example workflows
y, X =  unpack(channing,
               ==(:Exit),
               !=(:Time);
               rng=123);
scitype(y)
```

```@example workflows
schema(X)
```

Fixing wrong scientfic types in `X`:

```@example workflows
X = coerce(X, :Exit=>Continuous, :Entry=>Continuous, :Cens=>Multiclass)
schema(X)
```


Loading a built-in supervised dataset:

```@example workflows
table = load_iris();
schema(table)
```

Loading a built-in data set already split into `X` and `y`:

```@example workflows
X, y = @load_iris;
selectrows(X, 1:4) # selectrows works for any Tables.jl table
```

```@example workflows
y[1:4]
```

Splitting data vertically after row shuffling:

```@example workflows
channing_train, channing_test = partition(channing, 0.6, rng=123);
nothing # hide
```

Or, if already horizontally split:

```@example workflows
(Xtrain, Xtest), (ytrain, ytest)  = partition((X, y), 0.6, multi=true,  rng=123)
```


## Model search

*Reference:*   [Model Search](model_search.md)

Searching for a supervised model:

```@example workflows
X, y = @load_boston
ms = models(matching(X, y))
```

```@example workflows
ms[6]
```

```@example workflows
models("Tree");
```

A more refined search:

```@example workflows
models() do model
    matching(model, X, y) &&
    model.prediction_type == :deterministic &&
    model.is_pure_julia
end;
nothing # hide
```

Searching for an unsupervised model:

```@example workflows
models(matching(X))
```

Getting the metadata entry for a given model type:

```@example workflows
info("PCA")
info("RidgeRegressor", pkg="MultivariateStats") # a model type in multiple packages
```

## Instantiating a model

*Reference:*   [Getting Started](@ref), [Loading Model Code](@ref)

```@example workflows
Tree = @load DecisionTreeClassifier pkg=DecisionTree
tree = Tree(min_samples_split=5, max_depth=4)
```

or

```@julia
tree = (@load DecisionTreeClassifier)()
tree.min_samples_split = 5
tree.max_depth = 4
```

## Evaluating a model

*Reference:*   [Evaluating Model Performance](evaluating_model_performance.md)


```@example workflows
X, y = @load_boston
KNN = @load KNNRegressor
knn = KNN()
evaluate(knn, X, y,
         resampling=CV(nfolds=5),
         measure=[RootMeanSquaredError(), MeanAbsoluteError()])
```

Note `RootMeanSquaredError()` has alias `rms` and `MeanAbsoluteError()` has alias `mae`.

Do `measures()` to list all losses and scores and their aliases.


##  Basic fit/evaluate/predict by hand:

*Reference:*   [Getting Started](index.md), [Machines](machines.md),
[Evaluating Model Performance](evaluating_model_performance.md), [Performance Measures](performance_measures.md)

```@example workflows
crabs = load_crabs() |> DataFrames.DataFrame
schema(crabs)
```

```@example workflows
y, X = unpack(crabs, ==(:sp), !in([:index, :sex]); rng=123)


Tree = @load DecisionTreeClassifier pkg=DecisionTree
tree = Tree(max_depth=2) # hide
```

Bind the model and data together in a *machine* , which will
additionally store the learned parameters (*fitresults*) when fit:

```@example workflows
mach = machine(tree, X, y)
```

Split row indices into training and evaluation rows:

```@example workflows
train, test = partition(eachindex(y), 0.7); # 70:30 split
```

Fit on train and evaluate on test:

```@example workflows
fit!(mach, rows=train)
yhat = predict(mach, X[test,:])
mean(LogLoss(tol=1e-4)(yhat, y[test]))
```

Note `LogLoss()` has aliases `log_loss` and `cross_entropy`.

Run `measures()` to list all losses and scores and their aliases ("instances").

Predict on new data:

```@example workflows
Xnew = (FL = rand(3), RW = rand(3), CL = rand(3), CW = rand(3), BD =rand(3))
predict(mach, Xnew)      # a vector of distributions
```

```@example workflows
predict_mode(mach, Xnew) # a vector of point-predictions
```

## More performance evaluation examples

Evaluating model + data directly:

```@example workflows
evaluate(tree, X, y,
         resampling=Holdout(fraction_train=0.7, shuffle=true, rng=1234),
         measure=[LogLoss(), Accuracy()])
```

If a machine is already defined, as above:

```@example workflows
evaluate!(mach,
          resampling=Holdout(fraction_train=0.7, shuffle=true, rng=1234),
          measure=[LogLoss(), Accuracy()])
```

Using cross-validation:

```@example workflows
evaluate!(mach, resampling=CV(nfolds=5, shuffle=true, rng=1234),
          measure=[LogLoss(), Accuracy()])
```

With user-specified train/test pairs of row indices:

```@example workflows
f1, f2, f3 = 1:13, 14:26, 27:36
pairs = [(f1, vcat(f2, f3)), (f2, vcat(f3, f1)), (f3, vcat(f1, f2))];
evaluate!(mach,
          resampling=pairs,
          measure=[LogLoss(), Accuracy()])
```

Changing a hyperparameter and re-evaluating:

```@example workflows
tree.max_depth = 3
evaluate!(mach,
          resampling=CV(nfolds=5, shuffle=true, rng=1234),
          measure=[LogLoss(), Accuracy()])
```

##  Inspecting training results

Fit a ordinary least square model to some synthetic data:

```@example workflows
x1 = rand(100)
x2 = rand(100)

X = (x1=x1, x2=x2)
y = x1 - 2x2 + 0.1*rand(100);

OLS = @load LinearRegressor pkg=GLM
ols = OLS()
mach =  machine(ols, X, y) |> fit!
```

Get a named tuple representing the learned parameters,
human-readable if appropriate:

```@example workflows
fitted_params(mach)
```

Get other training-related information:

```@example workflows
report(mach)
```

##  Basic fit/transform for unsupervised models

Load data:

```@example workflows
X, y = @load_iris
train, test = partition(eachindex(y), 0.97, shuffle=true, rng=123)
```

Instantiate and fit the model/machine:

```@example workflows
PCA = @load PCA
pca = PCA(maxoutdim=2)
mach = machine(pca, X)
fit!(mach, rows=train)
```

Transform selected data bound to the machine:

```@example workflows
transform(mach, rows=test);
```

Transform new data:

```@example workflows
Xnew = (sepal_length=rand(3), sepal_width=rand(3),
        petal_length=rand(3), petal_width=rand(3));
transform(mach, Xnew)
```

##  Inverting learned transformations

```@example workflows
y = rand(100);
stand = Standardizer()
mach = machine(stand, y)
fit!(mach)
z = transform(mach, y);
@assert inverse_transform(mach, z) ≈ y # true
```

## Nested hyperparameter tuning

*Reference:*   [Tuning Models](tuning_models.md)

```@example workflows
X, y = @load_iris; nothing # hide
```

Define a model with nested hyperparameters:

```@example workflows
Tree = @load DecisionTreeClassifier pkg=DecisionTree
tree = Tree()
forest = EnsembleModel(model=tree, n=300)
```

Define ranges for hyperparameters to be tuned:

```@example workflows
r1 = range(forest, :bagging_fraction, lower=0.5, upper=1.0, scale=:log10)
```

```@example workflows
r2 = range(forest, :(model.n_subfeatures), lower=1, upper=4) # nested
```

Wrap the model in a tuning strategy:

```@example workflows
tuned_forest = TunedModel(model=forest,
                          tuning=Grid(resolution=12),
                          resampling=CV(nfolds=6),
                          ranges=[r1, r2],
                          measure=BrierLoss())
```

Bound the wrapped model to data:

```@example workflows
mach = machine(tuned_forest, X, y)
```

Fitting the resultant machine optimizes the hyperparameters specified
in `range`, using the specified `tuning` and `resampling` strategies
and performance `measure` (possibly a vector of measures), and
retrains on all data bound to the machine:

```@example workflows
fit!(mach)
```

Inspecting the optimal model:

```@example workflows
F = fitted_params(mach)
```

```@example workflows
F.best_model
```

Inspecting details of tuning procedure:

```@example workflows
r = report(mach);
keys(r)
```

```@example workflows
r.history[[1,end]]
```

Visualizing these results:

```julia
using Plots
plot(mach)
```

![](img/workflows_tuning_plot.png)

Predicting on new data using the optimized model:

```@example workflows
predict(mach, Xnew)
```

## Constructing linear pipelines

*Reference:*   [Composing Models](composing_models.md)

Constructing a linear (unbranching) pipeline with a *learned* target
transformation/inverse transformation:

```@example workflows
X, y = @load_reduced_ames
KNN = @load KNNRegressor
knn_with_target = TransformedTargetModel(model=KNN(K=3), target=Standardizer())
pipe = (X -> coerce(X, :age=>Continuous)) |> OneHotEncoder() |> knn_with_target
```

Evaluating the pipeline (just as you would any other model):

```@example workflows
pipe.one_hot_encoder.drop_last = true
evaluate(pipe, X, y, resampling=Holdout(), measure=RootMeanSquaredError(), verbosity=2)
```

Inspecting the learned parameters in a pipeline:

```@example workflows
mach = machine(pipe, X, y) |> fit!
F = fitted_params(mach)
F.transformed_target_model_deterministic.model
```

Constructing a linear (unbranching) pipeline with a *static* (unlearned)
target transformation/inverse transformation:

```@example workflows
Tree = @load DecisionTreeRegressor pkg=DecisionTree verbosity=0
tree_with_target = TransformedTargetModel(model=Tree(),
                                          target=y -> log.(y),
                                          inverse = z -> exp.(z))
pipe2 = (X -> coerce(X, :age=>Continuous)) |> OneHotEncoder() |> tree_with_target;
nothing # hide
```

## Creating a homogeneous ensemble of models

*Reference:* [Homogeneous Ensembles](homogeneous_ensembles.md)

```@example workflows
X, y = @load_iris
Tree = @load DecisionTreeClassifier pkg=DecisionTree
tree = Tree()
forest = EnsembleModel(model=tree, bagging_fraction=0.8, n=300)
mach = machine(forest, X, y)
evaluate!(mach, measure=LogLoss())
```

## Performance curves

Generate a plot of performance, as a function of some hyperparameter
(building on the preceding example)

Single performance curve:

```@example workflows
r = range(forest, :n, lower=1, upper=1000, scale=:log10)
curve = learning_curve(mach,
                       range=r,
                       resampling=Holdout(),
                       resolution=50,
                       measure=LogLoss(),
                       verbosity=0)
```

```julia
using Plots
plot(curve.parameter_values, curve.measurements, xlab=curve.parameter_name, xscale=curve.parameter_scale)
```

![](img/workflows_learning_curve.png)

Multiple curves:

```@example workflows
curve = learning_curve(mach,
                       range=r,
                       resampling=Holdout(),
                       measure=LogLoss(),
                       resolution=50,
                       rng_name=:rng,
                       rngs=4,
                       verbosity=0)
```

```julia
plot(curve.parameter_values, curve.measurements,
xlab=curve.parameter_name, xscale=curve.parameter_scale)
```

![](img/workflows_learning_curves.png)
# Acceleration and Parallelism

!!! warning "Experimental API"

    The acceleration API is experimental and may not work correctly in all
    cases, especially if trying to use an acceleration method that your
    version of Julia or installed packages cannot support. The API is also
    subject to breaking changes during minor or major releases without
    warning.

## User-facing interface

To enable composable, extensible acceleration of core MLJ methods,
[ComputationalResources.jl](https://github.com/timholy/ComputationalResources.jl)
is utilized to provide some basic types and functions to make implementing
acceleration easy. However, ambitious users or package authors have the option
to define their own types to be passed as resources to `acceleration`, which
must be `<:ComputationalResources.AbstractResource`.

Methods which support some form of acceleration support the `acceleration`
keyword argument, which can be passed a "resource" from
`ComputationalResources`. For example, passing `acceleration=CPUProcesses()`
will utilize `Distributed`'s multiprocessing functionality to accelerate the
computation, while `acceleration=CPUThreads()` will use Julia's PARTR
threading model to perform acceleration.

The default computational resource is `CPU1()`, which is simply serial
processing via CPU. The default resource can be changed as in this
example: `MLJ.default_resource(CPUProcesses())`. The argument must
always have type `<:ComputationalResource.AbstractResource`. To
inspect the current default, use `MLJ.default_resource()`.

!!! note

    The `CPUThreads()` resource is only available when running a version of
    Julia with `Threads.@spawn` available.

!!! note

    You cannot use `CPUThreads()` with models wrapping python code.
# Benchmarking

This feature not yet available.

[CONTRIBUTE.md](https://github.com/alan-turing-institute/MLJ.jl/blob/master/CONTRIBUTE.md)
# Machines

Recall from [Getting Started](@ref) that a machine binds a model
(i.e., a choice of algorithm + hyperparameters) to data (see more at
[Constructing machines](@ref) below). A machine is also the object
storing *learned* parameters.  Under the hood, calling `fit!` on a
machine calls either `MLJBase.fit` or `MLJBase.update`, depending on
the machine's internal state (as recorded in private fields
`old_model` and `old_rows`). These lower-level `fit` and `update`
methods, which are not ordinarily called directly by the user,
dispatch on the model and a view of the data defined by the optional
`rows` keyword argument of `fit!` (all rows by default). 

# Warm restarts

If a model `update` method has been implemented for the model, calls
to `fit!` will avoid redundant calculations for certain kinds of model
mutations. The main use-case is increasing an iteration parameter,
such as the number of epochs in a neural network. To test if
`SomeIterativeModel` supports this feature, check
`iteration_parameter(SomeIterativeModel)` is different from `nothing`.

```@example machines
using MLJ; color_off() # hide
tree = (@load DecisionTreeClassifier pkg=DecisionTree verbosity=0)()
forest = EnsembleModel(model=tree, n=10);
X, y = @load_iris;
mach = machine(forest, X, y)
fit!(mach, verbosity=2);
```

Generally, changing a hyperparameter triggers retraining on calls to
subsequent `fit!`:

```@repl machines
forest.bagging_fraction=0.5
fit!(mach, verbosity=2);
```

However, for this iterative model, increasing the iteration parameter
only adds models to the existing ensemble:

```@repl machines
forest.n=15
fit!(mach, verbosity=2);
```

Call `fit!` again without making a change and no retraining occurs:

```@repl machines
fit!(mach);
```

However, retraining can be forced:

```@repl machines
fit!(mach, force=true);
```

And is re-triggered if the view of the data changes:

```@repl machines
fit!(mach, rows=1:100);
```

```@repl machines
fit!(mach, rows=1:100);
```

If an iterative model exposes it's iteration parameter as a
hyper-parameter, and it implements the warm restart behaviour above,
then it can be wrapped in a "control strategy", like an early stopping
critetion. See [Controlling Iterative Models](@ref) for details.


## Inspecting machines

There are two methods for inspecting the outcomes of training in
MLJ. To obtain a named-tuple describing the learned parameters (in a
user-friendly way where possible) use `fitted_params(mach)`. All other
training-related outcomes are inspected with `report(mach)`.

```@example machines
X, y = @load_iris
pca = (@load PCA verbosity=0)()
mach = machine(pca, X)
fit!(mach)
```

```@repl machines
fitted_params(mach)
report(mach)
```

```@docs
fitted_params
report
```


## Constructing machines

A machine is constructed with the syntax `machine(model, args...)`
where the possibilities for `args` (called *training arguments*) are
summarized in table below. Here `X` and `y` represent inputs and
target, respectively, and `Xout` the output of a `transform` call.
Machines for supervised models may have additional training arguments,
such as a vector of per-observation weights (in which case
`supports_weights(model) == true`).

`model` supertype   | `machine` constructor calls | operation calls (first compulsory)
--------------------|-----------------------------|--------------------------------------
`Deterministic <: Supervised`    | `machine(model, X, y, extras...)` | `predict(mach, Xnew)`, `transform(mach, Xnew)`, `inverse_transform(mach, Xout)`
`Probabilistic <: Supervised`    | `machine(model, X, y, extras...)` | `predict(mach, Xnew)`, `predict_mean(mach, Xnew)`, `predict_median(mach, Xnew)`, `predict_mode(mach, Xnew)`, `transform(mach, Xnew)`, `inverse_transform(mach, Xout)`
`Unsupervised` (except `Static`) | `machine(model, X)` | `transform(mach, Xnew)`, `inverse_transform(mach, Xout)`, `predict(mach, Xnew)`
`Static`                        | `machine(model)`    | `transform(mach, Xnews...)`, `inverse_transform(mach, Xout)`

All operations on machines (`predict`, `transform`, etc) have exactly
one argument (`Xnew` or `Xout` above) after `mach`, the machine
instance. An exception is a machine bound to a `Static` model, which
can have any number of arguments after `mach`. For more on `Static`
transformers (which have no *training* arguments) see [Static
transformers](@ref).

A machine is reconstructed from a file using the syntax
`machine("my_machine.jlso")`, or `machine("my_machine.jlso", args...)`
if retraining using new data. See [Saving machines](@ref) below.


## Lowering memory demands

For large data sets you may be able to save memory by suppressing data
caching that some models perform to increase speed. To do this,
specify `cache=false`, as in

```julia
machine(model, X, y, cache=false)
```


### Constructing machines in learning networks

Instead of data `X`, `y`, etc,  the `machine` constructor is provided
`Node` or `Source` objects ("dynamic data") when building a learning
network. See [Composing Models](composing_models.md) for more on this
advanced feature. One also uses `machine` to wrap a machine
around a whole learning network; see [Learning network
machines](@ref).


## Saving machines

To save a machine to file, use the `MLJ.save` command:

```julia
tree = (@load DecisionTreeClassifier pkg=DecisionTree verbosity=0)()
mach = fit!(machine(tree, X, y))
MLJ.save("my_machine.jlso", mach)
```

To de-serialize, one uses the `machine` constructor:

```julia
mach2 = machine("my_machine.jlso")
predict(mach2, Xnew);
```

The machine `mach2` cannot be retrained; however, by providing data to
the constructor one can enable retraining using the saved model
hyperparameters (which overwrites the saved learned parameters):

```julia
mach3 = machine("my_machine.jlso", Xnew, ynew)
fit!(mach3)
```


## Internals

For a supervised machine the `predict` method calls a lower-level
`MLJBase.predict` method, dispatched on the underlying model and the
`fitresult` (see below). To see `predict` in action, as well as its
unsupervised cousins `transform` and `inverse_transform`, see
[Getting Started](index.md).

The fields of a `Machine` instance (which should not generally be
accessed by the user) are:

- `model` - the struct containing the hyperparameters to be used in
  calls to `fit!`

- `fitresult` - the learned parameters in a raw form, initially undefined

- `args` - a tuple of the data, each element wrapped in a source node;
  see [Learning Networks](@ref) (in the supervised learning example
  above, `args = (source(X), source(y))`)

- `report` - outputs of training not encoded in `fitresult` (eg, feature rankings)

- `old_model` - a deep copy of the model used in the last call to `fit!`

- `old_rows` -  a copy of the row indices used in last call to `fit!`

- `cache`

The interested reader can learn more on machine internals by examining
the simplified code excerpt in [Internals](internals.md).


## API Reference

```@docs
MLJBase.machine
fit!
fit_only!
MLJSerialization.save
```
# Weights

In machine learning it is possible to assign each observation an
independent significance, or *weight*, either in **training** or in
**performance evaluation**, or both. 

There are two kinds of weights in use in MLJ:

- *per observation weights* (also just called *weights*) refer to
  weight vectors of the same length as the number of observations

- *class weights* refer to dictionaries keyed on the target classes
  (levels) for use in classification problems
  

## Specifying weights in training

To specify weights in training you bind the weights to the model along
with the data when constructing a machine.  For supervised models the
weights are specified last:

```julia
KNNRegressor = @load KNNRegressor
model = KNNRegressor()
X, y = make_regression(10, 3)
w = rand(length(y))

mach = machine(model, X, y, w) |> fit!
```

Note that `model` supports per observation weights if
`supports_weights(model)` is `true`. To list all such models, do

```julia
models() do m
    m.supports_weights
end
```

The model `model` supports class weights if
`supports_class_weights(model)` is `true`.


## Specifying weights in performance evaluation

When calling an MLJ measure (metric) that supports weights, provide the
weights as the last argument, as in

```julia
_, y = @load_iris
ŷ = shuffle(y)
w = Dict("versicolor" => 1, "setosa" => 2, "virginica"=> 3)
macro_f1score(ŷ, y, w)
```

You can use `supports_weights` and `supports_class_weights` on
measures to check weight support. For example, to list all measures
supporting per observation weights, do

```julia 
measures() do m 
   m.supports_weights 
end 
```

See also [Evaluating Model Performance](@ref).

To pass weights to all the measures listed in an `evaluate!/evaluate`
call, use the keyword specifiers `weights=...` or
`class_weights=...`. For details, see [`evaluate!`](@ref).
# Loading Model Code

Once the name of a model, and the package providing that model, have
been identified (see [Model Search](@ref model_search)) one can either
import the model type interactively with `@iload`, as shown under
[Installation](@ref), or use `@load` as shown below. The `@load` macro
works from within a module, a package or a function, provided the
relevant package providing the MLJ interface has been added to your
package environment.

In general, the code providing core functionality for the model
(living in a packge you should consult for documentation) may be
different from the package providing the MLJ interface. Since the core
package is a dependency of the interface package, only the interface
package needs to be added to your environment.

For instance, suppose you have activated a Julia package environment
`my_env` that you wish to use for your MLJ project; for example you
have run:


```julia
using Pkg
Pkg.activate("my_env", shared=true)
```

And, furthermore, suppose you want to use `DecisionTreeClassifier`,
provided by the
[DecisionTree.jl](https://github.com/bensadeghi/DecisionTree.jl)
package. Then, to determine which package provides the MLJ interface
you call `load_path`:

```julia
julia> load_path("DecisionTreeClassifier", pkg="DecisionTree")
"MLJDecisionTreeInterface.DecisionTreeClassifier"
```

In this case we see that the package required is
MLJDecisionTreeInterface.jl. If this package is not in `my_env` (do
`Pkg.status()` to check) you add it by running

```julia
julia> Pkg.add("MLJDecisionTreeInterface");
```

So long as `my_env` is the active environment, this action need never
be repeated (unless you run `Pkg.rm("MLJDecisionTreeInterface")`). You
are now ready to instantiate a decision tree classifier:

```julia
julia> Tree = @load DecisionTree pkg=DecisionTree
julia> tree = Tree()
```

which is equivalent to

```julia
julia> import MLJDecisionTreeInterface.DecisionTreeClassifier
julia> Tree = MLJDecisionTreeInterface.DecisionTreeClassifier()
julia> tree = Tree()
```

*Tip.* The specification `pkg=...` above can be dropped for the many
models that are provided by only a single package.


## API

```@docs
load_path
@load
@iload
```
# [Model Search](@id model_search)

MLJ has a model registry, allowing the user to search models and their
properties, without loading all the packages containing model code. In
turn, this allows one to efficiently find all models solving a given
machine learning task. The task itself is specified with the help of
the `matching` method, and the search executed with the `models`
methods, as detailed below. 

For commonly encountered problems with model search, see also
[Preparing Data](@ref).

A table of all models is also given at [List of Supported
Models](@ref model_list).


## Model metadata

*Terminology.* In this section the word "model" refers to a metadata
entry in the model registry, as opposed to an actual model `struct`
that such an entry represents. One can obtain such an entry with the
`info` command:

```@setup tokai
using MLJ
MLJ.color_off()
```

```@repl tokai
info("PCA")
```

So a "model" in the present context is just a named tuple containing
metadata, and not an actual model type or instance. If two models with
the same name occur in different packages, the package name must be
specified, as in `info("LinearRegressor", pkg="GLM")`.


## General model queries

We list all models (named tuples) using `models()`, and list the models for which code is  already loaded with `localmodels()`:

```@repl tokai
localmodels()
localmodels()[2]
```

One can search for models containing specified strings or regular expressions in their `docstring` attributes, as in

```@repl tokai 
models("forest")
```

or by specifying a filter (`Bool`-valued function):

```@repl tokai
filter(model) = model.is_supervised &&
                model.input_scitype >: MLJ.Table(Continuous) &&
                model.target_scitype >: AbstractVector{<:Multiclass{3}} &&
                model.prediction_type == :deterministic
models(filter)
```

Multiple test arguments may be passed to `models`, which are applied
conjunctively.


## Matching models to data

Common searches are streamlined with the help of the `matching`
command, defined as follows:

- `matching(model, X, y) == true` exactly when `model` is supervised
   and admits inputs and targets with the scientific types of `X` and
   `y`, respectively

- `matching(model, X) == true` exactly when `model` is unsupervised
   and admits inputs with the scientific types of `X`.

So, to search for all supervised probabilistic models handling input
`X` and target `y`, one can define the testing function `task` by

```julia
task(model) = matching(model, X, y) && model.prediction_type == :probabilistic
```

And execute the search with

```julia
models(task)
```

Also defined are `Bool`-valued callable objects `matching(model)`,
`matching(X, y)` and `matching(X)`, with obvious behaviour. For example,
`matching(X, y)(model) = matching(model, X, y)`.

So, to search for all models compatible with input `X` and target `y`,
for example, one executes

```julia
models(matching(X, y))
```

while the preceding search can also be written

```julia
models() do model
    matching(model, X, y) &&
    model.prediction_type == :probabilistic
end
```

## API

```@docs
models
localmodels
```
# Quick-Start Guide to Adding Models

The definitive specification of the MLJ model interface is given in
[Adding Models for General Use](@ref). In the more informal and
condensed instructions below, we assume: (i) you have a Julia
**registered** package `YourPackage.jl` implementing some machine
learning models; (ii) that you would like to interface and register
these models with MLJ; and (iii) that you have a rough understanding
of how things work with MLJ.  In particular you are familiar with:

- what [scientific types](https://github.com/JuliaAI/ScientificTypes.jl) are

- what `Probabilistic`, `Deterministic` and `Unsupervised` models are

- the fact that MLJ generally works with tables rather than bare bone
  matrices. Here a *table* is a container satisfying the
  [Tables.jl](https://github.com/JuliaData/Tables.jl) API (e.g.,
  DataFrame, JuliaDB table, CSV file, named tuple of equi-length
  vectors)

- [CategoricalArrays.jl](https://github.com/JuliaData/CategoricalArrays.jl)
  (if working with finite discrete data, e.g., doing classification)

If you're not familiar with any one of these points, please refer to
relevant sections of this manual, and in particular [Getting
Started](@ref) and [Adding Models for General Use](@ref).

*But tables don't make sense for my model!* If a case can be made that
tabular input does not make sense for your particular model, then MLJ can
still handle this; you just need to define a non-tabular
`input_scitype` trait. However, you should probably open an issue to
clarify the appropriate declaration. The discussion below assumes
input data is tabular.

For simplicity, this document assumes no data front-end is to be
defined for your model. Adding a data front-end, which offers the MLJ
user some performances benefits, is easy to add post-facto, and is
described in [Implementing a data front-end](@ref).

### Overview

To write an interface create a file or a module in your package which
includes:

- a `using MLJModelInterface` or `import MLJModelInterface` statement

- MLJ-compatible model types and constructors,

- implementation of `fit`, `predict`/`transform` and optionally
  `fitted_params` for your models,

- metadata for your package and for each of your models

#### Important

[MLJModelInterface](https://github.com/JuliaAI/MLJModelInterface.jl)
is a very light-weight interface allowing you to *define* your
interface, but does not provide the functionality required to use or
test your interface; this requires
[MLJBase](https://github.com/JuliaAI/MLJBase.jl). So,
while you only need to add `MLJModelInterface` to your project's
[deps], for testing purposes you need to add
[MLJBase](https://github.com/JuliaAI/MLJBase.jl) to your
project's [extras] and [targets]. In testing, simply use `MLJBase` in
place of `MLJModelInterface`.

We give some details for each step below with, each time, a few
examples that you can mimic.  The instructions are intentionally
brief.


### Model type and constructor

MLJ-compatible constructors for your models need to meet the following requirements:

* be `mutable struct`,
* be subtypes of `MLJModelInterface.Probabilistic` or `MLJModelInterface.Deterministic` or `MLJModelInterface.Unsupervised`,
* have fields corresponding exclusively to hyperparameters,
* have a keyword constructor assigning default values to all
  hyperparameters.

You may use the `@mlj_model` macro from `MLJModelInterface` to declare
a (non parametric) model type:

```julia
MLJModelInterface.@mlj_model mutable struct YourModel <: MLJModelInterface.Deterministic
    a::Float64 = 0.5::(_ > 0)
    b::String  = "svd"::(_ in ("svd","qr"))
end
```

That macro specifies:

* A keyword constructor (here `YourModel(; a=..., b=...)`),
* Default values for the hyperparameters,
* Constraints on the hyperparameters where `_` refers to a value passed.

Further to the last point, `a::Float64 = 0.5::(_ > 0)` indicates that
the field `a` is a `Float64`, takes `0.5` as default value, and
expects its value to be positive.

Please see [this
issue](https://github.com/JuliaAI/MLJBase.jl/issues/68)
for a known issue and workaround relating to the use of `@mlj_model`
with negative defaults.

If you decide **not** to use the `@mlj_model` macro (e.g. in the case
of a parametric type), you will need to write a keyword constructor
and a `clean!` method:

```julia
mutable struct YourModel <: MLJModelInterface.Deterministic
    a::Float64
end
function YourModel(; a=0.5)
    model   = YourModel(a)
    message = MLJModelInterface.clean!(model)
    isempty(message) || @warn message
    return model
end
function MLJModelInterface.clean!(m::YourModel)
    warning = ""
    if m.a <= 0
        warning *= "Parameter `a` expected to be positive, resetting to 0.5"
        m.a = 0.5
    end
    return warning
end
```

**Additional notes**:

- Please annotate all fields with concrete types, if possible, using
  type parameters if necessary.

- Please prefer `Symbol` over `String` if you can (e.g. to pass the name of a solver).

- Please add constraints to your fields even if they seem obvious to you.

- Your model may have 0 fields, that's fine.

- Although not essential, try to avoid Union types for model
  fields. For example, a field declaration `features::Vector{Symbol}`
  with a default of `Symbol[]` (detected with the `isempty` method) is
  preferred to `features::Union{Vector{Symbol}, Nothing}` with a default
  of `nothing`.


**Examples**:

- [KNNClassifier](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/NearestNeighbors.jl#L62-L69) which uses `@mlj_model`,
- [XGBoostRegressor](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/XGBoost.jl#L17-L161) which does not.


### Fit

The implementation of `fit` will look like

```julia
function MLJModelInterface.fit(m::YourModel, verbosity, X, y, w=nothing)
    # body ...
    return (fitresult, cache, report)
end
```

where `y` should only be there for a supervised model and `w` for a
supervised model that supports sample weights.  You **must** type
`verbosity` to `Int` and you **must not** type `X`, `y` and `w` (MLJ
handles that).

#### Regressor

In the body of the `fit` function, you should assume that `X` is a
table and that `y` is an `AbstractVector` (for multitask regression it
may be a table).

Typical steps in the body of the `fit` function will be:

* forming a matrix-view of the data, possibly transposed if your model
  expects a `p x n` formalism (MLJ assumes columns are features by
  default i.e. `n x p`), use `MLJModelInterface.matrix` for this,

* passing the data to your model,

* returning the results as a tuple `(fitresult, cache, report)`.

The `fitresult` part should contain everything that is needed at the
`predict` or `transform` step, it should not be expected to be
accessed by users.  The `cache` should be left to `nothing` for now.
The `report` should be a `NamedTuple` with any auxiliary useful
information that a user would want to know about the fit (e.g.,
feature rankings). See more on this below.

**Example**: GLM's [LinearRegressor](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L95-L105)


#### Classifier

For a classifier, the steps are fairly similar to a regressor with these differences:

1. `y` will be a categorical vector and you will typically want to use
   the integer encoding of `y` instead of `CategoricalValue`s; use
   `MLJModelInterface.int` for this.
1.  You will need to pass the full pool of target labels (not just
   those observed in the training data) and additionally, in the
   `Deterministic` case, the encoding, to make these available to
   `predict`. A simple way to do this is to pass `y[1]` in the
   `fitresult`, for then `MLJModelInterface.classes(y[1])` is a complete list of
   possible categorical elements, and `d = MLJModelInterface.decoder(y[1])` is a
   method for recovering categorical elements from their integer
   representations (e.g., `d(2)` is the categorical element with `2`
   as encoding).
2. In the case of a *probabilistic* classifier you should pass all
   probabilities simultaneously to the `UnivariateFinite` constructor
   to get an abstract `UnivariateFinite` vector (type
   `UnivariateFiniteArray`) rather than use comprehension or
   broadcasting to get a vanilla vector. This is for performance
   reasons.
   
If implementing a classifier, you should probably consult the more
detailed instructions at [The predict method](@ref).

**Examples**:

-  GLM's [BinaryClassifier](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L119-L131) (`Probabilistic`)

- LIBSVM's [SVC](https://github.com/JuliaAI/MLJModels.jl/blob/master/src/LIBSVM.jl) (`Deterministic`)


#### Transformer

Nothing special for a transformer.

**Example**: [FillImputer](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/builtins/Transformers.jl#L54-L64)


### Fitted parameters

There is a function you can optionally implement which will return the
learned parameters of your model for purposes of user-inspection. For
instance, in the case of a linear regression, the user may want to get
direct access to the coefficients and intercept. This should be as human and
machine readable as practical (not a graphical representation) and the
information should be combined in the form of a named tuple.

The function will always look like:

```julia
function MLJModelInterface.fitted_params(model::YourModel, fitresult)
    # extract what's relevant from `fitresult`
    # ...
    # then return as a NamedTuple
    return (learned_param1 = ..., learned_param2 = ...)
end
```

**Example**: for [GLM models](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L133-L137)


### Summary of user interface points (or, What to put where?)

Recall that the `fitresult` returned as part of `fit` represents
everything needed by `predict` (or `transform`) to make new
predictions. It is not intended to be directly inspected by the
user. Here is a summary of the interface points for users that your
implementation creates:

- Use `fitted_params` to expose *learned parameters*, such as linear
  coefficients, to the user in a machine and human readable form (for
  re-use in another model, for example).
- Use the fields of your model struct for *hyperparameters*, i.e.,
  those parameters declared by the user ahead of time that generally
  affect the outcome of training. It is okay to add "control"
  parameters (such a specifying an `acceleration` parameter specifying
  computational resources, as
  [here](https://github.com/alan-turing-institute/MLJ.jl/blob/master/src/ensembles.jl#L193)).
- Use `report` to return *everything else*, including model-specific
  *methods* (or other callable objects). This includes: feature rankings,
  decision boundaries, SVM support vectors, clustering centres,
  methods for visualizing training outcomes, methods for saving
  learned parameters in a custom format, degrees of freedom, deviance,
  etc. If there is a performance cost to extra functionality you want
  to expose, the functionality can be toggled on/off through a
  hyperparameter, but this should otherwise be avoided. For, example,
  in a decision tree model `report.print_tree(depth)` might generate
  a pretty tree representation of the learned tree, up to the
  specified `depth`.

### Predict/Transform

The implementation of `predict` (for a supervised model) or
`transform` (for an unsupervised one) will look like:

```julia
function MLJModelInterface.predict(m::YourModel, fitresult, Xnew)
    # ...
end
```

Here `Xnew` is expected to be a table and part of the logic in `predict` or `transform` may be similar to that in `fit`.

The values returned should be:

model subtype   | return value of predict/transform
----------------|----------------------------------
`Deterministic` | vector of values (or table if multi-target)
`Probabilistic` | vector of `Distribution` objects, for classifiers in particular, a vector of `UnivariateFinite`
`Unsupervised`  | table

In the case of a `Probabilistic` model, you may further want to
implement a `predict_mean` or a `predict_mode`. However,
MLJModelInterface provides fallbacks, defined in terms of `predict`,
whose performance may suffice.


**Examples**

- Deterministic regression: [KNNRegressor](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/NearestNeighbors.jl#L124-L145)
- Probabilistic regression: [LinearRegressor](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L154-L158) and the [`predict_mean`](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L144-L147)
- Probabilistic classification: [LogisticClassifier](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L165-L168)

### Metadata

Adding metadata for your model(s) is crucial for the discoverability
of your package and its models and to make sure your model is used
with data it can handle.  You can individually overload a number of
trait functions that encode this metadata by following the instuctions
in [Adding Models for General Use](@ref)), which also explains these
traits in more detail. However, your most convenient option is to use
`metadata_model` and `metadata_pkg` functionalities from
`MLJModelInterface` to do this:

```julia
const ALL_MODELS = Union{YourModel1, YourModel2, ...}

MLJModelInterface.metadata_pkg.(ALL_MODELS
    name = "YourPackage",
    uuid = "6ee0df7b-...", # see your Project.toml
    url  = "https://...",  # URL to your package repo
    julia = true,          # is it written entirely in Julia?
    license = "MIT",       # your package license
    is_wrapper = false,    # does it wrap around some other package?
)

# Then for each model,
MLJModelInterface.metadata_model(YourModel1,
    input_scitype   = MLJModelInterface.Table(MLJModelInterface.Continuous),  # what input data is supported?
    target_scitype  = AbstractVector{MLJModelInterface.Continuous},           # for a supervised model, what target?
    output_scitype  = MLJModelInterface.Table(MLJModelInterface.Continuous),  # for an unsupervised, what output?
    supports_weights = false,                                                  # does the model support sample weights?
    descr   = "A short description of your model"
	load_path    = "YourPackage.SubModuleContainingModelStructDefinition.YourModel1"
    )
```

*Important.* Do not omit the `load_path` specification. Without a
correct `load_path` MLJ will be unable to import your model.

**Examples**:

- package metadata
  - [GLM](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L179-L186)
  - [MLJLinearModels](https://github.com/JuliaAI/MLJLinearModels.jl/blob/289a373a8357c4afc191711d0218aa1523e97f70/src/mlj/interface.jl#L91-L97)
- model metadata
  - [LinearRegressor](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/GLM.jl#L188-L193)
  - [DecisionTree](https://github.com/JuliaAI/MLJModels.jl/blob/3687491b132be8493b6f7a322aedf66008caaab1/src/DecisionTree.jl#L225-L229)
  - [A series of regressors](https://github.com/JuliaAI/MLJLinearModels.jl/blob/289a373a8357c4afc191711d0218aa1523e97f70/src/mlj/interface.jl#L105-L111)


### Adding a model to the model registry

See [here](https://github.com/JuliaAI/MLJModels.jl/tree/master#instructions-for-updating-the-mlj-model-registry).


# Adding Models for General Use

!!! note

    Models implementing the MLJ model interface according to the instructions given here should import MLJModelInterface version 1.0.0 or higher. This is enforced with a statement such as `MLJModelInterface = "^1" ` under `[compat]` in the Project.toml file of the package containing the implementation.

This guide outlines the specification of the MLJ model interface
and provides detailed guidelines for implementing the interface for
models intended for general use. See also the more condensed
[Quick-Start Guide to Adding Models](@ref).

For sample implementations, see
[MLJModels/src](https://github.com/JuliaAI/MLJModels.jl/tree/master/src/builtins).

Interface code can be hosted by the package providing the core machine
learning algorithm, or by a stand-alone "interface-only" package, using
the template
[MLJExampleInterface.jl](https://github.com/JuliaAI/MLJExampleInterface.jl)
(see [Where to place code implementing new models](@ref) below).

The machine learning tools provided by MLJ can be applied to the
models in any package that imports the package
[MLJModelInterface](https://github.com/JuliaAI/MLJModelInterface.jl) and
implements the API defined there, as outlined below. For a
quick-and-dirty implementation of user-defined models see [Simple User
Defined Models](simple_user_defined_models.md).  To make new models
available to all MLJ users, see [Where to place code implementing new
models](@ref).


#### Important

[MLJModelInterface](https://github.com/JuliaAI/MLJModelInterface.jl)
is a very light-weight interface allowing you to *define* your
interface, but does not provide the functionality required to use or
test your interface; this requires
[MLJBase](https://github.com/JuliaAI/MLJBase.jl).  So,
while you only need to add `MLJModelInterface` to your project's
[deps], for testing purposes you need to add
[MLJBase](https://github.com/JuliaAI/MLJBase.jl) to your
project's [extras] and [targets]. In testing, simply use `MLJBase` in
place of `MLJModelInterface`.

It is assumed the reader has read [Getting Started](index.md).
To implement the API described here, some familiarity with the
following packages is also helpful:

- [ScientificTypes.jl](https://github.com/JuliaAI/ScientificTypes.jl)
  (for specifying model requirements of data)

- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
  (for probabilistic predictions)

- [CategoricalArrays.jl](https://github.com/JuliaData/CategoricalArrays.jl)
  (essential if you are implementing a model handling data of
  `Multiclass` or `OrderedFactor` scitype; familiarity with
  `CategoricalPool` objects required)

- [Tables.jl](https://github.com/JuliaData/Tables.jl) (if your
  algorithm needs input data in a novel format).

In MLJ, the basic interface exposed to the user, built atop the model
interface described here, is the *machine interface*. After a first
reading of this document, the reader may wish to refer to [MLJ
Internals](internals.md) for context.


## Overview

A *model* is an object storing hyperparameters associated with some
machine learning algorithm, and that is all. In MLJ, hyperparameters
include configuration parameters, like the number of threads, and
special instructions, such as "compute feature rankings", which may or
may not affect the final learning outcome.  However, the logging level
(`verbosity` below) is excluded. *Learned parameters* (such as the
coefficients in a linear model) have no place in the model struct.

The name of the Julia type associated with a model indicates the
associated algorithm (e.g., `DecisionTreeClassifier`). The outcome of
training a learning algorithm is called a *fitresult*. For
ordinary multivariate regression, for example, this would be the
coefficients and intercept. For a general supervised model, it is the
(generally minimal) information needed to make new predictions.

The ultimate supertype of all models is `MLJModelInterface.Model`, which
has two abstract subtypes:

```julia
abstract type Supervised <: Model end
abstract type Unsupervised <: Model end
```

`Supervised` models are further divided according to whether they are
able to furnish probabilistic predictions of the target (which they
will then do by default) or directly predict "point" estimates, for each
new input pattern:

```julia
abstract type Probabilistic <: Supervised end
abstract type Deterministic <: Supervised end
```

Further division of model types is realized through [Trait declarations](@ref).

Associated with every concrete subtype of `Model` there must be a
`fit` method, which implements the associated algorithm to produce the
fitresult. Additionally, every `Supervised` model has a `predict`
method, while `Unsupervised` models must have a `transform`
method. More generally, methods such as these, that are dispatched on
a model instance and a fitresult (plus other data), are called
*operations*. `Probabilistic` supervised models optionally implement a
`predict_mode` operation (in the case of classifiers) or a
`predict_mean` and/or `predict_median` operations (in the case of
regressors) although MLJModelInterface also provides fallbacks that will suffice
in most cases. `Unsupervised` models may implement an
`inverse_transform` operation.


## New model type declarations and optional clean! method

Here is an example of a concrete supervised model type declaration,
for a model with a single hyper-parameter:

```julia
import MLJModelInterface
const MMI = MLJModelInterface

mutable struct RidgeRegressor <: MMI.Deterministic
    lambda::Float64
end
```

Models (which are mutable) should not be given internal
constructors. It is recommended that they be given an external lazy
keyword constructor of the same name. This constructor defines default values
for every field, and optionally corrects invalid field values by calling a
`clean!` method (whose fallback returns an empty message string):

```julia
function MMI.clean!(model::RidgeRegressor)
    warning = ""
    if model.lambda < 0
        warning *= "Need lambda ≥ 0. Resetting lambda=0. "
        model.lambda = 0
    end
    return warning
end

# keyword constructor
function RidgeRegressor(; lambda=0.0)
    model = RidgeRegressor(lambda)
    message = MMI.clean!(model)
    isempty(message) || @warn message
    return model
end
```

*Important.* The clean method must have the property that
`clean!(clean!(model)) == clean!(model)` for any instance `model`.

Although not essential, try to avoid `Union` types for model
fields. For example, a field declaration `features::Vector{Symbol}`
with a default of `Symbol[]` (detected with `isempty` method) is
preferred to `features::Union{Vector{Symbol}, Nothing}` with a default
of `nothing`.

### Hyper-parameters for parellizatioin options

The section [Acceleration and Parallelism](@ref) indicates how MLJ
models specify an option to run an algorithm using distributed
processing or multi-threading. A hyper-parameter specifying such an
option should be called `acceleration`. Its value `a` should satisfy 
`a isa AbstractResource` where `AbstractResource` is defined in the
ComputationalResources.jl package. An option to run on a GPU is
ordinarily indicated with the `CUDALibs()` resource.

### Hyper-parameter access and mutation

To support hyper-parameter optimization (see [Tuning Models](@ref))
any hyper-parameter to be individually controlled must be:

- property-accessible; nested property access allowed, as in
  `model.detector.K` 
  
- mutable 

For an un-nested hyper-parameter, the requirement is that
`getproperty(model, :param_name)` and `setproperty!(model,
:param_name, value)` have the expected behavior. (In hyper-parameter
tuning, recursive access is implemented using
[`MLJBase.recursive_getproperty`](@ref)` and
[`MLJBase.recursively_setproperty!`](@ref).)

Combining hyper-parameters in a named tuple does not generally
work, because, although property-accessible (with nesting), an
individual value cannot be mutated. 

For a suggested way to deal with hyper-parameters varying in number,
see the
[implementation](https://github.com/JuliaAI/MLJBase.jl/blob/dev/src/composition/models/stacking.jl)
of `Stack`, where the model struct stores a varying number of base
models internally as a vector, but components are named at
construction and accessed by overloading `getproperty/setproperty!`
appropriately.

### Macro shortcut

An alternative to declaring the model struct, clean! method and
keyword constructor, is to use the `@mlj_model` macro, as in the
following example:

```julia
@mlj_model mutable struct YourModel <: MMI.Deterministic
    a::Float64 = 0.5::(_ > 0)
    b::String  = "svd"::(_ in ("svd","qr"))
end
```

This declaration specifies:

* A keyword constructor (here `YourModel(; a=..., b=...)`),
* Default values for the hyperparameters,
* Constraints on the hyperparameters where `_` refers to a value
  passed.

For example, `a::Float64 = 0.5::(_ > 0)` indicates that
the field `a` is a `Float64`, takes `0.5` as default value, and
expects its value to be positive.

You cannot use the `@mlj_model` macro if your model struct has type
parameters.

#### Known issue with @mlj_macro

Defaults with negative values can trip up the `@mlj_macro` (see [this
issue](https://github.com/JuliaAI/MLJBase.jl/issues/68)). So,
for example, this does not work:

```julia
@mlj_model mutable struct Bar
    a::Int = -1::(_ > -2)
end
```

But this does:

```julia
@mlj_model mutable struct Bar
    a::Int = (-)(1)::(_ > -2)
end
```


## Supervised models

### Mathematical assumptions

At present, MLJ's performance estimate functionality (resampling using
`evaluate`/`evaluate!`) tacitly assumes that feature-label pairs of
observations `(X1, y1), (X2, y2), (X2, y2), ...` are being modelled as
identically independent random variables (i.i.d.), and constructs some
kind of representation of an estimate of the conditional probablility
`p(y | X)` (`y` and `X` single observations). It may be that a model
implementing the MLJ interface has the potential to make predictions
under weaker assumptions (e.g., time series forecasting
models). However the output of the compulsory `predict` method
described below should be the output of the model under the i.i.d
assumption.

In the future newer methods may be introduced to handle weaker
assumptions (see, e.g., [The predict_joint method](@ref) below).


### Summary of methods

The compulsory and optional methods to be implemented for each
concrete type `SomeSupervisedModel <: MMI.Supervised` are
summarized below. 

An `=` indicates the return value for a fallback version of the
method.

Compulsory:

```julia
MMI.fit(model::SomeSupervisedModel, verbosity, X, y) -> fitresult, cache, report
MMI.predict(model::SomeSupervisedModel, fitresult, Xnew) -> yhat
```

Optional, to check and correct invalid hyperparameter values:

```julia
MMI.clean!(model::SomeSupervisedModel) = ""
```

Optional, to return user-friendly form of fitted parameters:

```julia
MMI.fitted_params(model::SomeSupervisedModel, fitresult) = fitresult
```

Optional, to avoid redundant calculations when re-fitting machines
associated with a model:

```julia
MMI.update(model::SomeSupervisedModel, verbosity, old_fitresult, old_cache, X, y) =
   MMI.fit(model, verbosity, X, y)
```

Optional, to specify default hyperparameter ranges (for use in tuning):

```julia
MMI.hyperparameter_ranges(T::Type) = Tuple(fill(nothing, length(fieldnames(T))))
```

Optional, if `SomeSupervisedModel <: Probabilistic`:

```julia
MMI.predict_mode(model::SomeSupervisedModel, fitresult, Xnew) =
    mode.(predict(model, fitresult, Xnew))
MMI.predict_mean(model::SomeSupervisedModel, fitresult, Xnew) =
    mean.(predict(model, fitresult, Xnew))
MMI.predict_median(model::SomeSupervisedModel, fitresult, Xnew) =
    median.(predict(model, fitresult, Xnew))
```

Required, if the model is to be registered (findable by general users):

```julia
MMI.load_path(::Type{<:SomeSupervisedModel})    = ""
MMI.package_name(::Type{<:SomeSupervisedModel}) = "Unknown"
MMI.package_uuid(::Type{<:SomeSupervisedModel}) = "Unknown"
```

```julia
MMI.input_scitype(::Type{<:SomeSupervisedModel}) = Unknown
```

Strongly recommended, to constrain the form of target data passed to fit:

```julia
MMI.target_scitype(::Type{<:SomeSupervisedModel}) = Unknown
```

Optional but recommended:

```julia
MMI.package_url(::Type{<:SomeSupervisedModel})  = "unknown"
MMI.is_pure_julia(::Type{<:SomeSupervisedModel}) = false
MMI.package_license(::Type{<:SomeSupervisedModel}) = "unknown"
```

If `SomeSupervisedModel` supports sample weights or class weights,
then instead of the `fit` above, one implements

```julia
MMI.fit(model::SomeSupervisedModel, verbosity, X, y, w=nothing) -> fitresult, cache, report
```

and, if appropriate

```julia
MMI.update(model::SomeSupervisedModel, verbosity, old_fitresult, old_cache, X, y, w=nothing) =
   MMI.fit(model, verbosity, X, y, w)
```

Additionally, if `SomeSupervisedModel` supports sample weights, one must declare

```julia
MMI.supports_weights(model::Type{<:SomeSupervisedModel}) = true
```

Optionally, an implemenation may add a data front-end, for
transforming user data (such as a table) into some model-specific
format (such as a matrix), and for adding methods to specify how said
format is resampled. (This alters the meaning of `X`, `y` and `w` in
the signatures of `fit`, `update`, `predict`, etc; see [Implementing a
data front-end](@ref) for details). This can provide the MLJ user
certain performance advantages when fitting a machine.

```julia
MLJModelInterface.reformat(model::SomeSupervisedModel, args...) = args
MLJModelInterface.selectrows(model::SomeSupervisedModel, I, data...) = data
```

Optionally, to customized support for serialization of machines (see
[Serialization](@ref)), overload

```julia
MMI.save(filename, model::SomeModel, fitresult; kwargs...) = fitresult
```

and possibly

```julia
MMI.restore(filename, model::SomeModel, serializable_fitresult) -> serializable_fitresult
```

These last two are unlikely to be needed if wrapping pure Julia code.


### The form of data for fitting and predicting

The model implementer does not have absolute control over the types of
data `X`, `y` and `Xnew` appearing in the `fit` and `predict` methods
they must implement. Rather, they can specify the *scientific type* of
this data by making appropriate declarations of the traits
`input_scitype` and `target_scitype` discussed later under [Trait
declarations](@ref).

*Important Note.* Unless it genuinely makes little sense to do so, the
MLJ recommendation is to specify a `Table` scientific type for `X`
(and hence `Xnew`) and an `AbstractVector` scientific type (e.g.,
`AbstractVector{Continuous}`) for targets `y`. Algorithms requiring
matrix input can coerce their inputs appropriately; see below.


#### Additional type coercions

If the core algorithm being wrapped requires data in a different or
more specific form, then `fit` will need to coerce the table into the
form desired (and the same coercions applied to `X` will have to be
repeated for `Xnew` in `predict`). To assist with common cases, MLJ
provides the convenience method
[`MMI.matrix`](@ref). `MMI.matrix(Xtable)` has type `Matrix{T}` where
`T` is the tightest common type of elements of `Xtable`, and `Xtable`
is any table. (If `Xtable` is itself just a wrapped matrix,
`Xtable=Tables.table(A)`, then `A=MMI.table(Xtable)` will be returned
without any copying.)

Alternatively, a more performant option is to implement a data
front-end for your model; see [Implementing a data front-end](@ref).

Other auxiliary methods provided by MLJModelInterface for handling tabular data
are: `selectrows`, `selectcols`, `select` and `schema` (for extracting
the size, names and eltypes of a table's columns). See [Convenience
methods](@ref) below for details.


#### Important convention

It is to be understood that the columns of the table `X` correspond to
features and the rows to observations. So, for example, the predict
method for a linear regression model might look like `predict(model,
w, Xnew) = MMI.matrix(Xnew)*w`, where `w` is the vector of learned
coefficients.


### The fit method

A compulsory `fit` method returns three objects:

```julia
MMI.fit(model::SomeSupervisedModel, verbosity, X, y) -> fitresult, cache, report
```

1. `fitresult` is the fitresult in the sense above (which becomes an
    argument for `predict` discussed below).

2.  `report` is a (possibly empty) `NamedTuple`, for example,
    `report=(deviance=..., dof_residual=..., stderror=..., vcov=...)`.
    Any training-related statistics, such as internal estimates of the
    generalization error, and feature rankings, should be returned in
    the `report` tuple. How, or if, these are generated should be
    controlled by hyperparameters (the fields of `model`). Fitted
    parameters, such as the coefficients of a linear model, do not go
    in the report as they will be extractable from `fitresult` (and
    accessible to MLJ through the `fitted_params` method described below).

3.	The value of `cache` can be `nothing`, unless one is also defining
    an `update` method (see below). The Julia type of `cache` is not
    presently restricted.
	
!!! note

    The  `fit` (and `update`) methods should not mutate the `model`. If necessary, `fit` can create a `deepcopy` of `model` first. 


It is not necessary for `fit` to provide type or dimension checks on
`X` or `y` or to call `clean!` on the model; MLJ will carry out such
checks. 

The types of `X` and `y` are constrained by the `input_scitype` and
`target_scitype` trait declarations; see [Trait declarations](@ref)
below. (That is, unless a data front-end is implemented, in which case
these traits refer instead to the arguments of the overloaded
`reformat` method, and the types of `X` and `y` are determined by the
output of `reformat`.)

The method `fit` should never alter hyperparameter values, the sole
exception being fields of type `<:AbstractRNG`. If the package is able
to suggest better hyperparameters, as a byproduct of training, return
these in the report field.

The `verbosity` level (0 for silent) is for passing to learning
algorithm itself. A `fit` method wrapping such an algorithm should
generally avoid doing any of its own logging.

*Sample weight support.* If
`supports_weights(::Type{<:SomeSupervisedModel})` has been declared
`true`, then one instead implements the following variation on the
above `fit`:

```julia
MMI.fit(model::SomeSupervisedModel, verbosity, X, y, w=nothing) -> fitresult, cache, report
```


### The fitted_params method

A `fitted_params` method may be optionally overloaded. It's purpose is
to provide MLJ access to a user-friendly representation of the
learned parameters of the model (as opposed to the
hyperparameters). They must be extractable from `fitresult`.

```julia
MMI.fitted_params(model::SomeSupervisedModel, fitresult) -> friendly_fitresult::NamedTuple
```

For a linear model, for example, one might declare something like
`friendly_fitresult=(coefs=[...], bias=...)`.

The fallback is to return `(fitresult=fitresult,)`.


### The predict method

A compulsory `predict` method has the form

```julia
MMI.predict(model::SomeSupervisedModel, fitresult, Xnew) -> yhat
```

Here `Xnew` will have the same form as the `X` passed to
`fit`. 

Note that while `Xnew` generally consists of multiple observations
(e.g., has multiple rows in the case of a table) it is assumed, in view of
the i.i.d assumption recalled above, that calling `predict(..., Xnew)`
is equivalent to broadcasting some method `predict_one(..., x)` over
the individual observations `x` in `Xnew` (a method implementing the
probablility distribution `p(X |y)` above).


#### Prediction types for deterministic responses.

In the case of `Deterministic` models, `yhat` should have the same
scitype as the `y` passed to `fit` (see above). If `y` is a
`CategoricalVector` (classification) then elements of the predition
`yhat` **must have a pool == to the pool of the target `y` presented
in training**, even if not all levels appear in the training data or
prediction itself.

Unfortunately, code not written with the preservation of categorical
levels in mind poses special problems. To help with this,
MLJModelInterface provides some utilities:
[`MLJModelInterface.int`](@ref) (for converting a `CategoricalValue`
into an integer, the ordering of these integers being consistent with
that of the pool) and `MLJModelInterface.decoder` (for constructing a
callable object that decodes the integers back into `CategoricalValue`
objects). Refer to [Convenience methods](@ref) below for important
details.

Note that a decoder created during `fit` may need to be bundled with
`fitresult` to make it available to `predict` during re-encoding. So,
for example, if the core algorithm being wrapped by `fit` expects a
nominal target `yint` of type `Vector{<:Integer}` then a `fit` method
may look something like this:

```julia
function MMI.fit(model::SomeSupervisedModel, verbosity, X, y)
    yint = MMI.int(y)
    a_target_element = y[1]                # a CategoricalValue/String
    decode = MMI.decoder(a_target_element) # can be called on integers

    core_fitresult = SomePackage.fit(X, yint, verbosity=verbosity)

    fitresult = (decode, core_fitresult)
    cache = nothing
    report = nothing
    return fitresult, cache, report
end
```

while a corresponding deterministic `predict` operation might look like this:

```julia
function MMI.predict(model::SomeSupervisedModel, fitresult, Xnew)
    decode, core_fitresult = fitresult
    yhat = SomePackage.predict(core_fitresult, Xnew)
    return decode.(yhat)
end
```

For a concrete example, refer to the
[code](https://github.com/JuliaAI/MLJModels.jl/blob/master/src/ScikitLearn.jl)
for `SVMClassifier`.

Of course, if you are coding a learning algorithm from scratch, rather
than wrapping an existing one, these extra measures may be unnecessary.


#### Prediction types for probabilistic responses

In the case of `Probabilistic` models with univariate targets, `yhat`
must be an `AbstractVector` or table whose elements are distributions.
In the common case of a vector (single target), this means one
distribution per row of `Xnew`.

A *distribution* is some object that, at the least, implements
`Base.rng` (i.e., is something that can be sampled).  Currently, all
performance measures (metrics) defined in MLJBase.jl additionally
assume that a distribution is either:

- An instance of some subtype of `Distributions.Distribution`, an
  abstract type defined in the
  [`Distributions.jl`](https://juliastats.org/Distributions.jl/stable/)
  package; or
  
- An instance of `CategoricalDistributions.UnivariateFinite`, from the
  [CategoricalDistributions.jl](https://github.com/JuliaAI/CategoricalDistributions.jl)
  package, *which should be used for all probabilistic classifiers*,
  i.e., for predictors whose target has scientific type
  `<:AbstractVector{<:Finite}`.

All such distributions implement the probability mass or density
function `Distributions.pdf`. If your model's predictions cannot be
predict objects of this form, then you will need to implement
appropriate performance measures to buy into MLJ's performance
evaluation apparatus.

An implementation can avoid CategoricalDistributions.jl as a
dependency by using the "dummy" constructor
`MLJModelInterface.UnivariateFinite`, which is bound to the true one
when MLJBase.jl is loaded.

For efficiency, one should not construct `UnivariateFinite` instances
one at a time. Rather, once a probability vector, matrix, or
dictionary is known, construct an instance of `UnivariateFiniteVector
<: AbstractArray{<:UnivariateFinite},1}` to return. Both
`UnivariateFinite` and `UnivariateFiniteVector` objects are
constructed using the single `UnivariateFinite` function.

For example, suppose the target `y` arrives as a subsample of some
`ybig` and is missing some classes:

```julia
ybig = categorical([:a, :b, :a, :a, :b, :a, :rare, :a, :b])
y = ybig[1:6]
```

Your fit method has bundled the first element of `y` with the
`fitresult` to make it available to `predict` for purposes of tracking
the complete pool of classes. Let's call this `an_element =
y[1]`. Then, supposing the corresponding probabilities of the observed
classes `[:a, :b]` are in an `n x 2` matrix `probs` (where `n` the number of
rows of `Xnew`) then you return

```julia
yhat = MLJModelInterface.UnivariateFinite([:a, :b], probs, pool=an_element)
```

This object automatically assigns zero-probability to the unseen class
`:rare` (i.e., `pdf.(yhat, :rare)` works and returns a zero
vector). If you would like to assign `:rare` non-zero probabilities,
simply add it to the first vector (the *support*) and supply a larger
`probs` matrix.

In a binary classification problem it suffices to specify a single
vector of probabilities, provided you specify `augment=true`, as in
the following example, *and note carefully that these probablities are
associated with the* **last** *(second) class you specify in the
constructor:*

```julia
y = categorical([:TRUE, :FALSE, :FALSE, :TRUE, :TRUE])
an_element = y[1]
probs = rand(10)
yhat = MLJModelInterface.UnivariateFinite([:FALSE, :TRUE], probs, augment=true, pool=an_element)
```

The constructor has a lot of options, including passing a dictionary
instead of vectors. See
`CategoricalDistributions.UnivariateFinite`](@ref) for details.

See
[LinearBinaryClassifier](https://github.com/JuliaAI/MLJModels.jl/blob/master/src/GLM.jl)
for an example of a Probabilistic classifier implementation.

*Important note on binary classifiers.* There is no "Binary" scitype
distinct from `Multiclass{2}` or `OrderedFactor{2}`; `Binary` is just
an alias for `Union{Multiclass{2},OrderedFactor{2}}`. The
`target_scitype` of a binary classifier will generally be
`AbstractVector{<:Binary}` and according to the *mlj* scitype
convention, elements of `y` have type `CategoricalValue`, and *not*
`Bool`. See
[BinaryClassifier](https://github.com/JuliaAI/MLJModels.jl/blob/master/src/GLM.jl)
for an example.


### The predict_joint method

!!! warning "Experimental"

    The following API is experimental. It is subject to breaking changes during minor or major releases without warning.
	
```julia
MMI.predict_joint(model::SomeSupervisedModel, fitresult, Xnew) -> yhat
```

Any `Probabilistic` model type `SomeModel`may optionally implement a
`predict_joint` method, which has the same signature as `predict`, but
whose predictions are a single distribution (rather than a vector of
per-observation distributions). 

Specifically, the output `yhat` of `predict_joint` should be an
instance of `Distributions.Sampleable{<:Multivariate,V}`, where
`scitype(V) = target_scitype(SomeModel)` and samples have length `n`,
where `n` is the number of observations in `Xnew`.

If a new model type subtypes `JointProbablistic <: Probabilistic` then
implementation of `predict_joint` is compulsory.


### Trait declarations

Two trait functions allow the implementer to restrict the types of
data `X`, `y` and `Xnew` discussed above. The MLJ task interface uses
these traits for data type checks but also for model search. If they
are omitted (and your model is registered) then a general user may
attempt to use your model with inappropriately typed data.

The trait functions `input_scitype` and `target_scitype` take
scientific data types as values. We assume here familiarity with
[ScientificTypes.jl](https://github.com/JuliaAI/ScientificTypes.jl)
(see [Getting Started](index.md) for the basics).

For example, to ensure that the `X` presented to the
`DecisionTreeClassifier` `fit` method is a table whose columns all
have `Continuous` element type (and hence `AbstractFloat` machine
type), one declares

```julia
MMI.input_scitype(::Type{<:DecisionTreeClassifier}) = MMI.Table(MMI.Continuous)
```

or, equivalently,

```julia
MMI.input_scitype(::Type{<:DecisionTreeClassifier}) = Table(Continuous)
```

If, instead, columns were allowed to have either: (i) a mixture of `Continuous` and `Missing`
values, or (ii) `Count` (i.e., integer) values, then the
declaration would be

```julia
MMI.input_scitype(::Type{<:DecisionTreeClassifier}) = Table(Union{Continuous,Missing},Count)
```

Similarly, to ensure the target is an AbstractVector whose elements
have `Finite` scitype (and hence `CategoricalValue` machine type) we declare

```julia
MMI.target_scitype(::Type{<:DecisionTreeClassifier}) = AbstractVector{<:Finite}
```

#### Multivariate targets

The above remarks continue to hold unchanged for the case multivariate
targets.  For example, if we declare

```julia
target_scitype(SomeSupervisedModel) = Table(Continuous)
```

then this constrains the target to be any table whose columns have `Continous` element scitype (i.e., `AbstractFloat`), while

```julia
target_scitype(SomeSupervisedModel) = Table(Continuous, Finite{2})
```

restricts to tables with continuous or binary (ordered or unordered)
columns.

For predicting variable length sequences of, say, binary values
(`CategoricalValue`s) with some common size-two pool) we declare

```julia
target_scitype(SomeSupervisedModel) = AbstractVector{<:NTuple{<:Finite{2}}}
```

The trait functions controlling the form of data are summarized as follows:

method                   | return type       | declarable return values     | fallback value
-------------------------|-------------------|------------------------------|---------------
`input_scitype`          | `Type`            | some scientfic type          | `Unknown`
`target_scitype`         | `Type`            | some scientific type         | `Unknown`


Additional trait functions tell MLJ's `@load` macro how to find your
model if it is registered, and provide other self-explanatory metadata
about the model:

method                   | return type       | declarable return values           | fallback value
-------------------------|-------------------|------------------------------------|---------------
`load_path`              | `String`          | unrestricted                       | "unknown"
`package_name`           | `String`          | unrestricted                       | "unknown"
`package_uuid`           | `String`          | unrestricted                       | "unknown"
`package_url`            | `String`          | unrestricted                       | "unknown"
`package_license`        | `String`          | unrestricted                       | "unknown"
`is_pure_julia`          | `Bool`            | `true` or `false`                  | `false`
`supports_weights`       | `Bool`            | `true` or `false`                  | `false`

Here is the complete list of trait function declarations for
`DecisionTreeClassifier`, whose core algorithms are provided by
DecisionTree.jl, but whose interface actually lives at
[MLJDecisionTreeInterface.jl](https://github.com/JuliaAI/MLJDecisionTreeInterface.jl).

```julia
MMI.input_scitype(::Type{<:DecisionTreeClassifier}) = MMI.Table(MMI.Continuous)
MMI.target_scitype(::Type{<:DecisionTreeClassifier}) = AbstractVector{<:MMI.Finite}
MMI.load_path(::Type{<:DecisionTreeClassifier}) = "MLJDecisionTreeInterface.DecisionTreeClassifier"
MMI.package_name(::Type{<:DecisionTreeClassifier}) = "DecisionTree"
MMI.package_uuid(::Type{<:DecisionTreeClassifier}) = "7806a523-6efd-50cb-b5f6-3fa6f1930dbb"
MMI.package_url(::Type{<:DecisionTreeClassifier}) = "https://github.com/bensadeghi/DecisionTree.jl"
MMI.is_pure_julia(::Type{<:DecisionTreeClassifier}) = true
```

Alternatively these traits can also be declared using `MMI.metadata_pkg` and `MMI.metadata_model` helper functions as:

```julia
MMI.metadata_pkg(DecisionTreeClassifier,name="DecisionTree",
                     packge_uuid="7806a523-6efd-50cb-b5f6-3fa6f1930dbb",
                     package_url="https://github.com/bensadeghi/DecisionTree.jl",
                     is_pure_julia=true)

MMI.metadata_model(DecisionTreeClassifier,
                        input_scitype=MMI.Table(MMI.Continuous),
                        target_scitype=AbstractVector{<:MMI.Finite},
                        load_path="MLJDecisionTreeInterface.DecisionTreeClassifier")
```

*Important.* Do not omit the `load_path` specification. If unsure what
it should be, post an issue at
[MLJ](https://github.com/alan-turing-institute/MLJ.jl/issues).

```@docs
MMI.metadata_pkg
```

```@docs
MMI.metadata_model
```


### Iterative models and the update! method

An `update` method may be optionally overloaded to enable a call by
MLJ to retrain a model (on the same training data) to avoid repeating
computations unnecessarily.

```julia
MMI.update(model::SomeSupervisedModel, verbosity, old_fitresult, old_cache, X, y) -> fit
result, cache, report
MMI.update(model::SomeSupervisedModel, verbosity, old_fitresult, old_cache, X, y, w=nothing) -> fit
result, cache, report
```

Here the second variation applies if `SomeSupervisedModel` supports
sample weights.

If an MLJ `Machine` is being `fit!` and it is not the first time, then
`update` is called instead of `fit`, unless the machine `fit!` has
been called with a new `rows` keyword argument. However, `MLJModelInterface`
defines a fallback for `update` which just calls `fit`. For context,
see [MLJ Internals](internals.md).

Learning networks wrapped as models constitute one use-case (see
[Composing Models](index.md)): one would like each component model to
be retrained only when hyperparameter changes "upstream" make this
necessary. In this case MLJ provides a fallback (specifically, the
fallback is for any subtype of `SupervisedNetwork =
Union{DeterministicNetwork,ProbabilisticNetwork}`). A second more
generally relevant use-case is iterative models, where calls to
increase the number of iterations only restarts the iterative
procedure if other hyperparameters have also changed. (A useful method
for inspecting model changes in such cases is
`MLJModelInterface.is_same_except`. ) For an example, see
[MLJEnsembles.jl](https://github.com/JuliaAI/MLJEnsembles.jl).

A third use-case is to avoid repeating time-consuming preprocessing of
`X` and `y` required by some models.

In the event that the argument `fitresult` (returned by a preceding
call to `fit`) is not sufficient for performing an update, the author
can arrange for `fit` to output in its `cache` return value any
additional information required (for example, pre-processed versions
of `X` and `y`), as this is also passed as an argument to the `update`
method.

### Implementing a data front-end

!!! note

    It is suggested that packages implementing MLJ's model API, that later implement a data front-end, should tag their changes in a breaking release. (The changes will not break use of models for the ordinary MLJ user, who interacts with models exlusively through the machine interface. However, it will break usage for some external packages that have chosen to depend directly on the model API.)

```julia
MLJModelInterface.reformat(model, args...) -> data
MLJModelInterface.selectrows(::Model, I, data...) -> sampled_data
```

Models optionally overload `reformat` to define transformations of
user-supplied data into some model-specific representation (e.g., from
a table to a matrix). Computational overheads associated with multiple
`fit!`/`predict`/`transform` calls (on MLJ machines) are then avoided,
when memory resources allow. The fallback returns `args` (no
transformation). 

The `selectrows(::Model, I, data...)` method is overloaded to specify
how the model-specific data is to be subsampled, for some observation
indices `I` (a colon, `:`, or instance of
`AbstractVector{<:Integer}`). In this way, implementing a data
front-end also allow more efficient resampling of data (in user calls
to `evaluate!`).

After detailing formal requirments for implementing a data front-end,
we give a [Sample implementation](@ref). A simple implementation
[implementation](https://github.com/Evovest/EvoTrees.jl/blob/94b58faf3042009bd609c9a5155a2e95486c2f0e/src/MLJ.jl#L23)
also appears in the EvoTrees.jl package.

Here "user-supplied data" is what the MLJ user supplies when
constructing a machine, as in `machine(models, args...)`, which
coincides with the arguments expected by `fit(model, verbosity,
args...)` when `reformat` is not overloaded.

Implementing a `reformat` data front-end is permitted for any `Model`
subtype, except for subtypes of `Static`. Here is a complete list of
responsibilities for such an implementation, for some
`model::SomeModelType` (a sample implementation follows after):

- A `reformat(model::SomeModelType, args...) -> data` method must be
  implemented for each form of `args...` appearing in a valid machine
  construction `machine(model, args...)` (there will be one for each
  possible signature of `fit(::SomeModelType, ...)`).

- Additionally, if not included above, there must be a single argument
  form of reformat, `reformat(model::SommeModelType, arg) -> (data,)`,
  serving as a data front-end for operations like `predict`. It must
  always hold that `reformat(model, args...)[1] = reformat(model,
  args[1])`.

*Important.* `reformat(model::SomeModelType, args...)` must always
  return a tuple of the same length as `args`, even if this is one.

- `fit(model::SomeModelType, verbosity, data...)` should be
  implemented as if `data` is the output of `reformat(model,
  args...)`, where `args` is the data an MLJ user has bound to `model`
  in some machine. The same applies to any overloading of `update`.

- Each implemented operation, such as `predict` and `transform` - but
  excluding `inverse_transform` - must be defined as if its data
  arguments are `reformat`ed versions of user-supplied data. For
  example, in the supervised case, `data_new` in
  `predict(model::SomeModelType, fitresult, data_new)` is
  `reformat(model, Xnew)`, where `Xnew` is the data provided by the MLJ
  user in a call `predict(mach, Xnew)` (`mach.model == model`).

- To specify how the model-specific representation of data is to be
  resampled, implement `selectrows(model::SomeModelType, I, data...)
  -> resampled_data` for each overloading of `reformat(model::SomeModel,
  args...) -> data` above. Here `I` is an arbitrary abstract integer
  vector or `:` (type `Colon`).

*Important.* `selectrows(model::SomeModelType, I, args...)` must always
return a tuple of the same length as `args`, even if this is one.

The fallback for `selectrows` is described at [`selectrows`](@ref).


#### Sample implementation

Suppose a supervised model type `SomeSupervised` supports sample
weights, leading to two different `fit` signatures, and that it has a
single operation `predict`:

    fit(model::SomeSupervised, verbosity, X, y)
    fit(model::SomeSupervised, verbosity, X, y, w)

    predict(model::SomeSupervised, fitresult, Xnew)

Without a data front-end implemented, suppose `X` is expected to be a
table and `y` a vector, but suppose the core algorithm always converts
`X` to a matrix with features as rows (features corresponding to
columns in the table).  Then a new data-front end might look like
this:

    constant MMI = MLJModelInterface

    # for fit:
    MMI.reformat(::SomeSupervised, X, y) = (MMI.matrix(X, transpose=true), y)
    MMI.reformat(::SomeSupervised, X, y, w) = (MMI.matrix(X, transpose=true), y, w)
    MMI.selectrows(::SomeSupervised, I, Xmatrix, y) =
        (view(Xmatrix, :, I), view(y, I))
    MMI.selectrows(::SomeSupervised, I, Xmatrix, y, w) =
        (view(Xmatrix, :, I), view(y, I), view(w, I))

    # for predict:
    MMI.reformat(::SomeSupervised, X) = (MMI.matrix(X, transpose=true),)
    MMI.selectrows(::SomeSupervised, I, Xmatrix) = view(Xmatrix, I)

With these additions, `fit` and `predict` are refactored, so that `X`
and `Xnew` represent matrices with features as rows.


### Supervised models with a `transform` method

A supervised model may optionally implement a `transform` method,
whose signature is the same as `predict`. In that case the
implementation should define a value for the `output_scitype` trait. A
declaration

```julia
output_scitype(::Type{<:SomeSupervisedModel}) = T
```

is an assurance that `scitype(transform(model, fitresult, Xnew)) <: T`
always holds, for any `model` of type `SomeSupervisedModel`.

A use-case for a `transform` method for a supervised model is a neural
network that learns *feature embeddings* for categorical input
features as part of overall training. Such a model becomes a
transformer that other supervised models can use to transform the
categorical features (instead of applying the higher-dimensional one-hot
encoding representations).

## Models that learn a probability distribution


!!! warning "Experimental"

    The following API is experimental. It is subject to breaking changes during minor or major releases without warning. Models implementing this interface will not work with MLJBase versions earlier than 0.17.5.

Models that fit a probability distribution to some `data` should be
regarded as `Probablisitic <: Supervised` models with target `y = data`
and `X = nothing`. 

The `predict` method should return a single distribution. 

A working implementation of a model that fits a `UnivariateFinite`
distribution to some categorical data using [Laplace
smoothing](https://en.wikipedia.org/wiki/Additive_smoothing)
controlled by a hyper-parameter `alpha` is given
[here](https://github.com/JuliaAI/MLJBase.jl/blob/d377bee1198ec179a4ade191c11fef583854af4a/test/interface/model_api.jl#L36).


### Serialization 

!!! warning "Experimental"

    The following API is experimental. It is subject to breaking changes during minor or major releases without warning.

The MLJ user can serialize and deserialize a *machine*, which means
serializing/deserializing:

- the associated `Model` object (storing hyperparameters)
- the `fitresult` (learned parameters)
- the `report` generating during training

These are bundled into a single file or `IO` stream specified by the
user using the package `JLSO`. There are two scenarios in which a new
MLJ model API implementation will want to overload two additional
methods `save` and `restore` to support serialization:

1. The algorithm-providing package already has it's own serialization format for learned parameters and/or hyper-parameters, which users may want to access. In that case *the implementation overloads* `save`.
  
2. The `fitresult` is not a sufficiently persistent object; for example, it is a pointer passed from wrapped C code. In that case *the implementation overloads* `save` *and* `restore`.
  
In case 2, 1 presumably applies also, for otherwise MLJ serialization
is probably not going to be possible without changes to the
algorithm-providing package. An example is given below.

Note that in case 1, MLJ will continue to create it's own
self-contained serialization of the machine. Below `filename` refers
to the corresponding serialization file name, as specified by the
user, but with any final extension (e.g., ".jlso", ".gz") removed. If
the user has alternatively specified an `IO` object for serialization,
then `filename` is a randomly generated numeric string.


#### The save method

```julia
MMI.save(filename, model::SomeModel, fitresult; kwargs...) -> serializable_fitresult
```

Implement this method to serialize using a format specific to models
of type `SomeModel`. The `fitresult` is the first return value of
`MMI.fit` for such model types; `kwargs` is a list of keyword
arguments specified by the user and understood to relate to a some
model-specific serialization (cannot be `format=...` or
`compression=...`). The value of `serializable_fitresult` should be a
persistent representation of `fitresult`, from which a correct and
valid `fitresult` can be reconstructed using `restore` (see
below). 

The fallback of `save` performs no action and returns `fitresult`.


#### The restore method

```julia
MMI.restore(filename, model::SomeModel, serializable_fitresult) -> fitresult
```

Implement this method to reconstruct a `fitresult` (as returned by
`MMI.fit`) from a persistent representation constructed using
`MMI.save` as described above. 

The fallback of `restore` returns `serializable_fitresult`.

#### Example

Below is an example drawn from MLJ's XGBoost wrapper. In this example
the `fitresult` returned by `MMI.fit` is a tuple `(booster,
a_target_element)` where `booster` is the `XGBoost.jl` object storing
the learned parameters (essentially a pointer to some object created
by C code) and `a_target_element` is an ordinary `CategoricalValue`
used to track the target classes (a persistent object, requiring no
special treatment).

```julia
function MLJModelInterface.save(filename,
                                ::XGBoostClassifier,
                                fitresult;
                                kwargs...)
    booster, a_target_element = fitresult

    xgb_filename = string(filename, ".xgboost.model")
    XGBoost.save(booster, xgb_filename)
    persistent_booster = read(xgb_filename)
    @info "Additional XGBoost serialization file \"$xgb_filename\" generated. "
    return (persistent_booster, a_target_element)
end

function MLJModelInterface.restore(filename,
                                   ::XGBoostClassifier,
                                   serializable_fitresult)
    persistent_booster, a_target_element = serializable_fitresult

    xgb_filename = string(filename, ".tmp")
    open(xgb_filename, "w") do file
        write(file, persistent_booster)
    end
    booster = XGBoost.Booster(model_file=xgb_filename)
    rm(xgb_filename)
    fitresult = (booster, a_target_element)
    return fitresult
end
```

## Unsupervised models

Unsupervised models implement the MLJ model interface in a very
similar fashion. The main differences are:

- The `fit` method has only one training argument `X`, as in
  `MLJModelInterface.fit(model, verbosity, X)`. However, it has
  the same return value `(fitresult, cache, report)`. An `update`
  method (e.g., for iterative models) can be optionally implemented in
  the same way.

- A `transform` method is compulsory and has the same signature as
  `predict`, as in `MLJModelInterface.transform(model, fitresult, Xnew)`.

- Instead of defining the `target_scitype` trait, one declares an
  `output_scitype` trait (see above for the meaning).

- An `inverse_transform` can be optionally implemented. The signature
  is the same as `transform`, as in
  `MLJModelInterface.inverse_transform(model, fitresult, Xout)`, which:

   - must make sense for any `Xout` for which `scitype(Xout) <:
     output_scitype(SomeSupervisedModel)` (see below); and

   - must return an object `Xin` satisfying `scitype(Xin) <:
     input_scitype(SomeSupervisedModel)`.

- A `predict` method may be optionally implemented, and has the same
  signature as for supervised models, as in
  `MLJModelInterface.predict(model, fitresult, Xnew)`. A use-case is
  clustering algorithms that `predict` labels and `transform` new
  input features into a space of lower-dimension. See [Transformers
  that also predict](@ref) for an example.


## Convenience methods

```@docs
MLJBase.table
MLJBase.matrix
```

```@docs
MLJModelInterface.int
```

```@docs
CategoricalDistributions.UnivariateFinite
```

```@docs
CategoricalDistributions.classes
```

```@docs
MLJModelInterface.decoder
```

```@docs
MLJModelInterface.select
```

```@docs
MLJModelInterface.selectrows
```

```@docs
MLJModelInterface.selectcols
```

```@docs
MLJBase.recursive_getproperty
MLJBase.recursive_setproperty!
```


### Where to place code implementing new models

Note that different packages can implement models having the same name
without causing conflicts, although an MLJ user cannot simultaneously
*load* two such models.

There are two options for making a new model implementation available
to all MLJ users:

1. **Native implementations** (preferred option). The implementation
   code lives in the same package that contains the learning
   algorithms implementing the interface. An example is
   [`EvoTrees.jl`](https://github.com/Evovest/EvoTrees.jl/blob/master/src/MLJ.jl). In
   this case, it is sufficient to open an issue at
   [MLJ](https://github.com/alan-turing-institute/MLJ.jl) requesting
   the package to be registered with MLJ. Registering a package allows
   the MLJ user to access its models' metadata and to selectively load
   them.

2. **Separate interface package**. Implementation code lives in a
   separate *interface package*, which has the algorithm providing
   package as a dependency. See the template repository
   [MLJExampleInterface.jl](https://github.com/JuliaAI/MLJExampleInterface.jl).

Additionally, one needs to ensure that the implementation code defines
the `package_name` and `load_path` model traits appropriately, so that
`MLJ`'s `@load` macro can find the necessary code (see
[MLJModels/src](https://github.com/JuliaAI/MLJModels.jl/tree/master/src)
for examples).

### How to add models to the MLJ model registry?

The MLJ model registry is located in the [MLJModels.jl
repository](https://github.com/JuliaAI/MLJModels.jl). To
add a model, you need to follow these steps

- Ensure your model conforms to the interface defined above

- Raise an issue at
  [MLJModels.jl](https://github.com/JuliaAI/MLJModels.jl/issues)
  and point out where the MLJ-interface implementation is, e.g. by
  providing a link to the code.

- An administrator will then review your implementation and work with
  you to add the model to the registry

# Homogeneous Ensembles

Although an ensemble of models sharing a common set of hyperparameters
can defined using the learning network API, MLJ's `EnsembleModel`
model wrapper is preferred, for convenience and best performance.

```@docs
MLJEnsembles.EnsembleModel
```


# Modifying Behavior

To modify behaviour of MLJ you will need to clone the relevant
component package (e.g., MLJBase.jl) - or a fork thereof - and modify
your local julia environment to use your local clone in place of the
official release. For example, you might proceed something like this:

```julia
using Pkg
Pkg.activate("my_MLJ_enf", shared=true)
Pkg.develop("path/to/my/local/MLJBase")
```

To test your local clone, do

```julia
Pkg.test("MLJBase")
```

For more on package management, see [here](https://julialang.github.io/Pkg.jl/v1/).

# [Internals](@id internals_section)

## The machine interface, simplified

The following is simplified description of the `Machine` interface. It
predates the introduction of an optional data front-end for models
(see [Implementing a data front-end](@ref)). See also the
[Glossary](glossary.md)


### The Machine type

````julia
mutable struct Machine{M<Model}

    model::M
    fitresult
    cache
    args::Tuple    # e.g., (X, y) for supervised models
    report
    previous_rows # remember last rows used

    function Machine{M}(model::M, args...) where M<:Model
        machine = new{M}(model)
        machine.args = args
        machine.report = Dict{Symbol,Any}()
        return machine
    end

end
````

### Constructor

````julia
machine(model::M, Xtable, y) = Machine{M}(model, Xtable, y)
````

### fit! and predict/transform

````julia
function fit!(mach::Machine; rows=nothing, force=false, verbosity=1)

    warning = clean!(mach.model)
    isempty(warning) || verbosity < 0 || @warn warning

    if rows === nothing
        rows = (:)
    end

    rows_have_changed  = (!isdefined(mach, :previous_rows) ||
	    rows != mach.previous_rows)

    args = [MLJ.selectrows(arg, rows) for arg in mach.args]

    if !isdefined(mach, :fitresult) || rows_have_changed || force
        mach.fitresult, mach.cache, report =
            fit(mach.model, verbosity, args...)
    else # call `update`:
        mach.fitresult, mach.cache, report =
            update(mach.model, verbosity, mach.fitresult, mach.cache, args...)
    end

    if rows_have_changed
        mach.previous_rows = deepcopy(rows)
    end

    if report !== nothing
        merge!(mach.report, report)
    end

    return mach

end

function predict(machine::Machine{<:Supervised}, Xnew)
    if isdefined(machine, :fitresult)
        return predict(machine.model, machine.fitresult, Xnew))
    else
        throw(error("$machine is not trained and so cannot predict."))
    end
end

function transform(machine::Machine{<:Unsupervised}, Xnew)
    if isdefined(machine, :fitresult)
        return transform(machine.model, machine.fitresult, Xnew))
    else
        throw(error("$machine is not trained and so cannot transform."))
    end
end
````
!!! warning "Old post"

    This post is quite old. For a newer overview of the design of MLJ, see [here](https://github.com/alan-turing-institute/MLJ.jl/blob/master/paper/paper.md)


# Beyond machine learning pipelines with MLJ

Anthony Blaom, Diego Arenas, Franz Kiraly, Yiannis Simillides, Sebastian Vollmer

**May 1st, 2019.** Blog post also posted on the [Julia Language Blog](https://julialang.org/blog/2019/05/beyond-ml-pipelines-with-mlj)




![](img/learningcurves.png) | ![](img/heatmap.png)
------------------------|--------------------------
![](img/wrapped_ridge.png)  | ![](img/MLPackages.png)


## Introducing MLJ

[MLJ](https://github.com/alan-turing-institute/MLJ.jl) is an
open-source machine learning toolbox written in pure Julia. It
provides a uniform interface for interacting with supervised and
unsupervised learning models currently scattered in different Julia
packages.

Building on a earlier proof-of-concept, development began in earnest
at [The Alan Turing Institute](https://www.turing.ac.uk) in
December 2018. In a short time interest grew and the project is now
the Institute's most starred software repository.

After outlining MLJ's current functionality, this post introduces MLJ
**learning networks**, a super-charged pipelining feature for model
composition.

**Quick links:**

- [MLJ vs ScikitLearn.jl](https://alan-turing-institute.github.io/MLJ.jl/dev/frequently_asked_questions/)  

- Video from [London Julia User Group meetup in March 2019](https://www.youtube.com/watch?v=CfHkjNmj1eE) (skip to [demo at 21'39](https://youtu.be/CfHkjNmj1eE?t=21m39s)) &nbsp;

- [MLJ Tutorials](https://JuliaAI.github.io/MLJTutorials/)

- Implementing the MLJ interface for a [new model](https://alan-turing-institute.github.io/MLJ.jl/dev/adding_models_for_general_use/)

- How to [contribute](https://github.com/alan-turing-institute/MLJ.jl/blob/master/CONTRIBUTE.md)

- Julia [Slack](http://julialang.slack.com) channel: \#mlj.

- Star'ing us to show support for [MLJ](https://github.com/alan-turing-institute/MLJ.jl) would be greatly appreciated!


## MLJ features

MLJ already has substantial functionality:

- **Learning networks.** Flexible model composition beyond traditional
  pipelines (more on this below).

- **Automatic tuning.** Automated tuning of hyperparameters, including
  composite models. Tuning implemented as a model wrapper for
  composition with other meta-algorithms.

- **Homogeneous model ensembling.**

- **Registry for model metadata.** Metadata available without loading
  model code. Basis of a "task" interface and facilitates
  model composition.

- **Task interface.** Automatically match models to specified learning
  tasks, to streamline benchmarking and model selection.

- **Clean probabilistic API.** Improves support for Bayesian
  statistics and probabilistic graphical models.

- **Data container agnostic.** Present and manipulate data in your
  favorite Tables.jl format.

- **Universal adoption of categorical data types.** Enables model
  implementations to properly account for classes seen in training but
  not in evaluation.

Enhancements planned for the near future include integration of
Flux.jl **deep learning** models, and **gradient descent tuning** of
continuous hyperparameters using automatic differentiation.

While a relatively small number of machine learning models currently
implement the MLJ interface, work in progress aims to wrap models
supported by the popular python framework, scikit-learn, as a
temporary expedient. For a comparison of the MLJ's design with the
Julia wrap [ScitLearn.jl](https://github.com/cstjean/ScikitLearn.jl),
see this
[FAQ](https://github.com/alan-turing-institute/MLJ.jl/blob/master/docs/src/frequently_asked_questions.md).


## Learning networks

MLJ's model composition interface is flexible enough to implement, for
example, the [model
stacks](https://www.kdnuggets.com/2017/02/stacking-models-imropved-predictions.html)
popular in data science competitions. To treat examples of this kind,
the interface design must account for the fact that information flow
in prediction and training modes is different. This can be seen from
the following schematic of a simple two-model stack, viewed as a
network:

![](img/two_model_stack.png)

## Building a simple network

In MLJ, networks of models are built using a declarative syntax
already familiar from basic use of the package. For example, the
ordinary syntax for training a decision tree in MLJ, after one-hot
encoding the categorical features, looks like this:

```julia
using MLJ
@load DecisionTreeRegressor

# load some data:
task = load_reduced_ames();
X, y = task();

# one-hot encode the inputs, X:
hot_model = OneHotEncoder()
hot = machine(hot_model, X)
fit!(hot)
Xt = transform(hot, X)

# fit a decision tree to the transformed data:
tree_model = DecisionTreeRegressor()
tree = machine(tree_model, Xt, y)
fit!(tree, rows = 1:1300)
```

Note that a *model* in MLJ is just a struct containing
hyperparameters. Wrapping a model in data delivers a *machine* struct,
which will additionally record the results of training.

Without a pipeline, each time we want to present new data for
prediction we must first apply one-hot encoding:

```julia
Xnew = X[1301:1400,:];
Xnewt = transform(hot, Xnew);
yhat = predict(tree, Xnewt);
yhat[1:3]
 3-element Array{Float64,1}:
  223956.9999999999
  320142.85714285733
  161227.49999999994
```

To build a pipeline one simply wraps the supplied data in source nodes
and repeats similar declarations, omitting calls to
`fit!`. The difference now is that each "variable" (e.g., `Xt`,
`yhat`) is a node of our pipeline, instead of concrete data:

```julia
Xs = source(X)
ys = source(y)

hot = machine(hot_model, Xs)
Xt = transform(hot, Xs);

tree = machine(tree_model, Xt, ys)
yhat = predict(tree, Xt)
```

If we like, we can think of a node as *dynamic data* - "data" because
it can be called (indexed) on rows, but "dynamic" because the result
depends on the outcome of training events, which in turn depend on
hyperparameter values. For example, after fitting the completed pipeline,
we can make new predictions like this:

```julia
fit!(yhat, rows=1:1300)
 [ Info: Training NodalMachine @ 1…51.
 [ Info: Spawned 1300 sub-features to one-hot encode feature :Neighborhood.
 [ Info: Spawned 1300 sub-features to one-hot encode feature :MSSubClass.
 [ Info: Training NodalMachine @ 1…17.
 Node @ 1…79 = predict(1…17, transform(1…51, 1…07))

yhat(rows=1301:1302) # to predict on rows of source node
yhat(Xnew)           # to predict on new data
156-element Array{Float64,1}:
 223956.9999999999
 320142.85714285733
 ...
```


## Exporting and retraining

Once a pipeline like this has been built and tested on sample data, it
can be exported as a stand-alone model, ready to be trained on any
dataset. For details, see the MLJ
[documentation](https://alan-turing-institute.github.io/MLJ.jl/dev/learning_networks/). In
the future, Julia macros will allow common architectures (e.g., linear
pipelines) to be built in a couple of lines.

Finally, we mention that MLJ learning networks, and their exported
counterparts, are "smart" in the sense that changing a hyperparameter
does not trigger retraining of component models upstream of the
change:

```julia
tree_model.max_depth = 4
fit!(yhat, rows=1:1300)
 [ Info: Not retraining NodalMachine @ 1…51. It is up-to-date.
 [ Info: Updating NodalMachine @ 1…17.
 Node @ 1…79 = predict(1…17, transform(1…51, 1…07))
```


## Just "Write the math!"

Because of Julia's generic programming features, any kind of operation
you would normally apply to data (arithmetic, row selection, column
concatenation, etc) can be overloaded to work with nodes. In this way,
MLJ's network-building syntax is economical, intuitive and easy to
read. In this respect we have been inspired by [On Machine Learning
and Programming Languages](https://julialang.org/blog/2017/12/ml&pl).

## Invitation to the community
We now invite the community to try out our newly registered packages, [MLJ](https://github.com/alan-turing-institute/MLJ.jl)alongside [MLJModels](https://github.com/JuliaAI/MLJModels.jl), and provide any feedback or suggestions you may have going forward. We are also particularly interested in hearing how you would use our package, and what features it may be lacking.
# Getting Started

For an outline of MLJ's **goals** and **features**, see
[About MLJ](@ref).

This section introduces the most basic MLJ operations and concepts. It
assumes MJL has been successfully installed. See [Installation](@ref)
if this is not the case.


```@setup doda
import Random.seed!
using MLJ
using InteractiveUtils
MLJ.color_off()
seed!(1234)
```

## Choosing and evaluating a model

The following code loads Fisher's famous iris data set as a named tuple of
column vectors:

```@repl doda
using MLJ
iris = load_iris();
selectrows(iris, 1:3)  |> pretty
schema(iris)
```

Because this data format is compatible with
[Tables.jl](https://tables.juliadata.org/stable/), many MLJ methods
(such as `selectrows`, `pretty` and `schema` used above) as well as
many MLJ models, can work with it. However, as most new users are
already familiar with the access methods particular to
[DataFrames](https://dataframes.juliadata.org/stable/) (also
compatible with Tables.jl) we'll put our data into that format here:

```@example doda
import DataFrames
iris = DataFrames.DataFrame(iris);
nothing # hide
```

Next, let's split the data "horizontally" into input and target parts,
and specify an RNG seed, to force observations to be shuffled:

```@repl doda
y, X = unpack(iris, ==(:target); rng=123);
first(X, 3) |> pretty
```

This call to `unpack` splits off any column with name `==`
to `:target` into something called `y`, and all the remaining columns
into `X`.

To list *all* models available in MLJ's [model
registry](@ref model_search) do `models()`. Listing the models
compatible with the present data:

```@repl doda
models(matching(X,y))
```

In MLJ a *model* is a struct storing the hyperparameters of the
learning algorithm indicated by the struct name (and nothing
else). For common problems matching data to models, see [Model
Search](@ref model_search) and [Preparing Data](@ref).

Assuming the MLJDecisionTreeInterface.jl package is in your load path
(see [Installation](@ref)) we can use `@load` to import the
`DecisionTreeClassifier` model type, which we will bind to `Tree`:

```@repl doda
Tree = @load DecisionTreeClassifier pkg=DecisionTree
```

(In this case we need to specify `pkg=...` because multiple packages
provide a model type with name `DecisionTreeClassifier`.) Now we can
instantiate a model with default hyperparameters:

```@repl doda
tree = Tree()
```

*Important:* DecisionTree.jl and most other packages implementing
machine learning algorithms for use in MLJ are not MLJ
dependencies. If such a package is not in your load path you will
receive an error explaining how to add the package to your current
environment. Alternatively, you can use the interactive macro
`@iload`. For more on importing model types, see [Loading Model
Code](@ref).

Once instantiated, a model's performance can be evaluated with the
`evaluate` method. Our classifier is a *probabilistic* predictor (check
`prediction_type(tree) == :probabilistic`) which means we can specify
a probabilistic measure (metric) like `log_loss`, as well
deterministic measures like `accuracy` (which are applied after
computing the mode of each prediction):

```@repl doda
evaluate(tree, X, y,
         resampling=CV(shuffle=true),
                 measures=[log_loss, accuracy],
                 verbosity=0)
```

Under the hood, `evaluate` calls lower level functions `predict` or
`predict_mode` according to the type of measure, as shown in the
output. We shall call these operations directly below.

For more on performance evaluation, see [Evaluating Model
Performance](evaluating_model_performance.md) for details.


## A preview of data type specification in MLJ

The target `y` above is a categorical vector, which is appropriate
because our model is a decision tree *classifier*:

```@repl doda
typeof(y)
```

However, MLJ models do not actually prescribe the machine types for
the data they operate on. Rather, they specify a *scientific type*,
which refers to the way data is to be *interpreted*, as opposed to how
it is *encoded*:

```@repl doda
target_scitype(tree)
```

Here `Finite` is an example of a "scalar" scientific type with two
subtypes:

```@repl doda
subtypes(Finite)
```

We use the `scitype` function to check how MLJ is going to interpret
given data. Our choice of encoding for `y` works for
`DecisionTreeClassifier`, because we have:

```@repl doda
scitype(y)
```

and `Multiclass{3} <: Finite`. If we would encode with integers
instead, we obtain:

```@repl doda
yint = int.(y);
scitype(yint)
```

and using `yint` in place of `y` in classification problems will
fail. See also [Working with Categorical Data](@ref).

For more on scientific types, see [Data containers and scientific
types](@ref) below.


## Fit and predict

To illustrate MLJ's fit and predict interface, let's perform our
performance evaluations by hand, but using a simple holdout set,
instead of cross-validation.

Wrapping the model in data creates a *machine* which will store
training outcomes:

```@repl doda
mach = machine(tree, X, y)
```

Training and testing on a hold-out set:

```@repl doda
train, test = partition(eachindex(y), 0.7); # 70:30 split
fit!(mach, rows=train);
yhat = predict(mach, X[test,:]);
yhat[3:5]
log_loss(yhat, y[test]) |> mean
```

Note that `log_loss` and `cross_entropy` are aliases for `LogLoss()`
(which can be passed an optional keyword parameter, as in
`LogLoss(tol=0.001)`). For a list of all losses and scores, and their
aliases, run `measures()`.

Notice that `yhat` is a vector of `Distribution` objects, because
DecisionTreeClassifier makes probabilistic predictions. The methods
of the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
package can be applied to such distributions:

```@repl doda
broadcast(pdf, yhat[3:5], "virginica") # predicted probabilities of virginica
broadcast(pdf, yhat, y[test])[3:5] # predicted probability of observed class
mode.(yhat[3:5])
```

Or, one can explicitly get modes by using `predict_mode` instead of
`predict`:

```@repl doda
predict_mode(mach, X[test[3:5],:])
```

Finally, we note that `pdf()` is overloaded to allow the retrieval of
probabilities for all levels at once:

```@repl doda
L = levels(y)
pdf(yhat[3:5], L)
```

Unsupervised models have a `transform` method instead of `predict`,
and may optionally implement an `inverse_transform` method:

```@repl doda
v = Float64[1, 2, 3, 4]
stand = Standardizer() # this type is built-in
mach2 = machine(stand, v)
fit!(mach2)
w = transform(mach2, v)
inverse_transform(mach2, w)
```

[Machines](machines.md) have an internal state which allows them to
avoid redundant calculations when retrained, in certain conditions -
for example when increasing the number of trees in a random forest, or
the number of epochs in a neural network. The machine building syntax
also anticipates a more general syntax for composing multiple models,
an advanced feature explained in [Learning Networks](@ref).

There is a version of `evaluate` for machines as well as models. This
time we'll use a simple holdout strategy as above. (An exclamation
point is added to the method name because machines are generally
mutated when trained.)

```@repl doda
evaluate!(mach, resampling=Holdout(fraction_train=0.7),
                measures=[log_loss, accuracy],
                verbosity=0)
```
Changing a hyperparameter and re-evaluating:

```@repl doda
tree.max_depth = 3
evaluate!(mach, resampling=Holdout(fraction_train=0.7),
          measures=[log_loss, accuracy],
          verbosity=0)
```

## Next steps

To learn a little more about what MLJ can do, browse [Common MLJ
Workflows](common_mlj_workflows.md) or [Data Science Tutorials in
Julia](https://alan-turing-institute.github.io/DataScienceTutorials.jl/)
or try the [JuliaCon2020
Workshop](https://github.com/ablaom/MachineLearningInJulia2020) on MLJ
(recorded
[here](https://www.youtube.com/watch?time_continue=27&v=qSWbCn170HU&feature=emb_title))
returning to the manual as needed.

*Read at least the remainder of this page before considering serious
use of MLJ.*


## Data containers and scientific types

The MLJ user should acquaint themselves with some basic assumptions
about the form of data expected by MLJ, as outlined below. The basic
`machine` constructors look like this (see also [Constructing
machines](@ref)):

```
machine(model::Unsupervised, X)
machine(model::Supervised, X, y)
```

Each supervised model in MLJ declares the permitted *scientific type*
of the inputs `X` and targets `y` that can be bound to it in the first
constructor above, rather than specifying specific machine types (such
as `Array{Float32, 2}`). Similar remarks apply to the input `X` of an
unsupervised model.

Scientific types are julia types defined in the package
[ScientificTypesBase.jl](https://github.com/JuliaAI/ScientificTypesBase.jl);
the package
[ScientificTypes.jl](https://JuliaAI.github.io/ScientificTypes.jl/dev/)
implements the particular convention used in the MLJ universe for
assigning a specific scientific type (interpretation) to each julia
object (see the `scitype` examples below).

The basic "scalar" scientific types are `Continuous`, `Multiclass{N}`,
`OrderedFactor{N}`, `Count` and `Textual`. `Missing` and `Nothing` are
also considered scientific types. Be sure you read [Scalar scientific
types](@ref) below to guarantee your scalar data is interpreted
correctly. Tools exist to coerce the data to have the appropriate
scientific type; see
[ScientificTypes.jl](https://JuliaAI.github.io/ScientificTypes.jl/dev/)
or run `?coerce` for details.

Additionally, most data containers - such as tuples, vectors, matrices
and tables - have a scientific type parameterized by scitype of the
elements they contain.

![](img/scitypes_small.png)

*Figure 1. Part of the scientific type hierarchy in*
[ScientificTypesBase.jl](https://JuliaAI.github.io/ScientificTypes.jl/dev/).

```@repl doda
scitype(4.6)
scitype(42)
x1 = coerce(["yes", "no", "yes", "maybe"], Multiclass);
scitype(x1)
X = (x1=x1, x2=rand(4), x3=rand(4))  # a "column table"
scitype(X)
```

### Two-dimensional data

Generally, two-dimensional data in MLJ is expected to be *tabular*.
All data containers compatible with the
[Tables.jl](https://github.com/JuliaData/Tables.jl) interface (which
includes all source formats listed
[here](https://github.com/JuliaData/Tables.jl/blob/master/INTEGRATIONS.md))
have the scientific type `Table{K}`, where `K` depends on the
scientific types of the columns, which can be individually inspected
using `schema`:

```@repl doda
schema(X)
```

#### Matrix data

MLJ models expecting a table do not generally accept a matrix
instead. However, a matrix can be wrapped as a table, using
[`MLJ.table`](@ref):

```julia
matrix_table = MLJ.table(rand(2,3))
schema(matrix_table)
```

```
┌─────────┬─────────┬────────────┐
│ _.names │ _.types │ _.scitypes │
├─────────┼─────────┼────────────┤
│ x1      │ Float64 │ Continuous │
│ x2      │ Float64 │ Continuous │
│ x3      │ Float64 │ Continuous │
└─────────┴─────────┴────────────┘
_.nrows = 2

```

The matrix is *not* copied, only wrapped.  To manifest a table as a
matrix, use [`MLJ.matrix`](@ref).

### Observations correspond to rows, not columns

When supplying models with matrices, or wrapping them in tables, each
*row* should correspond to a different observation. That is, the
matrix should be `n x p`, where `n` is the number of observations and
`p` the number of features. However, *some models may perform better* if
supplied the *adjoint* of a `p x n` matrix instead, and observation
resampling is always more efficient in this case.


### Inputs

Since an MLJ model only specifies the scientific type of data, if that
type is `Table` - which is the case for the majority of MLJ models -
then any [Tables.jl](https://github.com/JuliaData/Tables.jl) format is
permitted.

Specifically, the requirement for an arbitrary model's input is `scitype(X)
<: input_scitype(model)`.

### Targets

The target `y` expected by MLJ models is generally an
`AbstractVector`. A multivariate target `y` will generally be a table.

Specifically, the type requirement for a model target is `scitype(y) <:
target_scitype(model)`.


### Querying a model for acceptable data types

Given a model instance, one can inspect the admissible scientific
types of its input and target, and without loading the code defining
the model;

```@setup doda
tree = @load DecisionTreeClassifier pkg=DecisionTree
```

```@repl doda
i = info("DecisionTreeClassifier", pkg="DecisionTree")
i.input_scitype
i.target_scitype
```

This output indicates that any table with `Continuous`, `Count` or
`OrderedFactor` columns is acceptable as the input `X`, and that any
vector with element scitype `<: Finite` is acceptable as the target
`y`.

For more on matching models to data, see [Model Search](@ref model_search).

### Scalar scientific types

Models in MLJ will always apply the `MLJ` convention described in
[ScientificTypes.jl](https://JuliaAI.github.io/ScientificTypes.jl/dev/)
to decide how to interpret the elements of your container types. Here
are the key features of that convention:

- Any `AbstractFloat` is interpreted as `Continuous`.

- Any `Integer` is interpreted as `Count`.

- Any `CategoricalValue` `x`, is interpreted as `Multiclass` or
  `OrderedFactor`, depending on the value of `isordered(x)`.

- `String`s and `Char`s are *not* interpreted as `Multiclass` or
  `OrderedFactor` (they have scitypes `Textual` and `Unknown`
  respectively).

- In particular, *integers* (including `Bool`s) *cannot be used to
  represent categorical data.* Use the preceding `coerce` operations
  to coerce to a `Finite` scitype.

- The scientific types of `nothing` and `missing` are `Nothing` and
  `Missing`, native types we also regard as scientific.

Use `coerce(v, OrderedFactor)` or `coerce(v, Multiclass)` to coerce a
vector `v` of integers, strings or characters to a vector with an
appropriate `Finite` (categorical) scitype.  See [Working with
Categorical Data](@ref).

For more on scitype coercion of arrays and tables, see [`coerce`](@ref),
[`autotype`](@ref) and [`unpack`](@ref) below and the examples at
[ScientificTypes.jl](https://JuliaAI.github.io/ScientificTypes.jl/dev/).



```@docs
scitype
coerce
autotype
```
# About MLJ

MLJ (Machine Learning in Julia) is a toolbox written in Julia
providing a common interface and meta-algorithms for selecting,
tuning, evaluating, composing and comparing [over 150 machine learning
models](@ref model_list) written in Julia and other languages. In
particular MLJ wraps a large number of
[scikit-learn](https://scikit-learn.org/stable/) models.

MLJ is released under the MIT licensed.

## Lightning tour

*For more elementary introductions to MLJ usage see [Basic
introductions](@ref) below.*

A self-contained notebook and julia script of this demonstration is
also available
[here](https://github.com/alan-turing-institute/MLJ.jl/tree/dev/examples/lightning_tour).

The first code snippet below creates a new Julia environment
`MLJ_tour` and installs just those packages needed for the tour. See
[Installation](@ref) for more on creating a Julia environment for use
with MLJ.

Julia installation instructions are
[here](https://julialang.org/downloads/).

```julia
using Pkg
Pkg.activate("MLJ_tour", shared=true)
Pkg.add("MLJ")
Pkg.add("MLJIteration")
Pkg.add("EvoTrees")
```

In MLJ a *model* is just a container for hyper-parameters, and that's
all. Here we will apply several kinds of model composition before
binding the resulting "meta-model" to data in a *machine* for
evaluation using cross-validation.

Loading and instantiating a gradient tree-boosting model:

```julia
using MLJ
Booster = @load EvoTreeRegressor # loads code defining a model type
booster = Booster(max_depth=2)   # specify hyper-parameter at construction
booster.nrounds=50               # or mutate afterwards
```

This model is an example of an iterative model. As is stands, the
number of iterations `nrounds` is fixed.


#### Composition 1: Wrapping the model to make it "self-iterating"

Let's create a new model that automatically learns the number of iterations,
using the `NumberSinceBest(3)` criterion, as applied to an
out-of-sample `l1` loss:

```julia
using MLJIteration
iterated_booster = IteratedModel(model=booster,
                                 resampling=Holdout(fraction_train=0.8),
                                 controls=[Step(2), NumberSinceBest(3), NumberLimit(300)],
                                 measure=l1,
                                 retrain=true)
```

#### Composition 2: Preprocess the input features

Combining the model with categorical feature encoding:

```julia
pipe = ContinuousEncoder() |> iterated_booster
```

#### Composition 3: Wrapping the model to make it "self-tuning"

First, we define a hyper-parameter range for optimization of a
(nested) hyper-parameter:

```julia
max_depth_range = range(pipe,
                        :(deterministic_iterated_model.model.max_depth),
                        lower = 1,
                        upper = 10)
```

Now we can wrap the pipeline model in an optimization strategy to make
it "self-tuning":

```julia
self_tuning_pipe = TunedModel(model=pipe,
                              tuning=RandomSearch(),
                              ranges = max_depth_range,
                              resampling=CV(nfolds=3, rng=456),
                              measure=l1,
                              acceleration=CPUThreads(),
                              n=50)
```

#### Binding to data and evaluating performance

Loading a selection of features and labels from the Ames
House Price dataset:

```julia
X, y = @load_reduced_ames;
```
Evaluating the "self-tuning" pipeline model's performance using 5-fold
cross-validation (implies multiple layers of nested resampling):

```julia
julia> evaluate(self_tuning_pipe, X, y,
                measures=[l1, l2],
                resampling=CV(nfolds=5, rng=123),
                acceleration=CPUThreads(),
                verbosity=2)
PerformanceEvaluation object with these fields:
  measure, measurement, operation, per_fold,
  per_observation, fitted_params_per_fold,
  report_per_fold, train_test_pairs
Extract:
┌───────────────┬─────────────┬───────────┬───────────────────────────────────────────────┐
│ measure       │ measurement │ operation │ per_fold                                      │
├───────────────┼─────────────┼───────────┼───────────────────────────────────────────────┤
│ LPLoss(p = 1) │ 17200.0     │ predict   │ [16500.0, 17100.0, 16300.0, 17500.0, 18900.0] │
│ LPLoss(p = 2) │ 6.83e8      │ predict   │ [6.14e8, 6.64e8, 5.98e8, 6.37e8, 9.03e8]      │
└───────────────┴─────────────┴───────────┴───────────────────────────────────────────────┘
```

Try out MLJ yourself in the following batteries-included Binder
[notebook](https://mybinder.org/v2/gh/alan-turing-institute/MLJ.jl/master?filepath=binder%2FMLJ_demo.ipynb). No
installation required.


## Key goals

* Offer a consistent way to use, compose and tune machine learning
  models in Julia,

* Promote the improvement of the Julia ML/Stats ecosystem by making it
  easier to use models from a wide range of packages,

* Unlock performance gains by exploiting Julia's support for
  parallelism, automatic differentiation, GPU, optimization etc.


## Key features

* Data agnostic, train models on any data supported by the
  [Tables.jl](https://github.com/JuliaData/Tables.jl) interface.

* Extensive, state-of-the art, support for model composition
  (*pipelines*, *stacks* and, more generally, *learning networks*). See more
  [below](#model-composability).

* Convenient syntax to tune and evaluate (composite) models.

* Consistent interface to handle probabilistic predictions.

* Extensible [tuning
  interface](https://github.com/JuliaAI/MLJTuning.jl),
  to support growing number of optimization strategies, and designed
  to play well with model composition.

* Options to accelerate model evaluation and tuning with
  multithreading and/or distributed processing.


## Model composability

The generic model composition API's provided by other toolboxes we
have surveyed share one or more of the following shortcomings, which
do not exist in MLJ:

- Composite models do not inherit all the behavior of ordinary
  models.

- Composition is limited to linear (non-branching) pipelines.

- Supervised components in a linear pipeline can only occur at the
  end of the pipeline.

- Only static (unlearned) target transformations/inverse
  transformations are supported.

- Hyper-parameters in homogeneous model ensembles cannot be coupled.

- Model stacking, with out-of-sample predictions for base learners,
  cannot be implemented (using the generic API alone).

- Hyper-parameters and/or learned parameters of component models are
  not easily inspected or manipulated (by tuning algorithms, for
  example)

- Composite models cannot implement multiple operations, for example,
  both a `predict` and `transform` method (as in clustering models) or
  both a `transform` and `inverse_transform` method.

Some of these features are demonstrated in [this
notebook](https://github.com/ablaom/MachineLearningInJulia2020/blob/master/wow.ipynb)

For more information see the [MLJ design
paper](https://doi.org/10.21105/joss.02704) or our detailed
[paper](https://arxiv.org/abs/2012.15505) on the composition
interface.


## Getting help and reporting problems

Users are encouraged to provide feedback on their experience using MLJ
and to report issues.

For a query to have maximum exposure to maintainers and users, start a
discussion thread at [Julia Discourse Machine
Learning](https://github.com/alan-turing-institute/MLJ.jl) and tag
your issue "mlj". Queries can also be posted as
[issues](https://github.com/alan-turing-institute/MLJ.jl/issues), or
on the `#mlj` slack workspace in the Julia Slack channel.

Bugs, suggestions, and feature requests can be posted
[here](https://github.com/alan-turing-institute/MLJ.jl/issues).

See also, [Known Issues](@ref)


## Installation

Initially it is recommended that MLJ and associated packages be
installed in a new
[environment](https://julialang.github.io/Pkg.jl/v1/environments/) to
avoid package conflicts. You can do this with

```julia
julia> using Pkg; Pkg.activate("my_MLJ_env", shared=true)
```

Installing MLJ is also done with the package manager:

```julia
julia> Pkg.add("MLJ")
```

**Optional:** To test your installation, run

```julia
julia> Pkg.test("MLJ")
```

It is important to note that MLJ is essentially a big wrapper
providing unified access to _model providing packages_. For this
reason, one generally needs to add further packages to your
environment to make model-specific code available. This
happens automatically when you use MLJ's interactive load command
`@iload`, as in

```julia
julia> Tree = @iload DecisionTreeClassifier # load type
julia> tree = Tree() # instance
```

where you will also be asked to choose a providing package, for more
than one provide a `DecisionTreeClassifier` model. For more on
identifying the name of an applicable model, see [Model Search](@ref model_search).
For non-interactive loading of code (e.g., from a
module or function) see [Loading Model Code](@ref).

It is recommended that you start with models from more mature
packages such as DecisionTree.jl, ScikitLearn.jl or XGBoost.jl.

MLJ is supported by a number of satellite packages (MLJTuning,
MLJModelInterface, etc) which the general user is *not* required to
install directly. Developers can learn more about these
[here](https://github.com/alan-turing-institute/MLJ.jl/blob/master/ORGANIZATION.md).

See also the alternative instalation instructions for [Modifying Behavior](@ref).


## Learning Julia

If you have experience in programming in another language but are new
to Julia, then we highly recommend Aaron Christinson's tutorial
[Dispatching Design
Patterns](https://github.com/ninjaaron/dispatching-design-patterns)
which is nicely compressed in his [half-hour
video presentation](https://www.youtube.com/watch?v=n-E-1-A_rZM).

However, one doesn't need to be able to program in Julia to start
using MLJ.

## Learning to use MLJ

The present document, although littered with examples, is primarily
intended as a complete reference. Resources for learning MLJ are:

### Basic introductions

- the [Getting Started](@ref) section of this manual

- an introductory [binder notebook](https://mybinder.org/v2/gh/alan-turing-institute/MLJ.jl/master?filepath=binder%2FMLJ_demo.ipynb) (no Julia/MLJ installation required)

To get direct help from maintainers and other users, see [Getting help
and reporting problems](@ref).


### In depth

- the MLJ JuliaCon2020 Workshop [materials](https://github.com/ablaom/MachineLearningInJulia2020) and [video recording](https://www.youtube.com/watch?time_continue=27&v=qSWbCn170HU&feature=emb_title)

- [Data Science Tutorials in Julia](https://alan-turing-institute.github.io/DataScienceTutorials.jl/)

Users are also welcome to join the `#mlj` Julia slack channel to ask
questions and make suggestions.

## Funding

MLJ was initially created as a Tools, Practices and Systems project at
the [Alan Turing Institute](https://www.turing.ac.uk/)
in 2019. Current funding is provided by a [New Zealand Strategic
Science Investment
Fund](https://www.mbie.govt.nz/science-and-technology/science-and-innovation/funding-information-and-opportunities/investment-funds/strategic-science-investment-fund/ssif-funded-programmes/university-of-auckland/)
awarded to the University of Auckland.


## Citing MLJ

An overview of MLJ design:


[![DOI](https://joss.theoj.org/papers/10.21105/joss.02704/status.svg)](https://doi.org/10.21105/joss.02704)

```bibtex
@article{Blaom2020,
  doi = {10.21105/joss.02704},
  url = {https://doi.org/10.21105/joss.02704},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {55},
  pages = {2704},
  author = {Anthony D. Blaom and Franz Kiraly and Thibaut Lienart and Yiannis Simillides and Diego Arenas and Sebastian J. Vollmer},
  title = {{MLJ}: A Julia package for composable machine learning},
  journal = {Journal of Open Source Software}
}
```

An in-depth view of MLJ's model composition design:

[![arXiv](https://img.shields.io/badge/arXiv-2012.15505-<COLOR>.svg)](https://arxiv.org/abs/2012.15505)

```bibtex
@misc{blaom2020flexible,
      title={Flexible model composition in machine learning and its implementation in {MLJ}},
      author={Anthony D. Blaom and Sebastian J. Vollmer},
      year={2020},
      eprint={2012.15505},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```
# Controlling Iterative Models

Iterative supervised machine learning models are usually trained until
an out-of-sample estimate of the performance satisfies some stopping
criterion, such as `k` consecutive deteriorations of the performance
(see [`Patience`](@ref EarlyStopping.Patience) below). A more
sophisticated kind of control might dynamically mutate parameters,
such as a learning rate, in response to the behavior of these
estimates.

Some iterative model implementations enable some form of automated
control, with the method and options for doing so varying from model
to model. But sometimes it is up to the user to arrange control, which
in the crudest case reduces to manually experimenting with the
iteration parameter.

In response to this ad hoc state of affairs, MLJ provides a uniform
and feature-rich interface for controlling any iterative model that
exposes its iteration parameter as a hyper-parameter, and which
implements the "warm restart" behavior described in [Machines](@ref).


## Basic use

As in [Tuning Models](@ref), iteration control in MLJ is implemented as
a model wrapper, which allows composition with other meta-algorithms.
Ordinarily, the wrapped model behaves just like the original model,
but with the training occurring on a subset of the provided data (to
allow computation of an out-of-sample loss) and with the iteration
parameter automatically determined by the controls specified in the
wrapper.

By setting `retrain=true` one can ask that the wrapped model retrain
on *all* supplied data, after learning the appropriate number of
iterations from the controlled training phase:

```@example gree
using MLJ

X, y = make_moons(100, rng=123, noise=0.5)
EvoTreeClassifier = @load EvoTreeClassifier verbosity=0

iterated_model = IteratedModel(model=EvoTreeClassifier(rng=123, η=0.005),
                               resampling=Holdout(),
                               measures=log_loss,
                               controls=[Step(5),
                                         Patience(2),
                                         NumberLimit(100)],
                               retrain=true)

mach = machine(iterated_model, X, y)
nothing # hide
```

```@repl gree
fit!(mach)
```

As detailed under [`IteratedModel`](@ref MLJIteration.IteratedModel)
below, the specified `controls` are repeatedly applied in sequence to
a *training machine*, constructed under the hood, until one of the
controls triggers a stop. Here `Step(5)` means "Compute 5 more
iterations" (in this case starting from none); `Patience(2)` means
"Stop at the end of the control cycle if there have been 2 consecutive
drops in the log loss"; and `NumberLimit(100)` is a safeguard ensuring
a stop after 100 control cycles (500 iterations). See [Controls
provided](@ref) below for a complete list.

Because iteration is implemented as a wrapper, the "self-iterating"
model can be evaluated using cross-validation, say, and the number of
iterations on each fold will generally be different:

```@example gree
e = evaluate!(mach, resampling=CV(nfolds=3), measure=log_loss, verbosity=0);
map(e.report_per_fold) do r
    r.n_iterations
end
```

Alternatively, one might wrap the self-iterating model in a tuning
strategy, using `TunedModel`; see [Tuning Models](@ref). In this way,
the optimization of some other hyper-parameter is realized
simultaneously with that of the iteration parameter, which will
frequently be more efficient than a direct two-parameter search.


## Controls provided

In the table below, `mach` is the *training machine* being iterated,
constructed by binding the supplied data to the `model` specified in
the `IteratedModel` wrapper, but trained in each iteration on a subset
of the data, according to the value of the `resampling`
hyper-parameter of the wrapper (using all data if `resampling=nothing`).

control                                                        | description                                                                             | can trigger a stop
---------------------------------------------------------------|-----------------------------------------------------------------------------------------|--------------------
[`Step`](@ref IterationControl.Step)`(n=1)`                    | Train model for `n` more iterations                                                     | no
[`TimeLimit`](@ref EarlyStopping.TimeLimit)`(t=0.5)`           | Stop after `t` hours                                                                    | yes
[`NumberLimit`](@ref EarlyStopping.NumberLimit)`(n=100)`       | Stop after `n` applications of the control                                              | yes
[`NumberSinceBest`](@ref EarlyStopping.NumberSinceBest)`(n=6)` | Stop when best loss occurred `n` control applications ago                               | yes
[`InvalidValue`](@ref IterationControl.InvalidValue)()         | Stop when `NaN`, `Inf` or `-Inf` loss/training loss encountered                         | yes 
[`Threshold`](@ref EarlyStopping.Threshold)`(value=0.0)`       | Stop when `loss < value`                                                                | yes
[`GL`](@ref EarlyStopping.GL)`(alpha=2.0)`                     | † Stop after the "generalization loss (GL)" exceeds `alpha`                             | yes
[`PQ`](@ref EarlyStopping.PQ)`(alpha=0.75, k=5)`               | † Stop after "progress-modified GL" exceeds `alpha`                                     | yes
[`Patience`](@ref EarlyStopping.Patience)`(n=5)`               | † Stop after `n` consecutive loss increases                                             | yes
[`Info`](@ref IterationControl.Info)`(f=identity)`             | Log to `Info` the value of `f(mach)`, where `mach` is current machine                   | no
[`Warn`](@ref IterationControl.Warn)`(predicate; f="")`        | Log to `Warn` the value of `f` or `f(mach)`, if `predicate(mach)` holds                 | no
[`Error`](@ref IterationControl.Error)`(predicate; f="")`      | Log to `Error` the value of `f` or `f(mach)`, if `predicate(mach)` holds and then stop  | yes
[`Callback`](@ref IterationControl.Callback)`(f=mach->nothing)`| Call `f(mach)`                                                                          | yes
[`WithNumberDo`](@ref IterationControl.WithNumberDo)`(f=n->@info(n))`                       | Call `f(n + 1)` where `n` is the number of complete control cycles so far | yes
[`WithIterationsDo`](@ref MLJIteration.WithIterationsDo)`(f=i->@info("iterations: $i"))`| Call `f(i)`, where `i` is total number of iterations           | yes
[`WithLossDo`](@ref IterationControl.WithLossDo)`(f=x->@info("loss: $x"))`                  | Call `f(loss)` where `loss` is the current loss            | yes
[`WithTrainingLossesDo`](@ref IterationControl.WithTrainingLossesDo)`(f=v->@info(v))`       | Call `f(v)` where `v` is the current batch of training losses | yes
[`WithEvaluationDo`](@ref MLJIteration.WithEvaluationDo)`(f->e->@info("evaluation: $e))`| Call `f(e)` where `e` is the current performance evaluation object | yes
[`WithFittedParamsDo`](@ref MLJIteration.WithFittedParamsDo)`(f->fp->@info("fitted_params: $fp))`| Call `f(fp)` where `fp` is fitted parameters of training machine | yes
[`WithReportDo`](@ref MLJIteration.WithReportDo)`(f->e->@info("report: $e))`| Call `f(r)` where `r` is the training machine report                    | yes
[`WithModelDo`](@ref MLJIteration.WithModelDo)`(f->m->@info("model: $m))`| Call `f(m)` where `m` is the model, which may be mutated by `f`             | yes
[`WithMachineDo`](@ref MLJIteration.WithMachineDo)`(f->mach->@info("report: $mach))`| Call `f(mach)` wher `mach` is the training machine in its current state    | yes
[`Save`](@ref MLJSerialization.Save)`(filename="machine.jlso")`| ⋆ Save current training machine to `machine1.jlso`, `machine2.jslo`, etc                         | yes

> Table 1. Atomic controls. Some advanced options omitted.

† For more on these controls see [Prechelt, Lutz
 (1998)](https://link.springer.com/chapter/10.1007%2F3-540-49430-8_3):
 "Early Stopping - But When?", in *Neural Networks: Tricks of the
 Trade*, ed. G. Orr, Springer.

⋆ If using `MLJIteration` without `MLJ`, then `Save` is not available
  unless one is also using `MLJSerialization`.

**Stopping option.** All the following controls trigger a stop if the
provided function `f` returns `true` and `stop_if_true=true` is
specified in the constructor: `Callback`, `WithNumberDo`,
`WithLossDo`, `WithTrainingLossesDo`.

There are also three control wrappers to modify a control's behavior:

wrapper                                                                    | description
---------------------------------------------------------------------------|-------------------------------------------------------------------------
[`IterationControl.skip`](@ref)`(control, predicate=1)`                    | Apply `control` every `predicate` applications of the control wrapper (can also be a function; see doc-string)
[`IterationControl.louder`](@ref IterationControl.louder)`(control, by=1)` | Increase the verbosity level of `control` by the specified value (negative values lower verbosity)
[`IterationControl.with_state_do`](@ref)`(control; f=...)`                 | Apply control *and* call `f(x)` where `x` is the internal state of control; useful for debugging. Default `f` logs state to `Info`. **Warning**: internal control state is not yet part of public API.
[`IterationControl.composite`](@ref)`(controls...)`                        | Apply each `control` in `controls` in sequence; used internally by IterationControl.jl

> Table 2. Wrapped controls


## Using training losses, and controlling model tuning

Some iterative models report a training loss, as a byproduct of a
`fit!` call, and these can be used in two ways:

1. To supplement an out-of-sample estimate of the loss in deciding when to stop, as in the `PQ` stopping criterion (see [Prechelt, Lutz (1998)](https://link.springer.com/chapter/10.1007%2F3-540-49430-8_3))); or

2. As a (generally less reliable) substitute for an out-of-sample loss, when wishing to train excusivley on all supplied data.

To have `IteratedModel` bind all data to the training machine and use
training losses in place of an out-of-sample loss, specify
`resampling=nothing`. To check if `MyFavoriteIterativeModel` reports
training losses, load the model code and inspect
`supports_training_losses(MyFavoriteIterativeModel)` (or do
`info("MyFavoriteIterativeModel")`)


### Controlling model tuning

An example of scenario 2 occurs when controlling hyper-parameter
optimization (model tuning). Recall that MLJ's [`TunedModel`](@ref
MLJTuning.TunedModel) wrapper is implemented as an iterative
model. Moreover, this wrapper reports, as a training loss, the lowest
value of the optimization objective function so far (typically the
lowest value of an out-of-sample loss, or -1 times an out-of-sample
score). One may want to simply end the hyper-parameter search when
this value meets the [`NumberSinceBest`](@ref
EarlyStopping.NumberSinceBest) stopping criterion discussed below,
say, rather than introduce an extra layer of resampling to first
"learn" the optimal value of the iteration parameter.

In the following example we conduct a [`RandomSearch`](@ref
MLJTuning.RandomSearch) for the optimal value of the regularization
parameter `lambda` in a `RidgeRegressor` using 6-fold
cross-validation. By wrapping our "self-tuning" version of the
regressor as an [`IteratedModel`](@ref MLJIteration.IteratedModel),
with `resampling=nothing` and `NumberSinceBest(20)` in the controls,
we terminate the search when the number of `lambda` values tested
since the previous best cross-validation loss reaches 20.

```@example gree
using MLJ

X, y = @load_boston;
RidgeRegressor = @load RidgeRegressor pkg=MLJLinearModels verbosity=0
model = RidgeRegressor()
r = range(model, :lambda, lower=-1, upper=2, scale=x->10^x)
self_tuning_model = TunedModel(model=model,
                               tuning=RandomSearch(rng=123),
                               resampling=CV(nfolds=6),
                               range=r,
                               measure=mae);
iterated_model = IteratedModel(model=self_tuning_model,
                               resampling=nothing,
                               control=[Step(1), NumberSinceBest(20), NumberLimit(1000)])
mach = machine(iterated_model, X, y)
nothing # hide
```

```@repl gree
fit!(mach)
```

```@repl gree
report(mach).model_report.best_model
```

We can use `mach` here to directly obtain predictions using the
optimal model (trained on all data), as in

```@repl gree
predict(mach, selectrows(X, 1:4))
```


## Custom controls

Under the hood, control in MLJIteration is implemented using
[IterationControl.jl](https://github.com/ablaom/IterationControl.jl). Rather
than iterating a training machine directly, we iterate a wrapped
version of this object, which includes other information that controls
may want to access, such the MLJ evaluation object. This information
is summarized under [The training machine wrapper](@ref) below.

Controls must implement two `update!` methods, one for initializing
the control's *state* on the first application of the control (this
state being external to the control `struct`) and one for all
subsequent control applications, which generally updates state
also. There are two optional methods: `done`, for specifying
conditions triggering a stop, and `takedown` for specifying actions to
perform at the end of controlled training.

We summarize the training algorithm, as it relates to controls, after
giving a simple example.


### Example 1 - Non-uniform iteration steps

Below we define a control, `IterateFromList(list)`, to train, on each
application of the control, until the iteration count reaches the next
value in a user-specified `list`, triggering a stop when the `list` is
exhausted. For example, to train on iteration counts on a log scale,
one might use `IterateFromList([round(Int, 10^x) for x in range(1, 2,
length=10)]`.

In the code, `wrapper` is an object that wraps the training machine
(see above). The variable `n` is a counter for control cycles (unused
in this example).

```julia

import IterationControl # or MLJ.IterationControl

struct IterateFromList
    list::Vector{<:Int} # list of iteration parameter values
    IterateFromList(v) = new(unique(sort(v)))
end

function IterationControl.update!(control::IterateFromList, wrapper, verbosity, n)
    Δi = control.list[1]
    verbosity > 1 && @info "Training $Δi more iterations. "
    MLJIteration.train!(wrapper, Δi) # trains the training machine
    return (index = 2, )
end

function IterationControl.update!(control::IterateFromList, wrapper, verbosity, n, state)
    index = state.positioin_in_list
    Δi = control.list[i] - wrapper.n_iterations
    verbosity > 1 && @info "Training $Δi more iterations. "
    MLJIteration.train!(wrapper, Δi)
    return (index = index + 1, )
end
```

The first `update` method will be called the first time the control is
applied, returning an initialized `state = (index = 2,)`, which is
passed to the second `update` method, which is called on subsequent
control applications (and which returns the updated `state`).

A `done` method articulates the criterion for stopping:

```julia
IterationControl.done(control::IterateFromList, state) =
    state.index > length(control.list)
```

For the sake of illustration, we'll implement a `takedown` method; its
return value is included in the `IteratedModel` report:

```julia
IterationControl.takedown(control::IterateFromList, verbosity, state)
    verbosity > 1 && = @info "Stepped through these values of the "*
                              "iteration parameter: $(control.list)"
    return (iteration_values=control.list, )
end
```


### The training machine wrapper

A training machine `wrapper` has these properties:

- `wrapper.machine` - the training machine, type `Machine`

- `wrapper.model`   - the mutable atomic model, coinciding with `wrapper.machine.model`

- `wrapper.n_cycles` - the number `IterationControl.train!(wrapper, _)` calls
  so far; generally the current control cycle count

- `wrapper.n_iterations` - the total number of iterations applied to the model so far

- `wrapper.Δiterations` - the number of iterations applied in the last
  `IterationControl.train!(wrapper, _)` call

- `wrapper.loss` - the out-of-sample loss (based on the first measure in `measures`)

- `wrapper.training_losses` - the last batch of training losses (if
  reported by `model`), an abstract vector of length
  `wrapper.Δiteration`.

- `wrapper.evaluation` - the complete MLJ performance evaluation
  object, which has the following properties: `measure`,
  `measurement`, `per_fold`, `per_observation`,
  `fitted_params_per_fold`, `report_per_fold` (here there is only one
  fold). For further details, see [Evaluating Model Performance](@ref).


### The training algorithm

Here now is a simplified description of the training of an
`IteratedModel`. First, the atomic `model` is bound in a machine - the
*training machine* above - to a subset of the supplied data, and then
wrapped in an object called `wrapper` below. To train the training
machine machine for `i` more iterations, and update the other data in
the wrapper, requires the call `MLJIteration.train!(wrapper, i)`. Only
controls can make this call (e.g., `Step(...)`, or
`IterateFromList(...)` above). If we assume for simplicity there is
only a single control, called `control`, then training proceeds as
follows:

```julia
n = 1 # initialize control cycle counter
state = update!(control, wrapper, verbosity, n)
finished = done(control, state)

# subsequent training events:
while !finished
    n += 1
    state = update!(control, wrapper, verbosity, n, state)
    finished = done(control, state)
end

# finalization:
return takedown(control, verbosity, state)
```


### Example 2 - Cyclic learning rates

The control below implements a triangular cyclic learning rate policy,
close to that described in [L. N. Smith
(2019)](https://ieeexplore.ieee.org/document/7926641): "Cyclical
Learning Rates for Training Neural Networks," 2017 IEEE Winter
Conference on Applications of Computer Vision (WACV), Santa Rosa, CA,
USA, pp. 464-472. [In that paper learning rates are mutated (slowly)
*during* each training iteration (epoch) while here mutations can occur
once per iteration of the model, at most.]

For the sake of illustration, we suppose the iterative model, `model`,
specified in the `IteratedModel` constructor, has a field called
`:learning_parameter`, and that mutating this parameter does not
trigger cold-restarts.

```julia
struct CycleLearningRate{F<:AbstractFloat}
    stepsize::Int
    lower::F
    upper::F
end

# return one cycle of learning rate values:
function one_cycle(c::CycleLearningRate)
    rise = range(c.lower, c.upper, length=c.stepsize + 1)
    fall = reverse(rise)
    return vcat(rise[1:end - 1], fall[1:end - 1])
end

function IterationControl.update!(control::CycleLearningRate,
                                  wrapper,
                                  verbosity,
                                  n,
                                  state = (learning_rates=nothing, ))
    rates = n == 0 ? one_cycle(control) : state.learning_rates
    index = mod(n, length(rates)) + 1
    r = rates[index]
    verbosity > 1 && @info "learning rate: $r"
    wrapper.model.iteration_control = r
    return (learning_rates = rates,)
end
```


## API Reference

```@docs
MLJIteration.IteratedModel
```

### Controls

```@docs
IterationControl.Step
EarlyStopping.TimeLimit
EarlyStopping.NumberLimit
EarlyStopping.NumberSinceBest
EarlyStopping.InvalidValue
EarlyStopping.Threshold
EarlyStopping.GL
EarlyStopping.PQ
EarlyStopping.Patience
IterationControl.Info
IterationControl.Warn
IterationControl.Error
IterationControl.Callback
IterationControl.WithNumberDo
MLJIteration.WithIterationsDo
IterationControl.WithLossDo
IterationControl.WithTrainingLossesDo
MLJIteration.WithEvaluationDo
MLJIteration.WithFittedParamsDo
MLJIteration.WithReportDo
MLJIteration.WithModelDo
MLJIteration.WithMachineDo
MLJSerialization.Save
```

### Control wrappers

```@docs
IterationControl.skip
IterationControl.louder
IterationControl.with_state_do
IterationControl.composite
```
# Linear Pipelines

In MLJ a *pipeline* is a composite model in which models are chained
together in a linear (non-branching) chain. For other arrangements,
including custom architectures via learning networks, see [Composing
Models](@ref).

For purposes of illustration, consider a supervised learning problem
with the following toy data:

```@setup 7
using MLJ
MLJ.color_off()
```

```@example 7
using MLJ
X = (age    = [23, 45, 34, 25, 67],
     gender = categorical(['m', 'm', 'f', 'm', 'f']));
y = [67.0, 81.5, 55.6, 90.0, 61.1]
     nothing # hide
```

We would like to train using a K-nearest neighbor model, but the
model type `KNNRegressor` assumes the features are all
`Continuous`. This can be fixed by first:

- coercing the `:age` feature to have `Continuous` type by replacing
  `X` with `coerce(X, :age=>Continuous)`
- standardizing continuous features and one-hot encoding the
  `Multiclass` features using the `ContinuousEncoder` model
  
However, we can avoid separately applying these preprocessing steps
(two of which require `fit!` steps) by combining them with the
supervised `KKNRegressor` model in a new *pipeline* model, using
Julia's `|>` syntax:

```@example 7
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels
pipe = (X -> coerce(X, :age=>Continuous)) |> ContinuousEncoder() |> KNNRegressor(K=2)
```

We see above that `pipe` is a model whose hyperparameters are
themselves other models or a function. (The names of these
hyper-parameters are automatically generated. To specify your own
names, use the explicit [`Pipeline`](@ref) constructor instead.)

The `|>` syntax can also be used to extend an existing pipeline or
concatenate two existing pipelines. So, we could instead have defined:

```julia
pipe_transformer = (X -> coerce(X, :age=>Continuous)) |> ContinuousEncoder()
pipe = pipe_transformer |> KNNRegressor(K=2)
```

A pipeline is just a model like any other. For example, we can
evaluate it's performance on the data above:

```@example 7
evaluate(pipe, X, y, resampling=CV(nfolds=3), measure=mae)
```

To include target transformations in a pipeline, wrap the supervised
component using [`TransformedTargetModel`](@ref).


```@docs
Pipeline
```
# Third Party Packages

A list of third party packages with integration with MLJ.

Last updated December 2020.

Pull requests to update this list are very welcome. Otherwise, you may
post an issue requesting this
[here](https://github.com/alan-turing-institute/MLJ.jl/issues).

## Packages providing models in the MLJ model registry

See [List of Supported Models](@ref model_list)


## Providing unregistered models:

- [SossMLJ.jl](https://github.com/cscherrer/SossMLJ.jl)
- [TimeSeriesClassification](https://github.com/alan-turing-institute/TimeSeriesClassification.jl)

## Packages providing other kinds of functionality:

- [MLJParticleSwarmOptimization.jl](https://github.com/JuliaAI/MLJParticleSwarmOptimization.jl) (hyper-parameter optimization strategy)
- [TreeParzen.jl](https://github.com/IQVIA-ML/TreeParzen.jl) (hyper-parameter optimization strategy)
- [Shapley.jl](https://gitlab.com/ExpandingMan/Shapley.jl) (feature ranking / interpretation)
- [ShapML.jl](https://github.com/nredell/ShapML.jl) (feature ranking / interpretation)
- [Fairness.jl](https://github.com/ashryaagr/Fairness.jl) (FAIRness metrics)
- [OutlierDetection.jl](https://github.com/OutlierDetectionJL/OutlierDetection.jl/blob/master/src/mlj_wrappers.jl) (provides the `ProbabilisticDetector` wrapper and other outlier detection meta-functionality)


# Simple User Defined Models


To quickly implement a new supervised model in MLJ, it suffices to:

- Define a `mutable struct` to store hyperparameters. This is either a subtype
  of `Probabilistic` or `Deterministic`, depending on
  whether probabilistic or ordinary point predictions are
  intended. This `struct` is the *model*.

- Define a `fit` method, dispatched on the model, returning
  learned parameters, also known as the *fitresult*.

- Define a `predict` method, dispatched on the model, and the
  fitresult, to return predictions on new patterns.

In the examples below, the training input `X` of `fit`, and the new
input `Xnew` passed to `predict`, are tables. Each training target `y`
is a `AbstractVector`.

The predictions returned by `predict` have the same form as `y` for
deterministic models, but are `Vector`s of distributions for
probabilistic models.

Advanced model functionality not addressed here includes: (i) optional
`update` method to avoid redundant calculations when calling `fit!` on
machines a second time; (ii) reporting extra training-related
statistics; (iii) exposing model-specific functionality; (iv) checking
the scientific type of data passed to your model in `machine`
construction; and (iv) checking validity of hyperparameter values. All
this is described in [Adding Models for General
Use](adding_models_for_general_use.md).

For an unsupervised model, implement `transform` and, optionally,
`inverse_transform` using the same signature at `predict` below.

## A simple deterministic regressor

Here's a quick-and-dirty implementation of a ridge regressor with no intercept:

```julia
import MLJBase
using LinearAlgebra

mutable struct MyRegressor <: MLJBase.Deterministic
    lambda::Float64
end
MyRegressor(; lambda=0.1) = MyRegressor(lambda)

# fit returns coefficients minimizing a penalized rms loss function:
function MLJBase.fit(model::MyRegressor, verbosity, X, y)
    x = MLJBase.matrix(X)                     # convert table to matrix
    fitresult = (x'x + model.lambda*I)\(x'y)  # the coefficients
    cache=nothing
    report=nothing
    return fitresult, cache, report
end

# predict uses coefficients to make new prediction:
MLJBase.predict(::MyRegressor, fitresult, Xnew) = MLJBase.matrix(Xnew) * fitresult
```

``` @setup regressor_example
using MLJ
import MLJBase
using LinearAlgebra
MLJBase.color_off()
mutable struct MyRegressor <: MLJBase.Deterministic
    lambda::Float64
end
MyRegressor(; lambda=0.1) = MyRegressor(lambda)
function MLJBase.fit(model::MyRegressor, verbosity, X, y)
    x = MLJBase.matrix(X)
    fitresult = (x'x + model.lambda*I)\(x'y)
    cache=nothing
    report=nothing
    return fitresult, cache, report
end
MLJBase.predict(::MyRegressor, fitresult, Xnew) = MLJBase.matrix(Xnew) * fitresult
```

After loading this code, all MLJ's basic meta-algorithms can be applied to `MyRegressor`:

```@repl regressor_example
X, y = @load_boston;
model = MyRegressor(lambda=1.0)
regressor = machine(model, X, y)
evaluate!(regressor, resampling=CV(), measure=rms, verbosity=0)

```

## A simple probabilistic classifier

The following probabilistic model simply fits a probability
distribution to the `MultiClass` training target (i.e., ignores `X`)
and returns this pdf for any new pattern:

```julia
import MLJBase
import Distributions

struct MyClassifier <: MLJBase.Probabilistic
end

# `fit` ignores the inputs X and returns the training target y
# probability distribution:
function MLJBase.fit(model::MyClassifier, verbosity, X, y)
    fitresult = Distributions.fit(MLJBase.UnivariateFinite, y)
    cache = nothing
    report = nothing
    return fitresult, cache, report
end

# `predict` returns the passed fitresult (pdf) for all new patterns:
MLJBase.predict(model::MyClassifier, fitresult, Xnew) =
    [fitresult for r in 1:nrows(Xnew)]
```

```julia
julia> X, y = @load_iris
julia> mach = fit!(machine(MyClassifier(), X, y))
julia> predict(mach, selectrows(X, 1:2))
2-element Array{UnivariateFinite{String,UInt32,Float64},1}:
 UnivariateFinite(setosa=>0.333, versicolor=>0.333, virginica=>0.333)
 UnivariateFinite(setosa=>0.333, versicolor=>0.333, virginica=>0.333)
```
