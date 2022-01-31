---
title: 'Grama: A Grammar of Model Analysis'
tags:
  - Python
  - Modeling
  - Uncertainty quantification
  - Functional programming
  - Pedagogy
  - Communication
  - Reproducibility
authors:
  - name: Zachary del Rosario
    orcid: 0000-0003-4676-1692
    affiliation: 1
affiliations:
 - name: Visiting Professor, Olin College of Engineering
   index: 1
date:
bibliography: paper.bib
---

# Summary

`Grama` is a Python package implementing a *functional grammar of model
analysis* emphasizing the quantification of uncertainties. In `Grama` a *model*
contains both a function mapping inputs to outputs as well as a distribution
characterizing uncertainties on those inputs. This conceptual object unifies the
engineer/scientist's definition of a model with that of a statistician. `Grama`
provides an *implementation* of this model concept, as well as *verbs* to carry
out model-building and model-analysis.

## Statement of Need

Uncertainty Quantification (UQ) is the science of analyzing uncertainty in
scientific problems and using those results to inform decisions. UQ has
important applications to building safety-critical engineering systems, and to
making high-consequence choices based on scientific models. However, UQ is
generally not taught at the undergraduate level: Many engineers leave their
undergraduate training with a purely deterministic view of their discipline,
which can lead to probabilistic design errors that negatively impact safety
[@delRosario2020design]. To that end, I have developed a grammar of model
analysis---`Grama`---to facilitate rapid model analysis, communication of
results, and the teaching of concepts, all with quantified uncertainties.
Intended users of `Grama` are scientists and engineers at the undergraduate
level and upward, seeking to analyze computationally-lightweight models.

## Differentiating Attributes

Packages similar to `Grama` exist, most notably Sandia National Lab's `Dakota`
[@adams2017sandia] and `UQLab` [@marelli2014uqlab] out of ETH Zurich. While both
of these packages are mature and highly featured, `Grama` has several
differentiating attributes. First, `Grama` emphasizes an explicit but flexible
*model object*: this object enables sharp decomposition of a UQ problem into a
model-building stage and a model-analysis stage. This logical decomposition
enables simplified syntax and a significant reduction in boilerplate code.
Second, `Grama` implements a functional programming syntax to emphasize
operations against the model object, improving readability of code. Finally,
`Grama` is designed from the ground-up as a pedagogical and communication tool.
For learnability: Its *verb-prefix* syntax is meant to remind the user how
functions are used based solely on their name, and the package is shipped with
fill-in-the-blank Jupyter notebooks [@kluyver2016jupyter] to take advantage of
the pedagogical benefits of active learning [@freeman2014active]. For
communication: The model object and functional syntax abstract away numerical
details for presentation in a notebook, while preserving tracability and
reproducibility of results through the inspection of source code.

## Inspiration and Dependencies

``Grama`` relies heavily on the SciKit package ecosystem for its numerical
backbone
[@scipy2020;@numpy2011;@matplotlib2007;@mckinney-proc-scipy-2010;@pedregosa2011scikit].
The functional design is heavily inspired by the `Tidyverse`
[@wickham2019tidyverse], while its implementation is built upon `dfply`
[@kiefer2019dfply]. Additional functionality for materials data via an optional
dependency on Matminer [@ward2018matminer].

# Acknowledgements

I acknowledge contributions from Richard W. Fenrich on the laminate plate model.

# References
# Contributing to `py_grama`

Thanks for your interest in contributing! Please follow the guidelines below, and contact one of the Maintainers (listed below) if you have any questions.

## Reporting an issue

If you find a issue with the software, please file a [new issue](https://github.com/zdelrosario/py_grama/issues). Please include a [reproducible example](https://stackoverflow.com/help/minimal-reproducible-example) in your issue. Note that if you have a question, please don't file an issue---the issue tracker is meant to document issues in the design and implementation of `py_grama`, not to answer questions.

## Contributing to the software

We welcome contributions to `py_grama`! To contribute, please determine what sort of contribution you plan to make. For more detailed information on forking and branching in the context of contributing, please see [this guide](https://opensource.com/article/19/7/create-pull-request-github).

### Bug-fix

If you find an issue with the software or want to fix an existing issue, please follow these steps:

1. If one does not yet exist, please report an Issue following the instructions in **Reporting an issue** above.
2. Fork `py_grama` and clone a local copy of the repository for your work.
3. Create a branch for your fix with the name `fix_name`, where `name` should be sensibly related to your fix.
4. If one does not already exist, create a unittest in the [tests](https://github.com/zdelrosario/py_grama/tree/master/tests) that captures the bug.
5. Implement your fix.
6. Verify the fix against the test suite; the `Makefile` in `py_grama` automates testing with the spell `make test`.
7. Create a pull request against `py_grama`; one of the Maintainers will review your contribution.

### Feature addition

If you wish to add a new feature to `py_grama`, please follow these steps:

1. Fork `py_grama` and clone a local copy of the repository for your work.
2. Create a branch for your fix with the name `dev_name`, where `name` should be sensibly related to your fix.
3. Create a an appropriate set of tests in the [tests](https://github.com/zdelrosario/py_grama/tree/master/tests) that verify the functionality of your feature.
4. Implement your feature.
5. Verify the feature against the test suite; the `Makefile` in `py_grama` automates testing with the spell `make test`.
6. Create a pull request against `py_grama`; one of the Maintainers will review your contribution.

### Design change

If you have a suggestion for a significant change to the design of `py_grama`, please reach out directly to one of the Maintainers (listed below).

## List of Maintainers

Updated 2020-06-05

- Zachary del Rosario (zdelrosario(at)outlook(doot)com)
# py_grama
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02462/status.svg)](https://doi.org/10.21105/joss.02462) [![PyPI version](https://badge.fury.io/py/py-grama.svg)](https://badge.fury.io/py/py-grama) [![Documentation Status](https://readthedocs.org/projects/py_grama/badge/?version=latest)](https://py_grama.readthedocs.io/en/latest/?badge=latest) ![Python package test](https://github.com/zdelrosario/py_grama/workflows/Python%20package%20test/badge.svg) [![codecov](https://codecov.io/gh/zdelrosario/py_grama/branch/master/graph/badge.svg)](https://codecov.io/gh/zdelrosario/py_grama) [![CodeFactor](https://www.codefactor.io/repository/github/zdelrosario/py_grama/badge/master)](https://www.codefactor.io/repository/github/zdelrosario/py_grama/overview/master) 

Implementation of a *grammar of model analysis* (*grama*). See the [documentation](https://py-grama.readthedocs.io/en/latest/) for more info.

# Overview

Grama is a *grammar of model analysis*---a Python package that supports building and analyzing models with quantified uncertainties. This language is heavily inspired by the [Tidyverse](https://www.tidyverse.org/). Grama provides convenient syntax for building a model (with functions and distributions), generating data, and visualizing results. The purpose of this language is to support scientists and engineers learning to handle uncertainty, and to improve documentation + reproducibility of results.

Uncertainty Quantification (UQ) is the science of analyzing uncertainty in scientific problems and using those results to inform decisions. UQ has important applications to building safety-critical engineering systems, and to making high-consequence choices based on scientific models. However, UQ is generally not taught at the undergraduate level: Many engineers leave their undergraduate training with a purely deterministic view of their discipline, which can lead to probabilistic design errors that [negatively impact safety](https://arc.aiaa.org/doi/abs/10.2514/6.2020-0414). To that end, Grama is designed to facilitate rapid model analysis, communication of results, and the teaching of concepts, all with quantified uncertainties. Intended users of `Grama` are scientists and engineers at the undergraduate level and upward, seeking to analyze computationally-lightweight models.

# Installation
Quick install:

```bash
$ pip install py-grama
```

For a manual install clone this repo, change directories and run the following to install dependencies. (Note: I recommend [Anaconda](https://www.anaconda.com/distribution/) as a Python distribution; it takes care of most of the dependencies.)

```bash
$ git clone git@github.com:zdelrosario/py_grama.git
$ cd py_grama/
$ pip install -r requirements.txt
$ pip install .
```

Run the following to check your install:

```bash
$ python
> import grama
```

# Quick Tour
`py_grama` has tools for both *building* and *analyzing* models. For a quick look at functionality, see the following notebooks:

- [video demo](https://youtu.be/jhyB-jQ7EC8)
- [model building demo](https://github.com/zdelrosario/py_grama/blob/master/examples/demo/builder_demo.ipynb)
- [model analysis demo](https://github.com/zdelrosario/py_grama/blob/master/examples/demo/analysis_demo.ipynb)

# Tutorials
The [tutorials](https://github.com/zdelrosario/py_grama/tree/master/tutorials) page has educational materials for learning to work with `py_grama`.

# Support and Contributing
If you are seeking support or want to contribute, please see [Contributing](https://github.com/zdelrosario/py_grama/blob/master/contributing.md).

# Cite As

If you find Grama useful in your work, we'd appreciate that you cite it as:

> del Rosario, Z., (2020). Grama: A Grammar of Model Analysis. Journal of Open Source Software, 5(51), 2462, https://doi.org/10.21105/joss.02462

Bibtex code:

```
@article{del Rosario2020,
  doi = {10.21105/joss.02462},
  url = {https://doi.org/10.21105/joss.02462},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {51},
  pages = {2462},
  author = {Zachary del Rosario},
  title = {Grama: A Grammar of Model Analysis},
  journal = {Journal of Open Source Software}
}
```
Copyright (C) 2019  Kiefer Katovich (kieferk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
# Tutorials
This directory contains self-paced Jupyter notebooks for learning `py_grama`. The `assignment` notebooks are intended to be worked through, and the `solution` notebooks are given for reference. See the [first notebook](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t01_introduction_assignment.ipynb) for advice on installation.

Curriculum:
- [t01_introduction_assignment](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t01_introduction_assignment.ipynb): Installation; intro to `py_grama` syntax and concepts
- [t02_explore_assignment](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t02_explore_assignment.ipynb): Intro to *exploratory model analysis* with `py_grama`
- [t03_building_assignment](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t03_building_assignment.ipynb): Building a fully-defined model with `py_grama`; using data to inform a model
- [t04_dag_assignment](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t04_dag_assignment.ipynb): Building hierarchical functions and debugging them with directed acyclic graphs
- [t05_data_summary_assignment](https://github.com/zdelrosario/py_grama/blob/master/tutorials/t05_data_summary_assignment.ipynb): Analyzing and aggregating data; estimating failure probabilities
# Language Details

---

The following is a fairly extensive introduction to the *grama* language. This
is (admittedly) not a great place to start learning how to *use* `py_grama`, but
is instead provided as a reference.

*grama* is a conceptual language, and `py_grama` is a code implementation of
that concept. This page is a description of both the concept and its
implementation.

## Running example

We'll use a running example throughout this page; the built-in `py_grama` Cantilever Beam model.

```python
import grama as gr
from grama.models import make_cantilever_beam

md_beam = make_cantilever_beam()
md_beam.printpretty()
```

```bash
model: Cantilever Beam

  inputs:
    var_det:
      w: [2, 4]
      t: [2, 4]
    var_rand:
      H: (+1) norm, {'loc': 500.0, 'scale': 100.0}
      V: (+1) norm, {'loc': 1000.0, 'scale': 100.0}
      E: (+0) norm, {'loc': 29000000.0, 'scale': 1450000.0}
      Y: (-1) norm, {'loc': 40000.0, 'scale': 2000.0}
  functions:
    cross-sectional area: ['w', 't'] -> ['c_area']
    limit state: stress: ['w', 't', 'H', 'V', 'E', 'Y'] -> ['g_stress']
    limit state: displacement: ['w', 't', 'H', 'V', 'E', 'Y'] -> ['g_disp']
```

## Objects

*grama* focuses on two categories of objects:

- **data** (`df`): observations on various quantities, implemented by the Python package [Pandas](https://pandas.pydata.org/)
- **models** (`md`): functions and a complete description of their inputs, implemented by [py_grama](https://github.com/zdelrosario/py_grama)

For readability, we suggest using prefixes `df_` and `md_` when naming DataFrames and models.

### Data

Data are observations on some quantities. Data often come from the real world,
but data are also used to inform models, and models can be used to generate new
data. `py_grama` uses the Pandas `DataFrame` implementation to represent data.
Since data operations are already well-handled by Pandas, `py_grama` uses the
existing Pandas infrastructure and focuses on providing tools to handle models
and their interface with data.

### Models

Models in *grama* have both functions and inputs. `py_grama` implements models
in the `Model` class, which in turn have three primary objects:

- `Domain`: Defines bounds on variables
- `Density`: Defines a joint density for variables
- List of`Function`s: Maps variables to outputs

#### Domain

The domain of the model defines bounds for all the variables. If a variable is
not included in the domain object, it is assumed to be unbounded. The model
above has bounds on `t, w`, both of which are `[2, 4]`.

#### Density

The density of the model defines a joint density for the *random* variables. If
a variable is included in the density it is random, otherwise it is
*deterministic*. The model above has a joint density on `H, V, E, Y`.
The model summary gives details on each marginal distribution.

#### Functions

A function has a set of variables, which map to a set of outputs; for instance,
the `cross-sectional area` function above maps `['w', 't'] -> ['c_area']`. The
other functions take more variables, all of which map to their respective
outputs.

#### Inputs

The full set of model inputs are organized into:

|            | Deterministic                         | Random     |
| ---------- | ------------------------------------- | ---------- |
| Variables  | `md.var_det`                          | `md.var_rand` |
| Parameters | `md.density.marginals[var].d_param`   | (Future*)  |

- **Variables** are inputs to the model's functions
  + **Deterministic** variables are chosen by the user; the model above has `w, t`
  + **Random** variables are not controlled; the model above has `H, V, E, Y`
- **Parameters** are inputs to the model's density
  + **Deterministic** parameters are currently implemented; these are listed under `var_rand` with their associated random variable
  + **Random** parameters* are not yet implemented

The full set of *variables* is determined by the domain, density, and functions.
Formally, the full set of variables is given (in pseudocode) by `domain.var +
[f.var for f in functions]`. The set of random variables is then given by
`domain.var + [f.var for f in functions] - density.marginals.keys()`, while the
deterministic variables are the remainder `var_det = var_full - var_rand`.

## Verbs

Verbs are used to take action on different *grama* objects. We use verbs to
generate data from models, build new models from data, and ultimately make sense
of the two.

The following table summarizes the categories of `py_grama` verbs. Verbs take
either data (`df`) or a model (`md`), and may return either object type. The
prefix of a verb immediately tells one both the input and output types. The
short prefix is used to denote the *pipe-enabled version* of a verb.

| Verb Type | Prefix (Short)  | In   | Out   |
| --------- | --------------- | ---- | ----- |
| Evaluate  | `eval_` (`ev_`) | `md` | `df`  |
| Fit       | `fit_`  (`ft_`) | `df` | `md`  |
| Transform | `tran_` (`tf_`) | `df` | `df`  |
| Compose   | `comp_` (`cp_`) | `md` | `md`  |
| Plot      | `plot_` (`pt_`) | `df` | (Plot) |

Since `py_grama` is focused on models, the majority of functions lie in the
Evaluate, Fit, and Compose categories, with only a few original Transform
utilities provided. Some shortcut plotting utilities are also provided
for covenience.

`py_grama` verbs are used to both *build* and *analyze* models.

### Model Building

The *recommended* way to build `py_grama` models is with *composition calls*.
Calling `Model()` creates an "empty" model, to which one can add.

```python
md = gr.Model()
md.printpretty()
```

```bash
model: None

  inputs:
    var_det:
    var_rand:
  functions:
```

We can then use Compose functions to build up a complete model step-by step. We
recommend starting with the functions, as those highlight the required
variables.

```python
md = gr.Model("Test") >> \
     gr.cp_function(
         fun=lambda x: [x[0], x[1]],
         var=["x0", "x1"],
         out=2,
         name="Identity"
     )
md.printpretty()
```

```bash
model: Test

  inputs:
    var_det:
      x0: (unbounded)
      x1: (unbounded)
    var_rand:
  functions:
    Identity: ['x0', 'x1'] -> ['y0', 'y1']
```

Note that by default all of the variables are assumed to be deterministic. We
can override this by adding marginal distributions for one or more of the
variables.

```python
md = gr.Model("Test") >> \
     gr.cp_function(
         fun=lambda x: [x[0], x[1]],
         var=["x0", "x1"],
         out=2,
         name="Identity"
     ) >> \
     gr.cp_marginals(
         x1=dict(dist="norm", loc=0, scale=1)
     )
md.printpretty()
```

```bash
model: Test

  inputs:
    var_det:
      x0: (unbounded)
    var_rand:
      x1: (+0) norm, {'loc': 0, 'scale': 1}
  functions:
    Identity: ['x0', 'x1'] -> ['y0', 'y1']
```

The marginals are implemented in terms of the Scipy [continuous
distributions](https://docs.scipy.org/doc/scipy/reference/stats.html); see the
variable `gr.valid_dist.keys()` for a list of implemented marginals. When
calling `gr.comp_marginals()`, we provide the target variable name as a keyword
argument, and the marginal information via dictionary. The marginal shape is
specified with the "dist" keyword; all distributions require the `loc, scale`
parameters, but some require additional keywords. See `gr.param_dist` for a
dictionary mapping between distributions and parameters.

Once we have constructed our model, we can analyze it with a number of tools.

### Model Analysis

One question in model analysis is to what degree the random variables affect the
outputs. A way to quantify this is with *Sobol' indices* (Sobol', 1999). We
can estimate Sobol' indices in `py_grama` with the following code.

```python
df_sobol = \
    md_beam >> \
    gr.ev_hybrid(n=1e3, df_det="nom", seed=101) >> \
    gr.tf_sobol()
print(df_sobol)
```

```bash
eval_hybrid() is rounding n...
     w    t  c_area     g_stress  g_disp  ind
0  NaN  NaN     NaN        -0.03    0.28  S_E
0  NaN  NaN     NaN         0.35    0.21  S_H
0  NaN  NaN     NaN         0.33    0.64  S_V
0  NaN  NaN     NaN         0.31    0.02  S_Y
0  0.0  0.0     0.0   -345263.39    0.01  T_E
0  0.0  0.0     0.0   4867712.30    0.01  T_H
0  0.0  0.0     0.0   4577175.85    0.02  T_V
0  0.0  0.0     0.0   4224965.01    0.00  T_Y
0  0.0  0.0     0.0  13758547.04    0.04  var
```

The normalized Sobol' indices are reported with `S_[var]` labels; they indicate
that `g_stress` is affected roughly equally by the inputs `H,V,Y`, while
`g_disp` is affected about twice as much by `V` as by `E` or `H`. Note that the
Sobol' indices are *only* defined for the random variables---since Sobol'
indices are defined in terms of fractional variances, they are only formally
valid for quantifying contributions from sources of randomness.

Under the hood `gr.eval_hybrid()` attaches metadata to its resulting DataFrame,
which `gr.tran_sobol()` detects and uses in post-processing the data.

`py_grama` also provides tools for constructing visual summaries of models. We
can construct a *sinew plot* with a couple lines of code. First we inspect the
design:

```python
md_beam >> \
    gr.ev_sinews(n_density=50, n_sweeps=10, df_det="nom", skip=True) >> \
    gr.pt_auto()
```

![beam sinew results](../images/ex_beam_sinews_doe.png)

The "sinews" are sweeps across random variable space which start at random
locations, and continue parallel to the variable axes. Evaluating these
samples allows us to construct a sinew plot:

```python
md_beam >> \
    gr.ev_sinews(n_density=50, n_sweeps=10, df_det="nom", skip=False) >> \
    gr.pt_auto()
```

![beam sinew results](../images/ex_beam_sinews_res.png)

Here we can see that inputs `H,E` tend to saturate in their effects on `g_disp`,
while `V` is linear over its domain. This may explain the difference in
contributed variance seen above via Sobol' indices.

By providing tools to quickly perform different analyses, one can quickly get a
sense of model behavior using `py_grama`.

## Layers and Defaults

`py_grama` is built around *layers and defaults*. As much as is possible
`py_grama` is designed to provide sensible defaults "out-of-the-box". We saw the
concept of *layers* above in the model building example. The following example
shows defaults in action.

### Example Defaults: `gr.eval_monte_carlo()`

Attempting to provide no arguments to `gr.eval_monte_carlo()` yields the
following error:

```python
df_res = md_beam >> gr.ev_monte_carlo()
```

```bash
...
ValueError: df_det must be DataFrame or 'nom'
```

One can sample over the random variables given their joint density, but this
tells us nothing about how to treat the deterministic variables. The error
message above tells us that we have to define the deterministic variable levels
through `df_det`. To perform simple studies, we can explicitly limit attention
to the nominal conditions.

```python
df_res = md_beam >> gr.ev_monte_carlo(df_det="nom")
print(df_res.describe())
```

```bash
                H           V             E  ...  c_area    g_stress    g_disp
count    1.000000    1.000000  1.000000e+00  ...     1.0     1.00000  1.000000
mean   505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044
std           NaN         NaN           NaN  ...     NaN         NaN       NaN
min    505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044
25%    505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044
50%    505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044
75%    505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044
max    505.743982  906.946892  2.824652e+07  ...     9.0  9758.06826  0.438044

[8 rows x 9 columns]
```

By default `gr.eval_monte_carlo()` will draw a single sample; this leads to the
`NaN` standard deviation (`std`) results. We can override this default behavior
by providing the `n` keyword.

```python
df_res = md_beam >> gr.ev_monte_carlo(df_det="nom", n=1e3)
print(df_res.describe())
```

```bash
eval_monte_carlo() is rounding n...
                 H            V             E  ...  c_area      g_stress       g_disp
count  1000.000000  1000.000000  1.000000e+03  ...  1000.0   1000.000000  1000.000000
mean    499.776481  1001.908814  2.899238e+07  ...     9.0   6589.731053     0.333387
std     100.878972   102.929434  1.523033e+06  ...     0.0   3807.292688     0.199576
min     168.466610   635.220875  2.427615e+07  ...     9.0  -5830.363269    -0.339594
25%     433.466349   933.698938  2.799035e+07  ...     9.0   3921.581411     0.200657
50%     500.201964  1002.514072  2.896188e+07  ...     9.0   6633.083925     0.339702
75%     569.159912  1065.489928  2.993924e+07  ...     9.0   9161.066426     0.465231
max     798.543698  1327.618914  3.378264e+07  ...     9.0  18682.531892     1.045742

[8 rows x 9 columns]
```

Formally `1e3` is a float, which is not a valid iteration count. The routine
`gr.eval_monte_carlo()` informs us that it first rounds the given value before
proceeding. Here we can see some variation in the inputs and outputs, though
`c_area` is clearly unaffected by the randomness.

We can also provide an explicit DataFrame to the `df_det` argument. The
`gr.eval_monte_carlo()` routine will automatically take an *outer product* of
the deterministic settings with the random samples; this will lead to a
multiplication in sample size of `df_det.shape[0] * n`. Since this can get quite
large, we should reduce `n` before proceeding. We can also *delay evaluation*
first with the `skip` keyword, and inspect the design first before evaluating
it.

```python
df_det = pd.DataFrame(dict(
    w=[3] * 10,
    t=[2.5 + i/10 for i in range(10)]
))

df_design = md_beam >> gr.ev_monte_carlo(df_det=df_det, n=1e2, skip=True)
print(df_design.describe())
```

```bash
eval_monte_carlo() is rounding n...
                 H            V             E             Y       w            t
count  1000.000000  1000.000000  1.000000e+03   1000.000000  1000.0  1000.000000
mean    485.021539   989.376268  2.890355e+07  40135.172902     3.0     2.950000
std     112.546373    92.381768  1.462286e+06   2135.788506     0.0     0.287372
min     137.606474   764.764234  2.586955e+07  33320.789307     3.0     2.500000
25%     411.652958   927.448139  2.763002e+07  39104.255172     3.0     2.700000
50%     490.765903  1002.513232  2.903287e+07  40121.853496     3.0     2.950000
75%     561.317283  1046.038620  2.994359e+07  41279.345157     3.0     3.200000
max     726.351812  1239.473097  3.176696e+07  45753.627265     3.0     3.400000
```

If we are happy with the design (possibly after visual inspection), we can
pass the input DataFrame to the straight evaluation routine

```python
df_res = md_beam >> gr.ev_df(df=df_design)
print(df_res.describe())
```

```bash
                 H            V             E  ...       c_area      g_stress       g_
disp
count  1000.000000  1000.000000  1.000000e+03  ...  1000.000000   1000.000000  1000.00
0000
mean    510.457048  1003.972079  2.887534e+07  ...     8.850000   4555.508097     0.12
0805
std     105.566743    91.443354  1.486829e+06  ...     0.862116   7027.137049     0.59
9538
min     269.102628   830.234813  2.553299e+07  ...     7.500000 -18370.732185    -1.84
8058
25%     437.342719   939.490354  2.794973e+07  ...     8.100000   -385.620306    -0.29
7274
50%     515.387982   995.468715  2.868415e+07  ...     8.850000   5347.337598     0.23
2133
75%     584.899669  1054.815972  2.982724e+07  ...     9.600000   9824.615141     0.60
9096
max     775.221819  1257.095025  3.369872e+07  ...    10.200000  21238.954122     1.23
6387

[8 rows x 9 columns]
```

## Functional Programming (Pipes)

[Functional programming](https://en.wikipedia.org/wiki/Functional_programming)
touches both the practical and conceptual aspects of the language. `py_grama`
provides tools to use functional programming patterns. Short-stem versions of
`py_grama` functions are *pipe-enabled*, meaning they can be used in functional
programming form with the pipe operator `>>`. These pipe-enabled functions are
simply aliases for the base functions, as demonstrated below:

```python
df_base = gr.eval_nominal(md_beam, df_det="nom")
df_functional = md_beam >> gr.ev_nominal(df_det="nom")

df_base.equals(df_functional)
```

```bash
True
```

Functional patterns enable chaining multiple commands, as demonstrated in the
following Sobol' index analysis. In nested form using base functions, this would
be:

```python
df_sobol = gr.tran_sobol(gr.eval_hybrid(md_beam, n=1e3, df_det="nom", seed=101))
```

From the code above, it is difficult to see that we first consider `md_beam`, perform a hybrid-point evaluation, then use those data to estimate Sobol' indices. With more chained functions, this only becomes more difficult. One could make the code significantly more readable by introducing intermediate variables:

```python
df_samples = gr.eval_hybrid(md_beam, n=1e3, df_det="nom", seed=101)
df_sobol = gr.tran_sobol(df_samples)
```

Conceptually, using *pipe-enabled* functions allows one to skip assigning intermediate variables, and instead pass results along to the next function. The pipe operator `>>` inserts the results of one function as the first argument of the next function. A pipe-enabled version of the code above would be:

```python
df_sobol = \
    md_beam >> \
    gr.ev_hybrid(n=1e3, df_det="nom", seed=101) >> \
    gr.tf_sobol()
```

## References

- I.M. Sobol', "Sensitivity Estimates for Nonlinear Mathematical Models" (1999) MMCE, Vol 1.
# Random Variable Modeling

---

TODO
# Overview

---

The `py_grama` package is an implementation of a *grammar of model analysis*
(*grama*)---a language for describing and analyzing models.

## What Models?

Statisticians often use "model" to refer to random variable models. Scientists
and engineers often use "model" to refer to simplified physics resulting in
function models. In *grama* we refer to a collection of random variables and
functions *together* as a model.

## Why *grama*?

Considering both the functional mapping between variables and the uncertainties
in those variables is of critical importance to a full understanding of a given
problem. Given the "split" perspective between statisticians and engineers,
unifying the perspectives is a conceptual challenge.

While much effort in the [uncertainty
quantification](https://en.wikipedia.org/wiki/Uncertainty_quantification) (UQ)
community has been made on merging the two perspectives on the *algorithmic*
side, relatively little work has been done to merge the two perspectives
*conceptually*. The aforementioned understanding of "models"---functions plus
random variables---is a step towards conceptually unifying these two
perspectives.

## Why `py_grama`?

Furthermore, virtually *no* work has been done to make UQ techniques easily
learnable and accessible. The `py_grama` package is heavily inspired by the
[Tidyverse](https://www.tidyverse.org/), partly in terms of functional
programming patterns, but primarily in terms of its *user-first perspective*.
`py_grama` is designed to help users learn and use UQ tools to analyze models.

## Why quantify uncertainty?

Uncertainty quantification is a relatively new scientific discipline, so the
motivation for doing UQ may not be immediately obvious. The following example notebooks demonstrate UQ in a number of settings:

- [structural safety: cable design example](https://github.com/zdelrosario/py_grama/blob/master/examples/tension/tension.ipynb)---failing to account for uncertainty can lead to unsafe structures. UQ enables safer design.

## What does it look like?

For a quick demonstration of `py_grama`, see the following demo notebooks:

- The [model building demo](https://github.com/zdelrosario/py_grama/blob/master/examples/demo/builder_demo.ipynb) shows how to build a *grama* model in a scientifically-reproducible way.
- The [model analysis demo](https://github.com/zdelrosario/py_grama/blob/master/examples/demo/analysis_demo.ipynb) shows how *grama* can be used to analyze an existing model, using compact syntax to probe how both functions and randomness affect model outputs.
.. py_grama documentation master file, created by
   sphinx-quickstart on Mon Dec 23 22:03:55 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to py_grama's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   source/overview
   source/language
   source/rv_modeling
   source/modules

.. autosummary::
   :toctree: stubs

.. automodule:: grama
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
grama.models package
====================

Submodules
----------

grama.models.cantilever\_beam module
------------------------------------

.. automodule:: grama.models.cantilever_beam
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.circuit\_RLC module
--------------------------------

.. automodule:: grama.models.circuit_RLC
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.ishigami module
----------------------------

.. automodule:: grama.models.ishigami
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.linear\_normal module
----------------------------------

.. automodule:: grama.models.linear_normal
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.plane\_laminate module
-----------------------------------

.. automodule:: grama.models.plane_laminate
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.plate\_buckling module
-----------------------------------

.. automodule:: grama.models.plate_buckling
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.poly module
------------------------

.. automodule:: grama.models.poly
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.test module
------------------------

.. automodule:: grama.models.test
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.time\_cantilever module
------------------------------------

.. automodule:: grama.models.time_cantilever
   :members:
   :undoc-members:
   :show-inheritance:

grama.models.trajectory\_linear\_drag module
--------------------------------------------

.. automodule:: grama.models.trajectory_linear_drag
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.models
   :members:
   :undoc-members:
   :show-inheritance:
grama.dfply package
===================

Submodules
----------

grama.dfply.base module
-----------------------

.. automodule:: grama.dfply.base
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.group module
------------------------

.. automodule:: grama.dfply.group
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.join module
-----------------------

.. automodule:: grama.dfply.join
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.mask\_helpers module
--------------------------------

.. automodule:: grama.dfply.mask_helpers
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.reshape module
--------------------------

.. automodule:: grama.dfply.reshape
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.select module
-------------------------

.. automodule:: grama.dfply.select
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.set\_ops module
---------------------------

.. automodule:: grama.dfply.set_ops
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.string\_helpers module
----------------------------------

.. automodule:: grama.dfply.string_helpers
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.subset module
-------------------------

.. automodule:: grama.dfply.subset
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.summarize module
----------------------------

.. automodule:: grama.dfply.summarize
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.summary\_functions module
-------------------------------------

.. automodule:: grama.dfply.summary_functions
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.transform module
----------------------------

.. automodule:: grama.dfply.transform
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.vector module
-------------------------

.. automodule:: grama.dfply.vector
   :members:
   :undoc-members:
   :show-inheritance:

grama.dfply.window\_functions module
------------------------------------

.. automodule:: grama.dfply.window_functions
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.dfply
   :members:
   :undoc-members:
   :show-inheritance:
grama.fit package
=================

Submodules
----------

grama.fit.fit\_lolo module
--------------------------

.. automodule:: grama.fit.fit_lolo
   :members:
   :undoc-members:
   :show-inheritance:

grama.fit.fit\_scikitlearn module
---------------------------------

.. automodule:: grama.fit.fit_scikitlearn
   :members:
   :undoc-members:
   :show-inheritance:

grama.fit.fit\_statsmodels module
---------------------------------

.. automodule:: grama.fit.fit_statsmodels
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.fit
   :members:
   :undoc-members:
   :show-inheritance:
grama.data package
==================

Submodules
----------

grama.data.datasets module
--------------------------

.. automodule:: grama.data.datasets
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.data
   :members:
   :undoc-members:
   :show-inheritance:
grama package
=============

Subpackages
-----------

.. toctree::

   grama.data
   grama.dfply
   grama.eval
   grama.fit
   grama.models
   grama.tran

Submodules
----------

grama.comp\_building module
---------------------------

.. automodule:: grama.comp_building
   :members:
   :undoc-members:
   :show-inheritance:

grama.comp\_metamodels module
-----------------------------

.. automodule:: grama.comp_metamodels
   :members:
   :undoc-members:
   :show-inheritance:

grama.core module
-----------------

.. automodule:: grama.core
   :members:
   :undoc-members:
   :show-inheritance:

grama.eval\_defaults module
---------------------------

.. automodule:: grama.eval_defaults
   :members:
   :undoc-members:
   :show-inheritance:

grama.eval\_opt module
----------------------

.. automodule:: grama.eval_opt
   :members:
   :undoc-members:
   :show-inheritance:

grama.eval\_random module
-------------------------

.. automodule:: grama.eval_random
   :members:
   :undoc-members:
   :show-inheritance:

grama.eval\_tail module
-----------------------

.. automodule:: grama.eval_tail
   :members:
   :undoc-members:
   :show-inheritance:

grama.fit\_synonyms module
--------------------------

.. automodule:: grama.fit_synonyms
   :members:
   :undoc-members:
   :show-inheritance:

grama.marginals module
----------------------

.. automodule:: grama.marginals
   :members:
   :undoc-members:
   :show-inheritance:

grama.mutate\_helpers module
----------------------------

.. automodule:: grama.mutate_helpers
   :members:
   :undoc-members:
   :show-inheritance:

grama.plot\_auto module
-----------------------

.. automodule:: grama.plot_auto
   :members:
   :undoc-members:
   :show-inheritance:

grama.string\_helpers module
----------------------------

.. automodule:: grama.string_helpers
   :members:
   :undoc-members:
   :show-inheritance:

grama.support module
--------------------

.. automodule:: grama.support
   :members:
   :undoc-members:
   :show-inheritance:

grama.tools module
------------------

.. automodule:: grama.tools
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran\_pivot module
------------------------

.. automodule:: grama.tran_pivot
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran\_shapley module
--------------------------

.. automodule:: grama.tran_shapley
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran\_summaries module
----------------------------

.. automodule:: grama.tran_summaries
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran\_tools module
------------------------

.. automodule:: grama.tran_tools
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama
   :members:
   :undoc-members:
   :show-inheritance:
grama
=====

.. toctree::
   :maxdepth: 4

   grama
grama.eval package
==================

Submodules
----------

grama.eval.eval\_pyDOE module
-----------------------------

.. automodule:: grama.eval.eval_pyDOE
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.eval
   :members:
   :undoc-members:
   :show-inheritance:
grama.tran package
==================

Submodules
----------

grama.tran.tran\_matminer module
--------------------------------

.. automodule:: grama.tran.tran_matminer
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran.tran\_scikitlearn module
-----------------------------------

.. automodule:: grama.tran.tran_scikitlearn
   :members:
   :undoc-members:
   :show-inheritance:

grama.tran.tran\_umap module
----------------------------

.. automodule:: grama.tran.tran_umap
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: grama.tran
   :members:
   :undoc-members:
   :show-inheritance:
