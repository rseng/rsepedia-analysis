# MF2: Multi-Fidelity-Functions

| Package Info | Status | Support |
|--------------|--------|---------|
| [![PyPI version](https://badge.fury.io/py/mf2.svg)](https://badge.fury.io/py/mf2) | [![Build Status](https://app.travis-ci.com/sjvrijn/mf2.svg?branch=main)](https://app.travis-ci.com/sjvrijn/mf2) | [![Documentation Status](https://readthedocs.org/projects/mf2/badge/?version=latest)][docs-badge] |
| [![Conda](https://img.shields.io/conda/v/conda-forge/mf2)](https://anaconda.org/conda-forge/mf2) | [![Coverage Status](https://coveralls.io/repos/github/sjvrijn/mf2/badge.svg?branch=master)](https://coveralls.io/github/sjvrijn/mf2?branch=master) | [![Gitter](https://badges.gitter.im/pymf2/community.svg)][gitter-badge] |
| ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mf2) | [![Codacy Badge](https://api.codacy.com/project/badge/Grade/54144e7d406b4558a14996b06a89adf8)](https://www.codacy.com/manual/sjvrijn/mf2?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=sjvrijn/mf2&amp;utm_campaign=Badge_Grade) | |
| [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) | |
| [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4540752.svg)](https://doi.org/10.5281/zenodo.4540752) | [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4231/badge)](https://bestpractices.coreinfrastructure.org/projects/4231) | |
| [![status](https://joss.theoj.org/papers/2575e93fc693c5c3bfa8736c60c35398/status.svg)](https://joss.theoj.org/papers/2575e93fc693c5c3bfa8736c60c35398) | | |

## Introduction

The `mf2` package provides consistent, efficient and tested Python
implementations of a variety of multi-fidelity benchmark functions. The goal is
to simplify life for numerical optimization researchers by saving time otherwise
spent reimplementing and debugging the same common functions, and enabling
direct comparisons with other work using the same definitions, improving
reproducibility in general.

A multi-fidelity function usually reprensents an objective which should be
optimized. The term 'multi-fidelity' refers to the fact that multiple versions
of the objective function exist, which differ in the accuracy to describe the
real objective. A typical real-world example would be the aerodynamic
efficiency of an airfoil, e.g., its drag value for a given lift value. The
different fidelity levels are given by the accuracy of the evaluation method
used to estimate the efficiency. Lower-fidelity versions of the objective
function refer to less accurate, but simpler approximations of the objective,
such as computational fluid dynamic simulations on rather coarse meshes,
whereas higher fidelity levels refer to more accurate but also much more
demanding evaluations such as prototype tests in wind tunnels. The hope of
multi-fildelity optimization approaches is that many of the not-so-accurate but
simple low-fidelity evaluations can be used to achieve improved results on the
realistic high-fidelity version of the objective where only very few
evaluations can be performed.

The only dependency of the mf2 package is the `numpy` package.

Documentation is available at [mf2.readthedocs.io][docs]

## Installation

The recommended way to install `mf2` is with Python's `pip`:
```
python3 -m pip install --user mf2
```
or alternatively using `conda`:
```
conda install -c conda-forge mf2
```

For the latest version, you can install directly from source:
```
python3 -m pip install --user https://github.com/sjvrijn/mf2/archive/master.zip
```

To work in your own version locally, it is best to clone the repository first,
and additionally install the dev-requirements:
```
git clone https://github.com/sjvrijn/mf2.git
cd mf2
python3 -m pip install --user -e .[dev]
```

## Example Usage

```python
import mf2
import numpy as np

# set numpy random seed for reproducibility
np.random.seed(42)
# generate 5 random samples in 2D as matrix
X = np.random.random((5, 2))

# print high fidelity function values
print(mf2.branin.high(X))
# Out: array([36.78994906 34.3332972  50.48149005 43.0569396  35.5268224 ])

# print low fidelity function values
print(mf2.branin.low(X))
# Out: array([-5.8762639  -6.66852889  3.84944507 -1.56314141 -6.23242223])
```

For more usage examples, please refer to the full documentation on
[readthedocs][docs].

## Contributing

Contributions to this project such as bug reports or benchmark function
suggestions are more than welcome! Please [refer to ``CONTRIBUTING.md``] for more
details.

## Contact

The [Gitter][gitter] channel is the preferred way to get in touch for any other
questions, comments or discussions about this package.

## Citation

Was this package useful to you? Great! If this leads to a publication, we'd
appreciate it if you would cite our [JOSS paper]:

```
@article{vanRijn2020,
  doi = {10.21105/joss.02049},
  url = {https://doi.org/10.21105/joss.02049},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2049},
  author = {Sander van Rijn and Sebastian Schmitt},
  title = {MF2: A Collection of Multi-Fidelity Benchmark Functions in Python},
  journal = {Journal of Open Source Software}
}
```


[docs]:               https://mf2.readthedocs.io/en/latest/
[docs-badge]:         https://mf2.readthedocs.io/en/latest/?badge=latest
[gitter]:             https://gitter.im/pymf2/community
[gitter-badge]:       https://gitter.im/pymf2/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[new-issue]:          https://github.com/sjvrijn/mf2/issues/new
[pytest-regressions]: https://github.com/ESSS/pytest-regressions
[JOSS paper]:         https://joss.theoj.org/papers/10.21105/joss.02049
[refer to ``CONTRIBUTING.md``]: https://github.com/sjvrijn/mf2/blob/master/CONTRIBUTING.md

# Contributing

## Bugs
If you've found a problem of some sort, please open an issue on
[GitHub][new-issue] with your code, the resulting output and expected output.

## Additions
Suggestions for new functions are welcome too. When suggesting a function,
please include the DOI of the source paper.

When submitting a PR to add new functions to this package, you can roughly
follow the following steps:

 1. Implement the function in a new file in the appropriate (sub)folder
 2. Add it to the tests:
    * Add the function in the `tests/property_test.py` and
    `tests/regression_test.py` files
    * Run the tests: `pytest tests`. It will fail the first time while the
    [`pytest-regressions`][pytest-regressions] package automatically creates
    the new output files.
    * Run the tests again to confirm that all tests now pass.
 3. Make sure to commit all new and updated files to git (Travis-CI will
    complain otherwise ;)
 4. Create a pull-request!

If you need any help with this process, please get in touch as specified in the
README under **Contact**.

[new-issue]:          https://github.com/sjvrijn/mf2/issues/new
[pytest-regressions]: https://github.com/ESSS/pytest-regressions<!-- Feel free to remove check-list items aren't relevant to your change -->

<!-- New function addition checklist -->
 - [ ] New function has summary docstring with source paper DOI
 - [ ] Summary added to `changelog.txt`
 - [ ] Function added to `tests\regression_test.py` and `tests\property_test.py`
 - [ ] Regression test files added
 - [ ] Tests pass and coverage is not reduced---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Code example to reproduce the behavior:
```python
import mf2

...
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Version information:**
 - MF2 version:
 - Numpy version:
 - Python version:

**Additional context**
Add any other context about the problem here.
---
name: Function suggestion
about: Suggest a new benchmark function to include
title: "[Suggestion] <function(s)> from <author et al.>"
labels: New function
assignees: ''

---

**In which paper is this function introduced?**
"Title", authors and DOI.

**What function(s) should be added from this paper**
Short name(s) of the function(s), as used in the paper.

**What are characteristics of these functions?**
Indicate briefly for each function their dimensionality, relevant landscape features, and any other peculiarities.

**Additional information**
Add any other relevant information here
---
name: Feature request
about: Suggest an idea for this project
title: "[Feature request]"
labels: enhancement
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
---
title: 'MF2: A Collection of Multi-Fidelity Benchmark Functions in Python'
tags:
  - Python
  - optimization
  - benchmarks
authors:
  - name: Sander van Rijn
    orcid: 0000-0001-6159-041X
    affiliation: 1
  - name: Sebastian Schmitt
    affiliation: 2
affiliations:
  - name: Leiden University, The Netherlands
    index: 1
  - name: Honda Research Institute Europe, Germany
    index: 2
date: 7 April 2020
bibliography: paper.bib

---


# Summary

The field of (evolutionary) optimization algorithms often works with expensive
black-box optimization problems. However, for the development of novel
algorithms and approaches, real-world problems are not feasible due to their
high computational cost. Instead, benchmark functions such as Sphere, Rastrigin,
and Ackley are typically used. These functions are not only fast to compute, but
also have known properties which are very helpful when examining the performance
of new algorithms.

As only a limited set of benchmark functions are typically used in the literature,
compiling a set of implementations for the most commonly used functions is
warranted. This ensures correctness of the functions, makes any results directly
comparable, and simply saves researchers time from not having to implement the functions
themselves. An example of a commonly used benchmark suite for optimizing
continuous problems is the COCO BBOB software by @nikolaus_hansen:2019.

As simulation-based problems in engineering are requiring increasingly more time
and computation power, a new sub-field of *multi-fidelity* optimization has
gained popularity. A multi-fidelity problem is characterised by having multiple
versions of an evaluation function that differ in their accuracy of describing
the real objective. A real-world example would be the aerodynamic efficiency of
an airfoil: A *low-fidelity* simulation would use a coarse mesh, and give lower
accuracy, but be fast to calculate, while a *high-fidelity* simulation would use
a much finer mesh and therefore be more accurate while taking longer to
calculate. Multi-fidelity methods aim to combine these multiple information
sources to obtain better results in equal or less time compared to only using a
single information source.

To this end, new multi-fidelity benchmark functions have been introduced in the
literature and are being adopted by other researchers. These multi-fidelity
problems naturally benefit from having the different fidelities combined into a
single 'problem'. The existing single-fidelity benchmark suites that exist
cannot be used for this field: no existing suite of benchmark functions
currently uses such a structure, or can easily accomodate it. Therefore, this
new class of benchmark problems is best served by introducing a new
implementation suite because their structure is inherently different from other
benchmarks. A new suite additionally gives more freedom to adapt to new multi
fidelity benchmarks as the field continues to evolve and new needs become
apparent.

The ``MF2`` package provides a consistent Python implementation of a collection
of these Multi-Fidelity Functions. It uses a standard interface that allows for
querying single vectors or multiple row-vectors as a single matrix, relying on
``numpy``'s optimized back-end to handle parallelization. It also offers a
simple factory pattern interface for functions with parameters for e.g.
correlation and dimensionality. A plot of how these implementations scale can
be seen in \autoref{fig:scalability}.

At this moment, ``MF2`` has collected functions
from the following previous works:

  * @forrester:2007 introduced a simple 1D bi-fidelity function for mostly
    illustrative purposes.
  * @simulationlib:2017 have previously collected a small collection of
    MATLAB/R implementations for the Borehole, Currin and Park91 A and B
    functions.
  * @dong_multi-fidelity:2015 introduced bi-fidelity versions of the
    Bohachevsky, Booth, Branin, Himmelblau and Six-hump Camelback functions.
  * @toal_considerations:2015 introduced correlation-adjustable multi-fidelity
    versions of the Branin, Paciorek, Hartmann3 and Trid functions.

As no convenient existing implementations written in Python could be found
during the authors' research on how the accuracy of multi-fidelity surrogate
models depends on the number of samples per fidelity, which required the
evaluation of many independent model training and test sets, the decision was
made to standardize the implementations and make them available for wider use.


![**Scalability plot** This plot shows how the evaluation time of high- and
low-fidelity functions scales with the number of points *N* being passed in
simultaneously. Running times were measured on a desktop PC with an Intel core
i7 5820k 6-core CPU, with Python 3.6.3 and Numpy 1.18.4. The times are divided
by the time needed for N=1 as a normalization. This is done independently for
each function and fidelity level. Results are grouped by function
dimensionality. If there are multiple functions, the mean is plotted with error
bars indicating the minimum and maximum time. Note that the 3D Hartmann3 and 6D
Hartmann6 function are significantly more computationally expensive than other
functions by definition, as they requires multiple matrix multiplications.
\label{fig:scalability}](../docs/_static/scalability.pdf)

# Acknowledgements

This work is part of the research program DAMIOSO with project number
628.006.002, which is partly financed by the Netherlands Organisation
for Scientific Research (NWO).

The first author would like to thank dr. Matthijs van Leeuwen, prof. dr. Thomas
Bäck, and dr. Markus Olhofer for their supervision and involvement in the
DAMIOSO project.

# References
# Matlab files for scalability comparison

This folder contains the matlab files that were used for the scalability
comparison plots:

 - Utility functions:
 
   - `myscale.m` scales input from [0, 1] to the given new range.

   - `mytiming.m` creates random input of desired size `N` and performs a timing
     experiment akin to Python's [`timeit`].
     
 - Test functions:
 
   - `test_function_output.m` confirms that the matlab code gives the same
     result as the Python code for the same input. Expects input
     `input_<N>d.mat` and output `output_<N>d_<Name>_<Fidelity>.mat`-files to
     be available. These can be created using [`scipy.io.savemat`] in Python.

   - `test_function_performance.m` calls `mytiming.m` for all listed
     function-fidelity combinations, with increasing evaluation sizes `N`. By
     default, values for `N` are successive powers of 10: 10^0 -- 10^6.
     
 - Matlab implementations of benchmark functions

   The following implementations were retrieved from
   ``https://www.sfu.ca/~ssurjano/multi.html`` on 2017-10-02, and are available
   here under their original GNU GPL v2.0 licenses:
   
   - [Borehole]: `borehole.m` & `boreholelc.m`

   - [Currin]: `curretal88exp.m` & `curretal88explc.m`

   - [Park91a]: `park91a.m` & `park91alc.m`

   - [Park91b]: `park91b.m` & `park91blc.m`


[`timeit`]: https://docs.python.org/3/library/timeit.html#timeit.Timer.autorange
[`scipy.io.savemat`]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html#scipy.io.savemat
[Borehole]: https://www.sfu.ca/~ssurjano/borehole.html
[Currin]: https://www.sfu.ca/~ssurjano/curretal88exp.html
[Park91a]: https://www.sfu.ca/~ssurjano/park91a.html
[Park91b]: https://www.sfu.ca/~ssurjano/park91b.html
