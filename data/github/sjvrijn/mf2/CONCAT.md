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
.. _getting_started:

Getting Started
===============

This page contains some explained examples to help get you started with using
the ``mf2`` package.

The Basics: What's in a MultiFidelityFunction?
----------------------------------------------

This package serves as a collection of functions with multiple fidelity levels.
The number of levels is at least two, but differs by function. Each function is
encoded as a :class:`~mf2.multiFidelityFunction.MultiFidelityFunction`
with the following attributes:

``.name``
    The *name* is simply a standardized format of the name as an attribute to
    help identify which function is being represented [#footnote_name]_ .

``.ndim``
    Number of dimensions. This is the dimensionality (i.e. length) of the input
    vector X of which the objective is evaluated.

``.fidelity_names``
    This is a list of the human-readable names given to each fidelity.

``.u_bound``, ``.l_bound``
    The upper and lower bounds of the search-space for the function.

``.functions``
    A list of the actual function references. You won't typically need this
    list though, as will be explained next in :ref:`accessing_functions`.



Simple Usage
------------

.. _accessing_functions:

Accessing the functions
^^^^^^^^^^^^^^^^^^^^^^^

As an example, we'll use the :mod:`~mf2.booth` function. As we can see using
``.ndim`` and the bounds, it is two-dimensional:

    >>> from mf2 import booth
    >>> print(booth.ndim)
    2
    >>> print(booth.l_bound, booth.u_bound)
    [-10. -10.] [10. 10.]

Most multi-fidelity functions in ``mf2`` are *bi-fidelity* functions, but a
function can have any number of fidelities. A bi-fidelity function has two
fidelity levels, which are typically called ``high`` and ``low``. You can
easily check the names of the fidelities by printing the ``fidelity_names``
attribute of a function:

    >>> print(len(booth.fidelity_names))
    2
    >>> print(booth.fidelity_names)
    ['high', 'low']

These are just the names of the fidelities. The functions they represent can be
accessed as an object-style *attribute*,

    >>> print(booth.high)
    <function booth_hf at 0x...>

as a dictionary-style *key*,

    >>> print(booth['low'])
    <function booth_lf at 0x...>

or with a list-style *index* (which just passes through to ``.functions``).

    >>> print(booth[0])
    <function booth_hf at 0x...>
    >>> print(booth[0] is booth.functions[0])
    True

The object-style notation ``function.fidelity()`` is recommended for explicit
access, but the other notations are available for more dynamic usage. With the
list-style access, the *highest* fidelity is always at index *0*.


Calling the functions
^^^^^^^^^^^^^^^^^^^^^

All functions in the ``mf2`` package assume *row-vectors* as input. To evaluate
the function at a  single point, it can be given as a simple Python list or 1D
numpy array. Multiple points can be passed to the function individually, or
combined into a 2D list/array. The output of the function will always be
returned as a 1D numpy array:

    >>> X1 = [0.0, 0.0]
    >>> print(booth.high(X1))
    [74.]
    >>> X2 = [
    ...     [ 1.0,  1.0],
    ...     [ 1.0, -1.0],
    ...     [-1.0,  1.0],
    ...     [-1.0, -1.0]
    ... ]
    >>> print(booth.high(X2))
    [ 20.  80.  72. 164.]


Using the bounds
^^^^^^^^^^^^^^^^

Each function also has a given upper and lower bound, stored as a 1D numpy
array. They will be of the same length, and exactly as long as the
dimensionality of the function [#footnote_ndim]_ .

Below is an example function to create a uniform sample within the bounds::

    import numpy as np

    def sample_in_bounds(func, n_samples):
        raw_sample = np.random.random((n_samples, func.ndim))

        scale = func.u_bound - func.l_bound
        sample = (raw_sample * scale) + func.l_bound

        return sample


Kinds of functions
------------------

Fixed Functions
^^^^^^^^^^^^^^^

The majority of multi-fidelity functions in this package are 'fixed' functions.
This means that everything about the function is fixed:

* dimensionality of the input
* number of fidelity levels
* relation between the different fidelity levels

Examples of these functions include the 2D :mod:`~mf2.booth` and 8D
:mod:`~mf2.borehole` functions.


Dynamic Dimensionality Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some functions are dynamic in the dimensionality of the input they accept. An
example of such a function is the ``forrester`` function. The regular 1D
function is included as ``mf2.forrester``, but a custom n-dimensional version
can be obtained by calling the factory::

    forrester_4d = mf2.Forrester(ndim=4)

This ``forrester_4d`` is then a regular fixed function as seen before.


Adjustable Functions
^^^^^^^^^^^^^^^^^^^^

Other functions have a tunable parameter that can be used to adjust the
correlation between the different high and low fidelity levels. For these too,
you can simply call a factory that will return a version of that function with
the parameter fixed to your specification::

    paciorek_high_corr = mf2.adjustable.paciorek(a2=0.1)

The exact relationship between the input parameter and resulting correlation
can be found in the documentation of the specific functions. See for example
:mod:`~mf2.adjustable.paciorek`.

Adding Your Own
---------------

Each function is stored as a ``MultiFidelityFunction``-object, which contains
the dimensionality, intended upper/lower bounds, and of course all fidelity
levels. This class can also be used to define your own multi-fidelity function.

To do so, first define regular functions for each fidelity. Then create the
``MultiFidelityFunction`` object by passing a name, the upper and lower bounds,
and a tuple of the functions for the fidelities.

The following is an example for a 1-dimensional multi-fidelity function named
``my_mf_sphere`` with three fidelities::

    import numpy as np
    from mf2 import MultiFidelityFunction

    def sphere_hf(x):
        return x*x

    def sphere_mf(x):
        return x * np.sqrt(x) * np.sign(x)

    def sphere_lf(x):
        return np.abs(x)

    my_mf_sphere = MultiFidelityFunction(
        name='sphere',
        u_bound=[1],
        l_bound=[-1],
        functions=(sphere_hf, sphere_mf, sphere_lf),
    )

These functions can be accessed using list-style *indices*, but as no names
are given, the object-style *attributes* or dict-style *keys* won't work:

    >>> print(my_mf_sphere[0])
    <function sphere_hf at 0x...>
    >>> print(my_mf_sphere['medium'])
    ---------------------------------------------------------------------------
    IndexError                                Traceback (most recent call last)
    ...
    IndexError: Invalid index 'medium'
    >>> print(my_mf_sphere.low)
    ---------------------------------------------------------------------------
    AttributeError                            Traceback (most recent call last)
    ...
    AttributeError: 'MultiFidelityFunction' object has no attribute 'low'
    >>> print(my_mf_sphere.fidelity_names)
    None

To enable access by attribute or key, a tuple containing a name for each fidelity
is required. Let's extend the previous example by adding
``fidelity_names=('high', 'medium', 'low')``::

    my_named_mf_sphere = MultiFidelityFunction(
        name='sphere',
        u_bound=[1],
        l_bound=[-1],
        functions=(sphere_hf, sphere_mf, sphere_lf),
        fidelity_names=('high', 'medium', 'low'),
    )

Now we the attribute and key access will work:

    >>> print(my_named_mf_sphere[0])
    <function sphere_hf at 0x...>
    >>> print(my_named_mf_sphere['medium'])
    <function sphere_mf at 0x...>
    >>> print(my_named_mf_sphere.low)
    <function sphere_lf at 0x...>
    >>> print(my_named_mf_sphere.fidelity_names)
    ('high', 'medium', 'low')




.. rubric:: Footnotes

.. [#footnote_name] This is as they're instances of MultiFidelityFunction instead
                    of separate classes.

.. [#footnote_ndim] In fact, ``.ndim`` is defined as ``len(self.u_bound)``
Performance
===========

Where possible, all functions are written using `numpy <https://numpy.org/>`_
to make use of optimized routines and vectorization. Evaluating a single
point typically takes less than 0.0001 seconds on a modern desktop system,
regardless of function. This page shows some more detailed information about the
performance, even though this library should not be a bottleneck in any
programs.

The scripts for generating following performance overviews can be found in the
`docs/scripts <https://github.com/sjvrijn/mf2/tree/master/docs/scripts>`_ folder
of the repository. Presented running times were measured on a desktop PC with an
Intel Core i7 5820k 6-core CPU, with Python 3.6.3 and Numpy 1.18.4.

Performance Scaling
-------------------

The image below shows how the runtime scales as ``N`` points are passed to the
functions simultaneously as a matrix of size ``(N, ndim)``. Performance for the
high- and low-fidelity formulations are shown separately to give a fair
comparison: many low-fidelities are defined as computations on top of the
high-fidelity definitions. As absolute performance will vary per system, the
runtime is divided by the time needed for ``N=1`` as a normalization. This is
done independently for each function and fidelity level.

Up to ``N=1_000``, the time required scales less than linearly thanks to
efficient and vectorized numpy routines.

.. image:: ../_static/scalability.png
  :width: 640


Performance Comparison
----------------------

The following image shows how the scaling for the ``mf2`` implementation of the
**Currin**, **Park91A**, **Park91B** and **Borehole** functions compares to the
*Matlab* implementations by `Surjanovic and Bingham
<https://www.sfu.ca/~ssurjano/multi.html>`_, which can only evaluate one point
at a time, so do not use any vectorization. Measurements were performed using
*Matlab* version R2020a (9.8.0.1323502).

.. image:: ../_static/scalability_comparison.png
  :width: 640
.. _example_usage:

Example Usage
=============

This example is a reproduction of Figure 1 from http://doi.org/10.1098/rspa.2007.1900 :

The original figure:

.. image:: https://royalsocietypublishing.org/cms/asset/efa57e07-5384-4503-8b2b-ccbe632ffe87/3251fig1.jpg
  :width: 640

Code to reproduce the above figure as close as possible:

.. literalinclude:: ../scripts/example-usage.py
   :language: python
   :emphasize-lines: 13,14,16,34,35
   :linenos:

Reproduced figure:

.. image:: ../_static/recreating-forrester-2007.png
  :width: 640
.. Multi-Fidelity Functions documentation master file, created by
   sphinx-quickstart on Thu Nov 14 00:14:31 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MF2: Multi-Fidelity Functions
=============================

This is the documentation for the ``mf2`` package. For a short introduction with
examples, have a look at the :ref:`getting_started` page. Otherwise, you can
look at the available functions in the package by category.

The ``mf2`` package provides consistent, efficient and tested Python
implementations of a variety of multi-fidelity benchmark functions. The goal is
to simplify life for numerical optimization researchers by saving time otherwise
spent reimplementing and debugging the same common functions, and enabling
direct comparisons with other work using the same definitions, improving
reproducibility in general.

A multi-fidelity function usually reprensents an objective which should be
optimized. The term 'multi-fidelity' refers to the fact, that multiple versions
of the objective function exist which differ in the accuray to describe the
real objective. A typical real-world example would be the aerodynamic
efficiency of an airfoil, e.g., its drag value for a given lift value. The
different fidelity levels are given by the accuracy of the evaluation method
used to estimate the efficiency. Lower-fidelity versions of the objective
function refer to less accurate, but simpler approximations of the objective,
such as computational fluid dynamic simulations on rather coarse meshes,
whereas higher fidelity levels refer to more accurate but also much more
demaning evaluations such as prototype tests in wind tunnels. The hope of
multi-fildelity optimization approaches is that many of the not-so-accurate but
simple low-fidelity evaluations can be used to achieve improved results on the
realistic high-fidelity version of the objective where only very few
evaluations can be performed.

The only dependency of the mf2 package is the `numpy <https://numpy.org/>`_
package.

The source for this package is hosted at `github.com/sjvrijn/mf2 <https://github.com/sjvrijn/mf2>`_.

Last updated: (|today|)


Contents
========

.. toctree::
   :maxdepth: 2

   install
   example-usage
   performance
   getting-started
   mf2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Installation
============

The recommended way to install `mf2` is with Python's `pip`::

    python3 -m pip install --user mf2

or alternatively using `conda`::

    conda install -c conda-forge mf2


For the latest version, you can install directly from source::

    python3 -m pip install --user https://github.com/sjvrijn/mf2/archive/master.zip


To work in your own version locally, it is best to clone the repository first::

    git clone https://github.com/sjvrijn/mf2.git
    cd mf2
    python3 -m pip install --user -e .[dev]
mf2 package
===========

Fixed Functions
---------------
.. toctree::

   functions/multiFidelityFunction
   functions/bohachevsky
   functions/booth
   functions/borehole
   functions/branin
   functions/currin
   functions/forrester
   functions/hartmann
   functions/himmelblau
   functions/park91a
   functions/park91b
   functions/six_hump_camelback


Adjustable Functions
--------------------
.. toctree::

   functions/adjustable/branin
   functions/adjustable/paciorek
   functions/adjustable/hartmann
   functions/adjustable/trid
Forrester
=========

.. automodule:: mf2.forrester
    :members:
    :undoc-members:
    :show-inheritance:
MultiFidelityFunction
=====================

.. automodule:: mf2.multiFidelityFunction
    :members:
    :special-members: __init__
    :undoc-members:
    :show-inheritance:
Park91 B
========

.. automodule:: mf2.park91b
    :members:
    :undoc-members:
    :show-inheritance:
Hartmann6
=========

.. automodule:: mf2.hartmann
    :members:
    :undoc-members:
    :show-inheritance:
Six-Hump Camelback
==================

.. automodule:: mf2.six_hump_camelback
    :members:
    :undoc-members:
    :show-inheritance:
Bohachevsky
===========

.. automodule:: mf2.bohachevsky
    :members:
    :undoc-members:
    :show-inheritance:
Booth
=====

.. automodule:: mf2.booth
    :members:
    :undoc-members:
    :show-inheritance:
Borehole
========

.. automodule:: mf2.borehole
    :members:
    :undoc-members:
    :show-inheritance:
Currin
======

.. automodule:: mf2.currin
    :members:
    :undoc-members:
    :show-inheritance:
Himmelblau
==========

.. automodule:: mf2.himmelblau
    :members:
    :undoc-members:
    :show-inheritance:
Park91 A
========

.. automodule:: mf2.park91a
    :members:
    :undoc-members:
    :show-inheritance:
Branin
======

.. automodule:: mf2.branin
    :members:
    :undoc-members:
    :show-inheritance:
Adjustable Trid
===============

.. automodule:: mf2.adjustable.trid
    :members:
    :undoc-members:
    :show-inheritance:
Adjustable Hartmann3
====================

.. automodule:: mf2.adjustable.hartmann
    :members:
    :undoc-members:
    :show-inheritance:
Adjustable Paciorek
===================

.. automodule:: mf2.adjustable.paciorek
    :members:
    :undoc-members:
    :show-inheritance:
Adjustable Branin
=================

.. automodule:: mf2.adjustable.branin
    :members:
    :undoc-members:
    :show-inheritance:
