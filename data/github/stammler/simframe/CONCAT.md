# Simframe

[![Documentation Status](https://readthedocs.org/projects/simframe/badge/?version=latest)](https://simframe.readthedocs.io/en/latest/?badge=latest) 
![GitHub](https://img.shields.io/github/license/stammler/simframe) 
[![status](https://joss.theoj.org/papers/0ef61e034c57445e846b2ec383c920a6/status.svg)](https://joss.theoj.org/papers/0ef61e034c57445e846b2ec383c920a6) 
![PyPI - Downloads](https://img.shields.io/pypi/dm/simframe?label=PyPI%20downloads)

### Framework for scientific simulations

`Simframe` is a Python framework to facilitate scientific simulations. The scope of the software is to provide a framework which can hold data fields, which can be used to integrate differential equations, and which can read and write data files.

Data fields are stored in modified `numpy.ndarray`s. Therefore, `Simframe` can only work with data, that can be stored in `NumPy` arrays.

## Installation

`pip install simframe`

## Documentation

[https://simframe.readthedocs.io/](https://simframe.readthedocs.io/)

* [1. Simple Integration](https://simframe.readthedocs.io/en/latest/1_simple_integration.html)
* [2. Advanced Integration](https://simframe.readthedocs.io/en/latest/2_advanced_integration.html)
* [3. Updating Groups and Fields](https://simframe.readthedocs.io/en/latest/3_updating.html)
* [4. Custom Integration Schemes](https://simframe.readthedocs.io/en/latest/4_custom_schemes.html)
* [5. Adaptive Integration Schemes](https://simframe.readthedocs.io/en/latest/5_adaptive_schemes.html)
* [6. Implicit Integration](https://simframe.readthedocs.io/en/latest/6_implicit_integration.html)
* [Example: Coupled Oscillators](https://simframe.readthedocs.io/en/latest/example_coupled_oscillators.html)
* [Example: Double Pendulum](https://simframe.readthedocs.io/en/latest/example_double_pendulum.html)
* [Example: Compartmental Models](https://simframe.readthedocs.io/en/latest/example_compartmental_models.html)

[Module Reference](https://simframe.readthedocs.io/en/latest/api.html)

## Contributing

To contribute to the software, please read the [contribution guidelines](https://github.com/stammler/simframe/blob/master/.github/CONTRIBUTING.md).

## Citation

If you use `Simframe`, please remember to cite [Stammler & Birnstiel (2022)](https://doi.org/10.21105/joss.03882).

```
@article{simframe,
  doi = {10.21105/joss.03882},
  url = {https://doi.org/10.21105/joss.03882},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {3882},
  author = {Sebastian M. Stammler and Tilman Birnstiel},
  title = {Simframe: A Python Framework for Scientific Simulations},
  journal = {Journal of Open Source Software}
}

```

## Ackowledgements

`simframe` has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme under grant agreement No 714769.

`simframe` was developed at the [University Observatory](https://www.usm.uni-muenchen.de/index_en.php) of the [Ludwig Maximilian University of Munich](https://www.en.uni-muenchen.de/index.html).
---
title: 'Simframe: A Python Framework for Scientific Simulations'
tags:
  - Python
  - Jupyter Notebook
  - NumPy
  - SciPy
authors:
  - name: Sebastian M. Stammler^[corresponding author]
    orcid: 0000-0002-1589-1796
    affiliation: 1
  - name: Tilman Birnstiel
    orcid: 0000-0002-1899-8783
    affiliation: 1, 2
affiliations:
  - name: University Observatory, Faculty of Physics, Ludwig-Maximilians-Universität München, Scheinerstr. 1, 81679 Munich, Germany
    index: 1
  - name: Exzellenzcluster ORIGINS, Boltzmannstr. 2, D-85748 Garching, Germany
    index: 2
date: 11 August 2021
bibliography: paper.bib
---

# Summary

`Simframe` is a Python framework to facilitate scientific simulations. The scope of the software is to provide a framework which can hold data fields, which can be used to integrate differential equations, and which can read and write data files.

Conceptually, upon initialization `Simframe` is an empty frame that can be filled with `Field`s containing data. `Field`s are derived from `numpy.ndarray`s [@Harris:2020], but with extended functionality. The user can then specify differential equations to those data fields and can set up an integrator which is integrating those fields according the given differential equations. Therefore, `Simframe` can only work with data, that can be stored in `NumPy` arrays.

Data fields that should not be integrated themselves, but are still required for the model, can have an update function assigned to them, according to which they will be updated once per integration step.

`Simframe` contains a number of integration schemes of different orders, both for explicit and implicit integration. Furthermore, `Simframe` includes methods to read and write output files.

Due to its modular structure, `Simframe` can be extended at will, for example, by implementing new integration schemes and/or user-defined output formats.

# Statement of need

Solving differential equations is part of the daily work of scientists. `Simframe` facilitates this by providing the infrastructure: Data structures, integration schemes, and methods to write and read output files.

On one hand, `Simframe` can be used to quickly solve small scientific problems, and, on the other hand, it can be easily extended to larger projects due to its versatility and modular structure.

Furthermore, `Simframe` is ideal for beginners without programming experience who are taking their first steps in solving differential equations. It can therefore be used to design lectures or practical courses at schools and universities, as it allows students to concentrate on the essentials without having to write larger programs on their own.

Plenty of ODE solver packages already exist for Python, like `solve_ivp` or `odeint` in `SciPy`'s `integrate` module, however, these do not provide data structures, nor input/output capabilities. `Simframe` offers a flexible framework to define, group, and describe data, define how it is updated, use existing integrators or define new ones, and to handle writing of data or serializing the entire simulation object, all in one modular package. Existing integrators like `solve_ivp` or `odeint` can be used within `Simframe` by simply adding them to an integration scheme.

# Features

## Data fields

The data fields of `Simframe` are subclassed `NumPy` `ndarray`s. The full `NumPy` functionality can therefore be used on `Simframe` data fields. The `ndarray`s have been extended to store additional information about differential equations or update functions and a string description of the field. The data fields can be arranged in groups to facilitate a clear structure within the data frame.

## Integration schemes

`Simframe` includes a number of basic integration schemes by default [@Hairer:1993]. All of the implemented schemes are Runge-Kutta methods of different orders. Some of the methods are adaptive, i.e., they are embedded Runge-Kutta methods, that return an optimal step size for the integration variable, such that the desired accuracy is achieved. The implicit methods require a matrix inversion that is either done directly by `NumPy` or by using the GMRES solver provided by `SciPy` [@Virtanen:2020].

Here is a list of all implemented integration schemes:

| Order | Scheme                      |          |          | solver |
| :---: | --------------------------- | :------: | :------: | ------ |
|   1   | Euler                       | explicit |          |        |
|   1   | Euler                       | implicit |          | direct |
|   1   | Euler                       | implicit |          | GMRES  |
|   2   | Fehlberg                    | explicit | adaptive |        |
|   2   | Heun                        | explicit |          |        |
|   2   | Heun-Euler                  | explicit | adaptive |        |
|   2   | midpoint                    | explicit |          |        |
|   2   | midpoint                    | implicit |          | direct |
|   2   | Ralston                     | explicit |          |        |
|   3   | Bogacki-Shampine            | explicit | adaptive |        |
|   3   | Gottlieb-Shu                | explicit | adaptive |        |
|   3   | Heun                        | explicit |          |        |
|   3   | Kutta                       | explicit |          |        |
|   3   | Ralston                     | explicit |          |        |
|   3   | Strong Stability Preserving | explicit |          |        |
|   4   | 3/8 rule                    | explicit |          |        |
|   4   | Ralston                     | explicit |          |        |
|   4   | Runge-Kutta                 | explicit |          |        |
|   5   | Cash-Karp                   | explicit | adaptive |        |
|   5   | Dormand-Prince              | explicit | adaptive |        |

## I/O

By default `Simframe` has two options for storing simulation results. One is by storing the data in a separate namespace within the `Simframe` object itself, useful for small simulations to access results without writing/reading data files. Another one is by storing the data in HDF5 data files using the `h5py` package [@Collette:2014].

If configured by the user, `Simframe` is writing dump files, from which the simulation can be resumed, in case the program crashed unexpectedly. These dump files are serialized `Simframe` objects using the `dill` package [@McKerns:2012].

# Acknowledgements

The authors acknowledge funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme under grant agreement No 714769 and funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under grants 361140270, 325594231, and Germany's Excellence Strategy - EXC-2094 - 390783311.

# References# How to contribute?

We invite everyone to contribute to the further development of `simframe`. Either by directly contributing to the code, by reporting bugs, or by requesting features.

## Contributing code

To contribute code please open a new pull request and describe the changes to the software your pull request introduces.

Please note, that we want to achieve a **code coverage of 100%**. Any addition to the software must therefore also come with unit tests. Additional features must also be described in the documentation.

## Reporting bugs

If you encountered a bug in `simframe`, please open a new [bug report issue](https://github.com/stammler/simframe/issues/new?template=bug_report.md&title=[BUG]+Descriptive+title+of+the+bug+report) and describe the bug, the expected behavior, and steps how to reproduce the bug.

## Requesting features

If you have an idea of a new feature, that is missing in `simframe`, or if you want to suggest an improvement, please open a new [feature request issue](https://github.com/stammler/simframe/issues/new?template=feature_request.md&title=[FEATURE]+Descriptive+title+of+the+feature+request).
# Pull Requests

Please have a look the [contribution guidelines](https://github.com/stammler/simframe/blob/master/.github/CONTRIBUTING.md).

## Reference issue

Please reference the issue this pull request is referring to, if applicable.

## Describe the changes

Please describe the changes or fixes this pull request introduces into the project.

## Additional information

Please provide additional information that is helpful for this pull request.
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] Descriptive title of the bug report"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Describe the steps to reproduce the behavior.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Software versions**
List the versions of the software packages involved in reproducing the bug.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE] Descriptive title of the feature request"
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
