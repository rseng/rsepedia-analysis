<img src="https://raw.githubusercontent.com/yaml2sbml-dev/yaml2sbml/main/doc/logo/logo_yaml2sbml_long.svg" alt="yaml2sbml logo"/>

[![CI](https://github.com/yaml2sbml-dev/yaml2sbml/workflows/CI/badge.svg)](https://github.com/yaml2sbml-dev/yaml2sbml/actions)
[![codecov](https://codecov.io/gh/yaml2sbml-dev/yaml2sbml/branch/master/graph/badge.svg)](https://codecov.io/gh/yaml2sbml-dev/yaml2sbml)
[![Documentation Status](https://readthedocs.org/projects/yaml2sbml/badge/?version=latest)](https://yaml2sbml.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/632acdc8d4ef4f50bf69892b8862fd24)](https://www.codacy.com/gh/yaml2sbml-dev/yaml2sbml/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=yaml2sbml-dev/yaml2sbml&amp;utm_campaign=Badge_Grade)
[![PyPI](https://badge.fury.io/py/yaml2sbml.svg)](https://badge.fury.io/py/yaml2sbml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03215/status.svg)](https://doi.org/10.21105/joss.03215)


## Table of contents

* [About](#about)

* [Installation](#installation)

* [Getting Started](#getting-started)

* [Basic Syntax](#basic-syntax)

* [How to cite](#how-to-cite)

* [Contact](#contact)


## About

`yaml2sbml` is a small package in Python to convert an ODE model specified in a YAML file into 
[**SBML**](http://www.sbml.org/) for ODE simulation and into 
[**PEtab**](https://github.com/PEtab-dev/PEtab/) for parameter fitting. `yaml2sbml` offers:


* a translator of ODE models specified in YAML into SBML/PEtab via a Python and a command-line interface;
* a format validator for the input YAML; and
* a model editor, which provides a simplified interface to generate, import and export YAML models.

## Installation

`yaml2sbml` can be installed from PyPI with

```shell
pip install yaml2sbml
```
For more info see the [docs](https://yaml2sbml.readthedocs.io/en/latest/).

## Getting Started

* The documentation of the tool is hosted on [Read the Docs](https://yaml2sbml.readthedocs.io/en/latest/).
* The [format documentation](https://yaml2sbml.readthedocs.io/en/latest/format_specification.html) describes the input YAML. 

* Jupyter notebooks containing examples can be found under [`doc/examples`](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples).  Most notably:
    * [Lotka_Volterra.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_python/Lotka_Volterra.ipynb) showing the Python package,
    * [Lotka_Volterra_CLI.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_CLI/Lotka_Volterra_CLI.ipynb) showing the command-line interface, and
    * [Lotka_Volterra_Model_Editor.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_Model_Editor/Lotka_Volterra_Model_Editor.ipynb) demonstrates the Model Editor.

## Basic Syntax

### Python

A YAML model can be translated to SBML/PEtab in Python via
```python
import yaml2sbml

# SBML conversion
yaml2sbml.yaml2sbml(yaml_dir, sbml_dir)

#PEtab conversion
yaml2sbml.yaml2petab(yaml_dir, 
                     output_dir,
                     sbml_name)
```
### Command-Line Interface

and in the command-line via 
```shell
# SBML conversion
yaml2sbml <yaml_dir> <sbml_dir>

#PEtab conversion
yaml2petab <yaml_dir> <output_dir> <sbml_name>
```

### Format Validation

Format validation is possible in Python via `yaml2sbml.validate_yaml` and in the command-line via `yaml2sbml_validate`.

## How to cite

`yaml2sbml` is published in the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03215). 

When using `yaml2sbml` in your project, please cite

* Vanhoefer J., Matos, M. R. A., Pathirana, D., Schälte, Y. and Hasenauer, J. (2021). yaml2sbml: Human-readable and -writable specification of ODE models and their conversion to SBML. Journal of Open Source Software, 6(61), 3215, https://doi.org/10.21105/joss.03215


```
@article{Vanhoefer2021,
  doi = {10.21105/joss.03215},
  url = {https://doi.org/10.21105/joss.03215},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {61},
  pages = {3215},
  author = {Jakob Vanhoefer and Marta R. A. Matos and Dilan Pathirana and Yannik Schälte and Jan Hasenauer},
  title = {yaml2sbml: Human-readable and -writable specification of ODE models and their conversion to SBML},
  journal = {Journal of Open Source Software}
}
```
## Contact
If you have a question regarding the tool: Please drop us an [issue](https://github.com/yaml2sbml-dev/yaml2sbml/issues/new) or a [mail](mailto:jakob.vanhoefer@uni-bonn.de), we are happy to help.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**What are the steps to reproduce the bug?**
A clear and concise description of how to reproduce the bug.

**Which output did you expect?**
A clear and concise description of what you expected to happen.

**What happened instead?**
A clear and concise description of what you observed instead.

**What could be a solution?**
If possible, a clear and concise description of a possible solution to the problem.

**Environment**
 - Operating system, e.g. Ubuntu 20.04, macOS 11.2.1  
 - `yaml2sbml` version (`pip show yaml2sbml`), e.g. 0.2.1
 - Python version (`python3 -V`)
 - If relevant: versions of other tools involved (e.g. `libsbml)

**Additional context**
Add any other context about the problem here.
# yaml2sbml examples

The [`yaml2sbml`](https://github.com/yaml2sbml-dev/yaml2sbml) package translates ODEs specified in YAML into [SBML](http://sbml.org/) for model simulation and [PEtab](https://github.com/PEtab-dev/PEtab) for parameter fitting. This folder contains several example notebooks for the usage of `yaml2sbml`.

## Scope of the Notebooks

*   _Lotka Volterra Notebooks_
    *   [Lotka Volterra (Python)](./Lotka_Volterra/Lotka_Volterra_python/Lotka_Volterra.ipynb)
        *   Introduces the input format, syntax & capabilities of `yaml2sbml`
        *   Showcases ODE simulation & parameter fitting using [AMICI](https://github.com/AMICI-dev/AMICI) & [pyPESTO](https://github.com/ICB-DCM/pyPESTO)
    *   [Lotka Volterra CLI](./Lotka_Volterra/Lotka_Volterra_CLI/Lotka_Volterra_CLI.ipynb)
        *   Shows the command line interface of `yaml2sbml`
    *   [Lotka Volterra Model Editor](./Lotka_Volterra/Lotka_Volterra_Model_Editor/Lotka_Volterra_Model_Editor.ipynb)
        *   Shows model construction using the model editor of `yaml2sbml`
*   [Sorensen](./Sorensen/yaml2sbml_Sorensen.ipynb)
    *   Application example, model of glucose and insulin metabolism
    *   Shows the extension of an existing YAML model by the model editor
*   [Finite State Projection](./Finite_State_Projection/Finite_State_Projection.ipynb)
    *   Application example of a stochastic gene transcription model having hundreds of states.
    *   Generates a complex and realistic model within a few lines of Python.

*   [Format Features](./Format_Features/Format_Features.ipynb)
    *   Several didactic examples, that show individual features of `yaml2sbml`:
        *   Time-dependent right hand sides
        *   Step functions in the right hand side
        *   Function definitions
                
## Running the examples

The notebooks come with additional dependencies, mainly for running the ODE-simulation and parameter fitting. Information on the installation of AMICI is given in its [installation guide](http://sbml.org/Special/Software/libSBML/docs/formatted/python-api/libsbml-math.html). Further dependencies ( e.g. [pyPESTO](https://github.com/ICB-DCM/pyPESTO)) can be installed via `pip install yaml2sbml[examples]`.
