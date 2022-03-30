---
title: 'PyNumDiff: A Python package for numerical differentiation of noisy time-series data'
tags:
  - Python
  - numerical differentiation
  - denoising
  - dynamics
  - time series
  - machine learning
authors:
  - name: Floris Van Breugel^[corresponding author]
    affiliation: 1
  - name: Yuying Liu
    affiliation: 2
  - name: Bingni W. Brunton
    affiliation: 3
  - name: J. Nathan Kutz
    affiliation: 2
affiliations:
 - name: Department of Mechanical Engineering, University of Nevada at Reno
   index: 1
 - name: Department of Applied Mathematics, University of Washington
   index: 2
 - name: Department of Biology, University of Washington
   index: 3
date: 10 July 2021
bibliography: paper.bib
---

# Statement of need

The numerical computation of derivatives is ubiquitous in every scientific discipline and engineering application because derivatives express fundamental relationships among many quantities of interest. As a result, a large number of diverse algorithms have been developed to differentiate numerical data.  These efforts are challenging because, in reality, practitioners  often have sparse and noisy measurements and data, which undermine the ability to estimate accurate derivatives.  Among the diversity of mathematical approaches that have been formulated, many are ad hoc in nature and require significant bespoke tuning of multiple parameters to produce reasonable results. Thus, at a practical level, it is often unclear which method should be used, how to choose parameters, and how to compare results from different methods. 

Regardless of application domain, scientists of various levels of mathematical expertise would benefit from a unified toolbox for differentiation techniques and parameter tuning. To address these needs, we built the open-source package `PyNumDiff`, with two primary goals in mind: (1) to develop a unified source for a diversity of differentiation methods using a common API, and (2) to provide an objective approach for choosing optimal parameters with a single universal hyperparameter (`gamma`) that functions similarly for all differentiation methods [@van2020numerical]. By filling these needs, `PyNumdiff` facilitates easy computations of derivatives on diverse time-series data sets.

# State of the field

Currently, practitioners in need of numerical differentiation tools must often implement a number of methods themselves, before selecting one that is appropriate for their application. High-quality data can leverage computationally efficient and algorithmically simple methods such as the finite-difference, as implemented by standard packages such as NumPy [@harris2020array], SciPy [@scipy; @Virtanen2020], or specialized packages like findiff [@findiff]. Data that are sparse and noisy, however, require more sophisticated algorithms that pracitioners must build themselves based on routines implemented across modules found in disparate packages such as SciPy, PyKalman [@pykalman], PyDMD [@demo18pydmd], or stand alone scripts such as these [implementations of total variation regularization](https://sites.google.com/site/dnartrahckcir/home/tvdiff-code) [@rudin1992nonlinear; @chartrand2011numerical]. At present, there is no centralized repository that offers a diverse range of vetted numerical differentiation tools under a unified API in Python, or other software languages. 

	
# Summary

`PyNumDiff` is a Python package that implements methods for computing numerical derivatives of noisy data. 
In this package, we implement four commonly used families of differentiation methods whose mathematical formulations have different 
underlying assumptions, including both global and local methods [@ahnert2007numerical]. The first family of methods usually start by 
applying a smoothing filter to the data, followed by a finite difference calculation [@butterworth1930theory]. 
The second family relies on building a local model of the data through linear regression, and then analytically 
calculating the derivative based on the model [@belytschko1996meshless; @schafer2011savitzky; @savitzky1964smoothing]. 
The third family we consider is the Kalman filter [@kalman1960new; @henderson2010fundamentals; @aravkin2017generalized; @crassidis2004optimal], 
with unknown noise and process characteristics.   The last family is an optimization approach based on total variation 
regularization (TVR) method [@rudin1992nonlinear; @chartrand2011numerical]. For more technical details, 
refer to @van2020numerical. Individual methods under each family are accessed through the API as `pynumdiff.family.method`. 

Applying `PyNumDiff` usually 
takes three steps: (i) pick a differentiation method, (ii) obtain optimized parameters, and (iii) apply the differentiation. 
Step (ii) can be skipped if one wants to manually assign the parameters, which is recommended when computation time is limited and the timeseries is long. Alternatively for long timeseries, optimal parameters can be chosen using a short but representative subset of the data. This optimization routine is provided as a sub-module (pynumdiff.optimize) with the same structure of differentiation families (i.e. `pynumdiff.optimize.family.method`). By default, the package performs the optimization using the open source CVXOPT package. Faster solutions can be achieved by using proprietary solvers such as MOSEK. 

The software package includes tutorials in the form of Jupyter notebooks. These tutorials demonstrate the usage of the aforementioned
features. For more detailed information, there is a more comprehensive Sphinx documentation associated with the repository.

# Acknowledgements

The work of J. Nathan Kutz was supported by the Air Force Office of Scientific Research under Grant FA9550-19-1-0011 and FA9550-19-1-0386. The work of F. van Breugel was supported by NIH grant P20GM103650, Air Force Research Lab award FA8651-20-1-0002 Airforce Office of Scientific Research FA9550-21-0122. BWB acknowledges support from the Air Force Office of Scientific Research award FA9550-19-1-0386.

# References

# PyNumDiff

Python methods for numerical differentiation of noisy data, including multi-objective optimization routines for automated parameter selection.

<p align="center">
  <a href="https://pynumdiff.readthedocs.io/en/master/" target="_blank" >
    <img alt="Python for Numerical Differentiation of noisy time series data" src="docs/source/_static/logo_PyNumDiff.png" width="300" height="200" />
  </a>
</p>

<p align="center">
    <a href="#travis" alt="Travis Build Status">
        <img src="https://travis-ci.com/florisvb/PyNumDiff.svg?branch=master"/></a>
    <a href='https://pynumdiff.readthedocs.io/en/master/?badge=master'>
        <img src='https://readthedocs.org/projects/pynumdiff/badge/?version=master' alt='Documentation Status' /></a>
    <a href="https://badge.fury.io/py/pynumdiff">
        <img src="https://badge.fury.io/py/pynumdiff.svg" alt="PyPI version" height="18"></a>
    <a href="https://zenodo.org/badge/latestdoi/159711175">
        <img src="https://zenodo.org/badge/159711175.svg" alt="DOI"></a>
    <a href="https://joss.theoj.org/papers/102257ee4b0142bf49bc18d7c810e9d5">
        <img src="https://joss.theoj.org/papers/102257ee4b0142bf49bc18d7c810e9d5/status.svg"></a>
</p>

## Table of contents
* [Introduction](#introduction)
* [Structure](#structure)
* [Getting Started](#getting-started)
    * [Prerequisite](#prerequisite)
    * [Installing](#installing)
* [Usage](#usage)
    * [Basic usages](#basic-usages)
    * [Notebook examples](#notebook-examples)
    * [Important notes](#important-notes)
    * [Running the tests](#running-the-tests)
* [Citation](#citation)
* [License](#license)

## Introduction

PyNumDiff is a Python package that implements various methods for computing numerical derivatives of noisy data, which 
can be a critical step in developing dynamic models or designing control. There are four different families of methods 
implemented in this repository: smoothing followed by finite difference calculation, local approximation with linear 
models, Kalman filtering based methods and total variation regularization methods. Most of these methods have multiple
parameters involved to tune. We take a principled approach and propose a multi-objective optimization framework for 
choosing parameters that minimize a loss function to balance the faithfulness and smoothness of the derivative estimate.
For more details, refer to [this paper](https://doi.org/10.1109/ACCESS.2020.3034077).

## Structure

    PyNumDiff/
      |- README.md
      |- pynumdiff/
         |- __init__.py
         |- __version__.py
         |- finite_difference/
         |- kalman_smooth/
         |- linear_model/
         |- smooth_finite_difference/
         |- total_variation_regularization/
         |- utils/
         |- optimize/
            |- __init__.py
            |- __optimize__.py
            |- finite_difference/
            |- kalman_smooth/
            |- linear_model/
            |- smooth_finite_difference/
            |- total_variation_regularization/
         |- tests/
      |- examples
         |- 1_basic_tutorial.ipynb
         |- 2a_optimizing_parameters_with_dxdt_known.ipynb
         |- 2b_optimizing_parameters_with_dxdt_unknown.ipynb
      |- docs/
         |- Makefile
         |- make.bat
         |- build/
         |- source/
            |- _static
            |- _summaries
            |- conf.py
            |- index.rst
            |- ...
      |- setup.py
      |- .gitignore
      |- .travis.yml
      |- LICENSE.txt
      |- requirements.txt

## Citation


#### PyNumDiff python package:

@article{PyNumDiff2022,
  doi = {10.21105/joss.04078},
  url = {https://doi.org/10.21105/joss.04078},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {71},
  pages = {4078},
  author = {Floris Van Breugel and Yuying Liu and Bingni W. Brunton and J. Nathan Kutz},
  title = {PyNumDiff: A Python package for numerical differentiation of noisy time-series data},
  journal = {Journal of Open Source Software}
}


#### Optimization algorithm:

@article{ParamOptimizationDerivatives2020, 
doi={10.1109/ACCESS.2020.3034077}
author={F. {van Breugel} and J. {Nathan Kutz} and B. W. {Brunton}}, 
journal={IEEE Access}, 
title={Numerical differentiation of noisy data: A unifying multi-objective optimization framework}, 
year={2020}
}

## Getting Started

### Prerequisite

PyNumDiff requires common packages like `numpy`, `scipy`, `matplotlib`, `pytest` (for unittests), `pylint` 
(for PEP8 style check). For a full list, you can check the file [requirements.txt](requirements.txt)

In addition, it also requires certain additional packages for select functions, though these are not required for a successful install of PyNumDiff:
* Total Variation Regularization methods: [`cvxpy`](http://www.cvxpy.org/install/index.html)
* Linear Model Chebychev: [`pychebfun`](https://github.com/pychebfun/pychebfun/)

When using `cvxpy`, our default solver is set to be `MOSEK` (highly recommended), you would need to download their 
free academic license from their [website](https://www.mosek.com/products/academic-licenses/). Otherwise, you can also 
use other solvers which are listed [here](https://www.cvxpy.org/tutorial/advanced/index.html).

### Installing

The code is compatible with >=Python 3.5. It can be installed using pip or directly from the source code. Basic installation options include:

* From PyPI using pip: `pip install pynumdiff`. May require pre-installing `numpy, scipy, matplotlib`. 
* From source using pip git+: `pip install git+https://github.com/florisvb/PyNumDiff`
* From local source code using setup.py: requires pre-installing `numpy, scipy, matplotlib`. Then run `python ./setup.py install` from inside this directory. See below for example.

Installation of the optional packages such as `cvxpy` can be tricky because `cvxpy` requires pythonX-dev packages. Depending on your version of Ubuntu it can be challenging to meet all the right requirements and installation options (e.g. it is difficult to install python3.6-dev on Ubuntu 16.04). Here are several tested example installation workflows:

###### Complete install on Ubuntu 16.04 using python3.5 in blank virtual environment using pip git+:

```console
sudo apt-get install python3.5-dev
python3.5 -m venv ~/PYNUMDIFF35
source ~/PYNUMDIFF35/bin/activate
pip install --upgrade pip
pip install --upgrade pip
pip install git+https://github.com/florisvb/PyNumDiff
pip install git+https://github.com/pychebfun/pychebfun
pip install cvxpy
pip install git+http://github.com/MOSEK/Mosek.pip
```

###### Complete install on Ubuntu 18.04 using python3.6 in blank virtual environment using pip git+:

```console
sudo apt-get install python3.6-dev
python3.6 -m venv ~/PYNUMDIFF36
source ~/PYNUMDIFF36/bin/activate
pip install --upgrade pip
pip install git+https://github.com/florisvb/PyNumDiff
pip install git+https://github.com/pychebfun/pychebfun
pip install cvxpy
pip install Mosek
```

###### Complete install on Ubuntu 16.04 using python3.5 in blank virtual environment using setup.py:

```console
sudo apt-get install python3.5-dev
python3.5 -m venv ~/PYNUMDIFF35
source ~/PYNUMDIFF35/bin/activate
pip install --upgrade pip
pip install --upgrade pip
pip install numpy scipy matplotlib
python ./setup.py install
pip install git+https://github.com/pychebfun/pychebfun
pip install cvxpy
pip install git+http://github.com/MOSEK/Mosek.pip
```
<em>Note: If using the optional MOSEK solver for cvxpy you will also need a [MOSEK license](https://www.mosek.com/products/academic-licenses/), free academic license.</em>


## Usage

**PyNumDiff** uses [Sphinx](http://www.sphinx-doc.org/en/stable/) for code documentation.
So you can see more details about the API usage [there](https://pynumdiff.readthedocs.io/en/latest/).

### Basic usages

* Basic Usage: you provide the parameters
```bash
        x_hat, dxdt_hat = pynumdiff.sub_module.method(x, dt, params, options)     
```
* Intermediate usage: automated parameter selection through multi-objective optimization
```bash
        params, val = pynumdiff.optimize.sub_module.method(x, dt, params=None, 
                                                           tvgamma=tvgamma, # hyperparameter
                                                           dxdt_truth=None, # no ground truth data
                                                           options={})
        print('Optimal parameters: ', params)
        x_hat, dxdt_hat = pynumdiff.sub_module.method(x, dt, params, options={'smooth': True})`
```
* Advanced usage: automated parameter selection through multi-objective optimization using a user-defined cutoff frequency
```bash
        # cutoff_freq: estimate by (a) counting the number of true peaks per second in the data or (b) look at power spectra and choose cutoff
        log_gamma = -1.6*np.log(cutoff_frequency) -0.71*np.log(dt) - 5.1 # see: https://ieeexplore.ieee.org/abstract/document/9241009
        tvgamma = np.exp(log_gamma) 

        params, val = pynumdiff.optimize.sub_module.method(x, dt, params=None, 
                                                           tvgamma=tvgamma, # hyperparameter
                                                           dxdt_truth=None, # no ground truth data
                                                           options={})
        print('Optimal parameters: ', params)
        x_hat, dxdt_hat = pynumdiff.sub_module.method(x, dt, params, options={'smooth': True})`
```

### Notebook examples

We will frequently update simple examples for demo purposes, and here are currently exisiting ones:
* Differentiation with different methods: [1_basic_tutorial.ipynb](examples/1_basic_tutorial.ipynb)
* Parameter Optimization with known ground truth (only for demonstration purpose):  [2a_optimizing_parameters_with_dxdt_known.ipynb](examples/2a_optimizing_parameters_with_dxdt_known.ipynb)
* Parameter Optimization with unknown ground truth:  [2b_optimizing_parameters_with_dxdt_unknown.ipynb](./examples/2b_optimizing_parameters_with_dxdt_unknown.ipynb)


### Important notes

* Larger values of `tvgamma` produce smoother derivatives
* The value of `tvgamma` is largely universal across methods, making it easy to compare method results
* The optimization is not fast. Run it on subsets of your data if you have a lot of data. It will also be much faster with faster differentiation methods, like savgoldiff and butterdiff, and probably too slow for sliding methods like sliding DMD and sliding LTI fit. 
* The following heuristic works well for choosing `tvgamma`, where `cutoff_frequency` is the highest frequency content of the signal in your data, and `dt` is the timestep: `tvgamma=np.exp(-1.6*np.log(cutoff_frequency)-0.71*np.log(dt)-5.1)`


### Running the tests

We are using Travis CI for continuous intergration testing. You can check out the current status 
[here](https://travis-ci.com/github/florisvb/PyNumDiff).

To run tests locally, type:
```bash
> pytest pynumdiff
```


## License

This project utilizes the [MIT LICENSE](LICENSE.txt).
100% open-source, feel free to utilize the code however you like. 
finite_difference
=================

.. toctree::
   :maxdepth: 1

   finite_difference/_finite_differenceAPI documentation
============================

.. toctree::
   :maxdepth: 1

   finite_difference
   kalman_smooth
   linear_model
   smooth_finite_difference
   total_variation_regularization
   optimize
   utils

Contact
=======

Feel free to contact the author for any information.

- Floris van Breugel: fvanbreugel@unr.edu
- Yuying Liu: liuyuyingufo@gmail.comtotal_variation_regularization
==============================

.. toctree::
   :maxdepth: 1

   total_variation_regularization/_total_variation_regularizationkalman_smooth
=============

.. toctree::
   :maxdepth: 1

   kalman_smooth/_kalman_smoothHow to contribute
===================

We'd love to accept your patches and contributions to this project. There are
just a few small guidelines you need to follow.

Submitting a patch:

  1. It's generally best to start by opening a new issue describing the bug or feature you're intending to fix. Even if you think it's relatively minor, it's helpful to know what people are working on. Mention in the initial issue that you are planning to work on that bug or feature so that it can be assigned to you.

  2. Follow the normal process of forking the project, and setup a new branch to work in. It's important that each group of changes be done in separate branches in order to ensure that a pull request only includes the commits related to that bug or feature.

  3. To ensure properly formatted code, please make sure to use 4 spaces to indent the code. You should also run pylint over your code. It's not strictly necessary that your code be completely "lint-free", but this will help you find common style issues.

  4. Any significant changes should almost always be accompanied by tests. The project already has good test coverage, so look at some of the existing tests if you're unsure how to go about it. We're using coveralls that is an invaluable tools for seeing which parts of your code aren't being exercised by your tests.

  5. Do your best to have well-formed commit messages for each change. This provides consistency throughout the project, and ensures that commit messages are able to be formatted properly by various git tools.

  6. Finally, push the commits to your fork and submit a pull request. Please, remember to rebase properly in order to maintain a clean, linear git history.utils
=====

.. toctree::
   :maxdepth: 1

   utils/_pi_cruise_control
   utils/evaluate
   utils/simulate
   utils/utilityPyNumDiff
=========

Python methods for numerical differentiation of noisy data, including
multi-objective optimization routines for automated parameter selection.

Table of contents
-----------------

-  `Introduction <#introduction>`__
-  `Structure <#structure>`__
-  `Getting Started <#getting-started>`__

   -  `Prerequisite <#prerequisite>`__
   -  `Installing <#installing>`__

-  `Usage <#usage>`__

   -  `Basic usages <#basic-usages>`__
   -  `Notebook examples <#notebook-examples>`__
   -  `Important notes <#important-notes>`__
   -  `Running the tests <#running-the-tests>`__

-  `Citation <#citation>`__
-  `License <#license>`__

Introduction
------------

PyNumDiff is a Python package that implements various methods for
computing numerical derivatives of noisy data, which can be a critical
step in developing dynamic models or designing control. There are four
different families of methods implemented in this repository: smoothing
followed by finite difference calculation, local approximation with
linear models, Kalman filtering based methods and total variation
regularization methods. Most of these methods have multiple parameters
involved to tune. We take a principled approach and propose a
multi-objective optimization framework for choosing parameters that
minimize a loss function to balance the faithfulness and smoothness of
the derivative estimate. For more details, refer to `this
paper <https://doi.org/10.1109/ACCESS.2020.3034077>`__.

Structure
---------

::

    PyNumDiff/
      |- README.md
      |- pynumdiff/
         |- __init__.py
         |- __version__.py
         |- finite_difference/
         |- kalman_smooth/
         |- linear_model/
         |- smooth_finite_difference/
         |- total_variation_regularization/
         |- utils/
         |- optimize/
            |- __init__.py
            |- __optimize__.py
            |- finite_difference/
            |- kalman_smooth/
            |- linear_model/
            |- smooth_finite_difference/
            |- total_variation_regularization/
         |- tests/
      |- examples
         |- 1_basic_tutorial.ipynb
         |- 2a_optimizing_parameters_with_dxdt_known.ipynb
         |- 2b_optimizing_parameters_with_dxdt_unknown.ipynb
      |- docs/
         |- Makefile
         |- make.bat
         |- build/
         |- source/
            |- _static
            |- _summaries
            |- conf.py
            |- index.rst
            |- ...
      |- setup.py
      |- .gitignore
      |- .travis.yml
      |- LICENSE.txt
      |- requirements.txt

Getting Started
---------------

Prerequisite
~~~~~~~~~~~~

PyNumDiff requires common packages like ``numpy``, ``scipy``,
``matplotlib``, ``pytest`` (for unittests), ``pylint`` (for PEP8 style
check). For a full list, you can check the file
`requirements.txt <https://github.com/florisvb/PyNumDiff/blob/master/requirements.txt>`__

In addition, it also requires certain additional packages for select
functions, though these are not required for a successful install of
PyNumDiff: 

-  Total Variation Regularization methods: `cvxpy <http://www.cvxpy.org/install/index.html>`__ 
-  Linear Model Chebychev: `pychebfun <https://github.com/pychebfun/pychebfun/>`__

When using ``cvxpy``, our default solver is set to be ``MOSEK`` (highly
recommended), you would need to download their free academic license
from their
`website <https://www.mosek.com/products/academic-licenses/>`__.
Otherwise, you can also use other solvers which are listed
`here <https://www.cvxpy.org/tutorial/advanced/index.html>`__.

Installing
~~~~~~~~~~

The code is compatible with >=Python 3.5. It can be installed using pip
or directly from the source code. Basic installation options include:

-  From PyPI using pip: ``pip install pynumdiff``. May require
   pre-installing ``numpy, scipy, matplotlib``.
-  From source using pip git+:
   ``pip install git+https://github.com/florisvb/PyNumDiff``
-  From local source code using setup.py: requires pre-installing
   ``numpy, scipy, matplotlib``. Then run ``python ./setup.py install``
   from inside this directory. See below for example.

Installation of the optional packages such as ``cvxpy`` can be tricky
because ``cvxpy`` requires pythonX-dev packages. Depending on your
version of Ubuntu it can be challenging to meet all the right
requirements and installation options (e.g. it is difficult to install
python3.6-dev on Ubuntu 16.04). Here are several tested example
installation workflows:

Complete install on Ubuntu 16.04 using python3.5 in blank virtual environment using pip git+:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: console

    sudo apt-get install python3.5-dev
    python3.5 -m venv ~/PYNUMDIFF35
    source ~/PYNUMDIFF35/bin/activate
    pip install --upgrade pip
    pip install --upgrade pip
    pip install git+https://github.com/florisvb/PyNumDiff
    pip install git+https://github.com/pychebfun/pychebfun
    pip install cvxpy
    pip install git+http://github.com/MOSEK/Mosek.pip

Complete install on Ubuntu 18.04 using python3.6 in blank virtual environment using pip git+:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: console

    sudo apt-get install python3.6-dev
    python3.6 -m venv ~/PYNUMDIFF36
    source ~/PYNUMDIFF36/bin/activate
    pip install --upgrade pip
    pip install git+https://github.com/florisvb/PyNumDiff
    pip install git+https://github.com/pychebfun/pychebfun
    pip install cvxpy
    pip install Mosek

Complete install on Ubuntu 16.04 using python3.5 in blank virtual environment using setup.py:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: console

    sudo apt-get install python3.5-dev
    python3.5 -m venv ~/PYNUMDIFF35
    source ~/PYNUMDIFF35/bin/activate
    pip install --upgrade pip
    pip install --upgrade pip
    pip install numpy scipy matplotlib
    python ./setup.py install
    pip install git+https://github.com/pychebfun/pychebfun
    pip install cvxpy
    pip install git+http://github.com/MOSEK/Mosek.pip

Note: If using the optional MOSEK solver for cvxpy you will also need a
`MOSEK license <https://www.mosek.com/products/academic-licenses/>`__,
free academic license.

Usage
-----

Basic usages
~~~~~~~~~~~~

-  Basic Usage: you provide the parameters

   .. code:: bash

           x_hat, dxdt_hat = pynumdiff.sub_module.method(x, dt, params, options)     

-  Advanced usage: automated parameter selection through multi-objective
   optimization

   .. code:: bash

           params, val = pynumdiff.optimize.sub_module.method(x, dt, params=None, 
                                                              tvgamma=tvgamma, # hyperparameter
                                                              dxdt_truth=None, # no ground truth data
                                                              options={})
           print('Optimal parameters: ', params)
           x_hat, dxdt_hat = pynumdiff.sub_module.method(x, dt, params, options={'smooth': True})`

Notebook examples
~~~~~~~~~~~~~~~~~

-  Differentiation with different methods: `1\_basic\_tutorial.ipynb <https://github.com/florisvb/PyNumDiff/tree/master/examples/1_basic_tutorial.ipynb>`__ 
-  Parameter Optimization with known ground truth (only for demonstration purpose): `2a\_optimizing\_parameters\_with\_dxdt\_known.ipynb <https://github.com/florisvb/PyNumDiff/tree/master/examples/2a_optimizing_parameters_with_dxdt_known.ipynb>`__
-  Parameter Optimization with unknown ground truth: `2b\_optimizing\_parameters\_with\_dxdt\_unknown.ipynb <https://github.com/florisvb/PyNumDiff/tree/master/examples/2b_optimizing_parameters_with_dxdt_unknown.ipynb>`__


Important notes
~~~~~~~~~~~~~~~

-  Larger values of ``tvgamma`` produce smoother derivatives
-  The value of ``tvgamma`` is largely universal across methods, making
   it easy to compare method results
-  The optimization is not fast. Run it on subsets of your data if you
   have a lot of data. It will also be much faster with faster
   differentiation methods, like savgoldiff and butterdiff, and probably
   too slow for sliding methods like sliding DMD and sliding LTI fit.
-  The following heuristic works well for choosing ``tvgamma``, where
   ``cutoff_frequency`` is the highest frequency content of the signal
   in your data, and ``dt`` is the timestep:
   ``tvgamma=np.exp(-1.6*np.log(cutoff_frequency)-0.71*np.log(dt)-5.1)``

Running the tests
~~~~~~~~~~~~~~~~~

We are using Travis CI for continuous intergration testing. You can
check out the current status
`here <https://travis-ci.com/github/florisvb/PyNumDiff>`__.

To run tests locally, type:

.. code:: bash

    > pytest pynumdiff

Citation
--------

@ARTICLE{9241009, author={F. {van Breugel} and J. {Nathan Kutz} and B.
W. {Brunton}}, journal={IEEE Access}, title={Numerical differentiation
of noisy data: A unifying multi-objective optimization framework},
year={2020}, volume={}, number={}, pages={1-1},
doi={10.1109/ACCESS.2020.3034077}}

Developer's Guide
-----------------

.. toctree::
   :maxdepth: 1

   code
   contact
   contributing
   LICENSEoptimize
=========

.. toctree::
   :maxdepth: 1

   optimize/__optimize__
   optimize/__finite_difference__
   optimize/__kalman_smooth__
   optimize/__linear_model__
   optimize/__smooth_finite_difference__
   optimize/__total_variation_regularization__
smooth_finite_difference
========================

.. toctree::
   :maxdepth: 1

   smooth_finite_difference/_smooth_finite_differencelinear_model
============

.. toctree::
   :maxdepth: 1

   linear_model/_linear_model_kalman_smooth
================

.. currentmodule:: pynumdiff.kalman_smooth._kalman_smooth

.. automodule:: pynumdiff.kalman_smooth._kalman_smooth
    :members:_linear_model
================

.. currentmodule:: pynumdiff.linear_model._linear_model

.. automodule:: pynumdiff.linear_model._linear_model
    :members:
utility
=======

.. currentmodule:: pynumdiff.utils.utility

.. automodule:: pynumdiff.utils.utility
    :members:simulate
========

.. currentmodule:: pynumdiff.utils.simulate

.. automodule:: pynumdiff.utils.simulate
    :members:_pi_cruise_control
==================

.. currentmodule:: pynumdiff.utils._pi_cruise_control

.. automodule:: pynumdiff.utils._pi_cruise_control
    :members:evaluate
========

.. currentmodule:: pynumdiff.utils.evaluate

.. automodule:: pynumdiff.utils.evaluate
    :members:_finite_difference
=====================

.. currentmodule:: pynumdiff.finite_difference._finite_difference

.. automodule:: pynumdiff.finite_difference._finite_difference
    :members:
__optimize__
==============

.. currentmodule:: pynumdiff.optimize.__optimize__

.. automodule:: pynumdiff.optimize.__optimize__
    :special-members: __optimize__
__kalman_smooth__
==================

.. currentmodule:: pynumdiff.optimize.kalman_smooth.__kalman_smooth__

.. automodule:: pynumdiff.optimize.kalman_smooth.__kalman_smooth__
    :members:__total_variation_regularization__
===================================

.. currentmodule:: pynumdiff.optimize.total_variation_regularization.__total_variation_regularization__

.. automodule:: pynumdiff.optimize.total_variation_regularization.__total_variation_regularization__
    :members:__smooth_finite_difference__
============================

.. currentmodule:: pynumdiff.optimize.smooth_finite_difference.__smooth_finite_difference__

.. automodule:: pynumdiff.optimize.smooth_finite_difference.__smooth_finite_difference__
    :members:__finite_difference__
======================

.. currentmodule:: pynumdiff.optimize.finite_difference.__finite_difference__

.. automodule:: pynumdiff.optimize.finite_difference.__finite_difference__
    :members:
__linear_model__
================

.. currentmodule:: pynumdiff.optimize.linear_model.__linear_model__

.. automodule:: pynumdiff.optimize.linear_model.__linear_model__
    :members:_total_variation_regularization
=================================

.. currentmodule:: pynumdiff.total_variation_regularization._total_variation_regularization

.. automodule:: pynumdiff.total_variation_regularization._total_variation_regularization
    :members:_smooth_finite_difference
===========================

.. currentmodule:: pynumdiff.smooth_finite_difference._smooth_finite_difference

.. automodule:: pynumdiff.smooth_finite_difference._smooth_finite_difference
    :members: