# Development notes

## Conda environment

Only required if you wish to install all the required prerequisites in [Miniconda](https://docs.conda.io/en/latest/miniconda.html)/[Anaconda](https://www.anaconda.com/) for local development and testing.

```
conda env create -f environment.yml
```

Then activate with `conda activate synthia`.


## Install synthia

During development:

```
pip install -e .[full]
pip install pytest
```


## Documentation

```
sphinx-build -v -b html docs/ docs/_build/
```

Use `SKIP_NB=1` to skip building the example notebooks (they take a while to build):
```
SKIP_NB=1 sphinx-build -v -b html docs/ docs/_build/
```

## Docstrings

Use [Google Style Python Docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/index.html#google-vs-numpy).

Note that type hints are intentionally not rendered as some of them become too complex and are better described in prose, following xarray's style.

## Testing

```
python -m pytest -s tests/
```

## Versioning

This project uses [semantic versioning](https://semver.org/).

## Deployment

Create and upload a new release with the following commands

```
python setup.py bdist_wheel
pip install --upgrade twine
python -m twine upload dist/*
```
<div align="center">
  <img src="assets/img/logo.png" alt="synthia" height="120"><br><br>

  [![PyPI](https://img.shields.io/pypi/v/synthia)](https://pypi.org/project/synthia) [![CI](https://github.com/dmey/synthia/workflows/CI/badge.svg)](https://github.com/dmey/synthia/actions) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02863/status.svg)](https://doi.org/10.21105/joss.02863)

  [Overview](#overview) | [Documentation](#documentation) | [How to cite](#how-to-cite) | [Contributing](#contributing) | [Development notes](#development-notes) | [Copyright and license](#copyright-and-license) | [Acknowledgements](#acknowledgements)
</div>

## Overview

Synthetic data need to preserve the statistical properties of real data in terms of their individual behavior and (inter-)dependences. [Copula](https://dmey.github.io/synthia/copula.html) and [functional Principle Component Analysis (fPCA)](https://dmey.github.io/synthia/fpca.html) are statistical models that allow these properties to be simulated ([Joe 2014](https://doi.org/10.1201/b17116)). As such, copula generated data have shown potential to improve the generalization of machine learning (ML) emulators ([Meyer et al. 2021](https://doi.org/10.5194/gmd-14-5205-2021)) or anonymize real-data datasets ([Patki et al. 2016](https://doi.org/10.1109/DSAA.2016.49)).

Synthia is an open source Python package to model univariate and multivariate data, parameterize data using empirical and parametric methods, and manipulate marginal distributions. It is designed to enable scientists and practitioners to handle labelled multivariate data typical of computational sciences. For example, given some vertical profiles of atmospheric temperature, we can use Synthia to generate new but statistically similar profiles in just three lines of code (Table 1).

Synthia supports three methods of multivariate data generation through: (i) fPCA, (ii) parametric (Gaussian) copula, and (iii) vine copula models for continuous (all), discrete (vine), and categorical (vine) variables. It has a simple and succinct API to natively handle [xarray](https://xarray.pydata.org)'s labelled arrays and datasets. It uses a pure Python implementation for fPCA and Gaussian copula, and relies on the fast and well tested C++ library [vinecopulib](https://github.com/vinecopulib/vinecopulib) through [pyvinecopulib](https://github.com/vinecopulib/pyvinecopulib)'s bindings for fast and efficient computation of vines. For more information, please see the website at https://dmey.github.io/synthia.


**Table 1**. *Example application of Gaussian and fPCA classes in Synthia. These are used to generate random profiles of atmospheric temperature similar to those included in the source data. The xarray dataset structure is maintained and returned by Synthia.*

| Source                                       | Synthetic with Gaussian Copula                           | Synthetic with fPCA                              |
| -------------------------------------------- | -------------------------------------------------------- | ------------------------------------------------ |
| `ds = syn.util.load_dataset()`               | `g = syn.CopulaDataGenerator()`                          | `g = syn.fPCADataGenerator()`                    |
|                                              | `g.fit(ds, syn.GaussianCopula())`                        | `g.fit(ds)`                                      |
|                                              | `g.generate(n_samples=500)`                              | `g.generate(n_samples=500)`                      |
|                                              |                                                          |                                                  |
| ![Source](./assets/img/temperature_true.png) | ![Gaussian](./assets/img/temperature_synth_gaussian.png) | ![fPCA](./assets/img/temperature_synth_fPCA.png) |


## Documentation

For installation instructions, getting started guides and tutorials, background information, and API reference summaries, please see the [website](https://dmey.github.io/synthia).


## How to cite

If you are using Synthia, please cite the following two papers using their respective Digital Object Identifiers (DOIs). Citations may be generated automatically using Crosscite's [DOI Citation Formatter](https://citation.crosscite.org/) or from the BibTeX entries below.

| Synthia Software                                                | Software Application                                                      |
| --------------------------------------------------------------- | ------------------------------------------------------------------------- |
| DOI: [10.21105/joss.02863](https://doi.org/10.21105/joss.02863) | DOI: [10.5194/gmd-14-5205-2021](https://doi.org/10.5194/gmd-14-5205-2021) |

```bibtex
@article{Meyer_and_Nagler_2021,
  doi = {10.21105/joss.02863},
  url = {https://doi.org/10.21105/joss.02863},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {2863},
  author = {David Meyer and Thomas Nagler},
  title = {Synthia: multidimensional synthetic data generation in Python},
  journal = {Journal of Open Source Software}
}

@article{Meyer_and_Nagler_and_Hogan_2021,
  doi = {10.5194/gmd-14-5205-2021},
  url = {https://doi.org/10.5194/gmd-14-5205-2021},
  year = {2021},
  publisher = {Copernicus {GmbH}},
  volume = {14},
  number = {8},
  pages = {5205--5215},
  author = {David Meyer and Thomas Nagler and Robin J. Hogan},
  title = {Copula-based synthetic data augmentation for machine-learning emulators},
  journal = {Geoscientific Model Development}
}
```

If needed, you may also cite the specific software version with [its corresponding Zendo DOI](https://doi.org/10.5281/zenodo.4701278). 

## Contributing

If you are looking to contribute, please read our [Contributors' guide](CONTRIBUTING.md) for details.


## Development notes

If you would like to know more about specific development guidelines, testing and deployment, please refer to our [development notes](DEVELOP.md).


## Copyright and license

Copyright 2020 D. Meyer and T. Nagler. Licensed under [MIT](LICENSE.txt).


## Acknowledgements

Special thanks to [@letmaik](https://github.com/letmaik) for his suggestions and contributions to the project.
# Contributing

Thank you for considering contributing to Synthia. You can contribute by reporting an issue or by directly contributing to the source code. For the latter, fork our repository, clone the fork, make your changes, and create a pull request (PR) with a clear description of your changes -- if you are unfamiliar about forking/creating PR, please see [this guide](https://guides.github.com/activities/forking/) first. If/when your changes are merged, you will appear as one of our [Contributors](https://github.com/dmey/synthia/graphs/contributors). If you would like to know more about specific development guidelines and testing, please refer to our [development notes](DEVELOP.md).

If you would like to report a bug, please check if similar issue have already been reported [here](https://github.com/dmey/synthia/issues). If none exist please create a new issue and include as many details as possible to reproduce the issue.
---
title: 'Synthia: multidimensional synthetic data generation in Python'
tags:
  - machine-learning
  - data-science
  - python
authors:
  - name: David Meyer
    orcid: 0000-0002-7071-7547
    affiliation: "1, 2"
  - name: Thomas Nagler
    orcid: 0000-0003-1855-0046
    affiliation: 3
affiliations:
 - name: Department of Meteorology, University of Reading, Reading, UK
   index: 1
 - name: Department of Civil and Environmental Engineering, Imperial College London, London, UK
   index: 2
 - name: Mathematical Institute, Leiden University, Leiden, The Netherlands
   index: 3
date: 22 October 2020
bibliography: paper.bib
---


# Summary and Statement of Need

Synthetic data -- artificially generated data that mimic the original (observed) data by preserving relationships between variables [@Nowok_2016] -- may be useful in several areas such as healthcare, finance, data science, and machine learning [@Dahmen_2019; @Kamthe_2021; @Nowok_2016; @Patki_2016]. As such, copula-based data generation models -- probabilistic models that allow for the statistical properties of observed data to be modelled in terms of individual behavior and (inter-)dependencies [@Joe_2014] -- have shown potential in several applications such as finance, data science, and meteorology [@Kamthe_2021; @Li_2020; @Meyer_2021a; @Patki_2016]. Although copula-based data generation tools have been developed for tabular data -- e.g. the Synthetic Data Vault project using Gaussian copulas and generative adversarial networks [@Patki_2016; @Xu_2018], or the Synthetic Data Generation via Gaussian Copula [@Li_2020] -- in computational sciences such as weather and climate, data often consist of large, labelled multidimensional datasets with complex dependencies.

Here we introduce Synthia, an open-source multidimensional synthetic data generator written in Python for xarray's [@Hoyer_2017] labelled arrays and datasets with support for parametric and vine copulas models and functional principal component analysis (fPCA) -- an extension of principal component analysis where data consist of functions instead of vectors [@Ramsay_2005] -- to allow for a wide range of data and dependent structures to be modelled. For efficiency, algorithms are implemented in NumPy [@Harris_2020] and SciPy [@SciPy_2020] for Gaussian (parametric) copula and fPCA classes and rely on the C++ library vinecopulib [@Nagler_2020a] through pyvinecopulib's [@Nagler_2020b] bindings for fast computation of vines.

Recent applications of Synthia include the generation of dependent [@Meyer_2021a] and independent [@Meyer_2021b] samples for improving the prediction of machine learning emulators in weather and climate. In this release we include examples and tutorials for univariate and multivariate synthetic data generation using copula and fPCA methods and look forward to enabling the generation of synthetic data in various scientific communities and for several applications.


# Acknowledgments

We thank Maik Riechert for his comments and contributions to the project.


# References
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
Steps to reproduce the behavior.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**System information (please complete the following information):**
 - OS name and version: [e.g. Ubuntu 18, macOS 10.14, Windows 10]
 - Python version
 - Synthia version
 - Dependencies version

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
# Copulas


## What copulas are

[Copulas](https://en.wikipedia.org/wiki/Copula_(probability_theory)) are probabilistic models for the dependence between random quantities $Z_1, \dots, Z_d$. The idea is to 
decouple the individual (or marginal) behavior of the quantities from the dependence. 
The key result, Sklar’s theorem  {cite}`sklar1959`, states that for any joint distribution $F$ on $d$ variables can be written as 

$$F(z_1,\dots, z_d) = C(F_1(z_1), \dots,F_d(z_d)), $$

where $F_1, \dots, F_d$ are the marginal distribution functions and $C$ is the *copula* function. The copula encodes all the information that is not contained in the marginals, i.e., the dependence between variables. 
The simplest example for a copula function is $C(u_1, \dots, u_d) = u_1 \times \cdots \times u_d$ which corresponds to independence: $F(z_1,\dots, z_d) = F_1(z_1) \times \cdots \times F_d(z_d)$. 


## The Gaussian copula

To induce dependence in $F$, it is common to consider sub-families of copulas that are conveniently parametrized. There's a variety of such parametric copula families. The most prominent one is the *Gaussian copula*. It is defined by inversion of Sklar's theorem:

$$C(u_1,\dots, u_d) = F(F_1^{-1}(u_1), \dots, F_d^{-1}(u_d)), $$

where $F$ is a multivariate Gaussian distribution and $F_1, \dots, F_d$ the corresponding marginals. The Gaussian copula is then parameterized by a correlation matrix and subsumes all possible dependence structure in a multivariate Gaussian distribution. The benefit comes from the fact that we can combine a given copula with any type of marginal distributions, not just the ones the copula was derived from. That way, we can build flexible models with arbitrary marginal distributions and Gaussian-like dependence. 


## Other copula families

The same principle applies to other multivariate distributions and many copula models have been derived, most prominently the Student t copula and Archimedean families. A comprehensive list can be found in {cite}`joe2014dependence`. When there are more than two variables ($d>2$) the types of dependence structures these models can generate is rather limited. For example, Gaussian and Student copulas only allow for symmetric dependencies. While Archimedean families allow for such asymmetries, they require all pairs of variables to have the same type and strength of dependence.

## Vine copulas

[Vine copula models](https://en.wikipedia.org/wiki/Vine_copula) ({cite}`aas2009pair`, {cite}`czado2019analyzing`) are a popular solution to this issue. The idea is to build a large dependence model from only two-dimensional building blocks. We can explain this with a simple example with just three variables $Z_1, Z_2, Z_3$.

We model the dependence between $Z_1$ and $Z_2$ by a two-dimensional copula 
$C_{1,2}$ and the dependence between $Z_2$ and $Z_3$ by another, possibly different, copula $C_{2,3}$. These two copulas already contain some information about the dependence between $Z_1$ and $Z_3$, the part of the dependence that is induced by $Z_2$. The missing piece is the dependence between $Z_1$ and $Z_3$ after the effect of $Z_2$ has been removed. Mathematically, this is the conditional dependence $Z_1$ and $Z_3$ given 
$Z_2$ and can be modeled by yet another two-dimensional copula 
$C_{1,3|2}$. The principle is easily extended to an arbitrary number of variables $Z_1, \dots, Z_d$. 

Because all two-dimensional copulas can be specified independently, such models are extremely flexible and allow for highly heterogeneous dependence structures. Algorithms for simulation and selecting the right conditioning order and parametric families for each (conditional) pair are given in {cite}`dissmann2013`. Using parametric models for pair-wise dependencies remain a limiting factor, however. If necessary, it is also possible to use nonparametric models for the two-dimensional building blocks, see {cite}`Nagler2017`.


```{bibliography} references.bib
:style: alpha
```
# Installation

## Required dependencies

- Python (3.8 or later)
- [numpy](http://www.numpy.org/)
- [scipy](https://www.scipy.org/scipylib/index.html)
- [xarray](http://xarray.pydata.org/)

## Optional dependencies

- [pyvinecopulib](https://github.com/vinecopulib/pyvinecopulib): For vine copula and quasirandom number support.

## Instructions

```
pip install synthia
```

or with optional [pyvinecopulib](https://github.com/vinecopulib/pyvinecopulib):

```
pip install synthia[full]
```
# How to cite

If you are using Synthia, please cite the following two papers using their respective Digital Object Identifiers (DOIs). Citations may be generated automatically using Crosscite's [DOI Citation Formatter](https://citation.crosscite.org/) or from the BibTeX entries below. If needed, you may also cite the specific software version with [its corresponding Zendo DOI](https://doi.org/10.5281/zenodo.4701278). 

| Synthia Software                   | Software Application                                              |
| ---------------------------------- | ----------------------------------------------------------------- |
| DOI: 10.21105/joss.02863 | DOI: [10.5194/gmd-2020-427](https://doi.org/10.5194/gmd-2020-427) |

```bibtex
@article{Meyer_Nagler_2021,
  title   = {Synthia: multidimensional synthetic data generation in Python},
  author  = {David Meyer and Thomas Nagler},
  year    = {2021},
  doi     = {10.21105/joss.02863},
  journal = {Journal of Open Source Software},
  note    = {Under review}
}

@article{Meyer_Nagler_Hogan_2021,
  title   = {Copula-Based Synthetic Data Generation for Machine Learning Emulators in Weather and Climate: Application to a Simple Radiation Model},
  author  = {David Meyer and Thomas Nagler and Robin J. Hogan},
  year    = {2021},
  volume  = {2021},
  doi     = {10.5194/gmd-2020-427},
  journal = {Geoscientific Model Development Discussions},
  note    = {Under review}
}
```

# Features

Synthia includes several software features for the generation and augmentation of multivariate data (Table 1). You can find out more about how to use these option in the [examples page](examples.rst).

**Table 1**. *Software features included in the latest version of Synthia for different fitting and generation methods. For correlated multivariate variables stretching and uniformization is available for copula models (^). All models support the handling of continuous variables (¹). Discrete (²) and categorical (³) variables are supported with vine copulas. All correlated multivariate methods support empirical and parametric fitting of marginal distributions*

| Variable relation                      | Method     | Option             |
| -------------------------------------- | ---------- | ------------------ |
| Univariate or Multivariate Independent |            |                    |
|                                        | Fitting    |                    |
|                                        |            | Empirical          |
|                                        |            | Parameterized      |
|                                        | Generation |                    |
|                                        |            | Pseudo-random      |
|                                        |            | Quasi-random       |
|                                        |            | Stretching         |
|                                        |            | Uniformization     |
| Multivariate                           |            |                    |
|                                        | Fitting    |                    |
|                                        |            | ^¹Copula: Gaussian |
|                                        |            | ^¹²³Copula: Vine   |
|                                        |            | ¹fPCA              |
|                                        |            | Empirical          |
|                                        |            | Parameterized      |
|                                        | Generation |                    |
|                                        |            | Pseudo-random      |
|                                        |            | Quasi-random       |
|                                        |            | ^Stretching        |
|                                        |            | ^Uniformization    |
# Functional Principal Component Analysis (fPCA)

## The general idea
[Principal component analysis (PCA)](https://en.wikipedia.org/wiki/Principal_component_analysis) is popular technique for dimensionality reduction and data compression. The idea is to find the directions in a high-dimensional space in which we find the most variation in the data. To reduce dimensionality (or compress the data), we then only look at projections on the first few principal directions. 

[Functional principal component analysis (fPCA)](https://en.wikipedia.org/wiki/Functional_principal_component_analysis) is an extension of this method to situations where data consist of functions, not vectors. Although mathematically more involved than regular PCA, the two are equivalent from a practical point of view. To keep it simple, we shall therefore explain the core ideas in terms of the regular PCA.

## Mathematical definition

Suppose each observation in the data is a vector $X$. The first principal component is defined as the direction $w$ in which there's the most variation:

$$v_1 = \arg\max_{\|v\| = 1} \mathrm{var}[X'v].$$

In practice, the variance of $X'v$ is approximated by the [sample variance](https://en.wikipedia.org/wiki/Variance#Sample_variance). The other principal components are defined similarly,

$$v_k = \arg\max_{\|v\| = 1} \mathrm{var}[X'v],$$

under the side-constraint that $v_{k}$ is orthogonal to all previous principal components $v_{1}, \dots, v_{k - 1}$. In practice, the collection of all principal components can be found in one go as the eigenvectors of the covariance matrix of $X$.

## PCA as a basis expansion

For a $d$-dimensional vector $X$, there are in total $d$ principal components to be found: $v^{1}, \dots, v_{d}$. The principal components form a basis of the space, i.e., every vector $X$ can be represented as 

$$X = \mu + \sum_{k = 1}^d a_k v_{k},$$

where $\mu$ is $E[X]$ and the coefficients $a_k = X'v_k$ are called *principal component scores*. When $d$ is very large, we can truncate the sum above at $K \ll d$ terms. This gives an approximation 

$$X \approx \mu + \sum_{k = 1}^K a_k v_{k}.$$

This allows us to represent the high-dimensional vector $X$  by only a few numbers $a_1, \dots, a_K$. The number $K$ determines the quality of the approximation.

## PCA for synthetic data generation

PCA can be used to generate synthetic data for the high-dimensional vector $X$. For every instance $X_i$ in the data set, we compute the principal component scores $a_{i, 1}, \dots, a_{i, K}$. The scores are uncorrelated by construction and we treat them as independent. We then fit a model $F_k$ for the marginal distribution of each score $a_k$. From these models, we can generate synthetic scores $\tilde a_k$ and transform them into a synthetic sample of $X$ via 

$$\tilde X = \mu + \sum_{k = 1}^K \tilde a_k v_{k}.$$
Examples
========


Generation
**********
.. toctree::
    :maxdepth: 1

    examples/univariate
    examples/multivariate-independent
    examples/multivariate-gaussian
    examples/multivariate-vine
    examples/fpca
    examples/discrete-and-categorical

Enhancement
***********
.. toctree::
    :maxdepth: 1

    examples/enhancement
Getting Started
---------------

* :doc:`overview`
* :doc:`features`
* :doc:`installation`
* :doc:`quickstart`
* :doc:`examples`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   overview
   features
   installation
   quickstart.ipynb
   examples

Background
---------------------

* :doc:`copula`
* :doc:`fpca`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Background

   copula
   fpca


Help & Reference
----------------

* :doc:`api`
* :doc:`citing`
* :doc:`contributing`
* :doc:`develop`
* :doc:`license`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Help & reference

   api
   citing
   contributing
   develop
   license.. currentmodule:: synthia

#############
API reference
#############

This page provides an auto-generated summary of synthia's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.


Data generators
===============

.. autoclass:: synthia.CopulaDataGenerator
   :members:
   :undoc-members:

.. autoclass:: synthia.FPCADataGenerator
   :members:
   :undoc-members:


Copulas
=======

.. autoclass:: synthia.GaussianCopula
   :members:
   :undoc-members:

.. autoclass:: synthia.VineCopula
   :members:
   :undoc-members:

Parameterizers
==============

.. autoclass:: synthia.ConstParameterizer
   :members:
   :undoc-members:

.. autoclass:: synthia.QuantileParameterizer
   :members:
   :undoc-members:

.. autoclass:: synthia.DistributionParameterizer
   :members:
   :undoc-members:

Transformers
============

.. autoclass:: synthia.ArcTanhTransformer
   :members:
   :undoc-members:

.. autoclass:: synthia.BoxCoxTransformer
   :members:
   :undoc-members:

.. autoclass:: synthia.CombinedTransformer
   :members:
   :undoc-members:

Utilities
=========

.. automodule:: synthia.util
   :members: load_dataset
   :undoc-members: