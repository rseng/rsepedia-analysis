# equadratures

*equadratures* is an open-source library for *uncertainty quantification*, *machine learning*, *optimisation*, *numerical integration* and *dimension reduction* -- all using orthogonal polynomials. It is particularly useful for models / problems where output quantities of interest are smooth and continuous; to this extent it has found widespread applications in computational engineering models (finite elements, computational fluid dynamics, etc). It is built on the latest research within these areas and has both deterministic and randomized algorithms. Effective Quadratures is actively being developed by researchers at the [University of Cambridge](https://www.cam.ac.uk), [Imperial College London](https://www.imperial.ac.uk), [Stanford University](https://www.stanford.edu), [The University of Utah](https://www.utah.edu), [The Alan Turing Institute](https://www.turing.ac.uk) and the [University of Cagliari](https://www.unica.it/unica/). *equadratures* is a NumFOCUS affiliated project.

**Key words associated with this code**: polynomial surrogates, polynomial chaos, polynomial variable projection, Gaussian quadrature, Clenshaw Curtis, polynomial least squares, compressed sensing, gradient-enhanced surrogates, supervised learning.

## Code

The latest version of the code is v9.1.0 *Narwhal*, released June 2021. 

![](https://travis-ci.com/equadratures/equadratures.svg?branch=master)
[![](https://coveralls.io/repos/github/Effective-Quadratures/Effective-Quadratures/badge.svg?branch=master)](https://coveralls.io/github/Effective-Quadratures/Effective-Quadratures)
[![](https://badge.fury.io/py/equadratures.svg)](https://pypi.org/project/equadratures/)
[![](https://joss.theoj.org/papers/10.21105/joss.00166/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00166)
[![](https://img.shields.io/pypi/pyversions/equadratures.svg)](https://pypi.python.org/pypi/equadratures)
![](https://img.shields.io/github/stars/Effective-Quadratures/Effective-Quadratures.svg?style=flat-square&logo=github&label=Stars&logoColor=white)
![](https://static.pepy.tech/badge/equadratures/week)
[![](https://img.shields.io/discourse/status?server=https%3A%2F%2Fdiscourse.equadratures.org)](https://discourse.equadratures.org)

If you use `pip` you can install the code with:

```python
pip install equadratures
```

or `pip` can be replaced with `python -m pip`, where `python` is the python version you wish to install *equadratures* for. Use of a virtual enviroment such as [virtualenv](https://pypi.org/project/virtualenv/) or [pyenv](https://github.com/pyenv/pyenv)/[pipenv](https://pypi.org/project/pipenv/) is also encouraged. Alternatively you can click either on the **Fork Code** button or **Clone**, and install from your local version of the code.

For issues with the code, please do *raise an issue* on our Github page; do make sure to add the relevant bits of code and specifics on package version numbers. We welcome contributions and suggestions from both users and folks interested in developing the code further.

Our code is designed to require minimal dependencies; current package requirements include ``numpy``, ``scipy`` and ``matplotlib``.

## Documentation, tutorials, Discourse

Code documentation and details on the syntax can be found [here](https://equadratures.org/index.html).

We've recently started a Discourse forum! Check it out [here](https://discourse.equadratures.org/).

## Code objectives

Specific goals of this code include:

* probability distributions and orthogonal polynomials
* supervised machine learning: regression and compressive sensing
* numerical quadrature and high-dimensional sampling
* transforms for correlated parameters
* computing moments from models and data-sets
* sensitivity analysis and Sobol' indices
* data-driven dimension reduction
* ridge approximations and neural networks
* surrogate-based design optimisation 

## Get in touch

Feel free to follow us via [Twitter](https://twitter.com/EQuadratures) or email us at contact@effective-quadratures.org. 


## Community guidelines

If you have contributions, questions, or feedback use either the Github repository, or get in touch. We welcome contributions to our code. In this respect, we follow the [NumFOCUS code of conduct](https://numfocus.org/code-of-conduct). 

## Acknowledgments

This work was supported by wave 1 of The UKRI Strategic Priorities Fund under the EPSRC grant EP/T001569/1, particularly the [Digital Twins in Aeronautics](https://www.turing.ac.uk/research/research-projects/digital-twins-aeronautics) theme within that grant, and [The Alan Turing Institute](https://www.turing.ac.uk).
---
title: 'Effective-Quadratures (EQ): Polynomials for Computational Engineering Studies'
tags:
  - effective quadrature subsampling
  - tensor and sparse grids
  - uncertainty quantification
  - surrogate modelling
  - polynomial interpolation and approximation
authors:
 - name: Pranay Seshadri
   orcid: 0000-0002-7351-012X
   affiliation: 1
 - name: Geoffrey Parks
   orcid: 0000-0001-8188-5047
   affiliation: 1
affiliations:
 - name: University of Cambridge
   index: 1
date: 23 November 2016
bibliography: paper.bib
---

# Summary

Effective-Quadratures (EQ) is a suite of tools for generating polynomials for parametric computational studies. These polynomials may be used to compute approximations to scientific models; tailored specifically for uncertainty quantification, approximation and numerical integration studies. This is a research code designed for the development and testing of new polynomial techniques in the areas listed above. EQ is written entirely in python and requires Numpy, Scipy and Matplotlib.

For a computational engineering problem, albeit analytical or a black-box model, EQ constructs a polynomial approximation (or interpolant) by sampling the model at specific points in its parameter space. This may be done using tensor or sparse grids—routines for which are available in EQ—based on the work of [@Xiu2002, @Constantine2012]. The code also has routines for effectively subsampling an existing tensor grid for computing least squares approximations, based on the work of [@Seshadri2016], which uses a QR column pivoting heuristic. Once these polynomials have been generated, they may be sampled, plotted, or used for computing moments. 

The need for this software is two fold: (1) to replicate state-of-the-art results in polynomial-based uncertainty quantification (for which the number of open source codes are limited), and (2) to provide an easy-to-use platform for carrying out further research in multivariate polynomial generation. 

# References
[@Xiu2002] Xiu, Dongbin, and George Em Karniadakis. "The Wiener--Askey polynomial chaos for stochastic differential equations." SIAM journal on scientific computing 24.2 (2002): 619-644.

[@Constantine2012] Constantine, Paul G., Michael S. Eldred, and Eric T. Phipps. "Sparse pseudospectral approximation method." Computer Methods in Applied Mechanics and Engineering 229 (2012): 1-12.

[@Seshadri2016] Seshadri, Pranay, Akil Narayan, and Sankaran Mahadevan. "Effectively Subsampled Quadratures For Least Squares Polynomial Approximations." arXiv preprint arXiv:1601.05470 (2016).
# Effective Quadratures
**Version 8.0**

## Instructions for adding your own sampling method.
Insert your python file in the ``sampling_methods/`` folder. In the folder above, open ``quadrature.py`` and insert the relevant header and add your method to the ``Quadrature`` constructor. To compute quadrature points, your method should only require a list of ``parameters`` and the ``basis``. 
