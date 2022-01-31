# GGLasso

[![PyPI version fury.io](https://badge.fury.io/py/gglasso.svg)](https://pypi.python.org/pypi/gglasso/)
[![PyPI license](https://img.shields.io/pypi/l/gglasso.svg)](https://pypi.python.org/pypi/gglasso/)
[![Python version](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8%20%7C%203.9-blue)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/gglasso/badge/?version=latest)](http://gglasso.readthedocs.io/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03865/status.svg)](https://doi.org/10.21105/joss.03865)
[![arXiv](https://img.shields.io/badge/arXiv-2011.00898-b31b1b.svg)](https://arxiv.org/abs/2110.10521)


This package contains algorithms for solving General Graphical Lasso (GGLasso) problems, including single, multiple, as well as latent 
Graphical Lasso problems. <br>

[Docs](https://gglasso.readthedocs.io/en/latest/) | [Examples](https://gglasso.readthedocs.io/en/latest/auto_examples/index.html)

## Getting started

### Install via pip

The package is available on pip and can be installed with

    pip install gglasso

### Install from source

Alternatively, you can install the package from source using the following commands:

    git clone https://github.com/fabian-sp/GGLasso.git
    pip install -r requirements.txt
    python setup.py

Test your installation with 

    pytest gglasso/ -v


### Advanced options

When installing from source, you can also install dependencies with `conda` via the command

	$ while read requirement; do conda install --yes $requirement || pip install $requirement; done < requirements.txt

If you wish to install `gglasso` in developer mode, i.e. not having to reinstall `gglasso` everytime the source code changes (either by remote or local changes), run

    python setup.py clean --all develop clean --all

## The `glasso_problem` class

`GGLasso` can solve multiple problem forumulations, e.g. single and multiple Graphical Lasso problems as well as with and without latent factors. Therefore, the main entry point for the user is the `glasso_problem` class which chooses automatically the correct solver and model selection functionality. See [our documentation](https://gglasso.readthedocs.io/en/latest/problem-object.html) for all the details.


## Algorithms

`GGLasso` contains algorithms for Single and Multiple Graphical Lasso problems. Moreover, it allows to model latent variables (Latent variable Graphical Lasso) in order to estimate a precision matrix of type **sparse - low rank**. The following algorithms are contained in the package.
<br>
1) ADMM for Single Graphical Lasso<br>

2) ADMM for Group and Fused Graphical Lasso<br>
The algorithm was proposed in [2] and [3]. To use this, import `ADMM_MGL` from `gglasso/solver/admm_solver`.<br>

3) A Proximal Point method for Group and Fused Graphical Lasso<br>
We implement the PPDNA Algorithm like proposed in [4]. To use this, import `warmPPDNA` from `gglasso/solver/ppdna_solver`.<br>

4) ADMM method for Group Graphical Lasso where the features/variables are non-conforming<br>
Method for problems where not all variables exist in all instances/datasets.  To use this, import `ext_ADMM_MGL` from `gglasso/solver/ext_admm_solver`.<br>

## Citation

If you use `GGLasso`, please consider the following citation

    @article{Schaipp2021,
      doi = {10.21105/joss.03865},
      url = {https://doi.org/10.21105/joss.03865},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {68},
      pages = {3865},
      author = {Fabian Schaipp and Oleg Vlasovets and Christian L. Müller},
      title = {GGLasso - a Python package for General Graphical Lasso computation},
      journal = {Journal of Open Source Software}
    }


## Community Guidelines

1)  Contributions and suggestions to the software are always welcome.
    Please, consult our [contribution guidelines](CONTRIBUTING.md) prior
    to submitting a pull request.
2)  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/fabian-sp/GGLasso/issues).
3)  Contributors must adhere to the [Code of
    Conduct](CODE_OF_CONDUCT.md).


## References
*  [1] Friedman, J., Hastie, T., and Tibshirani, R. (2007).  Sparse inverse covariance estimation with the Graphical Lasso. Biostatistics, 9(3):432–441.
*  [2] Danaher, P., Wang, P., and Witten, D. M. (2013). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2):373–397.
* [3] Tomasi, F., Tozzo, V., Salzo, S., and Verri, A. (2018). Latent Variable Time-varying Network Inference. InProceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining. ACM.
* [4] Zhang, Y., Zhang, N., Sun, D., and Toh, K.-C. (2020). A proximal point dual Newton algorithm for solving group graphical Lasso problems. SIAM J. Optim., 30(3):2197–2220.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at otorrent@mail.ru.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contributing to `GGLasso`

## Bug reports and feature requests

In case you have questions, possible bugs or feature requests, please, open an [issue](https://github.com/fabian-sp/GGLasso/issues) on Github. 

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the need for change. Please, describe the problem your PR addresses and show how your pull request solves it as concisely as possible.

Also, include this motivation in `NEWS` so that when a new release of
`gglasso` comes out it is easy for users to see what has changed. Add your
item at the top of the file and use markdown for formatting. The
news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please, check to make sure that you have NOT accidentally included any unrelated changes. These make it harder to see exactly what was changed, and to evaluate any unexpected side effects.

Each PR corresponds to a git branch, so if you expect to submit multiple changes
make sure to create multiple branches. If you have multiple changes that depend
on each other, start with the first one and do not submit any others until the
first one has been processed.

3. If you add new parameters or a new function, you will also need to document them.

## Contributing code

The procedure for contributing code is the following:

1. Fork the `GGLasso` repository into your own Github account (recommended: create a new branch for your feature).
2. Install the forked repository **from source** in developer mode (see the sections *Install from source* and *Advanced options* in our `README`).
3. Add your changes and commit them. 
4. Run `pytest gglasso/tests` to see whether all tests are successfull. If you add new functionalities, add unit tests for these.
5. Open a pull request for your feature branch on Github.

Also, have a look at [this guide on forking in Github](https://guides.github.com/activities/forking/).

Useful information on how to write bug reports and the principles of contributing to open-source code projects are described on [the scikit-learn webpage](https://scikit-learn.org/stable/developers/contributing.html).


## Fixing typos
You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the source file.

---
title: 'GGLasso - a Python package for General Graphical Lasso computation'
tags:
  - Python
  - graphical lasso
  - latent graphical model
  - structured sparsity
  - convex optimization
  - ADMM
authors:
  - name: Fabian Schaipp
    orcid: 0000-0002-0673-9944
    affiliation: "1"
  - name: Oleg Vlasovets
    affiliation: "2,3"
  - name: Christian L. Müller
    orcid: 0000-0002-3821-7083
    affiliation: "2,3,4"

affiliations:
  - name: Technische Universität München
    index: 1
  - name: Institute of Computational Biology, Helmholtz Zentrum München
    index: 2
  - name: Department of Statistics, Ludwig-Maximilians-Universität München
    index: 3
  - name: Center for Computational Mathematics, Flatiron Institute, New York
    index: 4
date: 17 November 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:

---

# Summary

We introduce `GGLasso`, a Python package for solving General Graphical Lasso problems. The Graphical Lasso scheme, introduced by @Friedman2007 (see also @Yuan2007, @Banerjee2008), estimates a sparse inverse covariance matrix $\Theta$ from multivariate Gaussian data $\mathcal{X} \sim \mathcal{N}(\mu, \Sigma) \in \mathbb{R}^p$. Originally proposed by @Dempster1972 under the name Covariance Selection, this estimation framework has been extended to include latent variables in @Chandrasekaran2012. Recent extensions also include the joint estimation of multiple inverse covariance matrices, see, e.g., in @Danaher2013, @Tomasi2018. The `GGLasso` package contains methods for solving the general problem formulation:

<div class="math">
\begin{align}
\label{eq:problem}
\min_{\Theta, L \in \mathbb{S}_{++}^K }\quad \sum_{k=1}^{K} \left(-\log\det(\Theta^{(k)} - L^{(k)}) + \langle S^{(k)},  \Theta^{(k)} - L^{(k)} \rangle \right)+ \mathcal{P}(\Theta) +\sum_{k=1}^{K} \mu_{1,k} \|L^{(k)}\|_{\star}.
\end{align}
</div>

Here, we denote with $\mathbb{S}_{++}^K$ the $K$-product of the space of symmetric, positive definite matrices. Moreover, we write $\Theta = (\Theta^{(1)},\dots,\Theta^{(K)})$ for the sparse component of the inverse covariances and $L = (L^{(1)},\dots,L^{(K)})$ for the low rank components, formed by potential latent variables. Here, $\mathcal{P}$ is a regularization function that induces a desired sparsity structure. The above problem formulation subsumes important special cases, including the single (latent variable) Graphical Lasso, the Group, and the Fused Graphical Lasso.

# Statement of need

Currently, there is no Python package available for solving general Graphical Lasso instances. The standard single Graphical Lasso problem (SGL) can be solved in `scikit-learn` [@Pedregosa2011]. The `skggm` package provides several algorithmic and model selection extensions for the single Graphical Lasso problem [@Laska2017]. The package `regain` [@Tomasi2018] comprises solvers for single and Fused Graphical Lasso problems, with and without latent variables. With `GGLasso`, we make the following contributions:

- Proposing a uniform framework for solving Graphical Lasso problems.
- Providing solvers for Group Graphical Lasso problems (with and without latent variables).
- Providing a solver for -- what we call -- *nonconforming GGL* problems where not all variables need to be present in every instance. We detail a use case of this novel extension on synthetic data.
- Implementing a block-wise ADMM solver for SGL problems following @Witten2011 as well as proximal point solvers for FGL and GGL problems [@Zhang2021; @Zhang2020].

In the table below we give an overview of existing functionalities and the `GGLasso` package.


|                                 | scikit-learn |  regain     |  GGLasso | comment |
| -----------                     | -----------  | ----------- | ----------- | ----------- |
| SGL                             | **yes**      | **yes**       | **yes**       | new: block-wise solver           |
| SGL + latent                    | **no**       | **yes**       | **yes**       |             |
| GGL                             | **no**       | **no**          | **yes**       |             |
| GGL + latent                    | **no**       | **no**          | **yes**       |             |
| FGL                             | **no**       | **yes**       | **yes**       | new: proximal point solver            |
| FGL + latent                    | **no**       | **yes**       | **yes**       |             |
| GGL nonconforming  (+latent)    | **no**       | **no**       | **yes**       |             |



# Functionalities

## Installation and problem instantiation

`GGLasso` can be installed via `pip`.

```shell
pip install gglasso
```

The central object of `GGLasso` is the class `glasso_problem`, which streamlines the solving or model selection procedure for SGL, GGL, and FGL problems with or without latent variables.

As an example, we instantiate a single Graphical Lasso problem (see the problem formulation below). We input the empirical covariance matrix `S` and the number of samples `N`. We can choose to model latent variables and set the regularization parameters via the other input arguments.

```python
# Import the main class of the package
from gglasso.problem import glasso_problem

# Define a SGL problem instance with given data S
problem  = glasso_problem(S, N, reg = None,
                          reg_params = {'lambda1': 0.01}, latent = False)
```

As a second example, we instantiate a Group Graphical Lasso problem with latent variables. Typically, the optimal choice of the regularization parameters are not known and are determined via model selection.

```python
# Define a GGL problem instance with given data S
problem  = glasso_problem(S, N, reg = "GGL", reg_params = None, latent = True)
```

Depending on the input arguments, `glasso_problem` comprises two main modes:

- if regularization parameters are specified, the problem-dependent default solver is called.
- if regularization parameters are *not* specified, `GGLasso` performs model selection via grid search and the extended BIC criterion [@Foygel2010]).

```python
problem.solve()
problem.model_selection()
```

For further information on the input arguments and methods, we refer to the [detailled documentation](https://gglasso.readthedocs.io/en/latest/problem-object.html).

![Illustration of the latent SGL: The estimated inverse covariance matrix $\hat \Omega$ decomposes into a sparse component $\hat \Theta$ (central) and a low-rank component $\hat L$ (right). \label{fig1}](../docs/source/pictures/SLRDecomp.pdf){width=90%}

## Problem formulation

We list important special cases of the problem formulation given in \autoref{eq:problem}. For a mathematical formulation of each special case, we refer to the [documentation](https://gglasso.readthedocs.io/en/latest/math-description.html).

### Single Graphical Lasso (*SGL*): {#SGL}
For $K=1$, the problem reduces to the single (latent variable) Graphical Lasso where
$$
\mathcal{P}(\Theta) = \lambda_1 \sum_{i \neq j} |\Theta_{ij}|.
$$
An illustration of the single latent variable Graphical Lasso model output is shown in \autoref{fig1}.

### Group Graphical Lasso (*GGL*): {#GGL}
For
$$
\mathcal{P}(\Theta) = \lambda_1 \sum_{k=1}^{K} \sum_{i \neq j} |\Theta_{ij}^{(k)}| + \lambda_2  \sum_{i \neq j} \left(\sum_{k=1}^{K} |\Theta_{ij}^{(k)}|^2 \right)^{\frac{1}{2}}
$$
we obtain the Group Graphical Lasso as formulated in @Danaher2013.

### Fused Graphical Lasso (*FGL*): {#FGL}
For
$$
\mathcal{P}(\Theta) = \lambda_1 \sum_{k=1}^{K} \sum_{i \neq j} |\Theta_{ij}^{(k)}| + \lambda_2  \sum_{k=2}^{K}   \sum_{i \neq j} |\Theta_{ij}^{(k)} - \Theta_{ij}^{(k-1)}|
$$
we obtain Fused (also called Time-Varying) Graphical Lasso [@Danaher2013; @Tomasi2018; @Hallac2017].

### Nonconforming GGL:

Consider the GGL case in a situation where not all variables are observed in every instance $k=1,\dots,K$. `GGLasso` is able to solve these problems and include latent variables. We provide the mathematical details in the [documentation](https://gglasso.readthedocs.io/en/latest/math-description.html#ggl-the-nonconforming-case) and give an [example](https://gglasso.readthedocs.io/en/latest/auto_examples/plot_nonconforming_ggl.html#sphx-glr-auto-examples-plot-nonconforming-ggl-py).


## Optimization algorithms

The `GGLasso` package implements several methods with provable convergence guarantees for solving the optimization problems formulated above.

- *ADMM*: for all problem formulations we implemented the ADMM algorithm [@Boyd2011]. ADMM is a flexible and efficient optimization scheme which is specifically suited for Graphical Lasso problems as it only relies on efficient computation of the proximal operators of the involved functions [@Danaher2013; @Tomasi2018; @Ma2013].  

- *PPDNA*: for GGL and FGL problems without latent variables, we included the proximal point solver proposed in @Zhang2021 and @Zhang2020. According to the numerical experiments in @Zhang2020, PPDNA can be an efficient alternative to ADMM especially for fast local convergence.

- *block-ADMM*: for SGL problems without latent variables, we implemented a method which solves the problem block-wise, following the proposal in @Witten2011. This wrapper simply applies the ADMM solver to all connected components of the empirical covariance matrix after thresholding.

![Runtime comparison for SGL problems of varying dimension and sample size at three different $\lambda_1$ values. The left column shows the runtime at low accuracy, the right column at high accuracy. \label{fig2}](../docs/source/pictures/runtime_merged.pdf){width=90%}


## Benchmarks and applications

In our example gallery, we included benchmarks comparing the solvers in `GGLasso` to state-of-the-art software as well as illustrative examples explaining the usage and functionalities of the package. We want to emphasize the following examples:

- [Benchmarks](https://gglasso.readthedocs.io/en/latest/auto_examples/plot_benchmarks.html#sphx-glr-auto-examples-plot-benchmarks-py) for SGL problems: our solver is competitive with `scikit-learn` and `regain`. The newly implemented block-wise solver is highly efficient for large sparse networks (see \autoref{fig2} for runtime comparison at [low and high accuracy](https://gglasso.readthedocs.io/en/latest/auto_examples/plot_benchmarks.html#calculating-the-accuracy), respectively).

- [Soil microbiome application](https://gglasso.readthedocs.io/en/latest/auto_examples/plot_soil_example.html#sphx-glr-auto-examples-plot-soil-example-py): following @Kurtz2019, we demonstrate how latent variables can be used to identify hidden confounders in microbial network inference.

- [Nonconforming GGL](https://gglasso.readthedocs.io/en/latest/auto_examples/plot_nonconforming_ggl.html#sphx-glr-auto-examples-plot-nonconforming-ggl-py): we illustrate how to use `GGLasso` for GGL problems with missing variables.


# Acknowledgements

We thank Prof. Dr. Michael Ulbrich, TU Munich, for supervising the Master's thesis of FS that led to the development of the software. We also thank Dr. Zachary D. Kurtz for helping with testing of the latent graphical model implementation.

# References
Getting started
======================

.. _Github: https://github.com/fabian-sp/GGLasso

Installation
^^^^^^^^^^^^^^^^

``GGLasso`` is available over Pypi or `Github`_. For installation with pip, simply run 

.. code-block::

     pip install gglasso

To install from source, clone the repository and make sure you have all requirements installed. Then move to the directory and run

.. code-block::

     python setup.py

This installs a package called ``gglasso`` in your Python environment. In case you want to edit the source code and use the ``gglasso`` package without re-installing, you can run instead

.. code-block::

     python setup.py clean --all develop clean --all

To make sure that everything works properly you can run unit tests in ``gglasso/tests``, for example

.. code-block::

     pytest gglasso/ -v

To import from ``GGLasso`` in Python, type for example

.. code-block:: python

     from gglasso.problem import glasso_problem 


Mathematical description
=============================

The ``GGLasso`` package can solve several problem formulations related to Graphical Lasso. On this page, we aim to define the exact formulation for each problem.

Single Graphical Lasso problems (SGL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider a multivariate Gaussian variable

.. math::
   \mathcal{X} \sim \mathcal{N}(\mu, \Sigma) \in \mathbb{R}^p

Zeros in the precision matrix, :math:`\Sigma^{-1}` correspond to conditional independence of two components of :math:`X` which motivates a sparse estimation of :math:`\Sigma^{-1}`.
Typically, we are given a sample of :math:`N` independent samples of :math:`X` for which we can compute the empirical covariance matrix :math:`S`.
Even though :math:`S^{-1}` (if exists) is the maximum-likelihood estimator of the precision matrix, it is not guaranteed to be sparse. 

This leads to the nonsmooth convex optimization problem, known under the name of Graphical Lasso [ref1]_, given by

.. math::
   \min_{\Theta \in \mathbb{S}^p_{++}} - \log \det \Theta + \mathrm{Tr}(S\Theta) + \lambda \|\Theta\|_{1,od}

where :math:`\mathbb{S}^p_{++}` is the cone of symmetric positive definite matrices of size :math:`p \times p`, :math:`\|\cdot\|_{1,od}` is the sum of absolute values of all off-diagonal elements and :math:`\mathrm{Tr}` is the trace of a matrix.

SGL with latent variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In presence of latent variables, the precision matrix of the marginal of observable variables turns out to have the structure *sparse - low rank* [ref4]_. 

.. image:: pictures/SLRDecomp.png
  :width: 800
  :alt: Illustration of the precision matrix

The problem can then be formulated as  

.. math::
   \min_{\Theta, L \in \mathbb{S}^p} - \log \det (\Theta -L) + \mathrm{Tr}(S(\Theta-L)) + \lambda_1 \|\Theta\|_{1,od} + \mu_1 \|L\|_{\star}

where :math:`\|\cdot\|_{\star}` is the nuclear norm (sum of singular values). :math:`\Theta` represents the sparse part while :math:`L` encodes the low rank component. In [ref12]_, latent variable Graphical Lasso has been studied in the context of microbial abundance analysis and microbiome interaction networks.

Multiple Graphical Lasso problems (MGL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, observations of the same variables but coming from different distributions or in a temporal context are available. For example, consider gene expression data coming from cancer tissue samples and normal tissue samples [ref2]_.
Hence, there has been an increased interest in estimating precision matrices for multiple instances jointly [ref2]_, [ref3]_. Mathematically, we consider :math:`K` Gaussians

.. math::
   \mathcal{X}^{(k)} \sim \mathcal{N}(\mu^{(k)}, \Sigma^{(k)})\in \mathbb{R}^{p}.


Group Graphical Lasso (GGL) describes the problem of estimating precision matrices across multiple instances of the same class under the assumption that the sparsity patterns are similar.
On the other hand, for time-varying data Fused Graphical Lasso was designed in order to get time-consistent estimates of the precision matrices, i.e. :math:`K` is the temporal index.

More generally, the problem formulation of Multiple Graphical Lasso is given by

.. math::
   \min_{\Theta}\quad \sum_{k=1}^{K} \left(-\log\det(\Theta^{(k)}) + \langle S^{(k)},  \Theta^{(k)} \rangle \right)+ \mathcal{P}(\Theta).

In the above, the feasible set is the :math:`K`-fold product of :math:`\mathbb{S}^p_{++}`, see an illustration below. We denote an element of this space :math:`\Theta =  (\Theta^{(1)}, \dots , \Theta^{(K)})`. As input, we have the empirical covariance matrices :math:`S =  (S^{(1)}, \dots , S^{(K)})` given.

.. image:: pictures/sparse_K.png
  :width: 600
  :alt: Illustration of the precision matrix

Group Graphical Lasso (GGL)
""""""""""""""""""""""""""""""""""""""""""""""""""  

For the Group Graphical Lasso problem, the regularization function is given by 

.. math::
   \mathcal{P}(\Theta) = \lambda_1 \sum_{k=1}^{K} \sum_{i \neq j} |\Theta_{ij}^{(k)}| + \lambda_2  \sum_{i \neq j} \left(\sum_{k=1}^{K} |\Theta_{ij}^{(k)}|^2 \right)^{\frac{1}{2}}

with positive numbers :math:`\lambda_1, \lambda_2`. The first term promotes off-diagonal sparsity of the estimator while the second term -- similar to the classical group penalty -- induces that the non-zero entries are present for all instances :math:`\Theta^{(k)}`.

Fused Graphical Lasso (FGL)
"""""""""""""""""""""""""""""""""""""""""""""""""" 

For the Fused Graphical Lasso problem, the regularization function is given by 

.. math::
   \mathcal{P}(\Theta) = \lambda_1 \sum_{k=1}^{K} \sum_{i \neq j} |\Theta_{ij}^{(k)}| + \lambda_2  \sum_{k=2}^{K}   \sum_{i \neq j} |\Theta_{ij}^{(k)} - \Theta_{ij}^{(k-1)}|

with positive numbers :math:`\lambda_1, \lambda_2`. The first term promotes off-diagonal sparsity of the estimator while the second term -- also known as total-variation penalty -- induces that subsequent estimates of :math:`\Theta^{(k)}` are similar.

MGL with latent variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Analogous to SGL, we can extend MGL problems with latent variables.  The problem formulation then becomes 

.. math::
   \min_{\Theta, L}\quad \sum_{k=1}^{K} \left(-\log\det(\Theta^{(k)}- L^{(k)}) + \langle S^{(k)},  \Theta^{(k)} - L^{(k)} \rangle \right)+ \mathcal{P}(\Theta) +\sum_{k=1}^{K} \mu_{1,k} \|L^{(k)}\|_{\star}.

Software is already available for FGL (with and without latent variables) in [ref3]_, however also the deviation of the low rank matrices is included in the penalty.

GGL - the nonconforming case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far, we have assumed that each component of :math:`\mathcal{X}^{(k)}` is present in each of the :math:`K` instances. However, in many practical situations this will not be the case. For example, assume that we have :math:`K` datasets of microbiome abundances but not every microbiome species (OTU) was measured in each dataset. Hence, we may want to estimate the association network but with a group sparsity penalty on all overlapping pairs of species. 

Consequently, assume that we have :math:`\mathcal{X}^{(k)} \sim \mathcal{N}(\mu^{(k)}, \Sigma^{(k)})\in \mathbb{R}^{p_k}` and that there exist groups of overlapping pairs of variables :math:`G_1, \dots, G_L` with 

.. math::
	G_l = \{(i_l^k, j_l^k) \in \mathbb{N}^2 \vert k \in K_l \}, \quad K_l \subset K

where :math:`k \in K_l` if and only if the pair of variables corresponding to :math:`G_l` exists in :math:`\mathcal{X}^{(k)}`. In that case :math:`(i_l^k, j_l^k)` are the indices of the relevant entry in :math:`\Theta^{(k)}` for group :math:`l`.

Now, the associated GGL regularizer becomes 

.. math::
	\mathcal{P}(\Theta) = \lambda_1 \sum_{i \neq j, k} |\Theta_{ij}^{(k)}| + \lambda_2 \sum_{l=1}^{L}\beta_l \|\Theta_{[l]}\|

where 
:math:`\Theta_{[l]}` is the vector with entries :math:`\{\Theta_{i_l^k j_l^k}^{(k)} \vert~ k \in K_l\} \in \mathbb{R}^{|K_l|}`. The scaling factor :math:`\beta_l > 0` is set to :math:`\beta_l = \sqrt{|K_l|}` in order to account for distinct group sizes.

In ``GGLasso`` we implemented an ADMM algorithm for the above described problem formulation, possibly extended with latent variables. Have a look at the :ref:`Nonconforming Group Graphical Lasso experiment` in our example gallery.

Optimization algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All of the above problem formulations are instances of nonlinear, convex and nonsmooth optimization problems. See :ref:`Algorithms` for an overview of solvers which we implemented for these problems and a short guide on how to use them.

References
^^^^^^^^^^^

.. [ref1]  Friedman, J., Hastie, T., and Tibshirani, R. (2007).  Sparse inverse covariance estimation with the Graphical Lasso. Biostatistics, 9(3):432–441.
.. [ref2]  Danaher, P., Wang, P., and Witten, D. M. (2013). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2):373–397.
.. [ref3] Tomasi, F., Tozzo, V., Salzo, S., and Verri, A. (2018). Latent Variable Time-varying Network Inference. InProceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining. ACM.
.. [ref4]  Chandrasekaran, V., Parrilo, P. A., and Willsky, A. S. (2012). Latent variable graphical model selection via convex optimization. The Annals of Statistics, 40(4):1935–1967.
.. [ref5] Ma,  S., Xue,  L., and Zou, H.  (2013). Alternating Direction Methods for Latent Variable Gaussian Graphical Model Selection. Neural Computation, 25(8):2172–2198.
.. [ref6] Zhang, Y., Zhang, N., Sun, D., and Toh, K.-C. (2020). A proximal point dual Newton algorithm for solving group graphical Lasso problems. SIAM J. Optim., 30(3):2197–2220.
.. [ref7] Zhang, N., Zhang, Y.,  Sun, D., and  Toh, K.-C. (2019). An efficient linearly convergent regularized proximal point algorithm for fused multiple graphical lasso problems.
.. [ref8] Boyd, S., Parikh, N., Chu, E., Peleato, B., and Eckstein, J. (2011). Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. Found. Trends Mach. Learn., 3(1):1–122.
.. [ref9] Witten, D. M., Friedman, J. H., and Simon, N. (2011). New Insights and Faster Computations for the Graphical Lasso. J. Comput. Graph. Statist., 20(4):892–900.
.. [ref10] Foygel, R. and Drton, M. (2010). Extended Bayesian Information Criteria for Gaussian Graphical Models. In Lafferty, J., Williams, C., Shawe-Taylor, J.,Zemel, R., and Culotta, A., editors, Advances in Neural Information Processing Systems, volume 23. Curran Associates, Inc.
.. [ref11] Condat, L. (2013). A Direct Algorithm for 1-D Total Variation Denoising, IEEE Signal Processing Letters, vol. 20, no. 11, pp. 1054-1057.
.. [ref12] Kurtz,  Z.  D.,  Bonneau,  R.,  and  Mueller,  C.  L.  (2019).   Disentangling microbial associations from hidden environmental and technical factors via latent graphical models.



Detailled solver documentation
=====================================

SGL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: gglasso.solver.single_admm_solver.ADMM_SGL
.. autofunction:: gglasso.solver.single_admm_solver.block_SGL


MGL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: gglasso.solver.admm_solver.ADMM_MGL

.. autofunction:: gglasso.solver.ppdna_solver.PPDNA


GGL nonconforming
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: gglasso.solver.ext_admm_solver.ext_ADMM_MGL


Using the problem object
=============================

If you want to solve a (Multiple) Graphical Lasso problem, you can of course use the solvers we list in :ref:`Algorithms` directly. However, in most situations it is not clear how to choose the regularization parameters a priori and thus model selection becomes necessary. Below, we describe the model selection functionalities implemented in ``GGLasso``. In order to make its usage as simple as possible, we implemented a class ``glasso_problem`` which calls the solvers/model selection procedures internally and returns an estimator of the precision matrix/matrices in ``sklearn``-style.

Class glasso_problem
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: gglasso.problem.glasso_problem
.. automethod:: gglasso.problem.glasso_problem.solve
.. automethod:: gglasso.problem.glasso_problem.model_selection

Other methods of glasso_problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: gglasso.problem.glasso_problem.set_modelselect_params
.. automethod:: gglasso.problem.glasso_problem.set_reg_params
.. automethod:: gglasso.problem.glasso_problem.set_start_point

Class GGLassoEstimator
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: gglasso.problem.GGLassoEstimator	


Model selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Choosing the regularization parameters :math:`\lambda_1` and :math:`\lambda_2` (and :math:`\mu_1` in the latent variable case) has crucial impact how well the Graphical Lasso solution recovers true edges/ non-zero entries of the precision matrix.

``GGLasso`` contains model selection functionalities for each of the problem described in :ref:`Mathematical description`. Model selection is done via grid searches on the regularization parameters where the quality of a solution is assessed either with the AIC (Akaike Information Criterion) or the eBIC (Extended Bayesian Information Criterion).

Typically, the eBIC chooses sparser solutions and thus leads to less false discoveries. For a single precision matrix estimate :math:`\hat{\Theta}` of dimension :math:`p` with :math:`N` samples it is given by

.. math::
   eBIC_\gamma(\lambda_1) = N \left(-\log \det \hat{\Theta} + \langle S, \hat{\Theta}\rangle\right) + E \log N + 4 E \gamma \log p

where :math:`E` is the number of non-zero off-diagonal entries on the upper triangular of :math:`\hat{\Theta}` and :math:`\gamma \in [0,1]`. The larger you choose :math:`\gamma`, the sparser the solution will become. For MGL problems, we extend this (according to [ref2]_) by

.. math::
	eBIC_\gamma(\lambda_1, \lambda_2) := \sum_{k=1}^{K} N_k \left(-\log \det \hat{\Theta}^{(k)} + \langle S^{(k)}, \hat{\Theta}^{(k)}\rangle\right) + E_k \log N_k + 4 E_k \gamma \log p_k

where :math:`E_k` is the analogous of :math:`E` for :math:`\hat{\Theta}^{(k)}`. We refer to [ref10]_ for details on the eBIC.
Algorithms
=============================

The ``GGLasso`` package contains solvers for several (Multiple) Graphical Lasso problem formulations. See :ref:`Mathematical description` for an overview of problem formulations.
This page aims to give an overview of the functions you need to call for solving a certain problem.

A popular algorithm for solving SGL and MGL problems is the ADMM [ref8]_, [ref5]_, [ref2]_, [ref3]_. Alternatively, a proximal point dual Newton algorithm (PPDNA) was proposed for GGL [ref6]_ and FGL [ref7]_.

The ``GGLasso`` package contains and ADMM solver for all problem formulations as well as the PPDNA solver for MGL problems without latent variables. For a detailled documentation of all of the solvers described below, see :ref:`Detailled solver documentation`.


**Note**: all solvers need as input an empirical covariance matrix (or a collection of matrices for MGL) and a starting point. Using the identity matrix as starting point, if no better guess is available, typically works fine. 

SGL solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For SGL problems, the standard ADMM is ``from gglasso.solver.single_admm_solver import ADMM_SGL``. However, it is shown in [ref9]_ that if the true precision matrix is block-sparse, it is sufficient to solve Graphical Lasso independently on each block. This gives a huge performance improvement if the number of nodes :math:`p` is large because ADMM scales in :math:`\mathcal{O}(p^3)`. We refer to [ref9]_ and our :ref:`Benchmarking` example for details. The resulting algorithm is implemented in ``from gglasso.solver.single_admm_solver import block_SGL``.

``ADMM_SGL`` can also solve latent variable Graphical Lasso problems via setting the option ``latent=True`` and specifying a positive value for the option ``mu1``, the penalty parameter for the nuclear norm.



MGL solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For MGL problems without latent variables two solvers are available, namely 

* ADMM: use ``from gglasso.solver.admm_solver import ADMM_MGL``
* PPDNA: use ``from gglasso.solver.ppdna_solver import PPDNA``. A two-stage version with ADMM for initial iterations and PPDNA for local convergence is implemented in ``warmPPDNA``.

Note that the formulation of the FGL regularizer differs sligtly in [ref2]_, [ref3]_ and our implementation.

Both ADMM and PPDNA have the option ``reg`` which can be set either to ``reg = 'GGL'`` for Group Graphical Lasso problems or ``reg = 'FGL'`` for Fused Graphical Lasso problems. 


For MGL problems with latent variables, only the ADMM solver is available. 

Nonconforming GGL 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the setting which we describe in :ref:`GGL - the nonconforming case`, we implemented an ADMM solver, see ``from gglasso.solver.ext_admm_solver import ext_ADMM_MGL``.
This solver is slightly more complicated to call as you have to tell the solver where the overlapping pairs of variables can be found in the repsective precision matrices. This can be done with the function argument ``G`` which can be seen as a bookeeping array: you should specify a ``(2,L,K)``-shaped array where :math:`L` is the number of groups. 

If your data samples are a list of ``pd.DataFrame`` objects where each Dataframe has the shape ``(n_variables,n_samples)`` and the index contains **unique identifiers** (preferably integers) for all variables, you can create ``G`` by simply calling the following two functions from ``gglasso.helper.ext_admm_helper``.

.. code-block:: python

     ix_exist, ix_location = construct_indexer(list_of_samples) 
     G = create_group_array(ix_exist, ix_location)

Here, ``list_of_samples`` stands for your list of data samples as described above.

We recommend to have a look at the :ref:`Nonconforming Group Graphical Lasso experiment` in our example gallery.


Further Remarks - proximal operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ADMM for Graphical Lasso relies on the efficient computation of proximal operators, in particular for the log-determinant function and the regularization function. For a convex function :math:`f` its proximal operator is given by

.. math::
	\mathrm{prox}_f(x) = \arg \min_z f(z) + \frac{1}{2} \|z-x\|^2.

For solving the FGL problem with ADMM, the proximal operator of the total variation norm (TV), i.e. :math:`x\mapsto \sum_{i=1}^{N-1} |x_{i+1} - x_i|`, needs to be computed. 
This is a case where the proximal operator is not known in closed form. An efficient algorithm for the TV-prox is Condat's algorithm [ref11]_. ``GGLasso`` contains an efficient implementation of this algorithm using ``numba``. This is implemented in ``prox_tv`` from ``gglasso.solver.ggl_helper``. Further, according to [ref7]_ and [ref6]_ we implemented the Clarke differential of this and other proximal operators (for example the :math:`\ell_1`-norm, :math:`\ell_2`-norm and a weighted sum of these norms). For details, see ``gglasso.solver.ggl_helper``... GGLasso documentation master file, created by
   sphinx-quickstart on Fri Apr 16 15:11:28 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GGLasso's documentation!
===================================

``GGLasso`` (stands for **G**\ eneral **G**\ raphical **Lasso**\ ) is a Python package for estimating sparse (or sparse plus low rank) inverse covariance matrices. Moreover, it contains solvers and model selection procedures for Multiple Graphical Lasso problems such as Group and Fused Graphical Lasso.

The package is available on `Github`_ .

.. _Github: https://github.com/fabian-sp/GGLasso

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting-started
   auto_examples/index
   math-description
   problem-object
   solvers-overview
   solvers-doc
   

