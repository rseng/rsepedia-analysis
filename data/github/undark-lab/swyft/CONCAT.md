# Community Guidelines and Contributing
[![codecov](https://codecov.io/gh/undark-lab/swyft/branch/master/graph/badge.svg?token=E253LRJWWE)](https://codecov.io/gh/undark-lab/swyft)
[![Tests](https://github.com/undark-lab/swyft/actions/workflows/tests.yml/badge.svg)](https://github.com/undark-lab/swyft/actions)
[![Syntax](https://github.com/undark-lab/swyft/actions/workflows/syntax.yml/badge.svg)](https://github.com/undark-lab/swyft/actions)
[![Documentation Status](https://readthedocs.org/projects/swyft/badge/?version=latest)](https://swyft.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

If you'd like to contribute to *swyft* or report a bug, please take a moment to read the guide on contributing to [open source software](https://opensource.guide/how-to-contribute/).

If you are a *swyft* user, we welcome a report of your experience. Simple bugs and minor changes could be directly added as [github issues](https://github.com/undark-lab/swyft/issues). Larger issues regarding your particular simulator or computational environment may require [more discussion](https://github.com/undark-lab/swyft/discussions).

## Report issues, bugs, or problems with *swyft*

Please be patient as *swyft* is research software but it is being actively developed. Small issues or bugs directly related to *swyft* code justify an issue. If the problem is specific to your simulator, then we will need a lot more detail about the simulator and its setup.

Please reference the version of *swyft* you are using in any issues / bug reports. Check out [issues on GitHub](https://github.com/undark-lab/swyft/issues).

## Contribute to *swyft*

We'd be happy to include your contribution! We operate by introducing issues and solving them with a corresponding pull request to address the issue.

### Setup Envrionment

Unless your contribution is purely documentation based, you will need to setup a development version of *swyft*. Create the environment, including pre-commit hooks.

```bash
git clone https://github.com/undark-lab/swyft.git
cd swyft
pip install -e .[dev]
pre-commit install
```

The :code:`-e` flag will install *swyft* in development mode such that your version of the code is used when *swyft* is imported.
The :code:`[dev]` flag installs the extra tools necessary to format and test your contribution.
`pre-commit` will enforce [black](https://github.com/psf/black),
[isort](https://github.com/timothycrosley/isort),
and a few other code format rules before every commit.

### Testing
Any code that ends up in the master branch must be tested. It is simple to test the software by running the following command from the root directory:

```bash
pytest tests/
```

If you're introducing new functionality, we ask that you create tests so that our code coverage percentage may remain high. The best tests are comprehensive, using `pytest.mark.parameterize` to cover various cases where the functionality may be called. New tests should be fast computationally, i.e., training a neural network should not be part of the test suite.

### Code Standards

#### Docstrings
- Please use [Google Style](http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for docstrings and comments.
- When we create a class, we put the docstring in the class rather than the `__init__` function, where appropriate.
- Use type annotations rather than putting the type in the docstring.
- When there is a default argument, do not repeat the default argument within the docstring since it is already visible when the user calls help on the function. The exception to this rule is when the default is `None`, then it may require an explanation which may include a default argument.

#### Types
Contributed code should have type hints. Some relevant types are defined within `swyft/types.py`, try to use them whenever possible.

#### Naming
- Integers which count the quantity of something should be proceded with an `n_*`, i.e. `n_parameters`.
- Although brevity is appreciated...
- Long names are generally more useful than unclear / single-use shortened versions. If you introduce a shortened version, please make sure it is consistent throughout your code and the existing codebase.
- When introducing a new variable, consider whether  it already has a name... see the table below.

For quick / common naming reference:

| Python Variable Name | Mathematical Object                           |
|----------------------|-----------------------------------------------|
| `v`                  | parameter vector in natural prior coordinates |
| `u`                  | parameter vector mapped to the hypercube      |
| `marginal_index`     | `tuple` of integers, often a key in a dict    |
| `marginal_indices`   | `tuple` of `marginal_index`                   |
| `observation`        | `dict` maps to simulated / observational data |
| `*_o`                | Either `v` or `observation` of interest       |

#### Converting between arrays and tensors
Please use the functions `array_to_tensor` and `tensor_to_array` when converting between arbitrary array data and pytorch tensors. This is to maintain consistency with default conversion types.

## Documentation

We have a [readthedocs](https://swyft.readthedocs.io/en/latest/) site and use the [sphinx_rtd_theme](https://github.com/readthedocs/sphinx_rtd_theme). The details of the configuration can be found in the [docs/source/conf.py](https://github.com/undark-lab/swyft/blob/master/docs/source/conf.py) file.

### Compiling documentation

To compile your own version of the documentation, run the following commands:

```bash
cd docs
make html
```

This will produce an html version of the documentation which can easily be viewed on a web browser by pointing to `swyft/docs/build/html/index.html`.
*swyft*
=======

.. image:: https://badge.fury.io/py/swyft.svg
   :target: https://badge.fury.io/py/swyft
   :alt: PyPI version


.. .. image:: https://github.com/undark-lab/swyft/actions/workflows/tests.yml/badge.svg
..    :target: https://github.com/undark-lab/swyft/actions
..    :alt: Tests


.. .. image:: https://github.com/undark-lab/swyft/actions/workflows/syntax.yml/badge.svg
..    :target: https://github.com/undark-lab/swyft/actions
..    :alt: Syntax


.. image:: https://codecov.io/gh/undark-lab/swyft/branch/master/graph/badge.svg?token=E253LRJWWE
   :target: https://codecov.io/gh/undark-lab/swyft
   :alt: codecov


.. .. image:: https://readthedocs.org/projects/swyft/badge/?version=latest
..    :target: https://swyft.readthedocs.io/en/latest/?badge=latest
..    :alt: Documentation Status


.. .. image:: https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat
..    :target: https://github.com/undark-lab/swyft/blob/master/CONTRIBUTING.md
..    :alt: Contributions welcome


.. .. image:: https://colab.research.google.com/assets/colab-badge.svg
..    :target: https://colab.research.google.com/github/undark-lab/swyft/blob/master/notebooks/Quickstart.ipynb
..    :alt: colab

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5752734.svg
   :target: https://doi.org/10.5281/zenodo.5752734

*swyft* is the official implementation of Truncated Marginal Neural Ratio Estimation (TMNRE),
a hyper-efficient, simulation-based inference technique for complex data and expensive simulators.

* **Documentation & installation**: https://swyft.readthedocs.io/en/latest/
* **Example usage**: https://swyft.readthedocs.io/en/latest/tutorial-notebooks.html
* **Source code**: https://github.com/undark-lab/swyft
* **Support & discussion**: https://github.com/undark-lab/swyft/discussions
* **Bug reports**: https://github.com/undark-lab/swyft/issues
* **Contributing**: https://swyft.readthedocs.io/en/latest/contributing-link.html
* **Citation**: https://swyft.readthedocs.io/en/latest/citation.html

*swyft*:

* estimates likelihood-to-evidence ratios for arbitrary marginal posteriors; they typically require fewer simulations than the corresponding joint.
* performs targeted inference by prior truncation, combining simulation efficiency with empirical testability.
* seamless reuses simulations drawn from previous analyses, even with different priors.
* integrates `dask <https://dask.org/>`_ and `zarr <https://zarr.readthedocs.io/en/stable/>`_ to make complex simulation easy.

*swyft* is designed to solve the Bayesian inverse problem when the user has access to a simulator that stochastically maps parameters to observational data.
In scientific settings, a cost-benefit analysis often favors approximating the posterior marginality; *swyft* provides this functionality.
The package additionally implements our prior truncation technique, routines to empirically test results by estimating the expected coverage,
and a `dask <https://dask.org/>`_ simulator manager with `zarr <https://zarr.readthedocs.io/en/stable/>`_ storage to simplify use with complex simulators.



Related
-------

* `tmnre <https://github.com/bkmi/tmnre>`_ is the implementation of the paper `Truncated Marginal Neural Ratio Estimation <https://arxiv.org/abs/2107.01214>`_.
* `v0.1.2 <https://github.com/undark-lab/swyft/releases/tag/v0.1.2>`_ is the implementation of the paper `Simulation-efficient marginal posterior estimation with swyft: stop wasting your precious time <https://arxiv.org/abs/2011.13951>`_.
* `sbi <https://github.com/mackelab/sbi>`_ is a collection of simulation-based inference methods. Unlike *swyft*, the repository does not include truncation nor marginal estimation of posteriors.

Educational videos
==================

Tutorial videos
---------------

These videos were created for an older version of swyft; however, they may remain interesting to users. Check back and we will update them to the new version and expand them!

* `1 - Linear regression with swyft <https://www.loom.com/share/cefac9e4e84d482c89c5281b90121974>`_
* `2 - ConvNets, parameter regression and swyft <https://www.loom.com/share/1fc4785159bf4f0081e59693133a5ad3>`_
* ...more to come...

Academic videos
---------------

* `Truncated Marginal Neural Ratio Estimation <https://www.youtube.com/watch?v=euUxDdB5XY8>`_
Why *swyft*?
============

Overview
--------

With *swyft* :cite:p:`miller2020simulation` our goal is to provide a general, flexible, reliable and practical
tool for solving hard Bayesian parameter inference problems in physics and
astronomy.  *swyft* uses a specific flavor of simulation-based neural inference
techniques called Truncated Marginal Neural Ratio Estimation :cite:p:`miller2021truncated`, that offers
multiple advantages over established Markov Chain based methods, or other
simulation-based neural approaches. It is based on the technique presented in :cite:t:`hermans2020likelihood`.

- *swyft* directly estimates marginal posteriors, which typically requires far
  less simulation runs than estimating the full joint posterior.
- *swyft* uses a simulation store that make re-use of simulations, even with
  different priors, efficient and seamless.
- *swyft* performs targeted inference by prior truncation, which combines the
  simulation efficiency of existing sequential methods with the testability of
  amortized methods.

Details
-------

Marginal posterior estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*swyft* can directly estimate marginal posteriors for parameters of interest
:math:`\boldsymbol{\vartheta}`, given some observation :math:`\mathbf{x}`. These are
formally obtained by integrating over all remaining nuisance parameters
:math:`\boldsymbol{\eta}`,

.. math::
   p(\boldsymbol{\vartheta}|\mathbf{x}) = \frac{\int d\boldsymbol{\eta}\,
   p(\mathbf{x}|\boldsymbol{\vartheta}, \boldsymbol{\eta}) p(\boldsymbol{\eta}, \boldsymbol{\vartheta})}
   {p(\mathbf{x})}\;.

Here, :math:`p(\mathbf{x}|\boldsymbol{\vartheta}, \boldsymbol{\eta})` is an abritrary
forward model that includes both the physics and detector simulator,
:math:`p(\boldsymbol{\vartheta}, \boldsymbol{\eta})` is the joint prior,
and :math:`p(\mathbf{x})` is the Bayesian evidence.


Nuisance parameters
^^^^^^^^^^^^^^^^^^^

*In the context of likelihood-based inference, nuisance parameters are an
integration problem.* Given the likelihood density
:math:`p(\mathbf{x}|\boldsymbol{\vartheta}, \boldsymbol{\eta})` for a particular
observation :math:`\mathbf{x}`, one attempts to solve the above integral over
:math:`\boldsymbol{\eta}`, e.g. through sampling based methods.  This becomes
increasingly challenging if the number of nuisance parameters grows.

*In the context of likelihood-free inference, nuisance parameters are noise.*
Posteriors are estimated based on a large number of training samples
:math:`\mathbf{x}, \boldsymbol{\vartheta}\sim p(\mathbf{x}|\boldsymbol{\vartheta},
\boldsymbol{\eta})p(\boldsymbol{\vartheta}, \boldsymbol{\eta})`, no matter the dimension
of the nuisance parameter space. For a given :math:`\boldsymbol{\vartheta}`, more nuisance
parameters just increase the variance of :math:`\mathbf{x}` (which oddly enough
can make the inference problem simpler rather than more difficult).


Simulation re-use
^^^^^^^^^^^^^^^^^

*Likelihood-based techniques often use Markov chains*, which require a simulation
for every link in the chain. Due to the properties of Markov chains, it is not
possible to utilize those simulations again for further analysis.
That effort has been lost.

*Likelihood-free inference can be based on simulations that sample the
(constrained) prior*. Reusing these simulations is allowed, we don’t
have to worry about breaking the Markov chain.


High precision through targeted inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Likelihood-based techniques are highly precise by focusing simulator
runs on parameter space regions that are consistent with a particular
observation.

Likelihood-free inference techniques can be less precise when there are
too few simulations in parameter regions that matter most.


Learned features
^^^^^^^^^^^^^^^^

*swyft* uses neural likelihood estimation. The package works out-of-the-box for
low-dimensional data.  Tackling complex and/or high-dimensional data (e.g.,
high-resolution images or spectra, combination of multiple data sets) is
possible through providing custom feature extractor networks in pytorch.

.. bibliography::
swyft API
=====================

.. toctree::
  :glob:

  module/*
.. swyft documentation master file, created by
   sphinx-quickstart on Thu Nov 12 11:23:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../../README.rst

.. toctree::
   :caption: Contents:
   :maxdepth: 2

   swyft
   install
   tutorial-notebooks
   videos
   autodoc
   contributing-link
   citation
Citation
-------------
If you use *swyft* in scientific publications, please include both of the following citations.

.. code-block:: bibtex

   @article{miller2021truncated,
      title={Truncated Marginal Neural Ratio Estimation},
      author={Miller, Benjamin Kurt and Cole, Alex and Forr{\'e}, Patrick and Louppe, Gilles and Weniger, Christoph},
      journal={Advances in Neural Information Processing Systems},
      volume={34},
      year={2021}
   }

   @article{miller2020simulation,
      title={Simulation-efficient marginal posterior estimation with swyft: stop wasting your precious time},
      author={Miller, Benjamin Kurt and Cole, Alex and Louppe, Gilles and Weniger, Christoph},
      journal={Machine Learning and the Physical Sciences: Workshop at the 34th Conference on Neural Information Processing Systems (NeurIPS)},
      year={2020}
   }
.. mdinclude:: ../../CONTRIBUTING.md
Installation
===============


pip
--------
**After installing** `pytorch <https://pytorch.org/get-started/locally/>`_, please run the command:

.. code-block:: bash

  pip install swyft



github
---------
If you want the lastest version, it is also possible to install *swyft* from github.
To do so, run the command:

.. code-block:: bash

  pip install git+https://github.com/undark-lab/swyft.git


develop
---------
If you're interested in contributing to swyft there is another procedure.
First clone the github repo, navigate to the repo in your terminal, from within that directory run the command:

.. code-block:: bash

  pip install -e .[dev]

The :code:`-e` flag will install *swyft* in development mode such that your version of the code is used when *swyft* is imported.
The :code:`[dev]` flag installs the extra tools necessary to format and test your contribution.


docs
---------
Compiling the docs requires an additional flag. Then the docs may be compiled by navigating to the docs folder.

.. code-block:: bash

  pip install -e .[docs]
  cd docs
  make html


.. _pytorch: https://pytorch.org/get-started/locally/
Tutorial notebooks
==================

You may find the tutorial notebooks here.

.. toctree::
   :caption: Notebooks:
   :maxdepth: 1
   :glob:

   notebooks/Quickstart
   notebooks/Examples*
swyft.bounds
---------------------------

.. automodule:: swyft.bounds
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.inference
-------------------

The primary marginal objects are contained here.

.. automodule:: swyft.inference
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.weightedmarginals
---------------------------

.. automodule:: swyft.weightedmarginals
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.utils
-----------------

.. automodule:: swyft.utils
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.prior
---------------------------

.. automodule:: swyft.prior
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.plot
-----------------

.. automodule:: swyft.plot
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
swyft.networks
--------------------

.. automodule:: swyft.networks
  :members:
  :special-members: __len__, __getitem__
swyft.store
------------

.. automodule:: swyft.store
  :members:
  :inherited-members:
  :special-members: __len__, __getitem__
Theoretical concepts
====================

Introduction
------------

Parametric stochastic simulators are ubiquitous in the physical sciences
 [1–3]. However, performing parameter inference based on simulator runs
using Markov chain Monte Carlo is inconvenient or even impossible if the
model parameter space is large or the likelihood function is
intractable. This problem is addressed by so-called likelihood-free
inference  [4] or simulation-based inference  [5] techniques. Deep
learning based likelihood-free inference algorithms were organized into
a taxonomy in Ref.  [6], where methods that estimated likelihood ratios
in a series of rounds were denoted Sequential Ratio Estimation
(SRE)  [7]. Our presented method is closely related.

We propose *Nested Ratio Estimation* (NRE), which approximates the
likelihood-to-evidence ratio in a sequence of rounds. Loosely inspired
by the contour sorting method of nested sampling  [8–10], the scheme
alternates between sampling from a constrained prior and estimating
likelihood-to-evidence ratios. It allows for efficient estimation of any
marginal posteriors of interest. Furthermore, we propose an algorithm
that we call *iP3 sample caching*, which facilitates simulator
efficiency by automatizing the reuse of previous simulator runs through
resampling of cached simulations.

The primary use case for these algorithms is the calculation of
arbitrary, low-dimensional marginal posteriors, typically in one or two
dimensions. In physics and astronomy, such marginals serve as the basis
for scientific conclusions by constraining individual model parameters
within uncertainty bounds. We implement a multi-target training regime
where all marginal posteriors of interest can be learned simultaneously.
We find that learning is simplified when one calculates each marginal
distribution directly rather than computing the full joint posterior and
marginalizing numerically. Furthermore, the method facilitates
effortless marginalization over arbitrary numbers of nuisance
parameters, increasing its utility in high-dimensional parameter
regimes–even to simulators with a tractable, yet high-dimensional,
likelihood  [11].

Nested Ratio Estimation (NRE)
-----------------------------

We operate in the context of simulation-based inference where our
simulator :math:`\mathbf{g}` is a nonlinear function mapping a vector of
parameters
:math:`\boldsymbol{\theta}= (\theta_{1}, \dots, \theta_{d}) \in \mathbb{R}^{d}`
and a stochastic latent state :math:`\mathbf{z}` to an observation
:math:`\mathbf{x}= \mathbf{g}(\boldsymbol{\theta}, \mathbf{z})`. The
likelihood function is therefore
:math:`p(\mathbf{x}\vert \boldsymbol{\theta}) = \int \delta(\mathbf{x}- \mathbf{g}(\boldsymbol{\theta}, \mathbf{z})) \, p(\mathbf{z}\vert \boldsymbol{\theta}) \, d\mathbf{z}`,
with :math:`\delta(\cdot)` denoting the Dirac delta. Consider a
factorizable prior
:math:`p(\boldsymbol{\theta}) = p(\theta_{1}) \cdots p(\theta_{d})` over
the parameters, the joint posterior is given via Bayes’ rule as
:math:`p(\boldsymbol{\theta}|\mathbf{x}) = p(\mathbf{x}|\boldsymbol{\theta})p(\boldsymbol{\theta})/p(\mathbf{x})`,
where :math:`p(\mathbf{x})` is the evidence.

Our goal is to compute the marginal posterior,
:math:`p(\boldsymbol{\vartheta}\vert \mathbf{x})`, where
:math:`\boldsymbol{\vartheta}` are the parameters of interest. We will
denote all other parameters by :math:`\boldsymbol{\eta}`, such that
:math:`\boldsymbol{\theta}= (\boldsymbol{\vartheta}, \boldsymbol{\eta})`.
The marginal posterior is obtained from the joint distribution
:math:`p(\boldsymbol{\vartheta}, \boldsymbol{\eta}|\mathbf{x}) \equiv p(\boldsymbol{\theta}|\mathbf{x})`
by integrating over all components of :math:`\boldsymbol{\eta}`,

.. math::

   \label{eqn:post}
   p(\boldsymbol{\vartheta}\vert \mathbf{x})  \equiv \int p(\boldsymbol{\vartheta}, \boldsymbol{\eta}| \mathbf{x}) d\boldsymbol{\eta}
   = \int \frac{p(\mathbf{x}| \boldsymbol{\vartheta}, \boldsymbol{\eta})}{p(\mathbf{x})}
   p(\boldsymbol{\theta})
   %\prod_{j \notin \texttt{idx}} d\theta_{j}
   d\boldsymbol{\eta}
   = \frac{p(\mathbf{x}|\boldsymbol{\vartheta})}{p(\mathbf{x})}p(\boldsymbol{\vartheta})\;,

where we used Bayes’ rule and defined the marginal likelihood
:math:`p(\mathbf{x}|\boldsymbol{\vartheta})` in the last step.

Just like in SRE, we focus on a specific observation of interest,
:math:`\mathbf{x}_0`. Only parameter values :math:`\boldsymbol{\theta}`
that could have plausibly generated observation :math:`\mathbf{x}_0`
will significantly contribute to the integrals in
Eq. `[eqn:post] <#eqn:post>`__. For implausible values the likelihood
:math:`p(\mathbf{x}_0|\boldsymbol{\theta})` will be negligible. We
denote priors that are suitably constrained to plausible parameter
values by :math:`\tilde{p}(\theta_1, \dots, \theta_d)`. Similarly,
:math:`\tilde{\square}` indicates quantities :math:`\square` that are
calculated using the constrained prior. Therefore, using a judiciously
chosen constrained prior, accurately approximates the marginal posterior
in place of our true prior beliefs,

.. math::

   p(\boldsymbol{\vartheta}| \mathbf{x}_0) =
   \frac{p(\mathbf{x}_0|\boldsymbol{\vartheta})}{p(\mathbf{x}_0)} p(\boldsymbol{\vartheta}) \simeq
   \frac{\tilde{p}(\mathbf{x}_0|\boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x}_0)} \tilde{p}(\boldsymbol{\vartheta})\;.

The increased probability that constrained priors assign to the
plausible parameter region cancels when dividing by the constrained
evidence :math:`\tilde p(\mathbf{x})`. We define the marginal
likelihood-to-evidence ratio

.. math::

   \label{eqn:likelihood_ratio}
       \tilde{r}(\mathbf{x}, \boldsymbol{\vartheta})
       \equiv \frac{\tilde{p}(\mathbf{x}\vert \boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x})}
       = \frac{\tilde{p}(\mathbf{x}, \boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x}) \tilde{p}(\boldsymbol{\vartheta})}
       = \frac{\tilde{p}(\boldsymbol{\vartheta}\vert\mathbf{x})}{\tilde{p}(\boldsymbol{\vartheta})}\;,

which is sufficient to evaluate the marginal posterior in
Eq. `[eqn:post] <#eqn:post>`__, and which we will now estimate. Under
the assumption of equal class population, it is known  [6,12] that one
can recover density ratios using binary classification to distinguish
between samples from two distributions. Our binary classification
problem is to distinguish positive samples,
:math:`(\mathbf{x}, \boldsymbol{\vartheta}) \sim \tilde{p}(\mathbf{x}, \boldsymbol{\vartheta}) = p(\mathbf{x}\vert \boldsymbol{\vartheta}) \tilde{p}(\boldsymbol{\vartheta})`,
drawn jointly, and negative samples,
:math:`(\mathbf{x}, \boldsymbol{\vartheta}) \sim \tilde{p}(\mathbf{x}) \tilde{p}(\boldsymbol{\vartheta})`,
drawn marginally. The binary classifier
:math:`\sigma(f_{\phi}(\mathbf{x}, \boldsymbol{\vartheta}))` performs
optimally when
:math:`f_{\phi}(\mathbf{x}, \boldsymbol{\vartheta}) = \log \tilde{r}(\mathbf{x}, \boldsymbol{\vartheta})`,
where :math:`\sigma(\cdot)` is the sigmoid function and :math:`f_{\phi}`
is a neural network parameterized by :math:`\phi`. The associated binary
cross-entropy loss function used to train the ratio
:math:`\tilde{r}(\boldsymbol{\vartheta}, \mathbf{x}_0)` via stochastic
gradient descent is given by

.. math:: \ell = -\int \left[ \tilde{p}(\mathbf{x}|\boldsymbol{\vartheta})\tilde{p}(\boldsymbol{\vartheta}) \ln \sigma(f_\phi(\mathbf{x}, \boldsymbol{\vartheta})) + \tilde{p}(\mathbf{x})\tilde{p}(\boldsymbol{\vartheta}) \ln \sigma(-f_\phi(\mathbf{x},\boldsymbol{\vartheta})) \right] d\mathbf{x}\, d\boldsymbol{\vartheta}\;.

We propose to iteratively improve marginal posterior estimates in
:math:`R` rounds by employing posterior estimates from previous rounds
to define constrained priors. In each round :math:`r`, we estimate *all*
1-dim marginal posteriors, using :math:`d` instances of the above
marginal likelihood-to-evidence ratio estimation in parallel by setting
:math:`\boldsymbol{\vartheta}= (\theta_i)` for :math:`i=1, \dots, d`. To
this end, we utilize the factorized constrained prior,
:math:`\tilde{p}_r(\theta) = \tilde{p}_r(\theta_1)\cdots\tilde{p}_r(\theta_d)`,
which is defined recursively by a cutoff criterion,

.. math::

   \tilde{p}_{r}(\theta_{i})
       \propto
       p(\theta_{i}) \Theta_{H} \left[ \frac{\tilde{r}_{r-1}(\theta_{i}, \mathbf{x})}{\max_{\theta_{i}} \tilde{r}_{r-1}(\theta_{i}, \mathbf{x})} - \epsilon \right],
       \label{eqn:it}

where :math:`\Theta_{H}` denotes the Heaviside step function and
:math:`\epsilon` denotes the minimum likelihood-ratio which passes
through the threshold. We use
:math:`\tilde{p}_1(\boldsymbol{\theta}) = p(\boldsymbol{\theta})` as an
initial prior in the iterative scheme.

In every round, each 1-dim posterior approximates a marginalization of
the same underlying constrained posterior, allowing us to effectively
reuse training data and train efficiently in a multi-target regime. The
inference network is therefore divided into a featurizer
:math:`\mathbf{F}(\mathbf{x})` with shared parameters and a set of
:math:`d` independent Multi-layer Perceptons
:math:`\{\textrm{MLP}_i(\cdot, \cdot)\}_{i=1}^{d}` which estimate
individual 1-dim marginal posteriors and do not share parameters, such
that
:math:`f_{\phi}(\mathbf{x}, \theta_i) = \textrm{MLP}_i(\mathbf{F}(\mathbf{x}), \theta_i)`.

This technique is valid as long as the excluded prior regions do not
significantly affect the integrals in Eq. `[eqn:post] <#eqn:post>`__.
For uncorrelated parameters, a sufficient criterion is that the impact
on the marginal posteriors is small, which we guarantee through the
iteration criterion Eq. `[eqn:it] <#eqn:it>`__. In the case of a very
large number of strongly correlated parameters the algorithm can
inadvertently cut away tails of the marginal posteriors. Decreasing
:math:`\epsilon` mitigates this effect. Discussion is left for future
study  [13].

With this design, the posteriors from the final round can be used to
approximate the true 1-dim marginal posteriors,
:math:`\tilde{p}_{R}(\theta_i \vert \mathbf{x}_{0}) \approx p(\theta_i\vert \mathbf{x}_{0})`,
while previous rounds were used to iteratively focus on relevant parts
of the parameter space. The key result and value of NRE lies in the
utility of our constrained prior from round :math:`R`. The final
constrainted prior, along with previously generated and cached samples,
allows for estimation of *any* higher dimensional marginal posterior
:math:`\tilde{p}_R(\boldsymbol{\vartheta}|\mathbf{x}_0) \approx p(\boldsymbol{\vartheta}|\mathbf{x}_0)`
of interest by doing likelihood-to-evidence ratio estimation, often
without further simulation.

Inhomogeneous Poisson Point Process (iP3) Sample Caching
--------------------------------------------------------

Simulating
:math:`(\mathbf{x}, \boldsymbol{\theta})\sim p(\mathbf{x}|\boldsymbol{\theta})p(\boldsymbol{\theta})`
can be extremely expensive. We develop a scheme to systematically reuse
appropriate subsets of previous simulator runs. Our method samples
:math:`N\sim \text{Pois}(\hat N)` parameter vectors from an arbitrary
distribution :math:`p(\boldsymbol{\theta})`, where :math:`\hat N` is the
expected number of samples. Taking :math:`N` samples from
:math:`p(\boldsymbol{\theta})` is equivalent to drawing a single sample
:math:`\Theta \equiv \{\boldsymbol{\theta}^{(n)}\}_{n=1}^{N}` from an
inhomogenous Poisson point process (PPP) with intensity function
:math:`\lambda_{r}(\boldsymbol{\theta}) = \hat{N} p(\boldsymbol{\theta})`.
In this context, :math:`\Theta` is known as a set of *points*. This
formulation provides convenient mathematical properties  [14], at the
low price of introducing variance in the number of samples drawn. The
precise number of samples doesn’t matter as long as
:math:`N \approx \hat{N}`, which is true in our regime of order
:math:`\geq 1000`.

We will need two properties of PPPs. *Superposition:* Given two
independent PPPs with intensity functions
:math:`\lambda_{1}(\boldsymbol{\theta})` and
:math:`\lambda_{2}(\boldsymbol{\theta})`, the sum yields another PPP
with intensity function
:math:`\lambda(\boldsymbol{\theta}) = \lambda_{1}(\boldsymbol{\theta}) + \lambda_{2}(\boldsymbol{\theta})`.
The union of two sets of points :math:`\Theta = \Theta_1 \cup \Theta_2`
from the individual PPPs is equivalent to a single set of points from
the combined PPP. *Thinning:* Consider a PPP with intensity function
:math:`\lambda(\boldsymbol{\theta})`, and an arbitrary function
:math:`q(\boldsymbol{\theta}): \mathbb{R}^{d} \to [0, 1]`. If we are
interested in drawing from a PPP with intensity function
:math:`\lambda_{q}(\boldsymbol{\theta}) = q(\boldsymbol{\theta}) \lambda(\boldsymbol{\theta})`,
we can achieve this by drawing a set of points :math:`\Theta`
distributed like :math:`\lambda(\boldsymbol{\theta})` and then rejecting
individual points :math:`\boldsymbol{\theta}^{(n)}` with probability
:math:`1 - q(\boldsymbol{\theta}^{(n)})`.

The parameter cache is defined by a set of points :math:`\Theta_{sc}`
drawn from a PPP with intensity function
:math:`\lambda_{sc}(\boldsymbol{\theta})`. For every point
:math:`\boldsymbol{\theta}\in\Theta_{sc}`, a corresponding observation
:math:`\mathbf{x}` is stored in an observation cache
:math:`\mathcal{X}_{sc}`. The iP3 cache sampling algorithm that is
responsible for maintaining the caches and sampling from a PPP with
target intensity function
:math:`\lambda_t(\boldsymbol{\theta}) = \hat{N} p(\boldsymbol{\theta})`
is written out in the supplementary material. It is summarized in two
steps: First, consider all points
:math:`\boldsymbol{\theta}\in \Theta_{sc}` from the cache and accept
them with probability
:math:`\min(1, \lambda_t(\boldsymbol{\theta})/\lambda_{sc}(\boldsymbol{\theta}))`.
The thinning operation yields a sample :math:`\Theta_1` from a PPP with
intensity function
:math:`\lambda_1(\boldsymbol{\theta}) = \min(\lambda_t(\boldsymbol{\theta}), \lambda_{sc}(\boldsymbol{\theta}))`.
Second, draw a new set of points :math:`\Theta_p` from
:math:`\lambda_t(\boldsymbol{\theta})`, and accept each
:math:`\boldsymbol{\theta}\in\Theta_p` with probability
:math:`\max(0, 1-\lambda_{sc}(\boldsymbol{\theta})/\lambda_t(\boldsymbol{\theta}))`.
This yields a sample :math:`\Theta_2` from a PPP with intensity function
:math:`\lambda_2(\boldsymbol{\theta}) = \max(0, \lambda_t(\boldsymbol{\theta}) - \lambda_{sc}(\boldsymbol{\theta}))`.
Thanks to superposition, the union
:math:`\Theta_1 \cup \Theta_2 = \Theta_t` yields a sample from the PPP
with intensity function :math:`\lambda_t(\boldsymbol{\theta})`–the
sample we were looking for. We only need to run simulations on points
from :math:`\Theta_1`. Points in :math:`\Theta_2` already have
corresponding observations in :math:`\mathcal{X}_{sc}` which we can
reuse. Finally, the new parameters are appended to the set of points in
the parameter cache, :math:`\Theta_{sc} \to \Theta_{sc} \cup \Theta_2`.
Similar for :math:`\mathcal{X}_{sc}`. On the basis of the superposition
principle, the intensity function of the :math:`\Theta_{sc}` cache is
updated
:math:`\lambda_{sc}(\boldsymbol{\theta}) \to \max(\lambda_{sc}(\boldsymbol{\theta}), \lambda_t(\boldsymbol{\theta}))`.

Storing and updating the parameter cache’s intensity function
:math:`\lambda_{sc}(\boldsymbol{\theta})` can pose challenges when it is
complex and high-dimensional. Our NRE implementation overcomes these
challenges by learning marginal 1-dim posteriors, guaranteeing that the
relevant target intensities always factorize,
:math:`\lambda_t(\boldsymbol{\theta}) = \lambda_t(\theta_1)\cdots \lambda_t(\theta_d)`.
Storage of and calculation with factorizable functions simplifies
matters.

.. container:: references csl-bib-body
   :name: refs

   .. container:: csl-entry
      :name: ref-Banik_2018

      [1]N. Banik, G. Bertone, J. Bovy, and N. Bozorgnia, Probing the
      Nature of Dark Matter Particles with Stellar Streams, Journal of
      Cosmology and Astroparticle Physics 2018, 061 (2018).

   .. container:: csl-entry
      :name: ref-Bartels_2016

      [2]R. Bartels, S. Krishnamurthy, and C. Weniger, Strong Support
      for the Millisecond Pulsar Origin of the Galactic Center GeV
      Excess, Physical Review Letters 116, (2016).

   .. container:: csl-entry
      :name: ref-Rodr_guez_Puebla_2016

      [3]A. Rodríguez-Puebla, P. Behroozi, J. Primack, A. Klypin, C.
      Lee, and D. Hellinger, Halo and Subhalo Demographics with Planck
      Cosmological Parameters: Bolshoi–Planck and MultiDark–Planck
      Simulations, Monthly Notices of the Royal Astronomical Society
      462, 893 (2016).

   .. container:: csl-entry
      :name: ref-sisson2018handbook

      [4]S. A. Sisson, Y. Fan, and M. Beaumont, Handbook of Approximate
      Bayesian Computation (CRC Press, 2018).

   .. container:: csl-entry
      :name: ref-Cranmer2020

      [5]K. Cranmer, J. Brehmer, and G. Louppe, The Frontier of
      Simulation-Based Inference, Proc. Natl. Acad. Sci. U. S. A.
      (2020).

   .. container:: csl-entry
      :name: ref-Durkan2020

      [6]C. Durkan, I. Murray, and G. Papamakarios, On Contrastive
      Learning for Likelihood-Free Inference, (2020).

   .. container:: csl-entry
      :name: ref-Hermans2019

      [7]J. Hermans, V. Begy, and G. Louppe, Likelihood-Free MCMC with
      Amortized Approximate Ratio Estimators, (2019).

   .. container:: csl-entry
      :name: ref-Skilling2006

      [8]J. Skilling, Nested Sampling for General Bayesian Computation,
      Bayesian Anal. 1, 833 (2006).

   .. container:: csl-entry
      :name: ref-Feroz2008

      [9]F. Feroz, M. P. Hobson, and M. Bridges, MultiNest: An Efficient
      and Robust Bayesian Inference Tool for Cosmology and Particle
      Physics, Mon. Not. Roy. Astron. Soc. 398: 1601-1614,2009 (2008).

   .. container:: csl-entry
      :name: ref-Handley2015

      [10]W. J. Handley, M. P. Hobson, and A. N. Lasenby, Polychord :
      Next-Generation Nested Sampling, Mon. Not. R. Astron. Soc. 453,
      4384 (2015).

   .. container:: csl-entry
      :name: ref-lensing

      [11]A. et. al., Precision Analysis of Gravitational Strong Lensing
      Images with Nested Likelihood-Free Inference, (2020).

   .. container:: csl-entry
      :name: ref-Cranmer2015

      [12]K. Cranmer, J. Pavez, and G. Louppe, Approximating Likelihood
      Ratios with Calibrated Discriminative Classifiers, (2015).

   .. container:: csl-entry
      :name: ref-swyft_future

      [13]A. et. al., Nested Ratio Estimation and iP3 Sample Caching,
      (2020).

   .. container:: csl-entry
      :name: ref-ppp

      [14]J. F. C. Kingman, Poisson Processes (Oxford University Press,
      1993).
Goals and plans
===============


Overview
--------

*swyft* aims at combining the convenience and precision of traditional
likelihood-based inference with the flexibility and power of modern neural
network based likelihood-free methods.  *swyft* is based on the concepts of
Amortized Approximate Likelihood Ratio Estimation (AALR), Nested Ratio
Estimation (NRE), and iP3 sample caching. In short, *swyft* = AALR + NRE +
iP3.


Our goals
---------

We are developing *swyft* to help solve our data analysis problems in
astroparticle physics, and provide a useful software tool for others simultanously.
We have the following goals in mind.

- Convenience: *swyft* utilizes neural networks for likelihood-free inference with
  the simplicity of traditional sampling tools.
- Precision: *swyft* is as precise as likelihood-based methods, where those can be
  applied, and flexible enough to cover the wide range of scenarios where they
  fail.
- Efficiency: *swyft* is simulator efficient, and automatically re-uses previous
  simulator runs in new analyses.

The current version of the code provides a prototype implementation of
algorithms that are capable of fullfilling the above goals.


Future plans
------------

- Automatized parallelization of simulator runs on computing clusters using
  `dask`.
- Automatized detection of optimal network structures by hyper-parameter
  optimization on GPU clusters.
- Automatized consistency and coverage checks of estimated posteriors.
- Convenient handling of cached simulations, e.g. changes of priors and
  parameter ranges.
- Flexible handling of high-dimensional constrained priors.

.. note::
   The future development of *swyft* happens in collaboration with engineers
   from the Dutch `eScience center <https://www.esciencecenter.nl/>`_ and
   `SURFsara <https://surf.nl>`_, as part of the eTEC-BIG grant `Dark
   Generators <https://www.esciencecenter.nl/projects/darkgenerators/>`_. Any
   feedback will be useful to shape the development of this tool. Stay tuned!
Quickstart with *swyft*
=======================
Run this example in a colab notebook here_.

..  _here: https://colab.research.google.com/github/undark-lab/swyft/blob/master/notebooks/Quickstart.ipynb

.. As a quick example, the following code defines a simple "simulator" and noise model and performs inference given a particular draw.
.. ::
..     import numpy as np
..     import pylab as plt
..     import swyft
..     import torch
..     from scipy import stats

..     DEVICE = 'cuda' #your gpu, or 'cpu' if a gpu is not available

..     #a simple simulator
..     def model(params):
..         a = params['a']
..         b = params['b']
..         x=np.array([a,2*(b-a)])
..         return dict(mu=x)

..     #a simple noise model
..     def noise(obs, params, noise = 0.01):
..         x = obs['mu']
..         n = np.random.randn(*x.shape)*noise
..         return dict(x=x + n)

..     #choose the "true" parameters for an inference problem
..     par0 = dict(a=0.55, b=0.45)
..     obs0 = model(par0) # using Asimov data

..     #give priors for model parameters
..     prior = swyft.Prior({"a": ["uniform", 0., 1.], "b": ["uniform",  0., 1.]})

..     #a simple inference
..     s = swyft.NestedRatios(model, prior, noise = noise, obs = obs0, device = DEVICE)
..     #train!
..     s.run(Ninit = 500)

.. The last line, which trains networks that estimate the 1-dimensional marginal posteriors, will output something like:
.. ::
..     Simulate:  14%|█▎        | 67/495 [00:00<00:00, 667.16it/s]

..     NRE ROUND 0

..     Simulate: 100%|██████████| 495/495 [00:00<00:00, 644.85it/s]

..     NRE ROUND 1

..     Simulate: 100%|██████████| 517/517 [00:00<00:00, 643.51it/s]

..     NRE ROUND 2

..     Simulate: 100%|██████████| 498/498 [00:00<00:00, 713.97it/s]

..     NRE ROUND 3

..     Simulate: 100%|██████████| 820/820 [00:01<00:00, 647.67it/s]

..     NRE ROUND 4

..     Simulate: 100%|██████████| 1598/1598 [00:02<00:00, 653.44it/s]

..     NRE ROUND 5

..     Simulate: 100%|██████████| 2745/2745 [00:04<00:00, 672.84it/s]

..     NRE ROUND 6

..     Simulate: 100%|██████████| 5027/5027 [00:07<00:00, 704.09it/s]

..     NRE ROUND 7
..     --> Posterior volume is converged. <--


.. This "zooms in" to the relevant region of parameter space. The resulting marginal posteriors can be plotted:
.. ::
..     #train 2d marginals
..     post = s.gen_2d_marginals(N = 15000)
..     #generate samples at which to evaluate posteriors
..     samples = post(obs0, 1000000);
..     #plot estimated posteriors
..     swyft.corner(samples, ["a", "b"], color='k', figsize = (15,15), truth=par0)

.. .. image:: images/quickstart-2d.png
..    :width: 600

For details on tweaking *swyft*, see the tutorial as a notebook on github_ or colab_.

.. _github: https://github.com/undark-lab/swyft/blob/master/notebooks/Tutorial.ipynb
.. _colab: https://colab.research.google.com/github/undark-lab/swyft/blob/master/notebooks/Tutorial.ipynb
Theoretical concepts
====================

Introduction
============

Parametric stochastic simulators are ubiquitous in the physical sciences
:raw-latex:`\cite{Banik_2018, Bartels_2016, Rodr_guez_Puebla_2016}`.
However, performing parameter inference based on simulator runs using
Markov chain Monte Carlo is inconvenient or even impossible if the model
parameter space is large or the likelihood function is intractable. This
problem is addressed by so-called likelihood-free inference
:raw-latex:`\cite{sisson2018handbook}` or simulation-based inference
:raw-latex:`\cite{Cranmer2020}` techniques. Deep learning based
likelihood-free inference algorithms were organized into a taxonomy in
Ref. :raw-latex:`\cite{Durkan2020}`, where methods that estimated
likelihood ratios in a series of rounds were denoted Sequential Ratio
Estimation (SRE) :raw-latex:`\cite{Hermans2019}`. Our presented method
is closely related.

We propose *Nested Ratio Estimation* (NRE), which approximates the
likelihood-to-evidence ratio in a sequence of rounds. Loosely inspired
by the contour sorting method of nested sampling
:raw-latex:`\cite{Skilling2006, Feroz2008, Handley2015}`, the scheme
alternates between sampling from a constrained prior and estimating
likelihood-to-evidence ratios. It allows for efficient estimation of any
marginal posteriors of interest. Furthermore, we propose an algorithm
that we call *iP3 sample caching*, which facilitates simulator
efficiency by automatizing the reuse of previous simulator runs through
resampling of cached simulations.

The primary use case for these algorithms is the calculation of
arbitrary, low-dimensional marginal posteriors, typically in one or two
dimensions. In physics and astronomy, such marginals serve as the basis
for scientific conclusions by constraining individual model parameters
within uncertainty bounds. We implement a multi-target training regime
where all marginal posteriors of interest can be learned simultaneously.
We find that learning is simplified when one calculates each marginal
distribution directly rather than computing the full joint posterior and
marginalizing numerically. Furthermore, the method facilitates
effortless marginalization over arbitrary numbers of nuisance
parameters, increasing its utility in high-dimensional parameter
regimes–even to simulators with a tractable, yet high-dimensional,
likelihood :raw-latex:`\cite{lensing}`.

Nested Ratio Estimation (NRE).
==============================

We operate in the context of simulation-based inference where our
simulator :math:`\mathbf{g}` is a nonlinear function mapping a vector of
parameters
:math:`\boldsymbol{\theta}= (\theta_{1}, \dots, \theta_{d}) \in \mathbb{R}^{d}`
and a stochastic latent state :math:`\mathbf{z}` to an observation
:math:`\mathbf{x}= \mathbf{g}(\boldsymbol{\theta}, \mathbf{z})`. The
likelihood function is therefore
:math:`p(\mathbf{x}\vert \boldsymbol{\theta}) = \int \delta(\mathbf{x}- \mathbf{g}(\boldsymbol{\theta}, \mathbf{z})) \, p(\mathbf{z}\vert \boldsymbol{\theta}) \, d\mathbf{z}`,
with :math:`\delta(\cdot)` denoting the Dirac delta. Consider a
factorizable prior
:math:`p(\boldsymbol{\theta}) = p(\theta_{1}) \cdots p(\theta_{d})` over
the parameters, the joint posterior is given via Bayes’ rule as
:math:`p(\boldsymbol{\theta}|\mathbf{x}) = p(\mathbf{x}|\boldsymbol{\theta})p(\boldsymbol{\theta})/p(\mathbf{x})`,
where :math:`p(\mathbf{x})` is the evidence.

Our goal is to compute the marginal posterior,
:math:`p(\boldsymbol{\vartheta}\vert \mathbf{x})`, where
:math:`\boldsymbol{\vartheta}` are the parameters of interest. We will
denote all other parameters by :math:`\boldsymbol{\eta}`, such that
:math:`\boldsymbol{\theta}= (\boldsymbol{\vartheta}, \boldsymbol{\eta})`.
The marginal posterior is obtained from the joint distribution
:math:`p(\boldsymbol{\vartheta}, \boldsymbol{\eta}|\mathbf{x}) \equiv p(\boldsymbol{\theta}|\mathbf{x})`
by integrating over all components of :math:`\boldsymbol{\eta}`,

.. math::

   \label{eqn:post}
   p(\boldsymbol{\vartheta}\vert \mathbf{x})  \equiv \int p(\boldsymbol{\vartheta}, \boldsymbol{\eta}| \mathbf{x}) d\boldsymbol{\eta}
   = \int \frac{p(\mathbf{x}| \boldsymbol{\vartheta}, \boldsymbol{\eta})}{p(\mathbf{x})}
   p(\boldsymbol{\theta})
   %\prod_{j \notin \texttt{idx}} d\theta_{j}
   d\boldsymbol{\eta}
   = \frac{p(\mathbf{x}|\boldsymbol{\vartheta})}{p(\mathbf{x})}p(\boldsymbol{\vartheta})\;,

where we used Bayes’ rule and defined the marginal likelihood
:math:`p(\mathbf{x}|\boldsymbol{\vartheta})` in the last step.

Just like in SRE, we focus on a specific observation of interest,
:math:`\mathbf{x}_0`. Only parameter values :math:`\boldsymbol{\theta}`
that could have plausibly generated observation :math:`\mathbf{x}_0`
will significantly contribute to the integrals in
Eq. `[eqn:post] <#eqn:post>`__. For implausible values the likelihood
:math:`p(\mathbf{x}_0|\boldsymbol{\theta})` will be negligible. We
denote priors that are suitably constrained to plausible parameter
values by :math:`\tilde{p}(\theta_1, \dots, \theta_d)`. Similarly,
:math:`\tilde{\square}` indicates quantities :math:`\square` that are
calculated using the constrained prior. Therefore, using a judiciously
chosen constrained prior, accurately approximates the marginal posterior
in place of our true prior beliefs,

.. math::

   p(\boldsymbol{\vartheta}| \mathbf{x}_0) =
   \frac{p(\mathbf{x}_0|\boldsymbol{\vartheta})}{p(\mathbf{x}_0)} p(\boldsymbol{\vartheta}) \simeq
   \frac{\tilde{p}(\mathbf{x}_0|\boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x}_0)} \tilde{p}(\boldsymbol{\vartheta})\;.

The increased probability that constrained priors assign to the
plausible parameter region cancels when dividing by the constrained
evidence :math:`\tilde p(\mathbf{x})`. We define the marginal
likelihood-to-evidence ratio

.. math::

   \label{eqn:likelihood_ratio}
       \tilde{r}(\mathbf{x}, \boldsymbol{\vartheta})
       \equiv \frac{\tilde{p}(\mathbf{x}\vert \boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x})}
       = \frac{\tilde{p}(\mathbf{x}, \boldsymbol{\vartheta})}{\tilde{p}(\mathbf{x}) \tilde{p}(\boldsymbol{\vartheta})}
       = \frac{\tilde{p}(\boldsymbol{\vartheta}\vert\mathbf{x})}{\tilde{p}(\boldsymbol{\vartheta})}\;,

which is sufficient to evaluate the marginal posterior in
Eq. `[eqn:post] <#eqn:post>`__, and which we will now estimate. Under
the assumption of equal class population, it is known
:raw-latex:`\cite{Durkan2020, Cranmer2015}` that one can recover density
ratios using binary classification to distinguish between samples from
two distributions. Our binary classification problem is to distinguish
positive samples,
:math:`(\mathbf{x}, \boldsymbol{\vartheta}) \sim \tilde{p}(\mathbf{x}, \boldsymbol{\vartheta}) = p(\mathbf{x}\vert \boldsymbol{\vartheta}) \tilde{p}(\boldsymbol{\vartheta})`,
drawn jointly, and negative samples,
:math:`(\mathbf{x}, \boldsymbol{\vartheta}) \sim \tilde{p}(\mathbf{x}) \tilde{p}(\boldsymbol{\vartheta})`,
drawn marginally. The binary classifier
:math:`\sigma(f_{\phi}(\mathbf{x}, \boldsymbol{\vartheta}))` performs
optimally when
:math:`f_{\phi}(\mathbf{x}, \boldsymbol{\vartheta}) = \log \tilde{r}(\mathbf{x}, \boldsymbol{\vartheta})`,
where :math:`\sigma(\cdot)` is the sigmoid function and :math:`f_{\phi}`
is a neural network parameterized by :math:`\phi`. The associated binary
cross-entropy loss function used to train the ratio
:math:`\tilde{r}(\boldsymbol{\vartheta}, \mathbf{x}_0)` via stochastic
gradient descent is given by

.. math:: \ell = -\int \left[ \tilde{p}(\mathbf{x}|\boldsymbol{\vartheta})\tilde{p}(\boldsymbol{\vartheta}) \ln \sigma(f_\phi(\mathbf{x}, \boldsymbol{\vartheta})) + \tilde{p}(\mathbf{x})\tilde{p}(\boldsymbol{\vartheta}) \ln \sigma(-f_\phi(\mathbf{x},\boldsymbol{\vartheta})) \right] d\mathbf{x}\, d\boldsymbol{\vartheta}\;.

We propose to iteratively improve marginal posterior estimates in
:math:`R` rounds by employing posterior estimates from previous rounds
to define constrained priors. In each round :math:`r`, we estimate *all*
1-dim marginal posteriors, using :math:`d` instances of the above
marginal likelihood-to-evidence ratio estimation in parallel by setting
:math:`\boldsymbol{\vartheta}= (\theta_i)` for :math:`i=1, \dots, d`. To
this end, we utilize the factorized constrained prior,
:math:`\tilde{p}_r(\theta) = \tilde{p}_r(\theta_1)\cdots\tilde{p}_r(\theta_d)`,
which is defined recursively by a cutoff criterion,

.. math::

   \tilde{p}_{r}(\theta_{i})
       \propto
       p(\theta_{i}) \Theta_{H} \left[ \frac{\tilde{r}_{r-1}(\theta_{i}, \mathbf{x})}{\max_{\theta_{i}} \tilde{r}_{r-1}(\theta_{i}, \mathbf{x})} - \epsilon \right],
       \label{eqn:it}

where :math:`\Theta_{H}` denotes the Heaviside step function and
:math:`\epsilon` denotes the minimum likelihood-ratio which passes
through the threshold. We use
:math:`\tilde{p}_1(\boldsymbol{\theta}) = p(\boldsymbol{\theta})` as an
initial prior in the iterative scheme.

In every round, each 1-dim posterior approximates a marginalization of
the same underlying constrained posterior, allowing us to effectively
reuse training data and train efficiently in a multi-target regime. The
inference network is therefore divided into a featurizer
:math:`\mathbf{F}(\mathbf{x})` with shared parameters and a set of
:math:`d` independent Multi-layer Perceptons
:math:`\{\textrm{MLP}_i(\cdot, \cdot)\}_{i=1}^{d}` which estimate
individual 1-dim marginal posteriors and do not share parameters, such
that
:math:`f_{\phi}(\mathbf{x}, \theta_i) = \textrm{MLP}_i(\mathbf{F}(\mathbf{x}), \theta_i)`.

This technique is valid as long as the excluded prior regions do not
significantly affect the integrals in Eq. `[eqn:post] <#eqn:post>`__.
For uncorrelated parameters, a sufficient criterion is that the impact
on the marginal posteriors is small, which we guarantee through the
iteration criterion Eq. `[eqn:it] <#eqn:it>`__. In the case of a very
large number of strongly correlated parameters the algorithm can
inadvertently cut away tails of the marginal posteriors. Decreasing
:math:`\epsilon` mitigates this effect. Discussion is left for future
study :raw-latex:`\cite{swyft_future}`.

With this design, the posteriors from the final round can be used to
approximate the true 1-dim marginal posteriors,
:math:`\tilde{p}_{R}(\theta_i \vert \mathbf{x}_{0}) \approx p(\theta_i\vert \mathbf{x}_{0})`,
while previous rounds were used to iteratively focus on relevant parts
of the parameter space. The key result and value of NRE lies in the
utility of our constrained prior from round :math:`R`. The final
constrainted prior, along with previously generated and cached samples,
allows for estimation of *any* higher dimensional marginal posterior
:math:`\tilde{p}_R(\boldsymbol{\vartheta}|\mathbf{x}_0) \approx p(\boldsymbol{\vartheta}|\mathbf{x}_0)`
of interest by doing likelihood-to-evidence ratio estimation, often
without further simulation.

Inhomogeneous Poisson Point Process (iP3) Sample Caching.
=========================================================

Simulating
:math:`(\mathbf{x}, \boldsymbol{\theta})\sim p(\mathbf{x}|\boldsymbol{\theta})p(\boldsymbol{\theta})`
can be extremely expensive. We develop a scheme to systematically reuse
appropriate subsets of previous simulator runs. Our method samples
:math:`N\sim \text{Pois}(\hat N)` parameter vectors from an arbitrary
distribution :math:`p(\boldsymbol{\theta})`, where :math:`\hat N` is the
expected number of samples. Taking :math:`N` samples from
:math:`p(\boldsymbol{\theta})` is equivalent to drawing a single sample
:math:`\Theta \equiv \{\boldsymbol{\theta}^{(n)}\}_{n=1}^{N}` from an
inhomogenous Poisson point process (PPP) with intensity function
:math:`\lambda_{r}(\boldsymbol{\theta}) = \hat{N} p(\boldsymbol{\theta})`.
In this context, :math:`\Theta` is known as a set of *points*. This
formulation provides convenient mathematical properties
:raw-latex:`\cite{ppp}`, at the low price of introducing variance in the
number of samples drawn. The precise number of samples doesn’t matter as
long as :math:`N \approx \hat{N}`, which is true in our regime of order
:math:`\geq 1000`.

We will need two properties of PPPs. *Superposition:* Given two
independent PPPs with intensity functions
:math:`\lambda_{1}(\boldsymbol{\theta})` and
:math:`\lambda_{2}(\boldsymbol{\theta})`, the sum yields another PPP
with intensity function
:math:`\lambda(\boldsymbol{\theta}) = \lambda_{1}(\boldsymbol{\theta}) + \lambda_{2}(\boldsymbol{\theta})`.
The union of two sets of points :math:`\Theta = \Theta_1 \cup \Theta_2`
from the individual PPPs is equivalent to a single set of points from
the combined PPP. *Thinning:* Consider a PPP with intensity function
:math:`\lambda(\boldsymbol{\theta})`, and an arbitrary function
:math:`q(\boldsymbol{\theta}): \mathbb{R}^{d} \to [0, 1]`. If we are
interested in drawing from a PPP with intensity function
:math:`\lambda_{q}(\boldsymbol{\theta}) = q(\boldsymbol{\theta}) \lambda(\boldsymbol{\theta})`,
we can achieve this by drawing a set of points :math:`\Theta`
distributed like :math:`\lambda(\boldsymbol{\theta})` and then rejecting
individual points :math:`\boldsymbol{\theta}^{(n)}` with probability
:math:`1 - q(\boldsymbol{\theta}^{(n)})`.

The parameter cache is defined by a set of points :math:`\Theta_{sc}`
drawn from a PPP with intensity function
:math:`\lambda_{sc}(\boldsymbol{\theta})`. For every point
:math:`\boldsymbol{\theta}\in\Theta_{sc}`, a corresponding observation
:math:`\mathbf{x}` is stored in an observation cache
:math:`\mathcal{X}_{sc}`. The iP3 cache sampling algorithm that is
responsible for maintaining the caches and sampling from a PPP with
target intensity function
:math:`\lambda_t(\boldsymbol{\theta}) = \hat{N} p(\boldsymbol{\theta})`
is written out in the supplementary material. It is summarized in two
steps: First, consider all points
:math:`\boldsymbol{\theta}\in \Theta_{sc}` from the cache and accept
them with probability
:math:`\min(1, \lambda_t(\boldsymbol{\theta})/\lambda_{sc}(\boldsymbol{\theta}))`.
The thinning operation yields a sample :math:`\Theta_1` from a PPP with
intensity function
:math:`\lambda_1(\boldsymbol{\theta}) = \min(\lambda_t(\boldsymbol{\theta}), \lambda_{sc}(\boldsymbol{\theta}))`.
Second, draw a new set of points :math:`\Theta_p` from
:math:`\lambda_t(\boldsymbol{\theta})`, and accept each
:math:`\boldsymbol{\theta}\in\Theta_p` with probability
:math:`\max(0, 1-\lambda_{sc}(\boldsymbol{\theta})/\lambda_t(\boldsymbol{\theta}))`.
This yields a sample :math:`\Theta_2` from a PPP with intensity function
:math:`\lambda_2(\boldsymbol{\theta}) = \max(0, \lambda_t(\boldsymbol{\theta}) - \lambda_{sc}(\boldsymbol{\theta}))`.
Thanks to superposition, the union
:math:`\Theta_1 \cup \Theta_2 = \Theta_t` yields a sample from the PPP
with intensity function :math:`\lambda_t(\boldsymbol{\theta})`–the
sample we were looking for. We only need to run simulations on points
from :math:`\Theta_1`. Points in :math:`\Theta_2` already have
corresponding observations in :math:`\mathcal{X}_{sc}` which we can
reuse. Finally, the new parameters are appended to the set of points in
the parameter cache, :math:`\Theta_{sc} \to \Theta_{sc} \cup \Theta_2`.
Similar for :math:`\mathcal{X}_{sc}`. On the basis of the superposition
principle, the intensity function of the :math:`\Theta_{sc}` cache is
updated
:math:`\lambda_{sc}(\boldsymbol{\theta}) \to \max(\lambda_{sc}(\boldsymbol{\theta}), \lambda_t(\boldsymbol{\theta}))`.

Storing and updating the parameter cache’s intensity function
:math:`\lambda_{sc}(\boldsymbol{\theta})` can pose challenges when it is
complex and high-dimensional. Our NRE implementation overcomes these
challenges by learning marginal 1-dim posteriors, guaranteeing that the
relevant target intensities always factorize,
:math:`\lambda_t(\boldsymbol{\theta}) = \lambda_t(\theta_1)\cdots \lambda_t(\theta_d)`.
Storage of and calculation with factorizable functions simplifies
matters.
What is *swyft*?
================

Marginal posterior estimation
-----------------------------

*swyft* directly estimates marginal posteriors for parameters of interest
:math:`\mathbf{z}`, given some observation :math:`\mathbf{x}`. These are
formally obtained by integrating over all remaining (nuisance) parameters
:math:`\boldsymbol{\eta}`,

.. math::
   p(\mathbf{z}|\mathbf{x}) = \frac{\int d\boldsymbol{\eta}\,
   p(\mathbf{x}|\mathbf{z}, \boldsymbol{\eta}) p(\boldsymbol{\eta}, \mathbf{z})}
   {p(\mathbf{x})}\;.

Here, :math:`p(\mathbf{x}|\mathbf{z}, \boldsymbol{\eta})` is an abritrary
forward model that includes both the physics and detector simulator,
:math:`p(\mathbf{z}, \boldsymbol{\eta})` is the joint prior,
and :math:`p(\mathbf{x})` is the Bayesian evidence.


Nuisance parameters — Yes please!
---------------------------------

*In the context of likelihood-based inference, nuisance parameters are an
integration problem.* Given the likelihood density
:math:`p(\mathbf{x}|\mathbf{z}, \boldsymbol{\eta})` for a particular
observation :math:`\mathbf{x}`, one attempts to solve the above integral over
:math:`\boldsymbol{\eta}`, e.g. through sampling based methods.  This becomes
increasingly challenging if the number of nuisance parameters grows.

*In the context of likelihood-free inference, nuisance parameters are noise.*
Posteriors are estimated based on a large number of training samples
:math:`\mathbf{x}, \mathbf{z}\sim p(\mathbf{x}|\mathbf{z},
\boldsymbol{\eta})p(\mathbf{z}, \boldsymbol{\eta})`, no matter the dimension
of the nuisance parameter space. For a given :math:`\mathbf{z}`, more nuisance
parameters just increase the variance of :math:`\mathbf{x}` (which oddly enough
can make the inference problem simpler rather than more difficult).

.. note::
   When performing likelihood-based modelling, it is sometimes necessary to
   simplify the model for computational purposes.
   *swyft* uses likelihood-free inference, which means that models can be as
   complex as they need to be to describe reality (more specifically, we use
   the effective AALR [1]).


Not a Markov Chain
------------------

*Likelihood-based techniques often use Markov chains*, which require a simulation
for every link in the chain. Due to the properties of Markov chains, it is not
possible to utilize those simulations again for further analysis.
That effort has been lost.

*Likelihood-free inference can be based on simulations that sample the
(constrained) prior*. Reusing these simulations is allowed, we don’t
have to worry about breaking the Markov chain.

.. note::
   *swyft* automatizes the re-use of simulator runs where appropriate, using a
   new resampling approach (iP3 sample caching [2]).


High precision
--------------

Likelihood-based techniques are highly precise by focusing simulator
runs on parameter space regions that are consistent with a particular
observation.

Likelihood-free inference techniques can be less precise when there are
too few simulations in parameter regions that matter most.

.. note::
   *swyft* uses a new nested sampling scheme to target parameter regions most
   relevant for a given observation. This allows similar precision to
   likelihood-based approaches, without the high number of simulator runs
   (nested ratio estimation, NRE [2]).


Where is the catch?
-------------------

*swyft* uses neural likelihood estimation. The package is supposed to work
out-of-the-box for simple low-dimensional data. However, tackling
complex and/or high-dimensional data (think of high-resolution images or
spectra, combination of multiple data sets) requires some basic skills
in writing neural networks using pytorch.

.. note::
   *swyft* SWYFT provides a simple gateway for spicing up your analysis
   with the power of neural network-based inference.


References
----------

[1] Joeri Hermans, Volodimir Begy, and Gilles Louppe. Likelihood-free mcmc
with amortized approximate ratio estimators. arXiv preprint arXiv:1903.04057, 2019.

[2] Benjamin Kurt Miller, Alex Cole, Gilles Louppe, and Christoph Weniger.
Simulation-efficient marginal posterior estimationwithswyft: stop wasting your freaking time.
arXiv preprint, 2020.
