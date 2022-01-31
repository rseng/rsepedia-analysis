Release Notes
=============

The list of changes for each statsmodels release can be found [here](https://www.statsmodels.org/devel/release/index.html). Full details are available in the [commit logs](https://github.com/statsmodels/statsmodels).
- [ ] closes #xxxx
- [ ] tests added / passed. 
- [ ] code/documentation is well formatted.  
- [ ] properly formatted commit message. See 
      [NumPy's guide](https://docs.scipy.org/doc/numpy-1.15.1/dev/gitwash/development_workflow.html#writing-the-commit-message). 

<details>


**Notes**:

* It is essential that you add a test when making code changes. Tests are not 
  needed for doc changes.
* When adding a new function, test values should usually be verified in another package (e.g., R/SAS/Stata).
* When fixing a bug, you must add a test that would produce the bug in main and
  then show that it is fixed with the new code.
* New code additions must be well formatted. Changes should pass flake8. If on Linux or OSX, you can
  verify you changes are well formatted by running 
  ```
  git diff upstream/main -u -- "*.py" | flake8 --diff --isolated
  ```
  assuming `flake8` is installed. This command is also available on Windows 
  using the Windows System for Linux once `flake8` is installed in the 
  local Linux environment. While passing this test is not required, it is good practice and it help 
  improve code quality in `statsmodels`.
* Docstring additions must render correctly, including escapes and LaTeX.

</details>
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

#### Describe the bug

[A clear and concise description of what the bug is. This should explain **why** the current behaviour is a problem and why the expected output is a better solution.]

#### Code Sample, a copy-pastable example if possible


```python
# Your code here that produces the bug
# This example should be self-contained, and so not rely on external data.
# It should run in a fresh ipython session, and so include all relevant imports.
```
<details>

**Note**: As you can see, there are many issues on our GitHub tracker, so it is very possible that your issue has been posted before. Please check first before submitting so that we do not have to handle and close duplicates.

**Note**: Please be sure you are using the latest released version of `statsmodels`, or a recent build of `main`. If your problem has been fixed in an unreleased version, you might be able to use `main` until a new release occurs. 

**Note**: If you are using a released version, have you verified that the bug exists in the main branch of this repository? It helps the limited resources if we know problems exist in the current main branch so that they do not need to check whether the code sample produces a bug in the next release.

</details>


If the issue has not been resolved, please file it in the issue tracker.

#### Expected Output

A clear and concise description of what you expected to happen.

#### Output of ``import statsmodels.api as sm; sm.show_versions()``

<details>

[paste the output of ``import statsmodels.api as sm; sm.show_versions()`` here below this line]

</details>
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

#### Is your feature request related to a problem? Please describe
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

#### Describe the solution you'd like
A clear and concise description of what you want to happen.

#### Describe alternatives you have considered
A clear and concise description of any alternative solutions or features you have considered.

#### Additional context
Add any other context about the feature request here.This directory holds files that were once part of statsmodels but
are no longer maintained.  They are retained here in order to have their
git histories readily available, but should *not* be considered usable.
# Continuous Integration Tools

These scripts are used to implement Continuous Integration on Travis, AppVeyor
and Azure.  They should not be removed without careful consideration of the 
consequences.
# Documentation Documentation

We use a combination of sphinx and Jupyter notebooks for the documentation.
Jupyter notebooks should be used for longer, self-contained examples demonstrating
a topic.
Sphinx is nice because we get the tables of contents and API documentation.

## Build Process

Building the docs requires a few additional dependencies. You can get most
of these with

```bash

   pip install -e .[docs]

```

From the root of the project.
Some of the examples rely on `rpy2` to execute R code from the notebooks.
It's not included in the setup requires since it's known to be difficult to
install.

To generate the HTML docs, run ``make html`` from the ``docs`` directory.
This executes a few distinct builds

1. datasets
2. notebooks
3. sphinx

# Notebook Builds

We're using `nbconvert` to execute the notebooks, and then convert them
to HTML. The conversion is handled by `statsmodels/tools/nbgenerate.py`.
The default python kernel (embedded in the notebook) is `python3`.
You need at least `nbconvert==4.2.0` to specify a non-default kernel,
which can be passed in the Makefile.
The code in this folder is based on the Federal Reserve Bank of New York code
found at https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting, which was
downloaded as of commit 19f365cab8269e3aac3faa11ad091d6e913c5c43. Only the
files from that repository which were required for generating the test results
are included here.

In additionm the following files from the original package have been modified
(use git diff against the above repository to see the changes)

- functions/dfm.m
- functions/update_nowcast.m

The following files are not a part of the original package:

- test_DFM_blocks.m
- test_DFM.m
- test_news_blocks.m
- test_news.m
- test_spec_blocks.xls
- test_spec.xls
Contributing guidelines
=======================

This page explains how you can contribute to the development of `statsmodels`
by submitting patches, statistical tests, new models, or examples.

`statsmodels` is developed on `Github <https://github.com/statsmodels/statsmodels>`_
using the `Git <https://git-scm.com/>`_ version control system.

Submitting a Bug Report
~~~~~~~~~~~~~~~~~~~~~~~

- Include a short, self-contained code snippet that reproduces the problem
- Specify the statsmodels version used. You can do this with ``sm.version.full_version``
- If the issue looks to involve other dependencies, also include the output of ``sm.show_versions()``

Making Changes to the Code
~~~~~~~~~~~~~~~~~~~~~~~~~~

For a pull request to be accepted, you must meet the below requirements. This greatly helps in keeping the job of maintaining and releasing the software a shared effort.

- **One branch. One feature.** Branches are cheap and github makes it easy to merge and delete branches with a few clicks. Avoid the temptation to lump in a bunch of unrelated changes when working on a feature, if possible. This helps us keep track of what has changed when preparing a release.
- Commit messages should be clear and concise. This means a subject line of less than 80 characters, and, if necessary, a blank line followed by a commit message body. We have an `informal commit format standard <https://www.statsmodels.org/devel/dev/maintainer_notes.html#commit-comments>`_ that we try to adhere to. You can see what this looks like in practice by ``git log --oneline -n 10``. If your commit references or closes a specific issue, you can close it by mentioning it in the `commit message <https://help.github.com/articles/closing-issues-via-commit-messages>`_.  (*For maintainers*: These suggestions go for Merge commit comments too. These are partially the record for release notes.)
- Code submissions must always include tests. See our `notes on testing <https://www.statsmodels.org/devel/dev/test_notes.html>`_.
- Each function, class, method, and attribute needs to be documented using docstrings. We conform to the `numpy docstring standard <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_.
- If you are adding new functionality, you need to add it to the documentation by editing (or creating) the appropriate file in ``docs/source``.
- Make sure your documentation changes parse correctly. Change into the top-level ``docs/`` directory and type::

   make clean
   make html

  Check that the build output does not have *any* warnings due to your changes.
- Finally, please add your changes to the release notes. Open the ``docs/source/release/versionX.X.rst`` file that has the version number of the next release and add your changes to the appropriate section.

Linting
~~~~~~~

Due to the way we have the CI builds set up, the linter will not do anything unless the environmental variable $LINT is set to a truthy value.

- On MacOS/Linux

    LINT=true ./lint.sh

- Dependencies: flake8, git

How to Submit a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So you want to submit a patch to `statsmodels` but are not too familiar with github? Here are the steps you need to take.

1. `Fork <https://help.github.com/articles/fork-a-repo>`_ the `statsmodels repository <https://github.com/statsmodels/statsmodels>`_ on Github.
2. `Create a new feature branch <https://git-scm.com/book/en/Git-Branching-Basic-Branching-and-Merging>`_. Each branch must be self-contained, with a single new feature or bugfix.
3. Make sure the test suite passes. This includes testing on Python 3. The easiest way to do this is to make a pull request and let the bot check for you. This can be slow, and if you are unsure about the fix or enhancement, it is best to run pytest locally.
4. Document your changes by editing the appropriate file in ``docs/source/``. If it is a big, new feature add a note and an example to the latest ``docs/source/release/versionX.X.rst`` file. See older versions for examples. If it's a minor change, it will be included automatically in our release notes.
5. Add an example. If it is a big, new feature please submit an example notebook by following `these instructions <https://www.statsmodels.org/devel/dev/examples.html>`_.
6. `Submit a pull request <https://help.github.com/articles/using-pull-requests>`_

Mailing List
~~~~~~~~~~~~

Conversations about development take place on the `statsmodels mailing list <https://groups.google.com/group/pystatsmodels?hl=en>`__.

Learn More
~~~~~~~~~~

The ``statsmodels`` documentation's `developer page <https://www.statsmodels.org/stable/dev/index.html>`_
offers much more detailed information about the process.

License
~~~~~~~

statsmodels is released under the
`Modified (3-clause) BSD license <https://www.opensource.org/licenses/BSD-3-Clause>`_.
|PyPI Version| |Conda Version| |License| |Azure CI Build Status|
|Codecov Coverage| |Coveralls Coverage| |PyPI downloads| |Conda downloads|

About statsmodels
=================

statsmodels is a Python package that provides a complement to scipy for
statistical computations including descriptive statistics and estimation
and inference for statistical models.


Documentation
=============

The documentation for the latest release is at

https://www.statsmodels.org/stable/

The documentation for the development version is at

https://www.statsmodels.org/dev/

Recent improvements are highlighted in the release notes

https://www.statsmodels.org/stable/release/

Backups of documentation are available at https://statsmodels.github.io/stable/
and https://statsmodels.github.io/dev/.


Main Features
=============

* Linear regression models:

  - Ordinary least squares
  - Generalized least squares
  - Weighted least squares
  - Least squares with autoregressive errors
  - Quantile regression
  - Recursive least squares

* Mixed Linear Model with mixed effects and variance components
* GLM: Generalized linear models with support for all of the one-parameter
  exponential family distributions
* Bayesian Mixed GLM for Binomial and Poisson
* GEE: Generalized Estimating Equations for one-way clustered or longitudinal data
* Discrete models:

  - Logit and Probit
  - Multinomial logit (MNLogit)
  - Poisson and Generalized Poisson regression
  - Negative Binomial regression
  - Zero-Inflated Count models

* RLM: Robust linear models with support for several M-estimators.
* Time Series Analysis: models for time series analysis

  - Complete StateSpace modeling framework

    - Seasonal ARIMA and ARIMAX models
    - VARMA and VARMAX models
    - Dynamic Factor models
    - Unobserved Component models

  - Markov switching models (MSAR), also known as Hidden Markov Models (HMM)
  - Univariate time series analysis: AR, ARIMA
  - Vector autoregressive models, VAR and structural VAR
  - Vector error correction model, VECM
  - exponential smoothing, Holt-Winters
  - Hypothesis tests for time series: unit root, cointegration and others
  - Descriptive statistics and process models for time series analysis

* Survival analysis:

  - Proportional hazards regression (Cox models)
  - Survivor function estimation (Kaplan-Meier)
  - Cumulative incidence function estimation

* Multivariate:

  - Principal Component Analysis with missing data
  - Factor Analysis with rotation
  - MANOVA
  - Canonical Correlation

* Nonparametric statistics: Univariate and multivariate kernel density estimators
* Datasets: Datasets used for examples and in testing
* Statistics: a wide range of statistical tests

  - diagnostics and specification tests
  - goodness-of-fit and normality tests
  - functions for multiple testing
  - various additional statistical tests

* Imputation with MICE, regression on order statistic and Gaussian imputation
* Mediation analysis
* Graphics includes plot functions for visual analysis of data and model results

* I/O

  - Tools for reading Stata .dta files, but pandas has a more recent version
  - Table output to ascii, latex, and html

* Miscellaneous models
* Sandbox: statsmodels contains a sandbox folder with code in various stages of
  development and testing which is not considered "production ready".  This covers
  among others

  - Generalized method of moments (GMM) estimators
  - Kernel regression
  - Various extensions to scipy.stats.distributions
  - Panel data models
  - Information theoretic measures

How to get it
=============
The main branch on GitHub is the most up to date code

https://www.github.com/statsmodels/statsmodels

Source download of release tags are available on GitHub

https://github.com/statsmodels/statsmodels/tags

Binaries and source distributions are available from PyPi

https://pypi.org/project/statsmodels/

Binaries can be installed in Anaconda

conda install statsmodels


Installing from sources
=======================

See INSTALL.txt for requirements or see the documentation

https://statsmodels.github.io/dev/install.html

Contributing
============
Contributions in any form are welcome, including:

* Documentation improvements
* Additional tests
* New features to existing models
* New models

https://www.statsmodels.org/stable/dev/test_notes

for instructions on installing statsmodels in *editable* mode.

License
=======

Modified BSD (3-clause)

Discussion and Development
==========================

Discussions take place on the mailing list

https://groups.google.com/group/pystatsmodels

and in the issue tracker. We are very interested in feedback
about usability and suggestions for improvements.

Bug Reports
===========

Bug reports can be submitted to the issue tracker at

https://github.com/statsmodels/statsmodels/issues

.. |Azure CI Build Status| image:: https://dev.azure.com/statsmodels/statsmodels-testing/_apis/build/status/statsmodels.statsmodels?branchName=main
   :target: https://dev.azure.com/statsmodels/statsmodels-testing/_build/latest?definitionId=1&branchName=main
.. |Codecov Coverage| image:: https://codecov.io/gh/statsmodels/statsmodels/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/statsmodels/statsmodels
.. |Coveralls Coverage| image:: https://coveralls.io/repos/github/statsmodels/statsmodels/badge.svg?branch=main
   :target: https://coveralls.io/github/statsmodels/statsmodels?branch=main
.. |PyPI downloads| image:: https://img.shields.io/pypi/dm/statsmodels?label=PyPI%20Downloads
   :alt: PyPI - Downloads
   :target: https://pypi.org/project/statsmodels/
.. |Conda downloads| image:: https://img.shields.io/conda/dn/conda-forge/statsmodels.svg?label=Conda%20downloads
   :target: https://anaconda.org/conda-forge/statsmodels/
.. |PyPI Version| image:: https://img.shields.io/pypi/v/statsmodels.svg
   :target: https://pypi.org/project/statsmodels/
.. |Conda Version| image:: https://anaconda.org/conda-forge/statsmodels/badges/version.svg
   :target: https://anaconda.org/conda-forge/statsmodels/
.. |License| image:: https://img.shields.io/pypi/l/statsmodels.svg
   :target: https://github.com/statsmodels/statsmodels/blob/main/LICENSE.txt
Developers
----------

This directory is only of interest to developers.  Many of the files are
substantially out of date and subject to clean-up.
.. currentmodule:: statsmodels.discrete.discrete_model


.. _discretemod:

Regression with Discrete Dependent Variable
===========================================

Regression models for limited and qualitative dependent variables. The module
currently allows the estimation of models with binary (Logit, Probit), nominal
(MNLogit), or count (Poisson, NegativeBinomial) data.

Starting with version 0.9, this also includes new count models, that are still
experimental in 0.9, NegativeBinomialP, GeneralizedPoisson and zero-inflated
models, ZeroInflatedPoisson, ZeroInflatedNegativeBinomialP and
ZeroInflatedGeneralizedPoisson.

See `Module Reference`_ for commands and arguments.

Examples
--------

.. ipython:: python
  :okwarning:

  # Load the data from Spector and Mazzeo (1980)
  import statsmodels.api as sm
  spector_data = sm.datasets.spector.load_pandas()
  spector_data.exog = sm.add_constant(spector_data.exog)

  # Logit Model
  logit_mod = sm.Logit(spector_data.endog, spector_data.exog)
  logit_res = logit_mod.fit()
  print(logit_res.summary())

Detailed examples can be found here:


* `Overview <examples/notebooks/generated/discrete_choice_overview.html>`__
* `Examples <examples/notebooks/generated/discrete_choice_example.html>`__

Technical Documentation
-----------------------

Currently all models are estimated by Maximum Likelihood and assume
independently and identically distributed errors.

All discrete regression models define the same methods and follow the same
structure, which is similar to the regression results but with some methods
specific to discrete models. Additionally some of them contain additional model
specific methods and attributes.


References
^^^^^^^^^^

General references for this class of models are::

    A.C. Cameron and P.K. Trivedi.  `Regression Analysis of Count Data`.
        Cambridge, 1998

    G.S. Madalla. `Limited-Dependent and Qualitative Variables in Econometrics`.
        Cambridge, 1983.

    W. Greene. `Econometric Analysis`. Prentice Hall, 5th. edition. 2003.

Module Reference
----------------

.. module:: statsmodels.discrete.discrete_model
   :synopsis: Models for discrete data

The specific model classes are:

.. autosummary::
   :toctree: generated/

   Logit
   Probit
   MNLogit
   Poisson
   NegativeBinomial
   NegativeBinomialP
   GeneralizedPoisson

.. currentmodule:: statsmodels.discrete.count_model
.. module:: statsmodels.discrete.count_model

.. autosummary::
   :toctree: generated/

   ZeroInflatedPoisson
   ZeroInflatedNegativeBinomialP
   ZeroInflatedGeneralizedPoisson

.. currentmodule:: statsmodels.discrete.conditional_models
.. module:: statsmodels.discrete.conditional_models

.. autosummary::
   :toctree: generated/

   ConditionalLogit
   ConditionalMNLogit
   ConditionalPoisson

The cumulative link model for an ordinal dependent variable is currently
in miscmodels as it subclasses GenericLikelihoodModel. This will change
in future versions.

.. currentmodule:: statsmodels.miscmodels.ordinal_model
.. module:: statsmodels.miscmodels.ordinal_model

.. autosummary::
   :toctree: generated/

   OrderedModel

The specific result classes are:

.. currentmodule:: statsmodels.discrete.discrete_model

.. autosummary::
   :toctree: generated/

   LogitResults
   ProbitResults
   CountResults
   MultinomialResults
   NegativeBinomialResults
   GeneralizedPoissonResults

.. currentmodule:: statsmodels.discrete.count_model

.. autosummary::
   :toctree: generated/

   ZeroInflatedPoissonResults
   ZeroInflatedNegativeBinomialResults
   ZeroInflatedGeneralizedPoissonResults
   
.. currentmodule:: statsmodels.miscmodels.ordinal_model

.. autosummary::
   :toctree: generated/

   OrderedResults


:class:`DiscreteModel` is a superclass of all discrete regression models. The
estimation results are returned as an instance of one of the subclasses of
:class:`DiscreteResults`. Each category of models, binary, count and
multinomial, have their own intermediate level of model and results classes.
This intermediate classes are mostly to facilitate the implementation of the
methods and attributes defined by :class:`DiscreteModel` and
:class:`DiscreteResults`.

.. currentmodule:: statsmodels.discrete.discrete_model

.. autosummary::
   :toctree: generated/

   DiscreteModel
   DiscreteResults
   BinaryModel
   BinaryResults
   CountModel
   MultinomialModel

.. currentmodule:: statsmodels.discrete.count_model

.. autosummary::
   :toctree: generated/

   GenericZeroInflated
:orphan:

.. module:: statsmodels.tsa.vector_ar.var_model
   :synopsis: Vector autoregressions

.. currentmodule:: statsmodels.tsa.vector_ar.var_model

.. _var:

Vector Autoregressions :mod:`tsa.vector_ar`
===========================================

:mod:`statsmodels.tsa.vector_ar` contains methods that are useful
for simultaneously modeling and analyzing multiple time series using
:ref:`Vector Autoregressions (VAR) <var>` and
:ref:`Vector Error Correction Models (VECM) <vecm>`.

.. _var_process:

VAR(p) processes
----------------

We are interested in modeling a :math:`T \times K` multivariate time series
:math:`Y`, where :math:`T` denotes the number of observations and :math:`K` the
number of variables. One way of estimating relationships between the time series
and their lagged values is the *vector autoregression process*:

.. math::

   Y_t = \nu + A_1 Y_{t-1} + \ldots + A_p Y_{t-p} + u_t

   u_t \sim {\sf Normal}(0, \Sigma_u)

where :math:`A_i` is a :math:`K \times K` coefficient matrix.

We follow in large part the methods and notation of `Lutkepohl (2005)
<https://www.springer.com/gb/book/9783540401728>`__,
which we will not develop here.

Model fitting
~~~~~~~~~~~~~

.. note::

    The classes referenced below are accessible via the
    :mod:`statsmodels.tsa.api` module.

To estimate a VAR model, one must first create the model using an `ndarray` of
homogeneous or structured dtype. When using a structured or record array, the
class will use the passed variable names. Otherwise they can be passed
explicitly:

.. ipython:: python
   :suppress:

   import pandas as pd
   pd.options.display.max_rows = 10
   import matplotlib
   import matplotlib.pyplot as plt
   matplotlib.style.use('ggplot')

.. ipython:: python
   :okwarning:

   # some example data
   import numpy as np
   import pandas
   import statsmodels.api as sm
   from statsmodels.tsa.api import VAR
   mdata = sm.datasets.macrodata.load_pandas().data

   # prepare the dates index
   dates = mdata[['year', 'quarter']].astype(int).astype(str)
   quarterly = dates["year"] + "Q" + dates["quarter"]
   from statsmodels.tsa.base.datetools import dates_from_str
   quarterly = dates_from_str(quarterly)

   mdata = mdata[['realgdp','realcons','realinv']]
   mdata.index = pandas.DatetimeIndex(quarterly)
   data = np.log(mdata).diff().dropna()

   # make a VAR model
   model = VAR(data)

.. note::

   The :class:`VAR` class assumes that the passed time series are
   stationary. Non-stationary or trending data can often be transformed to be
   stationary by first-differencing or some other method. For direct analysis of
   non-stationary time series, a standard stable VAR(p) model is not
   appropriate.

To actually do the estimation, call the `fit` method with the desired lag
order. Or you can have the model select a lag order based on a standard
information criterion (see below):

.. ipython:: python
   :okwarning:

   results = model.fit(2)
   results.summary()

Several ways to visualize the data using `matplotlib` are available.

Plotting input time series:

.. ipython:: python
   :okwarning:

   @savefig var_plot_input.png
   results.plot()


Plotting time series autocorrelation function:

.. ipython:: python

   @savefig var_plot_acorr.png
   results.plot_acorr()


Lag order selection
~~~~~~~~~~~~~~~~~~~

Choice of lag order can be a difficult problem. Standard analysis employs
likelihood test or information criteria-based order selection. We have
implemented the latter, accessible through the :class:`VAR` class:

.. ipython:: python

   model.select_order(15)

When calling the `fit` function, one can pass a maximum number of lags and the
order criterion to use for order selection:

.. ipython:: python

   results = model.fit(maxlags=15, ic='aic')

Forecasting
~~~~~~~~~~~

The linear predictor is the optimal h-step ahead forecast in terms of
mean-squared error:

.. math::

   y_t(h) = \nu + A_1 y_t(h − 1) + \cdots + A_p y_t(h − p)

We can use the `forecast` function to produce this forecast. Note that we have
to specify the "initial value" for the forecast:

.. ipython:: python

   lag_order = results.k_ar
   results.forecast(data.values[-lag_order:], 5)

The `forecast_interval` function will produce the above forecast along with
asymptotic standard errors. These can be visualized using the `plot_forecast`
function:

.. ipython:: python

   @savefig var_forecast.png
   results.plot_forecast(10)

Class Reference
~~~~~~~~~~~~~~~

.. module:: statsmodels.tsa.vector_ar
   :synopsis: Vector autoregressions and related tools

.. currentmodule:: statsmodels.tsa.vector_ar.var_model


.. autosummary::
   :toctree: generated/

   VAR
   VARProcess
   VARResults


Post-estimation Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several process properties and additional results after
estimation are available for vector autoregressive processes.

.. currentmodule:: statsmodels.tsa.vector_ar.var_model
.. autosummary::
   :toctree: generated/

   LagOrderResults

.. currentmodule:: statsmodels.tsa.vector_ar.hypothesis_test_results
.. autosummary::
   :toctree: generated/

   HypothesisTestResults
   NormalityTestResults
   WhitenessTestResults


Impulse Response Analysis
-------------------------

*Impulse responses* are of interest in econometric studies: they are the
estimated responses to a unit impulse in one of the variables. They are computed
in practice using the MA(:math:`\infty`) representation of the VAR(p) process:

.. math::

    Y_t = \mu + \sum_{i=0}^\infty \Phi_i u_{t-i}

We can perform an impulse response analysis by calling the `irf` function on a
`VARResults` object:

.. ipython:: python
   :okwarning:

   irf = results.irf(10)

These can be visualized using the `plot` function, in either orthogonalized or
non-orthogonalized form. Asymptotic standard errors are plotted by default at
the 95% significance level, which can be modified by the user.

.. note::

    Orthogonalization is done using the Cholesky decomposition of the estimated
    error covariance matrix :math:`\hat \Sigma_u` and hence interpretations may
    change depending on variable ordering.

.. ipython:: python

   @savefig var_irf.png
   irf.plot(orth=False)


Note the `plot` function is flexible and can plot only variables of interest if
so desired:

.. ipython:: python

   @savefig var_realgdp.png
   irf.plot(impulse='realgdp')

The cumulative effects :math:`\Psi_n = \sum_{i=0}^n \Phi_i` can be plotted with
the long run effects as follows:

.. ipython:: python

   @savefig var_irf_cum.png
   irf.plot_cum_effects(orth=False)


.. currentmodule:: statsmodels.tsa.vector_ar.irf
.. autosummary::
   :toctree: generated/

   IRAnalysis

Forecast Error Variance Decomposition (FEVD)
--------------------------------------------

Forecast errors of component j on k in an i-step ahead forecast can be
decomposed using the orthogonalized impulse responses :math:`\Theta_i`:

.. math::

    \omega_{jk, i} = \sum_{i=0}^{h-1} (e_j^\prime \Theta_i e_k)^2 / \mathrm{MSE}_j(h)

    \mathrm{MSE}_j(h) = \sum_{i=0}^{h-1} e_j^\prime \Phi_i \Sigma_u \Phi_i^\prime e_j

These are computed via the `fevd` function up through a total number of steps ahead:

.. ipython:: python

   fevd = results.fevd(5)
   fevd.summary()

They can also be visualized through the returned :class:`FEVD` object:

.. ipython:: python

   @savefig var_fevd.png
   results.fevd(20).plot()


.. currentmodule:: statsmodels.tsa.vector_ar.var_model
.. autosummary::
   :toctree: generated/

   FEVD

Statistical tests
-----------------

A number of different methods are provided to carry out hypothesis tests about
the model results and also the validity of the model assumptions (normality,
whiteness / "iid-ness" of errors, etc.).

Granger causality
~~~~~~~~~~~~~~~~~

One is often interested in whether a variable or group of variables is "causal"
for another variable, for some definition of "causal". In the context of VAR
models, one can say that a set of variables are Granger-causal within one of the
VAR equations. We will not detail the mathematics or definition of Granger
causality, but leave it to the reader. The :class:`VARResults` object has the
`test_causality` method for performing either a Wald (:math:`\chi^2`) test or an
F-test.

.. ipython:: python

   results.test_causality('realgdp', ['realinv', 'realcons'], kind='f')

Normality
~~~~~~~~~

As pointed out in the beginning of this document, the white noise component
:math:`u_t` is assumed to be normally distributed. While this assumption
is not required for parameter estimates to be consistent or asymptotically
normal, results are generally more reliable in finite samples when residuals
are Gaussian white noise. To test whether this assumption is consistent with
a data set, :class:`VARResults` offers the `test_normality` method.

.. ipython:: python

    results.test_normality()

Whiteness of residuals
~~~~~~~~~~~~~~~~~~~~~~

To test the whiteness of the estimation residuals (this means absence of
significant residual autocorrelations) one can use the `test_whiteness`
method of :class:`VARResults`.


.. currentmodule:: statsmodels.tsa.vector_ar.hypothesis_test_results
.. autosummary::
   :toctree: generated/

   HypothesisTestResults
   CausalityTestResults
   NormalityTestResults
   WhitenessTestResults

.. _svar:

Structural Vector Autoregressions
---------------------------------

There are a matching set of classes that handle some types of Structural VAR models.

.. module:: statsmodels.tsa.vector_ar.svar_model
   :synopsis: Structural vector autoregressions and related tools

.. currentmodule:: statsmodels.tsa.vector_ar.svar_model

.. autosummary::
   :toctree: generated/

   SVAR
   SVARProcess
   SVARResults

.. _vecm:

Vector Error Correction Models (VECM)
-------------------------------------

Vector Error Correction Models are used to study short-run deviations from
one or more permanent stochastic trends (unit roots). A VECM models the
difference of a vector of time series by imposing structure that is implied
by the assumed number of stochastic trends. :class:`VECM` is used to
specify and estimate these models.

A VECM(:math:`k_{ar}-1`) has the following form

.. math::

    \Delta y_t = \Pi y_{t-1} + \Gamma_1 \Delta y_{t-1} + \ldots
                   + \Gamma_{k_{ar}-1} \Delta y_{t-k_{ar}+1} + u_t

where

.. math::

    \Pi = \alpha \beta'

as described in chapter 7 of [1]_.

A VECM(:math:`k_{ar} - 1`) with deterministic terms has the form

.. math::

   \Delta y_t = \alpha \begin{pmatrix}\beta' & \eta'\end{pmatrix} \begin{pmatrix}y_{t-1} \\
                D^{co}_{t-1}\end{pmatrix} + \Gamma_1 \Delta y_{t-1} + \dots + \Gamma_{k_{ar}-1} \Delta y_{t-k_{ar}+1} + C D_t + u_t.

In :math:`D^{co}_{t-1}` we have the deterministic terms which are inside
the cointegration relation (or restricted to the cointegration relation).
:math:`\eta` is the corresponding estimator. To pass a deterministic term
inside the cointegration relation, we can use the `exog_coint` argument.
For the two special cases of an intercept and a linear trend there exists
a simpler way to declare these terms: we can pass ``"ci"`` and ``"li"``
respectively to the `deterministic` argument. So for an intercept inside
the cointegration relation we can either pass ``"ci"`` as `deterministic`
or `np.ones(len(data))` as `exog_coint` if `data` is passed as the
`endog` argument. This ensures that :math:`D_{t-1}^{co} = 1` for all
:math:`t`.

We can also use deterministic terms outside the cointegration relation.
These are defined in :math:`D_t` in the formula above with the
corresponding estimators in the matrix :math:`C`. We specify such terms by
passing them to the `exog` argument. For an intercept and/or linear trend
we again have the possibility to use `deterministic` alternatively. For
an intercept we pass ``"co"`` and for a linear trend we pass ``"lo"`` where
the `o` stands for `outside`.

The following table shows the five cases considered in [2]_. The last
column indicates which string to pass to the `deterministic` argument for
each of these cases.

====  ===============================  ===================================  =============
Case  Intercept                        Slope of the linear trend            `deterministic`
====  ===============================  ===================================  =============
I     0                                0                                    ``"nc"``
II    :math:`- \alpha \beta^T \mu`     0                                    ``"ci"``
III   :math:`\neq 0`                   0                                    ``"co"``
IV    :math:`\neq 0`                   :math:`- \alpha \beta^T \gamma`      ``"coli"``
V     :math:`\neq 0`                   :math:`\neq 0`                       ``"colo"``
====  ===============================  ===================================  =============

.. currentmodule:: statsmodels.tsa.vector_ar.vecm
.. autosummary::
   :toctree: generated/

   VECM
   coint_johansen
   JohansenTestResult
   select_order
   select_coint_rank
   VECMResults
   CointRankResults


References
----------
.. [1] Lütkepohl, H. 2005. *New Introduction to Multiple Time Series Analysis*. Springer.

.. [2] Johansen, S. 1995. *Likelihood-Based Inference in Cointegrated *
       *Vector Autoregressive Models*. Oxford University Press.
Pitfalls
========

This page lists issues which may arise while using statsmodels. These
can be the result of data-related or statistical problems, software design,
"non-standard" use of models, or edge cases.

statsmodels provides several warnings and helper functions for diagnostic
checking (see this `blog article
<http://jpktd.blogspot.ca/2012/01/anscombe-and-diagnostic-statistics.html>`_
for an example of misspecification checks in linear regression). The coverage
is of course not comprehensive, but more warnings and diagnostic functions will
be added over time.

While the underlying statistical problems are the same for all statistical
packages, software implementations differ in the way extreme or corner cases
are handled. Please report corner cases for which the models might not work, so
we can treat them appropriately.

Repeated calls to fit with different parameters
-----------------------------------------------

Result instances often need to access attributes from the corresponding model
instance. Fitting a model multiple times with different arguments can change
model attributes. This means that the result instance may no longer point to
the correct model attributes after the model has been re-fit.

It is therefore best practice to create separate model instances when we want
to fit a model using different fit function arguments.

For example, this works without problem because we are not keeping the results
instance for further use ::

  mod = AR(endog)
  aic = []
  for lag in range(1,11):
      res = mod.fit(maxlag=lag)
      aic.append(res.aic)


However, when we want to hold on to two different estimation results, then it
is recommended to create two separate model instances. ::

  mod1 = RLM(endog, exog)
  res1 = mod1.fit(scale_est='mad')
  mod2 = RLM(endog, exog)
  res2 = mod2.fit(scale_est=sm.robust.scale.HuberScale())


Unidentified Parameters
-----------------------

Rank deficient exog, perfect multicollinearity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Models based on linear models, GLS, RLM, GLM and similar, use a generalized
inverse. This means that:

+ Rank deficient matrices will not raise an error
+ Cases of almost perfect multicollinearity or ill-conditioned design matrices might produce numerically unstable results. Users need to manually check the rank or condition number of the matrix if this is not the desired behavior

Note: statsmodels currently fails on the NIST benchmark case for Filip if the
data is not rescaled, see `this blog <http://jpktd.blogspot.ca/2012/03/numerical-accuracy-in-linear-least.html>`_

Incomplete convergence in maximum likelihood estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, the maximum likelihood estimator might not exist, parameters
might be infinite or not unique (e.g. (quasi-)separation in models with binary
endogenous variable). Under the default settings, statsmodels will print
a warning if the optimization algorithm stops without reaching convergence.
However, it is important to know that the convergence criteria may sometimes
falsely indicate convergence (e.g. if the value of the objective function
converged but not the parameters). In general, a user needs to verify
convergence.

For binary Logit and Probit models, statsmodels raises an exception if perfect
prediction is detected. There is, however, no check for quasi-perfect
prediction.

Other Problems
--------------

Insufficient variation in the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible that there is insufficient variation in the data for small
datasets or for data with small groups in categorical variables. In these
cases, the results might not be identified or some hidden problems might occur.

The only currently known case is a perfect fit in robust linear model estimation.
For RLM, if residuals are equal to zero, then it does not cause an exception,
but having this perfect fit can produce NaNs in some results (scale=0 and 0/0
division) (issue #55).

True parameter outside domain of model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases domain restrictions for parameters in a model might be
inconsistent with the data. In those cases, estimation might stop at
parameter values close to the boundary of the parameter space, but can slso
fail with runtime errors or produce nans during optimization.

Two examples:

Estimating a negative binomial model when the data is not overdispersed
relative to Poisson, that is the true model has the same dispersion as Poisson
or is underdispersed, is inconsistent with the overdispersion assumption of
the Negative Binomial distribution. The loglikelihood and its derivatives, as
implemented in statsmodels, cannot be evaluated for dispersion parameter at
zero dispersion or in a positve neighborhood of zero dispersion.

Zero inflated models currently use either Logit or Probit as model for
inflation. This means that no inflation is at the boundary of the parameter
space and zero deflation is outside the valid parameter space.
When the data has no zero inflation or is zero deflated, then the zero
inflation model will often fail in the optimization, or end up close to the
no inflation boundary.
.. currentmodule:: statsmodels.regression.mixed_linear_model

.. _mixedlmmod:

Linear Mixed Effects Models
===========================

Linear Mixed Effects models are used for regression analyses involving
dependent data.  Such data arise when working with longitudinal and
other study designs in which multiple observations are made on each
subject.  Some specific linear mixed effects models are

* *Random intercepts models*, where all responses in a group are
  additively shifted by a value that is specific to the group.

* *Random slopes models*, where the responses in a group follow a
  (conditional) mean trajectory that is linear in the observed
  covariates, with the slopes (and possibly intercepts) varying by
  group.

* *Variance components models*, where the levels of one or more
  categorical covariates are associated with draws from distributions.
  These random terms additively determine the conditional mean of each
  observation based on its covariate values.

The statsmodels implementation of LME is primarily group-based,
meaning that random effects must be independently-realized for
responses in different groups.  There are two types of random effects
in our implementation of mixed models: (i) random coefficients
(possibly vectors) that have an unknown covariance matrix, and (ii)
random coefficients that are independent draws from a common
univariate distribution.  For both (i) and (ii), the random effects
influence the conditional mean of a group through their matrix/vector
product with a group-specific design matrix.

A simple example of random coefficients, as in (i) above, is:

.. math::

   Y_{ij} = \beta_0 + \beta_1X_{ij} + \gamma_{0i} + \gamma_{1i}X_{ij} + \epsilon_{ij}

Here, :math:`Y_{ij}` is the :math:`j^\rm{th}` measured response for subject
:math:`i`, and :math:`X_{ij}` is a covariate for this response.  The
"fixed effects parameters" :math:`\beta_0` and :math:`\beta_1` are
shared by all subjects, and the errors :math:`\epsilon_{ij}` are
independent of everything else, and identically distributed (with mean
zero).  The "random effects parameters" :math:`\gamma_{0i}` and
:math:`\gamma_{1i}` follow a bivariate distribution with mean zero,
described by three parameters: :math:`{\rm var}(\gamma_{0i})`,
:math:`{\rm var}(\gamma_{1i})`, and :math:`{\rm cov}(\gamma_{0i},
\gamma_{1i})`.  There is also a parameter for :math:`{\rm
var}(\epsilon_{ij})`.

A simple example of variance components, as in (ii) above, is:

.. math::

   Y_{ijk} = \beta_0 + \eta_{1i} + \eta_{2j} + \epsilon_{ijk}

Here, :math:`Y_{ijk}` is the :math:`k^\rm{th}` measured response under
conditions :math:`i, j`.  The only "mean structure parameter" is
:math:`\beta_0`.  The :math:`\eta_{1i}` are independent and
identically distributed with zero mean, and variance :math:`\tau_1^2`,
and the :math:`\eta_{2j}` are independent and identically distributed
with zero mean, and variance :math:`\tau_2^2`.

statsmodels MixedLM handles most non-crossed random effects models,
and some crossed models.  To include crossed random effects in a
model, it is necessary to treat the entire dataset as a single group.
The variance components arguments to the model can then be used to
define models with various combinations of crossed and non-crossed
random effects.

The statsmodels LME framework currently supports post-estimation
inference via Wald tests and confidence intervals on the coefficients,
profile likelihood analysis, likelihood ratio testing, and AIC.

Examples
--------

.. ipython:: python

  import statsmodels.api as sm
  import statsmodels.formula.api as smf

  data = sm.datasets.get_rdataset("dietox", "geepack").data

  md = smf.mixedlm("Weight ~ Time", data, groups=data["Pig"])
  mdf = md.fit()
  print(mdf.summary())

Detailed examples can be found here

* `Mixed LM <examples/notebooks/generated/mixed_lm_example.html>`__

There are some notebook examples on the Wiki:
`Wiki notebooks for MixedLM <https://github.com/statsmodels/statsmodels/wiki/Examples#linear-mixed-models>`_



Technical Documentation
-----------------------

The data are partitioned into disjoint groups.
The probability model for group :math:`i` is:

.. math::

    Y = X\beta + Z\gamma + Q_1\eta_1 + \cdots + Q_k\eta_k + \epsilon

where

* :math:`n_i` is the number of observations in group :math:`i`
* :math:`Y` is a :math:`n_i` dimensional response vector
* :math:`X` is a :math:`n_i * k_{fe}` dimensional matrix of fixed effects
  coefficients
* :math:`\beta` is a :math:`k_{fe}`-dimensional vector of fixed effects slopes
* :math:`Z` is a :math:`n_i * k_{re}` dimensional matrix of random effects
  coefficients
* :math:`\gamma` is a :math:`k_{re}`-dimensional random vector with mean 0
  and covariance matrix :math:`\Psi`; note that each group
  gets its own independent realization of gamma.
* :math:`Q_j` is a :math:`n_i \times q_j` dimensional design matrix for the
  :math:`j^\rm{th}` variance component.
* :math:`\eta_j` is a :math:`q_j`-dimensional random vector containing independent
  and identically distributed values with variance :math:`\tau_j^2`.
* :math:`\epsilon` is a :math:`n_i` dimensional vector of i.i.d normal
  errors with mean 0 and variance :math:`\sigma^2`; the :math:`\epsilon`
  values are independent both within and between groups

:math:`Y, X, \{Q_j\}` and :math:`Z` must be entirely observed.  :math:`\beta`,
:math:`\Psi`, and :math:`\sigma^2` are estimated using ML or REML estimation,
and :math:`\gamma`, :math:`\{\eta_j\}` and :math:`\epsilon` are
random so define the probability model.

The marginal mean structure is :math:`E[Y|X,Z] = X*\beta`.  If only
the marginal mean structure is of interest, GEE is a good alternative
to mixed models.

Notation:

* :math:`cov_{re}` is the random effects covariance matrix (referred
  to above as :math:`\Psi`) and :math:`scale` is the (scalar) error
  variance.  There is also a single estimated variance parameter
  :math:`\tau_j^2` for each variance component.  For a single group,
  the marginal covariance matrix of endog given exog is
  :math:`scale*I + Z * cov_{re} * Z`, where :math:`Z` is the design
  matrix for the random effects in one group.

References
^^^^^^^^^^

The primary reference for the implementation details is:

*   MJ Lindstrom, DM Bates (1988).  *Newton Raphson and EM algorithms for
    linear mixed effects models for repeated measures data*.  Journal of
    the American Statistical Association. Volume 83, Issue 404, pages 1014-1022.

See also this more recent document:

* http://econ.ucsb.edu/~doug/245a/Papers/Mixed%20Effects%20Implement.pdf

All the likelihood, gradient, and Hessian calculations closely follow
Lindstrom and Bates.

The following two documents are written more from the perspective of
users:

* https://r-forge.r-project.org/scm/viewvc.php/*checkout*/www/lMMwR/lrgprt.pdf?revision=949&root=lme4&pathrev=1781

* http://lme4.r-forge.r-project.org/slides/2009-07-07-Rennes/3Longitudinal-4.pdf

.. Class hierarchy: TODO

   General references for this class of models are

Module Reference
----------------

.. module:: statsmodels.regression.mixed_linear_model
   :synopsis: Mixed Linear Models


The model class is:

.. autosummary::
   :toctree: generated/

   MixedLM

The result class is:

.. autosummary::
   :toctree: generated/

   MixedLMResults



.. module:: statsmodels.miscmodels
.. currentmodule:: statsmodels.miscmodels


.. _miscmodels:


Other Models :mod:`miscmodels`
==============================

:mod:`statsmodels.miscmodels` contains model classes and that do not yet fit into
any other category, or are basic implementations that are not yet polished and will most
likely still change. Some of these models were written as examples for the generic
maximum likelihood framework, and there will be others that might be based on general
method of moments.

The models in this category have been checked for basic cases, but might be more exposed
to numerical problems than the complete implementation. For example, count.Poisson has
been added using only the generic maximum likelihood framework, the standard errors
are based on the numerical evaluation of the Hessian, while discretemod.Poisson uses
analytical Gradients and Hessian and will be more precise, especially in cases when there
is strong multicollinearity.
On the other hand, by subclassing GenericLikelihoodModel, it is easy to add new models,
another example can be seen in the zero inflated Poisson model, miscmodels.count.


Count Models :mod:`count`
--------------------------

.. module:: statsmodels.miscmodels.count
.. currentmodule:: statsmodels.miscmodels.count

.. autosummary::
   :toctree: generated/

   PoissonGMLE
   PoissonOffsetGMLE
   PoissonZiGMLE

Linear Model with t-distributed errors
--------------------------------------

This is a class that shows that a new model can be defined by only specifying the
method for the loglikelihood. All result statistics are inherited from the generic
likelihood model and result classes. The results have been checked against R for a
simple case.

.. module:: statsmodels.miscmodels.tmodel
.. currentmodule:: statsmodels.miscmodels.tmodel

.. autosummary::
   :toctree: generated/

   TLinearModel





.. module:: statsmodels.othermod
.. currentmodule:: statsmodels.othermod


.. _othermod:


Other Models :mod:`othermod`
==============================

:mod:`statsmodels.othermod` contains model classes that do not fit into
any other category, for example models for a response variable ``endog`` that
has support on the unit interval or is positive or non-negative.

:mod:`statsmodels.othermod` contains models that are, or will be fully developed
in contrast to :mod:`statsmodels.miscmodels` which contains mainly examples
for the use of the generic likelihood model setup.

Status is experimental. The api and implementation will need to adjust as we
support more types of models, for example models with multiple exog and
multiple link functions.


Interval Models :mod:`betareg`
------------------------------

Models for continuous dependent variables that are in the unit interval such
as fractions. These Models are estimated by full Maximum Likelihood. 
Dependent variables on the unit interval can also be estimate by 
Quasi Maximum Likelihood using models for binary endog, such as Logit and
GLM-Binomial. (The implementation of discrete.Probit assumes binary endog and
cannot estimate a QMLE for continuous dependent variable.)

.. module:: statsmodels.othermod.betareg
.. currentmodule:: statsmodels.othermod.betareg

.. autosummary::
   :toctree: generated/

   BetaModel
   BetaResults
.. module:: statsmodels.tsa
   :synopsis: Time-series analysis

.. currentmodule:: statsmodels.tsa


.. _tsa:


Time Series analysis :mod:`tsa`
===============================

:mod:`statsmodels.tsa` contains model classes and functions that are useful
for time series analysis. Basic models include univariate autoregressive models (AR),
vector autoregressive models (VAR) and univariate autoregressive moving average models
(ARMA). Non-linear models include Markov switching dynamic regression and
autoregression. It also includes descriptive statistics for time series, for example autocorrelation, partial
autocorrelation function and periodogram, as well as the corresponding theoretical properties
of ARMA or related processes. It also includes methods to work with autoregressive and
moving average lag-polynomials.
Additionally, related statistical tests and some useful helper functions are available.

Estimation is either done by exact or conditional Maximum Likelihood or conditional
least-squares, either using Kalman Filter or direct filters.

Currently, functions and classes have to be imported from the corresponding module, but
the main classes will be made available in the statsmodels.tsa namespace. The module
structure is within statsmodels.tsa is

- stattools : empirical properties and tests, acf, pacf, granger-causality,
  adf unit root test, kpss test, bds test, ljung-box test and others.
- ar_model : univariate autoregressive process, estimation with conditional
  and exact maximum likelihood and conditional least-squares
- arima.model : univariate ARIMA process, estimation with alternative methods
- statespace : Comprehensive statespace model specification and estimation. See
  the :ref:`statespace documentation <statespace>`.
- vector_ar, var : vector autoregressive process (VAR) and vector error correction
  models, estimation, impulse response analysis, forecast error variance decompositions,
  and data visualization tools. See the :ref:`vector_ar documentation <var>`.
- arma_process : properties of arma processes with given parameters, this
  includes tools to convert between ARMA, MA and AR representation as well as
  acf, pacf, spectral density, impulse response function and similar
- sandbox.tsa.fftarma : similar to arma_process but working in frequency domain
- tsatools : additional helper functions, to create arrays of lagged variables,
  construct regressors for trend, detrend and similar.
- filters : helper function for filtering time series
- regime_switching : Markov switching dynamic regression and autoregression models

Some additional functions that are also useful for time series analysis are in
other parts of statsmodels, for example additional statistical tests.

Some related functions are also available in matplotlib, nitime, and
scikits.talkbox. Those functions are designed more for the use in signal
processing where longer time series are available and work more often in the
frequency domain.


.. currentmodule:: statsmodels.tsa


Descriptive Statistics and Tests
""""""""""""""""""""""""""""""""

.. autosummary::
   :toctree: generated/

   stattools.acovf
   stattools.acf
   stattools.pacf
   stattools.pacf_yw
   stattools.pacf_ols
   stattools.pacf_burg
   stattools.ccovf
   stattools.ccf
   stattools.adfuller
   stattools.kpss
   stattools.range_unit_root_test
   stattools.zivot_andrews
   stattools.coint
   stattools.bds
   stattools.q_stat
   stattools.breakvar_heteroskedasticity_test
   stattools.grangercausalitytests
   stattools.levinson_durbin
   stattools.innovations_algo
   stattools.innovations_filter
   stattools.levinson_durbin_pacf
   stattools.arma_order_select_ic
   x13.x13_arima_select_order
   x13.x13_arima_analysis

Estimation
""""""""""

The following are the main estimation classes, which can be accessed through
statsmodels.tsa.api and their result classes

Univariate Autoregressive Processes (AR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The basic autoregressive model in Statsmodels is:

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   ar_model.AutoReg
   ar_model.AutoRegResults
   ar_model.ar_select_order

The `ar_model.AutoReg` model estimates parameters using conditional MLE (OLS),
and supports exogenous regressors (an AR-X model) and seasonal effects.

AR-X and related models can also be fitted with the `arima.ARIMA` class and the
`SARIMAX` class (using full MLE via the Kalman Filter).


Autoregressive Moving-Average Processes (ARMA) and Kalman Filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Basic ARIMA model and results classes are as follows:

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   arima.model.ARIMA
   arima.model.ARIMAResults

This model allows estimating parameters by various methods (including
conditional MLE via the Hannan-Rissanen method and full MLE via the Kalman
filter). It is a special case of the `SARIMAX` model, and it includes a large
number of inherited features from the :ref:`state space <statespace>` models
(including prediction / forecasting, residual diagnostics, simulation and
impulse responses, etc.).

Exponential Smoothing
~~~~~~~~~~~~~~~~~~~~~

Linear and non-linear exponential smoothing models are available:

.. currentmodule:: statsmodels.tsa.holtwinters

.. autosummary::
   :toctree: generated/

   ExponentialSmoothing
   SimpleExpSmoothing
   Holt
   HoltWintersResults

Separately, linear and non-linear exponential smoothing models have also been
implemented based on the "innovations" state space approach. In addition to the
usual support for parameter fitting, in-sample prediction, and out-of-sample
forecasting, these models also support prediction intervals, simulation, and
more.

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   exponential_smoothing.ets.ETSModel
   exponential_smoothing.ets.ETSResults

Finally, linear exponential smoothing models have also been separately
implemented as a special case of the general state space framework (this is
separate from the "innovations" state space approach described above). Although
this approach does not allow for the non-linear (multiplicative) exponential
smoothing models, it includes all features of :ref:`state space <statespace>`
models (including prediction / forecasting, residual diagnostics, simulation
and impulse responses, etc.).

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   statespace.exponential_smoothing.ExponentialSmoothing
   statespace.exponential_smoothing.ExponentialSmoothingResults

ARMA Process
""""""""""""

The following are tools to work with the theoretical properties of an ARMA
process for given lag-polynomials.

.. autosummary::
   :toctree: generated/

   arima_process.ArmaProcess
   arima_process.ar2arma
   arima_process.arma2ar
   arima_process.arma2ma
   arima_process.arma_acf
   arima_process.arma_acovf
   arima_process.arma_generate_sample
   arima_process.arma_impulse_response
   arima_process.arma_pacf
   arima_process.arma_periodogram
   arima_process.deconvolve
   arima_process.index2lpol
   arima_process.lpol2index
   arima_process.lpol_fiar
   arima_process.lpol_fima
   arima_process.lpol_sdiff

.. currentmodule:: statsmodels.sandbox.tsa.fftarma

.. autosummary::
   :toctree: generated/

   ArmaFft

.. currentmodule:: statsmodels.tsa

Autoregressive Distributed Lag (ARDL) Models
""""""""""""""""""""""""""""""""""""""""""""
Autoregressive Distributed Lag models span the space between
autoregressive models (:class:`~statsmodels.tsa.ar_model.AutoReg`)
and vector autoregressive models (:class:`~statsmodels.tsa.vector_ar.VAR`).

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   ardl.ARDL
   ardl.ARDLResults
   ardl.ardl_select_order
   ardl.ARDLOrderSelectionResults

The `ardl.ARDL` model estimates parameters using conditional MLE (OLS)
and allows for both simple deterministic terms (trends and seasonal
dummies) as well as complex deterministics using a
:class:`~statsmodels.tsa.deterministic.DeterministicProcess`.

AR-X and related models can also be fitted with
:class:`~statsmodels.tsa.statespace.sarimax.SARIMAX` class (using full MLE via
the Kalman Filter).

Error Correction Models (ECM)
"""""""""""""""""""""""""""""
Error correction models are reparameterizations of ARDL models that
regress the difference of the endogenous variable on the lagged levels
of the endogenous variables and optional lagged differences of the
exogenous variables.

.. currentmodule:: statsmodels.tsa

.. autosummary::
   :toctree: generated/

   ardl.UECM
   ardl.UECMResults
   ardl.BoundsTestResult


Statespace Models
"""""""""""""""""
See the :ref:`statespace documentation <statespace>`.


Vector ARs and Vector Error Correction Models
"""""""""""""""""""""""""""""""""""""""""""""
See the :ref:`vector_ar documentation. <var>`

Regime switching models
"""""""""""""""""""""""

.. currentmodule:: statsmodels.tsa.regime_switching.markov_regression
.. autosummary::
   :toctree: generated/

   MarkovRegression

.. currentmodule:: statsmodels.tsa.regime_switching.markov_autoregression
.. autosummary::
   :toctree: generated/

   MarkovAutoregression


Time Series Filters
"""""""""""""""""""

.. currentmodule:: statsmodels.tsa.filters.bk_filter
.. autosummary::
   :toctree: generated/

   bkfilter

.. currentmodule:: statsmodels.tsa.filters.hp_filter
.. autosummary::
   :toctree: generated/

   hpfilter

.. currentmodule:: statsmodels.tsa.filters.cf_filter
.. autosummary::
   :toctree: generated/

   cffilter

.. currentmodule:: statsmodels.tsa.filters.filtertools
.. autosummary::
   :toctree: generated/

   convolution_filter
   recursive_filter
   miso_lfilter
   fftconvolve3
   fftconvolveinv


.. currentmodule:: statsmodels.tsa.seasonal
.. autosummary::
   :toctree: generated/

   seasonal_decompose
   STL
   DecomposeResult

TSA Tools
"""""""""

.. currentmodule:: statsmodels.tsa.tsatools

.. autosummary::
   :toctree: generated/

   add_lag
   add_trend
   detrend
   lagmat
   lagmat2ds

VARMA Process
"""""""""""""

.. currentmodule:: statsmodels.tsa.varma_process
.. autosummary::
   :toctree: generated/

   VarmaPoly

Interpolation
"""""""""""""

.. currentmodule:: statsmodels.tsa.interp.denton
.. autosummary::
   :toctree: generated/

   dentonm


Deterministic Processes
"""""""""""""""""""""""

Deterministic processes simplify creating deterministic sequences with time
trend or seasonal patterns. They also provide methods to simplify generating
deterministic terms for out-of-sample forecasting. A
:class:`~statsmodels.tsa.deterministic.DeterministicProcess` can be directly
used with :class:`~statsmodels.tsa.ar_model.AutoReg` to construct complex
deterministic dynamics and to forecast without constructing exogenous trends.

.. currentmodule:: statsmodels.tsa.deterministic
.. autosummary::
   :toctree: generated/

   DeterministicProcess
   TimeTrend
   Seasonality
   Fourier
   CalendarTimeTrend
   CalendarSeasonality
   CalendarFourier
   DeterministicTerm
   CalendarDeterministicTerm
   FourierDeterministicTerm
   TimeTrendDeterministicTerm

Users who wish to write custom deterministic terms must subclass
:class:`~statsmodels.tsa.deterministic.DeterministicTerm`.

.. currentmodule:: statsmodels.tsa.deterministic
.. autosummary::
   :toctree: generated/

   DeterministicTerm

Forecasting Models
""""""""""""""""""
.. module:: statsmodels.tsa.forecasting
   :synopsis: Models designed for forecasting

.. currentmodule:: statsmodels.tsa.forecasting

The Theta Model
~~~~~~~~~~~~~~~
The Theta model is a simple forecasting method that combines a linear time
trend with a Simple Exponential Smoother (Assimakopoulos & Nikolopoulos).
An estimator for the parameters of the Theta model and methods to forecast
are available in:

.. module:: statsmodels.tsa.forecasting.theta
   :synopsis: Models designed for forecasting

.. currentmodule:: statsmodels.tsa.forecasting.theta

.. autosummary::
   :toctree: generated/

   ThetaModel
   ThetaModelResults

Forecasting after STL Decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`statsmodels.tsa.seasonal.STL` is commonly used to remove seasonal
components from a time series. The deseasonalized time series can then
be modeled using a any non-seasonal model, and forecasts are constructed
by adding the forecast from the non-seasonal model to the estimates of
the seasonal component from the final full-cycle which are forecast using
a random-walk model.

.. module:: statsmodels.tsa.forecasting.stl
   :synopsis: Models designed for forecasting

.. currentmodule:: statsmodels.tsa.forecasting.stl

.. autosummary::
   :toctree: generated/

   STLForecast
   STLForecastResults

Prediction Results
""""""""""""""""""
Most forecasting methods support a ``get_prediction`` method that return
a ``PredictionResults`` object that contains both the prediction, its
variance and can construct a prediction interval.

Results Class
~~~~~~~~~~~~~

.. module:: statsmodels.tsa.base.prediction
   :synopsis: Shared objects for predictive methods

.. currentmodule:: statsmodels.tsa.base.prediction

.. autosummary::
   :toctree: generated/

   PredictionResults

.. module:: statsmodels.sandbox.distributions
   :synopsis: Probability distributions

.. currentmodule:: statsmodels.sandbox.distributions

.. _distributions:


Distributions
=============

This section collects various additional functions and methods for statistical
distributions.

Empirical Distributions
-----------------------

.. module:: statsmodels.distributions.empirical_distribution
   :synopsis: Tools for working with empirical distributions

.. currentmodule:: statsmodels.distributions.empirical_distribution

.. autosummary::
   :toctree: generated/

   ECDF
   StepFunction
   monotone_fn_inverter

Count Distributions
-------------------

The `discrete` module contains classes for count distributions that are based
on discretizing a continuous distribution, and specific count distributions
that are not available in scipy.distributions like generalized poisson and
zero-inflated count models.

The latter are mainly in support of the corresponding models in
`statsmodels.discrete`. Some methods are not specifically implemented and will
use potentially slow inherited generic methods.

.. module:: statsmodels.distributions.discrete
   :synopsis: Support for count distributions

.. currentmodule:: statsmodels.distributions.discrete

.. autosummary::
   :toctree: generated/

   DiscretizedCount
   DiscretizedModel
   genpoisson_p
   zigenpoisson
   zinegbin
   zipoisson

Copula
------

The `copula` sub-module provides classes to model the dependence between
parameters. Copulae are used to construct a multivariate joint distribution and
provide a set of functions like sampling, PDF, CDF.

.. module:: statsmodels.distributions.copula.api
   :synopsis: Copula for modeling parameter dependence

.. currentmodule:: statsmodels.distributions.copula.api

.. autosummary::
   :toctree: generated/

   CopulaDistribution
   ArchimedeanCopula
   FrankCopula
   ClaytonCopula
   GumbelCopula
   GaussianCopula
   StudentTCopula
   ExtremeValueCopula
   IndependenceCopula

Distribution Extras
-------------------


.. module:: statsmodels.sandbox.distributions.extras
   :synopsis: Probability distributions and random number generators

.. currentmodule:: statsmodels.sandbox.distributions.extras

*Skew Distributions*

.. autosummary::
   :toctree: generated/

   SkewNorm_gen
   SkewNorm2_gen
   ACSkewT_gen
   skewnorm2

*Distributions based on Gram-Charlier expansion*

.. autosummary::
   :toctree: generated/

   pdf_moments_st
   pdf_mvsk
   pdf_moments
   NormExpan_gen

*cdf of multivariate normal* wrapper for scipy.stats


.. autosummary::
   :toctree: generated/

   mvstdnormcdf
   mvnormcdf

Univariate Distributions by non-linear Transformations
------------------------------------------------------

Univariate distributions can be generated from a non-linear transformation of an
existing univariate distribution. `Transf_gen` is a class that can generate a new
distribution from a monotonic transformation, `TransfTwo_gen` can use hump-shaped
or u-shaped transformation, such as abs or square. The remaining objects are
special cases.

.. module:: statsmodels.sandbox.distributions.transformed
   :synopsis: Experimental probability distributions and random number generators

.. currentmodule:: statsmodels.sandbox.distributions.transformed

.. autosummary::
   :toctree: generated/

   TransfTwo_gen
   Transf_gen

   ExpTransf_gen
   LogTransf_gen
   SquareFunc

   absnormalg
   invdnormalg

   loggammaexpg
   lognormalg
   negsquarenormalg

   squarenormalg
   squaretg


Helper Functions
----------------

.. module:: statsmodels.tools.rng_qrng
   :synopsis: Tools for working with random variable generation

.. currentmodule:: statsmodels.tools.rng_qrng

.. autosummary::
   :toctree: generated/

   check_random_state
.. currentmodule:: statsmodels.emplike


.. _emplike:


Empirical Likelihood :mod:`emplike`
====================================


Introduction
------------

Empirical likelihood is a method of nonparametric inference and estimation that lifts the
obligation of having to specify a family of underlying distributions.  Moreover, empirical
likelihood methods do not require re-sampling but still
uniquely determine confidence regions whose shape mirrors the shape of the data.
In essence, empirical likelihood attempts to combine the benefits of parametric
and nonparametric methods while limiting their shortcomings.  The main difficulties  of
empirical likelihood is the computationally intensive methods required to conduct inference.
:mod:`statsmodels.emplike` attempts to provide a user-friendly interface that allows the
end user to effectively conduct empirical likelihood analysis without having to concern
themselves with the computational burdens.

Currently, :mod:`emplike` provides methods to conduct hypothesis tests and form confidence
intervals for descriptive statistics.  Empirical likelihood estimation and inference
in a regression, accelerated failure time and instrumental variable model are
currently under development.

References
^^^^^^^^^^

The main reference for empirical likelihood is::

    Owen, A.B. "Empirical Likelihood." Chapman and Hall, 2001.



Examples
--------

.. ipython:: python

  import numpy as np
  import statsmodels.api as sm

  # Generate Data
  x = np.random.standard_normal(50)

  # initiate EL
  el = sm.emplike.DescStat(x)

  # confidence interval for the mean
  el.ci_mean()

  # test variance is 1
  el.test_var(1)


Module Reference
----------------

.. module:: statsmodels.emplike
   :synopsis: Empirical likelihood tools

.. autosummary::
   :toctree: generated/

   descriptive.DescStat
   descriptive.DescStatUV
   descriptive.DescStatMV
.. currentmodule:: statsmodels.robust


.. _rlm:

Robust Linear Models
====================

Robust linear models with support for the M-estimators listed under `Norms`_.

See `Module Reference`_ for commands and arguments.

Examples
--------

.. ipython:: python

    # Load modules and data
    import statsmodels.api as sm
    data = sm.datasets.stackloss.load()
    data.exog = sm.add_constant(data.exog)

    # Fit model and print summary
    rlm_model = sm.RLM(data.endog, data.exog, M=sm.robust.norms.HuberT())
    rlm_results = rlm_model.fit()
    print(rlm_results.params)

Detailed examples can be found here:

* `Robust Models 1 <examples/notebooks/generated/robust_models_0.html>`__
* `Robust Models 2 <examples/notebooks/generated/robust_models_1.html>`__

Technical Documentation
-----------------------

.. toctree::
   :maxdepth: 1

   rlm_techn1

References
^^^^^^^^^^

* PJ Huber. ‘Robust Statistics’ John Wiley and Sons, Inc., New York. 1981.
* PJ Huber. 1973, ‘The 1972 Wald Memorial Lectures: Robust Regression: Asymptotics, Conjectures, and Monte Carlo.’ The Annals of Statistics, 1.5, 799-821.
* R Venables, B Ripley. ‘Modern Applied Statistics in S’ Springer, New York,
* C Croux, PJ Rousseeuw, 'Time-efficient algorithms for two highly robust estimators of scale' Computational statistics. Physica, Heidelberg, 1992.

Module Reference
----------------

.. module:: statsmodels.robust

Model Classes
^^^^^^^^^^^^^

.. module:: statsmodels.robust.robust_linear_model
.. currentmodule:: statsmodels.robust.robust_linear_model

.. autosummary::
   :toctree: generated/

   RLM

Model Results
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   RLMResults

.. _norms:

Norms
^^^^^

.. module:: statsmodels.robust.norms
.. currentmodule:: statsmodels.robust.norms

.. autosummary::
   :toctree: generated/

   AndrewWave
   Hampel
   HuberT
   LeastSquares
   RamsayE
   RobustNorm
   TrimmedMean
   TukeyBiweight
   estimate_location


Scale
^^^^^

.. module:: statsmodels.robust.scale
.. currentmodule:: statsmodels.robust.scale

.. autosummary::
   :toctree: generated/

    Huber
    HuberScale
    mad
    hubers_scale
    iqr
    qn_scale
.. _sandbox:


Sandbox
=======

This sandbox contains code that is for various reasons not ready to be
included in statsmodels proper. It contains modules from the old stats.models
code that have not been tested, verified and updated to the new statsmodels
structure: cox survival model, mixed effects model with repeated measures,
generalized additive model and the formula framework. The sandbox also
contains code that is currently being worked on until it fits the pattern
of statsmodels or is sufficiently tested.

All sandbox modules have to be explicitly imported to indicate that they are
not yet part of the core of statsmodels. The quality and testing of the
sandbox code varies widely.


Examples
--------

There are some examples in the `sandbox.examples` folder. Additional
examples are directly included in the modules and in subfolders of
the sandbox.


Module Reference
----------------


Time Series analysis :mod:`tsa`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this part we develop models and functions that will be useful for time
series analysis. Most of the models and function have been moved to
:mod:`statsmodels.tsa`.

Moving Window Statistics
""""""""""""""""""""""""

Most moving window statistics, like rolling mean, moments (up to 4th order), min,
max, mean, and variance, are covered by the functions for `Moving (rolling)
statistics/moments <https://pandas.pydata.org/pandas-docs/stable/user_guide/computation.html#window-functions>`_ in Pandas.

.. module:: statsmodels.sandbox.tsa
   :synopsis: Experimental time-series analysis models

.. currentmodule:: statsmodels.sandbox.tsa

.. autosummary::
   :toctree: generated/

   movstat.movorder
   movstat.movmean
   movstat.movvar
   movstat.movmoment


Regression and ANOVA
^^^^^^^^^^^^^^^^^^^^

.. module:: statsmodels.sandbox.regression.anova_nistcertified
   :synopsis: Experimental ANOVA estimator

.. currentmodule:: statsmodels.sandbox.regression.anova_nistcertified

The following two ANOVA functions are fully tested against the NIST test data
for balanced one-way ANOVA. ``anova_oneway`` follows the same pattern as the
oneway anova function in scipy.stats but with higher precision for badly
scaled problems. ``anova_ols`` produces the same results as the one way anova
however using the OLS model class. It also verifies against the NIST tests,
with some problems in the worst scaled cases. It shows how to do simple ANOVA
using statsmodels in three lines and is also best taken as a recipe.


.. autosummary::
   :toctree: generated/

   anova_oneway
   anova_ols


The following are helper functions for working with dummy variables and
generating ANOVA results with OLS. They are best considered as recipes since
they were written with a specific use in mind. These function will eventually
be rewritten or reorganized.

.. module:: statsmodels.sandbox.regression
   :synopsis: Experimental regression tools

.. currentmodule:: statsmodels.sandbox.regression

.. autosummary::
   :toctree: generated/

   try_ols_anova.data2dummy
   try_ols_anova.data2groupcont
   try_ols_anova.data2proddummy
   try_ols_anova.dropname
   try_ols_anova.form2design

The following are helper functions for group statistics where groups are
defined by a label array. The qualifying comments for the previous group
apply also to this group of functions.


.. autosummary::
   :toctree: generated/

   try_catdata.cat2dummy
   try_catdata.convertlabels
   try_catdata.groupsstats_1d
   try_catdata.groupsstats_dummy
   try_catdata.groupstatsbin
   try_catdata.labelmeanfilter
   try_catdata.labelmeanfilter_nd
   try_catdata.labelmeanfilter_str

Additional to these functions, sandbox regression still contains several
examples, that are illustrative of the use of the regression models of
statsmodels.



Systems of Regression Equations and Simultaneous Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following are for fitting systems of equations models.  Though the returned
parameters have been verified as accurate, this code is still very
experimental, and the usage of the models will very likely change significantly
before they are added to the main codebase.

.. module:: statsmodels.sandbox.sysreg
   :synopsis: Experimental system regression models

.. currentmodule:: statsmodels.sandbox.sysreg

.. autosummary::
   :toctree: generated/

   SUR
   Sem2SLS

Miscellaneous
^^^^^^^^^^^^^
.. module:: statsmodels.sandbox.tools.tools_tsa
   :synopsis: Experimental tools for working with time-series

.. currentmodule:: statsmodels.sandbox.tools.tools_tsa


Descriptive Statistics Printing
"""""""""""""""""""""""""""""""

.. module:: statsmodels.sandbox
   :synopsis: Experimental tools that have not been fully vetted

.. currentmodule:: statsmodels.sandbox

.. autosummary::
   :toctree: generated/

   descstats.sign_test
   descstats.descstats




Original stats.models
^^^^^^^^^^^^^^^^^^^^^

None of these are fully working. The formula framework is used by cox and
mixed.

**Mixed Effects Model with Repeated Measures using an EM Algorithm**

:mod:`statsmodels.sandbox.mixed`


**Cox Proportional Hazards Model**

:mod:`statsmodels.sandbox.cox`

**Generalized Additive Models**

:mod:`statsmodels.sandbox.gam`

**Formula**

:mod:`statsmodels.sandbox.formula`


.. currentmodule:: statsmodels.nonparametric

.. _nonparametric:


Nonparametric Methods :mod:`nonparametric`
==========================================

This section collects various methods in nonparametric statistics. This
includes kernel density estimation for univariate and multivariate data,
kernel regression and locally weighted scatterplot smoothing (lowess).

sandbox.nonparametric contains additional functions that are work in progress
or do not have unit tests yet. We are planning to include here nonparametric
density estimators, especially based on kernel or orthogonal polynomials,
smoothers, and tools for nonparametric models and methods in other parts of
statsmodels.


Kernel density estimation
-------------------------

The kernel density estimation (KDE) functionality is split between univariate
and multivariate estimation, which are implemented in quite different ways.

Univariate estimation (as provided by `KDEUnivariate`) uses FFT transforms,
which makes it quite fast.  Therefore it should be preferred for *continuous,
univariate* data if speed is important.  It supports using different kernels;
bandwidth estimation is done only by a rule of thumb (Scott or Silverman).

Multivariate estimation (as provided by `KDEMultivariate`) uses product
kernels.   It supports least squares and maximum likelihood cross-validation
for bandwidth estimation, as well as estimating mixed continuous, ordered and
unordered data.  The default kernels (Gaussian, Wang-Ryzin and
Aitchison-Aitken) cannot be altered at the moment however.  Direct estimation
of the conditional density (:math:`P(X | Y) = P(X, Y) / P(Y)`) is supported
by `KDEMultivariateConditional`.

`KDEMultivariate` can do univariate estimation as well, but is up to two orders
of magnitude slower than `KDEUnivariate`.


Kernel regression
-----------------

Kernel regression (as provided by `KernelReg`) is based on the same product
kernel approach as `KDEMultivariate`, and therefore has the same set of
features (mixed data, cross-validated bandwidth estimation, kernels) as
described above for `KDEMultivariate`.  Censored regression is provided by
`KernelCensoredReg`.

Note that code for semi-parametric partial linear models and single index
models, based on `KernelReg`, can be found in the sandbox.


References
----------

* B.W. Silverman, "Density Estimation for Statistics and Data Analysis"
* J.S. Racine, "Nonparametric Econometrics: A Primer," Foundation and
  Trends in Econometrics, Vol. 3, No. 1, pp. 1-88, 2008.
* Q. Li and J.S. Racine, "Nonparametric econometrics: theory and practice",
  Princeton University Press, 2006.
* Hastie, Tibshirani and Friedman, "The Elements of Statistical Learning:
  Data Mining, Inference, and Prediction", Springer, 2009.
* Racine, J., Li, Q. "Nonparametric Estimation of Distributions
  with Categorical and Continuous Data." Working Paper. (2000)
* Racine, J. Li, Q. "Kernel Estimation of Multivariate Conditional
  Distributions Annals of Economics and Finance 5, 211-235 (2004)
* Liu, R., Yang, L. "Kernel estimation of multivariate
  cumulative distribution function." Journal of Nonparametric Statistics
  (2008)
* Li, R., Ju, G. "Nonparametric Estimation of Multivariate CDF
  with Categorical and Continuous Data." Working Paper
* Li, Q., Racine, J. "Cross-validated local linear nonparametric
  regression" Statistica Sinica 14(2004), pp. 485-512
* Racine, J.: "Consistent Significance Testing for Nonparametric
  Regression" Journal of Business & Economics Statistics
* Racine, J., Hart, J., Li, Q., "Testing the Significance of
  Categorical Predictor Variables in Nonparametric Regression
  Models", 2006, Econometric Reviews 25, 523-544


Module Reference
----------------

.. module:: statsmodels.nonparametric
   :synopsis: Nonparametric estimation of densities and curves

The public functions and classes are

.. currentmodule:: statsmodels.nonparametric.smoothers_lowess
.. autosummary::
   :toctree: generated/

   lowess

.. currentmodule:: statsmodels.nonparametric.kde
.. autosummary::
   :toctree: generated/

   KDEUnivariate

.. currentmodule:: statsmodels.nonparametric.kernel_density
.. autosummary::
   :toctree: generated/

   KDEMultivariate
   KDEMultivariateConditional
   EstimatorSettings

.. currentmodule:: statsmodels.nonparametric.kernel_regression
.. autosummary::
   :toctree: generated/

   KernelReg
   KernelCensoredReg

helper functions for kernel bandwidths

.. currentmodule:: statsmodels.nonparametric.bandwidths
.. autosummary::
   :toctree: generated/

   bw_scott
   bw_silverman
   select_bandwidth

There are some examples for nonlinear functions in
:mod:`statsmodels.nonparametric.dgp_examples`


Asymmetric Kernels
------------------

Asymmetric kernels like beta for the unit interval and gamma for positive
valued random variables avoid problems at the boundary of the support of the
distribution.

Statsmodels has preliminary support for estimating density and cumulative
distribution function using kernels for the unit interval, ``beta`` or the
positive real line, all other kernels.

Several of the kernels for the positive real line assume that the density at
the zero boundary is zero. The gamma kernel also allows the case of positive
or unbound density at the zero boundary.

There are currently no defaults and no support for choosing the bandwidth. the
user has to provide the bandwidth.

The functions to compute kernel density and kernel cdf are

.. currentmodule:: statsmodels.nonparametric.kernels_asymmetric
.. autosummary::
   :toctree: generated/

   pdf_kernel_asym
   cdf_kernel_asym

The available kernel functions for pdf and cdf are

.. autosummary::
   :toctree: generated/

   kernel_pdf_beta
   kernel_pdf_beta2
   kernel_pdf_bs
   kernel_pdf_gamma
   kernel_pdf_gamma2
   kernel_pdf_invgamma
   kernel_pdf_invgauss
   kernel_pdf_lognorm
   kernel_pdf_recipinvgauss
   kernel_pdf_weibull
   kernel_cdf_beta
   kernel_cdf_beta2
   kernel_cdf_bs
   kernel_cdf_gamma
   kernel_cdf_gamma2
   kernel_cdf_invgamma
   kernel_cdf_invgauss
   kernel_cdf_lognorm
   kernel_cdf_recipinvgauss
   kernel_cdf_weibull


The sandbox.nonparametric contains additional insufficiently tested classes
for testing functional form and for semi-linear and single index models.
.. currentmodule:: statsmodels.genmod.generalized_linear_model

.. _glm:

Generalized Linear Models
=========================

Generalized linear models currently supports estimation using the one-parameter
exponential families.

See `Module Reference`_ for commands and arguments.

Examples
--------

.. ipython:: python
   :okwarning:

   # Load modules and data
   import statsmodels.api as sm
   data = sm.datasets.scotland.load()
   data.exog = sm.add_constant(data.exog)

   # Instantiate a gamma family model with the default link function.
   gamma_model = sm.GLM(data.endog, data.exog, family=sm.families.Gamma())
   gamma_results = gamma_model.fit()
   print(gamma_results.summary())

Detailed examples can be found here:

* `GLM <examples/notebooks/generated/glm.html>`__
* `Formula <examples/notebooks/generated/glm_formula.html>`__

Technical Documentation
-----------------------

..   ..glm_techn1
..   ..glm_techn2

The statistical model for each observation :math:`i` is assumed to be

 :math:`Y_i \sim F_{EDM}(\cdot|\theta,\phi,w_i)` and
 :math:`\mu_i = E[Y_i|x_i] = g^{-1}(x_i^\prime\beta)`.

where :math:`g` is the link function and :math:`F_{EDM}(\cdot|\theta,\phi,w)`
is a distribution of the family of exponential dispersion models (EDM) with
natural parameter :math:`\theta`, scale parameter :math:`\phi` and weight
:math:`w`.
Its density is given by

 :math:`f_{EDM}(y|\theta,\phi,w) = c(y,\phi,w)
 \exp\left(\frac{y\theta-b(\theta)}{\phi}w\right)\,.`

It follows that :math:`\mu = b'(\theta)` and
:math:`Var[Y|x]=\frac{\phi}{w}b''(\theta)`. The inverse of the first equation
gives the natural parameter as a function of the expected value
:math:`\theta(\mu)` such that

 :math:`Var[Y_i|x_i] = \frac{\phi}{w_i} v(\mu_i)`

with :math:`v(\mu) = b''(\theta(\mu))`. Therefore it is said that a GLM is
determined by link function :math:`g` and variance function :math:`v(\mu)`
alone (and :math:`x` of course).

Note that while :math:`\phi` is the same for every observation :math:`y_i`
and therefore does not influence the estimation of :math:`\beta`,
the weights :math:`w_i` might be different for every :math:`y_i` such that the
estimation of :math:`\beta` depends on them.

================================================= ============================== ============================== ======================================== =========================================== ============================================================================ =====================
Distribution                                      Domain                         :math:`\mu=E[Y|x]`             :math:`v(\mu)`                           :math:`\theta(\mu)`                         :math:`b(\theta)`                                                            :math:`\phi`
================================================= ============================== ============================== ======================================== =========================================== ============================================================================ =====================
Binomial :math:`B(n,p)`                           :math:`0,1,\ldots,n`           :math:`np`                     :math:`\mu-\frac{\mu^2}{n}`              :math:`\log\frac{p}{1-p}`                   :math:`n\log(1+e^\theta)`                                                    1
Poisson :math:`P(\mu)`                            :math:`0,1,\ldots,\infty`      :math:`\mu`                    :math:`\mu`                              :math:`\log(\mu)`                           :math:`e^\theta`                                                             1
Neg. Binom. :math:`NB(\mu,\alpha)`                :math:`0,1,\ldots,\infty`      :math:`\mu`                    :math:`\mu+\alpha\mu^2`                  :math:`\log(\frac{\alpha\mu}{1+\alpha\mu})` :math:`-\frac{1}{\alpha}\log(1-\alpha e^\theta)`                             1
Gaussian/Normal :math:`N(\mu,\sigma^2)`           :math:`(-\infty,\infty)`       :math:`\mu`                    :math:`1`                                :math:`\mu`                                 :math:`\frac{1}{2}\theta^2`                                                  :math:`\sigma^2`
Gamma :math:`N(\mu,\nu)`                          :math:`(0,\infty)`             :math:`\mu`                    :math:`\mu^2`                            :math:`-\frac{1}{\mu}`                      :math:`-\log(-\theta)`                                                       :math:`\frac{1}{\nu}`
Inv. Gauss. :math:`IG(\mu,\sigma^2)`              :math:`(0,\infty)`             :math:`\mu`                    :math:`\mu^3`                            :math:`-\frac{1}{2\mu^2}`                   :math:`-\sqrt{-2\theta}`                                                     :math:`\sigma^2`
Tweedie :math:`p\geq 1`                           depends on :math:`p`           :math:`\mu`                    :math:`\mu^p`                            :math:`\frac{\mu^{1-p}}{1-p}`               :math:`\frac{\alpha-1}{\alpha}\left(\frac{\theta}{\alpha-1}\right)^{\alpha}` :math:`\phi`
================================================= ============================== ============================== ======================================== =========================================== ============================================================================ =====================

The Tweedie distribution has special cases for :math:`p=0,1,2` not listed in the
table and uses :math:`\alpha=\frac{p-2}{p-1}`.

Correspondence of mathematical variables to code:

* :math:`Y` and :math:`y` are coded as ``endog``, the variable one wants to
  model
* :math:`x` is coded as ``exog``, the covariates alias explanatory variables
* :math:`\beta` is coded as ``params``, the parameters one wants to estimate
* :math:`\mu` is coded as ``mu``, the expectation (conditional on :math:`x`)
  of :math:`Y`
* :math:`g` is coded as ``link`` argument to the ``class Family``
* :math:`\phi` is coded as ``scale``, the dispersion parameter of the EDM
* :math:`w` is not yet supported (i.e. :math:`w=1`), in the future it might be
  ``var_weights``
* :math:`p` is coded as ``var_power`` for the power of the variance function
  :math:`v(\mu)` of the Tweedie distribution, see table
* :math:`\alpha` is either

  * Negative Binomial: the ancillary parameter ``alpha``, see table
  * Tweedie: an abbreviation for :math:`\frac{p-2}{p-1}` of the power :math:`p`
    of the variance function, see table


References
^^^^^^^^^^

* Gill, Jeff. 2000. Generalized Linear Models: A Unified Approach. SAGE QASS Series.
* Green, PJ. 1984. “Iteratively reweighted least squares for maximum likelihood estimation, and some robust and resistant alternatives.” Journal of the Royal Statistical Society, Series B, 46, 149-192.
* Hardin, J.W. and Hilbe, J.M. 2007. “Generalized Linear Models and Extensions.” 2nd ed. Stata Press, College Station, TX.
* McCullagh, P. and Nelder, J.A. 1989. “Generalized Linear Models.” 2nd ed. Chapman & Hall, Boca Rotan.

Module Reference
----------------

.. module:: statsmodels.genmod.generalized_linear_model
   :synopsis: Generalized Linear Models (GLM)

Model Class
^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   GLM

Results Class
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   GLMResults
   PredictionResultsMean

.. _families:

Families
^^^^^^^^

The distribution families currently implemented are

.. module:: statsmodels.genmod.families.family
.. currentmodule:: statsmodels.genmod.families.family

.. autosummary::
   :toctree: generated/

   Family
   Binomial
   Gamma
   Gaussian
   InverseGaussian
   NegativeBinomial
   Poisson
   Tweedie


.. _links:

Link Functions
^^^^^^^^^^^^^^

The link functions currently implemented are the following. Not all link
functions are available for each distribution family. The list of
available link functions can be obtained by

::

    >>> sm.families.family.<familyname>.links

.. module:: statsmodels.genmod.families.links
.. currentmodule:: statsmodels.genmod.families.links

.. autosummary::
   :toctree: generated/

   Link
   CDFLink
   CLogLog
   LogLog
   Log
   Logit
   NegativeBinomial
   Power
   cauchy
   cloglog
   loglog
   identity
   inverse_power
   inverse_squared
   log
   logit
   nbinom
   probit

.. _varfuncs:

Variance Functions
^^^^^^^^^^^^^^^^^^

Each of the families has an associated variance function. You can access
the variance functions here:

::

    >>> sm.families.<familyname>.variance

.. module:: statsmodels.genmod.families.varfuncs
.. currentmodule:: statsmodels.genmod.families.varfuncs

.. autosummary::
   :toctree: generated/

   VarianceFunction
   constant
   Power
   mu
   mu_squared
   mu_cubed
   Binomial
   binary
   NegativeBinomial
   nbinom
.. currentmodule:: statsmodels.graphics.gofplots

.. _graphics:

Graphics
========

.. automodule:: statsmodels.graphics

Goodness of Fit Plots
---------------------

.. currentmodule:: statsmodels.graphics.gofplots
.. autosummary::
   :toctree: generated/

   qqplot
   qqline
   qqplot_2samples
   ProbPlot

Boxplots
--------

.. currentmodule:: statsmodels.graphics.boxplots

.. autosummary::
   :toctree: generated/

   violinplot
   beanplot

Correlation Plots
------------------

.. currentmodule:: statsmodels.graphics.correlation

.. autosummary::
   :toctree: generated/

   plot_corr
   plot_corr_grid

.. currentmodule:: statsmodels.graphics.plot_grids

.. autosummary::
   :toctree: generated/

   scatter_ellipse

Dot Plots
---------

.. currentmodule:: statsmodels.graphics.dotplots

.. autosummary::
   :toctree: generated/

   dot_plot

Functional Plots
----------------

.. currentmodule:: statsmodels.graphics.functional

.. autosummary::
   :toctree: generated/

   hdrboxplot
   fboxplot
   rainbowplot
   banddepth

Regression Plots
----------------

.. currentmodule:: statsmodels.graphics.regressionplots

.. autosummary::
   :toctree: generated/

   plot_fit
   plot_regress_exog
   plot_partregress
   plot_partregress_grid
   plot_ccpr
   plot_ccpr_grid
   plot_ceres_residuals
   abline_plot
   influence_plot
   plot_leverage_resid2

Time Series Plots
-----------------

.. currentmodule:: statsmodels.graphics.tsaplots

.. autosummary::
   :toctree: generated/

   plot_acf
   plot_pacf
   month_plot
   quarter_plot

Other Plots
-----------

.. currentmodule:: statsmodels.graphics.factorplots

.. autosummary::
   :toctree: generated/

   interaction_plot

.. currentmodule:: statsmodels.graphics.mosaicplot

.. autosummary::
   :toctree: generated/

   mosaic

.. currentmodule:: statsmodels.graphics.agreement

.. autosummary::
   :toctree: generated/

   mean_diff_plot
.. module:: statsmodels.multivariate
   :synopsis: Models for multivariate data

.. currentmodule:: statsmodels.multivariate

.. _multivariate:


Multivariate Statistics :mod:`multivariate`
===========================================

This section includes methods and algorithms from multivariate statistics.


Principal Component Analysis
----------------------------

.. module:: statsmodels.multivariate.pca
   :synopsis: Principal Component Analaysis

.. currentmodule:: statsmodels.multivariate.pca

.. autosummary::
   :toctree: generated/

   PCA
   pca


Factor Analysis
---------------

.. currentmodule:: statsmodels.multivariate.factor

.. autosummary::
   :toctree: generated/

   Factor
   FactorResults


Factor Rotation
---------------

.. currentmodule:: statsmodels.multivariate.factor_rotation

.. autosummary::
   :toctree: generated/

   rotate_factors
   target_rotation
   procrustes
   promax


Canonical Correlation
---------------------

.. currentmodule:: statsmodels.multivariate.cancorr

.. autosummary::
   :toctree: generated/

   CanCorr


MANOVA
------

.. currentmodule:: statsmodels.multivariate.manova

.. autosummary::
   :toctree: generated/

   MANOVA


MultivariateOLS
---------------

`_MultivariateOLS` is a model class with limited features. Currently it
supports multivariate hypothesis tests and is used as backend for MANOVA.

.. currentmodule:: statsmodels.multivariate.multivariate_ols

.. autosummary::
   :toctree: generated/

   _MultivariateOLS
   _MultivariateOLSResults
   MultivariateTestResults
.. currentmodule:: statsmodels.tools


.. _tools:

Tools
=====

Our tool collection contains some convenience functions for users and
functions that were written mainly for internal use.

Additional to this tools directory, several other subpackages have their own
tools modules, for example :mod:`statsmodels.tsa.tsatools`


Module Reference
----------------

.. module:: statsmodels.tools
   :synopsis: Tools for variable transformation and common numerical operations

Basic tools :mod:`tools`
^^^^^^^^^^^^^^^^^^^^^^^^

These are basic and miscellaneous tools. The full import path is
`statsmodels.tools.tools`.

.. autosummary::
   :toctree: generated/

   tools.add_constant

The next group are mostly helper functions that are not separately tested or
insufficiently tested.

.. autosummary::
   :toctree: generated/

   tools.clean0
   tools.fullrank
   tools.isestimable
   tools.recipr
   tools.recipr0
   tools.unsqueeze

.. currentmodule:: statsmodels.tools

.. _numdiff:

Numerical Differentiation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   numdiff.approx_fprime
   numdiff.approx_fprime_cs
   numdiff.approx_hess1
   numdiff.approx_hess2
   numdiff.approx_hess3
   numdiff.approx_hess_cs

Measure for fit performance :mod:`eval_measures`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first group of function in this module are standalone versions of
information criteria, aic bic and hqic. The function with `_sigma` suffix
take the error sum of squares as argument, those without, take the value
of the log-likelihood, `llf`, as argument.

The second group of function are measures of fit or prediction performance,
which are mostly one liners to be used as helper functions. All of those
calculate a performance or distance statistic for the difference between two
arrays. For example in the case of Monte Carlo or cross-validation, the first
array would be the estimation results for the different replications or draws,
while the second array would be the true or observed values.

.. currentmodule:: statsmodels.tools

.. autosummary::
   :toctree: generated/

   eval_measures.aic
   eval_measures.aic_sigma
   eval_measures.aicc
   eval_measures.aicc_sigma
   eval_measures.bic
   eval_measures.bic_sigma
   eval_measures.hqic
   eval_measures.hqic_sigma

   eval_measures.bias
   eval_measures.iqr
   eval_measures.maxabs
   eval_measures.meanabs
   eval_measures.medianabs
   eval_measures.medianbias
   eval_measures.mse
   eval_measures.rmse
   eval_measures.rmspe
   eval_measures.stde
   eval_measures.vare
.. module:: statsmodels.rlm
   :synopsis: Outlier robust linear models

.. currentmodule:: statsmodels.rlm


.. _rlm_techn1:

Weight Functions
----------------

Andrew's Wave

.. image:: images/aw.png

Hampel 17A

.. image:: images/hl.png

Huber's t

.. image:: images/ht.png

Least Squares

.. image:: images/ls.png

Ramsay's Ea

.. image:: images/re.png

Trimmed Mean

.. image:: images/tm.png

Tukey's Biweight

.. image:: images/tk.png


User Guide
==========

Background
----------

.. toctree::
   :maxdepth: 1

   endog_exog
   api-structure
   example_formulas
   pitfalls

Regression and Linear Models
----------------------------

.. toctree::
   :maxdepth: 1

   regression
   glm
   gee
   gam
   rlm
   mixed_linear
   discretemod
   mixed_glm
   anova
   other_models

Time Series Analysis
--------------------

.. toctree::
   :maxdepth: 1

   tsa
   statespace
   vector_ar

Other Models
------------

.. toctree::
   :maxdepth: 1

   duration
   nonparametric
   gmm
   miscmodels
   multivariate

Statistics and Tools
--------------------

.. toctree::
   :maxdepth: 1

   stats
   contingency_tables
   imputation
   emplike
   distributions
   graphics
   iolib
   tools
   large_data
   optimization

Data Sets
---------
.. toctree::
   :maxdepth: 1

   datasets/index

Sandbox
-------

.. toctree::
   :maxdepth: 1

   sandbox
.. module:: statsmodels.tsa.statespace
   :synopsis: Statespace models for time-series analysis

.. currentmodule:: statsmodels.tsa.statespace


.. _statespace:


Time Series Analysis by State Space Methods :mod:`statespace`
=============================================================

:mod:`statsmodels.tsa.statespace` contains classes and functions that are
useful for time series analysis using state space methods.

A general state space model is of the form

.. math::

  y_t & = Z_t \alpha_t + d_t + \varepsilon_t \\
  \alpha_{t+1} & = T_t \alpha_t + c_t + R_t \eta_t \\

where :math:`y_t` refers to the observation vector at time :math:`t`,
:math:`\alpha_t` refers to the (unobserved) state vector at time
:math:`t`, and where the irregular components are defined as

.. math::

  \varepsilon_t \sim N(0, H_t) \\
  \eta_t \sim N(0, Q_t) \\

The remaining variables (:math:`Z_t, d_t, H_t, T_t, c_t, R_t, Q_t`) in the
equations are matrices describing the process. Their variable names and
dimensions are as follows

Z : `design`          :math:`(k\_endog \times k\_states \times nobs)`

d : `obs_intercept`   :math:`(k\_endog \times nobs)`

H : `obs_cov`         :math:`(k\_endog \times k\_endog \times nobs)`

T : `transition`      :math:`(k\_states \times k\_states \times nobs)`

c : `state_intercept` :math:`(k\_states \times nobs)`

R : `selection`       :math:`(k\_states \times k\_posdef \times nobs)`

Q : `state_cov`       :math:`(k\_posdef \times k\_posdef \times nobs)`

In the case that one of the matrices is time-invariant (so that, for
example, :math:`Z_t = Z_{t+1} ~ \forall ~ t`), its last dimension may
be of size :math:`1` rather than size `nobs`.

This generic form encapsulates many of the most popular linear time series
models (see below) and is very flexible, allowing estimation with missing
observations, forecasting, impulse response functions, and much more.

**Example: AR(2) model**

An autoregressive model is a good introductory example to putting models in
state space form. Recall that an AR(2) model is often written as:

.. math::

   y_t = \phi_1 y_{t-1} + \phi_2 y_{t-2} + \epsilon_t,
   \quad \epsilon_t \sim N(0, \sigma^2)

This can be put into state space form in the following way:

.. math::

   y_t & = \begin{bmatrix} 1 & 0 \end{bmatrix} \alpha_t \\
   \alpha_{t+1} & = \begin{bmatrix}
      \phi_1 & \phi_2 \\
           1 &      0
   \end{bmatrix} \alpha_t + \begin{bmatrix} 1 \\ 0 \end{bmatrix} \eta_t

Where

.. math::

   Z_t \equiv Z = \begin{bmatrix} 1 & 0 \end{bmatrix}

and

.. math::

   T_t \equiv T & = \begin{bmatrix}
      \phi_1 & \phi_2 \\
           1 &      0
   \end{bmatrix} \\
   R_t \equiv R & = \begin{bmatrix} 1 \\ 0 \end{bmatrix} \\
   \eta_t \equiv \epsilon_{t+1} & \sim N(0, \sigma^2)

There are three unknown parameters in this model:
:math:`\phi_1, \phi_2, \sigma^2`.

Models and Estimation
---------------------

The following are the main estimation classes, which can be accessed through
`statsmodels.tsa.statespace.api` and their result classes.

Seasonal Autoregressive Integrated Moving-Average with eXogenous regressors (SARIMAX)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `SARIMAX` class is an example of a fully fledged model created using the
statespace backend for estimation. `SARIMAX` can be used very similarly to
:ref:`tsa <tsa>` models, but works on a wider range of models by adding the
estimation of additive and multiplicative seasonal effects, as well as
arbitrary trend polynomials.

.. autosummary::
   :toctree: generated/

   sarimax.SARIMAX
   sarimax.SARIMAXResults

For an example of the use of this model, see the
`SARIMAX example notebook <examples/notebooks/generated/statespace_sarimax_stata.html>`__
or the very brief code snippet below:


.. code-block:: python

   # Load the statsmodels api
   import statsmodels.api as sm

   # Load your dataset
   endog = pd.read_csv('your/dataset/here.csv')

   # We could fit an AR(2) model, described above
   mod_ar2 = sm.tsa.SARIMAX(endog, order=(2,0,0))
   # Note that mod_ar2 is an instance of the SARIMAX class

   # Fit the model via maximum likelihood
   res_ar2 = mod_ar2.fit()
   # Note that res_ar2 is an instance of the SARIMAXResults class

   # Show the summary of results
   print(res_ar2.summary())

   # We could also fit a more complicated model with seasonal components.
   # As an example, here is an SARIMA(1,1,1) x (0,1,1,4):
   mod_sarimax = sm.tsa.SARIMAX(endog, order=(1,1,1),
                                seasonal_order=(0,1,1,4))
   res_sarimax = mod_sarimax.fit()

   # Show the summary of results
   print(res_sarimax.summary())

The results object has many of the attributes and methods you would expect from
other statsmodels results objects, including standard errors, z-statistics,
and prediction / forecasting.

Behind the scenes, the `SARIMAX` model creates the design and transition
matrices (and sometimes some of the other matrices) based on the model
specification.

Unobserved Components
^^^^^^^^^^^^^^^^^^^^^

The `UnobservedComponents` class is another example of a statespace model.

.. autosummary::
   :toctree: generated/

   structural.UnobservedComponents
   structural.UnobservedComponentsResults

For examples of the use of this model, see the `example notebook <examples/notebooks/generated/statespace_structural_harvey_jaeger.html>`__ or a notebook on using the unobserved components model to `decompose a time series into a trend and cycle <examples/notebooks/generated/statespace_cycles.html>`__ or the very brief code snippet below:

.. code-block:: python

   # Load the statsmodels api
   import statsmodels.api as sm

   # Load your dataset
   endog = pd.read_csv('your/dataset/here.csv')

   # Fit a local level model
   mod_ll = sm.tsa.UnobservedComponents(endog, 'local level')
   # Note that mod_ll is an instance of the UnobservedComponents class

   # Fit the model via maximum likelihood
   res_ll = mod_ll.fit()
   # Note that res_ll is an instance of the UnobservedComponentsResults class

   # Show the summary of results
   print(res_ll.summary())

   # Show a plot of the estimated level and trend component series
   fig_ll = res_ll.plot_components()

   # We could further add a damped stochastic cycle as follows
   mod_cycle = sm.tsa.UnobservedComponents(endog, 'local level', cycle=True,
                                           damped_cycle=True,
                                           stochastic_cycle=True)
   res_cycle = mod_cycle.fit()

   # Show the summary of results
   print(res_cycle.summary())

   # Show a plot of the estimated level, trend, and cycle component series
   fig_cycle = res_cycle.plot_components()

Vector Autoregressive Moving-Average with eXogenous regressors (VARMAX)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `VARMAX` class is an example of a multivariate statespace model.

.. autosummary::
   :toctree: generated/

   varmax.VARMAX
   varmax.VARMAXResults

For an example of the use of this model, see the `VARMAX example notebook <examples/notebooks/generated/statespace_varmax.html>`__ or the very brief code snippet below:

.. code-block:: python

   # Load the statsmodels api
   import statsmodels.api as sm

   # Load your (multivariate) dataset
   endog = pd.read_csv('your/dataset/here.csv')

   # Fit a local level model
   mod_var1 = sm.tsa.VARMAX(endog, order=(1,0))
   # Note that mod_var1 is an instance of the VARMAX class

   # Fit the model via maximum likelihood
   res_var1 = mod_var1.fit()
   # Note that res_var1 is an instance of the VARMAXResults class

   # Show the summary of results
   print(res_var1.summary())

   # Construct impulse responses
   irfs = res_ll.impulse_responses(steps=10)

Dynamic Factor Models
^^^^^^^^^^^^^^^^^^^^^

Statsmodels has two classes that support dynamic factor models:
`DynamicFactorMQ` and `DynamicFactor`. Each of these models has strengths, but
in general the `DynamicFactorMQ` class is recommended. This is because it fits
parameters using the Expectation-Maximization (EM) algorithm, which is more
robust and can handle including hundreds of observed series. In addition, it
allows customization of which variables load on which factors. However, it does
not yet support including exogenous variables, while `DynamicFactor` does
support that feature.

.. autosummary::
   :toctree: generated/

   dynamic_factor_mq.DynamicFactorMQ
   dynamic_factor_mq.DynamicFactorMQResults

For an example of the `DynamicFactorMQ` class, see the very brief code snippet below:

.. code-block:: python

   # Load the statsmodels api
   import statsmodels.api as sm

   # Load your dataset
   endog = pd.read_csv('your/dataset/here.csv')

   # Create a dynamic factor model
   mod_dfm = sm.tsa.DynamicFactorMQ(endog, k_factors=1, factor_order=2)
   # Note that mod_dfm is an instance of the DynamicFactorMQ class

   # Fit the model via maximum likelihood, using the EM algorithm
   res_dfm = mod_dfm.fit()
   # Note that res_dfm is an instance of the DynamicFactorMQResults class

   # Show the summary of results
   print(res_ll.summary())

   # Show a plot of the r^2 values from regressions of
   # individual estimated factors on endogenous variables.
   fig_dfm = res_ll.plot_coefficients_of_determination()


The `DynamicFactor` class is suitable for models with a smaller number of
observed variables

.. autosummary::
   :toctree: generated/

   dynamic_factor.DynamicFactor
   dynamic_factor.DynamicFactorResults

For an example of the use of the `DynamicFactor` model, see the
`Dynamic Factor example notebook <examples/notebooks/generated/statespace_dfm_coincident.html>`__ 

Linear Exponential Smoothing Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `ExponentialSmoothing` class is an implementation of linear exponential
smoothing models using a state space approach.

**Note**: this model is available at `sm.tsa.statespace.ExponentialSmoothing`;
it is not the same as the model available at `sm.tsa.ExponentialSmoothing`.
See below for details of the differences between these classes.

.. autosummary::
   :toctree: generated/

   exponential_smoothing.ExponentialSmoothing
   exponential_smoothing.ExponentialSmoothingResults

A very brief code snippet follows:

.. code-block:: python

   # Load the statsmodels api
   import statsmodels.api as sm

   # Load your dataset
   endog = pd.read_csv('your/dataset/here.csv')

   # Simple exponential smoothing, denoted (A,N,N)
   mod_ses = sm.tsa.statespace.ExponentialSmoothing(endog)
   res_ses = mod_ses.fit()

   # Holt's linear method, denoted (A,A,N)
   mod_h = sm.tsa.statespace.ExponentialSmoothing(endog, trend=True)
   res_h = mod_h.fit()

   # Damped trend model, denoted (A,Ad,N)
   mod_dt = sm.tsa.statespace.ExponentialSmoothing(endog, trend=True,
                                                   damped_trend=True)
   res_dt = mod_dt.fit()

   # Holt-Winters' trend and seasonality method, denoted (A,A,A)
   # (assuming that `endog` has a seasonal periodicity of 4, for example if it
   # is quarterly data).
   mod_hw = sm.tsa.statespace.ExponentialSmoothing(endog, trend=True,
                                                   seasonal=4)
   res_hw = mod_hw.fit()

**Differences between Statsmodels' exponential smoothing model classes**

There are several differences between this model class, available at
`sm.tsa.statespace.ExponentialSmoothing`, and the model class available at
`sm.tsa.ExponentialSmoothing`.

- This model class only supports *linear* exponential smoothing models, while
  `sm.tsa.ExponentialSmoothing` also supports multiplicative models.
- This model class puts the exponential smoothing models into state space form
  and then applies the Kalman filter to estimate the states, while
  `sm.tsa.ExponentialSmoothing` is based on exponential smoothing recursions.
  In some cases, this can mean that estimating parameters with this model class
  will be somewhat slower than with `sm.tsa.ExponentialSmoothing`.
- This model class can produce confidence intervals for forecasts, based on an
  assumption of Gaussian errors, while `sm.tsa.ExponentialSmoothing` does not
  support confidence intervals.
- This model class supports concentrating initial values out of the objective
  function, which can improve performance when there are many initial states to
  estimate (for example when the seasonal periodicity is large).
- This model class supports many advanced features available to state space
  models, such as diagnostics and fixed parameters.

**Note**: this class is based on a "multiple sources of error" (MSOE) state
space formulation and not a "single source of error" (SSOE) formulation.

Custom state space models
^^^^^^^^^^^^^^^^^^^^^^^^^

The true power of the state space model is to allow the creation and estimation
of custom models. Usually that is done by extending the following two classes,
which bundle all of state space representation, Kalman filtering, and maximum
likelihood fitting functionality for estimation and results output.

.. autosummary::
   :toctree: generated/

   mlemodel.MLEModel
   mlemodel.MLEResults
   mlemodel.PredictionResults

For a basic example demonstrating creating and estimating a custom state space
model, see the `Local Linear Trend example notebook <examples/notebooks/generated/statespace_local_linear_trend.html>`__.
For a more sophisticated example, see the source code for the `SARIMAX` and
`SARIMAXResults` classes, which are built by extending `MLEModel` and
`MLEResults`.

In simple cases, the model can be constructed entirely using the MLEModel
class. For example, the AR(2) model from above could be constructed and
estimated using only the following code:

.. code-block:: python

   import numpy as np
   from scipy.signal import lfilter
   import statsmodels.api as sm

   # True model parameters
   nobs = int(1e3)
   true_phi = np.r_[0.5, -0.2]
   true_sigma = 1**0.5

   # Simulate a time series
   np.random.seed(1234)
   disturbances = np.random.normal(0, true_sigma, size=(nobs,))
   endog = lfilter([1], np.r_[1, -true_phi], disturbances)

   # Construct the model
   class AR2(sm.tsa.statespace.MLEModel):
       def __init__(self, endog):
           # Initialize the state space model
           super(AR2, self).__init__(endog, k_states=2, k_posdef=1,
                                     initialization='stationary')

           # Setup the fixed components of the state space representation
           self['design'] = [1, 0]
           self['transition'] = [[0, 0],
                                     [1, 0]]
           self['selection', 0, 0] = 1

       # Describe how parameters enter the model
       def update(self, params, transformed=True, **kwargs):
           params = super(AR2, self).update(params, transformed, **kwargs)

           self['transition', 0, :] = params[:2]
           self['state_cov', 0, 0] = params[2]

       # Specify start parameters and parameter names
       @property
       def start_params(self):
           return [0,0,1]  # these are very simple

   # Create and fit the model
   mod = AR2(endog)
   res = mod.fit()
   print(res.summary())

This results in the following summary table::

                              Statespace Model Results                           
   ==============================================================================
   Dep. Variable:                      y   No. Observations:                 1000
   Model:                            AR2   Log Likelihood               -1389.437
   Date:                Wed, 26 Oct 2016   AIC                           2784.874
   Time:                        00:42:03   BIC                           2799.598
   Sample:                             0   HQIC                          2790.470
                                  - 1000                                         
   Covariance Type:                  opg                                         
   ==============================================================================
                    coef    std err          z      P>|z|      [0.025      0.975]
   ------------------------------------------------------------------------------
   param.0        0.4395      0.030     14.730      0.000       0.381       0.498
   param.1       -0.2055      0.032     -6.523      0.000      -0.267      -0.144
   param.2        0.9425      0.042     22.413      0.000       0.860       1.025
   ===================================================================================
   Ljung-Box (Q):                       24.25   Jarque-Bera (JB):                 0.22
   Prob(Q):                              0.98   Prob(JB):                         0.90
   Heteroskedasticity (H):               1.05   Skew:                            -0.04
   Prob(H) (two-sided):                  0.66   Kurtosis:                         3.02
   ===================================================================================
   
   Warnings:
   [1] Covariance matrix calculated using the outer product of gradients (complex-step).

The results object has many of the attributes and methods you would expect from
other statsmodels results objects, including standard errors, z-statistics,
and prediction / forecasting.

More advanced usage is possible, including specifying parameter
transformations, and specifying names for parameters for a more informative
output summary.

Overview of usage
-----------------

All state space models follow the typical Statsmodels pattern:

1. Construct a **model instance** with an input dataset
2. Apply parameters to the model (for example, using `fit`) to construct a **results instance**
3. Interact with the results instance to examine the estimated parameters, explore residual diagnostics, and produce forecasts, simulations, or impulse responses.

An example of this pattern is as follows:

.. code-block:: python

  # Load in the example macroeconomic dataset
  dta = sm.datasets.macrodata.load_pandas().data
  # Make sure we have an index with an associated frequency, so that
  # we can refer to time periods with date strings or timestamps
  dta.index = pd.date_range('1959Q1', '2009Q3', freq='QS')

  # Step 1: construct an SARIMAX model for US inflation data
  model = sm.tsa.SARIMAX(dta.infl, order=(4, 0, 0), trend='c')

  # Step 2: fit the model's parameters by maximum likelihood
  results = model.fit()

  # Step 3: explore / use results

  # - Print a table summarizing estimation results
  print(results.summary())

  # - Print only the estimated parameters
  print(results.params)

  # - Create diagnostic figures based on standardized residuals:
  #   (1) time series graph
  #   (2) histogram
  #   (3) Q-Q plot
  #   (4) correlogram
  results.plot_diagnostics()

  # - Examine diagnostic hypothesis tests
  # Jarque-Bera: [test_statistic, pvalue, skewness, kurtosis]
  print(results.test_normality(method='jarquebera'))
  # Goldfeld-Quandt type test: [test_statistic, pvalue]
  print(results.test_heteroskedasticity(method='breakvar'))
  # Ljung-Box test: [test_statistic, pvalue] for each lag
  print(results.test_serial_correlation(method='ljungbox'))

  # - Forecast the next 4 values
  print(results.forecast(4))

  # - Forecast until 2020Q4
  print(results.forecast('2020Q4'))

  # - Plot in-sample dynamic prediction starting in 2005Q1
  #   and out-of-sample forecasts until 2010Q4 along with
  #   90% confidence intervals
  predict_results = results.get_prediction(start='2005Q1', end='2010Q4', dynamic=True)
  predict_df = predict_results.summary_frame(alpha=0.10)
  fig, ax = plt.subplots()
  predict_df['mean'].plot(ax=ax)
  ax.fill_between(predict_df.index, predict_df['mean_ci_lower'],
                  predict_df['mean_ci_upper'], alpha=0.2)

  # - Simulate two years of new data after the end of the sample
  print(results.simulate(8, anchor='end'))

  # - Impulse responses for two years
  print(results.impulse_responses(8))

Basic methods and attributes for estimation / filtering / smoothing
-------------------------------------------------------------------

The most-used methods for a state space model are:

- :py:meth:`fit <mlemodel.MLEModel.fit>` - estimate parameters via maximum
  likelihood and return a results object (this object will have also performed
  Kalman filtering and smoothing at the estimated parameters). This is the most
  commonly used method.
- :py:meth:`smooth <mlemodel.MLEModel.smooth>` - return a results object
  associated with a given vector of parameters after performing Kalman
  filtering and smoothing
- :py:meth:`loglike <mlemodel.MLEModel.loglike>` - compute the log-likelihood
  of the data using a given vector of parameters

Some useful attributes of a state space model are:

- :py:meth:`param_names <mlemodel.MLEModel.param_names>` - names of the
  parameters used by the model
- :py:meth:`state_names <mlemodel.MLEModel.state_names>` - names of the
  elements of the (unobserved) state vector
- :py:meth:`start_params <mlemodel.MLEModel.start_params>` - initial parameter
  estimates used a starting values for numerical maximum likelihood
  optimization

Other methods that are used less often are:

- :py:meth:`filter <mlemodel.MLEModel.filter>` - return a results object
  associated with a given vector of parameters after only performing Kalman
  filtering (but not smoothing)
- :py:meth:`simulation_smoother <mlemodel.MLEModel.simulation_smoother>` -
  return an object that can perform simulation smoothing

Output and postestimation methods and attributes
------------------------------------------------

Commonly used methods include:

- :py:meth:`summary <mlemodel.MLEResults.summary>` - construct a table that
  presents model fit statistics, estimated parameters, and other summary output
- :py:meth:`predict <mlemodel.MLEResults.predict>` - compute in-sample
  predictions and out-of-sample forecasts (point estimates only)
- :py:meth:`get_prediction <mlemodel.MLEResults.get_prediction>` - compute
  in-sample predictions and out-of-sample forecasts, including confidence
  intervals
- :py:meth:`forecast <mlemodel.MLEResults.forecast>` - compute out-of-sample
  forecasts (point estimates only) (this is a convenience wrapper around
  `predict`)
- :py:meth:`get_forecast <mlemodel.MLEResults.get_forecast>` - compute
  out-of-sample forecasts, including confidence intervals (this is a
  convenience wrapper around `get_prediction`)
- :py:meth:`simulate <mlemodel.MLEResults.simulate>` - simulate new data
  according to the state space model
- :py:meth:`impulse_responses <mlemodel.MLEResults.impulse_responses>` -
  compute impulse responses from the state space model

Commonly used attributes include:

- :py:meth:`params <mlemodel.MLEResults.params>` - estimated parameters
- :py:meth:`bse <mlemodel.MLEResults.bse>` - standard errors of estimated
  parameters
- :py:meth:`pvalues <mlemodel.MLEResults.pvalues>` - p-values associated with
  estimated parameters
- :py:meth:`llf <mlemodel.MLEResults.llf>` - log-likelihood of the data at
  the estimated parameters
- :py:meth:`sse <mlemodel.MLEResults.sse>`,
  :py:meth:`mse <mlemodel.MLEResults.mse>`, and
  :py:meth:`mae <mlemodel.MLEResults.mae>` - sum of squared errors,
  mean square error, and mean absolute error
- Information criteria, including: :py:meth:`aic <mlemodel.MLEResults.aic>`,
  :py:meth:`aicc <mlemodel.MLEResults.aicc>`,
  :py:meth:`bic <mlemodel.MLEResults.bic>`, and
  :py:meth:`hquc <mlemodel.MLEResults.hqic>`
- :py:meth:`fittedvalues <mlemodel.MLEResults.fittedvalues>` - fitted values
  from the model (note that these are one-step-ahead predictions)
- :py:meth:`resid <mlemodel.MLEResults.resid>` - residuals from the model (note
  that these are one-step-ahead prediction errors)

Estimates and covariances of the unobserved state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It can be useful to compute estimates of the unobserved state vector
conditional on the observed data. These are available in the results object
:py:meth:`states <mlemodel.MLEResults.states>`, which contains the following
elements:

- `states.filtered` - filtered (one-sided) estimates of the state vector. The
  estimate of the state vector at time `t` is based on the observed data up
  to and including time `t`.
- `states.smoothed` - smoothed (two-sided) estimates of the state vector. The
  estimate of the state vector at time `t` is based on all observed data in
  the sample.
- `states.filtered_cov` - filtered (one-sided) covariance of the state vector
- `states.smoothed_cov` - smoothed (two-sided) covariance of the state vector

Each of these elements are Pandas `DataFrame` objects.

As an example, in a "local level + seasonal" model estimated via the
`UnobservedComponents` components class we can get an estimates of the
underlying level and seasonal movements of a series over time.

.. code-block:: python

  fig, axes = plt.subplots(3, 1, figsize=(8, 8))

  # Retrieve monthly retail sales for clothing
  from pandas_datareader.data import DataReader
  clothing = DataReader('MRTSSM4481USN', 'fred', start='1992').asfreq('MS')['MRTSSM4481USN']

  # Construct a local level + seasonal model
  model = sm.tsa.UnobservedComponents(clothing, 'llevel', seasonal=12)
  results = model.fit()

  # Plot the data, the level, and seasonal
  clothing.plot(ax=axes[0])
  results.states.smoothed['level'].plot(ax=axes[1])
  results.states.smoothed['seasonal'].plot(ax=axes[2])

Residual diagnostics
^^^^^^^^^^^^^^^^^^^^

Three diagnostic tests are available after estimation of any statespace model,
whether built in or custom, to help assess whether the model conforms to the
underlying statistical assumptions. These tests are:

- :py:meth:`test_normality <mlemodel.MLEResults.test_normality>`
- :py:meth:`test_heteroskedasticity <mlemodel.MLEResults.test_heteroskedasticity>`
- :py:meth:`test_serial_correlation <mlemodel.MLEResults.test_serial_correlation>`

A number of standard plots of regression residuals are available for the same
purpose. These can be produced using the command
:py:meth:`plot_diagnostics <mlemodel.MLEResults.plot_diagnostics>`.

Applying estimated parameters to an updated or different dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three methods that can be used to apply estimated parameters from a
results object to an updated or different dataset:

- :py:meth:`append <mlemodel.MLEResults.append>` - retrieve a new results
  object with additional observations that follow after the end of the current
  sample appended to it (so the new results object contains both the current
  sample and the additional observations)
- :py:meth:`extend <mlemodel.MLEResults.extend>` - retrieve a new results
  object for additional observations that follow after end of the current
  sample (so the new results object contains only the new observations but NOT
  the current sample)
- :py:meth:`apply <mlemodel.MLEResults.apply>` - retrieve a new results object
  for a completely different dataset

One cross-validation exercise on time-series data involves fitting a model's
parameters based on a training sample (observations through time `t`) and
then evaluating the fit of the model using a test sample (observations `t+1`,
`t+2`, ...). This can be conveniently done using either `apply` or `extend`. In
the example below, we use the `extend` method.

.. code-block:: python

  # Load in the example macroeconomic dataset
  dta = sm.datasets.macrodata.load_pandas().data
  # Make sure we have an index with an associated frequency, so that
  # we can refer to time periods with date strings or timestamps
  dta.index = pd.date_range('1959Q1', '2009Q3', freq='QS')

  # Separate inflation data into a training and test dataset
  training_endog = dta['infl'].iloc[:-1]
  test_endog = dta['infl'].iloc[-1:]

  # Fit an SARIMAX model for inflation
  training_model = sm.tsa.SARIMAX(training_endog, order=(4, 0, 0))
  training_results = training_model.fit()

  # Extend the results to the test observations
  test_results = training_results.extend(test_endog)

  # Print the sum of squared errors in the test sample,
  # based on parameters computed using only the training sample
  print(test_results.sse)


Understanding the Impact of Data Revisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Statespace model results expose a :meth:`~mlemodel.MLEModel.news` method that
can be used to understand the impact of data revisions -- news -- on model
parameters.

.. autosummary::
   :toctree: generated/

   news.NewsResults


Additional options and tools
----------------------------

All state space models have the following options and tools:

Holding some parameters fixed and estimating the rest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :py:meth:`fit_constrained <mlemodel.MLEModel.fit_constrained>` method
allows fixing some parameters to known values and then estimating the rest via
maximum likelihood. An example of this is:

.. code-block:: python

  # Construct a model
  model = sm.tsa.SARIMAX(endog, order=(1, 0, 0))

  # To find out the parameter names, use:
  print(model.param_names)

  # Fit the model with a fixed value for the AR(1) coefficient:
  results = model.fit_constrained({'ar.L1': 0.5})

Alternatively, you can use the
:py:meth:`fix_params <mlemodel.MLEModel.fix_params>` context manager:

.. code-block:: python

  # Construct a model
  model = sm.tsa.SARIMAX(endog, order=(1, 0, 0))

  # Fit the model with a fixed value for the AR(1) coefficient using the
  # context manager
  with model.fix_params({'ar.L1': 0.5}):
      results = model.fit()

Low memory options
^^^^^^^^^^^^^^^^^^

When the observed dataset is very large and / or the state vector of the model
is high-dimensional (for example when considering long seasonal effects), the
default memory requirements can be too large. For this reason, the `fit`,
`filter`, and `smooth` methods accept an optional `low_memory=True` argument,
which can considerably reduce memory requirements and speed up model fitting.

Note that when using `low_memory=True`, not all results objects will be
available. However, residual diagnostics, in-sample (non-dynamic) prediction,
and out-of-sample forecasting are all still available.

Low-level state space representation and Kalman filtering
---------------------------------------------------------

While creation of custom models will almost always be done by extending
`MLEModel` and `MLEResults`, it can be useful to understand the superstructure
behind those classes.

Maximum likelihood estimation requires evaluating the likelihood function of
the model, and for models in state space form the likelihood function is
evaluated as a byproduct of running the Kalman filter.

There are two classes used by `MLEModel` that facilitate specification of the
state space model and Kalman filtering: `Representation` and `KalmanFilter`.

The `Representation` class is the piece where the state space model
representation is defined. In simple terms, it holds the state space matrices
(`design`, `obs_intercept`, etc.; see the introduction to state space models,
above) and allows their manipulation.

`FrozenRepresentation` is the most basic results-type class, in that it takes a
"snapshot" of the state space representation at any given time. See the class
documentation for the full list of available attributes.

.. autosummary::
   :toctree: generated/

   representation.Representation
   representation.FrozenRepresentation

The `KalmanFilter` class is a subclass of Representation that provides
filtering capabilities. Once the state space representation matrices have been
constructed, the :py:meth:`filter <kalman_filter.KalmanFilter.filter>`
method can be called, producing a `FilterResults` instance; `FilterResults` is
a subclass of `FrozenRepresentation`.

The `FilterResults` class not only holds a frozen representation of the state
space model (the design, transition, etc. matrices, as well as model
dimensions, etc.) but it also holds the filtering output, including the
:py:attr:`filtered state <kalman_filter.FilterResults.filtered_state>` and
loglikelihood (see the class documentation for the full list of available
results). It also provides a
:py:meth:`predict <kalman_filter.FilterResults.predict>` method, which allows
in-sample prediction or out-of-sample forecasting. A similar method,
:py:meth:`predict <kalman_filter.FilterResults.get_prediction>`, provides
additional prediction or forecasting results, including confidence intervals.

.. autosummary::
   :toctree: generated/

   kalman_filter.KalmanFilter
   kalman_filter.FilterResults
   kalman_filter.PredictionResults

The `KalmanSmoother` class is a subclass of `KalmanFilter` that provides
smoothing capabilities. Once the state space representation matrices have been
constructed, the :py:meth:`filter <kalman_smoother.KalmanSmoother.smooth>`
method can be called, producing a `SmootherResults` instance; `SmootherResults`
is a subclass of `FilterResults`.

The `SmootherResults` class holds all the output from `FilterResults`, but
also includes smoothing output, including the
:py:attr:`smoothed state <kalman_filter.SmootherResults.smoothed_state>` and
loglikelihood (see the class documentation for the full list of available
results). Whereas "filtered" output at time `t` refers to estimates conditional
on observations up through time `t`, "smoothed" output refers to estimates
conditional on the entire set of observations in the dataset.

.. autosummary::
   :toctree: generated/

   kalman_smoother.KalmanSmoother
   kalman_smoother.SmootherResults

The `SimulationSmoother` class is a subclass of `KalmanSmoother` that further
provides simulation and simulation smoothing capabilities. The
:py:meth:`simulation_smoother <simulation_smoother.SimulationSmoother.simulation_smoother>`
method can be called, producing a `SimulationSmoothResults` instance.

The `SimulationSmoothResults` class has a `simulate` method, that allows
performing simulation smoothing to draw from the joint posterior of the state
vector. This is useful for Bayesian estimation of state space models via Gibbs
sampling.

.. autosummary::
   :toctree: generated/

   simulation_smoother.SimulationSmoother
   simulation_smoother.SimulationSmoothResults
   cfa_simulation_smoother.CFASimulationSmoother

Statespace Tools
----------------

There are a variety of tools used for state space modeling or by the SARIMAX
class:

.. autosummary::
   :toctree: generated/

   tools.companion_matrix
   tools.diff
   tools.is_invertible
   tools.constrain_stationary_univariate
   tools.unconstrain_stationary_univariate
   tools.constrain_stationary_multivariate
   tools.unconstrain_stationary_multivariate
   tools.validate_matrix_shape
   tools.validate_vector_shape
.. _formula_examples:

Fitting models using R-style formulas
=====================================

Since version 0.5.0, ``statsmodels`` allows users to fit statistical
models using R-style formulas. Internally, ``statsmodels`` uses the
`patsy <https://patsy.readthedocs.io/en/latest/>`_ package to convert formulas and
data to the matrices that are used in model fitting. The formula
framework is quite powerful; this tutorial only scratches the surface. A
full description of the formula language can be found in the ``patsy``
docs:

-  `Patsy formula language description <https://patsy.readthedocs.io/en/latest/>`_

Loading modules and functions
-----------------------------

.. ipython:: python

    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    import numpy as np
    import pandas

Notice that we called ``statsmodels.formula.api`` in addition to the usual
``statsmodels.api``. In fact, ``statsmodels.api`` is used here only to load
the dataset. The ``formula.api`` hosts many of the same
functions found in ``api`` (e.g. OLS, GLM), but it also holds lower case
counterparts for most of these models. In general, lower case models
accept ``formula`` and ``df`` arguments, whereas upper case ones take
``endog`` and ``exog`` design matrices. ``formula`` accepts a string
which describes the model in terms of a ``patsy`` formula. ``df`` takes
a `pandas <https://pandas.pydata.org/>`_ data frame.

``dir(smf)`` will print a list of available models.

Formula-compatible models have the following generic call signature:
``(formula, data, subset=None, *args, **kwargs)``

OLS regression using formulas
-----------------------------

To begin, we fit the linear model described on the `Getting
Started <gettingstarted.html>`_ page. Download the data, subset columns,
and list-wise delete to remove missing observations:

.. ipython:: python

    df = sm.datasets.get_rdataset("Guerry", "HistData").data
    df = df[['Lottery', 'Literacy', 'Wealth', 'Region']].dropna()
    df.head()

Fit the model:

.. ipython:: python

    mod = smf.ols(formula='Lottery ~ Literacy + Wealth + Region', data=df)
    res = mod.fit()
    print(res.summary())

Categorical variables
---------------------

Looking at the summary printed above, notice that ``patsy`` determined
that elements of *Region* were text strings, so it treated *Region* as a
categorical variable. ``patsy``'s default is also to include an
intercept, so we automatically dropped one of the *Region* categories.

If *Region* had been an integer variable that we wanted to treat
explicitly as categorical, we could have done so by using the ``C()``
operator:

.. ipython:: python

    res = smf.ols(formula='Lottery ~ Literacy + Wealth + C(Region)', data=df).fit()
    print(res.params)


Examples more advanced features ``patsy``'s categorical variables
function can be found here: `Patsy: Contrast Coding Systems for
categorical variables <contrasts.html>`_

Operators
---------

We have already seen that "~" separates the left-hand side of the model
from the right-hand side, and that "+" adds new columns to the design
matrix.

Removing variables
~~~~~~~~~~~~~~~~~~

The "-" sign can be used to remove columns/variables. For instance, we
can remove the intercept from a model by:

.. ipython:: python

    res = smf.ols(formula='Lottery ~ Literacy + Wealth + C(Region) -1 ', data=df).fit()
    print(res.params)


Multiplicative interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

":" adds a new column to the design matrix with the product of the other
two columns. "\*" will also include the individual columns that were
multiplied together:

.. ipython:: python

    res1 = smf.ols(formula='Lottery ~ Literacy : Wealth - 1', data=df).fit()
    res2 = smf.ols(formula='Lottery ~ Literacy * Wealth - 1', data=df).fit()
    print(res1.params)
    print(res2.params)


Many other things are possible with operators. Please consult the `patsy
docs <https://patsy.readthedocs.io/en/latest/formulas.html>`_ to learn
more.

Functions
---------

You can apply vectorized functions to the variables in your model:

.. ipython:: python

    res = smf.ols(formula='Lottery ~ np.log(Literacy)', data=df).fit()
    print(res.params)


Define a custom function:

.. ipython:: python

    def log_plus_1(x):
        return np.log(x) + 1.0

    res = smf.ols(formula='Lottery ~ log_plus_1(Literacy)', data=df).fit()
    print(res.params)

.. _patsy-namespaces:

Namespaces
----------

Notice that all of the above examples use the calling namespace to look for the functions to apply. The namespace used can be controlled via the ``eval_env`` keyword. For example, you may want to give a custom namespace using the :class:`patsy:patsy.EvalEnvironment` or you may want to use a "clean" namespace, which we provide by passing ``eval_func=-1``. The default is to use the caller's namespace. This can have (un)expected consequences, if, for example, someone has a variable names ``C`` in the user namespace or in their data structure passed to ``patsy``, and ``C`` is used in the formula to handle a categorical variable. See the `Patsy API Reference <https://patsy.readthedocs.io/en/latest/API-reference.html>`_ for more information.

Using formulas with models that do not (yet) support them
---------------------------------------------------------

Even if a given ``statsmodels`` function does not support formulas, you
can still use ``patsy``'s formula language to produce design matrices.
Those matrices can then be fed to the fitting function as ``endog`` and
``exog`` arguments.

To generate ``numpy`` arrays:

.. ipython:: python

    import patsy
    f = 'Lottery ~ Literacy * Wealth'
    y, X = patsy.dmatrices(f, df, return_type='matrix')
    print(y[:5])
    print(X[:5])

``y`` and ``X`` would be instances of ``patsy.DesignMatrix`` which is a subclass of ``numpy.ndarray``.

To generate pandas data frames:

.. ipython:: python

    f = 'Lottery ~ Literacy * Wealth'
    y, X = patsy.dmatrices(f, df, return_type='dataframe')
    print(y[:5])
    print(X[:5])

.. ipython:: python

    print(sm.OLS(y, X).fit().summary())
.. module:: statsmodels.duration
   :synopsis: Models for durations

.. currentmodule:: statsmodels.duration


.. _duration:

Methods for Survival and Duration Analysis
==========================================

:mod:`statsmodels.duration` implements several standard methods for
working with censored data.  These methods are most commonly used when
the data consist of durations between an origin time point and the
time at which some event of interest occurred.  A typical example is a
medical study in which the origin is the time at which a subject is
diagnosed with some condition, and the event of interest is death (or
disease progression, recovery, etc.).

Currently only right-censoring is handled.  Right censoring occurs
when we know that an event occurred after a given time `t`, but we do
not know the exact event time.

Survival function estimation and inference
------------------------------------------

The :class:`statsmodels.api.SurvfuncRight` class can be used to
estimate a survival function using data that may be right censored.
``SurvfuncRight`` implements several inference procedures including
confidence intervals for survival distribution quantiles, pointwise
and simultaneous confidence bands for the survival function, and
plotting procedures.  The ``duration.survdiff`` function provides
testing procedures for comparing survival distributions.

Here we create a ``SurvfuncRight`` object using data from the
`flchain` study, which is available through the R datasets repository.
We fit the survival distribution only for the female subjects.


.. code-block:: python

   import statsmodels.api as sm

   data = sm.datasets.get_rdataset("flchain", "survival").data
   df = data.loc[data.sex == "F", :]
   sf = sm.SurvfuncRight(df["futime"], df["death"])

The main features of the fitted survival distribution can be seen by
calling the ``summary`` method:

.. code-block:: python

    sf.summary().head()

We can obtain point estimates and confidence intervals for quantiles
of the survival distribution.  Since only around 30% of the subjects
died during this study, we can only estimate quantiles below the 0.3
probability point:

.. code-block:: python

    sf.quantile(0.25)
    sf.quantile_ci(0.25)

To plot a single survival function, call the ``plot`` method:

.. code-block:: python

    sf.plot()

Since this is a large dataset with a lot of censoring, we may wish
to not plot the censoring symbols:

.. code-block:: python

    fig = sf.plot()
    ax = fig.get_axes()[0]
    pt = ax.get_lines()[1]
    pt.set_visible(False)

We can also add a 95% simultaneous confidence band to the plot.
Typically these bands only plotted for central part of the
distribution.

.. code-block:: python

    fig = sf.plot()
    lcb, ucb = sf.simultaneous_cb()
    ax = fig.get_axes()[0]
    ax.fill_between(sf.surv_times, lcb, ucb, color='lightgrey')
    ax.set_xlim(365, 365*10)
    ax.set_ylim(0.7, 1)
    ax.set_ylabel("Proportion alive")
    ax.set_xlabel("Days since enrollment")

Here we plot survival functions for two groups (females and males) on
the same axes:

.. code-block:: python

    gb = data.groupby("sex")
    ax = plt.axes()
    sexes = []
    for g in gb:
        sexes.append(g[0])
        sf = sm.SurvfuncRight(g[1]["futime"], g[1]["death"])
        sf.plot(ax)
    li = ax.get_lines()
    li[1].set_visible(False)
    li[3].set_visible(False)
    plt.figlegend((li[0], li[2]), sexes, "center right")
    plt.ylim(0.6, 1)
    ax.set_ylabel("Proportion alive")
    ax.set_xlabel("Days since enrollment")

We can formally compare two survival distributions with ``survdiff``,
which implements several standard nonparametric procedures.  The
default procedure is the logrank test:

.. code-block:: python

    stat, pv = sm.duration.survdiff(data.futime, data.death, data.sex)

Here are some of the other testing procedures implemented by survdiff:

.. code-block:: python

    # Fleming-Harrington with p=1, i.e. weight by pooled survival time
    stat, pv = sm.duration.survdiff(data.futime, data.death, data.sex, weight_type='fh', fh_p=1)

    # Gehan-Breslow, weight by number at risk
    stat, pv = sm.duration.survdiff(data.futime, data.death, data.sex, weight_type='gb')

    # Tarone-Ware, weight by the square root of the number at risk
    stat, pv = sm.duration.survdiff(data.futime, data.death, data.sex, weight_type='tw')


Regression methods
------------------

Proportional hazard regression models ("Cox models") are a regression
technique for censored data.  They allow variation in the time to an
event to be explained in terms of covariates, similar to what is done
in a linear or generalized linear regression model.  These models
express the covariate effects in terms of "hazard ratios", meaning the
the hazard (instantaneous event rate) is multiplied by a given factor
depending on the value of the covariates.


.. code-block:: python

   import statsmodels.api as sm
   import statsmodels.formula.api as smf

   data = sm.datasets.get_rdataset("flchain", "survival").data
   del data["chapter"]
   data = data.dropna()
   data["lam"] = data["lambda"]
   data["female"] = (data["sex"] == "F").astype(int)
   data["year"] = data["sample.yr"] - min(data["sample.yr"])
   status = data["death"].values

   mod = smf.phreg("futime ~ 0 + age + female + creatinine + "
                   "np.sqrt(kappa) + np.sqrt(lam) + year + mgus",
                   data, status=status, ties="efron")
   rslt = mod.fit()
   print(rslt.summary())


See :ref:`statsmodels-examples` for more detailed examples.


There are some notebook examples on the Wiki:
`Wiki notebooks for PHReg and Survival Analysis <https://github.com/statsmodels/statsmodels/wiki/Examples#survival-analysis>`_


.. todo::

   Technical Documentation

References
^^^^^^^^^^

References for Cox proportional hazards regression model::

    T Therneau (1996). Extending the Cox model. Technical report.
    http://www.mayo.edu/research/documents/biostat-58pdf/DOC-10027288

    G Rodriguez (2005). Non-parametric estimation in survival models.
    http://data.princeton.edu/pop509/NonParametricSurvival.pdf

    B Gillespie (2006). Checking the assumptions in the Cox proportional
    hazards model.
    http://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf


Module Reference
----------------

.. module:: statsmodels.duration.survfunc
   :synopsis: Models for Survival Analysis

.. currentmodule:: statsmodels.duration.survfunc

The class for working with survival distributions is:

.. autosummary::
   :toctree: generated/

   SurvfuncRight

.. module:: statsmodels.duration.hazard_regression
   :synopsis: Proportional hazards model for Survival Analysis

.. currentmodule:: statsmodels.duration.hazard_regression

The proportional hazards regression model class is:

.. autosummary::
   :toctree: generated/

   PHReg

The proportional hazards regression result class is:

.. autosummary::
   :toctree: generated/

   PHRegResults

The primary helper class is:

.. autosummary::
   :toctree: generated/

   rv_discrete_float:orphan:

.. _missing_data:

Missing Data
------------
All of the models can handle missing data. For performance reasons, the default is not to do any checking for missing data. If, however, you would like for missing data to be handled internally, you can do so by using the missing keyword argument. The default is to do nothing

.. ipython:: python

   import statsmodels.api as sm
   data = sm.datasets.longley.load()
   data.exog = sm.add_constant(data.exog)
   # add in some missing data
   missing_idx = np.array([False] * len(data.endog))
   missing_idx[[4, 10, 15]] = True
   data.endog[missing_idx] = np.nan
   ols_model = sm.OLS(data.endog, data.exog)
   ols_fit = ols_model.fit()
   print(ols_fit.params)

This silently fails and all of the model parameters are NaN, which is probably not what you expected. If you are not sure whether or not you have missing data you can use `missing = 'raise'`. This will raise a `MissingDataError` during model instantiation if missing data is present so that you know something was wrong in your input data.

.. ipython:: python
   :okexcept:

   ols_model = sm.OLS(data.endog, data.exog, missing='raise')

If you want statsmodels to handle the missing data by dropping the observations, use `missing = 'drop'`.

.. ipython:: python

   ols_model = sm.OLS(data.endog, data.exog, missing='drop')

We are considering adding a configuration framework so that you can set the option with a global setting.
.. currentmodule:: statsmodels.genmod.generalized_estimating_equations

.. _gee:

Generalized Estimating Equations
================================

Generalized Estimating Equations estimate generalized linear models for
panel, cluster or repeated measures data when the observations are possibly
correlated withing a cluster but uncorrelated across clusters. It supports
estimation of the same one-parameter exponential families as Generalized
Linear models (`GLM`).

See `Module Reference`_ for commands and arguments.

Examples
--------

The following illustrates a Poisson regression with exchangeable correlation
within clusters using data on epilepsy seizures.

.. ipython:: python

    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    data = sm.datasets.get_rdataset('epil', package='MASS').data

    fam = sm.families.Poisson()
    ind = sm.cov_struct.Exchangeable()
    mod = smf.gee("y ~ age + trt + base", "subject", data,
                  cov_struct=ind, family=fam)
    res = mod.fit()
    print(res.summary())


Several notebook examples of the use of GEE can be found on the Wiki:
`Wiki notebooks for GEE <https://github.com/statsmodels/statsmodels/wiki/Examples#generalized-estimating-equations-gee>`_


References
^^^^^^^^^^

* KY Liang and S Zeger. "Longitudinal data analysis using generalized
  linear models". Biometrika (1986) 73 (1): 13-22.
* S Zeger and KY Liang. "Longitudinal Data Analysis for Discrete and
  Continuous Outcomes". Biometrics Vol. 42, No. 1 (Mar., 1986),
  pp. 121-130
* A Rotnitzky and NP Jewell (1990). "Hypothesis testing of regression
  parameters in semiparametric generalized linear models for cluster
  correlated data", Biometrika, 77, 485-497.
* Xu Guo and Wei Pan (2002). "Small sample performance of the score test in
  GEE".
  http://www.sph.umn.edu/faculty1/wp-content/uploads/2012/11/rr2002-013.pdf
* LA Mancl LA, TA DeRouen (2001). A covariance estimator for GEE with improved
  small-sample properties.  Biometrics. 2001 Mar;57(1):126-34.


Module Reference
----------------

.. module:: statsmodels.genmod.generalized_estimating_equations
   :synopsis: Generalized estimating equations

Model Class
^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   GEE
   NominalGEE
   OrdinalGEE

.. module:: statsmodels.genmod.qif
   :synopsis: Quadratic inference functions

.. currentmodule:: statsmodels.genmod.qif

.. autosummary::
   :toctree: generated/

   QIF

Results Classes
^^^^^^^^^^^^^^^

.. currentmodule:: statsmodels.genmod.generalized_estimating_equations

.. autosummary::
   :toctree: generated/

   GEEResults
   GEEMargins

.. currentmodule:: statsmodels.genmod.qif

.. autosummary::
   :toctree: generated/

   QIFResults

Dependence Structures
^^^^^^^^^^^^^^^^^^^^^

The dependence structures currently implemented are

.. module:: statsmodels.genmod.cov_struct
   :synopsis: Covariance structures for Generalized Estimating Equations (GEE)

.. currentmodule:: statsmodels.genmod.cov_struct

.. autosummary::
   :toctree: generated/

   CovStruct
   Autoregressive
   Exchangeable
   GlobalOddsRatio
   Independence
   Nested


Families
^^^^^^^^

The distribution families are the same as for GLM, currently implemented are

.. module:: statsmodels.genmod.families.family
   :synopsis: Generalized Linear Model (GLM) families

.. currentmodule:: statsmodels.genmod.families.family

.. autosummary::
   :toctree: generated/

   Family
   Binomial
   Gamma
   Gaussian
   InverseGaussian
   NegativeBinomial
   Poisson
   Tweedie


Link Functions
^^^^^^^^^^^^^^

The link functions are the same as for GLM, currently implemented are the
following. Not all link functions are available for each distribution family.
The list of available link functions can be obtained by

::

    >>> sm.families.family.<familyname>.links

.. currentmodule:: statsmodels.genmod.families.links

.. autosummary::
   :toctree: generated/

   Link

   CDFLink
   CLogLog
   Log
   Logit
   NegativeBinomial
   Power
   cauchy
   cloglog
   identity
   inverse_power
   inverse_squared
   log
   logit
   nbinom
   probit
.. _endog_exog:

``endog``, ``exog``, what's that?
=================================

statsmodels is using ``endog`` and ``exog`` as names for the data, the
observed variables that are used in an estimation problem. Other names that
are often used in different statistical packages or text books are, for
example,

===================== ======================
endog                 exog
===================== ======================
y                     x
y variable            x variable
left hand side (LHS)  right hand side (RHS)
dependent variable    independent variable
regressand            regressors
outcome               design
response variable     explanatory variable
===================== ======================

The usage is quite often domain and model specific; however, we have chosen
to use `endog` and `exog` almost exclusively. A mnemonic hint to keep the two
terms apart is that exogenous has an "x", as in x-variable, in its name.

`x` and `y` are one letter names that are sometimes used for temporary
variables and are not informative in itself. To avoid one letter names we
decided to use descriptive names and settled on ``endog`` and ``exog``.
Since this has been criticized, this might change in future.

Background
----------

Some informal definitions of the terms are

`endogenous`: caused by factors within the system

`exogenous`: caused by factors outside the system

*Endogenous variables designates variables in an economic/econometric model
that are explained, or predicted, by that model.*
http://stats.oecd.org/glossary/detail.asp?ID=794

*Exogenous variables designates variables that appear in an
economic/econometric model, but are not explained by that model (i.e. they are
taken as given by the model).*  http://stats.oecd.org/glossary/detail.asp?ID=890

In econometrics and statistics the terms are defined more formally, and
different definitions of exogeneity (weak, strong, strict) are used depending
on the model. The usage in statsmodels as variable names cannot always be
interpreted in a formal sense, but tries to follow the same principle.


In the simplest form, a model relates an observed variable, y, to another set
of variables, x, in some linear or nonlinear form ::

   y = f(x, beta) + noise
   y = x * beta + noise

However, to have a statistical model we need additional assumptions on the
properties of the explanatory variables, x, and the noise. One standard
assumption for many basic models is that x is not correlated with the noise.
In a more general definition, x being exogenous means that we do not have to
consider how the explanatory variables in x were generated, whether by design
or by random draws from some underlying distribution, when we want to estimate
the effect or impact that x has on y, or test a hypothesis about this effect.

In other words, y is *endogenous* to our model, x is *exogenous* to our model
for the estimation.

As an example, suppose you run an experiment and for the second session some
subjects are not available anymore.
Is the drop-out relevant for the conclusions you draw for the experiment?
In other words, can we treat the drop-out decision as exogenous for our
problem.

It is up to the user to know (or to consult a text book to find out) what the
underlying statistical assumptions for the models are. As an example, ``exog``
in ``OLS`` can have lagged dependent variables if the error or noise term is
independently distributed over time (or uncorrelated over time). However, if
the error terms are autocorrelated, then OLS does not have good statistical
properties (is inconsistent) and the correct model will be ARMAX.
``statsmodels`` has functions for regression diagnostics to test whether some of
the assumptions are justified or not.

.. currentmodule:: statsmodels.stats.contingency_tables

.. _contingency_tables:


Contingency tables
==================

statsmodels supports a variety of approaches for analyzing contingency
tables, including methods for assessing independence, symmetry,
homogeneity, and methods for working with collections of tables from a
stratified population.

The methods described here are mainly for two-way tables.  Multi-way
tables can be analyzed using log-linear models.  statsmodels does not
currently have a dedicated API for loglinear modeling, but Poisson
regression in :class:`statsmodels.genmod.GLM` can be used for this
purpose.

A contingency table is a multi-way table that describes a data set in
which each observation belongs to one category for each of several
variables.  For example, if there are two variables, one with
:math:`r` levels and one with :math:`c` levels, then we have a
:math:`r \times c` contingency table.  The table can be described in
terms of the number of observations that fall into a given cell of the
table, e.g. :math:`T_{ij}` is the number of observations that have
level :math:`i` for the first variable and level :math:`j` for the
second variable.  Note that each variable must have a finite number of
levels (or categories), which can be either ordered or unordered.  In
different contexts, the variables defining the axes of a contingency
table may be called **categorical variables** or **factor variables**.
They may be either **nominal** (if their levels are unordered) or
**ordinal** (if their levels are ordered).

The underlying population for a contingency table is described by a
**distribution table** :math:`P_{i, j}`.  The elements of :math:`P`
are probabilities, and the sum of all elements in :math:`P` is 1.
Methods for analyzing contingency tables use the data in :math:`T` to
learn about properties of :math:`P`.

The :class:`statsmodels.stats.Table` is the most basic class for
working with contingency tables.  We can create a ``Table`` object
directly from any rectangular array-like object containing the
contingency table cell counts:

.. ipython:: python

    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    df = sm.datasets.get_rdataset("Arthritis", "vcd").data

    tab = pd.crosstab(df['Treatment'], df['Improved'])
    tab = tab.loc[:, ["None", "Some", "Marked"]]
    table = sm.stats.Table(tab)

Alternatively, we can pass the raw data and let the Table class
construct the array of cell counts for us:

.. ipython:: python

    data = df[["Treatment", "Improved"]]
    table = sm.stats.Table.from_data(data)


Independence
------------

**Independence** is the property that the row and column factors occur
independently. **Association** is the lack of independence.  If the
joint distribution is independent, it can be written as the outer
product of the row and column marginal distributions:

.. math::

    P_{ij} = \sum_k P_{ij} \cdot \sum_k P_{kj} \quad \text{for all} \quad  i, j

We can obtain the best-fitting independent distribution for our
observed data, and then view residuals which identify particular cells
that most strongly violate independence:

.. ipython:: python

    print(table.table_orig)
    print(table.fittedvalues)
    print(table.resid_pearson)

In this example, compared to a sample from a population in which the
rows and columns are independent, we have too many observations in the
placebo/no improvement and treatment/marked improvement cells, and too
few observations in the placebo/marked improvement and treated/no
improvement cells.  This reflects the apparent benefits of the
treatment.

If the rows and columns of a table are unordered (i.e. are nominal
factors), then the most common approach for formally assessing
independence is using Pearson's :math:`\chi^2` statistic.  It's often
useful to look at the cell-wise contributions to the :math:`\chi^2`
statistic to see where the evidence for dependence is coming from.

.. ipython:: python

    rslt = table.test_nominal_association()
    print(rslt.pvalue)
    print(table.chi2_contribs)

For tables with ordered row and column factors, we can us the **linear
by linear** association test to obtain more power against alternative
hypotheses that respect the ordering.  The test statistic for the
linear by linear association test is

.. math::

    \sum_k r_i c_j T_{ij}

where :math:`r_i` and :math:`c_j` are row and column scores.  Often
these scores are set to the sequences 0, 1, ....  This gives the
'Cochran-Armitage trend test'.

.. ipython:: python

    rslt = table.test_ordinal_association()
    print(rslt.pvalue)

We can assess the association in a :math:`r\times x` table by
constructing a series of :math:`2\times 2` tables and calculating
their odds ratios.  There are two ways to do this.  The **local odds
ratios** construct :math:`2\times 2` tables from adjacent row and
column categories.

.. ipython:: python

    print(table.local_oddsratios)
    taloc = sm.stats.Table2x2(np.asarray([[7, 29], [21, 13]]))
    print(taloc.oddsratio)
    taloc = sm.stats.Table2x2(np.asarray([[29, 7], [13, 7]]))
    print(taloc.oddsratio)

The **cumulative odds ratios** construct :math:`2\times 2` tables by
dichotomizing the row and column factors at each possible point.

.. ipython:: python

    print(table.cumulative_oddsratios)
    tab1 = np.asarray([[7, 29 + 7], [21, 13 + 7]])
    tacum = sm.stats.Table2x2(tab1)
    print(tacum.oddsratio)
    tab1 = np.asarray([[7 + 29, 7], [21 + 13, 7]])
    tacum = sm.stats.Table2x2(tab1)
    print(tacum.oddsratio)

A mosaic plot is a graphical approach to informally assessing
dependence in two-way tables.

.. ipython:: python

    from statsmodels.graphics.mosaicplot import mosaic
    fig, _ = mosaic(data, index=["Treatment", "Improved"])


Symmetry and homogeneity
------------------------

**Symmetry** is the property that :math:`P_{i, j} = P_{j, i}` for
every :math:`i` and :math:`j`.  **Homogeneity** is the property that
the marginal distribution of the row factor and the column factor are
identical, meaning that

.. math::

    \sum_j P_{ij} = \sum_j P_{ji} \forall i

Note that for these properties to be applicable the table :math:`P`
(and :math:`T`) must be square, and the row and column categories must
be identical and must occur in the same order.

To illustrate, we load a data set, create a contingency table, and
calculate the row and column margins.  The :class:`Table` class
contains methods for analyzing :math:`r \times c` contingency tables.
The data set loaded below contains assessments of visual acuity in
people's left and right eyes.  We first load the data and create a
contingency table.

.. ipython:: python

    df = sm.datasets.get_rdataset("VisualAcuity", "vcd").data
    df = df.loc[df.gender == "female", :]
    tab = df.set_index(['left', 'right'])
    del tab["gender"]
    tab = tab.unstack()
    tab.columns = tab.columns.get_level_values(1)
    print(tab)

Next we create a :class:`SquareTable` object from the contingency
table.

.. ipython:: python

    sqtab = sm.stats.SquareTable(tab)
    row, col = sqtab.marginal_probabilities
    print(row)
    print(col)


The ``summary`` method prints results for the symmetry and homogeneity
testing procedures.

.. ipython:: python

    print(sqtab.summary())

If we had the individual case records in a dataframe called ``data``,
we could also perform the same analysis by passing the raw data using
the ``SquareTable.from_data`` class method.

::

    sqtab = sm.stats.SquareTable.from_data(data[['left', 'right']])
    print(sqtab.summary())


A single 2x2 table
------------------

Several methods for working with individual 2x2 tables are provided in
the :class:`sm.stats.Table2x2` class.  The ``summary`` method displays
several measures of association between the rows and columns of the
table.

.. ipython:: python

    table = np.asarray([[35, 21], [25, 58]])
    t22 = sm.stats.Table2x2(table)
    print(t22.summary())

Note that the risk ratio is not symmetric so different results will be
obtained if the transposed table is analyzed.

.. ipython:: python

    table = np.asarray([[35, 21], [25, 58]])
    t22 = sm.stats.Table2x2(table.T)
    print(t22.summary())


Stratified 2x2 tables
---------------------

Stratification occurs when we have a collection of contingency tables
defined by the same row and column factors.  In the example below, we
have a collection of 2x2 tables reflecting the joint distribution of
smoking and lung cancer in each of several regions of China.  It is
possible that the tables all have a common odds ratio, even while the
marginal probabilities vary among the strata.  The 'Breslow-Day'
procedure tests whether the data are consistent with a common odds
ratio.  It appears below as the `Test of constant OR`.  The
Mantel-Haenszel procedure tests whether this common odds ratio is
equal to one.  It appears below as the `Test of OR=1`.  It is also
possible to estimate the common odds and risk ratios and obtain
confidence intervals for them.  The ``summary`` method displays all of
these results.  Individual results can be obtained from the class
methods and attributes.

.. ipython:: python

    data = sm.datasets.china_smoking.load_pandas()

    mat = np.asarray(data.data)
    tables = [np.reshape(x.tolist(), (2, 2)) for x in mat]

    st = sm.stats.StratifiedTable(tables)
    print(st.summary())


Module Reference
----------------

.. module:: statsmodels.stats.contingency_tables
   :synopsis: Contingency table analysis

.. currentmodule:: statsmodels.stats.contingency_tables

.. autosummary::
   :toctree: generated/

   Table
   Table2x2
   SquareTable
   StratifiedTable
   mcnemar
   cochrans_q

See also
--------

Scipy_ has several functions for analyzing contingency tables,
including Fisher's exact test which is not currently in statsmodels.

.. _Scipy: https://docs.scipy.org/doc/scipy-0.18.0/reference/stats.html#contingency-table-functions
.. currentmodule:: statsmodels.stats.anova

.. _anova:

ANOVA
=====

Analysis of Variance models containing anova_lm for ANOVA analysis with a
linear OLSModel, and AnovaRM for repeated measures ANOVA, within ANOVA for
balanced data.

Examples
--------

.. ipython:: python

    import statsmodels.api as sm
    from statsmodels.formula.api import ols

    moore = sm.datasets.get_rdataset("Moore", "carData",
                                     cache=True) # load data
    data = moore.data
    data = data.rename(columns={"partner.status":
                                "partner_status"}) # make name pythonic
    moore_lm = ols('conformity ~ C(fcategory, Sum)*C(partner_status, Sum)',
                    data=data).fit()

    table = sm.stats.anova_lm(moore_lm, typ=2) # Type 2 ANOVA DataFrame
    print(table)

A more detailed example for `anova_lm` can be found here:

*  `ANOVA <examples/notebooks/generated/interactions_anova.html>`__

Module Reference
----------------

.. module:: statsmodels.stats.anova
   :synopsis: Analysis of Variance

.. autosummary::
   :toctree: generated/

   anova_lm
   AnovaRM
.. currentmodule:: statsmodels.sandbox.regression.gmm


.. _gmm:


Generalized Method of Moments :mod:`gmm`
========================================

:mod:`statsmodels.gmm` contains model classes and functions that are based on
estimation with Generalized Method of Moments.
Currently the general non-linear case is implemented. An example class for the standard
linear instrumental variable model is included. This has been introduced as a test case, it
works correctly but it does not take the linear structure into account. For the linear
case we intend to introduce a specific implementation which will be faster and numerically
more accurate.

Currently, GMM takes arbitrary non-linear moment conditions and calculates the estimates
either for a given weighting matrix or iteratively by alternating between estimating
the optimal weighting matrix and estimating the parameters. Implementing models with
different moment conditions is done by subclassing GMM. In the minimal implementation
only the moment conditions, `momcond` have to be defined.

.. currentmodule:: statsmodels.sandbox.regression.gmm


Module Reference
""""""""""""""""

.. module:: statsmodels.sandbox.regression.gmm
   :synopsis: A framework for implementing Generalized Method of Moments (GMM)

.. autosummary::
   :toctree: generated/

   GMM
   GMMResults
   IV2SLS
   IVGMM
   IVGMMResults
   IVRegressionResults
   LinearIVGMM
   NonlinearIVGMM
.. currentmodule:: statsmodels.regression.linear_model


.. _regression:

Linear Regression
=================

Linear models with independently and identically distributed errors, and for
errors with heteroscedasticity or autocorrelation. This module allows
estimation by ordinary least squares (OLS), weighted least squares (WLS),
generalized least squares (GLS), and feasible generalized least squares with
autocorrelated AR(p) errors.

See `Module Reference`_ for commands and arguments.

Examples
--------

.. ipython:: python

    # Load modules and data
    import numpy as np
    import statsmodels.api as sm
    spector_data = sm.datasets.spector.load()
    spector_data.exog = sm.add_constant(spector_data.exog, prepend=False)

    # Fit and summarize OLS model
    mod = sm.OLS(spector_data.endog, spector_data.exog)
    res = mod.fit()
    print(res.summary())

Detailed examples can be found here:


* `OLS <examples/notebooks/generated/ols.html>`__
* `WLS <examples/notebooks/generated/wls.html>`__
* `GLS <examples/notebooks/generated/gls.html>`__
* `Recursive LS <examples/notebooks/generated/recursive_ls.html>`__
* `Rolling LS <examples/notebooks/generated/rolling_ls.html>`__

Technical Documentation
-----------------------

The statistical model is assumed to be

 :math:`Y = X\beta + \mu`,  where :math:`\mu\sim N\left(0,\Sigma\right).`

Depending on the properties of :math:`\Sigma`, we have currently four classes available:

* GLS : generalized least squares for arbitrary covariance :math:`\Sigma`
* OLS : ordinary least squares for i.i.d. errors :math:`\Sigma=\textbf{I}`
* WLS : weighted least squares for heteroskedastic errors :math:`\text{diag}\left  (\Sigma\right)`
* GLSAR : feasible generalized least squares with autocorrelated AR(p) errors
  :math:`\Sigma=\Sigma\left(\rho\right)`

All regression models define the same methods and follow the same structure,
and can be used in a similar fashion. Some of them contain additional model
specific methods and attributes.

GLS is the superclass of the other regression classes except for RecursiveLS,
RollingWLS and RollingOLS.

.. Class hierachy: TODO

.. yule_walker is not a full model class, but a function that estimate the
.. parameters of a univariate autoregressive process, AR(p). It is used in GLSAR,
.. but it can also be used independently of any models. yule_walker only
.. calculates the estimates and the standard deviation of the lag parameters but
.. not the additional regression statistics. We hope to include yule-walker in
.. future in a separate univariate time series class. A similar result can be
.. obtained with GLSAR if only the constant is included as regressors. In this
.. case the parameter estimates of the lag estimates are not reported, however
.. additional statistics, for example aic, become available.

References
^^^^^^^^^^

General reference for regression models:

* D.C. Montgomery and E.A. Peck. "Introduction to Linear Regression Analysis." 2nd. Ed., Wiley, 1992.

Econometrics references for regression models:

* R.Davidson and J.G. MacKinnon. "Econometric Theory and Methods," Oxford, 2004.
* W.Green. "Econometric Analysis," 5th ed., Pearson, 2003.

.. toctree::
..   :maxdepth: 1
..
..   regression_techn1

Attributes
^^^^^^^^^^

The following is more verbose description of the attributes which is mostly
common to all regression classes

pinv_wexog : array
    The `p` x `n` Moore-Penrose pseudoinverse of the whitened design matrix.
    It is approximately equal to
    :math:`\left(X^{T}\Sigma^{-1}X\right)^{-1}X^{T}\Psi`, where
    :math:`\Psi` is defined such that :math:`\Psi\Psi^{T}=\Sigma^{-1}`.
cholsimgainv : array
    The `n` x `n` upper triangular matrix :math:`\Psi^{T}` that satisfies
    :math:`\Psi\Psi^{T}=\Sigma^{-1}`.
df_model : float
    The model degrees of freedom. This is equal to `p` - 1, where `p` is the
    number of regressors. Note that the intercept is not counted as using a
    degree of freedom here.
df_resid : float
    The residual degrees of freedom. This is equal `n - p` where `n` is the
    number of observations and `p` is the number of parameters. Note that the
    intercept is counted as using a degree of freedom here.
llf : float
    The value of the likelihood function of the fitted model.
nobs : float
    The number of observations `n`
normalized_cov_params : array
    A `p` x `p` array equal to :math:`(X^{T}\Sigma^{-1}X)^{-1}`.
sigma : array
    The `n` x `n` covariance matrix of the error terms:
    :math:`\mu\sim N\left(0,\Sigma\right)`.
wexog : array
    The whitened design matrix :math:`\Psi^{T}X`.
wendog : array
    The whitened response variable :math:`\Psi^{T}Y`.

Module Reference
----------------

.. module:: statsmodels.regression.linear_model
   :synopsis: Least squares linear models

Model Classes
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   OLS
   GLS
   WLS
   GLSAR
   yule_walker
   burg

.. module:: statsmodels.regression.quantile_regression
   :synopsis: Quantile regression

.. currentmodule:: statsmodels.regression.quantile_regression

.. autosummary::
   :toctree: generated/

   QuantReg

.. module:: statsmodels.regression.recursive_ls
   :synopsis: Recursive least squares using the Kalman Filter

.. currentmodule:: statsmodels.regression.recursive_ls

.. autosummary::
   :toctree: generated/

   RecursiveLS

.. module:: statsmodels.regression.rolling
   :synopsis: Rolling (moving) least squares

.. currentmodule:: statsmodels.regression.rolling

.. autosummary::
   :toctree: generated/

   RollingWLS
   RollingOLS

.. module:: statsmodels.regression.process_regression
   :synopsis: Process regression

.. currentmodule:: statsmodels.regression.process_regression

.. autosummary::
   :toctree: generated/

   GaussianCovariance
   ProcessMLE

.. module:: statsmodels.regression.dimred
   :synopsis: Dimension reduction methods

.. currentmodule:: statsmodels.regression.dimred

.. autosummary::
   :toctree: generated/

    SlicedInverseReg
    PrincipalHessianDirections
    SlicedAverageVarianceEstimation


Results Classes
^^^^^^^^^^^^^^^

Fitting a linear regression model returns a results class. OLS has a
specific results class with some additional methods compared to the
results class of the other linear models.

.. currentmodule:: statsmodels.regression.linear_model

.. autosummary::
   :toctree: generated/

   RegressionResults
   OLSResults
   PredictionResults

.. currentmodule:: statsmodels.base.elastic_net

.. autosummary::
   :toctree: generated/

    RegularizedResults

.. currentmodule:: statsmodels.regression.quantile_regression

.. autosummary::
   :toctree: generated/

   QuantRegResults

.. currentmodule:: statsmodels.regression.recursive_ls

.. autosummary::
   :toctree: generated/

   RecursiveLSResults

.. currentmodule:: statsmodels.regression.rolling

.. autosummary::
   :toctree: generated/

   RollingRegressionResults

.. currentmodule:: statsmodels.regression.process_regression

.. autosummary::
   :toctree: generated/

   ProcessMLEResults

.. currentmodule:: statsmodels.regression.dimred

.. autosummary::
   :toctree: generated/

   DimReductionResults
.. currentmodule:: statsmodels.iolib

.. _iolib:

Input-Output :mod:`iolib`
=========================

``statsmodels`` offers some functions for input and output. These include a
reader for STATA files, a class for generating tables for printing in several
formats and two helper functions for pickling.

Users can also leverage the powerful input/output functions provided by :ref:`pandas.io <pandas:io>`. Among other things, ``pandas`` (a ``statsmodels`` dependency) allows reading and writing to Excel, CSV, and HDF5 (PyTables).

Examples
--------

    `SimpleTable: Basic example <examples/notebooks/generated/wls.html#ols-vs-wls>`__

Module Reference
----------------

.. module:: statsmodels.iolib
   :synopsis: Tools for reading datasets and producing summary output

.. autosummary::
   :toctree: generated/

   foreign.savetxt
   table.SimpleTable
   table.csv2st
   smpickle.save_pickle
   smpickle.load_pickle


The following are classes and functions used to return the summary of
estimation results, and mostly intended for internal use. There are currently
two versions for creating summaries.

.. autosummary::
   :toctree: generated/

   summary.Summary
   summary2.Summary
.. module:: statsmodels.base.distributed_estimation
.. currentmodule:: statsmodels.base.distributed_estimation

Working with Large Data Sets
============================

Big data is something of a buzzword in the modern world. While statsmodels
works well with small and moderately-sized data sets that can be loaded in
memory--perhaps tens of thousands of observations--use cases exist with
millions of observations or more. Depending your use case, statsmodels may or
may not be a sufficient tool.

statsmodels and most of the software stack it is written on operates in
memory. Resultantly, building models on larger data sets can be challenging
or even impractical. With that said, there are 2 general strategies for
building models on larger data sets with statsmodels.

Divide and Conquer - Distributing Jobs
--------------------------------------

If your system is capable of loading all the data, but the analysis you are
attempting to perform is slow, you might be able to build models on horizontal
slices of the data and then aggregate the individual models once fit.

A current limitation of this approach is that it generally does not support
`patsy <https://patsy.readthedocs.io/en/latest/>`_ so constructing your
design matrix (known as `exog`) in statsmodels, is a little challenging.

A detailed example is available
`here <examples/notebooks/generated/distributed_estimation.html>`_.

.. autosummary::
   :toctree: generated/

   DistributedModel
   DistributedResults

Subsetting your data
--------------------

If your entire data set is too large to store in memory, you might try storing
it in a columnar container like `Apache Parquet <https://parquet.apache.org/>`_
or `bcolz <http://bcolz.blosc.org/en/latest/>`_. Using the patsy formula
interface, statsmodels will use the `__getitem__` function (i.e. data['Item'])
to pull only the specified columns.

.. code-block:: python

    import pyarrow as pa
    import pyarrow.parquet as pq
    import statsmodels.formula.api as smf

    class DataSet(dict):
        def __init__(self, path):
            self.parquet = pq.ParquetFile(path)

        def __getitem__(self, key):
            try:
                return self.parquet.read([key]).to_pandas()[key]
            except:
                raise KeyError

    LargeData = DataSet('LargeData.parquet')

    res = smf.ols('Profit ~ Sugar + Power + Women', data=LargeData).fit()

Additionally, you can add code to this example `DataSet` object to return only
a subset of the rows until you have built a good model. Then, you can refit
your final model on more data.
.. currentmodule:: statsmodels.genmod.bayes_mixed_glm

Generalized Linear Mixed Effects Models
=======================================

Generalized Linear Mixed Effects (GLIMMIX) models are generalized
linear models with random effects in the linear predictors.
statsmodels currently supports estimation of binomial and Poisson
GLIMMIX models using two Bayesian methods: the Laplace approximation
to the posterior, and a variational Bayes approximation to the
posterior.  Both methods provide point estimates (posterior means) and
assessments of uncertainty (posterior standard deviation).

The current implementation only supports independent random effects.

Technical Documentation
-----------------------

Unlike statsmodels mixed linear models, the GLIMMIX implementation is
not group-based.  Groups are created by interacting all random effects
with a categorical variable.  Note that this creates large, sparse
random effects design matrices `exog_vc`.  Internally, `exog_vc` is
converted to a scipy sparse matrix.  When passing the arguments
directly to the class initializer, a sparse matrix may be passed.
When using formulas, a dense matrix is created then converted to
sparse.  For very large problems, it may not be feasible to use
formulas due to the size of this dense intermediate matrix.

References
^^^^^^^^^^

Blei, Kucukelbir, McAuliffe (2017).  Variational Inference: A review
for Statisticians https://arxiv.org/pdf/1601.00670.pdf

Module Reference
----------------

.. module:: statsmodels.genmod.bayes_mixed_glm
   :synopsis: Bayes Mixed Generalized Linear Models


The model classes are:

.. autosummary::
   :toctree: generated/

   BinomialBayesMixedGLM
   PoissonBayesMixedGLM

The result class is:

.. autosummary::
   :toctree: generated/

   BayesMixedGLMResults
.. currentmodule:: statsmodels.gam.api

.. _gam:

Generalized Additive Models (GAM)
=================================

Generalized Additive Models allow for penalized estimation of smooth terms
in generalized linear models.

See `Module Reference`_ for commands and arguments.

Examples
--------

The following illustrates a Gaussian and a Poisson regression where
categorical variables are treated as linear terms and the effect of
two explanatory variables is captured by penalized B-splines.
The data is from the automobile dataset
https://archive.ics.uci.edu/ml/datasets/automobile
We can load a dataframe with selected columns from the unit test module.

.. ipython:: python

    import statsmodels.api as sm
    from statsmodels.gam.api import GLMGam, BSplines

    # import data
    from statsmodels.gam.tests.test_penalized import df_autos

    # create spline basis for weight and hp
    x_spline = df_autos[['weight', 'hp']]
    bs = BSplines(x_spline, df=[12, 10], degree=[3, 3])

    # penalization weight
    alpha = np.array([21833888.8, 6460.38479])

    gam_bs = GLMGam.from_formula('city_mpg ~ fuel + drive', data=df_autos,
                                 smoother=bs, alpha=alpha)
    res_bs = gam_bs.fit()
    print(res_bs.summary())

    # plot smooth components
    res_bs.plot_partial(0, cpr=True)
    res_bs.plot_partial(1, cpr=True)

    alpha = np.array([8283989284.5829611, 14628207.58927821])
    gam_bs = GLMGam.from_formula('city_mpg ~ fuel + drive', data=df_autos,
                                 smoother=bs, alpha=alpha,
                                 family=sm.families.Poisson())
    res_bs = gam_bs.fit()
    print(res_bs.summary())

    # Optimal penalization weights alpha can be obtaine through generalized
    # cross-validation or k-fold cross-validation.
    # The alpha above are from the unit tests against the R mgcv package.

    gam_bs.select_penweight()[0]
    gam_bs.select_penweight_kfold()[0]


References
^^^^^^^^^^

* Hastie, Trevor, and Robert Tibshirani. 1986. Generalized Additive Models. Statistical Science 1 (3): 297-310.
* Wood, Simon N. 2006. Generalized Additive Models: An Introduction with R. Texts in Statistical Science. Boca Raton, FL: Chapman & Hall/CRC.
* Wood, Simon N. 2017. Generalized Additive Models: An Introduction with R. Second edition. Chapman & Hall/CRC Texts in Statistical Science. Boca Raton: CRC Press/Taylor & Francis Group.


Module Reference
----------------

.. module:: statsmodels.gam.generalized_additive_model
   :synopsis: Generalized Additive Models
.. currentmodule:: statsmodels.gam.generalized_additive_model

Model Class
^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   GLMGam
   LogitGam

Results Classes
^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   GLMGamResults

Smooth Basis Functions
^^^^^^^^^^^^^^^^^^^^^^

.. module:: statsmodels.gam.smooth_basis
   :synopsis: Classes for Spline and other Smooth Basis Function

.. currentmodule:: statsmodels.gam.smooth_basis

Currently there is verified support for two spline bases

.. autosummary::
   :toctree: generated/

   BSplines
   CyclicCubicSplines

`statsmodels.gam.smooth_basis` includes additional splines and a (global)
polynomial smoother basis but those have not been verified yet.



Families and Link Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distribution families in `GLMGam` are the same as for GLM and so are
the corresponding link functions.
Current unit tests only cover Gaussian and Poisson, and GLMGam might not
work for all options that are available in GLM.
:orphan:

.. currentmodule:: statsmodels

.. _faq:

Frequently Asked Question
=========================

What is statsmodels?
--------------------

statsmodels is a Python package that provides a collection of widely-used
statistical models. While statsmodels historically has an econometrics-heavy
user base, the package is designed to be useful to a large variety of
statistical use cases. In comparison with other Python-based modelling
tools, statsmodels focuses more heavily on the statistics and diagnostics
underlying the models than having the most cutting-edge or predictive models.

.. _endog-exog-faq:

What do endog and exog mean?
----------------------------

These are shorthand for endogenous and exogenous variables. You might be more
comfortable with the common ``y`` and ``X`` notation in linear models.
Sometimes the endogenous variable ``y`` is called a dependent variable.
Likewise, sometimes the exogenous variables ``X`` are called the independent
variables. You can read about this in greater detail at :ref:`endog_exog`

.. _missing-faq:

How does statsmodels handle missing data?
-----------------------------------------

Missing data can be handled via the ``missing`` keyword argument. Every model
takes this keyword. You can find more information in the docstring of
:class:`statsmodels.base.Model <statsmodels.base.model.Model>`.

.. _build-faq:

Why will not statsmodels build?
-------------------------------

Remember that to build, you must have:

- The appropriate dependencies (numpy, pandas, scipy, Cython) installed
- A suitable C compiler
- A working python installation

Please review our :ref:`installation instructions <install>` for details.

You might also try cleaning up your source directory by running:

.. code-block:: bash

    pip uninstall statsmodels
    python setup.py clean

And then attempting to re-compile. If you want to be more aggressive, you
could also reset git to a prior version by:

.. code-block:: bash

    git reset --hard
    git clean -xdf
    git checkout main

I'd like to contribute. Where do I start?
-----------------------------------------

Check out our :doc:`development pages <dev/index>` for a guide on how to
get involved. We accept Pull Requests on our GitHub page for bugfixes and
topics germane to statistics and statistical modeling. In addition, usability
and quality of life enhancements are greatly appreciated as well.

What if my question is not answered here?
-----------------------------------------

You may find answers for questions that have not yet been added here on GitHub
under the `FAQ issues tag <https://github.com/statsmodels/statsmodels/labels/FAQ>`_.
If not, please ask your question on stackoverflow using the
`statsmodels tag <https://stackoverflow.com/questions/tagged/statsmodels>`_ or
on the `mailing list <https://groups.google.com/forum/#!forum/pystatsmodels>`_.
.. image:: images/statsmodels-logo-v2-horizontal.svg
   :width: 50%
   :alt: statsmodels
   :align: left

:ref:`statsmodels <about:About statsmodels>` is a Python module that provides classes and functions for the estimation
of many different statistical models, as well as for conducting statistical tests, and statistical
data exploration. An extensive list of result statistics are available for each estimator.
The results are tested against existing statistical packages to ensure that they are correct. The
package is released under the open source Modified BSD (3-clause) license.
The online documentation is hosted at `statsmodels.org <https://www.statsmodels.org/>`__.

Introduction
============

``statsmodels`` supports specifying models using R-style formulas and ``pandas`` DataFrames.
Here is a simple example using ordinary least squares:

.. ipython:: python

    import numpy as np
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    # Load data
    dat = sm.datasets.get_rdataset("Guerry", "HistData").data

    # Fit regression model (using the natural log of one of the regressors)
    results = smf.ols('Lottery ~ Literacy + np.log(Pop1831)', data=dat).fit()

    # Inspect the results
    print(results.summary())

You can also use ``numpy`` arrays instead of formulas:

.. ipython:: python

    import numpy as np
    import statsmodels.api as sm

    # Generate artificial data (2 regressors + constant)
    nobs = 100
    X = np.random.random((nobs, 2))
    X = sm.add_constant(X)
    beta = [1, .1, .5]
    e = np.random.random(nobs)
    y = np.dot(X, beta) + e

    # Fit regression model
    results = sm.OLS(y, X).fit()

    # Inspect the results
    print(results.summary())

Have a look at `dir(results)` to see available results. Attributes are described in
`results.__doc__` and results methods have their own docstrings.

Citation
========

Please use following citation to cite statsmodels in scientific publications:


Seabold, Skipper, and Josef Perktold. "`statsmodels: Econometric and statistical modeling with
python. <http://conference.scipy.org/proceedings/scipy2010/pdfs/seabold.pdf>`_" *Proceedings
of the 9th Python in Science Conference.* 2010.

Bibtex entry::

  @inproceedings{seabold2010statsmodels,
    title={statsmodels: Econometric and statistical modeling with python},
    author={Seabold, Skipper and Perktold, Josef},
    booktitle={9th Python in Science Conference},
    year={2010},
  }

.. toctree::
   :maxdepth: 1

   install
   gettingstarted
   user-guide
   examples/index
   api
   about
   dev/index
   release/index


Index
=====

:ref:`genindex`

:ref:`modindex`
.. _importpaths:

Import Paths and Structure
--------------------------

We offer two ways of importing functions and classes from statsmodels:

1. `API import for interactive use`_

   + Allows tab completion

2. `Direct import for programs`_

   + Avoids importing unnecessary modules and commands

API Import for interactive use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For interactive use the recommended import is:

.. code-block:: python

    import statsmodels.api as sm

Importing `statsmodels.api` will load most of the public parts of statsmodels.
This makes most functions and classes conveniently available within one or two
levels, without making the "sm" namespace too crowded.

To see what functions and classes available, you can type the following (or use
the namespace exploration features of IPython, Spyder, IDLE, etc.):

.. code-block:: python

    >>> dir(sm)
    ['GLM', 'GLS', 'GLSAR', 'Logit', 'MNLogit', 'OLS', 'Poisson', 'Probit', 'RLM',
    'WLS', '__builtins__', '__doc__', '__file__', '__name__', '__package__',
    'add_constant', 'categorical', 'datasets', 'distributions', 'families',
    'graphics', 'iolib', 'nonparametric', 'qqplot', 'regression', 'robust',
    'stats', 'test', 'tools', 'tsa', 'version']

    >>> dir(sm.graphics)
    ['__builtins__', '__doc__', '__file__', '__name__', '__package__',
    'abline_plot', 'beanplot', 'fboxplot', 'interaction_plot', 'qqplot',
    'rainbow', 'rainbowplot', 'violinplot']

    >>> dir(sm.tsa)
    ['AR', 'ARMA', 'SVAR', 'VAR', '__builtins__', '__doc__',
    '__file__', '__name__', '__package__', 'acf', 'acovf', 'add_lag',
    'add_trend', 'adfuller', 'ccf', 'ccovf', 'datetools', 'detrend',
    'filters', 'grangercausalitytests', 'interp', 'lagmat', 'lagmat2ds', 'kpss',
    'pacf', 'pacf_ols', 'pacf_yw', 'periodogram', 'q_stat', 'range_unit_root_test',
    'stattools', 'tsatools', 'var']

Notes
^^^^^

The `api` modules may not include all the public functionality of statsmodels. If
you find something that should be added to the api, please file an issue on
github or report it to the mailing list.

The subpackages of statsmodels include `api.py` modules that are mainly
intended to collect the imports needed for those subpackages. The `subpackage/api.py`
files are imported into statsmodels api, for example ::

     from .nonparametric import api as nonparametric

Users do not need to load the `subpackage/api.py` modules directly.

Direct import for programs
~~~~~~~~~~~~~~~~~~~~~~~~~~

``statsmodels`` submodules are arranged by topic (e.g. `discrete` for discrete
choice models, or `tsa` for time series analysis). Our directory tree (stripped
down) looks something like this::

    statsmodels/
        __init__.py
        api.py
        discrete/
            __init__.py
            discrete_model.py
            tests/
                results/
        tsa/
            __init__.py
            api.py
            tsatools.py
            stattools.py
            arima_process.py
            vector_ar/
                __init__.py
                var_model.py
                tests/
                    results/
            tests/
                results/
        stats/
            __init__.py
            api.py
            stattools.py
            tests/
        tools/
            __init__.py
            tools.py
            decorators.py
            tests/

The submodules that can be import heavy contain an empty `__init__.py`, except
for some testing code for running tests for the submodules. The intention is to
change all directories to have an `api.py` and empty `__init__.py` in the next
release.

Import examples
^^^^^^^^^^^^^^^

Functions and classes::

    from statsmodels.regression.linear_model import OLS, WLS
    from statsmodels.tools.tools import rank, add_constant

Modules ::

    from statsmodels.datasets import macrodata
    import statsmodels.stats import diagnostic

Modules with aliases ::

    import statsmodels.regression.linear_model as lm
    import statsmodels.stats.diagnostic as smsdia
    import statsmodels.stats.outliers_influence as oi

We do not have currently a convention for aliases of submodules.
.. module:: statsmodels.stats
   :synopsis: Statistical methods and tests

.. currentmodule:: statsmodels.stats

.. _stats:


Statistics :mod:`stats`
=======================

This section collects various statistical tests and tools.
Some can be used independently of any models, some are intended as extension to the
models and model results.

API Warning: The functions and objects in this category are spread out in
various modules and might still be moved around. We expect that in future the
statistical tests will return class instances with more informative reporting
instead of only the raw numbers.


.. _stattools:


Residual Diagnostics and Specification Tests
--------------------------------------------

.. module:: statsmodels.stats.stattools
   :synopsis: Statistical methods and tests that do not fit into other categories

.. currentmodule:: statsmodels.stats.stattools

.. autosummary::
   :toctree: generated/

   durbin_watson
   jarque_bera
   omni_normtest
   medcouple
   robust_skewness
   robust_kurtosis
   expected_robust_kurtosis

.. module:: statsmodels.stats.diagnostic
   :synopsis: Statistical methods and tests to diagnose model fit problems

.. currentmodule:: statsmodels.stats.diagnostic

.. autosummary::
   :toctree: generated/

   acorr_breusch_godfrey
   acorr_ljungbox
   acorr_lm

   breaks_cusumolsresid
   breaks_hansen
   recursive_olsresiduals

   compare_cox
   compare_encompassing
   compare_j

   het_arch
   het_breuschpagan
   het_goldfeldquandt
   het_white
   spec_white

   linear_harvey_collier
   linear_lm
   linear_rainbow
   linear_reset


Outliers and influence measures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: statsmodels.stats.outliers_influence
   :synopsis: Statistical methods and measures for outliers and influence

.. currentmodule:: statsmodels.stats.outliers_influence

.. autosummary::
   :toctree: generated/

   OLSInfluence
   GLMInfluence
   MLEInfluence
   variance_inflation_factor

See also the notes on :ref:`notes on regression diagnostics <diagnostics>`

Sandwich Robust Covariances
---------------------------

The following functions calculate covariance matrices and standard errors for
the parameter estimates that are robust to heteroscedasticity and
autocorrelation in the errors. Similar to the methods that are available
for the LinearModelResults, these methods are designed for use with OLS.

.. currentmodule:: statsmodels.stats

.. autosummary::
   :toctree: generated/

   sandwich_covariance.cov_hac
   sandwich_covariance.cov_nw_panel
   sandwich_covariance.cov_nw_groupsum
   sandwich_covariance.cov_cluster
   sandwich_covariance.cov_cluster_2groups
   sandwich_covariance.cov_white_simple

The following are standalone versions of the heteroscedasticity robust
standard errors attached to LinearModelResults

.. autosummary::
   :toctree: generated/

   sandwich_covariance.cov_hc0
   sandwich_covariance.cov_hc1
   sandwich_covariance.cov_hc2
   sandwich_covariance.cov_hc3

   sandwich_covariance.se_cov


Goodness of Fit Tests and Measures
----------------------------------

some tests for goodness of fit for univariate distributions

.. module:: statsmodels.stats.gof
   :synopsis: Goodness of fit measures and tests

.. currentmodule:: statsmodels.stats.gof

.. autosummary::
   :toctree: generated/

   powerdiscrepancy
   gof_chisquare_discrete
   gof_binning_discrete
   chisquare_effectsize

.. currentmodule:: statsmodels.stats.diagnostic

.. autosummary::
   :toctree: generated/

   anderson_statistic
   normal_ad
   kstest_exponential
   kstest_fit
   kstest_normal
   lilliefors

Non-Parametric Tests
--------------------

.. module:: statsmodels.sandbox.stats.runs
   :synopsis: Experimental statistical methods and tests to analyze runs

.. currentmodule:: statsmodels.sandbox.stats.runs

.. autosummary::
   :toctree: generated/

   mcnemar
   symmetry_bowker
   median_test_ksample
   runstest_1samp
   runstest_2samp
   cochrans_q
   Runs

.. currentmodule:: statsmodels.stats.descriptivestats

.. autosummary::
   :toctree: generated/

   sign_test

.. currentmodule:: statsmodels.stats.nonparametric

.. autosummary::
   :toctree: generated/

   rank_compare_2indep
   rank_compare_2ordinal
   RankCompareResult
   cohensd2problarger
   prob_larger_continuous
   rankdata_2samp


Descriptive Statistics
----------------------

.. module:: statsmodels.stats.descriptivestats
   :synopsis: Descriptive statistics

.. currentmodule:: statsmodels.stats.descriptivestats

.. autosummary::
   :toctree: generated/

   describe
   Description

.. _interrater:

Interrater Reliability and Agreement
------------------------------------

The main function that statsmodels has currently available for interrater
agreement measures and tests is Cohen's Kappa. Fleiss' Kappa is currently
only implemented as a measures but without associated results statistics.

.. module:: statsmodels.stats.inter_rater
.. currentmodule:: statsmodels.stats.inter_rater

.. autosummary::
   :toctree: generated/

   cohens_kappa
   fleiss_kappa
   to_table
   aggregate_raters

Multiple Tests and Multiple Comparison Procedures
-------------------------------------------------

`multipletests` is a function for p-value correction, which also includes p-value
correction based on fdr in `fdrcorrection`.
`tukeyhsd` performs simultaneous testing for the comparison of (independent) means.
These three functions are verified.
GroupsStats and MultiComparison are convenience classes to multiple comparisons similar
to one way ANOVA, but still in development

.. module:: statsmodels.sandbox.stats.multicomp
   :synopsis: Experimental methods for controlling size while performing multiple comparisons


.. currentmodule:: statsmodels.stats.multitest

.. autosummary::
   :toctree: generated/

   multipletests
   fdrcorrection

.. currentmodule:: statsmodels.sandbox.stats.multicomp

.. autosummary::
   :toctree: generated/

   GroupsStats
   MultiComparison
   TukeyHSDResults

.. module:: statsmodels.stats.multicomp
   :synopsis: Methods for controlling size while performing multiple comparisons

.. currentmodule:: statsmodels.stats.multicomp

.. autosummary::
   :toctree: generated/

   pairwise_tukeyhsd

.. module:: statsmodels.stats.multitest
   :synopsis: Multiple testing p-value and FDR adjustments

.. currentmodule:: statsmodels.stats.multitest

.. autosummary::
   :toctree: generated/

   local_fdr
   fdrcorrection_twostage
   NullDistribution
   RegressionFDR

.. module:: statsmodels.stats.knockoff_regeffects
   :synopsis: Regression Knock-Off Effects

.. currentmodule:: statsmodels.stats.knockoff_regeffects

.. autosummary::
   :toctree: generated/

   CorrelationEffects
   OLSEffects
   ForwardEffects
   OLSEffects
   RegModelEffects

The following functions are not (yet) public

.. currentmodule:: statsmodels.sandbox.stats.multicomp

.. autosummary::
   :toctree: generated/

   varcorrection_pairs_unbalanced
   varcorrection_pairs_unequal
   varcorrection_unbalanced
   varcorrection_unequal

   StepDown
   catstack
   ccols
   compare_ordered
   distance_st_range
   ecdf
   get_tukeyQcrit
   homogeneous_subsets
   maxzero
   maxzerodown
   mcfdr
   qcrit
   randmvn
   rankdata
   rejectionline
   set_partition
   set_remove_subs
   tiecorrect

.. _tost:

Basic Statistics and t-Tests with frequency weights
---------------------------------------------------

Besides basic statistics, like mean, variance, covariance and correlation for
data with case weights, the classes here provide one and two sample tests
for means. The t-tests have more options than those in scipy.stats, but are
more restrictive in the shape of the arrays. Confidence intervals for means
are provided based on the same assumptions as the t-tests.

Additionally, tests for equivalence of means are available for one sample and
for two, either paired or independent, samples. These tests are based on TOST,
two one-sided tests, which have as null hypothesis that the means are not
"close" to each other.

.. module:: statsmodels.stats.weightstats
   :synopsis: Weighted statistics

.. currentmodule:: statsmodels.stats.weightstats

.. autosummary::
   :toctree: generated/

   DescrStatsW
   CompareMeans
   ttest_ind
   ttost_ind
   ttost_paired
   ztest
   ztost
   zconfint

weightstats also contains tests and confidence intervals based on summary
data

.. currentmodule:: statsmodels.stats.weightstats

.. autosummary::
   :toctree: generated/

   _tconfint_generic
   _tstat_generic
   _zconfint_generic
   _zstat_generic
   _zstat_generic2


Power and Sample Size Calculations
----------------------------------

The :mod:`power` module currently implements power and sample size calculations
for the t-tests, normal based test, F-tests and Chisquare goodness of fit test.
The implementation is class based, but the module also provides
three shortcut functions, ``tt_solve_power``, ``tt_ind_solve_power`` and
``zt_ind_solve_power`` to solve for any one of the parameters of the power
equations.


.. module:: statsmodels.stats.power
   :synopsis: Power and size calculations for common tests

.. currentmodule:: statsmodels.stats.power

.. autosummary::
   :toctree: generated/

   TTestIndPower
   TTestPower
   GofChisquarePower
   NormalIndPower
   FTestAnovaPower
   FTestPower
   normal_power_het
   normal_sample_size_one_tail
   tt_solve_power
   tt_ind_solve_power
   zt_ind_solve_power


.. _proportion_stats:

Proportion
----------

Also available are hypothesis test, confidence intervals and effect size for
proportions that can be used with NormalIndPower.

.. module:: statsmodels.stats.proportion
   :synopsis: Tests for proportions

.. currentmodule:: statsmodels.stats.proportion

.. autosummary::
   :toctree: generated

   proportion_confint
   proportion_effectsize

   binom_test
   binom_test_reject_interval
   binom_tost
   binom_tost_reject_interval

   multinomial_proportions_confint

   proportions_ztest
   proportions_ztost
   proportions_chisquare
   proportions_chisquare_allpairs
   proportions_chisquare_pairscontrol

   proportion_effectsize
   power_binom_tost
   power_ztost_prop
   samplesize_confint_proportion

Statistics for two independent samples
Status: experimental, API might change, added in 0.12

.. autosummary::
   :toctree: generated

   test_proportions_2indep
   confint_proportions_2indep
   power_proportions_2indep
   tost_proportions_2indep
   samplesize_proportions_2indep_onetail
   score_test_proportions_2indep
   _score_confint_inversion


Rates
-----

Statistical functions for rates. This currently includes hypothesis tests for
two independent samples.

Status: experimental, API might change, added in 0.12

.. module:: statsmodels.stats.rates
   :synopsis: Tests for Poisson rates

.. currentmodule:: statsmodels.stats.rates

.. autosummary::
   :toctree: generated

   test_poisson_2indep
   etest_poisson_2indep
   tost_poisson_2indep


Multivariate
------------

Statistical functions for multivariate samples.

This includes hypothesis test and confidence intervals for mean of sample
of multivariate observations and hypothesis tests for the structure of a
covariance matrix.

Status: experimental, API might change, added in 0.12

.. module:: statsmodels.stats.multivariate
   :synopsis: Statistical functions for multivariate samples.

.. currentmodule:: statsmodels.stats.multivariate

.. autosummary::
   :toctree: generated

   test_mvmean
   confint_mvmean
   confint_mvmean_fromstats
   test_mvmean_2indep
   test_cov
   test_cov_blockdiagonal
   test_cov_diagonal
   test_cov_oneway
   test_cov_spherical


.. _oneway_stats:

Oneway Anova
------------

Hypothesis test, confidence intervals and effect size for oneway analysis of
k samples.

Status: experimental, API might change, added in 0.12

.. module:: statsmodels.stats.oneway
   :synopsis: Statistical functions for oneway analysis, Anova.

.. currentmodule:: statsmodels.stats.oneway

.. autosummary::
   :toctree: generated


   anova_oneway
   anova_generic
   equivalence_oneway
   equivalence_oneway_generic
   power_equivalence_oneway
   _power_equivalence_oneway_emp

   test_scale_oneway
   equivalence_scale_oneway

   confint_effectsize_oneway
   confint_noncentrality
   convert_effectsize_fsqu
   effectsize_oneway
   f2_to_wellek
   fstat_to_wellek
   wellek_to_f2
   _fstat2effectsize

   scale_transform
   simulate_power_equivalence_oneway


.. _robust_stats:

Robust, Trimmed Statistics
--------------------------

Statistics for samples that are trimmed at a fixed fraction. This includes
class TrimmedMean for one sample statistics. It is used in `stats.oneway`
for trimmed "Yuen" Anova.

Status: experimental, API might change, added in 0.12

.. module:: statsmodels.stats.robust_compare
   :synopsis: Trimmed sample statistics.

.. currentmodule:: statsmodels.stats.robust_compare

.. autosummary::
   :toctree: generated

   TrimmedMean
   scale_transform
   trim_mean
   trimboth


Moment Helpers
--------------

When there are missing values, then it is possible that a correlation or
covariance matrix is not positive semi-definite. The following
functions can be used to find a correlation or covariance matrix that is
positive definite and close to the original matrix.
Additional functions estimate spatial covariance matrix and regularized
inverse covariance or precision matrix.

.. module:: statsmodels.stats.correlation_tools
   :synopsis: Procedures for ensuring correlations are positive semi-definite

.. currentmodule:: statsmodels.stats.correlation_tools

.. autosummary::
   :toctree: generated/

   corr_clipped
   corr_nearest
   corr_nearest_factor
   corr_thresholded
   cov_nearest
   cov_nearest_factor_homog
   FactoredPSDMatrix
   kernel_covariance

.. currentmodule:: statsmodels.stats.regularized_covariance

.. autosummary::
   :toctree: generated/

   RegularizedInvCovariance

These are utility functions to convert between central and non-central moments, skew,
kurtosis and cummulants.

.. module:: statsmodels.stats.moment_helpers
   :synopsis: Tools for converting moments

.. currentmodule:: statsmodels.stats.moment_helpers

.. autosummary::
   :toctree: generated/

   cum2mc
   mc2mnc
   mc2mvsk
   mnc2cum
   mnc2mc
   mnc2mvsk
   mvsk2mc
   mvsk2mnc
   cov2corr
   corr2cov
   se_cov


Mediation Analysis
------------------

Mediation analysis focuses on the relationships among three key variables:
an 'outcome', a 'treatment', and a 'mediator'. Since mediation analysis is a
form of causal inference, there are several assumptions involved that are
difficult or impossible to verify. Ideally, mediation analysis is conducted in
the context of an experiment such as this one in which the treatment is
randomly assigned. It is also common for people to conduct mediation analyses
using observational data in which the treatment may be thought of as an
'exposure'. The assumptions behind mediation analysis are even more difficult
to verify in an observational setting.

.. module:: statsmodels.stats.mediation
   :synopsis: Mediation analysis

.. currentmodule:: statsmodels.stats.mediation

.. autosummary::
   :toctree: generated/

   Mediation
   MediationResults


Oaxaca-Blinder Decomposition
----------------------------

The Oaxaca-Blinder, or Blinder-Oaxaca as some call it, decomposition attempts to explain
gaps in means of groups. It uses the linear models of two given regression equations to
show what is explained by regression coefficients and known data and what is unexplained
using the same data. There are two types of Oaxaca-Blinder decompositions, the two-fold
and the three-fold, both of which can and are used in Economics Literature to discuss
differences in groups. This method helps classify discrimination or unobserved effects.
This function attempts to port the functionality of the oaxaca command in STATA to Python.

.. module:: statsmodels.stats.oaxaca
   :synopsis: Oaxaca-Blinder Decomposition

.. currentmodule:: statsmodels.stats.oaxaca

.. autosummary::
   :toctree: generated/

   OaxacaBlinder
   OaxacaResults


Distance Dependence Measures
----------------------------

Distance dependence measures and the Distance Covariance (dCov) test.

.. module:: statsmodels.stats.dist_dependence_measures
   :synopsis: Distance Dependence Measures

.. currentmodule:: statsmodels.stats.dist_dependence_measures

.. autosummary::
   :toctree: generated/

   distance_covariance_test
   distance_statistics
   distance_correlation
   distance_covariance
   distance_variance


Meta-Analysis
-------------

Functions for basic meta-analysis of a collection of sample statistics.

Examples can be found in the notebook

 * `Meta-Analysis <examples/notebooks/generated/metaanalysis1.html>`__

Status: experimental, API might change, added in 0.12

.. module:: statsmodels.stats.meta_analysis
   :synopsis: Meta-Analysis

.. currentmodule:: statsmodels.stats.meta_analysis

.. autosummary::
   :toctree: generated/

   combine_effects
   effectsize_2proportions
   effectsize_smd
   CombineResults

The module also includes internal functions to compute random effects
variance.


.. autosummary::
   :toctree: generated/

   _fit_tau_iter_mm
   _fit_tau_iterative
   _fit_tau_mm
.. module:: statsmodels
   :synopsis: Statistical analysis in Python

.. currentmodule:: statsmodels

*****************
About statsmodels
*****************

Background
----------

The ``models`` module of ``scipy.stats`` was originally written by Jonathan
Taylor. For some time it was part of scipy but was later removed. During
the Google Summer of Code 2009, ``statsmodels`` was corrected, tested,
improved and released as a new package. Since then, the statsmodels
development team has continued to add new models, plotting tools, and
statistical methods.

Testing
-------

Most results have been verified with at least one other statistical package:
R, Stata or SAS. The guiding principle for the initial rewrite and for
continued development is that all numbers have to be verified. Some
statistical methods are tested with Monte Carlo studies. While we strive to
follow this test-driven approach, there is no guarantee that the code is
bug-free and always works. Some auxiliary function are still insufficiently
tested, some edge cases might not be correctly taken into account, and the
possibility of numerical problems is inherent to many of the statistical
models. We especially appreciate any help and reports for these kind of
problems so we can keep improving the existing models.

Code Stability
^^^^^^^^^^^^^^

The existing models are mostly settled in their user interface and we do not
expect many large changes going forward. For the existing code, although
there is no guarantee yet on API stability, we have long deprecation periods
in all but very special cases, and we try to keep changes that require
adjustments by existing users to a minimal level. For newer models we might
adjust the user interface as we gain more experience and obtain feedback.
These changes will always be noted in our release notes available in the
documentation.

Reporting Bugs
^^^^^^^^^^^^^^
If you encounter a bug or an unexpected behavior, please report it on
`the issue tracker <https://github.com/statsmodels/statsmodels/issues>`_.
Use the ``show_versions`` command to list the installed versions of
statsmodels and its dependencies.

.. autosummary::
   :toctree: generated/

   ~statsmodels.tools.print_version.show_versions


Financial Support
-----------------

We are grateful for the financial support that we obtained for the
development of statsmodels:

* Google `www.google.com <https://www.google.com/>`_ : Google Summer of Code
  (GSOC) 2009-2017.
* AQR `www.aqr.com <https://www.aqr.com/>`_ : financial sponsor for the work on
  Vector Autoregressive Models (VAR) by Wes McKinney

We would also like to thank our hosting providers, `github
<https://github.com/>`_ for the public code repository, `github.io
<https://www.statsmodels.org/stable/index.html>`_ for hosting our documentation
and `python.org <https://www.python.org/>`_ for making our downloads available
on PyPi.

We also thank our current and previous continuous integration providers, `Azure <https://azure.microsoft.com/>`_,
`Travis CI <https://travis-ci.org/>`_ and `AppVeyor <https://ci.appveyor.com>`_ for
unit testing, and `Codecov <https://codecov.io>`_ and `Coveralls <https://coveralls.io>`_ for
code coverage.

Brand Marks
-----------

Please make use of the statsmodels logos when preparing demonstrations involving
statsmodels code.

Color
^^^^^

+----------------+---------------------+
| Horizontal     | |color-horizontal|  |
+----------------+---------------------+
| Vertical       | |color-vertical|    |
+----------------+---------------------+
| Logo Only      | |color-notext|      |
+----------------+---------------------+

Monochrome (Dark)
^^^^^^^^^^^^^^^^^

+----------------+---------------------+
| Horizontal     | |dark-horizontal|   |
+----------------+---------------------+
| Vertical       | |dark-vertical|     |
+----------------+---------------------+
| Logo Only      | |dark-notext|       |
+----------------+---------------------+

Monochrome (Light)
^^^^^^^^^^^^^^^^^^

.. note::

   The light brand marks are light grey on transparent, and so are difficult to see on this
   page. They are intended for use on a dark background.


+----------------+---------------------+
| Horizontal     | |light-horizontal|  |
+----------------+---------------------+
| Vertical       | |light-vertical|    |
+----------------+---------------------+
| Logo Only      | |light-notext|      |
+----------------+---------------------+

.. |color-horizontal| image:: images/statsmodels-logo-v2-horizontal.svg
   :width: 50%

.. |color-vertical| image:: images/statsmodels-logo-v2.svg
   :width: 14%

.. |color-notext| image:: images/statsmodels-logo-v2-no-text.svg
   :width: 9%

.. |dark-horizontal| image:: images/statsmodels-logo-v2-horizontal-dark.svg
   :width: 50%

.. |dark-vertical| image:: images/statsmodels-logo-v2-dark.svg
   :width: 14%

.. |dark-notext| image:: images/statsmodels-logo-v2-no-text-dark.svg
   :width: 9%

.. |light-horizontal| image:: images/statsmodels-logo-v2-horizontal-light.svg
   :width: 50%

.. |light-vertical| image:: images/statsmodels-logo-v2-light.svg
   :width: 14%

.. |light-notext| image:: images/statsmodels-logo-v2-no-text-light.svg
   :width: 9%
:orphan:

Patsy: Contrast Coding Systems for categorical variables
===========================================================

.. note:: This document is based on `this excellent resource from UCLA <https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/>`__.

A categorical variable of K categories, or levels, usually enters a regression as a sequence of K-1 dummy variables. This amounts to a linear hypothesis on the level means. That is, each test statistic for these variables amounts to testing whether the mean for that level is statistically significantly different from the mean of the base category. This dummy coding is called Treatment coding in R parlance, and we will follow this convention. There are, however, different coding methods that amount to different sets of linear hypotheses.

In fact, the dummy coding is not technically a contrast coding. This is because the dummy variables add to one and are not functionally independent of the model's intercept. On the other hand, a set of *contrasts* for a categorical variable with `k` levels is a set of `k-1` functionally independent linear combinations of the factor level means that are also independent of the sum of the dummy variables. The dummy coding is not wrong *per se*. It captures all of the coefficients, but it complicates matters when the model assumes independence of the coefficients such as in ANOVA. Linear regression models do not assume independence of the coefficients and thus dummy coding is often the only coding that is taught in this context.

To have a look at the contrast matrices in Patsy, we will use data from UCLA ATS. First let's load the data.

.. ipython::
   :suppress:

   In [1]: import numpy as np
      ...: np.set_printoptions(precision=4, suppress=True)
      ...:
      ...: from patsy.contrasts import ContrastMatrix
      ...:
      ...: def _name_levels(prefix, levels):
      ...:     return ["[%s%s]" % (prefix, level) for level in levels]

   In [2]: class Simple(object):
      ...:     def _simple_contrast(self, levels):
      ...:         nlevels = len(levels)
      ...:         contr = -1./nlevels * np.ones((nlevels, nlevels-1))
      ...:         contr[1:][np.diag_indices(nlevels-1)] = (nlevels-1.)/nlevels
      ...:         return contr
      ...:
      ...:     def code_with_intercept(self, levels):
      ...:         contrast = np.column_stack((np.ones(len(levels)),
      ...:                                    self._simple_contrast(levels)))
      ...:         return ContrastMatrix(contrast, _name_levels("Simp.", levels))
      ...:
      ...:     def code_without_intercept(self, levels):
      ...:         contrast = self._simple_contrast(levels)
      ...:         return ContrastMatrix(contrast, _name_levels("Simp.", levels[:-1]))
      ...:

Example Data
------------

.. ipython:: python

   import pandas
   url = 'https://stats.idre.ucla.edu/stat/data/hsb2.csv'
   hsb2 = pandas.read_csv(url)

It will be instructive to look at the mean of the dependent variable, write, for each level of race ((1 = Hispanic, 2 = Asian, 3 = African American and 4 = Caucasian)).

.. ipython:: python

   hsb2.groupby('race')['write'].mean()

Treatment (Dummy) Coding
------------------------

Dummy coding is likely the most well known coding scheme. It compares each level of the categorical variable to a base reference level. The base reference level is the value of the intercept. It is the default contrast in Patsy for unordered categorical factors. The Treatment contrast matrix for race would be

.. ipython:: python

   from patsy.contrasts import Treatment
   levels = [1,2,3,4]
   contrast = Treatment(reference=0).code_without_intercept(levels)
   print(contrast.matrix)

Here we used `reference=0`, which implies that the first level, Hispanic, is the reference category against which the other level effects are measured. As mentioned above, the columns do not sum to zero and are thus not independent of the intercept. To be explicit, let's look at how this would encode the `race` variable.

.. ipython:: python

   contrast.matrix[hsb2.race-1, :][:20]

This is a bit of a trick, as the `race` category conveniently maps to zero-based indices. If it does not, this conversion happens under the hood, so this will not work in general but nonetheless is a useful exercise to fix ideas. The below illustrates the output using the three contrasts above

.. ipython:: python

   from statsmodels.formula.api import ols
   mod = ols("write ~ C(race, Treatment)", data=hsb2)
   res = mod.fit()
   print(res.summary())

We explicitly gave the contrast for race; however, since Treatment is the default, we could have omitted this.

Simple Coding
-------------

Like Treatment Coding, Simple Coding compares each level to a fixed reference level. However, with simple coding, the intercept is the grand mean of all the levels of the factors. See :ref:`user-defined` for how to implement the Simple contrast.


.. ipython:: python

   contrast = Simple().code_without_intercept(levels)
   print(contrast.matrix)

   mod = ols("write ~ C(race, Simple)", data=hsb2)
   res = mod.fit()
   print(res.summary())

Sum (Deviation) Coding
----------------------

Sum coding compares the mean of the dependent variable for a given level to the overall mean of the dependent variable over all the levels. That is, it uses contrasts between each of the first k-1 levels and level k In this example, level 1 is compared to all the others, level 2 to all the others, and level 3 to all the others.

.. ipython:: python

   from patsy.contrasts import Sum
   contrast = Sum().code_without_intercept(levels)
   print(contrast.matrix)

   mod = ols("write ~ C(race, Sum)", data=hsb2)
   res = mod.fit()
   print(res.summary())

This corresponds to a parameterization that forces all the coefficients to sum to zero. Notice that the intercept here is the grand mean where the grand mean is the mean of means of the dependent variable by each level.

.. ipython:: python

   hsb2.groupby('race')['write'].mean().mean()

Backward Difference Coding
--------------------------

In backward difference coding, the mean of the dependent variable for a level is compared with the mean of the dependent variable for the prior level. This type of coding may be useful for a nominal or an ordinal variable.

.. ipython:: python

   from patsy.contrasts import Diff
   contrast = Diff().code_without_intercept(levels)
   print(contrast.matrix)

   mod = ols("write ~ C(race, Diff)", data=hsb2)
   res = mod.fit()
   print(res.summary())

For example, here the coefficient on level 1 is the mean of `write` at level 2 compared with the mean at level 1. Ie.,

.. ipython:: python

   res.params["C(race, Diff)[D.1]"]
   hsb2.groupby('race').mean()["write"][2] - \
       hsb2.groupby('race').mean()["write"][1]

Helmert Coding
--------------

Our version of Helmert coding is sometimes referred to as Reverse Helmert Coding. The mean of the dependent variable for a level is compared to the mean of the dependent variable over all previous levels. Hence, the name 'reverse' being sometimes applied to differentiate from forward Helmert coding. This comparison does not make much sense for a nominal variable such as race, but we would use the Helmert contrast like so:

.. ipython:: python

   from patsy.contrasts import Helmert
   contrast = Helmert().code_without_intercept(levels)
   print(contrast.matrix)

   mod = ols("write ~ C(race, Helmert)", data=hsb2)
   res = mod.fit()
   print(res.summary())

To illustrate, the comparison on level 4 is the mean of the dependent variable at the previous three levels taken from the mean at level 4

.. ipython:: python

   grouped = hsb2.groupby('race')
   grouped.mean()["write"][4] - grouped.mean()["write"][:3].mean()

As you can see, these are only equal up to a constant. Other versions of the Helmert contrast give the actual difference in means. Regardless, the hypothesis tests are the same.

.. ipython:: python

   k = 4
   1./k * (grouped.mean()["write"][k] - grouped.mean()["write"][:k-1].mean())
   k = 3
   1./k * (grouped.mean()["write"][k] - grouped.mean()["write"][:k-1].mean())


Orthogonal Polynomial Coding
----------------------------

The coefficients taken on by polynomial coding for `k=4` levels are the linear, quadratic, and cubic trends in the categorical variable. The categorical variable here is assumed to be represented by an underlying, equally spaced numeric variable. Therefore, this type of encoding is used only for ordered categorical variables with equal spacing. In general, the polynomial contrast produces polynomials of order `k-1`. Since `race` is not an ordered factor variable let's use `read` as an example. First we need to create an ordered categorical from `read`.

.. ipython:: python

   _, bins = np.histogram(hsb2.read, 3)
   try: # requires numpy main
       readcat = np.digitize(hsb2.read, bins, True)
   except:
       readcat = np.digitize(hsb2.read, bins)
   hsb2['readcat'] = readcat
   hsb2.groupby('readcat').mean()['write']

.. ipython:: python

   from patsy.contrasts import Poly
   levels = hsb2.readcat.unique().tolist()
   contrast = Poly().code_without_intercept(levels)
   print(contrast.matrix)

   mod = ols("write ~ C(readcat, Poly)", data=hsb2)
   res = mod.fit()
   print(res.summary())

As you can see, readcat has a significant linear effect on the dependent variable `write` but not a significant quadratic or cubic effect.

.. _user-defined:

User-Defined Coding
-------------------

If you want to use your own coding, you must do so by writing a coding class that contains a code_with_intercept and a code_without_intercept method that return a `patsy.contrast.ContrastMatrix` instance.

.. ipython::

   In [1]: from patsy.contrasts import ContrastMatrix
      ...:
      ...: def _name_levels(prefix, levels):
      ...:     return ["[%s%s]" % (prefix, level) for level in levels]

   In [2]: class Simple(object):
      ...:     def _simple_contrast(self, levels):
      ...:         nlevels = len(levels)
      ...:         contr = -1./nlevels * np.ones((nlevels, nlevels-1))
      ...:         contr[1:][np.diag_indices(nlevels-1)] = (nlevels-1.)/nlevels
      ...:         return contr
      ...:
      ...:     def code_with_intercept(self, levels):
      ...:         contrast = np.column_stack((np.ones(len(levels)),
      ...:                                    self._simple_contrast(levels)))
      ...:         return ContrastMatrix(contrast, _name_levels("Simp.", levels))
      ...:
      ...:    def code_without_intercept(self, levels):
      ...:        contrast = self._simple_contrast(levels)
      ...:        return ContrastMatrix(contrast, _name_levels("Simp.", levels[:-1]))

   In [3]: mod = ols("write ~ C(race, Simple)", data=hsb2)
      ...: res = mod.fit()
      ...: print(res.summary())
:orphan:

.. _diagnostics:

Regression Diagnostics and Specification Tests
==============================================


Introduction
------------

In many cases of statistical analysis, we are not sure whether our statistical
model is correctly specified. For example when using ols, then linearity and
homoscedasticity are assumed, some test statistics additionally assume that
the errors are normally distributed or that we have a large sample.
Since our results depend on these statistical assumptions, the results are
only correct of our assumptions hold (at least approximately).

One solution to the problem of uncertainty about the correct specification is
to use robust methods, for example robust regression or robust covariance
(sandwich) estimators. The second approach is to test whether our sample is
consistent with these assumptions.

The following briefly summarizes specification and diagnostics tests for
linear regression.

Heteroscedasticity Tests
------------------------

For these test the null hypothesis is that all observations have the same
error variance, i.e. errors are homoscedastic. The tests differ in which kind
of heteroscedasticity is considered as alternative hypothesis. They also vary
in the power of the test for different types of heteroscedasticity.

:py:func:`het_breuschpagan <statsmodels.stats.diagnostic.het_breuschpagan>`
    Lagrange Multiplier Heteroscedasticity Test by Breusch-Pagan

:py:func:`het_white <statsmodels.stats.diagnostic.het_white>`
    Lagrange Multiplier Heteroscedasticity Test by White

:py:func:`het_goldfeldquandt <statsmodels.stats.diagnostic.het_goldfeldquandt>`
    test whether variance is the same in 2 subsamples


Autocorrelation Tests
---------------------

This group of test whether the regression residuals are not autocorrelated.
They assume that observations are ordered by time.

:py:func:`durbin_watson <statsmodels.stats.diagnostic.durbin_watson>`
  - Durbin-Watson test for no autocorrelation of residuals
  - printed with summary()

:py:func:`acorr_ljungbox <statsmodels.stats.diagnostic.acorr_ljungbox>`
  - Ljung-Box test for no autocorrelation of residuals
  - also returns Box-Pierce statistic

:py:func:`acorr_breusch_godfrey <statsmodels.stats.diagnostic.acorr_breusch_godfrey>`
  - Breusch-Pagan test for no autocorrelation of residuals


missing
  - ?


Non-Linearity Tests
-------------------

:py:func:`linear_harvey_collier <statsmodels.stats.diagnostic.linear_harvey_collier>`
  - Multiplier test for Null hypothesis that linear specification is
    correct

:py:func:`acorr_linear_rainbow <statsmodels.stats.diagnostic.acorr_linear_rainbow>`
  - Multiplier test for Null hypothesis that linear specification is
    correct.

:py:func:`acorr_linear_lm <statsmodels.stats.diagnostic.acorr_linear_lm>`
  - Lagrange Multiplier test for Null hypothesis that linear specification is
    correct. This tests against specific functional alternatives.

:py:func:`spec_white <statsmodels.stats.diagnostic.spec_white>`
  - White's two-moment specification test with null hypothesis of homoscedastic
    and correctly specified.

Tests for Structural Change, Parameter Stability
------------------------------------------------

Test whether all or some regression coefficient are constant over the
entire data sample.

Known Change Point
^^^^^^^^^^^^^^^^^^

OneWayLS :
  - flexible ols wrapper for testing identical regression coefficients across
    predefined subsamples (eg. groups)

missing
  - predictive test: Greene, number of observations in subsample is smaller than
    number of regressors


Unknown Change Point
^^^^^^^^^^^^^^^^^^^^

:py:func:`breaks_cusumolsresid <statsmodels.stats.diagnostic.breaks_cusumolsresid>`
  - cusum test for parameter stability based on ols residuals

:py:func:`breaks_hansen <statsmodels.stats.diagnostic.breaks_hansen>`
  - test for model stability, breaks in parameters for ols, Hansen 1992

:py:func:`recursive_olsresiduals <statsmodels.stats.diagnostic.recursive_olsresiduals>`
  Calculate recursive ols with residuals and cusum test statistic. This is
  currently mainly helper function for recursive residual based tests.
  However, since it uses recursive updating and does not estimate separate
  problems it should be also quite efficient as expanding OLS function.

missing
  - supLM, expLM, aveLM  (Andrews, Andrews/Ploberger)
  - R-structchange also has musum (moving cumulative sum tests)
  - test on recursive parameter estimates, which are there?


Multicollinearity Tests
--------------------------------

conditionnum (statsmodels.stattools)
  - -- needs test vs Stata --
  - cf Grene (3rd ed.) pp 57-8

numpy.linalg.cond
  - (for more general condition numbers, but no behind the scenes help for
    design preparation)

Variance Inflation Factors
  This is currently together with influence and outlier measures
  (with some links to other tests here: http://www.stata.com/help.cgi?vif)


Normality and Distribution Tests
--------------------------------

:py:func:`jarque_bera <statsmodels.stats.tools.jarque_bera>`
  - printed with summary()
  - test for normal distribution of residuals

Normality tests in scipy stats
  need to find list again

:py:func:`omni_normtest <statsmodels.stats.tools.omni_normtest>`
  - test for normal distribution of residuals
  - printed with summary()

:py:func:`normal_ad <statsmodels.stats.diagnostic.normal_ad>`
  - Anderson Darling test for normality with estimated mean and variance

:py:func:`kstest_normal <statsmodels.stats.diagnostic.kstest_normal>` :py:func:`lilliefors <statsmodels.stats.diagnostic.lilliefors>`
  Lilliefors test for normality, this is a Kolmogorov-Smirnov tes with for
  normality with estimated mean and variance. lilliefors is an alias for
  kstest_normal

qqplot, scipy.stats.probplot

other goodness-of-fit tests for distributions in scipy.stats and enhancements
  - kolmogorov-smirnov
  - anderson : Anderson-Darling
  - likelihood-ratio, ...
  - chisquare tests, powerdiscrepancy : needs wrapping (for binning)


Outlier and Influence Diagnostic Measures
-----------------------------------------

These measures try to identify observations that are outliers, with large
residual, or observations that have a large influence on the regression
estimates. Robust Regression, RLM, can be used to both estimate in an outlier
robust way as well as identify outlier. The advantage of RLM that the
estimation results are not strongly influenced even if there are many
outliers, while most of the other measures are better in identifying
individual outliers and might not be able to identify groups of outliers.

:py:class:`RLM <statsmodels.robust.robust_linear_model.RLM>`
    example from example_rlm.py ::

        import statsmodels.api as sm

        ### Example for using Huber's T norm with the default
        ### median absolute deviation scaling

        data = sm.datasets.stackloss.load()
        data.exog = sm.add_constant(data.exog)
        huber_t = sm.RLM(data.endog, data.exog, M=sm.robust.norms.HuberT())
        hub_results = huber_t.fit()
        print(hub_results.weights)

    And the weights give an idea of how much a particular observation is
    down-weighted according to the scaling asked for.

:py:class:`Influence <statsmodels.stats.outliers_influence.OLSInfluence>`
   Class in stats.outliers_influence, most standard measures for outliers
   and influence are available as methods or attributes given a fitted
   OLS model. This is mainly written for OLS, some but not all measures
   are also valid for other models.
   Some of these statistics can be calculated from an OLS results instance,
   others require that an OLS is estimated for each left out variable.

   - resid_press
   - resid_studentized_external
   - resid_studentized_internal
   - ess_press
   - hat_matrix_diag
   - cooks_distance - Cook's Distance `Wikipedia <https://en.wikipedia.org/wiki/Cook%27s_distance>`_ (with some other links)
   - cov_ratio
   - dfbetas
   - dffits
   - dffits_internal
   - det_cov_params_not_obsi
   - params_not_obsi
   - sigma2_not_obsi
:orphan:

.. _install:

Installing statsmodels
======================

The easiest way to install statsmodels is to install it as part of the `Anaconda <https://docs.continuum.io/anaconda/>`_
distribution, a cross-platform distribution for data analysis and scientific
computing. This is the recommended installation method for most users.

Instructions for installing from PyPI, source or a development version are also provided.


Python Support
--------------

statsmodels supports Python 3.8, 3.9, and 3.10.

Anaconda
--------
statsmodels is available through conda provided by
`Anaconda <https://www.anaconda.com/products/individual#Downloads>`__. The latest release can
be installed using:

.. code-block:: bash

   conda install -c conda-forge statsmodels

PyPI (pip)
----------

To obtain the latest released version of statsmodels using pip:

.. code-block:: bash

    pip install statsmodels

Follow `this link to our PyPI page <https://pypi.org/project/statsmodels/>`__ to directly
download wheels or source.

For Windows users, unofficial recent binaries (wheels) are occasionally
available `here <https://www.lfd.uci.edu/~gohlke/pythonlibs/#statsmodels>`__.

Obtaining the Source
--------------------

We do not release very often but the main branch of our source code is
usually fine for everyday use. You can get the latest source from our
`github repository <https://github.com/statsmodels/statsmodels>`__. Or if you
have git installed:

.. code-block:: bash

    git clone git://github.com/statsmodels/statsmodels.git

If you want to keep up to date with the source on github just periodically do:

.. code-block:: bash

    git pull

in the statsmodels directory.

Installation from Source
------------------------

You will need a C compiler installed to build statsmodels. If you are building
from the github source and not a source release, then you will also need
Cython. You can follow the instructions below to get a C compiler setup for
Windows.

If your system is already set up with pip, a compiler, and git, you can try:

.. code-block:: bash

    pip install git+https://github.com/statsmodels/statsmodels

If you do not have pip installed or want to do the installation more manually,
you can also type:

.. code-block:: bash

    python setup.py install

Or even more manually

.. code-block:: bash

    python setup.py build
    python setup.py install

statsmodels can also be installed in `develop` mode which installs statsmodels
into the current python environment in-place. The advantage of this is that
edited modules will immediately be re-interpreted when the python interpreter
restarts without having to re-install statsmodels.

.. code-block:: bash

    python setup.py develop

Compilers
~~~~~~~~~

Linux
^^^^^

If you are using Linux, we assume that you are savvy enough to install `gcc` on
your own. More than likely, it is already installed.

Windows
^^^^^^^

It is strongly recommended to use 64-bit Python if possible.

Getting the right compiler is especially confusing for Windows users. Over time,
Python has been built using a variety of different Windows C compilers.
`This guide <https://wiki.python.org/moin/WindowsCompilers>`_ should help
clarify which version of Python uses which compiler by default.

Mac
^^^

Installing statsmodels on MacOS requires installing `gcc` which provides
a suitable C compiler. We recommend installing Xcode and the Command Line
Tools.

Dependencies
------------

The current minimum dependencies are:

* `Python <https://www.python.org>`__ >= 3.8
* `NumPy <https://www.scipy.org/>`__ >= 1.18
* `SciPy <https://www.scipy.org/>`__ >= 1.4
* `Pandas <https://pandas.pydata.org/>`__ >= 1.0
* `Patsy <https://patsy.readthedocs.io/en/latest/>`__ >= 0.5.2

Cython is required to build from a git checkout but not to run or install from PyPI:

* `Cython <https://cython.org/>`__ >= 0.29.26 is required to build the code from
  github but not from a source distribution.

Given the long release cycle, statsmodels follows a loose time-based policy for
dependencies: minimal dependencies are lagged about one and a half to two
years. Our next planned update of minimum versions is expected in the first
half of 2020.

Optional Dependencies
---------------------

* `cvxopt <https://cvxopt.org/>`__ is required for regularized fitting of
  some models.
* `Matplotlib <https://matplotlib.org/>`__ >= 3 is needed for plotting
  functions and running many of the examples.
* If installed, `X-12-ARIMA <https://www.census.gov/srd/www/x13as/>`__ or
  `X-13ARIMA-SEATS <https://www.census.gov/srd/www/x13as/>`__ can be used
  for time-series analysis.
* `pytest <https://docs.pytest.org/en/latest/>`__ is required to run
  the test suite.
* `IPython <https://ipython.org>`__ >= 6.0 is required to build the
  docs locally or to use the notebooks.
* `joblib <https://joblib.readthedocs.io/>`__ >= 1.0can be used to accelerate distributed
  estimation for certain models.
* `jupyter <https://jupyter.org/>`__ is needed to run the notebooks.
.. module:: statsmodels.base.optimizer
.. currentmodule:: statsmodels.base.optimizer

Optimization
============

statsmodels uses three types of algorithms for the estimation of the parameters
of a model.

  1. Basic linear models such as :ref:`WLS and OLS <regression>` are directly
     estimated using appropriate linear algebra.
  2. :ref:`RLM <rlm>` and :ref:`GLM <glm>`, use iteratively re-weighted
     least squares. However, you can optionally select one of the scipy
     optimizers discussed below.
  3. For all other models, we use
     `optimizers <https://docs.scipy.org/doc/scipy/reference/optimize.html>`_
     from `scipy <https://docs.scipy.org/doc/scipy/reference/index.html>`_.

Where practical, certain models allow for the optional selection of a
scipy optimizer. A particular scipy optimizer might be default or an option.
Depending on the model and the data, choosing an appropriate scipy optimizer
enables avoidance of a local minima, fitting models in less time, or fitting a
model with less memory.

statsmodels supports the following optimizers along with keyword arguments
associated with that specific optimizer:

- ``newton`` - Newton-Raphson iteration. While not directly from scipy, we
  consider it an optimizer because only the score and hessian are required.

    tol : float
        Relative error in params acceptable for convergence.

- ``nm`` - scipy's ``fmin_nm``

    xtol : float
        Relative error in params acceptable for convergence
    ftol : float
        Relative error in loglike(params) acceptable for
        convergence
    maxfun : int
        Maximum number of function evaluations to make.

- ``bfgs`` - Broyden–Fletcher–Goldfarb–Shanno optimization, scipy's
  ``fmin_bfgs``.

      gtol : float
          Stop when norm of gradient is less than gtol.
      norm : float
          Order of norm (np.Inf is max, -np.Inf is min)
      epsilon
          If fprime is approximated, use this value for the step
          size. Only relevant if LikelihoodModel.score is None.

- ``lbfgs`` - A more memory-efficient (limited memory) implementation of
  ``bfgs``. Scipy's ``fmin_l_bfgs_b``.

      m : int
          The maximum number of variable metric corrections used to
          define the limited memory matrix. (The limited memory BFGS
          method does not store the full hessian but uses this many
          terms in an approximation to it.)
      pgtol : float
          The iteration will stop when
          ``max{|proj g_i | i = 1, ..., n} <= pgtol`` where pg_i is
          the i-th component of the projected gradient.
      factr : float
          The iteration stops when
          ``(f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr * eps``,
          where eps is the machine precision, which is automatically
          generated by the code. Typical values for factr are: 1e12
          for low accuracy; 1e7 for moderate accuracy; 10.0 for
          extremely high accuracy. See Notes for relationship to
          ftol, which is exposed (instead of factr) by the
          scipy.optimize.minimize interface to L-BFGS-B.
      maxfun : int
          Maximum number of iterations.
      epsilon : float
          Step size used when approx_grad is True, for numerically
          calculating the gradient
      approx_grad : bool
          Whether to approximate the gradient numerically (in which
          case func returns only the function value).

- ``cg`` - Conjugate gradient optimization. Scipy's ``fmin_cg``.

      gtol : float
          Stop when norm of gradient is less than gtol.
      norm : float
          Order of norm (np.Inf is max, -np.Inf is min)
      epsilon : float
          If fprime is approximated, use this value for the step
          size. Can be scalar or vector.  Only relevant if
          Likelihoodmodel.score is None.

- ``ncg`` - Newton conjugate gradient. Scipy's ``fmin_ncg``.

      fhess_p : callable f'(x, \*args)
          Function which computes the Hessian of f times an arbitrary
          vector, p.  Should only be supplied if
          LikelihoodModel.hessian is None.
      avextol : float
          Stop when the average relative error in the minimizer
          falls below this amount.
      epsilon : float or ndarray
          If fhess is approximated, use this value for the step size.
          Only relevant if Likelihoodmodel.hessian is None.

- ``powell`` - Powell's method. Scipy's ``fmin_powell``.

      xtol : float
          Line-search error tolerance
      ftol : float
          Relative error in loglike(params) for acceptable for
          convergence.
      maxfun : int
          Maximum number of function evaluations to make.
      start_direc : ndarray
          Initial direction set.

- ``basinhopping`` - Basin hopping. This is part of scipy's ``basinhopping``
  tools.

      niter : integer
          The number of basin hopping iterations.
      niter_success : integer
          Stop the run if the global minimum candidate remains the
          same for this number of iterations.
      T : float
          The "temperature" parameter for the accept or reject
          criterion. Higher "temperatures" mean that larger jumps
          in function value will be accepted. For best results
          `T` should be comparable to the separation (in function
          value) between local minima.
      stepsize : float
          Initial step size for use in the random displacement.
      interval : integer
          The interval for how often to update the `stepsize`.
      minimizer : dict
          Extra keyword arguments to be passed to the minimizer
          `scipy.optimize.minimize()`, for example 'method' - the
          minimization method (e.g. 'L-BFGS-B'), or 'tol' - the
          tolerance for termination. Other arguments are mapped from
          explicit argument of `fit`:
          - `args` <- `fargs`
          - `jac` <- `score`
          - `hess` <- `hess`

- ``minimize`` - Allows the use of any scipy optimizer.

  min_method : str, optional
      Name of minimization method to use.
      Any method specific arguments can be passed directly.
      For a list of methods and their arguments, see
      documentation of `scipy.optimize.minimize`.
      If no method is specified, then BFGS is used.

Model Class
-----------

Generally, there is no need for an end-user to directly call these functions
and classes. However, we provide the class because the different optimization
techniques have unique keyword arguments that may be useful to the user.

.. autosummary::
   :toctree: generated/

   Optimizer
   _fit_newton
   _fit_bfgs
   _fit_lbfgs
   _fit_nm
   _fit_cg
   _fit_ncg
   _fit_powell
   _fit_basinhopping
Getting started
===============

This very simple case-study is designed to get you up-and-running quickly with
``statsmodels``. Starting from raw data, we will show the steps needed to
estimate a statistical model and to draw a diagnostic plot. We will only use
functions provided by ``statsmodels`` or its ``pandas`` and ``patsy``
dependencies.

Loading modules and functions
-----------------------------

After `installing statsmodels and its dependencies <install.html>`_, we load a
few modules and functions:

.. ipython:: python

    import statsmodels.api as sm
    import pandas
    from patsy import dmatrices

`pandas <https://pandas.pydata.org/>`_ builds on ``numpy`` arrays to provide
rich data structures and data analysis tools. The ``pandas.DataFrame`` function
provides labelled arrays of (potentially heterogenous) data, similar to the
``R`` "data.frame". The ``pandas.read_csv`` function can be used to convert a
comma-separated values file to a ``DataFrame`` object.

`patsy <https://github.com/pydata/patsy>`_ is a Python library for describing
statistical models and building `Design Matrices
<https://en.wikipedia.org/wiki/Design_matrix>`_ using ``R``-like formulas.

.. note::

   This example uses the API interface.  See :ref:`importpaths` for information on
   the difference between importing the API interfaces (``statsmodels.api`` and
   ``statsmodels.tsa.api``) and directly importing from the module that defines
   the model.

Data
----

We download the `Guerry dataset
<https://vincentarelbundock.github.io/Rdatasets/doc/HistData/Guerry.html>`_, a
collection of historical data used in support of Andre-Michel Guerry's 1833
*Essay on the Moral Statistics of France*. The data set is hosted online in
comma-separated values format (CSV) by the `Rdatasets
<https://github.com/vincentarelbundock/Rdatasets/>`_ repository.
We could download the file locally and then load it using ``read_csv``, but
``pandas`` takes care of all of this automatically for us:

.. ipython:: python

    df = sm.datasets.get_rdataset("Guerry", "HistData").data

The `Input/Output doc page <iolib.html>`_ shows how to import from various
other formats.

We select the variables of interest and look at the bottom 5 rows:

.. ipython:: python

    vars = ['Department', 'Lottery', 'Literacy', 'Wealth', 'Region']
    df = df[vars]
    df[-5:]

Notice that there is one missing observation in the *Region* column. We
eliminate it using a ``DataFrame`` method provided by ``pandas``:

.. ipython:: python

    df = df.dropna()
    df[-5:]

Substantive motivation and model
--------------------------------

We want to know whether literacy rates in the 86 French departments are
associated with per capita wagers on the Royal Lottery in the 1820s. We need to
control for the level of wealth in each department, and we also want to include
a series of dummy variables on the right-hand side of our regression equation to
control for unobserved heterogeneity due to regional effects. The model is
estimated using ordinary least squares regression (OLS).


Design matrices (*endog* & *exog*)
----------------------------------

To fit most of the models covered by ``statsmodels``, you will need to create
two design matrices. The first is a matrix of endogenous variable(s) (i.e.
dependent, response, regressand, etc.). The second is a matrix of exogenous
variable(s) (i.e. independent, predictor, regressor, etc.). The OLS coefficient
estimates are calculated as usual:

.. math::

    \hat{\beta} = (X'X)^{-1} X'y

where :math:`y` is an :math:`N \times 1` column of data on lottery wagers per
capita (*Lottery*). :math:`X` is :math:`N \times 7` with an intercept, the
*Literacy* and *Wealth* variables, and 4 region binary variables.

The ``patsy`` module provides a convenient function to prepare design matrices
using ``R``-like formulas. You can find more information `here <https://patsy.readthedocs.io/en/latest/>`_.

We use ``patsy``'s ``dmatrices`` function to create design matrices:

.. ipython:: python

    y, X = dmatrices('Lottery ~ Literacy + Wealth + Region', data=df, return_type='dataframe')

The resulting matrices/data frames look like this:

.. ipython:: python

    y[:3]
    X[:3]

Notice that ``dmatrices`` has

* split the categorical *Region* variable into a set of indicator variables.
* added a constant to the exogenous regressors matrix.
* returned ``pandas`` DataFrames instead of simple numpy arrays. This is useful because DataFrames allow ``statsmodels`` to carry-over meta-data (e.g. variable names) when reporting results.

The above behavior can of course be altered. See the `patsy doc pages
<https://patsy.readthedocs.io/en/latest/>`_.

Model fit and summary
---------------------

Fitting a model in ``statsmodels`` typically involves 3 easy steps:

1. Use the model class to describe the model
2. Fit the model using a class method
3. Inspect the results using a summary method

For OLS, this is achieved by:

.. ipython:: python

    mod = sm.OLS(y, X)    # Describe model
    res = mod.fit()       # Fit model
    print(res.summary())   # Summarize model


The ``res`` object has many useful attributes. For example, we can extract
parameter estimates and r-squared by typing:


.. ipython:: python

    res.params
    res.rsquared

Type ``dir(res)`` for a full list of attributes.

For more information and examples, see the `Regression doc page <regression.html>`_

Diagnostics and specification tests
-----------------------------------

``statsmodels`` allows you to conduct a range of useful `regression diagnostics
and specification tests
<stats.html#residual-diagnostics-and-specification-tests>`_.  For instance,
apply the Rainbow test for linearity (the null hypothesis is that the
relationship is properly modelled as linear):

.. ipython:: python

    sm.stats.linear_rainbow(res)

Admittedly, the output produced above is not very verbose, but we know from
reading the `docstring <generated/statsmodels.stats.diagnostic.linear_rainbow.html>`_
(also, ``print(sm.stats.linear_rainbow.__doc__)``) that the
first number is an F-statistic and that the second is the p-value.

``statsmodels`` also provides graphics functions. For example, we can draw a
plot of partial regression for a set of regressors by:

.. ipython:: python

    @savefig gettingstarted_0.png
    sm.graphics.plot_partregress('Lottery', 'Wealth', ['Region', 'Literacy'],
                                 data=df, obs_labels=False)

Documentation
-------------
Documentation can be accessed from an IPython session
using :func:`~statsmodels.tools.web.webdoc`.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   ~statsmodels.tools.web.webdoc

More
----

Congratulations! You're ready to move on to other topics in the
`Table of Contents <index.html#table-of-contents>`_
.. module:: statsmodels.imputation.mice
   :synopsis: Multiple imputation for missing data

.. currentmodule:: statsmodels.imputation.mice

.. _imputation:


Multiple Imputation with Chained Equations
==========================================

The MICE module allows most statsmodels models to be fit to a dataset
with missing values on the independent and/or dependent variables, and
provides rigorous standard errors for the fitted parameters.  The
basic idea is to treat each variable with missing values as the
dependent variable in a regression, with some or all of the remaining
variables as its predictors.  The MICE procedure cycles through these
models, fitting each in turn, then uses a procedure called "predictive
mean matching" (PMM) to generate random draws from the predictive
distributions determined by the fitted models.  These random draws
become the imputed values for one imputed data set.

By default, each variable with missing variables is modeled using a
linear regression with main effects for all other variables in the
data set.  Note that even when the imputation model is linear, the PMM
procedure preserves the domain of each variable.  Thus, for example,
if all observed values for a given variable are positive, all imputed
values for the variable will always be positive.  The user also has
the option to specify which model is used to produce imputed values
for each variable.

.. code


Classes
-------

.. currentmodule:: statsmodels.imputation.mice

.. autosummary::
   :toctree: generated/

   MICE
   MICEData

.. currentmodule:: statsmodels.imputation.bayes_mi

.. autosummary::
   :toctree: generated/

   MI
   BayesGaussMI


Implementation Details
----------------------

Internally, this function uses
`pandas.isnull <https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html#working-with-missing-data>`_.
Anything that returns True from this function will be treated as missing data.
API Reference
=============
The main statsmodels API is split into models:

* ``statsmodels.api``: Cross-sectional models and methods. Canonically imported
  using ``import statsmodels.api as sm``.
* ``statsmodels.tsa.api``: Time-series models and methods. Canonically imported
  using ``import statsmodels.tsa.api as tsa``.
* ``statsmodels.formula.api``: A convenience interface for specifying models
  using formula strings and DataFrames. This API directly exposes the ``from_formula``
  class method of models that support the formula API. Canonically imported using
  ``import statsmodels.formula.api as smf``

The API focuses on models and the most frequently used statistical test, and tools.
:ref:`api-structure:Import Paths and Structure` explains the design of the two API modules and how
importing from the API differs from directly importing from the module where the
model is defined. See the detailed topic pages in the :ref:`user-guide:User Guide` for a complete
list of available models, statistics, and tools.

``statsmodels.api``
-------------------

Regression
~~~~~~~~~~
.. autosummary::

   ~statsmodels.regression.linear_model.OLS
   ~statsmodels.regression.linear_model.WLS
   ~statsmodels.regression.linear_model.GLS
   ~statsmodels.regression.linear_model.GLSAR
   ~statsmodels.regression.recursive_ls.RecursiveLS
   ~statsmodels.regression.rolling.RollingOLS
   ~statsmodels.regression.rolling.RollingWLS

Imputation
~~~~~~~~~~
.. autosummary::

   ~statsmodels.imputation.bayes_mi.BayesGaussMI
   ~statsmodels.imputation.bayes_mi.MI
   ~statsmodels.imputation.mice.MICE
   ~statsmodels.imputation.mice.MICEData

Generalized Estimating Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::

   ~statsmodels.genmod.generalized_estimating_equations.GEE
   ~statsmodels.genmod.generalized_estimating_equations.NominalGEE
   ~statsmodels.genmod.generalized_estimating_equations.OrdinalGEE

Generalized Linear Models
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::

   ~statsmodels.genmod.generalized_linear_model.GLM
   ~statsmodels.gam.generalized_additive_model.GLMGam
   ~statsmodels.genmod.bayes_mixed_glm.BinomialBayesMixedGLM
   ~statsmodels.genmod.bayes_mixed_glm.PoissonBayesMixedGLM

Discrete and Count Models
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::

   ~statsmodels.discrete.discrete_model.Logit
   ~statsmodels.discrete.discrete_model.Probit
   ~statsmodels.discrete.discrete_model.MNLogit
   ~statsmodels.miscmodels.ordinal_model.OrderedModel
   ~statsmodels.discrete.discrete_model.Poisson
   ~statsmodels.discrete.discrete_model.NegativeBinomial
   ~statsmodels.discrete.discrete_model.NegativeBinomialP
   ~statsmodels.discrete.discrete_model.GeneralizedPoisson
   ~statsmodels.discrete.count_model.ZeroInflatedPoisson
   ~statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialP
   ~statsmodels.discrete.count_model.ZeroInflatedGeneralizedPoisson
   ~statsmodels.discrete.conditional_models.ConditionalLogit
   ~statsmodels.discrete.conditional_models.ConditionalMNLogit
   ~statsmodels.discrete.conditional_models.ConditionalPoisson


Multivariate Models
~~~~~~~~~~~~~~~~~~~
.. autosummary::

   ~statsmodels.multivariate.factor.Factor
   ~statsmodels.multivariate.manova.MANOVA
   ~statsmodels.multivariate.pca.PCA

Other Models
~~~~~~~~~~~~
.. autosummary::

   ~statsmodels.regression.mixed_linear_model.MixedLM
   ~statsmodels.duration.survfunc.SurvfuncRight
   ~statsmodels.duration.hazard_regression.PHReg
   ~statsmodels.regression.quantile_regression.QuantReg
   ~statsmodels.robust.robust_linear_model.RLM
   ~statsmodels.othermod.betareg.BetaModel

Graphics
~~~~~~~~
.. autosummary::

   ~statsmodels.graphics.gofplots.ProbPlot
   ~statsmodels.graphics.gofplots.qqline
   ~statsmodels.graphics.gofplots.qqplot
   ~statsmodels.graphics.gofplots.qqplot_2samples

Statistics
~~~~~~~~~~
.. autosummary::

   ~statsmodels.stats.descriptivestats.Description
   ~statsmodels.stats.descriptivestats.describe

Tools
~~~~~
.. autosummary::

   ~statsmodels.__init__.test
   ~statsmodels.tools.tools.add_constant
   ~statsmodels.iolib.smpickle.load_pickle
   ~statsmodels.tools.print_version.show_versions
   ~statsmodels.tools.web.webdoc


``statsmodels.tsa.api``
-----------------------

Statistics and Tests
~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.stattools.acf
   ~statsmodels.tsa.stattools.acovf
   ~statsmodels.tsa.stattools.adfuller
   ~statsmodels.tsa.stattools.bds
   ~statsmodels.tsa.stattools.ccf
   ~statsmodels.tsa.stattools.ccovf
   ~statsmodels.tsa.stattools.coint
   ~statsmodels.tsa.stattools.kpss
   ~statsmodels.tsa.stattools.pacf
   ~statsmodels.tsa.stattools.pacf_ols
   ~statsmodels.tsa.stattools.pacf_yw
   ~statsmodels.tsa.stattools.q_stat
   ~statsmodels.tsa.stattools.range_unit_root_test
   ~statsmodels.tsa.stattools.zivot_andrews


Univariate Time-Series Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.ar_model.AutoReg
   ~statsmodels.tsa.ardl.ARDL
   ~statsmodels.tsa.arima.model.ARIMA
   ~statsmodels.tsa.statespace.sarimax.SARIMAX
   ~statsmodels.tsa.ardl.ardl_select_order
   ~statsmodels.tsa.stattools.arma_order_select_ic
   ~statsmodels.tsa.arima_process.arma_generate_sample
   ~statsmodels.tsa.arima_process.ArmaProcess
   ~statsmodels.tsa.ardl.UECM

Exponential Smoothing
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.holtwinters.ExponentialSmoothing
   ~statsmodels.tsa.holtwinters.Holt
   ~statsmodels.tsa.holtwinters.SimpleExpSmoothing
   ~statsmodels.tsa.statespace.exponential_smoothing.ExponentialSmoothing
   ~statsmodels.tsa.exponential_smoothing.ets.ETSModel


Multivariate Time Series Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.statespace.dynamic_factor.DynamicFactor
   ~statsmodels.tsa.statespace.dynamic_factor_mq.DynamicFactorMQ
   ~statsmodels.tsa.vector_ar.var_model.VAR
   ~statsmodels.tsa.statespace.varmax.VARMAX
   ~statsmodels.tsa.vector_ar.svar_model.SVAR
   ~statsmodels.tsa.vector_ar.vecm.VECM
   ~statsmodels.tsa.statespace.structural.UnobservedComponents

Filters and Decompositions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.seasonal.seasonal_decompose
   ~statsmodels.tsa.seasonal.STL
   ~statsmodels.tsa.filters.bk_filter.bkfilter
   ~statsmodels.tsa.filters.cf_filter.cffilter
   ~statsmodels.tsa.filters.hp_filter.hpfilter

Markov Regime Switching Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.regime_switching.markov_autoregression.MarkovAutoregression
   ~statsmodels.tsa.regime_switching.markov_regression.MarkovRegression

Forecasting
~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.forecasting.stl.STLForecast
   ~statsmodels.tsa.forecasting.theta.ThetaModel

Time-Series Tools
~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.tsatools.add_lag
   ~statsmodels.tsa.tsatools.add_trend
   ~statsmodels.tsa.tsatools.detrend
   ~statsmodels.tsa.tsatools.lagmat
   ~statsmodels.tsa.tsatools.lagmat2ds
   ~statsmodels.tsa.deterministic.DeterministicProcess

X12/X13 Interface
~~~~~~~~~~~~~~~~~

.. autosummary::

   ~statsmodels.tsa.x13.x13_arima_analysis
   ~statsmodels.tsa.x13.x13_arima_select_order

``statsmodels.formula.api``
---------------------------

Models
~~~~~~

The lower case names are aliases to the `from_formula` method of the
corresponding model class. The function descriptions of the methods exposed in
the formula API are generic. See the documentation for the parent model for
details.

.. autosummary::
   :toctree: generated/

   ~statsmodels.formula.api.gls
   ~statsmodels.formula.api.wls
   ~statsmodels.formula.api.ols
   ~statsmodels.formula.api.glsar
   ~statsmodels.formula.api.mixedlm
   ~statsmodels.formula.api.glm
   ~statsmodels.formula.api.gee
   ~statsmodels.formula.api.ordinal_gee
   ~statsmodels.formula.api.nominal_gee
   ~statsmodels.formula.api.rlm
   ~statsmodels.formula.api.logit
   ~statsmodels.formula.api.probit
   ~statsmodels.formula.api.mnlogit
   ~statsmodels.formula.api.poisson
   ~statsmodels.formula.api.negativebinomial
   ~statsmodels.formula.api.quantreg
   ~statsmodels.formula.api.phreg
   ~statsmodels.formula.api.glmgam
   ~statsmodels.formula.api.conditional_logit
   ~statsmodels.formula.api.conditional_mnlogit
   ~statsmodels.formula.api.conditional_poisson
{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block docstring %}
   {% endblock %}


:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}

:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}


:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :exclude-members: {% for item in methods %}{%- if not item.startswith('_') or item in ['__call__'] %}{{ item }},{% endif %}{%- endfor %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:

   {% for item in methods %}
   {%- if not item.startswith('_') or item in ['__call__'] %}   ~{{ name }}.{{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Properties

   .. autosummary::
      :toctree:

   {% for item in attributes %}
   {%- if not item.startswith('_') or item in ['__call__'] %}   ~{{ name }}.{{ item }}
   {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}
:orphan:

==============
Release 0.10.0
==============

Release summary
===============

statsmodels is using github to store the updated documentation. Two version are available:

* `Stable <https://www.statsmodels.org/stable/>`_, the latest release
* `Development <https://www.statsmodels.org/devel/>`_, the latest build of the main branch

**Warning**

API stability is not guaranteed for new features, although even in
this case changes will be made in a backwards compatible way if
possible. The stability of a new feature depends on how much time it
was already in statsmodels main and how much usage it has already
seen.  If there are specific known problems or limitations, then they
are mentioned in the docstrings.

Stats
-----
**Issues Closed**: 1052
**Pull Requests Merged**: 469

The list of pull requests for this release can be found on `github
<https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.10/>`_
(The list does not include some pull request that were merged before
the 0.9 release but not included in 0.9.)


The Highlights
==============

Generalized Additive Models
---------------------------

:class:`~statsmodels.gam.generalized_additive_model.GLMGam` adds support for Generalized additive models.

.. note::

    **Status: experimental**. This class has full unit test coverage for the core
    results with Gaussian and Poisson (without offset and exposure). Other
    options and additional results might not be correctly supported yet.
    (Binomial with counts, i.e. with n_trials, is most likely wrong in parts.
    User specified var or freq weights are most likely also not correct for
    all results.)

:class:`~statsmodels.gam.generalized_additive_model.LogitGam` adds a Logit version, although this is
unfinished. 

Conditional Models
------------------
Three conditional limited dependent variables models have been added:
:class:`~statsmodels.discrete.conditional_models.ConditionalLogit`,
:class:`~statsmodels.discrete.conditional_models.ConditionalMNLogit` and 
:class:`~statsmodels.discrete.conditional_models.ConditionalPoisson`. These are known
as fixed effect models in Econometrics. 

Dimension Reduction Methods
---------------------------
Three standard methods to perform dimension reduction when modeling data have been added:
:class:`~statsmodels.regression.dimred.SlicedInverseReg`,
:class:`~statsmodels.regression.dimred.PrincipalHessianDirections`, and
:class:`~statsmodels.regression.dimred.SlicedAverageVarianceEstimation`.

Regression using Quadratic Inference Functions (QIF)
----------------------------------------------------
Quadratic Inference Function, :class:`~statsmodels.genmod.qif.QIF`, improve the estimation of GEE models.

Gaussian Process Regression
---------------------------
:class:`~statsmodels.regression.process_regression.GaussianCovariance` implements Gaussian process
regression which is a nonparametric kernel-based method to model data.
:class:`~statsmodels.regression.process_regression.ProcessMLE` is a generic class that can be used
for other types of process regression. The results are returned in a
:class:`~statsmodels.regression.process_regression.ProcessMLEResults`.
:func:`~statsmodels.stats.correlation_tools.kernel_covariance`
provides a method that uses kernel averaging to estimate a multivariate covariance function.

Burg's Method
-------------
Burg's method, :func:`~statsmodels.regression.linear_model.burg`, provides an alternative estimator for the parameters
of AR models that is known to work well in small samples. It minimizes the forward and backward errors.

Time series Tools
-----------------
A number of common helper function for decomposing a time series have been added:
:func:`~statsmodels.tsa.innovations.arma_innovations.arma_innovations`, 
:func:`~statsmodels.tsa.stattools.innovations_algo`, and
:func:`~statsmodels.tsa.stattools.innovations_filter`. Two new PACF estimators have been added:
:func:`~statsmodels.tsa.stattools.levinson_durbin_pacf` and :func:`~statsmodels.tsa.stattools.pacf_burg`.

Other
-----
Knockoff effect estimation has been added for a many models:
:class:`~statsmodels.stats.knockoff_regeffects.RegModelEffects`,
:class:`~statsmodels.stats.knockoff_regeffects.CorrelationEffects`,
:class:`~statsmodels.stats.knockoff_regeffects.OLSEffects`,
:class:`~statsmodels.stats.knockoff_regeffects.ForwardEffects`, and 
:class:`~statsmodels.stats.knockoff_regeffects.OLSEffects`.

Influence functions are available for GLM and generic MLE models:
:class:`~statsmodels.stats.outliers_influence.GLMInfluence` and 
:class:`~statsmodels.stats.outliers_influence.MLEInfluence`.


What's new - an overview
========================

The following lists the main new features of statsmodels 0.10. In addition,
release 0.10 includes bug fixes, refactorings and improvements in many areas.

Submodules
----------

``base``
~~~~~~~~
- Add ``ModelWarning`` base class to avoid warning filter on standard UserWarning (:pr:`4712`)
- Add ultra-high screening with SCAD (:pr:`4683`)
- Add penalized mle scad (:pr:`4576`, :issue:`3677`, :issue:`2374`)
- Add score/LM conditional moment tests (:pr:`2096`)
- Fixed a bug which resulted in weights not being used in penalized models (:pr:`5762`, :issue:`4725`)
- Allow the constant index to be located even when ``hasconst=False`` (:pr:`5680`)
- Ensure ``mle_retvals`` is always set even when ``full_output=False`` (:pr:`5681`, :issue:`2752`)
- Fix a bug in Wald tests when testing a single constraint (:pr:`5684`, :issue:`5475`)
- Improve performance by skipping constant check when ``hasconst=True`` (:pr:`5698`)
- Deprecated ``scale`` parameter in the base model class (:pr:`5614`, :issue:`4598`)
- Fixed a bug that raised an error when a multi-index DataFrame was input into a model (:pr:`5634`, :issue:`5415`, :issue:`5414`)
- Fix bug in use of ``self.score`` in GenericLikelihoodModel (:pr:`5130`, :issue:`4453`)

``discrete``
~~~~~~~~~~~~
- Improve performance by only computing matrix_rank(exog) once in DiscreteModel.initialize (:pr:`4805`)
- Improve performance in discrete models by avoiding repeated calculations (:pr:`4515`)
- Add ``cov_type`` to summary of discrete models (:pr:`5672`, :issue:`4581`)
- Add conditional multinomial logit (:pr:`5510`)
- Add conditional logistic and Poisson regression (:pr:`5304`)

``genmod``
~~~~~~~~~~
- Fix arguments in poisson version of ``BayesMixedLM`` (:pr:`4809`)
- Ensure that column names are properly attached to the model (:pr:`4788`)
- Change ``cov_params`` in ``BayesMixedLM`` to act more like it does in other models (:pr:`4788`)
- Add missing predict and fit methods to ``BayesMixedGLM`` (:pr:`4702`)
- Add influence function support for GLM (:pr:`4732`, :issue:`4268`, :issue:`4257`)
- Fixed a bug in GEE where history was not saved (:pr:`5789`)
- Enable ``missing='drop'`` in GEE (:pr:`5771`)
- Improve score test to allow the submodel to be provided as a GEEResults object instead of as linear constraints (:pr:`5435`)
- Use GLM to get starting values for GEE (:pr:`5440`)
- Added regularized GEE (:pr:`5450`)
- Added Generalized Additive Models (GAM) (:pr:`5481`, :issue:`5370`, :issue:`5296`, :issue:`4575`, :issue:`2744`, :issue:`2435`)
- Added tweedie log-likelihood (:pr:`5521`)
- Added ridge regression by gradient for all GLM (:pr:`5521`)
- Added Tweedie EQL quasi-likelihood (:pr:`5543`)
- Allow ``dep_data`` to be specified using formula or names (:pr:`5345`)
- Fix a bug in stationary cov_struct for GEE (:pr:`5390`)
- Add QIC for GEE (:pr:`4909`)

``graphics``
~~~~~~~~~~~~
- Allow QQ plots using samples with different sizes (:pr:`5673`, :issue:`2896`, :issue:`3169`)
- Added examples of many graphics functions to the documentation (:pr:`5607`, :issue:`5309`)
- Fixed a bug in ``interaction_plot`` which lost information in a ``pd.Series`` index (:pr:`5548`)
- Remove change of global pickle method in functional plots (:pr:`4963`)

``imputation``
~~~~~~~~~~~~~~
- Add formula support to MI multiple imputation (:pr:`4722`)
- Saves the column names from ``pd.DataFrames`` and returns the imputed results as a DataFrame in ``BayesMI`` (:pr:`4722`)
- Fixed warnings in ``MICEData`` related to setting on copy (:pr:`5606`, :issue:`5431`)
- Allow results to be stored for multiple imputation (:pr:`5093`)
- Fixed a bug where MICEData sets initial imputation incorrectly (:pr:`5301`, :issue:`5254`)

``iolib``
~~~~~~~~~
- Deprecate ``StataReader``, ``StataWriter``, and ``genfromdta`` in favor of pandas equivalents (:pr:`5770`)
- Improve string escaping when exporting to LaTeX (:pr:`5683`, :issue:`5297`)
- Fixed a bug in ``summary2`` that ignored user float formatting  (:pr:`5655`, :issue:`1964`, :issue:`1965`)
- Remove ``$$`` from LaTeX output (:pr:`5588`,:issue:`5444`)

``multivariate``
~~~~~~~~~~~~~~~~
- Fixed a bug that only allowed ``MANOVA`` to work correctly when called using the formula interface (:pr:`5646`, :issue:`4903`, :issue:`5578`)
- Fix pickling bug in ``PCA`` (:pr:`4963`)

``nonparametric``
~~~~~~~~~~~~~~~~~
- Added input protection ``lowess` to ensure ``frac`` is always in bounds. (:pr:`5556`)
- Add check of inputs in ``KernelReg`` (:pr:`4968`, :issue:`4873`)

``regression``
~~~~~~~~~~~~~~
- Fix bug in  random effects covariance getter for ``MixedLM`` (:pr:`4704`)
- Add exact diffuse filtering for ``RecursiveLS`` (:pr:`4699`)
- Add Gaussian process regression (:pr:`4691`)
- Add linear restrictions to ``RecursiveLS`` (:pr:`4133`)
- Added regression with quadratic inference functions :class:`~statsmodels.genmod.qif.QIF` (:pr:`5803`)
- Allow mediation to be used with MixedLM as a mediator and/or outcome model (:pr:`5489`)
- Add square root LASSO (:pr:`5516`)
- Add dimension reduction regression methods: ``SlicedInverseReg``, ``PHD`` and ``SAVE`` (:pr:`5518`)
- Increased the number of methods available to optimize ``MixedLM`` models (:pr:`5551`)
- Added label to R2 when model is uncentered (:pr:`5083`, :issue:`5078`)
- Allow several optimizers to be tried in sequence for MixedLM (:pr:`4819`)
- Fix bug in Recursive LS with multiple constraints (:pr:`4826`)
- Fix a typo in ``ColinearityWarning`` (:pr:`4889`, :issue:`4671`)
- Add a finite check for ``_MinimalWLS`` (:pr:`4960`)
- Fix definition of R2 in ``GLS`` (:pr:`4967`, :issue:`1252`, :issue:`1171`)
- Add Burgs algorithm for estimating parameters of AR models (:pr:`5016`)

``sandbox``
~~~~~~~~~~~
- Add copulas (:pr:`5076`)

``stats``
~~~~~~~~~
- Implements a simple method of moments estimator of a spatial covariance in ``kernel_covariance`` (:pr:`4726`)
- Fixed a bug in multiple function in ``~statsmodels.stats.moment_helpers`` which prevents in-place modification of inputs (:pr:`5671`, :issue:`3362`, :issue:`2928`)
- Fixed a bug in contingency tables where shift was not correctly applied (:pr:`5654`, :issue:`3603`, :issue:`3579`)
- Added White's two-moment specification test with null hypothesis of homoskedastic and correctly specified(:pr:`5602`, :issue:`4721`)
- Added adjusted p-values for Tukey's HSD (:issue:`5418`, :pr:`5625`)
- Fixed a bug in ``medcouple`` that produced the incorrect estimate when there are ties in the data (:pr:`5397`, :issue:`5395`)
- Combine the real and knockoff features in init (:pr:`4920`)
- Modifying exog in-place leads to incorrect scaling (:pr:`4920`)
- Add Provide Knockoff+ (guaranteed to control FDR but slightly conservative) as well as Knockoff FDR (:pr:`4920`)
- Add RegModelEffects allows the user to specify which model is used for parameter estimation (:pr:`4920`)

``tools``
~~~~~~~~~
- Fixed a bug in ``group_sums`` that raised ``NameError`` (:pr:`5127`)

``tsa``
~~~~~~~
- Fix k_params in seasonal MAs (:pr:`4790`, :issue:`4789`)
- Fix prediction index in VAR predict (:pr:`4785`, :issue:`4784`)
- Standardized forecast error in state space when using Cholesky methods with partial missing data (:pr:`4770`)
- Add and fix VARMAX trend, exog. timing and polynomial trends (:pr:`4766`)
- Fix bug in exact diffuse filtering in complex data type case (:pr:`4743`)
- SARIMAX warns for non-stationary starting params (:pr:`4739`)
- Make arroots and maroots have consistent return type (:pr:`4559`)
- Add exact diffuse initialization to state space models (:pr:`4418`, :issue:`4042`)
- Allow concentrating scale out of log-likelihood in state space models (:pr:`3480`)
- Fixed a bug in ``coint_johansen`` that prevented it from running with 0 lags (:pr:`5783`)
- Improved performance in ``kpss`` using ``np.sum`` (:pr:`5774`)
- Enforce maximum number of lags in ``kpss`` (:pr:`5707`)
- Add ``arma_innovations`` to compute the innovations from an ARMA process (:pr:`5704`)
- Limit maximum lag length in ``adfuller`` so that model can always be estimated (:pr:`5699`, :issue:`5432`, :issue:`3330`)
- Added automatic data-dependent lag length selection in ``kpss`` (:pr:`5670`, :issue:`2781`, :issue:`5522`)
- Fixed a bug in ``VARMAX`` where the wrong form of the intercept was used when creating starting values (:pr:`5652`, :issue:`5651`)
- Fixed a bug ``sirf_errband_mc`` (:pr:`5641`, :issue:`5280`)
- Clarified error when input to ARMA is not a 1-d array (:pr:`5640`, :issue:`2575`)
- Improved the numerical stability of parameter transformation in ARIMA estimation (:pr:`5569`)
- Fixed a bug in the acf of a ``VAR`` which produced incorrect values (:pr:`5501`)
- Expose additional alternative estimation methods in ``pacf`` (:pr:`5153`, :issue:`3862`)
- Removed original implementation of Kalman Filter in favor of Cythonized version in ``statsmodels.tsa.statespace`` (:pr:`5171`)
- Issue warning when using ``VARResults.cov_params`` that it will become a method in the future (:pr:`5244`)
- Fix a bug in statespace models' ``predict`` that would fail when using row labels (:pr:`5250`)
- Allow ``summary`` even if filter_results=None, which happens after ``save`` and ``load`` (:pr:`5252`)
- Fixed a bug in sequential simulation in models with state_intercept (:pr:`5257`)
- Add an analytic version of ``arma_acovf`` (:pr:`5324`)
- Add a fast ARMA innovation algorithm and loglike computation (:pr:`5360`)
- Fix a bug in the Initialization of simulation smoother with exact diffuse initialization (:pr:`5383`)
- Fix bug in simulation smoothed measurement disturbance with FILTER_COLLAPSED (:pr:`4810`, :issue:`4800`)
- Improve SARIMAX for time series close to non-stationary (:pr:`4815`)
- Use Cython to improve speed of Exponential Smoothing models (:pr:`4845`)
- Fix a bug in ``arma_order_selection`` when data is passed in as a list (:pr:`4890`, :issue:`4727`)
- Add explicit exceptions in ARMA/ARIMA forecast with missing or wrong exog (:pr:`4915`, :issue:`3737`)
- Remove incorrect endog from results if constraints (:pr:`4921`)
- Add ``nlag`` argument to ``acovf`` (:pr:`4937`)
- Set reasonable default lags for acf/pacf plots (:pr:`4949`)
- Add innovations algorithm to convert acov to MA (:pr:`5042`)
- Add and innovations filter to filter for observations in a MA (:pr:`5042`)
- Fix a bug in initialization when simulating in state space models (:pr:`5043`)

``maintenance``
~~~~~~~~~~~~~~~
- Switch to standard setup.py so that ``pip install statsmodels`` can succeed in an empty virtual environment
- General compatibility fixes for recent versions of numpy, scipy and pandas
- Added new CI using Azure Pipelines (:pr:`5617`)
- Enable linting on travis to ensure code is up to standards (:pr:`4820`)
- Add coverage for Cython code (:pr:`4871`)
- Improve import speed (:pr:`5831`)
- Make all version of docs available (:pr:`5879`)

bug-wrong
---------

A new issue label `type-bug-wrong` indicates bugs that cause that incorrect
numbers are returned without warnings.
(Regular bugs are mostly usability bugs or bugs that raise an exception for
unsupported use cases.)
`see tagged issues <https://github.com/statsmodels/statsmodels/issues?q=is%3Aissue+label%3Atype-bug-wrong+is%3Aclosed+milestone%3A0.10/>`_

- :issue:`5475`
- :issue:`5316`


Major Bugs Fixed
================

See github issues for a list of bug fixes included in this release

- `Closed bugs <https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.10+label%3Atype-bug/>`_
- `Closed bugs (wrong result) <https://github.com/statsmodels/statsmodels/pulls?q=is%3Apr+is%3Amerged+milestone%3A0.10+label%3Atype-bug-wrong/>`_


Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Chad Fulton
* Brock Mendel
* Peter Quackenbush
* Kerby Shedden
* Kevin Sheppard

and the general maintainer and code reviewer

* Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.10 release (based on git log):


* Amir Masoud Abdol
* Andrew Davis
* Andrew Kittredge
* Andrew Theis
* bertrandhaut
* bksahu
* Brock Mendel
* Chad Fulton
* Chris Snow
* Chris Down
* Daniel Saxton
* donbeo
* Emlyn Price
* equinaut
* Eric Larson
* Evgeny Zhurko
* fourpoints
* Gabriel Reid
* Harry Moreno
* Hauke Jürgen Mönck
* Hugo
* hugovk
* Huize Wang
* JarnoRFB
* Jarrod Millman
* jcdang
* Jefferson Tweed
* Josef Perktold
* jtweeder
* Julian Taylor
* Kerby Shedden
* Kevin Sheppard
* Loknar
* Matthew Brett
* Max Ghenis
* Ming Li
* Mitch Negus
* Michael Handley
* Moritz Lotze
* Nathan Perkins
* Nathaniel J. Smith
* Niklas H
* Peter Quackenbush
* QuentinAndre
* Ralf Gommers
* Rebecca N. Palmer
* Rhys Ulerich
* Richard Barnes
* RonLek
* Stefaan Lippens
* Tad seldovia
* thequackdaddy
* Tom Augspurger
* Torsten Wörtwein
* Varanelli
* xrr
* Yichuan Liu
* zveryansky
* 郭飞

These lists of names are automatically generated based on git log, and may not
be complete.

Merged Pull Requests
--------------------

Thie following Pull Requests were merged since the last release:


* :pr:`2096`: Score/LM conditional moment tests
* :pr:`3480`: ENH: State space: allow concentrating scale out of log-likelihood
* :pr:`4048`: Remove redundant code for dropped Python 2.6
* :pr:`4133`: ENH: Add linear restrictions to RecursiveLS
* :pr:`4316`: ensure MultinomialResults has J, K.  Get rid of unnecessary lmap usage
* :pr:`4322`: Make DiscreteResults Unchanging
* :pr:`4371`: catch the correct exception, make assertions not-pointless
* :pr:`4418`: ENH: State space: Exact diffuse initialization
* :pr:`4458`: De-duplicate a bunch of identical code
* :pr:`4468`: remove unused resetlist
* :pr:`4487`: Get rid of non-standard imports and one-line functions
* :pr:`4494`: Fix imports math.foo -->np.foo in vecm
* :pr:`4501`: xfail test instead of commenting it out
* :pr:`4515`: PERF: Simplify algebra in discrete_model
* :pr:`4559`: REF: make arroots and maroots have consistent return type
* :pr:`4560`: Document and cleanup bits of cython code
* :pr:`4576`: Penalized mle scad rebased2
* :pr:`4593`: DOC:ArmaProcess class documentation typo fix
* :pr:`4594`: TEST/DOC: SMW linalg routines documentation and test
* :pr:`4640`: BF: DataTimeIndex.to_datetime removed in pandas
* :pr:`4648`: BUG/TEST: Make pattern order for multiple imputation deterministic
* :pr:`4650`: DISCUSS/BLD: Update minimum versions.
* :pr:`4653`: REF/MAINT: avoid dict with pandas
* :pr:`4658`: BLD: Use older version of Pandas for docbuild
* :pr:`4683`: ENH: add ultra-high screening with SCAD
* :pr:`4686`: TEST: Docstring edits and variable name changes for clarity
* :pr:`4689`: PERF: Declare temporary output for hessian
* :pr:`4691`: ENH: Gaussian process regression
* :pr:`4692`: DOC: Add GLM varfuncs and weights notebook to documentation
* :pr:`4696`: Configure doctr
* :pr:`4698`: REF: Remove compatibility mode for state space
* :pr:`4699`: ENH: Exact diffuse filtering for RecursiveLS
* :pr:`4702`: BUG: Add missing predict and fit methods to BayesMixedGLM
* :pr:`4704`: BUG Fix random effects covariance getter for MixedLM
* :pr:`4712`: BUG: add ModelWarning base class to avoid warning filter on standard UserWarning.
* :pr:`4717`: TST: allclose instead of exact match for floats and use machine precision
* :pr:`4720`: fix syntax-like error
* :pr:`4722`: ENH: Add formula support to MI multiple imputation
* :pr:`4726`: ENH Kernel covariance
* :pr:`4728`: TST: Openblas appveyor fixes
* :pr:`4732`: ENH: add GLMInfluence
* :pr:`4736`: DOC: Make custom function take effect
* :pr:`4739`: REF: SARIMAX: only warn for non stationary starting params
* :pr:`4743`: BUG: state space: exact diffuse filtering in complex data type case
* :pr:`4750`: DOC: Fix indentation of math formulas
* :pr:`4753`: DOC: Add notebook on concentrated scale in ssm
* :pr:`4758`: DOC: Added missing notebooks to examples
* :pr:`4760`: CLN: Provide better name for pooled risk ratio
* :pr:`4763`: replace copy/pasted code with import
* :pr:`4766`:  BUG/ENH: VARMAX Fix trend / exog. timing. Add polynomial trends.
* :pr:`4767`: MAINT: gitignore univariate_diffuse pyx files.
* :pr:`4770`: BUG: State space: standardized forecast error when using Cholesky methods with partial missing data
* :pr:`4777`: MAINT: conda specify numpy-base
* :pr:`4785`: BUG: Get prediction index in VAR predict.
* :pr:`4786`: CLEAN: fix indentation by four typos
* :pr:`4788`: BUG: bayes mixed GLM maintenance
* :pr:`4790`: BUG: k_params if seasonal MA
* :pr:`4805`: Only compute matrix_rank(exog) once in DiscreteModel.initialize
* :pr:`4809`: BUG: fix arguments in poisson mixed model
* :pr:`4810`: BUG: simulation smoothed measurement disturbance with FILTER_COLLAPSED
* :pr:`4814`: CLEAN: Removed unnecessary and non-informative print
* :pr:`4815`: ENH/BUG: Improve SARIMAX for time series close to non-stationary
* :pr:`4819`: ENH: Allow several optimizers to be tried in sequence for MixedLM
* :pr:`4820`: Implement basic linting for Travis
* :pr:`4823`: Fix deprecation warnings
* :pr:`4826`: BUG/ENH: Recursive LS: fix bug w/ multiple constraints
* :pr:`4834`: Implement full flake8 checking for a subset of files in good condition
* :pr:`4835`: CLEAN: Fix tab indentation, lint for it
* :pr:`4842`: CLN: Flake8 fixups and linting for statespace files (but not tests)
* :pr:`4844`: CLN: Fully lint regime_switching
* :pr:`4845`: ENH: Improve speed in Exponential Smoothing
* :pr:`4853`: CLN/REF: Remove recarrays from datasets
* :pr:`4855`: BUG: Attach vc_names for mixed Poisson models
* :pr:`4858`: MAINT: Delete migrate_issues_gh
* :pr:`4859`: Fix some NameErrors, do not delete unused [...]
* :pr:`4861`: DOC: Fix small doc errors
* :pr:`4864`: CLN: fix and lint for W391 blank line at end of file
* :pr:`4869`: Update setup.cfg
* :pr:`4871`: BLD: Refactor Setup
* :pr:`4872`: MAINT: Remove nose and related references
* :pr:`4879`: CLN: Fix documentation for Levinson-Durbin
* :pr:`4883`: CLN: remove empty __main__ sections
* :pr:`4886`: CLN: Fully lint recursive_ls.py
* :pr:`4889`: REF: Rename ColinearityWarning
* :pr:`4890`: BUG: Add check to ensure array in arma order selection
* :pr:`4891`: BLD: Fix linting and move coverage
* :pr:`4893`: TST: Restore incorrectly disabled test
* :pr:`4895`: CLN: Fix and lint for misleading indentation E125,E129
* :pr:`4896`: CLN: Fix and lint for potential double-negatives E713,E714
* :pr:`4897`: CLN: Fix and lint for multiple spaces after keyword E271
* :pr:`4900`: CLN: Lint for missing whitespace around modulo operator E228,E401
* :pr:`4901`: CLN: Fix and lint for E124 closing bracket does not match visual indentation
* :pr:`4909`: ENH: QIC for GEE
* :pr:`4910`: CLN: Blank Lines E301,E302,E303,E305,E306 in examples, tools, sm.base
* :pr:`4911`: MAINT: Remove future errors and warnings
* :pr:`4912`: BLD: Rebased ci improvements
* :pr:`4913`: TST: Add a fixture to close all plots
* :pr:`4914`: CLN: Blanks E301,E302,E303,E305,E306 in tsa
* :pr:`4915`: ENH: explicit exceptions in ARMA/ARIMA forecast with missing or wrong exog
* :pr:`4920`: BUG/ENH: Two bug fixes and several enhancements to knockoff filter (regression fdr)
* :pr:`4921`: BUG: remove faux endog from results if constraints
* :pr:`4924`: CLN: E242 space after tab, enforce all passing rules
* :pr:`4925`: CLN: Enforce E721, use isinstance
* :pr:`4926`: CLN: Enforce E306, blank lines in nested funcs
* :pr:`4927`: CLN: Enforce E272, multiple spaces
* :pr:`4929`: BLD: Add linting for any new files
* :pr:`4933`: Remove unused patsy import in quantile_regression.ipynb
* :pr:`4937`: ENH: Add nlag argument to acovf
* :pr:`4941`: MAINT: remove exact duplicate file datamlw.py
* :pr:`4943`: TST: Relax tolerance on failing test
* :pr:`4944`: BLD: Add pinned numpy on appveyor
* :pr:`4949`: BUG: Set default lags for acf/pacf plots
* :pr:`4950`: DOC: Fix small typo in unit root testing example
* :pr:`4953`: DOC: Fix nagging issues in docs
* :pr:`4954`: BUG: disallow use_self=False
* :pr:`4959`: DOC: Clean up tsa docs
* :pr:`4960`: BUG: Add finite check for _MinimalWLS
* :pr:`4963`: BUG: Remove change of global pickle method
* :pr:`4967`: BUG: Fix definition of GLS r2
* :pr:`4968`: BUG: Check inputs in KernelReg
* :pr:`4971`: DOC: Switch to https where used
* :pr:`4972`: MAINT/CLN Remove .bzrignore
* :pr:`4977`: [BUG/MAINT] Fix NameErrors caused by missing kwargs
* :pr:`4978`: [MAINT/Test] skip test instead of mangling name in test_generic_methods
* :pr:`4979`: [MAINT/TST] remove np.testing.dec unused imports (nose dependency)
* :pr:`4980`: [MAINT/TST] skip/xfail tests instead of mangling/commenting-out in genmod, regression
* :pr:`4981`: [MAINT] Remove info.py
* :pr:`4982`: DOC Fix typo Parameters-->Parameters
* :pr:`4983`: [TST] xfail/skip instead of commenting-out/mangling discrete tests
* :pr:`4984`: [TST/DOC] make commented-out code in tests/results into readable docs
* :pr:`4985`: [TST/DOC] Make test comments more readable
* :pr:`4986`: [MAINT/TST] turn commented-out code into readable docs in results_arma
* :pr:`4987`: [TST/MAINT] turn commented-out code into readable docs in results_ar,…
* :pr:`4988`: [TST/MAINT] de-duplicate get_correction_factor code
* :pr:`4989`: [MAINT/CLN] Remove code made unusable due to license issue
* :pr:`4990`: [MAINT/CLN] remove numdiff  __main__ section explicitly marked as scratch work
* :pr:`4993`: [TST/CLN] Turn decorators __main__ section into tests
* :pr:`4995`: [TST] make tools.linalg __main__ section into tests
* :pr:`4998`: [CLN/TST] Follow instructions to remove function
* :pr:`4999`: [MAINT] remove wrappers.py
* :pr:`5000`: [MAINT] update compat to remove unusable shims e.g. py26
* :pr:`5002`: [MAINT] add missing import
* :pr:`5003`: MAINT: fix invalid exception messages
* :pr:`5005`: [MAINT] remove unused imports in examples+tools
* :pr:`5007`: MAINT: unused imports in robust
* :pr:`5011`: [MAINT] remove text file relics from scikits/statsmodels
* :pr:`5012`: [MAINT/TST] move misplaced results files in regressions/tests
* :pr:`5013`: [MAINT] fix typo deprecated-->deprecated
* :pr:`5014`: [MAINT] typo in __init__ signature
* :pr:`5015`: [MAINT] move misplaced test_tsa_indexes
* :pr:`5016`: ENH: Burgs algorithm
* :pr:`5020`: MAINT: fix incorrect docstring summary-->summary2
* :pr:`5021`: MAINT: fix typo duplicated References in docstring
* :pr:`5024`: MAINT: silenced as_pandas warnings in documentation
* :pr:`5027`: MAINT: remove functions duplicated from scipy
* :pr:`5029`: MAINT: strict linting for sm.stats files _close_ to already passing
* :pr:`5040`: MAINT: clean up x13.py, delete main
* :pr:`5042`: ENH: Add innovations algorithm
* :pr:`5043`: BUG: Initialization when simulating
* :pr:`5045`: MAINT: strict linting for tsa.statespace.tests.results
* :pr:`5057`: BUG: Correct check for callable
* :pr:`5058`: BUG: Do not use mutable default values
* :pr:`5059`: BLD: Add line displaying CPU info to CI
* :pr:`5065`: TST: Fix incorrect assertion
* :pr:`5070`: MAINT: remove file that just says to remove it
* :pr:`5071`: MAINT: remove example file corresponding to removed module
* :pr:`5074`: MAINT: strict lint test_var.py
* :pr:`5075`: MAINT: strict linting test_univariate.py
* :pr:`5076`: ENH: more work on copula (deriv, classes)
* :pr:`5079`: MAINT: linting statespace tests
* :pr:`5080`: FIX failure caused by #5076
* :pr:`5083`: ENH: Add "(uncentered)" after rsquared label in .summary, .summary2 when appropriate
* :pr:`5086`: TST: parametrize tests instead of using for loops
* :pr:`5088`: DOC: Add javascript to link to other doc versions
* :pr:`5090`: MAINT: Chrome does not like having a secure link with an unsecure image
* :pr:`5093`: Allow results to be stored for multiple imputation
* :pr:`5096`: ENH remove unneeded restriction on QIC (GEE)
* :pr:`5099`: MAINT: fix and lint for W292 newline at end of file
* :pr:`5103`: BUG: fix missing new_branch_dir arg in upload_pdf
* :pr:`5105`: BUG/DOC: Description of k_posdef
* :pr:`5114`: MAINT: many but not all trailing whitespace
* :pr:`5119`: CLN: remove unused imports in tools, sm.tsa
* :pr:`5120`: BUG: Ensure internal tester exits with error if needed
* :pr:`5121`: MAINT: Avoid star imports
* :pr:`5122`: MAINT: Modernize R-->py script, lint output
* :pr:`5123`: CLN: Move results files to a location that is copied to install
* :pr:`5124`: MAINT: fix generated double whitespace
* :pr:`5127`: BUG: Fix NameError in grouputils, make __main__ into tests
* :pr:`5130`: BUG: incorrect self.score in GenericLikelihoodModel; closes #4453
* :pr:`5133`: TST: apply stacklevel to warning in Family.__init__
* :pr:`5135`: MAINT: Fix warnings
* :pr:`5136`: TST: improve testing util functions; de-duplicate
* :pr:`5138`: CLN: Use cache_readonly instead of OneTimeProperty
* :pr:`5141`: MAINT: Delete bspline source files
* :pr:`5143`: ENH/BUG Bootstrap clone rebased
* :pr:`5146`: Clean up the smf namespace
* :pr:`5148`: REF/TST: add seed to hdrboxplot, Use random order in pytest
* :pr:`5149`: TST: Theil test randomseed
* :pr:`5152`: REF: Use iterative cumsum_n
* :pr:`5153`: ENH: Add additional options for pacf ols
* :pr:`5156`: TST: Remove __main__ sections in tests
* :pr:`5162`: TST: Fix incorrect test closes #4325
* :pr:`5164`: BF: drop tolerance of a zero_constrained test
* :pr:`5165`: MAINT: Add decorator for tests that use matplotlib
* :pr:`5166`: DOC: Fix section title in QIC
* :pr:`5167`: TST/BUG: Fix missing SkipTest
* :pr:`5170`: DEPR: Remove items deprecated in previous versions
* :pr:`5171`: MAINT: Remove kalmanf StateSpace code supplanted by tsa.statespace
* :pr:`5176`:  TST: Fix random generation issue
* :pr:`5177`: DOC: Improve Holt Winters documentation
* :pr:`5178`: TST: Fix scale in test
* :pr:`5180`: TST: Change assert_approx_equal to assert_allclose
* :pr:`5184`: TST: parametrize tests in test_lme
* :pr:`5188`: BLD/TST: Add coverage for Cython files
* :pr:`5191`: MAINT: Remove selected __main__ sections
* :pr:`5192`: MAINT: Fix incorrect pass statements
* :pr:`5193`: MAINT: raise specific exceptions instead of just Exception
* :pr:`5194`: MAINT: fix incorrect TypeError --> ValueError
* :pr:`5195`: BLD: Include License in Wheel
* :pr:`5196`: TST: Set seed when using basin hopping
* :pr:`5198`: TST/CLN/BUG: Fix corr nearest factor
* :pr:`5200`: TST: Alter test condition due to parameter scale
* :pr:`5201`: TST/CLN: test_arima_exog_predict, Rescale data to avoid convergence issues
* :pr:`5203`: BUG: raise instead of return ValueError
* :pr:`5204`: MAINT: Avoid/Fix FutureWarnings
* :pr:`5207`: TST: Ensure random numbers are reproducible
* :pr:`5208`: TST/CLN: Tighten tol to reduce spurious test failure
* :pr:`5210`: BLD: Ensure main is available when linting
* :pr:`5211`: MAINT: Import instead of copy/pasting utils
* :pr:`5213`: MAINT: Move misplaced duration results files
* :pr:`5214`: MAINT: remove example-like file that could never run
* :pr:`5217`: MAINT: Remove outdated pandas compat shims
* :pr:`5218`: MAINT: Move misplaced genmod results files
* :pr:`5219`: fixed typo
* :pr:`5222`: MAINT: fully lint formula
* :pr:`5223`: MAINT: fully lint compat
* :pr:`5224`: REF: raise early on invalid method
* :pr:`5227`: MAINT: docstring and whitespace fixups
* :pr:`5228`: DOC: Fix many small errors in examples
* :pr:`5230`: DOC: Fix small doc build errors
* :pr:`5232`: TST: mark smoketests
* :pr:`5237`: TST: Add mac testing
* :pr:`5239`: BLD/TST: Add platform-specific skips to CI testing
* :pr:`5240`: MAINT: remove cythonize.py made unnecessary by #4871
* :pr:`5242`: DOC: Update release instructions [skip ci]
* :pr:`5244`: DEPR: warn that VARResults.cov_params will become method
* :pr:`5246`: DOC: Added documentation of elements in anova_lm
* :pr:`5248`: DOC: Revert incorrect docstring change [skip ci]
* :pr:`5249`: MAINT: Add Script to convert notebooks
* :pr:`5250`: BUG/TST: TSA models: _get_index_label_loc failed when using row labels.
* :pr:`5251`: DOC: Use correct “autoregressive” in docstring
* :pr:`5252`: BUG: Allow `summary` even if filter_results=None (e.g. after `save`, `load`
* :pr:`5257`: BUG: Sequential simulation in models with state_intercept
* :pr:`5260`: MAINT: avoid pandas FutureWarning by checking specific condition
* :pr:`5262`: MAINT: fix typos in pca, wrap long lines
* :pr:`5263`: BLD: Only unshallow when required
* :pr:`5265`: MAINT: Prefer signature over formatargspec
* :pr:`5267`: MAINT: implement _wrap_derivative_exog for de-duplication
* :pr:`5269`: MAINT: De-duplicate code in iolib.summary
* :pr:`5272`: WIP/MAINT: Identify defunct code in summary methods
* :pr:`5273`: Fix incorrect parameter name in docstring
* :pr:`5274`: MAINT: remove self.table pinning
* :pr:`5275`: ENH/BUG Modify GEE indexing to remove numpy warnings
* :pr:`5277`: DOC: Clarify/fix docs on GLM scale estimation for negative binomial
* :pr:`5292`: DOC: Remove only_directive
* :pr:`5295`: TST: Added random seed to test_gee and verified working
* :pr:`5300`: DOC fix docstring in stattools.py
* :pr:`5301`: BUG: MICEData sets initial imputation incorrectly
* :pr:`5304`: ENH: conditional logistic and Poisson regression
* :pr:`5306`: DOC: Workarounds to fix docbuild
* :pr:`5308`: REF: Collect covtype descriptions, de-duplicate normalization func
* :pr:`5314`: DOC: minor fix on documentation on Durbin Watson test
* :pr:`5322`: DOC: Move magic
* :pr:`5324`: ENH: analytic version of arma_acovf
* :pr:`5325`: BUG/TST: Fix innovations_filter, add test vs Kalman filter
* :pr:`5335`: MAINT: eliminate some pytest warnings
* :pr:`5345`: ENH: Allow dep_data to be specified using formula or names
* :pr:`5348`: Set python3 as interpreter for doc tools
* :pr:`5352`: CLN: Fix F901 and E306 mixups
* :pr:`5353`: CLN: W605 fixups in vector_ar
* :pr:`5359`: BUG: raise correct error
* :pr:`5360`: ENH: Fast ARMA innovation algorithm and loglike computation
* :pr:`5369`: MAINT: disable pytest minversion check (broken in pytest 3.10.0)
* :pr:`5383`: BUG: Initialization of simulation smoother with exact diffuse initialization
* :pr:`5390`: BUG/ENH: modify stationary cov_struct for GEE
* :pr:`5397`: BUG: Fix medcouple with ties
* :pr:`5399`: CLN: Fix some invalid escapes
* :pr:`5421`: CLN: informative names for test functions
* :pr:`5424`: MAINT: conda-forge use gcc7
* :pr:`5426`: Misspelling in the documentation proportions_ztest
* :pr:`5435`: ENH Score test enhancements for GEE
* :pr:`5440`: ENH: Use GLM to get starting values for GEE
* :pr:`5449`: ENH/DOC: Added linting instruction in CONTRIBUTING.rst
* :pr:`5450`: ENH: regularized GEE
* :pr:`5462`: Fixed broken link for Guerry Dataset
* :pr:`5471`: Fix broken link
* :pr:`5481`: ENH: Generalized Additive Models and splines (Gam 2744 rebased4)
* :pr:`5484`: DOC: fix gam.rst
* :pr:`5485`: MAINT: Travis fixes
* :pr:`5489`: ENH: Mediation for Mixedlm
* :pr:`5494`: BUG: Bad escapes
* :pr:`5497`: Fix typo in docstring
* :pr:`5501`: BUG: Correct error in VAR ACF
* :pr:`5510`: ENH Conditional multinomial logit
* :pr:`5513`: DOC: Fix spelling
* :pr:`5516`: ENH square root lasso
* :pr:`5518`: ENH dimension reduction regression
* :pr:`5521`: ENH: Tweedie log-likelihood (+ridge regression by gradient for all GLM)
* :pr:`5532`: DOC/ENH Docstring updates for clogit
* :pr:`5541`: DOC: Describe binomial endog formats
* :pr:`5542`: BUG/TEST: py27 needs slacker tolerances
* :pr:`5543`: BUG: Tweedie EQL quasi-likelihood
* :pr:`5548`: keep index of series when recoding a series
* :pr:`5551`: ENH: extend mixedlm optimizer attempts
* :pr:`5556`: Update _smoothers_lowess.pyx
* :pr:`5566`: Add project_urls to setup
* :pr:`5567`: Correct a spell mistake
* :pr:`5569`: ENH: Improve numerical stability of _ar_transparams, _ar_invtransparams
* :pr:`5582`: Jbrockmendel w605b
* :pr:`5583`: MAINT: Set language level for Cython
* :pr:`5584`: MAINT: Remov deprecation issues
* :pr:`5586`: DOC: Add issue and pr templates
* :pr:`5587`: MAINT: Resolve additional deprecations
* :pr:`5588`: BUG: Replace $$ in generated LaTeX
* :pr:`5589`: DOC: Updated the `for all i, j`
* :pr:`5590`: MAINT: Reorder travis so that legacy fails early
* :pr:`5591`: Jbrockmendel manywarns3
* :pr:`5592`: Jbrockmendel pversion
* :pr:`5593`: MAINT: remove never-needed callable and never-used compat functions
* :pr:`5594`: TST: Ensure test is identical in all runs
* :pr:`5595`: MAINT: Remove warnings from tests
* :pr:`5596`: TST: Explicitly set seed in basinhopping
* :pr:`5597`: MAINT: Remove unavailable imports
* :pr:`5599`: DOC: More emphasis and fix reference
* :pr:`5600`: TST: Relax tolerance for OpenBlas issue
* :pr:`5601`: Update mixed_linear.rst
* :pr:`5602`: ENH: White spec test (clean commit for PR 4721)
* :pr:`5604`: MAINT: Update template to encourage main check
* :pr:`5605`: Guofei9987 modify comments proportion confint
* :pr:`5606`: Mattwigway mice setting with copy warning
* :pr:`5607`: Jtweeder graphics addgraphics
* :pr:`5611`: BUG: Stop hardcoding parameters in results
* :pr:`5612`: MAINT: Ensure no warnings are produced by foreign
* :pr:`5613`: DOC: Improve PR template [skip ci]
* :pr:`5614`: MAINT: Deprecate scale in test function
* :pr:`5615`: Thequackdaddy docs
* :pr:`5616`: Bulleted list and minor typos in ttest_ind
* :pr:`5617`: CI: Implement azure-pipelines with multi-platform support
* :pr:`5621`: CLN: simplify lint configuration, fix some invalid escapes
* :pr:`5622`: DOC: Restore import
* :pr:`5625`: Andrew d davis tukey pvals
* :pr:`5626`: MAINT: Improve user-facing error message
* :pr:`5627`: BLD: Remove redundant travis config
* :pr:`5628`: MAINT: Relax tolerance on OSX only
* :pr:`5630`: MAINT: Enable xdist on azure
* :pr:`5631`: MAINT: Allow webuse fail
* :pr:`5633`: TST: change skip to xfail for test_compare_numdiff on OSX
* :pr:`5634`: Gabrielreid pandas multiindex handling bug
* :pr:`5635`: MAINT: Add a codecov config file
* :pr:`5636`: DOC: Update badges [skip ci]
* :pr:`5637`: CLN: strict linting for tools directory
* :pr:`5638`: MAINT: remove file with note to remove in 0.5.0
* :pr:`5640`: ENH: Improve error when ARMA endog is not 1d
* :pr:`5641`: Josef pkt svar irf errband 5280
* :pr:`5642`: TST: Relax tolerance on OSX for OpenBlas issues
* :pr:`5643`: MAINT: Consolidate platform checks
* :pr:`5644`: CLN/DOC: Remove unused module, vbench references
* :pr:`5645`: TST: Allow network failure in web tests
* :pr:`5646`: BUG: Fix MANOVA when not using formulas
* :pr:`5647`: TST: Adjust test_irf atol
* :pr:`5648`: BUG: Replace midrule with hline
* :pr:`5649`: CLN: strict linting for robust/tests directory
* :pr:`5650`: MAINT: Fix error in lint script
* :pr:`5652`: ENH/BUG: Use intercept form of trend / exog in VARMAX start params (not mean form)
* :pr:`5653`: MAINT: Reformat exceptions
* :pr:`5654`: Evgenyzhurko fix contingency table
* :pr:`5655`: BUG: summary2 use float_format when creating `_simple_tables` see #1964
* :pr:`5656`: BLD: Add linting to azure
* :pr:`5657`: TST: Protect multiprocess using test
* :pr:`5658`: BLD: Match requirements in setup and requirements
* :pr:`5659`: TST: Allow corr test to fail on Win32
* :pr:`5660`: MAINT: Fix make.bat [skip ci]
* :pr:`5661`: TST: Relax test tolerance on OSX
* :pr:`5662`: TST: Protect multiprocess on windows
* :pr:`5663`: MAINT: Add test runners
* :pr:`5664`: CLN: Fix and lint for E703 statements ending in semicolon
* :pr:`5666`: TST: Relax tolerance for irf test on windows
* :pr:`5667`: TST: Adjust tol and reset random state
* :pr:`5668`: TST: Adjust tolerance for test on windows
* :pr:`5669`: MAINT: Remove unused code
* :pr:`5670`: Jim varanelli issue2781
* :pr:`5671`: BUG: fix stats.moment_helpers inplace modification
* :pr:`5672`: ENH: Add cov type to summary for discrete models
* :pr:`5673`: ENH: Allow comparing two samples with different sizes
* :pr:`5675`: CLN: strict linting for emplike/tests
* :pr:`5679`: DOC: Clarify that predict expects arrays in dicts [skip ci]
* :pr:`5680`: ENH: Allow const idx to be found
* :pr:`5681`: BUG: Always set mle_retvals
* :pr:`5683`: BUG: Escape strings for latex output
* :pr:`5684`: BUG: fix df in summary for single constraint in wald_test_terms
* :pr:`5685`: Spelling
* :pr:`5686`: DOC: Fix parameter description in weightstats
* :pr:`5691`: MAINT: Near-duplicate example file, remove dominated version
* :pr:`5693`: CLN: Fix invalid escapes where possible
* :pr:`5694`: MAINT: Fix NameErrors in correlation_structures
* :pr:`5695`: MAINT: remove NameError-having version of levinson_durbin, just keep …
* :pr:`5696`: CLN: remove identical functions from garch
* :pr:`5697`: CLN: strict linting for examples/
* :pr:`5698`: PERF: Avoid implicit check when hasconst
* :pr:`5699`: BUG: Limit lag length in adf
* :pr:`5700`: MAINT: Update import of URLError
* :pr:`5701`: MAINT: missing imports, typos, fixes several NameErrors
* :pr:`5702`: MAINT: clean up docstring'd-out failure in __main__ block
* :pr:`5703`: MAINT: confirm that docstring'd-out traceback no longer raises; remove
* :pr:`5704`: ENH: expose innovations computation method to API.
* :pr:`5705`: WIP: TST: Sort dicts in test_multi
* :pr:`5707`: ENH: KPSS - detailed error message when lags > nobs
* :pr:`5709`: TST: Fix bad bash
* :pr:`5710`: CLN: clean up over/under indentation in tsa.tests.results, E12 codes
* :pr:`5712`: CLN: fix invalid escapes in test_stattools introduced in #5707
* :pr:`5713`: CLN/EX: Troubleshoot broken example, clean up now-working scratch paper
* :pr:`5715`: CLN: ellipses-out invalid escapes traceback
* :pr:`5716`: MAINT: Fix incorrect specification of loglike arg
* :pr:`5717`: MAINT: fix non-working example ex_pandas
* :pr:`5720`: CLN: remove impossible commented-out imports, close several
* :pr:`5721`: CLN: strict linting for dimred, processreg, and their tests.
* :pr:`5723`: Spelling fix in ValueError message
* :pr:`5724`: MAINT: close assorted small issues
* :pr:`5726`: DOC: Remove redundant attributes in GLM
* :pr:`5728`: CLN: remove and lint for unused imports
* :pr:`5729`: MAINT: use dummy_sparse func within method, see GH#5687
* :pr:`5730`: CLN: strict linting for discrete.tests.results
* :pr:`5732`: CLN: strict linting for genmod/tests/results
* :pr:`5734`: CLN: codes with only a few violations apiece
* :pr:`5736`: CLN: strict linting for regression/tests/results
* :pr:`5737`: CLN: strict linting for tsa.filters
* :pr:`5738`: CLN: strict linting for stats/tests/results/
* :pr:`5740`: CLN: strict linting for tsa.tests.results
* :pr:`5742`: CLN: strict linting for remaining results directories
* :pr:`5743`: CLN: strict linting for results files in sandbox/regression/tests/
* :pr:`5744`: CLN: Fix/lint for dangerous redefinitions and comparisons
* :pr:`5746`: MAINT: fix missing or redundant imports
* :pr:`5748`: CLN: clean up adfvalues, avoid using `eval`
* :pr:`5750`: CLN: E131 hanging indentation alignment
* :pr:`5758`: CLN: lint for ambiguous variable names
* :pr:`5760`: TST: test for intentionally emitted warnings, avoid some unintentional ones
* :pr:`5762`: BUG: rename wts to weights issue #4725
* :pr:`5765`: BUG/TST: Fix+test pieces of code that would raise NameError
* :pr:`5770`: DEPR: deprecate StataReader, StataWriter, genfromdta
* :pr:`5771`: ENH: improve missing data handling for GEE
* :pr:`5774`: PERF: use np.sum(...) instead of sum(...)
* :pr:`5778`: CLN: strict linting for test_varmax
* :pr:`5780`: TST: Protext against SSLError
* :pr:`5781`: CLN: Replace #5779
* :pr:`5783`: BUG: Ensure coint_johansen runs with 0 lags
* :pr:`5789`: BUG: GEE fit_history
* :pr:`5791`: Holder bunch
* :pr:`5792`: MAINT: matplotlib normed -> density
* :pr:`5793`: MAINT: Adjust tolerance for random fail on OSX
* :pr:`5796`: CLN: test_data.py
* :pr:`5798`: BUG: ignore bugs instead of fixing them
* :pr:`5801`: CI: Consolidate coveragerc spec
* :pr:`5803`: ENH: QIF regression
* :pr:`5805`: REF/CLN: collect imports at top of file, de-duplicate imports
* :pr:`5815`: CLN: test_gee.py
* :pr:`5816`: CLN: genmod/families/
* :pr:`5818`: CLN: qif
* :pr:`5825`: MAINT: use correct key name to check cov params presence
* :pr:`5830`: DOC: Add docstring for base class
* :pr:`5831`: PERF: Import speed
* :pr:`5833`: BUG: ARIMA fit with trend and constant exog
* :pr:`5834`: DOC: Fix small errors in release notes
* :pr:`5839`: MAINT: RangeIndex._start deprecated in pandas 0.25
* :pr:`5836`: CLN: over-indentation E117
* :pr:`5837`: CLN: invalid escapes in linear_model
* :pr:`5843`: MAINT: Catch intentional warnings
* :pr:`5846`: DOC: Update maintainer
* :pr:`5847`: BUG: Allow NumPy ints #
* :pr:`5848`: BUG: Warn rather than print
* :pr:`5850`: MAINT: Improve error message
* :pr:`5851`: BUG: Refactor method used to name variables
* :pr:`5853`: BUG: Add check for xnames length
* :pr:`5854`: BUG: Fix MNLogit summary with float values
* :pr:`5857`: BUG: Allow categorical to accept pandas dtype
* :pr:`5858`: BUG: Fix default alignment for SimpleTable
* :pr:`5859`: DOC: fix incorrect ARResults.predict docstring, closes #4498
* :pr:`5860`: Cdown gofplot typerror
* :pr:`5863`: MAINT: Use pd.Categorical() instead of .astype('categorical')
* :pr:`5868`: BUG: State space univariate smoothing w/ time-varying transition matrix: wrong transition matrix used
* :pr:`5869`: DOC: Improve ExponentialSmoothing docstring
* :pr:`5875`: DOC: Improve bug report template
* :pr:`5876`: BUG: Ensure keywords exist in partial reg plot
* :pr:`5879`: DOC: Update version dropdown javascript
:orphan:

==============
Release 0.10.2
==============

Release summary
===============
This is a bug release and adds compatibility with Python 3.8.

Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Chad Fulton
* Qingqing Mao
* Diego Mazon
* Brock Mendel
* Guglielmo Saggiorato
* Kevin Sheppard
* Tim Staley

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

* :pr:`5935`: CLN: port parts of #5220
* :pr:`5951`: BUG: Fix mosaic plot with missing category
* :pr:`5996`: BUG: Limit lags in KPSS
* :pr:`5998`: Replace alpha=0.05 with alpha=alpha
* :pr:`6030`: Turn relative import into an absolute import
* :pr:`6044`: DOC: Fix notebook due to pandas index change
* :pr:`6046`: DOC: Remove DynamicVAR
* :pr:`6091`: MAINT/SEC: Remove unnecessary pickle use
* :pr:`6092`: MAINT: Ensure r download cache works
* :pr:`6093`: MAINT: Fix new cache name
* :pr:`6117`: BUG: Remove extra LICENSE.txt and setup.cfg
* :pr:`6105`: Update correlation_tools.py
* :pr:`6050`: BUG: MLEModel now passes nobs to Representation
* :pr:`6205`: MAINT: Exclude pytest-xdist 1.30
* :pr:`6246`: TST: Add Python 3.8 environment
:orphan:

=============
Release 0.8.0
=============

See also changes in the unreleased 0.7

Release summary
---------------

The main features of this release are several new time series models based
on the statespace framework, multiple imputation using MICE as well as many
other enhancements. The codebase also has been updated to be compatible with
recent numpy and pandas releases.

statsmodels is using now github to store the updated documentation which
is available under
https://www.statsmodels.org/stable for the last release, and
https://www.statsmodels.org/devel/ for the development version.

This is the last release that supports Python 2.6.


**Warning**

API stability is not guaranteed for new features, although even in this case
changes will be made in a backwards compatible way if possible. The stability
of a new feature depends on how much time it was already in statsmodels main
and how much usage it has already seen.
If there are specific known problems or limitations, then they are mentioned
in the docstrings.


The following major new features appear in this version.

Statespace Models
-----------------

Building on the statespace framework and models added in 0.7, this release
includes additional models that build on it.
Authored by Chad Fulton largely during GSOC 2015

Kalman Smoother
^^^^^^^^^^^^^^^

The Kalman smoother (introduced in #2434) allows making inference on the
unobserved state vector at each point in time using data from the entire
sample. In addition to this improved inference, the Kalman smoother is required
for future improvements such as simulation smoothing and the expectation
maximization (EM) algorithm.

As a result of this improvement, all state space models now inherit a `smooth`
method for producing results with smoothed state estimates. In addition, the
`fit` method will return results with smoothed estimates at the maximum
likelihood estimates.

Postestimation
^^^^^^^^^^^^^^

Improved post-estimation output is now available to all state space models
(introduced in #2566). This includes the new methods `get_prediction` and
`get_forecast`, providing standard errors and confidence intervals as well
as point estimates, `simulate`, providing simulation of time series following
the given state space process, and `impulse_responses`, allowing computation
of impulse responses due to innovations to the state vector.

Diagnostics
^^^^^^^^^^^

A number of general diagnostic tests on the residuals from state space
estimation are now available to all state space models (introduced in #2431).
These include:

* `test_normality` implements the Jarque-Bera test for normality of residuals
* `test_heteroskedasticity` implements a test for homoskedasticity of
  residuals similar to the Goldfeld-Quandt test
* `test_serial_correlation` implements the Ljung-Box (or Box-Pierce) test for
  serial correlation of residuals

These test statistics are also now included in the `summary` method output. In
addition, a `plot_diagnostics` method is available which provides four plots
to visually assess model fit.

Unobserved Components
^^^^^^^^^^^^^^^^^^^^^

The class of univariate Unobserved Components models (also known as structural
time series models) are now available (introduced in #2432). This includes as
special cases the local level model and local linear trend model. Generically
it allows decomposing a time series into trend, cycle, seasonal, and
irregular components, optionally with exogenous regressors and / or
autoregressive errors.

Multivariate Models
^^^^^^^^^^^^^^^^^^^

Two standard multivariate econometric models - vector autoregressive
moving-average model with exogenous regressors (VARMAX) and Dynamic Factors
models - are now available (introduced in #2563). The first is a popular
reduced form method of exploring the covariance in several time series, and the
second is a popular reduced form method of extracting a small number of common
factors from a large dataset of observed series.

Recursive least squares
^^^^^^^^^^^^^^^^^^^^^^^

A model for recursive least squares, also known as expanding-window OLS, is
now available in `statsmodels.regression` (introduced in #2830).

Miscellaneous
^^^^^^^^^^^^^

Other improvements to the state space framework include:

* Improved missing data handling #2770, #2809
* Ongoing refactoring and bug fixes in fringes and corner cases


Time Series Analysis
--------------------

Markov Switching Models
^^^^^^^^^^^^^^^^^^^^^^^

Markov switching dynamic regression and autoregression models are now
available (introduced in #2980 by Chad Fulton). These models allow regression
effects and / or autoregressive dynamics to differ depending on an unobserved
"regime"; in Markov switching models, the regimes are assumed to transition
according to a Markov process.

Statistics
^^^^^^^^^^

* KPSS stationarity, unit root test #2775 (N-Wouda)
* The Brock Dechert Scheinkman (BDS) test for nonlinear dependence is now
  available (introduced in #934 by Chad Fulton)
* Augmented Engle/Granger cointegration test (refactor hidden function) #3146 (Josef Perktold)


New functionality in statistics
-------------------------------

Contingency Tables #2418 (Kerby Shedden)

Local FDR, multiple testing #2297 (Kerby Shedden)

Mediation Analysis #2352 (Kerby Shedden)

Confidence intervals for multinomial proportions #3162 (Sebastien Lerique, Josef Perktold)

other:

* weighted quantiles in DescrStatsW #2707 (Kerby Shedden)


Duration
--------

Kaplan Meier Survival Function #2614 (Kerby Shedden)

Cumulative incidence rate function #3016 (Kerby Shedden)

other:

* frequency weights in Kaplan-Meier #2992 (Kerby Shedden)
* entry times for Kaplan-Meier #3126 (Kerby Shedden)
* intercept handling for PHReg #3095 (Kerby Shedden)


Imputation
----------

new subpackage in `statsmodels.imputation`

MICE #2076  (Frank Cheng GSOC 2014 and Kerby Shedden)

Imputation by regression on Order Statistic  #3019 (Paul Hobson)


Penalized Estimation
--------------------

Elastic net: fit_regularized with L1/L2 penalization has been added to
OLS, GLM and PHReg (Kerby Shedden)


GLM
---

Tweedie is now available as new family #2872 (Peter Quackenbush, Josef Perktold)

other:

* frequency weights for GLM (currently without full support) #
* more flexible convergence options #2803 (Peter Quackenbush)


Multivariate
------------

new subpackage that currently contains PCA

PCA was added in 0.7 to statsmodels.tools and is now in statsmodels.multivariate


Documentation
-------------

New doc build with latest jupyter and Python 3 compatibility (Tom Augspurger)


Other important improvements
----------------------------

several existing functions have received improvements


* seasonal_decompose: improved periodicity handling #2987 (ssktotoro ?)
* tools add_constant, add_trend: refactoring and pandas compatibility #2240 (Kevin Sheppard)
* acf, pacf, acovf: option for missing handling #3020 (joesnacks ?)
* acf, pacf plots: allow array of lags #2989 (Kevin Sheppard)
* pickling support for ARIMA #3412 (zaemyung)
* io SimpleTable (summary): allow names with special characters #3015 (tvanessa ?)
* tsa tools lagmat, lagmat2ds: pandas support #2310 #3042 (Kevin Sheppard)
* CompareMeans: from_data, summary methods #2754 (Valery Tyumen)
* API cleanup for robust, sandwich covariances #3162 (Josef Perktold)
* influence plot used swapped arguments (bug) #3158



Major Bugs fixed
----------------

* see github issues

While most bugs are usability problems, there is now a new label `type-bug-wrong`
for bugs that cause that silently incorrect numbers are returned.
https://github.com/statsmodels/statsmodels/issues?q=label%3Atype-bug-wrong+is%3Aclosed



Backwards incompatible changes and deprecations
-----------------------------------------------

* predict now returns a pandas Series if the exog argument is a DataFrame,
  including missing/NaN values
* PCA moved to multivariate compared to 0.7


Development summary and credits
-------------------------------

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance came from

* Kevin Sheppard
* Pierre Barbier de Reuille
* Tom Augsburger

and the general maintainer and code reviewer

* Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.8 release (based on git log):

.. note::

   * Ashish
   * Brendan
   * Brendan Condon
   * BrianLondon
   * Chad Fulton
   * Chris Fonnesbeck
   * Christian Lorentzen
   * Christoph T. Weidemann
   * James Kerns
   * Josef Perktold
   * Kerby Shedden
   * Kevin Sheppard
   * Leoyzen
   * Matthew Brett
   * Niels Wouda
   * Paul Hobson
   * Pierre Barbier de Reuille
   * Pietro Battiston
   * Ralf Gommers
   * Roman Ring
   * Skipper Seabold
   * Soren Fuglede Jorgensen
   * Thomas Cokelaer
   * Tom Augspurger
   * ValeryTyumen
   * Vanessa
   * Yaroslav Halchenko
   * dhpiyush
   * joesnacks
   * kokes
   * matiumerca
   * rlan
   * ssktotoro
   * thequackdaddy
   * vegcev

Thanks to all of the contributors for the 0.7 release:

.. note::

   * Alex Griffing
   * Antony Lee
   * Chad Fulton
   * Christoph Deil
   * Daniel Sullivan
   * Hans-Martin von Gaudecker
   * Jan Schulz
   * Joey Stockermans
   * Josef Perktold
   * Kerby Shedden
   * Kevin Sheppard
   * Kiyoto Tamura
   * Louis-Philippe Lemieux Perreault
   * Padarn Wilson
   * Ralf Gommers
   * Saket Choudhary
   * Skipper Seabold
   * Tom Augspurger
   * Trent Hauck
   * Vincent Arel-Bundock
   * chebee7i
   * donbeo
   * gliptak
   * hlin117
   * jerry dumblauskas
   * jonahwilliams
   * kiyoto
   * neilsummers
   * waynenilsen

These lists of names are automatically generated based on git log, and may not be
complete.
:orphan:

==============
Release 0.12.0
==============

Release summary
===============

statsmodels is using github to store the updated documentation. Two version are available:

- `Stable <https://www.statsmodels.org/>`_, the latest release
- `Development <https://www.statsmodels.org/devel/>`_, the latest build of the main branch

**Warning**

API stability is not guaranteed for new features, although even in
this case changes will be made in a backwards compatible way if
possible. The stability of a new feature depends on how much time it
was already in statsmodels main and how much usage it has already
seen.  If there are specific known problems or limitations, then they
are mentioned in the docstrings.

Stats
-----
**Issues Closed**: 239

**Pull Requests Merged**: 221

The Highlights
==============

Statistics
----------

New functions for hypothesis tests return a `HolderTuple` instance which
allows tuple indexing and unpacking for ``(statistic, pvalue)``, but also has
attribute access for those and for additional results statistics.

Meta-Analysis
~~~~~~~~~~~~~

Functions for Meta-Analysis have been added in :mod:`~statsmodels.stats.meta_analysis`.
The function :func:`~statsmodels.stats.meta_analysis.combine_effects` performs
fixed effects and random effects analysis. Several methods such as Paule-Mandel
and DerSimonian-Laird are available to estimate the random effects variance.
The module also includes effect size functions for standardized mean difference
and for proportions that can be used with :func:`~statsmodels.stats.meta_analysis.combine_effects`.
A notebook illustrates the usage of the new features for meta-analysis.

New hypothesis test for 2 samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hypothesis tests, confidence intervals and power functions have been added
for proportions from two independent samples. Inferential statistics are
available for difference, ratio and odds-ratio of the two proportions.
Equivalence testing for two independent proportions is available based on
two one-sided tests TOST.

Hypothesis tests including equivalence test, for the ratio of two
independent Poisson rates are now available in
:func:`~statsmodels.stats.rates.test_poisson_2indep` and
:func:`~statsmodels.stats.rates.tost_poisson_2indep`

Oneway ANOVA-type analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

Several statistical methods for ANOVA-type analysis of k independent samples
have been added in module :mod:`~statsmodels.stats.oneway`. This includes
standard Anova, Anova for unequal variances (Welch, Brown-Forsythe for mean),
Anova based on trimmed samples (Yuen anova) and equivalence testing using
the method of Wellek.
Anova for equality of variances or dispersion are available for several
transformations. This includes Levene test and Browne-Forsythe test for equal
variances as special cases. It uses the `anova_oneway` function, so unequal
variance and trimming options are also available for tests on variances.
Several functions for effect size measures have been added, that can be used
for reporting or for power and sample size computation.

Multivariate statistics
~~~~~~~~~~~~~~~~~~~~~~~

The new module :mod:`~statsmodels.stats.multivariate` includes one and
two sample tests for multivariate means, Hotelling's t-tests',
:func:`~statsmodels.stats.multivariate.test_mvmean`,
:func:`~statsmodels.stats.multivariate.test_mvmean_2indep` and confidence
intervals for one-sample multivariate mean
:func:`~statsmodels.stats.multivariate.confint_mvmean`
Additionally, hypothesis tests for covariance patterns, and for oneway equality
of covariances are now available in several ``test_cov`` functions.


Time-Series Analysis
--------------------

New exponential smoothing model: ETS (Error, Trend, Seasonal)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Class implementing ETS models :class:`~statsmodels.tsa.exponential_smoothing.ets.ETSModel`.
- Includes linear and non-linear exponential smoothing models
- Supports parameter fitting, in-sample prediction and out-of-sample
  forecasting, prediction intervals, simulation, and more.
- Based on the innovations state space approach.

Statespace Models
-----------------

New dynamic factor model for large datasets and monthly / quarterly mixed frequency models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- New dynamic factor model :class:`~statsmodels.tsa.statespace.dynamic_factor_mq.DynamicFactorMQ`.
- Allows for hundreds of observed variables, by fitting with the EM algorithm
- Allows specifying factors that load only on a specific group of variables
- Allows for monthly / quarterly mixed frequency models. For example, this
  supports one popular approach to "Nowcasting" GDP

Decomposition of forecast updates based on the "news"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- New :meth:`~statsmodels.tsa.statespace.mlemodel.MLEResults.news` method for state space model results objects
- Links updated data to changes in forecasts
- Supports "nowcasting" exercises that progressively incorporate more and more
  information as time goes on

Sparse Cholesky Simulation Smoother
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- New option for simulation smoothing in state space models via the
  "Cholesky factor algorithm" (CFA) approach in
  :class:`~statsmodels.tsa.statespace.cfa_simulation_smoother.CFASimulationSmoother`
- Takes advantage of algorithms for sparse Cholesky factorization, rather than
  using the typical simulation smoother based on Kalman filtering and smoothing

Option to use Chadrasekhar recursions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- New option for state space models to use Chandrasekhar recursions rather than
  than the typical Kalman filtering recursions by setting ``filter_chandrasekhar=True``.
- Improved performance for some models with large state vectors

Forecasting Methods
~~~~~~~~~~~~~~~~~~~
Two popular method for forecasting time series, forecasting after STL decomposition
(:class:`~statsmodels.tsa.forecasting.stl.STLForecast`)
and the Theta model (:class:`~statsmodels.tsa.forecasting.theta.ThetaModel`) have
been added.

Complex Deterministic Terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`~statsmodels.tsa.deterministic.DeterministicProcess` can be used to generate
deterministic processes containing time trends, seasonal dummies and Fourier components.
A :class:`~statsmodels.tsa.deterministic.DeterministicProcess` can be used to produce
in-sample regressors or out-of-sample values suitable for forecasting.


What's new - an overview
========================

The following lists the main new features of statsmodels 0.12.0. In addition,
release 0.12.0 includes bug fixes, refactorings and improvements in many areas.

Submodules
----------


``Documentation``
~~~~~~~~~~~~~~~~~
- Fix the version that appears in the documentation  (:pr:`6452`)
- Send log to dev/null/  (:pr:`6456`)
- Correct spelling of various  (:pr:`6518`)
- Fix typos  (:pr:`6531`)
- Update interactions_anova.ipynb  (:pr:`6601`)
- Fix `true` type on statespace docs page  (:pr:`6616`)
- Minor fixes for holtwinters simulate  (:pr:`6631`)
- Change OLS example to use datasets  (:pr:`6656`)
- Fix AutoReg docstring  (:pr:`6662`)
- Fix `fdrcorrection` docstring missing `is_sorted` parameter  (:pr:`6680`)
- Add new badges  (:pr:`6704`)
- Fix number if notebook text  (:pr:`6709`)
- Improve Factor and related docstrings  (:pr:`6719`)
- Improve explantion of missing values in ACF and related  (:pr:`6726`)
- Notebook for quasibinomial regression  (:pr:`6732`)
- Improve "conservative" doc  (:pr:`6738`)
- Update broken link  (:pr:`6742`)
- Fix broken links with 404 error  (:pr:`6746`)
- Demonstrate variance components analysis  (:pr:`6758`)
- Make deprecations more visible  (:pr:`6775`)
- Numpydoc signatures  (:pr:`6825`)
- Correct reference in docs  (:pr:`6837`)
- Include dot_plot  (:pr:`6841`)
- Updated durbin_watson Docstring and Tests  (:pr:`6848`)
- Explain low df in cluster  (:pr:`6853`)
- Fix common doc errors  (:pr:`6862`)
- Small doc fixes  (:pr:`6874`)
- Fix issues in docs related to exponential smoothing  (:pr:`6879`)
- Spelling and other doc fixes  (:pr:`6902`)
- Correct spacing around colon in docstrings  (:pr:`6903`)
- Initial 0.12 Release Note  (:pr:`6923`)
- Fix doc errors and silence warning  (:pr:`6931`)
- Clarify deprecations  (:pr:`6932`)
- Document exceptions and warnings  (:pr:`6943`)
- Update pandas function in hp_filter example  (:pr:`6946`)
- Prepare docs  (:pr:`6948`)
- Fix final issues in release note  (:pr:`6951`)
- Final doc fixed for 0.12.0rc0  (:pr:`6965`)
- Update DeterministicProcess docs  (:pr:`6968`)
- Add docstring to string_like method  (:pr:`6972`)
- Fix LaTeX in seasonal notebook  (:pr:`6976`)
- Add new stats to release notes for 0.12  (:pr:`7001`)

``Performance``
~~~~~~~~~~~~~~~
- State space: add Chandrasekhar recursions  (:pr:`6411`)
- Speed up HC2/HC3 standard error calculation, using less memory  (:pr:`6664`)
- Sparse matrices in MixedLM  (:pr:`6766`)

``backport``
~~~~~~~~~~~~
- `MLEResults.states.predicted` has wrong index  (:pr:`6580`)
- State space: simulate with time-varying covariance matrices.  (:pr:`6607`)
- State space: error with collapsed observations when missing  (:pr:`6613`)
- Dataframe/series concatenation in statespace results append  (:pr:`6768`)
- Pass cov_type, cov_kwargs through ARIMA.fit  (:pr:`6770`)

``base``
~~~~~~~~
- Don't attach patsy constraint instance   (:pr:`6521`)
- Fix constraints and bunds when use scipy.optimize.minimize  (:pr:`6657`)
- Correct shape of fvalue and f_pvalue  (:pr:`6831`)
- Correct dimension when data removed  (:pr:`6888`)

``build``
~~~~~~~~~
- Use pip on Azure  (:pr:`6474`)
- Attempt to cache key docbuild files  (:pr:`6490`)
- Improve doc caching  (:pr:`6491`)
- Azure: Mac OSX 10.13 -> 10.14  (:pr:`6587`)

``discrete``
~~~~~~~~~~~~
- Don't attach patsy constraint instance   (:pr:`6521`)
- Sparse matrices in MixedLM  (:pr:`6766`)
- Catch warnings in discrete  (:pr:`6836`)
- Add improved .cdf() and .ppf() to discrete distributions  (:pr:`6938`)
- Remove k_extra from effects_idx  (:pr:`6939`)
- Improve count model tests  (:pr:`6940`)

``docs``
~~~~~~~~
- Fix doc errors and silence warning  (:pr:`6931`)
- Prepare docs  (:pr:`6948`)

``duration``
~~~~~~~~~~~~
- Allow more than 2 groups for survdiff in statsmodels.duration  (:pr:`6626`)

``gam``
~~~~~~~
- Fix GAM for 1-dim exog_linear   (:pr:`6520`)
- Fixed BSplines to match existing docs  (:pr:`6915`)

``genmod``
~~~~~~~~~~
- Change default optimizer for glm/ridge and make it user-settable  (:pr:`6438`)
- Fix exposure/offset handling in GEEResults  (:pr:`6475`)
- Use GLM starting values for QIF  (:pr:`6514`)
- Don't attach patsy constraint instance   (:pr:`6521`)
- Allow GEE weights to vary within clusters  (:pr:`6582`)
- Calculate AR covariance parameters for gridded data  (:pr:`6621`)
- Warn for non-convergence in elastic net  (:pr:`6697`)
- Gh 6627  (:pr:`6852`)
- Change of BIC formula in GLM  (:pr:`6941`)
- Make glm's predict function return numpy array even if exposure is a pandas series  (:pr:`6942`)
- Fix check for offset_exposure in null  (:pr:`6957`)
- Add test for offset exposure null  (:pr:`6959`)

``graphics``
~~~~~~~~~~~~
- Include figsize as parameter for IRF plot  (:pr:`6590`)
- Speed up banddepth calculations  (:pr:`6744`)
- Fix logic in labeling corr plot  (:pr:`6818`)
- Enable qqplot_2sample to handle uneven samples  (:pr:`6906`)
- Support frozen dist in ProbPlots  (:pr:`6910`)

``io``
~~~~~~
- Handle pathlib.Path objects  (:pr:`6654`)
- Added label option to summary.to_latex()  (:pr:`6895`)
- Fixed the shifted column names in summary.to_latex()  (:pr:`6900`)
- Removed additional hline between tabulars  (:pr:`6905`)

``maintenance``
~~~~~~~~~~~~~~~
- Special docbuild  (:pr:`6457`)
- Special docbuild"  (:pr:`6460`)
- Correcting typo  (:pr:`6461`)
- Avoid noise in f-pvalue  (:pr:`6465`)
- Replace Python 3.5 with 3.8 on Azure  (:pr:`6466`)
- Update supported versions  (:pr:`6467`)
- Fix future warnings  (:pr:`6469`)
- Fix issue with ragged array  (:pr:`6471`)
- Avoid future error  (:pr:`6473`)
- Silence expected visible deprecation warning  (:pr:`6477`)
- Remove Python 3.5 references  (:pr:`6492`)
- Avoid calling depr code  (:pr:`6493`)
- Use travis cache and optimize build times  (:pr:`6495`)
- Relax tolerance on test that occasionally fails  (:pr:`6534`)
- Relax tolerance on test that randomly fails  (:pr:`6588`)
- Fix appveyor/conda  (:pr:`6653`)
- Delete empty directory  (:pr:`6671`)
- Flake8 fixes  (:pr:`6710`)
- Remove deprecated keyword  (:pr:`6712`)
- Remove OrderedDict  (:pr:`6715`)
- Remove dtype np.integer for avoid Dep Warning  (:pr:`6728`)
- Update pip-pre links  (:pr:`6733`)
- Spelling and small fixes  (:pr:`6752`)
- Remove error on FutureWarning  (:pr:`6811`)
- Fix failing tests  (:pr:`6817`)
- Replace Warnings with Notes in regression summary  (:pr:`6828`)
- Numpydoc should work now  (:pr:`6842`)
- Deprecate categorical  (:pr:`6843`)
- Remove redundant definition  (:pr:`6845`)
- Relax tolerance on test that fails Win32  (:pr:`6849`)
- Fix error on nightly build  (:pr:`6850`)
- Correct debugging info  (:pr:`6855`)
- Mark VAR from_formula as NotImplemented  (:pr:`6865`)
- Allow skip if rdataset fails  (:pr:`6871`)
- Improve lint  (:pr:`6885`)
- Change default lag in serial correlation tests  (:pr:`6893`)
- Ensure setuptools is imported first  (:pr:`6894`)
- Remove FutureWarnings  (:pr:`6920`)
- Add tool to simplify documenting API in release notes  (:pr:`6922`)
- Relax test tolerance for future compat  (:pr:`6945`)
- Fixes for failures in wheel building  (:pr:`6952`)
- Fixes for wheel building  (:pr:`6954`)
- Remove print statements  (:pr:`6985`)
- Update Azure images  (:pr:`6992`)

``multivariate``
~~~~~~~~~~~~~~~~
- Multivariate mean tests and confint  (:pr:`4107`)
- Improve missing value handling in PCA  (:pr:`6705`)

``nonparametric``
~~~~~~~~~~~~~~~~~
- Fix #6511  (:pr:`6515`)
- Fix domain check  (:pr:`6547`)
- Ensure sigma estimate is positive in KDE  (:pr:`6713`)
- Fix access to normal_reference_constant  (:pr:`6806`)
- Add xvals param to lowess smoother  (:pr:`6908`)
- Return self from KDEUnivariate fit  (:pr:`6991`)
- Allow custom bandwidth functions in KDEUnivariate fit  (:pr:`7002`)

``regression``
~~~~~~~~~~~~~~
- Statsmodels.regression.linear_model.OLS.fit_regularized fails to generate correct answer (#6604)  (:pr:`6608`)
- Change OLS example to use datasets  (:pr:`6656`)
- Speed up HC2/HC3 standard error calculation, using less memory  (:pr:`6664`)
- Fix summary col R2 ordering  (:pr:`6714`)
- Insufficient input checks in QuantReg  (:pr:`6747`)
- Add expanding initialization to RollingOLS/WLS  (:pr:`6838`)
- Add  a note when R2 is uncentered  (:pr:`6844`)

``robust``
~~~~~~~~~~
- Add normalized iqr to robust.scales  (:pr:`6969`)
- Robust.scale.iqr does need centering, since quantiles are translation equivariant  (:pr:`6973`)
- Add robust qn scale  (:pr:`6990`)
- Fix bug where mad ignores center if center is not callable  (:pr:`7000`)

``stats``
~~~~~~~~~
- Multivariate mean tests and confint  (:pr:`4107`)
- Fix tukey-hsd for 1 pvalue   (:pr:`6470`)
- Add option for original Breusch-Pagan heteroscedasticity test  (:pr:`6508`)
- ENH Allow optional regularization in local fdr  (:pr:`6622`)
- Add meta-analysis (basic methods)  (:pr:`6632`)
- Add two independent proportion inference rebased  (:pr:`6675`)
- Rates, poisson means two-sample comparison  rebased  (:pr:`6677`)
- Stats.base, add HolderTuple, Holder class with indexing  (:pr:`6678`)
- Add covariance structure hypothesis tests  (:pr:`6693`)
- Raise exception when recursive residual is not well defined  (:pr:`6727`)
- Mediation support for PH regression  (:pr:`6782`)
- Stats robust rebased2  (:pr:`6789`)
- Hotelling's Two Sample Mean Test  (:pr:`6810`)
- Stats moment_helpers use random state in unit test  (:pr:`6835`)
- Updated durbin_watson Docstring and Tests  (:pr:`6848`)
- Add recent stats addition to docs  (:pr:`6859`)
- REF/DOC docs and refactor of recent stats  (:pr:`6872`)
- Api cleanup and improve docstrings in stats, round 3  (:pr:`6897`)
- Improve descriptivestats  (:pr:`6944`)
- Catch warning  (:pr:`6964`)

``tools``
~~~~~~~~~
- Return column information in add_constant  (:pr:`6830`)
- Add QR-based matrix rank  (:pr:`6834`)
- Add Root Mean Square Percentage Error  (:pr:`6926`)

``tsa``
~~~~~~~
- Fixes #6553, sliced predicted values according to predicted index  (:pr:`6556`)
- Holt-Winters simulations  (:pr:`6560`)
- Example notebook (r): stationarity and detrending (ADF/KPSS)  (:pr:`6614`)
- Ensure text comparison is lower  (:pr:`6628`)
- Minor fixes for holtwinters simulate  (:pr:`6631`)
- New exponential smoothing implementation  (:pr:`6699`)
- Improve warning message in KPSS  (:pr:`6711`)
- Change trend initialization in STL  (:pr:`6722`)
- Add check in test_whiteness  (:pr:`6723`)
- Raise on incorrectly sized exog  (:pr:`6730`)
- Add deterministic processes  (:pr:`6751`)
- Add Theta forecasting method  (:pr:`6767`)
- Automatic lag selection for Box-Pierce, Ljung-Box #6645  (:pr:`6785`)
- Fix missing str  (:pr:`6827`)
- Add support for PeriodIndex to AutoReg  (:pr:`6829`)
- Error in append for ARIMA model with trend  (:pr:`6832`)
- Add QR-based matrix rank  (:pr:`6834`)
- Rename unbiased to adjusted  (:pr:`6839`)
- Ensure PACF lag length is sensible  (:pr:`6846`)
- Allow Series as exog in predict  (:pr:`6847`)
- Raise on nonstationary parameters when attempting to use GLS  (:pr:`6854`)
- Relax test tolerance  (:pr:`6856`)
- Limit maxlags in VAR  (:pr:`6867`)
- Fix indexing with HoltWinters's forecast  (:pr:`6869`)
- Refactor Holt-Winters  (:pr:`6870`)
- Fix raise exception on granger causality test  (:pr:`6877`)
- Get_prediction method for ETS  (:pr:`6882`)
- Ets: test for simple exponential smoothing convergence  (:pr:`6884`)
- Added diagnostics test to ETS model  (:pr:`6892`)
- Stop transforming ES components  (:pr:`6904`)
- Fix extend in VARMAX with trend  (:pr:`6909`)
- Add STL Forecasting method  (:pr:`6911`)
- Dynamic is incorrect when not an int in statespace get_prediction  (:pr:`6917`)
- Correct IRF nobs with exog  (:pr:`6925`)
- Add get_prediction to AutoReg  (:pr:`6927`)
- Standardize forecast API  (:pr:`6933`)
- Fix small issues post ETS get_prediction merge  (:pr:`6934`)
- Modify failing test on Windows  (:pr:`6949`)
- Improve ETS / statespace documentation and highlights for v0.12   (:pr:`6950`)
- Remove FutureWarnings  (:pr:`6958`)

``tsa.statespace``
~~~~~~~~~~~~~~~~~~
- State space: add Chandrasekhar recursions  (:pr:`6411`)
- Use reset_randomstate  (:pr:`6433`)
- State space: add "Cholesky factor algorithm" simulation smoothing  (:pr:`6501`)
- Bayesian estimation of SARIMAX using PyMC3 NUTS  (:pr:`6528`)
- State space: compute smoothed state autocovariance matrices for arbitrary lags  (:pr:`6579`)
- `MLEResults.states.predicted` has wrong index  (:pr:`6580`)
- State space: simulate with time-varying covariance matrices.  (:pr:`6607`)
- State space: error with collapsed observations when missing  (:pr:`6613`)
- Notebook describing how to create state space custom models  (:pr:`6682`)
- Fix covariance estimation in parameterless models  (:pr:`6688`)
- Fix state space linting errors.  (:pr:`6698`)
- Decomposition of forecast updates in state space models due to the "news"  (:pr:`6765`)
- Dataframe/series concatenation in statespace results append  (:pr:`6768`)
- Pass cov_type, cov_kwargs through ARIMA.fit  (:pr:`6770`)
- Improve univariate smoother performance  (:pr:`6797`)
- Add `news` example notebook image.  (:pr:`6800`)
- Fix extend in VARMAX with trend  (:pr:`6909`)
- Dynamic is incorrect when not an int in statespace get_prediction  (:pr:`6917`)
- Add dynamic factor model with EM algorithm, option for monthly/quarterly mixed frequency model  (:pr:`6937`)
- Improve ETS / statespace documentation and highlights for v0.12   (:pr:`6950`)
- SARIMAX throwing different errors when length of endogenous var is too low  (:pr:`6961`)
- Fix start params computation with few nobs  (:pr:`6962`)
- Relax tolerance on random failure  (:pr:`6963`)

``tsa.vector.ar``
~~~~~~~~~~~~~~~~~
- Include figsize as parameter for IRF plot  (:pr:`6590`)
- Raise on incorrectly sized exog  (:pr:`6730`)
- Correct IRF nobs with exog  (:pr:`6925`)

bug-wrong
---------

A new issue label `type-bug-wrong` indicates bugs that cause that incorrect
numbers are returned without warnings.
(Regular bugs are mostly usability bugs or bugs that raise an exception for
unsupported use cases.)
`see tagged issues <https://github.com/statsmodels/statsmodels/issues?q=is%3Aissue+label%3Atype-bug-wrong+is%3Aclosed+milestone%3A0.12/>`_


Major Bugs Fixed
================

See github issues for a list of bug fixes included in this release

- `Closed bugs <https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.12+label%3Atype-bug/>`_
- `Closed bugs (wrong result) <https://github.com/statsmodels/statsmodels/pulls?q=is%3Apr+is%3Amerged+milestone%3A0.12+label%3Atype-bug-wrong/>`_


Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

- Chad Fulton
- Brock Mendel
- Peter Quackenbush
- Kerby Shedden
- Kevin Sheppard

and the general maintainer and code reviewer

- Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.12.0 release (based on git log):

- Alex Lyttle
- Amund Vedal
- Baran Karakus
- Batakrishna Sahu
- Chad Fulton
- Cinthia M. Tanaka
- Dorian Bivolaru
- Ezequiel Smucler
- Giulio Beseghi
- Haoyu Qi
- Hassan Kibirige
- He Yang
- Henning Blunck
- Jimmy2027
- Joon Ro
- Joonsuk Park
- Josef Perktold
- Kerby Shedden
- Kevin Rose
- Kevin Sheppard
- Manmeet Kumar Chaudhuri
- Markus Löning
- Martin Larralde
- Nolan Conaway
- Paulo Galuzio
- Peter Prescott
- Peter Quackenbush
- Samuel Scherrer
- Sean Lane
- Sebastian Pölsterl
- Skipper Seabold
- Thomas Brooks
- Thomas Marchand
- Tim Gates
- Victor Ananyev
- Wouter De Coster
- Zhiqing Xiao
- adrienpacifico
- aeturrell
- cd
- das-soham
- eirki
- pag
- partev
- tagoma
- w31ha0


These lists of names are automatically generated based on git log, and may not
be complete.

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

- :pr:`4107`: ENH: multivariate mean tests and confint
- :pr:`6411`: ENH: state space: add Chandrasekhar recursions
- :pr:`6433`: TST/BUG: use reset_randomstate
- :pr:`6438`: BUG: Change default optimizer for glm/ridge and make it user-settable
- :pr:`6452`: DOC: Fix the version that appears in the documentation
- :pr:`6456`: DOC: Send log to dev/null/
- :pr:`6457`: DOC: Special docbuild
- :pr:`6460`: Revert "DOC: Special docbuild"
- :pr:`6461`: MAINT: correcting typo
- :pr:`6465`: MAINT: Avoid noise in f-pvalue
- :pr:`6466`: MAINT: Replace Python 3.5 with 3.8 on Azure
- :pr:`6467`: MAINT: Update supported versions
- :pr:`6469`: MAINT: Fix future warnings
- :pr:`6470`: BUG: fix tukey-hsd for 1 pvalue 
- :pr:`6471`: MAINT: Fix issue with ragged array
- :pr:`6473`: MAINT: Avoid future error
- :pr:`6474`: BLD: Use pip on Azure
- :pr:`6475`: BUG: Fix exposure/offset handling in GEEResults
- :pr:`6477`: BUG: Silence expected visible deprecation warning
- :pr:`6490`: BLD: Attempt to cache key docbuild files
- :pr:`6491`: BLD: Improve doc caching
- :pr:`6492`: MAINT: Remove Python 3.5 references
- :pr:`6493`: MAINT: Avoid calling depr code
- :pr:`6495`: MAINT: Use travis cache and optimize build times
- :pr:`6501`: ENH: state space: add "Cholesky factor algorithm" simulation smoothing
- :pr:`6508`: ENH: Add option for original Breusch-Pagan heteroscedasticity test
- :pr:`6514`: ENH: use GLM starting values for QIF
- :pr:`6515`: BUG: fix #6511
- :pr:`6518`: DOC: Fix simple typo: various
- :pr:`6520`: BUG: fix GAM for 1-dim exog_linear
- :pr:`6521`: REF/BUG: don't attach patsy constraint instance
- :pr:`6528`: DOC: Bayesian estimation of SARIMAX using PyMC3 NUTS
- :pr:`6531`: DOC: fix typos
- :pr:`6534`: MAINT: Relax tolerance on test that occasionally fails
- :pr:`6547`: BUG: Fix domain check
- :pr:`6556`: BUG: fixes #6553, sliced predicted values according to predicted index
- :pr:`6560`: ENH: Holt-Winters simulations
- :pr:`6579`: ENH: state space: compute smoothed state autocovariance matrices for arbitrary lags
- :pr:`6580`: BUG: `MLEResults.states.predicted` has wrong index
- :pr:`6582`: ENH: Allow GEE weights to vary within clusters
- :pr:`6587`: BLD: Azure: Mac OSX 10.13 -> 10.14
- :pr:`6588`: MAINT: Relax tolerance on test that randomly fails
- :pr:`6590`: ENH: Include figsize as parameter for IRF plot
- :pr:`6601`: DOC: Update interactions_anova.ipynb
- :pr:`6607`: BUG: state space: simulate with time-varying covariance matrices.
- :pr:`6608`: BUG: statsmodels.regression.linear_model.OLS.fit_regularized fails to generate correct answer (#6604)
- :pr:`6613`: BUG: state space: error with collapsed observations when missing
- :pr:`6614`: DOC/ENH: example notebook (r): stationarity and detrending (ADF/KPSS)
- :pr:`6616`: DOC: Fix `true` type on statespace docs page
- :pr:`6621`: ENH: Calculate AR covariance parameters for gridded data
- :pr:`6622`: ENH Allow optional regularization in local fdr
- :pr:`6626`: ENH: allow more than 2 groups for survdiff in statsmodels.duration
- :pr:`6628`: BUG: Ensure text comparison is lower
- :pr:`6631`: DOC/TST: minor fixes for holtwinters simulate
- :pr:`6632`: ENH: add meta-analysis (basic methods)
- :pr:`6653`: MAINT: Fix appveyor/conda
- :pr:`6654`: ENH: Handle pathlib.Path objects
- :pr:`6656`: DOC: change OLS example to use datasets
- :pr:`6657`: BUG: fix constraints and bunds when use scipy.optimize.minimize
- :pr:`6662`: DOC: Fix AutoReg docstring
- :pr:`6664`: PERF: Speed up HC2/HC3 standard error calculation, using less memory
- :pr:`6671`: MAINT: Delete empty directory
- :pr:`6675`: ENH: add two independent proportion inference rebased
- :pr:`6677`: ENH: rates, poisson means two-sample comparison  rebased
- :pr:`6678`: ENH: stats.base, add HolderTuple, Holder class with indexing
- :pr:`6680`: DOC: Fix `fdrcorrection` docstring missing `is_sorted` parameter
- :pr:`6682`: ENH/DOC: Notebook describing how to create state space custom models
- :pr:`6688`: BUG: Fix covariance estimation in parameterless models
- :pr:`6693`: ENH: add covariance structure hypothesis tests
- :pr:`6697`: ENH: Warn for non-convergence in elastic net
- :pr:`6698`: CLN: Fix state space linting errors.
- :pr:`6699`: ENH: New exponential smoothing implementation
- :pr:`6704`: DOC: Add new badges
- :pr:`6705`: BUG\ENH: Improve missing value handling in PCA
- :pr:`6709`: DOC: Fix number if notebook text
- :pr:`6710`: MAINT: Flake8 fixes
- :pr:`6711`: ENH: Improve warning message in KPSS
- :pr:`6712`: MAINT: Remove deprecated keyword
- :pr:`6713`: BUG: Ensure sigma estimate is positive in KDE
- :pr:`6714`: BUG: Fix summary col R2 ordering
- :pr:`6715`: MAINT: Remove OrderedDict
- :pr:`6719`: DOC: Improve Factor and related docstrings
- :pr:`6722`: BUG: Change trend initialization in STL
- :pr:`6723`: ENH: Add check in test_whiteness
- :pr:`6726`: DOC: Improve explantion of missing values in ACF and related
- :pr:`6727`: ENH: Raise exception when recursive residual is not well defined
- :pr:`6728`: MAINT: Remove dtype np.integer for avoid Dep Warning
- :pr:`6730`: BUG: Raise on incorrectly sized exog
- :pr:`6732`: DOC: Notebook for quasibinomial regression
- :pr:`6733`: MAINT: Update pip-pre links
- :pr:`6738`: DOC: Improve "conservative" doc
- :pr:`6742`: Update broken link
- :pr:`6744`: ENH: Speed up banddepth calculations
- :pr:`6746`: DOC: Fix broken links with 404 error
- :pr:`6747`: BUG: Insufficient input checks in QuantReg
- :pr:`6751`: ENH: Add deterministic processes
- :pr:`6752`: MAINT: Spelling and small fixes
- :pr:`6758`: DOC: Demonstrate variance components analysis
- :pr:`6765`: ENH: Decomposition of forecast updates in state space models due to the "news"
- :pr:`6766`: PERF: Sparse matrices in MixedLM
- :pr:`6767`: ENH: Add Theta forecasting method
- :pr:`6768`: BUG: dataframe/series concatenation in statespace results append
- :pr:`6770`: BUG: pass cov_type, cov_kwargs through ARIMA.fit
- :pr:`6775`: DOC: Make deprecations more visible
- :pr:`6782`: ENH: Mediation support for PH regression
- :pr:`6785`: ENH: automatic lag selection for Box-Pierce, Ljung-Box #6645
- :pr:`6789`: ENH: Stats robust rebased2
- :pr:`6797`: ENH: improve univariate smoother performance
- :pr:`6800`: DOC: Add `news` example notebook image.
- :pr:`6806`: BUG: Fix access to normal_reference_constant
- :pr:`6810`: ENH: Hotelling's Two Sample Mean Test
- :pr:`6811`: MAINT: Remove error on FutureWarning
- :pr:`6817`: MAINT: Fix failing tests
- :pr:`6818`: BUG: Fix logic in labeling corr plot
- :pr:`6825`: DOC: Numpydoc signatures
- :pr:`6827`: BUG: Fix missing str
- :pr:`6828`: MAINT: Replace Warnings with Notes in regression summary
- :pr:`6829`: ENH: Add support for PeriodIndex to AutoReg
- :pr:`6830`: ENH: Return column information in add_constant
- :pr:`6831`: BUG: Correct shape of fvalue and f_pvalue
- :pr:`6832`: BUG: error in append for ARIMA model with trend
- :pr:`6834`: ENH: Add QR-based matrix rank
- :pr:`6835`: TST: stats moment_helpers use random state in unit test
- :pr:`6836`: MAINT: Catch warnings in discrete
- :pr:`6837`: DOC: Correct reference in docs
- :pr:`6838`: ENH: Add expanding initialization to RollingOLS/WLS
- :pr:`6839`: REF: Rename unbiased to adjusted
- :pr:`6841`: DOC: Include dot_plot
- :pr:`6842`: MAINT: numpydoc should work now
- :pr:`6843`: MAINT: Deprecate categorical
- :pr:`6844`: ENH: Add  a note when R2 is uncentered
- :pr:`6845`: MAINT: Remove redundant definition
- :pr:`6846`: BUG: Ensure PACF lag length is sensible
- :pr:`6847`: BUG: Allow Series as exog in predict
- :pr:`6848`: Updated durbin_watson Docstring and Tests
- :pr:`6849`: TST: Relax tolerance on test that fails Win32
- :pr:`6850`: MAINT: Fix error on nightly build
- :pr:`6852`: Gh 6627
- :pr:`6853`: DOC: Explain low df in cluster
- :pr:`6854`: BUG: Raise on nonstationary parameters when attempting to use GLS
- :pr:`6855`: MAINT: Correct debugging info
- :pr:`6856`: MAINT: Relax test tolerance
- :pr:`6859`: DOC: add recent stats addition to docs
- :pr:`6862`: DOC: Fix common doc errors
- :pr:`6865`: MAINT: Mark VAR from_formula as NotImplemented
- :pr:`6867`: BUG: Limit maxlags in VAR
- :pr:`6868`: TST: Refactor factor tests again
- :pr:`6869`: BUG: Fix indexing with HoltWinters's forecast
- :pr:`6870`: REF: Refactor Holt-Winters
- :pr:`6871`: MAINT: Allow skip if rdataset fails
- :pr:`6872`: REF/DOC docs and refactor of recent stats
- :pr:`6874`: DOC: Small doc fixes
- :pr:`6877`: BUG: fix raise exception on granger causality test
- :pr:`6879`: DOC: Fix issues in docs related to exponential smoothing
- :pr:`6882`: ENH: get_prediction method for ETS
- :pr:`6884`: TST: ets: test for simple exponential smoothing convergence
- :pr:`6885`: MAINT: Improve lint
- :pr:`6888`: BUG: Correct dimension when data removed
- :pr:`6892`: ENH: added diagnostics test to ETS model
- :pr:`6893`: MAINT: Change default lag in serial correlation tests
- :pr:`6894`: MAINT: Ensure setuptools is imported first
- :pr:`6895`: ENH: Added label option to summary.to_latex()
- :pr:`6897`: REF/DOC: api cleanup and improve docstrings in stats, round 3
- :pr:`6900`: ENH: Fixed the shifted column names in summary.to_latex()
- :pr:`6902`: DOC: Spelling and other doc fixes
- :pr:`6903`: DOC: Correct spacing around colon in docstrings
- :pr:`6904`: BUG: Stop transforming ES components
- :pr:`6905`: ENH: removed additional hline between tabulars
- :pr:`6906`: ENH: Enable qqplot_2sample to handle uneven samples
- :pr:`6908`: ENH: Add xvals param to lowess smoother
- :pr:`6909`: BUG: fix extend in VARMAX with trend
- :pr:`6910`: ENH: Support frozen dist in ProbPlots
- :pr:`6911`: ENH: Add STL Forecasting method
- :pr:`6915`: BUG: Fixed BSplines to match existing docs
- :pr:`6917`: BUG: dynamic is incorrect when not an int in statespace get_prediction
- :pr:`6920`: MAINT: Remove FutureWarnings
- :pr:`6922`: ENH: Add tool to simplify documenting API in release notes
- :pr:`6923`: DOC: Initial 0.12 Release Note
- :pr:`6925`: BUG: Correct IRF nobs with exog
- :pr:`6926`: ENH: Add Root Mean Square Percentage Error
- :pr:`6927`: ENH: Add get_prediction to AutoReg
- :pr:`6931`: DOC/MAINT: Fix doc errors and silence warning
- :pr:`6932`: DOC: Clarify deprecations
- :pr:`6933`: MAINT: Standardize forecast API
- :pr:`6934`: MAINT: Fix small issues post ETS get_prediction merge
- :pr:`6937`: ENH: Add dynamic factor model with EM algorithm, option for monthly/quarterly mixed frequency model
- :pr:`6938`: ENH: Add improved .cdf() and .ppf() to discrete distributions
- :pr:`6939`: BUG: remove k_extra from effects_idx
- :pr:`6940`: TST: Improve count model tests
- :pr:`6941`: REF: Change of BIC formula in GLM
- :pr:`6942`: BUG: Make glm's predict function return numpy array even if exposure is a pandas series
- :pr:`6943`: DOC: Document exceptions and warnings
- :pr:`6944`: ENH: improve descriptivestats
- :pr:`6945`: MAINT: Relax test tolerance for future compat
- :pr:`6946`: DOC: update pandas function in hp_filter example
- :pr:`6948`: Prepare docs
- :pr:`6949`: TST: Modify failing test on Windows
- :pr:`6950`: DOC: improve ETS / statespace documentation and highlights for v0.12
- :pr:`6951`: DOC: Fix final issues in release note
- :pr:`6952`: MAINT: Fixes for failures in wheel building
- :pr:`6954`: MAINT: Fixes for wheel building
- :pr:`6957`: BUG: Fix check for offset_exposure in null
- :pr:`6958`: MAINT: Remove FutureWarnings
- :pr:`6959`: TST: Add test for offset exposure null
- :pr:`6961`: BUG: SARIMAX throwing different errors when length of endogenous var is too low
- :pr:`6962`: BUG: Fix start params computation with few nobs
- :pr:`6963`: TST: Relax tolerance on random failure
- :pr:`6964`: MAINT: Catch warning
- :pr:`6965`: DOC: Final doc fixed for 0.12.0rc0
- :pr:`6968`: DOC: Update DeterministicProcess docs
- :pr:`6969`: ENH add normalized iqr to robust.scales
- :pr:`6972`: DOC: Add docstring to string_like method
- :pr:`6973`: ENH/BUG: robust.scale.iqr does need centering, since quantiles are translation equivariant
- :pr:`6976`: DOC: Fix LaTeX in seasonal notebook
- :pr:`6985`: MAINT: Remove print statements
- :pr:`6990`: ENH: Add robust qn scale
- :pr:`6991`: ENH: Return self from KDEUnivariate fit
- :pr:`6992`: MAINT: Update Azure images
- :pr:`7000`: BUG: fix bug where mad ignores center if center is not callable
- :pr:`7001`: DOC: add new stats to release notes for 0.12
- :pr:`7002`: ENH: Allow custom bandwidth functions in KDEUnivariate fit


API Changes
===========

Notable New Classes
-------------------
* :class:`statsmodels.stats.descriptivestats.Description`
* :class:`statsmodels.stats.meta_analysis.CombineResults`
* :class:`statsmodels.stats.robust_compare.TrimmedMean`
* :class:`statsmodels.tools.sm_exceptions.ParseError`
* :class:`statsmodels.tsa.base.prediction.PredictionResults`
* :class:`statsmodels.tsa.deterministic.CalendarDeterministicTerm`
* :class:`statsmodels.tsa.deterministic.CalendarFourier`
* :class:`statsmodels.tsa.deterministic.CalendarSeasonality`
* :class:`statsmodels.tsa.deterministic.CalendarTimeTrend`
* :class:`statsmodels.tsa.deterministic.DeterministicProcess`
* :class:`statsmodels.tsa.deterministic.DeterministicTerm`
* :class:`statsmodels.tsa.deterministic.Fourier`
* :class:`statsmodels.tsa.deterministic.FourierDeterministicTerm`
* :class:`statsmodels.tsa.deterministic.Seasonality`
* :class:`statsmodels.tsa.deterministic.TimeTrend`
* :class:`statsmodels.tsa.deterministic.TimeTrendDeterministicTerm`
* :class:`statsmodels.tsa.exponential_smoothing.ets.ETSModel`
* :class:`statsmodels.tsa.exponential_smoothing.ets.ETSResults`
* :class:`statsmodels.tsa.forecasting.stl.STLForecast`
* :class:`statsmodels.tsa.forecasting.stl.STLForecastResults`
* :class:`statsmodels.tsa.forecasting.theta.ThetaModel`
* :class:`statsmodels.tsa.forecasting.theta.ThetaModelResults`
* :class:`statsmodels.tsa.statespace.cfa_simulation_smoother.CFASimulationSmoother`
* :class:`statsmodels.tsa.statespace.dynamic_factor_mq.DynamicFactorMQ`
* :class:`statsmodels.tsa.statespace.dynamic_factor_mq.DynamicFactorMQResults`
* :class:`statsmodels.tsa.statespace.news.NewsResults`

Moved or Removed Classes
------------------------
* ``statsmodels.base._penalties.L2ContraintsPenalty``
* ``statsmodels.tools.docstring.ParseError``
* ``statsmodels.tsa.holtwinters.ExponentialSmoothing``
* ``statsmodels.tsa.holtwinters.Holt``
* ``statsmodels.tsa.holtwinters.HoltWintersResults``
* ``statsmodels.tsa.holtwinters.HoltWintersResultsWrapper``
* ``statsmodels.tsa.holtwinters.SimpleExpSmoothing``
* ``statsmodels.tsa.stattools.ResultsStore``

New Methods
-----------
* :meth:`statsmodels.genmod.generalized_linear_model.GLMResults.bic`
* :meth:`statsmodels.tsa.ar_model.AutoRegResults.forecast`
* :meth:`statsmodels.tsa.ar_model.AutoRegResults.get_prediction`
* :meth:`statsmodels.tsa.arima.model.ARIMAResults.append`
* :meth:`statsmodels.tsa.base.prediction.PredictionResults.predicted_mean`
* :meth:`statsmodels.tsa.base.prediction.PredictionResults.row_labels`
* :meth:`statsmodels.tsa.base.prediction.PredictionResults.t_test`
* :meth:`statsmodels.tsa.base.prediction.PredictionResults.tvalues`
* :meth:`statsmodels.tsa.base.prediction.PredictionResults.var_pred_mean`
* :meth:`statsmodels.tsa.statespace.kalman_smoother.SmootherResults.news`
* :meth:`statsmodels.tsa.statespace.kalman_smoother.SmootherResults.smoothed_state_autocovariance`
* :meth:`statsmodels.tsa.statespace.kalman_smoother.SmootherResults.smoothed_state_gain`
* :meth:`statsmodels.tsa.statespace.mlemodel.MLEResults.news`
* :meth:`statsmodels.tsa.statespace.representation.Representation.diff_endog`
* :meth:`statsmodels.tsa.vector_ar.var_model.VAR.from_formula`

Removed Methods
---------------
* ``statsmodels.base.model.GEEResults.remove_data``
* ``statsmodels.base.model.NominalGEEResults.remove_data``
* ``statsmodels.base.model.OrdinalGEEResults.remove_data``
* ``statsmodels.base.model.VAR.from_formula``
* ``statsmodels.genmod._prediction.PredictionResults.se_obs``
* ``statsmodels.genmod._prediction.PredictionResults.t_test``
* ``statsmodels.genmod._prediction.PredictionResults.tvalues``
* ``statsmodels.genmod.generalized_estimating_equations.GEE.predict``
* ``statsmodels.genmod.generalized_estimating_equations.NominalGEE.predict``
* ``statsmodels.genmod.generalized_estimating_equations.OrdinalGEE.predict``
* ``statsmodels.tsa.statespace.mlemodel.ARIMAResults.append``

Methods with New Arguments
--------------------------
* :meth:`statsmodels.discrete.discrete_model.BinaryModel`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.CountModel`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.DiscreteModel`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.GeneralizedPoisson`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.Logit`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.MNLogit`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.MultinomialModel`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.NegativeBinomial`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.NegativeBinomialP`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.Poisson`: ``check_rank``
* :meth:`statsmodels.discrete.discrete_model.Probit`: ``check_rank``
* :meth:`statsmodels.duration.hazard_regression.PHReg`: ``pred_only``
* :meth:`statsmodels.duration.hazard_regression.rv_discrete_float`: ``n``
* :meth:`statsmodels.genmod.cov_struct.Autoregressive`: ``grid``
* :meth:`statsmodels.iolib.summary2.Summary`: ``label``
* :meth:`statsmodels.regression.mixed_linear_model.MixedLMResults`: ``fit_kwargs``
* :meth:`statsmodels.regression.recursive_ls.RecursiveLSResults`: ``copy_initialization``
* :meth:`statsmodels.regression.rolling.RollingOLS`: ``expanding``
* :meth:`statsmodels.regression.rolling.RollingWLS`: ``expanding``
* :meth:`statsmodels.stats.mediation.Mediation`: ``outcome_predict_kwargs``
* :meth:`statsmodels.tsa.ar_model.AutoReg`: ``deterministic``, ``old_names``
* :meth:`statsmodels.tsa.arima.model.ARIMAResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.dynamic_factor.DynamicFactorResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.exponential_smoothing.ExponentialSmoothingResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.kalman_smoother.KalmanSmoother`: ``update_smoother``, ``update_filter``, ``update_representation``
* :meth:`statsmodels.tsa.statespace.mlemodel.MLEResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.sarimax.SARIMAXResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.simulation_smoother.SimulationSmoother`: ``update_smoother``, ``update_filter``, ``update_representation``
* :meth:`statsmodels.tsa.statespace.structural.UnobservedComponentsResults`: ``copy_initialization``
* :meth:`statsmodels.tsa.statespace.varmax.VARMAXResults`: ``truncate_endog_names``
* :meth:`statsmodels.tsa.vector_ar.irf.IRAnalysis`: ``figsize``

Methods with Changed Arguments
------------------------------
* :meth:`statsmodels.regression.mixed_linear_model.MixedLM`
   * New: ``MixedLM(start_params, reml, niter_sa, do_cg, fe_pen, cov_pen, free, full_output, method, fit_kwargs)``
   * Old: ``MixedLM(start_params, reml, niter_sa, do_cg, fe_pen, cov_pen, free, full_output, method, kwargs)``

New Functions
-------------
* :func:`statsmodels.robust.scale.iqr`
* :func:`statsmodels.robust.scale.qn_scale`
* :func:`statsmodels.stats.contrast.wald_test_noncent`
* :func:`statsmodels.stats.contrast.wald_test_noncent_generic`
* :func:`statsmodels.stats.descriptivestats.describe`
* :func:`statsmodels.stats.meta_analysis.combine_effects`
* :func:`statsmodels.stats.meta_analysis.effectsize_2proportions`
* :func:`statsmodels.stats.meta_analysis.effectsize_smd`
* :func:`statsmodels.stats.multivariate.confint_mvmean`
* :func:`statsmodels.stats.multivariate.confint_mvmean_fromstats`
* :func:`statsmodels.stats.multivariate.test_cov`
* :func:`statsmodels.stats.multivariate.test_cov_blockdiagonal`
* :func:`statsmodels.stats.multivariate.test_cov_diagonal`
* :func:`statsmodels.stats.multivariate.test_cov_oneway`
* :func:`statsmodels.stats.multivariate.test_cov_spherical`
* :func:`statsmodels.stats.multivariate.test_mvmean`
* :func:`statsmodels.stats.multivariate.test_mvmean_2indep`
* :func:`statsmodels.stats.oneway.anova_generic`
* :func:`statsmodels.stats.oneway.anova_oneway`
* :func:`statsmodels.stats.oneway.confint_effectsize_oneway`
* :func:`statsmodels.stats.oneway.confint_noncentrality`
* :func:`statsmodels.stats.oneway.convert_effectsize_fsqu`
* :func:`statsmodels.stats.oneway.effectsize_oneway`
* :func:`statsmodels.stats.oneway.equivalence_oneway`
* :func:`statsmodels.stats.oneway.equivalence_oneway_generic`
* :func:`statsmodels.stats.oneway.equivalence_scale_oneway`
* :func:`statsmodels.stats.oneway.f2_to_wellek`
* :func:`statsmodels.stats.oneway.fstat_to_wellek`
* :func:`statsmodels.stats.oneway.power_equivalence_oneway`
* :func:`statsmodels.stats.oneway.simulate_power_equivalence_oneway`
* :func:`statsmodels.stats.oneway.test_scale_oneway`
* :func:`statsmodels.stats.oneway.wellek_to_f2`
* :func:`statsmodels.stats.power.normal_power_het`
* :func:`statsmodels.stats.power.normal_sample_size_one_tail`
* :func:`statsmodels.stats.proportion.confint_proportions_2indep`
* :func:`statsmodels.stats.proportion.power_proportions_2indep`
* :func:`statsmodels.stats.proportion.samplesize_proportions_2indep_onetail`
* :func:`statsmodels.stats.proportion.score_test_proportions_2indep`
* :func:`statsmodels.stats.proportion.test_proportions_2indep`
* :func:`statsmodels.stats.proportion.tost_proportions_2indep`
* :func:`statsmodels.stats.rates.etest_poisson_2indep`
* :func:`statsmodels.stats.rates.test_poisson_2indep`
* :func:`statsmodels.stats.rates.tost_poisson_2indep`
* :func:`statsmodels.stats.robust_compare.scale_transform`
* :func:`statsmodels.stats.robust_compare.trim_mean`
* :func:`statsmodels.stats.robust_compare.trimboth`
* :func:`statsmodels.tools.eval_measures.rmspe`


Removed Functions
-----------------
* ``statsmodels.compat.python.iteritems``
* ``statsmodels.compat.python.iterkeys``
* ``statsmodels.compat.python.itervalues``
* ``statsmodels.stats.diagnostic.unitroot_adf``
* ``statsmodels.tools.decorators.nottest``
* ``statsmodels.tsa.stattools.periodogram``

Functions with New Arguments
----------------------------
* :func:`statsmodels.graphics.gofplots.qqline`: ``lineoptions``
* :func:`statsmodels.nonparametric.smoothers_lowess.lowess`: ``xvals``
* :func:`statsmodels.stats.diagnostic.acorr_ljungbox`: ``auto_lag``
* :func:`statsmodels.stats.diagnostic.het_breuschpagan`: ``robust``
* :func:`statsmodels.stats.diagnostic.linear_harvey_collier`: ``skip``
* :func:`statsmodels.stats.multitest.local_fdr`: ``alpha``
* :func:`statsmodels.tsa.ar_model.ar_select_order`: ``old_names``

Functions with Changed Arguments
--------------------------------
* :func:`statsmodels.graphics.tsaplots.plot_acf`
   * New: ``plot_acf(x, ax, lags, alpha, use_vlines, adjusted, fft, missing, title, zero, vlines_kwargs, kwargs)``
   * Old: ``plot_acf(x, ax, lags, alpha, use_vlines, unbiased, fft, missing, title, zero, vlines_kwargs, kwargs)``
* :func:`statsmodels.tsa.stattools.acf`
   * New: ``acf(x, adjusted, nlags, qstat, fft, alpha, missing)``
   * Old: ``acf(x, unbiased, nlags, qstat, fft, alpha, missing)``
* :func:`statsmodels.tsa.stattools.acovf`
   * New: ``acovf(x, adjusted, demean, fft, missing, nlag)``
   * Old: ``acovf(x, unbiased, demean, fft, missing, nlag)``
* :func:`statsmodels.tsa.stattools.ccf`
   * New: ``ccf(x, y, adjusted)``
   * Old: ``ccf(x, y, unbiased)``
* :func:`statsmodels.tsa.stattools.ccovf`
   * New: ``ccovf(x, y, adjusted, demean)``
   * Old: ``ccovf(x, y, unbiased, demean)``
* :func:`statsmodels.tsa.stattools.pacf_ols`
   * New: ``pacf_ols(x, nlags, efficient, adjusted)``
   * Old: ``pacf_ols(x, nlags, efficient, unbiased)``
:orphan:

==============
Release 0.12.1
==============

Release summary
===============
This is a bug fix release.

Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Kevin Sheppard
* Chad Fulton
* Josef Perktold
* Kerby Shedden
* Pratyush Sharan

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

* :pr:`7016`: BLD: avoid setuptools version 50, windows build problem
* :pr:`7017`: BUG: param names with higher order trend in VARMAX
* :pr:`7020`: REF: Don't validate specification in SARIMAX when cloning to get extended time varying matrices
* :pr:`7025`: BUG: Ensure bestlag is defined in autolag
* :pr:`7028`: BUG: Correct axis None case
* :pr:`7040`: Bug fix ets get prediction
* :pr:`7052`: ENH: handle/warn for singularities in MixedLM
* :pr:`7055`: BUG: Fix squeeze when nsimulation is 1
* :pr:`7073`: DOC: fix several doc issues in stats functions
* :pr:`7088`: DOC: augmented docstrings from statsmodels.base.optimizer
* :pr:`7090`: DOC: Fix contradicting KPSS-statistics interpretations in stationarity_detrending_adf_kpss.ipynb
* :pr:`7093`: BUG: Correct prediction intervals for ThetaModel
* :pr:`7109`: some fixes in the doc of grangercausalitytests
* :pr:`7116`: BUG: don't raise error in impacts table if no news.
* :pr:`7118`: MAINT: Fix issues in main branches of dependencies
:orphan:

=============
Release 0.9.0
=============

Release summary
---------------

statsmodels is using github to store the updated documentation which
is available under
https://www.statsmodels.org/stable for the last release, and
https://www.statsmodels.org/devel/ for the development version.


**Warning**

API stability is not guaranteed for new features, although even in
this case changes will be made in a backwards compatible way if
possible. The stability of a new feature depends on how much time it
was already in statsmodels main and how much usage it has already
seen.  If there are specific known problems or limitations, then they
are mentioned in the docstrings.


The list of pull requests for this release can be found on github
https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.9
(The list does not include some pull request that were merged before
the 0.8 release but not included in 0.8.)


The Highlights
--------------

- statespace refactoring, Markov Switching Kim smoother
- 3 Google summer of code (GSOC) projects merged
  - distributed estimation
  - VECM and enhancements to VAR (including cointegration test)
  - new count models: GeneralizedPoisson, zero inflated models
- Bayesian mixed GLM
- Gaussian Imputation
- new multivariate methods: factor analysis, MANOVA, repeated measures
  within ANOVA
- GLM var_weights in addition to freq_weights
- Holt-Winters and Exponential Smoothing


What's new - an overview
------------------------

The following lists the main new features of statsmodels 0.9. In addition,
release 0.9 includes bug fixes, refactorings and improvements in many areas.

**base**
 - distributed estimation #3396  (Leland Bybee GSOC, Kerby Shedden)
 - optimization option scipy minimize #3193 (Roman Ring)
 - Box-Cox #3477 (Niels Wouda)
 - t_test_pairwise #4365 (Josef Perktold)

**discrete**
 - new count models (Evgeny Zhurko GSOC, Josef Perktold)
    - NegativeBinomialP #3832 merged in #3874
    - GeneralizedPoisson #3727 merged in  #3795
    - zero-inflated count models #3755 merged in #3908

 - discrete optimization improvements #3921, #3928 (Josef Perktold)
 - extend discrete margin when extra params, NegativeBinomial #3811
   (Josef Perktold)

**duration**
 - dependent censoring in survival/duration #3090 (Kerby Shedden)
 - entry times for Kaplan-Meier #3126 (Kerby Shedden)

**genmod**
 - Bayesian GLMM #4189, #4540 (Kerby Shedden)
 - GLM add var_weights #3692 (Peter Quackenbush)
 - GLM: EIM in optimization #3646 (Peter Quackenbush)
 - GLM correction to scale handling, loglike #3856 (Peter Quackenbush)

**graphics**
 - graphics HDR functional boxplot #3876 merged in #4049 (Pamphile ROY)
 - graphics Bland-Altman or Tukey mean difference plot
   #4112 merged in #4200 (Joses W. Ho)
 - bandwidth options in violinplots #4510 (Jim Correia)

**imputation**
 - multiple imputation via Gaussian model #4394, #4520 (Kerby Shedden)
 - regularized fitting in MICE #4319 (Kerby Shedden)

**iolib**
 - improvements of summary_coll #3702 merged #4064 (Natasha Watkins,
   Kevin Sheppard)

**multivariate**
 - multivariate: MANOVA, CanCorr #3327 (Yichuan Liu)
 - Factor Analysis #4161, #4156, #4167, #4214 (Yichuan Liu, Kerby Shedden,
   Josef Perktold)
 - statsmodels now includes the rotation code by ....

**regression**
 - fit_regularized for WLS #3581 (Kerby Shedden)

**stats**
 - Knockoff FDR # 3204 (Kerby Shedden)
 - Repeated measures ANOVA #3303 merged in #3663, #3838 (Yichuan Liu, Richard
   Höchenberger)
 - lilliefors test for exponential distribution #3837 merged in #3936 (Jacob
   Kimmel, Josef Perktold)

**tools**
 - quasi-random, Halton sequences #4104 (Pamphile ROY)

**tsa**
 - VECM #3246 (Aleksandar Karakas GSOC, Josef Perktold)
 - exog support in VAR, incomplete for extra results, part of VECM
   #3246, #4538 (Aleksandar Karakas GSOC, Josef Perktold)
 - Markov switching, Kim smoother #3141 (Chad Fulton)
 - Holt-Winters #3817 merged in #4176 (tvanzyl)
 - seasonal_decompose: trend extrapolation and vectorized 2-D #3031
   (kernc, Josef Perktold)
 - add frequency domain seasonal components to UnobservedComponents #4250
   (Jordan Yoder)
 - refactoring of date handling in tsa #3276, #4457 (Chad Fulton)
 - SARIMAX without AR, MA #3383  (Chad Fulton)

**maintenance**
 - switch to pytest #3804 plus several other PRs (Kevin Sheppard)
 - general compatibility fixes for recent versions of numpy, scipy and pandas


`bug-wrong`
~~~~~~~~~~~

A new issue label `type-bug-wrong` indicates bugs that cause that incorrect
numbers are returned without warnings.
(Regular bugs are mostly usability bugs or bugs that raise an exception for
unsupported use cases.)
see https://github.com/statsmodels/statsmodels/issues?q=is%3Aissue+label%3Atype-bug-wrong+is%3Aclosed+milestone%3A0.9

- scale in GLM fit_constrained, #4193 fixed in #4195
  cov_params and bse were incorrect if scale is estimated as in Gaussian.
  (This did not affect families with scale=1 such as Poisson)
- incorrect `pearson_chi2` with binomial counts, #3612 fixed as part of #3692
- null_deviance and llnull in GLMResults were wrong if exposure was used and
  when offset was used with Binomial counts.
- GLM Binomial in the non-binary count case used incorrect endog in recreating
  models which is
  used by fit_regularized and fit_constrained #4599.
- GLM observed hessian was incorrectly computed if non-canonical link is used,
  fixed in #4620
  This fix improves convergence with gradient optimization and removes a usually
  numerically small error in cov_params.
- discrete predict with offset or exposure, #3569 fixed in #3696
  If either offset or exposure are not None but exog is None, then offset and
  exposure arguments in predict were ignored.
- discrete margins had wrong dummy and count effect if constant is prepended,
  #3695 fixed in #3696
- OLS outlier test, wrong index if order is True, #3971 fixed in #4385
- tsa coint ignored the autolag keyword, #3966 fixed in #4492
  This is a backwards incompatible change in default, instead of fixed maxlag
  it defaults now to 'aic' lag selection. The default autolag is now the same
  as the adfuller default.
- wrong confidence interval in contingency table summary, #3822 fixed in #3830
  This only affected the summary and not the corresponding attribute.
- incorrect results in summary_col if regressor_order is used,
  #3767 fixed in #4271


Description of selected new feature
-----------------------------------

The following provides more information about a selected set of new features.

Vector Error Correction Model (VECM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The VECM framework developed during GSOC 2016 by Aleksandar Karakas adds support
for non-stationary cointegrated VAR processes to statsmodels.
Currently, the following topics are implemented

* Parameter estimation for cointegrated VAR
* forecasting
* testing for Granger-causality and instantaneous causality
* testing for cointegrating rank
* lag order selection.

New methods have been added also to the existing VAR model, and VAR has now
limited support for user provided explanatory variables.


New Count Models
----------------

New count models have been added as part of GSOC 2017 by Evgeny Zhurko.
Additional models that are not yet finished will be added for the next release.

The new models are:

* NegativeBinomialP (NBP): This is a generalization of NegativeBinomial that
  allows the variance power parameter to be specified in the range between 1
  and 2. The current NegativeBinomial support NB1 and NB2 which are two special
  cases of NBP.
* GeneralizedPoisson (GPP): Similar to NBP this allows a large range of
  dispersion specification. GPP also allow some amount of under dispersion
* ZeroInflated Models: Based on a generic base class, zeroinflated models
  are now available for Poisson, GeneralizedPoisson and NegativeBinomialP.

Generalized linear mixed models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Limited support for GLIMMIX models is now included in the genmod
module.  Binomial and Poisson models with independent random effects
can be fit using Bayesian methods (Laplace and mean field
approximations to the posterior).

Multiple imputation
~~~~~~~~~~~~~~~~~~~

Multiple imputation using a multivariate Gaussian model is now
included in the imputation module.  The model is fit via Gibbs
sampling from the joint posterior of the mean vector, covariance
matrix, and missing data values.  A convenience function for fitting a
model to the multiply imputed data sets and combining the results is
provided.  This is an alternative to the existing MICE (Multiple
Imputation via Chained Equations) procedures.

Exponential smoothing models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exponential smoothing models are now available (introduced in #4176 by
Terence L van Zyl). These models are conceptually simple, decomposing a time
series into level, trend, and seasonal components that are constructed from
weighted averages of past observations. Nonetheless, they produce forecasts
that are competitive with more advanced models and which may be easier to
interpret.

Available models include:

- Simple exponential smoothing
- Holt's method
- Holt-Winters exponential smoothing

Improved time series index support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Handling of indexes for time series models has been overhauled (#3272) to
take advantage of recent improvements in Pandas and to shift to Pandas much of
the special case handling (especially for date indexes) that had previously been
done in statsmodels. Benefits include more consistent behavior, a reduced
number of bugs from corner cases, and a reduction in the maintenance burden.

Although an effort was made to maintain backwards compatibility with this
change, it is possible that some undocumented corner cases that previously
worked will now raise warnings or exceptions.

State space models
~~~~~~~~~~~~~~~~~~

The state space model infrastructure has been rewritten and improved (#2845).
New features include:

- Kalman smoother rewritten in Cython for substantial performance improvements
- Simulation smoother (Durbin and Koopman, 2002)
- Fast simulation of time series for any state space model
- Univariate Kalman filtering and smoothing (Koopman and Durbin, 2000)
- Collapsed Kalman filtering and smoothing (Jungbacker and Koopman, 2014)
- Optional computation of the lag-one state autocovariance
- Use of the Scipy BLAS functions for Cython interface if available
  (`scipy.linalg.cython_blas` for Scipy >= 0.16)

These features yield new features and improve performance for the existing
state space models (`SARIMAX`, `UnobservedComopnents`, `DynamicFactor`, and
`VARMAX`), and they also make Bayesian estimation by Gibbs-sampling possible.

**Warning**: this will be the last version that includes the original state
space code and supports Scipy < 0.16. The next release will only include the
new state space code.

Unobserved components models: frequency-domain seasonals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unobserved components models now support modeling seasonal factors from a
frequency-domain perspective with user-specified period and harmonics
(introduced in #4250 by Jordan Yoder). This not only allows for multiple
seasonal effects, but also allows the representation of seasonal components
with fewer unobserved states. This can improve computational performance and,
since it allows for a more parsimonious model, may also improve the
out-of-sample performance of the model.


Major Bugs fixed
----------------

* see github issues for a list of bug fixes included in this release
  https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.9+label%3Atype-bug
  https://github.com/statsmodels/statsmodels/pulls?q=is%3Apr+is%3Amerged+milestone%3A0.9+label%3Atype-bug-wrong

* Refitting elastic net regularized models using the `refit=True`
  option now returns the unregularized parameters for the coefficients
  selected by the regularized fitter, as documented. #4213

* In MixedLM, a bug that produced exceptions when calling
  `random_effects_cov` on models with variance components has been
  fixed.


Backwards incompatible changes and deprecations
-----------------------------------------------

* DynamicVAR and DynamicPanelVAR is deprecated and will be removed in
  a future version. It used rolling OLS from pandas which has been
  removed in pandas.

* In MixedLM, names for the random effects variance and covariance
  parameters have changed from, e.g. G RE to G Var or G x F Cov.  This
  impacts summary output, and also may require modifications to user
  code that extracted these parameters from the fitted results object
  by name.

* In MixedLM, the names for the random effects realizations for
  variance components have been changed.  When using formulas, the
  random effect realizations are named using the column names produced
  by Patsy when parsing the formula.


Development summary and credits
-------------------------------

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Kevin Sheppard
* Peter Quackenbush
* Brock Mendel

and the general maintainer and code reviewer

* Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.9 release (based on git log):

.. note::

    * Aleksandar Karakas
    * Alex Fortin
    * Alexander Belopolsky
    * Brock Mendel
    * Chad Fulton
    * ChadFulton
    * Christian Lorentzen
    * Dave Willmer
    * Dror Atariah
    * Evgeny Zhurko
    * Gerard Brunick
    * Greg Mosby
    * Jacob Kimmel
    * Jamie Morton
    * Jarvis Miller
    * Jasmine Mou
    * Jeroen Van Goey
    * Jim Correia
    * Joon Ro
    * Jordan Yoder
    * Jorge C. Leitao
    * Josef Perktold
    * Joses W. Ho
    * José Lopez
    * Joshua Engelman
    * Juan Escamilla
    * Justin Bois
    * Kerby Shedden
    * Kernc
    * Kevin Sheppard
    * Leland Bybee
    * Maxim Uvarov
    * Michael Kaminsky
    * Mosky Liu
    * Natasha Watkins
    * Nick DeRobertis
    * Niels Wouda
    * Pamphile ROY
    * Peter Quackenbush
    * Quentin Andre
    * Richard Höchenberger
    * Rob Klooster
    * Roman Ring
    * Scott Tsai
    * Soren Fuglede Jorgensen
    * Tom Augspurger
    * Tommy Odland
    * Tony Jiang
    * Yichuan Liu
    * ftemme
    * hugovk
    * kiwirob
    * malickf
    * tvanzyl
    * weizhongg
    * zveryansky

These lists of names are automatically generated based on git log, and may not
be complete.
:orphan:

.. _old_changes:

Pre 0.5.0 Release History
=========================

0.5.0
-----
*Main Changes and Additions*
* Add patsy dependency

*Compatibility and Deprecation*

* cleanup of import paths (lowess)
*

*Bug Fixes*

* input shapes of tools.isestimable
*

*Enhancements and Additions*

* formula integration based on patsy (new dependency)
* Time series analysis
  - ARIMA modeling
  - enhanced forecasting based on pandas datetime handling
* expanded margins for discrete models
* OLS outlier test

* empirical likelihood - Google Summer of Code 2012 project
  - inference for descriptive statistics
  - inference for regression models
  - accelerated failure time models

* expanded probability plots
* improved graphics
  - plotcorr
  - acf and pacf
* new datasets
* new and improved tools
  - numdiff numerical differentiation



0.4.3
-----

The only change compared to 0.4.2 is for compatibility with python 3.2.3
(changed behavior of 2to3)


0.4.2
-----

This is a bug-fix release, that affects mainly Big-Endian machines.

*Bug Fixes*

* discrete_model.MNLogit fix summary method
* tsa.filters.hp_filter do not use umfpack on Big-Endian machine (scipy bug)
* the remaining fixes are in the test suite, either precision problems
  on some machines or incorrect testing on Big-Endian machines.



0.4.1
-----

This is a backwards compatible (according to our test suite) release with
bug fixes and code cleanup.

*Bug Fixes*

* build and distribution fixes
* lowess correct distance calculation
* genmod correction CDFlink derivative
* adfuller _autolag correct calculation of optimal lag
* het_arch, het_lm : fix autolag and store options
* GLSAR: incorrect whitening for lag>1

*Other Changes*

* add lowess and other functions to api and documentation
* rename lowess module (old import path will be removed at next release)
* new robust sandwich covariance estimators, moved out of sandbox
* compatibility with pandas 0.8
* new plots in statsmodels.graphics
  - ABLine plot
  - interaction plot


0.4.0
-----

*Main Changes and Additions*

* Added pandas dependency.
* Cython source is built automatically if cython and compiler are present
* Support use of dates in timeseries models
* Improved plots
  - Violin plots
  - Bean Plots
  - QQ Plots
* Added lowess function
* Support for pandas Series and DataFrame objects. Results instances return
  pandas objects if the models are fit using pandas objects.
* Full Python 3 compatibility
* Fix bugs in genfromdta. Convert Stata .dta format to structured array
  preserving all types. Conversion is much faster now.
* Improved documentation
* Models and results are pickleable via save/load, optionally saving the model
  data.
* Kernel Density Estimation now uses Cython and is considerably faster.
* Diagnostics for outlier and influence statistics in OLS
* Added El Nino Sea Surface Temperatures dataset
* Numerous bug fixes
* Internal code refactoring
* Improved documentation including examples as part of HTML

*Changes that break backwards compatibility*

* Deprecated scikits namespace. The recommended import is now::

      import statsmodels.api as sm

* model.predict methods signature is now (params, exog, ...) where before
  it assumed that the model had been fit and omitted the params argument.
* For consistency with other multi-equation models, the parameters of MNLogit
  are now transposed.
* tools.tools.ECDF -> distributions.ECDF
* tools.tools.monotone_fn_inverter -> distributions.monotone_fn_inverter
* tools.tools.StepFunction -> distributions.StepFunction


0.3.1
-----

* Removed academic-only WFS dataset.
* Fix easy_install issue on Windows.

0.3.0
-----

*Changes that break backwards compatibility*

Added api.py for importing. So the new convention for importing is::

    import statsmodels.api as sm

Importing from modules directly now avoids unnecessary imports and increases
the import speed if a library or user only needs specific functions.

* sandbox/output.py -> iolib/table.py
* lib/io.py -> iolib/foreign.py (Now contains Stata .dta format reader)
* family -> families
* families.links.inverse -> families.links.inverse_power
* Datasets' Load class is now load function.
* regression.py -> regression/linear_model.py
* discretemod.py -> discrete/discrete_model.py
* rlm.py -> robust/robust_linear_model.py
* glm.py -> genmod/generalized_linear_model.py
* model.py -> base/model.py
* t() method -> tvalues attribute (t() still exists but raises a warning)

*Main changes and additions*

* Numerous bugfixes.
* Time Series Analysis model (tsa)

  - Vector Autoregression Models VAR (tsa.VAR)
  - Autoregressive Models AR (tsa.AR)
  - Autoregressive Moving Average Models ARMA (tsa.ARMA)
    optionally uses Cython for Kalman Filtering
    use setup.py install with option --with-cython
  - Baxter-King band-pass filter (tsa.filters.bkfilter)
  - Hodrick-Prescott filter (tsa.filters.hpfilter)
  - Christiano-Fitzgerald filter (tsa.filters.cffilter)

* Improved maximum likelihood framework uses all available scipy.optimize solvers
* Refactor of the datasets sub-package.
* Added more datasets for examples.
* Removed RPy dependency for running the test suite.
* Refactored the test suite.
* Refactored codebase/directory structure.
* Support for offset and exposure in GLM.
* Removed data_weights argument to GLM.fit for Binomial models.
* New statistical tests, especially diagnostic and specification tests
* Multiple test correction
* General Method of Moment framework in sandbox
* Improved documentation
* and other additions


0.2.0
-----

*Main changes*

 * renames for more consistency
   RLM.fitted_values -> RLM.fittedvalues
   GLMResults.resid_dev -> GLMResults.resid_deviance
 * GLMResults, RegressionResults:
   lazy calculations, convert attributes to properties with _cache
 * fix tests to run without rpy
 * expanded examples in examples directory
 * add PyDTA to lib.io -- functions for reading Stata .dta binary files
   and converting
   them to numpy arrays
 * made tools.categorical much more robust
 * add_constant now takes a prepend argument
 * fix GLS to work with only a one column design

*New*

 * add four new datasets

   - A dataset from the American National Election Studies (1996)
   - Grunfeld (1950) investment data
   - Spector and Mazzeo (1980) program effectiveness data
   - A US macroeconomic dataset
 * add four new Maximum Likelihood Estimators for models with a discrete
   dependent variables with examples

   - Logit
   - Probit
   - MNLogit (multinomial logit)
   - Poisson

*Sandbox*

 * add qqplot in sandbox.graphics
 * add sandbox.tsa (time series analysis) and sandbox.regression (anova)
 * add principal component analysis in sandbox.tools
 * add Seemingly Unrelated Regression (SUR) and Two-Stage Least Squares
   for systems of equations in sandbox.sysreg.Sem2SLS
 * add restricted least squares (RLS)


0.1.0b1
-------
 * initial release
:orphan:

=============
Release 0.7.0
=============

Release summary
---------------

**Note:** This version has never been officially released. Several models have
been refactored, improved or bugfixed in 0.8.


The following major new features appear in this version.

Principal Component Analysis
----------------------------

Author: Kevin Sheppard

A new class-based Principal Component Analysis has been added.  This
class replaces the function-based PCA that previously existed in the
sandbox.  This change bring a number of new features, including:

* Options to control the standardization (demeaning/studentizing)
* Scree plotting
* Information criteria for selecting the number of factors
* R-squared plots to assess component fit
* NIPALS implementation when only a small number of components are required and the dataset is large
* Missing-value filling using the EM algorithm

.. code-block:: python

   import statsmodels.api as sm
   from statsmodels.multivariate.pca import PCA

   data = sm.datasets.fertility.load_pandas().data

   columns = map(str, range(1960, 2012))
   data.set_index('Country Name', inplace=True)
   dta = data[columns]
   dta = dta.dropna()

   pca_model = PCA(dta.T, standardize=False, demean=True)
   pca_model.plot_scree()

*Note* : A function version is also available which is compatible with the
call in the sandbox.  The function version is just a thin wrapper around the
class-based PCA implementation.

Regression graphics for GLM/GEE
-------------------------------

Author: Kerby Shedden

Added variable plots, partial residual plots, and CERES residual plots
are available for GLM and GEE models by calling the methods
`plot_added_variable`, `plot_partial_residuals`, and
`plot_ceres_residuals` that are attached to the results classes.

State Space Models
------------------

Author: Chad Fulton

State space methods provide a flexible structure for the estimation and
analysis of a wide class of time series models. The statsmodels implementation
allows specification of state models, fast Kalman filtering, and built-in
methods to facilitate maximum likelihood estimation of arbitrary models. One of
the primary goals of this module is to allow end users to create and estimate
their own models. Below is a short example demonstrating the ease with which a
local level model can be specified and estimated:

.. code-block:: python

   import numpy as np
   import statsmodels.api as sm
   import pandas as pd

   data = sm.datasets.nile.load_pandas().data
   data.index = pd.DatetimeIndex(data.year.astype(int).astype(str), freq='AS')

   # Setup the state space representation
   class LocalLevel(sm.tsa.statespace.MLEModel):
       def __init__(self, endog):
           # Initialize the state space model
           super(LocalLevel, self).__init__(
               endog, k_states=1, initialization='approximate_diffuse')

           # Setup known components of state space representation matrices
           self.ssm['design', :] = 1.
           self.ssm['transition', :] = 1.
           self.ssm['selection', :] = 1.

       # Describe how parameters enter the model
       def update(self, params, transformed=True):
           params = super(LocalLevel, self).update(params, transformed)
           self.ssm['obs_cov', 0, 0] = params[0]
           self.ssm['state_cov', 0, 0] = params[1]

       def transform_params(self, params):
           return params**2  # force variance parameters to be positive

       # Specify start parameters and parameter names
       @property
       def start_params(self):
           return [np.std(self.endog)]*2

       @property
       def param_names(self):
           return ['sigma2.measurement', 'sigma2.level']

   # Fit the model with maximum likelihood estimation
   mod = LocalLevel(data['volume'])
   res = mod.fit()
   print res.summary()

The documentation and example notebooks provide further examples of how to
form state space models. Included in this release is a full-fledged
model making use of the state space infrastructure to estimate SARIMAX
models. See below for more details.

Time Series Models (ARIMA) with Seasonal Effects
------------------------------------------------

Author: Chad Fulton

A model for estimating seasonal autoregressive integrated moving average models
with exogenous regressors (SARIMAX) has been added by taking advantage of the
new state space functionality. It can be used very similarly to the existing
`ARIMA` model, but works on a wider range of specifications, including:

* Additive and multiplicative seasonal effects
* Flexible trend specification
* Regression with SARIMA errors
* Regression with time-varying coefficients
* Measurement error in the endogenous variables

Below is a short example fitting a model with a number of these components,
including exogenous data, a linear trend, and annual multiplicative seasonal
effects.

.. code-block:: python

   import statsmodels.api as sm
   import pandas as pd

   data = sm.datasets.macrodata.load_pandas().data
   data.index = pd.DatetimeIndex(start='1959-01-01', end='2009-09-01',
                                 freq='QS')
   endog = data['realcons']
   exog = data['m1']

   mod = sm.tsa.SARIMAX(endog, exog=exog, order=(1,1,1),
                        trend='t', seasonal_order=(0,0,1,4))
   res = mod.fit()
   print res.summary()


Generalized Estimating Equations GEE
------------------------------------

Author: Kerby Shedden

Enhancements and performance improvements for GEE:

* EquivalenceClass covariance structure allows covariances to be specified by
  arbitrary collections of equality constraints #2188
* add weights #2090
* refactored margins #2158


MixedLM
-------

Author: Kerby Shedden with Saket Choudhary

Enhancements to MixedLM (#2363): added variance components support for
MixedLM allowing a wider range of random effects structures to be specified;
also performance improvements from use of sparse matrices internally for
random effects design matrices.


Other important new features
----------------------------

* GLM: add scipy-based gradient optimization to fit #1961 (Kerby Shedden)
* wald_test_terms: new method of LikelihoodModels to compute wald tests (F or chi-square)
  for terms or sets of coefficients #2132  (Josef Perktold)
* add cov_type with fixed scale in WLS to allow chi2-fitting #2137 #2143
  (Josef Perktold, Christoph Deil)
* VAR: allow generalized IRF and FEVD computation #2067 (Josef Perktold)
* get_prediction new method for full prediction results (new API convention)



Major Bugs fixed
----------------

* see github issues for a full list
* bug in ARMA/ARIMA predict with `exog` #2470
* bugs in VAR
* x13: python 3 compatibility



Backwards incompatible changes and deprecations
-----------------------------------------------

* List backwards incompatible changes


Development summary and credits
-------------------------------



.. note::

  Thanks to all of the contributors for the 0.7 release:

.. note::

   * Alex Griffing
   * Antony Lee
   * Chad Fulton
   * Christoph Deil
   * Daniel Sullivan
   * Hans-Martin von Gaudecker
   * Jan Schulz
   * Joey Stockermans
   * Josef Perktold
   * Kerby Shedden
   * Kevin Sheppard
   * Kiyoto Tamura
   * Louis-Philippe Lemieux Perreault
   * Padarn Wilson
   * Ralf Gommers
   * Saket Choudhary
   * Skipper Seabold
   * Tom Augspurger
   * Trent Hauck
   * Vincent Arel-Bundock
   * chebee7i
   * donbeo
   * gliptak
   * hlin117
   * jerry dumblauskas
   * jonahwilliams
   * kiyoto
   * neilsummers
   * waynenilsen

These lists of names are automatically generated based on git log, and may not be
complete.
:orphan:

=============
Release 0.6.1
=============

statsmodels 0.6.1 is a bugfix release. All users are encouraged to upgrade to 0.6.1.

See the :ref:`list of fixed issues <issues_list_06>` for specific backported fixes.

=============
Release 0.6.0
=============

statsmodels 0.6.0 is another large release. It is the result of the work of 37 authors over the last year and includes over 1500 commits. It contains many new features, improvements, and bug fixes detailed below.


See the :ref:`list of fixed issues <issues_list_06>` for specific closed issues.

The following major new features appear in this version.

Generalized Estimating Equations
--------------------------------

Generalized Estimating Equations (GEE) provide an approach to handling
dependent data in a regression analysis.  Dependent data arise
commonly in practice, such as in a longitudinal study where repeated
observations are collected on subjects. GEE can be viewed as an
extension of the generalized linear modeling (GLM) framework to the
dependent data setting.  The familiar GLM families such as the
Gaussian, Poisson, and logistic families can be used to accommodate
dependent variables with various distributions.

Here is an example of GEE Poisson regression in a data set with four
count-type repeated measures per subject, and three explanatory
covariates.

.. code-block:: python

   import numpy as np
   import statsmodels.api as sm
   import statsmodels.formula.api as smf

   data = sm.datasets.get_rdataset("epil", "MASS").data

   md = smf.gee("y ~ age + trt + base", "subject", data,
                cov_struct=sm.cov_struct.Independence(), 
                family=sm.families.Poisson())
   mdf = md.fit()
   print mdf.summary()


The dependence structure in a GEE is treated as a nuisance parameter
and is modeled in terms of a "working dependence structure".  The
statsmodels GEE implementation currently includes five working
dependence structures (independent, exchangeable, autoregressive,
nested, and a global odds ratio for working with categorical data).
Since the GEE estimates are not maximum likelihood estimates,
alternative approaches to some common inference procedures have been
developed.  The statsmodels GEE implementation currently provides
standard errors, Wald tests, score tests for arbitrary parameter
contrasts, and estimates and tests for marginal effects.  Several
forms of standard errors are provided, including robust standard
errors that are approximately correct even if the working dependence
structure is misspecified.

Seasonality Plots
-----------------

Adding functionality to look at seasonality in plots. Two new functions are :func:`sm.graphics.tsa.month_plot` and :func:`sm.graphics.tsa.quarter_plot`. Another function :func:`sm.graphics.tsa.seasonal_plot` is available for power users.

.. code-block:: python

    import statsmodels.api as sm
    import pandas as pd

    dta = sm.datasets.elnino.load_pandas().data
    dta['YEAR'] = dta.YEAR.astype(int).astype(str)
    dta = dta.set_index('YEAR').T.unstack()
    dates = map(lambda x : pd.datetools.parse('1 '+' '.join(x)),
                                           dta.index.values)

    dta.index = pd.DatetimeIndex(dates, freq='M')
    fig = sm.tsa.graphics.month_plot(dta)

.. currentmodule:: statsmodels.tsa

Seasonal Decomposition
----------------------

We added a naive seasonal decomposition tool in the same vein as R's ``decompose``. This function can be found as :func:`sm.tsa.seasonal_decompose <tsa.seasonal.seasonal_decompose>`.


.. plot::
   :include-source:

    import statsmodels.api as sm

    dta = sm.datasets.co2.load_pandas().data
    # deal with missing values. see issue
    dta.co2.interpolate(inplace=True)

    res = sm.tsa.seasonal_decompose(dta.co2)
    res.plot()


Addition of Linear Mixed Effects Models (MixedLM)

Linear Mixed Effects Models
---------------------------

Linear Mixed Effects models are used for regression analyses involving
dependent data.  Such data arise when working with longitudinal and
other study designs in which multiple observations are made on each
subject.  Two specific mixed effects models are "random intercepts
models", where all responses in a single group are additively shifted
by a value that is specific to the group, and "random slopes models",
where the values follow a mean trajectory that is linear in observed
covariates, with both the slopes and intercept being specific to the
group.  The statsmodels MixedLM implementation allows arbitrary random
effects design matrices to be specified for the groups, so these and
other types of random effects models can all be fit.

Here is an example of fitting a random intercepts model to data from a
longitudinal study:

.. code-block:: python

    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    data = sm.datasets.get_rdataset('dietox', 'geepack', cache=True).data
    md = smf.mixedlm("Weight ~ Time", data, groups=data["Pig"])
    mdf = md.fit()
    print mdf.summary()

The statsmodels LME framework currently supports post-estimation
inference via Wald tests and confidence intervals on the coefficients,
profile likelihood analysis, likelihood ratio testing, and AIC.  Some
limitations of the current implementation are that it does not support
structure more complex on the residual errors (they are always
homoscedastic), and it does not support crossed random effects.  We
hope to implement these features for the next release.

Wrapping X-12-ARIMA/X-13-ARIMA
------------------------------

It is now possible to call out to X-12-ARIMA or X-13ARIMA-SEATS from statsmodels. These libraries must be installed separately.

.. code-block:: python

    import statsmodels.api as sm

    dta = sm.datasets.co2.load_pandas().data
    dta.co2.interpolate(inplace=True)
    dta = dta.resample('M').last()

    res = sm.tsa.x13_arima_select_order(dta.co2)
    print(res.order, res.sorder)

    results = sm.tsa.x13_arima_analysis(dta.co2)

    fig = results.plot()
    fig.set_size_inches(12, 5)
    fig.tight_layout()


Other important new features
----------------------------

* The Kalman filter Cython code underlying AR(I)MA estimation has been substantially optimized. You can expect speed-ups of one to two orders of magnitude.

* Added :func:`sm.tsa.arma_order_select_ic`. A convenience function to quickly get the information criteria for use in tentative order selection of ARMA processes.

* Plotting functions for timeseries is now imported under the ``sm.tsa.graphics`` namespace in addition to ``sm.graphics.tsa``.

* New `distributions.ExpandedNormal` class implements the Edgeworth expansion for weakly non-normal distributions.

* **New datasets**: Added new :ref:`datasets <datasets>` for examples. ``sm.datasets.co2`` is a univariate time-series dataset of weekly co2 readings. It exhibits a trend and seasonality and has missing values.

* Added robust skewness and kurtosis estimators in :func:`sm.stats.stattools.robust_skewness` and :func:`sm.stats.stattools.robust_kurtosis`, respectively.  An alternative robust measure of skewness has been added in :func:`sm.stats.stattools.medcouple`.

* New functions added to correlation tools: `corr_nearest_factor`
  finds the closest factor-structured correlation matrix to a given
  square matrix in the Frobenius norm; `corr_thresholded` efficiently
  constructs a hard-thresholded correlation matrix using sparse matrix
  operations.

* New `dot_plot` in graphics: A dotplot is a way to visualize a small dataset
  in a way that immediately conveys the identity of every point in the plot.
  Dotplots are commonly seen in meta-analyses, where they are known
  as "forest plots", but can be used in many other settings as well.
  Most tables that appear in research papers can be represented
  graphically as a dotplot.
* statsmodels has added custom warnings to ``statsmodels.tools.sm_exceptions``. By default all of these warnings will be raised whenever appropriate. Use ``warnings.simplefilter`` to turn them off, if desired.
* Allow control over the namespace used to evaluate formulas with patsy via the ``eval_env`` keyword argument. See the :ref:`patsy-namespaces` documentation for more information.


Major Bugs fixed
----------------

* NA-handling with formulas is now correctly handled. :issue:`805`, :issue:`1877`.
* Better error messages when an array with an object dtype is used. :issue:`2013`.
* ARIMA forecasts were hard-coded for order of integration with ``d = 1``. :issue:`1562`.

.. currentmodule:: statsmodels.tsa

Backwards incompatible changes and deprecations
-----------------------------------------------

* RegressionResults.norm_resid is now a readonly property, rather than a function.
* The function ``statsmodels.tsa.filters.arfilter`` has been removed. This did not compute a recursive AR filter but was instead a convolution filter. Two new functions have been added with clearer names :func:`sm.tsa.filters.recursive_filter <tsa.filters.filtertools.recursive_filter>` and :func:`sm.tsa.filters.convolution_filter <tsa.filters.filtertools.convolution_filter>`.

Development summary and credits
-------------------------------

The previous version (0.5.0) was released August 14, 2014. Since then we have closed a total of 528 issues, 276 pull requests, and 252 regular issues. Refer to the :ref:`detailed list<issues_list_06>` for more information.

This release is a result of the work of the following 37 authors who contributed a total of 1531 commits. If for any reason we have failed to list your name in the below, please contact us:

A blurb about the number of changes and the contributors list.

* Alex Griffing <argriffi-at-ncsu.edu>
* Alex Parij <paris.alex-at-gmail.com>
* Ana Martinez Pardo <anamartinezpardo-at-gmail.com>
* Andrew Clegg <andrewclegg-at-users.noreply.github.com>
* Ben Duffield <bduffield-at-palantir.com>
* Chad Fulton <chad-at-chadfulton.com>
* Chris Kerr <cjk34-at-cam.ac.uk>
* Eric Chiang <eric.chiang.m-at-gmail.com>
* Evgeni Burovski <evgeni-at-burovski.me>
* gliptak <gliptak-at-gmail.com>
* Hans-Martin von Gaudecker <hmgaudecker-at-uni-bonn.de>
* Jan Schulz <jasc-at-gmx.net>
* jfoo <jcjf1983-at-gmail.com>
* Joe Hand <joe.a.hand-at-gmail.com>
* Josef Perktold <josef.pktd-at-gmail.com>
* jsphon <jonathanhon-at-hotmail.com>
* Justin Grana <jg3705a-at-student.american.edu>
* Kerby Shedden <kshedden-at-umich.edu>
* Kevin Sheppard <kevin.sheppard-at-economics.ox.ac.uk>
* Kyle Beauchamp <kyleabeauchamp-at-gmail.com>
* Lars Buitinck <l.buitinck-at-esciencecenter.nl>
* Max Linke <max_linke-at-gmx.de>
* Miroslav Batchkarov <mbatchkarov-at-gmail.com>
* m <mngu2382-at-gmail.com>
* Padarn Wilson <padarn-at-gmail.com>
* Paul Hobson <pmhobson-at-gmail.com>
* Pietro Battiston <me-at-pietrobattiston.it>
* Radim Řehůřek <radimrehurek-at-seznam.cz>
* Ralf Gommers <ralf.gommers-at-googlemail.com>
* Richard T. Guy <richardtguy84-at-gmail.com>
* Roy Hyunjin Han <rhh-at-crosscompute.com>
* Skipper Seabold <jsseabold-at-gmail.com>
* Tom Augspurger <thomas-augspurger-at-uiowa.edu>
* Trent Hauck <trent.hauck-at-gmail.com>
* Valentin Haenel <valentin.haenel-at-gmx.de>
* Vincent Arel-Bundock <varel-at-umich.edu>
* Yaroslav Halchenko <debian-at-onerussian.com>

.. note::

   Obtained by running ``git log v0.5.0..HEAD --format='* %aN <%aE>' | sed 's/@/\-at\-/' | sed 's/<>//' | sort -u``.

.. _issues_list_06:

Issues closed in the 0.6.0 development cycle
============================================

Issues closed in 0.6.0
----------------------

GitHub stats for 2013/08/14 - 2014/10/15 (tag: v0.5.0)

We closed a total of 528 issues, 276 pull requests and 252 regular issues;
this is the full list (generated with the script :file:`tools/github_stats.py`):

This list is automatically generated and may be incomplete.

Pull Requests (276):

* :pr:`2044`: ENH: Allow unit interval for binary models. Closes #2040.
* :pr:`1426`: ENH: Import arima_process stuff into tsa.api
* :pr:`2042`: Fix two minor typos in contrast.py
* :pr:`2034`: ENH: Handle missing for extra data with formulas
* :pr:`2035`: MAINT: Remove deprecated code for 0.6
* :pr:`1325`: ENH: add the Edgeworth expansion based on the normal distribution
* :pr:`2032`: DOC: What it is what it is.
* :pr:`2031`: ENH: Expose patsy eval_env to users.
* :pr:`2028`: ENH: Fix numerical issues in links and families.
* :pr:`2029`: DOC: Fix versions to match other docs.
* :pr:`1647`: ENH: Warn on non-convergence.
* :pr:`2014`: BUG: Fix forecasting for ARIMA with d == 2
* :pr:`2013`: ENH: Better error message on object dtype
* :pr:`2012`: BUG: 2d 1 columns -> 1d. Closes #322.
* :pr:`2009`: DOC: Update after refactor. Use code block.
* :pr:`2008`: ENH: Add wrapper for MixedLM
* :pr:`1954`: ENH: PHReg formula improvements
* :pr:`2007`: BLD: Fix build issues
* :pr:`2006`: BLD: Do not generate cython on clean. Closes #1852.
* :pr:`2000`: BLD: Let pip/setuptools handle dependencies that are not installed at all.
* :pr:`1999`: Gee offset exposure 1994 rebased
* :pr:`1998`: BUG/ENH Lasso emptymodel rebased
* :pr:`1989`: BUG/ENH: WLS generic robust cov_type did not use whitened,
* :pr:`1587`: ENH: Wrap X12/X13-ARIMA AUTOMDL. Closes #442.
* :pr:`1563`: ENH: Add plot_predict method to ARIMA models.
* :pr:`1995`: BUG: Fix issue #1993
* :pr:`1981`: ENH: Add api for covstruct. Clear __init__. Closes #1917.
* :pr:`1996`: DEV: Ignore .venv file.
* :pr:`1982`: REF: Rename jac -> score_obs. Closes #1785.
* :pr:`1987`: BUG tsa pacf, base bootstrap
* :pr:`1986`: Bug multicomp 1927 rebased
* :pr:`1984`: Docs add gee.rst
* :pr:`1985`: Bug uncentered latex table 1929 rebased
* :pr:`1983`: BUG: Fix compat asunicode
* :pr:`1574`: DOC: Fix math.
* :pr:`1980`: DOC: Documentation fixes
* :pr:`1974`: REF/Doc beanplot change default color, add notebook
* :pr:`1978`: ENH: Check input to binary models
* :pr:`1979`: BUG: Typo
* :pr:`1976`: ENH: Add _repr_html_ to SimpleTable
* :pr:`1977`: BUG: Fix import refactor victim.
* :pr:`1975`: BUG: Yule walker cast to float
* :pr:`1973`: REF: Move and expose webuse
* :pr:`1972`: TST: Add testing against NumPy 1.9 and matplotlib 1.4
* :pr:`1939`: ENH: Binstar build files
* :pr:`1952`: REF/DOC: Misc
* :pr:`1940`: REF: refactor and speedup of mixed LME
* :pr:`1937`: ENH: Quick access to online documentation
* :pr:`1942`: DOC: Rename Change README type to rst
* :pr:`1938`: ENH: Enable Python 3.4 testing
* :pr:`1924`: Bug gee cov type 1906 rebased
* :pr:`1870`: robust covariance, cov_type in fit
* :pr:`1859`: BUG: Do not use negative indexing with k_ar == 0. Closes #1858.
* :pr:`1914`: BUG: LikelihoodModelResults.pvalues use df_resid_inference
* :pr:`1899`: TST: fix assert_equal for pandas index
* :pr:`1895`: Bug multicomp pandas
* :pr:`1894`: BUG fix more ix indexing cases for pandas compat
* :pr:`1889`: BUG: fix ytick positions closes #1561
* :pr:`1887`: Bug pandas compat asserts
* :pr:`1888`: TST test_corrpsd Test_Factor: add noise to data
* :pr:`1886`: BUG pandas 0.15 compatibility in grouputils labels
* :pr:`1885`: TST: corr_nearest_factor, more informative tests
* :pr:`1884`: Fix: Add compat code for pd.Categorical in pandas>=0.15
* :pr:`1883`: BUG: add _ctor_param to TransfGen distributions
* :pr:`1872`: TST: fix _infer_freq for pandas .14+ compat
* :pr:`1867`: Ref covtype fit
* :pr:`1865`: Disable tst distribution 1864
* :pr:`1856`: _spg_optim returns history of objective function values
* :pr:`1854`: BLD: Do not hard-code path for building notebooks. Closes #1249
* :pr:`1851`: MAINT: Cor nearest factor tests
* :pr:`1847`: Newton regularize
* :pr:`1623`: BUG Negbin fit regularized
* :pr:`1797`: BUG/ENH: fix and improve constant detection
* :pr:`1770`: TST: anova with `-1` noconstant, add tests
* :pr:`1837`: Allow group variable to be passed as variable name when using formula
* :pr:`1839`: BUG: GEE score
* :pr:`1830`: BUG/ENH Use t
* :pr:`1832`: TST error with scipy 0.14 location distribution class
* :pr:`1827`: fit_regularized for linear models   rebase 1674
* :pr:`1825`: Phreg 1312 rebased
* :pr:`1826`: Lme api docs
* :pr:`1824`: Lme profile 1695 rebased
* :pr:`1823`: Gee cat subclass 1694 rebase
* :pr:`1781`: ENH: Glm add score_obs
* :pr:`1821`: Glm maint #1734 rebased
* :pr:`1820`: BUG: revert change to conf_int in PR #1819
* :pr:`1819`: Docwork
* :pr:`1772`: REF: cov_params allow case of only cov_params_default is defined
* :pr:`1771`: REF numpy >1.9 compatibility, indexing into empty slice closes #1754
* :pr:`1769`: Fix ttest 1d
* :pr:`1766`: TST: TestProbitCG increase bound for fcalls closes #1690
* :pr:`1709`: BLD: Made build extensions more flexible
* :pr:`1714`: WIP: fit_constrained
* :pr:`1706`: REF: Use fixed params in test. Closes #910.
* :pr:`1701`: BUG: Fix faulty logic. Do not raise when missing='raise' and no missing data.
* :pr:`1699`: TST/ENH StandardizeTransform, reparameterize TestProbitCG
* :pr:`1697`: Fix for statsmodels/statsmodels#1689
* :pr:`1692`: OSL Example: redundant cell in example removed
* :pr:`1688`: Kshedden mixed rebased of #1398
* :pr:`1629`: Pull request to fix bandwidth bug in issue 597
* :pr:`1666`: Include pyx in sdist but do not install
* :pr:`1683`: TST: GLM shorten random seed closes #1682
* :pr:`1681`: Dotplot kshedden rebased of 1294
* :pr:`1679`: BUG: Fix problems with predict handling offset and exposure
* :pr:`1677`: Update docstring of RegressionModel.predict()
* :pr:`1635`: Allow offset and exposure to be used together with log link; raise except...
* :pr:`1676`: Tests for SVAR
* :pr:`1671`: ENH: avoid hard-listed bandwidths -- use present dictionary (+typos fixed)
* :pr:`1643`: Allow matrix structure in covariance matrices to be exploited
* :pr:`1657`: BUG: Fix refactor victim.
* :pr:`1630`: DOC: typo, "intercept"
* :pr:`1619`: MAINT: Dataset docs cleanup and automatic build of docs
* :pr:`1612`: BUG/ENH Fix negbin exposure #1611
* :pr:`1610`: BUG/ENH fix llnull, extra kwds to recreate model
* :pr:`1582`: BUG: wls_prediction_std fix weight handling, see 987
* :pr:`1613`: BUG: Fix proportions allpairs #1493
* :pr:`1607`: TST: adjust precision, CI Debian, Ubuntu testing
* :pr:`1603`: ENH: Allow start_params in GLM
* :pr:`1600`: CLN: Regression plots fixes
* :pr:`1592`: DOC: Additions and fixes
* :pr:`1520`: CLN: Refactored so that there is no longer a need for 2to3
* :pr:`1585`: Cor nearest 1384 rebased
* :pr:`1553`: Gee maint 1528 rebased
* :pr:`1583`: BUG: For ARMA(0,0) ensure 1d bse and fix summary.
* :pr:`1580`: DOC: Fix links. [skip ci]
* :pr:`1572`: DOC: Fix link title [skip ci]
* :pr:`1566`: BLD: Fix copy paste path error for >= 3.3 Windows builds
* :pr:`1524`: ENH: Optimize Cython code. Use scipy blas function pointers.
* :pr:`1560`: ENH: Allow ARMA(0,0) in order selection
* :pr:`1559`: MAINT: Recover lost commits from vbench PR
* :pr:`1554`: Silenced test output introduced in medcouple
* :pr:`1234`: ENH: Robust skewness, kurtosis and medcouple measures
* :pr:`1484`: ENH: Add naive seasonal decomposition function
* :pr:`1551`: COMPAT: Fix failing test on Python 2.6
* :pr:`1472`: ENH: using human-readable group names instead of integer ids in MultiComparison
* :pr:`1437`: ENH: accept non-int definitions of cluster groups
* :pr:`1550`: Fix test gmm poisson
* :pr:`1549`: TST: Fix locally failing tests.
* :pr:`1121`: WIP: Refactor optimization code.
* :pr:`1547`: COMPAT: Correct bit_length for 2.6
* :pr:`1545`: MAINT: Fix missed usage of deprecated tools.rank
* :pr:`1196`: REF: ensure O(N log N) when using fft for acf
* :pr:`1154`: DOC: Add links for build machines.
* :pr:`1546`: DOC: Fix link to wrong notebook
* :pr:`1383`: MAINT: Deprecate rank in favor of np.linalg.matrix_rank
* :pr:`1432`: COMPAT: Add NumpyVersion from scipy
* :pr:`1438`: ENH: Option to avoid "center" environment.
* :pr:`1544`: BUG: Travis miniconda
* :pr:`1510`: CLN: Improve warnings to avoid generic warnings messages
* :pr:`1543`: TST: Suppress RuntimeWarning for L-BFGS-B
* :pr:`1507`: CLN: Silence test output
* :pr:`1540`: BUG: Correct derivative for exponential transform.
* :pr:`1536`: BUG: Restores coveralls for a single build
* :pr:`1535`: BUG: Fixes for 2.6 test failures, replacing astype(str) with apply(str)
* :pr:`1523`: Travis miniconda
* :pr:`1533`: DOC: Fix link to code on github
* :pr:`1531`: DOC: Fix stale links with linkcheck
* :pr:`1530`: DOC: Fix link
* :pr:`1527`: DOCS: Update docs add FAQ page
* :pr:`1525`: DOC: Update with Python 3.4 build notes
* :pr:`1518`: DOC: Ask for release notes and example.
* :pr:`1516`: DOC: Update examples contributing docs for current practice.
* :pr:`1517`: DOC: Be clear about data attribute of Datasets
* :pr:`1515`: DOC: Fix broken link
* :pr:`1514`: DOC: Fix formula import convention.
* :pr:`1506`: BUG: Format and decode errors in Python 2.6
* :pr:`1505`: TST: Test co2 load_data for Python 3.
* :pr:`1504`: BLD: New R versions require NAMESPACE file. Closes #1497.
* :pr:`1483`: ENH: Some utility functions for working with dates
* :pr:`1482`: REF: Prefer filters.api to __init__
* :pr:`1481`: ENH: Add weekly co2 dataset
* :pr:`1474`: DOC: Add plots for standard filter methods.
* :pr:`1471`: DOC: Fix import
* :pr:`1470`: DOC/BLD: Log code exceptions from nbgenerate
* :pr:`1469`: DOC: Fix bad links
* :pr:`1468`: MAINT: CSS fixes
* :pr:`1463`: DOC: Remove defunct argument. Change default kw. Closes #1462.
* :pr:`1452`: STY: import pandas as pd
* :pr:`1458`: BUG/BLD: exclude sandbox in relative path, not absolute
* :pr:`1447`: DOC: Only build and upload docs if we need to.
* :pr:`1445`: DOCS: Example landing page
* :pr:`1436`: DOC: Fix auto doc builds.
* :pr:`1431`: DOC: Add default for getenv. Fix paths. Add print_info
* :pr:`1429`: MAINT: Use ip_directive shipped with IPython
* :pr:`1427`: TST: Make tests fit quietly
* :pr:`1424`: ENH: Consistent results for transform_slices
* :pr:`1421`: ENH: Add grouping utilities code
* :pr:`1419`: Gee 1314 rebased
* :pr:`1414`: TST temporarily rename tests probplot other to skip them
* :pr:`1403`: Bug norm expan shapes
* :pr:`1417`: REF: Let subclasses keep kwds attached to data.
* :pr:`1416`: ENH: Make handle_data overwritable by subclasses.
* :pr:`1410`: ENH: Handle missing is none
* :pr:`1402`: REF: Expose missing data handling as classmethod
* :pr:`1387`: MAINT: Fix failing tests
* :pr:`1406`: MAINT: Tools improvements
* :pr:`1404`: Tst fix genmod link tests
* :pr:`1396`: REF: Multipletests reduce memory usage
* :pr:`1380`: DOC :Update vector_ar.rst
* :pr:`1381`: BLD: Do not check dependencies on egg_info for pip. Closes #1267.
* :pr:`1302`: BUG: Fix typo.
* :pr:`1375`: STY: Remove unused imports and comment out unused libraries in setup.py
* :pr:`1143`: DOC: Update backport notes for new workflow.
* :pr:`1374`: ENH: Import tsaplots into tsa namespace. Closes #1359.
* :pr:`1369`: STY: Pep-8 cleanup
* :pr:`1370`: ENH: Support ARMA(0,0) models.
* :pr:`1368`: STY: Pep 8 cleanup
* :pr:`1367`: ENH: Make sure mle returns attach to results.
* :pr:`1365`: STY: Import and pep 8 cleanup
* :pr:`1364`: ENH: Get rid of hard-coded lbfgs. Closes #988.
* :pr:`1363`: BUG: Fix typo.
* :pr:`1361`: ENH: Attach mlefit to results not model.
* :pr:`1360`: ENH: Import adfuller into tsa namespace
* :pr:`1346`: STY: PEP-8 Cleanup
* :pr:`1344`: BUG: Use missing keyword given to ARMA.
* :pr:`1340`: ENH: Protect against ARMA convergence failures.
* :pr:`1334`: ENH: ARMA order select convenience function
* :pr:`1339`: Fix typos
* :pr:`1336`: REF: Get rid of plain assert.
* :pr:`1333`: STY: __all__ should be after imports.
* :pr:`1332`: ENH: Add Bunch object to tools.
* :pr:`1331`: ENH: Always use unicode.
* :pr:`1329`: BUG: Decode metadata to utf-8. Closes #1326.
* :pr:`1330`: DOC: Fix typo. Closes #1327.
* :pr:`1185`: Added support for pandas when pandas was installed directly from git trunk
* :pr:`1315`: MAINT: Change back to path for build box
* :pr:`1305`: TST: Update hard-coded path.
* :pr:`1290`: ENH: Add seasonal plotting.
* :pr:`1296`: BUG/TST: Fix ARMA forecast when start == len(endog). Closes #1295
* :pr:`1292`: DOC: cleanup examples folder and webpage
* :pr:`1286`: Make sure PeriodIndex passes through tsa. Closes #1285.
* :pr:`1271`: Silverman enhancement - Issue #1243
* :pr:`1264`: Doc work GEE, GMM, sphinx warnings
* :pr:`1179`: REF/TST: `ProbPlot` now uses `resettable_cache` and added some kwargs to plotting fxns
* :pr:`1225`: Sandwich mle
* :pr:`1258`: Gmm new rebased
* :pr:`1255`: ENH add GEE to genmod
* :pr:`1254`: REF: Results.predict convert to array and adjust shape
* :pr:`1192`: TST: enable tests for llf after change to WLS.loglike see #1170
* :pr:`1253`: Wls llf fix
* :pr:`1233`: sandbox kernels bugs uniform kernel and confint
* :pr:`1240`: Kde weights 1103 823
* :pr:`1228`: Add default value tags to adfuller() docs
* :pr:`1198`: fix typo
* :pr:`1230`: BUG: numerical precision in resid_pearson with perfect fit #1229
* :pr:`1214`: Compare lr test rebased
* :pr:`1200`: BLD: do not install \*.pyx \*.c  MANIFEST.in
* :pr:`1202`: MAINT: Sort backports to make applying easier.
* :pr:`1157`: Tst precision
* :pr:`1161`: add a fitting interface for simultaneous log likelihood and score, for lbfgs, tested with MNLogit
* :pr:`1160`: DOC: update scipy version from 0.7 to 0.9.0
* :pr:`1147`: ENH: add lbfgs for fitting
* :pr:`1156`: ENH: Raise on 0,0 order models in AR(I)MA. Closes #1123
* :pr:`1149`: BUG: Fix small data issues for ARIMA.
* :pr:`1092`: Fixed duplicate svd in RegressionModel
* :pr:`1139`: TST: Silence tests
* :pr:`1135`: Misc style
* :pr:`1088`: ENH: add predict_prob to poisson
* :pr:`1125`: REF/BUG: Some GLM cleanup. Used trimmed results in NegativeBinomial variance.
* :pr:`1124`: BUG: Fix ARIMA prediction when fit without a trend.
* :pr:`1118`: DOC: Update gettingstarted.rst
* :pr:`1117`: Update ex_arma2.py
* :pr:`1107`: REF: Deprecate stand_mad. Add center keyword to mad. Closes #658.
* :pr:`1089`: ENH: exp(poisson.logpmf()) for poisson better behaved.
* :pr:`1077`: BUG: Allow 1d exog in ARMAX forecasting.
* :pr:`1075`: BLD: Fix build issue on some versions of easy_install.
* :pr:`1071`: Update setup.py to fix broken install on OSX
* :pr:`1052`: DOC: Updating contributing docs
* :pr:`1136`: RLS: Add IPython tools for easier backporting of issues.
* :pr:`1091`: DOC: minor git typo
* :pr:`1082`: coveralls support
* :pr:`1072`: notebook examples title cell
* :pr:`1056`: Example: reg diagnostics
* :pr:`1057`: COMPAT: Fix py3 caching for get_rdatasets.
* :pr:`1045`: DOC/BLD: Update from nbconvert to IPython 1.0.
* :pr:`1026`: DOC/BLD: Add LD_LIBRARY_PATH to env for docs build.

Issues (252):

* :issue:`2040`: enh: fractional Logit, Probit
* :issue:`1220`: missing in extra data (example sandwiches, robust covariances)
* :issue:`1877`: error with GEE on missing data.
* :issue:`805`: nan with categorical in formula
* :issue:`2036`: test in links require exact class so Logit cannot work in place of logit
* :issue:`2010`: Go over deprecations again for 0.6.
* :issue:`1303`: patsy library not automatically installed
* :issue:`2024`: genmod Links numerical improvements
* :issue:`2025`: GEE requires exact import for cov_struct
* :issue:`2017`: Matplotlib warning about too many figures
* :issue:`724`: check warnings
* :issue:`1562`: ARIMA forecasts are hard-coded for d=1
* :issue:`880`: DataFrame with bool type not cast correctly.
* :issue:`1992`: MixedLM style
* :issue:`322`: acf / pacf do not work on pandas objects
* :issue:`1317`: AssertionError: attr is not equal [dtype]: dtype('object') != dtype('datetime64[ns]')
* :issue:`1875`: dtype bug object arrays (raises in clustered standard errors code)
* :issue:`1842`: dtype object, glm.fit() gives AttributeError: sqrt
* :issue:`1300`: Doc errors, missing
* :issue:`1164`: RLM cov_params, t_test, f_test do not use bcov_scaled
* :issue:`1019`: 0.6.0 Roadmap
* :issue:`554`: Prediction Standard Errors
* :issue:`333`: ENH tools: squeeze in R export file
* :issue:`1990`: MixedLM does not have a wrapper
* :issue:`1897`: Consider depending on setuptools in setup.py
* :issue:`2003`: pip install now fails silently
* :issue:`1852`: do not cythonize when cleaning up
* :issue:`1991`: GEE formula interface does not take offset/exposure
* :issue:`442`: Wrap x-12 arima
* :issue:`1993`: MixedLM bug
* :issue:`1917`: API: GEE access to genmod.covariance_structure through api
* :issue:`1785`: REF: rename jac -> score_obs
* :issue:`1969`: pacf has incorrect standard errors for lag 0
* :issue:`1434`: A small bug in GenericLikelihoodModelResults.bootstrap()
* :issue:`1408`: BUG test failure with tsa_plots
* :issue:`1337`: DOC: HCCM are now available for WLS
* :issue:`546`: influence and outlier documentation
* :issue:`1532`: DOC: Related page is out of date
* :issue:`1386`: Add minimum matplotlib to docs
* :issue:`1068`: DOC: keeping documentation of old versions on sourceforge
* :issue:`329`: link to examples and datasets from module pages
* :issue:`1804`: PDF documentation for statsmodels
* :issue:`202`: Extend robust standard errors for WLS/GLS
* :issue:`1519`: Link to user-contributed examples in docs
* :issue:`1053`: inconvenient: logit when endog is (1,2) instead of (0,1)
* :issue:`1555`: SimpleTable: add repr html for ipython notebook
* :issue:`1366`: Change default start_params to .1 in ARMA
* :issue:`1869`: yule_walker (from `statsmodels.regression`) raises exception when given an integer array
* :issue:`1651`: statsmodels.tsa.ar_model.ARResults.predict
* :issue:`1738`: GLM robust sandwich covariance matrices
* :issue:`1779`: Some directories under statsmodels dont have __init_.py
* :issue:`1242`: No support for (0, 1, 0) ARIMA Models
* :issue:`1571`: expose webuse, use cache
* :issue:`1860`: ENH/BUG/DOC: Bean plot should allow for separate widths of bean and violins.
* :issue:`1831`: TestRegressionNM.test_ci_beta2 i386 AssertionError
* :issue:`1079`: bugfix release 0.5.1
* :issue:`1338`: Raise Warning for HCCM use in WLS/GLS
* :issue:`1430`: scipy min version / issue
* :issue:`276`: memoize, last argument wins, how to attach sandwich to Results?
* :issue:`1943`: REF/ENH: LikelihoodModel.fit optimization, make hessian optional
* :issue:`1957`: BUG: Re-create OLS model using _init_keys
* :issue:`1905`: Docs: online docs are missing GEE
* :issue:`1898`: add python 3.4 to continuous integration testing
* :issue:`1684`: BUG: GLM NegativeBinomial: llf ignores offset and exposure
* :issue:`1256`: REF: GEE handling of default covariance matrices
* :issue:`1760`: Changing covariance_type on results
* :issue:`1906`: BUG: GEE default covariance is not used
* :issue:`1931`: BUG: GEE subclasses NominalGEE do not work with pandas exog
* :issue:`1904`: GEE Results does not have a Wrapper
* :issue:`1918`: GEE: required attributes missing, df_resid
* :issue:`1919`: BUG GEE.predict uses link instead of link.inverse
* :issue:`1858`: BUG: arimax forecast should special case k_ar == 0
* :issue:`1903`: BUG: pvalues for cluster robust, with use_t do not use df_resid_inference
* :issue:`1243`: kde silverman bandwidth for non-gaussian kernels
* :issue:`1866`: Pip dependencies
* :issue:`1850`: TST test_corr_nearest_factor fails on Ubuntu
* :issue:`292`: python 3 examples
* :issue:`1868`: ImportError: No module named compat  [ from statsmodels.compat import lmap ]
* :issue:`1890`: BUG tukeyhsd nan in group labels
* :issue:`1891`: TST test_gmm outdated pandas, compat
* :issue:`1561`: BUG plot for tukeyhsd, MultipleComparison
* :issue:`1864`: test failure sandbox distribution transformation with scipy 0.14.0
* :issue:`576`: Add contributing guidelines
* :issue:`1873`: GenericLikelihoodModel is not picklable
* :issue:`1822`: TST failure on Ubuntu pandas 0.14.0 , problems with frequency
* :issue:`1249`: Source directory problem for notebook examples
* :issue:`1855`: anova_lm throws error on models created from api.ols but not formula.api.ols
* :issue:`1853`: a large number of hardcoded paths
* :issue:`1792`: R² adjusted strange after including interaction term
* :issue:`1794`: REF: has_constant, k_constant, include implicit constant detection in base
* :issue:`1454`: NegativeBinomial missing fit_regularized method
* :issue:`1615`: REF DRYing fit methods
* :issue:`1453`: Discrete NegativeBinomialModel regularized_fit ValueError: matrices are not aligned
* :issue:`1836`: BUG Got an TypeError trying to import statsmodels.api
* :issue:`1829`: BUG: GLM summary show "t"  use_t=True for summary
* :issue:`1828`: BUG summary2 does not propagate/use use_t
* :issue:`1812`: BUG/ REF conf_int and use_t
* :issue:`1835`: Problems with installation using easy_install
* :issue:`1801`: BUG 'f_gen' missing in scipy 0.14.0
* :issue:`1803`: Error revealed by numpy 1.9.0r1
* :issue:`1834`: stackloss
* :issue:`1728`: GLM.fit maxiter=0  incorrect
* :issue:`1795`: singular design with offset ?
* :issue:`1730`: ENH/Bug cov_params, generalize, avoid ValueError
* :issue:`1754`: BUG/REF: assignment to slices in numpy >= 1.9 (emplike)
* :issue:`1409`: GEE test errors on Debian Wheezy
* :issue:`1521`: ubuntu failures: tsa_plot and grouputils
* :issue:`1415`: test failure test_arima.test_small_data
* :issue:`1213`: df_diff in anova_lm
* :issue:`1323`: Contrast Results after t_test summary broken for 1 parameter
* :issue:`109`: TestProbitCG failure on Ubuntu
* :issue:`1690`: TestProbitCG: 8 failing tests (Python 3.4 / Ubuntu 12.04)
* :issue:`1763`: Johansen method does not give correct index values
* :issue:`1761`: doc build failures: ipython version ? ipython directive
* :issue:`1762`: Unable to build
* :issue:`1745`: UnicodeDecodeError raised by get_rdataset("Guerry", "HistData")
* :issue:`611`: test failure foreign with pandas 0.7.3
* :issue:`1700`: faulty logic in missing handling
* :issue:`1648`: ProbitCG failures
* :issue:`1689`: test_arima.test_small_data: SVD fails to converge (Python 3.4 / Ubuntu 12.04)
* :issue:`597`: BUG: nonparametric: kernel, efficient=True changes bw even if given
* :issue:`1606`: BUILD from sdist broken if cython available
* :issue:`1246`: test failure test_anova.TestAnova2.test_results
* :issue:`50`: t_test, f_test, model.py for normal instead of t-distribution
* :issue:`1655`: newey-west different than R?
* :issue:`1682`: TST test failure on Ubuntu, random.seed
* :issue:`1614`: docstring for regression.linear_model.RegressionModel.predict() does not match implementation
* :issue:`1318`: GEE and GLM scale parameter
* :issue:`519`: L1 fit_regularized cleanup, comments
* :issue:`651`: add structure to example page
* :issue:`1067`: Kalman Filter convergence. How close is close enough?
* :issue:`1281`: Newton convergence failure prints warnings instead of warning
* :issue:`1628`: Unable to install statsmodels in the same requirements file as numpy, pandas, etc.
* :issue:`617`: Problem in installing statsmodels in Fedora 17 64-bit
* :issue:`935`: ll_null in likelihoodmodels discrete
* :issue:`704`: datasets.sunspot: wrong link in description
* :issue:`1222`: NegativeBinomial ignores exposure
* :issue:`1611`: BUG NegativeBinomial ignores exposure and offset
* :issue:`1608`: BUG: NegativeBinomial, llnul is always default 'nb2'
* :issue:`1221`: llnull with exposure ?
* :issue:`1493`: statsmodels.stats.proportion.proportions_chisquare_allpairs has hardcoded value
* :issue:`1260`: GEE test failure on Debian
* :issue:`1261`: test failure on Debian
* :issue:`443`: GLM.fit does not allow start_params
* :issue:`1602`: Fitting GLM with a pre-assigned starting parameter
* :issue:`1601`: Fitting GLM with a pre-assigned starting parameter
* :issue:`890`: regression_plots problems (pylint) and missing test coverage
* :issue:`1598`: Is "old" string formatting Python 3 compatible?
* :issue:`1589`: AR vs ARMA order specification
* :issue:`1134`: Mark knownfails
* :issue:`1259`: Parameterless models
* :issue:`616`: python 2.6, python 3 in single codebase
* :issue:`1586`: Kalman Filter errors with new pyx
* :issue:`1565`: build_win_bdist*_py3*.bat are using the wrong compiler
* :issue:`843`: UnboundLocalError When trying to install OS X
* :issue:`713`: arima.fit performance
* :issue:`367`: unable to install on RHEL 5.6
* :issue:`1548`: testtransf error
* :issue:`1478`: is sm.tsa.filters.arfilter an AR filter?
* :issue:`1420`: GMM poisson test failures
* :issue:`1145`: test_multi noise
* :issue:`1539`: NegativeBinomial   strange results with bfgs
* :issue:`936`: vbench for statsmodels
* :issue:`1153`: Where are all our testing machines?
* :issue:`1500`: Use Miniconda for test builds
* :issue:`1526`: Out of date docs
* :issue:`1311`: BUG/BLD 3.4 compatibility of cython c files
* :issue:`1513`: build on osx -python-3.4
* :issue:`1497`: r2nparray needs NAMESPACE file
* :issue:`1502`: coveralls coverage report for files is broken
* :issue:`1501`: pandas in/out in predict
* :issue:`1494`: truncated violin plots
* :issue:`1443`: Crash from python.exe using linear regression of statsmodels
* :issue:`1462`: qqplot line kwarg is broken/docstring is wrong
* :issue:`1457`: BUG/BLD: Failed build if "sandbox" anywhere in statsmodels path
* :issue:`1441`: wls function: syntax error "unexpected EOF while parsing" occurs when name of dependent variable starts with digits
* :issue:`1428`: ipython_directive does not work with ipython
* :issue:`1385`: SimpleTable in Summary (e.g. OLS) is slow for large models
* :issue:`1399`: UnboundLocalError: local variable 'fittedvalues' referenced before assignment
* :issue:`1377`: TestAnova2.test_results fails with pandas 0.13.1
* :issue:`1394`: multipletests: reducing memory consumption
* :issue:`1267`: Packages cannot have both pandas and statsmodels in install_requires
* :issue:`1359`: move graphics.tsa to tsa.graphics
* :issue:`356`: docs take up a lot of space
* :issue:`988`: AR.fit no precision options for fmin_l_bfgs_b
* :issue:`990`: AR fit with bfgs: large score
* :issue:`14`: arma with exog
* :issue:`1348`: reset_index + set_index with drop=False
* :issue:`1343`: ARMA does not pass missing keyword up to TimeSeriesModel
* :issue:`1326`: formula example notebook broken
* :issue:`1327`: typo in docu-code for "Outlier and Influence Diagnostic Measures"
* :issue:`1309`: Box-Cox transform (some code needed: lambda estimator)
* :issue:`1059`: sm.tsa.ARMA making ma invertibility
* :issue:`1295`: Bug in ARIMA forecasting when start is int len(endog) and dates are given
* :issue:`1285`: tsa models fail on PeriodIndex with pandas
* :issue:`1269`: KPSS test for stationary processes
* :issue:`1268`: Feature request: Exponential smoothing
* :issue:`1250`: DOCs error in var_plots
* :issue:`1032`: Poisson predict breaks on list
* :issue:`347`: minimum number of observations - document or check ?
* :issue:`1170`: WLS log likelihood, aic and bic
* :issue:`1187`:  sm.tsa.acovf fails when both unbiased and fft are True
* :issue:`1239`: sandbox kernels, problems with inDomain
* :issue:`1231`: sandbox kernels confint missing alpha
* :issue:`1245`: kernels cosine differs from Stata
* :issue:`823`: KDEUnivariate with weights
* :issue:`1229`: precision problems in degenerate case
* :issue:`1219`: select_order
* :issue:`1206`: REF: RegressionResults cov-HCx into cached attributes
* :issue:`1152`: statsmodels failing tests with pandas
* :issue:`1195`: pyximport.install() before import api crash
* :issue:`1066`: gmm.IV2SLS has wrong predict signature
* :issue:`1186`: OLS when exog is 1d
* :issue:`1113`: TST: precision too high in test_normality
* :issue:`1159`: scipy version is still >= 0.7?
* :issue:`1108`: SyntaxError: unqualified exec is not allowed in function 'test_EvalEnvironment_capture_flag
* :issue:`1116`: Typo in Example Doc?
* :issue:`1123`: BUG : arima_model._get_predict_out_of_sample, ignores exogenous of there is no trend ?
* :issue:`1155`: ARIMA - The computed initial AR coefficients are not stationary
* :issue:`979`: Win64 binary cannot find Python installation
* :issue:`1046`: TST: test_arima_small_data_bug on current main
* :issue:`1146`: ARIMA fit failing for small set of data due to invalid maxlag
* :issue:`1081`: streamline linear algebra for linear model
* :issue:`1138`: BUG: pacf_yw does not demean
* :issue:`1127`: Allow linear link model with Binomial families
* :issue:`1122`: no data cleaning for statsmodels.genmod.families.varfuncs.NegativeBinomial()
* :issue:`658`: robust.mad is not being computed correctly or is non-standard definition; it returns the median
* :issue:`1076`: Some issues with ARMAX forecasting
* :issue:`1073`: easy_install sandbox violation
* :issue:`1115`: EasyInstall Problem
* :issue:`1106`: bug in robust.scale.mad?
* :issue:`1102`: Installation Problem
* :issue:`1084`: DataFrame.sort_index does not use ascending when then value is a list with a single element
* :issue:`393`: marginal effects in discrete choice do not have standard errors defined
* :issue:`1078`: Use pandas.version.short_version
* :issue:`96`: deepcopy breaks on ResettableCache
* :issue:`1055`: datasets.get_rdataset   string decode error on python 3
* :issue:`46`: tsa.stattools.acf confint needs checking and tests
* :issue:`957`: ARMA start estimate with numpy main
* :issue:`62`: GLSAR incorrect initial condition in whiten
* :issue:`1021`: from_formula() throws error - problem installing
* :issue:`911`: noise in stats.power tests
* :issue:`472`: Update roadmap for 0.5
* :issue:`238`: release 0.5
* :issue:`1006`: update nbconvert to IPython 1.0
* :issue:`1038`: DataFrame with integer names not handled in ARIMA
* :issue:`1036`: Series no longer inherits from ndarray
* :issue:`1028`: Test fail with windows and Anaconda - Low priority
* :issue:`676`: acorr_breush_godfrey  undefined nlags
* :issue:`922`: lowess returns inconsistent with option
* :issue:`425`: no bse in robust with norm=TrimmedMean
* :issue:`1025`: add_constant incorrectly detects constant column
:orphan:

==============
Release 0.11.0
==============

Release summary
===============

statsmodels is using github to store the updated documentation. Two version are available:

- `Stable <https://www.statsmodels.org/>`_, the latest release
- `Development <https://www.statsmodels.org/devel/>`_, the latest build of the main branch

**Warning**

API stability is not guaranteed for new features, although even in
this case changes will be made in a backwards compatible way if
possible. The stability of a new feature depends on how much time it
was already in statsmodels main and how much usage it has already
seen.  If there are specific known problems or limitations, then they
are mentioned in the docstrings.

Stats
-----
**Issues Closed**: 335

**Pull Requests Merged**: 301


The Highlights
==============

Regression
----------
Rolling OLS and WLS are implemented in :class:`~statsmodels.regression.rolling.RollingOLS`
and :class:`~statsmodels.regression.rolling.RollingWLS`. These function similarly to the estimators
recently removed from pandas.

Statistics
----------
Add the Oaxaca-Blinder decomposition (:class:`~statsmodels.stats.oaxaca.OaxacaBlinder`) that
decomposes the difference in group means into with and between components.

Add the Distance dependence measures statistics
(:func:`~statsmodels.stats.dist_dependence_measures.distance_statistics`) and the Distance Covariance
test (:func:`~statsmodels.stats.dist_dependence_measures.distance_covariance_test`).

Statespace Models
-----------------

Linear exponential smoothing models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Model class: :class:`~statsmodels.tsa.statespace.ExponentialSmoothing`
- Alternative to :class:`~statsmodels.tsa.statespace.exponential_smoothing.ExponentialSmoothing`
- Only supports linear models
- Is part of the class of state space models and so inherits some additional
  functionality.

Methods to apply parameters fitted on one dataset to another dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Methods: ``extend``, ``append``, and ``apply``, for state space results classes
- These methods allow applying fitted parameters from a training dataset to a
  test dataset in various ways
- Useful for conveniently performing cross-validation exercises

Method to hold some parameters fixed at known values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Methods: :func:`statsmodels.tsa.statespace.mlemodel.MLEModel.fix_params` and
  :func:`statsmodels.tsa.statespace.mlemodel.MLEModel.fit_constrained` for state
  space model classes.
- These methods allow setting some parameters to known values and then
  estimating the remaining parameters

Option for low memory operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Argument: ``low_memory=True``, for ``fit``, ``filter``, and ``smooth``
- Only a subset of results are available when using this option, but it does
  allow for prediction, diagnostics, and forecasting
- Useful to speed up cross-validation exercises

Improved access to state estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Attribute: ``states``, for state space results classes
- Wraps estimated states (predicted, filtered, smoothed) as Pandas objects with
  the appropriate index.
- Adds human-readable names for states, where possible.

Improved simulation and impulse responses for time-varying models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Argument: ``anchor`` allows specifying the period after which to begin the simulation.
- Example: to simulate data following the sample period, use ``anchor='end'``

Time-Series Analysis
--------------------

STL Decomposition
~~~~~~~~~~~~~~~~~
- Class implementing the STL decomposition :class:`~statsmodels.tsa.seasonal.STL`.

New AR model
~~~~~~~~~~~~

- Model class: :class:`~statsmodels.tsa.ar_model.AutoReg`
- Estimates parameters using conditional MLE (OLS)
- Adds the ability to specify exogenous variables, include time trends,
  and add seasonal dummies.
- The function :class:`~statsmodels.tsa.ar_model.ar_select_order` performs lag length selection
  for AutoReg models.

New ARIMA model
~~~~~~~~~~~~~~~

- Model class: :class:`~statsmodels.tsa.arima.model.ARIMA`
- Incorporates a variety of SARIMA estimators
    - MLE via state space methods (SARIMA models)
    - MLE via innovations algorithm (SARIMA models)
    - Hannan-Rissanen (ARIMA models)
    - Burg's method (AR models)
    - Innovations algorithm (MA models)
    - Yule-Walker (AR models)
- Handles exogenous regressors via GLS or by MLE with state space methods.
- Is part of the class of state space models and so inherits some additional
  functionality.

Zivot-Andrews Test
~~~~~~~~~~~~~~~~~~
The Zivot-Andrews test for unit roots in the presence of structural breaks has
been added in :func:`~statsmodels.tsa.stattools.zivot_andrews`.

More robust regime switching models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Implementation of the Hamilton filter and Kim smoother in log space avoids
  underflow errors.


What's new - an overview
========================

The following lists the main new features of statsmodels 0.11. In addition,
release 0.11 includes bug fixes, refactorings and improvements in many areas.

Major Feature
-------------
- Allow fixing parameters in state space models  (:pr:`5735`)
- Add new version of ARIMA-type estimators (AR, ARIMA, SARIMAX)  (:pr:`5827`)
- Add STL decomposition for time series  (:pr:`5926`)
- Functional SIR  (:pr:`5963`)
- Zivot Andrews test  (:pr:`6014`)
- Added Oaxaca-Blinder Decomposition  (:pr:`6026`)
- Add rolling WLS and OLS  (:pr:`6028`)
- Replacement for AR  (:pr:`6087`)

Performance Improvements
------------------------
- Cythonize innovations algo and filter  (:pr:`5947`)
- Only perform required predict iterations in state space models  (:pr:`6064`)
- State space: Improve low memory usability; allow in fit, loglike  (:pr:`6071`)

Submodules
----------

``base``
~~~~~~~~
- Clarify xname length and purpose  (:pr:`5957`)
- Remove unnecessary pickle use  (:pr:`6091`)
- Fix accepting of eval environment for formula  (:pr:`6152`)
- Workaround NumPy ptp issue  (:pr:`6316`)


``discrete``
~~~~~~~~~~~~
- Test_constrained  (:pr:`5821`)
- Improve the cvxopt not found error  (:pr:`6163`)


``genmod``
~~~~~~~~~~
- Improvements to BayesMixedGLM docs, argument checking  (:pr:`5895`)
- Scale parameter handling in GEE  (:pr:`6208`)
- Add example notebook for GEE score tests  (:pr:`6299`)
- Fix bug in ridge for vector alpha  (:pr:`6442`)

``graphics``
~~~~~~~~~~~~
- Plot only unique censored points  (:pr:`6124`)
- Add missing keyword argument to plot_acf  (:pr:`6227`)
- And vlines option to plot_fit  (:pr:`6266`)
- Pass arguments through in plot_leverage_resid2  (:pr:`6281`)


``io``
~~~~~~
- Clarify summary2 documentation  (:pr:`6118`)


``nonparametric``
~~~~~~~~~~~~~~~~~
- Ensure BW is not 0  (:pr:`6292`)
- Check dtype in KDEUnivariate  (:pr:`6314`)
- Supporting custom kernel in local linear kernel regression  (:pr:`6375`)



``regression``
~~~~~~~~~~~~~~
- Test for anova_nistcertified  (:pr:`5797`)
- Remove no-longer-needed HC_se lookups  (:pr:`5841`)
- Dimension reduction for covariance matrices  (:pr:`5852`)
- Use class to define MixedLM variance components structure  (:pr:`5898`)
- Add rolling WLS and OLS  (:pr:`6028`)
- Prepare for Rolling Least Squares  (:pr:`6056`)
- Improve regression doc strings  (:pr:`6077`)
- Fix summary table header for mixedlm  (:pr:`6217`)


``robust``
~~~~~~~~~~
- Robust  (:pr:`5819`)
- Make mad function behave correctly when used on empty inputs  (:pr:`5968`)


``stats``
~~~~~~~~~
- Lilliefors min nobs not set  (:pr:`5610`)
- Replace alpha=0.05 with alpha=alpha  (:pr:`5998`)
- Added Oaxaca-Blinder Decomposition  (:pr:`6026`)
- Improve Ljung-Box  (:pr:`6079`)
- Correct thresholding in correlation tools  (:pr:`6105`)
- Use self.data consistently  (:pr:`6144`)
- Better argument checking for StratifiedTable  (:pr:`6294`)
- Restore multicomp  (:pr:`6320`)
- Improve Ljung Box diagnostics  (:pr:`6324`)
- Correct standardization in robust skewness  (:pr:`6374`)
- Distance dependence measures  (:pr:`6401`)
- Improve diagnostics   (:pr:`6410`)



``tools``
~~~~~~~~~
- Fix error introduced in isestimable  (:pr:`6081`)
- Fix axis in irq  (:pr:`6391`)



``tsa``
~~~~~~~
- Use cython fused types to simplify statespace code  (:pr:`5283`)
- Allow fixing parameters in state space models  (:pr:`5735`)
- Markov switching in log space: Hamilton filter / Kim smoother  (:pr:`5826`)
- Add new version of ARIMA-type estimators (AR, ARIMA, SARIMAX)  (:pr:`5827`)
- Exponential smoothing - damped trend gives incorrect param, predictions  (:pr:`5893`)
- State space: add methods to apply fitted parameters to new observations or new dataset  (:pr:`5915`)
- TVTP for Markov regression  (:pr:`5917`)
- Add STL decomposition for time series  (:pr:`5926`)
- Cythonize innovations algo and filter  (:pr:`5947`)
- Zivot Andrews test  (:pr:`6014`)
- Improve ARMA startparams  (:pr:`6018`)
- Fix ARMA so that it works with exog when trend=nc  (:pr:`6070`)
- Clean tsatools docs  (:pr:`6075`)
- Improve Ljung-Box  (:pr:`6079`)
- Replacement for AR  (:pr:`6087`)
- Incorrect TSA index if loc resolves to slice  (:pr:`6130`)
- Division by zero in exponential smoothing if damping_slope=0  (:pr:`6232`)
- Forecasts now ignore non-monotonic period index  (:pr:`6242`)
- Hannan-Rissanen third stage is invalid if non-stationary/invertible  (:pr:`6258`)
- Fix notebook  (:pr:`6279`)
- Correct VAR summary when model contains exog variables  (:pr:`6286`)
- Fix conf interval with MI  (:pr:`6297`)
- Ensure inputs are finite in granger causality test  (:pr:`6318`)
- Fix trend due to recent changes  (:pr:`6321`)
- Improve Ljung Box diagnostics  (:pr:`6324`)
- Documentation for release v0.11  (:pr:`6338`)
- Fix _get_index_loc with date strings  (:pr:`6340`)
- Use correct exog names  (:pr:`6389`)



``tsa.statespace``
~~~~~~~~~~~~~~~~~~
- Use cython fused types to simplify statespace code  (:pr:`5283`)
- Allow fixing parameters in state space models  (:pr:`5735`)
- Add new version of ARIMA-type estimators (AR, ARIMA, SARIMAX)  (:pr:`5827`)
- MLEModel now passes nobs to Representation  (:pr:`6050`)
- Only perform required predict iterations in state space models  (:pr:`6064`)
- State space: Improve low memory usability; allow in fit, loglike  (:pr:`6071`)
- State space: cov_params computation in fix_params context  (:pr:`6072`)
- Add conserve memory tests.  (:pr:`6073`)
- Improve cov_params in append, extend, apply  (:pr:`6074`)
- Seasonality in SARIMAX Notebook  (:pr:`6096`)
- Improve SARIMAX start_params if too few nobs  (:pr:`6102`)
- Fix score computation with fixed params  (:pr:`6104`)
- Add exact diffuse initialization as an option for SARIMAX, UnobservedComponents  (:pr:`6111`)
- Compute standardized forecast error in diffuse period if possible  (:pr:`6131`)
- Start_params for VMA model with exog.  (:pr:`6133`)
- Adds state space version of linear exponential smoothing models  (:pr:`6179`)
- State space: add wrapped states and, where possible, named states  (:pr:`6181`)
- Allow dynamic factor starting parameters computation with NaNs values  (:pr:`6231`)
- Dynamic factor model use AR model for error start params if error_var=False  (:pr:`6233`)
- SARIMAX index behavior with simple_differencing=True  (:pr:`6239`)
- Parameter names in DynamicFactor for unstructured error covariance matrix  (:pr:`6240`)
- SARIMAX: basic validation for order, seasonal_order  (:pr:`6241`)
- Update SARIMAX to use SARIMAXSpecification for more consistent input handling  (:pr:`6250`)
- State space: Add finer-grained memory conserve settings  (:pr:`6254`)
- Cloning of arima.ARIMA models.  (:pr:`6260`)
- State space: saving fixed params w/ extend, apply, append  (:pr:`6261`)
- State space: Improve simulate, IRF, prediction  (:pr:`6280`)
- State space: deprecate out-of-sample w/ unsupported index  (:pr:`6332`)
- State space: integer params can cause imaginary output  (:pr:`6333`)
- Append, extend check that index matches model  (:pr:`6334`)
- Fix k_exog, k_trend in arima.ARIMA; raise error when cloning a model with exog if no new exog given  (:pr:`6337`)
- Documentation for release v0.11  (:pr:`6338`)
- RecursiveLS should not allow `fix_params` method.  (:pr:`6415`)
- More descriptive error message for recursive least squares parameter constraints.  (:pr:`6424`)
- Diffuse multivariate case w/ non-diagonal observation innovation covariance matrix  (:pr:`6425`)



``tsa.vector.ar``
~~~~~~~~~~~~~~~~~
- Raise in GC test for VAR(0)  (:pr:`6285`)
- Correct VAR summary when model contains exog variables  (:pr:`6286`)
- Use correct exog names  (:pr:`6389`)


Build
-----
- Ignore warns on 32 bit linux  (:pr:`6005`)
- Travis CI: The sudo: tag is deprecated in Travis  (:pr:`6161`)
- Relax precision for ppc64el  (:pr:`6222`)

Documentation
-------------
- Remove orphaned docs files  (:pr:`5832`)
- Array-like -> array_like  (:pr:`5929`)
- Change some more links to https  (:pr:`5937`)
- Fix self-contradictory minimum dependency versions  (:pr:`5939`)
- Fix formula for log-like in WLS  (:pr:`5946`)
- Fix typo  (:pr:`5949`)
- Add parameters for CountModel predict  (:pr:`5986`)
- Fix many spelling errors  (:pr:`5992`)
- Small fixups after the spell check  (:pr:`5994`)
- Clarify that GARCH models are deprecated  (:pr:`6000`)
- Added content for two headings in VAR docs  (:pr:`6022`)
- Fix regression doc strings  (:pr:`6031`)
- Add doc string check to doc build  (:pr:`6036`)
- Apply documentation standardizations  (:pr:`6038`)
- Fix spelling  (:pr:`6041`)
- Merge pull request #6041 from bashtage/doc-fixes  (:pr:`6042`)
- Fix notebook due to pandas index change  (:pr:`6044`)
- Remove warning due to deprecated features  (:pr:`6045`)
- Remove DynamicVAR  (:pr:`6046`)
- Small doc site improvements  (:pr:`6048`)
- Small fix ups for modernized size  (:pr:`6052`)
- More small doc fixes  (:pr:`6053`)
- Small changes to doc building  (:pr:`6054`)
- Use the working branch of numpy doc  (:pr:`6055`)
- Fix spelling in notebooks  (:pr:`6057`)
- Fix missing spaces around colon  (:pr:`6058`)
- Continue fixing docstring formatting  (:pr:`6060`)
- Fix web font size  (:pr:`6062`)
- Fix web font size  (:pr:`6063`)
- Fix doc errors affecting build  (:pr:`6067`)
- Improve docs in tools and ar_model  (:pr:`6080`)
- Improve filter docstrings  (:pr:`6082`)
- Spelling and notebook link  (:pr:`6085`)
- Website fix  (:pr:`6089`)
- Changes summary_col's docstring to match variables  (:pr:`6106`)
- Update spelling in CONTRIBUTING.rst  (:pr:`6107`)
- Update link in CONTRIBUTING.rst  (:pr:`6108`)
- Update PR template Numpy guide link  (:pr:`6110`)
- Added interpretations to LogitResults.get_margeff  (:pr:`6113`)
- Improve docstrings  (:pr:`6116`)
- Switch doc theme  (:pr:`6119`)
- Add initial API doc  (:pr:`6120`)
- Small improvements to docs  (:pr:`6122`)
- Switch doc icon  (:pr:`6123`)
- Fix doc build failure  (:pr:`6125`)
- Update templates and add missing API functions  (:pr:`6126`)
- Add missing functions from the API  (:pr:`6134`)
- Restructure the documentation  (:pr:`6136`)
- Add a new logo  (:pr:`6142`)
- Fix validator so that it works  (:pr:`6143`)
- Add formula API  (:pr:`6145`)
- Fix sidebar TOC  (:pr:`6160`)
- Warn that only trusted files should be unpickled  (:pr:`6162`)
- Update pickle warning  (:pr:`6166`)
- Fix warning format  (:pr:`6167`)
- Clarify req for cvxopt  (:pr:`6198`)
- Spelling and Doc String Fixes  (:pr:`6204`)
- Fix a typo  (:pr:`6214`)
- Fix typos in install.rst  (:pr:`6215`)
- Fix a typo  (:pr:`6216`)
- Docstring fixes  (:pr:`6235`)
- Fix spelling in notebooks  (:pr:`6257`)
- Clarify patsy 0.5.1 is required  (:pr:`6275`)
- Fix notebook  (:pr:`6279`)
- Close issues  (:pr:`6283`)
- Doc string changes  (:pr:`6289`)
- Correct spells  (:pr:`6298`)
- Add simple, documented script to get github info  (:pr:`6303`)
- Update test running instructions  (:pr:`6317`)
- Restore test() autosummary  (:pr:`6319`)
- Fix alpha description for GLMGam  (:pr:`6322`)
- Move api docs  (:pr:`6327`)
- Update Release Note  (:pr:`6342`)
- Fix documentation errors  (:pr:`6343`)
- Fixes in preparation for release  (:pr:`6344`)
- Further doc fixes  (:pr:`6345`)
- Fix minor doc errors  (:pr:`6347`)
- Git notes  (:pr:`6348`)
- Finalize release notes for 0.11  (:pr:`6349`)
- Add version dropdown  (:pr:`6350`)
- Finalize release note  (:pr:`6353`)
- Change generated path  (:pr:`6363`)
- Doc updates  (:pr:`6368`)
- Improve doc strings  (:pr:`6369`)
- Clarify demeaning in ljungbox  (:pr:`6390`)
- Fix ridge regression formula in hpfilter  (:pr:`6398`)
- Fix link  (:pr:`6407`)
- Update release note for v0.11.0rc2  (:pr:`6416`)
- Replace array with ndarray (:pr:`6447`)
- Final release note change for 0.11 (:pr:`6450`)


Maintenance
-----------
- Implement cached_value, cached_data proof of concept  (:pr:`4421`)
- Use Appender pattern for docstrings  (:pr:`5235`)
- Remove sandbox.formula, supplanted by patsy  (:pr:`5692`)
- Remove docstring'd-out traceback for code that no longer raises  (:pr:`5757`)
- Enable/mark mangled/commented-out tests  (:pr:`5768`)
- Implement parts of #5220, deprecate ancient aliases  (:pr:`5784`)
- Catch warnings produced during tests  (:pr:`5799`)
- Parts of iolib  (:pr:`5814`)
- E701 multiple statements on one line (colon)  (:pr:`5842`)
- Remove ex_pairwise file dominated by try_tukey_hsd  (:pr:`5856`)
- Fix pandas compat  (:pr:`5892`)
- Use pytest.raises to check error message  (:pr:`5897`)
- Bump dependencies  (:pr:`5910`)
- Fix pandas imports  (:pr:`5922`)
- Remove Python 2.7 from Appveyor  (:pr:`5927`)
- Relax tol on test that randomly fails  (:pr:`5931`)
- Fix test that fails with positive probability  (:pr:`5933`)
- Port parts of #5220  (:pr:`5935`)
- Remove Python 2.7 from travis  (:pr:`5938`)
- Fix linting failures  (:pr:`5940`)
- Drop redundant travis configs  (:pr:`5950`)
- Mark MPL test as MPL  (:pr:`5954`)
- Deprecate periodogram  (:pr:`5958`)
- Ensure seaborn is available for docbuild  (:pr:`5960`)
- Cython cleanups  (:pr:`5962`)
- Remove PY3  (:pr:`5965`)
- Remove future and Python 2.7  (:pr:`5969`)
- Remove string_types in favor of str  (:pr:`5972`)
- Restore ResettableCache  (:pr:`5976`)
- Cleanup legacy imports  (:pr:`5977`)
- Follow-up to #5956  (:pr:`5982`)
- Clarify breusch_pagan is for scalars  (:pr:`5984`)
- Add W605 to lint codes  (:pr:`5987`)
- Follow-up to #5928  (:pr:`5988`)
- Add spell checking  (:pr:`5990`)
- Remove comment no longer relevant  (:pr:`5991`)
- Refactor X13 testing  (:pr:`6001`)
- Standardized on nlags for acf/pacf  (:pr:`6002`)
- Rename forecast years to forecast periods  (:pr:`6007`)
- Improve testing of seasonal decompose  (:pr:`6011`)
- Remove notes about incorrect test  (:pr:`6015`)
- Turn relative import into an absolute import  (:pr:`6030`)
- Change types for future changes in NumPy  (:pr:`6039`)
- Move garch to archive/  (:pr:`6059`)
- Fix small lint issue  (:pr:`6066`)
- Stop testing on old, buggy SciPy  (:pr:`6069`)
- Small fixes in preparation for larger changes  (:pr:`6088`)
- Add tools for programatically manipulating docstrings  (:pr:`6090`)
- Ensure r download cache works  (:pr:`6092`)
- Fix new cache name  (:pr:`6093`)
- Fix wrong test  (:pr:`6094`)
- Remove extra LICENSE.txt and setup.cfg  (:pr:`6117`)
- Be compatible with scipy 1.3  (:pr:`6164`)
- Don't assume that 'python' is Python 3  (:pr:`6165`)
- Exclude pytest-xdist 1.30  (:pr:`6205`)
- Add Python 3.8 environment  (:pr:`6246`)
- Ignore vscode  (:pr:`6255`)
- Update test tolerance  (:pr:`6288`)
- Remove open_help method  (:pr:`6290`)
- Remove deprecated code in preparation for release  (:pr:`6291`)
- Deprecate recarray support  (:pr:`6310`)
- Reduce test size to prevent 32-bit crash  (:pr:`6311`)
- Remove chain dot  (:pr:`6312`)
- Catch and fix warnings  (:pr:`6313`)
- Use NumPy's linalg when available  (:pr:`6315`)
- Pin xdist  (:pr:`6392`)
- Unify pandas testing import  (:pr:`6394`)
- Clarify codecov  (:pr:`6406`)
- Port non-diagnostic changes  (:pr:`6412`)
- Fixes for future SciPY and pandas  (:pr:`6414`)
- Fixes for rc2  (:pr:`6419`)
- Switch to bionic  (:pr:`6423`)
- Improve test that randomly fails  (:pr:`6426`)
- Fix future issues  (:pr:`6440`)
- Disable cvxopt for windows  (:pr:`6445`)
- Reduce tolerance on basin hopping test  (:pr:`6448`)
- Remove unused import  (:pr:`6449`)

bug-wrong
---------

A new issue label `type-bug-wrong` indicates bugs that cause that incorrect
numbers are returned without warnings.
(Regular bugs are mostly usability bugs or bugs that raise an exception for
unsupported use cases.)
`see tagged issues <https://github.com/statsmodels/statsmodels/issues?q=is%3Aissue+label%3Atype-bug-wrong+is%3Aclosed+milestone%3A0.11/>`_


Major Bugs Fixed
================

See github issues for a list of bug fixes included in this release

- `Closed bugs <https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.11+label%3Atype-bug/>`_
- `Closed bugs (wrong result) <https://github.com/statsmodels/statsmodels/pulls?q=is%3Apr+is%3Amerged+milestone%3A0.11+label%3Atype-bug-wrong/>`_


Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

- Chad Fulton
- Brock Mendel
- Peter Quackenbush
- Kerby Shedden
- Kevin Sheppard

and the general maintainer and code reviewer

- Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.10 release (based on git log):

- Atticus Yang
- Austin Adams
- Balazs Varga
- Brock Mendel
- Chad Fulton
- Christian Clauss
- Emil Mirzayev
- Graham Inggs
- Guglielmo Saggiorato
- Hassan Kibirige
- Ian Preston
- Jefferson Tweed
- Josef Perktold
- Keller Scholl
- Kerby Shedden
- Kevin Sheppard
- Lucas Roberts
- Mandy Gu
- Omer Ozen
- Padarn Wilson
- Peter Quackenbush
- Qingqing Mao
- Rebecca N. Palmer
- Ron Itzikovitch
- Samesh Lakhotia
- Sandu Ursu
- Tim Staley
- Varun Sriram
- Yasine Gangat
- comatrion
- luxiform
- partev
- vegcev
- 郭飞


These lists of names are automatically generated based on git log, and may not
be complete.

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

- :pr:`4421`: ENH: Implement cached_value, cached_data proof of concept
- :pr:`5235`: STY: use Appender pattern for docstrings
- :pr:`5283`: ENH: Use cython fused types to simplify statespace code
- :pr:`5610`: BUG: Lilliefors min nobs not set
- :pr:`5692`: MAINT: remove sandbox.formula, supplanted by patsy
- :pr:`5735`: ENH: Allow fixing parameters in state space models
- :pr:`5757`: MAINT: Remove docstring'd-out traceback for code that no longer raises
- :pr:`5768`: WIP/TST: enable/mark mangled/commented-out tests
- :pr:`5784`: MAINT: implement parts of #5220, deprecate ancient aliases
- :pr:`5797`: TST: test for anova_nistcertified
- :pr:`5799`: TST: Catch warnings produced during tests
- :pr:`5814`: CLN: parts of iolib
- :pr:`5819`: CLN: robust
- :pr:`5821`: CLN: test_constrained
- :pr:`5826`: ENH/REF: Markov switching in log space: Hamilton filter / Kim smoother
- :pr:`5827`: ENH: Add new version of ARIMA-type estimators (AR, ARIMA, SARIMAX)
- :pr:`5832`: DOC: remove orphaned docs files
- :pr:`5841`: MAINT: remove no-longer-needed HC_se lookups
- :pr:`5842`: CLN: E701 multiple statements on one line (colon)
- :pr:`5852`: ENH: Dimension reduction for covariance matrices
- :pr:`5856`: MAINT: remove ex_pairwise file dominated by try_tukey_hsd
- :pr:`5892`: BUG: fix pandas compat
- :pr:`5893`: BUG: exponential smoothing - damped trend gives incorrect param, predictions
- :pr:`5895`: DOC: improvements to BayesMixedGLM docs, argument checking
- :pr:`5897`: MAINT: Use pytest.raises to check error message
- :pr:`5898`: ENH: use class to define MixedLM variance components structure
- :pr:`5903`: BUG: Fix kwargs update bug in linear model fit_regularized
- :pr:`5910`: MAINT: Bump dependencies
- :pr:`5915`: ENH: state space: add methods to apply fitted parameters to new observations or new dataset
- :pr:`5917`: BUG: TVTP for Markov regression
- :pr:`5921`: BUG: Ensure exponential smoothers has continuous double data
- :pr:`5922`: MAINT: Fix pandas imports
- :pr:`5926`: ENH: Add STL decomposition for time series
- :pr:`5927`: MAINT: Remove Python 2.7 from Appveyor
- :pr:`5928`: ENH: Add array_like function to simplify input checking
- :pr:`5929`: DOC: array-like -> array_like
- :pr:`5930`: BUG: Limit lags in KPSS
- :pr:`5931`: MAINT: Relax tol on test that randomly fails
- :pr:`5933`: MAINT: Fix test that fails with positive probability
- :pr:`5935`: CLN: port parts of #5220
- :pr:`5937`: DOC: Change some more links to https
- :pr:`5938`: MAINT: Remove Python 2.7 from travis
- :pr:`5939`: DOC: Fix self-contradictory minimum dependency versions
- :pr:`5940`: MAINT: Fix linting failures
- :pr:`5946`: DOC: Fix formula for log-like in WLS
- :pr:`5947`: PERF: Cythonize innovations algo and filter
- :pr:`5948`: ENH: Normalize eigenvectors from coint_johansen
- :pr:`5949`: DOC: Fix typo
- :pr:`5950`: MAINT: Drop redundant travis configs
- :pr:`5951`: BUG: Fix mosaic plot with missing category
- :pr:`5952`: ENH: Improve RESET test stability
- :pr:`5953`: ENH: Add type checkers/converts for int, float and bool
- :pr:`5954`: MAINT: Mark MPL test as MPL
- :pr:`5956`: BUG: Fix multidimensional model cov_params when using pandas
- :pr:`5957`: DOC: Clarify xname length and purpose
- :pr:`5958`: MAINT: Deprecate periodogram
- :pr:`5960`: MAINT: Ensure seaborn is available for docbuild
- :pr:`5962`: CLN: cython cleanups
- :pr:`5963`: ENH: Functional SIR
- :pr:`5964`: ENH: Add start_params to RLM
- :pr:`5965`: MAINT: Remove PY3
- :pr:`5966`: ENH: Add JohansenResults class
- :pr:`5967`: BUG/ENH: Improve RLM in the case of perfect fit
- :pr:`5968`: BUG: Make mad function behave correctly when used on empty inputs
- :pr:`5969`: MAINT: Remove future and Python 2.7
- :pr:`5971`: BUG: Fix a future issue in ExpSmooth
- :pr:`5972`: MAINT: Remove string_types in favor of str
- :pr:`5976`: MAINT: Restore ResettableCache
- :pr:`5977`: MAINT: Cleanup legacy imports
- :pr:`5982`: CLN: follow-up to #5956
- :pr:`5983`: BUG: Fix return for RegressionResults
- :pr:`5984`: MAINT: Clarify breusch_pagan is for scalars
- :pr:`5986`: DOC: Add parameters for CountModel predict
- :pr:`5987`: MAINT: add W605 to lint codes
- :pr:`5988`: CLN: follow-up to #5928
- :pr:`5990`: MAINT/DOC: Add spell checking
- :pr:`5991`: MAINT: Remove comment no longer relevant
- :pr:`5992`: DOC: Fix many spelling errors
- :pr:`5994`: DOC: Small fixups after the spell check
- :pr:`5995`: ENH: Add R-squared and Adj. R_squared to summary_col
- :pr:`5996`: BUG: Limit lags in KPSS
- :pr:`5997`: ENH/BUG: Add check to AR instance to prevent bugs
- :pr:`5998`: BUG: Replace alpha=0.05 with alpha=alpha
- :pr:`5999`: ENH: Add summary to AR
- :pr:`6000`: DOC: Clarify that GARCH models are deprecated
- :pr:`6001`: MAINT: Refactor X13 testing
- :pr:`6002`: MAINT: Standardized on nlags for acf/pacf
- :pr:`6003`: BUG: Do not fit when fit=False
- :pr:`6004`: ENH/BUG: Allow ARMA predict to swallow typ
- :pr:`6005`: MAINT: Ignore warns on 32 bit linux
- :pr:`6006`: BUG/ENH: Check exog in ARMA and ARIMA predict
- :pr:`6007`: MAINT: Rename forecast years to forecast periods
- :pr:`6008`: ENH: Allow GC testing for specific lags
- :pr:`6009`: TST: Verify categorical is supported for MNLogit
- :pr:`6010`: TST: Improve test that is failing due to precision issues
- :pr:`6011`: MAINT/BUG/TST: Improve testing of seasonal decompose
- :pr:`6012`: BUG: Fix t-test and f-test for multidimensional params
- :pr:`6014`: ENH: Zivot Andrews test
- :pr:`6015`: CLN: Remove notes about incorrect test
- :pr:`6016`: TST: Add check for dtypes in Binomial
- :pr:`6017`: ENH: Set limit for number of endog variables when using formulas
- :pr:`6018`: ENH: Improve ARMA startparams
- :pr:`6019`: BUG: Fix ARMA cov_params
- :pr:`6020`: TST: Correct test to use trend not level
- :pr:`6022`: DOC: added content for two headings in VAR docs
- :pr:`6023`: TST: Verify missing exog raises in ARIMA
- :pr:`6026`: WIP: Added Oaxaca-Blinder Decomposition
- :pr:`6028`: ENH: Add rolling WLS and OLS
- :pr:`6030`: MAINT: Turn relative import into an absolute import
- :pr:`6031`: DOC: Fix regression doc strings
- :pr:`6036`: BLD/DOC: Add doc string check to doc build
- :pr:`6038`: DOC: Apply documentation standardizations
- :pr:`6039`: MAINT: Change types for future changes in NumPy
- :pr:`6041`: DOC: Fix spelling
- :pr:`6042`: DOC: Merge pull request #6041 from bashtage/doc-fixes
- :pr:`6044`: DOC: Fix notebook due to pandas index change
- :pr:`6045`: DOC/MAINT: Remove warning due to deprecated features
- :pr:`6046`: DOC: Remove DynamicVAR
- :pr:`6048`: DOC: Small doc site improvements
- :pr:`6050`: BUG: MLEModel now passes nobs to Representation
- :pr:`6052`: DOC: Small fix ups for modernized size
- :pr:`6053`: DOC: More small doc fixes
- :pr:`6054`: DOC: Small changes to doc building
- :pr:`6055`: DOC: Use the working branch of numpy doc
- :pr:`6056`: MAINT: Prepare for Rolling Least Squares
- :pr:`6057`: DOC: Fix spelling in notebooks
- :pr:`6058`: DOC: Fix missing spaces around colon
- :pr:`6059`: REF: move garch to archive/
- :pr:`6060`: DOC: Continue fixing docstring formatting
- :pr:`6062`: DOC: Fix web font size
- :pr:`6063`: DOC: Fix web font size
- :pr:`6064`: ENH/PERF: Only perform required predict iterations in state space models
- :pr:`6066`: MAINT: Fix small lint issue
- :pr:`6067`: DOC: Fix doc errors affecting build
- :pr:`6069`: MAINT: Stop testing on old, buggy SciPy
- :pr:`6070`: BUG: Fix ARMA so that it works with exog when trend=nc
- :pr:`6071`: ENH: state space: Improve low memory usability; allow in fit, loglike
- :pr:`6072`: BUG: state space: cov_params computation in fix_params context
- :pr:`6073`: TST: Add conserve memory tests.
- :pr:`6074`: ENH: Improve cov_params in append, extend, apply
- :pr:`6075`: DOC: Clean tsatools docs
- :pr:`6077`: DOC: Improve regression doc strings
- :pr:`6079`: ENH/DOC: Improve Ljung-Box
- :pr:`6080`: DOC: Improve docs in tools and ar_model
- :pr:`6081`: BUG: Fix error introduced in isestimable
- :pr:`6082`: DOC: Improve filter docstrings
- :pr:`6085`: DOC: Spelling and notebook link
- :pr:`6087`: ENH: Replacement for AR
- :pr:`6088`: MAINT: Small fixes in preparation for larger changes
- :pr:`6089`: DOC: Website fix
- :pr:`6090`: ENH/DOC: Add tools for programatically manipulating docstrings
- :pr:`6091`: MAINT/SEC: Remove unnecessary pickle use
- :pr:`6092`: MAINT: Ensure r download cache works
- :pr:`6093`: MAINT: Fix new cache name
- :pr:`6094`: TST: Fix wrong test
- :pr:`6096`: DOC: Seasonality in SARIMAX Notebook
- :pr:`6102`: ENH: Improve SARIMAX start_params if too few nobs
- :pr:`6104`: BUG: Fix score computation with fixed params
- :pr:`6105`: BUG: Correct thresholding in correlation tools
- :pr:`6106`: DOC: Changes summary_col's docstring to match variables
- :pr:`6107`: DOC: Update spelling in CONTRIBUTING.rst
- :pr:`6108`: DOC: Update link in CONTRIBUTING.rst
- :pr:`6110`: DOC: Update PR template Numpy guide link
- :pr:`6111`: ENH: Add exact diffuse initialization as an option for SARIMAX, UnobservedComponents
- :pr:`6113`: DOC: added interpretations to LogitResults.get_margeff
- :pr:`6116`: DOC: Improve docstrings
- :pr:`6117`: MAINT: Remove extra LICENSE.txt and setup.cfg
- :pr:`6118`: DOC: Clarify summary2 documentation
- :pr:`6119`: DOC: Switch doc theme
- :pr:`6120`: DOC: Add initial API doc
- :pr:`6122`: DOC: Small improvements to docs
- :pr:`6123`: DOC: Switch doc icon
- :pr:`6124`: ENH: Plot only unique censored points
- :pr:`6125`: DOC: Fix doc build failure
- :pr:`6126`: DOC: Update templates and add missing API functions
- :pr:`6130`: BUG: Incorrect TSA index if loc resolves to slice
- :pr:`6131`: ENH: Compute standardized forecast error in diffuse period if possible
- :pr:`6133`: BUG: start_params for VMA model with exog.
- :pr:`6134`: DOC: Add missing functions from the API
- :pr:`6136`: DOC: Restructure the documentation
- :pr:`6142`: DOC: Add a new logo
- :pr:`6143`: DOC: Fix validator so that it works
- :pr:`6144`: BUG: use self.data consistently
- :pr:`6145`: DOC: Add formula API
- :pr:`6152`: BUG: Fix accepting of eval environment for formula
- :pr:`6160`: DOC: fix sidebar TOC
- :pr:`6161`: BLD: Travis CI: The sudo: tag is deprecated in Travis
- :pr:`6162`: DOC/SEC: Warn that only trusted files should be unpickled
- :pr:`6163`: ENH: Improve the cvxopt not found error
- :pr:`6164`: MAINT: Be compatible with scipy 1.3
- :pr:`6165`: MAINT: Don't assume that 'python' is Python 3
- :pr:`6166`: DOC: Update pickle warning
- :pr:`6167`: DOC: Fix warning format
- :pr:`6179`: ENH: Adds state space version of linear exponential smoothing models
- :pr:`6181`: ENH: state space: add wrapped states and, where possible, named states
- :pr:`6198`: DOC: Clarify req for cvxopt
- :pr:`6204`: DOC: Spelling and Doc String Fixes
- :pr:`6205`: MAINT: Exclude pytest-xdist 1.30
- :pr:`6208`: ENH: Scale parameter handling in GEE
- :pr:`6214`: DOC: fix a typo
- :pr:`6215`: DOC: fix typos in install.rst
- :pr:`6216`: DOC: fix a typo
- :pr:`6217`: BUG: Fix summary table header for mixedlm
- :pr:`6222`: MAINT: Relax precision for ppc64el
- :pr:`6227`: ENH: Add missing keyword argument to plot_acf
- :pr:`6231`: BUG: allow dynamic factor starting parameters computation with NaNs values
- :pr:`6232`: BUG: division by zero in exponential smoothing if damping_slope=0
- :pr:`6233`: BUG: dynamic factor model use AR model for error start params if error_var=False
- :pr:`6235`: DOC: docstring fixes
- :pr:`6239`: BUG: SARIMAX index behavior with simple_differencing=True
- :pr:`6240`: BUG: parameter names in DynamicFactor for unstructured error covariance matrix
- :pr:`6241`: BUG: SARIMAX: basic validation for order, seasonal_order
- :pr:`6242`: BUG: Forecasts now ignore non-monotonic period index
- :pr:`6246`: TST: Add Python 3.8 environment
- :pr:`6250`: ENH: Update SARIMAX to use SARIMAXSpecification for more consistent input handling
- :pr:`6254`: ENH: State space: Add finer-grained memory conserve settings
- :pr:`6255`: MAINT: Ignore vscode
- :pr:`6257`: DOC: Fix spelling in notebooks
- :pr:`6258`: BUG: Hannan-Rissanen third stage is invalid if non-stationary/invertible
- :pr:`6260`: BUG: cloning of arima.ARIMA models.
- :pr:`6261`: BUG: state space: saving fixed params w/ extend, apply, append
- :pr:`6266`: ENH: and vlines option to plot_fit
- :pr:`6275`: MAINT/DOC: Clarify patsy 0.5.1 is required
- :pr:`6279`: DOC: Fix notebook
- :pr:`6280`: ENH: State space: Improve simulate, IRF, prediction
- :pr:`6281`: BUG: Pass arguments through in plot_leverage_resid2
- :pr:`6283`: MAINT/DOC: Close issues
- :pr:`6285`: BUG: Raise in GC test for VAR(0)
- :pr:`6286`: BUG: Correct VAR summary when model contains exog variables
- :pr:`6288`: MAINT: Update test tolerance
- :pr:`6289`: DOC: doc string changes
- :pr:`6290`: MAINT: Remove open_help method
- :pr:`6291`: MAINT: Remove deprecated code in preparation for release
- :pr:`6292`: BUG: Ensure BW is not 0
- :pr:`6294`: ENH: better argument checking for StratifiedTable
- :pr:`6297`: BUG: Fix conf interval with MI
- :pr:`6298`: DOC: Correct spells
- :pr:`6299`: DOC: Add example notebook for GEE score tests
- :pr:`6303`: DOC/MAINT: Add simple, documented script to get github info
- :pr:`6310`: MAINT: Deprecate recarray support
- :pr:`6311`: TST: Reduce test size to prevent 32-bit crash
- :pr:`6312`: MAINT: Remove chain dot
- :pr:`6313`: MAINT: Catch and fix warnings
- :pr:`6314`: BUG: Check dtype in KDEUnivariate
- :pr:`6315`: MAINT: Use NumPy's linalg when available
- :pr:`6316`: MAINT: Workaround NumPy ptp issue
- :pr:`6317`: DOC: Update test running instructions
- :pr:`6318`: BUG: Ensure inputs are finite in granger causality test
- :pr:`6319`: DOC: Restore test() autosummary
- :pr:`6320`: BUG: Restore multicomp
- :pr:`6321`: BUG: Fix trend due to recent changes
- :pr:`6322`: DOC: fix alpha description for GLMGam
- :pr:`6324`: ENH: Improve Ljung Box diagnostics
- :pr:`6327`: DOC: Move api docs
- :pr:`6332`: DEPR: state space: deprecate out-of-sample w/ unsupported index
- :pr:`6333`: BUG: state space: integer params can cause imaginary output
- :pr:`6334`: ENH: append, extend check that index matches model
- :pr:`6337`: BUG: fix k_exog, k_trend in arima.ARIMA; raise error when cloning a model with exog if no new exog given
- :pr:`6338`: DOC: Documentation for release v0.11
- :pr:`6340`: BUG: fix _get_index_loc with date strings
- :pr:`6342`: DOC: Update Release Note
- :pr:`6343`: DOC: Fix documentation errors
- :pr:`6344`: DOC: Fixes in preparation for release
- :pr:`6345`: DOC: Further doc fixes
- :pr:`6347`: DOC: Fix minor doc errors
- :pr:`6348`: DOC: git notes
- :pr:`6349`: DOC: Finalize release notes for 0.11
- :pr:`6350`: DOC: Add version dropdown
- :pr:`6353`: DOC: Finalize release note
- :pr:`6363`: DOC: Change generated path
- :pr:`6368`: Doc updates
- :pr:`6369`: DOC: Improve doc strings
- :pr:`6374`: BUG: Correct standardization in robust skewness
- :pr:`6375`: ENH: Supporting custom kernel in local linear kernel regression
- :pr:`6389`: BUG: Use correct exog names
- :pr:`6390`: DOC: Clarify demeaning in ljungbox
- :pr:`6391`: BUG: Fix axis in irq
- :pr:`6392`: MAINT: Pin xdist
- :pr:`6394`: MAINT: Unify pandas testing import
- :pr:`6398`: DOC: fix ridge regression formula in hpfilter
- :pr:`6401`: ENH: Distance dependence measures
- :pr:`6406`: MAINT: Clarify codecov
- :pr:`6407`: DOC: Fix link
- :pr:`6410`: ENH/CLN: Improve diagnostics
- :pr:`6412`: CLN/MAINT: Port non-diagnostic changes
- :pr:`6414`: CLN: Fixes for future SciPY and pandas
- :pr:`6415`: BUG: RecursiveLS should not allow `fix_params` method.
- :pr:`6416`: DOC: Update release note for v0.11.0rc2
- :pr:`6419`: MAINT: Fixes for rc2
- :pr:`6422`: BUG: Improve executable detection
- :pr:`6423`: MAINT: Switch to bionic
- :pr:`6424`: REF: More descriptive error message for recursive least squares parameter constraints.
- :pr:`6425`: BUG/ENH: Diffuse multivariate case w/ non-diagonal observation innovation covariance matrix
- :pr:`6426`: TST: Improve test that randomly fails
- :pr:`6440`: MAINT: Fix future issues
- :pr:`6442`: BUG: Fix bug in ridge for vector alpha
- :pr:`6445`: MAINT: Disable cvxopt for windows
- :pr:`6447`: DOC: Replace array with ndarray
- :pr:`6448`: MAINT: Reduce tolerance on basin hopping test
- :pr:`6449`: MAINT: Remove unused import
- :pr:`6450`: DOC: Final release note change for 0.11
:orphan:

==============
Release 0.13.0
==============

Release summary
===============

statsmodels is using github to store the updated documentation. Two version are available:

- `Stable <https://www.statsmodels.org/>`_, the latest release
- `Development <https://www.statsmodels.org/devel/>`_, the latest build of the main branch

**Warning**

API stability is not guaranteed for new features, although even in
this case changes will be made in a backwards compatible way if
possible. The stability of a new feature depends on how much time it
was already in statsmodels main and how much usage it has already
seen.  If there are specific known problems or limitations, then they
are mentioned in the docstrings.

Stats
-----
**Issues Closed**: 238

**Pull Requests Merged**: 165


The Highlights
==============

New cross-sectional models
--------------------------

Beta Regression
~~~~~~~~~~~~~~~
:class:`~statsmodels.othermod.betareg.BetaModel` estimates a regression model
for dependent variable in the unit interval such as fractions and proportions
based on the Beta distribution. The Model is parameterized by mean and
precision, where both can depend on explanatory variables through link
functions.

Ordinal Regression
~~~~~~~~~~~~~~~~~~
:class:`statsmodels.miscmodels.ordinal_model.OrderedModel` implements
cumulative link models for ordinal data, based on Logit, Probit or a
userprovided CDF link.

Distributions
-------------

Copulas
~~~~~~~
Statsmodels includes now basic support for mainly bivariate copulas.
Currently, 10 copulas are available, Archimedean, elliptical and asymmetric
extreme value copulas.
:class:`~statsmodels.distributions.copula.api.CopulaDistribution` combines a
copula with marginal distributions to create multivariate distributions.

Count distribution based on discretization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`~statsmodels.distributions.discrete.DiscretizedCount` provides count
distributions generated by discretizing continuous distributions available
in scipy. The parameters of the distribution can be estimated by maximum
likelihood with :class:`~statsmodels.distributions.discrete.DiscretizedModel`.

Bernstein Distribution
~~~~~~~~~~~~~~~~~~~~~~
:class:`~statsmodels.distributions.berstein.BernsteinDistribution` creates
nonparametric univariate and multivariate distributions using Bernstein
polynomials on a regular grid. This can be used to smooth histograms or
approximate distributions on the unit hypercube. When the marginal distributions
are uniform, then the BernsteinDistribution is a copula.

Statistics
----------

Brunner Munzel rank comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Brunner-Munzel test is nonparametric comparison of two samples and is an
extension of Wilcoxon-Mann-Whitney and Fligner-Policello tests that requires
only ordinal information without further assumption on
the distributions of the samples. Statsmodels provides the Brunner Munzel
hypothesis test for stochastic equality in
:func:`~statsmodels.stats.nonparametric.rank_compare_2indep` but also
confidence intervals and equivalence testing (TOST) for the stochastically
larger statistic, also known as Common Language effect size.

Nonparametric
-------------

Asymmetric kernels
~~~~~~~~~~~~~~~~~~

Asymmetric kernels can nonparametrically estimate density and cumulative
distribution function for random variables that have limited support, either
unit interval or positive or nonnegative real line. Beta kernels are available
for data in the unit interval. The available kernels for positive data are
“gamma”, “gamma2”, “bs”, “invgamma”, “invgauss”, “lognorm”, “recipinvgauss”
and “weibull”
:func:`~statsmodels.nonparametric.kernels_asymmetric.pdf_kernel_asym` estimates
a kernel density given a bandwidth parameter.
:func:`~statsmodels.nonparametric.kernels_asymmetric.cdf_kernel_asym`
estimates a kernel cdf.


Time series analysis
--------------------

Autoregressive Distributed Lag Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`~statsmodels.tsa.ardl.ARDL` adds support for specifying and estimating ARDL
models, and :class:`~statsmodels.tsa.ardl.UECM` support specifying models in error
correction form. :func:`~statsmodels.tsa.ardl.ardl_select_order` simplifies selecting
both AR and DL model orders. :func:`~statsmodels.tsa.ardl.UECM.bounds_test` implements
the bounds test of Peseran, Shin and Smith (2001) for testing whether there is a levels
relationship without knowing teh orders of integration of the variables.

.. ipython:: python

   from statsmodels.datasets import danish_data
   import statsmodels.tsa.api as tsa

   data = danish_data.load().data
   sel = tsa.ardl_select_order(data.lrm, 3, data[["lry", "ibo", "ide"]], 3, ic="aic")
   ardl = sel.model
   ardl.ardl_order

.. ipython:: python

   res = ardl.fit()
   print(res.summary())

.. ipython:: python

   uecm = tsa.UECM.from_ardl(ardl)
   uecm_res = uecm.fit()
   uecm_res.bounds_test(case=4)


Fixed parameters in ARIMA estimators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Allow fixing parameters in ARIMA estimator Hannan-Rissanen
  (:func:`~statsmodels.tsa.arima.estimators.hannan_rissanen`) through the new
  ``fixed_params`` argument


What's new - an overview
========================

The following lists the main new features of statsmodels 0.13.0. In addition,
release 0.13.0 includes bug fixes, refactorings and improvements in many areas.

Major Feature
-------------
- Allow fixing parameters in ARIMA estimator Hannan-Rissanen
  (:pr:`7497`, :pr:`7501`)
- OLS add "slim" option to summary method (:pr:`7693` based on :pr:`6880`)
- Add loglog link for use with GLM (:pr:`7594`)
- improved default derivatives in CDFLink (:pr:`7287`)
- GLM enhanced and corrected ``get_distribution`` (:pr:`7535`)
- GLMResults info_criteria, add ``dk_params`` option to include scale in
  parameter count (:pr:`7693`)
- GLMResults add pseudo R-squared, Cox-Snell and McFadden
  (:pr:`7682` based on :pr:`7367`)
- nonparametric: add tricube kernel (:pr:`7697` based on :pr:`7671`)


Submodules
----------


``Documentation``
~~~~~~~~~~~~~~~~~
- Port missed doc fix  (:pr:`7123`)
- Rls note  (:pr:`7293`)
- Minor updates to v0.12.2 release notes  (:pr:`7303`)
- Update doc for tweedie allowed links  (:pr:`7395`)
- Don't point to release version  (:pr:`7399`)
- Fixed error in linear mixed effects example  (:pr:`7402`)
- Fix for upstream changes in PyMC3 notebook  (:pr:`7416`)
- Remove redundant words in PCA docstring  (:pr:`7423`)
- Misc fixes in docstr of fdrcorrection  (:pr:`7426`)
- Small doc fixes  (:pr:`7434`)
- Typo, plats->plots  (:pr:`7458`)
- Specify impulse to impulse_responses in VARMAX notebook  (:pr:`7475`)
- Fix errors in theta notebook  (:pr:`7539`)
- Add github actions to build docs  (:pr:`7540`)
- Fix GH actions  (:pr:`7541`)
- Continue working on it  (:pr:`7552`)
- Continue working on push ability  (:pr:`7553`)
- Continue working on push ability  (:pr:`7554`)
- Finalize push ability  (:pr:`7555`)
- Finalize push ability  (:pr:`7556`)
- Finalize push ability  (:pr:`7557`)
- Finalize push ability  (:pr:`7558`)
- Get doc push to work  (:pr:`7559`)
- Get doc push to work  (:pr:`7560`)
- Get doc push to work  (:pr:`7561`)
- Improve rolling OLS notebook  (:pr:`7572`)
- Correct docstring  (:pr:`7587`)
- Copula in user guide and examples  (:pr:`7607`)
- Improve ARDL and documentation  (:pr:`7611`)
- Clarify which series is on x-axis  (:pr:`7612`)
- Small clean of example  (:pr:`7614`)
- Spelling error in docs fixed  (:pr:`7618`)
- Update dev page flake8 command to follow PULL_REQUEST_TEMPLATE.md  (:pr:`7644`)
- Improve copula notebook  (:pr:`7651`)
- Remove duplication methods section  (:pr:`7676`)
- Second try ixing duplicate methods  (:pr:`7677`)
- Fix a typo  (:pr:`7681`)
- Improve ARDL notebook  (:pr:`7699`)
- Update versions.json  (:pr:`7702`)
- Update versions file  (:pr:`7708`)
- Update release note  (:pr:`7714`)
- Update release note  (:pr:`7726`)
- Correct MultivariateTestResults doc string  (:pr:`7735`)
- Correct MultivariateTestResults doc string  (:pr:`7738`)
- Add missing function doc head  (:pr:`7740`)
- More 0.13  (:pr:`7757`)
- Fix lowess notebook  (:pr:`7770`)



``Performance``
~~~~~~~~~~~~~~~
- Added fft to ccovf and ccf  (:pr:`7721`)
- Improve Lowess  (:pr:`7768`)



``backport``
~~~~~~~~~~~~
- Backports  (:pr:`7222`)
- Backports  (:pr:`7291`)
- Forecast after extend w/ time varying matrix  (:pr:`7437`)



``base``
~~~~~~~~
- Use np.linalg.solve() instead of np.linalg.inv() in Newton-Raphson Algorithm  (:pr:`7429`)
- Allow remove_data to work when an attribute is not implemented  (:pr:`7511`)
- REF/BUG generic likelihood LLRMixin use df_resid instead of df_model for llr_pvalue  (:pr:`7586`)
- Raise when invalid optimization options passed to optimizer  (:pr:`7596`)



``datasets``
~~~~~~~~~~~~
- Add an error message for not found data  (:pr:`7490`)



``discrete``
~~~~~~~~~~~~
- Add discretized count distribution  (:pr:`7488`)
- ZI predict, fix offset default if None, allow exog_infl None if constant  (:pr:`7670`)



``distributions``
~~~~~~~~~~~~~~~~~
- Copula 7254 rebased  (:pr:`7408`)
- Add discretized count distribution  (:pr:`7488`)
- Random number generation wrapper for rng, qrng  (:pr:`7608`)
- BUG/REF copula another round for 0.13  (:pr:`7648`)
- Temporarily change the default RNG in check_random_state  (:pr:`7652`)
- More copula improvements for 0.13  (:pr:`7723`)



``docs``
~~~~~~~~
- Fix for upstream changes in PyMC3 notebook  (:pr:`7416`)
- Correct small typo in Theta model Notebook  (:pr:`7450`)
- Prevent indent running on None  (:pr:`7462`)
- Update versions file  (:pr:`7708`)
- Improve docs and docstrings, mainly for recent additions  (:pr:`7727`)
- Api.py, docstring improvements  (:pr:`7732`)
- Add to release notes, smaller doc fixes, references  (:pr:`7743`)



``genmod``
~~~~~~~~~~
- Change default derivative in CDFLink  (:pr:`7287`)
- Allow user to configure GEE qic  (:pr:`7471`)
- Score and Hessian for Tweedie models  (:pr:`7489`)
- BUG/ENH fix and enh GLM, family get_distribution  (:pr:`7535`)
- Enh glm loglog  (:pr:`7594`)
- McFadden and Cox&Snell Pseudo R squared to GLMResults  (:pr:`7682`)
- Add dk_params option to GLM info_criteria  (:pr:`7693`)
- Warn kwargs glm  (:pr:`7750`)
- GLM init invalid kwargs use ValueWarning  (:pr:`7751`)



``graphics``
~~~~~~~~~~~~
- Fix UserWarning: marker is redundantly defined (Matplotlib v 3.4.1)   (:pr:`7400`)
- Fix axis labels in qqplots  (:pr:`7413`)
- Remove typo in plot_pacf example  (:pr:`7514`)
- Start process of changing default in plot-pacf  (:pr:`7582`)
- Improve limit format in diff plot  (:pr:`7592`)
- Clarify which series is on x-axis  (:pr:`7612`)
- Graphics.plot_partregress add eval_env options  (:pr:`7673`)



``io``
~~~~~~
- Add support for pickling for generic path-like objects  (:pr:`7581`)
- Fix summary().as_latex, line in top table dropped  (:pr:`7748`)



``maintenance``
~~~~~~~~~~~~~~~
- V0.12.1 backports  (:pr:`7121`)
- Backport fixes for 0.12.2 compat release  (:pr:`7221`)
- Fix descriptive stats with extension dtypes  (:pr:`7404`)
- Fix pip pre test failures  (:pr:`7405`)
- Fix README badges  (:pr:`7406`)
- Silence warnings and future compat  (:pr:`7425`)
- Use loadscope to avoid rerunning setup  (:pr:`7432`)
- Remove cyclic import risks  (:pr:`7438`)
- Fit future and deprecation warnings  (:pr:`7474`)
- Avoid future issues in pandas  (:pr:`7495`)
- Remove 32-bit testing  (:pr:`7536`)
- Fix contrasts for Pandas changes  (:pr:`7546`)
- Correct example implementation  (:pr:`7547`)
- Check push ability  (:pr:`7551`)
- Remove deprecated functions  (:pr:`7575`)
- Remove additional deprecated features  (:pr:`7577`)
- Remove recarray  (:pr:`7578`)
- Remove deprecated code  (:pr:`7579`)
- Correct notebooks for deprecations  (:pr:`7580`)
- Fix spelling errors  (:pr:`7583`)
- Clarify minimum versions  (:pr:`7590`)
- Revert exception to warning  (:pr:`7599`)
- Silence future warnings  (:pr:`7617`)
- Avoid passing bad optimization param  (:pr:`7620`)
- Pin matplotlib  (:pr:`7641`)
- Modernize prediction in notebooks  (:pr:`7649`)
- Protect against changes in numeric indexes  (:pr:`7685`)
- Final issues in `__all__`  (:pr:`7742`)
- Fix hard to reach errors  (:pr:`7744`)



``multivariate``
~~~~~~~~~~~~~~~~
- Multivariate - Return E and H matrices in dict  (:pr:`5491`)
- Added the option `full_matrices=False` in the PCA method  (:pr:`7329`)
- Factor fit ml em resets seed (rebased)  (:pr:`7703`)
- Correct MultivariateTestResults doc string  (:pr:`7735`)
- Correct MultivariateTestResults doc string  (:pr:`7738`)
- Add missing function doc head  (:pr:`7740`)



``nonparametric``
~~~~~~~~~~~~~~~~~
- ENH add  tricube kernel  (:pr:`7697`)
- Fix lowess spikes/nans from epsilon values  (:pr:`7766`)
- Improve Lowess  (:pr:`7768`)



``othermod``
~~~~~~~~~~~~
- Betareg rebased3 Beta regression  (:pr:`7543`)
- REF/BUG generic likelihood LLRMixin use df_resid instead of df_model for llr_pvalue  (:pr:`7586`)
- Oaxaca Variance/Other Models  (:pr:`7713`)



``regression``
~~~~~~~~~~~~~~
- Allow remove_data to work when an attribute is not implemented  (:pr:`7511`)
- Fix scale parameter in elastic net  (:pr:`7571`)
- Regression, allow remove_data to remove wendog, wexog, wresid  (:pr:`7595`)
- Spelling error in docs fixed  (:pr:`7618`)
- Add dk_params option to GLM info_criteria  (:pr:`7693`)
- Quantile regression use dimension of x matrix rather than rank  (:pr:`7694`)
- Add option for slim summary in OLS results  (:pr:`7696`)
- Enable VIF to work with DataFrames  (:pr:`7704`)



``stats``
~~~~~~~~~
- Runs test numeric cutoff error  (:pr:`7422`)
- Resolve TODO in proportion.py  (:pr:`7515`)
- Improve sidak multipletest precision close to zero  (:pr:`7668`)
- Proportions_chisquare prevent integer overflow   (:pr:`7669`)
- Fix lilliefors results for single-column DataFrames  (:pr:`7698`)
- Describe / Description do not return percentiles  (:pr:`7710`)
- ENH: add options to meta-analysis plot_forest (:pr:`7772`)


``tools``
~~~~~~~~~
- Change default derivative in CDFLink  (:pr:`7287`)
- Fix style issue  (:pr:`7739`)



``tsa``
~~~~~~~
- Add Helper function to solve for polynomial coefficients from roots for ARIMA  (:pr:`6921`)
- Changed month abbreviations with localization  (:pr:`7409`)
- Add ARDL model  (:pr:`7433`)
- Fix typo in ets error  (:pr:`7435`)
- Add fixed_params to Hannan Rissanen (GH7202)  (:pr:`7497`)
- Enable ARIMA.fit(method='hannan_rissanen') with fixed parameters (GH7501)  (:pr:`7502`)
- Fix errors when making dynamic forecasts  (:pr:`7516`)
- Correct index location of seasonal  (:pr:`7545`)
- Handle non-date index with a freq    (:pr:`7574`)
- Start process of changing default in plot-pacf  (:pr:`7582`)
- Correct docstring  (:pr:`7587`)
- Let VAR results complete when model has perfect fit  (:pr:`7588`)
- Rename nc to n everywhere  (:pr:`7593`)
- Improve ARDL and documentation  (:pr:`7611`)
- Add RUR stationarity test to statsmodels.tsa.stattools  (:pr:`7616`)
- Improve ARDL and UECM  (:pr:`7619`)
- Improve error message in seasonal for bad freq  (:pr:`7643`)
- ENH Fixed Range Unit-Root critical values  (:pr:`7645`)
- Add SARIMAX FAQ  (:pr:`7656`)
- Add to the SARIMAX FAQ  (:pr:`7659`)
- Improve SARIMAX FAQ Notebook  (:pr:`7661`)
- Improve ARIMA documentation  (:pr:`7662`)
- Update TSA Api  (:pr:`7701`)
- Correct ArmaProcess.from_estimation  (:pr:`7709`)
- Added fft to ccovf and ccf  (:pr:`7721`)



``tsa.statespace``
~~~~~~~~~~~~~~~~~~
- Port missed doc fix  (:pr:`7123`)
- Forecast after extend w/ time varying matrix  (:pr:`7437`)
- Specify impulse to impulse_responses in VARMAX notebook  (:pr:`7475`)
- Column name can be passed as an argument in `impulse_responses` in `VARMAX`  (:pr:`7506`)
- Statespace MLEModel false validation error with nested fix_params (GH7507)  (:pr:`7508`)
- Ensure attributes exist  (:pr:`7538`)
- Ensure warning does not raise  (:pr:`7589`)
- Assert correct iloc dtypes  (:pr:`7737`)



``tsa.vector.ar``
~~~~~~~~~~~~~~~~~
- Fix float index usage in IRF error bands  (:pr:`7397`)
- Add error if too few values  (:pr:`7591`)





bug-wrong
---------

A new issue label `type-bug-wrong` indicates bugs that cause that incorrect
numbers are returned without warnings.
(Regular bugs are mostly usability bugs or bugs that raise an exception for
unsupported use cases.)
`see tagged issues <https://github.com/statsmodels/statsmodels/issues?q=is%3Aissue+label%3Atype-bug-wrong+is%3Aclosed+milestone%3A0.13/>`_


Major Bugs Fixed
================

See github issues for a list of bug fixes included in this release

- `Closed bugs <https://github.com/statsmodels/statsmodels/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Amerged+milestone%3A0.13+label%3Atype-bug/>`_
- `Closed bugs (wrong result) <https://github.com/statsmodels/statsmodels/pulls?q=is%3Apr+is%3Amerged+milestone%3A0.13+label%3Atype-bug-wrong/>`_


Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

- Chad Fulton
- Brock Mendel
- Peter Quackenbush
- Kerby Shedden
- Kevin Sheppard

and the general maintainer and code reviewer

- Josef Perktold

Additionally, many users contributed by participation in github issues and
providing feedback.

Thanks to all of the contributors for the 0.13.0 release (based on git log):

- Aidan Russell
- Alexander Stiebing
- Austin Adams
- Ben Greiner
- Brent Pedersen
- Chad Fulton
- Chadwick Boulay
- Edwin Rijgersberg
- Ezequiel Smucler
- G. D. Mcbain
- Graham Inggs
- Greg Mcmahan
- Helder Oliveira
- Hsiao Yi
- Jack Liu
- Jake Jiacheng Liu
- Jeremy Bejarano
- Joris Van Den Bossche
- Josef Perktold
- Juan Orduz
- Kerby Shedden
- Kevin Sheppard
- Luke Gregor
- Malte Zietlow
- Masanori Kanazu
- Max Mahlke
- Michele Fortunato
- Mike Ovyan
- Min Rk
- Natalie Heer
- Nikolai Korolev
- Omar Gutiérrez
- Oswaldo
- Pamphile Roy
- Pratyush Sharan
- Roberto Nunes Mourão
- Simardeep27
- Simon Høxbro Hansen
- Sin Kim
- Skipper Seabold
- Stefan Appelhoff
- Thomas Brooks
- Tomohiro Endo
- Wahram Andrikyan
- cxan96
- janosbiro
- partev
- w31ha0


These lists of names are automatically generated based on git log, and may not
be complete.

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

- :pr:`5491`: ENH: multivariate - Return E and H matrices in dict
- :pr:`6921`: ENH: Add Helper function to solve for polynomial coefficients from roots for ARIMA
- :pr:`7121`: MAINT: v0.12.1 backports
- :pr:`7123`: DOC: Port missed doc fix
- :pr:`7221`: MAINT: Backport fixes for 0.12.2 compat release
- :pr:`7222`: Backports
- :pr:`7287`: REF: change default derivative in CDFLink
- :pr:`7291`: Backports
- :pr:`7293`: Rls note
- :pr:`7303`: DOC: Minor updates to v0.12.2 release notes
- :pr:`7329`: ENH: Added the option `full_matrices=False` in the PCA method
- :pr:`7395`: DOC: update doc for tweedie allowed links
- :pr:`7397`: BUG: Fix float index usage in IRF error bands
- :pr:`7399`: DOC: Don't point to release version
- :pr:`7400`: MAINT: Fix UserWarning: marker is redundantly defined (Matplotlib v 3.4.1)
- :pr:`7402`: DOC: fixed error in linear mixed effects example
- :pr:`7404`: MAINT: Fix descriptive stats with extension dtypes
- :pr:`7405`: MAINT: Fix pip pre test failures
- :pr:`7406`: MAINT: Fix README badges
- :pr:`7408`: Copula 7254 rebased
- :pr:`7409`: ENH: changed month abbreviations with localization
- :pr:`7413`: BUG: Fix axis labels in qqplots
- :pr:`7416`: MAINT: Fix for upstream changes in PyMC3 notebook
- :pr:`7422`: BUG: Runs test numeric cutoff error
- :pr:`7423`: DOC/MAINT: Remove redundant words in PCA docstring
- :pr:`7425`: MAINT: Silence warnings and future compat
- :pr:`7426`: DOC: misc fixes in docstr of fdrcorrection
- :pr:`7429`: ENH: Use np.linalg.solve() instead of np.linalg.inv() in Newton-Raphson Algorithm
- :pr:`7432`: MAINT: Use loadscope to avoid rerunning setup
- :pr:`7433`: ENH: Add ARDL model
- :pr:`7434`: DOC: Small doc fixes
- :pr:`7435`: fix typo in ets error
- :pr:`7437`: BUG: forecast after extend w/ time varying matrix
- :pr:`7438`: MAINT: Remove cyclic import risks
- :pr:`7450`: Correct small typo in Theta model Notebook
- :pr:`7458`: DOC: typo, plats->plots
- :pr:`7462`: BUG: Prevent indent running on None
- :pr:`7471`: ENH: Allow user to configure GEE qic
- :pr:`7474`: MAINT: Fit future and deprecation warnings
- :pr:`7475`: Specify impulse to impulse_responses in VARMAX notebook
- :pr:`7488`: ENH: add discretized count distribution
- :pr:`7489`: BUG: score and Hessian for Tweedie models
- :pr:`7490`: ENH: Add an error message for not found data
- :pr:`7495`: MAINT: Avoid future issues in pandas
- :pr:`7497`: ENH: Add fixed_params to Hannan Rissanen (GH7202)
- :pr:`7502`: ENH: Enable ARIMA.fit(method='hannan_rissanen') with fixed parameters (GH7501)
- :pr:`7506`: ENH: Column name can be passed as an argument in `impulse_responses` in `VARMAX`
- :pr:`7508`: BUG: statespace MLEModel false validation error with nested fix_params (GH7507)
- :pr:`7511`: Allow remove_data to work when an attribute is not implemented
- :pr:`7514`: Remove typo in plot_pacf example
- :pr:`7515`: resolve TODO in proportion.py
- :pr:`7516`: BUG: Fix errors when making dynamic forecasts
- :pr:`7535`: BUG/ENH fix and enh GLM, family get_distribution
- :pr:`7536`: MAINT: Remove 32-bit testing
- :pr:`7538`: BUG: Ensure attributes exist
- :pr:`7539`: DOC: Fix errors in theta notebook
- :pr:`7540`: MAINT: Add github actions to build docs
- :pr:`7541`: MAINT: Fix GH actions
- :pr:`7543`: Betareg rebased3 Beta regression
- :pr:`7545`: BUG: Correct index location of seasonal
- :pr:`7546`: MAINT: Fix contrasts for Pandas changes
- :pr:`7547`: MAINT: Correct example implementation
- :pr:`7551`: MAINT: Check push ability
- :pr:`7552`: MAINT: Continue working on it
- :pr:`7553`: MAINT: Continue working on push ability
- :pr:`7554`: MAINT: Continue working on push ability
- :pr:`7555`: MAINT: Finalize push ability
- :pr:`7556`: MAINT: Finalize push ability
- :pr:`7557`: MAINT: Finalize push ability
- :pr:`7558`: MAINT: Finalize push ability
- :pr:`7559`: MAINT: Get doc push to work
- :pr:`7560`: MAINT: Get doc push to work
- :pr:`7561`: MAINT: Get doc push to work
- :pr:`7571`: BUG: Fix scale parameter in elastic net
- :pr:`7572`: DOC: Improve rolling OLS notebook
- :pr:`7574`: BUG: Handle non-date index with a freq
- :pr:`7575`: MAINT: Remove deprecated functions
- :pr:`7577`: MAINT: Remove additional deprecated features
- :pr:`7578`: MAINT: Remove recarray
- :pr:`7579`: MAINT: Remove deprecated code
- :pr:`7580`: MAINT: Correct notebooks for deprecations
- :pr:`7581`: ENH: Add support for pickling for generic path-like objects
- :pr:`7582`: ENH: Start process of changing default in plot-pacf
- :pr:`7583`: MAINT: Fix spelling errors
- :pr:`7586`: REF/BUG generic likelihood LLRMixin use df_resid instead of df_model for llr_pvalue
- :pr:`7587`: DOC: Correct docstring
- :pr:`7588`: BUG: Let VAR results complete when model has perfect fit
- :pr:`7589`: BUG: Ensure warning does not raise
- :pr:`7590`: MAINT: Clarify minimum versions
- :pr:`7591`: ENH: Add error if too few values
- :pr:`7592`: ENH: Improve limit format in diff plot
- :pr:`7593`: MAINT: Rename nc to n everywhere
- :pr:`7594`: Enh glm loglog
- :pr:`7595`: BUG: regression, allow remove_data to remove wendog, wexog, wresid
- :pr:`7596`: ENH: Raise when invalid optimization options passed to optimizer
- :pr:`7599`: MAINT: Revert exception to warning
- :pr:`7607`: DOC: copula in user guide and examples
- :pr:`7608`: ENH: random number generation wrapper for rng, qrng
- :pr:`7611`: ENH: Improve ARDL and documentation
- :pr:`7612`: BUG/DOC: Clarify which series is on x-axis
- :pr:`7614`: DOC: Small clean of example
- :pr:`7616`: ENH: Add RUR stationarity test to statsmodels.tsa.stattools
- :pr:`7617`: MAINT: Silence future warnings
- :pr:`7618`: DOC: spelling error in docs fixed
- :pr:`7619`: ENH: Improve ARDL and UECM
- :pr:`7620`: MAINT: Avoid passing bad optimization param
- :pr:`7641`: MAINT: Pin matplotlib
- :pr:`7643`: ENH: Improve error message in seasonal for bad freq
- :pr:`7644`: DOC: Update dev page flake8 command to follow PULL_REQUEST_TEMPLATE.md
- :pr:`7645`: ENH Fixed Range Unit-Root critical values
- :pr:`7648`: BUG/REF copula another round for 0.13
- :pr:`7649`: MAINT: Modernize prediction in notebooks
- :pr:`7651`: ENH: Improve copula notebook
- :pr:`7652`: MAINT: Temporarily change the default RNG in check_random_state
- :pr:`7656`: DOC: Add SARIMAX FAQ
- :pr:`7659`: DOC: Add to the SARIMAX FAQ
- :pr:`7661`: DOC: Improve SARIMAX FAQ Notebook
- :pr:`7662`: DOC: Improve ARIMA documentation
- :pr:`7668`: BUG: improve sidak multipletest precision close to zero
- :pr:`7669`: BUG: proportions_chisquare prevent integer overflow
- :pr:`7670`: BUG: ZI predict, fix offset default if None, allow exog_infl None if constant
- :pr:`7673`: ENH/BUG: graphics.plot_partregress add eval_env options
- :pr:`7676`: DOC: Remove duplication methods section
- :pr:`7677`: DOC: Second try ixing duplicate methods
- :pr:`7681`: fix a typo
- :pr:`7682`: ENH: McFadden and Cox&Snell Pseudo R squared to GLMResults
- :pr:`7685`: MAINT: Protect against changes in numeric indexes
- :pr:`7693`: ENH: add dk_params option to GLM info_criteria
- :pr:`7694`: ENH: quantile regression use dimension of x matrix rather than rank
- :pr:`7696`: ENH: add option for slim summary in OLS results
- :pr:`7697`: ENH add  tricube kernel
- :pr:`7698`: ENH: Fix lilliefors results for single-column DataFrames
- :pr:`7699`: DOC: Improve ARDL notebook
- :pr:`7701`: MAINT: Update TSA Api
- :pr:`7702`: DOC: Update versions.json
- :pr:`7703`: BUG: Factor fit ml em resets seed (rebased)
- :pr:`7704`: ENH: Enable VIF to work with DataFrames
- :pr:`7708`: MAINT: Update versions file
- :pr:`7709`: BUG: Correct ArmaProcess.from_estimation
- :pr:`7710`: BUG: describe / Description do not return percentiles
- :pr:`7713`: ENH: Oaxaca Variance/Other Models
- :pr:`7714`: DOC: Update release note
- :pr:`7721`: ENH: Added fft to ccovf and ccf
- :pr:`7723`: REF/ENH: more copula improvements for 0.13
- :pr:`7726`: DOC: Update release note
- :pr:`7727`: DOC: improve docs and docstrings, mainly for recent additions
- :pr:`7732`: DOC: api.py, docstring improvements
- :pr:`7735`: DOC: Correct MultivariateTestResults doc string
- :pr:`7737`: TST: Assert correct iloc dtypes
- :pr:`7738`: DOC: Correct MultivariateTestResults doc string
- :pr:`7739`: MAINT: Fix style issue
- :pr:`7740`: DOC: add missing function doc head
- :pr:`7742`: MAINT: Final issues in `__all__`
- :pr:`7743`: DOC: add to release notes, smaller doc fixes, references
- :pr:`7744`: MAINT: Fix hard to reach errors
- :pr:`7748`: BUG: fix summary().as_latex, line in top table dropped
- :pr:`7750`: ENH: Warn kwargs glm
- :pr:`7751`: REF: GLM init invalid kwargs use ValueWarning
- :pr:`7757`: BUG/MAINT/DOC: more 0.13
- :pr:`7766`: BUG: fix lowess spikes/nans from epsilon values
- :pr:`7768`: PERF/TST: Improve Lowess
- :pr:`7770`: DOC: Fix lowess notebook
- :pr:`7772`: ENH: add options to meta-analysis plot_forest
:orphan:

==============
Release 0.10.1
==============

Release summary
===============
This is a bug fix-only release

Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Chad Fulton
* Brock Mendel
* Peter Quackenbush
* Kerby Shedden
* Kevin Sheppard

and the general maintainer and code reviewer

* Josef Perktold

These lists of names are automatically generated based on git log, and may not
be complete.

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

* :pr:`5784`: MAINT: implement parts of #5220, deprecate ancient aliases
* :pr:`5892`: BUG: fix pandas compat
* :pr:`5893`: BUG: exponential smoothing - damped trend gives incorrect param, predictions
* :pr:`5895`: DOC: improvements to BayesMixedGLM docs, argument checking
* :pr:`5897`: MAINT: Use pytest.raises to check error message
* :pr:`5903`: BUG: Fix kwargs update bug in linear model fit_regularized
* :pr:`5917`: BUG: TVTP for Markov regression
* :pr:`5921`: BUG: Ensure exponential smoothers has continuous double data
* :pr:`5930`: BUG: Limit lags in KPSS
* :pr:`5933`: MAINT: Fix test that fails with positive probability
* :pr:`5935`: CLN: port parts of #5220
* :pr:`5940`: MAINT: Fix linting failures
* :pr:`5944`: BUG: Restore ResettableCache
* :pr:`5951`: BUG: Fix mosaic plot with missing category
* :pr:`5971`: BUG: Fix a future issue in ExpSmooth
.. _whatsnew_index:

=============
Release Notes
=============

.. toctree::
   :maxdepth: 1

   version0.13.0
   version0.12.1
   version0.12
   version0.11.1
   version0.11
   version0.10.2
   version0.10.1
   version0.10
   version0.9
   version0.8
   version0.7
   version0.6
   version0.5

For an overview of changes that occurred previous to the 0.5.0 release
see :ref:`old_changes`.

:orphan:

==============
Release 0.11.1
==============

Release summary
===============
This is a bug release.

Development summary and credits
===============================

Besides receiving contributions for new and improved features and for bugfixes,
important contributions to general maintenance for this release came from

* Kerby Shedden
* Josef Perktold
* Alex Lyttle
* Chad Fulton
* Kevin Sheppard
* Wouter De Coster

Merged Pull Requests
--------------------

The following Pull Requests were merged since the last release:

* :pr:`6433`: TST/BUG: use reset_randomstate
* :pr:`6438`: BUG: Change default optimizer for glm/ridge and make it user-settable
* :pr:`6453`: DOC: Fix the version that appears in the documentation
* :pr:`6456`: DOC: Send log to dev/null/
* :pr:`6461`: MAINT: correcting typo
* :pr:`6465`: MAINT: Avoid noise in f-pvalue
* :pr:`6469`: MAINT: Fix future warnings
* :pr:`6470`: BUG: fix tukey-hsd for 1 pvalue
* :pr:`6471`: MAINT: Fix issue with ragged array
* :pr:`6474`: BLD: Use pip on Azure
* :pr:`6515`: BUG: fix #6511
* :pr:`6520`: BUG: fix GAM for 1-dim exog_linear
* :pr:`6534`: MAINT: Relax tolerance on test that occasionally fails
* :pr:`6535`: MAINT: Restrict to Python 3.5+
=============
Release 0.5.0
=============

statsmodels 0.5 is a large and very exciting release that brings together a year of work done by 38 authors, including over 2000 commits. It contains many new features and a large amount of bug fixes detailed below.

See the :ref:`list of fixed issues <issues_list_05>` for specific closed issues.

The following major new features appear in this version.

Support for Model Formulas via Patsy
====================================

statsmodels now supports fitting models with a formula. This functionality is provided by `patsy <https://patsy.readthedocs.org/en/latest/>`_. Patsy is now a dependency for statsmodels. Models can be individually imported from the ``statsmodels.formula.api`` namespace or you can import them all as::

    import statsmodels.formula.api as smf

Alternatively, each model in the usual ``statsmodels.api`` namespace has a ``from_formula`` classmethod that will create a model using a formula. Formulas are also available for specifying linear hypothesis tests using the ``t_test`` and ``f_test`` methods after model fitting. A typical workflow can now look something like this.

.. code-block:: python

    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    url = 'https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/csv/HistData/Guerry.csv'
    data = pd.read_csv(url)

    # Fit regression model (using the natural log of one of the regressors)
    results = smf.ols('Lottery ~ Literacy + np.log(Pop1831)', data=data).fit()

See :ref:`here for some more documentation of using formulas in statsmodels <formula_examples>`

Empirical Likelihood (Google Summer of Code 2012 project)
---------------------------------------------------------

Empirical Likelihood-Based Inference for moments of univariate and multivariate variables is available as well as EL-based ANOVA tests. EL-based linear regression, including the regression through the origin model. In addition, the accelerated failure time model for inference on a linear regression model with a randomly right censored endogenous variable is available.

Analysis of Variance (ANOVA) Modeling
-------------------------------------

Support for ANOVA is now available including type I, II, and III sums of squares. See :ref:`anova`.

.. currentmodule:: statsmodels.nonparametric

Multivariate Kernel Density Estimators (GSoC 2012 project)
----------------------------------------------------------

Kernel density estimation has been extended to handle multivariate estimation as well via product kernels. It is available as :class:`sm.nonparametric.KDEMultivariate <kernel_density.KDEMultivariate>`. It supports least squares and maximum likelihood cross-validation for bandwidth estimation, as well as mixed continuous, ordered, and unordered categorical data. Conditional density estimation is also available via :class:`sm.nonparametric.KDEMUltivariateConditional <kernel_density.KDEMultivariateConditional>`.

nonparametric Regression (GSoC 2012 project)
---------------------------------------------
Kernel regression models are now available via :class:`sm.nonparametric.KernelReg <kernel_regression.KernelReg>`. It is based on the product kernel mentioned above, so it also has the same set of features including support for cross-validation as well as support for estimation mixed continuous and categorical variables. Censored kernel regression is also provided by `kernel_regression.KernelCensoredReg`.

Quantile Regression Model
-------------------------

.. currentmodule:: statsmodels.regression.quantile_regression

Quantile regression is supported via the :class:`sm.QuantReg <QuantReg>` class. Kernel and bandwidth selection options are available for estimating the asymptotic covariance matrix using a kernel density estimator.

Negative Binomial Regression Model
----------------------------------

.. currentmodule:: statsmodels.discrete.discrete_model

It is now possible to fit negative binomial models for count data via maximum-likelihood using the :class:`sm.NegativeBinomial <NegativeBinomial>` class. ``NB1``, ``NB2``, and ``geometric`` variance specifications are available.

l1-penalized Discrete Choice Models
-----------------------------------

A new optimization method has been added to the discrete models, which includes Logit, Probit, MNLogit and Poisson, that makes it possible to estimate the models with an l1, linear, penalization. This shrinks parameters towards zero and can set parameters that are not very different from zero to zero. This is especially useful if there are a large number of explanatory variables and a large associated number of parameters. `CVXOPT <https://cvxopt.org/>`_ is now an optional dependency that can be used for fitting these models.

New and Improved Graphics
-------------------------

.. currentmodule:: statsmodels.graphics

* **ProbPlot**: A new `ProbPlot` object has been added to provide a simple interface to create P-P, Q-Q, and probability plots with options to fit a distribution and show various reference lines. In the case of Q-Q and P-P plots, two different samples can be compared with the `other` keyword argument. :func:`sm.graphics.ProbPlot <gofplots.ProbPlot>`

.. code-block:: python
   
   import numpy as np
   import statsmodels.api as sm
   x = np.random.normal(loc=1.12, scale=0.25, size=37)
   y = np.random.normal(loc=0.75, scale=0.45, size=37)
   ppx = sm.ProbPlot(x)
   ppy =  sm.ProbPlot(y)
   fig1 = ppx.qqplot()
   fig2 = ppx.qqplot(other=ppy)

* **Mosaic Plot**: Create a mosaic plot from a contingency table. This allows you to visualize multivariate categorical data in a rigorous and informative way. Available with :func:`sm.graphics.mosaic <mosaicplot.mosaic>`.

* **Interaction Plot**: Interaction plots now handle categorical factors as well as other improvements. :func:`sm.graphics.interaction_plot <factorplots.interaction_plot>`.

* **Regression Plots**: The regression plots have been refactored and improved. They can now handle pandas objects and regression results instances appropriately. See :func:`sm.graphics.plot_fit <regressionplots.plot_fit>`, :func:`sm.graphics.plot_regress_exog <regressionplots.plot_regress_exog>`, :func:`sm.graphics.plot_partregress <regressionplots.plot_partregress>`, :func:`sm.graphics.plot_ccpr   <regressionplots.plot_ccpr>`, :func:`sm.graphics.abline_plot <regressionplots.abline_plot>`, :func:`sm.graphics.influence_plot <regressionplots.influence_plot>`, and :func:`sm.graphics.plot_leverage_resid2 <regressionplots.plot_leverage_resid2>`.

.. currentmodule:: statsmodels.stats.power

Power and Sample Size Calculations
----------------------------------

The power module (``statsmodels.stats.power``) currently implements power and sample size calculations for the t-tests (:class:`sm.stats.TTestPower <TTestPower>`, :class:`sm.stats.TTestIndPower <TTestIndPower>`), normal based test (:class:`sm.stats.NormIndPower <NormIndPower>`), F-tests (:class:`sm.stats.FTestPower <FTestPower>`, `:class:sm.stats.FTestAnovaPower <FTestAnovaPower>`) and Chisquare goodness of fit (:class:`sm.stats.GofChisquarePower <GofChisquarePower>`) test. The implementation is class based, but the module also provides three shortcut functions, :func:`sm.stats.tt_solve_power <tt_solve_power>`, :func:`sm.stats.tt_ind_solve_power <tt_ind_solve_power>` and :func:`sm.stats.zt_ind_solve_power <zt_ind_solve_power>` to solve for any one of the parameters of the power equations. See this `blog post <http://jpktd.blogspot.fr/2013/03/statistical-power-in-statsmodels.html>`_ for a more in-depth description of the additions.


Other important new features
----------------------------
* **IPython notebook examples**: Many of our examples have been converted or added as IPython notebooks now. They are available `here <https://www.statsmodels.org/devel/examples/index.html#notebook-examples>`_.

* **Improved marginal effects for discrete choice models**: Expanded options for obtaining marginal effects after the estimation of nonlinear discrete choice models are available. See :py:meth:`get_margeff <statsmodels.discrete.discrete_model.DiscreteResuls.get_margeff>`.

* **OLS influence outlier measures**: After the estimation of a model with OLS, the common set of influence and outlier measures and a outlier test are now available attached as methods ``get_influnce`` and ``outlier_test`` to the Results instance. See :py:class:`OLSInfluence <statsmodels.stats.outliers_influence.OLSInfluence>` and :func:`outlier_test <statsmodels.stats.outliers_influence.outlier_test>`.

* **New datasets**: New :ref:`datasets <datasets>` are available for examples.

* **Access to R datasets**: We now have access to many of the same datasets available to R users through the `Rdatasets project <https://vincentarelbundock.github.io/Rdatasets/>`_. You can access these using the :func:`sm.datasets.get_rdataset <statsmodels.datasets.get_rdataset>` function. This function also includes caching of these datasets.

* **Improved numerical differentiation tools**: Numerical differentiation routines have been greatly improved and expanded to cover all the routines discussed in::

    Ridout, M.S. (2009) Statistical applications of the complex-step method
        of numerical differentiation. The American Statistician, 63, 66-74

  See the :ref:`sm.tools.numdiff <numdiff>` module.

* **Consistent constant handling across models**: Result statistics no longer rely on the assumption that a constant is present in the model.

* **Missing value handling across models**: Users can now control what models do in the presence of missing values via the ``missing`` keyword available in the instantiation of every model. The options are ``'none'``, ``'drop'``, and ``'raise'``. The default is ``'none'``, which does no missing value checks. To drop missing values use ``'drop'``. And ``'raise'`` will raise an error in the presence of any missing data.

.. currentmodule:: statsmodels.iolib

* **Ability to write Stata datasets**: Added the ability to write Stata ``.dta`` files. See :class:`sm.iolib.StataWriter <foreign.StataWriter>`.

.. currentmodule:: statsmodels.tsa.arima.model

* **ARIMA modeling**: statsmodels now has support for fitting Autoregressive Integrated Moving Average (ARIMA) models. See :class:`ARIMA` and :class:`ARIMAResults` for more information.

* **Support for dynamic prediction in AR(I)MA models**: It is now possible to obtain dynamic in-sample forecast values in :class:`ARMA` and :class:`ARIMA` models.

* **Improved Pandas integration**: statsmodels now supports all frequencies available in pandas for time-series modeling. These are used for intelligent dates handling for prediction. These features are available, if you pass a pandas Series or DataFrame with a DatetimeIndex to a time-series model.

.. currentmodule:: statsmodels

* **New statistical hypothesis tests**: Added statistics for calculating interrater agreement including Cohen's kappa and Fleiss' kappa (See :ref:`interrater`), statistics and hypothesis tests for proportions (See :ref:`proportion stats <proportion_stats>`), Tukey HSD (with plot) was added as an enhancement to the multiple comparison tests (:class:`sm.stats.multicomp.MultiComparison <sandbox.stats.multicomp.MultiComparison>`, :func:`sm.stats.multicomp.pairwise_tukeyhsd <stats.multicomp.pairwise_tukeyhsd>`). Weighted statistics and t tests were enhanced with new options. Tests of equivalence for one sample and two independent or paired samples were added based on t tests and z tests (See :ref:`tost`).  


Major Bugs fixed
----------------

* Post-estimation statistics for weighted least squares that depended on the centered total sum of squares were not correct. These are now correct and tested. See :issue:`501`.

* Regression through the origin models now correctly use uncentered total sum of squares in post-estimation statistics. This affected the :math:`R^2` value in linear models without a constant. See :issue:`27`.

Backwards incompatible changes and deprecations
-----------------------------------------------

* Cython code is now non-optional. You will need a C compiler to build from source. If building from github and not a source release, you will also need Cython installed. See the :ref:`installation documentation <install>`.

* The ``q_matrix`` keyword to `t_test` and `f_test` for linear models is deprecated. You can now specify linear hypotheses using formulas.

.. currentmodule:: statsmodels.tsa

* The ``conf_int`` keyword to :func:`sm.tsa.acf <stattools.acf>` is deprecated.

* The ``names`` argument is deprecated in :class:`sm.tsa.VAR <vector_ar.var_model.VAR>` and `sm.tsa.SVAR <vector_ar.svar_model.SVAR>`. This is now automatically detected and handled.

.. currentmodule:: statsmodels.tsa

* The ``order`` keyword to :py:meth:`sm.tsa.ARMA.fit <ARMA.fit>` is deprecated. It is now passed in during model instantiation.

.. currentmodule:: statsmodels.distributions

* The empirical distribution function (:class:`sm.distributions.ECDF <ECDF>`) and supporting functions have been moved to ``statsmodels.distributions``. Their old paths have been deprecated.

* The ``margeff`` method of the discrete choice models has been deprecated. Use ``get_margeff`` instead. See above. Also, the vague ``resid`` attribute of the discrete choice models has been deprecated in favor of the more descriptive ``resid_dev`` to indicate that they are deviance residuals.

.. currentmodule:: statsmodels.nonparametric.kde

* The class ``KDE`` has been deprecated and renamed to :class:`KDEUnivariate` to distinguish it from the new ``KDEMultivariate``. See above.

Development summary and credits
-------------------------------

The previous version (statsmodels 0.4.3) was released on July 2, 2012. Since then we have closed a total of 380 issues, 172 pull requests and 208 regular issues. The :ref:`detailed list<issues_list_05>` can be viewed.

This release is a result of the work of the following 38 authors who contributed total of 2032 commits. If for any reason, we have failed to list your name in the below, please contact us:

* Ana Martinez Pardo <anamartinezpardo-at-gmail.com>
* anov <novikova.go.zoom-at-gmail.com>
* avishaylivne <avishay.livne-at-gmail.com>
* Bruno Rodrigues <rodrigues.bruno-at-aquitania.org>
* Carl Vogel <carljv-at-gmail.com>
* Chad Fulton <chad-at-chadfulton.com>
* Christian Prinoth <christian-at-prinoth.name>
* Daniel B. Smith <neuromathdan-at-gmail.com>
* dengemann <denis.engemann-at-gmail.com>
* Dieter Vandenbussche <dvandenbussche-at-axioma.com>
* Dougal Sutherland <dougal-at-gmail.com>
* Enrico Giampieri <enrico.giampieri-at-unibo.it>
* evelynmitchell <efm-github-at-linsomniac.com>
* George Panterov <econgpanterov-at-gmail.com>
* Grayson <graysonbadgley-at-gmail.com>
* Jan Schulz <jasc-at-gmx.net>
* Josef Perktold <josef.pktd-at-gmail.com>
* Jeff Reback <jeff-at-reback.net>
* Justin Grana <jg3705a-at-student.american.edu>
* langmore <ianlangmore-at-gmail.com>
* Matthew Brett <matthew.brett-at-gmail.com>
* Nathaniel J. Smith <njs-at-pobox.com>
* otterb <itoi-at-live.com>
* padarn <padarn-at-wilsonp.anu.edu.au>
* Paul Hobson <pmhobson-at-gmail.com>
* Pietro Battiston <me-at-pietrobattiston.it>
* Ralf Gommers <ralf.gommers-at-googlemail.com>
* Richard T. Guy <richardtguy84-at-gmail.com>
* Robert Cimrman <cimrman3-at-ntc.zcu.cz>
* Skipper Seabold <jsseabold-at-gmail.com>
* Thomas Haslwanter <thomas.haslwanter-at-fh-linz.at>
* timmie <timmichelsen-at-gmx-topmail.de>
* Tom Augspurger <thomas-augspurger-at-uiowa.edu>
* Trent Hauck <trent.hauck-at-gmail.com>
* tylerhartley <tyleha-at-gmail.com>
* Vincent Arel-Bundock <varel-at-umich.edu>
* VirgileFritsch <virgile.fritsch-at-gmail.com>
* Zhenya <evgeni-at-burovski.me>

.. note:: 

   Obtained by running ``git log v0.4.3..HEAD --format='* %aN <%aE>' | sed 's/@/\-at\-/' | sed 's/<>//' | sort -u``.

.. _issues_list_05:

Issues closed in the 0.5.0 development cycle
============================================

Issued closed in 0.5.0
-----------------------

GitHub stats for release 0.5.0 (07/02/2012/ - 08/14/2013/).

We closed a total of 380 issues, 172 pull requests and 208 regular issues. This is the full list (generated with the script  :file:`tools/github_stats.py`):

This list is automatically generated, and may be incomplete:

Pull Requests (172):

* :pr:`1015`: DOC: Bump version. Remove done tasks.
* :pr:`1010`: DOC/RLS: Update release notes workflow. Help Needed!
* :pr:`1014`: DOC: nbgenerate does not like the comment at end of line.
* :pr:`1012`: DOC: Add link to notebook and crosslink ref. Closes #924.
* :pr:`997`: misc, tests, diagnostic
* :pr:`1009`: MAINT: Add .mailmap file.
* :pr:`817`: Add 3 new unit tests for arima_process
* :pr:`1001`: BUG include_package_data for install closes #907
* :pr:`1005`: GITHUB: Contributing guidelines
* :pr:`1007`: Cleanup docs for release
* :pr:`1003`: BUG: Workaround for bug in sphinx 1.1.3. See #1002.
* :pr:`1004`: DOC: Update maintainer notes with branching instructions.
* :pr:`1000`: BUG: Support pandas 0.8.0.
* :pr:`996`: BUG: Handle combo of pandas 0.8.0 and dateutils 1.5.0
* :pr:`995`: ENH: Print dateutil version.
* :pr:`994`: ENH: Fail gracefully for version not found.
* :pr:`993`: More conservative error catching in TimeSeriesModel
* :pr:`992`: Misc fixes 12: adjustments to unit test
* :pr:`985`: MAINT: Print versions script.
* :pr:`986`: ENH: Prefer to_offset to get_offset. Closes #964.
* :pr:`984`: COMPAT: Pandas 0.8.1 compatibility. Closes #983.
* :pr:`982`: Misc fixes 11
* :pr:`978`: TST: generic mle pareto disable bsejac tests with estimated loc
* :pr:`977`: BUG python 3.3 fix for numpy str TypeError, see #633
* :pr:`975`: Misc fixes 10 numdiff
* :pr:`970`: BUG: array too long, raises exception with newer numpy closes #967
* :pr:`965`: Vincent summary2 rebased
* :pr:`933`: Update and improve GenericlikelihoodModel and miscmodels
* :pr:`950`: BUG/REF mcnemar fix exact pvalue, allow table as input
* :pr:`951`: Pylint emplike formula genmod
* :pr:`956`: Fix a docstring in KDEMultivariateConditional.
* :pr:`949`: BUG fix lowess sort when nans closes #946
* :pr:`932`: ENH: support basinhopping solver in LikelihoodModel.fit()
* :pr:`927`: DOC: clearer minimal example
* :pr:`919`: OLS summary crash
* :pr:`918`: Fixes10 emplike lowess
* :pr:`909`: Bugs in GLM pvalues, more tests, pylint
* :pr:`906`: ENH: No fmax with Windows SDK so define inline.
* :pr:`905`: MAINT more fixes
* :pr:`898`: Misc fixes 7
* :pr:`896`: Quantreg rebase2
* :pr:`895`: Fixes issue #832
* :pr:`893`: ENH: Remove unneeded restriction on low. Closes #867.
* :pr:`894`: MAINT: Remove broken function. Keep deprecation. Closes #781.
* :pr:`856`: Carljv improved lowess rebased2
* :pr:`884`: Pyflakes cleanup
* :pr:`887`: BUG: Fix kde caching
* :pr:`883`: Fixed pyflakes issue in discrete module
* :pr:`882`: Update predstd.py
* :pr:`871`: Update of sandbox doc
* :pr:`631`: WIP: Correlation positive semi definite
* :pr:`857`: BLD: apt get dependencies from Neurodebian, whitespace cleanup
* :pr:`855`: AnaMP issue 783 mixture rvs tests rebased
* :pr:`854`: Enrico multinear rebased
* :pr:`849`: Tyler tukeyhsd rebased
* :pr:`848`: BLD TravisCI use python-dateutil package
* :pr:`784`: Misc07 cleanup multipletesting and proportions
* :pr:`841`: ENH: Add load function to main API. Closes #840.
* :pr:`820`: Ensure that tuples are not considered as data, not as data containers
* :pr:`822`: DOC: Update for Cython changes.
* :pr:`765`: Fix build issues
* :pr:`800`: Automatically generate output from notebooks
* :pr:`802`: BUG: Use two- not one-sided t-test in t_test. Closes #740.
* :pr:`806`: ENH: Import formula.api in statsmodels.api namespace.
* :pr:`803`: ENH: Fix arima error message for bad start_params
* :pr:`801`: DOC: Fix ANOVA section titles
* :pr:`795`: Negative Binomial Rebased
* :pr:`787`: Origintests
* :pr:`794`: ENH: Allow pandas-in/pandas-out in tsa.filters
* :pr:`791`: Github stats for release notes
* :pr:`779`: added np.asarray call to durbin_watson in stattools
* :pr:`772`: Anova docs
* :pr:`776`: BUG: Fix dates_from_range with length. Closes #775.
* :pr:`774`: BUG: Attach prediction start date in AR. Closes #773.
* :pr:`767`: MAINT: Remove use of deprecated from examples and docs.
* :pr:`762`: ENH: Add new residuals to wrapper
* :pr:`754`: Fix arima predict
* :pr:`760`: ENH: Adjust for k_trend in information criteria. Closes #324.
* :pr:`761`: ENH: Fixes and tests sign_test. Closes #642.
* :pr:`759`: Fix 236
* :pr:`758`: DOC: Update VAR docs. Closes #537.
* :pr:`752`: Discrete cleanup
* :pr:`750`: VAR with 1d array
* :pr:`748`: Remove reference to new_t_test and new_f_test.
* :pr:`739`: DOC: Remove outdated note in docstring
* :pr:`732`: BLD: Check for patsy dependency at build time + docs
* :pr:`731`: Handle wrapped
* :pr:`730`: Fix opt fulloutput
* :pr:`729`: Get rid of warnings in docs build
* :pr:`698`: update url for hsb2 dataset
* :pr:`727`: DOC: Fix indent and add missing params to linear models. Closes #709.
* :pr:`726`: CLN: Remove unused method. Closes #694
* :pr:`725`: BUG: Should call anova_single. Closes #702.
* :pr:`723`: Rootfinding for Power
* :pr:`722`: Handle pandas.Series with names in make_lags
* :pr:`714`: Fix 712
* :pr:`668`: Allow for any pandas frequency to be used in TimeSeriesModel.
* :pr:`711`: Misc06 - bug fixes
* :pr:`708`: BUG: Fix one regressor case for conf_int. Closes #706.
* :pr:`700`: Bugs rebased
* :pr:`680`: BUG: Swap arguments in fftconvolve for scipy >= 0.12.0
* :pr:`640`: Misc fixes 05
* :pr:`663`: a typo in runs.py doc string for mcnemar test
* :pr:`652`: WIP: fixing pyflakes / pep8, trying to improve readability
* :pr:`619`: DOC: intro to formulas
* :pr:`648`: BF: Make RLM stick to Huber's description
* :pr:`649`: Bug Fix
* :pr:`637`: Pyflakes cleanup
* :pr:`634`: VAR DOC typo
* :pr:`623`: Slowtests
* :pr:`621`: MAINT: in setup.py, only catch ImportError for pandas.
* :pr:`590`: Cleanup test output
* :pr:`591`: Interrater agreement and reliability measures
* :pr:`618`: Docs fix the main warnings and errors during sphinx build
* :pr:`610`: nonparametric examples and some fixes
* :pr:`578`: Fix 577
* :pr:`575`: MNT: Remove deprecated scikits namespace
* :pr:`499`: WIP: Handle constant
* :pr:`567`: Remove deprecated
* :pr:`571`: Dataset docs
* :pr:`561`: Grab rdatasets
* :pr:`570`: DOC: Fixed links to Rdatasets
* :pr:`524`: DOC: Clean up discrete model documentation.
* :pr:`506`: ENH: Re-use effects if model fit with QR
* :pr:`556`: WIP:  L1 doc fix
* :pr:`564`: TST: Use native integer to avoid issues in dtype asserts
* :pr:`543`: Travis CI using M.Brett nipy hack
* :pr:`558`: Plot cleanup
* :pr:`541`: Replace pandas DataMatrix with DataFrame
* :pr:`534`: Stata test fixes
* :pr:`532`: Compat 323
* :pr:`531`: DOC: Add ECDF to distributions docs
* :pr:`526`: ENH: Add class to write Stata binary dta files
* :pr:`521`: DOC: Add abline plot to docs
* :pr:`518`: Small fixes: interaction_plot
* :pr:`508`: ENH: Avoid taking cholesky decomposition of diagonal matrix
* :pr:`509`: DOC: Add ARIMA to docs
* :pr:`510`: DOC: realdpi is disposable personal income. Closes #394.
* :pr:`507`: ENH: Protect numdifftools import. Closes #45
* :pr:`504`: Fix weights
* :pr:`498`: DOC: Add patys requirement to install docs
* :pr:`491`: Make _data a public attribute.
* :pr:`494`: DOC: Fix pandas links
* :pr:`492`: added intersphinx for pandas
* :pr:`422`: Handle missing data
* :pr:`485`: ENH: Improve error message for pandas objects without dates in index
* :pr:`428`: Remove other data
* :pr:`483`: Arima predict bug
* :pr:`482`: TST: Do array-array comparison when using numpy.testing
* :pr:`471`: Formula rename df -> data
* :pr:`473`: Vincent docs tweak rebased
* :pr:`468`: Docs 050
* :pr:`462`: El aft rebased
* :pr:`461`: TST: numpy 1.5.1 compatibility
* :pr:`460`: Emplike desc reg rebase
* :pr:`410`: Discrete model marginal effects
* :pr:`417`: Numdiff cleanup
* :pr:`398`: Improved plot_corr and plot_corr_grid functions.
* :pr:`401`: BUG: Finish refactoring margeff for dummy. Closes #399.
* :pr:`400`: MAINT: remove lowess.py, which was kept in 0.4.x for backwards compatibi...
* :pr:`371`: BF+TEST: fixes, checks and tests for isestimable
* :pr:`351`: ENH: Copy diagonal before write for upcoming numpy changes
* :pr:`384`: REF: Move mixture_rvs out of sandbox.
* :pr:`368`: ENH: Add polished version of acf/pacf plots with confidence intervals
* :pr:`378`: Infer freq
* :pr:`374`: ENH: Add Fair's extramarital affair dataset. From tobit-model branch.
* :pr:`358`: ENH: Add method to OLSResults for outlier detection
* :pr:`369`: ENH: allow predict to pass through patsy for transforms
* :pr:`352`: Formula integration rebased
* :pr:`360`: REF: Deprecate order in fit and move to ARMA init
* :pr:`366`: Version fixes
* :pr:`359`: DOC: Fix sphinx warnings

Issues (208):

* :issue:`1036`: Series no longer inherits from ndarray
* :issue:`1038`: DataFrame with integer names not handled in ARIMA
* :issue:`1028`: Test fail with windows and Anaconda - Low priority
* :issue:`676`: acorr_breush_godfrey  undefined nlags
* :issue:`922`: lowess returns inconsistent with option
* :issue:`425`: no bse in robust with norm=TrimmedMean
* :issue:`1025`: add_constant incorrectly detects constant column
* :issue:`533`: py3 compatibility ``pandas.read_csv(urlopen(...))``
* :issue:`662`: doc: install instruction: explicit about removing scikits.statsmodels
* :issue:`910`: test failure Ubuntu TestARMLEConstant.test_dynamic_predict
* :issue:`80`: t_model: f_test, t_test do not work
* :issue:`432`: GenericLikelihoodModel change default for score and hessian
* :issue:`454`: BUG/ENH: HuberScale instance is not used, allow user defined scale estimator
* :issue:`98`: check connection or connect summary to variable names in wrappers
* :issue:`418`: BUG: MNLogit loglikeobs, jac
* :issue:`1017`: nosetests warnings
* :issue:`924`: DOCS link in notebooks to notebook for download
* :issue:`1011`: power ttest endless loop possible
* :issue:`907`: BLD data_files for stats.libqsturng
* :issue:`328`: consider moving example scripts into IPython notebooks
* :issue:`1002`: Docs will not build with Sphinx 1.1.3
* :issue:`69`: Make methods like compare_ftest work with wrappers
* :issue:`503`: summary_old in RegressionResults
* :issue:`991`: TST precision of normal_power
* :issue:`945`: Installing statsmodels from github?
* :issue:`964`: Prefer to_offset not get_offset in tsa stuff
* :issue:`983`: bug: pandas 0.8.1 incompatibility
* :issue:`899`: build_ext inplace does not cythonize
* :issue:`923`: location of initialization code
* :issue:`980`: auto lag selection in  S_hac_simple
* :issue:`968`: genericMLE Ubuntu test failure
* :issue:`633`: python 3.3 compatibility
* :issue:`728`: test failure for solve_power with fsolve
* :issue:`971`: numdiff test cases
* :issue:`976`: VAR Model does not work in 1D
* :issue:`972`: numdiff: epsilon has no minimum value
* :issue:`967`: lowes test failure Ubuntu
* :issue:`948`: nonparametric tests: mcnemar, cochranq unit test
* :issue:`963`: BUG in runstest_2sample
* :issue:`946`: Issue with lowess() smoother in statsmodels
* :issue:`868`: k_vars > nobs
* :issue:`917`: emplike emplikeAFT stray dimensions
* :issue:`264`: version comparisons need to be made more robust (may be just use LooseVersion)
* :issue:`674`: failure in test_foreign, pandas testing
* :issue:`828`: GLMResults inconsistent distribution in pvalues
* :issue:`908`: RLM missing test for tvalues, pvalues
* :issue:`463`: formulas missing in docs
* :issue:`256`: discrete Nbin has zero test coverage
* :issue:`831`: test errors running bdist
* :issue:`733`: Docs: interrater cohens_kappa is missing
* :issue:`897`: lowess failure - sometimes
* :issue:`902`: test failure tsa.filters  precision too high
* :issue:`901`: test failure stata_writer_pandas, newer versions of pandas
* :issue:`900`: ARIMA.__new__   errors on python 3.3
* :issue:`832`: notebook errors
* :issue:`867`: Baxter King has unneeded limit on value for low?
* :issue:`781`: discreteResults margeff method not tests, obsolete
* :issue:`870`: discrete unit tests duplicates
* :issue:`630`: problems in regression plots
* :issue:`885`: Caching behavior for KDEUnivariate icdf
* :issue:`869`: sm.tsa.ARMA(..., order=(p,q)) gives "__init__() got an unexpected keyword argument 'order'" error
* :issue:`783`: statsmodels.distributions.mixture_rvs.py    no unit tests
* :issue:`824`: Multicomparison w/Pandas Series
* :issue:`789`: presentation of multiple comparison results
* :issue:`764`: BUG: multipletests incorrect reject for Holm-Sidak
* :issue:`766`: multipletests - status and tests of 2step FDR procedures
* :issue:`763`: Bug: multipletests raises exception with empty array
* :issue:`840`: sm.load should be in the main API namespace
* :issue:`830`: invalid version number
* :issue:`821`: Fail gracefully when extensions are not built
* :issue:`204`: Cython extensions built twice?
* :issue:`689`: tutorial notebooks
* :issue:`740`: why does t_test return one-sided p-value
* :issue:`804`: What goes in statsmodels.formula.api?
* :issue:`675`: Improve error message for ARMA SVD convergence failure.
* :issue:`15`: arma singular matrix
* :issue:`559`: Add Rdatasets to optional dependencies list
* :issue:`796`: Prediction Standard Errors
* :issue:`793`: filters are not pandas aware
* :issue:`785`: Negative R-squared
* :issue:`777`: OLS residuals returned as Pandas series when endog and exog are Pandas series
* :issue:`770`: Add ANOVA to docs
* :issue:`775`: Bug in dates_from_range
* :issue:`773`: AR model pvalues error with Pandas
* :issue:`768`: multipletests: numerical problems at threshold
* :issue:`355`: add draw if interactive to plotting functions
* :issue:`625`: Exog is not correctly handled in ARIMA predict
* :issue:`626`: ARIMA summary does not print exogenous variable coefficients
* :issue:`657`: order (0,1) breaks ARMA forecast
* :issue:`736`: ARIMA predict problem for ARMA model
* :issue:`324`: ic in ARResults, aic, bic, hqic, fpe inconsistent definition?
* :issue:`642`: sign_test   check
* :issue:`236`: AR start_params broken
* :issue:`235`: tests hang on Windows
* :issue:`156`: matplotlib deprecated legend ? var plots
* :issue:`331`: Remove stale tests
* :issue:`592`: test failures in datetools
* :issue:`537`: Var Models
* :issue:`755`: Unable to access AR fit parameters when model is estimated with pandas.DataFrame
* :issue:`670`: discrete: numerically useless clipping
* :issue:`515`: MNLogit residuals raise a TypeError
* :issue:`225`: discrete models only define deviance residuals
* :issue:`594`: remove skiptest in TestProbitCG
* :issue:`681`: Dimension Error in discrete_model.py When Running test_dummy_*
* :issue:`744`: DOC: new_f_test
* :issue:`549`: Ship released patsy source in statsmodels
* :issue:`588`: patsy is a hard dependency?
* :issue:`716`: Tests missing for functions if pandas is used
* :issue:`715`: statsmodels regression plots not working with pandas datatypes
* :issue:`450`: BUG: full_output in optimizers Likelihood model
* :issue:`709`: DOCstrings linear models do not have missing params
* :issue:`370`: BUG weightstats has wrong cov
* :issue:`694`: DiscreteMargins duplicate method
* :issue:`702`: bug, pylint stats.anova
* :issue:`423`: Handling of constant across models
* :issue:`456`: BUG: ARMA date handling incompatibility with recent pandas
* :issue:`514`: NaNs in Multinomial
* :issue:`405`: Check for existing old version of scikits.statsmodels?
* :issue:`586`: Segmentation fault with OLS
* :issue:`721`: Unable to run AR on named time series objects
* :issue:`125`: caching pinv_wexog broke iterative fit - GLSAR
* :issue:`712`: TSA bug with frequency inference
* :issue:`319`: Timeseries Frequencies
* :issue:`707`: .summary with alpha ignores parsed value
* :issue:`673`: nonparametric: bug in _kernel_base
* :issue:`710`: test_power failures
* :issue:`706`: .conf_int() fails on linear regression without intercept
* :issue:`679`: Test Baxter King band-pass filter fails with scipy 0.12 beta1
* :issue:`552`: influence outliers breaks when regressing on constant
* :issue:`639`: test folders not on python path
* :issue:`565`: omni_normtest does not propagate the axis argument
* :issue:`563`: error in doc generation for AR.fit
* :issue:`109`: TestProbitCG failure on Ubuntu
* :issue:`661`: from scipy import comb fails on the latest scipy 0.11.0
* :issue:`413`: DOC: example_discrete.py missing from 0.5 documentation
* :issue:`644`: FIX: factor plot + examples broken
* :issue:`645`: STY: pep8 violations in many examples
* :issue:`173`: doc sphinx warnings
* :issue:`601`: bspline.py dependency on old scipy.stats.models
* :issue:`103`: ecdf and step function conventions
* :issue:`18`: Newey-West sandwich covariance is missing
* :issue:`279`: cov_nw_panel not tests, example broken
* :issue:`150`: precision in test_discrete.TestPoissonNewton.test_jac ?
* :issue:`480`: rescale loglike for optimization
* :issue:`627`: Travis-CI support for scipy
* :issue:`622`: mark tests as slow in emplike
* :issue:`589`: OLS F-statistic error
* :issue:`572`: statsmodels/tools/data.py Stuck looking for la.py
* :issue:`580`: test errors in graphics
* :issue:`577`: PatsyData detection buglet
* :issue:`470`: remove deprecated features
* :issue:`573`: lazy imports are (possibly) very slow
* :issue:`438`: New results instances are not in online documentation
* :issue:`542`: Regression plots fail when Series objects passed to sm.OLS
* :issue:`239`: release 0.4.x
* :issue:`530`: l1 docs issues
* :issue:`539`: test for statawriter (failure)
* :issue:`490`: Travis CI on PRs
* :issue:`252`: doc: distributions.rst refers to sandbox only
* :issue:`85`: release 0.4
* :issue:`65`: MLE fit of AR model has no tests
* :issue:`522`: ``test`` does not propagate arguments to nose
* :issue:`517`: missing array conversion or shape in linear model
* :issue:`523`: test failure with ubuntu decimals too large
* :issue:`520`: web site documentation, source not updated
* :issue:`488`: Avoid cholesky decomposition of diagonal matrices in linear regression models
* :issue:`394`: Definition in macrodata NOTE
* :issue:`45`: numdifftools dependency
* :issue:`501`: WLS/GLS post estimation results
* :issue:`500`: WLS fails if weights is a pandas.Series
* :issue:`27`: add hasconstant indicator for R-squared and df calculations
* :issue:`497`: DOC: add patsy?
* :issue:`495`: ENH: add footer SimpleTable
* :issue:`402`: model._data -> model.data?
* :issue:`477`: VAR NaN Bug
* :issue:`421`: Enhancement: Handle Missing Data
* :issue:`489`: Expose model._data as model.data
* :issue:`315`: tsa models assume pandas object indices are dates
* :issue:`440`: arima predict is broken for steps > q and q != 1
* :issue:`458`: TST BUG?   comparing pandas and array in tests, formula
* :issue:`464`: from_formula signature
* :issue:`245`: examples in docs: make nicer
* :issue:`466`: broken example, pandas
* :issue:`57`: Unhelpful error from bad exog matrix in model.py
* :issue:`271`: ARMA.geterrors requires model to be fit
* :issue:`350`: Writing to array returned np.diag
* :issue:`354`: example_rst does not copy unchanged files over
* :issue:`467`: Install issues with Pandas
* :issue:`444`: ARMA example on stable release website not working
* :issue:`377`: marginal effects count and discrete adjustments
* :issue:`426`: "svd" method not supported for OLS.fit()
* :issue:`409`: Move numdiff out of the sandbox
* :issue:`416`: Switch to complex-step Hessian for AR(I)MA
* :issue:`415`: bug in kalman_loglike_complex
* :issue:`397`: plot_corr axis text labeling not working (with fix)
* :issue:`399`: discrete errors due to incorrect in-place operation
* :issue:`389`: VAR test_normality is broken with KeyError
* :issue:`388`: Add tsaplots to graphics.api as graphics.tsa
* :issue:`387`: predict date was not getting set with start = None
* :issue:`386`: p-values not returned from acf
* :issue:`385`: Allow AR.select_order to work without model being fit
* :issue:`383`: Move mixture_rvs out of sandbox.
* :issue:`248`: ARMA breaks with a 1d exog
* :issue:`273`: When to give order for AR/AR(I)MA
* :issue:`363`: examples folder -> tutorials folder
* :issue:`346`: docs in sitepackages
* :issue:`353`: PACF docs raise a sphinx warning
* :issue:`348`: python 3.2.3 test failure zip_longest
=======================
Exceptions and Warnings
=======================

Exceptions
----------

Errors derive from Exception or another custom error. Custom errors are
only needed if standard errors, for example ValueError or TypeError, are not
accurate descriptions of the cause for the error.

.. module:: statsmodels.tools.sm_exceptions
   :synopsis: Exceptions and Warnings

.. currentmodule:: statsmodels.tools.sm_exceptions

.. autosummary::
   :toctree: generated/

   ParseError
   PerfectSeparationError
   X13NotFoundError
   X13Error

Warnings
--------

Warnings derive from either an existing warning or another custom
warning, and are often accompanied by a string using the format
``warning_name_doc`` that services as a generic message to use when the
warning is raised.


.. currentmodule:: statsmodels.tools.sm_exceptions

.. autosummary::
   :toctree: generated/

   X13Warning
   IOWarning
   ModuleUnavailableWarning
   ModelWarning
   ConvergenceWarning
   CacheWriteWarning
   IterationLimitWarning
   InvalidTestWarning
   NotImplementedWarning
   OutputWarning
   DomainWarning
   ValueWarning
   EstimationWarning
   SingularMatrixWarning
   HypothesisTestWarning
   InterpolationWarning
   PrecisionWarning
   SpecificationWarning
   HessianInversionWarning
   CollinearityWarning
Get Involved
============

Where to Start?
---------------

Use grep or download a tool like `grin <https://pypi.python.org/pypi/grin>`__
to search the code for TODO notes::

    grin -i -I "*.py*" todo

This shows almost 700 TODOs in the code base right now. Feel free to inquire on
the `mailing list <https://groups.google.com/forum/#!forum/pystatsmodels>`_
about any of these.

Sandbox
-------

We currently have a large amount code in the :ref:`sandbox`. The medium term
goal is to move much of this to feature branches as it gets worked on and
remove the sandbox folder. Many of these models and functions are close to
done, however, and we welcome any and all contributions to complete them,
including refactoring, documentation, and tests. These models include
generalized additive models (GAM), information theoretic models such as
maximum entropy, survival models, systems of equation models, restricted least
squares, panel data models, and time series models such as (G)ARCH.

.. .. toctree::
..   :maxdepth: 4
..
..   ../sandbox

Contribute an Example
---------------------

Contribute an :ref:`example <examples>`, add some technical documentation, or
contribute a statistics tutorial.
.. _examples:

Examples
========

Examples are invaluable for new users who hope to get up and running quickly
with `statsmodels`, and they are extremely useful to those who wish to explore
new features of `statsmodels`. We hope to provide documentation and tutorials
for as many models and use-cases as possible! Please consider submitting an
example with any PR that introduces new functionality.

User-contributed examples/tutorials/recipes can be placed on the
`statsmodels examples wiki page <https://github.com/statsmodels/statsmodels/wiki/Examples>`_
That wiki page is freely editable. Please post your cool tricks,
examples, and recipes on there!

If you would rather have your example file officially accepted to the
`statsmodels` distribution and posted on this website, you will need to go
through the normal `patch submission process <index.html#submitting-a-patch>`_
and follow the instructions that follow.

File Format
-----------

Examples are best contributed as IPython notebooks. Save your notebook with all
output cells cleared in ``examples/notebooks``. From the notebook save the pure
Python output to ``examples/python``. The first line of the Notebook *must* be
a header cell that contains a title for the notebook, if you want the notebook
to be included in the documentation.


The Example Gallery
-------------------

We have a gallery of example notebooks available
`here <https://www.statsmodels.org/devel/examples/index.html>`_. If you would
like your example to show up in this gallery, add a link to the notebook in
``docs/source/examples/landing.json``. For the thumbnail, take a screenshot of
what you think is the best "hook" for the notebook. The image will be displayed
at 360 x 225 (W x H). It's best to save the image as a PNG with a resolution
that is some multiple of 360 x 225 (720 x 450 is preferred).

Please remember to shrink the PNG file, if you can.
`This website <https://tinypng.com>`_ can help with that.


Before submitting a PR
----------------------

To save you some time and to make the new examples nicely fit into the
existing ones consider the following points.

**Look at examples source code** to get a feel for how statsmodels examples
should look like.

**Build the docs** by running `make html` from the docs directory to see how
your example looks in the fully rendered html pages.
.. _model:



Internal Classes
================

The following summarizes classes and functions that are not intended to be
directly used, but of interest only for internal use or for a developer who
wants to extend on existing model classes.


Module Reference
----------------

Model and Results Classes
^^^^^^^^^^^^^^^^^^^^^^^^^

These are the base classes for both the estimation models and the results.
They are not directly useful, but layout the structure of the subclasses and
define some common methods.

.. module:: statsmodels.base.model
   :synopsis: Base classes that are inherited by models

.. currentmodule:: statsmodels.base.model

.. autosummary::
   :toctree: generated/

   Model
   LikelihoodModel
   GenericLikelihoodModel
   Results
   LikelihoodModelResults
   ResultMixin
   GenericLikelihoodModelResults

.. module:: statsmodels.stats.contrast
   :synopsis: Classes for statistical test

.. currentmodule:: statsmodels.stats.contrast

.. autosummary::
   :toctree: generated/

   ContrastResults

.. inheritance-diagram:: statsmodels.base.model statsmodels.discrete.discrete_model statsmodels.regression.linear_model statsmodels.miscmodels.count
   :parts: 3

.. inheritance-diagram:: statsmodels.regression.linear_model.GLS statsmodels.regression.linear_model.WLS statsmodels.regression.linear_model.OLS statsmodels.regression.linear_model.GLSAR
   :parts: 1

Linear Model
^^^^^^^^^^^^

.. inheritance-diagram:: statsmodels.regression.linear_model
   :parts: 1

Generalized Linear Model
^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: statsmodels.genmod.generalized_linear_model
   statsmodels.genmod.families.family statsmodels.genmod.families.links
   :parts: 1

Discrete Model
^^^^^^^^^^^^^^

.. inheritance-diagram:: statsmodels.discrete.discrete_model
   :parts: 1

Robust Model
^^^^^^^^^^^^

.. inheritance-diagram:: statsmodels.robust.robust_linear_model
   :parts: 1

Vector Autoregressive Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: statsmodels.tsa.vector_ar.var_model
   :parts: 3
Naming Conventions
------------------

File and Directory Names
~~~~~~~~~~~~~~~~~~~~~~~~
Our directory tree stripped down looks something like::

    statsmodels/
        __init__.py
        api.py
        discrete/
            __init__.py
            discrete_model.py
            tests/
                results/
        tsa/
            __init__.py
            api.py
            tsatools.py
            stattools.py
            arima_process.py
            vector_ar/
                __init__.py
                var_model.py
                tests/
                    results/
            tests/
                results/
        stats/
            __init__.py
            api.py
            stattools.py
            tests/
        tools/
            __init__.py
            tools.py
            decorators.py
            tests/

The submodules are arranged by topic, `discrete` for discrete choice models, or `tsa` for time series
analysis. The submodules that can be import heavy contain an empty __init__.py, except for some testing
code for running tests for the submodules. The namespace to be imported is in `api.py`. That way, we
can import selectively and do not have to import a lot of code that we do not need. Helper functions are
usually put in files named `tools.py` and statistical functions, such as statistical tests are placed
in `stattools.py`. Everything has directories for :ref:`tests <testing>`.

`endog` & `exog`
~~~~~~~~~~~~~~~~

Our working definition of a statistical model is an object that has
both endogenous and exogenous data defined as well as a statistical
relationship.  In place of endogenous and exogenous one can often substitute
the terms left hand side (LHS) and right hand side (RHS), dependent and
independent variables, regressand and regressors, outcome and design, response
variable and explanatory variable, respectively.  The usage is quite often
domain specific; however, we have chosen to use `endog` and `exog` almost
exclusively, since the principal developers of statsmodels have a background
in econometrics, and this feels most natural.  This means that all of the
models are objects with `endog` and `exog` defined, though in some cases
`exog` is None for convenience (for instance, with an autoregressive process).
Each object also defines a `fit` (or similar) method that returns a
model-specific results object.  In addition there are some functions, e.g. for
statistical tests or convenience functions.

See also the related explanation in :ref:`endog_exog`.

Variable Names
~~~~~~~~~~~~~~
All of our models assume that data is arranged with variables in columns. Thus, internally the data
is all 2d arrays. By convention, we will prepend a `k_` to variable names that indicate moving over
axis 1 (columns), and `n_` to variables that indicate moving over axis 0 (rows). The main exception to
the underscore is that `nobs` should indicate the number of observations. For example, in the
time-series ARMA model we have::

    `k_ar` - The number of AR lags included in the RHS variables
    `k_ma` - The number of MA lags included in the RHS variables
    `k_trend` - The number of trend variables included in the RHS variables
    `k_exog` - The number of exogenous variables included in the RHS variables excluding the trend terms
    `n_totobs` - The total number of observations for the LHS variables including the pre-sample values


Options
~~~~~~~
We are using similar options in many classes, methods and functions. They
should follow a standardized pattern if they recur frequently. ::

    `missing` ['none', 'drop', 'raise'] define whether inputs are checked for
        nans, and how they are treated
    `alpha` (float in (0, 1)) significance level for hypothesis tests and
        confidence intervals, e.g. `alpha=0.05`

patterns ::

    `return_xxx` : boolean to indicate optional or different returns
        (not `ret_xxx`)
.. _git_notes:

Working with the statsmodels Code
=================================

Github
------

The `statsmodels` code base is hosted on `Github <https://github.com/statsmodels/statsmodels>`_. To
contribute you will need to `sign up for a free Github account <https://github.com/>`_.

Version Control and Git
-----------------------

We use the `Git <https://git-scm.com/>`_ version control system for development.
Git allows many people to work together on the same project.  In a nutshell, it
allows you to make changes to the code independent of others who may also be
working on the code and allows you to easily contribute your changes to the
codebase. It also keeps a complete history of all changes to the code, so you
can easily undo changes or see when a change was made, by whom, and why.

To install and configure Git, and to setup SSH keys, see
`setting up git <https://help.github.com/articles/set-up-git/>`_.

To learn more about Git, you may want to visit:

+ `Git documentation (book and videos) <https://git-scm.com/documentation>`_
+ `Github help pages <https://help.github.com/>`_
+ `NumPy documentation <https://docs.scipy.org/doc/numpy/dev/index.html>`_
+ `Matthew Brett's Pydagogue <https://matthew-brett.github.io/pydagogue/>`_

Below, we describe the bare minimum git commands you need to contribute to
`statsmodels`.

statsmodels Git/Github Workflow
-------------------------------

Forking and cloning
~~~~~~~~~~~~~~~~~~~

After setting up git, you need to fork the main `statsmodels` repository. To do
this, visit the `statsmodels project page
<https://github.com/statsmodels/statsmodels>`_ and hit the fork button (see
instructions for
`forking a repo <https://help.github.com/articles/fork-a-repo/>`_ for details).
This should take you to your fork's page.

Then, you want to clone the fork to your machine::

    git clone https://github.com/your-user-name/statsmodels
    cd statsmodels
    git remote add upstream https://github.com/statsmodels/statsmodels
    git fetch --all

The third line sets-up a read-only connection to the upstream statsmodels
repository. This will allow you to periodically update your local code with
changes in the upstream.  The final command fetches both your repository and
the upstream statsmodels repository.

Create a Branch
~~~~~~~~~~~~~~~

All changes to the code should be made in a feature branch. To create a branch, type::

    git checkout main
    git rebase upstream/main
    git checkout -b shiny-new-feature

The first two lines ensure you are starting from an up-to-date version of the upstream
statsmodels repository.  The third creates and checkout a new branch.

Doing::

    git branch

will give something like::

    * shiny-new-feature
      main

to indicate that you are now on the `shiny-new-feature` branch.

Making changes
~~~~~~~~~~~~~~

Hack away! Make any changes that you want, but please keep the work in your
branch completely confined to one specific topic, bugfix, or feature
implementation. You can work across multiple files and have many commits, but
the changes should all be related to the feature of the feature branch,
whatever that may be.

Now imagine that you changed the file `foo.py`. You can see your changes by
typing::

    git status

This will print something like::

    # On branch shiny-new-feature
    # Changes not staged for commit:
    #   (use "git add <file>..." to update what will be committed)
    #   (use "git checkout -- <file>..." to discard changes in working directory)
    #
    #       modified:   relative/path/to/foo.py
    #
    no changes added to commit (use "git add" and/or "git commit -a")

Before you can commit these changes, you have to `add`, or `stage`, the
changes. You can do this by typing::

    git add path/to/foo.py

Then check the status to make sure your commit looks okay::

    git status

should give something like::

    # On branch shiny-new-feature
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #       modified:   /relative/path/to/foo.py
    #

Pushing your changes
~~~~~~~~~~~~~~~~~~~~

At any time you can push your feature branch (and any changes) to your github
(fork) repository by::

    git push

although the first time you will need to run

    git push --set-upstream origin shiny-new-feature

to instruct git to set the current branch to track its corresponding branch in
your github repository.

You can see the remote repositories by::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin  https://github.com/yourname/statsmodels.git (fetch)
    origin  https://github.com/yourname/statsmodels.git (push)
    upstream        https://github.com/statsmodels/statsmodels.git (fetch)
    upstream        https://github.com/statsmodels/statsmodels.git (push)

Before you push any commits, however, it is *highly* recommended that you make
sure what you are pushing makes sense and looks clean. You can review your
change history by::

    git log --oneline --graph

It pays to take care of things locally before you push them to github. So when
in doubt, do not push.  Also see the advice on keeping your history clean in
:ref:`merge-vs-rebase`.

.. _pull-requests:

Pull Requests
~~~~~~~~~~~~~

When you are ready to ask for a code review, we recommend that you file a pull
request. Before you do so you should check your changeset yourself. You can do
this by using `compare view
<https://github.com/blog/612-introducing-github-compare-view>`__ on github.

#. Navigate to your repository on github.
#. Click on `Branch List`
#. Click on the `Compare` button for your feature branch, `shiny-new-feature`.
#. Select the `base` and `compare` branches, if necessary. This will be `main` and
   `shiny-new-feature`, respectively.
#. From here you will see a nice overview of your changes. If anything is amiss, you can fix it.

If everything looks good you are read to make a `pull request <https://help.github.com/articles/about-pull-requests/>`__.

#. Navigate to your repository on github.
#. Click on the `Pull Request` button.
#. You can then click on `Commits` and `Files Changed` to make sure everything looks okay one last time.
#. Write a description of your changes in the `Preview Discussion` tab.
#. Click `Send Pull Request`.

Your request will then be reviewed. If you need to go back and make more
changes, you can make them in your branch and push them to github and the pull
request will be automatically updated.

One last thing to note. If there has been a lot of work in upstream/main
since you started your patch, you might want to rebase. However, you can
probably get away with not rebasing if these changes are unrelated to the work
you have done in the `shiny-new-feature` branch. If you can avoid it, then
do not rebase. If you have to, try to do it once and when you are at the end of
your changes. Read on for some notes on :ref:`merge-vs-rebase`.

Advanced Topics
---------------

.. _merge-vs-rebase:

Merging vs. Rebasing
~~~~~~~~~~~~~~~~~~~~

This is a topic that has been discussed at great length and with considerable
more expertise than we can offer here. This section will provide some resources
for further reading and some advice. The focus, though, will be for those who
wish to submit pull requests for a feature branch. For these cases rebase
should be preferred.

A rebase replays commits from one branch on top of another branch to preserve a
linear history. Recall that your commits were tested against a (possibly) older
version of main from which you started your branch, so if you rebase, you
could introduce bugs. However, if you have only a few commits, this might not
be such a concern. One great place to start learning about rebase is
:ref:`rebasing without tears <pydagogue:actual-rebase>`.  In particular, `heed
the warnings
<https://matthew-brett.github.io/pydagogue/rebase_without_tears.html#safety>`__.
Namely, **always make a new branch before doing a rebase**. This is good
general advice for working with git. I would also add **never use rebase on
work that has already been published**. If another developer is using your
work, do not rebase!!

As for merging, **never merge from trunk into your feature branch**. You will,
however, want to check that your work will merge cleanly into trunk. This will
help out the reviewers. You can do this in your local repository by merging
your work into your main branch (or any branch that tracks the remote main branch) and
:ref:`run-tests`.

Deleting Branches
~~~~~~~~~~~~~~~~~

Once your feature branch is accepted into upstream, you might want to get rid
of it. First you'll want to merge upstream main into your branch. That way
git will know that it can safely delete your branch::

    git fetch upstream
    git checkout main
    git merge upstream/main

Then you can just do::

    git branch -d shiny-new-feature

Make sure you use a lower-case -d. That way, git will complain if your feature
branch has not actually been merged. The branch will still exist on github
however. To delete the branch on github, do::

    git push origin :shiny-new-feature branch

.. Squashing with Rebase
.. ^^^^^^^^^^^^^^^^^^^^^

.. You have made a bunch of incremental commits, but you think they might be better off together as one
.. commit. You can do this with an interactive rebase. As usual, **only do this when you have local
.. commits. Do not edit the history of changes that have been pushed.**

.. see this reference http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html
.. _add_data:

Datasets
========

For a list of currently available datasets and usage instructions, see the
:ref:`datasets page <datasets>`.

License
-------

To be considered for inclusion in `statsmodels`, a dataset must be in the
public domain, distributed under a BSD-compatible license, or we must obtain
permission from the original author.

Adding a dataset: An example
----------------------------

The Nile River data measures the volume of the discharge of the Nile River at
Aswan for the years 1871 to 1970. The data are copied from the paper of Cobb
(1978).

**Step 1**: Create a directory `datasets/nile/`

**Step 2**: Add `datasets/nile/nile.csv` and  a new file `datasets/__init__.py` which contains ::

    from data import *

**Step 3**: If `nile.csv` is a transformed/cleaned version of the original data, create a `nile/src` directory and include the original raw data there. In the `nile` case, this step is not necessary.

**Step 4**: Copy `datasets/template_data.py` to `nile/data.py`. Edit `nile/data.py` by filling-in strings for COPYRIGHT, TITLE, SOURCE, DESCRSHORT, DESCLONG, and NOTE. ::

    COPYRIGHT   = """This is public domain."""
    TITLE       = """Nile River Data"""
    SOURCE      = """
    Cobb, G.W. 1978. The Problem of the Nile: Conditional Solution to a Changepoint
        Problem. Biometrika. 65.2, 243-251,
    """

    DESCRSHORT  = """Annual Nile River Volume at Aswan, 1871-1970""

    DESCRLONG   = """Annual Nile River Volume at Aswan, 1871-1970. The units of
    measurement are 1e9 m^{3}, and there is an apparent changepoint near 1898."""

    NOTE        = """
    Number of observations: 100
    Number of variables: 2
    Variable name definitions:
        year - Year of observation
        volume - Nile River volume at Aswan

    The data were originally used in Cobb (1987, See SOURCE). The author
    acknowledges that the data were originally compiled from various sources by
    Dr. Barbara Bell, Center for Astrophysics, Cambridge, Massachusetts. The data
    set is also used as an example in many textbooks and software packages.
    """

**Step 5:** Edit the docstring of the `load` function in `data.py` to specify
which dataset will be loaded. Also edit the path and the indices for the
`endog` and `exog` attributes. In the `nile` case, there is no `exog`, so
everything referencing `exog` is not used. The `year` variable is also not
used.

**Step 6:** Edit the `datasets/__init__.py` to import the directory.

That's it! The result can be found `here
<https://github.com/statsmodels/statsmodels/tree/main/statsmodels/datasets/nile>`_
for reference.
Developer Page
==============

This page explains how you can contribute to the development of `statsmodels`
by submitting patches, statistical tests, new models, or examples.

`statsmodels` is developed on `Github
<https://github.com/statsmodels/statsmodels>`_ using the `Git
<https://git-scm.com/>`_ version control system.

Submitting a Bug Report
-----------------------

- Include a short, self-contained code snippet that reproduces the problem
- Specify the statsmodels version used. You can do this with ``sm.version.full_version``
- If the issue looks to involve other dependencies, also include the output of ``sm.show_versions()``

Making Changes to the Code
--------------------------

First, review the :ref:`git_notes` section for an intro to the git version
control system.

For a pull request to be accepted, you must meet the below requirements. This
greatly helps the job of maintaining and releasing the software a shared effort.

- **One branch. One feature.** Branches are cheap and github makes it easy to
  merge and delete branches with a few clicks. Avoid the temptation to lump in a
  bunch of unrelated changes when working on a feature, if possible. This helps
  us keep track of what has changed when preparing a release.
- Commit messages should be clear and concise. This means a subject line of
  less than 80 characters, and, if necessary, a blank line followed by a commit
  message body. We have an
  `informal commit format standard <https://www.statsmodels.org/devel/dev/maintainer_notes.html#commit-comments>`_
  that we try to adhere to. You can see what this looks like in practice by
  ``git log --oneline -n 10``. If your commit references or closes a specific
  issue, you can close it by mentioning it in the
  `commit message <https://help.github.com/articles/closing-issues-via-commit-messages/>`_.
  (*For maintainers*: These suggestions go for Merge commit comments too. These are partially the record for release notes.)
- Code submissions must always include tests. See our notes on :ref:`testing`.
- Each function, class, method, and attribute needs to be documented using
  docstrings. We conform to the
  `numpy docstring standard <https://numpy.org/doc/stable/docs/howto_document.html#docstring-standard>`_.
- If you are adding new functionality, you need to add it to the documentation
  by editing (or creating) the appropriate file in ``docs/source``.
- Make sure your documentation changes parse correctly. Change into the
  top-level ``docs/`` directory and type::

    make clean
    make html

  Check that the build output does not have *any* warnings due to your changes.
- Follow `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style guidelines
  wherever possible. Compare your code to what's in main by running
  ``git diff upstream/main -u -- "*.py" | flake8 --diff --isolated`` prior to submitting.
- Finally, please add your changes to the release notes. Open the
  ``docs/source/release/versionX.X.rst`` file that has the version number of the
  next release and add your changes to the appropriate section.

How to Submit a Pull Request
----------------------------

So you want to submit a patch to `statsmodels` but are not too familiar with
github? Here are the steps you need to take.

1. `Fork <https://help.github.com/articles/fork-a-repo/>`_ the
   `statsmodels repository <https://github.com/statsmodels/statsmodels>`_ on Github.
2. `Create a new feature branch <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_.
   Each branch must be self-contained, with a single new feature or bugfix.
3. Make sure the test suite passes. This includes testing on Python 3. The
   easiest way to do this is to make a pull request and let the bot check for
   you. This can be slow, and if you are unsure about the fix or enhancement,
   it is best to run pytest locally.
4. `Submit a pull request <https://help.github.com/articles/about-pull-requests/>`_

Pull requests are thoroughly reviewed before being accepted into the codebase.
If your pull request becomes out of date, rebase your pull request on the
latest version in the central repository.


Mailing List
------------

Conversations about development take place on the
`statsmodels mailing list <https://groups.google.com/forum/?hl=en#!forum/pystatsmodels>`__.

License
-------

statsmodels is released under
the `Modified (3-clause) BSD license <https://opensource.org/licenses/BSD-3-Clause>`_.

Contents
--------

.. toctree::
   :maxdepth: 1

   git_notes
   maintainer_notes
   test_notes
   naming_conventions
   warnings-and-exceptions
   dataset_notes
   examples
   get_involved
   internal
   testing
.. _testing:

Testing
=======

Setting up development environment locally
------------------------------------------
Follow our :ref:`installation instructions <install>` and set up a suitable
environment to build statsmodels from source. We recommend that you develop
using a development install of statsmodels by running:

.. code-block:: bash

   pip install -e .

from the root directory of the git repository. The flag ``-e`` is for editable.

This command compiles the C code and add statsmodels to your activate python
environment by creating links from your python environment's libraries
to the statsmodels source code. Therefore, changes to pure python code will
be immediately available to the user without a re-install. Changes to C code or
Cython code require rerunning ``pip install -e .`` before these changes are
available.

Test Driven Development
-----------------------
We strive to follow a `Test Driven Development (TDD) <https://en.wikipedia.org/wiki/Test-driven_development>`_ pattern.
All models or statistical functions that are added to the main code base are to have
tests versus an existing statistical package, if possible.

Introduction to pytest
----------------------
Like many packages, statsmodels uses the `pytest testing system <https://docs.pytest.org/en/latest/contents.html>`__ and the convenient extensions in `numpy.testing <https://docs.scipy.org/doc/numpy/reference/routines.testing.html>`__.  Pytest will find any file, directory, function, or class name that starts with ``test`` or ``Test`` (classes only). Test function should start with ``test``, test classes should start with ``Test``. These functions and classes should be placed in files with names beginning with ``test`` in a directory called ``tests``.

.. _run-tests:

Running Tests
-------------
Test are run from the command line by calling ``pytest``. Directly running tests using
pytest requires that statsmodels is installed using ``pip install -e .`` as described
above.

Tests can be run at different levels of granularity:

* Project level, which runs all tests.  Running the entire test suite is slow
  and normally this would only be needed if making deep changes to statsmodels.

.. code-block:: bash

    pytest statsmodels

* Folder level, which runs all tests below a folder

.. code-block:: bash

    pytest statsmodels/regression/tests

* File level, which runs all tests in a file

.. code-block:: bash

    pytest statsmodels/regression/tests/test_regression.py

* Class level, which runs all tests in a class

.. code-block:: bash

    pytest statsmodels/regression/tests/test_regression.py::TestOLS

* Test level, which runs a single test.  The first example runs a test in a
  class.  The second runs a stand alone test.

.. code-block:: bash

    pytest statsmodels/regression/tests/test_regression.py::TestOLS::test_missing
    pytest statsmodels/regression/tests/test_regression.py::test_ridge

How To Write A Test
-------------------
NumPy provides a good introduction to unit testing with pytest and NumPy extensions `here <https://github.com/numpy/numpy/blob/main/doc/TESTS.rst.txt>`__. It is worth a read for some more details.
Here, we will document a few conventions we follow that are worth mentioning. Often we want to test
a whole model at once rather than just one function, for example. The following is a pared down
version test_discrete.py. In this case, several different models with different options need to be
tested. The tests look something like

.. code-block:: python

    from numpy.testing import assert_almost_equal
    import statsmodels.api as sm
    from results.results_discrete import Spector

    class CheckDiscreteResults(object):
        """
        res2 are the results. res1 are the values from statsmodels
        """

        def test_params(self):
            assert_almost_equal(self.res1.params, self.res2.params, 4)

        decimal_tvalues = 4
        def test_tvalues(self):
            assert_almost_equal(self.res1.params, self.res2.params, self.decimal_tvalues)

        # ... as many more tests as there are common results

    class TestProbitNewton(CheckDiscreteResults):
        """
        Tests the Probit model using Newton's method for fitting.
        """

        @classmethod
        def setup_class(cls):
            # set up model
            data = sm.datasets.spector.load()
            data.exog = sm.add_constant(data.exog)
            cls.res1 = sm.Probit(data.endog, data.exog).fit(method='newton', disp=0)

            # set up results
            res2 = Spector.probit
            cls.res2 = res2

            # set up precision
            cls.decimal_tvalues = 3

        def test_model_specifc(self):
            assert_almost_equal(self.res1.foo, self.res2.foo, 4)

The main workhorse is the `CheckDiscreteResults` class. Notice that we can set the level of precision
for `tvalues` to be different than the default in the subclass  `TestProbitNewton`. All of the test
classes have a ``@classmethod`` called ``setup_class``. Otherwise, pytest would reinstantiate the class
before every single test method. If the fitting of the model is time consuming, then this is clearly
undesirable. Finally, we have a script at the bottom so that we can run the tests should be running
the Python file.

Test Results
------------
The test results are the final piece of the above example. For many tests, especially those for the
models, there are many results against which you would like to test. It makes sense then to separate
the hard-coded results from the actual tests to make the tests more readable. If there are only a few
results it's not necessary to separate the results. We often take results from some other statistical
package. It is important to document where you got the results from and why they might differ from
the results that we get. Each tests folder has a results subdirectory. Consider the folder structure
for the discrete models::

    tests/
        __init__.py
        test_discrete.py
        results/
            __init__.py
            results_discrete.py
            nbinom_resids.csv

It is up to you how best to structure the results. In the discrete model example, you will notice
that there are result classes based around particular datasets with a method for loading different
model results for that dataset. You can also include text files that hold results to be loaded by
results classes if it is easier than putting them in the class itself.

Speeding up full runs
---------------------
Running the full test suite is slow. Fortunately it is only necessary to run the full suite when
making low-level changes (e.g., to ``statsmodels.base``) There are two methods available to
speed up runs of the full test suite when needed.

* Use the pytest-xdist package

.. code-block:: bash

   pip install pytest-xdist
   export MKL_NUM_THREADS=1
   export OMP_NUM_THREADS=1
   pytest -n auto statsmodels

* Skip slow tests using ``--skip-slow``

.. code-block:: bash

   pytest --skip-slow statsmodels


You can combine these two approaches for faster runs.

.. code-block:: bash

   export MKL_NUM_THREADS=1 && export OMP_NUM_THREADS=1
   pytest -n auto --skip-slow statsmodels


The ``test()`` method
---------------------
The root of statsmodels and all submodules expose a ``test()`` method which can
be used to run all tests either in the package (``statsmodels.test()``) or in
a module (``statsmodels.regression.test()``).  This method allows tests to be
run from an install copy of statsmodels even it is was not installed using the
*editable* flag as described above. This method is required for testing wheels in
release builds and is **not** recommended for development.

Using this method, all tests are run using:

.. code-block:: python

   import statsmodels.api as sm
   sm.test()

Submodules tests are run using:

.. code-block:: python

    sm.discrete.test()

.. autosummary::
   :toctree: generated/

   ~statsmodels.__init__.test
Maintainer Notes
================

This is for those with read-write access to upstream. It is recommended to name
the upstream remote something to remind you that it is read-write::

    git remote add upstream-rw git@github.com:statsmodels/statsmodels.git
    git fetch upstream-rw

Git Workflow
------------

Grabbing Changes from Others
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you need to push changes from others, you can link to their repository by doing::

    git remote add contrib-name git://github.com/contrib-name/statsmodels.git
    get fetch contrib-name
    git branch shiny-new-feature --track contrib-name/shiny-new-feature
    git checkout shiny-new-feature

The rest of the below assumes you are on your or someone else's branch with the changes you
want to push upstream.

.. _rebasing:

Rebasing
^^^^^^^^

If there are only a few commits, you can rebase to keep a linear history::

    git fetch upstream-rw
    git rebase upstream-rw/main

Rebasing will not automatically close the pull request however, if there is one,
so do not forget to do this.

.. _merging:

Merging
^^^^^^^

If there is a long series of related commits, then you'll want to merge. You
may ask yourself, :ref:`ff-no-ff`? See below for more on this choice. Once
decided you can do::

    git fetch upstream-rw
    git merge --no-ff upstream-rw/main

Merging will automatically close the pull request on github.

Check the History
^^^^^^^^^^^^^^^^^

This is very important. Again, any and all fixes should be made locally before
pushing to the repository::

    git log --oneline --graph

This shows the history in a compact way of the current branch. This::

    git log -p upstream-rw/main..

shows the log of commits excluding those that can be reached from
upstream-rw/main, and including those that can be reached from current HEAD.
That is, those changes unique to this branch versus upstream-rw/main. See
:ref:`Pydagogue <pydagogue:git-log-dots>` for more on using dots with log and
also for using :ref:`dots with diff <pydagogue:git-diff-dots>`.

Push Your Feature Branch
^^^^^^^^^^^^^^^^^^^^^^^^

All the changes look good? You can push your feature branch after
:ref:`merging` or :ref:`rebasing` by::

    git push upstream-rw shiny-new-feature:main

Cherry-Picking
^^^^^^^^^^^^^^

Say you are interested in some commit in another branch, but want to leave the
other ones for now. You can do this with a cherry-pick. Use `git log --oneline`
to find the commit that you want to cherry-pick. Say you want commit `dd9ff35`
from the `shiny-new-feature` branch. You want to apply this commit to main.
You simply do::

    git checkout main
    git cherry-pick dd9ff35

And that's all. This commit is now applied as a new commit in main.

.. _ff-no-ff:

Merging: To Fast-Forward or Not To Fast-Forward
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, `git merge` is a fast-forward merge. What does this mean, and when
do you want to avoid this?

.. figure:: images/git_merge.png
   :alt: git merge diagram
   :scale: 100%
   :align: center

   (source `nvie.com <http://nvie.com>`__, post `"A successful Git branching model" <http://nvie.com/posts/a-successful-git-branching-model/>`__)

The fast-forward merge does not create a merge commit. This means that the
existence of the feature branch is lost in the history. The fast-forward is the
default for Git basically because branches are cheap and, therefore, *usually*
short-lived. If on the other hand, you have a long-lived feature branch or are
following an iterative workflow on the feature branch (i.e. merge into main,
then go back to feature branch and add more commits), then it makes sense to
include only the merge in the main branch, rather than all the intermediate
commits of the feature branch, so you should use::

    git merge --no-ff

Handling Pull Requests
^^^^^^^^^^^^^^^^^^^^^^

You can apply a pull request through `fetch <https://www.kernel.org/pub/software/scm/git/docs/git-fetch.html>`__
and `merge <https://www.kernel.org/pub/software/scm/git/docs/git-merge.html>`__.
In your local copy of the main repo::

    git checkout main
    git remote add contrib-name git://github.com/contrib-name/statsmodels.git
    git fetch contrib-name
    git merge contrib-name/shiny-new-feature

Check that the merge applies cleanly and the history looks good. Edit the merge
message. Add a short explanation of what the branch did along with a
'Closes gh-XXX.' string. This will auto-close the pull request and link the
ticket and closing commit. To automatically close the issue, you can use any
of::

    gh-XXX
    GH-XXX
    #XXX

in the commit message. Any and all problems need to be taken care of locally
before doing::

    git push origin main

Releasing
---------

1. Checkout main::

    git checkout statsmodels/main

2. Clean the working tree with::

    git clean -xdf

   But you might want to do a dry-run first::

    git clean -xdfn

3. **Locally** tag the release. For a release candidate, for example::

    git tag -a v0.10.0rc1 -m "Version 0.10.0 Release Candidate 1" 7b2fb29

   or just::

    git tag -a v0.10.0rc1 -m "Version 0.10.0 Release Candidate 1"

   to use the last commit in main.

4. Checkout the tag::

    git checkout tags/v0.10.0rc1

5. Build a sdist to ensure that that the build is clean::

    python setup.py sdist --formats=gztar

   It is important that the build on the tar.gz file is the same as the tag. It must not be **dirty**

6. If on a new minor release (major.minor.micro format) start a new maintenance branch, for example::

    git checkout -b maintenance/0.10.x

   Any bug fixes and maintenance commits intended for the next micro release should be made
   against main as usual, but tagged with the milestone for the micro release it is intended
   for. Then merge into main as usual. When ready to do the backports, use the file
   ``tools/backport_pr.py`` to identify which PRs need to be backported and to apply them to the
   maintenance branch. The tag for the release should be made in the maintenance branch.

7. Upload the source distribution to PyPI::

    twine upload dist/*

   You might want to upload to test first::

    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

8. Go back to the main branch, and add an empty commit::

    git checkout statsmodels/main
    git commit --allow-empty -m "Start of 0.11.0 development"
    git tag -a v0.11.0.dev0 -m "Start of 0.11.0 development"

9. Push everything to statsmodels::

    git push --tags

   If a new branch was created::

    git push --set-upstream origin maintenance/0.10.x

10. Make an announcement, and inform maintainers of wheel builders.

11. Profit?

Releasing from Maintenance Branch
---------------------------------

Once any patches have been backported to a maintenance branch, the release steps are

1. Checkout the branch::

    git checkout maintenance/0.10.x

2. Clean up thoroughly::

    git clean -xdf


3. **Locally** tag the release::

    git tag -a v0.10.0 -m "Version 0.10.0"

4. Checkout the tag::

    git checkout tags/v0.10.0

5. Build a sdist to ensure that that the build is clean::

    python setup.py sdist --formats=gztar

   It is important that the build on the tar.gz file is the same as the tag. It must not be **dirty**.

6. Upload the source distribution to PyPI ot PyPI test::

    twine upload dist/*

   or::

    twine upload --repository-url https://test.pypi.org/legacy/ dist/*


7. Push the tag to statsmodels::

    git push --tags


8. Make an announcement, and inform maintainers of wheel builders.


Commit Comments
---------------
Prefix commit messages in the main branch of the main shared repository with
the following::

    ENH: Feature implementation
    BUG: Bug fix
    STY: Coding style changes (indenting, braces, code cleanup)
    DOC: Sphinx documentation, docstring, or comment changes
    CMP: Compiled code issues, regenerating C code with Cython, etc.
    REL: Release related commit
    TST: Change to a test, adding a test. Only used if not directly related to a bug.
    REF: Refactoring changes
Testing on Build Machines
-------------------------

There are currently several places that statsmodels is automatically built and
tested against different dependency and Python versions and architectures.
Check these logs periodically, make sure everything looks okay, and fix any
failures:

* `Azure Pipelines <https://dev.azure.com/statsmodels/statsmodels-testing/_build?definitionId=1&_a=summary>`_


The test coverage pages are here:

* `Coveralls <https://coveralls.io/github/statsmodels/statsmodels>`_
* `codecov <https://codecov.io/gh/statsmodels/statsmodels>`_

We use `Codacy <https://app.codacy.com/project/josef-pkt/statsmodels/dashboard>`_
to monitor code quality.
.. _datasets:

.. currentmodule:: statsmodels.datasets

.. ipython:: python
   :suppress:

   import numpy as np
   np.set_printoptions(suppress=True)

The Datasets Package
====================

``statsmodels`` provides data sets (i.e. data *and* meta-data) for use in
examples, tutorials, model testing, etc.

Using Datasets from Stata
-------------------------

.. autosummary::
   :toctree: ./

   webuse

Using Datasets from R
---------------------

The `Rdatasets project <https://vincentarelbundock.github.io/Rdatasets/>`__ gives access to the datasets available in R's core datasets package and many other common R packages. All of these datasets are available to statsmodels by using the :func:`get_rdataset` function. The actual data is accessible by the ``data`` attribute. For example:

.. ipython:: python

   import statsmodels.api as sm
   duncan_prestige = sm.datasets.get_rdataset("Duncan", "carData")
   print(duncan_prestige.__doc__)
   duncan_prestige.data.head(5)


R Datasets Function Reference
-----------------------------


.. autosummary::
   :toctree: ./

   get_rdataset
   get_data_home
   clear_data_home


Available Datasets
------------------

.. toctree::
   :maxdepth: 1
   :glob:

   generated/*

Usage
-----

Load a dataset:

.. ipython:: python

   import statsmodels.api as sm
   data = sm.datasets.longley.load_pandas()

The `Dataset` object follows the bunch pattern. The full dataset is available
in the ``data`` attribute.

.. ipython:: python

   data.data

Most datasets hold convenient representations of the data in the attributes `endog` and `exog`:

.. ipython:: python

   data.endog.iloc[:5]
   data.exog.iloc[:5,:]

Univariate datasets, however, do not have an `exog` attribute.

Variable names can be obtained by typing:

.. ipython:: python

   data.endog_name
   data.exog_name

If the dataset does not have a clear interpretation of what should be an
`endog` and `exog`, then you can always access the `data` or `raw_data`
attributes. This is the case for the `macrodata` dataset, which is a collection
of US macroeconomic data rather than a dataset with a specific example in mind.
The `data` attribute contains a record array of the full dataset and the
`raw_data` attribute contains an ndarray with the names of the columns given
by the `names` attribute.

.. ipython:: python

   type(data.data)
   type(data.raw_data)
   data.names

Loading data as pandas objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many users it may be preferable to get the datasets as a pandas DataFrame or
Series object. Each of the dataset modules is equipped with a ``load_pandas``
method which returns a ``Dataset`` instance with the data readily available as pandas objects:

.. ipython:: python

   data = sm.datasets.longley.load_pandas()
   data.exog
   data.endog

The full DataFrame is available in the ``data`` attribute of the Dataset object

.. ipython:: python

   data.data


With pandas integration in the estimation classes, the metadata will be attached
to model results:

.. ipython:: python
   :okwarning:

   y, x = data.endog, data.exog
   res = sm.OLS(y, x).fit()
   res.params
   res.summary()

Extra Information
^^^^^^^^^^^^^^^^^

If you want to know more about the dataset itself, you can access the
following, again using the Longley dataset as an example ::

    >>> dir(sm.datasets.longley)[:6]
    ['COPYRIGHT', 'DESCRLONG', 'DESCRSHORT', 'NOTE', 'SOURCE', 'TITLE']

Additional information
----------------------

* The idea for a datasets package was originally proposed by David Cournapeau.
* To add datasets, see the :ref:`notes on adding a dataset <add_data>`.
:orphan:

.. _statsmodels-examples:

Examples
========

This page provides a series of examples, tutorials and recipes to help you get
started with ``statsmodels``. Each of the examples shown here is made available
as an IPython Notebook and as a plain python script on the `statsmodels github
repository <https://github.com/statsmodels/statsmodels/tree/main/examples>`_.

We also encourage users to submit their own examples, tutorials or cool
`statsmodels` trick to the `Examples wiki page
<https://github.com/statsmodels/statsmodels/wiki/Examples>`_

.. toctree::


{# This content is white space sensitive. Do not reformat #}

{% for category in examples%}
{% set underscore = "-" * (category.header | length) %}

.. raw:: html

   <div class="example-header-clear">&nbsp;</div>

{{ category.header }}
{{ underscore }}

.. toctree::
   :maxdepth: 1
   :hidden:

{% for notebook in category.links  %}   {{ notebook.target | replace('.html','') }}
{% endfor %}

{%- for notebook in category.links  %}
{% set heading = "`" ~ notebook.text ~ " <" ~ notebook.target|e ~ ">`_" %}

.. container:: example

   {{ heading }}

   .. image:: {{ notebook.img }}
      :target: {{ notebook.target }}
      :width: 240px

{%- endfor %}

{% endfor %}
Copula tests
############

The reference results are coming from the R package Copula. The following
script is used:

    library(copula)
    sessionInfo()

    u <- matrix(c(0.33706249, 0.62232507,
                  0.2001457 , 0.77166391,
                  0.98534253, 0.72755898,
                  0.05943888, 0.0962475 ,
                  0.35496733, 0.44513594,
                  0.6075078 , 0.06241089,
                  0.54027684, 0.40610225,
                  0.99212789, 0.25913165,
                  0.61044613, 0.67585563,
                  0.79584436, 0.23050014),
                nrow=10)

    print("Vector u")
    print(u)

    gaussian <- normalCopula(0.8, dim = 2)
    student <- tCopula(0.8, dim = 2, df = 2)
    frank <- frankCopula(dim = 2, param = 3)
    clayton <- claytonCopula(dim = 2, param = 1.2)
    gumbel <- gumbelCopula(dim = 2, param = 1.5)

    # Compute the density and CDF
    pdf <- dCopula(u, gaussian)
    cdf <- pCopula(u, gaussian)
    print("Gaussian")
    print(pdf)
    print(cdf)

    pdf <- dCopula(u, student)
    cdf <- pCopula(u, student)
    print("Student")
    print(pdf)
    print(cdf)

    pdf <- dCopula(u, frank)
    cdf <- pCopula(u, frank)
    print("Frank")
    print(pdf)
    print(cdf)

    pdf <- dCopula(u, clayton)
    cdf <- pCopula(u, clayton)
    print("Clayton")
    print(pdf)
    print(cdf)

    pdf <- dCopula(u, gumbel)
    cdf <- pCopula(u, gumbel)
    print("Gumbel")
    print(pdf)
    print(cdf)

Which produces the following output:

    R version 4.0.3 (2020-10-10)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 20.04.1 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] copula_1.0-0

    loaded via a namespace (and not attached):
     [1] compiler_4.0.3      Matrix_1.2-18       ADGofTest_0.3
     [4] pspline_1.0-18      gsl_2.1-6           mvtnorm_1.1-1
     [7] grid_4.0.3          pcaPP_1.9-73        numDeriv_2016.8-1.1
    [10] stats4_4.0.3        lattice_0.20-41     stabledist_0.7-1
    [1] "Vector u"
                [,1]       [,2]
     [1,] 0.33706249 0.60750780
     [2,] 0.62232507 0.06241089
     [3,] 0.20014570 0.54027684
     [4,] 0.77166391 0.40610225
     [5,] 0.98534253 0.99212789
     [6,] 0.72755898 0.25913165
     [7,] 0.05943888 0.61044613
     [8,] 0.09624750 0.67585563
     [9,] 0.35496733 0.79584436
    [10,] 0.44513594 0.23050014
    [1] "Gaussian"
     [1]  1.03308741  0.06507279  0.72896012  0.65389439 16.45012399 0.34813218
     [7]  0.06768115  0.08168840  0.40521741  1.26723470
     [1] 0.31906854 0.06230196 0.19284669 0.39952707 0.98144792 0.25677003
     [7] 0.05932818 0.09605404 0.35211017 0.20885480
    [1] "Student"
     [1]  0.8303065  0.1359839  0.5157746  0.4776421 26.2173959  0.3070661
     [7]  0.1349173  0.1597064  0.3303230  1.0482301
     [1] 0.31140349 0.05942746 0.18548601 0.39143974 0.98347259 0.24894028
     [7] 0.05653947 0.09210693 0.34447385 0.20429882
    [1] "Frank"
     [1] 0.9646599 0.5627195 0.8941964 0.8364614 2.9570945 0.6665601 0.5779906
     [8] 0.5241333 0.7156741 1.1074024
     [1] 0.27467496 0.05492539 0.15995939 0.36750702 0.97782283 0.23412757
     [7] 0.05196265 0.08676979 0.32803721 0.16320730
    [1] "Clayton"
     [1] 1.0119836 0.2072728 0.8148839 0.9481976 2.1419659 0.6828507 0.2040454
     [8] 0.2838497 0.8197787 1.1096360
     [1] 0.28520375 0.06101690 0.17703377 0.36848218 0.97772088 0.24082057
     [7] 0.05811908 0.09343934 0.33012582 0.18738753
    [1] "Gumbel"
     [1]  1.0391696  0.6539579  0.9878446  0.8679504 16.6030932  0.7542073
     [7]  0.6668307  0.6275887  0.7477991  1.1564864
     [1] 0.27194634 0.05484380 0.15668190 0.37098420 0.98176346 0.23422865
     [7] 0.05188260 0.08659615 0.33086960 0.15803914
