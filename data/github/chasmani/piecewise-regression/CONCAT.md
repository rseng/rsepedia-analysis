---
title: 'piecewise-regression (aka segmented regression) in Python'
tags:
  - Python
  - regression
  - statistics
  - segmented regression
  - breakpoint analysis
authors:
  - name: Charlie Pilgrim
    orcid: 0000-0002-3800-677X
    affiliation: "1, 2" 
affiliations:
 - name: Centre for Doctoral Training in Mathematics for Real-World Systems, University of Warwick, Coventry, UK
   index: 1
 - name: The Alan Turing Institute, London, UK
   index: 2
date: 18 November 2021
bibliography: paper.bib

---

# Summary

Piecewise regression (also known as segmented regression, broken-line regression, or breakpoint analysis) fits a linear regression model to data that includes one or more breakpoints where the gradient changes. The `piecewise-regression` Python package uses the approach described by Muggeo [@muggeo2003estimating], where the breakpoint positions and the straight line models are simultaneously fit using an iterative method. This easy-to-use package includes an automatic comprehensive statistical analysis that gives confidence intervals for all model variables and hypothesis testing for the existence of breakpoints. 

# Statement of Need

A common problem in many fields is to fit a continuous straight line model to data that includes some change(s) in gradient known as breakpoint(s). Examples include investigating medical interventions [@wagner2002segmented], ecological thresholds [@toms2003piecewise], and geological phase transitions [@ryan2002defining]. Fitting such models involves the global problem of finding estimates for the breakpoint positions and the local problem of fitting line segments given breakpoints. Possible approaches involve using linear regression to fit line segments together with a global optimisation algorithm to find breakpoints—for example, an evolutionary algorithm as in the `pwlf` python package [@jekel2019pwlf]. Or one could take a non-linear least-squares approach using `scipy` [@virtanen2020scipy] or the `lmfit` python package [@newville2016lmfit]. Muggeo [@muggeo2003estimating] derived an alternative method whereby the breakpoint positions and the line segment models are fitted simultaneously using an iterative method, which is computationally efficient and allows for robust statistical analysis. Many R packages implement this method, including the `segmented` R package written by Muggeo himself [@muggeo2008segmented]. However, before the `piecewise-regression` package, there were not comparable resources in Python.

# Example

An example plot is shown in \autoref{fig:example}. Data was generated with 3 breakpoints and some noise, and a model was then fit to that data. The plot shows the maximum likelihood estimators for the straight line segments and breakpoint positions. The package automatically carries out a Davies hypothesis test [@davies1987hypothesis] for the existence of at least 1 breakpoint, in this example finding strong evidence for breakpoints with $p<0.001$.

![An example model fit (red line) to data (grey markers). The estimated breakpoint positions (blue lines) and confidence intervals (shaded blue regions) are shown. The data was generated using a piecewise linear model with a constant level of Gaussian noise. For example, this could represent observations with a sampling error of some physical process that undergoes phase transitions. \label{fig:example}](example.png)

# How It Works

We follow here the derivation by Muggeo [@muggeo2003estimating]. The general form of the model with one breakpoint is

\begin{equation}
    y = \alpha x + c + \beta (x-\psi) H(x-\psi) + \zeta \,,
\end{equation}

where given some data, $x$, $y$, we are trying to estimate the gradient of the first segment, $\alpha$, the intercept of the first segment, $c$, the change in gradient from the first to second segments, $\beta$, and the breakpoint position, $\psi$. $H$ is the Heaviside step function and $\zeta$ is a noise term. This cannot be solved directly through linear regression as the relationship is non-linear. We can take a linear approximation by a Taylor expansion around some initial guess for the breakpoint, $\psi^{(0)}$, 

\begin{equation}
    y \approx \alpha x + c + \beta (x - \psi^{(0)}) H (x - \psi^{(0)}) - \beta (\psi - \psi^{(0)}) H(x - \psi^{(0)}) + \zeta \,. \label{eqn:expansion}
\end{equation}


This is now a linear relationship and we can find a new breakpoint estimate, $\psi^{(1)}$, through ordinary linear regression using the `statsmodels` python package [@seabold2010statsmodels]. We iterate in this way until the breakpoint estimate converges, at which point we stop the algorithm. If considering multiple breakpoints, the same approach is followed using a multivariate Taylor expansion around an initial guess for each of the breakpoints. 

Muggeo's iterative algorithm is not guaranteed to converge on a globally optimal solution. Instead, it can converge to a local optimum or diverge. To address this limitation, we also implement bootstrap restarting [@wood2001minimizing], again following Muggeo's approach [@muggeo2008segmented]. The bootstrap restarting algorithm generates a non-parametric bootstrap of the data through resampling, which is then used to find new breakpoint values that may find a better global solution. This is repeated several times to escape local optima.

# Model Selection

The standard algorithm finds a good fit with a given number of breakpoints. In some instances we might not know how many breakpoints to expect in the data. We provide a tool to compare models with different numbers of breakpoints based on minimising the Bayesian Information Criterion [@wit2012all], which takes into account the value of the likelihood function while including a penalty for the number of model parameters, to avoid overfitting. When applied to the example in \autoref{fig:example}, a model with 3 breakpoints is the preferred choice.

# Features

The package includes the following features:

- Standard fit using the iterative method described by Muggeo.
- Bootstrap restarting to escape local optima.
- Bootstrap restarting with randomised initial breakpoint guesses. 
- Calculation of standard errors and confidence intervals.
- Davies hypothesis test for the existence of a breakpoint. 
- Customisable plots of fits.
- Customisable plots of algorithm iterations.
- Printable summary.
- Summary data output.
- Comprehensive tests.
- Model comparision with an unknown number of breakpoints, with the best fit based on the Bayesian information criterion.  

The package can be downloaded through the [Python Package Index](https://pypi.org/project/piecewise-regression/). The full code is publicly available on [github](https://github.com/chasmani/piecewise-regression). Documentation, including an API reference, can be found at [Read The Docs](https://piecewise-regression.readthedocs.io/en/latest/).

# Acknowledgements

I acknowledge support from Thomas Hills. The work was funded by the EPSRC grant for the Mathematics for Real-World Systems CDT at Warwick (grant number EP/L015374/1).

# References
.. piecewise-regression documentation master file, created by
   sphinx-quickstart on Tue Sep 14 14:20:43 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

piecewise-regression
================================================

.. include:: README.rst

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
==========================================================
piecewise-regression (aka segmented regression) in python
==========================================================
:piecewise-regression: fitting straight line models with breakpoints
:Author: Charlie Pilgrim
:Version: 1.2.1
:Github: https://github.com/chasmani/piecewise-regression
:Documentation: https://piecewise-regression.readthedocs.io/en/master/index.html
:Paper: https://joss.theoj.org/papers/10.21105/joss.03859

.. image:: https://github.com/chasmani/piecewise-regression/actions/workflows/python-package.yml/badge.svg
   :target: https://github.com/chasmani/piecewise-regression/actions/workflows/python-package.yml
   :alt: Build Status
.. image:: https://codecov.io/gh/chasmani/piecewise-regression/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/chasmani/piecewise-regression
   :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/piecewise-regression/badge/?version=latest
   :target: https://piecewise-regression.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://badge.fury.io/py/piecewise-regression.svg
   :target: https://badge.fury.io/py/piecewise-regression
   :alt: PyPi location
.. image:: https://joss.theoj.org/papers/b64e5e7d746efc5d91462a51b3fc5bf8/status.svg
   :target: https://joss.theoj.org/papers/b64e5e7d746efc5d91462a51b3fc5bf8
   :alt: Review Status
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/chasmani/piecewise-regresssion/blob/master/LICENSE
   :alt: License information

|

Easy-to-use piecewise regression (aka segmented regression) in Python. For fitting straight lines to data where there are one or more changes in gradient (known as breakpoints). Based on Muggeo's paper "Estimating regression models with unknown break-points" (2003). 

When using the package, please cite the `accompanying paper <https://joss.theoj.org/papers/10.21105/joss.03859>`_.

Example:

.. image:: https://raw.githubusercontent.com/chasmani/piecewise-regression/master/paper/example.png
    :alt: basic-example-plot-github

Code examples below, and more in this `Google Colab Jupyter Notebook <https://colab.research.google.com/drive/1Pwv6LqwZU8Zbl0VZH6cwOTwoRzm3CPPC#offline=true&sandboxMode=true/>`_.

Installation
========================

You can install piecewise-regression using python's `pip package index <https://pypi.org/project/piecewise-regression/>`_ ::

    pip install piecewise-regression

The package is tested on Python 3.7, 3.8 and 3.9.

Getting started
========================

The package requires some x and y data to fit. You need to specify either a) some initial breakpoint guesses as `start_values` or b) how many breakpoints you want to fit as `n_breakpoints` (or both). Here is an elementary example, assuming we already have some data `x` and `y`: ::

	import piecewise_regression
	pw_fit = piecewise_regression.Fit(x, y, n_breakpoints=2)
	pw_fit.summary()

Example
========================
*For demonstration purposes, substitute with your own data to fit.*

1. Start-off generating some data with a breakpoint: ::

	import piecewise_regression
	import numpy as np

	alpha_1 = -4    
	alpha_2 = -2
	constant = 100
	breakpoint_1 = 7
	n_points = 200
	np.random.seed(0)
	xx = np.linspace(0, 20, n_points)
	yy = constant + alpha_1*xx + (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)


2. Fit the model: ::

    # Given some data, fit the model
    pw_fit = piecewise_regression.Fit(xx, yy, start_values=[5], n_breakpoints=1)

    # Print a summary of the fit
    pw_fit.summary()

Example output: ::

	                    Breakpoint Regression Results                     
	====================================================================================================
	No. Observations                      200
	No. Model Parameters                    4
	Degrees of Freedom                    196
	Res. Sum of Squares               193.264
	Total Sum of Squares              46201.8
	R Squared                        0.995817
	Adjusted R Squared               0.995731
	Converged:                           True
	====================================================================================================
	====================================================================================================
	                    Estimate      Std Err            t        P>|t|       [0.025       0.975]
	----------------------------------------------------------------------------------------------------
	const                100.726        0.244       413.63     3.1e-290       100.25       101.21
	alpha1              -4.21998       0.0653      -64.605    4.37e-134      -4.3488      -4.0912
	beta1                2.18914       0.0689       31.788            -       2.0533        2.325
	breakpoint1          6.48706        0.137            -            -       6.2168       6.7573
	----------------------------------------------------------------------------------------------------
	These alphas(gradients of segments) are estimated from betas(change in gradient)
	----------------------------------------------------------------------------------------------------
	alpha2              -2.03084       0.0218      -93.068    3.66e-164      -2.0739      -1.9878
	====================================================================================================

	Davies test for existence of at least 1 breakpoint: p=5.13032e-295 (e.g. p<0.05 means reject null hypothesis of no breakpoints at 5% significance)

This includes estimates for all the model variables, along with confidence intervals. The Davies test is a hypothesis test for the existence of at least one breakpoint, against the null hypothesis of no breakpoints.  

3. Optional: Plotting the data and model results: ::

	import matplotlib.pyplot as plt

	# Plot the data, fit, breakpoints and confidence intervals
	pw_fit.plot_data(color="grey", s=20)
	# Pass in standard matplotlib keywords to control any of the plots
	pw_fit.plot_fit(color="red", linewidth=4) 
	pw_fit.plot_breakpoints()
	pw_fit.plot_breakpoint_confidence_intervals()
	plt.xlabel("x")
	plt.ylabel("y")
	plt.show()
	plt.close()

.. image:: https://raw.githubusercontent.com/chasmani/piecewise-regression/master/paper/example2.png
    :alt: fit-example-plot-github


You can extract data as well: ::

	# Get the key results of the fit 
	pw_results = pw_fit.get_results()
	pw_estimates = pw_results["estimates"]


How It Works
======================

The package implements Muggeo's iterative algorithm (Muggeo "Estimating regression models with unknown break-points" (2003)) to find breakpoints quickly. This method simultaneously fits breakpoint positions and the linear models for the different fit segments, and it gives confidence intervals for all the model estimates. See the accompanying paper for more details.

Muggeo's method doesn't always converge on the best solution - sometimes, it finds a locally optimal solution or doesn't converge at all. For this reason, the Fit method also implements a process called bootstrap restarting which involves taking a bootstrapped resample of the data to try to find a better solution. The number of times this process runs can be controlled with n_boot. To run the Fit without bootstrap restarting, set ``n_boot=0``.

If you do not have (or do not want to use) initial guesses for the number of breakpoints, you can set it to ``n_breakpoints=3``, and the algorithm will randomly generate start_values. With a 50% chance, the bootstrap restarting algorithm will either use the best currently converged breakpoints or randomly generate new ``start_values``, escaping the local optima in two ways in order to find better global optima. 

As is often the case with fitting non-linear models, even with these measures, the algorithm is not guaranteed to converge to a global optimum. However, increasing ``n_boot`` raises the probability of global convergence at the cost of computation time.


Model Selection
==========================

In addition to the main Fit tool, the package also offers a ModelSelection option based on the Bayesian Information Criterion (BIC). This additional tool is experimental and not as thorough as the main Fit function. In particular, the models are generated with random start_values, which can influence the model fit and give different values for the BIC. The tool can help explore other possible models but should not be used to choose the best model at this time. ::

	ms = piecewise_regression.ModelSelection(x, y, max_breakpoints=6)

Example output: ::

	                 Breakpoint Model Comparision Results                 
	====================================================================================================
	n_breakpoints            BIC    converged          RSS 
	----------------------------------------------------------------------------------------------------
	0                     421.09         True       1557.4 
	1                     14.342         True       193.26 
	2                     22.825         True       191.23 
	3                     24.169         True       182.59 
	4                     29.374         True       177.73 
	5                                   False              
	6                                   False              

	Minimum BIC (Bayesian Information Criterion) suggests the best model 

The data of the model fits can be accessed in ::

    ms.models 

For a robust comparision, you could run the ModelSelection tools many times and take the lowest BIC for each model. 


Testing
============

The package includes comprehensive tests.

To run all tests, from the main directory run (requires the pytest library): ::
	
	pytest

To get code coverage, run (requires pytest and pytest-cov libraries): ::

	pytest --cov=./

There are also a series of simulation tests that check the estimates have realistic confidence intervals, and the Davies test gives realistic p-values. These can be found in the folder "tests-manual". 

Requirements
=============

See requirements.txt for specific version numbers. Required packages, and their uses are:

- matplotlib for plotting.
- numpy for simple data handling and data transformations.  
- scipy for statistical tests including using t-distributions and Gaussians. 
- statsmodels for performing ordinary least squares.

Community Guidelines and Contributing
===================================================

We welcome community participation!

Sourced from `Open Source Guide: How to contribute. <https://opensource.guide/how-to-contribute/>`_

**Open an issue in the following situations:**

- Report an error you can’t solve yourself
- Discuss a high-level topic or idea (for example, community, vision or policies)
- Propose a new feature or other project ideas

**Tips for communicating on issues:**

- If you see an open issue that you want to tackle, comment on the issue to let people know you’re on it. That way, people are less likely to duplicate your work.
- If an issue was opened a while ago, it’s possible that it’s being addressed somewhere else, or has already been resolved, so comment to ask for confirmation before starting work.
- If you opened an issue, but figured out the answer later on your own, comment on the issue to let people know, then close the issue. Even documenting that outcome is a contribution to the project.

**Open a pull request in the following situations:**

- Submit trivial fixes (for example, a typo, a broken link or an obvious error)
- Start work on a contribution that was already asked for, or that you’ve already discussed, in an issue

**Tips for submitting PRs:** 

- Reference any relevant issues or supporting documentation in your PR (for example, “Closes #37.”)
- Include screenshots of the before and after if your changes include differences in HTML/CSS. Drag and drop the images into the body of your pull request.
- Test your changes by running them against any existing tests if they exist and create new ones when needed. Whether tests exist or not, make sure your changes don’t break the existing project.
- Contribute in the style of the project to the best of your abilities.

Installing From Source
===========================

To install from source: ::

	git clone https://github.com/chasmani/piecewise-regression
	cd piecewise_regression
	python3 setup.py install --user


Documentation
==============
`Full docs, including an API reference. <https://piecewise-regression.readthedocs.io/en/latest/>`_





API
================================================

Main
----------
The main module includes the Fit function, which runs the bootstrap restarting algorithm.

.. automodule:: piecewise_regression.main
	:members:

Model selection
---------------------
The model selection module is experimental. It compares models with different `n_breakpoints` using the Bayesian Information Criterion.

.. automodule:: piecewise_regression.model_selection
	:members:

Davies test
----------------------
Implements the Davies hypothesis test for existence of at least one breakpoint.

.. automodule:: piecewise_regression.davies
	:members:

