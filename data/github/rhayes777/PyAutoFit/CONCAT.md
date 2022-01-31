# Contributing

Contributions are welcome and greatly appreciated!

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/rhayes777/PyAutoFit/issues

If you are playing with the PyAutoFit library and find a bug, please
reporting it including:

* Your operating system name and version.
* Any details about your Python environment.
* Detailed steps to reproduce the bug.

### Propose New `NonLinearSearch` or Features

The best way to send feedback is to open an issue at
https://github.com/rhayes777/PyAutoFit/issues
with tag *enhancement*.

If you are proposing a new `NonLinearSearch` or a new feature:

* Explain in detail how it should work.
* Keep the scope as narrow as possible, to make it easier to implement.

### Implement `NonLinearSearch` or Features
Look through the Git issues for operator or feature requests.
Anything tagged with *enhancement* is open to whoever wants to
implement it.

### Add Examples or improve Documentation
Writing new features is not the only way to get involved and
contribute. Create examples with existing non-linear searches as well 
as improving the documentation of existing operators is as important
as making new non-linear searches and very much encouraged.


## Getting Started to contribute

Ready to contribute?

1. Follow the installation instructions for installing **PyAutoFit** from source root on our 
[readthedocs](https://pyautofit.readthedocs.io/en/latest/general/installation.html#forking-cloning>).

2. Create a branch for local development:
    ```
    git checkout -b name-of-your-branch
    ```
    Now you can make your changes locally.

3. When you're done making changes, check that old and new tests pass
succesfully:
    ```
    cd PyAutoFit/test_autofit
    python3 -m pytest
    ```

4. Commit your changes and push your branch to GitLab::
    ```
    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-branch
    ```
    Remember to add ``-u`` when pushing the branch for the first time.

5. Submit a pull request through the GitHub website.


### Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include new tests for all the core routines that have been developed.
2. If the pull request adds functionality, the docs should be updated accordingly.
---
title: "`PyAutoFit`: A Classy Probabilistic Programming Language for Model Composition and Fitting"
tags:
  - Python
  - statistics
  - Bayesian inference
  - probabilistic programming
  - model fitting
authors:
  - name: James. W. Nightingale
    orcid: 0000-0002-8987-7401
    affiliation: 1
  - name: Richard G. Hayes
    affiliation: 1
  - name: Matthew Griffiths
    orcid: 0000-0002-2553-2447
    affiliation: 2 
affiliations:
  - name: Institute for Computational Cosmology, Stockton Rd, Durham, United Kingdom, DH1 3LE
    index: 1
  - name: ConcR Ltd, London, UK
    index: 2    
date: 17 July 2020
codeRepository: https://github.com/rhayes777/PyAutoFit
license: MIT
bibliography: paper.bib
---

# Summary

A major trend in academia and data science is the rapid adoption of Bayesian statistics for data analysis and modeling, 
leading to the development of probabilistic programming languages (PPL). A PPL provides a framework that allows users 
to easily specify a probabilistic model and perform inference automatically. `PyAutoFit` is a Python-based PPL which 
interfaces with all aspects of the modeling (e.g., the model, data, fitting procedure, visualization, results) and
therefore provides complete management of every aspect of modeling. This includes composing high-dimensionality models 
from individual model components, customizing the fitting procedure and performing data augmentation before a model-fit. 
Advanced features include database tools for analysing large suites of modeling results and exploiting domain-specific 
knowledge of a problem via non-linear search chaining. Accompanying `PyAutoFit` is the [autofit workspace](https://github.com/Jammy2211/autofit_workspace), 
which includes example scripts and the `HowToFit` lecture series which introduces non-experts to model-fitting and 
provides a guide on how to begin a project using `PyAutoFit`. Readers can try `PyAutoFit` right now by 
going to [the introduction Jupyter notebook on Binder](https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD)
or checkout our [readthedocs](https://pyautofit.readthedocs.io/en/latest/) for a complete overview
of **PyAutoFit**'s features.

# Background of Probabilistic Programming

Probabilistic programming languages (PPLs) have enabled contemporary statistical inference techniques to be applied 
to a diverse range of problems across academia and industry. Packages such as PyMC3 [@Salvatier2016], 
Pyro [@Bingham2019] and STAN [@Carpenter2017] offer general-purpose frameworks where users can specify a generative 
model and fit it to data using a variety of non-linear fitting techniques. Each package is specialized to problems 
of a certain nature, with many focused on problems like generalized linear modeling or determining the 
distribution(s) from which the data was drawn. For these problems the model is typically composed of the equations and
distributions that are fitted to the data, which are easily expressed syntactically such that the PPL API offers an 
expressive way to define the model and extensions can be implemented in an intuitive and straightforward way.

# Statement of Need

`PyAutoFit` is a PPL whose core design is providing a direct interface with the model, data, fitting procedure and 
results, allowing it to provide comprehensive management of many different aspects of model-fitting. **PyAutoFit** began 
as an Astronomy project for fitting large imaging datasets of galaxies, after the developers found that existing PPLs 
were not suited to the type of model fitting problems Astronomers faced. This includes efficiently analysing large and 
homogenous datasets with an identical model fitting procedure, making it straight forward to fit many models to 
large datasets with streamlined model comparison and massively parallel support for problems where an expensive 
likelihood function means run-times can be of order days or longer. More recent development has generalized `PyAutoFit`, 
making it suitable to a broader range of model-fitting problems.

# Software Description

To compose a model with `PyAutoFit` model components are written as Python classes, allowing `PyAutoFit` to 
define the model and associated parameters in an expressive way that is tied to the modeling software's API. A 
model fit then requires that a `PyAutoFit` `Analysis` class is written, which combines the data, model and likelihood 
function and defines how the model-fit is performed using a `NonLinearSearch`. The `NonLinearSearch`
procedure is defined using an external inference library such as `dynesty` [@dynesty], `emcee` [@emcee]
or `PySwarms` [@pyswarms]. 

The `Analysis` class provides a model specific interface between `PyAutoFit` and the modeling software, allowing it 
to handle the 'heavy lifting' that comes with writing model-fitting software. This includes interfacing with the 
non-linear search, outputting results in a structured path format and model-specific visualization during and 
after the non-linear search. Results are output in a database structure that allows the `Aggregator` tool to load 
results post-analysis via a Python script or Jupyter notebook. This includes methods for summarizing the results of 
every fit, filtering results to inspect subsets of model fits and visualizing results. Results are loaded as `Python` 
generators, ensuring the `Aggregator` can be used to interpret large files in a memory efficient way. `PyAutoFit` is 
therefore suited to 'big data' problems where independent fits to large homogeneous data-sets using an identical 
model-fitting procedure are performed. 

# Model Abstraction and Composition

For many modeling problems the model comprises abstract model components representing objects or processes in a 
physical system. For example, galaxy morphology studies in astrophysics where model components represent the light 
profile of stars [@Haussler2013; @Nightingale2019]. For these problems the likelihood function is typically a 
sequence of numerical processes (e.g., convolutions, Fourier transforms, linear algebra) and extensions to the model 
often requires the addition of new model components in a way that is non-trivially included in the fitting process 
and likelihood function. Existing PPLs have tools for these problems, for example 'black-box' likelihood functions 
in PyMC3. However, these solutions decouple model composition from the data and fitting procedure, making the model 
less expressive, restricting model customization and reducing flexibility in how the model-fit is performed.

By writing model components as Python classes, the model and its associated parameters are defined in an expressive 
way that is tied to the modeling software’s API. Model composition with `PyAutoFit` allows complex models to be built 
from these individual components, abstracting the details of how they change model-fitting procedure from the user. 
Models can be fully customized, allowing adjustment of individual parameter priors, the fixing or coupling of 
parameters between model components and removing regions of parameter space via parameter assertions. Adding new model 
components to a `PyAutoFit` project is straightforward, whereby adding a new Python class means it works within 
the entire modeling framework. `PyAutoFit` is therefore ideal for problems where there is a desire to compose, fit and 
compare many similar (but slightly different) models to a single dataset, with the `Aggregator` including tools to 
facilitate this. 

For many model fitting problems, domain specific knowledge of the model can be exploited to speed up the non-linear 
search and ensure it locates the global maximum likelihood solution. For example, initial fits can be performed 
using simplified model parameterizations, augmented datasets and faster non-linear fitting techniques. Through 
experience users may know that certain model components share minimal covariance, meaning that separate fits to each 
model component (in parameter spaces of reduced dimensionality) can be performed before fitting them simultaneously. 
The results of these simplified fits can then be used to initialize fits using a higher dimensionality model. 
Breaking down a model-fit in this way uses `PyAutoFit`'s non-linear search chaining, which granularizes the non-linear 
fitting procedure into a series of linked non-linear searches. Initial model-fits are followed by fits that gradually 
increase the model complexity, using the information gained throughout the pipeline to guide each `NonLinearSearch` 
and thus enable accurate fitting of models of arbitrary complexity.

# History

`PyAutoFit` is a generalization of [PyAutoLens](https://github.com/Jammy2211/PyAutoLens), an Astronomy package 
developed to analyse images of gravitationally lensed galaxies. Modeling gravitational lenses historically requires 
large amounts of human time and supervision, an approach which does not scale to the incoming samples of 100000 objects. 
Domain exploitation enabled full automation of the lens modeling procedure [@Nightingale2015; @Nightingale2018], with 
model customization and the aggregator enabling one to fit large datasets with many different models. More 
recently, `PyAutoFit` has been applied to calibrating radiation damage to charge coupled imaging devices and a model 
of cancer tumour growth.  
 
# Workspace and HowToFit Tutorials

`PyAutoFit` is distributed with the [autofit workspace](https://github.com/Jammy2211/autofit_workspace), which 
contains example scripts for composing a model, performing a fit, using the `Aggregator` and `PyAutoFit`'s advanced 
statistical inference methods. Also included are the `HowToFit` tutorials, a series of Jupyter notebooks aimed at 
non-experts, introducing them to model-fitting and Bayesian inference. They teach users how to write model-components 
and `Analysis` classes in `PyAutoFit`, use these to fit a dataset and interpret the model-fitting results. The lectures 
are available on our [Binder](https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD) and may therefore be 
taken without a local `PyAutoFit` installation.

# Software Citations

`PyAutoFit` is written in Python 3.6+ [@python] and uses the following software packages:

- `corner.py` https://github.com/dfm/corner.py [@corner]
- `dynesty` https://github.com/joshspeagle/dynesty [@dynesty]
- `emcee` https://github.com/dfm/emcee [@emcee]
- `matplotlib` https://github.com/matplotlib/matplotlib [@matplotlib]
- `NumPy` https://github.com/numpy/numpy [@numpy]
- `PyMultiNest` https://github.com/JohannesBuchner/PyMultiNest [@multinest] [@pymultinest]
- `PySwarms` https://github.com/ljvmiranda921/pyswarms [@pyswarms]
- `Scipy` https://github.com/scipy/scipy [@scipy]

# Related Probabilistic Programming Languages

- `PyMC3` https://github.com/pymc-devs/pymc3 [@Salvatier2016]
- `Pyro` https://github.com/pyro-ppl/pyro [@Bingham2019]
- `STAN` https://github.com/stan-dev/stan [@Carpenter2017]
- `TensorFlow Probability` https://github.com/tensorflow/probability [@tensorflow]
- `uravu` https://github.com/arm61/uravu [@uravu]

# Acknowledgements

JWN and RJM are supported by the UK Space Agency, through grant ST/V001582/1, and by InnovateUK through grant TS/V002856/1. RGH is supported by STFC Opportunities grant ST/T002565/1.
This work used the DiRAC@Durham facility managed by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC capital grants ST/K00042X/1, ST/P002293/1, ST/R002371/1 and ST/S002502/1, Durham University and STFC operations grant ST/R000832/1. DiRAC is part of the National e-Infrastructure.

# References
PyAutoFit JOSS Paper
====================

Paper accompanying [PyAutoFit](https://github.com/rhayes777/PyAutoFit) for submission to the Journal of Open Source 
Software (JOSS).I**nesrt in the main body of the paper:**

We use the probabilistic programming language `PyAutoFit` https://github.com/rhayes777/PyAutoFit) [@pyautofit] to...

**At the end of the paper (delete as appropriate, see https://pyautofit.readthedocs.io/en/latest/general/citations.html):**

# Software Citations

This work uses the following software packages:

- `corner.py` https://github.com/dfm/corner.py [@corner]
- `dynesty` https://github.com/joshspeagle/dynesty [@dynesty]
- `emcee` https://github.com/dfm/emcee [@emcee]
- `matplotlib` https://github.com/matplotlib/matplotlib [@matplotlib]
- `NumPy` https://github.com/numpy/numpy [@numpy]
- `PyAutoFit` https://github.com/rhayes777/PyAutoFit [@pyautofit]
- `PyMultiNest` https://github.com/JohannesBuchner/PyMultiNest [@multinest] [@pymultinest]
- `PySwarms` https://github.com/ljvmiranda921/pyswarms [@pyswarms]
- `Python` https://www.python.org/ [@python]
- `Scipy` https://github.com/scipy/scipy [@scipy]
- `SQLite` https://www.sqlite.org/index.html [@sqlite]
- `UltraNest` https://github.com/JohannesBuchner/UltraNest [@ultranest]
- `Zeus` https://github.com/minaskar/zeus [@zeus1] [@zeus2]PyAutoFit: Classy Probabilistic Programming
===========================================

.. |binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.02550/status.svg
   :target: https://doi.org/10.21105/joss.02550

|binder| |JOSS|

`Installation Guide <https://pyautofit.readthedocs.io/en/latest/installation/overview.html>`_ |
`readthedocs <https://pyautofit.readthedocs.io/en/latest/index.html>`_ |
`Introduction on Binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=introduction.ipynb>`_ |
`HowToFit <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_

PyAutoFit is a Python based probabilistic programming language for the fully Bayesian analysis of extremely large
datasets which:

- Makes it simple to compose and fit mult-level models using a range of Bayesian inference libraries, such as `emcee <https://github.com/dfm/emcee>`_ and `dynesty <https://github.com/joshspeagle/dynesty>`_.

- Handles the 'heavy lifting' that comes with model-fitting, including model composition & customization, outputting results, model-specific visualization and posterior analysis.

- Is built for *big-data* analysis, whereby results are output as a sqlite database which can be queried after model-fitting is complete.

**PyAutoFit** supports advanced statistical methods such as `empirical Bayes <https://pyautofit.readthedocs.io/en/latest/features/empirical_bayes.html>`_, `sensitivity mapping <https://pyautofit.readthedocs.io/en/latest/features/sensitivity_mapping.html>`_ and `massively parallel model-fits <https://pyautofit.readthedocs.io/en/latest/features/search_grid_search.html>`_ .

Getting Started
---------------

The following links are useful for new starters:

- `The introduction Jupyter Notebook on Binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=introduction.ipynb>`_, where you can try **PyAutoFit** in a web browser (without installation).

- `The PyAutoFit readthedocs <https://pyautofit.readthedocs.io/en/latest>`_, which includes an `installation guide <https://pyautofit.readthedocs.io/en/latest/installation/overview.html>`_ and an overview of **PyAutoFit**'s core features.

- `The autofit_workspace GitHub repository <https://github.com/Jammy2211/autofit_workspace>`_, which includes example scripts and the `HowToFit Jupyter notebook tutorials <https://github.com/Jammy2211/autofit_workspace/tree/master/notebooks/howtofit>`_ which give new users a step-by-step introduction to **PyAutoFit**.

Why PyAutoFit?
--------------

**PyAutoFit** began as an Astronomy project for fitting large imaging datasets of galaxies after the developers found that existing PPLs
(e.g., `PyMC3 <https://github.com/pymc-devs/pymc3>`_, `Pyro <https://github.com/pyro-ppl/pyro>`_, `STAN <https://github.com/stan-dev/stan>`_)
were not suited to the model fitting problems many Astronomers faced. This includes:

- Efficiently analysing large and homogenous datasets with an identical model fitting procedure, with tools for processing the large libraries of results output.

- Problems where likelihood evaluations are expensive (e.g. run times of days per model-fit), necessitating highly customizable model-fitting pipelines with support for massively parallel computing.

- Fitting many different models to the same dataset with tools that streamline model comparison.

If these challenges sound familiar, then **PyAutoFit** may be the right software for your model-fitting needs!

API Overview
------------

To illustrate the **PyAutoFit** API, we'll use an illustrative toy model of fitting a one-dimensional Gaussian to
noisy 1D data. Here's the ``data`` (black) and the model (red) we'll fit:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/files/toy_model_fit.png
  :width: 400
  :alt: Alternative text

We define our model, a 1D Gaussian by writing a Python class using the format below:

.. code-block:: python

    class Gaussian:

        def __init__(
            self,
            centre=0.0,        # <- PyAutoFit recognises these
            normalization=0.1, # <- constructor arguments are
            sigma=0.01,        # <- the Gaussian's parameters.
        ):
            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

        """
        An instance of the Gaussian class will be available during model fitting.

        This method will be used to fit the model to data and compute a likelihood.
        """

        def profile_1d_via_xvalues_from(self, xvalues):

            transformed_xvalues = xvalues - self.centre

            return (self.normalization / (self.sigma * (2.0 * np.pi) ** 0.5)) * \
                    np.exp(-0.5 * (transformed_xvalues / self.sigma) ** 2.0)

**PyAutoFit** recognises that this Gaussian may be treated as a model component whose parameters can be fitted for via
a non-linear search like `emcee <https://github.com/dfm/emcee>`_.

To fit this Gaussian to the ``data`` we create an Analysis object, which gives **PyAutoFit** the ``data`` and a
``log_likelihood_function`` describing how to fit the ``data`` with the model:

.. code-block:: python

    class Analysis(af.Analysis):

        def __init__(self, data, noise_map):

            self.data = data
            self.noise_map = noise_map

        def log_likelihood_function(self, instance):

            """
            The 'instance' that comes into this method is an instance of the Gaussian class
            above, with the parameters set to values chosen by the non-linear search.
            """

            print("Gaussian Instance:")
            print("Centre = ", instance.centre)
            print("normalization = ", instance.normalization)
            print("Sigma = ", instance.sigma)

            """
            We fit the ``data`` with the Gaussian instance, using its
            "profile_1d_via_xvalues_from" function to create the model data.
            """

            xvalues = np.arange(self.data.shape[0])

            model_data = instance.profile_1d_via_xvalues_from(xvalues=xvalues)
            residual_map = self.data - model_data
            chi_squared_map = (residual_map / self.noise_map) ** 2.0
            log_likelihood = -0.5 * sum(chi_squared_map)

            return log_likelihood

We can now fit our model to the ``data`` using a non-linear search:

.. code-block:: python

    model = af.Model(Gaussian)

    analysis = Analysis(data=data, noise_map=noise_map)

    emcee = af.Emcee(nwalkers=50, nsteps=2000)

    result = emcee.fit(model=model, analysis=analysis)

The ``result`` contains information on the model-fit, for example the parameter samples, maximum log likelihood
model and marginalized probability density functions.

Support
-------

Support for installation issues, help with Fit modeling and using **PyAutoFit** is available by
`raising an issue on the GitHub issues page <https://github.com/rhayes777/PyAutoFit/issues>`_.

We also offer support on the **PyAutoFit** `Slack channel <https://pyautoFit.slack.com/>`_, where we also provide the 
latest updates on **PyAutoFit**. Slack is invitation-only, so if you'd like to join send 
an `email <https://github.com/Jammy2211>`_ requesting an invite.Probabilistic Programming
=========================

Probabilistic programming languages provide a framework that allows users to easily specify a probabilistic
model and perform inference automatically. PyAutoFit is a Python based probabilistic programming language for the
fully Bayesian analysis of extremely large datasets which:

- Makes it simple to compose and fit multi-level models using a range of Bayesian inference libraries, such as `emcee <https://github.com/dfm/emcee>`_ and `dynesty <https://github.com/joshspeagle/dynesty>`_.

- Handles the 'heavy lifting' that comes with model-fitting, including model composition & customization, outputting results, visualization and parameter inference.

- Is built for *big-data* analysis, whereby results are output as a sqlite database which can be queried after model-fitting is complete.

**PyAutoFit** supports advanced statistical methods such as `massively parallel non-linear search grid-searches <https://pyautofit.readthedocs.io/en/latest/features/search_grid_search.html>`_, `chaining together model-fits <https://pyautofit.readthedocs.io/en/latest/features/search_chaining.html>`_  and `sensitivity mapping <https://pyautofit.readthedocs.io/en/latest/features/sensitivity_mapping.html>`_.

Try it now
----------

You can try **PyAutoFit** now by going to the `introduction Jupyter Notebook on our
Binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=introduction.ipynb>`_, which runs
**PyAutoFit** in a web browser without installation.

Why PyAutoFit?
--------------

**PyAutoFit** is developed by Astronomers for fitting large imaging datasets of galaxies. We found that existing
probabilistic programming languages (e.g `PyMC3 <https://github.com/pymc-devs/pymc3>`_, `Pyro <https://github.com/pyro-ppl/pyro>`_,
`STAN <https://github.com/stan-dev/stan>`_) were not suited to the type of model fitting problems Astronomers faced,
for example:

- Fitting large and homogenous datasets with an identical model fitting procedure, with tools for processing the large libraries of results output.

- Problems where likelihood evaluations are expensive, leading to run times of days per fit and necessitating support for massively parallel computing.

- Fitting many different models to the same dataset with tools that streamline model comparison.

How does PyAutoFit Work?
========================

Model components are written as Python classes, allowing **PyAutoFit** to define the *model* and
associated *parameters* in an expressive way that is tied to the modeling software's API. Here is a simple example of
how a *model* representing a 1D Gaussian is written:

.. code-block:: python

    class Gaussian:

        def __init__(
            self,
            centre=0.0,     # <- PyAutoFit recognises these
            normalization=0.1,  # <- constructor arguments are
            sigma=0.01,     # <- the Gaussian's parameters.
        ):
            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

        """
        An instance of the Gaussian class will be available during model fitting.

        This method will be used to fit the model to ``data`` and compute a likelihood.
        """

        def profile_1d_via_xvalues_from(self, xvalues):

            transformed_xvalues = xvalues - self.centre

            return (self.normalization / (self.sigma * (2.0 * np.pi) ** 0.5)) * \
                    np.exp(-0.5 * transformed_xvalues / self.sigma)

A model-fit requires that a **PyAutoFit** ``Analysis`` class is written, which combines the data and model via
likelihood function:

.. code-block:: python

    class Analysis(af.Analysis):

        def __init__(self, data, noise_map):

            self.data = data
            self.noise_map = noise_map

        def log_likelihood_function(self, instance):

            """
            The 'instance' that comes into this method is an instance of the Gaussian class
            above, with the *parameters* set to values chosen by the non-linear search.
            """

            print("Gaussian Instance:")
            print("Centre = ", instance.centre)
            print("normalization = ", instance.normalization)
            print("Sigma = ", instance.sigma)

            """
            We fit the ``data`` with the Gaussian instance, using its
            "profile_1d_via_xvalues_from" function to create the model data.
            """

            xvalues = np.arange(self.data.shape[0])

            model_data = instance.profile_1d_via_xvalues_from(xvalues=xvalues)
            residual_map = self.data - model_data
            chi_squared_map = (residual_map / self.noise_map) ** 2.0
            log_likelihood = -0.5 * sum(chi_squared_map)

            return log_likelihood

The ``Analysis`` class provides a model specific interface between **PyAutoFit** and the modeling software, allowing
it to handle the 'heavy lifting' that comes with writing *model-fitting* software. This includes interfacing with the
non-linear search, model-specific visualization during and outputting results to a queryable sqlite database.

Performing a fit with a non-linear search, for example ``emcee``, is performed as follows:

.. code-block:: python

    model = af.Model(Gaussian)

    analysis = Analysis(data=data, noise_map=noise_map)

    emcee = af.Emcee(name="example_search", nwalkers=50, nsteps=2000)

    result = emcee.fit(model=model, analysis=analysis)

The ``result`` contains information on the model-fit, for example the parameter samples, maximum log likelihood
model and marginalized probability density functions.

Model Abstraction and Composition
=================================

For many model fitting problems the model comprises abstract *model components* representing objects or processes in a
physical system. For example, our child project `PyAutoLens <https://github.com/Jammy2211/PyAutoLens>`_,  where
*model components* represent the light and mass of galaxies. For these problems the likelihood function is typically a
sequence of numerical processes (e.g. convolutions, Fourier transforms, linear algebra) and extensions to the *model* 
often requires the addition of new *model components* in a way that is non-trivially included in the fitting process
and likelihood function. Existing PPLs have tools for these problems, however they decouple *model composition* from the
data and fitting procedure, making the *model* less expressive, restricting *model customization* and reducing
flexibility in how the *model-fit* is performed.

By writing *model components* as ``Python`` classes, the *model* and its associated *parameters* are defined in an
expressive way that is tied to the modeling software’s API. *Model composition* with **PyAutoFit** allows complex
*models* to be built from these individual components, abstracting the details of how they change *model-fitting*
procedure from the user. *Models* can be fully customized, allowing adjustment of individual parameter priors, the
fixing or coupling of parameters between *model components* and removing regions of parameter space via parameter
assertions. Adding new *model components* to a **PyAutoFit** project is straightforward, whereby adding a new
``Python`` class means it works within the entire modeling framework. **PyAutoFit** is therefore ideal for
problems where there is a desire to *compose*, *fit* and *compare* many similar (but slightly different) models to a
single dataset, with **Database** tools available to facilitate this.

The `overview section <https://pyautofit.readthedocs.io/en/latest/overview/model_fit.html>`_ gives a run-down of
**PyAutoFit**'s core features and the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_
provides new users with a more detailed introduction to **PyAutoFit**.

.. toctree::
   :caption: Overview:
   :maxdepth: 1
   :hidden:

   overview/model_fit
   overview/model_complex
   overview/non_linear_search
   overview/result
   overview/multi_level

.. toctree::
   :caption: Installation:
   :maxdepth: 1
   :hidden:

   installation/overview
   installation/conda
   installation/pip
   installation/source
   installation/troubleshooting

.. toctree::
   :caption: General:
   :maxdepth: 1
   :hidden:

   general/workspace
   general/cookbook
   general/adding_a_model_component
   general/configs
   general/roadmap
   general/software
   general/citations
   general/credits

.. toctree::
   :caption: Tutorials:
   :maxdepth: 1
   :hidden:

   howtofit/howtofit
   howtofit/chapter_1_introduction
   howtofit/chapter_database
   howtofit/chapter_graphical_models

.. toctree::
   :caption: API Reference:
   :maxdepth: 1
   :hidden:

   api/api

.. toctree::
   :caption: Graphical models:
   :maxdepth: 1
   :hidden:

   graphical/multiple_datasets
   graphical/expectation_propagation

.. toctree::
   :caption: Features:
   :maxdepth: 1
   :hidden:

   features/database
   features/search_grid_search
   features/empirical_bayes
   features/sensitivity_mapping.. _references:

Citations & References
======================

The bibtex entries for **PyAutoFit** and its affiliated software packages can be found
`here <https://github.com/rhayes777/PyAutoFit/blob/master/files/citations.bib>`_, with example text for citing **PyAutoFit**
in `.tex format here <https://github.com/rhayes777/PyAutoFit/blob/master/files/citation.tex>`_ format here and
`.md format here <https://github.com/rhayes777/PyAutoFit/blob/master/files/citations.md>`_. As shown in the examples, we
would greatly appreciate it if you mention **PyAutoFit** by name and include a link to our GitHub page!

**PyAutoFit** is published in the `Journal of Open Source Software <https://joss.theoj.org/papers/10.21105/joss.02550#>`_ and its
entry in the above .bib file is under the key ``pyautofit``... _workspace:

Workspace Tour
==============

You should have downloaded and configured the `autofit_workspace <https://github.com/Jammy2211/autofit_workspace>`_
when you installed **PyAutoFit**. If you didn't, checkout the
`installation instructions <https://pyautofit.readthedocs.io/en/latest/general/installation.html#installation-with-pip>`_
for how to downloaded and configure the workspace.

New users should begin by checking out the following parts of the workspace.

HowToFit
--------

The **HowToFit** lecture series are a collection of Jupyter notebooks describing how to build a **PyAutoFit** model
fitting project and giving illustrations of different statistical methods and techniiques.

Checkout the
`tutorials section <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_ for a
full description of the lectures and online examples of every notebook.

Scripts / Notebooks
-------------------

There are numerous example describing how perform model-fitting with **PyAutoFit** and providing an overview of its
advanced model-fitting features. All examples are provided as Python scripts and Jupyter notebooks.

A full description of the scripts available is given on
the `autofit workspace GitHub page <https://github.com/Jammy2211/autofit_workspace>`_.

Config
------

Here, you'll find the configuration files used by **PyAutoFit** which customize:

    - The default settings used by every non-linear search.
    - Example priors and notation configs which associate model-component with model-fitting.
    - The ``general.ini`` config which customizes other aspects of **PyAutoFit**.

Checkout the `configuration <https://pyautofit.readthedocs.io/en/latest/general/installation.html#installation-with-pip>`_
section of the ``readthedocs`` for a complete description of every configuration file.

Dataset
-------

This folder stores the example dataset's used in examples in the workspace.

Output
------

The folder where the model-fitting results of a non-linear search are stored... _cookbook:

Cookbook
========

This cookbook therefore provides an API reference for model composition, in particular:

- Creating individual model components via ``af.Model``.
- Creating collections of model components as a single model using ``af.Collection``.
- Creating multi-level models using hierarchies of Python classes.

Examples using different **PyAutoFit** API's for model composition are provided, which produce more concise and
readable code for different use-cases.

Python Class Template
---------------------

A model component is written as a Python class using the following format:

- The name of the class is the name of the model component, in this case, "Gaussian".

- The input arguments of the constructor are the parameters of the mode (here ``centre``, ``normalization`` and ``sigma``).

- The default values of the input arguments tell **PyAutoFit** whether a parameter is a single-valued ``float`` or a multi-valued ``tuple``.

.. code-block:: bash

    class Gaussian:
        def __init__(
            self,
            centre=30.0,  # <- **PyAutoFit** recognises these constructor arguments
            normalization=1.0,  # <- are the Gaussian``s model parameters.
            sigma=5.0,
        ):
            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

Model
-----

To instantiate a Python class as a model component using the ``af.Model`` object:

.. code-block:: bash

    model = af.Model(Gaussian)

To overwrite the priors of one or more parameters from the default value assumed via configuration files:

.. code-block:: bash

    model = af.Model(Gaussian)
    model.centre = af.UniformPrior(lower_limit=0.0, upper_limit=1.0)
    model.normalization = af.LogUniformPrior(lower_limit=1e-4, upper_limit=1e4)
    model.sigma = af.GaussianPrior(mean=0.0, sigma=1.0, lower_limit=0.0, upper_limit=1e5)

To fix a free parameter to a specific value (reducing the dimensionality of parameter space by 1):

.. code-block:: bash

    model = af.Model(Gaussian)
    model.centre = 0.0

To link two parameters together such they always assume the same value (reducing the dimensionality of parameter space by 1):

.. code-block:: bash

    model = af.Model(Gaussian)
    model.centre = model.normalization

Offsets between linked parameters are also possible:

.. code-block:: bash

    model = af.Model(Gaussian)
    model.centre = model.normalization - 1.0
    model.centre = model.normalization + model.sigma

Assertions remove regions of parameter space:

.. code-block:: bash

    model = af.Model(Gaussian)
    model.add_assertion(model.sigma > 5.0)
    model.add_assertion(model.centre > model.normalization)

Model (Alternative API)
-----------------------

The overwriting of priors shown above can be achieved via the following alternative API:

.. code-block:: bash

    model = af.Model(
        Gaussian,
        centre=af.UniformPrior(lower_limit=0.0, upper_limit=1.0),
        normalization=af.LogUniformPrior(lower_limit=1e-4, upper_limit=1e4),
        sigma=af.GaussianPrior(mean=0.0, sigma=1.0),
    )

This API can also be used for fixing a parameter to a certain value:

.. code-block:: bash

    model = af.Model(Gaussian, centre=0.0)

Collection
----------

To instantiate multiple Python classes into a combined model component using ``af.Collection`` and ``af.Model``:

.. code-block:: bash

    gaussian_0 = af.Model(Gaussian)
    gaussian_1 = af.Model(Gaussian)

    model = af.Collection(gaussian_0=gaussian_0, gaussian_1=gaussian_1)

By setting up each ``Model`` first the model can be customized using either of the ``af.Model`` API's shown above:

.. code-block:: bash

    gaussian_0 = af.Model(Gaussian)
    gaussian_0.normalization = 1.0
    gaussian_0.sigma = af.GaussianPrior(mean=0.0, sigma=1.0)

    gaussian_1 = af.Model(
        Gaussian,
        centre=af.UniformPrior(lower_limit=0.0, upper_limit=1.0),
        normalization=af.LogUniformPrior(lower_limit=1e-4, upper_limit=1e4),
        sigma=af.GaussianPrior(mean=0.0, sigma=1.0),
    )

    model = af.Collection(gaussian_0=gaussian_0, gaussian_1=gaussian_1)

Collection (Alternative API)
----------------------------

To create the ``Collection`` in one line of Python by not defining each ``Model`` beforehand:

.. code-block:: bash

    model = af.Collection(gaussian_0=af.Model(Gaussian), gaussian_1=af.Model(Gaussian))

Using this API, the ``af.Model()`` command can be omitted altogether (**PyAutoFit** will automatically determine
the ``Gaussian`` python classes should be set up as ``Model``'s):

.. code-block:: bash

    model = af.Collection(gaussian_0=Gaussian, gaussian_1=Gaussian)

To customize a model using this API the name of the model subcomponents (e.g. ``gaussian_0`` and ``gaussian_1``) are used
to access and customize the parameters.

.. code-block:: bash

    model = af.Collection(gaussian_0=Gaussian, gaussian_1=Gaussian)

    model.gaussian_0.normalization = 1.0
    model.gaussian_0.sigma = af.GaussianPrior(mean=0.0, sigma=1.0)

    model.gaussian_0.centre = model.gaussian_1.centre

    model.gaussian_1.add_assertion(model.gaussian_1.sigma > 5.0)
    model.gaussian_1.centre = model.gaussian_1.normalization - 1.0

Multi-level Models (Advanced)
-----------------------------

A multi-level model component is written as a Python class using the following format:

- The input arguments include one or more optional lists of Python classes that themselves are instantiated as model components.

- Addition parameters specific to the higher level of the model can be included in the constructor (in this example a parameter called the ``higher_level_parameter`` is used).

Like a normal model component, the name of the Python class is the name of the model component, input arguments are
the parameters of the model and default values tell **PyAutoFit** whether a parameter is a single-valued ``float`` or a
multi-valued ``tuple``.

.. code-block:: bash

    class MultiLevelGaussians:

        def __init__(
            self,
            higher_level_parameter=1.0,
            gaussian_list=None,  # This will optionally contain a list of ``af.Model(Gaussian)``'s
        ):

            self.higher_level_parameter = higher_level_parameter

            self.gaussian_list = gaussian_list

This multi-level model is instantiated via the ``af.Model()`` command, which is passed one or more ``Gaussian`` components:

.. code-block:: bash

    multi_level = af.Model(
        MultiLevelGaussians, gaussian_list=[af.Model(Gaussian), af.Model(Gaussian)]
    )

Again, if the ``af.Model()`` on the individual ``Gaussian``'s is omitted they are still created as model components:

.. code-block:: bash

    multi_level = af.Model(MultiLevelGaussians, gaussian_list=[Gaussian, Gaussian])

To customize the higher level parameters of a multi-level the usual ``Model`` API is used:

.. code-block:: bash

    multi_level = af.Model(MultiLevelGaussians, gaussian_list=[Gaussian, Gaussian])

    multi_level.higher_level_parameter = af.UniformPrior(lower_limit=0.0, upper_limit=1.0)

To customize a multi-level model instantiated via lists, each model component is accessed via its index:

.. code-block:: bash

    multi_level = af.Model(MultiLevelGaussians, gaussian_list=[Gaussian, Gaussian])

    multi_level.gaussian_list[0].centre = multi_level.gaussian_list[1].centre

Any combination of the API's shown above can be used for customizing this model:

.. code-block:: bash

    gaussian_0 = af.Model(Gaussian)
    gaussian_1 = af.Model(Gaussian)

    gaussian_0.centre = gaussian_1.centre

    multi_level = af.Model(
        MultiLevelGaussians, gaussian_list=[gaussian_0, gaussian_1, af.Model(Gaussian)]
    )

    multi_level.higher_level_parameter = 1.0
    multi_level.gaussian_list[2].centre = multi_level.gaussian_list[1].centre

Multi-level Models (Alternative API)
------------------------------------

A multi-level model can be instantiated where each model sub-component is setup using a name (as opposed to a list).

This means no list input parameter is required in the Python class of the model component:

.. code-block:: bash

    class MultiLevelGaussians:

        def __init__(self, higher_level_parameter=1.0):

            self.higher_level_parameter = higher_level_parameter

        multi_level = af.Model(MultiLevelGaussians, gaussian_0=Gaussian, gaussian_1=Gaussian)

Each model subcomponent can be customized using its name, analogous to the ``Collection`` API:

.. code-block:: bash

    multi_level = af.Model(MultiLevelGaussians, gaussian_0=Gaussian, gaussian_1=Gaussian)

    multi_level.gaussian_0.centre = multi_level.gaussian_1.centre

Multi-level Model Collections
-----------------------------

Models, multi-level models and collections can be combined to compose models of high complexity:

.. code-block:: bash

    multi_level_0 = af.Model(MultiLevelGaussians, gaussian_0=Gaussian, gaussian_1=Gaussian)

    multi_level_1 = af.Model(
        MultiLevelGaussians, gaussian_0=Gaussian, gaussian_1=Gaussian, gaussian_2=Gaussian
    )

    model = af.Collection(multi_level_0=multi_level_0, multi_level_1=multi_level_1)

    print(model.multi_level_0.gaussian_1.centre)
    print(model.multi_level_1.higher_level_parameter)

Wrap Up
-------

The API described here can be extended in all the ways one would expect.

For example, multi-level models composed of multiple levels are possible:

.. code-block:: bash

    multi_level_x2_model = af.Model(
        MultiLevelGaussians,
        multi_level_0=af.Model(MultiLevelGaussians, gaussian_0=Gaussian),
        multi_level_1=af.Model(MultiLevelGaussians, gaussian_0=Gaussian),
    )

    print(multi_level_x2_model.multi_level_0.gaussian_0.centre).. _credits:

Credits
-------

**Developers:**

`Richard Hayes <https://github.com/rhayes777>`_ - Lead developer

`James Nightingale <https://github.com/Jammy2211>`_ - Lead developer

`Matthew Griffiths <https://github.com/matthewghgriffiths>`_ - Graphical models guru.. _configs:

Configs
=======

The ``autofit_workspace`` includes configuration files that customize the behaviour of the non-linear search's,
visualization and other aspects of **PyAutoFit**. Here, we describe how to configure **PyAutoFit** to use the configs
and describe every configuration file complete with input parameters.

Setup
-----

By default, **PyAutoFit** looks for the config files in a ``config`` folder in the current working directory, which is
why we run autofit scripts from the ``autofit_workspace`` directory.

The configuration path can also be set manually in a script using **PyAutoConf** and the following command (the path
to the ``output`` folder where the results of a non-linear search are stored is also set below):

.. code-block:: bash

    from autoconf import conf

    conf.instance.push(
        config_path="path/to/config", output_path=f"path/to/output"
    )

general.ini
-----------

This config file is found at ``autofit_workspace/config/general.ini`` and contains the following sections and variables:

[output]
    log_file -> str
        The file name the logged output is written to (in the non-linear search output folder).
    log_level -> str
        The level of logging.
    model_results_decimal_places -> int
        The number of decimal places the estimated values and errors of all parameters in the ``model.results`` file are
        output to.
    remove_files -> bool
        If `True`, all output files of a non-linear search (e.g. samples, samples_backup, model.results, images, etc.)
        are deleted once the model-fit has completed.

        A .zip file of all output is always created before files are removed, thus results are not lost with this
        option turned on. If **PyAutoFit** does not find the output files of a model-fit (because they were removed) but
        does find this .zip file, it will unzip the contents and continue the analysis as if the files were
        there all along.

        This feature was implemented because super-computers often have a limit on the number of files allowed per
        user and the large number of files output by **PyAutoFit** can exceed this limit. By removing files the
        number of files is restricted only to the .zip files.
    grid_results_interval -> int
        For a ``GridSearch`` this interval sets after how many samples on the grid output is
        performed for. A ``grid_results_interval`` of -1 turns off output.

non_linear
----------

These config files are found at ``autofit_workspace/config/non_linear`` and they contain the default settings used by
every non-linear search. The ``[search]``, ``[settings]`` and ``[initialize]`` sections of the non-linear configs
contains settings specific to certain non-linear search's, and the documentation for these variables should be found
by inspecting the`API Documentation <https://pyautofit.readthedocs.io/en/latest/api/api.html>`_ of the relevent
non-linear search object.

The following config sections and variables are generic across all non-linear search configs (e.g.
``config/non_linear/nest/DynestyStatic.ini``, ``config/non_linear/mcmc/Emcee.ini``, etc.):

[updates]
   iterations_per_update -> int
        The number of iterations of the non-linear search performed between every 'update', where an update performs
        visualization of the maximum log likelihood model, backing-up of the samples, output of the ``model.results``
        file and logging.
   visualize_every_update -> int
        For every ``visualize_every_update`` updates visualization is performed and output to the hard-disk during the
        non-linear using the maximum log likelihood model. A ``visualization_interval`` of -1 turns off on-the-fly
        visualization.
   backup_every_update -> int
        For every ``backup_every_update`` the results of the non-linear search in the samples foler and backed up into the
        samples_backup folder. A ``backup_every_update`` of -1 turns off backups during the non-linear search (it is still
        performed when the non-linear search terminates).
   model_results_every_update -> int
        For every ``model_results_every_update`` the model.results file is updated with the maximum log likelihood model
        and parameter estimates with errors at 1 an 3 sigma confidence. A ``model_results_every_update`` of -1 turns off
        the model.results file being updated during the model-fit (it is still performed when the non-linear search
        terminates).
   log_every_update -> int
        For every ``log_every_update`` the log file is updated with the output of the Python interpreter. A
        ``log_every_update`` of -1 turns off logging during the model-fit.

[printing]
    silence -> bool
        If `True`, the default print output of the non-linear search is silenced and not printed by the Python
        interpreter.

[parallel]
    number_of_cores -> int
        For non-linear search's that support parallel procesing via the Python ``multiprocesing`` module, the number of
        cores the parallel run uses. If ``number_of_cores=1``, the model-fit is performed in serial omitting the use
        of the ``multiprocessing`` module.

visualize
---------

These config files are found at ``autofit_workspace/config/visualize`` and they contain the default settings used by
visualization in **PyAutoFit**. The ``general.ini`` config contains the following sections and variables:

[general]
    backend -> str
        The ``matploblib backend`` used for visualization (see
        https://gist.github.com/CMCDragonkai/4e9464d9f32f5893d837f3de2c43daa4 for a description of backends).

        If you use an invalid backend for your computer, **PyAutoFit** may crash without an error or reset your machine.
        The following backends have worked for **PyAutoFit** users:

        TKAgg (default)

        Qt5Agg (works on new MACS)

        Qt4Agg

        WXAgg

        WX

        Agg (outputs to .fits / .png but doesn't'display figures during a run on your computer screen)

priors
------

These config files are found at ``autofit_workspace/config/priors`` and they contain the default priors and related
variables for every model-component in a project, using ``.json`` format files (as opposed to ``.ini`` for most config files).

The autofit_workspace`` contains example ``prior`` files for the 1D ``data`` fitting problem. An example entry of the
json configs for the ``sigma`` parameter of the ``Gaussian`` class is as follows:

.. code-block:: bash

    "Gaussian": {
        "sigma": {
            "type": "Uniform",
            "lower_limit": 0.0,
            "upper_limit": 30.0,
            "width_modifier": {
                "type": "Absolute",
                "value": 0.2
            },
            "gaussian_limits": {
                "lower": 0.0,
                "upper": "inf"
            }
        },

The sections of this example config set the following:

json config
    type -> Prior
        The default prior given to this parameter which is used by the non-linear search. In the example above, a
        ``UniformPrior`` is used with ``lower_limit`` of 0.0 and ``upper_limit`` of 30.0. A ``GaussianPrior`` could be used by
        putting "``Gaussian``" in the "``type``" box, with "``mean``" and "``sigma``" used to set the default values. Any prior can be
        set in an analogous fashion (see the example configs).
    width_modifier
        When the results of a search are linked to a subsequent search to set up the priors of its non-linear search,
        this entry describes how the ``Prior`` is passed. For a full description of prior passing, checkout the examples
        in ``autofit_workspace/notebooks/features/search_chaining``.
    gaussian_limits
        When the results of a search are linked to a subsequent search, they are passed using a ``GaussianPrior``. The
        ``gaussian_limits`` set the physical lower and upper limits of this ``GaussianPrior``, such that parameter samples
        can not go beyond these limits.

notation
--------

The notation configs define the labels of every model-component parameter and its derived quantities, which are
used when visualizing results (for example labeling the axis of the PDF triangle plots output by a non-linear search).
Two examples using the 1D ``data`` fitting example for the config file **label.ini** are:

[label]
    centre -> str
        The label given to that parameter for non-linear search plots using that parameter, e.g. the PDF plots. For
        example, if centre=x, the plot axis will be labeled 'x'.

[superscript]
    Gaussian -> str
        The subscript used on certain plots that show the results of different model-components. For example, if
        Gaussian=g, plots where the Gaussian are plotted will have a subscript g.

The **label_format.ini** config file specifies the format certain parameters are output as in output files like the
*model.results* file... _roadmap:

Road Map
========

**PyAutoFit** is in active development and the road-map of features currently planned in the short and long term are
listed and described below:

**Non-Linear Searches:**

We are always striving to add new non-linear searches to **PyAutoFit*. In the short term, we aim to provide a wrapper to the many method available in the ``scipy.optimize`` library with support for outputting results to hard-disk.

If you would like to see a non-linear search implemented in **PyAutoFit** please `raise an issue on GitHub <https://github.com/rhayes777/PyAutoFit/issues>`_!

**Graphical Models**

Graphical models allow one to compose complex models that fit for global trends in many model-fits to individual
datasets. This feature is in development and described in the **Graphical Models** tab of the readthedocs,
however it is still in beta.

**Approximate Bayesian Computation**

Approximate Bayesian Computational (ABC) allows for one to infer parameter values for likelihood functions that are
intractable, by simulating many datasets and extracting from them a summary statistic that is compared to the
observed dataset.

ABC in **PyAutoFit** will be closely tied to the Database tools, ensuring that the simulation, fitting and extraction
of summary statistics can be efficiently scaled up to extremely large datasets... _software:

Software
--------

The following software projects use **PyAutoFit**:

`PyAutoLens <https://github.com/Jammy2211/PyAutoLens>`_ -
Astronomy software for modeling Strong Gravitational Lenses.

`PyAutoGalaxy <https://github.com/Jammy2211/PyAutoGalaxy>`_ -
Astronomy software for modeling galaxy light profiles and dynamics.

`PyAutoCTI <https://github.com/Jammy2211/PyAutoCTI>`_ -
Software for modeling Charge Transfer Inefficiency induced by radiation damage to CCDs... _adding_a_model_component:

Adding a Model Component
========================

Adding a class
--------------

The ``autofit_workspace`` comes ready for fitting 1D ``Gaussian`` and ``Exponential`` profiles, complete with configuration
files, analysis classes and example scripts.

However, once you're familiar with **PyAutoFit**, you will want to add your own model-components, specific to your
model-fitting task. There are a couple of extra steps that come with doing this, associated with configuration files,
that this brief guide explains.

The model-component we are going to add will perform a ``y = mx + c`` linear fit to noisy data drawn from a straight
line. We're only going to focus on the steps necessary to add this new model component, so we'll omit writing an
``Analysis`` class and performing the actual fit itself.

To perform a linear fit, we require a ``LinearFit`` model-component that fits the data with a
line ``y = mx + c`` or equivalently ``y = (gradient * x) + intercept``.

.. code-block:: bash

    class LinearFit:

        def __init__(self, gradient=1.0, intercept=0.0):

            self.gradient = gradient
            self.intercept = intercept

        def profile_1d_via_xvalues_from(self, xvalues):

            return (self.gradient * xvalues) + self.intercept

As should be clear on by now, the class ``LinearFit`` defines our model-component which has free parameters  ``gradient``
and ``intercept``.

However, if we tried to make this a ``Model`` PyAutoFit would raises an error, e.g.

.. code-block:: bash

    model = af.Model(LinearFit)

The error will read something like ``'KeyError: No prior config found for class LinearFit and path gradient in directories C:\\Users\\Jammy\\Code\\PyAuto\\autofit_workspace\\config\\priors'``.

**PyAutoFit** is informing us that it cannot find prior configuration files for the ``LinearFit`` model-component and that 
they are therefore missing from the folder ``autofit_workspace/config/priors``.


Every model-component must have a ``.json`` config file in the ``autofit_workspace/config/priors`` folder, so 
that **PyAutoFit** knows the default priors to associate with the model-component. If we do not manually override 
priors, these are the priors that will be used by default when a model-fit is performed.

Next, inspect the `TemplateObject.json  <https://github.com/Jammy2211/autofit_workspace/blob/master/config/priors/TemplateObject.json>`_ configuration file in ``autofit_workspace/config/priors``. You should see
the following ``.json`` text:

.. code-block:: bash

    {
        "parameter0": {
            "type": "Uniform",
            "lower_limit": 0.0,
            "upper_limit": 1.0
        },
        "parameter1": {
            "type": "Gaussian",
            "mean": 0.0,
            "sigma": 0.1,
            "lower_limit": "-inf",
            "upper_limit": "inf"
        }
    }

This specifies the default priors on two parameters, named ``parameter0`` and ``parameter1``. The ``type`` is the type of 
prior assumed by **PyAutoFit** by default for its corresponding parameter. 

In the example above: 

- ``parameter0`` is given a ``UniformPrior`` with limits between 0.0 and 1.0. 
- ``parameter1`` a ``GaussianPrior`` with mean 0.0 and sigma 1.0.

The ``lower_limit`` and ``upper_limit`` of a ``GaussianPrior`` define the boundaries of what parameter values are 
physically allowed. If a model-component is given a value outside these limits during model-fitting the model is
instantly resampled and discarded.
 
We can easily adapt this template for our ``LinearFit`` model component. First, copy and paste the `TemplateObject.json  <https://github.com/Jammy2211/autofit_workspace/blob/master/config/priors/TemplateObject.json>`_
file to create a new file called ``LinearFit.json``. 

**PyAutoFit** matches the name of the class to the name of the configuration file, therefore it is a requirement that 
the configuration file is named ``LinearFit.json``.

Next, rename ``parameter0`` to ``gradient``, ``parameter1`` to ``intercept`` and make it so both assume a ``UniformPrior`` 
between -10.0 to 10.0.

The ``.json`` file should read as follows:

.. code-block:: bash

    {
        "gradient": {
            "type": "Uniform",
            "lower_limit": -10.0,
            "upper_limit": 10.0
        },
        "intercept": {
            "type": "Uniform",
            "lower_limit": -10.0,
            "upper_limit": 10.0
        }
    }

We should now be able to make a ``Model`` of the ``LinearFit`` class.

.. code-block:: bash

    model = af.Model(LinearFit)

Adding a Module
---------------

For larger projects, it is not ideal to have to write all the model-component classes in a single Python script, 
especially as we may have many different model components. We instead would prefer them to be in their own dedicated 
Python module.

open the file:

- ``autofit_workspace/scripts/overview/adding_a_model_component/linear_fit.py``  OR
- ``autofit_workspace/notebooks/overview/adding_a_model_component/linear_fit.py``

Here, you will see the ``LinearFit`` class above is contained in the module ``linear_fit.py``. There is also a ``PowerFit`` 
class, fits the function ``y = m (x**p) + c``.

If we import this module and try to make a  ``Model`` of the ``linear_fit.LinearFit`` or ``linear_fit.PowerFit``
classes, we receive the same configuration error as before.

.. code-block:: bash

    import linear_fit
    
    model = af.Model(linear_fit.LinearFit)
    model = af.Model(linear_fit.PowerFit)

This is because if a model-component is contained in a Python module, the prior configuration file must be named after
that ``module`` and structured to contain Python class itself.

Open the file ``autofit_workspace/config/priors/template_module.json``, (https://github.com/Jammy2211/autofit_workspace/blob/master/config/priors/template_module.json) which reads as follows:

.. code-block:: bash
    
    {
        "ModelComponent0": {
            "parameter0": {
                "type": "Uniform",
                "lower_limit": 0.0,
                "upper_limit": 1.0
            },
            "parameter1": {
                "type": "LogUniform",
                "lower_limit": 1e-06,
                "upper_limit": 1e6
            },
            "parameter2": {
                "type": "Uniform",
                "lower_limit": 0.0,
                "upper_limit": 25.0
            }
        },
        "ModelComponent1": {
            "parameter0": {
                "type": "Uniform",
                "lower_limit": 0.0,
                "upper_limit": 1.0
            },
            "parameter1": {
                "type": "LogUniform",
                "lower_limit": 1e-06,
                "upper_limit": 1e6
            },
            "parameter2": {
                "type": "Uniform",
                "lower_limit": 0.0,
                "upper_limit": 1.0
            }
        }
    }

This looks very similar to ``TemplateObject``, the only differences are:

 - It now contains the model-component class name in the configuration file, e.g. ``ModelComponent0``, ``ModelComponent1``.
 - It includes multiple model-components, whereas ``TemplateObject.json`` corresponded to only one model component.
 
We can again easily adapt this template for our ``linear_fit.py`` module. Copy, paste and rename the ``.json`` file to
``linear_fit.json`` (noting again that **PyAutoFit** matches the module name to the configuration file) and update the
parameters as follows:

.. code-block:: bash
    
    {
        "LinearFit": {
            "gradient": {
                "type": "Uniform",
                "lower_limit": -10.0,
                "upper_limit": 10.0
            },
            "intercept": {
                "type": "Uniform",
                "lower_limit": -10.0,
                "upper_limit": 10.0
            }
        },
        "PowerFit": {
            "gradient": {
                "type": "Uniform",
                "lower_limit": -10.0,
                "upper_limit": 10.0
            },
            "intercept": {
                "type": "Uniform",
                "lower_limit": -10.0,
                "upper_limit": 10.0
            },
            "power": {
                "type": "Uniform",
                "lower_limit": 0.0,
                "upper_limit": 10.0
            }
        }
    }

We are now able to create both the ``linear_fit.LinearFit`` and ``linear_fit.PowerFit`` objects as ``Model``'s.

.. code-block:: bash

    model = af.Model(linear_fit.LinearFit)
    model = af.Model(linear_fit.PowerFit)

Optional Configs
----------------

There are a couple more configuration files you can optionally update, which change how results are output. Open the 
following configuration files:

``autofit_workspace/config/notation/label.ini``
``autofit_workspace/config/notation/label_format.ini``

These configuration files include the following additional settings for our model components:

``label_ini`` -> [label]: 
   This is a short-hand label for each parameter of each model-component used by certain **PyAutoFit** output files.

``label_ini`` -> [superscript]:
   A subscript for the model-component used by certain **PyAutoFit** output files.

``label_format.ini`` -> [format]
   The format that the values of a parameter appear in the ``model.results`` file.

For our ``LinearFit`` update the ``label.ini`` config as follows:

.. code-block:: bash

    [label]
    centre=x
    normalization=I
    sigma=sigma
    rate=\lambda
    gradient=m
    intercept=c
    power=p

.. code-block:: bash

    [superscript]
    Gaussian=g
    Exponential=e
    LinearFit=lin
    PowerFit=pow

and ``label_format.ini`` as:

.. code-block:: bash

    [format]
    centre={:.2f}
    normalization={:.2f}
    sigma={:.2f}
    rate={:.2f}
    gradient={:.2f}
    intercept={:.2f}
    power={:.2f}

You should now be able to add your own model-components to your **PyAutoFit** project!.. _troubleshooting:

Troubleshooting
===============

Current Working Directory
-------------------------

**PyAutoFit** scripts assume that the ``autofit_workspace`` directory is the Python working directory. This means
that, when you run an example script, you should run it from the ``autofit_workspace`` as follows:

.. code-block:: bash

    cd path/to/autofit_workspace (if you are not already in the autofit_workspace).
    python3 scripts/overview/simple/fit.py

The reasons for this are so that **PyAutoFit** can:

 - Load configuration settings from config files in the ``autofit_workspace/config`` folder.
 - Load example data from the ``autofit_workspace/dataset`` folder.
 - Output the results of models fits to your hard-disk to the ``autofit/output`` folder.
 - Import modules from the ``autofit_workspace``, for example ``from autofit_workspace.transdimensional import pipelines``.

If you have any errors relating to importing modules, loading data or outputting results it is likely because you
are not running the script with the ``autofit_workspace`` as the working directory!

Support
-------

If you are still having issues with installation or using **PyAutoFit** in general, please raise an issue on the
`autofit_workspace issues page <https://github.com/Jammy2211/autofit_workspace/issues>`_ with a description of the
problem and your system setup (operating system, Python version, etc.)... _source:

Building From Source
====================

Building from source means that you clone (or fork) the **PyAutoFit** GitHub repository and run **PyAutoFit** from
there. Unlike ``conda`` and ``pip`` this provides a build of the source code that you can edit and change, to
contribute the development **PyAutoFit** or experiment with yourself!

First, clone (or fork) the **PyAutoFit** GitHub repository:

.. code-block:: bash

    git clone https://github.com/Jammy2211/PyAutoFit

Next, install the **PyAutoFit** dependencies via pip:

.. code-block:: bash

   pip install -r PyAutoFit/requirements.txt

If you are using a `conda` environment, add the source repository as follows:

.. code-block:: bash

   conda-develop PyAutoFit

Alternatively, if you are using a Python environment include the **PyAutoFit** source repository in your PYTHONPATH
(noting that you must replace the text ``/path/to`` with the path to the **PyAutoFit** directory on your computer):

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:/path/to/PyAutoFit

Finally, check the **PyAutoFit** unit tests run and pass (you may need to install pytest via ``pip install pytest``):

.. code-block:: bash

   cd /path/to/PyAutoFit
   python3 -m pytest.. _pip:

Installation with pip
=====================

We strongly recommend that you install **PyAutoFit** in a
`Python virtual environment <https://www.geeksforgeeks.org/python-virtual-environment/>`_, with the link attached
describing what a virtual enviroment is and how to create one.

The latest version of **PyAutoFit** is installed via pip as follows (specifying the version as shown below ensures
the installation has clean dependencies):

.. code-block:: bash

    pip install autofit==2021.10.14.1

If this raises no errors **PyAutoFit** is installed! If there is an error check out
the `troubleshooting section <https://pyautofit.readthedocs.io/en/latest/installation/troubleshooting.html>`_.

Next, clone the ``autofit workspace`` (the line ``--depth 1`` clones only the most recent branch on
the ``autofit_workspace``, reducing the download size):

.. code-block:: bash

   cd /path/on/your/computer/you/want/to/put/the/autofit_workspace
   git clone https://github.com/Jammy2211/autofit_workspace --depth 1
   cd autofit_workspace

Run the ``welcome.py`` script to get started!

.. code-block:: bash

   python3 welcome.py.. _overview:

Overview
========

**PyAutoFit** requires Python 3.6+ and support the Linux, MacOS and Windows operating systems.

**PyAutoFit** can be installed via the Python distribution `Anaconda <https://www.anaconda.com/>`_ or using
`Pypi <https://pypi.org/>`_ to ``pip install`` **PyAutoFit** into your Python distribution.

We recommend Anaconda as it manages the installation of many major libraries used by **PyAutoFit** (e.g. numpy, scipy,
matplotlib, etc.) making installation more straight forward.

The installation guide for both approaches can be found at:

- `Anaconda installation guide <https://pyautofit.readthedocs.io/en/latest/installation/conda.html>`_

- `PyPI installation guide <https://pyautofit.readthedocs.io/en/latest/installation/pip.html>`_

Users who wish to build **PyAutoFit** from source (e.g. via a `git clone`) should follow
our `building from source installation guide <https://pyautofit.readthedocs.io/en/latest/installation/source.html>`_.

Known Issues
------------

There are currently no known issues with installing **PyAutoFit**.

Dependencies
------------

**PyAutoConf** https://github.com/rhayes777/PyAutoConf

**dynesty** https://github.com/joshspeagle/dynesty

**emcee** https://github.com/dfm/emcee

**PySwarms** https://github.com/ljvmiranda921/pyswarms

**astropy** https://www.astropy.org/

**corner.py** https://github.com/dfm/corner.py

**matplotlib** https://matplotlib.org/

**numpy** https://numpy.org/

**scipy** https://www.scipy.org/

And the following optional dependencies:

**PyMultiNest** http://johannesbuchner.github.io/pymultinest-tutorial/install.html.. _conda:

Installation with conda
=======================

Installation via a conda environment circumvents compatibility issues when installing certain libraries. This guide
assumes you have a working installation of `conda <https://conda.io/miniconda.html>`_.

First, create a conda environment (we name is ``autofit`` to signify it is for the **PyAutoFit** install).

The command below creates this environment with some of the bigger package requirements, the rest will be installed
with **PyAutoFit** via pip:

.. code-block:: bash

    conda create -n autofit numpy scipy

Activate the conda environment (you will have to do this every time you want to run **PyAutoFit**):

.. code-block:: bash

    conda activate autofit

The latest version of **PyAutoFit** is installed via pip as follows (specifying the version as shown below ensures
the installation has clean dependencies):

.. code-block:: bash

    pip install autofit==2021.10.14.1

Next, clone the ``autofit workspace`` (the line ``--depth 1`` clones only the most recent branch on
the ``autofit_workspace``, reducing the download size):

.. code-block:: bash

   cd /path/on/your/computer/you/want/to/put/the/autofit_workspace
   git clone https://github.com/Jammy2211/autofit_workspace --depth 1
   cd autofit_workspace

Run the `welcome.py` script to get started!

.. code-block:: bash

   python3 welcome.py.. _chapter_graphical_models:

Chapter: Graphical Models
=========================

NOTE: This is an in development feature and this chapter is incomplete.

In this chapter, we take you through how to compose and fit graphical models in **PyAutoFit**. Graphical models
can fit large datasets with a model that has 'local' parameters specific to each individual dataset and 'global'
parameters that fit for higher level parameters.

You can start the tutorials right now by going to `our binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD>`_
and navigating to the folder `notebooks/howtofit/chapter_graphical_models`. They are also on the `autofit_workspace <https://github.com/Jammy2211/autofit_workspace>`_.

The chapter contains the following tutorials:

`Tutorial 1: Global Model <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_graphical_models/tutorial_1_global_model.ipynb>`_
- An example of inferring global parameters from a dataset by fitting the model to each dataset one-by-one.

`Tutorial 2: Graphical Model <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_graphical_models/tutorial_2_graphical_model.ipynb>`_
- Fitting the dataset with a graphical model that fits all datasets simultaneously to infer the global parameters.. _howtofit:

HowToFit Lectures
=================

To learn how to use **PyAutoFit**, the best starting point is the **HowToFit** lecture series, which are found on
the `autofit_workspace <https://github.com/Jammy2211/autofit_workspace>`_ and at
our `binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD>`_.

The lectures are provided as *Jupyter notebooks* and currently consist of 3 chapters:

**Introduction**: How to perform model-fitting with **PyAutoFit** and analyse the results using the ``Aggregator``.

**Graphical Models**: How to compose and fit graphical models to large datasets (notebooks + feature are in development).

Config File Path
----------------

If, when running the first notebook, you get an error related to config files, this most likely means that
**PyAutoFit** is unable to find the config files in your autofit workspace. Checkout the
`configs section <https://pyautofit.readthedocs.io/en/latest/general/configs.html>`_ for a description of how to
fix this.

Jupyter Notebooks
-----------------

The tutorials are supplied as *Jupyter notebooks*, which come with a ``.ipynb`` suffix. For those new to
Python, *Jupyter notebooks* are a different way to write, view and use Python code. Compared to the
traditional Python scripts, they allow:

- Small blocks of code to be viewed and run at a time
- Images and visualization from a code to be displayed directly underneath it.
- Text script to appear between the blocks of code.

This makes them an ideal way for us to present the HowToFit lecture series, therefore I recommend you get
yourself a Jupyter notebook viewer (https://jupyter.org/) if you havent done so already.

If you *really* want to use Python scripts, all tutorials are supplied a ``.py`` python files in the ``scripts``
folder of each chapter.

For actual **PyAutoFit** use I recommend you use Python scripts. Therefore, as you go through the lecture
series you will notice that we will transition you to Python scripts.

Code Style and Formatting
-------------------------

You may notice the style and formatting of our Python code looks different to what you are used to. For
example, it is common for brackets to be placed on their own line at the end of function calls, the inputs
of a function or class may be listed over many separate lines and the code in general takes up a lot more
space then you are used to.

This is intentional, because we believe it makes the cleanest, most readable code possible. In fact, lots
of people do, which is why we use an auto-formatter to produce the code in a standardized format. If you're
interested in the style and would like to adapt it to your own code, check out the Python auto-code formatter
``black``.

https://github.com/python/black.. _chapter_1_introduction:

Chapter 1: Introduction
=======================

In chapter 1, we introduce you to the **PyAutoFit** and describe how to set up your own model, fit it to data via
a non-linear search and inspect and interpret the results.

You can start the tutorials right now by going to `our binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD>`_
and navigating to the folder ``notebooks/howtofit/chapter_1_introduction``. They are also on the `autofit_workspace <https://github.com/Jammy2211/autofit_workspace>`_.

The chapter contains the following tutorials:

`Tutorial 1: Model Composition <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_1_model_composition.ipynb>`_
- Composing a model in **PyAutoFit**.

`Tutorial 2: Fitting Data <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_2_fitting_data.ipynb>`_
- Fitting a model with an input set of parameters to data.

`Tutorial 3: Parameter Space and Priors <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_3_parameter_space_and_priors.ipynb>`_
- The Concepts of a parameter space and priors.

`Tutorial 4: Non Linear Search <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_4_non_linear_search.ipynb>`_
- Finding the model parameters that best-fit the data.

`Tutorial 5: Complex Models <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_5_complex_models.ipynb>`_
- Composing and fitting complex models.

`Tutorial 6: Results and Samples <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_1_introduction/tutorial_6_results_and_samples.ipynb>`_
- The results of a model-fit output by **PyAutoFit**.
.. _chapter_database:

Chapter: Database
=================

In this chapter 1, we introduce you to the **PyAutoFit**'s database functionality that allows one to output all model-fitting results
to a relational database. This means model-fits to large datasets can then be loaded and inspected in a Jupyter notebook.

You can start the tutorials right now by going to `our binder <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/HEAD>`_
and navigating to the folder ``notebooks/howtofit/chapter_database``. They are also on the `autofit_workspace <https://github.com/Jammy2211/autofit_workspace>`_.

The chapter contains the following tutorials:

`Tutorial 1: Fitting Multiple Datasets <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_database/tutorial_1_fitting_multiple_datasets.ipynb>`_
- Fitting a model to multiple similar datasets.

`Tutorial 2: Aggregator <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_database/tutorial_2_aggregator.ipynb>`_
- Loading large libraries of results using in-built database tools.

`Tutorial 3: Querying <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_database/tutorial_3_querying.ipynb>`_
- Filtering the loaded results to find specific results of interest.

`Tutorial 4: Data and Models <https://mybinder.org/v2/gh/Jammy2211/autofit_workspace/release?filepath=notebooks/howtofit/chapter_database/tutorial_4_data_and_models.ipynb>`_
- Replotting the model, data and fits of a set of results.
.. _database:

Database
========

The default behaviour of **PyAutoFit** is for model-fitting results to be output to hard-disc in folders, which are
straight forward to navigate and manually check the model-fitting results. For small model-fitting tasks this is
sufficient, however many users have a need to perform many model fits to very large datasets, making the manual
inspection of results time consuming.

PyAutoFit's database feature outputs all model-fitting results as a sqlite3 (https://docs.python.org/3/library/sqlite3.html)
relational database, such that all results can be efficiently loaded into a Jupyter notebook or Python script for
inspection, analysis and interpretation. This database supports advanced querying, so that specific
model-fits (e.g., which fit a certain model or dataset) can be loaded.

Session
-------

To make it so that results are output to an .sqlite database we simply open a database session and pass this session
to the non-linear search:

.. code-block:: bash

    session = af.db.open_database("database.sqlite")

    emcee = af.Emcee(
        path_prefix=path.join("features", "database"),
        session=session,  # This instructs the search to write to the .sqlite database.
    )

Unique Tag
----------

When a model-fit is performed, a unique identifier is generated based on the model and non-linear search. However,
if we were to fit many different datasets with the same model and non-linear search, they would all use the same
unique identifier and not be distinguishable by the database.

We can overcome this by using the name of the dataset as the ``unique_tag`` passed to the search, which is used alongside
the model and search to create the unique identifier:

.. code-block:: bash

    session = af.db.open_database("database.sqlite")

    dataset_name = "example_dataset_0"

    emcee = af.Emcee(
        path_prefix=path.join("features", "database"),
        unique_tag=dataset_name,  # This makes the unique identifier use the dataset name
        session=session,  # This instructs the search to write to the .sqlite database.
    )

Loading
-------

Lets suppose that we have performed 100 model-fits to 100 1D Gaussians, and when we ran **PyAutoFit** we told it
to write to the ``.sqlite`` database file. We can load these results in a Python script or Jupyter notebook using
the ``Aggregator``:

.. code-block:: bash

    agg = Aggregator.from_database("path/to/output/database.sqlite")

We can now use the ``Aggregator`` to inspect the results of all model-fits. For example, we can load the ``Samples``
object of all 100 model-fits, which contains information on the best-fit model, posterior, Bayesian evidence, etc.

Below, we use the samples generator to create a list of the maximum log likelihood of every model-fit and print it:

.. code-block:: bash

    for samples in agg.values("samples"):

        print(max(samples.log_likelihood))

This object (and all objects loaded by the ``Aggregator``) are returned as a generator (as opposed to a list,
dictionary or other Python type). This is because generators do not store large arrays or classes in memory until they
are used, ensuring that when we are manipulating large sets of results we do not run out of memory!

We can iterate over the samples to print the maximum log likelihood model ``Gaussian`` of every fit:

.. code-block:: bash

    for samps in agg.values("samples"):

        instance = samps.max_log_likelihood_instance

        print("Maximum Likelihood Model-fit \n")
        print("centre = ", instance.centre)
        print("normalization = ", instance.normalization)
        print("sigma = ", instance.sigma)


Queries
-------

The ``Aggregator`` contains tools for querying the database for certain results, for example to load subsets of
model-fits. This can be done in many different ways, depending on what information you want.

Below, we query based on the model fitted. For example, we can load all results which fitted a ``Gaussian``
model-component, which in this simple example is all 100 model-fits (note that when we performed the model fit, we
composed model using the name ``gaussian``):

.. code-block:: bash

    gaussian = agg.model.gaussian
    agg_query = agg.query(gaussian == m.Gaussian)

Queries using the results of model-fitting are also supported. Below, we query the database to find all fits where the
inferred value of ``sigma`` for the ``Gaussian`` is less than 3.0:

.. code-block:: bash

    agg_query = agg.query(gaussian.sigma < 3.0)

Advanced queries can be constructed using logic, for example we below we combine the two queries above to find all
results which fitted a ``Gaussian`` AND (using the & symbol) inferred a value of sigma less than 3.0.

The OR logical clause is also supported via the symbol |.

.. code-block:: bash

    agg_query = agg.query((gaussian == m.Gaussian) & (gaussian.sigma < 3.0))

We can query using the ``unique_tag`` to load the model-fit to a specific dataset:

.. code-block:: bash

    agg_query = agg.query(agg.unique_tag == "example_dataset_0")

Info
----

An ``info`` dictionary can be passed into a model-fit, which contains information on the model-fit. The example below
creates an ``info`` dictionary which is passed to the model-fit, which is then loaded via the database.

.. code-block:: bash

    info = {"example_key": "example_value"}

    emcee.fit(model=model, analysis=analysis, info=info)

    agg = Aggregator.from_database("path/to/output/database.sqlite")

    info_gen = agg.values("info")

Databases are an extremely powerful feature for users tasked with fitting extremely large datasets as well as fitting
many different models, where the scale of the problem can make the management of the large quantity of results produced
prohibitive. This is especially true on high performance computing facilities, which often have restrictions on the
number of files that a user can store on the machine.

Wrap Up
-------

If you'd like to see the ``Aggregator`` in action, checkout the
`database example <https://github.com/Jammy2211/autofit_workspace/blob/master/notebooks/features/database.ipynb>`_ on the
``autofit_workspace``.

The Database Chapter of the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_
provides more details, including how to visualize the results of a model fit fully... _sensitivity_mapping:

Sensitivity Mapping
===================

Bayesian model comparison allows us to take a dataset, fit it with multiple models and use the Bayesian evidence to
quantify which model objectively gives the best-fit following the principles of Occam's Razor.

However, a complex model may not be favoured by model comparison not because it is the 'wrong' model, but simply
because the dataset being fitted is not of a sufficient quality for the more complex model to be favoured. Sensitivity
mapping allows us to address what quality of data would be needed for the more complex model to be favoured or
alternatively for what sets of model parameter values it would be favoured for data of a given quality.

In order to do this, sensitivity mapping involves us writing a function that uses the model(s) to simulate a dataset.
We then use this function to simulate many datasets, for many different models, and fit each dataset using the same
model-fitting procedure we used to perform Bayesian model comparison. This allows us to infer how much of a Bayesian
evidence increase we should expect for datasets of varying quality and / or models with different parameters.

Data
----

To illustrate sensitivity mapping we will again use the example of fitting 1D Gaussian's in noisy data. This 1D data
includes a small feature to the right of the central ``Gaussian``, a second ``Gaussian`` centred on pixel 70.


.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_with_feature.png
  :width: 600
  :alt: Alternative text

Model Comparison
----------------

Before performing sensitivity mapping, we will quickly perform Bayesian model comparison on this data to get a sense
for whether the ``Gaussian`` feature is detectable and how much the Bayesian evidence increases when it is included in
the model.

We therefore fit the data using two models, one where the model is a single ``Gaussian``.

.. code-block:: bash

    model = af.Collection(gaussian_main=m.Gaussian)

    dynesty = af.DynestyStatic(
        path_prefix=path.join("features", "sensitivity_mapping", "single_gaussian"),
        nlive=100,
        iterations_per_update=500,
    )

    result_single = dynesty.fit(model=model, analysis=analysis)

For the second model it contains two ``Gaussians``. To avoid slow model-fitting and more clearly pronounce the results of
model comparison, we restrict the centre of the ``gaussian_feature`` to its true centre of 70 and sigma value of 0.5.

.. code-block:: bash

    model = af.Collection(gaussian_main=m.Gaussian, gaussian_feature=m.Gaussian)
    model.gaussian_feature.centre = 70.0
    model.gaussian_feature.sigma = 0.5

    dynesty = af.DynestyStatic(
        path_prefix=path.join("features", "sensitivity_mapping", "two_gaussians"),
        nlive=100,
        iterations_per_update=500,
    )

    result_multiple = dynesty.fit(model=model, analysis=analysis)

We can now print the ``log_evidence`` of each fit and confirm the model with two ``Gaussians`` was preferred to the model
with just one ``Gaussian``.

.. code-block:: bash

    print(result_single.samples.log_evidence)
    print(result_multiple.samples.log_evidence)

On my laptop, the increase in Bayesian evidence for the more compelx model is ~30, which is significant.

The model comparison above shows that in this dataset, the ``Gaussian`` feature was detectable and that it increased the
Bayesian evidence by ~25. Furthermore, the normalization of this ``Gaussian`` was ~0.3.

A lower value of normalization makes the ``Gaussian`` fainter and harder to detect. We will demonstrate sensitivity mapping
by answering the following question, at what value of normalization does the ``Gaussian`` feature become undetectable and
not provide us with a noticeable increase in Bayesian evidence?

Base Model
----------

To begin, we define the ``base_model`` that we use to perform sensitivity mapping. This model is used to simulate every
dataset. It is also fitted to every simulated dataset without the extra model component below, to give us the Bayesian
evidence of the every simpler model to compare to the more complex model.

The ``base_model`` corresponds to the ``gaussian_main`` above.

.. code-block:: bash

    base_model = af.Collection(gaussian_main=m.Gaussian)

Perturbation Model
------------------

We now define the ``perturbation_model``, which is the model component whose parameters we iterate over to perform
sensitivity mapping. Many instances of the ``perturbation_model`` are created and used to simulate the many datasets
that we fit. However, it is only included in half of the model-fits corresponding to the more complex models whose
Bayesian evidence we compare to the simpler model-fits consisting of just the ``base_model``.

The ``perturbation_model`` is therefore another ``Gaussian`` but now corresponds to the ``gaussian_feature`` above.

By fitting both of these models to every simulated dataset, we will therefore infer the Bayesian evidence of every
model to every dataset. Sensitivity mapping therefore maps out for what values of ``normalization`` in the ``gaussian_feature``
does the more complex model-fit provide higher values of Bayesian evidence than the simpler model-fit. We also fix the
values ot the ``centre`` and ``sigma`` of the ``Gaussian`` so we only map over its ``normalization``.

.. code-block:: bash

    perturbation_model = af.Model(m.Gaussian)
    perturbation_model.centre = 70.0
    perturbation_model.sigma = 0.5
    perturbation_model.normalization = af.UniformPrior(lower_limit=0.01, upper_limit=100.0)

Simulation
----------

We are performing sensitivity mapping to determine how bright the ``gaussian_feature`` needs to be in order to be
detectable. However, every simulated dataset must include the ``main_gaussian``, as its presence in the data will effect
the detectability of the ``gaussian_feature``.

We can pass the ``main_gaussian`` into the sensitivity mapping as the ``simulation_instance``, meaning that it will be used
in the simulation of every dataset. For this example we use the inferred ``main_gaussian`` from one of the model-fits
performed above.

.. code-block:: bash

    simulation_instance = result_single.instance

We now write the ``simulate_function``, which takes the ``instance`` of our model (defined above) and uses it to
simulate a dataset which is subsequently fitted.

Note that when this dataset is simulated, the quantity ``instance.perturbation`` is used in the ``simulate_function``.
This is an instance of the ``gaussian_feature``, and it is different every time the ``simulate_function`` is called.

In this example, this ``instance.perturbation`` corresponds to different ``gaussian_feature``'s with values of
``normalization`` ranging over 0.01 -> 100.0, such that our simulated datasets correspond to a very faint and very bright
gaussian features.

.. code-block:: bash

    def simulate_function(instance):

        """
        Specify the number of pixels used to create the xvalues on which the 1D line of the profile is generated using and
        thus defining the number of data-points in our data.
        """
        pixels = 100
        xvalues = np.arange(pixels)

        """
        Evaluate the ``Gaussian`` and Exponential model instances at every xvalues to create their model profile and sum
        them together to create the overall model profile.

        This print statement will show that, when you run ``Sensitivity`` below the values of the perturbation use fixed
        values of ``centre=70`` and ``sigma=0.5``, whereas the normalization varies over the ``step_size`` based on its prior.
        """

        print(instance.perturbation.centre)
        print(instance.perturbation.normalization)
        print(instance.perturbation.sigma)

        model_line = instance.gaussian_main.profile_1d_via_xvalues_from(xvalues=xvalues) + instance.perturbation.profile_1d_via_xvalues_from(xvalues=xvalues)

        """Determine the noise (at a specified signal to noise level) in every pixel of our model profile."""
        signal_to_noise_ratio = 25.0
        noise = np.random.normal(0.0, 1.0 / signal_to_noise_ratio, pixels)

        """
        Add this noise to the model line to create the line data that is fitted, using the signal-to-noise ratio to compute
        noise-map of our data which is required when evaluating the chi-squared value of the likelihood.
        """
        data = model_line + noise
        noise_map = (1.0 / signal_to_noise_ratio) * np.ones(pixels)

        return Dataset(data=data, noise_map=noise_map)

Here are what the two most extreme simulated datasets look like, corresponding to the highest and lowest normalization values

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/sensitivity_data_low.png
  :width: 600
  :alt: Alternative text

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/sensitivity_data_high.png
  :width: 600
  :alt: Alternative text

Summary
-------

We can now combine all of the objects created above and perform sensitivity mapping. The inputs to the ``Sensitivity``
object below are:

- ``simulation_instance``: This is an instance of the model used to simulate every dataset that is fitted. In this example it contains an instance of the ``gaussian_main`` model component.

- ``base_model``: This is the simpler model that is fitted to every simulated dataset, which in this example is composed of a single ``Gaussian`` called the ``gaussian_main``.

- ``perturbation_model``: This is the extra model component that alongside the ``base_model`` is fitted to every simulated dataset, which in this example  is composed of two ``Gaussians`` called the ``gaussian_main`` and ``gaussian_feature``.

- ``simulate_function``: This is the function that uses the ``instance`` and many instances of the ``perturbation_model`` to simulate many datasets that are fitted with the ``base_model`` and ``base_model`` + ``perturbation_model``.

- ``step_size``: The size of steps over which the parameters in the ``perturbation_model`` are iterated. In this example, normalization has a ``LogUniformPrior`` with lower limit 1e-4 and upper limit 1e2, therefore the ``step_size`` of 0.5 will simulate and fit just 2 datasets where the normalization is 1e-4 and 1e2.

- ``number_of_cores``: The number of cores over which the sensitivity mapping is performed, enabling parallel processing.

(Note that for brevity we have omitted a couple of extra inputs in this example, which can be found by going to the
full example script on the ``autofit_workspace``).

.. code-block:: bash

    sensitivity = s.Sensitivity(
        search=search,
        simulation_instance=simulation_instance,
        base_model=base_model,
        perturbation_model=perturbation_model,
        simulate_function=simulate_function,
        analysis_class=Analysis,
        step_size=0.5,
        number_of_cores=2,
    )

    sensitivity_result = sensitivity.run()

Here are what the fits to the two most extreme simulated datasets look like, for the models including the Gaussian
feature.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/sensitivity_data_low_fit.png
  :width: 600
  :alt: Alternative text

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/sensitivity_data_high_fit.png
  :width: 600
  :alt: Alternative text

The key point to note is that for every dataset, we now have a model-fit with and without the model ``perturbation``. By
compairing the Bayesian evidence of every pair of fits for every value of ``normalization`` we are able to determine when
our model was sensitivity to the ``Gaussian`` feature and therefore could detect it!.. _empirical_bayes:

Empirical Bayes
===============

In Empirical Bayes, the priors that are applied to the analysis are learnt from the data itself. This can manifest in
different ways, for example in MCMC analysis by fitting the data to determine where one initializes the MCMC walkers,
or in nested sampling fitting the data to derive priors where initial sampling occurs.

**PyAutoFit** has mature support for fitting models following an empirical Bayes approach, which typically first perform
efficient initial fits to data to estimate the priors on the model parameters. These priors are then used to fit more
complex models which require more computationally fitting procedures.

Search Chaining
---------------

To perform a model-fit, we typically compose one model and fit it to our data using one non-linear search.

Search chaining fits many different models to a dataset using a chained sequence of non-linear searches. Initial
fits are performed using simplified model parameterizations and faster non-linear fitting techniques. The results of
these simplified fits can then be used to initialize fits using a higher dimensionality model with more detailed
non-linear search.

To fit highly complex models our aim is therefore to **granularize** the fitting procedure into a series
of **bite-sized** searches which are faster and more reliable than fitting the more complex model straight away.

Our ability to construct chained non-linear searches that perform model fitting more accurately and efficiently relies
on our **domain specific knowledge** of the model fitting task. For example, we may know that our dataset contains
multiple features that can be accurately fitted separately before performing a joint fit, or that certain parameters
share minimal covariance such that allowing us to fix them to certain values in earlier model-fits.

We may also know tricks that can speed up the fitting of the initial model, for example reducing the size of the data
or changing making a likelihood evaluation faster (most likely at the expense of the quality of the fit itself). By
using chained searches speed-ups can be relaxed towards the end of the model-fitting sequence when we want the most
precise, most accurate model that best fits the dataset available.

Data
----

In this example we demonstrate search chaining using the example data where there are two ``Gaussians`` that are visibly
split:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x2_split.png
  :width: 600
  :alt: Alternative text

Approach
--------

Instead of fitting them simultaneously using a single non-linear search consisting of N=6 parameters, we break
this model-fit into a chained of three searches where:

1) The first model fits just the left ``Gaussian`` where N=3.
2) The second model fits just the right ``Gaussian`` where again N=3.
3) The final model is fitted with both ``Gaussians`` where N=6. Crucially, the results of the first two searches are used to initialize the search and tell it the highest likelihood regions of parameter space.

By initially fitting parameter spaces of reduced complexity we can achieve a more efficient and reliable model-fitting
procedure.

Search 1
--------

To fit the left ``Gaussian``, our first ``analysis`` receive only half data removing the right ``Gaussian``. Note that
this give a speed-up in log likelihood evaluation.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x2_left.png
  :width: 600
  :alt: Alternative text

We now create a search to fit this data. Given the simplicity of the model, we can use a low number of live points
to achieve a fast model-fit (had we fitted the more complex model right away we could not of done this).

.. code-block:: bash

    model = af.Collection(gaussian_left=m.Gaussian)

    dynesty = af.DynestyStatic(
        name=("search[1]__left_gaussian"),
        nlive=30,
    )

    search_2_result = dynesty.fit(model=model, analysis=analysis)

By plotting the result we can see we have fitted the left ``Gaussian`` reasonably well.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x2_left_fit.png
  :width: 600
  :alt: Alternative text

Search 2
--------

We now repeat the above process for the right ``Gaussian``.

We could remove the data on the left like we did the ``Gaussian`` above. However, we are instead going to fit the full
dataset. To fit the left Gaussian we use the maximum log likelihood model of the model inferred in search 1.

For search chaining, **PyAutoFit** has many convenient methods for passing the results of a search to a subsequence
search. Below, we achieve this by passing the result of the search above as an ``instance``.

.. code-block:: bash

    model = af.Collection(
        gaussian_left=search_1_result.instance.gaussian_left,
        gaussian_right=m.Gaussian
    )

We now run our second Dynesty search to fit the right ``Gaussian``. We can again exploit the simplicity of the model
and use a low number of live points to achieve a fast model-fit.

.. code-block:: bash

    dynesty = af.DynestyStatic(
        name=("search[2]__right_gaussian"),
        path_prefix=path.join("features", "search_chaining"),
        nlive=30,
        iterations_per_update=500,
    )

    search_2_result = dynesty.fit(model=model, analysis=analysis)

We can now see our model has successfully fitted both Gaussian's:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x2_right_fit.png
  :width: 600
  :alt: Alternative text

Search 3
--------

We now fit both ``Gaussians``'s simultaneously, using the results of the previous two searches to initialize where
the non-linear searches parameter space.

To pass the result in this way we use the command ``result.model``, which in contrast to ``result.instance`` above passes
the parameters not as the maximum log likelihood values but as ``GaussianPrior``'s that are fitted for by the
non-linear search.

The ``mean`` and ``sigma`` value of each parmeter's ``GaussianPrior`` are set using the results of searches 1 and
2 to ensure our model-fit only searches the high likelihood regions of parameter space.

.. code-block:: bash

    model = af.Collection(
        gaussian_left=search_1_result.model.gaussian_left,
        gaussian_right=search_2_result.model.gaussian_right
    )

    dynesty = af.DynestyStatic(
        name=("search[3]__both_gaussians"),
        path_prefix=path.join("features", "search_chaining"),
        nlive=100,
        iterations_per_update=500,
    )

    search_3_result = dynesty.fit(model=model, analysis=analysis)

We can now see our model has successfully fitted both Gaussians simultaneously:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x2_fit.png
  :width: 600
  :alt: Alternative text

Wrap Up
-------

This fit used a technique called 'prior passing' to pass results from searches 1 and 2 to search 3. Full details of how
prior passing works can be found in the ``search_chaining.ipynb`` feature notebook... _search_grid_search:

Search Grid-Search
==================

A classic method to perform model-fitting is a grid search, where the parameters of a model are divided on to a grid of
values and the likelihood of each set of parameters on this grid is sampled. For low dimensionality problems this
simple approach can be sufficient to locate high likelihood solutions, however it scales poorly to higher dimensional
problems.

**PyAutoFit** can perform a search grid search, which allows one to perform a grid-search over a subset of parameters
within a model, but use a non-linear search to fit for the other parameters. The parameters over which the grid-search
is performed are also included in the model fit and their values are simply confined to the boundaries of their grid
cell by setting these as the lower and upper limits of a ``UniformPrior``.

The benefits of using a search grid search are:

- For problems with complex and multi-model parameters spaces it can be difficult to robustly and efficiently perform model-fitting. If specific parameters are known to drive the multi-modality then sampling over a grid can ensure the parameter space of each individual model-fit is not multi-modal and therefore sampled more accurately and efficiently.

- It can provide a goodness-of-fit measure (e.g. the Bayesian evidence) of many model-fits over the grid. This can provide additional insight into where the model does and does not fit the data well, in a way that a standard non-linear search does not.

- The search grid search is embarrassingly parallel, and if sufficient computing facilities are available one can perform model-fitting faster in real-time than a single non-linear search. The **PyAutoFit** search grid search includes an option for parallel model-fitting via the Python ``multiprocessing`` module.

Data
----

In this example we will demonstrate the search grid search feature, again using the example of fitting 1D Gaussian's
in noisy data. This 1D data includes a small feature to the right of the central ``Gaussian``, a second ``Gaussian``
centred on pixel 70.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_with_feature.png
  :width: 600
  :alt: Alternative text

Basic Fit
---------

Without the search grid search we can fit this data as normal, by composing and fitting a model
containing two ``Gaussians``'s.

.. code-block:: bash

    model = af.Collection(gaussian_main=m.Gaussian, gaussian_feature=m.Gaussian)

    analysis = a.Analysis(data=data, noise_map=noise_map)

    dynesty = af.DynestyStatic(
        name="single_fit",
        nlive=100,
    )

    result = dynesty.fit(model=model, analysis=analysis)

For test runs on my laptop it is 'hit or miss' whether the feature is fitted correctly. This is because although models
including the feature corresponds to the highest likelihood solutions, they occupy a small volume in parameter space
which the non linear search may miss.

The image below shows a fit where we failed to detect the feature:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_with_feature_fit_no_feature.png
  :width: 600
  :alt: Alternative text

Grid Search
-----------

Lets now perform the search grid search using the ``SearchGridSearch`` object:

.. code-block:: bash

    dynesty = af.DynestyStatic(
        name="grid_fit",
        nlive=100,
    )

    grid_search = af.SearchGridSearch(
        search=dynesty,
        number_of_steps=5,
        number_of_cores=1,
    )


We specified two new inputs to the ``SearchGridSearch``:

``number_of_steps``: The number of steps in the grid search that are performed which is set to 5 below. Because the
prior on the parameter ``centre`` is a ``UniformPrior`` from 0.0 -> 100.0, this means the first grid search will
set the prior on the centre to be a ``UniformPrior`` from 0.0 -> 20.0. The second will run from 20.0 -> 40.0, the
third 40.0 -> 60.0, and so on.

``number_of_cores``: The number of cores the grid search will parallelize the run over. If ``number_of_cores=1``, the
search is run in serial. For > 1 core, 1 core is reserved as a farmer, e.g., if ``number_of_cores=4`` then up to 3
searches will be run in parallel.

We can now run the grid search, where we specify the parameter over which the grid search is performed, in this case
the ``centre`` of the ``gaussian_feature`` in our model.

.. code-block:: bash

    grid_search_result = grid_search.fit(
        model=model,
        analysis=analysis,
        grid_priors=[model.gaussian_feature.centre]
    )

Result
------

This returns a ``GridSearchResult``, which includes information on every model-fit performed on the grid. For example,
I can use it to print the ``log_evidence`` of all 5 model-fits.

.. code-block:: bash

    print(grid_search_result.log_evidence_values)

This shows a peak evidence value on the 4th cell of grid-search, where the ``UniformPrior`` on the ``centre`` ran from
60 -> 80 and therefore included the Gaussian feature. By plotting this model-fit we can see it has successfully
detected the feature.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_with_feature_fit_feature.png
  :width: 600
  :alt: Alternative text

A multi-dimensional grid search can be easily performed by adding more parameters to the ``grid_priors`` input.

The fit below belows performs a 5x5 grid search over the ``centres`` of both ``Gaussians``

.. code-block:: bash

    grid_search_result = grid_search.fit(
        model=model,
        analysis=analysis,
        grid_priors=[model.gaussian_feature.centre, model.gaussian_main.centre]
    ).. _result:

Results & Samples
=================

A non-linear search's fit function returns a ``Result`` object:

.. code-block:: bash

   analysis = Analysis(data=data, noise_map=noise_map)

   emcee = af.Emcee(number_of_cores=4)

   result = emcee.fit(model=model, analysis=analysis)

Here, we'll look in detail at what information is contained in the ``Result``.

Samples
-------

A result contains a ``Samples`` object, which contains information on the non-linear sampling, for example the parameters.
The parameters are stored as a list of lists, where the first entry corresponds to the sample index and second entry
the parameter index.

.. code-block:: bash

    samples = result.samples

    print("Final 10 Parameters:")
    print(samples.parameter_lists[-10:])

    print("Sample 10`s third parameter value (Gaussian -> sigma)")
    print(samples.parameter_lists[9][2], "\n")

The Samples class also contains the log likelihood, log prior, log posterior and weight_list of every accepted sample,
where:

- The log likelihood is the value evaluated from the likelihood function (e.g. -0.5 * chi_squared + the noise normalized).

- The log prior encodes information on how the priors on the parameters maps the log likelihood value to the log posterior value.

- The log posterior is log_likelihood + log_prior.

- The weight gives information on how samples should be combined to estimate the posterior. The weight values depend on the sampler used, for MCMC samples they are all 1 (e.g. all weighted equally).

Lets inspect the last 10 values of each for the analysis.

.. code-block:: bash

    print("Final 10 Log Likelihoods:")
    print(samples.log_likelihood_list[-10:])

    print("Final 10 Log Priors:")
    print(samples.log_prior_list[-10:])

    print("Final 10 Log Posteriors:")
    print(samples.log_posterior_list[-10:])

    print("Final 10 Sample Weights:")
    print(samples.weight_list[-10:], "\n")

Posterior
---------

The ``Result`` object therefore contains the full posterior information of our non-linear search, that can be used for
parameter estimation. The median pdf vector is readily available from the ``Samples`` object, which estimates the every
parameter via 1D marginalization of their PDFs.

.. code-block:: bash

    median_pdf_vector = samples.median_pdf_vector

The samples include methods for computing the error estimates of all parameters via 1D marginalization at an input sigma
confidence limit. This can be returned as the size of each parameter error:

.. code-block:: bash

    error_vector_at_upper_sigma = samples.error_vector_at_upper_sigma(sigma=3.0)
    error_vector_at_lower_sigma = samples.error_vector_at_lower_sigma(sigma=3.0)

    print("Upper Error values (at 3.0 sigma confidence):")
    print(error_vector_at_upper_sigma)

    print("lower Error values (at 3.0 sigma confidence):")
    print(error_vector_at_lower_sigma, "\n")

They can also be returned at the values of the parameters at their error values:

.. code-block:: bash

    vector_at_upper_sigma = samples.vector_at_upper_sigma(sigma=3.0)
    vector_at_lower_sigma = samples.vector_at_lower_sigma(sigma=3.0)

    print("Upper Parameter values w/ error (at 3.0 sigma confidence):")
    print(vector_at_upper_sigma)
    print("lower Parameter values w/ errors (at 3.0 sigma confidence):")
    print(vector_at_lower_sigma, "\n")

**PyAutoFit** includes many visualization tools for plotting the results of a non-linear search, for example we can
make a corner plot of the probability density function (PDF):

.. code-block:: bash

    emcee_plotter = aplt.EmceePlotter(samples=result.samples)
    emcee_plotter.corner()

Here is an example of how a PDF estimated for a lens model appears:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/images/corner.png
  :width: 600
  :alt: Alternative text

Other Vectors
-------------

The samples contain many useful vectors, including the samples with the highest likelihood and posterior values:

.. code-block:: bash

    max_log_likelihood_vector = samples.max_log_likelihood_vector
    max_log_posterior_vector = samples.max_log_posterior_vector

    print("Maximum Log Likelihood Vector:")
    print(max_log_likelihood_vector)

    print("Maximum Log Posterior Vector:")
    print(max_log_posterior_vector, "\n")


Labels
------

These vectors return the results as a list, which means you need to know the parameter ordering. The
list of ``parameter_names`` are available as a property of the ``Samples``, as are ``parameter_labels``
which can be used for labeling figures:

.. code-block:: bash

    samples.model.parameter_names
    samples.model.parameter_labels

Instances
---------

``Result``'s can instead be returned as an ``instance``, which is an instance of the model using the Python
classes used to compose it:

.. code-block:: bash

    max_log_likelihood_instance = samples.max_log_likelihood_instance

    print("Max Log Likelihood Gaussian Instance:")
    print("Centre = ", max_log_likelihood_instance.centre)
    print("normalization = ", max_log_likelihood_instance.normalization)
    print("Sigma = ", max_log_likelihood_instance.sigma)


For our example problem of fitting a 1D ``Gaussian`` profile, this makes it straight forward to plot
the maximum likelihood model:

.. code-block:: bash

    model_data = samples.max_log_likelihood_instance.profile_1d_via_xvalues_from(
        xvalues=np.arange(data.shape[0])
    )

    plt.plot(range(data.shape[0]), data)
    plt.plot(range(data.shape[0]), model_data)
    plt.title("Illustrative toy model fit to 1D Gaussian line profile data.")
    plt.xlabel("x values of line profile")
    plt.ylabel("Line profile normalization")
    plt.show()
    plt.close()

All methods above are available as an ``instance``:

.. code-block:: bash

    median_pdf_instance = samples.median_pdf_instance
    instance_at_upper_sigma = samples.instance_at_upper_sigma
    instance_at_lower_sigma = samples.instance_at_lower_sigma
    error_instance_at_upper_sigma = samples.error_instance_at_upper_sigma
    error_instance_at_lower_sigma = samples.error_instance_at_lower_sigma

An ``instance`` of any accepted sample can be created:

.. code-block:: bash

    instance = samples.instance_from_sample_index(sample_index=500)

Bayesian Evidence
-----------------

If a nested sampling non-linear search is used, the Bayesian evidence of the model is also
available which enables model comparison to be performed:

.. code-block:: bash

    log_evidence = samples.log_evidence

Result Extensions
-----------------

You might be wondering what else the results contains, as nearly everything we discussed above was a part of its
``samples`` property! The answer is, not much, however the result can be extended to include  model-specific results for
your project.

We detail how to do this in the **HowToFit** lectures, but for the example of fitting a 1D Gaussian we could extend
the result to include the maximum log likelihood profile:

.. code-block:: bash

    max_log_likelihood_profile = samples.max_log_likelihood_profile

Database
--------

For large-scaling model-fitting problems to large datasets, the results of the many model-fits performed can be output
and stored in a queryable sqlite3 database. The ``Result`` and ``Samples`` objects have been designed to streamline the
analysis and interpretation of model-fits to large datasets using the database.

The database is described `here <https://pyautofit.readthedocs.io/en/latest/features/database.html>`_

Wrap-Up
-------

More information on the ``Result`` class can be found at the
`results examples <https://github.com/Jammy2211/autofit_workspace/blob/master/notebooks/overview/simple/result.ipynb>`_ on
the ``autofit_workspace``. More details are provided in tutorial 7 or chapter 1 of
the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_.. _non_linear_search:

Non-linear Search
=================

**PyAutoFit** currently supports three types of non-linear search algorithms:

- **Optimizers**: ``PySwarms``.
- **MCMC**: ``emcee`` and ``Zeus``.
- **Nested Samplers**: ``dynesty``, ``UltraNest`` (optional) and ``PyMultiNest`` (optional).

Functionality
-------------

**PyAutoFit** extends the functionality of each non-linear search to ensure that they always perform the
following tasks, even if the original package does not:

- Stores the results of the non-linear search to the hard-disk, writing the results to human-readable files.

- Allows the non-linear search to be resumed if a previous run was finished.

- Can write results and associated metadata to an sqlite database for querying and inspection post model-fit.

- Extends the functionality of the non-linear search's, for example providing auto-correlation analysis and
  stopping criteria for MCMC algorithms.

Settings
--------

We've seen that we can call a non-linear search as follows:

.. code-block:: bash

   analysis = Analysis(data=data, noise_map=noise_map)

   emcee = af.Emcee(name="example_mcmc")

   result = emcee.fit(model=model, analysis=analysis)

However, ``Emcee`` has many settings associated with it (the number of walkers, the number of steps they take,
etc.). Above, we did not pass them to the ``Emcee`` constructor and they use the default values found in the
``autofit_workspace`` configuration files ``autofit_workspace/config/non_linear/mcmc/Emcee.ini``, which can be
viewed at this `link <https://github.com/Jammy2211/autofit_workspace/blob/master/config/non_linear/mcmc/Emcee.ini>`_.

Of course, we can manually specify all of the parameters instead:

.. code-block:: bash

   analysis = Analysis(data=data, noise_map=noise_map)

   emcee = af.Emcee(
       name="example_mcmc",
       nwalkers=50,
       nsteps=2000,
       initializer=af.InitializerBall(lower_limit=0.49, upper_limit=0.51),
       auto_correlations_settings=af.AutoCorrelationsSettings(
           check_for_convergence=True,
           check_size=100,
           required_length=50,
           change_threshold=0.01,
       ),
   )

   result = emcee.fit(model=model, analysis=analysis)

A number of these parameters are not part of the ``emcee`` package, but additional functionality added by
**PyAutoFit**:

- Initialization methods for the walkers are provided, including the strategy recommended at this `page <https://emcee.readthedocs.io/en/stable/user/faq/?highlight=ball#how-should-i-initialize-the-walkers>`_ where the walkers are initialized as a compact 'ball' in parameter space.

- Auto correlation lengths can be checked during sampling and used to determine whether the MCMC chains have converged, terminating ``emcee`` before all ``nwalkers`` have taken all ``nsteps``, as discussed at this `link <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/>`_.

The nested sampling algorithm ``dynesty`` has its own config file for default settings, which are at
this `link <https://github.com/Jammy2211/autofit_workspace/blob/master/config/non_linear/nest/Dynesty.ini>`_.
``DynestyStatic`` parameters can be manually specified as follows:

.. code-block:: bash

   analysis = Analysis(data=data, noise_map=noise_map)

   dynesty = af.DynestyStatic(
       name="example_nest",
       nlive=150,
       bound="multi",
       sample="auto",
       bootstrap=None,
       enlarge=None,
       update_interval=None,
       vol_dec=0.5,
       vol_check=2.0,
       walks=25,
       facc=0.5,
       slices=5,
       fmove=0.9,
       max_move=100,
   )

   result = dynesty.fit(model=model, analysis=analysis)

Output Paths
------------

We can also customize the output folder and path structure where results are output. The output folder is set
using the **PyAutoFit** parent project **PyAutoConf** and the following command:

.. code-block:: bash

   from autoconf import conf

   conf.instance.push(new_path="path/to/config", output_path="path/to/output")

The path structure within this folder of a given non-linear search is set using the ``path_prefix``.

Results are output to a folder which is a collection of random characters, which is the 'unique_identifier' of
the model-fit. This identifier is generated based on the model fitted and search used, such that an identical
combination of model and search generates the same identifier.

This ensures that rerunning an identical fit will use the existing results to resume the model-fit. In contrast, if
you change the model or search, a new unique identifier will be generated, ensuring that the model-fit results are
output into a separate folder.

The example code below would output the results to the
path ``/path/to/output/folder_0/folder_1/unique_tag/example_mcmc/sihfiuy838h``:

.. code-block:: bash

    from os import path

   emcee = af.Emcee(
       path_prefix=path.join("folder_0", "folder_1"),
       name="example_mcmc"
   )

Parallelization
---------------

Most searches support parallel analysis using the Python ``multiprocessing`` module. This distributes the
non-linear search analysis over multiple CPU's, speeding up the run-time roughly by the number of CPUs used.

The in-built parallelization of Libraries such as ``emcee`` and ``dynesty`` can be slow, because the default behaviour
is for them to pass the full likelihood function to every CPU. If this function includes a large dataset that is being
fitted, this can lead to long communication overheads and slow performance.

**PyAutoFit** implements *sneaky parallelization*, whereby the data is passed to every CPU before the model-fit. This
requires no extra user input and is performed by default. To perform a parallel search, you simply specify
the ``number_of_cores`` parameter (which is also found in the default config files):

.. code-block:: bash

   analysis = Analysis(data=data, noise_map=noise_map)

   emcee = af.Emcee(number_of_cores=4)

   result = emcee.fit(model=model, analysis=analysis)

Wrap-Up
-------

We are always looking to add more non-linear searches to **PyAutoFit**. If you are the developer of a package check out
our `contributions section <https://github.com/rhayes777/PyAutoFit/blob/master/CONTRIBUTING.md>`_ and please
contact us!.. _model_complex:

Model Composition
=================

Lets extend our example of fitting a 1D ``Gaussian`` profile to a problem where the data contains a signal from
two 1D profiles. Specifically, it contains signals from a 1D ``Gaussian`` signal and 1D symmetric ``Exponential``.

Data
----

The example ``data`` with errors (black), including the model-fit we'll perform (red) and individual
``Gaussian`` (blue dashed) and ``Exponential`` (orange dashed) components are shown below:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/images/toy_model_fit_x2.png
  :width: 600
  :alt: Alternative text

Model Composition
-----------------

we again define our 1D ``Gaussian`` profile as a *model component* in **PyAutoFit**:

.. code-block:: bash

    class Gaussian:
        def __init__(
            self,
            centre=0.0,     # <- PyAutoFit recognises these
            normalization=0.1,  # <- constructor arguments are
            sigma=0.01,     # <- the Gaussian's parameters.
        ):

            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

        def profile_1d_via_xvalues_from(self, xvalues):

            transformed_xvalues = xvalues - self.centre

            return np.multiply(
                np.divide(self.normalization, self.sigma * np.sqrt(2.0 * np.pi)),
                np.exp(-0.5 * np.square(np.divide(transformed_xvalues, self.sigma))),
            )

Now lets define a new *model component*, a 1D ``Exponential``, using the same Python class format:

.. code-block:: bash

    class Exponential:
        def __init__(
            self,
            centre=0.0,     # <- PyAutoFit recognises these
            normalization=0.1,  # <- constructor arguments are
            rate=0.01,      # <- the Exponential's parameters.
        ):

            self.centre = centre
            self.normalization = normalization
            self.rate = rate

        def profile_1d_via_xvalues_from(self, xvalues):

            transformed_xvalues = xvalues - self.centre

            return self.normalization * np.multiply(
                self.rate, np.exp(-1.0 * self.rate * abs(transformed_xvalues))
            )

We are now fitting multiple *model components*, therefore we create each component using the ``Model`` object we
used in the previous tutorial and put them together in a ``Collection`` to build the overall model.

.. code-block:: bash

    gaussian = af.Model(Gaussian)
    exponential = af.Model(Exponential)

    model = af.Collection(gaussian=gaussian, exponential=exponential)

The ``Collection`` allows us to *compose* models using multiple classes. This model is defined with 6 free
parameters (3 for the ``Gaussian``, 3 for the ``Exponential``), thus the dimensionality of non-linear parameter space
is 6.

Analysis
--------

The *model components* given to the ``Collection`` were also given names, in this case, ``gaussian`` and
``exponential``.

You are free to choose whichever names you want;  the names are used to pass the ``instance`` to the ``Analysis`` class:

.. code-block:: bash

    class Analysis(af.Analysis):

        def __init__(self, data, noise_map):

            super().__init__()

            self.data = data
            self.noise_map = noise_map

        def log_likelihood_function(self, instance):

            """
            The 'instance' that comes into this method is a Collection. It contains
            instances of every class we instantiated it with, where each instance is named
            following the names given to the Collection, which in this example is a
            Gaussian (with name 'gaussian) and Exponential (with name 'exponential'):
            """

            print("Gaussian Instance:")
            print("Centre = ", instance.gaussian.centre)
            print("normalization = ", instance.gaussian.normalization)
            print("Sigma = ", instance.gaussian.sigma)

            print("Exponential Instance:")
            print("Centre = ", instance.exponential.centre)
            print("normalization = ", instance.exponential.normalization)
            print("Rate = ", instance.exponential.rate)

            """
            Get the range of x-values the data is defined on, to evaluate the model of the
            line profiles.
            """

            xvalues = np.arange(self.data.shape[0])

            """
            The instance variable is a list of our model components. We can iterate over
            this list, calling their profile_1d_via_xvalues_from and summing the result to compute
            the summed line profile of our model.
            """

            model_data = sum([line.profile_1d_via_xvalues_from(xvalues=xvalues) for line in instance])

            """
            Fit the model line profile data to the observed data, computing the residuals and
            chi-squared.
            """

            residual_map = self.data - model_data
            chi_squared_map = (residual_map / self.noise_map) ** 2.0
            log_likelihood = -0.5 * sum(chi_squared_map)

            return log_likelihood

Model Fit
---------

Performing the *model-fit* uses the same steps as the previous example, whereby we  *compose* our *model* (now using a
``Collection``), instantiate the ``Analysis`` and pass them a non-linear search. In this example, we'll use
the nested sampling algorithm ``dynesty``, using the ``DynestyStatic`` sampler.

.. code-block:: bash

    model = af.Collection(gaussian=Gaussian, exponential=Exponential)

    analysis = Analysis(data=data, noise_map=noise_map)

    dynesty = af.DynestyStatic(name="example_search")

    result = dynesty.fit(model=model, analysis=analysis)

Model Priors
------------

Now, lets consider how we *customize* the models that we *compose*. To begin, lets *compose* a model using a single
``Gaussian`` with the ``Model`` object:

.. code-block:: bash

    gaussian = af.Model(Gaussian)

By default, the priors on the ``Gaussian``'s parameters are loaded from configuration files. If you have downloaded the
``autofit_workspace`` you can find these files at the path ``autofit_workspace/config/priors``. Alternatively,
you can check them out at this `link <https://github.com/Jammy2211/autofit_workspace/tree/master/config>`_.

Priors can be manually specified as follows:

.. code-block:: bash

    gaussian.centre = af.UniformPrior(lower_limit=0.0, upper_limit=100.0)
    gaussian.normalization = af.LogUniformPrior(lower_limit=0.0, upper_limit=1e2)
    gaussian.sigma = af.GaussianPrior(mean=10.0, sigma=5.0, lower_limit=0.0, upper_limit=np.inf)

These priors will be used by the non-linear search to determine how it samples parameter space. The ``lower_limit``
and ``upper_limit`` on the ``GaussianPrior`` set the physical limits of values of the parameter, specifying that the
``sigma`` value of the ``Gaussian`` cannot be negative.

We can fit this model, with all new priors, using a non-linear search as we did before:

.. code-block:: bash

    analysis = Analysis(data=data, noise_map=noise_map)

    emcee = af.Emcee(name="another_example_search")

    # The model passed here now has updated priors!

    result = emcee.fit(model=gaussian, analysis=analysis)

We can *compose* and *customize* the priors of multiple model components as follows:

.. code-block:: bash

    gaussian = af.Model(Gaussian)
    gaussian.normalization = af.UniformPrior(lower_limit=0.0, upper_limit=1e2)

    exponential = af.Model(Exponential)
    exponential.centre = af.UniformPrior(lower_limit=0.0, upper_limit=100.0)
    exponential.normalization = af.UniformPrior(lower_limit=0.0, upper_limit=1e2)
    exponential.rate = af.UniformPrior(lower_limit=0.0, upper_limit=10.0)

    model = af.Collection(gaussian=gaussian, exponential=exponential)

Model Customization
-------------------

The model can be *customized* to fix any *parameter* of the model to an input value:

.. code-block:: bash

    gaussian.sigma = 0.5

This fixes the ``Gaussian``'s ``sigma`` value to 0.5, reducing the number of free parameters and therefore
dimensionality of *non-linear parameter space* by 1.

We can also link two parameters, such that they always share the same value:

.. code-block:: bash

    model.gaussian.centre = model.exponential.centre

In this model, the ``Gaussian`` and ``Exponential`` will always be centrally aligned. Again, this reduces
the number of free *parameters* by 1.

Finally, assertions can be made on parameters that remove values that do not meet those assertions
from *non-linear parameter space*:

.. code-block:: bash

    gaussian.add_assertion(gaussian.sigma > 5.0)
    gaussian.add_assertion(gaussian.normalization > exponential.normalization)

Here, the ``Gaussian``'s ``sigma`` value must always be greater than 5.0 and its ``normalization`` is greater
than that of the ``Exponential``.

Wrap Up
-------

If you'd like to perform the fit shown in this script, checkout the
`complex examples <https://github.com/Jammy2211/autofit_workspace/tree/master/notebooks/overview/complex>`_ on the
``autofit_workspace``. We provide more details **PyAutoFit** works in the tutorials 5 and 6 of
the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_.. _composition:

Multi-level Models
==================

A **multi-level model** is a hierarchy of model components, where the different levels express the conditional
dependence between different parameters and model-components. Using hierarchies of Python classes **PyAutoFit** can
construct **multi-level models** via the ``Model`` and ``Collection`` objects, and these can be linked together to form
one over-arching model.

Astronomy Use Case
------------------

Multi-level models are a fairly abstract concept, and so to describe them we are going to introduce a real-world
model-fitting example. We will use an example from Astronomy; fitting images of gravitationally lensed galaxies.
This is the science case that sparked the development of **PyAutoFit** as a spin off of our astronomy software
`PyAutoLens <https://github.com/Jammy2211/PyAutoLens>`_.

The schematic below depicts a strong gravitational lens:

.. image:: https://raw.githubusercontent.com/Jammy2211/PyAutoLens/master/docs/overview/images/lensing/schematic.jpg
  :width: 600
  :alt: Alternative text

**Credit: F. Courbin, S. G. Djorgovski, G. Meylan, et al., Caltech / EPFL / WMKO**
https://www.astro.caltech.edu/~george/qsolens/

A strong gravitational lens is a system consisting of multiple galaxy's down the light-of-sight to earth. To model
a strong lens, we ray-trace the traversal of light throughout the Universe so as to fit it to imaging data of a strong
lens. The amount light is deflected by is defined by the distances between each galaxy, which is called their redshift.

Model Overview
--------------

We therefore need a model which contains separate model-components for every galaxy, and where each galaxy contains
separate model-components describing its light and mass. A multi-level representation of this model is as follows:

.. image:: https://github.com/rhayes777/PyAutoFit/blob/master/docs/overview/image/lens_model.png?raw=true
  :width: 600
  :alt: Alternative text

The image above shows that we need a model consisting of individual model-components for:

 1) The lens galaxy's *light* and *mass*.
 2) The source galaxy's *light*.

We also need each galaxy to be a **model-component** itself and for each of them to contain an additional parameter,
its ``redshift``. The galaxies can then be combined into an overall model for the strong lens system.

Model Example
-------------

To model the light of a galaxy, we define a ``LightProfile`` as a Python class, which behaves in the same way as
the ``Gaussian`` used in other **PyAutoFit** tutorials:

.. code-block:: bash

    class LightProfile:

        def __init__(
            self,
            centre: typing.Tuple[float, float] = (0.0, 0.0),
            normalization: float = 0.1,
            radius: float = 0.6,
        ):
            """
            A light profile used in Astronomy to represent the surface brightness distribution of galaxies.

            Parameters
            ----------
            centre
                The (y,x) coordinates of the profile centre.
            normalization
                Overall normalization normalisation of the light profile.
            radius
                The circular radius containing half the light of this profile.
            """

            self.centre = centre
            self.normalization = normalization
            self.effective_radius = effective_radius

        def image_from_grid(self, grid: np.ndarray) -> np.ndarray:
            """This function creates an image of the light profile, which is used in strong lens model-fitting"""
            ...

We have omitted the code that creates the image from the light profile as we want to focus purely on multi-level model
composition with **PyAutoFit**.

We also define a ``MassProfile``:

.. code-block:: bash

    class MassProfile:
        def __init__(
            self,
            centre: typing.Tuple[float, float] = (0.0, 0.0),
            mass: float = 1.0,
        ):
            """
            A mass profile used in Astronomy to represent the mass distribution of galaxies.

            Parameters
            ----------
            centre
                The (y,x) coordinates of the profile centre.
            mass
                The mass normalization of the profile.
            """

            self.centre = centre
            self.mass = mass

        def deflections_from_grid(self, grid: np.ndarray) -> np.ndarray:
            """This function describes the deflection of light due to the mass, which is used in strong lens model-fitting"""
            ...

We have again omitted the code which computes how this mass profile deflects the path of light.

We now define a ``Galaxy`` object, which contains instances of light and mass profiles and its redshift (e.g. distance 
from Earth):

.. code-block:: bash

    class Galaxy:

        def __init__(
            self,
            redshift: float,
            light_profile_list: Optional[List] = None,
            mass_profile_list: Optional[List] = None,
        ):
            """
            A galaxy, which contains light and mass profiles at a specified redshift.

            Parameters
            ----------
            redshift
                The redshift of the galaxy.
            light_profile_list
                A list of the galaxy's light profiles.
            mass_profile_list
                A list of the galaxy's mass profiles.
            """

            self.redshift = redshift
            self.light_profile_list = light_profile_list
            self.mass_profile_list = mass_profile_list

        def image_from_grid(self, grid: np.ndarray) -> np.ndarray:
            """Returns the image of all light profiles."""
            ...

        def deflections_from_grid(self, grid: np.ndarray) -> np.ndarray:
            """Returns the deflection angles of all mass profiles."""
            ...

If we were not composing a model, the code below shows how one would create an instance of the foreground lens galaxy,
which in the image above contains a light and mass profile:

.. code-block:: bash

    light = LightProfile(centre=(0.0, 0.0), normalization=10.0, radius=2.0)
    mass = MassProfile(centre=(0.0, 0.0), mass=0.5)

    lens = Galaxy(redshift=0.5, light_profile_list=[light], mass_profile_list=[mass])

The code creates instances of the ``LightProfile`` and ``MassProfile`` classes and uses them to create an
instance of the ``Galaxy`` class. This uses a **hierarchy of Python classes**.

Multi-level Model
-----------------

We can compose a multi-level model using this same hierarchy of classes, using the ``Model`` and ``Collection`` objects.

Lets first create a model of the lens galaxy:

.. code-block:: bash

    light = af.Model(LightProfile)
    mass = af.Model(MassProfile)

    lens = af.Model(
        cls=Galaxy,
        redshift=0.5,
        light_profile_list=[light],
        mass_profile_list=[mass]
    )

Lets consider what the code above is doing:

1) We use a ``Model`` to create the overall model component. The ``cls`` input is the ``Galaxy`` class, therefore the overall model that is created is a ``Galaxy``.

2) **PyAutoFit** next inspects whether the key word argument inputs to the ``Model`` match any of the ``__init__`` constructor arguments of the ``Galaxy`` class. This determine if these inputs are to be composed as **model sub-components** of the overall ``Galaxy`` model.

3) **PyAutoFit** matches the ``light_profile_list`` and  ``mass_profile_list`` inputs, noting they are passed as separate lists containing ``Model``'s of the ``LightProfile`` and ``MassProfile`` classes. They are both created as sub-components of the overall ``Galaxy`` model.

4) It also matches the ``redshift`` input, making it a fixed value of 0.5 for the model and not treating it as a free parameter.

We can confirm this by printing the ``prior_count`` of the lens, and noting it is 7 (4 parameters for
the ``LightProfile`` and 3 for the ``MassProfile``).

.. code-block:: bash

    print(lens.prior_count)
    print(lens.light_profile_list[0].prior_count)
    print(lens.mass_profile_list[0].prior_count)

The ``lens`` behaves exactly like the model-components we are used to previously. For example, we can unpack its
individual parameters to customize the model, where below we:

 1) Align the light profile centre and mass profile centre.
 2) Customize the prior on the light profile ``one``.
 3) Fix the ``one`` of the mass profile to 0.8.

.. code-block:: bash

    lens.light_profile_list[0].centre = lens.mass_profile_list[0].centre
    lens.light_profile_list[0].one = af.UniformPrior(lower_limit=0.7, upper_limit=0.9)
    lens.mass_profile_list[0].one = 0.8

We can now create a model of our source galaxy using the same API.

.. code-block:: bash

    source = af.Model(
        astro.Galaxy,
        redshift=1.0,
        light_profile_list=[af.Model(astro.lp.LightProfile)]
    )

We can now create our overall strong lens model, using a ``Collection`` in the same way we have seen previously. 

.. code-block:: bash

    model = af.Collection(galaxies=af.Collection(lens=lens, source=source))

The model contains both galaxies in the strong lens, alongside all of their light and mass profiles.

For every iteration of the non-linear search **PyAutoFit** generates an instance of this model, where all of the
``LightProfile``, ``MassProfile`` and ``Galaxy`` parameters of the are determined via their priors.

An example instance is show below:

.. code-block:: bash

    instance = model.instance_from_prior_medians()

    print("Strong Lens Model Instance:")
    print("Lens Galaxy = ", instance.galaxies.lens)
    print("Lens Galaxy Light = ", instance.galaxies.lens.light_profile_list)
    print("Lens Galaxy Light Centre = ", instance.galaxies.lens.light_profile_list[0].centre)
    print("Lens Galaxy Mass Centre = ", instance.galaxies.lens.mass_profile_list[0].centre)
    print("Source Galaxy = ", instance.galaxies.source)

This model can therefore be used in a **PyAutoFit** ``Analysis`` class and ``log_likelihood_function``.

Extensibility
-------------

This example highlights how multi-level models can make certain model-fitting problem fully extensible. For example:

 1) A ``Galaxy`` class can be created using any combination of light and mass profiles. Although this was not shown 
explicitly in this example, this is because it implements their ``image_from_grid`` and ``deflections_from_grid`` methods
as the sum of individual profiles.

 2) The overall strong lens model can contain any number of ``Galaxy``'s, as these methods and their redshifts are used
to implement the lensing calculations in the ``Analysis`` class and ``log_likelihood_function``.

Thus, for problems of this nature, we can design and write code in a way that fully utilizes **PyAutoFit**'s multi-level
modeling features to compose and fits models of arbitrary complexity and dimensionality.

To illustrate this further, consider the following dataset which is called a **strong lens galaxy cluster**:

.. image:: https://github.com/rhayes777/PyAutoFit/blob/master/docs/overview/image/cluster_example.jpg?raw=true
   :width: 600
   :alt: Alternative text

For this strong lens, there are many tens of strong lens galaxies as well as multiple background source galaxies.
However, despite it being a significantly more complex system than the single-galaxy strong lens we modeled above,
our use of graphical models ensures that we can model such datasets without any additional code development, for
example:

.. code-block:: bash

    lens_0 = af.Model(
        Galaxy,
        redshift=0.5,
        light_profile_list=[af.Model(LightProfile)],
        mass_profile_list=[af.Model(MassProfile)]
    )

    lens_1 = af.Model(
        Galaxy,
        redshift=0.5,
        light_profile_list=[af.Model(LightProfile)],
        mass_profile_list=[af.Model(MassProfile)]
    )

    source_0 = af.Model(
        astro.Galaxy,
        redshift=1.0,
        light_profile_list=[af.Model(LightProfile)]
    )

    # ... repeat for desired model complexity ...

    model = af.Collection(
        galaxies=af.Collection(
            lens_0=lens_0,
            lens_1=lens_1,
            source_0=source_0,
            # ... repeat for desired model complexity ...
        )
    )

Here is an illustration of this model's graph:

.. image:: https://github.com/rhayes777/PyAutoFit/blob/master/docs/overview/image/lens_model_cluster.png?raw=true
  :width: 600
  :alt: Alternative text

**PyAutoFit** therefore gives us full control over the composition and customization of high dimensional graphical
models.

Wrap-Up
-------

An example project on the **autofit_workspace** shows how to use **PyAutoFit** to set up code which fits strong
lensing data, using **multi-level model composition**.

If you'd like to perform the fit shown in this script, checkout the
`simple examples <https://github.com/Jammy2211/autofit_workspace/tree/master/notebooks/overview/simplee>`_ on the
``autofit_workspace``. We detail how **PyAutoFit** works in the first 3 tutorials of
the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_.

https://github.com/Jammy2211/autofit_workspace/tree/release/projects/astro.. _model_fit:

Fitting a Model
===============

To illustrate **PyAutoFit** we'll use the example modeling problem of fitting a 1D Gaussian profile to noisy data.

To begin, lets import ``autofit`` (and ``numpy``) using the convention below:

.. code-block:: bash

    import autofit as af
    import numpy as np

Data
----

The example ``data`` with errors (black) and the model-fit (red), are shown below:

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/images/toy_model_fit.png
  :width: 600
  :alt: Alternative text

Model
-----

We now need to define a 1D Gaussian profile as a **PyAutoFit** *model-component*, where a *model component* is a
component of the model we fit to the ``data``. It has associated with it a set of *parameters* which are varied for
during a *model-fit*, which is performed using a *non-linear search*.

*Model components* are defined using Python classes using the format below, where the class name is
the *model component* name and the constructor arguments are its *parameters*.

.. code-block:: bash

    class Gaussian:

        def __init__(
            self,
            centre=0.0,     # <- PyAutoFit recognises these
            normalization=0.1,  # <- constructor arguments are
            sigma=0.01,     # <- the Gaussian's parameters.
        ):

            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

The code above defines a **PyAutoFit** *model component* called a ``Gaussian``. When used for *model-fitting* it has
three parameters: ``centre``, ``normalization`` and ``sigma``.

When we fit the model to ``data`` and compute a likelihood an instance of the class above is accessible, with specific
values of ``centre``, ``normalization`` and ``sigma`` chosen by the non-linear search algorithm that fits the model to
the data.

This means that the class's functions are available to compute the likelihood, so lets add a ``profile_1d_via_xvalues_from``
function that generates the 1D profile from the ``Gaussian``.

.. code-block:: bash

    class Gaussian:
        def __init__(
            self,
            centre=0.0,     # <- PyAutoFit recognises these
            normalization=0.1,  # <- constructor arguments are
            sigma=0.01,     # <- the Gaussian's parameters.
        ):

            self.centre = centre
            self.normalization = normalization
            self.sigma = sigma

        def profile_1d_via_xvalues_from(self, xvalues):

            transformed_xvalues = xvalues - self.centre

            return np.multiply(
                np.divide(self.normalization, self.sigma * np.sqrt(2.0 * np.pi)),
                np.exp(-0.5 * np.square(np.divide(transformed_xvalues, self.sigma))),
            )

We use the ``Model`` object to compose the model, which in this case is a single ``Gaussian``.  The model is
defined with 3 free parameters, thus the dimensionality of non-linear parameter space is 3.

.. code-block:: bash

    model = af.Model(Gaussian)

Complex high dimensional models can be built from these individual model components, as described in
the `model composition overview page <https://pyautofit.readthedocs.io/en/latest/overview/model_complex.html>`_

Analysis
--------

Now we've defined our model, we need to tell **PyAutoFit** how to fit the model to data. This requires us to
define a **PyAutoFit** ``Analysis`` class:

.. code-block:: bash

    class Analysis(af.Analysis):

        def __init__(self, data, noise_map):

            super().__init__()

            self.data = data
            self.noise_map = noise_map

        def log_likelihood_function(self, instance):

            """
            The 'instance' that comes into this method is an instance of the Gaussian
            class, whose parameters were chosen by our non-linear search.

            The the print statements below will illustrate this when a model-fit is performed!
            """

            print("Gaussian Instance:")
            print("Centre = ", instance.centre)
            print("normalization = ", instance.normalization)
            print("Sigma = ", instance.sigma)

            """
            Get the range of x-values the data is defined on, to evaluate the model
            of the Gaussian.
            """

            xvalues = np.arange(self.data.shape[0])

            """
            Use these xvalues to create model_data of our Gaussian.
            """

            model_data = instance.profile_1d_via_xvalues_from(xvalues=xvalues)

            """
            Fit the model gaussian to the data, computing the residuals, chi-squareds
            and returning the log likelihood value to the non-linear search.
            """

            residual_map = self.data - model_data
            chi_squared_map = (residual_map / self.noise_map) ** 2.0
            log_likelihood = -0.5 * sum(chi_squared_map)

            return log_likelihood

Lets consider exactly what is happening in the ``Analysis`` class above.

- The ``data`` is passed into the constructor of the ``Analysis`` class. Above, only ``data`` and a ``noise_map`` are
  input, but the constructor can be easily extended to add other parts of the dataset.

- The ``log_likelihood_function`` receives an ``instance`` of the model, which in this example is an ``instance`` of the
  ``Gaussian`` class. This ``instance`` has values for its *parameters* (``centre``, ``normalization`` and ``sigma``) which
  are chosen by the non-linear search used to fit the model, as discussed next.

- The ``log_likelihood_function`` returns a log likelihood value, which the non-linear search uses evaluate the
  goodness-of-fit of a model to the data when sampling parameter space.

Non-Linear Search
-----------------

Next, we *compose* our model, set up our ``Analysis`` and fit the model to the ``data`` using a non-linear search:

.. code-block:: bash

    model = af.Model(Gaussian)
    analysis = Analysis(data=data, noise_map=noise_map)

    emcee = af.Emcee(name="example_search")

    result = emcee.fit(model=model, analysis=analysis)

We perform the fit using the non-linear search algorithm `emcee <https://github.com/dfm/emcee>`_. We cover
non-linear search's in more detail in the `non-linear search overview page <https://pyautofit.readthedocs.io/en/latest/overview/non_linear_search.html>`_.

Result
------

By running the code above **PyAutoFit** performs the model-fit, outputting all results into structured paths on you
hard-disk. It also returns a ``Result`` object in Python, which includes lists containing the non-linear search's
parameter samples, the maximum likelihood model, marginalized parameters estimates, errors are so on:

.. code-block:: bash

    print(result.samples.parameter_lists)
    print(result.samples.max_log_likelihood_vector)
    print(result.samples.median_pdf_vector)
    print(result.samples.error_vector_at_sigma)

It can even return *instances* of the ``Gaussian`` class using the values of the model results:

.. code-block:: bash

    instance = result.max_log_likelihood_instance

    print("Maximum Likelihood Gaussian Instance:")
    print("Centre = ", instance.centre)
    print("normalization = ", instance.normalization)
    print("Sigma = ", instance.sigma)

This can be used to straight forwardly plot the model fit to the data:

.. code-block:: bash

    instance = result.max_log_likelihood_instance

    model_data = instance.profile_1d_via_xvalues_from(xvalues=np.arange(data.shape[0]))

    plt.plot(range(data.shape[0]), data)
    plt.plot(range(data.shape[0]), model_data)

Results are covered in more detail in the `result overview page <https://pyautofit.readthedocs.io/en/latest/overview/result.html>`_.

Wrap-Up
-------

This completes our introduction to the **PyAutoFit** API. Next, we'll cover how to *compose* and *fit*
models using multiple *model components* and *customize* the model parameterization.

If you'd like to perform the fit shown in this script, checkout the
`simple examples <https://github.com/Jammy2211/autofit_workspace/tree/master/notebooks/overview/simplee>`_ on the
``autofit_workspace``. We detail how **PyAutoFit** works in the first 3 tutorials of
the `HowToFit lecture series <https://pyautofit.readthedocs.io/en/latest/howtofit/howtofit.html>`_.=============
API Reference
=============

.. currentmodule:: autofit

-------------------
Non-Linear Searches
-------------------

**Nested Samplers:**

.. autosummary::
   :toctree: generated/

   DynestyDynamic
   DynestyStatic
   MultiNest
   UltraNest

**MCMC:**

.. autosummary::
   :toctree: generated/

   Emcee
   Zeus

**Optimizers:**

.. autosummary::
   :toctree: generated/

   PySwarmsLocal
   PySwarmsGlobal

**GridSearch**:

.. autosummary::
   :toctree: generated/

   SearchGridSearch
   GridSearchResult

**Tools**:

.. autosummary::
   :toctree: generated/

   DirectoryPaths
   DatabasePaths
   Result
   InitializerBall
   InitializerPrior
   PriorPasser
   AutoCorrelationsSettings

--------
Plotters
--------

.. currentmodule:: autofit.plot
.. autosummary::
   :toctree: generated/

   DynestyPlotter
   UltraNestPlotter
   EmceePlotter
   ZeusPlotter
   PySwarmsPlotter


------
Models
------

.. currentmodule:: autofit

.. autosummary::
   :toctree: generated/

   PriorModel
   CollectionPriorModel


--------
Analysis
--------

.. currentmodule:: autofit

.. autosummary::
   :toctree: generated/

   Analysis


------
Priors
------

.. autosummary::
   :toctree: generated/

   UniformPrior
   GaussianPrior
   LogUniformPrior


-------
Samples
-------

.. autosummary::
   :toctree: generated/

   Samples
   PDFSamples
   MCMCSamples
   NestSamples
   StoredSamples

----------
Aggregator
----------

.. autosummary::
   :toctree: generated/

   Aggregator

-------
Backend
-------

.. autosummary::
   :toctree: generated/

   ModelMapper
   ModelInstance
.. _expectation_propagation:

Expectation Propagation
-----------------------

For large datasets, a graphical model may have hundreds, thousands, or *hundreds of thousands* of parameters. The
high dimensionality of such a parameter space can make it inefficient or impossible to fit the model.

Graphical models in **PyAutoFit** support the message passing framework below, which allows one to fit the local model
to every dataset individually and pass messages 'up and down' the graph to infer the global parameters efficiently.
Again, this feature is **still in beta** so contact us if you are interested in using this functionality ( https://github.com/Jammy2211 ).

https://arxiv.org/pdf/1412.4869.pdf.. _multiple_datasets:

Multiple Datasets
-----------------

**NOTE: Graphical models are an in-development feature. This example serves to illustrate what we currently developing , but the API is subject to change. If you are interested in using graphical models I recommend you contact me directly ( https://github.com/Jammy2211 ) so we can discuss how to implement **PyAutoFit** for your use-case.**

For graphical model composition, we saw how to compose a graphical model which fitted a single dataset. In this 
example, we will show how to build a graphical model that fits multiple datasets. 

It is common in statistical inference for us to have a large dataset and not be interested in how well a small aspect 
of the model fits each dataset individually. Instead, we want to fit the complete model to our full dataset and 
determine the global behaviour of the model's fit to every dataset.

Using graphical models, **PyAutoFit** can compose and fit models that have 'local' parameters specific to an individual 
dataset and higher-level model components that fit 'global' parameters. These higher level parameters will have 
conditional dependencies with the local parameters.  

The major selling point of **PyAutoFit**'s graphical modeling framework is the high level of customization it offers,
whereby: 

- Specific ``Analysis`` classes can be defined for fitting differnent local models to different datasets.
- Each pairing of a local model-fit to data can be given its own non-linear search.
- Graphical model networks of any topology can be defined and fitted.

In this example, we demonstrate the API for composing and fitting a graphical model to multiple-datasets, using the 
simple example of fitting noisy 1D Gaussians. In this example, I will explicitly write code that stores each dataset 
as its own Python variable (e.g. data_0, data_1, data_2, etc.), as opposed to a for loop or list. This is to make the 
API shown in this example clear, however examples in the ``autofit_workspace`` will use **PyAutoFit**'s bespoke API 
for setting up a graphical model.

We begin by loading noisy 1D data containing 3 Gaussian's.

.. code-block:: bash

    dataset_path = path.join("dataset", "example_1d")

    dataset_0_path = path.join(dataset_path, "gaussian_x1_0__low_snr")
    data_0 = af.util.numpy_array_from_json(file_path=path.join(dataset_0_path, "data.json"))
    noise_map_0 = af.util.numpy_array_from_json(
        file_path=path.join(dataset_0_path, "noise_map.json")
    )

    dataset_1_path = path.join(dataset_path, "gaussian_x1_1__low_snr")
    data_1 = af.util.numpy_array_from_json(file_path=path.join(dataset_1_path, "data.json"))
    noise_map_1 = af.util.numpy_array_from_json(
        file_path=path.join(dataset_1_path, "noise_map.json")
    )

    dataset_2_path = path.join(dataset_path, "gaussian_x1_2__low_snr")
    data_2 = af.util.numpy_array_from_json(file_path=path.join(dataset_2_path, "data.json"))
    noise_map_2 = af.util.numpy_array_from_json(
        file_path=path.join(dataset_2_path, "noise_map.json")
    )

This is what our three Gaussians look like. They are much lower signal-to-noise than the Gaussian's in other examples. 
We use lower signal-to-noise Gaussian's to demonstrate how fitting graphical models to lower quality data can still 
enable global parameters to be estimated precisely.

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_1__low_snr.png
  :width: 600
  :alt: Alternative text

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_2__low_snr.png
  :width: 600
  :alt: Alternative text

.. image:: https://raw.githubusercontent.com/rhayes777/PyAutoFit/master/docs/features/images/gaussian_x1_3__low_snr.png
  :width: 600
  :alt: Alternative text

For each dataset we now create a corresponding ``Analysis`` class. By associating each dataset with an ``Analysis``
class we are therefore associating it with a unique ``log_likelihood_function``. If our dataset had many different
formats (e.g. images) it would be straight forward to write customized ``Analysis`` classes for each dataset.

.. code-block:: bash

    analysis_0 = a.Analysis(data=data_0, noise_map=noise_map_0)
    analysis_1 = a.Analysis(data=data_1, noise_map=noise_map_1)
    analysis_2 = a.Analysis(data=data_2, noise_map=noise_map_2)

We now compose the graphical model we will fit using the ``Model`` and ``Collection`` objects. We begin by setting up a
shared prior for their ``centre`` using a single ``GaussianPrior``. This is passed to a unique ``Model`` for
each ``Gaussian`` and means that all three ``Gaussian``'s are fitted wih the same value of ``centre``. That is, we have
defined our graphical model to have a shared value of ``centre`` when it fits each dataset.

.. code-block:: bash

    centre_shared_prior = af.GaussianPrior(mean=50.0, sigma=30.0)

We now set up three ``Model`` objects, each of which contain a ``Gaussian`` that is used to fit each of the
datasets we loaded above. Because all three of these ``Model``'s use the ``centre_shared_prior`` the dimensionality of
parameter space is N=7, corresponding to three ``Gaussians`` with local parameters (``normalization`` and ``sigma``) and
a global parameter value of ``centre``.

.. code-block:: bash

    gaussian_0 = af.Model(m.Gaussian)
    gaussian_0.centre = centre_shared_prior
    gaussian_0.normalization = af.GaussianPrior(mean=10.0, sigma=10.0)
    gaussian_0.sigma = af.GaussianPrior(mean=10.0, sigma=10.0)  # This prior is used by all 3 Gaussians!

    gaussian_1 = af.Model(m.Gaussian)
    gaussian_1.centre = centre_shared_prior
    gaussian_1.normalization = af.GaussianPrior(mean=10.0, sigma=10.0)
    gaussian_1.sigma = af.GaussianPrior(mean=10.0, sigma=10.0)  # This prior is used by all 3 Gaussians!

    gaussian_2 = af.Model(m.Gaussian)
    gaussian_2.centre = centre_shared_prior
    gaussian_2.normalization = af.GaussianPrior(mean=10.0, sigma=10.0)
    gaussian_2.sigma = af.GaussianPrior(mean=10.0, sigma=10.0)  # This prior is used by all 3 Gaussians!

To build our graphical model which fits multiple datasets, we simply pair each model-component to each ``Analysis``
class, so that **PyAutoFit** knows that:

- ``gaussian_0`` fits ``data_0`` via ``analysis_0``.
- ``gaussian_1`` fits ``data_1`` via ``analysis_1``.
- ``gaussian_2`` fits ``data_2`` via ``analysis_2``.

The point where a ``Model`` and ``Analysis`` class meet is called a ``ModelFactor``.

This term is used to denote that we are composing a 'factor graph'. A factor defines a node on this graph where we have
some data, a model, and we fit the two together. The 'links' between these different factors then define the global
model we are fitting **and** the datasets used to fit it.

.. code-block:: bash

    model_factor_0 = g.ModelFactor(prior_model=prior_model_0, analysis=analysis_0)
    model_factor_1 = g.ModelFactor(prior_model=prior_model_1, analysis=analysis_1)
    model_factor_2 = g.ModelFactor(prior_model=prior_model_2, analysis=analysis_2)

We combine our ``ModelFactors`` into one, to compose the factor graph.

.. code-block:: bash

    factor_graph = g.FactorGraphModel(model_factor_0, model_factor_1, model_factor_2)

So, what does our factor graph looks like? Unfortunately, we haven't yet build visualization of this into **PyAutoFit**,
so you'll have to make do with a description for now. The factor graph above is made up of two components:

**Nodes:** these are points on the graph where we have a unique set of data and a model that is made up of a subset of
our overall graphical model. This is effectively the ``ModelFactor`` objects we created above.

**Links:** these define the model components and parameters that are shared across different nodes and thus retain the
same values when fitting different datasets.

.. code-block:: bash

    opt = g.optimise.LaplaceOptimiser(n_iter=3)
    model = factor_graph.optimise(opt)

**Road Map**

The example above which illustrated a simple graphical model where 3 datasets are fitted is fully functional in
**PyAutoFit** and can be ran at the following scripts:

https://github.com/Jammy2211/autofit_workspace/blob/release/scripts/features/graphical_models.py

https://github.com/Jammy2211/autofit_workspace/blob/release/scripts/howtofit/chapter_graphical_models/tutorial_2_graphical_model.py

However, graphical models **are still in beta testing** and I recommend you contact us if you wish to use the
functionality first (https://github.com/Jammy2211).