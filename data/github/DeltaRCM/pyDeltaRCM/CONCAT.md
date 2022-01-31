---
title: '*pyDeltaRCM*: a flexible numerical delta model'
tags:
  - Python
  - sedimentology
  - deltas
  - stratigraphy
authors:
  - name: Andrew J. Moodie
    orcid: 0000-0002-6745-036X
    affiliation: "1"
  - name: Jayaram Hariharan
    orcid: 0000-0002-1343-193X
    affiliation: "1"
  - name: Eric Barefoot
    orcid: 0000-0001-5770-2116
    affiliation: "2"
  - name: Paola Passalacqua
    affiliation: "1"
    orcid: 0000-0002-4763-7231
affiliations:
  - name: Department of Civil, Architectural, and Environmental Engineering, University of Texas at Austin, Austin, TX, USA
    index: 1
  - name: Department of Earth, Environmental and Planetary Sciences, Rice University, Houston, TX, USA
    index: 2
date: 08 June 2021
bibliography: paper.bib
---

# Summary

River deltas provide many societal benefits, and sustainability of these landforms may be impacted by human modification and global climate change.
Reduced-complexity numerical delta models incorporate limited physical processes, allowing researchers to assess the spatiotemporal evolution of landscape response to individual processes and environmental forcings.
Isolating individual processes is useful to understand, for example, shifting delta morphology due to sea-level rise, changing vegetal cover, or flooding intensity.
As a result, many numerical delta models have been proposed in the literature, and results from these studies are difficult to compare because of various design and implementation choices.
*pyDeltaRCM* (`v2.0`) delivers a computationally efficient and easy-to-customize implementation of the DeltaRCM numerical model [@liang_reduced_1_2015], enabling comparison and reproducibility in studies of delta change due to various environmental forcings.


# Statement of need

River deltas are societally important landforms because they provide arable land, deep inland ports, and are home to hundreds of millions of people globally [@edmonds_coastal_2020].
Existing at the interface between landmasses and water bodies, deltas are impacted by a multitude of processes arising in both of these domains.
For example, changes in sediment input to the delta modulate the rate at which new land is built; similarly, rising water levels in the downstream basin create flooded land.
In addition to natural processes, human landscape modification renders deltaic environments more sensitive to global climate change into the future [@paola_natural_2011].
Demand to understand natural delta processes, and how these processes will respond to various  environmental forcings, has led to a proliferation of numerical delta models in the literature [@overeem_three_2005].

The DeltaRCM delta model [@liang_reduced_1_2015] has gained popularity among geomorphologists due to an attractive balance of computational cost, realism, and interpretability [@larsen_appropriate_2016].
For example, studies have employed the DeltaRCM design to examine delta morphology and dynamism response to sea-level rise and regional subsidence [@liang_quantifying_2016; @liang_how_2016], as well as extended model design to simulate delta evolution with vegetation [@lauzon_comparing_2018] and ice and permafrost [@lauzon_ice_2019; @piliouras_unraveling_2021].
However, comparison among these studies is difficult, owing to disparate code bases, various implementation choices, lack of version control, and proprietary software dependencies.


# Background

Here, version 2.0 of *pyDeltaRCM* is introduced; *pyDeltaRCM* is a computationally efficient, free and open source, and easy-to-customize numerical delta model based on the original DeltaRCM design.
The original DeltaRCM framework is inspired by well-understood physical phenomena, and models mass movement as a probabilistic weighted random-walk process coupled with a set of hierarchical rules; the model is extensively described in @liang_reduced_1_2015 and @liang_reduced_2_2015.

This same framework is the basis for *pyDeltaRCM* v2.0, with a few modifications selected only to resolve known numerical instabilities, improve computational efficiency, and support reproducible simulations.
*PyDeltaRCM* depends only on common Python packages `numpy` [@harris2020], `matplotlib` [@hunter2007], `scipy` [@virtanen2020], `netCDF4`, `pyyaml`, and `numba` [@lam_numba_2015].

![Simulation with *pyDeltaRCM* v2.0, default parameter set, and random `seed: 10151919`. Simulation was run for 4000 timesteps, and assumes 10 days of bankfull discharge per year; computational time was \~2 hours. \label{fig:timeseries}](figures/timeseries.png)


# Flexible and easy to use

*pyDeltaRCM* is an object-oriented package, providing the central model class `DeltaModel`.
By creating custom model behavior as subclasses of `DeltaModel`, researchers can easily add, subtract, and modify model components without altering code that is not pertinent to the science objective.
Importantly, separating custom code from core model code makes clear how different studies can be compared.
The *pyDeltaRCM* documentation provides several examples for how to implement custom model behavior on top of the core `DeltaModel` object.

*pyDeltaRCM* also provides infrastructure to accelerate scientific exploration, such as the ability to configure multiple simulations from a single file.
Additionally, a preprocessor orchestrates `parallel` simulations for multi-core systems (optionally), and implements several tools to support simulations exploring a parameter space.
For example, `matrix` expansion converts lists of parameters into an n-dimensional set of simulations.
Similarly, replicate simulations can be created via an `ensemble` specification.

Reproducibility and computational efficiency were important priorities in *pyDeltaRCM* development.
For example, to-disk logging records all parameters, system-level and version data, and random-seed information to ensure that all runs can be recreated.
Additionally, "checkpoint" infrastructure has been added to the model, which records simulation progress during computation and can later resume model runs for further simulation.
Finally, *pyDeltaRCM* uses `numba` for computational optimization [@lam_numba_2015], and does not depend on any proprietary software.

*pyDeltaRCM* component units and integrations are thoroughly documented and tested.
Component-level documentation describes implementation notes, whereas narratives in "Guide" and "Example" documentation describes high-level model design and best practices for model use and development.
*pyDeltaRCM* also couples with other numerical models via the CSDMS Basic Model Interface 2.0 [@hutton_basic_2020; @BMI_pyDeltaRCM].


# Acknowledgments

We gratefully acknowledge Rebecca Lauzon and Mariela Perignon for developing an implementation of *DeltaRCM* in Python that was the basis for *pyDeltaRCM*. 
We also thank the National Science Foundation for supporting us in developing this software, by way of a Postdoctoral Fellowship to A.M. (EAR 1952772) and a grant to J.H. and P.P. (EAR 1719670).


# References***************
Code of Conduct
***************

Our Pledge
----------

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.


Our Standards
-------------

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members


Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting


Our Responsibilities
--------------------

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
-----

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at the issue tracker. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
-----------

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org/>`_, `version 1.4 <https://www.contributor-covenant.org/version/1/4/>`_,
available at https://www.contributor-covenant.org/version/1/4/.
**************
pyDeltaRCM
**************

.. image:: https://badge.fury.io/py/pyDeltaRCM.svg
    :target: https://badge.fury.io/py/pyDeltaRCM

.. image:: https://joss.theoj.org/papers/10.21105/joss.03398/status.svg
   :target: https://doi.org/10.21105/joss.03398

.. image:: https://github.com/DeltaRCM/pyDeltaRCM/actions/workflows/build.yml/badge.svg
    :target: https://github.com/DeltaRCM/pyDeltaRCM/actions
    
.. image:: https://codecov.io/gh/DeltaRCM/pyDeltaRCM/branch/develop/graph/badge.svg
  :target: https://codecov.io/gh/DeltaRCM/pyDeltaRCM

.. image:: https://app.codacy.com/project/badge/Grade/1c137d0227914741a9ba09f0b00a49a7
    :target: https://www.codacy.com/gh/DeltaRCM/pyDeltaRCM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeltaRCM/pyDeltaRCM&amp;utm_campaign=Badge_Grade    
    



*pyDeltaRCM* is a computationally efficient, free and open source, and easy-to-customize numerical delta model based on the original DeltaRCM model design (`Matlab deltaRCM <https://csdms.colorado.edu/wiki/Model:DeltaRCM>`_ model by Man Liang; `Liang et al., 2015 <https://doi.org/10.5194/esurf-3-67-2015>`_).
*pyDeltaRCM* delivers improved model stability and capabilities, infrastructure to support exploration with minimal boilerplate code, and establishes an approach to extending model capabilities that ensures reproducible and comparable studies.


.. figure:: https://deltarcm.org/pyDeltaRCM/_images/cover.png
    
    Weighted random walks for 20 water parcels, in a *pyDeltaRCM* model run with default parameters.


Documentation
#############

`Find the complete documentation here <https://deltarcm.org/pyDeltaRCM/index.html>`_.

Documentation includes an `installation guide <https://deltarcm.org/pyDeltaRCM/meta/installing.html>`_, a thorough `guide for users <https://deltarcm.org/pyDeltaRCM/guides/user_guide.html>`_, detailed `API documentation for developers <https://deltarcm.org/pyDeltaRCM/reference/index.html>`_, a `plethora of examples <https://deltarcm.org/pyDeltaRCM/examples/index.html>`_ to use and develop pyDeltaRCM in novel scientific experiments, and more!


Installation
############

See our complete `installation guide <https://deltarcm.org/pyDeltaRCM/meta/installing.html>`_, especially if you are a developer planning to modify or contribute code (`developer installation guide <https://deltarcm.org/pyDeltaRCM/meta/installing.html#developer-installation>`_), or if you are new to managing Python `venv` or `conda` environments.

For a quick installation into an existing Python 3.x environment:

.. code:: console

    $ pip install pyDeltaRCM


Executing the model
###################

We recommend you check out our `pyDeltaRCM in 10 minutes tutorial <https://deltarcm.org/pyDeltaRCM/guides/10min.html>`_, which is part of our documentation.

Beyond that breif tutorial, we have a comprehensive `User Documentation <https://deltarcm.org/pyDeltaRCM/index.html#user-documentation>`_ and `Developer Documentation <https://deltarcm.org/pyDeltaRCM/index.html#developer-documentation>`_ to check out.


Additional notes
################

This repository no longer includes the `Basic Model Interface (BMI) <https://bmi.readthedocs.io/en/latest/?badge=latest>`_ wrapper to the DeltaRCM model.
*pyDeltaRCM* maintains BMI compatibility through another repository (`the BMI_pyDeltaRCM model <https://deltarcm.org/BMI_pyDeltaRCM/>`_).
.. pyDeltaRCM documentation master file

Welcome to the PyDeltaRCM documentation
#########################################

*pyDeltaRCM* is a computationally efficient, free and open source, and easy-to-customize numerical delta model based on the original DeltaRCM model design (`Matlab DeltaRCM <https://csdms.colorado.edu/wiki/Model:DeltaRCM>`_ model by Man Liang).
*pyDeltaRCM* delivers improved model stability and capabilities, infrastructure to support exploration with minimal boilerplate code, and establishes an approach to extending model capabilities that ensures reproducible and comparable studies.

.. plot:: guides/cover.py

    The weighted random walk of 20 water parcels in a *pyDeltaRCM* model run with default parameters.


Project information
###################

.. toctree::
   :maxdepth: 1

   meta/installing
   meta/contributing
   meta/license
   meta/conduct

.. _user_documentation:

User documentation
##################

.. toctree::
   :maxdepth: 2

   guides/getting_started
   guides/user_guide
   info/index
   examples/index



Developer documentation
#######################

.. toctree::
   :maxdepth: 3

   guides/developer_guide
   reference/index
**********
User Guide
**********

.. rubric:: Preface

All of the documentation in this page, and elsewhere assumes you have imported the `pyDeltaRCM` package as ``pyDeltaRCM``:

.. code:: python

    >>> import pyDeltaRCM

Additionally, the documentation frequently refers to the `numpy` and `matplotlib` packages.
We may not always explicitly import these packages throughout the documentation; we may refer to these packages by their common shorthand as well.

.. code:: python

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt


=================
Running the model
=================

Running a `pyDeltaRCM` model is usually accomplished via a combination of configuration files and scripts.
In the simplest case, actually, none of these are required; the following command (executed at a console) runs a `pyDeltaRCM` model with the default set of model parameters for 10 timesteps:

.. code:: console

    $ pyDeltaRCM --timesteps 10

Internally, this console command calls a Python `function` that runs the model via a `method`.
Python code and scripts are the preferred and most common model use case.
Let's do the exact same thing (run a model for 10 timesteps), but do so with Python code.

Python is an object-oriented programming language, which means that the `pyDeltaRCM` model is set up to be manipulated *as an object*.
In technical jargon, the `pyDeltaRCM` model is built out as a Python `class`, specifically, it is a class named ``DeltaModel``; so, **the actual model is created by instantiating** (making an instance of a class) **the** :obj:`DeltaModel`.
Python objects are instantiated by calling the `class` name, followed by parentheses (and potentially any input arguments), from a python script (or at the Python console):

.. code:: python

    >>> mdl = pyDeltaRCM.DeltaModel()

.. plot::
    :context: reset

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     mdl = pyDeltaRCM.DeltaModel(out_dir=output_dir)

Here, ``mdl`` is an instance of the :obj:`DeltaModel` class; ``mdl`` is *an actual model that we can run*.
Based on the :doc:`default parameters </reference/model/yaml_defaults>`, the model domain is configured as an open basin, with an inlet centered along one wall.
Default parameters are set to sensible values, and generally follow those described in [1]_.

.. plot::
    :context:
    
    >>> fig, ax = plt.subplots(figsize=(4, 3))
    >>> ax.imshow(mdl.bed_elevation, vmax=-3)
    >>> plt.show()

To run the model, we use the :meth:`~pyDeltaRCM.model.DeltaModel.update` method of the :obj:`~pyDeltaRCM.model.DeltaModel` instance (i.e., call :meth:`~pyDeltaRCM.model.DeltaModel.update` on ``mdl``):

.. code:: python

    >>> for _ in range(0, 5):
    ...     mdl.update()

The :meth:`~pyDeltaRCM.model.DeltaModel.update` method manages a single timestep of the model run.
Check out the docstring of the method for a description of the routine, and additional steps that are called from :obj:`update`.
Additionally, it may be helpful to read through the :doc:`model informational guides </info/index>`.

.. note::

    The model timestep (:attr:`~pyDeltaRCM.model.DeltaModel.dt`) is determined by model stability criteria; see :doc:`the model stability guide </info/morphodynamics>` for more information.

Finally, we need to tie up some loose ends when we are done running the model.
We use the :obj:`~pyDeltaRCM.model.DeltaModel.finalize` method of the :obj:`DeltaModel` instance (i.e., call :meth:`~pyDeltaRCM.model.DeltaModel.finalize` on ``mdl``):

.. code:: python

    >>> mdl.finalize()

Specifically, calling :meth:`~pyDeltaRCM.model.DeltaModel.finalize` will ensure that any data output during the model run to the :doc:`output NetCDF4 file </info/outputfile>` is correctly saved to the disk and does not become corrupted while the script is exiting.
Note that the output file is periodically saved during the model run, so things might be okay if you forget to :meth:`finalize` at the end of your simulation, but it is considered a best practice to explicitly close the output file with :meth:`finalize`.

Putting the above snippets together gives a complete **minimum working example script**:

.. code-block:: python

    >>> mdl = pyDeltaRCM.DeltaModel()

    >>> for _ in range(0, 5):
    ...    mdl.update()

    >>> mdl.finalize()

Now, we probably don't want to just run the model with default parameters, so the remainder of the guide will cover other use cases with increasing complexity; first up is changing model parameter values via a ``YAML`` configuration file.


==============================
Configuring an input YAML file
==============================

The configuration for a pyDeltaRCM run is set up by a parameter set, usually described in the ``YAML`` markup format.
To configure a run, you should create a file called, for example, ``model_configuration.yml``.
Inside this file you can specify parameters for your run, with each parameter on a new line. For example, if ``model_configuration.yml`` contained the line:

.. code-block:: yaml

    S0: 0.005
    seed: 42

then a :obj:`~pyDeltaRCM.DeltaModel` model instance initialized with this file specified as ``input_file`` will have a slope of 0.005, and will use a random seed of 42.
Multiple parameters can be specified line by line.

Default values are substituted for any parameter not explicitly given in the ``input_file`` ``.yml`` file.
Default values of the YAML configuration are listed in the :doc:`/reference/model/yaml_defaults`.

.. important::

    The best practice for model configurations is to create a YAML file with only the settings you want to change specified. Hint: comment a line out with ``#`` so it will not be used in the model.

.. hint::

    Check the model log files to make sure your configuration was interpreted as you  expected!

===================
Starting model runs
===================

There are two API levels at which you can interact with the pyDeltaRCM model.
There is a "high-level" model API, which takes as argument a YAML configuration file, and will compose a list of jobs as indicated in the YAML file; the setup can be configured to automatically execute the job list, as well.
The "low-level" API consists of creating a model instance from a YAML configuration file and manually handling the timestepping, or optionally, augmenting operations of the model to implement new features.
Generally, if you are modifying the model source code or doing anything non-standard with your simulation runs, you will want to use the low-level API. The high-level API is extremely helpful for exploring parameter spaces.


.. _high_level_api:

High-level model API
====================

The high-level API is accessed via either a shell prompt or python script, and handles setting up the model configuration and running the model a specified duration.

For the following high-level API demonstrations, consider a YAML input file named ``model_configuration.yml`` which looks like:

.. code-block:: yaml

    S0: 0.005
    seed: 42
    timesteps: 500


Command line API
----------------

To invoke a model run from the command line using the YAML file ``model_configuration.yml`` defined above,
we would simply call:

.. code:: bash

    pyDeltaRCM --config model_configuration.yml

or equivalently:

.. code:: bash

    python -m pyDeltaRCM --config model_configuration.yml

These invocations will run the pyDeltaRCM :obj:`preprocessor <pyDeltaRCM.preprocessor.PreprocessorCLI>` with the parameters specified in the ``model_configuration.yml`` file.
If the YAML configuration indicates multiple jobs (:ref:`via matrix expansion or ensemble specification <configuring_multiple_jobs>`), the jobs will each be run automatically by calling :obj:`~pyDeltaRCM.DeltaModel.update` on the model 500 times.



Python API
----------

The Python high-level API is accessed via the :obj:`~pyDeltaRCM.Preprocessor` object.
First, the `Preprocessor` is instantiated with a YAML configuration file (e.g., ``model_configuration.yml``):

.. code:: python

    >>> pp = preprocessor.Preprocessor(p)

which returns an object containing the list of jobs to run.
Jobs are then run with:

.. code:: python

    >>> pp.run_jobs()


Model simulation duration in the high-level API
-----------------------------------------------

The duration of a model run configured with the high-level API can be set up with a number of different configuration parameters.

.. note:: see the complete description of "time" in the model: :doc:`../info/modeltime`.

Using the high-level API, you can specify the duration to run the model by two mechanisms: 1) the number of timesteps to run the model, or 2) the duration of time to run the model.

The former case is straightforward, insofar that the model determines the timestep duration and the high-level API simply iterates for the specified number of timestep iterations.
To specify the number of timesteps to run the model, use the argument ``--timesteps`` at the command line (or ``timesteps:`` in the configuration YAML file, or ``timesteps=`` with the Python :obj:`~pyDeltaRCM.Preprocessor`).

.. code:: bash
    
    $ pyDeltaRCM --config model_configuration.yml --timesteps 5000

The second case is more complicated, because the time specification is converted to model time according to a set of additional parameters.
In this case, the model run end condition is that the elapsed model time is *equal to or greater than* the specified input time.
Importantly, this means that the duration of the model run is unlikely to exactly match the input condition, because the model timestep is unlikely to be a factor of the specified time.
Again, refer to the complete description of model time :doc:`../info/modeltime` for more information.

To specify the duration of time to run the model in *seconds*, simply use the argument ``--time`` at the command line (or ``time:`` in the configuration YAML file, or ``time=`` with the Python :obj:`~pyDeltaRCM.Preprocessor`).
It is also possible to specify the input run duration in units of years with the similarly named argument ``--time_years`` (``time_years:``, ``time_years=``).

.. code:: bash
    
    $ pyDeltaRCM --config model_configuration.yml --time 31557600
    $ pyDeltaRCM --config model_configuration.yml --time_years 1

would each run a simulation for :math:`(86400 * 365.25)` seconds, or 1 year.

.. important::

    Do not specify both time arguments, or specify time arguments with the timesteps argument.
    In the case of multiple argument specification, precedence is given in the order `timesteps` > `time` > `time_years`.

When specifying the time to run the simulation, an additional parameter determining the intermittency factor (:math:`I_f`) may be specified ``--If`` at the command line (``If:`` in the YAML configuration file, ``If=`` with the Python :obj:`~pyDeltaRCM.Preprocessor`).
This argument will scale the specified time-to-model-time, such that the *scaled time* is equal to the input argument time.
Specifying the :math:`I_f`  value is essential when using the model duration run specifications.
See :doc:`../info/modeltime` for complete information on the scaling between model time and elapsed simulation time.

Running simulations in parallel
-------------------------------

The high-level API provides the ability to run simulations in parallel on Linux environments.
This option is only useful in the case where you are running multiple jobs with the :ref:`matrix expansion <matrix_expansion_tag>`, :ref:`ensemble expansion <ensemble_expansion_tag>`, or :ref:`set expansion <set_expansion_tag>` tools.

To run jobs in parallel simply specify the `--parallel` flag to the command line interface.
Optionally, you can specify the number of simulations to run at once by following the flag with a number.

.. code:: bash
    
    $ pyDeltaRCM --config model_configuration.yml --timesteps 5000 --parallel
    $ pyDeltaRCM --config model_configuration.yml --timesteps 5000 --parallel 6


Low-level model API
===================

The low-level API is the same as that described at the beginning of this guide.
Interact with the model by creating your own script, and manipulating model outputs at the desired level.
The simplest case to use the low-level API is to do

.. code::

    >>> delta = pyDeltaRCM.DeltaModel(input_file='model_configuration.yml')

    >>> for _ in range(0, 5000):
    ...    delta.update()

    >>> delta.finalize()

However, you can also inspect/modify the :obj:`~pyDeltaRCM.DeltaModel.update` method, and change the order of operations, or add operations, as desired; see the :ref:`guide to customizing the model <customize_the_model>` below.
If you are working with the low-level API, you can optionally pass any valid key in the YAML configuration file as a keyword argument during model instantiation. 
For example:

.. code::

    >>> delta = pyDeltaRCM.DeltaModel(input_file='model_configuration.yml',
    ...                    SLR=1e-9)


Keyword arguments supplied at this point will supersede values specified in the YAML configuration.
See :ref:`our guide for model customization <customize_the_model>` for a complete explanation and demonstration for how to modify model behavior.

..
    =============================
    Advanced model configurations
    =============================
    ** Advanced model configuration guide is imported from another file. **

.. include:: advanced_configuration_guide.inc


..
    =======================
    Working with Subsidence
    =======================
    ** Subsidence guide is imported from another file. **

.. include:: subsidence_guide.inc


..
    =======================
    dsfasf
    =======================
    ** Subsidence guide is imported from another file. **

.. include:: /info/outputfile.rst


Supporting documentation and files
==================================

Model reference:

    * :doc:`/reference/model/yaml_defaults`
    * :doc:`/reference/model/model_hooks`

Examples:

    * `Simple simulation in Jupyter Notebook <https://github.com/deltaRCM/pyDeltaRCM/blob/develop/docs/source/examples/simple_example.ipynb>`_
    * :doc:`/examples/slight_slope`
    * :doc:`/examples/subsidence_region`
    * :doc:`/examples/custom_saving`

References
==========

.. [1] A reduced-complexity model for river delta formation – Part 1: Modeling
       deltas with channel dynamics, M. Liang, V. R. Voller, and C. Paola, Earth
       Surf. Dynam., 3, 67–86, 2015. https://doi.org/10.5194/esurf-3-67-2015
***************
Developer Guide
***************

.. image:: https://github.com/DeltaRCM/pyDeltaRCM/workflows/actions/badge.svg
    :target: https://github.com/DeltaRCM/pyDeltaRCM/actions

This guide provides additional details for implementation that did not make it into the user guide.
If you have not yet read the user guide, that is the place to start, then visit and refer to this guide for details as necessary.


==============
Best Practices
==============

Reproducibility
---------------

.. important:: tl;dr: any random numbers must be generated within a jitted function. 

A major goal of the `pyDeltaRCM` project is to enable fully reproducible simulation results.
Variability in pyDeltaRCM arises from the weighted random walks of water and sediment parcels during model iteration, so the state of the random source is essential to reproducibility.
We encourage developers to ensure that their subclass models are also reproducible, and provide some information on how to do so in this section.

For many subclassing models, it will be straightforward to ensure models are reproducible.
Out of the box, models will use the core `pyDeltaRCM` "seed" functionality to make models reproducible, and checkpointing should easily integrate with most use cases.
However, models that implement processes or functions that rely on random numbers or samples from random distributions will need to take care to ensure models are reproducible.

Random numbers
~~~~~~~~~~~~~~

For `pyDeltaRCM`, the source of random values is a pseudo-random number generator (RNG) from `numpy`, which works by inputting the current state of the RNG to an algorithm, and returning a sample (an integer) that can be mapped to any probability density function.
Each time the RNG yields a sample, the RNG state is changed, such that repeated samples from the RNG appear to be random, but are actually a deterministic sequence for a given initial RNG state (i.e., for a given seed).

As a result of the deterministic RNG, model runs are *exactly* reproducible if they begin from the same initial RNG state.
However, this means that any change to the state of the RNG that is not recorded will make a run non-reproducible. 
Additionally, any "random" behavior implemented in the model *that does not also* modify the state of the RNG, is non-reproducible. 

With a default setup of `numpy`, the functions of `np.random` will utilize the same underlying random number generator, such that the following code always evaluates to the same results.
In this case, it would be relatively simple to keep runs reproducible, because any call to a `np.random` function would modify the state of the underlying generator.

.. doctest::

    >>> np.random.seed(10); np.random.uniform(); np.random.normal(); np.random.uniform();
    0.771320643266746
    0.03777261440227079
    0.7488038825386119

In `pyDeltaRCM` though, we use `numba` just-in-time compilation for several steps in the model routine to make execution faster, including getting the random values that drive parcel movement during model iteration.
Within a Python console, calls to the `numba` RNG will not affect the state of the `numpy` RNG, and vice versa; even though the pre-JIT compiled code appears to call `np.random.uniform` (e.g., :obj:`~pyDeltaRCM.shared_tools.get_random_uniform`).

.. important:: `pyDeltaRCM` only takes responsibility for the `numba` RNG!


A simple example
~~~~~~~~~~~~~~~~

So why does it matter?
Well, calling `np.random.normal(100, 10)` in a model subclass will give you a random value, and modify the state of the `numpy` RNG, but it will not affect the state of the `numba` RNG (the one `pyDeltaRCM` actually keeps track of).
Thus, the values returned from the call to `np.random.normal` are not known and are not reproducible.

In the following simple example, see how the reproducible model uses a random number generated from a jitted function. This ensures the `numba` RNG is used for random variability in the model, and runs are reproducible.


.. doctest:: 

    >>> from numba import njit

    >>> @njit
    ... def get_random_normal():
    ...     """Get a random number from standard normal distribution.
    ...     """
    ...     return np.random.normal(0, 1)
    

    >>> class BrokenAndNotReproducible(pyDeltaRCM.DeltaModel):
    ... 
    ...     def __init__(self, input_file=None, **kwargs):
    ...         """Initialize a model that can never be reproduced.
    ...         """
    ... 
    ...         super().__init__(input_file=input_file, **kwargs)
    ... 
    ...     def update(self):
    ...         """Reimplement update method for demonstration."""
    ... 
    ...         # the core pyDeltaRCM RNG is used in computations, e.g.,
    ...         _sample0 = pyDeltaRCM.shared_tools.get_random_uniform(1)
    ... 
    ...         # now, we do something custom in our subclass
    ...         _sample1 = np.random.normal(0, 1)
    ... 
    ...         # and write it out to view
    ...         print(_sample0, _sample1)
    

    >>> class BeautifulAndVeryReproducible(pyDeltaRCM.DeltaModel):
    ... 
    ...     def __init__(self, input_file=None, **kwargs):
    ...         """Initialize a reproducible model.
    ...         """
    ... 
    ...         super().__init__(input_file=input_file, **kwargs)
    ... 
    ...     def update(self):
    ...         """Reimplement update method for demonstration."""
    ... 
    ...         # the core pyDeltaRCM RNG is used in computations, e.g.,
    ...         _sample0 = pyDeltaRCM.shared_tools.get_random_uniform(1)
    ... 
    ...         # now, we do something custom in our subclass
    ...         _sample1 = get_random_normal()
    ... 
    ...         # and write it out to view
    ...         print(_sample0, _sample1)

Now, we will initialize and run each model for three timesteps, twice. Running each twice will allow us to see if the model is reproducible (i.e., are all of the numbers exactly the same between runs).
First, run the `Broken` model:

.. doctest::
    :hide:

    # generate this for good docs, but it is not shown
    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     broken = BrokenAndNotReproducible(
    ...         out_dir=output_dir, seed=10)

.. code::

    broken = BrokenAndNotReproducible(seed=10)

.. doctest::

    >>> for i in range(3):
    ...     broken.update() # doctest: +SKIP
    0.771320643266746 1.213653088541954
    0.0207519493594015 -0.40009453994985783
    0.6336482349262754 0.7719410676752912

.. doctest::
    :hide:

    # generate this for good docs, but it is not shown
    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     broken = BrokenAndNotReproducible(
    ...         out_dir=output_dir, seed=10)

.. code::

    broken = BrokenAndNotReproducible(seed=10)

.. doctest::

    >>> for i in range(3):
    ...     broken.update() # doctest: +SKIP
    0.771320643266746 -0.17926017697487434
    0.0207519493594015 -0.4421037872728855
    0.6336482349262754 -0.2725394596633578

Now, run the reproducible model:

.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     beautiful = BeautifulAndVeryReproducible(
    ...         out_dir=output_dir, seed=10)

.. code::

    beautiful = BeautifulAndVeryReproducible(seed=10)

.. doctest::

    >>> for i in range(3):
    ...     beautiful.update()
    0.771320643266746 0.03777261440227079
    0.7488038825386119 -0.1354484915560101
    0.4985070123025904 -0.6643797082723693


.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     beautiful = BeautifulAndVeryReproducible(
    ...         out_dir=output_dir, seed=10)

.. code::

    beautiful = BeautifulAndVeryReproducible(seed=10)

.. doctest::

    >>> for i in range(3):
    ...     beautiful.update()
    0.771320643266746 0.03777261440227079
    0.7488038825386119 -0.1354484915560101
    0.4985070123025904 -0.6643797082723693

From these results, we can see that the values returned from the built-in uniform RNG as the first sample of each iteration (i.e., the left column) is always deterministic (in `broken` and `beautiful`), whereas both the built-in and the custom RNG are deterministic (in `beautiful`).


.. important::
    
    Be sure to only generate random numbers inside jitted functions!

.. note::

    It is generally okay to not worry about reproducibility when you are developing your subclassing model and trying to work out how model mechanics will depend on randomness -- but once you start to do real simulations you may analyze, be sure to take the time to make your model reproducible.



Model development
-----------------

Slicing and neighbors 
~~~~~~~~~~~~~~~~~~~~~

Slicing an array to find the array values of neighbors is a common operation in the model.
The preferred way to slice is by 1) padding the array with :func:`~pyDeltaRCM.shared_tools.custom_pad`, and 2) looping through rows and columns to directly index. This approach makes for readable and reasonably fast code; for example, to find any cells that are higher than all neighbors:

.. code::
    
    pad_eta = shared_tools.custom_pad(self.eta)
    for i in range(self.L):
        for j in range(self.W):
            eta_nbrs = pad_eta[i - 1 + 1:i + 2 + 1, j - 1 + 1:j + 2 + 1]
            eta_nbrs[1, 1] = -np.inf
            
            np.all(self.eta[i, j] > eta_nbrs)

There are also several model attributes that may be helpful in development; we suggest using these builtins rather than creating your own whenever possible (see :meth:`~pyDeltaRCM.init_tools.init_tools.set_constants` and the model source code).
.. _getting_started:

===============
Getting Started
===============


It's easy to get started using *pyDeltaRCM*.
We recommend starting out with the 10-minute learn-by-example lesson:

.. toctree::
   :maxdepth: 3

   10min


Return to the main index :ref:`user_documentation` to see all the other documentation available to help you get started.

Some other helpful resources:

* An example `Jupyter Notebook <https://github.com/deltaRCM/pyDeltaRCM/blob/develop/docs/source/examples/simple_example.ipynb>`_ that shows a complete model runs and quick analysis of output file.
* A comprehensive :doc:`user_guide` on how to use and experiment with the model


And when you're ready to start developing, check out the :doc:`developer_guide` and the complete model :doc:`API documentation </reference/index>`.
******************
10-minute tutorial
******************

Use pyDeltaRCM in ten minutes!
This simple guide will show you the absolute basics of getting a `pyDeltaRCM` model running, and give you some direction on where to go from there.

.. important::

    If you haven't already, be sure to follow the :doc:`installation guide </meta/installing>` to get *pyDeltaRCM* set up properly on your computer.


A default model
---------------

You can get a model running with five simple lines of Python code. Note that you can run *pyDeltaRCM* in either a standalone script or part of an interactive session.
First, we instantiate the main :obj:`~pyDeltaRCM.DeltaModel` model object.

.. code:: python

    >>> import pyDeltaRCM

    >>> default_delta = pyDeltaRCM.DeltaModel()

Instantiating the :obj:`~pyDeltaRCM.model.DeltaModel()` without any arguments will use a set of :doc:`default parameters <../reference/model/yaml_defaults>` to configure the model run.
The default options are a reasonable set for exploring some controls of the model, and would work perfectly well for a simple demonstration here.

The delta model is run forward with a call to the :meth:`~pyDeltaRCM.DeltaModel.update()` method.
So, we simply create a `for` loop, and call the `update` function, and then wrap everything up with a call to :meth:`~pyDeltaRCM.DeltaModel.finalize` the model:

.. code:: python

    >>> for _ in range(0, 5):
    ...     default_delta.update()

    >>> default_delta.finalize()

.. note::

    Additional calls to update the model can be called up until the model is finalized.

That's it! You ran the pyDeltaRCM model for five timesteps, with just five lines of code.

We can visualize the delta bed elevation, though it's not very exciting after only five timesteps...

.. code:: python

    >>> import matplotlib.pyplot as plt

    >>> fig, ax = plt.subplots()
    >>> ax.imshow(default_delta.bed_elevation, vmax=-3)
    >>> plt.show()

.. plot:: guides/10min_demo.py


The model with set parameters
-----------------------------

To run a simulation with a non-default set of parameters, we use a configuration file written in the YAML markup language named `10min_tutorial.yaml`.
The markup file allows us to specify model boundary conditions and input and output settings, where anything set in the file will override the :doc:`default parameters <../reference/model/yaml_defaults>` for the model, and anything *not* specified will take the default value.

.. important::

    The best practice for model configurations is to create a YAML file with only the settings you want to change specified.

The YAML configuration file is central to managing *pyDeltaRCM* simulations, so we did not create this file for you; you will need to create the YAML file yourself.
To create the YAML file, open up your favorite plain-text editing application (e.g., gedit, notepad).
YAML syntax is pretty simple for basic configurations, essentially amounting to each line representing a parameter-value pair, separated by a colon.
For this example, let's specify three simulation controls: where we want the output file to be placed via the `out_dir` parameter, we will ensure that our simulation is easily reproducible by setting the random `seed` parameter, and we can examine what is the effect of a high fraction of bedload with the `f_bedload` parameter.
Enter the following in your text editor, and save the file as ``10min_tutorial.yaml``, making sure to place the file in a location accessible to your interpreter.

.. code:: yaml

    out_dir: '10min_tutorial'
    seed: 451220118313
    f_bedload: 0.9


Now, we can create a second instance of the :obj:`~pyDeltaRCM.model.DeltaModel()`, this time using the input yaml file.

.. code::

    >>> second_delta = pyDeltaRCM.DeltaModel(input_file='10min_tutorial.yaml')

and repeat the same `for` loop operation as above:

.. code:: python

    >>> for _ in range(0, 5):
    ...     second_delta.update()

    >>> second_delta.finalize()


Resources
---------

Consider reading through the :doc:`User Guide <user_guide>` as a first action, and determine how to set up the model to complete your experiment, including tutorials and examples for customizing the model to achieve any arbitrary behavior you need!

* :doc:`user_guide`
* :doc:`/reference/model/index`
************
Contributing
************

To contribute to this project, you must follow our :doc:`Code of Conduct </meta/conduct>` at all times.
If you are not familiar with our code of conduct policy, take a minute to read the policy before starting with your first contribution.


How to contribute
-----------------

We welcome contributions of many types, including but not limited to:

    * bug fixes
    * improvements to documentation
    * new features
    * additional examples for using and developing *pyDeltaRCM*

If you are interested in contributing, please submit a pull request or get in touch with the development team via the `Github issue tracker <https://github.com/DeltaRCM/pyDeltaRCM/issues>`_.


Issues and questions
--------------------

If you have identified a bug within the code, but aren't sure how to (or don't want to) fix it, we invite you to open an issue on our `Github issue tracker <https://github.com/DeltaRCM/pyDeltaRCM/issues>`_.
Please explain the problem as fully as possible, including any system information that may be relevant.

You can also ask any questions about the software or ask for help on the `Github issue tracker <https://github.com/DeltaRCM/pyDeltaRCM/issues>`_.
We will try our best to help you how we can!

.. include:: ../../../CODE_OF_CONDUCT.rst
************
Installing
************

We recommend installing *pyDeltaRCM* in a virtual environment.
That said, *pyDeltaRCM* depends on a small number of packages (:ref:`list of dependencies <dependencies-list>`), many of which are likely already in a Python user/developer's regular use, so it's probably safe to install *pyDeltaRCM* in your base environment, too.


Installing
==========

We describe installation flavors for both users and developers below.

.. hint::

    If you are looking to make any modifications to the model source code, you should follow the developer instructions.

We suggest using the Anaconda Python distribution, which you can obtain via `the project website <https://www.anaconda.com/products/individual>`_.

Before proceeding, you may wish to create a virtual environment for the *pyDeltaRCM* project.
With Anaconda on Linux:

.. code:: console

    $ conda create -n deltarcm python=3
    $ conda activate deltarcm

For more informtaion, see `this guide <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment>`_ for help on creating and activating a virtual environment with Anaconda on other platforms.
See `this helpful guide <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment>`_ for creating virtual environments with `venv` if you do not use Anaconda.


User installation
-----------------

For a user installation, simply install from the pypi package repository:

.. code:: console

    $ pip install pyDeltaRCM

.. note::

    You may need to `first install <https://pip.pypa.io/en/stable/installing/>`_ `pip`.


.. _dev-install:

Developer installation
----------------------

For a developer installation, you should first fork the repository on Github.
This will allow you to submit suggestions and contribute to *pyDeltaRCM*.

.. note::

    You do not *need* to create a fork if your are just testing, but it may save you time and headache down the road. If you choose not to, just use the main repository url below (https://github.com/DeltaRCM/pyDeltaRCM.git).

First, you will need to `install git <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_  if you do not already have it.
Then, download or clone your fork of the project:

.. code:: console

    $ git clone https://github.com/<your-username>/pyDeltaRCM.git

Then, with current working directory as the root of the repository (e.g., ``cd pyDeltaRCM``), run the following commands:

.. code:: console

    $ pip install -r requirements.txt
    $ pip install -r requirements-docs.txt
    $ pip install -r requirements-test.txt
    $ pip install -e .

To check installation, run the complete test suite with:

.. code:: console

    $ pytest --mpl --mpl-baseline-path=tests/imgs_baseline

Finally, add the `upstream` repository to your `remote` repository list:

.. code:: console

    $ git remote add upstream https://github.com/DeltaRCM/pyDeltaRCM.git

You can build a local copy of the documentation with:

.. code:: console

    $ (cd docs && make html)


Next steps
==========

Consider reading through the :doc:`10-minute tutorial </guides/10min>` or the :doc:`User Guide </guides/user_guide>` to help get you started using *pyDeltaRCM*.


.. _dependencies-list:

Dependencies
============

.. literalinclude:: ../../../requirements.txt
   :linenos:
******************
Time in pyDeltaRCM
******************

Time in the pyDeltaRCM model is simulated in units of seconds.
The duration of a timestep is determined during model initialization to maintain numerical stability; the calculation is based on model domain and configuration parameters.
The model is then iterated timestep by timestep, until an end-run condition is reached 

.. note:: 
    
    If you are using the :ref:`high-level API <high_level_api>`, the model run end condition is that the elapsed model time is *equal to or greater than* the specified input time. As a result, the model run duration is unlikely to exactly match the input time specification, because the model timestep is unlikely to be a factor of the specified time. 

    Please keep this in mind when evaluating model results, and especially when comparing  between different model runs.

Over the duration of the model, the water discharge (and thus sediment discharge) are assumed to be at bankfull.
This assumption is based on the concept of geomorphic work [WM60]_ and the strongly nonlinear relationship between water discharge and sediment transport.
To summarize the basis of this assumption, a very small proportion of sediment is in motion when discharge is at or near a low-flow "base flow" condition, relative to the amount of sediment in motion during flood conditions.
So while flooding only occurs for a small fraction of total time, most of the time-integrated sediment transport occurs during this fraction of time, and the remainder of time is deemed negligible in comparison.

For example, in the contrived river-delta hydrograph shown below, the discharge fluctuates around a base-flow condition for most of the year, with the exception of a 50 day period when a flood event occurs (from "start event" to "end event").
For 25 days of the flood event, the discharge is high enough to exceed the bankfull discharge ("flooding duration").

.. plot:: modeltime/one_year_plot.py

A numerical model could simulate every day of this hydrograph, but this would require substantial computing power.
To accelerate the computation time of models, we assume that the only important period of time is when the river is at-or-above bankfull discharge, and we collapse the hydrograph to a single value: the bankfull discharge.
The choice of bankfull discharge for the single-value discharge assumption is a matter of convenience, in that this is readily estimated from field data or measurements.

The duration of time represented by the single-value discharge assumption is also  determined by the hydrograph.
We arrive at the so-called *intermittency factor* (:math:`I_f`) by finding the fraction of unit-time per unit-time when the river is in flood.
For this example, the intermittency factor scales between days and years as 

.. math::

    \frac{25~\textrm{days}}{365.25~\textrm{days}} \approx 0.07 \equiv I_f

The intermittency factor thus gives a relationship to scale pyDeltaRCM model time (which is always computed in units of seconds elapsed) from seconds to years.
For example, for a model run that has elapsed :math:`4.32 \times 10^8` seconds and an assumed intermittency factor of 0.07:

.. math::

    \frac{\textrm{seconds~elapsed}}{I_f} = \frac{4.32 \times 10^8~\textrm{seconds}}{0.07} = 6.17 \times 10^9~\textrm{seconds} \approx 196~\textrm{years}

.. note:: A convenience function is supplied with the low-level API for these time conversions: :obj:`~pyDeltaRCM.shared_tools.scale_model_time`.

Applying the intermittency assumption to a four year model simulation reveals the computational advantage of this approach.
A simulation of every day in the simulation duration would require 1460 days of simulated time (:math:`1.26 \times 10^8` seconds).

.. plot:: modeltime/four_year_plot.py

In contrast, by condensing the simulated time down to a bankfull discharge only, the intermittency assumption allows us to simulate the same "four year" period in just under 59 days of simulated time (:math:`5.08 \times 10^6` seconds).


See also
--------

* :obj:`~pyDeltaRCM.preprocessor.scale_relative_sea_level_rise_rate`
* :obj:`~pyDeltaRCM.shared_tools.scale_model_time`


References
----------

.. [WM60] Wolman, M. G., & Miller, J. P. (1960). Magnitude and frequency of forces in geomorphic processes. The Journal of Geology, 68(1), 54–74. https://doi.org/10.1086/626637=================
Model Output File
=================

If configured to save any output data, model outputs are saved using the `netCDF4 <http://unidata.github.io/netcdf4-python/>`_ file format.


Gridded Variables
=================

In any given run, the saving parameters "save_<var>_grids" control whether or
not that 2-D grid variable (e.g. velocity) is saved to the netCDF4 file. In
the netCDF4 file, a 3-D array with the dimensions `time` :math:`\times`
`x` :math:`\times` `y` is created for each 2-D grid variable that is set to
be saved. Note that `x` is the *downstream* coordinate, rather than the
Cartesian `x` when displaying the grid. The appropriate units for all
variables are stored: for example "meters per second" for the *velocity*
grid.

.. note::
   
   The format of the output netCDF file coordinate changed in `v2.1.0`. The
   old format is documented
   in :attr:`~pyDeltaRCM.model.DeltaModel.legacy_netcdf`, and that input
   parameter `legacy_netcdf` can be used to create on output netcdf file with
   the old coordinate configuration.


Grid Coordinates
================

Grid coordinates are specified in the variables `time`, `x`, and `y` in the output netCDF4 file.
These arrays are 1D arrays, which specify the location of each cell in the domain in *dimensional* coordinates (e.g., meters).
In the downstream direction,  the distance of each cell from the inlet boundary is specified in `x` in meters.
Similarly, the cross-domain distance is specified in `y` in meters.
Lastly, the `time` variable is stored as a 1D array with model `time` in seconds.


Model Metadata
==============

In addition to the grid coordinates, model metadata is saved as a group of
1-D arrays (vectors) and 0-D arrays (floats and integers). The values that are
saved as metadata are the following:

- Length of the land surface: `L0`
- Width of the inlet channel: `N0`
- Center of the domain: `CTR`
- Length of cell faces: `dx`
- Depth of inlet channel: `h0`
- Sea level: `H_SL`
- Bedload fraction: `f_bedload`
- Sediment concentration: `C0_percent`
- Characteristic Velocity: `u0`
- If subsidence is enabled:
  - Subsidence start time: `start_subsidence`
  - Subsidence rate: `sigma`


Working with Model Outputs
==========================

The resulting netCDF4 output file can be read using any netCDF4-compatible
library. These libraries range from the
`netCDF4 Python package <https://github.com/Unidata/netcdf4-python>`_ itself,
to higher-level libraries such as
`xarray <https://github.com/pydata/xarray>`_. For deltas, and specifically
*pyDeltaRCM*, there is also a package under development called
`DeltaMetrics <https://github.com/DeltaRCM/DeltaMetrics>`_,
that is being designed to help post-process and analyze *pyDeltaRCM* outputs.


Here, we show how to read the output NetCDF file with Python package ``netCDF4``.

.. code::

   import netCDF4 as nc

   data = nc.Dataset('pyDeltaRCM_output.nc')  # the output file path!

This `data` object is a `Dataset` object that can be sliced the same was as a `numpy` array.
For example, we can slice the final bed elevation and velocity of a model run:

.. code::

   final_bed_elevation = data['eta'][-1, :, :]
   final_velocity = data['velocity'][-1, :, :]

These slices look like this, if we were to plot them.

.. plot:: guides/output_file.py
**************
Morphodynamics
**************

.. currentmodule:: pyDeltaRCM.sed_tools

pyDeltaRCM approximates sediment dispersal through the use of a weighted random
walk dictated by water flux.
In turn, sediment dispersal drives bed elevation change in the model domain by mass conservation.

See [1]_ for a complete description of morphodynamic assumptions in the DeltaRCM model.
In this documentation, we focus on the details of *model implementation*, rather than *model design*.


.. _sediment-transport:

==================
Sediment Transport
==================

Sediment transport in the model is computed according to an excess stress approach.
Conceptually, sand is routed as bed-material load, and mud is routed as fully suspended load.

For sand parcels, the *transport capacity* is determined by the scaling between sediment flux and flow velocity, and takes the form of the Meyer-Peter and Müller (1948) [3]_ formula:

.. math::

      q_{s\_cap} = q_{s0} \frac{u^\beta_{loc}}{u^\beta_0},

where :math:`u_{loc}` is the depth averaged flow velocity in the cell, :math:`beta` is an exponent set to 3 by default (:obj:`~pyDeltaRCM.model.DeltaModel.beta`), and :math:`q_{s0}` is the unit-width upstream sand flux input at the inlet channel.
At each step of the model domain, sand is either eroded or deposited to the bed depending on the local flow velocity :math:`u_{loc}` and local sediment transport :math:`q_{s\_loc}`. 
Sand is deposited where local transport (i.e., the sand put into that cell from upstream) is greater than the cell transport capacity :math:`q_{s\_loc} > q_{s\_cap}`.
Sand is eroded from the bed when the local velocity is greater than the threshold erosion velocity (:obj:`~pyDeltaRCM.model.DeltaModel.coeff_U_ero_sand`) **and** the local transport is less than the local transport capacity.


Mud parcels do not have any local capacity (i.e., fully suspended washload transport). 
At each parcel step, mud is either eroded or deposited (or neither), depending on the relative value of local flow velocity :math:`u_{loc}` and the threshold erosion and deposition values (:obj:`~pyDeltaRCM.model.DeltaModel.coeff_U_ero_mud` and :obj:`~pyDeltaRCM.model.DeltaModel.coeff_U_dep_mud`).

.. note:: 

      A complete conceptual description of sediment erosion and deposition routing rules can be found in the original DeltaRCM reference Liang et al., 2015 [1]_.


.. _sediment-routing-weighting:

Sediment routing weighting
--------------------------

Sediment routing probability for a given cell :math:`j` to neighbor cell :math:`i` is computed according to:

.. math::

    w_i = \frac{\frac{1}{R_i} \max(0, \mathbf{F}\cdot\mathbf{d_i})}{\Delta i},

where :math:`\mathbf{F}` is the local routing direction and :math:`\mathbf{d_i}` is a unit vector pointing to neighbor :math:`i` from cell :math:`j`, and :math:`\Delta_i` is the cellular distance to neighbor :math:`i` (:math:`1` for cells in main compass directions and :math:`\sqrt{2}` for corner cells.
:math:`R_i` is a resistance estimated as an inverse function of local water depth (:math:`h_i`):

.. math::

    R_i = \frac{1}{{h_i}^\theta}.

Here, :math:`\theta` takes the value of :obj:`~pyDeltaRCM.model.DeltaModel.coeff_theta_sand` for sand routing probabilities, and :obj:`~pyDeltaRCM.model.DeltaModel.coeff_theta_mud` for mud routing.


.. plot:: sed_tools/sediment_weights_examples.py


============================
Changes in the bed elevation
============================

Along the walk of a sediment parcel, the sediment parcel volume is modulated on each step, according to the sediment transport rules described above in :ref:`sediment-transport`.
As the volume of the sediment parcel changes, the channel bed elevation at the current parcel location is updated to reflect this volume change (:obj:`~sed_tools.BaseRouter._update_fields`), i.e., the bed is eroded or sediment is deposited on the bed.
The vertical change in the channel bed is dictated by sediment mass conservation (i.e., Exner equation) and is equal to:

.. math::

    \Delta \eta = \Delta V / dx^2

where :math:`\Delta V` is the volume of sediment to be eroded or deposited from the bed at a given cell along the parcel walk.

.. note::

    Total sediment mass is preserved during erosion, but individual categories of sand and mud are not. I.e., it is assumed that there is an infinite supply of sand and/or mud to erode and entrain at any location in the model domain.

Following a change in the bed elevation, the local flow depth is updated and then local flow velocity is updated according to fluid mass conservation (i.e., ``uw = qw / h``; :obj:`~sed_tools.BaseRouter._update_fields`; [1]_).

Sediment parcels are routed through the model domain step-by-step and in serial, such that changes in the bed elevation caused by one sediment parcel will affect the weighted random walks of all subsequent sediment parcels (:ref:`sediment-routing-weighting`), due to the updated flow field.

Sediment parcel routing is handled by first routing all sand parcels, applying a topographic diffusion (see below and :meth:`~sed_tools.topo_diffusion`), and then routing all mud parcels.
The impact of routing *all* sand and mud parcels on bed elevation is shown in the table below.

.. _sand-mud-route-comparison:

.. table::

    +-------------------------------------------+-----------------------------------------------+----------------------------------------------+
    | initial bed                               | :meth:`~sed_tools.route_all_sand_parcels`     | :meth:`~sed_tools.route_all_mud_parcels`     |
    +===========================================+===============================================+==============================================+
    | .. plot:: sed_tools/_initial_bed_state.py | .. plot:: sed_tools/route_all_sand_parcels.py | .. plot:: sed_tools/route_all_mud_parcels.py |
    +-------------------------------------------+-----------------------------------------------+----------------------------------------------+

.. _model-stability:

===============
Model Stability
===============

Model stability depends on a number of conditions.
At its core though, model stability depends on the bed elevation rate of change over space and time. 
Rapid and/or abrupt bed elevation change trigger numerical instability that can *occasionally* run-away and cause model runs to fail.
A number of processes are included in the DeltaRCM framework to help limit the possibility of failed runs.


.. _reference-volume:

Reference Volume
----------------

The reference volume if the foundational unit (:math:`V_0`) impacting model stability.
This value characterizes the volume on one inlet-channel cell, from the channel bed to the water surface:

.. math::

    V_0 = h_0 {\delta_c}^2

where :math:`h_0` is the inlet channel depth (meters) and :math:`\delta_c` is the cell length (meters).


Time stepping
-------------

Perhaps most important to model stability is the model timestep. 
Recall that for each iteration, the number of parcels and input sediment discharge (as `h0 u0 (c0_percent)/100`) are set by the user (fixed).
Therefore, to control stability, the duration of time represented by each iteration (i.e., the timestep) is determined such that changes in bed elevation per iteration are small.
The model timestep is determined as:

.. math::

    dt = dV_s / N_{p,sed}

where :math:`N_{p,sed}` is the number of sediment parcels, :math:`dV_s` is a characteristic sediment volume, based on the reference volume and inlet width as :math:`dV_s = 0.1 N_0^2 V_0`, where :math:`N_0` is the number of cells across the inlet.


Limiting bed elevation change
-----------------------------

At each sediment parcel step, bed elevation change is limited to 1/4 of the local flow depth. 
Additionally, an edge case where repeated channel bed deposition creates a local `depth` < 0 is restricted by enforcing zero deposition if the `depth` < 0.
These regulations are implemented in the `BaseRouter` class, as :obj:`~pyDeltaRCM.sed_tools.BaseRouter._limit_Vp_change`.


.. topographic-diffusion:

Topographic diffusion
---------------------

Abrupt change in bed elevation (i.e., steep local bed slope) may lead to numerical instability. 
To prevent this, a topographic diffusion is applied immediately following the routing of all sand parcels in the model sequence.

.. hint::

    Topographic diffusion is applied between routing sand parcels and routing mud parcels.

In implementation, topographic smoothing convolves topography with `3x3` cell kernels configured to a diffusive behavior.
The diffusion is repeated over the entire model domain :obj:`~pyDeltaRCM.DeltaModel.N_crossdiff` times.
In the following example, :obj:`~pyDeltaRCM.DeltaModel.N_crossdiff` takes the :doc:`default value </reference/model/yaml_defaults>`.

.. plot:: sed_tools/topo_diffusion.py

The impact of topographic diffusion is minor compared to the bed elevation change driven by parcel erosion or deposition (:ref:`sand and mud routing effects <sand-mud-route-comparison>`).


Notes for modeling best practices
=================================

* Stop simulations before the delta reaches the edge of the computational domain. Delta channel dynamics are changed when a distributary channel reaches the domain edge, because sediment is conveyed to outside the computational domain where it no longer feeds back on channel development. Channels become "locked" in place in this scenario [2]_, because the domain edge is an infinite sediment sink, and therefore rendering invalid any assumptions about stationarity of delta dynamics and/or stratigraphy. Moreover, the downstream water surface boundary condition (`H_SL`) will be violated if a channel reaches the domain edge. Generally, simulations that reach the edge of the domain should be discarded.
* Stop simulations before the delta reaches a condition where the topset slope is equal to the background slope parameter (`S0`). When the background slope is reached, the transport capacity of sediment through the delta is diminished such that channels "clog" up and trigger model instabilities. This is really only an issue for large domains run for long duration.
* Use a sufficient number of water and sediment parcels (> 2000). Too few parcels will result in a rough water surface, irregular sediment deposition, and a rough bed elevation.



References
==========

.. [1] A reduced-complexity model for river delta formation – Part 1: Modeling
       deltas with channel dynamics, M. Liang, V. R. Voller, and C. Paola, Earth
       Surf. Dynam., 3, 67–86, 2015. https://doi.org/10.5194/esurf-3-67-2015

.. [2] Liang, M., Kim, W., and Passalacqua, P. (2016), How much subsidence is
       enough to change the morphology of river deltas?, Geophysical Research Letters, 43, 10,266--10,276, doi:10.1002/2016GL070519.

.. [3] Meyer-Peter, E. and Müller, R.: Formulas for bed-load transport, in: 
       Proceedings of the 2nd Meeting of IAHSR, Stockholm, Sweden, 39–64, 1948.***************
YAML Parameters
***************

Configurable model parameters are listed in the
:doc:`../reference/model/yaml_defaults`.
Links to the API documentation for different types of YAML Parameters are
provided below. The YAML parameters are sorted by "type", for example,
:ref:`model-domain-parameters` are those parameters which control the
definition of the pyDeltaRCM model domain.

Model Settings
==============

:attr:`pyDeltaRCM.model.DeltaModel.out_dir`

:attr:`pyDeltaRCM.model.DeltaModel.verbose`

:attr:`pyDeltaRCM.model.DeltaModel.seed`


.. _model-domain-parameters:

Model Domain Parameters
=======================

:attr:`pyDeltaRCM.model.DeltaModel.Length`

:attr:`pyDeltaRCM.model.DeltaModel.Width`

:attr:`pyDeltaRCM.model.DeltaModel.dx`

:attr:`pyDeltaRCM.model.DeltaModel.L0_meters`

:attr:`pyDeltaRCM.model.DeltaModel.S0`

:attr:`pyDeltaRCM.model.DeltaModel.itermax`

:attr:`pyDeltaRCM.model.DeltaModel.Np_water`

:attr:`pyDeltaRCM.model.DeltaModel.u0`

:attr:`pyDeltaRCM.model.DeltaModel.N0_meters`

:attr:`pyDeltaRCM.model.DeltaModel.h0`

:attr:`pyDeltaRCM.model.DeltaModel.hb`

:attr:`pyDeltaRCM.model.DeltaModel.H_SL`

:attr:`pyDeltaRCM.model.DeltaModel.SLR`

:attr:`pyDeltaRCM.model.DeltaModel.Np_sed`

:attr:`pyDeltaRCM.model.DeltaModel.f_bedload`

:attr:`pyDeltaRCM.model.DeltaModel.active_layer_thickness`

:attr:`pyDeltaRCM.model.DeltaModel.C0_percent`

:attr:`pyDeltaRCM.model.DeltaModel.Csmooth`

:attr:`pyDeltaRCM.model.DeltaModel.toggle_subsidence`

:attr:`pyDeltaRCM.model.DeltaModel.subsidence_rate`

:attr:`pyDeltaRCM.model.DeltaModel.start_subsidence`


Output Settings
===============

:attr:`pyDeltaRCM.model.DeltaModel.save_eta_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_stage_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_depth_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_discharge_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_velocity_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_sedflux_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_sandfrac_figs`

:attr:`pyDeltaRCM.model.DeltaModel.save_figs_sequential`

:attr:`pyDeltaRCM.model.DeltaModel.save_eta_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_stage_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_depth_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_discharge_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_velocity_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_sedflux_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_sandfrac_grids`

:attr:`pyDeltaRCM.model.DeltaModel.save_discharge_components`

:attr:`pyDeltaRCM.model.DeltaModel.save_velocity_components`

:attr:`pyDeltaRCM.model.DeltaModel.save_dt`

:attr:`pyDeltaRCM.model.DeltaModel.checkpoint_dt`

:attr:`pyDeltaRCM.model.DeltaModel.save_checkpoint`

:attr:`pyDeltaRCM.model.DeltaModel.resume_checkpoint`

:attr:`pyDeltaRCM.model.DeltaModel.clobber_netcdf`

:attr:`pyDeltaRCM.model.DeltaModel.legacy_netcdf`


Reduced-Complexity Routing Parameters
=====================================

:attr:`pyDeltaRCM.model.DeltaModel.omega_sfc`

:attr:`pyDeltaRCM.model.DeltaModel.omega_flow`

:attr:`pyDeltaRCM.model.DeltaModel.Nsmooth`

:attr:`pyDeltaRCM.model.DeltaModel.theta_water`

:attr:`pyDeltaRCM.model.DeltaModel.coeff_theta_sand`

:attr:`pyDeltaRCM.model.DeltaModel.coeff_theta_mud`

:attr:`pyDeltaRCM.model.DeltaModel.beta`

:attr:`pyDeltaRCM.model.DeltaModel.sed_lag`

:attr:`pyDeltaRCM.model.DeltaModel.coeff_U_dep_mud`

:attr:`pyDeltaRCM.model.DeltaModel.coeff_U_ero_mud`

:attr:`pyDeltaRCM.model.DeltaModel.coeff_U_ero_sand`

:attr:`pyDeltaRCM.model.DeltaModel.alpha`

:attr:`pyDeltaRCM.model.DeltaModel.stepmax`
.. model_info:

=================
Model Information
=================

This section provides background information about the pyDeltaRCM model itself.
What are the model variables and parameters?
What do they mean?
How does the model simulate hydrodynamics?
How is morphodynamics handled?
All of these questions and more are answered here!

.. toctree::
   :maxdepth: 1

   yamlparameters
   modeltime
   hydrodynamics
   morphodynamics
   outputfile
*************
Hydrodynamics
*************

.. currentmodule:: pyDeltaRCM.water_tools

pyDeltaRCM approximates hydrodynamics through the use of a weighted random walk.
See [1]_ and [2]_ for a complete description of hydrodynamic assumptions in the DeltaRCM model.
In this documentation, we focus on the details of *model implementation*, rather than *model design*.

The water routing operations are orchestrated by :obj:`~water_tools.route_water`; the following sections narrate the sequence of events within a single call of this method.


Routing individual water parcels
================================

Probabilities for water parcel routing *to all neighbors and to self* for *each cell* are computed *once* at the beginning of the water routing routine (:obj:`~water_tools.get_water_weight_array` called from :obj:`~water_tools.run_water_iteration`).

Water routing probability for a given cell :math:`j` to neighbor cell :math:`i` is computed according to:

.. math::

    w_i = \frac{\frac{1}{R_i} \max(0, \mathbf{F}\cdot\mathbf{d_i})}{\Delta i},

where :math:`\mathbf{F}` is the local routing direction and :math:`\mathbf{d_i}` is a unit vector pointing to neighbor :math:`i` from cell :math:`j`, and :math:`\Delta_i` is the cellular distance to neighbor :math:`i` (:math:`1` for cells in main compass directions and :math:`\sqrt{2}` for corner cells.
:math:`R_i` is a flow resistance estimated as an inverse function of local water depth (:math:`h_i`):

.. math::

    R_i = \frac{1}{{h_i}^\theta}

The exponent :math:`\theta` takes a value of unity by default for water routing (:attr:`~pyDeltaRCM.DeltaModel.theta_water`, :doc:`/reference/model/yaml_defaults`), leading to routing weight for neighbor cell :math:`i`:

.. math::

    w_i = \frac{h_i \max(0, \mathbf{F}\cdot\mathbf{d_i})}{\Delta i},

These weights above are calculated only for wet neighbor cells; all dry neighbor cells take a weight value of 0 (:obj:`~water_tools._get_weight_at_cell_water`).
Finally, probability for routing from cell :math:`j` to cell :math:`i` is calculated as:

.. math::

    p_i = \frac{w_i}{\sum^8_{nb=1} w_{nb}}, i=1, 2, \ldots, 8

Weights are accumulated for 8 neighbors and a probability of 0 is assigned to moving from cell :math:`j` to cell :math:`j` (i.e., no movement).
These 9 probabilities are organized into an array ``self.water_weights`` with shape (:obj:`L`, :obj:`W`, 9)`. 

The following figure shows several examples of locations within the model domain, and the corresponding water routing weights determined for that location.

.. plot:: water_tools/water_weights_examples.py

Because probabilities are computed for all locations once at the beginning of water iteration, all water parcels can be routed *in parallel* step-by-step in :obj:`~water_tools.run_water_iteration`.
During iteration, the direction of the random walk is chosen for each parcel via :obj:`_choose_next_directions`, which internally uses :func:`~pyDeltaRCM.shared_tools.random_pick` for randomness.
For example, see the random walks of several parcels below:

.. plot::  water_tools/run_water_iteration.py


.. todo:: add sentence or two above about check_for_loops.

Water routing completes when all water parcels have either 1) reached the model domain boundary, 2) taken a number of steps exceeding :attr:`~pyDeltaRCM.model.DeltaModel.stepmax`, or 3) been removed from further routing via the :obj:`_check_for_loops` function.



Combining parcels into free surface
===================================

Following the routing of water parcels, these walks must be converted in some meaningful way to a model field representing a free surface (i.e., the water stage).
First, the :meth:`~water_tools.compute_free_surface` is called, which takes as input the current bed elevation, and the path of each water parcel (top row in figure below).

.. plot:: water_tools/compute_free_surface_inputs.py

The :meth:`~water_tools.compute_free_surface` method internally calls the :func:`_accumulate_free_surface_walks` function to determine 1) the number of times each cell has been visited by a water parcel
(``sfc_visit``), and 2) the *total sum of expected elevations* of the water surface at each cell (``sfc_sum``).
:func:`_accumulate_free_surface_walks` itself iterates through each water parcel, beginning from the end-point of the path, and working upstream; note that parcels that have been determined to "loop" (:func:`_check_for_loops` and described above) are excluded from computation in determining the free surface.
While downstream of the land-ocean boundary (determined by a depth-or-velocity threshold), the water surface elevation is assumed to be ``0``, whereas upstream of this boundary, the predicted elevation of the water surface is determined by the distance from the previously identified water surface elevation and the background land slope (:attr:`~pyDeltaRCM.DeltaModel.S0`), such the the water surface maintains an approximately constant slope for each parcel pathway.

.. plot:: water_tools/_accumulate_free_surface_walks.py

The algorithm tracks the number of times each cell has been visited by a water parcel (``sfc_visit``), and the *total sum of expected elevations* of the water surface at each cell (``sfc_sum``), by adding the predicted surface elevation of each parcel step while iterating through each step of each parcel.

Next, the output from :func:`_accumulate_free_surface_walks` is used to calculate a new stage surface (``H_new``) based only on the water parcel paths and expected water surface elevations, approximately as ``H_new = sfc_sum / sfc_visit``.
The updated water surface is combined with the previous timestep's water surface and an under-relaxation coefficient (:attr:`~pyDeltaRCM.DeltaModel.omega_sfc`).

.. plot:: water_tools/compute_free_surface_outputs.py

With a new free surface computed, a few final operations prepare the surface for boundary condition updates and eventually being passed to the sediment routing operations (inside :meth:`~water_tools.finalize_free_surface`).
A non-linear smoothing operation is applied to the free surface, whereby wet cells are iteratively averaged with neighboring wet cells to yield an overall smoother surface.
The smoothing is handled by :func:`_smooth_free_surface` and depends on the number of iterations (:attr:`~pyDeltaRCM.model.DeltaModel.Nsmooth`) and a weighting coefficient (:attr:`~pyDeltaRCM.model.DeltaModel.Csmooth`).

.. plot:: water_tools/_smooth_free_surface.py

Finally, a :meth:`~water_tools.flooding_correction` is applied to the domain.
In this correction, all "dry" cells (a cell where the flow depth is less than the `dry_depth`) are checked for any neighboring cells where the water surface elevation (`stage`) is higher than the bed elevation of the dry cell.
If this condition is met for a given dry cell, the dry cell is flooded: the stage of the dry cell is set to the maximum stage of neighboring cells.

.. plot:: water_tools/flooding_correction.py

Similar to :func:`~_smooth_free_surface` described above, :meth:`~water_tools.flooding_correction` acts to remove roughness in the water surface and contributes to :ref:`model stability <model-stability>`.


Finalizing and boundary conditions to sediment routing
======================================================

The final step of the water routing operations in each call to :obj:`~water_tools.route_water` is the updating of model fields, and application of a few corrections and boundary conditions.
These operations are handled within :obj:`~water_tools.finalize_water_iteration`.

First, the `stage` field is limited to elevations above or equal to sea level `H_SL`.
Then, the `depth` field is updated as the vertical distance between the `stage` field and the bed elevation `eta`.
These operations ensure that these fields are in agreement over the entire model domain.

Next, methods :obj:`~water_tools.update_flow_field` and :obj:`~water_tools.update_velocity_field` are called, which handle the updating of water discharge and water velocity fields, respectively.


.. note::
  
  The next step in the model `update` sequence is sediment routing (:obj:`~pyDeltaRCM.sed_tools.sed_tools.route_sediment`). For more information on the next stage, see :doc:`morphodynamics`.


References
==========

.. [1] A reduced-complexity model for river delta formation – Part 1: Modeling
       deltas with channel dynamics, M. Liang, V. R. Voller, and C. Paola, Earth
       Surf. Dynam., 3, 67–86, 2015. https://doi.org/10.5194/esurf-3-67-2015

.. [2] A reduced-complexity model for river delta formation – Part 2:
       Assessment of the flow routing scheme, M. Liang, N. Geleynse,
       D. A. Edmonds, and P. Passalacqua, Earth Surf. Dynam., 3, 87–104, 2015.
       https://doi.org/10.5194/esurf-3-87-2015
.. api:

=============
API reference
=============

The :obj:`DeltaModel` class, defined in `pyDeltaRCM.model`, is the main class of pyDeltaRCM, which provides the object that is manipulated to evolve the numerical delta model.
This class uses "mix-in" classes, which are defined in separate files (Python `modules`), that break out logically based on the various stages of model use, and components of the model iteration sequence.
Most model functionality is organized into the various mix-in classes, that are then inherited by the `DeltaModel`. 
Additionally, several routines of the model are organized into module-level functions, which are "jitted" via the `numba <https://numba.pydata.org/>`_  code optimizer library for Python.

This index lists the `pyDeltaRCM` organization, hopefully providing enough information to begin to determine where various components of the model are implemented.
The index includes model classes, methods, and attributes, as well as additionally utility classes and functions.

.. toctree::
   :maxdepth: 2

   model/index
   preprocessor/index
   iteration_tools/index
   init_tools/index
   water_tools/index
   sed_tools/index
   shared_tools/index
   hook_tools/index
   debug_tools/index


References
----------

* :ref:`modindex`
* :ref:`search`


Search the Index
==================

* :ref:`genindex`
.. api.sed_tools:

*********************************
sed_tools
*********************************

.. currentmodule:: pyDeltaRCM.sed_tools

.. todo:: add paragraph description of the module


The tools are defined in ``pyDeltaRCM.sed_tools``. 


Public API methods attached to model
------------------------------------

The following methods are defined in the ``sed_tools`` class, of the ``pyDeltaRCM.sed_tools`` module. 

.. autosummary::

    sed_tools

.. autoclass:: sed_tools


Router classes
--------------

The following classes are defined in the ``pyDeltaRCM.sed_tools`` module. These sediment routing classes are jitted for speed.

.. autosummary:: 
    :toctree: ../../_autosummary

    SandRouter
        :members:
        :inherited-members:
        :private-members:
    
    MudRouter
        :members:
        :inherited-members:
        :private-members:
    
    BaseRouter
        :members:
        :inherited-members:
        :private-members:

sed_tools helper functions
----------------------------

Additionally, the sediment parcel step-weighting function is defined at the module level in ``pyDeltaRCM.sed_tools``.

.. autofunction:: _get_weight_at_cell_sediment
.. api.preprocessor:

*********************************
Preprocessor
*********************************

.. currentmodule:: pyDeltaRCM.preprocessor

The high-level API is principally defined in ``pyDeltaRCM.preprocessor``. 

.. todo::

    add paragraph description of the module. What can we do with this? Why does it exist? Link to the relevant documentation in the user guide and examples.


Preprocessor classes and API
----------------------------

The following classes are defined in the ``pyDeltaRCM.preprocessor`` module and enable the high-level model API to work at both the command line and as a Python object.

.. autosummary:: 
    :toctree: ../../_autosummary

    PreprocessorCLI
        :members:
        :inherited-members:

    Preprocessor
        :members:
        :inherited-members:
    
    BasePreprocessor

Preprocessor function and utilities
-----------------------------------

.. todo:: add description, what are these?

.. autofunction:: preprocessor_wrapper
.. autofunction:: scale_relative_sea_level_rise_rate
.. autofunction:: _write_yaml_config_to_file
.. api.iteration_tools:

***************
iteration_tools
***************

.. currentmodule:: pyDeltaRCM.iteration_tools

.. todo::

    Add paragraph description of the module. What stages are defined here generally? Make a table with the main ones like in water tools?


Public API methods attached to model
------------------------------------

The following methods are defined in the ``iteration_tools`` class, of the ``pyDeltaRCM.iteration_tools`` module. 

.. autosummary::

    iteration_tools

.. autoclass:: iteration_tools
.. api.debug_tools:

***********
debug_tools
***********

.. currentmodule:: pyDeltaRCM.debug_tools

The debugging tools are defined in ``pyDeltaRCM.debug_tools``. 

.. todo::

    Add paragraph description of the module. What stages are defined here generally?


Public API methods attached to model
------------------------------------

The following methods are defined in the ``debug_tools`` class, of the ``pyDeltaRCM.debug_tools`` module.
They are then attached as methods of the `DeltaModel` and can be called at any time during run.

.. autosummary::

    debug_tools

.. autoclass:: debug_tools


Public plotting methods
-----------------------

The functions defined below are (generally) the actual workers that handle the plotting for the methods defined above and attached to model.
We expose these functions here because they are be useful in the documentation, where they are used extensively.

.. autofunction:: plot_domain
.. autofunction:: plot_ind
.. autofunction:: plot_line
.. api.shared_tools:

*********************************
shared_tools
*********************************

.. currentmodule:: pyDeltaRCM.shared_tools


.. todo:: add paragraph description of the module


The tools are defined in ``pyDeltaRCM.shared_tools``.


Shared functions
----------------

This module defines several functions that are used throughout the model, and so are organized here for convenience.

.. autofunction:: get_random_uniform
.. autofunction:: get_start_indices
.. autofunction:: get_steps
.. autofunction:: random_pick
.. autofunction:: custom_unravel
.. autofunction:: custom_ravel
.. autofunction:: custom_pad
.. autofunction:: get_weight_sfc_int
.. autofunction:: custom_yaml_loader


Time scaling functions
----------------------

Scaling of real-world time and model time is an important topic covered in detail in :doc:`/info/modeltime`.
Several functions are defined here which can help with scaling between model and real-world time.

.. autofunction:: scale_model_time
.. autofunction:: _scale_factor


Utilities
---------

Additionally, functions defined in ``pyDeltaRCM.shared_tools`` manage the random state of the model, and help with documentation and version management.

.. autofunction:: set_random_seed
.. autofunction:: get_random_state
.. autofunction:: set_random_state
.. autofunction:: _docs_temp_directory
.. autofunction:: _get_version
.. api.init_tools:

**********
init_tools
**********

.. currentmodule:: pyDeltaRCM.init_tools


The model initialization is managed by :obj:`~pyDeltaRCM.model.DeltaModel.__init__`, but the actual initialization is mostly handled by routines in `init_tools`.
The major steps of initialization are:

.. autosummary::

    init_tools.import_files
    init_tools.init_logger
    init_tools.process_input_to_model
    init_tools.determine_random_seed
    init_tools.create_other_variables
    init_tools.create_domain
    init_tools.init_sediment_routers
    init_tools.init_subsidence

and then depending on the checkpointing configuration, the following methods may be called:

.. autosummary::

    init_tools.load_checkpoint
    init_tools.init_output_file
    pyDeltaRCM.iteration_tools.iteration_tools.output_data
    pyDeltaRCM.iteration_tools.iteration_tools.output_checkpoint
    pyDeltaRCM.iteration_tools.iteration_tools.log_model_time


Public API methods attached to model
------------------------------------

The following methods are defined in the ``init_tools`` class, of the ``pyDeltaRCM.init_tools`` module. 

.. currentmodule:: pyDeltaRCM.init_tools

.. autosummary::

    init_tools

.. autoclass:: init_tools
*********************************
hook_tools
*********************************

.. currentmodule:: pyDeltaRCM.hook_tools

.. todo::

    Add paragraph description of the module. What stages are defined here generally? Link to relevant documentation about hooks in user guide and model list.


Public API methods attached to model
------------------------------------

The following methods are defined in the ``hook_tools`` class, of the ``pyDeltaRCM.hook_tools`` module. 

.. autosummary::

    hook_tools

.. autoclass:: hook_tools
Available Model Hooks
=====================

A complete list of hooks in the model follows:

.. currentmodule:: pyDeltaRCM.hook_tools

.. csv-table:: Available model hooks
    :header: "Initializing", "Updating"
    :widths: 40, 40

    :obj:`~hook_tools.hook_import_files`, :obj:`~hook_tools.hook_solve_water_and_sediment_timestep`
    :obj:`~hook_tools.hook_process_input_to_model`, :obj:`~hook_tools.hook_apply_subsidence`
    :obj:`~hook_tools.hook_create_other_variables`, :obj:`~hook_tools.hook_finalize_timestep`
    :obj:`~hook_tools.hook_create_domain`, :obj:`~hook_tools.hook_route_water`
    :obj:`~hook_tools.hook_load_checkpoint`, :obj:`~hook_tools.hook_init_water_iteration`
    :obj:`~hook_tools.hook_output_data`, :obj:`~hook_tools.hook_run_water_iteration`
    :obj:`~hook_tools.hook_output_checkpoint`, :obj:`~hook_tools.hook_compute_free_surface`
    :obj:`~hook_tools.hook_init_output_file`, :obj:`~hook_tools.hook_finalize_water_iteration`
    , :obj:`~hook_tools.hook_route_sediment`
    , :obj:`~hook_tools.hook_route_all_sand_parcels`
    , :obj:`~hook_tools.hook_topo_diffusion`
    , :obj:`~hook_tools.hook_route_all_mud_parcels`
    , :obj:`~hook_tools.hook_compute_sand_frac`
.. api.model:

*****
model
*****

This is the main model code.

This class is defined in ``pyDeltaRCM.model``. 


The DeltaModel class and all attributes and methods
---------------------------------------------------

.. currentmodule:: pyDeltaRCM.model

.. autosummary:: 
	:toctree: ../../_autosummary

	~DeltaModel




Model configuration reference
-----------------------------

.. toctree::
    :maxdepth: 1

    yaml_defaults
    model_hooks
Default Model Variable Values
=============================

Default model values are defined as:

.. literalinclude:: ../../../../pyDeltaRCM/default.yml
   :language: yaml
   :linenos: 
.. api.water_tools:

*********************************
water_tools
*********************************

.. currentmodule:: pyDeltaRCM.water_tools


The :obj:`~pyDeltaRCM.water_tools.water_tools.route_water` routine manages the water routing.
During :obj:`~pyDeltaRCM.water_tools.water_tools.route_water`, water iteration is repeated a total of :obj:`~pyDeltaRCM.model.DeltaModel.itermax` times.
During each of these iterations of the water routing, the following methods are called *in order*:

.. autosummary::

    water_tools.init_water_iteration
    water_tools.run_water_iteration
    water_tools.compute_free_surface
	water_tools.finalize_water_iteration


Public API methods attached to model
------------------------------------

The following methods are defined in the ``water_tools`` class, of the ``pyDeltaRCM.water_tools`` module.

.. autosummary::

    water_tools

.. autoclass:: water_tools


water_tools helper functions
----------------------------

The following routines are jitted for speed.
They generally take a large number of arrays and constants and return a new array(s) to continue with the model progression within the methods defined above.

.. autofunction:: _get_weight_at_cell_water
.. autofunction:: _choose_next_directions
.. autofunction:: _calculate_new_inds
.. autofunction:: _check_for_loops
.. autofunction:: _update_dirQfield
.. autofunction:: _update_absQfield
.. autofunction:: _accumulate_free_surface_walks
.. autofunction:: _smooth_free_surface
Time-varying bedload
====================

The following example demonstrates how to configure a subclassing model with a time-varying parameter.
In this example, the time-varying behavior arises by managing two "switches", ``self._changed`` and ``self._changed_back``, which change the state of the :obj:`f_bedload` parameters at a predetermined time in the mode sequence.

The following codes produce two runs (using `matrix` expansion from the Preprocessor), which has a baseline `f_bedload` value of either ``0.3`` or ``0.7``, and for a period in the middle of the run, the `f_bedload` values are mirrored by ``self.f_bedload = (1 - self.f_bedload)``, i.e., briefly switching the bedload values.


.. plot::
    :context: reset
    :include-source:

    class VariableBedloadModel(pyDeltaRCM.DeltaModel):

        def __init__(self, input_file=None, **kwargs):

            super().__init__(input_file, **kwargs)

            self._changed = False
            self._changed_back = False

        def hook_solve_water_and_sediment_timestep(self):
            """Change the state depending on the _time.
            """
            # check if the state has been changed, and time to change it
            if (not self._changed) and (self._time > 459909090):
                self.f_bedload = (1 - self.f_bedload)
                self._changed = True

                _msg = 'Bedload changed to {f_bedload}'.format(
                    f_bedload=self.f_bedload)
                self.log_info(_msg, verbosity=0)

            # check if the state has been changed back, and time to change it back
            if (not self._changed_back) and (self._time > 714545454):
                self.f_bedload = (1 - self.f_bedload)
                self._changed_back = True

                _msg = 'Bedload changed back to {f_bedload}'.format(
                    f_bedload=self.f_bedload)
                self.log_info(_msg, verbosity=0)


.. rubric:: Checking on the state-change effect

To demonstrate how this works, let's loop through time and check the model state.
Here, we will change the model time value directly, so that we can verify that the model is working as intended, but you should never do this in practice.

.. plot::
    :context:

    with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
        mdl_muddy = VariableBedloadModel(f_bedload=0.3,
                                         out_dir=output_dir)
        mdl_sandy = VariableBedloadModel(f_bedload=0.7)

.. code:: python

    mdl_muddy = VariableBedloadModel(f_bedload=0.3,
                                     out_dir=output_dir)
    mdl_sandy = VariableBedloadModel(f_bedload=0.7)

.. important::

    You should never modify the model time via ``self._time`` directly when working with the model.


.. plot::
    :context:
    :include-source:

    # create a figure
    fig, ax = plt.subplots()

    # set the "simulation" range
    _times = np.linspace(0, 1e9, num=100)
    fb_mdl_muddy = np.zeros_like(_times)
    fb_mdl_sandy = np.zeros_like(_times)

    # loop through time, change the model time and grab f_bedload values
    for i, _time in enumerate(_times):
        # change the model time directly
        mdl_muddy._time = _time  # you should never do this
        mdl_sandy._time = _time  # you should never do this

        # run the hooked method
        mdl_muddy.hook_solve_water_and_sediment_timestep()
        mdl_sandy.hook_solve_water_and_sediment_timestep()

        # grab the state of the `f_bedload` parameter
        fb_mdl_muddy[i] = mdl_muddy.f_bedload  # get the value
        fb_mdl_sandy[i] = mdl_sandy.f_bedload  # get the value

    # add it to the plot
    ax.plot(_times, fb_mdl_muddy, '-', c='saddlebrown', lw=2, label='muddy')
    ax.plot(_times, fb_mdl_sandy, '--', c='goldenrod', lw=2, label='sandy')
    ax.legend()

    # clean up
    ax.set_ylim(0, 1)
    ax.set_ylabel('f_bedload')
    ax.set_xlabel('model time (s)')

    plt.show()


.. rubric:: Running the model for real

Given a yaml file (``variable_bedload.yaml``):

.. code:: yaml

    Length: 5000.
    Width: 10000.
    dx: 50
    N0_meters: 500
    C0_percent: 0.05
    SLR: 1.5e-8
    h0: 2
    u0: 1.1
    coeff_U_dep_mud: 0.5
    parallel: True
    
    matrix:
      f_bedload:
        - 0.3
        - 0.7


and a script to run the code:

.. code:: python
    
    if __name__ == '__main__':

        # base yaml configuration
        base_yaml = 'variable_bedload.yaml'

        pp = pyDeltaRCM.Preprocessor(
                base_yaml,
                timesteps=12000)

        # run the jobs
        pp.run_jobs(DeltaModel=VariableBedloadModel)

Saving custom fields to the output file
=======================================

Given the flexibility of the pyDeltaRCM model to modifications via hooks and subclassing, it is necessary that the variables saved to the output netCDF file are similarly customizable.
Fortunately, the use of subclasses and hooks itself enables flexible setting of gridded variables to be saved as figures, as well as customization of the fields saved to the netCDF file as both variables, and as metadata.

To customize the figures and variables that are saved, the hook ``hook_init_output_file`` should be used to subclass the pyDeltaRCM model.

When adding a model attribute to the key-value pairs of grids to save as figures, the key indicates the name of the figure file that will be saved, and the value-pair can be a string representing the model attribute to be plotted, or a combination of model attributes, such as ``self.eta * self.depth``.
For example, ``self._save_fig_list['active_layer'] = ['active_layer']`` will properly indicate that figures of the active layer should be saved.

.. important::

    The built-in, on-the-fly, figure saving as the model runs is only supported for gridded variables with the shape ``L x W`` matching the model domain.
    Trying to set up figures to save that are not variables of that shape will result in an error.

When adding variables or metadata to be initialized and subsequently saved in the output netCDF, the key-value pair relationship is as follows.
The key added to ``self._save_var_list`` is the name of the variable as it will be recorded in the netCDF file, this *does not* have to correspond to the name of an attribute in the model.
To add a variable to the metadata, a key must be added to ``self._save_var_list['meta']``.
The expected value for a given key is a list containing strings indicating the model attribute to be saved, its units, the variable type, and lastly the variable dimensions (e.g., ``['active_layer', 'fraction', 'f4', ('time', 'x', 'y')]`` for the active layer).

.. important::

    The dimensions of the custom variable being specified must match *exactly* with one of the three standard dimensions: `x`, `y`, `time`.
    Use of an invalid dimension will result in an error.

An example of using the hook and creating a model subclass to customize the figures, gridded variables, and metadata being saved is provided below.

.. doctest::

    >>> import pyDeltaRCM

    >>> class CustomSaveModel(pyDeltaRCM.DeltaModel):
    ...     """A subclass of DeltaModel to save custom figures and variables.
    ...
    ...     This subclass modifies the list of variables and figures used to
    ...     initialize the netCDF file and save figures and grids before the
    ...     output file is setup and initial conditions are plotted.
    ...     """
    ...     def __init__(self, input_file=None, **kwargs):
    ...
    ...         # inherit base DeltaModel methods
    ...         super().__init__(input_file, **kwargs)
    ...
    ...     def hook_init_output_file(self):
    ...         """Add non-standard grids, figures and metadata to be saved."""
    ...         # save a figure of the active layer each save_dt
    ...         self._save_fig_list['active_layer'] = ['active_layer']
    ...
    ...         # save the active layer grid each save_dt w/ a short name
    ...         self._save_var_list['actlay'] = ['active_layer', 'fraction',
    ...                                          'f4', ('time',
    ...                                                 'x', 'y')]
    ...
    ...         # save number of water parcels w/ a long name
    ...         self._save_var_list['meta']['water_parcels'] = ['Np_water',
    ...                                                         'parcels',
    ...                                                         'i8', ()]

Next, we instantiate the model class.

.. code::

    >>> mdl = CustomSaveModel()


.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     mdl = CustomSaveModel(out_dir=output_dir)


This subclass has added the active layer as a figure and a grid to be saved, as well as the number of water parcels as metadata to be saved.
For simplicity we will just check that the appropriate parameters were added to the save figure and save variable lists, however please feel free to give this example a try on your local machine and examine the output figures and netCDF file.

.. doctest::

    >>> 'active_layer' in mdl._save_fig_list
    True

    >>> print(mdl._save_fig_list)
    {'active_layer': ['active_layer']}

    >>> print(mdl._save_var_list)
    {'meta': {'water_parcels': ['Np_water', 'parcels', 'i8', ()]}, 'actlay': ['active_layer', 'fraction', 'f4', ('time', 'x', 'y')]}
Constraining subsidence to part of the domain
=============================================

One case that has been explored in the literature with the DeltaRCM model is the case of subsidence limited to one region of the model domain [1]_.
This model configuration can be readily achieved with model subclassing.

Setting up the custom subclass
------------------------------

.. plot::
    :context: reset
    :include-source:

    class ConstrainedSubsidenceModel(pyDeltaRCM.DeltaModel):
        """A simple subclass of DeltaModel with subsidence region constrained.
    
        This subclass *overwrites* the `init_subsidence` method to
        constrain subsidence to only one region of the model domain.
        """
        def __init__(self, input_file=None, **kwargs):
    
             # inherit base DeltaModel methods
            super().__init__(input_file, **kwargs)

        def init_subsidence(self):
            """Initialize subsidence pattern constrained to a tighter region.

            Uses theta1 and theta2 to set the angular bounds for the
            subsiding region. theta1 and theta2 are set in relation to the
            inlet orientation. The inlet channel is at an angle of 0, if
            theta1 is -pi/3 radians, this means that the angle to the left of
            the inlet that will be included in the subsiding region is 30
            degrees. theta2 defines the right angular bounds for the subsiding
            region in a similar fashion.
            """
            _msg = 'Initializing custom subsidence field'
            self.log_info(_msg, verbosity=1)

            if self._toggle_subsidence:

                theta1 = -(np.pi / 3)
                theta2 = 0

                R1 = 0.3 * self.L  # radial limits (fractions of L)
                R2 = 0.8 * self.L

                Rloc = np.sqrt((self.y - self.L0)**2 + (self.x - self.W / 2.)**2)

                thetaloc = np.zeros((self.L, self.W))
                thetaloc[self.y > self.L0 - 1] = np.arctan(
                    (self.x[self.y > self.L0 - 1] - self.W / 2.)
                    / (self.y[self.y > self.L0 - 1] - self.L0 + 1))
                self.subsidence_mask = ((R1 <= Rloc) & (Rloc <= R2) &
                                        (theta1 <= thetaloc) &
                                        (thetaloc <= theta2))
                self.subsidence_mask[:self.L0, :] = False

                self.sigma = self.subsidence_mask * self.subsidence_rate * self.dt


Now, initialize the model and look at the field.
Note that the colorscale depicts the magnitude of subsidence in the model *per timestep* (`sigma`, which has units meters).

.. code:: python

    mdl = ConstrainedSubsidenceModel(toggle_subsidence=True)

.. plot::
    :context:

    with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
        mdl = ConstrainedSubsidenceModel(toggle_subsidence=True,
                                         out_dir=output_dir)

.. plot::
    :context:
    :include-source:

    fig, ax = plt.subplots()
    mdl.show_attribute('sigma', grid=False)
    plt.show()


Using the custom subclass with the preprocessor
-----------------------------------------------

We can configure a :obj:`Preprocessor` to handle a set of custom runs in conjunction with out custom `pyDeltaRCM` model subclass.
For example, in [1]_, the authors explore the impact of subsidence at various rates: 3 mm/yr, 6 mm/yr, 10 mm/yr, 25 mm/yr, 50 mm/yr, and 100 mm/yr.
We can scale these rates, assuming a model :doc:`intermittency factor </info/modeltime>` of 0.019, representing 7 of 365 days of flooding per year, by using the convenience function :obj:`~pyDeltaRCM.preprocessor.scale_relative_sea_level_rise_rate`:

.. plot::
    :context: close-figs
    :include-source:

    from pyDeltaRCM.preprocessor import scale_relative_sea_level_rise_rate

    subsidence_mmyr = np.array([3, 6, 10, 25, 50, 100])
    subsidence_scaled = scale_relative_sea_level_rise_rate(subsidence_mmyr, If=0.019)

Now, we use :ref:`matrix expansion <matrix_expansion_tag>` to set up the runs with a preprocessor.
For example, in a Python script, following the definition of the subclass above, define a dictionary with a `matrix` key and supply to the `Preprocessor`:

.. plot::
    :context:
    :include-source:

    # add a matrix with subsidence to the dict
    param_dict = {}
    param_dict['matrix'] = {'subsidence_rate': subsidence_scaled}

    # add other configurations
    param_dict.update(
        {'out_dir': 'liang_2016_reproduce',
         'toggle_subsidence': True,
         'parallel': 3})  # we can take advantage of parallel jobs

.. code::

    # create the preprocessor
    pp = pyDeltaRCM.Preprocessor(
        param_dict,
        timesteps=10000)

And finally run the jobs by specifying the model subclass as the class to use when instantiating the jobs with the preprocessor.

.. below, we overwrite the above, to make sure we only run for one timestep
.. plot::
    :context:

    with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
        param_dict['out_dir'] = output_dir
        pp = pyDeltaRCM.Preprocessor(
            param_dict,
            parallel=False,
            timesteps=1)
        pp.run_jobs(DeltaModel=ConstrainedSubsidenceModel)

.. code:: python

    # run the jobs
    pp.run_jobs(DeltaModel=ConstrainedSubsidenceModel)

We can check whether the runs were set up, as expected:

.. plot::
    :context:
    :include-source:

    from matplotlib.colors import Normalize

    fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(10, 4))
    norm = Normalize(vmin=3, vmax=100)

    for i, job in enumerate(pp.job_list):
        # first convert the field to a rate
        subsidence_rate_field = (job.deltamodel.sigma / job.deltamodel.dt)

        # now convert to mm/yr
        subsidence_rate_field = (subsidence_rate_field * 1000 *
            pyDeltaRCM.shared_tools._scale_factor(If=0.019, units='years'))

        # and display
        im = ax.flat[i].imshow(subsidence_rate_field, norm=norm)

    fig.colorbar(im, ax=ax.ravel().tolist())
    plt.show()


.. [1] Liang, M., Kim, W., and Passalacqua, P. (2016), How much subsidence is
   enough to change the morphology of river deltas?, Geophysical Research Letters, 43, 10,266--10,276, doi:10.1002/2016GL070519.
=============================================
Examples of using and working with pyDeltaRCM
=============================================


Modifying initial conditions
----------------------------

.. toctree::
   :maxdepth: 1

   slight_slope


Modifying boundary conditions
-----------------------------

.. toctree::
   :maxdepth: 1

   updating_boundary_conditions
   subsidence_region
   variable_bedload
   variable_velocity


Modifying internal computations
-------------------------------

* none


Modifying behavior of the Preprocessor
--------------------------------------

* none

Modifying Model Input/Output
----------------------------

.. toctree::
   :maxdepth: 1

   custom_yaml
   custom_saving
Updating model boundary conditions
==================================

In implementing custom model subclasses, it is common to want to change boundary conditions throughout the model run (see :doc:`variable_bedload`, :doc:`variable_velocity`).
In some situations, we want to change a single variable, and see the effect of changing *only this variable*, in essence, pushing the model out of a dynamic equilibrium. 
Another possibility is that we would want to change the boundary conditions of the inlet, *but maintain the dynamic equilibrium*.

Let's create a `DeltaModel` class to demonstrate.

.. code::

    mdl = DeltaModel()


.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     mdl = pyDeltaRCM.DeltaModel(out_dir=output_dir)

Now, we can see that by the default settings, after initialization, the model flow velocity is 1.0 m/s, flow depth is 5 m, and so the unit water discharge is 5 m2/s.

.. doctest::

    >>> mdl.u0
    1.0
    >>> mdl.h0
    5.0
    >>> mdl.qw0
    5.0

If after some number of model iterations, we wanted to change the inlet flow velocity to be 2.0 m/s, we could simply set this value directly.

.. doctest::

    >>> mdl.u0 = 2.0

But, now the model has been thrown out of equilibrium, where the unit water discharge no longer matches the product of the flow depth and flow velocity.

.. doctest::

    >>> mdl.u0
    2.0
    >>> mdl.h0
    5.0
    >>> mdl.qw0
    5.0

To remedy this, we need to use the :obj:`~pyDeltaRCM.init_tools.init_tools.create_boundary_conditions` method, which will reinitialize a number of fields, based on the current value of the inlet flow velocity.

.. doctest::

    >>> mdl.create_boundary_conditions()
    >>> mdl.qw0
    10.0

.. important::

    You are responsible for ensuring that boundary conditions are updated in the appropriate manner after changing certain model parameters. **You need to call the method to reinitialize boundary conditions yourself!**
Defining custom YAML parameters
===============================

.. currentmodule:: pyDeltaRCM.hook_tools

Custom subclasses for the standard ``DeltaModel`` may require additional or custom parameters not listed in the :doc:`/reference/model/yaml_defaults`.
By using the hook, :obj:`~hook_tools.hook_import_files`, it is straightforward to define custom YAML parameters along with expected types and default values for your subclassed model.

The following subclass model demonstrates this by defining a custom boolean parameter (could be used to toggle some custom functionality on/off), and a custom numeric parameter (could be required for the custom function).

.. doctest::

    >>> import pyDeltaRCM

    >>> class CustomParamsModel(pyDeltaRCM.DeltaModel):
    ...     """A subclass of DeltaModel with custom YAML parameters.
    ...
    ...     This subclass defines custom YAML parameters, their expected types
    ...     and default values.
    ...     """
    ...     def __init__(self, input_file=None, **kwargs):
    ...
    ...         # inherit base DeltaModel methods
    ...         super().__init__(input_file, **kwargs)
    ...
    ...     def hook_import_files(self):
    ...         """Define the custom YAML parameters."""
    ...         # custom boolean parameter
    ...         self.subclass_parameters['custom_bool'] = {
    ...             'type': 'bool', 'default': False
    ...         }
    ...
    ...         # custom numeric parameter
    ...         self.subclass_parameters['custom_number'] = {
    ...             'type': ['int', 'float'], 'default': 0
    ...         }

If the subclass model is loaded with a YAML configuration file that does not explicitly define these custom parameters, then the default values will be assigned as model attributes.

.. code::

    >>> defaults = CustomParamsModel()

.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     defaults = CustomParamsModel(out_dir=output_dir)

.. doctest::

    >>> print(defaults.custom_bool)
    False

    >>> print(defaults.custom_number)
    0

.. note::

   Since custom YAML parameters have expected types, ``TypeErrors`` are raised if the custom parameter type provided in the YAML does not agree with what is expected as defined in the subclass.


Once the custom parameters have been defined in the subclassed model they can be treated just like the default model parameters, and can be specified in the YAML or as keyword arguments (``**kwargs``).

.. code::

    >>> customized = CustomParamsModel(custom_bool=True, custom_number=15.3)


.. doctest::
    :hide:

    >>> with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
    ...     customized = CustomParamsModel(out_dir=output_dir,
    ...                                    custom_bool=True,
    ...                                    custom_number=15.3)

.. doctest::

    >>> print(customized.custom_bool)
    True

    >>> print(customized.custom_number)
    15.3
Slightly sloping basin
======================

Consider the case where we are a researcher seeking to explore the effects of a receiving basin that is sloped perpendicular to the channel outlet. 
This researcher asks: does this sloped basin cause channels to steer towards the deeper water, where compensation is higher?

The researcher can easily use subclassing and model hooks to achieve the desired effect.
Recall that anything added to the end of the the subclass' `__init__` method will be called during instantiation of the subclass.

.. plot::
    :context: reset
    :include-source:

    class SlightSlopeModel(pyDeltaRCM.DeltaModel):
        """A subclass of DeltaModel with sloping basin.
    
        This subclass simply modifies the basin geometry
        before any computation has occurred.
        """
        def __init__(self, input_file=None, **kwargs):
    
             # inherit base DeltaModel methods
            super().__init__(input_file, **kwargs)

             # modify the basin
            slope = 0.0005  # cross basin slope
            eta_line = slope * np.arange(0, self.Width,
                                          step=self.dx)
            eta_grid = np.tile(eta_line, (self.L - self.L0, 1))
            eta_grid = eta_grid - ((slope * self.Width)/2)  # center at inlet
            self.eta[self.L0:, :] += eta_grid

Next, we instantiate the model class.

.. code::

    mdl = SlightSlopeModel()


.. plot::
    :context:

    with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
        mdl = SlightSlopeModel(out_dir=output_dir)


And finally, make a plot of the initial condition using the :obj:`~pyDeltaRCM.debug_tools.debug_tools.show_attribute` method.

.. plot::
    :context:
    :include-source:

    fig, ax = plt.subplots()
    mdl.show_attribute('eta', grid=False)
    plt.show()

You can try this out for yourself, and even complete the model run.
Are the channels steered by the basin slope?

.. important:: 

    In this example, we did not take care to update the model `stage` or `depth` fields. In this simple case it works out fine, because after a single timestep, the fields are correctly computed relative to the modified bed. However, take caution when modifying `DeltaModel` fields directly, and be sure to change *all* relevant fields too.
Time-varying inlet velocity
===========================

It is straightforward to create models with time-varying boundary conditions. 
In this example, we set up an array of values to use for the inlet flow velocity, depending on the time of the model (i.e., a timeseries).

In implementation, because we don't know the timestep of the model before creating it, we interpolate the actual inlet flow velocity from the boundary condition timeseries.


Define the boundary condition
-----------------------------

First, let's set up the boundary condition array we want to use.
We'll set this up in a function, so that it can be called from the model subclass, as well as here for plotting.

.. plot::
    :context: reset
    :include-source:

    def create_velocity_array(end_time, a=1, b=5e4, h=3.2, k=2):
        """Create velocity timeseries.
        """
        _time = np.linspace(0, end_time, num=1000)

        _velocity = a * np.sin((_time - h)/b) + k
        return _time, _velocity

The function we define takes the `end_time` of the model run, as well as some shape parameters, and computes according to:

.. math::

    y = a \sin \left( \frac{ \left( t - h \right) }{ b } \right) + k

where :math:`t` is the time array in seconds.
We can inspect the result of this function, as it will be called in the model subclass below.

.. plot::
    :context:
    :include-source:

    end_time = 86400 * 100
    _time_array, _velocity_array = create_velocity_array(
                end_time)

    # make a plot of the boundary condition
    fig, ax = plt.subplots()
    ax.plot(_time_array, _velocity_array)
    ax.set_xlabel('time (seconds)')
    ax.set_ylabel('inlet velocity (meters/second)')
    plt.show()


Define the model subclass
-------------------------

We define a model subclass to handle the changing boundary condition:

.. plot::
    :context: close-figs
    :include-source:

    class ChangingVelocityModel(pyDeltaRCM.DeltaModel):
        """Model with changing flow velocity.

        Create a model that changes the inlet flow velocity throughout the run.
        In this example, the velocity is changed on each timestep, and the value
        it is set to is interpolated from a **predetermined timeseries** of
        velocities.
        """
        def __init__(self, input_file=None, end_time=86400, **kwargs):

            # inherit from the base model
            super().__init__(input_file, **kwargs)

            # set up the attributes for interpolation
            self._time_array, self._velocity_array = create_velocity_array(
                end_time)  # use default shape parameters for the array

        def hook_solve_water_and_sediment_timestep(self):
            """Change the velocity."""

            # find the new velocity and set it to the model
            self.u0 = np.interp(self._time, self._time_array, self._velocity_array)

            # update other boundary conditions using u0
            self.create_boundary_conditions()

            # log the new value
            _msg = 'Changed velocity value to {u0}'.format(
                u0=self.u0)
            self.log_info(_msg, verbosity=0)


and then simply run with:

.. plot::
    :context:

    # we create the model here, just to be sure it works (for good docs)
    with pyDeltaRCM.shared_tools._docs_temp_directory() as output_dir:
        mdl = ChangingVelocityModel(
            end_time=86400*100,
            out_dir=output_dir)

.. code::

    mdl = ChangingVelocityModel(end_time=end_time)

    while mdl.time < end_time:
        mdl.update()


.. note::

    For information on updating boundary conditions after changing certain model parameters see :doc:`updating_boundary_conditions`.
