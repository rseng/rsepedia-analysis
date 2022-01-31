---
title: 'OMG: A Scalable and Flexible Simulation and Testing Environment Toolbox for Intelligent Microgrid Control'
tags:
  - Python
  - OpenModelica
  - Microgrids
  - Reinforcement Learning
  - Energy Systems
  - Simulation
  - Testing
  - Control
authors:
  - name: Stefan Heid
    affiliation: 1
  - name: Daniel Weber
    affiliation: 2
  - name: Henrik Bode
    affiliation: 2
  - name: Eyke Hüllermeier
    affiliation: 1
  - name: Oliver Wallscheid
    orcid: 0000-0001-9362-8777
    affiliation: 2
affiliations:
 - name: Chair of Intelligent Systems and Machine Learning, Paderborn University, Paderborn, Germany
   index: 1
 - name: Chair of Power Electronics and Electrical Drives, Paderborn University, Paderborn, Germany
   index: 2
date: 25 May 2020
bibliography: paper.bib
---

# Summary

The OpenModelica Microgrid Gym (OMG) toolbox provides a transient simulation framework for local energy grids based on power electronic converters. OpenModelica is used as the backend, allowing users to set up arbitrary electric grid designs via its well-known graphical user interface in a plug-and-play fashion [@Fritzson2018]. Simulations can be configured using a python interface, making it easy to integrate software modules for the realization and testing of closed control loops. In addition, the OpenAI Gym interface is provided to connect data-driven reinforcement learning algorithms for investigating intelligent microgrid control approaches [@OpenAI2016]. 

![\label{fig:omg}](omg.png)
_Fig. 1:  Overview of the interconnections between the different parts of the  OMG  toolbox.  The  OpenModelica and  OpenAIGym logos are the property of their respective owners._


# Background on microgrids and their control

Micro- and smart-grids (MSG) play an important role for integrating renewable energy sources in conventional electricity grids and for providing power supply in remote areas [@Lund2017]. 
Due to their high efficiency and flexibility, power electronic converters are largely used to drive modern MSG.
Power electronics describes the application of solid-state electronics to the control and conversion of electric power, which is largely performed with semiconductor switching devices such as diodes or power transistors.
This includes energy conversion in terms of voltage and current amplitude, frequency and phase angle, as well as the number of phases between two or more electrical energy systems to be connected.


Controlling MSGs is a challenging task due to the high requirements on energy availability, safety, and voltage quality. 
This is particularly demanding due to the wide range of different MSG topologies depending on their field of application like industrial campuses, residential areas or remote off-grid electrification [@Kroposki2008].
This results in high demand for comprehensive testing of new control concepts during their development phase and comparisons with the state of the art to ensure their feasibility.
This applies in particular to data-driven control approaches such as reinforcement learning (RL), the stability and operating behavior of which cannot be evaluated a priori [@Garcia2015].


# State of field 
``OMG`` is a Python-based package for the modeling and simulation of microgrids based on power electronic energy conversion.
The OpenModelica [@Fritzson2018] library enables the user to define their microgrid (i.e. a local electricity grid containing arbitrary sources, storages and loads) in a flexible and scalable way or to use certain predefined example grids. 
Due to the component-oriented modeling framework based on OpenModelica, dynamic processes on small time scales are in focus, which allows for accurate control and test investigations during transients and steady-state.
This is an essential difference to already available open-source solutions for the simulation of electrical energy networks, which, in contrast, generally depict large-scale transmission networks with abstracted models in the (quasi)-stationary state (e.g. PyPSA [@Brown2018] or Pandapower [@Thurner2018]). In addition to the pure modeling and simulation of microgrids, basic building blocks for setting up a hierarchical control framework on the inner and primary level [@Guerrero2013] are provided with ``OMG``. 


# Interfaces for control and reinforcement learning  
The API is designed to provide a user-friendly interface to connect a modeled microgrid (the simulation environment) with a wide range of control methods such as classical linear feedback control or model predictive control techniques (cf. Fig. 1). Moreover, the standardized OpenAI Gym interface [@OpenAI2016] is also available for training data-driven control approaches like RL. 
This enables users who want to integrate contemporary open-source Python-based RL toolboxes such as ``Stable Baselines3`` [@stable-baselines3], ``TF-Agents`` [@TFAgents] or ``keras-rl`` [@plappert2016kerasrl].
Many auxiliary functionalities for the essential operation of microgrids are shipped with OMG such as coordinate transformations for basic controller classes, monitoring wrappers, and phase-locked loops for frequency and phase angle extraction. 
Following this structure, nearly every control approach, including data-driven RL, can be implemented and tested with ``OMG`` in a relatively short amount of time. 
To highlight the challenges of data-driven control approaches in safety-critical environments, application examples using safe Bayesian optimization [@Berkenkamp2020] for automated controller design are provided in the toolbox. 


# Intended use and targeted audience
``OMG`` is designed to be used by students, academics, and industrial researchers in the field of control and energy engineering and data science. The primary objective of the toolbox is to facilitate entry for new users into the modeling, control, and testing of microgrids and to provide a platform on which different control methods (including RL) can be compared under defined conditions (benchmarks).



# Features

The ``OMG`` toolbox provides the following key features:


* A library for the scalable and flexible design of local electricity grids in OpenModelica.
Users can select between a wide range of different grid components and connect them in a plug-and-play approach.

* Dynamic simulation of local electricity grids on component level including single and multi-phase systems as well as AC and DC operation. 

* Easy exchange of models between computing platforms and simulation of the models by using the FMI 2.0 standard [@FMI2020] with C++ code inside and PyFMI [@PyFMI2020] for access in Python. Appropriate numeric solvers for the underlying system of ordinary differential equations can be easily chosen within the usual Python packages (e.g. SciPy) due to the usage of co-simulation. 

* Calculation, evaluation and monitoring of every single time step covering states, action and auxiliary quantities provides an interface for manual or automated inspection. The latter is particularly useful for the automatic training of data-driven control approaches such as reinforcement learning.

* Large variety of predefined and parameterizable controllers (droop, voltage, current in multi- and singlephase) available, easy implementation of user-defined control structures possible.

* Monitoring tools to follow the live performance of the RL agent and to map the overall grid behaviour depending of each selected parameter set

* Interesting use cases applying safe data-driven learning to highlight the requirement of safety in a delicate control environment are available.

# Examples

Detailed examples are shown in the OMG whitepaper (https://arxiv.org/pdf/2005.04869.pdf) including the implementation and evaluation of 
a safe Bayesian controller [@Berkenkamp2020].
The SafeOpt learning algorithm is applied to an automatic controller tuning problem with safety-relevant state constraints in different microgrid topologies (e.g. different number of inverters, load characteristics). 
Furthermore, the provided evaluation tools enable users to compare the performance of different RL algorithms against each other and against manually tuned inverters. 



# Availability and installation

``OMG`` is supported and tested on Linux and Windows. Mac users are asked to run this toolbox on a Linux VM. 
The package should be installed in a conda environment. ``PyFMI`` can be installed via `conda install -c conda-forge pyfmi`, the ``OMG`` package by `pip` Python package manager using 
`pip install openmodelica_microgrid_gym` command. The source code, guide and 
datasets are available on the GitHub repository (https://github.com/upb-lea/openmodelica-microgrid-gym). 

# Individual contributions of the authors

Following are shown the main fields of each individual contributor of OMG: 

* S. Heid: Main software architecture, software module integration, unit tests 

* D. Weber: Application examples, control-related auxiliary features (e.g. basic expert controllers), unit tests

* H. Bode: Design of the specific OpenModelica library created for OMG, grid modelling and simulation, data transfer between OpenModelica and Python, unit tests

* O. Wallscheid: Concept design and idea generation, testing and technical feedback, administrative project management

* E. Hüllermeier: Administrative project management, concept-oriented feedback



# Acknowledgements

The authors kindly acknowledge the valuable contributions and advice regarding grid and controller design by Jarren Lange. The authors would also like to acknowledge the funding and support of this work by the Paderborn University research grant. 

# References

.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/upb-lea/openmodelica_microgrid_gym/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

OpenModelica Microgrid Gym could always use more documentation, whether as part of the
official OpenModelica Microgrid Gym docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/upb-lea/openmodelica_microgrid_gym/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `openmodelica_microgrid_gym` for local development.

1. Fork the `openmodelica_microgrid_gym` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/openmodelica_microgrid_gym.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv openmodelica_microgrid_gym
    $ cd openmodelica_microgrid_gym/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass pytest::

    $ pytest

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.5, 3.6, 3.7 and 3.8, and for PyPy. Check
   https://travis-ci.com/upb-lea/openmodelica_microgrid_gym/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

    $ pytest tests.test_openmodelica_microgrid_gym


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

    $ bump2version patch # possible: major / minor / patch
    $ git push
    $ git push --tags

Travis will then deploy to PyPI if tests pass.
=======
Credits
=======

Development Lead
----------------

* LEA - Uni Paderborn <upblea@mail.upb.de>

Contributors
------------

None yet. Why not be the first?
=======
History
=======

Next
-------

0.4.0 (2021-04-07)
------------------
Changes
^^^^^^^
* ModelicaEnv:
    - Introduced action clipping
    - model_params: None values are not passed to the OpenModelica env to allow initialization
    - model_params: negative time values are introduced for initialization (fix)
    - Introduced abort reward in env if episode is terminated
    - Introduced obs_output to define a subset of history given as observation to the agent

Fix
^^^
* omg.net.MasterInverter:
    - default values used to overwrite passed values

Add
^^^
* Random Process wrapper
* ObsTempl test
* reset test for initialized env




0.3.0 (2020-12-18)
------------------

API
^^^
* ModelicaEnv:
    - Uses Network
    - __init__:
      - removed: timestep, model_output, model_input
      - added: network
    - Delay buffer
* Network and Components:
    - Specify class structure using config file corresponding to fmu (see net-folder)
    - added noise
* SafeoptAgent:
    - __init__: Performance parameters and calculation
* aux_ctl.Contoller:
    - __init__: timestep and undersampling changed
    - added output clipping
* Plotmanager


Examples
^^^^^^^^
* updated to changed API

Experiments
^^^^^^^^^^^
* model validation:
    - experiment files
    - experiment environment managing testbench connection via SSH

Dependencies
^^^^^^^^^^^^
* Decreased Language Level to Python 3.7





0.2.0 (2020-05-27)
------------------


API
^^^
* ModelicaEnv:
   - reward function parameter
   - vis_cols now also supports Plotting templates

* EmptyHistory and descendant: update(), append()
* Agent: added properties
* StaticControlAgent and descendant: small changes in constructor params, specifically obs_template, added properties
* SafeOptAgent: added properties
* Runner: plotting can be disabled

Examples
^^^^^^^^
* added example for plotting

Performance
^^^^^^^^^^^
* 6.6× speedup

Dependencies
^^^^^^^^^^^^
* Increased Language Level to Python 3.8



0.1.3 (2020-05-13)
------------------

* best parameter set output after termination of SafeOpt agent (`#7`_)
* proper action and observation space (`#14`_)
* resolved problem related to environment :code:`model_params` (`#21`_)

|

* documentation improvements (more examples, installation)

.. _`#7`: https://github.com/upb-lea/openmodelica-microgrid-gym/issues/7
.. _`#14`: https://github.com/upb-lea/openmodelica-microgrid-gym/issues/14
.. _`#21`: https://github.com/upb-lea/openmodelica-microgrid-gym/issues/21


0.1.2 (2020-05-04)
------------------

* corrected pip install requirements


0.1.1 (2020-04-22)
------------------

* First release on PyPI.
==========================
OpenModelica Microgrid Gym
==========================

| |build| |cov| |nbsp| |nbsp| |python| |pypi| |download| |nbsp| |nbsp| |license|
| |doc| |whitepaper| |joss|

.. |nbsp|   unicode:: U+00A0 .. NO-BREAK SPACE

.. |build| image:: https://github.com/upb-lea/openmodelica-microgrid-gym/actions/workflows/build_and_test.yml/badge.svg
    :target: https://github.com/upb-lea/openmodelica-microgrid-gym/actions/workflows/build_and_test.yml

.. |cov| image:: https://codecov.io/gh/upb-lea/openmodelica-microgrid-gym/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/upb-lea/openmodelica-microgrid-gym

.. |license| image:: https://img.shields.io/github/license/upb-lea/openmodelica-microgrid-gym
    :target: LICENSE

.. |python| image:: https://img.shields.io/pypi/pyversions/openmodelica-microgrid-gym
    :target: https://pypi.python.org/pypi/openmodelica_microgrid_gym

.. |pypi| image:: https://img.shields.io/pypi/v/openmodelica_microgrid_gym
    :target: https://pypi.python.org/pypi/openmodelica_microgrid_gym

.. |download| image:: https://img.shields.io/pypi/dw/openmodelica-microgrid-gym
    :target: https://pypistats.org/packages/openmodelica-microgrid-gym

.. |doc| image:: https://img.shields.io/badge/doc-success-success
    :target: https://upb-lea.github.io/openmodelica-microgrid-gym

.. |whitepaper| image:: https://img.shields.io/badge/arXiv-whitepaper-informational
    :target: https://arxiv.org/pdf/2005.04869.pdf
    
.. |joss| image:: https://joss.theoj.org/papers/10.21105/joss.02435/status.svg
   :target: https://doi.org/10.21105/joss.02435



.. figure:: https://github.com/upb-lea/openmodelica-microgrid-gym/raw/develop/docs/pictures/omg_flow.png

**The OpenModelica Microgrid Gym (OMG) package is a software toolbox for the
simulation and control optimization of microgrids based on energy conversion by power electronic converters.**

The main characteristics of the toolbox are the plug-and-play grid design and simulation in OpenModelica as well as
the ready-to-go approach of intuitive reinfrocement learning (RL) approaches through a Python interface.

The OMG toolbox is built upon the `OpenAI Gym`_ environment definition framework.
Therefore, the toolbox is specifically designed for running reinforcement
learning algorithms to train agents controlling power electronic converters in microgrids. Nevertheless, also arbritary classical control approaches can be combined and tested using the OMG interface.

.. _OpenAI Gym: https://gym.openai.com/

* Free software: GNU General Public License v3
* Documentation: https://upb-lea.github.io/openmodelica-microgrid-gym


Video Tutorial
--------------

Following is a short YouTube video introduction, to get a fist impression how to use OMG.



- https://www.youtube.com/watch?v=rwBNFvCi_dY

Installation
------------


Install Python Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the short installation guide for Windows and Linux. OpenModelica_ is hardly supported for Mac, they suggest to install in a Linux VM. For this reason, running OMG in a Linux VM is strongly recommended for Mac users!

Since it is not possible to install PyFMI_, a package which is necessary for the communication between the python interface and the environment, via pip, we recommend to install this package in advance in a conda environment.
As of now, only Windows and Linux are supported officially.

- If conda is NOT installed on your PC, install miniconda_ for python 3.8
- Create a new conda environment (e.g. in PyCharm)
- Install PyFMI from the conda-forge channel in the terminal::

    $ conda install -c conda-forge pyfmi


- Install OpenModelica MicrogridGym from PyPI (recommended)::

    $ pip install openmodelica_microgrid_gym

.. _OpenModelica: https://openmodelica.org/download/download-mac
.. _miniconda: https://conda.io/en/latest/miniconda.html
.. _PyFMI: https://github.com/modelon-community/PyFMI

Installation of OpenModelica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OMG was create by using OMEdit_ v1.16

In case of installation issues you can resort to their pre-built `virtual machine`_.

.. _OMEdit: https://openmodelica.org/download/download-windows
.. _virtual machine: https://openmodelica.org/download/virtual-machine

Getting started
---------------

The environment is initialized and run like any other OpenAI Gym

.. code-block:: python

    import gym

    if __name__ == '__main__':
        env = gym.make('openmodelica_microgrid_gym:ModelicaEnv-v1',
                   max_episode_steps=None,
                   net='../net/net.yaml',
                   model_path='../omg_grid/grid.network.fmu')

        env.reset()
        for _ in range(1000):
            env.render()
            env.step(env.action_space.sample())  # take a random action
        env.close()




OMG uses the `FMI standard`_ for the exchange of the model between OpenModelica and Python.

.. _FMI standard: https://fmi-standard.org/

An example network consisting out of two inverters, three filters and an inductive load.

.. figure:: https://github.com/upb-lea/openmodelica-microgrid-gym/raw/master/docs/pictures/omedit.jpg

You can either use one of the provided FMUs (Windows and Linux, 64-bit, both included in the grid.network.fmu) or create your own by running::

    openmodelica_microgrid_gym\fmu> omc create_fmu.mos

Windows users might need to open the terminal out of OpenModelica by clicking 'tools' => 'OpenModelica Command Prompt' to make sure that the command 'omc' gets recognized.

Running the ``staticctrl.py`` starts a simulation with a manually tuned cascaded PIPI controller

.. figure:: https://github.com/upb-lea/openmodelica-microgrid-gym/raw/master/docs/pictures/control.jpg
    :scale: 70%
    :align: center

A save Bayesian approach of a reinforcement learning agent is provided under examples/berkamkamp.py.

.. figure:: https://github.com/upb-lea/openmodelica-microgrid-gym/raw/master/docs/pictures/kp_kp_J.png
    :figwidth: 60%
    :align: center

Using pytest
^^^^^^^^^^^^

OMG provides a big range of tests to ensure correct working toolbox after changes are done.
On some windows machines, the tests can only be started from the terminal via 'pytest'.

The standard test OS for the development is Linux. In some cases, we have noticed that the test_modelica.py on windows PCs might throw an error.
Since on Linux everything works fine, it seems to be a numerical issue connected with the FMUs.


Citation & white paper
----------------------

Please find a white paper on the OMG toolbox including an exemplary usage scenario here:

- https://arxiv.org/abs/2005.04869

Please use the following BibTeX entry for citing us::

    @article{OMG-code2020,
        title = {OMG: A Scalable and Flexible Simulation and Testing Environment Toolbox for Intelligent Microgrid Control},
        author = {Stefan Heid and Daniel Weber and Henrik Bode and Eyke Hüllermeier and Oliver Wallscheid},
        year = {2020},
        doi = {10.21105/joss.02435},
        url = {https://doi.org/10.21105/joss.02435},
        publisher = {The Open Journal},
        volume = {5},
        number = {54},
        pages = {2435},
        journal = {Journal of Open Source Software}
    }

    @article{OMG-whitepaper2020,
        title={Towards a Scalable and Flexible Simulation and
               Testing Environment Toolbox for Intelligent Microgrid Control},
        author={Henrik Bode and Stefan Heid and Daniel Weber and Eyke Hüllermeier and Oliver Wallscheid},
        year={2020},
        eprint={http://arxiv.org/abs/2005.04869},
        archivePrefix={arXiv},
        primaryClass={eess.SY}
    }


Contributing
------------

Please refer to the `contribution guide`_.

.. _`contribution guide`: https://github.com/upb-lea/openmodelica-microgrid-gym/blob/master/CONTRIBUTING.rst


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
Welcome to OpenModelica Microgrid Gym Toolbox Documentation!
=====================================================================

The OpenModelica Microgrid Gym (OMG) package is a software toolbox for the simulation of power electronics-driven microgrids and to
train and test reinforcement learning agents and to compare them with classical control approaches.

Content
*******

In the examples section all available use cases are presented with their default configuration.
For quick start, one of these can be selected and used out of the box.


The documentation of the base classes is important for the development of own modules like further reward functions or
reference generators. In this part, the basic interfaces of each module are specified.
For the creation of additional grid constellations, Openmodelica (nightly build recommended) can be used.


.. toctree::
   :maxdepth: 4
   :titlesonly:
   :caption: User Guide:

   parts/user_guide/getting_started
   parts/user_guide/OpenModelica
   parts/user_guide/fmu
   parts/user_guide/Pythoncode
   parts/user_guide/examples
   parts/user_guide/controller_tuning


.. toctree::
   :maxdepth: 4
   :titlesonly:
   :caption: API:

   api/omg.agents
   api/omg.aux_ctl
   api/omg.env
   api/omg.execution
   api/omg.net
   api/omg.util


.. GENERATE APIDOC
.. - sphinx-apidoc -o docs/api openmodelica_microgrid_gym/ -e
.. - delete module.rst
.. - remove package and module names:
..   - execute regex '.* package$' ''
..   - execute regex '.* module$' ''


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Examples
========

.. toctree::
   :maxdepth: 4
   :titlesonly:
   :hidden:

   examples/creating_env
   examples/basic_agent
   examples/plotting

   examples/single_inverter_current_control_safe_opt
   examples/two_inverter_static_droop_control




Basic Examples
--------------

To get an idea how the toolbox works and how to set-up an environment, an agent or use plotting, three very basic examples are shown below. They use the network presented above, but only the first inverter will be used.

- `Creating the Environment`_
- `Writing an simple Agent`_
- `Modify Plotting`_

.. _`Creating the Environment`: examples/creating_env
.. _`Writing an simple Agent`: examples/basic_agent
.. _`Modify Plotting`: examples/plotting



Advanced Examples
-----------------

- `Single Inverter Current Control`_
- `Two Inverter Static Droop Control`_

.. _`Single Inverter Current Control`: examples/single_inverter_current_control_safe_opt
.. _`Two Inverter Static Droop Control`: examples/two_inverter_static_droop_control
Toolbox Installation and General Remarks
================================


In the following, some installation hints and an introduction to the Python code written for the OMG
toolbox is presented.


Installation and Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the short installation guide for Windows and Linux. `OpenModelica <https://openmodelica.org/download/download-mac>`__ is hardly supported for Mac, they suggest to install in a Linux VM. For this reason, running OMG in a Linux `VM <https://openmodelica.org/download/virtual-machine>`__ is strongly recommended for Mac users!



It is recommended to install OMG via pip:

::

    pip install openmodelica_microgrid_gym

Alternatively, you can clone the GitHub repository. A list of
requirements is provided in the
home-directory.

.. literalinclude:: ../../../requirements.txt


**Hint:** If you are running a windows, PyFMI might throw some errors
while installing via pip. It can be installed via *conda* by running:

::

    conda install -c conda-forge pyfmi 

Simulation Settings
~~~~~~~~~~~~~~~~~~~

Heart of the program structure is the creation of the environment via
**gym.make()** in the main programm (in the folder example). Nearly
every simulation setting can be done directly in here. Some of the most
important ones are described in the following. For further information, see the `API-documentation <../../api/omg.env.modelica.html>`__.

-  **time\_step:** step size of the simulation in seconds. Too large
   timesteps may result in numerical issues, small timesteps result in a
   high simulation time. 1e-4 seems to be a good compromise as many real
   controllers operate in timesteps like this.

-  **reward\_fun:** Callable - Reward function for the RL-algorithm. Minimal value
   of rewards is for example used as lower bound for the safe Bayseian
   algorithm (see single_inverter_current_control_safe_opt.py). Has to be adjusted problem-specific.

-  **solver\_method:** Solver used for the ODE system. Every solver from
   `scipy.integrate.solve\_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`__
   can be selected. Though it is highly recommended to use the implicit
   ones. Default solver is the "LSODA" with its integrated stiffness
   detection. Non-stiff systems become solved with the faster "Adams"
   method, stiff systems with "BDF".


-  **model\_params:** Parameters for the simulation, which should be
   changed compared to the default values from the OpenModelica model.
   Also usable for loadsteps as replacement for
   `switches <OpenModelica.md>`__.

Example which increases the resistors in the load after 0.2 seconds from
20 Ohm to 40 Ohm:

::

    def f(t):
        return 20 if t < .2 else 40

    model_params={'rl.resistor1.R': f, 'rl.resistor.R': f, 'rl.resistor.R': f},

The function :code:`f` is passed as a callable (function reference).
The environment will evaluate this function in every timestep passing this timestep as parameter to the function
automatically.


Setting of v\_DC
~~~~~~~~~~~~~~~~

The DC supply voltage v\_DC can be set either directly in the
`OpenModelica model <OpenModelica.html#setting-of-v-dc>`__ or via
Python. The default value is 1000 V. It can be changed in the
environment creation with the line:

::

    model_params={'inverter1.v_DC': 700, 'inverter2.v_DC': 500}, 

It will be set for every of the three phases of the inverter. Take care
to set the param for every inverter which should no have the default
supply voltage of 1000 V.

Data Logging
~~~~~~~~~~~~

To enable logging, the root logger needs to be initialized in the
main function. To do so, call:

::

    import numpy as np

    logging.basicConfig()

    if __name__ == '__main__':
        ctrl = dict()

For further information about logging and the level see the `logging
standard library <https://docs.python.org/3/library/logging.html>`__.
Controller Tuning Hints
=======================

1. Current controller of primary inverter

-  With no droop, in other words constant mains frequency, apply a short
   circuit to inverter and tune Kp, Ki of the current controller. Can
   tune Kp then Ki, finally Kp again if both can’t be tuned at the same
   time.

   -  Try tune for 90-95% of max current. Not peak current limit, but
      maximum allowed during nominal operation.
   -  Ensure that the tuning does not allow the current to exceed the
      peak current limit. In real world scenarios this is when an
      inverter will shutoff or explode.

2. Voltage controller of primary inverter

-  With no droop and an open circuit load of just the inverter tune Kp,
   Ki of the voltage controller. Can tune Kp then Ki, finally Kp again
   if both can’t be tuned at the same time.

3. PLL of secondary inverter

-   With primary inverter providing a constant frequency start tuning
    Kp, Ki of the PLL. Noise on the frequency output of the PLL is
    acceptable, as long as the output phasors of the PLL accurately
    match the incoming voltage signal.
-   The secondary inverter power electronics should be
    disconnected/open-circuit for this.
-   If possible inject step changes to the frequency setpoint of the
    primary inverter, watching how accurately the PLL tracks the
    external voltage reference. Continue tuning if necessary.

4. Current controller of secondary inverter

-   With the PLL of the secondary inverter locked to the primary inverter
    tune the Kp, Ki of the current controller. No droops at this stage.

   -  Might need to create a step change for the setpoint for this to
      test accurate tracking.
   -  Try tune for 90-95% of max current. Not peak current limit, but
      maximum allowed during nominal operation.

5. Droop of the primary inverter

-  This isn’t really a tuning step as the parameters won’t affect any
   other system.

6. Droops of the secondary inverter

-  Firstly, the filter for the frequency and the voltage feedback to the
   droop controllers should be set before tuning the droop.
-  Tune the droop controllers of the secondary inverter.

   -  The frequency of the filter (for the frequency feedback from the
      PLL) affects the droop parameters and should thus be considered in
      the tuning.


Getting Started
===============

.. figure:: ../../pictures/omg_flow.png
   :alt: 

This user guide covers all of OpenModelica Microgrids Gym (OMG) toolbox features by
topic area. Each of the following steps is introduced in a more detailed
way in the different chapters of this users guide.

First step is to `create the microgrid <OpenModelica.html>`__, which is
the environment for training reinforcement learning agents in power electronic-driven microgrids.
In addition, the OMG toolbox can be used for pure simulation and classical control purpose using OpenModelica models with a Python interface.

Each microgrid model is built in the open source software
`OpenModelica <https://www.openmodelica.org/>`__ and can be easily adapted.

.. figure:: ../../pictures/network1.png
   :alt: 

For the transfer to Python, it needs to get exported as a `Functional
Mock-up Unit (FMU) <https://fmi-standard.org/>`__.

The creation process of the FMU is shown `here <fmu.html>`__. It is used to
build a gym environment like in the examples from `OpenAI
Gym <https://gym.openai.com/>`__. In OMG, the gym environment is defined
for example in (examples/two_inverter_static_droop_control.py).

After creating the environment, the network can be simulated in Python.
On the one hand, it is possible to test predefined, static controller designs
like described `here <examples.html#two-inverter-static-droop-control-py>`__.

.. figure::  ../../pictures/abc.png
   :alt: 

However, the main function of this toolbox is to apply reinforcement learning
approaches by utilizing the OMG interface for optimal microgrid control as shown in this
`example <examples.html#single-inverter-current-control-safe-opt-py>`__.

.. figure::  ../../pictures/kp_kp_J.png
   :alt: 

Basic Examples
~~~~~~~~~~~~~~

To get an idea how the toolbox works and how to set-up an environment or an agent, two very basic examples are shown below. Both use the network presented above, but only the first inverter will be used.


Creating an environment
-----------------------


Following is a minimal example how to set-up an run an environment.
Necessary is the definition of model inputs, in this case the three phases of the first inverter.
The model outputs will be shown as simulation results, and the model path is the relative location of the FMU file, which contains the network.
For any other simulation parameters, for example the step-size, default values will be used.

For the initialisation, the environment needs to be reseted, and env.render will define the output plots.
The simulation will perform 1000 steps. A different random number will be provided to every of the three previously defined model_inputs.
Afterwards, the inductor currents of the LC-filter "lc1"shown in the figure above will be plotted, which should result in three increasing and due to the random function noisy  lines.

.. literalinclude:: ../../../examples/basic_env.py
   :linenos:



Creating an agent and a runner
------------------------------

Additionally to the environment, an an agent will be created and a runner will be used. The runner class will take care of initializing and termination
of agents and environments, as well as the execution of multiple episodes. The class will handle all information
exchange between agent and environment like presented in the high level code architecture shown below:

.. figure:: ../../pictures/highlevel.png
   :width: 400
   :alt:

Since the inputs are used for both the agent and the environment, they are defined in advance. Although the Agent gets information of the environment, in this small example, its action is still a random number.

The environment is the same as above. Afterwards, the agent and the runner get defined, and the runner runs for one episode.

.. literalinclude:: ../../../examples/simple_agent.py
   :linenos:


Functional Mock-up Unit (FMU)
=============================

What is FMI and FMU?
^^^^^^^^^^^^^^^^^^^^

The Functional Mock-up Interface (`FMI <https://fmi-standard.org/>`__)
is a free standard that defines a container and an interface to exchange
dynamic models using a combination of XML files, binaries and C code
zipped into a single file. The generated files for the exchange are
called Functional Mock-up Units (FMU). The binaries are
platform-specific, so a linux user can´t run a FMU which is created on a
windows machine.

There are two versions of the FMI standard. The OMG toolbox uses **FMI
2.0**. Furthermore, there are two simulation types, co-simulation (CS) and
model exchange (ME). Co-simulation has an included solver in the FMU,
model exchange provides the possibility to use solvers in Python. Due to
a lack of implicit solvers in CS, the OMG toolbox uses **ME**.

Create FMU´s
^^^^^^^^^^^^

The best way to create a FMU for the OMG toolbox is to run the
create\_fmu.mos file in the folder "omg_grid" with the terminal. Make sure that
`OpenModelica <https://openmodelica.org/download/download-windows>`__ is
installed and your latest version of your OpenModelica package is in
the same folder as this script.

By running the create\_fmu.mos with the command *omc*, you can create a
FMU of the model *network* in the *Grids* folder.

::

    path\omg_grid> omc create_fmu.mos

Windows users might need to open the terminal out of OpenModelica by clicking 'tools' => 'OpenModelica Command Prompt' to make sure that the command *omc* gets recognized.

Creating FMUs of other models can be generated by changing
*omg_grid.Grids.Network* in the last line of the *create\_fmu.mos* to
*omg_grid.Grids.YourNetworksName*.

::

    OpenModelica.Scripting.loadFile("package.mo"); getErrorString();
    setCommandLineOptions("-d=newInst"); getErrorString();
    setCommandLineOptions("-d=initialization"); getErrorString();
    setCommandLineOptions("--simCodeTarget=Cpp"); getErrorString();
    setCommandLineOptions("-d=-disableDirectionalDerivatives"); getErrorString();
    OpenModelica.Scripting.translateModelFMU(omg_grid.Grids.Network, version="2.0", fmuType = "me"); getErrorString();

The lines in between are setting flags for the FMU-creation. The
target-language is C++ instead of the default C because of a missing
implementation of directional derivatives in C, which are necessary for the solvers in python.

It is possible to create FMUs directly in OpenModelica (File -> Export ->
FMU). Settings for this can be done in Tools -> Options -> Simulation
(Target language and Translation Flags) and Tools -> Options -> FMU for
the Version, type, name and path. This way is not recommended because of
the possibility to miss flags like the initialization. Furthermore,
problems with the providing of directional derivatives occurred during testing.

Merge FMUs
^^^^^^^^^^^

As mentioned above, FMUs only work on the operating system they are
created on. Though, there is a possibility to combine previously
generated FMUs to allow a usage on different platforms.

First, generate FMUs on different platforms from the **same** model.
Then run the Python file *merge*\_fmus.py\. Select the FMUs which
should be combined and choose a target file (does not need to exist
before). The script checks at several points if the FMUs were generated
from the same version of the library (careful: small changes like parameter variation
might not be detected) and if all FMUs contain binaries for different
platforms.

Annotation for Windows: Sometimes, the script throws an error because in
the temp-folder. The script still works and merges
the fmu. The resulting file just needs to be renamed with an ".fmu" at the end and can
directly be used.
OpenModelica
============

`OpenModelica <https://openmodelica.org/>`__ is an open-source
Modelica-based modeling and simulation environment intended for
industrial and academic usage.

Installation of OpenModelica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The models shown below were created by using
`OMEdit <https://openmodelica.org/download/download-windows>`__ v1.16.

Using a Linux OS, sometimes may lead to problems while trying to install
OpenModelica. In this case, try to download the pre-built `virtual
machine <https://openmodelica.org/download/virtual-machine>`__.

Mac users are strongly encouraged to run OpenModelica as well as the OMG toolbox in a Linux VM.

Creating Microgrids with OpenModelica
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The microgrids are created with a  user defined library. To start it, run the file *grid.mo* directly in the folder *omg_grid*.


This package contains all components and required for
creating the power electronics driven microgrids, as well as some example networks.

.. figure:: ../../pictures/library1.jpg
   :alt: 

It contains several folders with predefined components,
filters, loads etc. as well as microgrids for the FMU export to Python and some stand-alone examples which can be run directly in OpenModelica.



.. figure:: ../../pictures/omedit.jpg
   :alt: 

Main components of any microgrid are the three-phase inverters. They consist
of three input voltages controlled by a cascaded PI-PI controller in
Python. Default nomenclature for those inputs is i1p1 for
"**i**\ nverter 1 **p**\ hase 1" etc, but it can be changed in the
python code (model\_input=['one', 'two',...] in the env=gym.make()
call).

**Important**: By using filters/loads with inductors or capacitors,
leave the checkbox for the initialization in the provided settings. Changes are likely to
result in errors while creating the FMU or running the Python files.

The provided examples are designed for up to two inverters, but the underlying models can be
easily extended in order to investigate on more complex microgrid topologies.
Extended example showcases are also planed for future releases.

Power Losses
^^^^^^^^^^^^

In the default example "network", no power losses in the inverters or filters are included.
For the latter they can be added by using parts out of the "Filter.LossesFilter" package instead of
the "Filter.IdealFilter" package. Due to a big increase of components and
equations in the ODE-system, the simulation time will increase.
For modeling losses inside the power electronic converters, adding a model in the Python interface
scripts is recommending. Integrating, e.g switching losses, directly in the OpenModelia model will
require to reduce to simulation step size significantly.

For larger simulations with a demand of power loss modeling, it is recommended to
create user defined filters with only the resistors which are needed for
the calculation. To modify them, create a new package (ctrl+N,
specialization: package), duplicate the part which is to modify in the
new package (right-click on it, duplicate, select the previously created
package in path) and modify it there.

Switches / Transistors
^^^^^^^^^^^^^^^^^^^^^^

Modeling switching-like events inside the model,
e.g. triggering loadsteps by adding or removing loads, is
desirable, but difficult to implement with the possibilities of
OpenModelica. Switches in OpenModelica - like in many other free
modelling languages - are designed as resistors. A closed switch has a
low resistance, an open switch a high one. "Removed" loads are still
connected to the grid. Connections with resistors in such dimension
cause numerical issues while simulating as the ODE system becomes stiff.
There are solvers available for stiff equation systems like BDF and
Radau or ones with automatic stiffness detection, but using the switches
often runs into non-converging systems and execution errors.

The working alternative is a parameter-variation of the load. It is
possible to change the parameters of any load during a simulation and
apply loadsteps in this way (see the topic
`Python code <Pythoncode.html>`__).

Setting of v\_DC
^^^^^^^^^^^^^^^^

The DC supply voltage *v*\_DC can be set either directly in the
OpenModelica model or via `Python <Pythoncode.html#setting-of-v-dc>`__.
In the OM model, doubleclick in your network on the inverter, and change
the parameter *v*\_ DC to the demanded value. All three phases of the
inverter will be supplied with the same DC voltage. The default value is
1000 V. The default value can be changed with a right-click on an
inverter, *open class*, select *text view* on the top left corner of the
model canvas, and change the number in the following code line to
the demanded default value:

::

      parameter Real v_DC = 1000;
      

Phase-Locked Loop (PLL)
^^^^^^^^^^^^^^^^^^^^^^^

The PLL blocks are working for simulations in OpenModelica, but out of
structural reasons, the PLL is calculated in Python.
Creating an environment
=======================


Following is a minimal example how to set-up and run an environment.
Necessary is the definition of model inputs, in this case the three phases of the first inverter.
The model outputs will be shown as simulation results, and the model path is the relative location of the FMU file, which contains the network.
For any other simulation parameters, for example the step-size, default values will be used.

For the initialisation, the environment needs to be reseted, and env.render will define the output plots.
The simulation will perform 1000 steps. A different random number will be provided to every of the three previously defined model_inputs.
Afterwards, the inductor currents of the LC-filter "lc1"shown in the figure above will be plotted, which should result in three increasing and due to the random function noisy  lines.

.. literalinclude:: ../../../../examples/basic_env.py
   :linenos:Single Inverter Current Control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example a three phase inverter is supplying a load (rl1) via a filter (lc1)
like shown in the figure below. From that model a FMU is
built to create the environment.

.. figure:: ../../../pictures/Model.png

An optimization method developed by `Berkenkamp et al.`_ called Safe Controller Optimization (safeopt) is used which takes a Gaussian process and Bayesian
optimization to safely determine "optimal" controller parameters. The
goal of the standard PI current controller is to supply an exemplary 15 A d-current
to the load.

.. _`Berkenkamp et al.`: https://arxiv.org/abs/1509.01066

The `generated FMU <fmu.html>`__ is used in the environment to build up
a gym env like the examples from OpenAI Gym (https://gym.openai.com/).
The gym enviroment is defined in (examples/single\_inverter\_current\_control\_safe\_opt.py).
It generates a gym environment using

- a reward function,
- plotting the inductor values (current) from the lc1-filter (which should be controlled) like shown in the figure below,
- simulating 300 timesteps of delta\_t of the FMU grid.network\_singleInverter.fmu (generated from the model in the plot above),
- using the setpoints for the inverters (modulation indices) i1p{1,2,3} as inputs,
- and the inductor currents and capacitor voltages of lc1-filter as outputs.

.. figure:: ../../../pictures/i_abc_bk_kp15_Ki121.png

The agent used in this simple RL-example is taken from the class
:code:`SafeOptAgent`. It contains the controller a
:code:`MultiPhaseDQCurrentSourcingController`, which consists of multiphase
(3) PI controllers to control the current across the inductor of the
lc1-filter. There are also droop controllers implemented to calculate
e.g. the frequency drop due to load changes. The agent's task is to find better
parameters for the current controllers (Kp & Ki). Therefore, they are
defined as mutable\_params (e.g.
examples/single\_inverter\_current\_control\_safe\_opt.py) to
adopt them between the episodes. The SafeOpt algorithm uses a Gaussian
process to estimate the performance of the controller. Thus, the
bounds and the lengthscale (c.f. examples/single\_inverter\_current\_control\_safe\_opt.py) for
the gain parameters (Kp and Ki) have to be defined.

One can adjust one of the parameters (Kp or Ki) (1D case) or both of them
(2D case) using the algorithm. Therefore, the following flag parameters have to
be adjusted accoridngly:

- To adjust only Kp set :code:`adjust = 'Kp'`
- To adjust only Ki set :code:`adjust = 'Kp'`
- To adjust only Kp and Ki set :code:`adjust = 'Kpi'`

Due to SafeOpt the agent need a safe starting point (Kp and Ki). Then it
tries to calculate safely parameters with better performance. The
performance is calculated using the reward function from the environment.
There the mean-root-error (RME) from the measured currents and the setpoints are
calculated. Additionally a barrier function is used to penalize
over-currents. The barrier function can be adjusted using the parameter mu.

The safe threshold for the agent is set as safe\_threshold-times of
the initial performance (c.f. agents/safeopt.py). For example,
safe\_threshold = 1.2 and the initial reward is -10 the safe threshold
would be -12.

In the end of the script a :code:`Runner` is used to execute 10 episodes
using the agent to control the environment. For every episode the
controlled currents and the performance function as a function of Kp
and/or Ki are plotted.

Some exemplary results are shown below:

-  If :code:`adjust == 'Kp'`, the agent tries to
   find an optimal value for the proportional gain (Kp) of the
   controller in the range of [0, 0.03] with a
   lengthscale of 0.01. In the figure below on the x-axis is
   the value for Kp and on the y-axis the performance value calculated
   using the reward function mentioned above.

.. figure:: ../../../pictures/kp_J.png

-  If :code:`adjust == 'Ki'`, the agent tries to
   find an optimal value for the integral gain (Ki) of the controller in
   the range of [0, 300]  with a lengthscale of 50. In the figure below on the x-axis is the value for Ki and
   on the y-axis the performance value calculated using the reward
   function mentioned above.

.. figure:: ../../../pictures/ki_J.png

The - due to the algorithm - "unsafe" point on the right (for Kp as well
as for Ki) is not due to overcurrent but due to bad performance due to
permanent control error. The resulting currents for Kp = 0.01 and Ki = 0 ("unsafe" point on the right in the figure above)
is shown in the picture below. Due to the high error compared to the
reference value (15 A d-current), the performance is as bad as the
algorithm defines it as unsafe - in comparison to the performance
reached using the initial controller parameters.

.. figure:: ../../../pictures/i_abc_ki_J_bad.png

-  If :code:`adjust == 'Kpi'`, the agent tries to
   find an optimal value for the proportional gain (Kp) as well as for
   the integral gain (Ki) of the controller in the ranges of [0, 0.03]
   and a lengthscale of 0.01 for Kp and a range of [0, 300] and a
   lengthscale of 50 for Ki. In the figure below on the x-axis is the
   value for Kp, the y-axis the value for Ki and the z-axis the
   performance value calculated using the reward function.

.. figure:: ../../../pictures/kp_ki_J.png

The results of the algorithm are printed into the console in the form
like below:

Iteration, performance J, Params [Kp, Ki]

::

           J                                      Params
    0  -0.527522                                [0.01, 10.0]
    1  -0.442648    [0.01517286546392185, 14.85163114970222]
    2  -0.318154    [0.01426989823624961, 44.96747682456248]
    3  -0.296940   [0.007935547159879385, 63.12800825929393]
    4  -0.286636    [0.01482713453607815, 88.70170996759624]
    5  -0.286815  [0.006770598304777539, 108.12303673537075]
    6  -0.280167  [0.013261084415467694, 135.24448051372738]
    7  -0.313204   [0.02201710533671064, 56.446583269542394]
    8  -1.387003  [0.022868977920736434, 108.40140778199653]
    9  -0.304403   [0.002145673177669012, 55.14569829606201]
    10 -0.480421   [0.026197353734745858, 22.29566509028389]
    11 -1.097157  [0.0055262530542335535, 157.4879776902759]
    12 -0.391706                    [0.0, 17.86728037560901]
    13 -1.307038                    [0.0, 106.0724160092763]
    14 -1.561142                    [0.03, 42.1020413015999]

The best performance in this short example of -0.280167 produces the
parameter set of Kp = 0.0132... and Ki = 135.244...
Creating an agent and a runner
==============================

Additionally to the environment, an an agent will be created and a runner will be used. The runner class will take care of initializing and termination
of agents and environments, as well as the execution of multiple episodes. The class will handle all information
exchange between agent and environment like presented in the high level code architecture shown below:

.. figure:: ../../../pictures/highlevel.png
   :width: 400
   :alt:

Since the inputs are used for both the agent and the environment, they are defined in advance. Although the Agent gets information of the environment, in this small example, its action is still a random number.

The environment is the same as above. Afterwards, the agent and the runner get defined, and the runner runs for one episode.

.. literalinclude:: ../../../../examples/simple_agent.py
   :linenos:Two Inverter Static Droop Control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, a FMU generated by OpenModelica as gym environment containing two inverters, each connected via a
filter to supply in parallel a RC load is used which is shown in the figure below.
This example uses the controllers as defined in the auxiliaries. One inverter is set up as voltage forming inverter with a
direct droop controller which e.g. frequency drops due to the applied power. The other controller is used as current
sourcing inverter with an inverse droop controller which reacts on the frequency and voltage change due to its droop
control parameters by a power/reactive power change.
In the default settings, plots of the abc signal as well as the dq0 signals of
the master and slave are provided.

By default, the following small network will be simulated:

.. figure:: ../../../pictures/network.png

A short introduction to experimental controller tuning with some hints
can be found `here <controller_tuning.html>`__.

If the controller works fine, a three phase voltage similar to the
following one should be one of the plots.

.. figure:: ../../../pictures/abc.png

Any other demanded signal which is provided by the FMU or saved during
the simulating can be plotted by adding it to

::

    viz_cols=['*.m[dq0]', 'slave.freq', 'lcl1.*'],

in the gym.make() command. Make sure that demanded signal from the fmu
are listed as a model\_output:

::

    model_output={
                       'lc1': [
                           ['inductor1.i', 'inductor2.i', 'inductor3.i'],
                           ['capacitor1.v', 'capacitor2.v', 'capacitor3.v']],
                       'rl1': [f'inductor{i}.i' for i in range(1, 4)],
                       'lcl1':
                           [['inductor1.i', 'inductor2.i', 'inductor3.i'],
                            ['capacitor1.v', 'capacitor2.v', 'capacitor3.v']]},
                       )

Hint: Every possible variable which is provided by the FMU can be seen
the easiest in OpenModelica. Run the simulation without input signals,
so every result for voltages and currents should be 0. On the bottom right side, you can select
each component of the model in the tree structure. Clicking through the
components until reaching the variable will show the whole variable name
(for example lcl2.inductor2.i) on top of the plotting window.

The parameters of the controller like the control frequency delta\_t,
the voltage, frequency or droop characteristics can be set directly in
the main function.Creating plots of environment variables
=======================================

In this example is shown how the environment variables can be plotted and how the plots can be adjusted.
The environment gets an object of the class PlotTmpl, where can be chosen which environment variables are visualized.
For style, linestyle and color the grouping is detected in comparison to the environment variables to be plotted.
Ajdusting labels and saving the figure can be done via callable, like shown in the example below.

.. literalinclude:: ../../../../examples/plotting.py
   :linenos:omg.util.fastqueue
==================================================

.. automodule:: openmodelica_microgrid_gym.util.fastqueue
   :members:
   :undoc-members:
   :show-inheritance:
omg.env.pyfmi
=============================================

.. automodule:: openmodelica_microgrid_gym.env.pyfmi
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents.staticctrl
=======================================

.. automodule:: openmodelica_microgrid_gym.agents.staticctrl
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl.inverter\_controllers
================================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.inverter_controllers
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl.base
======================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.base
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl.filter
========================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.filter
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents
=============================

Submodules
----------

.. toctree::

   omg.agents.agent
   omg.agents.episodic
   omg.agents.safeopt
   omg.agents.staticctrl
   omg.agents.util

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.agents
   :members:
   :undoc-members:
   :show-inheritance:
omg.net
========================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   omg.net.base
   omg.net.components

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.net
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux\_ctl.droop\_controllers
===============================================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.droop_controllers
   :members:
   :undoc-members:
   :show-inheritance:
omg.env.modelica
==================================

.. automodule:: openmodelica_microgrid_gym.env.modelica
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl.pi_controllers
====================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.pi_controllers
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl
==================================

Submodules
----------

.. toctree::

   omg.aux_ctl.base
   omg.aux_ctl.filter
   omg.aux_ctl.params
   omg.aux_ctl.observers
   omg.aux_ctl.droop_controllers
   omg.aux_ctl.pi_controllers
   omg.aux_ctl.inverter_controllers

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.aux_ctl
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux_ctl.params
========================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.params
   :members:
   :undoc-members:
   :show-inheritance:
omg.execution
================================

Submodules
----------

.. toctree::

   omg.execution.runner
   omg.execution.callbacks

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.execution
   :members:
   :undoc-members:
   :show-inheritance:
omg.env
==========================

Submodules
----------

.. toctree::

   omg.env.modelica
   omg.env.pyfmi
   omg.env.plot
   omg.env.plotmanager

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.env
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents.util
=================================

.. automodule:: openmodelica_microgrid_gym.agents.util
   :members:
   :special-members:
   :exclude-members: __dict__,__module__,__repr__,__weakref__
   :undoc-members:
   :show-inheritance:
omg.util.randproc
==================================

.. automodule:: openmodelica_microgrid_gym.util.randproc
   :members:
   :undoc-members:
   :show-inheritance:
omg.execution.callbacks
=======================================================

.. automodule:: openmodelica_microgrid_gym.execution.callbacks
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents.agent
==================================

.. automodule:: openmodelica_microgrid_gym.agents.agent
   :members:
   :undoc-members:
   :show-inheritance:
omg.util.recorder
==================================

.. automodule:: openmodelica_microgrid_gym.util.recorder
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents.episodic
===================================================

.. automodule:: openmodelica_microgrid_gym.agents.episodic
   :members:
   :undoc-members:
   :show-inheritance:
omg.util.transforms
=======================================

.. automodule:: openmodelica_microgrid_gym.util.transforms
   :members:
   :undoc-members:
   :show-inheritance:
omg.net.base
============================================

.. automodule:: openmodelica_microgrid_gym.net.base
   :members:
   :undoc-members:
   :show-inheritance:
omg.util.itertools\_
========================================

.. automodule:: openmodelica_microgrid_gym.util.itertools_
   :members:
   :undoc-members:
   :show-inheritance:
omg.execution.runner
======================================

.. automodule:: openmodelica_microgrid_gym.execution.runner
   :members:
   :undoc-members:
   :show-inheritance:
omg.net.components
==================================================

.. automodule:: openmodelica_microgrid_gym.net.components
   :members:
   :undoc-members:
   :show-inheritance:
omg.agents.safeopt
====================================

.. automodule:: openmodelica_microgrid_gym.agents.safeopt
   :members:
   :undoc-members:
   :show-inheritance:
omg.util
=============================

Submodules
----------

.. toctree::

   omg.util.fastqueue
   omg.util.itertools_
   omg.util.transforms
   omg.util.randproc
   omg.util.recorder

Module contents
---------------

.. automodule:: openmodelica_microgrid_gym.util
   :members:
   :undoc-members:
   :show-inheritance:
omg.env.plotmanager
===================================================

.. automodule:: openmodelica_microgrid_gym.env.plotmanager
   :members:
   :undoc-members:
   :show-inheritance:
omg.aux\_ctl.observers
======================================================

.. automodule:: openmodelica_microgrid_gym.aux_ctl.observers
   :members:
   :undoc-members:
   :show-inheritance:
omg.env.plot
==================================

.. automodule:: openmodelica_microgrid_gym.env.plot
   :members:
   :undoc-members:
   :show-inheritance:
