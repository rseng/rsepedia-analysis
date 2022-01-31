[![CI general checks](https://github.com/tylunel/pvpumpingsystem/workflows/CI%20general%20checks/badge.svg)](https://github.com/tylunel/pvpumpingsystem/actions)
[![Coverage](https://codecov.io/gh/tylunel/pvpumpingsystem/branch/master/graph/badge.svg)](https://codecov.io/gh/tylunel/pvpumpingsystem)
[![Documentation Status](https://readthedocs.org/projects/pvpumpingsystem/badge/?version=latest)](https://pvpumpingsystem.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tylunel/pvpumpingsystem/master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02637/status.svg)](https://doi.org/10.21105/joss.02637)

![Logo](/docs/images/logo_pvpumpingsystem.jpg)
# pvpumpingsystem
*pvpumpingsystem* is a package providing tools for modeling and sizing
photovoltaic water pumping systems.

![Schema of a PV pumping system](/docs/images/schema_pvps.jpg)

It can model the whole functioning of such pumping system on an hourly basis
and eventually provide key financial and technical findings on a year.
Conversely it can help choose some elements of the pumping station
depending on output values wanted (like daily water consumption and
acceptable risk of water shortage). Further details are provided on the [scope page of the documentation](https://pvpumpingsystem.readthedocs.io/en/latest/package_overview.html).


# Documentation
The full package documentation is available on readthedocs:

[pvpumpingsystem docs](https://pvpumpingsystem.readthedocs.io/en/latest/?badge=latest)


# Installation
*pvpumpingsystem* works with Python 3.5 and superior only.

## With pip

For a rapid installation with pip, type in a command line interface:
```
python -m pip install pvpumpingsystem
```

Consult the docs for more information on installation:
https://pvpumpingsystem.readthedocs.io/en/latest/installation.html


# Hands-on start

Three examples of how the software can be used are in the folder
'docs/examples'.

For a given system, the first two show how to obtain the outflows,
probability of water shortage, life cycle cost and many other results:

[Basic usage example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_basic.ipynb)

[More advanced usage example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_advanced.ipynb)

The third shows how to optimize the selection of one or more component
on the pumping station based on user requirements:

[Sizing example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/sizing_example.ipynb)


# Contributions

All kind of contributions (documentation, testing, bug reports,
new features, suggestions...) are highly appreciated.
They can be reported as issues, pull requests, or simple message via
Github (prefered) or via mail of the maintainer.
﻿---
title: 'pvpumpingsystem: A Python package for modeling and sizing photovoltaic water pumping systems'
tags:
  - Python
  - sizing
  - modeling
  - water pumping
  - photovoltaics
  - solar energy
authors:
  - name: Tanguy R. Lunel
    orcid: 0000-0003-3008-1422
    affiliation: "1, 2"
  - name: Daniel R. Rousse
    orcid: 0000-0002-7247-5705
    affiliation: 1
affiliations:
 - name: Industrial Research Group In Technologies of Energy and Energy Efficiency (T3E), Department of Mechanical Sciences, Ecole de Technologie Supérieure Montreal
   index: 1
 - name: Department of Material Science and Engineering, Institut National des Sciences Appliquées Rennes
   index: 2
date: 14 April 2020
bibliography: paper.bib
---

# Summary

According to the World Health Organization, one tenth of the world population still lacks access to
basic water supply. One of the reasons for this is the remoteness of these populations from modern
water collection and distribution technologies, often coupled with an unfavorable socio-economic
situation. Photovoltaic (PV) pumping technology makes it possible to respond both to this problem
and to the criteria of sustainable development. However, these pumping systems must be carefully
modeled and sized in order to make the water supply cost efficient and reliable.

Pvpumpingsystem was conceived in order to tackle this issue. It is an open source package
providing various tools aimed at facilitating the modeling and the sizing of photovoltaic
powered water pumping systems. Even though the package is originally targeted at researchers
and engineers, three practical examples are provided in order to help anyone to use pvpumpingsystem.

Python is the programming language used in the software, and the code is structured with an
object-oriented approach. Continuous integration services allow for lint checking
and to test automation. Each class and function are documented with reference to the
literature when applicable. Pvpumpingsystem is released under a GPL-v3 license.

Pvpumpingsystem relies on already existing packages for photovoltaic and fluid mechanics modeling,
namely “pvlib-python” [@pvlib-python] and “fluids” [@fluids]. pvpumpingsystem's originality lies
in the implementation of various motor-pump models for finite power sources and in the coupling
of the distinct component models. In order to increase the understandability of the code,
each physical component of the PV pumping system corresponds to a class, like for example
the classes `Pump()`, `MPPT()`, `PipeNetwork()`, `Reservoir()`, and `PVGeneration()`. The previous objects
are then gathered in the class `PVPumpSystem()` which allows running a comprehensive model of
the pumping system.

The main inputs to the simulation are an hourly weather file, water source characteristics, expected water
consumption profile, and specifications of the photovoltaic array, motor-pump and water reservoir.
Typical outputs are hourly flow rates, unused electric power, efficiency of components, life
cycle cost and loss of load probability. The sizing module then builds on the
modeling tools, and uses them to provides functions to help choose
the best combination of components in order to minimize the total life cycle cost. Nevertheless,
sizing such complex systems is still an active field of research, and this module is subsequently
expected to be expanded with time.

Two software packages with similar scope already exist: PVsyst and online tool SISIFO, developed by
the MASLOWATEN consortium. Nevertheless, both are closed source, with restricted information
on models used internally, and no API is made accessible. Pvpumpingsystem also has the advantage
of providing ways to size PV pumping systems thanks to automation of pump and PV array choices.

Pvpumpingsystem is the second academic contribution of a broader research program on photovoltaic
water pumping launched in the Technologies of Energy and Energy Efficiency research group at Ecole de Technologie Supérieure Montreal, and is expected to grow with new
features and accuracy assessment provided by experimental studies. The authors also want to give
full access and help to anyone else interested in the use of the software.


# Acknowledgements

The authors would like to acknowledge Mr. Michel Trottier for his generous support to the T3E research group, as well as the NSERC and the FRQNT for their grants and subsidies.
The first author acknowledges the contributions and fruitful discussions with Louis Lamarche and Sergio Gualteros that inspired and helped with the current work.

# Reference
.. currentmodule:: pvpumpingsystem

API reference
=============

Classes
-------

The different classes of *pvpumpingsystem*.

.. autosummary::
   :toctree: generated/

   pvgeneration.PVGeneration
   mppt.MPPT
   pump.Pump
   pipenetwork.PipeNetwork
   reservoir.Reservoir
   consumption.Consumption
   pvpumpsystem.PVPumpSystem


Functions and methods
---------------------


Pump modeling
^^^^^^^^^^^^^

The core of the software's originality lies in the implementation of
different motor-pump models and in their coupling with the PV generator.

.. autosummary::
   :toctree: generated/
   
   pump.Pump.iv_curve_data
   pump.Pump.functIforVH
   pump.Pump.functIforVH_Arab
   pump.Pump.functIforVH_Kou
   pump.Pump.functIforVH_theoretical
   pump.Pump.functQforVH
   pump.Pump.functQforPH
   pump.Pump.functQforPH_Hamidat
   pump.Pump.functQforPH_Arab
   pump.Pump.functQforPH_Kou
   pump.Pump.functQforPH_theoretical
   pump.get_data_pump
   pump.specs_completeness
   pump._curves_coeffs_Arab06
   pump._curves_coeffs_Kou98
   pump._curves_coeffs_Hamidat08
   pump._curves_coeffs_theoretical
   pump._curves_coeffs_theoretical_variable_efficiency
   pump._curves_coeffs_theoretical_constant_efficiency
   pump._curves_coeffs_theoretical_basic
   pump._domain_V_H
   pump._domain_P_H
   pump._extrapolate_pow_eff_with_cst_efficiency
   pump.plot_Q_vs_P_H_3d
   pump.plot_I_vs_V_H_3d
   pump.plot_Q_vs_V_H_2d


Other components modeling
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   reservoir.Reservoir.change_water_volume

   consumption.adapt_to_flow_pumped

   pipenetwork.PipeNetwork.dynamichead

   pvgeneration.PVGeneration.run_model


Global modeling
^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   pvpumpsystem.PVPumpSystem.define_motorpump_model
   pvpumpsystem.PVPumpSystem.operating_point
   pvpumpsystem.PVPumpSystem.calc_flow
   pvpumpsystem.PVPumpSystem.calc_efficiency
   pvpumpsystem.PVPumpSystem.calc_reservoir
   pvpumpsystem.PVPumpSystem.run_model
   pvpumpsystem.function_i_from_v
   pvpumpsystem.operating_point
   pvpumpsystem.calc_flow_directly_coupled
   pvpumpsystem.calc_flow_mppt_coupled
   pvpumpsystem.calc_efficiency


Sizing tools
^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   sizing.shrink_weather_representative
   sizing.shrink_weather_worst_month
   sizing.subset_respecting_llp_direct
   sizing.size_nb_pv_direct
   sizing.subset_respecting_llp_mppt
   sizing.size_nb_pv_mppt
   sizing.sizing_minimize_npv


Ancillary functions
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   function_models.correlation_stats
   function_models.compound_polynomial_1_2
   function_models.compound_polynomial_1_3
   function_models.compound_polynomial_2_2
   function_models.compound_polynomial_2_3
   function_models.compound_polynomial_3_3
   function_models.polynomial_multivar_3_3_4
   function_models.polynomial_multivar_3_3_1
   function_models.polynomial_multivar_2_2_1
   function_models.polynomial_multivar_2_2_0
   function_models.polynomial_multivar_1_1_0
   function_models.polynomial_multivar_0_1_0
   function_models.polynomial_5
   function_models.polynomial_4
   function_models.polynomial_3
   function_models.polynomial_2
   function_models.polynomial_1
   function_models.polynomial_divided_2_1

   waterproperties.water_prop

   finance.initial_investment
   finance.net_present_value
.. _installation: pvpumpingsystem

Installation
============

Installing pvpumpingsystem can be done through different processes. Two of
them are detailled here, mainly thought for newcomers. Experienced users
can modify it to their liking.


    For people uncomfortable with package management, but who still plan on contributing or editing the code, follow the :ref:`anacondagit` instructions to install pvpumpingsystem along with Anaconda and Git.

    For people only interested in the use of the package, follow the :ref:`simple` instructions to install pvpumpingsystem alone.


Installing pvpumpingsystem is similar to installing most scientific python
packages, so in case of trouble see the :ref:`references` section
for further help.

Please see the :ref:`compatibility` section for information on the
optional packages that are needed for some pvpumpingsystem features.

.. _anacondagit:

Install pvpumpingsystem with Anaconda and Git
---------------------------------------------


- Anaconda:

The Anaconda distribution is an open source distribution providing Python
and others softwares and libraries useful for data science. Anaconda includes
many of the libraries needed for pvpumpingsystem (Pandas, NumPy, SciPy, etc).
Anaconda is especially recommended when using Windows.

Anaconda Python distribution is available at `<https://www.anaconda.com/download/>`_.

See `What is Anaconda? <https://www.anaconda.com/what-is-anaconda/>`_
and the `Anaconda Documentation <https://docs.anaconda.com/anaconda/>`_
for more information.


- Git:

Git is a version control system that widely help contribution and development
for open source softwares. Git should be native on most of Linux distribution,
but must be installed on Windows.

Git for Windows is available at `<https://gitforwindows.org/>`_.


- pvpumpingsystem:

Once you have Anaconda and git installed, open a command line interface
('Anaconda Prompt' on Windows, terminal in Linux and macOS), change
directory to the one you want to install pvpumpingsystem in, and type::

    pip install -e git+https://github.com/tylunel/pvpumpingsystem#egg=pvpumpingsystem



- Test pvpumpingsystem:

To ensure *pvpumpingsystem* and its dependencies are properly installed,
run the tests by going to the directory of pvpumpingsystem and by running
pytest::

    cd <relative/path/to/pvpumpingsystem/directory>
    pytest


.. _simple:

Install pvpumpingsystem alone
-----------------------------

.. note::

    Even if you decide not to use Anaconda or Git, you minimally need a Python
    version superior to 3.5, and to have pip and setuptools installed (installed
    by default with recent version of Python).

This second option simply uses pip::

    pip install pvpumpingsystem


If you have troubles with the use of pip, here is the
`pip documentation <https://pip.pypa.io/en/stable/user_guide/#installing-packages>`_
to help you.

To ensure *pvpumpingsystem* and its dependencies are properly installed,
you can consult the package information through pip::

    pip show pvpumpingsystem



.. _compatibility:

Compatibility
-------------

*pvpumpingsystem* is compatible with Python 3.5 and above.

Besides the libraries contained in Anaconda, *pvpumpingsystem* also requires:

* pvlib-python
* fluids
* numpy-financial

The full list of dependencies is detailled in
`setup.py <https://github.com/tylunel/pvpumpingsystem/docs/environment.rst>`_.


.. _references:

References
----------

.. note::

    This section was adapted from the pvlib-python documentation.
    Thanks to them for this useful listing!

Here are a few recommended references for installing Python packages:

* `Python Packaging Authority tutorial
  <https://packaging.python.org/tutorials/installing-packages/>`_
* `Conda User Guide
  <http://conda.pydata.org/docs/index.html>`_

Here are a few recommended references for git and GitHub:

* `The git documentation <https://git-scm.com/doc>`_:
  detailed explanations, videos, more links, and cheat sheets. Go here first!
* `Forking Projects <https://guides.github.com/activities/forking/>`_
* `Fork A Repo <https://help.github.com/articles/fork-a-repo/>`_
* `Cloning a repository
  <https://help.github.com/articles/cloning-a-repository/>`_


.. _getting_started:

Getting started
===============

To begin, a global view of how the code may be used is given here. Afterward,
it is recommended to go through the two first examples provided below
in order to get further step-by-step explanations.

General layout
--------------

The code make use of the possibilities provided by the object-oriented
paradigm of python. In particular, the code tries to match physical components
with objects (in the computer sense of the term) as much as possible.

Therefore, modeling a PV pumping system requires to start by defining each
object, i.e each component represented in below diagram.

.. image:: ../images/schema_pvps.jpg
  :width: 500
  :alt: Basic schema of a PV pumping system

Once all components have been declared in objects, they are gathered
in a parent object (an instance of class PVPumpSystem) where the type of
coupling between PV array and pump is declared.
Ultimately, the method ``run_model()`` launches the whole simulation.
It is also the step where the financial parameters can
be provided if a cost analysis is wanted.

Afterward, the results are contained in the object PVPumpSystem. At this point,
it is useful to have an IDE like Spyder or Pycharm to explore the values
internally contained. Otherwise, most of the results are actually stored in
the attributes ``flow``, ``efficiency``, ``water_stored``, ``npv``, ``llp``.

The examples below, in particular the jupyter notebook ones, provide further
details on each step of a standard simulation.


.. _examples:

Examples
--------

Three examples of how the software can be used are in the folder
``docs/examples``.
The examples are provided under two forms, as Jupyter Notebook files or as
Python files.


Jupyter Notebook
^^^^^^^^^^^^^^^^

Following examples can be run locally with Jupyter Notebook, or by clicking on the
corresponding icon in the upper right-hand corner of nbviewer pages, or by
accessing through the
`binder build <https://mybinder.org/v2/gh/tylunel/pvpumpingsystem/master>`_.


Simulation
""""""""""
The first two examples focus on simulation. These examples are important to understand because the modeling tools used here are the core of the software. These tools can be used later to get programs that fit a more particular use (for ex.: sizing process, parametric study, etc).
For a given system, the examples show how to obtain the values of interest for the user (output flow rates, total water pumped in a year, loss of load probability (llp), net present value (npv), efficiencies and others):

`Basic usage example <https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_basic.ipynb>`_

`More advanced usage example <https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_advanced.ipynb>`_

Once you went through these 2 examples, you are quite ready to dive into the code and adapt it to your needs.

Sizing
""""""
The third example shows how to use a sizing function written from the modeling
tools presented in the two examples above. This function aims at optimizing
the selection of the pump and the PV module, based on user requirements.

`Sizing example <https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/sizing_example.ipynb>`_


Python files
^^^^^^^^^^^^
These examples are also available in the form of python files in order to
freely adapt the code to your wishes. Directly check out in ``docs/examples``.
.. pvpumpingsystem documentation master file, created by sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to pvpumpingsystem's documentation!
===========================================

*pvpumpingsystem* is a package providing tools for modeling and sizing
offgrid photovoltaic water pumping systems. It is specially designed for
small to medium size systems, the type of pumping system typically used
for isolated communities.

It can model the whole functioning of such pumping system on an hourly basis
and eventually provide key financial and technical findings on a year.
Conversely it can help choose some elements of the pumping station
depending on output values wanted (like daily water consumption and
acceptable risk of water shortage). Find more on the scope of the software
in the section :ref:`package_overview`.

The source code for pvpumpingsystem is hosted on GitHub:
https://github.com/tylunel/pvpumpingsystem

The package was originally developped at T3E research group, in Ecole de
Technologie Superieure, Montreal, Qc, Canada, by Tanguy Lunel.

The software is published under the open source license GPL-v3.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   package_overview
   installation
   getting_started
   apireference



Citing pvpumpingsystem
======================
-----still to come------
Paper currently under review at Journal of Open Source Software.
.. _package_overview:

Package Overview
================

Introduction
------------

Scope
^^^^^
*Pvpumpingsystem* is an open source package providing various tools aimed
at facilitating the modeling and sizing of photovoltaic powered water
pumping systems.

This package helps users to model, test and validate different photovoltaic
pumping systems before actually installing it in situ. In order to guide the
designer in her/his choice, *pvpumpingsystem* provides both technical and
financial information on the system. Even though the package is originally
targeted at researchers and engineers, three practical examples are provided
in order to help anyone to use pvpumpingsystem

It models pumping systems minimally made of PV generator, DC motor-pump and
pipes. Each component can be precisely defined by the user in such a way
it corresponds closely to any actual system wanted.
User can choose to add a MPPT/DC-DC converter to increase the energy
yield of the PV array or to directly couple PV array and motor-pump.
The software also allows to add water tank to mitigate the effect of
intermittency.

.. image:: ../images/schema_pvps.jpg
  :width: 700
  :alt: Basic schema of a PV pumping system

The simulation eventually compute numerous outputs like hourly flow rates of
a given pump, efficiencies of components, risk of water shortage and
life cycle cost of the whole system.

.. image:: ../images/schema_simulation.jpg
  :width: 700
  :alt: Schema of a simulation on a PV pumping system

*Pvpumpingsystem* also offers to automate the process of sizing. In this case,
the user can provide a set of PV module, a set of motor-pumps and a
water needs file, and the software looks for the cheapest assembly while
making sure that it respects a minimum risk of water shortage.

.. image:: ../images/schema_sizing.jpg
  :width: 700
  :alt: Schema of standard sizing process

Nevertheless, the number of sizing processes can be infinite, and this module
is expected to significantly expand with time, welcoming new sizing process
based on different selection criteria or algorithms. In particular,
the reservoir size, the orientation of the PV array, the coupling strategy
or even the diameter of pipes are inputs that could ultimately become outputs of
the sizing process as well.


To better understand the possibilities of *pvpumpingsystem* and how it works,
you are invited to consult the examples available in the form of
Jupyter Notebook in :ref:`examples` or the corresponding python files in
``docs/examples``.



Code characteristics
^^^^^^^^^^^^^^^^^^^^

Python is the programming language used in the software, and the code is
structured within an object-oriented approach. Continuous integration
services allow checking for lint in the code and to automatize the tests.
Each class and function are documented in the docstring with reference to the
literature when applicable.

In *pvpumpingsystem*, in order to increase the understandability of the code,
the physical components of the PV pumping system corresponds to a class
when possible, like for example the classes Pump(), MPPT(), PipeNetwork(),
Reservoir() and PVGeneration().
Moreover, each of these classes are gathered into separate modules with
appropriate names (`pump.py`, `mppt.py`, etc).
The previous objects are then gathered in the class PVPumpSystem() which
allows running partial or comprehensive modeling of the pumping system.

A separate module `sizing.py` is dedicated to functions allowing to size these
systems. These functions are globally numerical methods, relying on numerous
simulations run according to an algorithm or to a factorial design.
`sizing.py` module can be expanded a lot as many strategies can be imagined to
size such a system.

Pvpumpingsystem relies on already existing packages for photovoltaic
and fluid mechanics modeling, namely *pvlib-python* and *fluids*.
*pvpumpingsystem*'s originality lies in the implementation of various
motor-pump models for finite power sources and in the coupling
of the distinct components models.

Pvpumpingsystem is released under a GPL-v3 license.


Databases accessible
--------------------

The PV module database of the California Energy Commission (CEC) is made
accessible through PVGeneration (being itself a wrapper of pvlib-python).
As this database is nearly comprehensive (more than 22,000 modules)
and regularly updated, it was considered that having a function to
define its own PV module was not relevant yet. Therefore, PV modules must
be declared by giving the reference in the corresponding attribute in
declaration of any PVGeneration instance.

Furthermore, the package also provide some pump and weather files
in the folder ``pvpumpingsystem/data``.

Concerning pump files, a template is provided in the folder in order to help
anyone fill it in with the specification of the pump they want to model.
A limited database coming from the company SunPumps is also accessible.
Nevertheless, it does not mean that the developers particularly encourage
their use, it rather reflects the difficulty to find other sources easily
accessible online. Any addition to the database is warmly welcomed here.

The weather files consist in a very restricted list of .epw files coming from
diverse climates and that users can exploit to learn and test the software.
Similar files for many location around the world are available at
`EnergyPlus website <https://energyplus.net/weather>`_, or can be constructed
using `PVGIS <https://re.jrc.ec.europa.eu/pvg_tools/en/#MR>`_.



Getting support and contribute
------------------------------

If you need help, you think you have discovered a bug, or if you would
like to edit *pvpumpingsystem*, then do not hesitate to open an issue on our
`GitHub issues page <https://github.com/tylunel/pvpumpingsystem/issues>`_
or on our
`GitHub pull request page <https://github.com/tylunel/pvpumpingsystem/pulls>`_.


Credits
-------
The T3E research group would like to acknowledge Mr. Michel Trottier for
his generous support, as well as the NSERC and the FRQNT for their grants
and subsidies. We also acknowledges the contributions and fruitful discussions
with Louis Lamarche and Sergio Gualteros that inspired and helped with the
current work.

pvpumpingsystem.pump.\_domain\_V\_H
===================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _domain_V_Hpvpumpingsystem.pump.Pump.functQforVH
=====================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforVHpvpumpingsystem.pump.Pump.functQforPH
=====================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforPHpvpumpingsystem.pvpumpsystem.PVPumpSystem.calc\_reservoir
=========================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.calc_reservoirpvpumpingsystem.consumption.adapt\_to\_flow\_pumped
===================================================

.. currentmodule:: pvpumpingsystem.consumption

.. autofunction:: adapt_to_flow_pumpedpvpumpingsystem.function\_models.compound\_polynomial\_1\_3
===========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: compound_polynomial_1_3pvpumpingsystem.pump.\_curves\_coeffs\_Arab06
=============================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_Arab06pvpumpingsystem.function\_models.polynomial\_divided\_2\_1
==========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_divided_2_1pvpumpingsystem.function\_models.compound\_polynomial\_3\_3
===========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: compound_polynomial_3_3pvpumpingsystem.pvpumpsystem.PVPumpSystem
=========================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autoclass:: PVPumpSystem

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~PVPumpSystem.__init__
      ~PVPumpSystem.calc_efficiency
      ~PVPumpSystem.calc_flow
      ~PVPumpSystem.calc_reservoir
      ~PVPumpSystem.define_motorpump_model
      ~PVPumpSystem.operating_point
      ~PVPumpSystem.run_model
   
   

   
   
   pvpumpingsystem.reservoir.Reservoir
===================================

.. currentmodule:: pvpumpingsystem.reservoir

.. autoclass:: Reservoir

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Reservoir.__init__
      ~Reservoir.change_water_volume
   
   

   
   
   pvpumpingsystem.finance.initial\_investment
===========================================

.. currentmodule:: pvpumpingsystem.finance

.. autofunction:: initial_investmentpvpumpingsystem.sizing.subset\_respecting\_llp\_mppt
====================================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: subset_respecting_llp_mpptpvpumpingsystem.pump.Pump.functQforPH\_Kou
==========================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforPH_Koupvpumpingsystem.pump.get\_data\_pump
====================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: get_data_pumppvpumpingsystem.pump.Pump
=========================

.. currentmodule:: pvpumpingsystem.pump

.. autoclass:: Pump

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Pump.__init__
      ~Pump.functIforVH
      ~Pump.functIforVH_Arab
      ~Pump.functIforVH_Kou
      ~Pump.functIforVH_theoretical
      ~Pump.functQforPH
      ~Pump.functQforPH_Arab
      ~Pump.functQforPH_Hamidat
      ~Pump.functQforPH_Kou
      ~Pump.functQforPH_theoretical
      ~Pump.functQforVH
      ~Pump.iv_curve_data
      ~Pump.starting_characteristics
   
   

   
   
   .. rubric:: Attributes

   .. autosummary::
   
      ~Pump.modeling_method
   
   pvpumpingsystem.pump.plot\_I\_vs\_V\_H\_3d
==========================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: plot_I_vs_V_H_3dpvpumpingsystem.pump.plot\_Q\_vs\_V\_H\_2d
==========================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: plot_Q_vs_V_H_2dpvpumpingsystem.pump.Pump.functQforPH\_theoretical
==================================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforPH_theoreticalpvpumpingsystem.function\_models.polynomial\_2
==============================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_2pvpumpingsystem.pipenetwork.PipeNetwork.dynamichead
===================================================

.. currentmodule:: pvpumpingsystem.pipenetwork

.. automethod:: PipeNetwork.dynamicheadpvpumpingsystem.pump.Pump.functIforVH\_theoretical
==================================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functIforVH_theoreticalpvpumpingsystem.function\_models.compound\_polynomial\_2\_3
===========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: compound_polynomial_2_3pvpumpingsystem.pvgeneration.PVGeneration
=========================================

.. currentmodule:: pvpumpingsystem.pvgeneration

.. autoclass:: PVGeneration

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~PVGeneration.__init__
      ~PVGeneration.run_model
   
   

   
   
   .. rubric:: Attributes

   .. autosummary::
   
      ~PVGeneration.pv_module_name
      ~PVGeneration.weather_data_and_metadata
   
   pvpumpingsystem.pump.\_curves\_coeffs\_theoretical
==================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_theoreticalpvpumpingsystem.function\_models.correlation\_stats
===================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: correlation_statspvpumpingsystem.pump.\_curves\_coeffs\_theoretical\_basic
=========================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_theoretical_basicpvpumpingsystem.pvpumpsystem.PVPumpSystem.operating\_point
==========================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.operating_pointpvpumpingsystem.pump.Pump.iv\_curve\_data
=========================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.iv_curve_datapvpumpingsystem.pvpumpsystem.calc\_flow\_mppt\_coupled
======================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autofunction:: calc_flow_mppt_coupledpvpumpingsystem.pump.Pump.functIforVH
=====================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functIforVHpvpumpingsystem.pvpumpsystem.function\_i\_from\_v
=================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autofunction:: function_i_from_vpvpumpingsystem.function\_models.polynomial\_multivar\_3\_3\_4
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_3_3_4pvpumpingsystem.pump.Pump.functQforPH\_Arab
===========================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforPH_Arabpvpumpingsystem.pump.\_curves\_coeffs\_Hamidat08
================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_Hamidat08pvpumpingsystem.function\_models.polynomial\_3
==============================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_3pvpumpingsystem.pump.\_curves\_coeffs\_theoretical\_constant\_efficiency
========================================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_theoretical_constant_efficiencypvpumpingsystem.function\_models.polynomial\_multivar\_1\_1\_0
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_1_1_0pvpumpingsystem.sizing.size\_nb\_pv\_direct
===========================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: size_nb_pv_directpvpumpingsystem.sizing.shrink\_weather\_worst\_month
====================================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: shrink_weather_worst_monthpvpumpingsystem.pump.specs\_completeness
========================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: specs_completenesspvpumpingsystem.function\_models.compound\_polynomial\_2\_2
===========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: compound_polynomial_2_2pvpumpingsystem.pvpumpsystem.calc\_efficiency
=============================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autofunction:: calc_efficiencypvpumpingsystem.sizing.subset\_respecting\_llp\_direct
======================================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: subset_respecting_llp_directpvpumpingsystem.pvpumpsystem.PVPumpSystem.run\_model
====================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.run_modelpvpumpingsystem.function\_models.polynomial\_1
==============================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_1pvpumpingsystem.pvpumpsystem.PVPumpSystem.calc\_flow
====================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.calc_flowpvpumpingsystem.pump.Pump.functQforPH\_Hamidat
==============================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functQforPH_Hamidatpvpumpingsystem.pvgeneration.PVGeneration.run\_model
====================================================

.. currentmodule:: pvpumpingsystem.pvgeneration

.. automethod:: PVGeneration.run_modelpvpumpingsystem.function\_models.polynomial\_multivar\_2\_2\_0
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_2_2_0pvpumpingsystem.mppt.MPPT
=========================

.. currentmodule:: pvpumpingsystem.mppt

.. autoclass:: MPPT

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~MPPT.__init__
   
   

   
   
   pvpumpingsystem.pump.\_extrapolate\_pow\_eff\_with\_cst\_efficiency
===================================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _extrapolate_pow_eff_with_cst_efficiencypvpumpingsystem.finance.net\_present\_value
===========================================

.. currentmodule:: pvpumpingsystem.finance

.. autofunction:: net_present_valuepvpumpingsystem.consumption.Consumption
=======================================

.. currentmodule:: pvpumpingsystem.consumption

.. autoclass:: Consumption

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Consumption.__init__
   
   

   
   
   pvpumpingsystem.pvpumpsystem.PVPumpSystem.define\_motorpump\_model
==================================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.define_motorpump_modelpvpumpingsystem.pvpumpsystem.operating\_point
=============================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autofunction:: operating_pointpvpumpingsystem.sizing.size\_nb\_pv\_mppt
=========================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: size_nb_pv_mpptpvpumpingsystem.function\_models.polynomial\_multivar\_2\_2\_1
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_2_2_1pvpumpingsystem.reservoir.Reservoir.change\_water\_volume
=========================================================

.. currentmodule:: pvpumpingsystem.reservoir

.. automethod:: Reservoir.change_water_volumepvpumpingsystem.waterproperties.water\_prop
===========================================

.. currentmodule:: pvpumpingsystem.waterproperties

.. autofunction:: water_proppvpumpingsystem.function\_models.compound\_polynomial\_1\_2
===========================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: compound_polynomial_1_2pvpumpingsystem.sizing.sizing\_minimize\_npv
============================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: sizing_minimize_npvpvpumpingsystem.function\_models.polynomial\_5
==============================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_5pvpumpingsystem.function\_models.polynomial\_4
==============================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_4pvpumpingsystem.function\_models.polynomial\_multivar\_0\_1\_0
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_0_1_0pvpumpingsystem.pump.Pump.functIforVH\_Kou
==========================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functIforVH_Koupvpumpingsystem.pvpumpsystem.calc\_flow\_directly\_coupled
==========================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. autofunction:: calc_flow_directly_coupledpvpumpingsystem.pipenetwork.PipeNetwork
=======================================

.. currentmodule:: pvpumpingsystem.pipenetwork

.. autoclass:: PipeNetwork

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~PipeNetwork.__init__
      ~PipeNetwork.dynamichead
   
   

   
   
   pvpumpingsystem.function\_models.polynomial\_multivar\_3\_3\_1
==============================================================

.. currentmodule:: pvpumpingsystem.function_models

.. autofunction:: polynomial_multivar_3_3_1pvpumpingsystem.pump.\_curves\_coeffs\_Kou98
============================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_Kou98pvpumpingsystem.sizing.shrink\_weather\_representative
======================================================

.. currentmodule:: pvpumpingsystem.sizing

.. autofunction:: shrink_weather_representativepvpumpingsystem.pvpumpsystem.PVPumpSystem.calc\_efficiency
==========================================================

.. currentmodule:: pvpumpingsystem.pvpumpsystem

.. automethod:: PVPumpSystem.calc_efficiencypvpumpingsystem.pump.Pump.functIforVH\_Arab
===========================================

.. currentmodule:: pvpumpingsystem.pump

.. automethod:: Pump.functIforVH_Arabpvpumpingsystem.pump.\_curves\_coeffs\_theoretical\_variable\_efficiency
========================================================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _curves_coeffs_theoretical_variable_efficiencypvpumpingsystem.pump.plot\_Q\_vs\_P\_H\_3d
==========================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: plot_Q_vs_P_H_3dpvpumpingsystem.pump.\_domain\_P\_H
===================================

.. currentmodule:: pvpumpingsystem.pump

.. autofunction:: _domain_P_H