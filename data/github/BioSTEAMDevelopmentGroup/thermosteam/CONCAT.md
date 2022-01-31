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
A script or a text file that runs the code to reproduce the error.

**Expected behavior**
A clear and concise description of what you expected to happen.

** Actual behavior**
Full trace back, as text or a screenshot

**Version**
 - PyPI  version or github

**Additional context**
Add any other context about the problem here.
---
title: "Thermosteam: BioSTEAM's Premier Thermodynamic Engine"
tags:
  - Python
  - BioSTEAM
  - thermodynamics
  - process modeling
  - mass and energy balance
  - chemical properties
  - mixture properties
  - phase equilibrium
authors:
  - name: Yoel Cortés-Peña
    orcid: 0000-0003-1742-5059
    affiliation: "1, 2"
affiliations:
  - name: Department of Civil and Environmental Engineering, University of Illinois at Urbana-Champaign
    index: 1
  - name: DOE Center for Advanced Bioenergy and Bioproducts Innovation (CABBI)
    index: 2
date: 12 May 2020
bibliography: paper.bib
---

# Summary

`Thermosteam` is a thermodynamic engine capable of solving mass and 
energy balances, estimating mixture properties, solving thermodynamic phase 
equilibria, and modeling stoichiometric reactions. All chemical data in 
`Thermosteam` is imported from the `chemicals` library [@chemicals], an 
open-source compilation of data and functions for the estimation of pure 
component chemical and mixture properties. `Thermosteam`'s fast and flexible
platform has enabled the evaluation of conceptual and emerging biochemical production
processes. The Biorefinery Simulation and Techno-Economic Analysis Modules 
(BioSTEAM) — capable of modeling reactors, distillation columns, heat exchangers, 
and other unit operations — has adopted `Thermosteam` as its premier 
thermodynamic engine [@BioSTEAM]. Published biorefinery designs modeled in 
BioSTEAM implement thermodynamic property packages created with `Thermosteam` 
[@Bioindustrial-Park], including a cornstover biorefinery for the production of 
cellulosic ethanol, a lipid-cane biorefinery for the co-production of ethanol 
and biodiesel, and a wheatstraw biorefinery for the production of cellulosic 
ethanol [@BioSTEAM;@Sanchis].

# Statement of Need

The overarching goal of `Thermosteam` is to aid the rigorous design and 
simulation of chemical production processes, whereby low value feedstocks are
converted to high value products via chemical reactions and thermodynamic-driven
separations. For example, modeling the separation of volatile chemicals from 
heavier ones in a distillation column (e.g., distilling ethanol from water),
requires vapor-liquid phase equilibrium calculations to predict how well volatile 
chemicals selectively partition into the vapor phase. Additionally, fluid 
viscosities, densities, and surface tensions are required to appropriately design 
a distillation column that can achieve a specified recovery of chemicals 
[@Perry]. 

Several open-source libraries in Python have comparable capabilities to 
`Thermosteam` in the estimation of fluid properties and phase equilibria: most 
notably `Cantera` and `CoolProp`. `Cantera` is a collection of software tools 
capable of modeling kinetic reactions, thermodynamic equilibrium, and mixture 
properties [@Cantera]. `Cantera`'s built-in chemicals are limited to 8, but new
chemicals can be defined by users with flexibility on the amount of detail. 
`Thermosteam` has yet to implement any features on kinetic reaction networks, 
but exposes a larger set of roughly 20,000 built-in chemicals from the `chemicals` 
library. Users may also define new models and pseudo-chemicals that are compatible 
with all of `Thermosteam`'s features. `CoolProp` offers fast and accurate 
thermodynamic and transport properties for 122 chemical components, and can 
estimate phase equilibrium and mixture properties [@CoolProp]. `CoolProp` 
also offers an interface to the NIST REFPROP software, which is considered the 
gold standard in thermophysical properties [@REFPROP]. It is within 
`Thermosteam`'s roadmap to use `CoolProp` as part of its built-in models. 
While `CoolProp` focuses on thermophysical chemical properties, `Thermosteam` 
also includes mass and energy balances and stoichiometric reactions as one of 
its central features.

# Roadmap

The main development  items in `Thermosteam`'s roadmap concerns the implementation
of fast, robust, and accurate algorithms for estimating mixture properties
and solving thermodynamic phase equilibria. Through `Thermosteam`, `BioSTEAM` is
able to evaluate a range of biofuels and bioproducts, but further efforts on these 
development items would enable the evaluation of a broader portfolio of 
potential bioproducts. 

In `Thermosteam`, Peng Robinson is the default equation of state 
for all pure components. However, the estimation of pure component chemical 
properties is not limited to solving the equation of state. Several models 
of thermodynamic properties (e.g., density, heat capacity, vapor pressure, 
heat of vaporization) are correlations that rely on fitted coefficients 
and key chemical properties (e.g., critical temperature and pressure). To 
facilitate the calculation of mixture properties, `Thermosteam`'s 
mixing rule estimates mixture properties by assuming a molar weighted average 
of the pure chemical properties. However, `Thermosteam` aims to implement rigorous 
equation of state (EOS) mixing rules for the estimation of mixture properties.

`Thermosteam` allows for fast estimation of thermodynamic equilibrium within 
hundreds of microseconds through the smart use of cache and Numba just-in-time 
(JIT) compiled functions [@numba]. The main vapor-liquid equilibrium (VLE) 
algorithm solves the modified Raoult’s law equation with activity coefficients
estimated through UNIQUAC Functional-group Activity Coefficients (UNIFAC) 
interaction parameters [@Gmehling]. Modified Raoult’s law is suitable to 
estimate VLE of nonideal mixtures under low to moderate pressures. At high to 
near-critical pressures, gaseous nonidealities become more significant. In a 
near future, `Thermosteam` may also implement the Predictive Soave–Redlich–Kwong
(PSRK) functional group model for estimating phase equilibrium of critical
mixtures. 

All of `Thermosteam`'s application program interface (API) is documented with 
examples. These examples also serve as preliminary tests that must pass before
accepting any changes to the software via continuous integration on Github.
Additionally, the online documentation includes a full tutorial that concludes 
with the creation of a property package. `Thermosteam`’s powerful features 
and extensive documentation encourage users to become a part of its
community-driven platform and help it become more industrially and academically 
relevant. 

# Acknowledgements

I would like to thank Caleb Bell for developing the open-source `chemicals` library
in Python, which has served as both groundwork and inspiration for developing `Thermosteam`. 
This material is based upon work supported by the National Science Foundation Graduate
Research Fellowship Program under Grant No. DGE—1746047. Any opinions, findings,
and conclusions or recommendations expressed in this publication are those of 
the authors and do not necessarily reflect the views of the National Science
Foundation. This work was funded by the DOE Center for Advanced Bioenergy and 
Bioproducts Innovation (U.S. Department of Energy, Office of Science, Office of 
Biological and Environmental Research under Award Number DE-SC0018420). Any 
opinions, findings, and conclusions or recommendations expressed in this 
publication are those of the author and do not necessarily reflect the views of
the U.S. Department of Energy.

# ReferencesThe Thermosteam Project Contributors is composed of:

* Yoel Cortes-Pena (Main Thermosteam author and maintainer)
* All other developers that have contributed to the thermosteam repository:
  
      https://github.com/BioSTEAMDevelopmentGroup/thermosteam/graphs/contributors

Additionally, some assets and code were originally sourced from third-party
authors or projects, including:

* Most chemical properties, including data and functions, are derived from `thermo <https://github.com/CalebBell/thermo>`_, by Caleb Bell.====================================================
Thermosteam: BioSTEAM's Premier Thermodynamic Engine 
====================================================

.. image:: http://img.shields.io/pypi/v/thermosteam.svg?style=flat
   :target: https://pypi.python.org/pypi/thermosteam
   :alt: Version_status
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://thermosteam.readthedocs.io/en/latest/
   :alt: Documentation
.. image:: https://img.shields.io/pypi/pyversions/thermosteam.svg
   :target: https://pypi.python.org/pypi/thermosteam
   :alt: Supported_versions
.. image:: http://img.shields.io/badge/license-UIUC-blue.svg?style=flat
   :target: https://github.com/BioSTEAMDevelopmentGroup/thermosteam/blob/master/LICENSE.txt
   :alt: license
.. image:: https://coveralls.io/repos/github/BioSTEAMDevelopmentGroup/thermosteam/badge.svg?branch=master
   :target: https://coveralls.io/github/BioSTEAMDevelopmentGroup/thermosteam?branch=master
   :alt: Coverage
.. image:: https://joss.theoj.org/papers/10.21105/joss.02814/status.svg
   :target: https://doi.org/10.21105/joss.02814

.. contents::

What is Thermosteam?
--------------------

Thermosteam is a standalone thermodynamic engine capable of estimating mixture 
properties, solving thermodynamic phase equilibria, and modeling stoichiometric 
reactions. Thermosteam builds upon `chemicals <https://github.com/CalebBell/chemicals>`_, 
the chemical properties component of the Chemical Engineering Design Library, 
with a robust and flexible framework that facilitates the creation of property packages.  
`The Biorefinery Simulation and Techno-Economic Analysis Modules (BioSTEAM) <https://biosteam.readthedocs.io/en/latest/>`_ 
is dependent on thermosteam for the simulation of unit operations.

Installation
------------

Get the latest version of Thermosteam from `PyPI <https://pypi.python.org/pypi/thermosteam/>`_.
If you have an installation of Python with pip, simple install it with::

    $ pip install thermosteam

To get the git version and install it, run::

    $ git clone --depth 100 git://github.com/BioSTEAMDevelopmentGroup/thermosteam
    $ cd thermosteam
    $ pip install .

We use the `depth` option to clone only the last 100 commits. Thermosteam has a 
long history, so cloning the whole repository (without using the depth option)
may take over 30 min.

If you would like to clone all branches, add the "--no-single-branch" flag as such::

    $ git clone --depth 100 --no-single-branch git://github.com/BioSTEAMDevelopmentGroup/thermosteam

Documentation
-------------

Thermosteam's documentation is available on the web:

    http://thermosteam.readthedocs.io/

Bug reports
-----------

To report bugs, please use the thermosteam's Bug Tracker at:

    https://github.com/BioSTEAMDevelopmentGroup/thermosteam


License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the thermosteam license, if it is convenient for you,
please cite Thermosteam if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


Citation
--------

To cite Thermosteam in publications use::

    Cortes-Pena, Y., (2020). Thermosteam: BioSTEAM's Premier Thermodynamic Engine. 
    Journal of Open Source Software, 5(56), 2814. doi.org/10.21105/joss.02814
Contributing to Thermosteam
===========================

General Process
---------------

Here’s the short summary of how to contribute using git bash:

#. If you are a first-time contributor:

   * Go to https://github.com/BioSTEAMDevelopmentGroup/thermosteam and click the “fork” button to create your own copy of the project.

   * Clone the project to your local computer::
    
        git clone --depth 100 https://github.com/your-username/thermosteam.git
    
   * Change the directory::
    
        cd thermosteam
    
   * Add the upstream repository::
    
        git remote add upstream https://github.com/BioSTEAMDevelopmentGroup/thermosteam.git
    
   * Now, git remote -v will show two remote repositories named "upstream" (which refers to the thermosteam repository), and "origin" (which refers to your personal fork).

#. Develop your contribution:

   * Pull the latest changes from upstream::

       git checkout master
       git pull upstream master

   * Create a branch for the feature you want to work on. Since the branch name will appear in the merge message, use a sensible name such as "Chemical-properties-enhancement"::

       git checkout -b Chemical-properties-enhancement

   * Commit locally as you progress (git add and git commit) Use a properly formatted commit message, write tests that fail before your change and pass afterward, run all the tests locally. Be sure to document any changed behavior in docstrings, keeping to the NumPy docstring standard.

#. To submit your contribution:

   * Push your changes back to your fork on GitHub::

       git push origin Chemical-properties-enhancement

   * Enter your GitHub username and password (repeat contributors or advanced users can remove this step by connecting to GitHub with SSH).

   * Go to GitHub. The new branch will show up with a green Pull Request button. Make sure the title and message are clear, concise, and self- explanatory. Then click the button to submit it.

   * If your commit introduces a new feature or changes functionality, post in https://github.com/BioSTEAMDevelopmentGroup/thermosteam/issues to explain your changes. For bug fixes, documentation updates, etc., this is generally not necessary, though if you do not get any reaction, do feel free to ask for a review.

Testing
-------

First install the developer version of thermosteam:

.. code-block:: bash

   $ cd thermosteam
   $ pip install -e .[dev]

This installs `pytest <https://docs.pytest.org/en/stable/>`__ and other
dependencies you need to run the tests locally. You can run tests by going
to your local thermosteam directory and running the following:

.. code-block:: bash
    
   $ pytest
    
This runs all the `doctests <https://docs.python.org/3.6/library/doctest.html>`__
in thermosteam, which covers most of the API. If any test is marked with a 
letter F, that test has failed. Pytest will point you to the location of the 
error, the values that were expected, and the values that were generated.


Documentation
-------------

Concise and thorough documentation is required for any contribution. Make sure to:

* Use NumPy style docstrings.
* Document all functions and classes.
* Document short functions in one line if possible.
* Mention and reference any equations or methods used and make sure to include the chapter and page number if it is a book or a long document.
* Preview the docs before making a pull request (open your cmd/terminal in the "docs" folder, run "make html", and open "docs/_build/html/index.html").


Authorship
----------

Authorship must be acknowledged for anyone contributing code, significant 
expertise, and/or other impactful efforts. Additionally, authorship should be 
included at the module-level, with a short description of the general contribution. 

If any code or implementation was copied from a third party, it should be rightfully
noted in the module-level documentation along with the corresponding copyright.

Any third-party code copied to the BioSTEAM software must be strictly open-source 
(not copy-left nor open-access). Additionally, if the license is different, 
the module should add the third-party license as an option (dual licensing is OK).


Best practices
--------------

Please refer to the following guides for best practices to make software designs more understandable, flexible, and maintainable:
    
* `PEP 8 style guide <https://www.python.org/dev/peps/pep-0008/>`__.
* `PEP 257 docstring guide <https://www.python.org/dev/peps/pep-0257/>`__.
* `Zen of Python philosophy <https://www.python.org/dev/peps/pep-0020/>`__.
* `SOLID programing principles <https://en.wikipedia.org/wiki/SOLID>`__.
