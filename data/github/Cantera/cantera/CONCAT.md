# Cantera Code of Conduct

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
reported to the community leaders responsible for enforcement at
conduct@cantera.org.
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
# Contributing to Cantera

* For significant changes, please start a discussion on the Cantera
  Users' Group or create an issue on the [Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose) repository
  on GitHub to plan your modifications so that they can be implemented
  efficiently and in a way that doesn't conflict with any other planned
  future development
* Fork the `Cantera/cantera` repository on Github
* Clone your new repository or add it as a remote to an existing repository
* Check out the existing `main` branch, then start a new feature branch for
  your work
* When making changes, write code that is consistent with the surrounding code
  (see the [style guidelines](#style-guidelines) below)
* Add tests for any new features that you are implementing to either the
  GoogleTest-based test suite or the Python test suite.
* Add examples that highlight new capabilities, or update existing
  examples to make use of new features.
* As you make changes, commit them to your feature branch
  * Configure Git with your name and e-mail address before making any commits
  * Use descriptive commit messages (summary line of no more than 72 characters,
    followed by a blank line and a more detailed summary, if any)
  * Make related changes in a single commit, and unrelated changes in separate
    commits
  * Make sure that your commits do not include any undesired files, e.g., files
    produced as part of the build process or other temporary files.
  * Use Git's history-rewriting features (i.e., `git rebase -i`; see
    https://help.github.com/articles/about-git-rebase/) to organize your commits
    and squash "fixup" commits and reversions.
  * Do not merge your branch with `main`. If needed, you should rebase your branch
    onto the most recent `HEAD` commit of `main`.
  * Periodically run the test suite (`scons test`) to make sure that your
    changes are not causing any test failures.
* Push the changes on your new feature branch to your forked copy of the
  `Cantera/cantera` repository on GitHub.

* Submit a Pull Request on Github, from your forked copy. Check the results
  of the continuous-integration tests run using GitHub Actions and resolve
  any issues that arise.
* Additional discussion of good Git & Github workflow is provided at
  http://matplotlib.org/devel/gitwash/development_workflow.html and
  https://docs.scipy.org/doc/numpy-1.15.0/dev/gitwash/development_workflow.html
* Cantera is licensed under a [BSD
  license](https://github.com/Cantera/cantera/blob/main/License.txt) which
  allows others to freely modify the code, and if your Pull Request is accepted,
  then that code will be release under this license as well. The copyright for
  Cantera is held collectively by the contributors. If you have made a
  significant contribution, please add your name to the `AUTHORS` file.

# Style Guidelines

* Try to follow the style of surrounding code, and use variable names that
  follow existing patterns. Pay attention to indentation and spacing.
* Configure your editor to use 4 spaces per indentation level, and **never to
  use tabs**.
* Avoid introducing trailing whitespace
* Limit line lengths to 88 characters when possible
* Write comments to explain non-obvious operations
* Use whitespaces to improve code readability (examples: after commas; before and
  after mathematical operators (`+`/`-`/`*`/`/` except `^`), binary operators
  (`&&`/`||`/...), and comparisons (`<`/`>`/`==`/...); before and after equality
  signs `=` unless used for the assignment of a default parameter)
* Do not go out of your way to change formatting in otherwise unmodified code

## C++

* All classes, member variables, and methods should have Doxygen-style comments
  (e.g., comment lines starting with `//!` or comment blocks starting with `/*!`)
* Avoid defining non-trivial functions in header files
* Header files should include an 'include guard'
* Protected and private member variable names are generally prefixed with
  `m_`. For most classes, member variables should not be public.
* Class names use `InitialCapsNames`
* Methods use `camelCaseNames`
* Do not indent the contents of namespaces
* Code should follow the C++11 standard, with minimum required compiler versions
  GCC 4.8, Clang 3.4, MSVC 14.0 (2015) and Intel 15.0.
* Avoid manual memory management (i.e. `new` and `delete`), preferring to use
  standard library containers, as well as `std::unique_ptr` and
  `std::shared_ptr` when dynamic allocation is required.
* Portions of Boost which are "header only" may be used. If possible, include
  Boost header files only within .cpp files rather than other header files to
  avoid unnecessary increases in compilation time. Boost should not be added
  to the public interface unless its existence and use is optional. This keeps
  the number of dependencies low for users of Cantera. In these cases,
  `CANTERA_API_NO_BOOST` should be used to conditionally remove Boost dependencies.
* While Cantera does not specifically follow these rules, the following style
  guides are useful references for possible style choices and the rationales behind them.
  * The Google C++ Style Guide: https://google.github.io/styleguide/cppguide.html
  * http://geosoft.no/development/cppstyle.html
* For any new code, do *not* use the `doublereal` and `integer` typedefs for the
  basic types `double` and `int`, but also do not go out of your way to change
  uses of these in otherwise unmodified code.

## Python

* Style generally follows PEP8 (https://www.python.org/dev/peps/pep-0008/)
* Code in `.py` and `.pyx` files needs to be written to work with Python 3
* The minimum Python version that Cantera supports is Python 3.6, so code should only use features added in Python 3.6 or earlier
* Please use double quotes in all new Python code
# How to get support

> This project has a [Code of Conduct](https://github.com/Cantera/cantera/blob/main/CODE_OF_CONDUCT.md).
> By interacting with this repository, organisation, or community you agree to
> abide by its terms.

For **help**, **support** and **questions** please create a post on the
**[Cantera Users' Group](https://groups.google.com/group/cantera-users)**.
Any discussion of Cantera functionality such as how to use certain function
calls, syntax problems, input files, etc. should be directed to the Users' Group.

Further, the **[Cantera Gitter Chat](https://gitter.im/Cantera/Lobby)** is an
infrequently monitored chat room that can be used to discuss tangentially-related
topics such as how to model the underlying physics of a problem, share cool
applications that you have developed, etc.

Please **_do not_** raise an issue on GitHub unless it is a bug report or a
feature request. Issues that do not fall into these categories will be closed.
If you're not sure, please make a post on the
[Users' Group](https://groups.google.com/group/cantera-users) and someone will
be able to help you out.

## Documentation

The [documentation](https://cantera.org/documentation)
offers a number of starting points:

- [Python tutorial](https://cantera.org/tutorials/python-tutorial.html)
- [Application Examples in Python (Jupyter)](https://github.com/Cantera/cantera-jupyter#cantera-jupyter)
- [A guide to Cantera's input file format](https://cantera.org/tutorials/input-files.html)
- [Information about the Cantera community](https://cantera.org/community.html)

Documentation for the [development version of
Cantera](https://cantera.org/documentation/dev-docs.html) is also available.

## Contributions

See [`CONTRIBUTING.md`](https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md) on how to contribute.
<!-- Thanks for contributing code! Please include a description of your change and check your pull request against the list below. For further questions, refer to the contributing guide (https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md). -->

**Changes proposed in this pull request**

<!-- Provide a clear and concise description of changes and/or features introduced in this pull request. -->

-
-
-

**If applicable, fill in the issue number this pull request is fixing**

<!-- Issues with issue number '<issue>' are referenced as #<issue>. To link to an issue in the enhancements repository, use Cantera/enhancements#<issue>. -->

Closes #

**If applicable, provide an example illustrating new features this pull request is introducing**

<!-- A minimal, complete, and reproducible example demonstrating features introduced by this pull request. See https://stackoverflow.com/help/minimal-reproducible-example for additional suggestions on how to create such an example. -->

**Checklist**

- [ ] The pull request includes a clear description of this code change
- [ ] Commit messages have short titles and reference relevant issues
- [ ] Build passes (`scons build` & `scons test`) and unit tests address code coverage
- [ ] Style & formatting of contributed code follows [contributing guidelines](https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md)
- [ ] The pull request is ready for review
---
name: Bug report
about: Report reproducible software issues so we can improve
title: ''
labels: ''
assignees: ''
---

<!-- Please fill in the following information to report a problem with Cantera. If you have a question about using Cantera, please post it on our Google Users' Group (https://groups.google.com/forum/#!forum/cantera-users). Feature enhancements should be discussed in the dedicated Cantera enhancements repository (https://github.com/Cantera/enhancements/new/choose) -->

**Problem description**

<!-- A clear and concise description of what the bug is. -->

**Steps to reproduce**

<!-- A minimal, complete, and reproducible example demonstrating the problem. See https://stackoverflow.com/help/minimal-reproducible-example for additional suggestions on how to create such an example. -->

1. Open '...'
2. Run '....'
3. See error '....'

**Behavior**

<!-- Describe the result of executing the above steps, and how this differs from what you expect to happen. -->

**System information**

- Cantera version: [for example, 2.5.0 or the git commit hash]
- OS: [for example, Windows 10]
- Python/MATLAB/other software versions:

**Attachments**

<!-- If applicable, attach scripts and/or input files to help explain your problem. Please do *not* attach screenshots of code or terminal output. -->

**Additional context**

<!-- Add any other context about the problem here. -->
---
name: Feature request
about: Suggest a new feature to enhance Cantera's capabilities
title: ''
labels: ''
assignees: ''
---

Feature requests have been moved to
[Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose) and should be
opened there. Thank you for your suggestions!
.. Cantera

|cantera|

|doi| |codecov| |ci| |release|


What is Cantera?
================

Cantera is an open-source collection of object-oriented software tools for
problems involving chemical kinetics, thermodynamics, and transport processes.
Among other things, it can be used to:

* Evaluate thermodynamic and transport properties of mixtures
* Compute chemical equilibrium
* Evaluate species chemical production rates
* Conduct kinetics simulations with large reaction mechanisms
* Simulate one-dimensional flames
* Conduct reaction path analysis
* Create process simulations using networks of stirred reactors
* Model non-ideal fluids

Cantera can be used from Python and Matlab, or in applications written in C++
and Fortran 90. A number of `examples of Cantera's capabilities
<https://github.com/Cantera/cantera-jupyter>`_ are available in the form of
Jupyter notebooks. These examples can be tried interactively, in the cloud by
using the following MyBinder link:

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/repo/cantera/cantera-jupyter

Installation
============

`Installation instructions for the current release of Cantera
<https://cantera.org/install/index.html>`_ are available from the main `Cantera
documentation site <https://cantera.org>`_. Installers are provided for Windows (MSI
packages), macOS, and Ubuntu. Anaconda packages containing the Cantera Python and Matlab
modules are also available for Windows, macOS, and Linux.

.. image:: https://anaconda.org/cantera/cantera/badges/installer/conda.svg
    :target: https://anaconda.org/Cantera/cantera

For other platforms, or for users wishing to install a development version of
Cantera, `compilation instructions <https://cantera.org/install/index.html>`_
are also available.

Documentation
=============

The `documentation <https://cantera.org/documentation>`_
offers a number of starting points:

- `Python tutorial
  <https://cantera.org/tutorials/python-tutorial.html>`_
- `Application Examples in Python
  <https://github.com/Cantera/cantera-jupyter#cantera-jupyter>`_
- `A guide to Cantera's input file format
  <https://cantera.org/tutorials/input-files.html>`_
- `Information about the Cantera community
  <https://cantera.org/community.html>`_
- `Affiliated packages
  <https://cantera.org/affiliated-packages.html>`_

`Documentation for the development version of Cantera
<https://cantera.org/documentation/dev-docs.html>`_ is also available.

Code of Conduct
===============

.. image:: https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg?style=flat-square
    :alt: conduct
    :target: https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

In order to have a more open and welcoming community, Cantera adheres to a
`code of conduct <CODE_OF_CONDUCT.md>`_ adapted from the `Contributor Covenant
code of conduct <https://contributor-covenant.org/>`_.

Please adhere to this code of conduct in any interactions you have in the
Cantera community. It is strictly enforced on all official Cantera
repositories, websites, users' group, and other resources. If you encounter
someone violating these terms, please `contact the code of conduct team
<mailto:conduct@cantera.org>`_ (`@speth <https://github.com/speth>`_,
`@bryanwweber <https://github.com/bryanwweber>`_, and `@kyleniemeyer
<https://github.com/kyleniemeyer>`_) and we will address it as soon as
possible.

Development Site
================

The current development version is 2.6.0a4. The current stable version is
2.5.1. The `latest Cantera source code <https://github.com/Cantera/cantera>`_,
the `issue tracker <https://github.com/Cantera/cantera/issues>`_ for bugs and
enhancement requests, `downloads of Cantera releases and binary installers
<https://github.com/Cantera/cantera/releases>`_ , and the `Cantera wiki
<https://github.com/Cantera/cantera/wiki>`_ can all be found on Github.

Users' Group
============

The `Cantera Users' Group <https://groups.google.com/group/cantera-users>`_ is a
message board/mailing list for discussions amongst Cantera users.

Continuous Integration Status
=============================

|ci|


NumFOCUS
========

Cantera is a fiscally-sponsored project of `NumFOCUS <https://numfocus.org>`__,
a non-profit dedicated to supporting the open source scientific computing
community. Please consider `making a donation
<https://numfocus.salsalabs.org/donate-to-cantera/index.html>`__ to support the
development of Cantera through NumFOCUS.

.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: https://numfocus.salsalabs.org/donate-to-cantera/index.html
    :alt: Powered by NumFOCUS

.. |cantera| image:: https://cantera.org/assets/img/cantera-logo.png
    :target: https://cantera.org
    :alt: cantera logo
    :width: 675px
    :align: middle

.. |ci| image:: https://github.com/Cantera/cantera/workflows/CI/badge.svg
    :target: https://github.com/Cantera/cantera/actions?query=workflow%3ACI+event%3Apush

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.170284.svg
   :target: https://doi.org/10.5281/zenodo.1174508

.. |codecov| image:: https://img.shields.io/codecov/c/github/Cantera/cantera/main.svg
   :target: https://codecov.io/gh/Cantera/cantera?branch=main

.. |release| image:: https://img.shields.io/github/release/cantera/cantera.svg
   :target: https://github.com/Cantera/cantera/releases
   :alt: GitHub release
.. Cantera documentation master file, created by
   sphinx-quickstart on Mon Mar 12 11:43:09 2012.

Documentation
=============

These are the detailed API documentation pages for the Python and Matlab
interfaces for Cantera. There is also documentation of the CTI input file
format.

.. toctree::
   :maxdepth: 2

   yaml/index
   cti/classes
   cython/index
   matlab/index

*******************
CTI Class Reference
*******************

.. py:module:: cantera.ctml_writer

.. py:currentmodule:: cantera.ctml_writer

Basic Classes & Functions
=========================

.. autofunction:: units

.. autofunction:: validate

.. autoclass:: state
   :no-undoc-members:

Phases of Matter
================

.. autoclass:: phase
   :no-members:

.. autoclass:: ideal_gas
   :no-undoc-members:

.. autoclass:: stoichiometric_solid
   :no-members:

.. autoclass:: stoichiometric_liquid
   :no-undoc-members:

.. autoclass:: metal
   :no-undoc-members:

.. autoclass:: lattice
   :no-undoc-members:

.. autoclass:: lattice_solid
   :no-undoc-members:

.. autoclass:: liquid_vapor
   :no-undoc-members:

.. autoclass:: ideal_interface
   :no-undoc-members:

.. autoclass:: edge
   :no-undoc-members:

Elements and Species
====================

.. autoclass:: element
   :no-undoc-members:

.. autoclass:: species
   :no-undoc-members:

Thermodynamic Properties
========================

.. autoclass:: Mu0_table
   :no-undoc-members:

.. autoclass:: NASA
   :no-undoc-members:

.. autoclass:: NASA9
   :no-undoc-members:

.. autoclass:: Shomate
   :no-undoc-members:

.. autoclass:: const_cp
   :no-undoc-members:

Transport Properties
====================

.. autoclass:: gas_transport
   :no-undoc-members:

Reactions
=========

.. autoclass:: reaction
   :no-undoc-members:

.. autoclass:: Arrhenius
   :no-undoc-members:

.. autoclass:: three_body_reaction
   :no-undoc-members:

.. autoclass:: falloff_reaction
   :no-undoc-members:

.. autoclass:: chemically_activated_reaction
   :no-undoc-members:

.. autoclass:: pdep_arrhenius
   :no-undoc-members:

.. autoclass:: chebyshev_reaction
   :no-undoc-members:

.. autoclass:: surface_reaction
   :no-undoc-members:

.. autoclass:: edge_reaction
   :no-undoc-members:

.. autoclass:: stick
   :no-undoc-members:

Falloff Parameterizations
-------------------------

.. autoclass:: Troe
   :no-undoc-members:

.. autoclass:: SRI
   :no-undoc-members:

.. autoclass:: Lindemann
   :no-undoc-members:
.. py:currentmodule:: cantera

Chemical Kinetics
=================

.. contents::
   :local:

Kinetics Managers
-----------------

Kinetics
^^^^^^^^
.. autoclass:: Kinetics(infile='', phaseid='', phases=())

InterfaceKinetics
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceKinetics

Reactions
---------

These classes contain the definition of a single reaction, independent of a specific
`Kinetics` object. For legacy objects (CTI/XML input), each class integrates associated
rate expressions, whereas for the new, YAML-based implementation, reaction rate
evaluation is handled by dedicated `ReactionRate` objects.

Reaction
^^^^^^^^
.. autoclass:: Reaction(reactants='', products='')
   :no-undoc-members:

ElementaryReaction
^^^^^^^^^^^^^^^^^^
.. autoclass:: ElementaryReaction(reactants='', products='')
   :no-undoc-members:

ThreeBodyReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ThreeBodyReaction(reactants='', products='')
   :no-undoc-members:

FalloffReaction
^^^^^^^^^^^^^^^
.. autoclass:: FalloffReaction(reactants='', products='')
   :no-undoc-members:

ChemicallyActivatedReaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ChemicallyActivatedReaction(reactants='', products='')
   :no-undoc-members:

PlogReaction
^^^^^^^^^^^^
.. autoclass:: PlogReaction(reactants='', products='')
   :no-undoc-members:

ChebyshevReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ChebyshevReaction(reactants='', products='')
   :no-undoc-members:

BlowersMaselReaction
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselReaction(reactants='', products='')
   :no-undoc-members:

InterfaceReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceReaction(reactants='', products='')
   :no-undoc-members:

BlowersMaselInterfaceReaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselInterfaceReaction(reactants='', products='')
   :no-undoc-members:

Reaction Rates
---------------------------------------

ReactionRate
^^^^^^^^^^^^
.. autoclass:: ReactionRate()

ArrheniusRate
^^^^^^^^^^^^^
.. autoclass:: ArrheniusRate(A, b, Ea)
   :no-undoc-members:

BlowersMaselRate
^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselRate(A, b, Ea, w)
   :no-undoc-members:

FalloffRate
^^^^^^^^^^^
.. autoclass:: FalloffRate()
   :no-undoc-members:

LindemannRate
^^^^^^^^^^^^^
.. autoclass:: LindemannRate(low, high, falloff_coeffs)
   :no-undoc-members:

TroeRate
^^^^^^^^
.. autoclass:: TroeRate(low, high, falloff_coeffs)
   :no-undoc-members:

SriRate
^^^^^^^
.. autoclass:: SriRate(low, high, falloff_coeffs)
   :no-undoc-members:

TsangRate
^^^^^^^^^
.. autoclass:: TsangRate(low, high, falloff_coeffs)
   :no-undoc-members:

PlogRate
^^^^^^^^
.. autoclass:: PlogRate(rates)
   :no-undoc-members:

ChebyshevRate
^^^^^^^^^^^^^
.. autoclass:: ChebyshevRate(temperature_range, pressure_range, data)
   :no-undoc-members:

CustomRate
^^^^^^^^^^
.. autoclass:: CustomRate(k)
   :no-undoc-members:

Auxilliary Reaction Data (legacy only)
--------------------------------------

Arrhenius
^^^^^^^^^
.. autoclass:: Arrhenius(A, b, E)

Falloff
^^^^^^^
.. autoclass:: Falloff(params=(), init=True)
   :no-undoc-members:

TroeFalloff
^^^^^^^^^^^
.. autoclass:: TroeFalloff(params=(), init=True)
   :no-undoc-members:

SriFalloff
^^^^^^^^^^
.. autoclass:: SriFalloff(params=(), init=True)
   :no-undoc-members:

BlowersMasel
^^^^^^^^^^^^
.. autoclass:: BlowersMasel(A, b, E0, w)

Reaction Path Analysis
----------------------

ReactionPathDiagram
^^^^^^^^^^^^^^^^^^^
.. autoclass:: ReactionPathDiagram(Kinetics kin, str element)
.. py:currentmodule:: cantera

.. _sec-cython-onedim:

One-dimensional Reacting Flows
==============================

.. contents::
   :local:

Composite Domains
-----------------

FreeFlame
^^^^^^^^^
.. autoclass:: FreeFlame(gas, grid=None, width=None)

BurnerFlame
^^^^^^^^^^^
.. autoclass:: BurnerFlame(gas, grid=None, width=None)

CounterflowDiffusionFlame
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowDiffusionFlame(gas, grid=None, width=None)

CounterflowPremixedFlame
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowPremixedFlame(gas, grid=None, width=None)

ImpingingJet
^^^^^^^^^^^^
.. autoclass:: ImpingingJet(gas, grid=None, width=None)

IonFreeFlame
^^^^^^^^^^^^
.. autoclass:: IonFreeFlame(gas, grid=None, width=None)

   .. autoattribute:: E
   .. autoattribute:: electric_field_enabled
   .. automethod:: solve

IonBurnerFlame
^^^^^^^^^^^^^^
.. autoclass:: IonBurnerFlame(gas, grid=None, width=None)

   .. autoattribute:: E
   .. autoattribute:: electric_field_enabled
   .. automethod:: solve

Flow Domains
------------

IdealGasFlow
^^^^^^^^^^^^
.. autoclass:: IdealGasFlow(thermo)
    :inherited-members:

IonFlow
^^^^^^^
.. autoclass:: IonFlow(thermo)

FreeFlow
^^^^^^^^
.. autoclass:: FreeFlow(thermo)

AxisymmetricStagnationFlow
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AxisymmetricStagnationFlow(thermo)

Boundaries
----------

Inlet1D
^^^^^^^
.. autoclass:: Inlet1D(phase, name=None)

Outlet1D
^^^^^^^^
.. autoclass:: Outlet1D(phase, name=None)

OutletReservoir1D
^^^^^^^^^^^^^^^^^
.. autoclass:: OutletReservoir1D(phase, name=None)

SymmetryPlane1D
^^^^^^^^^^^^^^^
.. autoclass:: SymmetryPlane1D(phase, name=None)

Surface1D
^^^^^^^^^
.. autoclass:: Surface1D(phase, name=None)

ReactingSurface1D
^^^^^^^^^^^^^^^^^
.. autoclass:: ReactingSurface1D(phase, name=None)


Base Classes
------------
Domain1D
^^^^^^^^
.. autoclass:: Domain1D(name=None)

Boundary1D
^^^^^^^^^^
.. autoclass:: Boundary1D(phase, name=None)

Sim1D
^^^^^
.. autoclass:: Sim1D(domains)

FlameBase
^^^^^^^^^
.. autoclass:: FlameBase(domains, gas, grid=None)
.. py:currentmodule:: cantera

.. _sec-cython-zerodim:

Zero-Dimensional Reactor Networks
=================================

.. contents::
   :local:

Defining Functions
------------------

.. autoclass:: Func1

Base Classes
------------

ReactorBase
^^^^^^^^^^^
.. autoclass:: ReactorBase(contents=None, name=None)

FlowDevice
^^^^^^^^^^
.. autoclass:: FlowDevice(upstream, downstream, *, name=None)

Reactor Networks
----------------

.. autoclass:: ReactorNet(reactors=())

Reactors
--------

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(contents=None, name=None)

Reactor
^^^^^^^
.. autoclass:: Reactor(contents=None, *, name=None, energy='on')

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(contents=None, *, name=None, energy='on')

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(contents=None, *, name=None, energy='on')

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(contents=None, *, name=None, energy='on')

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(contents=None, *, name=None, energy='on')

ExtensibleReactor
^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleReactor

ExtensibleIdealGasReactor
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasReactor(contents=None, *, name=None, energy='on')

ExtensibleConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureReactor(contents=None, *, name=None, energy='on')

ExtensibleIdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureReactor(contents=None, *, name=None, energy='on')

Walls
-----

Wall
^^^^
.. autoclass:: Wall(left, right, *, name=None, A=None, K=None, U=None, Q=None, velocity=None, kinetics=(None,None))

WallSurface
^^^^^^^^^^^
.. autoclass:: WallSurface(wall, side)

Surfaces
--------

ReactorSurface
^^^^^^^^^^^^^^
.. autoclass:: ReactorSurface(kin=None, r=None, *, A=None)

Flow Controllers
----------------

MassFlowController
^^^^^^^^^^^^^^^^^^
.. autoclass:: MassFlowController(upstream, downstream, *, name=None, mdot=None)
   :inherited-members:

Valve
^^^^^
.. autoclass:: Valve(upstream, downstream, *, name=None, K=None)
   :inherited-members:

PressureController
^^^^^^^^^^^^^^^^^^
.. autoclass:: PressureController(upstream, downstream, *, name=None, master=None, K=None)
   :inherited-members:
.. py:currentmodule:: cantera

Objects Representing Phases
===========================

.. contents::
   :local:

Composite Phase Objects
-----------------------

These classes are composite representations of a substance which has
thermodynamic, chemical kinetic, and (optionally) transport properties.

.. autoclass:: Solution(infile='', phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), **kwargs)

.. autoclass:: Interface(infile='', phaseid='', phases=(), thermo=None, species=(), kinetics=None, reactions=())

.. autoclass:: DustyGas(infile, phaseid='')

Pure Fluid Phases
-----------------

The following convenience functions can be used to create `PureFluid` objects
with the indicated equation of state:

.. autofunction:: CarbonDioxide
.. autofunction:: Heptane
.. autofunction:: Hfc134a
.. autofunction:: Hydrogen
.. autofunction:: Methane
.. autofunction:: Nitrogen
.. autofunction:: Oxygen
.. autofunction:: Water

Representing Quantities of Phases
---------------------------------

.. autoclass:: Quantity

Representing Multiple States
----------------------------

.. autoclass:: SolutionArray

Utility Functions
-----------------

.. autofunction:: add_directory
.. autofunction:: get_data_directories
.. py:currentmodule:: cantera


Thermodynamic Properties
========================

.. contents::
   :local:

Phases
------

These classes are used to describe the thermodynamic state of a system.

ThermoPhase
^^^^^^^^^^^
.. autoclass:: ThermoPhase(infile='', phaseid='')

InterfacePhase
^^^^^^^^^^^^^^
.. autoclass:: InterfacePhase(infile='', phaseid='')

PureFluid
^^^^^^^^^
.. autoclass:: PureFluid(infile='', phaseid='')

Mixture
-------

.. autoclass:: Mixture

Species
-------

.. autoclass:: Species

Species Thermodynamic Properties
--------------------------------

These classes are used to describe the reference-state thermodynamic properties
of a pure species.

SpeciesThermo
^^^^^^^^^^^^^
.. autoclass:: SpeciesThermo(T_low, T_high, P_ref, coeffs)

ConstantCp
^^^^^^^^^^
.. autoclass:: ConstantCp(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Mu0Poly
^^^^^^^
.. autoclass:: Mu0Poly(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

NasaPoly2
^^^^^^^^^
.. autoclass:: NasaPoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Nasa9PolyMultiTempRegion
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: Nasa9PolyMultiTempRegion(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

ShomatePoly2
^^^^^^^^^^^^
.. autoclass:: ShomatePoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Element
-------

.. autoclass:: Element
    :no-undoc-members:

    .. autoattribute:: num_elements_defined
    .. autoattribute:: element_symbols
    .. autoattribute:: element_names
.. py:currentmodule:: cantera

Transport Properties
====================

.. autoclass:: Transport(infile='', phaseid='')
.. autoclass:: DustyGasTransport(infile='', phaseid='')

Species Transport Properties
----------------------------

.. autoclass:: GasTransportData(geometry='', diameter=-1, well_depth=-1, dipole=0.0, polarizability=0.0, rotational_relaxation=0.0, acentric_factor=0.0)
.. _sec-cython-documentation:

Python Module Documentation
===========================

Contents:

.. toctree::
   :maxdepth: 2

   importing
   thermo
   kinetics
   transport
   zerodim
   onedim
   constants
.. py:currentmodule:: cantera

Physical Constants
==================

These values are the same as those in the Cantera C++ header file ct_defs.h.

.. data:: avogadro

   Avogadro's Number, kmol\ :sup:`-1`

.. data:: gas_constant

   The ideal gas constant, J kmol\ :sup:`-1` K\ :sup:`-1`

.. data:: one_atm

   One atmosphere, Pa

.. data:: boltzmann

   Boltzmann constant, m\ :sup:`2` kg s\ :sup:`-2` K\ :sup:`-1`

.. data:: planck

   Planck constant, J s

.. data:: stefan_boltzmann

   The Stefan-Boltzmann constant, W m\ :sup:`-2` K\ :sup:`-4`

.. data:: electron_charge

   The charge on an electron, C

.. data:: electron_mass

   The mass of an electron, kg

.. data:: faraday

   Faraday constant, C kmol\ :sup:`-1`

.. data:: light_speed

   Speed of Light, m s\ :sup:`-1`

.. data:: permeability_0

   Permeability of free space, m kg s\ :sup:`-2` A\ :sup:`-2`

.. data:: epsilon_0

   Permittivity of free space, s\ :sup:`4` A\ :sup:`2` m\ :sup:`-3` kg\ :sup:`-1`
***********************
CTML to YAML conversion
***********************

.. py:module:: cantera.ctml2yaml
.. py:currentmodule:: cantera.ctml2yaml

The script ``ctml2yaml.py`` will convert files from the legacy CTML format to YAML input
format. The documentation below describes the classes and functions in the script. Each
function/method is annotated with the Python types that the function accepts.

Most users will access the functionality of this module via the command line with the
``ctml2yaml`` entry-point script. For programmatic access, the `main` and/or `convert`
functions should be used. `main` should be used when command line arguments must be
processed, while `convert` takes an input filename or a string containing the CTML file
to be converted, and optionally the name of the output file.

Module-level functions
======================

.. autofunction:: float2string
.. autofunction:: represent_float
.. autofunction:: get_float_or_quantity
.. autofunction:: split_species_value_string
.. autofunction:: clean_node_text
.. autofunction:: create_species_from_data_node
.. autofunction:: create_reactions_from_data_node
.. autofunction:: create_phases_from_data_node
.. autofunction:: convert
.. autofunction:: main

Conversion classes
==================

.. autoclass:: Phase
   :no-undoc-members:
.. autoclass:: Species
   :no-undoc-members:
.. autoclass:: SpeciesThermo
   :no-undoc-members:
.. autoclass:: SpeciesTransport
   :no-undoc-members:
.. autoclass:: Reaction
   :no-undoc-members:

Exceptions
==========

.. autoexception:: MissingXMLNode
.. autoexception:: MissingXMLAttribute
.. autoexception:: MissingNodeText
.. highlight:: yaml

.. _sec-yaml-reactions:

*********
Reactions
*********

The fields common to all ``reaction`` entries are:

``equation``
    The stoichiometric equation for the reaction. Each term (i.e.,
    stoichiometric coefficient, species name, ``+`` or ``<=>``) in the equation
    must be separated by a space.

    Reversible reactions may be written using ``<=>`` or ``=`` to separate
    reactants and products. Irreversible reactions are written using ``=>``.

``type``
    A string specifying the type of reaction or rate coefficient
    parameterization. The default is ``elementary``. Reaction types are:

    - :ref:`elementary <sec-yaml-elementary>`
    - :ref:`three-body <sec-yaml-three-body>`
    - :ref:`falloff <sec-yaml-falloff>`
    - :ref:`chemically-activated <sec-yaml-chemically-activated>`
    - :ref:`pressure-dependent-Arrhenius <sec-yaml-pressure-dependent-Arrhenius>`
    - :ref:`Chebyshev <sec-yaml-Chebyshev>`
    - :ref:`Blowers-Masel <sec-yaml-Blowers-Masel>`
    - :ref:`surface-Blowers-Masel <sec-yaml-surface-Blowers-Masel>`

    Reactions without a specified ``type`` on surfaces or edges are
    automatically treated as :ref:`interface <sec-yaml-interface-reaction>`
    reactions, and reactions that involve charge transfer between phases are
    automatically treated as :ref:`electrochemical <sec-yaml-electrochemical-reaction>`
    reactions. Reactions on surfaces or edges specifying ``type`` as
    ``Blowers-Masel`` are treated as
    :ref:`surface-Blowers-Masel <sec-yaml-surface-Blowers-Masel>`.


``duplicate``
    Boolean indicating whether the reaction is a known duplicate of another
    reaction. The default is ``false``.

``orders``
    An optional mapping of species to explicit reaction orders to use. Reaction
    orders for reactant species not explicitly mentioned are taken to be their
    respective stoichiometric coefficients. See
    `Reaction orders <https://cantera.org/science/reactions.html#reaction-orders>`__
    for additional information.

``negative-orders``
    Boolean indicating whether negative reaction orders are allowed. The
    default is ``false``.

``nonreactant-orders``
    Boolean indicating whether orders for non-reactant species are allowed.
    The default is ``false``.

Depending on the reaction ``type``, other fields may be necessary to specify
the rate of the reaction.

.. _sec-yaml-arrhenius:

Arrhenius expression
====================

Arrhenius expressions can be specified as either a three-element list containing
the pre-exponential factor :math:`A`, the temperature exponent :math:`b`, and
the activation energy :math:`E_a`, or a mapping containing the fields ``A``,
``b``, and ``Ea``. The following are equivalent::

    {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
    [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol]


.. _sec-yaml-efficiencies:

Efficiencies
============

Some reaction types include parameters for the "efficiency" of different species
as third-body colliders. For these reactions, the following additional fields
are supported:

``efficiencies``
    A mapping of species names to efficiency values

``default-efficiency``
    The efficiency for use for species not included in the ``efficiencies``
    mapping. Defaults to 1.0.


Reaction types
==============

.. _sec-yaml-elementary:

``elementary``
--------------

A homogeneous reaction with a pressure-independent rate coefficient and mass
action kinetics, as
`described here <https://cantera.org/science/reactions.html#reactions-with-a-pressure-independent-rate>`__.

Additional fields are:

``rate-constant``
    An :ref:`Arrhenius-type <sec-yaml-arrhenius>` list or mapping.

``negative-A``
    A boolean indicating whether a negative value for the pre-exponential factor
    is allowed. The default is ``false``.

Example::

    equation: N + NO <=> N2 + O
    rate-constant: {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
    negative-A: true


.. _sec-yaml-three-body:

``three-body``
--------------

A three body reaction as
`described here <https://cantera.org/science/reactions.html#three-body-reactions>`__.

The reaction equation should include the third body collision partner ``M``.

Includes the fields of an ``elementary`` reaction, plus the fields for
specifying :ref:`efficiencies <sec-yaml-efficiencies>`.

Example::

    equation: 2 O + M = O2 + M
    type: three-body
    rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0]
    efficiencies: {AR: 0.83, H2O: 5}


.. _sec-yaml-falloff:

``falloff``
-----------

A falloff reaction as
`described here <https://cantera.org/science/reactions.html#falloff-reactions>`__.

The reaction equation should include the pressure-dependent third body collision
partner ``(+M)`` or ``(+name)`` where ``name`` is the name of a species. The
latter case is equivalent to setting the efficiency for ``name`` to 1 and the
efficiency for all other species to 0.

Includes field for specifying :ref:`efficiencies <sec-yaml-efficiencies>` as well
as:

``high-P-rate-constant``
    An :ref:`sec-yaml-arrhenius` expression for the high-pressure limit

``low-P-rate-constant``
    An :ref:`sec-yaml-arrhenius` expression for the low-pressure limit

``Troe``
    Parameters for the
    `Troe <https://cantera.org/science/reactions.html#the-troe-falloff-function>`__
    falloff function. A mapping containing the keys ``A``, ``T3``, ``T1`` and
    optionally ``T2``. The default value for ``T2`` is 0.

``SRI``
    Parameters for the
    `SRI <https://cantera.org/science/reactions.html#the-sri-falloff-function>`__
    falloff function. A mapping containing the keys ``A``, ``B``, ``C``, and
    optionally ``D`` and ``E``. The default values for ``D`` and ``E`` are 1.0
    and 0.0, respectively.

Example::

    equation: H + CH2 (+ N2) <=> CH3 (+N2)
    type: falloff
    high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0]
    low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600}
    Troe: {A: 0.562, T3: 91, T1: 5836}


.. _sec-yaml-chemically-activated:

``chemically-activated``
------------------------

A chemically activated reaction as
`described here <https://cantera.org/science/reactions.html#chemically-activated-reactions>`__.

The parameters are the same as for :ref:`sec-yaml-falloff` reactions.

Example::

    equation: CH3 + OH (+M) <=> CH2O + H2 (+M)
    type: chemically-activated
    high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
    low-P-rate-constant: [282320.078, 1.46878, -3270.56495]

.. _sec-yaml-pressure-dependent-Arrhenius:

``pressure-dependent-Arrhenius``
--------------------------------

A pressure-dependent reaction using multiple Arrhenius expressions as
`described here <https://cantera.org/science/reactions.html#pressure-dependent-arrhenius-rate-expressions-p-log>`__.

The only additional field in this reaction type is:

``rate-constants``
    A list of mappings, where each mapping is the mapping form of an
    :ref:`sec-yaml-arrhenius` expression with the addition of a pressure ``P``.

Example::

    equation: H + CH4 <=> H2 + CH3
    type: pressure-dependent-Arrhenius
    rate-constants:
    - {P: 0.039474 atm, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}
    - {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}
    - {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}
    - {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}


.. _sec-yaml-Chebyshev:

``Chebyshev``
-------------

A reaction parameterized as a bivariate Chebyshev polynomial as
`described here <https://cantera.org/science/reactions.html#chebyshev-reaction-rate-expressions>`__.

Additional fields are:

``temperature-range``
    A list of two values specifying the minimum and maximum temperatures at
    which the rate constant is valid

``pressure-range``
    A list of two values specifying the minimum and maximum pressures at
    which the rate constant is valid

``data``
    A list of lists containing the Chebyshev coefficients

Example::

    equation: CH4 <=> CH3 + H
    type: Chebyshev
    temperature-range: [290, 3000]
    pressure-range: [0.0098692326671601278 atm, 98.692326671601279 atm]
    data: [[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],
           [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],
           [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],
           [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],
           [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],
           [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]]

.. _sec-yaml-Blowers-Masel:

``Blowers-Masel``
-----------------

A reaction with parameters to calculate rate constant based on Blowers Masel
approximation as `described here <https://cantera.org/science/reactions.html#sec-blowers-masel>`__.

Additional fields are:

``rate-constant``
    A list of values containing the pre-exponential factor :math:`A`, the
    temperature exponent :math:`b`, the intrinsic activation energy :math:`E_{a0}`,
    and the average of the bond dissociation energy of the bond breaking and that
    being formed in the reaction :math:`w`.

``negative-A``
    A boolean indicating whether a negative value for the pre-exponential factor
    is allowed. The default is ``false``.

Example::

    equation: O + H2 <=> H + OH
    type: Blowers-Masel
    rate-constant: {A: 3.87e+04 cm^2/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}


.. _sec-yaml-interface-reaction:

``interface``
-------------

A reaction occurring on a surface between two bulk phases, or along an edge
at the intersection of two surfaces, as
`described here <https://cantera.org/science/reactions.html#sec-surface>`__.

Includes the fields of an :ref:`sec-yaml-elementary` reaction plus:

``sticking-coefficient``
    An :ref:`Arrhenius-type <sec-yaml-arrhenius>` expression for the sticking coefficient

``Motz-Wise``
    A boolean applicable to sticking reactions, indicating whether to use the
    Motz-Wise correction factor for sticking coefficients near unity. Defaults
    to ``false``.

``sticking-species``
    The name of the sticking species. Required for sticking reactions only if
    the reaction includes multiple non-surface species.

``coverage-dependencies``
    A mapping of species names to coverage dependence parameters, where these
    parameters are contained in a mapping with the fields:

    ``a``
        Coefficient for exponential dependence on the coverage

    ``m``
        Power-law exponent of coverage dependence

    ``E``
        Activation energy dependence on coverage

Example::

    equation: 2 H(s) => H2 + 2 Pt(s)
    rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea: 67400 J/mol}
    coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}


.. _sec-yaml-electrochemical-reaction:

``electrochemical``
-------------------

Interface reactions involving charge transfer between phases,
as `described here <https://cantera.org/documentation/dev/doxygen/html/d6/ddd/classCantera_1_1ElectrochemicalReaction.html#details>`__.

Includes the fields of an :ref:`sec-yaml-interface-reaction` reaction, plus:

``beta``
    The symmetry factor for the reaction. Default is 0.5.

``exchange-current-density-formulation``
    Set to ``true`` if the rate constant parameterizes the exchange current
    density. Default is ``false``.

Example::

    equation: LiC6 <=> Li+(e) + C6
    rate-constant: [5.74, 0.0, 0.0]
    beta: 0.4

.. _sec-yaml-surface-Blowers-Masel:

``surface-Blowers-Masel``
-------------------------

A reaction occurring on a surface between two bulk phases, or along an edge
at the intersection of two surfaces, which the rate constant can be calculated
by Blowers Masel Approximation with Arrhenius expression as
`described here <https://cantera.org/science/reactions.html#surface-blowers-masel-reactions>`__.

Includes the fields of a :ref:`sec-yaml-Blowers-Masel` reaction and
the fields of an :ref:`sec-yaml-interface-reaction` reaction.

Example::

    equation: 2 H(s) => H2 + 2 Pt(s)
    type: Blowers-Masel
    rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea0: 67400 J/mol, w: 1000000 J/mol}
    coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}
.. highlight:: yaml

.. _sec-yaml-species:

*******
Species
*******

The fields of a ``species`` entry are:

``name``
    String identifier used for the species. Required.

``composition``
    Mapping that specifies the elemental composition of the species,
    e.g., ``{C: 1, H: 4}``. Required.

``thermo``
    Mapping containing the reference state thermodynamic model specification
    and parameters. See :ref:`sec-yaml-species-thermo`.

``equation-of-state``
    A mapping or list of mappings. Each mapping contains an equation of state
    model specification for the species, any parameters for that model, and any
    parameters for interactions with other species. See
    :ref:`sec-yaml-species-eos`. If this field is absent and a model is
    required, the ``ideal-gas`` model is assumed.

``critical-parameters``
    Mapping containing parameters related to the critical state of a species. Used in
    models that incorporate "real gas" effects, such as
    :ref:`sec-yaml-eos-redlich-kwong`.
    See :ref:`sec-yaml-species-crit-props`.

``transport``
    Mapping containing the species transport model specification and
    parameters. See :ref:`sec-yaml-species-transport`.

``sites``
    The number of sites occupied by a surface or edge species. Default is 1.

``Debye-Huckel``
    Additional model parameters used in the Debye-Hückel model. See
    :ref:`sec-yaml-Debye-Huckel` for more information.


.. _sec-yaml-species-thermo:

Species thermo models
=====================

Fields of a species ``thermo`` entry used by all models are:

``model``
    String specifying the model to be used. Required. Supported model strings
    are:

    - :ref:`NASA7 <sec-yaml-nasa7>`
    - :ref:`NASA9 <sec-yaml-nasa9>`
    - :ref:`Shomate <sec-yaml-shomate>`
    - :ref:`constant-cp <sec-yaml-constcp>`
    - :ref:`piecewise-Gibbs <sec-yaml-piecewise-gibbs>`

``reference-pressure``
    The reference pressure at which the given thermodynamic properties apply.
    Defaults to 1 atm.


.. _sec-yaml-nasa7:

NASA 7-coefficient polynomials
------------------------------

The polynomial form `described here <https://cantera.org/science/science-species.html#the-nasa-7-coefficient-polynomial-parameterization>`__,
given for one or two temperature regions. Additional fields of a ``NASA7``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    For one temperature region, this list contains the minimum and maximum
    temperatures for the polynomial. For two temperature regions, this list
    contains the minimum, intermediate, and maximum temperatures.

``data``
    A list with one item per temperature region, where that item is a 7 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order. Note that this is different from the standard CHEMKIN
    formulation that uses two temperature regions listed in descending order.

Example::

    thermo:
      model: NASA7
      temperature-ranges: [300.0, 1000.0, 5000.0]
      data:
      - [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09,
        -2.444854e-12, -1020.8999, 3.950372]
      - [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10,
        -6.753351e-15, -922.7977, 5.980528]


.. _sec-yaml-nasa9:

NASA 9-coefficient polynomials
------------------------------

The polynomial form `described here <https://cantera.org/science/science-species.html#the-nasa-9-coefficient-polynomial-parameterization>`__,
given for any number of temperature regions. Additional fields of a ``NASA9``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    This list contains the minimum temperature, the intermediate temperatures
    between each set pair of regions, and the maximum temperature.

``data``
    A list with one item per temperature region, where that item is a 9 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order.

Example::

    thermo:
      model: NASA9
      temperature-ranges: [200.00, 1000.00, 6000.0, 20000]
      reference-pressure: 1 bar
      data:
      - [2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03,
         1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02,
         -1.076003744E+01]
      - [5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04,
         1.491806679E-07,  -1.923105485E-11, 1.061954386E-15, 1.283210415E+04,
         -1.586640027E+01]
      - [8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02,
         2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06,
         -1.672099740E+03]

.. _sec-yaml-shomate:

Shomate polynomials
-------------------

The polynomial form `described here <https://cantera.org/science/science-species.html#the-shomate-parameterization>`__,
given for one or two temperature regions. Additional fields of a ``Shomate``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    For one temperature region, this list contains the minimum and maximum
    temperatures for the polynomial. For two temperature regions, this list
    contains the minimum, intermediate, and maximum temperatures.

``data``
    A list with one item per temperature region, where that item is a 7 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order.

Example::

    thermo:
      model: Shomate
      temperature-ranges: [298, 1300, 6000]
      data:
      - [25.56759, 6.096130, 4.054656, -2.671301, 0.131021,
        -118.0089, 227.3665]
      - [35.15070, 1.300095, -0.205921, 0.013550, -3.282780,
        -127.8375, 231.7120]


.. _sec-yaml-constcp:

Constant heat capacity
----------------------

The constant heat capacity model `described here <https://cantera.org/science/science-species.html#constant-heat-capacity>`__.
Additional fields of a ``constant-cp`` thermo entry are:

``T0``
    The reference temperature. Defaults to 298.15 K.

``h0``
    The molar enthalpy at the reference temperature. Defaults to 0.0.

``s0``
    The molar entropy at the reference temperature. Defaults to 0.0.

``cp0``
    The heat capacity at constant pressure. Defaults to 0.0.

``T-min``
    The minimum temperature at which this thermo data should be used.
    Defaults to 0.0.

``T-max``
    The maximum temperature at which this thermo data should be used.
    Defaults to infinity.

Example::

    thermo:
      model: constant-cp
      T0: 1000 K
      h0: 9.22 kcal/mol
      s0: -3.02 cal/mol/K
      cp0: 5.95 cal/mol/K

.. _sec-yaml-piecewise-gibbs:

Piecewise Gibbs
---------------

A model based on piecewise interpolation of the Gibbs free energy as
`described here <https://cantera.org/documentation/dev/doxygen/html/d4/d9e/classCantera_1_1Mu0Poly.html#details>`__
Additional fields of a ``piecewise-Gibbs`` entry are:

``h0``
    The molar enthalpy at the reference temperature of 298.15 K. Defaults to
    0.0.

``dimensionless``
    A boolean flag indicating whether the values of the Gibbs free energy are
    given in a dimensionless form, i.e., divided by :math:`RT`. Defaults to
    ``false``.

``data``
    A mapping of temperatures to values of the Gibbs free energy. The Gibbs free
    energy can be either in molar units (if ``dimensionless`` is ``false``) or
    nondimensionalized by the corresponding temperature (if ``dimensionless`` is
    ``true``). A value must be provided at :math:`T^\circ = 298.15` K.

``T-min``
    The minimum temperature at which this thermo data should be used.
    Defaults to 0.0.

``T-max``
    The maximum temperature at which this thermo data should be used.
    Defaults to infinity.

Example::

    thermo:
      model: piecewise-Gibbs
      h0: -230.015 kJ/mol
      dimensionless: true
      data: {298.15: -91.50963, 333.15: -85.0}


.. _sec-yaml-species-eos:

Species equation of state models
================================

``model``
    String specifying the model to be used. Required. Supported model strings
    are:

    - :ref:`constant-volume <sec-yaml-eos-constant-volume>`
    - :ref:`density-temperature-polynomial <sec-yaml-eos-density-temperature-polynomial>`
    - :ref:`HKFT <sec-yaml-eos-hkft>`
    - :ref:`ideal-gas <sec-yaml-eos-ideal-gas>`
    - :ref:`ions-from-neutral-molecule <sec-yaml-eos-ions-from-neutral>`
    - :ref:`liquid-water-IAPWS95 <sec-yaml-eos-liquid-water-iapws95>`
    - :ref:`molar-volume-temperature-polynomial <sec-yaml-eos-molar-volume-temperature-polynomial>`
    - :ref:`Redlich-Kwong <sec-yaml-eos-redlich-kwong>`

.. _sec-yaml-species-crit-props:

Species critical state parameters
=================================

``critical-temperature``
    The critical temperature of the species.

``critical-pressure``
    The critical pressure of the species.

``acentric-factor``
    Pitzer's acentric factor :math:`omega`.

.. _sec-yaml-eos-constant-volume:

Constant volume
---------------

A constant volume model as
`described here <https://cantera.org/documentation/dev/doxygen/html/da/d33/classCantera_1_1PDSS__ConstVol.html#details>`__.

Any one of the following may be specified:

``molar-volume``
    The molar volume of the species.

``molar-density``
    The molar density of the species.

``density``
    The mass density of the species.

Example::

    equation-of-state:
      model: constant-volume
      molar-volume: 1.3 cm^3/mol


.. _sec-yaml-eos-density-temperature-polynomial:

Density temperature polynomial
------------------------------

A model in which the density varies with temperature as
`described here <https://cantera.org/documentation/dev/doxygen/html/d0/d2f/classCantera_1_1PDSS__SSVol.html#details>`__.

Additional fields:

``data``
    Vector of 4 coefficients for a cubic polynomial in temperature

Example::

    equation-of-state:
      model: density-temperature-polynomial
      units: {mass: g, length: cm}
      data: [0.536504, -1.04279e-4, 3.84825e-9, -5.2853e-12]


.. _sec-yaml-eos-hkft:

HKFT
----

The Helgeson-Kirkham-Flowers-Tanger model as
`described here <https://cantera.org/documentation/dev/doxygen/html/d9/d18/classCantera_1_1PDSS__HKFT.html#details>`__.

Additional fields:

``h0``
    Enthalpy of formation at the reference temperature and pressure

``s0``
    Entropy of formation at the reference temperature and pressure

``a``
    4-element vector containing the coefficients :math:`a_1, \ldots , a_4`

``c``
    2-element vector containing the coefficients :math:`c_1` and :math:`c_2`

``omega``
    The :math:`\omega` parameter at the reference temperature and pressure

Example::

    equation-of-state:
      model: HKFT
      h0: -57433. cal/gmol
      s0: 13.96 cal/gmol/K
      a: [0.1839 cal/gmol/bar, -228.5 cal/gmol,
         3.256 cal*K/gmol/bar, -27260. cal*K/gmol]
      c: [18.18 cal/gmol/K, -29810. cal*K/gmol]
      omega: 33060 cal/gmol


.. _sec-yaml-eos-ideal-gas:

Ideal gas
---------

A species using the ideal gas equation of state, as
`described here <https://cantera.org/documentation/dev/doxygen/html/df/d31/classCantera_1_1PDSS__IdealGas.html#details>`__.
This model is the default if no ``equation-of-state`` section is included.


.. _sec-yaml-eos-ions-from-neutral:

Ions from neutral molecule
--------------------------

A species equation of state model used with the ``ions-from-neutral-molecule``
phase model, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d5/df4/classCantera_1_1PDSS__IonsFromNeutral.html#details>`__.

Additional fields:

``special-species``
    Boolean indicating whether the species is the "special species" in the
    phase. Default is ``false``.

``multipliers``
    A dictionary mapping species to neutral species multiplier values.

Example::

    equation-of-state:
      model: ions-from-neutral-molecule
      multipliers: {KCl(l): 1.2}


.. _sec-yaml-eos-liquid-water-iapws95:

Liquid Water IAPWS95
--------------------

A detailed equation of state for liquid water as
`described here <https://cantera.org/documentation/dev/doxygen/html/de/d64/classCantera_1_1PDSS__Water.html#details>`__.


.. _sec-yaml-eos-molar-volume-temperature-polynomial:

Molar volume temperature polynomial
-----------------------------------

A model in which the molar volume varies with temperature as
`described here <https://cantera.org/documentation/dev/doxygen/html/d0/d2f/classCantera_1_1PDSS__SSVol.html#details>`__.

Additional fields:

``data``
    Vector of 4 coefficients for a cubic polynomial in temperature


.. _sec-yaml-eos-redlich-kwong:

Redlich-Kwong
-------------

A model where species follow the Redlich-Kwong equation of state as
`described here <https://cantera.org/documentation/dev/doxygen/html/d6/d29/classCantera_1_1RedlichKwongMFTP.html#details>`__.

Additional fields:

``a``
    Pure-species ``a`` coefficient. Scalar or list of two values for a
    temperature-dependent expression.

``b``
    Pure-species ``b`` coefficient.

``binary-a``
    Mapping where the keys are species and the values are the ``a``
    coefficients for binary interactions between the two species.


.. _sec-yaml-species-transport:

Species transport models
========================

``model``
    String specifying the model type. The only model that is specifically
    handled is ``gas``.

Gas transport
-------------

Species transport properties are a rare exception to Cantera's use of SI units,
and use the units in which these properties are customarily reported. No
conversions are supported.

The additional fields of a ``gas`` transport entry are:

``geometry``
    A string specifying the geometry of the molecule. One of ``atom``,
    ``linear``, or ``nonlinear``.

``diameter``
    The Lennard-Jones collision diameter [Å]

``well-depth``
    The Lennard-Jones well depth [K]

``dipole``
    The permanent dipole moment [Debye]. Default 0.0.

``polarizability``
    The dipole polarizability [Å^3]. Default 0.0.

``rotational-relaxation``
    The rotational relaxation collision number at 298 K [-]. Default 0.0.

``acentric-factor``
    Pitzer's acentric factor [-]. Default 0.0. This value may also be specified as part
    of the :ref:`critical-parameters <sec-yaml-species-crit-props>` field, in which case
    the value provided there supersedes this one.

``dispersion-coefficient``
    The dispersion coefficient, normalized by :math:`e^2` [Å^5]. Default 0.0.

``quadrupole-polarizability``
    The quadrupole polarizability [Å^5]. Default 0.0.

Example::

    transport:
      model: gas
      geometry: linear
      well-depth: 107.4
      diameter: 3.458
      polarizability: 1.6
      rotational-relaxation: 3.8
***********************
CTI to YAML conversion
***********************

.. automodule:: cantera.cti2yaml
.. highlight:: yaml

*****************
General Structure
*****************

Sections
--------

The top level of a Cantera `YAML <https://yaml.org/spec/1.2/spec.html#Introduction>`__
input file is a mapping that defines different input file sections. Each
section consists of a list of mappings that define objects of the same type,
e.g., reactions, species, phases, or elements. The ``phases`` section of an input
file contains all of the phase definitions. Multiple sections containing
reaction, species, or element definitions can be used. The specific names
``reactions``, ``species``, and ``elements`` are used as defaults when looking
for :ref:`sec-yaml-reactions`, :ref:`sec-yaml-species`, and
:ref:`sec-yaml-elements` to add to a phase. A simple input file has the
following structure::

    phases:
    - name: spam
      thermo: ideal-gas
      # Additional fields come after this
    - name: green-eggs
      thermo: model-name
      # Additional fields come after this

    species:
    - name: A
      # Additional fields come after this
    - name: B
      # Additional fields come after this
    - name: C
      # Additional fields come after this

    reactions:
    - equation: A + B <=> C + D
      # Additional fields come after this
    - equation: A + C <=> 2 D
      # Additional fields come after this

Units
-----

While Cantera generally works internally in SI units, input values can be
provided using a number of different units.

Compound units are written using the asterisk (``*``) to indicate
multiplication, the forward slash (``/``) to indicate division, and the caret
(``^``) to indicate exponentiation. Exponents can include negative and decimal
values. Standard one-letter metric prefixes can be applied to any unit.
Supported base units are:

- Mass: ``g``
- Length: ``m``, ``micron``, ``angstrom``, ``Å``
- Time: ``s``, ``min``, ``hr``
- Temperature: ``K``, ``C``
- Current: ``A``
- Quantity: ``mol`` (gram mole), ``gmol``, ``mole``, ``kmol``, ``kgmol``, ``molec``

Supported compound units are:

- Energy: ``J``, ``cal``, ``erg``, ``eV``
- Activation Energy: ``K``, or any unit of energy per quantity (``J/kmol``,
  ``cal/mol``, etc.)
- Force: ``N``, ``dyn``
- Pressure: ``Pa``, ``atm``, ``bar``, ``dyn/cm^2``
- Volume: ``m^3``, ``liter``, ``L``, ``l``, ``cc``
- Other electrical units: ``ohm``, ``V``, ``coulomb``

Units can be specified on individual input values by placing them after the
value, separated by a space::

    {A: 1.45e9 cm^3/kmol, b: 0.4, Ea: 21033 kJ/kmol}

or by using a ``units`` mapping::

    units: {mass: g, quantity: mol, pressure: atm, activation-energy: cal/mol}

A ``units`` mapping will set the default units for all values within the same
YAML list or mapping, including any nested lists and mappings. Units not
specified by a mapping use the values from higher level mappings, or the Cantera
defaults if no ``units`` mapping specifies applicable units. If a ``units``
mapping appears in a list, it must be the first item in that list.

Default units may be set for ``mass``, ``length``, ``time``, ``temperature``,
``current``, ``quantity``, ``pressure``, ``energy``, and ``activation-energy``.
The units ``pressure`` and ``energy`` are used when these units appear
explicitly in the units that a value is being converted to within Cantera. For
example, a conversion to ``N/m^2`` will use the default units for mass, length,
and time, while a conversion to ``Pa`` will use the default units for pressure.

Conversions of activation energies implicitly include scaling by the gas
constant where necessary. Setting default units for ``energy`` and ``quantity``
will determine the default units of ``activation-energy``, which can be
overridden by explicitly giving the desired units of ``activation-energy``.
.. highlight:: yaml

.. _sec-yaml-elements:

********
Elements
********

``element`` entries are needed only when defining custom elements that are not
standard chemical elements, or defining specific isotopes.

The fields of an ``element`` entry are:

``symbol``
    The symbol used for the element, as used when specifying the composition of
    species.

``atomic-weight``
    The atomic weight of the element, in unified atomic mass units (dalton).

``atomic-number``
    The atomic number of the element. Optional.

``entropy298``
    The standard molar entropy of the element at 298.15 K. Optional.

*************************
YAML Input File Reference
*************************

.. toctree::
   :maxdepth: 2

   general
   phases
   elements
   species
   reactions
   cti_conversion
   ctml_conversion
.. highlight:: yaml

*****************
Phase Definitions
*****************

A ``phase`` is a mapping that contains definitions for the elements, species,
and optionally reactions that can take place in that phase. The fields of a
``phase`` entry are:

``name``
    String identifier used for the phase. Required.

``elements``
    Specification for the elements present in the phase. This can be:

    - Omitted, in which case the standard elements will be added as needed by
      the species included in the phase.
    - A list of element symbols, which can be either defined in the ``elements``
      section of the file or taken from the standard elements.
    - A list of single-key mappings of section names to lists of element
      symbols. These sections can be in the same file as the phase definition,
      or from another file if written as ``file-path/sectionname``. If a
      relative path is specified, the directory containing the current file is
      searched first, followed by the Cantera data path. Standard elements can
      be included by referencing the fictitious section ``default``.

``species``
    Specification for the species present in the phase. This can be:

    - a list of species that appear in the ``species`` section of the file.
    - The string ``all``, to indicate that all species in the ``species``
      section should be included. This is the default if no ``species`` entry
      is present.
    - A list of single-key mappings of section names to either the string
      ``all`` or a list of species names. These sections can be in the same
      file as the phase definition, or from another file if written as
      ``file-path/sectionname``. If a relative path is specified, the directory
      containing the current file is searched first, followed by the Cantera
      data path.

    Species may be skipped depending on the setting of the
    ``skip-undeclared-elements`` option.

``skip-undeclared-elements``
    If set to ``true``, do not add species that contain elements that are not
    explicitly included in the phase. The default is ``false``, where the
    presence of such species is considered an error.

``skip-undeclared-third-bodies``
   If set to ``true``, ignore third body efficiencies for species that are not
   defined in the phase. The default is ``false``, where the presence of
   such third body specifications is considered an error.

``state``
    A mapping specifying the thermodynamic state. See
    :ref:`sec-yaml-setting-state`.

``adjacent-phases``
    For interface phases, specification of adjacent phases that participate in reactions
    on the interface. This can be:

    - a list of phase names that appear in the ``phases`` section of the file.
    - A list of single-key mappings of section names to a list of phase names. These
      sections can be in the same file as the current phase definition, or from another
      file if written as ``file-path/section-name``. If a relative path is specified,
      the directory containing the current file is searched first, followed by the
      Cantera data path.

``thermo``
    String specifying the phase thermodynamic model to be used. Supported model
    strings are:

    - :ref:`binary-solution-tabulated <sec-yaml-binary-solution-tabulated>`
    - :ref:`compound-lattice <sec-yaml-compound-lattice>`
    - :ref:`constant-density <sec-yaml-constant-density>`
    - :ref:`Debye-Huckel <sec-yaml-Debye-Huckel>`
    - :ref:`edge <sec-yaml-edge>`
    - :ref:`fixed-stoichiometry <sec-yaml-fixed-stoichiometry>`
    - :ref:`HMW-electrolyte <sec-yaml-HMW-electrolyte>`
    - :ref:`ideal-gas <sec-yaml-ideal-gas>`
    - :ref:`ideal-molal-solution <sec-yaml-ideal-molal-solution>`
    - :ref:`ideal-condensed <sec-yaml-ideal-condensed>`
    - :ref:`ideal-solution-VPSS <sec-yaml-ideal-solution-VPSS>`
    - :ref:`ideal-surface <sec-yaml-ideal-surface>`
    - :ref:`ions-from-neutral-molecule <sec-yaml-ions-from-neutral-molecule>`
    - :ref:`lattice <sec-yaml-lattice>`
    - :ref:`liquid-water-IAPWS95 <sec-yaml-liquid-water-IAPWS95>`
    - :ref:`Margules <sec-yaml-Margules>`
    - :ref:`Maskell-solid-solution <sec-yaml-Maskell-solid-solution>`
    - :ref:`electron-cloud <sec-yaml-electron-cloud>`
    - :ref:`pure-fluid <sec-yaml-pure-fluid>`
    - :ref:`Redlich-Kister <sec-yaml-Redlich-Kister>`
    - :ref:`Redlich-Kwong <sec-yaml-Redlich-Kwong>`

``kinetics``
    String specifying the kinetics model to be used. Supported model strings
    are:

    - none
    - `gas <https://cantera.org/documentation/dev/doxygen/html/de/dae/classCantera_1_1GasKinetics.html#details>`__
    - `surface <https://cantera.org/documentation/dev/doxygen/html/d1/d72/classCantera_1_1InterfaceKinetics.html#details>`__
    - `edge <https://cantera.org/documentation/dev/doxygen/html/d0/df0/classCantera_1_1EdgeKinetics.html#details>`__

``reactions``
    Source of reactions to include in the phase, if a kinetics model has been
    specified. This can be:

    - The string ``all``, which indicates that all reactions from the
      ``reactions`` section of the file should be included. This is the default
      if no ``reactions`` entry is present.
    - The string ``declared-species``, which indicates that all reactions from
      the ``reactions`` section involving only species present in the phase
      should be included.
    - The string ``none``, which indicates that no reactions should be added.
      This can be used if reactions will be added programmatically after
      the phase is constructed.
    - A list of sections from which to include reactions. These sections can be
      in the same file as the phase definition, or from another file if written
      as ``file-path/sectionname``. If a relative path is specified, the
      directory containing the current file is searched first, followed by the
      Cantera data path.
    - A list of single-key mappings of section names to rules for adding
      reactions, where for each section name, that rule is either ``all`` or
      ``declared-species`` and is applied as described above.

``Motz-Wise``
    Boolean indicating whether the Motz-Wise correction should be applied to
    sticking reactions. Applicable only to interface phases. The default is
    ``false``. The value set at the phase level may be overridden on individual
    reactions.

``transport``
    String specifying the transport model to be used. Supported model strings
    are:

    - none
    - `high-pressure <https://cantera.org/documentation/dev/doxygen/html/d9/d63/classCantera_1_1HighPressureGasTransport.html#details>`__
    - `ionized-gas <https://cantera.org/documentation/dev/doxygen/html/d4/d65/classCantera_1_1IonGasTransport.html#details>`__
    - `mixture-averaged <https://cantera.org/documentation/dev/doxygen/html/d9/d17/classCantera_1_1MixTransport.html#details>`__
    - `mixture-averaged-CK <https://cantera.org/documentation/dev/doxygen/html/d9/d17/classCantera_1_1MixTransport.html#details>`__
    - `multicomponent <https://cantera.org/documentation/dev/doxygen/html/df/d7c/classCantera_1_1MultiTransport.html#details>`__
    - `multicomponent-CK <https://cantera.org/documentation/dev/doxygen/html/df/d7c/classCantera_1_1MultiTransport.html#details>`__
    - `unity-Lewis-number <https://cantera.org/documentation/dev/doxygen/html/d3/dd6/classCantera_1_1UnityLewisTransport.html#details>`__
    - `water <https://cantera.org/documentation/dev/doxygen/html/df/d1f/classCantera_1_1WaterTransport.html#details>`__



.. _sec-yaml-setting-state:

Setting the state
=================

The state of a ``phase`` can be set using two properties to set the
thermodynamic state, plus the composition.

The composition can be set using one of the following fields, depending on the
phase type. The composition is specified as a mapping of species names to
values. Where necessary, the values will be automatically normalized.

- ``mass-fractions`` or ``Y``
- ``mole-fractions`` or ``X``
- ``coverages``
- ``molalities`` or ``M``

The thermodynamic state can be set using the following property pairs, with some
exceptions for phases where setting that property pair is not implemented. All
properties are on a per unit mass basis where relevant:

- ``T`` and ``P``
- ``T`` and ``D``
- ``T`` and ``V``
- ``H`` and ``P``
- ``U`` and ``V``
- ``S`` and ``V``
- ``S`` and ``P``
- ``S`` and ``T``
- ``P`` and ``V``
- ``U`` and ``P``
- ``V`` and ``H``
- ``T`` and ``H``
- ``S`` and ``H``
- ``D`` and ``P``

The following synonyms are also implemented for use in any of the pairs:

- ``temperature``, ``T``
- ``pressure``, ``P``
- ``enthalpy``, ``H``
- ``entropy``, ``S``
- ``int-energy``, ``internal-energy``, ``U``
- ``specific-volume``, ``V``
- ``density``, ``D``


.. _sec-yaml-phase-thermo-models:

Phase thermodynamic models
==========================

.. _sec-yaml-binary-solution-tabulated:

``binary-solution-tabulated``
-----------------------------

A phase implementing tabulated standard state thermodynamics for one species in
a binary solution, as `described here <https://cantera.org/documentation/dev/doxygen/html/de/ddf/classCantera_1_1BinarySolutionTabulatedThermo.html#details>`__.

Includes the fields of :ref:`sec-yaml-ideal-molal-solution`, plus:

``tabulated-species``
    The name of the species to which the tabulated enthalpy and entropy is
    added.

``tabulated-thermo``
    A mapping containing three lists of equal lengths:

    ``mole-fractions``
        A list of mole fraction values for the tabulated species.

    ``enthalpy``
        The extra molar enthalpy to be added to the tabulated species at these
        mole fractions.

    ``entropy``
        The extra molar entropy to be added to the tabulated species at these
        mole fractions.


.. _sec-yaml-compound-lattice:

``compound-lattice``
--------------------

A phase that is comprised of a fixed additive combination of other lattice
phases, as `described here <https://cantera.org/documentation/dev/doxygen/html/de/de1/classCantera_1_1LatticeSolidPhase.html#details>`__.

Additional fields:

``composition``
    A mapping of component phase names to their relative stoichiometries.

Example::

    thermo: compound-lattice
    composition: {Li7Si3(s): 1.0, Li7Si3-interstitial: 1.0}


.. _sec-yaml-constant-density:

``constant-density``
--------------------

An incompressible phase with constant density, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d9/de4/classCantera_1_1ConstDensityThermo.html#details>`__.

Additional fields:

``density``
    The density of the phase

Example::

    thermo: constant-density
    density: 0.7 g/cm^3


.. _sec-yaml-Debye-Huckel:

``Debye-Huckel``
----------------

The Debye-Hückel model as
`described here <https://cantera.org/documentation/dev/doxygen/html/d8/d9a/classCantera_1_1DebyeHuckel.html#details>`__.

Additional parameters for this model are contained in the ``activity-data``
field:

``activity-data``
    The activity data field contains the following fields:

    ``model``
        One of ``dilute-limit``, ``B-dot-with-variable-a``,
        ``B-dot-with-common-a``, ``beta_ij``, or ``Pitzer-with-beta_ij``

    ``A_Debye``
        The value of the Debye "A" parameter, or the string ``variable`` to use
        a calculation based on the water equation of state.

    ``B_Debye``
        The Debye "B" parameter

    ``max-ionic-strength``
        The maximum ionic strength

    ``use-Helgeson-fixed-form``
        Boolean, ``true`` or ``false``

    ``default-ionic-radius``
        Ionic radius to use for species where the ionic radius has not been
        specified.

    ``B-dot``
        The value of B-dot.

    ``beta``
        List of mappings providing values of :math:`\beta_{ij}` for different
        species pairs. Each mapping contains a ``species`` key that contains a
        list of two species names, and a ``beta`` key that contains the
        corresponding value of :math:`\beta_{ij}`.

Example::

    thermo: Debye-Huckel
    activity-data:
      model: beta_ij
      max-ionic-strength: 3.0
      use-Helgeson-fixed-form: true
      default-ionic-radius: 3.042843 angstrom
      beta:
      - species: [H+, Cl-]
        beta: 0.27
      - species: [Na+, Cl-]
        beta: 0.15
      - species: [Na+, OH-]
        beta: 0.06

In addition, the Debye-Hückel model uses several species-specific properties
which may be defined in the ``Debye-Huckel`` field of the *species* entry. These
properties are:

``ionic-radius``
    Size of the species.

``electrolyte-species-type``
    One of ``solvent``, ``charged-species``, ``weak-acid-associated``,
    ``strong-acid-associated``, ``polar-neutral``, or ``nonpolar-neutral``.
    The type ``solvent`` is the default for the first species in the phase. The
    type ``charged-species`` is the default for species with a net charge.
    Otherwise, the default is and ``nonpolar-neutral``.

``weak-acid-charge``
    Charge to use for species that can break apart into charged species.

Example::

    name: NaCl(aq)
    composition: {Na: 1, Cl: 1}
    thermo:
      model: piecewise-Gibbs
      h0: -96.03E3 cal/mol
      dimensionless: true
      data: {298.15: -174.5057463, 333.15: -174.5057463}
    equation-of-state:
      model: constant-volume
      molar-volume: 1.3
    Debye-Huckel:
      ionic-radius: 4 angstrom
      electrolyte-species-type: weak-acid-associated
      weak-acid-charge: -1.0


.. _sec-yaml-edge:

``edge``
--------

A one-dimensional edge between two surfaces, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d9/d17/classCantera_1_1EdgePhase.html#details>`__.

Additional fields:

``site-density``
    The molar density of sites per unit length along the edge

Example::

    thermo: edge
    site-density: 5.0e-17 mol/cm


.. _sec-yaml-fixed-stoichiometry:

``fixed-stoichiometry``
-----------------------

A phase with fixed composition, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d3/d50/classCantera_1_1StoichSubstance.html#details>`__.


.. _sec-yaml-HMW-electrolyte:

``HMW-electrolyte``
-------------------

A dilute or concentrated liquid electrolyte phase that obeys the Pitzer
formulation for nonideality, as
`described here <https://cantera.org/documentation/dev/doxygen/html/de/d1d/classCantera_1_1HMWSoln.html#details>`__.

Additional parameters for this model are contained in the ``activity-data``
field:

``activity-data``
    The activity data field contains the following fields:

    ``temperature-model``
        The form of the Pitzer temperature model. One of ``constant``,
        ``linear`` or ``complex``.

    ``A_Debye``
        The value of the Debye "A" parameter, or the string ``variable`` to use
        a calculation based on the water equation of state.

    ``max-ionic-strength``
        The maximum ionic strength

    ``interactions``
        A list of mappings, where each mapping describes a binary or ternary
        interaction among species. Fields of this mapping include:

        ``species``
            A list of one to three species names

        ``beta0``
            The :math:`\beta^{(0)}` parameters for an cation/anion interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

        ``beta1``
            The :math:`\beta^{(1)}` parameters for an cation/anion interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

        ``beta2``
            The :math:`\beta^{(2)}` parameters for an cation/anion interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

        ``Cphi``
            The :math:`C^\phi` parameters for an cation/anion interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

        ``alpha1``
            The :math:`\alpha^{(1)}` parameter for an cation/anion interaction.

        ``alpha2``
            The :math:`\alpha^{(2)}` parameter for an cation/anion interaction.

        ``theta``
            The :math:`\theta` parameters for a like-charged binary interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

        ``lambda``
            The :math:`\lambda` parameters for binary interactions involving at
            least one neutral species. 1, 2, or 5 values depending on the value
            of ``temperature-model``.

        ``psi``
            The :math:`\Psi` parameters for ternary interactions involving three
            charged species. 1, 2, or 5 values depending on the value of
            ``temperature-model``.

        ``zeta``
            The :math:`\zeta` parameters for ternary interactions involving one
            neutral species. 1, 2, or 5 values depending on the value of
            ``temperature-model``.

        ``mu``
            The :math:`\mu` parameters for a neutral species self-interaction.
            1, 2, or 5 values depending on the value of ``temperature-model``.

    ``cropping-coefficients``

        ``ln_gamma_k_min``
            Default -5.0.

        ``ln_gamma_k_max``
            Default 15.0.

        ``ln_gamma_o_min``
            Default -6.0.

        ``ln_gamma_o_max``
            Default 3.0.

Example::

    thermo: HMW-electrolyte
    activity-data:
      temperature-model: complex
      A_Debye: 1.175930 kg^0.5/gmol^0.5
      interactions:
      - species: [Na+, Cl-]
        beta0: [0.0765, 0.008946, -3.3158E-6, -777.03, -4.4706]
        beta1: [0.2664, 6.1608E-5, 1.0715E-6, 0.0, 0.0]
        beta2: [0.0, 0.0, 0.0, 0.0, 0.0]
        Cphi: [0.00127, -4.655E-5, 0.0, 33.317, 0.09421]
        alpha1: 2.0
      - species: [H+, Cl-]
        beta0: [0.1775]
        beta1: [0.2945]
        beta2: [0.0]
        Cphi: [0.0008]
        alpha1: 2.0
      - species: [Na+, OH-]
        beta0: 0.0864
        beta1: 0.253
        beta2: 0.0
        Cphi: 0.0044
        alpha1: 2.0
        alpha2: 0.0
      - {species: [Cl-, OH-], theta: -0.05}
      - {species: [Na+, Cl-, OH-], psi: -0.006}
      - {species: [Na+, H+], theta: 0.036}
      - {species: [Cl-, Na+, H+], psi: [-0.004]}


.. _sec-yaml-ideal-gas:

``ideal-gas``
-------------

The ideal gas model as
`described here <https://cantera.org/documentation/dev/doxygen/html/d7/dfa/classCantera_1_1IdealGasPhase.html#details>`__.


.. _sec-yaml-ideal-molal-solution:

``ideal-molal-solution``
------------------------

A phase based on the mixing-rule assumption that all molality-based activity
coefficients are equal to one, as
`described here <https://cantera.org/documentation/dev/doxygen/html/da/d5c/classCantera_1_1IdealMolalSoln.html#details>`__.

Additional fields:

``standard-concentration-basis``
    A string specifying the basis for the standard concentration. One of
    ``unity``, ``species-molar-volume``, or ``solvent-molar-volume``.

``cutoff``
    Parameters for cutoff treatments of activity coefficients

    ``model``
        ``poly`` or ``polyExp``

    ``gamma_o``
        gamma_o value for the cutoff process at the zero solvent point

    ``gamma_k``
        gamma_k minimum for the cutoff process at the zero solvent point

    ``X_o``
        value of the solute mole fraction that centers the cutoff polynomials
        for the cutoff = 1 process

    ``c_0``
        Parameter in the polyExp cutoff treatment having to do with rate of
        exponential decay

    ``slope_f``
        Parameter in the ``polyExp`` cutoff treatment

    ``slope_g``
        Parameter in the ``polyExp`` cutoff treatment

Example::

    thermo: ideal-molal-solution
    standard-concentration-basis: solvent-molar-volume
    cutoff:
      model: polyexp
      gamma_o: 0.0001
      gamma_k: 10.0
      X_o: 0.2
      c_0: 0.05
      slope_f: 0.6
      slope_g: 0.0


.. _sec-yaml-ideal-condensed:

``ideal-condensed``
-------------------

A condensed phase ideal solution as
`described here <https://cantera.org/documentation/dev/doxygen/html/d3/d4c/classCantera_1_1IdealSolidSolnPhase.html#details>`__.

Additional fields:

``standard-concentration-basis``
    A string specifying the basis for the standard concentration. One of
    ``unity``, ``species-molar-volume``, or ``solvent-molar-volume``.


.. _sec-yaml-ideal-solution-VPSS:

``ideal-solution-VPSS``
-----------------------

An ideal solution model using variable pressure standard state methods as
`described here <https://cantera.org/documentation/dev/doxygen/html/dc/ddb/classCantera_1_1IdealSolnGasVPSS.html#details>`__.

Additional fields:

``standard-concentration-basis``
    A string specifying the basis for the standard concentration. One of
    ``unity``, ``species-molar-volume``, or ``solvent-molar-volume``.


.. _sec-yaml-ideal-surface:

``ideal-surface``
-----------------

An ideal surface phase, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d2/d95/classCantera_1_1SurfPhase.html#details>`__.

Additional fields:

``site-density``
    The molar density of surface sites


.. _sec-yaml-ions-from-neutral-molecule:

``ions-from-neutral-molecule``
------------------------------

A model that handles the specification of the chemical potentials for ionic
species, given a specification of the chemical potentials for the same phase
expressed in terms of combinations of the ionic species that represent neutral
molecules, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d7/d4a/classCantera_1_1IonsFromNeutralVPSSTP.html#details>`__.

Additional fields:

``neutral-phase``
    The ``name`` of the phase definition for the phase containing the neutral
    molecules.

Example::

    - name: KCl-ions
      thermo: ions-from-neutral-molecule
      neutral-phase: KCl-neutral
      species: [K+, Cl-]
    - name: KCl-neutral
      species: [KCl(l)]
      thermo: Margules


.. _sec-yaml-lattice:

``lattice``
-----------

A simple thermodynamic model for a bulk phase, assuming a lattice of solid
atoms, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d1/da0/classCantera_1_1LatticePhase.html#details>`__.

Additional fields:

``site-density``
    The molar density of lattice sites


.. _sec-yaml-liquid-water-IAPWS95:

``liquid-water-IAPWS95``
------------------------

An equation of state for liquid water, as
`described here <https://cantera.org/documentation/dev/doxygen/html/dc/d86/classCantera_1_1WaterSSTP.html#details>`__.


.. _sec-yaml-Margules:

``Margules``
------------

A phase employing the Margules approximation for the excess Gibbs free energy, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d7/dfe/classCantera_1_1MargulesVPSSTP.html#details>`__.

Additional fields:

``interactions``
    A list of mappings, where each mapping has the following fields:

    ``species``
        A list of two species names

    ``excess-enthalpy``
        A list of two values specifying the first and second excess enthalpy
        coefficients for the interaction of the specified species. Defaults to
        [0, 0].

    ``excess-entropy``
        A list of two values specifying the first and second excess entropy
        coefficients for the interaction of the specified species. Defaults to
        [0, 0].

    ``excess-volume-enthalpy``
        A list of two values specifying the first and second enthalpy
        coefficients for the excess volume interaction of the specified species.
        Defaults to [0, 0].

    ``excess-volume-entropy``
        A list of two values specifying the first and second entropy
        coefficients for the excess volume interaction of the specified species.
        Defaults to [0, 0].

Example::

  thermo: Margules
  interactions:
  - species: [KCl(l), LiCl(l)]
    excess-enthalpy: [-17570, -377]
    excess-entropy: [-7.627, 4.958]


.. _sec-yaml-Maskell-solid-solution:

``Maskell-solid-solution``
--------------------------

A condensed phase non-ideal solution with two species, as
`described here <https://cantera.org/documentation/dev/doxygen/html/dd/d3a/classCantera_1_1MaskellSolidSolnPhase.html#details>`__.

Additional fields:

``excess-enthalpy``
    The molar excess enthalpy

``product-species``
    String specifying the "product" species

Example::

    thermo: Maskell-solid-solution
    excess-enthalpy: 5 J/mol
    product-species: H(s)


.. _sec-yaml-electron-cloud:

``electron-cloud``
------------------

A phase representing an electron cloud, such as conduction electrons in a metal,
as `described here <https://cantera.org/documentation/dev/doxygen/html/d9/d13/classCantera_1_1MetalPhase.html#details>`__.

Additional fields:

``density``
    The density of the bulk metal


.. _sec-yaml-pure-fluid:

``pure-fluid``
--------------

A phase representing a pure fluid equation of state for one of several species,
as `described here <https://cantera.org/documentation/dev/doxygen/html/d1/d29/classCantera_1_1PureFluidPhase.html#details>`__.

Additional fields:

``pure-fluid-name``
    Name of the pure fluid model to use:
    - ``carbon-dioxide``
    - ``heptane``
    - ``HFC-134a``
    - ``hydrogen``
    - ``methane``
    - ``nitrogen``
    - ``oxygen``
    - ``water``


.. _sec-yaml-Redlich-Kister:

``Redlich-Kister``
------------------

A phase employing the Redlich-Kister approximation for the excess Gibbs free
energy, as
`described here <https://cantera.org/documentation/dev/doxygen/html/d0/d23/classCantera_1_1RedlichKisterVPSSTP.html#details>`__.

Additional fields:

``interactions``
    A list of mappings, where each mapping has the following fields:

    ``species``
        A list of two species names

    ``excess-enthalpy``
        A list of polynomial coefficients for the excess enthalpy of the
        specified binary interaction

    ``excess-entropy``
        A list of polynomial coefficients for the excess entropy of the
        specified binary interaction

Example::

  thermo: Redlich-Kister
  interactions:
  - species: [Li(C6), V(C6)]
    excess-enthalpy: [-3.268e+06, 3.955e+06, -4.573e+06, 6.147e+06, -3.339e+06,
                      1.117e+07, 2.997e+05, -4.866e+07, 1.362e+05, 1.373e+08,
                      -2.129e+07, -1.722e+08, 3.956e+07, 9.302e+07, -3.280e+07]
    excess-entropy: [0.0]


.. _sec-yaml-Redlich-Kwong:

``Redlich-Kwong``
-----------------

A multi-species Redlich-Kwong phase as
`described here <https://cantera.org/documentation/dev/doxygen/html/d6/d29/classCantera_1_1RedlichKwongMFTP.html#details>`__.

The parameters for each species are contained in the corresponding species
entries.

*****************************
Matlab Interface User's Guide
*****************************

.. toctree::
   :maxdepth: 2

   importing
   thermodynamics
   kinetics
   transport
   zero-dim
   one-dim
   data
   utilities
