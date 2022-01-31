---
title: 'OpenSCM Two Layer Model: A Python implementation of the two-layer climate model'
tags:
  - Python
  - climate science
  - temperature projections
  - simple climate model
  - energy balance
  - reduced complexity climate model
authors:
  - name: Zebedee Nicholls
    orcid: 0000-0002-4767-2723
    affiliation: "1, 2, 3"
  - name: Jared Lewis
    orcid: 0000-0002-8155-8924
    affiliation: "1, 2, 3"
affiliations:
 - name: Australian-German Climate & Energy College, The University of Melbourne, Parkville, Victoria, Australia
   index: 1
 - name: School of Geography, Earth and Atmosphere Sciences, The University of Melbourne, Parkville, Victoria, Australia
   index: 2
 - name: Climate Resource, Northcote, Victoria, Australia
   index: 3
date: 7 October 2020
bibliography: paper.bib
---

# Summary

The evolution of the climate is controlled by highly complex physical dynamics.
However, simplified representations are surprisingly powerful tools for understanding these dynamics [@held_2010_two_layer] and making climate projections [@meinshausen_2011_rcp].
The field of simple climate modelling is widely used, in particular for assessing the climatic implications of large numbers of different emissions scenarios, a task that cannot be performed with more complex models because of computational constraints.

One of the most commonly used models of the climate's response to changes in the "Earth's energy balance"
(energy input compared to energy output of the earth system) is the two-layer model originally introduced by @held_2010_two_layer.
While this model must be given energy imbalances (more precisely, radiative forcing) rather than emissions, it is nonetheless a widely used tool within climate science.
Approximately speaking, the model represents the Earth System as an ocean with two-layers.
The upper layer absorbs excess heat in the Earth System and then exchanges heat with the deep layer.
As a result, in response to a perturbation, the model responds with a distinctive two-timescale response, commonly referred to as the "fast" and "slow" warming components.
Since @held_2010_two_layer, the model has been extended to include updated representations of the efficiency of ocean heat uptake [@geoffroy_2013_two_layer2] as well as a state-dependent response to radiative forcing [@bloch_johnson_2015_feedback_dependence; @rohrschneider_2019_simple].

There are many simple climate models in the scientific literature [@rcmip_phase_1].
Given the context of this paper, we provide below a table of openly accessible models, their programming language, and their approach.
These models are conceptually similar to the two-layer model implemented here except they use different parameterisations for ocean heat uptake and the relationship between ocean heat uptake and warming.
On top of the relationship between ocean heat uptake and warming, these models also implement many other components of the climate system, e.g., carbon cycle, methane cycle, and the relationship between changes in atmospheric greenhouse gas concentrations and atmospheric energy fluxes.
The exception is the FaIR model [@smith_2018_fairv1_3], which uses the two-layer model as its thermal core.

OpenSCM Two Layer Model is an object-oriented and open-source implementation of the two-layer model.
It is written in Python, a user-friendly open-source language that is popular in the climate sciences, and uses the Pint package [@pint], a widely used units library, for unit handling.
It provides an extensible interface for the two-layer model, which could then be coupled with other modules as researchers see fit.
The implementation also provides an easy way to convert between the two-layer model of @held_2010_two_layer and the mathematically equivalent two-timescale impulse response model, used most notably as the thermal core of the FaIR model [@smith_2018_fairv1_3].
The conversion between the two is an implementation of the proof by @geoffroy_2013_two_layer1.

| Model | Brief description and URL | Language |
|-------|---------------------------|----------|
| [FaIR](https://github.com/OMS-NetZero/FAIR) | Modified impulse response [@smith_2018_fairv1_3], [github.com/OMS-NetZero/FAIR](https://github.com/OMS-NetZero/FAIR) | Python |
| [GREB](https://github.com/christianstassen/greb-official) | Coarse grid energy balance [@Dommenget_2011_greb], [github.com/christianstassen/greb-official](https://github.com/christianstassen/greb-official) | Fortran 90 |
| [Hector](https://github.com/JGCRI/hector) | Upwelling-diffusion ocean energy balance [@hartin_2015_hector], [github.com/JGCRI/hector](https://github.com/JGCRI/hector) | C++  |
| [MAGICC](http://magicc.org) | Upwelling-diffusion ocean four-box (hemispheric land/ocean) energy balance [@Meinshausen_2011_magicc], [live.magicc.org](http://live.magicc.org) (Pymagicc [@Gieseke_2018_pymagicc] provides a Python wrapper in [github.com/openclimatedata/pymagicc](https://github.com/openclimatedata/pymagicc)) | Fortran 90 |
| [OSCAR](https://github.com/tgasser/OSCAR) | Energy balance with book-keeping land carbon cycle [@Gasser_2020_asdfjk], [github.com/tgasser/OSCAR](https://github.com/tgasser/OSCAR) | Python |
| [WASP](https://github.com/WASP-ESM/WASP_Earth_System_Model) | Energy balance with 8-box carbon cycle [@Goodwin_2019_ggfp6s], [github.com/WASP-ESM/WASP_Earth_System_Model](https://github.com/WASP-ESM/WASP_Earth_System_Model) | C++ |

  : Brief overview of other simple climate models available in the scientific literature. Shown is the model name, a brief description and relevant URL(s), and the programming language in which the model is written. The programming language shown is the one used for the model's core; other languages might be used in the development repositories for, e.g., plotting. For a more extensive list of simple climate models and references which describe the models in detail, see Table 1 of @rcmip_phase_1.

# Statement of need

OpenSCM Two Layer Model was designed to provide a clean, modularised, extensible interface for one of the most commonly used simple climate models.
It was used in Phase 1 of the Reduced Complexity Model Intercomparison Project [@rcmip_phase_1] as a point of comparison for the other participating models.

The FaIR model [@fair_repo] implements a mathematically equivalent model (under certain assumptions) but does not provide as clear a conversion between the two-layer model and the two-timescale response as is provided here.
We hope that this implementation could interface with other simple climate models like FaIR to allow simpler exploration of the combined behaviour of interacting climate components with minimal coupling headaches.

As implemented here, the OpenSCM Two Layer Model interface is intended to be used in research or education.

# Acknowledgements

We thank Robert Gieseke for comments on the manuscript and for all of his efforts within the OpenSCM project.

# References
<Please replace this paragraph by a description of what this PR does and alter the TODO-list below to what applies here. Please make a 'draft PR' to begin and then only mark it as ready for review once all the relevant points are done.>

- [ ] Tests added
- [ ] Documentation added
- [ ] Example added (in the documentation, to an existing notebook, or in a new notebook)
- [ ] Description in ``CHANGELOG.rst`` added (single line such as: ``(`#XX <https://github.com/openscm/openscm-twolayermodel/pull/XX>`_) Added feature which does something``)
---
name: Bug Report
about: Write a report to help us improve

---

**Describe the bug**

A clear and concise description of what the bug is, in particular:

- What did you do?
- What actually happened?
- What did you expect to happen?

**Failing Test**

Please put code (ideally in the form of a unit test) which fails below

**Screenshots**

If applicable, add screenshots to help explain your problem.

**System (please complete the following information):**

 - OS: [e.g. Windows, Linux, macOS]
 - Python and openscm commit/version [e.g. Python 3.6]

**Additional context**

Add any other context about the problem here.
---
name: Feature Request
about: Suggest an idea for this project

---

**Is your feature request related to a problem? Please describe.**

A clear and concise description of what the problem is. E.g. It's annoying that I always have to [...]

**Describe the solution you'd like**

A clear and concise description of the solution you would like to see.

**Describe alternatives you've considered**

A clear and concise description of any alternative solutions or features you've considered.

**Additional context**

Add any other context or screenshots about the feature request here.
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_, and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

The changes listed in this file are categorised as follows:

    - Added: new features
    - Changed: changes in existing functionality
    - Deprecated: soon-to-be removed features
    - Removed: now removed features
    - Fixed: any bug fixes
    - Security: in case of vulnerabilities.

v0.2.3 - 2021-04-27
-------------------

Fixed
~~~~~

- (`#34 <https://github.com/openscm/openscm-twolayermodel/pull/34>`_, `#35 <https://github.com/openscm/openscm-twolayermodel/pull/35>`_, `#36 <https://github.com/openscm/openscm-twolayermodel/pull/36>`_) Final tweaks to JOSS paper

v0.2.2 - 2021-04-27
-------------------

Added
~~~~~

- (`#33 <https://github.com/openscm/openscm-twolayermodel/pull/33>`_) Information in README and testing for conda install

Changed
~~~~~~~

- (`#32 <https://github.com/openscm/openscm-twolayermodel/pull/32>`_) Include ``LICENSE``, ``README.rst`` and ``CHANGELOG`` in package
- (`#30 <https://github.com/openscm/openscm-twolayermodel/pull/30>`_) Require ``scmdata>=0.9``
- (`#27 <https://github.com/openscm/openscm-twolayermodel/pull/27>`_) Fixed the discussion (in the relevant notebook) of how a one-layer model can be made from the two-layer implementation here

Fixed
~~~~~

- (`#30 <https://github.com/openscm/openscm-twolayermodel/pull/30>`_) Incorrect call to :meth:`scmdata.ScmRun` in tests

v0.2.1 - 2020-12-23
-------------------

Added
~~~~~

- (`#20 <https://github.com/openscm/openscm-twolayermodel/pull/20>`_) Statement of need to the README following `JOSS review <https://github.com/openjournals/joss-reviews/issues/2766>`_ (closes `#18 <https://github.com/openscm/openscm-twolayermodel/issues/18>`_)

Changed
~~~~~~~

- (`#26 <https://github.com/openscm/openscm-twolayermodel/pull/26>`_) Updated to new scmdata version (and hence new openscm-units API)
- (`#25 <https://github.com/openscm/openscm-twolayermodel/pull/25>`_) JOSS paper following `JOSS review 1 <https://github.com/openjournals/joss-reviews/issues/2766#issuecomment-718025503>`_
- (`#23 <https://github.com/openscm/openscm-twolayermodel/pull/23>`_) Moved notebooks into full documentation following `JOSS review <https://github.com/openjournals/joss-reviews/issues/2766>`_ (closes `#17 <https://github.com/openscm/openscm-twolayermodel/issues/17>`_)
- (`#21 <https://github.com/openscm/openscm-twolayermodel/pull/21>`_) Quoted pip install instructions to ensure cross-shell compatibility following `JOSS review <https://github.com/openjournals/joss-reviews/issues/2766>`_ (closes `#16 <https://github.com/openscm/openscm-twolayermodel/issues/16>`_)
- (`#20 <https://github.com/openscm/openscm-twolayermodel/pull/20>`_) Option to remove tqdm progress bar by passing ``progress=False``

v0.2.0 - 2020-10-09
-------------------

Added
~~~~~

- (`#7 <https://github.com/openscm/openscm-twolayermodel/pull/7>`_) JOSS paper draft

Changed
~~~~~~~

- (`#7 <https://github.com/openscm/openscm-twolayermodel/pull/7>`_) Require ``scmdata>=0.7``

v0.1.2 - 2020-31-07
-------------------

Changed
~~~~~~~

- (`#12 <https://github.com/openscm/openscm-twolayermodel/pull/12>`_) Upgrade to ``scmdata>=0.6.2`` so that package can be installed

v0.1.1 - 2020-06-29
-------------------

Added
~~~~~

- (`#8 <https://github.com/openscm/openscm-twolayermodel/pull/8>`_) Add notebook showing how to run a single-layer model

Changed
~~~~~~~

- (`#11 <https://github.com/openscm/openscm-twolayermodel/pull/11>`_) Re-wrote the getting started notebook
- (`#10 <https://github.com/openscm/openscm-twolayermodel/pull/10>`_) Re-wrote CHANGELOG
- (`#9 <https://github.com/openscm/openscm-twolayermodel/pull/9>`_) Update to scmdata 0.5.Y

v0.1.0 - 2020-05-15
-------------------

Added
~~~~~

- (`#3 <https://github.com/openscm/openscm-twolayermodel/pull/3>`_) Add first implementation of the models
- (`#1 <https://github.com/openscm/openscm-twolayermodel/pull/1>`_) Setup repository
OpenSCM Two Layer Model
=======================

+-------------------+----------------+--------------+
| Repository health |    |CI CD|     |  |Coverage|  |
+-------------------+----------------+--------------+

+---------------+--------+--------+
| Documentation | |Docs| | |JOSS| |
+---------------+--------+--------+

+------+------------------+----------------+------------------+
| PyPI |  |PyPI Install|  |     |PyPI|     |  |PyPI Version|  |
+------+------------------+----------------+------------------+

+-------+-----------------+-------------------+-----------------+
| Conda | |conda install| | |conda platforms| | |conda version| |
+-------+-----------------+-------------------+-----------------+

+-----------------+----------------+---------------+-----------+
|   Other info    | |Contributors| | |Last Commit| | |License| |
+-----------------+----------------+---------------+-----------+


Brief summary
+++++++++++++

.. sec-begin-long-description
.. sec-begin-index

OpenSCM two layer model contains implementations of the two layer radiative forcing driven models by `Held et al. <https://journals.ametsoc.org/doi/full/10.1175/2009JCLI3466.1>`_, `Geoffroy et al. <https://journals.ametsoc.org/doi/pdf/10.1175/JCLI-D-12-00195.1>`_ and `Bloch-Johnson et al. <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2015GL064240>`_

OpenSCM Two Layer Model was designed to provide a clean, modularised, extensible interface for one of the most commonly used simple climate models.
It was used in `Phase 1 of the Reduced Complexity Model Intercomparison Project <https://doi.org/10.5194/gmd-13-5175-2020>`_ as a point of comparison for the other participating models.

The `FaIR model <https://github.com/OMS-NetZero/FAIR>`_ implements a mathematically equivalent model (under certain assumptions) but does not provide as clear a conversion between the two-layer model and the two-timescale response as is provided here.
We hope that this implementation could interface with other simple climate models like FaIR to allow simpler exploration of the combined behaviour of interacting climate components with minimal coupling headaches.

As implemented here, the "OpenSCM Two Layer Model" interface is target at researchers and as an education tool.
Users from other fields are of course encouraged to use it if they find it helpful.

.. sec-end-index

License
-------

.. sec-begin-license

OpenSCM two layer model is free software under a BSD 3-Clause License, see
`LICENSE <https://github.com/openscm/openscm-twolayermodel/blob/master/LICENSE>`_.

.. sec-end-license
.. sec-end-long-description

.. sec-begin-installation

Installation
------------

OpenSCM two layer model has only two dependencies:

.. begin-dependencies

- scmdata>=0.9
- tqdm

.. end-dependencies

OpenSCM two layer model can be installed with pip

.. code:: bash

    pip install openscm-twolayermodel

If you also want to run the example notebooks install additional
dependencies using

.. code:: bash

    pip install "openscm-twolayermodel[notebooks]"

**Coming soon** OpenSCM two layer model can also be installed with conda

.. code:: bash

    conda install -c conda-forge openscm-twolayermodel

We do not ship our tests with the packages.
If you wish to run the tests, you must install from source (or download the tests separately and run them on your installation).

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

To install from source, simply clone the repository and then install it using pip e.g. ``pip install ".[dev]"``.
Having done this, the tests can be run with ``pytest tests`` or using the ``Makefile`` (``make test`` will run only the code tests, ``make checks`` will run the code tests and all other tests e.g. linting, notebooks, documentation).

.. sec-end-installation

For more details, see the `development section of the docs <https://openscm-two-layer-model.readthedocs.io/en/latest/development.html>`_.

Documentation
-------------

Documentation can be found at our `documentation pages <https://openscm-two-layer-model.readthedocs.io/en/latest/>`_
(we are thankful to `Read the Docs <https://readthedocs.org/>`_ for hosting us).

Getting help
------------

.. sec-begin-getting-help

If you have any issues or would like to discuss a feature request, please raise them in the `OpenSCM Two Layer Model issue tracker <https://github.com/openscm/openscm-twolayermodel/issues>`_.
If your issue is a feature request or a bug, please use the templates available, otherwise, simply open a normal issue.

.. sec-end-getting-help

Contributing
------------

Please see the `Development section of the docs <https://openscm-two-layer-model.readthedocs.io/en/latest/development.html>`_.

.. sec-begin-links

.. |CI CD| image:: https://github.com/openscm/openscm-twolayermodel/workflows/OpenSCM%20Two%20Layer%20Model%20CI-CD/badge.svg
    :target: https://github.com/openscm/openscm-twolayermodel/actions?query=workflow%3A%22OpenSCM+Two+Layer+Model+CI-CD%22
.. |Coverage| image:: https://codecov.io/gh/openscm/openscm-twolayermodel/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/openscm/openscm-twolayermodel
.. |Docs| image:: https://readthedocs.org/projects/openscm-two-layer-model/badge/?version=latest
    :target: https://openscm-two-layer-model.readthedocs.io/en/latest/?badge=latest
.. |JOSS| image:: https://joss.theoj.org/papers/94a3759c9ea117499b90c56421ef4857/status.svg
    :target: https://joss.theoj.org/papers/94a3759c9ea117499b90c56421ef4857
.. |PyPI Install| image:: https://github.com/openscm/openscm-twolayermodel/workflows/Test%20PyPI%20install/badge.svg
    :target: https://github.com/openscm/openscm-twolayermodel/actions?query=workflow%3A%22Test+PyPI+install%22
.. |PyPI| image:: https://img.shields.io/pypi/pyversions/openscm-twolayermodel.svg
    :target: https://pypi.org/project/openscm-twolayermodel/
.. |PyPI Version| image:: https://img.shields.io/pypi/v/openscm-twolayermodel.svg
    :target: https://pypi.org/project/openscm-twolayermodel/
.. |conda install| image:: https://github.com/openscm/openscm-twolayermodel/workflows/Test%20conda%20install/badge.svg
    :target: https://github.com/openscm/openscm-twolayermodel/actions?query=workflow%3A%22Test+conda+install%22
.. |conda platforms| image:: https://img.shields.io/conda/pn/conda-forge/openscm-twolayermodel.svg
    :target: https://anaconda.org/conda-forge/openscm-twolayermodel
.. |conda version| image:: https://img.shields.io/conda/vn/conda-forge/openscm-twolayermodel.svg
    :target: https://anaconda.org/conda-forge/openscm-twolayermodel
.. |Contributors| image:: https://img.shields.io/github/contributors/openscm/openscm-twolayermodel.svg
    :target: https://github.com/openscm/openscm-twolayermodel/graphs/contributors
.. |Last Commit| image:: https://img.shields.io/github/last-commit/openscm/openscm-twolayermodel.svg
    :target: https://github.com/openscm/openscm-twolayermodel/commits/master
.. |License| image:: https://img.shields.io/github/license/openscm/openscm-twolayermodel.svg
    :target: https://github.com/openscm/openscm-twolayermodel/blob/master/LICENSE

.. sec-end-links
.. _development:

Development
===========

If you're interested in contributing to OpenSCM Two Layer Model, we'd love to have you on board!
This section of the docs details how to get setup to contribute and how best to communicate.

.. contents:: :local:

Contributing
------------

All contributions are welcome, some possible suggestions include:

- tutorials (or support questions which, once solved, result in a new tutorial :D)
- blog posts
- improving the documentation
- bug reports
- feature requests
- pull requests

Please report issues or discuss feature requests in the `OpenSCM Two Layer Model issue tracker`_.
If your issue is a feature request or a bug, please use the templates available, otherwise, simply open a normal issue.

As a contributor, please follow a couple of conventions:

- Create issues in the `OpenSCM Two Layer Model issue tracker`_ for changes and enhancements, this ensures that everyone in the community has a chance to comment
- Be welcoming to newcomers and encourage diverse new contributors from all backgrounds: see the `Python Community Code of Conduct <https://www.python.org/psf/codeofconduct/>`_
- Only push to your own branches, this allows people to force push to their own branches as they need without fear or causing others headaches
- Start all pull requests as draft pull requests and only mark them as ready for review once they've been rebased onto master, this makes it much simpler for reviewers
- Try and make lots of small pull requests, this makes it easier for reviewers and faster for everyone as review time grows exponentially with the number of lines in a pull request


Getting setup
-------------

To get setup as a developer, we recommend the following steps (if any of these tools are unfamiliar, please see the resources we recommend in `Development tools`_):

#. Install conda and make
#. Run ``make virtual-environment``, if that fails you can try doing it manually

    #. Change your current directory to OpenSCM Two Layer Model's root directory (i.e. the one which contains ``README.rst``), ``cd openscm-twolayermodel``
    #. Create a virtual environment to use with OpenSCM Two Layer Model ``python3 -m venv venv``
    #. Activate your virtual environment ``source ./venv/bin/activate``
    #. Upgrade pip ``pip install --upgrade pip``
    #. Install the development dependencies (very important, make sure your virtual environment is active before doing this) ``pip install -e .[dev]``

#. Make sure the tests pass by running ``make checks``, if that fails the commands can be read out of the ``Makefile``


Getting help
~~~~~~~~~~~~

Whilst developing, unexpected things can go wrong (that's why it's called 'developing', if we knew what we were doing, it would already be 'developed').
Normally, the fastest way to solve an issue is to contact us via the `issue tracker <https://github.com/openscm/openscm-twolayermodel/issues>`_.
The other option is to debug yourself.
For this purpose, we provide a list of the tools we use during our development as starting points for your search to find what has gone wrong.

Development tools
+++++++++++++++++

This list of development tools is what we rely on to develop OpenSCM Two Layer Model reliably and reproducibly.
It gives you a few starting points in case things do go inexplicably wrong and you want to work out why.
We include links with each of these tools to starting points that we think are useful, in case you want to learn more.

- `Git <http://swcarpentry.github.io/git-novice/>`_

- `Make <https://swcarpentry.github.io/make-novice/>`_

- `Conda virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_

- `Pip and pip virtual environments <https://www.dabapps.com/blog/introduction-to-pip-and-virtualenv-python/>`_

- `Tests <https://semaphoreci.com/community/tutorials/testing-python-applications-with-pytest>`_

    - we use a blend of `pytest <https://docs.pytest.org/en/latest/>`_ and the inbuilt Python testing capabilities for our tests so checkout what we've already done in ``tests`` to get a feel for how it works

- `Continuous integration (CI) <https://help.github.com/en/actions>`_ (also `brief intro blog post <https://gabrieltanner.org/blog/an-introduction-to-github-actions>`_ and a `longer read <https://dev.to/bnb/an-unintentionally-comprehensive-introduction-to-github-actions-ci-blm>`_)

    - we use GitHub CI for our CI but there are a number of good providers

- `Jupyter Notebooks <https://medium.com/codingthesmartway-com-blog/getting-started-with-jupyter-notebook-for-python-4e7082bd5d46>`_

    - Jupyter is automatically included in your virtual environment if you follow our `Getting setup`_ instructions

- Sphinx_


Other tools
+++++++++++

We also use some other tools which aren't necessarily the most familiar.
Here we provide a list of these along with useful resources.

.. _regular-expressions:

- `Regular expressions <https://www.oreilly.com/ideas/an-introduction-to-regular-expressions>`_

    - we use `regex101.com <regex101.com>`_ to help us write and check our regular expressions, make sure the language is set to Python to make your life easy!

Formatting
----------

To help us focus on what the code does, not how it looks, we use a couple of automatic formatting tools.
These automatically format the code for us and tell use where the errors are.
To use them, after setting yourself up (see `Getting setup`_), simply run ``make format``.
Note that ``make format`` can only be run if you have committed all your work i.e. your working directory is 'clean'.
This restriction is made to ensure that you don't format code without being able to undo it, just in case something goes wrong.


Buiding the docs
----------------

After setting yourself up (see `Getting setup`_), building the docs is as simple as running ``make docs`` (note, run ``make -B docs`` to force the docs to rebuild and ignore make when it says '... index.html is up to date').
This will build the docs for you.
You can preview them by opening ``docs/build/html/index.html`` in a browser.

For documentation we use Sphinx_.
To get ourselves started with Sphinx, we started with `this example <https://pythonhosted.org/an_example_pypi_project/sphinx.html>`_ then used `Sphinx's getting started guide <http://www.sphinx-doc.org/en/master/usage/quickstart.html>`_.


Gotchas
~~~~~~~

To get Sphinx to generate pdfs (rarely worth the hassle), you require `Latexmk <https://mg.readthedocs.io/latexmk.html>`_.
On a Mac this can be installed with ``sudo tlmgr install latexmk``.
You will most likely also need to install some other packages (if you don't have the full distribution).
You can check which package contains any missing files with ``tlmgr search --global --file [filename]``.
You can then install the packages with ``sudo tlmgr install [package]``.


Docstring style
~~~~~~~~~~~~~~~

For our docstrings we use numpy style docstrings.
For more information on these, `here is the full guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_ and `the quick reference we also use <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.


Releasing
---------

First step
~~~~~~~~~~

#. Test installation with dependencies ``make test-install``
#. Update ``CHANGELOG.rst``

    - add a header for the new version between ``master`` and the latest bullet point
    - this should leave the section underneath the master header empty

#. ``git add .``
#. ``git commit -m "Prepare for release of vX.Y.Z"``
#. ``git tag vX.Y.Z``
#. Test version updated as intended with ``make test-install``


Push to repository
~~~~~~~~~~~~~~~~~~

To do the release, push the tags and the repository state.

#. ``git push``
#. ``git push --tags``

Assuming all the checks pass, this automatically triggers a release on PyPI via the ``.github/workflows/ci-cd-workflow.yml`` action.


Why is there a ``Makefile`` in a pure Python repository?
--------------------------------------------------------

Whilst it may not be standard practice, a ``Makefile`` is a simple way to automate general setup (environment setup in particular).
Hence we have one here which basically acts as a notes file for how to do all those little jobs which we often forget e.g. setting up environments, running tests (and making sure we're in the right environment), building docs, setting up auxillary bits and pieces.

.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _OpenSCM Two Layer Model issue tracker: https://github.com/openscm/openscm-twolayermodel/issues
.. _`OpenSCM Two Layer Model's PyPI`: https://pypi.org/project/openscm-twolayermodel/
Usage
=====

Here we provide examples of OpenSCM two layer model's behaviour and usage.
The source code of these usage examples is available in the folder
`docs/source/usage`_ of the `GitHub repository`_.

.. _`docs/source/usage`:
   https://github.com/openscm/openscm-twolayermodel/tree/master/docs/source/usage

.. _`GitHub repository`:
   https://github.com/openscm/openscm-twolayermodel

Basic demos
+++++++++++

.. toctree::
   :maxdepth: 1

   usage/getting-started.ipynb
   usage/running-scenarios.ipynb

More detail
+++++++++++

.. toctree::
   :maxdepth: 1

   usage/impulse-response-equivalence.ipynb
   usage/one-layer-model.ipynb
.. _utils-reference:

Utils API
---------

.. automodule:: openscm_twolayermodel.utils
.. _impulse_response_model-reference:

Impulse Response Model API
--------------------------

.. automodule:: openscm_twolayermodel.impulse_response_model
.. _two_layer_model-reference:

Two Layer Model API
-------------------

.. automodule:: openscm_twolayermodel.two_layer_model
.. include:: ../../README.rst
    :start-after: sec-begin-installation
    :end-before: sec-end-installation

For more details, see the :ref:`development section of the docs <development>`.
.. _base-reference:

Base API
---------

.. automodule:: openscm_twolayermodel.base
.. OpenSCM Two Layer Model documentation master file, created by
   sphinx-quickstart on Thu Mar 19 12:35:34 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

OpenSCM Two Layer Model
=======================

.. include:: ../../README.rst
    :start-after: sec-begin-index
    :end-before: sec-end-index

.. include:: ../../README.rst
    :start-after: sec-begin-license
    :end-before: sec-end-license

.. include:: ../../README.rst
    :start-after: sec-begin-getting-help
    :end-before: sec-end-getting-help

.. toctree::
    :maxdepth: 2
    :caption: Documentation

    installation
    usage
    development

.. toctree::
    :maxdepth: 2
    :caption: API reference

    base
    impulse_response_model
    two_layer_model
    constants
    errors
    utils

.. toctree::
    :maxdepth: 2
    :caption: Versions

    changelog

Index
-----

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`
.. _constants-reference:

Constants API
-------------

.. automodule:: openscm_twolayermodel.constants
.. include:: ../../CHANGELOG.rst
.. _errors-reference:

Errors API
----------

.. automodule:: openscm_twolayermodel.errors
