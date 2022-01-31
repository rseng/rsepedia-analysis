QuTiP: Quantum Toolbox in Python
================================

[A. Pitchford](https://github.com/ajgpitch),
[C. Granade](https://github.com/cgranade),
[A. Grimsmo](https://github.com/arnelg),
[N. Shammah](https://github.com/nathanshammah),
[S. Ahmed](https://github.com/quantshah),
[N. Lambert](https://github.com/nwlambert),
[E. Giguère](https://github.com/ericgig),
[B. Li](https://github.com/boxili),
[J. Lishman](https://github.com/jakelishman),
[S. Cross](https://github.com/hodgestar),
[P. D. Nation](https://github.com/nonhermitian),
and [J. R. Johansson](https://github.com/jrjohansson)

[![Build Status](https://github.com/qutip/qutip/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/qutip/qutip/actions/workflows/tests.yml)
[![Coverage Status](https://img.shields.io/coveralls/qutip/qutip.svg?logo=Coveralls)](https://coveralls.io/r/qutip/qutip)
[![Maintainability](https://api.codeclimate.com/v1/badges/df502674f1dfa1f1b67a/maintainability)](https://codeclimate.com/github/qutip/qutip/maintainability)
[![license](https://img.shields.io/badge/license-New%20BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)  
[![PyPi Downloads](https://img.shields.io/pypi/dm/qutip?label=downloads%20%7C%20pip&logo=PyPI)](https://pypi.org/project/qutip)
[![Conda-Forge Downloads](https://img.shields.io/conda/dn/conda-forge/qutip?label=downloads%20%7C%20conda&logo=Conda-Forge)](https://anaconda.org/conda-forge/qutip)

QuTiP is open-source software for simulating the dynamics of closed and open quantum systems.
The QuTiP library uses the excellent Numpy, Scipy, and Cython packages as numerical backend, and graphical output is provided by Matplotlib.
QuTiP aims to provide user-friendly and efficient numerical simulations of a wide variety of quantum mechanical problems, including those with Hamiltonians and/or collapse operators with arbitrary time-dependence, commonly found in a wide range of physics applications.
QuTiP is freely available for use and/or modification, and it can be used on all Unix-based platforms and on Windows.
Being free of any licensing fees, QuTiP is ideal for exploring quantum mechanics in research as well as in the classroom.

Support
-------

[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=flat)](https://unitary.fund)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)

We are proud to be affiliated with [Unitary Fund](https://unitary.fund) and [numFOCUS](https://numfocus.org).
QuTiP development is supported by [Nori's lab](https://dml.riken.jp/) at RIKEN, by the University of Sherbrooke, and by Aberystwyth University, [among other supporting organizations](https://qutip.org/#supporting-organizations).


Installation
------------

[![Pip Package](https://img.shields.io/pypi/v/qutip?logo=PyPI)](https://pypi.org/project/qutip)
[![Conda-Forge Package](https://img.shields.io/conda/vn/conda-forge/qutip?logo=Conda-Forge)](https://anaconda.org/conda-forge/qutip)

QuTiP is available on both `pip` and `conda` (the latter in the `conda-forge` channel).
You can install QuTiP from `pip` by doing

```bash
pip install qutip
```

to get the minimal installation.
You can instead use the target `qutip[full]` to install QuTiP with all its optional dependencies.
For more details, including instructions on how to build from source, see [the detailed installation guide in the documentation](https://qutip.org/docs/latest/installation.html).

All back releases are also available for download in the [releases section of this repository](https://github.com/qutip/qutip/releases), where you can also find per-version changelogs.
For the most complete set of release notes and changelogs for historic versions, see the [changelog](https://qutip.org/docs/latest/changelog.html) section in the documentation.


Documentation
-------------

The documentation for official releases, in HTML and PDF formats, can be found in the [documentation section of the QuTiP website](https://qutip.org/documentation.html).
The latest development documentation is available in this repository in the `doc` folder.

A [selection of demonstration notebooks is available](https://qutip.org/tutorials.html), which demonstrate some of the many features of QuTiP.
These are stored in the [qutip/qutip-notebooks repository](https://github.com/qutip/qutip-notebooks) here on GitHub.
You can run the notebooks online using myBinder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qutip/qutip-notebooks/master?filepath=index.ipynb)

Contribute
----------

You are most welcome to contribute to QuTiP development by forking this repository and sending pull requests, or filing bug reports at the [issues page](https://github.com/qutip/qutip/issues).
You can also help out with users' questions, or discuss proposed changes in the [QuTiP discussion group](https://groups.google.com/g/qutip).
All code contributions are acknowledged in the [contributors](https://qutip.org/docs/latest/contributors.html) section in the documentation.

For more information, including technical advice, please see the ["contributing to QuTiP development" section of the documentation](https://qutip.org/docs/latest/development/contributing.html).


Citing QuTiP
------------

If you use QuTiP in your research, please cite the original QuTiP papers that are available [here](https://dml.riken.jp/?s=QuTiP).
# Contributor Covenant Code of Conduct

As contributors and maintainers of this project, and in the interest of fostering an open and welcoming community, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic addresses, without explicit permission
* Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. By adopting this Code of Conduct, project maintainers commit themselves to fairly and consistently applying these principles to every aspect of managing this project. Project maintainers who do not follow or enforce the Code of Conduct may be permanently removed from the project team.

This code of conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. 

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers. 

This Code of Conduct is adapted from the Contributor Covenant , version 1.2.0, available at https://www.contributor-covenant.org/version/1/2/0/code-of-conduct.html 

[homepage]: https://contributor-covenant.org
[version]: https://contributor-covenant.org/version/1/2/
# Contributing to QuTiP Development

You are most welcome to contribute to QuTiP development by forking this repository and sending pull requests, or filing bug reports at the [issues page](https://github.com/qutip/qutip/issues).
You can also help out with users' questions, or discuss proposed changes in the [QuTiP discussion group](https://groups.google.com/g/qutip).
All code contributions are acknowledged in the [contributors](https://qutip.org/docs/latest/contributors.html) section in the documentation.

For more information, including technical advice, please see the ["contributing to QuTiP development" section of the documentation](https://qutip.org/docs/latest/development/contributing.html).
**Checklist**
Thank you for contributing to QuTiP! Please make sure you have finished the following tasks before opening the PR.

- [ ] Please read [Contributing to QuTiP Development](http://qutip.org/docs/latest/development/contributing.html)
- [ ] Contributions to qutip should follow the [pep8 style](https://www.python.org/dev/peps/pep-0008/).
You can use [pycodestyle](http://pycodestyle.pycqa.org/en/latest/index.html) to check your code automatically
- [ ] Please add tests to cover your changes if applicable.
- [ ] If the behavior of the code has changed or new feature has been added, please also update the documentation in the `doc` folder, and the [notebook](https://github.com/qutip/qutip-notebooks). Feel free to ask if you are not sure.

Delete this checklist after you have completed all the tasks. If you have not finished them all, you can also open a [Draft Pull Request](https://github.blog/2019-02-14-introducing-draft-pull-requests/) to let the others know this on-going work and keep this checklist in the PR description.

**Description**
Describe here the proposed change.

**Related issues or PRs**
Please mention the related issues or PRs here. If the PR fixes an issue, use the keyword fix/fixes/fixed followed by the issue id, e.g. fix #1184

**Changelog**
Give a short description of the PR in a few words. This will be shown in the QuTiP change log after the PR gets merged.
For example: 
Fixed error checking for null matrix in essolve.
Added option for specifying resolution in Bloch.save function.
Repository for QuTiP documentation
==================================

This repository contains the source files for the QuTiP documentation.

For pre-built documentation, see https://www.qutip.org/documentation.html

Building
--------

The main Python requirements for the documentation are `sphinx`, `sphinx-gallery`, `sphinx_rtd_theme`, `numpydoc` and `ipython`.
You should build or install the version of QuTiP you want to build the documentation against in the same environment.
You will also need a sensible copy of `make`, and if you want to build the LaTeX documentation then also a `pdflatex` distribution.
As of 2021-04-20, the `conda` recipe for `sphinx_rtd_theme` is rather old compared to the `pip` version, so it's recommended to use a mostly `pip`-managed environment to do the documentation build.

The simplest way to get a functional build environment is to use the `requirements.txt` file in this repository, which completely defines a known-good `pip` environment (tested on Python 3.8, but not necessarily limited to it).
If you typically use conda, the way to do this is
```bash
$ conda create -n qutip-doc-build python=3.8
$ conda activate qutip-doc-build
$ pip install -r /path/to/qutip/doc/requirements.txt
```
You will also need to build or install the main QuTiP library in the same environment.
If you simply want to build the documentation without editing the main library, you can install a release version of QuTiP with `pip install qutip`.
Otherwise, refer to [the main repository](https://github.com/qutip/qutip) for the current process to build from source.
You need to have the optional QuTiP dependency `Cython` to build the documentation, but this is included in this repository's `requirements.txt` so you do not need to do anything separately.

After you have done this, you can effect the build with `make`.
The targets you might want are `html`, `latexpdf` and `clean`, which build the HTML pages, build the PDFs, and delete all built files respectively.
For example, to build the HTML files only, use
```bash
$ make html
```

*Note (2021-04-20):* the documentation build is currently broken on Windows due to incompatibilities in the main library in multiprocessing components.

Writing User Guides
-------------------

The user guide provides an overview of QuTiP's functionality. The guide is composed of individual reStructuredText (`.rst`) files which each get rendered as a webpage. Each page typically tackles one area of functionality. To learn more about how to write `.rst` files, it is useful to follow the [Sphinx Guide](https://www.sphinx-doc.org/en/master/usage/index.html).

The documentation build also utilizes a number of [Sphinx Extensions](https://www.sphinx-doc.org/en/master/usage/extensions/index.html) including but not limited to
[doctest](https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html), [autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html), [sphinx gallery](https://sphinx-gallery.github.io/stable/index.html), [plot](https://matthew-brett.github.io/nb2plots/nbplots.html#module-nb2plots.nbplots). Additional extensions can be configured in the `conf.py` file.

Tests can also be run on examples in the documentation using the doctest extension
and plots are generated using the `plot` directive. For more specific
guidelines on how to incorporate code examples into the guide, refer to (insert reference).
.. _developers:

************
Developers
************


.. _developers-lead:

Lead Developers
===============

- `Alex Pitchford <https://github.com/ajgpitch>`_
- `Nathan Shammah <https://nathanshammah.com/>`_
- `Shahnawaz Ahmed <http://sahmed.in/>`_
- `Neill Lambert <https://github.com/nwlambert>`_
- `Eric Giguère <https://github.com/Ericgig>`_
- `Boxi Li <https://github.com/BoxiLi>`_
- `Jake Lishman <https://binhbar.com>`_
- `Simon Cross <http://hodgestar.za.net/>`_

Past Lead Developers
====================

- `Robert Johansson <https://jrjohansson.github.io/research.html>`_ (RIKEN)
- `Paul Nation <http://nqdl.korea.ac.kr>`_ (Korea University)
- `Chris Granade <https://www.cgranade.com>`_
- `Arne Grimsmo <https://www.sydney.edu.au/science/about/our-people/academic-staff/arne-grimsmo.html>`_


.. _developers-contributors:

Contributors
============

.. note::

	Anyone is welcome to contribute to QuTiP.
        If you are interested in helping, please let us know!

- Abhisek Upadhyaya
- Adriaan
- Alexander Pitchford
- Alexios-xi
- Amit
- Anubhav Vardhan
- Arie van Deursen
- Arne Grimsmo
- Arne Hamann
- Asier Galicia Martinez
- Ben Bartlett
- Ben Criger
- Ben Jones
- Bo Yang
- Boxi Li
- Canoming
- Christoph Gohlke
- Christopher Granade
- Craig Gidney
- Denis Vasilyev
- Dominic Meiser
- Drew Parsons
- Eric Giguère
- Eric Hontz
- Felipe Bivort Haiek
- Florestan Ziem
- Gilbert Shih
- Harry Adams
- Ivan Carvalho
- Jake Lishman
- Jevon Longdell
- Johannes Feist
- Jonas Hoersch
- Jonas Neergaard-Nielsen
- Jonathan A. Gross
- Julian Iacoponi
- Kevin Fischer
- Laurence Stant
- Louis Tessler
- Lucas Verney
- Marco David
- Marek Narozniak
- Markus Baden
- Martín Sande
- Mateo Laguna
- Matthew O'Brien
- Michael Goerz
- Michael V. DePalatis
- Moritz Oberhauser
- Nathan Shammah
- Neill Lambert
- Nicolas Quesada
- Nikolas Tezak
- Nithin Ramu
- Paul Nation
- Peter Kirton
- Philipp Schindler
- Piotr Migdal
- Rajiv-B
- Ray Ganardi
- Reinier Heeres
- Richard Brierley
- Robert Johansson
- Sam Griffiths
- Samesh Lakhotia
- Sebastian Krämer
- Shahnawaz Ahmed
- Sidhant Saraogi
- Simon Cross
- Simon Humpohl
- Simon Whalen
- Stefan Krastanov
- Tarun Raheja
- Thomas Walker
- Viacheslav Ostroukh
- Vlad Negnevitsky
- Wojciech Rzadkowski
- Xiaodong Qi
- Xiaoliang Wu
- Yariv Yanay
- YouWei Zhao
- alex
- eliegenois
- essence-of-waqf
- fhenneke
- gecrooks
- jakobjakobson13
- maij
- sbisw002
- yuri@FreeBSD
- Élie Gouzien
.. _copyright:

***********************
Copyright and Licensing
***********************

The text of this documentation is licensed under the `Creative Commons Attribution 3.0 Unported License <https://creativecommons.org/licenses/by/3.0/>`_.
Unless specifically indicated otherwise, all code samples, the source code of QuTiP, and its reproductions in this documentation, are licensed under the terms of the 3-clause BSD license, reproduced below.

License Terms for Documentation Text
====================================

The canonical form of this license is available at `https://creativecommons.org/licenses/by/3.0/ <https://creativecommons.org/licenses/by/3.0/>`_, which should be considered the binding version of this license.
It is reproduced here for convenience.

.. include:: LICENSE_cc-by-3.0.txt


License Terms for Source Code of QuTiP and Code Samples
=======================================================

.. include:: ../LICENSE.txt
.. _frontmatter:

*************
Frontmatter
*************

.. _about-docs:

About This Documentation
==========================

This document contains a user guide and automatically generated API documentation for QuTiP. A PDF version of this text is available at the `documentation page <https://www.qutip.org/documentation.html>`_.

**For more information see the** `QuTiP project web page`_.

.. _QuTiP project web page: https://www.qutip.org


:Author: J.R. Johansson

:Author: P.D. Nation

:Author: Alexander Pitchford

:Author: Arne Grimsmo

:Author: Chris Grenade

:Author: Nathan Shammah

:Author: Shahnawaz Ahmed

:Author: Neill Lambert

:Author: Eric Giguere

:Author: Boxi Li

:Author: Jake Lishman

:Author: Simon Cross

:release: |release|

:copyright:
   The text of this documentation is licensed under the Creative Commons Attribution 3.0 Unported License.
   All contained code samples, and the source code of QuTiP, are licensed under the 3-clause BSD licence.
   Full details of the copyright notices can be found on the `Copyright and Licensing <copyright>`_ page of this documentation.

.. _citing-qutip:

Citing This Project
==========================

If you find this project useful, then please cite:

.. centered:: J. R. Johansson, P.D. Nation, and F. Nori, "QuTiP 2: A Python framework for the dynamics of open quantum systems", Comp. Phys. Comm. **184**, 1234 (2013).

or

.. centered:: J. R. Johansson, P.D. Nation, and F. Nori, "QuTiP: An open-source Python framework for the dynamics of open quantum systems", Comp. Phys. Comm. **183**, 1760 (2012).


which may also be downloaded from https://arxiv.org/abs/1211.6518 or https://arxiv.org/abs/1110.0573, respectively.

.. _funding-qutip:

Funding
=======
QuTiP is developed under the auspice of the non-profit organizations:

.. _image-numfocus:

.. figure:: figures/NumFocus_logo.png
   :width: 3in
   :figclass: align-center

.. _image-unitaryfund:

.. figure:: figures/unitaryfund_logo.png
   :width: 3in
   :figclass: align-center

QuTiP was partially supported by

.. _image-jsps:

.. figure:: figures/jsps.jpg
   :width: 2in
   :figclass: align-center

.. _image-riken:

.. figure:: figures/riken-logo.png
	:width: 1.5in
	:figclass: align-center

.. _image-korea:

.. figure:: figures/korea-logo.png
	:width: 2in
	:figclass: align-center

.. figure:: figures/inst_quant_sher.png
	:width: 2in
	:figclass: align-center

.. _about:

About QuTiP
===========

Every quantum system encountered in the real world is an open quantum system. For although much care is taken experimentally to eliminate the unwanted influence of external interactions, there remains, if ever so slight, a coupling between the system of interest and the external world. In addition, any measurement performed on the system necessarily involves coupling to the measuring device, therefore introducing an additional source of external influence. Consequently, developing the necessary tools, both theoretical and numerical, to account for the interactions between a system and its environment is an essential step in understanding the dynamics of practical quantum systems.

In general, for all but the most basic of Hamiltonians, an analytical description of the system dynamics is not possible, and one must resort to numerical simulations of the equations of motion. In absence of a quantum computer, these simulations must be carried out using classical computing techniques, where the exponentially increasing dimensionality of the underlying Hilbert space severely limits the size of system that can be efficiently simulated. However, in many fields such as quantum optics, trapped ions, superconducting circuit devices, and most recently nanomechanical systems, it is possible to design systems using a small number of effective oscillator and spin components, excited by a limited number of quanta, that are amenable to classical simulation in a truncated Hilbert space.

The Quantum Toolbox in Python, or QuTiP, is an open-source framework written in the Python programming language, designed for simulating the open quantum dynamics of systems such as those listed above. This framework distinguishes itself from other available software solutions in providing the following advantages:

* QuTiP relies entirely on open-source software.  You are free to modify and use it as you wish with no licensing fees or limitations.

* QuTiP is based on the Python scripting language, providing easy to read, fast code generation without the need to compile after modification.

* The numerics underlying QuTiP are time-tested algorithms that run at C-code speeds, thanks to the `Numpy <https://numpy.org>`_, `Scipy <https://scipy.org>`_, and `Cython <https://cython.org>`_ libraries, and are based on many of the same algorithms used in propriety software.

* QuTiP allows for solving the dynamics of Hamiltonians with (almost) arbitrary time-dependence, including collapse operators.

* Time-dependent problems can be automatically compiled into C++-code at run-time for increased performance.

* Takes advantage of the multiple processing cores found in essentially all modern computers.

* QuTiP was designed from the start to require a minimal learning curve for those users who have experience using the popular quantum optics toolbox by Sze M. Tan.

* Includes the ability to create high-quality plots, and animations, using the excellent `Matplotlib <https://matplotlib.org>`_ package.


For detailed information about new features of each release of QuTiP, see the :ref:`changelog`.

.. _plugin-qutip:

QuTiP Plugins
=============

Several libraries depend on QuTiP heavily making QuTiP a super-library

:Matsubara: `Matsubara <https://matsubara.readthedocs.io/en/latest/>`_ is a plugin to study the ultrastrong coupling regime with structured baths

:QNET: `QNET <https://qnet.readthedocs.io/en/latest/readme.html>`_ is a computer algebra package for quantum mechanics and photonic quantum networks

.. _libraries:

Libraries Using QuTiP
=====================

Several libraries rely on QuTiP for quantum physics or quantum information processing. Some of them are:

:Krotov: `Krotov <https://qucontrol.github.io/krotov/v1.2.0/01_overview.html>`_ focuses on the python implementation of Krotov's method for quantum optimal control

:pyEPR: `pyEPR <https://pyepr-docs.readthedocs.io/en/latest/index.html>`_ interfaces classical distributed microwave analysis with that of quantum structures and hamiltonians by providing easy to use analysis function and automation for the design of quantum chips

:scQubits: `scQubits <https://scqubits.readthedocs.io/en/latest/>`_ is a Python library which provides a convenient way to simulate superconducting qubits by providing an interface to QuTiP

:SimulaQron: `SimulaQron <https://softwarequtech.github.io/SimulaQron/html/index.html>`_ is a distributed simulation of the end nodes in a quantum internet with the specific goal to explore application development

:QInfer: `QInfer <http://qinfer.org/>`_ is a library for working with sequential Monte Carlo methods for parameter estimation in quantum information

:QPtomographer: `QPtomographer <https://qptomographer.readthedocs.io/en/latest/>`_ derive quantum error bars for quantum processes in terms of the diamond norm to a reference quantum channel

:QuNetSim: `QuNetSim <https://tqsd.github.io/QuNetSim/_build/intro.html>`_ is a quantum networking simulation framework to develop and test protocols for quantum networks

:qupulse: `qupulse <https://qupulse.readthedocs.io/en/latest/>`_ is a toolkit to facilitate experiments involving pulse driven state manipulation of physical qubits

:Pulser: `Pulser <https://pulser.readthedocs.io/en/latest/>`_ is a framework for composing, simulating and executing pulse sequences for neutral-atom quantum devices.



Contributing to QuTiP
=====================

We welcome anyone who is interested in helping us make QuTiP the best package for simulating quantum systems.
There are :ref:`detailed instructions on how to contribute code and documentation <development-contributing>` in the developers' section of this guide.
You can also help out our users by answering questions in the `QuTiP discussion mailing list <https://groups.google.com/g/qutip>`_, or by raising issues in `the main GitHub repository <https://github.com/qutip/qutip>`_ if you find any bugs.
Anyone who contributes code will be duly recognized.
Even small contributions are noted.
See :ref:`developers-contributors` for a list of people who have helped in one way or another.
.. This file can be edited using retext 6.1 https://github.com/retext-project/retext

.. _install:

**************
Installation
**************

.. _quick-start:

Quick Start
===========

From QuTiP version 4.6 onwards, you should be able to get a working version of QuTiP with the standard

.. code-block:: bash

   pip install qutip

It is not recommended to install any packages directly into the system Python environment; consider using ``pip`` or ``conda`` virtual environments to keep your operating system space clean, and to have more control over Python and other package versions.

You do not need to worry about the details on the rest of this page unless this command did not work, but do also read the next section for the list of optional dependencies.
The rest of this page covers `installation directly from conda <install-with-conda_>`_, `installation from source <install-from-source_>`_, and `additional considerations when working on Windows <install-on-windows_>`_.


.. _install-requires:

General Requirements
=====================

QuTiP depends on several open-source libraries for scientific computing in the Python programming language.
The following packages are currently required:

.. cssclass:: table-striped

+----------------+--------------+-----------------------------------------------------+
| Package        | Version      | Details                                             |
+================+==============+=====================================================+
| **Python**     | 3.6+         |                                                     |
+----------------+--------------+-----------------------------------------------------+
| **NumPy**      | 1.16+        |                                                     |
+----------------+--------------+-----------------------------------------------------+
| **SciPy**      | 1.0+         | Lower versions may have missing features.           |
+----------------+--------------+-----------------------------------------------------+


In addition, there are several optional packages that provide additional functionality:

.. cssclass:: table-striped

+--------------------------+--------------+-----------------------------------------------------+
| Package                  | Version      | Details                                             |
+==========================+==============+=====================================================+
| ``matplotlib``           | 1.2.1+       | Needed for all visualisation tasks.                 |
+--------------------------+--------------+-----------------------------------------------------+
| ``cython``               | 0.29.20+     | Needed for compiling some time-dependent            |
|                          |              | Hamiltonians.                                       |
+--------------------------+--------------+-----------------------------------------------------+
| ``cvxpy``                | 1.0+         | Needed to calculate diamond norms.                  |
+--------------------------+--------------+-----------------------------------------------------+
| C++                      | GCC 4.7+,    | Needed for compiling Cython files, made when        |
| Compiler                 | MS VS 2015   | using string-format time-dependence.                |
+--------------------------+--------------+-----------------------------------------------------+
| ``pytest``,              | 5.3+         | For running the test suite.                         |
| ``pytest-rerunfailures`` |              |                                                     |
+--------------------------+--------------+-----------------------------------------------------+
| LaTeX                    | TeXLive 2009+| Needed if using LaTeX in matplotlib figures, or for |    
|                          |              | nice circuit drawings in IPython.                   |
+--------------------------+--------------+-----------------------------------------------------+

In addition, there are several additional packages that are not dependencies, but may give you a better programming experience.
`IPython <https://ipython.org/>`_ provides an improved text-based Python interpreter that is far more full-featured that the default interpreter, and runs in a terminal.
If you prefer a more graphical set-up, `Jupyter <https://jupyter.org/>`_ provides a notebook-style interface to mix code and mathematical notes together.
Alternatively, `Spyder <https://www.spyder-ide.org/>`_ is a free integrated development environment for Python, with several nice features for debugging code.
QuTiP will detect if it is being used within one of these richer environments, and various outputs will have enhanced formatting.

.. _install-with-conda:

Installing with conda
=====================

QuTiP is designed to work best when using the `Anaconda <https://www.anaconda.com/products/individual>`_ or `Intel <https://software.intel.com/en-us/python-distribution>`_ Python distributions that support the conda package management system.
It is still possible to use ``pip`` to install QuTiP while using conda, but uniformly using conda will make complete dependency management easier.

If you already have your conda environment set up, and have the ``conda-forge`` channel available, then you can install QuTiP using:

.. code-block:: bash

   conda install qutip

This will install the minimum set of dependences, but none of the optional packages.

.. _adding-conda-forge:

Adding the conda-forge channel
------------------------------

To install QuTiP from conda, you will need to add the conda-forge channel.
The following command adds this channel with lowest priority, so conda will still try and install all other packages normally:

.. code-block:: bash

   conda config --append channels conda-forge

If you want to change the order of your channels later, you can edit your ``.condarc`` (user home folder) file manually, but it is recommended to keep ``defaults`` as the highest priority.


.. _building-conda-environment:

New conda environments
----------------------

The default Anaconda environment has all the Python packages needed for running QuTiP installed already, so you will only need to add the ``conda-forge`` channel and then install the package.
If you have only installed Miniconda, or you want a completely clean virtual environment to install QuTiP in, the ``conda`` package manager provides a convenient way to do this.

To create a conda environment for QuTiP called ``qutip-env``:

.. code-block:: bash

   conda create -n qutip-env python qutip

This will automatically install all the necessary packages, and none of the optional packages.
You activate the new environment by running

.. code-block:: bash

   conda activate qutip-env

You can also install any more optional packages you want with ``conda install``, for example ``matplotlib``, ``ipython`` or ``jupyter``.

.. _install-from-source:

Installing from Source
======================

Official releases of QuTiP are available from the download section on `the project's web pages <https://qutip.org/download.html>`_, and the latest source code is available in `our GitHub repository <https://github.com/qutip/qutip>`_.
In general we recommend users to use the latest stable release of QuTiP, but if you are interested in helping us out with development or wish to submit bug fixes, then use the latest development version from the GitHub repository.

You can install from source by using the `Python-recommended PEP 517 procedure <build-pep517_>`_, or if you want more control or to have a development version, you can use the `low-level build procedure with setuptools <build-setuptools_>`_.

.. _build-pep517:

PEP 517 Source Builds
---------------------

The easiest way to build QuTiP from source is to use a PEP-517-compatible builder such as the ``build`` package available on ``pip``.
These will automatically install all build dependencies for you, and the ``pip`` installation step afterwards will install the minimum runtime dependencies.
You can do this by doing (for example)

.. code-block:: bash

   pip install build
   python -m build <path to qutip>
   pip install <path to qutip>/dist/qutip-<version>.whl

The first command installs the reference PEP-517 build tool, the second effects the build and the third uses ``pip`` to install the built package.
You will need to replace ``<path to qutip>`` with the actual path to the QuTiP source code.
The string ``<version>`` will depend on the version of QuTiP, the version of Python and your operating system.
It will look something like ``4.6.0-cp39-cp39-manylinux1_x86_64``, but there should only be one ``.whl`` file in the ``dist/`` directory, which will be the correct one.


.. _build-setuptools:

Direct Setuptools Source Builds
-------------------------------

This is the method to have the greatest amount of control over the installation, but it the most error-prone and not recommended unless you know what you are doing.
You first need to have all the runtime dependencies installed.
The most up-to-date requirements will be listed in ``pyproject.toml`` file, in the ``build-system.requires`` key.
As of the 4.6.0 release, the build requirements can be installed with

.. code-block:: bash

   pip install setuptools wheel packaging 'cython>=0.29.20' 'numpy>=1.16.6,<1.20' 'scipy>=1.0'

or similar with ``conda`` if you prefer.
You will also need to have a functional C++ compiler installed on your system.
This is likely already done for you if you are on Linux or macOS, but see the `section on Windows installations <install-on-windows_>`_ if that is your operating system.

To install QuTiP from the source code run:

.. code-block:: bash

   python setup.py install

To install OpenMP support, if available, run:

.. code-block:: bash

   python setup.py install --with-openmp

This will attempt to load up OpenMP libraries during the compilation process, which depends on you having suitable C++ compiler and library support.
If you are on Linux this is probably already done, but the compiler macOS ships with does not have OpenMP support.
You will likely need to refer to external operating-system-specific guides for more detail here, as it may be very non-trivial to correctly configure.
   
If you wish to contribute to the QuTiP project, then you will want to create your own fork of `the QuTiP git repository <https://github.com/qutip/qutip>`_, clone this to a local folder, and install it into your Python environment using:

.. code-block:: bash

   python setup.py develop

When you do ``import qutip`` in this environment, you will then load the code from your local fork, enabling you to edit the Python files and have the changes immediately available when you restart your Python interpreter, without needing to rebuild the package.
Note that if you change any Cython files, you will need to rerun the build command.

You should not need to use ``sudo`` (or other superuser privileges) to install into a personal virtual environment; if it feels like you need it, there is a good chance that you are installing into the system Python environment instead.


.. _install-on-windows:

Installation on Windows
=======================

As with other operating systems, the easiest method is to use ``pip install qutip``, or use the ``conda`` procedure described above.
If you want to build from source or use runtime compilation with Cython, you will need to have a working C++ compiler.

You can `download the Visual Studio IDE from Microsoft <https://visualstudio.microsoft.com/downloads/>`_, which has a free Community edition containing a sufficient C++ compiler.
This is the recommended compiler toolchain on Windows.
When installing, be sure to select the following components:

- Windows "X" SDK (where "X" stands for your version: 7/8/8.1/10)
- Visual Studio C++ build tools

You can then follow the `installation from source <install-from-source_>`_ section as normal.

.. important::

   In order to prevent issues with the ``PATH`` environment variable not containing the compiler and associated libraries, it is recommended to use the developer command prompt in the Visual Studio installation folder instead of the built-in command prompt.

The Community edition of Visual Studio takes around 10GB of disk space.
If this is prohibitive for you, it is also possible to install `only the build tools and necessary SDKs <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ instead, which should save about 2GB of space.


.. _install-verify:

Verifying the Installation
==========================

QuTiP includes a collection of built-in test scripts to verify that an installation was successful.
To run the suite of tests scripts you must also have the ``pytest`` testing library.
After installing QuTiP, leave the installation directory, run Python (or IPython), and call:

.. code-block:: python

   import qutip.testing
   qutip.testing.run()

This will take between 10 and 30 minutes, depending on your computer.
At the end, the testing report should report a success; it is normal for some tests to be skipped, and for some to be marked "xfail" in yellow.
Skips may be tests that do not run on your operating system, or tests of optional components that you have not installed the dependencies for.
If any failures or errors occur, please check that you have installed all of the required modules.
See the next section on how to check the installed versions of the QuTiP dependencies.
If these tests still fail, then head on over to the `QuTiP Discussion Board <https://groups.google.com/group/qutip>`_ or `the GitHub issues page <https://github.com/qutip/qutip/issues>`_ and post a message detailing your particular issue.

.. _install-about:

Checking Version Information
============================

QuTiP includes an "about" function for viewing information about QuTiP and the important dependencies installed on your system.
To view this information:

.. code-block:: python

   import qutip
   qutip.about()
.. figure:: figures/logo.png
   :align: center
   :width: 7in


QuTiP: Quantum Toolbox in Python
================================

.. toctree::
   :maxdepth: 3

   frontmatter.rst
   installation.rst
   guide/guide.rst
   gallery/build/index.rst
   apidoc/apidoc.rst

   changelog.rst
   contributors.rst
   development/development.rst
   biblio.rst
   copyright.rst


Indices and tables
====================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _biblo:
    
Bibliography
============

.. [BCSZ08]
    W. Bruzda, V. Cappellini, H.-J. Sommers, K. Życzkowski, *Random Quantum Operations*, Phys. Lett. A **373**, 320-324 (2009). :doi:`10.1016/j.physleta.2008.11.043`.

.. [Hav03]
    Havel, T. *Robust procedures for converting among Lindblad, Kraus and matrix representations of quantum dynamical semigroups*. Journal of Mathematical Physics **44** 2, 534 (2003). :doi:`10.1063/1.1518555`.

.. [Wat13]
    Watrous, J. |theory-qi|_, lecture notes.

..  The trick with |text|_ is to get an italic link, and is described in the
    Docutils FAQ at https://docutils.sourceforge.net/FAQ.html#is-nested-inline-markup-possible.
    
.. |theory-qi| replace:: *Theory of Quantum Information*
.. _theory-qi: https://cs.uwaterloo.ca/~watrous/CS766/

.. [Mez07]
    F. Mezzadri, *How to generate random matrices from the classical compact groups*, Notices of the AMS **54** 592-604 (2007). :arxiv:`math-ph/0609050`.

.. [Moh08]
    M. Mohseni, A. T. Rezakhani, D. A. Lidar, *Quantum-process tomography: Resource analysis of different strategies*, Phys. Rev. A **77**, 032322 (2008). :doi:`10.1103/PhysRevA.77.032322`.

.. [Gri98]
    M. Grifoni, P. Hänggi, *Driven quantum tunneling*, Physics Reports **304**, 299 (1998). :doi:`10.1016/S0370-1573(98)00022-2`.

.. [Gar03]
    Gardineer and Zoller, *Quantum Noise* (Springer, 2004).

.. [Bre02]
    H.-P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems* (Oxford, 2002).

.. [Coh92]
    C. Cohen-Tannoudji, J. Dupont-Roc, G. Grynberg, *Atom-Photon Interactions: Basic Processes and Applications*, (Wiley, 1992).

.. [WBC11]
    C. Wood, J. Biamonte, D. G. Cory, *Tensor networks and graphical calculus for
    open quantum systems*. :arxiv:`1111.6950`
    
.. [dAless08]
    D. d’Alessandro, *Introduction to Quantum Control and Dynamics*, (Chapman & Hall/CRC, 2008).
    
.. [Byrd95]
    R. H. Byrd, P. Lu, J. Nocedal, and C. Zhu, *A Limited Memory Algorithm for Bound Constrained Optimization*, SIAM J. Sci. Comput. **16**, 1190 (1995). :doi:`10.1137/0916069`

.. [Flo12]
    F. F. Floether, P. de Fouquieres, and S. G. Schirmer, *Robust quantum gates for open systems via optimal control: Markovian versus non-Markovian dynamics*, New J. Phys. **14**, 073023 (2012). :doi:`10.1088/1367-2630/14/7/073023`

.. [Lloyd14]
    S. Lloyd and S. Montangero, *Information theoretical analysis of quantum optimal control*, Phys. Rev. Lett. **113**, 010502 (2014). :doi:`10.1103/PhysRevLett.113.010502`
    
.. [Doria11]
    P. Doria, T. Calarco & S. Montangero, *Optimal Control Technique for Many-Body Quantum Dynamics*, Phys. Rev. Lett. **106**, 190501 (2011). :doi:`10.1103/PhysRevLett.106.190501`
    
.. [Caneva11]
    T. Caneva, T. Calarco, & S. Montangero, *Chopped random-basis quantum optimization*, Phys. Rev. A **84**, 022326 (2011). :doi:`10.1103/PhysRevA.84.022326`
    
.. [Rach15]
    N. Rach, M. M. Müller, T. Calarco, and S. Montangero, *Dressing the chopped-random-basis optimization: A bandwidth-limited access to the trap-free landscape*, Phys. Rev. A. **92**, 062343 (2015). :doi:`10.1103/PhysRevA.92.062343`

.. [DYNAMO]
    S. Machnes, U. Sander, S. J. Glaser, P. De Fouquieres, A. Gruslys, S. Schirmer, and T. Schulte-Herbrueggen, *Comparing, Optimising and Benchmarking Quantum Control Algorithms in a Unifying Programming Framework*, Phys. Rev. A. **84**, 022305 (2010). :arxiv:`1011.4874`

.. [Wis09]

    Wiseman, H. M. & Milburn, G. J. *Quantum Measurement and Control*, (Cambridge University Press, 2009).
.. _changelog:

**********
Change Log
**********

Version 4.6.2 (June 2, 2021)
++++++++++++++++++++++++++++

This minor release adds a function to calculate the quantum relative entropy, fixes a corner case in handling time-dependent Hamiltonians in ``mesolve`` and adds back support for a wider range of matplotlib versions when plotting or animating Bloch spheres.

It also adds a section in the README listing the papers which should be referenced while citing QuTiP. 
 

Improvements
------------
- Added a "Citing QuTiP" section to the README, containing a link to the QuTiP papers. (`#1554 <https://github.com/qutip/qutip/pull/1554>`_)
- Added ``entropy_relative`` which returns the quantum relative entropy between two density matrices. (`#1553 <https://github.com/qutip/qutip/pull/1553>`_)

Bug Fixes
---------
- Fixed Bloch sphere distortion when using Matplotlib >= 3.3.0. (`#1496  <https://github.com/qutip/qutip/pull/1496>`_)
- Removed use of integer-like floats in math.factorial since it is deprecated as of Python 3.9. (`#1550 <https://github.com/qutip/qutip/pull/1550>`_)
- Simplified call to ffmpeg used in the the Bloch sphere animation tutorial to work with recent versions of ffmpeg. (`#1557 <https://github.com/qutip/qutip/pull/1557>`_)
- Removed blitting in Bloch sphere FuncAnimation example. (`#1558 <https://github.com/qutip/qutip/pull/1558>`_)
- Added a version checking condition to handle specific functionalities depending on the matplotlib version. (`#1556 <https://github.com/qutip/qutip/pull/1556>`_)
- Fixed ``mesolve`` handling of time-dependent Hamiltonian with a custom tlist and ``c_ops``. (`#1561 <https://github.com/qutip/qutip/pull/1561>`_)

Developer Changes
-----------------
- Read documentation version and release from the VERSION file.


Version 4.6.1 (May 4, 2021)
+++++++++++++++++++++++++++

This minor release fixes bugs in QIP gate definitions, fixes building from
the source tarball when git is not installed and works around an MKL
bug in versions of SciPy <= 1.4.

It also adds the ``[full]`` pip install target so that ``pip install qutip[full]``
installs qutip and all of its optional and developer dependencies.

Improvements
------------
- Add the ``[full]`` pip install target (by **Jake Lishman**)

Bug Fixes
---------
- Work around pointer MKL eigh bug in SciPy <= 1.4 (by **Felipe Bivort Haiek**)
- Fix berkeley, swapalpha and cz gate operations (by **Boxi Li**)
- Expose the CPHASE control gate (by **Boxi Li**)
- Fix building from the sdist when git is not installed (by **Jake Lishman**)

Developer Changes
-----------------
- Move the qutip-doc documentation into the qutip repository (by **Jake Lishman**)
- Fix warnings in documentation build (by **Jake Lishman**)
- Fix warnings in pytest runs and make pytest treat warnings as errors (by **Jake Lishman**)
- Add Simon Cross as author (by **Simon Cross**)


Version 4.6.0 (April 11, 2021)
++++++++++++++++++++++++++++++

This release brings improvements for qubit circuits, including a pulse scheduler, measurement statistics, reading/writing OpenQASM and optimisations in the circuit simulations.

This is the first release to have full binary wheel releases on pip; you can now do ``pip install qutip`` on almost any machine to get a correct version of the package without needing any compilers set up.
The support for Numpy 1.20 that was first added in QuTiP 4.5.3 is present in this version as well, and the same build considerations mentioned there apply here too.
If building using the now-supported PEP 517 mechanisms (e.g. ``python -mbuild /path/to/qutip``), all build dependencies will be correctly satisfied.

Improvements
------------
- **MAJOR** Add saving, loading and resetting functionality to ``qutip.settings`` for easy re-configuration. (by **Eric Giguère**)
- **MAJOR** Add a quantum gate scheduler in ``qutip.qip.scheduler``, to help parallelise the operations of quantum gates.  This supports two scheduling modes: as late as possible, and as soon as possible. (by **Boxi Li**)
- **MAJOR** Improved qubit circuit simulators, including OpenQASM support and performance optimisations. (by **Sidhant Saraogi**)
- **MAJOR** Add tools for quantum measurements and their statistics. (by **Simon Cross** and **Sidhant Saraogi**)
- Add support for Numpy 1.20.  QuTiP should be compiled against a version of Numpy ``>= 1.16.6`` and ``< 1.20`` (note: does _not_ include 1.20 itself), but such an installation is compatible with any modern version of Numpy.  Source installations from ``pip`` understand this constraint.
- Improve the error message when circuit plotting fails. (by **Boxi Li**)
- Add support for parsing M1 Mac hardware information. (by **Xiaoliang Wu**)
- Add more single-qubit gates and controlled gates. (by **Mateo Laguna** and **Martín Sande Costa**)
- Support decomposition of ``X``, ``Y`` and ``Z`` gates in circuits. (by **Boxi Li**)
- Refactor ``QubitCircuit.resolve_gate()`` (by **Martín Sande Costa**)

Bug Fixes
---------
- Fix ``dims`` in the returns from ``Qobj.eigenstates`` on superoperators. (by **Jake Lishman**)
- Calling Numpy ufuncs on ``Qobj`` will now correctly raise a ``TypeError`` rather than returning a nonsense ``ndarray``. (by **Jake Lishman**)
- Convert segfault into Python exception when creating too-large tensor products. (by **Jake Lishman**)
- Correctly set ``num_collapse`` in the output of ``mesolve``. (by **Jake Lishman**)
- Fix ``ptrace`` when all subspaces are being kept, or the subspaces are passed in order. (by **Jake Lishman**)
- Fix sorting bug in ``Bloch3d.add_points()``. (by **pschindler**)
- Fix invalid string literals in docstrings and some unclosed files. (by **Élie Gouzien**)
- Fix Hermicity tests for matrices with values that are within the tolerance of 0. (by **Jake Lishman**)
- Fix the trace norm being incorrectly reported as 0 for small matrices. (by **Jake Lishman**)
- Fix issues with ``dnorm`` when using CVXPy 1.1 with sparse matrices. (by **Felipe Bivort Haiek**)
- Fix segfaults in ``mesolve`` when passed a bad initial ``Qobj`` as the state. (by **Jake Lishman**)
- Fix sparse matrix construction in PIQS when using Scipy 1.6.1. (by **Drew Parsons**)
- Fix ``zspmv_openmp.cpp`` missing from the pip sdist. (by **Christoph Gohlke**)
- Fix correlation functions throwing away imaginary components. (by **Asier Galicia Martinez**)
- Fix ``QubitCircuit.add_circuit()`` for SWAP gate. (by **Canoming**)
- Fix the broken LaTeX image conversion. (by **Jake Lishman**)
- Fix gate resolution of the FREDKIN gate. (by **Bo Yang**)
- Fix broken formatting in docstrings. (by **Jake Lishman**)

Deprecations
------------
- ``eseries``, ``essolve`` and ``ode2es`` are all deprecated, pending removal in QuTiP 5.0.  These are legacy functions and classes that have been left unmaintained for a long time, and their functionality is now better achieved with ``QobjEvo`` or ``mesolve``.

Developer Changes
-----------------
- **MAJOR** Overhaul of setup and packaging code to make it satisfy PEP 517, and move the build to a matrix on GitHub Actions in order to release binary wheels on pip for all major platforms and supported Python versions. (by **Jake Lishman**)
- Default arguments in ``Qobj`` are now ``None`` rather than mutable types. (by **Jake Lishman**)
- Fixed comsumable iterators being used to parametrise some tests, preventing the testing suite from being re-run within the same session. (by **Jake Lishman**)
- Remove unused imports, simplify some floats and remove unnecessary list conversions. (by **jakobjakobson13**)
- Improve Travis jobs matrix for specifying the testing containers. (by **Jake Lishman**)
- Fix coverage reporting on Travis. (by **Jake Lishman**)
- Added a ``pyproject.toml`` file. (by **Simon Humpohl** and **Eric Giguère**)
- Add doctests to documentation. (by **Sidhant Saraogi**)
- Fix all warnings in the documentation build. (by **Jake Lishman**)



Version 4.5.3 (February 19, 2021)
+++++++++++++++++++++++++++++++++

This patch release adds support for Numpy 1.20, made necessary by changes to how array-like objects are handled. There are no other changes relative to version 4.5.2.

Users building from source should ensure that they build against Numpy versions >= 1.16.6 and < 1.20 (not including 1.20 itself), but after that or for those installing from conda, an installation will support any current Numpy version >= 1.16.6.

Improvements
------------
- Add support for Numpy 1.20.  QuTiP should be compiled against a version of Numpy ``>= 1.16.6`` and ``< 1.20`` (note: does _not_ include 1.20 itself), but such an installation is compatible with any modern version of Numpy.  Source installations from ``pip`` understand this constraint.



Version 4.5.2 (July 14, 2020)
+++++++++++++++++++++++++++++

This is predominantly a hot-fix release to add support for Scipy 1.5, due to changes in private sparse matrix functions that QuTiP also used.

Improvements
------------
- Add support for Scipy 1.5. (by **Jake Lishman**)
- Improved speed of ``zcsr_inner``, which affects ``Qobj.overlap``. (by **Jake Lishman**)
- Better error messages when installation requirements are not satisfied. (by **Eric Giguère**)

Bug Fixes
---------
- Fix ``zcsr_proj`` acting on matrices with unsorted indices.  (by **Jake Lishman**)
- Fix errors in Milstein's heterodyne. (by **Eric Giguère**)
- Fix datatype bug in ``qutip.lattice`` module. (by **Boxi Li**)
- Fix issues with ``eigh`` on Mac when using OpenBLAS.  (by **Eric Giguère**)

Developer Changes
-----------------
- Converted more of the codebase to PEP 8.
- Fix several instances of unsafe mutable default values and unsafe ``is`` comparisons.



Version 4.5.1 (May 15, 2020)
++++++++++++++++++++++++++++

Improvements
------------
- ``husimi`` and ``wigner`` now accept half-integer spin (by **maij**)
- Better error messages for failed string coefficient compilation. (issue raised by **nohchangsuk**)

Bug Fixes
---------
- Safer naming for temporary files. (by **Eric Giguère**)
- Fix ``clebsch`` function for half-integer (by **Thomas Walker**)
- Fix ``randint``'s dtype to ``uint32`` for compatibility with Windows. (issue raised by **Boxi Li**)
- Corrected stochastic's heterodyne's m_ops (by **eliegenois**)
- Mac pool use spawn. (issue raised by **goerz**)
- Fix typos in ``QobjEvo._shift``. (by **Eric Giguère**)
- Fix warning on Travis CI. (by **Ivan Carvalho**)

Deprecations
------------
- ``qutip.graph`` functions will be deprecated in QuTiP 5.0 in favour of ``scipy.sparse.csgraph``.

Developer Changes
-----------------
- Add Boxi Li to authors. (by **Alex Pitchford**)
- Skip some tests that cause segfaults on Mac. (by **Nathan Shammah** and **Eric Giguère**)
- Use Python 3.8 for testing on Mac and Linux. (by **Simon Cross** and **Eric Giguère**)



Version 4.5.0 (January 31, 2020)
++++++++++++++++++++++++++++++++

Improvements
------------
- **MAJOR FEATURE**: Added `qip.noise`, a module with pulse level description of quantum circuits allowing to model various types of noise and devices (by **Boxi Li**).

- **MAJOR FEATURE**: Added `qip.lattice`, a module for the study of lattice dynamics in 1D (by **Saumya Biswas**).

- Migrated testing from Nose to PyTest (by **Tarun Raheja**).

- Optimized testing for PyTest and removed duplicated test runners (by **Jake Lishman**).

- Deprecated importing `qip` functions to the qutip namespace (by **Boxi Li**).

- Added the possibility to define non-square superoperators relevant for quantum circuits (by **Arne Grimsmo** and **Josh Combes**).

- Implicit tensor product for `qeye`, `qzero` and `basis` (by **Jake Lishman**).

- QObjEvo no longer requires Cython for string coefficient (by **Eric Giguère**).

- Added marked tests for faster tests in `testing.run()` and made faster OpenMP benchmarking in CI (by **Eric Giguère**).

- Added entropy and purity for Dicke density matrices, refactored into more general dicke_trace (by **Nathan Shammah**).

- Added option for specifying resolution in Bloch.save function (by **Tarun Raheja**).

- Added information related to the value of hbar in `wigner` and `continuous_variables` (by **Nicolas Quesada**).

- Updated requirements for `scipy 1.4` (by **Eric Giguère**).

- Added previous lead developers to the qutip.about() message (by **Nathan Shammah**).

- Added improvements to `Qobj` introducing the `inv` method and making the partial trace, `ptrace`, faster, keeping both sparse and dense methods (by **Eric Giguère**).

- Allowed general callable objects to define a time-dependent Hamiltonian (by **Eric Giguère**).

- Added feature so that `QobjEvo` no longer requires Cython for string coefficients (by **Eric Giguère**).

- Updated authors list on Github and added `my binder` link (by **Nathan Shammah**).


Bug Fixes
---------

- Fixed `PolyDataMapper` construction for `Bloch3d` (by **Sam Griffiths**).

- Fixed error checking for null matrix in essolve (by **Nathan Shammah**).

- Fixed name collision for parallel propagator (by **Nathan Shammah**).

- Fixed dimensional incongruence in `propagator` (by **Nathan Shammah**)

- Fixed bug by rewriting clebsch function based on long integer fraction (by **Eric Giguère**).

- Fixed bugs in QobjEvo's args depending on state and added solver tests using them (by **Eric Giguère**).

- Fixed bug in `sesolve` calculation of average states when summing the timeslot states (by **Alex Pitchford**).

- Fixed bug in `steadystate` solver by removing separate arguments for MKL and Scipy (by **Tarun Raheja**).

- Fixed `Bloch.add_ponts` by setting `edgecolor = None` in `plot_points` (by **Nathan Shammah**).

- Fixed error checking for null matrix in `essolve` solver affecting also `ode2es` (by **Peter Kirton**).

- Removed unnecessary shebangs in .pyx and .pxd files (by **Samesh Lakhotia**).

- Fixed `sesolve` and  import of `os` in `codegen` (by **Alex Pitchford**).

- Updated `plot_fock_distribution` by removing the offset value 0.4 in the plot (by **Rajiv-B**).


Version 4.4.1 (August 29, 2019)
+++++++++++++++++++++++++++++++

Improvements
------------

- QobjEvo do not need to start from 0 anymore (by **Eric Giguère**).

- Add a quantum object purity function (by **Nathan Shammah** and **Shahnawaz Ahmed**).

- Add step function interpolation for array time-coefficient (by **Boxi Li**).

- Generalize expand_oper for arbitrary dimensions, and new method for cyclic permutations of given target cubits (by **Boxi Li**).


Bug Fixes
---------

- Fixed the pickling but that made solver unable to run in parallel on Windows (Thank **lrunze** for reporting)

- Removed warning when mesolve fall back on sesolve (by **Michael Goerz**).

- Fixed dimension check and confusing documentation in random ket (by **Yariv Yanay**).

- Fixed Qobj isherm not working after using Qobj.permute (Thank **llorz1207** for reporting).

- Correlation functions call now properly handle multiple time dependant functions (Thank **taw181** for reporting).

- Removed mutable default values in mesolve/sesolve (by **Michael Goerz**).

- Fixed simdiag bug (Thank **Croydon-Brixton** for reporting).

- Better support of constant QobjEvo (by **Boxi Li**).

- Fixed potential cyclic import in the control module (by **Alexander Pitchford**).


Version 4.4.0 (July 03, 2019)
+++++++++++++++++++++++++++++

Improvements
------------

- **MAJOR FEATURE**: Added methods and techniques to the stochastic solvers (by **Eric Giguère**) which allows to use a much broader set of solvers and much more efficiently.

- **MAJOR FEATURE**: Optimization of the montecarlo solver (by **Eric Giguère**). Computation are faster in many cases. Collapse information available to time dependant information.

- Added the QObjEvo class and methods (by **Eric Giguère**), which is used behind the scenes by the dynamical solvers, making the code more efficient and tidier. More built-in function available to string coefficients.

- The coefficients can be made from interpolated array with variable timesteps and can obtain state information more easily. Time-dependant collapse operator can have multiple terms.

- New wigner_transform and plot_wigner_sphere function. (by **Nithin Ramu**).

- ptrace is faster and work on bigger systems, from 15 Qbits to 30 Qbits.

- QIP module: added the possibility for user-defined gates, added the possibility to remove or add gates in any point of an already built circuit, added the molmer_sorensen gate, and fixed some bugs (by **Boxi Li**).

- Added the quantum Hellinger distance to qutip.metrics (by **Wojciech Rzadkowski**).

- Implemented possibility of choosing a random seed (by **Marek Marekyggdrasil**).

- Added a code of conduct to Github.


Bug Fixes
---------

- Fixed bug that made QuTiP incompatible with SciPy 1.3.


Version 4.3.0 (July 14, 2018)
+++++++++++++++++++++++++++++

Improvements
------------

- **MAJOR FEATURE**: Added the Permutational Invariant Quantum Solver (PIQS) module (by **Nathan Shammah** and **Shahnawaz Ahmed**) which allows the simluation of large TLSs ensembles including collective and local Lindblad dissipation. Applications range from superradiance to spin squeezing.

- **MAJOR FEATURE**: Added a photon scattering module (by **Ben Bartlett**) which can be used to study scattering in arbitrary driven systems coupled to some configuration of output waveguides.

- Cubic_Spline functions as time-dependent arguments for the collapse operators in mesolve are now allowed.

- Added a faster version of bloch_redfield_tensor, using components from the time-dependent version. About 3x+ faster for secular tensors, and 10x+ faster for non-secular tensors.

- Computing Q.overlap() [inner product] is now ~30x faster.

- Added projector method to Qobj class.

- Added fast projector method, ``Q.proj()``.

- Computing matrix elements, ``Q.matrix_element`` is now ~10x faster.

- Computing expectation values for ket vectors using ``expect`` is now ~10x faster.

- ``Q.tr()`` is now faster for small Hilbert space dimensions.

- Unitary operator evolution added to sesolve

- Use OPENMP for tidyup if installed.


Bug Fixes
---------

- Fixed bug that stopped simdiag working for python 3.

- Fixed semidefinite cvxpy Variable and Parameter.

- Fixed iterative lu solve atol keyword issue.

- Fixed unitary op evolution rhs matrix in ssesolve.

- Fixed interpolating function to return zero outside range.

- Fixed dnorm complex casting bug.

- Fixed control.io path checking issue.

- Fixed ENR fock dimension.

- Fixed hard coded options in propagator 'batch' mode

- Fixed bug in trace-norm for non-Hermitian operators.

- Fixed bug related to args not being passed to coherence_function_g2

- Fixed MKL error checking dict key error


Version 4.2.0 (July 28, 2017)
+++++++++++++++++++++++++++++

Improvements
------------

- **MAJOR FEATURE**: Initial implementation of time-dependent Bloch-Redfield Solver.

- Qobj tidyup is now an order of magnitude faster.

- Time-dependent codegen now generates output NumPy arrays faster.

- Improved calculation for analytic coefficients in coherent states (Sebastian Kramer).

- Input array to correlation FFT method now checked for validity.

- Function-based time-dependent mesolve and sesolve routines now faster.

- Codegen now makes sure that division is done in C, as opposed to Python.

- Can now set different controls for a each timeslot in quantum optimization.
  This allows time-varying controls to be used in pulse optimisation.


Bug Fixes
---------

- rcsolve importing old Odeoptions Class rather than Options.

- Non-int issue in spin Q and Wigner functions.

- Qobj's should tidyup before determining isherm.

- Fixed time-dependent RHS function loading on Win.

- Fixed several issues with compiling with Cython 0.26.

- Liouvillian superoperators were hard setting isherm=True by default.

- Fixed an issue with the solver safety checks when inputing a list with Python functions as time-dependence.

- Fixed non-int issue in Wigner_cmap.

- MKL solver error handling not working properly.



Version 4.1.0 (March 10, 2017)
++++++++++++++++++++++++++++++

Improvements
------------

*Core libraries*

- **MAJOR FEATURE**: QuTiP now works for Python 3.5+ on Windows using Visual Studio 2015.

- **MAJOR FEATURE**: Cython and other low level code switched to C++ for MS Windows compatibility.

- **MAJOR FEATURE**: Can now use interpolating cubic splines as time-dependent coefficients.

- **MAJOR FEATURE**: Sparse matrix - vector multiplication now parallel using OPENMP.

- Automatic tuning of OPENMP threading threshold.

- Partial trace function is now up to 100x+ faster.

- Hermitian verification now up to 100x+ faster.

- Internal Qobj objects now created up to 60x faster.

- Inplace conversion from COO -> CSR sparse formats (e.g. Memory efficiency improvement.)

- Faster reverse Cuthill-Mckee and sparse one and inf norms.



Bug Fixes
---------

- Cleanup of temp. Cython files now more robust and working under Windows.



Version 4.0.2 (January 5, 2017)
+++++++++++++++++++++++++++++++

Bug Fixes
---------
- td files no longer left behind by correlation tests
- Various fast sparse fixes



Version 4.0.0 (December 22, 2016)
+++++++++++++++++++++++++++++++++

Improvements
------------
*Core libraries*

- **MAJOR FEATURE**: Fast sparse: New subclass of csr_matrix added that overrides commonly used methods to avoid certain checks that incurr execution cost. All Qobj.data now fast_csr_matrix
- HEOM performance enhancements
- spmv now faster
- mcsolve codegen further optimised

*Control modules*

- Time dependent drift (through list of pwc dynamics generators)
- memory optimisation options provided for control.dynamics

Bug Fixes
---------

- recompilation of pyx files on first import removed
- tau array in control.pulseoptim funcs now works

Version 3.2.0 (Never officially released)
+++++++++++++++++++++++++++++++++++++++++

New Features
------------

*Core libraries*

- **MAJOR FEATURE**: Non-Markovian solvers: Hierarchy (**Added by Neill Lambert**), Memory-Cascade, and Transfer-Tensor methods.
- **MAJOR FEATURE**: Default steady state solver now up to 100x faster using the Intel Pardiso library under the Anaconda and Intel Python distributions.
- The default Wigner function now uses a Clenshaw summation algorithm to evaluate a polynomial series that is applicable for any number of exciations (previous limitation was ~50 quanta), and is ~3x faster than before. (**Added by Denis Vasilyev**)
- Can now define a given eigen spectrum for random Hermitian and density operators.
- The Qobj ``expm`` method now uses the equivilent SciPy routine, and performs a much faster ``exp`` operation if the matrix is diagonal.
- One can now build zero operators using the ``qzero`` function.

*Control modules*

- **MAJOR FEATURE**: CRAB algorithm added
  This is an alternative to the GRAPE algorithm, which allows for analytical control functions, which means that experimental constraints can more easily be added into optimisation.
  See tutorial notebook for full information.


Improvements
------------
*Core libraries*

- Two-time correlation functions can now be calculated for fully time-dependent Hamiltonians and collapse operators. (**Added by Kevin Fischer**)
- The code for the inverse-power method for the steady state solver has been simplified.
- Bloch-Redfield tensor creation is now up to an order of magnitude faster. (**Added by Johannes Feist**)
- Q.transform now works properly for arrays directly from sp_eigs (or eig).
- Q.groundstate now checks for degeneracy.
- Added ``sinm`` and ``cosm`` methods to the Qobj class.
- Added ``charge`` and ``tunneling`` operators.
- Time-dependent Cython code is now easier to read and debug.


*Control modules*

- The internal state / quantum operator data type can now be either Qobj or ndarray
  Previous only ndarray was possible. This now opens up possibility of using Qobj methods in fidelity calculations
  The attributes and functions that return these operators are now preceded by an underscore, to indicate that the data type could change depending on the configuration options.
  In most cases these functions were for internal processing only anyway, and should have been 'private'.
  Accessors to the properties that could be useful outside of the library have been added. These always return Qobj. If the internal operator data type is not Qobj, then there could be signicant overhead in the conversion, and so this should be avoided during pulse optimisation.
  If custom sub-classes are developed that use Qobj properties and methods (e.g. partial trace), then it is very likely that it will be more efficient to set the internal data type to Qobj.
  The internal operator data will be chosen automatically based on the size and sparsity of the dynamics generator. It can be forced by setting ``dynamics.oper_dtype = <type>``
  Note this can be done by passing ``dyn_params={'oper_dtype':<type>}`` in any of the pulseoptim functions.

  Some other properties and methods were renamed at the same time. A full list is given here.

  - All modules
    - function: ``set_log_level`` -> property: ``log_level``

  - dynamics functions

    - ``_init_lists`` now ``_init_evo``
    - ``get_num_ctrls`` now property: ``num_ctrls``
    - ``get_owd_evo_target`` now property: ``onto_evo_target``
    - ``combine_dyn_gen`` now ``_combine_dyn_gen`` (no longer returns a value)
    - ``get_dyn_gen`` now ``_get_phased_dyn_gen``
    - ``get_ctrl_den_gen`` now ``_get_phased_ctrl_dyn_gen``
    - ``ensure_decomp_curr`` now ``_ensure_decomp_curr``
    - ``spectral_decomp`` now ``_spectral_decomp``

  - dynamics properties

    - ``evo_init2t`` now ``_fwd_evo`` (``fwd_evo`` as Qobj)
    - ``evo_t2end`` now ``_onwd_evo`` (``onwd_evo`` as Qobj)
    - ``evo_t2targ`` now ``_onto_evo`` (``onto_evo`` as Qobj)

  - fidcomp properties

    - ``uses_evo_t2end`` now ``uses_onwd_evo``
    - ``uses_evo_t2targ`` now ``uses_onto_evo``
    - ``set_phase_option`` function now property ``phase_option``

  - propcomp properties

    - ``grad_exact`` (now read only)

  - propcomp functions

    - ``compute_propagator`` now ``_compute_propagator``
    - ``compute_diff_prop`` now ``_compute_diff_prop``
    - ``compute_prop_grad`` now ``_compute_prop_grad``

  - tslotcomp functions

    - ``get_timeslot_for_fidelity_calc`` now ``_get_timeslot_for_fidelity_calc``


*Miscellaneous*

- QuTiP Travis CI tests now use the Anaconda distribution.
- The ``about`` box and ipynb ``version_table`` now display addition system information.
- Updated Cython cleanup to remove depreciation warning in sysconfig.
- Updated ipynb_parallel to look for ``ipyparallel`` module in V4 of the notebooks.


Bug Fixes
---------
- Fixes for countstat and psuedo-inverse functions
- Fixed Qobj division tests on 32-bit systems.
- Removed extra call to Python in time-dependent Cython code.
- Fixed issue with repeated Bloch sphere saving.
- Fixed T_0 triplet state not normalized properly. (**Fixed by Eric Hontz**)
- Simplified compiler flags (support for ARM systems).
- Fixed a decoding error in ``qload``.
- Fixed issue using complex.h math and np.kind_t variables.
- Corrected output states mismatch for ``ntraj=1`` in the mcf90 solver.
- Qobj data is now copied by default to avoid a bug in multiplication. (**Fixed by Richard Brierley**)
- Fixed bug overwriting ``hardware_info`` in ``__init__``. (**Fixed by Johannes Feist**)
- Restored ability to explicity set Q.isherm, Q.type, and Q.superrep.
- Fixed integer depreciation warnings from NumPy.
- Qobj * (dense vec) would result in a recursive loop.
- Fixed args=None -> args={} in correlation functions to be compatible with mesolve.
- Fixed depreciation warnings in mcsolve.
- Fixed neagtive only real parts in ``rand_ket``.
- Fixed a complicated list-cast-map-list antipattern in super operator reps. (**Fixed by Stefan Krastanov**)
- Fixed incorrect ``isherm`` for ``sigmam`` spin operator.
- Fixed the dims when using ``final_state_output`` in ``mesolve`` and ``sesolve``.



Version 3.1.0 (January 1, 2015)
+++++++++++++++++++++++++++++++

New Features
------------

- **MAJOR FEATURE**: New module for quantum control (qutip.control).
- **NAMESPACE CHANGE**: QuTiP no longer exports symbols from NumPy and matplotlib, so those modules must now be explicitly imported when required.
- New module for counting statistics.
- Stochastic solvers now run trajectories in parallel.
- New superoperator and tensor manipulation functions
  (super_tensor, composite, tensor_contract).
- New logging module for debugging (qutip.logging).
- New user-available API for parallelization (parallel_map).
- New enhanced (optional) text-based progressbar (qutip.ui.EnhancedTextProgressBar)
- Faster Python based monte carlo solver (mcsolve).
- Support for progress bars in propagator function.
- Time-dependent Cython code now calls complex cmath functions.
- Random numbers seeds can now be reused for successive calls to mcsolve.
- The Bloch-Redfield master equation solver now supports optional Lindblad type collapse operators.
- Improved handling of ODE integration errors in mesolve.
- Improved correlation function module (for example, improved support for time-dependent problems).
- Improved parallelization of mcsolve (can now be interrupted easily, support for IPython.parallel, etc.)
- Many performance improvements, and much internal code restructuring.

Bug Fixes
---------

- Cython build files for time-dependent string format now removed automatically.
- Fixed incorrect solution time from inverse-power method steady state solver.
- mcsolve now supports `Options(store_states=True)`
- Fixed bug in `hadamard` gate function.
- Fixed compatibility issues with NumPy 1.9.0.
- Progressbar in mcsolve can now be suppressed.
- Fixed bug in `gate_expand_3toN`.
- Fixed bug for time-dependent problem (list string format) with multiple terms in coefficient to an operator.

Version 3.0.1 (Aug 5, 2014)
+++++++++++++++++++++++++++

Bug Fixes
---------

- Fix bug in create(), which returned a Qobj with CSC data instead of CSR.
- Fix several bugs in mcsolve: Incorrect storing of collapse times and collapse
  operator records. Incorrect averaging of expectation values for different
  trajectories when using only 1 CPU.
- Fix bug in parsing of time-dependent Hamiltonian/collapse operator arguments
  that occurred when the args argument is not a dictionary.
- Fix bug in internal _version2int function that cause a failure when parsingthe version number of the Cython package.
-


Version 3.0.0 (July 17, 2014)
+++++++++++++++++++++++++++++

New Features
------------

- New module `qutip.stochastic` with stochastic master equation and stochastic
  Schrödinger equation solvers.

- Expanded steady state solvers. The function ``steady`` has been deprecated in
  favor of ``steadystate``. The steadystate solver no longer use umfpack by
  default. New pre-processing methods for reordering and balancing the linear
  equation system used in direct solution of the steady state.

- New module `qutip.qip` with utilities for quantum information processing,
  including pre-defined quantum gates along with functions for expanding
  arbitrary 1, 2, and 3 qubit gates to N qubit registers, circuit
  representations, library of quantum algorithms, and basic physical models for
  some common QIP architectures.

- New module `qutip.distributions` with unified API for working with
  distribution functions.

- New format for defining time-dependent Hamiltonians and collapse operators,
  using a pre-calculated numpy array that specifies the values of the
  Qobj-coefficients for each time step.

- New functions for working with different superoperator representations,
  including Kraus and Chi representation.

- New functions for visualizing quantum states using Qubism and Schimdt plots:
  ``plot_qubism`` and ``plot_schmidt``.

- Dynamics solver now support taking argument ``e_ops`` (expectation value
  operators) in dictionary form.

- Public plotting functions from the ``qutip.visualization`` module are now
  prefixed with ``plot_`` (e.g., ``plot_fock_distribution``). The
  ``plot_wigner`` and ``plot_wigner_fock_distribution`` now supports 3D views
  in addition to contour views.

- New API and new functions for working with spin operators and states,
  including for example ``spin_Jx``, ``spin_Jy``, ``spin_Jz`` and
  ``spin_state``, ``spin_coherent``.

- The ``expect`` function now supports a list of operators, in addition to the
  previously supported list of states.

- Simplified creation of qubit states using ``ket`` function.

- The module ``qutip.cyQ`` has been renamed to ``qutip.cy`` and the sparse
  matrix-vector functions ``spmv`` and ``spmv1d`` has been combined into one
  function ``spmv``. New functions for operating directly on the underlaying
  sparse CSR data have been added (e.g., ``spmv_csr``). Performance
  improvements. New and improved Cython functions for calculating expectation
  values for state vectors, density matrices in matrix and vector form.

- The ``concurrence`` function now supports both pure and mixed states. Added
  function for calculating the entangling power of a two-qubit gate.

- Added function for generating (generalized) Lindblad dissipator
  superoperators.

- New functions for generating Bell states, and singlet and triplet states.

- QuTiP no longer contains the demos GUI. The examples are now available on the
  QuTiP web site. The ``qutip.gui`` module has been renamed to ``qutip.ui`` and
  does no longer contain graphical UI elements. New text-based and HTML-based
  progressbar classes.

- Support for harmonic oscillator operators/states in a Fock state basis that
  does not start from zero (e.g., in the range [M,N+1]). Support for
  eliminating and extracting states from Qobj instances (e.g., removing one
  state from a two-qubit system to obtain a three-level system).

- Support for time-dependent Hamiltonian and Liouvillian callback functions that
  depend on the instantaneous state, which for example can be used for solving
  master equations with mean field terms.

Improvements
------------

- Restructured and optimized implementation of Qobj, which now has
  significantly lower memory footprint due to avoiding excessive copying of
  internal matrix data.

- The classes ``OdeData``, ``Odeoptions``, ``Odeconfig`` are now called
  ``Result``, ``Options``, and ``Config``, respectively, and are available in
  the module `qutip.solver`.

- The ``squeez`` function has been renamed to ``squeeze``.

- Better support for sparse matrices when calculating propagators using the
  ``propagator`` function.

- Improved Bloch sphere.

- Restructured and improved the module ``qutip.sparse``, which now only
  operates directly on sparse matrices (not on Qobj instances).

- Improved and simplified implement of the ``tensor`` function.

- Improved performance, major code cleanup (including namespace changes),
  and numerous bug fixes.

- Benchmark scripts improved and restructured.

- QuTiP is now using continuous integration tests (TravisCI).

Version 2.2.0 (March 01, 2013)
++++++++++++++++++++++++++++++


New Features
------------

- **Added Support for Windows**

- New Bloch3d class for plotting 3D Bloch spheres using Mayavi.

- Bloch sphere vectors now look like arrows.

- Partial transpose function.

- Continuos variable functions for calculating correlation and covariance
  matrices, the Wigner covariance matrix and the logarithmic negativity for
  for multimode fields in Fock basis.

- The master-equation solver (mesolve) now accepts pre-constructed Liouvillian
  terms, which makes it possible to solve master equations that are not on
  the standard Lindblad form.

- Optional Fortran Monte Carlo solver (mcsolve_f90) by Arne Grimsmo.

- A module of tools for using QuTiP in IPython notebooks.

- Increased performance of the steady state solver.

- New Wigner colormap for highlighting negative values.

- More graph styles to the visualization module.


Bug Fixes
---------

- Function based time-dependent Hamiltonians now keep the correct phase.

- mcsolve no longer prints to the command line if ntraj=1.


Version 2.1.0 (October 05, 2012)
++++++++++++++++++++++++++++++++


New Features
------------

- New method for generating Wigner functions based on Laguerre polynomials.

- coherent(), coherent_dm(), and thermal_dm() can now be expressed using analytic values.

- Unittests now use nose and can be run after installation.

- Added iswap and sqrt-iswap gates.

- Functions for quantum process tomography.

- Window icons are now set for Ubuntu application launcher.

- The propagator function can now take a list of times as argument, and returns a list of corresponding propagators.


Bug Fixes
---------

- mesolver now correctly uses the user defined rhs_filename in Odeoptions().

- rhs_generate() now handles user defined filenames properly.

- Density matrix returned by propagator_steadystate is now Hermitian.

- eseries_value returns real list if all imag parts are zero.

- mcsolver now gives correct results for strong damping rates.

- Odeoptions now prints mc_avg correctly.

- Do not check for PyObj in mcsolve when gui=False.

- Eseries now correctly handles purely complex rates.

- thermal_dm() function now uses truncated operator method.

- Cython based time-dependence now Python 3 compatible.

- Removed call to NSAutoPool on mac systems.

- Progress bar now displays the correct number of CPU's used.

- Qobj.diag() returns reals if operator is Hermitian.

- Text for progress bar on Linux systems is no longer cutoff.


Version 2.0.0 (June 01, 2012)
+++++++++++++++++++++++++++++

The second version of QuTiP has seen many improvements in the performance of the original code base, as well as the addition of several new routines supporting a wide range of functionality.  Some of the highlights of this release include:

New Features
------------

- QuTiP now includes solvers for both Floquet and Bloch-Redfield master equations.

- The Lindblad master equation and Monte Carlo solvers allow for time-dependent collapse operators.

- It is possible to automatically compile time-dependent problems into c-code using Cython (if installed).

- Python functions can be used to create arbitrary time-dependent Hamiltonians and collapse operators.

- Solvers now return Odedata objects containing all simulation results and parameters, simplifying the saving of simulation results.

.. important:: This breaks compatibility with QuTiP version 1.x.

- mesolve and mcsolve can reuse Hamiltonian data when only the initial state, or time-dependent arguments, need to be changed.

- QuTiP includes functions for creating random quantum states and operators.

- The generation and manipulation of quantum objects is now more efficient.

- Quantum objects have basis transformation and matrix element calculations as built-in methods.

- The quantum object eigensolver can use sparse solvers.

- The partial-trace (ptrace) function is up to 20x faster.

- The Bloch sphere can now be used with the Matplotlib animation function, and embedded as a subplot in a figure.

- QuTiP has built-in functions for saving quantum objects and data arrays.

- The steady-state solver has been further optimized for sparse matrices, and can handle much larger system Hamiltonians.

- The steady-state solver can use the iterative bi-conjugate gradient method instead of a direct solver.

- There are three new entropy functions for concurrence, mutual information, and conditional entropy.

- Correlation functions have been combined under a single function.

- The operator norm can now be set to trace, Frobius, one, or max norm.

- Global QuTiP settings can now be modified.

- QuTiP includes a collection of unit tests for verifying the installation.

- Demos window now lets you copy and paste code from each example.


Version 1.1.4 (May 28, 2012)
++++++++++++++++++++++++++++

Bug Fixes
---------

- Fixed bug pointed out by Brendan Abolins.

- Qobj.tr() returns zero-dim ndarray instead of float or complex.

- Updated factorial import for scipy version 0.10+


Version 1.1.3 (November 21, 2011)
+++++++++++++++++++++++++++++++++

New Functions
-------------

- Allow custom naming of Bloch sphere.

Bug Fixes
---------
- Fixed text alignment issues in AboutBox.

- Added fix for SciPy V>0.10 where factorial was moved to scipy.misc module.

- Added tidyup function to tensor function output.

- Removed openmp flags from setup.py as new Mac Xcode compiler does not recognize them.

- Qobj diag method now returns real array if all imaginary parts are zero.

- Examples GUI now links to new documentation.

- Fixed zero-dimensional array output from metrics module.


Version 1.1.2 (October 27, 2011)
++++++++++++++++++++++++++++++++

Bug Fixes
---------

- Fixed issue where Monte Carlo states were not output properly.


Version 1.1.1 (October 25, 2011)
++++++++++++++++++++++++++++++++

**THIS POINT-RELEASE INCLUDES VASTLY IMPROVED TIME-INDEPENDENT MCSOLVE AND ODESOLVE PERFORMANCE**

New Functions
-------------

- Added linear entropy function.

- Number of CPU's can now be changed.

Bug Fixes
---------

- Metrics no longer use dense matrices.

- Fixed Bloch sphere grid issue with matplotlib 1.1.

- Qobj trace operation uses only sparse matrices.

- Fixed issue where GUI windows do not raise to front.


Version 1.1.0 (October 04, 2011)
++++++++++++++++++++++++++++++++

**THIS RELEASE NOW REQUIRES THE GCC COMPILER TO BE INSTALLED**

New Functions
-------------

- tidyup function to remove small elements from a Qobj.

- Added concurrence function.

- Added simdiag for simultaneous diagonalization of operators.

- Added eigenstates method returning eigenstates and eigenvalues to Qobj class.

- Added fileio for saving and loading data sets and/or Qobj's.

- Added hinton function for visualizing density matrices.

Bug Fixes
---------

- Switched Examples to new Signals method used in PySide 1.0.6+.

- Switched ProgressBar to new Signals method.

- Fixed memory issue in expm functions.

- Fixed memory bug in isherm.

- Made all Qobj data complex by default.

- Reduced ODE tolerance levels in Odeoptions.

- Fixed bug in ptrace where dense matrix was used instead of sparse.

- Fixed issue where PyQt4 version would not be displayed in about box.

- Fixed issue in Wigner where xvec was used twice (in place of yvec).


Version 1.0.0 (July 29, 2011)
+++++++++++++++++++++++++++++

- **Initial release.**
.. _classes:

***************
Classes
***************

.. _classes-qobj:

Qobj
--------------

.. autoclass:: qutip.Qobj
    :members:

.. _classes-qobjevo:

QobjEvo
--------------

.. autoclass:: qutip.QobjEvo
    :members:

.. _classes-eseries:

eseries
-----------------

.. autoclass:: qutip.eseries
    :members:

.. _classes-bloch:

Bloch sphere
---------------

.. autoclass:: qutip.bloch.Bloch
    :members:


Distributions
-------------

.. autoclass:: qutip.QFunc
    :members:


Cubic Spline
---------------

.. autoclass:: qutip.interpolate.Cubic_Spline
    :members:


.. _classes-non_markov:

Non-Markovian Solvers
---------------------

.. autoclass:: qutip.nonmarkov.heom.HEOMSolver
    :members:

.. autoclass:: qutip.nonmarkov.heom.HSolverDL
    :members:

.. autoclass:: qutip.nonmarkov.heom.BathExponent
    :members:

.. autoclass:: qutip.nonmarkov.heom.Bath
    :members:

.. autoclass:: qutip.nonmarkov.heom.BosonicBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.DrudeLorentzBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.DrudeLorentzPadeBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.UnderDampedBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.FermionicBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.LorentzianBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.LorentzianPadeBath
    :members:

.. autoclass:: qutip.nonmarkov.heom.HierarchyADOs
    :members:

.. autoclass:: qutip.nonmarkov.heom.HierarchyADOsState
    :members:

.. autoclass:: qutip.nonmarkov.dlheom_solver.HSolverDL
    :members:

.. autoclass:: qutip.nonmarkov.dlheom_solver.HEOMSolver
    :members:

.. autoclass:: qutip.nonmarkov.memorycascade.MemoryCascade
    :members:

.. autoclass:: qutip.nonmarkov.transfertensor.TTMSolverOptions
    :members:


.. _classes-odeoptions:

Solver Options and Results
---------------------------

.. autoclass:: qutip.solver.ExpectOps
    :members:

.. autoclass:: qutip.solver.Options
    :members:

.. autoclass:: qutip.solver.Result
    :members:

.. autoclass:: qutip.solver.SolverConfiguration
    :members:

.. autoclass:: qutip.solver.Stats
    :members:

.. autoclass:: qutip.stochastic.StochasticSolverOptions
    :members:

.. _classes-piqs:

Permutational Invariance
------------------------

.. autoclass:: qutip.piqs.Dicke
    :members:

.. autoclass:: qutip.piqs.Pim
    :members:

.. _classes-distributions:

One-Dimensional Lattice
-----------------------

.. autoclass:: qutip.lattice.Lattice1d
    :members:

Distribution functions
----------------------------

.. autoclass:: qutip.distributions.Distribution
    :members:

.. autoclass:: qutip.distributions.WignerDistribution
    :members:

.. autoclass:: qutip.distributions.QDistribution
    :members:

.. autoclass:: qutip.distributions.TwoModeQuadratureCorrelation
    :members:

.. autoclass:: qutip.distributions.HarmonicOscillatorWaveFunction
    :members:

.. autoclass:: qutip.distributions.HarmonicOscillatorProbabilityFunction
    :members:

.. _classes-qip:

Quantum information processing
------------------------------

.. autoclass:: qutip.qip.Gate
    :members:

.. autoclass:: qutip.qip.circuit.Measurement
    :members:

.. autoclass:: qutip.qip.circuit.QubitCircuit
    :members:

.. autoclass:: qutip.qip.circuit.CircuitResult
    :members:

.. autoclass:: qutip.qip.circuit.CircuitSimulator
    :members:

.. autoclass:: qutip.qip.device.Processor
    :members:

.. autoclass:: qutip.qip.device.OptPulseProcessor
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.device.ModelProcessor
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.device.SpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.device.LinearSpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.device.CircularSpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.device.DispersiveCavityQED
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.noise.Noise
    :members:

.. autoclass:: qutip.qip.noise.DecoherenceNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.noise.RelaxationNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.noise.ControlAmpNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.noise.RandomNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.pulse.Pulse
    :members:

.. autoclass:: qutip.qip.compiler.GateCompiler
    :members:

.. autoclass:: qutip.qip.compiler.CavityQEDCompiler
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.compiler.SpinChainCompiler
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.compiler.Scheduler
    :members:

.. autoclass:: qutip.qip.compiler.Instruction
    :members:

.. _classes-control:

Optimal control
---------------

.. autoclass:: qutip.control.optimizer.Optimizer
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerBFGS
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerLBFGSB
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerCrab
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerCrabFmin
    :members:

.. autoclass:: qutip.control.optimizer.OptimIterSummary
    :members:

.. autoclass:: qutip.control.termcond.TerminationConditions
    :members:

.. autoclass:: qutip.control.optimresult.OptimResult
    :members:

.. autoclass:: qutip.control.dynamics.Dynamics
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsGenMat
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsUnitary
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsSymplectic
    :members:

.. autoclass:: qutip.control.propcomp.PropagatorComputer
    :members:

.. autoclass:: qutip.control.propcomp.PropCompApproxGrad
    :members:

.. autoclass:: qutip.control.propcomp.PropCompDiag
    :members:

.. autoclass:: qutip.control.propcomp.PropCompFrechet
    :members:

.. autoclass:: qutip.control.fidcomp.FidelityComputer
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompUnitary
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompTraceDiff
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompTraceDiffApprox
    :members:

.. autoclass:: qutip.control.tslotcomp.TimeslotComputer
    :members:

.. autoclass:: qutip.control.tslotcomp.TSlotCompUpdateAll
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGen
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenRandom
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenZero
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenLinear
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenPeriodic
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSine
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSquare
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSaw
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenTriangle
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenGaussian
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenGaussianEdge
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenCrab
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenCrabFourier
    :members:

.. autoclass:: qutip.control.stats.Stats
    :members:

.. autoclass:: qutip.control.dump.Dump
    :members:

.. autoclass:: qutip.control.dump.OptimDump
    :members:

.. autoclass:: qutip.control.dump.DynamicsDump
    :members:

.. autoclass:: qutip.control.dump.DumpItem
    :members:

.. autoclass:: qutip.control.dump.EvoCompDumpItem
    :members:

.. autoclass:: qutip.control.dump.DumpSummaryItem
    :members:
.. _apidoc:

*****************
API documentation
*****************

This chapter contains automatically generated API documentation, including a
complete list of QuTiP's public classes and functions.

.. toctree::
   :maxdepth: 3

   classes.rst 
   functions.rst
.. _functions:

***************
Functions
***************

Manipulation and Creation of States and Operators
=================================================

Quantum States
--------------

.. automodule:: qutip.states
    :members: basis, bell_state, bra, coherent, coherent_dm, enr_state_dictionaries, enr_thermal_dm, enr_fock, fock, fock_dm, ghz_state, maximally_mixed_dm, ket, ket2dm, phase_basis, projection, qutrit_basis, singlet_state, spin_state, spin_coherent, state_number_enumerate, state_number_index, state_index_number, state_number_qobj, thermal_dm, triplet_states, w_state, zero_ket


Quantum Operators
-----------------

.. automodule:: qutip.operators
    :members: charge, commutator, create, destroy, displace, enr_destroy, enr_identity, jmat, num, qeye, identity, momentum, phase, position, qdiags, qutrit_ops, qzero, sigmam, sigmap, sigmax, sigmay, sigmaz, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, squeeze, squeezing, tunneling


.. _functions-rand:

Quantum Objects
---------------

.. automodule:: qutip.qobj
    :members: qobj_list_evaluate, ptrace, dag, isequal, issuper, isoper, isoperket, isoperbra, isket, isbra, isherm, shape, dims


Random Operators and States
---------------------------

.. automodule:: qutip.random_objects
    :members: rand_dm, rand_dm_ginibre, rand_dm_hs, rand_herm, rand_ket, rand_ket_haar, rand_stochastic, rand_unitary, rand_unitary_haar, rand_super, rand_super_bcsz


Three-Level Atoms
-----------------

.. automodule:: qutip.three_level_atom
    :members: three_level_basis, three_level_ops
    :undoc-members:


Superoperators and Liouvillians
-------------------------------

.. automodule:: qutip.superoperator
    :members: operator_to_vector, vector_to_operator, liouvillian, spost, spre, sprepost, lindblad_dissipator

Superoperator Representations
-----------------------------

.. automodule:: qutip.superop_reps
    :members: super_to_choi, choi_to_super, choi_to_kraus, kraus_to_choi, kraus_to_super, choi_to_chi, chi_to_choi, to_choi, to_chi, to_super, to_kraus, to_stinespring
    :undoc-members:

Operators and Superoperator Dimensions
--------------------------------------

.. automodule:: qutip.dimensions
    :members: is_scalar, is_vector, is_vectorized_oper, type_from_dims, flatten, deep_remove, unflatten, collapse_dims_oper, collapse_dims_super, enumerate_flat, deep_map, dims_to_tensor_perm, dims_to_tensor_shape, dims_idxs_to_tensor_idxs


Functions acting on states and operators
========================================

Expectation Values
------------------

.. automodule:: qutip.expect
    :members: expect, variance


Tensor
------

.. automodule:: qutip.tensor
    :members: tensor, super_tensor, composite, tensor_contract



Partial Transpose
-----------------

.. automodule:: qutip.partial_transpose
    :members: partial_transpose


.. _functions-entropy:

Entropy Functions
-----------------

.. automodule:: qutip.entropy
    :members: concurrence, entropy_conditional, entropy_linear, entropy_mutual, entropy_relative, entropy_vn


Density Matrix Metrics
----------------------

.. automodule:: qutip.metrics
    :members: fidelity, tracedist, bures_dist, bures_angle, hilbert_dist, average_gate_fidelity, process_fidelity


Continuous Variables
--------------------

.. automodule:: qutip.continuous_variables
    :members: correlation_matrix, covariance_matrix, correlation_matrix_field, correlation_matrix_quadrature, wigner_covariance_matrix, logarithmic_negativity


Measurement
===========

Measurement of quantum states
-----------------------------

.. automodule:: qutip.measurement
    :members: measure, measure_observable, measurement_statistics, measurement_statistics_observable


Dynamics and Time-Evolution
===========================

Schrödinger Equation
--------------------

.. automodule:: qutip.sesolve
    :members: sesolve

Master Equation
---------------

.. automodule:: qutip.mesolve
    :members: mesolve

Monte Carlo Evolution
---------------------

.. automodule:: qutip.mcsolve
    :members: mcsolve

.. ignore f90 stuff for now
    .. automodule:: qutip.fortran.mcsolve_f90
        :members: mcsolve_f90


Exponential Series
------------------

.. automodule:: qutip.essolve
    :members: essolve, ode2es


Bloch-Redfield Master Equation
------------------------------

.. automodule:: qutip.bloch_redfield
    :members: brmesolve, bloch_redfield_tensor, bloch_redfield_solve


Floquet States and Floquet-Markov Master Equation
-------------------------------------------------

.. automodule:: qutip.floquet
    :members: fmmesolve, floquet_modes, floquet_modes_t, floquet_modes_table, floquet_modes_t_lookup, floquet_states, floquet_states_t, floquet_wavefunction, floquet_wavefunction_t, floquet_state_decomposition, fsesolve, floquet_master_equation_rates, floquet_master_equation_steadystate, floquet_basis_transform, floquet_markov_mesolve


Stochastic Schrödinger Equation and Master Equation
---------------------------------------------------

.. automodule:: qutip.stochastic
    :members: ssesolve, photocurrent_sesolve, smepdpsolve, smesolve, photocurrent_mesolve, ssepdpsolve, stochastic_solvers, general_stochastic


Correlation Functions
---------------------

.. automodule:: qutip.correlation
    :members: correlation, correlation_ss, correlation_2op_1t, correlation_2op_2t, correlation_3op_1t, correlation_3op_2t, correlation_4op_1t, correlation_4op_2t, spectrum, spectrum_ss, spectrum_pi, spectrum_correlation_fft, coherence_function_g1, coherence_function_g2


Steady-state Solvers
--------------------

.. automodule:: qutip.steadystate
    :members: steadystate, build_preconditioner
    :undoc-members:

Propagators
-----------

.. automodule:: qutip.propagator
    :members: propagator, propagator_steadystate
    :undoc-members:


Time-dependent problems
-----------------------

.. automodule:: qutip.rhs_generate
    :members: rhs_generate, rhs_clear

Scattering in Quantum Optical Systems
-------------------------------------

.. automodule:: qutip.scattering
    :members: temporal_basis_vector, temporal_scattered_state, scattering_probability
    :undoc-members:

Permutational Invariance
------------------------

.. automodule:: qutip.piqs
    :members: num_dicke_states, num_dicke_ladders, num_tls, isdiagonal, dicke_blocks, dicke_blocks_full, dicke_function_trace, purity_dicke, entropy_vn_dicke, state_degeneracy, m_degeneracy, energy_degeneracy, ap, am, spin_algebra, jspin, collapse_uncoupled, dicke_basis, dicke, excited, superradiant, css, ghz, ground, identity_uncoupled, block_matrix, tau_column,


Lattice
=======

Lattice Properties
------------------
.. automodule:: qutip.lattice
    :members: cell_structures


Topology
--------
.. automodule:: qutip.topology
    :members: berry_curvature, plot_berry_curvature


Visualization
===============

Pseudoprobability Functions
---------------------------

.. automodule:: qutip.wigner
    :members: qfunc, spin_q_function, spin_wigner, wigner


Graphs and Visualization
------------------------

.. automodule:: qutip.visualization
    :members: hinton, matrix_histogram, matrix_histogram_complex, plot_energy_levels, plot_fock_distribution, plot_wigner_fock_distribution, plot_wigner, sphereplot, plot_schmidt, plot_qubism, plot_expectation_values, plot_spin_distribution_2d, plot_spin_distribution_3d, plot_wigner_sphere
    :undoc-members:

.. automodule:: qutip.orbital
    :members: orbital

.. automodule:: qutip.matplotlib_utilities
   :members: wigner_cmap, complex_phase_cmap


Quantum Process Tomography
--------------------------

.. automodule:: qutip.tomography
    :members: qpt, qpt_plot, qpt_plot_combined
    :undoc-members:



.. _functions-qip:

Quantum Information Processing
==============================

Gates
-----

.. automodule:: qutip.qip.operations.gates
    :members: rx, ry, rz, sqrtnot, snot, phasegate, cphase, cnot, csign, berkeley, swapalpha, swap, iswap, sqrtswap, sqrtiswap, fredkin, toffoli, rotation, controlled_gate, globalphase, hadamard_transform, gate_sequence_product, gate_expand_1toN, gate_expand_2toN, gate_expand_3toN, expand_operator

Qubits
------

.. automodule:: qutip.qip.qubits
    :members: qubit_states

Algorithms
----------

.. automodule:: qutip.qip.algorithms.qft
    :members: qft, qft_steps, qft_gate_sequence


Circuit
-------

.. automodule:: qutip.qip.qasm
    :members: read_qasm, save_qasm, print_qasm, circuit_to_qasm_str

.. _functions-non_markov:

Non-Markovian Solvers
=====================

.. automodule:: qutip.nonmarkov.transfertensor
    :members: ttmsolve

.. _functions-control:

Optimal control
===============

.. automodule:: qutip.control.pulseoptim
    :members: optimize_pulse, optimize_pulse_unitary, create_pulse_optimizer, opt_pulse_crab, opt_pulse_crab_unitary

.. automodule:: qutip.control.pulsegen
    :members: create_pulse_gen

Utility Functions
=================

.. _functions-graph:

Graph Theory Routines
---------------------

.. automodule:: qutip.graph
    :members: breadth_first_search, graph_degree, reverse_cuthill_mckee, maximum_bipartite_matching, weighted_bipartite_matching


.. _functions-utilities:

Utility Functions
-----------------

.. automodule:: qutip.utilities
    :members: n_thermal, clebsch, convert_unit


.. _functions-fileio:

File I/O Functions
------------------

.. automodule:: qutip.fileio
    :members: file_data_read, file_data_store, qload, qsave


.. _functions-parallel:

Parallelization
---------------

.. automodule:: qutip.parallel
    :members: parfor, parallel_map, serial_map


.. _functions-ipython:

Semidefinite Programming
------------------------

.. automodule:: qutip.semidefinite
    :members: complex_var, herm, pos_noherm, pos, dens, kron, conj, bmat, bmat, memoize, qudit_swap, dnorm_problem


.. _functions-semidefinite:

IPython Notebook Tools
----------------------

.. automodule:: qutip.ipynbtools
    :members: parfor, parallel_map, version_table

.. _functions-misc:

Miscellaneous
-------------

.. automodule:: qutip
    :members: about, simdiag
.. _visual:


.. plot::
   :include-source: False

   import numpy as np

   from qutip import *

   import pylab as plt

   from warnings import warn

   plt.close("all")


*********************************************
Visualization of quantum states and processes
*********************************************

Visualization is often an important complement to a simulation of a quantum
mechanical system. The first method of visualization that come to mind might be
to plot the expectation values of a few selected operators. But on top of that,
it can often be instructive to visualize for example the state vectors or
density matices that describe the state of the system, or how the state is
transformed as a function of time (see process tomography below). In this
section we demonstrate how QuTiP and matplotlib can be used to perform a few
types of  visualizations that often can provide additional understanding of
quantum system.

.. _visual-fock:

Fock-basis probability distribution
===================================

In quantum mechanics probability distributions plays an important role, and as
in statistics, the expectation values computed from a probability distribution
does not reveal the full story. For example, consider an quantum harmonic
oscillator mode with Hamiltonian :math:`H = \hbar\omega a^\dagger a`, which is
in a state described by its density matrix :math:`\rho`, and which on average
is occupied by two photons, :math:`\mathrm{Tr}[\rho a^\dagger a] = 2`. Given
this information we cannot say whether the oscillator is in a Fock state,
a thermal state, a coherent state, etc. By visualizing the photon distribution
in the Fock state basis important clues about the underlying state can be
obtained.

One convenient way to visualize a probability distribution is to use histograms.
Consider the following histogram visualization of the number-basis probability
distribution, which can be obtained from the diagonal of the density matrix,
for a few possible oscillator states with on average occupation of two photons.

First we generate the density matrices for the coherent, thermal and fock states.

.. plot::
    :context: reset

    N = 20

    rho_coherent = coherent_dm(N, np.sqrt(2))

    rho_thermal = thermal_dm(N, 2)

    rho_fock = fock_dm(N, 2)


Next, we plot histograms of the diagonals of the density matrices:

.. plot::
    :context:

    fig, axes = plt.subplots(1, 3, figsize=(12,3))

    bar0 = axes[0].bar(np.arange(0, N)-.5, rho_coherent.diag())

    lbl0 = axes[0].set_title("Coherent state")

    lim0 = axes[0].set_xlim([-.5, N])

    bar1 = axes[1].bar(np.arange(0, N)-.5, rho_thermal.diag())

    lbl1 = axes[1].set_title("Thermal state")

    lim1 = axes[1].set_xlim([-.5, N])

    bar2 = axes[2].bar(np.arange(0, N)-.5, rho_fock.diag())

    lbl2 = axes[2].set_title("Fock state")

    lim2 = axes[2].set_xlim([-.5, N])

    plt.show()


All these states correspond to an average of two photons, but by visualizing
the photon distribution in Fock basis the differences between these states are
easily appreciated.

One frequently need to visualize the Fock-distribution in the way described
above, so QuTiP provides a convenience function for doing this, see
:func:`qutip.visualization.plot_fock_distribution`, and the following example:

.. plot::
    :context: close-figs

    fig, axes = plt.subplots(1, 3, figsize=(12,3))

    plot_fock_distribution(rho_coherent, fig=fig, ax=axes[0], title="Coherent state");

    plot_fock_distribution(rho_thermal, fig=fig, ax=axes[1], title="Thermal state");

    plot_fock_distribution(rho_fock, fig=fig, ax=axes[2], title="Fock state");

    fig.tight_layout()

    plt.show()

.. _visual-dist:

Quasi-probability distributions
===============================

The probability distribution in the number (Fock) basis only describes the
occupation probabilities for a discrete set of states. A more complete
phase-space probability-distribution-like function for harmonic modes are
the Wigner and Husumi Q-functions, which are full descriptions of the
quantum state (equivalent to the density matrix). These are called
quasi-distribution functions because unlike real probability distribution
functions they can for example be negative. In addition to being more complete descriptions
of a state (compared to only the occupation probabilities plotted above),
these distributions are also great for demonstrating if a quantum state is
quantum mechanical, since for example a negative Wigner function
is a definite indicator that a state is distinctly nonclassical.


Wigner function
---------------

In QuTiP, the Wigner function for a harmonic mode can be calculated with the
function :func:`qutip.wigner.wigner`. It takes a ket or a density matrix as
input, together with arrays that define the ranges of the phase-space
coordinates (in the x-y plane). In the following example the Wigner functions
are calculated and plotted for the same three states as in the previous section.

.. plot::
    :context: close-figs

    xvec = np.linspace(-5,5,200)

    W_coherent = wigner(rho_coherent, xvec, xvec)

    W_thermal = wigner(rho_thermal, xvec, xvec)

    W_fock = wigner(rho_fock, xvec, xvec)

    # plot the results

    fig, axes = plt.subplots(1, 3, figsize=(12,3))

    cont0 = axes[0].contourf(xvec, xvec, W_coherent, 100)

    lbl0 = axes[0].set_title("Coherent state")

    cont1 = axes[1].contourf(xvec, xvec, W_thermal, 100)

    lbl1 = axes[1].set_title("Thermal state")

    cont0 = axes[2].contourf(xvec, xvec, W_fock, 100)

    lbl2 = axes[2].set_title("Fock state")

    plt.show()

.. _visual-cmap:

Custom Color Maps
~~~~~~~~~~~~~~~~~

The main objective when plotting a Wigner function is to demonstrate that the underlying
state is nonclassical, as indicated by negative values in the Wigner function.  Therefore,
making these negative values stand out in a figure is helpful for both analysis and publication
purposes.  Unfortunately, all of the color schemes used in Matplotlib (or any other plotting software)
are linear colormaps where small negative values tend to be near the same color as the zero values, and
are thus hidden.  To fix this dilemma, QuTiP includes a nonlinear colormap function :func:`qutip.matplotlib_utilities.wigner_cmap`
that colors all negative values differently than positive or zero values.  Below is a demonstration of how to use
this function in your Wigner figures:

.. plot::
    :context: close-figs

    import matplotlib as mpl

    from matplotlib import cm

    psi = (basis(10, 0) + basis(10, 3) + basis(10, 9)).unit()

    xvec = np.linspace(-5, 5, 500)

    W = wigner(psi, xvec, xvec)

    wmap = wigner_cmap(W)  # Generate Wigner colormap

    nrm = mpl.colors.Normalize(-W.max(), W.max())

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    plt1 = axes[0].contourf(xvec, xvec, W, 100, cmap=cm.RdBu, norm=nrm)

    axes[0].set_title("Standard Colormap");

    cb1 = fig.colorbar(plt1, ax=axes[0])

    plt2 = axes[1].contourf(xvec, xvec, W, 100, cmap=wmap)  # Apply Wigner colormap

    axes[1].set_title("Wigner Colormap");

    cb2 = fig.colorbar(plt2, ax=axes[1])

    fig.tight_layout()

    plt.show()



Husimi Q-function
-----------------

The Husimi Q function is, like the Wigner function, a quasiprobability
distribution for harmonic modes. It is defined as

.. math::

    Q(\alpha) = \frac{1}{\pi}\left<\alpha|\rho|\alpha\right>

where :math:`\left|\alpha\right>` is a coherent state and
:math:`\alpha = x + iy`. In QuTiP, the Husimi Q function can be computed given
a state ket or density matrix using the function :func:`.qfunc`, as
demonstrated below.

.. plot::
    :context: close-figs

    Q_coherent = qfunc(rho_coherent, xvec, xvec)
    Q_thermal = qfunc(rho_thermal, xvec, xvec)
    Q_fock = qfunc(rho_fock, xvec, xvec)
    fig, axes = plt.subplots(1, 3, figsize=(12,3))
    cont0 = axes[0].contourf(xvec, xvec, Q_coherent, 100)
    lbl0 = axes[0].set_title("Coherent state")
    cont1 = axes[1].contourf(xvec, xvec, Q_thermal, 100)
    lbl1 = axes[1].set_title("Thermal state")
    cont0 = axes[2].contourf(xvec, xvec, Q_fock, 100)
    lbl2 = axes[2].set_title("Fock state")
    plt.show()

If you need to calculate the Q function for many states with the same
phase-space coordinates, it is more efficient to use the :obj:`.QFunc` class.
This stores various intermediary results to achieve an order-of-magnitude
improvement compared to calling :obj:`.qfunc` in a loop.

.. code-block:: python

   xs = np.linspace(-1, 1, 101)
   qfunc_calculator = qutip.QFunc(xs, xs)
   q_state1 = qfunc_calculator(qutip.rand_dm(5))
   q_state2 = qfunc_calculator(qutip.rand_ket(100))


.. _visual-oper:

Visualizing operators
=====================

Sometimes, it may also be useful to directly visualizing the underlying matrix
representation of an operator. The density matrix, for example, is an operator
whose elements can give insights about the state it represents, but one might
also be interesting in plotting the matrix of an Hamiltonian to inspect the
structure and relative importance of various elements.

QuTiP offers a few functions for quickly visualizing matrix data in the
form of histograms, :func:`qutip.visualization.matrix_histogram` and
:func:`qutip.visualization.matrix_histogram_complex`, and as Hinton diagram of weighted
squares, :func:`qutip.visualization.hinton`. These functions takes a
:class:`qutip.Qobj.Qobj` as first argument, and optional arguments to, for
example, set the axis labels and figure title (see the function's documentation
for details).

For example, to illustrate the use of :func:`qutip.visualization.matrix_histogram`,
let's visualize of the Jaynes-Cummings Hamiltonian:

.. plot::
    :context: close-figs

    N = 5

    a = tensor(destroy(N), qeye(2))

    b = tensor(qeye(N), destroy(2))

    sx = tensor(qeye(N), sigmax())

    H = a.dag() * a + sx - 0.5 * (a * b.dag() + a.dag() * b)

    # visualize H

    lbls_list = [[str(d) for d in range(N)], ["u", "d"]]

    xlabels = []

    for inds in tomography._index_permutations([len(lbls) for lbls in lbls_list]):
       xlabels.append("".join([lbls_list[k][inds[k]] for k in range(len(lbls_list))]))

    fig, ax = matrix_histogram(H, xlabels, xlabels, limits=[-4,4])

    ax.view_init(azim=-55, elev=45)

    plt.show()


Similarly, we can use the function :func:`qutip.visualization.hinton`, which is
used below to visualize the corresponding steadystate density matrix:

.. plot::
    :context: close-figs

    rho_ss = steadystate(H, [np.sqrt(0.1) * a, np.sqrt(0.4) * b.dag()])

    hinton(rho_ss)

    plt.show()

.. _visual-qpt:

Quantum process tomography
==========================

Quantum process tomography (QPT) is a useful technique for characterizing experimental implementations of quantum gates involving a small number of qubits. It can also be a useful theoretical tool that can give insight in how a process transforms states, and it can be used for example to study how noise or other imperfections deteriorate a gate. Whereas a fidelity or distance measure can give a single number that indicates how far from ideal a gate is, a quantum process tomography analysis can give detailed information about exactly what kind of errors various imperfections introduce.

The idea is to construct a transformation matrix for a quantum process (for example a quantum gate) that describes how the density matrix of a system is transformed by the process. We can then decompose the transformation in some operator basis that represent well-defined and easily interpreted transformations of the input states.

To see how this works (see e.g. [Moh08]_ for more details), consider a process that is described by quantum map :math:`\epsilon(\rho_{\rm in}) = \rho_{\rm out}`, which can be written

.. math::
    :label: qpt-quantum-map

    \epsilon(\rho_{\rm in}) = \rho_{\rm out} = \sum_{i}^{N^2} A_i \rho_{\rm in} A_i^\dagger,

where :math:`N` is the number of states of the system (that is, :math:`\rho` is represented by an :math:`[N\times N]` matrix). Given an orthogonal operator basis of our choice :math:`\{B_i\}_i^{N^2}`, which satisfies :math:`{\rm Tr}[B_i^\dagger B_j] = N\delta_{ij}`, we can write the map as

.. math::
    :label: qpt-quantum-map-transformed

    \epsilon(\rho_{\rm in}) = \rho_{\rm out} = \sum_{mn} \chi_{mn} B_m \rho_{\rm in} B_n^\dagger.

where :math:`\chi_{mn} = \sum_{ij} b_{im}b_{jn}^*` and :math:`A_i = \sum_{m} b_{im}B_{m}`. Here, matrix :math:`\chi` is the transformation matrix we are after, since it describes how much :math:`B_m \rho_{\rm in} B_n^\dagger` contributes to :math:`\rho_{\rm out}`.

In a numerical simulation of a quantum process we usually do not have access to the quantum map in the form Eq. :eq:`qpt-quantum-map`. Instead, what we usually can do is to calculate the propagator :math:`U` for the density matrix in superoperator form, using for example the QuTiP function :func:`qutip.propagator.propagator`. We can then write

.. math::

    \epsilon(\tilde{\rho}_{\rm in}) = U \tilde{\rho}_{\rm in} = \tilde{\rho}_{\rm out}

where :math:`\tilde{\rho}` is the vector representation of the density matrix :math:`\rho`. If we write Eq. :eq:`qpt-quantum-map-transformed` in superoperator form as well we obtain

.. math::

    \tilde{\rho}_{\rm out} = \sum_{mn} \chi_{mn} \tilde{B}_m \tilde{B}_n^\dagger \tilde{\rho}_{\rm in} = U \tilde{\rho}_{\rm in}.

so we can identify

.. math::

    U = \sum_{mn} \chi_{mn} \tilde{B}_m \tilde{B}_n^\dagger.

Now this is a linear equation systems for the :math:`N^2 \times N^2` elements in :math:`\chi`. We can solve it by writing :math:`\chi` and the superoperator propagator as :math:`[N^4]` vectors, and likewise write the superoperator product :math:`\tilde{B}_m\tilde{B}_n^\dagger` as a :math:`[N^4\times N^4]` matrix :math:`M`:

.. math::

    U_I = \sum_{J}^{N^4} M_{IJ} \chi_{J}

with the solution

.. math::

    \chi = M^{-1}U.

Note that to obtain :math:`\chi` with this method we have to construct a matrix :math:`M` with a size that is the square of the size of the superoperator for the system. Obviously, this scales very badly with increasing system size, but this method can still be a very useful for small systems (such as system comprised of a small number of coupled qubits).

Implementation in QuTiP
-----------------------

In QuTiP, the procedure described above is implemented in the function :func:`qutip.tomography.qpt`, which returns the :math:`\chi` matrix given a density matrix propagator. To illustrate how to use this function, let's consider the :math:`i`-SWAP gate for two qubits. In QuTiP the function :func:`qutip.qip.operations.iswap` generates the unitary transformation for the state kets:


.. plot::
    :context: close-figs

    from qutip.qip.operations import iswap

    U_psi = iswap()

To be able to use this unitary transformation matrix as input to the function :func:`qutip.tomography.qpt`, we first need to convert it to a transformation matrix for the corresponding density matrix:

.. plot::
    :context:

    U_rho = spre(U_psi) * spost(U_psi.dag())


Next, we construct a list of operators that define the basis :math:`\{B_i\}` in the form of a list of operators for each composite system. At the same time, we also construct a list of corresponding labels that will be used when plotting the :math:`\chi` matrix.

.. plot::
    :context:

    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
    op_label = [["i", "x", "y", "z"]] * 2


We are now ready to compute :math:`\chi` using :func:`qutip.tomography.qpt`, and to plot it using :func:`qutip.tomography.qpt_plot_combined`.

.. plot::
    :context:

    chi = qpt(U_rho, op_basis)

    fig = qpt_plot_combined(chi, op_label, r'$i$SWAP')

    plt.show()



For a slightly more advanced example, where the density matrix propagator is calculated from the dynamics of a system defined by its Hamiltonian and collapse operators using the function :func:`qutip.propagator.propagator`, see notebook "Time-dependent master equation: Landau-Zener transitions" on the tutorials section on the QuTiP web site.
.. _control:

*********************************************
Quantum Optimal Control
*********************************************


Introduction
=============

In quantum control we look to prepare some specific state, effect some state-to-state transfer, or effect some transformation (or gate) on a quantum system. For a given quantum system there will always be factors that effect the dynamics that are outside of our control. As examples, the interactions between elements of the system or a magnetic field required to trap the system. However, there may be methods of affecting the dynamics in a controlled way, such as the time varying amplitude of the electric component of an interacting laser field. And so this leads to some questions; given a specific quantum system with known time-independent dynamics generator (referred to as the *drift* dynamics generators) and set of externally controllable fields for which the interaction can be described by *control* dynamics generators:

1. What states or transformations can we achieve (if any)?

2. What is the shape of the control pulse required to achieve this?

These questions are addressed as *controllability* and *quantum optimal control* [dAless08]_. The answer to question of *controllability* is determined by the commutability of the dynamics generators and is formalised as the *Lie Algebra Rank Criterion* and is discussed in detail in [dAless08]_. The solutions to the second question can be determined through optimal control algorithms, or control pulse optimisation.

.. figure:: figures/quant_optim_ctrl.png
   :align: center
   :width: 2.5in

   Schematic showing the principle of quantum control.

Quantum Control has many applications including NMR, *quantum metrology*, *control of chemical reactions*, and *quantum information processing*.

To explain the physics behind these algorithms we will first consider only finite-dimensional, closed quantum systems.

Closed Quantum Systems
======================
In closed quantum systems the states can be represented by kets, and the transformations on these states are unitary operators. The dynamics generators are Hamiltonians. The combined Hamiltonian for the system is given by

.. math::

    H(t) = H_0 + \sum_{j=1} u_j(t) H_j

where :math:`H_0` is the drift Hamiltonian and the :math:`H_j` are the control Hamiltonians. The :math:`u_j` are time varying amplitude functions for the specific control.

The dynamics of the system are governed by *Schrödingers equation*.

.. math::

    \tfrac{d}{dt} \ket{\psi} = -i H(t)\ket{\psi}

Note we use units where :math:`\hbar=1` throughout. The solutions to Schrödinger's equation are of the form:

.. math::

    \ket{\psi(t)} = U(t)\ket{\psi_0}

where :math:`\psi_0` is the state of the system at :math:`t=0` and :math:`U(t)` is a unitary operator on the Hilbert space containing the states. :math:`U(t)` is a solution to the *Schrödinger operator equation*

.. math::

    \tfrac{d}{dt}U = -i H(t)U ,\quad U(0) = \mathbb{1}

We can use optimal control algorithms to determine a set of :math:`u_j` that will drive our system from :math:`\ket{\psi_0}` to :math:`\ket{\psi_1}`, this is state-to-state transfer, or drive the system from some arbitary state to a given state :math:`\ket{\psi_1}`, which is state preparation, or effect some unitary transformation :math:`U_{target}`, called gate synthesis. The latter of these is most important in quantum computation.


The GRAPE algorithm
===================
The **GR**\ adient **A**\ scent **P**\ ulse **E**\ ngineering was first proposed in [2]. Solutions to Schrödinger's equation for a time-dependent Hamiltonian are not generally possible to obtain analytically. Therefore, a piecewise constant approximation to the pulse amplitudes is made. Time allowed for the system to evolve :math:`T` is split into :math:`M` timeslots (typically these are of equal duration), during which the control amplitude is assumed to remain constant. The combined Hamiltonian can then be approximated as:

.. math::

    H(t) \approx H(t_k) = H_0 + \sum_{j=1}^N u_{jk} H_j\quad

where :math:`k` is a timeslot index, :math:`j` is the control index, and :math:`N` is the number of controls. Hence :math:`t_k` is the evolution time at the start of the timeslot, and :math:`u_{jk}` is the amplitude of control :math:`j` throughout timeslot :math:`k`. The time evolution operator, or propagator, within the timeslot can then be calculated as:

.. math::

    X_k:=e^{-iH(t_k)\Delta t_k}

where :math:`\Delta t_k` is the duration of the timeslot. The evolution up to (and including) any timeslot :math:`k` (including the full evolution :math:`k=M`) can the be calculated as

.. math::

    X(t_k):=X_k X_{k-1}\cdots X_1 X_0

If the objective is state-to-state transfer then :math:`X_0=\ket{\psi_0}` and the target :math:`X_{targ}=\ket{\psi_1}`, for gate synthesis :math:`X_0 = U(0) = \mathbb{1}` and the target :math:`X_{targ}=U_{targ}`.

A *figure of merit* or *fidelity* is some measure of how close the evolution is to the target, based on the  control amplitudes in the timeslots. The typical figure of merit for unitary systems is the normalised overlap of the evolution and the target.

.. math::

    f_{PSU} = \tfrac{1}{d} \big| \tr \{X_{targ}^{\dagger} X(T)\} \big|

where :math:`d` is the system dimension. In this figure of merit the absolute value is taken to ignore any differences in global phase, and :math:`0 \le f \le 1`. Typically the fidelity error (or *infidelity*) is more useful, in this case defined as :math:`\varepsilon = 1 - f_{PSU}`.  There are many other possible objectives, and hence figures of merit.

As there are now :math:`N \times M` variables (the :math:`u_{jk}`) and one
parameter to minimise :math:`\varepsilon`, then the problem becomes a finite
multi-variable optimisation problem, for which there are many established
methods, often referred to as 'hill-climbing' methods. The simplest of these to
understand is that of steepest ascent (or descent). The gradient of the
fidelity with respect to all the variables is calculated (or approximated) and
a step is made in the variable space in the direction of steepest ascent (or
descent). This method is a first order gradient method. In two dimensions this
describes a method of climbing a hill by heading in the direction where the
ground rises fastest. This analogy also clearly illustrates one of the main
challenges in multi-variable optimisation, which is that all methods have a
tendency to get stuck in local maxima. It is hard to determine whether one has
found a global maximum or not - a local peak is likely not to be the highest
mountain in the region. In quantum optimal control we can typically define an
infidelity that has a lower bound of zero. We can then look to minimise the
infidelity (from here on we will only consider optimising for infidelity
minima). This means that we can terminate any pulse optimisation when the
infidelity reaches zero (to a sufficient precision). This is however only
possible for fully controllable systems; otherwise it is hard (if not
impossible) to know that the minimum possible infidelity has been achieved. In
the hill walking analogy the step size is roughly fixed to a stride, however,
in computations the step size must be chosen. Clearly there is a trade-off here
between the number of steps (or iterations) required to reach the minima and
the possibility that we might step over a minima. In practice it is difficult
to determine an efficient and effective step size.

The second order differentials of the infidelity with respect to the variables
can be used to approximate the local landscape to a parabola. This way a step
(or jump) can be made to where the minima would be if it were parabolic. This
typically vastly reduces the number of iterations, and removes the need to
guess a step size. The method where all the second differentials are calculated
explicitly is called the *Newton-Raphson* method. However, calculating the
second-order differentials (the Hessian matrix) can be computationally
expensive, and so there are a class of methods known as *quasi-Newton* that
approximate the Hessian based on successive iterations. The most popular of
these (in quantum optimal control) is the Broyden–Fletcher–Goldfarb–Shanno
algorithm (BFGS). The default method in the QuTiP Qtrl GRAPE implementation is
the L-BFGS-B method in Scipy, which is a wrapper to the implementation
described in [Byrd95]_. This limited memory and bounded method does not need to
store the entire Hessian, which reduces the computer memory required, and
allows bounds to be set for variable values, which considering these are field
amplitudes is often physical.

The pulse optimisation is typically far more efficient if the gradients can be
calculated exactly, rather than approximated. For simple fidelity measures such
as :math:`f_{PSU}` this is possible. Firstly the propagator gradient for each
timeslot with respect to the control amplitudes is calculated. For closed
systems, with unitary dynamics, a method using the eigendecomposition is used,
which is efficient as it is also used in the propagator calculation (to
exponentiate the combined Hamiltonian). More generally (for example open
systems and symplectic dynamics) the Frechet derivative (or augmented matrix)
method is used, which is described in [Flo12]_. For other optimisation goals it
may not be possible to calculate analytic gradients. In these cases it is
necessary to approximate the gradients, but this can be very expensive, and can
lead to other algorithms out-performing GRAPE.


The CRAB Algorithm
===================
It has been shown [Lloyd14]_, the dimension of a quantum optimal control
problem is a polynomial function of the dimension of the manifold of the
time-polynomial reachable states, when allowing for a finite control precision
and evolution time. You can think of this as the information content of the
pulse (as being the only effective input) being very limited e.g. the pulse is
compressible to a few bytes without loosing the target.

This is where the **C**\ hopped **RA**\ ndom **B**\ asis (CRAB) algorithm
[Doria11]_, [Caneva11]_ comes into play: Since the pulse complexity is usually
very low, it is sufficient to transform the optimal control problem to a few
parameter search by introducing a physically motivated function basis that
builds up the pulse. Compared to the number of time slices needed to accurately
simulate quantum dynamics (often equals basis dimension for Gradient based
algorithms), this number is lower by orders of magnitude, allowing CRAB to
efficiently optimize smooth pulses with realistic experimental constraints. It
is important to point out, that CRAB does not make any suggestion on the basis
function to be used. The basis must be chosen carefully considered, taking into
account a priori knowledge of the system (such as symmetries, magnitudes of
scales,...) and solution (e.g. sign, smoothness, bang-bang behavior,
singularities, maximum excursion or rate of change,....). By doing so, this
algorithm allows for native integration of experimental constraints such as
maximum frequencies allowed, maximum amplitude, smooth ramping up and down of
the pulse and many more. Moreover initial guesses, if they are available, can
(however not have to) be included to speed up convergence.

As mentioned in the GRAPE paragraph, for CRAB local minima arising from
algorithmic design can occur, too. However, for CRAB a 'dressed' version has
recently been introduced [Rach15]_ that allows to escape local minima.

For some control objectives and/or dynamical quantum descriptions, it is either
not possible to derive the gradient for the cost functional with respect to
each time slice or it is computationally expensive to do so. The same can apply
for the necessary (reverse) propagation of the co-state. All this trouble does
not occur within CRAB as those elements are not in use here. CRAB, instead,
takes the time evolution as a black-box where the pulse goes as an input and
the cost (e.g. infidelity) value will be returned as an output. This concept,
on top, allows for direct integration in a closed loop experimental environment
where both the preliminarily open loop optimization, as well as the final
adoption, and integration to the lab (to account for modeling errors,
experimental systematic noise, ...) can be done all in one, using this
algorithm.

Optimal Quantum Control in QuTiP
================================

There are two separate implementations of optimal control inside QuTiP. The
first is an implementation of first order GRAPE, and is not further described
here, but there are the example notebooks. The second is referred to as Qtrl
(when a distinction needs to be made) as this was its name before it was
integrated into QuTiP. Qtrl uses the Scipy optimize functions to perform the
multi-variable optimisation, typically the L-BFGS-B method for GRAPE and
Nelder-Mead for CRAB. The GRAPE implementation in Qtrl was initially based on
the open-source package  DYNAMO, which is a MATLAB implementation, and is
described in [DYNAMO]_. It has since been restructured and extended for
flexibility and compatibility within QuTiP.

The rest of this section describes the Qtrl implementation and how to use it.

Object Model
  The Qtrl code is organised in a hierarchical object model in order to try and maximise configurability whilst maintaining some clarity. It is not necessary to understand the model in order to use the pulse optimisation functions, but it is the most flexible method of using Qtrl. If you just want to use a simple single function call interface, then jump to :ref:`pulseoptim-functions`

.. figure:: figures/qtrl-code_object_model.png
   :align: center
   :width: 3.5in

   Qtrl code object model.

The object's properties and methods are described in detail in the documentation, so that will not be repeated here.

OptimConfig
  The OptimConfig object is used simply to hold configuration parameters used by all the objects. Typically this is the subclass types for the other objects and parameters for the users specific requirements. The ``loadparams`` module can be used read parameter values from a configuration file.

Optimizer
  This acts as a wrapper to the ``Scipy.optimize`` functions that perform the work of the pulse optimisation algorithms. Using the main classes the user can specify which of the optimisation methods are to be used. There are subclasses specifically for the BFGS and L-BFGS-B methods. There is another subclass for using the CRAB algorithm.

Dynamics
  This is mainly a container for the lists that hold the dynamics generators, propagators, and time evolution operators in each timeslot. The combining of dynamics generators is also complete by this object. Different subclasses support a range of types of quantum systems, including closed systems with unitary dynamics, systems with quadratic Hamiltonians that have Gaussian states and symplectic transforms, and a general subclass that can be used for open system dynamics with Lindbladian operators.

PulseGen
  There are many subclasses of pulse generators that generate different types of pulses as the initial amplitudes for the optimisation. Often the goal cannot be achieved from all starting conditions, and then typically some kind of random pulse is used and repeated optimisations are performed until the desired infidelity is reached or the minimum infidelity found is reported.
  There is a specific subclass that is used by the CRAB algorithm to generate the pulses based on the basis coefficients that are being optimised.

TerminationConditions
  This is simply a convenient place to hold all the properties that will determine when the single optimisation run terminates. Limits can be set for number of iterations, time, and of course the target infidelity.

Stats
  Performance data are optionally collected during the optimisation. This object is shared to a single location to store, calculate and report run statistics.

FidelityComputer
  The subclass of the fidelity computer determines the type of fidelity measure. These are closely linked to the type of dynamics in use. These are also the most commonly user customised subclasses.

PropagatorComputer
  This object computes propagators from one timeslot to the next and also the propagator gradient. The options are using the spectral decomposition or Frechet derivative, as discussed above.

TimeslotComputer
  Here the time evolution is computed by calling the methods of the other computer objects.

OptimResult
  The result of a pulse optimisation run is returned as an object with properties for the outcome in terms of the infidelity, reason for termination, performance statistics, final evolution, and more.

.. _pulseoptim-functions:

Using the pulseoptim functions
==============================
The simplest method for optimising a control pulse is to call one of the functions in the ``pulseoptim`` module. This automates the creation and configuration of the necessary objects, generation of initial pulses, running the optimisation and returning the result. There are functions specifically for unitary dynamics, and also specifically for the CRAB algorithm (GRAPE is the default). The ``optimise_pulse`` function can in fact be used for unitary dynamics and / or the CRAB algorithm, the more specific functions simply have parameter names that are more familiar in that application.

A semi-automated method is to use the ``create_optimizer_objects`` function to generate and configure all the objects, then manually set the initial pulse and call the optimisation. This would be more efficient when repeating runs with different starting conditions.
.. _dynamics:

******************************************
Time Evolution and Quantum System Dynamics
******************************************

.. toctree::
   :maxdepth: 2

   dynamics/dynamics-data.rst
   dynamics/dynamics-master.rst
   dynamics/dynamics-monte.rst
   dynamics/dynamics-photocurrent.rst
   dynamics/dynamics-stochastic.rst
   dynamics/dynamics-time.rst
   dynamics/dynamics-bloch-redfield.rst
   dynamics/dynamics-floquet.rst
   dynamics/dynamics-piqs.rst
   dynamics/dynamics-options.rst
.. _heom:

********************************
Hierarchical Equations of Motion
********************************

.. toctree::
   :maxdepth: 2

   heom/intro.rst
   heom/bosonic.rst
   heom/fermionic.rst
   heom/history.rst
   heom/references.rst
.. _saving:

**********************************
Saving QuTiP Objects and Data Sets
**********************************


With time-consuming calculations it is often necessary to store the results to files on disk, so it can be post-processed and archived. In QuTiP there are two facilities for storing data: Quantum objects can be stored to files and later read back as python pickles, and numerical data (vectors and matrices) can be exported as plain text files in for example CSV (comma-separated values), TSV (tab-separated values), etc. The former method is preferred when further calculations will be performed with the data, and the latter when the calculations are completed and data is to be imported into a post-processing tool (e.g. for generating figures).

Storing and loading QuTiP objects
=================================

To store and load arbitrary QuTiP related objects (:class:`qutip.Qobj`, :class:`qutip.solver.Result`, etc.) there are two functions: :func:`qutip.fileio.qsave` and :func:`qutip.fileio.qload`. The function :func:`qutip.fileio.qsave` takes an arbitrary object as first parameter and an optional filename as second parameter (default filename is `qutip_data.qu`). The filename extension is always `.qu`. The function :func:`qutip.fileio.qload` takes a mandatory filename as first argument and loads and returns the objects in the file.

To illustrate how these functions can be used, consider a simple calculation of the steadystate of the harmonic oscillator ::

    >>> a = destroy(10); H = a.dag() * a
    >>> c_ops = [np.sqrt(0.5) * a, np.sqrt(0.25) * a.dag()]
    >>> rho_ss = steadystate(H, c_ops)

The steadystate density matrix `rho_ss` is an instance of :class:`qutip.Qobj`. It can be stored to a file `steadystate.qu` using ::

    >>> qsave(rho_ss, 'steadystate')
    >>> !ls *.qu
    density_matrix_vs_time.qu  steadystate.qu

and it can later be loaded again, and used in further calculations ::

    >>> rho_ss_loaded = qload('steadystate')
    Loaded Qobj object:
    Quantum object: dims = [[10], [10]], shape = (10, 10), type = oper, isHerm = True
    >>> a = destroy(10)
    >>> np.testing.assert_almost_equal(expect(a.dag() * a, rho_ss_loaded), 0.9902248289345061)

The nice thing about the :func:`qutip.fileio.qsave` and :func:`qutip.fileio.qload` functions is that almost any object can be stored and load again later on. We can for example store a list of density matrices as returned by :func:`qutip.mesolve` ::

    >>> a = destroy(10); H = a.dag() * a ; c_ops = [np.sqrt(0.5) * a, np.sqrt(0.25) * a.dag()]
    >>> psi0 = rand_ket(10)
    >>> times = np.linspace(0, 10, 10)
    >>> dm_list = mesolve(H, psi0, times, c_ops, [])
    >>> qsave(dm_list, 'density_matrix_vs_time')

And it can then be loaded and used again, for example in an other program ::

    >>> dm_list_loaded = qload('density_matrix_vs_time')
    Loaded Result object:
    Result object with mesolve data.
    --------------------------------
    states = True
    num_collapse = 0
    >>> a = destroy(10)
    >>> expect(a.dag() * a, dm_list_loaded.states) # doctest: +SKIP
    array([4.63317086, 3.59150315, 2.90590183, 2.41306641, 2.05120716,
       1.78312503, 1.58357995, 1.4346382 , 1.32327398, 1.23991233])


Storing and loading datasets
============================

The :func:`qutip.fileio.qsave` and :func:`qutip.fileio.qload` are great, but the file format used is only understood by QuTiP (python) programs. When data must be exported to other programs the preferred method is to store the data in the commonly used plain-text file formats. With the QuTiP functions :func:`qutip.fileio.file_data_store` and :func:`qutip.fileio.file_data_read` we can store and load **numpy** arrays and matrices to files on disk using a deliminator-separated value format (for example comma-separated values CSV). Almost any program can handle this file format.

The :func:`qutip.fileio.file_data_store` takes two mandatory and three optional arguments:

>>> file_data_store(filename, data, numtype="complex", numformat="decimal", sep=",") # doctest: +SKIP

where `filename` is the name of the file, `data` is the data to be written to the file (must be a *numpy* array), `numtype` (optional) is a flag indicating numerical type that can take values `complex` or `real`, `numformat` (optional) specifies the numerical format that can take the values `exp` for the format `1.0e1` and `decimal` for the format `10.0`, and `sep` (optional) is an arbitrary single-character field separator (usually a tab, space, comma, semicolon, etc.).

A common use for the :func:`qutip.fileio.file_data_store` function is to store the expectation values of a set of operators for a sequence of times, e.g., as returned by the :func:`qutip.mesolve` function, which is what the following example does

.. plot::
    :context:

    >>> a = destroy(10); H = a.dag() * a ; c_ops = [np.sqrt(0.5) * a, np.sqrt(0.25) * a.dag()]
    >>> psi0 = rand_ket(10)
    >>> times = np.linspace(0, 100, 100)
    >>> medata = mesolve(H, psi0, times, c_ops, [a.dag() * a, a + a.dag(), -1j * (a - a.dag())])
    >>> np.shape(medata.expect)
    (3, 100)
    >>> times.shape
    (100,)
    >>> output_data = np.vstack((times, medata.expect))   # join time and expt data
    >>> file_data_store('expect.dat', output_data.T) # Note the .T for transpose!
    >>> with open("expect.dat", "r") as f:
    ...    print('\n'.join(f.readlines()[:10]))
    # Generated by QuTiP: 100x4 complex matrix in decimal format [',' separated values].
    0.0000000000+0.0000000000j,3.2109553666+0.0000000000j,0.3689771549+0.0000000000j,0.0185002867+0.0000000000j
    1.0101010101+0.0000000000j,2.6754598872+0.0000000000j,0.1298251132+0.0000000000j,-0.3303672956+0.0000000000j
    2.0202020202+0.0000000000j,2.2743186810+0.0000000000j,-0.2106241300+0.0000000000j,-0.2623894277+0.0000000000j
    3.0303030303+0.0000000000j,1.9726633457+0.0000000000j,-0.3037311621+0.0000000000j,0.0397330921+0.0000000000j
    4.0404040404+0.0000000000j,1.7435892209+0.0000000000j,-0.1126550232+0.0000000000j,0.2497182058+0.0000000000j
    5.0505050505+0.0000000000j,1.5687324121+0.0000000000j,0.1351622725+0.0000000000j,0.2018398581+0.0000000000j
    6.0606060606+0.0000000000j,1.4348632045+0.0000000000j,0.2143080535+0.0000000000j,-0.0067820038+0.0000000000j
    7.0707070707+0.0000000000j,1.3321818015+0.0000000000j,0.0950352763+0.0000000000j,-0.1630920429+0.0000000000j
    8.0808080808+0.0000000000j,1.2533244850+0.0000000000j,-0.0771210981+0.0000000000j,-0.1468923919+0.0000000000j


In this case we didn't really need to store both the real and imaginary parts, so instead we could use the ``numtype="real"`` option

.. plot::
   :context:

    >>> file_data_store('expect.dat', output_data.T, numtype="real")
    >>> with open("expect.dat", "r") as f:
    ...    print('\n'.join(f.readlines()[:5]))
    # Generated by QuTiP: 100x4 real matrix in decimal format [',' separated values].
    0.0000000000,3.2109553666,0.3689771549,0.0185002867
    1.0101010101,2.6754598872,0.1298251132,-0.3303672956
    2.0202020202,2.2743186810,-0.2106241300,-0.2623894277
    3.0303030303,1.9726633457,-0.3037311621,0.0397330921

and if we prefer scientific notation we can request that using the ``numformat="exp"`` option

.. plot::
    :context:

    >>> file_data_store('expect.dat', output_data.T, numtype="real", numformat="exp")

Loading data previously stored using :func:`qutip.fileio.file_data_store` (or some other software) is a even easier. Regardless of which deliminator was used, if data was stored as complex or real numbers, if it is in decimal or exponential form, the data can be loaded using the :func:`qutip.fileio.file_data_read`, which only takes the filename as mandatory argument.

.. plot::
    :context:

    input_data = file_data_read('expect.dat')
    plt.plot(input_data[:,0], input_data[:,1]);  # plot the data


(If a particularly obscure choice of deliminator was used it might be necessary to use the optional second argument, for example ``sep="_"`` if ``_`` is the deliminator).
.. _qip:

******************************
Quantum Information Processing
******************************

.. toctree::
   :maxdepth: 2

   qip/qip-basics.rst
   qip/qip-simulator.rst
   qip/qip-processor.rst
.. _measurement:

******************************
Measurement of Quantum Objects
******************************

.. note::
   New in QuTiP 4.6

.. _measurement-intro:

Introduction
------------

Measurement is a fundamental part of the standard formulation of quantum
mechanics and is the process by which classical readings are obtained from
a quantum object. Although the interpretation of the procedure is at times
contentious, the procedure itself is mathematically straightforward and is
described in many good introductory texts.

Here we will show you how to perform simple measurement operations on QuTiP
objects. The same functions :func:`~qutip.measurement.measure` and
:func:`~qutip.measurement.measurement_statistics` can be used
to handle both observable-style measurements and projective style measurements.

.. _measurement-basic:

Performing a basic measurement (Observable)
-------------------------------------------

First we need to select some states to measure. For now, let us create an *up*
state and a *down* state:

.. testcode::

   up = basis(2, 0)

   down = basis(2, 1)

which represent spin-1/2 particles with their spin pointing either up or down
along the z-axis.

We choose what to measure (in this case) by selecting a **measurement operator**.
For example,
we could select :func:`~qutip.sigmaz` which measures the z-component of the
spin of a spin-1/2 particle, or :func:`~qutip.sigmax` which measures the
x-component:

.. testcode::

   spin_z = sigmaz()

   spin_x = sigmax()

How do we know what these operators measure? The answer lies in the measurement
procedure itself:

* A quantum measurement transforms the state being measured by projecting it into
  one of the eigenvectors of the measurement operator.

* Which eigenvector to project onto is chosen probabilistically according to the
  square of the amplitude of the state in the direction of the eigenvector.

* The value returned by the measurement is the eigenvalue corresponding to the
  chosen eigenvector.

.. note::

   How to interpret this "random choosing" is the famous
   "quantum measurement problem".

The eigenvectors of `spin_z` are the states with their spin pointing either up
or down, so it measures the component of the spin along the z-axis.

The eigenvectors of `spin_x` are the states with their spin pointing either
left or right, so it measures the component of the spin along the x-axis.

When we measure our `up` and `down` states using the operator `spin_z`, we
always obtain:

.. testcode::

   from qutip.measurement import measure, measurement_statistics

   measure(up, spin_z) == (1.0, up)

   measure(down, spin_z) == (-1.0, down)

because `up` is the eigenvector of `spin_z` with eigenvalue `1.0` and `down`
is the eigenvector with eigenvalue `-1.0`. The minus signs are just an
arbitrary global phase -- `up` and `-up` represent the same quantum state.

Neither eigenvector has any component in the direction of the other (they are
orthogonal), so `measure(spin_z, up)` returns the state `up` 100% percent of the
time and `measure(spin_z, down)` returns the state `down` 100% of the time.

Note how :func:`~qutip.measurement.measure` returns a pair of values. The
first is the measured value, i.e. an eigenvalue of the operator (e.g. `1.0`),
and the second is the state of the quantum system after the measurement,
i.e. an eigenvector of the operator (e.g. `up`).

Now let us consider what happens if we measure the x-component of the spin
of `up`:

.. testcode::

   measure(up, spin_x)

The `up` state is not an eigenvector of `spin_x`. `spin_x` has two eigenvectors
which we will call `left` and `right`. The `up` state has equal components in
the direction of these two vectors, so measurement will select each of them
50% of the time.

These `left` and `right` states are:

.. testcode::

   left = (up - down).unit()

   right = (up + down).unit()

When `left` is chosen, the result of the measurement will be `(-1.0, -left)`.

When `right` is chosen, the result of measurement with be `(1.0, right)`.

.. note::

  When :func:`~qutip.measurement.measure` is invoked with the second argument
  being an observable, it acts as an alias to
  :func:`~qutip.measurement.measure_observable`.

Performing a basic measurement (Projective)
-------------------------------------------

We can also choose what to measure by specifying a *list of projection operators*. For
example, we could select the projection operators :math:`\ket{0} \bra{0}` and
:math:`\ket{1} \bra{1}` which measure the state in the :math:`\ket{0}, \ket{1}`
basis. Note that these projection operators are simply the projectors determined by
the eigenstates of the :func:`~qutip.sigmaz` operator.

.. testcode::

   Z0, Z1 = ket2dm(basis(2, 0)), ket2dm(basis(2, 1))

The probabilities and respective output state
are calculated for each projection operator.

.. testcode::

   measure(up, [Z0, Z1]) == (0, up)

   measure(down, [Z0, Z1]) == (1, down)

In this case, the projection operators are conveniently eigenstates corresponding
to subspaces of dimension :math:`1`. However, this might not be
the case, in which case it is not possible to have unique eigenvalues for each
eigenstate. Suppose we want to measure only the first
qubit in a two-qubit system. Consider the two qubit state :math:`\ket{0+}`

.. testcode::

   state_0 = basis(2, 0)

   state_plus = (basis(2, 0) + basis(2, 1)).unit()

   state_0plus = tensor(state_0, state_plus)

Now, suppose we want to measure only the first qubit in the computational basis.
We can do that by measuring with the projection operators
:math:`\ket{0}\bra{0} \otimes I` and  :math:`\ket{1}\bra{1} \otimes I`.

.. testcode::

   PZ1 = [tensor(Z0, identity(2)), tensor(Z1, identity(2))]

   PZ2 = [tensor(identity(2), Z0), tensor(identity(2), Z1)]

Now, as in the previous example, we can measure by supplying a list of projection operators
and the state.

.. testcode::

    measure(state_0plus, PZ1) == (0, state_0plus)

The output of the measurement is the index of the measurement outcome as well
as the output state on the full Hilbert space of the input state. It is crucial to
note that we do not discard the measured qubit after measurement (as opposed to
when measuring on quantum hardware).

.. note::

  When :func:`~qutip.measurement.measure` is invoked with the second argument
  being a list of projectors, it acts as an alias to
  :func:`~qutip.measurement.measure_povm`.

The :func:`~qutip.measurement.measure` function can perform measurements on
density matrices too. You can read about these and other details at
:func:`~qutip.measurement.measure_povm` and :func:`~qutip.measurement.measure_observable`.

Now you know how to measure quantum states in QuTiP!

.. _measurement-statistics:

Obtaining measurement statistics(Observable)
--------------------------------------------

You've just learned how to perform measurements in QuTiP, but you've also
learned that measurements are probabilistic. What if instead of just making
a single measurement, we want to determine the probability distribution of
a large number of measurements?

One way would be to repeat the measurement many times -- and this is what
happens in many quantum experiments. In QuTiP one could simulate this using:

.. testcode::
    :hide:

    np.random.seed(42)

.. testcode::

   results = {1.0: 0, -1.0: 0}  # 1 and -1 are the possible outcomes
   for _ in range(1000):
      value, new_state = measure(up, spin_x)
      results[round(value)] += 1
   print(results)

**Output**:

.. testoutput::

   {1.0: 497, -1.0: 503}

which measures the x-component of the spin of the `up` state `1000` times and
stores the results in a dictionary. Afterwards we expect to have seen the
result `1.0` (i.e. left) roughly 500 times and the result `-1.0` (i.e. right)
roughly 500 times, but, of course, the number of each will vary slightly
each time we run it.

But what if we want to know the distribution of results precisely? In a
physical system, we would have to perform the measurement many many times,
but in QuTiP we can peak at the state itself and determine the probability
distribution of the outcomes exactly in a single line:

.. doctest::
    :hide:

    >>> np.random.seed(42)

.. doctest::

   >>> eigenvalues, eigenstates, probabilities = measurement_statistics(up, spin_x)

   >>> eigenvalues # doctest: +NORMALIZE_WHITESPACE
   array([-1., 1.])

   >>> eigenstates # doctest: +NORMALIZE_WHITESPACE
   array([Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[ 0.70710678]
    [-0.70710678]],
          Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[0.70710678]
    [0.70710678]]], dtype=object)

   >>> probabilities  # doctest: +NORMALIZE_WHITESPACE
   [0.5000000000000001, 0.4999999999999999]

The :func:`~qutip.measurement.measurement_statistics` function then returns three values
when called with a single observable:

- `eigenvalues` is an array of eigenvalues of the measurement operator, i.e.
  a list of the possible measurement results. In our example
  the value is `array([-1., -1.])`.

- `eigenstates` is an array of the eigenstates of the measurement operator, i.e.
  a list of the possible final states after the measurement is complete.
  Each element of the array is a :obj:`~qutip.Qobj`.

- `probabilities` is a list of the probabilities of each measurement result.
  In our example the value is `[0.5, 0.5]` since the `up` state has equal
  probability of being measured to be in the left (`-1.0`) or
  right (`1.0`) eigenstates.

All three lists are in the same order -- i.e. the first eigenvalue is
`eigenvalues[0]`, its corresponding eigenstate is `eigenstates[0]`, and
its probability is `probabilities[0]`, and so on.

.. note::

   When :func:`~qutip.measurement.measurement_statistics`
   is invoked with the second argument
   being an observable, it acts as an alias to
   :func:`~qutip.measurement.measurement_statistics_observable`.


Obtaining measurement statistics(Projective)
--------------------------------------------

Similarly, when we want to obtain measurement statistics for projection operators,
we can use the `measurement_statistics` function with the second argument being a list of projectors.
Consider again, the state :math:`\ket{0+}`.
Suppose, now we want to obtain the measurement outcomes for the second qubit. We
must use the projectors specified earlier by `PZ2` which allow us to measure only
on the second qubit. Since the second qubit has the state :math:`\ket{+}`, we get
the following result.

.. testcode::

   collapsed_states, probabilities = measurement_statistics(state_0plus, PZ2)

   print(collapsed_states)

**Output**:

.. testoutput::
   :options: +NORMALIZE_WHITESPACE

   [Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
    Qobj data =
    [[1.]
     [0.]
     [0.]
     [0.]], Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
    Qobj data =
    [[0.]
     [1.]
     [0.]
     [0.]]]

.. testcode::

   print(probabilities)

**Output**:

.. testoutput::
   :options: +NORMALIZE_WHITESPACE

   [0.4999999999999999, 0.4999999999999999]

The function :func:`~qutip.measurement.measurement_statistics` then returns two values:

* `collapsed_states` is an array of the possible final states after the
  measurement is complete. Each element of the array is a :obj:`~qutip.Qobj`.

* `probabilities` is a list of the probabilities of each measurement outcome.

Note that the collapsed_states are exactly :math:`\ket{00}` and :math:`\ket{01}`
with equal probability, as expected. The two lists are in the same order.

.. note::

   When :func:`~qutip.measurement.measurement_statistics`
   is invoked with the second argument
   being a list of projectors, it acts as an alias to
   :func:`~qutip.measurement.measurement_statistics_povm`.

The :func:`~qutip.measurement.measurement_statistics` function can provide statistics for measurements
of density matrices too.
You can read about these and other details at
:func:`~qutip.measurement.measurement_statistics_observable`
and :func:`~qutip.measurement.measurement_statistics_povm`.

Furthermore, the :func:`~qutip.measurement.measure_povm`
and :func:`~qutip.measurement.measurement_statistics_povm` functions can
handle POVM measurements which are more general than projective measurements.
.. _settings:

*********************************
Modifying Internal QuTiP Settings
*********************************

.. _settings-params:

User Accessible Parameters
==========================

In this section we show how to modify a few of the internal parameters used by QuTiP. The settings that can be modified are given in the following table:

.. tabularcolumns:: | p{3cm} | p{5cm} | p{5cm} |

.. cssclass:: table-striped

+-------------------------------+-------------------------------------------+-----------------------------+
| Setting                       | Description                               | Options                     |
+===============================+===========================================+=============================+
| `auto_herm`                   | Automatically calculate the hermicity of  | True / False                |
|                               | quantum objects.                          |                             |
+-------------------------------+-------------------------------------------+-----------------------------+
| `auto_tidyup`                 | Automatically tidyup quantum objects.     | True / False                |
+-------------------------------+-------------------------------------------+-----------------------------+
| `auto_tidyup_atol`            | Tolerance used by tidyup                  | any `float` value > 0       |
+-------------------------------+-------------------------------------------+-----------------------------+
| `atol`                        | General tolerance                         | any `float` value > 0       |
+-------------------------------+-------------------------------------------+-----------------------------+
| `num_cpus`                    | Number of CPU's used for multiprocessing. | `int` between 1 and # cpu's |
+-------------------------------+-------------------------------------------+-----------------------------+
| `debug`                       | Show debug printouts.                     | True / False                |
+-------------------------------+-------------------------------------------+-----------------------------+
| `openmp_thresh`               | NNZ matrix must have for OPENMP.          | Int                         |
+-------------------------------+-------------------------------------------+-----------------------------+

.. _settings-usage:

Example: Changing Settings
==========================

The two most important settings are ``auto_tidyup`` and ``auto_tidyup_atol`` as they control whether the small elements of a quantum object should be removed, and what number should be considered as the cut-off tolerance. Modifying these, or any other parameters, is quite simple::

>>> qutip.settings.auto_tidyup = False

These settings will be used for the current QuTiP session only and will need to be modified again when restarting QuTiP.  If running QuTiP from a script file, then place the `qutip.settings.xxxx` commands immediately after `from qutip import *` at the top of the script file.  If you want to reset the parameters back to their default values then call the reset command::

>>> qutip.settings.reset()

Persistent Settings
===================

When QuTiP is imported, it looks for a file named ``qutiprc`` in a folder called ``.qutip`` user's home directory. If this file is found, it will be loaded and overwrite the QuTiP default settings, which allows for persistent changes in the QuTiP settings to be made. A sample ``qutiprc`` file is show below. The syntax is a simple key-value format, where the keys and possible values are described in the table above::

    [qutip]
    auto_tidyup=True
    auto_herm=True
    auto_tidyup_atol=1e-12
    num_cpus=4
    debug=False

Note that the ``openmp_thresh`` value is automatically generatd by QuTiP.  It is also possible to set a specific compiler for QuTiP to use when generating runtime Cython code for time-dependent problems.  For example, the following section in the ``qutiprc`` file will set the compiler to be ``clang-3.9``::

    [compiler]
    cc = clang-3.9
    cxx = clang-3.9


.. _correlation:

******************************
Two-time correlation functions
******************************

With the QuTiP time-evolution functions (for example :func:`qutip.mesolve` and :func:`qutip.mcsolve`), a state vector or density matrix can be evolved from an initial state at :math:`t_0` to an arbitrary time :math:`t`, :math:`\rho(t)=V(t, t_0)\left\{\rho(t_0)\right\}`, where :math:`V(t, t_0)` is the propagator defined by the equation of motion. The resulting density matrix can then be used to evaluate the expectation values of arbitrary combinations of *same-time* operators.

To calculate *two-time* correlation functions on the form :math:`\left<A(t+\tau)B(t)\right>`, we can use the quantum regression theorem (see, e.g., [Gar03]_) to write

.. math::

    \left<A(t+\tau)B(t)\right> = {\rm Tr}\left[A V(t+\tau, t)\left\{B\rho(t)\right\}\right]
                               = {\rm Tr}\left[A V(t+\tau, t)\left\{BV(t, 0)\left\{\rho(0)\right\}\right\}\right]

We therefore first calculate :math:`\rho(t)=V(t, 0)\left\{\rho(0)\right\}` using one of the QuTiP evolution solvers with :math:`\rho(0)` as initial state, and then again use the same solver to calculate :math:`V(t+\tau, t)\left\{B\rho(t)\right\}` using :math:`B\rho(t)` as initial state.

Note that if the initial state is the steady state, then :math:`\rho(t)=V(t, 0)\left\{\rho_{\rm ss}\right\}=\rho_{\rm ss}` and

.. math::

    \left<A(t+\tau)B(t)\right> = {\rm Tr}\left[A V(t+\tau, t)\left\{B\rho_{\rm ss}\right\}\right]
                               = {\rm Tr}\left[A V(\tau, 0)\left\{B\rho_{\rm ss}\right\}\right] = \left<A(\tau)B(0)\right>,

which is independent of :math:`t`, so that we only have one time coordinate :math:`\tau`.

QuTiP provides a family of functions that assists in the process of calculating two-time correlation functions. The available functions and their usage is shown in the table below. Each of these functions can use one of the following evolution solvers: Master-equation, Exponential series and the Monte-Carlo. The choice of solver is defined by the optional argument ``solver``.

.. cssclass:: table-striped

+----------------------------------------------+--------------------------------------------------+
| QuTiP function                               | Correlation function                             |
+==============================================+==================================================+
|                                              | :math:`\left<A(t+\tau)B(t)\right>` or            |
| :func:`qutip.correlation.correlation_2op_2t` | :math:`\left<A(t)B(t+\tau)\right>`.              |
+----------------------------------------------+--------------------------------------------------+
|                                              | :math:`\left<A(\tau)B(0)\right>` or              |
| :func:`qutip.correlation.correlation_2op_1t` | :math:`\left<A(0)B(\tau)\right>`.                |
+----------------------------------------------+--------------------------------------------------+
| :func:`qutip.correlation.correlation_3op_1t` | :math:`\left<A(0)B(\tau)C(0)\right>`.            |
+----------------------------------------------+--------------------------------------------------+
| :func:`qutip.correlation.correlation_3op_2t` | :math:`\left<A(t)B(t+\tau)C(t)\right>`.          |
+----------------------------------------------+--------------------------------------------------+

The most common use-case is to calculate correlation functions of the kind :math:`\left<A(\tau)B(0)\right>`, in which case we use the correlation function solvers that start from the steady state, e.g., the :func:`qutip.correlation.correlation_2op_1t` function. These correlation function solvers return a vector or matrix (in general complex) with the correlations as a function of the delays times.

.. _correlation-steady:

Steadystate correlation function
================================

The following code demonstrates how to calculate the :math:`\left<x(t)x(0)\right>` correlation for a leaky cavity with three different relaxation rates.

.. plot::
    :context:

    times = np.linspace(0,10.0,200)
    a = destroy(10)
    x = a.dag() + a
    H = a.dag() * a

    corr1 = correlation_2op_1t(H, None, times, [np.sqrt(0.5) * a], x, x)
    corr2 = correlation_2op_1t(H, None, times, [np.sqrt(1.0) * a], x, x)
    corr3 = correlation_2op_1t(H, None, times, [np.sqrt(2.0) * a], x, x)

    plt.figure()
    plt.plot(times, np.real(corr1), times, np.real(corr2), times, np.real(corr3))
    plt.legend(['0.5','1.0','2.0'])
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Correlation $\left<x(t)x(0)\right>$')
    plt.show()


Emission spectrum
=================

Given a correlation function :math:`\left<A(\tau)B(0)\right>` we can define the corresponding power spectrum as

.. math::

    S(\omega) = \int_{-\infty}^{\infty} \left<A(\tau)B(0)\right> e^{-i\omega\tau} d\tau.

In QuTiP, we can calculate :math:`S(\omega)` using either :func:`qutip.correlation.spectrum_ss`, which first calculates the correlation function using one of the time-dependent solvers and then performs the Fourier transform semi-analytically, or we can use the function :func:`qutip.correlation.spectrum_correlation_fft` to numerically calculate the Fourier transform of a given correlation data using FFT.

The following example demonstrates how these two functions can be used to obtain the emission power spectrum.

.. plot:: guide/scripts/spectrum_ex1.py
   :width: 5.0in
   :include-source:

.. _correlation-spectrum:


Non-steadystate correlation function
====================================

More generally, we can also calculate correlation functions of the kind :math:`\left<A(t_1+t_2)B(t_1)\right>`, i.e., the correlation function of a system that is not in its steady state. In QuTiP, we can evaluate such correlation functions using the function :func:`qutip.correlation.correlation_2op_2t`. The default behavior of this function is to return a matrix with the correlations as a function of the two time coordinates (:math:`t_1` and :math:`t_2`).

.. plot:: guide/scripts/correlation_ex2.py
   :width: 5.0in
   :include-source:

However, in some cases we might be interested in the correlation functions on the form :math:`\left<A(t_1+t_2)B(t_1)\right>`, but only as a function of time coordinate :math:`t_2`. In this case we can also use the :func:`qutip.correlation.correlation_2op_2t` function, if we pass the density matrix at time :math:`t_1` as second argument, and `None` as third argument. The :func:`qutip.correlation.correlation_2op_2t` function then returns a vector with the correlation values corresponding to the times in `taulist` (the fourth argument).

Example: first-order optical coherence function
-----------------------------------------------

This example demonstrates how to calculate a correlation function on the form :math:`\left<A(\tau)B(0)\right>` for a non-steady initial state. Consider an oscillator that is interacting with a thermal environment. If the oscillator initially is in a coherent state, it will gradually decay to a thermal (incoherent) state. The amount of coherence can be quantified using the first-order optical coherence function :math:`g^{(1)}(\tau) = \frac{\left<a^\dagger(\tau)a(0)\right>}{\sqrt{\left<a^\dagger(\tau)a(\tau)\right>\left<a^\dagger(0)a(0)\right>}}`. For a coherent state :math:`|g^{(1)}(\tau)| = 1`, and for a completely incoherent (thermal) state :math:`g^{(1)}(\tau) = 0`. The following code calculates and plots :math:`g^{(1)}(\tau)` as a function of :math:`\tau`.

.. plot:: guide/scripts/correlation_ex3.py
   :width: 5.0in
   :include-source:

For convenience, the steps for calculating the first-order coherence function have been collected in the function :func:`qutip.correlation.coherence_function_g1`.

Example: second-order optical coherence function
------------------------------------------------

The second-order optical coherence function, with time-delay :math:`\tau`, is defined as

.. math::

    \displaystyle g^{(2)}(\tau) = \frac{\langle a^\dagger(0)a^\dagger(\tau)a(\tau)a(0)\rangle}{\langle a^\dagger(0)a(0)\rangle^2}

For a coherent state :math:`g^{(2)}(\tau) = 1`, for a thermal state :math:`g^{(2)}(\tau=0) = 2` and it decreases as a function of time (bunched photons, they tend to appear together), and for a Fock state with :math:`n` photons :math:`g^{(2)}(\tau = 0) = n(n - 1)/n^2 < 1` and it increases with time (anti-bunched photons, more likely to arrive separated in time).

To calculate this type of correlation function with QuTiP, we can use :func:`qutip.correlation.correlation_3op_1t`, which computes a correlation function on the form :math:`\left<A(0)B(\tau)C(0)\right>` (three operators, one delay-time vector).
We first have to combine the central two operators into one single one as they are evaluated at the same time, e.g. here we do :math:`a^\dagger(\tau)a(\tau) = (a^\dagger a)(\tau)`.

The following code calculates and plots :math:`g^{(2)}(\tau)` as a function of :math:`\tau` for a coherent, thermal and Fock state.

.. plot:: guide/scripts/correlation_ex4.py
   :width: 5.0in
   :include-source:

For convenience, the steps for calculating the second-order coherence function have been collected in the function :func:`qutip.correlation.coherence_function_g2`.
.. _random:

********************************************
Generating Random Quantum States & Operators
********************************************

.. testsetup:: [random]

   from qutip import rand_herm, rand_dm, rand_super_bcsz, rand_dm_ginibre

QuTiP includes a collection of random state, unitary and channel generators for simulations, Monte Carlo evaluation, theorem evaluation, and code testing.
Each of these objects can be sampled from one of several different distributions including the default distributions
used by QuTiP versions prior to 3.2.0.

For example, a random Hermitian operator can be sampled by calling `rand_herm` function:

.. doctest:: [random]
    :hide:

    >>> np.random.seed(42)

.. doctest:: [random]

  >>> rand_herm(5) # doctest: +NORMALIZE_WHITESPACE
  Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
  Qobj data =
  [[-0.25091976+0.j          0.        +0.j          0.        +0.j
    -0.21793701+0.47037633j -0.23212846-0.61607187j]
   [ 0.        +0.j         -0.88383278+0.j          0.836086  -0.23956218j
    -0.09464275+0.45370863j -0.15243356+0.65392096j]
   [ 0.        +0.j          0.836086  +0.23956218j  0.66488528+0.j
    -0.26290446+0.64984451j -0.52603038-0.07991553j]
   [-0.21793701-0.47037633j -0.09464275-0.45370863j -0.26290446-0.64984451j
    -0.13610996+0.j         -0.34240902-0.2879303j ]
   [-0.23212846+0.61607187j -0.15243356-0.65392096j -0.52603038+0.07991553j
    -0.34240902+0.2879303j   0.        +0.j        ]]



.. tabularcolumns:: | p{2cm} | p{3cm} | c |

.. cssclass:: table-striped

+-------------------------------+--------------------------------------------+------------------------------------------+
| Random Variable Type          | Sampling Functions                         | Dimensions                               |
+===============================+============================================+==========================================+
| State vector (``ket``)        | `rand_ket`, `rand_ket_haar`                | :math:`N \times 1`                       |
+-------------------------------+--------------------------------------------+------------------------------------------+
| Hermitian operator (``oper``) | `rand_herm`                                | :math:`N \times 1`                       |
+-------------------------------+--------------------------------------------+------------------------------------------+
| Density operator (``oper``)   | `rand_dm`, `rand_dm_hs`, `rand_dm_ginibre` | :math:`N \times N`                       |
+-------------------------------+--------------------------------------------+------------------------------------------+
| Unitary operator (``oper``)   | `rand_unitary`, `rand_unitary_haar`        | :math:`N \times N`                       |
+-------------------------------+--------------------------------------------+------------------------------------------+
| CPTP channel (``super``)      | `rand_super`, `rand_super_bcsz`            | :math:`(N \times N) \times (N \times N)` |
+-------------------------------+--------------------------------------------+------------------------------------------+

In all cases, these functions can be called with a single parameter :math:`N` that specifies the dimension of the relevant Hilbert space. The optional
``dims`` keyword argument allows for the dimensions of a random state, unitary or channel to be broken down into subsystems.

.. doctest:: [random]

    >>> rand_super_bcsz(7).dims
    [[[7], [7]], [[7], [7]]]
    >>> rand_super_bcsz(6, dims=[[[2, 3], [2, 3]], [[2, 3], [2, 3]]]).dims
    [[[2, 3], [2, 3]], [[2, 3], [2, 3]]]

Several of the distributions supported by QuTiP support additional parameters as well, namely *density* and *rank*. In particular,
the `rand_herm` and `rand_dm` functions return quantum objects such that a fraction of the elements are identically equal to zero.
The ratio of nonzero elements is passed as the ``density`` keyword argument. By contrast, the `rand_dm_ginibre` and
`rand_super_bcsz` take as an argument the rank of the generated object, such that passing ``rank=1`` returns a random
pure state or unitary channel, respectively. Passing ``rank=None`` specifies that the generated object should be
full-rank for the given dimension.

For example,

.. doctest:: [random]
    :hide:

    >>> np.random.seed(42)

.. doctest:: [random]

   >>> rand_dm(5, density=0.5)
   Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
   Qobj data =
   [[ 0.05157906+0.j          0.04491736+0.01043329j  0.06966148+0.00344713j
      0.        +0.j          0.04031493-0.01886791j]
    [ 0.04491736-0.01043329j  0.33632352+0.j         -0.08046093+0.02954712j
      0.0037455 +0.03940256j -0.05679126-0.01322392j]
    [ 0.06966148-0.00344713j -0.08046093-0.02954712j  0.2938209 +0.j
      0.0029377 +0.04463531j  0.05318743-0.02817689j]
    [ 0.        +0.j          0.0037455 -0.03940256j  0.0029377 -0.04463531j
      0.22553181+0.j          0.01657495+0.06963845j]
    [ 0.04031493+0.01886791j -0.05679126+0.01322392j  0.05318743+0.02817689j
      0.01657495-0.06963845j  0.09274471+0.j        ]]

   >>> rand_dm_ginibre(5, rank=2)
   Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
   Qobj data =
   [[ 0.07318288+2.60675616e-19j  0.10426866-6.63115850e-03j
     -0.05377455-2.66949369e-02j -0.01623153+7.66824687e-02j
     -0.12255602+6.11342416e-02j]
    [ 0.10426866+6.63115850e-03j  0.30603789+1.44335373e-18j
     -0.03129486-4.16194216e-03j -0.09832531+1.74110000e-01j
     -0.27176358-4.84608761e-02j]
    [-0.05377455+2.66949369e-02j -0.03129486+4.16194216e-03j
      0.07055265-8.76912454e-19j -0.0183289 -2.72720794e-02j
      0.01196277-1.01037189e-01j]
    [-0.01623153-7.66824687e-02j -0.09832531-1.74110000e-01j
     -0.0183289 +2.72720794e-02j  0.14168414-1.51340961e-19j
      0.07847628+2.07735199e-01j]
    [-0.12255602-6.11342416e-02j -0.27176358+4.84608761e-02j
      0.01196277+1.01037189e-01j  0.07847628-2.07735199e-01j
      0.40854244-6.75775934e-19j]]




See the API documentation: :ref:`functions-rand` for details.

.. warning::

    When using the ``density`` keyword argument, setting the density too low may result in not enough diagonal elements to satisfy trace
    constraints.

Random objects with a given eigen spectrum
==========================================

It is also possible to generate random Hamiltonian (``rand_herm``) and densitiy matrices (``rand_dm``) with a given eigen spectrum.  This is done by passing an array of eigenvalues as the first argument to either function.  For example,

.. doctest:: [random]
    :hide:

    >>> np.random.seed(42)

.. doctest:: [random]

   >>> eigs = np.arange(5)

   >>> H = rand_herm(eigs, density=0.5)

   >>> H # doctest: +NORMALIZE_WHITESPACE
   Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
   Qobj data =
   [[ 2.51387054-5.55111512e-17j  0.81161447+2.02283642e-01j
      0.        +0.00000000e+00j  0.875     +3.35634092e-01j
      0.81161447+2.02283642e-01j]
    [ 0.81161447-2.02283642e-01j  1.375     +0.00000000e+00j
      0.        +0.00000000e+00j -0.76700198+5.53011066e-01j
      0.375     +0.00000000e+00j]
    [ 0.        +0.00000000e+00j  0.        +0.00000000e+00j
      2.        +0.00000000e+00j  0.        +0.00000000e+00j
      0.        +0.00000000e+00j]
    [ 0.875     -3.35634092e-01j -0.76700198-5.53011066e-01j
      0.        +0.00000000e+00j  2.73612946+0.00000000e+00j
     -0.76700198-5.53011066e-01j]
    [ 0.81161447-2.02283642e-01j  0.375     +0.00000000e+00j
      0.        +0.00000000e+00j -0.76700198+5.53011066e-01j
      1.375     +0.00000000e+00j]]


   >>> H.eigenenergies() # doctest: +NORMALIZE_WHITESPACE
   array([7.70647994e-17, 1.00000000e+00, 2.00000000e+00, 3.00000000e+00,
       4.00000000e+00])


In order  to generate a random object with a given spectrum QuTiP applies a series of random complex Jacobi rotations.  This technique requires many steps to build the desired quantum object, and is thus suitable only for objects with Hilbert dimensionality :math:`\lesssim 1000`.



Composite random objects
========================

In many cases, one is interested in generating random quantum objects that correspond to composite systems generated using the :func:`qutip.tensor.tensor` function.  Specifying the tensor structure of a quantum object is done using the `dims` keyword argument in the same fashion as one would do for a :class:`qutip.Qobj` object:

.. doctest:: [random]
    :hide:

    >>> np.random.seed(42)

.. doctest:: [random]

   >>> rand_dm(4, 0.5, dims=[[2,2], [2,2]]) # doctest: +NORMALIZE_WHITESPACE
   Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
   Qobj data =
   [[ 0.13622928+0.j          0.        +0.j          0.01180807-0.01739166j
      0.        +0.j        ]
    [ 0.        +0.j          0.14600238+0.j          0.10335328+0.21790786j
     -0.00426027-0.02193627j]
    [ 0.01180807+0.01739166j  0.10335328-0.21790786j  0.57566072+0.j
     -0.0670631 +0.04124094j]
    [ 0.        +0.j         -0.00426027+0.02193627j -0.0670631 -0.04124094j
      0.14210761+0.j        ]]
.. _tensor:

******************************************
Using Tensor Products and Partial Traces
******************************************


.. _tensor-products:

Tensor products
===============

To describe the states of multipartite quantum systems - such as two coupled qubits, a qubit coupled to an oscillator, etc. - we need to expand the Hilbert space by taking the tensor product of the state vectors for each of the system components. Similarly, the operators acting on the state vectors in the combined Hilbert space (describing the coupled system) are formed by taking the tensor product of the individual operators.

In QuTiP the function :func:`qutip.tensor.tensor` is used to accomplish this task. This function takes as argument a collection::

>>> tensor(op1, op2, op3) # doctest: +SKIP

or a ``list``::

>>> tensor([op1, op2, op3]) # doctest: +SKIP

of state vectors *or* operators and returns a composite quantum object for the combined Hilbert space. The function accepts an arbitrary number of states or operators as argument. The type returned quantum object is the same as that of the input(s).

For example, the state vector describing two qubits in their ground states is formed by taking the tensor product of the two single-qubit ground state vectors:

.. testcode:: [tensor]

    print(tensor(basis(2, 0), basis(2, 0)))

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
    Qobj data =
    [[1.]
     [0.]
     [0.]
     [0.]]

or equivalently using the ``list`` format:

.. testcode:: [tensor]

    print(tensor([basis(2, 0), basis(2, 0)]))

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
    Qobj data =
    [[1.]
     [0.]
     [0.]
     [0.]]

This is straightforward to generalize to more qubits by adding more component state vectors in the argument list to the :func:`qutip.tensor.tensor` function, as illustrated in the following example:

.. testcode:: [tensor]

    print(tensor((basis(2, 0) + basis(2, 1)).unit(), (basis(2, 0) + basis(2, 1)).unit(), basis(2, 0)))

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
    Qobj data =
    [[0.5]
     [0. ]
     [0.5]
     [0. ]
     [0.5]
     [0. ]
     [0.5]
     [0. ]]


This state is slightly more complicated, describing two qubits in a superposition between the up and down states, while the third qubit is in its ground state.

To construct operators that act on an extended Hilbert space of a combined system, we similarly pass a list of operators for each component system to the :func:`qutip.tensor.tensor` function. For example, to form the operator that represents the simultaneous action of the :math:`\sigma_x` operator on two qubits:

.. testcode:: [tensor]

    print(tensor(sigmax(), sigmax()))

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
    Qobj data =
    [[0. 0. 0. 1.]
     [0. 0. 1. 0.]
     [0. 1. 0. 0.]
     [1. 0. 0. 0.]]

To create operators in a combined Hilbert space that only act on a single component, we take the tensor product of the operator acting on the subspace of interest, with the identity operators corresponding to the components that are to be unchanged. For example, the operator that represents :math:`\sigma_z` on the first qubit in a two-qubit system, while leaving the second qubit unaffected:

.. testcode:: [tensor]

    print(tensor(sigmaz(), identity(2)))

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
    Qobj data =
    [[ 1.  0.  0.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0. -1.  0.]
     [ 0.  0.  0. -1.]]


.. _tensor-product-example:

Example: Constructing composite Hamiltonians
============================================

The :func:`qutip.tensor.tensor` function is extensively used when constructing Hamiltonians for composite systems. Here we'll look at some simple examples.

.. _tensor-product-example-2qubits:

Two coupled qubits
------------------

First, let's consider a system of two coupled qubits. Assume that both the qubits have equal energy splitting, and that the qubits are coupled through a :math:`\sigma_x\otimes\sigma_x` interaction with strength g = 0.05 (in units where the bare qubit energy splitting is unity). The Hamiltonian describing this system is:

.. testcode:: [tensor]

    H = tensor(sigmaz(), identity(2)) + tensor(identity(2), sigmaz()) + 0.05 * tensor(sigmax(), sigmax())

    print(H)

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
    Qobj data =
    [[ 2.    0.    0.    0.05]
     [ 0.    0.    0.05  0.  ]
     [ 0.    0.05  0.    0.  ]
     [ 0.05  0.    0.   -2.  ]]

.. _tensor-product-example-3qubits:

Three coupled qubits
--------------------

The two-qubit example is easily generalized to three coupled qubits:

.. testcode:: [tensor]

    H = (tensor(sigmaz(), identity(2), identity(2)) + tensor(identity(2), sigmaz(), identity(2)) + tensor(identity(2), identity(2), sigmaz()) + 0.5 * tensor(sigmax(), sigmax(), identity(2)) + 0.25 * tensor(identity(2), sigmax(), sigmax()))

    print(H)

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2, 2], [2, 2, 2]], shape = (8, 8), type = oper, isherm = True
    Qobj data =
    [[ 3.    0.    0.    0.25  0.    0.    0.5   0.  ]
     [ 0.    1.    0.25  0.    0.    0.    0.    0.5 ]
     [ 0.    0.25  1.    0.    0.5   0.    0.    0.  ]
     [ 0.25  0.    0.   -1.    0.    0.5   0.    0.  ]
     [ 0.    0.    0.5   0.    1.    0.    0.    0.25]
     [ 0.    0.    0.    0.5   0.   -1.    0.25  0.  ]
     [ 0.5   0.    0.    0.    0.    0.25 -1.    0.  ]
     [ 0.    0.5   0.    0.    0.25  0.    0.   -3.  ]]


.. _tensor-product-example-jcmodel:

A two-level system coupled to a cavity: The Jaynes-Cummings model
-------------------------------------------------------------------

The simplest possible quantum mechanical description for light-matter interaction is encapsulated in the Jaynes-Cummings model, which describes the coupling between a two-level atom and a single-mode electromagnetic field (a cavity mode). Denoting the energy splitting of the atom and cavity ``omega_a`` and ``omega_c``, respectively, and the atom-cavity interaction strength ``g``, the Jaynes-Cummings Hamiltonian can be constructed as:

.. testcode:: [tensor]

    N = 10

    omega_a = 1.0

    omega_c = 1.25

    g = 0.05

    a = tensor(identity(2), destroy(N))

    sm = tensor(destroy(2), identity(N))

    sz = tensor(sigmaz(), identity(N))

    H = 0.5 * omega_a * sz + omega_c * a.dag() * a + g * (a.dag() * sm + a * sm.dag())

    print(H)

**Output**:

.. testoutput:: [tensor]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 10], [2, 10]], shape = (20, 20), type = oper, isherm = True
    Qobj data =
    [[ 0.5         0.          0.          0.          0.          0.
       0.          0.          0.          0.          0.          0.
       0.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          1.75        0.          0.          0.          0.
       0.          0.          0.          0.          0.05        0.
       0.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          3.          0.          0.          0.
       0.          0.          0.          0.          0.          0.07071068
       0.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          4.25        0.          0.
       0.          0.          0.          0.          0.          0.
       0.08660254  0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          5.5         0.
       0.          0.          0.          0.          0.          0.
       0.          0.1         0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          6.75
       0.          0.          0.          0.          0.          0.
       0.          0.          0.1118034   0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       8.          0.          0.          0.          0.          0.
       0.          0.          0.          0.12247449  0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          9.25        0.          0.          0.          0.
       0.          0.          0.          0.          0.13228757  0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.         10.5         0.          0.          0.
       0.          0.          0.          0.          0.          0.14142136
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.          0.         11.75        0.          0.
       0.          0.          0.          0.          0.          0.
       0.15        0.        ]
     [ 0.          0.05        0.          0.          0.          0.
       0.          0.          0.          0.         -0.5         0.
       0.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.07071068  0.          0.          0.
       0.          0.          0.          0.          0.          0.75
       0.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.08660254  0.          0.
       0.          0.          0.          0.          0.          0.
       2.          0.          0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.1         0.
       0.          0.          0.          0.          0.          0.
       0.          3.25        0.          0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.1118034
       0.          0.          0.          0.          0.          0.
       0.          0.          4.5         0.          0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.12247449  0.          0.          0.          0.          0.
       0.          0.          0.          5.75        0.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.13228757  0.          0.          0.          0.
       0.          0.          0.          0.          7.          0.
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.          0.14142136  0.          0.          0.
       0.          0.          0.          0.          0.          8.25
       0.          0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.          0.          0.15        0.          0.
       0.          0.          0.          0.          0.          0.
       9.5         0.        ]
     [ 0.          0.          0.          0.          0.          0.
       0.          0.          0.          0.          0.          0.
       0.          0.          0.          0.          0.          0.
       0.         10.75      ]]


Here ``N`` is the number of Fock states included in the cavity mode.

.. _tensor-ptrace:

Partial trace
=============

The partial trace is an operation that reduces the dimension of a Hilbert space by eliminating some degrees of freedom by averaging (tracing). In this sense it is therefore the converse of the tensor product. It is useful when one is interested in only a part of a coupled quantum system.  For open quantum systems, this typically involves tracing over the environment leaving only the system of interest.  In QuTiP the class method  :func:`qutip.Qobj.ptrace` is used to take partial traces. :func:`qutip.Qobj.ptrace` acts on the :class:`qutip.Qobj` instance for which it is called, and it takes one argument ``sel``, which is a ``list`` of integers that mark the component systems that should be **kept**. All other components are traced out.

For example, the density matrix describing a single qubit obtained from a coupled two-qubit system is obtained via:

.. doctest:: [tensor]
  :options: +NORMALIZE_WHITESPACE

  >>> psi = tensor(basis(2, 0), basis(2, 1))

  >>> psi.ptrace(0)
  Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
  Qobj data =
  [[1. 0.]
   [0. 0.]]

  >>> psi.ptrace(1)
  Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
  Qobj data =
  [[0. 0.]
   [0. 1.]]

Note that the partial trace always results in a density matrix (mixed state), regardless of whether the composite system is a pure state (described by a state vector) or a mixed state (described by a density matrix):

.. doctest:: [tensor]
  :options: +NORMALIZE_WHITESPACE

  >>> psi = tensor((basis(2, 0) + basis(2, 1)).unit(), basis(2, 0))

  >>> psi
  Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
  Qobj data =
  [[0.70710678]
   [0.        ]
   [0.70710678]
   [0.        ]]

  >>> psi.ptrace(0)
  Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
  Qobj data =
  [[0.5 0.5]
   [0.5 0.5]]

  >>> rho = tensor(ket2dm((basis(2, 0) + basis(2, 1)).unit()), fock_dm(2, 0))

  >>> rho
  Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[0.5 0.  0.5 0. ]
   [0.  0.  0.  0. ]
   [0.5 0.  0.5 0. ]
   [0.  0.  0.  0. ]]

  >>> rho.ptrace(0)
  Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
  Qobj data =
  [[0.5 0.5]
   [0.5 0.5]]

Superoperators and Tensor Manipulations
=======================================

As described in :ref:`states-super`, *superoperators* are operators
that act on Liouville space, the vectorspace of linear operators.
Superoperators can be represented
using the isomorphism
:math:`\mathrm{vec} : \mathcal{L}(\mathcal{H}) \to \mathcal{H} \otimes \mathcal{H}` [Hav03]_, [Wat13]_.
To represent superoperators acting on :math:`\mathcal{L}(\mathcal{H}_1 \otimes \mathcal{H}_2)` thus takes some tensor rearrangement to get the desired ordering
:math:`\mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \mathcal{H}_1 \otimes \mathcal{H}_2`.

In particular, this means that :func:`qutip.tensor` does not act as
one might expect on the results of :func:`qutip.to_super`:

.. doctest:: [tensor]

  >>> A = qeye([2])

  >>> B = qeye([3])

  >>> to_super(tensor(A, B)).dims
  [[[2, 3], [2, 3]], [[2, 3], [2, 3]]]

  >>> tensor(to_super(A), to_super(B)).dims
  [[[2], [2], [3], [3]], [[2], [2], [3], [3]]]

In the former case, the result correctly has four copies
of the compound index with dims ``[2, 3]``. In the latter
case, however, each of the Hilbert space indices is listed
independently and in the wrong order.

The :func:`qutip.super_tensor` function performs the needed
rearrangement, providing the most direct analog to :func:`qutip.tensor` on
the underlying Hilbert space. In particular, for any two ``type="oper"``
Qobjs ``A`` and ``B``, ``to_super(tensor(A, B)) == super_tensor(to_super(A), to_super(B))`` and
``operator_to_vector(tensor(A, B)) == super_tensor(operator_to_vector(A), operator_to_vector(B))``. Returning to the previous example:

.. doctest:: [tensor]

  >>> super_tensor(to_super(A), to_super(B)).dims
  [[[2, 3], [2, 3]], [[2, 3], [2, 3]]]

The :func:`qutip.composite` function automatically switches between
:func:`qutip.tensor` and :func:`qutip.super_tensor` based on the ``type``
of its arguments, such that ``composite(A, B)`` returns an appropriate Qobj to
represent the composition of two systems.

.. doctest:: [tensor]

  >>> composite(A, B).dims
  [[2, 3], [2, 3]]

  >>> composite(to_super(A), to_super(B)).dims
  [[[2, 3], [2, 3]], [[2, 3], [2, 3]]]

QuTiP also allows more general tensor manipulations that are
useful for converting between superoperator representations [WBC11]_.
In particular, the :func:`tensor_contract` function allows for
contracting one or more pairs of indices. As detailed in
the `channel contraction tutorial`_, this can be used to find
superoperators that represent partial trace maps.
Using this functionality, we can construct some quite exotic maps,
such as a map from :math:`3 \times 3` operators to :math:`2 \times 2`
operators:

.. doctest:: [tensor]

  >>> tensor_contract(composite(to_super(A), to_super(B)), (1, 3), (4, 6)).dims
  [[[2], [2]], [[3], [3]]]



.. _channel contraction tutorial: https://nbviewer.ipython.org/github/qutip/qutip-notebooks/blob/master/examples/superop-contract.ipynb
.. _states:

*************************************
Manipulating States and Operators
*************************************

.. _states-intro:

Introduction
=================

In the previous guide section :ref:`basics`, we saw how to create states and operators, using the functions built into QuTiP. In this portion of the guide, we will look at performing basic operations with states and operators.  For more detailed demonstrations on how to use and manipulate these objects, see the examples on the `tutorials <https://qutip.org/tutorials.html>`_ web page.


.. _states-vectors:

State Vectors (kets or bras)
==============================

Here we begin by creating a Fock :func:`qutip.states.basis` vacuum state vector :math:`\left|0\right>` with in a Hilbert space with 5 number states, from 0 to 4:

.. testcode:: [states]

    vac = basis(5, 0)

    print(vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[1.]
     [0.]
     [0.]
     [0.]
     [0.]]




and then create a lowering operator :math:`\left(\hat{a}\right)` corresponding to 5 number states using the :func:`qutip.operators.destroy` function:

.. testcode:: [states]

    a = destroy(5)

    print(a)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = False
    Qobj data =
    [[0.         1.         0.         0.         0.        ]
     [0.         0.         1.41421356 0.         0.        ]
     [0.         0.         0.         1.73205081 0.        ]
     [0.         0.         0.         0.         2.        ]
     [0.         0.         0.         0.         0.        ]]


Now lets apply the destruction operator to our vacuum state ``vac``,


.. testcode:: [states]

    print(a * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [0.]
     [0.]
     [0.]]

We see that, as expected, the vacuum is transformed to the zero vector.  A more interesting example comes from using the adjoint of the lowering operator, the raising operator :math:`\hat{a}^\dagger`:

.. testcode:: [states]

    print(a.dag() * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
    [1.]
    [0.]
    [0.]
    [0.]]

The raising operator has in indeed raised the state `vec` from the vacuum to the :math:`\left| 1\right>` state.  Instead of using the dagger ``Qobj.dag()`` method to raise the state, we could have also used the built in :func:`qutip.operators.create` function to make a raising operator:

.. testcode:: [states]

    c = create(5)

    print(c * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [1.]
     [0.]
     [0.]
     [0.]]

which does the same thing.  We can raise the vacuum state more than once by successively apply the raising operator:

.. testcode:: [states]

    print(c * c * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.        ]
     [0.        ]
     [1.41421356]
     [0.        ]
     [0.        ]]

or just taking the square of the raising operator :math:`\left(\hat{a}^\dagger\right)^{2}`:

.. testcode:: [states]

    print(c ** 2 * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.        ]
     [0.        ]
     [1.41421356]
     [0.        ]
     [0.        ]]

Applying the raising operator twice gives the expected :math:`\sqrt{n + 1}` dependence.  We can use the product of :math:`c * a` to also apply the number operator to the state vector ``vac``:

.. testcode:: [states]

    print(c * a * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [0.]
     [0.]
     [0.]]

or on the :math:`\left| 1\right>` state:

.. testcode:: [states]

    print(c * a * (c * vac))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [1.]
     [0.]
     [0.]
     [0.]]

or the :math:`\left| 2\right>` state:

.. testcode:: [states]

    print(c * a * (c**2 * vac))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.        ]
     [0.        ]
     [2.82842712]
     [0.        ]
     [0.        ]]

Notice how in this last example, application of the number operator does not give the expected value :math:`n=2`, but rather :math:`2\sqrt{2}`.  This is because this last state is not normalized to unity as :math:`c\left| n\right> = \sqrt{n+1}\left| n+1\right>`.  Therefore, we should normalize our vector first:

.. testcode:: [states]

    print(c * a * (c**2 * vac).unit())

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [2.]
     [0.]
     [0.]]

Since we are giving a demonstration of using states and operators, we have done a lot more work than we should have.  For example, we do not need to operate on the vacuum state to generate a higher number Fock state.  Instead we can use the :func:`qutip.states.basis` (or :func:`qutip.states.fock`) function to directly obtain the required state:

.. testcode:: [states]

    ket = basis(5, 2)

    print(ket)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [1.]
     [0.]
     [0.]]

Notice how it is automatically normalized.  We can also use the built in :func:`qutip.operators.num` operator:

.. testcode:: [states]

    n = num(5)

    print(n)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0. 0. 0. 0. 0.]
     [0. 1. 0. 0. 0.]
     [0. 0. 2. 0. 0.]
     [0. 0. 0. 3. 0.]
     [0. 0. 0. 0. 4.]]

Therefore, instead of ``c * a * (c ** 2 * vac).unit()`` we have:

.. testcode:: [states]

    print(n * ket)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [2.]
     [0.]
     [0.]]

We can also create superpositions of states:

.. testcode:: [states]

    ket = (basis(5, 0) + basis(5, 1)).unit()

    print(ket)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.70710678]
     [0.70710678]
     [0.        ]
     [0.        ]
     [0.        ]]

where we have used the :func:`qutip.Qobj.unit` method to again normalize the state. Operating with the number function again:

.. testcode:: [states]

    print(n * ket)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.        ]
     [0.70710678]
     [0.        ]
     [0.        ]
     [0.        ]]

We can also create coherent states and squeezed states by applying the :func:`qutip.operators.displace` and :func:`qutip.operators.squeeze` functions to the vacuum state:

.. testcode:: [states]

    vac = basis(5, 0)

    d = displace(5, 1j)

    s = squeeze(5, np.complex(0.25, 0.25))

    print(d * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[ 0.60655682+0.j        ]
     [ 0.        +0.60628133j]
     [-0.4303874 +0.j        ]
     [ 0.        -0.24104351j]
     [ 0.14552147+0.j        ]]

.. testcode:: [states]

    print(d * s * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[ 0.65893786+0.08139381j]
     [ 0.10779462+0.51579735j]
     [-0.37567217-0.01326853j]
     [-0.02688063-0.23828775j]
     [ 0.26352814+0.11512178j]]

Of course, displacing the vacuum gives a coherent state, which can also be generated using the built in :func:`qutip.states.coherent` function.


.. _states-dm:

Density matrices
=================

One of the main purpose of QuTiP is to explore the dynamics of **open** quantum systems, where the most general state of a system is no longer a state vector, but rather a density matrix.  Since operations on density matrices operate identically to those of vectors, we will just briefly highlight creating and using these structures.

The simplest density matrix is created by forming the outer-product :math:`\left|\psi\right>\left<\psi\right|` of a ket vector:

.. testcode:: [states]

    ket = basis(5, 2)

    print(ket * ket.dag())

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 1. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]]

A similar task can also be accomplished via the :func:`qutip.states.fock_dm` or :func:`qutip.states.ket2dm` functions:

.. testcode:: [states]

    print(fock_dm(5, 2))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 1. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]]

.. testcode:: [states]

    print(ket2dm(ket))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 1. 0. 0.]
     [0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0.]]

If we want to create a density matrix with equal classical probability of being found in the :math:`\left|2\right>` or :math:`\left|4\right>` number states we can do the following:

.. testcode:: [states]

    print(0.5 * ket2dm(basis(5, 4)) + 0.5 * ket2dm(basis(5, 2)))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0.  0.  0.  0.  0. ]
     [0.  0.  0.  0.  0. ]
     [0.  0.  0.5 0.  0. ]
     [0.  0.  0.  0.  0. ]
     [0.  0.  0.  0.  0.5]]

or use ``0.5 * fock_dm(5, 2) + 0.5 * fock_dm(5, 4)``. There are also several other built-in functions for creating predefined density matrices, for example :func:`qutip.states.coherent_dm` and :func:`qutip.states.thermal_dm` which create coherent state and thermal state density matrices, respectively.


.. testcode:: [states]

    print(coherent_dm(5, 1.25))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0.20980701 0.26141096 0.23509686 0.15572585 0.13390765]
     [0.26141096 0.32570738 0.29292109 0.19402805 0.16684347]
     [0.23509686 0.29292109 0.26343512 0.17449684 0.1500487 ]
     [0.15572585 0.19402805 0.17449684 0.11558499 0.09939079]
     [0.13390765 0.16684347 0.1500487  0.09939079 0.0854655 ]]

.. testcode:: [states]

    print(thermal_dm(5, 1.25))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0.46927974 0.         0.         0.         0.        ]
     [0.         0.26071096 0.         0.         0.        ]
     [0.         0.         0.14483942 0.         0.        ]
     [0.         0.         0.         0.08046635 0.        ]
     [0.         0.         0.         0.         0.04470353]]

QuTiP also provides a set of distance metrics for determining how close two density matrix distributions are to each other. Included are the trace distance :func:`qutip.metrics.tracedist`, fidelity :func:`qutip.metrics.fidelity`, Hilbert-Schmidt distance :func:`qutip.metrics.hilbert_dist`, Bures distance :func:`qutip.metrics.bures_dist`, Bures angle :func:`qutip.metrics.bures_angle`, and quantum Hellinger distance :func:`qutip.metrics.hellinger_dist`.

.. testcode:: [states]

    x = coherent_dm(5, 1.25)

    y = coherent_dm(5, np.complex(0, 1.25))  # <-- note the 'j'

    z = thermal_dm(5, 0.125)

    np.testing.assert_almost_equal(fidelity(x, x), 1)

    np.testing.assert_almost_equal(hellinger_dist(x, y), 1.3819080728932833)

We also know that for two pure states, the trace distance (T) and the fidelity (F) are related by :math:`T = \sqrt{1 - F^{2}}`, while the quantum Hellinger distance (QHE) between two pure states :math:`\left|\psi\right>` and :math:`\left|\phi\right>` is given by :math:`QHE = \sqrt{2 - 2\left|\left<\psi | \phi\right>\right|^2}`.

.. testcode:: [states]

    np.testing.assert_almost_equal(tracedist(y, x), np.sqrt(1 - fidelity(y, x) ** 2))

For a pure state and a mixed state, :math:`1 - F^{2} \le T` which can also be verified:

.. testcode:: [states]

    assert 1 - fidelity(x, z) ** 2 < tracedist(x, z)

.. _states-qubit:

Qubit (two-level) systems
=========================

Having spent a fair amount of time on basis states that represent harmonic oscillator states, we now move on to qubit, or two-level quantum systems (for example a spin-1/2). To create a state vector corresponding to a qubit system, we use the same :func:`qutip.states.basis`, or :func:`qutip.states.fock`, function with only two levels:


.. testcode:: [states]

    spin = basis(2, 0)

Now at this point one may ask how this state is different than that of a harmonic oscillator in the vacuum state truncated to two energy levels?

.. testcode:: [states]

    vac = basis(2, 0)

At this stage, there is no difference.  This should not be surprising as we called the exact same function twice.  The difference between the two comes from the action of the spin operators :func:`qutip.operators.sigmax`, :func:`qutip.operators.sigmay`, :func:`qutip.operators.sigmaz`, :func:`qutip.operators.sigmap`, and :func:`qutip.operators.sigmam` on these two-level states.  For example, if ``vac`` corresponds to the vacuum state of a harmonic oscillator, then, as we have already seen, we can use the raising operator to get the :math:`\left|1\right>` state:

.. testcode:: [states]

    print(vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[1.]
     [0.]]

.. testcode:: [states]

    c = create(2)

    print(c * vac)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[0.]
     [1.]]

For a spin system, the operator analogous to the raising operator is the sigma-plus operator :func:`qutip.operators.sigmap`.  Operating on the ``spin`` state gives:

.. testcode:: [states]

    print(spin)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[1.]
     [0.]]

.. testcode:: [states]

    print(sigmap() * spin)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[0.]
     [0.]]

Now we see the difference!  The :func:`qutip.operators.sigmap` operator acting on the ``spin`` state returns the zero vector.  Why is this?  To see what happened, let us use the :func:`qutip.operators.sigmaz` operator:

.. testcode:: [states]

    print(sigmaz())

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[ 1.  0.]
     [ 0. -1.]]

.. testcode:: [states]

    print(sigmaz() * spin)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[1.]
     [0.]]

.. testcode:: [states]

    spin2 = basis(2, 1)

    print(spin2)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
    Qobj data =
    [[0.]
     [1.]]

.. testcode:: [states]

    print(sigmaz() * spin2)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
  Qobj data =
  [[ 0.]
   [-1.]]

The answer is now apparent.  Since the QuTiP :func:`qutip.operators.sigmaz` function uses the standard z-basis representation of the sigma-z spin operator, the ``spin`` state corresponds to the :math:`\left|\uparrow\right>` state of a two-level spin system while ``spin2`` gives the :math:`\left|\downarrow\right>` state.  Therefore, in our previous example ``sigmap() * spin``, we raised the qubit state out of the truncated two-level Hilbert space resulting in the zero state.

While at first glance this convention might seem somewhat odd, it is in fact quite handy. For one, the spin operators remain in the conventional form. Second, when the spin system is in the :math:`\left|\uparrow\right>` state:

.. testcode:: [states]

    print(sigmaz() * spin)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
  Qobj data =
  [[1.]
   [0.]]

the non-zero component is the zeroth-element of the underlying matrix (remember that python uses c-indexing, and matrices start with the zeroth element).  The :math:`\left|\downarrow\right>` state therefore has a non-zero entry in the first index position. This corresponds nicely with the quantum information definitions of qubit states, where the excited :math:`\left|\uparrow\right>` state is label as :math:`\left|0\right>`, and the :math:`\left|\downarrow\right>` state by :math:`\left|1\right>`.

If one wants to create spin operators for higher spin systems, then the :func:`qutip.operators.jmat` function comes in handy.

.. _states-expect:

Expectation values
===================

Some of the most important information about quantum systems comes from calculating the expectation value of operators, both Hermitian and non-Hermitian, as the state or density matrix of the system varies in time.  Therefore, in this section we demonstrate the use of the :func:`qutip.expect` function.  To begin:

.. testcode:: [states]

    vac = basis(5, 0)

    one = basis(5, 1)

    c = create(5)

    N = num(5)

    np.testing.assert_almost_equal(expect(N, vac), 0)

    np.testing.assert_almost_equal(expect(N, one), 1)

    coh = coherent_dm(5, 1.0j)

    np.testing.assert_almost_equal(expect(N, coh), 0.9970555745806597)

    cat = (basis(5, 4) + 1.0j * basis(5, 3)).unit()

    np.testing.assert_almost_equal(expect(c, cat), 0.9999999999999998j)


The :func:`qutip.expect` function also accepts lists or arrays of state vectors or density matrices for the second input:

.. testcode:: [states]

    states = [(c**k * vac).unit() for k in range(5)]  # must normalize

    print(expect(N, states))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [0. 1. 2. 3. 4.]

.. testcode:: [states]

    cat_list = [(basis(5, 4) + x * basis(5, 3)).unit() for x in [0, 1.0j, -1.0, -1.0j]]

    print(expect(c, cat_list))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [ 0.+0.j  0.+1.j -1.+0.j  0.-1.j]

Notice how in this last example, all of the return values are complex numbers.  This is because the :func:`qutip.expect` function looks to see whether the operator is Hermitian or not.  If the operator is Hermitian, then the output will always be real.  In the case of non-Hermitian operators, the return values may be complex.  Therefore, the :func:`qutip.expect` function will return an array of complex values for non-Hermitian operators when the input is a list/array of states or density matrices.

Of course, the :func:`qutip.expect` function works for spin states and operators:


.. testcode:: [states]

    up = basis(2, 0)

    down = basis(2, 1)

    np.testing.assert_almost_equal(expect(sigmaz(), up), 1)

    np.testing.assert_almost_equal(expect(sigmaz(), down), -1)


as well as the composite objects discussed in the next section :ref:`tensor`:

.. testcode:: [states]

    spin1 = basis(2, 0)

    spin2 = basis(2, 1)

    two_spins = tensor(spin1, spin2)

    sz1 = tensor(sigmaz(), qeye(2))

    sz2 = tensor(qeye(2), sigmaz())

    np.testing.assert_almost_equal(expect(sz1, two_spins), 1)

    np.testing.assert_almost_equal(expect(sz2, two_spins), -1)

.. _states-super:

Superoperators and Vectorized Operators
=======================================

In addition to state vectors and density operators, QuTiP allows for
representing maps that act linearly on density operators using the Kraus,
Liouville supermatrix and Choi matrix formalisms. This support is based on the
correspondence between linear operators acting on a Hilbert space, and vectors
in two copies of that Hilbert space,
:math:`\mathrm{vec} : \mathcal{L}(\mathcal{H}) \to \mathcal{H} \otimes \mathcal{H}`
[Hav03]_, [Wat13]_.

This isomorphism is implemented in QuTiP by the
:obj:`~qutip.superoperator.operator_to_vector` and
:obj:`~qutip.superoperator.vector_to_operator` functions:

.. testcode:: [states]

    psi = basis(2, 0)

    rho = ket2dm(psi)

    print(rho)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[1. 0.]
     [0. 0.]]

.. testcode:: [states]

    vec_rho = operator_to_vector(rho)

    print(vec_rho)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[[2], [2]], [1]], shape = (4, 1), type = operator-ket
  Qobj data =
  [[1.]
   [0.]
   [0.]
   [0.]]

.. testcode:: [states]

    rho2 = vector_to_operator(vec_rho)

    np.testing.assert_almost_equal((rho - rho2).norm(), 0)

The :attr:`~qutip.Qobj.type` attribute indicates whether a quantum object is
a vector corresponding to an operator (``operator-ket``), or its Hermitian
conjugate (``operator-bra``).

Note that QuTiP uses the *column-stacking* convention for the isomorphism
between :math:`\mathcal{L}(\mathcal{H})` and :math:`\mathcal{H} \otimes \mathcal{H}`:

.. testcode:: [states]

    A = Qobj(np.arange(4).reshape((2, 2)))

    print(A)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = False
    Qobj data =
    [[0. 1.]
     [2. 3.]]

.. testcode:: [states]

    print(operator_to_vector(A))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [1]], shape = (4, 1), type = operator-ket
    Qobj data =
    [[0.]
     [2.]
     [1.]
     [3.]]

Since :math:`\mathcal{H} \otimes \mathcal{H}` is a vector space, linear maps
on this space can be represented as matrices, often called *superoperators*.
Using the :obj:`~qutip.Qobj`, the :obj:`~qutip.superoperator.spre` and :obj:`~qutip.superoperator.spost` functions, supermatrices
corresponding to left- and right-multiplication respectively can be quickly
constructed.

.. testcode:: [states]

    X = sigmax()

    S = spre(X) * spost(X.dag()) # Represents conjugation by X.

Note that this is done automatically by the :obj:`~qutip.superop_reps.to_super` function when given
``type='oper'`` input.

.. testcode:: [states]

    S2 = to_super(X)

    np.testing.assert_almost_equal((S - S2).norm(), 0)

Quantum objects representing superoperators are denoted by ``type='super'``:

.. testcode:: [states]

  print(S)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True
  Qobj data =
  [[0. 0. 0. 1.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]
   [1. 0. 0. 0.]]

Information about superoperators, such as whether they represent completely
positive maps, is exposed through the :attr:`~qutip.Qobj.iscp`, :attr:`~qutip.Qobj.istp`
and :attr:`~qutip.Qobj.iscptp` attributes:

.. testcode:: [states]

    print(S.iscp, S.istp, S.iscptp)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    True True True

In addition, dynamical generators on this extended space, often called
*Liouvillian superoperators*, can be created using the :func:`~qutip.superoperator.liouvillian` function. Each of these takes a Hamiltonian along with
a list of collapse operators, and returns a ``type="super"`` object that can
be exponentiated to find the superoperator for that evolution.

.. testcode:: [states]

    H = 10 * sigmaz()

    c1 = destroy(2)

    L = liouvillian(H, [c1])

    print(L)

    S = (12 * L).expm()

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = False
    Qobj data =
    [[ 0.  +0.j  0.  +0.j  0.  +0.j  1.  +0.j]
     [ 0.  +0.j -0.5+20.j  0.  +0.j  0.  +0.j]
     [ 0.  +0.j  0.  +0.j -0.5-20.j  0.  +0.j]
     [ 0.  +0.j  0.  +0.j  0.  +0.j -1.  +0.j]]

For qubits, a particularly useful way to visualize superoperators is to plot them in the Pauli basis,
such that :math:`S_{\mu,\nu} = \langle\!\langle \sigma_{\mu} | S[\sigma_{\nu}] \rangle\!\rangle`. Because
the Pauli basis is Hermitian, :math:`S_{\mu,\nu}` is a real number for all Hermitian-preserving superoperators
:math:`S`,
allowing us to plot the elements of :math:`S` as a `Hinton diagram <https://matplotlib.org/examples/specialty_plots/hinton_demo.html>`_. In such diagrams, positive elements are indicated by white squares, and negative elements
by black squares. The size of each element is indicated by the size of the corresponding square. For instance,
let :math:`S[\rho] = \sigma_x \rho \sigma_x^{\dagger}`. Then :math:`S[\sigma_{\mu}] = \sigma_{\mu} \cdot \begin{cases} +1 & \mu = 0, x \\ -1 & \mu = y, z \end{cases}`. We can quickly see this by noting that the :math:`Y` and :math:`Z` elements
of the Hinton diagram for :math:`S` are negative:

.. plot::

    from qutip import *
    settings.colorblind_safe = True

    import matplotlib.pyplot as plt
    plt.rcParams['savefig.transparent'] = True

    X = sigmax()
    S = spre(X) * spost(X.dag())

    hinton(S)

Choi, Kraus, Stinespring and :math:`\chi` Representations
=========================================================

In addition to the superoperator representation of quantum maps, QuTiP
supports several other useful representations. First, the Choi matrix
:math:`J(\Lambda)` of a quantum map :math:`\Lambda` is useful for working with
ancilla-assisted process tomography (AAPT), and for reasoning about properties
of a map or channel. Up to normalization, the Choi matrix is defined by acting
:math:`\Lambda` on half of an entangled pair. In the column-stacking
convention,

.. math::

    J(\Lambda) = (\mathbb{1} \otimes \Lambda) [|\mathbb{1}\rangle\!\rangle \langle\!\langle \mathbb{1}|].

In QuTiP, :math:`J(\Lambda)` can be found by calling the :func:`~qutip.superop_reps.to_choi`
function on a ``type="super"`` :obj:`~Qobj`.

.. testcode:: [states]

    X = sigmax()

    S = sprepost(X, X)

    J = to_choi(S)

    print(J)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True, superrep = choi
    Qobj data =
    [[0. 0. 0. 0.]
     [0. 1. 1. 0.]
     [0. 1. 1. 0.]
     [0. 0. 0. 0.]]

.. testcode:: [states]

  print(to_choi(spre(qeye(2))))

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True, superrep = choi
  Qobj data =
  [[1. 0. 0. 1.]
   [0. 0. 0. 0.]
   [0. 0. 0. 0.]
   [1. 0. 0. 1.]]

If a :obj:`~Qobj` instance is already in the Choi :attr:`~Qobj.superrep`, then calling :func:`~qutip.superop_reps.to_choi`
does nothing:

.. testcode:: [states]

    print(to_choi(J))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True, superrep = choi
    Qobj data =
    [[0. 0. 0. 0.]
     [0. 1. 1. 0.]
     [0. 1. 1. 0.]
     [0. 0. 0. 0.]]

To get back to the superoperator representation, simply use the :func:`~qutip.superop_reps.to_super` function.
As with :func:`~qutip.superop_reps.to_choi`, :func:`~qutip.superop_reps.to_super` is idempotent:

.. testcode:: [states]

    print(to_super(J) - S)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True
    Qobj data =
    [[0. 0. 0. 0.]
     [0. 0. 0. 0.]
     [0. 0. 0. 0.]
     [0. 0. 0. 0.]]

.. testcode:: [states]

    print(to_super(S))

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True
    Qobj data =
    [[0. 0. 0. 1.]
     [0. 0. 1. 0.]
     [0. 1. 0. 0.]
     [1. 0. 0. 0.]]

We can quickly obtain another useful representation from the Choi matrix by taking its eigendecomposition.
In particular, let :math:`\{A_i\}` be a set of operators such that
:math:`J(\Lambda) = \sum_i |A_i\rangle\!\rangle \langle\!\langle A_i|`.
We can write :math:`J(\Lambda)` in this way
for any hermicity-preserving map; that is, for any map :math:`\Lambda` such that :math:`J(\Lambda) = J^\dagger(\Lambda)`.
These operators then form the Kraus representation of :math:`\Lambda`. In particular, for any input :math:`\rho`,

.. math::

    \Lambda(\rho) = \sum_i A_i \rho A_i^\dagger.

Notice using the column-stacking identity that :math:`(C^\mathrm{T} \otimes A) |B\rangle\!\rangle = |ABC\rangle\!\rangle`,
we have that

.. math::

      \sum_i (\mathbb{1} \otimes A_i) (\mathbb{1} \otimes A_i)^\dagger |\mathbb{1}\rangle\!\rangle \langle\!\langle\mathbb{1}|
    = \sum_i |A_i\rangle\!\rangle \langle\!\langle A_i| = J(\Lambda).

The Kraus representation of a hermicity-preserving map can be found in QuTiP
using the :func:`~qutip.superop_reps.to_kraus` function.

.. testcode:: [states]

    del sum # np.sum overwrote sum and caused a bug.


.. testcode:: [states]

    I, X, Y, Z = qeye(2), sigmax(), sigmay(), sigmaz()

.. testcode:: [states]

    S = sum([sprepost(P, P) for P in (I, X, Y, Z)]) / 4
    print(S)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True
    Qobj data =
    [[0.5 0.  0.  0.5]
     [0.  0.  0.  0. ]
     [0.  0.  0.  0. ]
     [0.5 0.  0.  0.5]]

.. testcode:: [states]

    J = to_choi(S)
    print(J)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True, superrep = choi
    Qobj data =
    [[0.5 0.  0.  0. ]
     [0.  0.5 0.  0. ]
     [0.  0.  0.5 0. ]
     [0.  0.  0.  0.5]]

.. testcode:: [states]

    print(J.eigenstates()[1])

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [Quantum object: dims = [[[2], [2]], [1, 1]], shape = (4, 1), type = operator-ket
    Qobj data =
    [[1.]
     [0.]
     [0.]
     [0.]]
     Quantum object: dims = [[[2], [2]], [1, 1]], shape = (4, 1), type = operator-ket
    Qobj data =
    [[0.]
     [1.]
     [0.]
     [0.]]
     Quantum object: dims = [[[2], [2]], [1, 1]], shape = (4, 1), type = operator-ket
    Qobj data =
    [[0.]
     [0.]
     [1.]
     [0.]]
     Quantum object: dims = [[[2], [2]], [1, 1]], shape = (4, 1), type = operator-ket
    Qobj data =
    [[0.]
     [0.]
     [0.]
     [1.]]]

.. testcode:: [states]

    K = to_kraus(S)
    print(K)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0.70710678 0.        ]
     [0.         0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = False
    Qobj data =
    [[0.         0.        ]
     [0.70710678 0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = False
    Qobj data =
    [[0.         0.70710678]
     [0.         0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0.         0.        ]
     [0.         0.70710678]]]

As with the other representation conversion functions, :func:`~qutip.superop_reps.to_kraus`
checks the :attr:`~Qobj.superrep` attribute of its input, and chooses an appropriate
conversion method. Thus, in the above example, we can also call :func:`~qutip.superop_reps.to_kraus`
on ``J``.

.. testcode:: [states]

    KJ = to_kraus(J)
    print(KJ)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0.70710678 0.        ]
     [0.         0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = False
    Qobj data =
    [[0.         0.        ]
     [0.70710678 0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = False
    Qobj data =
    [[0.         0.70710678]
     [0.         0.        ]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0.         0.        ]
     [0.         0.70710678]]]

.. testcode:: [states]

    for A, AJ in zip(K, KJ):
      print(A - AJ)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0. 0.]
     [0. 0.]]
    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0. 0.]
     [0. 0.]]
    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0. 0.]
     [0. 0.]]
    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0. 0.]
     [0. 0.]]

The Stinespring representation is closely related to the Kraus representation,
and consists of a pair of operators :math:`A` and :math:`B` such that for
all operators :math:`X` acting on :math:`\mathcal{H}`,

.. math::

    \Lambda(X) = \operatorname{Tr}_2(A X B^\dagger),

where the partial trace is over a new index that corresponds to the
index in the Kraus summation. Conversion to Stinespring
is handled by the :func:`~qutip.superop_reps.to_stinespring`
function.

.. testcode:: [states]

    a = create(2).dag()

    S_ad = sprepost(a * a.dag(), a * a.dag()) + sprepost(a, a.dag())
    S = 0.9 * sprepost(I, I) + 0.1 * S_ad

    print(S)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = False
    Qobj data =
    [[1.  0.  0.  0.1]
     [0.  0.9 0.  0. ]
     [0.  0.  0.9 0. ]
     [0.  0.  0.  0.9]]

.. testcode:: [states]

    A, B = to_stinespring(S)
    print(A)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 3], [2]], shape = (6, 2), type = oper, isherm = False
    Qobj data =
    [[-0.98845443  0.        ]
     [ 0.          0.31622777]
     [ 0.15151842  0.        ]
     [ 0.         -0.93506452]
     [ 0.          0.        ]
     [ 0.         -0.16016975]]

.. testcode:: [states]

    print(B)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 3], [2]], shape = (6, 2), type = oper, isherm = False
    Qobj data =
    [[-0.98845443  0.        ]
     [ 0.          0.31622777]
     [ 0.15151842  0.        ]
     [ 0.         -0.93506452]
     [ 0.          0.        ]
     [ 0.         -0.16016975]]

Notice that a new index has been added, such that :math:`A` and :math:`B`
have dimensions ``[[2, 3], [2]]``, with the length-3 index representing the
fact that the Choi matrix is rank-3 (alternatively, that the map has three
Kraus operators).

.. testcode:: [states]

    to_kraus(S)
    print(to_choi(S).eigenenergies())

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [0.         0.04861218 0.1        1.85138782]

Finally, the last superoperator representation supported by QuTiP is
the :math:`\chi`-matrix representation,

.. math::

    \Lambda(\rho) = \sum_{\alpha,\beta} \chi_{\alpha,\beta} B_{\alpha} \rho B_{\beta}^\dagger,

where :math:`\{B_\alpha\}` is a basis for the space of matrices acting
on :math:`\mathcal{H}`. In QuTiP, this basis is taken to be the Pauli
basis :math:`B_\alpha = \sigma_\alpha / \sqrt{2}`. Conversion to the
:math:`\chi` formalism is handled by the :func:`~qutip.superop_reps.to_chi`
function.

.. testcode:: [states]

    chi = to_chi(S)
    print(chi)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = True, superrep = chi
    Qobj data =
    [[3.7+0.j  0. +0.j  0. +0.j  0.1+0.j ]
     [0. +0.j  0.1+0.j  0. +0.1j 0. +0.j ]
     [0. +0.j  0. -0.1j 0.1+0.j  0. +0.j ]
     [0.1+0.j  0. +0.j  0. +0.j  0.1+0.j ]]


One convenient property of the :math:`\chi` matrix is that the average
gate fidelity with the identity map can be read off directly from
the :math:`\chi_{00}` element:

.. testcode:: [states]

    np.testing.assert_almost_equal(average_gate_fidelity(S), 0.9499999999999998)

    print(chi[0, 0] / 4)

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    (0.925+0j)

Here, the factor of 4 comes from the dimension of the underlying
Hilbert space :math:`\mathcal{H}`. As with the superoperator
and Choi representations, the :math:`\chi` representation is
denoted by the :attr:`~Qobj.superrep`, such that :func:`~qutip.superop_reps.to_super`,
:func:`~qutip.superop_reps.to_choi`, :func:`~qutip.superop_reps.to_kraus`,
:func:`~qutip.superop_reps.to_stinespring` and :func:`~qutip.superop_reps.to_chi`
all convert from the :math:`\chi` representation appropriately.

Properties of Quantum Maps
==========================

In addition to converting between the different representations of quantum maps,
QuTiP also provides attributes to make it easy to check if a map is completely
positive, trace preserving and/or hermicity preserving. Each of these attributes
uses :attr:`~Qobj.superrep` to automatically perform any needed conversions.

In particular, a quantum map is said to be positive (but not necessarily completely
positive) if it maps all positive operators to positive operators. For instance, the
transpose map :math:`\Lambda(\rho) = \rho^{\mathrm{T}}` is a positive map. We run into
problems, however, if we tensor :math:`\Lambda` with the identity to get a partial
transpose map.

.. testcode:: [states]

    rho = ket2dm(bell_state())
    rho_out = partial_transpose(rho, [0, 1])
    print(rho_out.eigenenergies())

**Output**:

.. testoutput:: [states]
    :options: +NORMALIZE_WHITESPACE

    [-0.5  0.5  0.5  0.5]

Notice that even though we started with a positive map, we got an operator out
with negative eigenvalues. Complete positivity addresses this by requiring that
a map returns positive operators for all positive operators, and does so even
under tensoring with another map. The Choi matrix is very useful here, as it
can be shown that a map is completely positive if and only if its Choi matrix
is positive [Wat13]_. QuTiP implements this check with the :attr:`~Qobj.iscp`
attribute. As an example, notice that the snippet above already calculates
the Choi matrix of the transpose map by acting it on half of an entangled
pair. We simply need to manually set the ``dims`` and ``superrep`` attributes to reflect the
structure of the underlying Hilbert space and the chosen representation.

.. testcode:: [states]

    J = rho_out
    J.dims = [[[2], [2]], [[2], [2]]]
    J.superrep = 'choi'
    print(J.iscp)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  False

This confirms that the transpose map is not completely positive. On the other hand,
the transpose map does satisfy a weaker condition, namely that it is hermicity preserving.
That is, :math:`\Lambda(\rho) = (\Lambda(\rho))^\dagger` for all :math:`\rho` such that
:math:`\rho = \rho^\dagger`. To see this, we note that :math:`(\rho^{\mathrm{T}})^\dagger
= \rho^*`, the complex conjugate of :math:`\rho`. By assumption, :math:`\rho = \rho^\dagger
= (\rho^*)^{\mathrm{T}}`, though, such that :math:`\Lambda(\rho) = \Lambda(\rho^\dagger) = \rho^*`.
We can confirm this by checking the :attr:`~Qobj.ishp` attribute:

.. testcode:: [states]

    print(J.ishp)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  True

Next, we note that the transpose map does preserve the trace of its inputs, such that
:math:`\operatorname{Tr}(\Lambda[\rho]) = \operatorname{Tr}(\rho)` for all :math:`\rho`.
This can be confirmed by the :attr:`~Qobj.istp` attribute:

.. testcode:: [states]

    print(J.istp)

**Output**:

.. testoutput:: [states]
  :options: +NORMALIZE_WHITESPACE

  False

Finally, a map is called a quantum channel if it always maps valid states to valid
states. Formally, a map is a channel if it is both completely positive and trace preserving.
Thus, QuTiP provides a single attribute to quickly check that this is true.

.. doctest:: [states]

    >>> print(J.iscptp)
    False

    >>> print(to_super(qeye(2)).iscptp)
    True
.. _basics:

************************************
Basic Operations on Quantum Objects
************************************

.. _basics-first:

First things first
==================

.. warning:: Do not run QuTiP from the installation directory.

To load the qutip modules, we must first call the import statement:

.. code-block:: Python

   from qutip import *

that will load all of the user available functions. Often, we also need to import the NumPy and Matplotlib libraries with:

.. code-block:: Python

   import numpy as np

   import matplotlib.pyplot as plt

Note that, in the rest of the documentation, functions are written using `qutip.module.function()` notation which links to the corresponding function in the QuTiP API: :ref:`functions`. However, in calling `import *`, we have already loaded all of the QuTiP modules. Therefore, we will only need the function name and not the complete path when calling the function from the interpreter prompt, Python script, or Jupyter notebook.

.. _basics-qobj:

The quantum object class
========================

.. _basics-qobj-intro:

Introduction
---------------

The key difference between classical and quantum mechanics lies in the use of operators instead of numbers as variables. Moreover, we need to specify state vectors and their properties. Therefore, in computing the dynamics of quantum systems we need a data structure that is capable of encapsulating the properties of a quantum operator and ket/bra vectors. The quantum object class, :func:`qutip.Qobj`, accomplishes this using matrix representation.

To begin, let us create a blank ``Qobj``:

.. testcode:: [basics]

    print(Qobj())

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
    Qobj data =
    [[0.]]


where we see the blank ``Qobj`` object with dimensions, shape, and data. Here the data corresponds to a 1x1-dimensional matrix consisting of a single zero entry.

.. Hint:: By convention, Class objects in Python such as ``Qobj()`` differ from functions in the use of a beginning capital letter.

We can create a ``Qobj`` with a user defined data set by passing a list or array of data into the ``Qobj``:

.. testcode:: [basics]

    print(Qobj([[1],[2],[3],[4],[5]]))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[1.]
    [2.]
    [3.]
    [4.]
    [5.]]

.. testcode:: [basics]

    x = np.array([[1, 2, 3, 4, 5]])
    print(Qobj(x))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[1], [5]], shape = (1, 5), type = bra
    Qobj data =
    [[1. 2. 3. 4. 5.]]

.. testcode:: [basics]
    :hide:

    np.random.seed(42)

.. testcode:: [basics]

    r = np.random.rand(4, 4)
    print(Qobj(r))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[4], [4]], shape = (4, 4), type = oper, isherm = False
    Qobj data =
    [[0.37454012 0.95071431 0.73199394 0.59865848]
     [0.15601864 0.15599452 0.05808361 0.86617615]
     [0.60111501 0.70807258 0.02058449 0.96990985]
     [0.83244264 0.21233911 0.18182497 0.18340451]]

Notice how both the dims and shape change according to the input data.  Although dims and shape appear to have the same function, the difference will become quite clear in the section on :ref:`tensor products and partial traces <tensor>`.

.. note:: If you are running QuTiP from a python script you must use the :func:`print` function to view the Qobj attributes.

.. _basics-qobj-states:

States and operators
---------------------

Manually specifying the data for each quantum object is inefficient. Even more so when most objects correspond to commonly used types such as the ladder operators of a harmonic oscillator, the Pauli spin operators for a two-level system, or state vectors such as Fock states. Therefore, QuTiP includes predefined objects for a variety of states:

.. cssclass:: table-striped

+--------------------------+----------------------------------+----------------------------------------+
| States                   | Command (# means optional)       | Inputs                                 |
+==========================+==================================+========================================+
| Fock state ket vector    | ``basis(N,#m)``/``fock(N,#m)``   | N = number of levels in Hilbert space, |
|                          |                                  | m = level containing excitation        |
|                          |                                  | (0 if no m given)                      |
+--------------------------+----------------------------------+----------------------------------------+
| Fock density matrix      | ``fock_dm(N,#p)``                | same as basis(N,m) / fock(N,m)         |
| (outer product of basis) |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+
| Coherent state           | ``coherent(N,alpha)``            | alpha = complex number (eigenvalue)    |
|                          |                                  | for requested coherent state           |
+--------------------------+----------------------------------+----------------------------------------+
| Coherent density matrix  | ``coherent_dm(N,alpha)``         | same as coherent(N,alpha)              |
| (outer product)          |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+
| Thermal density matrix   | ``thermal_dm(N,n)``              | n = particle number expectation value  |
| (for n particles)        |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+

and operators:

.. cssclass:: table-striped

+--------------------------+----------------------------+----------------------------------------+
| Operators                | Command (# means optional) | Inputs                                 |
+==========================+============================+========================================+
| Charge operator          | ``charge(N,M=-N)``         | Diagonal operator with entries         |
|                          |                            | from M..0..N.                          |
+--------------------------+----------------------------+----------------------------------------+
| Commutator               | ``commutator(A, B, kind)`` | Kind = 'normal' or 'anti'.             |
+--------------------------+----------------------------+----------------------------------------+
| Diagonals operator       | ``qdiags(N)``              | Quantum object created from arrays of  |
|                          |                            | diagonals at given offsets.            |
+--------------------------+----------------------------+----------------------------------------+
| Displacement operator    | ``displace(N,alpha)``      | N=number of levels in Hilbert space,   |
| (Single-mode)            |                            | alpha = complex displacement amplitude.|
+--------------------------+----------------------------+----------------------------------------+
| Higher spin operators    | ``jmat(j,#s)``             | j = integer or half-integer            |
|                          |                            | representing spin, s = 'x', 'y', 'z',  |
|                          |                            | '+', or '-'                            |
+--------------------------+----------------------------+----------------------------------------+
| Identity                 | ``qeye(N)``                | N = number of levels in Hilbert space. |
+--------------------------+----------------------------+----------------------------------------+
| Lowering (destruction)   | ``destroy(N)``             | same as above                          |
| operator                 |                            |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Momentum operator        | ``momentum(N)``            | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Number operator          | ``num(N)``                 | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Phase operator           | ``phase(N, phi0)``         | Single-mode Pegg-Barnett phase         |
| (Single-mode)            |                            | operator with ref phase phi0.          |
+--------------------------+----------------------------+----------------------------------------+
| Position operator        | ``position(N)``            | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Raising (creation)       | ``create(N)``              | same as above                          |
| operator                 |                            |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Squeezing operator       | ``squeeze(N, sp)``         | N=number of levels in Hilbert space,   |
| (Single-mode)            |                            | sp = squeezing parameter.              |
+--------------------------+----------------------------+----------------------------------------+
| Squeezing operator       | ``squeezing(q1, q2, sp)``  | q1,q2 = Quantum operators (Qobj)       |
| (Generalized)            |                            | sp = squeezing parameter.              |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-X                  | ``sigmax()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-Y                  | ``sigmay()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-Z                  | ``sigmaz()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma plus               | ``sigmap()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma minus              | ``sigmam()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Tunneling operator       | ``tunneling(N,m)``         | Tunneling operator with elements of the|
|                          |                            | form :math:`|N><N+m| + |N+m><N|`.      |
+--------------------------+----------------------------+----------------------------------------+

As an example, we give the output for a few of these functions:

.. doctest:: [basics]
  :options: +NORMALIZE_WHITESPACE

   >>> basis(5,3)
   Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
   Qobj data =
   [[0.]
    [0.]
    [0.]
    [1.]
    [0.]]

   >>> coherent(5,0.5-0.5j)
   Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
   Qobj data =
   [[ 0.7788017 +0.j        ]
    [ 0.38939142-0.38939142j]
    [ 0.        -0.27545895j]
    [-0.07898617-0.07898617j]
    [-0.04314271+0.j        ]]

   >>> destroy(4)
   Quantum object: dims = [[4], [4]], shape = (4, 4), type = oper, isherm = False
   Qobj data =
   [[0.         1.         0.         0.        ]
    [0.         0.         1.41421356 0.        ]
    [0.         0.         0.         1.73205081]
    [0.         0.         0.         0.        ]]

   >>> sigmaz()
   Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
   Qobj data =
   [[ 1.  0.]
    [ 0. -1.]]

   >>> jmat(5/2.0,'+')
   Quantum object: dims = [[6], [6]], shape = (6, 6), type = oper, isherm = False
   Qobj data =
   [[0.         2.23606798 0.         0.         0.         0.        ]
    [0.         0.         2.82842712 0.         0.         0.        ]
    [0.         0.         0.         3.         0.         0.        ]
    [0.         0.         0.         0.         2.82842712 0.        ]
    [0.         0.         0.         0.         0.         2.23606798]
    [0.         0.         0.         0.         0.         0.        ]]

.. _basics-qobj-props:

Qobj attributes
---------------

We have seen that a quantum object has several internal attributes, such as data, dims, and shape.  These can be accessed in the following way:

.. doctest:: [basics]
  :options: +NORMALIZE_WHITESPACE

   >>> q = destroy(4)

   >>> q.dims
   [[4], [4]]

   >>> q.shape
   (4, 4)

In general, the attributes (properties) of a ``Qobj`` object (or any Python class) can be retrieved using the `Q.attribute` notation.  In addition to the attributes shown with the ``print`` function, the ``Qobj`` class also has the following:

.. cssclass:: table-striped

+---------------+---------------+----------------------------------------+
| Property      | Attribute     | Description                            |
+===============+===============+========================================+
| Data          | ``Q.data``    | Matrix representing state or operator  |
+---------------+---------------+----------------------------------------+
| Dimensions    | ``Q.dims``    | List keeping track of shapes for       |
|               |               | individual components of a             |
|               |               | multipartite system (for tensor        |
|               |               | products and partial traces).          |
+---------------+---------------+----------------------------------------+
| Shape         | ``Q.shape``   | Dimensions of underlying data matrix.  |
+---------------+---------------+----------------------------------------+
| is Hermitian? | ``Q.isherm``  | Is the operator Hermitian or not?      |
+---------------+---------------+----------------------------------------+
| Type          | ``Q.type``    | Is object of type 'ket, 'bra',         |
|               |               | 'oper', or 'super'?                    |
+---------------+---------------+----------------------------------------+

.. figure:: quide-basics-qobj-box.png
   :align: center
   :width: 3.5in

   The ``Qobj`` Class viewed as a container for the properties need to characterize a quantum operator or state vector.


For the destruction operator above:

.. doctest:: [basics]
  :options: +NORMALIZE_WHITESPACE

    >>> q.type
    'oper'

    >>> q.isherm
    False

    >>> q.data
    <4x4 sparse matrix of type '<class 'numpy.complex128'>'
	   with 3 stored elements in Compressed Sparse Row format>



The data attribute returns a message stating that the data is a sparse matrix. All ``Qobj`` instances store their data as a sparse matrix to save memory. To access the underlying dense matrix one needs to use the :func:`qutip.Qobj.full` function as described below.

.. _basics-qobj-math:

Qobj Math
----------

The rules for mathematical operations on ``Qobj`` instances are similar to standard matrix arithmetic:

.. doctest:: [basics]
  :options: +NORMALIZE_WHITESPACE

    >>> q = destroy(4)

    >>> x = sigmax()

    >>> q + 5
    Quantum object: dims = [[4], [4]], shape = (4, 4), type = oper, isherm = False
    Qobj data =
    [[5.         1.         0.         0.        ]
     [0.         5.         1.41421356 0.        ]
     [0.         0.         5.         1.73205081]
     [0.         0.         0.         5.        ]]

    >>> x * x
    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[1. 0.]
     [0. 1.]]

    >>> q ** 3
    Quantum object: dims = [[4], [4]], shape = (4, 4), type = oper, isherm = False
    Qobj data =
    [[0.         0.         0.         2.44948974]
     [0.         0.         0.         0.        ]
     [0.         0.         0.         0.        ]
     [0.         0.         0.         0.        ]]

    >>> x / np.sqrt(2)
    Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
    Qobj data =
    [[0.         0.70710678]
     [0.70710678 0.        ]]


Of course, like matrices, multiplying two objects of incompatible shape throws an error:

.. doctest:: [basics]
  :options: +SKIP

    >>> print(q * x)
    ------------------------------------------------------------------
    TypeError                        Traceback (most recent call last)
    <ipython-input-33-0b599f41213e> in <module>
    ----> 1 print(q * x)

    ~/Documents/qutip_dev/qutip/qutip/qobj.py in __mul__(self, other)
        553
        554             else:
    --> 555                 raise TypeError("Incompatible Qobj shapes")
        556
        557         elif isinstance(other, np.ndarray):

    TypeError: Incompatible Qobj shapes


In addition, the logic operators is equal `==` and is not equal `!=` are also supported.

.. _basics-functions:

Functions operating on Qobj class
==================================

Like attributes, the quantum object class has defined functions (methods) that operate on ``Qobj`` class instances. For a general quantum object ``Q``:

.. cssclass:: table-striped

+-----------------+-------------------------------+----------------------------------------+
| Function        | Command                       | Description                            |
+=================+===============================+========================================+
| Check Hermicity | ``Q.check_herm()``            | Check if quantum object is Hermitian   |
+-----------------+-------------------------------+----------------------------------------+
| Conjugate       | ``Q.conj()``                  | Conjugate of quantum object.           |
+-----------------+-------------------------------+----------------------------------------+
| Cosine          | ``Q.cosm()``                  | Cosine of quantum object.              |
+-----------------+-------------------------------+----------------------------------------+
| Dagger (adjoint)| ``Q.dag()``                   | Returns adjoint (dagger) of object.    |
+-----------------+-------------------------------+----------------------------------------+
| Diagonal        | ``Q.diag()``                  | Returns the diagonal elements.         |
+-----------------+-------------------------------+----------------------------------------+
| Diamond Norm    | ``Q.dnorm()``                 | Returns the diamond norm.              |
+-----------------+-------------------------------+----------------------------------------+
| Eigenenergies   | ``Q.eigenenergies()``         | Eigenenergies (values) of operator.    |
+-----------------+-------------------------------+----------------------------------------+
| Eigenstates     | ``Q.eigenstates()``           | Returns eigenvalues and eigenvectors.  |
+-----------------+-------------------------------+----------------------------------------+
| Eliminate States| ``Q.eliminate_states(inds)``  | Returns quantum object with states in  |
|                 |                               | list inds removed.                     |
+-----------------+-------------------------------+----------------------------------------+
| Exponential     | ``Q.expm()``                  | Matrix exponential of operator.        |
+-----------------+-------------------------------+----------------------------------------+
| Extract States  | ``Q.extract_states(inds)``    | Qobj with states listed in inds only.  |
+-----------------+-------------------------------+----------------------------------------+
| Full            | ``Q.full()``                  | Returns full (not sparse) array of     |
|                 |                               | Q's data.                              |
+-----------------+-------------------------------+----------------------------------------+
| Groundstate     | ``Q.groundstate()``           | Eigenval & eigket of Qobj groundstate. |
+-----------------+-------------------------------+----------------------------------------+
| Matrix Element  | ``Q.matrix_element(bra,ket)`` | Matrix element <bra|Q|ket>             |
+-----------------+-------------------------------+----------------------------------------+
| Norm            | ``Q.norm()``                  | Returns L2 norm for states,            |
|                 |                               | trace norm for operators.              |
+-----------------+-------------------------------+----------------------------------------+
| Overlap         | ``Q.overlap(state)``          | Overlap between current Qobj and a     |
|                 |                               | given state.                           |
+-----------------+-------------------------------+----------------------------------------+
| Partial Trace   | ``Q.ptrace(sel)``             | Partial trace returning components     |
|                 |                               | selected using 'sel' parameter.        |
+-----------------+-------------------------------+----------------------------------------+
| Permute         | ``Q.permute(order)``          | Permutes the tensor structure of a     |
|                 |                               | composite object in the given order.   |
+-----------------+-------------------------------+----------------------------------------+
| Projector       | ``Q.proj()``                  | Form projector operator from given     |
|                 |                               | ket or bra vector.                     |
+-----------------+-------------------------------+----------------------------------------+
| Sine            | ``Q.sinm()``                  | Sine of quantum operator.              |
+-----------------+-------------------------------+----------------------------------------+
| Sqrt            | ``Q.sqrtm()``                 | Matrix sqrt of operator.               |
+-----------------+-------------------------------+----------------------------------------+
| Tidyup          | ``Q.tidyup()``                | Removes small elements from Qobj.      |
+-----------------+-------------------------------+----------------------------------------+
| Trace           | ``Q.tr()``                    | Returns trace of quantum object.       |
+-----------------+-------------------------------+----------------------------------------+
| Transform       | ``Q.transform(inpt)``         | A basis transformation defined by      |
|                 |                               | matrix or list of kets 'inpt' .        |
+-----------------+-------------------------------+----------------------------------------+
| Transpose       | ``Q.trans()``                 | Transpose of quantum object.           |
+-----------------+-------------------------------+----------------------------------------+
| Truncate Neg    | ``Q.trunc_neg()``             | Truncates negative eigenvalues         |
+-----------------+-------------------------------+----------------------------------------+
| Unit            | ``Q.unit()``                  | Returns normalized (unit)              |
|                 |                               | vector Q/Q.norm().                     |
+-----------------+-------------------------------+----------------------------------------+

.. doctest:: [basics]
  :options: +NORMALIZE_WHITESPACE

    >>> basis(5, 3)
    Quantum object: dims = [[5], [1]], shape = (5, 1), type = ket
    Qobj data =
    [[0.]
     [0.]
     [0.]
     [1.]
     [0.]]

    >>> basis(5, 3).dag()
    Quantum object: dims = [[1], [5]], shape = (1, 5), type = bra
    Qobj data =
    [[0. 0. 0. 1. 0.]]

    >>> coherent_dm(5, 1)
    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = True
    Qobj data =
    [[0.36791117 0.36774407 0.26105441 0.14620658 0.08826704]
     [0.36774407 0.36757705 0.26093584 0.14614018 0.08822695]
     [0.26105441 0.26093584 0.18523331 0.10374209 0.06263061]
     [0.14620658 0.14614018 0.10374209 0.05810197 0.035077  ]
     [0.08826704 0.08822695 0.06263061 0.035077   0.0211765 ]]

    >>> coherent_dm(5, 1).diag()
    array([0.36791117, 0.36757705, 0.18523331, 0.05810197, 0.0211765 ])

    >>> coherent_dm(5, 1).full()
    array([[0.36791117+0.j, 0.36774407+0.j, 0.26105441+0.j, 0.14620658+0.j,
            0.08826704+0.j],
           [0.36774407+0.j, 0.36757705+0.j, 0.26093584+0.j, 0.14614018+0.j,
            0.08822695+0.j],
           [0.26105441+0.j, 0.26093584+0.j, 0.18523331+0.j, 0.10374209+0.j,
            0.06263061+0.j],
           [0.14620658+0.j, 0.14614018+0.j, 0.10374209+0.j, 0.05810197+0.j,
            0.035077  +0.j],
           [0.08826704+0.j, 0.08822695+0.j, 0.06263061+0.j, 0.035077  +0.j,
            0.0211765 +0.j]])

    >>> coherent_dm(5, 1).norm()
    1.0000000175063126

    >>> coherent_dm(5, 1).sqrtm()
    Quantum object: dims = [[5], [5]], shape = (5, 5), type = oper, isherm = False
    Qobj data =
    [[0.36791117+3.66778589e-09j 0.36774407-2.13388761e-09j
      0.26105441-1.51480558e-09j 0.14620658-8.48384618e-10j
      0.08826704-5.12182118e-10j]
     [0.36774407-2.13388761e-09j 0.36757705+2.41479965e-09j
      0.26093584-1.11446422e-09j 0.14614018+8.98971115e-10j
      0.08822695+6.40705133e-10j]
     [0.26105441-1.51480558e-09j 0.26093584-1.11446422e-09j
      0.18523331+4.02032413e-09j 0.10374209-3.39161017e-10j
      0.06263061-3.71421368e-10j]
     [0.14620658-8.48384618e-10j 0.14614018+8.98971115e-10j
      0.10374209-3.39161017e-10j 0.05810197+3.36300708e-10j
      0.035077  +2.36883273e-10j]
     [0.08826704-5.12182118e-10j 0.08822695+6.40705133e-10j
      0.06263061-3.71421368e-10j 0.035077  +2.36883273e-10j
      0.0211765 +1.71630348e-10j]]

    >>> coherent_dm(5, 1).tr()
    1.0

    >>> (basis(4, 2) + basis(4, 1)).unit()
    Quantum object: dims = [[4], [1]], shape = (4, 1), type = ket
    Qobj data =
    [[0.        ]
     [0.70710678]
     [0.70710678]
     [0.        ]]
.. _guide:

*******************
Users Guide
*******************

.. toctree::
   :maxdepth: 2

   guide-overview.rst
   guide-basics.rst
   guide-states.rst
   guide-tensor.rst
   guide-dynamics.rst
   guide-heom.rst
   guide-steady.rst
   guide-correlation.rst
   guide-control.rst
   guide-bloch.rst
   guide-visualization.rst
   guide-parfor.rst
   guide-saving.rst
   guide-random.rst
   guide-settings.rst
   guide-qip.rst
   guide-measurement.rst
.. _overview:

******************
Guide Overview
******************

The goal of this guide is to introduce you to the basic structures and functions that make up QuTiP. This guide is divided up into several sections, each highlighting a specific set of functionalities. In combination with the examples that can be found on the project web page `https://qutip.org/tutorials.html <https://qutip.org/tutorials.html>`_, this guide should provide a more or less complete overview. In addition, the :ref:`apidoc` for each function is located at the end of this guide.


.. _overview-org:

Organization
=============

QuTiP is designed to be a general framework for solving quantum mechanics problems such as systems composed of few-level quantum systems and harmonic oscillators. To this end, QuTiP is built from a large (and ever growing) library of functions and classes; from :func:`qutip.states.basis` to :func:`qutip.wigner`.  The general organization of QuTiP, highlighting the important API available to the user, is shown in the figure below.


.. _figure-qutip-org:

.. figure:: figures/qutip_tree.png
   :align: center
   :figwidth: 100%

   Tree-diagram of the 468 user accessible functions and classes in QuTiP 4.6. A vector image of the code tree is in :download:`qutip_tree.pdf <doc/qutip_tree.pdf>`.

.. _parfor:

******************************************
Parallel computation
******************************************

Parallel map and parallel for-loop
----------------------------------

Often one is interested in the output of a given function as a single-parameter is varied. For instance, we can calculate the steady-state response of our system as the driving frequency is varied.  In cases such as this, where each iteration is independent of the others, we can speedup the calculation by performing the iterations in parallel. In QuTiP, parallel computations may be performed using the :func:`qutip.parallel.parallel_map` function or the :func:`qutip.parallel.parfor` (parallel-for-loop) function.

To use the these functions we need to define a function of one or more variables, and the range over which one of these variables are to be evaluated. For example:


.. doctest::
  :skipif: not os_nt
  :options: +NORMALIZE_WHITESPACE

   >>> def func1(x): return x, x**2, x**3

   >>> a, b, c = parfor(func1, range(10))

   >>> print(a)
   [0 1 2 3 4 5 6 7 8 9]

   >>> print(b)
   [ 0  1  4  9 16 25 36 49 64 81]

   >>> print(c)
   [  0   1   8  27  64 125 216 343 512 729]

or

.. doctest::
  :skipif: not os_nt
  :options: +NORMALIZE_WHITESPACE

   >>> result = parallel_map(func1, range(10))

   >>> result_array = np.array(result)

   >>> print(result_array[:, 0])  # == a
   [0 1 2 3 4 5 6 7 8 9]

   >>> print(result_array[:, 1])  # == b
   [ 0  1  4  9 16 25 36 49 64 81]

   >>> print(result_array[:, 2])  # == c
   [  0   1   8  27  64 125 216 343 512 729]


Note that the return values are arranged differently for the :func:`qutip.parallel.parallel_map` and the :func:`qutip.parallel.parfor` functions, as illustrated below. In particular, the return value of :func:`qutip.parallel.parallel_map` is not enforced to be NumPy arrays, which can avoid unnecessary copying if all that is needed is to iterate over the resulting list:


.. doctest::
  :skipif: not os_nt
  :options: +NORMALIZE_WHITESPACE

   >>> result = parfor(func1, range(5))

   >>> print(result)
   [array([0, 1, 2, 3, 4]), array([ 0,  1,  4,  9, 16]), array([ 0,  1,  8, 27, 64])]

   >>> result = parallel_map(func1, range(5))

   >>> print(result)
   [(0, 0, 0), (1, 1, 1), (2, 4, 8), (3, 9, 27), (4, 16, 64)]

The :func:`qutip.parallel.parallel_map` and :func:`qutip.parallel.parfor` functions are not limited to just numbers, but also works for a variety of outputs:

.. doctest::
  :skipif: not os_nt
  :options: +NORMALIZE_WHITESPACE

   >>> def func2(x): return x, Qobj(x), 'a' * x

   >>> a, b, c = parfor(func2, range(5))

   >>> print(a)
   [0 1 2 3 4]

   >>> print(b)
   [Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
   Qobj data =
   [[0.]]
    Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
   Qobj data =
   [[1.]]
    Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
   Qobj data =
   [[2.]]
    Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
   Qobj data =
   [[3.]]
    Quantum object: dims = [[1], [1]], shape = (1, 1), type = bra
   Qobj data =
   [[4.]]]

   >>>print(c)
   ['' 'a' 'aa' 'aaa' 'aaaa']


One can also define functions with **multiple** input arguments and even keyword arguments. Here the :func:`qutip.parallel.parallel_map` and :func:`qutip.parallel.parfor` functions behaves differently:
While :func:`qutip.parallel.parallel_map` only iterate over the values `arguments`, the :func:`qutip.parallel.parfor` function simultaneously iterates over all arguments:

.. doctest::
  :skipif: not os_nt
  :options: +NORMALIZE_WHITESPACE

    >>> def sum_diff(x, y, z=0): return x + y, x - y, z

    >>> parfor(sum_diff, [1, 2, 3], [4, 5, 6], z=5.0)
    [array([5, 7, 9]), array([-3, -3, -3]), array([5., 5., 5.])]

    >>> parallel_map(sum_diff, [1, 2, 3], task_args=(np.array([4, 5, 6]),), task_kwargs=dict(z=5.0))
    [(array([5, 6, 7]), array([-3, -4, -5]), 5.0),
     (array([6, 7, 8]), array([-2, -3, -4]), 5.0),
     (array([7, 8, 9]), array([-1, -2, -3]), 5.0)]

Note that the keyword arguments can be anything you like, but the keyword values are **not** iterated over. The keyword argument *num_cpus* is reserved as it sets the number of CPU's used by parfor. By default, this value is set to the total number of physical processors on your system. You can change this number to a lower value, however setting it higher than the number of CPU's will cause a drop in performance. In :func:`qutip.parallel.parallel_map`, keyword arguments to the task function are specified using `task_kwargs` argument, so there is no special reserved keyword arguments.

The :func:`qutip.parallel.parallel_map` function also supports progressbar, using the keyword argument `progress_bar` which can be set to `True` or to an instance of :class:`qutip.ui.progressbar.BaseProgressBar`. There is a function called :func:`qutip.parallel.serial_map` that works as a non-parallel drop-in replacement for :func:`qutip.parallel.parallel_map`, which allows easy switching between serial and parallel computation.

.. doctest::
  :options: +SKIP

   >>> import time

   >>> def func(x): time.sleep(1)

   >>> result = parallel_map(func, range(50), progress_bar=True)

   10.0%. Run time:   3.10s. Est. time left: 00:00:00:27
   20.0%. Run time:   5.11s. Est. time left: 00:00:00:20
   30.0%. Run time:   8.11s. Est. time left: 00:00:00:18
   40.0%. Run time:  10.15s. Est. time left: 00:00:00:15
   50.0%. Run time:  13.15s. Est. time left: 00:00:00:13
   60.0%. Run time:  15.15s. Est. time left: 00:00:00:10
   70.0%. Run time:  18.15s. Est. time left: 00:00:00:07
   80.0%. Run time:  20.15s. Est. time left: 00:00:00:05
   90.0%. Run time:  23.15s. Est. time left: 00:00:00:02
   100.0%. Run time:  25.15s. Est. time left: 00:00:00:00
   Total run time:  28.91s

Parallel processing is useful for repeated tasks such as generating plots corresponding to the dynamical evolution of your system, or simultaneously simulating different parameter configurations.


IPython-based parallel_map
--------------------------

When QuTiP is used with IPython interpreter, there is an alternative parallel for-loop implementation in the QuTiP  module :func:`qutip.ipynbtools`, see :func:`qutip.ipynbtools.parallel_map`. The advantage of this parallel_map implementation is based on IPython's powerful framework for parallelization, so the compute processes are not confined to run on the same host as the main process.
.. _bloch:

******************************
Plotting on the Bloch Sphere
******************************

.. _bloch-intro:

Introduction
============

When studying the dynamics of a two-level system, it is often convenient to visualize the state of the system by plotting the state-vector or density matrix on the Bloch sphere.  In QuTiP, we have created two different classes to allow for easy creation and manipulation of data sets, both vectors and data points, on the Bloch sphere.  The :class:`qutip.Bloch` class, uses Matplotlib to render the Bloch sphere, where as :class:`qutip.Bloch3d` uses the Mayavi rendering engine to generate a more faithful 3D reconstruction of the Bloch sphere.

.. _bloch-class:

The Bloch and Bloch3d Classes
=============================

In QuTiP, creating a Bloch sphere is accomplished by calling either:

.. plot::
    :context:

    b = qutip.Bloch()

which will load an instance of the :class:`qutip.Bloch` class, or using ::

   >>> b3d = qutip.Bloch3d()

that loads the :class:`qutip.Bloch3d` version.  Before getting into the details of these objects, we can simply plot the blank Bloch sphere associated with these instances via:

.. plot::
    :context:

    b.make_sphere()

or

.. _image-blank3d:

.. figure:: figures/bloch3d-blank.png
    :width: 3.5in
    :figclass: align-center

In addition to the ``show`` command, see the API documentation for :class:`~Bloch` for a full list of other available functions.
As an example, we can add a single data point:

.. plot::
    :context: close-figs

    pnt = [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]
    b.add_points(pnt)
    b.render()

and then a single vector:

.. plot::
    :context: close-figs

    b.fig.clf()
    vec = [0, 1, 0]
    b.add_vectors(vec)
    b.render()

and then add another vector corresponding to the :math:`\left|\rm up \right>` state:

.. plot::
    :context: close-figs

    up = qutip.basis(2, 0)
    b.add_states(up)
    b.render()

Notice that when we add more than a single vector (or data point), a different color will automatically be applied to the later data set (mod 4).
In total, the code for constructing our Bloch sphere with one vector, one state, and a single data point is:

.. plot::
    :context: close-figs

    b = qutip.Bloch()

    pnt = [1./np.sqrt(3), 1./np.sqrt(3), 1./np.sqrt(3)]
    b.add_points(pnt)
    vec = [0, 1, 0]
    b.add_vectors(vec)
    up = qutip.basis(2, 0)
    b.add_states(up)
    b.render()

where we have removed the extra ``show()`` commands.  Replacing ``b=Bloch()`` with ``b=Bloch3d()`` in the above code generates the following 3D Bloch sphere.

.. _image-bloch3ddata:

.. figure:: figures/bloch3d+data.png
    :width: 3.5in
    :figclass: align-center


We can also plot multiple points, vectors, and states at the same time by passing list or arrays instead of individual elements.  Before giving an example, we can use the `clear()` command to remove the current data from our Bloch sphere instead of creating a new instance:

.. plot::
    :context: close-figs

    b.clear()
    b.render()


Now on the same Bloch sphere, we can plot the three states associated with the x, y, and z directions:

.. plot::
    :context: close-figs

    x = (qutip.basis(2, 0) + (1+0j)*qutip.basis(2, 1)).unit()
    y = (qutip.basis(2, 0) + (0+1j)*qutip.basis(2, 1)).unit()
    z = (qutip.basis(2, 0) + (0+0j)*qutip.basis(2, 1)).unit()

    b.add_states([x, y, z])
    b.render()

a similar method works for adding vectors:

.. plot::
    :context: close-figs

    b.clear()
    vec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    b.add_vectors(vec)
    b.render()

Adding multiple points to the Bloch sphere works slightly differently than adding multiple states or vectors.  For example, lets add a set of 20 points around the equator (after calling `clear()`):

.. plot::
    :context: close-figs

    b.clear()

    th = np.linspace(0, 2*np.pi, 20)
    xp = np.cos(th)
    yp = np.sin(th)
    zp = np.zeros(20)

    pnts = [xp, yp, zp]
    b.add_points(pnts)
    b.render()

Notice that, in contrast to states or vectors, each point remains the same color as the initial point.  This is because adding multiple data points using the ``add_points`` function is interpreted, by default, to correspond to a single data point (single qubit state) plotted at different times.  This is very useful when visualizing the dynamics of a qubit.  An example of this is given in the example .  If we want to plot additional qubit states we can call additional ``add_points`` functions:

.. plot::
    :context: close-figs

    xz = np.zeros(20)
    yz = np.sin(th)
    zz = np.cos(th)
    b.add_points([xz, yz, zz])
    b.render()

The color and shape of the data points is varied automatically by the Bloch class.  Notice how the color and point markers change for each set of data.  Again, we have had to call ``add_points`` twice because adding more than one set of multiple data points is *not* supported by the ``add_points`` function.

What if we want to vary the color of our points.  We can tell the :class:`qutip.Bloch` class to vary the color of each point according to the colors listed in the ``b.point_color`` list (see :ref:`bloch-config` below).  Again after ``clear()``:

.. plot::
    :context: close-figs

    b.clear()

    xp = np.cos(th)
    yp = np.sin(th)
    zp = np.zeros(20)
    pnts = [xp, yp, zp]
    b.add_points(pnts, 'm')  # <-- add a 'm' string to signify 'multi' colored points
    b.render()


Now, the data points cycle through a variety of predefined colors.  Now lets add another set of points, but this time we want the set to be a single color, representing say a qubit going from the :math:`\left|\rm up\right>` state to the :math:`\left|\rm down\right>` state in the y-z plane:

.. plot::
    :context: close-figs

    xz = np.zeros(20)
    yz = np.sin(th)
    zz = np.cos(th)

    b.add_points([xz, yz, zz])  # no 'm'
    b.render()

Again, the same plot can be generated using the :class:`qutip.Bloch3d` class by replacing ``Bloch`` with ``Bloch3d``:

.. figure:: figures/bloch3d+points.png
    :width: 3.5in
    :figclass: align-center

A more slick way of using this 'multi' color feature is also given in the example, where we set the color of the markers as a function of time.

Differences Between Bloch and Bloch3d
-------------------------------------
While in general the ``Bloch`` and ``Bloch3d`` classes are interchangeable, there are some important differences to consider when choosing between them.

- The ``Bloch`` class uses Matplotlib to generate figures.  As such, the data plotted on the sphere is in reality just a 2D object.  In contrast the ``Bloch3d`` class uses the 3D rendering engine from VTK via mayavi to generate the sphere and the included data.  In this sense the ``Bloch3d`` class is much more advanced, as objects are rendered in 3D leading to a higher quality figure.

- Only the ``Bloch`` class can be embedded in a Matplotlib figure window.  Thus if you want to combine a Bloch sphere with another figure generated in QuTiP, you can not use ``Bloch3d``.  Of course you can always post-process your figures using other software to get the desired result.

- Due to limitations in the rendering engine, the ``Bloch3d`` class does not support LaTeX for text.  Again, you can get around this by post-processing.

- The user customizable attributes for the ``Bloch`` and ``Bloch3d`` classes are not identical.  Therefore, if you change the properties of one of the classes, these changes will cause an exception if the class is switched.


.. _bloch-config:

Configuring the Bloch sphere
============================

Bloch Class Options
--------------------

At the end of the last section we saw that the colors and marker shapes of the data plotted on the Bloch sphere are automatically varied according to the number of points and vectors added.  But what if you want a different choice of color, or you want your sphere to be purple with different axes labels? Well then you are in luck as the Bloch class has 22 attributes which one can control.  Assuming ``b=Bloch()``:

.. tabularcolumns:: | p{3cm} | p{7cm} |  p{7cm} |

.. cssclass:: table-striped

+---------------+---------------------------------------------------------+-------------------------------------------------+
| Attribute     | Function                                                | Default Setting                                 |
+===============+=========================================================+=================================================+
| b.axes        | Matplotlib axes instance for animations. Set by ``axes``| ``None``                                        |
|               | keyword arg.                                            |                                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.fig         | User supplied Matplotlib Figure instance. Set by ``fig``| ``None``                                        |
|               | keyword arg.                                            |                                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.font_color  | Color of fonts                                          | 'black'                                         |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.font_size   | Size of fonts                                           | 20                                              |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.frame_alpha | Transparency of wireframe                               | 0.1                                             |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.frame_color | Color of wireframe                                      | 'gray'                                          |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.frame_width | Width of wireframe                                      | 1                                               |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.point_color | List of colors for Bloch point markers to cycle through | ``['b', 'r', 'g', '#CC6600']``                  |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.point_marker| List of point marker shapes to cycle through            | ``['o', 's', 'd', '^']``                        |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.point_size  | List of point marker sizes (not all markers look the    | ``[55, 62, 65, 75]``                            |
|               | same size when plotted)                                 |                                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.sphere_alpha| Transparency of Bloch sphere                            | 0.2                                             |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.sphere_color| Color of Bloch sphere                                   | ``'#FFDDDD'``                                   |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.size        | Sets size of figure window                              | ``[7, 7]`` (700x700 pixels)                     |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.vector_color| List of colors for Bloch vectors to cycle through       | ``['g', '#CC6600', 'b', 'r']``                  |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.vector_width| Width of Bloch vectors                                  | 4                                               |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.view        | Azimuthal and Elevation viewing angles                  | ``[-60,30]``                                    |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.xlabel      | Labels for x-axis                                       | ``['$x$', '']`` +x and -x (labels use LaTeX)    |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.xlpos       | Position of x-axis labels                               | ``[1.1, -1.1]``                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.ylabel      | Labels for y-axis                                       | ``['$y$', '']`` +y and -y (labels use LaTeX)    |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.ylpos       | Position of y-axis labels                               | ``[1.2, -1.2]``                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.zlabel      | Labels for z-axis                                       | ``['$\left|0\right>$', '$\left|1\right>$']``    |
|               |                                                         | +z and -z (labels use LaTeX)                    |
+---------------+---------------------------------------------------------+-------------------------------------------------+
| b.zlpos       | Position of z-axis labels                               | ``[1.2, -1.2]``                                 |
+---------------+---------------------------------------------------------+-------------------------------------------------+

Bloch3d Class Options
---------------------

The Bloch3d sphere is also customizable.  Note however that the attributes for the ``Bloch3d`` class are not in one-to-one
correspondence to those of the ``Bloch`` class due to the different underlying rendering engines. Assuming ``b=Bloch3d()``:

.. tabularcolumns:: | p{3cm} | p{7cm} |  p{7cm} |

.. cssclass:: table-striped

+---------------+---------------------------------------------------------+---------------------------------------------+
| Attribute     | Function                                                | Default Setting                             |
+===============+=========================================================+=============================================+
| b.fig         | User supplied Mayavi Figure instance. Set by ``fig``    | ``None``                                    |
|               | keyword arg.                                            |                                             |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.font_color  | Color of fonts                                          | ``'black'``                                 |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.font_scale  | Scale of fonts                                          | 0.08                                        |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.frame       | Draw wireframe for sphere?                              | ``True``                                    |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.frame_alpha | Transparency of wireframe                               | 0.05                                        |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.frame_color | Color of wireframe                                      | ``'gray'``                                  |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.frame_num   | Number of wireframe elements to draw                    | 8                                           |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.frame_radius| Radius of wireframe lines                               | 0.005                                       |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.point_color | List of colors for Bloch point markers to cycle through | ``['r', 'g', 'b', 'y']``                    |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.point_mode  | Type of point markers to draw                           | ``'sphere'``                                |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.point_size  | Size of points                                          | 0.075                                       |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.sphere_alpha| Transparency of Bloch sphere                            | 0.1                                         |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.sphere_color| Color of Bloch sphere                                   | ``'#808080'``                               |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.size        | Sets size of figure window                              | ``[500, 500]`` (500x500 pixels)             |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.vector_color| List of colors for Bloch vectors to cycle through       | ``['r', 'g', 'b', 'y']``                    |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.vector_width| Width of Bloch vectors                                  | 3                                           |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.view        | Azimuthal and Elevation viewing angles                  | ``[45, 65]``                                |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.xlabel      | Labels for x-axis                                       | ``['|x>', '']`` +x and -x                   |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.xlpos       | Position of x-axis labels                               | ``[1.07, -1.07]``                           |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.ylabel      | Labels for y-axis                                       | ``['$y$', '']`` +y and -y                   |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.ylpos       | Position of y-axis labels                               | ``[1.07, -1.07]``                           |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.zlabel      | Labels for z-axis                                       | ``['|0>', '|1>']`` +z and -z                |
+---------------+---------------------------------------------------------+---------------------------------------------+
| b.zlpos       | Position of z-axis labels                               | ``[1.07, -1.07]``                           |
+---------------+---------------------------------------------------------+---------------------------------------------+

These properties can also be accessed via the print command:

.. doctest::

    >>> b = qutip.Bloch()

    >>> print(b) # doctest: +NORMALIZE_WHITESPACE
    Bloch data:
    -----------
    Number of points:  0
    Number of vectors: 0
    <BLANKLINE>
    Bloch sphere properties:
    ------------------------
    font_color:      black
    font_size:       20
    frame_alpha:     0.2
    frame_color:     gray
    frame_width:     1
    point_color:     ['b', 'r', 'g', '#CC6600']
    point_marker:    ['o', 's', 'd', '^']
    point_size:      [25, 32, 35, 45]
    sphere_alpha:    0.2
    sphere_color:    #FFDDDD
    figsize:         [5, 5]
    vector_color:    ['g', '#CC6600', 'b', 'r']
    vector_width:    3
    vector_style:    -|>
    vector_mutation: 20
    view:            [-60, 30]
    xlabel:          ['$x$', '']
    xlpos:           [1.2, -1.2]
    ylabel:          ['$y$', '']
    ylpos:           [1.2, -1.2]
    zlabel:          ['$\\left|0\\right>$', '$\\left|1\\right>$']
    zlpos:           [1.2, -1.2]
    <BLANKLINE>

.. _bloch-animate:

Animating with the Bloch sphere
===============================

The Bloch class was designed from the outset to generate animations.  To animate a set of vectors or data points the basic idea is: plot the data at time t1, save the sphere, clear the sphere, plot data at t2,... The Bloch sphere will automatically number the output file based on how many times the object has been saved (this is stored in b.savenum).  The easiest way to animate data on the Bloch sphere is to use the ``save()`` method and generate a series of images to convert into an animation.  However, as of Matplotlib version 1.1, creating animations is built-in.  We will demonstrate both methods by looking at the decay of a qubit on the bloch sphere.

.. _bloch-animate-decay:

Example: Qubit Decay
--------------------

The code for calculating the expectation values for the Pauli spin operators of a qubit decay is given below.  This code is common to both animation examples.

.. literalinclude:: scripts/ex_bloch_animation.py

.. _bloch-animate-decay-images:

Generating Images for Animation
++++++++++++++++++++++++++++++++

An example of generating images for generating an animation outside of Python is given below::

     import numpy as np
     b = qutip.Bloch()
     b.vector_color = ['r']
     b.view = [-40, 30]
     for i in range(len(sx)):
         b.clear()
         b.add_vectors([np.sin(theta), 0, np.cos(theta)])
         b.add_points([sx[:i+1], sy[:i+1], sz[:i+1]])
         b.save(dirc='temp')  # saving images to temp directory in current working directory

Generating an animation using FFmpeg (for example) is fairly simple::

   ffmpeg -i temp/bloch_%01d.png bloch.mp4

.. _bloch-animate-decay-direct:

Directly Generating an Animation
++++++++++++++++++++++++++++++++

.. important::
   Generating animations directly from Matplotlib requires installing either MEncoder or FFmpeg.
   While either choice works on linux, it is best to choose FFmpeg when running on the Mac.
   If using macports just do: ``sudo port install ffmpeg``.

The code to directly generate an mp4 movie of the Qubit decay is as follows ::

   from matplotlib import pyplot, animation
   from mpl_toolkits.mplot3d import Axes3D

   fig = pyplot.figure()
   ax = Axes3D(fig, azim=-40, elev=30)
   sphere = qutip.Bloch(axes=ax)

   def animate(i):
      sphere.clear()
      sphere.add_vectors([np.sin(theta), 0, np.cos(theta)])
      sphere.add_points([sx[:i+1], sy[:i+1], sz[:i+1]])
      sphere.make_sphere()
      return ax

   def init():
      sphere.vector_color = ['r']
      return ax

   ani = animation.FuncAnimation(fig, animate, np.arange(len(sx)),
                                 init_func=init, blit=False, repeat=False)
   ani.save('bloch_sphere.mp4', fps=20)

The resulting movie may be viewed here: `bloch_decay.mp4 <https://raw.githubusercontent.com/qutip/qutip/master/doc/figures/bloch_decay.mp4>`_
.. _steady:

*************************************
Solving for Steady-State Solutions
*************************************

.. _steady-intro:

Introduction
============

For time-independent open quantum systems with decay rates larger than the corresponding excitation rates, the system will tend toward a steady state as :math:`t\rightarrow\infty` that satisfies the equation

.. math::
    \frac{d\hat{\rho}_{ss}}{dt}=\mathcal{L}\hat{\rho}_{ss}=0.

Although the requirement for time-independence seems quite resitrictive, one can often employ a transformation to the interaction picture that yields a time-independent Hamiltonian.  For many these systems, solving for the asymptotic density matrix :math:`\hat{\rho}_{ss}` can be achieved using direct or iterative solution methods faster than using master equation or Monte Carlo simulations.  Although the steady state equation has a simple mathematical form, the properties of the Liouvillian operator are such that the solutions to this equation are anything but straightforward to find.

Steady State solvers in QuTiP
=============================

In QuTiP, the steady-state solution for a system Hamiltonian or Liouvillian is given by :func:`qutip.steadystate.steadystate`.  This function implements a number of different methods for finding the steady state, each with their own pros and cons, where the method used can be chosen using the ``method`` keyword argument.

.. cssclass:: table-striped

.. list-table::
   :widths: 10 15 30
   :header-rows: 1

   * - Method
     - Keyword
     - Description
   * - Direct (default)
     - 'direct'
     - Direct solution solving :math:`Ax=b` via sparse LU decomposition.
   * - Eigenvalue
     - 'eigen'
     - Iteratively find the zero eigenvalue of :math:`\mathcal{L}`.
   * - Inverse-Power
     - 'power'
     - Solve using the inverse-power method.
   * - GMRES
     - 'iterative-gmres'
     - Solve using the GMRES method and optional preconditioner.
   * - LGMRES
     - 'iterative-lgmres'
     - Solve using the LGMRES method and optional preconditioner.
   * - BICGSTAB
     - 'iterative-bicgstab'
     - Solve using the BICGSTAB method and optional preconditioner.
   * - SVD
     - 'svd'
     - Steady-state solution via the **dense** SVD of the Liouvillian.


The function :func:`qutip.steadystate.steadystate` can take either a Hamiltonian and a list of collapse operators as input, generating internally the corresponding Liouvillian super operator in Lindblad form, or alternatively, a Liouvillian passed by the user. When possible, we recommend passing the Hamiltonian and collapse operators to :func:`qutip.steadystate.steadystate`, and letting the function automatically build the Liouvillian (in Lindblad form) for the system.

As of QuTiP 3.2, the ``direct`` and ``power`` methods can take advantage of the Intel Pardiso LU solver in the Intel Math Kernel library that comes with the Anacoda (2.5+) and Intel Python distributions.  This gives a substantial increase in performance compared with the standard SuperLU method used by SciPy.  To verify that QuTiP can find the necessary libraries, one can check for ``INTEL MKL Ext: True`` in the QuTiP about box (:func:`qutip.about`).


.. _steady-usage:

Using the Steadystate Solver
=============================

Solving for the steady state solution to the Lindblad master equation for a general system with :func:`qutip.steadystate.steadystate` can be accomplished using::

>>> rho_ss = steadystate(H, c_ops)

where ``H`` is a quantum object representing the system Hamiltonian, and ``c_ops`` is a list of quantum objects for the system collapse operators. The output, labeled as ``rho_ss``, is the steady-state solution for the systems.  If no other keywords are passed to the solver, the default 'direct' method is used, generating a solution that is exact to machine precision at the expense of a large memory requirement.  The large amount of memory need for the direct LU decomposition method stems from the large bandwidth of the system Liouvillian and the correspondingly large fill-in (extra nonzero elements) generated in the LU factors.  This fill-in can be reduced by using bandwidth minimization algorithms such as those discussed in :ref:`steady-args`.  However, in most cases, the default fill-in reducing algorithm is nearly optimal.  Additional parameters may be used by calling the steady-state solver as:

.. code-block:: python

   rho_ss = steadystate(H, c_ops, method='power', use_rcm=True)

where ``method='power'`` indicates that we are using the inverse-power solution method, and ``use_rcm=True`` turns on a bandwidth minimization routine.


Although it is not obvious, the ``'direct'``, ``eigen``, and ``'power'`` methods all use an LU decomposition internally and thus suffer from a large memory overhead.  In contrast, iterative methods such as the ``'iterative-gmres'``, ``'iterative-lgmres'``, and ``'iterative-bicgstab'`` methods do not factor the matrix and thus take less memory than these previous methods and allowing, in principle, for extremely large system sizes. The downside is that these methods can take much longer than the direct method as the condition number of the Liouvillian matrix is large, indicating that these iterative methods require a large number of iterations for convergence.  To overcome this, one can use a preconditioner :math:`M` that solves for an approximate inverse for the (modified) Liouvillian, thus better conditioning the problem, leading to faster convergence.  The use of a preconditioner can actually make these iterative methods faster than the other solution methods.  The problem with precondioning is that it is only well defined for Hermitian matrices.  Since the Liouvillian is non-Hermitian, the ability to find a good preconditioner is not guaranteed.  And moreover, if a preconditioner is found, it is not guaranteed to have a good condition number. QuTiP can make use of an incomplete LU preconditioner when using the iterative ``'gmres'``, ``'lgmres'``, and ``'bicgstab'`` solvers by setting ``use_precond=True``. The preconditioner optionally makes use of a combination of symmetric and anti-symmetric matrix permutations that attempt to improve the preconditioning process.  These features are discussed in the :ref:`steady-args` section.  Even with these state-of-the-art permutations, the generation of a successful preconditoner for non-symmetric matrices is currently a trial-and-error process due to the lack of mathematical work done in this area.  It is always recommended to begin with the direct solver with no additional arguments before selecting a different method.

Finding the steady-state solution is not limited to the Lindblad form of the master equation. Any time-independent Liouvillian constructed from a Hamiltonian and collapse operators can be used as an input::

>>> rho_ss = steadystate(L)

where ``L`` is the Louvillian.  All of the additional arguments can also be used in this case.


.. _steady-args:

Additional Solver Arguments
=============================

The following additional solver arguments are available for the steady-state solver:

.. cssclass:: table-striped

.. list-table::
   :widths: 10 30 60
   :header-rows: 1

   * - Keyword
     - Options (default listed first)
     - Description
   * - method
     - 'direct', 'eigen', 'power', 'iterative-gmres','iterative-lgmres', 'svd'
     - Method used for solving for the steady-state density matrix.
   * - sparse
     - True, False
     - Use sparse version of direct solver.
   * - weight
     - None
     - Allows the user to define the weighting factor used in the ``'direct'``, ``'GMRES'``, and ``'LGMRES'`` solvers.
   * - permc_spec
     - 'COLAMD', 'NATURAL'
     - Column ordering used in the sparse LU decomposition.
   * - use_rcm
     - False, True
     - Use a Reverse Cuthill-Mckee reordering to minimize the bandwidth of the modified Liouvillian used in the LU decomposition.  If ``use_rcm=True`` then the column ordering is set to ``'Natural'`` automatically unless explicitly set.
   * - use_precond
     - False, True
     - Attempt to generate a preconditioner when using the ``'iterative-gmres'`` and ``'iterative-lgmres'`` methods.
   * - M
     - None, sparse_matrix, LinearOperator
     - A user defined preconditioner, if any.
   * - use_wbm
     - False, True
     - Use a Weighted Bipartite Matching algorithm to attempt to make the modified Liouvillian more diagonally dominate, and thus for favorable for preconditioning.  Set to ``True`` automatically when using a iterative method, unless explicitly set.
   * - tol
     - 1e-9
     - Tolerance used in finding the solution for all methods expect ``'direct'`` and ``'svd'``.
   * - maxiter
     - 10000
     - Maximum number of iterations to perform for all methods expect ``'direct'`` and ``'svd'``.
   * - fill_factor
     - 10
     - Upper-bound on the allowed fill-in for the approximate inverse preconditioner.  This value may need to be set much higher than this in some cases.
   * - drop_tol
     - 1e-3
     - Sets the threshold for the relative magnitude of preconditioner elements that should be dropped.  A lower number yields a more accurate approximate inverse at the expense of fill-in and increased runtime.
   * - diag_pivot_thresh
     - None
     - Sets the threshold between :math:`[0,1]` for which diagonal elements are considered acceptable pivot points when using a preconditioner.
   * - ILU_MILU
     - 'smilu_2'
     - Selects the incomplete LU decomposition method algorithm used.

Further information can be found in the :func:`qutip.steadystate.steadystate` docstrings.


.. _steady-example:

Example: Harmonic Oscillator in Thermal Bath
============================================

A simple example of a system that reaches a steady state is a harmonic oscillator coupled to a thermal environment.  Below we consider a harmonic oscillator, initially in the :math:`\left|10\right>` number state, and weakly coupled to a thermal environment characterized by an average particle expectation value of :math:`\left<n\right>=2`.  We calculate the evolution via master equation and Monte Carlo methods, and see that they converge to the steady-state solution.  Here we choose to perform only a few Monte Carlo trajectories so we can distinguish this evolution from the master-equation solution.

.. plot:: guide/scripts/ex_steady.py
   :include-source:
####################
Bosonic Environments
####################

In this section we consider a simple two-level system coupled to a
Drude-Lorentz bosonic bath. The system Hamiltonian, :math:`H_{sys}`, and the bath
spectral density, :math:`J_D`, are

.. math::

    H_{sys} &= \frac{\epsilon \sigma_z}{2} + \frac{\Delta \sigma_x}{2}

    J_D &= \frac{2\lambda \gamma \omega}{(\gamma^2 + \omega^2)},

We will demonstrate how to describe the bath using two different expansions
of the spectral density correlation function (Matsubara's expansion and
a Padé expansion), how to evolve the system in time, and how to calculate
the steady state.

First we will do this in the simplest way, using the built-in implementations of
the two bath expansions, :class:`~qutip.nonmarkov.heom.DrudeLorentzBath` and
:class:`~qutip.nonmarkov.heom.DrudeLorentzPadeBath`. We will do this both with a
truncated expansion and show how to include an approximation to all of the
remaining terms in the bath expansion.

Afterwards, we will show how to calculate the bath expansion coefficients and to
use those coefficients to construct your own bath description so that you can
implement your own bosonic baths.

Finally, we will demonstrate how to simulate a system coupled to multiple
independent baths, as occurs, for example, in certain photosynthesis processes.

A notebook containing a complete example similar to this one implemented in
BoFiN can be found in
`example notebook 1a <https://github.com/tehruhn/bofin/blob/main/examples/example-1a-Spin-bath-model-basic.ipynb>`__.


Describing the system and bath
------------------------------

First, let us construct the system Hamiltonian, :math:`H_{sys}`, and the initial
system state, ``rho0``:

.. plot::
    :context: reset
    :nofigs:

    from qutip import basis, sigmax, sigmaz

    # The system Hamiltonian:
    eps = 0.5  # energy of the 2-level system
    Del = 1.0  # tunnelling term
    H_sys = 0.5 * eps * sigmaz() + 0.5 * Del * sigmax()

    # Initial state of the system:
    rho0 = basis(2,0) * basis(2,0).dag()

Now let us describe the bath properties:

.. plot::
    :context:
    :nofigs:

    # Bath properties:
    gamma = 0.5  # cut off frequency
    lam = 0.1  # coupling strength
    T = 0.5  # temperature

    # System-bath coupling operator:
    Q = sigmaz()

where :math:`\gamma` (``gamma``), :math:`\lambda` (``lam``) and :math:`T` are
the parameters of a Drude-Lorentz bath, and ``Q`` is the coupling operator
between the system and the bath.

We may the pass these parameters to either
:class:`~qutip.nonmarkov.heom.DrudeLorentzBath` or
:class:`~qutip.nonmarkov.heom.DrudeLorentzPadeBath` to construct an expansion of
the bath correlations:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import DrudeLorentzBath
    from qutip.nonmarkov.heom import DrudeLorentzPadeBath

    # Number of expansion terms to retain:
    Nk = 2

    # Matsubara expansion:
    bath = DrudeLorentzBath(Q, lam, gamma, T, Nk)

    # Padé expansion:
    bath = DrudeLorentzPadeBath(Q, lam, gamma, T, Nk)

Where ``Nk`` is the number of terms to retain within the expansion of the
bath.


.. _heom-bosonic-system-and-bath-dynamics:

System and bath dynamics
------------------------

Now we are ready to construct a solver:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import HEOMSolver
    from qutip import Options

    max_depth = 5  # maximum hierarchy depth to retain
    options = Options(nsteps=15_000)

    solver = HEOMSolver(H_sys, bath, max_depth=max_depth, options=options)

and to calculate the system evolution as a function of time:

.. code-block:: python

    tlist = [0, 10, 20]  # times to evaluate the system state at
    result = solver.run(rho0, tlist)

The ``max_depth`` parameter determines how many levels of the hierarchy to
retain. As a first approximation hierarchy depth may be thought of as similar
to the order of Feynman Diagrams (both classify terms by increasing number
of interactions).

The ``result`` is a standard QuTiP results object with the attributes:

- ``times``: the times at which the state was evaluated (i.e. ``tlist``)
- ``states``: the system states at each time
- ``expect``: the values of each ``e_ops`` at each time
- ``ado_states``: see below (an instance of
  :class:`~qutip.nonmarkov.heom.HierarchyADOsState`)

If ``ado_return=True`` is passed to ``.run(...)`` the full set of auxilliary
density operators (ADOs) that make up the hierarchy at each time will be
returned as ``.ado_states``. We will describe how to use these to determine
other properties, such as system-bath currents, later in the fermionic guide
(see :ref:`heom-determining-currents`).

If one has a full set of ADOs from a previous call of ``.run(...)`` you may
supply it as the initial state of the solver by calling
``.run(result.ado_states[-1], tlist, ado_init=True)``.

As with other QuTiP solvers, if expectation operators or functions are supplied
using ``.run(..., e_ops=[...])`` the expectation values are available in
``result.expect``.

Below we run the solver again, but use ``e_ops`` to store the expectation
values of the population of the system states and the coherence:

.. plot::
    :context:

    # Define the operators that measure the populations of the two
    # system states:
    P11p = basis(2,0) * basis(2,0).dag()
    P22p = basis(2,1) * basis(2,1).dag()

    # Define the operator that measures the 0, 1 element of density matrix
    # (corresonding to coherence):
    P12p = basis(2,0) * basis(2,1).dag()

    # Run the solver:
    tlist = np.linspace(0, 20, 101)
    result = solver.run(rho0, tlist, e_ops={"11": P11p, "22": P22p, "12": P12p})

    # Plot the results:
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,8))
    axes.plot(result.times, result.expect["11"], 'b', linewidth=2, label="P11")
    axes.plot(result.times, result.expect["12"], 'r', linewidth=2, label="P12")
    axes.set_xlabel(r't', fontsize=28)
    axes.legend(loc=0, fontsize=12)


Steady-state
------------

Using the same solver, we can also determine the steady state of the
combined system and bath using:

.. plot::
    :context:
    :nofigs:

    steady_state, steady_ados = solver.steady_state()

where ``steady_state`` is the steady state of the system and ``steady_ados``
if the steady state of the full hierarchy. The ADO states are
described more fully in :ref:`heom-determining-currents` and
:class:`~qutip.nonmarkov.heom.HierarchyADOsState`.


Matsubara Terminator
--------------------

When constructing the Drude-Lorentz bath we have truncated the expansion at
``Nk = 2`` terms and ignore the remaining terms.

However, since the coupling to these higher order terms is comparatively weak,
we may consider the interaction with them to be Markovian, and construct an
additional Lindbladian term that captures their interaction with the system and
the lower order terms in the expansion.

This additional term is called the ``terminator`` because it terminates the
expansion.

The :class:`~qutip.nonmarkov.heom.DrudeLorentzBath` and
:class:`~qutip.nonmarkov.heom.DrudeLorentzPadeBath` both provide a means of
calculating the terminator for a given expansion:

.. plot::
    :context:
    :nofigs:

    # Matsubara expansion:
    bath = DrudeLorentzBath(Q, lam, gamma, T, Nk)

    # Padé expansion:
    bath = DrudeLorentzPadeBath(Q, lam, gamma, T, Nk)

    # Add terminator to the system Liouvillian:
    delta, terminator = bath.terminator()
    HL = liouvillian(H_sys) + terminator

    # Construct solver:
    solver = HEOMSolver(HL, bath, max_depth=max_depth, options=options)

This captures the Markovian effect of the remaining terms in the expansion
without having to fully model many more terms.

The value ``delta`` is an approximation to the strength of the effect of
the remaining terms in the expansion (i.e. how strongly the terminator is
coupled to the rest of the system).


Matsubara expansion coefficients
--------------------------------

So far we have relied on the built-in
:class:`~qutip.nonmarkov.heom.DrudeLorentzBath` to construct the Drude-Lorentz
bath expansion for us. Now we will calculate the coefficients ourselves and
construct a :class:`~qutip.nonmarkov.heom.BosonicBath` directly. A similar
procedure can be used to apply :class:`~qutip.nonmarkov.heom.HEOMSolver` to any
bosonic bath for which we can calculate the expansion coefficients.

The real and imaginary parts of the correlation function, :math:`C(t)`, for the
bosonic bath is expanded in an expontential series:

.. math::

      C(t) &= C_{real}(t) + i C_{imag}(t)

      C_{real}(t) &= \sum_{k=0}^{\infty} c_{k,real} e^{- \nu_{k,real} t}

      C_{imag}(t) &= \sum_{k=0}^{\infty} c_{k,imag} e^{- \nu_{k,imag} t}

In the specific case of Matsubara expansion for the Drude-Lorentz bath, the
coefficients of this expansion are, for the real part, :math:`C_{real}(t)`:

.. math::

    \nu_{k,real} &= \begin{cases}
        \gamma                & k = 0\\
        {2 \pi k} / {\beta }  & k \geq 1\\
    \end{cases}

    c_{k,real} &= \begin{cases}
        \lambda \gamma [\cot(\beta \gamma / 2) - i]             & k = 0\\
        \frac{4 \lambda \gamma \nu_k }{ (\nu_k^2 - \gamma^2)\beta}    & k \geq 1\\
    \end{cases}

and the imaginary part, :math:`C_{imag}(t)`:

.. math::

    \nu_{k,imag} &= \begin{cases}
        \gamma                & k = 0\\
        0                     & k \geq 1\\
    \end{cases}

    c_{k,imag} &= \begin{cases}
        - \lambda \gamma      & k = 0\\
        0                     & k \geq 1\\
    \end{cases}

And now the same numbers calculated in Python:

.. plot::
    :context:
    :nofigs:

    # Convenience functions and parameters:

    def cot(x):
        return 1. / np.tan(x)

    beta = 1. / T

    # Number of expansion terms to calculate:
    Nk = 2

    # C_real expansion terms:
    ck_real = [lam * gamma / np.tan(gamma / (2 * T))]
    ck_real.extend([
        (8 * lam * gamma * T * np.pi * k * T /
            ((2 * np.pi * k * T)**2 - gamma**2))
        for k in range(1, Nk + 1)
    ])
    vk_real = [gamma]
    vk_real.extend([2 * np.pi * k * T for k in range(1, Nk + 1)])

    # C_imag expansion terms (this is the full expansion):
    ck_imag = [lam * gamma * (-1.0)]
    vk_imag = [gamma]

After all that, constructing the bath is very straight forward:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import BosonicBath

    bath = BosonicBath(Q, ck_real, vk_real, ck_imag, vk_imag)

And we're done!

The :class:`~qutip.nonmarkov.heom.BosonicBath` can be used with the
:class:`~qutip.nonmarkov.heom.HEOMSolver` in exactly the same way as the baths
we constructed previously using the built-in Drude-Lorentz bath expansions.


Multiple baths
--------------

The :class:`~qutip.nonmarkov.heom.HEOMSolver` supports having a system interact
with multiple environments. All that is needed is to supply a list of baths
instead of a single bath.

In the example below we calculate the evolution of a small system where
each basis state of the system interacts with a separate bath. Such
an arrangement can model, for example, the Fenna–Matthews–Olson (FMO)
pigment-protein complex which plays an important role in photosynthesis (
for a full FMO example see the notebook
https://github.com/tehruhn/bofin/blob/main/examples/example-2-FMO-example.ipynb
).

For each bath expansion, we also include the terminator in the system
Liouvillian.

At the end, we plot the populations of the system states as a function of
time, and show the long-time beating of quantum state coherence that
occurs:

.. plot::
    :context: close-figs

    # The size of the system:
    N_sys = 3

    def proj(i, j):
        """ A helper function for creating an interaction operator. """
        return basis(N_sys, i) * basis(N_sys, j).dag()

    # Construct one bath for each system state:
    baths = []
    for i in range(N_sys):
        Q = proj(i, i)
        baths.append(DrudeLorentzBath(Q, lam, gamma, T, Nk))

    # Construct the system Liouvillian from the system Hamiltonian and
    # bath expansion terminators:
    H_sys = sum((i + 0.5) * eps * proj(i, i) for i in range(N_sys))
    H_sys += sum(
      (i + j + 0.5) * Del * proj(i, j)
      for i in range(N_sys) for j in range(N_sys)
      if i != j
    )
    HL = liouvillian(H_sys) + sum(bath.terminator()[1] for bath in baths)

    # Construct the solver (pass a list of baths):
    solver = HEOMSolver(HL, baths, max_depth=max_depth, options=options)

    # Run the solver:
    rho0 = basis(N_sys, 0) * basis(N_sys, 0).dag()
    tlist = np.linspace(0, 5, 200)
    e_ops = {
        f"P{i}": proj(i, i)
        for i in range(N_sys)
    }
    result = solver.run(rho0, tlist, e_ops=e_ops)

    # Plot populations:
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,8))
    for label, values in result.expect.items():
        axes.plot(result.times, values, label=label)
    axes.set_xlabel(r't', fontsize=28)
    axes.set_ylabel(r"Population", fontsize=28)
    axes.legend(loc=0, fontsize=12)


.. plot::
    :context: reset
    :include-source: false
    :nofigs:

    # reset the context at the end
########################
Previous implementations
########################

The current HEOM implementation in QuTiP is the latest in a succession of HEOM
implementations by various contributors:


HSolverDL
---------

The original HEOM solver was implemented by Neill Lambert, Anubhav Vardhan,
and Alexander Pitchford and is still available as
:class:`qutip.nonmarkov.dlheom_solver.HSolverDL` and only directly provided
support for the Drude-Lorentz bath although there was the possibility of
sub-classing the solver to implement other baths.

A compatible interface using the current implementation is available under the
same name in :class:`qutip.nonmarkov.heom.HSolverDL`.


BoFiN-HEOM
----------

BoFiN-HEOM (the bosonic and fermionic HEOM solver) was a much more
flexible re-write of the original QuTiP ``HSolverDL`` that added support for
both bosonic and fermionic baths and for baths to be specified directly via
their correlation function expansion coefficients. Its authors were
Neill Lambert, Tarun Raheja, Shahnawaz Ahmed, and Alexander Pitchford.

BoFiN was written outside of QuTiP and is can still be found in its original
repository at https://github.com/tehruhn/bofin.

The construction of the right-hand side matrix for BoFiN was slow, so
BoFiN-fast, a hybrid C++ and Python implementation, was written that performed
the right-hand side construction in C++. It was otherwise identical to the
pure Python version. BoFiN-fast can be found at
https://github.com/tehruhn/bofin_fast.

BoFiN also came with an extensive set of example notebooks that are available
at https://github.com/tehruhn/bofin/tree/main/examples.


Current implementation
----------------------

The current implementation is a rewrite of BoFiN in pure Python. It's
right-hand side construction has similar speed to BoFiN-fast, but is written
in pure Python. Built-in implementations of a variety of different baths
are provided, and a single solver is used for both fermionic and bosonic baths.
Multiple baths of the same kind (either fermionic or bosonic) may be
specified in a single problem, and there is good support for working with
the auxiliary density operator (ADO) state and extracting information from it.

The code was written by Neill Lambert and Simon Cross.
######################
Fermionic Environments
######################

Here we model a single fermion coupled to two electronic leads or reservoirs
(e.g.,  this can describe a single quantum dot, a molecular transistor, etc).
The system hamiltonian, :math:`H_{sys}`, and bath spectral density, :math:`J_D`,
are

.. math::

    H_{sys} &= c^{\dagger} c

    J_D &= \frac{\Gamma W^2}{(w - \mu)^2 + W^2},

We will demonstrate how to describe the bath using two different expansions
of the spectral density correlation function (Matsubara's expansion and
a Padé expansion), how to evolve the system in time, and how to calculate
the steady state.

Since our fermion is coupled to two reservoirs, we will construct two baths --
one for each reservoir or lead -- and call them the left (:math:`L`) and right
(:math:`R`) baths for convenience. Each bath will have a different chemical
potential :math:`\mu` which we will label :math:`\mu_L` and :math:`\mu_R`.

First we will do this using the built-in implementations of the bath expansions,
:class:`~qutip.nonmarkov.heom.LorentzianBath` and
:class:`~qutip.nonmarkov.heom.LorentzianPadeBath`.

Afterwards, we will show how to calculate the bath expansion coefficients and to
use those coefficients to construct your own bath description so that you can
implement your own fermionic baths.

Our implementation of fermionic baths primarily follows the definitions used by
Christian Schinabeck in his dissertation (
https://opus4.kobv.de/opus4-fau/files/10984/DissertationChristianSchinabeck.pdf
) and related publications.

A notebook containing a complete example similar to this one implemented in
BoFiN can be found in `example notebook 4b
<https://github.com/tehruhn/bofin/blob/main/examples/example-4b-fermions-single-impurity-model.ipynb>`__.


Describing the system and bath
------------------------------

First, let us construct the system Hamiltonian, :math:`H_{sys}`, and the initial
system state, ``rho0``:

.. plot::
    :context: reset
    :nofigs:

    from qutip import basis, destroy

    # The system Hamiltonian:
    e1 = 1.  # site energy
    H_sys = e1 * destroy(2).dag() * destroy(2)

    # Initial state of the system:
    rho0 = basis(2,0) * basis(2,0).dag()

Now let us describe the bath properties:

.. plot::
    :context:
    :nofigs:

    # Shared bath properties:
    gamma = 0.01   # coupling strength
    W = 1.0  # cut-off
    T = 0.025851991  # temperature
    beta = 1. / T

    # Chemical potentials for the two baths:
    mu_L = 1.
    mu_R = -1.

    # System-bath coupling operator:
    Q = destroy(2)

where :math:`\Gamma` (``gamma``), :math:`W` and :math:`T` are the parameters of
an Lorentzian bath, :math:`\mu_L` (``mu_L``) and :math:`\mu_R` (``mu_R``) are
the chemical potentials of the left and right baths, and ``Q`` is the coupling
operator between the system and the baths.

We may the pass these parameters to either ``LorentzianBath`` or
``LorentzianPadeBath`` to construct an expansion of the bath correlations:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import LorentzianBath
    from qutip.nonmarkov.heom import LorentzianPadeBath

    # Number of expansion terms to retain:
    Nk = 2

    # Matsubara expansion:
    bath_L = LorentzianBath(Q, gamma, W, mu_L, T, Nk, tag="L")
    bath_R = LorentzianBath(Q, gamma, W, mu_R, T, Nk, tag="R")

    # Padé expansion:
    bath_L = LorentzianPadeBath(Q, gamma, W, mu_L, T, Nk, tag="L")
    bath_R = LorentzianPadeBath(Q, gamma, W, mu_R, T, Nk, tag="R")

Where ``Nk`` is the number of terms to retain within the expansion of the
bath.

Note that we haved labelled each bath with a tag (either "L" or "R") so that
we can identify the exponents from individual baths later when calculating
the currents between the system and the bath.


System and bath dynamics
------------------------

Now we are ready to construct a solver:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import HEOMSolver
    from qutip import Options

    max_depth = 5  # maximum hierarchy depth to retain
    options = Options(nsteps=15_000)
    baths = [bath_L, bath_R]

    solver = HEOMSolver(H_sys, baths, max_depth=max_depth, options=options)

and to calculate the system evolution as a function of time:

.. code-block:: python

    tlist = [0, 10, 20]  # times to evaluate the system state at
    result = solver.run(rho0, tlist)

As in the bosonic case, the ``max_depth`` parameter determines how many levels
of the hierarchy to retain.

As in the bosonic case, we can specify ``e_ops`` in order to retrieve the
expectation values of operators at each given time. See
:ref:`heom-bosonic-system-and-bath-dynamics` for a fuller description of
the returned ``result`` object.

Below we run the solver again, but use ``e_ops`` to store the expectation
values of the population of the system states:

.. plot::
    :context:

    # Define the operators that measure the populations of the two
    # system states:
    P11p = basis(2,0) * basis(2,0).dag()
    P22p = basis(2,1) * basis(2,1).dag()

    # Run the solver:
    tlist = np.linspace(0, 500, 101)
    result = solver.run(rho0, tlist, e_ops={"11": P11p, "22": P22p})

    # Plot the results:
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,8))
    axes.plot(result.times, result.expect["11"], 'b', linewidth=2, label="P11")
    axes.plot(result.times, result.expect["22"], 'r', linewidth=2, label="P22")
    axes.set_xlabel(r't', fontsize=28)
    axes.legend(loc=0, fontsize=12)

The plot above is not very exciting. What we would really like to see in
this case are the currents between the system and the two baths. We will plot
these in the next section using the auxiliary density operators (ADOs)
returned by the solver.


.. _heom-determining-currents:

Determining currents
--------------------

The currents between the system and a fermionic bath may be calculated from the
first level auxiliary density operators (ADOs) associated with the exponents
of that bath.

The contribution to the current into a given bath from each exponent in that
bath is:

.. math::

    \mathrm{Contribution from Exponent} = \pm i \mathrm{Tr}(Q^\pm \cdot A)

where the :math:`\pm` sign is the sign of the exponent (see the
description later in :ref:`heom-fermionic-pade-expansion-coefficients`) and
:math:`Q^\pm` is :math:`Q` for ``+`` exponents and :math:`Q^{\dagger}` for
``-`` exponents.

The first-level exponents for the left bath are retrieved by calling
``.filter(tags=["L"])`` on ``ado_state`` which is an instance of
:class:`~qutip.nonmarkov.heom.HierarchyADOsState` and also provides access to
the methods of :class:`~qutip.nonmarkov.heom.HierarchyADOs` which describes the
structure of the hierarchy for a given problem.

Here the tag "L" matches the tag passed when constructing ``bath_L`` earlier
in this example.

Similarly, we may calculate the current to the right bath from the exponents
tagged with "R".

.. plot::
    :context:
    :nofigs:

    def exp_current(aux, exp):
        """ Calculate the current for a single exponent. """
        sign = 1 if exp.type == exp.types["+"] else -1
        op = exp.Q if exp.type == exp.types["+"] else exp.Q.dag()
        return 1j * sign * (op * aux).tr()

    def heom_current(tag, ado_state):
        """ Calculate the current between the system and the given bath. """
        level_1_ados = [
            (ado_state.extract(label), ado_state.exps(label)[0])
            for label in ado_state.filter(tags=[tag])
        ]
        return np.real(sum(exp_current(aux, exp) for aux, exp in level_1_ados))

    heom_left_current = lambda t, ado_state: heom_current("L", ado_state)
    heom_right_current = lambda t, ado_state: heom_current("R", ado_state)

Once we have defined functions for retrieving the currents for the
baths, we can pass them to ``e_ops`` and plot the results:

.. plot::
    :context: close-figs

    # Run the solver (returning ADO states):
    tlist = np.linspace(0, 100, 201)
    result = solver.run(rho0, tlist, e_ops={
        "left_currents": heom_left_current,
        "right_currents": heom_right_current,
    })

    # Plot the results:
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,8))
    axes.plot(
        result.times, result.expect["left_currents"], 'b',
        linewidth=2, label=r"Bath L",
    )
    axes.plot(
        result.times, result.expect["right_currents"], 'r',
        linewidth=2, label="Bath R",
    )
    axes.set_xlabel(r't', fontsize=28)
    axes.set_ylabel(r'Current', fontsize=20)
    axes.set_title(r'System to Bath Currents', fontsize=20)
    axes.legend(loc=0, fontsize=12)

And now we have a more interesting plot that shows the currents to the
left and right baths decaying towards their steady states!

In the next section, we will calculate the steady state currents directly.


Steady state currents
---------------------

Using the same solver, we can also determine the steady state of the
combined system and bath using:

.. plot::
    :context:
    :nofigs:

    steady_state, steady_ados = solver.steady_state()

and calculate the steady state currents to the two baths from ``steady_ados``
using the same ``heom_current`` function we defined previously:

.. plot::
    :context:
    :nofigs:

    steady_state_current_left = heom_current("L", steady_ados)
    steady_state_current_right = heom_current("R", steady_ados)

Now we can add the steady state currents to the previous plot:

.. plot::
    :context: close-figs

    # Plot the results and steady state currents:
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,8))
    axes.plot(
        result.times, result.expect["left_currents"], 'b',
        linewidth=2, label=r"Bath L",
    )
    axes.plot(
        result.times, [steady_state_current_left] * len(result.times), 'b:',
        linewidth=2, label=r"Bath L (steady state)",
    )
    axes.plot(
        result.times, result.expect["right_currents"], 'r',
        linewidth=2, label="Bath R",
    )
    axes.plot(
        result.times, [steady_state_current_right] * len(result.times), 'r:',
        linewidth=2, label=r"Bath R (steady state)",
    )
    axes.set_xlabel(r't', fontsize=28)
    axes.set_ylabel(r'Current', fontsize=20)
    axes.set_title(r'System to Bath Currents (with steady states)', fontsize=20)
    axes.legend(loc=0, fontsize=12)

As you can see, there is still some way to go beyond ``t = 100`` before the
steady state is reached!


.. _heom-fermionic-pade-expansion-coefficients:

Padé expansion coefficients
---------------------------

We now look at how to calculate the correlation expansion coefficients for the
Lorentzian spectral density ourselves. Once we have calculated the coefficients
we can construct a :class:`~qutip.nonmarkov.heom.FermionicBath` directly from
them. A similar procedure can be used to apply
:class:`~qutip.nonmarkov.heom.HEOMSolver` to any fermionic bath for which we can
calculate the expansion coefficients.

In the fermionic case we must descriminate between the order in which
excitations are created within the bath, so we define two different correlation
functions, :math:`C_{+}(t)`, and :math:`C_{-}(t)`:

.. math::

    C^{\sigma}(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} d\omega e^{\sigma i \omega t} J(\omega) f_F[\sigma\beta(\omega - \mu)]

where :math:`\sigma` is either ``+`` or ``-`` and, :math:`f_F` is the Fermi
distribution function, and :math:`J(\omega)` is the Lorentzian spectral density
we defined at the start.

The Fermi distribution function is:

.. math::

    f_F (x) = (\exp(x) + 1)^{-1}

As in the bosonic case we can approximate this integral with a Matsubara or
Padé expansion. For the Lorentzian bath the Padé expansion converges much
more quickly, so we will calculate the Padé expansion coefficients here.

The Padé decomposition approximates the Fermi distribution as:

.. math::

    f_F(x) \approx f_F^{\mathrm{approx}}(x) = \frac{1}{2} - \sum_{l=0}^{Nk} \frac{2k_l x}{x^2 + \epsilon_l^2}

where :math:`k_l` and :math:`\epsilon_l` are coefficients defined in
`J. Chem Phys 133, "Efficient on the fly calculation of time correlation functions in computer simulations" <https://doi.org/10.1063/1.3491098>`_,
and :math:`Nk` specifies the cut-off in the expansion.

Evaluating the integral for the correlation functions gives:

.. math::

    C^{\sigma}(t) \approx \sum_{l=0}^{Nk} \eta^{\sigma,l} e^{-\gamma_{\sigma,l}t}

where:

.. math::

    \eta_{\sigma, l} &= \begin{cases}
        \frac{\Gamma W}{2} f_F^{approx}(i\beta W)  & l = 0\\
        -i\cdot \frac{k_l}{\beta} \cdot \frac{\Gamma W^2}{-\frac{\epsilon^2_l}{\beta^2} + W^2}  & l \neq 0\\
    \end{cases}

    \gamma_{\sigma,l} &= \begin{cases}
        W - \sigma i\mu  & l = 0\\
        \frac{\epsilon_l}{\beta} - \sigma i \mu  & l \neq 0\\
    \end{cases}

and :math:`\beta = \frac{1}{T}`.

And now we calculate the same numbers in Python:

.. plot::
    :context:
    :nofigs:

    # Imports
    from numpy.linalg import eigvalsh

    # Convenience functions and parameters:
    def deltafun(j, k):
        """ Kronecker delta function. """
        return 1.0 if j == k else 0.

    def f_approx(x, Nk):
        """ Padé approxmation to Fermi distribution. """
        f = 0.5
        for ll in range(1, Nk + 1):
            # kappa and epsilon are calculated further down
            f = f - 2 * kappa[ll] * x / (x**2 + epsilon[ll]**2)
        return f

    def kappa_epsilon(Nk):
        """ Calculate kappa and epsilon coefficients. """

        alpha = np.zeros((2 * Nk, 2 * Nk))
        for j in range(2 * Nk):
            for k in range(2 * Nk):
                alpha[j][k] = (
                    (deltafun(j, k + 1) + deltafun(j, k - 1))
                    / np.sqrt((2 * (j + 1) - 1) * (2 * (k + 1) - 1))
                )

        eps = [-2. / val for val in eigvalsh(alpha)[:Nk]]

        alpha_p = np.zeros((2 * Nk - 1, 2 * Nk - 1))
        for j in range(2 * Nk - 1):
            for k in range(2 * Nk - 1):
                alpha_p[j][k] = (
                    (deltafun(j, k + 1) + deltafun(j, k - 1))
                    / np.sqrt((2 * (j + 1) + 1) * (2 * (k + 1) + 1))
                )

        chi = [-2. / val for val in eigvalsh(alpha_p)[:Nk - 1]]

        eta_list = [
            0.5 * Nk * (2 * (Nk + 1) - 1) * (
                np.prod([chi[k]**2 - eps[j]**2 for k in range(Nk - 1)]) /
                np.prod([
                    eps[k]**2 - eps[j]**2 + deltafun(j, k) for k in range(Nk)
                ])
            )
            for j in range(Nk)
        ]

        kappa = [0] + eta_list
        epsilon = [0] + eps

        return kappa, epsilon

    kappa, epsilon = kappa_epsilon(Nk)

    # Phew, we made it to function that calculates the coefficients for the
    # correlation function expansions:

    def C(sigma, mu, Nk):
        """ Calculate the expansion coefficients for C_\sigma. """
        beta = 1. / T
        ck = [0.5 * gamma * W * f_approx(1.0j * beta * W, Nk)]
        vk = [W - sigma * 1.0j * mu]
        for ll in range(1, Nk + 1):
            ck.append(
                -1.0j * (kappa[ll] / beta) * gamma * W**2
                / (-(epsilon[ll]**2 / beta**2) + W**2)
            )
            vk.append(epsilon[ll] / beta - sigma * 1.0j * mu)
        return ck, vk

    ck_plus_L, vk_plus_L = C(1.0, mu_L, Nk)  # C_+, left bath
    ck_minus_L, vk_minus_L = C(-1.0, mu_L, Nk)  # C_-, left bath

    ck_plus_R, vk_plus_R = C(1.0, mu_R, Nk)  # C_+, right bath
    ck_minus_R, vk_minus_R = C(-1.0, mu_R, Nk)  # C_-, right bath

Finally we are ready to construct the
:class:`~qutip.nonmarkov.heom.FermionicBath`:

.. plot::
    :context:
    :nofigs:

    from qutip.nonmarkov.heom import FermionicBath

    # Padé expansion:
    bath_L = FermionicBath(Q, ck_plus_L, vk_plus_L, ck_minus_L, vk_minus_L)
    bath_R = FermionicBath(Q, ck_plus_R, vk_plus_R, ck_minus_R, vk_minus_R)

And we're done!

The :class:`~qutip.nonmarkov.heom.FermionicBath` can be used with the
:class:`~qutip.nonmarkov.heom.HEOMSolver` in exactly the same way as the baths
we constructed previously using the built-in Lorentzian bath expansions.


.. plot::
    :context: reset
    :include-source: false
    :nofigs:

    # reset the context at the end
References
==========

.. bibliography:: heom.bib
    :all:
############
Introduction
############

The Hierarchical Equations of Motion (HEOM) method was originally developed by
Tanimura and Kubo :cite:`Tanimura_1989` in the context of physical chemistry to
''exactly'' solve a quantum system in contact with a bosonic environment,
encapsulated in the Hamiltonian:

.. math::

	H = H_s + \sum_k \omega_k a_k^{\dagger}a_k + \hat{Q} \sum_k g_k \left(a_k + a_k^{\dagger}\right).

As in other solutions to this problem, the properties of the bath are
encapsulated by its temperature and its spectral density,

.. math::

    J(\omega) = \pi \sum_k g_k^2 \delta(\omega-\omega_k).

In the HEOM, for bosonic baths, one typically chooses a Drude-Lorentz spectral
density:

.. math::

    J_D = \frac{2\lambda \gamma \omega}{(\gamma^2 + \omega^2)},

or an under-damped Brownian motion spectral density:

.. math::

    J_U = \frac{\alpha^2 \Gamma \omega}{[(\omega_c^2 - \omega^2)^2 + \Gamma^2 \omega^2]}.

Given the spectral density, the HEOM requires a decomposition of the bath
correlation functions in terms of exponentials. In :doc:`bosonic` we describe
how this is done with code examples, and how these expansions are passed to the
solver.

In addition to support for bosonic environments, QuTiP also provides support for
feriomic environments which is described in :doc:`fermionic`.

Both bosonic and fermionic environments are supported via a single solver,
:class:`HEOMSolver`, that supports solving for both dynamics and steady-states.
.. _qip_processor:

******************************
Pulse-level circuit simulation
******************************

Modelling quantum hardware with Processor
-----------------------------------------

Based on the open system solver, :class:`~qutip.qip.device.Processor` in the :mod:`qutip.qip` module simulates quantum circuits at the level of time evolution. One can consider the processor as a simulator of a quantum device, on which the quantum circuit is to be implemented. 

The procedure is illustrated in the figure below.
It first compiles circuit into a Hamiltonian model, adds noisy dynamics and then uses the QuTiP open time evolution solvers to simulation the evolution.

.. image:: /figures/qip/illustration.png

Like a real quantum device, the processor is determined by a list of Hamiltonians, i.e. the control pulses driving the evolution. Given the intensity of the control pulses and the corresponding time slices for each pulse, the evolution is then computed. A control pulse is characterized by :class:`~qutip.qip.pulse.Pulse`, consisting of the control Hamiltonian, the targets qubit, the pulse coefficients and the time sequence. We can either use the coefficients as a step function or with cubic spline. For step function, ``tlist`` specifies the start and the end of each pulse and thus is one element longer the ``coeffs``. One example of defining the control pulse coefficients and the time array is as follows:

.. testcode::

    import numpy as np
    from qutip import sigmaz
    from qutip.qip.device import Processor

    processor = Processor(2)
    processor.add_control(sigmaz(), cyclic_permutation=True)  # sigmaz for all qubits
    processor.pulses[0].coeffs = np.array([[1.0, 1.5, 2.0], [1.8, 1.3, 0.8]])
    processor.pulses[0].tlist = np.array([0.1, 0.2, 0.4, 0.5])

It defines a :math:`\sigma_z` operator on both qubits and a pulse that acts on the first qubit.
An equivalent approach is using the :meth:`~qutip.qip.device.Processor.add_pulse` method.

.. testcode::

    from qutip.qip.pulse import Pulse

    processor = Processor(2)
    coeff=np.array([0.1, 0.2, 0.4, 0.5])
    tlist=np.array([[1.0, 1.5, 2.0], [1.8, 1.3, 0.8]])
    pulse = Pulse(sigmaz(), targets=0, coeff=coeff, tlist=tlist)
    processor.add_pulse(pulse)

One can also use choose the ``pulse_mode`` attribute of :class:`~qutip.qip.device.Processor`
between ``"discrete"`` and ``"continuous"``.

.. note::

   If the coefficients represent dicrete pulse, the length of each array is 1 element shorter than ``tlist``. If it is supposed to be a continuous function, the length should be the same as ``tlist``.


The above example shows the framework and the most essential part of the simulator's API. So far, it looks like just a wrapper for the open system solvers. However, based on this, we can implement different physical realizations. They differ mainly in how to find the control pulse for a quantum circuit, which gives birth to different sub-classes:

| Processor
| ├── ModelProcessor
| │   ├── DispersiveCavityQED
| │   └── SpinChain
| └── OptPulseProcessor

In general, there are two ways to find the control pulses. The first one, :class:`~qutip.qip.device.ModelProcessor`, is more experiment-oriented and based on physical models. A universal set of
gates is defined in the processor as well as the pulse implementing them in this particular physical model. This is usually the case where control pulses realizing those gates are well known and can be concatenated to realize the whole quantum circuits. Two realizations have already been implemented: the spin chain and the Cavity QED model for quantum computing. In those models, the driving Hamiltonians are predefined. Another approach, based on the optimal control module in QuTiP (see :ref:`control`), is called :class:`~qutip.qip.device.OptPulseProcessor`. In this subclass, one only defines the available Hamiltonians in their system. The processor then uses algorithms to find the optimal control pulses that realize the desired unitary evolution.

Despite this difference, the logic behind all processors is the same:

* One defines a processor by a list of available Hamiltonians and, as explained later, hardware-dependent noise. In model based processors, the Hamiltonians are predefined and one only needs to give the device parameters like frequency and interaction strength.

* The control pulse coefficients and time slices are either specified by the user or calculated by the method :meth:`~qutip.qip.device.Processor.load_circuit`, which takes a :class:`~qutip.qip.circuit.QubitCircuit` and find the control pulse for this evolution.

* The processor calculates the evolution using the QuTiP solvers. Collapse operators can be added to simulate decoherence. The method :meth:`~qutip.qip.device.Processor.run_state` returns a object :class:`qutip.solver.Result`.

It is also possible to calculate the evolution analytically with matrix exponentiation by setting ``analytical=True``. A list of the matrices representing the gates is returned just like for :meth:`~qutip.qip.circuit.QubitCircuit.propagators`. However, this does not consider the collapse operators or other noise. As the system size gets larger, this approach will become very inefficient.

In the following we describe the predefined subclasses for :class:`~qutip.qip.device.Processor`:

**SpinChain**

:class:`~qutip.qip.device.LinearSpinChain` and :class:`~qutip.qip.device.CircularSpinChain` are quantum computing models base on the spin chain realization. The control Hamiltonians are :math:`\sigma_x`, :math:`\sigma_z` and :math:`\sigma_x \sigma_x + \sigma_y \sigma_y`. This processor will first decompose the gate into the universal gate set with ISWAP or SQRTISWAP as two-qubit gates, resolve them into quantum gates of adjacent qubits and then calculate the pulse coefficients.

An example of simulating a simple circuit is shown below:

.. testcode::

    from qutip import basis
    from qutip.qip.circuit import QubitCircuit
    from qutip.qip.device import LinearSpinChain

    qc = QubitCircuit(2)
    qc.add_gate("X", targets=0)
    qc.add_gate("X", targets=1)
    processor = LinearSpinChain(2)
    processor.load_circuit(qc)
    result = processor.run_state(basis([2,2], [0,0]))
    print(result.states[-1].tidyup(1.0e-6))

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket
    Qobj data =
    [[ 0.]
    [ 0.]
    [ 0.]
    [-1.]]

We can also visualize the pulses implementing this circuit:

.. plot::

    from qutip import basis
    from qutip.qip.circuit import QubitCircuit
    from qutip.qip.device import LinearSpinChain

    qc = QubitCircuit(2)
    qc.add_gate("X", targets=0)
    qc.add_gate("X", targets=1)
    processor = LinearSpinChain(2)
    processor.load_circuit(qc)
    fig, axis = processor.plot_pulses()
    fig.show()

**DispersiveCavityQED**

Same as above, :class:`~qutip.qip.device.DispersiveCavityQED` is a simulator based on Cavity Quantum Electrodynamics. The workflow is similar to the one for the spin chain, except that the component systems are a multi-level cavity and a qubits system. The control Hamiltonians are the single-qubit rotation together with the qubits-cavity interaction :math:`a^{\dagger} \sigma^{-} + a \sigma^{+}`. The device parameters including the cavity frequency, qubits frequency, detuning and interaction strength etc.

.. note::

   The :meth:`~qutip.qip.device.DispersiveCavityQED.run_state` method of :class:`~qutip.qip.device.DispersiveCavityQED`
   returns the full simulation result of the solver,
   hence including the cavity.
   To obtain the circuit result, one needs to first trace out the cavity state.

**OptPulseProcessor**

The :class:`~qutip.qip.device.OptPulseProcessor` uses the function in :func:`~qutip.control.pulseoptim.optimize_pulse_unitary` in the optimal control module to find the control pulses. The Hamiltonian includes a drift part and a control part and only the control part will be optimized. The unitary evolution follows

.. math::

   U(\Delta t)=\exp(\rm{i} \cdot \Delta t [H_d  + \sum_j u_j H_j] )

To let it find the optimal pulses, we need to give the parameters for :func:`~qutip.control.pulseoptim.optimize_pulse_unitary` as keyword arguments to :meth:`~qutip.qip.device.OptPulseProcessor.load_circuit`. Usually, the minimal requirements are the evolution time ``evo_time`` and the number of time slices ``num_tslots`` for each gate. Other parameters can also be given in the keyword arguments. For available choices, see :func:`~qutip.control.pulseoptim.optimize_pulse_unitary`. It is also possible to specify different parameters for different gates, as shown in the following example:

.. testcode::

      from qutip.qip.device import OptPulseProcessor
      from qutip.operators import sigmaz, sigmax, sigmay
      from qutip.tensor import tensor

      # Same parameter for all the gates
      qc = QubitCircuit(N=1)
      qc.add_gate("SNOT", 0)

      num_tslots = 10
      evo_time = 10
      processor = OptPulseProcessor(N=1, drift=sigmaz())
      processor.add_control(sigmax())
      # num_tslots and evo_time are two keyword arguments
      tlist, coeffs = processor.load_circuit(
      qc, num_tslots=num_tslots, evo_time=evo_time)

      # Different parameters for different gates
      qc = QubitCircuit(N=2)
      qc.add_gate("SNOT", 0)
      qc.add_gate("SWAP", targets=[0, 1])
      qc.add_gate('CNOT', controls=1, targets=[0])

      processor = OptPulseProcessor(N=2, drift=tensor([sigmaz()]*2))
      processor.add_control(sigmax(), cyclic_permutation=True)
      processor.add_control(sigmay(), cyclic_permutation=True)
      processor.add_control(tensor([sigmay(), sigmay()]))

      setting_args = {"SNOT": {"num_tslots": 10, "evo_time": 1},
                      "SWAP": {"num_tslots": 30, "evo_time": 3},
                      "CNOT": {"num_tslots": 30, "evo_time": 3}}

      tlist, coeffs = processor.load_circuit(
                      qc, setting_args=setting_args, merge_gates=False)

Compiler and scheduler
----------------------

.. note::

   New in QuTiP 4.6

In order to simulate quantum circuits at the level of time evolution.
We need to first compile the circuit into the Hamiltonian model, i.e.
the control pulses.
Hence each :class:`~qutip.qip.device.Processor` has a corresponding 
:class:`~qutip.qip.compiler.GateCompiler` class.
The compiler takes a :class:`~qutip.qip.circuit.QubitCircuit`
and returns the compiled ``tlist`` and ``coeffs``.
It is called implicitly when calling the method
:class:`~qutip.qip.device.Processor.run_state`.

.. testcode::

    from qutip.qip.compiler import SpinChainCompiler
    qc = QubitCircuit(2)
    qc.add_gate("X", targets=0)
    qc.add_gate("X", targets=1)

    processor = LinearSpinChain(2)
    compiler = SpinChainCompiler(
        2, params=processor.params, pulse_dict=processor.pulse_dict)
    resolved_qc = qc.resolve_gates(["RX", "RZ", "ISWAP"])
    tlists, coeffs = compiler.compile(resolved_qc)
    print(tlists)
    print(coeffs)

**Output**

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    [array([0., 1.]), array([0., 1., 2.]), None, None, None]
    [array([1.57079633]), array([0.        , 1.57079633]), None, None, None]

Here we first use :meth:`~qutip.qip.circuit.QubitCircuit.resolve_gates`
to decompose the X gate to its natural gate on Spin Chain model,
the rotation over X-axis.
We pass the hardware parameters of the :class:`~qutip.qip.device.SpinChain `` model, ``processor.params``, as well as a map between the pulse name and pulse index ``pulse_dict`` to the compiler.
The later one allows one to address the pulse more conveniently in the compiler.

The compiler returns a list of ``tlist`` and ``coeff``, corresponding to each pulse.
The first pulse starts from ``t=0`` and ends at ``t=1``, with the strengh :math:`\pi/2`.
The second one is turned on from ``t=1`` to ``t=2`` with the same strength.
The compiled pulse here is different from what is shown in the plot
in the previous subsection because the scheduler is turned off by default.

The scheduler is implemented in the class :class:`~qutip.qip.compiler.Scheduler`,
based on the idea of https://doi.org/10.1117/12.666419.
It schedules the order of quantum gates and instructions for the
shortest execution time.
It works not only for quantum gates but also for pulse implementation of gates
(:class:`~qutip.qip.compiler.Instruction`) with varying pulse duration.

The scheduler first generates a quantum gates dependency graph,
containing information about which gates have to be executed before some other gates.
The graph preserves the mobility of the gates,
i.e. commuting gates are not dependent on each other, even if they use the same qubits.
Next, it computes the longest distance of each node to the start and end nodes.
The distance for each dependency arrow is defined by the execution time of the instruction
(By default, it is 1 for all gates).
This is used as a priority measure in the next step.
The gate with a longer distance to the end node and a shorter distance to the start node has higher priority.
In the last step, it uses a list-schedule algorithm with hardware constraint and
priority and returns a list of cycles for gates/instructions.
Since the algorithm is heuristics, sometimes it does not find the optimal solution.
Hence, we offer an option that randomly shuffles the commuting gates and
repeats the scheduling a few times to get a better result.

.. testcode::

    from qutip.qip.circuit import QubitCircuit
    from qutip.qip.compiler import Scheduler
    circuit = QubitCircuit(7)
    circuit.add_gate("SNOT", 3)  # gate0
    circuit.add_gate("CZ", 5, 3)  # gate1
    circuit.add_gate("CZ", 4, 3)  # gate2
    circuit.add_gate("CZ", 2, 3)  # gate3
    circuit.add_gate("CZ", 6, 5)  # gate4
    circuit.add_gate("CZ", 2, 6)  # gate5
    circuit.add_gate("ISWAP", [0, 2])  # gate6
    scheduler = Scheduler("ASAP")
    result = scheduler.schedule(circuit, gates_schedule=True)
    print(result)

**Output**

.. testoutput::

    [0, 1, 3, 2, 2, 3, 4]

The result shows the scheduling order of each gate in the original circuit.

For pulse schedule, or scheduling gates with different duration,
one will need to wrap the :class:`qutip.qip.circuit.Gate` object with :class:`qutip.qip.compiler.instruction` object,
with a parameter `duration`.
The result will then be the start time of each instruction.

Noise Simulation
----------------

In the common way of QIP simulation, where evolution is carried out by gate matrix product, the noise is usually simulated with bit flipping and sign flipping errors.
The typical approaches are either applying bit/sign flipping gate probabilistically
or applying Kraus operators representing different noisy channels (e.g. amplitude damping, dephasing) after each unitary gate evolution. In the case of a single qubit, they have the same effect and the parameters in the Kraus operators are exactly the probability of a flipping error happens during the gate operation time.

Since the processor simulates the state evolution at the level of the driving Hamiltonian, there is no way to apply an error operator to the continuous-time evolution. Instead, the error is added to the pulses (coherent control error) or the collapse operators (Lindblad error) contributing to the evolution. Mathematically, this is no different from adding error channel probabilistically (it is actually how :func:`qutip.mcsolve` works internally). The collapse operator for single-qubit amplitude damping and dephasing are exactly the destroying operator and the sign-flipping operator. One just needs to choose the correct coefficients for them to simulate the noise, e.g. the relaxation time T1 and dephasing time T2. Because it is based on the open system evolution instead of abstract operators, this simulation is closer to the physical implementation and requires less pre-analysis of the system.

Compared to the approach of Kraus operators, this way of simulating noise is more computationally expensive. If you only want to simulate the decoherence of single-qubit relaxation and the relaxation time is much longer than the gate duration, there is no need to go through all the calculations. However, this simulator is closer to the real experiment and, therefore, more convenient in some cases, such as when coherent noise or correlated noise exist. For instance, a pulse on one qubit might affect the neighbouring qubits, the evolution is still unitary but the gate fidelity will decrease. It is not always easy or even possible to define a noisy gate matrix. In our simulator, it can be done by defining a :class:`~qutip.qip.noise.ControlAmpNoise` (Control Amplitude Noise).

In the simulation, noise can be added to the processor at different levels:

* The decoherence time T1 and T2 can be defined for the processor or for each qubit. When calculating the evolution, the corresponding collapse operators will be added automatically to the solver.

* The noise of the physical parameters (e.g. detuned frequency) can be simulated by changing the parameters in the model, e.g. laser frequency in cavity QED. (This can only be time-independent since QuTiP open system solver only allows varying coefficients, not varying Hamiltonian operators.)

* The noise of the pulse intensity can be simulated by modifying the coefficients of the Hamiltonian operators or even adding new Hamiltonians.

To add noise to a processor, one needs to first define a noise object :class:`~qutip.qip.noise.Noise`. The simplest relaxation noise can be defined directly in the processor with relaxation time. Other pre-defined noise can be found as subclasses of  :class:`~qutip.qip.noise.Noise`. We can add noise to the simulator with the method :meth:`~qutip.qip.device.Processor.add_noise`.

Below, we show two examples.

The first example is a processor with one qubit under rotation around the z-axis and relaxation time :math:`T_2=5`. We measure the population of the :math:`\left| + \right\rangle` state and observe the Ramsey signal:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from qutip import sigmaz, destroy, basis
    from qutip.qip.device import Processor
    from qutip.qip.operations import snot

    a = destroy(2)
    Hadamard = snot()
    plus_state = (basis(2,1) + basis(2,0)).unit()
    tlist = np.arange(0.00, 20.2, 0.2)

    T2 = 5
    processor = Processor(1, t2=T2)
    processor.add_control(sigmaz())
    processor.pulses[0].coeff = np.ones(len(tlist))
    processor.pulses[0].tlist = tlist
    result = processor.run_state(
        plus_state, e_ops=[a.dag()*a, Hadamard*a.dag()*a*Hadamard])

    fig, ax = plt.subplots()
    # detail about length of tlist needs to be fixed
    ax.plot(tlist[:-1], result.expect[1][:-1], '.', label="simulation")
    ax.plot(tlist[:-1], np.exp(-1./T2*tlist[:-1])*0.5 + 0.5, label="theory")
    ax.set_xlabel("t")
    ax.set_ylabel("Ramsey signal")
    ax.legend()
    ax.set_title("Relaxation T2=5")
    ax.grid()
    fig.tight_layout()
    fig.show()

The second example demonstrates a biased Gaussian noise on the pulse amplitude. For visualization purposes, we plot the noisy pulse intensity instead of the state fidelity. The three pulses can, for example, be a zyz-decomposition of an arbitrary single-qubit gate:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from qutip import sigmaz, sigmay
    from qutip.qip.device import Processor
    from qutip.qip.noise import RandomNoise

    # add control Hamiltonians
    processor = Processor(N=1)
    processor.add_control(sigmaz(), targets=0)

    # define pulse coefficients and tlist for all pulses
    processor.pulses[0].coeff = np.array([0.3, 0.5, 0. ])
    processor.set_all_tlist(np.array([0., np.pi/2., 2*np.pi/2, 3*np.pi/2]))

    # define noise, loc and scale are keyword arguments for np.random.normal
    gaussnoise = RandomNoise(
                dt=0.01, rand_gen=np.random.normal, loc=0.00, scale=0.02)
    processor.add_noise(gaussnoise)

    # Plot the ideal pulse
    fig1, axis1 = processor.plot_pulses(title="Original control amplitude", figsize=(5,3))

    # Plot the noisy pulse
    qobjevo, _ = processor.get_qobjevo(noisy=True)
    noisy_coeff = qobjevo.to_list()[1][1] + qobjevo.to_list()[2][1]
    fig2, axis2 = processor.plot_pulses(title="Noisy control amplitude", figsize=(5,3))
    axis2[0].step(qobjevo.tlist, noisy_coeff)


Customize the simulator
-----------------------
The number of predefined physical models and compilers are limited.
However, it is designed for easy customization and one can easily build customized model and compiling routines.
For guide and examples, please refer to the tutorial notebooks
at https://qutip.org/tutorials.html

The workflow of the simulator
-------------------------------
The following plot demonstrates the workflow of the simulator.

.. image:: /figures/qip/workflow.png

The core of the simulator is :class:`~qutip.qip.device.Processor`,
which characterizes the quantum hardware of interest,
containing the information such as the non-controllable drift Hamiltonian and
the control Hamiltonian.
Apart from the ideal system representing the qubits, one can also define
hardware-dependent or pulse-dependent noise in :class:`~qutip.qip.noise.Noise`.
It describes how noisy terms such as imperfect control
and decoherence can be added once the ideal control pulse is defined.
When loading a quantum circuit, a :class:`~qutip.qip.compiler.GateCompiler` compiles the circuit into a sequence of control pulse signals and schedule the pulse for parallel execution.
For each control Hamiltonian, a :class:`~qutip.qip.pulse.Pulse` instance is created that including the ideal evolution and associated noise.
They will then be sent to the QuTiP solvers for the computation.
.. _qip_simulator:

*********************************
Operator-level circuit simulation
*********************************

.. note::

   New in QuTiP 4.6

Run a quantum circuit
---------------------

Let's start off by defining a simple circuit which we use to demonstrate a few
examples of circuit evolution.
We take `a circuit from OpenQASM 2 <https://github.com/Qiskit/openqasm/blob/OpenQASM2.x/examples/W-state.qasm>`_

.. testcode::

    from qutip.qip.circuit import QubitCircuit, Gate
    from qutip.qip.operations import controlled_gate, hadamard_transform
    def controlled_hadamard():
        # Controlled Hadamard
        return controlled_gate(
            hadamard_transform(1), 2, control=0, target=1, control_value=1)
    qc = QubitCircuit(N=3, num_cbits=3)
    qc.user_gates = {"cH": controlled_hadamard}
    qc.add_gate("QASMU", targets=[0], arg_value=[1.91063, 0, 0])
    qc.add_gate("cH", targets=[0,1])
    qc.add_gate("TOFFOLI", targets=[2], controls=[0, 1])
    qc.add_gate("X", targets=[0])
    qc.add_gate("X", targets=[1])
    qc.add_gate("CNOT", targets=[1], controls=0)

It corresponds to the following circuit:

.. image:: /figures/qip/quantum_circuit_w_state.png

We will add the measurement gates later. This circuit prepares the W-state :math:`(\ket{001} + \ket{010} + \ket{100})/\sqrt{3}`.
The simplest way to carry out state evolution through a quantum circuit is
providing a input state to the :meth:`~qutip.qip.circuit.QubitCircuit.run`
method.

.. testcode::

  from qutip import tensor
  zero_state = tensor(basis(2, 0), basis(2, 0), basis(2, 0))
  result = qc.run(state=zero_state)
  wstate = result

  print(wstate)

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
  Qobj data =
  [[0.        ]
   [0.57734961]
   [0.57734961]
   [0.        ]
   [0.57735159]
   [0.        ]
   [0.        ]
   [0.        ]]


As expected, the state returned is indeed the required W-state.

As soon as we introduce measurements into the circuit, it can lead to multiple outcomes
with associated probabilities.  We can also carry out circuit evolution in a manner such that it returns all the possible state
outputs along with their corresponding probabilities. Suppose, in the previous circuit,
we measure each of the three qubits at the end.

.. testcode::

  qc.add_measurement("M0", targets=[0], classical_store=0)
  qc.add_measurement("M1", targets=[1], classical_store=1)
  qc.add_measurement("M2", targets=[2], classical_store=2)

To get all the possible output states along with the respective probability of observing the
outputs, we can use the :meth:`~qutip.qip.circuit.QubitCircuit.run_statistics` function:

.. testcode::

    result = qc.run_statistics(state=tensor(basis(2, 0), basis(2, 0), basis(2, 0)))
    states = result.get_final_states()
    probabilities = result.get_probabilities()

    for state, probability in zip(states, probabilities):
        print("State:\n{}\nwith probability {}".format(state, probability))

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

    State:
    Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
    Qobj data =
    [[0.]
    [1.]
    [0.]
    [0.]
    [0.]
    [0.]
    [0.]
    [0.]]
    with probability 0.33333257054168813
    State:
    Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
    Qobj data =
    [[0.]
    [0.]
    [1.]
    [0.]
    [0.]
    [0.]
    [0.]
    [0.]]
    with probability 0.33333257054168813
    State:
    Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
    Qobj data =
    [[0.]
    [0.]
    [0.]
    [0.]
    [1.]
    [0.]
    [0.]
    [0.]]
    with probability 0.33333485891662384

The function returns a :class:`~qutip.qip.Result` object which contains
the output states.
The method :meth:`~qutip.qip.Result.get_results` can be used to obtain the
possible states and probabilities.
Since the state created by the circuit is the W-state, we observe the states
:math:`\ket{001}`,  :math:`\ket{010}` and :math:`\ket{100}` with equal probability.


Circuit simulator
-----------------

.. _simulator_class:

The :meth:`~qutip.qip.circuit.QubitCircuit.run` and :meth:`~qutip.qip.circuit.QubitCircuit.run_statistics` functions
make use of the :class:`~qutip.qip.circuit.CircuitSimulator` which enables exact simulation with more
granular options. The simulator object takes a quantum circuit as an argument. It can optionally
be supplied with an initial state. There are two modes in which the exact simulator can function. The default mode is the
"state_vector_simulator" mode. In this mode, the state evolution proceeds maintaining the ket state throughout the computation.
For each measurement gate, one of the possible outcomes is chosen probabilistically
and computation proceeds. To demonstrate, we continue with our previous circuit:


.. testcode::

  from qutip.qip.circuit import CircuitSimulator

  sim = CircuitSimulator(qc, state=zero_state)

This initializes the simulator object and carries out any pre-computation
required. There are two ways to carry out state evolution with the simulator.
The primary way is to use the :meth:`~qutip.qip.circuit.CircuitSimulator.run` and
:meth:`~qutip.qip.circuit.CircuitSimulator.run_statistics` functions just like before (only
now with the :class:`~qutip.qip.circuit.CircuitSimulator` class).

The :class:`~qutip.qip.circuit.CircuitSimulator` class also enables stepping through the circuit:

.. testcode::

  print(sim.step())

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2, 2, 2], [1, 1, 1]], shape = (8, 1), type = ket
  Qobj data =
  [[0.57735159]
   [0.        ]
   [0.        ]
   [0.        ]
   [0.81649565]
   [0.        ]
   [0.        ]
   [0.        ]]

This only excutes one gate in the circuit and
allows for a better understanding of how the state evolution takes place.
The method steps through both the gates and the measurements.

Precomputing the unitary
------------------------

By default, the :class:`~qutip.qip.circuit.CircuitSimulator` class is initialized such that
the circuit evolution is conducted by applying each unitary to the state interactively.
However, by setting the argument ``precompute_unitary=True``, :class:`~qutip.qip.circuit.CircuitSimulator`
precomputes the product of the unitaries (in between the measurements):

.. testcode::

  sim = CircuitSimulator(qc, precompute_unitary=True)

  print(sim.ops)

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  [Quantum object: dims = [[2, 2, 2], [2, 2, 2]], shape = (8, 8), type = oper, isherm = False
    Qobj data =
    [[ 0.          0.57734961  0.         -0.57734961  0.          0.40824922
       0.         -0.40824922]
     [ 0.57734961  0.         -0.57734961  0.          0.40824922  0.
      -0.40824922  0.        ]
     [ 0.57734961  0.          0.57734961  0.          0.40824922  0.
       0.40824922  0.        ]
     [ 0.          0.57734961  0.          0.57734961  0.          0.40824922
       0.          0.40824922]
     [ 0.57735159  0.          0.          0.         -0.81649565  0.
       0.          0.        ]
     [ 0.          0.57735159  0.          0.          0.         -0.81649565
       0.          0.        ]
     [ 0.          0.          0.57735159  0.          0.          0.
      -0.81649565  0.        ]
     [ 0.          0.          0.          0.57735159  0.          0.
       0.         -0.81649565]],
       Measurement(M0, target=[0], classical_store=0),
       Measurement(M1, target=[1], classical_store=1),
       Measurement(M2, target=[2], classical_store=2)]


Here, ``sim.ops`` stores all the circuit operations that are going to be applied during
state evolution. As observed above, all the unitaries of the circuit are compressed into
a single unitary product with the precompute optimization enabled.
This is more efficient if one runs the same circuit one multiple initial states.
However, as the number of qubits increases, this will consume more and more memory
and become unfeasible.

Density Matrix Simulation
-------------------------

By default, the state evolution is carried out in the "state_vector_simulator" mode
(specified by the **mode** argument) as described before.
In the "density_matrix_simulator" mode, the input state can be either a ket or a density
matrix. If it is a ket, it is converted into a density matrix before the evolution is
carried out. Unlike the "state_vector_simulator" mode, upon measurement, the state
does not collapse to one of the post-measurement states. Rather, the new state is now
the density matrix representing the ensemble of post-measurement states.
In this sense, we measure the qubits and forget all the results.

To demonstrate this consider the original W-state preparation circuit which is followed
just by measurement on the first qubit:

.. testcode::

    qc = QubitCircuit(N=3, num_cbits=3)
    qc.user_gates = {"cH": controlled_hadamard}
    qc.add_gate("QASMU", targets=[0], arg_value=[1.91063, 0, 0])
    qc.add_gate("cH", targets=[0,1])
    qc.add_gate("TOFFOLI", targets=[2], controls=[0, 1])
    qc.add_gate("X", targets=[0])
    qc.add_gate("X", targets=[1])
    qc.add_gate("CNOT", targets=[1], controls=0)
    qc.add_measurement("M0", targets=[0], classical_store=0)
    qc.add_measurement("M0", targets=[1], classical_store=0)
    qc.add_measurement("M0", targets=[2], classical_store=0)
    sim = CircuitSimulator(qc, mode="density_matrix_simulator")
    print(sim.run(zero_state).get_final_states()[0])

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

    Quantum object: dims = [[2, 2, 2], [2, 2, 2]], shape = (8, 8), type = oper, isherm = True
    Qobj data =
    [[0.         0.         0.         0.         0.         0.
      0.         0.        ]
     [0.         0.33333257 0.         0.         0.         0.
      0.         0.        ]
     [0.         0.         0.33333257 0.         0.         0.
      0.         0.        ]
     [0.         0.         0.         0.         0.         0.
      0.         0.        ]
     [0.         0.         0.         0.         0.33333486 0.
      0.         0.        ]
     [0.         0.         0.         0.         0.         0.
      0.         0.        ]
     [0.         0.         0.         0.         0.         0.
      0.         0.        ]
     [0.         0.         0.         0.         0.         0.
      0.         0.        ]]

We are left with a mixed state.

Import and export quantum circuits
----------------------------------

QuTiP supports importation and exportation of quantum circuit in the `OpenQASM 2 format <https://github.com/Qiskit/openqasm/tree/OpenQASM2.x>`_
through the functions :func:`~qutip.qip.qasm.read_qasm` and :func:`~qutip.qip.qasm.save_qasm`.
We demonstrate this using the w-state generation circuit.
The following code is in OpenQASM format:

.. code-block::

    // Name of Experiment: W-state v1

    OPENQASM 2.0;
    include "qelib1.inc";


    qreg q[4];
    creg c[3];
    gate cH a,b {
    h b;
    sdg b;
    cx a,b;
    h b;
    t b;
    cx a,b;
    t b;
    h b;
    s b;
    x b;
    s a;
    }

    u3(1.91063,0,0) q[0];
    cH q[0],q[1];
    ccx q[0],q[1],q[2];
    x q[0];
    x q[1];
    cx q[0],q[1];

    measure q[0] -> c[0];
    measure q[1] -> c[1];
    measure q[2] -> c[2];

One can save it in a ``.qasm`` file and import it using the following code:

.. testcode::

  from qutip.qip.qasm import read_qasm
  qc = read_qasm("guide/qip/w-state.qasm")
.. _qip_intro:

******************************
Quantum Information Processing
******************************

Introduction
============

The Quantum Information Processing (QIP) module aims at providing basic tools for quantum computing simulation both for simple quantum algorithm design and for experimental realization. It offers two different approaches, one with :class:`~qutip.qip.QubitCircuit` calculating unitary evolution under quantum gates by matrix product, another called :class:`~qutip.qip.device.Processor` using open system solvers in QuTiP to simulate noisy quantum device.

.. _quantum_circuits:

Quantum Circuit
===============

The most common model for quantum computing is the quantum circuit model.
In QuTiP, we use :class:`~qutip.qip.QubitCircuit` to represent a quantum circuit.
The circuit is characterized by registers and gates:

- **Registers**: The argument ``N`` specifies the number of qubit registers in the circuit
  and the argument ``num_cbits`` (optional) specifies the number of classical bits available for measurement
  and control.

- **Gates**: Each quantum gate is saved as a class object :class:`~qutip.qip.Gate`
  with information such as gate name, target qubits and arguments.
  Gates can also be controlled on a classical bit by specifying the register number
  with the argument ``classical_controls``.

- **Measurements**: We can also carry out measurements on individual qubit (both in the middle and at the end of the circuit).
  Each measurement is saved as a class object :class:`~qutip.qip.Measurement` with parameters such as `targets`,
  the target qubit on which the measurement will be carried out, and `classical_store`,
  the index of the classical register which stores the result of the measurement.

A circuit with the various gates and registers available is demonstrated below:

.. testcode::

  from qutip.qip.circuit import QubitCircuit, Gate
  from qutip import tensor, basis

  qc = QubitCircuit(N=2, num_cbits=1)
  swap_gate = Gate(name="SWAP", targets=[0, 1])

  qc.add_gate(swap_gate)
  qc.add_measurement("M0", targets=[1], classical_store=0) # measurement gate
  qc.add_gate("CNOT", controls=0, targets=1)
  qc.add_gate("X", targets=0, classical_controls=[0]) # classically controlled gate
  qc.add_gate(swap_gate)

  print(qc.gates)

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  [Gate(SWAP, targets=[0, 1], controls=None, classical controls=None, control_value=None),
   Measurement(M0, target=[1], classical_store=0),
   Gate(CNOT, targets=[1], controls=[0], classical controls=None, control_value=None),
   Gate(X, targets=[0], controls=None, classical controls=[0], control_value=None),
   Gate(SWAP, targets=[0, 1], controls=None, classical controls=None, control_value=None)]

Unitaries
=========

There are a few useful functions associated with the circuit object. For example,
the :meth:`~qutip.qip.circuit.QubitCircuit.propagators` method returns a list of the unitaries associated
with the sequence of gates in the circuit. By default, the unitaries are expanded to the
full dimension of the circuit:

.. testcode::

  U_list = qc.propagators()
  print(U_list)

**Output**:

.. testoutput::

  [Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]], Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]
   [0. 0. 1. 0.]], Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[0. 0. 1. 0.]
   [0. 0. 0. 1.]
   [1. 0. 0. 0.]
   [0. 1. 0. 0.]], Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]]]

Another option is to only return the unitaries in their original dimension. This
can be achieved with the argument ``expand=False`` specified to the
:meth:`~qutip.qip.circuit.QubitCircuit.propagators`.

.. testcode::

  U_list = qc.propagators(expand=False)
  print(U_list)

**Output**:

.. testoutput::

  [Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]], Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]
   [0. 0. 1. 0.]], Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
  Qobj data =
  [[0. 1.]
   [1. 0.]], Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
  Qobj data =
  [[1. 0. 0. 0.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]
   [0. 0. 0. 1.]]]

.. _quantum_gates:

Gates
=====

The pre-defined gates for the class :class:`qutip.qip.Gate` are shown in the table below:

====================  ========================================
Gate name                           Description
====================  ========================================
"RX"                  Rotation around x axis
"RY"                  Rotation around y axis
"RZ"                  Rotation around z axis
"X"                   Pauli-X gate
"Y"                   Pauli-Y gate
"Z"                   Pauli-Z gate
"S"                   Single-qubit rotation or Z90
"T"                   Square root of S gate
"SQRTNOT"             Square root of NOT gate
"SNOT"                Hardmard gate
"PHASEGATE"           Add a phase one the state 1
"CRX"                 Controlled rotation around x axis
"CRY"                 Controlled rotation around y axis
"CRZ"                 Controlled rotation around z axis
"CX"                  Controlled X gate
"CY"                  Controlled Y gate
"CZ"                  Controlled Z gate
"CS"                  Controlled S gate
"CT"                  Controlled T gate
"CPHASE"              Controlled phase gate
"CNOT"                Controlled NOT gate
"CSIGN"               Same as CPHASE
"QASMU"               U rotation gate used as a primitive in the QASM standard
"BERKELEY"            Berkeley gate
"SWAPalpha"           SWAPalpha gate
"SWAP"                Swap the states of two qubits
"ISWAP"               Swap gate with additional phase for 01 and 10 states
"SQRTSWAP"            Square root of the SWAP gate
"SQRTISWAP"           Square root of the ISWAP gate
"FREDKIN"             Fredkin gate
"TOFFOLI"             Toffoli gate
"GLOBALPHASE"         Global phase
====================  ========================================

For some of the gates listed above, :class:`~qutip.qip.QubitCircuit` also has a primitive :func:`~qutip.qip.QubitCircuit.resolve_gates()` method that decomposes them into elementary gate sets such as CNOT or SWAP with single-qubit gates (RX, RY and RZ). However, this method is not fully optimized. It is very likely that the depth of the circuit can be further reduced by merging quantum gates. It is required that the gate resolution be carried out before the measurements to the circuit are added.

**Custom Gates**

In addition to these pre-defined gates, QuTiP also allows the user to define their own gate.
The following example shows how to define a customized gate.
The key step is to define a
gate function returning a :class:`qutip.Qobj` and save it in the attribute ``user_gates``.

.. testcode::

      from qutip.qip.circuit import Gate
      from qutip.qip.operations import rx

      def user_gate1(arg_value):
           # controlled rotation X
           mat = np.zeros((4, 4), dtype=np.complex)
           mat[0, 0] = mat[1, 1] = 1.
           mat[2:4, 2:4] = rx(arg_value)
           return Qobj(mat, dims=[[2, 2], [2, 2]])


      def user_gate2():
           # S gate
           mat = np.array([[1.,   0],
                           [0., 1.j]])
           return Qobj(mat, dims=[[2], [2]])

      qc = QubitCircuit(2)
      qc.user_gates = {"CTRLRX": user_gate1,
                       "S"     : user_gate2}

      # qubit 0 controls qubit 1
      qc.add_gate("CTRLRX", targets=[0,1], arg_value=np.pi/2)

      # qubit 1 controls qubit 0
      qc.add_gate("CTRLRX", targets=[1,0], arg_value=np.pi/2)

      # we also add a gate using a predefined Gate object
      g_T = Gate("S", targets=[1])
      qc.add_gate(g_T)
      props = qc.propagators()

      print(props[0])

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
  Qobj data =
  [[1.        +0.j         0.        +0.j         0.        +0.j
    0.        +0.j        ]
   [0.        +0.j         1.        +0.j         0.        +0.j
    0.        +0.j        ]
   [0.        +0.j         0.        +0.j         0.70710678+0.j
    0.        -0.70710678j]
   [0.        +0.j         0.        +0.j         0.        -0.70710678j
    0.70710678+0.j        ]]

.. testcode::

      print(props[1])

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE


  Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
  Qobj data =
  [[1.        +0.j         0.        +0.j         0.        +0.j
    0.        +0.j        ]
   [0.        +0.j         0.70710678+0.j         0.        +0.j
    0.        -0.70710678j]
   [0.        +0.j         0.        +0.j         1.        +0.j
    0.        +0.j        ]
   [0.        +0.j         0.        -0.70710678j 0.        +0.j
    0.70710678+0.j        ]]


.. testcode::

      print(props[2])

**Output**:

.. testoutput::
  :options: +NORMALIZE_WHITESPACE

  Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
  Qobj data =
  [[1.+0.j 0.+0.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+1.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 1.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 0.+0.j 0.+1.j]]

.. _quantum_circuit_plots:

Plotting a Quantum Circuit
===================================

A quantum circuit (described above) can directly be plotted using the QCircuit library (https://github.com/CQuIC/qcircuit).
QCiruit is a quantum circuit drawing application and is implemented directly into QuTiP.

The circuit image visualization requires LaTeX and ImageMagick for display.
The module automatically generates the LaTeX code for plotting the circuit,
produces the pdf and converts it to the png format. On Mac and Linux,
ImageMagick can be easily installed with the command conda install imagemagick if you have conda installed.
Otherwise, please follow the installation instructions on the ImageMagick documentation.

On windows, you need to download and install ImageMagick installer.
In addition, you also need perl (for ``pdfcrop``) and
Ghostscript (additional dependency of ImageMagick for png conversion).

If you want to check whether all dependencies are installed,
see if the following three commands work correctly:
``pdflatex``, ``pdfcrop`` and ``magick anypdf.pdf anypdf.png``,
where ``anypdf.pdf`` is any pdf file you have.

An example code for plotting the example quantum circuit from above is given:

.. code-block:: python

    from qutip.qip.circuit import QubitCircuit, Gate
    # create the quantum circuit
    qc = QubitCircuit(2, num_cbits=1)
    qc.add_gate("CNOT", controls=0, targets=1)
    qc.add_gate("H", targets=1)
    qc.add_gate("ISWAP", targets=[0,1])
    qc.add_measurement("M0", targets=1, classical_store=0)
    # plot the quantum circuit
    qc.png

.. image:: /figures/qip/quantum_circuit_example.png

..
   _This: is a comment, do not test the png generation as it requires additional installation!


Circuit simulation
==================

There are two different ways to simulate the action of quantum circuits using QuTiP:

- The first method utilizes unitary application through matrix products on the input states.
  This method simulates circuits exactly in a deterministic manner. This is achieved through
  :class:`~qutip.qip.CircuitSimulator`. A short guide to exact simulation can be
  found at :ref:`qip_simulator`. The teleportation notebook is also useful as an example.

- A different method of circuit simulation employs driving Hamiltonians with the ability to
  simulate circuits in the presence of noise. This can be achieved through the various classes
  in :class:`~qutip.qip.device`.A short guide to processors for QIP simulation can be found at :ref:`qip_processor`.
.. _stochastic:

*******************************************
Stochastic Solver
*******************************************

.. _stochastic-intro:

When a quantum system is subjected to continuous measurement, through homodyne detection for example, it is possible to simulate the conditional quantum state using stochastic Schrodinger and master equations. The solution of these stochastic equations are quantum trajectories, which represent the conditioned evolution of the system given a specific measurement record.

In general, the stochastic evolution of a quantum state is calculated in
QuTiP by solving the general equation

.. math::
    :label: general_form

    d \rho (t) = d_1 \rho dt + \sum_n d_{2,n} \rho dW_n,

where :math:`dW_n` is a Wiener increment, which has the expectation values :math:`E[dW] = 0` and :math:`E[dW^2] = dt`. Stochastic evolution is implemented with the :func:`qutip.stochastic.general_stochastic` function.

Stochastic Schrodinger Equation
===============================

.. _sse-solver:

The stochastic Schrodinger equation is given by (see section 4.4, [Wis09]_)

.. math::
    :label: jump_ssesolve

    d \psi(t) = - i H \psi(t) dt
                     - \sum_n \left( \frac{S_n^\dagger S_n}{2} -\frac{e_n}{2} S_n
                     + \frac{e_n^2}{8} \right) \psi(t) dt
                     + \sum_n \left( S_n - \frac{e_n}{2} \right) \psi(t) dW_n,

where :math:`H` is the Hamiltonian, :math:`S_n` are the stochastic collapse operators, and :math:`e_n` is

.. math::
   :label: jump_matrix_element

   e_n = \left<\psi(t)|S_n + S_n^\dagger|\psi(t)\right>

In QuTiP, this equation can be solved using the function :func:`qutip.stochastic.ssesolve`, which is implemented by defining :math:`d_1` and :math:`d_{2,n}` from Equation :eq:`general_form` as

.. math::
    :label: d1_def

    d_1 = -iH -  \frac{1}{2} \sum_n \left(S_n^\dagger S_n - e_n S_n + \frac{e_i^2}{4}  \right),

and

.. math::
    :label: d2_def

    d_{2, n} = S_n - \frac{e_n}{2}.

The solver :func:`qutip.stochastic.ssesolve` will construct the operators :math:`d_1` and :math:`d_{2,n}` once the user passes the Hamiltonian (``H``) and the stochastic operator list (``sc_ops``). As with the :func:`qutip.mcsolve`, the number of trajectories and the seed for the noise realisation can be fixed using the arguments: ``ntraj`` and ``noise``, respectively. If the user also requires the measurement output, the argument ``store_measurement=True`` should be included.

Additionally, homodyne and heterodyne detections can be easily simulated by passing the arguments ``method='homodyne'`` or ``method='heterodyne'`` to :func:`qutip.stochastic.ssesolve`.

Examples of how to solve the stochastic Schrodinger equation using QuTiP can be found in this `development notebook <https://nbviewer.ipython.org/github/qutip/qutip-notebooks/blob/master/development/development-ssesolve-tests.ipynb>`_.

Stochastic Master Equation
==========================

.. Stochastic Master equation

When the initial state of the system is a density matrix :math:`\rho`, the stochastic master equation solver :func:`qutip.stochastic.smesolve` must be used. The stochastic master equation is given by (see section 4.4, [Wis09]_)

.. math::
   :label: stochastic_master

    d \rho (t) = -i[H, \rho(t)] dt + D[A]\rho(t) dt + \mathcal{H}[A]\rho dW(t)

where

.. math::
    :label: dissipator

    D[A] \rho = \frac{1}{2} \left[2 A \rho A^\dagger
               - \rho A^\dagger A - A^\dagger A \rho \right],

and

.. math::
    :label: h_cal

    \mathcal{H}[A]\rho = A\rho(t) + \rho(t) A^\dagger - \tr[A\rho(t) + \rho(t) A^\dagger].


In QuTiP, solutions for the stochastic master equation are obtained using the solver :func:`qutip.stochastic.smesolve`. The implementation takes into account 2 types of collapse operators. :math:`C_i` (``c_ops``) represent the dissipation in the environment, while :math:`S_n` (``sc_ops``) are monitored operators. The deterministic part of the evolution, described by the :math:`d_1` in Equation :eq:`general_form`, takes into account all operators :math:`C_i` and :math:`S_n`:

.. math::
    :label: liouvillian

    d_1 = - i[H(t),\rho(t)]
                 + \sum_i D[C_i]\rho
                 + \sum_n D[S_n]\rho,



The stochastic part, :math:`d_{2,n}`, is given solely by the operators :math:`S_n`

.. math::
    :label: stochastic_smesolve

    d_{2,n} = S_n \rho(t) + \rho(t) S_n^\dagger - \tr \left(S_n \rho (t)
                     + \rho(t) S_n^\dagger \right)\rho(t).

As in the stochastic Schrodinger equation, the detection method can be specified using the ``method`` argument.

Example
-------

Below, we solve the dynamics for an optical cavity at 0K whose output is monitored using homodyne detection. The cavity decay rate is given by :math:`\kappa` and the :math:`\Delta` is the cavity detuning with respect to the driving field. The measurement operators can be passed using the option ``m_ops``. The homodyne current :math:`J_x` is calculated using

.. math::
    :label: measurement_result

    J_x = \langle x \rangle + dW,

where :math:`x` is the operator passed using ``m_ops``. The results are available in ``result.measurements``.

.. plot::
    :context: close-figs

    import numpy as np
    import matplotlib.pyplot as plt
    import qutip as qt

    # parameters
    DIM = 20             # Hilbert space dimension
    DELTA = 5*2*np.pi    # cavity detuning
    KAPPA = 2            # cavity decay rate
    INTENSITY = 4        # intensity of initial state
    NUMBER_OF_TRAJECTORIES = 500

    # operators
    a = qt.destroy(DIM)
    x = a + a.dag()
    H = DELTA*a.dag()* a

    rho_0 = qt.coherent(DIM, np.sqrt(INTENSITY))
    times = np.arange(0, 1, 0.0025)

    stoc_solution = qt.smesolve(H, rho_0, times,
                                c_ops=[],
                                sc_ops=[np.sqrt(KAPPA) * a],
                                e_ops=[x],
                                ntraj=NUMBER_OF_TRAJECTORIES,
                                nsubsteps=2,
                                store_measurement=True,
                                dW_factors=[1],
                                method='homodyne')

    fig, ax = plt.subplots()
    ax.set_title('Stochastic Master Equation - Homodyne Detection')
    ax.plot(times, np.array(stoc_solution.measurement).mean(axis=0)[:].real,
            'r', lw=2, label=r'$J_x$')
    ax.plot(times, stoc_solution.expect[0], 'k', lw=2,
            label=r'$\langle x \rangle$')
    ax.set_xlabel('Time')
    ax.legend()


For other examples on :func:`qutip.stochastic.smesolve`, see the `following notebook <https://nbviewer.ipython.org/github/qutip/qutip-notebooks/blob/master/development/development-smesolve-tests.ipynb>`_, as well as these notebooks available at `QuTiP Tutorials page <https://qutip.org/tutorials.html>`_: `heterodyne detection <https://nbviewer.ipython.org/github/qutip/qutip-notebooks/blob/master/examples/smesolve-heterodyne.ipynb>`_, `inneficient detection <https://nbviewer.ipython.org/github/qutip/qutip-notebooks/blob/master/examples/smesolve-inefficient-detection.ipynb>`_, and `feedback control <https://nbviewer.ipython.org/github/jrjohansson/reproduced-papers/blob/master/Reproduce-SIAM-JCO-46-445-2007-Mirrahimi.ipynb>`_.
.. _master:

*********************************
Lindblad Master Equation Solver
*********************************

.. _master-unitary:

Unitary evolution
====================
The dynamics of a closed (pure) quantum system is governed by the Schrödinger equation

.. math::
   :label: schrodinger

	i\hbar\frac{\partial}{\partial t}\Psi = \hat H \Psi,

where :math:`\Psi` is the wave function, :math:`\hat H` the Hamiltonian, and :math:`\hbar` is Planck's constant. In general, the Schrödinger equation is a partial differential equation (PDE) where both :math:`\Psi` and :math:`\hat H` are functions of space and time. For computational purposes it is useful to expand the PDE in a set of basis functions that span the Hilbert space of the Hamiltonian, and to write the equation in matrix and vector form

.. math::

   i\hbar\frac{d}{dt}\left|\psi\right> = H \left|\psi\right>

where :math:`\left|\psi\right>` is the state vector and :math:`H` is the matrix representation of the Hamiltonian. This matrix equation can, in principle, be solved by diagonalizing the Hamiltonian matrix :math:`H`. In practice, however, it is difficult to perform this diagonalization unless the size of the Hilbert space (dimension of the matrix :math:`H`) is small. Analytically, it is a formidable task to calculate the dynamics for systems with more than two states. If, in addition, we consider dissipation due to the inevitable interaction with a surrounding environment, the computational complexity grows even larger, and we have to resort to numerical calculations in all realistic situations. This illustrates the importance of numerical calculations in describing the dynamics of open quantum systems, and the need for efficient and accessible tools for this task.

The Schrödinger equation, which governs the time-evolution of closed quantum systems, is defined by its Hamiltonian and state vector. In the previous section, :ref:`tensor`, we showed how Hamiltonians and state vectors are constructed in QuTiP. Given a Hamiltonian, we can calculate the unitary (non-dissipative) time-evolution of an arbitrary state vector :math:`\left|\psi_0\right>` (``psi0``) using the QuTiP function :func:`qutip.mesolve`. It evolves the state vector and evaluates the expectation values for a set of operators ``expt_ops`` at the points in time in the list ``times``, using an ordinary differential equation solver.

For example, the time evolution of a quantum spin-1/2 system with tunneling rate 0.1 that initially is in the up state is calculated, and the  expectation values of the :math:`\sigma_z` operator evaluated, with the following code

.. plot::
    :context:

    >>> H = 2*np.pi * 0.1 * sigmax()
    >>> psi0 = basis(2, 0)
    >>> times = np.linspace(0.0, 10.0, 20)
    >>> result = sesolve(H, psi0, times, [sigmaz()])


The brackets in the fourth argument is an empty list of collapse operators, since we consider unitary evolution in this example. See the next section for examples on how dissipation is included by defining a list of collapse operators.

The function returns an instance of :class:`qutip.solver.Result`, as described in the previous section :ref:`solver_result`. The attribute ``expect`` in ``result`` is a list of expectation values for the operators that are included in the list in the fifth argument. Adding operators to this list results in a larger output list returned by the function (one array of numbers, corresponding to the times in times, for each operator)

.. plot::
    :context:

    >>> result = sesolve(H, psi0, times, [sigmaz(), sigmay()])
    >>> result.expect # doctest: +NORMALIZE_WHITESPACE
    [array([ 1.        ,  0.78914057,  0.24548559, -0.40169513, -0.8794735 ,
        -0.98636142, -0.67728219, -0.08258023,  0.54694721,  0.94581685,
         0.94581769,  0.54694945, -0.08257765, -0.67728015, -0.98636097,
        -0.87947476, -0.40169736,  0.24548326,  0.78913896,  1.        ]),
     array([ 0.00000000e+00, -6.14212640e-01, -9.69400240e-01, -9.15773457e-01,
            -4.75947849e-01,  1.64593874e-01,  7.35723339e-01,  9.96584419e-01,
             8.37167094e-01,  3.24700624e-01, -3.24698160e-01, -8.37165632e-01,
            -9.96584633e-01, -7.35725221e-01, -1.64596567e-01,  4.75945525e-01,
             9.15772479e-01,  9.69400830e-01,  6.14214701e-01,  2.77159958e-06])]

The resulting list of expectation values can easily be visualized using matplotlib's plotting functions:

.. plot::
    :context:

    >>> H = 2*np.pi * 0.1 * sigmax()
    >>> psi0 = basis(2, 0)
    >>> times = np.linspace(0.0, 10.0, 100)
    >>> result = sesolve(H, psi0, times, [sigmaz(), sigmay()])
    >>> fig, ax = plt.subplots()
    >>> ax.plot(result.times, result.expect[0]) # doctest: +SKIP
    >>> ax.plot(result.times, result.expect[1]) # doctest: +SKIP
    >>> ax.set_xlabel('Time') # doctest: +SKIP
    >>> ax.set_ylabel('Expectation values') # doctest: +SKIP
    >>> ax.legend(("Sigma-Z", "Sigma-Y")) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

If an empty list of operators is passed as fifth parameter, the :func:`qutip.mesolve` function returns a :class:`qutip.solver.Result` instance that contains a list of state vectors for the times specified in ``times``

.. plot::
    :context: close-figs

    >>> times = [0.0, 1.0]
    >>> result = mesolve(H, psi0, times, [], [])
    >>> result.states # doctest: +NORMALIZE_WHITESPACE
    [Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
     Qobj data =
     [[1.]
      [0.]], Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
     Qobj data =
     [[0.80901699+0.j        ]
      [0.        -0.58778526j]]]

.. _master-nonunitary:

Non-unitary evolution
=======================

While the evolution of the state vector in a closed quantum system is deterministic, open quantum systems are stochastic in nature. The effect of an environment on the system of interest is to induce stochastic transitions between energy levels, and to introduce uncertainty in the phase difference between states of the system. The state of an open quantum system is therefore described in terms of ensemble averaged states using the density matrix formalism. A density matrix :math:`\rho` describes a probability distribution of quantum states :math:`\left|\psi_n\right>`, in a matrix representation :math:`\rho = \sum_n p_n \left|\psi_n\right>\left<\psi_n\right|`, where :math:`p_n` is the classical probability that the system is in the quantum state :math:`\left|\psi_n\right>`. The time evolution of a density matrix :math:`\rho` is the topic of the remaining portions of this section.

.. _master-master:

The Lindblad Master equation
=============================

The standard approach for deriving the equations of motion for a system interacting with its environment is to expand the scope of the system to include the environment. The combined quantum system is then closed, and its evolution is governed by the von Neumann equation

.. math::
   :label: neumann_total

   \dot \rho_{\rm tot}(t) = -\frac{i}{\hbar}[H_{\rm tot}, \rho_{\rm tot}(t)],

the equivalent of the Schrödinger equation :eq:`schrodinger` in the density matrix formalism. Here, the total Hamiltonian

.. math::

 	H_{\rm tot} = H_{\rm sys} + H_{\rm env} + H_{\rm int},

includes the original system Hamiltonian :math:`H_{\rm sys}`, the Hamiltonian for the environment :math:`H_{\rm env}`, and a term representing the interaction between the system and its environment :math:`H_{\rm int}`. Since we are only interested in the dynamics of the system, we can at this point perform a partial trace over the environmental degrees of freedom in Eq. :eq:`neumann_total`, and thereby obtain a master equation for the motion of the original system density matrix. The most general trace-preserving and completely positive form of this evolution is the Lindblad master equation for the reduced density matrix :math:`\rho = {\rm Tr}_{\rm env}[\rho_{\rm tot}]`

.. math::
	:label: lindblad_master_equation

	\dot\rho(t)=-\frac{i}{\hbar}[H(t),\rho(t)]+\sum_n \frac{1}{2} \left[2 C_n \rho(t) C_n^\dagger - \rho(t) C_n^\dagger C_n - C_n^\dagger C_n \rho(t)\right]

where the :math:`C_n = \sqrt{\gamma_n} A_n` are collapse operators, and :math:`A_n` are the operators through which the environment couples to the system in :math:`H_{\rm int}`, and :math:`\gamma_n` are the corresponding rates.  The derivation of Eq. :eq:`lindblad_master_equation` may be found in several sources, and will not be reproduced here.  Instead, we emphasize the approximations that are required to arrive at the master equation in the form of Eq. :eq:`lindblad_master_equation` from physical arguments, and hence perform a calculation in QuTiP:

- **Separability:** At :math:`t=0` there are no correlations between the system and its environment such that the total density matrix can be written as a tensor product :math:`\rho^I_{\rm tot}(0) = \rho^I(0) \otimes \rho^I_{\rm env}(0)`.

- **Born approximation:** Requires: (1) that the state of the environment does not significantly change as a result of the interaction with the system;  (2) The system and the environment remain separable throughout the evolution. These assumptions are justified if the interaction is weak, and if the environment is much larger than the system. In summary, :math:`\rho_{\rm tot}(t) \approx \rho(t)\otimes\rho_{\rm env}`.

- **Markov approximation** The time-scale of decay for the environment :math:`\tau_{\rm env}` is much shorter than the smallest time-scale of the system dynamics :math:`\tau_{\rm sys} \gg \tau_{\rm env}`. This approximation is often deemed a "short-memory environment" as it requires that environmental correlation functions decay on a time-scale fast compared to those of the system.

- **Secular approximation** Stipulates that elements in the master equation corresponding to transition frequencies satisfy :math:`|\omega_{ab}-\omega_{cd}| \ll 1/\tau_{\rm sys}`, i.e., all fast rotating terms in the interaction picture can be neglected. It also ignores terms that lead to a small renormalization of the system energy levels. This approximation is not strictly necessary for all master-equation formalisms (e.g., the Block-Redfield master equation), but it is required for arriving at the Lindblad form :eq:`lindblad_master_equation` which is used in :func:`qutip.mesolve`.


For systems with environments satisfying the conditions outlined above, the Lindblad master equation :eq:`lindblad_master_equation` governs the time-evolution of the system density matrix, giving an ensemble average of the system dynamics. In order to ensure that these approximations are not violated, it is important that the decay rates :math:`\gamma_n` be smaller than the minimum energy splitting in the system Hamiltonian. Situations that demand special attention therefore include, for example, systems strongly coupled to their environment, and systems with degenerate or nearly degenerate energy levels.


For non-unitary evolution of a quantum systems, i.e., evolution that includes
incoherent processes such as relaxation and dephasing, it is common to use
master equations. In QuTiP, the same function (:func:`qutip.mesolve`) is used for
evolution both according to the Schrödinger equation and to the master equation,
even though these two equations of motion are very different. The :func:`qutip.mesolve`
function automatically determines if it is sufficient to use the Schrödinger
equation (if no collapse operators were given) or if it has to use the
master equation (if collapse operators were given). Note that to calculate
the time evolution according to the Schrödinger equation is easier and much
faster (for large systems) than using the master equation, so if possible the
solver will fall back on using the Schrödinger equation.

What is new in the master equation compared to the Schrödinger equation are
processes that describe dissipation in the quantum system due to its interaction
with an environment. These environmental interactions are defined by the
operators through which the system couples to the environment, and rates that
describe the strength of the processes.

In QuTiP, the product of the square root of the rate and the operator that
describe the dissipation process is called a collapse operator. A list of
collapse operators (``c_ops``) is passed as the fourth argument to the
:func:`qutip.mesolve` function in order to define the dissipation processes in the master
equation. When the ``c_ops`` isn't empty, the :func:`qutip.mesolve` function will use
the master equation instead of the unitary Schrödinger equation.

Using the example with the spin dynamics from the previous section, we can
easily add a relaxation process (describing the dissipation of energy from the
spin to its environment), by adding ``np.sqrt(0.05) * sigmax()`` to
the previously empty list in the fourth parameter to the :func:`qutip.mesolve` function:


.. plot::
    :context:

    >>> times = np.linspace(0.0, 10.0, 100)
    >>> result = mesolve(H, psi0, times, [np.sqrt(0.05) * sigmax()], [sigmaz(), sigmay()])
    >>> fig, ax = plt.subplots()
    >>> ax.plot(times, result.expect[0]) # doctest: +SKIP
    >>> ax.plot(times, result.expect[1]) # doctest: +SKIP
    >>> ax.set_xlabel('Time') # doctest: +SKIP
    >>> ax.set_ylabel('Expectation values') # doctest: +SKIP
    >>> ax.legend(("Sigma-Z", "Sigma-Y"))  # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP


Here, 0.05 is the rate and the operator :math:`\sigma_x` (:func:`qutip.operators.sigmax`) describes the dissipation
process.

Now a slightly more complex example: Consider a two-level atom coupled to a leaky single-mode cavity through a dipole-type interaction, which supports a coherent exchange of quanta between the two systems. If the atom initially is in its groundstate and the cavity in a 5-photon Fock state, the dynamics is calculated with the lines following code

.. plot::
    :context: close-figs

    >>> times = np.linspace(0.0, 10.0, 200)
    >>> psi0 = tensor(fock(2,0), fock(10, 5))
    >>> a  = tensor(qeye(2), destroy(10))
    >>> sm = tensor(destroy(2), qeye(10))
    >>> H = 2 * np.pi * a.dag() * a + 2 * np.pi * sm.dag() * sm + 2 * np.pi * 0.25 * (sm * a.dag() + sm.dag() * a)
    >>> result = mesolve(H, psi0, times, [np.sqrt(0.1)*a], [a.dag()*a, sm.dag()*sm])
    >>> plt.figure() # doctest: +SKIP
    >>> plt.plot(times, result.expect[0]) # doctest: +SKIP
    >>> plt.plot(times, result.expect[1]) # doctest: +SKIP
    >>> plt.xlabel('Time') # doctest: +SKIP
    >>> plt.ylabel('Expectation values') # doctest: +SKIP
    >>> plt.legend(("cavity photon number", "atom excitation probability")) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP
.. _floquet:

*****************
Floquet Formalism
*****************

.. _floquet-intro:

Introduction
============

Many time-dependent problems of interest are periodic. The dynamics of such systems can be solved for directly by numerical integration of the Schrödinger or Master equation, using the time-dependent Hamiltonian. But they can also be transformed into time-independent problems using the Floquet formalism. Time-independent problems can be solve much more efficiently, so such a transformation is often very desirable.

In the standard derivations of the Lindblad and Bloch-Redfield master equations the Hamiltonian describing the system under consideration is assumed to be time independent. Thus, strictly speaking, the standard forms of these master equation formalisms should not blindly be applied to system with time-dependent Hamiltonians. However, in many relevant cases, in particular for weak driving, the standard master equations still turns out to be useful for many time-dependent problems. But a more rigorous approach would be to rederive the master equation taking the time-dependent nature of the Hamiltonian into account from the start. The Floquet-Markov Master equation is one such a formalism, with important applications for strongly driven systems (see e.g., [Gri98]_).

Here we give an overview of how the Floquet and Floquet-Markov formalisms can be used for solving time-dependent problems in QuTiP. To introduce the terminology and naming conventions used in QuTiP we first give a brief summary of quantum Floquet theory.

.. _floquet-unitary:

Floquet theory for unitary evolution
====================================

The Schrödinger equation with a time-dependent Hamiltonian :math:`H(t)` is

.. math::
   :label: eq_td_schrodinger

	H(t)\Psi(t) = i\hbar\frac{\partial}{\partial t}\Psi(t),

where :math:`\Psi(t)` is the wave function solution. Here we are interested in problems with periodic time-dependence, i.e., the Hamiltonian satisfies :math:`H(t) = H(t+T)` where :math:`T` is the period. According to the Floquet theorem, there exist solutions to :eq:`eq_td_schrodinger` on the form

.. math::
   :label: eq_floquet_states

    \Psi_\alpha(t) = \exp(-i\epsilon_\alpha t/\hbar)\Phi_\alpha(t),

where :math:`\Psi_\alpha(t)` are the *Floquet states* (i.e., the set of wave function solutions to the Schrödinger equation), :math:`\Phi_\alpha(t)=\Phi_\alpha(t+T)` are the periodic *Floquet modes*, and :math:`\epsilon_\alpha` are the *quasienergy levels*. The quasienergy levels are constants in time, but only uniquely defined up to multiples of :math:`2\pi/T` (i.e., unique value in the interval :math:`[0, 2\pi/T]`).

If we know the Floquet modes (for :math:`t \in [0,T]`) and the quasienergies for a particular :math:`H(t)`, we can easily decompose any initial wavefunction :math:`\Psi(t=0)` in the Floquet states and immediately obtain the solution for arbitrary :math:`t`

.. math::
   :label: eq_floquet_wavefunction_expansion

    \Psi(t) = \sum_\alpha c_\alpha \Psi_\alpha(t) = \sum_\alpha c_\alpha \exp(-i\epsilon_\alpha t/\hbar)\Phi_\alpha(t),

where the coefficients :math:`c_\alpha` are determined by the initial wavefunction :math:`\Psi(0) = \sum_\alpha c_\alpha \Psi_\alpha(0)`.

This formalism is useful for finding :math:`\Psi(t)` for a given :math:`H(t)` only if we can obtain the Floquet modes :math:`\Phi_a(t)` and quasienergies :math:`\epsilon_\alpha` more easily than directly solving :eq:`eq_td_schrodinger`. By substituting :eq:`eq_floquet_states` into the Schrödinger equation :eq:`eq_td_schrodinger` we obtain an eigenvalue equation for the Floquet modes and quasienergies

.. math::
   :label: eq_floquet_eigen_problem

    \mathcal{H}(t)\Phi_\alpha(t) = \epsilon_\alpha\Phi_\alpha(t),

where :math:`\mathcal{H}(t) = H(t) - i\hbar\partial_t`. This eigenvalue problem could be solved analytically or numerically, but in QuTiP we use an alternative approach for numerically finding the Floquet states and quasienergies [see e.g. Creffield et al., Phys. Rev. B 67, 165301 (2003)]. Consider the propagator for the time-dependent Schrödinger equation :eq:`eq_td_schrodinger`, which by definition satisfies

.. math::

    U(T+t,t)\Psi(t) = \Psi(T+t).

Inserting the Floquet states from :eq:`eq_floquet_states` into this expression results in

.. math::
    U(T+t,t)\exp(-i\epsilon_\alpha t/\hbar)\Phi_\alpha(t) = \exp(-i\epsilon_\alpha(T+t)/\hbar)\Phi_\alpha(T+t),

or, since :math:`\Phi_\alpha(T+t)=\Phi_\alpha(t)`,

.. math::
    U(T+t,t)\Phi_\alpha(t) = \exp(-i\epsilon_\alpha T/\hbar)\Phi_\alpha(t) = \eta_\alpha \Phi_\alpha(t),

which shows that the Floquet modes are eigenstates of the one-period propagator. We can therefore find the Floquet modes and quasienergies :math:`\epsilon_\alpha = -\hbar\arg(\eta_\alpha)/T` by numerically calculating :math:`U(T+t,t)` and diagonalizing it. In particular this method is useful to find :math:`\Phi_\alpha(0)` by calculating and diagonalize :math:`U(T,0)`.

The Floquet modes at arbitrary time :math:`t` can then be found by propagating :math:`\Phi_\alpha(0)` to :math:`\Phi_\alpha(t)` using the wave function propagator :math:`U(t,0)\Psi_\alpha(0) = \Psi_\alpha(t)`, which for the Floquet modes yields

.. math::

    U(t,0)\Phi_\alpha(0) = \exp(-i\epsilon_\alpha t/\hbar)\Phi_\alpha(t),

so that :math:`\Phi_\alpha(t) = \exp(i\epsilon_\alpha t/\hbar) U(t,0)\Phi_\alpha(0)`. Since :math:`\Phi_\alpha(t)` is periodic we only need to evaluate it for :math:`t \in [0, T]`, and from :math:`\Phi_\alpha(t \in [0,T])` we can directly evaluate :math:`\Phi_\alpha(t)`, :math:`\Psi_\alpha(t)` and :math:`\Psi(t)` for arbitrary large :math:`t`.

Floquet formalism in QuTiP
--------------------------

QuTiP provides a family of functions to calculate the Floquet modes and quasi energies, Floquet state decomposition, etc., given a time-dependent Hamiltonian on the *callback format*, *list-string format* and *list-callback format* (see, e.g., :func:`qutip.mesolve` for details).

Consider for example the case of a strongly driven two-level atom, described by the Hamiltonian

.. math::
   :label: eq_driven_qubit

    H(t) = -\frac{1}{2}\Delta\sigma_x - \frac{1}{2}\epsilon_0\sigma_z + \frac{1}{2}A\sin(\omega t)\sigma_z.

In QuTiP we can define this Hamiltonian as follows:

.. plot::
   :context:

   >>> delta = 0.2 * 2*np.pi
   >>> eps0 = 1.0 * 2*np.pi
   >>> A = 2.5 * 2*np.pi
   >>> omega = 1.0 * 2*np.pi
   >>> H0 = - delta/2.0 * sigmax() - eps0/2.0 * sigmaz()
   >>> H1 = A/2.0 * sigmaz()
   >>> args = {'w': omega}
   >>> H = [H0, [H1, 'sin(w * t)']]

The :math:`t=0` Floquet modes corresponding to the Hamiltonian :eq:`eq_driven_qubit` can then be calculated using the :func:`qutip.floquet.floquet_modes` function, which returns lists containing the Floquet modes and the quasienergies

.. plot::
   :context:

   >>> T = 2*np.pi / omega
   >>> f_modes_0, f_energies = floquet_modes(H, T, args)
   >>> f_energies # doctest: +NORMALIZE_WHITESPACE
   array([-2.83131212,  2.83131212])
   >>> f_modes_0 # doctest: +NORMALIZE_WHITESPACE
   [Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[ 0.72964231+0.j      ]
    [-0.39993746+0.554682j]],
   Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[0.39993746+0.554682j]
    [0.72964231+0.j      ]]]

For some problems interesting observations can be draw from the quasienergy levels alone. Consider for example the quasienergies for the driven two-level system introduced above as a function of the driving amplitude, calculated and plotted in the following example. For certain driving amplitudes the quasienergy levels cross. Since the quasienergies can be associated with the time-scale of the long-term dynamics due that the driving, degenerate quasienergies indicates a "freezing" of the dynamics (sometimes known as coherent destruction of tunneling).

.. plot::
   :context:

   >>> delta = 0.2 * 2*np.pi
   >>> eps0  = 0.0 * 2*np.pi
   >>> omega = 1.0 * 2*np.pi
   >>> A_vec = np.linspace(0, 10, 100) * omega
   >>> T = (2*np.pi)/omega
   >>> tlist = np.linspace(0.0, 10 * T, 101)
   >>> spsi0 = basis(2,0)
   >>> q_energies = np.zeros((len(A_vec), 2))
   >>> H0 = delta/2.0 * sigmaz() - eps0/2.0 * sigmax()
   >>> args = {'w': omega}
   >>> for idx, A in enumerate(A_vec): # doctest: +SKIP
   >>>   H1 = A/2.0 * sigmax() # doctest: +SKIP
   >>>   H = [H0, [H1, lambda t, args: np.sin(args['w']*t)]] # doctest: +SKIP
   >>>   f_modes, f_energies = floquet_modes(H, T, args, True) # doctest: +SKIP
   >>>   q_energies[idx,:] = f_energies # doctest: +SKIP
   >>> plt.figure() # doctest: +SKIP
   >>> plt.plot(A_vec/omega, q_energies[:,0] / delta, 'b', A_vec/omega, q_energies[:,1] / delta, 'r') # doctest: +SKIP
   >>> plt.xlabel(r'$A/\omega$') # doctest: +SKIP
   >>> plt.ylabel(r'Quasienergy / $\Delta$') # doctest: +SKIP
   >>> plt.title(r'Floquet quasienergies') # doctest: +SKIP
   >>> plt.show() # doctest: +SKIP

Given the Floquet modes at :math:`t=0`, we obtain the Floquet mode at some later time :math:`t` using the function :func:`qutip.floquet.floquet_mode_t`:

.. plot::
   :context: close-figs

   >>> f_modes_t = floquet_modes_t(f_modes_0, f_energies, 2.5, H, T, args)
   >>> f_modes_t # doctest: +SKIP
   [Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[-0.89630512-0.23191946j]
    [ 0.37793106-0.00431336j]],
   Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
   Qobj data =
   [[-0.37793106-0.00431336j]
    [-0.89630512+0.23191946j]]]

The purpose of calculating the Floquet modes is to find the wavefunction solution to the original problem :eq:`eq_driven_qubit` given some initial state :math:`\left|\psi_0\right>`. To do that, we first need to decompose the initial state in the Floquet states, using the function :func:`qutip.floquet.floquet_state_decomposition`

.. plot::
   :context:

   >>> psi0 = rand_ket(2)
   >>> f_coeff = floquet_state_decomposition(f_modes_0, f_energies, psi0)
   >>> f_coeff # doctest: +SKIP
   [(-0.645265993068382+0.7304552549315746j),
   (0.15517002114250228-0.1612116102238258j)]

and given this decomposition of the initial state in the Floquet states we can easily evaluate the wavefunction that is the solution to :eq:`eq_driven_qubit` at an arbitrary time :math:`t` using the function :func:`qutip.floquet.floquet_wavefunction_t`

.. plot::
   :context:

   >>> t = 10 * np.random.rand()
   >>> psi_t = floquet_wavefunction_t(f_modes_0, f_energies, f_coeff, t, H, T, args)

The following example illustrates how to use the functions introduced above to calculate and plot the time-evolution of :eq:`eq_driven_qubit`.

.. plot:: guide/scripts/floquet_ex1.py
   :width: 4.0in
   :include-source:

Pre-computing the Floquet modes for one period
----------------------------------------------

When evaluating the Floquet states or the wavefunction at many points in time it is useful to pre-compute the Floquet modes for the first period of the driving with the required resolution. In QuTiP the function :func:`qutip.floquet.floquet_modes_table` calculates a table of Floquet modes which later can be used together with the function :func:`qutip.floquet.floquet_modes_t_lookup` to efficiently lookup the Floquet mode at an arbitrary time. The following example illustrates how the example from the previous section can be solved more efficiently using these functions for pre-computing the Floquet modes.

.. plot:: guide/scripts/floquet_ex2.py
   :width: 4.0in
   :include-source:

Note that the parameters and the Hamiltonian used in this example is not the same as in the previous section, and hence the different appearance of the resulting figure.

For convenience, all the steps described above for calculating the evolution of a quantum system using the Floquet formalisms are encapsulated in the function :func:`qutip.floquet.fsesolve`. Using this function, we could have achieved the same results as in the examples above using

.. code-block:: python

    output = fsesolve(H, psi0=psi0, tlist=tlist, e_ops=[qutip.num(2)], args=args)
    p_ex = output.expect[0]

.. _floquet-dissipative:

Floquet theory for dissipative evolution
========================================

A driven system that is interacting with its environment is not necessarily well described by the standard Lindblad master equation, since its dissipation process could be time-dependent due to the driving. In such cases a rigorious approach would be to take the driving into account when deriving the master equation. This can be done in many different ways, but one way common approach is to derive the master equation in the Floquet basis. That approach results in the so-called Floquet-Markov master equation, see Grifoni et al., Physics Reports 304, 299 (1998) for details.


The Floquet-Markov master equation in QuTiP
-------------------------------------------

The QuTiP function :func:`qutip.floquet.fmmesolve` implements the Floquet-Markov master equation. It calculates the dynamics of a system given its initial state, a time-dependent Hamiltonian, a list of operators through which the system couples to its environment and a list of corresponding spectral-density functions that describes the environment. In contrast to the :func:`qutip.mesolve` and :func:`qutip.mcsolve`, and the :func:`qutip.floquet.fmmesolve` does characterize the environment with dissipation rates, but extract the strength of the coupling to the environment from the noise spectral-density functions and the instantaneous Hamiltonian parameters (similar to the Bloch-Redfield master equation solver :func:`qutip.bloch_redfield.brmesolve`).

.. note::

    Currently the :func:`qutip.floquet.fmmesolve` can only accept a single environment coupling operator and spectral-density function.

The noise spectral-density function of the environment is implemented as a Python callback function that is passed to the solver. For example:


.. code-block:: python

    gamma1 = 0.1
    def noise_spectrum(omega):
        return 0.5 * gamma1 * omega/(2*pi)

The other parameters are similar to the :func:`qutip.mesolve` and :func:`qutip.mcsolve`, and the same format for the return value is used :class:`qutip.solver.Result`. The following example extends the example studied above, and uses :func:`qutip.floquet.fmmesolve` to introduce dissipation into the calculation

.. plot:: guide/scripts/floquet_ex3.py
   :width: 4.0in
   :include-source:

Alternatively, we can let the :func:`qutip.floquet.fmmesolve` function transform the density matrix at each time step back to the computational basis, and calculating the expectation values for us, but using::

    output = fmmesolve(H, psi0, tlist, [sigmax()], [num(2)], [noise_spectrum], T, args)
    p_ex = output.expect[0]
.. _stochastic_photo:

********************************
Stochastic Solver - Photocurrent
********************************

.. _photocurrent-intro:

Photocurrent method, like monte-carlo method, allows for simulating an
individual realization of the system evolution under continuous measurement.

Closed system
-------------

.. photocurent_Schrodinger_equation

Photocurrent evolution have the state evolve deterministically between quantum jumps.
During the deterministic part, the system evolve by schrodinger equation with a
non-hermitian, norm conserving effective Hamiltonian.

.. math::
	:label: pssesolve_heff

	H_{\rm eff}=H_{\rm sys}+
	\frac{i\hbar}{2}\left( -\sum_{n}C^{+}_{n}C_{n}+ |C_{n} \psi |^2\right).

With :math:`C_{n}`, the collapse operators.
This effective Hamiltonian is equivalent to the monte-carlo effective
Hamiltonian with an extra term to keep the state normalized.
At each time step of :math:`\delta t`, the wave function has a probability

.. math::
	:label: pssesolve_jump_prob

	\delta p_{n} = \left<\psi(t)|C_{n}^{+}C_{n}|\psi(t)\right>  \delta t

of making a quantum jump. :math:`\delta t` must be chosen small enough to keep
that probability small :math:`\delta p << 1`. *If multiple jumps happen at the
same time step, the state become unphysical.*
Each jump result in a sharp variation of the state by,

.. math::
	:label: pssesolve_jump

	\delta \psi = \left( \frac{C_n \psi} {\left| C_n \psi  \right|} - \psi \right)

The basic photocurrent method directly integrates these equations to the first-order.
Starting from a state :math:`\left|\psi(0)\right>`, it evolves the state according to

.. math::
	:label: pssesolve_sde

	\delta \psi(t) = - i H_{\rm sys} \psi(t) \delta t + \sum_n \left(
					 -\frac{C_n^{+} C_n}{2}  \psi(t) \delta t
					 + \frac{ \left| C_n \psi  \right| ^2}{2} \delta t
	                 +  \delta N_n \left( \frac{C_n \psi}
					 {\left| C_n \psi  \right|} - \psi \right)\right),

for each time-step.
Here :math:`\delta N = 1` with a probability of :math:`\delta \omega` and
:math:`\delta N_n = 0` with a probability of :math:`1-\delta \omega`.

Trajectories obtained with this algorithm are equivalent to those obtained with
monte-carlo evolution (up to :math:`O(\delta t^2)`).
In most cases, :func:`qutip.mcsolve` is more efficient than
:func:`qutip.photocurrent_sesolve`.

Open system
-----------
.. photocurent_Master_equation

Photocurrent approach allows to obtain trajectories for a system with
both measured and dissipative interaction with the bath.
The system evolves according to the master equation between jumps with a modified
liouvillian

.. math::
	:label: master_equation

	L_{\rm eff}(\rho(t)) = L_{\rm sys}(\rho(t)) +
	                      \sum_{n}\left(
						  \rm{tr} \left(C_{n}^{+}C_{n}  \rho C_{n}^{+}C_{n} \right)
						  - C_{n}^{+}C_{n}  \rho C_{n}^{+}C_{n} \right),

with the probability of jumps in a time step :math:`\delta t` given by

.. math::
	:label: psmesolve_rate

	\delta p = \rm{tr} \left( C \rho C^{+} \right) \delta t.

After a jump, the density matrix become

.. math::

	\rho' = \frac{C \rho C^{+}}{\rm{tr} \left( C \rho C^{+} \right)}.

The evolution of the system at each time step if thus given by

.. math::
	:label: psmesolve_sde

	\rho(t + \delta t) = \rho(t) + L_{\rm eff}(\rho) \delta t + \delta N
	\left(\frac{C \rho C^{+}}{\rm{tr} \left( C \rho C^{+} \right)} - \rho \right).
.. _time:

*************************************************
Solving Problems with Time-dependent Hamiltonians
*************************************************


Methods for Writing Time-Dependent Operators
============================================

In the previous examples of quantum evolution,
we assumed that the systems under consideration were described by time-independent Hamiltonians.
However, many systems have explicit time dependence in either the Hamiltonian,
or the collapse operators describing coupling to the environment, and sometimes both components might depend on time.
The time-evolutions  solvers
:func:`qutip.mesolve`, :func:`qutip.mcsolve`, :func:`qutip.sesolve`, :func:`qutip.brmesolve`
:func:`qutip.ssesolve`, :func:`qutip.photocurrent_sesolve`, :func:`qutip.smesolve`, and :func:`qutip.photocurrent_mesolve`
are all capable of handling time-dependent Hamiltonians and collapse terms.
There are, in general, three different ways to implement time-dependent problems in QuTiP:


1. **Function based**: Hamiltonian / collapse operators expressed using [qobj, func] pairs, where the time-dependent coefficients of the Hamiltonian (or collapse operators) are expressed using Python functions.

2. **String (Cython) based**: The Hamiltonian and/or collapse operators are expressed as a list of [qobj, string] pairs, where the time-dependent coefficients are represented as strings.  The resulting Hamiltonian is then compiled into C code using Cython and executed.

3. **Array Based**: The Hamiltonian and/or collapse operators are expressed as a list of [qobj, np.array] pairs. The arrays are 1 dimensional and dtype are complex or float. They must contain one value for each time in the tlist given to the solver. Cubic spline interpolation will be used between the given times.

4. **Hamiltonian function (outdated)**: The Hamiltonian is itself a Python function with time-dependence.  Collapse operators must be time independent using this input format.


Give the multiple choices of input style, the first question that arrises is which option to choose?
In short, the function based method (option #1) is the most general,
allowing for essentially arbitrary coefficients expressed via user defined functions.
However, by automatically compiling your system into C++ code,
the second option (string based) tends to be more efficient and will run faster
[This is also the only format that is supported in the :func:`qutip.brmesolve` solver].
Of course, for small system sizes and evolution times, the difference will be minor.
Although this method does not support all time-dependent coefficients that one can think of,
it does support essentially all problems that one would typically encounter.
Time-dependent coefficients using any of the following functions,
or combinations thereof (including constants) can be compiled directly into C++-code::

  'abs', 'acos', 'acosh', 'arg', 'asin', 'asinh', 'atan', 'atanh', 'conj',
   'cos', 'cosh','exp', 'erf', 'zerf', 'imag', 'log', 'log10', 'norm', 'pi',
   'proj', 'real', 'sin', 'sinh', 'sqrt', 'tan', 'tanh'

In addition, QuTiP supports cubic spline based interpolation functions [:ref:`time-interp`].

If you require mathematical functions other than those listed above,
it is possible to call any of the functions in the NumPy library using the prefix ``np.``
before the function name in the string, i.e ``'np.sin(t)'`` and  ``scipy.special`` imported as ``spe``.
This includes a wide range of functionality, but comes with a small overhead created by going from C++->Python->C++.

Finally option #4, expressing the Hamiltonian as a Python function,
is the original method for time dependence in QuTiP 1.x.
However, this method is somewhat less efficient then the previously mentioned methods.
However, in contrast to the other options
this method can be used in implementing time-dependent Hamiltonians that cannot be
expressed as a function of constant operators with time-dependent coefficients.

A collection of examples demonstrating the simulation of time-dependent problems can be found on the `tutorials <https://qutip.org/tutorials.html>`_ web page.

.. _time-function:

Function Based Time Dependence
==============================

A very general way to write a time-dependent Hamiltonian or collapse operator is by using Python functions as the time-dependent coefficients.  To accomplish this, we need to write a Python function that returns the time-dependent coefficient.  Additionally, we need to tell QuTiP that a given Hamiltonian or collapse operator should be associated with a given Python function.  To do this, one needs to specify operator-function pairs in list format: ``[Op, py_coeff]``, where ``Op`` is a given Hamiltonian or collapse operator and ``py_coeff`` is the name of the Python function representing the coefficient.  With this format, the form of the Hamiltonian for both ``mesolve`` and ``mcsolve`` is:

>>> H = [H0, [H1, py_coeff1], [H2, py_coeff2], ...] # doctest: +SKIP

where ``H0`` is a time-independent Hamiltonian, while ``H1``,``H2``, are time dependent. The same format can be used for collapse operators:

>>> c_ops = [[C0, py_coeff0], C1, [C2, py_coeff2], ...] # doctest: +SKIP

Here we have demonstrated that the ordering of time-dependent and time-independent terms does not matter.  In addition, any or all of the collapse operators may be time dependent.

.. note:: While, in general, you can arrange time-dependent and time-independent terms in any order you like, it is best to place all time-independent terms first.

As an example, we will look at an example that has a time-dependent Hamiltonian of the form :math:`H=H_{0}-f(t)H_{1}` where :math:`f(t)` is the time-dependent driving strength given as :math:`f(t)=A\exp\left[-\left( t/\sigma \right)^{2}\right]`.  The follow code sets up the problem

.. plot::
    :context:

    ustate = basis(3, 0)
    excited = basis(3, 1)
    ground = basis(3, 2)

    N = 2 # Set where to truncate Fock state for cavity
    sigma_ge = tensor(qeye(N), ground * excited.dag())  # |g><e|
    sigma_ue = tensor(qeye(N), ustate * excited.dag())  # |u><e|
    a = tensor(destroy(N), qeye(3))
    ada = tensor(num(N), qeye(3))

    c_ops = []  # Build collapse operators
    kappa = 1.5 # Cavity decay rate
    c_ops.append(np.sqrt(kappa) * a)
    gamma = 6  # Atomic decay rate
    c_ops.append(np.sqrt(5*gamma/9) * sigma_ue) # Use Rb branching ratio of 5/9 e->u
    c_ops.append(np.sqrt(4*gamma/9) * sigma_ge) # 4/9 e->g

    t = np.linspace(-15, 15, 100) # Define time vector
    psi0 = tensor(basis(N, 0), ustate) # Define initial state

    state_GG = tensor(basis(N, 1), ground) # Define states onto which to project
    sigma_GG = state_GG * state_GG.dag()
    state_UU = tensor(basis(N, 0), ustate)
    sigma_UU = state_UU * state_UU.dag()

    g = 5  # coupling strength
    H0 = -g * (sigma_ge.dag() * a + a.dag() * sigma_ge)  # time-independent term
    H1 = (sigma_ue.dag() + sigma_ue)  # time-dependent term

Given that we have a single time-dependent Hamiltonian term, and constant collapse terms, we need to specify a single Python function for the coefficient :math:`f(t)`.  In this case, one can simply do

.. plot::
    :context:

    def H1_coeff(t, args):
        return 9 * np.exp(-(t / 5.) ** 2)

In this case, the return value dependents only on time.  However, when specifying Python functions for coefficients, **the function must have (t,args) as the input variables, in that order**.  Having specified our coefficient function, we can now specify the Hamiltonian in list format and call the solver (in this case :func:`qutip.mesolve`)

.. plot::
    :context:

    H = [H0,[H1, H1_coeff]]
    output = mesolve(H, psi0, t, c_ops, [ada, sigma_UU, sigma_GG])

We can call the Monte Carlo solver in the exact same way (if using the default ``ntraj=500``):


..
  Hacky fix because plot has complicated conditional code execution

.. doctest::
    :skipif: True

    output = mcsolve(H, psi0, t, c_ops, [ada, sigma_UU, sigma_GG])

The output from the master equation solver is identical to that shown in the examples, the Monte Carlo however will be noticeably off, suggesting we should increase the number of trajectories for this example.  In addition, we can also consider the decay of a simple Harmonic oscillator with time-varying decay rate

.. plot::
    :context:

    kappa = 0.5

    def col_coeff(t, args):  # coefficient function
        return np.sqrt(kappa * np.exp(-t))

    N = 10  # number of basis states
    a = destroy(N)
    H = a.dag() * a  # simple HO
    psi0 = basis(N, 9)  # initial state
    c_ops = [[a, col_coeff]]  # time-dependent collapse term
    times = np.linspace(0, 10, 100)
    output = mesolve(H, psi0, times, c_ops, [a.dag() * a])


Using the args variable
------------------------
In the previous example we hardcoded all of the variables, driving amplitude :math:`A` and width :math:`\sigma`, with their numerical values.  This is fine for problems that are specialized, or that we only want to run once.  However, in many cases, we would like to change the parameters of the problem in only one location (usually at the top of the script), and not have to worry about manually changing the values on each run.  QuTiP allows you to accomplish this using the keyword ``args`` as an input to the solvers.  For instance, instead of explicitly writing 9 for the amplitude and 5 for the width of the gaussian driving term, we can make us of the args variable

.. plot::
    :context:

    def H1_coeff(t, args):
        return args['A'] * np.exp(-(t/args['sigma'])**2)

or equivalently,

.. plot::
    :context:

    def H1_coeff(t, args):
          A = args['A']
          sig = args['sigma']
          return A * np.exp(-(t / sig) ** 2)

where args is a Python dictionary of ``key: value`` pairs ``args = {'A': a, 'sigma': b}`` where ``a`` and ``b`` are the two parameters for the amplitude and width, respectively.  Of course, we can always hardcode the values in the dictionary as well ``args = {'A': 9, 'sigma': 5}``, but there is much more flexibility by using variables in ``args``.  To let the solvers know that we have a set of args to pass we append the ``args`` to the end of the solver input:

.. plot::
    :context:

    output = mesolve(H, psi0, times, c_ops, [a.dag() * a], args={'A': 9, 'sigma': 5})

or to keep things looking pretty

.. plot::
    :context:

    args = {'A': 9, 'sigma': 5}
    output = mesolve(H, psi0, times, c_ops, [a.dag() * a], args=args)

Once again, the Monte Carlo solver :func:`qutip.mcsolve` works in an identical manner.

.. _time-string:

String Format Method
=====================

.. note:: You must have Cython installed on your computer to use this format.  See :ref:`install` for instructions on installing Cython.

The string-based time-dependent format works in a similar manner as the previously discussed Python function method.  That being said, the underlying code does something completely different.  When using this format, the strings used to represent the time-dependent coefficients, as well as Hamiltonian and collapse operators, are rewritten as Cython code using a code generator class and then compiled into C code.  The details of this meta-programming will be published in due course.  however, in short, this can lead to a substantial reduction in time for complex time-dependent problems, or when simulating over long intervals.

Like the previous method, the string-based format uses a list pair format ``[Op, str]`` where ``str`` is now a string representing the time-dependent coefficient.  For our first example, this string would be ``'9 * exp(-(t / 5.) ** 2)'``.  The Hamiltonian in this format would take the form:

.. plot::
   :context:

   ustate = basis(3, 0)
   excited = basis(3, 1)
   ground = basis(3, 2)

   N = 2 # Set where to truncate Fock state for cavity

   sigma_ge = tensor(qeye(N), ground * excited.dag())  # |g><e|
   sigma_ue = tensor(qeye(N), ustate * excited.dag())  # |u><e|
   a = tensor(destroy(N), qeye(3))
   ada = tensor(num(N), qeye(3))

   c_ops = []  # Build collapse operators
   kappa = 1.5 # Cavity decay rate
   c_ops.append(np.sqrt(kappa) * a)
   gamma = 6  # Atomic decay rate
   c_ops.append(np.sqrt(5*gamma/9) * sigma_ue) # Use Rb branching ratio of 5/9 e->u
   c_ops.append(np.sqrt(4*gamma/9) * sigma_ge) # 4/9 e->g

   t = np.linspace(-15, 15, 100) # Define time vector
   psi0 = tensor(basis(N, 0), ustate) # Define initial state
   state_GG = tensor(basis(N, 1), ground) # Define states onto which to project
   sigma_GG = state_GG * state_GG.dag()
   state_UU = tensor(basis(N, 0), ustate)
   sigma_UU = state_UU * state_UU.dag()

   g = 5  # coupling strength
   H0 = -g * (sigma_ge.dag() * a + a.dag() * sigma_ge)  # time-independent term
   H1 = (sigma_ue.dag() + sigma_ue)  # time-dependent term


.. plot::
    :context:

    H = [H0, [H1, '9 * exp(-(t / 5) ** 2)']]

Notice that this is a valid Hamiltonian for the string-based format as ``exp`` is included in the above list of suitable functions. Calling the solvers is the same as before:

.. plot::
   :context:

   output = mesolve(H, psi0, t, c_ops, [a.dag() * a])

We can also use the ``args`` variable in the same manner as before, however we must rewrite our string term to read: ``'A * exp(-(t / sig) ** 2)'``

.. plot::
    :context:

    H = [H0, [H1, 'A * exp(-(t / sig) ** 2)']]
    args = {'A': 9, 'sig': 5}
    output = mesolve(H, psi0, times, c_ops, [a.dag()*a], args=args)


.. important:: Naming your ``args`` variables ``exp``, ``sin``, ``pi`` etc. will cause errors when using the string-based format.

Collapse operators are handled in the exact same way.


.. _time-interp:

Modeling Non-Analytic and/or Experimental Time-Dependent Parameters using Interpolating Functions
=================================================================================================

Sometimes it is necessary to model a system where the time-dependent parameters are non-analytic functions, or are derived from experimental data (i.e. a collection of data points).  In these situations, one can use interpolating functions as an approximate functional form for input into a time-dependent solver.  QuTiP includes it own custom cubic spline interpolation class :class:`qutip.interpolate.Cubic_Spline` to provide this functionality.  To see how this works, lets first generate some noisy data:

.. plot::
    :context:

    t = np.linspace(-15, 15, 100)
    func = lambda t: 9*np.exp(-(t / 5)** 2)
    noisy_func = lambda t: func(t)+(0.05*func(t))*np.random.randn(t.shape[0])
    noisy_data = noisy_func(t)

    plt.figure()
    plt.plot(t, func(t))
    plt.plot(t, noisy_data, 'o')
    plt.show()


To turn these data points into a function we call the QuTiP :class:`qutip.interpolate.Cubic_Spline` class using the first and last domain time points, ``t[0]`` and ``t[-1]``, respectively, as well as the entire array of data points:


.. plot::
    :context: close-figs

    S = Cubic_Spline(t[0], t[-1], noisy_data)

    plt.figure()
    plt.plot(t, func(t))
    plt.plot(t, noisy_data, 'o')
    plt.plot(t, S(t), lw=2)
    plt.show()


Note that, at present, only equally spaced real or complex data sets can be accommodated.  This cubic spline class ``S`` can now be pasted to any of the ``mesolve``, ``mcsolve``, or ``sesolve`` functions where one would normally input a time-dependent function or string-representation.  Taking the problem from the previous section as an example.  We would make the replacement:

.. code-block:: python

    H = [H0, [H1, '9 * exp(-(t / 5) ** 2)']]

to

.. code-block:: python

    H = [H0, [H1, S]]


When combining interpolating functions with other Python functions or strings, the interpolating class will automatically pick the appropriate method for calling the class.  That is to say that, if for example, you have other time-dependent terms that are given in the string-format, then the cubic spline representation will also be passed in a string-compatible format.  In the string-format, the interpolation function is compiled into c-code, and thus is quite fast.  This is the default method if no other time-dependent terms are present.


.. _time-dynargs:

Accesing the state from solver
==============================

New in QuTiP 4.4

The state of the system, the ket vector or the density matrix,
is available to time-dependent Hamiltonian and collapse operators in ``args``.
Some keys of the argument dictionary are understood by the solver to be values
to be updated with the evolution of the system.
The state can be obtained in 3 forms: ``Qobj``, vector (1d ``np.array``), matrix (2d ``np.array``),
expectation values and collapse can also be obtained.

+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
|                   | Preparation             | usage                | Notes                                                            |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
| state as Qobj     | ``name+"=Qobj":psi0``   | ``psi_t=args[name]`` | The ket or density matrix as a Qobj with ``psi0``'s dimensions   |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
| state as matrix   | ``name+"=mat":psi0``    | ``mat_t=args[name]`` | The state as a matrix, equivalent to ``state.full()``            |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
| state as vector   | ``name+"=vec":psi0``    | ``vec_t=args[name]`` | The state as a vector, equivalent to ``state.full().ravel('F')`` |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
| expectation value | ``name+"=expect":O``    | ``e=args[name]``     | Expectation value of the operator ``O``, either                  |
|                   |                         |                      | :math:`\left<\psi(t)|O|\psi(t)\right>`                           |
|                   |                         |                      | or :math:`\rm{tr}\left(O \rho(t)\right)`                         |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+
| collpases         | ``name+"=collapse":[]`` | ``col=args[name]``   | List of collapse,                                                |
|                   |                         |                      | each collapse is a tuple of the pair ``(time, which)``           |
|                   |                         |                      | ``which`` being the indice of the collapse operator.             |
|                   |                         |                      | ``mcsolve`` only.                                                |
+-------------------+-------------------------+----------------------+------------------------------------------------------------------+

Here ``psi0`` is the initial value used for tests before the evolution begins.
:func:`qutip.brmesolve` does not support these arguments.

Reusing Time-Dependent Hamiltonian Data
=======================================

.. note:: This section covers a specialized topic and may be skipped if you are new to QuTiP.

When repeatedly simulating a system where only the time-dependent variables, or initial state change, it is possible to reuse the Hamiltonian data stored in QuTiP and there by avoid spending time needlessly preparing the Hamiltonian and collapse terms for simulation.  To turn on the the reuse features, we must pass a :class:`qutip.Options` object with the ``rhs_reuse`` flag turned on.  Instructions on setting flags are found in :ref:`Options`.  For example, we can do

.. plot::
    :context: close-figs

    H = [H0, [H1, 'A * exp(-(t / sig) ** 2)']]
    args = {'A': 9, 'sig': 5}
    output = mcsolve(H, psi0, times, c_ops, [a.dag()*a], args=args)
    opts = Options(rhs_reuse=True)
    args = {'A': 10, 'sig': 3}
    output = mcsolve(H, psi0, times, c_ops, [a.dag()*a], args=args, options=opts)

The second call to :func:`qutip.mcsolve` does not reorganize the data, and in the case of the string format, does not recompile the Cython code.  For the small system here, the savings in computation time is quite small, however, if you need to call the solvers many times for different parameters, this savings will obviously start to add up.


.. _time-parallel:

Running String-Based Time-Dependent Problems using Parfor
==========================================================

.. note:: This section covers a specialized topic and may be skipped if you are new to QuTiP.

In this section we discuss running string-based time-dependent problems using the :func:`qutip.parfor` function.  As the :func:`qutip.mcsolve` function is already parallelized, running string-based time dependent problems inside of parfor loops should be restricted to the :func:`qutip.mesolve` function only. When using the string-based format, the system Hamiltonian and collapse operators are converted into C code with a specific file name that is automatically genrated, or supplied by the user via the ``rhs_filename`` property of the :class:`qutip.Options` class. Because the :func:`qutip.parfor` function uses the built-in Python multiprocessing functionality, in calling the solver inside a parfor loop, each thread will try to generate compiled code with the same file name, leading to a crash.  To get around this problem you can call the :func:`qutip.rhs_generate` function to compile simulation into C code before calling parfor.  You **must** then set the :class:`qutip.Odedata` object ``rhs_reuse=True`` for all solver calls inside the parfor loop that indicates that a valid C code file already exists and a new one should not be generated.  As an example, we will look at the Landau-Zener-Stuckelberg interferometry example that can be found in the notebook "Time-dependent master equation: Landau-Zener-Stuckelberg inteferometry" in the tutorials section of the QuTiP web site.

To set up the problem, we run the following code:

.. plot::
   :context:

   delta = 0.1  * 2 * np.pi  # qubit sigma_x coefficient
   w = 2.0  * 2 * np.pi      # driving frequency
   T = 2 * np.pi / w         # driving period
   gamma1 = 0.00001          # relaxation rate
   gamma2 = 0.005            # dephasing  rate

   eps_list = np.linspace(-10.0, 10.0, 51) * 2 * np.pi  # epsilon
   A_list = np.linspace(0.0, 20.0, 51) * 2 * np.pi	# Amplitude

   sx = sigmax(); sz = sigmaz(); sm = destroy(2); sn = num(2)

   c_ops = [np.sqrt(gamma1) * sm, np.sqrt(gamma2) * sz]  # relaxation and dephasing
   H0 = -delta / 2.0 * sx
   H1 = [sz, '-eps / 2.0 + A / 2.0 * sin(w * t)']
   H_td = [H0, H1]
   Hargs = {'w': w, 'eps': eps_list[0], 'A': A_list[0]}


where the last code block sets up the problem using a string-based Hamiltonian, and ``Hargs`` is a dictionary of arguments to be passed into the Hamiltonian.  In this example, we are going to use the :func:`qutip.propagator` and :func:`qutip.propagator.propagator_steadystate` to find expectation
values for different values of :math:`\epsilon` and :math:`A` in the
Hamiltonian :math:`H = -\frac{1}{2}\Delta\sigma_x -\frac{1}{2}\epsilon\sigma_z- \frac{1}{2}A\sin(\omega t)`.

We must now tell the :func:`qutip.mesolve` function, that is called by :func:`qutip.propagator` to reuse a
pre-generated Hamiltonian constructed using the :func:`qutip.rhs_generate` command:

.. plot::
   :context:

   opts = Options(rhs_reuse=True)
   rhs_generate(H_td, c_ops, Hargs, name='lz_func')

Here, we have given the generated file a custom name ``lz_func``, however this is not necessary as a generic name will automatically be given.  Now we define the function ``task`` that is called by :func:`qutip.parallel.parfor` with the m-index parallelized in loop over the elements of ``p_mat[m,n]``:

.. code-block:: python

   def task(args):
      m, eps = args
      p_mat_m = np.zeros(len(A_list))
      for n, A in enumerate(A_list):
          # change args sent to solver, w is really a constant though.
          Hargs = {'w': w, 'eps': eps,'A': A}
          U = propagator(H_td, T, c_ops, Hargs, opts) #<- IMPORTANT LINE
          rho_ss = propagator_steadystate(U)
          p_mat_m[n] = expect(sn, rho_ss)
      return [m, p_mat_m]

Notice the Options ``opts`` in the call to the :func:`qutip.propagator` function.  This is tells the :func:`qutip.mesolve` function used in the propagator to call the pre-generated file ``lz_func``. If this were missing then the routine would fail.
.. _options:

*********************************************
Setting Options for the Dynamics Solvers
*********************************************

.. testsetup:: [dynamics_options]

   from qutip import Options

   import numpy as np

Occasionally it is necessary to change the built in parameters of the dynamics solvers used by for example the :func:`qutip.mesolve` and :func:`qutip.mcsolve` functions.  The options for all dynamics solvers may be changed by using the Options class :class:`qutip.solver.Options`.

.. testcode:: [dynamics_options]

   options = Options()

the properties and default values of this class can be view via the `print` function:

.. testcode:: [dynamics_options]

   print(options)

**Output**:

.. testoutput:: [dynamics_options]
  :options: +NORMALIZE_WHITESPACE

  Options:
  -----------
  atol:              1e-08
  rtol:              1e-06
  method:            adams
  order:             12
  nsteps:            1000
  first_step:        0
  min_step:          0
  max_step:          0
  tidy:              True
  num_cpus:          2
  norm_tol:          0.001
  norm_steps:        5
  rhs_filename:      None
  rhs_reuse:         False
  seeds:             0
  rhs_with_state:    False
  average_expect:    True
  average_states:    False
  ntraj:             500
  store_states:      False
  store_final_state: False

These properties are detailed in the following table.  Assuming ``options = Options()``:

.. cssclass:: table-striped

+-----------------------------+-----------------+----------------------------------------------------------------+
| Property                    | Default setting | Description                                                    |
+=============================+=================+================================================================+
| options.atol                | 1e-8            | Absolute tolerance                                             |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.rtol                | 1e-6            | Relative tolerance                                             |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.method              | 'adams'         | Solver method.  Can be 'adams' (non-stiff) or 'bdf' (stiff)    |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.order               | 12              | Order of solver.  Must be <=12 for 'adams' and <=5 for 'bdf'   |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.nsteps              | 1000            | Max. number of steps to take for each interval                 |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.first_step          | 0               | Size of initial step.  0 = determined automatically by solver. |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.min_step            | 0               | Minimum step size.  0 = determined automatically by solver.    |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.max_step            | 0               | Maximum step size.  0 = determined automatically by solver.    |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.tidy                | True            | Whether to run tidyup function on time-independent Hamiltonian.|
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.store_final_state   | False           | Whether or not to store the final state of the evolution.      |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.store_states        | False           | Whether or not to store the state vectors or density matrices. |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.rhs_filename        | None            | RHS filename when using compiled time-dependent Hamiltonians.  |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.rhs_reuse           | False           | Reuse compiled RHS function.  Useful for repetitive tasks.     |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.rhs_with_state      | False           | Whether or not to include the state in the Hamiltonian         |
|                             |                 | function callback signature.                                   |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.num_cpus            | installed num   | Integer number of cpus used by mcsolve.                        |
|                             | of processors   |                                                                |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.seeds               | None            | Array containing random number seeds for mcsolver.             |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.norm_tol            | 1e-6            | Tolerance used when finding wavefunction norm in mcsolve.      |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.norm_steps          | 5               | Max. number of steps used to find wavefunction's norm to within|
|                             |                 | norm_tol in mcsolve.                                           |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.steady_state_average| False           | Include an estimation of the steady state  in mcsolve.         |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.ntraj               | 500             | Number of trajectories in stochastic solvers.                  |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.average_expect      | True            | Average expectation values over trajectories.                  |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.average_states      | False           | Average of the states over trajectories.                       |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.openmp_threads      | installed num   | Number of OPENMP threads to use.                               |
|                             | of processors   |                                                                |
+-----------------------------+-----------------+----------------------------------------------------------------+
| options.use_openmp          | None            | Use OPENMP for sparse matrix vector multiplication.            |
+-----------------------------+-----------------+----------------------------------------------------------------+

As an example, let us consider changing the number of processors used, turn the GUI off, and strengthen the absolute tolerance.  There are two equivalent ways to do this using the Options class.  First way,

.. testcode:: [dynamics_options]

    options = Options()
    options.num_cpus = 3
    options.atol = 1e-10

or one can use an inline method,

.. testcode:: [dynamics_options]

    options = Options(num_cpus=4, atol=1e-10)

Note that the order in which you input the options does not matter.  Using either method, the resulting `options` variable is now:

.. testcode:: [dynamics_options]

  print(options)

**Output**:

.. testoutput:: [dynamics_options]
  :options: +NORMALIZE_WHITESPACE

  Options:
  -----------
  atol:              1e-10
  rtol:              1e-06
  method:            adams
  order:             12
  nsteps:            1000
  first_step:        0
  min_step:          0
  max_step:          0
  tidy:              True
  num_cpus:          4
  norm_tol:          0.001
  norm_steps:        5
  rhs_filename:      None
  rhs_reuse:         False
  seeds:             0
  rhs_with_state:    False
  average_expect:    True
  average_states:    False
  ntraj:             500
  store_states:      False
  store_final_state: False



To use these new settings we can use the keyword argument ``options`` in either the func:`qutip.mesolve` and :func:`qutip.mcsolve` function.  We can modify the last example as::

    >>> mesolve(H0, psi0, tlist, c_op_list, [sigmaz()], options=options)
    >>> mesolve(hamiltonian_t, psi0, tlist, c_op_list, [sigmaz()], H_args, options=options)

or::

    >>> mcsolve(H0, psi0, tlist, ntraj,c_op_list, [sigmaz()], options=options)
    >>> mcsolve(hamiltonian_t, psi0, tlist, ntraj, c_op_list, [sigmaz()], H_args, options=options)
.. _monte:

*******************************************
Monte Carlo Solver
*******************************************

.. _monte-intro:

Introduction
=============

Where as the density matrix formalism describes the ensemble average over many identical realizations of a quantum system, the Monte Carlo (MC), or quantum-jump approach to wave function evolution, allows for simulating an individual realization of the system dynamics.  Here, the environment is continuously monitored, resulting in a series of quantum jumps in the system wave function, conditioned on the increase in information gained about the state of the system via the environmental measurements.  In general, this evolution is governed by the Schrödinger equation with a **non-Hermitian** effective Hamiltonian

.. math::
	:label: heff

	H_{\rm eff}=H_{\rm sys}-\frac{i\hbar}{2}\sum_{i}C^{+}_{n}C_{n},

where again, the :math:`C_{n}` are collapse operators, each corresponding to a separate irreversible process with rate :math:`\gamma_{n}`.  Here, the strictly negative non-Hermitian portion of Eq. :eq:`heff` gives rise to a reduction in the norm of the wave function, that to first-order in a small time :math:`\delta t`, is given by :math:`\left<\psi(t+\delta t)|\psi(t+\delta t)\right>=1-\delta p` where

.. math::
	:label: jump

	\delta p =\delta t \sum_{n}\left<\psi(t)|C^{+}_{n}C_{n}|\psi(t)\right>,

and :math:`\delta t` is such that :math:`\delta p \ll 1`.  With a probability of remaining in the state :math:`\left|\psi(t+\delta t)\right>` given by :math:`1-\delta p`, the corresponding quantum jump probability is thus Eq. :eq:`jump`.  If the environmental measurements register a quantum jump, say via the emission of a photon into the environment, or a change in the spin of a quantum dot, the wave function undergoes a jump into a state defined by projecting :math:`\left|\psi(t)\right>` using the collapse operator :math:`C_{n}` corresponding to the measurement

.. math::
	:label: project

	\left|\psi(t+\delta t)\right>=C_{n}\left|\psi(t)\right>/\left<\psi(t)|C_{n}^{+}C_{n}|\psi(t)\right>^{1/2}.

If more than a single collapse operator is present in Eq. :eq:`heff`, the probability of collapse due to the :math:`i\mathrm{th}`-operator :math:`C_{i}` is given by

.. math::
	:label: pcn

	P_{i}(t)=\left<\psi(t)|C_{i}^{+}C_{i}|\psi(t)\right>/\delta p.

Evaluating the MC evolution to first-order in time is quite tedious.  Instead, QuTiP uses the following algorithm to simulate a single realization of a quantum system.  Starting from a pure state :math:`\left|\psi(0)\right>`:

- **Ia:** Choose a random number :math:`r_1` between zero and one, representing the probability that a quantum jump occurs.

- **Ib:** Choose a random number :math:`r_2` between zero and one, used to select which collapse operator was responsible for the jump.

- **II:** Integrate the Schrödinger equation, using the effective Hamiltonian :eq:`heff` until a time :math:`\tau` such that the norm of the wave function satisfies :math:`\left<\psi(\tau)\right.\left|\psi(\tau)\right> = r_1`, at which point a jump occurs.

- **III:** The resultant jump projects the system at time :math:`\tau` into one of the renormalized states given by Eq. :eq:`project`.  The corresponding collapse operator :math:`C_{n}` is chosen such that :math:`n` is the smallest integer satisfying:

.. math::
    :label: mc3

    \sum_{i=1}^{n} P_{n}(\tau) \ge r_2

where the individual :math:`P_{n}` are given by Eq. :eq:`pcn`.  Note that the left hand side of Eq. :eq:`mc3` is, by definition, normalized to unity.

- **IV:** Using the renormalized state from step III as the new initial condition at time :math:`\tau`, draw a new random number, and repeat the above procedure until the final simulation time is reached.


.. _monte-qutip:

Monte Carlo in QuTiP
====================

In QuTiP, Monte Carlo evolution is implemented with the :func:`qutip.mcsolve` function. It takes nearly the same arguments as the :func:`qutip.mesolve`
function for master-equation evolution, except that the initial state must be a ket vector, as oppose to a density matrix, and there is an optional keyword parameter ``ntraj`` that defines the number of stochastic trajectories to be simulated.  By default, ``ntraj=500`` indicating that 500 Monte Carlo trajectories will be performed.

To illustrate the use of the Monte Carlo evolution of quantum systems in QuTiP, let's again consider the case of a two-level atom coupled to a leaky cavity. The only differences to the master-equation treatment is that in this case we invoke the :func:`qutip.mcsolve` function instead of :func:`qutip.mesolve`

.. plot::
    :context:

    times = np.linspace(0.0, 10.0, 200)
    psi0 = tensor(fock(2, 0), fock(10, 5))
    a  = tensor(qeye(2), destroy(10))
    sm = tensor(destroy(2), qeye(10))
    H = 2*np.pi*a.dag()*a + 2*np.pi*sm.dag()*sm + 2*np.pi*0.25*(sm*a.dag() + sm.dag()*a)
    data = mcsolve(H, psi0, times, [np.sqrt(0.1) * a], [a.dag() * a, sm.dag() * sm])

    plt.figure()
    plt.plot(times, data.expect[0], times, data.expect[1])
    plt.title('Monte Carlo time evolution')
    plt.xlabel('Time')
    plt.ylabel('Expectation values')
    plt.legend(("cavity photon number", "atom excitation probability"))
    plt.show()

.. guide-dynamics-mc1:

The advantage of the Monte Carlo method over the master equation approach is that only the state vector is required to be kept in the computers memory, as opposed to the entire density matrix. For large quantum system this becomes a significant advantage, and the Monte Carlo solver is therefore generally recommended for such systems. For example, simulating a Heisenberg spin-chain consisting of 10 spins with random parameters and initial states takes almost 7 times longer using the master equation rather than Monte Carlo approach with the default number of trajectories running on a quad-CPU machine.  Furthermore, it takes about 7 times the memory as well. However, for small systems, the added overhead of averaging a large number of stochastic trajectories to obtain the open system dynamics, as well as starting the multiprocessing functionality, outweighs the benefit of the minor (in this case) memory saving. Master equation methods are therefore generally more efficient when Hilbert space sizes are on the order of a couple of hundred states or smaller.

Like the master equation solver :func:`qutip.mesolve`, the Monte Carlo solver returns a :class:`qutip.solver.Result` object consisting of expectation values, if the user has defined expectation value operators in the 5th argument to ``mcsolve``, or state vectors if no expectation value operators are given.  If state vectors are returned, then the :class:`qutip.solver.Result` returned by :func:`qutip.mcsolve` will be an array of length ``ntraj``, with each element containing an array of ket-type qobjs with the same number of elements as ``times``.  Furthermore, the output :class:`qutip.solver.Result` object will also contain a list of times at which collapse occurred, and which collapse operators did the collapse, in the ``col_times`` and ``col_which`` properties, respectively.


.. _monte-ntraj:

Changing the Number of Trajectories
-----------------------------------

As mentioned earlier, by default, the ``mcsolve`` function runs 500 trajectories.  This value was chosen because it gives good accuracy, Monte Carlo errors scale as :math:`1/n` where :math:`n` is the number of trajectories, and simultaneously does not take an excessive amount of time to run.  However, like many other options in QuTiP you are free to change the number of trajectories to fit your needs.  If we want to run 1000 trajectories in the above example, we can simply modify the call to ``mcsolve`` like:

.. plot::
    :context: close-figs

    data = mcsolve(H, psi0, times, [np.sqrt(0.1) * a], [a.dag() * a, sm.dag() * sm], ntraj=1000)

where we have added the keyword argument ``ntraj=1000`` at the end of the inputs.  Now, the Monte Carlo solver will calculate expectation values for both operators, ``a.dag() * a, sm.dag() * sm`` averaging over 1000 trajectories.  Sometimes one is also interested in seeing how the Monte Carlo trajectories converge to the master equation solution by calculating expectation values over a range of trajectory numbers.  If, for example, we want to average over 1, 10, 100, and 1000 trajectories, then we can input this into the solver using:

.. plot::
    :context:

    ntraj = [1, 10, 100, 1000]

Keep in mind that the input list must be in ascending order since the total number of trajectories run by ``mcsolve`` will be calculated using the last element of ``ntraj``.  In this case, we need to use an extra index when getting the expectation values from the :class:`qutip.solver.Result` object returned by ``mcsolve``.  In the above example using:

.. plot::
    :context:

    data = mcsolve(H, psi0, times, [np.sqrt(0.1) * a], [a.dag() * a, sm.dag() * sm], ntraj=[1, 10, 100, 1000])

we can extract the relevant expectation values using:

.. plot::
    :context:

    expt1 = data.expect[0]
    expt10 = data.expect[1]
    expt100 = data.expect[2]
    expt1000 = data.expect[3]

The Monte Carlo solver also has many available options that can be set using the :func:`qutip.solver.Options` class as discussed in :ref:`options`.


.. _monte-reuse:

Reusing Hamiltonian Data
------------------------

.. note:: This section covers a specialized topic and may be skipped if you are new to QuTiP.

In order to solve a given simulation as fast as possible, the solvers in QuTiP take the given input operators and break them down into simpler components before passing them on to the ODE solvers.  Although these operations are reasonably fast, the time spent organizing data can become appreciable when repeatedly solving a system over, for example, many different initial conditions. In cases such as this, the Hamiltonian and other operators may be reused after the initial configuration, thus speeding up calculations.  Note that, unless you are planning to reuse the data many times, this functionality will not be very useful.

To turn on the "reuse" functionality we must set the ``rhs_reuse=True`` flag in the :func:`qutip.solver.Options`:

.. plot::
    :context:

    options = Options(rhs_reuse=True)

A full account of this feature is given in :ref:`options`.  Using the previous example, we will calculate the dynamics for two different initial states, with the Hamiltonian data being reused on the second call

.. plot::
    :context:

    times = np.linspace(0.0, 10.0, 200)
    psi0 = tensor(fock(2, 0), fock(10, 5))
    a  = tensor(qeye(2), destroy(10))
    sm = tensor(destroy(2), qeye(10))

    H = 2*np.pi*a.dag()*a + 2*np.pi*sm.dag()*sm + 2*np.pi*0.25*(sm*a.dag() + sm.dag()*a)
    data1 = mcsolve(H, psi0, times, [np.sqrt(0.1) * a], [a.dag() * a, sm.dag() * sm])
    psi1 = tensor(fock(2, 0), coherent(10, 2 - 1j))
    opts = Options(rhs_reuse=True) # Run a second time, reusing RHS
    data2 = mcsolve(H, psi1, times, [np.sqrt(0.1) * a], [a.dag() * a, sm.dag() * sm], options=opts)

    plt.figure()
    plt.plot(times, data1.expect[0], times, data1.expect[1], lw=2)
    plt.plot(times, data2.expect[0], '--', times, data2.expect[1], '--', lw=2)
    plt.title('Monte Carlo time evolution')
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Expectation values', fontsize=14)
    plt.legend(("cavity photon number", "atom excitation probability"))
    plt.show()

.. guide-dynamics-mc2:

In addition to the initial state, one may reuse the Hamiltonian data when changing the number of trajectories ``ntraj`` or simulation times ``times``.  The reusing of Hamiltonian data is also supported for time-dependent Hamiltonians.  See :ref:`time` for further details.
.. _solver_result:

********************************************************
Dynamics Simulation Results
********************************************************

.. _solver_result-class:

The solver.Result Class
=======================

Before embarking on simulating the dynamics of quantum systems, we will first look at the data structure used for returning the simulation results to the user. This object is a :func:`qutip.solver.Result` class that stores all the crucial data needed for analyzing and plotting the results of a simulation.  Like the :func:`qutip.Qobj` class, the ``Result`` class has a collection of properties for storing information.  However, in contrast to the ``Qobj`` class, this structure contains no methods, and is therefore nothing but a container object.  A generic ``Result`` object ``result`` contains the following properties for storing simulation data:

.. cssclass:: table-striped

+------------------------+-----------------------------------------------------------------------+
| Property               | Description                                                           |
+========================+=======================================================================+
| ``result.solver``      | String indicating which solver was used to generate the data.         |
+------------------------+-----------------------------------------------------------------------+
| ``result.times``       | List/array of times at which simulation data is calculated.           |
+------------------------+-----------------------------------------------------------------------+
| ``result.expect``      | List/array of expectation values, if requested.                       |
+------------------------+-----------------------------------------------------------------------+
| ``result.states``      | List/array of state vectors/density matrices calculated at ``times``, |
|                        | if requested.                                                         |
+------------------------+-----------------------------------------------------------------------+
| ``result.num_expect``  | The number of expectation value operators in the simulation.          |
+------------------------+-----------------------------------------------------------------------+
| ``result.num_collapse``| The number of collapse operators in the simulation.                   |
+------------------------+-----------------------------------------------------------------------+
| ``result.ntraj``       | Number of Monte Carlo trajectories run.                               |
+------------------------+-----------------------------------------------------------------------+
| ``result.col_times``   | Times at which state collapse occurred. Only for Monte Carlo solver.  |
+------------------------+-----------------------------------------------------------------------+
| ``result.col_which``   | Which collapse operator was responsible for each collapse in          |
|                        | in ``col_times``. Only used by Monte Carlo solver.                    |
+------------------------+-----------------------------------------------------------------------+
| ``result.seeds``       | Seeds used in generating random numbers for Monte Carlo solver.       |
+------------------------+-----------------------------------------------------------------------+


.. _odedata-access:

Accessing Result Data
======================

To understand how to access the data in a Result object we will use an example as a guide, although we do not worry about the simulation details at this stage.  Like all solvers, the Monte Carlo solver used in this example returns an Result object, here called simply ``result``.  To see what is contained inside ``result`` we can use the print function:

.. doctest::
  :options: +SKIP

  >>> print(result)
  Result object with mcsolve data.
  ---------------------------------
  expect = True
  num_expect = 2, num_collapse = 2, ntraj = 500

The first line tells us that this data object was generated from the Monte Carlo solver ``mcsolve`` (discussed in :ref:`monte`).  The next line (not the ``---`` line of course) indicates that this object contains expectation value data.  Finally, the last line gives the number of expectation value and collapse operators used in the simulation, along with the number of Monte Carlo trajectories run.  Note that the number of trajectories ``ntraj`` is only displayed when using the Monte Carlo solver.

Now we have all the information needed to analyze the simulation results.
To access the data for the two expectation values one can do:


.. testcode::
  :skipif: True

  expt0 = result.expect[0]
  expt1 = result.expect[1]

Recall that Python uses C-style indexing that begins with zero (i.e., [0] => 1st collapse operator data). Together with the array of times at which these expectation values are calculated:

.. testcode::
  :skipif: True

  times = result.times

we can plot the resulting expectation values:

.. testcode::
  :skipif: True

  plot(times, expt0, times, expt1)
  show()


State vectors, or density matrices, as well as ``col_times`` and ``col_which``, are accessed in a similar manner, although typically one does not need an index (i.e [0]) since there is only one list for each of these components.  The one exception to this rule is if you choose to output state vectors from the Monte Carlo solver, in which case there are ``ntraj`` number of state vector arrays.

.. _odedata-saving:

Saving and Loading Result Objects
==================================

The main advantage in using the Result class as a data storage object comes from the simplicity in which simulation data can be stored and later retrieved. The :func:`qutip.fileio.qsave` and :func:`qutip.fileio.qload` functions are designed for this task.  To begin, let us save the ``data`` object from the previous section into a file called "cavity+qubit-data" in the current working directory by calling:

.. testcode::
  :skipif: True

  qsave(result, 'cavity+qubit-data')

All of the data results are then stored in a single file of the same name with a ".qu" extension.  Therefore, everything needed to later this data is stored in a single file.  Loading the file is just as easy as saving:

.. doctest::
  :options: +SKIP

  >>> stored_result = qload('cavity+qubit-data')
  Loaded Result object:
  Result object with mcsolve data.
  ---------------------------------
  expect = True
  num_expect = 2, num_collapse = 2, ntraj = 500

where ``stored_result`` is the new name of the Result object.  We can then extract the data and plot in the same manner as before:

.. testcode::
    :skipif: True

    expt0 = stored_result.expect[0]
    expt1 = stored_result.expect[1]
    times = stored_result.times
    plot(times, expt0, times, expt1)
    show()

Also see :ref:`saving` for more information on saving quantum objects, as well as arrays for use in other programs.
.. _bloch_redfield:

******************************
Bloch-Redfield master equation
******************************


.. plot::
      :include-source: False

      import pylab as plt
      from scipy import *
      from qutip import *
      import numpy as np

.. _bloch-redfield-intro:

Introduction
============

The Lindblad master equation introduced earlier is constructed so that it describes a physical evolution of the density matrix (i.e., trace and positivity preserving), but it does not provide a connection to any underlaying microscopic physical model. The Lindblad operators (collapse operators) describe phenomenological processes, such as for example dephasing and spin flips, and the rates of these processes are arbitrary parameters in the model. In many situations the collapse operators and their corresponding rates have clear physical interpretation, such as dephasing and relaxation rates, and in those cases the Lindblad master equation is usually the method of choice.

However, in some cases, for example systems with varying energy biases and eigenstates and that couple to an environment in some well-defined manner (through a physically motivated system-environment interaction operator), it is often desirable to derive the master equation from more fundamental physical principles, and relate it to for example the noise-power spectrum of the environment.

The Bloch-Redfield formalism is one such approach to derive a master equation from a microscopic system. It starts from a combined system-environment perspective, and derives a perturbative master equation for the system alone, under the assumption of weak system-environment coupling. One advantage of this approach is that the dissipation processes and rates are obtained directly from the properties of the environment. On the downside, it does not intrinsically guarantee that the resulting master equation unconditionally preserves the physical properties of the density matrix (because it is a perturbative method). The Bloch-Redfield master equation must therefore be used with care, and the assumptions made in the derivation must be honored. (The Lindblad master equation is in a sense more robust -- it always results in a physical density matrix -- although some collapse operators might not be physically justified). For a full derivation of the Bloch Redfield master equation, see e.g. [Coh92]_ or [Bre02]_. Here we present only a brief version of the derivation, with the intention of introducing the notation and how it relates to the implementation in QuTiP.

.. _bloch-redfield-derivation:


Brief Derivation and Definitions
================================

The starting point of the Bloch-Redfield formalism is the total Hamiltonian for the system and the environment (bath): :math:`H = H_{\rm S} + H_{\rm B} + H_{\rm I}`, where :math:`H` is the total system+bath Hamiltonian, :math:`H_{\rm S}` and :math:`H_{\rm B}` are the system and bath Hamiltonians, respectively, and :math:`H_{\rm I}` is the interaction Hamiltonian.

The most general form of a master equation for the system dynamics is obtained by tracing out the bath from the von-Neumann equation of motion for the combined system (:math:`\dot\rho = -i\hbar^{-1}[H, \rho]`). In the interaction picture the result is

.. math::
   :label: br-nonmarkovian-form-one

    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^t d\tau\;  {\rm Tr}_B [H_I(t), [H_I(\tau), \rho_S(\tau)\otimes\rho_B]],

where the additional assumption that the total system-bath density matrix can be factorized as :math:`\rho(t) \approx \rho_S(t) \otimes \rho_B`. This assumption is known as the Born approximation, and it implies that there never is any entanglement between the system and the bath, neither in the initial state nor at any time during the evolution. *It is justified for weak system-bath interaction.*

The master equation :eq:`br-nonmarkovian-form-one` is non-Markovian, i.e., the change in the density matrix at a time :math:`t` depends on states at all times :math:`\tau < t`, making it intractable to solve both theoretically and numerically. To make progress towards a manageable master equation, we now introduce the Markovian approximation, in which :math:`\rho(s)` is replaced by :math:`\rho(t)` in Eq. :eq:`br-nonmarkovian-form-one`. The result is the Redfield equation

.. math::
   :label: br-nonmarkovian-form-two

    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^t d\tau\; {\rm Tr}_B [H_I(t), [H_I(\tau), \rho_S(t)\otimes\rho_B]],

which is local in time with respect the density matrix, but still not Markovian since it contains an implicit dependence on the initial state. By extending the integration to infinity and substituting :math:`\tau \rightarrow t-\tau`, a fully Markovian master equation is obtained:

.. math::
   :label: br-markovian-form

    \frac{d}{dt}\rho_S(t) = - \hbar^{-2}\int_0^\infty d\tau\; {\rm Tr}_B [H_I(t), [H_I(t-\tau), \rho_S(t)\otimes\rho_B]].

The two Markovian approximations introduced above are valid if the time-scale with which the system dynamics changes is large compared to the time-scale with which correlations in the bath decays (corresponding to a "short-memory" bath, which results in Markovian system dynamics).

The master equation :eq:`br-markovian-form` is still on a too general form to be suitable for numerical implementation. We therefore assume that the system-bath interaction takes the form :math:`H_I = \sum_\alpha A_\alpha \otimes B_\alpha` and where :math:`A_\alpha` are system operators and :math:`B_\alpha` are bath operators. This allows us to write master equation in terms of system operators and bath correlation functions:

.. math::

    \frac{d}{dt}\rho_S(t) =
    -\hbar^{-2}
    \sum_{\alpha\beta}
    \int_0^\infty d\tau\;
    \left\{
    g_{\alpha\beta}(\tau) \left[A_\alpha(t)A_\beta(t-\tau)\rho_S(t) - A_\alpha(t-\tau)\rho_S(t)A_\beta(t)\right]
    \right. \nonumber\\
    \left.
    g_{\alpha\beta}(-\tau) \left[\rho_S(t)A_\alpha(t-\tau)A_\beta(t) - A_\alpha(t)\rho_S(t)A_\beta(t-\tau)\right]
    \right\},

where :math:`g_{\alpha\beta}(\tau) = {\rm Tr}_B\left[B_\alpha(t)B_\beta(t-\tau)\rho_B\right] = \left<B_\alpha(\tau)B_\beta(0)\right>`, since the bath state :math:`\rho_B` is a steady state.

In the eigenbasis of the system Hamiltonian, where :math:`A_{mn}(t) = A_{mn} e^{i\omega_{mn}t}`, :math:`\omega_{mn} = \omega_m - \omega_n` and :math:`\omega_m` are the eigenfrequencies corresponding the eigenstate :math:`\left|m\right>`, we obtain in matrix form in the Schrödinger picture

.. math::

    \frac{d}{dt}\rho_{ab}(t)
    =
    -i\omega_{ab}\rho_{ab}(t)
    -\hbar^{-2}
    \sum_{\alpha,\beta}
    \sum_{c,d}^{\rm sec}
    \int_0^\infty d\tau\;
    \left\{
    g_{\alpha\beta}(\tau)
    \left[\delta_{bd}\sum_nA^\alpha_{an}A^\beta_{nc}e^{i\omega_{cn}\tau}
    -
    A^\alpha_{ac} A^\beta_{db} e^{i\omega_{ca}\tau}
    \right]
    \right. \nonumber\\
    +
    \left.
    g_{\alpha\beta}(-\tau)
    \left[\delta_{ac}\sum_n A^\alpha_{dn}A^\beta_{nb} e^{i\omega_{nd}\tau}
    -
    A^\alpha_{ac}A^\beta_{db}e^{i\omega_{bd}\tau}
    \right]
    \right\} \rho_{cd}(t),
    \nonumber\\

where the "sec" above the summation symbol indicate summation of the secular terms which satisfy :math:`|\omega_{ab}-\omega_{cd}| \ll \tau_ {\rm decay}`. This is an almost-useful form of the master equation. The final step before arriving at the form of the Bloch-Redfield master equation that is implemented in QuTiP, involves rewriting the bath correlation function :math:`g(\tau)` in terms of the noise-power spectrum of the environment :math:`S(\omega) = \int_{-\infty}^\infty d\tau e^{i\omega\tau} g(\tau)`:

.. math::
   :label: br-nonmarkovian-form-four

    \int_0^\infty d\tau\; g_{\alpha\beta}(\tau) e^{i\omega\tau} = \frac{1}{2}S_{\alpha\beta}(\omega) + i\lambda_{\alpha\beta}(\omega),

where :math:`\lambda_{ab}(\omega)` is an energy shift that is neglected here. The final form of the Bloch-Redfield master equation is


.. math::
    :label: br-final

    \frac{d}{dt}\rho_{ab}(t)
    =
    -i\omega_{ab}\rho_{ab}(t)
    +
    \sum_{c,d}^{\rm sec}R_{abcd}\rho_{cd}(t),

where

.. math::
   :label: br-nonmarkovian-form-five

    R_{abcd} =  -\frac{\hbar^{-2}}{2} \sum_{\alpha,\beta}
    \left\{
    \delta_{bd}\sum_nA^\alpha_{an}A^\beta_{nc}S_{\alpha\beta}(\omega_{cn})
    -
    A^\alpha_{ac} A^\beta_{db} S_{\alpha\beta}(\omega_{ca})
    \right. \nonumber\\
    +
    \left.
    \delta_{ac}\sum_n A^\alpha_{dn}A^\beta_{nb} S_{\alpha\beta}(\omega_{dn})
    -
    A^\alpha_{ac}A^\beta_{db} S_{\alpha\beta}(\omega_{db})
    \right\},

is the Bloch-Redfield tensor.

The Bloch-Redfield master equation in the form Eq. :eq:`br-final` is suitable for numerical implementation. The input parameters are the system Hamiltonian :math:`H`, the system operators through which the environment couples to the system :math:`A_\alpha`, and the noise-power spectrum :math:`S_{\alpha\beta}(\omega)` associated with each system-environment interaction term.

To simplify the numerical implementation we assume that :math:`A_\alpha` are Hermitian and that cross-correlations between different environment operators vanish, so that the final expression for the Bloch-Redfield tensor that is implemented in QuTiP is

.. math::
   :label: br-tensor

    R_{abcd} =  -\frac{\hbar^{-2}}{2} \sum_{\alpha}
    \left\{
    \delta_{bd}\sum_nA^\alpha_{an}A^\alpha_{nc}S_{\alpha}(\omega_{cn})
    -
    A^\alpha_{ac} A^\alpha_{db} S_{\alpha}(\omega_{ca})
    \right. \nonumber\\
    +
    \left.
    \delta_{ac}\sum_n A^\alpha_{dn}A^\alpha_{nb} S_{\alpha}(\omega_{dn})
    -
    A^\alpha_{ac}A^\alpha_{db} S_{\alpha}(\omega_{db})
    \right\}.


.. _bloch-redfield-qutip:

Bloch-Redfield master equation in QuTiP
=======================================



In QuTiP, the Bloch-Redfield tensor Eq. :eq:`br-tensor` can be calculated using the function :func:`qutip.bloch_redfield.bloch_redfield_tensor`. It takes two mandatory arguments: The system Hamiltonian :math:`H`, a nested list of operator  :math:`A_\alpha`, spectral density functions :math:`S_\alpha(\omega)` pairs that characterize the coupling between system and bath. The spectral density functions are Python callback functions that takes the (angular) frequency as a single argument.

To illustrate how to calculate the Bloch-Redfield tensor, let's consider a two-level atom

.. math::
   :label: qubit

    H = -\frac{1}{2}\Delta\sigma_x - \frac{1}{2}\epsilon_0\sigma_z


.. testcode:: [dynamics-br]

    delta = 0.2 * 2*np.pi
    eps0 = 1.0 * 2*np.pi
    gamma1 = 0.5

    H = - delta/2.0 * sigmax() - eps0/2.0 * sigmaz()

    def ohmic_spectrum(w):
      if w == 0.0: # dephasing inducing noise
        return gamma1
      else: # relaxation inducing noise
        return gamma1 / 2 * (w / (2 * np.pi)) * (w > 0.0)


    R, ekets = bloch_redfield_tensor(H, [[sigmax(), ohmic_spectrum]])

    print(R)

**Output**:

.. testoutput:: [dynamics-br]

    Quantum object: dims = [[[2], [2]], [[2], [2]]], shape = (4, 4), type = super, isherm = False
    Qobj data =
    [[ 0.        +0.j         0.        +0.j         0.        +0.j
       0.24514517+0.j       ]
     [ 0.        +0.j        -0.16103412-6.4076169j  0.        +0.j
       0.        +0.j       ]
     [ 0.        +0.j         0.        +0.j        -0.16103412+6.4076169j
       0.        +0.j       ]
     [ 0.        +0.j         0.        +0.j         0.        +0.j
      -0.24514517+0.j       ]]

Note that it is also possible to add Lindblad dissipation superoperators in the Bloch-Refield tensor by passing the operators via the ``c_ops`` keyword argument like you would in the :func:`qutip.mesolve` or :func:`qutip.mcsolve` functions. For convenience, the function :func:`qutip.bloch_redfield.bloch_redfield_tensor` also returns a list of eigenkets `ekets`, since they are calculated in the process of calculating the Bloch-Redfield tensor `R`, and the `ekets` are usually needed again later when transforming operators between the computational basis and the eigenbasis.


.. plot::
    :context:
    :include-source: False

    delta = 0.2 * 2*np.pi
    eps0 = 1.0 * 2*np.pi
    gamma1 = 0.5

    H = - delta/2.0 * sigmax() - eps0/2.0 * sigmaz()

    def ohmic_spectrum(w):
      if w == 0.0: # dephasing inducing noise
        return gamma1
      else: # relaxation inducing noise
        return gamma1 / 2 * (w / (2 * np.pi)) * (w > 0.0)

    R, ekets = bloch_redfield_tensor(H, [[sigmax(), ohmic_spectrum]])

The evolution of a wavefunction or density matrix, according to the Bloch-Redfield master equation :eq:`br-final`, can be calculated using the QuTiP function :func:`qutip.bloch_redfield.bloch_redfield_solve`. It takes five mandatory arguments: the Bloch-Redfield tensor ``R``, the list of eigenkets ``ekets``, the initial state ``psi0`` (as a ket or density matrix), a list of times ``tlist`` for which to evaluate the expectation values, and a list of operators ``e_ops`` for which to evaluate the expectation values at each time step defined by `tlist`. For example, to evaluate the expectation values of the :math:`\sigma_x`, :math:`\sigma_y`, and :math:`\sigma_z` operators for the example above, we can use the following code:

.. plot::
    :context:

    tlist = np.linspace(0, 15.0, 1000)

    psi0 = rand_ket(2)

    e_ops = [sigmax(), sigmay(), sigmaz()]

    expt_list = bloch_redfield_solve(R, ekets, psi0, tlist, e_ops)

    sphere = Bloch()

    sphere.add_points([expt_list[0], expt_list[1], expt_list[2]])

    sphere.vector_color = ['r']

    sphere.add_vectors(np.array([delta, 0, eps0]) / np.sqrt(delta ** 2 + eps0 ** 2))

    sphere.make_sphere()

The two steps of calculating the Bloch-Redfield tensor and evolving according to the corresponding master equation can be combined into one by using the function :func:`qutip.bloch_redfield.brmesolve`, which takes same arguments as :func:`qutip.mesolve` and :func:`qutip.mcsolve`, save for the additional nested list of operator-spectrum pairs that is called ``a_ops``.

.. plot::
    :context: close-figs

    output = brmesolve(H, psi0, tlist, a_ops=[[sigmax(),ohmic_spectrum]], e_ops=e_ops)

where the resulting `output` is an instance of the class :class:`qutip.solver.Result`.


.. _td-bloch-redfield:

Time-dependent Bloch-Redfield Dynamics
=======================================

.. warning::

    It takes ~3-5 seconds (~30 if using Visual Studio) to compile a time-dependent Bloch-Redfield problem.  Therefore,
    if you are doing repeated simulations by varying parameters, then it is best to pass
    ``options = Options(rhs_reuse=True)`` to the solver.

If you have not done so already, please read the section: :ref:`time`.

As we have already discussed, the Bloch-Redfield master equation requires transforming into the eigenbasis of the system Hamiltonian.  For time-independent systems, this transformation need only be done once.  However, for time-dependent systems, one must move to the instantaneous eigenbasis at each time-step in the evolution, thus greatly increasing the computational complexity of the dynamics.  In addition, the requirement for computing all the eigenvalues severely limits the scalability of the method.  Fortunately, this eigen decomposition occurs at the Hamiltonian level, as opposed to the super-operator level, and thus, with efficient programming, one can tackle many systems that are commonly encountered.


The time-dependent Bloch-Redfield solver in QuTiP relies on the efficient numerical computations afforded by the string-based time-dependent format, and Cython compilation.  As such, all the time-dependent terms, and noise power spectra must be expressed in the string format.  To begin, lets consider the previous example, but formatted to call the time-dependent solver:


.. plot::
    :context:

    ohmic = "{gamma1} / 2.0 * (w / (2 * pi)) * (w > 0.0)".format(gamma1=gamma1)

    output = brmesolve(H, psi0, tlist, a_ops=[[sigmax(),ohmic]], e_ops=e_ops)


Although the problem itself is time-independent, the use of a string as the noise power spectrum tells the solver to go into time-dependent mode.  The string is nearly identical to the Python function format, except that we replaced ``np.pi`` with ``pi`` to avoid calling Python in our Cython code, and we have hard coded the ``gamma1`` argument into the string as limitations prevent passing arguments into the time-dependent Bloch-Redfield solver.


For actual time-dependent Hamiltonians, the Hamiltonian itself can be passed into the solver like any other string-based Hamiltonian, as thus we will not discuss this topic further.  Instead, here the focus is on time-dependent bath coupling terms.  To this end, suppose that we have a dissipative harmonic oscillator, where the white-noise dissipation rate decreases exponentially with time :math:`\kappa(t) = \kappa(0)\exp(-t)`.  In the Lindblad or monte-carlo solvers, this could be implemented as a time-dependent collapse operator list ``c_ops = [[a, 'sqrt(kappa*exp(-t))']]``.  In the Bloch-Redfield solver, the bath coupling terms must be Hermitian.  As such, in this example, our coupling operator is the position operator ``a+a.dag()``.  In addition, we do not need the ``sqrt`` operation that occurs in the ``c_ops`` definition.  The complete example, and comparison to the analytic expression is:


.. plot::
    :context: close-figs

    N = 10  # number of basis states to consider

    a = destroy(N)

    H = a.dag() * a

    psi0 = basis(N, 9)  # initial state

    kappa = 0.2  # coupling to oscillator

    a_ops = [[a+a.dag(), '{kappa}*exp(-t)*(w>=0)'.format(kappa=kappa)]]

    tlist = np.linspace(0, 10, 100)

    out = brmesolve(H, psi0, tlist, a_ops, e_ops=[a.dag() * a])

    actual_answer = 9.0 * np.exp(-kappa * (1.0 - np.exp(-tlist)))

    plt.figure()

    plt.plot(tlist, out.expect[0])

    plt.plot(tlist, actual_answer)

    plt.show()


In many cases, the bath-coupling operators can take the form :math:`A = f(t)a + f(t)^* a^{+}`.  In this case, the above format for inputting the ``a_ops`` is not sufficient. Instead, one must construct a nested-list of tuples to specify this time-dependence.  For example consider a white-noise bath that is coupled to an operator of the form ``exp(1j*t)*a + exp(-1j*t)* a.dag()``.  In this example, the ``a_ops`` list would be:

.. plot::
    :context: close-figs

    a_ops = [ [ (a, a.dag()), ('{0} * (w >= 0)'.format(kappa), 'exp(1j*t)', 'exp(-1j*t)') ] ]


where the first tuple element ``(a, a.dag())`` tells the solver which operators make up the full Hermitian coupling operator.  The second tuple ``('{0} * (w >= 0)'.format(kappa), 'exp(1j*t)', 'exp(-1j*t)')``, gives the noise power spectrum, and time-dependence of each operator.  Note that the noise spectrum must always come first in this second tuple. A full example is:

.. plot::
    :context:

    N = 10

    w0 = 1.0 * 2 * np.pi

    g = 0.05 * w0

    kappa = 0.15

    times = np.linspace(0, 25, 1000)

    a = destroy(N)

    H = w0 * a.dag() * a + g * (a + a.dag())

    psi0 = ket2dm((basis(N, 4) + basis(N, 2) + basis(N, 0)).unit())

    a_ops = [[ (a, a.dag()), ('{0} * (w >= 0)'.format(kappa), 'exp(1j*t)', 'exp(-1j*t)') ]]

    e_ops = [a.dag() * a, a + a.dag()]

    res_brme = brmesolve(H, psi0, times, a_ops, e_ops)

    plt.figure()

    plt.plot(times,res_brme.expect[0], label=r'$a^{+}a$')

    plt.plot(times,res_brme.expect[1], label=r'$a+a^{+}$')

    plt.legend()

    plt.show()


Further examples on time-dependent Bloch-Redfield simulations can be found in the online tutorials.
.. _master-piqs:

*********************************
Permutational Invariance
*********************************

.. _master-unitary-piqs:

Permutational Invariant Quantum Solver (PIQS)
=============================================
The *Permutational Invariant Quantum Solver (PIQS)* is a QuTiP module that allows to study the dynamics of an open quantum system consisting of an ensemble of identical qubits that can dissipate through local and collective baths according to a Lindblad master equation.

The Liouvillian of an ensemble of :math:`N` qubits, or two-level systems (TLSs), :math:`\mathcal{D}_{TLS}(\rho)`, can be built using only polynomial – instead of exponential – resources.
This has many applications for the study of realistic quantum optics models of many TLSs and in general as a tool in cavity QED.

Consider a system evolving according to the equation

.. math::
    \dot{\rho} = \mathcal{D}_\text{TLS}(\rho)=-\frac{i}{\hbar}\lbrack H,\rho \rbrack
    +\frac{\gamma_\text{CE}}{2}\mathcal{L}_{J_{-}}[\rho]
    +\frac{\gamma_\text{CD}}{2}\mathcal{L}_{J_{z}}[\rho]
    +\frac{\gamma_\text{CP}}{2}\mathcal{L}_{J_{+}}[\rho]

    +\sum_{n=1}^{N}\left(
    \frac{\gamma_\text{E}}{2}\mathcal{L}_{J_{-,n}}[\rho]
    +\frac{\gamma_\text{D}}{2}\mathcal{L}_{J_{z,n}}[\rho]
    +\frac{\gamma_\text{P}}{2}\mathcal{L}_{J_{+,n}}[\rho]\right)


where :math:`J_{\alpha,n}=\frac{1}{2}\sigma_{\alpha,n}` are SU(2) Pauli spin operators, with :math:`{\alpha=x,y,z}` and :math:`J_{\pm,n}=\sigma_{\pm,n}`. The collective spin operators are :math:`J_{\alpha} = \sum_{n}J_{\alpha,n}` . The Lindblad super-operators are :math:`\mathcal{L}_{A} = 2A\rho A^\dagger - A^\dagger A \rho - \rho A^\dagger A`.

The inclusion of local processes in the dynamics lead to using a Liouvillian space of dimension :math:`4^N`. By exploiting the permutational invariance of identical particles [2-8], the Liouvillian :math:`\mathcal{D}_\text{TLS}(\rho)` can be built as a block-diagonal matrix in the basis of Dicke states :math:`|j, m \rangle`.

The system under study is defined by creating an object of the
:code:`Dicke` class, e.g. simply named
:code:`system`, whose first attribute is

- :code:`system.N`, the number of TLSs of the system :math:`N`.

The rates for collective and local processes are simply defined as

- :code:`collective_emission` defines :math:`\gamma_\text{CE}`, collective (superradiant) emission
- :code:`collective_dephasing` defines :math:`\gamma_\text{CD}`, collective dephasing
- :code:`collective_pumping` defines :math:`\gamma_\text{CP}`, collective pumping.
- :code:`emission` defines :math:`\gamma_\text{E}`, incoherent emission (losses)
- :code:`dephasing` defines :math:`\gamma_\text{D}`, local dephasing
- :code:`pumping`  defines :math:`\gamma_\text{P}`, incoherent pumping.

Then the :code:`system.lindbladian()` creates the total TLS Lindbladian superoperator matrix. Similarly, :code:`system.hamiltonian` defines the TLS hamiltonian of the system :math:`H_\text{TLS}`.

The system's Liouvillian can be built using :code:`system.liouvillian()`. The properties of a Piqs object can be visualized by simply calling
:code:`system`. We give two basic examples on the use of *PIQS*. In the first example the incoherent emission of N driven TLSs is considered.

.. code-block:: python

    from piqs import Dicke
    from qutip import steadystate
    N = 10
    system = Dicke(N, emission = 1, pumping = 2)
    L = system.liouvillian()
    steady = steadystate(L)

For more example of use, see the "Permutational Invariant Lindblad Dynamics" section in the tutorials section of the website, `https://qutip.org/tutorials.html <https://qutip.org/tutorials.html>`_.

.. list-table:: Useful PIQS functions.
   :widths: 25 25 50
   :header-rows: 1

   * - Operators
     - Command
     - Description
   * - Collective spin algebra :math:`J_x,\ J_y,\ J_z`
     - ``jspin(N)``
     - The collective spin algebra  :math:`J_x,\ J_y,\ J_z` for :math:`N` TLSs
   * - Collective spin :math:`J_x`
     - ``jspin(N, "x")``
     - The collective spin operator :math:`Jx`. Requires :math:`N` number of TLSs
   * - Collective spin :math:`J_y`
     - ``jspin(N, "y")``
     - The collective spin operator :math:`J_y`. Requires :math:`N` number of TLSs
   * - Collective spin :math:`J_z`
     - ``jspin(N, "z")``
     - The collective spin operator :math:`J_z`. Requires :math:`N` number of TLSs
   * - Collective spin :math:`J_+`
     - ``jspin(N, "+")``
     - The collective spin operator :math:`J_+`.
   * - Collective spin :math:`J_-`
     - ``jspin(N, "-")``
     - The collective spin operator :math:`J_-`.
   * - Collective spin :math:`J_z` in uncoupled basis
     - ``jspin(N, "z", basis='uncoupled')``
     - The collective spin operator :math:`J_z` in the uncoupled basis of dimension :math:`2^N`.
   * - Dicke state :math:`|j,m\rangle` density matrix
     - ``dicke(N, j, m)``
     - The density matrix for the Dicke state given by :math:`|j,m\rangle`
   * - Excited-state density matrix  in Dicke basis
     - ``excited(N)``
     - The excited state in the Dicke basis
   * - Excited-state density matrix  in uncoupled basis
     - ``excited(N, basis="uncoupled")``
     - The excited state in the uncoupled basis
   * - Ground-state density matrix  in Dicke basis
     - ``ground(N)``
     - The ground state in the Dicke basis
   * - GHZ-state density matrix in the Dicke basis
     - ``ghz(N)``
     - The GHZ-state density matrix in the Dicke (default) basis for N number of TLS
   * - Collapse operators of the ensemble
     - ``Dicke.c_ops()``
     - The collapse operators for the ensemble can be called by the `c_ops` method of the Dicke class.

Note that the mathematical object representing the density matrix of the full system that is manipulated (or obtained from `steadystate`) in the Dicke-basis formalism used here is a *representative of the density matrix*. This *representative object* is of linear size N^2, whereas the full density matrix is defined over a 2^N Hilbert space. In order to calculate nonlinear functions of such density matrix, such as the Von Neumann entropy or the purity, it is necessary to take into account the degeneracy of each block of such block-diagonal density matrix. Note that as long as one calculates expected values of operators, being Tr[A*rho] a *linear* function of `rho`, the *representative density matrix* give straightforwardly the correct result. When a *nonlinear* function of the density matrix needs to be calculated, one needs to weigh each degenerate block correctly; this is taken care by the `dicke_function_trace` in `qutip.piqs`, and the user can use it to define general nonlinear functions that can be described as the trace of a Taylor expandable function. Two nonlinear functions that use `dicke_function_trace` and are already implemented are `purity_dicke`, to calculate the purity of a density matrix in the Dicke basis, and `entropy_vn_dicke`, which can be used to calculate the Von Neumann entropy.

More functions relative to the `qutip.piqs` module can be found at :ref:`apidoc`. Attributes to the :class:`qutip.piqs.Dicke` and :class:`qutip.piqs.Pim` class can also be found there.
Gallery
=======

This is the gallery for QuTiP examples, you can click on the image to see the source code.
.. _development:

*************************
Development Documentation
*************************

This chapter covers the development of QuTiP and its subpackages, including
a roadmap for upcoming releases and ideas for future improvements.

.. toctree::
   :maxdepth: 3

   contributing.rst
   roadmap.rst
   ideas.rst
   docs.rst
   release_distribution.rst
.. _development-contributing:

*********************************
Contributing to QuTiP Development
*********************************

Quick Start
===========

QuTiP is developed through wide collaboration using the ``git`` version-control system, with the main repositories hosted in the `qutip organisation on GitHub <https://github.com/qutip>`_.
You will need to be familiar with ``git`` as a tool, and the `GitHub Flow <https://guides.github.com/introduction/flow/>`_ workflow for branching and making pull requests.
The exact details of environment set-up, build process and testing vary by repository and are discussed below, however in overview, the steps to contribute are:

#. Consider creating an issue on the GitHub page of the relevant repository, describing the change you think should be made and why, so we can discuss details with you and make sure it is appropriate.
#. (If this is your first contribution.) Make a fork of the relevant repository on GitHub and clone it to your local computer.  Also add our copy as a remote (``git remote add qutip https://github.com/qutip/<repo>``)
#. Begin on the ``master`` branch (``git checkout master``), and pull in changes from the main QuTiP repository to make sure you have an up-to-date copy (``git pull qutip master``).
#. Switch to a new ``git`` branch (``git checkout -b <branch-name>``).
#. Make the changes you want to make, then create some commits with short, descriptive names (``git add <files>`` then ``git commit``).
#. Follow the build process for this repository to build the final result so you can check your changes work sensibly.
#. Run the tests for the repository (if it has them).
#. Push the changes to your fork (``git push -u origin <branch-name>``).  You won't be able to push to the main QuTiP repositories directly.
#. Go to the GitHub website for the repository you are contributing to, click on the "Pull Requests" tab, click the "New Pull Request" button, and follow the instructions there.

Once the pull request is created, some members of the QuTiP admin team will review the code to make sure it is suitable for inclusion in the library, to check the programming, and to ensure everything meets our standards.
For some repositories, several automated tests will run whenever you create or modify a pull request; in general these will be the same tests you can run locally, and all tests are required to pass online before your changes are merged.
There may be some feedback and possibly some requested changes.
You can add more commits to address these, and push them to the relevant branch of your fork to update the pull request.

The rest of this document covers programming standards, and particular considerations for some of the more complicated repositories.


.. _contributing-qutip:

Core Library: qutip/qutip
=========================

The core library is in the `qutip/qutip repository on GitHub <https://github.com/qutip/qutip>`_.

Building
--------

Building the core library from source is typically a bit more difficult than simply installing the package for regular use.
You will most likely want to do this in a clean Python environment so that you do not compromise a working installation of a release version, for example by starting from ::

   conda create -n qutip-dev python

:ref:`Complete instructions for the build <install>` are elsewhere in this guide, however beware that you will need to follow the :ref:`installation from source using setuptools section <build-setuptools>`, not the general installation.
You will need all the *build* and *tests* "optional" requirements for the package.
The build requirements can be found in the |pyproject.toml file|_, and the testing requirements are in the ``tests`` key of the ``options.extras_require`` section of |setup.cfg|_.
You will also need the requirements for any optional features you want to test as well.

.. |pyproject.toml file| replace:: ``pyproject.toml`` file
.. _pyproject.toml file: https://github.com/qutip/qutip/blob/master/pyproject.toml
.. |setup.cfg| replace:: ``setup.cfg``
.. _setup.cfg: https://github.com/qutip/qutip/blob/master/setup.cfg

Refer to the main instructions for the most up-to-date version, however as of version 4.6 the requirements can be installed into a conda environment with ::

   conda install setuptools wheel 'numpy>=1.16.6,<1.20' 'scipy>=1.0' 'cython>=0.29.20' packaging 'pytest>=5.2' pytest-rerunfailures

Note that ``qutip`` should *not* be installed with ``conda install``.

.. note::
   If you prefer, you can also use ``pip`` to install all the dependencies.
   We typically recommend ``conda`` when doing main-library development because it is easier to switch low-level packages around like BLAS implementations, but if this doesn't mean anything to you, feel free to use ``pip``.

You will need to make sure you have a functioning C++ compiler to build QuTiP.
If you are on Linux or Mac, this is likely already done for you, however if you are on Windows, refer to the :ref:`Windows installation <install-on-windows>` section of the installation guide.

The command to build QuTiP in editable mode is ::

   python setup.py develop

from the repository directory.
If you now load up a Python interpreter, you should be able to ``import qutip`` from anywhere as long as the correct Python environment is active.
Any changes you make to the Python files in the git repository should be immediately present if you restart your Python interpreter and re-import ``qutip``.

On the first run, the setup command will compile many C++ extension modules built from Cython sources (files ending ``.pxd`` and ``.pyx``).
Generally the low-level linear algebra routines that QuTiP uses are written in these files, not in pure Python.
Unlike Python files, changes you make to Cython files will not appear until you run ``python setup.py develop`` again; you will only need to re-run this if you are changing Cython files.
Cython will detect and compile only the files that have been changed, so this command will be faster on subsequent runs.

.. note::

   When undertaking Cython development, the reason we use ``python setup.py develop`` instead of ``pip install -e .`` is because Cython's changed-file detection does not reliably work in the latter.
   ``pip`` tends to build in temporary virtual environments, which often makes Cython think its core library files have been updated, triggering a complete, slow rebuild of everything.

Code Style
----------

The biggest concern you should always have is to make it easy for your code to be read and understood by the person who comes next.

All new contributions must follow `PEP 8 style <https://www.python.org/dev/peps/pep-0008/>`_; all pull requests will be passed through a linter that will complain if you violate it.
You should use the ``pycodestyle`` package locally (available on ``pip``) to test you satisfy the requirements before you push your commits, since this is rather faster than pushing 10 different commits trying to fix minor niggles.
Keep in mind that there is quite a lot of freedom in this style, especially when it comes to line breaks.
If a line is too long, consider the *best* way to split it up with the aim of making the code readable, not just the first thing that doesn't generate a warning.

Try to stay consistent with the style of the surrounding code.
This includes using the same variable names, especially if they are function arguments, even if these "break" PEP 8 guidelines.
*Do not* change existing parameter, attribute or method names to "match" PEP 8; these are breaking user-facing changes, and cannot be made except in a new major release of QuTiP.

Other than this, general "good-practice" Python standards apply: try not to duplicate code; try to keep functions short, descriptively-named and side-effect free; provide a docstring for every new function; and so on.

Documenting
-----------

When you make changes in the core library, you should update the relevant documentation if needed.
If you are making a bug fix, or other relatively minor changes, you will probably only need to make sure that the docstrings of the modified functions and classes are up-to-date; changes here will propagate through to the documentation the next time it is built.
Be sure to follow the |numpydoc|_ when writing docstrings.
All docstrings will be parsed as reStructuredText, and will form the API documentation section of the documentation.

.. |numpydoc| replace:: Numpy documentation standards (``numpydoc``)
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html

Testing
-------

We use ``pytest`` as our test runner.
The base way to run every test is ::

   pytest /path/to/repo/qutip/tests

This will take around 10 to 30 minutes, depending on your computer and how many of the optional requirements you have installed.
It is normal for some tests to be marked as "skip" or "xfail" in yellow; these are not problems.
True failures will appear in red and be called "fail" or "error".

While prototyping and making changes, you might want to use some of the filtering features of ``pytest``.
Instead of passing the whole ``tests`` directory to the ``pytest`` command, you can also pass a list of files.
You can also use the ``-k`` selector to only run tests whose names include a particular pattern, for example ::

   pytest qutip/tests/test_qobj.py -k "expm"

to run the tests of :meth:`Qobj.expm`.


.. _contributing-docs:

Documentation: qutip/qutip (doc directory)
==========================================

The core library is in the `qutip/qutip repository on GitHub, inside the doc directory <https://github.com/qutip/qutip>`_.

Building
--------

The documentation is built using ``sphinx``, ``matplotlib`` and ``numpydoc``, with several additional extensions including ``sphinx-gallery`` and ``sphinx-rtd-theme``.
The most up-to-date instructions and dependencies will be in the ``README.md`` file of the documentation directory.
You can see the rendered version of this file simply by going to the `documentation GitHub page <https://github.com/qutip/qutip/tree/master/doc>`_ and scrolling down.

Building the documentation can be a little finnicky on occasion.
You likely will want to keep a separate Python environment to build the documentation in, because some of the dependencies can have tight requirements that may conflict with your favourite tools for Python development.
We recommend creating an empty ``conda`` environment containing only Python with ::

   conda create -n qutip-doc python=3.8

and install all further dependencies with ``pip``.
There is a ``requirements.txt`` file in the repository root that fixes all package versions exactly into a known-good configuration for a completely empty environment, using ::

   pip install -r requirements.txt

This known-good configuration was intended for Python 3.8, though in principle it is possible that other Python versions will work.

.. note::

   We recommend you use ``pip`` to install dependencies for the documentation rather than ``conda`` because several necessary packages can be slower to update their ``conda`` recipes, so suitable versions may not be available.

The documentation build includes running many components of the main QuTiP library to generate figures and to test the output, and to generate all the API documentation.
You therefore need to have a version of QuTiP available in the same Python environment.
If you are only interested in updating the users' guide, you can use a release version of QuTiP, for example by running ``pip install qutip``.
If you are also modifying the main library, you need to make your development version accessible in this environment.
See the `above section on building QuTiP <contributing-qutip_>`_ for more details, though the ``requirements.txt`` file will have already installed all the build requirements, so you should be able to simply run ::

   python setup.py develop

in the main library repository.

The documentation is built by running the ``make`` command.
There are several targets to build, but the most useful will be ``html`` to build the webpage documentation, ``latexpdf`` to build the PDF documentation (you will also need a full ``pdflatex`` installation), and ``clean`` to remove all built files.
The most important command you will want to run is ::

   make html

You should re-run this any time you make changes, and it should only update files that have been changed.

.. important::
   The documentation build includes running almost all the optional features of QuTiP.
   If you get failure messages in red, make sure you have installed all of the optional dependencies for the main library.

The HTML files will be placed in the ``_build/html`` directory.
You can open the file ``_build/html/index.html`` in your web browser to check the output.

Code Style
----------

All user guide pages and docstrings are parsed by Sphinx using reStructuredText.
There is a general `Sphinx usage guide <https://www.sphinx-doc.org/en/master/usage/index.html>`_, which has a lot of information that can sometimes be a little tricky to follow.
It may be easier just to look at other ``.rst`` files already in the documentation to copy the different styles.

.. note::
   reStructuredText is a very different language to the Markdown that you might be familiar with.
   It's always worth checking your work in a web browser to make sure it's appeared the way you intended.

Testing
-------

There are unfortunately no automated tests for the documentation.
You should ensure that no errors appeared in red when you ran ``make html``.
Try not to introduce any new warnings during the build process.
The main test is to open the HTML pages you have built (open ``_build/html/index.html`` in your web browser), and click through to the relevant pages to make sure everything has rendered the way you expected it to.
.. This file was created using retext 6.1 https://github.com/retext-project/retext

.. _release_distribution:

************************
Release and Distribution
************************

Preamble
++++++++

This document covers the process for managing updates to the current minor release and making new releases.
Within this document, the git remote ``upstream`` refers to the main QuTiP organsiation repository, and ``origin`` refers to your personal fork.

In short, the steps you need to take are:

1. Prepare the release branch (see git_).
2. Run the "Build wheels, optionally deploy to PyPI" GitHub action to build binary and source packages and upload them to PyPI (see deploy_).
3. Retrieve the built documentation from GitHub (see docbuild_).
4. Create a GitHub release and uploaded the built files to it (see github_).
5. Update `qutip.org <https://qutip.org/>`_ with the new links and documentation (web_).
6. Update the conda feedstock, deploying the package to ``conda`` (cforge_).



.. _git:

Setting Up The Release Branch
+++++++++++++++++++++++++++++

In this step you will prepare a git branch on the main QuTiP repository that has the state of the code that is going to be released.
This procedure is quite different if you are releasing a new minor or major version compared to if you are making a bugfix patch release.
For a new minor or major version, do update-changelog_ and then jump to release_.
For a bug fix to an existing release, do update-changelog_ and then jump to bugfix_.

Changes that are not backwards-compatible may only be made in a major release.
New features that do not affect backwards-compatibility can be made in a minor release.
Bug fix releases should be small, only fix bugs, and not introduce any new features.

There are a few steps that *should* have been kept up-to-date during day-to-day development, but might not be quite accurate.
For every change that is going to be part of your release, make sure that:

- The user guide in the documentation is updated with any new features, or changes to existing features.
- Any new API classes or functions have entries in a suitable RST file in ``doc/apidoc``.
- Any new or changed docstrings are up-to-date and render correctly in the API documentation.

Please make a normal PR to ``master`` correcting anything missing from these points and have it merged before you begin the release, if necessary.

.. _update-changelog:

Updating the Changelog
----------------------

This needs to be done no matter what type of release is being made.

#. Create a new branch to use to make a pull request.
#. Write the changelog for this version in ``doc/changelog.rst``.
   Look at recent entries in that file to get a feel for the style.
   In general, the format is one or two paragraphs written in regular prose describing the major new features of the version, and anything that needs special attention.
   After that, in suitable headings, list all the changes and who made them.
   Headings you may want to have include "Features", "Improvements", "Bug Fixes", "Deprecations", "Removals" and "Developer Changes", but feel free to use anything sensible.
#. Make a pull request on the main ``qutip/qutip`` repository with this changelog, and get other members of the admin team to approve it.
#. Merge this into ``master``.

Now jump to release_ if you are making a major or minor release, or bugfix_ if you are only fixing bugs in a previous release.

.. _release:

Create a New Minor or Major Release
-----------------------------------

This involves making a new branch to hold the release and adding some commits to set the code into "release" mode.
This release should be done by branching directly off the ``master`` branch at its current head.

#. On your machine, make sure your copy of ``master`` is up-to-date (``git checkout master; git pull upstream master``).
   This should at least involve fetching the changelog PR that you just made.
   Now create a new branch off a commit in ``master`` that has the state of the code you want to release.
   The command is ``git checkout -b qutip-<major>.<minor>.X``, for example ``qutip-4.7.X``.
   This branch name will be public, and must follow this format.
#. Push the new branch (with no commits in it relative to ``master``) to the main ``qutip/qutip`` repository (``git push upstream qutip-4.7.X``).
   Creating a branch is one of the only situations in which it is ok to push to ``qutip/qutip`` without making a pull request.
#. Create a second new branch, which will be pushed to your fork and used to make a pull request against the ``qutip-<major>.<minor>.X`` branch on ``qutip/qutip`` you just created.
   You can call this branch whatever you like because it is not going to the main repository, for example ``git checkout -b prepare-qutip-4.7.0``.
#. - Change the ``VERSION`` file to contain the new version number exactly, removing the ``.dev`` suffix.
     For example, if you are releasing the first release of the minor 4.7 track, set ``VERSION`` to contain the string ``4.7.0``.
     (*Special circumstances*: if you are making an alpha, beta or release candidate release, append a ``.a<n>``, ``.b<n>`` or ``.rc<n>`` to the version string, where ``<n>`` is an integer starting from 0 that counts how many of that pre-release track there have been.)
   - Edit ``setup.cfg`` by changing the "Development Status" line in the ``classifiers`` section to ::

        Development Status :: 5 - Production/Stable

   Commit both changes (``git add VERSION setup.cfg; git commit -m "Set release mode for 4.7.0"``), and then push them to your fork (``git push -u origin prepare-qutip-4.7.0``)
#. Using GitHub, make a pull request to the release branch (e.g. ``qutip-4.7.X``) using this branch that you just created.
   You will need to change the "base branch" in the pull request, because GitHub will always try to make the PR against ``master`` at first.
   When the tests have passed, merge this in.
#. Finally, back on ``master``, make a new pull request that changes the ``VERSION`` file to be ``<next-expected-version>.dev``, for example ``4.8.0.dev``.
   The "Development Status" in ``setup.cfg`` on ``master`` should not have changed, and should be ::

       Development Status :: 2 - Pre-Alpha

   because ``master`` is never directly released.

You should now have a branch that you can see on the GitHub website that is called ``qutip-4.7.X`` (or whatever minor version), and the state of the code in it should be exactly what you want to release as the new minor release.
If you notice you have made a mistake, you can make additional pull requests to the release branch to fix it.
``master`` should look pretty similar, except the ``VERSION`` will be higher and have a ``.dev`` suffix, and the "Development Status" in ``setup.cfg`` will be different.

You are now ready to actually perform the release.
Go to deploy_. 



.. _bugfix:

Create a Bug Fix Release
------------------------

In this you will modify an already-released branch by "cherry-picking" one or more pull requests that have been merged to ``master`` (including your new changelog), and bump the "patch" part of the version number.

#. On your machine, make sure your copy of ``master`` is up-to-date (``git checkout master; git pull upstream master``).
   In particular, make sure the changelog you wrote in the first step is visible.
#. Find the branch of the release that you will be modifying.
   This should already exist on the ``qutip/qutip`` repository, and be called ``qutip-<major>.<minor>.X`` (e.g. ``qutip-4.6.X``).
   If you cannot see it, run ``git fetch upstream`` to update all the branch references from the main repository.
   Checkout a new private branch, starting from the head of the release branch (``git checkout -b prepare-qutip-4.6.1 upstream/qutip-4.6.X``).
   You can call this branch whatever you like (in the example it is ``prepare-qutip-4.6.1``), because it will only be used to make a pull request.
#. Cherry-pick all the commits that will be added to this release in order, including your PR that wrote the new changelog entries (this will be the last one you cherry-pick).
   You will want to use ``git log`` to find the relevant commits, going from **oldest to newest** (their "age" is when they were merged into ``master``, not when the PR was first opened).
   The command is slightly different depending on which merge strategy was used for a particular PR:

   - "merge": you only need to find one commit though the log will have included several; there will be an entry in ``git log`` with a title such as "Merge pull request #1000 from <...>".
     Note the first 7 characters of its hash.
     Cherry-pick this by ``git cherry-pick --mainline 1 <hash>``.
   - "squash and merge": there will only be a single commit for the entire PR.
     Its name will be "<Name of the pull request> (#1000)".
     Note the first 7 characters of its hash.
     Cherry-pick this by ``git cherry-pick <hash>``.
   - "rebase and merge": this is the most difficult, because there will be many commits that you will have to find manually, and cherry-pick all of them.
     Go to the GitHub page for this PR, and go to the "Commits" tab.
     Using your local ``git log`` (you may find ``git log --oneline`` useful), find the hash for every single commit that is listed on the GitHub page, in order from **oldest to newest** (top-to-bottom in the GitHub view, which is bottom-to-top in ``git log``).
     You will need to use the commit message to do this; the hashes that GitHub reports will probably not be the same as how they appear locally.
     Find the first 7 characters of each of the hashes.
     Cherry-pick these all in one go by ``git cherry-pick <hash1> <hash2> ... <hash10>``, where ``<hash1>`` is the oldest.

   If any of the cherry-picks have merge conflicts, first verify that you are cherry-picking in order from oldest to newest.
   If you still have merge conflicts, you will either need to manually fix them (if it is a *very* simple fix), or else you will need to find which additional PR this patch depends on, and restart the bug fix process including this additional patch.
   This generally should not happen if you are sticking to very small bug fixes; if the fixes had far-reaching changes, a new minor release may be more appropriate.
#. Change the ``VERSION`` file by bumping the last number up by one (double-digit numbers are fine, so ``4.6.10`` comes after ``4.6.9``), and commit the change.
#. Push this branch to your fork, and make a pull request against the release branch.
   On GitHub in the PR screen, you will need to change the "Base" branch to ``qutip-4.6.X`` (or whatever version), because GitHub will default to making it against ``master``.
   It should be quite clear if you have forgotten to do this, because there will probably be many merge conflicts.
   Once the tests have passed and you have another admin's approval, merge the PR.

You should now see that the ``qutip-4.6.X`` (or whatever) branch on GitHub has been updated, and now includes all the changes you have just made.
If you have made a mistake, feel free to make additonal PRs to rectify the situation.

You are now ready to actually perform the release.
Go to deploy_. 


.. _deploy:

Build Release Distribution and Deploy
+++++++++++++++++++++++++++++++++++++

This step builds the source (sdist) and binary (wheel) distributions, and uploads them to PyPI (pip).
You will also be able to download the built files yourself in order to upload them to the QuTiP website.

Build and Deploy
----------------

This is handled entirely by a GitHub Action.
Go to the `"Actions" tab at the top of the QuTiP code repository <https://github.com/qutip/qutip/actions>`_.
Click on the "Build wheels, optionally deploy to PyPI" action in the left-hand sidebar.
Click the "Run workflow" dropdown in the header notification; it should look like the image below.

.. image:: /figures/release_guide_run_build_workflow.png

- Use the drop-down menu to choose the branch or tag you want to release from.
  This should be called ``qutip-4.5.X`` or similar, depending on what you made earlier.
  This must *never* be ``master``.
- To make the release to PyPI, type the branch name (e.g. ``qutip-4.5.X``) into the "Confirm chosen branch name [...]" field.
  You *may* leave this field blank to skip the deployment and only build the package.
- (Special circumstances) If for some reason you need to override the version number (for example if the previous deployment to PyPI only partially succeeded), you can type a valid Python version identifier into the "Override version number" field.
  You probably do not need to do this.
  The mechanism is designed to make alpha-testing major upgrades with nightly releases easier.
  For even a bugfix release, you should commit the change to the ``VERSION`` file.
- Click the lower "Run workflow" to perform the build and deployment.

At this point, the deployment will take care of itself.
It should take between 30 minutes and an hour, after which the new version will be available for install by ``pip install qutip``.
You should see the new version appear on `QuTiP's PyPI page <https://pypi.org/project/qutip>`_.

Download Built Files
--------------------

When the build is complete, click into its summary screen.
This is the main screen used to both monitor the build and see its output, and should look like the below image on a success.

.. image:: /figures/release_guide_after_workflow.png

The built binary wheels and the source distribution are the "build artifacts" at the bottom.
You need to download both the wheels and the source distribution.
Save them on your computer, and unzip both files; you should have many wheel ``qutip-*.whl`` files, and two sdist files: ``qutip-*.tar.gz`` and ``qutip-*.zip``.
These are the same files that have just been uploaded to PyPI.


Monitoring Progress (optional)
------------------------------

While the build is in progress, you can monitor its progress by clicking on its entry in the list below the "Run workflow" button.
You should see several subjobs, like the completed screen, except they might not yet be completed.

The "Verify PyPI deployment confirmation" should get ticked, no matter what.
If it fails, you have forgotten to choose the correct branch in the drop-down menu or you made a typo when confirming the correct branch, and you will need to restart this step.
You can check that the deployment instruction has been understood by clicking the "Verify PyPI deployment confirmation" job, and opening the "Compare confirmation to current reference" subjob.
You will see a message saying "Built wheels will be deployed" if you typed in the confirmation, or "Only building wheels" if you did not.
If you see "Only building wheels" but you meant to deploy the release to PyPI, you can cancel the workflow and re-run it after typing the confirmation.


.. _docbuild:

Getting the Built Documentation
+++++++++++++++++++++++++++++++

The documentation will have been built automatically for you by a GitHub Action when you merged the final pull request into the release branch before building the wheels.
You do not need to re-release the documentation on either GitHub or the website if this is a patch release, unless there were changes within it.

Go to the "Actions" tab at the top of the ``qutip/qutip`` repository, and click the "Build HTML documentation" heading in the left column.
You should see a list of times this action has run; click the most recent one whose name is exactly "Build HTML documentation", with the release branch name next to it (e.g. ``qutip-4.6.X``).
Download the ``qutip_html_docs`` artifact to your local machine and unzip it somewhere safe.
These are all the HTML files for the built documentation; you should be able to open ``index.html`` in your own web browser and check that everything is working.


.. _github:

Making a Release on GitHub
++++++++++++++++++++++++++

This is all done through `the "Releases" section <https://github.com/qutip/qutip/releases>`_ of the ``qutip/qutip`` repository on GitHub.

- Click the "Draft a new release" button.
- Choose the correct branch for your release (e.g. ``qutip-4.5.X``) in the drop-down.
- For the tag name, use ``v<your-version>``, where the version matches the contents of the ``VERSION`` file.
  In other words, if you are releasing a micro version 4.5.3, use ``v4.5.3`` as the tag, or if you are releasing major version 5.0.0, use ``v5.0.0``.
- The title is "QuTiP <your-version>", e.g. "QuTiP 4.6.0".
- For the description, write a short (~two-line for a patch release) summary of the reason for this release, and note down any particular user-facing changes that need special attention.
  Underneath, put the changelog you wrote when you did the documentation release.
  Note that there may be some syntax differences between the ``.rst`` file of the changelog and the Markdown of this description field (for example, GitHub's markdown typically maintains hard-wrap linebreaks, which is probably not what you wanted).
- Drag-and-drop all the ``qutip-*.whl``, ``qutip-*.tar.gz`` and ``qutip-*.zip`` files you got after the build step into the assets box.
  You may need to unzip the files ``wheels.zip`` and ``sdist.zip`` to find them if you haven't already; **don't** upload those two zip files.

Click on the "Publish release" button to finalise.


.. _web:

Website
+++++++

This assumes that qutip.github.io has already been forked and familiarity with the website updating workflow.
The documentation need not be updated for every patch release.

Copying New Files
-----------------

You only need to copy in new documentation to the website repository.
Do not copy the ``.whl``, ``.tar.gz`` or ``.zip`` files into the git repository, because we can access the public links from the GitHub release stage, and this keeps the website ``.git`` folder a reasonable size.

For all releases move (no new docs) or copy (for new docs) the ``qutip-doc-<MAJOR>.<MINOR>.pdf`` into the folder ``downloads/<MAJOR>.<MINOR>.<MICRO>``.

The legacy html documentation should be in a subfolder like ::

    docs/<MAJOR>.<MINOR>
    
For a major or minor release the previous version documentation should be moved into this folder. 

The latest version HTML documentation should be the folder ::

    docs/latest
    
For any release which new documentation is included
- copy the contents ``qutip/doc/_build/html`` into this folder. **Note that the underscores at start of the subfolder names will need to be removed, otherwise Jekyll will ignore the folders**. There is a script in the ``docs`` folder for this. 
https://github.com/qutip/qutip.github.io/blob/master/docs/remove_leading_underscores.py


HTML File Updates
-----------------

- Edit ``download.html``

    * The 'Latest release' version and date should be updated.
    * The tar.gz and zip links need to have their micro release numbers updated in their filenames, labels and trackEvent javascript.
      These links should point to the "Source code" links that appeared when you made in the GitHub Releases section.
      They should look something like ``https://github.com/qutip/qutip/archive/refs/tags/v4.6.0.tar.gz``.
    * For a minor or major release links to the last micro release of the previous version will need to be moved (copied) to the 'Previous releases' section.

- Edit ``_includes/sidebar.html``

    * The 'Latest release' version should be updated. The gztar and zip file links will need the micro release number updating in the traceEvent and file name.
    * The link to the documentation folder and PDF file (if created) should be updated.

- Edit ``documentation.html``

    * The previous release tags should be moved (copied) to the 'Previous releases' section.

.. _cforge:

Conda Forge
+++++++++++

If not done previously then fork the `qutip-feedstock <https://github.com/conda-forge/qutip-feedstock_>`_.

Checkout a new branch on your fork, e.g. ::

    $ git checkout -b version-4.0.2

Find the sha256 checksum for the tarball that the GitHub web interface generated when you produced the release called "Source code".
This is *not* the sdist that you downloaded earlier, it's a new file that GitHub labels "Source code".
When you download it, though, it will have a name that *looks* like it's the sdist ::

    $ openssl sha256 qutip-4.0.2.tar.gz

Edit the ``recipe/meta.yaml`` file.
Change the version at the top of the file, and update the sha256 checksum.
Check that the recipe package version requirements at least match those in ``setup.cfg``, and that any changes to the build process are reflected in ``meta.yml``.
Also ensure that the build number is reset ::

    build:
        number: 0

Push changes to your fork, e.g. ::

    $ git push --set-upstream origin version-4.0.2

Make a Pull Request.
This will trigger tests of the package build process.

If (when) the tests pass, the PR can be merged, which will trigger the upload of the packages to the conda-forge channel.
To test the packages, add the conda-forge channel with lowest priority ::

    $ conda config --append channels conda-forge

This should mean that the prerequistes come from the default channel, but the qutip packages are found in conda-forge.
.. _development_roadmap:

*************************
QuTiP Development Roadmap
*************************

Preamble
========

This document outlines plan and ideas for the current and future development of
QuTiP. The document is maintained by the QuTiP Admim team. Contributuions from
the QuTiP Community are very welcome.

In particular this document outlines plans for the next major release of qutip,
which will be version 5. And also plans and dreams beyond the next major
version.

There is lots of development going on in QuTiP that is not recorded in here.
This a just an attempt at coordinated stragetgy and ideas for the future.

.. _what-is-qutip:

What is QuTiP?
--------------

The name QuTiP refers to a few things. Most famously, qutip is a Python library
for simulating quantum dynamics. To support this, the library also contains
various software tools (functions and classes) that have more generic
applications, such as linear algebra components and visualisation utilities, and
also tools that are specifically quantum related, but have applications beyond
just solving dynamics (for instance partial trace computation).

QuTiP is also an organisation, in the Github sense, and in the sense of a group
of people working collaboratively towards common objectives, and also a web
presence `qutip.org <https://qutip.org/>`_. The QuTiP Community includes all the
people who have supported the project since in conception in 2010, including
manager, funders, developers, maintainers and users.

These related, and overlapping, uses of the QuTiP name are of little consequence
until one starts to consider how to organise all the software packages that are
somehow related to QuTiP, and specifically those that are maintained by the
QuTiP Admim Team. Herin QuTiP will refer to the project / organisation and qutip
to the library for simulating quantum dyanmics.

Should we be starting again from scratch, then we would probably chose another
name for the main qutip library, such as qutip-quantdyn. However, qutip is
famous, and the name will stay.


Library package structure
=========================

With a name as general as Quantum Toolkit in Python, the scope for new code
modules to be added to qutip is very wide. The library was becoming increasingly
difficult to maintain, and in c. 2020 the QuTiP Admim Team decided to limit the
scope of the 'main' (for want of a better name) qutip package. This scope is
restricted to components for the simulation (solving) of the dynamics of quantum
systems. The scope includes utilities to support this, including analysis and
visualisation of output.

At the same time, again with the intention of easing maintence, a decision to
limit dependences was agreed upon. Main qutip runtime code components should
depend only upon Numpy and Scipy. Installation (from source) requires Cython,
and some optional components also require Cython at runtime. Unit testing
requires Pytest. Visualisation (optional) components require Matplotlib.

Due to the all encompassing nature of the plan to abstract the linear algebra
data layer, this enhancement (developed as part of a GSoC project) was allowed
the freedom (potential for non-backward compatibility) of requiring a major
release. The timing of such allows for a restructuring of the qutip compoments,
such that some that could be deemed out of scope could be packaged in a
different way -- that is, not installed as part of the main qutip package. Hence
the proposal for different types of package described next. With reference to
the :ref:`discussion above <what-is-qutip>` on the name QuTiP/qutip, the planned
restructuring suffers from confusing naming, which seems unavoidable without
remaining either the organisation or the main package (neither of which are
desirable).

QuTiP family packages
  The main qutip package already has sub-packages,
  which are maintained in the main qutip repo. Any packages maitained by the
  QuTiP organisation will be called QuTiP 'family' packages. Sub-packages within
  qutip main will be called 'integrated' sub-packages. Some packages will be
  maintained in their own repos and installed separately within the main qutip
  folder structure to provide backwards compatibility, these are (will be)
  called qutip optional sub-packages. Others will be installed in their own
  folders, but (most likely) have qutip as a dependency -- these will just be
  called 'family' packages.

QuTiP affilliated packages
  Other packages have been developed by others
  outside of the QuTiP organisation that work with, and are complementary to,
  qutip. The plan is to give some recognition to those that we deem worthy of
  such [this needs clarification]. These packages will not be maintained by the
  QuTiP Team.

.. todo::

   Do we really need optional subpackages? It seems that are a bit fiddly and as
   we are side-stepping bw compat with a major version, we could just make QIP a
   separate package.

Family packages
---------------

.. _qmain:

qutip main
^^^^^^^^^^

* **current package status**: family package `qutip`
* **planned package status**: family package `qutip`

The in-scope components of the main qutip package all currently reside in the
base folder. The plan is to move some components into integrated subpackages as
follows:

- `core` quantum objects and operations
- `solver` quantum dynamics solvers

What will remain in the base folder will be miscellaneous modules. There may be
some opportunity for grouping some into a `visualisation` subpackage. There is
also some potential for renaming, as some module names have underscores, which
is unconventional.

Qtrl
^^^^

* **current package status**: integrated sub-package `qutip.control`
* **planned package status**: family package `qtrl`

There are many OSS Python packages for quantum control optimisation. There are
also many different algorithms. The current `control` integrated subpackage
provides the GRAPE and CRAB algorithms. It is too ambitious for QuTiP to attempt
(or want) to provide for all options. Control optimisation has been deemed out
of scope and hence these components will be separated out into a family package
called Qtrl.

Potentially Qtrl may be replaced by separate packages for GRAPE and CRAB, based
on the QuTiP Control Framework.

QIP
^^^

* **current package status**: integrated sub-package `qutip.qip`
* **planned package status**: optional sub-package `qutip.qip`

.. todo::

   Is it really necessary for this to be a sub-package? It could just be a
   separate package.

The QIP subpackage has been deemed out of scope (feature-wise). It also depends
on `qutip.control` and hence would be out of scope for dependency reasons. A
separate repository has already been made for qutip-qip.

qutip-symbolic
^^^^^^^^^^^^^^

* **current package status**: independent package `sympsi`
* **planned package status**: family package `qutip-symbolic`

Long ago Robert Johansson and Eunjong Kim developed Sympsi. It is a fairly
coomplete library for quantum computer algebra (symbolic computation). It is
primarily a quantum wrapper for `Sympy <https://www.sympy.org>`_.

It has fallen into unmaintained status. The latest version on the `sympsi repo
<https://github.com/sympsi/sympsi>`_ does not work with recent versions of
Sympy. Alex Pitchford has a `fork <https://github.com/ajgpitch/sympsi>`_ that
does 'work' with recent Sympy versions -- unit tests pass, and most examples
work. However, some (important) examples fail, due to lack of respect for
non-commuting operators in Sympy simplifcation functions (note this was true as
of Nov 2019, may be fixed now).

There is a [not discussed with RJ & EK] plan to move this into the QuTiP family
to allow the Admin Team to maintain, develop and promote it. The 'Sympsi' name
is cute, but a little abstract, and qutip-symbolic is proposed as an
alternative, as it is plainer and more distinct from Sympy.


Affilliated packages
--------------------

qucontrol-krotov
^^^^^^^^^^^^^^^^

* **code repository**: https://github.com/qucontrol/krotov

A package for quantum control optimisation using Krotov, developed mainly by
Michael Goerz.

Generally accepted by the Admin Team as well developed and maintained. A solid
candiate for affilliation.


Development Projects
====================

.. _dl-abs:

data layer abstraction
----------------------

:tag: dl-abs
:status: majority of development completed.
:admin lead: `Eric <https://github.com/Ericgig>`_
:main dev: `Jake Lishman <https://github.com/jakelishman>`_

Development completed as a GSoC project. Fully implemented in the dev.major
branch. Currently being used by some research groups.

Abstraction of the linear algebra data from code qutip components, allowing
for alternatives, such as sparse, dense etc. Difficult to summarize. Almost
every file in qutip affected in some way. A major milestone for qutip.
Significant performance improvements throughout qutip.

Some developments tasks remain, including providing full control over how the
data-layer dispatchers choose the most appropriate output type.

.. _qmain-reorg:

qutip main reorganization
-------------------------

:tag: qmain-reorg
:status: development [pretty much] complete
:admin lead: `Eric <https://github.com/Ericgig>`_
:main dev: `Jake Lishman <https://github.com/jakelishman>`_

Reorganise qutip main components to the structure :ref:`described above <qmain>`.

.. _qmain-docs:

qutip user docs migration
-------------------------

:tag: qmain-docs
:status: conceptualised
:admin lead: TBA; `Shahnawaz <https://github.com/quantshah>`_?
:main dev: TBA

The qutip user documentation build files are to be moved to the qutip/qutip
repo. This is more typical for an OSS package.

As part of the move, the plan is to reconstruct the Sphinx structure from
scratch. Historically, there have been many issues with building the docs.
Sphinx has come a long way since qutip docs first developed. The main source
(rst) files will remain [pretty much] as they are, although there is a lot of
scope to improve them.

The qutip-doc repo will afterwards just be used for documents, such as this one,
pertaining to the QuTiP project.

.. _solve-dl:

Solver data layer integration
-----------------------------

:tag: solve-dl
:status: development ongoing
:admin lead: `Eric <https://github.com/Ericgig>`_
:main dev: `Eric <https://github.com/Ericgig>`_

The new data layer gives opportunity for significantly improving performance of
the qutip solvers. Eric has been revamping the solvers by deploying `QobjEvo`
(the time-dependent quantum object) that he developed. `QobjEvo` will exploit
the data layer, and the solvers in turn exploit `QobjEvo`.

.. _qip-mig:

QIP migration
-------------

:tag: qip-mig
:status: development [pretty much] complete
:admin lead: `Boxi <https://github.com/BoxiLi>`_
:main dev: `Sidhant Saraogi <https://github.com/sarsid>`_

A separate package for qutip-qip was created during Sidhant's GSoC project.
There is some fine tuning required, especially after qutip.control is migrated.

.. _heom-revamp:

HEOM revamp
-----------

:tag: heom-revamp
:status: development [pretty much] complete
:admin lead: `Neill <https://github.com/nwlambert>`_
:main dev: `Tarun Raheja <https://github.com/tehruhn>`_

An overhaul of the HEOM solver. C++ components used to speed up construction of
the hierarchy.

.. _qtrl-mig:

Qtrl migration
--------------

:tag: qtrl-mig
:status: conceptualised
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

The components currently packaged as an integrated subpackage of qutip main will
be moved to separate package called Qtrl. This is the original codename of the
package before it was integrated into qutip. Also changes to exploit the new
data layer will be implemented.

.. _ctrl-fw:

QuTiP control framework
-----------------------

:tag: ctrl-fw
:status: conceptualised
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

Create new package qutip-ctrlfw "QuTiP Control Framework". The aim is provide a
common framework that can be adopted by control optimisation packages, such that
different packages (algorithms) can be applied to the same problem.

Classes for defining a controlled system:

- named control parameters. Scalar and n-dim. Continuous and discrete variables
- mapping of control parameters to dynamics generator args
- masking for control parameters to be optimised

Classes for time-dependent variable parameterisation

- piecewise constant
- piecewise linear
- Fourier basis
- more

Classes for defining an optimisation problem:

- single and multiple objectives

.. _qutip-optim:

QuTiP optimisation
------------------

:tag: qutip-optim
:status: conceptualised
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

A wrapper for multi-variable optimisation functions. For instance those in
`scipy.optimize` (Nelder-Mead, BFGS), but also others, such as Bayesian
optimisation and other machine learning based approaches. Initially just
providing a common interface for quantum control optimisation, but applicable
more generally.

.. _sympsi-mig:

Sympsi migration
----------------

:tag: sympsi-mig
:status: conceptualised
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

Create a new family package qutip-symbolic from ajgpitch fork of Sympy. Must
gain permission from Robert Johansson and Eunjong Kim. Extended Sympy simplify
to respect non-commuting operators. Produce user documentation.

.. _status-mig:

Status messaging and recording
------------------------------

:tag: status-msg
:status: conceptualised
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

QuTiP has various ways of recording and reporting status and progress.

- `ProgressBar` used by some solvers
- Python logging used in qutip.control
- `Dump` used in qutip.control
- heom records `solver.Stats`

Some consolidation of these would be good.

Some processes (some solvers, correlation, control optimisation) have many
stages and many layers. `Dump` was initially developed to help with debugging,
but it is also useful for recording data for analysis. qutip.logging_utils has
been criticised for the way it uses Python logging. The output goes to stderr
and hence the output looks like errors in Jupyter notebooks.

Clearly, storing process stage data is costly in terms of memory and cpu time,
so any implementation must be able to be optionally switched on/off, and avoided
completely in low-level processes (cythonized components).

Required features:

- optional recording (storing) of process stage data (states, operators etc)
- optionally write subsets to stdout
- maybe other graphical representations
- option to save subsets to file
- should ideally replace use of `ProgressBar`, Python logging, `control.Dump`, `solver.Stats`

.. _qutip-gui:

qutip Interactive
-----------------

:status: conceptualised
:tag: qutip-gui
:admin lead: `Alex <https://github.com/ajgpitch>`_
:main dev: TBA

QuTiP is pretty simple to use at an entry level for anyone with basic Python
skills. However, *some* Python skills are necessary. A graphical user interface
(GUI) for some parts of qutip could help make qutip more accessible. This could
be particularly helpful in education, for teachers and learners.

This would make an good GSoC project. It is independent and the scope is
flexible.

The scope for this is broad and flexible. Ideas including, but not limited to:

Interactive Bloch sphere
^^^^^^^^^^^^^^^^^^^^^^^^

Matplotlib has some interactive features (sliders, radio buttons, cmd buttons)
that can be used to control parameters. They are a bit clunky to use, but they
are there. Could maybe avoid these and develop our own GUI. An interactive Bloch
sphere could have sliders for qubit state angles. Buttons to add states, toggle
state evolution path.

Interactive solvers
^^^^^^^^^^^^^^^^^^^

Options to configure dynamics generators (Lindbladian / Hamiltonian args etc)
and expectation operators. Then run solver and view state evolution.

Animated circuits
^^^^^^^^^^^^^^^^^

QIP circuits could be animated. Status lights showing evolution of states during
the processing. Animated Bloch spheres for qubits.


QuTiP major release roadmap
===========================

QuTiP v.5
---------

These Projects need to be completed for the qutip v.5 release.

- :ref:`dl-abs`
- :ref:`qmain-reorg`
- :ref:`qmain-docs`
- :ref:`solve-dl`
- :ref:`qip-mig`
- :ref:`qtrl-mig`
- :ref:`heom-revamp`

The planned timeline for the release is:

- **alpha version, April 2021**. Core features packaged and available for
  experienced users to test.
- **beta version, July 2021**. All required features and documentation complete,
  packaged and ready for community testing.
- **full release, September 2021**. Full tested version released.
.. _user_guide.rst:

************************************
Working with the QuTiP Documentation
************************************


The user guide provides an overview of QuTiP's functionality.
The guide is composed of individual reStructuredText (``.rst``) files which each get rendered as a webpage.
Each page typically tackles one area of functionality.
To learn more about how to write ``.rst`` files, it is useful to follow the `sphinx guide <https://www.sphinx-doc.org/en/master/usage/index.html>`_.

The documentation build also utilizes a number of
`Sphinx Extensions <https://www.sphinx-doc.org/en/master/usage/extensions/index.html>`_
including but not limited to
`doctest <https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html>`_,
`autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_,
`sphinx gallery <https://sphinx-gallery.github.io/stable/index.html>`_ and
`plot <https://matplotlib.org/3.1.1/devel/plot_directive.html>`_.
Additional extensions can be configured in the `conf.py <https://github.com/qutip/qutip/blob/master/doc/conf.py>`_ file.

.. _directives.rst:

Directives
==========

There are two Sphinx directives that can be used to write code examples in the user guide:

- `Doctest <https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html>`_
- `Plot <https://matplotlib.org/3.1.1/devel/plot_directive.html>`_

For a more comprehensive account of the usage of each directive, please refer to their individual pages. Here we outline some general guidelines on how to these directives while making a user guide.

Doctest
-------

The doctest directive enables tests on interactive code examples.
The simplest way to do this is by specifying a prompt along with its respective output: ::

    .. doctest::

        >>> a = 2
        >>> a
        2

This is rendered in the documentation as follows:

.. doctest::

    >>> a = 2
    >>> a
    2


While specifying code examples under the ``.. doctest::`` directive, either all statements must be specified by the ``>>>`` prompt or without it.
For every prompt, any potential corresponding output must be specified immediately after it.
This directive is ideally used when there are a number of examples that need to be checked in quick succession.

A different way to specify code examples (and test them) is using the associated ``.. testcode::`` directive which is effectively a code block: ::

    .. testcode::

        a = 2
        print(a)

followed by its results.
The result can be specified with the ``.. testoutput::`` block: ::

    .. testoutput::

        2

The advantage of the ``testcode`` directive is that it is a lot simpler to
specify and amenable to copying the code to clipboard. Usually, tests are
more easily specified with this directive as the input and output are
specified in different blocks. The rendering is neater too.

.. note::
    The ``doctest`` and ``testcode`` directives should not be assumed to
    have the same namespace.

**Output:**

.. testcode::

    a = 2
    print(a)

.. testoutput::

    2

A few notes on using the doctest extension:

- By default, each ``testcode`` and ``doctest`` block is run in a fresh namespace.
  To share a common namespace, we can specify a common group across the blocks
  (within a single ``.rst`` file). For example, ::

        .. doctest:: [group_name]

          >>> a = 2

  can be followed by some explanation followed by another code block
  sharing the same namespace ::

        .. doctest:: [group_name]

          >>> print(a)
          2

- To only print the code blocks (or the output), use the option ``+SKIP`` to
  specify the block without the code being tested when running ``make doctest``.

- To check the result of a ``Qobj`` output, it is useful to make sure that
  spacing irregularities between the expected and actual output are ignored.
  For that, we can use the option ``+NORMALIZE_WHITESPACE``.

Plot
----

Since the doctest directive cannot render matplotlib figures, we use Matplotlib's
`Plot <https://matplotlib.org/3.1.1/devel/plot_directive.html>`_
directive when rendering to LaTeX or HTML.

The plot directive can also be used in the doctest format. In this case,
when running doctests (which is enabled by specifying all statements with the
``>>>`` prompts), tests also include those specified under the plot directive.

**Example:**
::

    First we specify some data:

    .. plot::

      >>> import numpy as np
      >>> x = np.linspace(0, 2 * np.pi, 1000)
      >>> x[:10] # doctest: +NORMALIZE_WHITESPACE
      array([ 0.        ,  0.00628947,  0.01257895,  0.01886842,  0.0251579 ,
              0.03144737,  0.03773685,  0.04402632,  0.0503158 ,  0.05660527])


    .. plot::
      :context:

      >>> import matplotlib.pyplot as plt
      >>> plt.plot(x, np.sin(x))
      [...]

Note the use of the ``NORMALIZE_WHITESPACE`` option to ensure that the
multiline output matches.

**Render:**


.. plot::

    >>> import numpy as np
    >>> x = np.linspace(0, 2 * np.pi, 1000)
    >>> x[:10] # doctest: +SKIP
    array([ 0.        ,  0.00628947,  0.01257895,  0.01886842,  0.0251579 ,
            0.03144737,  0.03773685,  0.04402632,  0.0503158 ,  0.05660527])
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, np.sin(x))
    [...]

A few notes on using the plot directive:

- A useful argument to specify in plot blocks is that of ``context`` which ensures
  that the code is being run in the namespace of the previous plot block within the
  same file.

- By default, each rendered figure in one plot block (when using ``:context:``)
  is carried over to the next block.

- When the ``context`` argument is specified with the ``reset`` option
  as ``:context: reset``, the namespace is reset to a new one and all figures are
  erased.

- When the ``context`` argument is specified with the ``close-figs`` option
  as ``:context: reset``, the namespace is reset to a new one and all figures are
  erased.


The Plot directive cannot be used in conjunction with Doctest because they do not
share the same namespace when used in the same file.
Since Plot can also be used in doctest mode, in
the case where code examples require both testing and rendering figures, it is
easier to use the Plot directive. To learn more about each directive, it is useful
to refer to their individual pages.
.. _development_ideas:

**********************************
Ideas for future QuTiP development
**********************************

This chapter covers the development of QuTiP and its subpackages, including
a roadmap for upcoming releases and ideas for future improvements.

.. toctree::
   :maxdepth: 1

   ideas/qutip-interactive.rst
   ideas/pulse-level-quantum-circuits.rst
   ideas/tensorflow-data-backend.rst
   ideas/quantum-error-mitigation.rst
   ideas/heom-gpu.rst
***********************
TensorFlow Data Backend
***********************

.. contents:: Contents
    :local:
    :depth: 3

QuTiP's data layer provides the mathematical operations needed to work with
quantum states and operators, i.e. ``Qobj``, inside QuTiP. As part of Google
Summer of Code 2020, the data layer was rewritten to allow new backends to
be added more easily and for different backends to interoperate with each
other. Backends using in-memory spares and dense matrices already exist,
and we would like to add a backend that implements the necessary operations
using TensorFlow [1]_.

Why a TensorFlow backend?
-------------------------

TensorFlow supports distributing matrix operations across multiple GPUs and
multiple machines, and abstracts away some of the complexities of doing so
efficiently. We hope that by using TensorFlow we might enable QuTiP to scale
to bigger quantum systems (e.g. more qubits) and decrease the time taken to
simulate them.

There is particular interest in trying the new backend with the
BoFiN HEOM (Hierarchical Equations of Motion) solver [2]_.

Challenges
----------

TensorFlow is a very different kind of computational framework to the existing
dense and sparse matrix backends. It uses flow graphs to describe operations,
and to work efficiently. Ideally large graphs of operations need to be
executed together in order to efficiently compute results.

The QuTiP data layer might need to be adjusted to accommodate these
differences, and it is possible that this will prove challenging or even
that we will not find a reasonable way to achieve the desired performance.

Expected outcomes
=================

* Add a ``qutip.core.data.tensorflow`` data type.
* Implement specialisations for some important operations (e.g. ``add``,
  ``mul``, ``matmul``, ``eigen``, etc).
* Write a small benchmark to show how ``Qobj`` operations scale on the new
  backend in comparison to the existing backends. Run the benchmark both
  with and without using a GPU.
* Implement enough for a solver to run on top of the new TensorFlow data
  backend and benchmark that (stretch goal).

Skills
======

* Git, Python and familiarity with the Python scientific computing stack
* Familiarity with TensorFlow (beneficial, but not required)
* Familiarity with Cython (beneficial, but not required)

Difficulty
==========

* Medium

Mentors
=======

* Simon Cross (hodgestar@gmail.com)
* Jake Lishman (jake@binhbar.com)
* Alex Pitchford (alex.pitchford@gmail.com)

References
==========

.. [1] https://www.tensorflow.org/
.. [2] https://github.com/tehruhn/bofin
**********************************************************
GPU implementation of the Hierarchical Equations of Motion
**********************************************************

.. contents:: Contents
    :local:
    :depth: 3

The Hierarchical Equations of Motion (HEOM) method is a non-perturbative
approach to simulate the evolution of the density matrix of dissipative quantum
systems. The underlying equations are a system of coupled ODEs which can be run
on a GPU. This will allow the study of larger systems as discussed in [1]_. The
goal of this project would be to extend QuTiP's HEOM method [2]_ and implement
it on a GPU.

Since the method is related to simulating large, coupled ODEs, it can also be
quite general and extended to other solvers.

Expected outcomes
=================

* A version of HEOM which runs on a GPU.
* Performance comparison with the CPU version.
* Implement dynamic scaling.

Skills
======

* Git, python and familiarity with the Python scientific computing stack
* CUDA and OpenCL knowledge

Difficulty
==========

* Hard

Mentors
=======

* Neill Lambert (nwlambert@gmail.com)
* Alex Pitchford (alex.pitchford@gmail.com)
* Shahnawaz Ahmed (shahnawaz.ahmed95@gmail.com)
* Simon Cross (hodgestar@gmail.com)

References
==========

.. [1] https://pubs.acs.org/doi/abs/10.1021/ct200126d?src=recsys&journalCode=jctcce
.. [2] https://arxiv.org/abs/2010.10806
*****************
QuTiP Interactive
*****************

.. contents:: Contents
    :local:
    :depth: 3

QuTiP is pretty simple to use at an entry level for anyone with basic Python
skills. However, *some* Python skills are necessary. A graphical user interface
(GUI) for some parts of qutip could help make qutip more accessible. This could
be particularly helpful in education, for teachers and learners.

Ideally, interactive components could be embedded in web pages. Including, but
not limited to, Jupyter notebooks.

The scope for this is broad and flexible. Ideas including, but not limited to:

Interactive Bloch sphere
------------------------

QuTiP has a Bloch sphere virtualisation for qubit states. This could be made
interactive through sliders, radio buttons, cmd buttons etc. An interactive
Bloch sphere could have sliders for qubit state angles. Buttons to add states,
toggle state evolution path. Potential for recording animations. Matplotlib has
some interactive features (sliders, radio buttons, cmd buttons) that can be used
to control parameters. that could potentially be used.

Interactive solvers
-------------------

Options to configure dynamics generators (Lindbladian / Hamiltonian args etc)
and expectation operators. Then run solver and view state evolution.

Animated circuits
-----------------

QIP circuits could be animated. Status lights showing evolution of states during
the processing. Animated Bloch spheres for qubits.

Expected outcomes
=================

* Interactive graphical components for demonstrating quantum dynamics
* Web pages for qutip.org or Jupyter notebooks introducing quantum dynamics
  using the new components

Skills
======

* Git, Python and familiarity with the Python scientific computing stack
* elementary understanding of quantum dynamics

Difficulty
==========

* Variable

Mentors
=======

* Nathan Shammah (nathan.shammah@gmail.com)
* Alex Pitchford (alex.pitchford@gmail.com)
* Simon Cross (hodgestar@gmail.com)
* Boxi Li (etamin1201@gmail.com) [QuTiP GSoC 2019 graduate]
************************
Quantum Error Mitigation
************************

.. contents:: Contents
    :local:
    :depth: 3

From the QuTiP 4.5 release, the qutip.qip module now contains the noisy quantum
circuit simulator (which was a GSoC project) providing enhanced features for a
pulse-level description of quantum circuits and noise models. A new class
`Processor` and several subclasses are added to represent different platforms
for quantum computing. They can transfer a quantum circuit into the
corresponding control sequence and simulate the dynamics with QuTiP solvers.
Different noise models can be added to `qutip.qip.noise` to simulate noise in a
quantum device.

This module is still young and many features can be improved, including new
device models, new noise models and integration with the existing general
framework for quantum circuits (`qutip.qip.circuit`). There are also possible
applications such as error mitigation techniques ([1]_, [2]_, [3]_).

The tutorial notebooks can be found at https://qutip.org/tutorials.html#nisq. A
recent presentation on the FOSDEM conference may help you get an overview
(https://fosdem.org/2020/schedule/event/quantum_qutip/). See also the Github
Project page for a collection of related issues and ongoing Pull Requests.

Expected outcomes
=================

- Make an overview of existing libraries and features in error mitigation,
  similarly to a literature survey for a research article, but for a code
  project (starting from Refs. [4]_, [5]_). This is done in order to best
  integrate the features in QuTiP with existing libraries and avoid
  reinventing the wheel.
- Features to perform error mitigation techniques in QuTiP, such as zero-noise
  extrapolation by pulse stretching.
- Tutorials implementing basic quantum error mitigation protocols
- Possible integration with Mitiq [6]_

Skills
======

* Background in quantum physics and quantum circuits.
* Git, python and familiarity with the Python scientific computing stack

Difficulty
==========

* Medium

Mentors
=======

* Nathan Shammah (nathan.shammah@gmail.com)
* Alex Pitchford (alex.pitchford@gmail.com)
* Eric Giguère (eric.giguere@usherbrooke.ca)
* Neill Lambert (nwlambert@gmail.com)
* Boxi Li (etamin1201@gmail.com) [QuTiP GSoC 2019 graduate]

References
==========

.. [1] Kristan Temme, Sergey Bravyi, Jay M. Gambetta, **Error mitigation for short-depth quantum circuits**, Phys. Rev. Lett. 119, 180509 (2017)

.. [2] Abhinav Kandala, Kristan Temme, Antonio D. Corcoles, Antonio Mezzacapo, Jerry M. Chow, Jay M. Gambetta,
 **Extending the computational reach of a noisy superconducting quantum processor**, Nature *567*, 491 (2019)

.. [3] S. Endo, S.C. Benjamin, Y. Li, **Practical quantum error mitigation for near-future applications**, Physical Review X *8*, 031027 (2018)

.. [4] Boxi Li's blog on the GSoC 2019 project on pulse-level control, https://gsoc2019-boxili.blogspot.com/

.. [5] Video of a recent talk on the GSoC 2019 project, https://fosdem.org/2020/schedule/event/quantum_qutip/

.. [6] `Mitiq <https://mitiq.readthedocs.io/>`_
*******************************************
Pulse level description of quantum circuits
*******************************************

.. contents:: Contents
    :local:
    :depth: 3

The aim of this proposal is to enhance QuTiP quantum-circuit compilation
features with regard to quantum information processing. While QuTiP core modules
deal with dynamics simulation, there is also a module for quantum circuits
simulation. The two subsequent Google Summer of Code projects, in 2019 and 2020,
enhanced them in capabilities and features, allowing the simulation both at the
level of gates and at the level of time evolution. To connect them, a compiler
is implemented to compile quantum gates into the Hamiltonian model. We would
like to further enhance this feature in QuTiP and the connection with other
libraries.

Expected outcomes
=================

* APIs to import and export pulses to other libraries. Quantum compiler is a
  current research topic in quantum engineering. Although QuTiP has a simple
  compiler, many may want to try their own compiler which is more compatible
  with their quantum device. Allowing importation and exportation of control
  pulses will make this much easier. This will include a study of existing
  libraries, such as `qiskit.pulse` and `OpenPulse` [1]_, comparing them with
  `qutip.qip.pulse` module and building a more general and comprehensive
  description of the pulse.

* More examples of quantum system in the `qutip.qip.device` module. The circuit
  simulation and compilation depend strongly on the physical system. At the
  moment, we have two models: spin chain and cavity QED. We would like to
  include some other commonly used planform such as Superconducting system [2]_,
  Ion trap system [3]_ or silicon system. Each model will need a new set of
  control Hamiltonian and a compiler that finds the control pulse of a quantum
  gate. More involved noise models can also be added based on the physical
  system. This part is going to involve some physics and study of commonly used
  hardware platforms. The related code can be found in `qutip.qip.device` and
  `qutip.qip.compiler`.

Skills
======

* Git, Python and familiarity with the Python scientific computing stack
* quantum information processing and quantum computing (quantum circuit formalism)

Difficulty
==========

* Medium

Mentors
=======

* Boxi Li (etamin1201@gmail.com) [QuTiP GSoC 2019 graduate]
* Nathan Shammah (nathan.shammah@gmail.com)
* Alex Pitchford (alex.pitchford@gmail.com)

References
==========

.. [1] McKay D C, Alexander T, Bello L, et al. Qiskit backend specifications for openqasm and openpulse experiments[J]. arXiv preprint arXiv:1809.03452, 2018.

.. [2] Häffner H, Roos C F, Blatt R, **Quantum computing with trapped ions**, Physics reports, 2008, 469(4): 155-203.

.. [3] Krantz P, Kjaergaard M, Yan F, et al. **A quantum engineer's guide to superconducting qubits**, Applied Physics Reviews, 2019, 6(2): 021318.
