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
