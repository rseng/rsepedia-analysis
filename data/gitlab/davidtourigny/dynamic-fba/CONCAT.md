# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at  USAdst2156@cumc.columbia.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

---
title: 'dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python'
tags:
 - Mathematical optimization
 - Linear programming
 - Numerical simulation
 - Multi-scale metabolic modeling
 - C++
 - Python
authors:
 - name: David S. Tourigny
   orcid: 0000-0002-3987-8078
   affiliation: 1
 - name: Jorge Carrasco Muriel
   orcid: 0000-0001-7365-0299
   affiliation: 2
 - name: Moritz E. Beber
   orcid: 0000-0003-2406-1978
   affiliation: 2
affiliations:
 - name: Columbia University Irving Medical Center, 630 West 168th Street, New York, NY 10032 USA
   index: 1
 - name: Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark, Building 220, Kemitorvet, 2800 Kongens Lyngby, Denmark
   index: 2
date: 13 March 2020
bibliography: paper.bib
---

# Summary

Flux-balance analysis (FBA) is a computational method based on linear programming (LP) that has had enormous success modeling the metabolic behaviour of organisms and cellular systems existing in steady state with their environment [@Varma94; @Orth10]. Representing the corresponding biological model as an LP problem means that FBA can be used to study metabolism at genome-scale. Unfortunately, the underlying assumption of an unchanging environment means that FBA is not immediately applicable to systems with dynamics where, for example, environmental conditions may vary in time. Extensions of FBA, including dynamic FBA (DFBA) [@Mahadevan02], have therefore been developed in order to accommodate temporal dynamics into the framework of genome-scale metabolic modeling.

Although DFBA is well-defined mathematically as an LP problem embedded in a system of ordinary differential equations (ODEs), numerical simulation of DFBA models proves particularly challenging (as described in @Harwood16). Consequently, @Harwood16 proposed an algorithm for efficiently simulating DFBA models based on reformulating the ODE system as a differential algebraic equation (DAE) with root detection and representing solutions of the LP problem using an optimal basis formulation. An initial implementation of this algorithm has been provided in the software package DFBAlab [@Gomez14] written in MATLAB and using commercial LP solvers.

Increasingly, researchers engaged in metabolic modeling prefer open source software. Python is quickly becoming their platform of choice [@Carey20]. Among other packages, open source resources for building and simulating FBA models using Python can be found in the COBRApy [@Ebrahim13] package which is part of the openCOBRA organization [@opencobra]. Until now, COBRApy lacked an efficient implementation of DFBA using the DAE formulation.

## Statement of need:
_Researchers wanting to build and simulate specific models of interest often lack the background in numerical analysis or high-performance computing required to overcome the numerical challenges of DFBA_.

We have solved this issue by developing a software package based on open source libraries GLPK [@glpk] and SUNDIALS [@Hindmarsh05] that implements the most efficient algorithms in a compiled programming language that is made accessible to users through a simple and intuitive pybind11 [@pybind11] interface with pandas [@McKinney11] and the openCOBRA Python module COBRApy [@Ebrahim13]. ODEs are constructed using the symbolic expression enabled by optlang [@Jensen16], SymPy, and SymEngine [@Meurer17].

# Acknowledgements

DST is a Simons Foundation Fellow of the Life Sciences Research Foundation. MEB
received funding from the European Union’s Horizon 2020 research and innovation
programme under grant agreement 686070 (DD-DeCaF). We thank Peter St. John and Christian Diener for discussions and suggestions.

# References
## *EMBLP* source files

This directory contains the following content

* [`./emblp_direct.cpp`](./emblp_direct.cpp): contains member functions for
  direct method derived class
* [`./emblp_harwood.cpp`](./emblp_harwood.cpp): contains member functions for
  Harwood et al. derived class
* [`./emblp.cpp`](./emblp.cpp): contains member functions for embedded
  LP problem base class
* [`./emblp_scott.h`](./emblp_scott.h): contains member functions for Scott et
  al. derived class
* [`./emblp.h`](./emblp.h): contains class declarations for embedded LP problems
* [`./README.md`](./README.md): this document

## *METHODS* source files

This directory contains the following content

* [`./methods_direct.cpp`](./methods_direct.cpp): contains source code for model integration
  using direct method
* [`./methods_harwood.cpp`](./methods_harwood.cpp): contains source code for model integration
  using Harwood et al. algorithm
* [`./methods.h`](./methods.h): contains *SUNDIALS* includes and declarations
  for integration methods
* [`./methods.cpp`](./methods.cpp): contains functions used by methods
* [`./methods_scott.cpp`](./methods_scott.cpp): contains source code for model integration using
  Scott et al. algorithm
* [`./README.md`](./README.md): this document

#### Problem description

Please explain:
* **what** you tried to achieve,
* **how** you went about it (referring to the code sample), and
* **why** the current behaviour is a problem and what output
  you expected instead.

#### Code Sample

Create a [minimal, complete, verifiable example
](https://stackoverflow.com/help/mcve).

```python
# Paste your code here or link to a gist.
```

```
# If there was a crash, please include the traceback here.
```

### Context

Please run the following code and paste the output inside the details
block.

```
python -c "import dfba;dfba.show_versions()"
```

<details>

</details>
* [ ] fix #(issue number)
* [ ] description of feature/fix
* [ ] tests added/passed
* [ ] add an entry to the [next release](../../CHANGELOG.rst)
## Examples

This directory contains the following content

* [`./example1.py`](./example1.py): anaerobic growth of *E. coli* on glucose and
  xylose
* [`./example2.py`](./example2.py): aerobic growth of *E. coli* on glucose and
  xylose
* [`./example3.py`](./example3.py): aerobic growth of *S. cerevisiae* on glucose
  with switch to anaerobic conditions at *t=7.7h*
* [`./example4.py`](./example4.py): aerobic growth of *S. cerevisiae* on glucose and
  xylose
* [`./example5.py`](./example5.py): anaerobic growth of *E. coli* on glucose and
  xylose simulated using direct method
* [`./README.md`](./README.md): this document

## Build scripts

This directory contains the following content

* [`./build_pybind11.sh`](./build_pybind11.sh): script for downloading, building and
  installing specified pybind11
* [`./build_glpk.sh`](./build_glpk.sh): script for downloading, building and
  installing specified GLPK version
* [`./build_sundials.sh`](./build_sundials.sh): script for downloading, building
  and installing specified SUNDIALS version
* [`./jupyterlab_plotly.sh`](./jupyterlab_plotly.sh): script for building jupyterlab extensions
* [`./README.md`](./README.md): this document

## cmake modules

* [`./FindGLPK.cmake`](./FindGLPK.cmake) adapted from [ycollet/coinor-cmake](https://github.com/ycollet/coinor-cmake/blob/master/Clp/cmake/FindGlpk.cmake)
* [`./FindSUNDIALS.cmake`](./FindSUNDIALS.cmake) adapted from [casadi/casadi](https://github.com/casadi/casadi/blob/master/cmake/FindSUNDIALS.cmake)
* [`./README.md`](./README.md): this document

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Report Bugs
===========

Report bugs at https://gitlab.com/davidtourigny/dynamic-fba/issues.

If you are reporting a bug, please follow the template guide lines. The more 
detailed your report, the easier and thus faster we can help you.

Fix Bugs
========

Look through the GitLab issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
==================

Look through the GitLab issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
===================

dynamic-fba could always use more documentation, whether as part of the
official documentation, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
===============

The best way to send feedback is to file an issue at
https://gitlab.com/davidtourigny/dynamic-fba/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
============

Ready to contribute? Here's how to set up dynamic-fba for
local development.

1. Fork the https://gitlab.com/davidtourigny/dynamic-fba
   repository on GitLab.
2. Clone your fork locally

   .. code-block:: console
   
       git clone git@gitlab.com:your_name_here/dynamic-fba.git

3. Install your local copy into a a Python virtual environment.
   You can `read this guide to learn more
   <https://realpython.com/python-virtual-environments-a-primer/>`_
   about them and how to create one. Alternatively, particularly if you are a 
   Windows or Mac user, you can also use
   `Anaconda <https://docs.anaconda.com/anaconda/>`_. Assuming you have 
   virtualenvwrapper installed, this is how you set up your fork for local development

   .. code-block:: console
   
       mkvirtualenv my-env
       cd dynamic-fba/
       pip install -e .[development]

4. Create a branch for local development using the ``devel`` branch as a 
   starting point. Use ``fix`` or ``feat`` as a prefix

   .. code-block:: console
   
       git checkout devel
       git checkout -b fix-name-of-your-bugfix

   Now you can make your changes locally.

5. When you're done making changes, apply the quality assurance tools and check 
   that your changes pass our test suite. This is all included with tox

   .. code-block:: console
   
       make qa
       tox

   You can run all tests in parallel using detox. To get detox, just
   pip install it into your virtualenv.

6. Commit your changes and push your branch to GitLab. Please use `semantic
   commit messages <http://karma-runner.github.io/2.0/dev/git-commit-msg.html>`_.

   .. code-block:: console
   
       git add .
       git commit -m "fix: Your summary of changes"
       git push origin fix-name-of-your-bugfix

7. Open the link displayed in the message when pushing your new branch 
   in order to submit a pull request.

Pull Request Guidelines
=======================

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
3. The pull request should work for Python 3.6 and 3.7. This is also ensured 
   by our Travis CI.
=======
History
=======

Next Release
------------
* Upcoming features and fixes

0.1.0 (2019-07-19)
------------------
* First release of the dfba package with five examples applying dynamic FBA.
=======
Support
=======

* dynamic-fba `gitter chat <https://gitter.im/opencobra/dynamic-fba>`_

=============================
Dynamic Flux Balance Analysis
=============================

.. image:: https://img.shields.io/pypi/v/dfba.svg
   :target: https://pypi.org/project/dfba/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/dfba.svg
   :target: https://pypi.org/project/dfba/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/dfba.svg
   :target: http://www.gnu.org/licenses/
   :alt: GPLv3+

.. image:: https://gitlab.com/davidtourigny/dynamic-fba/badges/master/pipeline.svg
   :target: https://travis-ci.org/davidtourigny/dynamic-fba/commits/master
   :alt: Pipeline Status

.. image:: https://gitlab.com/davidtourigny/dynamic-fba/badges/master/coverage.svg
   :target: https://gitlab.com/davidtourigny/dynamic-fba/commits/master
   :alt: Coverage Report

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. image:: https://joss.theoj.org/papers/10.21105/joss.02342/status.svg
   :target: https://doi.org/10.21105/joss.02342

.. _`Harwood et al., 2016`: https://link.springer.com/article/10.1007/s00211-015-0760-3
.. _GLPK: https://www.gnu.org/software/glpk/
.. _SUNDIALS: https://computation.llnl.gov/projects/sundials
.. _Python: https://www.python.org/
.. _cobrapy: https://github.com/opencobra/cobrapy
.. _optlang: https://github.com/biosustain/optlang
.. _symengine: https://github.com/symengine/symengine

This project provides an object-oriented software package for dynamic
flux-balance analysis (DFBA) simulations using implementations of the direct
method or Algorithm 1 described in the paper `Harwood et al., 2016`_. The main
algorithms for solving embedded LP problems are written in *C++* and use the GNU
Linear Programming Kit (GLPK_) and the Suite of Nonlinear and
Differential/Algebraic Equation Solvers (SUNDIALS_) CVODE or IDA. Extension
modules to cobrapy_ are provided for easy generation and simulation of DFBA
models.

Installation
============

.. _GLPK: https://www.gnu.org/software/glpk/
.. _SUNDIALS: https://computation.llnl.gov/projects/sundials
.. _Python: https://www.python.org/
.. _cobrapy: https://github.com/opencobra/cobrapy
.. _optlang: https://github.com/biosustain/optlang
.. _symengine: https://github.com/symengine/symengine

Currently, we do not provide Python wheels for this package and therefore `Installing from source`_ is a bit more involved. However, the package is now available on `conda-forge <https://anaconda.org/conda-forge/dfba>`_ and can be installed in a conda environment. First, create an environment with Python_ (e.g. version 3.7):

.. code-block:: console

    conda create --name dfba python=3.7
    

Next, activate the environment using

.. code-block:: console

    conda activate dfba
    
and install `dfba <https://anaconda.org/conda-forge/dfba>`_ using

.. code-block:: console

     conda install -c conda-forge dfba
     
You can then use the software within the environment on the examples in the repository as well as your own DFBA models.


The quickest way to run the software without building anything locally
is from the provided `Docker <https://docs.docker.com/>`_ image:

.. code-block:: console

    docker run --rm -it davidtourigny/dfba:latest


Installing from source
----------------------

Currently this package is compatible with most UNIX-like operating systems.
Provided the following `Dependencies`_ are installed, the module
can be installed from source using the command:

.. code-block:: console

    pip install dfba

Dependencies
~~~~~~~~~~~~

.. _`build_glpk.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/master/scripts/build_glpk.sh
.. _`build_pybind11.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/master/scripts/build_pybind11.sh
.. _`build_sundials.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/master/scripts/build_sundials.sh
.. _Dockerfile: https://gitlab.com/davidtourigny/dynamic-fba/tree/master/Dockerfile
.. _`pybind11`: https://github.com/pybind/pybind11


* A version of Python_ 3.6 or higher is required
* You need `cmake <https://cmake.org/>`_ for the build process
* You will need `git <https://git-scm.com/>`_ to clone this repository to access
  the scripts and build files
* You need a working compiler with C++11 support, for example, by installing
  ``build-essential`` on Debian-derived Linux systems
* GLPK_ version 4.65 is required or can be installed using `build_glpk.sh`_
* SUNDIALS_ version 5.0.0 or higher is required or can be installed using `build_sundials.sh`_
* pybind11_ is required or can be installed using `build_pybind11.sh`_

Be aware that some of these packages have their own dependencies that must
therefore be installed also (e.g. GLPK_ depends on `GMP <https://gmplib.org/>`_
and pybind11_ requires `pytest <https://docs.pytest.org/en/latest/>`_).


Alternatively, a Dockerfile_ is provided for building a `Docker <https://docs.docker.com/>`_
image to run the software from an interactive container. The `Docker <https://docs.docker.com/>`_ image can be
built in one step by issuing the command:

.. code-block:: console

    make build

from the root of this repository. It can then be started using:

.. code-block:: console

    make run

Documentation
=============

Documentation for dfba is provided at `readthedocs <https://dynamic-fba.readthedocs.io>`_

Authors
=======

* David S. Tourigny
* Moritz E. Beber

Additional contributors
=======================

* Jorge Carrasco Muriel (visualization and documentation)

Citation
========

Tourigny DS, Carrasco Muriel J, Beber ME (2020). dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python. `Journal of Open Source Software, 5(52), 2342 <https://doi.org/10.21105/joss.02342>`_

Copyright
=========

* Copyright © 2018,2019 Columbia University Irving Medical Center, New York, USA
* Copyright © 2019 Novo Nordisk Foundation Center for Biosustainability,
  Technical University of Denmark
* Free software distributed under the `GNU General Public License v3 or later
  (GPLv3+) <http://www.gnu.org/licenses/>`_.
============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Report Bugs
===========

Report bugs at https://gitlab.com/davidtourigny/dynamic-fba/issues.

If you are reporting a bug, please follow the template guide lines. The more 
detailed your report, the easier and thus faster we can help you.

Fix Bugs
========

Look through the GitLab issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
==================

Look through the GitLab issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
===================

dfba could always use more documentation, whether as part of the
official documentation, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
===============

The best way to send feedback is to file an issue at
https://gitlab.com/davidtourigny/dynamic-fba/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
============

Ready to contribute? Here's how to set up dfba for
local development.

1. Fork the https://gitlab.com/davidtourigny/dynamic-fba
   repository on GitLab.
2. Clone your fork locally

   .. code-block:: console
   
       git clone git@gitlab.com:your_name_here/dynamic-fba.git

3. Install your local copy into a a Python virtual environment.
   You can `read this guide to learn more
   <https://realpython.com/python-virtual-environments-a-primer/>`_
   about them and how to create one. Alternatively, particularly if you are a 
   Windows or Mac user, you can also use
   `Anaconda <https://docs.anaconda.com/anaconda/>`_. Assuming you have 
   virtualenvwrapper installed, this is how you set up your fork for local development

   .. code-block:: console
   
       mkvirtualenv my-env
       cd dynamic-fba/
       pip install -e .[development]

4. Create a branch for local development using the ``devel`` branch as a 
   starting point. Use ``fix`` or ``feat`` as a prefix

   .. code-block:: console
   
       git checkout devel
       git checkout -b fix-name-of-your-bugfix

   Now you can make your changes locally.

5. When you're done making changes, apply the quality assurance tools and check 
   that your changes pass our test suite. This is all included with tox

   .. code-block:: console
   
       make qa
       tox

   You can run all tests in parallel using detox. To get detox, just
   pip install it into your virtualenv.

6. Commit your changes and push your branch to GitLab. Please use `semantic
   commit messages <http://karma-runner.github.io/2.0/dev/git-commit-msg.html>`_.

   .. code-block:: console
   
       git add .
       git commit -m "fix: Your summary of changes"
       git push origin fix-name-of-your-bugfix

7. Open the link displayed in the message when pushing your new branch 
   in order to submit a pull request.

Pull Request Guidelines
=======================

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
3. The pull request should work for Python 3.6 and 3.7. This is also ensured 
   by our Travis CI.
======
Citing
======

If this package has contributed to your own work you can use the following citation:

Tourigny DS, Carrasco Muriel J, Beber ME (2020). dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python. `Journal of Open Source Software, 5(52), 2342 <https://doi.org/10.21105/joss.02342>`_
Example DFBA models
===================

.. _example1: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example1.py
.. _example2: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example2.py
.. _example3: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example3.py
.. _example4: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example4.py
.. _example5: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example5.py

The current version is distributed with several examples related to Examples
6.2.1 and 6.3 in `Harwood et al., 2016 <https://link.springer.com/article/10.1007/s00211-015-0760-3>`_.

Examples example1_ and example2_ are based on `Hanly & Henson, 2011 <https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.22954>`_ and also Example 1 in
`DFBAlab <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0409-8>`_.
Example example5_ implements the same model as example1_,
but uses the direct method in place of Algorithm 1 from
`Harwood et al., 2016 <https://link.springer.com/article/10.1007/s00211-015-0760-3>`_. Conversely, example3_ and example4_ are based on
`Hanly & Henson, 2011 <https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.22954>`_.

Genome-scale Metabolic Networks
-------------------------------

.. _sbml-models: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/sbml-models

The genome-scale metabolic models used in these examples are provided in the directory sbml-models_

* `iJR904 <http://bigg.ucsd.edu/models/iJR904>`_: *Escherichia coli* bacterium
  iJR904 contains 761 metabolites and 1075 reaction fluxes.
* `iND750 <http://bigg.ucsd.edu/models/iND750/>`_: *Saccharomyces cerevisiae*
  strain S288C iND750 contains 1059 metabolites and 1266 reaction fluxes.
Installation
============

.. _GLPK: https://www.gnu.org/software/glpk/
.. _SUNDIALS: https://computation.llnl.gov/projects/sundials
.. _Python: https://www.python.org/
.. _cobrapy: https://github.com/opencobra/cobrapy
.. _optlang: https://github.com/biosustain/optlang
.. _symengine: https://github.com/symengine/symengine
.. _pandas: https://pandas.pydata.org
.. _example1: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples/example1.py

Currently, we do not provide Python wheels for this package and therefore :ref:`Installing from source` is a bit more involved. However, the package is now available on `conda-forge <https://anaconda.org/conda-forge/dfba>`_ and can be installed in a conda environment. First, create an environment with Python_ (e.g. version 3.7):

.. code-block:: console

    conda create --name dfba python=3.7
    

Next, activate the environment using

.. code-block:: console

    conda activate dfba
    
and install `dfba <https://anaconda.org/conda-forge/dfba>`_ using

.. code-block:: console

     conda install -c conda-forge dfba
     
You can then use the software within the environment on the examples in the repository as well as your own DFBA models.


The quickest way to run the software without building anything locally
is from the provided `Docker <https://docs.docker.com/>`_ image. There are several ways to do this; we suggest using the following command to run Example 1 from the root of this repsitory:

.. code-block:: console

    docker run -it -v ${PWD}:/opt/examples davidtourigny/dfba python3 examples/example1.py

Provided `Docker <https://docs.docker.com/>`_ is installed correctly, the first time you run this command the latest image will be pulled from dockerhub (unless you have previously pulled or built it yourself using the instructions at the bottom of this page). All subsequent runs will use this same image, and you can replace  ``examples/example1.py`` with the path to a Python_ script containing your own DFBA model. To access results of each simulation, we recommend including the pandas_ method  ``to_csv()`` inside your script in order to write the data frame of results to a csv file (e.g., see :doc:`example1` or uncomment lines 107 and 108 in script example1_).

.. _Installing from source:

Installing from source
----------------------------------------

Currently this package is compatible with most UNIX-like operating systems.
Provided the following :ref:`Dependencies` are installed, the module
can be installed using the command:

.. code-block:: console

    pip install dfba

.. _Dependencies:

Dependencies
~~~~~~~~~~~~

.. _`build_glpk.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/scripts/build_glpk.sh
.. _`build_pybind11.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/scripts/build_pybind11.sh
.. _`build_sundials.sh`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/scripts/build_sundials.sh
.. _Dockerfile: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/Dockerfile
.. _`pybind11`: https://github.com/pybind/pybind11


* A version of Python_ 3.6 or higher is required
* You need `cmake <https://cmake.org/>`_ for the build process
* You will need `git <https://git-scm.com/>`_ to clone this repository to access
  the scripts and build files
* You need a working compiler with C++11 support, for example, by installing
  ``build-essential`` on Debian-derived Linux systems
* GLPK_ version 4.65 is required or can be installed using `build_glpk.sh`_
* SUNDIALS_ version 5.0.0 or higher is required or can be installed using `build_sundials.sh`_
* pybind11_ is required or can be installed using `build_pybind11.sh`_

Be aware that some of these packages have their own dependencies that must
therefore be installed also (e.g. GLPK_ depends on `GMP <https://gmplib.org/>`_
and pybind11_ requires `pytest <https://docs.pytest.org/en/latest/>`_).


Alternatively, a Dockerfile_ is provided for building a `Docker <https://docs.docker.com/>`_
image to run the software from an interactive container. The `Docker <https://docs.docker.com/>`_ image can be
built in one step by issuing the command:

.. code-block:: console

    make build

from the root of the `repository <https://gitlab.com/davidtourigny/dynamic-fba>`_. It can then be started using:

.. code-block:: console

    make run
.. dfba documentation main file, created by
   sphinx-quickstart on Thu Nov 28 10:05:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dfba's documentation!
=======================================

.. image:: https://img.shields.io/pypi/v/dfba.svg
   :target: https://pypi.org/project/dfba/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/dfba.svg
   :target: https://pypi.org/project/dfba/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/dfba.svg
   :target: http://www.gnu.org/licenses/
   :alt: GPLv3+

.. image:: https://gitlab.com/davidtourigny/dynamic-fba/badges/main/pipeline.svg
   :target: https://travis-ci.org/davidtourigny/dynamic-fba/commits/main
   :alt: Pipeline Status

.. image:: https://gitlab.com/davidtourigny/dynamic-fba/badges/main/coverage.svg
   :target: https://gitlab.com/davidtourigny/dynamic-fba/commits/main
   :alt: Coverage Report

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v1.4%20adopted-ff69b4.svg
   :target: https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
   :alt: Contributor Covenant

.. image:: https://joss.theoj.org/papers/10.21105/joss.02342/status.svg
   :target: https://doi.org/10.21105/joss.02342

.. _GLPK: https://www.gnu.org/software/glpk/
.. _SUNDIALS: https://computation.llnl.gov/projects/sundials
.. _Python: https://www.python.org/
.. _cobrapy: https://github.com/opencobra/cobrapy
.. _optlang: https://github.com/biosustain/optlang
.. _symengine: https://github.com/symengine/symengine

This project provides an object-oriented software package for dynamic
flux-balance analysis (DFBA) simulations using implementations of the direct
method or Algorithm 1 described in the paper `Harwood et al., 2016 <https://link.springer.com/article/10.1007/s00211-015-0760-3>`_. The main
algorithms for solving embedded LP problems are written in C++ and use the GNU
Linear Programming Kit (GLPK_) and the Suite of Nonlinear and
Differential/Algebraic Equation Solvers (SUNDIALS_) CVODE or IDA. Extension
modules to cobrapy_ are provided for easy generation and simulation of DFBA
models.

Overview
========

.. _`pybind11`: https://github.com/pybind/pybind11
.. _examples: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/examples

This software package implements the fastest and most robust DFBA simulation algorithms in a compiled programming language that is made accessible to users through a simple and intuitive interface with cobrapy_. The target audience is researchers wanting to build and efficiently simulate specific DFBA models of interest, but not spend their time wrestling with the numerical challenges of DFBA. Users are not expected to interact directly with the lower-level *C++* interface
and once installed the package should ideally remain untouched. Instead, the
classes and functions for solving embedded LP problems have been exposed to
Python_ using `pybind11`_. Combined with the provided cobrapy_ extension
modules, this provides the user with the ability to build their own DFBA model
exclusively in Python_.  

The Python_ class :class:`dfba.DfbaModel` intuitively encapsulates
all the data required for a full definition of a DFBA model by combining an
underlying cobrapy_ object with instances of the :class:`dfba.KineticVariable` and
:class:`dfba.ExchangeFlux` classes. The :class:`dfba.DfbaModel` class instance ensures all user data are
consistent with the initialization and simulation requirements of an embedded LP
problem. User data are passed directly to the algorithms and symbolic functions
are dynamically compiled and loaded prior to simulation.  

The directory examples_ also contains scripts for the examples described in the :doc:`examples` section,
and details on how the user can adapt these to build and simulate their own model are outlined in the section :doc:`example1`.

.. toctree::
   :maxdepth: 2
   :hidden:

   installation.rst
   plotting.rst
   examples.rst
   example1.ipynb
   example3.ipynb
   example5.ipynb
   contributing.rst
   citing.rst
   source_files.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Source Files
============

.. _src: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src
.. _extension: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension
.. _`dfba_utils.cpp`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension/dfba_utils.cpp
.. _emblp: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension/emblp
.. _methods: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension/methods
.. _`solver_data.h`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension/solver_data.h
.. _`user_data.h`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/extension/user_data.h
.. _dfba: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba
.. _`control.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/control.py
.. _`exchange.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/exchange.py
.. _`helpers.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/helpers.py
.. _`jit.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/jit.py
.. _`model.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/model.py
.. _`library.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/library.py
.. _`variables.py`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/variables.py
.. _`plot`: https://gitlab.com/davidtourigny/dynamic-fba/tree/main/src/dfba/plot


Source files contained within the directory src_ are split between two
sub-directories separated by their language of implementation.

*C++*
-----

.. _Python: https://www.python.org/

The sub-directory extension_ contains the following content:

* `dfba_utils.cpp`_: contains source code for exposing the extension module to Python_
* emblp_: contains class and function declarations for embedded LP problems
* methods_: contains algorithms for integration of embedded LP problems
* `solver_data.h`_: struct exposed to Python_ for solver options
* `user_data.h`_: struct exposed to Python_ for model specification

*Python*
--------

The directory dfba_ contains the following content:

* `control.py`_: definition of class :class:`dfba.ControlParameter`
* `exchange.py`_: definition of class :class:`dfba.ExchangeFlux`
* `helpers.py`_: general helper functions
* `jit.py`_: tools for JIT compilation of dynamic library
* `model.py`_: definition of class :class:`dfba.DfbaModel`
* `library.py`_: methods for writing dynamic library
* `variables.py`_: definition of class :class:`dfba.KineticVariable`
* `plot`_: directory for additional visualization dependency
Visualization
-------------

Visualization tools are available as an extra dependency, optionally installed 
from the root of the `repository <https://gitlab.com/davidtourigny/dynamic-fba>`_ using the commands.

.. code-block:: console

    pip install .[plotly]

or

.. code-block:: console

    pip install .[matplotlib]

Usage of these libraries in the package is exemplified in :doc:`example1`.

It is also possible to run `Jupyterlab <https://jupyterlab.readthedocs.io/en/stable/>`_
in an interactive `Docker <https://docs.docker.com/>`_ container and access the visualization tools
through your browser:

.. code-block:: console

    docker run --rm -it -p 8888:8888 -v $PWD:/opt davidtourigny/dfba:latest-interactive bash

Once the image is running, start `Jupyterlab <https://jupyterlab.readthedocs.io/en/stable/>`_ using the command:

.. code-block:: console

    jupyter lab --ip 0.0.0.0 --no-browser --allow-root

This will display a URL that should be used to open the notebook in your browser. From
there you will be able to interact with the examples and build your own DFBA model.
