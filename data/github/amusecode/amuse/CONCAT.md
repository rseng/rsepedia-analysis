# AMUSE: The Astrophysical Multipurpose Software Environment
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3466805.svg)](https://doi.org/10.5281/zenodo.3466805)
[![PyPI version](https://badge.fury.io/py/amuse.svg)](https://badge.fury.io/py/amuse)

This repository contains the AMUSE software. With AMUSE you can write
scripts to simulate astrophysical problems in different domains.

The project website is:

* www.amusecode.org

and the documentation can be found at:

* https://amuse.readthedocs.io

Getting Started
===============

In short, most probably

```bash
pip install amuse
```
should get you going if you have a linux or Mac were you compile 
codes on (HDF5 and an MPI libraries must be installed). 

Below are some hints for a quick install, if these fail please 
look for options at the detailed descriptions of the installation 
procedure in the documents in the 'doc/install' directory.

Compilers
=========

To build AMUSE from source you need to have a working build
environment. The AMUSE build system needs C/C++ and fortan 90
compilers, we recommend a recent version of GCC. 

In Ubuntu you can setup the environment with (as root):

```bash
apt-get install build-essential curl g++ gfortran gettext zlib1g-dev
```

Other distributions have similar package or package groups available.

In macOS you can use the homebrew or macports package manager (both
require the Apple Developer Tools and Xcode to be installed).

For a Windows 10 machine, AMUSE can be installed in the Windows Subsystem
for linux (WSL), and installing e.g. Ubuntu from the Microsoft store. 
Its recommended to use WSL 2. For further installation instructions, see the 
Linux install instructions.

Python
======

AMUSE needs Python 3 version >=3.5 installed preferably with pip and 
virtualenv. It may be necessary to update pip to a recent version.
If you cannot use Python 3, legacy support for Python 2 is available in the 
AMUSE 12 release and the python2 branch.

Installing Prerequisites
========================

The following libraries need to be installed:

* HDF (version 1.6.5 - 1.8.x)
* MPI (OpenMPI or MPICH)

The following are needed for some codes:
* FFTW (version >= 3.0)
* GSL
* CMake (version >= 2.4)
* GMP (version >= 4.2.1)
* MPFR (version >= 2.3.1)

Installing+building AMUSE
=========================

AMUSE can be installed through pip:

```bash
pip install [--user] amuse
```

This will build and install AMUSE with an extensive set of codes.
If necessary this will also install some required Python packages:

* Numpy (version >= 1.3.0)
* h5py (version >= 1.2.0)
* mpi4py (version >= 1.0)
* pytest (version >= 5.0)
* docutils (version >= 0.6)

If you are not using pip these must be installed by hand.

It is possible to install the minimal framework by:

```bash
pip install [--user] amuse-framework
```

This does not include any codes. These can be added
```bash
pip install [--user] amuse-<code name>
```

AMUSE Development 
=================

An AMUSE development install can also be handled through pip by executing (in the root of a clone of the 
repository)

```bash
pip install -e .
```

after this the codes need to be build:

```bash
python setup.py develop_build
```

Running the tests
=================
AMUSE comes with a large set of tests, most can be run automatically.
To run these tests start the py.test command from the main
amuse directory (directory this README file lives in).

To run these tests do:

1. install the tests

```bash
pip install [--user] amuse-tests
```
(this will install all tests whether or not you have installed the full amuse package)

2. Run the automatic tests

```bash
pytest --pyargs -v amuse.test.suite
```
you can also just run the tests for the specific packages you have installed e.g.
```bash
pytest --pyargs amuse.test.suite.codes_tests.test_huayno
```
you may have to prefix ```mpiexec -n 1 --oversubscribe``` to the pytest command.
Contributing to AMUSE
=====================

So you're interested in contributing code to AMUSE? Excellent! 

Reporting Issues
----------------

When opening an issue to report a problem, please try and provide a minimal
code example that reproduces the issue, and also include details of the
operating system, compiler, and the Python, Numpy, and AMUSE versions you are using.

Contributing code
-----------------

Most contributions to AMUSE are done via pull requests from GitHub users'
forks of the [amuse repository](https://github.com/amusecode/amuse).

Once you open a pull request (which should be opened against the ``master``
branch, not against any of the other branches), please make sure that you
include the following:

- **Code**: the code you are adding

- **Tests**: these are usually tests to ensure that code that previously
  failed now works (regression tests) or tests that cover as much as possible
  of the new functionality to make sure it doesn't break in future, and also
  returns consistent results on all platforms (since we run these tests on many
  platforms/configurations). 


Checklist for Contributed Code
------------------------------

A pull request for a new feature will be reviewed to see if it meets the
following requirements.  For any pull request, an AMUSE maintainer can
help to make sure that the pull request meets the requirements for inclusion
in the package.

**Scientific Quality**
(when applicable)
  * Is the submission relevant to AMUSE?
  * Are references included to the origin paper for the simulation code?
  * Does the code perform as expected?
  * Has the code been tested against previously existing codes in the same domain?

**Code Quality**
  * Is the code compatible with Python >=2.7?
  * Are there dependencies other than AMUSE, MPI, the Python Standard
    Library, and NumPy?
  * For compatibility reasons we prefer code that also works on older 
    versions of Numpy, matplotlib etc.
  * Are additional dependencies handled appropiately? If possible, factor out 
    additional dependencies or make them optional.
  * Does the code follow the AMUSE Style Guide (http://www.amusecode.org/doc/reference/style_guide.html)?

**Testing**
  * Are the inputs to the functions sufficiently tested?
  * Are there tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Does python setup.py test run without failures?

**Documentation**
  * Is there any information needed to be added to the docs to describe the code?

**License**
  * Does the code require a specific license other than the AMUSE license?
  * Are there any conflicts with this code and AMUSE?

**AMUSE requirements**
  * Can you checkout the pull request and repeat the examples and tests?
asterisk
========

Astrophysics visualization for AMUSE, based on eSight


Getting started:
================

Before running the code, please add the root directory of this project to the classpath.
### Bonsai2

This is a more optimized version of Bonsai1, however it also contains fewer features.
This version does not contain stopping conditions or individual/dynamic time-steps.

To run the code, make sure to type the following in your console otherwise the code hangs/crashes.

`$ ulimit -s unlimited`


AMUSE-Distributed
=================

Distributed code for AMUSE project
ETICS
=====

This is ETICS (Expansion Techniques In Collisionless Systems), a GPU (currently
**CUDA only**) N-body code which uses series expansion to calculate the
gravitational field. See more details in this publication:

Meiron, Y., Li, B., Holley-Bockelmann, K., & Spurzem, R. 2014, ApJ, 792, 98

http://adsabs.harvard.edu/abs/2014ApJ...792...98M


What's inside
-------------

- ETICS standalone program
- ETICS static library
- ETICS module for AMUSE


Prerequisites
-------------

- CUDA (>= 6; mandatory)
- HDF5 (optional)
- Boost (optional)
- AMUSE (mandatory only for AMUSE module)

To disable HDF5 [insert explanation here]
To disable Boost [insert explanation here]


Compilation
-----------

### Standalone program

    make standalone

builds in `src/`.


### Static library

    make library

builds in `src/`.


### AMUSE module

    make

builds in top level directory. The whole `etics` directory has to be placed (or
linked) inside:

    $AMUSE_DIR/src/amuse/community/


How to use
----------

The `file.ini` is a self-explanatory input file; if compiled without Boost, fill
in the relevant variables in file `noboost.inc` which is compiled into the
executable (any change of parameters requires re-compilation). Start simulation
with:

    ./etics file.ini

Any input file name is acceptable.


Known issues
------------

* No MEX

The MEX (Multipole Expansion) method is not available in this version; thus, the
SCF (Self-Consistent Field) method is the only expansion technique available.
The ETICS program has been heavily restructured and the MEX routines are no
longer compatible. Hopefully this will be fixed.

* Hardcoded launch configuration

For at least one of the CUDA kernels, for various reasons, it seems that "brute
force" search is needed to find the optimal launch configuration. Currently it
is hardcoded, and a primitive search routine is provided.

* Problem for particles with |θ| << 1

Due to using an unstable recursion relation to calculate the Associated Legendre
polynomials, particles in a narrow cone around the z-axis cannot be considered
accurately. This means that they give an erroneous contribution to the
gravitational field and also are assigned erroneous force and potential. The
size of this cone increases with the angular term of the expansion. To partly
solve this, the current (ugly) fix is to only consider particles with cos(θ) >
0.999 at the monopole level. This is not so bad because the monopole is always
the most dominant term (and is error free) and the number of particles in this
cone is small and they are transient (i.e. they come out of it usually in a
small number of steps). A somewhat better solution is to make this arbitrary
cutoff of 0.999 l-dependent, and an even better solution would be to use an
asymptotic expression for the Associated Legendre polynomials.
---
name: Other issue
about: Report an issue that is not a feature request or bug report
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:

**Expected behavior**
A clear and concise description of what you expected to happen.

**Logs**
If applicable, add logfiles to help explain your problem.

**Environment (please complete the following information):**
 - OS and version: [e.g. macOS High Sierra; Ubuntu Linux 17.10]
 - Compiler: [e.g. gcc6]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest a feature or improvement for AMUSE
title: ''
labels: feature request
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or logfiles about the feature request here.
---
name: Question
about: Ask an AMUSE-related question
title: ''
labels: question
assignees: ''

---


This package installs the Mameclot community code for AMUSE.
This package installs the Astrophysical Multipurpose Software Environment (AMUSE).
This package installs the Gadget2 community code for AMUSE.
This package installs the Mercury community code for AMUSE.
This package installs the ph4 community code for AMUSE.
This package installs the Capreole community code for AMUSE.
This package installs the Mikkola community code for AMUSE.
This package installs the FI community code for AMUSE.
This package installs the Athena community code for AMUSE.
This package installs the Hermite community code for AMUSE.
This package installs the EVTwin community code for AMUSE.
This package installs the SecularMultiple community code for AMUSE.
This package installs the Astrophysical Multipurpose Software Environment (AMUSE).
This package installs the Hop community code for AMUSE.
This package installs the Huyano community code for AMUSE.
This package installs the SPHRay community code for AMUSE.
This package installs the phiGRAPE community code for AMUSE.
This package installs the SSE community code for AMUSE.
This package installs the fractalcluster community code for AMUSE.
This package installs the Kepler community code for AMUSE.
This package installs the twobody community code for AMUSE.
This package installs the BSE community code for AMUSE.
This package installs the Brutus community code for AMUSE.
This package installs the smalln community code for AMUSE.
This package installs the Simplex community code for AMUSE.
This package installs the PeTar community code for AMUSE.
This package installs the AarsethZare community code for AMUSE.
This package installs the framework for the Astrophysical Multipurpose Software Environment (AMUSE).
This package installs the BHTree community code for AMUSE.
This package installs the Kepler-orbiters community code for AMUSE.
This package installs the Halogen community code for AMUSE.
This package installs the MMAMS community code for AMUSE.
This package installs the fastkick community code for AMUSE.
This package installs the distributed community code for AMUSE needed for
the distributed communication channel. 
This package installs the MOSSE community code for AMUSE.
This package installs the SeBa community code for AMUSE.
This package installs the tutorial for the Astrophysical Multipurpose Software Environment (AMUSE).
This package installs the Galaxia community code for AMUSE.
This package installs the Galactics community code for AMUSE.
This package installs the MOBSE community code for AMUSE.
