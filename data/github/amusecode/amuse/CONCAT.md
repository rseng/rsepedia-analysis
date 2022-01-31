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
.. AMUSE documentation master file, created by
   sphinx-quickstart on Tue Sep 29 13:22:44 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Astrophysical Multipurpose Software Environment
===================================================

Contents:

.. toctree::
   :maxdepth: 2
   
   install/index
   tutorial/index
   interactive_tutorial/index
   reference/index
   design/index
   
.. htmlonly::
   - `Examples <examples/index.html>`_
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

===================
Plotting with amuse
===================

matplotlib
==========

`Matplotlib <https://matplotlib.org/index.html>`_ is a
python plotting library capable of working with many  graphical user
interface toolkits. It is not required by AMUSE, but  if installed
then AMUSE provides extended plot functionality.  If a plot is made,
axis labels will be made automatically yielding  the concerning
units. 

To use matplotlib within AMUSE with the extended plot functionality
you need the following import:

.. code-block:: python

    >>> from amuse.plot import *

The native matplotlib plot functions are still available in the
native_plot namespace, e.g.:

.. code-block:: python

    >>> native_plot.subplot(2,2,1)

`matplotlib documentation <https://matplotlib.org/contents.html>`_

install matplotlib
------------------

Either use a pre-packaged version or install from source. If you
install from  source and you have installed the prerequisites in a
user directory make sure the ``PATH`` settings are correct. 

The source can be found `here
<https://github.com/matplotlib/matplotlib>`_

Installation instructions can be found `here
<https://matplotlib.org/users/installing.html>`_

Latex support
-------------
 
Latex support for labels can be enabled by issuing:

.. code-block:: python

    >>> latex_support()

This command will temporarily change the matplotlibrc settings:

.. code-block:: python

    rc('text', usetex=True)

Mathtext
~~~~~~~~

Latex support, while flexible and rich in features, might be slow and
requires latex, dvipng and Ghostscript to be installed. You can use
matplotlib's builtin TeX expression parser, mathtext instead. This is
a subset TeX markup usable in any matplotlib text string by placing it
inside a pair of dollar signs ($).  See `mathtext
<https://matplotlib.org/users/usetex.html>`_ for details.

Supported functions
-------------------

* plot
* semilogx
* semilogy
* loglog
* scatter
* hist
* xlabel
* ylabel

Example code
------------

.. literalinclude:: ../../examples/applications/test_plot.py

.. image:: plots/plot1.png
   :width: 18cm
   :align: left
=================================================================
Estimating the typical mass of the most massive star in a cluster
=================================================================

In this tutorial we will estimate the typical mass of the most 
massive star for a cluster with a Salpeter initial mass function (IMF).
This tutorial also illustrates that the units won't get in the way of 
calculations.

First we import numpy and set the seed of its random number 
generator (RNG). All random numbers used in AMUSE are drawn from this 
RNG. Seeding it is of course not necessary, but it will make the 
results reproducible.

.. code-block:: python

    >>> import numpy
    >>> numpy.random.seed(123456)

We will also need units and the Salpeter IMF generator. The easiest 
way to import almost everything from the AMUSE toolbox is to do:

.. code-block:: python

    >>> from amuse.lab import *

But here we will import the required module and function manually for clarity:

.. code-block:: python

    >>> from amuse.units import units
    >>> from amuse.ext.salpeter import new_salpeter_mass_distribution

Now, we can calculate the maximum stellar mass for each of 
*number_of_cluster_realizations* realizations of an 
*number_of_stars* star cluster. The Salpeter IMF runs from 0.1 to 
125 solar mass with a slope of -2.35, by default. Suppose the 
maximum of 125 solar mass is a bit too high for our taste, so we set 
it to 100 solar mass.

.. code-block:: python

    >>> number_of_stars = 1000
    >>> number_of_cluster_realizations = 100
    >>> maxmasses = [] | units.MSun
    >>> for i in range(number_of_cluster_realizations):
    ...     maxmasses.append(max(new_salpeter_mass_distribution(
    ...         number_of_stars, 
    ...         mass_max = 100. | units.MSun
    ...     )))
    ...

Note that the way we initialize *maxmasses*, with the solar mass 
unit, forces it to be a VectorQuantity instead of a Python list of 
ScalarQuantities. This is always recommended, because it is much 
faster, and it will make sure that AMUSE always recognizes it as a 
Quantity.  

If we want to know the mean mass of the *maxmasses* VectorQuantity, 
we simply use the *mean* function of VectorQuantities:

.. code-block:: python

    >>> print ("mean:  ", maxmasses.mean())
    mean:   27.4915750164 MSun

The same works for the median or the standard deviation of *maxmasses*.

.. code-block:: python

    >>> print ("median:", maxmasses.median())
    >>> print ("stddev:", maxmasses.std())
    median: 21.0983403429 MSun
    stddev: 19.7149800906 MSun

Slightly slower, but giving the same result, we can use the numpy functions:

.. code-block:: python

    >>> print ("mean:  ", numpy.mean( maxmasses))
    >>> print ("median:", numpy.median(maxmasses))
    >>> print ("stddev:", numpy.std(maxmasses))
    mean:   27.4915750164 MSun
    median: 21.0983403429 1.98892e+30 * kg
    stddev: 19.7149800906 MSun

Something weird has happened to the unit of the median 
mass. The result is still correct but the unit is converted to SI 
units. This is usually caused by a multiplication of a Quantity, 
where AMUSE tries to simplify the result, cancelling out for example 
factors of kg / kg. There's no need to bother, but if it annoys you, 
it can easily be fixed by:

.. code-block:: python

    >>> print ("median:", numpy.median(maxmasses).in_(units.MSun))
    median: 21.0983403429 MSun



===========================
Integrate a Fortran 90 code
===========================

In this tutorial we will create an AMUSE interface to a fortran 90 
code. We will first define the legacy interface, then implement the 
code and finally build an object oriented interface on top of the 
legacy interface.

The legacy code interface supports methods that transfer values to 
and from the code. The values do not have any units and no error 
handling is provided by the interface. We can add error handling, unit 
handling and more functions to the legacy interface by defining a 
subclass of a InCodeComponentImplementation (this is the objected oriented interface).

.. graphviz::

   digraph layers4 {
      fontsize=10.0;
        rankdir="LR";
        node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow];
        
        "Legacy Code" -> "Legacy Interface" -> "Object Oriented Interface" -> "Script";
    }

The legacy code in this tutorial will be a very simple and naive
implementation to find 3 nearest neighbors of a particle.

Two paths
---------
When defining the interface will walk 2 paths:

1. Management of particles in AMUSE (python)
2. Management of particles in the code (C or Fortran)

The first path makes sense for legacy codes that perform a 
transformation on the particles, or analyse the particles state or 
do not store any internal state between function calls (all data is 
external). For every function of the code, data of every particle is 
send to the code. If we expect multiple calls, the code would incur a 
high communication overhead and we are better of choosing path 2.

The second path makes sense for codes that already have management 
of a particles (or grid) or were we want to call multiple functions 
of the code and need to send the complete model to code for every 
function call. The code is first given the data, then calls are made 
to the code to evolve it's model or perform reduction steps on the 
data, finally the updated data is retrieved from the code. 

Procedure
---------

The suggested procedure for creating a new interface is as follows:

0. **Legacy Interface.** Start with creating the legacy 
   interface. Define functions on the interface to input and
   output relevant data.
   The InCodeComponentImplementation code depends on the legacy interface code.   
1. **Make a Class.** Create a subclass of the InCodeComponentImplementation class
2. **Define methods.** In the legacy interface we have defined functions
   with parameters. In the code interface we need to define the
   units of the parameters and if a parameter or return value
   is used as an errorcode.
3. **Define properties.** Some functions in the legacy interface can
   be better described as a property of the code. These are read only 
   variables, like the current model time.
4. **Define parameters.** Some functions in the legacy interface provide
   access to parameters of the code. Units and default values
   need to be defined for the parameters in this step
5. **Define sets or grids.** A code usually handles objects or gridpoints with
   attributes. In this step a generic interface is defined for these
   objects so that the interoperability between codes increases.

Before we start
---------------

This tutorial assumes you have a working amuse environment. Please 
ensure that amuse is setup correctly by running 'nosetests' in the 
amuse directory.

Environment variables
~~~~~~~~~~~~~~~~~~~~~
To simplify the work in the coming sections, we first define the 
environment variable 'AMUSE_DIR'. This environment variable must 
point to the root directory of AMUSE (this is the directory 
containing the build.py script).

.. code-block:: bash

    > export AMUSE_DIR=<path to the amuse root directory>
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR <path to the amuse root directory>

After building the code, we want to run and test the code. Check if 
amuse is available in your python path by running the following code 
on the command line.

.. code-block:: bash

    > python -c "import amuse"
    Traceback (most recent call last):
    File "<string>", line 1, in <module>
    ImportError: No module named amuse
    
If this code ends in a "ImportError" as shown in the example, the 
PYTHONPATH environment variable must be extended with the src directory
in AMUSE_DIR. 
We can do so by using one of the following commands.

.. code-block:: bash

    > export PYTHONPATH=${PYTHONPATH}:${AMUSE_DIR}/src
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR ${PYTHONPATH}:${AMUSE_DIR}/src
    
    
The name of our project
~~~~~~~~~~~~~~~~~~~~~~~
We will be writing a code to find the nearest neighbors of a particle, 
so let's call our project 'NearestNeighbor'.

Creating the initial directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First we need to create a directory for our project and put some 
files in it to help build the code. The fastest method to setup the 
directory is by using the build.py script.

.. code-block:: bash

    > $AMUSE_DIR/build.py --type=f90 --mode=dir NearestNeighbor

The script will generate a directory with all the files needed to 
start our project. It has also generates a very small legacy code 
with only one function ```echo_int```. We can build and test our new 
module::

    > cd nearestneighbor/
    > make all
    > $AMUSE_DIR/amuse.sh -c 'from interface import NearestNeighbor; print NearestNeighbor().echo_int(10)' 
    OrderedDictionary({'int_out':10, '__result':0})
    > nosetests -v
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.556s

    OK
    
    
.. note::

    The build.py script can be used to generate a range of files. To
    see what this file can do you can run the script with a --help
    parameter, like so::

        > $AMUSE_DIR/build.py --help

The Legacy Code
---------------
Normally the legacy code already exists and our task is limited to 
defining and implementing an interface so that AMUSE scripts can 
access the code. For this tutorial we will implement our legacy code.

When a legacy code is integrated all interface code is put in one 
directory and all the legacy code is put in a **src** directory 
placed under this directory. The build.py script created a **src** 
directory for us, and we will put the nearest neighbor algorithm in 
this directory.

Go to the **src** directory and create a **code.f90** file, open 
this file in your favorite editor and copy and paste this code into 
it:
    
.. literalinclude:: nearestneighbor/code.f90
    :language: fortran
    

.. note::

    This algorithm is un-optimized and has N*N order. It is not meant
    as very efficient code but as a readable example.

Before we can continue we also need to alter the **Makefile** in the 
**src** directory, so that our **code.cc** file is included in the 
build. To do so, open an editor on the Makefile and change the line::

    CODEOBJS = test.o

to::

    CODEOBJS = test.o code.o
    

Test if the code builds. As we have not coupled our algorithm to the 
interface we (we have not even defined an interface) we do not have 
any new functionality. In the legacy interface directory (not the 
**src** directory) do:

.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works, if the test fails for any reason please check that the fortran
code is correct and that ``worker_code`` exists in your directory.

Path 1
------

Defining the legacy interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will first define a legacy interface so that we can call the 
**find_nearest_neighbors** function from python. AMUSE can interact 
with 2 classes of functions:

1. A function with all scalar input and output variables. All variables
   are simple, non-composite variables (like INTEGER or DOUBLE PRECISION).
   For example:
   
   .. code-block:: fortran
    
        FUNCTION example1(input, output)
            IMPLICIT NONE
            INTEGER example1
            DOUBLE PRECISION input, output
            output = input
            example1 = 0
        END FUNCTION
        
2. A function with all vector (or list) input and output variables and
   a length variable. The return value is a scalar value.
   For example:
   
   .. code-block:: fortran
    
        FUNCTION example2(input, output, N)
            IMPLICIT NONE
            INTEGER example2, N, i
            DOUBLE PRECISION input(N), output(N)
            DO i = 1, N
                output(i) = input(i)
            END DO
            example2 = 0
        END FUNCTION

If you have functions that don't follow this pattern you need to 
define a convert function in fortran that provides an interface 
following one of the two patterns supported by AMUSE.

In our case the **find_nearest_neighbors** complies to pattern 2 and
we do not have to write any code in fortran to convert the function
to a compliant interface. We only have to specify the function in 
python. We do so by adding a **find_nearest_neighbors** method to
the **NearestNeighborInterface** class in the **interface.py** file.
Open and editor on the interfaces.py file and add the following method
to the NearestNeighborInterface class:


.. code-block:: python

    class NearestNeighborInterface(InCodeComponentImplementation):    
        #...
        
        @legacy_function
        def find_nearest_neighbors():
            function = LegacyFunctionSpecification()  
            function.must_handle_array = True 
            function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
            function.addParameter('x', dtype='float64', direction=function.IN)
            function.addParameter('y', dtype='float64', direction=function.IN)
            function.addParameter('z', dtype='float64', direction=function.IN)
            function.addParameter('n0', dtype='int32', direction=function.OUT)
            function.addParameter('n1', dtype='int32', direction=function.OUT)
            function.addParameter('n2', dtype='int32', direction=function.OUT)
            function.result_type = 'int32'
            return function
            
In the *find_nearest_neighbors* method we specify every parameter of 
the fortran function and the result type. For each parameter we need 
to define a name, data type and wether we will input, output (or 
input and output) data using this parameter. AMUSE knows only a 
limited amount of data types for parameters: float64, float32, int32 
and string. We also have a special parameter, with LENGTH as 
direction. This parameter is needed for all functions that follow 
pattern 2, it will be filled with the length of the input arrays. We 
also must specify that the function follows pattern 2 by setting 
```function.must_handle_array = True```.

Save the file and recompile the code.

.. code-block:: bash

    > make clean
    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works! But, how do we know the find_nearest_neighbors method really works?
Let's write a test and find out. Open an editor on the *test_nearestneighbor.py*
file and add the following method to the **NearestNeighborInterfaceTests** class:

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        x = [0.0, 1.0, 4.0, 7.5]
        y = [0.0] * len(x)
        z = [0.0] * len(x)
        
        n0, n1, n2, error = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(error[0], 0)
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()

This test calls the *find_nearest_neighbors* method with 4 positions
and checks if the nearest neighbors are determined correctly.
Let's run the test, and see if everything is working:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 test in 0.491s

    OK

We now have a simple interface that works, but we have to do our own 
indexing after the call and we could send data of any unit to the
method, also we have to do our own error checking after the method. Let's
define a object oriented interface to solve these problems

Defining the Object Oriented Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The object oriented interface sits on top of the legacy interface.
It decorates this interface with sets, unit handling, state engine
and more. We start creating the object oriented interface by inheriting
from InCodeComponentImplementation and writing the __init__ function. The build script
has added this class to the interface.py file for us. Open an editor
on interface.py and make sure this code is in the file (at the end
of the file):

.. code-block:: python

    
    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)

Configuring the handlers
~~~~~~~~~~~~~~~~~~~~~~~~

We configure the object oriented interface by implementing several 
methods. The object oriented interface is implement by several 
"handlers". Each handler provides support for a specific aspect of 
the interface. AMUSE defines a handler for the unit conversion, a 
handler for the interfacing with sets of particles, a handler to 
ensure the methods are called in the right order, etc. Each handler 
is very generic and needs to be configured before use. The handler 
are configured using the "Visitor" pattern. The following 
pseudo-code shows how the handlers are configured

.. code-block:: python

    class InCodeComponentImplementation(object):
        #...
        
        def configure_handlers(self):
            #...
            for handler in self.get_all_handlers():
                handler.configure(self)
        
        def define_converter(self, handler):
            """ configure the units converter handler """
            
            handler.set_nbody_converter(...)
            
        def define_particle_sets(self, handler):
            """ configure sets of particles """
            
            handler.define_incode_particle_set(...)
            handler.set_getter(...)
            
    class HandleConvertUnits(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_converter(self)
            
    class HandleParticles(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_particle_sets(self)
            
            
Configuration of the handlers is optional, we only have to define
those handler that we need in our interface. In our example we need
to configure the "HandleMethodsWithUnits" handler (to define units and
error handling) and the "HandleParticles" to define a particle set.

Defining methods with units
+++++++++++++++++++++++++++

We first want to add units and error handling to the 
**find_nearest_neighbors**. We do this by creating a 
**define_methods** function on the **NearestNeighbor** class. Open 
an editor on *interface.py* and add this method to the class:

.. code-block:: python

    def define_methods(self, handler):
        
        handler.add_method(
            "find_nearest_neighbors",
            (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                handler.INDEX,
                handler.INDEX,
                handler.INDEX,
                handler.ERROR_CODE
            )
        )
        
The **add_method** call expects the name of the function in the
legacy interface as it's first, next it expects a list of the units
of the input parameters and a list of the units of the output parameters.
The return value of a function is always the last item in the list
of output parameters. We specify a **generic_unit_system.length** unit
for the x, y and z parameters. The output parameters are indices and 
an errorcode. The errorcode will be handled by the AMUSE interface (0
means success and < 0 means raise an exception).

Let's write a test to see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test3(self):
        instance = NearestNeighbor()
        x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        y = [0.0] * len(x) | generic_unit_system.length
        z = [0.0] * len(x) | generic_unit_system.length
        
        n0, n1, n2 = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()

.. note::
    This test looks a lot like test2, but we now have to
    define a unit and we do not need to handle the errorcode.
    
Now build and test the code:

.. code-block:: bash
    
    > make clean; make all
    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.650s

    OK
    
.. note::
    
    Although we only edited python code we still need to run make. The
    code will check if the "worker_code" executable is up to date on 
    every run. It cannot detect if the update broke the code but it
    will still demand that the code is rebuilt.
    

Defining the particle set
+++++++++++++++++++++++++

Particle sets in AMUSE can be handled by python (we call these 
"inmemory") and by the community code (we call these "incode"). In our 
case the code does not handle the particles and we need to 
configure the particles handler to manage an inmemory particle set. 
Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    def define_particle_sets(self, object):
        object.define_inmemory_set('particles')

That's all we now have a "particles" attribute on the class and 
we can add, remove, delete particles from this set. But we are
still missing a connection between the particles and the nearest 
neighbors. AMUSE provides no handler for this, instead, we 
will write a method to run the find_nearest_neighbors function and
set the indices on the particles set.

Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    
    def run(self):
        indices0, indices1, indices2 = self.find_nearest_neighbors(
            self.particles.x,
            self.particles.y,
            self.particles.z
        )
        self.particles.neighbor0 = list(self.particles[indices0 - 1])
        self.particles.neighbor1 = list(self.particles[indices1 - 1])
        self.particles.neighbor2 = list(self.particles[indices2 - 1])
        
This function gets the "x", "y" and "z" attributes from the particles
set and sends these to the "find_nearest_neighbors" method. This methods
returns 3 lists of indices and we need to find the particles with
these indices. As this is fortran code (indices start with 1) and we use
python (indices start with 0) we need to subtract 1 from the array of
indices and use this to find the particles.

.. note::
    
    Particle sets have no given sequence, deletion and addition of 
    particles will change the order of the particles in the set. It
    is therefor never a good idea to use th index of the particle in the
    set as a reference to that particle. However, in the "run" method
    we "own" the particle set, it cannot change between the find_nearest_neighbor
    call and the moment we find the particles in the set by index (using
    self.particles[indices0-1]), and in this case it is save to use
    index as a valid reference.

Let's write a test and see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test4(self):
        instance = NearestNeighbor()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()

Now, make and run the tests:

.. code-block:: bash

    > make clean; make all
    > nosetests
    ....
    ----------------------------------------------------------------------
    Ran 4 tests in 0.797s

    OK
 
We are done, we have defined an object oriented interface on the
legacy interface. Only, if we look at our tests, the code seems
to be more rather than less complex. But, remember we now
have units and we are compatible with other parts of amuse. And
we can make more complex scripts easier.

Let's make a plummer model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash

    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png



Path 2
------

Defining the legacy interface
-----------------------------
We define our code interface so that a user can add, update and 
delete particles, start the nearest neighbors finding
algorithm and retrieve the ids of the nearest neighbors.

To define the interface, open interface.py with your favorite
editor and replace the contents of this file with:

.. literalinclude:: nearestneighbor/nn2.py

We can generate a stub from the interface code with::

    > $AMUSE_DIR/build.py --type=f90 --mode=stub interface.py NearestNeighborInterface -o interface.f90

The generated **interface.f90** replaces the original file generated in
the previous section. 

We will be implementing the interface.f90 file as a module and we 
need to add the module definition to this file. We do 
this by adding a ```MODULE NN``` line to the beginning of the file 
and a ```END MODULE``` line to the end of the file. To do so open 
the interface.f90 file and append/prepend the following code:

.. code-block:: fortran

    MODULE NN

    CONTAINS

    !.. original code
    
    END MODULE
    
.. note::
    We need to remove the original USE NN line at the start
    of the file.

The code builds, but does not have any functionality yet::

    > make clean
    > make all
    
.. note::
    
    Compiling the interface code will result in a lot of warnings about
    unused dummy arguments. These warnings can be savely ignored for now.
    
The tests are broken (the echo_int function has been removed):

.. code-block:: bash

    > nosetests
    E
    ======================================================================
    ERROR: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "../src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "nearestneighbor/test_nearestneighbor.py", line 11, in test1
        result,error = instance.echo_int(12)
    AttributeError: 'NearestNeighborInterface' object has no attribute 'echo_int'

    ----------------------------------------------------------------------
    Ran 1 test in 0.315s

    FAILED (errors=1)

Let's create a working test by calling the new_particle method, open an editor
on the test_nearestneighbor.py file and replace the test1 method with:

.. code-block:: python

    def test1(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()

As this is python code we do not need to rebuild the code, instead
we can run the tests right after saving the code. Unfortunately, when
we run the test, it still fails.

.. code-block:: bash

    > nosetests
    F
    ======================================================================
    FAIL: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "/src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "/nearestneighbor/test_nearestneighbor.py", line 13, in test1
        self.assertEquals(result, 1)
      File "/src/amuse/test/amusetest.py", line 62, in failUnlessEqual
        self._raise_exceptions_if_any(failures, first, second, '{0} != {1}', msg)
      File "/src/amuse/test/amusetest.py", line 49, in _raise_exceptions_if_any
        raise self.failureException(msg or err_fmt_string.format(first, second, *args))
    AssertionError: 0 != 1
    -------------------- >> begin captured logging << --------------------
    legacy: INFO: start call 'NearestNeighborInterface.new_particle'
    legacy: INFO: end call 'NearestNeighborInterface.new_particle'
    --------------------- >> end captured logging << ---------------------

    ----------------------------------------------------------------------
    Ran 1 test in 0.319s

When you look closely at the output of the test you see that the result from the
method is 0 and not the excepted 1. We need to edit the fortran code to make this 
test work. Open an editor on interface.f90 and go to the *new_particle* function.

.. code-block:: fortran

    FUNCTION new_particle(index_of_the_particle, x, y, z)
      INTEGER :: index_of_the_particle
      DOUBLE PRECISION :: x, y, z
      index_of_the_particle = 1
      new_particle=0
    END FUNCTION

.. note ::
    
    In AMUSE all interface functions return an errorcode. Any other return
    values must be passed through the arguments of the functions.


We need to rebuild the code, and after building run the tests.


.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

The are tests work again! Only, we do not have any real working legacy code.

Filling the stubs
-----------------
The implementation of the algorithm does not match the interface
we defined and created. We need to write some glue code to
connect the code with the interface. To do so we fill in the
stubs generated earlier.

Open the **interface.f90** file in your favorite editor and
change its contents to:

.. literalinclude:: nearestneighbor/interface1.f90
    :language: fortran

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all
    
The tests will still fail as we need to set the 
"maximum_number_of_particles" before allocating any new particles.
Let's update the test, open an editor on test_nearestneighbor.py 
file and update the test1 method to:

.. code-block:: python

    def test1(self):
        instance = NearestNeighborInterface()
        instance.set_maximum_number_of_particles(100)
        instance.commit_parameters()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()

Now the tests should succeed:

.. code-block:: bash

    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.311s

    OK
    
Let's check some more functionaly by adding another test

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        instance.set_maximum_number_of_particles(2)
        instance.commit_parameters()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 2)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, -1)
        error = instance.delete_particle(1)
        self.assertEquals(error, 0)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()
        

Now the tests should succeed:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 tests in 0.448s

    OK


We now have done everything in Step 0 *Legacy Interface*. We have
a legacy code and can access it in our python script. But, our
interface is not very friendly to work with. We have to think about
errorcodes and we have not information about units. To make our
interface easier to works with we start defining methods, properties
and parameters.

Defining methods
----------------
The object oriented interface is also defined in the **interface.py**.
So, we continue by opening an editor on this file. We will be 
writing methods for the **NearestNeighbor** class, in your editor
seek this code (at the end of the file)::

    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)
            
We will start by defining methods, we will do this by implementing
the **define_methods** function, like so::

    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            
            builder.add_method(
                "new_particle", 
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.INDEX, builder.ERROR_CODE)
            )
            
            builder.add_method(
                "delete_particle", 
                (builder.INDEX,),
                (builder.ERROR_CODE)
            )
            
            builder.add_method(
                "get_state", 
                (builder.INDEX,),
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, builder.ERROR_CODE),
                public_name = "get_position"
            )
            
            builder.add_method(
                "set_state", 
                (builder.INDEX, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.ERROR_CODE),
                public_name = "set_position"
            )
            
            builder.add_method(
                "run", 
                (),
                (builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_close_neighbors", 
                (builder.INDEX,),
                (builder.LINK('particles'), builder.LINK('particles'), builder.LINK('particles'), builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_nearest_neighbor", 
                (builder.INDEX,),
                (builder.LINK('particles'), generic_unit_system.length, builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_number_of_particles", 
                (),
                (builder.NO_UNIT, builder.ERROR_CODE),
            )

        
With this code, we define the methods and specify how to interpret the
arguments and return values. We get a special object (the **builder**
object) that provides us with the **add_method** function to be able
to this. The definition of the **add_method** function is as follows::
    
    add_method(
        name of the original function in the legacy interface,
        list of arguments (unit or type),
        list of return values,
        public_name = name for the user of the class (optional)
    )
    
For every argument or return value we can specify if it has a unit or
if it is special. The special arguments are:

=========================== ===================================
definition                  description
=========================== ===================================
builder.ERROR_CODE          the value is interpreted as an errorcode,
                            zero means no error for all other values and
                            Exception will be raise, only valid for one return value.
builder.NO_UNIT             the value has not unit (for example for the
                            number of items in a list)
builder.INDEX               the value is interpreted as an index for
                            object identifiers
builder.LINK('particles')   the value is interpreted as a link to
                            objects in the set with the given name
=========================== ===================================

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test::

    def test3(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(2)
        instance.commit_parameters()
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length, 
            3.0 | generic_unit_system.length
        )
        self.assertEquals(result, 1)
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length
        )
        self.assertEquals(result, 2)
        x,y,z = instance.get_position(1)
        self.assertEquals(1.0 | generic_unit_system.length, x)
        self.assertEquals(2.0 | generic_unit_system.length, y)
        self.assertEquals(3.0 | generic_unit_system.length, z)
        instance.stop()

And run the tests:

.. code-block:: bash

    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.664s

    OK

As you can see our script is now a little simpler and we support units.
We do not have to think about the errorcodes in this script, AMUSE will
interpret the errorcodes and raise the right exceptions if needed. The 
units are also automatically converted to the right units for the 
code. But the script is still not very easy and we have to manage
all the ids we get from the code. To make our code even easier to
handle we will continue by defining a **set**.

.. note::

    We skip defining parameters and properties, we will come back to 
    this later in this tutorial.

Defining a set
--------------
We have made our interface a little easier but we still have to
do a some management work in our script. We would like to work
with objects and adding or removing these objects from the code.
AMUSE supports this by defining **sets**. Each set is capable of
storing specific attributes of the objects in the set. Our code
is capable of storing the x, y and z position of an object. An object
in AMUSE is called a *Particle* and the sets that contain these
particles are called *ParticleSets* or shorter *Particles*.

We define our particle set by implementing a **define_particle_sets** 
function on our **NearestNeighbor** class like so::

    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            ...
            
        def define_particle_sets(self, builder):
            builder.define_set('particles', 'index_of_the_particle')
            builder.set_new('particles', 'new_particle')
            builder.set_delete('particles', 'delete_particle')
            builder.add_setter('particles', 'set_position')
            builder.add_getter('particles', 'get_position')
            builder.add_getter('particles', 'get_close_neighbors', names=('neighbor0', 'neighbor1', 'neighbor2') )


That's all, we now have defined a set called **particles**. Again, we 
get a builder object to use in defining our set. All methods have
the name of the set as their first argument, this name can be any
name you want, but in AMUSE most codes provide a set called 
**particles**. For the **add_setter**, **add_getter**, **set_new** 
and **set_delete** functions, the second argument is the name of
the method we defined in the previous step. Finally you can set
the name of the attribute in the particles set with the **names**
argument. This is optional for legacy functions, if not given the names 
of the attributes will be derived from the names of the arguments in
the original calls. For example, the **get_position** call we specified
earlier has parameter name **x**, **y** and **z**, these names are 
also used in the particles set.

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test::
    
    def test4(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(100)
        instance.commit_parameters()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()
    
    
The support for a particle set means we can now also interact
with other parts of AMUSE. Let's make a plummer
model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash

    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png

    

    
    









===========================================
Using Blender for visualising Amuse results
===========================================

.. image:: blender.jpeg
   :width: 2cm

Blender
=======

Blender is a free open source 3D content creation suite, *GPL-ed*. As it is specialised in visualisation of 3D content *and* 
scriptable in python it provides handy a tool for visualisation of Amuse data.

More information on Blender can be found on https://www.blender.org/

Starting Blender
================
We start blender from the command line, assuming all environment variables are set for amuse, and make sure mpd is running.

.. code-block:: bash

    >>blender &

After blender opened we press **CTRL-RIGHTCURSOR** three times to open a text panel next to our 3D view. In the text panel we use the Text 
menu item to open a text file:

.. image:: blenderopenscript.png
   :width: 800
   :align: center

This way we can retrieve scripts we had written before. Now, we are going to write a new script, but use our favourite editor, and once finished 
we will load the script using the method described above.

**Note:** By default, blender places a cube in the scene. You can delete it else it will obscure the sun in our example. 

The script is based on the sun-earth test in Hermite. We start with that amuse script and add a couple of lines to communicate with blender, which 
are marked by comments in the example below.

Amuse blender API
=================

To simplify communication with blender, amuse has the module amuse.ext.blender, which contains a very basic API. In our example we start with importing 
it on top of the regular stellar dynamics imports:

.. code-block:: python

    from amuse.ext.blender import blender 

To create a sphere, e.g. the sun,  with 32 segments, 32 rings and a radius of 1.0 in the current scene we type:

.. code-block:: python

    sun = blender.Primitives.sphere(segments = 32, rings = 32, radius = 1.0)

We can move and rotate our object using:

.. code-block:: python

    x,y,z = 1.0, 0.0, 0.0
    alpha = 0.2
    sun.loc = (x, y, z)
    sun.rotZ = alpha


The Hermite Sun-Earth model with blender visualisation will become like this:

.. literalinclude:: ../../examples/applications/test_blender.py

The path to this example code is {amusedir}/trunk/examples/applications/test_blender.py
      
We save this file as *myname.py*, open it in blender and run it by typing **ALT-P** (cursor in editor):

.. image:: blenderrunscript.png
   :width: 800
   :align: center 

=======================================
Create an low-level Interface to a Code
=======================================

In this tutorial we will create an interface to a code. 

.. warning::
    
    This tutorial does not use the build script provided with 
    AMUSE. All files and directories need to be created "by hand".
    Use this tutorial if you want to get a deeper understanding of
    how the build process works and which steps are involved in
    creating a low level interface. To learn how to create
    an interface to a code we recommend :doc:`c_code` and 
    :doc:`fortran_code`.


Work directory
~~~~~~~~~~~~~~
We start by making a directory for our code. This directory should
be a subdirectory of the "src/amuse/community" directory. It also will be
a python package and we need to create the file "__init__.py" in 
this directory. So, let's open a shell and go to the AMUSE 
root directory. To create the code directory we then do:

.. code-block:: bash
    
    >> cd src/amuse/community
    >> mkdir mycode
    >> cd mycode
    >> touch __init__.py
    >> pwd
    ../src/amuse/community/mycode
    >> ls
    __init__.py
    

The code
~~~~~~~~
To create an interface we first need the code. For
this example we will use a very simple code do some calculations 
on two numbers.

The contents of the code.c file:

.. code-block:: cpp

  #include "code.h"
   
  double sum(double x, double y) {
    return x + y;
  }
  
  int divide(double x, double y, double * result) {
    if(y == 0.0) {
        return -1;
    } else {
        *result = x / y;
        return 0;
    }
  }

We need to access these function from another C file, so we need to
define a header file. 

The contents of the code.h file:

.. code-block:: cpp

    double sum(double x, double y);
    
    int divide(double x, double y, double * result);


The interface code
~~~~~~~~~~~~~~~~~~
Now we can define the interface class for our code in python. An 
interface needs to inherit from the class "CodeInterface".

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
In this example we first import names from the ``amuse.community`` 
module on line 1. The ``amuse.community`` module defines the typical 
classes and function needed to write a legacy interface. On line 3 
we define our class and inherit from ``CodeInterface``. The class 
will be used to generate a C++ file. In this file we will need the 
definitions of our legacy functions. So, on line 4 we specify the 
necessary include files in a array of strings. Each string will be
converted to an include statement.

Building the code
~~~~~~~~~~~~~~~~~
Building the code takes a couple of steps, first generating the C file
and then compiling the code. We will construct a makefile to automate
the build process.

.. code-block:: make
    :linenos:
    
    ifndef AMUSE_DIR
        AMUSE_DIR=../../../..
    endif

    CODE_GENERATOR = $(AMUSE_DIR)/build.py

    CXXFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
    LDFLAGS = -lm $(MUSE_LD_FLAGS)

    OBJS = code.o

    all: worker_code

    cleanall: clean
        $(RM) worker_code *~
        
    clean:
        rm -f *.so *.o *.pyc worker_code.cc

    worker_code.cc: interface.py
        $(CODE_GENERATOR) --type=c interface.py MyCode -o $@

    worker_code: worker_code.cc $(OBJS)
        mpicxx $@.cc $(OBJS) -o $@

    .cc.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<
        
    .c.o: $<
        g++ $(CXXFLAGS) -c -o $@ $<

.. compound:

    You need to convert the spaces into tabs, 
    if you copy the above text to a new file.
    

Let's start make and build the ``worker_code`` application

.. code-block:: bash
    
    >> make clean
    >> make
    ...
    mpicxx worker_code.cc code.o -o worker_code
    >> ls
    ... worker_code ...
    
Running the code
~~~~~~~~~~~~~~~~
Before we run the code we need to start the MPI daemon process ''mpd''. 
This daemon process manages the start of child processes.

.. code-block:: bash
    
    >> mpd &

We can use ``amuse.sh`` and try the interface.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance
    <amuse.community.mycode.interface.MyCode object at 0x7f57abfb2550>
    >>> del instance
    
We have not defined any methods and our interface class is not
very useful. We can only create an instance of the code. When we 
create this instance the "worker_code" application will start 
to handle all the function calls. We can see the application 
running when we do "ps x | grep worker_code"

Implementing a method
~~~~~~~~~~~~~~~~~~~~~~
Now we will define the ``sum`` method. We will add the definition to
the MyCode class.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function
            
The new code starts from line 9. On line 9 we specify a decorator. This
decorator will convert the following function into a specification that
can be used to call the function and generate the C++ code. On line 10
we give the function the same name as the function in our code. This
function may not have any arguments. On line 11 we create an instance of
the "LegacyFunctionSpecification" class, this class has methods to specify the intercase.
On line 12 and 13 we specify the parameters for out functions. Parameters have
a name, type and direction. The type is specified with a single character
*type code*. The following type codes are defined:
            
            
=========  ===========  ================
Type code  C type       Fortran type  
=========  ===========  ================
'i'        int          integer
'd'        double       double precision
'f'        float        single precision
=========  ===========  ================

The direction of the parameter can be ``IN``, ``OUT`` or ``INOUT``. On line 
14 we define the return type, this can be a *type code* or ``None``. The default
value is ``None``, specifying no return value (void function).

Let's rebuild the code.

.. code-block:: bash
    
    >> make clean
    >> make
    ...
    mpicxx worker_code.cc code.o -o worker_code

We can now start ```amuse.sh``` again and do a simple sum

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> instance.sum(40.5, 10.3)
    50.799999999999997
    >>> 40.5 + 10.3
    50.799999999999997
    >>> del instance

And we see that our interface correctly sums two numbers.

A method with an OUT parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We can complete out interface by defining the ``divide`` function.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            return function

        @legacy_function
        def divide():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.addParameter('result', 'd', function.OUT)
            function.result_type = 'i'
            return function
            
On line 22 we define the parameter "result" as an OUT parameter. In python
we do not have to provide this parameter as an argument to our function. It
After rebuilding we can try this new function.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> (result, error) =  instance.divide(10.2, 30.2)
    >>> result
    0.33774834437086093
    >>> error
    0
    >>> del instance

We see that the function returns two values, the OUT parameter and also
the return value of the function.

Working with arrays
~~~~~~~~~~~~~~~~~~~
Some functions will be called to perform on the elements of an array. 
For example:

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> x_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    >>> y_values = [10.3, 20.3, 30.4 , 40.4, 50.6]
    >>> results = []
    >>> for x , y in map(None, x_values, y_values):
    ...     results.append(instance.sum(x,y))
    ...
    >>> print results
    [11.300000000000001, 22.300000000000001, 33.399999999999999, 
    44.399999999999999, 55.600000000000001]
    
    
The MPI message passing overhead is incurred for every call on 
the method. We can change this by specifying the function to be able
to handle arrays.

.. code-block:: python
    :linenos:
    
    from amuse.community import *
    
    class MyCode(CodeInterface):
        include_headers = ['code.h']
        
        def __init__(self):
             CodeInterface.__init__(self)
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.result_type = 'd'
            function.can_handle_array = True
            return function

On line 15 we specify that the function can be called with an array of
values. The function will be called for every element of the array. The 
array will be send in one MPI message, reducing the overhead.

Let's rebuild the code and run an example.

.. code-block:: pycon

    >>> from amuse.community.mycode import interface
    >>> instance = interface.MyCode()
    >>> x_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    >>> y_values = [10.3, 20.3, 30.4 , 40.4, 50.6]
    >>> results = instance.sum(x_values, y_values)
    >>> print results
    [ 11.3  22.3  33.4  44.4  55.6]
    

Other interfaces
~~~~~~~~~~~~~~~~
The community codes directory contains a number of codes. Please look at
these codes to see how the interfaces are defined.


             


.. _working_with_units:

==================
Working with Units
==================

The AMUSE framework provides a unit handling library. This library
is used throughout the AMUSE framework. When interacting with a code
all data has a unit, even scaled systems have units.

==========
Quantities
==========

The basic data object is a quantity, a quantity is made up of 
a value and a unit. The value can be a single number (a scalar 
quantity) or a multi-dimensional array (a vector quantity). 

Quantities are created by the ```|(bar)``` operator. All 
quantities can be used like numbers and a lot of numpy 
functions also work on quantities.

.. code-block:: python
    
    >>> from amuse.units.si import *
    >>> from amuse.units.core import named_unit
    >>>
    >>> weigth = 80 | kg
    >>> persons = 10
    >>> print ("Total weight: ", persons * weigth)
    Total weight:  800 kg
    
    >>> day = named_unit("day", "d", s * 60 * 60 * 24 )
    >>> weight_loss = (0.1 | kg) / (1 | day)
    >>> print ("Weight loss: ", weight_loss)
    Weight loss:  0.1 1.15740740741e-05 * kg * s**-1
    >>> print ("Weight loss: ", weight_loss.as_quantity_in(kg/day))
    Weight loss:  0.1 kg / d
    
===================
Working with arrays
===================

A vector quantity can be used like a python list. Take care to only 
put quantities into a vector quantity. These vector quantities behave
more like numpy arrays than like python lists in that numeric operators
operate on them element-wise; they can also be initialized from numpy
arrays.

.. code-block:: python
    
    >>> from amuse.units.units import MSun
    >>>
    >>> masses = [] | MSun
    >>> for i in range(10):
    ...     masses.append(i**2 | MSun)
    >>> print ("Masses:", masses)
    Masses: [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0] MSun
    
.. note::

    When working with arrays, some care must be taken to ensure that
    vector quantities are created and not arrays of quantities. The
    following code will create an array of quantities
    
    >>> from amuse.units.units import MSun
    >>>
    >>> masses = []
    >>> for i in range(2):
    ...     masses.append(i**2 | MSun)
    >>> print ("Masses:", masses)
    Masses: [quantity<0 MSun>, quantity<1 MSun>]

    
====================
Integrate a C++ code
====================


In this tutorial we will create an AMUSE interface to a C++
code. We will first define the legacy interface, then implement the 
code and finally build an object oriented interface on top of the 
legacy interface.

The legacy code will be a very simple and naive implementation to find
3 nearest neighbors of a particle.

The legacy code interface supports methods that transfer values to 
and from the code. The values do not have any units and no error 
handling is provided by the interface. We can add error handling, unit 
handling and more functions to the legacy interface by defining a 
subclass of a InCodeComponentImplementation (this is the objected oriented interface).

.. graphviz::

   digraph layers4 {
      fontsize=10.0;
        rankdir="LR";
        node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow];
        
        "Legacy Code" -> "Legacy Interface" -> "Object Oriented Interface" -> "Script";
    }

The legacy code in this tutorial will be a very simple and naive
implementation to find 3 nearest neighbors of a particle.

Two paths
---------
When defining the interface will walk 2 paths:

1. Management of particles in AMUSE (python)
2. Management of particles in the code (C or Fortran)

The first path makes sense for legacy codes that perform a 
transformation on the particles, or analyse the particles state or 
do not store any internal state between function calls (all data is 
external). For every function of the code, data of every particle is 
send to the code. If we expect multiple calls, the code would incur a 
high communication overhead and we are better of choosing path 2.

The second path makes sense for codes that already have management 
of a particles (or grid) or were we want to call multiple functions 
of the code and need to send the complete model to code for every 
function call. The code is first given the data, then calls are made 
to the code to evolve it's model or perform reduction steps on the 
data, finally the updated data is retrieved from the code. 

Procedure
---------

The suggested procedure for creating a new interface is as follows:

0. **Legacy Interface.** Start with creating the legacy 
   interface. Define functions on the interface to input and
   output relevant data.
   The InCodeComponentImplementation code depends on the legacy interface code.   
1. **Make a Class.** Create a subclass of the InCodeComponentImplementation class
2. **Define methods.** In the legacy interface we have defined functions
   with parameters. In the code interface we need to define the
   units of the parameters and if a parameter or return value
   is used as an errorcode.
3. **Define properties.** Some functions in the legacy interface can
   be better described as a property of the code. These are read only 
   variables, like the current model time.
4. **Define parameters.** Some functions in the legacy interface provide
   access to parameters of the code. Units and default values
   need to be defined for the parameters in this step
5. **Define sets or grids.** A code usually handles objects or gridpoints with
   attributes. In this step a generic interface is defined for these
   objects so that the interoperability between codes increases.

Before we start
---------------

This tutorial assumes you have a working amuse development build. Please 
ensure that amuse is setup correctly by running 'nosetests' in the 
amuse directory.


Environment variables
~~~~~~~~~~~~~~~~~~~~~
To simplify the work in the coming sections, we first define the 
environment variable 'AMUSE_DIR'. This environment variable must 
point to the root directory of AMUSE (this is the directory 
containing the build.py script).

.. code-block:: bash

    > export AMUSE_DIR=<path to the amuse root directory>
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR <path to the amuse root directory>

After building the code, we want to run and test the code. Check if 
amuse is available in your python path by running the following code 
on the command line.

.. code-block:: bash

    > python -c "import amuse"
    Traceback (most recent call last):
    File "<string>", line 1, in <module>
    ImportError: No module named amuse
    
If this code ends in a "ImportError" as shown in the example, the 
PYTHONPATH environment variable must be extended with the src directory
in AMUSE_DIR. 
We can do so by using one of the following commands.

.. code-block:: bash

    > export PYTHONPATH=${PYTHONPATH}:${AMUSE_DIR}/src
    
or in a c shell:

.. code-block:: csh

    > setenv AMUSE_DIR ${PYTHONPATH}:${AMUSE_DIR}/src
    
    
The name of our project
~~~~~~~~~~~~~~~~~~~~~~~
We will be writing a code to find the nearest neighbors of a particle, 
so let's call our project 'NearestNeighbor'.

Creating the initial directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First we need to create a directory for our project and put some 
files in it to help build the code. The fastest method to setup the 
directory is by using the build.py script.

.. code-block:: bash

    > $AMUSE_DIR/build.py --type=c --mode=dir NearestNeighbor

The script will generate a directory with all the files needed to 
start our project. It also generates a very small example legacy code 
with only one function ```echo_int```. We can build and test our new 
module::

    > cd nearestneighbor/
    > make all
    > $AMUSE_DIR/amuse.sh -c 'from interface import NearestNeighbor; print NearestNeighbor().echo_int(10)' 
    OrderedDictionary({'int_out':10, '__result':0})
    > nosetests -v
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.556s

    OK
    
    
.. note::

    The build.py script can be used to generate a range of files. To
    see what this file can do you can run the script with a ```--help```
    parameter, like so::

        > $AMUSE_DIR/build.py --help

The Legacy Code
---------------
Normally the legacy code already exists and our task is limited to 
defining and implementing an interface so that AMUSE scripts can 
access the code. For this tutorial we will implement our legacy code.

When a legacy code is integrated all interface code is put in one 
directory and all the legacy code is put in a **src** directory 
placed under this directory. The build.py script created a **src** 
directory for us, and we will put the nearest neighbor algorithm in 
this directory.

Go to the **src** directory and create a **code.cc** file, open 
this file in your favorite editor and copy and paste this code into 
it:
    
.. literalinclude:: nearestneighbor/code.cc
    :language: c++
    

.. note::

    This algorithm is un-optimized and has N*N order. It is not meant
    as very efficient code but as a readable example.

Before we can continue we also need to alter the **Makefile** in the 
**src** directory, so that our **code.cc** file is included in the 
build. To do so, open an editor on the Makefile and change the line::

    CODEOBJS = test.o

to::

    CODEOBJS = test.o code.o
    

Test if the code builds. As we have not coupled our algorithm to the 
interface we (we have not even defined an interface) we do not have 
any new functionality. In the legacy interface directory (not the 
**src** directory) do:

.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works, if the test fails for any reason please check that the C++
code is correct and that ``worker_code`` exists in your directory.

Path 1
------

Defining the legacy interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will first define a legacy interface so that we can call the 
**find_nearest_neighbors** function from python. AMUSE can interact 
with 2 classes of functions:


1. A function with all scalar input and output variables. All variables
   are simple, non-composite variables (like int or double).
   For example:
   
   .. code-block:: c++
    
        int example1(double input, double *output)
        {
            *output = input;
            return 0;
        }
        
2. A function with all vector (or list) input and output variables and
   a length variable. The return value is a scalar value.
   For example:
   
   .. code-block:: c++
    
        int example2(double *input, double *output, int n)
        {
            for(int i = 0; i < n; i++)
            {
                output[i] = input[i];
            }
            return 0;
        }

If you have functions that don't follow this pattern you need to 
define a convert function in C++ that provides an interface 
following one of the two patterns supported by AMUSE.

In our case the **find_nearest_neighbors** complies to pattern 2 and
we do not have to write any code in C++ to convert the function
to a compliant interface. We only have to specify the function in 
python. We do so by adding a **find_nearest_neighbors** method to
the **NearestNeighborInterface** class in the **interface.py** file.
Open and editor on the interfaces.py file and add the following method
to the NearestNeighborInterface class:

.. code-block:: python

    class NearestNeighborInterface(InCodeComponentImplementation):    
        #...
        
        @legacy_function
        def find_nearest_neighbors():
            function = LegacyFunctionSpecification()  
            function.must_handle_array = True 
            function.addParameter(
                'npoints', 
                dtype='int32', 
                direction=function.LENGTH)
            function.addParameter(
                'x', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'y', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'z', 
                dtype='float64', 
                direction=function.IN)
            function.addParameter(
                'n0', 
                dtype='int32', 
                direction=function.OUT)
            function.addParameter(\
                'n1', 
                dtype='int32', 
                direction=function.OUT)
            function.addParameter(
                'n2', 
                dtype='int32', 
                direction=function.OUT)
            function.result_type = 'int32'
            return function
            
In the *find_nearest_neighbors* method we specify every parameter of 
the C++ function and the result type. For each parameter we need 
to define a name, data type and whether we will input, output (or 
input and output) data using this parameter. AMUSE knows only a 
limited amount of data types for parameters: float64, float32, int32, \
string, bool and int64. We also have a special parameter, with LENGTH as 
direction. This parameter is needed for all functions that follow 
pattern 2, it will be filled with the length of the input arrays. We 
also must specify that the function follows pattern 2 by setting 
```function.must_handle_array = True```.

Save the file and recompile the code.

.. code-block:: bash

    > make clean
    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

It works! But, how do we know the find_nearest_neighbors method really works?
Let's write a test and find out. Open an editor on the *test_nearestneighbor.py*
file and add the following method to the **NearestNeighborInterfaceTests** class:

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        x = [0.0, 1.0, 4.0, 7.5]
        y = [0.0] * len(x)
        z = [0.0] * len(x)
        
        n0, n1, n2, error = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(error[0], 0)
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()


This test calls the *find_nearest_neighbors* method with 4 positions
and checks if the nearest neighbors are determined correctly.
Let's run the test, and see if everything is working:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 test in 0.491s

    OK

We now have a simple interface that works, but we have to do our own 
indexing after the call and we could send data of any unit to the
method, also we have to do our own error checking after the method. Let's
define a object oriented interface to solve these problems

Defining the Object Oriented Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The object oriented interface sits on top of the legacy interface.
It decorates this interface with sets, unit handling, state engine
and more. We start creating the object oriented interface by inheriting
from InCodeComponentImplementation and writing the __init__ function. The build script
has added this class to the interface.py file for us. Open an editor
on interface.py and make sure this code is in the file (at the end
of the file):

.. code-block:: python

    
    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)

Configuring the handlers
~~~~~~~~~~~~~~~~~~~~~~~~

We configure the object oriented interface by implementing several 
methods. The object oriented interface is implement by several 
"handlers". Each handler provides support for a specific aspect of 
the interface. AMUSE defines a handler for the unit conversion, a 
handler for the interfacing with sets of particles, a handler to 
ensure the methods are called in the right order, etc. Each handler 
is very generic and needs to be configured before use. The handler 
are configured using the "Visitor" pattern. The following 
pseudo-code shows how the handlers are configured

.. code-block:: python

    class InCodeComponentImplementation(object):
        #...
        
        def configure_handlers(self):
            #...
            for handler in self.get_all_handlers():
                handler.configure(self)
        
        def define_converter(self, handler):
            """ configure the units converter handler """
            
            handler.set_nbody_converter(...)
            
        def define_particle_sets(self, handler):
            """ configure sets of particles """
            
            handler.define_incode_particle_set(...)
            handler.set_getter(...)
            
    class HandleConvertUnits(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_converter(self)
            
    class HandleParticles(AbstractHandler):
        #...
        
        def configure(self, interface):
            interface.define_particle_sets(self)
            
            
Configuration of the handlers is optional, we only have to define
those handler that we need in our interface. In our example we need
to configure the "HandleMethodsWithUnits" handler (to define units and
error handling) and the "HandleParticles" to define a particle set.

Defining methods with units
+++++++++++++++++++++++++++

We first want to add units and error handling to the 
**find_nearest_neighbors**. We do this by creating a 
**define_methods** function on the **NearestNeighbor** class. Open 
an editor on *interface.py* and add this method to the class:

.. code-block:: python

    def define_methods(self, handler):
        
        handler.add_method(
            "find_nearest_neighbors",
            (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
            ),
            (
                handler.INDEX,
                handler.INDEX,
                handler.INDEX,
                handler.ERROR_CODE
            )
        )
        
The **add_method** call expects the name of the function in the
legacy interface as it's first, next it expects a list of the units
of the input parameters and a list of the units of the output parameters.
The return value of a function is always the last item in the list
of output parameters. We specify a **generic_unit_system.length** unit
for the x, y and z parameters. The output parameters are indices and 
an errorcode. The errorcode will be handled by the AMUSE interface (0
means success and < 0 means raise an exception).

Let's write a test to see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test3(self):
        instance = NearestNeighbor()
        x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        y = [0.0] * len(x) | generic_unit_system.length
        z = [0.0] * len(x) | generic_unit_system.length
        
        n0, n1, n2 = instance.find_nearest_neighbors(x,y,z)
        
        self.assertEquals(n0, [2,1,2,3])
        self.assertEquals(n1, [3,3,4,2])
        self.assertEquals(n2, [4,4,1,1])
        
        instance.stop()

.. note::
    This test looks a lot like test2, but we now have to
    define a unit and we do not need to handle the errorcode.
    
Now build and test the code:

.. code-block:: bash
    
    > make clean; make all
    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.650s

    OK
    
.. note::
    
    Although we only edited python code we still need to run make. The
    code will check if the "worker_code" executable is up to date on 
    every run. It cannot detect if the update broke the code but it
    will still demand that the code is rebuilt.
    

Defining the particle set
+++++++++++++++++++++++++

Particle sets in AMUSE can be handled by python (we call these 
"inmemory") and by the legacy code (we call these "incode"). In our 
case the code does not handle the particles and we need to 
configure the particles handler to manage an inmemory particle set. 
Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    def define_particle_sets(self, object):
        object.define_inmemory_set('particles')

That's all we now have a "particles" attribute on the class and 
we can add, remove, delete particles from this set. But we are
still missing a connection between the particles and the nearest 
neighbors. AMUSE provides no handler for this, instead, we 
will write a method to run the find_nearest_neighbors function and
set the indices on the particles set.

Open an editor on *interface.py* and add this method to the 
**NearestNeighbor** class:

.. code-block:: python

    
    def run(self):
        indices0, indices1, indices2 = self.find_nearest_neighbors(
            self.particles.x,
            self.particles.y,
            self.particles.z
        )
        self.particles.neighbor0 = list(self.particles[indices0])
        self.particles.neighbor1 = list(self.particles[indices1])
        self.particles.neighbor2 = list(self.particles[indices2])
        
This function gets the "x", "y" and "z" attributes from the particles
set and sends these to the "find_nearest_neighbors" method. This methods
returns 3 lists of indices and we need to find the particles with
these indices. 

.. note::
    
    Particle sets have no given sequence, deletion and addition of 
    particles will change the order of the particles in the set. It
    is therefor never a good idea to use the index of the particle in the
    set as a reference to that particle. However, in the "run" method
    we "own" the particle set, it cannot change between the find_nearest_neighbor
    call and the moment we find the particles in the set by index (using
    self.particles[indices0]), and in this case it is save to use
    index as a valid reference.

Let's write a test and see if it works, open an editor on the 
test_nearestneighbor.py class and add this method:

.. code-block:: python

    def test4(self):
        instance = NearestNeighbor()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()

Now, make and run the tests:

.. code-block:: bash

    > make clean; make all
    > nosetests
    ....
    ----------------------------------------------------------------------
    Ran 4 tests in 0.797s

    OK
 
We are done, we have defined an object oriented interface on the
legacy interface. Only, if we look at our tests, the code seems
to be more rather than less complex. But, remember we now
have units and we are compatible with other parts of amuse. And
we can make more complex scripts easier.

Let's make a plummer model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash

    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png



Path 2
------

Defining the legacy interface
-----------------------------
We define our code interface so that a user can add, update and 
delete particles, start the nearest neighbors finding
algorithm and retrieve the ids of the nearest neighbors.

To define the interface, open interface.py with your favorite
editor and replace the contents of this file with:

.. literalinclude:: nearestneighbor/nn1.py

We can generate a stub from the interface code with::

    > $AMUSE_DIR/build.py --type=c --mode=stub interface.py NearestNeighborInterface -o interface.cc

The generated **interface.cc** replaces the original file generated in
the previous section.     

The code builds, but does not have any functionality yet::

    > make clean
    > make all
    
.. note::
    
    Compiling the interface code will result in a lot of warnings about
    unused dummy arguments. These warnings can be safely ignored for now.
    
The tests are broken (the echo_int function has been removed):

.. code-block:: bash

    > nosetests
    E
    ======================================================================
    ERROR: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "../src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "nearestneighbor/test_nearestneighbor.py", line 11, in test1
        result,error = instance.echo_int(12)
    AttributeError: 'NearestNeighborInterface' object has no attribute 'echo_int'

    ----------------------------------------------------------------------
    Ran 1 test in 0.315s

    FAILED (errors=1)

Let's create a working test by calling the new_particle method, open an editor
on the test_nearestneighbor.py file and replace the test1 method with:

.. code-block:: python

    def test1(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()

As this is python code we do not need to rebuild the code, instead
we can run the tests right after saving the code. Unfortunately, when
we run the test, it still fails.

.. code-block:: bash

    > nosetests
    F
    ======================================================================
    FAIL: test1 (nearestneighbor.test_nearestneighbor.NearestNeighborInterfaceTests)
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "/src/amuse/test/amusetest.py", line 146, in run
        testMethod()
      File "/nearestneighbor/test_nearestneighbor.py", line 13, in test1
        self.assertEquals(result, 1)
      File "/src/amuse/test/amusetest.py", line 62, in failUnlessEqual
        self._raise_exceptions_if_any(failures, first, second, '{0} != {1}', msg)
      File "/src/amuse/test/amusetest.py", line 49, in _raise_exceptions_if_any
        raise self.failureException(msg or err_fmt_string.format(first, second, *args))
    AssertionError: 0 != 1
    -------------------- >> begin captured logging << --------------------
    legacy: INFO: start call 'NearestNeighborInterface.new_particle'
    legacy: INFO: end call 'NearestNeighborInterface.new_particle'
    --------------------- >> end captured logging << ---------------------

    ----------------------------------------------------------------------
    Ran 1 test in 0.319s

When you look closely at the output of the test you see that the result from the
method is 0 and not the expected 1. We need to edit the c code to make this 
test work. Open an editor on interface.cc and go to the *new_particle* function.

.. code-block:: c++

    int new_particle(int * index_of_the_particle, double x, double y, 
        double z)
    {
        *index_of_the_particle = 1;
        return 0;
    }

.. note ::
    
    In AMUSE all interface functions return an errorcode. Any other return
    values must be passed through the arguments of the functions.

We need to rebuild the code, and after building run the tests.

.. code-block:: bash

    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.427s

    OK

The are tests work again! Only, we do not have any real working legacy code.

Filling the stubs
-----------------
The implementation of the algorithm does not match the interface
we defined and created. We need to write some glue code to
connect the code with the interface. To do so we fill in the
stubs generated earlier.

Open the **interface.cc** file in your favorite editor and
change its contents to:

.. literalinclude:: nearestneighbor/interface1.cc
    :language: c++

Test if the code builds and try it out. In the legacy interface
directory do:
    
.. code-block:: bash

    > make clean
    > make all
    > nosetests
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.311s

    OK
    
Let's check some more functionality by adding another test

.. code-block:: python

    def test2(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 2)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, -1)
        error = instance.delete_particle(1)
        self.assertEquals(error, 0)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()
        

The tests should succeed:

.. code-block:: bash

    > nosetests
    ..
    ----------------------------------------------------------------------
    Ran 2 tests in 0.448s

    OK


We now have done everything in Step 0 *Legacy Interface*. We have
a legacy code and can access it in our python script. But, our
interface is not very friendly to work with. We have to think about
errorcodes and we have not information about units. To make our
interface easier to works with we start defining methods, properties
and parameters.

Defining methods
----------------
The object oriented interface is also defined in the **interface.py**.
So, we continue by opening an editor on this file. We will be 
writing methods for the **NearestNeighbor** class, in your editor
seek this code (at the end of the file)::

    class NearestNeighbor(CodeInterface):

        def __init__(self, **options):
            CodeInterface.__init__(self,  NearestNeighborInterface(), **options)
            
We will start by defining methods, we will do this by implementing
the **define_methods** function, like so::

    class NearestNeighbor(CodeInterface):

        def __init__(self, **options):
            CodeInterface.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            
            builder.add_method(
                "new_particle", 
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.INDEX, builder.ERROR_CODE)
            )
            
            builder.add_method(
                "delete_particle", 
                (builder.INDEX,),
                (builder.ERROR_CODE)
            )
            
            builder.add_method(
                "get_state", 
                (builder.INDEX,),
                (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, builder.ERROR_CODE),
                public_name = "get_position"
            )
            
            builder.add_method(
                "set_state", 
                (builder.INDEX, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
                (builder.ERROR_CODE),
                public_name = "set_position"
            )
            
            builder.add_method(
                "run", 
                (),
                (builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_close_neighbors", 
                (builder.INDEX,),
                (builder.LINK('particles'), builder.LINK('particles'), builder.LINK('particles'), builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_nearest_neighbor", 
                (builder.INDEX,),
                (builder.LINK('particles'), generic_unit_system.length, builder.ERROR_CODE),
            )
            
            builder.add_method(
                "get_number_of_particles", 
                (),
                (builder.NO_UNIT, builder.ERROR_CODE),
            )

        
With this code, we define the methods and specify how to interpret the
arguments and return values. We get a special object (the **builder**
object) that provides us with the **add_method** function to be able
to this. The definition of the **add_method** function is as follows::
    
    add_method(
        name of the original function in the legacy interface,
        list of arguments (unit or type),
        list of return values,
        public_name = name for the user of the class (optional)
    )
    
For every argument or return value we can specify if it has a unit or
if it is special. The special arguments are:

=========================== ===================================
definition                  description
=========================== ===================================
builder.ERROR_CODE          the value is interpreted as an errorcode,
                            zero means no error for all other values and
                            Exception will be raise, only valid for one return value.
builder.NO_UNIT             the value has not unit (for example for the
                            number of items in a list)
builder.INDEX               the value is interpreted as an index for
                            object identifiers
builder.LINK('particles')   the value is interpreted as a link to
                            objects in the set with the given name
=========================== ===================================

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test:

.. code-block:: python

    def test3(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(2)
        instance.commit_parameters()
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length, 
            3.0 | generic_unit_system.length
        )
        self.assertEquals(result, 1)
        result = instance.new_particle(
            1.0 | generic_unit_system.length, 
            1.0 | generic_unit_system.length, 
            2.0 | generic_unit_system.length
        )
        self.assertEquals(result, 2)
        x,y,z = instance.get_position(1)
        self.assertEquals(1.0 | generic_unit_system.length, x)
        self.assertEquals(2.0 | generic_unit_system.length, y)
        self.assertEquals(3.0 | generic_unit_system.length, z)
        instance.stop()

And run the tests:

.. code-block:: bash

    > nosetests
    ...
    ----------------------------------------------------------------------
    Ran 3 tests in 0.664s

    OK

As you can see our script is now a little simpler and we support units.
We do not have to think about the errorcodes in this script, AMUSE will
interpret the errorcodes and raise the right exceptions if needed. The 
units are also automatically converted to the right units for the 
code. But the script is still not very easy and we have to manage
all the ids we get from the code. To make our code even easier to
handle we will continue by defining a **set**.

.. note::

    We skip defining parameters and properties, we will come back to 
    this in the next tutorial.

Defining a set
--------------
We have made our interface a little easier but we still have to
do a some management work in our script. We would like to work
with objects and adding or removing these objects from the code.
AMUSE supports this by defining **sets**. Each set is capable of
storing specific attributes of the objects in the set. Our code
is capable of storing the x, y and z position of an object. An object
in AMUSE is called a *Particle* and the sets that contain these
particles are called *ParticleSets* or shorter *Particles*.

We define our particle set by implementing a **define_particle_sets** 
function on our **NearestNeighbor** class like so::

    class NearestNeighbor(InCodeComponentImplementation):

        def __init__(self, **options):
            InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)
            
        
        def define_methods(self, builder):
            ...
            
        def define_particle_sets(self, builder):
            builder.define_set('particles', 'index_of_the_particle')
            builder.set_new('particles', 'new_particle')
            builder.set_delete('particles', 'delete_particle')
            builder.add_setter('particles', 'set_position')
            builder.add_getter('particles', 'get_position')
            builder.add_getter('particles', 'get_close_neighbors', names=('neighbor0', 'neighbor1', 'neighbor2') )


That's all, we now have defined a set called **particles**. Again, we 
get a builder object to use in defining our set. All methods have
the name of the set as their first argument, this name can be any
name you want, but in AMUSE most codes provide a set called 
**particles**. For the **add_setter**, **add_getter**, **set_new** 
and **set_delete** functions, the second argument is the name of
the method we defined in the previous step. Finally you can set
the name of the attribute in the particles set with the **names**
argument. This is optional for legacy functions, if not given the names 
of the attributes will be derived from the names of the arguments in
the original calls. For example, the **get_position** call we specified
earlier has parameter name **x**, **y** and **z**, these names are 
also used in the particles set.

Test if the code builds and try it out. In the legacy interface
directory do::
    
    > make clean
    > make all

Let's add another test::
    
    def test4(self):
        instance = NearestNeighbor()
        instance.set_maximum_number_of_particles(100)
        instance.commit_parameters()
        
        particles = datamodel.Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length
        
        instance.particles.add_particles(particles)
        instance.run()
        
        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])
        
        instance.stop()
    
    
The support for a particle set means we can now also interact
with other parts of AMUSE. Let's make a plummer
model and find the nearest neighbors in this model.

First make a file with the following contents, let's call this file
**plummer2.py**:

.. literalinclude:: nearestneighbor/plummer2.py
    :language: python
    
We can run this file with python::

.. code-block:: bash
    
    $AMUSE_DIR/amuse.sh plummer2.py
    
It will create an **output.txt** file and we can show this file
with gnuplot.

.. code-block:: gnuplot

    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    gnuplot> #we can zoom into the center
    gnuplot> set xr[-0.5:0.5]
    gnuplot> set yr[-0.5:0.5]
    gnuplot> set zr[-0.5:0.5]
    gnuplot> splot 'output.txt' using 1:2:3:4:5:6 with vectors nohead, 'output.txt' using 1:2:3
    

.. image:: nearestneighbor/plummer1.png

.. _tutorials-label:

Tutorials
=========

.. toctree::
   :maxdepth: 2
   
   units
   most_massive_from_salpeter
   particle_sets
   c_code
   fortran_code
   gravitational_dynamics_code
   legacy_code
   plot
   grid_boundary
   documenting
   

==========================
Working with Particle Sets
==========================

.. image:: sets.png
   :width: 2cm

There are Sets, Subsets, Particles, Children and Parents...
===========================================================
The fundamental data structures used in AMUSE are particle sets. Based 
on attributes of the elements in the sets (particles), selections can
be made using the selection method which return subsets. These subsets
are views or scopes on the set and do not hold values of their own.

It is also possible to *add* structure to the set by defining 
parent-child relationships between **particles in a set**. 
These structures exist only in the set and are a property
of it and have no meaning with respect to the particles
outside the set.


Selecting stars in a plummer model
==================================
In this tutorial we generate a set of stars, a plummer model, and 
two subsets of stars. The first subset contains the stars within
a certain radius from the center of mass of the plummer system and
the second subset contains all the other stars.

The plummer module generates a particle set of the form datamodel.Stars(*N*).

Imports
-------
We start with some AMUSE imports in a python shell:
The core module contains the set objects, and the units and nbody_system
are needed for the correct assignment of properties to
the particles, which **need to be quantities**. 

.. code-block:: python

    >>> from amuse import datamodel
    >>> from amuse.units import nbody_system
    >>> from amuse.units import units
    >>> from amuse.ic.plummer import new_plummer_sphere
    >>> import numpy


Model
-----
.. note ::

    A quick way to look at the sets we are going to make
    is by using gnuplot. If you have gnuplot, you can install
    the **gnuplot-py** package to control gnuplot directly
    from your script.
    
    To install **gnuplot-py**, open a shell and do::
        
        easy_install gnuplot-py
        
        
Let's generate a plummer model consisting of 1000 stars, 

.. code-block:: python

    >>> convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
    >>> plummer = new_plummer_sphere(1000, convert_nbody)

We can work on the new **plummer** particle set, but we want to keep this
set unchanged for now. So, we copy all data to a working set:

.. code-block:: python
   
    >>> stars = plummer.copy()
    
.. note ::

    To look at the stars in gnuplot do::
    
        >>> plotter = Gnuplot.Gnuplot()
        >>> plotter.splot(plummer.position.value_in(units.parsec))
    
    .. image:: sets/plummer1.png
    
    
    
    
Selection
---------
At this stage we select the subsets based on the distance of the
individual stars with respect to the center of mass, being 
**1 parsec** in this example.

We need the center of mass, of course:

.. code-block:: python

    >>> center_of_mass = stars.center_of_mass()   

and the selection of the sets:
 
.. code-block:: python
 
    >>> innersphere = stars.select(lambda r: (center_of_mass-r).length()<1.0 | units.parsec,["position"])  
    >>> outersphere = stars.select(lambda r: (center_of_mass-r).length()>=1.0 | units.parsec,["position"])  

.. note ::

    To look at the stars in gnuplot do::
    
        >>> plotter = Gnuplot.Gnuplot()
        >>> plotter.splot(outersphere.position.value_in(units.parsec), innersphere.position.value_in(units.parsec),)
    
    .. image:: sets/plummer2.png
    
    

We can achieve the same result in another way by using the fact that outersphere is the difference of the innersphere set 
and the stars set:

.. code-block:: python

    >>> outersphere_alt = stars.difference(innersphere)

or using the particle subtraction '-' operator:

.. code-block:: python

    >>> outersphere_alt2 = stars - innersphere

The selections are all subsets as we can verify:

.. code-block:: python

    >>> outersphere #doctest:+ELLIPSIS 
    <amuse.datamodel.particles.ParticlesSubset object at ...>

Set operations
--------------

The result should be the same, but we'll check:

.. code-block:: python

    >>> hopefully_empty = outersphere.difference(outersphere_alt)
    >>> hopefully_empty.is_empty()
    True
    >>> len(outersphere - outersphere_alt2)==0
    True

From our selection criteria we would expect to have selected all stars, to check this we could do something like this:

.. code-block:: python

    >>> len(innersphere)+len(outersphere) == 1000
    True
    >>> len(innersphere)+len(outersphere_alt) == 1000
    True
    >>> len(innersphere)+len(outersphere_alt2) == 1000
    True

The union of the innersphere and outersphere set should give the stars set, we can check:

.. code-block:: python

    >>> like_stars = innersphere.union(outersphere)
    >>> stars.difference(like_stars).is_empty()
    True
    >>> (innersphere + outersphere_alt2 - stars).is_empty()
    True

Iteration
---------

We can iterate over sets and we will put that to use here to check whether we really selected the right stars in the outersphere:

.. code-block:: python

    >>> should_not_be_there_stars = 0
    >>> for star in outersphere:
    ...     if (center_of_mass-star.position).length()<1.0|units.parsec:
    ...         should_not_be_there_stars += 1
    >>> should_not_be_there_stars
    0

Indexation
----------

Not very set like, but can-do.


Using codes
===========  

Channels
--------

.. image:: channels.png
   :width: 2cm
   :align: center

Imagine we evaluate stars in some MPI bound legacy code, and want to know what happened to the stars in our subsets. Each query for attributes in de 
code set invokes one MPI call, which is inefficient if we have many queries. Copying the entire set to the model, however, costs only one MPI call too
and we can query from the model set at will without the MPI overhead. Copying the data from one existing set to another set *or* subset 
can be done via channels.

First we define a simple dummy legacy code object, with some typical gd interface methods:

.. code-block:: python

    >>> class DummyLegacyCode(object):
    ...     def __init__(self):
    ...         self.particles = datamodel.Particles()
    ...     def add_particles(self, new_particles):
    ...         self.particles.add_particles(new_particles)
    ...     def update_particles(self, particles):
    ...         self.particles.copy_values_of_attributes_to(['x', 'y', 'z', 'mass'], particles)
    ...         particles = self.particles.copy()           
    ...     def evolve(self):
    ...         self.particles.position *= 1.1

We instantiate the code and use it to evolve (expand) our plummer model. We add a channel to our innersphere subset to track the changes:

.. code-block:: python

    >>> code = DummyLegacyCode()
    >>> channel_to_innersphere = code.particles.new_channel_to(innersphere)
    >>> code.add_particles(stars)
    >>> code.evolve()
    >>> r_inner_init = innersphere.position
    >>> r_outer_init = outersphere.position
    >>> channel_to_innersphere.copy()
    >>> r_inner_fini = innersphere.position
    >>> r_outer_fini = outersphere.position

Checking the changes by looking at the positions (all changed after the evolve), we will see that only the innersphere particles are updated:

.. code-block:: python

    >>> numpy.all(r_inner_init[0]==r_inner_fini[0])
    False
    >>> numpy.all(r_outer_init[0]==r_outer_fini[0])
    True

If we want to update all particles in our model, we can use:

.. code-block:: python

    >>> code.update_particles(stars)
    
and check:

.. code-block:: python

    >>> r_outer_fini = outersphere.position
    >>> numpy.all(r_outer_init[0] == r_outer_fini[0])
    False

Particle Hierarchies
====================

.. image:: children.png
   :width: 2cm

Let us suppose that the zero-th star, stars[0] has a child and a grandchild star who do not belong to the plummer model:

.. code-block:: python

    >>> child_star = datamodel.Particle()
    >>> grandchild_star = datamodel.Particle()
    >>> child_star.mass = 0.001|units.MSun
    >>> child_star.position = [0,0,0]|units.AU
    >>> child_star.velocity = [0,0,0]|units.AUd
    >>> grandchild_star.mass = 0.0001|units.MSun
    >>> grandchild_star.position = [0,0.1,0]|units.AU
    >>> grandchild_star.velocity = [0,0,0]|units.AUd

We can add them as child and grandchild etc. to the set of plummer stars. But first we have to add them to the set as regular stars:

.. code-block:: python

    >>> child_star_in_set = stars.add_particle(child_star)
    >>> grandchild_star_in_set = stars.add_particle(grandchild_star)

Now we can define the hierarchy:

.. code-block:: python

    >>> stars[0].add_child(child_star_in_set)
    >>> child_star_in_set.add_child(grandchild_star_in_set)

The descendants of star 0 form a subset:

.. code-block:: python

    >>> stars[0].descendents() #doctest:+ELLIPSIS
    <amuse.datamodel.ParticlesSubset object at ...>
    >>> stars[0].children().mass.value_in(units.MSun)
    array([ 0.001])
    >>> stars[0].descendents().mass
    quantity<[1.98892e+27, 1.98892e+26] kg>


Methods to retreive physical properties of the particles set
============================================================

Particle sets have a number of functions for calculating physical
properties that apply to the sets. Although some of these might be implemented
in the legacy codes as well, using them via particle sets guarantees
the applicability to *all* particles when multiple legacy codes
are used. Furthermore, the particle set functions provide a uniform
way of doing the calculations. Speed might be the downside.

.. automodule:: amuse.datamodel.particle_attributes
       
.. autofunction:: center_of_mass 
    :noindex:
.. autofunction:: center_of_mass_velocity
    :noindex:
.. autofunction:: kinetic_energy
    :noindex:
.. autofunction:: potential_energy
    :noindex: 
.. autofunction:: particle_specific_kinetic_energy
    :noindex:
.. autofunction:: particle_potential
    :noindex:

====================================
Adding a Gravitational Dynamics Code
====================================

In the previous tutorial we have seen how to create an AMUSE interface to
a C++ or a Fortran code from scratch. In this tutorial, we will expand on this
with a number of additional features. We will do this through the implementation
of a simple gravitational dynamics code.

As this tutorial deals only with features of the Python interface, we restrict our
legacy code example to C++. The concepts described can just as easily be applied to
an interface to a Fortran code. 

We will also assume that we are working with the same environment settings as in
the previous example. Thus, we can create the initial directory structure of our new code 'SimpleGrav' with the following terminal command:

.. code-block:: bash

    > $AMUSE_DIR/build.py --type=c --mode=dir SimpleGrav


The Legacy Code
---------------
For this example we will use a very simple forward Eulerian integration routine as our
legacy code. We simply have a function that takes the initial condition of a set of particles (in the form of seven dynamic arrays), an integer containing the number of particles, and double precision scalars containing the time step, smoothing length, and gravitational constant. It outputs the new particle positions and velocities by updating the input arrays.

.. literalinclude:: simplegrav/SimpleGrav.cc
    :language: C++


.. note::

    Just like in the previous tutorial, the example algorithm is by no means
    particularly good. In fact, forward Eulerian integration might be the worst 
    method of integrating gravitational systems due to the tendency of the energy
    to diverge. Again, we use it as a readable example.

Simply put this script in the **src** directory as **code.cc**, and again change
the Makefile line::

    CODEOBJS = test.o

to::

    CODEOBJS = test.o code.o


Interface Templates
-------------------
In order to promote uniformity among interfaces, and to make it easier to create
an interface, a number of types of codes have interface templates. These classes
have a number of functions common to that type of code already defined. They can be
found in the folder **src/amuse/community/interface** of the main AMUSE github
repository. 

Gravitational dynamics is one type of code having a template interface, defined
in **gd.py**. We can use this template by inheriting from it:

.. code-block:: python

    from amuse.community.interface.gd import GravitationalDynamics
    from amuse.community.interface.gd import GravitationalDynamicsInterface

    class SimpleGravInterface(CodeInterface, GravitationalDynamicsInterface):
        ...

    class SimpleGrav(GravitationalDynamics):
        ...

The legacy interface now contains all legacy functions defined in **gd.py**.
This includes adding and deleting particles, getting and setting particle properties
(mass, position, velocity, acceleration, and radius), getting and setting 
parameters (smoothing length and begin time), and a host of other functions. Our
definition of the ``SimpleGravInterface`` is then simply:

.. code-block:: python

    class SimpleGravInterface(CodeInterface, GravitationalDynamicsInterface):

        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="SimpleGrav_worker",
                **keyword_arguments)

Note that when compiling an interface inheriting from
``GravitationalDynamicsInterface``, all legacy functions defined there
(and in its parent class, ``CommonCodeInterface``) must be
defined in the C++ interface. For now, you can simply let them return an error 
value (-1 or -2).

The object oriented interface similarly contains definitions for methods, properties,
particle sets, parameters, and states (the latter we will discuss later). Its definition
is simply:

.. code-block:: python

    class SimpleGrav(GravitationalDynamics):

        def __init__(self, convert_nbody = None, **keyword_arguments):
            
            GravitationalDynamics.__init__(self,
                                           SimpleGravInterface(**keyword_arguments),
                                           convert_nbody,
                                           **keyword_arguments)


Note the ``convert_nbody`` keyword that appears in the initializer. This is used
to handle the conversion of units between the legacy code and the object-oriented
interface. The gravitational dynamics interface template assumes that the legacy
codes work in dimensionless N-body units, but by passing a ``convert_nbody``
instance the unit conversions can be handled automatically.

The only part left to write for a functional version is the C++ interface. As mentioned before,
this must contain the corresponding functions of all legacy functions. Also note
that ``GravitationalDynamicsInterface`` inherits from yet another interface,
``CommonCodeInterface``, which contains four more legacy functions to be defined.
These are very general functions, mostly concerned with code logistics. Their 
implementation, and whether they are even necessary, depends on the legacy code. If
they are not needed they can simply return 0 (in the case of our legacy code, only
``cleanup_code`` has actual functionality). They are often called automatically,
when needed, through the state model (discussed below). The four functions from
``CommonCode``, plus three similar functions from
``GravitationalDynamics``, are described in the table below:

=========================== ===================================
definition                  description
=========================== ===================================
initialize_code             any code that must be executed when the code is
                            initialized, but before parameters are set, particles
                            are added, a grid is initialized, etc. 
cleanup_code                any code that must be executed when the code is
                            stopped. This typically frees up memory allocated
                            by the code.
commit_parameters           any initialization that must be done that is dependent
                            on parameters, for example the allocation of memory
                            for a grid.
recommit_parameters         similar to commit_parameters, but after 
                            commit_parameters has been called once. This could
                            be a redefinition of a grid.
commit_particles            any initialization that must be done that is dependent
                            on particles, for example the construction of a tree
                            of some kind from a set of particles.
recommit_particles          similar to commit_particles, but after 
                            commit_particles has been called once. This could
                            the addition or removal of particles from a tree.
synchronize_model           evolve all particles to the same time, if they
                            are at different times.
=========================== ===================================

The following code contains definitions for all legacy functions, although some
non-essential functions only return error values:


.. literalinclude:: simplegrav/interface_1.cc
    :language: C++



Properties & Parameters
-----------------------
A property of an AMUSE code is a read-only quantity of the code that contains
information about the (often physical) state of the code. Examples include the
current model time and the potential and kinetic energies of the total system.

A parameter of an AMUSE code is a quantity that is used in the internal workings
of the code. Examples include the time step in most types of evolution codes and
the metallicity in stellar evolution codes. These are often, but not always,
writable. 

Both are defined in the object-oriented interface, in the ``define_properties``
and ``define_parameters`` functions. 

Properties are added as follows:

.. code-block:: python

    def define_properties(self, handler):
        handler.add_property("get_property", public_name="property_name")

Here, ``get_property`` is a legacy interface function that returns a single value.
``public_name`` is an optional argument that defines what the property is called;
if it is not given, the name of the parameter of the legacy function is used.

Properties are accessed directly as properties of the code object:

.. code-block:: python

    gravity = SimpleGrav()
    the_current_time = gravity.model_time

Parameters are added as follows:

.. code-block:: python

   def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_parameter",
            "set_parameter",
            "parameter",
            "what is this parameter?",
            default_value = 1. | some_unit,
            is_vector = BOOL,
            must_set_before_get = BOOL,
        )

Here, ``get_parameter`` and ``set_parameter`` are two legacy interface functions
that respectively have one output and one input. ``parameter`` is the name of the
parameter in the code. The fourth parameter is a documentation string describing
what the parameter signifies. Finally, there are three optional parameters. 
``default_value`` is the value the parameter is set to by default.
If ``is_vector == True``, the parameter is a vector of a fixed length. This
length is defined as the number of parameters of the getters and setters
(scalar parameters have only one parameter). 
If ``must_set_before_get == True``, the parameter must be set before it can be returned. 
In most cases, this can be False.


Note that there are methods to define other types of parameters:

=========================== ============================
function                    purpose
=========================== ============================
add_boolean_parameter       the parameter is a boolean,
                            instead of an integer/float.
add_caching_parameter       the parameter is only set
                            once ``commit_parameters`` is
                            called.
add_array_parameter         the parameter is an array;
                            in contrast with a normal
                            parameter with ``is_vector=True``,
                            the getters and setters have
                            only a single argument, which
                            must be able to handle arrays.
                            Between the setter and the name,
                            an additional function must be
                            passed that specifies the size
                            of the array.
=========================== ============================

``add_alias_parameter`` is also available to add an alias to an existing parameter.
It is of the following form: ``add_alias_parameter("alias_parameter", "original_parameter", 
"parameter documentation")``.


Parameters are accessed through the ``parameters`` property of the code object:

.. code-block:: python

    gravity = SimpleGrav()
    param = gravity.parameters.parameter




.. note::

    With properties, particles, and parameters, we have three ways of communicating
    data between the legacy code and AMUSE. With complex codes that are not
    necessarily similar to other codes in AMUSE, it might be nontrivial in which
    of these manners data is to be communicated. As a rule of thumb, any quantity 
    relating to individual particles should be communicated as particle properties, 
    even if it is conceptually closer to a parameter. To demonstrate the thin line
    between these, if a stellar evolution code has a single metallicity for all
    stars, it will be parameter, whereas if different stars can have different
    metallicities, it will be a particle property.

    Parameters and properties, on the other hand, both apply to the system as a whole.
    Of these, parameters should not be altered by the legacy code itself, whereas
    properties are read-only. 

In our SimpleGrav example we have introduced a number of parameters that the template
interface does not have. For one, we want the time step to be a settable parameter
(note that the definition in **gd.py** has ``None`` as a setter function). Additionally,
we want to be able to set a smoothing length, and for didactic reasons we include an
entirely new parameter, ``gravity_strength``, which functions as a scaling of the
magnitude of the gravitational force. 

We have to define new legacy functions to set the time step and get and set the gravity
strength (the others are defined in **gd.py**). The legacy interface then looks as follows:

.. code-block:: python

    class SimpleGravInterface(CodeInterface,
                       GravitationalDynamicsInterface):

        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="SimpleGrav_worker",
                **keyword_arguments)

        @legacy_function
        def get_grav_fac ():
            function = LegacyFunctionSpecification()
            function.addParameter('grav_fac', dtype='float64', direction=function.OUT)
            function.result_type = 'int32'
            return function

        @legacy_function
        def set_grav_fac ():
            function = LegacyFunctionSpecification()
            function.addParameter('grav_fac', dtype='float64', direction=function.IN)
            function.result_type = 'int32'
            return function

        @legacy_function
        def set_time_step ():
            function = LegacyFunctionSpecification()
            function.addParameter('time_step', dtype='float64', direction=function.IN)
            function.result_type = 'int32'
            return function

In the object-oriented interface, we have to overload the ``define_methods`` and the 
``define_parameters`` functions:

.. code-block:: python

    def define_methods (self, handler):

        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            "get_grav_fac",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE)
        )

        handler.add_method(
            "set_grav_fac",
            (handler.NO_UNIT),
            (handler.ERROR_CODE)
        )

        handler.add_method(
            "set_time_step",
            (nbody_system.time),
            (handler.ERROR_CODE)
        )

        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE)
        )

        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length),
            (handler.ERROR_CODE)
        )


    def define_parameters (self, handler):

        GravitationalDynamics.define_parameters(self, handler)

        handler.add_method_parameter(
            "get_grav_fac",
            "set_grav_fac",
            "gravity_strength",
            "constant factor by which G is multiplied",
            default_value = 1.
        )

        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "time_step",
            "constant integrator timestep",
            default_value = 0.01 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length


We also call the equivalent functions of **gd.py** so we still inherit
from the template interface. 




State Model
-----------
Many legacy codes have pieces of logistical code that must be executed between
physically relevant pieces of code. For example, in a gravitational tree code,
the tree must be constructed after all particles have been added, but before
the system is evolved, and not necessarily between evolve calls. In order to automate this, a state model can be defined
for the code. This defines what functions can be run in what state of the code,
but also how to transition between states, and what functions trigger that
transition. A particularly convenient function is that it allows a transition, 
and thus
the function associated with that, to be triggered automatically. This allows the
aforementioned tree code to automatically build its tree between the definition
of the particle set and the evolution of the system. 

The state of the code can be handled automatically by defining a state model that 
describes the states the code can be in as a graph with the nodes as the state and 
transitions mediated by interface calls. The state model of a code is defined in
the ``define_state`` interface function by calling the following methods on its 
handler argument: 

- ``set_initial_state(state)``  this defines the initial entry state, given as a string label (usually 'UNINITIALIZED').

- ``add_method(state, method_name)`` this defines that the given method is allowed in state ``state``. Again the state is a string, the string can be prepended with '!' to indicate a method is forbidden in the state. If state is not yet defined it will be added to the state model.

- ``add_transition(state1, state2, method_name, is_auto=True)`` this adds a transition between states state1 and state2 (and adds these states if not previously defined) which is triggered by a call to the interface function method_name. The is_auto argument determines whether this transition is allowed to  be triggered automatically. 

A method that is not mentioned in any add_method, is allowed in any state (and
doesn't trigger any transitions). If the state model detects that an interface call 
needs a state change it tries to hop to a state where the interface call is allowed
by calling transitions that are added with is_auto=True. (this is only possible if
they don't have non-default arguments)

The state model can be built up by the above methods. The state model of a code can be printed:

- ``interface.state_machine.to_table_string()`` or 

- ``interface.state_machine.to_plantuml_string()``. 

Note that function arguments in the above methods are strings! They are evaluated later to methods of the low level code interface class!
(so they must be defined there, they need not be remote function but can be ordinary class methods) 

Note that it can be convenient to use a number of sentinel methods which do not (by default) do anything, these are:

- before_get_parameter, before_set_parameter, before_set_interface_parameter, 

- before_new_set_instance, before_get_data_store_names.

(e.g. before_get_parameter is called before each parameter is retrieved.)

State models have been defined for our code's parent class, ``GravitationalDynamics``, and its grandparent class, ``CommonCode``,
and our code is simple enough that we do not need to expand this. Thus we do not need to define it.

To complete this example we will take a look at the state model of ``GravitationalDynamics``:

.. code-block:: python

    def define_state(self, handler): 
        common.CommonCode.define_state(self, handler)   
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        ...
        
        
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_method('UPDATE', 'new_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        handler.add_method('EVOLVED', 'evolve_model')
        handler.add_transition('EVOLVED','RUN', 'synchronize_model')
        ...

Note that the entry state is defined in ``CommonCode``, and we have shortened it
for readability.

Let's focus on what happens when we add a particle after evolving for some time.
We can see that ``evolve_model`` is only allowed in the ``EVOLVED`` state, which
means that that is the state after evolving. ``new_particle``, then, is allowed in
the ``EDIT`` and ``UPDATE`` states, but also when transitioning from ``RUN`` to 
``UPDATE``. From ``EVOLVED`` to ``RUN``, we can transition automatically with
``synchronize_model``. Thus, after evolving, if we add a particle, the function
``synchronize_model`` is automatically called, and we end up in the ``UPDATE`` 
state. If we call ``evolve_model`` again, we first go from ``UPDATE`` to ``RUN`` 
with an automatic call of ``recommit_particles``, and from there ``evolve_model``
leads us from ``RUN`` back to ``EVOLVED``.




Literature References
---------------------
When adding your code to AMUSE you of course want your work to be recognized.
AMUSE actively provides the references of every code used in an AMUSE script,
at the end of every run of the script. The references are defined in the Python
interface, as in the following code snippet:

.. code-block:: python

    class SimpleGravInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface):
        """
        SimpleGrav is a forward Euler code to dynamically evolve a Newtonian
        gravity N-body system. 

        .. [#] Alvarez, Etienne; Monthly Notices of the Astronomical and Astrophysical Review Letters Z, Vol. 42 (2020)
        .. [#] Adams, Douglas; Hitchhiker's Guide to the Galaxy (1979)
        """
        include_headers = ['worker_code.h']

        def __init__ (self, **keyword_arguments):
            CodeInterface.__init__(self, 
                name_of_the_worker="simplegrav_worker",
                **keyword_arguments)
            LiteratureReferencesMixIn.__init__(self)

Upon finishing a script using ``SimpleGrav`` we will get the following warning:

.. code-block::
    /home/.../amuse/src/amuse/support/literature.py:78: AmuseWarning: 

    You have used the following codes, which contain literature references:

	    "SimpleGravInterface"
		    Alvarez, Etienne; Monthly Notices of the Astronomical and Astrophysical Review Letters Z, Vol. 42 (2020)
		    Adams, Douglas; Hitchhiker's Guide to the Galaxy (1979)


	    "AMUSE"
		    Portegies Zwart, S. & McMillan, S.L.W., 2018, "Astrophysical Recipes: the art of AMUSE", 
                AAS IOP Astronomy publishing (411 pages) [2018araa.book.....P]
		    ** Portegies Zwart, S. et al., 2013, Multi-physics Simulations Using a Hierarchical Interchangeable 
                Software Interface, Computer Physics Communications 183, 456-468 [2013CoPhC.183..456P]
		    ** Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software Environment, 
                Astronomy and Astrophysics 557, 84 [2013A&A...557A..84P]
		    Portegies Zwart, S. et al., 2009, A multiphysics and multiscale software environment 
                for modeling astrophysical systems, *New Astronomy*, **Volume 14**, **Issue 4**, 369-378 [2009NewA...14..369P]

      warnings.warn(prefix + self.all_literature_references_string(), exceptions.AmuseWarning)

Note how the ``.. [#]`` denotes each entry for the literature list.
==================================
Setting values for grid boundaries
==================================

.. highlight:: python

.. note::

    At the time of writing of this tutorial (september 2012), only the
    Athena and Capreole hydrodynamics codes support the user
    defined boundaries.
    
In AMUSE, you can specify how a hydrodynamics grid code handles the 
boundary conditions. Hydrodynamic grid codes define a number of 
special grid cells at all ends of the grid. These cells are often 
called boundary or ghost cells and are handled differently by the 
code. The flow-variables in these cells are not evolved but 
calculated from other cells or from a user defined function.

Supported boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each code in amuse will need to provide at least 3 kinds of boundary 
conditions:

periodic
    To simulate periodically repeating problems. Periodic boundary 
    conditions simulates as if the grid is connected to a copy of 
    itself, often implemented by copying over the data at the 
    opposite site of the grid to the boundary cells.
reflect
    Mirrors all grid variables in the boundary cells, often 
    implemented by copying the data of the cells connected to the 
    boundary into the boundary cells.
outflow
    Sets the derivative of the flow variables with respect to the 
    direction of the boundary to zero.

The following picture shows all three boundary conditions for a one
dimensional grid with 4 boundary cells:

.. image:: boundaries1.png

Boundary conditions are defined by setting a boundary parameter for
each boundary of the grid of a code. These boundaries are defined:

xbound1
    The x-axis boundary on the left side of the grid, often called
    the inflow boundary
xbound2
    The x-axis boundary on the right side of the grid, often called
    the outflow boundary
ybound1
    The y-axis boundary on the bottom side of the grid, often called
    the inflow boundary
ybound2
    The y-axis boundary on the top side of the grid, often called
    the outflow boundary
zbound1
    The z-axis boundary on the front side of the grid, often called
    the inflow boundary
zbound2
    The z-axis boundary on the back side of the grid, often called
    the outflow boundary
    
.. image:: boundaries2.png

.. code-block:: python
    
    from amuse.lab import *
    
    code = Athena()
    code.parameters.xbound1 = "periodic"
    code.parameters.ybound1 = "reflect"

You can also set the boundary per axis as an ``(inflow, outflow)`` parameter
with the ``x_, y_, z_boundary_condition`` parameters:


.. code-block:: python

    from amuse.lab import *
    
    code = Athena()
    
    instance.parameters.x_boundary_conditions = ("periodic","periodic")
    instance.parameters.y_boundary_conditions = ("periodic","periodic")
    instance.parameters.z_boundary_conditions = ("periodic","periodic")
    
    
Custom boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can specify a custom boundary condition by setting the boundary 
condition parameters to ``interface``. After these parameters have 
been set and committed you can change the boundary values by filling 
in the correct boundary grid. These boundary grids do not overlap 
but patch together to form a larger cuboid around the grid cuboid. 
The x-boundary conditions have N boundary cells in the x direction 
and the same amount of cells in the y and z direction as the grid. 
The y-boundary conditions have N boundary cells in the y direction 
and in the x direction it  has 2 * N + the number of cells in the x 
direction of the grid, in the z direction it has the same amount of 
cells as the grid. Finally the z-boundary conditions have N boundary 
cells in the z direction and in the x and y directions it  has 2 * N 
+ the number of cells in the x or y direction of the grid.

.. image:: boundaries3.png

You can get the boundary grid of a specific boundary by calling 
`get_boundary_grid` on the code. 

.. note::

    You fill the returned grid by copying over data from a memory 
    grid, as you can only set all attributes of the grid in one go 
    (you cannot set the individual attributes yet at we did not 
    implement the required methods, september 2012)

En example of filling the boundary grid:

.. code-block:: python

    from amuse.lab import *
    
    instance=Athena()
    instance.parameters.mesh_size = (10 , 20, 10)
    instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
    instance.parameters.x_boundary_conditions = ("interface", "outflow")
    
    # request the boundary grid
    # it will have shape of 4 x 20 x 10 (x, y, z)
    # as Athena has a 4 cell boundary depth
    xbound1_grid = instance.get_boundary_grid('xbound1')
    
    # copy the grid to memory, so we can manipulate it easier
    memxbound1 = xbound1_grid.copy()
    
    # just set all cells in the grid to the same values
    memxbound1.rho = 0.02 | density
    memxbound1.rhovx = 0.2 | momentum
    memxbound1.rhovy = 0.0 | momentum
    memxbound1.rhovz = 0.0 | momentum
    memxbound1.energy =  p / (instance.parameters.gamma - 1)
    memxbound1.energy += 0.5 * (
            memxbound1.rhovx ** 2 +
            memxbound1.rhovy ** 2 + 
            memxbound1.rhovz ** 2 
        ) / memxbound1.rho
    
    # copy over the data so that the code has the correct
    # boundary values
    channel = memxbound1.new_channel_to(xbound1_grid)
    channel.copy()
    
A boundary grid connects to the grid in a very specific way. The 
general rule for the *1* or *inflow* boundaries is that the last 
cell of the boundary will connect to the first cell of the grid.  
The general rule for the *2* or *outflow* boundaries is that the 
first cell of the boundary will connect to the last cell of the grid.

For example for xbound1:

xbound1
    The left, inflow boundary on the x-axis. Will have cells ``[0, 
    # boundary cells - 1]`` in the x direction. The last cell of 
    the boundary will connect to the first cell of the grid.
ybound1. 
    The right, outflow boundary on the x-axis. Will have cells 
    ``[0, # boundary cells - 1]`` in the x direction. The first 
    cell of the boundary will connect to the last cell of the grid.
    
A boundary grid will have correct x, y and z positions that match the
grid. So, if a 1D grid starts at 0.0 and the cell size is 1.0, 
the last grid cell on the left(the inflow boundary) will have position -0.5.
 
.. code-block:: python

    from amuse.lab import *
    
    instance=Athena()
    instance.parameters.mesh_size = (10 , 1, 1)
    instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
    instance.parameters.x_boundary_conditions = ("interface", "interface")
    
    # request the boundary grid on both sides
    # it will have shape of 4 x 1 x 1 (x, y, z)
    # as Athena has a 4 cell boundary depth
    xbound1_grid = instance.get_boundary_grid('xbound1')
    xbound2_grid = instance.get_boundary_grid('xbound2')
    
    print xbound1_grid[...,0,0].x
    print instance.grid[...,0,0].x
    print xbound2_grid[...,0,0].x
    
    
Interactive tutorial
====================

This is the static description of the interactive tutorial for AMUSE.

The interactive tutorial is always delivered with the binary distribution
and you can start it in that distribution with:

.. code-block:: sh

    ./amuse-tutorial
    
For the source distribution you can also run the interactive tutorial 
but you will need to install IPython notebook first (`ipythonnb`_). 
Once you have installed ipython you can run the notebook with:

.. code-block:: sh

    export PYTHONPATH=$PWD/src:$PYTHONPATH
    
    cd doc/interactive_tutorial/
    
    ipython notebook --pylab inline 

.. note::
    
    You will need to run the ipython command from the 
    interactive tutorial directory.
    You also need to have amuse in your python path.

.. note::
    
    The ipython notebook code depends on the `tornado` and the `pyzmq` 
    modules. You will need to install these before 
    installing the ipython notebook.
    
.. toctree::
    :maxdepth: 1

    01-Loading_AMUSE
    02-Quantities_with_units
    03-Generic_units
    04-Collections_of_Particles
    05-Attributes_and_functions_on_particle_collections
    06-Using_a_community_code
    07-Channels
    08-Grids


.. _ipythonnb: http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html

=====================
Architecture Overview
=====================

Layers
------
The AMUSE architecture is based on a layered design with three layers. 
The highest layer is a python script, written for a single problem
or set of problems. The next layer contains the AMUSE code, this layer
provides a library of objects and function to use in the python script. 
The last layer contains all the existing or legacy codes. In this layer
the physical models are implemented.

Each layer builds upon a lower layer, adding functionality or ease of
use to the previous layer:

.. raw:: latex

   \begin{figure}[htp] \centering
   
.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [shape=box, style=filled, fillcolor=lightyellow, width=3];
      subgraph cluster0 {
            style=filled;
            color=azure2;
            labeljust="l";
            label="Layer 1";
            level1 [label = "User Script"];
      } 
      subgraph cluster1 {
            style=filled;
            color=azure2;
            labeljust="l";
            label="Layer 2";
            level2 [label = "AMUSE Library"];
      }
    
      subgraph cluster2 {
            style=filled;
            
            color=azure2;
            labeljust="l";
            label="Layer 3";
            level3 [label = "Community Codes"];
      }
      level1 -> level2 -> level3
    }

.. raw:: latex

   \caption{The 3 layers in AMUSE}
   \end{figure} 


Each layer has a different role in the AMUSE architecture:

1. **User Script layer**. The code in this layer implements a specific 
   physical problem or set of problems. This layer contains the 
   example scripts and scripts written by the user. This layer is 
   conceptually comparable to a User Interface layer in applications 
   with a GUI. Coupling two or more codes happens in this layer (with
   the help of support classes from the *AMUSE Library Layer*.

2. **AMUSE Library layer**. This layer provides an object oriented 
   interface on top of the legacy codes. It also provides a library
   of functionalities, such as unit handling and data conversion.
   The role of this layer is very generic, it is not specific 
   for one problem or for one physical domain.
   
3. **Community Codes layer**.  This layer defines the interfaces
   to the community codes and contains the actual codes. It provides 
   process management for the community codes and functional
   interfaces to these. The code in this layer is generic in 
   respect to problems, but specific for different physical domains.

The following sections contain a detailed explanation of the layers,
starting with the lowest layer to the highest. Some details are further
worked out in other chapters or in the reference manual.

Community codes layer
*********************

The **Community Codes layer** contains the actual applications and 
the functionality to communicate with these applications. This
layer exposes every community code as a set of functions. These 
functions are grouped in one class per code.

The AMUSE framework code and the community codes are designed to be run 
as separate applications. The AMUSE framework code consists of a python
script and the AMUSE library. The community codes consist of 
the original code-base of a scientific code extended with a 
new main application that handles messages send to it from 
the python library. 
Function calls into the community codes are send via a *message passing framework* 
to the actual running codes. 

.. raw:: latex

   \begin{figure}[htp] \centering

.. graphviz::

   digraph layers4 {
      fontsize=10.0;
        rankdir="LR";
        node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow];
        subgraph cluster0 {
            style=filled;
            color=azure2;
            label="Application";
            
            "Python Interfaces";
        }
        "Message Passing Framework";
        subgraph cluster1 {
            style=filled;
            color=azure2;
            label="Application";
            "Community Code";
        }
        "Python Interfaces" -> "Message Passing Framework";
        "Community Code" -> "Message Passing Framework";
        
        "Message Passing Framework" -> "Community Code";
        "Message Passing Framework" -> "Python Interfaces";
    }

.. raw:: latex

   \caption{The AMUSE script and community codes are separate applications. The application
   communicate using a message passing framework}
   \end{figure} 

The number of applications started and the machines on which these
run can all be set dynamically in AMUSE. Depending on the problem
a researcher can run all of AMUSE on a single desktop computer
or in a mixed environment with clusters of computers. Every AMUSE run
starts with one python script. This script can in turn start
a number of different community codes (as separate applications).
A complete run can consist of multiple applications running in parallel
or in sequence and managed by one python script. 


.. graphviz::

   digraph multiples {
      fontsize=8.0;
        rankdir="LR";
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];
        subgraph cluster0 {
            style=filled;
            color=azure2;
            label="Application";
            
            "Python Script";
        }
        subgraph cluster1 {
            style=filled;
            color=azure2;
            label="Application, running on a GPU";
            "Gravitational Dynamics";
        }
        subgraph cluster2 {
            style=filled;
            color=azure2;
            label="Application, running on a cluster";
            "Hydrodynamics";
            "Hydrodynamics 1";
            "Hydrodynamics 2";
            "Hydrodynamics 3";
            "Hydrodynamics 4";
            "Hydrodynamics" -> "Hydrodynamics 1"
            "Hydrodynamics" -> "Hydrodynamics 2"
            "Hydrodynamics" -> "Hydrodynamics 3"
            "Hydrodynamics" -> "Hydrodynamics 4"
        }
        subgraph cluster3 {
            style=filled;
            color=azure2;
            label="Application";
            "Stellar Evolution 1";
        }
        subgraph cluster4 {
            style=filled;
            color=azure2;
            label="Application";
            "Stellar Evolution 2";
        }
        subgraph cluster5 {
            style=filled;
            color=azure2;
            label="Application";
            "Stellar Evolution 3";
        }
        subgraph cluster6 {
            style=filled;
            color=azure2;
            label="Application";
            "Stellar Evolution 4";
        }
        "Python Script" -> "Gravitational Dynamics";
        "Python Script" -> "Hydrodynamics";
        "Python Script" -> "Stellar Evolution 1";
        "Python Script" -> "Stellar Evolution 2";
        "Python Script" -> "Stellar Evolution 3";
        "Python Script" -> "Stellar Evolution 4";
        
    }

Message passing
~~~~~~~~~~~~~~~
The amuse framework interacts with legacy codes via a message passing
framework. Function calls in the python scripts are translated
to messages and these messages are send to the community codes 
using the message passing framework . The community codes wait 
for message events and will decode the message upon arrival and 
perform the requested function. The results will be send back 
using a similar message.

.. image:: message_passing.png


AMUSE Library layer
*******************

The **Library layer** is responsible for providing an object
oriented interface to the community codes. It also provides extra 
functionality to help write a user script, such as file input 
and output of common file formats and unit conversions. These
extra functionalities can be used independent of the community codes.

Every community code has a *low-level* interface (defined in the community
interface layer) and an *object-oriented* interface. The *low-level*
interface is defined as as set of functions. The *object-oriented* interface
uses these functions and combines these with models for state-transitions,
units and data sets to provide an interface that is easier to use (less error
prone) and easier to couple with other codes.

.. graphviz::

    digraph amcode_0{
      fontsize=10.0;
        compound=true;
        ranksep=1;
        
        node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow];
        subgraph cluster0 {
            style=filled;
            color=azure2;
            label="Object Oriented Interface";
            "Unit Conversion";
            "Code Interface";
            "State";
            "Exceptions";
        }
        
        subgraph cluster1 {
            style=filled;
            color=azure2;
            label="Data Model";
            "Particles" ;
            "Grid Points";
        }
        "Legacy Interface";
        "Code Interface" -> "Community Code Interface"[ ltail=cluster0];
        subgraph cluster2 {
            style=filled;
            color=azure2;
            label="Support";
            
            "Input/Ouput";
            "Units";
        }
        
        "Code Interface" -> "Particles"[lhead=cluster1, ltail=cluster0]; 
        "Code Interface" -> "Units"[lhead=cluster2, ltail=cluster0]; 
    }

Model of a community code
~~~~~~~~~~~~~~~~~~~~~~~~~

The community codes of every module in all physical domains are modelled using
the same template. The template consists of attributes and wrappers. **Attributes**
provide a common interface for sub-parts of the code, for example the *particles*
attribute provides an interface to add, update and remove the particles in
a code. Attributes combine several functions in a legacy interface into
one object.  **Wrappers** are defined on top of the community functions and
add functionality to existing methods. For example for every method 
the units of the arguments and return values can be defined in a filter.
Wrappers add functionality to individual methods.


Attributes
++++++++++
The template divides the interface object of a code into
a number of attributes. Each attribute refers to an object implementing
a specific sub-interface of the code. For example a code can have a
*parameter* attribute, this attribute implements the *ParameterSet* 
sub-interface. The *ParameterSet* sub-interface defines how to
interact with the parameters of a code (in this case each parameter
can be set or queried from the set by name using normal python attribute access).

The *template* for all codes is divided into the following sub-interfaces:

parameters
    Parameters influence how the code works. Parameters are usually
    set just after creating a code. Parameters should be read-write
    or write-only.
    
properties
    Properties inform the user about the state of the code. The
    current model time is a property. Properties are always read-only.

particle sets
    Particle sets provide a common interface for a set of particles
    in the code. A code can have multiple particle sets defined under
    different names (for example gas, stars and dark matter)

grids
    Grids provide access to multi-dimensional data. A code can
    have multiple grids defined in a hierarchy (for AMR or SMR codes)


Wrappers
++++++++

Wrappers decorate a method. Wrappers can do pre- and post-processing of
the arguments or decide if a method can safely be called.

units and error code
    Defines a unit for each argument of the wrapped method. When called
    the arguments will be converted to numbers in the correct unit. The
    return values will be converted to quantities (numbers with a unit).
    
state
    The state of a code determines which functions are valid to call
    and how the code can transfer from one state into another. For example, 
    a code might give incorrect answers if the potential energy is requested before
    the particles are entered into the code, the state model will raise an error to
    inform the script writer of this problem.
    
Implementation
++++++++++++++

The implementation of the *object-oriented* interface is based on the adaptor 
pattern. A *Community Code Interface* class is adapted to create a
class which provides *"parameters"*, *"particle sets/gridpoints"* , 
*"methods with units"* , *"properties with units"* , 
*"state control"* and *"Unit conversions for incompatible unit systems"*. 
Each functionality has the same interface for all codes in the system.

.. raw:: latex

   \begin{figure}[htp] \centering
   
.. graphviz::

    digraph amcode_1{
      fontsize=10.0;
        compound=true;
        ranksep=1.0;
        rankdir="LR";
        node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow];
        subgraph cluster1 {
            style=filled;
            color=azure2;
            label="Adaptor";
            labelloc="b";
            labeljust="r";
            "Particles or Gridpoints"
            "Parameters"
            "Methods with Units"
            "Properties with Units"
            "State control"
            "Unit conversions for incompatible unit systems"
        }
        "Community Code Interface"
        "Community Code Interface"->"Particles or Gridpoints"[lhead=cluster1];
    }


.. raw:: latex

   \caption{A legacy interface is adapted to provide an object
   oriented interface and more functionality.}
   \end{figure} 


User Script
***********

The final layer is the **User Script Layer** this layer contains 
all the scripts written by a researcher for a specific 
problem or set of problems. These scripts are always written
in *python* and can use all the functionality provided by the two
lower layers in the AMUSE framework. The scripts don't need to 
follow a fixed design.

.. image:: user_script_sequence.png




=============
MPI interface
=============

The interface between the AMUSE python core and the legacy-codes is 
based on the MPI framework. Choosing MPI and not SWIG (or any other direct
coupling method) has several advantages:

* MPI is a well-known framework in the astrophysics community. 
  Other coupling methods are less well known (like SWIG)
* Legacy code does not run in the python space (memory usage, names)
* Multiple instances of the same legacy code can easily be supported (not so
  in SWIG / f2py couplings)
* Multi-process support taken into account at the start of 
  the project.
* Coupling is much looser.

There are also be some disadvantages:

* Need to define a protocol over MPI
* More "hand-work" needed to couple code. Other frameworks, like SWIG and f2py,
  generate an interface based on the application code.
* More overhead for every call, slower calls

These disadvantages are mitigated by creating a library that handles
most of the coupling details. This library has a Python, C++ and
Fortran version. It implements the protocol and generates
hooks to connect with the legacy codes.

The overhead per call may be an important factor in the speed of the
framework. This will be tested during development of the first codes. It
should be possible to limit the overhead by sending a lot of data per call. For
example, setting the properties of a lot of stars in one call. Calling a lot of methods
with limited data will be compared to sending one method with a lot of data.


==============
Coupling Codes
==============

The design for coupling codes in AMUSE is based on providing the same set of functions for
every community code and using these to devise different coupling methods. As the 
coupling methods are not fixed and can change on a per problem basis the functions
to be very generic. The AMUSE library defines three sets of functions to support coupling
codes:

particle or gridpoint manipulation
    Most properties of particles (or gridpoints) can be queried and updated during the run, 
    providing a direct method of manipulating the data of a community code. Further most
    codes support removing and adding particles during the run.

stopping conditions
    Stopping conditions are designed to interrupt a code during model evolutions. Stopping
    conditions are triggered when a code encounters a predefined state (for example
    a particle escaping out of the bounding box). 
    
services
    Services are functions added to a code that use the model of that code to provide
    data for other codes. For example a smoothed particle hydrodynamics code can provide
    the state of the model at any random point (not just on the particles) which can be
    used to create a grid from an particle model.
    
    .. _design-label:

Design documentation
====================

.. toctree::
   :maxdepth: 2
   
   introduction
   architecture
   coupling
   packages
   sets
   mpi
   

Datamodel
=========


All data is stored in sets
--------------------------
In the datamodel of AMUSE all data are stored in sets. The sets
store objects and the values of attributes belonging to the objects. 
All objects in a set have the same attributes, but not same values for
these attributes.

.. graphviz::
    
    digraph set0{
      fontsize=10.0;
      node [shape=record, fontsize=10.0,style=filled, fillcolor=lightyellow, width=3];
      struct1 [label="{set | object 1 | object 2 | .. | object n} | {attribute a | value 1.a| value 2.a | .. | value n.a} | {attribute b | value 1.b| value 2.b | .. | value n.b}| {.. | ..| .. | .. | ..} |  {attribute z | value 1.z| value 2.z | .. | value n.z}"];
    }
    

Like the relation database model
--------------------------------
For every object in a set, the set will store the values of the
attributes of the object. This model is like a relation 
database with Tables (sets in AMUSE), Records (an object in the set)
, Columns (an attribute of an object) and Fields (the value of an attribute of an object). 

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow, width=1.5];
      subgraph cluster0 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="AMUSE Model";
            "Set";
            "Object";
            "Attribute";
            "Value";
            "Set" -> "Attribute" [fontsize=10.0,label="defines"];
            "Object" -> "Value" [fontsize=10.0,label="has"];
            "Set" -> "Object" [fontsize=10.0,label="contains"];
            "Set" -> "Value"  [fontsize=10.0,label="stores"];
      } 
      subgraph cluster1 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Relation Database Model";
            "Table";
            "Record";
            "Column";
            "Field";
            "Table" -> "Column" [fontsize=10.0,label="defines"];
            "Table" -> "Record" [fontsize=10.0,label="stores"];
            "Record" -> "Field" [fontsize=10.0,label="has"];
      }
    
    }

Objects are views
------------------
Objects from a set do not store any values, instead they defer
to the set to provide their attribute values. In a sense these
objects are pointers to a location in the set. When comparing
to the relational database model an object is like
a cursor. It can be used to access the values of the attributes
belonging to the object stored in the set.

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow, width=3];
      style=filled;
      color=azure2;
      "Set";
      "Object";
      "Object" -> "Set" [fontsize=10.0,label=" pointer to a location in"];
    }

When a user asks an object for its mass the object will query the
set for the stored value and return the answer to the user.

.. image:: objects_in_set.png

Objects have a key
------------------
The objects in a set can be identified with a unique key. All objects 
having the same key are seen as the same object by the AMUSE system.
The same object can exist in multiple sets. In each set this object
can have a different value of an attribute or different attributes. 
Or, in each set a different version of the object can exist. 

Sets use Storage Models
-----------------------
The actual storage of attribute values in a set is provided by a storage
model. The set provides the interface to the script writer, the 
storage model manages the data. Each storage model must decide how and
where to store the data. All data can be stored in the memory area
of the script or in the mememory area of the code or on a file or in
a relational database.

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow, width=1.5];
      subgraph cluster0 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="AMUSE Model";
            "Set";
      } 
      subgraph cluster1 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Storage Models";
            "In Memory";
            "In Code";
            "In File";
            "In Database";
      }
      "Set" -> "In Memory";
      "Set" -> "In Code";
      "Set" -> "In File";
      "Set" -> "In Database";
    }


Selections on the set
---------------------
The datamodel provides subsets to handle a selection of the objects in a set.
When comparing to the relational database model an subset is like a view.
The subset does not store any data, all the data is stored in the original
set. When an attribute is updated in a subset, the attribute is also
updated in the original data.

 .. graphviz::
 
    digraph set0{
        fontsize=10.0;
        node [shape=record, fontsize=10.0,style=filled, fillcolor=lightyellow, width=4];
        subgraph cluster0 {
            fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Subset";
            struct1 [label="{set |<here>  object 2 | .. | object m}"];
        }
 
        subgraph cluster1 {
            fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Original Set";
            struct2 [label="{set | object 1 |<there> object 2 | .. | object n} | {attribute a | value 1.a| value 2.a | .. | value n.a} | {attribute b | value 1.b| value 2.b | .. | value n.b}| {.. | ..| .. | .. | ..} |  {attribute z | value 1.z| value 2.z | .. | value n.z}"];
        }
        
        struct1:here:e -> struct2:there:w
    
    }
    
=============
Introduction
=============


In this document we will describe the high level design of AMUSE. During
the development period of AMUSE this document will be a 
"Work in Progress". It will be updated to state latests ideas about
the design and reflect the current implementation. More detailed 
documentation can be found in the reference documentation.

AMUSE
-----
AMUSE combines existing astrophysical numerical codes
into a single system. 

Goal
~~~~

To develop a toolbox for performing numerical astrophysical
experiments.  The toolbox provides:

* A standard way of input and output of astrophysical data.
* Support for set-up and management of numerical experiments.
* A unique method to couple a wide variety of physical codes
* A legacy set of standard, proven codes. These codes will be integrated 
  into AMUSE as modules. Each module can be used stand-alone or in
  combination with other modules
* A standard way for adding new modules to AMUSE.
* Examples to show the use of each module and possible couplings 
  between the modules.
* Documentation containing introduction, howtos and reference documents.


Development
~~~~~~~~~~~
AMUSE is originally developed at the Leiden Observatory. The Leiden Observatory is
a faculty of the Leiden University in the Netherlands. Funding is provided
by a NOVA grant.

.. image:: ../logos/universiteit_leiden_logo.png
   
.. image:: ../logos/strw_logo.png
   :width: 2.5cm

   
.. image:: ../logos/nova_logo.jpg
   :width: 2.5cm

===============
Python Packages
===============

Like all large python projects, the AMUSE source-code is structured 
in packages, subpackages and modules. Packages and subpackages are
directories in the filesystem and modules are python files (A detailed
explanation can be found in 
`Modules <http://docs.python.org/tutorial/modules.html>`). In this
chapter the most important packages are named. 

The source code of the **AMUSE Code** and **Community Codes** layer is 
combined in one package. The package name is **amuse**.

amuse
    Root package of the AMUSE code, does not contain test files,
    the build system or the test system.

The **amuse** and the package is further divided into three
subpackages:

amuse.support
    Contains the code of the **AMUSE Code** layer. The units system,
    data io and model and all base classes can be found in 
    this package.
    
amuse.community
    Contains the code the **Community Codes** layer. The community codes
    can be found in this package as well as support code for generating
    the script to legacy code messaging framework.
    
amuse.ext
    Contains **extra** and/or **extension** code. For example, making
    initial data models is not one of the main functionalities of AMUSE,
    but it is very useful to include this into the codebase.

===============================================
Support code for the community code interfaces
===============================================

.. automodule:: amuse.rfi.core
    :members:
===================================================
Message protocol between python and community codes
===================================================

Introduction
------------

The implementation of the interfaces of the community codes is based
on sending messages. These messages are encoded as MPI primitive datatypes
and send to each code using MPI. In this section we will describe the
overall operation of the interface implementation and specify the message
format.

Overall Operation
-----------------
The method call interface is a request/response protocol. For every method
call a message is send to the code. This message is decoded 
by the code into a function-id and arguments. The code will call 
the function with the function-id using the decoded arguments. After the
function completes a result message is returned containing all the
results and the function-id of the called function. This process is depicted
in the following table.

============================= === =============
Python script                 ..  Code
(client)                      ..  (server)
============================= === =============
start of function call
encode arguments
send MPI messages              >  
..                                receive MPI messages
..                                decode messages
..                                handle (setting data, evolving the solution)
..                                encode response messages
..                             <  send response messages
receive MPI response messages 
decode response message
return result to script
end of function call
============================= === =============


Message format
--------------
Every message sent between the python script and a code has the
same format. This format consists of a header and zero or more
content arrays. The header contains the function-id, 
the number of calls to the same function and 
the number of values (arguments or results) per primitive type for
one function call. Every content array contains the sent values of 
a primitive datatype. For example, a content array with 
all the integer values in the arguments of the function. 

.. code-block:raw
    HEADER
        function id
        number of calls
        number of values of type 1
        ...
        number of values of type n

    CONTENT-ARRAY 1
        int 1
        ...
        int m
    ...
    CONTENT-ARRAY n
        double 1
        ...
        double m

Message header
~~~~~~~~~~~~~~
The header is sent with a MPI broadcast message. The header consists
of an array of ``n + 2`` integers. ``1`` integer to specify the function, ``1`` 
integer to specify the number of calls and ``n`` integers to 
specify the number of arguments of each type. Version 0.2 of the
interface contains support for 4 types (float64, int32, float32 and string)
The header is a 5 integer long array, the specification of each integer is given
in the following table:

*Message header*

=========  ===============
Position   Description
=========  ===============
0          function-id
1          number of calls
2          number of arguments/results of type float64
3          number of arguments/results of type int32
4          number of arguments/results of type float32
5          number of arguments/results of type string
=========  ===============


Content array
~~~~~~~~~~~~~
The arguments are sent with a MPI broadcast message, the results are
sent with a MPI send message. The arguments or results are sent when
1 or more values are needed for a function. If no values are needed for 
a type, no message is sent for that type.

Encoding of the content arrays
******************************
To sent the arguments or results values of a function, the values must
be encoded in arrays (each of a different primitive type).

The arguments of a function are normally not sorted by type. The first argument
may be an integer, the second a double and the last an integer. The message 
format does specify a fixed sorting, all float64 values are sent first, then 
the integers, then the other types. To sent the arguments or result, the
values are encoded following a fixed scheme.

The arguments are encoded by extracting all values of a primitive type in the 
order they occur in the function definition. This is done for every type. For 
example when the first argument of a function is an integer, the second
a double and the last an integer, two content arrays will be sent.
One for the two integers and one for the single double. 
The integer array has at the first position the first argument 
and at the second position the last argument. The double array has at the
first position the second double.


For this function::

  int example_function(int id, double x, int type)

Two content arrays are sent::

  int[id, type] 
         (the first argument and the last argument 
          to the method are integers)
          
  double[x] 
         (the second, argument to the method is a double)

The arguments are encoded in order, going from left to right in C or
fortran function definition. A content array is as long as the number of arguments 
or results of a primitive type, the specification of member in a content array
is given in the following table:

*Content Array*

=========  ===============
Position   Description
=========  ===============
0          first argument of type X
1          second argument of type X
...
n          last argument of type X
=========  ===============



Multiple calls to the same function
***********************************
The MPI messaging has a significant overhead. To reduce this overhead
the arguments and results of a number of calls to the same method
can be encoded in one set of messages (header and content-arrays).

Creating arrays of values:
    
    x = [i for i in range(1000)]
    y = [i * 2 for i in range(1000)]
    z = [i * 3 for i in range(1000)]

Calling the same function multiple times with different values for the arguments::

    for i in range(100):
        instance.add_position(x[i], y[i], z[i])

Can be converted to calling the function once with an array of arguments::
    
    instance.add_position(x, y, z)

The encoding of the arguments for the call with arrays
follows the same strategy as the call with single value. The values of
the first argument are encoded first, the values of the second found argument 
of a type are encoded second etc. 

*Content Array format, when message is encoded for multiple calls to the same function*

=========  ===============
Position   Description
=========  ===============
0          first value in the array of the first argument of type X
1          second value in the array of the first argument of type X
...
m - 1      last value in the array of the first argument of type X
m          first value in the array of the second argument of type X
m + 1      second value in the array of the second argument of type X
...
n * m      last value in the array of the last argument of type X
=========  ===============

To get the value of the Mth value of the Nth argument of a type (starting 
to count at zero, n = 0 is the first argument, m = 0 is the first value of 
the argument)::

    value = array[ n * len + m ]
    
This encoding degrades into the case for the single call when len = 1 (m = 0, as
the array contains only one value)::

    value = array[ n ]


Examples
--------
TBD


Interface Specifications
========================

.. toctree::
    :maxdepth: 2

    introduction_interface_specification
    stellar_dynamics_interface_specification
    stellar_evolution_interface_specification
    hydrodynamics_interface_specification
    radiative_transfer_interface_specification
   

==========================
Datamodel
==========================


Introduction
------------

Containers
~~~~~~~~~~
The AMUSE datamodel is based on objects stored in containers. The 
containers store all the relevant data of the objects. One can work 
with the data in a container as a whole or with an individual object 
in a container. When working with an individual object all data will 
be retrieved from and stored in the container of that object. *This 
is different from normal python objects and lists where the lists 
store references to the objects and the data is stored on the 
objects*.

Amuse containers::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(2)
    >>> stars.mass = [2, 3] | units.MSun
    >>> sun = stars[0]
    >>> sun.mass = 1 |  units.MSun
    >>> print stars.mass
    [1.0, 3.0] MSun
    >>> stars.mass = [0.5, 5] | units.MSun
    >>> print sun.mass
    0.5 MSun
    
Python lists::

    >>> from amuse.datamodel import Particle
    >>> sun = Particle(mass = 1 | units.MSun)
    >>> star = Particle(mass = 3 | units.MSun)
    >>> stars = [sun, star]
    >>> print stars.mass # cannot access attributes on the list
    Traceback (most recent call last):
    AttributeError: 'list' object has no attribute 'mass'
    

Set or grid
~~~~~~~~~~~
AMUSE supports two kinds of containers. One container kind 
implements a one dimensional set of objects. You can add and remove 
objects from the set and combine sets in different ways. This 
container kind is called **'Particles'** and the individual objects 
stored are called **'Particle'**. The other container kind 
implements a multidimensional grid of objects. You cannot change the 
shape of a grid and objects cannot be added or removed from a 
grid. This container kind is called **'Grid'** and individual 
objects stored in the container are called **'GridPoint'**.  

With particles::
    
    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(2)
    >>> stars.mass = [2, 3] | units.MSun
    >>> print stars.mass
    [2.0, 3.0] MSun
    
With grids::

    >>> from amuse.datamodel import Grid
    >>> from amuse.units import units
    >>> field = Grid(2,2)
    >>> field.density = [[1, 2], [3, 4]] | units.g / units.cm ** 3
    >>> point = field[1,0]
    >>> point.density = 5  | units.g / units.cm ** 3
    >>> print field.density
    [[ 1.  2.], [ 5.  4.]] g / (cm**3)
    

Memory, code or file
~~~~~~~~~~~~~~~~~~~~
The containers in AMUSE can save the data for individual objects in 
different kinds of **'Stores'**. AMUSE currently implements 3 kinds 
of stores. One kind of store saves all data in memory, this kind is 
used by default and provides the fastest way to work with the 
containers. Another kind of store saves and retrieves all data from 
a community code. It uses MPI messages to communicate this data. 
This kind is used by default for containers provided by a code and 
is the primary way by which you can interact with the data inside a 
running code. The last kind of container stores and retrieves data from
an HDF5 file. This kind is only used when saving or loading a container.

Varying attributes
~~~~~~~~~~~~~~~~~~
In memory containers can support a user defined set of attributes 
for the contained objects. You can define a new attribute for the 
container by assigning a value to the attribute. In code containers 
support a pre-defined, per code set of attributes for the contained 
objects. You cannot define a new attribute on these containers. Also 
as a code may not allow some attributes to be set individually, the 
container also cannot set some attributes individually. For example 
you often cannot set the X position of a particle, you must set the 
X, Y, and Z position in one go.

Adding attributes to a set or object::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(keys = [1,2])
    >>> # you can add an attribute by assigning to a name on a set
    >>> stars.mass = [1, 3] | units.MSun 
    >>> sun = stars[0]
    >>> print sun
    Particle(1, mass=1.0 MSun)
    >>> sun.radius = 1 | units.RSun
    >>> print sun
    Particle(1, mass=1.0 MSun, radius=1.0 RSun)

Objects not classes
~~~~~~~~~~~~~~~~~~~
Different kinds of particles or gridpoints are not modelled by 
implementing subclasses. Instead, different kinds of particles are 
defined ad-hoc, by variable naming (for example planets, stars or 
cloud_particles) and varying attributes (stars have position and 
mass, cloud_particles have position and density). This allows you to 
model your problem with the names and attributes that best fit your 
model. Unfortunately, this way of working does remove a level of 
self description in the system. To mitigate this problem the 
containers defined in community codes and in example scripts all 
follow the same conventions. We describe these conventions in a 
later section in this chapter.

Different attributes and names for different kinds::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> # stars and planets share some attributes (radius)
    >>> # but also have attributes that make sense only for
    >>> # the specific kind (population for planets, 
    >>> # luminosity for stars)
    >>> stars = Particles(keys = [1,2])
    >>> stars.luminosity = [1, 3] | units.LSun 
    >>> stars.radius = 1 | units.RSun
    >>> planets = Particles(keys = [3,4])
    >>> planets.radius = 6371  | units.km
    >>> planets.population = [7000000000, 0]
    
Identity
~~~~~~~~
Each object in a container has a unique identity, no two objects in 
a container can have the same identity. In a **'Particles'** 
container this identity is provided by a unique 64-bit key. In a 
**'Grid'** container this identity is provided by the n-dimensional 
index in the grid. Objects with the same identity can exists in 
multiple containers. These objects are considered as the same 
*conceptual* object in AMUSE. Different containers will provide 
different information about that object. For example the same *star* 
could live in a container in memory, in a container of a stellar 
evolution code and in a container of a stellar dynamics code. AMUSE 
provides several functions to link these objects and to transfer 
data between them.


Different sets store information on the same object::

    >>> from amuse.datamodel import Particles
    >>> from amuse.units import units
    >>> stars = Particles(keys = [1,2])
    >>> stars.luminosity = [1, 3] | units.LSun 
    >>> bodies = Particles(keys = [1,2])
    >>> bodies.mass = [1, 3] | units.MSun 
    >>> print bodies[0] == stars[0] # the same 'conceptual' object in different sets
    True
    >>> print bodies[0] == stars[1] # not the same object
    False
    >>> print bodies[0]
    Particle(1, mass=1.0 MSun)
    >>> # finding the coresponding particle in the stars
    >>> print bodies[0].as_particle_in_set(stars)
    Particle(1, luminosity=1.0 LSun)
    
    
Ownership
~~~~~~~~~
Objects in a container are owned by that container, the container 
controls the data and the life-cycle of each object. 




Particle keys
-------------
All particles have a unique 64-bit key. This key is created using 
a random number generator. The chances of duplicate keys using 
64-bit integers are finite but very low. The chance of a duplicate 
key can be determined by a generalization of the birthday problem.

Duplicate keys::

    >>> # given n random integers drawn from a discrete uniform distribution 
    >>> # with range [1,highest_integer], what is the probability
    >>> # p(n;highest_integer) that at least two numbers are the same?
    >>> import math
    >>> number_of_bits = 64
    >>> highest_integer = 2**number_of_bits
    >>> number_of_particles = 1000000.0 # one million
    >>> probability = 1.0 - math.exp( (-number_of_particles * (number_of_particles - 1.0))/ (2.0* highest_integer) )
    >>> print probability
    2.71050268896e-08
    >>> # can also set the probablity and determine the set size
    >>> probability = 0.00001 # 0.001 percent
    >>> number_of_particles = math.sqrt(2 * highest_integer * math.log(1 / (1.0 - probability)))
    >>> print number_of_particles
    19207725.6894
    
If you use large sets or want to load a lot of simulations with 
different particles into a script the probability of encountering a 
duplicate may be too high. You can check for duplicates in a set of 
particles by calling ``has_duplicates`` on a set. You can also change
the key generator to better match your problem.
    

Sets of particles
------------------
The AMUSE datamodel assumes all particles come in sets. The
data of a particle is stored in the set.

.. automodule:: amuse.datamodel

    .. autoclass:: AbstractParticleSet
        :members:
        
        .. automethod:: __add__
        .. automethod:: __sub__
        
    .. autoclass:: Particles
        :members:
        
    .. autoclass:: ParticlesSubset
        :members:
        
    .. autoclass:: ParticlesWithUnitsConverted
        :members:

    object
    ---------------

    .. autoclass:: Particle
        :members:
        
        .. automethod:: __add__
        .. automethod:: __sub__
    
    Methods to retreive physical properties of the particles set
    ------------------------------------------------------------
    
.. automodule:: amuse.datamodel.particle_attributes

    .. autofunction:: center_of_mass 
        :noindex:
    .. autofunction:: center_of_mass_velocity
        :noindex:
    .. autofunction:: kinetic_energy
        :noindex:
    .. autofunction:: potential_energy 
        :noindex:
    .. autofunction:: particle_specific_kinetic_energy
        :noindex:
    .. autofunction:: particle_potential
        :noindex:

    
Implementation
--------------


.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];
        
        "AttributeStorage" -> "AbstractSet";
        "AbstractSet" -> "AbstractParticlesSet";
        "AbstractSet" -> "AbstractGrid";
        
        "AbstractParticlesSet" -> "Particles";
        "AbstractParticlesSet" -> "ParticlesSubSet";
        "AbstractParticlesSet" -> "ParticlesSuperSet";
        "AbstractGrid" -> "Grid";
        "AbstractGrid" -> "SubGrid";
        
        "AttributeStorage" -> "InMemoryStorage" [label = "implements"];
        "AttributeStorage" -> "InCodeStorage";
    }



.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];
        edge [color="dodgerblue2", fontcolor="dodgerblue2"];

        "Particles" -> "AttributeStorage" [headlabel=" 1", taillabel="1",label = "store"];
        "Grid" -> "AttributeStorage" [headlabel=" 1", taillabel="1", label = "store"];
        "SubGrid" -> "Grid" [headlabel=" 1", taillabel="1", label="view on", color="dodgerblue2"];
        "ParticlesSubSet" -> "Particles" [headlabel=" 1", taillabel="1", label="view on", color="dodgerblue2"];
        "ParticlesSuperSet" -> "Particles" [headlabel=" *", taillabel="1", label="view on", color="dodgerblue2"];
    }


.. graphviz::

   digraph multiples {
      fontsize=8.0;
        node [fontsize=8.0,shape=box, style=filled, fillcolor=lightyellow];

        "AbstractSet" -> "AttributeStorage" [headlabel=" 1", taillabel="1", label = "stored_attributes", color="dodgerblue2" , fontcolor="dodgerblue2"];
        "AbstractSet" -> "DerivedAttribute" [headlabel=" *", taillabel="1", label = "derived_attributes", color="dodgerblue2", fontcolor="dodgerblue2"]; 
        "DerivedAttribute" -> "CalculatedAttribue";
        "DerivedAttribute" -> "VectorAttribute";
        "DerivedAttribute" -> "FunctionAttribute";
        "DerivedAttribute" -> "DomainAttribute";
    }

.. image:: derived_attributes.png


.. autoclass:: amuse.datamodel.AttributeStorage
    :members:



.. _cuda-setup-label:

==========================
Setting up GPGPU with CUDA
==========================

Introduction
~~~~~~~~~~~~

Here we provide help for setting up general-purpose computing on graphics processing units (GPGPU)
using CUDA. Performing (part of) the calculations on a graphics card can result 
in a significant speed-up. Several codes in AMUSE support or require GPGPU: 
phi-GRAPE (using Sapporo), ph4 (using Sapporo), Octgrav, Bonsai, HiGPUs.


Self-help script
~~~~~~~~~~~~~~~~

In the AMUSE root directory a self-help script can be found. If building or testing any of the 
codes mentioned above fails and you wonder why, it will hopefully provide you with helpful suggestions.
From a command line run the bash script `cuda_self_help`:

.. code-block:: sh

    > ./amuse-x.0/cuda_self_help


Step-by-step
~~~~~~~~~~~~

* :ref:`CUDA-capable`
* :ref:`CUDA_SDK`
* :ref:`env-vars`
* :ref:`configure-with-cuda`
* :ref:`test`


.. _CUDA-capable:

Check that your computer has a CUDA-capable Nvidia graphics card
-----------------------------------------------------------------

First determine the model of your GPU.

On Linux:

.. code-block:: sh

    > nvidia-settings -q gpus


On Mac:

   1. Click on “Apple Menu”
   2. Click on “About this Mac”
   3. Click on “More Info”
   4. Select “Graphics/Displays” under Contents list

Check whether your GPU model is listed among 
`Nvidia's CUDA-enabled GPUs <https://www.nvidia.com/object/cuda_gpus>`_.


.. _CUDA_SDK:

Check that you have installed the CUDA Toolkit (TK) and software development kit (SDK)
--------------------------------------------------------------------------------------

If not, download and install it from `CUDA Development Tools <https://developer.nvidia.com/cuda-downloads>`_.


.. _env-vars:

Set the CUDA_TK and CUDA_SDK environment variables
--------------------------------------------------

After installing the CUDA TK and SDK, make sure the environment variables CUDA_TK and CUDA_SDK are set correctly.
For shell (bash) you need to do:

.. code-block:: sh

   export CUDA_TK=/path/to/cuda_tk
   export CUDA_SDK=/path/to/cuda_sdk

'/path/to/cuda_tk' should hold directories named 'include', 'lib', and 'lib64' (where libcudart.so is located)

'/path/to/cuda_sdk' should hold a directory named 'common/inc' (where various header files are located)

We recommend you add these lines to your '.bashrc' file so that
the variables are set correctly for all sessions. If you have a
C shell you need to do a *setenv* and edit the '.cshrc file.


.. _configure-with-cuda:

Configure AMUSE with CUDA enabled
---------------------------------

AMUSE needs to be configured with the option ``--enable-cuda``. See :ref:`configuration-gpu-label`.


.. _test:

Testing
-------

Now try building for example Octgrav and run the nosetests (from AMUSE root directory),
but first re-initialize mpd (or it will remember its original environment):


.. code-block:: sh

   mpdallexit
   mpd &
   make octgrav.code
   nosetests ./test/codes_tests/test_octgrav.py

If this fails, please contact us through the `'amusecode' google group <http://groups.google.com/group/amusecode>`_, 
or on IRC at the #amuse channel on irc.freenode.net. 
==========================
Supported File Formats
==========================

.. automodule:: amuse.io

    Introduction
    ------------
    The AMUSE framework provides a generic way of reading and writing
    sets of entities (these entities can be particles, gasclouds or
    gridpoints). AMUSE provides a function to write a set to a file
    and a function to read a set from a file. These functions need
    the name of the file and the format of the data in the file. We
    will describe these functions first in this chapter. The functions
    can throw 3 kinds of exceptions, these are described next.
    
    Use
    ---
    
    To write a data set to a space separated text file do:: 
        
        >>> from amuse.io import write_set_to_file
        >>> from amuse.datamodel.core import Particles
        >>> from amuse.units import units
        >>> x = Particles(100)
        >>> x.mass = 10.0 | units.MSun
        >>> x.radius = range(1.0,101.0) | units.RSun
        >>> write_set_to_file(x, "test.csv","txt", attribute_types = [units.MSun, units.RSun])
        
    .. autofunction:: write_set_to_file
    .. autofunction:: read_set_from_file
    .. autofunction:: get_options_for_format
    
    Exceptions
    ----------
    
    .. autoclass:: UnsupportedFormatException
    .. autoclass:: CannotSaveException
    .. autoclass:: CannotLoadException
    
    
    
    Amuse HDF5 files
    ----------------
    
    AMUSE provides it's own data format based on hdf5 files. The HDF5 file storage
    provides support for storing particle set and grids in the same file. It also
    retains the particle references and can store multiple versions of the
    same set.
    
        >>> from amuse.support import io
        >>> from amuse.ic.plummer import new_plummer_model
        >>> particles = new_plummer_model(1000)
        >>> io.write_set_to_file(
                    particles,
                   'plummer.hdf5',
                   'hdf5',
            )
        >>> io.read_set_from_file('plummer.hdf5', 'hdf5')
    
    
    .. iooptions:: hdf5
    
    
    
    Text files
    ----------
    
    AMUSE has support for reading and writing text (``.txt``, ``.csv``) files.
    You can specify the txt format by entering "txt" (for space separated files) 
    or "csv" (for comma separated files)
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file(
                   'plummer.txt',
                   'txt',
                   attribute_types = (units.m, units.m, units.kms),
                   attribute_names= ('x', 'y', 'vx')
            )
        >>> io.write_set_to_file(particles, 'plummer.txt', 'txt')
        
    .. iooptions:: txt
    
    By default, text files are stored with some unit information on 
    a comment line in the header. Unfortunately, amuse cannot depend 
    on this line to read the file. When reading a text file you 
    always need to specify the units of every column in the file. For
    example, to read a file with particle positions where all positions
    are stored in parsec, do::
    
        particles = io.read_set_from_file(
               'positions.txt',
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names= ('x', 'y', 'z')
        )
    
    When writing a text file, you do not need to specify the units. If
    you don't specify the units the ``write_set_to_file`` function will
    default to the units of the quantities being saved. For reliability 
    and reproducibility we suggest to always specify the units and
    the names of the attributes to save when saving a text file::
        
        particles = io.write_set_from_file(
               particles,
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names= ('x', 'y', 'z')
        )
        
    The amuse text file routines are very generic and since version 6.0,
    you can also specify a list of vector quantities instead of a particle
    set::
    
        particles = io.write_set_from_file(
               None,
               'txt',
               attribute_types = (units.parsec, units.parsec, units.parsec),
               attribute_names = ('x', 'y', 'z'),
               quantities = [
                    [0.1,0.2] | units.parsec, 
                    [0.3,0.4] | units.parsec, 
                    [0.5,0.6] | units.parsec, 
               ]
        )
    
    Starlab
    -------
    
    AMUSE has support for reading and writing starlab (``.dyn``) files.
    You can specify the starlab format by entering "dyn" or "starlab"
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.dyn','starlab')
        >>> io.write_set_to_file(particles, 'output.dyn', 'dyn')
    
    The starlab format support sevaral options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.dyn','starlab', must_scale = False, return_children = False)
    
    .. iooptions:: starlab

    The units of the values in the star (stellar properties)
    section of a starlab file are always in derived S.I units (solar
    mass, million years, solar luminocity etc.).
    
    The masses given in de the dynamics section of a starlab file
    are usually in *nbody* units. Some starlab tools will set the
    mass values in Solar mass units (for example the ``makemass``
    tool will return the masses in solar mass units). To read these
    files you need to set the ``dynamics_mass_units``.
    
    .. code-block:: bash
    
        > makeplummer -n 1000 > plummer1.dyn
        > cat plummer1.dyn | makemass -f 1 -x -2.0 -l 0.1 -u 20 > plummer2.dyn
        > cat plummer2.dyn | add_star -Q 0.5 -R 5 > plummer3.dyn
        > cat plummer3.dyn | scale -s > plummer4.dyn
        > cat plummer4.dyn | kira -S > plummer5.dyn
        
    The ``plummer1.dyn``, ``plummer4.dyn`` and ``plummer5.dyn`` files 
    will provide masses (and all other dynamical properties) in scaled
    *nbody* units. The ``plummer2.dyn`` and ``plummer3.dyn`` files
    will have masses in solar masses. To read each file in AMUSE, and
    return the particles with S.I. units, you need to do::
    
    
        >>> from amuse.support import io
        >>> from amuse.units import nbody_system, units
        >>> converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        >>> particles1 = io.read_set_from_file('plummer1.dyn','starlab', nbody_to_si_converter = converter)
        >>> particles2 = io.read_set_from_file('plummer2.dyn','starlab', dynamics_mass_units = units.MSun, nbody_to_si_converter = converter)
        >>> particles3 = io.read_set_from_file('plummer3.dyn','starlab', dynamics_mass_units = units.MSun, nbody_to_si_converter = converter)
        >>> particles4 = io.read_set_from_file('plummer4.dyn','starlab')
        >>> particles5 = io.read_set_from_file('plummer5.dyn','starlab')
    
    .. note::
        
        No nbody converter object is needed for the last files, as the scale 
        factors given in the files will be used.
    
    The ``plummer1.dyn``, ``plummer4.dyn`` and ``plummer5.dyn`` can also
    be read in nbody units. In the following example the returned 
    particles have dynamic attributes (mass, radius, velocity, acceleration)
    in *nbody* units::
    
        >>> from amuse.support import io
        >>> particles1 = io.read_set_from_file('plummer1.dyn','starlab')
        >>> particles4 = io.read_set_from_file('plummer4.dyn','starlab', must_scale = False)
        >>> particles5 = io.read_set_from_file('plummer5.dyn','starlab', must_scale = False)
    
    
    NEMO
    ----
    
    AMUSE has support for reading and writing nemo (``.tsf``) files.
    You can specify the starlab format by entering "nemo" or "tsf"
    as the format::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.tsf','nemo')
        >>> io.write_set_to_file(particles, 'output.tsf', 'tsf')
        
    The nemo format support several options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example::
    
        >>> from amuse.support import io
        >>> from amuse.units import nbody_system, units
        >>> converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
        >>> particles = io.read_set_from_file('plummer.nemo','tsf', nbody_to_si_converter = converter)
        
    .. iooptions:: nemo
        
        
    Gadget
    ------
    
    AMUSE has support for reading and writing gadget 2 files.
    You can specify the gadget format by entering "gadget"
    as the format::
    
        >>> from amuse.support import io
        >>> gas, halo, disk, bulge, stars, bndry = io.read_set_from_file('plummer.dat','gadget')
        >>> io.write_set_to_file(particles, 'output.dat', 'gadget')
    
    The gadget file format reader will return a tuple of gas, halo, disk, 
    bulge, stars and boundary particles. This is different form other
    readers which only return one set. The write_set_to_file can take
    a particles set (will be saved as halo particles) or a tuple of 
    gas, halo, disk, bulge, stars and boundary particles.
    
    The gadget format support several options, listed below. You can
    use these options by adding additional keyword arguments to
    the :func:`read_set_from_file` or :func:`write_set_to_file` functions. 
    For example (will read a gadget file with timestep information)::
    
        >>> from amuse.support import io
        >>> particles = io.read_set_from_file('plummer.nemo','gadget', has_timestep = True)
        
    .. iooptions:: gadget
        
======================
Source Code Management
======================




Committing Code
--------------


Git repository location
-----------------------

The AMUSE source code repository can be found at::

    https://github.com/amusecode/amuse/

========================================
Radiative Transfer Interface Definition
========================================

=========== ============ ========= =========
Date        Author(s)    Version   State
=========== ============ ========= =========
02-11-2009  AvE          0.0       TBD
=========== ============ ========= =========

Introduction
~~~~~~~~~~~~
TBD
================================
Support code for AMUSE framework
================================

.. automodule:: amuse.support.core


    .. autoclass:: late
        :members:
        
        
    .. autoclass:: print_out
        :members: 
    
        .. automethod:: __add__
        
    .. autoclass:: OrderedDictionary
        :members:
        
        
====================
Quantities and Units
====================

Introduction
------------

We want to be able to use physical quantities rather than just
measures (represented by e.g. floats or integers on computers) in
order to raise the ambiguity caused by the implicit choice of
units. Serving this purpose, AMUSE comes with a **quantity** class. Like
particle sets (fundamental data structure in AMUSE), quantities are
fundamental variables. When interacting with code all data has units,
even scaled systems. In handling quantities we regard units as being
integral part of the mathematical description of the variable
[D.C. Ipsen, Units, Dimension, and Dimensionless Numbers, 1960,
McGraw-Hill Book Company], i.e. we can state things like:

1 AU = 149597870.691 km

.. code-block:: python

    >>> from amuse.units import units
    >>> q1 = 1.0|units.MSun
    >>> q2 = 1.98892e30|units.kg
    >>> q1 == q2
    True

Quantity objects have basic conversion ability, when different units
belonging to the same dimension exist the quantity object will convert
from one to the other. For more elaborate conversion facilities, like
interfacing with a natural unit code, AMUSE provides the
generic_unit_converter and the derived nbody_system modules.  


Quantities
----------

.. inheritance-diagram:: amuse.units.quantities

.. automodule:: amuse.units.quantities
    
    .. autoclass:: Quantity
        :members:
   
    .. autoclass:: ScalarQuantity
        :members:

    .. autoclass:: VectorQuantity
        :members:

        .. automethod:: __getitem__
       
        .. automethod:: __setitem__
    
    .. autoclass:: NonNumericQuantity
        :members:
        

Units
-----
.. inheritance-diagram:: amuse.units.core

.. automodule:: amuse.units.core
    :members:

Unit systems and converters
---------------------------

The amuse framework gives you the ability to choose a unit system for
your model through the 'generic_unit_converter' module.  This enables
you to work with e.g. natural units or n-body units within AMUSE.

The implementation makes use of a dimension-space, which is a
vector-space where the chosen units form a base. For a detailed
description of the method see e.g.: Maksymowicz, A, *American Journal
of Physics*, Vol.44, No.3, **1976**.         
    
Generic units
~~~~~~~~~~~~~

.. automodule:: amuse.units.generic_unit_system

N-body units
~~~~~~~~~~~~

.. automodule:: amuse.units.nbody_system

Converter
~~~~~~~~~

.. automodule:: amuse.units.generic_unit_converter

     .. autoclass:: ConvertBetweenGenericAndSiUnits
         :members:
                     
Use with codes
~~~~~~~~~~~~~~

For convenience, the gravitational dynamics interface works with a
unit-converter which converts between the units used by the code and
the preferred units in the script (user's choice).

We show two examples. The first one uses the (derived from the generic
units converter) nbody_system converter, which is the logical choice
for dynamics (n-body) codes.

The second one uses the generic unit converter directly, this is just
an example.

~~~~~~~~~~~~~~~~
Usage: example 1
~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from amuse.community.hermite0.interface import Hermite
    >>> from amuse.units import nbody_system
    >>> from amuse.units import constants, units
    >>> convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    >>> hermite = Hermite(convert_nbody)

~~~~~~~~~~~~~~~~
Usage: example 2
~~~~~~~~~~~~~~~~

.. code-block:: python
    
    >>> from amuse.community.hermite0.interface import Hermite
    >>> from amuse.units import generic_unit_converter as gc
    >>> from amuse.units import constants, units
    >>> converter = gc.ConvertBetweenGenericAndSiUnits(constants.G, 1.0 | units.MSun, 1 | units.AU) 
    >>> hermite = Hermite(converter)

~~~~~~~
Example
~~~~~~~

More examples can be found in the tutorial, :ref:`working_with_units`

.. automodule:: amuse.community.interface.gd

    .. autoclass:: GravitationalDynamics
        :members: __init__

=============================================
Running AMUSE on The Cartesius supercomputer
=============================================

The Cartesius is the Dutch national supercomputer. For more information on e.g. obtaining an account, see the SURFSara website: https://userinfo.surfsara.nl/systems/cartesius.

Using AMUSE on the Cartesius is relatively straightforward. Below is a tested method for installing and using AMUSE using the prerequisites, though other options (e.g. using only software pre-installed on the machine) should also be possible. For a generic description of the installation of prerequisites, see the :ref:`prerequisite-label` section.

Obtaining AMUSE
---------------

We assume a copy of AMUSE has been downloaded to an `amuse` folder in the users home. For instance using git:

.. code-block:: sh

	> cd /home/USERNAME
	
	> git clone https://github.com/amusecode/amuse.git

Environment settings
--------------------

Since we will be using the AMUSE prerequisites software, we need to set some environment variables. Make sure the following lines are present in your .bashrc file:

.. code-block:: sh

	#file: ~/.bashrc

	#load java 1.7 needed for some codes
	module load java/oracle

	export PREFIX=~/amuse/prerequisites
	export PATH=${PREFIX}/bin:${PATH}
	export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

	#use gfortran  	
	export FC=gfortran
	export F77=gfortran
  	
Also make sure that .bashrc is loaded in your .bash_profile. This way, the enviroment is always set correctly, both in interactive and non-interactive mode.
 
.. code-block:: sh

	#file: ~/.bash_profile

	if [ -f ~/.bashrc ]; then
	    . ~/.bashrc
	fi

**Note: be sure to re-connect to the machine for these changes to take effect.**

Install AMUSE prerequisites
---------------------------

Next, we will install all prerequisites of amuse using the AMUSE supplied scripts.

.. code-block:: sh

	# create a directory for the prerequisites
	> mkdir ~/amuse/prerequisites

	# go to the <doc/install> directory
	> cd ~/amuse/doc/install
	
	# Start the installation script for Python.
	> ./install-python.sh

	# Download the prerequisite packages.
	> ./install.py download
	
	# Install prerequisites. Use hydra as the default MPI process manager. May take a while...
  	> ./install.py --hydra install
  	
  	# Optionally also install matplotlib
  	> ./install.py --matplotlib install
 
Configure and build AMUSE
-------------------------

Configuring and building amuse is now as normal.

.. code-block:: sh

	# go to the amuse directory
	> cd ~/amuse

	# configure amuse
	> ./configure MPIEXEC=mpiexec.hydra
	
	# build amuse
	> make

        # optionally also install codes requiring downloading files
	> make DOWNLOAD_CODES=1
	
	
Test the installation
---------------------

To test your AMUSE installation, run nosetests.

**Note: do not run simulations on the frontend of the cartesius. This is not allowed!**

.. code-block:: sh

	# go to the amuse directory
	> cd ~/amuse
	
	> mpiexec.hydra -n 1 nosetests -v tests

	
Running on a Cartesius node
---------------------------

Running on the Cartesius is typically done by submitting a slurm script. See the surfsare site for more info:

https://userinfo.surfsara.nl/systems/cartesius/usage/batch-usage
 
Below is a simple example script for running amuse on Cartesius.

.. code-block:: sh
	
	#!/bin/bash
	#SBATCH -N 1
	#SBATCH -n 1
	#SBATCH -p short
	#SBATCH -t 10

	cd ~/amuse

	mpiexec.hydra -n 1 ./amuse.sh examples/syllabus/gravity_simple.py
	
Submit using sbatch, get status using squeue, cancel using scancel
	
.. code-block:: sh

	#submit a script
	> sbatch example-script
	
	#list jobs of current user
	> squeue -u $USER

	#cancel job 505224
	scancel 505224

	#cancel all jobs of the current user
	> scancel -u $USER


Using multiple nodes should also work, by specifying this to slurm. MPI will automatically pickup on this and spread workers over all nodes.

.. code-block:: sh
	
	#!/bin/bash
	#SBATCH -N 2
	#SBATCH -n 10
	#SBATCH -p short
	#SBATCH -t 10

	cd ~/amuse

	#note that this simple example uses only a single worker
	mpiexec.hydra -n 1 ./amuse.sh examples/syllabus/gravity_simple.py
	
==================================
Hydrodynamics Interface Definition
==================================


Introduction
~~~~~~~~~~~~

In this chapter we describe the common interface for all 
hydrodynamics, grid based codes. For particle based, SPH codes see
the gravitational dynamics specifications.

Parameters
~~~~~~~~~~
All parameters have to be accessed with functions following
the template of the ``get_timestep`` and ``set_timestep`` functions.
A parameter access function may only retrieve or 
update the value of a single parameter. After all parameters have been set, the 
``commit_parameters`` function should be called,
this gives the code the opportunity prepare the model.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_timestep, get_timestep, commit_parameters
   

Grid Management
~~~~~~~~~~~~~~~~~
Most hydrodynamical codes work on grids or a hierarchy of grids. The following 
methods define the functionality to setup and query the grids. 

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_index_range_inclusive, set_boundary, setup_mesh, get_position_of_index,  get_index_of_position
    
Grid state
~~~~~~~~~~~~
Grid points in a hydrodynamics code have a well known, *minimal* state. This state is is defined
by a density, momentum density and energy density. The state can 
be retrieved and updated with the following functions.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_grid_state, get_grid_state
   

Grid State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a grid point can be transferred with the fill_grid_state and get_grid_state functions. To
support other properties (like pressure or MHD properties), the code can define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)


.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_density, set_grid_density, get_momentum_density, set_grid_momentum_density, get_energy_density, set_grid_energy_density


Model evolution
~~~~~~~~~~~~~~~
The hydrodynamics codes evolve the properties all grid cells in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: initialize_grid, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. The following 
function can be used to query the exact model time.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: get_time

External fields
~~~~~~~~~~~~~~~
Some Hydrodynamics codes support external acceleration or potential fields. The following functions
can be used to enter a gravitational potential field.

.. autoclass:: amuse.community.interface.hydro.HydrodynamicsInterface
   :members: set_has_external_gravitational_potential, get_has_external_gravitational_potential, get_index_range_for_potential, set_potential, get_potential



Glossary
--------

.. glossary::
    :sorted:
    
    
    Physical Domain
      A region in space, time and objects (stars / gas clouds) that can be 
      described by a set of physical models. A physical domain is a 
      limited description of the physics of a stellar system. For example 
      in the gravity domain the gravitational forces are taken into 
      account but not the evolution of stars or the hydrodynamics of th 
      interstellar medium
      
    Module
      Interface to a code in AMUSE. Also, file containing python definitions. 
      These definitions can be imported into the main Python program

    Unit test
      Automatic test case for a part of the code. For example one can 
      define a unit tests that checks the result of one method given a 
      specified argument.

    Application
      A script that combines modules into a multi-physical simulator.
      
    Run
      During a run the physical models are evolved to a wanted end state. 
      A run can end when a certain model time has passed. For example, one 
      can run a script until the age is 12 million years.
      
    Evolve
      Work out / develop the physical model equations. Usually, small steps 
      in time are taken, the new state of the system is calculated from the 
      previous state.

    Step
      An identifiable stage in a run or evolve. A code can take multiple 
      steps during one evolve step. A step can be the integrator step, 
      force calculation step.
      
    Data-model
      A single representation model for data of different physical codes.

    Set
      A collection of identifyable objects (for example all stars in a 
      globular cluster). Every value in the set can be seen as belonging
      to an object.
      
    Grid
      A collection of identifyable locations or regions (for example a box 
      of interstellar medium with a density, molecular composition). Every 
      value in the collection can be seen as belongin to a specific
      location (x,y,z point).

    Code
      A computer program to solve or approximate the equations of state of 
      a physical domain. Codes are usually written in Fortran or C/C++. 
      For example 'sse' is a code to perform stellar evolution.
      
    Production Code
      A code used for modelling stellar systems and described in published 
      articles.
      
    Toy Code
      A code with limited accurancy but that is easy to use.
      
    Numerical Experiment
      
    SPH, smoothed particle hydrodynamics
      A computational method used for simulating fluid flows

    AMR, adaptive mesh refinement
      
    Input / output format

    Binary star
      Two stars in close orbit
      
      
    Star Cluster
      A group of stars. The mean distance between the stars in smaller 
      than the mean distance between stars in the galaxy disks. Two types 
      of star clusters exist: globular clusters, open clusters

    Dense stellar systems
     
    Interstellar Medium
      The matter between the stars in a galaxy. 

    Interstellar Cloud
      Generic name for an accumulation of gas, plasma and dust. A denser 
      region of the insterstellar medium
      
    Molecular Cloud
      A type of insterstellar cloud whose density and size permits the
      formation of molecules,
     
    Gravity Domain

    Hydrodynamics

    Star  


    Stellar Evolution
      
    ZAMS, zero age main sequence
      Start of the main sequence in the evolution of a star. Stable stage
      after formation of the star from a collapsing gas cloud
      
    Main sequence star
      A star that derives its energy from the conversion of hydrogen into 
      helium in its core.
        
    Star cluster
      
      
    Metallicity
      Proportion of star matter made up of chemical elements other than 
      hydrogen and helium.
      
    Galaxy

    Star Core

    Star Envelope





      
      
===========================================
Simplified interface function specification
===========================================

In addition to the specification of remote interface functions by using 
the ``@legacy_function`` decorator, a simplified method is available, 
which can handle most cases of interface functions. Here we will 
describe the simplified interface function specification.

Example
~~~~~~~
Let's start with a simple example specification:

.. code-block:: python
             
        @legacy_function
        def sum():
            function = LegacyFunctionSpecification()
            function.addParameter('x', 'd', function.IN)
            function.addParameter('y', 'd', function.IN)
            function.addParameter('sum', 'd', function.OUT)
            function.result_type = 'i'
            function.can_handle_array = True
            return function

This can be converted to:

.. code-block:: python
             
        @remote_function(can_handle_array=True)
        def sum(x='d',y='d'):
            returns (sum='d')

As can be seen the parameters are specified in keyword/value style. 
``can_handle_array`` and ``must_handle_array`` are supplied, if 
necessary, as keywords to the decorator. A default integer return value 
(for the error code) is implied (but can be overridden, see below). The 
following table lists the options for the parameter specifications:

=========  ================  ===============================
data type  no default value  default value (for input) 
=========  ================  ===============================
boolean    "b", "bool"       True, False
integer    "i", "int32"      <int>, numpy.int32(<value>)
long       "l", "int64"      numpy.int64(<value>)
float      "f", "float32"    numpy.float32(<value>)
double     "d", "float64"    <float>, numpy.float64(<value>)
string     "s", "string"     "any other string"
=========  ================  ===============================

A unit specification can be added. Remember that parameters without 
default cannot follow parameters with default. So the following will 
generate an error:

.. code-block:: python
             
        @remote_function(can_handle_array=True)
        def sum(x=0.,y='d'):
            returns (sum='d')

Below are some more examples of valid specifications:

.. code-block:: python
             
        @remote_function
        def initialize_code():
            pass

        @remote_function
        def get_time():
            returns (time=0. | units.s)

        @remote_function(must_handle_array=True)
        def inout_sum(x='d',y='d', sum=0.):
            returns (sum=0.)

        @remote_function
        def function_with_float_return():
            return (__result=0.)

One limitation of this type of specification is that they won't
work if generated dynamically (so don't try, for example, to generate
a bunch of functions based on a list of parameter names).
Bridge 
======



Code
~~~~

.. automodule:: amuse.couple.bridge
    
    .. autoclass:: Bridge
        :members:
        
    .. autoclass:: AbstractCalculateFieldForCodes
        :members:
        
    .. autoclass:: CalculateFieldForCodes
        :members:
        
    .. autoclass:: CalculateFieldForCodesUsingReinitialize
        :members:
        
    .. autoclass:: CalculateFieldForCodesUsingRemove
        :members:
    
    .. autoclass:: CalculateFieldForParticles
        :members:

    .. autoclass:: GravityCodeInField
        :members:
        ==================
Distributed AMUSE
==================

It is possible to run AMUSE on multiple machines simultaneously. 
The AMUSE script itself always runs on a users' local machine, while workers for codes can be "send out" to remote machines such as workstations, clusters, etc.


Installation
------------

Deploying workers one remote machines requires a full installation of AMUSE on each machine. For each code "sockets" support needs to be present. This is done by default, and should be available for all codes. Note that older versions of Fortran lack the necessary features to support sockets workers.

On each machine, the distributed code also needs to be build. Distributed AMUSE requires a Java Development Kit (JDK), preferably Oracle Java version 7 or 8. The ``configure`` script tries to locate the JDK, but you may need to specify it by hand. For details, see:

.. code-block:: sh

	> ./configure --help

To build distributed amuse run the following at the amuse root:
	
.. code-block:: sh

	> make distributed.code

To check if the installation is set-up properly, run all the tests related to the worker interface:

.. code-block:: sh

	> cd $AMUSE_DIR
	> nosetests -v test/codes_tests/test*implementation.py
	
Note that Distributed AMUSE is mostly tested with the version of MPI includes in the amuse "prerequisites". 
If you get MPI errors while running remote (parallel) workers, try using the install.py script included in AMUSE to install the prerequisites.  

Overview
--------

Usage of Distributed Amuse is (by design) very close to the usage of any other code in AMUSE. 
The main difference being it contains resources, pilots, and jobs, instead of particles.

Resource
	Description of a machine capable of running jobs. 
    For each resource distributed AMUSE will launch a support process (HUB) to facilitate communication and coordination between the workers and AMUSE
	
Pilot
	Process running on a resource (often within a reservation on a resource) waiting for jobs to run on behalf of AMUSE. 
    Can consist of multiple machines.
	
Job
	Worker process.
    Will search for a suitable pilot to run on.

In general, a user will first define resources, then deploy pilots on these resources, and finally create codes that make use of the machines offered by the pilots.


Initializing the Distributed AMUSE code
---------------------------------------

Distributed Amuse can be initialized like any other code:

    >>> from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
    >>> from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
    >>> 
    >>> #initialize code, print output of code to console
    >>> instance = DistributedAmuse(redirection='none')


Parameters
----------

Distributed AMUSE supports a few parameters to adjust settings. 
All parameters need to be set before any resource, pilot or job is made to have effect.

Overview of settings:

debug
	Boolean parameters, defaults to False. 
    If true/enabled, will output additional debugging information and logs, both in the code output, and in a `distributed-amuse-logs` folder on any target machine used.
webinterface_port
	Port on which a simple webinterface is available for monitoring. 
    Defaults to "0", for a port determined automatically.
start_hubs
	To facilitate communication across different networks (with for instance firewalls), as hub is by default started on each resource. 
    This can be turned off if needed, for instance if all resources are within the same network.
worker_queue_timeout
	The user is responsible for making sure enough slots are available to run a worker. 
    If not, it will end up in the queue. 
    The time the worker will wait before giving up can be set using this parameter.
worker_startup_timeout
	The distributed code starts AMUSE workers running the actual codes. 
    This can take a while on some machines. 
    If needed, this parameter can be used to increase the time waited.
	
    >>> instance.parameters.debug = True
    >>> instance.parameters.webinterface_port = 5555
    >>> instance.commit_parameters()
    >>>
    >>> print instance.parameters.webinterface_port
    

Monitoring
----------

Distributed Amuse has a small build-in webinterface for monitoring. 
A utility function is available to get the url:

    >>> import webbrowser
    >>>
    >>> webbrowser.open(instance.get_webinterface_url())

Specifying resources
--------------------

In order to use a remote machine, AMUSE needs to have some information about this resource such as the host name, type of machine, username to gain access, etc.
This can be specified by creating a "Resource" in Distributed AMUSE. 
As a side effect, a communication hub is also started on the (frontend of) the resource.

    >>> resource = Resource()
    >>> resource.name = "some.resource"
    >>> resource.location = "user@machine.example.com"
    >>> resource.scheduler = "ssh"
    >>> resource.amuse_dir = "/home/user/amuse"
    >>>
    >>> instance.resources.add_resource(resource)

Overview of all options:

name
	Some user chosen name for the resource
location
	Address of the resource. Usually a hostname (e.g. somehost.somedomain.com). Could also be an IP address
amuse_dir
	Location of amuse on the remote machine (e.g. /home/user/amuse-svn)
tmp_dir
	Where all temporary files will be put on the remote machine
gateway
	Sometimes a machine is not reachable directly due to firewalls and such. Use this setting to provide an intermediate resource to route traffic via. This resource should already have been created.
scheduler_type
	The type of scheduler present on the remote machine. Defaults to 'ssh' useful for single machines. Current supported scheduler types: 'ssh', 'sge', 'slurm'
hub_queue_name
	Normally the support process is started on the front end. However, it can also be submitted to a queue by specifying it here.
hub_time_minutes
	When a hub is submitted, this option denotes the time the hub will be available.


Starting Pilots
---------------

The next step in running jobs remotely is to start a so-called pilot job on the resource specified previously. This pilot will submit a job to the resource, create necessary communication channels with the main amuse application, and wait for jobs to be started (currently mostly workers)

Note that pilots may not be started for a while. A function is available to wait until all created pilots have started.

    >>> pilot = Pilot()
    >>> pilot.resource_name='local'
    >>> pilot.node_count=1
    >>> pilot.time= 2|units.hour
    >>> pilot.slots_per_node=22
    >>> pilot.label='local'
    >>>
    >>> instance.pilots.add_pilot(pilot)
    >>> 
    >>> print "Pilots:"
    >>> print instance.pilots
    >>> 
    >>> print "Waiting for pilots"
    >>> instance.wait_for_pilots()

Overview of all options:

resource_name
	name of the resource to start the pilot on
queue_name
	queue to use to run the pilot (cluster specific, not used in case of ssh)
node_count
	number of nodes to start the pilot on
time
	time to keep the pilot active
slots_per_node
	number of workers to start on a node. Usually the number of cores, but could be less if memory is a limiting factor, or workers are multi-core capable
label
	label to attach to the pilot. Can be used when starting workers to run workers on specific pilots
options
	Additional options. Usually not required.


Starting jobs
-------------

When running remote workers, they can be started as normal. 
However, AMUSE needs to be signalled to use the distributed code to start them instead of the normal process.
A function is available to enable and disable this.

    >>> print "starting all workers using the distributed code"
    >>> instance.use_for_all_workers()

    >>> print "not using distributed workers any longer"
    >>> instance.use_for_all_workers(enable=False)

Alternatively, you can also explicitly enable the distributed code per worker

    >>> print "using this distributed instance for all distributed workers"
    >>> instance.use_for_all_distributed_workers(enable=True)
    >>> worker = Hermite(channel_type='distributed')

Or, even pass the instance of the distributed code you would like to use, in the rare case you have multiple distributed codes

    >>> worker = Hermite(channel_type='distributed', distributed_instance=instance)

Worker options
--------------

This section lists all the relevant worker options for Distributed AMUSE. 
Most are new, some are also supported in the other channel implementations.
You are normally not required to use any options.

number_of_workers
	Number of worker processes started (thus working as normally the case). 
    Each worker takes up a slot of the pilot (see above)
label
	Label of the pilot to use. By default any pilot with enough free slots found will be used to start this worker. 
    Using the labels an explicit selection can be done.
number_of_threads
	Number of threads used in the process. 
    This can be used to explicitly set the OMP_NUM_THREADS environment variable in the worker
channel_type
	Set this to "distributed" to start workers using the distributed code. 
    Alternatively, use the use_for_all_workers functions as described above to set this by default
distributed_instance
	This is a reference to the distributed instance used to start the worker, in the rare case you have multiple distributed codes.
dynamic_python_code
	Boolean option stating if this code is a dynamic python code. 
    If so, all .py files in the worker directory will be copied to the remote machine before starting the code.


Labels
------

By default workers are started on any available pilot with enough slots available. 
However, sometimes you would like to have more control over which worker is started where, for instance if special hardware is present on some machines.

The concept of labels can be used within Distributed AMUSE to get this functionality.
If a label is attached to a worker (one of the parameters when starting a worker, see above), only pilots with exactly the same label (specified when the pilot is started) are considered candidates for running the worker. 
The name of labels is completely up to the user.

For instance, say a simulation uses a number of workers running on a CPU, and a single GPU worker.
The following code will put all the cpu workers on one machine, and the single gpu worker on another.

    >>> cpu_pilot = Pilot()
    >>> cpu_pilot.resource_name='machine1'
    >>> cpu_pilot.node_count=1
    >>> cpu_pilot.time= 2|units.hour
    >>> cpu_pilot.slots_per_node=30
    >>> cpu_pilot.label='CPU'
    >>> instance.pilots.add_pilot(cpu_pilot)
    >>>
    >>> gpu_pilot = Pilot()
    >>> gpu_pilot.resource_name='machine2'
    >>> gpu_pilot.node_count=1
    >>> gpu_pilot.time= 2|units.hour
    >>> gpu_pilot.slots_per_node=1
    >>> gpu_pilot.label='GPU'
    >>> instance.pilots.add_pilot(gpu_pilot)
    >>>
    >>> ...
    >>> worker1 = Hermite(label='CPU')
    >>> worker2 = Bonsai(label='GPU')
    >>>
    >>> #will not start due to a lack of slots.
    >>> worker3 = Bonsai(label='GPU')
 

Examples
--------

AMUSE contains a number of examples for the distributed code. See examples/applications/

Gateways
--------

Gateways can be used in case of connectivity problems between machines, such as firewalls and private IP addresses. 
This is for instance the case at the LGM. 
A gateway is started like any other resource (and thus require a valid installation of AMUSE on each gateway). 
This resource can then be specified to be a "gateway" to another resource. 
In this case all ssh connections will be made via the gateway, so make sure you can login from the gateway to the target machine without using a password, as well as from your local machine.

Commonly Encountered Problems
-----------------------------

Most initial setup problems with the Distributed AMUSE code can be solved by checking:

- Can you login to each machine you plan to use using ssh without using a password? 
  See for instance here on how to set this up: https://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/
- Did you configure a Java JDK version 1.7 or higher using ./configure? 
  Check the content of config.mk to see which java is used, and what version was detected. 
  Make sure to do a "make clean" and "make" in case you make any changes. This should also be done on all machines.
- Is AMUSE configured properly on each and every machine? 
  Running the code implementation tests is a good way of spotting issues:

    >>> nosetests -v test/codes_tests/test_*_implementation.py

- Are the settings provided for each resource correct (username, amuse location, etc)
- Have you set the correct mpiexec in ./configure? This setting is normally not used by AMUSE, so you may only now notice it is misconfigured

In case this does not help, it is probably best to check the output for any errors. 
Normally worker output is discarded by most scripts. 
Use 'redirect=none' to see the output of the workers, a lot of errors show up in this output only. 
There is also a "debug" parameter in Distributed Amuse.
If enabled, output for each pilot will be in a "distributed-amuse-logs" folder in the home of each remote machine used, and additional information is printed to the log from the local AMUSE script.

=======================
Sets and Grids in Codes
=======================

.. automodule:: amuse.datamodel.incode_storage
    :members:
    :undoc-members:
    :show-inheritance:
.. _options-label:

Options
=======

.. automodule:: amuse.support.options

    Introduction
    ------------
    The AMUSE framework provides a generic way to handle optional
    attributes for any class. When an optional attribute is requested
    the class can read the value from an configuration file or
    return a default value when nothing is found in the file. 
    The values for all options are stored in one configuration file. 
    
    Configuration file location
    ----------------------------
    The AMUSE framework will search for the configuration file in 
    following locations:

    1. First, the framework tries to find a file named **amuserc** 
       in the current working directory
    2. Second, the framework tries to find a hidden file named 
       **.amuserc** in the home directory.
    3. Last, the framework tries to find a file named **amuserc**
       in the amuse installation directory
    
    When no file is found all options will refer to the default-values.
    
    
    Configuration file format
    -------------------------
    The configuration file is formatted similar to Microsoft Windows 
    INI files. The configuration file consists of sections, led by a 
    [section] header and followed by name: value entries.
        
    For example::
        
        [channel]
        redirection=none
        debugger=xterm
    
    
    Option values
    -------------
    Values for optional attributes are determined in four
    different ways:
    
    * the attribute is set on an instance of the object::
        
        channel.debugger = "ddd"
            
    * given when creating an object of the class::
    
        channel = MessageChannnel(debugger = "ddd")
     
    * given in the ini file::
    
        [channel]
        debugger = ddd
    
    * set as default value::
    
        @option
        def debugger(self):
            return "ddd"
    
    
    Sections lookup
    ---------------
    Options can be set in different sections. To determine
    the value of an option the framework first searches
    the sections of a class and then the sections of the option.
    
    
    Use
    ---
    
    .. autoclass:: option
        :members:
        
    .. autoclass:: OptionalAttributes
        :members:
    
    
    
=====================================
Stellar Dynamics Interface Definition
=====================================

Introduction
~~~~~~~~~~~~

In this chapter we describe the common interface for all stellar dynamics codes.

Parameters
~~~~~~~~~~
Gravity dynamics codes have at least one specified parameter. Other parameters need
to be specified on a per code basis. All parameters have to be accessed with functions following
the template of the ``get_eps`` and ``set_eps`` functions. A parameter access function may only
retrieve or update the value of a single parameter. After all parameters have been set, the 
``initialize_code`` function should be called, this gives the code the opportunity prepare the
model.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_eps2, set_eps2, initialize_code
   

Object Management
~~~~~~~~~~~~~~~~~
Most gravitational dynamics codes work on particles (stars, black holes or gas). The following 
methods define the functionality to create, remove and query the particles in the code. *Currently 
the interface does not handle different types of particles*

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: new_particle, delete_particle, get_number_of_particles, get_index_of_first_particle, get_index_of_next_particle

    
Object state
~~~~~~~~~~~~
Particles in gravitational dynamics have a well known, *minimal* state. This state is is defined
by a location, velocity and mass and radius. The state can be retrieved and updated with 
the following functions.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_state, set_state
   

Object State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a particle can be transferred with the get_state and set_state functions. To
support other properties (like acceleration), the code can define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)


.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_mass, set_mass, get_position, set_position, set_acceleration, get_acceleration, get_potential


Model evolution
~~~~~~~~~~~~~~~
The gravitational dynamics codes evolve the properties of the particles in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: commit_particles, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. The following 
function can be used to query the exact model time, the total energies and colliding particles.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_time, get_kinetic_energy, get_potential_energy, get_indices_of_colliding_particles, get_center_of_mass_velocity, get_center_of_mass_position, get_total_mass, get_total_radius

Services
~~~~~~~~
Some Gravitational Dynamics codes can provide services for other codes. Currently calculating
the gravity at a given point is the only specified function.

.. autoclass:: amuse.community.interface.gd.GravitationalDynamicsInterface
   :members: get_gravity_at_point
==================
AMUSE Style Guide
==================


This Style Guide covers the source code written in the AMUSE project. The
source of the existing codes, integrated into AMUSE, do not have to
follow this guide. Existing codes are not rewritten.

Python Code
-----------
All python code should be consistent with the 
`Python Style Guide <https://www.python.org/dev/peps/pep-0008/>`_.
We have defined some deviations from the python style guide. 
These are listed in this document.

Maximum Line Length
    Limit *most* lines to a maximum of 79 characters.
    
    Some lines may be longer than 79 characters, especially when
    this increased readability. When you need to define a class with
    a long name, having a superclass with a long name, a line with
    more than 79 characters is OK.
    
    For example::
    
        class NeedsALongNameToDescribeItsFunction (AComplicatedNameForTheSuperClass):
        
    Is better than::
    
        class NeedsALongNameToDescribeItsFunction \
            (AComplicatedNameForTheSuperClass):
    
    Or::
    
        class NeedsALongNameToDescribeItsFunction (
            AComplicatedNameForTheSuperClass):
    
    
Naming Conventions
    Naming the classes and methods in the code is very important. However,
    do not spend too much time on naming. Use long names. Use a longer
    name when a shorter name cannot be determined quickly. When you've 
    defined a very long name, it's use may be harder. When it is used 
    a lot, a refactoring to shorter name will present itself. When it is
    not used a lot, at least the class or method is described better
    in the longer name.
    
    Do not use abbreviations in names. Especially, do not use the 
    first three or four letters of a word (for example, *dir* for directory
    or *ref* for reference). Also, do not remove the vowels to shorten a
    name (for example *tmp* for temp).
    
    The following one letter variables are allowed.
    
     * Loop variables, when used as ``for x in ...``,
       names allowed are 'x' and 'y'. Example::
       
            for x in stars:
                print x.mass
        
     * Loop index variables, when used as ``for i in range(10):``, 
       only name allowed is 'i'. Thus when writing a loop in a loop, you need
       to use longer names.
       For example, to iterate over every item in a 10 by 5 table, do:s:
            
            for row in range(10):
                for column in range(5):
                    print table[row][column]
        
    These names are only allowed in short loops of maximum 6 lines. Longer
    loops need better names. Also, for looping over row indices and column
    indices use *row* and *column*.
    
    
C / C++ code
------------

All C and C++ code should be consistent with the 
`Google C++ style guid <https://google.github.io/styleguide/cppguide.html>`_

File Names
    All C++ files should have a ``.cc`` suffix.
     
Naming Conventions
    We follow the Python style for naming function names. Function names \
    are all lowercase, with underscores between words.
    

Fortran code
------------
All Fortran code written in AMUSE should be Fortran 90.  We will follow
the coding standard given in 
`Fortran Coding Standard for the Community Atmospheric Model <http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html>`_

File Names
    All fortran 90 files should have a ``.f90`` suffix.
    
Layout
    Source code should follow the "free form" fortran 90. Do not use
    the "fixed form" fortran 77 layout.
    
Indentation
    Use indentation to highlight the block structure of your
    code::
    
        FUNCTION example(arg1)
            REAL  :: arg1, example
            example = arg1 * 2.0
        END FUNCTION
        
    Use 4 character to indent the code.
    

    
    
    

    
    
==================
From Codes to Data
==================

Introduction
------------

The framework contains two distinct levels on which interaction with the codes takes place.
The first level interacts directly with the code and is implemented as a set of functions
on a class. The second level is built on top of the first level and abstracts the function calls
to handling objects and datasets. Often multiple function-calls on the first level can be 
abstracted to a single statement (assignment or query) on the second level.

The first level
---------------
The first level is a direct interface to the code, in this chapter the kind of functions
supported by the first level interface code will be briefly described. All functions that can
be defined on the first level fall in two categories, those that handle scalars (a single
value for each parameter) or those that handle 1-D vectors (a list of values for each
parameters). For functions that handle 1-D vectors each vector must be of the same length. 

.. note:: 

    Not every function will fit in the two categories, but it is usually possible
    to rewrite a function or create a interfacing function that do fit into
    one of the two categories. Supporting only these two categories keeps the 
    communication layer simpler and allows for some optimizations in the 
    communication between python and C/Fortran codes.

An example of using the first level with scalars and vectors::

    from amuse.community.codes.athena.interface import AthenaInterface
    
    # create an instance of the code (will start an application 
    # in the background to handle all requests)
    hydro = AthenaInterface()
    
    # set parameters needed by the code
    # these are functions handling one scalar input 
    # parameter
    hydro.set_gamma(1.6666666666666667)
    hydro.set_courant_friedrichs_lewy_number(0.8)
    
    # define a grid having 5 cells in all directions and 
    # with the total length of the grid in each direction
    # is 1.0
    # this is a function handling multiple
    # scalar input parameters
    hydro.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
    
    # setup boundary conditions 
    # (can be periodic, reflective, outflow)
    hydro.set_boundary(
        "periodic","periodic",
        "periodic","periodic",
        "periodic","periodic"
    )
    
    # let the code do some work
    # (athena will allocate the grid)
    # this is a function handling no 
    # scalar input parameters and having 
    # 1 scalar output parameter
    hydro.commit_parameters()
    
    # lets print the center position
    # of one grid point
    print hydro.get_position_of_index(1,2,3)
    
    # all calls so far have been to functions handling scalar values
    # the next calls will be to functions handling vectors of values
    
    # lets print the center positions
    # of all grid points on one line
    print hydro.get_position_of_index(range(0,5), [0] * 5, [0] * 5)

In the previous example we used functions with scalar parameters and
vector parameters. The functions handling vectors often can also
handle scalars, the framework will take care of the necessary 
conversions.

All first level functions are not actual python functions, these
functions are instances of a special Python class that
implements function call handling. To continue our example::

    # let's take a look at the kind of functions
    # on the first level
    print hydro.get_position_of_index
    
    # you can ask the specificition of a
    # first level function
    print hydro.get_position_of_index.specification


Adding units
~~~~~~~~~~~~
The first step after defining a first level function is to
specify the units of the in- and output-parameters of the first
level function. 

    
============
Introduction
============

In this document we will specify the low-level interfaces of the comminity codes. 
These interfaces should follow a pattern. This pattern is described bellow. This
pattern can be used to keep the interfaces of the codes consistent
between codes of different physical domains and to help in learning about
the functionality of the interfaces. In later chapters we define the specific
interfaces for the modules of different physical domains. The actual interface of
a code should derive from these interfaces.

These interfaces only specify the low-level interfaces. The low-level interfaces
do not support unit handling, data conversion and attributes. This functionality
is handled by the next-level interfaces described in nextlevel_.

Data Types
----------
The exact size (in bytes) and type of the values sent between the python 
layer and the community codes is very important. In this document we will 
specify the type of each argument and return value. 
Currently AMUSE supports 4 types, these are described in the following table. 

========== ======== ================ ================
Type name  Size     Fortran          C               
..         bytes    type             type             
========== ======== ================ ================ 
int32      4        integer          long             
float64    8        double precision double
float32    4        real             float
string     n                         char *
========== ======== ================ ================ 

Function template
------------------
All functions in the interface should follow the same template.
Each function returns an error or status code. Results are returned through 
the function arguments. In C these arguments need to be pointers
to valid memory locations.

.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: example_function
        
The error codes all have the same general form. Zero stands for no error, a negative
value indicates some error happened, a positive value is returned when the function
ends in a special state.

0 - OK
    Function encountered no error or special state
<0 - ERROR
    Something went wrong in the execution of the function
>0 - STATE
    Function has encountered an expected special state. For example the code
    has detected a collision between two stars.
    

Function categories
-------------------

Parameters
~~~~~~~~~~
Codes can have a number of parameters. Some parameters are domain specific, these
parameters are found for all codes in a specific domain (for example the smoothing
length in gravitational dynamics). Other parameters are only defined and used by
a specific code. The domain specific parameters are defined on the domain specific
interfaces, when a code supports the parameter it should implement the
specified functions. Other parameters have to be accessed with functions following
the template of the :meth:`~ amuse.community.interface.example.ExampleInterface.get_example_parameter` 
and :meth:`~amuse.community.interface.example.ExampleInterface.set_example_parameter` functions. 

.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: get_example_parameter, set_example_parameter
   

A function used to access (set or get) a parameter may only retrieve
or update the value of a single parameter. Functions setting two or more
parameters in one go are not supported by the next-level interfaces. After all 
parameters have been set, the  :meth:`~amuse.community.interface.example.ExampleInterface.initialize_code`
function should be called, this gives the code the opportunity prepare the model.


.. autoclass:: amuse.community.interface.example.ExampleInterface
   :members: initialize_code

   
Object Management
~~~~~~~~~~~~~~~~~
Codes can work on particles or grids (stars, black holes or gas). The methods 
in the *Object Management* category define the functionality to create, remove 
and query the particles or gridpoints in the codes. 

When a code supports objects, the code is responsible for managing these objects. 
The code needs to assign a unique index to a particle so that the particle can be 
referred to in other function calls. This is a *major* difference with the 
``MUSE`` code. Where the user of the code was responsible for assigning unique ids
to the particles. This change makes the implementation of the code simpler and
allows the code to support creation of new objects during simulation. For example
a hydrocode can add or delete gridpoints during the evolution of the model.

    
Object state
~~~~~~~~~~~~
Particles in the same physical domain can have a well known, *minimal* state. 
For example, in the gravitational dynamics domain the state of a particle can be
defined by a location vector, velocity vector, a mass and a radius. The methods
in the *Object State* category provide a way to access this state in one function.
  

Object State, Extension Mechanism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Not all information of a particle can be transfered with the get_state
and set_state functions. Some codes may support other properties of a particle,
the code can define ``get_`` and ``set_`` functions for these properties. These 
functions follow the pattern defined in the *Parameters* category. The functions
must either get or set a scalar property (1 argument) or
a vector property (3 arguments).

Model evolution
~~~~~~~~~~~~~~~
The main function of a code is often evolving a model in time or solving a steady
state solution. The methods that control model evolution or start and stop the 
model calculations all belong to the *Model evolution* category. At this time, 
no pattern is defined for the functions in this category.

Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model
calculations. All the functions in this category follow 
the 'get_name' pattern. The state of code should not change during a function 
call to a function in this category. The functions must either get
a scalar property (1 argument) or a vector property (3 arguments).

Services
~~~~~~~~
Some codes can provide services for other codes in the same or other physical
domains. For example, gravitational dynamics code might provide a function to
calculate the gravity force at a point. The methods that provide these
services all belong to this category.
====================
Directory Structure
====================

The amuse source-code is separated into 3 directories:

* ``src`` - source code, implementation of the environment.
* ``test`` - applications, examples and unittests.
* ``support`` - build system, test system.

Under the ``src`` directories all code needed to run AMUSE can be found.
One can view this code as an *library* that can be used to create
*applications* to do numerical astrophysical experiments. This code will
contain the building blocks needed to interface with codes,
import and export data, do unit conversions, and all other AMUSE
functionality.

Under the ``test`` directories all application and test code can be
found. This directory tree will contain scripts to do a complete
astrophysical experiment. Also all unit-tests can be found here. These
unit tests each cover only a small part (unit) of the functionality of
AMUSE. For example a  test to check the import of a file to AMUSE data
format.

Under the ``support`` directories all support code for the building
system can be found.

The ``src`` directories
~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``src`` directory are further split into:


* ``community`` - contains the source code of existing astrophysical
  applications and *glue* code to the AMUSE interface classes. In other
  words this directory contains the implementation of the interfaces.
* ``support`` - contains the AMUSE generic code, defines the data
  representation and input/output routines and also provides the generic
  unit handling code. Code in the ``interface`` and ``community`` directories
  use these functions and classes to provide their functionality.


The ``test`` directories
~~~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``test`` directory are further split into:

* ``unit_tests`` - All unit testing code. These tests are coded using
  the standard unit testing framework that is included in the Python
  distribution (``unittest``). See python module documentation for further
  information: http://docs.python.org/library/unittest.html.
* ``application`` - contains the source code of published applications.
* ``examples`` - contains documented example codes.

The ``support`` directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The directories under the ``support`` directory are further split into:

* ``test`` - Scripts to support the testing of AMUSE code.
* ``build`` - Scripts used by the building system of AMUSE.



==========================
Input / Output Framework
==========================


Introduction
------------
The reading and writing of the files is done by subclasses
of the :class:`FileFormatProcessor` class. 

    
Extending
---------

.. autoclass:: amuse.io.FileFormatProcessor
    :members:
    
.. autoclass:: amuse.io.FullTextFileFormatProcessor
    :members:
        
        
        
        
====================
Stopping Conditions
====================

Introduction
------------
Codes in AMUSE evolve their models in three loops, the "inner", the 
"outer" loop, and the "script" loop. The "inner" loop is controlled 
by the code and evolves the model in a sufficiently small steps to 
limit errors and be able to simulate important physics. The "outer" 
loop invokes the "inner" loop until a condition is met specified by 
the AMUSE script. The AMUSE script (or the "script" loop) interacts 
with the "outer" loop.

 .. graphviz::

    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow, width=1.5];
      subgraph cluster0 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="AMUSE Script";
            "script loop";
      }
      subgraph cluster1 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Community Code";
            "outer loop";
            "inner loop";
      }
      "script loop" -> "outer loop"[minlen=2];
      "outer loop" -> "inner loop";

      "script loop" -> "script loop"  [dir="back", constraint=false,minlen=2, headport="se", tailport="ne"];
      "outer loop" -> "outer loop"[dir="back", constraint=false, headport="se", tailport="ne"];
      "inner loop" -> "inner loop"[dir="back", constraint=false, headport="se", tailport="ne"];
    }

In AMUSE the outer loop is always limited by the model time. Every
:func:`evolve_model` call has one argument, the model time. When the
simulation is about to go beyond this time (or is at this time),
control is returned to the caller.
This process is shown in the following code example.

.. code-block:: c

    // example outer loop
    void outer_loop(double end_time)
    {
        while(model_time < end_time)
        {
            ...
            inner_loop()
            ...
            model_time += delta_time_for_this_step
        }
        return
    }

The model time is often enough to control the loop but some
situations call for an earlier break in the loop when the model
time has not reached the end time. To be able to break the outer loop
before the end time is reached, the framework supports
**stopping conditions**. The following code example
shows how these conditions might work in a C code.

.. code-block:: c

    // example outer loop with stopping conditions
    void outer_loop(double end_time)
    {
        int special_condition_was_met = 0;

        while(model_time < end_time)
        {
            ...
            inner_loop()
            ...
            check_if_special_conditions_are_met();
            ...
            model_time += delta_time_for_this_step
            ...
            if(special_condition_was_met) {
                break;
            }
        }
        return
    }

For example, a code might be able to detect a collision of two stars 
and the user wants to handle this event. If the code would proceed 
to the end time, the collision is long over and the user has no way 
to know about the collision and cannot handle it.

.. image:: stopping_conditions_collision.png

The **stopping conditions** framework in AMUSE provides 
functionality for using and implementing stopping conditions. This 
framework uses simple tables in the codes and translates these 
tables to objects in the python scripts. In the next section the use 
of this framework is described. In the final section the 
implementation of the framework in a community code is described.

Using stopping conditions
-------------------------

All interaction with the *stopping conditions* of a *code*
is via the :attr:`stopping_conditions` attribute of a :class:`CodeInterface`
class. When the :attr:`stopping_conditions` attribute is
queried it returns a :class:`StoppingConditions` object. All
:class:`StoppingConditions` objects have an attribute for each
stopping condition supported by AMUSE. The
attributes have the following names:

collision_detection
    The outer loop breaks when two stars connect
pair_detection
    The outer loop breaks when two stars pair, possibly
    creating a binary star.
escaper_detection
    The outer loop breaks when a star escapes the
    simulation (is no longer gravitationally bound
    to the other stars)
timeout_detection
    The outer loop breaks when a certain amount
    of computer (wall-clock) time has elapsed.
number_of_steps_detection
    The outer loop breaks when a certain amount
    of steps (iterations) is reached.

.. note ::

    :class:`CodeInterface` objects provide a object oriented
    interface to each code.


All :class:`StoppingConditions` objects provide a string
representation describing the supported conditions and their
states.

.. code-block:: python

    >>> from amuse.lab import *
    >>> code = Hermite()
    >>> print code.stopping_conditions
    Stopping conditions of a 'Hermite' object
    * supported conditions: collision_detection, pair_detection, timeout_detection
    * enabled conditions: none
    * set conditions: none
    >>> code.stop()

For every attribute the user can determine if it is
supported, if it was enabled or if it was hit during
the :func:`evolve` loop.

.. code-block:: python

    >>> from amuse.lab import *
    >>> code = Hermite()
    >>> print code.stopping_conditions.collision_detection.is_supported()
    True
    >>> print code.stopping_conditions.collision_detection.is_enabled()
    False
    >>> print code.stopping_conditions.collision_detection.is_set()
    False
    >>> code.stop()

When a stopping condition was hit during the inner or outer loop 
evaluation the user can query which particles were involved. The 
particles involved are stored in columns. For pair wise detections 
(collection and pair) the condition provides two columns. The first 
column contains all the first particles of the pairs and the second 
column contains all the second particles of the pairs. These columns 
can be queried with the :meth:`particles` method on a stopping 
condition attribute. This method takes one argument, the column 
index and returns a particle subset with the involved particles.

========  ============  ============
Pair      particles(0)  particles(1)
========  ============  ============
1 + 2     1              2
5 + 11    5              11
12 + 10   12             10
========  ============  ============

.. code-block:: python

    >>> from amuse.lab import *
    >>> code = Hermite()
    >>> sc = code.stopping_conditions.collision_detection
    >>> sc.enable()
    >>> code.evolve_model(0.1 | nbody_system.time)
    >>> print sc.is_set()
    True
    >>> print sc.particles(0)
    ...
    >>> pair = sc.particles(0)[0], sc.particles(1)[0]
    >>> code.stop()

Implementing stopping conditions
--------------------------------
All codes supporting stopping conditions must implement the 
:class:`StoppingConditionInterface` interface. This can be easily 
done by inheriting from the :class:`StoppingConditionInterface` and 
implementing all the functions defined in the interface:

.. code-block:: python

    class MyCodeInterface(CodeInterface, StoppingConditionInterface):
        ...

This :class:`StoppingConditionInterface` interface models the interface to a
simple data model, best described as two tables:

Defined stopping conditions table
    This table list the types of stopping conditions
    and if these are supported by the code, and if so, if
    these are enabled by the user.

    ======  ========= ========
    type    supported enabled
    int     int       int
            read only read/write
    ======  ========= ========
    0       0         0
    1       1         0
    2       0         0
    3       1         0
    ======  ========= ========

    type
        The type of the stopping condition, integer between
        0 and N. Defined by the framework
    supported
        1 if the stopping condition is supported, 0 otherwise.
        Needs to be set by the implementer
    enabled
        1 if enabled by the user, 0 otherwise.
        Implementer must provide a mechanism to set and unset
        this field.

Set stopping conditions table
    This table has a row for each stopping condition
    set during the outer and inner loop evaluations.
    Every row lists the type of condition set and which
    particles were involved (if any). The id of the
    last particle in the list must be -1.

    ======  =========== =========== ======= ============
    type    particle[0] particle(1) ...     particle(n)
    int     int         int         int     int
    ======  =========== =========== ======= ============
    1       1           2           -1      -1
    1       10          12          -1      -1
    3       -1          ...
    1       11          20          ...
    ...
    3       -1
    ======  =========== =========== ======= ============

    type
        The type of the stopping condition, integer between
        0 and N. Defined by the framework. Set by the code
    particle(0...n)
        Index of the particle involved in the stopping condition.
        The index of the last particle must be -1

.. autoclass:: amuse.support.codes.stopping_conditions.StoppingConditionInterface
    :members:

This interface contains a lot of functions and implementing it
might be a time consuming task. To help implementing the
interface, AMUSE provides a C library that does most of the
table management, this library is described in the next section.

Using the stopping conditions library
-------------------------------------
The stopping conditions tables are implemented as a
static C library. This library is part of the AMUSE distribution
and can be found in the ```lib/stopcond``` directory.

The library implements *all* functions of the :class:`StoppingConditionInterface`
interface and provides a number functions for the
code implementer to actually set and support
stopping conditions. The library does not implement any
stopping condition detection routines. These have to
be implemented on a per code basis. The functionality
of the library is limited to managing the stopping conditions
tables and the library has to be controlled by the code
so that stopping conditions are set.

The library implements the defined stopping conditions table as
three long integers. These integers are used as bitmaps (the state
of individual bits is used to determine if a condition is true or
false). The individual bits can be tested with the following
macros:

.. data:: COLLISION_DETECTION_BITMAP

.. data:: PAIR_DETECTION_BITMAP

.. data:: ESCAPER_DETECTION_BITMAP

.. data:: TIMEOUT_DETECTION_BITMAP

.. data:: NUMBER_OF_STEPS_BITMAP

The three variables are not accessible directly from FORTRAN
codes. The stopping condition library has special functions
for FORTRAN to access them. 

The three global variables are:

.. c:var:: long enabled_conditions

Individual bits are set when the user enables a stopping condition.
Access in FORTRAN through:

.. code-block:: fortran

    INCLUDE "../../lib/stopcond/stopcond.inc"
    INTEGER:: is_stopping_condition_enabled
    INTEGER:: is_collision_detection_enabled
    INTEGER:: error
    error = is_stopping_condition_enabled(COLLISION_DETECTION, is_collision_detection_enabled)

.. c:var:: long set_conditions

Individual bits are set when a condition is detected by the
code. (the :c:func:`set_stopping_condition_info` will fill
this bitmap). Access in FORTRAN:

.. code-block:: fortran

    INCLUDE "../../lib/stopcond/stopcond.inc"
    INTEGER:: is_stopping_condition_set
    INTEGER:: is_collision_detection_set
    INTEGER:: error
    error = is_stopping_condition_set(COLLISION_DETECTION, is_collision_detection_set)

.. c:var:: long supported_conditions

.. code-block:: fortran

     set_supported_conditions()
 
Individual bits are set for supported conditions. This global
must be defined by the code. 

For a code that support
collision detection and timeout detection the variable can
be defined in c as:

    .. code-block:: c

        set_support_for_condition(COLLISION_DETECTION);
        set_support_for_condition(PAIR_DETECTION);

or in FORTRAN as:
    .. code-block:: fortran

        INCLUDE "../../lib/stopcond/stopcond.inc"
        
        set_support_for_condition(COLLISION_DETECTION)
        set_support_for_condition(PAIR_DETECTION)

To refer to the type of a stopping condition, the library defines the
following macros:

.. data:: COLLISION_DETECTION

.. data:: PAIR_DETECTION

.. data:: ESCAPER_DETECTION

.. data:: TIMEOUT_DETECTION

.. data:: NUMBER_OF_STEPS_DETECTION

The following four functions can be used to set and reset
the stopping conditions.

.. c:function:: int reset_stopping_conditions()

    Resets all stopping conditions. Needs to be called just
    before entering the outer_loop

    .. code-block:: c

        void outer_loop(double end_time)
        {
            reset_stopping_conditions()
            while(model_time < end_time)
            {
                ...
            }
        }

.. c:function:: int next_index_for_stopping_condition()

    Reserves a row on the stopping conditions table.
    Needs to be called after detecting a stopping condition so
    that the code can fill in the stopping condition information.

    .. code-block:: c

        if(COLLISION_DETECTION_BITMAP & enabled_conditions) {
            int row = next_index_for_stopping_condition();
            set_stopping_condition_info(row, COLLISION_DETECTION);
            ...
        }

.. c:function:: int set_stopping_condition_info(int index, int type)

    Sets the type of stopping condition encountered for
    the given row index. Use one of the predefined type macros
    for the type argument

    .. code-block:: c

        set_stopping_condition_info(row, COLLISION_DETECTION);


.. c:function:: int set_stopping_condition_particle_index(int index, int column_index, int index_of_particle)

    Sets the id of a particles involved in the stopping
    condition (only needed for stopping conditions with particles)

    .. code-block:: c

        set_stopping_condition_info(row, COLLISION_DETECTION);
        set_stopping_condition_particle_index(stopping_index, 0, ident[i]);
        set_stopping_condition_particle_index(stopping_index, 1, ident[j]);

Stopping condtions with parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *timeout* and *number of steps* stopping conditions are parameter
based, i.e. the conditions are met if the outer loop took a certain time
or a reached a certain number of steps (iterations). The parameters are
accessed through the following procedures:

.. c:function:: int set_stopping_condition_timeout_parameter(double value); 

    Sets the computer (wall-clock) time parameter.

.. c:function:: get_stopping_condition_timeout_parameter(double *value);

    Gets the computer (wall-clock) time parameter.

.. c:function:: int set_stopping_condition_number_of_steps_parameter(int value);

    Sets the number of steps parameter.

.. c:function:: int get_stopping_condition_number_of_steps_parameter(int *value);

    Gets the number of steps parameter.

For MPI enabled codes, the stopping conditions must be distributed
and collected. The following three functions can be used to
manage the MPI communication.

.. c:function:: int mpi_setup_stopping_conditions()

    Setups some globals for MPI communication, this needs
    to be called once per code execution (for example
    in a :func:`initialize_code` function)

.. c:function:: int mpi_distribute_stopping_conditions()

    Distributes the state of the conditions, so each process
    knows which conditions are enabled by the user. Needs
    to be called when the parallel code is entered (and all
    other data is distributed).


.. c:function:: int mpi_collect_stopping_conditions()

    Collects the individual tables of the set conditions and
    combines these to one table on the root process (MPI rank
    == 0). Needs to be called at the end of the parallel part
    of the code.

The procedure for adding stopping conditions is explained by
the following steps.

1.  Add compiler and linker flags to the build script or
    Makefile.

    To use the ```stopcond``` library, the compile and link steps
    of the community code must be changed. The following compiler flags
    need to be added::

        -I$(AMUSE_DIR)/lib/stopcond

    For the link step the following linker flags are needed::

        -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi

2.  Include the header file and defining the globals

    The name of the header file is "stopcond.h" and
    this header file needs to be included in all
    files referring to the stopping conditions.

    .. code-block:: c

        #include <stopcond.h>

    The code should take care of the setting
    of the :c:data:`supported_conditions` bitmap
    in its initialization, e.g.:

    .. code-block:: c

        set_support_for_condition(COLLISION_DETECTION);
        set_support_for_condition(PAIR_DETECTION);

    In FORTRAN include "stopcond.inc" in subroutines that
    require the constants.
    
    .. code-block:: fortran
    
        INCLUDE "../../lib/stopcond/stopcond.inc"
        set_support_for_condition(COLLISION_DETECTION)
        set_support_for_condition(PAIR_DETECTION)

3.  Check for a condition when it is enabled and set the condition
    if found. In c:

    .. code-block:: c

        int is_collision_detection_enabled;
        is_stopping_condition_enabled(COLLISION_DETECTION, &is_collision_detection_enabled);
    
        if(is_collision_detection_enabled)
        {
            real sum_of_the_radii = radius_of_star[i] + radius_of_star[j];
            if (distance_between_stars <= sum_of_the_radii)
            {
                int stopping_index  = next_index_for_stopping_condition();
                set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                set_stopping_condition_particle_index(
                    stopping_index,
                    0,
                    identifier_of_star[i]
                );
                set_stopping_condition_particle_index(
                    stopping_index,
                    1,
                    identifier_of_star[j]
                );
            }
        }

    and in FORTRAN:

    .. code-block:: fortran

        INCLUDE "../../lib/stopcond/stopcond.inc"
	INTEGER::is_stopping_conditions_enabled
        INTEGER::is_collision_detection_enabled, stopping_index
        INTEGER::error
        
        error = is_stopping_condition_enabled(COLLISION_DETECTION, is_collision_detection_enabled)
        !...
        IF is_collision_detection_enabled.EQ.1 THEN
          IF there_is_a_collision THEN
            stopping_index = next_index_for_stopping_conditions()
            error = set_stopping_condition_info(stopping_index, COLLISION_DETECTION)
        !...
          ENDIF
        ENDIF          

4.  Break the outer loop when any condition is set. The
    :c:data:`set_conditions` can be tested:

    .. code-block:: c

        while(...) {
            ...
            if(set_conditions)
            {
                break;
            }
        }

    .. code-block:: fortran

        INCLUDE "../../lib/stopcond/stopcond.inc"
        INTEGER::is_any_condition_set
        DO WHILE ()
        !...   
          IF (is_any_condition_set().EQ.1) EXIT
          
        END DO 


A good example of a code implementing the stopping conditions
is the "hermite" community code. The code can be found
in the 'src/amuse/community/hermite' directory.
.. _reference-label:

Reference documentation
=======================

.. toctree::
   :maxdepth: 2
   
   quantities_and_units
   particles
   particle_attributes
   bridge
   available-codes
   core_support
   legacy_support
   fileformat
   report
   stopping_conditions
   from_codes_to_data
   ioframework
   options
   directory-structure
   incode_storage
   cuda-setup
   simplified_function_interface
   interface_specification
   message-protocol
   distributed
   cartesius
   slurm
   style_guide
   code_management
   glossary
    
===================
Particle Attributes
===================

Particle attributes are defined in the 
:py:mod:`amuse.datamodel.particle_attributes`. These attributes
can be accessed in two ways, on the particle(s) or as a function
imported from the module. When accessed on the particles, the
first parameter (usually called ```particles```) must not
be given.

To access these functions directly do:

.. code-block:: python

    from amuse.datamodel.particle_attributes import *
    from amuse.lab import new_plummer_sphere
    
    particles = new_plummer_sphere(100)
    print kinetic_enery(particles)

To access these functions as an attribute do:

.. code-block:: python

    from amuse.lab import new_plummer_sphere
    
    particles = new_plummer_sphere(100)
    
    print particles.kinetic_energy()
    

.. automodule:: amuse.datamodel.particle_attributes
    :members:
=============================================
Running AMUSE through slurm
=============================================

Many supercomputers use slurm (see
https://slurm.schedmd.com/documentation.html) as a job-submission
environment.  The Cartesius supercomputer, for example, uses slurm,
but also the Leiden university supercomputer ALICE
https://wiki.alice.universiteitleiden.nl/index.php?title=ALICE_User_Documentation_Wiki&_ga=2.254956621.431957536.1583247665-1261213469.1568893554

The parallel AMUSE script
-------------------------

Slurm operates via batch scripts which have to be written and
subsequently submitted via slurm to the job scheduler.
We use a simple amuse example code as an example.
.. code-block:: sh

	> cd ${AMUSE_DIR}/examples/tutorial/parallel_amuse_script.py

The slurm batch script
----------------------
	
The slurm script to run this code with 6 processors can be callded run_example.sh

.. code-block:: sh

        #!/bin/bash
	sbatch <<EOT
	#!/bin/sh
	#SBATCH --mem=1000
	#SBATCH -p cpu-short
	#SBATCH -N 1
	#SBATCH -n 6

	export OMP_NUM_THREADS=6
	export OMPI_MCA_rmaps_base_oversubscribe=yes 
	export OMPI_MCA_mpi_warn_on_fork=0
	export OMPI_MCA_rmaps_base_oversubscribe=yes

	module load AMUSE/12.0.0-foss-2018a-Python-2.7.14
	module load openmpi/gcc/64/1.10.7
	module load matplotlib/3.1.1-foss-2019b-Python-3.7.4

	mpiexec -n 6 python -u parallel_amuse_script.py
	EOT

The various ingredients of this script identifies the amount of memory
requested, the queue (here called cpu-short), the number of nodes and
processors.

In the next block the environment variable are set.

Then the various modules are loaded.

And finally, the mpiexec command starts the python scrip on 6
processors and runs the script called parallel_amuse_script.py


Running the script
------------------

A slurm script is generally started with some command line, for example:	

.. code-block:: sh

	> sbatch run_example.sh

The output can subsequently be downloaded via scp	

.. code-block:: sh

	> scp user@remote_computer_address:file .
======================
Reporting during a run
======================

.. autoclass:: amuse.io.ReportTable
    :members:
=======================================
Stellar Evolution Interface Definition
=======================================

Introduction
~~~~~~~~~~~~
In this chapter we describe the common interface for stellar evolution codes.
Currently the interface for stellar evolutions codes that store state of
a star is specified. The Stellar Evolution codes that get all the needed
state from the function interface (are stateless), are not yet described.     

Parameters
~~~~~~~~~~
Stellar Evolution codes have at least one specified parameter. Other parameters need
to be specified on a per code basis. All parameters have to be accessed with functions following
the template of the ``get_metallicity`` and ``set_metallicity`` functions. A parameter access function may only
retrieve or update the value of a single parameter. After all parameters have been set, the 
``initialize_code`` function should be called, this gives the code the opportunity prepare the
model.

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: get_metallicity, set_metallicity, initialize_code
   

Object Management
~~~~~~~~~~~~~~~~~
A number of stellar evolution codes work on star objects. The following 
methods define the functionality to create, remove and query the particles in the code. 
*Currently the interface does not specify query function for stellar evolution, see stellar evolution for possible direction*

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: new_particle, delete_star

Object State
~~~~~~~~~~~~~~
To support properties (like acceleration), the code must define ``get_`` and ``set_`` functions. These
functions must get or set one scalar property (1 argument) or a vector property (3 arguments)
*Currently only get functions are specified*

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: get_mass, get_radius, get_luminosity, get_temperature, get_age, get_stellar_type, get_stellar_type


Model evolution
~~~~~~~~~~~~~~~
The stellar evolution codes evolve the properties of the star in time. The following functions
are needed to control the evolution in the code.

.. autoclass:: amuse.community.interface.se.StellarEvolution
   :members: commit_particles, evolve


Diagnostics
~~~~~~~~~~~
The state of the code can be queried, before, during and after the model calculations. 
*Currently no specific stellar evolution diagnostics functions have been defined*


Services
~~~~~~~~
Some stellar evolution codes can provide services for other codes. 
*Currently no specific stellar evolution service functions have been defined*
.. _supported-codes-label:

===================================
Currently supported Community Codes
===================================

Introduction
~~~~~~~~~~~~

.. note::
  This document is not up-to-date!

Here we provide an overview of some currently supported community codes in AMUSE,
along with a concise explanation of how each code works.
This document serves as an initial guide in finding the code with the highest
applicability to a specific astrophysical problem. The supported codes
have been sorted according to their astrophysical domain:

* :ref:`Dynamics`
* :ref:`Evolution`
* :ref:`Hydrodynamics`
* :ref:`Radiative-Transfer`

.. _Dynamics:

Stellar Dynamics
~~~~~~~~~~~~~~~~

* BHtree_
* hermite_
* HiGPUs_
* huayno_
* mercury_
* mmc_
* octgrav_
* phiGRAPE_
* smallN_
* twobody_

General
-------

The general parameters and methods for the gravitational dynamics are described in
:doc:`../reference/stellar_dynamics_interface_specification`.
Here we describe the exceptions for the specific codes under "Specifics".

Code comparison table
---------------------

Note: holds for AMUSE implementation.

+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+ 
|name     |approximation scheme |timestep scheme  |CPU     |GPU      |GRAPE |language |stopcond (1)  |parallel (2) |
|         |                     |                 |        |         |      |         |              |             |
+=========+=====================+=================+========+=========+======+=========+==============+=============+
|bhtree   |tree                 |shared/fixed     |Y       |N        |N     |C/C++    |CST           |N            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|bonsai   |tree                 |                 |        |         |      |         |              |             |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|fi       |tree                 |block/variable   |Y       |N        |N     |FORTRAN  |S             |N            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|hermite  |direct               |shared/variable  |Y       |N        |N     |C/C++    |CSOPT         |Y            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|HiGPUs   |direct               |block time steps |N       |Y        |N     |C/C++    |              |Y (on gpus   |
|         |                     |                 |        |         |      |         |              |cluster)     |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|huayno   |Approx symplectic    |                 |Y       |Y(opencl)|N     |C        |              |N            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|gadget   |tree                 |individual       |Y       |N        |N     |C/C++    |S             |Y            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|mercury  |MVS symplectic       |                 |Y       |N        |N     |FORTRAN  |              |N            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|octgrav  |tree                 |shared           |N       |Y        |N     |C/C++    |S             |N            | 
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|rebound  |                     |                 |        |         |      |         |              |N            |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|phigrape |direct               |block/variable   |Y g6    |Y        |Y     |FORTAN   |CSPT          |N            | 
|         |                     |                 |        |sapporo  |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|smallN   |direct Hermite 4th   |individual       |Y       |N        |N     |C/C++    |              |N            |
|         |order                |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+
|twobody  |universal variables, |none, exact      |Y       |N        |N     |Python   |              |N            |
|         |Kepler eq.           |                 |        |         |      |         |              |             |
|         |                     |                 |        |         |      |         |              |             |
+---------+---------------------+-----------------+--------+---------+------+---------+--------------+-------------+



(1) stopping conditions
   
    ==== ==========================
    code name of stopping condition
    ==== ========================== 
    C    Collision detection
    E    Escaper detection
    S    Number of steps detection
    O    Out of box detection
    P    Pair detection
    T    Timeout detection
    ==== ==========================

(2) Parallel in the following sense: AMUSE uses MPI to communicate with the codes, but
    for some codes it can be used to parallelize the calculations. Some codes (GPU) 
    are already parallel, however in this table we *do not* refer to that.
  
    Codes designated *Y* for parallel can set the number of (parallel) workers, e.g. to set
    10 workers for hermite do:

    .. code-block:: python

       >>> instance = Hermite(number_of_workers=10)

.. _BHtree:

BHtree
------

N-body integration module. An implementation of the Barnes & Hut tree code [#bh]_
by Jun Makino  `BHTree-code <http://jun.artcompsci.org//softwares/C++tree/index.html>`_.

Specifics
#########

Parameters
^^^^^^^^^^

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  |unit           |description        |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|use_self_gravity                    |1              |none           |flag for usage of  |
|                                    |               |               |self gravity, 1 or |
|                                    |               |               |0 (true or false)  |
+------------------------------------+---------------+---------------+-------------------+
|timestep                            |0.015625       |time           |time step          |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.125          |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |>0!                |
+------------------------------------+---------------+---------------+-------------------+
|ncrit_for_tree                      |1024           |none           |maximum number of p|
|                                    |               |               |articles sharing an|
|                                    |               |               |interaction list   |
+------------------------------------+---------------+---------------+-------------------+
|opening_angle                       |0.75           |none           |opening angle,     |
|                                    |               |               |theta, for building|
|                                    |               |               |the tree: between 0|
|                                    |               |               |and 1              |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_number_of_steps |1              |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_timeout         |4.0            |seconds        |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_out_of_box_size |0.0            |length         |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|time                                |0.0            |time           |current            |
|                                    |               |               |simulation         |
|                                    |               |               |time               |
+------------------------------------+---------------+---------------+-------------------+
|dt_dia                              | 1.0           |time           |time interval      |
|                                    |               |               |between            |
|                                    |               |               |diagnostics        |
|                                    |               |               |output             |
+------------------------------------+---------------+---------------+-------------------+

.. automodule:: amuse.community.bhtree.interface
    
    .. autoclass:: BHTree
    
        .. autoparametersattribute:: parameters
        
example

.. code-block:: python

     >>> from amuse.community.bhtree.interface import BHTreeInterface, BHTree
     >>> from amuse.units import nbody_system
     >>> instance = BHTree(BHTree.NBODY)
     >>> instance.parameters.epsilon_squared = 0.00001 | nbody_system.length**2

.. [#bh] Barnes, J. & Hut, P. 1986. *Nature* **324**, 446.

.. _hermite:

hermite
--------

Time-symmetric N-body integration module with shared but variable time step
(the same for all particles but its size changing in time),
using the Hermite integration scheme [#hutmakino]_. See also : `ACS <http://www.artcompsci.org/>`_

Specifics
#########

Parameters
^^^^^^^^^^

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|pair factor                         |1.0            | none          |radius factor      |
|                                    |               |               |for pair detec     |
|                                    |               |               |tion               |
+------------------------------------+---------------+---------------+-------------------+
|dt_param                            |0.03           |none           |timestep           |
|                                    |               |               |scaling factor     |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.0            |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_number_of_steps |1              |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_timeout         |4.0            |seconds        |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_out_of_box_size |0.0            |length         |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|time                                |0.0            |time           |current            |
|                                    |               |               |simulation         |
|                                    |               |               |time               |
+------------------------------------+---------------+---------------+-------------------+
|dt_dia                              | 1.0           |time           |time interval      |
|                                    |               |               |between            |
|                                    |               |               |diagnostics        |
|                                    |               |               |output             |
+------------------------------------+---------------+---------------+-------------------+

.. automodule:: amuse.community.hermite.interface
    
    .. autoclass:: Hermite
    
        .. autoparametersattribute:: parameters

.. [#hutmakino] Hut, P., Makino, J. & McMillan, S., 1995, *ApJL* **443**, L93.

.. _phiGRAPE:

phiGRAPE
--------

phiGRAPE is a direct N-body code optimized for running on a parallel
GRAPE cluster. See Harfst et al. [#harfst]_ for more details.  The
Amusean version is capable of working on other platforms as well by
using interfaces that mimic GRAPE hardware.

* Sapporo
    Sapporo is a library that mimics the behaviour of
    GRAPE hardware and uses the GPU to execute the force calculations [#gaburov]_.

* Sapporo-light
    This version of Sapporo is without
    multi-threading support and does not need
    C++. This makes it easier to integrate into
    fortran codes, but beware, it can only use one
    GPU device per application!

* g6
    Library which mimics the behavior of GRAPE and uses the
    CPU. Lowest on hw requirements.

Specifics
#########

Hardware modes
^^^^^^^^^^^^^^

Parameters
^^^^^^^^^^

Amuse tries to build all implementations at compile time. In the phiGRAPE interface
module the preferred mode can be selected whith the mode parameter:

* MODE_G6LIB = 'g6lib'
    Just make it work, no optimizations, no special hw requirements
* MODE_GPU   = 'gpu'
    Using sapporo, CUDA needed.
* MODE_GRAPE = 'grape'
    Using GRAPE hw.
* MODE_PG    = 'pg'
    Phantom grape, optimized for x86_64 processors

.. code-block:: python

   >>> from amuse.community.phigrape.interface import PhiGRAPEInterface, PhiGRAPE
   >>> instance = PhiGRAPE(PhiGRAPE.NBODY, PhiGRAPEInterface.MODE_GPU)

The default is **MODE_G6LIB**.

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|initialize_gpu_once                 |0              |none           |set to 1 if the gpu|
|                                    |               |               |must only be       |
|                                    |               |               |initialized once, 0|
|                                    |               |               |if it can be       |
|                                    |               |               |initialized for    |
|                                    |               |               |every call\nIf you |
|                                    |               |               |want to run        |
|                                    |               |               |multiple instances |
|                                    |               |               |of the code on the |
|                                    |               |               |same gpu this      |
|                                    |               |               |parameter needs to |
|                                    |               |               |be 0 (default)     |
+------------------------------------+---------------+---------------+-------------------+
|initial_timestep_parameter          |0.0            |none           |parameter to       |
|                                    |               |               |determine the      |
|                                    |               |               |initial timestep   |
+------------------------------------+---------------+---------------+-------------------+
|timestep_parameter                  |0.0            |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.0            |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |>0!                |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_number_of_steps |1              |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_timeout         |4.0            |seconds        |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_out_of_box_size |0.0            |length         |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+

.. automodule:: amuse.community.phigrape.interface
    
    .. autoclass:: PhiGRAPE
    
        .. autoparametersattribute:: parameters


.. code-block:: python

    >>> instance.timestep_parameter = 0.1 |nbody_system.time


.. [#harfst]  Harfst, S., Gualandris, A., Merritt, D., Spurzem, R., Portegies Zwart, S., & Berczik, P. 2006, *NewAstron.* **12**, 357-377.
.. [#gaburov]  Gaburov, E., Harfst, S., Portegies Zwart, S. 2009, *NewAstron.* **14** 630-637.

.. _twobody:

twobody
-------

Semi analytical code based on Kepler [#bate]_. The particle set provided has length one or two. If one particle is given, the mass is assigned to
a particle in the origin and the phase-coordinates are assigned to the other particle. This is usefull when *m1* >> *m2*.

Specifics
#########

Parameters
^^^^^^^^^^


.. [#bate]  Bate, R.R, Mueller, D.D., White, J.E. "FUNDAMENTALS OF ASTRODYNAMICS" *Dover* 0-486-60061-0

.. _smallN:

smallN
------

Interface to the Kira Small-N Integrator and Kepler modules from
Starlab. https://www.sns.ias.edu/~starlab/ 

You will need to download Starlab from the above site, make it, install
it, and then set the STARLAB_INSTALL_PATH variable to be equal to the
installation directory (typically something like ~/starlab/usr).

Starlab is available under the GNU General Public Licence (version 2),
and is developed by:

    * Piet Hut
    * Steve McMillan
    * Jun Makino
    * Simon Portegies Zwart

Other Starlab Contributors:

    * Douglas Heggie
    * Kimberly Engle
    * Peter Teuben

Specifics
#########

Parameters
^^^^^^^^^^
+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.0            |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |>0!                |
+------------------------------------+---------------+---------------+-------------------+
|number_of_particles                 |0.0            |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+


Methods
^^^^^^^

.. _octgrav:

Octgrav
-------

Tree-code which runs on GPUs with NVIDIA CUDA architecture. [#oct]_

Specifics
#########

Parameters
^^^^^^^^^^

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|opening_angle                       |0.8            |none           |opening angle for  |
|                                    |               |               |building the tree  |
|                                    |               |               |between 0 and 1    |
+------------------------------------+---------------+---------------+-------------------+
|timestep                            |0.01           |time           |constant timestep  |
|                                    |               |               |for iteration      |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.01           |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_number_of_steps |1              |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_timeout         |4.0            |seconds        |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_out_of_box_size |0.0            |length         |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+

.. automodule:: amuse.community.octgrav.interface
    
    .. autoclass:: Octgrav
    
        .. autoparametersattribute:: parameters


.. [#oct] Gaburov, E., Bedorf, J., Portegies Zwart S., 2010, "Gravitational tree-code on graphics processing units: implementations in CUDA", *ICCS*


Specifics
#########

Parameters
^^^^^^^^^^

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|opening_angle                       |0.8            |none           |opening angle for  |
|                                    |               |               |building the tree  |
|                                    |               |               |between 0 and 1    |
+------------------------------------+---------------+---------------+-------------------+
|timestep                            |0.01           |time           |constant timestep  |
|                                    |               |               |for iteration      |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|epsilon squared                     |0.01           |length*length  |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_number_of_steps |1              |none           |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_timeout         |4.0            |seconds        |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|stopping_conditions_out_of_box_size |0.0            |length         |                   |
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+

.. _mercury:

mercury 
-------

Mercury is a general-purpose N-body integration package for problems
in celestial mechanics.

This package contains some subroutines taken from the Swift
integration package by H.F.Levison and M.J.Duncan (1994) Icarus, vol
108, pp18.  Routines taken from Swift have names beginning with
`drift` or `orbel`.

The standard symplectic (MVS) algorithm is described in J.Widsom and
M.Holman (1991) Astronomical Journal, vol 102, pp1528.

The hybrid symplectic algorithm is described in J.E.Chambers (1999)
Monthly Notices of the RAS, vol 304, pp793.

Currently Mercury has an interface that differs from the other grav
dyn interfaces. The class is called MercuryWayWard and will loose its
predicate once the interface is standardized (work in progress). It
handles two kinds of particles: centre_particle and orbiters. The
centre_particle is restricted to contain only one particle and it
should be the heaviest and much heavier than the orbiters. Its at the
origin in phase space.

Appart from the usual phase space coordinate, particles have a spin in
mercury and the centre particle has oblateness parameters expressed in
moments 2, 4 and 6. Orbiters have density.

Furthermore, mercury does not use nbody units but instead units as listed in table.

Central particle
################

+--------------+-----------+------------------+
|mass          |           |MSun              |
+--------------+-----------+------------------+
|radius        |           |AU                |
+--------------+-----------+------------------+
|oblateness    |j2, j4, j6 |AU^2 etc.         |
|              |           |                  |
|              |           |                  |
+--------------+-----------+------------------+
|angluar       |           |MSun AU ^2 day ^-1|
|momentum      |           |                  |
+--------------+-----------+------------------+

Orbiters
########

+--------------+-----------+-----------+------------------+
|mass          |           |           |MSun              |
+--------------+-----------+-----------+------------------+
|density       |           |           |g/cm^3            |
+--------------+-----------+-----------+------------------+
|position      |           |x, y, z    |AU                |
|              |           |           |                  |
|              |           |           |                  |
+--------------+-----------+-----------+------------------+
|velocity      |           |vx, vy, vz |AU/day            |
|              |           |           |                  |
+--------------+-----------+-----------+------------------+
|celimit       |close      |           |Hill radii, but   |
|              |encounters |           |units.none in     |
|              |           |           |amuse             |
+--------------+-----------+-----------+------------------+
|angluar       |           |Lx, Ly, Lz |MSun AU ^2 day ^-1|
|momentum      |           |           |                  |
+--------------+-----------+-----------+------------------+

.. automodule:: amuse.community.mercury.interface
    
    .. autoclass:: MercuryWayWard
    
        .. autoparametersattribute:: parameters
        
.. _huayno:

Huayno
------

Hierarchically split-Up Approximately sYmplectic N-body sOlver 
(HUAYNO)

Inti Pelupessy - january 2011

short description
#################

  HUAYNO is a code to solve the astrophysical N-body problem. It uses
  recursive Hamiltonian splitting to generate multiple-timestep integrators
  which conserve momentum to machine precision. A number of different 
  integrators are available. The code has been developed within the 
  AMUSE environment. It can make use of GPUs - for this an OpenCL 
  version can be compiled.

Use
###

  Use of the code is the same as any gravity the code within AMUSE. There are
  three parameters for the code: Smoothing length (squared) eps2, timestep
  parameter eta and a parameter to select the integrator:

  timestep_parameter( eta ): recommended values eta=0.0001 - 0.05
  epsilon_squared (eps2) : eps2 can be zero or non-zero
  inttype_parameter( inttype ): possible values for inttype are described below. 

  Miscellaneous:

  - The code assumes G=1, 
  - Collisions are not implemented (needs rewrite),
  - workermp option can be used for OpenMP parallelization,
  - The floating point precision of the calculations can be changed by setting
    the FLOAT and DOUBLE definitions in evolve.h. FLOAT sets the precision of
    the calculations, DOUBLE the precision of the position and velocity
    reductions. They can be set to e.g. float, double, long double or __float128
    It is advantageous to choose set DOUBLE at a higher precision than FLOAT.
    recommended is the combination: double/ long double  
  - the AMUSE interface uses double precision
  - It is unlikely the integer types in evolve.h would need to be changed
    (INT should be able to hold particle number, LONG should be able to hold
    interaction counts) 

  OpenCL operation:

  By compiling worker_cl it is possible to offload the force and timestep
  loops to the GPU.  The implementation is based on the STDCL library
  (www.browndeertechnology.com) so this library should be compiled first. 
  In the Makefile the corresponding directories should point to the
  installation directory of STDCL. The define in evolce_cl.h should be set
  appropiately for the OpenCL configuration: CLCONTEXT:  stdgpu or stdcpu
  NTHREAD:  64 for GPU, 2 for CPU (actually 64 will also work for CPU just
  fine) BLOCKSIZE: number of particles stored in local memory (64 good
  starting guess) evolve_kern.cl contains the OpenCL kernels. Precision of
  the calculation is controlled by FLOAT(4) defines in evolve_kern.cl and
  CLFLOAT(4) in evolve.h. They should agree with each other (i.e. float and
  cl_float or  double and cl_double)

.. _mmc:
   
Monte Carlo
-----------


Giersz, M. 1998, MNRAS, 298, 1239
Giersz, M. 2001, MNRAS, 324, 218
Giersz, M. 2006, MNRAS, 371, 484


Specifics
#########

parameters
^^^^^^^^^^

+--------+------------------------------------------+----+
|name    | description                              |    |
+========+==========================================+====+
|irun    |initial sequence of random numbers        |    |
|        |                                          |    |
|        |                                          |    |
|        |                                          |    |
+--------+------------------------------------------+----+
|nt      |total number o f objects (sta rs and      |    |
|        |binarie s) at T=0 ns - number of single   |    |
|        |tars, nb - number o f binaries, nt ns+nb, |    |
|        |nss - number of star s(nss = nt+nb)       |    |
+--------+------------------------------------------+----+
|istart  |1 - initial model, .ne.1 - restart        |    |
+--------+------------------------------------------+----+
|ncor    |number of stars to calculate the central  |    |
|        |parameters                                |    |
+--------+------------------------------------------+----+
|nmin    |minimum number of stars to calculate the  |    |
|        |central parameters                        |    |
+--------+------------------------------------------+----+
|nz0     |number of stars in each zone at T=0       |    |
+--------+------------------------------------------+----+
|nzonc   |minimum number of zones in the core       |    |
+--------+------------------------------------------+----+
|nminzo  |minimum number of stars in a zone         |    |
+--------+------------------------------------------+----+
|ntwo    |maximum index of 2                        |    |
+--------+------------------------------------------+----+
|imodel  |initial model: 1- uniform & isotropic, 2- |    |
|        |Plummer, 3- King, 4 - M67                 |    |
+--------+------------------------------------------+----+
|iprint  |0- full diagnostic information, 1-        |    |
|        |diagnostic info.  suppressed              |    |
+--------+------------------------------------------+----+
|ib3f    |1 - Spitzer's, 2 - Heggie's formula for   |    |
|        |three-body binary interaction with field  |    |
|        |stars, 3 - use Pmax for interaction *     |    |
|        |probability 4 - three- and four-body      |    |
|        |numerical integration                     |    |
|        |                                          |    |
+--------+------------------------------------------+----+
|iexch   |0 - no exchange in any interactions, 1 -  |    |
|        |exchange only in binary field star        |    |
|        |interacions, 2 - exchange in all          |    |
|        |interactions (binary - field and binary - |    |
|        |binary)                                   |    |
+--------+------------------------------------------+----+
|tcrit   |termination time in units of the crossing |    |
|        |time                                      |    |
+--------+------------------------------------------+----+
|tcomp   |maximum computing time in hours           |    |
+--------+------------------------------------------+----+
|qe      |energy tolerance                          |    |
+--------+------------------------------------------+----+
|alphal  |power-law index for initial mass function |    |
|        |for masses smaller than breake mass: -1 - |    |
|        |equal mass case                           |    |
+--------+------------------------------------------+----+
|alphah  |power-law index for initial mass function |    |
|        |for masses greater than breake mass. If   |    |
|        |alphal=alphah the IMF does not have a     |    |
|        |break                                     |    |
+--------+------------------------------------------+----+
|brakem  |the mass in which the IMF is broken. If   |    |
|        |brakem is smaller * than the minimum mass |    |
|        |(bodyn) than the break mass is as for *   |    |
|        |the Kroupa mass function (brakem = 0.5 Mo)|    |
|        |                                          |    |
+--------+------------------------------------------+----+
|body1   |maximum particle mass before scaling      |    |
|        |(solar mass)                              |    |
+--------+------------------------------------------+----+
|bodyn   |minimum particle mass before scaling      |    |
|        |(solar mass)                              |    |
+--------+------------------------------------------+----+
|fracb   |primordial binary fraction by number. nb =|    |
|        |fracb*nt, * ns = (1 - fracb)*nt, nss = (1 |    |
|        |+ fracb)*nt * fracb > 0 - primordial      |    |
|        |binaries * fracb = 0 - only dynamical     |    |
|        |binaries                                  |    |
+--------+------------------------------------------+----+
|amin    |minimum semi-major axis of binaries (in   |    |
|        |sollar units) * = 0 then amin = 2*(R1+R2),|    |
|        |> 0 then amin = amin                      |    |
+--------+------------------------------------------+----+
|amax    |maximum semi-major axis of binaries (in   |    |
|        |sollar units)                             |    |
+--------+------------------------------------------+----+
|qvir    |virial ratio (qvir = 0.5 for equilibrium) |    |
+--------+------------------------------------------+----+
|rbar    |tidal radius in pc, halfmass radius in pc |    |
|        |for isolated * cluster. No scaling - rbar |    |
|        |= 1                                       |    |
+--------+------------------------------------------+----+
|zmbar   |total mass of the cluster in sollar mass, |    |
|        |* no scaling zmbar = 1                    |    |
+--------+------------------------------------------+----+
|w0      |king model parameter                      |    |
+--------+------------------------------------------+----+
|bmin    |minimum value of sin(beta^2/2)            |    |
+--------+------------------------------------------+----+
|bmax    |maximum value of sin(beta^2/2)            |    |
+--------+------------------------------------------+----+
|tau0    |time step for a complite cluster model    |    |
+--------+------------------------------------------+----+
|gamma   |parameter in the Coulomb logarithm        |    |
|        |(standard value = 0.11)                   |    |
+--------+------------------------------------------+----+
|xtid    |coeficient in the front of cluster tidal  |    |
|        |energy: * -xtid*smt/rtid                  |    |
+--------+------------------------------------------+----+
|rplum   |for M67 rtid = rplum*rsplum (rsplum -     |    |
|        |scale radius for * plummer model)         |    |
+--------+------------------------------------------+----+
|dttp    |time step (Myr) for profile output        |    |
+--------+------------------------------------------+----+
|dtte    |time step (Myr) for mloss call for all    |    |
|        |objects                                   |    |
+--------+------------------------------------------+----+
|dtte0   |time step (Myr) for mloss call for all    |    |
|        |objects for tphys * less then tcrevo. For |    |
|        |tphys greater then tcrevo time step * is  |    |
|        |eqiual to dtte                            |    |
+--------+------------------------------------------+----+
|tcrevo  |critical time for which time step for     |    |
|        |mloss call changes from * dtte0 to dtte   |    |
+--------+------------------------------------------+----+
|xtau    |call mloss for a particlular object when *|    |
|        |(uptime(im1) - olduptime(im1))/tau/tscale |    |
|        |< xtau                                    |    |
+--------+------------------------------------------+----+
|ytau    |multiplication of tau0 (tau = ytau*tau0)  |    |
|        |after time * greater than tcrevo          |    |
+--------+------------------------------------------+----+
|ybmin   |multiplication of bmin0 (bmin =           |    |
|        |ybmin*bmin0) after time * greater than    |    |
|        |tcrevo                                    |    |
+--------+------------------------------------------+----+
|zini    |initial metalicity (solar z = 0.02,       |    |
|        |globular clusters * M4 - z = 0.002,       |    |
|        |NGC6397 - z = 0.0002)                     |    |
+--------+------------------------------------------+----+
|ikroupa |0 - the initial binary parameters are     |    |
|        |picked up * according Kroupa's            |    |
|        |eigenevolution and feeding algorithm *    |    |
|        |(Kroupa 1995, MNRAS 277, 1507) * 1 - the  |    |
|        |initial binary parameters are picked as   |    |
|        |for M67 * model (Hurley et al. 2005)      |    |
+--------+------------------------------------------+----+
|iflagns |0 - no SN natal kiks for NS formation, 1 -|    |
|        |SN natal kicks only for single NS         |    |
|        |formation, 2 - SN natal kick for single NS|    |
|        |formation and NS formation in binaries    |    |
+--------+------------------------------------------+----+
|iflagbh |0 - no SN natal kiks for BH formation, 1 -|    |
|        |SN natal kicks * only for single BH       |    |
|        |formation, 2 - SN natal kick for single * |    |
|        |BH formation and BH formation in binaries |    |
+--------+------------------------------------------+----+
|nitesc  |0 - no iteration of the tidal radius and  |    |
|        |induced mass loss * due to stellar        |    |
|        |evolution, 1 - iteration of the tidal     |    |
|        |radius * and induced mass loss due to     |    |
|        |stellar evolution                         |    |
+--------+------------------------------------------+----+


.. automodule:: amuse.community.mmc.interface
    
    .. autoclass:: mmc
    

HiGPUs
--------

HiGPUs is a parallel direct N-body code based on a 6th order Hermite 
integrator. The code has been developed by Capuzzo-Dolcetta, Punzo 
and Spera (Dep. of Physics, Sapienza, Univ. di Roma; 
see astrowww.phys.uniroma1.it/dolcetta/HPCcodes/HiGPUs.html) and uses, at 
the same time, MPI, OpenMP and CUDA libraries to fully exploit all 
the capabilities offered by hybrid supercomputing platforms. 
Moreover, it is implemented using block time steps (individual time stepping) 
such to be able to deal with stiff problems like highly collisional gravitational 
N-body problems.

Specifics
#########

Parameters
^^^^^^^^^^

+------------------------------------+---------------+---------------+-------------------+
|name                                |default value  | unit          | description       |
|                                    |               |               |                   |
|                                    |               |               |                   |
+====================================+===============+===============+===================+
|eta_6                               |0.4            |none           |eta parameter for  |
|                                    |               |               |determining stars  |
|                                    |               |               |time steps         |
+------------------------------------+---------------+---------------+-------------------+
|eta_4                               |0.01           |none           |eta parameter for  |
|                                    |               |               |initializing blocks|
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|eps                                 |0.001          |length         |smoothing parameter|
|                                    |               |               |for gravity        |
|                                    |               |               |calculations       |
+------------------------------------+---------------+---------------+-------------------+
|r_scale_galaxy                      |0.0            |length         |scale radius for   |
|                                    |               |               |analytical galaxy  |
|                                    |               |               |potential          |
+------------------------------------+---------------+---------------+-------------------+
|mass_galaxy                         |0.0            |mass           |total mass for     |
|                                    |               |               |analytical galaxy  |
|                                    |               |               |potential          |

+------------------------------------+---------------+---------------+-------------------+
|r_core_plummer                      |0.0            |length         |core radius for    |
|                                    |               |               |analytical plummer |
|                                    |               |               |potential          |
+------------------------------------+---------------+---------------+-------------------+
|mass_plummer                        |0.0            |mass           |total mass for     |
|                                    |               |               |analytical plummer |
|                                    |               |               |potential          |
+------------------------------------+---------------+---------------+-------------------+
|start_time                          |0.0            |time           |initial simulation |
|                                    |               |               |time               |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|min_step                            |-30.0          |none           |exponent which defi|
|                                    |               |               |nes the minimum    |
|                                    |               |               |time step allowed  |
|                                    |               |               |for stars          |
|                                    |               |               |(2^exponent)       |
+------------------------------------+---------------+---------------+-------------------+
|max_step                            |-3.0           |none           |exponent which defi|
|                                    |               |               |nes the maximum    |
|                                    |               |               |time step allowed  |
|                                    |               |               |for stars          | 
|                                    |               |               |(2^exponent)       |
+------------------------------------+---------------+---------------+-------------------+
|n_Print                             |1000000        |none           |maximum number of  | 
|                                    |               |               |snapshots          |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|dt_Print                            |1.0            |time           |time interval      |
|                                    |               |               |between            |
|                                    |               |               |diagnostics        |
|                                    |               |               |output             |
+------------------------------------+---------------+---------------+-------------------+
|n_gpu                               |2              |none           |number of GPUs per | 
|                                    |               |               |node               |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|gpu_name                            |GeForce GTX 480|none           |GPUs to use        | 
|                                    |               |               |                   |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|Threads                             |128            |none           |number of gpus     | 
|                                    |               |               |threads per block  |
|                                    |               |               |                   |
+------------------------------------+---------------+---------------+-------------------+
|output_path_name                    |../../test_re  |none           |path where HiGPUs  |
|                                    |sults/         |               |output will be     |
|                                    |               |               |stored             |
+------------------------------------+---------------+---------------+-------------------+



For more information about parameters check the readme file in the docs folder.
These are the maximum performance (in Gflops) reached using different single GPUs installed, one at a time, 
on a workstation equipped with 2 CPUs Intel Xeon X5650, 12 GB of ECC RAM memory 1333 MHz, 
Ubuntu Lucid 10.04 x86_64, motherboard Supermicro X8DTG-QF:
  
    * TESLA C1060  :  107
    * TESLA C2050  :  395
    * TESLA M2070  :  391
    * GeForce GTX 480  :  265
    * GeForce GTX 580  :  311

To use the code with AMP GPUs you can download the OpenCL version from the website

.. automodule:: amuse.community.higpus.interface
    
    .. autoclass:: HiGPUs
    
        .. autoparametersattribute:: parameters       



.. _Evolution:

Stellar Evolution
~~~~~~~~~~~~~~~~~

* bse_
* evtwin_
* mesa_
* seba_
* sse_


.. _sse:

sse
---

Stellar evolution is performed by the **rapid** single-star evolution (SSE) algorithm. This is a package of **analytical formulae** fitted to the detailed models of Pols et al. (1998)  that covers all phases of evolution from the zero-age main-sequence up to and including remnant phases. It is valid for masses in the range 0.1-100 Msun and metallicity can be varied. The SSE package contains a prescription for mass loss by stellar winds. It also follows the evolution of rotational angular momentum for the star. Full details can be found in the SSE paper:

* "Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity"
    Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543


=========== ====== ==== =========================
..           min    max  unit
=========== ====== ==== =========================
Mass        0.1    100  Msun
Metallicity 0.0001 0.03 fraction (0.02 is solar)
=========== ====== ==== =========================

.. automodule:: amuse.community.sse.interface
    
    .. autoclass:: SSE
    
        .. autoparametersattribute:: parameters
        
.. _bse:

bse
---

Binary evolution is performed by the **rapid** binary-star evolution (BSE) algorithm. Circularization of eccentric orbits and synchronization of stellar rotation with the orbital motion owing to tidal interaction is modelled in detail. Angular momentum loss mechanisms, such as gravitational radiation and magnetic braking, are also modelled. Wind accretion, where the secondary may accrete some of the material lost from the primary in a wind, is allowed with the necessary adjustments made to the orbital parameters in the event of any mass variations. Mass transfer also occurs if either star fills its Roche lobe and may proceed on a nuclear, thermal or dynamical time-scale. In the latter regime, the radius of the primary increases in response to mass-loss at a faster rate than the Roche-lobe of the star. Stars with deep surface convection zones and degenerate stars are unstable to such dynamical time-scale mass loss unless the mass ratio of the system is less than some critical value. The outcome is a common-envelope event if the primary is a giant star. This results in merging or formation of a close binary, or a direct merging if the primary is a white dwarf or low-mass main-sequence star. On the other hand, mass transfer on a nuclear or thermal time-scale is assumed to be a steady process. Prescriptions to determine the type and rate of mass transfer, the response of the secondary to accretion and the outcome of any merger events are in place in BSE and the details can be found in the BSE paper:

* "Evolution of binary stars and the effect of tides on binary populations"
    Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
* "Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity"
    Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543


============ ====== ==== =========================
..           min    max  unit
============ ====== ==== =========================
Mass         0.1    100  Msun
Metallicity  0.0001 0.03 fraction (0.02 is solar)
Period       all    all
Eccentricity 0.0    1.0
============ ====== ==== =========================

.. automodule:: amuse.community.bse.interface
    
    .. autoclass:: BSE
    
        .. autoparametersattribute:: parameters

.. _seba:

seba
----

Single-star and binary evolution is performed by the **rapid** stellar evolution algorithm SeBa. Stars are evolved from
the zero-age main sequence until remnant formation and beyond.
Single stellar evolution is modeled with **analytical formulae** based on fits to detailed single star tracks at different metallicities (Hurley, Pols & Tout, 2000, 315, 543). Stars are parametrised by mass, radius, luminosity, core mass, etc. as functions of time and initial mass. Mass
loss from winds, which is substantial e.g. for massive stars and
post main-sequence stars, is included.

Furthermore, the SeBa package contains an algorithm for **rapid** binary evolution calculations. 
Binary interactions such as wind accretion, tidal interaction and angular momentum loss through (wind) mass loss, magnetic braking, or gravitational radiation are taken into account at every timestep with appropriate recipes. 
The stability and rate of mass transfer are dependent on the reaction to mass
change of the stellar radii and the corresponding Roche lobes. If the mass transfer takes place on the dynamical timescale of the donor star, the mass transfer becomes quickly unstable, and a common-envelope phase follows. If mass transfer occurs in a stable way, SeBa models the response of the companion star and possible mass loss and angular momentum loss at every timestep.  
After mass transfer ceases in a binary system, the donor star turns into a remnant or a helium-burning
star without a hydrogen envelope. When instead, the mass transfer leads to a merger between the binary stars, the resulting stellar product is estimated and the subsequent evolution is followed. 
  
More information on SeBa can be found in the papers:    

Relevant papers:

* "Population synthesis of high-mass binaries"
  Portegies Zwart, S.F., Verbunt, F. 1996, 309, 179P
* "Population synthesis for double white dwarfs . I. Close detached systems"
  Nelemans, G., Yungelson, L.R., Portegies Zwart, S.F., Verbunt, F. 2001, 365, 491N
* "Supernova Type Ia progenitors from merging double white dwarfs. Using a new population synthesis model"
  Toonen, S., Nelemans, G., Portegies Zwart, S. 2012, 546A, 70T  
* "The effect of common-envelope evolution on the visible population of post-common-envelope binaries"
  Toonen, S., Nelemans, G. 2013, 557A, 87T

============ ====== ==== =========================
..           min    max  unit
============ ====== ==== =========================
Mass         0.1    100  Msun
Metallicity  0.0001 0.03 fraction (0.02 is solar)
Period       all    all
Eccentricity 0.0    1.0
============ ====== ==== =========================

.. automodule:: amuse.community.seba.interface
    
    .. autoclass:: SeBa
    
        .. autoparametersattribute:: parameters

.. _evtwin:

evtwin
------

Evtwin is based on Peter Eggleton's stellar evolution code, and actually solves the differential equations that apply to the interior of a star. Therefore it is more accurate, but also much slower than the analytic fits-based sse_ and seba_ algorithm explained above.
Binaries are not yet supported in the AMUSE interface to evtwin, neither is the work-around for the helium flash. Currently only solar metallicity.

Relevant papers:

* "The evolution of low mass stars"
   Eggleton, P.P. 1971, MNRAS, 151, 351
* "Composition changes during stellar evolution"
   Eggleton, P.P. 1972, MNRAS, 156, 361
* "A numerical treatment of double shell source stars"
   Eggleton, P.P. 1973, MNRAS, 163, 279
* "An Approximate Equation of State for Stellar Material"
   Eggleton, P.P., Faulkner, J., & Flannery, B.P. 1973, A&A, 23, 325
* "A Possible Criterion for Envelope Ejection in Asymptotic Giant Branch or First Giant Branch Stars"
   Han, Z., Podsiadlowski, P., & Eggleton, P.P. 1994, MNRAS, 270, 121
* "Approximate input physics for stellar modelling"
   Pols, O.R., Tout, C.A., Eggleton, P.P., & Han, Z. 1995, MNRAS, 274, 964
* "The Braking of Wind"
   Eggleton, P.P. 2001, Evolution of Binary and Multiple Star Systems, 229, 157
* "A Complete Survey of Case A Binary Evolution with Comparison to Observed Algol-type Systems"
   Nelson, C.A., & Eggleton, P.P. 2001, ApJ, 552, 664
* "The Evolution of Cool Algols"
   Eggleton, P.P., & Kiseleva-Eggleton, L. 2002, ApJ, 575, 461
* For thermohaline mixing:
   Stancliffe, Glebbeek, Izzard & Pols, 2007 A&A
* For the OPAL 1996 opacity tables:
   Eldridge & Tout, 2004 MNRAS 348
* For enhancements to the solver:
   Glebbeek, Pols & Hurley, 2008 A&A


.. automodule:: amuse.community.evtwin.interface
    
    .. autoclass:: EVtwin
    
        .. autoparametersattribute:: parameters


.. _mesa:

mesa
----

The software project MESA (Modules for Experiments in Stellar 
Astrophysics, `<http://mesa.sourceforge.net/>`_), aims to provide 
state-of-the-art, robust, and efficient open source modules, usable 
singly or in combination for a wide range of applications in stellar 
astrophysics. Since the package is rather big (about 800 MB 
download, >2 GB built), this community code is optional and does not 
install automatically. Set the environment variable DO_INSTALL_MESA 
and run `make` to download and install it. The AMUSE interface to 
MESA can create and evolve stars using the MESA/STAR module. If you 
order a metallicity you haven't used before, starting models will be 
computed automatically and saved in the 
`mesa/src/data/star_data/starting_models` directory (please be 
patient...). All metallicities are supported, even the interesting 
case of Z=0. The supported stellar mass range is from about 0.1 to 
100 Msun.

References:

* Paxton, Bildsten, Dotter, Herwig, Lesaffre & Timmes 2010, ApJS submitted, arXiv:1009.1622
* `<http://mesa.sourceforge.net/>`_


.. automodule:: amuse.community.mesa.interface
    
    .. autoclass:: MESA
    
        .. autoparametersattribute:: parameters



.. _Hydrodynamics:

Hydrodynamics
~~~~~~~~~~~~~

* athena_ (grid code)
* capreole_ (grid code)
* fi_ (N-body/SPH code)
* gadget2_ (N-body/SPH code)


.. _athena:

athena
--------

Athena is a grid-based code for astrophysical hydrodynamics. Athena can solve magnetohydrodynamics (MHD) as well, but this is currently not supported from AMUSE. It was developed primarily for studies of the interstellar medium, star formation, and accretion flows.

The current version (Athena v4.0) implements algorithms for the following physics:
  * compressible hydrodynamics and MHD in 1D, 2D, and 3D,
  * ideal gas equation of state with arbitrary γ (including γ = 1, an isothermal EOS),
  * an arbitrary number of passive scalars advected with the flow,
  * self-gravity, and/or a static gravitational potential,
  * Ohmic resistivity, ambipolar diffusion, and the Hall effect,
  * both Navier-Stokes and anisotropic (Braginskii) viscosity,
  * both isotropic and anisotropic thermal conduction,
  * optically-thin radiative cooling. 

In addition, Athena allows for the following grid and parallelization options:
  * Cartesian or cylindrical coordinates,
  * static (fixed) mesh refinement,
  * shearing-box source terms, and an orbital advection algorithm for MHD,
  * parallelization using domain decomposition and  MPI. 

A variety of choices are also available for the numerical algorithms, such as different Riemann solvers and spatial reconstruction methods.

The relevant references are:

* Gardiner & Stone 2005, JCP, 205, 509  (2D JCP Method)
* Gardiner & Stone 2007, JCP, 227, 4123 (3D JCP Method)
* Stone et al. 2008, ApJS, 178, 137 (Method)
* Stone & Gardiner 2009, NewA, 14, 139 (van Leer Integrator)
* Skinner & Ostriker 2010, ApJ, 188, 290 (Cylindrical Integrator)
* Stone & Gardiner 2010, ApJS, 189, 142 (Shearing Box Method)


.. automodule:: amuse.community.athena.interface
    
    .. autoclass:: Athena
    
        .. autoparametersattribute:: parameters



.. _capreole:

capreole
--------

Capreole is a grid-based astrophysical hydrodynamics code developed by Garrelt Mellema. 
It works in one, two dimensions, and three spatial dimensions and is programmed in 
Fortran 90. It is parallelized with MPI. For the hydrodynamics it relies on the 
Roe-Eulderink-Mellema (REM) solver, which is an approximate Riemann solver for arbitrary
metrics. It can solve different hydrodynamics problems. Capreole has run on single 
processors, but also on massively parallel systems (e.g. 512 processors on a BlueGene/L).

The reference for Capreole (original version):

* Mellema, Eulderink & Icke 1991, A&A 252, 718



.. automodule:: amuse.community.capreole.interface
    
    .. autoclass:: Capreole
    
        .. autoparametersattribute:: parameters


.. _fi:

fi
--

FI is a parallel TreeSPH code for galaxy simulations. Extensively
rewritten, extended and parallelized, it is a development from code from
Jeroen Gerritsen and Roelof Bottema, which itself goes back to Treesph.

The relevant references are:

* Hernquist \& Katz 1989, ApJS 70, 419
* Gerritsen \& Icke 1997, A&A 325, 972
* Pelupessy, van der Werf & Icke 2004, A&A 422, 55
* Pelupessy, PhD thesis 2005, Leiden Observatory


.. automodule:: amuse.community.fi.interface
    
    .. autoclass:: Fi
    
        .. autoparametersattribute:: parameters


.. _gadget2:

gadget2
-------

GADGET-2 computes gravitational forces with a hierarchical tree 
algorithm (optionally in combination with a particle-mesh 
scheme for long-range gravitational forces, currently not 
supported from the AMUSE interface) and represents fluids by 
means of smoothed particle hydrodynamics (SPH). The code can 
be used for studies of isolated systems, or for simulations 
that include the cosmological expansion of space, both with 
or without periodic boundary conditions. In all these types 
of simulations, GADGET follows the evolution of a self-
gravitating collisionless N-body system, and allows gas 
dynamics to be optionally included. Both the force computation 
and the time stepping of GADGET are fully adaptive, with a 
dynamic range which is, in principle, unlimited. 

The relevant references are:

* Springel V., 2005, MNRAS, 364, 1105  (GADGET-2)
* Springel V., Yoshida N., White S. D. M., 2001, New Astronomy, 6, 51  (GADGET-1)


.. automodule:: amuse.community.gadget2.interface
    
    .. autoclass:: Gadget2
    
        .. autoparametersattribute:: parameters


.. _Radiative-Transfer:

Radiative Transfer
~~~~~~~~~~~~~~~~~~

* SimpleX_ (Delaunay triangulation based)


SimpleX
-------

SimpleX computes the transport of radiation on an irregular grid composed of
the Delaunay triangulation of a particle set. Radiation is transported along
the vertices of the triangulation. The code can be considered as a particle
based radiative transfer code: in this case particles sample the gas
density, but can be both absorbers and sources of radiation. Calculation
time with SimpleX scales linearly with the number of particles. At the moment
the code calculates the transport of ionizing radiation in the grey (one
frequency) approximation. It is especially well suited to couple with SPH
codes. 


Specifics
#########

* particle sets send to SimpleX must have attributes x [pc], y[pc], z[pc], 
  rho [amu/cm**3], flux [s**-1] and xion [none].
* care must be taken that the particle sets fit in the box_size
* the default hilbert_order should work for most particle distributions

.. automodule:: amuse.community.simplex.interface
    
    .. autoclass:: SimpleX
    
        .. autoparametersattribute:: parameters


References:

* Paardekooper J.-P., 2010, PhD thesis, University of Leiden
* Paardekooper J.-P., Kruip, C. J. H., Icke V., 2010, A&A, 515, 79 (SimpleX2)
* Ritzerveld, J., & Icke, V. 2006, Phys. Rev. E, 74, 26704 (SimpleX)


Work in progress
Installing the prerequisites
============================

For a full AMUSE installation, you will need to install some further dependencies that can be installed via your package manager - e.g. apt or yum on Linux; macports or homebrew on macOS.

Ubuntu
******

You can choose between openmpi and mpich as desired, both work with AMUSE. Please do not install both!
In the examples below we choose GCC-7 as the compiler, but more recent versions of GCC will also work.

* For openmpi:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      libopenmpi-dev openmpi-bin \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      libblas-dev liblapack-dev \
      python3-venv python3-pip git

* For mpich:

.. code-block:: sh

    sudo apt-get install build-essential gfortran python3-dev \
      mpich libmpich-dev \
      libgsl-dev cmake libfftw3-3 libfftw3-dev \
      libgmp3-dev libmpfr6 libmpfr-dev \
      libhdf5-serial-dev hdf5-tools \
      libblas-dev liblapack-dev \
      python3-venv python3-pip git


macOS
*****


On macOS, you will first need to install Xcode. You can do so via the app store.
In macOS Big Sur and later, you may have to add the following line to your .bashrc or .zshrc profile:

.. code-block:: sh

    export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

In this section we assume a default macOS installation (up to Big Sur) with MacPorts, but other methods (such as Homebrew) will also work.

You can choose between openmpi and mpich as desired, both work with AMUSE. 
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.
In the examples below we choose GCC 9 as the compiler, but other versions of GCC should also work.

* For openmpi:

.. code-block:: sh

    sudo port install gcc9 openmpi-gcc9 hdf5 gsl cmake gmp mpfr fftw-3 +gcc9 openblas lapack
    sudo port install python39
    sudo port select --set mpi openmpi-gcc9-fortran
    sudo port select --set gcc mp-gcc9
    sudo port select --set python3 python39

* For mpich:

.. code-block:: sh

    sudo port install gcc9 mpich-gcc9 hdf5 gsl cmake gmp mpfr fftw-3 +gcc9 openblas lapack
    sudo port install python39
    sudo port select --set mpi mpich-gcc9
    sudo port select --set gcc mp-gcc9
    sudo port select --set python3 python39



Installing AMUSE
================


After installing the prerequisites, you can install AMUSE.
Optionally, first create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don’t need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory):

.. code-block:: sh

    python3 -m venv Amuse-env

When the environment is created, you can activate it with:

.. code-block:: sh

    . Amuse-env/bin/activate

You may want to make an alias for this, e.g.:

.. code-block:: sh

    alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'

From this point, your prompt will have ‘Amuse-env’ in front of it, so you will always know when you’re in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE:

.. code-block:: sh

    pip install --upgrade pip

    pip install numpy docutils mpi4py h5py wheel

Probably, you’ll want to install these Python modules too:

.. code-block:: sh

    pip install scipy astropy jupyter pandas seaborn matplotlib

Now we can finally install AMUSE itself.
This is done easiest via pip:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse

If you only require a subset of AMUSE, you can install any of the individual packages as such:

.. code-block:: sh

    pip install amuse-framework
    pip install amuse-$(community_code_name)



Re-installation notes and troubleshooting pip installs
******************************************************

The packages installed with pip are distributed as source packages that must be compiled against the libraries
installed on your local machine. After compilation pip saves a binary package version in its cache.
In case of problems with the AMUSE installation using pip or if the environment changes it may be necessary to clean the pip cache (e.g. at ```~/.cache/pip```). In addition, the cache can be disabled using the ```--no-cache-dir``` option. the ```--no-build-isolation``` may also be tried in case the virtualenv has all the prerequisites, but the build still fails.
The ```--no-clean``` pip install option preserves the build directory for debugging purposes (The actual directory is reported 
in verbose mode ```-v```). 



Development build
*****************

Alternatively, you can install amuse as a development build, which allows you to modify the source code. It is potentially also more convenient when encountering issues with installation of specific codes as the build.log file in the root directory of the repository contains the error logs of the installation process.

Installation can also be handled through pip by executing (in the root of a clone of the repository)

.. code-block:: sh

    pip install -e .

after this the codes need to be build:

.. code-block:: sh

    python setup.py develop_build

individual codes can be build with:

.. code-block:: sh

    make {code}.code

with {code} the name of the code in lower case. 
Getting started with AMUSE
==========================

AMUSE is based on python, so if you’re new to Python, you’ll find the official `Python documentation <https://docs.python.org/3/>`_ a valuable resource. Like with Python, there are basically two ways to use AMUSE. Firstly, directly via the interactive (Python) command line:

.. code-block:: sh

    > python
    Python 3.8.0 (default, Nov  3 2019, 10:55:54) 
    [Clang 11.0.0 (clang-1100.0.33.8)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 
    >>> quit()

Secondly, by writing (Python) scripts. Suppose you wrote the following script myscript.py, and saved it in the current working directory:

.. code-block:: python

    from amuse.units.units import *
    from amuse.units import constants

    def convert_to_freq(wavelengths = [355.1, 468.6, 616.5, 748.1, 893.1] | nano(m)):
        """
        This function converts wavelength to frequency, using the speed of
        light in vacuum.
        """
        print("The speed of light in vacuum:", constants.c)
        print("wavelength -->  frequency")
        for wavelength in wavelengths:
            print(wavelength, "  --> ", (constants.c/wavelength).as_quantity_in(giga(Hz)))

Then this script can be executed from the AMUSE interactive command line:

.. code-block:: python

    >>> import myscript
    >>> help(myscript) # Tells you what myscript can do, ...
    >>>    # ... for example that it has a function to convert wavelength to frequency.
    >>> myscript.convert_to_freq()
    The speed of light in vacuum: 299792458.0 m * s**-1
    wavelength -->  frequency
    355.1 nm   -->  844247.98085 GHz
    468.6 nm   -->  639761.967563 GHz
    616.5 nm   -->  486281.359286 GHz
    748.1 nm   -->  400738.481486 GHz
    893.1 nm   -->  335676.24902 GHz
    >>> from amuse.units.units import *
    >>> myscript.convert_to_freq([21.0, 18.0, 6.0] | cm)
    The speed of light in vacuum: 299792458.0 m * s**-1
    wavelength -->  frequency
    21.0 cm   -->  1.42758313333 GHz
    18.0 cm   -->  1.66551365556 GHz
    6.0 cm   -->  4.99654096667 GHz
    >>> quit()

You can also run scripts directly from the terminal prompt. Calling python with a file name argument will execute the file. For this you need to add the following line to your script, telling the script which of its functions to call when executed:

.. code-block:: python

    if __name__ == '__main__':
        convert_to_freq()

Your script can now be executed directly from the terminal prompt:

.. code-block:: sh

> python myscript.py

    The speed of light in vacuum: 299792458.0 m \* s\*\*-1 wavelength --\>
    frequency 355.1 nm --\> 844247.98085 GHz 468.6 nm --\> 639761.967563
    GHz 616.5 nm --\> 486281.359286 GHz 748.1 nm --\> 400738.481486 GHz
    893.1 nm --\> 335676.24902 GHz



Example interactive session
===========================

This is an example of an interactive session with AMUSE, showing how the interface to a typical (gravitational dynamics) legacy code works. Using the Barnes & Hut Tree code, the dynamics of the Sun-Earth system is solved. This two-body problem is chosen for simplicity, and is, of course, not exactly what a Tree code normally is used for. First we import the necessary AMUSE modules.

.. code-block:: python

    >>> from amuse.community.bhtree.interface import BHTree
    >>> from amuse.datamodel import Particles
    >>> from amuse.units import nbody_system
    >>> from amuse.units import units

Gravitational dynamics legacy codes usually work with `N-body <https://en.wikipedia.org/wiki/N-body_units>` units internally. We have to tell the code how to convert these to the natural units of the specific system, when creating an instance of the legacy code class.

.. code-block:: python

    >>> convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    >>> instance = BHTree(convert_nbody)

Now we can tell the instance to change one of its parameters, before it initializes itself:

.. code-block:: python

    >>> instance.parameters.epsilon_squared = 0.001 | units.AU**2

Then we create two particles, with properties set to those of the Sun and the Earth, and hand them over to the BHTree instance.

.. code-block:: python

    >>> stars = Particles(2)
    >>> sun = stars[0]
    >>> sun.mass = 1.0 | units.MSun
    >>> sun.position = [0.0,0.0,0.0] | units.m
    >>> sun.velocity = [0.0,0.0,0.0] | units.m / units.s
    >>> sun.radius = 1.0 | units.RSun
    >>> earth = stars[1]
    >>> earth.mass = 5.9736e24 | units.kg
    >>> earth.radius = 6371.0 | units.km 
    >>> earth.position = [1.0, 0.0, 0.0] | units.AU
    >>> earth.velocity = [0.0, 29783, 0.0] | units.m / units.s
    >>> instance.particles.add_particles(stars)

We need to setup a channel to copy values from the code to our model in python:

.. code-block:: python

    >>> channel = instance.particles.new_channel_to(stars)

Now the model can be evolved up to a specified end time. The current values of the particles are retieved from the legacy code by using copy from the channel.

.. code-block:: python

    >>> print(earth.position[0])
    149597870691.0 m
    >>> print(earth.position.in_(units.AU)[0])
    1.0 AU
    >>> instance.evolve_model(1.0 | units.yr)
    >>> print(earth.position.in_(units.AU)[0])  # This is the outdated value! (should update_particles first)
    1.0 AU
    >>> channel.copy()
    >>> print(earth.position.in_(units.AU)[0])
    0.999843742682 AU
    >>> instance.evolve_model(1.5 | units.yr)
    >>> channel.copy()
    >>> print(earth.position.in_(units.AU)[0])
    -1.0024037469 AU

It’s always a good idea to clean up after you’re finished:

.. code-block:: python

    >>> instance.stop()

Example scripts
===============

In the `test/examples <https://github.com/amusecode/amuse/tree/master/examples>` subdirectory several example scripts are included. They show how the different legacy codes can be used. One such example is `test_HRdiagram_cluster.py <https://github.com/amusecode/amuse/blob/master/examples/applications/test_HRdiagram_cluster.py>`. It has several optional arguments. The example script can be executed from the AMUSE command line as well as from the terminal prompt (in the latter case use -h to get a list of the available command line options):

.. code-block:: python

    >>> import test_HRdiagram_cluster
    >>> test_HRdiagram_cluster.simulate_stellar_evolution()
    The evolution of  1000  stars will be  simulated until t= 1000.0 Myr ...
    Using SSE legacy code for stellar evolution.
    Deriving a set of  1000  random masses following a Salpeter IMF between 0.1 and 125 MSun (alpha = -2.35).
    Initializing the particles
    Start evolving...
    Evolved model successfully.
    Plotting the data...
    All done!
    >>> from amuse.units.units import *
    >>> test_HRdiagram_cluster.simulate_stellar_evolution(end_time=5000 | Myr)
    The evolution of  1000  stars will be  simulated until t= 5000 Myr ...
    ...

.. code-block:: python

    > python test_HRdiagram_cluster.py -h
    Usage: test_HRdiagram_cluster.py [options]

    This script will generate HR diagram for an 
    evolved cluster of stars with a Salpeter mass 
    distribution.

    Options:
      -h, --help            show this help message and exit
    ...
    > python test_HRdiagram_cluster.py
    The evolution of  1000  stars will be  simulated until t= 1000.0 Myr ...
    ...



If instead of “Plotting the data…” the script printed “Unable to produce plot: couldn’t find matplotlib.”, this probably means you do not have Matplotlib installed. See the subsection on Matplotlib_ below.

Matplotlib
**********

Matplotlib is a python plotting library which produces publication quality figures. Many of the AMUSE example scripts use this library to produce graphical output. If you would like to take advantage of this library, get it from `https://matplotlib.org/ <https://matplotlib.org/>` and install it in the Python site-packages directory. For your own work, it is of course also possible to print the required output to the terminal and use your favourite plotting tool to make the figures.
Installing on Suse Linux
========================

Installing on OpenSuse 11
~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a normal desktop install of OpenSuse 11. Not
all packages are available in the default OpenSuse package repository.
We recommend to add the **Packman Repository** to the list of 
configured software repositories (To do so, open Yast and go to 
*Software Repositories*).

Python
------
OpenSuse comes with python2.6 pre-installed, you can check if
python is installed by doing:

.. code-block:: sh

	> python --version
	Python 2.6

If this failes with an error or a version before 2.6, please install 
python first(the package is called ``python``). You also need 
the ``python-devel`` development package.
To install it, do::

    > sudo zypper install python-devel
    

GCC
---
By default, OpenSuse does not install a fortran 90 or a C++ compiler. We
suggest using gfortran and g++. These compilers are installed with
the ``gcc``, ``gcc-c++`` and the ``gcc-fortran`` packages. 
To install these, do::

    > sudo zypper install gcc gcc-c++ gcc-fortran

MPI2
----
The Packman Repository provides an OpenMPI package.
To install the openmpi packages, do::

    > sudo zypper install openmpi openmpi-devel

Unfortunately the openmpi installation does not work out
of the box, you need to set the  **LD_LIBRARY_PATH** variable
and edit a configuration file first.

Setting the LD_LIBRARY_PATH
****************************

The LD_LIBRARY_PATH must be set so that mpi4py can find the
openmpi libraries. To set the variable we must first find out
where the openmpi libs can be found, to do so execute::

    > mpicxx -showme:link
    -pthread -L/usr/lib/mpi/gcc/openmpi/lib -lmpi_cxx -lmpi 
    -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
    

We need to set LD_LIBRARY_PATH variable to the path after the **-L**
in the output (so in this example case '/usr/lib/mpi/gcc/openmpi/lib',
this may be a different path if you system is 64-bits or if the
opensuse version is different).

In bash do::
    
    > export LD_LIBRARY_PATH=/usr/lib/mpi/gcc/openmpi/lib
    
We recommend you add this line to your '.bashrc' file so that
the variable is set correctly for all sessions. If you have a
C shell you need to do a *setenv* and edit the .cshrc file.

Editing the configuration file
*******************************

It seems that the default openmpi installation has some problems
with loading an LDAP library. To check if your installation has 
this problem do::

    > python -c "from mpi4py import MPI; print MPI.get_vendor()"
    ...
    WARNING: ....
    ...
    DAT: library load failure: libdaplscm.so.2: cannot open shared object file: No such file or directory
    ...

If you get a long list of warnings about DAT providers not found, you
need to edit the configuration file and turn off ldap. To do so, 
open an editor (as root) on the file 
**/etc/openmpi-mca-params.conf** 
and add this line to the bottom of the file::

    btl = ^udapl

After saving the file, you can rerun the python statement::

    > python -c "from mpi4py import MPI; print MPI.get_vendor()"
    ('Open MPI', (1, 2, 8))
    
    
HDF5
----
Amuse can work with HDF5 versions 1.6.* and 1.8.*. The Packman Repository
has a package with HDF5 version 1.8.1. To install it, do::

    > sudo zypper install hdf5 hdf5-devel

FFTW
-------
Some codes in AMUSE need FFTW 3, FFTW can be installed with::

    > sudo zypper install fftw3 fftw3-devel

GSL
-------
On OpenSuse (10.2 and newer), GSL can be installed with::

    > sudo zypper install gsl gsl-devel

CMake
-------
CMake is used to build EVTwin. On OpenSuse, CMake can be installed with::

    > sudo zypper install cmake

GMP
-------
GMP is required for Adaptb. On OpenSuse, GMP can be installed with::

    > sudo zypper install gmp-devel

MPFR
-------
MPFR is required for Adaptb. On OpenSuse, MPFR can be installed with::

    > sudo zypper install libmpfr4 mpfr-devel

Python packages in Fedora
-------------------------
Fedora comes with python packages for numpy. You also need 
the setuptools package to be able to install the other python
packages. To install these, do::

    > sudo zypper install python-numpy \
        python-setuptools python-setuptools-devel

Python packages with easy_install
---------------------------------
The  ``nose``, ``mpi4py``, ``h5py`` and ``docutils`` can be 
installed with the ``easy_install`` command::

    > sudo easy_install nose
    > sudo easy_install mpi4py
    > sudo easy_install h5py
    > sudo easy_install docutils 


Installing on Fedora 18
~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a basic install of Fedora 18 installation.


All in One
----------
The prerequisites can be installed with a couple of commands
on Fedora.

For mpich2 do::
	
	> sudo yum install make gcc gcc-c++ gcc-gfortran\
		cmake zlib-devel\
		mpich2 mpich2-devel\
		hdf5 hdf5-devel\
		fftw fftw-devel\
		gsl gsl-devel\
		gmp gmp-devel\
		mpfr mpfr-devel\
		python-nose numpy numpy-f2py\
		h5py\
		python-setuptools python-setuptools-devel\
		mpi4py-mpich2\
		python-matplotlib
    


.. note::
    This line will also install `matplotlib`, this package is
    used for all plotting in AMUSE. If you do not need any plotting
    you can leave it out.

After installing mpich2, you need to activate it using the 'module'
command::

	> module load mpi/mpich2-$(uname -i)
	
.. note::

    We recommend to put the module activation script
    in your .bashrc or .cshrc file.

For openmpi do::

	> sudo yum install make gcc gcc-c++ gcc-gfortran\
		cmake zlib-devel\
		openmpi openmpi-devel\
		hdf5 hdf5-devel\
		fftw fftw-devel\
		gsl gsl-devel\
		gmp gmp-devel\
		mpfr mpfr-devel\
		python-nose numpy numpy-f2py\
		h5py\
		python-setuptools python-setuptools-devel\
		mpi4py-openmpi\
		python-matplotlib
	
After installing openmpi, you need to activate it using the 'module'
command::

	> module load mpi/openmpi-$(uname -i)


.. note::

    On Fedora you can install both mpich2 and openmpi, the module
    command will keep manage these separate installation, so
    no conflict will exists. If you change between implementation, you
    will need to recompile the amuse community codes with::
    
	> make clean; make
	


Installing on Fedora 11
~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a live-cd install of Fedora 11 installation.

Python
------
Fedora comes with python2.6 pre-installed, you can check if
python is installed by doing:

.. code-block:: sh

	> python --version
	Python 2.6.2

If this fails with an error or a version before 2.6, please install 
python first(the package is called ``python``). You also need 
the ``python-devel`` development package.
To install it, do::

    > sudo yum install python-devel
    

GCC
---
By default, Fedora does not install a fortran 90 or a C++ compiler. We
suggest using gfortran and g++. These compilers are installed with
the ``gcc``, ``gcc-c++`` and the ``gcc-gfortran`` packages. 
To install these, do::

    > sudo yum install gcc gcc-c++ gcc-gfortran

MPI2
----
Fedora comes with packages for MPICH2 and Openmpi.

To install MPICH2, do::
    
    > sudo yum install mpich2 mpich2-devel

If you prefer OpenMpi over MPICH2, you can install openmpi
from the Fedora yum database. 
To install the openmpi packages, do::

     > sudo yum install openmpi openmpi-devel

HDF5
----
Amuse can work with HDF5 versions 1.6.* and 1.8.3. Fedora 11 has a package
with HDF5 version 1.8.3. To install it, do::

    > sudo yum install hdf5 hdf5-devel

FFTW
-------
On Fedora, FFTW can be installed with::

    > sudo yum install fftw fftw-devel

GSL
-------
On Fedora, GSL can be installed with::

    > sudo yum install gsl gsl-devel

CMake
-------
CMake is used to build EVTwin. On Fedora, CMake can be installed with::

    > sudo yum install cmake

GMP
-------
GMP is required for Adaptb. On Fedora, GMP can be installed with::

    > sudo yum install gmp

MPFR
-------
MPFR is required for Adaptb. On Fedora, MPFR is currently included in the gmp 
package. So, if you have not already done so, MPFR can be installed with::

    > sudo yum install gmp

Python packages in Fedora
-------------------------
Fedora comes with python packages for nose and numpy. You also need 
the setuptools package to be able to install the ``mpi4py`` and ``h5py`` 
software. To install these , do::

    > sudo yum install python-nose numpy numpy-f2py \
        python-setuptools python-setuptools-devel

Python packages with easy_install
---------------------------------
The ``mpi4py``, ``h5py`` and ``docutils`` can be 
installed with the ``easy_install`` command::

    > sudo easy_install mpi4py
    > sudo easy_install h5py
    > sudo easy_install docutils 
Installing on Ubuntu
====================

Installing on Ubuntu (up to date for version 18.04)
---------------------------------------------------

In this section we assume a default Ubuntu desktop installation.

All
---
The prerequisites can be installed with a couple of commands
on Ubuntu. The only choice to make is between openmpi and mpich2. 

For openmpi do::

	> sudo apt-get install build-essential gfortran python-dev \
	  libopenmpi-dev openmpi-bin \
	  libgsl-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr6 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	>  [sudo] pip install mpi4py
	or alternatively setuptools easy_install (deprecated):
	>  [sudo] easy_install mpi4py


For mpich do::
	
	> sudo apt-get install build-essential gfortran python-dev \
	  mpich libmpich-dev \
	  libgsl-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr6 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  python-nose python-numpy python-setuptools python-docutils \
	  python-h5py python-setuptools git
	
	>  [sudo] pip install mpi4py
	or alternatively setuptools easy_install (deprecated):
	>  [sudo] easy_install mpi4py

.. note::
	
	Please make sure not to install mpich2 and openmpi together. 
	When both openmpi and mpich2 are installed strange errors
	will occur and AMUSE will not work. If you see both installed
	please remove both and install one.
Installing on macOS High Sierra with Python 3.6
===============================================

In this section we assume a default macOS installation with macports installed.

Installing prerequisites
------------------------
These prerequisites are essential for building AMUSE and the community codes.
They can be installed system-wide with the commands below.
You can choose between openmpi and mpich as desired, both work with AMUSE.
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.

For openmpi do::

  > sudo port install gcc7 openmpi-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python36 py36-virtualenv
  > sudo port select --set mpi openmpi-gcc7-fortran
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python3 python36
  > sudo port select --set virtualenv virtualenv36
  
For mpich do::
	
  > sudo port install gcc7 mpich-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python36 py36-virtualenv
  > sudo port select --set mpi mpich-gcc7
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python2 python27
  > sudo port select --set virtualenv virtualenv36

.. note:
  Please make sure not to install mpich and openmpi together. 
  When both are installed strange errors will occur and AMUSE will not work.
  If you have both installed please first remove both and then install one.
  
Installing AMUSE
----------------

First, create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don't need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory)::

  > python3 -m venv Amuse-env
  
When the environment is created, you can activate it with::

  > . Amuse-env/bin/activate

You may want to make an alias for this, e.g.::

  > alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'
  
From this point, your prompt will have 'Amuse-env' in front of it, so you will always know when you're in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE::

  > pip install numpy nose docutils mpi4py h5py
  
Probably, you'll want to install these Python modules too::

  > pip install scipy astropy jupyter pandas seaborn
  
Now we can finally install AMUSE itself.
First, download AMUSE or preferably make a git clone (in a desired directory)::

  > git clone https://github.com/amusecode/amuse.git

Then, change to the AMUSE directory and run configure, enabling optional GPU features if present/required::

  > cd amuse
  > ./configure [--enable-cuda] [--enable-sapporo2]

Finally, build and install AMUSE, with optionally downloaded codes if desired::

  > [export DOWNLOAD_CODES=1]
  > python setup.py install
 
.. note:
  The part below does not currently work in a Python 3 environment. Please skip it for now.

Optionally, to test if your setup was successful, run (this will take a long time)::

  > python setup.py test
.. _prerequisite-label:


Installation of the prerequisite software
=========================================


.. toctree::
   :maxdepth: 1
   
   install-prerequisites-ubuntu
   install-prerequisites-osx
   install-prerequisites-arch
   install-prerequisites-fedora
   install-prerequisites-centos
   install-prerequisites-suse

   

Before installing AMUSE several software packages must be installed. These 
software packages can be installed manually or with two prepared installation
scripts. The installation scripts will install python and the 
other prerequisites in a user directory. No "root" access is required.

These are the packages AMUSE needs:

* Python (version >= 2.6)
* Numpy (version >= 1.3.0)
* HDF (version 1.6.5 - 1.8.x)
* h5py (version >= 1.2.0)
* MPI (OpenMPI or MPICH)
* mpi4py (version >= 1.0)
* nose (version >= 0.11)
* docutils (version >= 0.6)
* FFTW (version >= 3.0)
* GSL
* CMake (version >= 2.4)
* GMP (version >= 4.2.1)
* MPFR (version >= 2.3.1)

In the first two sections (compilers_ and installation_scripts_) we explain how to use the two
installation scripts to install AMUSE. In the last section (manual_) 
we have specified the required packages with the needed version for each.

.. _compilers:

Compilers
*********

To build AMUSE from source you need to have a working  build environment.
The AMUSE build system needs a C++ and fortan 90 compiler. Please check first if you
have a working build environment on your system.

In Ubuntu you can setup the environment with (as root):

.. code-block:: sh

	apt-get install build-essential curl g++ gfortran gettext zlib1g-dev



In Fedora you can setup the environment with (as root):

.. code-block:: sh

	yum groupinstall "Development Tools" "Development Libraries"

.. _installation_scripts:

Installation scripts
~~~~~~~~~~~~~~~~~~~~

We have created two installation scripts to automate the installation of
the required packages on a LINUX and OS.X system. These scripts will
install these packages in a user directory. One script downloads and
installs python while the other script downloads and installs the libraries
and python packages. As everything is installed in a user directory these
packages can be installed even if a version of the software is already 
installed on your system. 

The scripts will download and install the software in a user directory. This
user directory must be specified with the ``PREFIX`` environment variable. Before
running the installation scripts you must set the ``PREFIX`` environment 
variable and update the path and library path. For shell (bash) you need to do:

.. code-block:: sh

	export PREFIX=~/amuse/prerequisites
	export PATH=${PREFIX}/bin:${PATH}
  	export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}


One script will download, build and install python on your system. The other 
script is written in Python and will download and install the other packages. 
Both scripts can be found in the ``doc/install`` directory. 

To start the installation do:

.. code-block:: sh

	# 1. Open a shell and go to the <doc/install> directory
	>

	# 2. Set the PREFIX, PATH and LD_LIBRARY_PATH environment variables:
  	> export PREFIX=~/amuse/prerequisites
  	> export PATH=${PREFIX}/bin:${PATH}
  	> export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

	# 3. Start the installation script for python
	> ./install-python.sh

	# 4. Start the installation script for the prerequisite packages
	> ./install.py download
  	> ./install.py install

	# 5. Update your PATH variable in your profile. 
	# Make sure the `${PREFIX}/bin` directory is the first entry in the PATH!

You should now be able to install AMUSE.

Using the installation scripts on macOS
---------------------------------------

For macOS you need to install XCode and a gfortran compiler first.
The XCode development package is available at https://developer.apple.com/xcode/ or in the App Store. 

The standard XCode release does not come with a gfortran compiler. 
Go to the `HPC Mac OS X site <http://hpc.sourceforge.net/index.php>`_ 
for a recent gfortran compiler, compatible with the XCode tools.

After installing XCode and gfortan, follow the steps described in the
previous paragraph.

.. _manual:

Manually installing the prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python
------
Python is probably already installed on your system. To check the version of python do:

.. code-block:: sh

	> python --version
	Python 2.6.2

You can download python from https://www.python.org. 

Numpy
-----
To check if numpy is installed on your system do:

.. code-block:: sh

	> python -c 'import numpy; print numpy.version.version'
	1.3.0

If this fails with an error or a version before 1.3 you need to install numpy.
You can download numpy from http://www.numpy.org. 

HDF5 library
------------
HDF5 is a data format specification. The HDF group provides a C library 
to write and access HDF files.

To check if the HDF library is installed on your system do:

.. code-block:: sh

	> h5ls -V
	h5ls: Version 1.8.3

If this fails with an error or a version before 1.6.5 you 
need to install the HDF library. 
You can download HDF from https://www.hdfgroup.org/. 

h5py
----
To access HDF5 files from python we use the ``h5py`` library.

To check if the h5py library is installed on your system do:

.. code-block:: sh

	> python -c 'import h5py; print h5py.version.version'
	1.2.0

If this fails with an error or a version before 1.2.0 you need to install h5py.
You can download h5py from https://www.h5py.org. 

docutils
--------
To check if the python docutils are installed on your system do:

.. code-block:: sh

        > python -c 'import docutils; print docutils.__version__'
        0.6

If this fails with an error or a version before 0.6 you need to install docutils.
You can download docutils from http://docutils.sourceforge.net/

MPI
---
The installed MPI framework must be MPI 2 compatible. AMUSE will work with 
MPICH or OpenMPI

MPICH
^^^^^
MPICH is a portable implementation of the MPI 2 standard.

To check if MPICH is installed on your system do:

.. code-block:: sh
    
    > mpdhelp
    
    The following mpd commands are available.  For usage of any specific one,
    invoke it with the single argument --help .

    mpd           start an mpd daemon
    mpdtrace      show all mpds in ring
    mpdboot       start a ring of daemons all at once
    mpdringtest   test how long it takes 
    ...
    
If this fails with an error you need to install MPICH or check for OpenMPI
support.
You can download MPICH from http://www.mpich.org.

OpenMPI
^^^^^^^
OpenMPI is another portable implementation of the MPI 2 standard

To check if OpenMPI is installed on your system do:

.. code-block:: sh
    
    > mpicxx -v 
    

If this fails with an error you need to install MPICH or OpenMPI
support. Most examples in the dopcumentation assume OpenMPI.
You can download OpenMPI from https://www.open-mpi.org/.

MPI4PY
------
To access MPI from python we use the ``mpi4py`` software.
To check if the mpi4py library is installed on your system do:

.. code-block:: sh

	> python -c 'import mpi4py; print mpi4py.__version__'
	1.0.0

If this fails with an error or a version before 1.0 you need to install mpi4py.
You can find mpi4py at https://mpi4py.readthedocs.io/. 

Nose
----
Nose is an extension of the python testing framework. It is used for all
unit testing in AMUSE.


To check if Nose is installed on your system do:

.. code-block:: sh
    
    > nosetests --version
    nosetests version 0.11.1
    ...
    
If this fails with an error or a version before 0.11 you need to install nose.
You can download nose from http://somethingaboutorange.com/mrl/projects/nose/. 

FFTW
----
FFTW is a C subroutine library for computing discrete Fourier transforms. To 
check for the availability of fftw on your system, you can use ``fftw-wisdom``:

.. code-block:: sh

   > fftw-wisdom --version
   fftw-wisdom tool for FFTW version 3.2.1.


You can download the FFTW library from http://www.fftw.org. 

GSL
-------
The GNU Scientific Library (GSL) is a numerical library for C and C++ 
programmers. It is free software under the GNU General Public License.
To check for the availability of GSL on your system, you can use ``gsl-config``:

.. code-block:: sh

   > gsl-config --version
   1.14


You can download GSL from http://www.gnu.org/software/gsl/. 

CMake
-------
CMake is a cross-platform, open-source build system. CMake is used to control 
the software compilation process using simple platform and compiler independent
configuration files. CMake generates native makefiles and workspaces that can 
be used in the compiler environment of your choice.
CMake is used to build EVTwin. 
To check whether you have CMake installed on your system:

.. code-block:: sh

   > cmake --version
   cmake version 2.8.2


You can download CMake from https://cmake.org/download/. 

GMP
-------
GNU MP is a library for arbitrary precision arithmetic (ie, a bignum package). 
It can operate on signed integer, rational, and floating point numeric types.
GMP is required for Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda 
Boekholt). 
The best way to check whether you have the right version of GMP installed on your 
system depends on the package manager you use, but this should always work (note 
that the library numbers do not match the release version):

.. code-block:: sh

   > locate libgmp
   /usr/lib64/libgmp.so
   /usr/lib64/libgmp.so.10
   /usr/lib64/libgmp.so.10.0.3
   
   > locate gmp.h
   /usr/include/gmp.h
   
   > grep GNU_MP_VERSION /usr/include/gmp.h
   #define __GNU_MP_VERSION 5
   #define __GNU_MP_VERSION_MINOR 0
   #define __GNU_MP_VERSION_PATCHLEVEL 3


You can download GMP from https://gmplib.org. 

MPFR
-------
The MPFR library is a C library for multiple-precision floating-point 
computations with correct rounding.
MPFR is required for Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda 
Boekholt). 
The best way to check whether you have the right version of MPFR installed on your 
system depends on the package manager you use, but this should always work (note 
that the library numbers do not match the release version):

.. code-block:: sh

   > locate libmpfr
   /usr/lib64/libmpfr.so
   /usr/lib64/libmpfr.so.4
   /usr/lib64/libmpfr.so.4.1.0
   
   > locate mpfr.h
   /usr/include/mpfr.h
   
   > grep MPFR_VERSION /usr/include/mpfr.h
   #define MPFR_VERSION_MAJOR 3
   #define MPFR_VERSION_MINOR 1
   #define MPFR_VERSION_PATCHLEVEL 0

You can download MPFR from https://www.mpfr.org. 





    
    





    
    

    
author: Arjen van Elteren (vanelteren@strw.leidenuniv.nl)
date: 2010/09/22





.. _configuration-label:

=================
Configuring AMUSE
=================

Introduction
~~~~~~~~~~~~
In AMUSE a configuration script is used to for two purposes; to run 
on different operating systems and to set compile time options. The 
AMUSE framework has been built on Linux, AIX and OS X systems, it 
also runs on Windows. AMUSE can be configured to run with or without 
MPI, GPU (CUDA) and openmp. In this document we will provide a short 
overview of the configuration options and their effects.


Basic
~~~~~
The basic configuration of AMUSE uses MPI as the communication 
channel, does not build any GPU enabled codes (or GPU enabled 
versions) and uses openmp if available.  The configuration script
can be run as::

    > ./configure
    
To get a list of options and important environment variables run 
```configure``` with the help flag::

    > ./configure --help
    
A very important variable for the configuration script is the 
location of the python executable. The python executable is searched 
for in the PATH and you can override it by setting the ```PYTHON``` 
environment variable::

    > ./configure PYTHON=/path/to/libraries/python
 
The configuration script will look for dependent libraries in 
default locations of the system and, if defined, also in directories 
under the ```PREFIX``` environment variable. If you installed the 
prerequisites with the AMUSE installation scripts (see 
:doc:`howto-install-prerequisites`), the configuration script should 
find all the packages installed. For most libraries the 
```PREFIX/lib``` or ```PREFIX/lib64``` is searched before the system 
path. You can override the ```PREFIX``` environment variable::
   
    > ./configure PREFIX=/path/to/libraries/root

.. _configuration-gpu-label:

GPU
~~~

Currently all codes in AMUSE capable of using the GPU are based 
on CUDA. To run these codes you will need CUDA libraries and drivers.
Once these have been installed you can run configure like so::

    > ./configure --enable-cuda

The configuration script will look for the ```nvcc``` compiler and the 
cuda libraries, in well known paths. Unfortunately it often will not find
the cuda tools and you have to specify the some environment variables or
configuration options. 

If the configuration script cannot find the nvcc compiler (or if it 
finds the wrong one) you can specify the nvcc compiler with the 
```NVCC``` environment variable::

    > ./configure --enable-cuda NVCC=/path/to/nvcc
    
The configure script also searches for the nvcc compiler in the 
```$CUDA_TK/bin``` directory::

    > ./configure --enable-cuda CUDA_TK=/opt/nvidia
    
The configure script looks for cuda and cudart libraries in 
```$NVCC/../lib``` or ```$NVCC/../lib64```, if your libraries cannot 
be found there you can override the library path with::

    > ./configure --enable-cuda --with-cuda-libdir=/path/to/cuda/lib
    
Using ```--with-cuda-libdir``` will always override the local 
search paths and should also work if you have an old version of cuda 
in ```/usr/lib```.

Finally, if all else fails, you can edit the ```config.mk``` file 
after configure has finished. The important variables in the file are:

 * CUDA_ENABLED, valid values are "yes" or "no".
 * NVCC, absolute path to the nvcc executable.
 * CUDA_TK, directory of the cuda toolkit installation
 * CUDA_LIBS, library flags the add in the linking stage (-L/path -lcuda -lcudart)
 
Please remember that the ```config.mk``` file is overwritten 
every time configure is run.

Sapporo library version
-----------------------
 
The Sapporo library will be build when CUDA is enabled. The Sapporo 
library implements the GRAPE6 API on GPU hardware. AMUSE is shipped 
with two versions of the Sapporo library:

 * An older version ```sapporo_light``` that runs on most CUDA devices but is not maintained any longer
 * The latests version ```sapporo``` that runs on modern GPU hardware. This version should also run on 
   OpenCL devices but this is still a work in progress.

By default AMUSE will use the older ```sapporo_light``` version, to enable
the latests version do:

    > ./configure --enable-cuda --enable-sapporo2
Installing on Ubuntu version 18.04 with Python 3.6
==================================================

In this section we assume a default Ubuntu desktop installation.

Installing prerequisites
------------------------
These prerequisites are essential for building AMUSE and the community codes.
They can be installed system-wide with the commands below.
You can choose between openmpi and mpich as desired, both work with AMUSE.

For openmpi do::

  > sudo apt-get install build-essential gfortran python3-dev \
	  libopenmpi-dev openmpi-bin \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  git

For mpich do::
	
  > sudo apt-get install build-essential gfortran python3-dev \
	  mpich libmpich-dev \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  git

.. note:
  Please make sure not to install mpich and openmpi together. 
  When both are installed strange errors will occur and AMUSE will not work.
  If you have both installed please first remove both and then install one.

  
Installing AMUSE
----------------

First, create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don't need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory)::

  > python3 -m venv Amuse-env
  
When the environment is created, you can activate it with::

  > . Amuse-env/bin/activate

You may want to make an alias for this, e.g.::

  > alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'
  
From this point, your prompt will have 'Amuse-env' in front of it, so you will always know when you're in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE::

  > pip install numpy nose docutils mpi4py h5py
  
Probably, you'll want to install these Python modules too::

  > pip install scipy astropy jupyter pandas seaborn
  
Now we can finally install AMUSE itself.
First, download AMUSE or preferably make a git clone (in a desired directory)::

  > git clone https://github.com/amusecode/amuse.git

Then, change to the AMUSE directory and run configure, enabling optional GPU features if present/required::

  > cd amuse
  > ./configure [--enable-cuda] [--enable-sapporo2]

Finally, build and install AMUSE, with optionally downloaded codes if desired::

  > [export DOWNLOAD_CODES=1]
  > python setup.py install

.. note:
  The part below does not currently work in a Python 3 environment. Please skip it for now.

Optionally, to test if your setup was successful, run (this will take a long time)::

  > python setup.py test
Installing on macOS High Sierra with Python 2.7
===============================================

In this section we assume a default macOS installation with macports installed.

Installing prerequisites
------------------------
These prerequisites are essential for building AMUSE and the community codes.
They can be installed system-wide with the commands below.
You can choose between openmpi and mpich as desired, both work with AMUSE.
Please make sure to set the compilers installed here as default, as it will greatly simplify things later on.

For openmpi do::

  > sudo port install gcc7 openmpi-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python27 py27-virtualenv
  > sudo port select --set mpi openmpi-gcc7-fortran
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python2 python27
  > sudo port select --set virtualenv virtualenv27
  
For mpich do::
	
  > sudo port install gcc7 mpich-gcc7 hdf5 gsl cmake gmp mpfr fftw-3 +gcc7
  > sudo port install python27 py27-virtualenv
  > sudo port select --set mpi mpich-gcc7
  > sudo port select --set gcc mp-gcc7
  > sudo port select --set python2 python27
  > sudo port select --set virtualenv virtualenv27

.. note:
  Please make sure not to install mpich and openmpi together. 
  When both are installed strange errors will occur and AMUSE will not work.
  If you have both installed please first remove both and then install one.
  
Installing AMUSE
----------------

First, create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don't need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory)::

  > virtualenv Amuse-env
  
When the environment is created, you can activate it with::

  > . Amuse-env/bin/activate

You may want to make an alias for this, e.g.::

  > alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'
  
From this point, your prompt will have 'Amuse-env' in front of it, so you will always know when you're in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE::

  > pip install numpy nose docutils mpi4py h5py
  
Probably, you'll want to install these Python modules too::

  > pip install scipy astropy jupyter pandas seaborn
  
Now we can finally install AMUSE itself.
First, download AMUSE or preferably make a git clone (in a desired directory)::

  > git clone https://github.com/amusecode/amuse.git

Then, change to the AMUSE directory and run configure, enabling optional GPU features if present/required::

  > cd amuse
  > ./configure [--enable-cuda] [--enable-sapporo2]

Finally, build and install AMUSE, with optionally downloaded codes if desired::

  > [export DOWNLOAD_CODES=1]
  > python setup.py install
  
Optionally, to test if your setup was successful, run (this will take a long time)::

  > python setup.py test
Installing on macOS
*******************

Installing on macOS with MacPorts 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a clean MacPorts installation. The MacPorts build
system will build most packages from source so installation may take a while.
The packages in MacPorts support different *variants*, each *variant* is built
differently. Below, with all installation commands we will specify the variant
where needed. AMUSE is tested with gcc versions 4.3 up to 8. Below, we will use
gcc 7.

.. note::
    
    If you want to use a different fortran compiler (ifort), you are better 
    off using the **install.py** script in the **doc/install** directory.

.. note::

    Make sure you have a recent MacPorts installation. Do this by running 
    
    > sudo port selfupdate

    before installing.
    
.. note::
    
    If you are unsure of your installation you can uninstall and clear the 
    packages with::
    
        port uninstall py27-docutils py27-nose py27-mpi4py py27-h5py py27-numpy hdf5 fftw-3 gsl openmpi python27
    
    To make a clean install of MacPorts, please remove the MacPorts directory
    and read the guide at:
    https://guide.macports.org/
    
All in one
----------

You can install all the packages described below in one go with::

    > sudo port install gcc7
    
    > sudo port install python27
    > sudo port install openmpi-gcc7
    > sudo port install fftw-3 +gcc7
    > sudo port install hdf5 gsl cmake gmp mpfr
    > sudo port install py27-numpy py27-h5py py27-nose py27-docutils 
    > sudo port install py27-mpi4py +openmpi
    > sudo port install py27-matplotlib

To make sure the right MacPorts compilers and python are set as default, do the
following::

    > sudo port select --set mpi openmpi-gcc7-fortran
    > sudo port select --set gcc mp-gcc7
    > sudo port select --set python python27
    > sudo port select --set nosetests nosetests27
    
After installing you will need to configure the code with the following line::

    > ./configure --with-fftw=/opt/local

.. note::

    The ``--with-fftw`` option will ensure that fftw is found in 
    ``/opt/local`` and not in any other location on macOS. Sometimes, incompatible versions of
    fftw may be installed in ``/usr/include`` or ``/usr/local/inlude``. These versions may
    have the wrong processor type (32 vs 64bits) or not contain a fortran API. For both cases
    compiling ``Fi`` will fail.      
    In case the configure script does not pick the wanted fftw directories, you
    can edit the ``config.mk`` file to point to the right version.
    

.. note::

    The ``PREFIX`` variable will make sure some support libraries for 
    community codes (gsl, gmp and mpfr) are found in ``/opt/local``.

.. note::

    Please, make sure you have no other compile variables specified
    (like CC or CXX or CFLAGS), unless you have customized MacPorts in
    any way. Some variable settings will likely give compile errors
    for the community codes. 
    
    For example, BHTree is compiled with openmpicxx and $CXX. 
    The command in the CXX variable must be compatible 
    with openmpicxx (you can do ``openmpicxx --show`` to get 
    the command openmpicxx is using)


GCC
---
By default MacPorts uses the XCode compilers, these compilers have no support
for fortran, a MacPorts gcc compiler set needs to be installed. We suggest
installing gcc 7 as the most reliable option:

.. code-block:: sh
    
    > sudo port install gcc7
    
.. note::
    
    If you have installed a different version of gcc, you need to select
    a different variant of the packages below. To select a different variant
    replace **+gcc7** with **+gcc6**, **+gcc8** or any other version
    matching your gcc installation. Note, apple-gcc versions will not work,
    these do not support fortran.

Python
------
MacPorts supports several python versions in different variants, we will install
the python27 versions

.. code-block:: sh

    > sudo port install python27
   
MPI2
----
MacPorts provides packages for mpich and openmpi. Although you can
probably install both, this is not recommended. We suggest you install
openmpi.

To install openmpi, do::

     > sudo port install openmpi +gcc7

HDF5
----
Amuse can work with HDF5 versions 1.6.*, 1.8.3 and higher. MacPorts comes
with HDF5 version 1.8.* and 1.10.* To install the most recent version, do::

    > sudo port install hdf5 

FFTW-3
------
MacPorts comes with a FFTW and FFTW-3 package, for AMUSE we need FFTW-3.
FFTW-3 can be installed with::

    > sudo port install fftw-3 +gcc7

GSL
---
GSL is used to build Gadget2, GSL can be installed with::

    > sudo port install gsl

CMake
-----
CMake is used to build EVTwin, CMake can be installed with::

    > sudo port install cmake

GMP
-------
GMP is required for Adaptb. With MacPorts, GMP can be installed with::

    > sudo port install gmp

MPFR
-------
MPFR is required for Adaptb. With MacPorts, MPFR can be installed with::

    > sudo port install mpfr


Python packages
---------------
By this point all libraries and frameworks are installed. We can now
install python packages (some depend on the installed libraries)::

    > sudo port install py27-numpy py27-h5py py27-nose py27-docutils

If you installed openmpi in the MPI2 step you need to set the
"openmpi" variant for "py27-mpi4py"::

    > sudo port install py27-mpi4py +openmpi


Matplotlib
----------
Matplotlib is not required but is highly recommended for creating graphics, 
you can install it with::

    > sudo port install py27-matplotlib
    

.. note::

    Macports will install the compilers under non standard names.
    To make sure the right MacPorts compilers and python are set as default, do
    the following::

    > sudo port select --set mpi openmpi-gcc7-fortran
    > sudo port select --set gcc mp-gcc7
    > sudo port select --set python python27
    > sudo port select --set nosetests nosetests27

    Alternatively, to use the right compilers you can specify these during the
    configure stage of AMUSE.
    
    See the output for ``configure --help`` for a list of all 
    environment variables you can set.
    
    If you installed openmpi you need to specify the mpi compilers 
    like so (replacing OPENMPICXX OPENMPICC and OPENMPIF90 with your installed compiler)::
    
        ./configure MPICXX=OPENMPICXX MPICC=OPENMPICC MPIFC=OPENMPIF90

Installation
============

.. toctree::
   :maxdepth: 1
   
   howto-install-AMUSE
   getting-started
   

===============
Obtaining AMUSE
===============

Download
--------

Go to the `Amusecode github <https://github.com/amusecode/amuse>`_ page.

Getting started
---------------

The first step in getting AMUSE to work is obtaining the AMUSE
source code. We advice you to do this even before installation of the 
prerequisite software (:ref:`prerequisite-label`). In the following 
installation instructions we assume that you will install AMUSE in a 
directory ``/amuse``. 

Releases
--------

For the official releases we provide tarballs via https://github.com/amusecode/amuse/releases.

Tarball
~~~~~~~

Obtain the tarball (e.g. release-11.2.tar.gz) from the download site and unpack it 
in the amuse directory using:

.. code-block:: sh

    > tar -xf release-11.2.tar.gz

this will make an amuse sub-directory ``amuse-release-11.2``, which we will be referring to as
the AMUSE root directory.

From here proceed by reading the  :ref:`prerequisite-label` section.

Bleeding edge
-------------

The current development version is available via git repository access 
by issuing the following command:

.. code-block:: sh

    > git clone https://github.com/amusecode/amuse.git

This will make an AMUSE root directory with the name "amuse".  
Installing on Ubuntu version 18.04 with Python 2.7
==================================================

In this section we assume a default Ubuntu desktop installation.

Installing prerequisites
------------------------
These prerequisites are essential for building AMUSE and the community codes.
They can be installed system-wide with the commands below.
You can choose between openmpi and mpich as desired, both work with AMUSE.

For openmpi do::

  > sudo apt-get install build-essential gfortran python-dev \
	  libopenmpi-dev openmpi-bin \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  git

For mpich do::
	
  > sudo apt-get install build-essential gfortran python-dev \
	  mpich libmpich-dev \
	  libgsl0-dev cmake libfftw3-3 libfftw3-dev \
	  libgmp3-dev libmpfr4 libmpfr-dev \
	  libhdf5-serial-dev hdf5-tools \
	  git

.. note:
  Please make sure not to install mpich and openmpi together. 
  When both are installed strange errors will occur and AMUSE will not work.
  If you have both installed please first remove both and then install one.

  
Installing AMUSE
----------------

First, create a virtual environment to install AMUSE and other desired Python packages in.
This ensures that you don't need root privileges and that your AMUSE environment is isolated from other system-installed packages.

To create the virtual environment, do (from a desired directory)::

  > virtualenv Amuse-env
  
When the environment is created, you can activate it with::

  > . Amuse-env/bin/activate

You may want to make an alias for this, e.g.::

  > alias amuse-env='. ~/virtualenvironments/Amuse-env/bin/activate'
  
From this point, your prompt will have 'Amuse-env' in front of it, so you will always know when you're in this virtual environment.

Now you can use pip to install the prerequisite python modules for AMUSE::

  > pip install numpy nose docutils mpi4py h5py
  
Probably, you'll want to install these Python modules too::

  > pip install scipy astropy jupyter pandas seaborn
  
Now we can finally install AMUSE itself.
First, download AMUSE or preferably make a git clone (in a desired directory)::

  > git clone https://github.com/amusecode/amuse.git

Then, change to the AMUSE directory and run configure, enabling optional GPU features if present/required::

  > cd amuse
  > ./configure [--enable-cuda] [--enable-sapporo2]

Finally, build and install AMUSE, with optionally downloaded codes if desired::

  > [export DOWNLOAD_CODES=1]
  > python setup.py install
  
Optionally, to test if your setup was successful, run (this will take a long time)::

  > python setup.py test
Installing on Arch Linux
========================

In this section we assume a default Arch Linux installation.

All
---
The prerequisites can be installed with a couple of commands
on Arch Linux. 

To install the prerequisites do (for base-devel select *all* members)::

    > sudo pacman -Syu base-devel curl gcc-fortran gettext zlib

Install python and dependencies::

    > sudo pacman -Syu python2 python2-numpy \
      hdf5 docutils openmpi
      python2-mpi4py python2-nose\
      fftw gsl cmake gmp mpfr

To install h5py, first install distribute and then run easy_install::

    > sudo pacman -Syu python2-distribute
    
    > sudo easy_install-2.7 h5py
.. _documenting:

======================
Writing documentation
======================

Getting started
===============

The documentation for AMUSE is generated from ReStructured Text 
using the Sphinx_ documentation generation tool. Sphinx version 1.0 
or later is required. You might still run into problems, so most 
developers work from the sphinx source repository (Mercurial based) 
because it is a rapidly evolving project:

.. code-block:: bash

  > hg clone http://bitbucket.org/birkenfeld/sphinx/
  > cd sphinx
  > python setup.py install

.. _Sphinx: http://www.sphinx-doc.org/en/master/

The documentation sources are found in the :file:`doc/` directory in the trunk.
To build the AMUSE documentation in html format, cd into :file:`doc/` and
do::

  make html

You can also pass a ``pdflatex`` flag to make to build a pdf, or pass no
arguments to show help information.

The output produced by Sphinx can be configured by editing the :file:`conf.py`
file located in the :file:`doc/`.


Organization of the AMUSE documentation
==========================================

The actual ReStructured Text files are kept in :file:`doc/install`,
:file:`doc/design`, :file:`doc/tutorial`. The main entry point is
:file:`doc/index.rst`. The documentation suite is
built as a single document in order to make the most effective use of cross
referencing, we want to make navigating the AMUSE documentation as easy as
possible.

Additional files can be added to the various sections by including their base
file name (the .rst extension is not necessary) in the table of contents.
Installing on RedHat (CentOS)
=============================

Installing on CentOS 6
~~~~~~~~~~~~~~~~~~~~~~~

In this section we assume a minimal CentOS 6 installation.

All
---
The prerequisites can be installed with a couple of commands
on CentOS 6. 

To install the prerequisites do (for base-devel select *all* members)::

    > sudo yum install make gcc gcc-c++ gcc-gfortran \
	cmake zlib-devel\
	openmpi openmpi-devel \
	fftw fftw-devel \
	gsl gsl-devel gmp
	
After installing openmpi, you need to activate it using the 'module'
command::
    
    > module load openmpi-$(uname -i)

.. note::

    We recommend to put the openmpi module activation script
    in your .bashrc or .cshrc file.

Install python and dependencies::

    > sudo yum install python-devel \
	docutils python-nose \
	numpy numpy-f2py\
	python-docutils

To install hdf5 and docutils first install an additional rpm forge.
For documentation see http://wiki.centos.org/AdditionalResources/Repositories/RPMForge

After installing an rpm forge do::

    > sudo yum install hdf5 hdf5-devel
    
To install h5py do::

    > sudo easy_install h5py
    
Last, you need to install mpi4py with::

    > su -
    > module load openmpi-$(uname -i)
    > easy_install mpi4py

.. note::
    
    The default CentOS sudo policy resets the environments variables and
    thereby removes the openmpi settings. So for the last step
    you cannot use ```sudo easy_install mpi4py``` but must install
    under root directly.

