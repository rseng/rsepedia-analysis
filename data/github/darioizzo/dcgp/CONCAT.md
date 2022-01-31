[![windows](https://github.com/darioizzo/dcgp/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/darioizzo/dcgp/actions/workflows/windows.yml)
[![linux](https://github.com/darioizzo/dcgp/actions/workflows/linux.yml/badge.svg?branch=master)](https://github.com/darioizzo/dcgp/actions/workflows/linux.yml)
[![osx](https://github.com/darioizzo/dcgp/actions/workflows/macos.yml/badge.svg)](https://github.com/darioizzo/dcgp/actions/workflows/macos.yml)

[![PyPI](https://img.shields.io/pypi/v/dcgpy?style=for-the-badge)](https://pypi.python.org/pypi?:action=display&name=dcgpy&version=1.0.1)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/dcgp-python?style=for-the-badge)](https://anaconda.org/conda-forge/dcgp-python)

[![Gitter](https://img.shields.io/gitter/room/esa/dcgp?logo=gitter-white&style=for-the-badge)](https://gitter.im/esa/dcgp)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/darioizzo/dcgp?style=for-the-badge)

[![DOI](https://zenodo.org/badge/38923383.svg)](https://zenodo.org/badge/latestdoi/38923383)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02290/status.svg)](https://doi.org/10.21105/joss.02290)

# dCGP
Implementation of differentiable Cartesian Genetic Programming (dCGP)

The dCGP is a development in the field of Genetic Programming that adds the information about the derivatives of the program output with respect to the input values and various parameters (weights, biases, etc..). In doing so, it enables a number of new applications currently the subject of active research.

 * The representation of deep neural networks using a dCGPANN allows the whole network to be encoded and evolved, including the connection topology, the weights, the biases, etc, .... as well as to learn the weights and biases using backpropagation.
 * The solution to boundary values problems, differential equations etc. can be encoded in a dCGP and evolved against different boundary conditions or initial conditions
 * Prime integrals of motion can be represented by a dCGP and learned
 * Symbolic regression tasks can learn ephemeral constants using backprop or even higher order methods.
 
Please quote the use of differentiable Cartesian Programming as:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable Genetic Programming." arXiv preprint arXiv:1611.04766 (2016).

Please quote the use of the software here provided as:

Izzo, Dario, Francesco Biscani (2020). dcgp: Differentiable Cartesian Genetic Programming made easy. Journal of Open Source Software, 5(51), 2290, https://doi.org/10.21105/joss.02290

Preliminary documentation can be found at http://darioizzo.github.io/dcgp/

A web based version of an early dCGP version can be found here: https://esa.github.io/dcgp-web/ thanks to Mike Heddes!

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* [Boost](http://www.boost.org/), various C++ utilities (>=1.72).
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), linear algebra library (>=3.3.0)
* [Pagmo](https://github.com/esa/pagmo2), parallel optimization library (>=2.15).
* [tbb](https://github.com/intel/tbb), lets you easily write parallel C++ programs that take full advantage of multicore performance (>=2020.1).
* [AuDi](http://darioizzo.github.io/audi/), high order automated differentiation library (>=1.8).
* [Symengine](https://github.com/symengine/symengine), symbolic manipulation of math expressions (>=0.6).

In case you plan to use the package manager ``conda`` to install the necessary dependencies, the following packages will be required: 

```
boost-cpp, eigen, pagmo-devel, tbb, audi, symengine, obake-devel
```
while we also suggest to install:
```
cmake, cxx-compiler
```

------------------------------------------------------------

After making sure the dependencies above are installed and found in your system, you may download
the latest dCGP code via git:

```bash
git clone https://github.com/darioizzo/dcgp.git
```
and configure your build using CMake. 

When done, type (in your build directory):

```bash
 make 
```
When finished, to run the tests type:

```bash
make test
```

If succesfull, you may now install cgp:

```bash
make install
```

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 

Python
------
The main functionalities of dCGP are exposed into a Python module called ``dcgpy`` which
can be installed from conda (OSx, linux and Win), pip (only linux) or by building the module.

### Installing with conda

``dcgpy`` is available in the `conda <https://conda.io/en/latest/>`__ package manager
from the `conda-forge <https://conda-forge.org/>`__ channel. A single package is available:

* `dcgp-python <https://anaconda.org/conda-forge/dcccgp-python>`__, which contains the ``dcgpy`` python module.

In order to install ``dcgpy`` via conda, you just need
to add ``conda-forge`` to the channels:

```
conda config --add channels conda-forge
conda install dcgp-python
```

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

You may test the successfull installation by running the python tests typing:

```
python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"
```

### Installing with pip (deprecated)

We also provide the pip packages (mainly for linux 64 bit architectures and versions <= 1.4.1).
Check on the `PyPi dcgpy page <https://pypi.org/project/dcgpy/>`_ if the needed package is provided.

```
pip install dcgpy
```

### Building


To build the module you need to have the Boost Python libraries installed and to activate the BUILD_DCGPY option from within CMake (and deselect BUILD_DCGP)

Check carefully what Python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the Python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and Python module after install.

When done, type (in your build directory):

```
make install
```

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.


## Other CGP libraries
### Comparison to the CGP-Library
If all below statements are true:
 * You do not care about knowing derivatives of your encoded program
 * You do not care about run-time capabilities
 * You do not care about the Python module
 * You do not care about the possibility of defining your kernel functions as complex functors (e.g. CGP expressions.)
 * You do not care about thread-safety

then you should consider using, instead, Andrew Turner's CGP-Library (http://www.cgplibrary.co.uk/files2/About-txt.html) which is, roughly, twice as fast to compute a CGP expression as it makes use of function pointers rather than a type-erased function wrapper to define the kernel functions.

## Community guidelines

### Report bugs

Please submit any bugs as an issue to https://github.com/darioizzo/dcgp/issues.

It would be greatly appreciated if the issue templates we provide can be used.

### Contributing

If anyone wishes to contribute to this project, please submit early on a pull request marked as "(WIP) - pull request title". We will be then discussing 
how to implement the desired functionality directly on the PR.

A guide on how to submit pull requests can be found on
https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request.
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment (please complete the following information):**
 - OS: [e.g. Linux]
 - Installation method: [e.g. conda, source]
 - Version: [e.g. 2.15]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]"
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# yacma

Yet another CMake modules archive.
.. _credits:

Credits 
============

Who is behind all of this? Once again, the elusive and yet so resourceful "pagmo foundation",
which is, basically, a bunch of people fooling around, not taking themselves too seriously and,
in the meantime trying to make some good progress in science.

dcgp/dcgpy is a software developed and offered to the users in a full FLOSS philosophy. We have 
created it, mostly, in our spare time and as a tool to explore some research ideas we started
nourishing since our early investigations into the evolutionary computations / optimization field.

We use dcgp in the context of our work for the European Space Agency to perform research on genetic programming, 
but its capabilities go well beyond that of a purely scientific/academic tool.

Over the years, our work attracted a wider interest. We acknowledge and thanks the support of
The Dow Corporation during some phases of the development and of the European Space Agency
for having its Advanced Concepts Team to look into Genetic Programming as a potentially game changing technology.

Main developers
^^^^^^^^^^^^^^^
Dario Izzo - European Space Agency

Francesco Biscani - Max Planck Institute for Astronomy

--------------------------------------------------

.. image:: ../sphinx/_static/dow.png
   :width: 40%

.. image:: ../sphinx/_static/esa.png
   :width: 40%



.. quickstart examples


Quick start examples
====================

C++
---

After following the :ref:`installationguide` you will be able to compile and run your first C++ dCGP program, 
put the following text into a ``getting_started.cpp`` file:

.. _getting_started_c++:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

To compile it, create also, in the same folder, a ``CmakeLists.txt`` file with the following content:

.. code-block:: cmake

    project(getting_started)

    cmake_minimum_required(VERSION 3.8)

    find_package(dcgp REQUIRED)

    add_executable(getting_started getting_started.cpp)
    target_link_libraries(getting_started Dcgp::dcgp)

    set_property(TARGET getting_started PROPERTY CXX_STANDARD 17)
    set_property(TARGET getting_started PROPERTY CXX_STANDARD_REQUIRED YES)
    set_property(TARGET getting_started PROPERTY CXX_EXTENSIONS NO)

then:

.. code-block:: console

    $ mkdir build
    $ cd build
    $ cmake ../
    $ make
    $ ./getting_started

-----------------------------------------------------------------------

Python
------

If you have successfully compiled and installed dcgpy following the :ref:`installationguide` you will be able to test its
use by running the following script:

.. _getting_started_py:

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a ``getting_started.py`` text file and run it with:

.. code-block:: console

   $ python getting_started.py

We recommend the use of `Jupyter <https://jupyter.org/>`_ or `Ipython <https://github.com/ipython/ipython>`_ to enjoy ``dcgpy`` the most.

.. _installationguide:

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* `Boost <http://www.boost.org/>`_, various C++ utilities. (>=1.72).
* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, linear algebra library. (>=3.3.0)
* `Pagmo <https://github.com/esa/pagmo2>`_, parallel optimization library. (>=2.15).
* `tbb <https://github.com/intel/tbb>`_, lets you easily write parallel C++ programs that take full advantage of multicore performance. (>=2020.1).
* `AuDi <http://darioizzo.github.io/audi/>`_, high order automated differentiation library. (>=1.8).
* `obake <https://github.com/bluescarni/obake>`,  symbolic manipulation of sparse polynomials. (>=0.6.0). 
* `Symengine <https://github.com/symengine/symengine>`_, symbolic manipulation of math expressions. (>=0.6).

In case you are familiar with the conda package manager an environent ready for the dcgp installation can be created via the single command

.. code-block:: console

   $ conda create -n build_dcgp cmake cxx-compiler boost-cpp pagmo-devel tbb-devel audi symengine obake-devel

After making sure the dependencies above are installed and found in your system, you may download
the latest dCGP code via git:

.. code-block:: console

   $ git clone https://github.com/darioizzo/dcgp.git

and configure your build using CMake. For example, if you are using conda, make a build directory and type:

.. code-block:: console

   $ cmake -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/build_dcgp/ -DCMAKE_PREFIX_PATH=~/miniconda3/envs/build_dcgp/ ../

where the conda environment has been assumed in ```~/miniconda3/envs/build_dcgp/```.

When done, type (in your build directory):

.. code-block:: console

   $ make 

When finished, to run the tests type:

.. code-block:: console

   $ make test

If successfull, you may now install cgp:

.. code-block:: console

   $ make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 
To check that all went well compile the :ref:`quick-start example <getting_started_c++>`.



-----------------------------------------------------------------------

Python
------
The main functionalities of dCGP are exposed into a Python module called ``dcgpy`` which
can be installed from conda (OSx, linux and Win), pip (only linux) or by building the module.

The following third party dependencies are required to have full access to the ``dcgpy`` API:

* `numpy <https://numpy.org/>`_, The fundamental package for scientific computing with Python. (>=1.18)
* `matplotlib <https://matplotlib.org/>`_,  A comprehensive library for creating static, animated, and interactive visualizations in Python. (>=3.2)
* `pyaudi <http://darioizzo.github.io/audi/>`_, A library that implements the differential algebra of Taylor truncated polynomials. (>=1.8)
* `sympy <https://www.sympy.org/en/index.html>`_, A Python library for symbolic mathematics. (>=1.6)
* `graphviz <https://graphviz.readthedocs.io/en/stable/>`_, A simple pure-Python interface for the Graphviz graph-drawing software. (>=2.42)


Installing with conda
^^^^^^^^^^^^^^^^^^^^^
``dcgpy`` is available in the `conda <https://conda.io/en/latest/>`__ package manager
from the `conda-forge <https://conda-forge.org/>`__ channel. A single package is available:

* `dcgp-python <https://anaconda.org/conda-forge/dcccgp-python>`__, which contains the ``dcgpy`` python module.

In order to install ``dcgpy`` via conda, you just need
to add ``conda-forge`` to the channels:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda install dcgp-python

note that all the required dependencies will be installed automatically, as well as the C++ ```dcgp``` headers.

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

You may test the successfull installation by running the python tests typing:

.. code-block:: console

   $ python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"


Building
^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the python module you need to first install the dcgp C++ header library and its dependencies (see above) as well as the additional dependency:

* `pybind11 <https://github.com/pybind/pybind11>`_, Seamless operability between C++11 and Python. (>=2.5.0). 

In case you are familiar with the conda package manager an environent ready for the python module installation can be created via the single command

.. code-block:: console

   $ conda create -n build_dcgp cmake cxx-compiler boost-cpp pagmo-devel tbb-devel audi symengine obake-devel pybind11

Install the latest dCGP code via git:

.. code-block:: console

   $ git clone https://github.com/darioizzo/dcgp.git

After installing the C++ dcgp library (see above) and making sure your environment is correctly set up to find all dependencies, you can 
configure your python module build using CMake. For example, if you are using conda, in the build directory type:

.. code-block:: console

   $ cmake -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/build_dcgp/ -DCMAKE_PREFIX_PATH=~/miniconda3/envs/build_dcgp/ -DDCGP_BUILD_DCGP=OFF -DDCGP_BUILD_DCGPY=ON ../

where the conda environment has been assumed in ```~/miniconda3/envs/build_dcgp/```.

When done, type (in your build directory):

.. code-block:: console

   $ make install

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.
Differentiable Cartesian Genetic Programming
============================================

Differentiable Cartesian Genetic Programming (dCGP) is a recent development in the field of Genetic Programming
that adds the information about the derivatives of the output nodes (the programs, or expressions encoded) with
respect to the input nodes (the input values) and/or weights. In doing so, it enables a number of new applications
currently the subject of active research.

The evolution of the genetic program can now be supported by using the information on the derivatives, hence enabling
for the equivalent of back-propagation in Neural Networks. The fitness function can be defined in terms of the
derivatives, allowing to go beyond simple regression tasks and, additionally, solve differential equations, learn
differential models, capture conserved quantities in dynamical systems.

Beyond the standard **CGP** we provide two new encodings called **dCGP-W** and **dCGPANN**. 
The first one adds weights to the standard CGP nodes connections, while the second one allows to encode and evolve 
Artificial Neural Networks, while learning its parameters too, using an underlying CGP representation.

.. figure:: _static/expression_home.png 
   :alt: dCGP expression
   :align: center
   :width: 800px

   A dCGP expression and its differential expansion in a specified point

----------------------------------------------------------------

Table of contents:
^^^^^^^^^^^^^^^^^^^^

.. toctree::
  :maxdepth: 1

  installation
  quickstart
  docs/index
  tutorials/index
  theory/index
  credits
  changelog

---------------------------------------------------------------

References
^^^^^^^^^^

.. [dCGP1] Dario Izzo, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." In European Conference on Genetic Programming, pp. 35-51. Springer, 2017. `Arxiv version (2016) <https://arxiv.org/pdf/1611.04766v1.pdf>`_ 

.. [dCGP2] Marcus Märtens and Dario Izzo. Neural network architecture search with differentiable cartesian genetic programming for regression. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (pp. 181-182). ACM. `Arxiv version (2018) <https://arxiv.org/abs/1907.01939>`_ 
.. _changelog:

Changelog
=========

1.6 (Unreleased)
------------------

* All UDAs, :class:`dcgpy.es4cgp`, :class:`dcgpy.mes4cgp`, 
  :class:`dcgpy.moes4cgp`, :class:`dcgpy.momes4cgp` use random mutations rather than active mutations
  to evolve expressions. This improves significantly the algorithm behaviour.
* Bug fix when using a large number of ephemeral constants ... now udas crash is avoided
* New :func:`dcgpy.enable_threading` and :func:`dcgpy.disable_threading` utilities 
  allowing to switch off multithreding entirely. 

1.5 (6/7/2020)
-------------------

New 
~~~

* A new UDA is introduced to solve symbolic regression problems. 
  Its called moes (Multi-Objective Evolutionary Startegy) and completes the 
  evolutionary approaches in dcgp which can now be selected to be memetic or not
  and single objective or multi-objective.

Changes
~~~~~~~

- BREAKING: the API has been made uniform for the four UDAs: :class:`dcgpy.es4cgp`, :class:`dcgpy.mes4cgp`, 
  :class:`dcgpy.moes4cgp`, :class:`dcgpy.momes4cgp` as well as the mutation mechanism. 
  Named parameters have thus changed and default values too. Note that, for example, what
  was *n_mut* in some algos, is now *max_mut*.

- The underlying computations of the symbolic regression optimization problem (UDP) 
  is now performed by obake using a vectorized type. Speed improvements are observed
  of magnitudes between x4 and x100.

- The problem on nans appearing and exceptions being thrown has been solved 
  for :class:`dcgpy.symbolic_regression` by guarding against symengine exceptions
  and by discarding zero columns and rows when inverting hessians for the Newton step of memetic algorithms.

- The UDA :class:`dcgpy.es4cgp` is no longer using a thread bfe to compute the loss. This avoids crashes when pythonic, 
  non thread-safe kernels are used. A bfe can still be set by the user (deprecated in python) after
  the UDA has been instantiated.
  
- Documentation has been improved and all tutorials and examples updated to the new API.

Differentiable Cartesian Genetic Programming (dCGP)
===================================================
Differentiable Cartesian Genetic Programming (dCGP) is a recent development of CGP
that adds the information about the derivatives of the output nodes (the programs,
or expressions encoded) with respect to the input nodes (the input values) and/or the
weights (that can be associated with every connection in the graph, similarly to Neural Networks).

.. figure:: ../_static/expression_theory.png
   :alt: weighted dCGP expression
   :align: center
   :width: 800px

   A weighted dCGP expression

The derivatives are obtained by means of Automatic Differentiation through the
`AuDi <http://darioizzo.github.io/audi/>`_ library, that implements the Taylor truncated
polynomial algebra. 

The use of the derivatives of the outputs (and hence of any fitness function that is a
combination of these) with respect to inputs/weights enables a number of new applications
currently the subject of active research.

The evolution of the genetic program can now be supported by using the information
on the derivatives, hence enabling for the equivalent of back-propagation in Neural Networks.

Furthermore, the fitness function itself can be defined in terms of the derivatives, allowing 
dCGP to encode solutions to initial value problems, boundary value problems, Lyapunov functions, etc.

.. figure:: ../_static/expression_ann.png
   :alt: weighted dCGP expression
   :align: center
   :width: 500px

   A dCGP-ANN expression

Possible applications of a dCGP are:

* Symbolic Regression.
* Solving differential equations.
* Searching for Lyapunov functions (in non-linear control problems).
* Searching for conserved quantities in dynamical systems.
* Concurrent learning of network topologies and parameters. 

In [dCGP1]_ we introduce the ideas behind dCGP and study a few interesting applications. we
present an extension to the normal, canonical CGP, which introduces weights on the node connections
and allows for a variable number of model parameters, albeit complicating the learning process.

In [dCGP2]_ we introduce a further new variant of a canonical CGP, the dCGP-ANN, which makes use
of backpropagation to learn the many model parameters (biases and weights) of an artificial neural network 
represented by the CGP encoding.
Background
====================

Cartesian Genetic Programming (CGP) is a tree based representation of a computer program introduced by 
Julian F. Miller and Peter Thomson in 1997. Several open source codes exist that implement some
form of CGP, Miller himself maintain a resource page `here <https://www.cartesiangp.com/resources>`_.

The dcgpy library here introduced, differs in several aspects from all existing implementations:
it allows to perform any-order derivatives on the represented programs (hence implementing a
differentiable form of genetic programming), it allows to use Python to use runtime scripting, 
it is thread-safe and it allows to define kernels as generic functors.

.. toctree::
  :maxdepth: 1

  cgp
  dcgp
Cartesian Genetic Programming (CGP)
====================================

.. figure:: ../_static/cgp.jpg
   :alt: The original CGP encoding
   :align: center
   :width: 800px

Cartesian Genetic Programming (CGP) is a form of Genetic Programming where a computer program
is represented by an acyclic graph and is evolved using a standard evolutionary strategy.

For more information visit `Cartesian Genetic Programming <http://www.cartesiangp.co.uk/>`_.
Tutorials
====================
.. toctree::
  :maxdepth: 1

  cpp/index
  py/index
Python
=================

Core
^^^^^^^

.. toctree::
  :maxdepth: 1

  ../../notebooks/custom_kernels.ipynb


Symbolic Regression 
^^^^^^^^^^^^^^^^^^^
.. toctree::
  :maxdepth: 1

  ../../notebooks/symbolic_regression_1.ipynb
  ../../notebooks/symbolic_regression_2.ipynb
  ../../notebooks/symbolic_regression_3.ipynb
  ../../notebooks/real_world1.ipynb
  ../../notebooks/real_world2.ipynb
  ../../notebooks/real_world3.ipynb
  ../../notebooks/evo_in_parallel.ipynb

Phenotype correction
^^^^^^^^^^^^^^^^^^^^

.. toctree::
  :maxdepth: 1

  ../../notebooks/phenotype_correction_ex.ipynb


Neural Networks
^^^^^^^^^^^^^^^
 
.. toctree::
  :maxdepth: 1

  ../../notebooks/An_intro_to_dCGPANNs.ipynb
  ../../notebooks/dCGPANNs_for_classification.ipynb
  ../../notebooks/dCGPANNs_for_function_approximation.ipynb

Other
^^^^^^^^

.. toctree::
  :maxdepth: 1

  ../../notebooks/solving_odes.ipynb
  ../../notebooks/weighted_symbolic_regression.ipynb
  ../../notebooks/finding_prime_integrals.ipynb


Deprecated
^^^^^^^^^^

.. toctree::
  :maxdepth: 1

  ../../notebooks/learning_constants.ipynb
  ../../notebooks/learning_constants2.ipynb

Find a model including one parameter (using a purely evolutionary approach)
===========================================================================

In this second tutorial we show how to find a model for our input data when we also
want to learn some constants. 

Constants can, in general, be learned via two main techniques:
 1. evolutionary (common and standard practice in GP)
 2. memetic (original with dCGP)

The difference is that the evolutionary approach cannot find the correct and exact values for constants, only approximations. 
In this tutorial we follow the evolutionary approach 1. In the next tutorial we will follow a memetic approach 2.
We use the problem P1 from the dcgp::gym, that is x^5 - pi*x^3 + x

Code:
^^^^^^^^
.. literalinclude:: ../../../../examples/symbolic_regression_2.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output will be different on your computers as its non deterministic.

.. code-block:: python

   Gen:        Fevals:          Best:	Constants:	Formula:
      0              0        3803.05	[1.11378]	[x0 + c1**2] ...
    500           2000         8.2487	[0.959988]	[(-c1 + x0)*x0**4 - x0**2] ...
   1000           4000        2.70799	[1.14114]	[-x0**2*c1 + (-c1 + x0)*x0**4*c1] ...
   1500           6000        2.70797	[1.14053]	[-x0**2*c1 + (-c1 + x0)*x0**4*c1] ...
   2000           8000        2.70796	[1.14064]	[-x0**2*c1 + (-c1 + x0)*x0**4*c1] ...
   2500          10000      0.0901456	[1.0711]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   3000          12000      0.0597183	[1.07233]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   3500          14000      0.0597183	[1.07233]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   4000          16000      0.0597183	[1.07233]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   4500          18000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   5000          20000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   5500          22000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   6000          24000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   6500          26000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   7000          28000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   7500          30000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   8000          32000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   8500          34000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   9000          36000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
   9500          38000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
  10000          40000      0.0596577	[1.0723]	[-x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3] ...
  Exit condition -- generations = 10000

  Best fitness: [0.0596577]
  Chromosome: [1.0723, 1, 0, 0, 2, ... ]
  Pretty Formula: [((((x0*c1)*((x0*c1)*x0))*((x0*c1)*(x0-c1)))-((x0*c1)*((x0*c1)*x0)))]
  Prettier Formula: -x0**3*c1**2 + (-c1 + x0)*x0**4*c1**3
  Expanded Formula: -x0**3*c1**2 - x0**4*c1**4 + x0**5*c1**3
  
Solving the NLODE3 problem from Tsoulos
================================================

.. literalinclude:: ../../../../examples/tsoulos_nlode3.cpp
   :language: c++
   :linenos:Search for first integrals of Kepler's problem
================================================

.. literalinclude:: ../../../../examples/hamiltonian_kepler.cpp
   :language: c++
   :linenos:Solving the ODE1 problem from Tsoulos
================================================

.. literalinclude:: ../../../../examples/tsoulos_ode1.cpp
   :language: c++
   :linenos:Find an exact model for the Koza quintic problem
==================================================

In this first tutorial we show how to find an exact formula for some input data that do not require
any real valued constant. This is the easiest case for a symbolic regression task and thus makes it for a perfect entry tutorial.

We use the classic problem Koza quintic polynomial, that is x - 2x^3 + x^5.

Code:
^^^^^^^^
.. literalinclude:: ../../../../examples/symbolic_regression_1.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output will be different on your computers as its non deterministic.

.. code-block:: python

    Gen:        Fevals:          Best:  Constants:    Formula:
       0              0        3898.35           []    [2*x0**3] ...
     500           2000        638.426           []    [x0**5] ...
    1000           4000        138.482           []    [(-x0**2 + x0**4)*x0] ...
    1500           6000        101.734           []    [-x0 + (-x0**2 + x0**4)*x0] ...
    1698           6792     5.2071e-30           []    [x0*(1 - x0**2) - x0**3*(1 - x0**2)] ...
    Exit condition -- ftol < 1e-08
    
    Best fitness: [5.2071e-30]
    Chromosome: [2, 0, 0, 3, 1, ... ]
    Pretty Formula: [(((((x0*x0)/(x0*x0))-(x0*x0))*x0)-(((((x0*x0)/(x0*x0))-(x0*x0))*x0)*(x0*x0)))]
    Prettier Formula: x0*(1 - x0**2) - x0**3*(1 - x0**2)
    Expanded Formula: x0 - 2*x0**3 + x0**5Finding the prime integral for a spring mass system using Lipson method
=======================================================================

.. literalinclude:: ../../../../examples/hamiltonian_spring_mass_lipson.cpp
   :language: c++
   :linenos:Find an exact model inclduing one parameter using a memetic approach
====================================================================

In this second tutorial we show how to find a model for our input data when we also
want to learn some constants. 

Constants can, in general, be learned via two main techniques:
 1. evolutionary (common and standard practice in GP)
 2. memetic (original with dCGP)

The difference is that the evolutionary approach cannot find the correct and exact values for constants, only approximations. 
In this tutorial we follow the evolutionary approach 2. In the next tutorial we will follow a memetic approach 2.
We use the problem P1 from the dcgp::gym, that is x^5 - pi*x^3 + x

Code:
^^^^^^^^
.. literalinclude:: ../../../../examples/symbolic_regression_3.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output is non deterministic. Sometimes, with a bit of luck :)
the problem is solved exaclty (loss goes to zero). The following output reports one of these occurences.
Note that in the traditional evolutionary approach this result is incredibly hard to obtain.

.. code-block:: python

   Gen:        Fevals:          Best:	Constants:	Formula:
      0              0        4009.59	[-2.41031]	[(c1 + x0)*x0] ...
     50            200       0.978909	[0.238557]	[x0**2*c1*(-x0 + x0**4) - (c1 + x0 + x0* ...
    100            400        0.84565	[0.240548]	[x0**2*c1*(-x0 + x0**4) - (c1 + 2*x0)] ...
    150            600       0.761757	[0.240032]	[-2*x0 + x0**2*c1*(-x0 + x0**4)] ...
    200            800      0.0170582	[-1.16484]	[(-x0 + x0**3)*(c1 + x0**2) - x0**3] ...
    250           1000      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    300           1200      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    350           1400      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    357           1428    2.17578e-29	[-1.14159]	[x0**3*(c1 + x0**2) - (-x0 + 2*x0**3)] ...
   Exit condition -- ftol < 1e-08
   
   Best fitness: [2.17578e-29]
   Chromosome: [-1.14159, 2, 0, 0, 2, ... ]
   Pretty Formula: [(((c1+(x0*x0))*((x0*x0)*x0))-(((x0*x0)*x0)+(((x0*x0)*x0)-x0)))]
   Prettier Formula: x0**3*(c1 + x0**2) - (-x0 + 2*x0**3)
   Expanded Formula: x0 + x0**3*c1 - 2*x0**3 + x0**5
C++ 
=================

Symbolic Regression 
^^^^^^^^^^^^^^^^^^^
.. toctree::
  :maxdepth: 1

  cpp_symbolic_regression_1
  cpp_symbolic_regression_2
  cpp_symbolic_regression_3
  cpp_symbolic_regression_4

Others
^^^^^^^^^^^^^^^^^^^
.. toctree::
  :maxdepth: 1

  cpp_tsoulos_ode1
  cpp_tsoulos_nlode3
  cpp_hamiltonian_spring_mass
  cpp_hamiltonian_spring_mass_lipson
  cpp_kepler_first_integralFinding an entire non dominated front of formulas.
==================================================

In this fourth tutorial on symbolic regression we solve the multiobjective symbolic regression problem.
The Mean Squared Error (i.e the loss) of our model is considered next to the model complexity to determine how good 
a certain model is. The result is thus a whole non-dominated front of models.

This case is arguably the most complete and useful among  symbolic regression tasks. We use here
the problem vladi6 from the dcgp::gym, that is: 6*cos(x*sin(y))

Code:
^^^^^^^^
.. literalinclude:: ../../../../examples/symbolic_regression_4.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output will be different on your computers as its non deterministic.

.. code-block:: python

  Gen:        Fevals:     Best loss: Ndf size:  Compl.:
      0              0        4.10239        11        62
     10           1000        1.42703         7        74
     20           2000       0.828554         8        45
     30           3000       0.374803        13        78
     40           4000       0.164032        16        66
     50           5000        0.03552        14        48
     60           6000      0.0200792        11        45
     70           7000              0        10        19
     80           8000              0        10        19
     90           9000              0        10        19
    100          10000              0        11        19
    110          11000              0        10        18
    120          12000              0        11        18
    130          13000              0        11        18
    140          14000              0        11        18
    150          15000              0        11        18
    160          16000              0        11        18
    170          17000              0        11        18
    180          18000              0        11        18
    190          19000              0        11        18
    200          20000              0        11        18
    210          21000              0        11        18
    220          22000              0        11        18
    230          23000              0        11        18
    240          24000              0        10        18
    250          25000              0        10        18
    Exit condition -- max generations = 250
    
    Non dominated Front at the end:
    1  - Loss: 0            Complexity:    18   Formula:  c1*cos(x0*sin(x1))
    2  - Loss: 0.844272     Complexity:    17   Formula:  3*c1 - 2*x0 - 2*x0*x1
    3  - Loss: 1.10197      Complexity:    15   Formula:  4*c1 - 3*x0 - x0*x1
    4  - Loss: 1.17331      Complexity:    14   Formula:  3*c1 - 3*x0 - 2*x1
    5  - Loss: 1.27379      Complexity:    10   Formula:  c1 - 3*x0 - x1
    6  - Loss: 1.92403      Complexity:    9    Formula:  2*c1 - 3*x0
    7  - Loss: 1.92403      Complexity:    7    Formula:  c1 - 3*x0
    8  - Loss: 3.22752      Complexity:    5    Formula:  c1 - x0
    9  - Loss: 4.74875      Complexity:    2    Formula:  c1
    10 - Loss: 4.8741       Complexity:    1    Formula:  4
Finding the prime integral for a spring mass system
=====================================================

.. literalinclude:: ../../../../examples/hamiltonian_spring_mass.cpp
   :language: c++
   :linenos:.. docs tab

API documentation
=================

The documentation is split into two parts. The C++ documentation and the Python documentation. 
We try to keep the APIs as similar as possible, but since templates are not available in python, 
we expose an arbitrary number of template instantiations indicating in the class name
after an underscore, the type. For example, the Python class *expression_double* corresponds, 
to the C++ class *expression<double>*.

-------------------------------------------------

.. toctree::
  :maxdepth: 1

  cpp/index
  python/index
 Multi-Objective Memetic Evolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. figure:: ../../_static/multiple_objectives.png
   :alt: MES

.. autoclass:: dcgpy.momes4cgp  
   :members:Kernel set
-------------

kernel_set_double
^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_double

    .. automethod:: dcgpy.kernel_set_double.push_back()

kernel_set_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_gdual_double

    .. automethod:: dcgpy.kernel_set_gdual_double.push_back()

kernel_set_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_gdual_vdouble

    .. automethod:: dcgpy.kernel_set_gdual_vdouble.push_back()
NIST problems (StRD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
These problems contain real observation of physical processes
collected in the NIST web pages (https://www.nist.gov/itl) 

.. autofunction:: dcgpy.generate_chwirut1
.. autofunction:: dcgpy.generate_chwirut2
.. autofunction:: dcgpy.generate_daniel_wood
.. autofunction:: dcgpy.generate_gauss1
.. autofunction:: dcgpy.generate_kirby2
.. autofunction:: dcgpy.generate_lanczos2
.. autofunction:: dcgpy.generate_misra1b
expression_weighted (dCGP-W)
----------------------------

.. figure:: ../../_static/expression_weighted.png
   :alt: weighted dCGP expression
   :align: center

Three types of expression_weighted are currently available in python. 
They differ from the type they operate upon: double, gdual, vectorized gdual.

--------------------------------------------------------

expression_weighted_double
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_double
    :members:

--------------------------------------------------------

expression_weighted_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_gdual_double
    :members:

--------------------------------------------------------

expression_weighted_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_gdual_vdouble
    :members:Problems P1-P7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
These seven problems are introduced studied in the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

.. autofunction:: dcgpy.generate_P1
.. autofunction:: dcgpy.generate_P2
.. autofunction:: dcgpy.generate_P3
.. autofunction:: dcgpy.generate_P4
.. autofunction:: dcgpy.generate_P5
.. autofunction:: dcgpy.generate_P6
.. autofunction:: dcgpy.generate_P7Kernel
-------------

kernel_double
^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_double
    :members:

kernel_gdual_double
^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_gdual_double
    :members:

kernel_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_gdual_vdouble
    :members:
Memetic Evolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. figure:: ../../_static/neo-darwinism.jpg
   :alt: MES

.. autoclass:: dcgpy.mes4cgp  
   :members:Problems from Vladislavleva
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These problems are introduced and studied in the paper:

Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

.. autofunction:: dcgpy.generate_kotanchek
.. autofunction:: dcgpy.generate_salutowicz
.. autofunction:: dcgpy.generate_salutowicz2d
.. autofunction:: dcgpy.generate_uball5d
.. autofunction:: dcgpy.generate_ratpol3d
.. autofunction:: dcgpy.generate_sinecosine
.. autofunction:: dcgpy.generate_ripple
.. autofunction:: dcgpy.generate_ratpol2dEvolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. figure:: ../../_static/EvolutionaryStrategy.jpg
   :alt: ES

.. autoclass:: dcgpy.es4cgp  
   :members:.. python docs

Python Documentation
====================

Utilities
---------

Assorted functionalities.

.. toctree::
   :maxdepth: 1

   threading


Kernels
---------

.. figure:: ../../_static/kernels.png
   :alt: dCGP kernels
   :align: right
   :scale: 40 %

**Kernels**, (also called **non-linearities** in the ANN literature) describe the fundamental computational units of a CGP. 
Things like addition, multiplication, trigonomoetric functions are all kernels. The classes :class:`dcgpy.kernel_double`, 
:class:`dcgpy.kernel_gdual_double` and :class:`dcgpy.kernel_gdual_vdouble` allow the user
the definition of their own kernels able to operate on the choosen type. 

The most popular kernels are already coded and shipped with dcgpy. A python list containing multiple kernels can be easily instantiated 
via the classes :class:`dcgpy.kernel_set_double`, :class:`dcgpy.kernel_set_gdual_double` and :class:`dcgpy.kernel_set_gdual_vdouble`.

.. toctree::
  :maxdepth: 1

  kernel
  kernel_set
  kernel_list

----------------------------------------------------------------------------------

Types of dCGPs
--------------

Several types of **Cartesian Genetic Program** are provided in dcgpy. Since a dCGP is some kind of a mathematical expression,
we use the term expression to name the related, different, classes. 

We essentially provide three types of CGPs:

* **expression**: this is the original CGP as introduced by Miller in 1999.
* **expression_weighted**: this adds to the original CGP formulation weights on each of the graph edges. (original with dCGP - 2016)
* **expression_ann**: this represents an Artificial Neural Network inclusive of biases and weights, via a CGP and allows to learn the model parameters using backproagation. (original with dCGP - 2018)

Each of the above CGPs can operate over different numerical types. For example an **expression** can operate over floats (:class:`dcgpy.expression_double`), in which case the 
result of evaluating the inner computational graph will be a float, but also on gduals (:class:`dcgpy.expression_gdual_double`), in which case, the result of evaluating the
inner computational graph will be a gdual (hence it will contain all the program derivatives up to the chosen order and hence is referred to as a dCGP.)

Another important type some CGPs can operate upon is the vectorized gdual (:class:`dcgpy.expression_gdual_vdoubles`). This type is the same as the gdual type, but its vectorized,
allowing order of magnitude speed ups when a CGP needs to be evaluated over several points (such as in the case of a loss evaluation)

.. toctree::
  :maxdepth: 1

  expression
  expression_weighted
  expression_ann

----------------------------------------------------------------------------------

 
Symbolic Regression
---------------------

.. figure:: ../../_static/non_linear_regression.png
   :alt: non-linear regression
   :align: right
   :scale: 70 %

Mathematically, a symbolic regression problem is a global optimization problem. In order to facilitate its solution,
a number of classes have been designed to interface a dcgpy ``expression`` to the optimisation suite pygmo. 
In particular we provide UDPs (in pagmo's jargon user defined problems) that can be used to build :class:`pygmo.problem` objects 
(representing symbolic regression problems) and UDAs (in pagmo's jargon user defined algorithms) that can be used
to build :class:`pygmo.algorithm` objects.

.. toctree::
  :maxdepth: 1

  symbolic_regression
  es4cgp
  gd4cgp
  mes4cgp
  momes4cgp


We also make available, as a gym to test the capabilities of various proposed methodologies, a number of data sets
coming from different scientific sources. These are collected in what we call the ``dcgpy gym``.

.. toctree::
  :maxdepth: 1

  koza_quintic
  p_problems
  vladislavleva_problems
  nist_problems

Available kernels
----------------------------------

When constructing a :class:``dcgpy.kernel_set_double``, :class:``dcgpy.kernel_set_gdual_double `` 
or :class:``dcgpy.kernel_set_gdual_vdouble`` we can use the following names to add the corresponding
kernels to the set.

---------------------------------------------------------------------------

.. cssclass:: table-bordered table-striped

   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |Kernel name     |  Function              |  Definition                                                                           |
   +================+========================+=======================================================================================+
   | **Basic operations**                                                                                                            |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sum"           |addition                | :math:`\sum_i x_i`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"diff"          |subtraction             | :math:`x_1 - \sum_{i=2} x_i`                                                          |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"mul"           |multiplication          | :math:`\prod_i x_i`                                                                   |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"div"           |division                | :math:`x_1 / \prod_{i=2} x_i`                                                         |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"pdiv"          |protected division      | :math:`x_1 / \prod_{i=2} x_i` or 1 if NaN                                             |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   | **Unary non linearities (ignoring all inputs after the first one)**                                                             |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sin"           |sine                    | :math:`\sin x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"cos"           |cosine                  | :math:`\cos x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"log"           |natural logarithm       | :math:`\log x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"exp"           |exponential             | :math:`e^{x_1}`                                                                       |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"gaussian"      |gaussian                | :math:`e^{-x_1^2}`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sqrt"          |square root             | :math:`\sqrt{x_1}`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"psqrt"         |protected square root   | :math:`\sqrt{|x_1|}`                                                                  |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   | **Non linearities suitable also for dCGPANN**                                                                                   |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sig"           |sigmoid                 | :math:`\frac{1}{1 + e^{-\sum_i x_i}}`                                                 |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"tanh"          |hyperbolic tangent      | :math:`\tanh \left(\sum_i x_i\right)`                                                 |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ReLu"          |rectified linear unit   | :math:`\sum_i x_i` if positive, 0 otherwise                                           |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ELU"           |exp linear unit         | :math:`\sum_i x_i` if positive, :math:`e^{\sum_i x_i} - 1` otherwise                  |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ISRU"          |Inverse square root unit| :math:`\frac{\sum_i x_i}{\sqrt{1 + \left(\sum_i x_i\right)^2}}`                       |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sum"           |addition                | :math:`\sum_i x_i`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
expression_ann (dCGP-ANN)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. important::
   This Cartesian Genetic Program is able to encode an Artificial Neural Network. Weights and biases are added to the acyclic graph
   as well as extra methods to allow to perform backpropagation (in parallel), to visualize the network and more ...

.. autoclass:: dcgpy.expression_ann
    :members:
Gradient Descent (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. figure:: ../../_static/gd.png
   :alt: gradient descent

.. autoclass:: dcgpy.gd4cgp  
   :members:Symbolic Regression (UDP)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../../_static/symbolic_regression.jpg
   :alt: weighted dCGP expression

.. autoclass:: dcgpy.symbolic_regressionKoza Quintic
^^^^^^^^^^^^^

.. autofunction:: dcgpy.generate_koza_quinticThreading
---------

disable_threading
^^^^^^^^^^^^^^^^^

.. autofunction:: dcgpy.disable_threading

enable_threading
^^^^^^^^^^^^^^^^^

.. autofunction:: dcgpy.enable_threading

expression (dCGP)
-----------------

.. figure:: ../../_static/expression.png
   :alt: dCGP expression
   :align: center

Three types of expression are currently available in python. 
They differ from the type they operate upon: double, gdual, vectorized gdual.

.. contents::

---------------------------------------------------------------

expression_double
^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_double
    :members:

---------------------------------------------------------------

expression_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_gdual_double
    :members:

---------------------------------------------------------------

expression_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_gdual_vdouble
    :members:Multi-Objective Memetic Evolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: c++
 
.. doxygenclass:: dcgp::momes4cgp
   :project: dCGP
   :members:kernel_set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Classes such as :cpp:class:`dcgp::expression`, :cpp:class:`dcgp::symbolic_regression` and others are constructed 
from an ``std::vector`` of :cpp:class:`dcgp::kernel`. Assembling such a ``std::vector`` is an operation that is greatly 
facilitated by the class :cpp:class:`dcgp::kernel_set`.

Intended use of the class is:

.. highlight:: c++

.. code-block:: c++

   // (optional) - Say we have a custom kernel named f 
   kernel<double> f(my_sum<double>, print_my_sum, "my_sum");
   // We can construct a kernel set (and fill it in with provided kernels)
   kernel_set<double> kernels({"sum", "div", "mul","diff"});
   // We can get a vector with these four basic kernels ...
   auto kernel_vector = kernels();
   // (optional) - ... or add our own custom kernel to the set ...
   basic_set.push_back(f);
   // (optional) - ... and get a vector with all five kernels
   auto kernel_vector2 = kernels();

---------------------------------------------------------------------------

.. doxygenclass:: dcgp::kernel_set
   :project: dCGP
   :members:
expression_weighted (dCGP-W)
--------------------------------

This class represents a **Weighted Cartesian Genetic Program**. Each node connection is associated to a weight so that more generic mathematical expressions
can be represented. When instantiated with the type *gdual<T>*, also the weights are defined as gduals, hence the program output can be expanded also with respect to the weights
thus allowing to train the weights using algorithms such as stochastic gradient descent, while the rest of the expression remains fixed. 


The class template can be instantiated using the types *double* or *gdual<T>*. 

.. figure:: ../../_static/expression_weighted.png
   :alt: weighted dCGP expression
   :align: center

   A weighted dCGP expression

.. doxygenclass:: dcgp::expression_weighted
   :project: dCGP
   :members:kernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. highlight:: c++
 
.. doxygenclass:: dcgp::kernel
   :project: dCGP
   :members:Memetic Evolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: c++
 
.. doxygenclass:: dcgp::mes4cgp
   :project: dCGP
   :members:Evolutionary Strategy (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: c++
 
.. doxygenclass:: dcgp::es4cgp
   :project: dCGP
   :members:.. cpp docs

C++ Documentation
=================

.. contents::

Kernels
--------- 

.. figure:: ../../_static/kernels.png
   :alt: dCGP kernels
   :align: right
   :scale: 40 %

**Kernels**, (also called **non-linearities** in the ANN literature) describe the fundamental computational units of a CGP. 
Things like addition, multiplication, trigonomoetric functions are all kernels. The templated class :cpp:class:`dcgp::kernel` allow 
the user the definition of their own kernels able to operate on the choosen type. 

The most popular kernels are already coded and shipped with dcgpy. A python list containing multiple kernels can be easily 
instantiated via the :cpp:class:`dcgp::kernel_set`.

.. toctree::
  :maxdepth: 1

  kernel
  kernel_set
  kernel_list

----------------------------------------------------------------------------------

Types of dCGPs
--------------

Several types of **Cartesian Genetic Program** are provided in dcgpy. Since a dCGP is some kind of a mathematical expression,
we use the term expression to name the related, different, classes. 

We essentially provide three types of CGPs:

* **expression**: this is the original CGP as introduced by Miller in 1999.
* **expression_weighted**: this adds to the original CGP formulation weights on each of the graph edges. (original with dCGP - 2016)
* **expression_ann**: this represents an Artificial Neural Network inclusive of biases and weights, via a CGP and allows to learn the model parameters using backproagation. (original with dCGP - 2018)

Each of the above CGPs can operate over different numerical types, hence the corresponding classes are templated. 
For example a :cpp:class:`dcgp::expression <dcgp::expression>` can operate over doubles (``T`` = ``double``), in which case the 
result of evaluating the inner computational graph is a ``double``, but also on a ``gdual`` (``T`` = :cpp:class:`audi::gdual <gdual>` ``<Cf>``)
with ``Cf`` = ``double``, in which case, the result of evaluating the inner computational graph will
be a :cpp:class:`audi::gdual <gdual>` ``<Cf>`` and hence it will contain all the program derivatives
up to the chosen order and is thus referred to as a dCGP. 

Another important type some CGPs can operate upon is the :cpp:class:`audi::gdual <gdual>` ``<Cf>`` with ``Cf`` = ``vectorized_double``.
This type offers order of magnitude speed ups whenever a CGP needs derivatives and to be evaluated over several points 
(such as in the case of a loss evaluation).

.. toctree::
  :maxdepth: 1

  expression
  expression_weighted
  expression_ann

----------------------------------------------------------------------------------

 
Symbolic Regression
---------------------

.. figure:: ../../_static/non_linear_regression.png
   :alt: non-linear regression
   :align: right
   :scale: 70 %

Mathematically, a symbolic regression problem is a global optimization problem. In order to facilitate its solution,
a number of classes have been designed to interface a :cpp:class:`dcgp::expression` to the optimisation suite pygmo. 
In particular we provide UDPs (in pagmo's jargon user defined problems) that can be used to build
:cpp:class:`pagmo::problem <pagmo::problem>` objects (representing symbolic regression problems) and UDAs 
(in pagmo's jargon user defined algorithms) that can be used to build :cpp:class:`pagmo::algorithm <pagmo::algorithm>` objects.

.. toctree::
  :maxdepth: 1

  symbolic_regression
  es4cgp
  gd4cgp
  mes4cgp
  momes4cgp




  Available kernels
----------------------------------

When constructing a :cpp:class:`dcgp::kernel_set` we can use the following names to add the corresponding
kernels to the set.

---------------------------------------------------------------------------

.. cssclass:: table-bordered table-striped

   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |Kernel name     |  Function              |  Definition                                                                           |
   +================+========================+=======================================================================================+
   | **Basic operations**                                                                                                            |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sum"           |addition                | :math:`\sum_i x_i`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"diff"          |subtraction             | :math:`x_1 - \sum_{i=2} x_i`                                                          |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"mul"           |multiplication          | :math:`\prod_i x_i`                                                                   |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"div"           |division                | :math:`x_1 / \prod_{i=2} x_i`                                                         |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"pdiv"          |protected division      | :math:`x_1 / \prod_{i=2} x_i` or 1 if NaN                                             |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   | **Unary non linearities (ignoring all inputs after the first one)**                                                             |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sin"           |sine                    | :math:`\sin x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"cos"           |cosine                  | :math:`\cos x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"log"           |natural logarithm       | :math:`\log x_1`                                                                      |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"exp"           |exponential             | :math:`e^{x_1}`                                                                       |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"gaussian"      |gaussian                | :math:`e^{-x_1^2}`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sqrt"          |square root             | :math:`\sqrt{x_1}`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"psqrt"         |protected square root   | :math:`\sqrt{|x_1|}`                                                                  |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   | **Non linearities suitable also for dCGPANN**                                                                                   |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sig"           |sigmoid                 | :math:`\frac{1}{1 + e^{-\sum_i x_i}}`                                                 |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"tanh"          |hyperbolic tangent      | :math:`\tanh \left(\sum_i x_i\right)`                                                 |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ReLu"          |rectified linear unit   | :math:`\sum_i x_i` if positive, 0 otherwise                                           |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ELU"           |exp linear unit         | :math:`\sum_i x_i` if positive, :math:`e^{\sum_i x_i} - 1` otherwise                  |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"ISRU"          |Inverse square root unit| :math:`\frac{\sum_i x_i}{\sqrt{1 + \left(\sum_i x_i\right)^2}}`                       |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sum"           |addition                | :math:`\sum_i x_i`                                                                    |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"sin_nu"        |sine (non unary)        | :math:`\sin(\sum_i x_i)`                                                              |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"cos_nu"        |cosine (non unary)      | :math:`\cos(\sum_i x_i)`                                                              |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"gaussian_nu"   |gaussian (non unary)    | :math:`e^{\sum_i x_i}`                                                                |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"abs"           |absolute value          | :math:`\vert\sum_i x_i\vert`                                                          |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"inv_sum"       |inverse sum             | :math:`- \sum_i x_i`                                                                  |
   +----------------+------------------------+---------------------------------------------------------------------------------------+
   |"step"          |step function           | 1 if positive, 0 otherwise                                                            |
   +----------------+------------------------+---------------------------------------------------------------------------------------+expression_ann (dCGP-ANN)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class represents a **Artificial Neural Network Cartesian Genetic Program**. Each node connection is associated to a weight and each node to a bias. Only a subset of the kernel functions
is allowed, including the most used nonlinearities in ANN research: *tanh*, *sig*, *ReLu*, *ELU* and *ISRU*. The resulting expression can represent any feed forward neural network but also other
less obvious architectures. Weights and biases of the expression can be trained using the efficient backpropagation algorithm (gduals are not allowed for this class, they correspond to forward mode
automated differentiation which is super inefficient for deep networks ML.)

.. figure:: ../../_static/expression_ann.png
   :alt: weighted dCGP expression
   :align: center

   A, small, artificial neural network as using the dCPP-ANN approach.

.. doxygenclass:: dcgp::expression_ann
   :project: dCGP
   :members:
Gradient Descent (UDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: c++
 
.. doxygenclass:: dcgp::gd4cgp
   :project: dCGP
   :members:Symbolic Regression (UDP)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: c++
 
.. doxygenclass:: dcgp::symbolic_regression
   :project: dCGP
   :members:expression (dCGP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class represents a **Cartesian Genetic Program**. Since that is, essentially, an artificial genetic encoding for a mathematical expression, we named the templated class *expression*.
The class template can be instantiated using the types *double* or *gdual<T>*. In the case of *double*, the class would basically reproduce a canonical CGP expression. In the case of *gdual<T>*
the class would operate in the differential algebra of truncated Taylor polynomials with coefficients in *T*, and thus provide also any order derivative information on the program 
(i.e. the Taylor expansion of the program output with respect to its inputs).

.. figure:: ../../_static/expression.png
   :alt: dCGP expression
   :align: center

   A dCGP expression

.. doxygenclass:: dcgp::expression
   :project: dCGP
   :members: