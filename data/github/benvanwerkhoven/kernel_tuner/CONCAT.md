# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added
- new optimization strategies: dual annealing, greedly ILS, ordered greedy MLS, greedy MLS
- support for constant memory in cupy backend

### Removed
- Alternative Bayesian Optimization strategies that could not be used directly
- C++ wrapper module that was too specific and hardly used

## [0.4.1] - 2021-09-10
### Added
- support for PyTorch Tensors as input data type for kernels
- support for smem_args in run_kernel
- support for (lambda) function and string for dynamic shared memory size
- a new Bayesian Optimization strategy

### Changed
- optionally store the kernel_string with store_results
- improved reporting of skipped configurations

## [0.4.0] - 2021-04-09
### Added
- support for (lambda) function instead of list of strings for restrictions
- support for (lambda) function instead of list for specifying grid divisors
- support for (lambda) function instead of tuple for specifying problem_size
- function to store the top tuning results
- function to create header file with device targets from stored results
- support for using tuning results in PythonKernel
- option to control measurements using observers
- support for NVML tunable parameters
- option to simulate auto-tuning searches from existing cache files
- Cupy backend to support C++ templated CUDA kernels
- support for templated CUDA kernels using PyCUDA backend
- documentation on tunable parameter vocabulary

## [0.3.2] - 2020-11-04
### Added
- support loop unrolling using params that start with loop_unroll_factor
- always insert "define kernel_tuner 1" to allow preprocessor ifdef kernel_tuner
- support for user-defined metrics
- support for choosing the optimization starting point x0 for most strategies

### Changed
- more compact output is printed to the terminal
- sequential runner runs first kernel in the parameter space to warm up device
- updated tutorials to demonstrate use of user-defined metrics

## [0.3.1] - 2020-06-11
### Added
- kernelbuilder functionality for including kernels in Python applications
- smem_args option for dynamically allocated shared memory in CUDA kernels

### Changed
- bugfix for Nvidia devices without internal current sensor

## [0.3.0] - 2019-12-20
### Changed
- fix for output checking, custom verify functions are called just once
- benchmarking now returns multiple results not only time
- more sophisticated implementation of genetic algorithm strategy
- how the "method" option is passed, now use strategy_options

### Added
- Bayesian Optimizaton strategy, use strategy="bayes_opt"
- support for kernels that use texture memory in CUDA
- support for measuring energy consumption of CUDA kernels
- option to set strategy_options to pass strategy specific options
- option to cache and restart from tuned kernel configurations cachefile

### Removed
- Python 2 support, it may still work but we no longer test for Python 2
- Noodles parallel runner

## [0.2.0] - 2018-11-16
### Changed
- no longer replacing kernel names with instance strings during tuning
- bugfix in tempfile creation that lead to too many open files error

### Added
- A minimal Fortran example and basic Fortran support
- Particle Swarm Optimization strategy, use strategy="pso" 
- Simulated Annealing strategy, use strategy="simulated_annealing" 
- Firefly Algorithm strategy, use strategy="firefly_algorithm" 
- Genetic Algorithm strategy, use strategy="genetic_algorithm" 

## [0.1.9] - 2018-04-18
### Changed
- bugfix for C backend for byte array arguments
- argument type mismatches throw warning instead of exception

### Added
- wrapper functionality to wrap C++ functions
- citation file and zenodo doi generation for releases

## [0.1.8] - 2017-11-23
### Changed
- bugfix for when using iterations smaller than 3
- the install procedure now uses extras, e.g. [cuda,opencl]
- option quiet makes tune_kernel completely quiet
- extensive updates to documentation

### Added
- type checking for kernel arguments and answers lists
- checks for reserved keywords in tunable paramters
- checks for whether thread block dimensions are specified
- printing units for measured time with CUDA and OpenCL
- option to print all measured execution times

## [0.1.7] - 2017-10-11
### Changed
- bugfix install when scipy not present
- bugfix for GPU cleanup when using Noodles runner
- reworked the way strings are handled internally

### Added
- option to set compiler name, when using C backend

## [0.1.6] - 2017-08-17
### Changed
- actively freeing GPU memory after tuning
- bugfix for 3D grids when using OpenCL

### Added
- support for dynamic parallelism when using PyCUDA
- option to use differential evolution optimization
- global optimization strategies basinhopping, minimize

## [0.1.5] - 2017-07-21
### Changed
- option to pass a fraction to the sample runner
- fixed a bug in memset for OpenCL backend

### Added
- parallel tuning on single node using Noodles runner
- option to pass new defaults for block dimensions
- option to pass a Python function as code generator
- option to pass custom function for output verification

## [0.1.4] - 2017-06-14
### Changed
- device and kernel name are printed by runner
- tune_kernel also returns a dict with environment info
- using different timer in C vector add example

## [0.1.3] - 2017-04-06
### Changed
- changed how scalar arguments are handled internally

### Added
- separate install and contribution guides

## [0.1.2] - 2017-03-29
### Changed
- allow non-tuple problem_size for 1D grids
- changed default for grid_div_y from None to block_size_y
- converted the tutorial to a Jupyter Notebook
- CUDA backend prints device in use, similar to OpenCL backend
- migrating from nosetests to pytest
- rewrote many of the examples to save results to json files

### Added
- full support for 3D grids, including option for grid_div_z
- separable convolution example

## [0.1.1] - 2017-02-10
### Changed
- changed the output format to list of dictionaries

### Added
- option to set compiler options

## [0.1.0] - 2016-11-02
### Changed
- verbose now also prints debug output when correctness check fails
- restructured the utility functions into util and core
- restructured the code to prepare for different strategies
- shortened the output printed by the tune_kernel
- allowing numpy integers for specifying problem size

### Added
- a public roadmap
- requirements.txt
- example showing GPU code unit testing with the Kernel Tuner
- support for passing a (list of) filenames instead of kernel string
- runner that takes a random sample of 10 percent
- support for OpenCL platform selection
- support for using tuning parameter names in the problem size

## [0.0.1] - 2016-06-14
### Added
- A function to type check the arguments to the kernel
- Example (convolution) that tunes the number of streams 
- Device interface to C functions, for tuning host code
- Correctness checks for kernels during tuning
- Function for running a single kernel instance
- CHANGELOG file
- Compute Cartesian product and process restrictions before main loop
- Python 3.5 compatible code, thanks to Berend
- Support for constant memory arguments to CUDA kernels
- Use of mocking in unittests
- Reporting coverage to codacy
- OpenCL support
- Documentation pages with Convolution and Matrix Multiply examples
- Inspecting device properties at runtime
- Basic Kernel Tuning functionality


# Roadmap for Kernel Tuner

This roadmap presents an overview of the features we are currently planning to
implement. Please note that this is a living document that will evolve as
priorities grow and shift.

## version 0.5.0

 * Allow strategies to tune for a metric other than time

## version 0.9.0

 * Multi-objective optimization

## version 1.0.0

These functions are to be implemented by version 1.0.0, but may already be
implemented in earlier versions.

 * Tuning kernels in parallel on a set of nodes in a GPU cluster
 * Functionality for including auto-tuned kernels in applications

## Wish list

These are the things that we would like to implement, but we currently have no
immediate demand for it. If you are interested in any of these, let us know!

 * Provide API for analysis of tuning results
 * Tuning compiler options in combination with other parameters
 * Example that tunes a kernel using thread block re-indexing
 * Extend Fortran support, no more warnings on data types or missing block size parameter etc.
 * Turn the C backend into a more general compiler backend
 * A get_parameterized_kernel_source function to return the parameterized kernel source for inspection
 * Function to generate wrapper kernels for directly calling/testing device functions

These notebooks are part of the [Kernel Tuner documentation pages](https://benvanwerkhoven.github.io/kernel_tuner/). 

For the materials belonging to the instructor-led tutorial on Kernel Tuner, please see the separate [Kernel Tuner Tutorial repository](https://github.com/benvanwerkhoven/kernel_tuner_tutorial).
Contribution guide
==================
Thank you for considering to contribute to Kernel Tuner!

Reporting Issues
----------------
Not all contributions are code, creating an issue also helps us to improve. When you create an issue about a problem, please ensure the following:

* Describe what you expected to happen.
* If possible, include a minimal example to help us reproduce the issue.
* Describe what actually happened, including the output of any errors printed.
* List the version of Python, CUDA or OpenCL, and C compiler, if applicable. 

Contributing Code
-----------------
For contributing code to Kernel Tuner please select an issue to work on or create a new issue to propose a change or addition. For significant changes, it is required to first create an issue and discuss the proposed changes. Then fork the repository, create a branch, one per change or addition, and create a pull request.

Kernel Tuner follows the Google Python style guide, with Sphinxdoc docstrings for module public functions. Please use `pylint` to check your Python changes.

Before creating a pull request please ensure the following:

* You have written unit tests to test your additions and all unit tests pass
* The examples still work and produce the same (or better) results
* The code is compatible with Python 3.5 or newer
* You have run `pylint` to check your code
* An entry about the change or addition is created in CHANGELOG.md
* Any matching entries in the roadmap.md are updated/removed

If you are in doubt on where to put your additions to the Kernel Tuner, please
have look at the `design documentation
<http://benvanwerkhoven.github.io/kernel_tuner/design.html>`__, or discuss it in the issue regarding your additions.

Development setup
-----------------
You can install the packages required to run the tests using:

.. code-block:: bash

    pip install -e .[dev]

After this command you should be able to run the tests and build the documentation.
See below on how to do that. The ``-e`` flag installs the package in *development mode*.
This means files are not copied, but linked to, such that your installation tracks
changes in the source files.

Running tests
-------------
To run the tests you can use ``pytest -v test/`` in the top-level directory.

Note that tests that require PyCuda and/or a CUDA capable GPU will be skipped if these
are not installed/present. The same holds for tests that require PyOpenCL.

Contributions you make to the Kernel Tuner should not break any of the tests
even if you cannot run them locally.

The examples can be seen as *integration tests* for the Kernel Tuner. Note that
these will also use the installed package.

Building documentation
----------------------
Documentation is located in the ``doc/`` directory. This is where you can type
``make html`` to generate the html pages in the ``doc/build/html`` directory.

The source files used for building the documentation are located in
``doc/source``. The tutorials should be included in the ``tutorials/`` directory
and a symlink can be used to add them to the source file directory before building
documentation.

To update the documentation pages hosted on the GitHub the generated contents of
``doc/build/html`` should be copied to the top-level directory of the
``gh-pages`` branch.
Installation Guide
==================

The Kernel Tuner requires several packages to be installed. First of all, you need a 
working Python version, several Python packages, and optionally CUDA and/or OpenCL 
installations. All of this is explained in detail in this guide.


Python
------

You need a Python installation. I recommend using Python 3 and 
installing it with `Miniconda <https://conda.io/miniconda.html>`__.

Linux users could type the following to download and install Python 3 using Miniconda:

.. code-block:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

You are of course also free to use your own Python installation, and the Kernel Tuner
is developed to be fully compatible with Python 3.6 and newer.

Installing Python Packages
--------------------------

Note that when you are using a native Python installation, the `pip` command used 
Kernel Tuner and its dependencies require `sudo` rights for system wide installation. 

Sudo rights are typically not required when using Miniconda or virtual environments.
You could also use e.g. the `--user` or `--prefix` option of `pip` to install into 
your home directory,
this requires that your home directory is on your `$PYTHONPATH` environment variable
(see for further details the pip documentation).

The following command will install Kernel Tuner together with the required dependencies:

.. code-block:: bash

    pip install kernel_tuner

There are also optional dependencies, explained below.

CUDA and PyCUDA
---------------

Installing CUDA and PyCUDA is optional, because you may want to only use Kernel 
Tuner for tuning OpenCL or C kernels. 

If you want to use the Kernel Tuner to tune 
CUDA kernels you will first need to install the CUDA toolkit 
(https://developer.nvidia.com/cuda-toolkit). A recent version of the 
CUDA toolkit (and the PyCUDA Python bindings for CUDA) are 
recommended (older version may work, but may not support all features of 
Kernel Tuner). 

It's very important that you install the CUDA toolkit before trying to install PyCuda.

You can install PyCuda manually using:

.. code-block:: bash

    pip install pycuda

Or you could install Kernel Tuner and PyCUDA together if you haven't done so already:

.. code-block:: bash

    pip install kernel_tuner[cuda]

If you run into trouble with installing PyCuda, make sure you have CUDA installed first.
Also make sure that the Python package Numpy is already installed, e.g. using `pip install numpy`.

If you retry the ``pip install pycuda`` command, you may need to use the 
``--no-cache-dir`` option to ensure the pycuda installation really starts over and not continues
from an installation that is failing.

If this fails, I recommend to see the PyCuda installation guide (https://wiki.tiker.net/PyCuda/Installation)


OpenCL and PyOpenCL
-------------------

Before we can install PyOpenCL you'll need an OpenCL compiler. There are several 
OpenCL compilers available depending on the OpenCL platform you want to your 
code to run on.

* `AMD APP SDK <http://developer.amd.com/tools-and-sdks/opencl-zone/amd-accelerated-parallel-processing-app-sdk/>`__
* `Intel OpenCL <https://software.intel.com/en-us/iocl_rt_ref>`__
* `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`__
* `Apple OpenCL <https://developer.apple.com/opencl/>`__
* `Beignet <https://www.freedesktop.org/wiki/Software/Beignet/>`__

You can also look at this `OpenCL Installation Guide <https://wiki.tiker.net/OpenCLHowTo>`__ for PyOpenCL.

As with the CUDA toolkit, recent versions of one or more of the above OpenCL SDK's and 
PyOpenCL are recommended to support all features of the Kernel Tuner.

After you've installed your OpenCL compiler of choice you can install PyOpenCL using:

.. code-block:: bash

    pip install pyopencl

Or you could install Kernel Tuner and PyOpenCL together if you haven't done so already:

.. code-block:: bash

    pip install kernel_tuner[opencl]

If this fails, please see the PyOpenCL installation guide (https://wiki.tiker.net/PyOpenCL/Installation)


Installing Kernel Tuner
-----------------------

You can also install from the git repository. This way you also get the 
examples and the tutorials.

.. code-block:: bash

    git clone https://github.com/benvanwerkhoven/kernel_tuner.git
    cd kernel_tuner
    pip install .

You can install Kernel Tuner with several optional dependencies, the full list is:

- `cuda`: install pycuda along with kernel_tuner
- `opencl`: install pycuda along with kernel_tuner
- `doc`: installs packages required to build the documentation
- `tutorial`: install packages required to run the tutorials
- `dev`: install everything you need to start development on Kernel Tuner

For example, use:
```
pip install .[dev,cuda,opencl]
```
To install Kernel Tuner along with all the packages required for development.


Dependencies for the Tutorial
-----------------------------

Some addition Python packages are required to run the tutorial. These packages are
actually very commonly used and chances are that you already have these installed.

However, to install Kernel Tuner along with the dependencies to run the tutorials,
you could use:

.. code-block:: bash

    pip install kernel_tuner[tutorial,cuda]

Or if you have already installed Kernel Tuner and PyCUDA, just use:

.. code-block:: bash

    pip install jupyter matplotlib pandas




Kernel Tuner
============

|Build Status| |CodeCov Badge| |PyPi Badge| |Zenodo Badge| |SonarCloud Badge| |FairSoftware Badge|

Kernel Tuner simplifies the software development of optimized and auto-tuned GPU programs, by enabling Python-based unit testing of GPU code and making it easy to develop scripts for auto-tuning GPU kernels. This also means no extensive changes and no new dependencies are required in the kernel code. The kernels can still be compiled and used as normal from any host programming language.

Kernel Tuner provides a comprehensive solution for auto-tuning GPU programs, supporting auto-tuning of user-defined parameters in both host and device code, supporting output verification of all benchmarked kernels during tuning, as well as many optimization strategies to speed up the tuning process.

Documentation
-------------

The full documentation is available
`here <http://benvanwerkhoven.github.io/kernel_tuner/index.html>`__.

Installation
------------

The easiest way to install the Kernel Tuner is using pip:

To tune CUDA kernels:

- First, make sure you have the `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`_ installed
- Then type: ``pip install kernel_tuner[cuda]``

To tune OpenCL kernels:

- First, make sure you have an OpenCL compiler for your intended OpenCL platform
- Then type: ``pip install kernel_tuner[opencl]``

Or both:

- ``pip install kernel_tuner[cuda,opencl]``

More information about how to install Kernel Tuner and its
dependencies can be found in the `installation guide 
<http://benvanwerkhoven.github.io/kernel_tuner/install.html>`__

Example usage
-------------

The following shows a simple example for tuning a CUDA kernel:

.. code:: python

    kernel_string = """
    __global__ void vector_add(float *c, float *a, float *b, int n) {
        int i = blockIdx.x * block_size_x + threadIdx.x;
        if (i<n) {
            c[i] = a[i] + b[i];
        }
    }
    """

    size = 10000000

    a = numpy.random.randn(size).astype(numpy.float32)
    b = numpy.random.randn(size).astype(numpy.float32)
    c = numpy.zeros_like(b)
    n = numpy.int32(size)
    args = [c, a, b, n]

    tune_params = dict()
    tune_params["block_size_x"] = [32, 64, 128, 256, 512]

    tune_kernel("vector_add", kernel_string, size, args, tune_params)

The exact same Python code can be used to tune an OpenCL kernel:

.. code:: python

    kernel_string = """
    __kernel void vector_add(__global float *c, __global float *a, __global float *b, int n) {
        int i = get_global_id(0);
        if (i<n) {
            c[i] = a[i] + b[i];
        }
    }
    """

The Kernel Tuner will detect the kernel language and select the right compiler and 
runtime. For every kernel in the parameter space, the Kernel Tuner will insert C 
preprocessor defines for the tunable parameters, compile, and benchmark the kernel. The 
timing results will be printed to the console, but are also returned by tune_kernel to 
allow further analysis. Note that this is just the default behavior, what and how 
tune_kernel does exactly is controlled through its many `optional arguments 
<http://benvanwerkhoven.github.io/kernel_tuner/user-api.html#kernel_tuner.tune_kernel>`__.

You can find many - more extensive - example codes, in the
`examples directory <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/>`__
and in the `Kernel Tuner
documentation pages <http://benvanwerkhoven.github.io/kernel_tuner/index.html>`__.

Search strategies for tuning
----------------------------

Kernel Tuner supports many optimization algorithms to accelerate the auto-tuning process. Currently 
implemented search algorithms are: Brute Force (default), Nelder-Mead, Powell, CG, BFGS, L-BFGS-B, TNC, 
COBYLA, SLSQP, Random Search, Basinhopping, Differential Evolution, a Genetic Algorithm, Particle Swarm 
Optimization, the Firefly Algorithm, and Simulated Annealing.

.. image:: doc/gemm-amd-summary.png
    :width: 100%
    :align: center

Using a search strategy is easy, you only need to specify to ``tune_kernel`` which strategy and method 
you would like to use, for example ``strategy="genetic_algorithm"`` or ``strategy="basinhopping"``. 
For a full overview of the supported search strategies and methods please see the `user 
api documentation <http://benvanwerkhoven.github.io/kernel_tuner/user-api.html>`__.

Tuning host and kernel code
---------------------------

It is possible to tune for combinations of tunable parameters in
both host and kernel code. This allows for a number of powerfull things,
such as tuning the number of streams for a kernel that uses CUDA Streams
or OpenCL Command Queues to overlap transfers between host and device
with kernel execution. This can be done in combination with tuning the
parameters inside the kernel code. See the `convolution\_streams example
code <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/>`__
and the
`documentation <http://benvanwerkhoven.github.io/kernel_tuner/hostcode.html>`__
for a detailed explanation of the kernel tuner Python script.


Correctness verification
------------------------

Optionally, you can let the kernel tuner verify the output of every
kernel it compiles and benchmarks, by passing an ``answer`` list. This
list matches the list of arguments to the kernel, but contains the
expected output of the kernel. Input arguments are replaced with None.

.. code:: python

    answer = [a+b, None, None]  # the order matches the arguments (in args) to the kernel
    tune_kernel("vector_add", kernel_string, size, args, tune_params, answer=answer)

Contributing
------------

Please see the `Contributions Guide <http://benvanwerkhoven.github.io/kernel_tuner/contributing.html>`__.

Citation
--------
If you use Kernel Tuner in research or research software, please cite the most relevant among the following publications:

.. code:: latex

    @article{kerneltuner,
      author  = {Ben van Werkhoven},
      title   = {Kernel Tuner: A search-optimizing GPU code auto-tuner},
      journal = {Future Generation Computer Systems},
      year = {2019},
      volume  = {90},
      pages = {347-358},
      url = {https://www.sciencedirect.com/science/article/pii/S0167739X18313359},
      doi = {https://doi.org/10.1016/j.future.2018.08.004}
    }

    @article{willemsen2021bayesian,
      author = {Willemsen, Floris-Jan and Van Nieuwpoort, Rob and Van Werkhoven, Ben},
      title = {Bayesian Optimization for auto-tuning GPU kernels},
      journal = {International Workshop on Performance Modeling, Benchmarking and Simulation
         of High Performance Computer Systems (PMBS) at Supercomputing (SC21)},
      year = {2021},
      url = {https://arxiv.org/abs/2111.14991}
    }


Related work
------------

You may also like `CLTune <https://github.com/CNugteren/CLTune>`__ by
Cedric Nugteren. CLTune is a C++ library for kernel tuning.


.. |Build Status| image:: https://github.com/benvanwerkhoven/kernel_tuner/actions/workflows/python-app.yml/badge.svg
   :target: https://github.com/benvanwerkhoven/kernel_tuner/actions/workflows/python-app.yml
.. |CodeCov Badge| image:: https://codecov.io/gh/benvanwerkhoven/kernel_tuner/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/benvanwerkhoven/kernel_tuner
.. |PyPi Badge| image:: https://img.shields.io/pypi/v/kernel_tuner.svg?colorB=blue 
   :target: https://pypi.python.org/pypi/kernel_tuner/
.. |Zenodo Badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1220113.svg
   :target: https://doi.org/10.5281/zenodo.1220113
.. |SonarCloud Badge| image:: https://sonarcloud.io/api/project_badges/measure?project=benvanwerkhoven_kernel_tuner&metric=alert_status
   :target: https://sonarcloud.io/dashboard?id=benvanwerkhoven_kernel_tuner
.. |FairSoftware Badge| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu
Kernel Tuner Examples
=====================

Most of the examples show how to use Kernel Tuner to tune a
CUDA, OpenCL, or C kernel, while demonstrating a particular usecase of Kernel Tuner.

Except for `test\_vector\_add.py <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/test_vector_add.py>`__  and 
`test\_vector\_add_parameterized.py <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/test_vector_add_parameterized.py>`__,
which show how to write tests for GPU kernels with Kernel Tuner.

Below we list the example applications and the features they illustrate.

Vector Add
----------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/vector_add.py>`__] [`CUDA-C++ <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda-c++/vector_add.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/vector_add.py>`__] [`C <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/c/vector_add.py>`__] [`Fortran <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/fortran/vector_add.py>`__]
 - use Kernel Tuner to tune a simple kernel

Stencil
-------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/stencil.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/stencil.py>`__]
 -  use a 2-dimensional problem domain with 2-dimensional thread blocks in a simple and clean example

Matrix Multiplication
---------------------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/matmul.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/matmul.py>`__]
 -  pass a filename instead of a string with code
 -  use 2-dimensional thread blocks and tiling in both dimensions
 -  tell Kernel Tuner to compute the grid dimensions for 2D thread blocks with tiling
 -  use the restrictions option to limit the search to only valid configurations
 -  use a user-defined performance metric like GFLOP/s

Convolution
-----------
There are several different examples centered around the convolution
kernel [`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/convolution.cu>`__]
[`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/convolution.cl>`__]

convolution.py
~~~~~~~~~~~~~~
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/convolution.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/convolution.py>`__]
 - use tunable parameters for tuning for multiple input sizes
 - pass constant memory arguments to the kernel
 - write output to a json file

sepconv.py
~~~~~~~~~~
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/sepconv.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/sepconv.py>`__]
 - use the convolution kernel for separable filters
 - write output to a csv file using Pandas

convolution\_correct.py
~~~~~~~~~~~~~~~~~~~~~~~
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/convolution_correct.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/convolution_correct.py>`__]
 - use run\_kernel to compute a reference answer
 - verify the output of every benchmarked kernel

convolution\_streams.py
~~~~~~~~~~~~~~~~~~~~~~~
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/convolution_streams.py>`__]
 - allocate page-locked host memory from Python
 - overlap transfers to and from the GPU with computation
 - tune parameters in the host code in combination with those in the kernel
 - use the lang="C" option and set compiler options
 - pass a list of filenames instead of strings with kernel code

Reduction
---------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/reduction.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/reduction.py>`__]
 - use vector types and shuffle instructions (shuffle is only available in CUDA)
 - tune the number of thread blocks the kernel is executed with
 - tune the partial loop unrolling factor of a for-loop
 - tune pipeline that consists of two kernels
 - tune with custom output verification function

Sparse Matrix Vector Multiplication
-----------------------------------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/spmv.py>`__]
 -  use scipy to compute a reference answer and verify all benchmarked kernels
 -  express that the number of thread blocks depends on the values of tunable parameters

Point-in-Polygon
----------------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/pnpoly.py>`__]
 -  overlap transfers with device mapped host memory
 -  tune on different implementations of an algorithm

ExpDist
-------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/expdist.py>`__]
 -  in-thread block 2D reduction using CUB library
 -  C++ in CUDA kernel code
 -  tune multiple kernels in pipeline

Code Generator
--------------
[`CUDA <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/cuda/vector_add_codegen.py>`__] [`OpenCL <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/opencl/vector_add_codegen.py>`__]
 - use a Python function as a code generator

.. toctree::
   :maxdepth: 2

.. _contributing:


.. include:: ../../CONTRIBUTING.rst

.. toctree::
   :maxdepth: 2

.. _details:

API Documentation 
=================

This file provides all the details you need about how to call the Kernel Tuner's functions, including all the optional arguments.

.. autofunction:: kernel_tuner.tune_kernel

.. autofunction:: kernel_tuner.run_kernel

.. autofunction:: kernel_tuner.store_results

.. autofunction:: kernel_tuner.create_device_targets
.. toctree::
   :maxdepth: 2

.. _examples:


.. include:: ../../examples/README.rst

.. highlight:: python
    :linenothreshold: 50


Parameter Vocabulary
--------------------

There are certain tunable parameters that have special meaning in Kernel Tuner.
This document specifies which parameters are special and what there uses are when auto-tuning GPU kernels.

In general, it is best to avoid using these parameter names for purposes other than the ones indicated in this document.

.. code-block:: python

    kernel_tuner #is inserted by Kernel Tuner to signal the code is compiled using the tuner

    block_size_* #reserved for thread block dimensions
    grid_size_* #reserved for grid dimensions, if you want to tune these use problem_size

    compiler_opt_* #reserved for future support for tuning compiler options

    loop_unroll_factor_* #reserved for tunable parameters that specify loop unrolling factors

    nvml_* #reserved for tunable parameters and outputs related to NVML
    nvml_pwr_limit #use NVML to set power limit
    nvml_gr_clock #use NVML to set graphics clock
    nvml_mem_clock #use NVML to set memory clock


There are also a number of names that Kernel Tuner uses for reporting benchmarking results. 
Because these are reported along with the tunable parameters, it is generally a good idea to not use these names for any tunable parameters.

.. code-block:: python

    time* #reserved for time measurements

    Information that can be observed using kernel_tuner.nvml.NVMLObserver:
    nvml_energy
    nvml_power
    power_readings
    core_freq
    mem_freq
    temperature

    ps_energy  # Energy as measured by PowerSensor
    ps_power   # Power as measured by PowerSensor


.. highlight:: python
    :linenothreshold: 5



Tuning Host Code
----------------

With the Kernel Tuner it is also possible to tune the host code of your GPU programs, or even just any C function for that matter.
Tuning host code can be useful when it contains parameters that have impact on the performance of kernel on the GPU, such as the number of
streams to use when executing a kernel across multiple streams. Another example is when you want to include the data transfers between
host and device into your tuning setup, or tune for different methods of moving data between host and device.

There are few differences with tuning just a single CUDA or OpenCL kernel, to list them:  
 * You have to specify the lang="C" option
 * The C function should return a ``float``
 * You have to do your own timing and error handling in C

You have to specify the language as "C" because the Kernel Tuner will be calling a host function. This means that the Kernel
Tuner will have to interface with C and in fact uses a different backend. This also means you can use this way of tuning
without having PyCuda installed, because the C functions interface calls the CUDA compiler directly.

The C function should return a float, this is the convention used by the Kernel Tuner. The returned float is also the number
that you are tuning for. Meaning that this does not necessarily needs to be time, you could also optimize a program for
a different quality, as long as you can express that quality in a single floating-point value. When benchmarking an instance
of the parameter space the returned floats will be averaged for the multiple runs in the same way as with direct CUDA or OpenCL kernel tuning.

By itself the C language does not provide any very precise timing functions. If you are tuning the host code of a CUDA program you can use
CUDA Events to do the timing for you. However, if you are using plain C then you have to supply your own timing function.
In the `C vector add example <https://github.com/benvanwerkhoven/kernel_tuner/blob/master/examples/c/vector_add.py>`__ we are using the ``omp_get_wtime()`` function from OpenMP to measure time on the CPU.

Tuning the number of streams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following describes the example in ``examples/cuda/convolution_streams.py``.
In this example, the same convolution kernel is used as with correctness checking and convolution application example.

What is different is that we also supply the host code, which you can find in ``examples/cuda/convolution_streams.cu``. It is a bit
too long and complex to include here, but we will explain what it does. The idea behind the host code is that the kernel computation
is spread across a number of CUDA streams. In this way, it is possible to overlap the data transfers from host to device with kernel execution, and with
transfers from the device back to the host.

The way we split the computation across streams is by dividing the problem in the y-dimension into chunks. The data transferred by the first stream is slightly 
larger to account for the overlapping border between the data needed by different streams. Before the kernel in stream `n` can start executing the data transfers 
in streams `n` and `n-1` have to be finished. To ensure the latter, we use CUDA Events and in particular cudaStreamWaitEvent(), which halts stream `n` until the 
transfer in stream `n-1` has finished.

The way you use the Kernel Tuner to tune this CUDA program is very similar to when you are tuning a CUDA kernel directly, as you can see below:

.. code-block:: python

    with open('convolution_streams.cu', 'r') as f:
        kernel_string = f.read()

    problem_size = (4096, 4096)
    size = numpy.prod(problem_size)
    input_size = (problem_size[0]+16) * (problem_size[0]+16)

    output = numpy.zeros(size).astype(numpy.float32, order='C')
    input = numpy.random.randn(input_size).astype(numpy.float32, order='C')
    filter = numpy.random.randn(17*17).astype(numpy.float32)
    args = [output, input, filter]

    tune_params = dict()
    tune_params["block_size_x"] = [16*i for i in range(1,9)]
    tune_params["block_size_y"] = [2**i for i in range(6)]

    tune_params["tile_size_x"] = [2**i for i in range(3)]
    tune_params["tile_size_y"] = [2**i for i in range(3)]

    tune_params["num_streams"] = [2**i for i in range(6)]

    grid_div_x = ["block_size_x", "tile_size_x"]
    grid_div_y = ["block_size_y", "tile_size_y", "num_streams"]

    kernel_tuner.tune_kernel("convolution_streams", kernel_string,
        problem_size, args, tune_params,
        grid_div_y=grid_div_y, grid_div_x=grid_div_x, verbose=True, lang="C")

In fact, the only differences with the simple convolution example are:  
 * The source file also contains host code 
 * "num_streams" is added to the tuning parameters
 * "num_streams" is added to the "grid_div_y" list
 * The kernel_name "convolution_streams" is a C function
 * lang="C" is used to tell this is a C function
 * ``filter`` is not passed as a constant memory argument

Most differences have been explained, but we clarify a few things below.

The function that we are tuning is a C function that launches the CUDA kernel by itself, yet we supply the grid_div_x and 
grid_div_y lists. We are, however, not required to do so. The C function could just compute the grid dimensions in whatever way it sees fit. Using grid_div_y 
and grid_div_x at this point is matter of choice. To support this convenience, the values grid_size_x and grid_size_y are inserted by the Kernel Tuner into the 
compiled C code. This way, you don't have to compute the grid size in C, you can just use the grid dimensions as computed by the Kernel Tuner.

The filter is not passed separately as a constant memory argument, because the CudaMemcpyToSymbol operation is now performed by the C host function. Also, 
because the code is compiled differently, we have no direct reference to the compiled module that is uploaded to the device and therefore we can not perform this 
operation directly from Python. If you are tuning host code, you have to perform all memory allocations, frees, and memcpy operations inside the C host code, 
that's the purpose of host code after all. That is also why you have to do the timing yourself in C, as you may not want to include the time spent on memory 
allocations and other setup into your time measurements.





.. highlight:: python
    :linenothreshold: 5



Kernel Correctness Verification
-------------------------------

Whenever you optimize a program for performance it is very important to
ensure that the program is still producing the correct output. What good
is a program that is fast but not correct?

Therefore an important feature of the kernel tuner is to verify the output
of every kernel instance in the parameter space. To use the kernel tuner
with correctness checking you need to pass the ``answer`` option to
``tune_kernel()``. Answer is a list that should match the order and types of
the kernel arguments. However, if an argument to the kernel is input-only
you may insert ``None`` at that location in the list.

After kernel compilation, but before benchmarking the kernel, the kernel
tuner runs the kernel once to verify the output it produces. For each
argument in the ``answer`` list that is not None, it will check the results
produced by the current kernel against the expected result specified in
``answer``. The comparison is currently implemented using numpy.allclose()
with an maximum allowed absolute error of 1e-6. If you want to use a 
difference tolerance value, use the optional argument ``atol``.

The example in ``examples/cuda/convolution_correct.py`` demonstrates how
to use the ``answer`` option of ``tune_kernel()``:

.. code-block:: python

    import numpy
    import kernel_tuner

    with open('convolution.cu', 'r') as f:
        kernel_string = f.read()

    problem_size = (4096, 4096)
    size = numpy.prod(problem_size)
    input_size = ((problem_size[0]+16) * (problem_size[1]+16))

    output = numpy.zeros(size).astype(numpy.float32)
    input = numpy.random.randn(input_size).astype(numpy.float32)

    filter = numpy.random.randn(17*17).astype(numpy.float32)
    cmem_args= {'d_filter': filter }

    args = [output, input, filter]
    tune_params = dict()
    tune_params["block_size_x"] = [16*i for i in range(1,9)]
    tune_params["block_size_y"] = [2**i for i in range(6)]

    tune_params["tile_size_x"] = [2**i for i in range(3)]
    tune_params["tile_size_y"] = [2**i for i in range(3)]

    grid_div_x = ["block_size_x", "tile_size_x"]
    grid_div_y = ["block_size_y", "tile_size_y"]

    #compute the answer using a naive kernel
    params = { "block_size_x": 16, "block_size_y": 16 }
    results = kernel_tuner.run_kernel("convolution_naive", kernel_string,
        problem_size, args, params,
        grid_div_y=["block_size_y"], grid_div_x=["block_size_x"])

    #set non-output fields to None
    answer = [results[0], None, None]

    #start kernel tuning with correctness verification
    kernel_tuner.tune_kernel("convolution_kernel", kernel_string,
        problem_size, args, tune_params,
        grid_div_y=grid_div_y, grid_div_x=grid_div_x,
        verbose=True, cmem_args=cmem_args, answer=answer)

This example uses the ``run_kernel()`` function of the kernel tuner
to run a single kernel and return its results, with almost the same
interface as ``tune_kernel()``. In this example we run a naive CUDA
kernel whose results are trusted to be correct.

The ``answer`` list is constructed out of the results from the naive
kernel, but only includes the kernel arguments that are actually outputs.
The arguments that are input are replaced by a ``None`` value in the
``answer`` list before the list is passed to ``tune_kernel()``.

There are cases, however, where simply comparing the results computed on the device to precomputed values is not enough,
and more flexibility is necessary.
In this case, it is possible to use the ``verify`` option of ``tune_kernel()`` and specify a ``callable`` object that
implements a user-defined correctness check.
This function should accept three parameters: ``cpu_result``, ``gpu_result``, and ``atol``.
Although the name of the parameters can be different, their semantic is position dependent and reflected in the names
used in the documentation.

The example in ``examples/cuda/reduction.py`` demonstrates how to use the ``verify`` option of ``tune_kernel()``;
what follows is a snippet from the example:

.. code-block:: python

    # gpu_result
    args = [sum_x, x, n]
    # cpu_result
    reference = [numpy.sum(x), None, None]
    # custom verify function
    def verify_partial_reduce(cpu_result, gpu_result, atol=None):
        return numpy.isclose(cpu_result, numpy.sum(gpu_result), atol=atol)
    # call to tune_kernel()
    first_kernel, _ = tune_kernel("sum_floats", kernel_string, problem_size,
        args, tune_params, grid_div_x=[], verbose=True, answer=reference, verify=verify_partial_reduce)

The first argument, ``cpu_result``, is mapped to the NumPy array provided to the ``answer`` option; in this example it
is mapped to ``reference``.
The second argument, ``gpu_result``, is mapped to the NumPy array provided to the ``arguments`` option of
``tune_kernel()``; in this example it is mapped to ``args``.
The third argument, ``atol``, is set to ``None``; the default maximum allowed absolute error of 1e-6 is then used.

In the example, the user-defined ``verify`` function is used to compare the partial results, computed on the GPU,
to the final result, computed on the CPU.
The same could not be achieved just by using the ``answer`` option, because the number of elements in ``args[0]`` does
not necessarily match the number of elements in ``reference[0]`` in this example... toctree::
   :maxdepth: 2


Design documentation
====================

This section provides detailed information about the design and internals 
of the Kernel Tuner. **This information is mostly relevant for developers.**

The Kernel Tuner is designed to be extensible and support 
different search and execution strategies. The current architecture of 
the Kernel Tuner can be seen as:

.. image:: design.png
   :width: 500pt

At the top we have the kernel code and the Python script that tunes it, 
which uses any of the main functions exposed in the user interface.

The strategies are responsible for iterating over and searching through 
the search space. The default strategy is ``brute_force``, which 
iterates over all valid kernel configurations in the search space. 
``random_sample`` simply takes a random sample of the search space. More 
advanced strategies currently implemented in Kernel Tuner are 
``minimize``, ``basinhopping``, and differential evolution 
(``diff_evo``). How to use these is explained in the :doc:`user-api`,
see the options ``strategy`` and ``strategy_options``.

The runners are responsible for compiling and benchmarking the kernel 
configurations selected by the strategy. The sequential runner is currently
the only supported runner, which does exactly what its name says. It compiles 
and benchmarks configurations using a single sequential Python process.
Other runners are foreseen in future releases.

The runners are implemented on top of a high-level *Device Interface*,
which wraps all the functionality for compiling and benchmarking
kernel configurations based on the low-level *Device Function Interface*.
Currently, we have 
four different implementations of the device function interface, which 
basically abstracts the different backends into a set of simple 
functions such as ``ready_argument_list`` which allocates GPU memory and 
moves data to the GPU, and functions like ``compile``, ``benchmark``, or 
``run_kernel``. The functions in the core are basically the main 
building blocks for implementing runners.

At the bottom, three of the backends are shown. 
PyCUDA and PyOpenCL are for tuning either CUDA or OpenCL kernels.
A relatively new addition is the Cupy backend based on Cupy for tuning
CUDA kernels using the NVRTC compiler.
The C 
Functions implementation can actually call any compiler, typically NVCC 
or GCC is used. This backend was created not just to be able to tune C 
functions, but mostly to tune C functions that in turn launch GPU kernels.

The rest of this section contains the API documentation of the modules 
discussed above. For the documentation of the user API see the 
:doc:`user-api`.



Strategies
----------

kernel_tuner.strategies.brute_force
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.brute_force
    :members:

kernel_tuner.strategies.random_sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.random_sample
    :members:

kernel_tuner.strategies.minimize
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.minimize
    :members:

kernel_tuner.strategies.basinhopping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.basinhopping
    :members:

kernel_tuner.strategies.diff_evo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.diff_evo
    :members:

kernel_tuner.strategies.genetic_algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.genetic_algorithm
    :members:

kernel_tuner.strategies.pso
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.pso
    :members:

kernel_tuner.strategies.firefly_algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.firefly_algorithm
    :members:

kernel_tuner.strategies.simulated_annealing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.strategies.simulated_annealing
    :members:



Runners
-------

kernel_tuner.runners.sequential.SequentialRunner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.runners.sequential.SequentialRunner
    :special-members: __init__
    :members:

kernel_tuner.runners.sequential.SimulationRunner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.runners.simulation.SimulationRunner
    :special-members: __init__
    :members:


Device Interfaces
-----------------

kernel_tuner.core.DeviceInterface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.core.DeviceInterface
    :special-members: __init__
    :members:

kernel_tuner.cuda.CudaFunctions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.cuda.CudaFunctions
    :special-members: __init__
    :members:

kernel_tuner.cupy.CupyFunctions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.cupy.CupyFunctions
    :special-members: __init__
    :members:

kernel_tuner.opencl.OpenCLFunctions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.opencl.OpenCLFunctions
    :special-members: __init__
    :members:

kernel_tuner.c.CFunctions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kernel_tuner.c.CFunctions
    :special-members: __init__
    :members:


Util Functions
--------------

kernel_tuner.util
~~~~~~~~~~~~~~~~~
.. automodule:: kernel_tuner.util
    :members:

.. highlight:: python
    :linenothreshold: 5



Templated kernels
-----------------

It is quite common in CUDA programming to write kernels that use C++ templates. This can be very useful when writing code that can work for several types, for example floats and doubles. However, the use of C++ templates makes it slightly more difficult to directly 
integrate the CUDA kernel into applications that are not written in C++, for example Matlab, Fortran, or Python. And since Kernel Tuner is written in Python, we needed to take a few extra steps to provide support for templated CUDA kernels. Let's first look at an 
example of what it's like to tune a templated kernel with Kernel Tuner.

Example
~~~~~~~

Say we have a templated CUDA kernel in a file called vector_add.cu:

.. code-block:: cuda

    template<typename T>
    __global__ void vector_add(T *c, T *a, T *b, int n) {
        auto i = blockIdx.x * block_size_x + threadIdx.x;
        if (i<n) {
            c[i] = a[i] + b[i];
        }
    }

Then the Python script to tune this kernel would be as follows:

.. code-block:: python

    import numpy
    from kernel_tuner import tune_kernel

    size = 1000000

    a = numpy.random.randn(size).astype(numpy.float32)
    b = numpy.random.randn(size).astype(numpy.float32)
    c = numpy.zeros_like(b)
    n = numpy.int32(size)

    args = [c, a, b, n]

    tune_params = dict()
    tune_params["block_size_x"] = [128+64*i for i in range(15)]

    tune_kernel("vector_add<float>", "vector_add.cu", size, args, tune_params)

What you can see is that in the Python code we specify the template instantiation to use. Kernel Tuner will detect the use of templated kernels when the kernel_name positional argument to tune_kernel contains a template argument.

This feature also allows use to auto-tune template parameters to the kernel. We could for example define a tunable parameter:

.. code-block:: python

    tune_params["my_type"] = ["float", "double"]

and call tune_kernel using a tunable parameter inside the template arguments:

.. code-block:: python

    tune_kernel("vector_add<my_type>", "vector_add.cu", size, args, tune_params)

Selecting a backend
~~~~~~~~~~~~~~~~~~~

Kernel Tuner supports multiple backends, for CUDA these are based on PyCUDA and Cupy. The following explains how to enable tuning of templated kernels with either backend.

The PyCuda backend is the default backend in Kernel Tuner and is selected if the user does not supply the 'lang' option and CUDA code is detected in the kernel source, or when lang is set to "CUDA" by the user. PyCuda requires CUDA kernels to have extern C linkage, 
which means that C++ templated kernels are not supported. To support templated kernels regardless of this limitation Kernel Tuner attempts to wrap the templated CUDA kernel by inserting a compile-time template instantiation statement and a wrapper kernel that calls 
the templated CUDA kernel, which is actually demoted to a __device__ function in the process. These automatic code rewrites have a real risk of breaking the code. To minimize the chance of errors due to Kernel Tuner's automatic code rewrites, it's best to isolate the 
templated kernel in a single source file and include it where needed in the larger application.

The Cupy backend provides much more advanced support for C++ templated kernels, because it internally uses NVRTC, the Nvidia runtime compiler. NVRTC does come with some restrictions however, for example NVRTC does not allow any host code to be inside code that
is passed. So, like with the PyCuda backend it helps to separate the source code of device and host functions into seperate files. You can force Kernel Tuner to use the Cupy backend by passing the lang="cupy" option to tune_kernel. 






.. kernel_tuner documentation master file, created by
   sphinx-quickstart on Tue Mar 29 15:46:32 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Kernel Tuner documentation
==============================

Contents:

.. toctree::
   :maxdepth: 1

   Introduction <self>
   install
   convolution
   diffusion
   examples
   matrix_multiplication
   correctness
   hostcode
   templates
   user-api
   vocabulary
   design
   contributing

Introduction
============

.. include:: ../../README.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2

.. _install:


.. include:: ../../INSTALL.rst

