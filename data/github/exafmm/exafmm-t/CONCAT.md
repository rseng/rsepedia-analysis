---
title: 'ExaFMM: a high-performance fast multipole method library with C++ and Python interfaces'
tags:
  - C++
  - Python
  - fast multipole method
  - low-rank approximation
  - high performance computing
authors:
 - name: Tingyu Wang
   orcid: 0000-0003-2520-0511
   affiliation: 1
 - name: Rio Yokota
   orcid: 0000-0001-7573-7873
   affiliation: 2
 - name: Lorena A. Barba
   orcid: 0000-0001-5812-2711
   affiliation: 1
affiliations:
 - name: The George Washington University
   index: 1
 - name: Tokyo Institute of Technology
   index: 2
date: 21 July 2020
bibliography: paper.bib
---

# Summary

ExaFMM is an open-source library for fast multipole algorithms, providing high-performance evaluation of N-body problems in three dimensions, with C++ and Python interfaces.
This new implementation is the most recent of many across a decade of work in our research group.
Our goal for all these years has been to produce reusable, standard code for this intricate algorithm. 
The new header-only C/C++ implementation is easier to extend, still competitive with state-of-the-art codes, and includes a `pybind11` Python interface [@pybind11].

The fast multipole method (FMM) was introduced more than 30 years ago [@GreengardRokhlin1987]
as a means of reducing 
the complexity of N-body problems from $\mathcal{O}(N^2)$ to $\mathcal{O}(N)$
using hierarchical approximations of long-range interactions.
Two major variants of hierarchical N-body algorithms exist: treecodes and FMM.
(Algebraic analogues also exist, such as H-matrix methods.)
Both were originally developed for fast evaluation of the gravitational potential field, but now have found many applications in different fields.
For example, the integral formulation of problems modeled by elliptic partial differential equations can be reinterpreted as N-body interactions.
In this way, N-body algorithms are applicable to acoustics, electromagenetics, fluid dynamics, and more.
The present version of ExaFMM implements the so-called kernel-independent variant of FMM [@yingKernelindependentAdaptiveFast2004].
It supports computing both the potential and forces of Laplace, low-frequency Helmholtz and modified Helmholtz (a.k.a. Yukawa) kernels.
Users can add other non-oscillatory kernels with modest programming effort.

# Statement of Need

Over the past three decades, a plethora of fast N-body implementations have emerged.
We will mention a few notable ones for context.
`Bonsai` [@bedorfSparseOctreeGravitational2012] is a gravitational treecode that runs entirely on GPU hardware.
`ChaNGa` [@jetleyMassivelyParallelCosmological2008] is also a treecode that uses `Charm++` to automate dynamic load balancing.
`ScalFMM` [@agullo2014task] implements the black-box FMM, a kernel-independent variant based on interpolation.
It comes with an option to use `StarPU` runtime system to handle heterogeneous task scheduling.
`TBFMM` [@bramasTBFMMGenericParallel2020] is a task-based FMM library that features a generic C++ design to support various types of tree structures and kernels, through heavy use of C++ templates.
`PVFMM` [@malhotraPVFMMParallelKernel2015] can compute both particle and volume potentials using a kernel-independent FMM, KIFMM [@yingKernelindependentAdaptiveFast2004].

The first version of ExaFMM focused on low-accuracy optimizations and implemented a dual-tree traversal [@BarbaYokota2012-figshare; @YokotaBarba2011a; @yokotaFMMBasedDual2013].
It was GPU-enabled using CUDA, parallel with MPI and exploited multithreading using OpenMP.
Despite all these efforts, it has remained a challenge in the FMM community to have a well-established open-source software package, analogous to FFTW for the fast Fourier transform,
delivering compelling performance with a standard and easy-to-use interface.

The "alpha" version of ExaFMM is long and complex, and hard to maintain.
Its length is partly due to a focus on fast execution, which led to specialized treatment of low-$p$ evaluations and extensive hand-coded optimizations.
This emphasis on achieving high performance led to poor reusability and maintainability.
The current version uses the kernel-independent variant of the method for higher performance at high accuracy (higher truncation order $p$).
Keeping focus on maintainability, it stays close to standards, with clean function interfaces, shallow inheritance and conservative use of classes.
Instead of a fixation with being the fastest of all, we focus on "ballpark" competitive performance.
Above all, the Python API and ability to call it from Jupyter notebooks should make this software valuable for many applications.

# Features of the software design

`exafmm-t` is designed to be standard and lean.
First, it only uses C++ STL containers and depends on mature math libraries: BLAS, LAPACK and FFTW3.
Second, `exafmm-t` is moderately object-oriented, namely, the usage of encapsulation, inheritance and polymorphism is conservative or even minimal in the code.
The core library consists of around 6,000 lines of code, which is an order of magnitude shorter than many other FMM packages.

`exafmm-t` is concise but highly optimized.
To achieve competitive performance, our work combines techniques and optimizations from several past efforts.
On top of multi-threading using OpenMP, we further speed up the P2P operator (near-range interactions) using SIMD vectorization with SSE/AVX/AVX-512 compatibility;
we apply the cache optimization proposed in PVFMM [@malhotraPVFMMParallelKernel2015] to improve the performance of M2L operators (far-range interactions).
In addition, `exafmm-t` also allows users to pre-compute and store translation operators, which benefits applications that require iterative FMM evaluations.
The single-node performance of `exafmm-t` is on par with the state-of-the-art packages that we mentioned above.
We ran a benchmark that solves a Laplace N-body problem with 1 million randomly distributed particles on a workstation with a 14-core Intel i9-7940X CPU.
It took 0.95 and 1.48 seconds to obtain 7 and 10 digits of accuracy on the potential, respectively.

`exafmm-t` is also relatively easy to extend.
Adding a new kernel only requires users to create a derived `FMM` class and provide the kernel function.
Last but not least, it offers high-level Python APIs to support Python applications.
Thanks to `pybind11`, most STL containers can be automatically converted to Python-native data structures.
Since Python uses duck typing, we have to expose overloaded functions to different Python objects.
To avoid naming collisions and keep a clean interface, we chose to create a Python module for each kernel under `exafmm-t`'s Python package, instead of adding suffixes to function and class names to identify types.

# Application

We have recently integrated `exafmm-t` with `Bempp-cl`, an open-source boundary element method (BEM) package in Python [@Betcke2021],
whose predecessor, `BEM++` [@smigajSolvingBoundaryIntegral2015], has enabled many acoustic and electromagnetic applications.
In BEM applications, computations are dominated by the dense matrix-vector multiplication (mat-vec) in each iteration.
`exafmm-t` reduces both time and memory cost of mat-vec to a linear complexity, thus makes `Bempp-cl` feasible to solve large-scale problems.
In an upcoming paper, we demonstrate the capabilities and performance of Bempp-Exafmm on biomolecular electrostatics simulations, including solving problems at the scale of a virus [@wangETal2021].
The showcase calculation in that paper (submitted) obtains the surface electrostatic potential of a Zika virus, modeled with 1.6 million atoms, 10 million boundary elements (30M points), at a runtime of 1.5 hours on 1 CPU node.

# Acknowledgements

This software builds on efforts over more than a decade. It directly received support from grants to LAB in both the UK and the US, including EPSRC Grant EP/E033083/1, and NSF Grants OCI-0946441, NSF CAREER OAC-1149784, and CCF-1747669.
Other support includes faculty start-up funds at Boston University and George Washington University, and NVIDIA via hardware donations. 

# References
# exafmm-t 

[![Build Status](https://travis-ci.com/exafmm/exafmm-t.svg?branch=master)](https://travis-ci.com/exafmm/exafmm-t)
[![status](https://joss.theoj.org/papers/0faabca7e0ef645b42d7dd72cc924ecc/status.svg)](https://joss.theoj.org/papers/0faabca7e0ef645b42d7dd72cc924ecc)

## Cite as

"ExaFMM: a high-performance fast multipole method library with C++ and Python interfaces", Tingyu Wang, Rio Yokota, Lorena A. Barba. The Journal of Open Source Software, 6(61):3145 (2021). doi:10.21105/joss.03145 

**exafmm-t** is a kernel-independent fast multipole method library for solving N-body problems.
It provides both C++ and Python APIs.
We use [pybind11](https://github.com/pybind/pybind11) to create Python bindings from C++ code.
exafmm-t aims to deliver compelling performance with a simple code design and a user-friendly interface.
It currently supports both potential and force calculation of Laplace, low-frequency Helmholtz and modified Helmholtz (Yukawa) kernel in 3D.
In addition, users can easily add other non-oscillatory kernels under exafmm-t's framework.

## Documentation

The full documentation is available [here](https://exafmm.github.io/exafmm-t).

Please use [GitHub issues](https://github.com/exafmm/exafmm-t/issues) for tracking bugs and requests.
To contribute to exafmm-t, please review [CONTRIBUTING](https://github.com/exafmm/exafmm-t/blob/master/CONTRIBUTING.md).
# Contributing guidelines

To contribute to exafmm-t, please first fork the repository and push changes on your fork and then submit a pull request (PR).
For minor fixes, please make sure that your code passes all [current tests](https://exafmm.github.io/exafmm-t/compile.html#install-exafmm-t) before submitting a PR.
If your contribution introduces new features, please also go through the checklist below:

- add unit test source files in `tests` folder
- add compilation instructions to `tests/Makefile.am` (since we are using autotools)
- new functions and classes should have a Doxygen style docstring

Once your branch is merged into `master`, Travis CI will automatically generate new documentation and push to `gh-pages` branch.
You could follow [these intructions](https://exafmm.github.io/exafmm-t/documentation.html) to build documentation locally.# History

As a bit of history of this project, it started in 2008 with [`PyFMM`](https://github.com/barbagroup/pyfmm), a 2D serial prototype in Python; then followed `PetFMM` in 2009, a PETSc-based parallel code with heavy templating [@cruz2011petfmm]; the third effort was [`GemsFMM`](https://github.com/barbagroup/gemsfmm) in 2010, a serial code with CUDA kernels for execution on GPUs [@yokota2011gems].
Another student in the group felt that [`exafmm-alpha`](https://github.com/exafmm/exafmm-alpha) overused class inheritance and was inadequate for his application, so he and a collaborator re-used only the kernels and implemented another tree construction in [`fmmtl`](https://github.com/ccecka/fmmtl).
The first author of this paper began working with `exafmm-alpha` in 2017, but at this point the mission was clear: simplify and aim for reusability.
That work led to the sixth implementation of FMM in the group, still called [`exafmm`](https://github.com/exafmm/exafmm), but it is not what we present here and we haven't published about it. 
With this new version of ExaFMM (called `exafmm-t`), we aim to bring the FMM to a broader audience and to increased scientific impact.

### Versions

#### exafmm-t
History: 2017/12/11 - Now  
Branch: gpu, vanilla-m2l  
Kernel: LaplaceKI, HelmholtzKI, YukawaKI 
Periodic: no  
SIMD: vec.h  
Thread: OpenMP loops and OpenMP tasks
MPI: none  
GPU: P2P, M2L  
Build: autoconf 
Wrapper: none  
Plotting: none  

#### exafmm
History: 2017/03/03 - Now  
Branch: dev, learning  
Kernel: Laplace, LaplaceKI, Helmholtz, Stokes  
Periodic: yes  
SIMD: vec.h  
Thread: OpenMP tasks  
MPI: HOT (global histogram sort)  
GPU: no  
Build: autoconf  
Wrapper: none  
Plotting: Python  

#### exafmm-alpha > exafmm-beta
History: 2012/07/21 - 2017/03/01  
Branch: develop, sc16  
Kernel: Laplace, Helmholtz, BiotSavart, Van der Waals  
Periodic: yes  
SIMD: Agner's vectorclass  
Thread: OpenMP, Cilk, TBB, Mthreads  
MPI: ORB (bisection, octsection)  
GPU: separate code (Bonsai hack)  
Build: autoconf & CMake  
Wrapper: CHARMM, GROMACS, general MD, PetIGA  
Plotting: none  

#### exafmm-alpha/old + vortex_method > exafmm-alpha
History: 2010/12/22 - 2012/07/21  
Branch: none  
Kernel: Laplace, Van der Waals, Biot Savart, Stretching, Gaussian  
Periodic: yes  
SIMD: none  
Thread: QUARK  
MPI: ORB (global nth_element)  
GPU: offload all kernels  
Build: Makefile  
Wrapper: MR3 compatible MD  
Plotting: VTK  

#### old_fmm_bem
History: mid 2010 - late 2010  
Branch: none  
Kernel: Laplace, Laplace Gn, Helmholtz  
Periodic: no  
SIMD: none  
Thread: none  
MPI: Allgather LET  
GPU: offload all kernels  
Build: Makefile  
Wrapper: none  
Plotting: none  

#### old_fmm_vortex
History: early 2010 - mid 2010  
Branch: none  
Kernel: Laplace, Biot Savart, Stretching (transpose,mixed), Gaussian  
Periodic: yes (shear, channel)  
SIMD: none  
Thread: none  
MPI: Allgather LET  
GPU: offload all kernels  
Build: Makefile  
Wrapper: none  
Plotting: none  

#### Other FMM codes
[ASKIT](http://padas.ices.utexas.edu/libaskit/)  
[Bosnsai](https://github.com/treecode/Bonsai)  
[ChaNGa](https://github.com/N-BodyShop/changa/wiki/ChaNGa)  
[FDPS](https://github.com/FDPS/FDPS)  
[KIFMM](https://cs.nyu.edu/~harper/kifmm3d/documentation/index.html)  
[KIFMM new](https://github.com/jeewhanchoi/kifmm--hybrid--double-only)  
[Modylas](https://github.com/rioyokotalab/modylas)  
[PKIFMM](https://github.com/roynalnaruto/FMM_RPY_BROWNIAN/tree/master/pkifmm)  
[PEPC](http://www.fz-juelich.de/ias/jsc/EN/AboutUs/Organisation/ComputationalScience/Simlabs/slpp/SoftwarePEPC/_node.html)  
[pfalcON](https://pfalcon.lip6.fr)  
[PVFMM](https://github.com/dmalhotra/pvfmm)  
[Salmon treecode](https://github.com/rioyokotalab/salmon_treecode)  
[ScaFaCos](http://www.scafacos.de)  
[ScalFMM](http://people.bordeaux.inria.fr/coulaud/Softwares/scalFMM.html)
###TODO
-------------
- [ ] M2L cache optimization for non-uniform distribution (Tingyu)
- [ ] Revive vanilla m2l branch (Tingyu)
- [ ] M2M, L2L, L2P on GPU (Elket)
- [ ] Compare exafmm vs. exafmm-t (Rio)

###LONG TERM
-------------
- [ ] GPU kernels
- [ ] MPI
- [ ] Stokes
========
Examples
========

C++ Examples
------------

There are three major classes in exafmm-t:

- ``Body<T>``: The class for bodies (particles).
- ``Node<T>``: The class for nodes in the octree.
- ``Fmm``: The FMM class.

The choice of template parameter ``T`` depends on the data type of the potential:
``T`` should be set to ``real_t`` for real-valued kernels (ex. Laplace and modified Helmholtz),
and set to ``complex_t`` for complex-valued kernels (ex. Helmholtz).

exafmm-t uses double precision by default, i.e., ``real_t`` and ``complex_t`` are mapped to ``double`` and ``std::complex<double>`` respectively.
If you want to use single precision, you should still use ``real_t`` and ``complex_t`` in your code,
and add ``-DFLOAT`` to your compiler flags which predefines the macro ``FLOAT`` as true.

All exafmm-t's types, classes and functions are in ``exafmm_t`` namespace.
API documentation can be found in the last section.

Let's solve a Laplace N-body problem as an example, we first need to create ``sources`` and ``targets``.
Here we create 100,000 sources and targets that are randomly distributed in a cube from -1 to 1.
Their type ``Bodies`` is a STL vector of ``Body``.

.. code-block:: cpp
   
   using exafmm_t::real_t;
   std::random_device rd;
   std::mt19937 gen(rd());  // random number generator
   std::uniform_real_distribution<> dist(-1.0, 1.0);
   int ntargets = 100000;
   int nsources = 100000;
   
   exafmm_t::Bodies<real_t> sources(nsources);
   for (int i=0; i<nsources; i++) {
     sources[i].ibody = i;
     sources[i].q = dist(gen);        // charge
     for (int d=0; d<3; d++)
       sources[i].X[d] = dist(gen);   // location
   }

   exafmm_t::Bodies<real_t> targets(ntargets);
   for (int i=0; i<ntargets; i++) {
     targets[i].ibody = i;
     for (int d=0; d<3; d++)
       targets[i].X[d] = dist(gen);   // location
   }

Next, we need to create an FMM instance ``fmm`` for Laplace kernel, and set the order of expansion and ncrit.
We use the former to control the accuracy and the latter to balance the workload between near-field and far-field.

.. code-block:: cpp

   int P = 8;         // expansion order
   int ncrit = 400;   // max number of bodies per leaf
   exafmm_t::LaplaceFmm fmm(P, ncrit);

We can then build and balance the octree. The variable ``nodes`` represents the tree, whose type is ``Nodes``, a STL vector of ``Node``.
To facilitate creating lists and evaluation, we also store a vector of leaf nodes - ``leafs`` and a vector of non-leaf nodes - ``nonleafs``.
Their type ``NodePtrs`` is a STL vector of ``Node*``.

.. code-block:: cpp

   exafmm_t::get_bounds(sources, targets, fmm.x0, fmm.r0);
   exafmm_t::NodePtrs<real_t> leafs, nonleafs;
   exafmm_t::Nodes<real_t> nodes = exafmm_t::build_tree<real_t>(sources, targets, leafs, nonleafs, fmm);
   exafmm_t::balance_tree<real_t>(nodes, sources, targets, leafs, nonleafs, fmm);

Next, we can build lists and pre-compute invariant matrices.

.. code-block:: cpp

   exafmm_t::init_rel_coord();        // compute all possible relative positions of nodes for each FMM operator
   exafmm_t::set_colleagues(nodes);   // find colleague nodes
   exafmm_t::build_list(nodes, fmm);  // create list for each FMM operator
   fmm.M2L_setup(nonleafs);           // an extra setup for M2L operator

Finally, we can use FMM to evaluate potentials and gradients

.. code-block:: cpp

   fmm.upward_pass(nodes, leafs);
   fmm.downward_pass(nodes, leafs);
   
After the downward pass, the calculated potentials and gradients are stored in the leaf nodes of the tree.
You can compute the error in L2 norm by comparing with direct summation:

.. code-block:: cpp

   std::vector<real_t> error = fmm.verify(leafs);
   std::cout << "potential error: " << error[0] << "\n"
             << "gradient error:  " << error[1] << "\n";

Other examples can be found in ``examples/cpp`` folder.

Python Examples
---------------

For simplicity, the name of our Python package is just ``exafmm``.
It has a separate module for each kernel: ``exafmm.laplace``, ``exafmm.helmholtz`` and ``exafmm.modified_helmholtz``.

Compare with C++ interface, exafmm-t's Python interface only exposes high-level APIs.
Now, the steps for tree construction, list construction and pre-computation are merged into one function called ``setup()``.
Also, the evaluation now only requires to call one function ``evalute()``.
Below are Python examples on Jupyter notebooks.

- `Laplace <https://nbviewer.jupyter.org/github/exafmm/exafmm-t/blob/master/examples/python/laplace.ipynb>`__
- `Helmholtz <https://nbviewer.jupyter.org/github/exafmm/exafmm-t/blob/master/examples/python/helmholtz.ipynb>`__
- `Modified Helmholtz <https://nbviewer.jupyter.org/github/exafmm/exafmm-t/blob/master/examples/python/modified_helmholtz.ipynb>`__
===================
Build Documentation
===================

exafmm-t depends on **doxygen**, **sphinx** and **breathe** to generate this documentation. 
To faciliate generating C++ API documentation, we use **exhale**, a Sphinx extension,
to automate launching Doxygen and calling Sphinx to create documentation based on Doxygen xml output.

To build this documentation locally, you need to install doxygen with your package manager, install other dependencies using
``pip install -r docs/requirements.txt``, and then use the following commands:

.. code-block:: bash

   $ cd docs
   $ make html

The HTML documentation will be generated in ``docs/_build/html`` directory.

We also have set up Travis CI to automatically deploy the documentation to Github Pages... exafmm-t documentation master file, created by
   sphinx-quickstart on Mon Apr  8 18:01:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: intro.rst

.. toctree::
   :caption: Contents
   :maxdepth: 3

   compile
   examples
   documentation

.. toctree::
   :caption: API Reference
   :maxdepth: 2

   api/library_root
============
Installation
============

Dependencies
------------
* a C++ compiler that supports OpenMP and C++11 standard (or newer).
* GNU Make
* BLAS
* LAPACK
* `FFTW3 <http://www.fftw.org/download.html>`_
* GFortran

Notes:

* GNU and Intel compilers have been tested. Compilers that use LLVM, such as clang, are not supported yet.
* We recommend to install OpenBLAS. The standard build also includes a full LAPACK library.
* There is no Fortran code in exafmm-t, but gfortran is required in the autotools macros that help configure BLAS and LAPACK libraries.

You can use the following commands to install these dependencies on Ubuntu:

.. code-block:: bash

   $ apt-get update
   $ apt-get -y install libopenblas-dev libfftw3-dev gfortran

Modify these commands accordingly if you are running other Linux distributions.


Install exafmm-t
----------------
This section is only necessary for the users who want to use exafmm-t in C++ applications.
Python users can skip to next section: :ref:`Install exafmm-t's Python package`.

exafmm-t uses **autotools** as the build-system. Go to the root directory of exafmm-t and configure the build:

.. code-block:: bash

   $ cd exafmm-t
   $ ./configure

By default, the configure script will use the most advanced SIMD instruction set 
available on the CPU and enable double precision option. Use ``./configure --help`` to see all available options.

After configuration, you can compile and run exafmm-t's tests with:

.. code-block:: bash

   $ make check

at the root directory of the repo.

Optionally, you can install the headers to the configured location:

.. code-block:: bash

   $ make install


Install exafmm-t's Python package
---------------------------------
exafmm-t relies on `pybind11 <https://github.com/pybind/pybind11>`_ to generate Python bindings.
It requires Python 2.7 or 3.x. To install the Python package, you need first to install OpenBLAS (as the choice of BLAS library),
in addition to the aforementioned dependencies.

Then install exafmm-t to your Python environment using **pip**:

.. code-block:: bash

   $ pip install git+https://github.com/exafmm/exafmm-t.git========
exafmm-t
========

exafmm-t is an open-source fast multipole method (FMM) library to simulate N-body interactions.
It implements the kernel-independent FMM and provides both C++ and Python APIs.
We use `pybind11 <https://github.com/pybind/pybind11>`__ to create Python bindings from the C++ source code.

Exafmm-t currently is a shared-memory implementation using OpenMP.
It aims to deliver competitive performance with a simple code design.
It also has the following features:

- offer high-level APIs in Python
- only use C++ STL containers
- support both single- and double-precision
- vectorization on near-range interactions
- cache optimization on far-range interactions