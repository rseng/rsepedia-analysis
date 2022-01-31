[![Documentation Status](https://readthedocs.org/projects/turtlefsi2/badge/?version=latest)](https://turtlefsi2.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/KVSlab/turtleFSI.svg?branch=master)](https://travis-ci.org/KVSlab/turtleFSI)
[![status](https://joss.theoj.org/papers/b7febdaa2709205d40b51227091c3b0b/status.svg)](https://joss.theoj.org/papers/b7febdaa2709205d40b51227091c3b0b)

# turtleFSI - a Fluid-Structure Interaction Solver

<p align="center">
    <img src="figs/turtleFSI_swim.gif" width="288" height="200" alt="turtleFSI_swim"/>
    <img src="figs/turek_benchmark.gif" width="404" height="200" alt="turtleFSI_swim"/>
</p>
<p align="center">
  To the left we show a turtle swimming (in turtleFSI), and to the right, the classical Turek benchmark (FSI2).
</p>


Description
-----------
turtleFSI is a monolithic fluid-structure interaction solver written in FEniCS, and has out-of-the-box high performance capabilities. The goal of turtleFSI is to provide research groups, and other individuals, with a simple, but robust solver to investigate fluid structure interaction problems.


Authors
-------
turtleFSI is developed by:

  * Andreas Slyngstad
  * Sebastian Gjertsen
  * Aslak W. Bergersen
  * Alban Souche
  * Kristian Valen-Sendstad


Licence
-------
turtleFSI is licensed under the GNU GPL, version 3 or (at your option) any
later version. turtleFSI is Copyright (2016-2019) by the authors.


Documentation
-------------
For an introduction to turtleFSI, and tutorials, please refer to the [documentation](https://turtlefsi2.readthedocs.io/en/latest/).

If you wish to use turtleFSI for journal publications, please refer to the [JOSS publication](https://joss.theoj.org/papers/10.21105/joss.02089#):

Bergersen et al., (2020). turtleFSI: A Robust and Monolithic FEniCS-based Fluid-Structure Interaction Solver. Journal of Open Source Software, 5(50), 2089, https://doi.org/10.21105/joss.02089

Installation
------------
turtleFSI is build upon the open source Finite Elements FEniCS project (version 2018.1.0 or 2019.1.0).
Please refer to the respective FEniCS documentation for installing the dependencies on your system.

However, if you are using Linux or MaxOSX you can install turtleFSI through anaconda::

        conda create -n your_environment -c conda-forge turtleFSI

You can then activate your environment by runing ``source activate your_environment``.
You are now all set, and can start running fluid-structure interaction simulations. If 
you would like to use ``save_deg > 1``, then please install fenicstools by runningthe following::

	pip install git+https://github.com/mikaem/fenicstools

Use
---
Run turtleFSI with all the default parameters::
   ``turtleFSI``

See all the command line parameters run the following command::
  ``turtleFSI -h``

Run a specific problem file::
  ``turtleFSI --problem [path_to_problem]``

When calling a specific problem file, turtleFSI will first look for the file name locally, then check if the file name is present in the directory "/turtleFSI/problems/".
Please refere to the [documentation](https://turtlefsi2.readthedocs.io/en/latest/) to learn how to define a new problem file and for a more complete description of usage.


Contact
-------
The latest version of this software can be obtained from

  https://github.com/KVSlab/turtleFSI

Please report bugs and other issues through the issue tracker at:

  https://github.com/KVSlab/turtleFSI/issues
# Contributing to turtleFSI

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.

### Testing

Please provide unit tests for the new code you create, testing the main functionality or feature to be submitted. We can always use more test coverage!

### Submitting changes

In order to submit you changes, please send a [GitHub Pull Request to turtleFSI](https://github.com/KVSlab/turtleFSI/pull/new/master) with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)). Please follow our coding conventions (below) and make sure all of your commits are atomic (one feature per commit).

Always write a clear commit message when submitting your changes. One-line messages are fine for small changes, but bigger changes should look like this:

    > $ git commit -m "A brief summary of the commit
    >
    > A paragraph describing what changed and its impact."

### Coding conventions

#### Formatting
* Avoid inline comments.
* Break long lines after 120 characters.
* Delete trailing whitespace.
* Don't include spaces after `(`, `[` or before `]`, `)`.
* Don't misspell.
* Use 4 space indentation.
* Use an empty line between methods.
* Use spaces around operators, except for unary operators, such as `!`.
* Use spaces after commas, after colons and semicolons, around `{` and before
  `}`.

#### Naming

* Avoid abbreviations.
* Use snake case for variables and methods.
* Use camel case for classes.
* Name variables, methods, and classes to reveal intent.

#### Organization

* Order methods so that caller methods are earlier in the file than the methods
  they call.
* Place methods receiving command line arguments at the bottom of the file, but above the top-level script environment check.
* Separate local and global imports of modules.

### Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

### Our Standards

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

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at alban@simula.no. All complaints will be reviewed and
investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'turtleFSI: A Robust and Monolithic FEniCS-based Fluid-Structure Interaction Solver'
tags:
  - fluid-structure interaction
  - FEniCS
  - finite-elements
  - numerical methods
authors:
  - name: Aslak W. Bergersen
    orcid: 0000-0001-5063-3680
    affiliation: 1
  - name: Andreas Slyngstad
    affiliation: 1
  - name: Sebastian Gjertsen
    affiliation: 1
  - name: Alban Souche
    orcid: 0000-0001-7547-7979
    affiliation: 1
  - name: Kristian Valen-Sendstad
    orcid: 0000-0002-2907-0171
    affiliation: 1
affiliations:
  - name: Department of Computational Physiology, Simula Research Laboratory, Fornebu, Norway
    index: 1
date: 15 June 2020
bibliography: paper.bib
---

# Summary

It is often sufficient to study fluids [@Moin:1998] and solids [@Holzapfel:2002] in isolation to gain fundamental insights into a physical problem, as other factors may play a secondary role and can be neglected.  On the other hand, there are certain phenomena or situations where the stresses on or by a fluid or a solid can lead to large deformations, and the interaction between fluids and solids are essential [@LeTallec:2001]. Computational fluid-structure interaction (FSI) is an active field of research with much focus on numerical accuracy, stability,  and convergence rates. At the same time, there is also a sweet spot in between these areas of research where there is a need to experiment with FSI without having an in-depth, bottom-up mathematical understanding of the problem, but where a physical insight might suffice. Therefore, the aim was to develop a fully monolithic and robust entry-level research code with ease-of-use targeted towards students, educators, and researchers.

FEniCS [@Logg:2012] has emerged as one of the leading platforms for development of scientific software due to the close connection between mathematical notation and compact computer implementation, where highly efficient C++ code is compiled during execution of a program. Combined with the out-of-the-box entry-level high-performance computing capabilities, FEniCS was a natural choice of computing environment. Compared to other open-source FSI solvers [@Malinen:2013; @Heil:2006; @Jasak:2007], turtleFSI is written in only a couple of hundred lines of high-level Python code, in contrast to tens of thousands of lines of low-level C++ code. This provides full transparency and a unique opportunity for researchers and educators to modify and experiment with the code, while still providing out of the box entry-level high-performance computing capabilities. Furthermore, because of the close resemblance between mathematics and code in FEniCS, users can make additions or modifications with ease.

The turtleFSI solver relies on a fully monolithic approach in the classical arbitrary Lagrangian-Eulerian formulation, and we used the generalized theta scheme for temporal discretization and P2P1P2 elements for velocity, pressure, and displacement, respectively. We implemented and evaluated four different mesh lifting operators, ranging from a simple and efficient second-order Laplace equation, most suitable for small deformations, to more sophisticated and computationally expensive 4th order bi-harmonic equations that can handle larger mesh deformations. We used The Method of Manufactured Solutions to verify the implementation. The obtained results are formally second-order accurate (L2) in space and time [@Wick:2011], respectively, and we demonstrate that all building blocks of code exhibit desired properties. The solver's validity was confirmed using the classical Turek Flag benchmark case [@Turek:2006] with a good agreement – including a diverged numerical solution for long term evolution under certain conditions, as expected. For a complete justification of computational approaches and further details, we refer to [@Slyngstad:2017; @Gjertsen:2017]. We demonstrate adequate strong scaling up to 64 cores (from one cluster node), although the latter is problem size-dependent. In the online documentation, we provide benchmarks, tutorials, and simple demos. The naive FEniCS implementation provides full transparency with compact code, which can easily be adapted to other 2D or 3D FSI problems.

In conclusion, turtleFSI is not a superior FSI solver in terms of speed, but it is a robust entry-level FSI solver and performs exactly as designed and intended; ‘slow and steady wins the race’.


# turtleFSI in Action

turtleFSI comes with several problem files, found under /turtleFSI/problems/, to illustrate the usage and document the Turek flag benchmarks used to validate the implementation of the solver. Here are some illustrations of the execution and outputs expected from the solver.

![Fluid_Turek*='#center'](./cfd_illu.png){ width=100% }\
**Figure 1:**
  Fluid dynamics benchmark snapshot. Simulation executed with the command:
  ```
  turtleFSI --problem TF_cfd
  ```

![Solid_Turek*='#center'](./csm_illu.png){ width=100% }\
**Figure 2:**
  Solid mechanics benchmark snapshots. Simulation executed with the command:
  ```
  turtleFSI --problem TF_csm
  ```

![FSI_Turek*='#center'](./fsi_illu.png){ width=100% }\
**Figure 3:**
  Full fluid-structure interaction benchmark snapshot. Simulation executed with the command:
  ```
  turtleFSI --problem TF_fsi
  ```

# Acknowledgements
The study was supported by The Research Council of Norway through the Center for Biomedical Computing (grant 179578), the Centre for Cardiological Innovation (grant number 203489), and the SIMMIS project (grant number 262827). Simulations were performed on the Abel Cluster (University of Oslo and the Norwegian metacenter for High Performance Computing (NOTUR), project nn9316k), and the Experimental Infrastructure for Exploration of Exascale Computing (eX3) cluster (Norwegian Research Council grant 270053).

# References
.. title:: Acknowledgements and references

.. _acknow_ref:

===============================
Acknowledgements and references
===============================

We would like to acknowledge the open-source project `FEniCS <https://www.fenicsproject.org>`_,
with is the basis of turtleFSI.

The numerical schemes in turtleFSI is presented and tested in Slyngstad [1]_ and Gjertsen [2]_.
The input problem set up for the TF_cfd, TF_csm, and TF_fsi is taken from the Turek et al. [3]_ benchmark
paper.

.. [1] Slyngstad, Andreas S. Verification and Validation of a Monolithic Fluid-Structure Interaction Solver in FEniCS. A comparison of mesh lifting operators. MS thesis. 2017.
.. [2] Gjertsen, Sebastian. Development of a Verified and Validated Computational Framework for Fluid-Structure Interaction: Investigating Lifting Operators and Numerical Stability. MS thesis. 2017.
.. [3] Turek, Stefan, and Jaroslav Hron. "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow." Fluid-structure interaction. Springer, Berlin, Heidelberg, 2006. 371-385.
.. title:: New features

.. _new_features:

============
New features
============

The existing methods provide many degrees of freedom, however, if you need a specific method
or functionality, please do not hesitate to propose enhancements in the
`issue tracker <https://github.com/KVSlab/turtleFSI/issues/>`_, or create a pull request with new features.
Our only request is that you follow our
`guidelines for contributing <https://github.com/KVSlab/turtleFSI/blob/master/CONTRIBUTING.md>`_.
.. title:: Solver verification and performance

.. _verif_perf:

=========================================
Newton solver convergence and performance
=========================================

We illustrate the solver convergence and performance of turtleFSI with a 3D FSI pipe flow problem
(animation below). We impose a constant Dirichlet plug flow at one end of the fluid domain and
compute the resulting structure deformation and fluid flow along the pipe.

.. figure:: ../../figs/movie_36_tstep.gif
   :width: 600px
   :align: center

   Animation of the problem used for benchmarking convergence and performance.

Solver convergence
~~~~~~~~~~~~~~~~~~
The robustness of turtleFSI relies on the use of a direct solver (here, MUMPS) within each Newton
iterations. Using a direct solver for large 3D problems can rapidly become computationally
demanding. Indeed, most of the computational cost of turtleFSI is spent to perform the factorization
of the linear system at each Newton iterations. To mitigate this limitation, we can reuse the same
factorization over several iterations by not updating the Jacobian matrix of the problem. This can
typically be done in simulations where the Newton iterations exhibit “good” converge behavior. In
turtleFSI, the reuse of the Jacobian matrix can be set by the flags ``recompute``, which takes
an integer and controls how many iterations reusing the same Jacobian matrix in the same timestep. 
``recompute_tstep``, does the same, but controls how many time steps to take reusing the same
Jacobian matrix.

.. note::
   Any increase of the relative or absolute residual will trigger turtleFSI to recompute the
   Jacobian matrix, irrespective of the prescribed user values set by the ``recompute`` and
   ``recompute_tstep``.  

Figure 1 illustrates the convergence of the solver using the full Newton procedure, updating the
Jacobian matrix at each iteration step versus reusing the Jacobian matrix over the iterations
(and the time steps). Reusing the Jacobian matrix leads to a larger number of iterations per time
steps, typically ~10 iterations instead of ~5 when updating the Jacobian matrix, but the compute
time is drastically reduced from ca. 20 min. to only 25 s. per time step. The results were produced
with an AMD CPU Ryzen 5 1600, using 12 threads.

.. figure:: ../../figs/reuse_jac_iterations.png
    :width: 600px
    :align: center

    **Figure 1**: Comparison of the convergence behavior of the Newton procedure when updating or
    reusing the Jacobian matrix. The residuals are plotted for the time steps 31 to 36 of the 3D
    FSI pipe flow problem (see above animation). The average execution time for each time step with
    updated Jacobian is 1170 seconds and 25 seconds for the reused Jacobian. 


HPC performance
~~~~~~~~~~~~~~~
turtleFSI benefits from the high-performance computing (HPC) functionality of FEniCS and the solver
can be executed with MPI parallel tasks as follow without any addition to the code::

  mpirun -np 4 turtleFSI

We performed a strong scaling of a 3D FSI pipe flow to illustrate the behavior of the solver using a
relatively large number of parallel MPI tasks. We present the results obtained at the second time step
of the simulation starting from initial rest. We demonstrate an adequate scaling using up to 64 cores
of one cluster node, both executing turtleFSI from a module installation or within a docker container
(Figure 2c). A direct consequence of splitting the geometry in several MPI domains (Figure 2b) is an
increase of the system size associated with the handling of the degree of freedoms along the inner
split boundaries. We illustrate this effect in Figure d where the total memory usage is monitored as
function of the number of MPI tasks used to solve the problem. In our example, we reuse the Jacobian
matrix and the factorization of the direct solver over five Newton’s iterations. As shown in Figure 2c,
the total execution time for computing the Jacobian matrix once and factorization of the system is about
two orders of magnitude larger than solving five iteration steps by reusing the factorization.

.. figure:: ../../figs/figure_hpc2.png
    :width: 600px
    :align: center

    **Figure 2**: Strong scaling of a 3D FSI pipe flow problem. a) Meshing of the inner fluid domain, and
    outer solid pipe with a total of 63 thousand elements. b) Split of the geometry in 16 MPI domains.
    c) Total time spent for one time step (jac.: Jacobian matrix evaluation, fact.: direct solver
    factorization step, it.: direct solver solve steps) as function of the number of MPI tasks. d)
    System memory usage as function of the number of MPI tasks for three different mesh discretizations
    of the problem illustrated in panel a).

.. title:: Installation

.. _installation:

============
Installation
============

Compatibility and Dependencies
==============================
The dependencies of turtleFSI are:

* FEniCS 2019.1.0
* Numpy >1.1X
* Python >=3.5
* FEniCStools 2019.1.0 (optional)

Basic Installation
==================
If you have a MacOX or Linux operating system we recommend that you
install turtleFSI through Anaconda. First, install `Anaconda or Miniconda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda>`_,
depending on your need. For just installing turtleFSI we recommend Miniconda.
Then execute the following command in a terminal window::

    $ conda create -n your_environment -c conda-forge turtleFSI

You can then activate your environment by running ``source activate your_environment``.
Now you are all set, and can start using turtleFSI. A detailed explanation for usage of
turtleFSI can be found `here <https://turtlefsi2.readthedocs.io/en/latest/using_turtleFSI.html>`_.

If you are using turtleFSI on a high performance computing (HPC) cluster we always
recommend that you build from source, as described below. This is in accordance
with the guidelines provided by the `FEniCS project <https://fenicsproject.org/download/>`_
users to install FEniCS from source when on a HPC cluster. To have a fully featured version
of turtleFSI you would also need to install fenicstools manually with::

     $ pip install git+https://github.com/mikaem/fenicstools


Development version
===================

Downloading
~~~~~~~~~~~
The latest development version of turtleFSI can be found on the official
`turtleFSI git repository <https://github.com/KVSlab/turtleFSI>`_ on Github.
To clone the turtleFSI repository, open a terminal, navigate to the directory where you wish
turtleFSI to be stored, type the following command, and press Enter::

    $ git clone https://github.com/KVSlab/turtleFSI

After the source distribution has been downloaded, all the files will be located
in the newly created ``turtleFSI`` folder.

Building
~~~~~~~~
In order to build and install turtleFSI, navigate into the ``turtleFSI`` folder, where a ``setup.py``
file will be located. First, make sure that all dependencies are installed.
Then, you can install turtleFSI be executing the following::

    $ python setup.py install

If you are installing turtleFSI somewhere you do not have root access, typically on a cluster, you can add
``--user`` to install locally. If you are using fenicstools to use the ``save_deg`` functionality, follow the
same procedure, but clone the `fenicstools git repository <https://github.com/mikaem/fenicstools>`_ instead.
.. turtleFSI documentation master file, created by
   sphinx-quickstart on Tue Mar 26 14:24:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: turtleFSI

=====================================
Welcome to turtleFSI's documentation!
=====================================

.. toctree::
   :maxdepth: 3

   installation
   using_turtleFSI
   verif_perf
   new_features
   acknow_ref
.. title:: Using turtleFSI

.. _using_turtleFSI:

===============
Using turtleFSI
===============


Execute a run
=============

turtleFSI is aimed to be user friendly, where all parameters can be controlled from the command line.
We here provide an example on how to create your own problem file. First a quick recap on how to
execute turtleFSI.

To run turtleFSI with all the default parameters you can execute::

 turtleFSI

which will execute the default problem file ``turtle_demo.py`` with all the default parameter. In the folder where
you executed the command there will be a folder ``turtle_demo_results/1`` where you can find ``Visualization`` files,
and under ``Checkpoint`` you can find a list of all parameters used in ``default_variables.pickle``.

.. note::
   The default number of time steps is set to 30, which takes a couple of minutes to execute on a normal desktop.
   To re-create the video shown in readme set the end time (T) to 15 seconds. This will take some time to run,
   but you can look at the intermediate results in ``turtle_demo_results/1`` in paraview while running the simulation.


To run a specific problem file, run::

 turtleFSI --problem [path_to_problem]

To get an overview of all parameters, please run::

 turtleFSI -h


Built-in functionality
======================
TurtleFSI is designed to be lightweight and generic, but we have added to features generally of interest: checkpointing/restart, and
storing files for visualization.

For checkpointing you can set the variable ``--checkpoint-step`` to set the checkpoint frequency, i.e., how often a
checkpoint should be stored. To restart from a previous checkpoint, use the command
``--restart-folder [folder]/[sub-folder]``. Note that the variables from the previous simulation will overwrite any
parameters defined in ``set_problem_parameters`` or on the commandline. If you need to change a parameter from the
previous checkpoint file (for instance, end time ``T``), you can still do it by explicitly redefining the variable
in the ``initiate`` function.

To set how often you save files for visualization you can set ``--save-step``. Note that the default is ``10``.


Setting parameters
==================
All the default parameters are set in the ``problem/__init__.py`` file. Problem specific parameters
are then overwritten in the problem file under ``set_problem_parameters`` or from the command line.
In summary; the parameters will be updated in the following order: default parameters < ``set_problem_parameters`` <
command line (< checkpointing) < other problem specific functions.


Create your own problem file
============================

We have created a step-by-step explanation, see the below, on how you can create your own problem file.

For all numerical problems we have to specify a set parameters, provide a mesh and boundary conditions,
how to solve the equations, and finally specify which metric we are interested in measuring.
In turtleFSI problem file you can define up to seven functions which provides the solver with
the above mentioned information. Listed in the order they are first executed:

- ``set_problem_parameters``
- ``get_mesh_domain_and_boundaries``
- ``initiate`` (optional)
- ``create_bcs``
- ``pre_solve`` (optional)
- ``post_solve`` (optional)
- ``finished`` (optional)


set_problem_parameters
~~~~~~~~~~~~~~~~~~~~~~
This function is for defining parameters of the problem like: time step size ``dt``, end time ``T``, and
physical parameters of the problem. To see a full list of the default parameters, we refer to the
``default_variables`` defined in ``turtleFSI/problems/__init__.py``. Please note that you can, and should, also
define other problem specific variables in ``default_variables``, for instance, geometry information like hight
and length of the problem is of interest to the problem you are solving.

In the ``set_problem_parameters`` function, you have the possibility to set your problem parameters within the
``default_variables`` dictionary. However, keep in mind that any command line arguments will overwrite the ``default_variables``.

A simple example of this function can look like this::

    def set_problem_parameters(default_variables, **namespace):
        # Overwrite default values
        default_variables.update(dict(
            T=15,                          # End time [s]
            dt=0.005,                      # Time step [s]
            theta=0.505,                   # theta value (0.5 + dt), shifted Crank-Nicolson scheme
            Um=1.0,                        # Max. velocity inlet [m/s]
            rho_f=1.0E3,                   # Fluid density [kg/m3]
            mu_f=1.0,                      # Fluid dynamic viscosity [Pa.s]
            rho_s=1.0E3,                   # Solid density [kg/m3]
            mu_s=5.0E4,                    # Solid shear modulus or 2nd Lame Coef. [Pa]
            lambda_s=4.5E5,                # Solid 1st Lame Coef. [Pa]
            nu_s=0.45,                     # Solid Poisson ratio [-]
            dx_f_id=1,                     # ID of marker in the fluid domain
            dx_s_id=2,                     # ID of marker in the solid domain
            extrapolation="biharmonic",    # laplace, elastic, biharmonic, no-extrapolation
            extrapolation_sub_type="constrained_disp",  # ["constant", "small_constant", "volume", "volume_change", "constrained_disp", "constrained_disp_vel"]
            recompute=15,                  # recompute the Jacobian matrix every "recompute" Newton iterations
            folder="turtle_demo_results"), # name of the folder to save the data
            save_step=1                    # frequency of data saving
         )
         return default_variables


.. note::
    The parameter extrapolation here refers to mesh lifting operator, i.e., how to de extrapolate the deformation of
    the solid into the fluid to create an ALE frame of reference. Laplace is the cheapest in terms of computational
    cost, but is less robust for large deformations and sharp edges. In contrast, biharmonic is very robust, but
    at the cost of computational efficiency and increased memory load.

get_mesh_domain_and_boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The function is essential and unique for every problem, and has to be provided for the problem file to run.
Here you read or define your mesh, domain markers, and boundary markers. In ``turtle_demo.py`` we
read the mesh data from pre-existing ".xdmf" mesh files. In contrast to defining the domain and
boundary markers using FEniCS functions, like in the the 'Turek flag'-example (``TF_fsi.py``).

``pygmsh`` and ``meshio`` are relevant tools to create geometries. For any question regarding meshing,
we refer to the FEniCS documentation and `discourse group <https://fenicsproject.discourse.group/>`_.


In ``turtle_demo.py``, the function looks like this::


    def get_mesh_domain_and_boundaries(args, **namespace):
        mesh_folder = path.join(path.dirname(path.abspath(__file__)), "..", "mesh", "turtle_demo")

        # In this example, the mesh and markers are stored in the 3 following files
        mesh_path = path.join(mesh_folder, "turtle_mesh.xdmf")     # mesh geometry
        domains_marker_path = path.join(mesh_folder, "mc.xdmf")    # marker over the elements (domains)
        boundaries_marker_path = path.join(mesh_folder, "mf.xdmf") # markers of the segments (boundaries)

        # "mesh" collects the mesh geometry of the entire domain (fluid + solid).
        # In this example, we import a mesh stored in a .xdmf file, but other formats
        # are supported such as .xml files.
        mesh = Mesh()
        xdmf = XDMFFile(MPI.comm_world, mesh_path)
        xdmf.read(mesh)

        # "domains" collects the element markers of the fluid domain (marked as 1)
        # and the solid domain (marked as 2).
        domains = MeshFunction("size_t", mesh, mesh.geometry().dim())
        xdmf = XDMFFile(MPI.comm_world, domains_marker_path)
        xdmf.read(domains)

        # "boundaries" collects the boundary markers that are used to apply the
        # Dirichlet boundary conditions on both the fluid and solid domains.
        # Marker values ranging from 11 to 15.
        mesh_collection = MeshValueCollection("size_t", mesh, mesh.geometry().dim() - 1)
        xdmf = XDMFFile(MPI.comm_world, boundaries_marker_path)
        xdmf.read(mesh_collection)
        boundaries = cpp.mesh.MeshFunctionSizet(mesh, mesh_collection)

        return mesh, domains, boundaries

.. figure:: ../../figs/Turtle_boundaries.png
   :width: 600px
   :align: center

   Domain boundaries.


initiate
~~~~~~~~
This function is not strictly necessary, but can be used to initiate variables before
entering the time loop of the simulation. Here we have no need for that, and have therefore
not included it. See ``TF_fsi.py`` for an example.


create_bcs
~~~~~~~~~~
The function ``create_bcs`` is used to define the boundary conditions of the problem to be solved,
and is required for the problem file to run. In ``turtle_demo.py``, the inlet boundary condition
is defined the ``Inlet`` class, which inherits the FEniCS ``UserExpression`` class.
This class is then used in the function ``create_bcs`` to prescribe Dirichlet boundary condition to the
inlet velocity. When defining the boundary conditions to specific domain regions or boundaries, make sure to
be consistent with the markers provided in ``get_mesh_domain_and_boundaries``::

    class Inlet(UserExpression):
        def __init__(self, Um, **kwargs):
            self.t = 0.0
            self.t_ramp = 0.5  # time to ramp-up to max inlet velocity (from 0 to Um)
            self.Um = Um       # Max. velocity inlet [m/s]
            super().__init__(**kwargs)

        def update(self, t):
            self.t = t
            if self.t < self.t_ramp:
                self.value = self.Um * np.abs(np.cos(self.t/self.t_ramp*np.pi)-1)/2  # ramp-up the inlet velocity
                print(self.value)
            else:
                Um_min = self.Um/6  # lower velocity during oscillations
                self.value = (self.Um-Um_min) * np.abs(np.cos(self.t/self.t_ramp*np.pi)-1)/2 + Um_min
                print(self.value)

        def eval(self, value, x):
            value[0] = self.value
            value[1] = 0

        def value_shape(self):
            return (2,)


    def create_bcs(DVP, boundaries, Um, v_deg, extrapolation_sub_type, **namespace):
        if MPI.rank(MPI.comm_world) == 0:
            print("Create bcs")

        inlet = Inlet(Um, degree=v_deg)
        noslip = ((0.0, 0.0))

        # Segments indices (make sure of the consistency with the boundary file)
        bottom_id = 11  # segments at the bottom of the model
        outlet_id = 12  # segments at the outlet (right wall) of the model
        top_id = 13     # segments at the top (right wall) of the model
        inlet_id = 14   # segments at the inlet (left wall) of the model
        turtle_head_tail_id = 15   # segments along the head and tail of the turtle

        # Fluid velocity boundary conditions
        u_inlet = DirichletBC(DVP.sub(1), inlet, boundaries, inlet_id)
        u_bot = DirichletBC(DVP.sub(1).sub(1), (0.0), boundaries, bottom_id)  # slip in x-direction
        u_top = DirichletBC(DVP.sub(1).sub(1), (0.0), boundaries, top_id)     # slip in x-direction
        u_head_tail = DirichletBC(DVP.sub(1), noslip, boundaries, turtle_head_tail_id)

        # Pressure boundary conditions
        p_outlet = DirichletBC(DVP.sub(2), (0.0), boundaries, outlet_id)

        # List boundary conditions for the fluid
        bcs = [u_bot, u_top, u_inlet, p_outlet, u_head_tail]

        # Mesh uplifting boundary conditions
        d_inlet = DirichletBC(DVP.sub(0), noslip, boundaries, inlet_id)
        d_bot = DirichletBC(DVP.sub(0), noslip, boundaries, bottom_id)
        d_top = DirichletBC(DVP.sub(0), noslip, boundaries, top_id)
        d_outlet = DirichletBC(DVP.sub(0), noslip, boundaries, outlet_id)
        d_head_tail = DirichletBC(DVP.sub(0), noslip, boundaries, turtle_head_tail_id)

        # Add boundary conditions for the structure
        bcs += [d_bot, d_top, d_outlet, d_inlet, d_head_tail]:

        return dict(bcs=bcs, inlet=inlet)

.. figure:: ../../figs/Turtle_boundaries_zoom.png
    :width: 600px
    :align: center

    Boundaries between the fluid and structures and fixed boundaries.

.. figure:: ../../figs/Turtle_inlet_vel.png
   :width: 600px
   :align: center

   Inlet velocity amplitude variation with time as defined by the class Inlet().



pre_solve
~~~~~~~~~
This function is called within the time loop of the simulation before calling the solver
at the given time step. In ``turtle_demo.py``, we used this function to update the time variable of the
``Inlet`` expression used for the inlet boundary conditions::

    def pre_solve(t, inlet, **namespace):
        # Update the time variable used for the inlet boundary condition
        inlet.update(t)


post_solve
~~~~~~~~~~~
This function is called within the time loop of the simulation after
calling the solver at the given time step. In ``turtle_demo.py``, we do not have any use for
this function, but see ``TF_fsi.py`` for an example.


finished
~~~~~~~~
Function called once at the end of the time loop. An example of use is given in the
``TF_fsi.py`` where text file are saved to store informations from the simulation::

    def finished(folder, dis_x, dis_y, Drag_list, Lift_list, Time_list, **namespace):
        if MPI.rank(MPI.comm_world) == 0:
            np.savetxt(path.join(folder, 'Lift.txt'), Lift_list, delimiter=',')
            np.savetxt(path.join(folder, 'Drag.txt'), Drag_list, delimiter=',')
            np.savetxt(path.join(folder, 'Time.txt'), Time_list, delimiter=',')
            np.savetxt(path.join(folder, 'dis_x.txt'), dis_x, delimiter=',')
            np.savetxt(path.join(folder, 'dis_y.txt'), dis_y, delimiter=',')


Visualizing the result
======================
Given that the parameter ``--save-step`` not was set larger than the number of time steps, there will
be a folder: ``[folder]/[sub-folder]/Visualization`` with ``xdmf`` files that can be opened in a
visualization probrem, for instance ParaView. Below we have visualized the pressure and velocity at 2.5 s.

.. figure:: ../../figs/Turtle_Flow_Pressure_Fields_t_2.5s.png
   :width: 600px
   :align: center

   Pressure and velocity fields at 2.5 s. obtained by running the turtle_demo.py problem file.
