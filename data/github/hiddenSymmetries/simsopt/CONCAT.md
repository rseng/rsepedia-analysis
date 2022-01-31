# simsopt

![GitHub](https://img.shields.io/github/license/hiddensymmetries/simsopt)
[![codecov](https://codecov.io/gh/hiddenSymmetries/simsopt/branch/master/graph/badge.svg?token=ltN6qonZ5p)](https://codecov.io/gh/hiddenSymmetries/simsopt)
[![DOI](https://zenodo.org/badge/247710081.svg)](https://zenodo.org/badge/latestdoi/247710081)

![SIMSOPT](docs/source/logo.png)
![SIMSOPT](docs/source/coils_and_surfaces.png)

`simsopt` is a framework for optimizing
[stellarators](https://en.wikipedia.org/wiki/Stellarator).
The high-level routines of `simsopt` are in python, with calls to C++
or fortran where needed for performance. Several types of components
are included:

- Interfaces to physics codes, e.g. for MHD equilibrium.
- Tools for defining objective functions and parameter spaces for
  optimization.
- Geometric objects that are important for stellarators - surfaces and
  curves - with several available parameterizations.
- Efficient implementations of the Biot-Savart law and other magnetic
  field representations, including derivatives.
- Tools for parallelized finite-difference gradient calculations.

The design of `simsopt` is guided by several principles:

- Thorough unit testing, regression testing, and continuous
  integration.
- Extensibility: It should be possible to add new codes and terms to
  the objective function without editing modules that already work,
  i.e. the [open-closed principle](https://en.wikipedia.org/wiki/Open%E2%80%93closed_principle).
  This is because any edits to working code can potentially introduce bugs.
- Modularity: Physics modules that are not needed for your
  optimization problem do not need to be installed. For instance, to
  optimize SPEC equilibria, the VMEC module need not be installed.
- Flexibility: The components used to define an objective function can
  be re-used for applications other than standard optimization. For
  instance, a `simsopt` objective function is a standard python
  function that can be plotted, passed to optimization packages
  outside of `simsopt`, etc.

`simsopt` is fully open-source, and anyone is welcome to use it, make
suggestions, and contribute.

Several methods are available for installing `simsopt`. One
recommended approach is to use pip:

    pip install simsopt

For detailed installation instructions on some specific systems, see
[the wiki](https://github.com/hiddenSymmetries/simsopt/wiki).
Also, a Docker container is available with `simsopt` and its components pre-installed, which
can be started using

    docker run -it --rm hiddensymmetries/simsopt

More [installation
options](https://simsopt.readthedocs.io/en/latest/getting_started.html#),
[instructions for the Docker
container](https://simsopt.readthedocs.io/en/latest/containers.html), and
other information can be found in the [main simsopt documentation
here.](https://simsopt.readthedocs.io)

Some of the physics modules with compiled code reside in separate
repositories. These separate modules include

- [VMEC](https://github.com/hiddenSymmetries/VMEC2000), for MHD
  equilibrium.
- [SPEC](https://github.com/PrincetonUniversity/SPEC), for MHD
  equilibrium.
- [booz_xform](https://hiddensymmetries.github.io/booz_xform), for
  Boozer coordinates.
  
If you use `simsopt` in your research, kindly cite the code using
[this reference](https://doi.org/10.21105/joss.03525):

[1] M Landreman, B Medasani, F Wechsung, A Giuliani, R Jorge, and C Zhu,
    "SIMSOPT: A flexible framework for stellarator optimization",
    *J. Open Source Software* **6**, 3525 (2021).

See also [the simsopt publications page](https://simsopt.readthedocs.io/en/latest/publications.html).

We gratefully acknowledge funding from the [Simons Foundation's Hidden
symmetries and fusion energy
project](https://hiddensymmetries.princeton.edu). 
# Simsopt testing framework

## Overview

This directory contains integrated/regression tests. Source code for unit tests of each component is stored in the subdirectory for that component.

The layout of the subfolders within **tests** nearly mimics that of the simsopt code in **src/simsopt** folder. The test files (inputs, outputs or any other data files) are all collected into **tests/test_files**.

## Running tests

To run the tests, you must first install simsopt with `pip install .` from the main simsopt directory.
Then, change to the `tests` directory, and run `python -m unittest` to run all tests.
See the [python unittest documentation](https://docs.python.org/3/library/unittest.html) for more options.
# Background on files in this directory

    1DOF_Garabedian.sp

This SPEC input file describes a circular-cross-section torus with
some torsion of the magnetic axis, for a 1-volume vacuum field.

---

    2DOF_targetIotaAndVolume.sp

This SPEC input file describes a classical stellarator with rotating
elliptical cross-section, for a 1-volume vacuum field.

---

    QH-residues.sp

This SPEC input file describes the WISTELL-A configuration in Bader et
al, Journal of Plasma Physics 86, 905860506 (2020). It uses a 1-volume
vacuum field model.

---

    input.circular_tokamak
    boozmn_circular_tokamak.nc

This VMEC input file and BOOZ_XFORM output file are for an
axisymmetric tokamak with circular cross-section, using ITER-like
parameters (but not the ITER poloidal shaping). This example is useful
for testing simsopt functions for the case of axisymmetry.

---

    input.li383_low_res
    wout_li383_low_res_reference.nc
    boozmn_li383_low_res.nc
    
These VMEC input and output files and BOOZ_XFORM output file refer to
the LI383 configuration of NCSX. In contrast to the "official" version
of NCSX, the VMEC resolution parameters `mpol`, `ntor`, and `ns` have
been lowered so the tests run quickly.

---

    input.simsopt_nfp2_QA_20210328-01-020_000_000251
    input.20210406-01-002-nfp4_QH_000_000240

These VMEC input files come from previous simsopt optimizations for
quasi-axisymmetry and quasi-helical symmetry, respectively. These
files are useful for testing measures of quasisymmetry.

---

    input.LandremanSengupta2019_section5.4_B2_A80

This VMEC input file is for the quasi-helically-symmetric
configuration in section 5.4 of Landreman & Sengupta, Journal of
Plasma Physics 85, 815850601 (2019), except that the aspect ratio has
been raised to 80, and the mean field has been increased to 2
Tesla. This configuration is useful for comparing to analytic results
from the near-axis expansion.

---

    input.LandremanSenguptaPlunk_section5p3
    wout_LandremanSenguptaPlunk_section5p3_reference.nc

These VMEC input and output files refer to the configuration in
section 5.3 of Landreman, Sengupta, and Plunk, Journal of Plasma
Physics 85, 905850103 (2019). This configuration is not
stellarator-symmetric, so these files are useful for testing that
simsopt functions work for non-stellarator-symmetric geometries.

---

    input.NuhrenbergZille_1988_QHS

This VMEC input file is for the first quasisymmetric configuration to
be found, from Nuhrenberg and Zille (1988).

---

    input.W7-X_standard_configuration

This VMEC input file is for the W7-X standard configuration with no
plasma pressure or current.

---

    input.cfqs_2b40

This VMEC input file is for a configuration of the CFQS experiment.

---

    input.rotating_ellipse

This VMEC input file is for a classical stellarator with rotating
elliptical cross-section.

---

    tf_only_half_tesla.plasma

This FOCUS input file gives the NCSX plasma boundary, as well as the
component of B normal to this boundary due to a 0.5 Tesla purely
toroidal field.

# Examples

The examples are divided in four folders: **1_Simple**, **2_Intermediate**, **3_Advanced** and **stellarator_benchmarks**. The majority of the examples can be run using the *run_examples* script, which is also called during the continuous integration tests. The files that are generated by running the examples can be easily deleted using the *cleanup* script. The VMEC/SPEC input files needed for each script are inside a subfolder called *inputs*

---

## 1_Simple

Examples where SIMSOPT takes as objective function one/several simple geometric measure such as length, are or volume. There is no need for external libraries or dependencies.

### just_a_quadratic
Minimize f(x,y,z) = ((x-1)/1)^2 + ((y-2)/2)^2 + ((z-3)/3)^2.
### logger_example
Example file for transparently logging both MPI and serial jobs
### minimize_curve_length
Minimize the length of a curve, holding the 0-frequency Fourier mode fixed resulting in a circle.
### surface_volume_and_area
Optimize the minor radius and elongation of an axisymmetric torus to obtain a desired volume and area.
### graph_surf_vol_area
Optimize the minor radius and elongation of an axisymmetric torus to obtain a desired volume and area using the graph optimizable objects.

## 2_Intermediate

Examples where SIMSOPT specifically optimizes for an objective function associated with a stellarator magnetic field such as quasi-symmetry or rotational transform. These scripts need external dependencies, such as VMEC, SPEC or QSC.

### boozer
How to compute surfaces in Boozer coordinates for a magnetic field induced by coils.
### eliminate_magnetic_islands
Show how the shape of a boundary magnetic
surface can be adjusted to eliminate magnetic islands inside it,
considering a vacuum field. The SPEC code is used with a single radial domain.
### QAS
Perform several runs with the VMEC python wrapper while changing a particular surface Fourier coefficient.
### QH_fixed_resolution
Optimize for quasi-helical symmetry (M=1, N=1) at a given radius.
### QSC
Optimize an axis shape and the first-order shape of the flux surface
at first order near the magnetic axis for a target iota and low elongation
using the Stellarator Quasisymmetry Construction code https://github.com/landreman/pyQSC
### resolution_increase
Show how to increase the size of the parameter space and refine the resolution of the calculations during an optimization. The objective function targets quasi-axisymmetry and the iota profile.

## 3_Advanced

Examples where SIMSOPT takes several external libraries together (such as VMEC+SPEC) to optimize for an objective function associated with a stellarator magnetic field, such as quasi-symmetry and magnetic islands.

### optimize_qs_and_islands_simultaneously
simultaneously optimize for quasisymmetry and the elimination of magnetic islands, with both VMEC and SPEC called in the objective function.

## stellarator_benchmarks

This folder contains the scenarios present in https://github.com/landreman/stellopt_scenarios . These are several benchmark problems for stellarator optimization that may be useful for comparing optimization algorithms and for testing new codes. The examples often need either VMEC, SPEC or both.
# Building documentation 

This is a guide to building documentation locally. Online documentation gets
updated whenever master branch is updated.

## Prerequsites

Install *sphinx*, *sphinx-rtd-theme*, *sphinxcontrib-napoleon*,
*sphinx-autodoc-napoleon-typehints* with pip.

## Build
1. Use sphinx-build

```bash
cd docs
sphinx-build -b html source build
```
The documentation is generated in html format and is in the **docs/build**
folder. Start with index.html file.

2. Use the supplied makefile

```bash
cd docs
make html 
```
The documentation is generated in html format and is in the **docs/build/html**
folder. Start with index.html file.

## Update the documentation

Whenever the code is updated, repopulate the code tree and run either step-1  or step-2.

```bash
cd docs
sphinx-apidoc -f -o source ../src/simsopt 
```
# Simsopt and Docker Containers


This document explains how to build docker container for simsopt. It also provides
instructions on how to run simsopt docker container


## Build the container

0. Install docker
1. Build the docker image by running the `docker build` command:

   ```bash
   docker build -t <[user_name/]image_name:[tag]> -f Dockerfile.ubuntu .
   # In the above command, user_name is typically used for uploading the image to docker hub
   # For local builds, you can omit the user_name/ and just use your choice of image_name
   # Tag is optional. If not given, default is latest
   ```
2. Upload the image to docker hub if you want to distribute it
   ```bash
   docker push <image_name:tag>
   ```
   

## Run the container
There are two cases here: 1) docker image is built locally and 2) docker image downloaded from a repo

### For local builds:
0. Identify the local image using the command
   ```bash
    docker image ls
    ```
    After identifying the image built locally, you can use either <image_name:tag> or the image id, which 
    is an hexadecimal number, to run the docker image in the next command

### Remote repo
0. Assuming image is uploaded to docker hub, use the <user_name/image_name:tag> to automatically download
   the image in the next command

1. For interactive runs, run the docker container with the command:
   ```bash
   docker run -it --rm  <image>
   ```

2. To execute a python script
   ```bash
   docker run -v$PWD:$PWD <image> python3 $PWD/<python_script>
   ```
   **This command works only if no additional files are needed to run the script**



simsopt.\_core package
======================

Submodules
----------

simsopt.\_core.derivative module
--------------------------------

.. automodule:: simsopt._core.derivative
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

simsopt.\_core.dofs module
--------------------------

.. automodule:: simsopt._core.dofs
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

simsopt.\_core.graph\_optimizable module
----------------------------------------

.. automodule:: simsopt._core.graph_optimizable
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

simsopt.\_core.optimizable module
---------------------------------

.. automodule:: simsopt._core.optimizable
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

simsopt.\_core.util module
--------------------------

.. automodule:: simsopt._core.util
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

Module contents
---------------

.. automodule:: simsopt._core
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
Contributing to Simsopt
=======================

We are glad you decided to contribute to ``simsopt``! During development please
follow the contribution guidlines posted here. 


Types of Contribution
^^^^^^^^^^^^^^^^^^^^^

Both big and small contributions to ``simsopt`` are welcome. Some ways you can contribute to 
``simsopt`` are:

- Submit feedback
- Report bugs
- Fix bugs
- Improve documentation
- Request new features
- Contribute new algorithms
- Implement new features

Submit Feedback
---------------

You can give your feedback on ``simsopt``  by opening an `issue <https://github.com/hiddensymmetries/simsopt/issues>`_.

If you are proposing/requesting a new feature or new algorithm:

- Explain in detail why the feature is needed with a demo problem that couldn't be implemented with existing code.
- Features that are ambitious will take time or may not be implemented at all. So, keep the scope of the feature as narrow as possible, to make it easier to implement.
- Your contributions are always welcome!


Bug Reports
-----------

To report a bug in the package, open an `issue <https://github.com/hiddensymmetries/simsopt/issues>`_.

Please include in your bug report:

* Your operating system type (mac or linux) and version.
* Pertinent details about your local setup (such as MPI and compiler info) that might be helpful in troubleshooting.
* Steps to reproduce the bug.

Fix Bugs
--------

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement New Features or Algorithms
------------------------------------

First become familiar with simsopt code (at the least, the subpackage you are contributing to).
Quickest way is to reach out to the developers and ask for help. Simsopt enforces some code quality
throught unit-tests. Have the unit tests ready to go with your code. Having few examples showcasing
a problem or two solved with the new feature would be even better.

Improve Documentation
---------------------

If you feel the documentation is lagging at any place, please feel
free to submit a PR focused on fixing or improving the 
documentation. Steps to build documentation locally can be found `here <https://github.com/hiddenSymmetries/simsopt/tree/contributing/docs>`_.


Code development Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


1. Log in to GitHub and fork ``simsopt`` repository. 
   This will create a duplicate at 'github.com/<your_account>/simsopt' 
   which is your personal fork (or just 'your fork'), as opposed to this repository
   (github.com/hiddensymmetries/simsopt), which is conventionally referred to as upstream repository in git parlance.

2. Clone your fork to your machine by opening a console and executing

   .. code-block::

        git clone https://github.com/<your_account>/simsopt.git

   Make sure to clone your fork, not the upstream repo. This will create a
   directory called ``simsopt``. Navigate to it and execute

   .. code-block::

        git remote add upstream https://github.com/hiddensymmetries/simsopt.git

   In this way, your machine will know of both your fork (which git calls
   `origin`) and the upstream repository (`upstream`).

3. During development, you make changes to the code in your fork.
   code. To prevent frequent reinstallation of simsopt after each modification, 
   and to reflect the changes immediately, install ``simsopt`` as editable.

   .. code-block::
	
        pip install -e .

4. Using a new branch to start implementing your changes would be a good idea.

   .. code-block::

        git checkout -b <your_branch_name>

5. Use git add and commit commands to save your changes.
    
   .. code-block::

        git add <your_new_file_or_modified_file>
        git commit -m "Brief message highlighting the changes implemented"

6. Make sure run_tests, run_tests_mpi, and examples/run_examples all pass. Running these locally will help you to catch bugs while developing.

7. Before submitting your changes, run ``autopep8`` to fix formatting issues using the supplied ``run_autopep`` script.
   Don't forget to run step 5, once again.

9. Push the changes to github. 

    .. code block::
        git push

10. Once the changes are in your fork you can submit a pull-request. PRs will only be merged if run_tests, 
    run_tests_mpi, and examples/run_examples all pass. Request at least 1 review from the ``simsopt`` 
    developers (Bharat Medasani, Matt Landreman, Florian Wechshung). The last reviewer to approve will be in charge of merging.
    Your contributions will be reviewed and merged to upstream repository if ``simsopt`` developers are 
    satisfied with the code quality. Please note this is not a full tutorial on using git and you need to know additional
    git commands to be more efficient or not to get stuck with git conflicts.
simsopt.geo package
===================

Submodules
----------

simsopt.geo.boozersurface module
--------------------------------

.. automodule:: simsopt.geo.boozersurface
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.config module
-------------------------

.. automodule:: simsopt.geo.config
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curve module
------------------------

.. automodule:: simsopt.geo.curve
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curvehelical module
-------------------------------

.. automodule:: simsopt.geo.curvehelical
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curveperturbed module
---------------------------------

.. automodule:: simsopt.geo.curveperturbed
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curveobjectives module
----------------------------------

.. automodule:: simsopt.geo.curveobjectives
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curverzfourier module
---------------------------------

.. automodule:: simsopt.geo.curverzfourier
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.curvexyzfourier module
----------------------------------

.. automodule:: simsopt.geo.curvexyzfourier
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.jit module
----------------------

.. automodule:: simsopt.geo.jit
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.plot module
-----------------------

.. automodule:: simsopt.geo.plot
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.qfmsurface module
-----------------------------

.. automodule:: simsopt.geo.qfmsurface
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surface module
--------------------------

.. automodule:: simsopt.geo.surface
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfacegarabedian module
------------------------------------

.. automodule:: simsopt.geo.surfacegarabedian
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfacehenneberg module
-----------------------------------

.. automodule:: simsopt.geo.surfacehenneberg
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfaceobjectives module
------------------------------------

.. automodule:: simsopt.geo.surfaceobjectives
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfacerzfourier module
-----------------------------------

.. automodule:: simsopt.geo.surfacerzfourier
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfacexyzfourier module
------------------------------------

.. automodule:: simsopt.geo.surfacexyzfourier
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.geo.surfacexyztensorfourier module
------------------------------------------

.. automodule:: simsopt.geo.surfacexyztensorfourier
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.geo
   :members:
   :undoc-members:
   :show-inheritance:
Field line and particle tracing
===============================

The :obj:`simsopt.field.tracing` module provides tools for tracing
field lines, and for following particle motion in magnetic fields.
For the latter, full orbits (including gyromotion) or guiding center
trajectories can be followed in cylindrical coordinates, or guiding
center motion can be followed in Boozer coordinates.  Examples of
these various tracing features can be found in
``examples/1_Simple/tracing_fieldline.py``,
``examples/1_Simple/tracing_particle.py``, and
``examples/2_Intermediate/tracing_boozer.py``,



Particle motion in cylindrical coordinates
------------------------------------------

The main function to use in this case is
:obj:`simsopt.field.tracing.trace_particles` (click the link for more
information on the input and output parameters) and it is able to use
two different sets of equations depending on the input parameter
``mode``:

- In the case of ``mode='full'`` it solves

.. math::

  [\ddot x, \ddot y, \ddot z] = \frac{q}{m}  [\dot x, \dot y, \dot z] \times \mathbf B

- In the case of ``mode='gc_vac'`` it solves the guiding center equations under the assumption :math:`\nabla p=0`, that is

.. math::

  [\dot x, \dot y, \dot z] &= v_{||}\frac{\mathbf B}{B} + \frac{m}{q|B|^3}  \left(\frac{v_\perp^2}{2} + v_{||}^2\right)  \mathbf B\times \nabla B\\
  \dot v_{||}    &= -\mu  \mathbf B \cdot \nabla B

where :math:`v_\perp^2 = 2\mu B`.
See equations (12) and (13) of `Guiding Center Motion, H.J. de Blank <https://doi.org/10.13182/FST04-A468>`_.

Below is an example of the vertical drift experienced by two particles in a simple toroidal magnetic field.

.. code-block::

    from simsopt.field.magneticfieldclasses import ToroidalField
    from simsopt.geo.curvexyzfourier import CurveXYZFourier
    from simsopt.util.constants import PROTON_MASS, ELEMENTARY_CHARGE, ONE_EV
    from simsopt.field.tracing import trace_particles_starting_on_curve

    bfield = ToroidalField(B0=1., R0=1.)
    start_curve = CurveXYZFourier(300, 1)
    start_curve.set_dofs([0, 0, 1.01, 0, 1.01, 0., 0, 0., 0.])
    nparticles = 2
    m = PROTON_MASS
    q = ELEMENTARY_CHARGE
    tmax = 1e-6
    Ekin = 100*ONE_EV
    gc_tys, gc_phi_hits = trace_particles_starting_on_curve(
        start_curve, bfield, nparticles, tmax=tmax, seed=1, mass=m, charge=q,
        Ekin=Ekin, umin=0.2, umax=0.5, phis=[], mode='gc_vac', tol=1e-11)
    z_particle_1 = gc_tys[0][:][2]
    z_particle_2 = gc_tys[1][:][2]
    print(z_particle_1)
    print(z_particle_2)


We note that SIMSOPT draws initial data for particles consisting of
the guiding center position, the parallel velocity, and the
perpendicular speed.

* To compute the speed, the user specifies the total kinetic energy and an interval of pitch-angles :math:`[u_\min, u_\max]`. Given a pitch angle :math:`u` and total speed :math:`v` (computed from the kinetic energy and the particle mass), the parallel speed is given by :math:`v_{||} = u v` and, and :math:`v_\perp = \sqrt{v^2-v_{||}^2}`.
* In the case of full orbit simulations, we need velocity initial data and hence this only defines the initial data up to the phase. To specify the phase, one can pass the ``phase_angle`` variable to the tracing functions. A value in :math:`[0, 2\pi]` is expected.
Optimizing for quasisymmetry
============================

In this tutorial it is shown how to optimize the boundary shape of a
VMEC configuration to achieve quasisymmetry.  Several variants of this
problem are shown here, including both quasi-helical symmetry and
quasi-axisymmetry.  There are two objective functions available in
simsopt to achieve quasisymmetry, and both are demonstrated here.  One
objective function is
:obj:`simsopt.mhd.vmec_diagnostics.QuasisymmetryRatioResidual`, based
on the fact that :math:`(\vec{B}\times\nabla B
\cdot\nabla\psi)/(\vec{B}\cdot\nabla B)` is a flux function in
quasisymmetry, demonstrated in the paper `arXiv:2108.03711
<https://arxiv.org/pdf/2108.03711>`__.  The other objective function,
using :obj:`simsopt.mhd.boozer.Quasisymmetry`, is a more traditional
one based on the Fourier modes in Boozer coordinates :math:`B_{m,n}`
that break the symmetry. The recommended approach is to use
:obj:`~simsopt.mhd.vmec_diagnostics.QuasisymmetryRatioResidual`, since
it has the modest advantage of not requiring a transformation to
Boozer coordinates, although both methods are effective. The first two
examples on this page do not require the ``booz_xform`` package,
whereas the third example does.

For these optimization examples, you probably want to write a script
that is submitted to the queue on a computing cluster, where they are
likely to take minutes or tens of minutes on one node.



Fixed resolution
----------------

..
   This example was run on IPP-Cobra in /ptmp/mlan/20211217-01-simsopt_docs_tutorials/20211217-01-001_QH_fixed_resolution
   The final configuration is also available at
   ~/Box Sync/work21/wout_20211217-01-001_simsopt_docs_tutorials_nfp4_QH_warm_start_000_000038.nc

This example is also available as ``QH_fixed_resolution.py`` in the
``examples/2_Intermediate`` directory.  As usual, a driver script begins with
imports of the classes and functions we will need::

  import numpy as np
  from simsopt.util.mpi import MpiPartition
  from simsopt.mhd.vmec import Vmec
  from simsopt.mhd.vmec_diagnostics import QuasisymmetryRatioResidual
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_mpi import least_squares_mpi_solve

For this problem we will want MPI for parallelized finite difference
gradients. As explained below, this particular problem has 24
independent variables, so we can take advantage of 24 + 1 = 25 concurrent
function evaluations for one-sided differences. It is therefore
optimal to divide the MPI processes into 25 worker groups. For more
information about MPI and worker groups, see :doc:`mpi` or
:obj:`~simsopt.util.mpi.MpiPartition`.  In our script, we therefore
use::

  mpi = MpiPartition(25)

It is not necessary to use 25 groups - the code will function properly
for any choice.  If you select more than 25 groups, the groups beyond
the 25th will sit idle. If you select fewer than 25 groups, the 25
function evaluations needed for a finite difference gradient will be
distributed among the groups that are available.  There is no need to
make the number of groups a multiple of the number of available MPI
processes, although there cannot be more groups than processes.

We next initialize a VMEC object from an input file::

  vmec = Vmec("input.nfp4_QH_warm_start", mpi=mpi)

This file can be found in the ``examples/2_Intermediate/inputs``
directory. The file describes a configuration that has already been
partially optimized for quasi-helical symmetry in a very small
parameter space, keeping poloidal and toroidal mode numbers (the
latter divided by the number of field periods) only up through 1:

.. image:: example_quasisymmetry_QH_before.png
   :width: 400

In the present example we will refine the configuration by optimizing
in a larger parameter space, with poloidal and toroidal mode numbers
up through 2. To define this parameter space, we set the ``fixed``
property of the boundary's Fourier modes as follows::

  # Define parameter space:
  surf = vmec.boundary
  surf.fix_all()
  max_mode = 2
  surf.fixed_range(mmin=0, mmax=max_mode,
                   nmin=-max_mode, nmax=max_mode, fixed=False)
  surf.fix("rc(0,0)") # Major radius

The above code first fixes all modes of the boundary, since we want
the mode numbers greater than 2 to all be fixed. Then the desired
range of modes is set to be not fixed. This range includes the m=n=0
mode which is essentially the mean major radius. We don't need or want
to vary the overall scale of the configuration, so it is convenient to
remove this mode from the parameter space. To confirm which modes will
be varied, you can print out the names of the free degrees of freedom
(dofs)::

  print('Parameter space:', surf.dof_names)

The result is

.. code-block::

   Parameter space: ['SurfaceRZFourier1:rc(0,1)', 'SurfaceRZFourier1:rc(0,2)',
   'SurfaceRZFourier1:rc(1,-2)', 'SurfaceRZFourier1:rc(1,-1)',
   'SurfaceRZFourier1:rc(1,0)', 'SurfaceRZFourier1:rc(1,1)',
   'SurfaceRZFourier1:rc(1,2)', 'SurfaceRZFourier1:rc(2,-2)',
   'SurfaceRZFourier1:rc(2,-1)', 'SurfaceRZFourier1:rc(2,0)',
   'SurfaceRZFourier1:rc(2,1)', 'SurfaceRZFourier1:rc(2,2)',
   'SurfaceRZFourier1:zs(0,1)', 'SurfaceRZFourier1:zs(0,2)',
   'SurfaceRZFourier1:zs(1,-2)', 'SurfaceRZFourier1:zs(1,-1)',
   'SurfaceRZFourier1:zs(1,0)', 'SurfaceRZFourier1:zs(1,1)',
   'SurfaceRZFourier1:zs(1,2)', 'SurfaceRZFourier1:zs(2,-2)',
   'SurfaceRZFourier1:zs(2,-1)', 'SurfaceRZFourier1:zs(2,0)',
   'SurfaceRZFourier1:zs(2,1)', 'SurfaceRZFourier1:zs(2,2)']

Next, we need to configure a term in the objective function to
represent the departure from quasisymmetry. This can be done as
follows::

  # Configure quasisymmetry objective:
  qs = QuasisymmetryRatioResidual(vmec,
                                  np.arange(0, 1.01, 0.1),  # Radii to target
				  helicity_m=1, helicity_n=-1)  # (M, N) you want in |B|

There are several adjustable options, the details of which can be
found in the API documentation for
:obj:`~simsopt.mhd.vmec_diagnostics.QuasisymmetryRatioResidual`.
There you can also find the mathematical expression for the objective
function.  The second argument to
:obj:`~simsopt.mhd.vmec_diagnostics.QuasisymmetryRatioResidual` above
sets the quasisymmetry objective to be evaluated at a uniform grid of
11 surfaces ``[0, 0.1, 0.2, ..., 1]`` in the normalized toroidal flux
:math:`s`, with the result that quasisymmetry is targeted throughout
the volume.  You are free to provide different values, or a single
float if you only want to target quasisymmetry on a single
surface. There is also an optional argument ``weights`` if you wish to
have different weights in the objective function for quasisymmetry on
different surfaces. The ``helicity_n`` argument can also be set to
``+1`` rather than ``-1`` for quasi-helical symmetry, amounting to a
mirror-reversal, though the initial configuration used here is
consistent with the ``-1`` handedness.

We are now ready to define the total objective function. Here we will
include quasisymmetry and aspect ratio. Aspect ratio must be included
because otherwise quasisymmetry can be made arbitrarily good by
increasing the aspect ratio to infinity. The simsopt objective
function is defined as follows::

  # Define objective function
  prob = LeastSquaresProblem.from_tuples([(vmec.aspect, 7, 1),
                                          (qs.residuals, 0, 1)])

It can be seen that we are targeting an aspect ratio of 7. This
objective function will be a sum of 44,353 least-squares terms, 44,352
of which correspond to the quasisymmetry residual on 63x64 grid points
on the 11 flux surfaces targeted, plus one additional term
``(vmec.aspect - 7) ** 2``. (The 63x64 resolution is a default in
:obj:`~simsopt.mhd.vmec_diagnostics.QuasisymmetryRatioResidual`.)  This
large number of residual terms is no problem - it introduces
negligible computational cost compared to the cost of the equilibrium
calculations, so we may as well use this high resolution.

You can check the value of the objective functions before the
optimization. Rather than print each residual term, the scalar total
for the quasisymmetry term can be obtained with the ``.total()``
method.

.. code-block::

   print("Quasisymmetry objective before optimization:", qs.total())
   print("Total objective before optimization:", prob.objective())

The results are both 0.304, since the aspect ratio term is negligible.

Finally, we solve the optimization problem::

  least_squares_mpi_solve(prob, mpi, grad=True)

Suppose you have written the above commands in a file named
``simsopt_driver``.  Depending on your computing system, the script
can be run using a command like ``srun python simsopt_driver`` (for
SLURM systems) or ``mpirun -n 25 simsopt_driver``.

Since this objective function has multiple local minima, the final
result of the optimization can be sensitive to small changes in
simsopt, VMEC, or the packages they depend on. Therefore you will not
necessarily obtain exactly the result shown here. But one result
produced by this optimization script is the following configuration:

.. image:: example_quasisymmetry_QH_after.png
   :width: 400
.. image:: example_quasisymmetry_QH_after_3D.png
   :width: 400
..
   Figure produced by ~/Box Sync/MATLAB/m20210207_plotVMECWout.m
.. image:: example_quasisymmetry_QH_after_Boozer.png
   :width: 400
..
   Figure produced by ~/Box Sync/work21/boozPlotHalfFluxUnfilled wout_20211217-01-001_simsopt_docs_tutorials_nfp4_QH_warm_start_000_000038.nc

This last figure shows that reasonably good quasisymmetry has been
achieved. The quality of quasisymmetry can be improved significantly
by further refining the configuration using one or more rounds of
optimization with more Fourier modes in the parameter space. Printing
``qs.total()`` or ``prob.objective()`` at the end of the optimization,
it can be seen that both have been reduced significantly, to 0.00794
for the result shown here.


Dynamic resolution
------------------
..
   This example was run on IPP-Cobra in /ptmp/mlan/20211217-01-simsopt_docs_tutorials/20211217-01-003_QA_dynamic_resolution
   The final configuration is also available at
   ~/Box Sync/work21/wout_20211217-01-003_simsopt_docs_tutorials_QA_dynamic_resolution_000_000205.nc

Since simsopt optimization problems are defined using a python script,
you are free to add other scripting in your problem definition. Here
we show how this capability can be used to increase the numerical
resolution of codes such as VMEC during the optimization. At the same
time, we will increase the number of Fourier modes in the parameter
space during the optimization. This example can also be found in the
``examples/2_Intermediate`` directory as
``resolution_increase.py``. This example is very similar to the
quasi-axisymmetry optimization shown in `arXiv:2108.03711
<https://arxiv.org/pdf/2108.03711>`__.

As usual, we begin with the necessary imports::

  import numpy as np
  from simsopt.util.mpi import MpiPartition
  from simsopt.mhd.vmec import Vmec
  from simsopt.mhd.vmec_diagnostics import QuasisymmetryRatioResidual
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_mpi import least_squares_mpi_solve

We again split the pool of MPI processes into worker groups. Here, for
simplicity, we make each process its own worker group, by omitting the
argument::

  mpi = MpiPartition()

We initialize a VMEC configuration from an input file. This starting
configuration is axisymmetric with a circular cross-section, so we are
starting "from scratch"::

  vmec = Vmec("input.nfp2_QA", mpi=mpi)

This input file can be found in the ``examples/2_Intermediate/inputs``
directory. We define the quasisymmetry objective as in the previous
section, except that we specify a helicity of (1,0) instead of (1,1)
or (1,-1) to get quasi-axisymmetry instead of quasi-helical symmetry::

  # Configure quasisymmetry objective:
  qs = QuasisymmetryRatioResidual(vmec,
                                  np.arange(0, 1.01, 0.1),  # Radii to target
				  helicity_m=1, helicity_n=0)  # (M, N) you want in |B|
				  
We now define the total objective function. For this example, it is
necessary to include a nonzero target value for the rotational
transform in the objective, to prevent the optimum from being truly
axisymmetric::

  # Define objective function
  prob = LeastSquaresProblem.from_tuples([(vmec.aspect, 6, 1),
                                          (vmec.mean_iota, 0.42, 1),
                                          (qs, 0, 1)])

It can be seen here that we are seeking a configuration with aspect
ratio 6, and average iota of 0.42, slightly above the resonance at 2 /
5 = 0.4. The function :func:`simsopt.mhd.vmec.Vmec.mean_iota()` used
here returns :math:`\int_0^1 \iota\, ds` where :math:`s` is the
toroidal flux normalized by its value at the VMEC boundary.

Now, we set up a loop over several optimization steps. At each step,
the resolution parameters ``mpol`` and ``ntor`` for VMEC increase. At
the same time, in each optimization step a larger range of poloidal
and toroidal mode numbers are set to be varied in the optimization::

  for step in range(4):
      max_mode = step + 1
    
      # VMEC's mpol & ntor will be 3, 4, 5, 6:
      vmec.indata.mpol = 3 + step
      vmec.indata.ntor = vmec.indata.mpol
    
      if mpi.proc0_world:
          print("Beginning optimization with max_mode =", max_mode, \
                ", vmec mpol=ntor=", vmec.indata.mpol, \
                ". Previous vmec iteration = ", vmec.iter)

      # Define parameter space:
      surf.fix_all()
      surf.fixed_range(mmin=0, mmax=max_mode, 
                       nmin=-max_mode, nmax=max_mode, fixed=False)
      surf.fix("rc(0,0)") # Major radius

      # Carry out the optimization for this step:
      least_squares_mpi_solve(prob, mpi, grad=True)

      if mpi.proc0_world:
          print("Done optimization with max_mode =", max_mode, \
                ". Final vmec iteration = ", vmec.iter)

If you like, other parameters could be adjusted at each step too, such
as the radial resolution or number of iterations in VMEC, the solver
tolerances, or the maximum number of iteration of the optimization
algorithm.

As in the previous section, the final result of this optimization can
be sensitive to small changes in simsopt, VMEC, or the packages they
depend on. Therefore you will not necessarily obtain exactly the
result shown here. But one result produced by this optimization script
is the following configuration:

.. image:: example_quasisymmetry_QA_after.png
   :width: 400
.. image:: example_quasisymmetry_QA_after_3D.png
   :width: 400
..
   Figure produced by ~/Box Sync/MATLAB/m20210207_plotVMECWout.m
.. image:: example_quasisymmetry_QA_after_Boozer.png
   :width: 400
..
   Figure produced by ~/Box Sync/work21/boozPlotHalfFluxUnfilled wout_20211217-01-003_simsopt_docs_tutorials_QA_dynamic_resolution_000_000205.nc


Bmn objective
-------------

Here we show an alternative method of quasisymmetry optimization using
a different objective function,
:obj:`simsopt.mhd.boozer.Quasisymmetry`, based on the
symmetry-breaking Fourier mode aplitudes :math:`B_{m,n}` in Boozer
coordinates.  This example can also be found in the
``examples/2_Intermediate`` directory as
``resolution_increase_boozer.py``.

In this case, the imports needed are::

  from simsopt.util.mpi import MpiPartition
  from simsopt.mhd.vmec import Vmec
  from simsopt.mhd.boozer import Boozer, Quasisymmetry
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_mpi import least_squares_mpi_solve

We again split the pool of MPI processes into worker groups and
initialize a ``Vmec`` object as in the previous example::

  mpi = MpiPartition()
  vmec = Vmec("input.nfp2_QA", mpi=mpi)

This input file, corresponding to an axisymmetric torus with circular
cross-section, can be found in the ``examples/2_Intermediate/inputs``
directory. Next, this alternative quasisymmetry objective can be
created as follows::

  # Configure quasisymmetry objective:
  boozer = Boozer(vmec)
  qs = Quasisymmetry(boozer,
                     0.5, # Radius to target
                     1, 0) # (M, N) you want in |B|

There are several adjustable options, the details of which can be
found in the API documentation for :obj:`~simsopt.mhd.boozer.Boozer`
and :obj:`~simsopt.mhd.boozer.Quasisymmetry`. The numerical resolution
of the Boozer-coordinate transformation can be adjusted by passing
parameters to the :obj:`~simsopt.mhd.boozer.Boozer` constructor, as in
``Boozer(vmec, mpol=64, ntor=32)``. The second argument to
``Quasisymmetry`` above sets the quasisymmetry objective to be
evaluated at normalized toroidal flux of 0.5, but you are free to
provide different values.  Or, a list of values can be provided to
target quasisymmetry on multiple surfaces. The
:obj:`~simsopt.mhd.boozer.Quasisymmetry` also has optional arguments
to adjust the normalization and weighting of different Fourier modes.

We now define the total objective function. As with the previous
quasi-axisymmetry example, it is necessary to include a nonzero target
value for the rotational transform in the objective, to prevent the
optimum from being truly axisymmetric. Here we will constrain iota
at the edge and magnetic axis, in order to prescribe the magnetic shear::

  # Define objective function
  prob = LeastSquaresProblem.from_tuples([(vmec.aspect, 6, 1),
                                          (vmec.iota_axis, 0.465, 1),
                                          (vmec.iota_edge, 0.495, 1),
                                          (qs, 0, 1)])

It can be seen here that we are seeking a configuration with aspect
ratio 6, and iota slightly below 0.5.

Now, we set up a loop over several optimization steps. At each step,
the resolution parameters ``mpol`` and ``ntor`` for VMEC increase, as
do the the Fourier resolution parameters for ``booz_xform``. At the
same time, in each optimization step a larger range of poloidal and
toroidal mode numbers are set to be varied in the optimization::

  for step in range(4):
      max_mode = step + 1
    
      # VMEC's mpol & ntor will be 3, 4, 5, 6:
      vmec.indata.mpol = 3 + step
      vmec.indata.ntor = vmec.indata.mpol

      # booz_xform's mpol & ntor will be 16, 24, 32, 40:
      boozer.mpol = 16 + step * 8
      boozer.ntor = boozer.mpol
    
      if mpi.proc0_world:
          print("Beginning optimization with max_mode =", max_mode, \
                ", vmec mpol=ntor=", vmec.indata.mpol, \
                ", boozer mpol=ntor=", boozer.mpol, \
                ". Previous vmec iteration = ", vmec.iter)

      # Define parameter space:
      surf.fix_all()
      surf.fixed_range(mmin=0, mmax=max_mode, 
                       nmin=-max_mode, nmax=max_mode, fixed=False)
      surf.fix("rc(0,0)") # Major radius

      # Carry out the optimization for this step:
      least_squares_mpi_solve(prob, mpi, grad=True)

      if mpi.proc0_world:
          print("Done optimization with max_mode =", max_mode, \
                ". Final vmec iteration = ", vmec.iter)

If you like, other parameters could be adjusted at each step too, such
as the radial resolution or number of iterations in VMEC, the solver
tolerances, or the maximum number of iteration of the optimization
algorithm.

As with the previous examples, the final result of this optimization
can be sensitive to small changes in simsopt, VMEC, or the packages
they depend on. Therefore you will not necessarily obtain exactly the
result shown here. But one result produced by this optimization script
is the following configuration:

.. image:: example_quasisymmetry_QA_Bmn_after.png
   :width: 400
.. image:: example_quasisymmetry_QA_Bmn_after_3D.png
   :width: 400
..
   Figure produced by ~/Box Sync/MATLAB/m20210207_plotVMECWout.m
.. image:: example_quasisymmetry_QA_Bmn_after_Boozer.png
   :width: 400
..
   Figure produced by ~/Box Sync/work21/boozPlotHalfFluxUnfilled simsopt_nfp2_QA_20210328-01-020_000_000251/wout_simsopt_nfp2_QA_20210328-01-020_000_000251_scaled.nc
simsopt.objectives package
==========================

Submodules
----------

simsopt.objectives.fluxobjective module
---------------------------------------

.. automodule:: simsopt.objectives.fluxobjective
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.objectives.functions module
-----------------------------------

.. automodule:: simsopt.objectives.functions
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.objectives.graph\_functions module
------------------------------------------

.. automodule:: simsopt.objectives.graph_functions
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.objectives.graph\_least\_squares module
-----------------------------------------------

.. automodule:: simsopt.objectives.graph_least_squares
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.objectives.least\_squares module
----------------------------------------

.. automodule:: simsopt.objectives.least_squares
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.objectives
   :members:
   :undoc-members:
   :show-inheritance:
simsopt package
===============

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   simsopt._core
   simsopt.field
   simsopt.geo
   simsopt.mhd
   simsopt.objectives
   simsopt.solve
   simsopt.util
   simsoptpp

Module contents
---------------

.. automodule:: simsopt
   :members:
   :undoc-members:
   :show-inheritance:
simsopt.util package
====================

Submodules
----------

simsopt.util.constants module
-----------------------------

.. automodule:: simsopt.util.constants
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.dev module
-----------------------

.. automodule:: simsopt.util.dev
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.log module
-----------------------

.. automodule:: simsopt.util.log
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.mpi module
-----------------------

.. automodule:: simsopt.util.mpi
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.mpi\_logger module
-------------------------------

.. automodule:: simsopt.util.mpi_logger
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.types module
-------------------------

.. automodule:: simsopt.util.types
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.util.zoo module
-----------------------

.. automodule:: simsopt.util.zoo
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.util
   :members:
   :undoc-members:
   :show-inheritance:
`Source code on GitHub <https://github.com/hiddenSymmetries/booz_xform>`_
=========================================================================

The source code for ``simsopt`` can be found at https://github.com/hiddenSymmetries/simsopt.
MPI partitions and worker groups
--------------------------------

To use MPI parallelization with simsopt, you must install the
`mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_ package.

MPI simsopt programs can generally be launched using ``mpirun`` or
``mpiexec``, e.g.

.. code-block::

   mpiexec -n 32 ./name_of_your_python_driver_script

where 32 can be replaced by the number of processors you want to use.
There may be specific instructions to run MPI programs on your HPC system,
so consult the documentation for your system.
   
Given a set of :math:`M` MPI processes, an optimization problem can be
parallelized in different ways.  At one extreme, all :math:`M`
processes could work together on each evaluation of the objective
function, with the optimization algorithm itself only initiating a
single evaluation of the objective function :math:`f` at a time.  At
the other extreme, the optimization algorithm could request :math:`M`
evaluations of :math:`f` all at once, with each evaluation performed
by a single processor.  This type of parallelization can be useful for
instance when finite differences are used to evaluate gradients.  Or,
both types of parallelization could be used at the same time. To allow
all of these types of parallelization, simsopt uses the concepts of
*worker groups* and *MPI partitions*.

A *worker group* is a set of MPI processes that works together to
evaluate the objective function.  An MPI partition is a division of
the total set of :math:`M` MPI processors into one or more worker
groups.  If there are :math:`W` worker groups, each group has
approximately :math:`M/W` processors (approximate because one may need to
round up or down.)  Simsopt has a class
:obj:`simsopt.util.mpi.MpiPartition` to manage this division of
processors into worker groups.  Each MPI-aware simsopt object, such as
:obj:`simsopt.mhd.vmec.Vmec`, keeps an instance of :obj:`~simsopt.util.mpi.MpiPartition`
which can be set from its constructor.  The number of worker
groups can be set by an argument to :obj:`~simsopt.util.mpi.MpiPartition`.
For instance, to tell a Vmec object to use four worker groups, one could write

.. code-block::

   from simsopt.mhd.vmec import Vmec
   from simsopt.util.mpi import MpiPartition
   
   mpi = MpiPartition(4)
   equil = Vmec("input.li383_low_res", mpi=mpi)

The same :obj:`~simsopt.util.mpi.MpiPartition` instance should be passed to the solver::

  # ... code to define an optimization problem "prob" ...
  
  from simsopt.solve.graph_mpi import least_squares_mpi_solve
  
  least_squares_mpi_solve(prob, mpi, grad=True)

Many optimization algorithms that do not use derivatives do not
support concurrent evaluations of the objective.  In this case, the
number of worker groups, :math:`W`, should be equal to 1.  Any algorithm that
uses derivatives, such as Levenberg-Marquardt, can take advantage of
multiple worker groups to evaluate derivatives by finite
differences. If the number of parameters (i.e. independent variables)
is :math:`N`, you ideally want to set :math:`W=N+1` when using 1-sided
finite differences, and set :math:`W=2N+1` when using centered
differences.  These ideal values are not required however -
``simsopt`` will evaluate finite difference derivatives for any value
of :math:`W`, and results should be exactly independent of :math:`W`.
Other derivative-free algorithms intrinsically support
parallelization, such as HOPSPACK, though no such algorithm is
available in ``simsopt`` yet.

An MPI partition is associated with three MPI communicators, "world",
"groups", and "leaders". The "world" communicator
represents all :math:`M` MPI processors available to the program. (Normally
this communicator is the same as ``MPI_COMM_WORLD``, but it could be a
subset thereof if you wish.)  The "groups" communicator also
contains the same :math:`M` processors, but grouped into colors, with
a different color representing each worker group. Therefore
operations such as ``MPI_Send`` and ``MPI_Bcast`` on this communicator
exchange data only within one worker group.  This "groups"
communicator is therefore the one that must be used for evaluation of
the objective function.  Finally, the "leaders" communicator
only includes the :math:`W` processors that have rank 0 in the
"groups" communicator.  This communicator is used for
communicating data within a parallel optimization *algorithm* (as
opposed to within a parallelized objective function).

Given an instance of :obj:`simsopt.util.mpi.MpiPartition` named
``mpi``, the number of worker groups is available as ``mpi.ngroups``,
and the index of a given processor's group is ``mpi.group``.  The
communicators are available as ``mpi.comm_world``,
``mpi.comm_groups``, and ``mpi.comm_leaders``.  The number of
processors within the communicators can be determined from
``mpi.nprocs_world``, ``mpi.nprocs_groups``, and
``mpi.nprocs_leaders``.  The rank of the present processor within the
communicators is available as ``mpi.rank_world``, ``mpi.rank_groups``,
and ``mpi.rank_leaders``.  To determine if the present processor has
rank 0 in a communicator, you can use the variables
``mpi.proc0_world`` or ``mpi.proc0_groups``, which have type ``bool``.

Publications
============

If you use simsopt in your research, kindly cite simsopt using
reference 1 below.  It is recommended that you also give a reference
for the simsopt release version used in your work. Each release has a
Zenodo archive with a unique DOI. The DOI for the latest release is

   .. image:: https://zenodo.org/badge/247710081.svg
        :target: https://zenodo.org/badge/latestdoi/247710081

Clicking this badge will take you to Zenodo, and on the right of the
page you will find a list of DOIs for previous releases.
   
Here is a list of publications in which simsopt results appear:

.. # The | symbols below are used to put a blank line between each item.

#. | M Landreman, B Medasani, F Wechsung, A Giuliani, R Jorge, and C Zhu,
     "SIMSOPT: A flexible framework for stellarator optimization",
     *J. Open Source Software* **6**, 3525 (2021).
     `[journal version] <https://doi.org/10.21105/joss.03525>`__
   | 

#. | F Wechsung, A Giuliani, M Landreman, A Cerfon, and G Stadler,
     "Single-stage gradient-based stellarator coil design: stochastic optimization",
     *Submitted*, (2021).
     `[arXiv version] <https://arxiv.org/pdf/2106.12137>`__
   |
   
#. | M Landreman, B Medasani, and C Zhu,
     "Stellarator optimization for good magnetic surfaces at the same time as quasisymmetry",
     *Phys. Plasmas* **28**, 092505 (2021).
     `[journal version] <https://doi.org/10.1063/5.0061665>`__
     `[arXiv version] <https://arxiv.org/pdf/2106.14930>`__
   |

#. | A Bader, D T Anderson, M Drevlak, B J Faber, C C Hegna, S Henneberg, M Landreman, J C Schmitt, Y Suzuki, and A Ware,
     "Energetic particle transport in optimized stellarators",
     *Nuclear Fusion* **61**, 116060 (2021).
     `[journal version] <https://doi.org/10.1088/1741-4326/ac2991>`__
     `[arXiv version] <https://arxiv.org/pdf/2106.00716>`__
   |
   
#. | M Landreman and E Paul,
     "Magnetic fields with precise quasisymmetry for plasma confinement",
     *Physical Review Letters* **128**, 035001 (2022).
     `[journal version] <https://doi.org/10.1103/PhysRevLett.128.035001>`__
     `[arXiv version] <https://arxiv.org/pdf/2108.03711>`__
   |

#. | A Baillod, J Loizu, J P Graves, and M Landreman,
     "Stellarator optimization for good magnetic surfaces at finite β and toroidal current",
     *Submitted*, (2021).
     `[arXiv version] <https://arxiv.org/pdf/2111.15564>`__
   |
simsopt.solve package
=====================

Submodules
----------

simsopt.solve.graph\_mpi module
-------------------------------

.. automodule:: simsopt.solve.graph_mpi
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.solve.graph\_serial module
----------------------------------

.. automodule:: simsopt.solve.graph_serial
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.solve.mpi module
------------------------

.. automodule:: simsopt.solve.mpi
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.solve.serial module
---------------------------

.. automodule:: simsopt.solve.serial
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.solve
   :members:
   :undoc-members:
   :show-inheritance:
simsopt.field package
=====================

Submodules
----------

simsopt.field.biotsavart module
-------------------------------

.. automodule:: simsopt.field.biotsavart
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.boozermagneticfield module
----------------------------------------

.. automodule:: simsopt.field.boozermagneticfield
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.coil module
-------------------------

.. automodule:: simsopt.field.coil
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.magneticfield module
----------------------------------

.. automodule:: simsopt.field.magneticfield
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.magneticfieldclasses module
-----------------------------------------

.. automodule:: simsopt.field.magneticfieldclasses
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.sampling module
-----------------------------

.. automodule:: simsopt.field.sampling
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.field.tracing module
----------------------------

.. automodule:: simsopt.field.tracing
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.field
   :members:
   :undoc-members:
   :show-inheritance:
Installation
============

This page provides general information on installation.  Detailed
installation instructions for some specific systems can also be found
on `the wiki <https://github.com/hiddenSymmetries/simsopt/wiki>`_.

Requirements
^^^^^^^^^^^^

``simsopt`` is a python package focused on stellarator optimization and requires
python version 3.7 or higher.  ``simsopt``
also requires some mandatory python packages, listed in
``requirements.txt``.  These packages are all installed automatically
when you install using ``pip`` or another python package manager such as ``conda``, as discussed below.  If you prefer to
install via ``python setup.py install`` or ``python setup.py
develop``, you will need to install these python packages manually
using ``pip`` or ``conda``.

Mandatory Packages
------------------
- numpy
- jax
- jaxlib
- scipy
- nptyping
- ruamel.yaml

Optional Packages
-----------------
- For MPI support:
    * mpi4py
- For SPEC support:
    * py_spec
    * pyoculus
    * h5py
    * f90nml
- For VMEC support:
    * https://github.com/hiddenSymmetries/vmec2000
- For computing Boozer coordinates:
    * `booz_xform <https://hiddensymmetries.github.io/booz_xform/>`_

For requirements of separate physics modules like VMEC, see the
documentation of the module you wish to use.


Virtual Environments
^^^^^^^^^^^^^^^^^^^^


This is an optional step, but users are strongly encouraged to use a python virtual environment
to install simsopt. There are two popular ways to create a python virtual environment using 
either ``venv`` module supplied with python or the conda virtual environment.

venv
----

A python virtual envionment can be created with venv using

.. code-block::

    python3 -m venv <path/to/new/virtual/environment>

Activate the newly created virtual environmnet (for bash shell)

.. code-block::
   
    . <path/to/new/virtual/environment>/bin/activate

If you are on a different shell, use the ``activate`` file with an appropriate extension reflecting the shell type.
For more information, please refer to `venv official documentation <https://https://docs.python.org/3/library/venv.html>`_.

conda
-----
Install either `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `anaconda <https://www.anaconda.com/>`_.
If you are on a HPC system, anaconda is either available by default or via a module.

A conda python virtual environment can be created by running

.. code-block::

    conda create -n <your_virtual_env_name> python=3.8

For the new virtual environment, python version 3.8 was chosen in the above command, but you are free to choose any version you want. 
The newly created virtual environment can be activated with a simple command

.. code-block::

    conda activate <your_virtual_env_name>

After activating the conda virtual environment, the name of the environment should appear in the shell prompt.

Installation methods
^^^^^^^^^^^^^^^^^^^^

PyPi
----

This works for both venv and conda virtual environments.

.. code-block::

    pip install simsopt

Running the above command will install simsopt and all of its mandatory dependencies. To install
optional dependencies related to SPEC and MPI, run the following command:

.. code-block::

    pip install simsopt[MPI,SPEC]
    
On some systems, you may not have permission to install packages to
the default location. In this case, add the ``--user`` flag to ``pip``
so the package can be installed for your user only::

    pip install --user simsopt
    
conda
-----

A pre-compiled conda package for simsopt is available. This
installation approach works only with conda virtual environments.
First we need to add conda-forge as one of the channels.

.. code-block::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

Then install simsopt by running

.. code-block::

    conda install -c hiddensymmetries simsopt


From source
-----------

This approach works for both venv and conda virtual environments.
First, install ``git`` if not already installed. Then clone the repository using

.. code-block::

    git clone https://github.com/hiddenSymmetries/simsopt.git

Then install the package to your local python environment with

.. code-block::

    cd simsopt
    pip install -e .

The ``-e`` flag makes the installation "editable", meaning that the
installed package is a pointer to your local repository rather than
being a copy of the source files at the time of installation. Hence,
edits to code in your local repository are immediately reflected in
the package you can import.

Again, if you do not have permission to install python packages to the
default location, add the ``--user`` flag to ``pip`` so the package
can be installed for your user only::

    pip install --user -e .
    
.. warning::
    Installation from local source creates a directory called **build**. If you are reinstalling simsopt from source after updating the code by making local changes or by git pull, remove the directory **build** before reinstalling.

If you want to build SIMSOPT locally with the optional dependencies,
you can run

.. code-block::

    pip install --user -e .[MPI,SPEC]

However, if you're using a zsh terminal (example: latest Macbook versions),
you'll need to run instead

.. code-block::

    pip install --user -e ".[MPI,SPEC]"


Docker container
----------------

A docker image with simsopt along with its dependencies, VMEC, SPEC,
and BOOZ_XFORM pre-installed is available from docker hub. This
container allows you to use simsopt without having to compile any code
yourself.  After `installing docker
<https://docs.docker.com/get-docker/>`_, you can run the simsopt
container directly from the docker image uploaded to Docker Hub.

.. code-block::

   docker run -it --rm hiddensymmetries/simsopt python

The above command should load the python shell that comes with the
simsopt docker container. When you run it first time, the image is
downloaded automatically, so be patient. More information about using
simsopt with Docker can be found :doc:`here <containers>`.

Post-Installation
^^^^^^^^^^^^^^^^^

If the installation is successful, ``simsopt`` will be added to your
python environment. You should now be able to import the module from
python::

  >>> import simsopt

Defining optimization problems
==============================

The Optimizable class
---------------------

A basic tool for defining optimization problems in simsopt is the
class :obj:`~simsopt._core.graph_optimizable.Optimizable`. Many
classes in simsopt are subclasses of this class.  This parent class
provides several functions.  First, it allows for the parameters of an
object to be either fixed or varied in an optimization, for a useful
name string to be associated with each such degree of freedom, and for
box constraints on each parameter to be set.  Second, the
:obj:`~simsopt._core.graph_optimizable.Optimizable` class manages
dependencies between objects.  For example, if an MHD equilibrium
depends on a :obj:`~simsopt.geo.surface.Surface` object representing
the boundary, the equilibrium object will know it needs to recompute
the equilibrium if the :obj:`~simsopt.geo.surface.Surface` changes.
Third, when a set of objects with dependencies is combined into an
objective function, the
:obj:`~simsopt._core.graph_optimizable.Optimizable` class
automatically combines the non-fixed degrees of freedom into a global
state vector, which can be passed to numerical optimization
algorithms.

Users can create their own optimizable objects in two ways. One method
is to create a standard python function, and apply the
:obj:`simsopt.make_optimizable()` function to it, as explained
below. Or, you can directly subclass
:obj:`simsopt._core.graph_optimizable.Optimizable`.


Optimizable degrees of freedom
------------------------------

..
    A notebook containing the example in this section can be found in
    ~/Box Sync/work21/20211219-01 Simsopt optimizable demo.ipynb
    
The parameters of an object that can potentially be included in the
parameter space for optimization are referred to in simsopt as "dofs"
(degrees of freedom). Each dof is a float; integers are not optimized
in simsopt.  Each dof has several properties: it can be either "fixed"
or "free", it has a string name, and it has upper and lower bounds.
The free dofs are varied in an optimization, whereas the fixed ones
are not.

To demonstrate these functions, we can use a
:obj:`simsopt.geo.curvexyzfourier.CurveXYZFourier` object as a
concrete example::

  >>> from simsopt.geo.curvexyzfourier import CurveXYZFourier
  >>> c = CurveXYZFourier(quadpoints=30, order=1)

This object provides a Fourier representation of a closed curve, as is
commonly used for coil optimization.  The dofs for this object are the
Fourier series amplitudes of the Cartesian coordinates
:math:`(X(\theta), Y(\theta), Z(\theta))` with respect to the
parameter :math:`\theta`. By choosing ``order=1``, only mode numbers 0
and 1 are included.

Each dof has a string name, which can be queried using the
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_dof_names`
property::

  >>> c.local_dof_names

  ['xc(0)', 'xs(1)', 'xc(1)', 'yc(0)', 'ys(1)', 'yc(1)', 'zc(0)', 'zs(1)', 'zc(1)']

Evidently there are nine dofs in this case. For each, the number in
parentheses is the mode number :math:`m`. The values of the dofs can
be read or written to using the
:obj:`~simsopt._core.graph_optimizable.Optimizable.x` property::

  >>> c.x

  array([0., 0., 0., 0., 0., 0., 0., 0., 0.])

  >>> c.x = [1, 0.1, 0, -2, 0, 0.3, 3, 0, 0.2]
  >>> c.x

  array([ 1. ,  0.1,  0. , -2. ,  0. ,  0.3,  3. ,  0. ,  0.2])

  >>> c.x[3]

  -2.0

Although you can use indices to retrieve selected elements of
:obj:`~simsopt._core.graph_optimizable.Optimizable.x`, as in the last
line, you *cannot* assign values to individual elements of
:obj:`~simsopt._core.graph_optimizable.Optimizable.x`, i.e. ``c.x[2] =
7.0`` will not work -- you can only assign an entire array to
:obj:`~simsopt._core.graph_optimizable.Optimizable.x`. You can get or
set individual dofs using their index or string name with the
:obj:`~simsopt._core.graph_optimizable.Optimizable.get()` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.set()` methods::

  >>> c.get(5)

  0.3
  
  >>> c.get('xs(1)')

  0.1

  >>> c.set(7, -0.5)
  >>> c.x
  
  array([ 1. ,  0.1,  0. , -2. ,  0. ,  0.3,  3. , -0.5,  0.2])

  >>> c.set('zc(1)', 0.4)
  >>> c.x

  array([ 1. ,  0.1,  0. , -2. ,  0. ,  0.3,  3. , -0.5,  0.4])

Sometimes we may want to vary a particular dof in an optimization, and
other times we may want to hold that same dof fixed. Some common use
cases for fixing dofs are fixing the major radius or minor radius of a
surface, fixing the high-mode-number modes of a surface, or fixing the
current in a coil.  All dofs in our
:obj:`~simsopt.geo.curvexyzfourier.CurveXYZFourier` object are free by
default. We can fix a dof using the
:obj:`~simsopt._core.graph_optimizable.Optimizable.fix()` method.
When a dof is fixed, it is excluded from the state vector
:obj:`~simsopt._core.graph_optimizable.Optimizable.x`, but you can
still access its value either by name, or with the
:obj:`~simsopt._core.graph_optimizable.Optimizable.full_x` property
(which gives both the free and fixed dofs)::

  >>> c.fix('xc(0)')
  >>> c.x

  array([ 0.1,  0. , -2. ,  0. ,  0.3,  3. , -0.5,  0.4])

  >>> c.full_x

  array([ 1. ,  0.1,  0. , -2. ,  0. ,  0.3,  3. , -0.5,  0.4])

  >>> c.get('xc(0)')

  1.0

To check which dofs are free, you can use the
:obj:`~simsopt._core.graph_optimizable.Optimizable.dofs_free_status`
property. The status of individual dofs can also be checked using
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_fixed` or
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_free`, specify
the dof either using its index or string name ::

  >>> c.dofs_free_status

  array([False,  True,  True,  True,  True,  True,  True,  True,  True])

  >>> c.is_fixed(0)

  True

  >>> c.is_fixed('xc(0)')

  True

  >>> c.is_free('xc(0)')

  False

In addition to
:obj:`~simsopt._core.graph_optimizable.Optimizable.fix()`, you can
also manipulate the fixed/free status of dofs using the functions
:obj:`~simsopt._core.graph_optimizable.Optimizable.unfix()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.fix_all()`, and
:obj:`~simsopt._core.graph_optimizable.Optimizable.unfix_all()`::

  >>> c.fix_all()
  >>> c.x

  array([], dtype=float64)

  >>> c.unfix('yc(0)')
  >>> c.x

  array([-2.])

  >>> c.unfix_all()
  >>> c.x

  array([ 1. ,  0.1,  0. , -2. ,  0. ,  0.3,  3. , -0.5,  0.4])
  
Dependencies
------------

A collection of optimizable objects with dependencies is represented
in simsopt as a directed acyclic graph (DAG): each vertex in the graph
is an instance of an
:obj:`~simsopt._core.graph_optimizable.Optimizable` object, and the
direction of each edge indicates dependency.  An
:obj:`~simsopt._core.graph_optimizable.Optimizable` object can depend
on the dofs of other objects, which are called its parents. The
orignal object is considered a child of the parent objects. An
object's "ancestors" are the an object's parents, their parents, and
so on, i.e. all the objects it depends on.  Note that each dof is
"owned" by only one object, even if multiple objects depend on the
value of that dof.

Many of the functions and properties discussed in the previous section
each have two variants: one that applies just to the dofs owned
directly by an object, and another that applies to the dofs of an
object together with its ancestors. The version that applies just to
the dofs directly owned by an object has a name beginning ``local_``.
For example, analogous to the properties
:obj:`~simsopt._core.graph_optimizable.Optimizable.x` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.dof_names`, which
include all ancestor dofs, there are also properties
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_x` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_dof_names`.
To demonstrate these features, we can consider the following small
collection of objects: a :obj:`simsopt.field.coil.Coil`, which is a
pairing of a :obj:`simsopt.field.coil.Current` with a
:obj:`simsopt.geo.curve.Curve`.  For the latter, we can use the
subclass :obj:`simsopt.geo.curvexyzfourier.CurveXYZFourier` as in the
previous section.  These objects can be created as follows::

  >>> from simsopt.field.coil import Current, Coil
  >>> from simsopt.geo.curvexyzfourier import CurveXYZFourier
  >>>
  >>> current = Current(1.0e4)
  >>> curve = CurveXYZFourier(quadpoints=30, order=1)
  >>> coil = Coil(curve, current)

Here, ``coil`` is a child of ``curve`` and ``current``, and ``curve``
and ``current`` are parents of ``coil``. The corresponding graph looks
as follows:

..
    The original vector graphics for the following figure are on Matt's laptop in
    ~/Box Sync/work21/20211220-01 Simsopt optimizable docs graphs.pptx

.. image:: graph1.png
   :width: 400

(Arrows point from children to parents.) You can access a list of the
parents or ancestors of an object with the ``parents`` or
``ancestors`` attributes::

  >>> coil.parents

  [<simsopt.geo.curvexyzfourier.CurveXYZFourier at 0x1259ac630>,
   <simsopt.field.coil.Current at 0x1259a2040>]

The object ``coil`` does not own any dofs of its own, so its
``local_`` properties return empty arrays, whereas its non-``local_``
properties include the dofs of both of its parents::

  >>> coil.local_dof_names

  []

  >>> coil.dof_names

  ['Current1:x0', 'CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)',
   'CurveXYZFourier1:xc(1)', 'CurveXYZFourier1:yc(0)', 'CurveXYZFourier1:ys(1)',
   'CurveXYZFourier1:yc(1)', 'CurveXYZFourier1:zc(0)', 'CurveXYZFourier1:zs(1)',
   'CurveXYZFourier1:zc(1)']

Note that the names returned by
:obj:`~simsopt._core.graph_optimizable.Optimizable.dof_names` have the
name of the object and a colon prepended, to distinguish which
instance owns the dof. This unique name for each object instance can
be accessed by
:obj:`~simsopt._core.graph_optimizable.Optimizable.name`. For the ``current`` and ``curve`` objects,
since they have no ancestors, their
:obj:`~simsopt._core.graph_optimizable.Optimizable.dof_names` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_dof_names` are the same, except
that the non-``local_`` versions have the object name prepended::

  >>> curve.local_dof_names

  ['xc(0)', 'xs(1)', 'xc(1)', 'yc(0)', 'ys(1)', 'yc(1)', 'zc(0)', 'zs(1)', 'zc(1)']

  >>> curve.dof_names

  ['CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)', 'CurveXYZFourier1:xc(1)',
   'CurveXYZFourier1:yc(0)', 'CurveXYZFourier1:ys(1)', 'CurveXYZFourier1:yc(1)',
   'CurveXYZFourier1:zc(0)', 'CurveXYZFourier1:zs(1)', 'CurveXYZFourier1:zc(1)']

  >>> current.local_dof_names

  ['x0']

  >>> current.dof_names

  ['Current1:x0']

The :obj:`~simsopt._core.graph_optimizable.Optimizable.x` property
discussed in the previous section includes dofs from ancestors. The
related property
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_x` applies
only to the dofs directly owned by an object. When the dofs of a
parent are changed, the
:obj:`~simsopt._core.graph_optimizable.Optimizable.x` property of
child objects is automatically updated::

  >>> curve.x = [1.7, -0.2, 0.1, -1.1, 0.7, 0.3, 1.3, -0.6, 0.5]
  >>> curve.x

  array([ 1.7, -0.2,  0.1, -1.1,  0.7,  0.3,  1.3, -0.6,  0.5])

  >>> curve.local_x

  array([ 1.7, -0.2,  0.1, -1.1,  0.7,  0.3,  1.3, -0.6,  0.5])

  >>> current.x

  array([10000.])

  >>> current.local_x

  array([10000.])

  >>> coil.x

  array([ 1.0e+04,  1.7e+00, -2.0e-01,  1.0e-01, -1.1e+00,  7.0e-01,
        3.0e-01,  1.3e+00, -6.0e-01,  5.0e-01])

  >>> coil.local_x

  array([], dtype=float64)

Above, you can see that
:obj:`~simsopt._core.graph_optimizable.Optimizable.x` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_x`
give the same results for ``curve`` and ``current`` since these objects have no ancestors.
For ``coil``,
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_x`
returns an empty array because ``coil`` does not
own any dofs itself, while
:obj:`~simsopt._core.graph_optimizable.Optimizable.x`
is a concatenation of the dofs of its ancestors.

The functions
:obj:`~simsopt._core.graph_optimizable.Optimizable.get()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.set()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.fix()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.unfix()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_fixed()`, and
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_free()` refer
only to dofs directly owned by an object. If an integer index is
supplied to these functions it must be the local index, and if a
string name is supplied to these functions, it does not have the
object name and colon prepended. So for instance,
``curve.fix('yc(0)')`` works, but
``curve.fix('CurveXYZFourier3:yc(0)')``, ``coil.fix('yc(0)')``, and
``coil.fix('CurveXYZFourier3:yc(0)')`` do not.

When some dofs are fixed in parent objects, these dofs are
automatically removed from the global state vector
:obj:`~simsopt._core.graph_optimizable.Optimizable.x` of a child
object::

  >>> curve.fix_all()
  >>> curve.unfix('zc(0)')
  >>> coil.x

  array([1.0e+04, 1.3e+00])

  >>> coil.dof_names

  ['Current1:x0', 'CurveXYZFourier1:zc(0)']

Thus, the :obj:`~simsopt._core.graph_optimizable.Optimizable.x`
property of a child object is convenient to use as the state vector
for numerical optimization packages, as it automatically combines the
selected degrees of freedom that you wish to vary from all objects
that are involved in the optimization problem. If you wish to get or
set the state vector *including* the fixed dofs, you can use the
properties :obj:`~simsopt._core.graph_optimizable.Optimizable.full_x`
(which includes ancestors) or
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_full_x`
(which does not). The corresponding string labels including the fixed
dofs can be accessed using
:obj:`~simsopt._core.graph_optimizable.Optimizable.full_dof_names` and
:obj:`~simsopt._core.graph_optimizable.Optimizable.local_full_dof_names`::
       
  >>> coil.full_x

  array([ 1.0e+04,  1.7e+00, -2.0e-01,  1.0e-01, -1.1e+00,  7.0e-01,
        3.0e-01,  1.3e+00, -6.0e-01,  5.0e-01])

  >>> coil.full_dof_names

  ['CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)', 'CurveXYZFourier1:xc(1)',
   'CurveXYZFourier1:yc(0)', 'CurveXYZFourier1:ys(1)', 'CurveXYZFourier1:yc(1)',
   'CurveXYZFourier1:zc(0)', 'CurveXYZFourier1:zs(1)', 'CurveXYZFourier1:zc(1)']
  
Realistic optimization problems can have significantly more complicated graphs.
For example, here is the graph for the problem described in the paper
`"Stellarator optimization for good magnetic surfaces at the same time as quasisymmetry",
M Landreman, B Medasani, and C Zhu,
Phys. Plasmas 28, 092505 (2021). <https://doi.org/10.1063/5.0061665>`__

.. image:: graph2.png
   :width: 400


   
Function reference
------------------

The following tables provide a reference for many of the properties
and functions of :obj:`~simsopt._core.graph_optimizable.Optimizable`
objects. Many come in a set of 2x2 variants:

.. list-table:: State vector
   :widths: 20 20 20
   :header-rows: 1
   :stub-columns: 1

   * -
     - Excluding ancestors
     - Including ancestors
   * - Both fixed and free
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_full_x`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.full_x`
   * - Free only
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_x`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.x`

.. list-table:: Number of elements in the state vector
   :widths: 20 20 20
   :header-rows: 1
   :stub-columns: 1

   * -
     - Excluding ancestors
     - Including ancestors
   * - Both fixed and free
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_full_dof_size`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.full_dof_size`
   * - Free only
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_dof_size`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.dof_size`

.. list-table:: String names
   :widths: 20 20 20
   :header-rows: 1
   :stub-columns: 1

   * -
     - Excluding ancestors
     - Including ancestors
   * - Both fixed and free
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_full_dof_names`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.full_dof_names`
   * - Free only
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_dof_names`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.dof_names`

.. list-table:: Whether dofs are free
   :widths: 20 20 20
   :header-rows: 1
   :stub-columns: 1

   * -
     - Excluding ancestors
     - Including ancestors
   * - Both fixed and free
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.local_dofs_free_status`
     - :obj:`~simsopt._core.graph_optimizable.Optimizable.dofs_free_status`
   * - Free only
     - N/A
     - N/A

Other attributes: ``name``, ``parents``, ``ancestors``

Other functions:
:obj:`~simsopt._core.graph_optimizable.Optimizable.get()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.set()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.fix()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.unfix()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_fixed()`,
:obj:`~simsopt._core.graph_optimizable.Optimizable.is_free()`.

       
Caching
-------

Optimizable objects may need to run a relatively expensive
computation, such as computing an MHD equilibrium.  As long as no dofs
change, results can be re-used without re-running the computation.
However if any dofs change, either dofs owned locally or by an
ancestor object, this computation needs to be re-run. Many Optimizable
objects in simsopt therefore implement caching: results are saved,
until the cache is cleared due to changes in dofs.  The
:obj:`~simsopt._core.graph_optimizable.Optimizable` base class
provides a function
:obj:`~simsopt._core.graph_optimizable.Optimizable.recompute_bell()`
to assist with caching. This function is called automatically whenever
dofs of an object or any of its ancestors change. Subclasses of
:obj:`~simsopt._core.graph_optimizable.Optimizable` can overload the
default (empty)
:obj:`~simsopt._core.graph_optimizable.Optimizable.recompute_bell()`
function to manage their cache in a customized way.


Specifying least-squares objective functions
--------------------------------------------

A common use case is to minimize a nonlinear least-squares objective
function, which consists of a sum of several terms. In this case the
:obj:`simsopt.objectives.graph_least_squares.LeastSquaresProblem`
class can be used.  Suppose we want to solve a least-squares
optimization problem in which an
:obj:`~simsopt._core.graph_optimizable.Optimizable` object ``obj`` has
some dofs to be optimized. If ``obj`` has a function ``func()``, we
can define the objective function ``weight * ((obj.func() - goal) **
2)`` as follows::

  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  prob = LeastSquaresProblem.from_tuples([(obj.func, goal, weight)])

Note that the problem was defined using a 3-element tuple of the form
``(function_handle, goal, weight)``.  In this example, ``func()``
could return a scalar, or it could return a 1D numpy array. In the
latter case, ``sum(weight * ((obj.func() - goal) ** 2))`` would be
included in the objective function, and ``goal`` could be either a
scalar or a 1D numpy array of the same length as that returned by
``func()``.  Similarly, we can define least-squares problems with
additional terms with a list of multiple tuples::

  prob = LeastSquaresProblem.from_tuples([(obj1.func1, goal1, weight1),
                                          (obj2.func2, goal2, weight2)])

The corresponding objective funtion is then ``weight1 *
((obj1.func1() - goal1) ** 2) + weight2 * ((obj2.func2() - goal2) **
2)``. The list of tuples can include any mixture of terms defined by
scalar functions and by 1D numpy array-valued functions.  Note that
the function handles that are specified should be members of an
:obj:`~simsopt._core.graph_optimizable.Optimizable` object.  As
:obj:`~simsopt.objectives.graph_least_squares.LeastSquaresProblem` is
a subclass of :obj:`~simsopt._core.graph_optimizable.Optimizable`, the
free dofs of all the objects that go into the objective function are
available in the global state vector ``prob.x``. The overall scalar
objective function is available from
:func:`simsopt.objectives.graph_least_squares.LeastSquaresProblem.objective`.
The vector of residuals before scaling by the ``weight`` factors
``obj.func() - goal`` is available from
:func:`simsopt.objectives.graph_least_squares.LeastSquaresProblem.unweighted_residuals`.
The vector of residuals after scaling by the ``weight`` factors,
``sqrt(weight) * (obj.func() - goal)``, is available from
:func:`simsopt.objectives.graph_least_squares.LeastSquaresProblem.residuals`.

Least-squares problems can also be defined in an alternative way::
  
  prob = LeastSquaresProblem([goal1, goal2, goal3],
                             [weight1, weight2, weight3],
                             [obj1.fn1, obj2.fn2, obj3.fn3])

If you prefer, you can specify
``sigma = 1 / sqrt(weight)`` rather than ``weight`` and use the
``LeastSquaresProblem.from_sigma``  as::

  prob = LeastSquaresProblem.from_sigma([goal1, goal2, goal3],
                                        [sigma1, sigma2, sigma3],
                                        [obj1.fn1, obj2.fn2, obj3.fn3])

Custom objective functions and optimizable objects
--------------------------------------------------

You may wish to use a custom objective function.  The recommended
approach for this is to use
:func:`simsopt._core.graph_optimizable.make_optimizable()`, which can
be imported from the top-level ``simsopt`` module. In this approach,
you first define a standard python function which takes as arguments
any :obj:`~simsopt._core.graph_optimizable.Optimizable` objects that
the function depends on. This function can return a float or 1D numpy
array.  You then apply
:func:`~simsopt._core.graph_optimizable.make_optimizable()` to the
function handle, including the parent objects as additional
arguments. The newly created
:obj:`~simsopt._core.graph_optimizable.Optimizable` object will have a
function ``.J()`` that returns the function you created.

For instance, suppose we wish to minimize the objective function
``(m - 0.1)**2``, where ``m`` is the value of VMEC's ``DMerc`` array
(for Mercier stability) at the outermost available grid point. This
can be accomplished as follows::

  from simsopt import make_optimizable
  from simsopt.mhd.vmec import Vmec
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem

  def myfunc(v):
     v.run()  # Ensure VMEC has run with the latest dofs.
     return v.wout.DMerc[-2]

  vmec = Vmec('input.extension')
  myopt = make_optimizable(myfunc, vmec)
  prob = LeastSquaresProblem.from_tuples([(myopt.J, 0.1, 1)])
      
In this example, the new
:obj:`~simsopt._core.graph_optimizable.Optimizable` object did not own
any dofs.  However the
:func:`~simsopt._core.graph_optimizable.make_optimizable()` can also
create :obj:`~simsopt._core.graph_optimizable.Optimizable` objects
with their own dofs and other parameters. For this syntax, see the API documentation for
:func:`~simsopt._core.graph_optimizable.make_optimizable()`.

An alternative to using
:func:`~simsopt._core.graph_optimizable.make_optimizable()` is to
write your own subclass of
:obj:`~simsopt._core.graph_optimizable.Optimizable`.  In this
approach, the above example looks as follows::
  
  from simsopt._core.graph_optimizable import Optimizable
  from simsopt.mhd.vmec import Vmec
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem

  class Myopt(Optimizable):
      def __init__(self, v):
          self.v = v
	  Optimizable.__init__(self, depends_on=[v])

      def J(self):
          self.v.run()  # Ensure VMEC has run with the latest dofs.
	  return self.v.wout.DMerc[-2]

  vmec = Vmec('input.extension')
  myopt = Myopt(vmec)
  prob = LeastSquaresProblem.from_tuples([(myopt.J, 0.1, 1)])

  
Derivatives
-----------

Simsopt can be used for both derivative-free and derivative-based
optimization. Examples are included in which derivatives are computed
analytically, with adjoint methods, or with automatic differentiation.
Generally, the objects in the :obj:`simsopt.geo` and
:obj:`simsopt.field` modules provide derivative information, while
objects in :obj:`simsopt.mhd` do not, aside from several adjoint
methods in the latter.  For problems with derivatives, the class
:obj:`simsopt._core.derivative.Derivative` is used.  See the API
documentation of this class for details.  This class provides the
chain rule, and automatically masks out rows of the gradient
corresponding to fixed dofs. The chain rule is computed with "reverse
mode", using vector-Jacobian products, which is efficient for cases in
which the objective function is a scalar or a vector with fewer
dimensions than the number of dofs.  For objects that return a
gradient, the gradient function is typically named ``.dJ()``.
Eliminating magnetic islands
============================

In this example, we show how the shape of a boundary magnetic surface
can be adjusted to eliminate magnetic islands inside it, considering a
vacuum field. For this example we will use the SPEC code with a single
radial domain. The geometry comes from a quasi-helically symmetric
configuration developed at the University of Wisconsin.  We will
eliminate the islands by minimizing an objective function involving
Greene's residue for several O-points and X-points, similar to the
approach of Hanson and Cary (1984).

The initial configuration is defined in the SPEC input file
``QH-residues.sp``, which can be found in the ``examples/2_Intermediate/inputs``
directory. If you generate a Poincare plot for this boundary shape by
running standalone SPEC, it can be seen that there is an island chain
corresponding to the :math:`\iota = -8/7` resonance:

..
   Figure generated by Matt with ~/Box Sync/work20/20201231-01-AtenAndSimsopt/aten_poincare_redBoundary

.. image:: QH-residues-before.png
   :width: 400

This island chain can be eliminated by a slight change to the shape of
the boundary magnetic surface. To do this with simsopt, the first step
is to import some items we will need::

  import numpy as np
  from simsopt.mhd.spec import Spec, Residue
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_serial import least_squares_serial_solve

We then create a Spec object based on the input file::

  s = Spec('QH-residues.sp')

Many Fourier modes of the boundary affect the island width, so we can
choose nearly any subset or all of the boundary modes to vary in the
optimization. To keep this example fast, we will pick out just two
modes to vary rather than all of the modes. Modes with high poloidal
mode number m couple particularly strongly to the islands. Here we
will choose to vary two modes with m=6. Since the original
configuration contains only modes up to m=3, we must increase the
number of modes in the boundary shape in order to have m=6 modes
available to vary::

  s.boundary.change_resolution(6, s.boundary.ntor)

Now we can pick out a few modes of the boundary shape to vary in the
optimization::

  s.boundary.fix_all()
  s.boundary.unfix('zs(6,1)')
  s.boundary.unfix('zs(6,2)')

Next, let us define the objective function. The objective function
will be based on the residue defined by Greene (1979). A residue is a
property of a periodic trajectory, such as the field line at the
O-point or X-point of a magnetic island. If a good magnetic surface
exists at a rational value of :math:`\iota`, the residues will be
zero, so minimizing the squares of residues will promote good surface
quality. In simsopt, we can define a Residue object based on a Spec
object, together with the rational number of the island chain::

  # Resonant surface is iota = p / q:
  p = -8
  q = 7
  # Guess for radial location of the island chain:
  s_guess = 0.9
  residue1 = Residue(s, p, q, s_guess=s_guess)

The initial guess for the radial coordinate s is not critical; the
Newton method to find the periodic field line is robust in this
case. By default, the residue will be computed at a poloidal angle of
zero, corresponding to the O-point of the island chain. We can define
a second residue for the X-point, by specifying the poloidal angle to
be :math:`\pi` instead of zero::

  residue2 = Residue(s, p, q, s_guess=s_guess, theta=np.pi)

To get a numerical value for the residues, you can call the ``J()``
function of these objects. Doing so will automatically cause SPEC to
run::

  initial_r1 = residue1.J()
  initial_r2 = residue2.J()
  print(f"Initial residues: {initial_r1}, {initial_r2}")

After SPEC runs, the output should be close to the following::

  Initial residues: 0.02331532869136166, -0.02287637681580268
  
You are free to also define residues for other rational surfaces in
the confinement region. In this example it turns out to be useful to
also control the residues of the :math:`\iota=-12/11` surface, since
boundary adjustments to eliminate the -8/7 island can open up -12/11
islands. Therefore let us define Residue objects for the O and X
points at this second rational surface::

  p = -12
  q = 11
  s_guess = -0.1

  residue3 = Residue(s, p, q, s_guess=s_guess)
  residue4 = Residue(s, p, q, s_guess=s_guess, theta=np.pi)

We now combine the four residues into a least-squares objective
function, by summing the squares of the residues::

  # Objective function is \sum_j residue_j ** 2
  prob = LeastSquaresProblem.from_tuples([(residue1, 0, 1),
                                          (residue2, 0, 1),
                                          (residue3, 0, 1),
                                          (residue4, 0, 1)])

If you wanted an island to be present instead of absent, which might
be the case when designing an island divertor, a value other than zero
could be used for the goal values above, e.g. ``(residue1, 0.1, 1)``.

Finally, let us solve the optimization problem::

  least_squares_serial_solve(prob)

The solution takes about 18 function evaluations, which likely will
take a minute or two.  Afterward, we can examine the optimum::

  final_r1 = residue1.J()
  final_r2 = residue2.J()
  print(f"Final residues: {final_r1}, {final_r2}")

The residues have been reduced::
  
  Final residues: 2.9093984016959062e-06, 2.5974339906698063e-06

Generating a Poincare plot of the final configuration using standalone
SPEC, the island chain has been eliminated:

..
   Figure generated by Matt with ~/Box Sync/work20/20201231-01-AtenAndSimsopt/aten_poincare_optimized

.. image:: QH-residues-after.png
   :width: 400

(Note that to make Poincare plots like this with SPEC, you can
increase the values of ``nppts`` and ``nptrj`` in the SPEC input
file.)
simsoptpp package
=================


Module contents
---------------

.. automodule:: simsoptpp
   :members:
   :undoc-members:
   :show-inheritance:
Overview
========

Ways to use simsopt
-------------------

Simsopt is a collection of classes and functions that can be used in
several ways.  One application is to solve optimization problems
involving stellarators, similar to STELLOPT.  You could also define an
objective function using simsopt, but use an optimization library
outside simsopt to solve it.  Or, you could use the simsopt
optimization infrastructure to optimize your own objects, which may or
may not have any connection to stellarators.  Alternatively, you can
use the stellarator-related objects in a script for some purpose other
than optimization, such as plotting how some code output varies as an
input parameter changes, or evaluating the finite-difference gradient
of some code outputs.  Or, you can manipulate the objects
interactively, at the python command line or in a Jupyter notebook.


Input files
-----------

Simsopt does not use input data files to define optimization problems,
in contrast to ``STELLOPT``. Rather, problems are specified using a
python driver script, in which objects are defined and
configured. However, objects related to specific physics codes may use
their own input files. In particular, a :obj:`simsopt.mhd.vmec.Vmec` object
can be initialized using a standard VMEC ``input.*`` input file, and a
:obj:`simsopt.mhd.spec.Spec` object can be initialized using a standard
SPEC ``*.sp`` input file.


Optimization stages
-------------------

Recent optimized stellarators have been designed using two stages,
both of which can be performed using simsopt. In the first stage, the
parameter space is the shape of a toroidal boundary flux
surface. Coils are not considered explicitly in this stage.  The
objective function involves surrogates for confinement and stability
in the plasma inside the boundary surface.  In the second optimization
stage, coil shapes are optimized to produce the plasma shape that
resulted from stage 1.  The parameter space for stage 2 represents the
space of coil shapes. The objective function for stage 2 usually
involves several terms.  One term is the deviation between the
magnetic field produced by the coils and the magnetic field desired at
the plasma boundary, given the stage 1 solution. Other terms in the
objective function introduce regularization on the coil shapes, such
as the coil length and/or curvature, and reflect other engineering
considerations such as the distance between coils. In the future, we
aim to introduce alternative optimization strategies in simsopt
besides this two-stage approach, such as combined single-stage
methods.



Optimization
------------

To do optimization using simsopt, there are four basic steps:

1. Define the physical entities in the optimization problem (coils, MHD equilibria, etc.) by creating instances of the relevant simsopt classes.
2. Define the independent variables for the optimization, by choosing which degrees of freedom of these objects are free vs fixed.
3. Define an objective function.
4. Solve the optimization problem that has been defined.

This pattern is evident in the tutorials in this documentation
and in the ``examples`` directory of the repository.

Some typical objects are a MHD equilibrium represented by the VMEC or
SPEC code, or some electromagnetic coils. To define objective
functions, a variety of additional objects can be defined that depend
on the MHD equilibrium or coils, such as a
:obj:`simsopt.mhd.boozer.Boozer` object for Boozer-coordinate
transformation, a :obj:`simsopt.mhd.spec.Residue` object to represent
Greene's residue of a magnetic island, or a
:obj:`simsopt.geo.objectives.LpCurveCurvature` penalty on coil
curvature.

More details about setting degrees of freedom and defining
objective functions can be found on the :doc:`optimizable` page.

For the solution step, two functions are provided presently,
:meth:`simsopt.solve.graph_serial.least_squares_serial_solve` and
:meth:`simsopt.solve.graph_mpi.least_squares_mpi_solve`.  The first
is simpler, while the second allows MPI-parallelized finite differences
to be used in the optimization.


Modules
-------

Classes and functions in simsopt are organized into several modules:

- :obj:`simsopt.geo` contains several representations of curves and surfaces.
- :obj:`simsopt.field` contains machinery for the Biot-Savart law and other magnetic field representations.
- :obj:`simsopt.mhd` contains interfaces to MHD equilibrium codes and tools for diagnosing their output.
- :obj:`simsopt.objectives` contains tools for some common objective functions.
- :obj:`simsopt.solve` contains wrappers for some optimization algorithms.
- :obj:`simsopt.util` contains other utility functions.
- :obj:`simsopt._core` defines the ``Optimizable`` class and other tools used internally in simsopt.
Simsopt documentation
=====================

``simsopt`` is a framework for optimizing `stellarators
<https://en.wikipedia.org/wiki/Stellarator>`_.  The high-level
routines are in python, with calls to C++ or fortran where needed for
performance. Several types of components are included:

- Interfaces to physics codes, e.g. for MHD equilibrium.
- Tools for defining objective functions and parameter spaces for
  optimization.
- Geometric objects that are important for stellarators -- surfaces and
  curves -- with several available parameterizations.
- Efficient implementations of the Biot-Savart law and other magnetic
  field representations, including derivatives.
- Tools for parallelized finite-difference gradient calculations.

The design of ``simsopt`` is guided by several principles:

- Thorough unit testing, regression testing, and continuous
  integration.
- Extensibility: It should be possible to add new codes and terms to
  the objective function without editing modules that already work,
  i.e. the `open-closed principle
  <https://en.wikipedia.org/wiki/Open%E2%80%93closed_principle>`_
  . This is because any edits to working code can potentially introduce bugs.
- Modularity: Physics modules that are not needed for your
  optimization problem do not need to be installed. For instance, to
  optimize SPEC equilibria, the VMEC module need not be installed.
- Flexibility: The components used to define an objective function can
  be re-used for applications other than standard optimization. For
  instance, a ``simsopt`` objective function is a standard python
  function that can be plotted, passed to optimization packages
  outside of ``simsopt``, etc.

``simsopt`` is fully open-source, and anyone is welcome to use it,
make suggestions, and contribute.

Some of the physics modules with compiled code reside in separate
repositories. These separate modules include

- `VMEC <https://github.com/hiddenSymmetries/VMEC2000>`_, for MHD
  equilibrium.
- `SPEC <https://github.com/PrincetonUniversity/SPEC>`_, for MHD
  equilibrium.
- `booz_xform <https://hiddensymmetries.github.io/booz_xform/>`_, for
  Boozer coordinates and quasisymmetry.
  
We gratefully acknowledge funding from the `Simons Foundation's Hidden
symmetries and fusion energy project
<https://hiddensymmetries.princeton.edu>`_.

``simsopt`` is one of several available systems for stellarator
optimization.  Others include `STELLOPT
<https://github.com/PrincetonUniversity/STELLOPT>`_, `ROSE
<https://doi.org/10.1088/1741-4326/aaed50>`_, and `LASSO
<https://gitlab.com/wistell>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   overview
   installation
   containers
   optimizable
   geo
   fields
   tracing
   mpi
   testing
   source
   publications
   contributing

.. toctree::
   :maxdepth: 3
   :caption: Tutorials

   example_vmec_only
   example_quasisymmetry
   example_islands
   example_coils

.. toctree::
   :maxdepth: 3
   :caption: API

   simsopt

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Optimizing an equilibrium code
==============================

Here we walk through an example of a small optimization problem in
which a fixed-boundary MHD equilibrium code is used, but no other
codes are required.  Either VMEC or SPEC can be used, and both
variants of the problem are shown here.  There are two independent
variables for the optimization, controlling the shape of a toroidal
boundary, and the objective function involves the rotational transform
and volume.  This same example is documented in the `stellopt
scenarios collection
<https://github.com/landreman/stellopt_scenarios/tree/master/2DOF_vmecOnly_targetIotaAndVolume>`_.
You can also find the source code for this example in the `examples
directory
<https://github.com/hiddenSymmetries/simsopt/blob/master/examples/stellarator_benchmarks/2DOF_vmecOnly_targetIotaAndVolume.py>`_.

For this problem the two independent variables are ``RBC(1,1)`` and
``ZBS(1,1)``, which control the shape of the plasma boundary::
  
  R(phi) = 1 + 0.1 * cos(theta) + RBC(1,1) * cos(theta - 5 * phi),
  Z(phi) =     0.1 * sin(theta) + ZBS(1,1) * sin(theta - 5 * phi).

Note that this boundary has five field periods. We consider the vacuum
field inside this boundary magnetic surface, i.e. there is no plasma
current or pressure.  The objective function is

.. code-block::
   
   f = (iota - iota_target)^2 + (volume - volume_target)^2,
   
where
``iota`` is the rotational transform on the magnetic axis,
``iota_target = 0.41``,
``volume_target = 0.15 m^3``,
and ``volume`` is the volume enclosed by the plasma boundary.

Here is what the objective function landscape looks like:

.. image:: 2DOF_vmecOnly_targetIotaAndVolume_landscape.png
   :width: 500

It can be seen that the total objective function has two long narrow
valleys that are fairly straight.  There are two symmetric optima, one
at the bottom of each valley.  When either of the two independent
variables is +/- 0.1m, the boundary surface becomes infinitesmally
thin, so equilibrium codes are likely to fail.
	   
The optimum configuration looks as follows, with the color indicating
the magnetic field strength:

.. image:: 2DOF_vmecOnly_targetIotaAndVolume_optimum.png
   :width: 500

	   
SPEC version
------------

Here we show the SPEC version of this example.
For simplicity, MPI parallelization will not be used for now.
To start, we must import several classes::

  from simsopt.mhd import Spec
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_serial import least_squares_serial_solve

Then we create the equilibrium object, starting from an input file::

  equil = Spec('2DOF_targetIotaAndVolume.sp')

This file can be found in the ``examples/stellarator_benchmarks/inputs`` directory; you can prepend
the path to the filename if needed.

Next, we define the independent variables for the optimization, by
choosing which degrees of freedom are fixed or not fixed. In this
case, the independent variables are two of the Fourier amplitudes
defining the boundary. We choose the independent variables as follows::

  surf = equil.boundary
  surf.fix_all()
  surf.unfix('rc(1,1)')
  surf.unfix('zs(1,1)')

The next step is to define the objective function. For details of how to define
least-squares objective functions, see the :doc:`optimizable` page. For the present problem, we use

.. code-block::

   desired_volume = 0.15
   volume_weight = 1
   term1 = (equil.volume, desired_volume, volume_weight)

   desired_iota = -0.41
   iota_weight = 1
   term2 = (equil.iota, desired_iota, iota_weight)

   prob = LeastSquaresProblem.from_tuples([term1, term2])

Finally, we solve the optimization problem::

  least_squares_serial_solve(prob)

SPEC will then run many times; it will likely take a bit less than a
minute to find the optimum.  Once the problem is solved, we can
examine some properties of the optimum::

  print("At the optimum,")
  print(" rc(m=1,n=1) = ", surf.get_rc(1, 1))
  print(" zs(m=1,n=1) = ", surf.get_zs(1, 1))
  print(" volume, according to SPEC    = ", equil.volume())
  print(" volume, according to Surface = ", surf.volume())
  print(" iota on axis = ", equil.iota())
  print(" objective function = ", prob.objective())

The results are

.. code-block::
   
   At the optimum,
    rc(m=1,n=1) =  0.03136534181915223
    zs(m=1,n=1) =  -0.03127549335108014
    volume, according to SPEC    =  0.17802858467026614
    volume, according to Surface =  0.1780285846702657
    iota on axis =  -0.41148381548239504
    objective function =  0.0007878032670040736

These numbers match the solution found using stellopt and VMEC in
`stellopt_scenarios
<https://github.com/landreman/stellopt_scenarios/tree/master/2DOF_vmecOnly_targetIotaAndVolume>`_

    
VMEC version
------------

To use VMEC instead of SPEC, the only essential change is to use a
:obj:`simsopt.mhd.vmec.Vmec` object for the equilibrium instead of the
Spec object.

Here we can also show how to add MPI to the example.  MPI can be used
for parallelized finite-difference gradients, within each VMEC
computation, or both at the same time.  To introduce MPI we first
initialize an :obj:`simsopt.util.mpi.MpiPartition` object and choose
the number of worker groups.  The instance is then passed as an
argument to the Vmec object and to the
:meth:`simsopt.solver.mpi_solve.least_squares_mpi_solve` function.
For more details about MPI, see :doc:`mpi`.

The complete example is then as follows::

  from simsopt.util.mpi import MpiPartition
  from simsopt.mhd import Vmec
  from simsopt.objectives.graph_least_squares import LeastSquaresProblem
  from simsopt.solve.graph_mpi import least_squares_mpi_solve

  # In the next line, we can adjust how many groups the pool of MPI
  # processes is split into.
  mpi = MpiPartition(ngroups=3)

  # Initialize VMEC from an input file:
  equil = Vmec('input.2DOF_vmecOnly_targetIotaAndVolume', mpi)
  surf = equil.boundary

  # You can choose which parameters are optimized by setting their 'fixed' attributes.
  surf.fix_all()
  surf.unfix('rc(1,1)')
  surf.unfix('zs(1,1)')

  # Each Target is then equipped with a shift and weight, to become a
  # term in a least-squares objective function.  A list of terms are
  # combined to form a nonlinear-least-squares problem.
  desired_volume = 0.15
  volume_weight = 1
  term1 = (equil.volume, desired_volume, volume_weight)

  desired_iota = 0.41
  iota_weight = 1
  term2 = (equil.iota_axis, desired_iota, iota_weight)

  prob = LeastSquaresProblem.from_tuples([term1, term2])

  # Solve the minimization problem:
  least_squares_mpi_solve(prob, mpi, grad=True)

The VMEC input file used here can be found in the ``examples``
directory of the repository.
Geometric objects
-----------------

Simsopt contains implementations of two types of geometric objects
that are important for stellarators: curves and surfaces. Each are
available in a variety of parameterizations, and several derivatives
are provided.  The curve and surface objects can all be found in the
:obj:`simsopt.geo` module.

.. _curves:
     
Curves
~~~~~~

Curves are useful for representing electromagnetic coils and the
magnetic axis.  Curves are represented in simsopt with subclasses of
the base class :obj:`simsopt.geo.curve.Curve`.  A simsopt curve is
modelled as a function :math:`\Gamma:[0, 1] \to \mathbb{R}^3`.  Note
that the curve parameter goes up to 1, not to :math:`2\pi`.  Curves in
simsopt are assumed to be periodic in the parameter. A curve object
stores a list of :math:`n_\theta` "quadrature points" :math:`\{\theta_1,
\ldots, \theta_{n_\theta}\} \subset [0, 1]`.  A variety of class methods
return information about the curve at these quadrature points. Some of
the available methods are the following:

- ``Curve.gamma()``: returns a ``(n_theta, 3)`` array containing :math:`\Gamma(\theta_i)` for :math:`i\in\{1, \ldots, n_\theta\}`, i.e. returns a list of XYZ coordinates along the curve.
- ``Curve.gammadash()``: returns a ``(n_theta, 3)`` array containing :math:`\Gamma'(\theta_i)` for :math:`i\in\{1, \ldots, n_\theta\}`, i.e. returns the (non-unit-length) tangent vector along the curve.
- ``Curve.kappa()``: returns a ``(n_theta, 1)`` array containing the curvature :math:`\kappa` of the curve at the quadrature points.
- ``Curve.torsion()``: returns a ``(n_theta, 1)`` array containing the torsion :math:`\tau` of the curve at the quadrature points.
- ``Curve.frenet_frame()``: returns a 3-element tuple. The leading element is a ``(n_theta, 3)`` array containing the Cartesian components of the unit tangent vector at the quadrature points. Similarly, the remaining two entries of the tuple give the unit normal and binormal vectors.

The different curve classes, such as
:obj:`simsopt.geo.curverzfourier.CurveRZFourier` and
:obj:`simsopt.geo.curvexyzfourier.CurveXYZFourier` differ in the way
curves are discretized.  Each of these take an array of "dofs"
(parameters, e.g. Fourier coefficients) and turn these into a function
:math:`\Gamma:[0, 1] \to \mathbb{R}^3`.  These dofs can be queried and
set via the ``.x`` or ``.full_x`` properties; the former gives just
the non-fixed dofs, whereas the latter gives all dofs including those
that are fixed. You can also set an individual dof using its string
name.  For example, for a
:obj:`simsopt.geo.curvexyzfourier.CurveXYZFourier` object ``c``, the
zero-frequency Fourier mode of the y Cartesian component can be set to
5.0 using ``c.set("yc(0)", 5.0)``.  Changing the dofs will change the
shape of the curve. Simsopt is able to compute derivatives of all
relevant quantities with respect to the discretization parameters.
For example, to compute the derivative of the coordinates of the curve
at quadrature points, one calls ``Curve.dgamma_by_dcoeff()``.  One
obtains a numpy array of shape ``(n_theta, 3, n_dofs)``, containing the
derivative of the position at every quadrature point with respect to
every degree of freedom of the curve.  In the same way one can compute
the derivative of quantities such as curvature (via
``Curve.dkappa_by_dcoeff()``) or torsion (via
``Curve.dtorsion_by_dcoeff()``).

A number of quantities are implemented in
:obj:`simsopt.geo.curveobjectives` and are computed on a
:obj:`simsopt.geo.curve.Curve`:

- ``CurveLength``: computes the length of the ``Curve``.
- ``LpCurveCurvature``: computes a penalty based on the :math:`L_p` norm of the curvature on a curve.
- ``LpCurveTorsion``: computes a penalty based on the :math:`L_p` norm of the torsion on a curve.
- ``MinimumDistance``: computes a penalty term on the minimum distance between a set of curves.

The value of the quantity and its derivative with respect to the curve
dofs can be obtained by calling e.g., ``CurveLength.J()`` and
``CurveLength.dJ()``.

.. _surfaces:

Surfaces
~~~~~~~~

Surfaces are used to represent flux surfaces, particularly for the
boundary of MHD equilibria, and for the target surface in stage-2 coil
optimization.  Surfaces are represented in simsopt using subclasses of
the base class :obj:`simsopt.geo.surface.Surface`.  A surface is
modelled in simsopt as a function :math:`\Gamma:[0, 1] \times [0, 1]
\to \mathbb{R}^3` and is evaluated at quadrature points
:math:`\{\phi_1, \ldots, \phi_{n_\phi}\}\times\{\theta_1, \ldots,
\theta_{n_\theta}\}`.  Here, :math:`\phi` is the toroidal angle and
:math:`\theta` is the poloidal angle. Note that :math:`\phi` and
:math:`\theta` go up to 1, not up to :math:`2 \pi`! Surfaces in
simsopt are assumed to be periodic in both angles.

In practice, you almost never use the base
:obj:`~simsopt.geo.surface.Surface` class.  Rather, you typically use
one of the subclasses corresponding to a specific parameterization.
Presently, the available subclasses are
:obj:`~simsopt.geo.surfacerzfourier.SurfaceRZFourier`,
:obj:`~simsopt.geo.surfacegarabedian.SurfaceGarabedian`,
:obj:`~simsopt.geo.surfacehenneberg.SurfaceHenneberg`,
:obj:`~simsopt.geo.surfacexyzfourier.SurfaceXYZFourier`,
and
:obj:`~simsopt.geo.surfacexyztensorfourier.SurfaceXYZTensorFourier`.
In many cases you can convert a surface from one type to another by going through
:obj:`~simsopt.geo.surfacerzfourier.SurfaceRZFourier`, as most surface types have
``to_RZFourier()`` and ``from_RZFourier()`` methods.
Note that :obj:`~simsopt.geo.surfacerzfourier.SurfaceRZFourier`
corresponds to the surface parameterization used internally in the VMEC and SPEC codes.
However when using these codes in simsopt, any of the available surface subclasses
can be used to represent the surfaces, and simsopt will automatically handle the conversion
to :obj:`~simsopt.geo.surfacerzfourier.SurfaceRZFourier` when running the code.

The points :math:`\phi_j` and :math:`\theta_j` are used for evaluating
the position vector and its derivatives, for computing integrals, and
for plotting, and there are several available methods to specify these
points.  For :math:`\theta_j`, you typically specify a keyword
argument ``ntheta`` to the constructor when instantiating a surface
class. This results in a grid of ``ntheta`` uniformly spaced points
between 0 and 1, with no endpoint at 1. Alternatively, you can specify
a list or array of points to the ``quadpoints_theta`` keyword argument
when instantiating a surface class, specifying the :math:`\theta_j`
directly.  If both ``ntheta`` and ``quadpoints_theta`` are specified,
an exception will be raised.  For the :math:`\phi` coordinate, you
sometimes want points up to 1 (the full torus), sometimes up to
:math:`1/n_{fp}` (one field period), and sometimes up to :math:`1/(2
n_{fp})` (half a field period). These three cases can be selected by
setting the ``range`` keyword argument of the surface subclasses to
``"full torus"``, ``"field period"``, or ``"half period"``.
Equivalently, you can set ``range`` to the constants
``S.RANGE_FULL_TORUS``, ``S.RANGE_FIELD_PERIOD``, or
``S.RANGE_HALF_PERIOD``, where ``S`` can be
:obj:`simsopt.geo.surface.Surface` or any of its subclasses.  Note
that the :math:`\phi` grid points begin at 0 for ``"full torus"`` and
``"field period"``, whereas for ``"half period"`` the :math:`\phi`
grid is shifted by half of the grid spacing to preserve spectral
accuracy of integration.  For all three cases, the ``nphi`` keyword
argument can be set to the desired number of :math:`\phi` grid
points. Alternatively, you can pass a list or array to the
``quadpoints_phi`` keyword argument of the constructor for any Surface
subclass to specify the :math:`\phi_j` points directly.  An exception
will be raised if both ``nphi`` and ``quadpoints_phi`` are specified.
For more information about these arguments, see the
:obj:`~simsopt.geo.surfacerzfourier.SurfaceRZFourier` API
documentation.

The methods available to each surface class are similar to those of
the :obj:`~simsopt.geo.curve.Curve` class:

- ``Surface.gamma()``: returns a ``(n_phi, n_theta, 3)`` array containing :math:`\Gamma(\phi_i, \theta_j)` for :math:`i\in\{1, \ldots, n_\phi\}, j\in\{1, \ldots, n_\theta\}`, i.e. returns a list of XYZ coordinates on the surface.
- ``Surface.gammadash1()``: returns a ``(n_phi, n_theta, 3)`` array containing :math:`\partial_\phi \Gamma(\phi_i, \theta_j)` for :math:`i\in\{1, \ldots, n_\phi\}, j\in\{1, \ldots, n_\theta\}`.
- ``Surface.gammadash2()``: returns a ``(n_phi, n_theta, 3)`` array containing :math:`\partial_\theta \Gamma(\phi_i, \theta_j)` for :math:`i\in\{1, \ldots, n_\phi\}, j\in\{1, \ldots, n_\theta\}`.
- ``Surface.normal()``: returns a ``(n_phi, n_theta, 3)`` array containing :math:`\partial_\phi \Gamma(\phi_i, \theta_j)\times \partial_\theta \Gamma(\phi_i, \theta_j)` for :math:`i\in\{1, \ldots, n_\phi\}, j\in\{1, \ldots, n_\theta\}`.
- ``Surface.area()``: returns the surface area.
- ``Surface.volume()``: returns the volume enclosed by the surface.

A number of quantities are implemented in :obj:`simsopt.geo.surfaceobjectives` and are computed on a :obj:`simsopt.geo.surface.Surface`:

- ``ToroidalFlux``: computes the flux through a toroidal cross section of a ``Surface``.

The value of the quantity and its derivative with respect to the surface dofs can be obtained by calling e.g., ``ToroidalFlux.J()`` and ``ToroidalFlux.dJ_dsurfacecoefficients()``.


Caching
~~~~~~~

The quantities that Simsopt can compute for curves and surfaces often
depend on each other.  For example, the curvature or torsion of a
curve both rely on ``Curve.gammadash()``; to avoid repeated
calculation, geometric objects contain a cache that is automatically
managed.  If a quantity for the curve is requested, the cache is
checked to see whether it was already computed.  This cache can be
cleared manually by calling ``Curve.invalidate_cache()``.  This
function is called every time values are assigned to ``Curve.x``
(meaning the shape of the curve changes).

Graphics
~~~~~~~~

Some basic graphics functions are provided for curve and surface
objects.  To plot a single curve or surface, you can call the
``.plot()`` function of the object.  Presently, three graphics engines
are supported: matplotlib, mayavi, and plotly.  You can select the
plotting engine by passing the ``engine`` keyword argument, e.g. if
``c`` is a Curve object you can call ``c.plot(engine="mayavi")``. You
can use the ``close`` argument to control whether segments are drawn
between the last quadrature point and the first. For these and other
options, see the API documentation for
:func:`simsopt.geo.curve.Curve.plot()` and
:func:`simsopt.geo.surface.Surface.plot()`.

If you have multiple curve and/or surface objects, a convenient way to
plot them together on the same axes is the function
:func:`simsopt.geo.plot.plot()`, which accepts a list of objects as
its argument. Any keywords passed to this function are passed to the
``.plot()`` methods of the individual objects, so you may wish to pass
keywords such as ``engine`` or ``close``.  Alternatively, you can also
use the ``ax`` and ``show`` arguments of the ``.plot()`` methods for
individual curve and surface objects to put them on shared axes.

It is also possible to export curve and surface objects in VTK format,
so they can be viewed in Paraview.  This functionality requires the
python package ``pyevtk``, which can be installed via ``pip install
pyevtk``. A list of curve objects can be exported using the function
:func:`simsopt.geo.curve.curves_to_vtk()`. To export a VTK file for a
surface, call the ``.to_vtk(filename)`` function of the object.  See
:func:`simsopt.geo.surface.Surface.to_vtk()` for more details.


Coil optimization
=================

Here we work through an example of "stage 2" coil optimization.  In
this approach, we consider the target plasma boundary shape to have
been already chosen in a "stage 1", and the present task is to
optimize the shapes of coils to produce this target field. The
tutorial here is similar to the example
``examples/2_Intermediate/stage_two_optimization.py``.  Note that for
this stage-2 problem, no MHD codes like VMEC or SPEC are used, so you
do not need to have them installed.

The objective function we will minimize is

.. math::

   J = \frac{1}{2} \int |\vec{B}\cdot\vec{n}|^2 ds + \alpha \sum_j L_j  + \beta J_{dist}

The first right-hand-side term is the "quadratic flux", the area
integral over the target plasma surface of the square of the magnetic
field normal to the surface. If the coils exactly produce a flux
surface of the target shape, this term will vanish.  Next, :math:`L_j`
is the length of the :math:`j`-th coil.  The scalar regularization
parameter :math:`\alpha` is chosen to balance a trade-off: large
values will give smooth coils at the cost of inaccuracies in producing
the target field; small values of :math:`\alpha` will give a more
accurate match to the target field at the cost of complicated coils.
Finally, :math:`J_{dist}` is a penalty that prevents coils from
becoming too close.  For its precise definition, see
:obj:`simsopt.geo.curveobjectives.MinimumDistance`.  The constant
:math:`\beta` is selected to balance this minimum distance penalty
against the other objectives.  Analytic derivatives are used for the
optimization.

In this tutorial we will consider vacuum fields, so the magnetic field
due to current in the plasma does not need to be subtracted in the
quadratic flux term. The configuration considered is the "precise QA"
case from `arXiv:2108.03711 <http://arxiv.org/pdf/2108.03711.pdf>`_,
which has two field periods and is quasi-axisymmetric.

The stage-2 optimization problem is automatically parallelized in
simsopt using OpenMP and vectorization, but MPI is not used, so the
``mpi4py`` python package is not needed. This example can be run in a
few seconds on a laptop.

To begin solving this optimization problem in simsopt, we first import
some classes that will be used::

  import numpy as np
  from scipy.optimize import minimize
  from simsopt.geo.surfacerzfourier import SurfaceRZFourier
  from simsopt.objectives.fluxobjective import SquaredFlux
  from simsopt.geo.curve import curves_to_vtk, create_equally_spaced_curves
  from simsopt.field.biotsavart import BiotSavart
  from simsopt.field.coil import Current, coils_via_symmetries
  from simsopt.geo.curveobjectives import CurveLength, MinimumDistance
  from simsopt.geo.plot import plot

The target plasma surface is given in the VMEC input file ``tests/test_files/input.LandremanPaul2021_QA``.
We load the surface using

.. code-block::

  nphi = 32
  ntheta = 32
  filename = "tests/test_files/input.LandremanPaul2021_QA"
  s = SurfaceRZFourier.from_vmec_input(filename, range="half period", nphi=nphi, ntheta=ntheta)

You can adjust the directory in ``"filename"`` as appropriate for your
system. The surface could also be read in from a VMEC wout file using
:obj:`simsopt.geo.surfacerzfourier.SurfaceRZFourier.from_wout()`.
(VMEC does not need to be installed to initialize a surface from a
VMEC input or output file.)  Note that surface objects carry a grid of
"quadrature points" at which the position vector is evaluated, and in
different circumstances we may want these points to cover different
ranges of the toroidal angle. For this problem with stellarator
symmetry and field-period symmetry, we need only consider half of a
field period in order to evaluate integrals over the entire
surface. For this reason, the ``range`` parameter of the surface is
set to ``"half period"`` here.

Next, we set the initial condition for the coils, which will be equally spaced circles.
Here we will consider a case with four unique coil shapes, each of which is repeated four times due to
stellarator symmetry and two-field-period symmetry, giving a total of 16 coils.
The four unique coil shapes are called the "base coils". Each copy of a coil also carries the same current,
but we will allow the unique coil shapes to have different current from each other,
as is allowed in W7-X. For this tutorial we will consider the coils to be infinitesimally thin filaments.
In simsopt, such a coil is represented with the :obj:`simsopt.field.coil.Coil` class,
which is essentially a curve paired with a current, represented using
:obj:`simsopt.geo.curve.Curve` and :obj:`simsopt.field.coil.Current` respectively.
The initial conditions are set as follows::

  # Number of unique coil shapes:
  ncoils = 4

  # Major radius for the initial circular coils:
  R0 = 1.0
  
  # Minor radius for the initial circular coils:
  R1 = 0.5

  # Number of Fourier modes describing each Cartesian component of each coil:
  order = 5

  base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=R0, R1=R1, order=order)
  base_currents = [Current(1e5) for i in range(ncoils)]

One detail of optimizing coils for a vacuum configuration is that the
optimizer can "cheat" by making all the currents go to zero, which
makes the quadratic flux vanish. To close this loophole, we can fix
the current of the first base coil::

  base_currents[0].fix_all()

(A ``Current`` object only has one degree of freedom, hence we can use
``fix_all()``.)  If you wish, you can fix the currents in all the
coils to force them to have the same value. Now the full set of 16
coils can be obtained using stellarator symmetry and field-period
symmetry::

  coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

It is illuminating to look at the non-fixed degrees of freedom that
each coil depends on. This can be done by printing the ``dof_names``
property::

  >>> print(coil[0].dof_names)

  ['CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)', 'CurveXYZFourier1:xc(1)', ...

  >>> print(coil[1].dof_names)

  ['Current2:x0', 'CurveXYZFourier2:xc(0)', 'CurveXYZFourier2:xs(1)', 'CurveXYZFourier2:xc(1)', ...

  >>> print(coil[4].dof_names)

  ['CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)', 'CurveXYZFourier1:xc(1)', ...

Notice that the current appears in the list of dofs for ``coil[1]``
but not for ``coil[0]``, since we fixed the current for
``coil[0]``. Also notice that ``coil[4]`` has the same degrees of
freedom (owned by ``CurveXYZFourier1``) as ``coil[0]``, because coils
0 and 4 refer to the same base coil shape.

There are several ways to view the objects we have created so far. One
approach is the function :obj:`simsopt.geo.plot.plot()`, which accepts
a list of Coil, Curve, and/or Surface objects::

  plot(coils + [s], engine="mayavi", close=True)

.. image:: coils_init.png
   :width: 500
	
Instead of ``"mayavi"`` you can select ``"matplotlib"`` or
``"plotly"`` as the graphics engine, although matplotlib has problems
with displaying multiple 3D objects in the proper
order. Alternatively, you can export the objects in VTK format and
open them in Paraview::

  curves = [c.curve for c in coils]
  curves_to_vtk(curves, "curves_init")
  s.to_vtk("surf_init")
  
To evaluate the magnetic field on the target surface, we create
:obj:`simsopt.field.biotsavart.BiotSavart` object based on the coils,
and instruct it to evaluate the field on the surface::

  bs = BiotSavart(coils)
  bs.set_points(s.gamma().reshape((-1, 3)))

(The surface position vector ``gamma()`` returns an array of size
``(nphi, ntheta, 3)``, which we reshaped here to
``(number_of_evaluation_points, 3)`` for the
:obj:`~simsopt.field.biotsavart.BiotSavart` object.)  Let us check the
size of the field normal to the target surface before optimization::

  B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
  print('Initial max B dot n:', np.max(B_dot_n))

The result is 0.19 Tesla. We now define the objective function::

  # Weight on the curve lengths in the objective function:
  ALPHA = 1e-6
  # Threshhold for the coil-to-coil distance penalty in the objective function:
  MIN_DIST = 0.1
  # Weight on the coil-to-coil distance penalty term in the objective function:
  BETA = 10
  
  Jf = SquaredFlux(s, bs)
  Jls = [CurveLength(c) for c in base_curves]
  Jdist = MinimumDistance(curves, MIN_DIST)
  # Scale and add terms to form the total objective function:
  objective = Jf + ALPHA * sum(Jls) + BETA * Jdist

In the last line, we have used the fact that the Optimizable objects
representing the individual terms in the objective can be scaled by a
constant and added.  (This feature applies to Optimizable objects that
have a function ``J()`` returning the objective and, if gradients are
used, a function ``dJ()`` returning the gradient.)

You can check the degrees of freedom that will be varied in the
optimization by printing the ``dof_names`` property of the objective::

  >>> print(objective.dof_names)

  ['Current2:x0', 'Current3:x0', 'Current4:x0', 'CurveXYZFourier1:xc(0)', 'CurveXYZFourier1:xs(1)', ...
   'CurveXYZFourier1:zc(5)', 'CurveXYZFourier2:xc(0)', 'CurveXYZFourier2:xs(1)', ...
   'CurveXYZFourier4:zs(5)', 'CurveXYZFourier4:zc(5)']

As desired, the Fourier amplitudes of all four base coils appear, as
do three of the four currents.  Next, to interface with scipy's
minimization routines, we write a small function::

  def fun(dofs):
    objective.x = dofs
    return objective.J(), objective.dJ()

Note that when the ``dJ()`` method of the objective is called to
compute the gradient, simsopt automatically applies the chain rule to
assemble the derivatives from the various terms in the objective, and
entries in the gradient corresponding to degrees of freedom that are
fixed (such as the current in the first coil) are automatically
removed.  We can now run the optimization using the `L-BFGS-B algorithm
from scipy
<https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html#optimize-minimize-lbfgsb>`_::

  res = minimize(fun, objective.x, jac=True, method='L-BFGS-B',
                 options={'maxiter': 200, 'iprint': 5}, tol=1e-15)
  
The optimization takes a few seconds, and the output will look like

.. code-block:: none
   
   RUNNING THE L-BFGS-B CODE

           * * *

  Machine precision = 2.220D-16
   N =          135     M =           10
   This problem is unconstrained.

  At X0         0 variables are exactly at the bounds

  At iterate    0    f=  3.26880D-02    |proj g|=  5.14674D-02

  At iterate    5    f=  6.61538D-04    |proj g|=  2.13561D-03

  At iterate   10    f=  1.13772D-04    |proj g|=  6.27872D-04

  ...
  At iterate  195    f=  1.81723D-05    |proj g|=  4.18583D-06

  At iterate  200    f=  1.81655D-05    |proj g|=  6.31030D-06

           * * *

  Tit   = total number of iterations
  Tnf   = total number of function evaluations
  Tnint = total number of segments explored during Cauchy searches
  Skip  = number of BFGS updates skipped
  Nact  = number of active bounds at final generalized Cauchy point
  Projg = norm of the final projected gradient
  F     = final function value

           * * *

   N    Tit     Tnf  Tnint  Skip  Nact     Projg        F
  135    200    234      1     0     0   6.310D-06   1.817D-05
  F =   1.8165520700970273E-005

  STOP: TOTAL NO. of ITERATIONS REACHED LIMIT                 

You can adjust parameters such as the tolerance and number of
iterations. Let us check the final :math:`\vec{B}\cdot\vec{n}` on the surface::

  B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
  print('Final max B dot n:', np.max(B_dot_n))

The final value is 0.0017 Tesla, reduced two orders of magnitude from
the initial state.  As with the initial conditions, you can plot the
optimized coil shapes directly from simsopt using

.. code-block::

  plot(coils + [s], engine="mayavi", close=True)
  
or you can export the objects in VTK format and open them in
Paraview. For this latter option, we can also export the final
:math:`\vec{B}\cdot\vec{n}` on the surface using the following
syntax::

  curves = [c.curve for c in coils]
  curves_to_vtk(curves, "curves_opt")
  s.to_vtk("surf_opt", extra_data={"B_N": B_dot_n[:, :, None]})

.. image:: coils_final.png
   :width: 500
	
The optimized value of the current in coil ``j`` can be obtained using
``coils[j].current.get_value()``. The optimized Fourier coefficients
for coil ``j`` can be obtained from ``coils[j].curve.x``, where the
meaning of each array element can be seen from
``coils[j].curve.dof_names``.  The position vector for coil ``j`` in
Cartesian coordinates can be obtained from ``coils[j].curve.gamma()``.
Magnetic fields
---------------

Simsopt provides several representations of magnetic fields, including
the Biot-Savart field associated with coils, as well as others.
Magnetic fields are represented as subclasses of the base class
:obj:`simsopt.field.magneticfield.MagneticField`.  The various field
types can be scaled and summed together; for instance you can add a
purely toroidal field to the field from coils.  Each field object has
an associated set of evaluation points.  At these points, you can
query the field B, the vector potential A, and their derivatives with
respect to position or the field parameters.


Field types
^^^^^^^^^^^

Coils and BiotSavart
~~~~~~~~~~~~~~~~~~~~

In simsopt, a filamentary coil is represented by the class
:obj:`simsopt.field.coil.Coil`. A ``Coil`` is a pairing of a
:obj:`~simsopt.geo.curve.Curve` with a current magnitude. The latter
is represented by a :obj:`simsopt.field.coil.Current` object.  For
information about Curve objects see :ref:`the page on geometric
objects <curves>`. If you wish for several filamentary curves to carry
the same current, you can reuse a :obj:`~simsopt.field.coil.Current`
object in multiple :obj:`~simsopt.field.coil.Coil` objects.

Once :obj:`~simsopt.field.coil.Coil` objects are defined, the magnetic
field they produce can be evaluated using the class
:obj:`simsopt.field.biotsavart.BiotSavart`. This class provides an
implementation of the Biot-Savart law

.. math::

  B(\mathbf{x}) = \frac{\mu_0}{4\pi} \sum_{k=1}^{n_\mathrm{coils}} I_k \int_0^1 \frac{(\Gamma_k(\phi)-\mathbf{x})\times \Gamma_k'(\phi)}{\|\Gamma_k(\phi)-\mathbf{x}\|^3} d\phi

where :math:`\mu_0=4\pi \times 10^{-7}` is the vacuum permitivity,
:math:`\Gamma_k` is the position vector of coil :math:`k`, and :math:`I_k`
indicates the electric currents.

Example::

  import numpy as np
  from simsopt.geo.curvexyzfourier import CurveXYZFourier
  from simsopt.field.coil import Current, Coil
  from simsopt.field.biotsavart import BiotSavart

  curve = CurveXYZFourier(100, 1)  # 100 = Number of quadrature points, 1 = max Fourier mode number
  curve.x = [0, 0, 1., 0., 1., 0., 0., 0., 0.]  # Set Fourier amplitudes
  coil = Coil(curve, Current(1.0e4))  # 10 kAmpere-turns
  field = BiotSavart([coil])  # Multiple coils can be included in the list 
  field.set_points(np.array([[0.5, 0.5, 0.1], [0.1, 0.1, -0.3]]))
  print(field.B())

For a more complex example of a
:obj:`~simsopt.field.biotsavart.BiotSavart` object used in coil
optimization, see
``examples/2_Intermediate/stage_two_optimization.py``.

ToroidalField
~~~~~~~~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.ToroidalField` class
represents a purely toroidal magnetic field. The field is given by
:math:`\mathbf B = B_0 \frac{R_0}{R} \mathbf e_\phi`, where
:math:`R_0` and :math:`B_0` are input scalar quantities, with
:math:`R_0` representing a reference major radius, :math:`B_0` the
magnetic field at :math:`R_0`, and :math:`R` is the radial coordinate
of the cylindrical coordinate system :math:`(R,Z,\phi)`.  The vector
:math:`\mathbf e_\phi` is a unit vector pointing in the direction of
increasing :math:`\phi`, with :math:`\phi` the standard azimuthal
angle. Given Cartesian coordinates :math:`(x,y,z)`, :math:`\mathbf e_\phi`
is calculated as :math:`\mathbf e_\phi=-\sin \phi \mathbf e_x+\cos
\phi \mathbf e_y` with :math:`\phi=\arctan(y/x)`.

PoloidalField
~~~~~~~~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.PoloidalField` class
represents a poloidal magnetic field acording to the formula
:math:`\mathbf B = B_0 \frac{r}{q R_0} \mathbf e_\theta`, where
:math:`R_0, q` and :math:`B_0` are input scalar
quantities. :math:`R_0` represents the major radius of the magnetic
axis, :math:`B_0` the magnetic field at :math:`r=R_0 q` and :math:`q`
the safety factor associated with the sum of a poloidal magnetic field
and a toroidal magnetic field with major radius :math:`R_0` and
magnetic field on-axis :math:`B_0`. :math:`r` is the radial coordinate
of the simple toroidal coordinate system
:math:`(r,\phi,\theta)`. Given a set of points :math:`(x,y,z)`,
:math:`r` is calculated as
:math:`r=\sqrt{(\sqrt{x^2+y^2}-R_0)^2+z^2}`. The vector :math:`\mathbf
e_\theta` is a unit vector pointing in the direction of increasing
:math:`\theta`, with :math:`\theta` the poloidal angle in the simple
toroidal coordinate system :math:`(r,\phi,\theta)`. Given a set of
points :math:`(x,y,z)`, :math:`\mathbf e_\theta` is calculated as
:math:`\mathbf e_\theta=-\sin \theta \cos \phi \mathbf e_x+\sin \theta
\sin \phi \mathbf e_y+\cos \theta \mathbf e_z` with
:math:`\phi=\arctan(y/x)` and
:math:`\theta=\arctan(z/(\sqrt{x^2+y^2}-R_0))`.

ScalarPotentialRZMagneticField
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The
:obj:`simsopt.field.magneticfieldclasses.ScalarPotentialRZMagneticField`
class initializes a vacuum magnetic field :math:`\mathbf B = \nabla
\Phi` defined via a scalar potential :math:`\Phi` in cylindrical
coordinates :math:`(R,Z,\phi)`. The field :math:`\Phi` is specified as
an analytical expression via a string argument. Simsopt performs the
necessary partial derivatives in order find :math:`\mathbf B` and its
derivatives. For example, the function
``ScalarPotentialRZMagneticField("2*phi")`` represents a toroidal
magnetic field :math:`\mathbf B = \nabla (2\phi)=2/R \mathbf e_\phi`.
Note: this functions needs the library ``sympy`` for the analytical
derivatives.

CircularCoil
~~~~~~~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.CircularCoil` class
represents a magnetic field created by a single circular coil. It
takes four input quantities: :math:`a`, the radius of the coil,
:math:`\mathbf c=[c_x,c_y,c_z]`, the center of the coil, :math:`I`,
the current flowing through the coil and :math:`\mathbf n`, the normal
vector to the plane of the coil centered at the coil radius, which
could be specified either with its three Cartesian components
:math:`\mathbf n=[n_x,n_y,n_z]` or as :math:`\mathbf n=[\theta,\phi]`
with the spherical angles :math:`\theta` and :math:`\phi`.

The magnetic field is calculated analitically using the following
expressions (`reference
<https://ntrs.nasa.gov/citations/20010038494>`_)

- :math:`B_x=\frac{\mu_0 I}{2\pi}\frac{x z}{\alpha^2 \beta \rho^2}\left[(a^2+r^2)E(k^2)-\alpha^2 K(k^2)\right]`
- :math:`B_y=\frac{y}{x}B_x`
- :math:`B_z=\frac{\mu_0 I}{2\pi \alpha^2 \beta}\left[(a^2-r^2)E(k^2)+\alpha^2 K(k^2)\right]`

where :math:`\rho^2=x^2+y^2`, :math:`r^2=x^2+y^2+z^2`, :math:`\alpha^2=a^2+r^2-2a\rho`, :math:`\beta^2=a^2+r^2+2 a \rho`, :math:`k^2=1-\alpha^2/\beta^2`.

Dommaschk
~~~~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.Dommaschk` class
represents a vacuum magnetic field :math:`\mathbf B = \nabla \Phi`
with basis functions for the scalar potential :math:`\Phi` described
in `W. Dommaschk (1986), Computer Physics Communications 40, 203-218
<https://www.sciencedirect.com/science/article/pii/0010465586901098>`_. This
representation provides explicit analytic formulae for vacuum fields
with a mixture of flux surfaces, islands, and chaos. Following the
original reference, a toroidal field with :math:`B_0=R_0=1` is already
included in the definition. As input parameters, it takes two arrays:

- The first array is an :math:`N\times2` array :math:`[(m_1,n_1),(m_2,n_2),...]` specifying which harmonic coefficients :math:`m` and :math:`n` are non-zero.
- The second array is an :math:`N\times2` array :math:`[(b_1,c_1),(b_2,c_2),...]` with :math:`b_i=b_{m_i,n_i}` and :math:`c_i=c_{m_i,n_i}` the coefficients used in the Dommaschk representation.

Reiman
~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.Reiman` provides the
magnetic field model in section 5 of `Reiman and Greenside, Computer
Physics Communications 43 (1986) 157—167
<https://www.sciencedirect.com/science/article/pii/0010465586900597>`_.
It is an analytical magnetic field representation that allows the
explicit calculation of the width of the magnetic field islands.

InterpolatedField
~~~~~~~~~~~~~~~~~

The :obj:`simsopt.field.magneticfieldclasses.InterpolatedField`
function takes an existing field and interpolates it on a regular grid
in :math:`r,\phi,z`. This resulting interpolant can then be evaluated
very quickly. This is useful for efficiently tracing field lines and
particle trajectories.

Scaling and summing fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

Magnetic field objects can be added together, either by using the
``+`` operator, or by creating an instance of the class
:obj:`simsopt.field.magneticfield.MagneticFieldSum`. (The ``+``
operator creates the latter.)

Magnetic fields can also be scaled by a constant. This can be accomplished either using the ``*`` operator,
or by creating an instance of the class
:obj:`simsopt.field.magneticfield.MagneticFieldMultiply`. (The ``*``
operator creates the latter.)

Example::

   from simsopt.field.magneticfieldclasses import ToroidalField, CircularCoil
   
   field1 = CircularCoil(I=1.e7, r0=1.)
   field2 = ToroidalField(R0=1., B0=1.)
   total_field = field1 + 2.5 * field2

Common operations
^^^^^^^^^^^^^^^^^

Magnetic field objects have a large number of functions available. Before evaluating the field, you must
set the evaluation points. This can be done using either Cartesian or cylindrical coordinates.
Let ``m`` be a :obj:`~simsopt.field.magneticfield.MagneticField` object, and suppose there are ``n`` points
at which you wish to evaluate the field.

- ``m.set_points_cart()`` takes a numpy array of size ``(n, 3)`` with the Cartesian coordinates ``(x, y, z)`` of the points.
- ``m.set_points_cyl()`` takes a numpy array of size ``(n, 3)`` with the cylindrical coordinates ``(r, phi, z)`` of the points.
- ``m.set_points()`` is shorthand for ``m.set_points_cart()``.
- ``m.get_points_cart()`` returns a numpy array of size ``(n, 3)`` with the Cartesian coordinates ``(x, y, z)`` of the points.
- ``m.get_points_cyl()`` returns a numpy array of size ``(n, 3)`` with the cylindrical coordinates ``(r, phi, z)`` of the points.

A variety of functions are available to return the magnetic field
:math:`B`, vector potential :math:`A`, and their gradients.  The most
commonly used ones are the following:

- ``m.B()`` returns an array of size ``(n, 3)`` with the Cartesian coordinates of :math:`B`.
- ``m.B_cyl()`` returns an array of size ``(n, 3)`` with the cylindrical ``(r, phi, z)`` coordinates of :math:`B`.
- ``m.A()`` returns an array of size ``(n, 3)`` with the Cartesian coordinates of :math:`A`.
- ``m.AbsB()`` returns an array of size ``(n, 1)`` with the field magnitude :math:`|B|`.
- ``m.dB_by_dX()`` returns an array of size ``(n, 3, 3)`` with the Cartesian coordinates of :math:`\nabla B`. Denoting the indices
  by :math:`(i,j,l)`, the result contains  :math:`\partial_j B_l(x_i)`.
- ``m.d2B_by_dXdX()`` returns an array of size ``(n, 3, 3, 3)`` with the Cartesian coordinates of :math:`\nabla\nabla B`. Denoting the indices
  by :math:`(i,j,k,l)`, the result contains  :math:`\partial_k \partial_j B_l(x_i)`.
- ``m.dA_by_dX()`` returns an array of size ``(n, 3, 3)`` with the Cartesian coordinates of :math:`\nabla A`. Denoting the indices
  by :math:`(i,j,l)`, the result contains  :math:`\partial_j A_l(x_i)`.
- ``m.d2A_by_dXdX()`` returns an array of size ``(n, 3, 3, 3)`` with the Cartesian coordinates of :math:`\nabla\nabla A`. Denoting the indices
  by :math:`(i,j,k,l)`, the result contains  :math:`\partial_k \partial_j A_l(x_i)`.
- ``m.GradAbsB()`` returns an array of size ``(n, 3)`` with the Cartesian components of :math:`\nabla |B|`.

Example:

.. code-block::

   import numpy as np
   from simsopt.field.magneticfieldclasses import CircularCoil
   
   field = CircularCoil(I=1.e7, r0=1.)
   points = np.array([[0.5, 0.5, 0.1], [0.1, 0.1, -0.3]])
   field.set_points(points)
   print(field.B())
   print(field.dB_by_dX())

Containers
**********

.. _docker_doc:

Docker container
================

A Docker container for simsopt is available, allowing you to use
simsopt without having to compile any code yourself.  The container
includes VMEC, SPEC, and BOOZ_XFORM.

.. warning::

   Docker is not generally allowed to run on computers at HPC centers due to security issues.
   For those wishing to run simsopt on NERSC machines, please refer to :ref:`shifter_doc`.

Requirements
^^^^^^^^^^^^
Docker needs to be installed before running the docker container. Docker
software can be obtained from `docker website <https://docs.docker.com/get-docker/>`_.
Check the `docker get started webpage <https://docs.docker.com/get-started/>`_ for installation instructions 
as well as for tutorials to get a feel for docker containers. On linux, you may need to start the docker daemon
before proceeding further.

.. warning::

   On Mac, the default 2 GB memory per container assigned by Docker Desktop is not sufficient. Increase the memory of
   the container to at least 3 GB to run simsopt much faster.

Install From Docker Hub
^^^^^^^^^^^^^^^^^^^^^^^
The easiest way to get simsopt docker image which comes with simsopt and all of its dependencies such as
SPEC and VMEC pre-installed is to use Docker Hub. After 
`installing docker <https://docs.docker.com/get-started/>`_, you can run
the simsopt container directly from the simsopt docker image uploaded to
Docker Hub.

.. code-block::

   docker run -it --rm hiddensymmetries/simsopt python # Linux users, prefix the command with sudo

The above command should load the python shell that comes with the simsopt
docker container. When you run it first time, the image is downloaded
automatically, so be patient.  You should now be able to import the module from
python::

  >>> import simsopt

Ways to use simsopt docker container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IPython Shell
-------------

Easiest way is to start ipython shell and import the simsopt
library. But this approach is only useful if a few commands need to be
executed or to load a python module and execute it.

.. code-block::

    docker run -it --rm hiddensymmetries/simsopt ipython

Bash Shell
----------

In this approach, you write a simsopt based driver script for your optimization problem. One
needs to mount the host directory in the container before running scripts. Use the ``-v`` flag 
to mount the current directory

.. code-block:: 

    docker run -it --rm -v $PWD:/my_mount hiddensymmetries/simsopt 
    <container ###> cd /my_mount
    <container ###> python <driver_script>

Jupyter notebook
----------------

The simsopt container comes with jupyter notebook preloaded. You can launch the jupyter from
the container using the command:

.. code-block::
   
    docker run -it --rm -v $PWD:/my_mount -p 8888:8888 hiddensymmetries/simsopt 
    <container ###> cd /my_mount
    <container ###> jupyter notebook --ip 0.0.0.0 --no-browser --allow-root 

Running the above command, a link will be printed that you can use to
open the jupyter console in a web browser. The link typically starts
with ``http://127.0.0.1:8888/?token=``. (Several links are printed,
but only the last one listed will work for browsers running outside
the container.) Copy the full link and paste it into any browser in
your computer.


Persistent containers
^^^^^^^^^^^^^^^^^^^^^

Using the intructions above will create a fresh container each time and delete the container after exiting.
If you would like to create a persistent container (e.g. because you are installing additional pip packages inside) that you can reuse at any time,
you can do so by removing the ``--rm`` command and specifying a container name via ``--name=``

.. code-block::

    docker run --name=mycontainer -it -v $PWD:/my_mount hiddensymmetries/simsopt
    <container ###> cd /my_mount
    <container ###> python <driver_script>

And to restart and rejoin the container:

.. code-block::

    docker start mycontainer
    docker exec -it mycontainer /bin/bash
    <container ###> source /venv/bin/activate



.. _shifter_doc:

Shifter container
=================

`Shifter <https://docs.nersc.gov/development/shifter/>`_ is the
container technology deployed at NERSC to circumvent the security
issues associated with Docker containers. Shifter allows to you use
the simsopt Docker image files hosted on Docker Hub.  Detailed
instructions for using Shifter can be found at the `NERSC page on the
simsopt wiki
<https://github.com/hiddenSymmetries/simsopt/wiki/NERSC-Cori>`_.

Testing
^^^^^^^

``simsopt`` includes unit and regression tests, and continuous integration.

Python test suite
*****************

The main test suite is based on the standard ``unittest`` python module.
Source code for the python tests is located in the ``tests`` directory.
These tests will use the installed version of the ``simsopt`` python package,
which may differ from the code in your local repository if you did not
make an editable install (see :doc:`installation`).

To run all of the tests in the test suite on one processor, you can type

.. code-block::

    ./run_tests

from the command line in the repository's home directory. Equivalently,
you can run

.. code-block::

    python -m unittest

from the ``tests`` directory.

For some of the tests involving MPI, it is useful to execute the tests
for various numbers of processes to make sure the tests pass in each
case. For this purpose, it is not necessary to run the entire test
suite, only the tests for which MPI is involved.  For convenience, you
can run the script

.. code-block::

    ./run_tests_mpi

in the repository's home directory. This script runs only the tests
that have ``mpi`` in the name, and the tests are run on 1, 2, and 3
processors.

The tests make use of data files in the ``tests/test_files`` directory.


Longer examples
***************

For convenience, the main test suite is designed to run in no more than a few minutes.
This means that some more complicated integrated and regression tests that require substantial time
are not included. You may wish to run some of these more complicated tests by hand during development.
A number of such examples can be found in the ``examples`` subdirectory.
Also, ``simsopt`` and ``stellopt`` have been benchmarked for several problems in the
`stellopt_scenarios collection <https://github.com/landreman/stellopt_scenarios>`_,
which includes the corresponding ``simsopt`` input files.


Continuous integration
**********************

The serial and MPI tests are automatically run after every commit to
the repository.  This automation is handled by GitHub Actions, and
controlled by the script ``.github/workflows/ci.yml``.
To view the results of the continuous integration runs, you can click on the "Actions"
link from the `GitHub repository page <https://github.com/hiddenSymmetries/simsopt>`_,
or you can directly visit `<https://github.com/hiddenSymmetries/simsopt/actions>`_.
simsopt.mhd package
===================

Submodules
----------

simsopt.mhd.boozer module
-------------------------

.. automodule:: simsopt.mhd.boozer
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.mhd.spec module
-----------------------

.. automodule:: simsopt.mhd.spec
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.mhd.vmec module
-----------------------

.. automodule:: simsopt.mhd.vmec
   :members:
   :undoc-members:
   :show-inheritance:

simsopt.mhd.vmec\_diagnostics module
------------------------------------

.. automodule:: simsopt.mhd.vmec_diagnostics
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: simsopt.mhd
   :members:
   :undoc-members:
   :show-inheritance:
