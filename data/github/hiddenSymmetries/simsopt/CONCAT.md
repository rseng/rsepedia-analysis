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
container](https://simsopt.readthedocs.io/en/latest/docker.html), and
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



