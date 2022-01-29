# Changelog

## 0.1.0 
- Initial release
# NLMech 

[![CircleCI](https://circleci.com/gh/nonlocalmodels/NLMech.svg?style=shield)](https://circleci.com/gh/nonlocalmodels/NLMech) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/118379d7d745464584b73e9e06f60462)](https://www.codacy.com/gh/nonlocalmodels/NLMech?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=nonlocalmodels/NLMech&amp;utm_campaign=Badge_Grade) [![Coverage Status](https://coveralls.io/repos/github/nonlocalmodels/NLMech/badge.svg?branch=main)](https://coveralls.io/github/nonlocalmodels/NLMech?branch=main) [![GitHub issues](https://img.shields.io/github/issues/nonlocalmodels/nlmech.svg)](https://github.com/nonlocalmodels/NLMech/issues) ![GitHub release](https://img.shields.io/github/release/nonlocalmodels/NLMech.svg) [![GitHub license](https://img.shields.io/github/license/nonlocalmodels/nonlocalmodels.github.io.svg)](https://github.com/nonlocalmodels/nonlocalmodels.github.io/blob/main/LICENSE) [![status](https://joss.theoj.org/papers/271dd66ea91b7fbfdccb4b10a7ba462c/status.svg)](https://joss.theoj.org/papers/271dd66ea91b7fbfdccb4b10a7ba462c) [![DOI](https://zenodo.org/badge/178625941.svg)](https://zenodo.org/badge/latestdoi/178625941)

 

<p style="text-align:center;"><img src="https://github.com/nonlocalmodels/NLMech/blob/main/assets/logo/logo_sim.png?raw=true" alt="logo" width="400"/></p>

## Introduction
Welcome to the NLMech repository. In this project, we implement Peridynamics model of fracture using meshfree and finite element discretizations. A brief overview of the equations is available [here](https://nonlocalmodels.github.io/documentation/md_content_equations.html). NLMech primarily served as a code for academic research (e.g., [1,2]), however, we plan to improve it further for a large-scale usage. The plan also includes development of fully distributed solver using HPX library for asynchronous computation to its maximum potential. In [3], we discuss the structure of NLMech and use HPX library for multi-threading computation. For further list of publications based on this library, we refer to the [publication list](https://nonlocalmodels.github.io/publications/).

At present, the library consists of multiple material models such as 
- **RNP** - Regularized Nonlinear Potential. This is implemented in class [RNPBond](src/material/pd/rnpBond.h).
- **State-based peridynamics** - State-based peridynamics model. This is implemented in class [ElasticState](src/material/pd/ElasticState.h).

One of the main features of NLMech is the implementation of both explicit time discretization for dynamic problems (see [FDModel](src/model/fd/fDModel.h)) and implicit discretization for quasi-static problems (see [QuasiStaticModel](src/model/quasistatic/QuasiStaticModel.h)). 

## Documentation and getting started
All source and header files are fairly well documented. We used doxygen to automatically generate the documentation of methods, classes, etc. For complete documentation follow this [link](https://nonlocalmodels.github.io/documentation/).

We provide shell scripts to help with the installation of dependencies and the NLMech itself. We also provide Docker images for a quick test of the library and to run the examples. In section `Installation`, we describe the dependencies, installation of dependencies, and building NLMech code. In section `Running NLMech`, we discuss how to run the examples.


## Installation

### Build tools
The following build tools are needed to compile the NLMech and its dependencies:
  * GCC compiler collection (gcc) > 4.9, however, gcc >= 8 is recommended
  * [autoconf](https://www.gnu.org/software/autoconf/)
  * [wget](https://www.gnu.org/software/wget/)
  * [cmake](https://cmake.org/)
  * [git](https://git-scm.com/)

   On Ubuntu you might install all dependencies using the papackage manager:

  ```bash
  apt-get install build-essential git wget cmake libssl-dev libblas-dev liblapack-dev autoconf freeglut3-dev

  ```

  On Fedora you might install all dependencies  using the package manager

  ```bash
  dnf install @development-tools cmake git wget blas-devel lapack-devel freeglut-devel
  ```

### Dependencies
We use cmake to build the code. We list the dependencies and how they are used in the code below:
  * [CMake](https://cmake.org/) 3.16
    - To generate makefiles
  * [Boost](https://www.boost.org/) 1.75
    - Required to build HPX and YAML libraries
  * [HPX](https://github.com/STEllAR-GROUP/hpx) 1.6.0
    - Provides methods for multithreading computation
  * [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/) 3.8
    - Required to build the BlazeIterative library
  * [Blaze_Iterative](https://github.com/STEllAR-GROUP/BlazeIterative) master
    - Provides linear algebra support such as storage and inversion of stiffness matrix
  * [gmsh](https://gmsh.info/) 4.7
    - Our code directly interfaces with gmsh to use input-ouput functions of gmsh
  * [VTK](https://www.vtk.org) 9.0
    - For read-write operations on `.vtu` type files
  * [YAML-CPP](https://github.com/jbeder/yaml-cpp) 0.6
    - To parse `.yaml` input files

  On Ubuntu you might install all dependencies using the papackage manager:

  ```bash
  apt-get install libyaml-cpp-dev libvtk7-dev gmsh boost-devel
  ```

  Note that on Ubuntu you need to install HPX, Blaze, and Blaze_Iterative since
  there is no package available.

  On Fedora you might install all dependencies  using the package manager

  ```bash
  dnf install hpx-devel cmake blaze-devel vtk-devel yaml-cpp-devel gmsh-devel
  ```

  Note on Fedora, you need only to install Blaze_Itertaive.

Following dependencies are optional, but are recommended for the large simulations:
  * [PCL](https://github.com/PointCloudLibrary/pcl) 1.11 
    - Used in neighbor search calculation (KDTree)
  * [Flann](https://github.com/flann-lib/flann)  1.9
    - Required to build PCL library

### Building dependencies
Building above dependencies is quite a challenge. To help with this, we provide the bash scripts for Ubuntu-20.04 and Fedor operating systems: [Bash scripts](https://github.com/nonlocalmodels/buildscripts/tree/main/bash)).

Further, we provide various docker files
* to build the code on Fedora using Fedora packages, see [Using Fedora Packages](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/Fedora) 
* to build all the dependencies and the code on Fedora, see [Build Dependencies and Code on Fedora](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/FedoraAll)
* to build dependencies and the code on Ubuntu-20.04, see [Build Dependencies and Code on Ubuntu](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/FedoraAll).

In [Scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure), bash scripts to build individual dependencies such as blaze, vtk, hpx, etc, on HPC systems is provided.

A more detailed version of the build instruction is available [here](https://nonlocalmodels.github.io/documentation/md_content_install_instructions.html).

> :exclamation: We recommend to use the same CMake version and same compiler version to build the HPX and NLMech.

### Compiling library
Assuming all the dependencies are installed at standard paths such as `/usr/local`, we build NLMech as follows
```sh
git clone https://github.com/nonlocalmodels/NLMech.git 
cd NLMech
mkdir build && cd build 

# turn building documentation and tools, dependency on PCL off by using 'OFF' instead of 'ON'
cmake -DCMAKE_BUILD_TYPE=Release \
      -DEnable_Documentation=ON \
      -DEnable_Tools=ON \
      -DEnable_PCL=ON \
      .. 
      
# compile the code
make -j $(cat /proc/cpuinfo | grep processor | wc -l) VERBOSE=1
```

In case certain libraries such as `HPX, PCL, Blaze_Iterative` are installed at custom paths, one would need to provide correct paths to `cmake` as follows:
```sh
cmake -DCMAKE_BUILD_TYPE=Release \
      -DEnable_Documentation=ON \
      -DEnable_Tools=ON \
      -DEnable_PCL=ON \
      -DHPX_DIR=<hpx-path>/lib/cmake/HPX \
      -DPCL_DIR=<pcl-path> \
      -DYAML_CPP_DIR=<yaml-cpp-path> \
      -DBLAZEITERATIVE_DIR=<blaze-iterative path> \
      -DGMSH_DIR=<gmsh-path> \
      .. 
      
make -j $(cat /proc/cpuinfo | grep processor | wc -l) VERBOSE=1
```

## Running NLMech
To quickly run the tests and examples, you may use Docker image with the [latest Successful build](https://hub.docker.com/r/diehlpk/nlmech/tags?page=1&ordering=last_updated) of the main branch. 

```sh
podman/docker pull diehlpk/nlmech:latest
podman/docker run -it docker.io/diehlpk/nlmech /bin/bash
cd /app/NLMech/examples/qsModel/1D
# Generate the mesh file
/app/NLMech/build/bin/mesh -i input_mesh.yaml -d 1
# Run the simulation
/app/NLMech/build/bin/NLMech -i input.yaml --hpx:threads=2
```

In [examples](https://nonlocalmodels.github.io/examples/), we provide information on how to prepare a simulation setup input file using [YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html).

Asume, you have build NLMech on your own, you can go the the `build` folder and run the executable as below

```sh
cd build
./bin/NLMech -i input.yaml --hpx:threads=n
```

with the first argument `-i` the `input.yaml` file is specified and the second argument `--hpx:threads` the amount
of CPU cores HPX is allowed to use is specified. If you do not specify any number there all coes of the CPU are used
to run the simulation. Note that in the current version only shared memory parallism is provided, however, we plan to add
dsitributed memory to the code in the near future. 

The one-dimensional quasi-static example is computational inexpesive, therfore, we used it in the Docker example to finish the
simulation soon. For scaling test, we recommend to use any of the two-dimensional examples. 

## Trouble, issues, bugs
In case you found a bug in the library, want to contribute, or need a feature, please create a new [issue](https://github.com/nonlocalmodels/NLMech/issues). 


## Releases
The current stable version is [![GitHub release](https://img.shields.io/github/release/nonlocalmodels/NLMech.svg)](https://GitHub.com/nonlocalmodels/NLMech/releases/). Main development branch is the [main branch](https://github.com/nonlocalmodels/NLMech). For more details, we refer to the [Changelog](https://github.com/nonlocalmodels/NLMech/blob/main/CHANGELOG.md) file.

## Code of conduct
We have adopted a [code of conduct](https://github.com/nonlocalmodels/NLMech/blob/main/CODE_OF_CONDUCT.md) for this project. 

## Contributing
The source code is released under the [![GitHub license](https://img.shields.io/github/license/nonlocalmodels/nonlocalmodels.github.io.svg)](https://github.com/nonlocalmodels/nonlocalmodels.github.io/blob/main/LICENSE) license. If you like to contribute, we only accept your pull request using the same license. Please feel free to add your name to license header of the files you added or contributed to. If possible please add a test for your new feature using [CTest](https://gitlab.kitware.com/cmake/community/-/wikis/doc/ctest/Testing-With-CTest). We adapted the Google C++ [Style Guide](https://google.github.io/styleguide/cppguide.html) for this project. We use the [clang-format](https://clang.llvm.org/docs/ClangFormat.html) tool to format the source code with respect to this style guide. Please run the `format.sh` script before your do the pull request.

## Citing
In publications, please use our paper as the main citation for NLMech: 
* Diehl, P., Jha, P. K., Kaiser, H., Lipton, R., & Lévesque, M. (2020). An asynchronous and task-based implementation of peridynamics utilizing HPX—the C++ standard library for parallelism and concurrency. SN Applied Sciences, 2(12), 1-21.
In addition, please use our [JOSS](https://joss.theoj.org/) for the referencing the code

* Jha et al., (2021). NLMech: Implementation of finite difference/meshfree discretization of nonlocal fracture models. Journal of Open Source Software, 6(65), 3020, [10.21105/joss.03020](https://doi.org/10.21105/joss.03020)

For more references, we refer to NLMech's [publication list](https://nonlocalmodels.github.io/publications/).

## Acknowledgments
NLMech has been funded by:
* Army Research Office Grant # W911NF-16-1-0456 to PI Dr. Robert Lipton (Professor at Louisiana State University). This grant supported Prashant K. Jha on a postdoctoral position from October 2016 - July 2019.
*  Canada Research Chairs Program under the Canada Research Chair in Multiscale Modelling of Advanced Aerospace Materials held by M. Lévesque; Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants Program under Discovery Grant RGPIN-2016-06412.
* We are grateful for the support of the [Google Summer of Code program](https://summerofcode.withgoogle.com/) funding internships.

## References
[1] Jha, P. K., & Lipton, R. (2019). Numerical convergence of finite difference approximations for state based peridynamic fracture models. Computer Methods in Applied Mechanics and Engineering, 351, 184-225.

[2] Jha, P. K., & Lipton, R. P. (2020). Kinetic relations and local energy balance for LEFM from a nonlocal peridynamic model. International Journal of Fracture, 226(1), 81-95.

[3] Diehl, P., Jha, P. K., Kaiser, H., Lipton, R., & Lévesque, M. (2020). An asynchronous and task-based implementation of peridynamics utilizing HPX—the C++ standard library for parallelism and concurrency. SN Applied Sciences, 2(12), 1-21.
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
reported by contacting the project team at nlm-fellows@cct.lsu.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
## Expected Behavior

... Please describe the behavior you would have expected.

## Actual Behavior

... Please describe the behavior you actually observed.


## Steps to Reproduce the Problem

... Please be as specific as possible while describing how to reproduce your problem.

  1.
  1.
  1.

## Specifications

... Please describe your environment

  - HPX Version:
  - Platform (compiler, OS):


Fixes #

## Proposed Changes

  -
  -
  -

## Any background context you want to provide?
### Crack propagation in Glass material

Fore more details, see [here](https://nonlocalmodels.github.io/examples/fd-crack-glass-material.html)


### Example files for restarting a simualiton

Fore more details, see [here](https://nonlocalmodels.github.io/examples/restart.html)


### Example files for the code's logo

Fore more details, see [here](https://nonlocalmodels.github.io/examples/fd-logo-soft-material-2.html)


### Example files for the code's logo

Fore more details, see [here](https://nonlocalmodels.github.io/examples/fd-logo-soft-material.html)


### Postprocessing simulation results
We provide various methods in [Compute](../../../../../tools/pp/src/compute.cpp) to perform postprocessing on simulation results. In this example, we will compute damage at nodes and strain and stress.

We assume that input `.yaml` file used in obtaining simulatio is available and is in the same path where simulation results `.vtu` files are.

### Input file for postprocessing

1. Provide path where simulation results and input file is located.

1. Specifiy the prefix the simulation filenames have. E.g. if files are `output_1.vtu, output_2.vtu,...` then prefix is `output`.

1. If simulation input file does not have elastic material properties, then provide them.

1. Specify output path.

1. We can specify multiply compute sets. In this example we have two compute sets. One for damage, and another for strain and stress. Thus `Sets: 2`.

1. For each set give the required instruction. 

#### Damage 
Damage calculation is performed in `Set_1`.  

- Provide filename for this set. We have `Tag_Filename: damage_Z`.

- Specify (optional) start and end index of simulation file, and specify (optional) interval between processing each file. 

	- For example, we want to skip first 10 files, we will have `Dt_Start: 11`. 

	- If we want to skip 1 file each time, i.e. consider `output_11.vtu, output_13.vtu, ...` then we will have `Dt_Interval: 2`.

- Set `Damage_Z: true` which is responsible for computing damage. 

#### Strain and stress
This calculation is performed in `Set_2`.  

- Description of some fields are same as before. However, we can have different values for them in different sets. 
	
	- We take `Tag_Filename: strain`, `Dt_Start: 20`, `Dt_Interval: 4`.

- `Output_Only_Nodes: true` means we do not write element-node connectivity in the output `.vtu` file. As strain/stress are cell data, we set it to `false`.

- Set `Strain_Stress: true` which is responsible for computing strain and stress. 

- Additionally, we can compute the magnitude of strain tensor. 

	- If we want to compute magnitude of tensor, we only need to have 'Magnitude_Strain_Tensor:' with no entry. 

	- If we want to compute magnitude of `yy` component, we can specify `Component: yy`. 
### Mesh
We describe the input file required to run [Mesh](../../../../../tools/mesh/mesh.cpp). 

- We consider rectangle domain of size `0.1 m x 0.1 m`. 

- Horizon is specified to be `horizon = 0.02 m`. 

- Ratio of horizon to mesh size `r = 4`. Thus mesh size is `horizon / 4`.

- Type of mesh is `uniform_tri` which produces uniform mesh of triangles. Other options are `uniform_square` and `uniform_tri_sym`.

### Results
Mesh looks like

<p id="result" align="center">
	<img src="result.png" alt="setup" width="400" height="400" />
</p>
---
title: 'NLMech: Implementation of finite difference/meshfree discretization of nonlocal fracture models'
tags:
  - Peridynamics
  - Finite difference
  - Finite element
  - HPX
  - Asynchronous many-task systems
authors:
  - name: Prashant K. Jha^[co-first author]
    orcid: 0000-0003-2158-364X
    affiliation: 1
  - name: Patrick Diehl^[co-first author]
    orcid: 0000-0003-0872-7098
    affiliation: 2
affiliations:
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, TX, United States of America
   index: 1
 - name: Center for Computation \& Technology, Louisiana State University, Baton Rouge, LA, United States of America
   index: 2
date: 25 November 2020
bibliography: paper.bib

---

![NLMech's logo which shows the obtained damage of a peridynamic simulation.\label{fig:logo}](../assets/logo/logo_joss.png){ width=40% }

# Summary

The open source code *NLMech* is an implementation of finite difference approximation of nonlocal models, \emph{e.g.}\ peridynamic. Peridynamic (PD) [@silling2007peridynamic;@silling2005meshfree] is a nonlocal formulation of classical continuum mechanics that is particularly robust in mechanical deformations involving crack (discontinuous displacement) and damage. The model seamlessly handles the two regimes of deformation: elastic/smooth deformation and fracture. The constitutive laws describing the material behavior are simple to conceptualize and implement. Particularly, in numerical implementation, no special care for the modeling of cracks is required. Successful comparison of PD against a variety of experiments has been done [@diehl2019review,diehl2021comparative]. 

Unlike classical continuum mechanics, where the internal force is written in terms of the stress, in PD, the internal force at a given material point is due to the sum of the pairwise forces of the neighboring points. \emph{i.e.}\ the force in PD is expressed as the integral of the pairwise force density between the given point and another point in the neighborhood. The neighborhood of point $x$ is typically defined as all points in the sphere of radius $\delta$, centered at $x$, where $\delta$ is the nonlocal length scale and is referred to as the \textit{horizon}. PD is often divided in two classes: bond-based and state-based models. In bond-based models, the two material points interact via a pairwise force law, and the forces between the material points do not depend on the deformation state of surrounding points. In contrast, in the state-based models the volumetric deformation in the neighborhood of two points plays a role in the pairwise force. The governing equation of motion for the bond-based PD [@silling2005meshfree] reads as

$$ \varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) = \int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X}) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D$$

and the governing equation for the state-based PD [@silling2007peridynamic] reads as 

$$  \varrho (\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) =  \int\limits_{B_\delta(\mathbf{X})} (T[\mathbf{X},t]\langle \mathbf{X}' - \mathbf{X} \rangle - T[\mathbf{X}',t]\langle \mathbf{X} - \mathbf{X}' \rangle) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D \text{.} $$
Here $\varrho$ denotes density of the material, $\mathbf{u}$ displacement field in the material, $\ddot{\mathbf{u}}$ acceleration, and $\mathbf{b}$ external force density. The constitutive law, relating bond strain with bond force, is prescribed using either the pairwise force function $\mathbf{f}$ or the PD state $T$ [@silling2007peridynamic]. In the NLMech library, the following material models are implemented:

* Elastic state-based PD model [@silling2007peridynamic],
* Prototype micro-elastic brittle bond-based PD model [@silling2005meshfree],
* Nonlinear bond-based PD model [@lipton2014dynamic;@lipton2016cohesive], and
* Nonlocal double-well state-based peridynamic model [@Lipton2018;@jha2019numerical].

Examples for these types of model implementations are provided in the [documentation](https://nonlocalmodels.github.io/examples/).

Currently, the library supports finite difference (or more generally meshfree) discretization. Using the triangulation of an arbitrary domain, the library can create a meshfree discretization. The library is equipped with necessary modules, such as FE elements and quadrature integration rules, for finite element discretization of PD. Next, we briefly discuss the finite difference/meshfree discretization of PD. \autoref{fig:discrete} shows the domain D discretized with the nodes $X = \{ X_i \in \mathbb{R}^3 \vert i=1,\ldots,n\}$. Each node $X_i$ represents a small area/volume denoted by $V_i$. In PD, as previously mentioned, each point $X_i$ interacts with neighboring points in the sphere (discrete)  $B_\delta(X_i) = \{X_j: |X_i - X_j| < \delta \}$. 

![ Adapted from [@Diehl2020].\label{fig:discrete}](discrete.pdf){ width=30% }

The discrete equation of motion is written as, for the bond-based PD,

$$ \varrho(X_i)\ddot{\mathbf{u}}(t,X_i) = \sum\limits_{j \in B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i) V_j + \mathbf{b}(t,X_i) \text{ in } D,$$

and, the state-based PD,

$$  \varrho (X_i)\ddot{\mathbf{u}}(t,X_i) =  \sum\limits_{j \in B_\delta(X_i)} (T[X_i,t]\langle X_j - X_i \rangle - T[X_j,t]\langle X_i - X_j \rangle) V_j + \mathbf{b}(t,X_i) \text{ in } D \text{.} $$

Here $\mathbf{u}(t,X_i)$ denotes the displacement of node $X_i$ at time $0 \leq t\leq T$. For the time discretization, we can consider: \textit{1)} implicit time integration and \textit{2)} explicit time integration using either a central difference or velocity verlet scheme.

## Software Implementation and Applications

NLMech relies on the following open source softwares: HPX [@Kaiser2020], Blaze [@iglberger2012high], Blaze_Iterative, Gmsh [@geuzaine2009gmsh], VTK [@schroeder2004visualization], and yaml-cpp. For details 
about the specific version, we refer to NLMech's [documentation](https://github.com/nonlocalmodels/NLMech#building).

NLMech was used for the following applications/publications:

* Numerical convergence of finite difference approximations for state based perdidynamic fracture
models [@jha2019numerical] 
* Complex fracture nucleation and evolution with nonlocal elastodynamics [@lipton2019complex]
* Free damage propagation with memory [@lipton2018free] 
* Kinetic relations and local energy balance for linear elastic fracture mechanics from a
nonlocal peridynamic model [@jha2020kinetic]
* A Fracture Multiscale Model for Peridynamic enrichment within the Partition of Unity Method [@DBLP:journals/corr/abs-2108-02336]
* Peridynamics for Quasistatic Fracture Modeling [@bhattacharya2021peridynamics]

For an updated list of applications/publications, we refer to corresponding [NLMech documentation](https://nonlocalmodels.github.io/publications/).

# Statement of need

Nonlocal models, like peridynamic, are computationally expensive. Several 
publications on GPU-based implementations [@mossaiby2017opencl;@diehl2012implementierung;@diehl2015efficient] and one commercial implementation in LS-DYNA [@ren20173d] can be found in literature. However, 
from an open source perspective, only two other peridynamic implementations are available: [Peridigm](https://github.com/peridigm/peridigm) [@littlewood2015roadmap] and [PDLammps](https://lammps.sandia.gov/doc/pair_peri.html) [@parks2008implementing]. Both of these codes rely on the Message Passing Interface (MPI). On modern supercomputers, many core architectures where the threads per computational node increase, it is more and more important to focus on the fine-grain parallelism with increasing cores per computational nodes. NLMech is based on the C++ standard library for parallelism and concurrency (HPX) [@Kaiser2020]. For more details on utilization of asynchronous many-task systems, we refer to @diehl2018implementation. The library implements the experimental nonlinear bond-based and state-based models, and the process of adding new material models is simple following the existing templates.

# Future directions

We are interested in extending/improving the library with

- implementation of new material models,
- higher order time discretization schemes,
- local-nonlocal coupling methods, and
- further optimization of nonlocal computation. 
 
If you are interested in contributing, please read our [guidelines](https://github.com/nonlocalmodels/NLMech#contributing) and our [code of conduct](https://github.com/nonlocalmodels/NLMech/blob/master/CODE_OF_CONDUCT.md) before doing so. 

# Acknowledgments

NLMech has been funded by:

*  Army Research Office Grant # W911NF-16-1-0456 to PI Dr. Robert Lipton (Professor at Louisiana State University). This grant supported Prashant K. Jha on a postdoctoral position from October 2016 - July 2019.
*  Canada Research Chairs Program under the Canada Research Chair in Multiscale Modelling of Advanced Aerospace Materials held by M. Lévesque; Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants Program under Discovery Grant RGPIN-2016-06412.
* We are grateful for the support of the Google Summer of Code program funding internships.

For a updated list of previous and current funding, we refer to the corresponding [NLMech website](https://github.com/nonlocalmodels/NLMech#acknowledgements).

# References
# CMake options

## NLMech 

* Enable_Documentation : Generates target for generating the documentation (Default = False)
* Enable_Tools : Enables the tools to the build target (Default = False)
* Enable_Expensive_Tests : Enables the computation intense tests (Default = False)
* Enable_RPM: Enables to generate RPM packages (Default = False)
* Enable_PCL: Add [PCL](https://github.com/PointCloudLibrary/pcl) for faster neighborseach. Recommend for large amount of nodes (Default=False)

## General options

* CMAKE_BUILD_TYPE : Specifies the build type (Default = Release)
* CMAKE_INSTALL_PREFIX : Specifies the install directory (Default = /usr/local)

## Dependencies

* VTK_DIR : Path the the VTK installation
* BLAZEITERATIVE_INCLUDE : Path to the Blaze_Iterative installation 
* GMSH_DIR : Path to the Gmsh installation 
* YAML_CPP_DIR : Path to the YAML CPP installation 
* HPX_DIR : Path to the HPX installation 

Note that HPX's CMake will load the same Boost version as HPX was compiled with. Note that it is important 
to use the same Boost version for this code. 

# Installation 

This tutorial guides you how to install [NLMech](https://github.com/nonlocalmodels/NLMech) from scratch. In addition, we provide some [Docker files](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/) to build the code on [Fedora](https://getfedora.org/), since we use the same OS on [Circle-CI](https://app.circleci.com/pipelines/github/nonlocalmodels/NLMech). In addition, we provide some [bash script](https://github.com/nonlocalmodels/buildscripts/tree/main/bash) to install the code in one step. 

In this installation guide, we will use a set of [scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure) to build the code on HPC clusters.

## Prerequisites 

### Tools

To compile NLMech and its dependencies the following tools are needed:
* GCC compiler collection (gcc) > 4.9, however, gcc >= 8 is recommended
* [autoconf](https://www.gnu.org/software/autoconf/)
* [wget](https://www.gnu.org/software/wget/)
* [cmake](https://cmake.org/)
* [git](https://git-scm.com/)

### Dependencies

* [Lapack](http://www.netlib.org/lapack/)
* [BLAS](http://www.netlib.org/blas/)
* [Openssl](https://www.openssl.org/)
* [OpenGL](https://www.opengl.org/)

We recommend to install these libraries using the package manager. Bewlow you can find some
examples to install them using apt and dnf

```bash
apt-get install build-essential git wget cmake libssl-dev libblas-dev liblapack-dev autoconf freeglut3-dev
```

```bash
dnf install @development-tools cmake git wget blas-devel lapack-devel freeglut-devel
```

## Using the HPCBuildInfrastructure

First, we clone the [HPCBuildInfrastructure](https://github.com/nonlocalmodels/HPCBuildInfrastructure) 

```bash
git clone https://github.com/nonlocalmodels/HPCBuildInfrastructure.git
cd HPCBuildInfrastructure/
```
The uses version of each library is defined in the `config.sh`, if you need to change them. 

### Build cmake

To build all dependencies, the script `build-all.sh` is used and the last argument is always the dependency to be build.

```bash
./build-all.sh Release without-gcc cmake
```
The first command `Release` specifies the `CMAKE_BUILD_TYPE` of the build. Note that you should use the same build type for all of the dependencies. The second command specifies if you use the gcc from the system in `/usr/bin` by using `without-gcc`. If you do not have access to gcc > 4.9 on your system, we provide the option `with-gcc` to use your own build of gcc. Note you have to run `./build-all.sh Release with-gcc gcc` first to build your gcc.

### Build HPX

```bash
./build-all.sh Release without-gcc hwloc
./build-all.sh Release without-gcc jemalloc
./build-all.sh Release without-gcc boost
./build-all.sh Release without-gcc hpx
```

### Build Blaze and blaze_iterative

```bash
./build-all.sh Release without-gcc blaze
./build-all.sh Release without-gcc blazeIterative
```

### Build VTK

```bash
./build-all.sh Release without-gcc vtk
```

### Build yampl-cpp

```bash
./build-all.sh Release without-gcc yamlcpp
```

### Build NLMech

```bash
./build-all.sh Release without-gcc NLMech
```

All these instructions are available in a single [bash script](https://github.com/nonlocalmodels/buildscripts/blob/main/bash/buildAll.sh) and a [Dockerfile](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/FedoraAll).
# Brief overview of the equations

The governing equations and the discretization is briefly introduced on this page. For more details we
refer to the following two references [1,2].

## Governing equations

Unlike classical continuum mechanics, where the internal force is written in terms of the stress, in PD, the internal force at a given material point is due to the sum of the pairwise forces of the neighboring points. i.e. the force in PD is expressed as the integral of the pairwise force density between the given point and another point in the neighborhood. Neighborhood of point <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;x" title="x" /></a> is typically taken as the ball of radius <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\delta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\delta" title="\delta" /></a>, centered at <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;x" title="x" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\delta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\delta" title="\delta" /></a> is the nonlocal length scale and is referred to as horizon. PD is often divided in two classes: bond-based and state-based models. In bond-based models, the two material points interact via a pairwise force law and the forces between the material points do not depend on the deformation state of surrounding points. In contrast, in the state-based models the volumetric deformation in the neighborhood of two points plays a role in the pariwise force. The governing equation of motion for the bond-based PD [3] reads as

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X})&space;=&space;\int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X})&space;d\mathbf{X}'&space;&plus;&space;\mathbf{b}(t,\mathbf{X})&space;\text{&space;in&space;}&space;D" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X})&space;=&space;\int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X})&space;d\mathbf{X}'&space;&plus;&space;\mathbf{b}(t,\mathbf{X})&space;\text{&space;in&space;}&space;D" title="\varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) = \int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X}) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D" /></a>

and the governing equation for the state-based PD [4] reads as 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\varrho&space;(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X})&space;=&space;\int\limits_{B_\delta(\mathbf{X})}&space;(T[\mathbf{X},t]\langle&space;\mathbf{X}'&space;-&space;\mathbf{X}&space;\rangle&space;-&space;T[\mathbf{X}',t]\langle&space;\mathbf{X}&space;-&space;\mathbf{X}'&space;\rangle)&space;d\mathbf{X}'&space;&plus;&space;\mathbf{b}(t,\mathbf{X})&space;\text{&space;in&space;}&space;D&space;\text{.}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\varrho&space;(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X})&space;=&space;\int\limits_{B_\delta(\mathbf{X})}&space;(T[\mathbf{X},t]\langle&space;\mathbf{X}'&space;-&space;\mathbf{X}&space;\rangle&space;-&space;T[\mathbf{X}',t]\langle&space;\mathbf{X}&space;-&space;\mathbf{X}'&space;\rangle)&space;d\mathbf{X}'&space;&plus;&space;\mathbf{b}(t,\mathbf{X})&space;\text{&space;in&space;}&space;D&space;\text{.}" title="\varrho (\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) = \int\limits_{B_\delta(\mathbf{X})} (T[\mathbf{X},t]\langle \mathbf{X}' - \mathbf{X} \rangle - T[\mathbf{X}',t]\langle \mathbf{X} - \mathbf{X}' \rangle) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D \text{.}" /></a>

Here <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\varrho" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\varrho" title="\varrho" /></a> denotes density of the material, <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\mathbf{u}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\mathbf{u}" title="\mathbf{u}" /></a> displacement field in the material, <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\ddot{\mathbf{u}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\ddot{\mathbf{u}}" title="\ddot{\mathbf{u}}" /></a> acceleration, and <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\mathbf{b}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\mathbf{b}" title="\mathbf{b}" /></a> external force density. The constitutive law, relating bond strain with bond force, is prescribed using either the pairwise force function <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\mathbf{f}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\mathbf{f}" title="\mathbf{f}" /></a> or the PD state <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;T" title="T" /></a> [4].


## Discretization 

Currently, the library supports finite difference (or more generally meshfree) discretization. Using the triangulation of arbitrary domain, the library can create a meshfree discretization. The library is equipped with necessary modules, such as FE elements and quadrature integration rules, for finite element discretization of PD. Next, we briefly discuss the finite difference/meshfree discretization of PD.  The domain <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;D" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;D" title="D" /></a> is discretized with the nodes <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;X&space;=&space;\{&space;X_i&space;\in&space;\mathbb{R}^3&space;\vert&space;i=1,\ldots,n\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;X&space;=&space;\{&space;X_i&space;\in&space;\mathbb{R}^3&space;\vert&space;i=1,\ldots,n\}" title="X = \{ X_i \in \mathbb{R}^3 \vert i=1,\ldots,n\}" /></a>. Each node <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;X_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;X_i" title="X_i" /></a> represents a small area/volume denoted by <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;V_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;V_i" title="V_i" /></a>. In PD, as previously mentioned, each point <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;X_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;X_i" title="X_i" /></a> interacts with neighboring points in ball (discrete) <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;B_\delta(X_i)&space;=&space;\{X_j:&space;\vert&space;X_i&space;-&space;X_j|&space;<&space;\delta&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;B_\delta(X_i)&space;=&space;\{X_j:&space;\vert&space;X_i&space;-&space;X_j|&space;<&space;\delta&space;\}" title="B_\delta(X_i) = \{X_j: \vert X_i - X_j| < \delta \}" /></a>. 

The discrete equation of motion is written as, for the bond-based PD,

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\varrho(X_i)\ddot{\mathbf{u}}(t,X_i)&space;=&space;\sum\limits_{j&space;\in&space;B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i)&space;V_j&space;&plus;&space;\mathbf{b}(t,X_i)&space;\text{&space;in&space;}&space;D," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\varrho(X_i)\ddot{\mathbf{u}}(t,X_i)&space;=&space;\sum\limits_{j&space;\in&space;B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i)&space;V_j&space;&plus;&space;\mathbf{b}(t,X_i)&space;\text{&space;in&space;}&space;D," title="\varrho(X_i)\ddot{\mathbf{u}}(t,X_i) = \sum\limits_{j \in B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i) V_j + \mathbf{b}(t,X_i) \text{ in } D," /></a>

and, the state-based PD,

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\varrho&space;(X_i)\ddot{\mathbf{u}}(t,X_i)&space;=&space;\sum\limits_{j&space;\in&space;B_\delta(X_i)}&space;(T[X_i,t]\langle&space;X_j&space;-&space;X_i&space;\rangle&space;-&space;T[X_j,t]\langle&space;X_i&space;-&space;X_j&space;\rangle)&space;V_j&space;&plus;&space;\mathbf{b}(t,X_i)&space;\text{&space;in&space;}&space;D&space;\text{.}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\varrho&space;(X_i)\ddot{\mathbf{u}}(t,X_i)&space;=&space;\sum\limits_{j&space;\in&space;B_\delta(X_i)}&space;(T[X_i,t]\langle&space;X_j&space;-&space;X_i&space;\rangle&space;-&space;T[X_j,t]\langle&space;X_i&space;-&space;X_j&space;\rangle)&space;V_j&space;&plus;&space;\mathbf{b}(t,X_i)&space;\text{&space;in&space;}&space;D&space;\text{.}" title="\varrho (X_i)\ddot{\mathbf{u}}(t,X_i) = \sum\limits_{j \in B_\delta(X_i)} (T[X_i,t]\langle X_j - X_i \rangle - T[X_j,t]\langle X_i - X_j \rangle) V_j + \mathbf{b}(t,X_i) \text{ in } D \text{.}" /></a>

Here <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\mathbf{u}(t,X_i)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\mathbf{u}(t,X_i)" title="\mathbf{u}(t,X_i)" /></a> denotes the displacement of node <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;X_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;X_i" title="X_i" /></a> at time <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;0&space;\leq&space;t\leq&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;0&space;\leq&space;t\leq&space;T" title="0 \leq t\leq T" /></a>. For the time discretization, we can consider: *1)* implicit time integration and *2)* explicit time integration using either central difference or velocity verlet scheme.


## References

1. P. Diehl, P. K. Jha, H. Kaiser, R. Lipton, and M. Lévesque. An asynchronous and task-based implementation of peridynamics utilizing hpx—the C++ standard library for parallelism and concurrency. SN Applied Sciences, 2(12):2144, 2020, [10.1007/s42452-020-03784-x]({https://doi.org/10.1007/s42452-020-03784-x), [Preprint](https://arxiv.org/abs/1806.06917). 
2. Jha, PK, Lipton R. "Numerical convergence of finite difference approximations for state based peridynamic fracture models."  Computer Methods in Applied Mechanics and Engineering, 1 July 2019, 351(1), 184 - 225. [Link](https://doi.org/10.1016/j.cma.2019.03.024)
3. Silling, S. A., & Askari, E. (2005). A meshfree method based on the peridynamic model of164solid mechanics.Computers & Structures,83(17-18), 1526–1535. [https://doi.org/10.1651016/j.compstruc.2004.11.026](https://doi.org/10.1651016/j.compstruc.2004.11.026)
4. Silling, S. A., Epton, M., Weckner, O., Xu, J., & Askari. (2007). Peridynamic states and167constitutive modeling.Journal of Elasticity,88(2), 151–184.[https://doi.org/10.1007/168s10659-007-9125-1](https://doi.org/10.1007/168s10659-007-9125-1)# Configuration file

We use the [YAML](https://de.wikipedia.org/wiki/YAML) file format for the configuration files. For more details on the YAML file format, we refer to this [tutorial](https://gettaurus.org/docs/YAMLTutorial/). Here, we introduce the option of within the YAMl file to control the simulation.  

## Mesh generation

To generate a mesh the tool [mesh]() is provided and it is used as

```sh
mesh -i input_mesh.yaml -d 1
```

where `-d` specify the dimension of the problem and `-i` specifies the YAML fine with the mesh information. Let us look at the [example](https://github.com/nonlocalmodels/NLMech/tree/main/examples/qsModel/1D) for the one-dimensional implicit time integration example:


```yaml
Output:
  Path: ./
  Mesh: mesh
  File_Format: vtu
Domain:
  - 0.
  - 16
Horizon: 2
Horizon_h_Ratio: 4 
Compress_Type: zlib
```

### Output

The tag `Output` describes details about the generate `.vtu` files using the following attributes:

* `Path:` Describes the path where the generated mesh is written. Using the value `./` means the file will be written in the folder where the executable `mesh` was executed. 
* `Mesh:` Describes the name of the generated file. 
* `File_Format`: Describes the file format of the generated file, e.g. `vtu` and `gmsh`. For the `vtu` file format we refer to the [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) file format documentation and for the `gmsh` file format to the [Gmsh](http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php) documentation.
* `Compress_Type` Describes the compression type of the generated file. 

### Domain

The tag `Domain` describes the simulation domain using the following attributes:

* The dashes `-` list the first point and the last point of the domain
* `Horizon` Describes the length of the horizon as an integer
* `Horizon_h_Ratio` Describes the number of nodes within the horizon

For more complex geometries, the [Gmsh](https://gmsh.info/) tool is recommended. Further examples are available [here](https://nonlocalmodels.github.io/examples/).

## Simulation

Here, all tags and attributes of the YAML files are listed and described. A detailed set of examples is available [here](https://nonlocalmodels.github.io/examples/).

### Model deck

Example of a model deck:

```yaml
Model:
  Dimension: 2
  Discretization_Type:
    Spatial: finite_difference
    Time: velocity_verlet
  Final_Time: 0.001000
  Time_Steps: 50000
  Horizon: 0.002000
  Horizon_h_Ratio: 4
```

The tag `Model` describes the model with following attributes:

* `Dimension` Defines the dimension of the problem, e.g. one-dimensional 
* `Discretization_Type` Describes the discretization in time and space
  * `Time` Defines if a central difference scheme `central_difference` or a velocity verlet `velocity_verlet` scheme is used for the discretization in time
  * `Spatial` Defines if a finite difference `finite_difference` scheme is used
* `Horizon` Defines the horizon
* `Horizon_h_Ratio` Defines the ratio of nodes with the distance of h in the horizon
* `Mesh_Size` Defines the mesh size h
* `Final_Time` Defines the final time 
* `Time_Steps` Defines the amount of time steps

### Restart

Example of a `Restart` deck:

```yaml
Restart:
  File: out/output_3.vtu
  Step: 1500
```

The tag `Restart` describes a restart of the simulation using following attributes:

* `File` The path and file name of the last Successful simulation step
* `Step` The time step to restart from

One example of the restart of one simulation is available [here](https://nonlocalmodels.github.io/examples/restart.html).

### Mesh

Example of a `Mesh` deck:

```yaml
Mesh:
  File: mesh.vtu
  Keep_Element_Conn: true
```

The tag `Mesh` describes the mesh of the simulation using following attributes:

* `File` Path and file name to the mesh either in the gmsh or vtu file format.
* `Keep_Element_Conn` Keep the mesh information available

### Output

Example of a `Output` deck

```yaml
Output:
  Path: out_top_bottom_const/
  Tags:
    - Displacement
    - Velocity
    - Force
    - Damage_Z
  Output_Interval: 500
  Compress_Type: zlib
  Perform_FE_Out: true
```

The tag `Output` describes the mesh of the simulation using following attributes:

* `File_Format` Specifies if the gmsh or vtu file format is used
* `Path` Defines the path were the output is written to
* `Compress_Type` Defines the compression type for the vtu file (`zlib` or `ascii`)
* `Output_Interval` Defines the interval a output file is written
* `Tag` Defines the field appended to the out file. The mesh is always added and all other simulation data needs to be listed.
* `Perform_FE_Out` Store the mesh information in the output

### Boundary conditions

#### Displacement boundary conditions

Example of a `Displacement_BC` deck

```yaml
Displacement_BC:
  Sets: 1
  Set_1:
    Location:
      Rectangle: [0.000000e+00, 9.900000e-02, 1.000000e-01, 1.000000e-01]
    Direction: [1,2]
    Time_Function:
      Type: constant
      Parameters:
        - 0.0
    Spatial_Function:
      Type: constant
      Parameter: [1]
```


The tag `Displacement_BC` describes the displacement boundary conditions and the attribute `Sets` specifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Location` Defines the location of the boundary condition
  *  `Line: [ a, b ]` Applies the boundary condition along the line from `a` to `b`
  *  `Rectangle: [a, b, c, d]` Applies the boundary condition on all nodes inside the rectangle given by the two points `(a,b)` and `c,d)`  
* `Direction: [n]` Defines the direction (x=1, y=2, or z=3) of the boundary condition 
* `Time_Function` Specifies how the boundary condition is applied over time
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load
* `Spatial_Function` Specifies how the boundary condition is applied in space
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load

#### Force boundary conditions

Example of a `Force_BC` deck:

```yaml
Force_BC:
  Sets: 2
  Set_1:
    Location:
      Rectangle: [0, 0, 0.1, 0.002]
    Direction: [2]
    Time_Function:
      Type: linear
      Parameters:
        - -5000000000.0
    Spatial_Function:
      Type: constant
      Parameters: [1]
  Set_2:
    Location:
      Rectangle: [0, 0.098, 0.1, 0.1]
    Direction: [2]
    Time_Function:
      Type: linear
      Parameters:
        - 5000000000.0
    Spatial_Function:
      Type: constant
      Parameters: [1]
```

The tag `Force_BC` describes the displacement boundary conditions and the attribute `Sets` specifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Location` Defines the location of the boundary condition
  *  `Line: [ a, b ]` Applies the boundary condition along the line from `a` to `b`
  *  `Rectangle: [a, b, c, d]` Applies the boundary condition on all nodes inside the rectangle given by the two points `(a,b)` and `c,d)`  
* `Direction: [n]` Defines the direction (x=1, y=2, or z=3) of the boundary condition 
* `Time_Function` Specifies how the boundary condition is applied over time
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load
* `Spatial_Function` Specifies how the boundary condition is applied in space
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load

### Fracture


Example of a `Fracture` deck:

```yaml
Fracture:
  Cracks:
    Sets: 1
    Set_1:
      Orientation: -1
      Line: [5.000500e-02, 0.000000e+00, 5.000500e-02, 2.000000e-02]
```

The tag `Fracture` describes the displacement boundary conditions and the attribute `Sets` specifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Orientation` Describes the orientation of the crack
* `Line` Describes a line and all bonds interesting this line are initially broken

#### References

* C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009. 
* W. Schroeder, Ken Martin, and Bill Lorensen, Visualization Toolkit: An Object-Oriented Approach to 3D Graphics, 4th Edition

