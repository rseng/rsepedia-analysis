Material Point Method (mpm) {#mainpage}
===========================

An explicit Material Point Method code.
# CB-Geo Material Point Method contributors

> ⚠ This document reflects the contributions of team members as of August 2020.

## MPM Framework
* Krishna Kumar, UT Austin, USA. [krishnak@utexas.edu](krishnak@utexas.edu)
* Shyamini Kularathna, UC Berkeley, USA.
* Ezra Tjung, UC Berkeley, USA. [ezrayst@gmail.com](ezrayst@gmail.com)
* Christopher Wilkes, University of Cambridge, UK.
* Tianchi Zhao, UC Berkeley, USA. [ztc@tongji.edu.cn](ztc@tongji.edu.cn)
* Thiago Araujo, UT Austin, USA. [thiago.araujo@utexas.edu](thiago.araujo@utexas.edu)
* Bodhinanda Chandra, UC Berkeley, USA.
* Weijian Liang, UC Berkeley, USA and HKUST, HK. [wliangab@ust.hk](wliangab@ust.hk)
* Jeffrey Salmond, University of Cambridge, UK.
* Wentao Xu, University of Pennsylvania, USA. [wentaoxu@seas.upenn.edu](wentaoxu@seas.upenn.edu)
* Xinyi Qian, UC Berkeley, USA.
* Joel Given, UC Berkeley, USA. [joelgiven@berkeley.edu](joelgiven@berkeley.edu)
* Yong Liang, UC Berkeley, USA. [yliang_sn@berkeley.edu](yliang_sn@berkeley.edu)

## Generalized Interpolation Material Point Method
* Christopher Wilkes

## Domain decomposition / Parallelization schemes
* Krishna Kumar
* Jeffrey Salmond
* Wentao Xu

## Two-phase MPM formulation
* Shyamini Kularathna
* Tianchi Zhao
* Weijian Liang
* Bodhinanda Chandra

## Constitutive models
* Tianchi Zhao (Mohr-Coulomb, Modified Cam Clay)
* Ezra Tjung (Bingham, NorSand)
* Shyamini Kularathna
* Krishna Kumar (Linear Elastic, Newtonian)
* Bodhinanda Chandra
* Christopher Wilkes

## Visualization
* Krishna Kumar

## Professors
* Kenichi Soga, UC Berkeley, USA.
* Giovanna Biscontin, University of Cambridge, UK.


# Developer Certificate of Origin
> Version 1.1

By contributing to CB-Geo, You accept and agree to the following 
terms and conditions for Your present and future Contributions 
submitted to CB-Geo Except for the license granted herein to 
CB-Geo and recipients of software distributed by CB-Geo, You 
reserve all right, title, and interest in and to Your Contributions.

Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
1 Letterman Drive
Suite D4700
San Francisco, CA, 94129

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.

## Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.

# High-Performance Material Point Method (CB-Geo mpm)
> [CB-Geo Computational Geomechanics Research Group](https://www.cb-geo.com)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/mpm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/mpm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://mpm.cb-geo.com/)
[![CircleCI](https://circleci.com/gh/cb-geo/mpm.svg?style=svg)](https://circleci.com/gh/cb-geo/mpm)
[![codecov](https://codecov.io/gh/cb-geo/mpm/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-geo/mpm)
[![](https://img.shields.io/github/issues-raw/cb-geo/mpm.svg)](https://github.com/cb-geo/mpm/issues)
[![Coverity](https://scan.coverity.com/projects/14389/badge.svg)](https://scan.coverity.com/projects/14389/badge.svg)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/cb-geo/mpm.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cb-geo/mpm/context:cpp)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/orgs/cb-geo/projects/1)
[![Discourse forum](https://img.shields.io/badge/forum-mpm-blueviolet.svg)](https://cb-geo.discourse.group/c/mpm/)

## Documentation

Please refer to [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for information on compiling, and running the code. The documentation also include the MPM theory.

If you have any issues running or compiling the MPM code please open a issue on the [CB-Geo Discourse forum](https://cb-geo.discourse.group/c/mpm/).

## Running code on Docker

* Docker image for CB-Geo mpm code [https://hub.docker.com/r/cbgeo/mpm](https://hub.docker.com/r/cbgeo/mpm)

* Instructions for running mpm docker container: [https://github.com/cb-geo/docker-mpm/blob/master/README.md](https://github.com/cb-geo/mpm-container/blob/master/README.md).

## Running code locally

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Boost](http://www.boost.org/)
* [Eigen](http://eigen.tuxfamily.org/)
* [HDF5](https://support.hdfgroup.org/HDF5/)

#### Optional
* [MKL](https://software.intel.com/en-us/mkl)
* [MPI](https://www.open-mpi.org/)
* [OpenMP 5.0](https://www.openmp.org/specifications/)
* [KaHIP](https://github.com/schulzchristian/KaHIP)
* [Partio](https://github.com/wdas/partio)
* [VTK](https://www.vtk.org/)

### Fedora installation (recommended)

Please run the following command:

```shell
dnf install -y boost boost-devel clang clang-analyzer clang-tools-extra cmake cppcheck dnf-plugins-core \
                   eigen3-devel findutils freeglut freeglut-devel gcc gcc-c++ git hdf5 hdf5-devel \
                   kernel-devel lcov libnsl make ninja-build openmpi openmpi-devel tar \
                   valgrind vim vtk vtk-devel wget
```

### Ubuntu installation

Please run the following commands to install dependencies:

```
sudo apt update
sudo apt upgrade
sudo apt install -y gcc git libboost-all-dev libeigen3-dev libhdf5-serial-dev libopenmpi-dev libomp-dev
```

If you are running Ubuntu 18.04 or below, you may want to update the GCC version to 9 to have OpenMP 5 specifications
support.

```
sudo apt install software-properties-common
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install gcc-9 g++-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90 --slave /usr/bin/g++ g++ /usr/bin/g++-9 --slave /usr/bin/gcov gcov /usr/bin/gcov-9

```

To install other dependencies:
> CMake 3.15
```
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt update
sudo apt upgrade
```

> OpenGL and X11:Xt
```
sudo apt-get install freeglut3-dev libxt-dev
```

> VTK
```
git clone https://gitlab.kitware.com/vtk/vtk.git VTK
cd VTK && mkdir build && cd build/
cmake -DCMAKE_BUILD_TYPE:STRING=Release ..
make -j
sudo make install
```

### Partio for Houdini SFX Visualization

```shell
mkdir -p ~/workspace && cd ~/workspace/ && git clone https://github.com/wdas/partio.git && \
    cd partio && cmake . && make
```

Houdini supported (*.bgeo) files will be generated. These can be rendered using the non-commercial [Houdini Apprentice](https://www.sidefx.com/download/).

### KaHIP installation for domain decomposition

```shell
cd ~/workspace/ && git clone https://github.com/schulzchristian/KaHIP && \
   cd KaHIP && sh ./compile_withcmake.sh
```

## Compile
> See [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for more detailed instructions.

0. Run `mkdir build && cd build && cmake -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `make clean && make -jN` (where N is the number of cores).

> To compile without KaHIP partitioning use `cmake -DNO_KAHIP=True ..`

### Compile mpm or mpmtest

* To compile either `mpm` or `mpmtest` alone, run `make mpm -jN` or `make mpmtest -jN` (where N is the number of cores).

### Compile without tests [Editing CMake options]

To compile without tests run: `mkdir build && cd build && cmake -DMPM_BUILD_TESTING=Off  -DCMAKE_CXX_COMPILER=g++ ..`.

## Compile with MPI (Running on a cluster)

The CB-Geo mpm code can be compiled with `MPI` to distribute the workload across compute nodes in a cluster.

Additional steps to load `OpenMPI` on Fedora:

```
source /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/usr/share/modulefiles
module load mpi/openmpi-x86_64
```

Compile with OpenMPI (with halo exchange):

```
mkdir build && cd build
export CXX_COMPILER=mpicxx
cmake -DCMAKE_BUILD_TYPE=Release -DKAHIP_ROOT=~/workspace/KaHIP/ -DHALO_EXCHANGE=On ..
make -jN
```

To enable halo exchange set `-DHALO_EXCHANGE=On` in `CMake`. Halo exchange is a better MPI communication protocol, however, use this only for larger number of MPI tasks (> 4).

### Compile with Ninja build system [Alternative to Make]

0. Run `mkdir build && cd build && cmake -GNinja -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `ninja`

### Compile with Partio viz support

Please include `-DPARTIO_ROOT=/path/to/partio/` in the cmake command. A typical cmake command would look like `cmake -DCMAKE_BUILD_TYPE=Release -DPARTIO_ROOT=~/workspace/partio/ ..`

## Run tests

0. Run `./mpmtest -s` (for a verbose output) or `ctest -VV`.

## Run MPM
> See [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for more detailed instructions.

The CB-Geo MPM code uses a `JSON` file for input configuration. To run the mpm code:

```
   ./mpm  [-p <parallel>] [-i <input_file>] -f <working_dir> [--]
          [--version] [-h]
```

For example:

```
export OMP_SCHEDULE="static,4"
./mpm -f /path/to/input-dir/ -i mpm-usf-3d.json
```

Where:

```
   -p <parallel>,  --parallel <parallel>
     Number of parallel threads

   -i <input_file>,  --input_file <input_file>
     Input JSON file [mpm.json]

   -f <working_dir>,  --working_dir <working_dir>
     (required)  Current working folder

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

### Running the code with MPI

To run the CB-Geo mpm code on a cluster with MPI:

```
mpirun -N <#-MPI-tasks> ./mpm -f /path/to/input-dir/ -i mpm.json
```

For example to run the code on 4 compute nodes (MPI tasks):

```
mpirun -N 4 ./mpm -f ~/benchmarks/3d/uniaxial-stress -i mpm.json
```

## Authors

Please refer to the [list of contributors to the CB-Geo MPM code](AUTHORS.md).

## Citation

If you publish results using our code, please acknowledge our work by quoting the following paper:

Kumar, K., Salmond, J., Kularathna, S., Wilkes, C., Tjung, E., Biscontin, G., & Soga, K. (2019). Scalable and modular material point method for large scale simulations. 2nd International Conference on the Material Point Method. Cambridge, UK. [https://arxiv.org/abs/1909.13380](https://arxiv.org/abs/1909.13380)
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at info@cb-geo.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
## Styleguides

### Developer Certificate of Origin

Please read the [Developer Certificate of Origin](./dco.md) before contributing.
The project is licensed under MIT.

We appreciate your contribution and suggest following our [developer guidelines](https://mpm.cb-geo.com/#/code/developer)# References

## Volume calculations
[Hexahedron volume](1985_davies.pdf)

## Inverse isoparametric mapping
[Newton Raphson for Hexahedron: Thermomechanical Modeling of Additive Manufacturing Large Parts (Delinger et al 2014](https://manufacturingscience.asmedigitalcollection.asme.org/article.aspx?articleID=1910535)

[Dealii: Affine transformation for initial guess](https://github.com/dealii/dealii/blob/6d75a550b12999a4167b372b51f4affaa80133bb/source/grid/tria_accessor.cc#L1625-L1748)

[A consistent point‐searching algorithm for solution interpolation in unstructured meshes consisting of 4‐node bilinear quadrilateral elements (Zhao et al., 1999)](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0207%2819990810%2945%3A10%3C1509%3A%3AAID-NME643%3E3.0.CO%3B2-1)
**Describe the PR**
A clear and concise description of the PR.

**Related Issues/PRs**
A clear and concise description of what you expected to happen.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project

---

**Describe the feature**
A clear and concise description of what you want and the reasoning behind this request.

**Describe alternatives**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Request for comments
about: Request for comments discussion

---
# (RFC title goes here)

## Summary

> One paragraph explanation of the change.

## Motivation

> Why are we doing this? What use cases does it support? What is the expected
outcome?

## Design Detail

> This is the bulk of the RFC. Include code snippets, class outlines, inheritance schemes,
feature design outlines, etc. Describe how a new feature / function / class will be implemented
in detail. Elaborate on the feature on how it will be implemented with examples.
RFCs should be used as a design and discussion board to have full details of the implementation,
before a feature is implemented. 

> Explain the design in enough detail for somebody
familiar with the infrastructure to understand. This should get into specifics and corner-cases,
and include examples of how the service is used. Any new terminology should be
defined here.

## Drawbacks

> Why should we *not* do this? Please consider the impact on users,
on the integration of this change with other existing and planned features etc.

> There are tradeoffs to choosing any path, please attempt to identify them here.

## Rationale and Alternatives

> Why is this design the best in the space of possible designs?

> What other designs have been considered and what is the rationale for not choosing them?

> What is the impact of not doing this?


## Prior Art

Discuss prior art, both the good and the bad, in relation to this proposal. A few examples of what this can include are:

> How other services / infrastructures in the same domain have solved this problem.

> Are there any published papers or great posts that discuss this? If you have some relevant papers to refer to, this can serve as a more detailed theoretical background.

This section is intended to encourage you as an author to think about the lessons from other organisations, provide readers of your RFC with a fuller picture. If there is no prior art, that is fine - your ideas are interesting whether they are brand new or if it is an adaptation from other services.

## Unresolved questions

> What parts of the design do you expect to resolve through the RFC process before this gets merged?

> What parts of the design do you expect to resolve through the implementation of this feature before stabilisation?

> What related issues do you consider out of scope for this RFC that could be addressed in the future independently of the solution that comes out of this RFC?

## Changelog

> Add edit summaries to the current RFC here.---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of the bug.

**To Reproduce**
Steps to reproduce the behavior:
1. Compile '...'
2. Run on '....'
3. On condition '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Runtime environment (please complete the following information):**
 - OS/Docker image:
 - Branch [e.g. develop]

**Additional context**
Add any other context about the problem here.
