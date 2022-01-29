geoclaw-landspill
=================

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/barbagroup/geoclaw-landspill/raw/master/LICENSE)
[![Travis CI](https://img.shields.io/travis/com/barbagroup/geoclaw-landspill/master?label=Travis%20CI)](https://travis-ci.com/barbagroup/geoclaw-landspill)
[![GitHub Action CI](https://img.shields.io/github/workflow/status/barbagroup/geoclaw-landspill/CI/master?label=GitHub%20Action%20CI)](https://github.com/barbagroup/geoclaw-landspill/actions?query=workflow%3ACI)
[![status](https://joss.theoj.org/papers/fb7b012799a70c9b4c55eb4bb0f36f97/status.svg)](https://joss.theoj.org/papers/fb7b012799a70c9b4c55eb4bb0f36f97)
[![Conda](https://anaconda.org/barbagroup/geoclaw-landspill/badges/installer/conda.svg)](https://anaconda.org/barbagroup/geoclaw-landspill)

***Note: if looking for content of `geoclaw-landspill-cases`, please checkout tag
`v0.1`. This repository has been converted to a fully working solver package.***

*geoclaw-landspill* is a package for running oil overland flow simulations for
applications in pipeline risk management. It includes a numerical solver and
some pre-/post-processing utilities.

<center><img src="./doc/sample.gif" /></center>

The numerical solver is a modified version of
[GeoClaw](http://www.clawpack.org/geoclaw.html).
GeoClaw solves full shallow-water equations. We added several new features and
utilities to it and make it usable to simulate the overland flow from pipeline
ruptures. These features include:

* adding point sources to mimic the rupture points
* adding evaporation models
* adding Darcy-Weisbach bottom friction models with land roughness
* adding temperature-dependent viscosity
* recording detail locations and time of oil flowing into in-land waterbodies
* downloading topography and hydrology data automatically (the US only)
* generating CF-1.7 compliant NetCDF files

## Documentation
1. [Dependencies, installation, and tests](doc/deps_install_tests.md)
2. [Usage](doc/usage.md)
3. [Configuration file: `setrun.py`](doc/configuration.md)
4. [Example cases](cases/README.md)
5. [Containers: Docker and Singularity](doc/container.md)

------------------------------------------------------------------------
## Quick start

We only maintain compatibility with Linux. Though using `pip` or building from
source may still work in Mac OS or Windows (e.g., through WSL), we are not able
to help with the installation issues on these two systems.

Beyond this quick start, to see more details, please refer to the
[documentation](#documentation) section.

### 1. Installation

The fast way to install *geoclaw-landspill* is through
[Anaconda](https://www.anaconda.com/)'s `conda` command. The following command
creates a conda environment (called `landspill`) and installs the package and
dependencies:

```
$ conda create \
    -n landspill -c barbagroup -c conda-forge \
    python=3.8 geoclaw-landspill
```

Then use `conda activate landspill` or
`source <conda installation prefix>/bin/activate landspill` to activate the
environment. Type `geoclaw-landspill --help` in the terminal to see if
*geoclaw-landspill* is correctly installed.

### 2. Running an example case

To run an example case under the folder `cases`, users have to clone this
repository. We currently don't maintain another repository for cases. After
cloning this repository, run
```
$ geoclaw-landspill run <path to an example case folder>
```
For example, to run `utal-flat-maya`:
```
$ geoclaw-landspill run ./cases/utah-flat-maya
```
Users can use environment variable `OMP_NUM_THREADS` to control how many CPU
threads the simulation should use for OpenMP parallelization.

### 3. Creating a CF-compliant NetCDF raster file

After a simulation is done, users can convert flow depth in raw simulation data
into a CF-compliant NetCDF raster file. For example,
```
$ geoclaw-landspill createnc ./case/utah-flat-maya
```
Replace `./cases/utah-flat-maya` with the path to another desired case.

QGIS and ArcGIS should be able to read the resulting NetCDF raster file.

------------------------------------------------------------------------
## Third-party codes and licenses

* amrclaw: https://github.com/clawpack/amrclaw
  ([BSD 3-Clause License](https://github.com/clawpack/amrclaw/blob/ee85c1fe178ec319a8403503e779d3f8faf22840/LICENSE))
* geoclaw: https://github.com/clawpack/geoclaw
  ([BSD 3-Clause License](https://github.com/clawpack/geoclaw/blob/3593cb1b418fd52739c186a8845a288037c8f575/LICENSE))
* pyclaw: https://github.com/clawpack/pyclaw
  ([BSD 3-Clause License](https://github.com/clawpack/pyclaw/blob/a85a01a5f20be1a18dde70b7bb37dc1cdcbd0b26/LICENSE))
* clawutil: https://github.com/clawpack/clawutil
  ([BSD 3-Clause License](https://github.com/clawpack/clawutil/blob/116ffb792e889fbf0854d7ac599657039d7b1f3e/LICENSE))
* riemann: https://github.com/clawpack/riemann
  ([BSD 3-Clause License](https://github.com/clawpack/riemann/blob/597824c051d56fa0c8818e00d740867283329b24/LICENSE))

------------------------------------------------------------------------
## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

------------------------------------------------------------------------
## Contact

Pi-Yueh Chuang: pychuang@gwu.edu

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
pychuang@gwu.edu or labarba@gwu.edu.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations

Contributing
============

Thank you for considering contributing to this project. *geoclaw-landspil* is an
open-source project hosted on GitHub at
https://github.com/barbagroup/geoclaw-landspill.
All contributions are welcome, including but not limited to:

* bug reports,
* bug fixes,
* documentation improvement,
* more examples or tests,
* new features, and
* performance enhancement.

Don't hesitate to ask any questions! Questions help us know what can be improved
and make this project more helpful to people.

--------------------
## How to contribute

Please open an issue at our
[GitHub repository](https://github.com/barbagroup/geoclaw-landspill) if there
are any questions, bugs, or suggestions.

Make a pull request if you want to add new features/modifications, fix bugs, or
add more examples/tests. As we try to keep the project simple, all pull requests
should make the `master` branch as the base branch. A pull request against the
`master` branch triggers Travis CI and tests defined with GitHub Actions. Please
make sure the newly added code passes all tests.

------------------------------------------------------------------------
## How to install *geoclaw-landspill* for development and run tests

In addition to the standard installation methods described in
[`deps_install_tests.md`](doc/deps_install_tests.md), developers should consider
an editable installation. With an editable installation, code modifications take
effect immediately without re-installation. Assuming currently under the
top-level folder, do

```
$ pip install --editable .
```

You may need to uninstall previously installed *geoclaw-landspill* with
`$ pip uninstall geoclaw-landspill`.

To run tests, see the last section in
[`deps_install_tests.md`](doc/deps_install_tests.md).
---
title: "geoclaw-landspill: an oil land-spill and overland flow simulator for pipeline rupture events"
tags:
  - shallow-water-equations
  - overland-flow
  - land-spill
  - pipeline
  - geoclaw
authors:
  - name: Pi-Yueh Chuang
    orcid: 0000-0001-6330-2709
    affiliation: 1
  - name: Tracy Thorleifson
    affiliation: 2
  - name: Lorena A. Barba
    orcid: 0000-0001-5812-2711
    affiliation: 1
affiliations:
  - name: Department of Mechanical and Aerospace Engineering, the George Washington University, Washington, DC, USA
    index: 1
  - name: G2 Integrated Solutions, Houston, TX, USA
    index: 2
date: 15 January 2021
bibliography: paper.bib
---

# Summary

The package *geoclaw-landspill* builds on the *geoclaw* shallow-water solver to numerically simulate oil land-spills and overland flows that occur during pipeline accidents.
It helps understanding how oil flows above ground after accidental release and to study a pipeline's impact on the environment and economy.
By understanding the impact of a to-be-constructed pipeline, one can choose a pipeline route with the least loss if accidents were to happen.
On the other hand, understanding how oil flows also helps develop rescue teams' deployment strategies and remedy plans for potential accidents.

The package provides a numerical solver for the full shallow-water equations, and post-processing utilities.
The solver is an expanded version of GeoClaw [@Berger2011],
a parallel shallow-water equation solver for tsunami simulations using adaptive mesh refinement (AMR) and finite-volume methods.
We added several new features and modifications to simulate the overland flow of pipeline rupture events, including (citations refer to the details of the adopted models):

* point sources with multi-stage inflow rates to mimic rupture points along a pipeline;
* Lewis Squires Correlation for temperature-dependent flow viscosity [@mehrotra_generalized_1991];
* Darcy-Weisbach friction model with multi-regime coefficient models (laminar, transient, and turbulent regimes) [@Yen2002] and Churchill's model [@churchill_friction-factor_1977];
* inland waterbody interactions;
* Fingas' evaporation models [@fingas_modeling_2004]; and
* optimizations to improve performance in overland flow simulations.

In addition to the numerical solver, *geoclaw-landspill* is also able to:

* automatically download high-resolution topography and hydrology data, and
* create CF-compliant NetCDF raster files for mainstream GIS software (e.g., QGIS, ArcGIS).

*geoclaw-landspill* includes several examples to showcase its capability to simulate the overland flow on different terrain.

We implemented the core numerical solver and new features using Fortran 2008 to better integrate into the original GeoClaw code.
Other utilities are in Python.
Users can find *geoclaw-landspill* as a Python package on PyPI and install it through pip.
The only dependency that pip does not manage is the Fortran compiler.
Docker images and Singularity images are also available.
They ease the deployment of the solver and simulations to cloud-based high-performance clusters.

# Statement of need

In the US, between 2010 and 2017, an average of 388 hazardous liquid pipeline accidents happened per year.
Half of accidents contaminate soil, and 41% of accidents affect areas with high consequences in either ecology or economy.
Moreover, 85%, on average, of the released oil was not recovered and kept damaging the environment
[@belvederesi_statistical_2018].
From the perspective of risk management, while pipelines are unavoidable in modern days, it is necessary to understand how a pipeline may impact the environment if any accidental release happens.
*geoclaw-landspill* serves this purpose.
It provides a free and open-source simulation tool to researchers investigating the danger, risk, and loss posed by potential pipeline accidents.

To our knowledge, *geoclaw-landspill* is the only open-source high-fidelity flow simulator for oil pipeline rupture events.
High fidelity means the results provide more details and accuracy because of high-resolution digital elevation data, fine spatial discretization, and full shallow-water equations.
Commercial products with a similar capability to *geoclaw-landspill* are available [@Zuczek2008; @RPSGroup; @Hydronia; @Gin2012].
Other non-commercial software more or less serving a similar purpose usually relies on simplified models, such as 1D open-channel models, diffusive wave approximation, gravity current models, and gradient-based route selection models [@Hussein2002; @Simmons2003; @Ronnie2004; @farrar_gis_2005; @Guo2006; @Su2017].
Moreover, these non-commercial codes are no longer available or are not open-source.

Another value of *geoclaw-landspill* is that it provides a platform for scholars who study oil flow modeling to implement and test their models.
As the main flow solver is under the BSD 3-Clause License, scholars can add their models to *geoclaw-landspill* freely.

# Past or ongoing research projects using the software

The following conference presentations and posters used previous versions of
*geoclaw-landspill*:

1. Chuang, P.-Y., Thorleifson, T., & Barba, L. A. (2019a, May). GeoClaw-ArcGIS Integration for Advanced Modeling of Overland Hydrocarbon Plumes. *2019 Petroleum GIS Conference Proceedings. 2019 Esri Petroleum GIS Conference, Houston, TX, USA*.

2. Chuang, P.-Y., Thorleifson, T., & Barba, L. A. (2019b, July). Python Workflow for High-Fidelity Modeling of Overland Hydrocarbon Flows with GeoClaw and Cloud Computing. *Proceedings of the 18th Python in Science Conference. 18th Python in Science Conference (SciPy 2019), Austin, TX, USA*.


# Acknowledgments

The project received funding from G2 Integrated Solutions, Houston, TX, USA.

# Author Contributions

Pi-Yueh Chuang developed and is the maintainer of this software.
Tracy Thorleifson was involved in software design, model selection, and software testing.
Lorena A. Barba designed the software specification, development roadmap, and software framework.
She also oversaw the progress of the project.

# Reference
Examples
========

Currently, all examples are from the nearby region of Salt Lake City, Utah, USA.
The naming of examples follows
**{location}-{terrain type}-{oil type}-{additional setting}**.

There are three files Under each example folder:

* `setrun.py`: it serves as the configuration file of the simulation.
* `roughness.txt`: all examples use a Darcy-Weisbach model that considers
  surface roughness. This file provides roughness information to the model,
  though the values are simply `0.1` everywhere.
* `animation.gif`: this is an animation of the flow to give users a sense of how
  their results should look.

For the details of `setrun.py`, please refer to the documentation.

To run an example, do

```
$ OMP_NUM_THREADS={number of CPU threads to use} geoclaw-landspill run {path to a case folder}
```

To create the `animation.gif`s shown below, do

```
$ geoclaw-landspill plotdepth --use-sat {path to a case folder}
```

This command only creates frames used by the animations. The frames are saved
in the folder `_plots/sat/level02/` under each case. Finally, we use
[`ffmpeg`](https://ffmpeg.org/) to combine frames into a GIF animation file:

```
$ ffmpeg \
    -i {case folder}/_plots/sat/level02/frame00%03d.png \
    -vf "fps=20,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" \
    {filename to save the resulting animation}
```

For a more serious analysis, users should consider creating a temporal NetCDF
raster file for each case with

```
$ geoclaw-landspill createnc {path to a case folder}
```

The resulting NetCDF file will be at `_output/{case name}-depth-lvl02.nc`.

----------
## List of cases and results

### 1. utah-flat-gasoline-no-evap

Gasoline without evaporation on flat terrain.

![utah-flat-gasoline-no-evap](./utah-flat-gasoline-no-evap/animation.gif)

### 2. utah-flat-gasoline

Gasoline with evaporation on flat terrain. File `_output/evaporated_fluid.dat`
stores the total volume evaporated at the end of the simulation.

![utah-flat-gasoline](./utah-flat-gasoline/animation.gif)

### 3. utah-flat-maya-no-evap

Maya crude oil without evaporation on flat terrain.

![utah-flat-maya-no-evap](./utah-flat-maya-no-evap/animation.gif)

### 4. utah-flat-maya

Maya crude oil with evaporation on flat terrain. File `_output/evaporated_fluid.dat`
stores the total volume evaporated at the end of the simulation.

![utah-flat-maya](./utah-flat-maya/animation.gif)

### 5. utah-hill-maya-no-hydro

Maya crude oil on hilly terrain. Evaporation is on, but the in-land waterbody
(a creek in this case) is turned off to see how the solver works with a
drainage system.

![utah-hill-maya-no-hydro](./utah-hill-maya-no-hydro/animation.gif)

### 6. utah-hill-maya

Maya crude oil on hilly terrain. Evaporation is on. The in-land waterbody
(a creek in this case) is on. The creek catches all oil.
`_output/removed_fluid.csv` records contact points and the fluid volumes flowing
into waterbodies at the contact points.


![utah-hill-maya](./utah-hill-maya/animation.gif)

### 7. utah-structures-maya

Maya crude oil on flat terrain but with some structures in the area. Both
evaporation and in-land waterbodies are on.

![utah-structures-maya](./utah-structures-maya/animation.gif)

### 8. utah-waterbody-gasoline-no-evap

Gasoline above land surrounded by in-land waterbodies.
`_output/removed_fluid.csv` records contact points and the fluid volumes flowing
into waterbodies at the contact points. The evaporation is off.

![utah-waterbody-gasoline-no-evap](./utah-waterbody-gasoline-no-evap/animation.gif)

### 9. utah-waterbody-gasoline

Gasoline above land surrounded by in-land waterbodies.
`_output/removed_fluid.csv` records contact points and the fluid volumes flowing
into waterbodies at the contact points. The evaporation is on.

![utah-waterbody-gasoline](./utah-waterbody-gasoline/animation.gif)

### 10. utah-waterbody-maya-no-evap

Maya crude oil above land surrounded by in-land waterbodies.
`_output/removed_fluid.csv` records contact points and the fluid volumes flowing
into waterbodies at the contact points. The evaporation is off.

![utah-waterbody-maya-no-evap](./utah-waterbody-maya-no-evap/animation.gif)

### 11. utah-waterbody-maya

Maya crude oil above land surrounded by in-land waterbodies.
`_output/removed_fluid.csv` records contact points and the fluid volumes flowing
into waterbodies at the contact points. The evaporation is on.

![utah-waterbody-maya](./utah-waterbody-maya/animation.gif)
# Usage

`geoclaw-landspill` is the top-level executable. It has several subcommands. To
see the subcommands:
```
$ geoclaw-landspill --help
```

-----------------
## Run a case

To run a case, execute

```
$ geoclaw-landspill run <path to a case folder>
```

Interested users can find example cases in directory `cases`.

A case folder must have at least a `setrun.py` that configures the simulation.
The usage of `setrun.py` follows the convention of GeoClaw and Clawpack-related
projects.

The command `run` automatically downloads topography (in the US) and
hydrology data from the USGS database if it can not find specified files.
For example, running the `utah-flat-maya` case under `cases` downloads
topography and hydrology files, `utah-flat.asc` and `utah-flat-hydro.asc`,
into the folder `common-files`.

The solver runs with OpenMP-parallelization. By default, the number of threads
involved in a run is system-dependent. Users can explicitly control how many
threads to use for a simulation by environment variable `OMP_NUM_THREADS`. See
OpenMP's documentation. For example, to run `utah-flat-maya` with 4 threads:

```
$ OMP_NUM_THREADS=4 geoclaw-landspill run <path to utah-flat-maya>
```

Raw simulation results are under folder `<case folder>/_output`. If running a
case multiple times, old `_output` folders are renamed automatically to
`_output.<timestamp>` to avoid losing old results.

Use the following commands to generate GIS-software-readable raster files or to
produce quick visualizations.

------------------------------------------------
## Create NetCDF raster files with CF convention

```
$ geoclaw-landspill createnc <path to a case>
```

This subcommand converts raw simulation results to a CF-compliant temporal NetCDF
raster file. At least QGIS and ArcGIS are able to read this file format. The
resulting NetCDF file will be at `<case folder>/_output/<case name>-level<XX>.nc`.
Only the results at one AMR level are used. The default level is the finest AMR
level. Use flag `--level=<number>` to create raster files for other AMR levels.

To see more parameters to control the conversion, use
```
$ geoclaw-landspill createnc --help
```

--------------------------------------
## Create quick Visualization of depth

```
$ geoclaw-landspill plotdepth <path to a case>
```

This subcommand creates flow depth contours for all time frames. Figures are
saved to `<case folder>/_plots/depth/level<xx>`. By default, only the results on
the finest AMR grid level are plotted.

There are several parameters to fine-tune the output images, including plotting
the depth data on a satellite image (such as the one shown in README). See
`--help` for more arguments.

-----------------------------------------------------
## Visualize elevation data used by runtime AMR grids

```
$ geoclaw-landspill plottopo <path to a case>
```

`plottopo` plots the runtime topography on AMR grids at each time frame of a 
simulation. This is different from plotting the topography using a topography
file. Runtime topography means the actual elevation values defined on AMR grid
cells during a simulation.

The output figures are saved to `<case folder>/_plots/topo`.

See `--help` for more arguments.

## Calculate total volume above ground

```
$ geoclaw-landspill volumes <path to a case>
```

This subcommand calculates total fluid volumes at each AMR grid level at each
time frame. The results are saved to a CSV file `<case folder>/_output/volumes.csv`.
The main use case of this volume data is to check the mass conservation.
# Containers: Docker and Singularity

We provide Docker images on Docker Hub and Singularity images on Singularity Hub.

---------------
## Docker usage

Pull the Docker image through:
```
$ docker pull barbagroup/landspill:<version>
```

Replace `<version>` with a desired version string. For example, `v1.0.dev1`. To
see all tags/versions, check the
[Docker Hub page](https://hub.docker.com/repository/docker/barbagroup/landspill).

To get into the shell of a Docker container:
```
$ docker run -it --name landspill barbagroup/landspill:<version>
```
Once getting into the shell, the example cases are under
`~/geoclaw-landspill-cases`. *geoclaw-landspill* is ready to use. For example,
in the shell of the Docker container, run the `utah-flat-maya` case with
```
$ geoclaw-landspill run ~/geoclaw-landspill-cases/utah-flat-maya
```

Interested users can build Docker image locally using the Dockerfile in folder
`Dockerfiles`. For example, to build an image with version `1.0.dev1` (the
corresponding tag is `v1.0.dev1`):
```
$ cd Dockerfiles
$ docker build \
    --tag "my_image" --build-arg "VER=v1.0.dev1" --target production
    --file Dockerfile .
```
`my_image` will be the image name in the local Docker registry.

------------------------------------------------------------------------
## Singularity usage (Singularity version >= 3.4)

On many HPC clusters or supercomputers, Docker is not available due to
security concerns, and [Singularity](https://www.sylabs.io/singularity/) is the 
only container technology available. 

To pull a Singularity image with tag `<version>` (e.g., `v1.0.dev1`) and save to
a local image file, do
```
$ singularity pull lanspill-<version>.sif shub://barbagroup/geoclaw-landspill:<version>
```

One can get into the shell of a Singularity container and run cases as described
in the section of Docker usage. Also, the Singularity image can be used as an
executable, which is equivalent to the executable `geoclaw-landspill`. For
example, to run a case:
```
$ ./landspill-<version>.sif run <case on the host>
```

Or, to see the help of subcommand `plotdepth`:
```
$ ./landspill-<version>.sif plotdepth --help
```

Alternatively, to run a Singularity image with the regular approach and, for
example, to create a NetCDF file:
```
$ singularity run landspill-<version>.sif createnc <path to a case>
```

Singularity automatically maps some folders on the host into a container. For
example, a user's `$HOME`. So if there's a case at a user's `$HOME`, such as
`$HOME/smaple_case`, then the user can do
```
$ ./landspill-<version>.sif run $HOME/sample_case
```
without copying case data into the container. Please refer to Singularity's
documentation for more details.

Interested users can build Singularity images locally with the configuration
files in folder `Singularityfiles`. Currently, there is one Singularity file
per release version:
```
# singularity build <the name of the image> Singularity.<version>
```
`<the name of the image>` is whatever name a user would like to use for his/her
image.
# Dependencies, installation, and tests

The only operating system officially supported is Linux. We are not maintaining
compatibility with other systems, though they may still work.

---------------
## Dependencies

Build time and runtime dependencies are described in `requirements-build.txt`
and `requirements.txt`, respectively. `gfortran >= 8.0` is the only build time
dependency not included in the file, and, correspondingly, `libgfortran5 >= 8.0`
is the only runtime dependency not included in the file.

Anaconda users do not need to worry about dependencies at all.

On the other hand, `pip` users have to install `gfortran` or/and `libgfortran5`
in advance using the package managers of their Linux distributions. For example,
in Arch Linux, use:
```
# pacman -S gcc-fortran
```
And in Ubuntu 20.04:
```
# apt install gfortran
```

Alternatively, though not recommended, one can use `conda` to get `gfortran` and
then continue using `pip` for other dependencies. The command to get `gfortran`
from Anaconda is
```
$ conda install -c conda-forge "gfortran_linux-64>=8.0"
```
`conda` renames the compiler executable to `x86_64-conda-linux-gnu-gfortran`.
However, this should not concern users because CMake should be able to find it
automatically.

After installing `gfortran` manually, `pip` users can continue on the
installation of *geoclaw-landspill* (the next section).

---------------
## Installation

### Option 1: use `conda` and install binary files from Anaconda

As described in README, the following command creates an environment called
`landspill`, and it has *geoclaw-landspill* installed:
```
$ conda create \
    -n landspill -c barbagroup -c conda-forge \
    python=3.8 geoclaw-landspill
```
Once activate the environment, the executable `geoclaw-landspill` should already
be available.

### Option 2: use `pip` to install from PyPI or from source

Note, when using the `pip` command, users can always add the `--user` flag to
install to users' local paths and avoid root privilege. However, if using the
`--user` flag, users should make sure `pip`'s local `bin` path is in `PATH`.

#### Option 2.1: install from PyPI

To install the package from PyPI. 
```
$ pip install geoclaw-landspill
```

We only distribute source tarballs on PyPI due to the requirement of a Fortran
compiler. Wheels or binary releases of this package are not available. `pip`
will download the source tarball, compile/build the package, and then install
it. `gfortran` has to be installed in advance as described in the previous
section.

#### Option 2.2: install with a source tarball from GitHub

Download a release tarball from the repository's
[release page](https://github.com/barbagroup/geoclaw-landspill/releases) on GitHub,
and install the package directly with pip and the tarball:
```
$ pip install <tarball name>.tar.gz
```

#### Option 2.3: install with the repository in developer mode

Clone/pull the repository from GitHub:
```
$ git clone --recurse-submodules https://github.com/barbagroup/geoclaw-landspill.git
```

Go into the folder, and then install the Python dependencies with:
```
$ pip install -r requirements.txt -r requirements-build.txt
```
Then, install *geoclaw-landspill*:
```
$ pip install --editable .
```

Under the developer mode, installation is just a link referencing the source
directory, so any changes in the source Python files take effect immediately (
but not the Fortran files because they have to be re-compiled).

--------
## Tests

To run tests without installation, users can use
[`tox`](https://tox.readthedocs.io/en/latest/) to do so:

```
$ tox
```

End-users can run tests against the installed package if they install
*geoclaw-landspill* through `pip` and using the source tarball or code
repository. Use `pytest`:
```
$ pytest -v tests
```

Currently, the number and the coverage of the tests are limited. It's still a
WIP.
Configuration file: `setrun.py`
===============================

A simulation case folder must at least contains a `setrun.py`, which describes
the configurations of a simulation. `setrun.py` defines a function called
`setrun` that accepts no argument and returns an instance of
`gclandspill.data.ClawRunData`. For those familiar with Clawpack's ecosystem,
`gclandspill.data.ClawRunData` is a modified version of
`clawpack.clawutil.data.ClawRunData`. Our modified `ClawRunData` has some
default settings configured to meet *geoclaw-landspill*'s needs. Also, it has
additional settings for landspill simulations.

The basic working function `setrun` looks like:

```python
import gclandspill

def setrun():
    """Whatever docstring."""

    rundata = gclandspill.data.ClawRunData()

    <...setting up rundata...>

    return rundata
```

For those unfamiliar with Clawpack/GeoClaw, roughly speaking, there are four
parts in the `rundata` object that needs to be configured: core solver, AMR
(adaptive mesh refinement), GeoClaw-specific, and landspill-specific
configurations. The core solver, AMR, and GeoClaw configurations are described
in the official Clawpack documentation:

* Core solver: [Specifying classic run-time parameters in setrun.py](http://www.clawpack.org/setrun.html)
* AMR: [Specifying AMRClaw run-time parameters in setrun.py](http://www.clawpack.org/setrun_amrclaw.html)
* GeoClaw: [Specifying GeoClaw parameters in setrun.py](http://www.clawpack.org/setrun_geoclaw.html)

This document covers only required settings and also settings specific for
landspill. To better understand how a `setrun.py` looks like, please see the
`setrun.py` in examples under the `cases` folder.

--------------------------------------------
## I. Settings in core solver, AMR, and GeoClaw

Most settings in the core solver, AMR, and GeoClaw are tuned to work with
*geoclaw-landspill*, so end-users usually don't need to configure them. Here we
list some settings in these three parts that users must configure:

* `rundata.clawdata.lower[:]`: A list of two `float` representing *xmin* and
  *ymin*, respectively. In other words, lower boundaries in x and y.

* `rundata.clawdata.upper[:]`: A list of two `float` representing *xmax* and
  *ymax*, respectively. In other words, upper boundaries in x and y.

* `rundata.clawdata.num_cells[:]`: A list of two `int` representing the numbers
  of cells at the coarsest AMR level in the x and y direction.

* `rundata.topo_data.topofiles`: A list of two-element lists. For each of the
  two-element lists, the first element is an `int` indicating the topography
  file's format. The topography file is specified in the second element. If a
  relative path is used for the topography file, it must be relative to the
  case's folder. Usually, the format is `3`, which means it's an
  [Esri ASCII format raster](https://en.wikipedia.org/wiki/Esri_grid) file.
  See [Topography data](https://www.clawpack.org/topo.html#topo) for other
  acceptable raster formats.

* `rundata.clawdata.output_style`: An `int` specifying which output style to
  use for saving temporal results. Each style comes with different follow-up
  settings that need to be set:

  * `1`: Output a fixed number of frames at equally spaced times up to a final
    time.
    * `rundata.clawdata.num_output_times`: an `int`; total number of frames
    * `rundata.clawdata.tfinal`: a `float`; time of the last frame
    * `rundata.clawdata.output_t0`: a `bool`; whether to additionally output the
      initial state
  * `2`: Specify a list of times.
    * `rundata.clawdata.output_times`: a list of `float`; output at these times
  * `3`: Specify a number of time steps, and output every this number of time
    steps up to a given number of total steps.
    * `rundata.clawdata.output_step_interval`: an `int`; output every this
      number of steps
    * `rundata.clawdata.total_steps`: an `int`; output until the total number of
      steps reaches this value
    * `rundata.clawdata.output_t0`: a `bool` indicating whether to also output
      the initial state

--------------------------------------
## II. Settings specific for landspill

### i. Basic

* `rundata.landspill_data.ref_mu`: Reference dynamic viscosity at
  `ref_temperature`. The default is the viscosity of Maya crude oil. Unit: cP
  (i.e., mPa-s, or 1e-3 kg/s-m). (default: 332.0)

* `rundata.landspill_data.ref_temperature`: The temperature at which the
  `ref_mu` is measured. Unit: Celsius. (default: 15.0)

* `rundata.landspill_data.ambient_temperature`: The ambient temperature in
  simulation. The solver adjusts the viscosity and density based on this
  temperature. Unit: Celsius. (default: 25.0)

* `rundata.landspill_data.density`: The density measured at `ref_temperature`.
  The default is the density of Maya crude oil. Unit: kg/m^3. (default: 926.6)

### ii. Point sources

A point source is used to mimic a rupture point along a pipeline. It provides
fluid inflow into a computational domain. The inflow profile can be a
multi-stage constant rate.

* `rundata.landspill_data.point_source.n_point_sources`: An `int` indicating how
  many point sources. ***NOTE: currently only cases with one point source are
  well tested.*** (default: 0)

* `rundata.landspill_data.point_source.point_sources`: A list of four-element
  lists. Each four-element list describes a point source. (default: [])
  * The first element is a length-2 list of x and y coordinates of the point
    source.
  * The second element is the number of stages in the multi-stage inflow rate.
  * The third element is a list indicating the end time of each stage.
  * The fourth element is a list providing the volumetric rate of each stage.

  For example, if there is one single point source located at x=10, y=11. and if
  it has three stages of inflow profile (excluding when the rate is zero):
  * 1 m^3/sec when time &lt; 60 seconds
  * 0.5 m^3/sec when time &ge; 60 and &lt; 1800 seconds
  * 0.1 m^3/sec when time &ge; 1800 and &lt; 7200 seconds
  * 0 m^3/sec after 7200 seconds

  Then the corresponding configuration is

  ```python
  rundata.landspill_data.point_source.n_point_sources = 1
  rundata.landspill_data.point_source.point_sources = [
     [[10., 11.], 3, [60., 1800., 7200.], [1., 0.5, 0.1]]
  ]
  ```

### iii. Darcy-Weisbach friction model

In Darcy-Weisbach friction model, the friction coefficient can be calculated
by several different models. We implemented some of them. The type of models to
use is controlled by:

* `rundata.landspill_data.darcy_weisbach_friction.type`: (default: 0)

Each type has its own set of follow-up parameters.

* type is `0`: No Darcy-Weisbach friction.

* type is `1`: Use a constant friction coefficient everywhere. The only
follow-up parameter is
`rundata.landspill_data.darcy_weisbach_friction.coefficient` (default: 0.25).

* type is `2`: Block-wise constant coefficients. Users specify several blocks in
  a computational domain, and each block has a constant coefficient.  Available
  parameters are:

    * `rundata.landspill_data.darcy_weisbach_friction.n_blocks`: Number of blocks.
      (default: 0) 
    * `rundata.landspill_data.darcy_weisbach_friction.xlowers`: The lower boundaries
      in x direction of all blocks. (default: [])
    * `rundata.landspill_data.darcy_weisbach_friction.ylowers`: The lower boundaries
      in y direction of all blocks. (default: [])
    * `rundata.landspill_data.darcy_weisbach_friction.xuppers`: The upper boundaries
      in x direction of all blocks. (default: [])
    * `rundata.landspill_data.darcy_weisbach_friction.yuppers`: The upper boundaries
      in y direction of all blocks. (default: [])
    * `rundata.landspill_data.darcy_weisbach_friction.coefficients`: The
      coefficients in all blocks. (default: [])
    * `rundata.landspill_data.darcy_weisbach_friction.default_coefficient`: Regions
      that are not covered by blocks will have this value. (default: 0.25)

* type is `3`: Cell-wise coefficients. The coefficients are set through a raster
  file in Esri ASCII format.

    * `rundata.landspill_data.darcy_weisbach_friction.filename`: the name of the
      raster file for friction coefficients. (default: "")
    * `rundata.landspill_data.darcy_weisbach_friction.default_coefficient`: the
      coefficient for regions not covered by the raster file. (default: 0.25)

* type is `4`: Three-regime coefficient model. The coefficient is determined
in each cell based on whether the cell is laminar, transient, or turbulent
flow. See reference [1]. 

    * `rundata.landspill_data.darcy_weisbach_friction.filename`: the name of the
      raster file for surface roughness. (default: "")
    * `rundata.landspill_data.darcy_weisbach_friction.default_roughness`: the
      surface roughness for regions not covered by the raster file. (default:
      0.0)

* type is `5`: Churchill's coefficient model. See reference [2].

    * `rundata.landspill_data.darcy_weisbach_friction.filename`: the name of the
      raster file for surface roughness. (default: "")
    * `rundata.landspill_data.darcy_weisbach_friction.default_roughness`: the
      surface roughness for regions not covered by the raster file. (default:
      0.0)

* type is `6`: Two-regime coefficient. Similar to the three-regime model but
ignore transient flow regime.

    * `rundata.landspill_data.darcy_weisbach_friction.filename`: the name of the
      raster file for surface roughness. (default: "")
    * `rundata.landspill_data.darcy_weisbach_friction.default_roughness`: the
      surface roughness for regions not covered by the raster file. (default:
      0.0)

Usually, type `4` or `5` is used because they consider the surface roughness.

### iv. In-land waterbodies

*geoclaw-landspill* is able to records how much fluid volume flowing into
in-land waterbodies and at what locations. Users have to provide rasterized
hydrology data to enable this feature.

* `rundata.landspill_data.hydro_features.files`: A list of filenames. If
  relative paths are provided, they are assumed to be relative to the case
  folder. (default: [])

***Note: using multiple files has not been widely tested. It's recommended to
combine all data into one single raster file.***

Usually, hydrology data are vector data, such as those obtained from the USGS
database. Users can burn the vector data into raster files using
[`gdal_rasterize`](https://gdal.org/programs/gdal_rasterize.html).
Alternatively, users can specify a filename without actually providing the file.
*geoclaw-landspill* automatically downloads and rasterize hydrology data in this
case.

### v. Evaporation model

*geoclaw-landspill* implements Fingas' 1996 model, including the natural log
model and square root model. See reference [3].

* `rundata.landspill_data.evaporation.type`: The type of evaporation models to
  use. (default: 0)

  * `0`: No evaporation.
  * `1`: Fingas' natural log model
  * `2`: Fingas' square root model

* `rundata.landspill_data.evaporation.coefficients`: a list of `float` for model
  coefficients. Currently, we only have Fingas' models, so this variable should
  be a list of two `float`. See reference [3] for coefficients of a variety of
  oils.

-------------------------------------------------------------------
## III. Optional settings affecting stability, accuracy, and performance

The following optional settings are commonly used to control numerical
stability, accuracy, or/and performance:

* `rundata.clawdata.dt_initial`: This determines the size of the first time
  step. By default, it is automatically determined by the following
  formula:

  ```
  dt_initial = 0.3 * (cell size at the finest AMR grid) / (inflow rate)
  ```

  This makes the cell containing the first point source to have a depth of 0.3
  meters at the end of the first time step. Currently, the coefficient, 0.3, is
  hard-coded. This may give some trouble to simulations under small scales.
  For example, if the inflow rate of a point source is 1e-6 m^3 / sec, and if
  the finest cell size is 1e-4 m^2, this formula returns a `dt_initial` of 30
  seconds, which is apparently too big for a problem at this scale. In this
  case, users may want to manually set `dt_initial`.

* `rundata.clawdata.dt_max`: The maximum time step size allowed during a
  simulation. By default, the solver adjusts time step sizes based on flow
  conditions. Users sometimes may want to limit or increase how big a step size
  can be. (default: 5.0)

* `rundata.clawdata.cfl_desired`: The desired Courant–Friedrichs–Lewy (CFL)
  number. The solver adjusts time step sizes to make the resulting CFL numbers
  below this value. (default: 0.9)

* `rundata.clawdata.cfl_max`: The maximum allowed CFL number. If a CFL number is
  larger than this value, it triggers the solver to adjust time step sizes.
  (default: 0.95)

* `rundata.refinement_data.variable_dt_refinement_ratios`: Whether to use
  GeoClaw's adaptive time-step refinement mechanism. When AMR refines a coarse
  cell into several smaller cells, the smaller cells need smaller time step
  sizes. The time step sizes' refinement ratios can be fixed values or
  automatically determined by the GeoClaw solver. (default: True)

* `rundata.amrdata.amr_levels_max`: The number of AMR levels to be used. By
  default, *geoclaw-landspill* uses two levels: the coarse level for dry
  regions and the fine level for wet regions. (default: 2)

* `rundata.amrdata.refinement_ratios_x`: The refinement ratio between each two
  consecutive AMR levels in the x direction. (default: [4])

* `rundata.amrdata.refinement_ratios_y`: The refinement ratio between each two
  consecutive AMR levels in the y direction. (default: [4])

* `rundata.amrdata.refinement_ratios_t`: The refinement ratio between each two
  consecutive AMR levels in time step sizes. This only have effect when
  `rundata.refinement_data.variable_dt_refinement_ratios` is `False`. (default:
  [4])

* `rundata.geo_data.dry_tolerance`: The solver considers cells with depth below
  this value to be dry cells. (default: 1.e-4)

* `rundata.landspill_data.update_tol`: Affect how a coarse cell should
  distribute the depth to its children cells. This setting affects the
  conservation of mass. (default: the same as the `dry_tolerance`)

* `rundata.landspill_data.refine_tol`: By default, as long as there's fluid in a
  cell, it is refined by AMR (even if the depth is smaller than `dry_tolerance`.
  The behavior can be changed with this variable. (default: 0.0)

------------
## Reference

[1] B. C. Yen, “Open Channel Flow Resistance,” Journal of Hydraulic Engineering,
vol. 128, no. 1, pp. 20–39, Jan. 2002, doi: 10.1061/(ASCE)0733-9429(2002)128:1(20).

[2] S. W. Churchill, “Friction-factor equation spans all fluid flow regimes,”
Chem. Eng., vol. 84, no. 24, pp. 91–92, 1977

[3] M. F. Fingas, “Modeling evaporation using models that are not boundary-layer
regulated,” Journal of Hazardous Materials, vol. 107, no. 1–2, pp. 27–36, Feb.
2004, doi: 10.1016/j.jhazmat.2003.11.007.
