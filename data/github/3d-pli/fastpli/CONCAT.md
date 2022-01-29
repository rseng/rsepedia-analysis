<!--
________             ___________________________
___  __/_____ _________  /___  __ \__  /____  _/
__  /_ _  __ `/_  ___/  __/_  /_/ /_  /  __  /
_  __/ / /_/ /_(__  )/ /_ _  ____/_  /____/ /
/_/    \__,_/ /____/ \__/ /_/     /_____/___/
-->

# Fiber Architecture Simulation Toolbox for 3D-PLI

![fastpli-logo](logo.svg)

The [Fiber Architecture Simulation Toolbox for 3D-PLI (fastpli)](https://github.com/3d-pli/fastpli) is a toolbox for [polarized light imaging (PLI)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) with three main purposes:

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/fiber.png" alt="fiber" width="125">

- [`Sandbox` - designing of nerve fiber models](https://github.com/3d-pli/fastpli/wiki/NerveFiber):
  The first module allows the user to create different types of nerve fiber bundles and additionally fill them with individual nerve fibers.

  - [Details](https://github.com/3d-pli/fastpli/wiki/NerveFiber)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-sandbox)

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/solver_1_cropped.gif" alt="fiber" width="125">

- [`Solver` - generating collision free models](https://github.com/3d-pli/fastpli/wiki/Solver):
  The second module takes as input a configuration of nerve fibers and checks them for spatial collisions.
  Since nerve fibers cannot overlap in reality, one must ensure that the models follow the same rules.
  The solver module implements a simple algorithm that checks for collisions and, if it finds any, pushes the colliding segments of the fibers slightly apart.
  This is repeated until all collisions are solved.

  - [Details](https://github.com/3d-pli/fastpli/wiki/Solver)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-solver)

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/optic_chiasm_10_0.png" alt="fiber" width="125">

- [`Simulation` - simulation of 3D-Polarized Light Imaging](https://github.com/3d-pli/fastpli/wiki/Simulation):
  The simulation module enables the simulation of [3D Polarized Light Imaging (3D-PLI)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html).
  This is a microscopic technique that allows the polarization change of light moving through a brain section to be measured.
  Due to the birefringence property of the myelin surrounding the nerve fibers, the polarization state changes.
  This change enables the calculation of the 3d orientation of the nerve fibers in the brain slice.

  - [Details](https://github.com/3d-pli/fastpli/wiki/Simulation)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-simulation)

## Wiki

<https://github.com/3d-pli/fastpli/wiki>

## Example

As an example, a simplified model of the optic chiasm is presented.
This structure in the brain allows nerve fibers from the eyes to cross each other and connect to the opposite side of the brain.
In addition, a certain portion remains on the same side of the brain.

- [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-optic_chiasm)

![png](https://raw.githubusercontent.com/wiki/3d-pli/fastpli/optic_chiasm_files/optic_chiasm_18_1.png)
![png](https://raw.githubusercontent.com/wiki/3d-pli/fastpli/optic_chiasm_files/optic_chiasm_19_0.png)

## Module lists

[API documentation](https://3d-pli.github.io/fastpli/)

| module                                                                                               | information                                                   |
| ---------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| [`fastpli.analysis` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.analysis.html)           | analysis of 3D-PLI results                                    |
| [`fastpli.io` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.io.html)                       | input/output functions, e.g. to read/save fiber_bundles data  |
| [`fastpli.model.sandbox` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.model.sandbox.html) | building of simple 3d nerve fiber models                      |
| [`fastpli.model.solver` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.model.solver.html)   | generation of non intersection nerve fiber models             |
| [`fastpli.objects` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.objects.html)             | manipulation of fastpli objects (e.g. rotation)               |
| [`fastpli.tools` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.tools.html)                 | mathematical tools and helper function                        |
| [`fastpli.simulation` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.simulation.html)       | simulation of fiber models inside a virtual 3D-PLI microscope |

# Installation

### Note:

> The current version of `fastpli` can only be run under Linux as operating system due to dependencies.
> If you want to use `fastpli` under Windows, please use the Windows subsystem for Linux.
> To enable graphical output, you must install an X server.
> For more information, see <https://wiki.ubuntu.com/WSL>.
> Support for macOS is planned for the future.

## Dependencies

### Requirements

- C++17
- Make
- CMake
- Python3
- MPI
- OpenGL (optional, recommended)

### Submodules

- pybind11

## Install instructions

### Packages

Install all necessary packages.

For Ubuntu:

```sh
sudo apt install gcc g++ cmake make git
sudo apt install python3-dev python3-venv
sudo apt install libopenmpi-dev freeglut3-dev
```

### Clone repository

```sh
git clone --recursive https://github.com/3d-pli/fastpli.git
cd fastpli
```

### Compilation

Use your favorite environment e. g. `python3 -m venv env` and `source env/bin/activate`.
Probably you also have to update your pip version `pip3 install pip -U`.

```sh
make fastpli
pip3 install .
```

# Examples

## Tutorials

```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

jupyter-notebook examples/sandbox.ipynb
jupyter-notebook examples/solver.ipynb
jupyter-notebook examples/simulation.ipynb
jupyter-notebook examples/optic_chiasm.ipynb
```

## Scripts

```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

# run examples
python3 examples/sandbox.py
python3 examples/solver.py
python3 examples/simulation.py
python3 examples/optic_chiasm.py
```

# Tests

```sh
python3 setup.py test
```

# About this Project

## Libraries

All computationally intensive calculations are optimized either with **numba** on the Python side or with multithreading **C++**, which can be accessed via **pybind11**.
Additionally the simulation module supports the **Message Passing Interface (MPI)**.

## Contributions and Bug Reports

Please submit [issues](https://github.com/3d-pli/fastpli/issues) on GitHub to report
problems or suggest features. [Pull requests](https://github.com/3d-pli/fastpli/pulls)
are also welcome to add features or correct problems.
Please run the local env-CI environment `./CI/run-all.sh` or docker container `make docker` in advance.

## Literature

- [3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034)
- [dense fiber modeling](https://arxiv.org/abs/1901.10284)
- [MEDUSA](https://doi.org/10.1016/j.neuroimage.2019.02.055)
- [simulation](https://doi.org/10.1016/j.neuroimage.2015.02.020)
- [tilting analysis](https://doi.org/10.3389/fnana.2018.00075)

## Authors

- **Felix Matuschke**

## References

[fastPLI](https://github.com/3d-pli/fastpli) is an open source toolbox for modeling nerve fibers, simulating them in a [3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034) microscope and the signal processing developed by the [fiber architecture group](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) at the [Forschungszentrum Jülich](https://www.fz-juelich.de) - [INM1](https://www.fz-juelich.de/inm/inm-1/EN/Home/home_node.html).
This project has received funding from the European Union’s Horizon 2020 Research and Innovation Programme under Grant Agreement No. 7202070 ([Human Brain Project](https://www.humanbrainproject.eu/en/) SGA2).

|                                                                                                                                                                                  |                                                                                                                                                              |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
|                  [![Forschungszentrum Jülich](https://www.fz-juelich.de/SharedDocs/Bilder/INM/INM-1/EN/FZj_Logo.jpg?__blob=normal)](https://www.fz-juelich.de)                   | [Forschungszentrum Jülich](https://www.fz-juelich.de)                                                                                                        |
| [![FA-INM-1](https://avatars2.githubusercontent.com/u/51479655?s=200&v=4)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) | [Fiber Architecture - INM1 - Forschungszentrum Jülich](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) |
|                                 [![HBP](https://sos-ch-dk-2.exo.io/public-website-production/img/HBP.png)](https://www.humanbrainproject.eu/en/)                                 | [Human Brain Project](https://www.humanbrainproject.eu/en/)                                                                                                  |

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/3d-pli/fastpli/blob/main/LICENSE) file for details
# TODO

## Issues

* mpi results and non mpi results can differ around eps
* multiprocessing has initial high cpu values with model.solver
* simpli.pixel_size does not call property
* nan values detected in light signal (for 90degree fibers?)
* simpli simulation has high initial copying and allocation time
* doc generation warnings about Classes: stub file not found

## simpli

* [ ] overlap warning only with multiple fiber bundles
* [ ] fiber_bundles have 3 copies, 1 in python, and 2 in cpp (fbs_, fbs_org_)
* [ ] interpolation for mu, dn
* [ ] stokes vs jones
* [ ] polarization value px, py
* [ ] test cells
* [ ] test mpi on multiple nodes

## VCS

* [ ] warning if fiber radius, segment length and r_min are questionable
* [ ] warning if collision close to 0 but "never ending"
* [ ] memory warning for fiber_bundles (gets quite big for small segment lengths)
  * [ ] non linear splitting and merging
* [ ] segment length per fiber, automatic determination for each fiber fun(f_radius)

---

[back](README.md)
---
title: "fastPLI: A Fiber Architecture Simulation Toolbox for 3D-PLI"
tags:
  - 3D-PLI
  - Python
  - microscopy
  - simulation
authors:
  - name: Felix Matuschke
    orcid: 0000-0001-6435-5351
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Katrin Amunts
    orcid: 0000-0001-5828-0867
    affiliation: 1, 2
  - name: Markus Axer
    orcid: 0000-0001-5871-9331
    affiliation: 1
affiliations:
  - name: Institute of Neuroscience and Medicine (INM-1), Forschungszentrum Jülich GmbH, 52425, Jülich, Germany
    index: 1
  - name: Cécile and Oskar Vogt Institute for Brain Research, University Hospital Düsseldorf, University of Düsseldorf, 40204, Düsseldorf, Germany.
    index: 2
date: 18. December 2020
bibliography: paper.bib
---

# Statement of need

3D Polarized Light Imaging (3D-PLI) is a microscopic neuroimaging technique used to study the nerve fiber architecture in unstained histological brain sections at the micrometer scale [@Axer2011].
It provides image contrast for fibers and fiber tracts, and ultimately enables reconstruction of 3D nerve fiber orientations.
The physical effect behind 3D-PLI is the optical property of the nerve fibers called birefringence.
Due to this intrinsic birefringence, it is possible to use polarized light, pass it through a thin brain section and observe the change of the polarization state of light.
This change is directly related to the 3D orientation of the fibers and also provides strong contrasts between fibers and other tissue components.

To understand the influence of the underlying fiber structure, simulations are an essential tool. They allow testing different hypotheses, knowing the ground truth. It has already been shown that simulations with scattered light within tissue sections require models with irregularities to mimic the behavior of the scattered light. This knowledge can now be used to understand structures such as fiber crossings that have been difficult to interpret in 3D-PLI [@Menzel2020].

In addition, the generated nerve fiber models can be used in other imaging simulation techniques such as diffusion magnetic resonance imaging (dMRI).

In recent years, various software tools have been developed to design fibre models.
Nerve fiber modeling is commonly used in dMRI.
But many modeling techniques do not use volumetric representations, especially collision-free ones.
In the last decade, an increasing number of algorithms for non-trivial overlapping structures have been developed [@Altendorf2011; @Chapelle2015; @Mingasson2017; @Ginsburger2019].
While these algorithms are specialized in their field, _fastPLI_ also provides a dedicated tool for 3D-PLI simulation based on linear optics.
Here, the focus is on the previously developed algorithm [@Matuschke2019], which provides a fast method to generate collision-free results for white matter structures in the brain.

Different types of simulations for polarized light are for example described in [RamellaRoman2005; @vanTurnhout2009; @Jiang2020].
However, to our knowledge, none of these techniques have been used to simulate the effects of polarized light on nerve fibers, except `simPLI` [@Dohmen2015], which is included in this toolbox.

# Summary

_fastPLI_ is an open source toolbox based on Python and C++ for modeling myelinated axons, i.e. nerve fibers and simulating the results of measurement of fiber orientations with a polarization microscope using 3D-PLI.

The _fastPLI_ package includes the following modules:

1. **Fiber Modelling Modules:**
   A detailed 3D modelling of nerve fibers at the micrometer level is essential as input for the measurement simulation.
   In order to recreate biological tissue as a model, it is important that the nerve fibers do not spatially overlap.
   We have decided to implement a solver module that takes any configuration of fiber bundles as input and converts it over several iterations into a collision-free configuration.
   In order to generate collision free fiber arrangements, a dedicated algorithm to prohibit such overlaps has been developed [@Matuschke2019].

2. **Simulation Module:**
   The 3D-PLI simulation is based on Stokes vector and Müller matrix approaches as described in [@Dohmen2015; @Menzel2015].
   For the simulation the polarimetric setup can be equipped with a tiltable specimen stage [@Axer2011; @Schmitz2018].
   By this means the brain section can be scanned from oblique views which adds important information to unambiguously analyze the 3D fiber orientation.

3. **Analysis Module:**
   The resulting simulated measurements (i.e., image stacks of a section acquired at different polarizing filter rotation angles and, optionally, at different oblique views) can be processed similarly to the real, experimental 3D-PLI [@Axer2011; @Schmitz2018].

All computationally intensive calculations are optimized either with _numba_ on the Python side or with multithreading _C++_ algorithms, which can be accessed via _pybind11_ inside the Python package [@Lam2015;@pybind11].
Additionally, the simulation module supports the Message Passing Interface (MPI) to facilitate the simulation of very large volumes on multiple computer nodes.

# Installation

The _fastPLI_ package has to be built with a _C++17_ compiler.
A Makefile allows a simple local installation.
It generates the necessary libraries inside a build folder and a matching setup.py file, which can be used for a second installation process with pip.

```sh
make fastpli
pip3 install .
```

All necessary configurations are handled in the background with CMake.
All required software libraries are listed in the [software repository](https://github.com/3d-pli/fastpli).

# Usage & Examples

A more detailed description of _fastPLI_ including examples, jupyter notebooks and tutorials can be found in the software package [^1], at the wiki pages [^2] and the API documentation [^3].

[^1]: [https://github.com/3d-pli/fastpli](https://github.com/3d-pli/fastpli/tree/main/examples)
[^2]: [https://github.com/3d-pli/fastpli/wiki](https://github.com/3d-pli/fastpli/wiki)
[^3]: [https://3d-pli.github.io/fastpli/](https://3d-pli.github.io/fastpli/)

```sh
python3 examples/sandbox.py
python3 examples/solver.py
python3 examples/simulation.py
python3 examples/optic_chiasm.py
```

## Fiber Modelling

Two modules exist to allow the user to build non-colliding white matter nerve fiber models: `fastpli.model.sandbox` and `fastpli.model.solver`.

The `fastpli.model.sandbox` module contains multiple functions to build from coordinates single individual fibers up to fiber bundles composed of individual fibers.

The `fastpli.model.solver` module contains a `Solver` class to generate non-colliding fibers constellations from an initial input dataset [@Matuschke2019].

An example of the solving process is shown in fig. \ref{generation} a).
The red colored segments indicate that a collision with this segment is detected.
With further iterations, the number of collisions decreases.
At the end, a collision-free configurations results.
Figure \ref{generation} b) shows the cross-section through the resulting fiber configuration.
Each fiber bundle has a different gray value.

![a) Solving process of two crossing fiber bundles. Individual colliding segments are colored in red. b) Cross-section of the collision solving process.\label{generation}](generation.png){width=100%}

The resulting orientations for each segment can be then visualized in a polar histogram (see fig. \ref{orientation}).

![orientation distribution a) initial, b) resulting fiber configuration. \label{orientation}](orientation_plot.png){width=100%}

## Simulation and Analysis

The module `fastpli.simulation` contains the class `Simpli`.
This class contains all tools required to simulate the 3D-PLI setup [@Dohmen2015; @Menzel2015] and to analyze the generated images.
Final results include the stack of images as well as the derived modalities referred to as transmittance, direction, retardation, inclination, relative thickness, and 3D fiber orientation (FOM) maps [@Schmitz2018].

An example of the simulation results is shown in fig. \ref{crossingmodalities}.
a) indicates the resulting transmittance map.
It represents the mean light intensity transmission through the cut tissue.
b) represents the fiber direction in the tissue plane.
c) shows the resulting retardation, which is a measure of the amplitude of the resulting sinusoidal signal from the simulation in each pixel.

![PLI modalities: a) transmittance, b) direction, c) retardation\label{crossingmodalities}](crossing_modalities.png){width=100%}

The data can be further analyzed with an advanced tilting simulation and analysis \ref{crossingrofl}.
a) represents the fiber direction in the tissue plane.
b) visualizes the fiber inclination in the tissue plane.
c) shows the relative thickness of the fibers within a pixel volume.
d) represents the resulting FOM. The colors encode the 3D orientation.

![Tilting analysis results: a) direction, b) inclination, c) relative thickness, d) FOM\label{crossingrofl}](crossing_rofl.png){width=100%}

The analysis is accessible inside the simulation module via a _pipeline_ inside the `fastpli.simulation.Simpli` class. For further information see `examples/simulation_pipeline.py`.

# Acknowledgements

This study was funded by the European Union's Horizon 2020 Research and Innovation Programme under Grant Agreement No. 785907 (Human Brain Project SGA2) and No. 945539 (Human Brain Project SGA3).
We gratefully acknowledge the computing time granted through JARA-HPC on the supercomputer JURECA [@jureca] at Forschungszentrum Jülich, Germany.

# Author Contributions

F.M. developed the software and wrote the documentation.
M.A. provided guidance throughout the project.
K.A. and M.A. supervised the project.
F.M. wrote the manuscript with input from all authors.

# References
# CI - Scripts

These scripts provide all CI functions.
They also allow developers to run them localy or the provided dockerfile

**All scripts have to be run from the repository path!**

## Environment

To be able to check all functionalities, a virtual python environment `env-CI` will be created.

```sh
./CI/run-env.sh
```

## Checking Code Format

The `C++` formatation is checked via clang-format.
The style is define in the `.clang-format` file n the repository.

For the python formation `flake8` is used.

```sh
./CI/run-code-check.sh
```

## Checking Build

The `fastpli` package will be installed inside the `env-CI` environment.

```sh
./CI/run-build.sh
```

## Checking Tests

All tests inside the `.\tests\` path will be run.

```sh
./CI/run-pytest.sh
```

## Checking Examples

All examples inside the `./examples` path will be checked for runnability.

```sh
./CI/run-examples.sh
```

## API and WIKI:

- make docs
- make wiki
- switch to ssh credential
- push new version
- notebook md tutorials have to be formated by `pretty`
