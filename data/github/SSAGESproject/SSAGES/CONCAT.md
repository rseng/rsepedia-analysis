SSAGES v0.9 Release Notes
=========================

## v0.9.2
- Critical bug fix for patching Eigen with GROMACS
- Correction to gradient of AngleCV
- Additional checks and docs for multiwalker loggers
- Support newer LAMMPS versions
- Minor source cleanup

## v0.9.1
- Critical bug fix for neural-network based methods
- Updated Eigen dependency to v3.3.7
- Many documentation improvements
- Support for newer HOOMD versions (2.6+)
- Build process fixes

## v0.9.0
- New Combined Force–Frequency sampling method
- Addition of non-weighted internal center of mass calculation
- Documentation additions
- New ANN-based collective variable

SSAGES v0.8 Release Notes
=========================

## v0.8.6
- New RMSD Collective Variable
- Improved Qbox examples
- Standardized test suite
- Improved documentation
- Addition of Travis-CI testing
- General code cleanup and bugfixes (See commits)

## v0.8.5
- Improved Elastic Band Sampling
- Temporary removal of GROMACS 2019
- Additional Documentation
- General code cleanup and bugfixes (See commits)

## v0.8.4
- Major documentation updates!
- Improved HOOMD-blue support
- ANN restarts
- Python script for Metadynamics example
- Support for newer GROMACS versions (up to 2019.1)
- Several bugfixes (See commits)

## v0.8.3
- HOOMD-blue support!
- Support for newer GROMACS and LAMMPS versions
- CV definition checking for Methods
- Documentation updates
- Several bugfixes (See commits)

## v0.8.2
- Grid internal updates
- ABF Integrator handles interpolation in each direction independently
- ABF restarts now handled with JSON member
- GyrationTensorCV can be projected into any number of dimensions
- Documentation updates
- Several bugfixes (See commits)

## v0.8.1
- GROMACS support for all 5.1.x, 2016.x, 2018.x!
- Better CMake handling for Hooks
- Handling of LAMMPS line continuations (&)
- Correct handling of multiprocessor ABF method
- BFS method cleanup
- Minor documentation updates
- Eigen include update (3.3.4)
- googletest include update (1.8.0)
- jsoncpp include update

## v0.8.0
- Added ANN sampling!
- More documentation updates
- Updates to examples
- Added Fourier and Chebyshev basis sets to BFS
- Added 3D ABF integrator
- Improved periodicity handling on grid
- Secondary structure CV bug fixes
- Improved Qbox integration


SSAGES v0.7 Release Notes
=========================

## v0.7.5
- Major documentation update!
- Major examples update!
- New generalized pairwise CV
- New secondary structure CVs (alpha and anti/parallel beta sheet RMSD)
- New CV logging capability
- CV selection by name in methods
- Updated unit tests
- Updated ABF integrators for periodic and non-periodic CVs
- Gromacs 2016.3 support!
- Many bug fixes!

## v0.7.0
- New simplified JSON syntax
- Support for multiple simultaneous methods!
- Eliminated boost dependency!
- CV selector for methods
- Argument forwarding for Gromacs
- Updated forward flux examples
- Significant under-the-hood improvements
- Fixed Gromacs auto-download


SSAGES v0.6 Release Notes
=========================

## v0.6.0
- Support for QBox first-principles MD engine
- Support for OpenMD engine
- Coordination number CV
- Polymer Rouse modes CV
- Box volume CV
- Virial contribution (NPT support) for some CVs and methods
- Updated examples and documentation
- New backend grid
- Grid-based metadynamics
- Updated forward flux sampling
- Fixed regression with string methods
- Performance and other improvements!


SSAGES v0.5 Release Notes
=========================

## v0.5.0
- Gromacs restart support
- New gyration tensor CVs
- Updated examples and documentation
- Metadynamics optimizations
- Better engine error handling
- More! (See commit log)
<div align="center">
  <a href="http://ssagesproject.github.io" target="_blank">
    <img src="doc/assets/ssages-logo.png" alt="SSAGES" height="200">
  </a>
</div>

<h2 align="center">
<p align="center">
  <a href="http://ssagesproject.github.io/docs/index.html" target="_blank">
    <img src="https://img.shields.io/badge/docs-v0.9-blue.svg" alt="Documentation">
  </a>
  &nbsp;
  <a href="https://doi.org/10.1063/1.5008853" target="_blank">
    <img src="https://img.shields.io/badge/doi-10.1063%2F1.5008853-blue" alt="Cite SSAGES">
  </a>
</p>
</h2>

**SSAGES** (**S**oftware **S**uite for **A**dvanced **G**eneral **E**nsemble
**S**imulations) is an open-source, engine agnostic, C++11 based advanced
sampling package.  It is designed to be easy to use, extendable and extremely
versatile. It is currently pre-beta, meaning that there are many rough edges,
but we are working rapidly to expand its features and fix any bugs. Keep an eye
on this page for future updates and see below on how to contribute!

## What's New (v0.9.3)
- Numerous improvements to the GROMACS hook
- New Qbox examples
- Support for newer LAMMPS versions
- Engine version is now logged
- ANN and CFF schemas now include a CVs section
- Source cleanup: Fetching external dependencies at build time
- Other minor source cleanup changes

To view the full changelog history, refer to [HISTORY](HISTORY.md).

<a id="features"></a>
## Features
**SSAGES** currently works with multiple molecular dynamics engines. It contains a variety of collective variables (CVs) and advanced sampling methods.

### Highlights
- Engine agnostic framework
- Simple JSON input file syntax
- Easy to add new CVs
- Easy to add new methods
- Much more!

### Engines
- GROMACS 5.1.x, 2016.x, 2018.x
- LAMMPS (Most recent versions)
- OpenMD (2.5+)
- QBox (1.63+)

### CVs
- Artificial Neural Network (as a function of group positions)
- Atom group coordinate
- Atom group position
- Atom group separation
- Bend angle
- Box volume
- Components of gyration tensor
- Pairwise kernel (coordination number, nearest neighbors)
- Polymer Rouse modes
- Root-mean-square deviation (RMSD)
- Secondary structure (alpha, anti/parallel beta sheet) RMSD
- Torsional angle

### Methods
- Adaptive biasing force
- Artificial neural network sampling
- Basis function sampling
- Combined Force–Frequency Sampling
- Metadynamics
- Umbrella sampling
- Finite temperature string
- Nudged elastic band
- Swarm of trajectories
- Forward flux sampling

<a id="installation"></a>
## Installation
The first step is to clone the repository locally.

```bash
$ git clone https://github.com/SSAGESproject/SSAGES.git
```
**SSAGES** uses a CMake build system. It also requires the use of a support MD engine.
For example, to compile with LAMMPS, execute the following

```bash
$ cd SSAGES
$ mkdir build && cd build
$ cmake -DLAMMPS_SRC=/path/to/lammps/src ..
$ make
```

This will build a SSAGES executable which will reside in the build directory.

If you want to use a specific compiler (or if your default compiler is not supported),
set the C and C++ compilers with `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`, respectively.
For example, to use gcc/g++, replace the CMake command with

```bash
$ cmake -DLAMMPS_SRC=/path/to/lammps/src -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
```

If you want to compile and run unit and integration tests, replace the cmake command
in the example above with

```bash
$ cmake -DLAMMPS_SRC=/path/to/lammps/src -DBUILD_TESTS=ON ..
```

### MPI

A requisite underlying MPI library also required to run SSAGES.
On recent Debian based systems using OpenMPI, the requirement
can be installed via:

```bash
$ sudo apt-get install libopenmpi-dev openmpi-bin
```

For more detail on the build system, please check the documentation.

To build the documentation, refer to [Documentation README](doc/README.md).

## Known issues
**SSAGES** is currently in pre-beta. That means there may be known issues that are not yet resolved. Major issues are listed here.

- Restarts are not fully functioning for all methods.

## Contributing
Feel free to fork this project on GitHub. Any pull-requests, feature requests or other form of contributions are welcome.
# HOOMD-blue Umbrella Sampling Example

*Author: Bradley Dice (bdice@umich.edu)*

Please refer to the SSAGES documentation on Umbrella Sampling for instructions on how to use this example code.
This folder contains sample outputs from the SSAGES HOOMD-blue Umbrella Sampling example code.

`cv_vs_time.png` plots the collective variable (dihedral angle) over time. This helps check that enough autocorrelation times have passed.

![CV vs Time](cv_vs_time.png "CV vs Time")

`histogram_trajectories.png` shows a histogram from each of the trajectories and the regions of the CV that were sampled.

![Histogram Trajectories](histogram_trajectories.png "Histogram Trajectories")

`histogram_combined.png` shows a histogram summed over all trajectories to ensure that the entire range of angles were sampled.

![Histograms Combined](histogram_combined.png "Histograms Combined")

`wham_free_energy.png` is the free energy as a function of the dihedral angle.

![WHAM Free Energy](wham_free_energy.png "WHAM Free Energy")
This example demonstrates the CFF method by calculating the potential of mean force
alanine dipeptide in water, using Gromacs.

Inside the directory `1walker`, run CFF using a single walker.

To run this example, execute:

```
{SSAGES_path}/build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f npt.mdp -p topol.top -c npt.gro -o adp.tpr
mpirun -np ${nprocessors} {SSAGES_path}/build/ssages CFF.json
```

This run will produce the main output file called `CFF.out`, similar to the file
`CFF.out1000` located in this folder.

Please see the SSAGES documentation under `Methods/CFF` for more information on definition
of input variables.

**Additional notes**

The Gromacs input files can be modified to change any underlying properties of the
simulation. After changing any of these parameters, the `.tpr` file must be regenerated
using `gmx grompp`, as follows:

```
{SSAGES_path}/build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f npt.mdp -p topol.top -c npt.gro -o adp.tpr
```

Free energy surface or the potential of mean force can be plotted from one of the CFF
output data files, `CFF.out`. First two columns in `CFF.out` are the peptide torsional
angles, phi and psi (radians), and the last column is the potential of mean force (kJ/mol)
#RUN PARALLEL TEMPERING
1) lammps input file
replace temperature values with keyword ptemp in Butane_SSAGES.in

2) temeplate json file
set the parameters in the Template_Input.json
MDSteps = total timesteps to run
frequency = timestep frequency for swaping configurations
CVS = it is a dummpy CVs and does have any influence on simulation

3) generate json file
set the parameters in the ParallelTemp_Input_Generator.py
nwalker = number of walkers
min_temp = minimum temperature
max_temp = maximum temperature

run
python ParallelTemp_Input_Generator.py

5) run SSAGES
mpirun -np nwalker ../../../build/ssages ParallelTemp.json

This example demonstrates the free energy surface (FES) of the Cl-C-C-Cl
torsional angle in 1,2-dichloroethane (CH2Cl-CH2Cl) in vacuum, using Qbox as
the engine and ABF as the method.

This run will produce three output files: F_out, Fworld_cv0, Nworld.
The FES can be created by running the script from the Tools folder:

    python /path/to/ssages/Tools/ABF_integrator.py -i F_out -o G -p True

**Additional notes**

In the example input.json, `"md_iterations"` is set to be 80000, which can
create a rough estimate of the FES with respect to the torsional angle.
To get a more accurate estimate of the FES, increase `"md_iterations"` to
improve sampling.



This example calculates the MFEP of the isomerization of alanine dipeptide, using LAMMPS.  To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in basis.csv) by executing:
```
python plotter.py 22 bfs
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using LAMMPS.  This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
To run the example, execute:

```
mpirun -np 44 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using Gromacs.
This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 44 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follow:

```
gmx_mpi grompp -f nvt.mdp -c adp_water.gro -p topol_water.top -o adp_H2O.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf

```
This example calculates the MFEP of a particle moving on a 2D energy surface with two wells (at [-1, -1] and [1, 1]) and a barrier at [0,0].
By default, this example runs with 16 images.  

```
mpirun -np 16 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.
This example calculates the MFEP of the isomerization of alanine dipeptide, using LAMMPS.  This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
To run the example, execute:

```
mpirun -np 44 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in basis.csv) by executing:
```
python plotter.py 22 bfs
```
This example calculates the MFEP of the isomerization of alanine dipeptide, using Gromacs.  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follows:

```
gmx_mpi grompp -f nvt.mdp -c adp.gro -p topol.top -o adp.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide, using Gromacs.  
This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 44 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follows:

```
gmx_mpi grompp -f nvt.mdp -c adp.gro -p topol.top -o adp.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using Gromacs.  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follow:

```
gmx_mpi grompp -f nvt.mdp -c adp_water.gro -p topol_water.top -o adp_H2O.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using LAMMPS.  
To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide, using LAMMPS.  To run the example, execute:

```
mpirun -np 22 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in basis.csv) by executing:
```
python plotter.py 22 bfs
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using LAMMPS.  This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
To run the example, execute:

```
mpirun -np 44 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using Gromacs.
This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 44 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follow:

```
gmx_mpi grompp -f nvt.mdp -c adp_water.gro -p topol_water.top -o adp_H2O.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of a particle moving on a 2D energy surface with two wells (at [-1, -1] and [1, 1]) and a barrier at [0,0].
By default, this example runs with 16 images.  

```
mpirun -np 16 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.
This example calculates the MFEP of the isomerization of alanine dipeptide, using LAMMPS.  This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
To run the example, execute:

```
mpirun -np 44 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in basis.csv) by executing:
```
python plotter.py 22 bfs
```
This example calculates the MFEP of the isomerization of alanine dipeptide, using Gromacs.  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 22 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follows:

```
gmx_mpi grompp -f nvt.mdp -c adp.gro -p topol.top -o adp.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide, using Gromacs.  
This example is configured to run on two processor per walker.
Some parameters of the inputs may differ from the single processor per walker example.
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 44 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follows:

```
gmx_mpi grompp -f nvt.mdp -c adp.gro -p topol.top -o adp.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using Gromacs.  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 22 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follow:

```
gmx_mpi grompp -f nvt.mdp -c adp_water.gro -p topol_water.top -o adp_H2O.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
This example calculates the MFEP of the isomerization of alanine dipeptide in water, using LAMMPS.  
To run the example, execute:

```
mpirun -np 22 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf

```
SSAGES Documentation
====================

SSAGES's documentation can be found in the [online
manual](https://ssagesproject.github.io/docs/index.html). Alternatively, after you
have built SSAGES you can further build the documentation from the `build/`
directory.

## Requirements

Make sure you have the following additional tools installed:

- [Doxygen]              – generate documentation from annotated C++ source code.
- [Graphviz]             – `dot` tool used to draw graph visualizations.
- [Sphinx]               – documentation builder.
- [sphinx-rtd-theme]     – Read the Docs Sphinx theme.
- [sphinxcontrib-bibtex] – Sphinx extensions for BibTeX style citations.

[Doxygen]:              https://www.doxygen.nl/index.html
[Graphviz]:             https://graphviz.org
[Sphinx]:               https://sphinx-doc.org
[sphinx-rtd-theme]:     https://sphinx-rtd-theme.readthedocs.io
[sphinxcontrib-bibtex]: https://sphinxcontrib-bibtex.readthedocs.io

On Debian-based distributions (e.g. Ubuntu), you can easily install them with
the following commands:
```
$ sudo apt install doxygen graphviz python3-sphinx python3-pip
$ pip install sphinx_rtd_theme sphinxcontrib-bibtex
```

## Building

You can build the documentation with
```
$ make doc
```
You can also build the API-references and the User Manual separately with
```
$ make apiref
```
and
```
$ make manual
```

If you have `pdflatex` installed, you can also build a PDF file for the
documentation. To compile the API-reference into a PDF file do
```
$ cd doc/API-doc/latex/
$ make
```
The PDF will be called `refman.pdf`

Similarly, you can build a PDF version of the Manual with
```
$ cd doc/Manual/
$ make
```
The PDF will be called `SSAGES.pdf`

## Viewing the documentation

Once you have built the documentation you will find it in the `doc/API-doc/`
and `doc/Manual` directories. To view the documentation in a browser just open
the files `doc/Manual/index.html` or `doc/API-doc/html/index.html`.
.. _metadynamics:

Metadynamics
------------

Introduction
^^^^^^^^^^^^

Metadynamics defines a class of flat-histogram methods useful for
molecular dynamics simulations. Within metadynamics, a
history-dependent bias is accrued through the periodic application of
elemental Gaussian biases to the collective variables (CVs) of
interest. The form of the individual biases is

.. math::

	g(\vec{\xi},\vec{s}_i) = W_i
	e^{-\frac{\left[\vec{\xi}-\vec{s}_i\right]^2}{2\sigma_i^2}}\;.

Here, :math:`\vec{\xi}` is the collective variable, :math:`\vec{s}_i` is
the location of the :math:`i^\text{th}` hill, :math:`W_i` is the weight and
:math:`\sigma_i` the width of the hill. This bias acts to push the
system away from previously visited states. As this bias accrues to
the height of nearby features in the free energy surface, the system's
trajectory will begin to explore an expanded area in CV space. After a
sufficient number of biases have been applied, the total applied bias
within a given region will begin to oscillate around the negative of
the free energy in that region.

The free energy surface (FES) at any time :math:`t` may be
reconstructed by summing over all applied biases thusly

.. math::

	F(\vec{\xi},t) = -V(\vec{\xi}) = -\sum_{i < n(t)} W_i
	e^{-\frac{\left[\vec{\xi}-\vec{s}_i\right]^2}{2\sigma_i^2}}\;.

Here, :math:`n(t)` refers to the number of biases applied before time
:math:`t`. The time-dependent FES can be used to determine whether or
not a simulation has reached a converged FES by computing block
averages of the applied bias (see :cite:`SINGH2012369`).

Metadynamics can be applied to unbounded regions of CV space, or to
bounded regions. For the case of bounded, non-periodic CVs, boundary
corrections must be applied which alter the structure of the hills
:math:`g(\vec{\xi},\vec{s}_i)` (see :cite:`MCGOVERN2013084102`)

Metadynamics exists in many flavors, in which :math:`W_i` and
:math:`\sigma_i` are altered in a time- or trajectory- dependent
fashion. These include Well-Tempered Metadynamics, Transition-Tempered
Metadynamics and Metadynamics with Adaptive Gaussians. Each of these
methods has advantages and drawbacks. Currently, SSAGES includes only
standard Metadynamics with fixed-shape hills and no tempering
algorithm. Details on how to use these algorithms within SSAGES are
given in the following sections.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Metadynamics is selected by defining ``type" : "Metadynamics"`` as the
method in the JSON input file. It supports the following options:

widths
	*array of doubles (length: number of CVs)*.
	This array defines the width of the Gaussians being deposited over time
	in each CV dimension.

height
	*double*.
	This value defines the height of the Gaussians being deposited over time,
	in units of energy used by the MD engine.

hill_frequency
	*double*
	This value defines the frequency in iterations with which the Gaussians
	are deposited.

.. note::

	The Metadynamics method does not yet support CV bounds.

.. note::

	The Metadynamics method does not yet support restarts.

.. warning::

	Metadynamics will run for the duration specified in the input files. It is
	the user's responsibility to ensure that enough time is given for
	Metadynamics to obtain an satisfactory representation of the free energy
	surface of interest. A simple way to prevent Metadynamics from terminating
	before convergence is to define a very large number of timesteps and to
	check the output file periodically for convergence.

Example Input
^^^^^^^^^^^^^

.. code-block:: javascript

	"method": {
		"type": "Metadynamics",
		"widths": [0.1, 0.1],
		"height": 0.1,
		"hill_frequency": 1
	}

Output
^^^^^^

The output of Metadynamics is stored in a file called "hills.out". This file
contains the location, width and height of each Gaussian that has been deposited
over time. Each time Gaussian is deposited, a new line is written to the file
in the following format:

*center1 center2 ... width1 width2 ... height*

The centers denote the locations in CV space where each Gaussian has been
deposited, listed in the order in which the CVs appear within the SSAGES JSON
input file. The widths denote the corresponding Gaussian widths for each CV
dimension. The height is the Gaussian height, which should match the parameter
height defined in the JSON input file.

.. note::

	Although the widths and height of the Gaussian currently do not change in
	time, future additions to the Metadynamics method will allow for adaptive
	Gaussians.

Example MATLAB scripts are provided in the ``Examples/User/Meta`` directory.
These scripts sum the Gaussians and generate a free energy surface from
the "hills.out" file.

.. _metadynamics-tutorial:

Tutorial
^^^^^^^^

Two Metadynamics examples are included in the ``Examples/User/Meta`` directory.
In the first example, Metadynamics is used to sample the free energy surface of
a two-dimensional particle undergoing Langevin dynamics. This example is found in
the ``Single_Atom`` directory and uses LAMMPS.
The files included are described below:

* ``in.LAMMPS_Meta_Test``: LAMMPS input file describing the Langevin particle
  and underlying free energy surface to be sampled. The free energy surface
  consists of two Gaussian wells at (0.98, 0.98) and (-0.98, -0.98)
  respectively, and one Gaussian barrier at the origin.
* ``Meta.json``: SSAGES JSON input file specifying Metadynamics and CVs to be
  sampled. In this case the CVs are the *x* and *y* coordinates of the particle.
* ``analysis.m``: MATLAB script that analyzes the output of the Metadynamics
  method.
* ``Movie.m``: MATLAB script that generates a movie of the free energy
  surface estimate over time.

To run this example:

1. Either copy or create a symbolic link to the SSAGES executable in the
   examples directory.

.. code-block:: bash

	ln -s /path/to/SSAGES/build/ssages

2. Run the example by issuing the command below. Please note that in this
   example, two walkers are used to explore the system more efficiently. If
   you would like to use more walkers (1 processor per walker), simply include
   more drivers in the ``Meta.json`` input file.

.. code-block:: bash

	mpirun -np 2 ./ssages Meta.json

3. After the run is complete use the provided ``analysis.m`` or ``analysis.py``
   script to generate a representation of the underlying free energy surface using
   MatLab or python (Matplotlib and numpy required).

Developer
^^^^^^^^^

* Hythem Sidky
.. _references:

.. Developer's note: This file is named "zReferences.rst" to appear at the
   end of the file list, as sphinxcontrib-bibtex must have the references
   file compiled *last*.

References
==========

.. bibliography:: ../references.bib
	:style: plain
.. _basis-function-sampling:

Basis Function Sampling
-----------------------

Introduction
^^^^^^^^^^^^

The Basis Function Sampling method is a variant of the Continuous
Wang-Landau Sampling method developed by
Whitmer *et al.* :cite:`WHITMER2014190602`, which biases a PMF
through the summation of Kronecker deltas. In this method, the Kronecker delta
is approximated by projection of a locally biased histogram to a truncated set
of orthogonal basis functions.

.. math::

	\int_\Xi f_{i}(\vec{\xi}) f_{j}(\vec{\xi}) w(\vec{\xi}) d\vec{\xi} =
	\delta_{ij}c_{i}

By projecting a basis set, the system resolves the same properties as the
Kronecker deltas, but in a continuous and differentiable manner that lends well
towards MD simulations. The current version of SSAGES has support for Chebyshev,
Fourier, and Legendre polynomials. Each of these has their defined weight
function :math:`w(\xi)` implemented specific to the method. Additionally, any
combination of implemented basis sets can be used for any system. It is advised
that a periodic basis set (e.g. Fourier) be used with a periodic CV, but it is
not required.

The BFS method applies its bias in sweeps of :math:`N` through a histogram
(:math:`H_{i}`) that is updated at every :math:`j` microstate or timestep.
This histogram is then modified to an unbiased partition function estimate
(:math:`\tilde{H_{i}}`) by exponentiation with the current bias potential
(:math:`\Phi_{i}`).

.. math::

	\tilde{H}_{i}(\xi) = H_{i}(\xi) e^{\beta \Phi_{i}}

A weight function has been added into this implementation (:math:`W(t_{j})`) so that
the user can define the effective strength of the applied bias. If not chosen, the weight is normalized to
the length of the interval.

.. math::

	Z_{i}(\xi) = \sum_{j} W(t_{j})\tilde{H_{j}}(\xi)

This final estimate is then projected to the truncated basis set. After this set
is evaluated, the coefficients of the basis set are evaluated. This process is
iterated until the surface converges, which is determined by the overall update
of the coefficients.

.. math::

	\beta \Phi_{i+1}(\xi) &= \sum_j^N \alpha^i_j L_j(\xi)\\
	\alpha^i_j &= \frac{2j + 1}{2} \int_{-1}^1 \log(Z_i(\xi))L_j(\xi)d\xi

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

These are all the options that SSAGES provides for running Basis Function
Sampling. In order to add BFS to the JSON file, the method should be labeled as
``"BFSMethod"``.

Basis Function Sampling requires the use of a basis set. These are defined
by defining an object of "basis_functions". These have the following properties:

type
	Currently can either be Chebyshev, Fourier, or Legendre.

polynomial_order
	Order of the polynomial. In the case of Chebyshev or Legendre, this results
	in an order of input value + 1 as the method takes the 0th order internally.
	For a Fourier series, the order is the total number of coefficients
	including the sine and cosine series.

upper_bound
	Only exists for Chebyshev and Fourier series. This is the upper bound of the CV.

lower_bound
	Only exists for Chebyshev and Fourier series. This is the lower bound of the CV.

CV_restraint_spring_constants
	The strength of the springs keeping the system in bounds in a non-periodic
	system.

CV_restraint_maximums
	The upper bounds of each CV in a non-periodic system.

CV_restraint_minimums
	The lower bounds of each CV in a non-periodic system.

cycle_frequency
	The frequency of updating the projection bias.

frequency
	The frequency of each integration step. This should almost always be set to 1.

weight
	The weight of each visited histogram step. Should be kept around the same value
	as the ``cycle_frequency`` (usually 0.1 times that).

.. note::

	The system has a higher chance of exploding at higher weight values.

basis_filename
	A suffix to name the output file. If not specified, the output will be
	``basis.out``.

temperature
	The temperature of the simulation.

tolerance
	Convergence criteria. The sum of the difference in subsequent updates of the
	coefficients squared must be less than this for convergence to work.

convergence_exit
	A boolean option to let the user choose if the system should exit once the
	convergence is met.

Required to Run BFS
^^^^^^^^^^^^^^^^^^^

In order to use the method properly a few things must be put in the JSON file. A
grid is required to run Basis Function Sampling. Refer to the Grid section in
order to understand options available for the grid implementation.
The only inputs required to run the method:

* ``cycle_frequency``
* ``frequency``
* ``basis_functions``
* ``temperature``

Example Input
^^^^^^^^^^^^^
.. code-block:: javascript

	"methods": [{
		"type": "BFSMethod",
		"basis_functions": [
		{
				"type": "Fourier",
				"polynomial_order": 30,
				"upper_bound": 3.14,
				"lower_bound": -3.14
			},
			{
				"type": "Fourier",
				"polynomial_order": 30,
				"upper_bound": 3.14,
				"lower_bound": -3.14
			}
		],
		"cvs": [0, 1],
		"cycle_frequency": 100000,
		"basis_filename": "example",
		"frequency": 1,
		"temperature": 300.0,
		"weight": 1.0,
		"tolerance": 1e-3,
		"convergence_exit": true,
		"grid": {
			"lower": [-3.14, -3.14],
			"upper": [3.14, 3.14],
			"number_points": [100, 100],
			"periodic": [true, true]
		}
	}]

Guidelines for Running BFS
^^^^^^^^^^^^^^^^^^^^^^^^^^

* It is generally a good idea to choose a lower order polynomial initially.
  Excessive number of polynomials may create an unwanted
  `"ringing" effect <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_
  that could result in much slower convergence.
* For higher order polynomials, the error in projection is less, but the number
  of bins must increase in order to accurately project the surface. This may
  also create an unwanted
  `"ringing" effect <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.
* A good rule of thumb for these simulations is to do at least one order of
  magnitude more bins than polynomial order.

If the system that is to be used requires a non-periodic boundary condition,
then it is typically a good idea to place the bounds approximately 0.1--0.2
units outside the grid boundaries.

The ``convergence_exit`` option is available if the user chooses to continue
running past convergence, but a good heuristic for tolerance is around 0.001.

.. _BFS-tutorial:

Tutorial
^^^^^^^^

This tutorial will provide a reference for running BFS in SSAGES. There are
multiple examples provided in the ``Examples/User/BasicFunc`` directory of
SSAGES, but this tutorial will cover the Alanine Dipeptide example.

In the ``Examples/User/BasicFunc/ADP`` subdirectory, there should be two
LAMMPS input files (titled ``in.ADP_BFS_example{0,1}``) and two JSON input files.
Both of these files will work for SSAGES, but the one titled
``ADP_BFS_2walkers.json`` makes use of multiple walkers.

For LAMMPS to run the example, it must be made with ``rigid`` and
``molecule`` packages. In order to do so, issue the following commands from
your build directory:

.. code-block:: bash

	make yes-rigid
	make yes-molecule
	make

Use the following command to run the example:

.. code-block:: bash

	mpiexec -np 2 ./ssages ADP_BFS_2walkers.json

This should prompt SSAGES to begin the simulation. If the run is
successful, the console will output the current sweep number on each node.
At this point, the user can elect to read the output information after
each sweep.

**basis.out**

The ``basis.out`` file outputs in at least 3 columns. These columns refer to the
CV values, the projected PMF from the basis set, and the log of the histogram.
Depending on the number of CVs chosen for a simulation, the
number of CV columns will also correspond. Only the first CV column should be
labeled.

The important line for graphing purposes is the projected PMF, which is the
basis set projection from taking the log of the biased histogram. The biased
histogram is printed so that it can be read in for doing restart runs (subject to
change). For plotting the PMF, a simple plotting tool over the CV value and
projected PMF columns will result in the free energy surface of the simulation.
The free energy surface will return a crude estimate within the first few
sweeps, and then will take a longer period of time to retrieve the fully
converged surface. A reference image of the converged  alanine dipeptide example
is provided in the same directory as the LAMMPS and JSON input files.

**restart.out**

This holds all the coefficient values after each bias projection update, as well
as the biased histogram. This file is entirely used for restart runs.

Developers
^^^^^^^^^^

* Joshua Moller
* Julian Helfferich
.. _combined-force-frequency-sampling:

Combined Force-Frequency Sampling
---------------------------------

Introduction
^^^^^^^^^^^^

Combined Force-Frequency (CFF) sampling is a free energy sampling method which uses
artificial neural networks to generate an on-the-fly adaptive bias capable of rapidly
resolving free energy landscapes. It is a recent method proposed in :cite:`SEVGEN2020`,
and its extension and application to ab initio MD has been demonstrated in
:cite:`LEE2020`. The method's main strength resides in its ability to learn both from the
frequency of visits to distinct states and the generalized force estimates that arise in a
system as it evolves in phase space. This is accomplished by introducing a
self-integrating artificial neural network, which generates an estimate of the free energy
directly from its derivatives.

CFF algorithm proceeds in sweeps. Statistics are collected over a user specified interval
(``sweep``), which is a very flexible choice. At the end of each sweep artificial neural
networks are fit to the statistics to determine the optimal bias which is then applied in
the subsequent sweep. This proceeds until the free energy landscape has converged.

Example Input
^^^^^^^^^^^^^

.. code-block:: javascript

    "methods" : [
        {
            "type" : "CFF",
            "topology" : [12,8],
            "nsweep" : 10000,
            "temperature" : 298.15,
            "grid" : {
                "lower" : [-3.14159,-3.14159],
                "upper" : [3.14159,3.14159],
                "number_points" : [30,30],
                "periodic" : [false,false]
            },
            "lower_bounds" : [-5,-5],
            "upper_bounds" : [5,5],
            "lower_bound_restraints" : [0.0,0.0],
            "upper_bound_restraints" : [0.0,0.0],
            "timestep" : 0.001,
            "unit_conversion" : 1,
            "minimum_count" : 3000,
            "overwrite_output": true
        }
    ]

.. warning::

    Be sure to follow correct JSON syntax for your input, with a comma after every line
    except the last within each bracket.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

**Define CFF**

In the methods block, define the CFF method through the syntax:

.. code-block:: javascript

    "type" : "CFF",

To define neural network topology:

.. code-block:: javascript

    "topology" : [12,8],

Array of integers (length: number of hidden layers). This array defines the architecture
of the neural network. Like ANN Sampling, good heuristics is to use a network with a
single hidden layer if you are biasing on one CV, and two hidden layers if you are biasing
on two or more.

.. code-block:: javascript

    "temperature" : 298.15,

The temperature of the simulation must be specified.

To define the length of each sweep:

.. code-block:: javascript

    "nsweep" : 10000,

Typical values range from 1,000 to 10,000 depending on the size of the system. The slower
the system dynamics, the longer the sweep. This is not going to heavily affect convergence
time, and the method is generally quite robust to choice of sweep length. The main
consequence of this choice is that the neural network training time becomes relatively
expensive if the system's free energy is very cheap to evaluate.


**Define the grid**

To define the bounds:

.. code-block:: javascript

    "lower" : [-3.14159,-3.14159],
    "upper" : [3.14159,3.14159],

These are arrays of doubles whose length is the number of CVs used. This defines the
minimum and maximum values for the CVs for the range in which the method will be used in
order.

To define the number of CV bins used:

.. code-block:: javascript

    "number_points" : [30,30],

This array of integers defines the number of histogram bins in each CV dimension in order.

.. code-block:: javascript

    "periodic" : [false,false],

This array defines whether a given CV is periodic for restraint purposes. This is only
used to apply minimum image convention to CV restraints. The value can be safely set to
``false`` *even for periodic CVs* if no restraints are being used.

**Define the restraints**

.. code-block:: javascript

    "lower_bounds" : [-5,-5],
    "upper_bounds" : [5,5],

These arrays define the minimum and maximum values for the CV restraints in order.

.. code-block:: javascript

    "lower_bound_restraints" : [0,0],
    "upper_bound_restraints" : [0,0],

These arrays define the spring constant for the lower and upper bounds.

**Define time and unit parameters**

.. code-block:: javascript

    "timestep" : 0.001,

The timestep of the simulation. Units depend on the conversion factor that follows. This
must be entered correctly, otherwise the generalized force estimate will be incorrect.

.. code-block:: javascript

    "unit_conversion" : 1,

Defines the unit conversion from d(momentum)/d(time) to force for the simulation. For
LAMMPS using units real, this is 2390.06 (gram.angstrom/mole.femtosecond^2 ->
kcal/mole.angstrom). For GROMACS, this is 1.

.. code-block:: javascript

    "minimum_count" : 3000,

This is the number of hits required to a bin in the general histogram before the full
biasing force is active. Below this value, the bias linearly decreases to zero at
``hits = 0``. Default value is ``200``, but user should provide a reasonable value for
their system.

**Handle outputs (Optional)**

.. code-block:: javascript

    "overwrite_output" : [true],

If this is enabled, output files are overwritten at each sweep such.  Otherwise, output
files are saved at each sweep. Default, ``true``.

Output
^^^^^^

There are six output files from this method: ``CFF.out``, ``F_out``, ``netstate.net``,
``netstate.net2``, ``CFF.out_gamma``, and ``traintime.out``.

The main output of this method is stored in ``CFF.out``. Each column corresponds
respectively to the CVs, visit frequencies (histogram), bias based on frequency-based ANN,
bias based on force-based ANN, average bias, and the average free energy estimate (which
is the negative of the average bias with a constant shift). The format is as follows:

``cv1 cv2 ... hist bias(freq_only) bias(force_only) bias(avg) free_energy(avg)``

A file called ``CFF.out_gamma`` outputs network complexity term ``gamma`` for each neural
net and the ratio of gammas from both neural nets. (See :cite:`SEVGEN2020` for more
information.) The format is as follows:

``sweep_iter gamma(freq_only) gamma(force_only) gamma_ratio``

File ``traintime.out`` contains CPU wall time (in seconds) taken for training of neural
networks during each sweep.

Files called ``netstate.dat`` and ``netstate2.dat`` contain the neural network parameters
for each of the two neural networks used to apply biases during the sampling
(``netstate.dat`` stores frequency-based parameters, while ``netstate2.dat`` stores
force-based ones). For more information about these files, see the ANN Sampling Method.

A file called ``F_out`` contains the Generalized Force vector field, which is the same
output file as the Adaptive Biasing Force method. Vectors defined on each point on a grid
that goes from ``CV_lower_bounds`` to ``CV_upper_bounds`` of each CV in its dimension, with
``CV_bins`` of grid points in each dimension. The printout is in the following format: 2*N
number of columns, where N is the number of CVs. First N columns are coordinates in CV
space, the N+1 to 2N columns are components of the Generalized Force vectors. See the ABF
Sampling Method for more information.

An example of CFF (e.g., alanine dipeptide in water) is located in ``Examples/User/CFF/ADP``.

Developers
^^^^^^^^^^

* Elizabeth M.Y. Lee
* Emre Sevgen
* Boyuan Yu


.. warning::

    Please make sure to cite :cite:`SEVGEN2020` and :cite:`LEE2020` if you use this
    method!
Acknowledgments
================

We are grateful to Argonne National Laboratory for initiating this project and
their continued support. Julian Helfferich acknowledges financial support from
the DFG research fellowship program, grant No. HE 7429/1.

Some important core functionality of SSAGES comes from SAPHRON.
SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
https://github.com/hsidky/SAPHRON

Project Supervisors
-------------------

* Juan de Pablo
* Jonathan Whitmer

Project Leads
-------------

* Michael Quevillon
* Emre Sevgen

Past Project Leads
------------------

* Yamil J. Colón (2015-2018)
* Hythem Sidky (2015-2018)

Core Development
----------------

* Julian Helfferich
* Hythem Sidky
* Benjamin Sikora

Methods
-------

* Cody Bezik
* Yamil J. Colón
* Federico Giberti
* Ashley Guo
* Joshua Lequieu
* Jiyuan Li
* Joshua Moller
* Hadi Ramezani-Dakhel
* Emre Sevgen
* Hythem Sidky
* Benjamin Sikora
* Jonathan Whitmer

Collective Variables
--------------------

* Yamil J. Colón
* Ashley Guo
* Michael Quevillon
* Hythem Sidky
* Benjamin Sikora
* Mike Webb

Documentation
-------------

* Cody Bezik
* Yamil J. Colón
* Federico Giberti
* Ashley Guo
* Julian Helfferich
* Joshua Moller
* Michael Quevillon
* Hadi Ramezani-Dakhel
* Emre Sevgen
* Hythem Sidky

Engine Adapters
---------------

* Carl Simon Adorf
* Federico Giberti
* Bradley Dice
* Michael Quevillon
* Emre Sevgen
* Hythem Sidky
* Benjamin Sikora
.. _inputfiles:

Input Files
============

A SSAGES input file contains multiple sections that define the Collective
Variables (CVs), Methods, and various other components that go into an advanced
sampling simulation. There is a brief primer below on JSON, the format used by
SSAGES for input files. The remaining topics describe the basic syntax and
requirements for each section of an input file. Detailed information for
particular methods or collective variables can be found in their respective
locations in this manual.

JSON
----

SSAGES is run using input files written in the JSON_ file format. JSON is a
lightweight text format which is easy for humans to read and write and for
machines to parse and generate. Almost every programming language offers some
level of native JSON support, making it particularly convenient to script or
automate input file generation using, say, Python. If you've never used JSON
before, don't worry. Throughout the documentation we make no assumptions about
the user's knowledge of JSON and provide clear easy-to-follow examples.

A SSAGES input file is a valid JSON document. Here, we will define a bit of
terminology relating to JSON. Take the following JSON structure as an example,
obtained from Wikipedia_:

.. code-block:: javascript

	{
		"firstName": "John",
		"lastName": "Smith",
		"age": 25,
		"address": {
			"streetAddress": "21 2nd Street",
			"city": "New York",
			"state": "NY",
			"postalCode": "10021"
		},
		"phoneNumber": [
			{
				"type": "home",
				"number": "212 555-1234"
			},
			{
				"type": "fax",
				"number": "646 555-4567"
			}
		],
		"gender": {
			"type": "male"
		}
	}

The first pair of curly brackets define the *root* section, we will signify
this using ``#``. An item in the hierarchy, such as the street address, can be
referenced like this:  ``#/address/streetAddress``.

Square brackets ``[]`` in JSON refer to *arrays*, while curly brackets refer
to  *objects*. They can be thought of as Python lists and dictionaries
respectively. That would make ``#/phoneNumber`` an array of phone number
objects, each containing a type and a number. The fax number can be referenced
by ``#/phoneNumber/1/number``, where ``1`` is the array index beginning from
zero.

Items in a JSON object (Python dictionary) are unique. In the example above,
``#/age`` can only be defined once - it is a key in the root tree.  Defining
``#/age`` again will not throw an error, but instead the last definition will
override any previous definitions. This is actually very powerful behavior.
It means a user can import a general template JSON file and override whatever
parameters they wish. The exact behavior of the merging process is described in
detail in the user guide.

Types matter in JSON. Notice how ``#/age`` is specified by a number that is not
surrounded in quotes. This is a number, more specifically an integer. On the
other hand, ``#/address/postalCode`` is a string, even though the contents of
the string are all numbers. Certain fields in a SSAGES input file may be
required to be a string, integer, or number. The user should be aware of this
and take care to format their input file appropriately.

Simulation Properties
---------------------

A SSAGES build is compiled with support for a particular MD engine, and the
requirements for each engine vary slightly. For detailed information on
specific engines and their options check the :ref:`Engines <engines>` section.
The following parameters are needed to define a simulation in the JSON root.

.. warning::

	The properties specified below are case-sensitive. Please be sure to check
	that you have defined it according to the documentation.

Input
~~~~~

The ``"input"`` property specifies the name of the input file used by the
simulation engine.

.. code-block:: javascript

	"input": "in.system"

.. code-block:: javascript

	"input": ["in.system1","in.system2","in.system3"]

The first syntax is used if there is a single input file. For multi-walker
simulations, it is possible to use a single file for all walkers (though this
may not be recommended depending on the method) or specify a separate input
file for each walker.

.. note::

	This property not used by GROMACS (see ``"args"`` property).

Args
~~~~

.. warning::

	This property is *exclusively* for GROMACS and HOOMD-blue.

The ``"args"`` property specifies additional command line arguments to be
passed to the engine.

.. code-block:: javascript

	"args": ["-v", "-deffnm", "runfile"]

.. code-block:: javascript

	"args": "-v -deffnm runfile"

For GROMACS, a standard simulation can be invoked using
``gmx mdrun -deffnm runfile`` to execute a ``runfile.tpr`` binary, the
equivalent arguments must be specified in the ``"args"`` property. This
provides the user with the flexibility of calling command-line arguments in the
same fashion as the standard **mdrun** utility. The only exception is in the
case of multi-walker simulations. If a user wishes to use the multi-walker
capabilities, then ``"args"`` is invoked in the same fashion as a single-walker
simulation. **Do not specify the** ``-multi`` **option. This will be done
automatically.** If ``-deffnm`` is called, GROMACS expects the ``.tpr`` files
for each walker to  be named according to the walker ID starting from zero. In
the example above, if there were three walkers, then GROMACS will look for the
files "runfile0.tpr", "runfile1.tpr", and "runfile2.tpr".

Walkers
~~~~~~~

The ``"walkers"`` property specifies the number of walkers (independent instances
of the simulation engine) to run with SSAGES.

.. code-block:: javascript

	"walkers": 5

Many advanced sampling methods support multi-walker simulations which improve
the convergence of many algorithms. Typically, each walker has an independent
system configuration in a separate input file. It is important to note that
when specifying more than a single walker, the number of processors passed to
``mpiexec`` must be divisible by the number of walkers requested. Otherwise,
SSAGES will terminate with an error.

.. note::

	It is not possible to allocate a different number of processors to each
	walker, at this time.

Collective Variables
~~~~~~~~~~~~~~~~~~~~

The ``"CVs"`` property specifies the collective variables on which SSAGES
will perform its advanced sampling.

.. code-block:: javascript

	"CVs":
	[
		{
			"type": "Torsional",
			"name": "mytorsion_1",
			"atom_ids": [5,7,9,15]
		},
		{
			"type": "ParticleCoordinate",
			"atom_ids": [1],
			"dimension": "x"
		}
	]

Collective variables are specified in an array, where each element is a CV
object. Collective variables can be assigned names or referenced
by index, beginning with zero.

Methods
~~~~~~~

The ``"methods"`` property specifies the advanced sampling algorithms to which
SSAGES will apply to the system.

.. code-block:: javascript

	"methods":
	[
		{
			"type": "Umbrella",
			"ksprings": [100],
			"output_file": "ulog.dat",
			"output_frequency": 10,
			"centers": [1.0],
			"cvs": ["mytorsion_1"]
		},
		{
			"type": "Metadynamics",
			"widths": [0.3],
			"height": 1.0,
			"hill_frequency": 500,
			"lower_bounds": [0.2],
			"upper_bounds": [1.4],
			"lower_bound_restraints": [100],
			"upper_bound_restraints": [100],
			"cvs": [1]
		}
	]

Methods are specified in an array, since it is possible to run multiple methods
simultaneously. This is useful if a user is interested in performing advanced
sampling on a system subject to some restraint, typically applied via an
umbrella. Each method can selectively operate on a subset of CVs by referencing
them either by name or index, as shown above.

Logger
~~~~~~

The ``"logger"`` property specifies an output file to track any or all CVs
as the simulation proceeds.

.. code-block:: javascript

	"logger": {
		"frequency": 100,
		"output_file": "cvs.dat",
		"cvs": [0, 3]
	}

If your simulation is using multiple walkers, you must define an array of
``"output_file"`` that has the same number of filenames as number of walkers.
For instance, with two walkers, use the syntax below.

.. code-block:: javascript

	"logger": {
		"frequency": 100,
		"output_file": ["cvs_w0.dat","cvs_w1.dat"],
		"cvs": [0, 3]
	}

The logger is useful in tracking the evolution of the CVs over the course of an
advanced sampling calculation. Logging CVs can allow for post-simulation
reweighting, or indicate if there are sampling problems in the system being
studied. The frequency of logging the CVs can be specified and each walker in a
multi-walker simulation will have a separate output file. A user can choose to
selectively log individual CVs as well.

Putting It All Together
~~~~~~~~~~~~~~~~~~~~~~~

Combining the previous sections into a single input file yields the following
(purely hypothetical) example input for a LAMMPS simulation.

.. code-block:: javascript

	{
		"walkers": 2,
		"input": ["in.first", "in.second"],
		"CVs":
		[
			{
				"type": "Torsional",
				"name": "mytorsion_1",
				"atom_ids": [5,7,9,15]
			},
			{
				"type": "ParticleCoordinate",
				"atom_ids": [1],
				"dimension": "x"
			}
		],
		"methods":
		[
			{
				"type": "Umbrella",
				"ksprings": [100],
				"output_file": ["walker1.dat", "walker2.dat"],
				"output_frequency": 10,
				"centers": [1.0],
				"cvs": ["mytorsion_1"]
			},
			{
				"type": "Metadynamics",
				"widths": [0.3],
				"height": 1.0,
				"hill_frequency": 500,
				"lower_bounds": [0.2],
				"upper_bounds": [1.4],
				"lower_bound_restraints": [100],
				"upper_bound_restraints": [100],
				"cvs": [1]
			}
		]
	}

To execute this input file, assigning two processors per walker, one would call the command below.

.. code-block:: bash

	mpirun -np 4 ./ssages inputfile.json

.. _JSON: https://www.json.org/json-en.html
.. _Wikipedia: https://en.wikipedia.org/wiki/JSON#JSON_sample
The SSAGES Cookbook
===================

.. A collection of short solutions to common problems. Just like a FAQ.

A collection of resources relating to SSAGES.

MICCoM Summer School 2017
-------------------------

On July 17th--19th, 2017, `MICCoM <http://miccom-center.org/>`_ held a summer
school session focused on the codes developed within the center, which includes
SSAGES. The main page for this event is here:
`MICCoM Computational School 2017 <http://miccom-center.org/summer-school-2017/index.html>`_.

The slides from the relevant talks are included here:

* `Setting up SSAGES <http://miccom-center.org/images/summer_school_setting_up_your_ssages.pdf>`_
* `How to Add a New Collective Variable (CV) <http://miccom-center.org/images/MiCCoM_CV_talk_revealed.pdf>`_
* `Classical Molecular Dynamics and Sampling Methods <http://miccom-center.org/images/Classical%20MD%20and%20Sampling-v4.pdf>`_
* `Enhanced Sampling: String Methods & Forward Flux <http://miccom-center.org/images/enhanced_sampling-string_methods.pdf>`_
* `Free energy methods <http://miccom-center.org/images/free_energy_methods.pdf>`_
Contribute to SSAGES
====================

The SSAGES project is built by an inclusive and welcoming group of physicists,
chemists, and chemical engineers working on complex Molecular Dynamics (MD)
simulations employing advanced sampling techniques. Advanced sampling is an
exciting and rapidly developing field. Similarly, this project is designed to
facilitate the usage and implementation of a wide array of sampling methods.
We welcome you heartily to join us as we embark on this great adventure.

There are many ways to contribute to SSAGES, and you do not necessarily need
programming skills to be part of this project (even though they surely help).
But, if you decide to work on the code base, you will be happy to find that
SSAGES is designed to be easy to use and is just as easy to extend. We put a
high priority on maintaining a readable and clearly structured code base as
well as an inclusive community welcoming new ideas and contributions.

Here is a short summary of ideas how you can become part of SSAGES:

Reporting, Triaging, and Fixing Bugs
    No software is without errors, inconsistencies, and strange behaviors. Even
    with zero programming knowledge, you can help tremendously by reporting
    bugs or confirming issued bugs.
    :ref:`Read more... <report_bugs>`

Improving the SSAGES documentation
    We strive for SSAGES to have a detailed, yet comprehensive, documentation
    on what it does and how it does it. This should include concise
    introductions to methods, quick to learn tutorials, complete coverage of
    the intricacies of each method, and helpful pointers in case you run into
    errors. While the documentation is already expansive, improvements on it
    never go unappreciated.
    :ref:`Read more... <improve_documentation>`

Including your Method and CV in SSAGES
    You have developed a new sampling method or collective variable and
    want to make it available to the community via SSAGES? Great!
    :ref:`Read more... <add_your_method>`

Working on the core SSAGES system
    If you would like to climb into the heart of SSAGES and get your hands
    dirty, this task is for you! :ref:`Read more... <work_on_core_ssages>`

.. _report_bugs:

Reporting Bugs and Requesting Features
--------------------------------------

SSAGES is an open-source project maintained by a small team of scientific
developers. Because of this, we appreciate any help that the community can
provide!

If you have been using SSAGES and have come across an error that might affect
other users or the core code base, please report it on the `GitHub Issues`_
page within the repository. When reporting bugs, it is helpful to include as
much detail as possible. This can include, but is not limited to, input files
necessary for replicating the bug, the output you expected to see versus what
you got from SSAGES, and any relevant error messages. Before submitting a bug,
browse the current issues to make sure that it has not already been reported.

If you feel that there is a feature missing from SSAGES, you can either submit
a feature request (also via `GitHub Issues`_) or you can add the feature
yourself and submit a pull request on the `GitHub Pull Requests`_ page. When
submitting a feature request, please include a reference to relevant literature
for new methods/CVs/etc., as well as a generalized use case (or more) that can
motivate adding in the requested feature.

.. _GitHub Issues: https://github.com/SSAGESproject/SSAGES/issues
.. _GitHub Pull Requests: https://github.com/SSAGESproject/SSAGES/pulls

.. _improve_documentation:

Improving the Documentation
---------------------------

| Great documentation and great code produces great software.
| ---SSAGE advice
|

Improvements on the documentation are always highly appreciated.
The SSAGES documentation is split into two parts: The User Manual (which you
are reading right now) and the API documentation. While the User Manual uses
the `Sphinx`_ documentation and contains all information necessary to use the
program, the API docs are built on `Doxygen`_ to describe the usage of the
underlying classes and functions for everyone willing to extend and improve
SSAGES.

.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _Doxygen: https://www.doxygen.nl/index.html

Here are a few ideas on how you can help:

* Fix typos: Even though we are regularly checking, there are certainly still
  a few hidden somewhere.
* Check links: Verify all internal and external links are working.
* Bring up to date: Make sure that the documentation is current, i.e. that it
  reflects the usage of the latest version.
* Add examples: An example on how to use a method, avoid a common problem, etc.
  is more helpful than a hundred pages of dry descriptions.
* Write a tutorial: If there are aspects of the code that are missing a
  tutorial, and one might be useful, consider adding a short tutorial with a
  simple toy system.

Building the Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Before you can work on the documentation, you first have to build it. The
documentation is part of the SSAGES source code. It is assumed that you have
already downloaded and built the source code as described in the
:ref:`Getting Started <getting-started>` section. You will find a
collection of ``.rst`` files comprising the User Manual under ``doc/source/``
where the file ending ``.rst`` stands for ReStructured Text. The API
documentation, on the other hand, resides directly in the header files, right
next to the classes and functions they describe.

Assuming you have already built SSAGES, building the documentation is as easy as
typing

.. code:: bash

	``make doc``

within your build directory. In order to make the documentation, you must have
the following programs installed:

* Sphinx (with PyPI via ``pip install Sphinx`` for example)
* Doxygen
* dot (in Ubuntu this is part of the graphViz package)
* Sphinx "`Read the docs`_" theme (via ``pip install sphinx_rtd_theme``)

.. _Read the docs: https://github.com/readthedocs/sphinx_rtd_theme

Once you have successfully built the documentation, you will find the User Manual
under ``doc/Manual/`` and the API documentation under ``doc/API-doc/html/``
(relative to your build directory - do not confuse it with the ``doc/`` folder
in the main directory of the project). To view it in your favorite web
browser (using FireFox as an example) just type

``firefox doc/Manual/index.html``

for the User Manual or

``firefox doc/API-doc/html/index.html``

for the API documentation.

Writing Documentation
^^^^^^^^^^^^^^^^^^^^^

Here are a few pointers on how to write helpful documentation, before we dive
into the details of **Sphinx** and **Doxygen** for the User Manual and the API
documentation:

* Write documentation "along the way". Do not code first and write the
  documentation later.
* Use helpful error messages. These are considered part of the documentation and
  probably are the part that is read most frequently.
* Do everything you can to structure the text. Let's face it: most people will
  just skim the documentation. Feel encouraged to use any and all techniques that
  help to spot the relevant information, for example:

  * Format your text **bold**, *italic*, ``code``, etc.
  * Use effective headers
  * Write in short paragraphs
  * Use lists, code blocks, tables, etc.

  .. note::

    These Note blocks are extremely helpful for example.

  .. warning::

    Warnings work great, too!

  .. seealso::

    Here you can find more examples for helpful Sphinx markup:
    
    https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

* Use examples, a lot of them.
* In the initial stages: Don't be a perfectionist. Missing documentation is
  the worst kind of documentation.

| "It is better to have written and coded than to have never written at all."
| ---SSAGE advice
|

Documenting with Sphinx
~~~~~~~~~~~~~~~~~~~~~~~

The **Sphinx** documentation system uses reStructuredText which is loosely
based on the Markdown format. Examples for documentations written with Sphinx
include:

* `LAMMPS`_
* `HOOMD-blue`_
* Virtually all of the `Python`_ Documentation

The following tutorials are extremely helpful:

* https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
* https://docutils.sourceforge.io/docs/user/rst/quickref.html
* http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html

.. _LAMMPS: https://lammps.sandia.gov/doc/Manual.html
.. _HOOMD-blue: https://hoomd-blue.readthedocs.io/en/stable/index.html
.. _Python: https://docs.python.org/3/

One of the great things of Sphinx is that most documentations have a "view page
source" link at the top of the page, where you can take a look at the Sphinx
source code. Thus, the best way to learn Sphinx is to click on this link right
now and look at the source code of this page. But here is a short summary of
the most important commands:

* Markup: You can use \*italic*, \**bold**, and \``code`` for *italic*, **bold**
  and ``code``.
* Headers: Underline your headers with at least three ``===`` for titles,
  ``---`` for subtitles, ``^^^`` for subsubtitles and ``~~~`` for paragraphs.
* Lists: Bulleted lists are indicated by lines beginning with ``*``.

.. note::

    These highlighted blocks can be created with ``.. note::``. The content of
    this block needs to be indented. You can also use ``warning`` and
    ``seealso``. Even more can be found
    `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.

Documenting with Doxygen
~~~~~~~~~~~~~~~~~~~~~~~~

**Doxygen** follows a very different philosophy compared to Sphinx and is more
steered towards API documentation, exactly what we use it for in SSAGES.
Instead of maintaining the documentation separately from the source code, the
classes and functions are documented in the same place where they are declared:
the header files. Doxygen then reads the source code and automatically builds
the documentation. Examples for documentation created with Doxygen include:

* `PLUMED`_
* `Root`_

.. _PLUMED: https://plumed.github.io/doc-v2.3/user-doc/html/index.html
.. _Root: https://root.cern.ch/doc/master/index.html

The mainpage of the Doxygen documentation is written in a separate header file,
in our case ``doc/mainpage.h``. A good introduction to the Doxygen syntax can
be found at

* https://www.doxygen.nl/manual/docblocks.html

The basic rule is that Doxygen comments start with ``//!`` or ``/*!`` and
document the class, namespace or function that directly follows it. Let's start
with a short example:

.. code-block:: cpp

    //! Function taking the square of a value
    /*!
     * \param val Input value
     * \returns Square of the input value
     *
     * This function calculates the square of a given value.
     */
    double square(double val)
    {
        return val*val;
    }

This example documents the function ``square()`` which simply calculates the
square of a number. The first line, starting with ``//!``, is the brief
description and should not be longer than one line. The second comment block,
starting with ``/*!`` is the full description. Here, two special commands
are used:

\\param
    This command documents one parameter of the function.

\\returns
    This command documents the return value of the function.

There are many special Doxygen commands. They all start with a backslash and
the most important, apart from the two mentioned above, are:

\\tparam
    Used to document a template parameter.

\\ingroup
    This class is part of a group, such as Methods or Core. The groups are
    defined in ``doc/mainpage.h``.

There are also helpful boxes that highlight a given aspect of the function, such as:

\\attention
    Puts the following text in a raised box. A blank line ends the attention box.

\\note
    Starts a highlighted block. A blank line ends the note block.

\\remark
    Starts a paragraph where remarks may be entered.

\\see
    Paragraph for "See also".

\\deprecated
    The documented class or function is deprecated and only kept for backwards
    compatibility.

\\todo
    Leave a "to do" note with this command.

You can also highlight your text:

\\em
    For an *italic* word. To apply to multiple words, use <em> *Italic text* </em>.

\\b
    For a **bold** word. To apply to multiple words, use <b> **Bold text** </b>.

\\c
    For a ``code`` word (in typewriter font). To apply to multiple words, use
    <tt> ``Typewriter`` </tt>.

\\code
    Starts a ``code`` block. The block ends with **\\endcode**.

\\li
    A line starting with **\\li** is an entry in a bullet list.

Another benefit of Doxygen is that you can use some LaTeX syntax. For example:

\\f$
    Starts and ends an inline math equation, similar to $ in Latex.

\\f[ and \\f]
    Start and end a display-style LaTeX equation.

\\cite <label>
    Cite a reference. The references are listed in ``doc/references.bib`` and
    follow BibTex syntax.

Doxygen is very clever in producing automatic links. For example, there
exists a class ``Method`` in SSAGES. Thus, Doxygen automatically creates a
link to the documentation of this class where the word "Method" appears. This,
however, does not work for the plural, "Methods". Instead, you can write
``\link Method Methods \endlink``. On the other hand, if you want to prevent
Doxygen from creating an autolink, put a ``%`` in front of the word.

What to Document
^^^^^^^^^^^^^^^^

We strive for a comprehensive documentation of all the methods available in
SSAGES, as well as the core features. Thus, for each method, the documentation
should include:

* An introduction to the method, what it does, and how it does it.
* A short tutorial based on one of the working examples. The reader should be
  able to complete the tutorial in ~30 min and should leave with a sense of
  accomplishment, e.g. a nice energy profile or a picture of a folded protein.
* A detailed description on how to use the method, the parameters, constraints,
  requirements, etc.

.. _add_your_method:

Adding your Method to SSAGES
----------------------------

.. seealso::

    See :ref:`here <Write-your-own-method>` for an introduction to how to
    develop your own method.

So, you have developed a new sampling method or new collective variable (CV)?
Great! SSAGES is about collaboration and integrating your new CV or method is
a priority. But before we do that, make sure you check the following boxes:

* Your code needs to compile and run (obviously).
* If you have implemented a new method, this method should have been published
  in a peer reviewed journal and the publication should be cited in the
  documentation of the method (see next point). If you have implemented a CV,
  please give a small example of usage. In which case(s) does the new CV come
  in handy?
* Your method needs to come with the necessary documentation. For others to
  be able to use your method, you will have to explain how it works. You can
  take a look at the section :ref:`"Improving the Documentation"
  <improve_documentation>` for a starter on how to write good documentation.
* Please provide an example system. This could be the folding of an
  Alanine Dipeptide molecule, a NaCl system, or just a toy model with a simple
  energy landscape. As long as the system is small and the method can easily
  complete within a few hours, it will be fine.

Once these boxes have been checked, our team of friendly developers will review
your source code to help you meet the standards of the SSAGES code.

.. _work_on_core_ssages:

Working on the Core Classes
---------------------------

.. todo::

    Describe SSAGES code hierarchy and development workflow.
.. _advanced-sampling-methods:

Advanced Sampling Methods
=========================

List of Advanced Sampling Methods (alphabetical order):

.. toctree::
   :maxdepth: 1

   Adaptive Biasing Force Algorithm
   Artificial Neural Network Sampling
   Basis Function Sampling
   Combined Force-Frequency Sampling
   Elastic Band
   Finite Temperature String
   Forward-Flux
   Metadynamics
   Swarm of Trajectories
   Umbrella Sampling

.. The following are :orphan: pages
..   Image Method
..   Replica Exchange
.. _adaptive-biasing-force:

Adaptive Biasing Force Algorithm
--------------------------------

Introduction
^^^^^^^^^^^^

Relative to most other free energy methods, which are based on adding biases
to the Hamiltonian of a system to obtain a flat histogram in converged
sampling, the adaptive biasing force method (ABF) is unique, as it seeks a
flattening of the generalized force. Like flat histogram methods, this bias
force is applied until the landscape seen by a simulation is flat in
collective variable (CV) space and freely diffusive sampling is achieved.

In practice, ABF functions by partitioning relevant CV space into histogram-
like bins with a grid, and keeping a running tabulation of the number of
visits to each bin. Concurrently, a running sum of the instantaneous
generalized force experienced by the system. Putting these together gives an
estimate for the mean generalized force, and integrating this quantity yields
the free energy. Note that this means ABF returns a vector field rather than a
free energy surface (which is the output of most methods.) Ref.
:cite:`COMER20141129` contains
an excellent write-up on the method in general. The specific implementation in
SSAGES is discussed in Ref. :cite:`DARVE2008144120`.

.. note::

    An integrator for 1D, 2D and 3D surfaces are provided in SSAGES/Tools/ABF_integrator (requires numpy, scipy and matplotlib). Syntax is below; this is illustrated further in the tutorial found on this page.

.. code-block:: bash

    ./ABF_integrator.py -i <inputfile> -o <outputname> -p <bool> (<bool> <bool>) --interpolate <integer> (<integer> <integer>) --scale <float>

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

ABF calculates a generalized force on CVs at each timestep, and biases
simulations using the negative of the estimated generalized force. This
requires a grid, which requires that you define a CV range. Outside of the CV
range, the simulation will continue to run, but the grid does not extend
there, so there is no bias applied and no histogram hits are collected.
histogram hits will be collected.


ABF can restart from a previous run. Simply include Fworld_cvX and Nworld
outputs generated by the previous run in your working directory, and set the
``"restart"`` option to ``true``.

.. warning::

  ABF keeps only a single backup of old files, and overwrites older backups
  with new data if they already exist. If restarting, a backup will not be
  created, instead ABF will read from and update the newest files. If you want
  to keep a copy of the output from a previous run after restart, be sure to
  rename the output or place a copy in a different directory.


If you are using the multiple walkers option, they will read from and write to
the same histogram and force estimate during runtime. The resulting histogram
and force data are saved in the same way as single-walker simulations.

ABF can optionally define a restraint range, which biases simulations back
toward the region of interest using a harmonic restraint with user-chosen
spring constant(s). To disable restraints, enter a spring constant k equal to
or less than zero.

.. warning::

  The restraint range should be WIDER than the CV range by at least one bin size
  in each direction. 

If restraints are used on a periodic system, one can define the periodic
boundaries, so that minimum image convention to CVs can be applied using the
commands ``CV_periodic_boundary_upper_bounds`` and
``CV_periodic_boundary_lower_bounds``. For example, on a :math:`-\pi` to
:math:`\pi` CV, if the CV is restrained between -3.14 to -2.36 and the system
crosses the periodic boundary, setting this will ensure the restraint is
applied correctly back towards -3.14 rather than applying an incorrect a large
force to push it toward -2.36.

Example Input
^^^^^^^^^^^^^

.. code-block:: javascript

    "methods" : [
        {
            "type" : "ABF",
            "cvs" : [0,1],
            "CV_lower_bounds" : [-3.14, -3.14],
            "CV_upper_bounds" : [3.14, 3.14],
            "CV_bins" : [21,21],
            "CV_restraint_minimums" : [-5,-5],
            "CV_restraint_maximums" : [5,5],
            "CV_restraint_spring_constants" : [0,0],
            "CV_isperiodic" : [false,false],
            "timestep" : 0.002,
            "minimum_count" : 50,
            "output_file" : "F_out",
            "output_frequency" : 1000,
            "unit_conversion" : 1,
            "frequency" : 1
        }
    ]

.. warning:: 

    Be sure to follow correct JSON syntax for your input, with a comma after every line except the last within each bracket.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

**Define ABF**:

In the methods block, define the ABF method through the syntax: 

.. code-block:: javascript

    "type" : "ABF"

**Define CVs**

To define the collective variables:

.. code-block:: javascript 

   "cvs" : [0,1]

In the example input, this defines a two-dimensional CV-space to be sampled by ABF, with indices [0,1]. The argument to this must be a list of integers defining the CVs to be operated on by ABF. 

**Define the grid**

To define the bounds:

.. code-block:: javascript

    "CV_lower_bounds" : [-3.14, -3.14] 
    "CV_upper_bounds" : [3.14, 3.14]

Thee are arrays of doubles whose length is the number of CVs used. This
defines the minimum and maximum values for the CVs for the range in which the
method will be used in order.

To define the number of CV bins used:

.. code-block:: javascript

    "CV_bins" : [21,21]

This array of integers defines the number of histogram bins in each CV dimension in order.


**Define the restraints**

.. code-block:: javascript

    "CV_restraint_minimums" : [-5,-5],
    "CV_restraint_maximums" : [5,5],

These arrays define the minimum and maximum values for the CV restraints in order. 

.. code-block:: javascript  

    "CV_restraint_spring_constants" : [0,0],

This array defines the spring constants for the CV restraints in order.
Enter a value equal to or less than zero to turn restraints off.

.. code-block:: javascript  

    "CV_isperiodic" : [false,false],

    This array defines whether a given CV is periodic for restraint purposes.
    This is only used to apply minimum image convention to CV restraints. The
    value can be safely set to ``false`` *even for periodic CVs* if no
    restraints are being used.

.. warning::

    If ANY CV is set to periodic, then ``CV_periodic_boundary_lower_bounds``
    and ``CV_periodic_boundary_upper_bounds`` must be provided for ALL CVs.
    Values entered for non-periodic CVs are not used.

.. code-block:: javascript  
    
    "CV_periodic_boundary_lower_bounds" : [-3.14, -3.14],
    "CV_periodic_boundary_upper_bounds" : [3.14, 3.14],

These arrays define the lower and upper end of the period. This only matters if
``CV_isperiodic`` is ``true`` for the CV.


**Define time and unit parameters**

.. code-block:: javascript

    "timestep" : 0.002,

The timestep of the simulation. Units depend on the conversion factor that
follows. This must be entered correctly, otherwise the generalized force estimate
will be incorrect.

.. code-block:: javascript

    "unit_conversion" : 1,

Defines the unit conversion from d(momentum)/d(time) to force for the
simulation. For LAMMPS using units real, this is
2390.06 (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom).
For GROMACS, this is 1.

.. code-block:: javascript

    "minimum_count" : 50,

This is the number of hits required to a bin in the general histogram before
the full bias is active. Below this value, the bias linearly decreases to
equal 0 at hits = 0. Default = 200, but user should provide a reasonable
value for their system. See :cite:`COMER20141129` and :cite:`DARVE2008144120`
for more details.

**Output parameters**

.. code-block:: javascript

    "output_frequency" : 1000,

*Optional*: This defines how many timesteps pass in between output of the generalized force.

.. code-block:: javascript

    "output_file" : "F_out",
    
This is a string value defining the file name for the adaptive vector force
field that is acquired. The default name is "F_out".

.. code-block:: javascript
    
    "Fworld_output_file" : "Fworld_cv"

*Optional*: This is the name of the file to backup raw Fworld force output for
use in restarts. There will be separate outputs for each CV. The default
filename is Fworld_cv, which saves each CV's output to Fworld_cvX.

.. code-block:: javascript
    
    "Nworld_output_file" : "Nworld"

*Optional*: This is name of the file which backs up the raw histogram data for
restart purposes. The default filename is "Nworld".

**Optional Parameters**

.. code-block:: javascript

    "mass_weighting" : false,

Turns on/off mass weighing of the adaptive force. The default is ``false``,
which turns off the weighting.

.. warning::

    Leave this off if your system has massless sites such as in TIP4P water.


.. code-block:: javascript

    "restart" : false

This boolean determines whether the simulation is a restart. The default value
is ``false``. If set to ``true``, ABF will attempt to load a previous state
from Nworld and Fworld files.

.. code-block:: javascript

    "frequency" : 1  

Leave at 1. 


Output
^^^^^^

The main output of the method is stored in a file specified in 'filename'. This 
file will contain the Adaptive Force vector field printed out every 
'backup_frequency' steps and at the end of a simulation. The method outputs a vector 
field, with vectors defined on each point on a grid that goes from 
(CV_lower_bounds) to (CV_upper_bounds) of each CV in its dimension, with (CV_bins) of grid points 
in each dimension. For example, for 2 CVs defined from (-1,1) and (-1,0) with 3 and
2 bins respectively would be a 3x2 grid (6 grid points). The printout is in the
following format: 2*N number of columns, where N is the number of CVs. First N columns 
are coordinates in CV space, the N+1 to 2N columns are components of the Adaptive Force 
vectors. An example for N=2 is:

+-----------+-----------+-------------+-------------+
| CV1 Coord | CV2 Coord | d(A)/d(CV1) | d(A)/d(CV2) |
+===========+===========+=============+=============+
| -1        | -1        | -1          | 1           |
+-----------+-----------+-------------+-------------+
| -1        | 0         | 2           | 1           |
+-----------+-----------+-------------+-------------+
| 0         | -1        | 1           | 2           |
+-----------+-----------+-------------+-------------+
| 0         | 0         | 2           | 3           |
+-----------+-----------+-------------+-------------+
| 1         | -1        | 2           | 4           |
+-----------+-----------+-------------+-------------+
| 1         | 0         | 3           | 5           |
+-----------+-----------+-------------+-------------+

.. _ABF-tutorial:

Tutorial
^^^^^^^^

Alanine Dipeptide

For LAMMPS (must be built with RIGID and MOLECULE packages)
To build RIGID and MOLECULE: 

1) Go to LAMMPS src folder (/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/ for -DLAMMPS=YES)
2) Do:

.. code-block:: bash

   make yes-RIGID
   make yes-MOLECULE

3) Go to your build folder and make.

Find the following input files in Examples/User/ABF/Example_AlanineDipeptide:

* ``in.ADP_ABF_Example(0-1)`` (2 files)
* ``example.input``
* ``ADP_ABF_1walker.json``
* ``ADP_ABF_2walkers.json``

1) Put the contents of ABF_ADP_LAMMPS_Example folder in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ADP_ABF_1walker.json.json
    
For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ADP_ABF_2walkers.json

For GROMACS:

Optional:

* ``adp.gro``
* ``topol.top``
* ``nvt.mdp``

Required:

* ``example_adp(0-1).tpr`` (2 files)
* ``ADP_ABF_1walker.json``
* ``ADP_ABF_2walkers.json``

1) Put the contents of ABF_ADP_Gromacs_Example in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ABF_ADP_1walker.json

For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ABF_ADP_2walkers.json

These will run using the pre-prepared input files in .tpr format. If you wish to
prepare the input files yourself using GROMACS tools (if compiled with -DGROMACS=YES):

.. code-block:: bash

    /build/hooks/gromacs/gromacs/bin/gmx_ssages grompp -f nvt.mdp -p topol.top -c adp.gro -o example_adp0.tpr
    /build/hooks/gromacs/gromacs/bin/gmx_ssages grompp -f nvt.mdp -p topol.top -c adp.gro -o example_adp1.tpr

Be sure to change the seed in .mdp files for random velocity generation, 
so walkers can explore different places on the free energy surface.

Multiple walkers initiated from different seeds will
explore different regions and will all contribute to the same adaptive force.

After the run is finished, you can check that your output matches the sample
outputs given in the examples folders:

1) Copy ABF_integrator.py (requires numpy, scipy and matplotlib) into your build folder.
2) Run the integrator:

.. code-block:: bash

    python ABF_integrator.py --periodic1 True --periodic2 True --interpolate 200

3) This will output a contour map, a gradient field and a heatmap. Compare these to the sample outputs.




Sodium Chloride

For LAMMPS (must be built with KSPACE and MOLECULE packages)
To build RIGID and MOLECULE: 

1) Go to LAMMPS src folder (/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/ for -DLAMMPS=YES)
2) Do:

.. code-block:: bash

   make yes-KSPACE
   make yes-MOLECULE

3) Go to your build folder and make.

Find the following input files in Examples/User/ABF/Example_NaCl/ABF_NaCl_LAMMPS_Example:

* ``in.NaCl_ADP_example(0-1)`` (2 files)
* ``data.spce``
* ``ADP_NaCl_1walker.json``
* ``ADP_NaCl_2walkers.json``

1) Put the contents of ABF_NaCl_LAMMPS_Example folder in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ADP_NaCl_1walker.json.json
    
For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ADP_NaCl_2walkers.json

For GROMACS:

Optional:

* ``NaCl.gro``
* ``topol.top``
* ``npt.mdp``

Required:

* ``example_NaCl(0-1).tpr`` (2 files)
* ``ADP_NaCl_1walker.json``
* ``ADP_NaCl_2walkers.json``

1) Put the contents of ABF_NaCl_Gromacs_Example in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ABF_NaCl_1walker.json

For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ABF_NaCl_2walkers.json

These will run using the pre-prepared input files in .tpr format. If you wish to
prepare the input files yourself using GROMACS tools (if compiled with -DGROMACS=YES):

.. code-block:: bash

    /build/hooks/gromacs/gromacs/bin/gmx_ssages grompp -f npt.mdp -p topol.top -c NaCl.gro -o example_NaCl0.tpr
    /build/hooks/gromacs/gromacs/bin/gmx_ssages grompp -f npt.mdp -p topol.top -c NaCl.gro -o example_NaCl1.tpr

Be sure to change the seed in .mdp files for random velocity generation, 
so walkers can explore different places on the free energy surface.

Multiple walkers initiated from different seeds will
explore different regions and will all contribute to the same adaptive force.

After the run is finished, you can check that your output matches the sample
outputs given in the examples folders:

1) Copy ABF_integrator.py (requires numpy, scipy and matplotlib) into your build folder.
2) Run the integrator:

.. code-block:: bash

    python ABF_integrator.py

3) This will output a Potential of Mean Force graph. Compare this to the sample output.


Developers
^^^^^^^^^^

* Emre Sevgen
* Hythem Sidky
.. swarm:

Swarm of Trajectories
---------------------

Introduction
^^^^^^^^^^^^

Like all string methods in general, the **String Method with Swarms of
Trajectories** (often abbreviated to "Swarm of Trajectories" or even more simply
"SoT") is a method to identify a transition pathway in an arbitrarily
high-dimensional collective variable space between two metastable states of a
system. This pathway (the string) is a parametrized curve discretized into a
set of images, each of which is itself a molecular system. The classical *String
Method in Collective Variables* :cite:`MARAGLIANO2006024106`
evolves each image by estimating a mean force and
metric tensor at each image with restrained molecular dynamics simulations. In
the SoT method, the string is instead evolved by launching a large number (a
swarm) of **unrestrained** trajectories from each image and estimating the average
drift of the collective variables over the swarm.

The mathematical background of the method can be expressed in a few relatively
straightforward equations, with further detail available in the original work of
Benoît Roux and collaborators :cite:`PAN20083432`.
First, consider a path :math:`z(\alpha)`
constructed between two metastable states, such that :math:`\alpha=0` represents
the starting state and :math:`\alpha=1` is the final state. The "Most Probable
Transition Pathway" (MPTP) is defined such that a molecular system started from
anywhere on the path will most probably evolve while staying on the path. It is
shown in the original work that a mathematical definition for such a path is
given when the collective variables evolve according to:

.. math::

	z_{i}(\alpha) = z_{i}(\alpha') + \sum\limits_{j}\left(
	\beta D_{ij}\left[ z(0) \right] F_{j}\left[z(0)\right] +
	\frac{\partial}{\partial z_{j}}\left( D_{ij}\left[z(0)\right]\right)
	\right)\delta\tau

Where the following notation is used: :math:`z_{i}` represents the collective
variables belonging to the string, :math:`\alpha` represents the parameter
identifying that point on the string, :math:`\beta` represents the temperature,
:math:`D_{ij}` represents the diffusion tensor, :math:`F_{j}` represents the
mean force, :math:`z` represents the collective variables constructed from the
molecular system at a given moment in time, and :math:`\delta\tau` represents
the time step of the evolution of the dynamics. The SoT method approximates
this equation using the average drift evaluated from a large number of unbiased
trajectories, each of length :math:`\delta\tau`, launched from each image:

.. math::

	\overline{\Delta z_{i}(\delta\tau)} =
	\overline{z_{i}(\delta\tau) - z_{i}(0)} \equiv
	\sum\limits_{j} \left( \beta D_{ij}\left[z(0)\right] F_{j}\left[z(0)]\right] +
	\frac{\partial}{\partial z_{j}}\left( D_{ij}\left[
	z(0)\right]\right)\right)\delta\tau

Like all string methods, there is an additional step beyond evolving the
collective variables - after one iteration of evolution, the images along the
path must be reparameterized such that they lie (for example) an equal arc length
apart. This step is necessary to ensure that all images do not fall into one
metastable basin or the other.

Algorithmically, the SoT method is implemented as follows:

1. An initial string is defined between the two states of interest. This can be
   defined however one wishes; often it is simply a linear interpolation through
   the space of the collective variables. In fact, the ends of the string need
   not necessarily be in the basins of interest; the dynamic nature of the
   method should allow the ends to naturally fall into nearby metastable basins.

2. For each image of the string, a molecular system with atomic coordinates that
   roughly correspond to the collective variables of that image is constructed.

3. A set of equilibrium trajectories are generated from that system by performing
   restrained sampling around the image’s collective variables.

4. That set of equilibrium trajectories is used as the starting point of a large
   number of short unbiased trajectories; the resulting average displacement of
   each collective variable is used to update the positions of the images.

5. A reparameterization scheme is enforced to ensure that, for example, the
   string images are equally distant in collective variable space.

Steps 2--5 are iterated upon, leading to convergence of the method
and the MPTP.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

To construct a Swarm input file, the following options are available. A
complete Swarm JSON file will inherit some of its inputs from the String
schema (for parameters common to all string methods).
The options unique to Swarm are:

initial_steps
	For each iteration of the method, this is the number of steps to spend
	doing restrained sampling and not harvesting trajectories. This time is
	important to ensure the underlying molecular system’s CV values are
	close to the string CV values.

harvest_length
	After the initial restraining is finished, a trajectory is harvested for
	later use in launching an unrestrained trajectory every so often - harvest
	length specifies how often this will be done. Harvest length multiplied by
	number of trajectories (see below) will determine overall how many more
	steps will be taken under restrained sampling.

number_of_trajectories
	The total number of unrestrained trajectories to be included in each swarm.

swarm_length
	The length of each unrestrained trajectory in the swarm. Swarm length
	multiplied by number of trajectories specifies how many total steps will be
	spent doing unrestrained sampling.

From the String schema, the options are:

type
	This parameter identifies that a String-type method is being used, and
	thus should be set to "String".

flavor
	This parameter identifies the specific kind of string-type method
	being used; for swarm, it should be set to "Swarm".

centers
	The initial values of each CV for each image on the string are specified
	under "centers", which is an array of size equal to the total number of
	images, with each entry consisting of an array with size equal to the
	number of CVs used for the string method. In this way, the initial string
	is defined.

tolerance
	This is a tolerance threshold that can be set to trigger the end of
	the method; it is a percentage by which, if no node CV changes by this
	percentage, the method will end. It must be specified as an array with
	one entry for each CV desired.

max_iterations
	A complementary stopping criterion can be specified; the method will
	stop if it undergoes this many iterations of the string method.

ksprings
	A unique spring constant must be defined for each CV; its purpose is
	described above.

frequency
	The frequency of each integration step. This should almost always be
	set to 1.

.. _Swarm_tutorial:

Tutorial
^^^^^^^^

This tutorial will walk you step by step through the user example provided with
the SSAGES source code that runs the SoT method on the alanine dipeptide using
LAMMPS. First, be sure you have compiled SSAGES with LAMMPS. Then, navigate to
the ``SSAGES/Examples/User/Swarm/ADP`` subdirectory. Now, take a moment to
observe the ``in.ADP_Test and data.input`` files. In general, these should be
the same as what you would use for any other method, but for the SoT method, it
is important to define a larger skin distance than one normally would in the
neighbor command in LAMMPS. This is because, under the hood, each unrestrained
trajectory in the swarm is started by manually resetting the positions of each
atom in the LAMMPS simulation to the start of a new trajectory. From the
perspective of LAMMPS, this is a huge amount of distance to move in a single
time step; this move triggers neighbor list rebuilding, but LAMMPS considers it
a "dangerous build" which threatens to crash the simulation. Thus, we increase
the skin distance, which forces LAMMPS to keep track of more pairs in the
neighbor lists, and thus reduces the number of dangerous builds. Keep this in
mind for future runs of the SoT method.

The next two files of interest are the ``Template_Input.json`` input file and
the ``Input_Generator.py`` script. Both of these files can be modified in your
text editor of choice to customize the inputs, but for this tutorial, simply
observe them and leave them be. ``Template_Input.json`` contains all the
information necessary to fully specify one driver; ``Input_Generator.py`` copies
this information a number of times specified within the script (for this
tutorial, 22 times) while also linearly interpolating through the start and end
states defined in the script and substituting the correct values into the
"centers" portion of the method definition. Execute this script as follows:

.. code-block:: bash

	python Input_Generator.py

You will produce a file called ``Swarm.json``. You can also open this file to
verify for yourself that the script did what it was supposed to do. Now, with
your JSON input and your SSAGES binary, you have everything you need to perform
a simulation. Simply run:

.. code-block:: bash

	mpiexec -np 22 ./ssages Swarm.json

Soon, the simulation will produce a ``node-X.log`` file for each driver, where
X is the number specifying the driver (in this case, 0-21 for our 22 drivers).
Each one will report the following information, in order: the node number, the
iteration number, and for each CV, the current value of the string CV as well as
the current value of the CV calculated from the molecular system.

Allow your system to run for the desired number of MD steps, but keep an eye on
it - the system should exit once one driver reaches the maximum number of MD
steps, but it is possible that instead one driver will exit and the rest will
get stuck. Check in on your node files and see if they’ve been updated recently - if
not, the simulation has likely finished. Once this is done, you can execute the
included plotter.py function in a directory containing the node files with the
command line argument of how many images your string had. The script also
accepts an argument to plot a free energy surface alongside the string
(generated with another method), but that
goes beyond the scope of this tutorial. Thus, simply execute:

.. code-block:: bash

	python plotter.py 22 none

And in a moment you should have a graph of your converged string. Thus concludes
this tutorial.

Developer
^^^^^^^^^

* Cody Bezik
Introduction
============

Welcome to SSAGES, our extensive advanced sampling package. You might be
wondering---what is SSAGES and what can it do for my research?

Over the past several decades, molecular simulation has emerged as a powerful
tool for investigating a wide range of physical phenomena. Molecular
simulation is, in essence, a computational "microscope" whereby computers are
used to "look at" the properties of a system that are difficult to observe or
measure through traditional experimental setups. The comparison between
simulations and the corresponding experimental systems can sometimes be
challenging, usually due to factors such as the length and time scales
explored. In simulation, a molecular model must have sufficient temporal and
spatial accuracy to resolve the fastest time scales and shortest length scales
within a system. Unfortunately, due to computational constraints, this
detailed resolution has limited the length of time and number of particles
that a model can simulate, typically simulating systems that are smaller than
analogous experimental setups in laboratory environments for much shorter
times than the duration of the experiments. Recent advancements in processing
power, including custom-built computer architectures and GPU-based computing,
have continued to increase the time and length scales accessible by molecular
simulation, with current state-of-the-art simulations able to analyze systems
for milliseconds (:math:`10^{-3}` s) or more :cite:`PIANA2013201218321`.

However, other challenges arise in obtaining good statistics from molecular
simulations.  Thermal fluctuations dominate motion at the nano-scale and
result in motion that appears random (i.e. Brownian), with no two molecular
trajectories being identical. As a result, statistically meaningful averages
are necessary in order to calculate thermodynamic and kinetic quantities of
interest in these systems :cite:`FRENKELANDSMIT`. An incredibly powerful
thermodynamic quantity referred to as the relative free energy of a system can
be calculated in this way. The relative free energy can characterize
underlying system behavior in the presence of thermally-induced random
noise. Performing this necessary averaging within simulations is challenging.
In essence, the requirement of averaging compounds the issue of time scales
described previously; not only must long simulations be performed, but they
must be performed a prohibitively large number of times in order to extract
sufficient statistics. It is therefore necessary to develop efficient
techniques to calculate meaningful averages from simulations.

Advanced sampling methods represent a class of simulation techniques that seek
to improve this improper averaging and accelerate the extraction of useful
properties (e.g. free energies, transition paths) from simulations.  At the
heart of all advanced sampling methods is statistical mechanics, a field of
physics that relates microscopic phenomena (e.g. the motion of particles) to
macroscopic observables (e.g. temperature and pressure). By taking advantage of
statistical mechanics, advanced sampling methods are used to apply a systematic
bias to a simulation to speed convergence, and then mathematically remove this
bias to extract the true underlying behavior. Throughout the past decade,
advanced sampling methods have become wildly successful, and have now become an
essential component in the toolbox of molecular simulation.

Despite the demonstrated utility of advanced sampling techniques, they have only
been adopted by a fraction of the scientists working in the field. One
explanation for this slow adoption is technical: advanced sampling methods are
complicated, and not all research groups have the expertise required in order to
implement these methods themselves. In the worst case, this leads to long stages
of code development, possibly leading to unknown implementation errors or
insufficient validation. Even in cases when advanced sampling methods are
implemented, they are typically done so with a specific problem in mind and are
custom-built for a certain model or application. This specificity necessitates
modification of the custom-built advanced sampling code when studying new
systems. This prevents the distribution of code between researches in the field.
As a result, the same methods are implemented again and again by different
members of the community. Sadly, in molecular simulation, it is quite common to
"reinvent the wheel".

SSAGES is an answer to this problem :cite:`SSAGES`. SSAGES (Software Suite for
Advanced General Ensemble Simulations) is a free, open-source software package
that allows users to easily apply advanced sampling techniques to any
molecular system of interest. Simply put, SSAGES is a wrapper that converts a
molecular simulation engine (e.g. LAMMPS, GROMACS) into an advanced sampling
engine. SSAGES contains a library of widely used advanced sampling methods
that can be used to calculate everything from free energies to transition
pathways. Importantly, SSAGES works with many widely used simulation
packages, and can simply be added on top of the simulations a researcher is
already running. SSAGES is implemented in a highly modular way, and is easily
extended to incorporate a new method or to modify an existing one and has been
rigorously tested to ensure the accuracy of its calculations.

In short, SSAGES makes advanced sampling methods easy. We hope that it will do
just that for your research.
:orphan:

.. _replica-exchange:

Replica Exchange
----------------

Introduction
^^^^^^^^^^^^

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Tutorial
^^^^^^^^

Developer
^^^^^^^^^

.. _getting-started:

Getting Started
===============

Prerequisites
-------------

Before you try to build SSAGES, make sure that you have the following packages
installed:

+------------+------------------+-----------------------------------+
| Package    | Required version | Package name in Ubuntu repository |
+============+==================+===================================+
| `openmpi`_ | 1.8 or higher    | openmpi-common, libopenmpi-dev    |
+------------+------------------+-----------------------------------+
| `gcc`_     | 4.9 or higher    | gcc-4.9                           |
+------------+------------------+-----------------------------------+
| `cmake`_   | 2.8 or higher    | cmake                             |
+------------+------------------+-----------------------------------+
| `python`_  | 2.7              | python2.7                         |
+------------+------------------+-----------------------------------+

.. _openmpi: https://www.open-mpi.org/
.. _gcc: https://gcc.gnu.org/
.. _cmake: https://cmake.org/
.. _python: https://www.python.org/

Get the Source Code
-------------------

There are two ways of getting the source code for SSAGES: Download a ZIP file
from GitHub or clone the Git repository. We strongly recommend the second
method, as it allows you to easily stay up-to-date with the latest version.

To clone the Git repository, call

.. code-block:: bash

    git clone https://github.com/SSAGESproject/SSAGES.git

Build SSAGES
------------

SSAGES supports a number of simulation engines. If you have already downloaded
an engine's source files, the general steps for building SSAGES with this
engine are as follows. (The example here assumes the engine is called "Engine"
and uses the respective CMake flag ``-DENGINE_SRC``.)

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DENGINE_SRC=/path/to/engine
    make

The different supported engines have respective ``_SRC`` flags to specify
their source code path. For engine-specific steps and options, consult the
:ref:`Engines <engines>` page, which provides more in-depth directions. SSAGES
supports most CMake flags, including ``-DCMAKE_BUILD_TYPE`` or specifying the
compilers with ``-DCMAKE_C_COMPILER=`` and ``-DCMAKE_CXX_COMPILER=``.

It is possible to have SSAGES auto-download the source codes for LAMMPS and
GROMACS. This is done by providing the option ``-DLAMMPS=YES`` for LAMMPS and
``-DGROMACS=YES`` for GROMACS to CMake. However, in many cases, it will be
necessary to build SSAGES using your local copy of the MD engine source code.
For example, if you have modified it to fit a special need, LAMMPS or GROMACS
does not support natively.

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DLAMMPS=YES
    make

or

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DGROMACS=YES
    make

This set of commands will automatically download LAMMPS/GROMACS and build it together
with SSAGES.

For these and other supported engines, you can build with a local copy of the
MD engine source code, as described in the :ref:`Engines <engines>` section.
This can be useful if you are building on a system that has limited or no
connectivity or if you want to add user modifications to the engine's
source code before compiling.

Run SSAGES
----------

In order to run SSAGES, you call the executable followed by the input file.
For example, with an input file called input.json, simple single-core jobs
can call

.. code-block:: bash

    ./ssages input.json

while jobs running on multiple threads can call

.. code-block:: bash

    mpiexec -np 6 ./ssages input.json

Here, the ``-np`` flag dictates the total number of processors on which the
simulation will run and input.json is the input file. For more information,
consult the :ref:`Input Files <inputfiles>` section.

Advanced Options
----------------

In case these simple steps do not meet your need, you can find engine-specific
options and advanced information on building and running SSAGES in the
:ref:`Engines <engines>` section.
:orphan:

.. image_method:

Image Method
------------

Introduction
^^^^^^^^^^^^

Surface charging or polarization can strongly affect the nature of interactions
between charged dielectric objects, particularly when sharp dielectric
discontinuities are involved. However, not any efficient and accurate
computation tools are publicly available especially for the description of
polarization effects in many-body systems. 

For this purpose, Image Method, an analytic perturbative approach we recently
developed for evaluating the polarization energy of a many-body collection of
charged dielectric spheres embedded in a dielectric medium becomes particularly
suitable [1]_.

The polarization-induced interactions between these spheres depend on the ratio
of dielectric constants for the spheres and the medium, and the ratio of the
distance between particles and the radii of the particles. We have shown that,
in some cases, polarization completely alters the qualitative behavior, and in
some other cases, polarization leads to stable configurations that otherwise
could not occur in its absence. 

We think it is helpful to include Image Method into SSAGES for users to include
polarization corrections properly in their systems, and meanwhile, to couple
with advanced sampling methods to accelerate their simulations. 

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

SSAGES Image method is implemented in a way that is as easy as conducting a
simulation using LAMMPS that only includes pairwise Coulombic interactions into
electrostatic interactions. To achieve this, we update the electrostatic forces
acting on all objects by adding up the polarization corrections using SSAGES
engine and then pass the modified snapshot back to LAMMPS engine at each time
step. The JSON file needed for SSAGES engine should include:

einner
    The relative dielectric permittivity of polarizable object. 

ion-type-start
    For cases that you have both polarizable objects and non-polarizable objects
    in you system, for example, in which colloids and ions are treated as
    polarizable and non-polarizable, respectively. This parameter controls where
    the non-polarizable typos start. 

atom type radius
    Radius of all types of objects. 

Guidelines
^^^^^^^^^^

It is very similar as running a simulation including electrostatic interactions
using LAMMPS. Referring to the exampled LAMMPS INPUTFILE and DATAFILE, you need
to double check you have declared the following variables that are particularly
necessary for Image Method to compute polarization corrections: 

* charges
* dielectric (relative dielectric permittivity of the surrounding continuum)

Method Output
^^^^^^^^^^^^^

There are not special outputs files generated for Image method since it only
provides an updated electrostatic forces by including polarization corrections.
Nevertheless, we provided options of dumping trajectories and printing out
force-distance data in the LAMMPS INPUTFILE examples for users to visualize how
significant the polarization effects are in some cases more conveniently. 

.. _IM_tutorial:

Tutorial
^^^^^^^^

.. todo::

    Write a tutorial. 

Developer
^^^^^^^^^

* Jiyuan Li

References
^^^^^^^^^^

.. [1] J. Qin, J. Li, V. Lee, H. Jaeger, J. J. de Pablo, and K. Freed,
       *A theory of interactions between polarizable dielectric spheres*,
       J. Coll. Int. Sci. **469**, 237 - 241 (2016)
.. SSAGES documentation master file, created by
   sphinx-quickstart on Mon Jun  6 16:05:35 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SSAGES!
==================

In their simplest form, particle-based simulations are limited to generating
ensembles of configurations (in Monte Carlo [MC] simulations) or trajectories
in time (in molecular dynamics [MD] or Brownian dynamics [BD]). One can then
extract mechanical variables such as the potential energy or pressure and
perform ensemble or time averages. There are two important limitations to such
calculations:

1. For complex materials, the time scales available to standard
   MD simulations are often insufficient to sample relevant regions of phase
   space
2. In order to develop a fundamental understanding of materials,
   researchers are primarily interested in calculating the free energy, the
   entropy, and their derivatives with respect to various thermodynamic
   quantities (which lead to material properties such as elastic moduli, heat
   capacity, and various other susceptibilities).

These quantities are difficult to obtain or intractable in standard MC and MD
simulations. To overcome these limitations, MC and MD simulations must be
supplemented with advanced sampling techniques. These methods are critical for
the efficient simulation of complex assembly processes.

SSAGES (Software Suite for Advanced General Ensemble Simulations) is
designed to perform these calculations. The framework is designed to treat
molecular simulation routines as a black box, using the coordinates of the
system as evolved by an MD engine to compute collective variables which
permit a meaningful reduced-dimensionality representation of the phase space
within a system. This information is then used to define evolving reactive
pathways or to bias the statistics of a simulation for the purposes of
computing free energies. The internal structure of the code has been designed
to be simple and extensible to new sampling methods and engines. For further
details on examples and capabilities of SSAGES, peruse the documentation
for :ref:`specific methods <advanced-sampling-methods>`.

Contents:

.. toctree::
   :maxdepth: 2

   Introduction
   Getting Started
   Input Files
   Engines
   Collective Variables
   Advanced Sampling Methods
   Write your own Methods and CVs
   Contribute to SSAGES
   The SSAGES Cookbook
   Acknowledgements
   zReferences
   Copyright and License

.. Developer's note: References is named "zReferences.rst" to appear at the
   end of the file list, as sphinxcontrib-bibtex must have the references
   file compiled *last*.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _Finite-temperature-string:

Finite Temperature String
-------------------------

Introduction
^^^^^^^^^^^^

Along with Nudged Elastic Band and Swarm of Trajectories, Finite Temperature
String Method (FTS) is a chain-of-states method. As in other chain-of-states
methods, multiple copies of a system are simulated, with each copy ("image")
corresponding to a different state of the system along some proposed transition
pathway. In FTS, each image is associated with a node along a smooth curve
through in collective variable space, representing a transition pathway.

The goal of FTS is to evolve the path of this smooth curve or "string" until it
approximates a transition pathway by finding the **principal curve**, which by
definition intersects each of the perpendicular hyperplanes that it passes
through at the expected value of each hyperplane. As such, a principal curve is
often referred to as being its own expectation. Rather than sampling along each
hyperplane belonging to each node along the string, we use the Voronoi
approximation introduced by Vanden-Eijnden and Venturoli in 2009
:cite:`VANDENEIJNDEN200905B605`. We
associate each node along the string with a corresponding Voronoi cell,
consisting of the region in state space where any point is closer to its origin
node than any other node along the string. Each image is free to explore within
the bounds of its associated Voronoi cell. To evolve the string toward its own
expectation, the string is evolved toward the running averages in CV space for
each image along the string.

The evolution of the string can be broken down into the following steps:

1. Evolve the individual images with some dynamics scheme, using the location of
   the initial image as a starting point. Only keep the new update at each time
   step if it falls within the Voronoi cell of its associated image; if the
   updated position leaves the Voronoi cell, the system is returned back to the
   state at the previous timestep. 

2. Keep track of a running average of locations visited in CV space for each
   image.

3. Update each node on the string toward the running average while keeping the
   path smooth; specific equations can be found in
   :cite:`VANDENEIJNDEN200905B605`.

4. Enforce parametrization (e.g. interpolate a smooth curve through the new node
   locations, and redistribute the nodes to new locations along the smooth curve
   such that there is equal arc length between any two adjacent nodes).

5. After images have been moved, their respective Voronoi cells have also
   changed. Check that each image still falls within the new Voronoi cell of its
   associated image. If the image is no longer in the correct Voronoi cell, the
   system must be returned to the Voronoi cell.

6. Return to step 1 and repeat until convergence (i.e. until change in the string
   falls below some tolerance criteria or stop iterating after a certain number
   of string method iterations)

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

The following parameters need to be set under ``"method"`` in the JSON input file:

.. code-block:: javascript

	"type" : "String"
	"flavor" : "FTS"

The following options are available as FTS inputs: 

centers (required)
	Array containing this image's coordinates in CV space.

ksprings (required)
	Array of spring constants corresponding to each CV.
	Used to ensure that each simulation remains within its own respective
	Voronoi cell.

block_iterations (required)
	Number of integration steps to perform before updating the string.
	(Default: 2000)

time_step (required)
	Parameter used for updating the string.
	(:math:`\Delta\tau` in :cite:`VANDENEIJNDEN200905B605`).
	(Default: 0.1)

kappa (required)
	Parameter used for smoothing the string.
	(:math:`\kappa` in :cite:`VANDENEIJNDEN200905B605`).
	(Default: 0.1)

frequency (required)
	Frequency to perform integration; should almost always be set to 1.

max_iterations (required)
	Maximum number of string method iterations to perform.

tolerance (required)
	Array of tolerance values corresponding to each CV. Simulation will stop
	after tolerance criteria has been met for all CVs.

iteration
	Value of initial string method iterator.
	(Default: 0)

.. _FTS_tutorial:

Tutorial
^^^^^^^^

Two examples for running FTS can be found in the ``Examples/User/FTS``
directory. This tutorial will go through running FTS on a 2D single particle
system, using LAMMPS as the MD engine. The necessary files are found in
``Examples/User/FTS/2D_Particle``, which should contain the following:

``in.LAMMPS_2DParticle``
	LAMMPS input file; sets up 1 particle on a 2D surface with two Gaussian
	wells of different depths (at :math:`(-0.98, -0.98)` and
	:math:`(0.98, 0.98)`) and one Gaussian barrier at the origin.

``Template_Input.json``
	Template JSON input containing information for one image on the string. We
	are looking at two CVs: x and y coordinates. We will use
	``Input_Generator.py`` to use this template to create a JSON input file
	containing information for all string images.

``Input_Generator.py``
	Python script for creating FTS JSON input file.

After compiling SSAGES with LAMMPS, we will use ``Input_Generator.py`` to
create a JSON input file for FTS. Run this script

.. code-block:: bash

	python Input_Generator.py

to create a file called ``FTS.json``. A string with 16 images is initalized on
the 2D surface, evenly spaced on a straight line from :math:`(-0.98, -0.68)` to
:math:`(0.98, 1.28)`. If you take a look at ``FTS.json``, you will see that the
location of each image along the string has been appended to the ``"centers"``
field. These center locations are listed from one end of the string to the
other; the first center listed corresponds to one end of the string, and the 
final center listed corresponds to the opposite end of the string.

Once ``FTS.json`` has been generated, we can run the example with the following
command: 

.. code-block:: bash

	mpiexec -np 16 ./ssages FTS.json

As SSAGES runs, a series of output files are generated: 

``log.lammps``
	Output from LAMMPS.

``node-00xx.log``
	FTS output for each of the 16 nodes on the string. The first column contains
	the image number (0-15). The second column contains the iteration number. The
	remaining columns list the location of the image and the instantaneous value
	for each of the CVs. For this example we have two CVs (x coordinate and y
	coordinate), so the remaining columns are (from left to right): x coordinate
	of the string node, instantaneous x coordinate of the particle, y coordinate
	of the string node, instantaneous y coordinate of the particle.

To visualize the string, we can plot the appropriate values from the last line
of each ``node-00xx.log`` file. For example, one can quickly plot the final
string using gnuplot with the command

.. code-block:: gnuplot

	plot "< tail -n 1 node*" u 3:5

The following image shows the initial string in blue, compared with the final
string plotted in green: 

.. figure:: images/2dsingle.png
	:align: center

The two ends of the string have moved to the two energy minima (at
:math:`(-0.98, -0.98)` and :math:`(0.98, 0.98)`), and the center of the string
has curved away from the energy barrier at the origin. 

Developers
^^^^^^^^^^

* Ashley Guo
* Benjamin Sikora
* Yamil J. Colón
.. _Forward-flux:

Forward-Flux
------------

Forward Flux Sampling (FFS) is an advanced sampling method to simulate
"rare events" in non-equilibrium and equilibrium systems. Several review
articles in the literature present a comprehensive perspective on the basics,
applications, implementations, and recent advances of FFS. Here, we provide a
brief general introduction to FFS, and describe the Rosenbluth-like variant of
forward flux method. We also explain various options and variables to setup
and run an efficient FFS simulation using SSAGES.

Introduction
^^^^^^^^^^^^

Rare events are ubiquitous in nature. Important examples include crystal
nucleation, earthquake formation, slow chemical reactions, protein
conformational changes, switching in biochemical networks, and translocation
through pores. The activated/rare process from a stable/metastable region A to
a stable/metastable region B is characterized by a long waiting time between
events, which is several orders of magnitude longer than the transition process
itself. This long waiting time typically arises due to the presence of a large
free energy barrier that the system has to overcome to make the transition from
one region to another. The outcomes of rare events are generally substantial
and thereby it is essential to obtain a molecular-level understanding of the
mechanisms and kinetics of these events.

"Thermal fluctuations" commonly drive the systems from an initial state to a
final state over an energy barrier :math:`\Delta E`. The transition frequency
from state A to state B is proportional to :math:`e^{\frac{-\Delta E}{k_{B}T}}`,
where :math:`k_{B}T` is the thermal energy of the system. Accordingly, the time
required for an equilibrated system in state A to reach state B grows
exponentially (at a constant temperature) as the energy barrier :math:`\Delta E`
become larger. Eventually, none or only a few transitions may occur within the
typical timescale of molecular simulations. In FFS method, several intermediate
states or so-called interfaces (:math:`\lambda_{i}`) are placed along a
"reaction coordinate" or an "order parameter" between the initial state A and
the final state B (Figure 1). These intermediate states are chosen such that the
energy barrier between adjacent interfaces are readily surmountable using
typical simulations. Using the stored configurations at an interface, several
attempts are made to arrive at the next interface in the forward direction (the
order parameter must change monotonically when going from A to B). This
incremental progress makes it more probable to observe a full transition path
from state A to state B. FFS uses positive flux expression to calculate rate
constant. The system dynamics are integrated forward in time and therefore
detailed balance is not required.

.. figure:: images/forward_flux_image1.png
	:align: center
	:alt: Forward Flux Sampling setup

	In the Forward Flux Sampling method, several intermediate states are placed
	along the order parameter to link the initial state A and the final
	state B. Incremental progress of the system is recorded and analyzed to
	obtain relevant kinetic and thermodynamic properties.

Several protocols of Forward Flux Sampling have been adopted in the literature to

1. generate the intermediate configurations,
2. calculate the conditional probability of reaching state B starting from
   state A, :math:`P(\lambda_{B} = \lambda_{n} | \lambda_{A} = \lambda_{0})`,
3. compute various thermodynamic properties, and
4. optimize overall efficiency of the method :cite:`ALLEN2006194111`.
   The following are the widely-used variants of Forward Flux Sampling method:

* Direct FFS (DFFS) (currently implemented in SSAGES)
* Branched Growth FFS (BGFFS)
* Rosenbluth-like FFS (RBFFS)
* Restricted Branched Growth FFS (RBGFFS)
* FFS Least-Squares Estimation (FFS-LSE)
* FF Umbrella Sampling (FF-US)

Rate Constant and Initial Flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall rate constant or the frequency of going from state A to state B is
computed using the following equation:

.. math::

	k_{AB} = \Phi_{A,0} \cdot P\left(\lambda_{N} \vert \lambda_{0}\right)

here, :math:`\Phi_{A,0}` is the initial forward flux or the flux at the initial
interface, and :math:`P\left(\lambda_{N} \vert \lambda_{0}\right)` is the
conditional probability of the trajectories that initiated from A and reached B
before returning to A. In practice, :math:`\Phi_{A,0}` can be obtained by
simulating a single trajectory in State A for a certain amount of time
:math:`t_{A}`, and counting the number of crossings of the initial interface
:math:`\lambda_{0}`. Alternatively, a simulation may be carried out around state
A for a period of time until :math:`N_{0}` number of accumulated
configurations is stored:

.. math::

	\Phi_{A,0} = \frac{N_{0}}{t_{A}}

here, :math:`N_{0}` is the number of instances in which :math:`\lambda_{0}` is
crossed in forward direction, and :math:`t_{A}` is the simulation time that the system was run around
state A. Note that

1. :math:`\lambda_{0}` can be crossed in either forward
   (:math:`\lambda_{t} < \lambda_{0}`) or backward
   (:math:`\lambda_{t} > \lambda_{0}`) directions, but only "forward crossing"
   marks a checkpoint (see Figure 2) and
2. :math:`t_{A}` should only include the simulation time around state A and
   thereby the portion of time spent around state B must be excluded, if any.

In general, the conditional probability is computed using the following expression:

.. math::

	P\left(\lambda_{n} \vert \lambda_{0}\right) =
	\prod\limits_{i=0}^{n-1} P\left(\lambda_{i+1} \vert \lambda_{i}\right) =
	P\left(\lambda_{1}\vert\lambda_{0}\right) \cdot
	P\left(\lambda_{2}\vert\lambda_{1}\right) \dots
	P\left(\lambda_{n}\vert\lambda_{n-1}\right)

:math:`P\left(\lambda_{i+1}\vert\lambda_{i}\right)` is computed by initiating a
large number of trials from the current interface and recording the number of
successful trials that reaches the next interface. The successful trials in
which the system reaches the next interface are stored and used as
checkpoints in the next interface. The failed trajectories that go all the way
back to state A are terminated. Different flavors of forward flux method use
their unique protocol to select checkpoints to initiate trials at a given
interface, compute final probabilities, create transitions paths, and analyze
additional statistics.

.. figure:: images/forward_flux_image2.png
	:align: center
	:alt: Forward Flux Sampling initial flux

	A schematic representation of computation of initial flux using a single
	trajectory initiated in state A. The simulation runs for a certain period of
	time :math:`t_{A}` and number of forward crossing is recorded. Alternatively,
	we can specify the number of necessary checkpoints :math:`N_{0}` and run a
	simulation until desired number of checkpoints are collected. In this figure,
	green circles show the checkpoints that can be used to generate transition
	paths.

Rosenbluth-like Forward Flux Sampling (RBFFS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rosenbluth-like Forward Flux Sampling (RBFFS) method is an adaptation of
Rosenbluth method in polymer sampling to simulate rare events
:cite:`ROSENBLUTH1955356`.
The RBFFS is comparable to Branched Growth Forward Flux (BGFFS)
:cite:`ALLEN2009463102,ESCOBEDO2009333101` but,
in contrast to BGFFS, a single checkpoint is randomly selected at a non-initial
interface instead of initiation of trials from all checkpoints at a given
interface (Figure 3). In RBFFS, first a checkpoint at :math:`\lambda_{0}` is
selected and :math:`k_{0}` trials are initiated. The successful runs that reach
:math:`\lambda_{1}` are stored. Next, one of the checkpoints at :math:`\lambda_{1}` is randomly chosen (in
contrast to Branched Growth where all checkpoints are involved), and
:math:`k_{1}` trials are initiated to :math:`\lambda_{2}`. Last, this procedure
is continued for the following interfaces until state B is reached or all trials
fail. This algorithm is then repeated for the remaining checkpoints at
:math:`\lambda_{0}` to generate multiple "transition paths".

.. figure:: images/forward_flux_image3.png
	:align: center
	:alt: Rosenbluth-like Forward Flux Sampling

	Rosenbluth-like Forward Flux Sampling (RBFFS) involves sequential generation
	of unbranched transition paths from all available checkpoints at the first
	interface :math:`\lambda_{0}`. A single checkpoint at the interface
	:math:`\lambda_{i > 0}`  is randomly marked and :math:`k_{i}` trials are
	initiated from that checkpoint which may reach to the next interface
	:math:`\lambda_{i+1}` (successful trials) or may return to state A (failed
	trial).

In Rosenbluth-like forward flux sampling, we choose one checkpoint from each
interface independent of the number of successes. The number of available
checkpoints at an interface are not necessarily identical for different
transition paths :math:`p`. This implies that more successful transition paths
are artificially more depleted than less successful paths. Therefore, we need to
enhance those extra-depleted paths by reweighting them during post-processing.
The weight of path :math:`p` at the interface :math:`\lambda_{i}` is given by:

.. math::

	w_{i,b} = \prod\limits_{j=0}^{i-1} \frac{S_{j,p}}{k_{j}}

where :math:`S_{j,p}` is the number of successes at the interface :math:`j` for
path :math:`p`. The conditional probability is then computed using the following
expression:

.. math::

	P\left(\lambda_{n}\vert\lambda_{0}\right) =
	\prod\limits_{i=0}^{n-1} P\left(\lambda_{i+1} \vert \lambda_{i}\right) =
	frac{ \prod_{i=0}^{n-1}\sum_{p} w_{i,p} S_{i,p} / k_{i} }{ \sum_{p} w_{i,p} }

Here, the summation runs over all transition paths in the simulation.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

The notation used in SSAGES implementation of the FFS is mainly drawn from
Ref. :cite:`ALLEN2009463102`. We recommend referring to this review article if
the user is unfamiliar with the terminology.

To run a DFFS simulation using SSAGES, an input file in JSON format is required
along with a general input file designed for your choice of molecular dynamics
engine (MD engine). For your convenience, two files ``Template_Input.json`` and
``FF_Input_Generator.py`` are provided to assist you in generating the JSON
file. Here we describe the parameters and the options that should be set in
``Template_Input.json`` file in order to successfully generate an input file and
run a DFFS simulation.

.. warning::

	The current implementation of FFS only accepts one CV.

The following parameters need to be set under ``"method"`` in the JSON input file:

.. code-block:: javascript

	"type": "ForwardFlux"

The following options are available for Forward Flux Sampling:

flavor (required)
	Specifies the flavor of the FFS method that SSAGES should run.
	Available options: "DirectForwardFlux"

.. note::

	Currently, only DFFS has been implemented in SSAGES.
	RBFFS and BGFFS will be available in the future releases.

trials (required)
	Array of number of trials to be spawned from each interface. The length of
	this array should match the length of the array of ``interfaces``, or can
	be left blank (``[]``) if defined in ``FF_Input_Generator.py``.

interfaces (required)
	Array of intermediate interfaces linking the initial state A to the final
	state B. This array can either be defined in ``Template_Input.json``
	or ``FF_Input_Generator.py``. In the latter case, the values of
	**interfaces** is left blank in the ``Template_Input.json`` file.

nInterfaces (optional)
	Total number of interfaces connecting the initial state A to
	the final state B, inclusive. (Default: 5)

.. warning::

	Minimum of two interfaces must be defined.

N0Target (optional)
	Number of configurations to be generated (or provided by user) at the
	first interface. (Default: 100)

computeInitialFlux (optional)
	Specifies whether a calculation of the initial flux should be performed.
	If this parameter is set to ``true``, SSAGES would also generate the
	user-specified number of initial configurations (``N0Target``) at the first
	interface. To compute the initial flux, user must provide an initial
	configuration in state A, otherwise SSAGES would issue an error. If this
	parameter is set to ``false``, the user must provide the necessary number
	of the initial configurations in separate files. The files name and the
	files content should follow a specific format. The format of the filenames
	should be ``l0-n<n>.dat`` where <n> is the configuration number
	(i.e. ``1``, ``2``, ..., ``N0Target``).
	
	The first line of the configuration files includes three numbers
	``<l> <n> <a>``, where ``<l>`` is the interface number (set to zero here),
	``<n>`` is the configuration number, and ``<a>`` is the attempt
	number (set to zero here). The rest of the lines include the atoms IDs and
	their corresponding values of positions and velocities, in the format
	``<atom ID> <x> <y> <z> <vx> <vy> <vz>``
	where <atom ID> is the ID of an atoms,
	<x>, <y>, <z> are the coordinates of that atom,
	and <vx>, <vy>, and <vz> are the components of the velocity in
	the x, y, and z directions. Please note that the stored configurations at
	other interfaces follow a similar format. (Default: ``true``)

saveTrajectories (optional)
	This flag determines if the FFS trajectories should be saved.
	(Default: ``true``)

.. warning::

	Saving trajectories of thousands of atoms may require large amount
	of storage space.

currentInterface (optional)
	Specifies the interface from which the calculations should start
	(or continue). This parameter is helpful in restarting a FFS calculation
	from interfaces other than the initial state A. (Default: 0)

outputDirectoryName (optional)
	Specifies the directory name that contains the output of the FFS
	calculations including the initial flux, the successful and failed
	configurations, commitor probabilities, and the trajectories. The output
	data related to the computation of the initial flux is stored in the file
	``initial_flux_value.dat``, and the data related to transition probabilities
	is stored in the file ``commitor_probabilities.dat``. (Default: "FFSoutput")

.. _FFS_tutorial:

Tutorial
^^^^^^^^

This tutorial will walk you step-by-step through the user example provided with
the SSAGES source code that runs the forward flux method on a Langevin particle
in a two-dimensional potential energy surface using LAMMPS. This example shows
how to prepare a multi-walker simulation (here we use 2 walkers). First, be
sure you have compiled SSAGES with LAMMPS. Then, navigate to the
``Examples/User/ForwardFlux/LAMMPS/Langevin`` subdirectory. Now, take a moment
to observe the ``in.LAMMPS_FF_Test_1d`` file to familiarize
yourself with the system being simulated.

The next two files of interest are the ``Template_Input.json`` input file and
the ``FF_Input_Generator.py`` script. These files are provided to help setup
sophisticated simulations. Both of these files can be modified in your
text editor of choice to customize your input files, but for this tutorial,
simply observe them and leave them be. ``FF_Template.json`` contains all
information necessary to fully specify a walker; ``FF_Input_Generator.py``
uses the information in this file and generates a new JSON along with necessary
LAMMPS input files. Issue the following command to generate the files:

.. code-block:: bash

	python FF_Input_Generator.py

You will produce a file called ``Input-2walkers.json`` along with
``in.LAMMPS_FF_Test_1d-0`` and ``in.LAMMPS_FF_Test_1d-1``. You can also open
these files to verify for yourself that the script did what it was supposed to
do. Now, with your JSON input and your SSAGES binary, you have everything you
need to perform a simulation. Simply run:

.. code-block:: bash

	mpiexec -np 2 ./ssages Input-2walkers.json

This should run a quick FFS calculation and generate the necessary output.

Developers
^^^^^^^^^^

* Joshua Lequieu
* Hadi Ramezani-Dakhel
* Benjamin Sikora
* Vikram Thapar
.. _engines:

Engines
=======

SSAGES supports multiple molecular dynamics engines. However, supported
features between MD engines. This may be the result of engine limitations
or work in progress. The table below summarizes the main features that
vary between supported engines.

+---------------+-----------------------+--------------+------------+
| Engine        | Supported versions    | Multi-walker | NPT virial |
+===============+=======================+==============+============+
| `LAMMPS`_     | 2010 or newer         | yes          | yes        |
+---------------+-----------------------+--------------+------------+
| `GROMACS`_    | 5.1.x, 2016.x, 2018.x | yes          | yes        |
+---------------+-----------------------+--------------+------------+
| `OpenMD`_     | 2.5, 2.6              | no           | no         |
+---------------+-----------------------+--------------+------------+
| `Qbox`_       | 1.60 or newer         | yes          | no         |
+---------------+-----------------------+--------------+------------+
| `HOOMD-blue`_ | 2.4.0 or newer        | yes          | no         |
+---------------+-----------------------+--------------+------------+

Special instructions on how to use SSAGES with a particular engine are
listed under the appropriate section.

LAMMPS
^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of LAMMPS. To compile SSAGES with a
compatible version of LAMMPS, either ``-DLAMMPS=YES`` or
``-DLAMMPS_SRC=/path/to/LAMMPS`` must be specified in the CMake command.
For example,

.. code:: bash

	cmake -DLAMMPS=YES ..
	make

will automatically download LAMMPS (version 30 Jul 2016, tagged ``r15407``)
and compile SSAGES. If a user is interested in using a different version of
LAMMPS, SSAGES can download a specific stable release (must be supported) with

.. code:: bash

	cmake -DLAMMPS="22 Aug 2018" ..

If a user is interested in using an already-downloaded source or one with
personal modifications, then SSAGES can be pointed to that particular source
repository.

.. code:: bash

	cmake -DLAMMPS_SRC=/path/to/lammps ..

Because many users may take advantage of optional LAMMPS
packages, SSAGES forwards the make commands necessary to do so. To enable or
disable these packages, you can call

.. code-block:: bash

	make yes-package

or

.. code-block:: bash

	make no-package

For more information on optional packages for LAMMPS, refer to
the `LAMMPS User Manual <https://lammps.sandia.gov/doc/Packages.html>`_.

.. warning::

	Once you link SSAGES to a particular LAMMPS source, you will be
	**unable** to compile that LAMMPS source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following stable versions of LAMMPS have been tested extensively, but we
are confident that SSAGES will also work with most other LAMMPS versions.

* 10 Aug 2015
* 7 Dec 2015
* 16 Feb 2016
* 14 May 2016
* 30 Jul 2016
* 5 Nov 2016
* 17 Nov 2016
* 31 Mar 2017
* 11 Aug 2017
* 16 Mar 2018
* 22 Aug 2018

Running
~~~~~~~

SSAGES integrates with LAMMPS though the flexible *fix* API offered
by LAMMPS. It is therefore necessary to define a SSAGES fix within
the LAMMPS input file as follows.

.. code-block:: none

	fix ssages all ssages

This directive ensures that SSAGES is able to locate the appropriate
adapter and interface with the LAMMPS library. **It is very important to
name the fix "ssages" as shown above. Otherwise, SSAGES will not work
properly**. It is highly recommended that the SSAGES fix command be placed
after all integrator fixes. Also, make sure that the fix is specified before
the run command, which will begin the advanced sampling simulation.

.. note::

	Due to the nature of how SSAGES forwards commands to LAMMPS, the use
	of ``include`` and ``label/jump`` within a LAMMPS input script is
	currently not supported.

SSAGES is compatible with typical LAMMPS workflows that include equilibration
or energy minimization steps before production. So long as the SSAGES fix is not
declared, LAMMPS will run without any modification.

The only LAMMPS-specific property required in a SSAGES input file is the ``input``
property which points to the LAMMPS input script. Details can be found on the
:ref:`input files page <inputfiles>`.

GROMACS
^^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of GROMACS. To compile SSAGES with a
compatible version of GROMACS, either ``-DGROMACS=YES`` or
``-DGROMACS_SRC=/path/to/GROMACS`` must be specified in the CMake command.
For example,

.. code:: bash

	cmake -DGROMACS=YES ..
	make

will automatically download GROMACS 5.1.3 and compile SSAGES.
If a user is interested in using a different version of GROMACS, SSAGES can
download a specific release (must be supported) with

.. code:: bash

	cmake -DGROMACS=2018.3 ..

If a user is interested in using an already-downloaded source or one with
personal modifications, then SSAGES can be pointed to that particular source
repository.

.. code:: bash

	cmake -DGROMACS_SRC=/path/to/gromacs ..

Common options for building GROMACS will be passed through to the GROMACS
compilation step. For instance, options such as ``-DGMX_BUILD_OWN_FFTW=ON``,
``-DGMX_GPU=ON``, and ``-DGMX_DOUBLE=ON`` are supported. With newer versions
of the ``hwloc`` package, GROMACS may not be able to automatically detect the
SIMD level. A workaround is to manually specify this with
``-DGMX_SIMD=AVX2_256``, which will turn off automatic detection using the
``hwloc`` package.

.. warning::

	Once you link SSAGES to a particular GROMACS source, you will be
	**unable** to compile that GROMACS source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following versions of GROMACS are supported by SSAGES, but we are very
confident that SSAGES will *not* work with other versions of GROMACS out of
the box, in contrast to LAMMPS. We are working hard to make SSAGES
compatible with new versions of GROMACS, as they are released.

* 5.1.x
* 2016.x
* 2018.x

Setup
~~~~~

After compiling GROMACS with SSAGES, you can use all of GROMACS’s available
tools to set up systems and generate input files. The executable is located at
``hooks/gromacs/gromacs/bin/gmx_ssages`` within the build directory.

.. note::

	Note that, the ``gmx_ssages`` executable in the SSAGES folder will NOT
	function normally for running regular GROMACS simulations via
	``gmx_ssages mdrun``.

As GROMACS has in-depth
`Documentation <http://manual.gromacs.org/documentation/current/user-guide/>`_
and a helpful
`Getting Started <http://manual.gromacs.org/documentation/current/user-guide/getting-started.html>`_
section, we will not dwell much on how to use these tools to generate systems.

Briefly, generating a GROMACS binary input file (``.tpr``) requires the
following three files:

1. A "box" of particle coordinates to simulate (``.gro`` file)
2. A topology that describes the forcefield and connectivity (``.top`` file,
   optionally ``.itp`` files)
3. A simulation detail file that sets parameters such as which thermostat
   and barostat to use, number and length of time steps, integrator, saving
   frequency, and many more (``.mdp`` file)

For example, one can convert a protein ``.pdb`` file from an
`online database <https://www.rcsb.org/>`_ using GROMACS tools to generate
a ``.gro`` and a ``.top`` file. To generate an input file, use the preprocessor
``gmx_ssages grompp`` command:

.. code-block:: bash

	gmx_ssages grompp -f npt.mdp -p topol.top -c conf.gro -o input.tpr

There are example ``.gro``, ``.mdp``, ``.top``, ``.tpr``, and ``.json`` inputs
available in the Examples folder.

After an energy minimization and brief NVT and NPT equilibration runs, you
should be ready to use SSAGES with your system. First, generate a ``.json``
file for your SSAGES input. If using a single walker, the ``inputfile`` should
be the same as your ``.tpr`` file name. If using multiple walkers, you should
number your input files right before the extension, include a numberless file,
and set the “inputfile” to be the same as the numberless. For example, if using
four walkers, you should set your “inputfile” to ``input.tpr`` and have the
following in your folder:

* ``input.tpr``
* ``input0.tpr``
* ``input1.tpr``
* ``input2.tpr``
* ``input3.tpr``

Finally, define your CV(s) and Methods, as detailed in the
:ref:`input files <inputfiles>` page.

Running
~~~~~~~

SSAGES forwards arguments to the GROMACS **mdrun** library. The
``args`` property must specified in the SSAGES input file as
described on the :ref:`input files <inputfiles>` page.

You can start your simulation by calling the SSAGES executable:

.. code-block:: bash

	mpiexec -np N ./ssages input.json

where `N` is the total number of MPI processes. For example, for three walkers
using 2 processors each, set :math:`N = 3*2 = 6`.

OpenMD
^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of OpenMD. To compile SSAGES with a
compatible version of OpenMD, the location of the already-downloaded source
must be specified in the CMake command.

.. code:: bash

	cmake -DOPENMD_SRC=/path/to/OpenMD ..

.. warning::

	Once you link SSAGES to a particular OpenMD source, you will be
	**unable** to compile that OpenMD source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following versions of OpenMD are supported by SSAGES, but we are very
confident that SSAGES will *not* work with other versions of OpenMD out of
the box, in contrast to LAMMPS.

* 2.5
* 2.6

Running
~~~~~~~

The only OpenMD-specific property required in a SSAGES input file is the ``input``
property which points to the OpenMD input script. Details can be found on the
:ref:`input files page <inputfiles>`.

Qbox
^^^^

Building
~~~~~~~~

SSAGES and Qbox can be run together to use advanced sampling methods in
*ab initio* molecular dynamics simulations. The coupling with Qbox is performed
in a server--driver mode, with SSAGES acting as the driver and Qbox as the
server. This means that if you have access to a version of Qbox (minimum 1.60)
you do not need to recompile SSAGES and Qbox together. However, it is necessary
to configure SSAGES to be used with Qbox, so that it will compile the correct
Hook and Driver. To do so, add the following flag during the configuration of
SSAGES:

.. code:: bash

	cmake -DQBOX=YES ..

It is important to remark that in this case, **SSAGES will not automatically
download Qbox**, it will be simply configured so to communicate with it. You
are required to have access to a Qbox executable. If you do not have access to
a pre-compiled version, then you will need to
`download and compile it yourself <http://qboxcode.org/build/>`_.

Setup
~~~~~
As for other engines, there are two input scripts necessary to run a
Qbox--SSAGES calculation composed of ``N`` walkers:

1. A JSON input file, specifying the methods and CVs that you want to use.
   Also, it specifies the Qbox input file names and the number of MD, density,
   and wavefunction steps that you want to use.
2. A number ``N`` of Qbox input files, that will be used in the first step of
   the calculation to obtain the ground state density in the first step.

The JSON file contains the same field that would usually have (CVs, methods, logger, etc.) with three additional options:

.. code:: javascript

	{
		"walkers": N,
		"input": ["md.1", "md.2", ..., "md.N"],
		"md_iterations" : 10,
		"qm_iterations" : 30,
		"wf_iterations" : 1,
	}

The keywords ``walkers`` and ``input`` are the standard SSAGES keywords to
declare the number of walkers and the starting input file of each walker. The
keywords ``md_iterations``, ``qm_iterations`` and ``wf_iterations``  are the
respectively the number of MD steps to perform, the number of `scf` to perform
per MD step, and the number of wave-function optimization per `scf` steps.
These parameters correspond to the first, second and third number in the
command ``run 20 10 0``. Consult the
`Qbox Documentation <http://qboxcode.org/doc/html/>`_ for more information.

The Qbox input file of each walker, specifies the parameters to be used in the
DFT calculations (``xc``, ``ecut``, ``T``, etc.). This file will be parsed by
Qbox **at the first time step of the simulations** to set up the calculations.
If the file contains a command such as ``run 200 10``, the 200 MD steps that
Qbox will perform **will be unbiased**. If wanted, this feature can be used to
equilibrate the system. After this first step, the command
``run 1 qm_iterations wf_iterations`` will be repeated for ``md_iterations``.
**The QBox input file must contain at least 1 MD step in order to run with SSAGES.**
Thus always include the ``run 1`` command.

An example of ``input.json`` and ``md.i`` is present in the ``Examples/User/ABF/NaCl-Qbox`` directory.

Running
~~~~~~~

As previously reported, Qbox and SSAGES communicate in a server--driver mode.
To launch Qbox in a server mode is sufficient to use the proper keyword and
specify its input and output file:

.. code:: bash

	mpirun -n X qb -server ssages_in_0 ssages_out_0

for a single walker or

.. code:: bash

	mpirun -n X qb -server ssages_in_0 ssages_out_0
	mpirun -n X qb -server ssages_in_1 ssages_out_1
	....
	mpirun -n X qb -server ssages_in_N ssages_out_N

for multiple walkers. At the moment, the name ``ssages_in_`` and
``ssages_out_`` are **mandatory** and cannot be changed. When launched in this
way, Qbox creates ``N`` files called ``ssages_in_N.lock``, and then wait for
input. When the files ``ssages_in_N.lock`` are deleted from disk, Qbox will
execute the commands contained in the files ``ssages_in_N``, write the result
of the calculation in ``ssages_out_N``, and create N ``ssages_in_N.lock``
files. Without the deletion of the ``.lock`` files, Qbox will not execute any
command and will remain idle.

After Qbox has started the server mode run (so it is idling and the ``.lock``
files are present on disk), we can launch SSAGES to drive the calculations:

.. code::

	mpirun -n N ssages input.json

After SSAGES has started, the two codes will alternate with each other in the
following way:

1. SSAGES will write the script ``md.i`` to file ``ssages_in_i``, which will
   initialize the DFT parameters of the calculations. Then, it will trigger
   Qbox execution by deleting the ``.lock`` files.
2. Qbox will perform the DFT calculation specified in ``ssages_in_i`` and write
   the output in ``ssages_out_i`` and will recreate the ``.lock`` files.
3. SSAGES will read the Qbox output, calculate the CVs and the bias, and write
   to the files ``ssages_in_i``, containing the external forces and the
   position of the atoms, as well as the command
   ``run 1 qm_iterations wf_iterations``. It will then delete the ``.lock``
   file, triggering another MD step calculation in Qbox.
4. Steps 2 and 3 will be repeated for ``md_iterations`` number of time.
5. After the last iteration, SSAGES will write an input file that will instruct
   Qbox to save a ``restart_i.xml`` file that can be used to restart the
   calculations, and terminate the Qbox instance.
6. Qbox and SSAGES will then finish the execution.

Normally, Qbox overwrites the output ``ssages_out_i`` in server mode. To
preserve the trajectory and avoid the loss of data, SSAGES will append the
``ssages_out_i`` file to a ``ssages_out_i_run_j.xml`` file. In the latter, the
``i`` index identifies the walker, while the ``j`` index identifies the number
of runs. (For example, if you restarted two times, you would have
``_run_1.xml``, ``_run_2.xml``, and ``_run_3.xml``.) We suggest using the
``restart_i.xml`` files to avoid discontinuities in the trajectories; when
restarting, create a ``md.i`` file that contains the ``load restart_i.xml``
instructions.

There are useful scripts to analyze and plot Qbox trajectories, which are available in the `Qbox tools webpage <http://qboxcode.org/tools//>`_. To run any of these scripts, first reformat ``ssages_out_i_run_j.xml`` file by running a python script ``Qbox-xml-cleaning.py`` present in ``Tools/`` directory. For example, if the original ssages-qbox output file is ``ssages_out_0_run_0.xml`` then the command line to reformat this xml file is

.. code:: bash

	python3 Qbox-xml-cleaning.py ssages_out_0_run_0.xml ssages_out_0_run_0_cleaned.xml

where first arugment the name of the original xml output file, and the second argument is the reformatted xml file. Now these files can be analyzed using the scripts in `Qbox tools webpage <http://qboxcode.org/tools//>`_. For example, to create xyz trajectory file from the reformatted output  ``ssages_out_0_run_0_cleaned.xml``, run the command

.. code:: bash

	python2 qbox_xyz.py -all ssages_out_0_run_0_cleaned.xml > out_0_run_0.xyz
	
Running on Clusters
~~~~~~~~~~~~~~~~~~~

Most likely, you are going to launch this calculation on a cluster or a
supercomputer, where you will need to prepare a submission scripts and then
launch through a job scheduler. Given the fact that SSAGES needs to start
*after* Qbox, it is better to either separate the scripts that submit the two
different calculations, or use a syntax that ensure that the submission occurs
in the right order. For example, on `Slurm <https://slurm.schedmd.com/>`_, we
can use one script:

.. code:: bash

	srun -n X  -N 1 qb -server ssages_in0 ssages_out0  &
	srun -n X  -N 1 qb -server ssages_in1 ssages_out1  &
	srun -n 2  -N 1 ssages input.json &
	wait

which ensures that the scripts are executed in the right way.

If you want to have different scripts for Qbox and SSAGES:

In the Qbox scripts, ``qb.sh``

.. code:: bash

	srun -n X  -N 1 qb -server ssages_in0 ssages_out0
	srun -n X  -N 1 qb -server ssages_in1 ssages_out1

In the SSAGES script, ``ssages.sh``

.. code:: bash

	srun -n 2  -N 1 ssages input.json

Then you will need to submit both of them with a third script, ``launch.sh``

.. code:: bash

	#!/bin/bash

	j_qb=`sbatch qb.sh | awk '{print $4}'`

	sbatch --dependency=after:${j_qb} ssages.sh


The advantage of the latter method, with three scripts, is that it will avoid
conflict between modules, which may be present in the first example, depending
on how you have compiled Qbox and SSAGES.

HOOMD-blue
^^^^^^^^^^

Building
~~~~~~~~

HOOMD-blue supports SSAGES in releases v2.4.0 and later. HOOMD-blue must be built
with **MPI support enabled**, using the HOOMD CMake flag ``-DENABLE_MPI=ON``.
It is recommended to use a Python virtual environment with HOOMD-blue, as
described in `the HOOMD compilation instructions <https://hoomd-blue.readthedocs.io/en/stable/installation.html#compile-hoomd-blue>`_.
With the virtual environment active, building HOOMD-blue and SSAGES from source
can use Python to appropriately set relevant CMake variables.
For example (using 8 ``make`` job slots),

.. code:: bash

    cd /path/to/HOOMD
    mkdir build
    cd build/
    python3 -m venv /path/to/new/virtual/environment --system-site-packages
    source /path/to/new/virtual/environment/bin/activate
    cmake .. -DCMAKE_INSTALL_PREFIX=$(python3 -c "import site; print(site.getsitepackages()[0])") -DENABLE_MPI=ON
    make install -j 8

    cd /path/to/SSAGES
    mkdir build
    cd build/
    cmake .. -DHOOMD=ON
    make -j 8

Alternatively, HOOMD can be built using ``make`` (without installation) using
the HOOMD CMake flag ``-DCOPY_HEADERS=ON``.
For more information on HOOMD CMake flags, see `the HOOMD documentation <https://hoomd-blue.readthedocs.io/en/stable/installation.html>`_.
In this case, you must specify the SSAGES CMake flag
``-DHOOMD_ROOT=/path/to/hoomd_installation/hoomd``,
pointing to the ``hoomd`` directory within the HOOMD-blue installation prefix.

.. code:: bash

    cd /path/to/HOOMD
    mkdir build
    cd build/
    cmake .. -DCOPY_HEADERS=ON -DENABLE_MPI=ON
    make -j 8

    cd /path/to/SSAGES
    mkdir build
    cd build/
    cmake .. -DHOOMD_ROOT=/path/to/HOOMD/build/hoomd
    make -j 8

It may be necessary to put the HOOMD build directory in your ``PYTHONPATH``,
so that your Python interpreter can find the ``hoomd`` module. Set this
variable with the following command.

.. code:: bash

    export PYTHONPATH="${PYTHONPATH}:/path/to/hoomd/build"

Running
~~~~~~~

HOOMD-blue offers a "half-step hook" API to support SSAGES.
This feature is automatically configured when SSAGES launches
the HOOMD-blue simulation.

The HOOMD-blue SSAGES user script is written in Python. This script should
contain necessary ``import`` statements, configure the simulation, and set
types, interactions, integrators, log outputs and so forth. However, the
simulation user script should *not* call ``hoomd.run(steps)`` as this will
be called within SSAGES.

To set
`HOOMD-blue command-line options <https://hoomd-blue.readthedocs.io/en/v2.5.0/command-line-options.html>`_,
use the SSAGES JSON input file. Set the key ``"args"`` with a string of command
line options that will be passed to ``hoomd.context.initialize()``.
.. _Write-your-own-method:

Write Your Own CVs and Methods
==============================

One of the basic design goals of SSAGES is that it should be easily extensible.
To this end, it provides intuitive and simple tools to implement new collective
variables (CVs) and new advanced sampling methods. This section covers the
basic steps to implement a new CV or a new Method. Let us first start with the
implementation of a new CV. The techniques to implement a new Method are
covered :ref:`below <write_new_method>`.

.. _write_new_CV:

How to Write a New CV
---------------------

Each CV consists of two components: a header file and a schema file. The header
file contains the source code for the calculation of the CV and the schema file
describes the properties of the CV in a simple JSON format. Finally, you will
have to make SSAGES aware of the new CV by incorporating it into the core
classes.

The CV Header File
^^^^^^^^^^^^^^^^^^

Each CV in SSAGES is implemented as a child of the class ``CollectiveVariable``.
The header file should be placed in the directory ``src/CVs/`` and has to
(re)implement the following functions:

:code:`void Initialize(const Snapshot&)` (optional)
    Called during the pre-simulation phase. It is typically used
    to allocate or reserve memory.

:code:`void Evaluate(const Snapshot&)`
    Evaluation of the CV based on a simulation snapshot. This function should
    calculate both the value and gradient of the CV. The gradient should be a
    vector of length :math:`n`, where :math:`n` is the number of atoms
    in the Snapshot. Each element in the vector is the derivative of the CV with
    respect to the corresponding atom's coordinates. This method is called
    in the post-integration phase of every iteration.

:code:`double GetValue() const`
    Return the current value of the CV, as calculated in
    ``Evaluate(const Snapshot&)``.

:code:`double GetPeriodicValue() const`
    Return the current value of the CV, taking periodic boundary conditions into
    account. An example would be an angular CV which is bound to the region
    :math:`(-\pi,\pi]`. In this case, ``GetValue()`` could return *any* angle,
    while ``GetPeriodicValue()`` should return the angle mapped back into the
    region :math:`(-\pi,\pi]`. If the CV does not use periodic boundaries, this
    function should return the same value as ``GetValue()``.

:code:`const std::vector<Vector3>& GetGradient() const`
    Return the gradient of the CV (see ``Evaluate(const Snapshot&)`` for how the
    gradient is defined).

:code:`const std::array<double, 2>& GetBoundaries() const`
    Return a two-element array containing the lower and the upper boundary for
    the CV.

:code:`double GetDifference(double Location) const`
    Return the distance of the current value of the CV from a specified
    location, taking periodic boundary conditions into account. If the CV does
    not use periodic boundary conditions, the return value should simply be
    ``GetValue() - Location``.

The CV Schema File
^^^^^^^^^^^^^^^^^^

Together with the header file that contains the source code of the CV, you will
have to provide a schema file to make the CV accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/CVs/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of ``type`` should be set to ``object``.

varname
    The name of your new CV, of the form ``ExampleCV``.

properties
    The properties contain the ``type`` which is the internal name of the CV and
    a set of other properties that have to be supplied to the constructor of the
    CV.

required
    A list containing the required properties. Optional parameters to the CV
    constructor are not listed here.

additionalProperties
    A boolean value allowing unlisted JSON members to be included in input files.

Integrate the New CV into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header and schema files, there are two more
steps in order to make SSAGES aware of the newly included CV.

To include your new CV, you have to edit the file
``src/CVs/CollectiveVariable.cpp``, and

1. ``#include`` your CV header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildCV()``. The if-test checks for the
   type of CV selected and returns the respective ``Build(json, path)``
   function from the new CV.

.. _write_new_method:

How to Write a New Method
-------------------------

Each method consists of three components: a header file (``.h``),
a source file (``.cpp``), and a schema file (``.json``).
The header and source files contain the main code for the method
and the schema file describes the properties of the method in a simple JSON
format. Finally, you will have to make SSAGES aware of the new method by
incorporating it into the core classes.

The Method Header File
^^^^^^^^^^^^^^^^^^^^^^

Each method in SSAGES is implemented as a child of the class ``Methods``.
The header file should be placed in the directory ``src/Methods`` and has to
declare the following functions:

:code:`void PreSimulation(Snapshot* snapshot, const CVList& cvs)`
    Called before the method actually runs. It is typically used
    to allocate or reserve memory.

:code:`void PostIntegration(Snapshot* snapshot, const CVList& cvs)`
    Called after each MD integration step. This is where the heart of your
    method should go. By using Snapshot and the CVs, this function modifies
    the forces, positions, velocities, etc. appropriated by the new method.

:code:`void PostSimulation(Snapshot* snapshot, const CVList& cvs)`
    Called at the end of the simulation run. Use it to close files
    your method opened, to write out data that the method is storing, etc.

:code:`void Build`
	Called upon instantiation of SSAGES to build your method. Standard parts
	include reading the JSON input and setting variables.

Any other variables or functions that your new method will need to use
must be declared here (or defined, if simple enough).

The Method Source File
^^^^^^^^^^^^^^^^^^^^^^

The source file should be placed in the directory ``src/Methods`` and has to
(re)implement the functions defined above: ``PreSimulation()``,
``PostIntegration()``, ``PostSimulation``, and ``Build()``. Any functions
that were declared in the header file need to be defined, as well.
For String Method variants, edit the logic framework in
``src/Methods/StringMethod.cpp`` for the ``Build()`` function, instead of
your new source file.

The Method Schema File
^^^^^^^^^^^^^^^^^^^^^^

Together with the source code of the method, you will
have to provide a schema file to make the method accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/Methods/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of ``type`` should be set to ``object``.

varname
    The name of your new method, of the form ``ExampleMethod``.

properties
    The properties contain the ``type`` which is the internal name of the
    method and a set of other properties that have to be supplied to the
    constructor of the method.

required
    A list containing the required properties. Optional parameters to the method
    constructor are not listed here.

additionalProperties
    A boolean value allowing unlisted JSON members to be included in input files.

If your new method uses the String Method framework, you simply need to add a
new "flavor" of this method, defined in ``schema/Methods/string.method.json``.

Integrate the New Method into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header, source, and schema files, there are two more
steps in order to make SSAGES aware of the newly included method.

To include your new method, you have to edit the file
``src/Methods/Method.cpp``, and

1. ``#include`` your method header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildMethod()``. The if-test checks for
   the type of method selected and calls the respective
   ``Build(json, world, comm, path)`` function from the new method. A pointer
   to the newly created object should be stored in the variable named ``method``.

Finally, add the method ``.cpp`` file to CMakeLists.txt as a source.
.. _umbrella-sampling:

Umbrella Sampling
-----------------

Introduction
^^^^^^^^^^^^

Calculations of thermodynamic data and other properties rely on proper sampling of the configurational space.
However, the presence of energy barriers can prevent certain configurations from being sampled properly or even sampled at
all. Umbrella sampling is a simulation technique that helps to overcome those barriers and improve sampling
by applying a bias along a specified collective variable. The bias takes the form of a harmonic potential and is typically constant throughout a simulation.
Usually, a series of umbrella-sampled simulations are performed and analyzed together using the weighted histogram analysis method
(WHAM) :cite:`KUMAR19921011`.

The functional form of the artificial bias is

.. math::

	U_\text{umbrella} = \frac{1}{2} k \left(\xi - \xi_0\right)^2

where :math:`k` is the spring constant, :math:`\xi` is the current value of the collective variable and :math:`\xi_0` is the desired value of the collective variable.
The success of the umbrella sampling method depends on the correct choice of :math:`k` and :math:`\xi_0` for the different simulations. Suitable values of :math:`k` and :math:`\xi_0` depend on the distances between adjacent umbrellas, and the gradient of the free energy surface at :math:`\xi_0`.

For more information about the umbrella sampling method, the interested reader is referred to Ref. :cite:`KASTNER2011932`.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

The following parameters need to be set under ``"method"`` in the JSON input file:

.. code-block:: javascript

	"type" : "Umbrella"

The following options are available for Umbrella Sampling:

centers (required)
	Array of target values for each CV. This can either be an array of numbers, or
	an array of an array of numbers, for multiple walkers. If only a single array is
	defined, it is used across all walkers.

ksprings (required)
	Array of spring constants for each CV. This can either be an array of numbers, or
	an array of an array of numbers, for multiple walkers. If only a single array is
	defined, it is used across all walkers.

output_file (required)
	Output file name for umbrella sampling data. This can be a string or
	an array of strings for multiple walkers.

output_frequency (optional)
	Frequency of writing to output file. (Default: 1)

append (optional)
	Boolean value which will cause umbrella sampling to either append to
	an existing file or override it entirely. (Default: ``false``)

The umbrella sampling method can also be run with time dependent centers.
This is equivalent to the "steered MD" method. To do so, omit ``"centers"``
and use the following options instead:

centers0 (required)
	Array of initial target values for each CV. This can either be an array of numbers, or
	an array of an array of numbers, for multiple walkers. If only a single array is
	defined, it is used across all walkers.

centers1 (required)
	Array of final target values for each CV. This can either be an array of numbers, or
	an array of an array of numbers, for multiple walkers. If only a single array is
	defined, it is used across all walkers.

timesteps (required)
	The number of timesteps over which to scale the umbrella centers

By setting the above options, the umbrella method will linearly interpolate
between ``centers0`` and ``centers1`` in ``timesteps`` iterations.

.. _Umbrella_tutorial:

LAMMPS Tutorial
^^^^^^^^^^^^^^^

This tutorial will go through running Umbrella Sampling on an atomistic model
of butane using LAMMPS as the MD engine.
Umbrella sampling will be performed on the torsional CV of the butane C atoms.

The butane implementation in LAMMPS requires several modules to be added
before being linked to SSAGES.
To do this, return to your build directory and issue the following commands:

.. code-block:: bash

	make yes-molecule
	make yes-kspace
	make

which add the ``MOLECULE`` and ``KSPACE`` LAMMPS packages.
For more information on optional packages for LAMMPS, refer to
the `LAMMPS User Manual <https://lammps.sandia.gov/doc/Packages.html>`_.

The files for running this example can
be found in ``Examples/User/Umbrella/LAMMPS`` and consist of the following files:

``Butane_SSAGES.in``
	LAMMPS input file

``Butane.data``
	LAMMPS data file describing butane molecule.

``umbrella_input.json``
	Template JSON input containing information for one Umbrella Sampling simulation.

``umbrella_multiwalker_gen.py``
	Python script for creating ``multiwalker_umbrella.json`` input file. The total number of
	simulations and the ``centers`` values are controlled in this file.

Once in the directory, the appropriate ``.json`` file needs to be generated. A ``.json`` file
is already in the directory, ``umbrella_input.json``, which contains the CV information
and specifies the LAMMPS input files to be used. A single-walker umbrella simulation can
be run directly using

.. code-block:: bash

	ssages umbrella_input.json

The simulation will create an output file named ``umbrella.dat1`` containing the value of
the CV and the target value (the center) every 100 timesteps. From this histogram, the
local free energy can be calculated.

While it is possible to run Umbrella Sampling using a single walker, typically multiple
walkers (multiple umbrellas) are simulated. To run multiwalker Umbrella Sampling of butane,
you can generate an input file using the ``umbrella_multiwalker_gen.py`` script via

.. code-block:: bash

	python umbrella_multiwalker_gen.py

This will generate an input file called ``multiwalker_umbrella.json`` containing the
information from ``umbrella_input.json`` duplicated 12 times with varying values of
``centers``. These values correspond to the target values of the torsional CV.

To run multiwalker SSAGES issue the command:

.. code-block:: bash

	mpiexec -np 12 /path/to/SSAGES/build/ssages multiwaler_umbrella.json

This will run 12 different umbrella sampling simulations simultaneously.
Ideally, this example will be run in computing environment where each process can run
on a different processor. The example will still work if run on a users local desktop
or laptop machine, but the runtime of the code will be very large.

During the simulation 12 different output files will be generated, each containing the
iteration, target value of the corresponding 'center' CV, and the value of the CV at
the iteration number.

These output files can then be used to construct a complete free energy surface using
the WHAM algorithm :cite:`KUMAR19921011`. Though SSAGES does not currently contain its own implementation
of WHAM, there are many implementations available, such as that provided by the
Grossfield Lab :cite:`WHAM`.

HOOMD-blue Tutorial
^^^^^^^^^^^^^^^^^^^

This example uses the HOOMD-blue engine to run parallel simulations of a butane molecule.
The free energy is measured as a function of the dihedral angle between the terminal carbons.
The butane molecule has a backbone of four carbon atoms that `rotates into different conformations <https://chem.libretexts.org/Textbook_Maps/Organic_Chemistry/Supplemental_Modules_(Organic_Chemistry)/Chirality/Stereoisomers/Butane_Conformers>`_ (*anti*, *gauche*, and *eclipsed*).
We wish to extract the free energy of this rotation, to know the energy cost of any angle between -180 degrees and 180 degrees.

This example uses Umbrella Sampling with the weighted histogram analysis method (WHAM).
The WHAM tool developed by Alan Grossfield :cite:`WHAM` is used to determine the free energy from the biased sampling we perform.
Disclaimer: The parameters of this simulation (number of walkers, strength of bias potential springs, length of run, etc.) may not provide ideal sampling for this example problem, and improvements to this code are welcomed.

The files for running this example can be found in ``Examples/User/Umbrella/HOOMD``.

Sample output files from this example code are in the ``Examples/User/Umbrella/HOOMD/sample_outputs`` folder.

**Running the Example Script:**

1. Modify HOOMD-blue script: Set desired parameters (e.g. ``kT``) in ``Butane_SSAGES.py``

2. Modify input generator: Set the parameters in ``umbrella_multiwalker_gen.py`` and ``umbrella_input.json``. Important parameters:

``umbrella_multiwalker_gen.py``:

	* ``nwalkers`` is the number of walkers, determining how many independent simulations will be run.

``umbrella_input.json``:

	* ``ksprings`` gives the bias potential spring strength.
	* ``hoomd_steps`` changes the length of the run.

	Most of the other parameters are used to define the system and collective variables and should not be changed.

3. Generate inputs:

.. code-block:: bash

	python umbrella_multiwalker_gen.py

4. Run SSAGES, replacing "nwalker" with the number of walkers specified previously:

.. code-block:: bash

	mpiexec -np nwalker /path/to/SSAGES/build/ssages multiwalker_umbrella_input.json

5. Analyze data:

	Download the WHAM code available here :cite:`WHAM`.
	Compile the program using the instructions and documentation provided.
	It is recommended to read `this talk about the theory and practice of WHAM <http://membrane.urmc.rochester.edu/sites/default/files/wham/wham_talk.pdf>`_.

	a) Call wham: The script ``wham_analysis.sh`` contains a set of parameters for calling ``wham``.
		This requires that the executable ``wham`` is in this directory.

	.. code-block:: bash

		./wham_analysis.sh

	b) Run visualization script:

	.. code-block:: bash

		python wham_visualization.py

	The script ``wham_visualization.py`` will read the output data files from SSAGES and the ``wham`` software to produce sets of figures similar to those in the talk linked above.
	The visualization outputs include:

	* ``cv_vs_time.png`` plots the collective variable (dihedral angle) over
	  time. This helps check that enough autocorrelation times have passed.
	* ``histogram_trajectories.png`` shows a histogram from each of the
	  trajectories and the regions of the CV that were sampled.
	* ``histogram_combined.png`` shows a histogram summed over all trajectories
	  to ensure that the entire range of angles were sampled.
	* ``wham_free_energy.png`` is the free energy as a function of the dihedral
	  angle.

Developers
^^^^^^^^^^

* Hythem Sidky
* Benjamin Sikora
.. _elastic-band:

Elastic Band
------------

Introduction
^^^^^^^^^^^^

There are many methods, several of which are included in SSAGES, to calculate
transition pathways between metastable states. One kind of pathway between
states in the *Minimum Energy Pathway* (MEP), quite simply the lowest energy
pathway a system can take between these states. An MEP has the condition that
the force everywhere along the pathway points only along the path, that is, it
has no perpendicular component. By finding the MEP, one also finds the saddle
points of the potential energy surface, as they are by definition the maxima of
the MEP. The *Nudged Elastic Band* (NEB) method is a popular and efficient
method to calculate the MEP between the initial and final state of a transition
:cite:`HENKELMAN20009978,HENKELMAN20009901`.

The method involves the evolution of a series of images connected by a spring
interaction (hence the "elastic" nature of the band). The force acting on the
images (a combination of the spring force along the band and the true force
acting perpendicular to the band) is minimized to ensure convergence to the MEP.
The **nudged** nature of NEB refers to a force projection that ensures the
spring forces do not interfere with the elastic band converging to the MEP, as
well as that the true force does not alter the distribution of images along the
band (that is, it ensures all the images do not fall into the metastable states).
This projection is accomplished by using the parallel portion of the spring
force and the perpendicular portion of the true force. In this way, the spring
forces act similarly to reparameterization schemes common to the string method.

Full mathematical background is available in :cite:`HENKELMAN20009978` and
:cite:`HENKELMAN20009901`, but a brief overview is given here. The
band is discretized as a series of N+1 images, and the force on each image is
given by:

.. math::

	F_{i} = F_{i,\parallel}^{s} - \nabla E(R_{i})_{\perp}

Where :math:`F_{i}` is the total force on the image, :math:`F_{i,\parallel}^{s}`
refers to the parallel component of the spring force on the
:math:`i^\text{th}` image, and
:math:`\nabla E(R_{i})_{\perp}` is the perpendicular component of the gradient
of the energy evaluated at each image :math:`R_{i}`. The second term on the
right hand side is the "true force" and is evaluated as:

.. math::

	\nabla E(R_{i})_{\perp} = \nabla E(R_{i}) - \nabla E(R_{i})\cdot\hat{\tau_{i}}

The term :math:`\hat{\tau_{i}}` represents the normalized local tangent at the
:math:`i^\text{th}` image, and thus this equation states simply that the perpendicular component
of the gradient is the full gradient minus the parallel portion of the gradient.
There are different schemes available in literature to evaluate the tangent
vector :cite:`HENKELMAN20009978`. The "spring force" is calculated as:

.. math::

	F_{i,\parallel}^{s} =
	k \left( \lvert R_{i+1} - R_{i} \rvert -
	         \lvert R_{i} - R_{i-1} \rvert \right) \hat{\tau_{i}}

Where :math:`k` is the spring constant, which can be different for each image of
the band. One can evolve the images with these forces according to any number
of schemes - a straightforward Verlet integration scheme is used in the SSAGES
implementation, described below.

Algorithmically, the NEB method is implemented in SSAGES as follows:

1. An initial band is defined between the two states of interest. This can be
   defined however one wishes; often it is simply a linear interpolation through
   the space of the collective variables. In fact, the ends of the band need
   not necessarily be in the basins of interest; the method should allow the
   ends to naturally fall into nearby metastable basins.

2. For each image of the band, a molecular system with atomic coordinates that
   roughly correspond to the collective variables of that image is constructed.
   A period of equilibration is performed to ensure that the underlying systems'
   CVs match their respective band images.

3. The gradient is sampled over a user-defined period of time and intervals,
   this being the only quantity with statistical variance that needs to be
   averaged over.

4. When sufficient sampling of the gradient is done, the band is updated one
   time-step forward with a simple Verlet scheme.

Steps 2--4 are iterated upon, leading to convergence of the method and the MEP.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

To construct an EB input file, the following options are available. A
complete EB JSON file will inherit some of its inputs from the String
schema (for parameters common to all string methods).
The options unique to EB are:

equilibration_steps
	The number of MD steps to simply perform umbrella sampling without
	invoking the NEB method. A sufficiently long number of steps ensures
	that the underlying molecular systems have CVs close to the CVs of their
	associated image on the band.

evolution_steps
	The number of steps to perform the NEB over; the band is updated after
	evolution steps times the number of samples total MD steps. A new value
	of the gradient is harvested every time the number of MD steps taken is
	an integer multiple of evolution steps.

kstring
	The constant used in calculating the spring force at each image. It
	can be specified uniquely for each image. Please notice its difference
	from ksprings.

From the String schema, the options are:

type
	This parameter identifies that a String-type method is being used, and
	thus should be set to ``"String"``.

flavor
	This parameter identifies the specific kind of string-type method
	being used; for EB, it should be set to ``"ElasticBand"``.

centers
	This parameter assigns a location in CV space for each individual image
	along the string/elastic band. This should be an array with size equal to
	the total number of images, with each entry consisting of an array with size
	equal to the number of CVs used for the elastic band method.

tolerance
	This is a tolerance threshold that can be set to trigger the end of
	the method; it is a percentage by which, if no node CV changes by this
	percentage, the method will end. It must be specified as an array with
	one entry for each CV desired.

max_iterations
	A complementary stopping criterion can be specified; the method will
	stop if it undergoes this many iterations of the string method.

ksprings
	A unique spring constant must be defined for each CV; its purpose is
	described above.

frequency
	The frequency of each integration step. This should almost always be
	set to 1.

.. _EB_tutorial:

Tutorial
^^^^^^^^

This tutorial will walk you step by step through the user example provided with
the SSAGES source code that runs the NEB method on the alanine dipeptide using
LAMMPS. First, be sure you have compiled SSAGES with LAMMPS. Then, navigate to
the ``Examples/User/ElasticBand/ADP`` subdirectory. Now, take a moment
to observe the ``in.ADP_Test and data.input`` files in order to familiarize
yourself with the system being simulated.

The next two files of interest are the ``Template_Input.json`` input file and the
``Input_Generator.py`` script. Both of these files can be modified in your
text editor of choice to customize the inputs, but for this tutorial, simply
observe them and leave them be. ``Template_Input.json`` contains all the information
necessary to fully specify one driver; ``Input_Generator.py`` copies this
information a number of times specified within the script (for this tutorial,
22 times) while also linearly interpolating through the start and end states
defined in the script and substituting the correct values into the "centers"
portion of the method definition. Execute this script as follows:

.. code-block:: bash

	python Input_Generator.py

You will produce a file called ``ElasticBand.json``. You can also open this file to
verify for yourself that the script did what it was supposed to do. Now, with
your JSON input and your SSAGES binary, you have everything you need to perform
a simulation. Simply run:

.. code-block:: bash

	mpiexec -np 22 ./ssages ElasticBand.json

Soon, the simulation will produce a ``node-X.log`` file for each driver, where
X is the number specifying the driver (in this case, 0-21 for our 22 drivers).
Each one will report the following information, in order: the node number, the
iteration number, and for each CV, the current value of the band CV as well as
the current value of the CV calculated from the molecular system.

Allow your system to run for the specified number of iterations (1000 for this
tutorial). The last line of every node file can be analyzed to view the last
positions of each image of the elastic band.

Developer
^^^^^^^^^

* Benjamin Sikora
.. _cvs:

Collective Variables
====================

Collective variables (CVs) are arbitrary differentiable functions of the
:math:`3N` Cartesian coordinates of the atoms in a simulation. These
usually represent some structurally, thermodynamically, or chemically
meaningful quantity along which advanced sampling can be performed. Listed
below are the collective variables currently supported in SSAGES, along with a
brief description and information on how the syntax is implemented in SSAGES.
In addition to specific properties for each CV, a name property may be
defined for any CV which is used to reference the CV within the methods of SSAGES.

.. code-block:: javascript

	"name" : "mycvname"

The name specified for a CV must be unique. It is possible, however, to omit
this and simply reference a CV by its numerical index.

.. note::

	We tacitly assume we are working with spherical atoms that join together to form molecules.

Angle
-----

Description
^^^^^^^^^^^

This CV calculates the bend angle, in radians, formed between three selected atoms :math:`i,j,k`,

.. math::

	\xi = \cos^{-1}\left(\frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{kj}}{\Vert \mathbf{r}_{ij} \Vert \Vert \mathbf{r}_{kj} \Vert} \right).

This can be helpful when probing the conformations of a molecule to understand its stable and metastable states. Angles do not have to be defined between bonded atoms.

Example
^^^^^^^

.. code-block:: javascript

	{
		"type" : "Angle",
		"atom_ids" : [0, 1, 2]
	}

.. warning::

	The angle must be between three individual particles rather than the centers-of-mass of particle groups.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"Angle"``.


.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must contain three integers consisting of the atom ID forming the angle of interest.

ANNCV
-----

Description
^^^^^^^^^^^

This CV takes scaled (specified by ``scaling_factor``) Cartesian coordinates of a group of atoms (specified by ``atomids``) as inputs to a neural network (its number of nodes, connection weights, and activation functions are specified by ``num_nodes``, ``coeff_file``, ``activations``, respectively), computes one component (specified by ``out_index``) of the final neural network outputs as the CV value. The coefficients of the neural network can be obtained from any feed-forward neural networks trained with Cartesian coordinates as inputs.  Examples of the neural networks include the encoders of the autoencoders in the `MESA framework <https://github.com/weiHelloWorld/accelerated_sampling_with_autoencoder>`_, or `State-free Reversible VAMPnets <https://github.com/hsidky/srv>`_.

In the following example, we define an ANN CV which takes Cartesian coordinates of atoms ``[2, 5, 7, 9, 15, 17, 19]``, scaled by factor 0.5, as inputs to a neural network with node numbers ``[21, 40, 2]`` and activation functions ``["Tanh", "Tanh"]``, and weights defined in file ``autoencoder_info_1.txt`` (which stores weights for the neural network), and outputs two components (marked as index 0 and 1) as CVs.

Example
^^^^^^^

.. code-block:: javascript

	"CVs": [
				{
					"type": "ANNCV",
					"atom_ids": [2, 5, 7, 9, 15, 17, 19],
					"scaling_factor": 0.5,
					"num_nodes": [21, 40, 2],
					"activations": ["Tanh", "Tanh"],
					"index": 0,
					"coeff_file": "autoencoder_info_1.txt"
				},
				{
					"type": "ANNCV",
					"atom_ids": [2, 5, 7, 9, 15, 17, 19],
					"scaling_factor": 0.5,
					"num_nodes": [21, 40, 2],
					"activations": ["Tanh", "Tanh"],
					"index": 1,
					"coeff_file": "autoencoder_info_1.txt"
				}
			]
			

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type": "ANNCV"

Selects this collective variable.


.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must contain integers consisting of the atom ID for the inputs of ANN.

.. code-block:: javascript

	"scaling_factor"

Property ``scaling_factor`` is the scaling factor of the inputs.

.. code-block:: javascript

	"num_nodes"

Property ``num_nodes`` defines the number of nodes for each layer of the neural network.

.. code-block:: javascript

	"activations"

Property ``activations`` defines the activation functions for each layer of the neural network.

.. code-block:: javascript

	"coeff_file"

Property ``coeff_file`` defines the file which stores weights for the neural network.

.. code-block:: javascript

	"index"

Property ``index`` defines the output index we want to use for CV.

Box Volume
----------

Description
^^^^^^^^^^^

The current volume of a simulation box is an important parameter determining the thermodynamic state. Constant-pressure simulations where volume information is recorded may be reweighted according to standard methods :cite:`CONRAD199851`. This CV calculates the box volume as the determinant of the Parrinello-Rahman matrix :math:`\mathbf{H}`,

.. math::

	\xi = \det\left( H_{ij} \right)

Example
^^^^^^^

.. code-block:: javascript
	
	{
		"type" : "BoxVolume"
	}

.. warning::

	Non-orthorhombic boxes are currently not supported. Only Gromacs and LAMMPS are currently supported

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"BoxVolume"``.

Gyration Tensor
---------------

Description
^^^^^^^^^^^

This CV calculates quantities derived from the symmetric *mass-weighted*
gyration tensor of a group of :math:`N` atoms defined as,

.. math::

	\mathbf{S} = \frac{1}{\sum_{i=1}^{N}{m_i}}\sum_{i=1}^{N}{m_i \left( \mathbf{r}_i - \mathbf{r}_\mathrm{COM}\right) \otimes \left( \mathbf{r}_i - \mathbf{r}_\mathrm{COM}\right)}

where :math:`m_i` is the mass and :math:`\mathbf{r}_i` is the vector of
coordinates of the :math:`i^{\mathrm{th}}` atom, :math:`\mathbf{r}_\mathrm{COM}`
is the vector of the center of mass of all :math:`N` atoms in the group, and
:math:`\otimes` is the outer, or tensor, product.

The eigenvalues of the radius of gyration tensor are particularly useful as
collective variables which quantify the conformation of a molecule (such as a
long polymer) or the shape of a given assembly of molecules. With eigenvalues
of :math:`\lambda_x^2,~\lambda_y^2,~\lambda_z^2` (in increasing order)
defined in the frame of the principal axes of inertia,
the following quantities may be computed:

Radius of Gyration (Squared)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

	R_g^2 = \lambda_x^2 + \lambda_y^2 + \lambda_z^2

Principal Moment
~~~~~~~~~~~~~~~~

.. math::

	\lambda_i^2,\ i \in \{x,y,z\}

Asphericity
~~~~~~~~~~~

.. math::

	b = \lambda_z^2 - \frac{1}{2}\left(\lambda_x^2 + \lambda_y^2 \right)

Acylindricity
~~~~~~~~~~~~~

.. math::

	c = \lambda_y^2 - \lambda_x^2

Shape Anisotropy
~~~~~~~~~~~~~~~~

.. math::

	\kappa^2 = \frac{3}{2}\frac{\lambda_x^4+\lambda_y^4+\lambda_z^4}{\left(\lambda_x^2+\lambda_y^2+\lambda_z^2\right)^2}-\frac{1}{2}



Example
^^^^^^^

This example computes the shape anisotropy of a ten-atom group.

.. code-block:: javascript

	"type" : "GyrationTensor",
	"atom_ids" : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
	"component" : "shapeaniso"

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript
	
	"type"

Property ``type`` must be set to string ``"GyrationTensor"``.

.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must be an array of integers containing the atom IDs which will enter the calculation.

.. code-block:: javascript

	"component"

Property ``component`` must be a string defining the gyration tensor component of interest.
Valid options are ``"Rg"``, ``"principal1"``, ``"principal2"``, ``"principal3"``, ``"asphericity"``,
``"acylindricity"``, or ``"shapeaniso"``.

Optional
~~~~~~~~

.. code-block:: javascript

    "dimension"

Property ``dimension`` is a 3-element array of booleans specifying which
Cartesian components to include in the calculation. If left unspecified, all
three xyz components will be used.

Particle Coordinate
-------------------

Description
^^^^^^^^^^^

This CV calculates the :math:`x`, :math:`y` or :math:`z` position of the center of mass for a
group of atoms.

.. math::

	\xi = \frac{1}{\sum_i{m^i}}\sum_{i=1}^{N}{r_\alpha^i}\ \ \ \alpha \in {x,y,z}

Example
^^^^^^^

.. code-block:: javascript

	{
		"type" : "ParticleCoordinate",
		"atom_ids" : [1, 5, 6, 10],
		"dimension" : "x"
	}


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript
	
	"type"

Property ``type`` must be set to string ``"ParticleCoordinate"``.

.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must be an array of integers containing the atom IDs which will enter the calculation.

.. code-block:: javascript

	"dimension"

Property ``dimension`` must be a string defining the Cartesian component of interest ``"x"``, ``"y"``, or ``"z"``.

Pairwise
--------

Description
^^^^^^^^^^^

This CV calculates a variety of pairwise properties. The functions (kernels) used are continous analogs for otherwise discontinuous CVs. If parameters are chosen judiciously, these kernels can be used in place of some standard, discontinuous CVs. A Gaussian kernel can emulate a count of nearest neighbors; a switching function kernel can emulate a coordination number.

.. math::

	\xi = \sum_{i \in A}\sum_{i \in B}{f_{ij}}

where :math:`f_{ij}` is a pairwise function for atoms :math:`i` and :math:`j`. are at a distance of the center of the Gaussian, :math:`r_{ij}=\mu`, and decreases to zero as the distance deviates away from :math:`\mu`.

Example
^^^^^^^

This example uses a Gaussian pairwise kernel to compute contributions from contact-type interactions between two atoms of size 1.0.

.. code-block:: javascript
	
	{
		"type" : "Pairwise",
		"group1" : [1, 5],
		"group2" : [2, 3, 4, 6, 7, 8],
		"kernel" : {
			"type" : "gaussian",
			"mu" : 1.0,
			"sigma" : 0.2
		}
	}



Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript
	
	"type"

Property ``type`` must be set to string ``"Pairwise"``.

.. code-block:: javascript
	
	"group1"

Property ``group1`` must be an array of integers containing the atom IDs in the first set.

.. code-block:: javascript
	
	"group2"

Property ``group2`` must be an array of integers containing the atom IDs in the second set.

.. note::

	Atoms can exist in both ``group1`` and ``group2`` simultaneously. Contacts are automatically
	skipped if :math:`i = j`.

.. code-block:: javascript
	
	"kernel"

Property ``kernel`` must be an object defining the properties of the pairwise kernel function and its associated properties.

Pairwise Kernels
~~~~~~~~~~~~~~~~

Gaussian Function
*****************

The Gaussian function is defined as:

.. math::

	g_{ij} = e^{-\frac{\left(r_{ij} - \mu\right)^2}{2\sigma^2}}.

This type of kernel is useful to select between conformations which have a different position of (e.g.) neighbors and next nearest neighbors in a particle cluster. Selection of particle separations approximates a math:`\delta` distribution.

Properties
++++++++++

.. code-block:: javascript
	
	"mu"

Property ``mu`` is required and must be numeric.

.. code-block:: javascript
	
	"sigma"

Property ``sigma`` is required and must be numeric.

Rational Switching Function
***************************

The rational switching function is defined as:

.. math::

	s_{ij} = \frac{1-\left(\frac{r_{ij} - d_0}{r_0}\right)^n}{1-\left(\frac{r_{ij} - d_0}{r_0}\right)^m}.

This quantity is useful for measuring how many atoms in group 2 occupy a spherical shell around atoms in group 1. The form is chosen so that the variable is continuous and differentiable. Through tuning :math:`n` and :math:`m` this can be made arbitrarily close to a Heaviside switching function.

Properties
++++++++++

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"rationalswitch"``.

.. code-block:: javascript

	"d0"

Property ``d0`` is required and must be numeric.

.. code-block:: javascript

	"r0"

Property ``r0`` is required and must be numeric.

.. code-block:: javascript

	"n"

Property ``n`` is required and must be an integer.

.. code-block:: javascript

	"m"

Property ``m`` is required and must be an integer.

Particle Position
-------------------

Example
^^^^^^^

.. code-block:: javascript

	{
		"type" : "ParticlePosition",
		"atom_ids" : [1, 5, 6, 10],
		"dimension" : [true, false, true],
		"position" : [3.51, 6.66, 2.14]
	}

Description
^^^^^^^^^^^

This CV calculates the distance of the center of mass of a group of atoms
from a particular point in Cartesian space.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript
	
	"type"

Property ``type`` must be set to string ``"ParticlePosition"``.

.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must be an array of integers containing the atom IDs which
will enter the calculation.

.. code-block:: javascript

	"position"

Property ``position`` must be a 3-element array of numbers defining the reference
point in the simulation box.

Optional
~~~~~~~~

.. code-block:: javascript

    "dimension"

Property ``dimension`` is a 3-element array of booleans specifying which
Cartesian components to include in the calculation. If left unspecified, all
three xyz components will be used.

Particle Separation
-------------------


Description
^^^^^^^^^^^

This CV calculates the distance between the centers of mass of two groups of
atoms. The variable is unsigned, as the distance is a magnitude.

Example
^^^^^^^

.. code-block:: javascript

    {
        "type" : "ParticleSeparation",
        "group1" : [1],
        "group2" : [5, 6, 10]
    }


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

    "type"

Property ``type`` must be set to string ``"ParticleSeparation"``.

.. code-block:: javascript

    "group1"

Property ``group1`` must be an array of integers containing the atom ID(s) which
make up the first group of atoms. The CV will calculate the distance between
the center of mass of this group and the group defined by property ``group2``.

.. code-block:: javascript

    "group2"

Property ``group2`` must be an array of integers containing the atom ID(s) which
make up the second group of atoms. The CV will calculate the distance between
the center of mass of this group and the group defined by property ``group1``.

Optional
~~~~~~~~

.. code-block:: javascript

    "dimension"

Property ``dimension`` is a 3-element array of booleans specifying which
Cartesian components to include in the calculation. If left unspecified, all
three xyz components will be used.

Polymer Rouse Modes
-------------------

Description
^^^^^^^^^^^

This CV calculates the magnitude of a given Rouse mode for a set of atoms as

.. math::

    X_p = \sqrt{\mathbf{X}_p\cdot\mathbf{X}_p},

with the :math: `p` th Rouse mode defined as

.. math::

    \mathbf{X}_p = \sqrt{\frac{c_p}{N}}\sum_{i=1}^N \mathbf{R}_i \cos \Bigl[\frac{p\pi}{N}\bigl(i-\frac{1}{2}\bigr) \Bigr],

where :math: `N` is the number of groups or beads comprising the polymer, :math: `\mathbf{R}_i` is the center-of-mass of the :math: `i` th bead, and :math: `c_p` is a constant equal to 1 for :math: `p=0` and equal to 2 for :math: `p=1,\cdots,N-1`. This CV can be helpful to bias the conformations of both moderate-size and long-chain proteins and polymers.


Example
^^^^^^^

.. code-block:: javascript

    {
        "type": "RouseMode",
        "mode": 1,
        "groups":  [
                    [ 1, 2, 3, 4, 5],
                    [ 6, 7, 8, 9,10],
                    [11,12,13,14,15],
                    [16,17,18,19,20],
                    [21,22,23,24,25],
                    [26,27,28,29,30],
                    [31,32,33,34,35],
                    [36,37,38,39,40],
                    [41,42,43,44,45],
                    [46,47,48,49,50]
                   ]
    }


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
^^^^^^^^

.. code-block:: javascript

    "type"

Property ``mode`` must be set to string ``"RouseMode"``.

.. code-block:: javascript

    "groups"

Property ``groups`` is an array of arrays containing the atom IDs (as integers) that comprise the discretized polymer beads. The number of groups provided implicitly defines :math: `N`, the number of polymer beads.

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer indicating the index of the desired Rouse mode. Valid values range from 0 up to one less than the number of groups, or `0,\cdots, N-1`.

Torsional Angle
---------------

Description
^^^^^^^^^^^

This CV calculates the dihedral angle, in radians, formed by four atoms :math:`i,j,k,l`.
It is computed as in :cite:`BLONDEL19961132`,

.. math::

	\xi = \tan^{-1}\left( \frac{\left[(r_{lk} \times r_{jk}) \times (r_{ij} \times r_{jk}) \right] \cdot \frac{r_{jk}}{\Vert r_{jk}\Vert}}{(r_{lk} \times r_{jk}) \cdot (r_{ij} \times r_{jk}) } \right).

Specifically, the function ``atan2`` is used for the inverse tangent calculation to yield a four-quadrant angle.

.. warning::

	The torsional angle can only be defined between four atoms rather than four groups of atoms.


Example
^^^^^^^

.. code-block:: javascript

	{
		"type" : "Torsional",
		"atom_ids" : [1, 5, 6, 10]
	}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript
	
	"type"

Property ``type`` must be set to string ``"Torsional"``.

.. code-block:: javascript

	"atom_ids"

Property ``atom_ids`` must be an array of 4 integers containing the atom IDs which
form the dihedral.

Alpha Helix RMSD
----------------

Description
^^^^^^^^^^^

This CV calculates alpha helix character by comparision to an "ideal" alpha
helix structure composed of 6 amino acids. This is computed by performing a
summation over all possible sequences of 6 consecutive amino acids in the
segment of interest:

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of switching functions applied to the RMSD rather than to coordinates; in future versions, the user will be able to choose custom parameters for the switching function.

.. note::

	Unlike the simpler CVs discussed above, this one takes atomic labels in the form of a reference PDB structure. This is true of all protein-like CVs below which compare to a reference structure.

.. warning::

	Since the definition of this CV uses nanometers as a unit length, you must specify the ``unitconv`` parameter, as outlined below, in order to apply this CV when that is not the base unit of length.

Example
^^^^^^^

.. code-block:: javascript

	{
            "type" : "AlphaRMSD",
            "residue_ids" : [3, 21],
            "reference" : "reference_structure.pdb",
            "unitconv" : 10
	}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"AlphaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for alpha helix character. For now, all
residues in the system must be numbered in increasing order, even if they belong
to separate chains. For example, if your system has two chains of 20 amino acids
each, the first amino acid in the second chain should be numbered 21.

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal alpha helix
reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10.

Anti Beta RMSD
----------------

Description
^^^^^^^^^^^

This CV calculates anti beta-sheet character by comparision to an "ideal" anti
beta-sheet structure composed of 6 amino acids. This is computed by performing a
summation over all possible sequences of 6 amino acids, consisting of two
segments of 3 consecutive amino acids each, in the region of interest.

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of switching functions applied to the RMSD rather than to coordinates; in future versions, the user will be able to choose custom parameters for the switching function.

.. note::

	Unlike the simpler CVs discussed above, this one takes atomic labels in the form of a reference PDB structure. This is true of all protein-like CVs below which compare to a reference structure.

.. warning::

	Since the definition of this CV uses nanometers as a unit length, you must specify the ``unitconv`` parameter, as outlined below, in order to apply this CV when that is not the base unit of length.

Example
^^^^^^^

.. code-block:: javascript

	{
            "type" : "AntiBetaRMSD",
            "residue_ids" : [3, 21],
            "reference" : "reference_structure.pdb",
            "unitconv" : 10,
            "mode" : 0
	}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"AntiBetaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for anti beta-sheet character. For now, all
residues in the system must be numbered in increasing order, even if they belong
to separate chains. For example, if your system has two chains of 20 amino acids
each, the first amino acid in the second chain should be numbered 21.

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal anti
beta-sheet reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10.

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer specifying whether to calculate beta-sheets
formed only between residues on the same chain (intra) or only between residues
on separate chains (inter). If ``mode`` is set to 0, both modes will be used.
A value of 1 selects for the intra mode; a value of 2 selects for inter mode.

Parallel Beta RMSD
------------------

Description
^^^^^^^^^^^

This CV calculates anti beta-sheet character by comparision to an "ideal"
parallel beta-sheet structure composed of 6 amino acids. This is computed by
performing a summation over all possible sequences of 6 amino acids, consisting
of two segments of 3 consecutive amino acids each, in the region of interest.

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of switching functions applied to the RMSD rather than to coordinates; in future versions, the user will be able to choose custom parameters for the switching function.

.. note::

	Unlike the simpler CVs discussed above, this one takes atomic labels in the form of a reference PDB structure. This is true of all protein-like CVs below which compare to a reference structure.

.. warning::

	Since the definition of this CV uses nanometers as a unit length, you must specify the ``unitconv`` parameter, as outlined below, in order to apply this CV when that is not the base unit of length.

Example
^^^^^^^

.. code-block:: javascript

	{
		"type" : "ParallelBetaRMSD",
		"residue_ids" : [3, 21],
		"reference" : "reference_structure.pdb",
		"unitconv" : 10,
		"mode" : 0
	}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"ParallelBetaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for parallel beta-sheet character. For now,
all residues in the system must be numbered in increasing order, even if they
belong to separate chains. For example, if your system has two chains of 20
amino acids each, the first amino acid in the second chain should be numbered
21.

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal parallel
beta-sheet reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10.

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer specifying whether to calculate beta-sheets
formed only between residues on the same chain (intra) or only between residues
on separate chains (inter). If ``mode`` is set to 0, both modes will be used.
A value of 1 selects for the intra mode; a value of 2 selects for inter mode.
.. _artificial-neural-network-sampling: 

Artificial Neural Network Sampling
----------------------------------

Introduction
^^^^^^^^^^^^

Artificial Neural Network (ANN) sampling is a free energy sampling method which uses
ANNs to generate an on-the-fly adaptive bias capable of rapidly resolving free
energy landscapes. It is a recent method proposed by Sidky & Whitmer
:cite:`SIDKY2018104111` which is demonstrated to be robust to user inputs and
requires a minimal number of parameters. It is quite simple to use as described below.

Like Basis Function Sampling, the algorithm proceeds in sweeps. Statistics are
collected over a user specified interval, which is a very flexible choice. At
the end of each sweep an ANN is fit to the statistics to determine the optimal
bias which is then applied in the subsequent sweep. This proceeds until the
free energy landscape has converged.

Detailed information on the method and implementation can be found in
the publication :cite:`SIDKY2018104111`.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

ANN sampling is selected by defining ``"type": "ANN"`` as the
method in the JSON input file. It supports the following options:

topology
	*Array of integers (length: number of hidden layers)*.
	This array defines the architecture of the neural network. ANN sampling
	is quite robust to the choice of network architecture. Nonetheless, there are
	very good heuristics that will ensure you have a network that trains quickly
	and performs well. The key rule is to use a network with a single hidden layer
	if you are biasing on one CV, and two hidden layers if you are biasing on two or
	more. **Please take a look at the example section for more details on network
	architectures**.

grid
	This is a grid object which defines how CV space will be discretized.

nsweep
	*integer*
	Defines the length of a sweep in iterations. Typical values range from
	1,000 to 10,000 depending on the size of the system. The slower the system
	dynamics, the longer the sweep. This is not going to heavily affect
	convergence time, and the method is generally quite robust to choice of
	sweep length. The main consequence of this choice is that the ANN training
	time become relatively expensive if the system is very cheap to evaluate.

weight
	*double*
	(Default = 1)
	This defines how much relative weight is assigned to the statistics
	collected during a sweep. The default value works fine in all cases.
	However, if you know the free energy barriers are very large
	or very small, you can adjust this to speed up convergence.

temperature
	*double*
	The temperature of the simulation must be specified.

output_file
	*string*
	(Default = ``ann.out``)
	This specifies the name of the output file containing the bias and
	reweighted histogram.

overwrite_output
	*boolean*
	(Default = ``true``)
	Setting this to ``false`` will append the sweep number to the output file name
	so that a series of files will be created over time. This is useful if a
	user is interested in plotting the evolution of the bias over time.

lower_bounds
	*array of doubles (nr of CVs) long*.
	This array defines the minimum values for the CV restraints.

upper_bounds
	*array of doubles (nr of CVs) long*.
	This array defines the maximum values for the CV restraints.

lower_bound_restraints
	*array of doubles (nr of CVs) long*.
	This array defines the spring constants for the lower bound CV restraints.
	If you do not wish to restrain the CVs to a particular interval, you can
	set this to zero.

upper_bound_restraints
	*array of doubles (nr of CVs) long*.
	This array defines the spring constants for the upper bound CV restraints.
	If you do not wish to restrain the CVs to a particular interval, you can
	set this to zero.

max_iters
	*integer*
	(Default = 1000)
	Defines the maximum number of training iterations/epochs per sweep. If you
	are running a very large network or a large number of CVs then this can be
	made quite small (10 to 100) since the network will continually improve as
	the system proceeds. Unless you have a very large network, or a very cheap
	system, you can leave this alone.

prev_weight
	*double*
	(Default = 1)
	This is a special feature that allows ANN sampling to "forget" previous
	history specified by a fraction. For example, if you want to retain 90% of
	accumulated data each sweep, set this to 0.9. This is useful if
	you are sampling along a CV you know to be bad, and there is clear
	quasi-nonergodicity in the system. By allowing the ANN method to forget
	accumulated data, you allow it to adapt more efficiently to newly
	accessible regions of phase space. Note that unless you really know why
	you need to be using this, you should leave it alone.

Example Input
^^^^^^^^^^^^^

There is an online repository containing complete examples and analysis for
systems of up to 4 CVs. Along with the original publication, this should
provide you a good idea of how to specify the network architecture. The
repository can be found `here <https://github.com/hsidky/ann_sampling>`_.

For one CV, such as Na--Cl distance, the input is as follows:

.. code-block:: javascript

	"methods" : [
		{
			"type" : "ANN",
			"topology" : [15],
			"nsweep" : 10000,
			"temperature" : 298.15,
			"grid" : {
			"lower" : [0.24],
			"upper" : [0.80],
			"number_points" : [300],
			"periodic": [false]
			},
			"lower_bounds" : [0.23],
			"upper_bounds" : [0.81],
			"lower_bound_restraints" : [0],
			"upper_bound_restraints" : [1000]
		}
	]

for two simultaneous CVs, such as the dihedral angles of ADP, the input is
as follows:

.. code-block:: javascript

	"methods" : [
		{
			"type" : "ANN",
			"topology" : [10, 6],
			"nsweep" : 5000,
			"overwrite_output" : false,
			"temperature" : 298.15,
			"grid" : {
				"lower" : [-3.141592653589793, -3.141592653589793],
				"upper" : [3.141592653589793, 3.141592653589793],
				"number_points" : [30, 30],
				"periodic" : [true, true]
			},
			"lower_bounds" : [-4, -4],
			"upper_bounds" : [4, 4],
			"lower_bound_restraints" : [0, 0],
			"upper_bound_restraints" : [0, 0]
		}
	]

For more examples, and higher dimensions, please check out the repository linked above.

Output
^^^^^^

ANN sampling writes either a single output file or a series of output files
over time. Each file contains columns corresponding to the CVs, a column
containing the unbiased histogram estimate and a final column containing the
bias. The format is as follows:

*cv1 cv2 ... histogram bias*

This file can be loaded and visualized easily in many scripting languages, such
as Python and MATLAB. An example of how to load data in Python for a 2D CV is
shown below.

.. code-block:: python

	# Load data.
	X = np.loadtxt("ann.dat")
	xg = np.reshape(X[:,0], (61, 61))
	yg = np.reshape(X[:,1], (61, 61))
	zg = np.reshape(-X[:,3], (61, 61))
	zg = zg - np.max(zg)

	# Plot data.
	fig = plt.figure(figsize=(5,5))
	plt.contour(xg, yg, zg, linewidths=0.5, colors="k")
	plt.contourf(xg, yg, zg)

A file called "netstate.dat" is also written out which contains the neural
network parameters. This network can be evaluated in Python using a ANN library
such as `TensorFlow <https://www.tensorflow.org/>`_ or
`Keras <https://keras.io/>`_.

.. code-block:: python

	from keras.models import Sequential
	from keras.layers import Dense, Activation

	# Import and define Keras network.
	params = []
	xshift = []
	xscale = []
	yshift = []
	yscale = []
	net = Sequential()
	with open("netstate.dat", "r") as f:
		# Topology.
		layers = int(f.readline())
		arch = [int(x) for x in f.readline().split()]

		# Scaling and shifting.
		xscale = [float(x) for x in f.readline().split()]
		xshift = [float(x) for x in f.readline().split()]
		yscale = [float(x) for x in f.readline().split()]
		yshift = [float(x) for x in f.readline().split()]

		# Weights and biases.
		for i in range(1, layers):
			b = []
			for j in range(arch[i]):
				b.append(float(f.readline()))
			b = np.array(b)

			w = []
			for j in range(arch[i]*arch[i-1]):
				w.append(float(f.readline()))
			w = np.array(w).reshape(arch[i-1], arch[i])

			params.append(w)
			params.append(b)

			if i==1:
				net.add(Dense(arch[i], activation="tanh", input_dim=arch[i-1]))
			elif i==layers-1:
				net.add(Dense(arch[i], activation="linear"))
			else:
				net.add(Dense(arch[i], activation="tanh"))

	net.set_weights(params)

The network can then be evaluated on a high-resolution grid and plotted.

.. code-block:: python

	# Define new high-resolution grid.
	x = np.linspace(-np.pi, np.pi, 500, endpoint=True)
	y = np.linspace(-np.pi, np.pi, 500, endpoint=True)
	xg, yg = np.meshgrid(x, y)

	# Scale data.
	xs = np.vstack((xg.flatten(), yg.flatten())).T
	xs = (xs - xshift)*xscale

	# Evaluate network. Unscale data.
	ys = net.predict(xs)
	ys = ys/yscale + yshift
	zg = -ys.reshape(500, 500)

	# Plot data.
	plt.figure(figsize=(12,10))
	zg = zg - np.max(zg)
	plt.contour(xg, yg, zg, linewidths=0.5, colors="k")
	plt.contourf(xg, yg, zg)
	cb = plt.colorbar()
	cb.set_label("G (kJ/mol)")
	plt.xlabel("$\phi$")
	plt.ylabel("$\psi$")

These examples and more are also found in the `online repository <https://github.com/hsidky/ann_sampling>`_.

Developer
^^^^^^^^^

* Hythem Sidky

.. warning::

	Please make sure to cite the paper :cite:`SIDKY2018104111` if you use
	this method!
