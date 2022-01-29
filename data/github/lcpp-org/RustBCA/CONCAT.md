# `RustBCA`

`RustBCA` is a general-purpose, high-performance code for simulating
ion-material interactions including sputtering, reflection, and implantation
using the binary collision approximation ([BCA]), written in [Rust]!
RustBCA consists of a standalone code and libraries for including
ion-material interactions in simulations written in C/C++, Python,
and Fortran.

By discretizing the collision cascade into a sequence of binary collisions,
[BCA] codes can accurately and efficiently model the prompt interaction
between an energetic ion and a target material. This includes reflection,
implantation, and transmission of the incident ion, s well as sputtering
and displacement damage of the target. Generally, [BCA] codes can be
valid for incident ion energies between approximately ~1 eV/nucleon 
to <1 GeV/nucleon.

Check out the `RustBCA` [Wiki] for detailed information, installation
instructions, use cases, examples, and more. See the RustBCA paper at the
Journal of Open Source Software by clicking the badge below:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03298/status.svg)](https://doi.org/10.21105/joss.03298)

## Getting started

For those eager to get started, try running one of the examples in the
`rustBCA` directory. Note that to automatically manipulate input files and reproduce the plots located on the [Wiki], these require several optional, but common,
[Python] packages (`matplotlib`, `numpy`, `scipy`, `shapely`, and `toml`).

### H trajectories and collision cascades in a boron nitride dust grain

First, run the example using:

```bash
cargo run --release examples/boron_nitride.toml
```

Afterwords, fire up your favourite [Python] interpreter
(e.g., [IPython]) and execute:

```python
from scripts.rustbca import *
do_trajectory_plot("boron_dust_grain_")
```

### He implantation into a layered TiO<sub>2</sub>/Al/Si target

First, run the example using:

```bash
cargo run --release examples/layered_geometry.toml
```

Afterwords, fire up your favourite [Python] interpreter
(e.g., [IPython]) and execute:

```python
import numpy as np
import matplotlib.pyplot as plt

deposited_ions = np.genfromtxt(
    "2000.0eV_0.0001deg_He_TiO2_Al_Sideposited.output",
    delimiter=",",
    names=["M", "Z", "x", "y", "z", "collisions"],
)

plt.hist(deposited_ions["x"], bins=100)

plt.show()
```

## Features

The following features are implemented in `rustBCA`:

* Ion-material interactions for all combinations of incident ion and target species.
* Infinite, homogeneous targets (Mesh0D), Layered, finite-depth inhomogeneous targets (Mesh1D), arbitrary 2D geometry composition through a triangular mesh (Mesh2D), homogeneous spherical geometry (Sphere) and homogeneous, arbitrary triangular mesh geometry (TriMesh).
* Amorphous Solid/Liquid targets, Gaseous targets, and targets with both solid/liquid and gaseous elements
* Low energy (< 25 keV/nucleon) electronic stopping modes including:
  * local (Oen-Robinson),
  * nonlocal (Lindhard-Scharff),
  * and equipartition forms.
* Biersack-Varelas interpolation is also included for electronic stopping up to ~1 GeV/nucleon. Note that high energy physics beyond electronic stopping are not included.
* Optionally, the Biersack-Haggmark treatment of high-energy free-flight paths between collisions can be included to greatly speed up high-energy simulations (i.e., by neglecting very small angle scattering).
* A wide range of interaction potentials are provided, including:
  * the Kr-C, ZBL, Lenz-Jensen, and Moliere universal, screened-Coulomb potentials.
  * the Lennard-Jones 12-6, Lennard-Jones 6.5-6, and Morse attractive-repulsive potentials.
* Solving the distance-of-closest-approach problem is achieved using:
  * the Newton-Raphson method for simple root-finding,
  * or, for attractive-repulsive potentials, an Adaptive Chebyshev Proxy Rootfinder with Automatic Subdivision algorithm and a Polynomial root-finding algorithm are provided through the [rcpr] crate.
* Multiple interaction potentials can be used in a single simulation for any number of potentials/species.
  * For example, the He-W interaction can be specified using a Lennard-Jones 12-6 potential, while the W-W interaction can be defined using a Kr-C potential.
* The scattering integral can be calculated using:
  * Gauss-Mehler quadrature,
  * Gauss-Legendre quadrature,
  * Mendenall-Weller quadrature,
  * or the MAGIC algorithm.
* Input files use the [TOML] format, making them both human-readable and easily parsable.
* RustBCA generates user-friendly, context-providing error messages, which help pinpoint the cause of errors and provide suggested fixes to the user.
* The simulation results are comma-delimited (`csv` format) and include:
  * the energies and directions of emitted particles (reflected ions and sputtered atoms),
  * the final positions of implanted ions,
  * full trajectory tracking for both the incident ions and target atoms,
  * and many other parameters such as position of origin of sputtered particles and energy loss along trajectories.
* Optionally, the code can produce energy-angle and implantation distributions when built with the `--features distributions` flag and disable space-intensive particle list output with `--features no_list_output`.
* Library functions for modeling ion reflection, implantation, and sputtering in C++/C, Python, and Fortran codes.

## Installation

Without optional features, `rustBCA` should compile with `cargo` alone on
Windows, MacOS, and Linux systems.

[HDF5] for particle list input has been tested on Windows, but version 1.10.6 must be used.
[rcpr], the adaptive Chebyshev Proxy Rootfinder with automatic subdivision and
polynomial rootfinder package for [Rust], has not yet been successfully compiled
on Windows.
However, it can be compiled on the Windows Subsystem for Linux (WSL) and, likely,
on Ubuntu for Windows or Cygwin.

#### Manual Dependences

* [rustup], the [Rust] toolchain (includes `cargo`, the [Rust] package manager, `rustc`, the [Rust] compiler, and more).

#### Automatic Dependencies

* see [Cargo.toml](https://github.com/lcpp-org/RustBCA/blob/master/Cargo.toml) for a complete list.

#### Optional Dependencies

* [HDF5] libraries
* [rcpr], a CPR and polynomial rootfinder, required for using attractive-repulsive interaction potentials such as Lennard-Jones or Morse. It may require additional software (see below).
* For manipulating input files and running associated scripts, the following are required:
  * [Python] 3.6+
  * The [Python] libraries: `numpy`, `matplotlib`, `toml` (must build from source), `shapely`, and `scipy`.

### Detailed instructions for Ubuntu 18.04 LTS

1. (Optional) Install Python 3.6+ (this comes natively in Ubuntu 18.04)
2. Install `curl`:
```bash
sudo apt-get install curl
```
3. Install [rustup], the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/ by running the following command and following on-screen instructions:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
4. (Optional) Install `pip` for [Python]:
```bash
sudo apt-get install python3-pip
```
5. (Optional) Install [Python] libraries for making input files:
```bash
python3 -m pip install numpy matplotlib shapely scipy
```
6. (Optional) Install [Python] [TOML] library from source:
```bash
git clone https://github.com/uiri/toml.git
cd toml
python3 setup.py install
```
7. (Optional) Install software for [rcpr]:
```bash
sudo apt-get install gcc gfortran build-essential cmake liblapack-dev libblas-dev liblapacke-dev
```
8. Install `cargo`:
```bash
sudo apt-get install cargo
```
9. Build `rustBCA`:
```bash
git clone https://github.com/lcpp-org/rustBCA
cd rustBCA
cargo build --release
```
10. (Optional) Build `rustBCA` with optional dependencies, `hdf5` and/or `rcpr` (with your choice of backend: `openblas`, `netlib`, or `intel-mkl`):
```bash
cargo build --release --features cpr_rootfinder_netlib,hdf5_input
cargo build --release --features cpr_rootfinder_openblas,hdf5_input
cargo build --release --features cpr_rootfinder_intel_mkl,hdf5_input
 ```
11. `input.toml` is the input file - see [Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files) for more information
12. Run the required tests using:
```bash
cargo test
```
13. (Optional) Run the required and optional tests for the desired backend(s):
```bash
cargo test --features cpr_rootfinder_netlib
cargo test --features cpr_rootfinder_openblas
cargo test --features cpr_rootfinder_intel_mkl
```

### Detailed instructions for Fedora 33

Most of the ingredients for building `rustBCA` and running the [Python] helper
scripts are available natively in the Fedora software repository, so the setup
is relatively painless.

The [Rust] toolchain can be aquired using:

```bash
sudo dnf install rust rust-src rust-std-static rust-analysis rust-gdb rust-lldb rustfmt
```

The (optional) [Python] packages can be obtained using:

```bash
sudo dnf install python3-numpy python3-scipy python3-matplotlib python3-toml python3-shapely
```

or, alternatively, using `pip3`.

If the [rcpr] is desired, it's probably also a good idea to install the following:

```bash
sudo dnf install gcc gcc-gfortran cmake lapack lapack-devel blas blas-devel
```

Building `rustBCA` is straightforward and can be done using:

```bash
git clone https://github.com/lcpp-org/rustBCA
cd RustBCA
cargo build --release
```

with all of the explicit dependencies listed in `Cargo.toml` handled
automatically during the build.

## Usage

To use `rustBCA`, modify the `input.toml` file, which is used to configure each
simulation.
To run a simulation, execute:

```bash
./rustBCA
```

with `input.toml` in the same directory as `rustBCA`.
Alternatively, `rustBCA` accepts the name of a`.toml` input file as a single
command line argument:

```bash
./rustBCA /path/to/input.toml
```

For further details, have a look at
[Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files)
on the `rustBCA` [Wiki] for usage instructions.
Also have a look at the examples on the [Wiki] for writing `.toml` input files.

[BCA]: https://en.wikipedia.org/wiki/Binary_collision_approximation
[HDF5]: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
[IPython]: https://en.wikipedia.org/wiki/IPython
[Python]: https://en.wikipedia.org/wiki/Python_(programming_language)
[rcpr]: https://github.com/drobnyjt/rcpr
[rustup]: https://rustup.rs
[Rust]: https://en.wikipedia.org/wiki/Rust_(programming_language)
[TOML]: https://en.wikipedia.org/wiki/TOML
[Wiki]: https://github.com/lcpp-org/RustBCA/wiki
---
title: 'RustBCA: A High-Performance Binary-Collision-Approximation Code for Ion-Material Interactions'
tags:
  - Rust
  - plasma material interactions
  - binary collision approximation
  - ion solid interactions
  - ion material interactions
  - sputtering
  - reflection
  - implantation
authors:
  - name: Jon. T Drobny
    orcid: 0000-0002-9733-6058
    affiliation: 1
  - name: Davide Curreli
    affiliation: 1
affiliations:
  - name: Department of Nuclear, Plasma, and Radiological Engineering, University of Illinois at Urbana-Champaign
    index: 1
date: 24 October 2020
bibliography: paper.bib
---

# Summary

Ion-material interactions are of vital importance in industrial applications, the study and design of nuclear fusion devices, the engineering of survivable spacecraft components, and more. In particular, plasma-material interactions are typically dominated by ion-material interactions, including the phenomena of sputtering, reflection, and implantation. These phenomena are difficult to model analytically, and many such models rely on empirical or semi-empirical formulas, such as the Yamamura sputtering yield formula [@Yamamura1982], or the Thomas reflection coefficient [@Thomas1992]. However, such models are inherently limited, and of little applicability to complex geometry, multi-component surfaces, or for coupling to plasma or material dynamics codes. Since ion-material interactions span a range of energies from sub-eV to GeV and beyond, n-body approaches such as molecular dynamics can be computationally infeasible for many applications where the characteristic ion range exceeds the limits of reasonable molecular dynamics domains. Instead, approximations to the full n-body problem are used; the most common of these is the Binary Collision Approximation (BCA), a set of simplifying assumptions to the full n-body problem. RustBCA is a high-performance, general purpose, ion-material interactions BCA code, built for scientific flexibility and ease of use. RustBCA features include:

 - electronic stopping formulations for low energy (up to 25 keV/nucleon) and high energy (up to 1 GeV/nucleon)
 - Kr-C, ZBL, Moliere, and Lenz-Jensen screened coulomb interatomic potentials
 - Lennard-Jones and Morse attractive-repulsive potentials
 - the unique capability of using multiple interatomic potentials in a single simulation
 - choice of Gaussian quadrature or the approximate MAGIC algorithm for determining scattering angles
 - full trajectory tracking of ions and material atoms, including local nuclear and electronic energy losses
 - a human- and machine-readable configuration file
 - full 6D output of all particles that leave the simulation (via sputtering or reflection)
 - multiple geometry types

# Binary Collision Approximation Codes

RustBCA is an amorphous-material BCA code, following the TRIM [@Biersack1980] family of codes, which includes Tridyn [@Möller1988], SDTrimSP [@Mutzke2019], F-TRIDYN [@Drobny2017], and SRIM [@Ziegler2010]; this has historically been the most popular implementation of the BCA. Based on the number of citations recorded in Google Scholar at the time of writing, SRIM is the most popular amorphous-material BCA code, likely due to its being free to download, available on Windows, and having a graphical user interface. It is followed by the original TRIM code, upon which SRIM was based, then Tridyn, and finally SDTrimSP. Crystalline-material BCA codes have also been developed, such as MARLOWE [@Robinson1974], OKSANA [@Shulga1984],  and some versions of ACAT [@Yamamura1996], but are not as widely used. The BCA itself is a set of simplifying assumptions for the ion-material interaction problem; the assumptions used in the amorphous-material BCA can be summarized as follows:

* Particles in the code, ions and material atoms both, are "superparticles" that represent many real ions or atoms each
* Energetic particles interact with initially stationary atoms in the material through elastic, binary collisions
* Collisions occur at mean-free-path lengths, or exponentially distributed path lengths for gaseous targets
* Particle trajectories are approximated by the classical asymptotic trajectories
* Electronic interactions are separated from the nuclear, elastic interactions
* Local electronic energy losses occur at each collision
* Nonlocal electronic energy losses occur along each segment of the asymptotic trajectories
* Material atoms are mobile and transfer momentum following collisions
* Particles are stopped when their energy drops below a threshold, cutoff energy, or when they leave the simulation as sputtered or reflected/transmitted particles
* Particles that leave a surface experience reflection by or refraction through a locally planar surface binding potential
* When simulating radiation damage, only material atoms given an energy larger than the threshold displacement energy will be considered removed from their original location

For detailed summaries of the history and theory of binary collision approximation codes, see the review by Robinson [@Robinson1994] and the text by Eckstein [@Eckstein1991].

# Statement of Need

Ion-material interactions have been historically modeled using analytical and semi-empirical formulas, such as Sigmund's sputtering theory [@Sigmund1987], the Bohdansky formula [@Bohdansky1980; @Bohdansky1984], the Yamamura formula [@Yamamura1982; @Yamamura1983; @Yamamura1984], and the Thomas et al. reflection coefficient [@Thomas1992]. However, for any physical situation beyond the regimes of validity of these formulas (e.g., non-normal angles of incidence), or for complex geometry, or for inhomogeneous composition, straightforward empirical formulas cannot be reliably used. Many BCA codes have been developed to provide computationally efficient solutions to these problems, including SRIM [@Ziegler2010], Tridyn [@Möller1988], F-TRIDYN [@Drobny2017], SDTrimSP [@Mutzke2019] and its derivatives, which are based on the original TRIM [@Biersack1980] code. However, each has limitations that prevent widespread adoption across a broad range of applications. In particular, SRIM, which is free-use but closed-source, suffers from relatively poor computational performance and significant anomalies in sputtered atom angular distributions and light ion sputtering yields [@Shulga2019; @Shulga2018; @Hofsass2014; @Wittmaack2004]. Tridyn and F-TRIDYN, which are not open source, are limited to low ion energy, specific screened-coulomb potentials, mono-angular ion beams, atomically flat and atomically rough surfaces respectively, and are single-threaded. SDTrimSP, although significantly more advanced than the preceding codes, is built on the original TRIM source code and is not open-source.

As far as the authors are aware, there is no widely-used open-source BCA code suitable for simulating plasma-material interactions. Iradina is an open source BCA that has been used for ion-material interactions in a semicondcutor manufacturing context [@HollandMoritz2017; @Johannes2014], but sputtering yields, reflection coefficients, or other key quantities of interest from iradina have not, to the knowledge of the authors, been reported for a wide range of ions, targets, energies, or angles. Additionally, those BCA codes that are available, through licensing agreements or as closed-source software, are not well suited to a wide range of physical problems. Particularly, the direct integration of BCA codes to particle or subsurface dynamics codes, such as those performed using F-TRIDYN for ITER divertor simulations [@Lasa2020], requires costly external wrappers to manage simulations and process output files to perform file-based coupling. RustBCA, as part of the  [Plasma Surface Interactions 2 SciDAC Project](https://collab.cels.anl.gov/display/PSIscidac2/Plasma+Surface+Interactions+2) suite of codes, has been developed to fill that gap and expand upon the feature set included in currently available BCA codes.

Features unique to RustBCA include the ability to handle attractive-repulsive interatomic potentials, use multiple interatomic potentials in one simulation, handle high-energy incident ions and multiple geometry types, use large file input of incident particles to facilitate coupling to other codes via HDF5, output pre-binned distributions without post-processing of text-based particle lists, and use a human- and machine-readable configuration file. RustBCA has been designed with modern programming techniques, robust error-handling, and multi-threading capability. RustBCA is being developed as both a standalone code and as a library code that may be used to add BCA routines to other high-performance codes to avoid file-based code coupling entirely. Additionally, the TRIM family of codes typically relies on the MAGIC algorithm to approximate the scattering integral with 5 fitting coefficients. RustBCA includes not only an implementation of the MAGIC algorithm, but also Mendenhall-Weller, Gauss-Mehler, and Gauss-Legendre quadrature, the three of which are significantly more accurate than the MAGIC algorithm. We hope that giving users direct access to a user-friendly, flexible, high-performance, open-source BCA will encourage and enable heretofore unexplored research in ion-materials interactions.

![Figure showing sputtering yields of silicon from SRIM, RustBCA, F-TRIDYN, Yamamura's formula for Q=0.33-0.99, and a smooth analytical fit to experimental data by Wittmaack [@Wittmaack2004], for an incident energy of 1 keV and for many different projectiles.](corrected_yields.png)

Quantities of interest from RustBCA, including sputtering yields, have been benchmarked against F-TRIDYN, SRIM, empirical formulas, and experiments. This summary figure shows the sputtering yields of silicon by 1 keV helium, beryllium, oxygen, neon, aluminum, silicon, argon, titanium, copper, krypton, xenon, ytterbium, tungsten, gold, lead and uranium ions. SRIM's unphysical Z1 dependence is clearly visible, as is the divergence of Yamamura's formula (for Q = 0.66, the reported value for silicon, and +/- 0.33) at high mass ratios (M1 >> M2) from the experimental data collected by Wittmaack [@Wittmaack2004]. RustBCA and F-TRIDYN both reproduce the correct Z1 dependence of the sputtering yield, and correctly model the magnitude of the yield for all projectiles. It should be noted that, for this simulation, F-TRIDYN uses corrected MAGIC coefficients [@Ziegler2010], that differ from those originally included in the Tridyn source code; Tridyn's orignal MAGIC coefficients underestimated the sputtering yield for high mass ratios. A soft grey line depicts the point of silicon on silicon sputtering. Reflection coefficients, although very low for mass ratios above one, are also shown, with F-TRIDYN and RustBCA agreeing with the semi-empirical Thomas reflection coefficient formula.

# Examples

RustBCA includes multiple example input files, under the examples/ folder on the directory, as well as discussion of each on the RustBCA [Wiki](https://github.com/lcpp-org/RustBCA/wiki) page. Three examples will be summarized here.

## Example 1: Layered Targets

First, an example of 2 keV helium ions at normal incidence on a layered titanium dioxide, aluminum, and silicon target can be run in 2D with:

```
cargo run --release examples/layered_geometry.toml
```

 The same example using the 1D layered geometry can be run with:

```
cargo run --release 1D examples/layered_geometry_1D.toml
```

 ![Helium implantation depth distributions at 2 keV in a layered TiO2-Al-Si target.](layered_target.png)

 The depth distribution, compared to F-TRIDYN [@Drobny2017], clearly shows the effect of layer composition and sharp interfaces on the combined nuclear and electronic stopping of helium.

## Example 2: 2D Geometry

 Second, as an example of the capability of RustBCA to handle 2D geometry, the trajectories of 1 keV hydrogen on a circular cross-section of boron-nitride can be simulated.

```
cargo run --release examples/boron_nitride.toml
```

 ![Trajectories of hydrogen and mobile boron and nitrogen resulting from 10 1 keV hydrogen ions impacting on a circular cross-section boron-nitride target.](H_B_N.png)

## Example 3: Spherical geometry

Third, the 2D boron nitride example can be run as a spherical boron nitride dust grain, by running the following command:

```
cargo run --release SPHERE examples/boron_nitride_sphere.toml
```

The trajectories can be plotted in 3D with mayavi using `do_trajectory_plot_3d()` or in 2D with matplotlib using `do_trajectory_plot()` in `scripts/rustbca.py`.

![Trajectories of hydrogen, boron, and nitrogen in a 3D boron target. Hydrogen is medium blue, nitrogen yellow, and boron light blue.](sphere_trajectories_bordered.png)

# Acknowledgements
This work was funded by the U.S. Department of Energy, Office of Fusion Energy Sciences through the Scientific Discovery through Advanced Computing (SciDAC) project on Plasma Surface Interactions 2 (Grant No. DE-SC0018141).

# References
Thank you for your contribution to rustBCA!

Before submitting this PR, please make sure you have:

- [ ] Opened an issue
- [ ] Referenced the relevant issue number(s) below
- [ ] Provided a description of the changes below
- [ ] Ensured all tests pass and added any necessary tests for new code

Fixes # (issue)

## Description
Please include a concise description of the change and how it addresses a relevant issue, including references, input files, and figures as relevant.

## Tests
Please describe how the changes in this pull request have been tested, including system information such as OS.
---
name: Bug report
about: Report a problem with the code
title: "[bug]"
labels: bug
assignees: ''

---

**Description**
A clear and concise description of what the bug is.

**To Reproduce**
Please attach the relevant TOML input file and list the command line arguments used.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Error messages, output files, or figures**
Please include error messages, output files or figures that demonstrate the error if applicable.

**System (please complete the following information):**
OS: 

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for rustBCA
title: "[feature]"
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is.

**Proposed solution**
A clear and concise description of what you want to happen.

**Alternative solution(s)**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context, including but not limited to references, figures, or examples.
---
name: Question
about: Question about implementation, physics, or features
title: "[question]"
labels: question
assignees: ''

---


